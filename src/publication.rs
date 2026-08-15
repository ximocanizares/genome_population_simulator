//! Ensamblado provisional y publicación Linux atómica sin sobrescritura.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File, OpenOptions};
use std::io::{Read, Write};
use std::os::fd::OwnedFd;
use std::os::unix::ffi::OsStrExt;
use std::path::{Component, Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use rustix::fs::{
    AtFlags, CWD, Dir, FileType, Mode, OFlags, RenameFlags, fstat, openat, renameat_with, statat,
};

use crate::diagnostics::Diagnostic;
use crate::integrity::{ArtifactIdentity, Sha256Writer, identify_file, verify_identity};
use crate::results::{
    ArtifactRecord, UnitRecord, UnitState, ValidatedGeneralManifest,
    bind_general_manifest_to_workspace, parse_general_manifest_bytes,
    validate_general_manifest_bytes,
};

const INCOMPLETE_MARKER: &str = "INCOMPLETE";
static WORKSPACE_COUNTER: AtomicU64 = AtomicU64::new(0);

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ArtifactStage {
    Open,
    Closed,
    Validated,
    Cataloged,
    Failed,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ArtifactLifecycle {
    id: String,
    stage: ArtifactStage,
}

impl ArtifactLifecycle {
    #[must_use]
    pub fn new(id: impl Into<String>) -> Self {
        Self {
            id: id.into(),
            stage: ArtifactStage::Open,
        }
    }

    #[must_use]
    pub const fn stage(&self) -> ArtifactStage {
        self.stage
    }

    pub fn close(&mut self) -> Result<(), Diagnostic> {
        self.transition(ArtifactStage::Open, ArtifactStage::Closed)
    }

    pub fn validate(&mut self) -> Result<(), Diagnostic> {
        self.transition(ArtifactStage::Closed, ArtifactStage::Validated)
    }

    pub fn catalog(&mut self) -> Result<(), Diagnostic> {
        self.transition(ArtifactStage::Validated, ArtifactStage::Cataloged)
    }

    pub fn fail(&mut self) -> Result<(), Diagnostic> {
        if matches!(
            self.stage,
            ArtifactStage::Open | ArtifactStage::Closed | ArtifactStage::Validated
        ) {
            self.stage = ArtifactStage::Failed;
            return Ok(());
        }
        Err(invalid_unit_state(
            &self.id,
            format!("a {:?} artifact cannot transition to failed", self.stage),
        ))
    }

    fn transition(
        &mut self,
        required: ArtifactStage,
        next: ArtifactStage,
    ) -> Result<(), Diagnostic> {
        if self.stage == required {
            self.stage = next;
            return Ok(());
        }
        Err(invalid_unit_state(
            &self.id,
            format!("expected {required:?}, found {:?}", self.stage),
        ))
    }
}

#[derive(Clone, Debug)]
pub struct PublishableUnit {
    id: String,
    kind: String,
    required_artifacts: BTreeMap<String, ArtifactStage>,
    state: UnitState,
    diagnostic_code: Option<String>,
}

impl PublishableUnit {
    #[must_use]
    pub fn new(
        id: impl Into<String>,
        kind: impl Into<String>,
        required_artifacts: impl IntoIterator<Item = String>,
    ) -> Self {
        Self {
            id: id.into(),
            kind: kind.into(),
            required_artifacts: required_artifacts
                .into_iter()
                .map(|id| (id, ArtifactStage::Open))
                .collect(),
            state: UnitState::Incomplete,
            diagnostic_code: None,
        }
    }

    pub fn record_artifact(&mut self, id: &str, stage: ArtifactStage) -> Result<(), Diagnostic> {
        if self.state != UnitState::Incomplete {
            return Err(invalid_unit_state(
                &self.id,
                format!("a {:?} unit cannot change its artifact catalog", self.state),
            ));
        }
        let known = self
            .required_artifacts
            .get_mut(id)
            .ok_or_else(|| invalid_unit_state(&self.id, format!("unexpected artifact {id}")))?;
        if *known == ArtifactStage::Cataloged || *known == ArtifactStage::Failed {
            return Err(invalid_unit_state(
                &self.id,
                format!("artifact {id} is already terminal at {known:?}"),
            ));
        }
        if stage != ArtifactStage::Failed && stage_rank(stage) < stage_rank(*known) {
            return Err(invalid_unit_state(
                &self.id,
                format!("artifact {id} cannot regress from {known:?} to {stage:?}"),
            ));
        }
        *known = stage;
        Ok(())
    }

    pub fn complete(&mut self) -> Result<(), Diagnostic> {
        if self.state != UnitState::Incomplete {
            return Err(invalid_unit_state(
                &self.id,
                format!("unit is already {:?}", self.state),
            ));
        }
        if self
            .required_artifacts
            .values()
            .all(|stage| *stage == ArtifactStage::Cataloged)
        {
            self.state = UnitState::Complete;
            return Ok(());
        }
        Err(invalid_unit_state(
            &self.id,
            "not every required artifact is closed, validated and cataloged",
        ))
    }

    pub fn fail(&mut self, diagnostic: &Diagnostic) -> Result<(), Diagnostic> {
        if self.state != UnitState::Incomplete {
            return Err(invalid_unit_state(
                &self.id,
                format!("unit is already {:?}", self.state),
            ));
        }
        self.state = UnitState::Failed;
        self.diagnostic_code = Some(diagnostic.code().to_owned());
        Ok(())
    }

    #[must_use]
    pub fn record(&self) -> UnitRecord {
        UnitRecord {
            id: self.id.clone(),
            kind: self.kind.clone(),
            state: self.state,
            exportable: self.state == UnitState::Complete,
            diagnostic_code: self.diagnostic_code.clone(),
        }
    }
}

pub(crate) trait AtomicPublisher {
    /// Hook controlado que se ejecuta antes de la última revalidación física.
    ///
    /// La implementación de producción no lo modifica. Los dobles de prueba
    /// pueden alterar un artefacto aquí para demostrar que la ventana entre una
    /// comprobación anterior y el `rename` queda cubierta.
    fn before_final_validation(&self, _provisional: &Path) -> Result<(), Diagnostic> {
        Ok(())
    }

    /// Publica con una única operación atómica no-clobber. Si devuelve error,
    /// debe dejar el provisional en su ruta original; nunca puede retirar ni
    /// interpretar el marcador `INCOMPLETE`.
    fn publish_no_replace(&self, provisional: &Path, destination: &Path) -> Result<(), Diagnostic>;
}

#[derive(Debug, Default)]
pub(crate) struct LinuxAtomicPublisher;

impl AtomicPublisher for LinuxAtomicPublisher {
    fn publish_no_replace(&self, provisional: &Path, destination: &Path) -> Result<(), Diagnostic> {
        renameat_with(CWD, provisional, CWD, destination, RenameFlags::NOREPLACE).map_err(|error| {
            if error == rustix::io::Errno::EXIST {
                Diagnostic::io(
                    "GPS-IO-003",
                    "El destino final ya existe y no puede sobrescribirse.",
                )
                .with("path", destination.display().to_string())
            } else if matches!(
                error,
                rustix::io::Errno::NOSYS | rustix::io::Errno::INVAL | rustix::io::Errno::OPNOTSUPP
            ) {
                unsupported_publication(destination, error.to_string())
            } else {
                Diagnostic::io("GPS-IO-001", "Una operación de lectura o escritura falló.")
                    .with("path", destination.display().to_string())
                    .with("operation", "renameat2(RENAME_NOREPLACE)")
                    .with("detail", error.to_string())
            }
        })
    }
}

/// Doble determinista para demostrar el rechazo de sistemas sin no-clobber.
#[derive(Debug, Default)]
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) struct UnsupportedAtomicPublisher;

impl AtomicPublisher for UnsupportedAtomicPublisher {
    fn publish_no_replace(
        &self,
        _provisional: &Path,
        destination: &Path,
    ) -> Result<(), Diagnostic> {
        Err(unsupported_publication(
            destination,
            "simulated_unsupported",
        ))
    }
}

#[derive(Debug)]
pub struct ExecutionWorkspace {
    destination: PathBuf,
    provisional: PathBuf,
    artifacts: Vec<ArtifactRecord>,
    token: u64,
}

impl ExecutionWorkspace {
    pub fn create(destination: &Path) -> Result<Self, Diagnostic> {
        let parent = destination.parent().unwrap_or_else(|| Path::new("."));
        let name = destination
            .file_name()
            .and_then(|name| name.to_str())
            .filter(|name| !name.is_empty())
            .ok_or_else(|| invalid_unit_state("execution", "destination has no file name"))?;
        fs::create_dir_all(parent).map_err(|error| io_error(parent, "create_dir_all", error))?;
        for _ in 0..128 {
            let counter = WORKSPACE_COUNTER.fetch_add(1, Ordering::Relaxed);
            let provisional = parent.join(format!(
                ".{name}.gps-incomplete-{}-{counter}",
                std::process::id()
            ));
            match fs::create_dir(&provisional) {
                Ok(()) => {
                    fs::write(provisional.join(INCOMPLETE_MARKER), b"incomplete\n")
                        .map_err(|error| io_error(&provisional, "write_marker", error))?;
                    return Ok(Self {
                        destination: destination.to_path_buf(),
                        provisional,
                        artifacts: Vec::new(),
                        token: counter,
                    });
                }
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(error) => return Err(io_error(&provisional, "create_provisional", error)),
            }
        }
        Err(invalid_unit_state(
            "execution",
            "could not allocate a unique provisional sibling",
        ))
    }

    #[must_use]
    pub(crate) fn root(&self) -> &Path {
        &self.provisional
    }

    #[must_use]
    pub fn destination(&self) -> &Path {
        &self.destination
    }

    #[must_use]
    pub fn artifacts(&self) -> &[ArtifactRecord] {
        &self.artifacts
    }

    /// Valida el manifiesto y lo liga a esta instancia y a su catálogo actual.
    pub fn validate_general_manifest(
        &self,
        path: &Path,
        bytes: &[u8],
    ) -> Result<ValidatedGeneralManifest, Diagnostic> {
        let manifest_path = self.provisional.join(path);
        let manifest = bind_general_manifest_to_workspace(
            &manifest_path,
            bytes,
            self.token,
            self.artifacts.clone(),
        )?;
        self.validate_manifest_workspace_artifacts(&manifest)?;
        self.validate_physical_tree(false)?;
        Ok(manifest)
    }

    pub(crate) fn mark_existing_artifacts_failed(&mut self, ids: &BTreeSet<String>) {
        for artifact in &mut self.artifacts {
            if ids.contains(&artifact.id) {
                artifact.state = UnitState::Failed;
                artifact.validated = false;
            }
        }
    }

    pub fn write_artifact(
        &mut self,
        id: &str,
        kind: &str,
        relative_path: &Path,
        bytes: &[u8],
    ) -> Result<ArtifactRecord, Diagnostic> {
        if kind == "manifest_general" {
            return Err(invalid_unit_state(
                id,
                "a general manifest must use the semantic validation gate",
            ));
        }
        self.write_artifact_bytes(id, kind, relative_path, bytes)
    }

    /// Valida, escribe y relee el manifiesto general antes de catalogarlo.
    ///
    /// Esta es la única vía que admite `manifest_general`; evita que un nuevo
    /// productor pueda omitir la segunda puerta de enlaces dinámicos.
    pub fn write_general_manifest(
        &mut self,
        manifest: &ValidatedGeneralManifest,
    ) -> Result<ArtifactRecord, Diagnostic> {
        if manifest.workspace_token() != self.token
            || manifest.workspace_snapshot() != self.artifacts
        {
            return Err(invalid_unit_state(
                "manifest-general",
                "the general-manifest capability belongs to another workspace or stale catalog",
            ));
        }
        self.validate_manifest_workspace_artifacts(manifest)?;
        self.validate_physical_tree(false)?;
        let relative_path = Path::new("manifest.json");
        let bytes = manifest.bytes();
        let record = self.write_artifact_bytes(
            "manifest-general",
            "manifest_general",
            relative_path,
            bytes,
        )?;
        let serialized =
            read_workspace_artifact_nofollow(&self.provisional, relative_path, &record)?;
        let manifest_path = self.provisional.join(relative_path);
        validate_general_manifest_bytes(&manifest_path, &serialized)?;
        let (_, reloaded_artifacts) = parse_general_manifest_bytes(&manifest_path, &serialized)?;
        compare_artifact_catalogs(&reloaded_artifacts, manifest.artifacts())?;
        self.validate_physical_tree(true)?;
        Ok(record)
    }

    fn write_artifact_bytes(
        &mut self,
        id: &str,
        kind: &str,
        relative_path: &Path,
        bytes: &[u8],
    ) -> Result<ArtifactRecord, Diagnostic> {
        validate_relative_path(relative_path)?;
        let final_path = self.provisional.join(relative_path);
        let parent = final_path
            .parent()
            .ok_or_else(|| invalid_unit_state(id, "artifact path has no provisional parent"))?;
        fs::create_dir_all(parent).map_err(|error| io_error(parent, "create_dir_all", error))?;
        let file_name = final_path
            .file_name()
            .and_then(|name| name.to_str())
            .ok_or_else(|| invalid_unit_state(id, "artifact file name is not UTF-8"))?;
        let temporary = parent.join(format!(".{file_name}.gps-tmp"));
        let mut lifecycle = ArtifactLifecycle::new(id);
        let mut file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary)
            .map_err(|error| io_error(&temporary, "create_artifact_temporary", error))?;
        file.write_all(bytes)
            .map_err(|error| io_error(&temporary, "write", error))?;
        file.flush()
            .map_err(|error| io_error(&temporary, "flush", error))?;
        drop(file);
        lifecycle.close()?;
        let identity = identify_file(&temporary, relative_path)?;
        verify_identity(&temporary, &identity)?;
        lifecycle.validate()?;
        fs::rename(&temporary, &final_path)
            .map_err(|error| io_error(&final_path, "finalize_artifact", error))?;
        lifecycle.catalog()?;
        let record = ArtifactRecord::from_identity(id, kind, &identity)?;
        self.artifacts.push(record.clone());
        Ok(record)
    }

    pub fn write_streamed_artifact(
        &mut self,
        id: &str,
        kind: &str,
        relative_path: &Path,
        write: impl FnOnce(&mut File) -> Result<(), Diagnostic>,
    ) -> Result<ArtifactRecord, Diagnostic> {
        if kind == "manifest_general" {
            return Err(invalid_unit_state(
                id,
                "a general manifest must use the semantic validation gate",
            ));
        }
        validate_relative_path(relative_path)?;
        let final_path = self.provisional.join(relative_path);
        let parent = final_path
            .parent()
            .ok_or_else(|| invalid_unit_state(id, "artifact path has no provisional parent"))?;
        fs::create_dir_all(parent).map_err(|error| io_error(parent, "create_dir_all", error))?;
        let file_name = final_path
            .file_name()
            .and_then(|name| name.to_str())
            .ok_or_else(|| invalid_unit_state(id, "artifact file name is not UTF-8"))?;
        let temporary = parent.join(format!(".{file_name}.gps-tmp"));
        let mut lifecycle = ArtifactLifecycle::new(id);
        let mut file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary)
            .map_err(|error| io_error(&temporary, "create_artifact_temporary", error))?;
        write(&mut file)?;
        file.flush()
            .map_err(|error| io_error(&temporary, "flush", error))?;
        drop(file);
        lifecycle.close()?;
        let identity = identify_file(&temporary, relative_path)?;
        verify_identity(&temporary, &identity)?;
        lifecycle.validate()?;
        fs::rename(&temporary, &final_path)
            .map_err(|error| io_error(&final_path, "finalize_artifact", error))?;
        lifecycle.catalog()?;
        let record = ArtifactRecord::from_identity(id, kind, &identity)?;
        self.artifacts.push(record.clone());
        Ok(record)
    }

    /// Copia una entrada exacta al temporal del artefacto, calculando de nuevo
    /// tamaño y SHA-256 antes de catalogarla. Una mutación entre admisión y
    /// copia se rechaza y nunca alcanza la ruta final del artefacto.
    pub fn copy_verified_artifact(
        &mut self,
        id: &str,
        kind: &str,
        source: &Path,
        relative_path: &Path,
        expected: &ArtifactIdentity,
    ) -> Result<ArtifactRecord, Diagnostic> {
        if kind == "manifest_general" {
            return Err(invalid_unit_state(
                id,
                "a general manifest must use the semantic validation gate",
            ));
        }
        validate_relative_path(relative_path)?;
        let final_path = self.provisional.join(relative_path);
        let parent = final_path
            .parent()
            .ok_or_else(|| invalid_unit_state(id, "artifact path has no provisional parent"))?;
        fs::create_dir_all(parent).map_err(|error| io_error(parent, "create_dir_all", error))?;
        let file_name = final_path
            .file_name()
            .and_then(|name| name.to_str())
            .ok_or_else(|| invalid_unit_state(id, "artifact file name is not UTF-8"))?;
        let temporary = parent.join(format!(".{file_name}.gps-tmp"));
        let mut input = File::open(source).map_err(|error| input_error(source, error))?;
        let mut output = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary)
            .map_err(|error| io_error(&temporary, "create_artifact_temporary", error))?;
        let mut digest = Sha256Writer::default();
        let mut buffer = [0_u8; 64 * 1024];
        loop {
            let read = input
                .read(&mut buffer)
                .map_err(|error| input_error(source, error))?;
            if read == 0 {
                break;
            }
            output
                .write_all(&buffer[..read])
                .map_err(|error| io_error(&temporary, "write", error))?;
            digest.update(&buffer[..read])?;
        }
        output
            .flush()
            .map_err(|error| io_error(&temporary, "flush", error))?;
        drop(output);

        let observed = digest.finish(relative_path);
        if observed.size_bytes() != expected.size_bytes() || observed.digest() != expected.digest()
        {
            return Err(Diagnostic::input(
                "GPS-INPUT-002",
                "Una entrada cambió durante su copia o no coincide con su identidad declarada.",
            )
            .with("path", source.display().to_string())
            .with("expected_digest", expected.digest())
            .with("observed_digest", observed.digest()));
        }
        verify_identity(&temporary, &observed)?;
        fs::rename(&temporary, &final_path)
            .map_err(|error| io_error(&final_path, "finalize_artifact", error))?;
        let record = ArtifactRecord::from_identity(id, kind, &observed)?;
        self.artifacts.push(record.clone());
        Ok(record)
    }

    pub fn create_truncated_temporary(
        &self,
        relative_path: &Path,
        bytes: &[u8],
    ) -> Result<PathBuf, Diagnostic> {
        validate_relative_path(relative_path)?;
        let final_path = self.provisional.join(relative_path);
        let parent = final_path
            .parent()
            .ok_or_else(|| invalid_unit_state("artifact", "temporary path has no parent"))?;
        if !parent.exists() {
            fs::create_dir_all(parent)
                .map_err(|error| io_error(parent, "create_dir_all", error))?;
        }
        let file_name = final_path
            .file_name()
            .and_then(|name| name.to_str())
            .ok_or_else(|| invalid_unit_state("artifact", "temporary file name is not UTF-8"))?;
        let path = parent.join(format!(".{file_name}.gps-tmp"));
        let mut file = File::create(&path).map_err(|error| io_error(&path, "create", error))?;
        file.write_all(bytes)
            .map_err(|error| io_error(&path, "write", error))?;
        Ok(path)
    }

    /// Publica mediante la única implementación de producción, que mantiene la
    /// última revalidación y el rename no-clobber dentro del mismo componente.
    ///
    /// Los consumidores no pueden sustituir el publicador ni insertar trabajo
    /// entre la última revalidación y el `rename`:
    ///
    /// ```compile_fail
    /// use genome_population_simulator::publication::AtomicPublisher;
    /// ```
    pub fn publish(self) -> Result<PathBuf, Diagnostic> {
        self.publish_with(&LinuxAtomicPublisher)
    }

    pub(crate) fn publish_with(
        self,
        publisher: &impl AtomicPublisher,
    ) -> Result<PathBuf, Diagnostic> {
        self.validate_publication_snapshot()?;
        let provisional_marker = self.provisional.join(INCOMPLETE_MARKER);
        preserve_or_restore_marker(&provisional_marker)?;
        if let Err(diagnostic) = publisher.before_final_validation(&self.provisional) {
            preserve_or_restore_marker(&provisional_marker).map_err(|restoration| {
                marker_invariant_error(
                    &provisional_marker,
                    "restore_marker_after_prepublication_hook_failure",
                    &diagnostic,
                    &restoration,
                )
            })?;
            return Err(diagnostic);
        }
        preserve_or_restore_marker(&provisional_marker)?;
        self.validate_publication_snapshot()?;
        match publisher.publish_no_replace(&self.provisional, &self.destination) {
            Ok(()) => {
                let published_marker = self.destination.join(INCOMPLETE_MARKER);
                preserve_or_restore_marker(&published_marker)?;
                match fs::remove_file(&published_marker) {
                    Ok(()) => Ok(self.destination),
                    Err(error) => {
                        let removal = io_error(&published_marker, "remove_published_marker", error);
                        preserve_or_restore_marker(&published_marker).map_err(|restoration| {
                            marker_invariant_error(
                                &published_marker,
                                "remove_published_marker",
                                &removal,
                                &restoration,
                            )
                        })?;
                        Err(removal)
                    }
                }
            }
            Err(diagnostic) => {
                preserve_or_restore_marker(&provisional_marker).map_err(|restoration| {
                    marker_invariant_error(
                        &provisional_marker,
                        "restore_marker_after_publish_failure",
                        &diagnostic,
                        &restoration,
                    )
                })?;
                Err(diagnostic)
            }
        }
    }

    fn validate_publication_snapshot(&self) -> Result<(), Diagnostic> {
        // La frontera física se comprueba antes de interpretar cualquier enlace
        // del manifiesto. Así, una mutación de árbol conserva su diagnóstico de
        // E/S y no queda enmascarada por la relectura semántica no-follow.
        self.validate_physical_tree(true)?;
        let manifests = self
            .artifacts
            .iter()
            .filter(|artifact| artifact.kind == "manifest_general")
            .collect::<Vec<_>>();
        let [manifest_record] = manifests.as_slice() else {
            return Err(invalid_unit_state(
                "execution",
                "exactly one validated general manifest is required before publication",
            ));
        };
        if manifest_record.id != "manifest-general"
            || manifest_record.path != "manifest.json"
            || manifest_record.state != UnitState::Complete
            || !manifest_record.validated
        {
            return Err(invalid_unit_state(
                "manifest-general",
                "the general manifest catalog record is not canonical and complete",
            ));
        }
        let manifest_bytes = read_workspace_artifact_nofollow(
            &self.provisional,
            Path::new("manifest.json"),
            manifest_record,
        )?;
        let manifest_path = self.provisional.join("manifest.json");
        validate_general_manifest_bytes(&manifest_path, &manifest_bytes)?;
        let (_, manifest_artifacts) =
            parse_general_manifest_bytes(&manifest_path, &manifest_bytes)?;
        let observed = self
            .artifacts
            .iter()
            .filter(|artifact| artifact.kind != "manifest_general")
            .cloned()
            .collect::<Vec<_>>();
        compare_artifact_catalogs(&manifest_artifacts, &observed)?;
        self.validate_physical_tree(true)
    }

    fn validate_manifest_workspace_artifacts(
        &self,
        manifest: &ValidatedGeneralManifest,
    ) -> Result<(), Diagnostic> {
        let observed = self
            .artifacts
            .iter()
            .filter(|artifact| artifact.kind != "manifest_general")
            .cloned()
            .collect::<Vec<_>>();
        compare_artifact_catalogs(manifest.artifacts(), &observed)?;
        Ok(())
    }

    fn validate_physical_tree(&self, manifest_expected: bool) -> Result<(), Diagnostic> {
        let mut expected_files = BTreeMap::new();
        let mut expected_directories = BTreeSet::new();
        let mut manifest_records = 0_usize;
        for artifact in &self.artifacts {
            validate_relative_path(Path::new(&artifact.path))?;
            if artifact.kind == "manifest_general" {
                manifest_records += 1;
            }
            if expected_files
                .insert(artifact.path.clone(), Some(artifact))
                .is_some()
            {
                return Err(invalid_unit_state(
                    &artifact.id,
                    "two artifact records resolve to the same canonical path",
                ));
            }
            register_parent_directories(Path::new(&artifact.path), &mut expected_directories)?;
        }
        if manifest_records != usize::from(manifest_expected) {
            return Err(invalid_unit_state(
                "manifest-general",
                "the workspace has an unexpected number of general-manifest records for this stage",
            ));
        }
        if expected_files
            .insert(INCOMPLETE_MARKER.to_owned(), None)
            .is_some()
        {
            return Err(invalid_unit_state(
                INCOMPLETE_MARKER,
                "the incomplete marker path is reserved and cannot be an artifact",
            ));
        }

        let root = openat(
            CWD,
            &self.provisional,
            OFlags::RDONLY | OFlags::DIRECTORY | OFlags::NOFOLLOW | OFlags::CLOEXEC,
            Mode::empty(),
        )
        .map_err(|error| {
            rustix_io_error(&self.provisional, "open_workspace_root_nofollow", error)
        })?;
        let mut observed_files = BTreeSet::new();
        let mut observed_directories = BTreeSet::new();
        let mut observed_file_objects = BTreeSet::new();
        inspect_directory_nofollow(
            &root,
            Path::new(""),
            &expected_files,
            &expected_directories,
            &mut observed_files,
            &mut observed_directories,
            &mut observed_file_objects,
        )?;
        let expected_file_paths = expected_files.keys().cloned().collect::<BTreeSet<_>>();
        if observed_files != expected_file_paths || observed_directories != expected_directories {
            return Err(invalid_unit_state(
                "execution",
                "the physical provisional tree differs from the exact cataloged layout",
            ));
        }
        Ok(())
    }
}

fn register_parent_directories(
    path: &Path,
    directories: &mut BTreeSet<String>,
) -> Result<(), Diagnostic> {
    let mut parent = path.parent();
    while let Some(path) = parent {
        if path.as_os_str().is_empty() {
            break;
        }
        validate_relative_path(path)?;
        directories.insert(path.to_string_lossy().into_owned());
        parent = path.parent();
    }
    Ok(())
}

fn inspect_directory_nofollow(
    directory_fd: &OwnedFd,
    relative_directory: &Path,
    expected_files: &BTreeMap<String, Option<&ArtifactRecord>>,
    expected_directories: &BTreeSet<String>,
    observed_files: &mut BTreeSet<String>,
    observed_directories: &mut BTreeSet<String>,
    observed_file_objects: &mut BTreeSet<(u64, u64)>,
) -> Result<(), Diagnostic> {
    let mut directory = Dir::read_from(directory_fd).map_err(|error| {
        rustix_io_error(
            relative_directory,
            "read_workspace_directory_nofollow",
            error,
        )
    })?;
    for entry in &mut directory {
        let entry = entry
            .map_err(|error| rustix_io_error(relative_directory, "read_workspace_entry", error))?;
        let name = entry.file_name();
        if matches!(name.to_bytes(), b"." | b"..") {
            continue;
        }
        let name_path = Path::new(std::ffi::OsStr::from_bytes(name.to_bytes()));
        let relative = relative_directory.join(name_path);
        validate_relative_path(&relative)?;
        let relative_text = relative
            .to_str()
            .ok_or_else(|| invalid_unit_state("execution", "workspace path is not UTF-8"))?
            .to_owned();
        let stat = statat(directory_fd, name, AtFlags::SYMLINK_NOFOLLOW)
            .map_err(|error| rustix_io_error(&relative, "stat_workspace_entry_nofollow", error))?;
        match FileType::from_raw_mode(stat.st_mode) {
            file_type if file_type.is_dir() => {
                if !expected_directories.contains(&relative_text)
                    || !observed_directories.insert(relative_text)
                {
                    return Err(invalid_unit_state(
                        &relative.display().to_string(),
                        "unexpected or duplicate directory in provisional tree",
                    ));
                }
                let child = openat(
                    directory_fd,
                    name,
                    OFlags::RDONLY | OFlags::DIRECTORY | OFlags::NOFOLLOW | OFlags::CLOEXEC,
                    Mode::empty(),
                )
                .map_err(|error| {
                    rustix_io_error(&relative, "open_workspace_directory_nofollow", error)
                })?;
                if !FileType::from_raw_mode(
                    fstat(&child)
                        .map_err(|error| rustix_io_error(&relative, "fstat_directory", error))?
                        .st_mode,
                )
                .is_dir()
                {
                    return Err(invalid_unit_state(
                        &relative.display().to_string(),
                        "workspace directory changed type while being inspected",
                    ));
                }
                inspect_directory_nofollow(
                    &child,
                    &relative,
                    expected_files,
                    expected_directories,
                    observed_files,
                    observed_directories,
                    observed_file_objects,
                )?;
            }
            file_type if file_type.is_file() => {
                let expected = expected_files.get(&relative_text).ok_or_else(|| {
                    invalid_unit_state(
                        &relative.display().to_string(),
                        "uncataloged physical file in provisional tree",
                    )
                })?;
                if !observed_files.insert(relative_text) {
                    return Err(invalid_unit_state(
                        &relative.display().to_string(),
                        "duplicate physical path in provisional tree",
                    ));
                }
                let file_fd = openat(
                    directory_fd,
                    name,
                    OFlags::RDONLY | OFlags::NOFOLLOW | OFlags::CLOEXEC,
                    Mode::empty(),
                )
                .map_err(|error| {
                    rustix_io_error(&relative, "open_workspace_file_nofollow", error)
                })?;
                let opened_stat = fstat(&file_fd)
                    .map_err(|error| rustix_io_error(&relative, "fstat_file", error))?;
                if !FileType::from_raw_mode(opened_stat.st_mode).is_file() {
                    return Err(invalid_unit_state(
                        &relative.display().to_string(),
                        "workspace file changed type while being inspected",
                    ));
                }
                if opened_stat.st_nlink != 1 {
                    return Err(invalid_unit_state(
                        &relative.display().to_string(),
                        "every publishable file must have exactly one physical hard link",
                    ));
                }
                if !observed_file_objects.insert((opened_stat.st_dev, opened_stat.st_ino)) {
                    return Err(invalid_unit_state(
                        &relative.display().to_string(),
                        "two canonical paths alias the same physical file",
                    ));
                }
                if let Some(record) = expected {
                    verify_open_artifact(file_fd, &relative, record)?;
                }
            }
            _ => {
                return Err(invalid_unit_state(
                    &relative.display().to_string(),
                    "symlinks and special files are forbidden in the provisional tree",
                ));
            }
        }
    }
    Ok(())
}

fn verify_open_artifact(
    file_fd: OwnedFd,
    relative: &Path,
    artifact: &ArtifactRecord,
) -> Result<(), Diagnostic> {
    let mut file = File::from(file_fd);
    let mut digest = Sha256Writer::default();
    let mut buffer = [0_u8; 64 * 1024];
    loop {
        let read = file
            .read(&mut buffer)
            .map_err(|error| io_error(relative, "read_workspace_file_nofollow", error))?;
        if read == 0 {
            break;
        }
        digest.update(&buffer[..read])?;
    }
    let observed = digest.finish(relative);
    if observed.size_bytes() != artifact.size_bytes || observed.digest() != artifact.digest {
        return Err(invalid_unit_state(
            &artifact.id,
            "the no-follow physical file differs in size or SHA-256 from its record",
        ));
    }
    Ok(())
}

fn read_workspace_artifact_nofollow(
    root: &Path,
    relative: &Path,
    artifact: &ArtifactRecord,
) -> Result<Vec<u8>, Diagnostic> {
    validate_relative_path(relative)?;
    let components = relative
        .to_str()
        .ok_or_else(|| invalid_unit_state(&artifact.id, "artifact path is not UTF-8"))?
        .split('/')
        .collect::<Vec<_>>();
    let mut directory = openat(
        CWD,
        root,
        OFlags::RDONLY | OFlags::DIRECTORY | OFlags::NOFOLLOW | OFlags::CLOEXEC,
        Mode::empty(),
    )
    .map_err(|error| rustix_io_error(root, "open_workspace_root_nofollow", error))?;
    for component in &components[..components.len() - 1] {
        directory = openat(
            &directory,
            *component,
            OFlags::RDONLY | OFlags::DIRECTORY | OFlags::NOFOLLOW | OFlags::CLOEXEC,
            Mode::empty(),
        )
        .map_err(|error| rustix_io_error(relative, "open_workspace_parent_nofollow", error))?;
    }
    let file_fd = openat(
        &directory,
        components[components.len() - 1],
        OFlags::RDONLY | OFlags::NOFOLLOW | OFlags::CLOEXEC,
        Mode::empty(),
    )
    .map_err(|error| rustix_io_error(relative, "open_workspace_file_nofollow", error))?;
    let stat = fstat(&file_fd)
        .map_err(|error| rustix_io_error(relative, "fstat_workspace_file", error))?;
    if !FileType::from_raw_mode(stat.st_mode).is_file() {
        return Err(invalid_unit_state(
            &artifact.id,
            "the no-follow manifest path is not a regular file",
        ));
    }
    if stat.st_nlink != 1 {
        return Err(invalid_unit_state(
            &artifact.id,
            "the no-follow manifest file must have exactly one physical hard link",
        ));
    }
    let mut file = File::from(file_fd);
    let mut bytes = Vec::new();
    file.read_to_end(&mut bytes)
        .map_err(|error| io_error(relative, "read_workspace_file_nofollow", error))?;
    let mut digest = Sha256Writer::default();
    digest.update(&bytes)?;
    let observed = digest.finish(relative);
    if observed.size_bytes() != artifact.size_bytes || observed.digest() != artifact.digest {
        return Err(invalid_unit_state(
            &artifact.id,
            "the no-follow file bytes differ in size or SHA-256 from their record",
        ));
    }
    Ok(bytes)
}

fn compare_artifact_catalogs(
    expected: &[ArtifactRecord],
    observed: &[ArtifactRecord],
) -> Result<(), Diagnostic> {
    let mut expected_by_id = BTreeMap::new();
    let mut expected_paths = BTreeSet::new();
    for artifact in expected {
        if artifact.kind == "manifest_general" {
            return Err(invalid_unit_state(
                &artifact.id,
                "the manifest cannot catalog itself as a payload artifact",
            ));
        }
        validate_relative_path(Path::new(&artifact.path))?;
        if expected_by_id.insert(&artifact.id, artifact).is_some()
            || !expected_paths.insert(&artifact.path)
        {
            return Err(invalid_unit_state(
                &artifact.id,
                "manifest artifact IDs and paths must be unique",
            ));
        }
    }

    let mut observed_by_id = BTreeMap::new();
    let mut observed_paths = BTreeSet::new();
    for artifact in observed {
        validate_relative_path(Path::new(&artifact.path))?;
        if observed_by_id.insert(&artifact.id, artifact).is_some()
            || !observed_paths.insert(&artifact.path)
        {
            return Err(invalid_unit_state(
                &artifact.id,
                "workspace artifact IDs and paths must be unique",
            ));
        }
    }
    if expected_by_id.len() != observed_by_id.len() {
        return Err(invalid_unit_state(
            "manifest-general",
            "manifest and workspace contain a different number of payload artifacts",
        ));
    }
    for (id, expected) in expected_by_id {
        let observed = observed_by_id.get(id).ok_or_else(|| {
            invalid_unit_state(id, "a manifest artifact is absent from this workspace")
        })?;
        if *observed != expected {
            return Err(invalid_unit_state(
                id,
                "manifest and workspace disagree on artifact class, path, state, validation, size or SHA-256",
            ));
        }
    }
    Ok(())
}

const fn stage_rank(stage: ArtifactStage) -> u8 {
    match stage {
        ArtifactStage::Open | ArtifactStage::Failed => 0,
        ArtifactStage::Closed => 1,
        ArtifactStage::Validated => 2,
        ArtifactStage::Cataloged => 3,
    }
}

fn require_regular_marker(marker: &Path) -> Result<(), Diagnostic> {
    let metadata =
        fs::metadata(marker).map_err(|error| io_error(marker, "inspect_marker", error))?;
    if metadata.is_file() {
        return Ok(());
    }
    Err(invalid_unit_state(
        "execution",
        format!("{} is not a regular INCOMPLETE marker", marker.display()),
    ))
}

fn preserve_or_restore_marker(marker: &Path) -> Result<(), Diagnostic> {
    match fs::metadata(marker) {
        Ok(metadata) if metadata.is_file() => Ok(()),
        Ok(_) => Err(invalid_unit_state(
            "execution",
            format!("{} is not a regular INCOMPLETE marker", marker.display()),
        )),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
            let mut file = OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(marker)
                .map_err(|error| io_error(marker, "restore_marker", error))?;
            file.write_all(b"incomplete\n")
                .map_err(|error| io_error(marker, "restore_marker_write", error))?;
            file.flush()
                .map_err(|error| io_error(marker, "restore_marker_flush", error))?;
            drop(file);
            require_regular_marker(marker)
        }
        Err(error) => Err(io_error(marker, "inspect_marker", error)),
    }
}

fn marker_invariant_error(
    marker: &Path,
    operation: &str,
    primary: &Diagnostic,
    restoration: &Diagnostic,
) -> Diagnostic {
    Diagnostic::io(
        "GPS-IO-005",
        "Una unidad incompleta o inválida no puede publicarse como completa.",
    )
    .with("unit", "execution")
    .with("state", "incomplete_marker_unrecoverable")
    .with("path", marker.display().to_string())
    .with("operation", operation)
    .with("primary_diagnostic", primary.code())
    .with("restoration_diagnostic", restoration.code())
}

fn validate_relative_path(path: &Path) -> Result<(), Diagnostic> {
    let text = path.to_str().ok_or_else(|| {
        invalid_unit_state(
            &path.display().to_string(),
            "artifact path must use canonical UTF-8",
        )
    })?;
    let canonical_components = text.split('/').collect::<Vec<_>>();
    if text.is_empty()
        || path.is_absolute()
        || text.contains(['\\', '\0'])
        || text.chars().any(char::is_control)
        || canonical_components
            .iter()
            .any(|component| component.is_empty() || matches!(*component, "." | ".."))
        || canonical_components.join("/") != text
        || path
            .components()
            .any(|component| !matches!(component, Component::Normal(_)))
    {
        return Err(invalid_unit_state(
            &path.display().to_string(),
            "artifact path must be one unique canonical relative UTF-8 representation",
        ));
    }
    Ok(())
}

fn invalid_unit_state(unit: &str, detail: impl Into<String>) -> Diagnostic {
    Diagnostic::io(
        "GPS-IO-005",
        "Una unidad incompleta o inválida no puede publicarse como completa.",
    )
    .with("unit", unit)
    .with("state", "incomplete")
    .with("detail", detail.into())
}

fn unsupported_publication(path: &Path, filesystem: impl Into<String>) -> Diagnostic {
    Diagnostic::io(
        "GPS-IO-004",
        "El sistema de ficheros no ofrece publicación atómica no-clobber.",
    )
    .with("path", path.display().to_string())
    .with("filesystem", filesystem.into())
    .with("operation", "renameat2(RENAME_NOREPLACE)")
}

fn io_error(path: &Path, operation: &str, error: std::io::Error) -> Diagnostic {
    Diagnostic::io("GPS-IO-001", "Una operación de lectura o escritura falló.")
        .with("path", path.display().to_string())
        .with("operation", operation)
        .with("detail", error.to_string())
}

fn rustix_io_error(path: &Path, operation: &str, error: rustix::io::Errno) -> Diagnostic {
    Diagnostic::io("GPS-IO-001", "Una operación de lectura o escritura falló.")
        .with("path", path.display().to_string())
        .with("operation", operation)
        .with("detail", error.to_string())
}

fn input_error(path: &Path, error: std::io::Error) -> Diagnostic {
    Diagnostic::input(
        "GPS-INPUT-001",
        "No se puede abrir o leer una entrada declarada.",
    )
    .with("path", path.display().to_string())
    .with("detail", error.to_string())
}

#[cfg(test)]
mod tests {
    use std::cell::Cell;
    use std::fs;
    use std::os::unix::fs::{MetadataExt, symlink};
    use std::path::{Path, PathBuf};
    use std::sync::atomic::{AtomicU64, Ordering};

    use super::{
        AtomicPublisher, ExecutionWorkspace, UnsupportedAtomicPublisher, compare_artifact_catalogs,
        validate_relative_path,
    };
    use crate::application::load_and_validate_config;
    use crate::diagnostics::Diagnostic;
    use crate::execution::execute_validated_config;
    use crate::results::{ArtifactRecord, UnitState};

    static NONCE: AtomicU64 = AtomicU64::new(0);

    fn temporary(label: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "gps-publication-unit-{}-{}-{label}",
            std::process::id(),
            NONCE.fetch_add(1, Ordering::Relaxed)
        ))
    }

    fn fixture(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/config")
            .join(name)
    }

    fn matching_workspace(destination: &Path) -> (ExecutionWorkspace, PathBuf) {
        let source = temporary("source-result");
        let effective = load_and_validate_config(&fixture("general-seeded.toml"))
            .expect("fixture configuration");
        execute_validated_config(&effective, &source).expect("fixture execution");
        let manifest_bytes = fs::read(source.join("manifest.json")).expect("fixture manifest");
        let manifest: serde_json::Value =
            serde_json::from_slice(&manifest_bytes).expect("fixture JSON");
        let mut workspace = ExecutionWorkspace::create(destination).expect("workspace");
        let mut first_payload = None;
        for artifact in manifest["artifacts"].as_array().expect("artifact catalog") {
            let id = artifact["id"].as_str().expect("artifact ID");
            let kind = artifact["kind"].as_str().expect("artifact kind");
            let relative = PathBuf::from(artifact["path"].as_str().expect("artifact path"));
            let bytes = fs::read(source.join(&relative)).expect("artifact bytes");
            workspace
                .write_artifact(id, kind, &relative, &bytes)
                .expect("catalog fixture artifact");
            first_payload.get_or_insert(relative);
        }
        let capability = workspace
            .validate_general_manifest(Path::new("manifest.json"), &manifest_bytes)
            .expect("bind fixture manifest");
        workspace
            .write_general_manifest(&capability)
            .expect("write fixture manifest");
        fs::remove_dir_all(source).expect("source fixture cleanup");
        (workspace, first_payload.expect("at least one payload"))
    }

    fn record() -> ArtifactRecord {
        ArtifactRecord {
            id: "artifact-a".to_owned(),
            kind: "reference_fasta".to_owned(),
            path: "resources/a.fa".to_owned(),
            size_bytes: 4,
            digest: "a".repeat(64),
            state: UnitState::Complete,
            validated: true,
        }
    }

    #[test]
    fn workspace_catalog_requires_every_artifact_field_exactly() {
        let expected = [record()];
        compare_artifact_catalogs(&expected, &expected).unwrap();

        let mutations: [fn(&mut ArtifactRecord); 7] = [
            |artifact| artifact.id = "artifact-b".to_owned(),
            |artifact| artifact.kind = "reference_reduction_map".to_owned(),
            |artifact| artifact.path = "resources/b.fa".to_owned(),
            |artifact| artifact.size_bytes += 1,
            |artifact| artifact.digest = "b".repeat(64),
            |artifact| artifact.state = UnitState::Failed,
            |artifact| artifact.validated = false,
        ];
        for mutate in mutations {
            let mut observed = record();
            mutate(&mut observed);
            assert!(compare_artifact_catalogs(&expected, &[observed]).is_err());
        }
        assert!(compare_artifact_catalogs(&expected, &[]).is_err());
        assert!(compare_artifact_catalogs(&[], &expected).is_err());

        let mut second = record();
        second.id = "artifact-b".to_owned();
        assert!(compare_artifact_catalogs(&[record(), second], &expected).is_err());
    }

    #[derive(Clone, Copy)]
    enum Mutation {
        Extra,
        FileSymlink,
        Content,
    }

    struct MutatingPublisher {
        mutation: Mutation,
        payload: PathBuf,
        external: PathBuf,
        rename_called: Cell<bool>,
    }

    impl AtomicPublisher for MutatingPublisher {
        fn before_final_validation(&self, provisional: &Path) -> Result<(), Diagnostic> {
            match self.mutation {
                Mutation::Extra => fs::write(provisional.join("uncataloged.bin"), b"extra")
                    .expect("inject extra file"),
                Mutation::FileSymlink => {
                    let path = provisional.join(&self.payload);
                    fs::copy(&path, &self.external).expect("copy external alias target");
                    fs::remove_file(&path).expect("remove cataloged file");
                    symlink(&self.external, &path).expect("inject symlink");
                }
                Mutation::Content => fs::write(provisional.join(&self.payload), b"changed")
                    .expect("inject content mutation"),
            }
            Ok(())
        }

        fn publish_no_replace(
            &self,
            _provisional: &Path,
            _destination: &Path,
        ) -> Result<(), Diagnostic> {
            self.rename_called.set(true);
            panic!("rename must not be reached after a controlled mutation")
        }
    }

    #[test]
    fn internal_hooks_reject_extra_file_symlink_and_content_before_rename() {
        for (label, mutation) in [
            ("extra", Mutation::Extra),
            ("file-symlink", Mutation::FileSymlink),
            ("content", Mutation::Content),
        ] {
            let root = temporary(label);
            fs::create_dir_all(&root).expect("test root");
            let destination = root.join("published");
            let (workspace, payload) = matching_workspace(&destination);
            let provisional = workspace.root().to_path_buf();
            let publisher = MutatingPublisher {
                mutation,
                payload,
                external: root.join("external.bin"),
                rename_called: Cell::new(false),
            };
            let error = workspace
                .publish_with(&publisher)
                .expect_err("final revalidation must reject mutation");
            assert_eq!(error.code(), "GPS-IO-005");
            assert!(!publisher.rename_called.get());
            assert!(provisional.join("INCOMPLETE").is_file());
            assert!(!destination.exists());
            fs::remove_dir_all(root).expect("test cleanup");
        }
    }

    struct ParentSymlinkPublisher {
        payload: PathBuf,
        external_directory: PathBuf,
        rename_called: Cell<bool>,
    }

    impl AtomicPublisher for ParentSymlinkPublisher {
        fn before_final_validation(&self, provisional: &Path) -> Result<(), Diagnostic> {
            let payload = provisional.join(&self.payload);
            let parent = payload
                .parent()
                .expect("fixture payload must have an intermediate directory");
            assert_ne!(parent, provisional);
            fs::rename(parent, &self.external_directory)
                .expect("move real intermediate directory outside the workspace");
            symlink(&self.external_directory, parent)
                .expect("replace intermediate directory with a symlink");
            assert!(
                payload.is_file(),
                "cataloged bytes remain reachable through the injected parent symlink"
            );
            Ok(())
        }

        fn publish_no_replace(
            &self,
            _provisional: &Path,
            _destination: &Path,
        ) -> Result<(), Diagnostic> {
            self.rename_called.set(true);
            panic!("publisher must not run after an intermediate-directory symlink")
        }
    }

    #[test]
    fn intermediate_directory_symlink_is_rejected_before_publisher_and_preserves_incomplete() {
        let root = temporary("parent-symlink");
        fs::create_dir_all(&root).expect("test root");
        let destination = root.join("published");
        let (workspace, payload) = matching_workspace(&destination);
        let provisional = workspace.root().to_path_buf();
        let publisher = ParentSymlinkPublisher {
            payload,
            external_directory: root.join("external-parent"),
            rename_called: Cell::new(false),
        };

        let error = workspace
            .publish_with(&publisher)
            .expect_err("a symbolic intermediate directory must fail final validation");
        assert_eq!(error.code(), "GPS-IO-005");
        assert!(!publisher.rename_called.get());
        assert!(provisional.join("INCOMPLETE").is_file());
        assert!(!destination.exists());
        fs::remove_dir_all(root).expect("test cleanup");
    }

    #[test]
    fn canonical_relative_path_rejects_every_textual_alias_directly() {
        for valid in ["a", "a/b", "resources/references/reference.fa"] {
            validate_relative_path(Path::new(valid)).expect("canonical relative path");
        }
        for alias in ["a//b", "./a/b", "a/./b", "a/../b", "a/b/", "a\\b"] {
            let error = validate_relative_path(Path::new(alias))
                .expect_err("textual path alias must be rejected directly");
            assert_eq!(error.code(), "GPS-IO-005", "alias={alias}");
            assert_eq!(
                error.context().get("detail").map(String::as_str),
                Some("artifact path must be one unique canonical relative UTF-8 representation"),
                "alias={alias} must fail for canonical form, not content integrity"
            );
        }
    }

    #[test]
    fn external_hardlink_with_valid_bytes_is_rejected_before_rename() {
        let root = temporary("hardlink");
        fs::create_dir_all(&root).expect("test root");
        let destination = root.join("published");
        let (workspace, payload) = matching_workspace(&destination);
        let provisional = workspace.root().to_path_buf();
        let external = root.join("external-hardlink");
        fs::hard_link(provisional.join(payload), &external).expect("external hardlink");
        let error = workspace
            .publish_with(&UnsupportedAtomicPublisher)
            .expect_err("multiple physical links must fail before publisher invocation");
        assert_eq!(error.code(), "GPS-IO-005");
        assert!(provisional.join("INCOMPLETE").is_file());
        assert!(!destination.exists());
        fs::remove_dir_all(root).expect("test cleanup");
    }

    #[test]
    fn regular_files_with_one_link_publish_even_when_directories_have_multiple_links() {
        let root = temporary("positive-links");
        fs::create_dir_all(&root).expect("test root");
        let destination = root.join("published");
        let (workspace, _) = matching_workspace(&destination);
        let provisional = workspace.root().to_path_buf();
        assert!(
            fs::metadata(&provisional)
                .expect("provisional metadata")
                .nlink()
                > 1,
            "directory link counts are not constrained by the file invariant"
        );
        workspace
            .publish()
            .expect("legitimate workspace publication");
        assert!(destination.join("manifest.json").is_file());
        assert!(!destination.join("INCOMPLETE").exists());
        fs::remove_dir_all(root).expect("test cleanup");
    }
}
