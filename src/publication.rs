//! Ensamblado provisional y publicación Linux atómica sin sobrescritura.

use std::collections::BTreeMap;
use std::fs::{self, File, OpenOptions};
use std::io::{Read, Write};
use std::path::{Component, Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use rustix::fs::{CWD, RenameFlags, renameat_with};

use crate::diagnostics::Diagnostic;
use crate::integrity::{ArtifactIdentity, Sha256Writer, identify_file, verify_identity};
use crate::results::{ArtifactRecord, UnitRecord, UnitState};

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

pub trait AtomicPublisher {
    /// Publica con una única operación atómica no-clobber. Si devuelve error,
    /// debe dejar el provisional en su ruta original; nunca puede retirar ni
    /// interpretar el marcador `INCOMPLETE`.
    fn publish_no_replace(&self, provisional: &Path, destination: &Path) -> Result<(), Diagnostic>;
}

#[derive(Debug, Default)]
pub struct LinuxAtomicPublisher;

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
pub struct UnsupportedAtomicPublisher;

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
    pub fn root(&self) -> &Path {
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

    pub fn write_artifact(
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

    pub fn publish(self, publisher: &impl AtomicPublisher) -> Result<PathBuf, Diagnostic> {
        if contains_artifact_temporary(&self.provisional)? {
            return Err(invalid_unit_state(
                "execution",
                "an incomplete per-artifact temporary remains in the provisional directory",
            ));
        }
        if !self.artifacts.iter().any(|artifact| {
            artifact.kind == "manifest_general"
                && artifact.path == "manifest.json"
                && artifact.state == UnitState::Complete
                && artifact.validated
        }) {
            return Err(invalid_unit_state(
                "execution",
                "a validated general manifest is required before publication",
            ));
        }
        let provisional_marker = self.provisional.join(INCOMPLETE_MARKER);
        preserve_or_restore_marker(&provisional_marker)?;
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
}

const fn stage_rank(stage: ArtifactStage) -> u8 {
    match stage {
        ArtifactStage::Open | ArtifactStage::Failed => 0,
        ArtifactStage::Closed => 1,
        ArtifactStage::Validated => 2,
        ArtifactStage::Cataloged => 3,
    }
}

fn contains_artifact_temporary(directory: &Path) -> Result<bool, Diagnostic> {
    for entry in fs::read_dir(directory).map_err(|error| io_error(directory, "read_dir", error))? {
        let entry = entry.map_err(|error| io_error(directory, "read_dir_entry", error))?;
        let file_type = entry
            .file_type()
            .map_err(|error| io_error(&entry.path(), "file_type", error))?;
        if file_type.is_dir() {
            if contains_artifact_temporary(&entry.path())? {
                return Ok(true);
            }
        } else if entry
            .file_name()
            .to_str()
            .is_some_and(|name| name.ends_with(".gps-tmp"))
        {
            return Ok(true);
        }
    }
    Ok(false)
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
    if path.as_os_str().is_empty()
        || path.is_absolute()
        || path
            .components()
            .any(|component| !matches!(component, Component::Normal(_)))
    {
        return Err(invalid_unit_state(
            &path.display().to_string(),
            "artifact path must be a non-empty relative path below the provisional root",
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

fn input_error(path: &Path, error: std::io::Error) -> Diagnostic {
    Diagnostic::input(
        "GPS-INPUT-001",
        "No se puede abrir o leer una entrada declarada.",
    )
    .with("path", path.display().to_string())
    .with("detail", error.to_string())
}
