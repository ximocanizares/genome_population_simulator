//! Provisional privado y publicación Linux atómica sin sobrescritura.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, DirBuilder, File, OpenOptions};
use std::io::{Read, Write};
use std::os::unix::fs::{DirBuilderExt, PermissionsExt};
use std::path::{Component, Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use rustix::fs::{CWD, RenameFlags, renameat_with};

use crate::diagnostics::Diagnostic;
use crate::integrity::{ArtifactIdentity, Sha256Writer, identify_file, verify_identity};
use crate::results::{
    ArtifactRecord, parse_general_manifest_bytes, validate_general_manifest_bytes,
    validate_published_resources, validate_relative_link,
};

const INCOMPLETE_MARKER: &str = "INCOMPLETE";
static WORKSPACE_COUNTER: AtomicU64 = AtomicU64::new(0);
static TEMPORARY_COUNTER: AtomicU64 = AtomicU64::new(0);

#[derive(Debug)]
pub struct ExecutionWorkspace {
    destination: PathBuf,
    /// Wrapper privado que conserva `INCOMPLETE` fuera del árbol publicable.
    provisional: PathBuf,
    /// Árbol que se renombra atómicamente al destino y nunca contiene marcador.
    payload: PathBuf,
    artifacts: Vec<ArtifactRecord>,
    manifest_written: bool,
}

impl ExecutionWorkspace {
    pub fn create(destination: &Path) -> Result<Self, Diagnostic> {
        let parent = destination.parent().unwrap_or_else(|| Path::new("."));
        let name = destination
            .file_name()
            .and_then(|name| name.to_str())
            .filter(|name| !name.is_empty())
            .ok_or_else(|| invalid_state("execution", "destination has no file name"))?;
        fs::create_dir_all(parent).map_err(|error| io_error(parent, "create_dir_all", error))?;
        if destination.exists() {
            return Err(destination_exists(destination));
        }
        for _ in 0..128 {
            let counter = WORKSPACE_COUNTER.fetch_add(1, Ordering::Relaxed);
            let provisional = parent.join(format!(
                ".{name}.gps-incomplete-{}-{counter}",
                std::process::id()
            ));
            let mut builder = DirBuilder::new();
            builder.mode(0o700);
            match builder.create(&provisional) {
                Ok(()) => {
                    let initialized = (|| {
                        fs::write(provisional.join(INCOMPLETE_MARKER), b"incomplete\n")
                            .map_err(|error| io_error(&provisional, "write_marker", error))?;
                        let payload = provisional.join("payload");
                        let mut payload_builder = DirBuilder::new();
                        payload_builder.mode(0o700);
                        payload_builder
                            .create(&payload)
                            .map_err(|error| io_error(&payload, "create_payload", error))?;
                        fs::set_permissions(&provisional, fs::Permissions::from_mode(0o700))
                            .map_err(|error| {
                                io_error(&provisional, "set_private_permissions", error)
                            })
                    })();
                    if let Err(diagnostic) = initialized {
                        let _ = fs::remove_dir_all(&provisional);
                        return Err(diagnostic);
                    }
                    return Ok(Self {
                        destination: destination.to_path_buf(),
                        payload: provisional.join("payload"),
                        provisional,
                        artifacts: Vec::new(),
                        manifest_written: false,
                    });
                }
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(error) => return Err(io_error(&provisional, "create_provisional", error)),
            }
        }
        Err(invalid_state(
            "execution",
            "could not allocate a unique provisional sibling",
        ))
    }

    #[must_use]
    pub(crate) fn root(&self) -> &Path {
        &self.payload
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
        self.write_streamed_artifact(id, kind, relative_path, |file| {
            file.write_all(bytes)
                .map_err(|error| io_error(relative_path, "write", error))
        })
    }

    pub fn write_streamed_artifact(
        &mut self,
        id: &str,
        kind: &str,
        relative_path: &Path,
        write: impl FnOnce(&mut File) -> Result<(), Diagnostic>,
    ) -> Result<ArtifactRecord, Diagnostic> {
        self.check_new_artifact(id, relative_path)?;
        let final_path = self.payload.join(relative_path);
        let parent = final_path
            .parent()
            .ok_or_else(|| invalid_state(id, "artifact path has no provisional parent"))?;
        fs::create_dir_all(parent).map_err(|error| io_error(parent, "create_dir_all", error))?;
        let temporary = temporary_path(parent, &final_path, id)?;
        let result = (|| {
            let mut file = OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&temporary)
                .map_err(|error| io_error(&temporary, "create_artifact_temporary", error))?;
            write(&mut file)?;
            file.flush()
                .map_err(|error| io_error(&temporary, "flush", error))?;
            drop(file);
            let identity = identify_file(&temporary, relative_path)?;
            verify_identity(&temporary, &identity)?;
            fs::rename(&temporary, &final_path)
                .map_err(|error| io_error(&final_path, "finalize_artifact", error))?;
            ArtifactRecord::from_identity(id, kind, &identity)
        })();
        if result.is_err() {
            let _ = fs::remove_file(&temporary);
        }
        let record = result?;
        self.artifacts.push(record.clone());
        Ok(record)
    }

    /// Copia una entrada admitida y exige que sus bytes sigan siendo los mismos.
    pub fn copy_verified_artifact(
        &mut self,
        id: &str,
        kind: &str,
        source: &Path,
        relative_path: &Path,
        expected: &ArtifactIdentity,
    ) -> Result<ArtifactRecord, Diagnostic> {
        self.check_new_artifact(id, relative_path)?;
        let final_path = self.payload.join(relative_path);
        let parent = final_path
            .parent()
            .ok_or_else(|| invalid_state(id, "artifact path has no provisional parent"))?;
        fs::create_dir_all(parent).map_err(|error| io_error(parent, "create_dir_all", error))?;
        let temporary = temporary_path(parent, &final_path, id)?;
        let result = (|| {
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
            if observed.size_bytes() != expected.size_bytes()
                || observed.digest() != expected.digest()
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
            ArtifactRecord::from_identity(id, kind, &observed)
        })();
        if result.is_err() {
            let _ = fs::remove_file(&temporary);
        }
        let record = result?;
        self.artifacts.push(record.clone());
        Ok(record)
    }

    /// Valida, escribe y relee el manifiesto general. `manifest.json` no entra
    /// en el catálogo porque su auto-hash no tendría solución finita.
    pub fn write_general_manifest(&mut self, bytes: &[u8]) -> Result<(), Diagnostic> {
        if self.manifest_written {
            return Err(invalid_state(
                "manifest-general",
                "the general manifest is already finalized",
            ));
        }
        let manifest_path = self.payload.join("manifest.json");
        let declared = parse_general_manifest_bytes(&manifest_path, bytes)?;
        compare_catalogs(&declared, &self.artifacts)?;
        self.verify_catalog()?;
        self.validate_effective_config(&manifest_path, bytes)?;
        let temporary = self.payload.join(format!(
            ".manifest.json.gps-tmp-{}-{}",
            std::process::id(),
            TEMPORARY_COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        let result = (|| {
            let mut file = OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&temporary)
                .map_err(|error| io_error(&temporary, "create_manifest_temporary", error))?;
            file.write_all(bytes)
                .map_err(|error| io_error(&temporary, "write", error))?;
            file.flush()
                .map_err(|error| io_error(&temporary, "flush", error))?;
            drop(file);
            fs::rename(&temporary, &manifest_path)
                .map_err(|error| io_error(&manifest_path, "finalize_manifest", error))?;
            let reloaded = fs::read(&manifest_path)
                .map_err(|error| io_error(&manifest_path, "read_manifest", error))?;
            validate_general_manifest_bytes(&manifest_path, &reloaded)?;
            let reparsed = parse_general_manifest_bytes(&manifest_path, &reloaded)?;
            compare_catalogs(&reparsed, &self.artifacts)?;
            self.validate_effective_config(&manifest_path, &reloaded)
        })();
        if result.is_err() {
            let _ = fs::remove_file(&temporary);
        }
        result?;
        self.manifest_written = true;
        Ok(())
    }

    /// Publica con `renameat2(RENAME_NOREPLACE)`. El marcador permanece en el
    /// wrapper exterior; el árbol renombrado nunca contiene `INCOMPLETE`.
    pub fn publish(self) -> Result<PathBuf, Diagnostic> {
        if !self.manifest_written {
            return Err(invalid_state(
                "manifest-general",
                "a validated general manifest is required before publication",
            ));
        }
        self.verify_catalog()?;
        let marker = self.provisional.join(INCOMPLETE_MARKER);
        if fs::read(&marker).map_err(|error| io_error(&marker, "read_marker", error))?
            != b"incomplete\n"
        {
            return Err(invalid_state(
                "execution",
                "the provisional INCOMPLETE marker is invalid",
            ));
        }
        let manifest_path = self.payload.join("manifest.json");
        let manifest = fs::read(&manifest_path)
            .map_err(|error| io_error(&manifest_path, "read_manifest", error))?;
        let declared = parse_general_manifest_bytes(&manifest_path, &manifest)?;
        compare_catalogs(&declared, &self.artifacts)?;
        self.validate_effective_config(&manifest_path, &manifest)?;
        self.rename_payload()?;
        // Desde este punto el destino ya es un árbol final limpio. El wrapper
        // exterior solo registra un residuo de limpieza si el proceso cae.
        // El rename es el punto de commit. Una interrupción o un error de
        // limpieza posterior solo puede dejar el wrapper privado exterior; no
        // convierte un destino final ya íntegro en una ejecución fallida.
        let _ = fs::remove_file(&marker);
        let _ = fs::remove_dir(&self.provisional);
        Ok(self.destination)
    }

    fn rename_payload(&self) -> Result<(), Diagnostic> {
        rename_no_replace(&self.payload, &self.destination)
    }

    fn check_new_artifact(&self, id: &str, path: &Path) -> Result<(), Diagnostic> {
        validate_relative_path(path)?;
        if path == Path::new(INCOMPLETE_MARKER) {
            return Err(invalid_state(
                id,
                "INCOMPLETE is reserved for the provisional wrapper marker",
            ));
        }
        if id.is_empty() {
            return Err(invalid_state("artifact", "artifact ID is empty"));
        }
        let path_text = path.to_string_lossy().replace('\\', "/");
        if self
            .artifacts
            .iter()
            .any(|artifact| artifact.id == id || artifact.path == path_text)
        {
            return Err(invalid_state(
                id,
                "artifact ID or path is already cataloged",
            ));
        }
        if self.manifest_written {
            return Err(invalid_state(id, "the artifact catalog is already closed"));
        }
        Ok(())
    }

    pub(crate) fn discard_last_artifact(
        &mut self,
        expected: &ArtifactRecord,
    ) -> Result<(), Diagnostic> {
        if self.manifest_written || self.artifacts.last() != Some(expected) {
            return Err(invalid_state(
                &expected.id,
                "only the last artifact may be discarded before the manifest is finalized",
            ));
        }
        let path = self.payload.join(&expected.path);
        fs::remove_file(&path).map_err(|error| io_error(&path, "discard_artifact", error))?;
        self.artifacts.pop();
        Ok(())
    }

    pub(crate) fn retain_only_artifact(
        &mut self,
        retained_id: Option<&str>,
    ) -> Result<(), Diagnostic> {
        if self.manifest_written {
            return Err(invalid_state(
                "artifact-catalog",
                "cannot reduce a finalized artifact catalog",
            ));
        }
        for artifact in self
            .artifacts
            .iter()
            .filter(|artifact| Some(artifact.id.as_str()) != retained_id)
        {
            let path = self.payload.join(&artifact.path);
            fs::remove_file(&path)
                .map_err(|error| io_error(&path, "discard_failed_artifact", error))?;
        }
        self.artifacts
            .retain(|artifact| Some(artifact.id.as_str()) == retained_id);
        Ok(())
    }

    #[cfg(test)]
    pub(crate) fn replace_artifact_bytes_for_validation_test(
        &mut self,
        id: &str,
        bytes: &[u8],
    ) -> Result<(), Diagnostic> {
        let index = self
            .artifacts
            .iter()
            .position(|artifact| artifact.id == id)
            .ok_or_else(|| invalid_state(id, "test artifact is absent"))?;
        let previous = self.artifacts[index].clone();
        let path = self.payload.join(&previous.path);
        fs::write(&path, bytes).map_err(|error| io_error(&path, "replace_test_artifact", error))?;
        let identity = identify_file(&path, Path::new(&previous.path))?;
        self.artifacts[index] =
            ArtifactRecord::from_identity(previous.id, previous.kind, &identity)?;
        Ok(())
    }

    fn verify_catalog(&self) -> Result<(), Diagnostic> {
        let mut ids = BTreeSet::new();
        let mut paths = BTreeSet::new();
        for artifact in &self.artifacts {
            if !ids.insert(&artifact.id) || !paths.insert(&artifact.path) {
                return Err(invalid_state(
                    &artifact.id,
                    "artifact ID or path is duplicated",
                ));
            }
            let path = self.payload.join(&artifact.path);
            let identity = identify_file(&path, Path::new(&artifact.path))?;
            if identity.size_bytes() != artifact.size_bytes || identity.digest() != artifact.digest
            {
                return Err(Diagnostic::io(
                    "GPS-IO-002",
                    "El tamaño o SHA-256 de un artefacto no coincide con el valor esperado.",
                )
                .with("path", path.display().to_string())
                .with("artifact_id", &artifact.id));
            }
        }
        Ok(())
    }

    fn validate_effective_config(
        &self,
        manifest_path: &Path,
        manifest_bytes: &[u8],
    ) -> Result<(), Diagnostic> {
        let document: serde_json::Value =
            serde_json::from_slice(manifest_bytes).map_err(|error| {
                invalid_state(
                    "manifest-general",
                    format!("manifest JSON changed: {error}"),
                )
            })?;
        let effective_id = document
            .get("effective_config_artifact_id")
            .and_then(serde_json::Value::as_str);
        let Some(effective_id) = effective_id else {
            return Ok(());
        };
        let artifact = self
            .artifacts
            .iter()
            .find(|artifact| artifact.id == effective_id)
            .ok_or_else(|| {
                invalid_state(effective_id, "effective config is absent from catalog")
            })?;
        let effective_path = self.payload.join(&artifact.path);
        let effective_bytes = fs::read(&effective_path)
            .map_err(|error| io_error(&effective_path, "read_effective_config", error))?;
        validate_published_resources(
            &self.payload,
            &effective_path,
            &effective_bytes,
            manifest_path,
            manifest_bytes,
        )
    }
}

fn temporary_path(parent: &Path, final_path: &Path, id: &str) -> Result<PathBuf, Diagnostic> {
    let name = final_path
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| invalid_state(id, "artifact file name is not UTF-8"))?;
    Ok(parent.join(format!(
        ".{name}.gps-tmp-{}-{}",
        std::process::id(),
        TEMPORARY_COUNTER.fetch_add(1, Ordering::Relaxed)
    )))
}

fn compare_catalogs(
    declared: &[ArtifactRecord],
    observed: &[ArtifactRecord],
) -> Result<(), Diagnostic> {
    let declared = declared
        .iter()
        .map(|record| (record.id.as_str(), record))
        .collect::<BTreeMap<_, _>>();
    let observed = observed
        .iter()
        .map(|record| (record.id.as_str(), record))
        .collect::<BTreeMap<_, _>>();
    if declared == observed {
        Ok(())
    } else {
        Err(invalid_state(
            "manifest-general",
            "the manifest artifact catalog differs from the workspace catalog",
        ))
    }
}

fn rename_no_replace(provisional: &Path, destination: &Path) -> Result<(), Diagnostic> {
    renameat_with(CWD, provisional, CWD, destination, RenameFlags::NOREPLACE)
        .map_err(|error| rename_error(destination, error))
}

fn rename_error(destination: &Path, error: rustix::io::Errno) -> Diagnostic {
    if error == rustix::io::Errno::EXIST {
        destination_exists(destination)
    } else if matches!(
        error,
        rustix::io::Errno::NOSYS | rustix::io::Errno::INVAL | rustix::io::Errno::OPNOTSUPP
    ) {
        Diagnostic::io(
            "GPS-IO-004",
            "El sistema de ficheros no ofrece publicación atómica no-clobber.",
        )
        .with("path", destination.display().to_string())
        .with("detail", error.to_string())
    } else {
        io_error(destination, "renameat2(RENAME_NOREPLACE)", error.into())
    }
}

fn validate_relative_path(path: &Path) -> Result<(), Diagnostic> {
    let text = path
        .to_str()
        .ok_or_else(|| invalid_state("artifact", "artifact path is not UTF-8"))?;
    validate_relative_link(text)?;
    if path
        .components()
        .any(|component| !matches!(component, Component::Normal(_)))
    {
        return Err(invalid_state("artifact", "artifact path is not canonical"));
    }
    Ok(())
}

fn destination_exists(destination: &Path) -> Diagnostic {
    Diagnostic::io(
        "GPS-IO-003",
        "El destino final ya existe y no puede sobrescribirse.",
    )
    .with("path", destination.display().to_string())
}

fn invalid_state(id: &str, detail: impl Into<String>) -> Diagnostic {
    Diagnostic::internal(
        "GPS-INTERNAL-001",
        "Una invariante interna considerada imposible se ha roto.",
    )
    .with("component", "publication")
    .with("id", id)
    .with("detail", detail.into())
}

fn input_error(path: &Path, error: std::io::Error) -> Diagnostic {
    Diagnostic::input("GPS-INPUT-001", "No se pudo leer una entrada declarada.")
        .with("path", path.display().to_string())
        .with("detail", error.to_string())
}

fn io_error(path: &Path, operation: &str, error: std::io::Error) -> Diagnostic {
    Diagnostic::io("GPS-IO-001", "Una operación de lectura o escritura falló.")
        .with("path", path.display().to_string())
        .with("operation", operation)
        .with("detail", error.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::results::{
        ConfigurationInputRecord, GeneralManifest, ManifestReferenceLink, UnitState,
        general_manifest_document, json_document_bytes,
    };

    fn temporary(label: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "gps-publication-unit-{}-{}-{label}",
            std::process::id(),
            TEMPORARY_COUNTER.fetch_add(1, Ordering::Relaxed)
        ))
    }

    fn ready_failed_workspace(destination: &Path) -> ExecutionWorkspace {
        let mut workspace = ExecutionWorkspace::create(destination).expect("workspace");
        let diagnostic = Diagnostic::format("GPS-FMT-001", "fallo controlado");
        let manifest = GeneralManifest {
            execution_id: "execution-publication-unit",
            state: UnitState::Failed,
            configuration_input: None,
            effective_config: None,
            references: &[],
            reduced_references: &[],
            artifacts: &[],
            diagnostics: std::slice::from_ref(&diagnostic),
        };
        workspace
            .write_general_manifest(
                &json_document_bytes(&general_manifest_document(&manifest)).expect("JSON"),
            )
            .expect("manifest");
        workspace
    }

    #[test]
    fn crash_after_payload_rename_never_exposes_incomplete_in_final_tree() {
        let root = temporary("rename-boundary");
        fs::create_dir(&root).expect("root");
        let destination = root.join("result");
        let workspace = ready_failed_workspace(&destination);
        let wrapper = workspace.provisional.clone();

        // Simula exactamente el corte entre el rename publicador y la limpieza.
        workspace.rename_payload().expect("rename payload");

        assert!(destination.join("manifest.json").is_file());
        assert!(!destination.join(INCOMPLETE_MARKER).exists());
        assert!(wrapper.join(INCOMPLETE_MARKER).is_file());
        fs::remove_dir_all(root).expect("cleanup");
    }

    #[test]
    fn unsupported_no_clobber_rename_maps_to_gps_io_004_and_exit_five() {
        let diagnostic = rename_error(Path::new("result"), rustix::io::Errno::OPNOTSUPP);
        assert_eq!(diagnostic.code(), "GPS-IO-004");
        assert_eq!(diagnostic.exit_code(), 5);
    }

    #[test]
    fn incomplete_marker_name_is_reserved_outside_the_payload() {
        let root = temporary("reserved-incomplete");
        fs::create_dir(&root).expect("root");
        let destination = root.join("result");
        let mut workspace = ExecutionWorkspace::create(&destination).expect("workspace");
        let error = workspace
            .write_artifact(
                "forbidden-marker",
                "test_artifact",
                Path::new(INCOMPLETE_MARKER),
                b"not a marker\n",
            )
            .expect_err("root INCOMPLETE is reserved");
        assert_eq!(error.code(), "GPS-FMT-006");
        assert!(!destination.exists());
        fs::remove_dir_all(root).expect("cleanup");
    }

    #[test]
    fn publication_rereads_and_rejects_invalid_effective_config() {
        let root = temporary("invalid-effective");
        fs::create_dir(&root).expect("root");
        let destination = root.join("result");
        let mut workspace = ExecutionWorkspace::create(&destination).expect("workspace");
        let config = workspace
            .write_artifact(
                "config-general",
                "config_general_toml",
                Path::new("inputs/config/general.toml"),
                b"config_kind = \"general\"\n",
            )
            .expect("config");
        let fasta = workspace
            .write_artifact(
                "reference-fasta",
                "reference_fasta",
                Path::new("resources/references/main.fa"),
                b">chr1\nA\n",
            )
            .expect("FASTA");
        let effective = workspace
            .write_artifact(
                "effective-config",
                "effective_config",
                Path::new("effective-config.json"),
                br#"{"effective_config_kind":"unknown"}"#,
            )
            .expect("invalid effective bytes are still a physical artifact");
        let configuration = ConfigurationInputRecord {
            kind: "general".to_owned(),
            original_path: "config/general.toml".to_owned(),
            copy: config,
        };
        let references = [ManifestReferenceLink {
            id: "main".to_owned(),
            fasta_artifact_id: fasta.id.clone(),
        }];
        let artifacts = workspace.artifacts().to_vec();
        let manifest = GeneralManifest {
            execution_id: "execution-invalid-effective",
            state: UnitState::Complete,
            configuration_input: Some(&configuration),
            effective_config: Some(&effective),
            references: &references,
            reduced_references: &[],
            artifacts: &artifacts,
            diagnostics: &[],
        };
        let error = workspace
            .write_general_manifest(
                &json_document_bytes(&general_manifest_document(&manifest)).expect("JSON"),
            )
            .expect_err("publication must reject the effective document");
        assert_eq!(error.code(), "GPS-FMT-006");
        assert_eq!(error.category(), "FMT");
        assert_eq!(error.exit_code(), 1);
        assert!(!destination.exists());
        fs::remove_dir_all(root).expect("cleanup");
    }
}
