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
    validate_published_event_manifests, validate_published_resources, validate_relative_link,
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
    /// El manifiesto que cerró el catálogo es el general v5 de F5.
    ///
    /// v5 no es una ampliación de v4: su raíz no tiene `artifacts`, ni
    /// `references`, ni `effective_config_artifact_id`, de modo que el lector
    /// v4 lo rechazaría por versión antes de mirar nada. Esta marca elige, en
    /// el cierre y otra vez en la publicación, la puerta física que
    /// corresponde a la familia escrita.
    f5_manifest: bool,
    /// Documentos raíz del paquete que ningún catálogo enumera.
    ///
    /// `manifest.json` nunca estuvo en el catálogo —su auto-hash no tendría
    /// solución finita— y `effective-config.json` v5 tampoco puede estarlo: el
    /// objeto `artifacts` del manifiesto v5 es una lista cerrada de
    /// exactamente los siete TSV. Se registran aquí para impedir que dos
    /// escrituras reclamen el mismo nombre o que un documento raíz pise la
    /// ruta de un artefacto.
    root_documents: BTreeSet<String>,
    /// Identificadores y rutas de artefactos abiertos en flujo y todavía sin
    /// publicar. El catálogo no los contiene aún, de modo que sin este
    /// registro dos salidas simultáneas podrían reclamar el mismo nombre.
    reserved: BTreeSet<(String, String)>,
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
                        f5_manifest: false,
                        root_documents: BTreeSet::new(),
                        reserved: BTreeSet::new(),
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

    /// Raíz de la carga útil que se renombrará al destino.
    ///
    /// La expone para que un lector independiente pueda juzgar los artefactos
    /// ya escritos **antes** de que exista el destino final: es lo que el
    /// recuento MAN-70 de F5.4.3 necesita, y hacerlo sobre el provisional es
    /// justamente lo que impide publicar un paquete que su propio oráculo
    /// rechazaría.
    #[must_use]
    pub fn payload_root(&self) -> &Path {
        &self.payload
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
        self.write_streamed_artifact_with_validation(id, kind, relative_path, write, |_| Ok(()))
    }

    /// Escribe un artefacto por flujo y valida el temporal ya cerrado antes de
    /// calcular su identidad y hacerlo visible dentro del provisional.
    ///
    /// El validador recibe el temporal, nunca la ruta final. Por tanto, un
    /// fallo de formato o de contadores no deja un artefacto parcial en el
    /// catálogo ni en el árbol de carga útil.
    pub fn write_streamed_artifact_with_validation(
        &mut self,
        id: &str,
        kind: &str,
        relative_path: &Path,
        write: impl FnOnce(&mut File) -> Result<(), Diagnostic>,
        validate: impl FnOnce(&Path) -> Result<(), Diagnostic>,
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
            validate(&temporary)?;
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

    /// Abre un artefacto cuyo temporal permanece disponible entre llamadas.
    ///
    /// A diferencia de `write_streamed_artifact_with_validation`, que abre,
    /// valida y publica en una sola llamada, esto permite tener varias salidas
    /// escribiéndose a la vez. Es lo que necesitan los FASTQ finales: R1 y R2
    /// de todas las bibliotecas se llenan en un único recorrido de las parejas,
    /// de modo que ninguna pareja tiene que conservarse en memoria.
    ///
    /// El artefacto no entra en el catálogo hasta `finish_streamed_artifact`.
    /// Abandonar el handle elimina su temporal y no publica nada.
    pub(crate) fn begin_streamed_artifact(
        &mut self,
        id: &str,
        kind: &str,
        relative_path: &Path,
    ) -> Result<OpenArtifact, Diagnostic> {
        self.check_new_artifact(id, relative_path)?;
        let final_path = self.payload.join(relative_path);
        let parent = final_path
            .parent()
            .ok_or_else(|| invalid_state(id, "artifact path has no provisional parent"))?;
        fs::create_dir_all(parent).map_err(|error| io_error(parent, "create_dir_all", error))?;
        let temporary = temporary_path(parent, &final_path, id)?;
        let file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary)
            .map_err(|error| io_error(&temporary, "create_artifact_temporary", error))?;
        let path_text = relative_path.to_string_lossy().replace('\\', "/");
        self.reserved.insert((id.to_owned(), path_text));
        Ok(OpenArtifact {
            id: id.to_owned(),
            kind: kind.to_owned(),
            relative_path: relative_path.to_path_buf(),
            final_path,
            temporary,
            file: Some(file),
            published: false,
        })
    }

    /// Cierra un artefacto abierto en flujo con las mismas comprobaciones que
    /// la vía de una sola llamada: validar el temporal cerrado, identificarlo,
    /// verificar la identidad y solo entonces hacerlo visible y catalogarlo.
    pub(crate) fn finish_streamed_artifact(
        &mut self,
        mut artifact: OpenArtifact,
        validate: impl FnOnce(&Path) -> Result<(), Diagnostic>,
    ) -> Result<ArtifactRecord, Diagnostic> {
        if self.manifest_written {
            return Err(invalid_state(
                &artifact.id,
                "the artifact catalog is already closed",
            ));
        }
        let result = (|| {
            let mut file = artifact
                .file
                .take()
                .ok_or_else(|| invalid_state(&artifact.id, "artifact stream is already closed"))?;
            file.flush()
                .map_err(|error| io_error(&artifact.temporary, "flush", error))?;
            drop(file);
            validate(&artifact.temporary)?;
            let identity = identify_file(&artifact.temporary, &artifact.relative_path)?;
            verify_identity(&artifact.temporary, &identity)?;
            fs::rename(&artifact.temporary, &artifact.final_path)
                .map_err(|error| io_error(&artifact.final_path, "finalize_artifact", error))?;
            ArtifactRecord::from_identity(&artifact.id, &artifact.kind, &identity)
        })();
        let path_text = artifact.relative_path.to_string_lossy().replace('\\', "/");
        self.reserved.remove(&(artifact.id.clone(), path_text));
        let record = result?;
        artifact.published = true;
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

    /// Escribe un documento raíz del paquete que ningún catálogo enumera.
    ///
    /// Es la vía de `effective-config.json`: el manifiesto v5 no puede
    /// enlazarlo —su lista `artifacts` está cerrada a los siete TSV— y por
    /// tanto tampoco puede compararse con el catálogo. Sigue escribiéndose
    /// como todos los demás: a un temporal propio que se renombra ya completo,
    /// de modo que una interrupción no deja un documento raíz a medias.
    ///
    /// # Errors
    ///
    /// `GPS-INTERNAL-001` si el nombre no es un componente relativo simple, si
    /// ya lo reclama un artefacto u otro documento raíz o si el catálogo ya
    /// está cerrado; `GPS-IO-001` ante un fallo de escritura.
    pub fn write_root_document(&mut self, name: &str, bytes: &[u8]) -> Result<(), Diagnostic> {
        if self.manifest_written {
            return Err(invalid_state(
                name,
                "the artifact catalog is already closed",
            ));
        }
        let relative = Path::new(name);
        validate_relative_path(relative)?;
        if relative.components().count() != 1 {
            return Err(invalid_state(
                name,
                "a root document is a single path component",
            ));
        }
        if name == "manifest.json" {
            return Err(invalid_state(
                name,
                "manifest.json is written by its own closing gate",
            ));
        }
        if self.artifacts.iter().any(|artifact| artifact.path == name)
            || self.reserved.iter().any(|(_, path)| path == name)
            || !self.root_documents.insert(name.to_owned())
        {
            return Err(invalid_state(
                name,
                "root document path is already cataloged, reserved or written",
            ));
        }
        let final_path = self.payload.join(relative);
        let temporary = temporary_path(&self.payload, &final_path, name)?;
        let result = (|| {
            let mut file = OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&temporary)
                .map_err(|error| io_error(&temporary, "create_root_document_temporary", error))?;
            file.write_all(bytes)
                .map_err(|error| io_error(&temporary, "write", error))?;
            file.flush()
                .map_err(|error| io_error(&temporary, "flush", error))?;
            drop(file);
            fs::rename(&temporary, &final_path)
                .map_err(|error| io_error(&final_path, "finalize_root_document", error))
        })();
        if result.is_err() {
            let _ = fs::remove_file(&temporary);
            self.root_documents.remove(name);
        }
        result
    }

    /// Cierra el catálogo con el manifiesto general v5 de F5.
    ///
    /// Escribe el temporal, lo renombra, **relee los bytes publicados** y los
    /// juzga con la puerta física de F5, que vuelve a abrir cada artefacto
    /// enlazado y recomputa su tamaño y SHA-256. No comparte ninguna
    /// comprobación con la familia v4: ni el catálogo `artifacts` de la raíz,
    /// ni los manifiestos de evento, ni el efectivo v4 existen aquí.
    ///
    /// # Errors
    ///
    /// `GPS-INTERNAL-001` si el manifiesto ya está cerrado, `GPS-FMT-006` si
    /// el documento publicado no satisface el contrato v5 o no corresponde al
    /// paquete físico, y `GPS-IO-001`/`GPS-IO-002` ante fallos de escritura o
    /// de identidad.
    pub fn write_f5_manifest<T: serde::Serialize>(
        &mut self,
        document: &T,
    ) -> Result<(), Diagnostic> {
        if self.manifest_written {
            return Err(invalid_state(
                "manifest-general",
                "the general manifest is already finalized",
            ));
        }
        self.verify_catalog()?;
        let manifest_path = self.payload.join("manifest.json");
        let temporary = self.payload.join(format!(
            ".manifest.json.gps-tmp-{}-{}",
            std::process::id(),
            TEMPORARY_COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        let result = (|| {
            let file = OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&temporary)
                .map_err(|error| io_error(&temporary, "create_manifest_temporary", error))?;
            // El documento se serializa **directamente al fichero**: el
            // manifiesto v5 crece con `B*N_pop*S` y materializar su texto
            // costaba, solo él, más que todo el pico que la admisión sella.
            let mut writer = std::io::BufWriter::new(file);
            serde_json::to_writer_pretty(&mut writer, document).map_err(|error| {
                invalid_state(
                    "manifest-general",
                    format!("manifest is not serializable: {error}"),
                )
            })?;
            writer
                .write_all(b"\n")
                .map_err(|error| io_error(&temporary, "write", error))?;
            let mut file = writer
                .into_inner()
                .map_err(|error| io_error(&temporary, "flush", error.into_error()))?;
            file.flush()
                .map_err(|error| io_error(&temporary, "flush", error))?;
            drop(file);
            fs::rename(&temporary, &manifest_path)
                .map_err(|error| io_error(&manifest_path, "finalize_manifest", error))
        })();
        if result.is_err() {
            let _ = fs::remove_file(&temporary);
        }
        result
    }

    /// Ejecuta la puerta física sobre el manifiesto v5 ya escrito y cierra el
    /// catálogo.
    ///
    /// Está separada de [`Self::write_f5_manifest`] para que el llamador pueda
    /// soltar el documento que lo produjo **antes** de juzgarlo: la puerta
    /// relee el fichero publicado y no debe coexistir con la copia en memoria.
    /// Sin esa separación, el objeto `f5` y los resúmenes MAN-70 seguían vivos
    /// durante la verificación y sumaban su tamaño al pico.
    ///
    /// # Errors
    ///
    /// `GPS-INTERNAL-001` si el manifiesto ya está cerrado o si no se ha
    /// escrito, y `GPS-FMT-006`/`GPS-IO-001` si la puerta lo rechaza.
    pub fn close_f5_manifest(&mut self) -> Result<(), Diagnostic> {
        if self.manifest_written {
            return Err(invalid_state(
                "manifest-general",
                "the general manifest is already finalized",
            ));
        }
        let manifest_path = self.payload.join("manifest.json");
        if !manifest_path.is_file() {
            return Err(invalid_state(
                "manifest-general",
                "the F5 manifest has not been written yet",
            ));
        }
        crate::f5_publication::verify_published_f5_package(
            &self.payload,
            &manifest_path,
            &self.artifacts,
        )?;
        self.manifest_written = true;
        self.f5_manifest = true;
        Ok(())
    }

    /// Escribe, relee y valida físicamente el manifiesto general. `manifest.json`
    /// no entra en el catálogo porque su auto-hash no tendría solución finita.
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
            validate_published_event_manifests(&self.payload, &manifest_path, &reloaded)?;
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
        // Revalida todos los enlaces físicos sobre los bytes que realmente se
        // van a publicar. Entre el cierre del manifiesto y este punto todavía
        // puede haberse alterado el árbol provisional. La familia v5 lee el
        // manifiesto por flujo y no lo materializa; la v4 conserva su lectura
        // completa, porque sus lectores la exigen.
        if self.f5_manifest {
            crate::f5_publication::verify_published_f5_package(
                &self.payload,
                &manifest_path,
                &self.artifacts,
            )?;
        } else {
            let manifest = fs::read(&manifest_path)
                .map_err(|error| io_error(&manifest_path, "read_manifest", error))?;
            let declared = parse_general_manifest_bytes(&manifest_path, &manifest)?;
            compare_catalogs(&declared, &self.artifacts)?;
            validate_published_event_manifests(&self.payload, &manifest_path, &manifest)?;
            self.validate_effective_config(&manifest_path, &manifest)?;
        }
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
        if self
            .reserved
            .iter()
            .any(|(held_id, held_path)| held_id == id || *held_path == path_text)
        {
            return Err(invalid_state(
                id,
                "artifact ID or path is already reserved by an open stream",
            ));
        }
        // El invariante de `root_documents` es simétrico: un documento raíz no
        // puede pisar la ruta de un artefacto **ni al revés**. Sin esta mitad,
        // un artefacto podía reclamar `effective-config.json` ya escrito,
        // sobrescribirlo y entrar además en el catálogo.
        if self.root_documents.contains(&path_text) {
            return Err(invalid_state(
                id,
                "artifact path is already written as a root document",
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

    /// Conserva únicamente artefactos independientes ya validados al cerrar
    /// un fallo controlado. No reabre ni completa temporales: los demás se
    /// eliminan antes de formar el manifiesto final fallido.
    pub(crate) fn retain_only_artifacts(
        &mut self,
        retained_ids: &BTreeSet<String>,
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
            .filter(|artifact| !retained_ids.contains(&artifact.id))
        {
            let path = self.payload.join(&artifact.path);
            fs::remove_file(&path)
                .map_err(|error| io_error(&path, "discard_failed_artifact", error))?;
        }
        self.artifacts
            .retain(|artifact| retained_ids.contains(&artifact.id));
        Ok(())
    }

    /// Retira **todos** los artefactos ya escritos y sus ficheros.
    ///
    /// Es la puerta de un fallo controlado F5: DTC-089 prohíbe que el paquete
    /// fallido liste ninguno de los siete TSV «aunque existieran temporales»,
    /// de modo que el catálogo se vacía y los ficheros se borran antes de
    /// escribir el diagnóstico. Reutiliza `retain_only_artifacts` con el
    /// conjunto vacío para no tener dos borradores del mismo invariante.
    ///
    /// # Errors
    ///
    /// `GPS-INTERNAL-001` si el catálogo ya está cerrado y `GPS-IO-001` si un
    /// fichero no puede retirarse.
    pub fn discard_all_artifacts(&mut self) -> Result<(), Diagnostic> {
        self.retain_only_artifacts(&BTreeSet::new())
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

/// Artefacto en curso cuyo temporal vive entre llamadas.
///
/// Solo `ExecutionWorkspace::finish_streamed_artifact` lo publica. Si el
/// handle se destruye antes —por un error o un desenrollado— el temporal se
/// elimina, de modo que un fallo a mitad de escritura no deja restos en el
/// provisional ni un artefacto parcial en el catálogo.
#[derive(Debug)]
pub(crate) struct OpenArtifact {
    id: String,
    kind: String,
    relative_path: PathBuf,
    final_path: PathBuf,
    temporary: PathBuf,
    file: Option<File>,
    published: bool,
}

impl OpenArtifact {
    /// Cede el fichero para envolverlo en un codificador. El handle queda sin
    /// fichero hasta que `restore` lo devuelva, de manera que un cierre
    /// olvidado se detecta en `finish_streamed_artifact`.
    pub(crate) fn take_file(&mut self) -> Result<File, Diagnostic> {
        self.file
            .take()
            .ok_or_else(|| invalid_state(&self.id, "artifact stream has no open file"))
    }

    pub(crate) fn restore_file(&mut self, file: File) {
        self.file = Some(file);
    }
}

impl Drop for OpenArtifact {
    fn drop(&mut self) {
        if !self.published {
            self.file = None;
            let _ = fs::remove_file(&self.temporary);
        }
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
    use crate::application::build_effective_config;
    use crate::configuration::ExpectedConfigKind;
    use crate::event_manifest::{
        EventManifest, EventManifestArtifact, EventManifestDemultiplexing, EventManifestReference,
        event_manifest_bytes, event_manifest_relative_path, validate_event_manifest_bytes,
    };
    use crate::identity::ReferenceId;
    use crate::randomness::RootSeed;
    use crate::results::{
        ConfigurationInputRecord, GeneralManifest, ManifestEventLink, ManifestReferenceLink,
        UnitState, general_manifest_document, json_document_bytes,
    };
    use crate::technical_identity::derive_technical_identity_plan;

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
            population_snvs: None,
            vcf_outputs: &[],
            events: None,
            authorizations: crate::execution::ExecutionAuthorizations::default(),
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
            population_snvs: None,
            vcf_outputs: &[],
            events: None,
            authorizations: crate::execution::ExecutionAuthorizations::default(),
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

    #[test]
    fn publication_rejects_an_artifact_mutated_after_manifest_finalization() {
        let root = temporary("artifact-mutated-before-publish");
        fs::create_dir(&root).expect("root");
        let destination = root.join("result");
        let mut workspace = ExecutionWorkspace::create(&destination).expect("workspace");
        let artifact = workspace
            .write_artifact(
                "config-general",
                "config_general_toml",
                Path::new("inputs/config/general.toml"),
                b"config_kind = \"general\"\n",
            )
            .expect("artifact");
        let configuration = ConfigurationInputRecord {
            kind: "general".to_owned(),
            original_path: "config/general.toml".to_owned(),
            copy: artifact,
        };
        let diagnostic = Diagnostic::format("GPS-FMT-001", "fallo controlado");
        let artifacts = workspace.artifacts().to_vec();
        let manifest = GeneralManifest {
            execution_id: "execution-mutated-artifact",
            state: UnitState::Failed,
            configuration_input: Some(&configuration),
            effective_config: None,
            references: &[],
            reduced_references: &[],
            population_snvs: None,
            vcf_outputs: &[],
            events: None,
            authorizations: crate::execution::ExecutionAuthorizations::default(),
            artifacts: &artifacts,
            diagnostics: std::slice::from_ref(&diagnostic),
        };
        workspace
            .write_general_manifest(
                &json_document_bytes(&general_manifest_document(&manifest)).expect("JSON"),
            )
            .expect("manifest");

        fs::write(
            workspace.root().join("inputs/config/general.toml"),
            b"config_kind = \"mutated\"\n",
        )
        .expect("mutate provisional artifact");

        let error = workspace
            .publish()
            .expect_err("a modified artifact must fail the lightweight catalog check");
        assert_eq!(error.code(), "GPS-IO-002");
        assert!(!destination.exists());
        fs::remove_dir_all(root).expect("cleanup");
    }

    #[test]
    #[ignore = "F4.1 retira la publicación FASTQ SE; F4.2/F4.5 cubrirán la pareja publicada"]
    fn publication_revalidates_event_manifests_immediately_before_rename() {
        let root = temporary("event-tampered-before-publish");
        fs::create_dir(&root).expect("root");
        let destination = root.join("result");
        let mut workspace = ExecutionWorkspace::create(&destination).expect("workspace");
        let reference = workspace
            .write_artifact(
                "reference-main-fasta",
                "reference_fasta",
                Path::new("resources/references/reference-main.fa"),
                b">chr1\nACGT\n",
            )
            .expect("reference");
        let effective = effective_technical_event();
        let event = effective.technical_event().expect("technical event");
        let reference_id = ReferenceId::parse("reference-main", "test").expect("reference ID");
        let identities = derive_technical_identity_plan(
            RootSeed::parse(
                "root",
                "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f",
            )
            .expect("seed"),
            &reference_id,
            event.event_id().clone(),
            event.individual_ids(),
        )
        .expect("identities");
        let libraries = identities.libraries().collect::<Vec<_>>();
        let fastqs = libraries
            .iter()
            .map(|library| {
                workspace.write_artifact(
                    &format!("fastq-{}", library.library_id().as_str()),
                    "fastq",
                    &PathBuf::from("outputs/fastq")
                        .join(event.event_id().as_str())
                        .join(format!("{}.fastq.gz", library.individual_id().as_str())),
                    b"@read\nA\n+\nI\n",
                )
            })
            .collect::<Result<Vec<_>, _>>()
            .expect("complete FASTQ artifacts");
        let complete_artifacts = libraries
            .iter()
            .zip(&fastqs)
            .map(|(library, artifact)| EventManifestArtifact {
                id: &artifact.id,
                kind: "fastq",
                state: UnitState::Complete,
                individual_id: library.individual_id(),
                library_id: library.library_id(),
                mate: crate::technical_identity::TechnicalMate::R1,
                artifact: Some(artifact),
                diagnostic: None,
            })
            .collect::<Vec<_>>();
        let complete_event = EventManifest {
            event,
            state: UnitState::Complete,
            identities: &identities,
            reference: EventManifestReference {
                reference_id: &reference_id,
                fasta: &reference,
            },
            artifacts: &complete_artifacts,
            coverage: &[],
            final_fastq_summary: &[],
            error_quality: crate::event_manifest::EventManifestErrorQuality::neutral(),
            demultiplexing: EventManifestDemultiplexing::disabled(),
            chimera_rate: event.chimera_rate(),
            realized_chimera_pairs: 0,
            biological_contamination: &[],
            diagnostics: &[],
        };
        let event_path = event_manifest_relative_path(event.event_id());
        let initial_event_bytes = event_manifest_bytes(&complete_event).expect("event bytes");
        let event_record = workspace
            .write_streamed_artifact_with_validation(
                "event-manifest-reads-main",
                "event_manifest",
                &event_path,
                |file| {
                    file.write_all(&initial_event_bytes)
                        .map_err(|error| io_error(&event_path, "write_test_event_manifest", error))
                },
                |path| {
                    let bytes = fs::read(path)
                        .map_err(|error| io_error(path, "read_test_event_manifest", error))?;
                    validate_event_manifest_bytes(&event_path, &bytes)
                },
            )
            .expect("event manifest");
        let root_failure = Diagnostic::format("GPS-FMT-001", "fallo controlado");
        let events = [ManifestEventLink {
            event_id: event.event_id().as_str().to_owned(),
            state: UnitState::Complete,
            manifest_artifact_id: event_record.id.clone(),
            chimera_rate: event.chimera_rate(),
            realized_chimera_pairs: 0,
        }];
        let references = [ManifestReferenceLink {
            id: reference_id.as_str().to_owned(),
            fasta_artifact_id: reference.id.clone(),
        }];
        let catalog = workspace.artifacts().to_vec();
        let general = GeneralManifest {
            execution_id: "event-tamper-regression",
            state: UnitState::Failed,
            configuration_input: None,
            effective_config: None,
            references: &references,
            reduced_references: &[],
            population_snvs: None,
            vcf_outputs: &[],
            events: Some(&events),
            authorizations: crate::execution::ExecutionAuthorizations::default(),
            artifacts: &catalog,
            diagnostics: std::slice::from_ref(&root_failure),
        };
        workspace
            .write_general_manifest(
                &json_document_bytes(&general_manifest_document(&general)).expect("general bytes"),
            )
            .expect("finalized general manifest");

        // Simula una manipulación con acceso al provisional: reescribe ambos
        // catálogos de hash, pero deja el enlace raíz declarado como complete.
        let event_failure = Diagnostic::internal_invariant("test", "event-tampered");
        let failed_artifacts = [
            EventManifestArtifact {
                id: &fastqs[0].id,
                kind: "fastq",
                state: UnitState::Complete,
                individual_id: libraries[0].individual_id(),
                library_id: libraries[0].library_id(),
                mate: crate::technical_identity::TechnicalMate::R1,
                artifact: Some(&fastqs[0]),
                diagnostic: None,
            },
            EventManifestArtifact {
                id: "fastq-unclosed-library",
                kind: "fastq",
                state: UnitState::Failed,
                individual_id: libraries[1].individual_id(),
                library_id: libraries[1].library_id(),
                mate: crate::technical_identity::TechnicalMate::R1,
                artifact: None,
                diagnostic: Some(&event_failure),
            },
        ];
        let failed_event = EventManifest {
            event,
            state: UnitState::Failed,
            identities: &identities,
            reference: EventManifestReference {
                reference_id: &reference_id,
                fasta: &reference,
            },
            artifacts: &failed_artifacts,
            coverage: &[],
            final_fastq_summary: &[],
            error_quality: crate::event_manifest::EventManifestErrorQuality::neutral(),
            demultiplexing: EventManifestDemultiplexing::disabled(),
            chimera_rate: event.chimera_rate(),
            realized_chimera_pairs: 0,
            biological_contamination: &[],
            diagnostics: std::slice::from_ref(&event_failure),
        };
        let replacement = event_manifest_bytes(&failed_event).expect("tampered event bytes");
        let physical_event_path = workspace.root().join(&event_path);
        fs::write(&physical_event_path, replacement).expect("tamper event manifest");
        let replacement_identity =
            identify_file(&physical_event_path, &event_path).expect("identity");
        let replacement_record = ArtifactRecord::from_identity(
            event_record.id.clone(),
            event_record.kind.clone(),
            &replacement_identity,
        )
        .expect("replacement record");
        let index = workspace
            .artifacts
            .iter()
            .position(|artifact| artifact.id == event_record.id)
            .expect("event catalog record");
        workspace.artifacts[index] = replacement_record.clone();
        let manifest_path = workspace.root().join("manifest.json");
        let mut root_document: serde_json::Value =
            serde_json::from_slice(&fs::read(&manifest_path).expect("general bytes"))
                .expect("general JSON");
        *root_document["artifacts"]
            .as_array_mut()
            .expect("catalog")
            .iter_mut()
            .find(|artifact| artifact["id"] == event_record.id)
            .expect("event catalog JSON") = replacement_record.to_json();
        fs::write(
            &manifest_path,
            json_document_bytes(&root_document).expect("rewritten general bytes"),
        )
        .expect("tamper general manifest");

        let error = workspace
            .publish()
            .expect_err("event link altered after finalization must prevent rename");
        assert_eq!(error.code(), "GPS-FMT-006");
        assert!(!destination.exists());
        fs::remove_dir_all(root).expect("cleanup");
    }

    fn effective_technical_event() -> crate::application::EffectiveConfig {
        const CONFIG: &str = r#"config_kind = "general"
schema_version = "2"
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"
[population_snvs]
mode = "generate"
[population_snvs.generate]
count_mode = "exact"
locus_count = 1
spatial_mode = "uniform"
substitution_profile = "uniform_snv12_v1"
[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.2
max_alt_af = 0.8
[reference]
kind = "simulated"
id = "reference-main"
[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 200
[[initial_populations]]
population_id = "P1"
individual_count = 2
individual_ids = ["sample-A", "sample-B"]
ploidy = 2
genotype_model = "hwe"
[technical_event]
event_id = "reads-main"
individual_ids = ["sample-A", "sample-B"]
target_coverage_wgs = 1.0
[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
[[technical_event.libraries]]
individual_id = "sample-B"
i7 = "AAAACCCC"
i5 = "GGGGTTTT"
"#;
        let admitted = crate::adapters::admit_toml_header(
            Path::new("publication-event-regression.toml"),
            CONFIG.as_bytes().to_vec(),
            ExpectedConfigKind::General,
        )
        .expect("header");
        build_effective_config(&crate::adapters::parse_declared_config(admitted).expect("declared"))
            .expect("effective")
    }

    #[test]
    fn a_root_document_claims_its_name_and_never_enters_the_catalog() {
        let root = temporary("root-document");
        fs::create_dir(&root).expect("root");
        let mut workspace = ExecutionWorkspace::create(&root.join("published")).expect("workspace");
        workspace
            .write_artifact(
                "sites",
                "f5_sites",
                Path::new("outputs/f5-sites.tsv"),
                b"header\n",
            )
            .expect("artifact");

        workspace
            .write_root_document("effective-config.json", b"{}\n")
            .expect("root document");
        assert_eq!(
            fs::read(workspace.root().join("effective-config.json")).expect("reread"),
            b"{}\n"
        );
        // No entra en el catálogo: el manifiesto v5 no puede enlazarlo.
        assert_eq!(workspace.artifacts().len(), 1);
        assert_eq!(workspace.artifacts()[0].id, "sites");

        // Nombre repetido, ruta de un artefacto, `manifest.json` y ruta
        // compuesta se rechazan por igual.
        for (name, reason) in [
            ("effective-config.json", "repetido"),
            ("outputs/f5-sites.tsv", "ruta de artefacto"),
            ("manifest.json", "manifiesto"),
            ("outputs/extra.json", "ruta compuesta"),
            ("INCOMPLETE", "marcador"),
        ] {
            assert!(
                workspace.write_root_document(name, b"{}\n").is_err(),
                "{name} debe rechazarse: {reason}"
            );
        }
        // La dirección contraria del mismo invariante: un artefacto no puede
        // reclamar la ruta de un documento raíz ya escrito. Antes de la
        // revisión la sobrescribía devolviendo `Ok` y se catalogaba.
        let collision =
            workspace.write_artifact("colision", "f5", Path::new("effective-config.json"), b"x\n");
        assert!(
            collision.is_err(),
            "un artefacto no puede pisar un documento raíz"
        );
        assert_eq!(workspace.artifacts().len(), 1);
        assert_eq!(
            fs::read(workspace.root().join("effective-config.json")).expect("reread"),
            b"{}\n",
            "el documento raíz conserva sus bytes"
        );
        // Y ningún temporal sobrevive a los rechazos.
        let residue = fs::read_dir(workspace.root())
            .expect("payload")
            .filter_map(Result::ok)
            .any(|entry| entry.file_name().to_string_lossy().contains(".gps-tmp-"));
        assert!(!residue);
        let _ = fs::remove_dir_all(&root);
    }

    #[test]
    fn discarding_all_artifacts_empties_the_catalog_and_the_tree() {
        let root = temporary("discard-all");
        fs::create_dir(&root).expect("root");
        let mut workspace = ExecutionWorkspace::create(&root.join("published")).expect("workspace");
        for (id, path) in [
            ("sites", "outputs/f5-sites.tsv"),
            ("truth", "outputs/truth.jsonl"),
        ] {
            workspace
                .write_artifact(id, "f5", Path::new(path), b"x\n")
                .expect("artifact");
        }
        assert_eq!(workspace.artifacts().len(), 2);
        workspace.discard_all_artifacts().expect("discard");
        assert!(workspace.artifacts().is_empty());
        assert!(!workspace.root().join("outputs/f5-sites.tsv").exists());
        assert!(!workspace.root().join("outputs/truth.jsonl").exists());
        // El mismo identificador vuelve a estar libre.
        workspace
            .write_artifact("sites", "f5", Path::new("outputs/f5-sites.tsv"), b"y\n")
            .expect("reuse");
        let _ = fs::remove_dir_all(&root);
    }

    #[test]
    fn an_open_stream_reserves_its_name_and_leaves_no_temporary_when_abandoned() {
        let root = temporary("open-artifact");
        fs::create_dir(&root).expect("root");
        let destination = root.join("published");
        let mut workspace = ExecutionWorkspace::create(&destination).expect("workspace");
        let path = Path::new("outputs/stream.txt");
        let mut open = workspace
            .begin_streamed_artifact("stream", "text", path)
            .expect("open the artifact");

        // Mientras esté abierto, su nombre y su ruta están tomados aunque el
        // catálogo todavía no lo contenga.
        let repeated_id =
            workspace.begin_streamed_artifact("stream", "text", Path::new("other.txt"));
        assert!(repeated_id.is_err(), "a repeated ID must be rejected");
        let repeated_path = workspace.begin_streamed_artifact("other", "text", path);
        assert!(repeated_path.is_err(), "a repeated path must be rejected");
        assert!(workspace.artifacts().is_empty());

        // Un segundo artefacto distinto sí puede convivir con el primero.
        let mut second = workspace
            .begin_streamed_artifact("second", "text", Path::new("outputs/second.txt"))
            .expect("open a second artifact");
        let mut file = second.take_file().expect("second file");
        file.write_all(b"second\n").expect("write");
        second.restore_file(file);
        let record = workspace
            .finish_streamed_artifact(second, |_| Ok(()))
            .expect("publish the second artifact");
        assert_eq!(workspace.artifacts(), &[record]);

        // Abandonar el primero elimina su temporal y no publica nada.
        let mut file = open.take_file().expect("first file");
        file.write_all(b"abandoned\n").expect("write");
        open.restore_file(file);
        drop(open);
        assert!(!workspace.root().join(path).exists());
        let remaining = fs::read_dir(workspace.root().join("outputs"))
            .expect("outputs directory")
            .filter_map(Result::ok)
            .map(|entry| entry.file_name())
            .collect::<Vec<_>>();
        assert_eq!(
            remaining,
            vec![std::ffi::OsString::from("second.txt")],
            "only the published artifact may remain"
        );
        fs::remove_dir_all(root).expect("cleanup");
    }
}
