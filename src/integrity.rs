//! Integridad SHA-256 incremental sobre todos los bytes exactos.

use std::fs::File;
use std::io::{self, Read};
use std::path::{Path, PathBuf};

use sha2::{Digest, Sha256};

use crate::diagnostics::Diagnostic;

pub const INTEGRITY_ALGORITHM: &str = "sha256";

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ArtifactIdentity {
    relative_path: PathBuf,
    size_bytes: u64,
    digest: String,
}

impl ArtifactIdentity {
    #[must_use]
    pub fn relative_path(&self) -> &Path {
        &self.relative_path
    }

    #[must_use]
    pub const fn size_bytes(&self) -> u64 {
        self.size_bytes
    }

    #[must_use]
    pub const fn algorithm(&self) -> &'static str {
        INTEGRITY_ALGORITHM
    }

    #[must_use]
    pub fn digest(&self) -> &str {
        &self.digest
    }
}

#[derive(Debug)]
pub struct Sha256Writer {
    hasher: Sha256,
    size_bytes: u64,
}

impl Default for Sha256Writer {
    fn default() -> Self {
        Self {
            hasher: Sha256::new(),
            size_bytes: 0,
        }
    }
}

impl Sha256Writer {
    pub fn update(&mut self, bytes: &[u8]) -> Result<(), Diagnostic> {
        self.size_bytes = self
            .size_bytes
            .checked_add(u64::try_from(bytes.len()).map_err(|_| size_error())?)
            .ok_or_else(size_error)?;
        self.hasher.update(bytes);
        Ok(())
    }

    pub fn finish(self, relative_path: impl Into<PathBuf>) -> ArtifactIdentity {
        ArtifactIdentity {
            relative_path: relative_path.into(),
            size_bytes: self.size_bytes,
            digest: format!("{:x}", self.hasher.finalize()),
        }
    }
}

pub fn identify_file(
    path: &Path,
    relative_path: impl Into<PathBuf>,
) -> Result<ArtifactIdentity, Diagnostic> {
    let mut file = File::open(path).map_err(|error| io_error(path, "open", error))?;
    let mut writer = Sha256Writer::default();
    let mut buffer = [0_u8; 64 * 1024];
    loop {
        let read = file
            .read(&mut buffer)
            .map_err(|error| io_error(path, "read", error))?;
        if read == 0 {
            break;
        }
        writer.update(&buffer[..read])?;
    }
    Ok(writer.finish(relative_path))
}

pub fn verify_identity(path: &Path, expected: &ArtifactIdentity) -> Result<(), Diagnostic> {
    let observed = identify_file(path, expected.relative_path.clone())?;
    if observed.size_bytes == expected.size_bytes && observed.digest == expected.digest {
        return Ok(());
    }
    Err(Diagnostic::io(
        "GPS-IO-002",
        "El tamaño o SHA-256 de un artefacto no coincide con el valor esperado.",
    )
    .with("path", path.display().to_string())
    .with("expected_size", expected.size_bytes.to_string())
    .with("observed_size", observed.size_bytes.to_string())
    .with("expected_digest", &expected.digest)
    .with("observed_digest", observed.digest))
}

fn size_error() -> Diagnostic {
    Diagnostic::internal(
        "GPS-INTERNAL-001",
        "Una invariante interna considerada imposible se ha roto.",
    )
    .with("component", "integrity")
    .with("operation", "count_bytes")
    .with("detail", "artifact size exceeds u64")
}

fn io_error(path: &Path, operation: &str, error: io::Error) -> Diagnostic {
    Diagnostic::io("GPS-IO-001", "Una operación de lectura o escritura falló.")
        .with("path", path.display().to_string())
        .with("operation", operation)
        .with("detail", error.to_string())
}
