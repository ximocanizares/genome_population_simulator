//! Guarda previa de los lectores BGZF de los jueces (DT-030).
//!
//! `noodles-bgzf` 0.47 no acota el tamaño descomprimido que declara el pie de
//! cada bloque: uno que declara más de 65 536 bytes hace **entrar en pánico** al
//! lector («range end index 65537 out of range for slice of length 65536») en
//! vez de devolver un error. Un juez no puede abortar ante un fichero
//! manipulado, y subir de versión exige subir todo el grafo de `noodles`
//! (`noodles-csi` y `noodles-tabix` dependen de la 0.47; DTC-043 y DTC-048 fijan
//! las versiones exactas). Esta guarda recorre los bloques **antes** de
//! entregarlos al lector y rechaza el fichero.
//!
//! # Qué hace y qué no
//!
//! Recorre solo las cabeceras y los pies: lee 18 bytes de cabecera y los 4
//! últimos de cada bloque (el `ISIZE`) y salta el resto con un desplazamiento,
//! sin descomprimir nada. Usa el mismo encuadre que el lector —cabecera BGZF
//! de 18 bytes y `BSIZE` en sus dos últimos—, de modo que no se salta ningún
//! bloque que el lector fuera a abrir.
//!
//! **Solo rechaza un `ISIZE` imposible.** Una cabecera que no es BGZF o un
//! bloque más corto de lo que dice su cabecera **detienen el recorrido sin
//! rechazar**: es el lector quien las diagnostica como errores, y lo hace antes
//! de llegar a ningún bloque posterior. La suma de comprobación no se lee: la
//! guarda sigue adelante, y es el lector quien la verifica al descomprimir.

use std::io::{self, Read, Seek, SeekFrom};

/// Tamaño máximo descomprimido de un bloque BGZF (especificación SAM, § 4.1).
pub(crate) const BGZF_MAX_BLOCK_ISIZE: u32 = 1 << 16;

const HEADER_SIZE: usize = 18;
const TRAILER_SIZE: u64 = 8;
const MIN_BLOCK_SIZE: u64 = HEADER_SIZE as u64 + TRAILER_SIZE;

#[derive(Debug)]
pub(crate) enum BgzfBlockFault {
    Io(io::Error),
    /// El bloque que empieza en `offset` declara `declared` bytes descomprimidos.
    OversizedBlock {
        offset: u64,
        declared: u32,
    },
}

impl std::fmt::Display for BgzfBlockFault {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(error) => write!(formatter, "no se pudo recorrer el fichero: {error}"),
            Self::OversizedBlock { offset, declared } => write!(
                formatter,
                "el bloque BGZF que empieza en el byte {offset} declara {declared} bytes \
descomprimidos y el máximo es {BGZF_MAX_BLOCK_ISIZE}"
            ),
        }
    }
}

/// La cabecera fija de un bloque BGZF, tal y como la exige el lector.
fn is_bgzf_header(header: &[u8; HEADER_SIZE]) -> bool {
    header[0..2] == [0x1f, 0x8b]
        && header[2] == 0x08
        && header[3] == 0x04
        && header[10..12] == [0x06, 0x00]
        && header[12..14] == *b"BC"
        && header[14..16] == [0x02, 0x00]
}

/// Rechaza un fichero BGZF con algún bloque que declare más de 65 536 bytes
/// descomprimidos. Recorre desde el byte cero y **deja la posición en el
/// principio**, listo para entregárselo al lector.
pub(crate) fn reject_oversized_bgzf_blocks<R: Read + Seek>(
    source: &mut R,
) -> Result<(), BgzfBlockFault> {
    let result = scan(source);
    source
        .seek(SeekFrom::Start(0))
        .map_err(BgzfBlockFault::Io)?;
    result
}

fn scan<R: Read + Seek>(source: &mut R) -> Result<(), BgzfBlockFault> {
    let mut offset = source
        .seek(SeekFrom::Start(0))
        .map_err(BgzfBlockFault::Io)?;
    loop {
        let mut header = [0_u8; HEADER_SIZE];
        match source.read_exact(&mut header) {
            Ok(()) => {}
            Err(error) if error.kind() == io::ErrorKind::UnexpectedEof => return Ok(()),
            Err(error) => return Err(BgzfBlockFault::Io(error)),
        }
        if !is_bgzf_header(&header) {
            return Ok(());
        }
        let block_size = u64::from(u16::from_le_bytes([header[16], header[17]])) + 1;
        if block_size < MIN_BLOCK_SIZE {
            return Ok(());
        }
        // El `ISIZE` son los cuatro últimos bytes del bloque.
        let isize_position = offset + block_size - 4;
        source
            .seek(SeekFrom::Start(isize_position))
            .map_err(BgzfBlockFault::Io)?;
        let mut isize_bytes = [0_u8; 4];
        match source.read_exact(&mut isize_bytes) {
            Ok(()) => {}
            Err(error) if error.kind() == io::ErrorKind::UnexpectedEof => return Ok(()),
            Err(error) => return Err(BgzfBlockFault::Io(error)),
        }
        let declared = u32::from_le_bytes(isize_bytes);
        if declared > BGZF_MAX_BLOCK_ISIZE {
            return Err(BgzfBlockFault::OversizedBlock { offset, declared });
        }
        offset += block_size;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    /// Un bloque BGZF mínimo con el `ISIZE` que se pida, relleno hasta
    /// `block_size` bytes.
    fn block(block_size: usize, declared: u32) -> Vec<u8> {
        let mut bytes = vec![0_u8; block_size];
        bytes[..16].copy_from_slice(&[
            0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00, b'B', b'C', 0x02, 0x00,
        ]);
        bytes[16..18]
            .copy_from_slice(&u16::try_from(block_size - 1).expect("tamaño").to_le_bytes());
        bytes[block_size - 4..].copy_from_slice(&declared.to_le_bytes());
        bytes
    }

    #[test]
    fn accepts_blocks_up_to_the_maximum_and_rewinds() {
        let mut bytes = block(30, 10);
        bytes.extend(block(40, BGZF_MAX_BLOCK_ISIZE));
        bytes.extend(block(26, 0));
        let mut cursor = Cursor::new(bytes);
        reject_oversized_bgzf_blocks(&mut cursor).expect("un fichero válido se acepta");
        assert_eq!(cursor.position(), 0);
    }

    #[test]
    fn rejects_the_first_oversized_block_with_its_offset() {
        let mut bytes = block(30, 10);
        bytes.extend(block(40, BGZF_MAX_BLOCK_ISIZE + 1));
        let fault = reject_oversized_bgzf_blocks(&mut Cursor::new(bytes))
            .expect_err("un bloque de más de 65536 bytes se rechaza");
        assert!(
            matches!(
                fault,
                BgzfBlockFault::OversizedBlock {
                    offset: 30,
                    declared: 65_537
                }
            ),
            "{fault:?}"
        );
    }

    #[test]
    fn leaves_irregular_framing_to_the_reader() {
        // Una cabecera que no es BGZF, un bloque truncado y un fichero vacío no
        // son un `ISIZE` imposible: los diagnostica el lector.
        for bytes in [Vec::new(), vec![0_u8; 18], block(30, 10)[..25].to_vec()] {
            reject_oversized_bgzf_blocks(&mut Cursor::new(bytes)).expect("no es de esta guarda");
        }
    }

    /// Una cabecera BGZF válida cuyo `BSIZE` declara un bloque más corto que
    /// cabecera más pie (26 bytes). Sin la comprobación del tamaño mínimo, la
    /// posición del `ISIZE` (`offset + block_size - 4`) retrocede o desborda:
    /// con `BSIZE` 0 y desplazamiento 0 es una resta que da la vuelta.
    #[test]
    fn leaves_a_block_shorter_than_its_header_and_trailer_to_the_reader() {
        for block_size in [1_usize, 5, 25] {
            let mut bytes = block(30, 10)[..18].to_vec();
            bytes[16..18]
                .copy_from_slice(&u16::try_from(block_size - 1).expect("tamaño").to_le_bytes());
            bytes.extend_from_slice(&[0_u8; 16]);
            reject_oversized_bgzf_blocks(&mut Cursor::new(bytes))
                .unwrap_or_else(|fault| panic!("BSIZE {}: {fault:?}", block_size - 1));
        }
    }
}
