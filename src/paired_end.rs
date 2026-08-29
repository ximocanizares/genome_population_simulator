//! Fragmentos y mates *paired-end* transitorios de F4.1.
//!
//! Esta frontera no publica FASTQ ni decide cuántas parejas corresponden a una
//! cobertura. Recibe un número ya decidido por su llamador (F4.2 será el único
//! que podrá derivarlo) y materializa una pareja completa desde una plantilla
//! y un fragmento único.

use crate::diagnostics::Diagnostic;
use crate::haplotype::HaplotypeCopyId;
use crate::technical_identity::{
    PreHopTechnicalPairIdentity, TechnicalMate, TechnicalTemplateIdentity,
};

/// Fragmento biológico semiabierto antes de orientar sus dos mates.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PairedEndFragment {
    template: TechnicalTemplateIdentity,
    haplotype_copy_id: HaplotypeCopyId,
    start: u64,
    end: u64,
    sequence: Vec<u8>,
}

impl PairedEndFragment {
    pub fn new(
        template: TechnicalTemplateIdentity,
        haplotype_copy_id: HaplotypeCopyId,
        start: u64,
        sequence: Vec<u8>,
    ) -> Result<Self, Diagnostic> {
        if sequence.is_empty() {
            return Err(pair_error("fragment", "el inserto no puede estar vacío"));
        }
        if !sequence
            .iter()
            .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T' | b'N'))
        {
            return Err(pair_error(
                "fragment",
                "el inserto contiene una base no canónica",
            ));
        }
        let length = u64::try_from(sequence.len())
            .map_err(|_| pair_error("fragment", "la longitud del inserto desborda u64"))?;
        let end = start
            .checked_add(length)
            .ok_or_else(|| pair_error("fragment", "el extremo derecho del inserto desborda u64"))?;
        Ok(Self {
            template,
            haplotype_copy_id,
            start,
            end,
            sequence,
        })
    }

    #[must_use]
    pub const fn template(&self) -> &TechnicalTemplateIdentity {
        &self.template
    }

    #[must_use]
    pub const fn haplotype_copy_id(&self) -> &HaplotypeCopyId {
        &self.haplotype_copy_id
    }

    #[must_use]
    pub const fn start(&self) -> u64 {
        self.start
    }

    #[must_use]
    pub const fn end(&self) -> u64 {
        self.end
    }

    #[must_use]
    pub fn sequence(&self) -> &[u8] {
        &self.sequence
    }
}

/// Un mate PE verdadero. La calidad es transitoria y solo expresa la regla
/// determinista vigente antes de F4.2: A/C/G/T con Phred 40 y N con Phred 0.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PairedEndMate {
    mate: TechnicalMate,
    sequence: Vec<u8>,
    phred: Vec<u8>,
}

impl PairedEndMate {
    #[must_use]
    pub const fn mate(&self) -> TechnicalMate {
        self.mate
    }

    #[must_use]
    pub fn sequence(&self) -> &[u8] {
        &self.sequence
    }

    #[must_use]
    pub fn phred(&self) -> &[u8] {
        &self.phred
    }
}

/// Pareja normal FR. No hay representación para un mate huérfano.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PairedEndPair {
    pre_hop_identity: PreHopTechnicalPairIdentity,
    fragment: PairedEndFragment,
    r1: PairedEndMate,
    r2: PairedEndMate,
}

impl PairedEndPair {
    #[must_use]
    pub const fn pre_hop_identity(&self) -> &PreHopTechnicalPairIdentity {
        &self.pre_hop_identity
    }

    #[must_use]
    pub const fn fragment(&self) -> &PairedEndFragment {
        &self.fragment
    }

    #[must_use]
    pub const fn r1(&self) -> &PairedEndMate {
        &self.r1
    }

    #[must_use]
    pub const fn r2(&self) -> &PairedEndMate {
        &self.r2
    }
}

/// Emite R1 desde el extremo izquierdo y R2 como complemento inverso del
/// extremo derecho. El solapamiento se conserva: ambos cortes se hacen sobre
/// el mismo inserto íntegro y no se fusionan.
pub fn emit_paired_end(
    fragment: PairedEndFragment,
    pre_hop_identity: PreHopTechnicalPairIdentity,
    common_read_length: u16,
) -> Result<PairedEndPair, Diagnostic> {
    if pre_hop_identity.template() != fragment.template() {
        return Err(pair_error(
            "pre_hop_identity",
            "la identidad PH no corresponde a la plantilla del fragmento",
        ));
    }
    let read_length = usize::from(common_read_length);
    if read_length == 0 {
        return Err(pair_error(
            "read_length",
            "la longitud común debe ser al menos una base",
        ));
    }
    if fragment.sequence.len() < read_length {
        return Err(pair_error(
            "read_length",
            "la longitud común no puede superar el inserto",
        ));
    }
    let r1_sequence = fragment.sequence[..read_length].to_vec();
    let mut r2_sequence = fragment.sequence[fragment.sequence.len() - read_length..].to_vec();
    reverse_complement(&mut r2_sequence);
    let r1 = PairedEndMate {
        mate: TechnicalMate::R1,
        phred: phred_for(&r1_sequence),
        sequence: r1_sequence,
    };
    let r2 = PairedEndMate {
        mate: TechnicalMate::R2,
        phred: phred_for(&r2_sequence),
        sequence: r2_sequence,
    };
    Ok(PairedEndPair {
        pre_hop_identity,
        fragment,
        r1,
        r2,
    })
}

fn phred_for(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .map(|base| if *base == b'N' { 0 } else { 40 })
        .collect()
}

fn reverse_complement(sequence: &mut [u8]) {
    sequence.reverse();
    for base in sequence {
        *base = match *base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'N' => b'N',
            _ => unreachable!("PairedEndFragment validates its base alphabet"),
        };
    }
}

fn pair_error(field: &'static str, detail: &'static str) -> Diagnostic {
    Diagnostic::internal_invariant("paired_end", field).with("detail", detail)
}
