//! Sorteo de los alelos de repetición por copia haplotípica y su aplicación.
//!
//! # Qué decide este módulo
//!
//! El usuario declara, por locus, un puñado de alelos con sus **frecuencias**;
//! aquí se sortea cuál le toca a cada copia de cada individuo. Es la decisión
//! S-04 del plan de repeticiones: declarar a mano cien individuos por dos
//! copias y por locus no es generar, es transcribir.
//!
//! # El ámbito de azar, y por qué sus componentes son los que son
//!
//! Cada sorteo tiene subsemilla propia, derivada de **identidades estables**:
//! escenario, referencia, población, individuo, **índice de copia** e
//! **identificador del locus**. No entran ni la posición del locus en la
//! referencia ni el orden en que individuos o loci aparecen en la
//! configuración. La consecuencia observable es la que importa: insertar un
//! individuo al principio de la lista, reordenarla, o añadir un locus en un
//! cromosoma que ordene antes que los existentes **no cambia ni un genotipo ya
//! sorteado**.
//!
//! # Cuándo no se consume azar
//!
//! Cuando un locus tiene un solo alelo con peso positivo, todas las copias lo
//! reciben y **no se pide ni una palabra al generador**. No es un atajo de
//! rendimiento: es lo que hace que añadir un locus monomórfico a un escenario
//! no desplace la secuencia de azar de nada. El cortocircuito vive en
//! `categorical_f64`, que ya lo tenía.

use crate::application::{EffectiveConfig, EffectiveRepeatLocus};
use crate::diagnostics::Diagnostic;
use crate::haplotype::{HaplotypeCopy, HaplotypeStateError, Individual};
use crate::identity::{EventId, IndividualId, PopulationId};
use crate::randomness::{
    CATEGORICAL_F64_ALGORITHM, RandomWordSource, RootSeed, RuntimeSubseedScope,
    categorical_f64_with_source, derive_subseed,
};
use crate::reference::Reference;
use crate::segments::{
    CurrentIndelEvent, CurrentRepeatAllelePiece, ModificationProvenance, SegmentError,
};

/// Ámbito de subsemilla del sorteo, tal y como lo publica el catálogo.
pub const REPEAT_ALLELE_SCOPE: &str = "initial_repeat_allele";
/// Algoritmo normativo del sorteo (DTC-072). No estrena primitiva.
pub const REPEAT_ALLELE_ALGORITHM: &str = CATEGORICAL_F64_ALGORITHM;
/// Escenario de los ámbitos de azar de la población inicial, igual que el de
/// la composición y la fase genotípicas.
const REPEAT_ALLELE_SCENARIO: &str = "direct";

/// El alelo que le tocó a una copia concreta de un locus concreto.
///
/// Es lo que permite publicar la **frecuencia realizada** sin volver a
/// recorrer las copias: el recuento de esta lista es el numerador.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RepeatAlleleDraw {
    locus_id: String,
    individual_id: IndividualId,
    copy_index: u8,
    allele_index: usize,
}

impl RepeatAlleleDraw {
    #[must_use]
    pub fn locus_id(&self) -> &str {
        &self.locus_id
    }

    #[must_use]
    pub const fn individual_id(&self) -> &IndividualId {
        &self.individual_id
    }

    #[must_use]
    pub const fn copy_index(&self) -> u8 {
        self.copy_index
    }

    /// Índice del alelo dentro de la lista publicada del locus.
    #[must_use]
    pub const fn allele_index(&self) -> usize {
        self.allele_index
    }
}

/// Individuo con sus alelos de repetición ya aplicados, y lo que le tocó.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AppliedRepeatAlleles {
    individual: Individual,
    draws: Vec<RepeatAlleleDraw>,
}

impl AppliedRepeatAlleles {
    #[must_use]
    pub const fn individual(&self) -> &Individual {
        &self.individual
    }

    #[must_use]
    pub fn draws(&self) -> &[RepeatAlleleDraw] {
        &self.draws
    }

    #[must_use]
    pub fn into_parts(self) -> (Individual, Vec<RepeatAlleleDraw>) {
        (self.individual, self.draws)
    }
}

/// Sortea el alelo de un locus para una copia, sin tocar ningún estado.
pub fn draw_repeat_allele_index(
    root_seed: RootSeed,
    reference_id: &str,
    population_id: &PopulationId,
    individual_id: &IndividualId,
    copy_index: u8,
    locus: &EffectiveRepeatLocus,
) -> Result<usize, Diagnostic> {
    let mut rng = derive_subseed(
        root_seed,
        RuntimeSubseedScope::InitialRepeatAllele {
            scenario_id: REPEAT_ALLELE_SCENARIO,
            reference_id,
            population_id: population_id.as_str(),
            individual_id: individual_id.as_str(),
            copy_index,
            locus_id: locus.id(),
        },
    )?
    .local_rng();
    draw_with_source(&mut rng, locus)
}

/// El mismo sorteo contra una fuente de palabras cualquiera.
///
/// Existe para que una prueba pueda **contar** cuántas palabras consume: con un
/// solo alelo de peso positivo tienen que ser cero, y sin poder contarlas esa
/// afirmación no sería comprobable.
pub(crate) fn draw_with_source(
    rng: &mut impl RandomWordSource,
    locus: &EffectiveRepeatLocus,
) -> Result<usize, Diagnostic> {
    // Los pesos son los **normalizados** por la suma compensada de DTC-049,
    // que es la misma que consume el sorteo categórico.
    let weights: Vec<f64> = locus
        .alleles()
        .iter()
        .map(crate::application::EffectiveRepeatAllele::frequency)
        .collect();
    categorical_f64_with_source(rng, &weights).map_err(|error| {
        Diagnostic::internal(
            "GPS-INTERNAL-001",
            "Una invariante interna considerada imposible se ha roto.",
        )
        .with("component", "repeat_allele_draw")
        .with("operation", REPEAT_ALLELE_ALGORITHM)
        .with("detail", format!("locus {}: {error:?}", locus.id()))
    })
}

/// Sortea y aplica a un individuo el alelo de cada locus declarado.
///
/// # Lo que **no** se aplica, y por qué importa
///
/// Cuando el alelo sorteado es byte a byte el tracto que ya hay en la
/// referencia, no se aplica nada: no hay variante que publicar. Es lo que hace
/// que un escenario cuyo único alelo es el de la referencia produzca una verdad
/// idéntica a la de antes de que existieran los microsatélites.
pub fn apply_repeat_alleles(
    individual: &Individual,
    reference: &Reference,
    config: &EffectiveConfig,
    population_id: &PopulationId,
) -> Result<AppliedRepeatAlleles, Diagnostic> {
    let mut updated = individual.clone();
    let mut draws = Vec::new();
    if config.repeat_loci().is_empty() {
        return Ok(AppliedRepeatAlleles {
            individual: updated,
            draws,
        });
    }
    let root_seed = config.reproducibility().root_seed();
    for locus in config.repeat_loci() {
        let interval = locus.declaration().interval();
        let sequence_id = interval.sequence_id();
        let reference_tract = reference_tract(reference, locus)?;
        let event_id = repeat_event_id(locus)?;
        let copies = individual
            .copies_for(sequence_id)
            .ok_or_else(|| {
                repeat_genotype_error(
                    locus.id(),
                    format!("el individuo no tiene copias de {sequence_id}"),
                )
            })?
            .iter()
            .map(|copy| (copy.phase(), copy.phase().copy_index()))
            .collect::<Vec<_>>();
        for (phase, copy_index) in copies {
            let allele_index = draw_repeat_allele_index(
                root_seed,
                config.reference_id(),
                population_id,
                individual.id(),
                copy_index,
                locus,
            )?;
            let allele = locus.alleles().get(allele_index).ok_or_else(|| {
                repeat_genotype_error(locus.id(), "el sorteo devolvió un alelo inexistente")
            })?;
            draws.push(RepeatAlleleDraw {
                locus_id: locus.id().to_owned(),
                individual_id: individual.id().clone(),
                copy_index,
                allele_index,
            });
            let copy = updated
                .copies_for(sequence_id)
                .and_then(|copies| copies.iter().find(|copy| copy.phase() == phase))
                .ok_or_else(|| {
                    repeat_genotype_error(locus.id(), "la copia objetivo desapareció")
                })?;
            // **La condición se comprueba siempre, se aplique o no el alelo.**
            // Cuando el alelo sorteado coincide con el tracto de la referencia
            // no hay nada que escribir, pero la verdad va a publicar que esta
            // copia conserva ese tracto, y eso sería falso si algo lo hubiera
            // tocado. Comprobarlo solo al aplicar hacía que el rechazo
            // dependiera del sorteo y, con él, de la semilla.
            copy.segments()
                .intact_repeat_interval(interval)
                .map_err(|error| invaded_locus_error(locus, copy_index, copy, error))?;
            if allele.tract().bases() == reference_tract.as_slice() {
                continue;
            }
            let segments = copy
                .segments()
                .with_current_repeat_allele(CurrentRepeatAllelePiece::new(
                    locus.declaration().id().clone(),
                    interval.clone(),
                    allele.tract().clone(),
                    ModificationProvenance::random(event_id.clone()),
                ))
                .map_err(|error| invaded_locus_error(locus, copy_index, copy, error))?;
            updated = updated
                .with_copy_segments(sequence_id, phase, segments)
                .map_err(|error| {
                    repeat_genotype_error(
                        locus.id(),
                        format!(
                            "no se pudo sustituir la copia {copy_index}: {}",
                            describe_haplotype_error(error)
                        ),
                    )
                })?;
        }
    }
    Ok(AppliedRepeatAlleles {
        individual: updated,
        draws,
    })
}

/// Las bases que el locus ocupa hoy en la referencia.
pub fn reference_tract(
    reference: &Reference,
    locus: &EffectiveRepeatLocus,
) -> Result<Vec<u8>, Diagnostic> {
    let interval = locus.declaration().interval();
    let sequence = reference
        .sequences()
        .iter()
        .find(|sequence| sequence.id() == interval.sequence_id())
        .ok_or_else(|| {
            repeat_genotype_error(
                locus.id(),
                format!("la referencia no contiene {}", interval.sequence_id()),
            )
        })?;
    let bases = sequence.sequence().as_bytes().map_err(|error| {
        repeat_genotype_error(
            locus.id(),
            format!("las bases de la referencia no están materializadas: {error:?}"),
        )
    })?;
    let start = usize::try_from(interval.start())
        .map_err(|_| repeat_genotype_error(locus.id(), "el inicio del intervalo desborda"))?;
    let end = usize::try_from(interval.end())
        .map_err(|_| repeat_genotype_error(locus.id(), "el final del intervalo desborda"))?;
    bases
        .get(start..end)
        .map(<[u8]>::to_vec)
        .ok_or_else(|| repeat_genotype_error(locus.id(), "el intervalo rebasa la secuencia"))
}

/// El evento que coloca un alelo de repetición **es el propio locus**.
///
/// No se inventa un identificador nuevo ni se le pone prefijo: el prefijo
/// podría rebasar los 64 caracteres que admite un identificador aportado, y un
/// identificador automático nuevo contradiría la decisión S-11, que dice
/// expresamente que este trabajo no estrena clase de identificador.
fn repeat_event_id(locus: &EffectiveRepeatLocus) -> Result<EventId, Diagnostic> {
    EventId::parse(locus.id(), "repeat_locus")
}

/// Traduce a castellano llano el rechazo del estado cuando el tramo del locus
/// no llega intacto.
///
/// **Es la comprobación exacta de la decisión S-32 para las inserciones y las
/// deleciones**, no un remate: la frontera de configuración solo juzga lo que
/// se declara en coordenadas de la referencia —cambios de base y filas de la
/// matriz—, y un indel no lo es. Aquí el intervalo se busca sobre el estado ya
/// compuesto, de modo que no hay deriva de coordenadas posible.
///
/// El mensaje **enumera los cambios que tocaron esta copia**. No puede decir
/// cuál de ellos invadió el tracto —lo que el estado sabe es que el tramo ya no
/// llega entero—, pero dar la lista es lo que convierte el rechazo en algo
/// accionable en vez de en un callejón.
fn invaded_locus_error(
    locus: &EffectiveRepeatLocus,
    copy_index: u8,
    copy: &HaplotypeCopy,
    error: SegmentError,
) -> Diagnostic {
    let cause = match error {
        SegmentError::RepeatIntervalNotIntact => {
            "una inserción o una deleción anterior ha partido, borrado o sustituido ese tramo"
        }
        SegmentError::RepeatIntervalAlreadyModified => {
            "un cambio de base anterior cae dentro de ese tramo"
        }
        // Las demás variantes no las produce esta operación; se nombran en vez
        // de absorberse con un comodín para que una nueva obligue a decidir.
        SegmentError::RepeatAlleleIsNotDivisible => {
            "una inserción o una deleción posterior caería dentro del tracto"
        }
        _ => "el estado de la copia no admite el tracto en ese tramo",
    };
    // Solo las inserciones y las deleciones: son la única clase que puede
    // partir, borrar o desplazar el tramo. Un cambio de base dentro del
    // intervalo lo rechaza antes la frontera de configuración, que sí puede
    // nombrarlo, de modo que listarlos aquí sería ruido — y la matriz aporta
    // identificadores derivados que el usuario no escribió.
    let culprits = copy
        .segments()
        .small_event_history()
        .iter()
        .filter(|event| {
            matches!(
                event,
                CurrentIndelEvent::Insertion { .. } | CurrentIndelEvent::Deletion { .. }
            )
        })
        .map(|event| event.event_id().as_str())
        .collect::<Vec<_>>();
    let listed = if culprits.is_empty() {
        String::new()
    } else {
        format!(
            " Las inserciones y deleciones declaradas sobre esta copia son: {}.",
            culprits.join(", ")
        )
    };
    let interval = locus.declaration().interval();
    repeat_genotype_error(
        locus.id(),
        format!(
            "el tracto declarado ({}:{} a {}) no llega intacto a la copia {copy_index}: {cause}. Un alelo de repetición describe el tracto que hay en la referencia, de modo que ninguna otra modificación puede tocarlo.{listed}",
            interval.sequence_id(),
            interval.start() + 1,
            interval.end()
        ),
    )
}

/// Las tres formas en que una copia puede no admitir una cadena nueva, dichas
/// en castellano en vez de con el nombre de su variante de Rust.
fn describe_haplotype_error(error: HaplotypeStateError) -> &'static str {
    match error {
        HaplotypeStateError::MissingSequence => "el individuo no tiene esa secuencia",
        HaplotypeStateError::MissingCopyPhase => "el individuo no tiene esa copia",
        HaplotypeStateError::SegmentSourceMismatch => "los segmentos son de otra referencia",
        HaplotypeStateError::SegmentSourceLengthMismatch { .. } => {
            "los segmentos declaran otra longitud de secuencia"
        }
        HaplotypeStateError::PhaseOutOfRange { .. } => "la fase queda fuera del rango de ploidía",
    }
}

fn repeat_genotype_error(locus_id: &str, detail: impl Into<String>) -> Diagnostic {
    Diagnostic::configuration(
        "GPS-CONF-008",
        "La configuración contiene una combinación de campos incompatible.",
    )
    .with("path", "effective-config")
    .with("fields", "reference.repeat_loci")
    .with("detail", format!("locus {locus_id}: {}", detail.into()))
}

#[cfg(test)]
mod tests {
    use super::draw_with_source;
    use crate::adapters::{admit_toml_header, parse_declared_config};
    use crate::application::{EffectiveConfig, build_effective_config};
    use crate::configuration::ExpectedConfigKind;
    use crate::randomness::RandomWordSource;
    use std::path::Path;

    /// Fuente que **cuenta** las palabras que se le piden.
    ///
    /// Sin ella la afirmación «con un solo alelo no se consume aleatoriedad»
    /// no sería comprobable: el resultado del sorteo es el mismo se consuma o
    /// no, de modo que solo el contador la distingue.
    struct Counting {
        consumed: usize,
    }

    impl RandomWordSource for Counting {
        fn next_word(&mut self) -> u64 {
            self.consumed += 1;
            0
        }
    }

    fn config(alleles: &str) -> EffectiveConfig {
        let source = format!(
            r#"config_kind = "general"
schema_version = "4"
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"

[reference]
kind = "simulated"
id = "ref-main"

[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 200

[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = 21
end = 80
motif = "CA"
{alleles}
"#
        );
        let admitted = admit_toml_header(
            Path::new("repeat-genotype-draw-unit.toml"),
            source.into_bytes(),
            ExpectedConfigKind::General,
        )
        .expect("header");
        build_effective_config(&parse_declared_config(admitted).expect("declared"))
            .expect("effective")
    }

    /// **Con un solo alelo de peso positivo no se pide ni una palabra.**
    ///
    /// La pareja es lo que la hace discriminante: el mismo contador sobre un
    /// locus de dos alelos **sí** consume, de modo que un cero no puede venir
    /// de que la fuente no se use nunca.
    #[test]
    fn a_single_positive_weight_consumes_no_randomness() {
        let single = config(
            r#"
[[reference.repeat_loci.alleles]]
units = 25
frequency = 1.0
"#,
        );
        let mut counter = Counting { consumed: 0 };
        assert_eq!(
            draw_with_source(&mut counter, &single.repeat_loci()[0]),
            Ok(0)
        );
        assert_eq!(counter.consumed, 0);

        let pair = config(
            r#"
[[reference.repeat_loci.alleles]]
units = 25
frequency = 0.5

[[reference.repeat_loci.alleles]]
units = 20
frequency = 0.5
"#,
        );
        let mut counter = Counting { consumed: 0 };
        draw_with_source(&mut counter, &pair.repeat_loci()[0]).expect("sorteo de dos alelos");
        assert_eq!(
            counter.consumed, 1,
            "con dos alelos el sorteo sí consume, de modo que el cero de arriba significa algo"
        );
    }
}
