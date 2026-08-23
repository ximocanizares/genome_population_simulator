//! Núcleo del simulador de genomas y poblaciones.
//!
//! La biblioteca mantiene separadas las fronteras de dominio, aplicación y
//! adaptadores. Las fases posteriores añadirán comportamiento detrás de estas
//! fronteras sin filtrar tipos de dependencias externas al dominio.

pub mod adapters;
pub mod application;
pub mod configuration;
pub mod coordinates;
pub mod diagnostics;
pub mod domain;
pub mod exact_decimal;
pub mod execution;
pub mod haplotype;
pub mod identity;
pub mod indel_application;
pub mod initial_genotypes;
pub mod integrity;
pub mod materialization;
pub mod population_snvs;
pub mod publication;
pub mod randomness;
pub mod reference;
pub mod reference_fasta;
pub mod reference_generation;
pub mod reference_reduction;
pub mod reference_reduction_io;
pub mod results;
pub mod segments;
pub mod small_event_composition;
pub mod snp_application;
