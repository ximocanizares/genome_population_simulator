//! Núcleo del simulador de genomas y poblaciones.
//!
//! ## `unsafe` prohibido en la biblioteca, y por qué la prohibición vive aquí
//!
//! `#![forbid(unsafe_code)]` cubre **todo** `src/`: ningún módulo de la
//! biblioteca publicada puede escribir `unsafe`, y un intento de hacerlo es
//! error de compilación, no un aviso. La prohibición estuvo en
//! `[lints.rust]` de `Cargo.toml` hasta la revisión de F5.4.1; desde allí
//! alcanzaba también a los objetivos de prueba, donde `forbid` no puede
//! levantarse, y eso impedía instalar el `#[global_allocator]` contador que
//! mide el pico real del escritor de artefactos F5 contra el presupuesto que
//! la admisión sella. Mover el atributo aquí deja la garantía de la
//! biblioteca **idéntica** y libera únicamente a `tests/`. No lo devuelvas a
//! `Cargo.toml`: no se relajó nada, se acotó al lugar que debía cubrir.
//!
//! La biblioteca mantiene separadas las fronteras de dominio, aplicación y
//! adaptadores. Las fases posteriores añadirán comportamiento detrás de estas
//! fronteras sin filtrar tipos de dependencias externas al dominio.

#![forbid(unsafe_code)]

pub mod adapters;
pub mod application;
pub mod biological_contamination;
pub mod chimera;
pub mod configuration;
pub mod coordinates;
pub mod demultiplexing;
pub mod diagnostics;
pub mod domain;
pub mod error_quality;
pub mod event_manifest;
pub mod exact_decimal;
pub mod execution;
pub mod f5_admission;
pub mod f5_artifacts;
pub mod f5_closure;
pub mod f5_founders;
pub mod f5_haplotypes;
pub mod f5_layout;
pub mod f5_ld;
pub mod f5_marginals;
pub mod f5_oracle;
pub mod f5_preflight;
pub mod f5_profile;
pub mod f5_publication;
pub mod f5_search;
pub mod fastq_output;
pub mod fastq_validation;
pub mod haplotype;
pub mod identity;
pub mod indel_application;
pub mod initial_genotypes;
pub mod initial_individual_integration;
pub mod initial_individual_materialization;
pub mod integrity;
pub mod materialization;
pub mod paired_end;
pub mod population_snvs;
pub mod publication;
pub mod randomness;
pub mod reference;
pub mod reference_fasta;
pub mod reference_generation;
pub mod reference_reduction;
pub mod reference_reduction_io;
pub mod reference_tract_insertion;
pub mod repeat_genotypes;
pub mod results;
pub mod segments;
pub mod small_event_composition;
pub mod snp_application;
pub mod tandem_repeat;
pub mod technical_identity;
pub mod technical_profile;
pub mod truth_projection;
pub mod truth_schema;
pub mod truth_validation;
pub mod vcf_contract;
pub mod vcf_grouping;
pub mod vcf_projection;
pub mod vcf_publication;
pub mod vcf_representation;
pub mod vcf_semantic_validation;
pub mod wgs_sampling;
