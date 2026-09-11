//! Núcleo del simulador de genomas y poblaciones.
//!
//! La biblioteca mantiene separadas las fronteras de dominio, aplicación y
//! adaptadores. Las fases posteriores añadirán comportamiento detrás de estas
//! fronteras sin filtrar tipos de dependencias externas al dominio.

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
pub mod f5_haplotypes;
pub mod f5_layout;
pub mod f5_ld;
pub mod f5_marginals;
pub mod f5_preflight;
pub mod f5_profile;
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
pub mod results;
pub mod segments;
pub mod small_event_composition;
pub mod snp_application;
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
