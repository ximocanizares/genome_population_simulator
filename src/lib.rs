//! Núcleo del simulador de genomas y poblaciones.
//!
//! La biblioteca mantiene separadas las fronteras de dominio, aplicación y
//! adaptadores. Las fases posteriores añadirán comportamiento detrás de estas
//! fronteras sin filtrar tipos de dependencias externas al dominio.

pub mod adapters;
pub mod application;
pub mod configuration;
pub mod diagnostics;
pub mod domain;
pub mod execution;
pub mod identity;
pub mod integrity;
pub mod numeric;
pub mod publication;
pub mod randomness;
pub mod results;
pub mod rules;
