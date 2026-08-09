//! Núcleo del simulador de genomas y poblaciones.
//!
//! La biblioteca mantiene separadas las fronteras de dominio, aplicación y
//! adaptadores. Las fases posteriores añadirán comportamiento detrás de estas
//! fronteras sin filtrar tipos de dependencias externas al dominio.

pub mod adapters;
pub mod application;
pub mod domain;
