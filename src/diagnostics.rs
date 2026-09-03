//! Diagnósticos públicos estables independientes de errores de crates.

use std::collections::BTreeMap;
use std::sync::OnceLock;

/// Diagnóstico estructurado que alimenta todos los renderizados públicos.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Diagnostic {
    code: &'static str,
    category: &'static str,
    message: &'static str,
    exit_code: u8,
    context: BTreeMap<String, String>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct CatalogContract {
    code: &'static str,
    category: &'static str,
    meaning: &'static str,
    exit_code: u8,
}

impl Diagnostic {
    /// Construye el diagnóstico estable para una invariante interna rota en
    /// una frontera pública, sin exponer el constructor general del catálogo.
    #[must_use]
    pub fn internal_invariant(component: impl Into<String>, operation: impl Into<String>) -> Self {
        Self::internal(
            "GPS-INTERNAL-001",
            "Una invariante interna considerada imposible se ha roto.",
        )
        .with("component", component.into())
        .with("operation", operation.into())
    }

    pub(crate) fn configuration(code: &'static str, _message: &'static str) -> Self {
        Self::from_catalog(code)
    }

    /// Factoría F5: no reutiliza el contrato histórico GPS-CONF-008 y siempre
    /// materializa el contexto público obligatorio del catálogo F5.
    #[must_use]
    pub fn f5_configuration(
        path: impl Into<String>,
        fields: impl Into<String>,
        detail: impl Into<String>,
    ) -> Self {
        Self::from_catalog("GPS-CONF-101")
            .with("path", path.into())
            .with("fields", fields.into())
            .with("detail", detail.into())
    }

    #[must_use]
    pub fn into_f5_configuration(self) -> Self {
        let mut context = self.context;
        let path = context
            .remove("path")
            .unwrap_or_else(|| "<unknown>".to_owned());
        let fields = context
            .remove("fields")
            .or_else(|| context.remove("field"))
            .unwrap_or_else(|| "document".to_owned());
        let detail = context
            .remove("detail")
            .unwrap_or_else(|| self.message.to_owned());
        let mut diagnostic = Self::f5_configuration(path, fields, detail);
        for (key, value) in context {
            diagnostic = diagnostic.with(key, value);
        }
        diagnostic
    }

    pub(crate) fn input(code: &'static str, _message: &'static str) -> Self {
        Self::from_catalog(code)
    }

    pub(crate) fn format(code: &'static str, _message: &'static str) -> Self {
        Self::from_catalog(code)
    }

    pub(crate) fn resource(code: &'static str, _message: &'static str) -> Self {
        Self::from_catalog(code)
    }

    pub(crate) fn domain(code: &'static str, _message: &'static str) -> Self {
        Self::from_catalog(code)
    }

    pub(crate) fn io(code: &'static str, _message: &'static str) -> Self {
        Self::from_catalog(code)
    }

    pub(crate) fn internal(code: &'static str, _message: &'static str) -> Self {
        Self::from_catalog(code)
    }

    fn from_catalog(code: &'static str) -> Self {
        let contract = catalog_contract(code)
            .unwrap_or_else(|| panic!("diagnostic code absent from runtime catalog: {code}"));
        Self {
            code: contract.code,
            category: contract.category,
            message: contract.meaning,
            exit_code: contract.exit_code,
            context: BTreeMap::new(),
        }
    }

    #[must_use]
    pub fn with(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.context.insert(key.into(), value.into());
        self
    }

    #[must_use]
    pub const fn code(&self) -> &'static str {
        self.code
    }

    #[must_use]
    pub const fn category(&self) -> &'static str {
        self.category
    }

    #[must_use]
    pub const fn message(&self) -> &'static str {
        self.message
    }

    #[must_use]
    pub const fn exit_code(&self) -> u8 {
        self.exit_code
    }

    #[must_use]
    pub fn context(&self) -> &BTreeMap<String, String> {
        &self.context
    }

    #[must_use]
    pub fn render_human(&self) -> String {
        let context = self
            .context
            .iter()
            .map(|(key, value)| format!("{key}={value:?}"))
            .collect::<Vec<_>>()
            .join(" ");
        if context.is_empty() {
            format!("{} [{}] {}", self.code, self.category, self.message)
        } else {
            format!(
                "{} [{}] {} {}",
                self.code, self.category, self.message, context
            )
        }
    }

    /// Documento RFC 8259 único generado desde el mismo objeto que el modo
    /// humano. No depende de un modelo de error de terceros.
    #[must_use]
    pub fn render_json_document(&self) -> String {
        let context = self
            .context
            .iter()
            .map(|(key, value)| format!("\"{}\":\"{}\"", escape_json(key), escape_json(value)))
            .collect::<Vec<_>>()
            .join(",");
        format!(
            concat!(
                "{{\"diagnostic_kind\":\"gps_diagnostics\",",
                "\"schema_version\":\"1\",\"exit_code\":{},",
                "\"diagnostics\":[{{\"code\":\"{}\",\"category\":\"{}\",",
                "\"severity\":\"error\",\"message\":\"{}\",\"context\":{{{}}}}}]}}"
            ),
            self.exit_code,
            escape_json(self.code),
            escape_json(self.category),
            escape_json(self.message),
            context
        )
    }
}

fn catalog_contract(code: &'static str) -> Option<CatalogContract> {
    diagnostic_catalog().get(code).copied()
}

/// Catálogo activo para los diagnósticos que producen las rutas vigentes.
/// Los documentos ya publicados no se validan contra él: su versión fija el
/// catálogo histórico aplicable.
fn diagnostic_catalog() -> &'static BTreeMap<&'static str, CatalogContract> {
    static CATALOG: OnceLock<BTreeMap<&'static str, CatalogContract>> = OnceLock::new();
    CATALOG.get_or_init(|| {
        let source = include_str!("../catalogs/diagnostics-v4.tsv");
        assert!(
            source
                .lines()
                .any(|line| line == "# diagnostic_catalog_version=4"),
            "diagnostics-v4.tsv does not declare version 4"
        );
        let mut catalog = BTreeMap::new();
        for line in source
            .lines()
            .filter(|line| !line.is_empty() && !line.starts_with('#'))
        {
            let fields = line.split('|').collect::<Vec<_>>();
            assert_eq!(fields.len(), 6, "invalid diagnostics-v4.tsv row: {line}");
            let code = fields[0];
            let contract = CatalogContract {
                code,
                category: fields[1],
                meaning: fields[2],
                exit_code: fields[4]
                    .parse::<u8>()
                    .expect("invalid diagnostics-v4.tsv exit code"),
            };
            assert!(
                catalog.insert(code, contract).is_none(),
                "duplicate diagnostics-v4.tsv code: {code}"
            );
        }
        catalog
    })
}

/// Comprueba un diagnóstico contra el catálogo congelado por un documento
/// versionado. Una versión posterior puede añadir códigos, pero no los vuelve
/// válidos retroactivamente en un manifiesto anterior.
pub(crate) fn diagnostic_matches_catalog_version(
    catalog_version: &str,
    code: &str,
    message: &str,
) -> bool {
    diagnostic_catalog_for_version(catalog_version)
        .is_some_and(|catalog| catalog.get(code).is_some_and(|meaning| *meaning == message))
}

fn diagnostic_catalog_for_version(
    catalog_version: &str,
) -> Option<BTreeMap<&'static str, &'static str>> {
    let source = match catalog_version {
        "3" => include_str!("../catalogs/diagnostics-v3.tsv"),
        "4" => include_str!("../catalogs/diagnostics-v4.tsv"),
        _ => return None,
    };
    Some(
        source
            .lines()
            .filter(|line| !line.is_empty() && !line.starts_with('#'))
            .map(|line| {
                let fields = line.split('|').collect::<Vec<_>>();
                assert_eq!(fields.len(), 6, "invalid diagnostics catalog row: {line}");
                (fields[0], fields[2])
            })
            .collect(),
    )
}

fn escape_json(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\u{08}' => escaped.push_str("\\b"),
            '\u{0c}' => escaped.push_str("\\f"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            control if control <= '\u{1f}' => {
                escaped.push_str(&format!("\\u{:04x}", u32::from(control)));
            }
            other => escaped.push(other),
        }
    }
    escaped
}

#[cfg(test)]
mod tests {
    use std::collections::{BTreeMap, BTreeSet};
    use std::fs;
    use std::path::Path;

    use super::{Diagnostic, catalog_contract, diagnostic_matches_catalog_version};

    #[test]
    fn every_active_runtime_diagnostic_matches_active_catalog_v4() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"));
        let catalog = include_str!("../catalogs/diagnostics-v4.tsv")
            .lines()
            .filter(|line| !line.is_empty() && !line.starts_with('#'))
            .map(|line| {
                let fields = line.split('|').collect::<Vec<_>>();
                assert_eq!(fields.len(), 6, "invalid catalog row: {line}");
                (
                    fields[0].to_owned(),
                    (
                        fields[1].to_owned(),
                        fields[2].to_owned(),
                        fields[4].parse::<u8>().expect("catalog exit code"),
                    ),
                )
            })
            .collect::<BTreeMap<_, _>>();
        let mut active = BTreeSet::new();
        collect_codes(&root.join("src"), &mut active);
        assert!(!active.is_empty());
        for code in active {
            let static_code = catalog
                .get_key_value(&code)
                .map(|(code, _)| code.as_str())
                .unwrap_or_else(|| panic!("active code absent from diagnostics-v4.tsv: {code}"));
            // Los códigos del productor son literales estáticos; esta conversión
            // solo presta esa vida al test sin ampliar la API de producción.
            let static_code: &'static str = Box::leak(static_code.to_owned().into_boxed_str());
            let contract = catalog_contract(static_code)
                .unwrap_or_else(|| panic!("active code absent from runtime catalog: {code}"));
            let declared = &catalog[&code];
            assert_eq!(contract.category, declared.0, "category for {code}");
            assert_eq!(contract.meaning, declared.1, "meaning for {code}");
            assert_eq!(contract.exit_code, declared.2, "exit code for {code}");
            let diagnostic = Diagnostic::from_catalog(static_code);
            assert_eq!(diagnostic.code(), code);
            assert_eq!(diagnostic.category(), declared.0);
            assert_eq!(diagnostic.message(), declared.1);
            assert_eq!(diagnostic.exit_code(), declared.2);
            assert!(diagnostic_matches_catalog_version("4", &code, &declared.1));
            assert!(!diagnostic_matches_catalog_version(
                "4",
                &code,
                "significado alterado"
            ));
        }
        let unknown_code = format!("{}FMT-999", "GPS-");
        assert!(!diagnostic_matches_catalog_version(
            "4",
            &unknown_code,
            "significado inventado"
        ));
    }

    #[test]
    fn versioned_catalog_validation_does_not_backport_f5_diagnostics() {
        let message = "El escenario F5, su objetivo, F/P o sus límites no satisfacen el contrato de admisión.";
        assert!(!diagnostic_matches_catalog_version(
            "3",
            "GPS-CONF-101",
            message
        ));
        assert!(diagnostic_matches_catalog_version(
            "4",
            "GPS-CONF-101",
            message
        ));
    }

    fn collect_codes(directory: &Path, output: &mut BTreeSet<String>) {
        for entry in fs::read_dir(directory).expect("source directory") {
            let path = entry.expect("source entry").path();
            if path.is_dir() {
                collect_codes(&path, output);
            } else if path.extension().and_then(|value| value.to_str()) == Some("rs") {
                let source = fs::read_to_string(path).expect("UTF-8 Rust source");
                for token in source.split(|character: char| {
                    !(character.is_ascii_uppercase()
                        || character.is_ascii_digit()
                        || character == '-')
                }) {
                    if is_diagnostic_code(token) {
                        output.insert(token.to_owned());
                    }
                }
            }
        }
    }

    fn is_diagnostic_code(value: &str) -> bool {
        let parts = value.split('-').collect::<Vec<_>>();
        parts.len() == 3
            && parts[0] == "GPS"
            && !parts[1].is_empty()
            && parts[1].bytes().all(|byte| byte.is_ascii_uppercase())
            && parts[2].len() == 3
            && parts[2].bytes().all(|byte| byte.is_ascii_digit())
    }
}
