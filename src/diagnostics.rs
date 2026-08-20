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

fn diagnostic_catalog() -> &'static BTreeMap<&'static str, CatalogContract> {
    static CATALOG: OnceLock<BTreeMap<&'static str, CatalogContract>> = OnceLock::new();
    CATALOG.get_or_init(|| {
        let source = include_str!("../catalogs/diagnostics-v2.tsv");
        assert!(
            source
                .lines()
                .any(|line| line == "# diagnostic_catalog_version=2"),
            "diagnostics-v2.tsv does not declare version 2"
        );
        let mut catalog = BTreeMap::new();
        for line in source
            .lines()
            .filter(|line| !line.is_empty() && !line.starts_with('#'))
        {
            let fields = line.split('|').collect::<Vec<_>>();
            assert_eq!(fields.len(), 6, "invalid diagnostics-v2.tsv row: {line}");
            let code = fields[0];
            let contract = CatalogContract {
                code,
                category: fields[1],
                meaning: fields[2],
                exit_code: fields[4]
                    .parse::<u8>()
                    .expect("invalid diagnostics-v2.tsv exit code"),
            };
            assert!(
                catalog.insert(code, contract).is_none(),
                "duplicate diagnostics-v2.tsv code: {code}"
            );
        }
        catalog
    })
}

pub(crate) fn diagnostic_matches_catalog(code: &str, message: &str) -> bool {
    diagnostic_catalog()
        .get(code)
        .is_some_and(|contract| message == contract.meaning)
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

    use super::{Diagnostic, catalog_contract, diagnostic_matches_catalog};

    #[test]
    fn every_active_runtime_diagnostic_matches_catalog_v2() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"));
        let catalog = include_str!("../catalogs/diagnostics-v2.tsv")
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
                .unwrap_or_else(|| panic!("active code absent from diagnostics-v2.tsv: {code}"));
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
            assert!(diagnostic_matches_catalog(&code, &declared.1));
            assert!(!diagnostic_matches_catalog(&code, "significado alterado"));
        }
        let unknown_code = format!("{}FMT-999", "GPS-");
        assert!(!diagnostic_matches_catalog(
            &unknown_code,
            "significado inventado"
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
