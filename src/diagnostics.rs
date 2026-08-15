//! Diagnósticos públicos estables independientes de errores de crates.

use std::collections::BTreeMap;
use std::sync::OnceLock;

const DIAGNOSTIC_CATALOG_V2: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/catalogs/diagnostics-v2.tsv"
));

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct DiagnosticContract {
    pub category: &'static str,
    pub meaning: &'static str,
    pub required_context: Vec<&'static str>,
}

static DIAGNOSTIC_CONTRACTS: OnceLock<Result<BTreeMap<&'static str, DiagnosticContract>, String>> =
    OnceLock::new();

pub(crate) fn diagnostic_contracts()
-> Result<&'static BTreeMap<&'static str, DiagnosticContract>, String> {
    DIAGNOSTIC_CONTRACTS
        .get_or_init(parse_diagnostic_contracts)
        .as_ref()
        .map_err(Clone::clone)
}

fn parse_diagnostic_contracts() -> Result<BTreeMap<&'static str, DiagnosticContract>, String> {
    let mut contracts = BTreeMap::new();
    for (index, line) in DIAGNOSTIC_CATALOG_V2.lines().enumerate() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut fields = line.split('|');
        let code = fields
            .next()
            .ok_or_else(|| format!("diagnostics-v2 línea {} sin código", index + 1))?;
        let category = fields
            .next()
            .ok_or_else(|| format!("diagnostics-v2 línea {} sin categoría", index + 1))?;
        let meaning = fields
            .next()
            .ok_or_else(|| format!("diagnostics-v2 línea {} sin significado", index + 1))?;
        let required_context = fields
            .next()
            .ok_or_else(|| format!("diagnostics-v2 línea {} sin contexto", index + 1))?;
        let exit_code = fields
            .next()
            .ok_or_else(|| format!("diagnostics-v2 línea {} sin salida", index + 1))?;
        let planned_steps = fields
            .next()
            .ok_or_else(|| format!("diagnostics-v2 línea {} sin pasos", index + 1))?;
        if fields.next().is_some()
            || meaning.is_empty()
            || required_context.is_empty()
            || exit_code.is_empty()
            || planned_steps.is_empty()
            || code
                .strip_prefix("GPS-")
                .and_then(|value| value.split_once('-'))
                .map(|pair| pair.0)
                != Some(category)
        {
            return Err(format!("diagnostics-v2 línea {} incoherente", index + 1));
        }
        let contract = DiagnosticContract {
            category,
            meaning,
            required_context: required_context.split(',').collect(),
        };
        if contracts.insert(code, contract).is_some() {
            return Err(format!("diagnostics-v2 repite el código {code}"));
        }
    }
    if contracts.is_empty() {
        return Err("diagnostics-v2 no contiene contratos".to_owned());
    }
    Ok(contracts)
}

/// Diagnóstico estructurado que alimenta todos los renderizados públicos.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Diagnostic {
    code: &'static str,
    category: &'static str,
    message: &'static str,
    exit_code: u8,
    context: BTreeMap<String, String>,
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

    pub(crate) fn configuration(code: &'static str, message: &'static str) -> Self {
        Self {
            code,
            category: "CONF",
            message,
            exit_code: 2,
            context: BTreeMap::new(),
        }
    }

    pub(crate) fn input(code: &'static str, message: &'static str) -> Self {
        Self {
            code,
            category: "INPUT",
            message,
            exit_code: 2,
            context: BTreeMap::new(),
        }
    }

    pub(crate) fn format(code: &'static str, message: &'static str) -> Self {
        Self {
            code,
            category: "FMT",
            message,
            exit_code: 2,
            context: BTreeMap::new(),
        }
    }

    pub(crate) fn resource(code: &'static str, message: &'static str) -> Self {
        Self {
            code,
            category: "RES",
            message,
            exit_code: 3,
            context: BTreeMap::new(),
        }
    }

    pub(crate) fn domain(code: &'static str, message: &'static str) -> Self {
        Self {
            code,
            category: "DOMAIN",
            message,
            exit_code: 1,
            context: BTreeMap::new(),
        }
    }

    pub(crate) fn io(code: &'static str, message: &'static str) -> Self {
        Self {
            code,
            category: "IO",
            message,
            exit_code: 5,
            context: BTreeMap::new(),
        }
    }

    pub(crate) fn internal(code: &'static str, message: &'static str) -> Self {
        Self {
            code,
            category: "INTERNAL",
            message,
            exit_code: 70,
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
