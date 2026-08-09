//! Diagnósticos públicos estables independientes de errores de crates.

use std::collections::BTreeMap;

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

    pub(crate) fn resource(code: &'static str, message: &'static str) -> Self {
        Self {
            code,
            category: "RES",
            message,
            exit_code: 3,
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
