# Cierre local F2 — verdad estructurada y VCF

## Estado y alcance

La implementación de F2 queda **cerrada localmente de forma condicional** para
la vía directa: una configuración general
produce la referencia, el estado haplotípico, `outputs/truth.jsonl` y, solo
cuando se solicita, una vista VCF 4.5 multimuestra. `truth.jsonl` permanece
obligatorio e independiente; VCF nunca reconstruye, corrige ni alimenta la
verdad.

F2.5.7 completa la frontera Marimo de DTC-051/DTC-079 sin añadir semántica
biológica, servicios ni infraestructura de evidencia. La pestaña `Salidas`
ofrece `none` (valor por defecto), `text` y `bgzf`. La omisión no serializa
`[vcf_output]`; `text` serializa un VCF sin índice y `bgzf` delega en Rust la
elección del único TBI o CSI. Tras una ejecución correcta, Marimo lee solo el
manifiesto, comprueba ruta/tamaño/catálogo de los enlaces VCF y muestra un
resumen conciso; no abre ni carga el VCF.

## Trazabilidad comprobada

- **DTC-051 / REQ-UI vigente:** `ui/config_editor.py` importa y exporta la
  tabla opcional cerrada `[vcf_output]`; `ui/gps_app.py` la expone con los tres
  únicos estados autorizados. Rust conserva admisión, valores efectivos,
  publicación y validación.
- **REQ-MAN-27, REQ-MAN-57 y REQ-SAL-21 (alcance F2):**
  `tests/vcf_publication.rs` verifica el enlace del VCF a referencia y muestras,
  el catálogo y estados completos, modos/rutas, SHA-256 y exclusividad del
  índice. El resumen UI verifica precisamente esos IDs y artefactos físicos sin
  interpretar bytes VCF.
- **REQ-VAL-28 y verdad independiente:** `tests/vcf_semantic_validation.rs`
  ejerce el oráculo independiente FASTA+`truth.jsonl`+VCF para REF/ALT/GT,
  fase, dosis, muestras y `GPS_CASE`, antes de abrir VCF cuando la verdad no es
  válida. Exige una vez `##reference=<reference_id efectivo>`,
  `##gps_execution_id=<execution_id derivado>`, `##gps_qual=.`,
  `##FILTER=<ID=PASS,...>` y `##FORMAT=<ID=GT,...>`; omisiones o intercambios
  se rechazan. `tests/truth_validation.rs` conserva la reconstrucción sin VCF.
- **REC-26 vigente:** `tests/vcf_semantic_validation.rs::fully_deleted_snp_is_absent_from_truth_and_vcf`
  recorre verdad actual→proyección→VCF y comprueba que un SNP cubierto por una
  deleción posterior permanece solo en el historial, sin genotipo SNP ni
  registro VCF SNP.
- **Recorrido F2.5.7:** `tests_ui/test_ui_contract.py` compara TOML manual y
  TOML de UI con `vcf_output = "text"`, valida ambos y comprueba que ambos
  producen `truth.jsonl`, `truth.vcf` y la misma configuración efectiva.
  `tests_ui/test_f2_5_vcf_ui.py` cubre la ruta UI
  configuración→`truth.jsonl`→BGZF+índice. `tests_ui/test_vcf_summary.py`
  cubre texto, BGZF+TBI, omisión y enlaces incompletos.

## Validación ejecutada

Todos estos comandos terminaron con código 0:

```text
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --all-features --locked -- -D warnings
cargo test --workspace --all-targets --all-features --locked
cargo test --workspace --doc --all-features --locked
cargo doc --workspace --no-deps --all-features --locked
uv run --frozen python -m unittest discover -s tests_ui
timeout 60s env MARIMO_SKIP_UPDATE_CHECK=1 uv run --frozen marimo check --strict ui/gps_app.py
git diff --check
```

La suite Cargo dejó dos pruebas informativas omitidas: el oráculo de
normalización externo DTC-066 y el auxiliar de subproceso. No son una
aceptación externa superada.

## DTC-066 externo pendiente: condición del cierre

La puerta `scripts/validate_vcf_external.sh` exige exactamente `bcftools 1.24`
y `tabix (htslib) 1.24`. Este entorno tiene `bcftools 1.21` y
`tabix (htslib) 1.19`; por tanto no se ejecutó ni se acredita la validación
externa. Es una condición pendiente, no una sustitución por los lectores,
índices y oráculo semántico locales. Tampoco se generó OCI, servicios remotos,
escenarios ni formatos de fases futuras.

## Pausa F2.5

Antes de iniciar F3 corresponde la revisión independiente de las invariantes
genómicas, la trazabilidad y la ruta local completa. Esa revisión y la puerta
externa DTC-066 1.24 son condiciones posteriores al cierre local; ninguna de
las dos queda declarada como superada aquí.
