# Validación de F1.5

- Fecha: 10 de agosto de 2026.
- Plataforma: Ubuntu 24.04 sobre WSL2, x86-64, ext4.
- Estado: E1–E12 ejecutados; E13 reauditado con 25/27 casos superados y dos
  parciales. F1.5 llega a la pausa sin cierre total.

## Pasos y pruebas específicas

- E1: `noodles-fasta = 0.66.0` aceptado; grafo, features, licencias, código
  nativo y cuatro `build.rs` auditados. `cargo audit --no-fetch`: código 0.
- E2: `cargo test --locked --test reference_domain`: 8 superadas.
- E3: `cargo test --locked --test fasta_validation`: 5 superadas; incluye
  vacío, IDs, descripción/no UTF-8, CRLF, gzip, gap, IUPAC y posición.
- E4/E8: `reference_copy` (4) y `results_publication` (19) superadas; copia
  exacta, mutación concurrente, carga autoritativa desde la copia, digest,
  enlaces y JSON Schema externo.
- E5/E6: `reference_generation` (4), `subseed_derivation` y
  `normative_catalogs` superadas. Catálogo v1 intacto; v2 añade
  `reference_generation`. FASTA canónico v1: 80 bases, LF, mayúsculas.
- E7: `reference_index_acceptance` superada; cero índices normales, FAI bajo
  demanda y consulta contrastada con `samtools faidx`.
- E9: `f1_5_acceptance` (4) superadas; CLI simulado repetido, CLI suministrado
  con N/cabecera/caja exactas y CLI inválido sin destino. Oráculos `cmp`,
  `sha256sum`, Python jsonschema y samtools terminaron correctamente.
- E10/E11, reabiertos: release portable; A (100 kbp) y B (1/10/100 Mb), tres repeticiones
  cada tamaño. Se miden preflight y recorrido completo, residuo temporal y se
  declara `not_measured` para el pico temporal no instrumentado. C de 1 Gb
  queda programado y no ejecutado.
- E12: tar listado/extraído y verificado; OCI por digest construido,
  inspeccionado y arrancado. Digests en `distribution.md`.
- E13: 27 casos reaudidatos en `requirements.md`: 25 superados y 2 parciales.

## Regresiones alternativas solicitadas

Ejecutadas individualmente antes de las puertas generales:

- `cargo test --locked --test f1_5_acceptance
  partial_ac_composition_is_completed_with_zero_g_t_before_cli_execution --
  --exact`: 1 superada. Recorre CLI `run`, admisión, normalización A/C→A/C/G/T,
  generación y publicación completa; el FASTA solo contiene A/C.
- `cargo test --locked --test reference_index_acceptance
  normal_f1_run_has_zero_indexes_and_acceptance_query_creates_fai_on_demand --
  --exact`: 1 superada. El adaptador y samtools usan copias y FAI separados y
  ambos devuelven la verdad literal `CGTA`.
- `cargo test --locked --test results_publication
  public_effective_document_builder_is_schema_valid_before_one_nested_mutation
  -- --exact`: 1 superada. Valida la base y rechaza una sola adición anidada.
- `cargo test --locked --test results_publication
  schemas_reject_unknown_fields_in_every_nested_f1_4_contract -- --exact`: 1
  superada. Cada ayudante negativo valida inmediatamente su documento base y
  después aplica exactamente una corrupción.
- `cargo test --locked --test reference_copy
  authoritative_representation_is_loaded_from_the_copy_and_keeps_original_provenance
  -- --exact`: 1 superada. Tras validar y copiar, muta el original antes de
  construir la representación autoritativa; cabecera y bases proceden de la
  copia, mientras `ReferenceOrigin` conserva la ruta original. Producción usa
  el mismo adaptador que separa `read_path` de `origin_path`.

## Puertas completas

- `cargo fmt --all -- --check`: código 0.
- `cargo clippy --workspace --all-targets --all-features --locked -- -D warnings`:
  código 0.
- `cargo test --workspace --all-targets --all-features --locked`: 138 pruebas
  superadas. La primera pasada detectó una expectativa F1.3 del catálogo 1; se
  corrigió para esperar el catálogo 2 sin cambiar sus vectores antiguos.
- `cargo test --workspace --doc --all-features --locked`: 1 `compile_fail`
  superada.
- `cargo doc --workspace --no-deps --all-features --locked`: código 0.
- `cargo build --release --workspace --all-features --locked`: código 0.
- `cargo audit --no-fetch`: 1199 avisos cargados, 65 dependencias examinadas,
  código 0; se conserva el aviso no bloqueante por lock del índice local.
- `git diff --check`: código 0.

## Pendientes de la pausa

- TEST-REQ-REF-05 no demuestra «varios» índices porque F1 solo tiene FAI y no
  existe un segundo tipo con consumidor real aprobado.
- TEST-REQ-EVO-63 enlaza referencia, configuración, subsemilla y manifiesto por
  evento, pero los resultados biológicos/técnicos reales siguen en F5.6.
- E10/E11 permanecen abiertos: el residuo temporal sí se mide, pero el pico de
  bytes temporales se declara `not_measured` y requiere instrumentación directa.
- La evidencia nativa Ubuntu 24.04 de esta implementación F1.5 deberá provenir
  de CI; WSL2 y OCI no se presentan como ejecución nativa independiente.
