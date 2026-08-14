# Validación de F1.5

- Fecha de actualización: 14 de agosto de 2026.
- Plataforma: Ubuntu 24.04 sobre WSL2, x86-64, ext4.
- Estado: E1–E13 ejecutados y 27/27 casos de requisitos superados; pausa F1.5
  aprobada el 14 de agosto de 2026.

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
- E10/E11: release portable; A (100 kbp) y B (1/10/100 Mb), tres repeticiones
  cada tamaño. Se miden preflight y recorrido completo, pico temporal directo y
  residuo. La sonda gemela usa `strace` 6.8 con filtro seccomp y retraso de
  `renameat2`, exige salida idéntica y no contamina los tiempos ordinarios. C de
  1 Gb queda programado y no ejecutado conforme a DTC-046. La sonda se ejecutó
  fuera del sandbox restringido porque necesita `ptrace`; dentro de él falla de
  forma explícita y nunca registra cero como sustituto.
- E12: tar listado/extraído y verificado; OCI por digest construido,
  inspeccionado y arrancado. Digests en `distribution.md`.
- E13: 27 casos reaudidatos en `requirements.md`: 27 superados; puertas
  generales repetidas después de cerrar E10/E11.

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

## Cierre de la pausa

- TEST-REQ-REF-05 queda superado: v1 exige cero índices normales o un único FAI
  bajo demanda; DTF-048 reserva la coexistencia de índices múltiples.
- TEST-REQ-EVO-63 queda superado con directorio, manifiesto y enlaces propios de
  evento; REQ-EVO-78 planifica en F5.6 los enlaces a resultados reales.
- E10/E11 quedan cerrados con instrumentación y resultados completos.
- El commit `b5dc2207a36561fa447e9015d7be72e33f7235ac` superó el trabajo
  `Cargo gates` de GitHub Actions sobre Ubuntu 24.04 nativo en la ejecución
  [31783456974](https://github.com/ximocanizares/genome_population_simulator/actions/runs/31783456974):
  formato, Clippy, pruebas objetivo, doctests y compilación release portable.
- Con la evidencia local WSL2/ext4, OCI y CI nativo no queda ningún pendiente
  de la pausa F1.5. Su aprobación no cierra F1, que continúa con F1.6.
