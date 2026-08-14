# Revisión de F1.5

## Aplicación de la revisión independiente H-F15-R01–R05

| ID | Gravedad | Decisión inicial | Estado final | Prueba principal | Seguimiento |
|---|---|---|---|---|---|
| H-F15-R01 | Grave | Aplicar: admisión y ejecución discrepaban para composiciones parciales. | Aplicado | `f1_5_acceptance::{partial_ac_composition_is_completed_with_zero_g_t_before_cli_execution,generated_composition_rejects_nonexact_ids_categories_and_duplicates_at_admission}` | Ninguno. |
| H-F15-R02 | Grave | Aplicar; la segunda revisión confirmó que el mero enlace no cambiaba la representación ya parseada desde el original. | Aplicado | `reference_copy::authoritative_representation_is_loaded_from_the_copy_and_keeps_original_provenance` | Producción reparsea la copia catalogada con la misma función probada y conserva por separado la procedencia original. |
| H-F15-R03 | Grave | Aplicar: el constructor público sin artefacto no podía satisfacer el esquema. | Aplicado con adaptación | `results_publication::public_effective_document_builder_is_schema_valid_before_one_nested_mutation` | La base incompleta pasó a privada; el único constructor público exige `ReferenceOutput` autoritativo y retorna `Result`. |
| H-F15-R04 | Menor | Aplicar sin inventar el pico temporal ausente. | Aplicado | `scripts/run-f1-reference-benchmarks.sh A 3` y `B 3` | Una ejecución gemela detenida antes de `renameat2` permite medir el provisional con `du -sb`; salida idéntica, pico no nulo y residuo cero son puertas obligatorias. |
| H-F15-R05 | Menor | Aplicar: samtools reutilizaba el FAI interno. | Aplicado | `reference_index_acceptance::normal_f1_run_has_zero_indexes_and_acceptance_query_creates_fai_on_demand` | FASTA/FAI externos separados y verdad literal `CGTA`. |

Balance: cinco hallazgos contabilizados, cinco atendidos; cuatro `Aplicado` y uno
`Aplicado con adaptación`. Ningún hallazgo desapareció ni se marcó como
corregido sin regresión.

La instrumentación posterior produjo dos residuos menores, H-F15-R04-R1 y
H-F15-R04-R2. Ambos quedaron aplicados y conservan prueba ejecutable; no hay
residuos de esa corrección pendientes.

La pasada alternativa posterior confirmó explícitamente los cuatro bordes
pedidos: CLI A/C hasta publicación, FAI externos sin índice compartido, base
positiva antes de cada corrupción de `effective-config` y carga autoritativa
desde la copia después de mutar el original. La segunda revisión mantuvo R02
abierto porque la primera regresión mutaba el origen demasiado tarde; la nueva
regresión cubre el orden temporal exigido y la ruta de producción usa el mismo
adaptador `read_path`/`origin_path`. No se produjo un sexto hallazgo.

## Hallazgos recibidos y corregidos

- **H-F15-E1-01 — resuelto con decisión explícita:** los cuatro `build.rs`
  transitivos son admisibles tras auditar efectos; ninguno ejecuta procesos,
  usa red, modifica fuentes ni genera/compila código nativo.
- **H-F15-E2-01 — aplicado:** `ReferenceSequence` conserva la `DnaSequence`
  validada completa; la regresión comprueba las bases.
- **H-F15-E2-02 — aplicado:** la descripción conserva los bytes de frontera y
  su codificación reversible.

## Hallazgos de implementación E3–E12

- **H-F15-E3-01 — aplicado:** Clippy detectó un literal de slice innecesario en
  la prueba; se corrigió antes de cerrar E3.
- **H-F15-E8-01 — aplicado:** el primer esquema seguía describiendo solo la
  declaración F1.4 y rechazaba metadatos reales de referencia. Se amplió con
  contratos cerrados para origen, cabecera/descripcion reversible, FASTA,
  formato canónico, índices y enlaces de evento; las pruebas externas pasan.
- **H-F15-E13-03 — aplicado:** la puerta completa reveló que una prueba
  preparatoria F1.3 esperaba todavía `catalog_version = 1`. Se actualizó solo la
  expectativa al catálogo 2; las huellas y los cuatro vectores del catálogo 1
  siguen protegidos sin cambio.
- **H-F15-PERF-01 — medio, medido y aceptado como línea base:** el escritor
  propietario de noodles duplica temporalmente una secuencia generada; a 100 Mb
  el RSS es ~194 MiB. No existe umbral aprobado y cambiar la frontera sin perfil
  adicional sería especulativo. Queda documentado, no oculto.
- **H-F15-R04-R1 — resuelto durante la instrumentación:** el primer observador
  por sondeo podía perder por completo el provisional de A y devolver cero. La
  versión final usa `strace --seccomp-bpf`, retrasa 100 ms `renameat2`, exige un
  pico no nulo y compara los árboles ordinario e instrumentado. La ejecución
  temporizada no incluye la sonda.
- **H-F15-R04-R2 — menor, aplicado:** la ausencia de `strace` y una divergencia
  entre árboles podían terminar sin diagnóstico propio. El script identifica
  ambas causas antes de abandonar la medición.

No quedan hallazgos bloqueantes o graves de código abiertos en E3–E12.

## Hallazgos de reauditación E13

- **H-F15-E13-01 — resuelto por decisión funcional:** REQ-REF-05 limita v1 a
  cero índices normales o un único FAI bajo demanda, que F1 demuestra con un
  consumidor y un oráculo independiente. DTF-048 reserva los índices múltiples;
  no se inventa un segundo índice.
- **H-F15-E13-02 — resuelto por separación de obligaciones:** REQ-EVO-63 queda
  satisfecho en F1.5 con la infraestructura y los enlaces propios del evento.
  REQ-EVO-78 exige en F5.6 enlazar los resultados reales cuando existan, sin
  inventarlos en F1.

La revisión marca 27/27 casos de requisitos superados y E10/E11 cerrados. E13
queda completo y F1.5 alcanza su pausa obligatoria; esto no aprueba todavía la
pausa ni declara F1 cerrada.
