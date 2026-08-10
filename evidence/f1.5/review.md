# Revisión de F1.5

## Aplicación de la revisión independiente H-F15-R01–R05

| ID | Gravedad | Decisión inicial | Estado final | Prueba principal | Seguimiento |
|---|---|---|---|---|---|
| H-F15-R01 | Grave | Aplicar: admisión y ejecución discrepaban para composiciones parciales. | Aplicado | `f1_5_acceptance::{partial_ac_composition_is_completed_with_zero_g_t_before_cli_execution,generated_composition_rejects_nonexact_ids_categories_and_duplicates_at_admission}` | Ninguno. |
| H-F15-R02 | Grave | Aplicar; la segunda revisión confirmó que el mero enlace no cambiaba la representación ya parseada desde el original. | Aplicado | `reference_copy::authoritative_representation_is_loaded_from_the_copy_and_keeps_original_provenance` | Producción reparsea la copia catalogada con la misma función probada y conserva por separado la procedencia original. |
| H-F15-R03 | Grave | Aplicar: el constructor público sin artefacto no podía satisfacer el esquema. | Aplicado con adaptación | `results_publication::public_effective_document_builder_is_schema_valid_before_one_nested_mutation` | La base incompleta pasó a privada; el único constructor público exige `ReferenceOutput` autoritativo y retorna `Result`. |
| H-F15-R04 | Menor | Aplicar sin inventar el pico temporal ausente. | Aplicado con adaptación | `scripts/run-f1-reference-benchmarks.sh A 3` y `B 3` | Preflight/total y residuo quedan medidos; pico `not_measured`; E10/E11 reabiertos hasta instrumentarlo. |
| H-F15-R05 | Menor | Aplicar: samtools reutilizaba el FAI interno. | Aplicado | `reference_index_acceptance::normal_f1_run_has_zero_indexes_and_acceptance_query_creates_fai_on_demand` | FASTA/FAI externos separados y verdad literal `CGTA`. |

Balance: cinco hallazgos contabilizados, cinco atendidos; tres `Aplicado` y dos
`Aplicado con adaptación`. Ningún hallazgo desapareció ni se marcó como
corregido sin regresión.

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

No quedan hallazgos bloqueantes o graves de código abiertos en E3–E12.

## Hallazgos de reauditación E13

- **H-F15-E13-01 — bloquea el cierre total:** REQ-REF-05 exige poder representar
  cero, uno o varios índices. F1 demuestra cero y un FAI real, pero no tiene un
  segundo tipo ni consumidor autorizado; no se inventa uno.
- **H-F15-E13-02 — cierre compartido conservado:** REQ-EVO-63 sigue asignada a
  F1.5/F5.6. F1.5 completa los enlaces de referencia, pero no existen todavía
  resultados biológicos/técnicos reales de evento; permanece parcial.

La revisión marca 25/27 casos superados. E13 se ejecutó como auditoría, pero no
se aprueba ni se declara F1 cerrado mientras esos dos casos no tengan destino o
criterio explícitamente aprobado.
