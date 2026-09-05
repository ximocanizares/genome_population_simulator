# Validación integral F5 — 5 de septiembre de 2026

Estado revisado: HEAD `cab105c`, con las modificaciones locales de F5.1.2.3
y su segunda corrección aún sin confirmar. Estas comprobaciones no constituyen
aprobación de F5.1.2.3. No se modificó código para obtener estos resultados.

## Resultados

| Comprobación | Resultado observado |
|---|---|
| `cargo test --lib f5_profile::tests::` | 11 superadas |
| `cargo test --test f5_bounded_input --test f5_bootstrap_dispatch --test f5_resource_preflight --test cli_config_diagnostics` | Código 101: CLI 7/7; bootstrap 7/9; Cargo detuvo los targets posteriores |
| `cargo test --all-targets --all-features --no-fail-fast` | Código 101: 663 superadas, 7 fallidas, 12 ignoradas; 92 resultados de targets |
| `cargo fmt --check` | Código 0 |
| `git diff --check` | Código 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | Código 101: `too_many_arguments` en `src/f5_profile.rs:606`, `parse_integer_intervals` (8 argumentos, máximo 7) |
| `cargo doc --no-deps` | Código 0 |
| `cargo test --test publication_process complete_failed_and_incomplete_are_real_minimal_trees -- --exact` | Código 101: reproduce el timeout de publicación por separado |

La ejecución completa usó `--no-fail-fast` para recoger los fallos de todos los
targets. Su salida de trabajo está en
`/tmp/gps-f5-review-20260905-tests.log`; ese archivo temporal no es evidencia
persistente ni está versionado. Los resultados relevantes se conservan aquí.

## Fallos de la ejecución completa

1. `f5_layout::tests::production_v5_instrumentation_covers_fragmentation_interior_last_residues_and_retry`:
   fallo H10 previamente reservado a F5.1.2.6; `operation_units=4957/4903`.
2. `dispatcher_opens_each_bootstrap_toml_once_and_records_identity_from_written_bytes`:
   `GPS-INPUT-001`, no existe el perfil del fixture al continuar desde bootstrap.
3. `fixed_width_bootstrap_memory_limit_is_inclusive_and_minus_one_stops_before_profile_and_rng`:
   `GPS-RES-101` al entrar en metadata del perfil, con acumulado y límite
   `131784`; informa que no queda presupuesto para abrir la entrada.
4. `cli_admits_or_rejects_profile_document_bytes_before_materializing_them`:
   comparación fallida, observado `0`, esperado `8184`.
5. `fragmented_mask_changes_only_tract_terms_and_is_bounded_inclusively`:
   comparación fallida, observado `71489`, esperado `73661`.
6. `retained_sequence_id_peak_is_persisted_in_the_published_memory_bound`:
   comparación fallida, observado `104257`, esperado `106429`.
7. `complete_failed_and_incomplete_are_real_minimal_trees`:
   `timeout waiting for provisionals`, reproducido también aisladamente.

La clasificación entre defectos del código, expectativas desactualizadas y
contradicciones del contrato corresponde a la revisión adversarial independiente.
No se presenta la suite como aprobada ni se atribuyen los siete fallos al H10.

## Estado posterior al primer ciclo A1–A7

Los resultados anteriores describen la revisión que descubrió los hallazgos y
se conservan como evidencia histórica. Tras aplicar A1–A7 en el árbol local se
obtuvieron estos resultados, que tampoco acreditaron por sí solos un cierre:

| Comprobación posterior | Resultado observado |
|---|---|
| `cargo test --lib f5_profile` | 13/13 |
| `cargo test --lib f5_admission::tests::` | 2/2 |
| Targets focales de bootstrap, entrada acotada, preflight y diagnósticos CLI | 57/57 |
| Targets de identidad, catálogo RNG, admisión y esquemas externos | 20/20 |
| Regresión aislada de publicación | 1/1 |
| `cargo check --all-targets --all-features` | Código 0 |
| `cargo fmt --all -- --check` | Código 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | Código 0 |
| `cargo doc --no-deps` | Código 0 |
| `cargo test --all-targets --all-features --no-fail-fast` | Código 101: un único fallo en `--lib`, H10; todos los targets posteriores terminaron verdes |

El fallo residual es
`f5_layout::tests::production_v5_instrumentation_covers_fragmentation_interior_last_residues_and_retry`.
Su fixture conserva un límite histórico de operaciones y el preflight corregido
informa `operation_units=5005/4951`, de los que 48 corresponden a la nueva cota
del soporte de gaps. H10 sigue reservado a F5.1.2.6; no se afirma que haya
pasado ni se amplió el alcance para corregirlo.

Una revisión independiente posterior encontró R1–R3: versión externa v6 para
un cálculo v7, falta de reserva prospectiva de la memoria TSV conocida y una
dimensión diagnóstica incorrecta. Esos hallazgos y los gates repetidos después
de aplicarlos están en
`evidence/f5.1.2/aplicacion-correcciones-integrales-2026-09-05.md`.
F5.1.2.3 permanece abierto y requiere otra revisión independiente de la pausa
A; esta evidencia histórica no afirma que R1–R3 hayan sido aprobados.
