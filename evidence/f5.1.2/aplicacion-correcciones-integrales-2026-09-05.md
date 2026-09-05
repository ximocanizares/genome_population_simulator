# Aplicación de correcciones F5.1.2-A — 5 de septiembre de 2026

Estado: A1–A7 se aplicaron en un primer ciclo, pero la revisión independiente
posterior encontró dos lagunas graves (R1/R2) y una menor (R3); por tanto, aquel
ciclo no constituyó una corrección integral ni aprobó la pausa. Una revisión de
cierre posterior volvió a rechazar R1 por una clave pública v7 heredada de v6.
La corrección recurrente está aplicada localmente con regresión rojo→verde y
queda pendiente de otra revisión independiente junto con R2/R3. F5.1.2.3
permanece abierto. H10 sigue reservado a F5.1.2.6 y no se presenta como
resuelto.

## Contabilidad de hallazgos

| ID | Decisión inicial | Estado aplicado | Corrección y evidencia de regresión |
|---|---|---|---|
| A1, bloqueante | Aplicar | Corregido | `F5BlockLengthPlan::from_preflight` conserva el `scenario_id` configurado. Una prueba por la API pública comprueba el plan y, con oráculos independientes, el sub-seed y el identificador de bloque; renombrar el escenario cambia ambos sin introducir el literal `scenario_id`. |
| A2, grave | Aplicar con adaptación contractual | Corregido | El preflight calcula con aritmética comprobada la cardinalidad exacta e inclusiva del soporte discontinuo de gaps y carga la cota conservadora `2 × L × B × S_bound × G`. La cota cubre los dos recorridos literales del selector y el gap terminal sin consumir RNG; no cambia la distribución ni la secuencia RNG. Hay pruebas de soporte discontinuo, exacto/menos uno y rango hasta `10^12`. La derivación queda registrada en DTC-089 y arquitectura §15. |
| A3, grave | Aplicar | Corregido | El TSV usa un arreglo fijo de campos y rechaza la primera columna sobrante, sin `Vec` proporcional al número de columnas. Las retenciones de auditoría se cargan antes de normalizar o conservar. Las regresiones cubren 200 000 campos bajo límite reducido y rechazo previo a normalización. |
| A4, grave | Aplicar con reconciliación | Corregido | El presupuesto de entrada resta memoria actualmente retenida, no el máximo histórico. El ledger separa retención viva y pico duro; bootstrap libera su buffer al cambiar de fase. La transición común de aplicación evita dos construcciones divergentes. Las pruebas cubren casos dominados por bootstrap y por TSV, exacto/menos uno y retención viva no nula, con oráculos independientes. |
| A5, grave de validación | Aplicar | Corregido | Se reabrió F5.1.2.3 y se eliminó la afirmación de cierre verde. El fixture de publicación usa la barrera existente y un TOML pequeño para producir el fallo controlado, sin depender de un documento superior a 65 536 bytes. La prueba aislada vuelve a ser verde. |
| A6, grave | Aplicar | Corregido | El layout compara el gap con la distancia restante antes de sumarlo. Un gap singleton `u64::MAX` termina por agotamiento normal `GPS-DOMAIN-106`, no por overflow interno. |
| A7, menor | Aplicar | Corregido | Un fallo inyectado de reserva deja los totales atómicamente intactos y el ledger terminalmente rechazado; una transición posterior no puede leer entrada. Los diagnósticos de recursos emplean el límite de la misma dimensión —operaciones, memoria o salida— que se agotó. |

Además, `parse_integer_intervals` agrupa el estado común en un contexto de
parseo pequeño; así desaparece `clippy::too_many_arguments` sin permitir el
lint.

## Derivación conservadora de A2

Sea `G` la suma comprobada de las anchuras inclusivas de todos los intervalos
disjuntos de gaps. En cada intento de un bloque ya materializado pueden caber a
lo sumo `S_bound - 1` gaps y existir una selección terminal que ya no cabe, por
lo que se acotan `S_bound` selecciones. El selector literal recorre una vez todo
el soporte para sumar pesos y puede recorrerlo otra vez para localizar el valor
seleccionado. La cota de recorridos obligatorios es por tanto
`2 × L × B × S_bound × G`, incluso cuando una ejecución concreta termine antes.
Esta sobreaproximación es deliberada, se calcula antes de crear o consumir el
RNG y conserva exactamente el algoritmo publicado.

## Hallazgos de la revisión independiente posterior

| ID | Gravedad | Decisión inicial | Estado final | Archivos | Pruebas y seguimiento |
|---|---|---|---|---|---|
| R1 | Grave | Aplicar | Aplicado con adaptación | `src/f5_preflight.rs`, `schemas/effective-config-v5.schema.json`, `tests/f5_0_1_external_schemas.rs`, `tests/f5_resource_preflight.rs` | La primera aplicación fue incompleta: productor y esquema publicaban `input_parser_peak_bytes`, ausente de la lista contractual v7. La corrección recurrente conserva ese valor solo para contabilidad interna y lo elimina de la proyección, de `required` y de `properties`; el manifiesto hereda la corrección mediante su `$ref`. Dos oráculos literales independientes del productor/esquema comprueban el conjunto escalar publicable por .3; efectivo y manifiesto aceptan v7 y rechazan v5, v6, cada ausencia y cualquier extra. `input_identities` y `phases` siguen explícitamente reservados a .5. Pendiente: revisión independiente de cierre; R1 no se declara cerrado. |
| R2 | Grave | Aplicar | Aplicado con adaptación | `src/f5_admission.rs`, `src/f5_profile.rs` | Tras metadata se calcula con aritmética comprobada el máximo `2*size_bytes` de los cinco TSV y se reserva antes de abrir el primero. La reserva eleva el pico admitido sin convertirlo en retención viva; las filas siguen cargándose incrementalmente. Una tabla posterior sobredimensionada y la primera ausente prueban rechazo `GPS-RES-101` en exacto−1 y apertura sólo con el límite exacto. Las regresiones previas de A4 siguen verdes. Pendiente: revisión independiente de cierre. |
| R3 | Menor | Aplicar | Aplicado | `src/f5_admission.rs` | La dimensión se elige por las capacidades ya divididas; operaciones gana el empate. Se prueban operaciones `3/2`, el caso simétrico de memoria y el empate, todos terminales. Pendiente: revisión independiente de cierre. |

## Corrección recurrente de R1

La comparación directa con DTC-089 (líneas 1400–1408 y tabla de la línea 1622),
arquitectura §15 y el reparto del plan confirmó que
`input_parser_peak_bytes` es contabilidad interna y no una clave pública v7.
Antes de editar se clasificó el único hallazgo de esta revisión como **Aplicar**,
grave y de confianza alta. La prueba exacta
`public_v7_preflight_serializes_gap_terms_and_validates_both_contracts` falló en
rojo: el conjunto observado contenía esa clave adicional frente al oráculo
literal.

La corrección mínima retiró la clave de `to_contract_json`, de `required` y de
`properties`. No retiró el campo interno, su accesor ni las pruebas de la cuenta
de memoria. El manifiesto no duplica la forma: consume la misma definición
corregida mediante
`effective-config-v5.schema.json#/$defs/preflight`. El oráculo del productor y
el oráculo separado de documentos completos enumeran literalmente las 27 claves
escalares publicables en .3; el segundo las contrasta por separado en efectivo
y manifiesto y no deriva la expectativa de `preflight()`, el productor o el
esquema. La forma final con `input_identities` y `phases` sigue siendo propiedad
de F5.1.2.5 y no se adelantó.

| Comprobación posterior a la corrección recurrente | Resultado observado |
|---|---|
| Regresión exacta del productor y ambas referencias de esquema | 1 superada; antes fallaba por `input_parser_peak_bytes` adicional |
| Casos externos exactos de efectivo y manifiesto | 2 superadas; v7 aceptada, v5/v6, cada clave ausente y clave adicional rechazadas en ambas rutas |
| `cargo test --test f5_resource_preflight --test f5_0_1_external_schemas --test cli_config_diagnostics --test f5_bootstrap_dispatch --test f5_bounded_input --test f5_0_1_config_admission --test f5_identity --test f5_rng_catalog --no-fail-fast` | 78 superadas, 0 fallidas |
| `cargo fmt --all -- --check` | Código 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | Código 0 |
| `cargo doc --no-deps` | Código 0 |
| `git diff --check` | Código 0 antes y después de redactar esta evidencia |
| `cargo test --all-targets --all-features --no-fail-fast` | Código 101: 666 superadas, 1 fallida y 12 ignoradas; único fallo H10 en `production_v5_instrumentation_covers_fragmentation_interior_last_residues_and_retry` (`5005/4951`) |

No correspondía una línea base de rendimiento: el cambio afecta únicamente a
una clave serializada, su esquema y los oráculos de contrato. No quedan
preguntas funcionales abiertas para esta corrección. Su estado es **Aplicado
con adaptación**, pendiente de revisión independiente; no equivale al cierre de
R1 ni de F5.1.2.3.

## Validación del primer ciclo A1–A7 (histórica)

| Comprobación | Resultado |
|---|---|
| `cargo test --lib f5_profile` | 13 superadas |
| `cargo test --lib f5_admission::tests::` | 2 superadas |
| `cargo test --test f5_bootstrap_dispatch --test f5_bounded_input --test f5_resource_preflight --test cli_config_diagnostics --no-fail-fast` | 57 superadas, 0 fallidas |
| `cargo test --test f5_identity --test f5_rng_catalog --test f5_0_1_config_admission --test f5_0_1_external_schemas --no-fail-fast` | 20 superadas, 0 fallidas |
| `cargo test --test publication_process complete_failed_and_incomplete_are_real_minimal_trees -- --exact` | 1 superada |
| `cargo check --all-targets --all-features` | Código 0 |
| `cargo fmt --all -- --check` | Código 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | Código 0 |
| `cargo doc --no-deps` | Código 0 |
| `git diff --check` | Código 0 antes de redactar esta evidencia |
| `cargo test --all-targets --all-features --no-fail-fast` | Código 101; único fallo: H10 en `production_v5_instrumentation_covers_fragmentation_interior_last_residues_and_retry`; todos los demás targets terminaron verdes |

La prueba H10 falla ahora en el preflight de su fixture por el límite histórico
de operaciones (`5005/4951`, incluidos 48 recorridos de soporte de gaps). No se
adaptó porque H10 está expresamente reservado a F5.1.2.6. Tampoco se ejecutaron
acciones de red, commits ni publicación remota.

## Estado de la pausa

Tras R1–R3 se repitieron las comprobaciones sobre el árbol corregido:

| Comprobación | Resultado |
|---|---|
| Regresiones rojo→verde de R1, R2 y R3 | Las tres fallaron antes de la corrección por el mecanismo descrito y pasan después |
| `cargo test f5_profile::tests` | 14 superadas |
| `cargo test f5_admission::tests` | 3 superadas |
| `cargo test --test f5_resource_preflight` | 34 superadas |
| `cargo test --test f5_0_1_external_schemas` | 6 superadas |
| Bootstrap, lector acotado, CLI, admisión, identidad y RNG F5 | 38 superadas |
| Regresión aislada de publicación | 1 superada |
| `cargo fmt --all -- --check` | Código 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | Código 0 |
| `cargo doc --no-deps` | Código 0 |
| `git diff --check` | Código 0 antes de actualizar esta evidencia |
| `cargo test --all-targets --all-features --no-fail-fast` | Código 101: 666 superadas, 1 fallida y 12 ignoradas; único fallo H10 en `production_v5_instrumentation_covers_fragmentation_interior_last_residues_and_retry` (`5005/4951`) |

F5.1.2.3 vuelve a ser candidato a revisión independiente de la pausa A, pero
no está aprobado. Los estados R1–R3 describen aplicación local, no conformidad
independiente; el checkbox permanece abierto.
