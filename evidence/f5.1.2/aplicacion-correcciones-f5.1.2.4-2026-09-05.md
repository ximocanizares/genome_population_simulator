# Aplicación de correcciones F5.1.2.4 — 5 de septiembre de 2026

Estado: los hallazgos C1–C6, L1–L3 y C8–C9 están corregidos localmente y
cuentan con regresiones focales. La suite completa conserva un único fallo,
H10, que C7 identifica como trabajo ya reservado a F5.1.2.6. Este informe no
constituye la revisión independiente ni aprueba la pausa F5.1.2-B; F5.1.2.4
permanece sin marcar.

## Contabilidad de hallazgos

| ID | Gravedad | Decisión inicial | Estado final | Archivos | Pruebas | Seguimiento |
|---|---|---|---|---|---|---|
| C1 | Grave | Aplicar | Aplicado con adaptación | `src/application.rs`, `tests/f5_resource_preflight.rs` | La configuración provisional carga perfil y conserva ledger sin semilla; referencia, máscara y fórmula se admiten antes de resolver la entropía. El espía comprueba cero llamadas ante FASTA/máscara ausente, inválido, con SHA incorrecto o presupuesto menos uno, y una llamada tras el éxito. | La semilla declarada conserva su comportamiento. El estado todavía no es un receipt sellado. |
| C2 | Grave | Aplicar | Aplicado con adaptación | `src/f5_admission.rs`, `src/f5_preflight.rs`, `tests/f5_resource_preflight.rs` | Cada `Read::read` reserva primero el bloque prospectivo y liquida después los bytes observados, devolviendo la capacidad no usada incluso en EOF. La máscara usa una lectura acotada sin read-ahead y reserva su línea antes de abrir el handle. Los límites exacto/menos uno y el lector espía prueban que no se pide un byte no admitido. | La segunda visita FASTA se contabiliza, pero su lectura física sigue reservada a .5. |
| C3 | Grave | Aplicar | Aplicado con adaptación | `src/application.rs`, `src/f5_preflight.rs`, `tests/f5_resource_preflight.rs` | El ledger de bootstrap/perfil continúa en el escaneo y en `FinalFormula`. `other` queda ligado a esa admisión inicial antes de referencia/máscara; la ruta pública solo recupera el mismo resultado con el mismo valor y rechaza cualquier sustitución tardía sin E/S ni entropía. | El ledger y resultado finales quedan retenidos internamente para que .5 pueda sellarlos; no se recalculan desde cero. |
| C4 | Grave | Aplicar | Aplicado | `src/f5_admission.rs`, `src/f5_preflight.rs`, `tests/f5_resource_preflight.rs` | El mismo handle que entrega bytes al escáner calcula tamaño y SHA; el resultado no sellado conserva las identidades físicas de referencia y máscara en orden de rol. La prueba contrasta rol, ruta absoluta normalizada, tamaño y SHA, y ningún resultado se publica tras rechazo. | Proyección pública y sellado quedan expresamente para .5. |
| C5 | Grave | Aplicar | Aplicado con adaptación | `src/f5_preflight.rs`, `tests/f5_resource_preflight.rs` | El oráculo independiente suma retención viva de perfil y coexistencia de buffers FASTA/máscara. IDs, tramos y líneas se cargan antes del crecimiento; el trabajo aprendido de IDs/secuencias/tramos se carga incrementalmente. El único backing `Vec` de tramos se comparte en las recuperaciones públicas, sin una clonación proporcional a la fragmentación. | No se añadió una segunda pasada; se modeló la coexistencia real de .4. |
| C6 | Grave | Aplicar | Aplicado | `src/application.rs`, `src/main.rs`, `tests/cli_config_diagnostics.rs`, `tests/f5_bootstrap_dispatch.rs` | `gps validate` ejecuta la admisión física .4. Las regresiones cubren caso declarado válido, FASTA ausente o malformado y máscara inválida, sin crear destino ni resolver entropía. | La CLI continúa sin ejecutar layout ni publicar artefactos F5. |
| C7 | Grave, solo prueba fuera de .4 | Aplazar | Aplazado | `src/f5_layout.rs` (sin cambios) | La suite completa reproduce únicamente H10: la reserva prospectiva alcanza `4961` con límite `4952`, delta +9. C8 retiró una unidad EOF no contractual del límite anterior, sin tocar el recorrido H10. | Destino: F5.1.2.6. No se alteró el fixture ni la fórmula para ocultar el fallo. |

### Corrección residual tras revisión independiente

| ID | Gravedad | Decisión inicial | Estado final | Evidencia focal | Seguimiento |
|---|---|---|---|---|---|
| L1 | Grave | Aplicar | Aplicado con adaptación | La API conserva la firma de `preflight_f5_resources`, pero ya no crea un ledger nuevo: recupera el resultado admitido inicialmente. Mismo `other` funciona tras retirar FASTA/máscara; otro distinto devuelve `GPS-CONF-101`, sin nueva E/S ni RNG. | Revisión independiente pendiente; .5 sellará el objeto ya ligado. |
| L2 | Grave | Aplicar | Aplicado con adaptación | `covered_bases` deriva y valida `B=ceil(T/Lmax)` antes del FASTA; `S_bound`, `H*S` y los términos triangulares conocidos también preceden la apertura. La frontera conocida exacta alcanza el sentinel FASTA ausente y exacta−1 rechaza por operaciones. Los costes de ID/secuencia/tramo se cargan antes de retener y las sondas EOF conservadoras mantienen inclusivo el presupuesto exacto. | `coverage_fraction` sigue resolviendo B después de conocer E, como exige su contrato. |
| L3 | Grave | Aplicar | Aplicado | El escáner transfiere el `Vec` al preflight sin `scanned.clone()` y las recuperaciones comparten su backing mediante `Arc`; 4096 tramos conservan el mismo puntero en dos resultados públicos. | El sellado de .5 podrá consumir la misma geometría compartida. |

### Corrección residual final

| ID | Gravedad | Decisión inicial | Estado final | Archivos | Pruebas | Seguimiento |
|---|---|---|---|---|---|---|
| C8 | Grave | Aplicar | Aplicado con adaptación | `src/f5_admission.rs`, `src/f5_preflight.rs`, `tests/f5_resource_preflight.rs` | La fórmula publicada elimina `1 + mask_present`. El ledger respalda cada sonda EOF con una parte del trabajo futuro ya cargado: cero bytes no cambia el total y un byte real consume respaldo y exige margen adicional en la fórmula final. Un oráculo desplegado desde los bytes literales de general, escenario, metadata, cinco TSV, FASTA y máscara, sin llamar a la fórmula productiva, coincide término por término con y sin máscara; el límite exacto admite y exacto−1 rechaza. El lector espía demuestra que cada byte entregado usa headroom o trabajo ya cargado. | Nueva revisión independiente de cierre; `.4` permanece abierto. |
| C9 | Grave | Aplicar | Aplicado | `src/f5_preflight.rs`, `tests/f5_resource_preflight.rs` | Dos máscaras sintácticamente válidas, una más corta y otra más larga que el contenido original cuyo SHA sigue declarado, producen `GPS-INPUT-002` con `expected_digest`/`observed_digest`, cero entropía y cero destino. Una cabecera o fila malformada sigue produciendo `GPS-FMT-101`. | Nueva revisión independiente de cierre; `.4` permanece abierto. |

La observación menor sobre dos accesores provisionales `#[doc(hidden)] pub` no
se aplicó: esta corrección no amplía ni cambia la interfaz pública. Su
visibilidad puede reevaluarse cuando .5 sustituya la proyección provisional por
el sellado interno.

## Explicación y límites

La corrección introduce una configuración F5 provisional interna: el
bootstrap y el perfil dejan un único ledger en `ReferenceAndMaskScan`, pero no
resuelven la semilla. El escaneo continúa ese ledger, cobra antes de cada
lectura o retención, calcula identidades desde los handles observados y alcanza
`FinalFormula`. Solo entonces la ruta de biblioteca resuelve una semilla
omitida. La ruta `validate` realiza la misma admisión, pero no necesita semilla
porque no existe consumidor RNG.

La reserva de lectura es transaccional. Se reserva el máximo bloque que el
presupuesto restante permite antes de llamar a `Read::read`; después se
conserva únicamente el coste de los bytes realmente entregados. La máscara no
usa un lector con read-ahead: carga prospectivamente su buffer acotado antes de
abrir y consume bytes mediante la misma frontera. La identidad física se
termina únicamente al alcanzar EOF válido y comprobar el SHA declarado.

No se creó ni publicó `F5AdmissionReceipt`, no se añadió la segunda pasada
FASTA y no se modificó H10. Tampoco hubo red, commits, publicación externa ni
agentes hijos. No hubo preguntas funcionales abiertas ni contradicciones que
obligaran a ampliar el alcance.

## Comandos y resultados reales

| Comprobación | Resultado observado |
|---|---|
| `cargo test --test f5_resource_preflight --no-fail-fast` | 43 superadas, 0 fallidas |
| `cargo test --test f5_bootstrap_dispatch --no-fail-fast` | 9 superadas, 0 fallidas |
| `cargo test --lib f5_admission --no-fail-fast` | 4 superadas, 0 fallidas |
| `cargo test --test cli_config_diagnostics --no-fail-fast` | 7 superadas, 0 fallidas |
| `cargo fmt --all -- --check` | Código 0 después de aplicar formato |
| `cargo clippy --all-targets --all-features -- -D warnings` | Código 0 |
| `cargo doc --no-deps` | Código 0 |
| `git diff --check` | Código 0 antes de redactar esta evidencia |
| `cargo test --all-targets --all-features --no-fail-fast` | Código 101: 671 superadas, 1 fallida y 12 ignoradas; único fallo H10 en `production_v5_instrumentation_covers_fragmentation_interior_last_residues_and_retry` (`4961/4952`) |
| `cargo test --lib f5_layout::tests::production_v5_instrumentation_covers_fragmentation_interior_last_residues_and_retry -- --exact` | Código 101: H10 reproducido aisladamente con `4961/4952`; al retirar la unidad EOF no contractual de la magnitud publicada, el déficit observado pasa de +8 a +9 sin modificar F5.1.2.6. |

No correspondía guardar una línea base de rendimiento: el trabajo corrige
orden de admisión, cotas prospectivas y contabilidad, y no modifica un
recorrido con gate de rendimiento establecido. No se reconstruyó una línea
base artificial.

## Pendientes graves

No queda un hallazgo C1–C6, L1–L3 o C8–C9 conocido sin corregir dentro de
F5.1.2.4. La conformidad final sigue pendiente de revisión independiente.
C7/H10 permanece deliberadamente pendiente con destino F5.1.2.6 y por ello la
suite completa no se declara verde.
