# Revisión adversarial integral de F5 — 5 de septiembre de 2026

## Veredicto

**No aprobar el cierre de F5 ni la pausa F5.1.2-A en este estado.** La revisión
encuentra un defecto silencioso en la identidad del escenario que alimenta RNG
e IDs, una contradicción entre el coste del muestreo de gaps y su presupuesto,
y fallos actuales en la admisión de memoria. No son solo expectativas antiguas
de pruebas ni ausencia de pasos futuros.

No se ha reparado código. Se ha utilizado `revisar-codigo-rust`, con todas sus
referencias, y `redactar-texto-tecnico` para separar defectos demostrados,
trabajo pendiente y propuestas. Esta revisión independiente se realizó con
gpt-6-astra, esfuerzo high, solicitado expresamente.

## 1. Alcance y autoridad

Base: `cab105c222795d1724d33d12a6675a052a91c4a4`, más el árbol no confirmado de
F5.1.2.3 disponible durante la revisión. Se inspeccionaron los cuatro módulos
`src/f5_*.rs`, sus entradas públicas en `application.rs` y `main.rs`,
validación TOML en `adapters.rs`, tipos de configuración, fronteras RNG/ID,
contratos y pruebas F5, y efectos sobre la suite general. El historial relevante
es `307995f`, `c73a14e`, `08376f3` y `cab105c`.

Autoridades: [DTC-089](../../design/decisiones-tecnicas-0.1.md),
[arquitectura §15](../../design/arquitectura-0.1.md) y
[plan F5](../../design/plan-implementacion-f5-0.1.md), con sus requisitos y
criterios de aceptación trazados; se contrastaron las
[evidencias previas](redisenio-admision-v7.md).
La autoridad exige identidad y RNG reproducibles, cotas prospectivas y una
admisión por fases; no autoriza cambiar la distribución biológica para abaratar
el algoritmo.

F5.1.2.4, .5 y .6 están pendientes. Por tanto, no se presentan como defectos
nuevos de .3 la ausencia del receipt sellado, la integración v7 de FASTA/máscara,
la segunda lectura física acotada, la sustitución final de esquemas v6 ni H10.
Sí se revisa el layout ya implementado: sus defectos científicos no quedan
excluidos porque el cableado v7 futuro esté pendiente. F5.1.3 y F5.2–F5.6 no se
evalúan como funcionalidades terminadas.

Se preservaron los cambios ajenos en `TO_DO.md`, plan VE, `.claude/` y
`benchmarks/`. No hubo red, agentes hijos, commits ni publicación externa.
Las reproducciones escribieron únicamente fixtures en `/tmp`.

## 2. Rutas públicas y alcance real de las pruebas

| Ruta | Qué ejecuta hoy | Límite de la evidencia existente |
|---|---|---|
| CLI → `load_and_validate_cli_config` | general y escenario acotados; valida perfil con un ledger nuevo; descarta perfil y devuelve bootstrap | El observador público instrumenta solo aperturas TOML, no los handles del perfil; las expectativas de “bootstrap sin perfil” quedaron obsoletas al añadir .3 |
| `load_and_validate_f5_config[_with_entropy]` | bootstrap, otro montaje de ledger, perfil, configuración efectiva y resolución de semilla | No conserva el ledger para el preflight; la unificación completa pertenece a .4/.5 |
| `preflight_f5_resources` | escáner heredado y fórmulas finales | Sigue anunciando v6 aunque el diff ya cambia operaciones y memoria; no equivale a admisión v7 sellada |
| `F5BlockLengthPlan::from_preflight` → `F5LayoutRequest` → layout | camino público productivo de geometría, RNG, IDs y recuperación REF | Los vectores unitarios RNG/ID no verifican qué escenario les entrega este camino |
| Pruebas externas de JSON Schema | validan objetos JSON construidos por el fixture | No reconstruyen un paquete F5 publicado ni las fórmulas; por ejemplo el fixture `preflight()` admite escalares artificiales iguales a 1 |
| Pruebas privadas de layout | mayoría construye `F5LayoutInputs` / `from_test_parts` | Son útiles para geometría, pero omiten el adaptador productivo que contiene A1; H10 sí usa producción, pero falla antes de alcanzar su instrumentación |

## 3. Compilación y pruebas

El agente raíz ejecutó los gates sobre el mismo árbol. Leí la salida completa
relevante en `/tmp/gps-f5-review-20260905-tests.log`; no repetí la suite.
El registro detallado queda en
[validación integral](validacion-integral-2026-09-05.md).

- `cargo fmt --check`, `git diff --check` y `cargo doc --no-deps`: correctos.
- Clippy estricto: falla en `src/f5_profile.rs:606`,
  `parse_integer_intervals`, 8 argumentos frente a 7.
- Suite `--all-targets --all-features --no-fail-fast`: **663 correctas,
  7 fallidas, 12 ignoradas**, 92 targets; fallan cuatro targets.
- Perfil: 11/11; CLI diagnostics: 7/7. Eso no acredita cierre integral.
- El fallo de publicación se repitió aislado por el agente raíz.

Clasificación de los siete fallos:

| Fallo | Interpretación |
|---|---|
| H10 productivo v5 | Preexistente y reservado a .6; ahora observa operaciones 4957/4903. Sigue impidiendo alcanzar el oráculo productivo |
| Bootstrap “una apertura e identidad” | Expectativa vieja: deja perfil ausente y llama a un dispatcher que ahora sí lo valida; no demuestra una reapertura TOML |
| Bootstrap memoria exacta | Además de cambiar el alcance de la prueba, revela A4: el límite 131784 ya admitido se convierte en headroom cero |
| Cita de perfil: delta de memoria 8184 | Expectativa v6 de sumar picos que no coexisten; con máximo de fases el delta puede ser 0. No corregir el código para recuperar esa suma |
| Fragmentación: memoria 73661 | Literal antiguo; resultado actual 71489. Debe reconstruirse la fase dominante, no copiar el nuevo literal |
| IDs retenidos: memoria 106429 | Literal antiguo; resultado actual 104257. Mismo problema de oráculo de memoria |
| Publicación: timeout de provisional | `tests/publication_process.rs:345` añade 750000 líneas al general v2; la cota universal de 65536 lo rechaza antes de crear el provisional. Interacción con el nuevo dispatcher, no un timeout casual |

## 4. Tres prioridades

1. Restituir la identidad de escenario en el camino público y comprobar sus
   consecuencias sobre subsemillas e IDs (A1).
2. Resolver explícitamente la contradicción de admisión del soporte de gaps
   (A2), sin cambiar silenciosamente la distribución ni los vectores RNG.
3. Cerrar las invariantes de memoria y sus oráculos independientes (A3/A4)
   antes de volver a marcar .3 aprobado.

## 5. Hallazgos

### A1 — Bloqueante: todos los layouts productivos usan el escenario literal `scenario_id`

- Ubicación: [`src/f5_layout.rs`](../../src/f5_layout.rs):257.
- Confianza: alta. Preexistente desde `c73a14e` (F5.1.1), demostrado con
  `git blame` y recorrido público ejecutado.
- Contrato: escenario real forma parte del material RNG y de la identidad F5
  (DTC-037/039/089; catálogos F5).
- Mecanismo: `ScenarioId::parse(value, scope)`, definido en
  `src/identity.rs:85`, recibe los argumentos invertidos:
  `ScenarioId::parse("scenario_id", config.scenario_id())`. El plan retiene la
  constante. En `src/f5_layout.rs:942` esa misma identidad entra en
  `F5IdentityContext`; las decisiones de longitud, posición y gaps consumen
  también el escenario del plan. La comprobación interna compara dos valores
  que ya comparten el error y por eso pasa.
- Evidencia: cambiar únicamente el ID en general y escenario, manteniendo
  semilla, referencia y ExecutionId, produjo:
  `scenario-config=scenario scenario-plan=scenario_id
  renamed-config=different-scenario renamed-plan=scenario_id
  layouts-with-ids-equal=true`.
  La igualdad incluye IDs de bloques/sitios. No se afirma que toda geometría
  deba cambiar al cambiar la semilla; sí deben cambiar el material de identidad
  y los ámbitos RNG contractuales.
- Corrección mínima: pasar el valor real como primer argumento; revisar las
  expectativas productivas que hayan congelado IDs derivados de la constante.
  No modificar catálogos para legitimar ese resultado.
- Prueba faltante: entrada pública con un ID que no sea `scenario_id`;
  igualdad exacta entre escenario configurado y plan; derivación independiente
  de subsemilla e ID desde el literal de entrada; segundo escenario con mismos
  restantes componentes. Los tests privados actuales construyen directamente
  un `ScenarioId` correcto y no detectan el error del adaptador.

### A2 — Grave: el soporte de gaps puede ejecutar un trabajo arbitrariamente mayor que el admitido

- Ubicación: [`src/f5_layout.rs`](../../src/f5_layout.rs):517.
- Confianza: alta. Camino de gaps incorporado en `08376f3`; defecto
  preexistente respecto al diff .3.
- Contrato: DTC-089 exige la selección literal creciente de enteros con
  Neumaier, también para gaps, y rechazar todo coste obligatorio antes del RNG.
- Mecanismo: `choose_weighted_integer_with_source` recorre todos los enteros
  de cada intervalo (primera pasada en línea 1956; segunda en 2000). El
  preflight solo incorpora `length_support_cardinality`, no la cardinalidad
  de gaps ni su multiplicación por los sorteos de sitios/intentos. La propia
  fórmula normativa también omite ese término: aquí hay una contradicción
  entre algoritmo y cota, no una simple línea Rust incorrecta.
- Evidencia pública, sin ejecutar el bucle enorme: longitud fija 2, B=L=A=Q=1,
  H=2, referencia ACGT. Gap `1..1`: 4380 operaciones y 71489 bytes; gap
  `1..1000000000000`: **4404 operaciones y los mismos 71489 bytes**. Se
  admite en ambos casos. El segundo requiere 10^12 visitas solo en el primer
  total Neumaier, muy por encima de 2^28.
- Corrección mínima propuesta: decisión técnica acotada sobre cómo reconciliar
  la fórmula con el selector literal. La opción que conserva su semántica es
  presupuestar conservadoramente los sorteos y la cardinalidad de soporte de
  gaps con aritmética comprobada. Cambiar a selección por bandas no es una
  corrección autorizada equivalente de RNG/binary64.
- Prueba faltante: perfiles con iguales filas, gap mínimo y dimensiones pero
  distinto ancho de soporte, incluidos soportes discontinuos; contador de
  visitas del selector y frontera exacta/−1 antes de abrir sus RNG.
- Detención: **requiere confirmar una enmienda técnica de la fórmula**. Este
  informe no la da por aprobada ni cambia diseño o código.

### A3 — Grave: el parser TSV puede asignar más memoria que el presupuesto antes de rechazar el formato

- Ubicación: [`src/f5_profile.rs`](../../src/f5_profile.rs):921.
- Confianza: alta. La materialización no acotada de columnas procede del
  parser previo y sigue dentro del nuevo camino de admisión .3.
- Contrato: arquitectura §15 y DTC-089 exigen comprobar crecimiento de
  búfer/vector antes de reservar; un TSV solo admite 2, 3 o 5 columnas.
- Mecanismo: `line.split('\t').collect::<Vec<_>>()` materializa todas las
  columnas antes de comprobar su número. En esta plataforma cada `&str`
  ocupa 16 bytes, mientras que la fase reserva solo `2*table_bytes` más
  filas persistentes. El número de separadores no está limitado por las
  columnas contractuales durante esa asignación.
- Evidencia pública: primer TSV de 400037 bytes con 200001 campos, tamaño y
  SHA correctos, límite de memoria 1000000. La carga llega a
  `GPS-FMT-101: número de columnas ... inválido` después de construir un
  vector cuyo almacenamiento mínimo es **3200016 bytes**. Esta cifra es
  tamaño de elementos, no una estimación de RSS ni una medición del allocator.
- Corrección mínima: descomponer con una estructura de tamaño fijo o iterador
  que solo acepte el número contractual y detecte el campo siguiente sin
  acumularlo. Reservar prospectivamente las estructuras variables restantes.
- Superficie de la misma invariante: líneas 641/703/748 construyen
  `normalize_weights/normalize_slice` antes de cobrar las filas de auditoría
  mediante `retain_profile_rows`; el rechazo puede llegar después de las
  asignaciones de los vectores normalizados. El cargo debe precederlas.
- Prueba faltante: primera fila con muchas columnas y presupuesto reducido,
  observador de reservas que demuestre tamaño acotado y cero asignación de
  auditoría si su cargo no cabe. La prueba no debe limitarse al código final
  del diagnóstico.

### A4 — Grave: el límite del metadata resta memoria histórica ya liberada

- Ubicación: [`src/f5_admission.rs`](../../src/f5_admission.rs):262.
- Confianza: alta. Nuevo método del diff .3; el problema persiste después del
  último informe de corrección.
- Contrato: memoria global es máximo de picos; headroom de la fase usa
  estructuras retenidas vivas, no el máximo histórico del proceso.
- Mecanismo: `remaining_input_bytes` resta
  `totals.working_memory_peak` del límite. `charge_memory_phase` calcula
  correctamente un máximo, pero ese arreglo no se propagó al derivador de
  límite físico. La prueba de bootstrap exacto deja 131784 como pico y límite;
  después de liberar los temporales, metadata sigue viendo cero.
- Evidencia: la suite falla con
  `phase=profile_metadata accumulated=131784 limit=131784`. La función solo
  almacena el pico, no la retención actual necesaria para esa resta.
- Corrección mínima: conservar la retención de la fase y derivar el headroom
  desde ella. No permitir superar el límite ni eliminar el coste del bootstrap.
- Prueba faltante: fase anterior dominante, transitorios liberados y perfil
  pequeño; la admisión exacta debe pasar y pico−1 debe rechazar. Añadir una fase
  con retención real no nula para impedir un arreglo que reste siempre cero.
- Límite relacionado todavía visible: `src/f5_preflight.rs:1102` aproxima
  “pico de tabla” como `2*mayor_TSV + todas_las_filas_finales`. La tabla mayor
  y las filas de tablas posteriores no necesariamente coexisten. Para cerrar
  v7 debe copiar/reconstruir cada fase real, no fabricar otro máximo a partir
  de escalares que perdieron el orden. Esto pertenece a la reconciliación de
  .3 con .4/.5, no justifica adelantar todo ese trabajo en esta corrección.

### A5 — Grave de validación: .3 figura cerrado con evidencia que no describe el árbol actual

- Ubicación: [plan F5](../../design/plan-implementacion-f5-0.1.md):188 y 200.
- Confianza: alta. Estado del diff actual.
- Contrato: pausa A exige pruebas exactas, cero apertura posterior y gates
  antes de continuar. “Completado” no sustituye la revisión aprobada.
- Evidencia: el texto mantiene ✅, Clippy correcto y suite roja “solo por H10”.
  Clippy falla y hay seis fallos adicionales. Además, los tests de memoria
  `src/f5_profile.rs:1335` y `:1414` obtienen el presupuesto llamando al
  productor sin límite y después usan ese mismo resultado como esperado.
  Comprueban autoconsistencia, no reconstrucción independiente. El test de
  identidad del perfil sí deriva tamaño/SHA desde bytes literales y es útil.
- Mecanismo: cambió el alcance del dispatcher y cambió la fórmula, pero el
  punto de observación, los fixtures y la acreditación de cierre permanecieron
  en versiones diferentes. Reescribir los literales fallidos hasta conseguir
  verde ocultaría esta causa.
- Corrección mínima: reabrir el estado documental de .3; clasificar y adaptar
  cada prueba según su contrato; registrar un cierre nuevo solo tras los
  gates sobre un árbol identificado. En publicación sustituir el retardo por
  TOML gigante por una sincronización que no viole el nuevo límite de entrada.
- Prueba faltante: oráculo de memoria con bytes/filas/retención por fase
  reconstruidos desde entradas, y trazado de todas las aperturas reales del
  camino público, no solo del opener TOML.

### A6 — Menor: un gap válido que no cabe se informa como fallo interno por overflow

- Ubicación: [`src/f5_layout.rs`](../../src/f5_layout.rs):1640.
- Confianza: alta; preexistente respecto a .3.
- Contrato: el siguiente gap que no cabe termina la lista de SNP; si quedan
  menos de dos, se descarta el intento y el agotamiento es `GPS-DOMAIN-106`.
- Evidencia: gap singleton `u64::MAX..u64::MAX`, longitud 2, B=L=1:
  perfil y preflight aceptan; layout produce `GPS-INTERNAL-101`, salida 70,
  “la posición SNP desborda u64”, sin mejor intento.
- Mecanismo: suma coordenada+gap antes de comparar con el fin. Aquí no es
  preciso representar la suma para saber que el gap no cabe.
- Corrección mínima: comparar gap con distancia restante comprobada y sumar
  solo cuando cabe. No añadir un límite biológico de gap inexistente.
- Prueba faltante: este singleton extremo por el camino público; no ejecuta
  soporte enorme y debe finalizar como intento insuficiente/agotado.

### A7 — Menor: el fallo de reserva deja utilizable el ledger

- Ubicación: [`src/f5_admission.rs`](../../src/f5_admission.rs):835.
- Confianza: alta; comportamiento preexistente de .2.
- Contrato: rechazo de recursos terminal, sin aperturas posteriores.
- Evidencia: `injected_buffer_reservation_failure_restores_ledger_totals`
  devuelve `GPS-RES-101` por reserva fallida y exige expresamente
  `phase == ProfileMetadata` en línea 992. El ledger puede avanzar después.
- Mecanismo: restaura totales y devuelve diagnóstico sin `reject`.
  Restaurar los totales y dejar la fase en `Rejected` son compatibles.
- Corrección mínima: conservar atomicidad de totales y estado terminal.
- Prueba faltante: el mismo fallo inyectado seguido de intento de transición
  y de lectura con espía; ambos deben rechazarse sin consumir el lector.
  Los llamadores actuales propagan el error, por lo que no se ha demostrado
  apertura posterior en la CLI; la gravedad queda limitada al contrato público.

Nota de diagnóstico de la misma frontera: `reject_resource` en línea 441
siempre usa `working_memory_peak` como límite, incluso cuando recibe
`OperationUnits` para overflow de tamaños de TSV. Debe escoger el límite de
la dimensión indicada; las unidades tipadas no bastan si el formatter las
vuelve a mezclar. No se eleva a un defecto de admisión porque sí rechaza.

## 6. Por qué se acumulan errores entre ciclos

Los siguientes son mecanismos observables, no juicios sobre la capacidad de
un modelo ni reconstrucciones inventadas de revisiones no conservadas.

**La unidad que se revisa ha sido menor que la invariante real.** El defecto A1
entró en F5.1.1 y sobrevivió a varias correcciones de admisión porque las pruebas
de RNG/ID reciben identidades correctas directamente y las de geometría crean
inputs privados. La unión “TOML → identidad → consumidor RNG” quedó sin
oráculo. No faltaba otra prueba del hash: faltaba comprobar el argumento que
llegaba al hash.

**Se ha corregido el lugar donde apareció el síntoma, pero no todos los
consumidores del concepto.** A4 muestra tres interpretaciones de memoria:
derivación del restante por pico histórico, cargo por fase y recomposición
final de escalares. Cambiar una no repara las otras. En el árbol actual también
hay dos funciones de aplicación que reconstruyen separadamente el ledger:
`validate_f5_profile_after_bootstrap` (1048) y
`materialize_f5_config_after_bootstrap` (1216). No es necesaria una arquitectura
nueva, pero sí una única transición compartida y propiedad clara del estado.

**Algunas pruebas son espejos del productor.** Los dos tests de presupuesto
del perfil ejecutan primero el productor para obtener el “esperado”. Eso deja
pasar una omisión si ambos recorridos la repiten. El test de reserva fallida
incluso congela el estado no terminal. En cambio, tamaño/SHA reconstruidos desde
bytes literales sí constituyen un oráculo útil. “Independiente” debe aplicarse
al resultado y a su reconstrucción, no solo a ejecutar otro test.

**El contrato también contiene una omisión.** A2 cumple la regla literal de
muestreo y puede seguir la fórmula escrita, pero esa fórmula no cubre su coste.
Una revisión limitada a comparar fórmula Rust con fórmula Markdown concluye
igualdad y pierde el problema. Hay que vincular cada término de coste a las
visitas efectivas del algoritmo, también en los rangos extremos.

**La migración atraviesa capas y pruebas aparentemente ajenas.** El dispatcher
de 65 KiB cambia la entrada v2 además de F5; el test de publicación usaba un
TOML grande para generar una ventana temporal. Añadir la validación de perfil
a la CLI invalida fixtures diseñados para tener perfil ausente. No demuestra
que el nuevo dispatcher sea conceptualmente incorrecto: demuestra que el
cierre exigía recorrer sus consumidores, no solo el target con prefijo F5.

**La evidencia mezcla generaciones de implementación.** `08376f3` concentra
5526 inserciones y 352 eliminaciones en 20 archivos, junto con evidencia de
v4/v6 y rediseño v7. `cierre-v4.md` reconoce que los informes detallados
H1–H3 y H6–H9 no están conservados. No se puede atribuirles causas o afirmar
que esta revisión redescubre exactamente esos fallos. El plan vigente sí
contiene una afirmación de cierre refutada por los gates actuales. Esa falta
de un estado validado identificable alimenta los ciclos.

### Complejidad necesaria y complejidad evitable

Son necesarias las tres unidades de recursos, identidades físicas por rol,
SHA sobre el mismo handle, separación de transitorios/retención, selección
RNG normativa y admisión antes de acciones posteriores. Quitarlas para reducir
código eliminaría protecciones explícitamente acordadas.

Es evitable mantener varias reconstrucciones del ledger y de los picos, usar
parámetros sueltos allí donde ya existe un estado de fase, mantener pruebas de
bootstrap conectadas a una función que hace más trabajo y tratar cada nueva
fórmula como un parche aislado. El selector literal caro no es por sí solo
“sobrediseño”: está prescrito. Su coste debe incorporarse al contrato o
modificarse mediante decisión explícita.

La evidencia no justifica concluir que todo F5 necesite reescribirse.
Sí justifica sustituir las reparaciones por síntoma por un cierre pequeño,
trazado por invariantes y con un único recorrido productivo por transición.

## 7. Estrategia acotada de cierre propuesta

Esta estrategia no autoriza implementación en este turno.

1. Identificar un árbol de revisión y reabrir .3 documentalmente. Mantener una
   tabla por hallazgo con estado, cambio mínimo, regresión y gate. No crear una
   nueva numeración que pierda estos mecanismos.
2. Corregir A1 como unidad independiente, con oráculo público de escenario,
   subsemilla e ID. No mezclarlo con memoria ni con el receipt.
3. Decidir el término de coste de gaps de A2, conservar semántica RNG y añadir
   contador instrumental. Corregir A6 en esa misma frontera geométrica.
4. Cerrar A3/A4/A7 alrededor de una única propiedad del ledger: “qué está vivo,
   qué se reserva ahora y qué ya terminó”. Reservar antes de asignar; derivar
   headroom desde retención; reusar la transición común en CLI y biblioteca.
   Esto no exige implementar FASTA/receipt por adelantado.
5. Construir una matriz corta de aceptación: bootstrap dominante, metadata
   dominante, cada tabla dominante, retención dominante, error de identidad,
   exceso de columnas y reserva fallida. Para cada caso: resultado independiente,
   frontera exacta/−1, primera apertura/primera asignación y cero acciones
   posteriores. Congelar estos oráculos antes de aceptar las correcciones.
6. Adaptar las pruebas obsoletas según la clasificación de §3, incluida
   publicación. Exigir Clippy y suite dirigida verdes en A; conservar únicamente
   H10 como deuda explícita de .6, nunca afirmar suite completa verde.
7. Una revisión independiente de A sobre ese estado. Si reaparece una invariante,
   corregir su propietario y su oráculo; no reabrir módulos sin relación ni
   iniciar otra revisión integral por cada literal cambiado.
8. Solo después de A, continuar .4/.5 y sus pausas ya aprobadas. .6 debe cerrar
   el camino CLI/biblioteca → v7 sellado → layout con fórmulas y trazado
   independientes, JSON Schema coherente y suite completa verde.

Condición objetiva de salida: cada invariante tiene una implementación dueña,
un caso negativo que falla antes del arreglo, un oráculo independiente y una
evidencia sobre el mismo árbol. No se promete ausencia absoluta de defectos.

## 8. Aspectos correctamente resueltos y límites

El dispatcher lee y clasifica después de EOF bajo 65536+byte de exceso; las
pruebas de identidad TOML reconstruyen hashes desde bytes escritos. El perfil
usa un único handle para lectura/hash/medida, reserva juntos los tamaños TSV
declarados, libera los cinco buffers secuencialmente y conserva tipos y
auditorías. El parser valida continuidad/disjunción, dominios de probabilidades,
rangos de haplotipos y cobertura LD por población. La geometría mantiene
transacciones de intento, orden canónico, prefijos globales de inicios y
selección ALT distinta de REF. Los contratos heredados cuentan con pruebas de
inmutabilidad. Nada de ello compensa A1–A5, pero permite acotar la reparación.

No se ejecutó el bucle de 10^12 gaps, ni se midió RSS, ni se completaron fases
futuras. Las comparaciones de memoria del TSV se basan en su número de elementos
y `size_of::<&str>()`, más la ruta de asignación leída. Las reproducciones
siguientes usan APIs públicas de producción; reutilizan únicamente constructores
de fixtures de los tests, no sus oráculos de resultados.

## 9. Reproducciones adicionales conservadas

Se creó el siguiente programa en `/tmp/gps_f5_adversarial_probe.rs`. El
`include!` reutiliza constructores de fixtures; al compilar sin `--test` no
ejecuta los tests anotados ni habilita los inputs privados de layout. Los
archivos de prueba se crean bajo `/tmp`; no modifica el repositorio.

Compilación y ejecución efectivamente realizadas, después de los gates
compartidos y con Cargo libre:

```text
rustc -Awarnings --edition=2024 /tmp/gps_f5_adversarial_probe.rs -L dependency=target/debug/deps --extern genome_population_simulator=target/debug/deps/libgenome_population_simulator-862cc19d977f35af.rlib --extern sha2=target/debug/deps/libsha2-9e34213456fdf3c2.rlib --extern serde_json=target/debug/deps/libserde_json-e9c890ee917c05ca.rlib -o /tmp/gps_f5_adversarial_probe
/tmp/gps_f5_adversarial_probe
```

Los sufijos de `rlib` pertenecen a este árbol/toolchain; tras recompilar se
localizan con `rg --files target/debug/deps`. No usar una biblioteca antigua
para verificar una corrección. Resultado observado:

```text
singleton: gap_max=1 operation_units=4380 memory=71489
scenario-config=scenario scenario-plan=scenario_id
renamed-config=different-scenario renamed-plan=scenario_id
layouts-with-ids-equal=true
huge: gap_max=1000000000000 operation_units=4404 memory=71489
extra-columns: table_bytes=400037 field_count=200001 slice_size=16
fields_min_bytes=3200016 declared_memory_limit=1000000 error=GPS-FMT-101
singleton-u64-max-gap: GPS-INTERNAL-101, la posición SNP desborda u64
```

Fuente completa:

```rust
include!("/home/ximo/devel/genome_population_simulator/tests/f5_resource_preflight.rs");

fn main() {
    let root = temporary("adversarial-gap-support");
    for (label, gap) in [("singleton", "1"), ("huge", "1000000000000")] {
        let general = custom_configuration_with_gap(
            &root.join(label),
            "",
            "ACGT",
            "length_min_bp\tlength_max_bp\tweight\n2\t2\t1\n",
            &format!("gap_min_bp\tgap_max_bp\tweight\n1\t{gap}\t1\n"),
            "haplotype_rank\tweight\n0\t0.5\n1\t0.5\n",
            "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t1\t0.5\t0.05\n",
        );
        let config = load_and_validate_f5_config(&general).expect("config admits gap support");
        let preflight = preflight_f5_resources(&config, F5OtherResources::new(0, 0, 0)).expect("preflight admits gap support");
        println!("label={label} gap_max={gap} operation_units={} memory={} general={}", preflight.operation_units(), preflight.working_memory_bytes(), general.display());
        if label == "singleton" {
            let plan = F5BlockLengthPlan::from_preflight(&config, &preflight).unwrap();
            let execution = ExecutionId::parse("RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4", "test").unwrap();
            let first = build_f5_physical_layout(F5LayoutRequest::new(&config,&preflight,&execution).unwrap()).unwrap();
            let scenario_path = general.parent().unwrap().join("scenarios/scenario.toml");
            for file in [&general, &scenario_path] {
                let original = fs::read_to_string(file).unwrap();
                fs::write(file, original.replace("scenario_id = \"scenario\"", "scenario_id = \"different-scenario\"")).unwrap();
            }
            let renamed = load_and_validate_f5_config(&general).unwrap();
            let renamed_preflight = preflight_f5_resources(&renamed,F5OtherResources::new(0,0,0)).unwrap();
            let renamed_plan = F5BlockLengthPlan::from_preflight(&renamed,&renamed_preflight).unwrap();
            let second = build_f5_physical_layout(F5LayoutRequest::new(&renamed,&renamed_preflight,&execution).unwrap()).unwrap();
            println!("scenario-config={} scenario-plan={} renamed-config={} renamed-plan={} layouts-with-ids-equal={}",config.scenario_id(),plan.scenario_id().as_str(),renamed.scenario_id(),renamed_plan.scenario_id().as_str(),first==second);
        }
    }
    let malicious_row = format!("{}1\n", "1\t".repeat(200_000));
    let block_table = format!("length_min_bp\tlength_max_bp\tweight\n{malicious_row}");
    let general = custom_configuration_with_gap(
        &root.join("extra-columns"),
        "max_working_memory_bytes = 1000000\n",
        "ACGT",
        &block_table,
        "gap_min_bp\tgap_max_bp\tweight\n1\t1\t1\n",
        "haplotype_rank\tweight\n0\t0.5\n1\t0.5\n",
        "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t1\t0.5\t0.05\n",
    );
    let error = load_and_validate_f5_config(&general).expect_err("extra columns");
    println!("extra-columns table_bytes={} field_count=200001 slice_size={} fields_min_bytes={} declared_memory_limit=1000000 error={}", block_table.len(), std::mem::size_of::<&str>(), 200001 * std::mem::size_of::<&str>(), error.render_human());
    let general = custom_configuration_with_gap(
        &root.join("gap-overflow"), "", "ACGT",
        "length_min_bp\tlength_max_bp\tweight\n2\t2\t1\n",
        "gap_min_bp\tgap_max_bp\tweight\n18446744073709551615\t18446744073709551615\t1\n",
        "haplotype_rank\tweight\n0\t0.5\n1\t0.5\n",
        "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t1\t0.5\t0.05\n",
    );
    let config = load_and_validate_f5_config(&general).expect("singleton maximum gap accepted");
    let preflight = preflight_f5_resources(&config, F5OtherResources::new(0,0,0)).expect("singleton maximum gap preflight");
    let execution = ExecutionId::parse("RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4", "test").unwrap();
    let request = F5LayoutRequest::new(&config, &preflight, &execution).unwrap();
    println!("singleton-u64-max-gap layout={:?}", build_f5_physical_layout(request));
}

```

Gates para volver a verificar, una vez autorizadas y aplicadas las correcciones:

```text
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --test f5_bootstrap_dispatch --test f5_bounded_input --test f5_resource_preflight --test cli_config_diagnostics
cargo test --lib f5_profile
cargo test --test f5_identity --test f5_rng_catalog --test f5_0_1_config_admission --test f5_0_1_external_schemas
cargo test --test publication_process complete_failed_and_incomplete_are_real_minimal_trees -- --exact --nocapture
cargo test --all-targets --all-features --no-fail-fast
cargo doc --no-deps
git diff --check
```

