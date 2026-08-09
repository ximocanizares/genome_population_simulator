# Revisión de F1.2

- Fecha: 9 de agosto de 2026.
- Alcance: B1–B12, sin C1/F1.3.
- Método: revisión completa del diff conforme a `revisar-codigo-rust`, separada
  en sintaxis, semántica, dominancia y diagnósticos; correcciones conforme a
  `aplicar-correcciones-rust` y repetición de todas las puertas.

## Resultado por frente

- **Sintaxis y admisión:** la sintaxis TOML se comprueba antes de clase y
  versión; clase y versión preceden a la semántica. Se conservan bytes, ruta y
  rangos. Se revisaron desconocidos, duplicados, cabeceras cruzadas y los
  envoltorios válidos directo, evolutivo y de familias. Sin hallazgos pendientes.
- **Semántica:** se revisaron identidades, referencias, escenarios, eventos
  explícitos, ploidía, sumas comprobadas, formas de frecuencias, reglas TOML y
  `EffectiveConfig`. La configuración efectiva es poseída, no vuelve a leer el
  TOML y conserva valores normalizados, ganador, selectores y procedencia. Sin
  hallazgos pendientes.
- **Dominancia:** se inspeccionó el retículo manual, la ausencia de prioridad
  entre ejes, las reglas permutadas, la regla combinada y las máximas
  incomparables aun con valores iguales. Duplicados y solapes cerrados se
  rechazan antes de resolver. Sin hallazgos pendientes.
- **Diagnósticos:** humano y JSON parten del mismo `Diagnostic`; el JSON es un
  documento único, escapa controles, satisface el esquema independiente y
  conserva código, categoría, mensaje, contexto y salida 2. Un rechazo no crea
  el destino solicitado. Sin hallazgos pendientes.

## Hallazgos y correcciones

| ID | Prioridad | Hallazgo confirmado | Corrección y regresión |
| --- | --- | --- | --- |
| R1 | media | Las fases declaradas dentro de escenarios no alimentaban el catálogo de identidades de los selectores. | Se incorporaron las fases de escenarios a la validación; el fixture evolutivo contiene una regla por `phase-1` y pasa por CLI. |
| R2 | media | El resolvedor exigía regla base incluso cuando la ploidía declarada ya era un valor general válido, y una sustitución ganadora no actualizaba la secuencia efectiva. | Se añadió fallback explícito solo donde existe valor declarado, resolución automática por secuencia y prueba de ploidía efectiva `4` con ganador y procedencia. |
| R3 | baja | Algunos `GPS-CONF-008` de validación cruzada omitían la ruta requerida por el catálogo. | Todos se construyen ahora con la ruta general o de escenario; la integración humano/JSON comprueba su conservación. |
| R4 | baja | La escritura con `println!` podía producir `panic` ante una tubería cerrada. | La frontera CLI usa escritura comprobada y devuelve 70 si falla el canal público. |
| R5 | media | El catálogo materializado omitía la familia normativa de recombinación. | Se añadió `recombination_rate`, familia, ejes grupo/secuencia, prueba exhaustiva y nuevo fingerprint previo a publicar v1. |
| R6 | baja | Los ejemplos no ejercitaban de extremo a extremo las dos rutas con escenario ni una referencia suministrada válida. | Se añadieron tres fixtures generales y se validan por CLI sin publicar salida. |

Balance: **6 confirmados, 6 corregidos, 0 descartados y 0 pendientes**.

## Segunda revisión: hallazgos H-F12

La revisión posterior cuestionó de forma adversarial la admisión completa de
intersecciones, los dominios de valores de regla y la evidencia de los 36
`TEST-REQ-*`. Los tres hallazgos se confirmaron antes de editar y se aplicaron
con regresiones previas a la corrección.

| ID | Prioridad | Decisión inicial | Corrección y regresión | Estado |
| --- | --- | --- | --- | --- |
| H-F12-01 | grave | Aplicar | La validación previa recorre las intersecciones compatibles de reglas declaradas y resuelve contra el conjunto completo. Se añadieron fixtures E2E permutados para población+fase, una regla combinada dominante y pruebas directas de individuo+región y región+subpoblación; el diagnóstico conserva parámetro, ámbito y todas las máximas de la intersección resuelta. | Corregido |
| H-F12-02 | grave | Aplicar | El catálogo fija un dominio por `ParameterKey`; la admisión convierte solo enteros decimales admitidos como reales y `Rule` conserva campos privados para que ningún consumidor salte `Rule::new`. Ploidía acepta enteros 1–10, las tres probabilidades aceptan `f64` finitos en `[0,1]`, salida exige booleano y tasa de variantes, recombinación y cobertura aceptan `f64` finitos no negativos. Doce casos inválidos se ejecutan en CLI humano/JSON y otra regresión cubre todos los bordes inclusivos. | Corregido |
| H-F12-03 | grave | Aplicar con reasignación | Se auditaron literalmente las 36 filas. Solo 16 REQ completos pertenecen y están demostrados en F1.2; 20 se reasignaron a F2.2/F2.3/F4.2/F6.2/F7.6/F9.7. La prueba de matriz verifica conjunto, estados, destinos existentes en la hoja de ruta, toda evidencia enlazada y coherencia con la trazabilidad maestra. | Corregido |

No se inventó prioridad entre ejes ni un límite superior para tasas o cobertura
que la especificación no fija. La tasa de variantes, la recombinación y la
cobertura solo reciben el dominio mínimo aprobado de número real finito no
negativo. La revisión completa posterior de sintaxis, semántica, dominancia,
errores, pruebas, dependencias y diff no encontró hallazgos nuevos.

## Reapertura de H-F12-03 y encapsulación

La revisión independiente rechazó la primera corrección de H-F12-03 porque aún
quedaban una cifra antigua y evidencias parciales. Se confirmaron y corrigieron
los tres residuos. Además se confirmó un hallazgo menor nuevo en H-F12-02:

| ID | Prioridad | Decisión inicial | Estado final | Evidencia |
| --- | --- | --- | --- | --- |
| H-F12-03-R1 | grave | Aplicar | Aplicado | El cierre de F1 usa 72 requisitos vigentes y conserva 20 reasignados; no queda la cifra 92 como asignación vigente. |
| H-F12-03-R2 | grave | Aplicar | Aplicado | Nueva regresión resuelve generación>fase, subpoblación>población, locus>región>secuencia, individuo/evento/grupo/bloque>ausencia y demuestra que dos ejes siguen incomparables. |
| H-F12-03-R3 | grave | Aplicar con adaptación | Aplicado con adaptación | Reauditoría literal de las 23 filas antes marcadas `pass`: 16 permanecen superadas y 7 se reasignan porque una prueba de admisión no puede cerrar consumo, no herencia, consolidación familiar o tolerancia 0,5/C. |
| H-F12-04 | menor | Aplicar | Aplicado | Los cinco campos de `Rule` son privados y solo se exponen mediante accesores inmutables; toda construcción pasa por `Rule::new`. |

Resultado de la reauditoría literal:

| REQ antes `pass` | Resultado | Razón o evidencia principal |
| --- | --- | --- |
| CONF-01 | `pass` | Admisión del TOML estructurado y conservación exacta de bytes. |
| CONF-04 | `pass` | Prueba nueva de cantidad, longitudes, clases y composición de referencia simulada. |
| CONF-05 | `pass` | Suma conocida y desbordamiento comprobado. |
| CONF-07 | `pass` | Tres vías válidas por CLI. |
| CONF-09 | `pass` | Enlaces de escenario evolutivo y familiar; ausencia en vía directa. |
| CONF-10 | `pass` | Ficheros de escenario estructurados con identidad y orden. |
| CONF-14 | `pass` | Editar únicamente TOML cambia el valor efectivo. |
| CONF-17 | F9.7 | Falta cierre conjunto de individuos directos, fundadores y parentales. |
| CONF-18 | F9.7 | Existencia/categoría está probada, pero no consumo, no herencia y no retroactividad sobre descendientes reales. |
| CONF-29 | `pass` | Jerarquía espacial literal en contexto aplicable y sin prioridad transversal. |
| CONF-30 | `pass` | Solapes cerrados parciales, totales y de borde con igual valor. |
| CONF-31 | `pass` | Intersección espacial+individual conflictiva y regla combinada dominante bajo permutación. |
| CONF-65 | `pass` | Ausente, no entero, menor que uno y mayor que diez con referencia y secuencia. |
| CONF-106 | `pass` | Borde inclusivo, normalización/rechazo y registro exacto antes/después. |
| CONF-116 | F9.7 | Solo puede auditarse completamente cuando existan todos los parámetros v1. |
| CONF-117 | `pass` | Todas las jerarquías internas literales y ausencia de prioridad entre ejes. |
| CONF-118 | `pass` | Máximas incomparables, diagnóstico, combinación dominante y orden permutado. |
| CONF-119 | `pass` | Duplicados exactos con igual o distinto valor. |
| CONF-120 | F9.7 | El rechazo de selector existe, pero no el cierre sobre descendientes y retroactividad reales. |
| CONF-121 | F9.7 | Requiere el catálogo completo de parámetros de todos los hitos v1. |
| CONF-122 | F9.7 | Requiere demostrar que la técnica no altera verdad o trayectoria y que eventos biológicos no se reinterpretan. |
| VAL-40 | F6.2 | Falta la tolerancia específica 0,5/C y su separación del normalizador general. |
| VAL-51 | `pass` | Evidencia múltiple real para base, jerarquías, combinada, conflictos, duplicados, solapes, prohibiciones y frontera biológica/técnica. |

Balance de la reapertura: **4 confirmados, 4 aplicados, 0 descartados y 0
pendientes**.

## Aprobación de la pausa

La revisión independiente final del 9 de agosto de 2026 repitió formato,
Clippy, las 63 pruebas, pruebas documentales, documentación, construcción
release, metadatos, grafo bloqueado, JSON Schema, `cargo audit` y
`git diff --check`, sin hallar defectos bloqueantes o graves pendientes. F1.2
queda aprobada y se autoriza comenzar F1.3 en un cambio posterior.
