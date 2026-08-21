# Revisión adversarial del plan de implementación F2 0.3

> **Informe histórico.** Revisa el plan de 130 pasos sustituido por F2 0.2 y no
> es autoridad para reintroducir expedientes, protocolos o hardening local.

## 1. Alcance y método

- **Documento revisado:**
  [plan-implementacion-f2-0.1.md](./plan-implementacion-f2-0.1.md), en el estado
  del árbol compartido del 15 de agosto de 2026, después de DTC-070 a DTC-075
  y de las correcciones declaradas para RAF2-028 a RAF2-044.
- **Tipo de revisión:** adversarial, repetida desde cero; sus correcciones
  técnicas se incorporan a los contratos y al plan, sin modificar código.
- **Fuentes contrastadas:** especificaciones, arquitectura, decisiones
  técnicas, hoja de ruta, trazabilidad, registro de ampliaciones futuras,
  plan F1, punto de reanudación, informes adversariales F2 0.1/0.2, catálogo
  RNG y fronteras actuales de código/CLI cuando afectan a los prerrequisitos.
- **Escala:** un hallazgo `bloqueante` impide ejecutar el plan con un contrato
  único; uno `mayor` obliga a corregir el plan antes de aprobarlo; uno `menor`
  no invalida por sí solo el recorrido, pero debe corregirse para evitar una
  evidencia o mantenimiento ambiguos.

La revisión distingue el diseño que F2 todavía debe implementar del estado
actual del producto. Por ello, que el repositorio conserve por ahora el
catálogo RNG v3 y una CLI de F1 no se considera un defecto: F2.2.8 y los pasos
posteriores son precisamente quienes deben incorporar v4 y las órdenes nuevas.

## 2. Veredicto ejecutivo

**No aprobar todavía el plan F2.**

El plan corregido supera sus comprobaciones mecánicas: contiene 128 pasos, todos tienen
los cinco campos obligatorios, existen seis pausas y el inventario separa 79
requisitos exclusivos, nueve compartidos y nueve parciales directos. También
han quedado cerrados de forma convincente `GPS_CASE`, el modo `exact`, el
catálogo de capacidades, la evidencia real de plataforma y la receta compacta
de perfiles humanos.

La revisión encontró cuatro defectos propios del plan, ya corregidos:

1. DTC-072, el esquema y el perfil SFS expresaban reglas distintas para 0/1 y
   el cálculo lineal podía redondear a cero.
2. El plan asignaba a indels el selector estático, incompatible con reenumerar
   un estado que puede ganar o perder oportunidades.
3. F2.2.2 y F2.2.8 compartían propiedad de `log-uniform-v1`, con una dependencia
   oculta y demasiadas responsabilidades en una sola casilla.

Los dos primeros eran bloqueantes porque dos implementaciones razonables podían
aceptar entradas distintas o producir resultados distintos con la misma
semilla. El tercero es mayor porque deja una dependencia oculta, un productor
compartido sin dueño único y una unidad que no se puede revisar aisladamente.
Ninguno requirió una decisión biológica nueva. DTC-075, DTC-068/071/072
corregidas y F2.2.8.1–F2.2.8.6 les dan un contrato único y verificable.

F1 continúa siendo el único bloqueo abierto: la revisión independiente no la
ha aprobado. Aunque estos hallazgos de F2 ya estén corregidos, F2.0.1 no puede
iniciarse hasta disponer de aprobación y commit base F1 identificable.

## 3. Comprobaciones mecánicas y de inventario

| Comprobación | Resultado | Valoración |
|---|---:|---|
| Pasos jerárquicos F2 | 128 | Correcto tras dividir F2.2.8 |
| Pasos con exactamente `Depende de`, `Procede de`, `Cambia`, `Se comprueba con` y `Evidencia ejecutable` | 128/128 | Correcto |
| Referencias simples o extremos de rango F2 inexistentes en `Depende de` | 0 | Correcto |
| Pausas `PAUSA-F2.0` a `PAUSA-F2.5` | 6 | Correcto |
| REQ de cierre exclusivo | 79 = 14+4+22+12+2+25 | Correcto |
| REQ compartidos | 9 | Correcto |
| REQ parciales directos | 9 | Correcto |

Las nueve obligaciones parciales directas son `REQ-SAL-02`, `REQ-SAL-12`,
`REQ-FMT-07`, `REQ-VER-23`, `REQ-VER-24`, `REQ-VER-25`, `REQ-VER-26`,
`REQ-MAN-31` y `REQ-MAN-58`. F2.4.8 cubre la primera y F2.5.12 las otras ocho,
siempre con estado `partial` y destino posterior. Las nueve compartidas tienen
productores y consumidores reales señalados en §15.1; no se confunden con los
79 cierres exclusivos.

El primer paso de cada hito depende de la pausa anterior. Los pasos finales de
Marimo dependen de todo el hito biológico correspondiente, por lo que su
referencia adicional a `PAUSA-F2.0` no permite eludir las pausas intermedias.

## 4. Revisión de RAF2-028 a RAF2-040

| Hallazgo previo | Estado tras esta revisión | Evidencia principal |
|---|---|---|
| RAF2-028, catálogo RNG incompleto/circular | Cerrado en su causa original | DTC-071 fija v4 como superconjunto de v3, once ámbitos, componentes ordenados y ALT anterior a `variant_id`. Los defectos de integración encontrados se registran aparte en RAF2-043. |
| RAF2-029, asignación `exact` indeterminada | Cerrado | DTC-070 separa cuotas, permutación por población y flujo de valor por población×variante; prohíbe asociar por coordenada. |
| RAF2-030, `GPS_CASE` ambiguo | Cerrado | DTC-065/066, verdad y F2.5 fijan `Number=.`, IDs originales únicos/canónicos y relación ALT↔casos sin caso compuesto. |
| RAF2-031, parciales omitidos | Cerrado | §15.2 contiene los nueve expedientes, pasos productores y destino final. |
| RAF2-032, pausas no representadas | Cerrado | Las seis pausas existen y el grafo obliga a atravesarlas en orden. |
| RAF2-033, algoritmos RNG no cerrados | Cerrado | DTC-072/075 fijan fronteras, selección dinámica y productores únicos; RAF2-041/042/043 quedan resueltos. |
| RAF2-034, capacidades anticipadas | Cerrado | DTC-073 separa `catalogued` y `available`; cada flag se habilita solo junto con parser, validación y semántica. |
| RAF2-035, plataforma autoetiquetada | Cerrado | DTC-074 exige hechos detectados, rechaza discordancias y separa WSL2, Ubuntu nativo y actas manuales. |
| RAF2-036, pasos demasiado grandes | Cerrado | F2.2.8 conserva su ID padre y se divide en seis hijos revisables con evidencia propia. |
| RAF2-037, productores compartidos incorrectos | Cerrado para la matriz §15.1 | Los nueve compartidos enlazan consumidores efectivos. La propiedad interna de `log-uniform-v1` es un nuevo defecto de planificación, no una fila REQ compartida. |
| RAF2-038, estado contradictorio RAF2-001 | Cerrado | Cabecera, pausa y cierre mantienen RAF2-001 resuelto. |
| RAF2-039, DTPD-103 resuelta/pendiente | Cerrado | El plan consume DTC-062 a DTC-066 sin reabrir la decisión. |
| RAF2-040, «fichero D» ambiguo | Cerrado | El plan nombra ya los registros de modificaciones/procedencia de `truth.jsonl`. |

## 5. Hallazgos bloqueantes propios de F2

### RAF2-041 — Los extremos del espectro alélico tienen contratos incompatibles

- **Estado:** resuelto el 15 de agosto de 2026 mediante DTC-068/072/075 y F2.2.8.4.
- **Gravedad original:** bloqueante técnico.
- **Ámbito:** DTC-068, DTC-072, especificaciones de espectro y pasos
  F2.2.12, F2.2.14.2 y F2.2.24.
- **Evidencia:** el contrato general define intervalos `[min_af,max_af)` y
  afirma que solo una fila `fixed` produce exactamente 0 o 1. DTC-072 repite
  que un intervalo continuo nunca alcanza su extremo superior. En cambio, las
  especificaciones, DTPD-121 y F2.2.14.2 escriben el último tramo del perfil
  como `[0.05,1]` y dicen que «puede alcanzar 1». Es una notación cerrada y una
  afirmación incompatible con el algoritmo normativo.
- **Defecto numérico adicional:** `linear-uniform-v1` admite cualquier
  `[min,max)` finito con `min < max`. Si `min=0` y el intervalo es
  suficientemente pequeño, `u × (max-min)` puede redondear a 0 aunque
  `open-unit-f64-v1` produzca `u>0`. La corrección definida solo protege el
  extremo superior; por tanto, la regla «solo `fixed` produce exactamente 0»
  no se cumple para todo el dominio admitido.
- **Riesgo:** un implementador puede tratar el último tramo como cerrado,
  exclusivo o especial; otro puede rechazar el perfil por contradecir el
  esquema. Los vectores no pueden certificar simultáneamente ambas reglas.
- **Resolución aplicada:** todos los intervalos son `[min,max)`, incluido
  `[0.05,1)`; `max=1` cubre sin producir 1 y 0/1 exactos exigen `fixed`. La
  corrección inferior de `linear-uniform-v1` devuelve el menor `f64` positivo
  si el redondeo diera cero. F2.2.8.4 congela vectores para cero, subnormales,
  `next_down(1)` y ambos redondeos; DTC-068 conserva pesos 64:12:8.
- **Decisión del usuario:** no necesaria. La intención aprobada —0 y 1 solo en
  filas `fixed`— ya permite escoger la corrección técnica compatible.

### RAF2-042 — La selección estática sin reemplazo no define los *indels* sobre estado vigente

- **Estado:** resuelto el 15 de agosto de 2026 mediante `dynamic-indel-selection-v1`.
- **Gravedad original:** bloqueante técnico.
- **Ámbito:** DTC-020/072, ACE-107, arquitectura del motor de mutaciones y
  F2.2.30–F2.2.31.
- **Evidencia:** `weighted-without-replacement-v1` recibe una colección,
  comprueba antes de consumir que contiene suficientes candidatos positivos y
  después solo retira los elegidos. El contrato de *indels* exige lo contrario:
  sortear una sola vez el recuento, elegir un evento, modificar la secuencia y
  volver a enumerar las oportunidades sobre el nuevo estado antes del evento
  siguiente. Una inserción puede crear anclas nuevas; una deleción puede
  eliminar o desplazar candidatas. F2.2.30 exige esa actualización, mientras
  F2.2.31 vuelve a nombrar el consumidor estático y una precomprobación de
  insuficiencia sin definir cómo se combinan ambas reglas.
- **Riesgo:** preflightar el espacio inicial puede rechazar una realización que
  habría ganado oportunidades, o aceptar una cuya colección deja de ser
  válida tras la primera edición. Volver a llamar con `k=1` tampoco es una
  solución normativa mientras no se definan la identidad canónica de
  candidatas entre estados, el consumo continuo del flujo y el momento exacto
  del error por agotamiento.
- **Resolución aplicada:** `weighted-without-replacement-v1` queda reservado al
  conjunto SNV estático. `dynamic-indel-selection-v1` preflighta solo
  invariantes estáticas, reenumera/pondera/aplica por ordinal con el mismo flujo
  `indel_position` y falla transaccionalmente al agotarse, sin nuevo Poisson ni
  prefijo publicado. Oráculos independientes cubren creación y pérdida de
  oportunidades, orden, consumo cero y agotamiento temprano.
- **Decisión del usuario:** no necesaria. Es la materialización técnica de
  ACE-107 y de la decisión ya aprobada de actualizar candidatos evento a
  evento.

## 6. Hallazgos mayores

### RAF2-043 — La familia RNG tiene propiedad duplicada, una dependencia oculta y un paso sobredimensionado

- **Estado:** resuelto el 15 de agosto de 2026 con F2.2.8.1–F2.2.8.6.
- **Gravedad original:** mayor.
- **Ámbito:** F2.2.2 y F2.2.8, con consumidores F2.2.24 y F2.3.5–F2.3.8.
- **Evidencia:** F2.2.2 «implementa [...] la matemática compartida por
  `log-uniform-v1`» y verifica sus vectores. F2.2.8 dice que «implementa la
  familia DTC-072 antes de cualquier consumidor», lo que vuelve a incluir
  `log-uniform-v1`, pero no depende de F2.2.2. A su vez, F2.2.24 usa
  `log-uniform-v1` sin depender de F2.2.2. El grafo permite, por tanto, aprobar
  la familia y ejecutar el consumidor antes de que exista la matemática
  compartida que el propio plan asigna a otro paso.
- **Tamaño:** F2.2.8 publica el catálogo v4, conserva seis filas/vectores
  anteriores, añade once ámbitos, deriva BLAKE3, implementa diez algoritmos,
  fija consumos/fronteras y aporta oráculos de aislamiento. Esta unidad es más
  amplia que varias agrupaciones que RAF2-036 obligó a subdividir.
- **Riesgo:** dos módulos pueden terminar implementando el mismo algoritmo; la
  evidencia de un paso puede depender de código aún no aprobado; una casilla
  verde no permite saber qué primitiva, consumidor o vector fue revisado.
- **Resolución aplicada:** F2.2.8 es padre de seis pasos para catálogo,
  uniforme/entero, selectores, continuas, genotipo/fase y oráculos. F2.2.2 solo
  posee Poisson y evalúa `libm`; F2.2.8.4 es el único productor de
  `log-uniform-v1` y depende de esa evaluación. Los consumidores citan los
  hijos necesarios y F2.2.8.6 verifica que no haya productores duplicados.
- **Decisión del usuario:** no necesaria; es una corrección de propiedad,
  tamaño y orden de implementación.

## 7. Hallazgos menores

### RAF2-044 — El paso que cuenta candidatos no nombra su ámbito RNG normativo

- **Estado:** resuelto el 15 de agosto de 2026 en F2.2.11/25/39 y F2.3.14.
- **Gravedad original:** menor.
- **Ámbito:** F2.2.11 y DTC-071.
- **Evidencia:** DTC-071 reserva `matrix_candidate_count` para el Poisson de
  `mean_loci`/densidad, pero F2.2.11 solo dice que «sortean Poisson» y su campo
  `Procede de` no cita DTC-071. La relación se puede deducir de F2.2.8 y del
  catálogo, pero no queda afirmada en el paso consumidor.
- **Riesgo:** la prueba local puede cubrir el algoritmo Poisson sin demostrar
  que usa el ámbito v4 correcto ni que `exact_loci` no consume ese flujo.
- **Resolución aplicada:** F2.2.11 cita DTC-071/075, usa
  `matrix_candidate_count` solo para media/densidad y exige consumo cero en
  `exact_loci`; publicación, resultado F2.2 y manifiesto registran algoritmo,
  ámbito/subsemilla/consumo o su ausencia normativa. El oráculo prueba
  aislamiento frente a selección, ALT y frecuencias.
- **Decisión del usuario:** no necesaria.

## 8. Decisiones que requieren al usuario

**Ninguna.**

Esta revisión no descubre una alternativa biológica razonable que deba volver
a decidirse. Las reglas aprobadas ya indican que:

- los extremos exactos 0/1 proceden de filas `fixed`;
- los eventos aleatorios actualizan su espacio después de cada aceptación;
- un algoritmo o artefacto compartido debe tener un productor único y una
  dependencia explícita.

Las cuatro RAF2 nuevas pueden corregirse sin ampliar el alcance ni cambiar esas
decisiones.

## 9. Aspectos comprobados como correctos

- **Catálogo RNG v4:** DTC-071 conserva las seis filas v3 y añade exactamente
  once ámbitos F2. Los componentes están ordenados, los enteros tienen
  codificación canónica y `matrix_alternate_allele` no depende de ALT ni de
  `variant_id`.
- **Modo `exact`:** cuotas por mayores restos, permutación independiente por
  población y flujo separado de valor; `fixed` y cuota cero tienen el consumo
  nulo especificado.
- **Capabilities:** `available` es subconjunto de `catalogued`; F2.0 no
  anticipa opciones y F2.2.21/F2.3.9 habilitan cada autorización con su
  semántica completa.
- **Plataforma:** la evidencia procede de detección y mantiene separados WSL2,
  Ubuntu 24.04 nativo y revisiones manuales.
- **Productores compartidos de requisitos:** las nueve filas de §15.1 enlazan
  autoría y consumidores reales, incluidos individuos creados en F2.3 y
  manifiestos/verdad/VCF posteriores.
- **`GPS_CASE` y verdad:** cardinalidad multivaluada, unicidad, orden canónico,
  ALT↔casos y ausencia de complemento artificial en F2 son coherentes entre
  especificación, DTC, arquitectura, pasos y oráculo.
- **Perfiles compactos:** fuente pequeña, licencia, tamaño, SHA-256, receta
  64:12:8, transformación Harris 96→6→12, productor y oráculo independiente
  están definidos sin materializar el VCF de 1000 Genomes; el último intervalo
  es ahora inequívocamente `[0.05,1)`.
- **Evidencia:** cada paso tiene una orden concreta, destino JSON y regla de
  éxito/bloqueo; las pausas exigen aprobación explícita y los oráculos
  principales no usan al productor como autoridad.

## 10. Prerrequisito externo F1

F1 no se reevalúa ni se convierte en un defecto propio de F2 en este informe.
La revisión independiente ha emitido veredicto de **no aprobación**, por lo que
RAF2-027 continúa bloqueando el inicio aunque el plan F2 se corrija. F2.0.1
debe conservar la exigencia de revisión F1 superada, aprobación formal y commit
base publicado e identificable.

## 11. Condiciones para aprobar el plan F2

RAF2-041 a RAF2-044 están corregidos y los recuentos 128×5, seis pausas y
79+9+9 son coherentes. Para aprobar el plan queda repetir una revisión
adversarial sobre estos cambios y confirmar que no aparecen nuevos bloqueantes
o mayores. Por separado, F1 debe aprobarse y fijar su commit base antes de
ejecutar F2.0.1; RAF2-027 sigue siendo la única causa abierta conocida.

## 12. Validaciones realizadas

- Lectura completa del plan F2 y contraste con las fuentes de diseño vigentes.
- Recuento mecánico de 128 pasos y comprobación de exactamente cinco campos en
  cada uno.
- Comprobación de seis pausas y de ausencia de identificadores F2 inexistentes
  en referencias simples o extremos de rango de `Depende de`.
- Revisión manual de los inventarios 79 exclusivos, nueve compartidos y nueve
  parciales directos, incluidos sus productores y destinos.
- Inspección del catálogo RNG físico y de la CLI actual solo para distinguir
  prerrequisitos F1 de trabajo futuro F2; no se ejecutaron pruebas de producto
  ni se modificó código.
- Comprobación final de Markdown y espacios con `git diff --check` sobre este
  archivo nuevo.
