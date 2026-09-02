# Revisión adversarial del plan de implementación F2 0.2

> **Informe histórico.** La numeración del informe no corresponde al plan F2
> 0.2 vigente del 21 de agosto. Sus conclusiones sobre el plan anterior quedan
> sustituidas por DTC-079 y requieren reevaluación proporcional.

## 1. Cabecera, alcance y método

- **Documento revisado:**
  [plan-implementacion-f2-0.1.md](./plan-implementacion-f2-0.1.md), en el estado
  del árbol compartido del 15 de agosto de 2026.
- **Tipo de revisión:** adversarial, documental y de prerrequisitos; sus
  correcciones técnicas se incorporan al plan y documentos de diseño, sin
  modificar código.
- **Fuentes contrastadas:** especificación, arquitectura, decisiones técnicas,
  hoja de ruta, trazabilidad, registro de ampliaciones futuras, plan F1, punto
  de reanudación, revisión adversarial F2 0.1 y fronteras actuales del código,
  CLI y catálogo RNG.
- **Cambios incorporados desde la revisión 0.1:** DTC-058 a DTC-074,
  DTPD-121 resuelta mediante DTC-068, perfiles humanos compactos, DTC-069 para
  importación/exportación Marimo, resolución de RAF2-029/030 y correcciones
  RAF2 previas.
- **Escala de gravedad:** `bloqueante` impide autorizar o ejecutar con un
  contrato no ambiguo; `mayor` exige corregir el plan antes de aprobarlo, aunque
  no cambie por sí solo la biología; `menor` es una incoherencia editorial o de
  mantenimiento que no invalida el recorrido por separado.

Tras las correcciones, la revisión comprobó que el plan contiene **122 pasos**, que
cada paso tiene exactamente una línea `Depende de`, `Procede de`, `Cambia`,
`Se comprueba con` y `Evidencia ejecutable`, y que existen **seis pausas**. La
matriz contiene 79 cierres exclusivos, nueve requisitos compartidos y nueve
obligaciones parciales directas con cierre posterior. Estos recuentos son
correctos, pero no autorizan F2 mientras F1 siga sin aprobarse.

## 2. Veredicto ejecutivo

**No aprobar F2 para ejecución todavía.**

El plan corregido ha mejorado de forma sustancial desde la revisión 0.1: ya no depende de
materializar el VCF de 1000 Genomes, la receta compacta DTC-068 es ejecutable,
la importación Marimo tiene un límite técnico explícito, todas las decisiones
DTPD propias de F2 figuran resueltas y los comandos de evidencia son mucho más
concretos. Tras incorporar las correcciones técnicas, permanece **un único
bloqueante externo**, RAF2-027: F1 aún no está aprobada. No quedan hallazgos
mayores ni menores abiertos en esta revisión.

Las dos decisiones funcionales detectadas por esta revisión ya están resueltas:
DTC-070 fija la asignación `exact` sin sesgo espacial y DTC-065 fija
`GPS_CASE` multivaluado con relaciones ALT↔casos en `truth.jsonl`. Los bloqueos
restantes se corrigieron mediante DTC-071 a DTC-074, inventario completo,
pausas explícitas y subpasos jerárquicos. No debe iniciarse F2.0
hasta que F1 tenga revisión y aprobación explícitas y el árbol base esté
publicado e identificable.

## 3. Hallazgos bloqueantes

### RAF2-027 — F1 no es todavía una base aprobada ni inmutable

- **Estado:** bloqueado únicamente por la revisión, aprobación y commit base de F1.
- **Gravedad:** bloqueante de ejecución; no es una nueva decisión funcional de
  F2.
- **Evidencia:** la hoja de ruta y el punto de reanudación declaran F1.6
  implementada pero pendiente de revisión independiente, y el árbol contiene
  cambios F1 sin commit. El propio paso F2.0.1 exige «cierre y commit base de
  F1». La revisión independiente concurrente ha reproducido además un defecto
  grave en `src/results.rs::validate_effective_values`: acepta para `ind-1` el
  valor `variant_rate=0.01` de una regla base aunque una regla individual con
  `0.02` la domina, y acepta una fila de un objeto ajeno `zzz-foreign`. El
  supuesto oráculo independiente de `src/execution.rs` repite la misma lógica y
  no lo detecta.
- **Riesgo:** empezar Marimo, catálogos o esquemas F2 contra una frontera F1 aún
  corregible hace que la evidencia de compatibilidad se calcule sobre una base
  móvil. También impediría atribuir con precisión qué commit fue validado.
- **Corrección necesaria:** terminar la revisión independiente de F1.6, reparar
  el defecto, añadir regresiones con un oráculo que reconstruya contexto,
  identidad y regla ganadora sin reutilizar producción, retirar mientras tanto
  las afirmaciones 22/22 y 83/83, repetir sus oráculos, aprobar formalmente F1 y
  fijar el commit base antes de ejecutar F2.0.1. No basta con que las pruebas
  locales estén verdes.

### RAF2-028 — El contrato RNG de F2 es incompleto y el sorteo de ALT es circular

- **Estado:** resuelto el 15 de agosto de 2026 mediante DTC-071 y los pasos F2.2.8/28 y F2.3.3.
- **Gravedad original:** bloqueante técnico.
- **Evidencia:** F2.2.8 define
  `matrix/alternate_allele/<variant_id>`, pero F2.2.16 define `variant_id` a
  partir de referencia, secuencia, posición, `REF` y **`ALT`**. Por tanto, la
  subsemilla necesaria para sortear `ALT` depende del valor que todavía se debe
  sortear. Además, F2.2.8, F2.2.28 y F2.3.3 expresan ámbitos como rutas, pero no
  fijan el nombre canónico del ámbito, los componentes nombrados y ordenados,
  la siguiente versión del catálogo ni la versión de derivación. El catálogo
  actual v3 y `SubseedRequest` solo contienen los seis ámbitos de F1.
- **Riesgo:** un implementador tendría que inventar la identidad del flujo y la
  versión del catálogo. Dos implementaciones conformes al texto podrían
  producir matrices distintas, y el flujo de `ALT` no puede construirse
  literalmente.
- **Resolución aplicada:** `subseed-scopes-v4.tsv` es superconjunto estricto de
  v3, conserva derivación v1 y añade once ámbitos F2 con componentes/orden y
  consumidores exactos. `matrix_alternate_allele` usa la identidad previa
  `(scenario_id, reference_id, sequence_id, position_1based, ref)` y
  `variant_id` se crea después. Vectores y oráculos demuestran que v1–v3 no
  cambian y que matriz, indels, genotipo y fase están aislados.

### RAF2-029 — El modo `exact` fija cantidades por intervalo, pero no qué loci las reciben

- **Estado:** resuelto el 15 de agosto de 2026 mediante DTC-070.
- **Gravedad original:** bloqueante funcional.
- **Evidencia:** F2.2.13 reparte el número de loci entre filas mediante mayores
  restos, pero F2.2.24 vuelve a indicar que para cada
  `población × variant_id` se «elige fila y valor». No se establece cómo se
  asignan a loci concretos las cantidades exactas, ni si el orden canónico de
  variantes, una permutación aleatoria u otro mecanismo decide esa asociación.
- **Riesgo:** asignar consecutivamente las filas del espectro a variantes
  ordenadas por coordenada introduce una asociación artificial entre frecuencia
  y posición. Sortearlas sin una regla normativa rompe la reproducibilidad y el
  aislamiento por población. Los oráculos de F2.2.38 tampoco pueden calcular un
  resultado literal único.
- **Resolución aplicada:** primero calcula por mayores restos las cantidades
  exactas. Después, cada población permuta de forma reproducible e independiente
  todos los `variant_id` canónicos con una subsemilla exclusiva y aplica las
  filas siguiendo esa permutación; el valor del intervalo usa otro flujo por
  población×variante. `fixed` no consume RNG de valor y una cuota cero no
  consume RNG. Queda prohibida la asignación por coordenada. F2.2.13/24/38 y
  sus oráculos recogen cuotas, reordenación, ausencia de correlación espacial,
  aislamiento poblacional y separación de flujos. DTC-071 fija los nombres,
  componentes, orden y versión literales de ambos ámbitos.

### RAF2-030 — La agrupación VCF no define `GPS_CASE` cuando confluyen varios casos

- **Estado:** resuelto el 15 de agosto de 2026 mediante DTC-065/066 revisadas.
- **Gravedad original:** bloqueante funcional.
- **Evidencia:** DTC-064 y F2.5.5 exigen exactamente un registro por
  `CHROM+POS`, integrando todos los alelos finales compatibles. DTC-065 define
  `GPS_CASE` como el identificador local estable de «su caso» en singular. El
  plan no define qué ocurre si cambios o alelos de casos de verdad distintos se
  normalizan a la misma posición y deben compartir registro.
- **Riesgo:** seleccionar un único caso pierde trazabilidad; repetir la línea
  viola la decisión de un registro por posición; concatenar IDs sin declarar
  cardinalidad, orden y correspondencia ALT↔caso deja un VCF formal o
  semánticamente ambiguo. El oráculo independiente de F2.5.10 no puede fijar la
  cardinalidad caso↔registro.
- **Resolución aplicada:** `GPS_CASE` usa `Number=.`/`Type=String` y contiene
  los `case_id` originales únicos en orden canónico, sin caso compuesto.
  `truth.jsonl` conserva por cada ALT la lista ordenada completa de casos
  fuente. El corpus cubre dos ALT/dos casos, un ALT/dos casos, un caso/varios
  registros, reordenación, duplicados e ID inexistente. En F2 no aparecen
  `GPS_LOSS`, `GPS_COMPLEMENT` ni complemento porque todos sus casos son
  representables; F6.5 fijará la cardinalidad positiva multicaso únicamente
  sobre una pérdida real.

### RAF2-031 — El inventario de cierre omite siete obligaciones parciales asignadas a F2

- **Estado:** resuelto el 15 de agosto de 2026 en §15.2 y F2.5.12.
- **Gravedad original:** bloqueante de trazabilidad y aprobación.
- **Evidencia:** la sección 15 solo cuenta los 79 REQ exclusivos, nueve con F2
  en el campo `Hito` y dos parciales directos. Sin embargo, la matriz maestra
  asigna expresamente trabajo parcial F2 a otros siete REQ cuyo cierre final
  está en otra fase:

  | REQ | Trabajo que la trazabilidad exige ya en F2 | Cierre final |
  |---|---|---|
  | REQ-FMT-07 | lectura completa, normalización, consulta y oráculo VCF | F9.5 |
  | REQ-VER-23 | declaraciones exactas de `GPS_CASE/GPS_LOSS/GPS_COMPLEMENT` | F6.5 |
  | REQ-VER-24 | ausencia de complemento para casos representables | F6.5 |
  | REQ-VER-25 | mecanismo y ausencia condicional sin fabricar pérdida | F6.5 |
  | REQ-VER-26 | límite explícito de la reconstrucción positiva aplazada | F6.5 |
  | REQ-MAN-31 | enlace condicional de ruta, tamaño y SHA-256 del complemento | F6.5 |
  | REQ-MAN-58 | versión, normalización, modo, índice, rutas y validación VCF | F9.5 |

- **Riesgo:** F2 puede aprobarse con un recuento 79+9+2 aun faltando una de
  esas siete evidencias. A la inversa, un paso puede ejecutarlas sin producir
  expediente `partial` reutilizable para su cierre posterior.
- **Resolución aplicada:** §15.2 enumera nueve parciales directos —los dos ya
  presentes y estos siete— con pasos, destino y `status = "partial"`; F2.5.12
  exige ocho de ellos y F2.4.8 el restante. Ninguno se cuenta entre los 79
  cierres F2.

## 4. Hallazgos mayores

### RAF2-032 — Las dependencias formales no representan todas las pausas obligatorias

- **Estado:** resuelto el 15 de agosto de 2026 con PAUSA-F2.0 a PAUSA-F2.5.
- **Gravedad original:** mayor.
- **Evidencia:** la pausa F2.0 prohíbe iniciar F2.1, pero F2.1.1 solo depende de
  F1 y de la aprobación del plan; no declara la pausa F2.0. La pausa F2.1 dice
  específicamente «No iniciar F2.2.27 ni F2.3/F2.4», aunque F2.2.1 y F2.2.2 sí
  dependen de su aprobación y el orden general afirma cerrar F2.1 antes de
  F2.2. La política de bloqueo mezcla así una barrera total con una parcial.
- **Riesgo:** un ejecutor que siga `Depende de` puede adelantar pasos que el
  texto de pausa prohíbe, mientras una persona puede interpretar lo contrario.
- **Resolución aplicada:** las seis pausas tienen ID estable, el primer paso de
  cada hito depende de la anterior y PAUSA-F2.1 es una barrera total coherente;
  los pasos Marimo finales también nombran PAUSA-F2.0.

### RAF2-033 — Los algoritmos aleatorios, salvo Poisson, no están cerrados al nivel exigido

- **Estado:** resuelto el 15 de agosto de 2026 mediante DTC-072.
- **Gravedad original:** mayor.
- **Evidencia:** el plan menciona selección ponderada sin reemplazo, selección
  categórica de filas, `linear_uniform`, `log_uniform`, muestreo HWE/F y
  permutación de fase, pero no fija para todos ellos la transformación exacta de
  `next_u64`, los límites abiertos/cerrados, el número de consumos, el manejo de
  cero ni una versión de algoritmo. DTC-067 sí aporta ese nivel para Poisson y
  demuestra que el proyecto lo considera parte de la reproducibilidad.
- **Riesgo:** cambiar de una función auxiliar a otra o usar matemática de
  plataforma distinta altera matriz, genotipos o fase con la misma semilla. Un
  oráculo que elija su propia interpretación no prueba el contrato del
  producto.
- **Resolución aplicada:** DTC-072 fija conversión abierta, entero acotado,
  categórico, Fisher-Yates, selección ponderada sin reemplazo, uniformes
  lineal/log, HWE, F diploide y fase, con versión, fronteras, consumo cero,
  `libm`, vectores y oráculos ajenos. Solo `fixed` puede producir 0 o 1.

### RAF2-034 — `capabilities` anuncia autorizaciones antes de definir su disponibilidad efectiva

- **Estado:** resuelto el 15 de agosto de 2026 mediante DTC-073.
- **Gravedad original:** mayor.
- **Evidencia:** F2.0.8 obliga a publicar desde el Marimo base la unión
  `--allow-empty-snv-matrix` y `--allow-incompatible-genotypes`, pero sus
  semánticas y aceptación CLI no se implementan hasta F2.2.21 y F2.3.9. La CLI
  actual solo admite `validate`/`run`, diagnósticos y salida.
- **Riesgo:** el binario puede declarar una capacidad que aún rechaza, o F2.0
  tendría que anticipar opciones F2.2/F2.3 sin que el plan lo indique. Esto
  invalida la comparación aplicación↔binario usada para habilitar acciones.
- **Resolución aplicada:** `capabilities` distingue `catalogued_authorizations`
  de `available_authorizations`; Marimo solo usa la intersección. F2.2.21 y
  F2.3.9 activan cada flag únicamente junto con parser, validación y semántica.

### RAF2-035 — La evidencia de plataforma puede autoetiquetarse sin demostrar la plataforma real

- **Estado:** resuelto el 15 de agosto de 2026 mediante DTC-074 y F2.0.7/17.
- **Gravedad original:** mayor.
- **Evidencia:** F2.0.7 propone
  `verify-ui-platform --platform ubuntu-24.04` y `--platform wsl2`. El parámetro
  describe la plataforma esperada, pero el paso no exige comprobar
  `/etc/os-release`, kernel, arquitectura, marcador WSL, sistema de ficheros ni
  procedencia de CI. El entorno actual es WSL2 y no puede aportar por sí mismo
  evidencia nativa Ubuntu.
- **Riesgo:** ejecutar ambas órdenes en el mismo host produciría dos carpetas
  correctamente nombradas pero una prueba falsa. La revisión manual de lector
  de pantalla tiene el mismo problema si no conserva entorno y responsable.
- **Resolución aplicada:** el verificador registra OS, kernel, arquitectura,
  triple, WSL, montaje y procedencia CI, y rechaza expectativas discordantes.
  Ubuntu nativo enlaza un job real; WSL2 y las actas manuales quedan separados.

### RAF2-036 — Varias unidades siguen siendo demasiado grandes para una revisión aislada

- **Estado:** resuelto el 15 de agosto de 2026 con quince subpasos jerárquicos.
- **Gravedad original:** mayor.
- **Evidencia:** F2.0.13 agrupa inventario de grafo, apertura segura, defensa
  TOCTOU, copia privada, preservación opaca, límites y exportación atómica;
  F2.0.17 reúne navegador, accesibilidad automática y dos revisiones manuales;
  F2.2.14 combina herramienta, dos recursos, licencia y doble oráculo; F2.5.4
  reúne composición haplotípica, intervalo común, recorte y alineamiento a la
  izquierda.
- **Riesgo:** una sola casilla puede quedar verde aunque una subresponsabilidad
  carezca de prueba específica; las reparaciones posteriores serán difíciles de
  atribuir y revisar.
- **Resolución aplicada:** F2.0.13 y F2.0.17 tienen cuatro hijos cada uno,
  F2.2.14 otros cuatro y F2.5.4 tres. Los padres conservan ID, dependen de todos
  sus hijos y cada casilla tiene los cinco campos y evidencia propia.

### RAF2-037 — La matriz de pasos compartidos señala productores que no crean la evidencia indicada

- **Estado:** resuelto el 15 de agosto de 2026 en §15.1 y trazabilidad maestra.
- **Gravedad original:** mayor.
- **Evidencia:** §15.1 enlaza REQ-CONF-08 con F2.2.4, aunque ese paso solo fija
  metadatos de matriz; la autoría de cambios está en F2.2.5 y su consumo en
  pasos posteriores. REQ-CONF-21 apunta a F2.2.8, que solo versiona RNG y no
  demuestra la distribución espacial. REQ-IND-02 apunta a F2.2.27 y
  REQ-IND-03 a F2.2.18, antes de crear un individuo o conservar su referencia
  por muestra. F2.2.40 dice además «22 REQ previos más CONF-128/129 e
  IND-60/61», aunque esos cuatro ya forman parte de los 22 exclusivos.
- **Riesgo:** los expedientes pueden acreditar configuración, espacio o
  referencia individual sin haber ejecutado el consumidor real, y el texto
  puede interpretarse como 26 cierres en vez de 22.
- **Resolución aplicada:** CONF-08/21 enlazan autoría y consumidores reales;
  IND-02/03 pasan a creación/materialización F2.3 y su hito compuesto se corrige
  en trazabilidad. F2.2.40 declara 22 cierres más solo dos parciales distintos;
  F2.3.16 declara doce cierres y cinco compartidos.

## 5. Hallazgos menores

### RAF2-038 — El estado de RAF2-001 se contradice dentro del mismo plan

- **Estado:** resuelto el 15 de agosto de 2026.
- **Gravedad original:** menor.
- **Evidencia:** la cabecera y el cierre afirman que RAF2-001 queda resuelto,
  mientras la pausa F2.0 dice literalmente «RAF2-001 continúa abierto».
- **Resolución aplicada:** PAUSA-F2.0 declara RAF2-001 resuelto y conserva solo
  RAF2-027 como bloqueo externo por F1.

### RAF2-039 — DTPD-103 aparece a la vez resuelta y pendiente de completar

- **Estado:** resuelto el 15 de agosto de 2026.
- **Gravedad original:** menor.
- **Evidencia:** cabecera, arquitectura, DTC-066 y F2.5.1 la declaran resuelta,
  pero §4.2 todavía dice que «DTPD-103 debe completar» normalización, INFO,
  complemento y comandos.
- **Resolución aplicada:** §4.2 expresa consumo/verificación de DTC-062 a
  DTC-066 y mantiene DTPD-103 cerrada.

### RAF2-040 — Un nombre de artefacto vuelve a ser ambiguo en F2.5.6

- **Estado:** resuelto el 15 de agosto de 2026.
- **Gravedad original:** menor.
- **Evidencia:** la comprobación de F2.5.6 dice que «el fichero D conserva todo
  el historial», pero no existe un artefacto normativo llamado «fichero D».
- **Resolución aplicada:** F2.5.6 nombra los registros de
  modificaciones/procedencia de `truth.jsonl` que conservan el historial.

## 6. Decisiones que requieren al usuario

No queda una decisión funcional nueva de esta revisión pendiente del usuario.
RAF2-029 y RAF2-030 han sido resueltas mediante DTC-070 y DTC-065/066,
respectivamente.

No requieren decisión del usuario: el uso de `site_id` previo a `ALT`, la
versión del siguiente catálogo RNG, la detección real de plataforma, la
corrección de pasos/evidencias, la división jerárquica, la disponibilidad de
autorizaciones y la limpieza de estados documentales.

## 7. Revisión del estado de los hallazgos 0.1

| Hallazgo 0.1 | Estado actual |
|---|---|
| RAF2-004, RAF2-007 y RAF2-012 | Resueltos mediante evidencia parcial y destinos F5/F6/F7/F8 explícitos. |
| RAF2-013 | Sigue abierto como puerta externa: F1 aún no está aprobada; se actualiza como RAF2-027. |
| RAF2-014 y RAF2-015 | Resueltos en protocolo de autorizaciones, DTC-073 y regla DTC-049. |
| RAF2-016 | Resuelto: DTC-071/072 fijan ámbitos de indels y consumidores. |
| RAF2-017 | Resuelto por DTC-069 y F2.0.13.1–F2.0.13.4. |
| RAF2-018 | Resuelto por DTC-068 y la puerta compacta F2.2.14. |
| RAF2-019 a RAF2-021 | Resueltos por DTC-060, el entorno Playwright bloqueado y la separación F2.0.1/F2.0.2. |
| RAF2-022 a RAF2-025 | Resueltos: IDs jerárquicos, comandos, drenaje concurrente, pasos divididos e inventario parcial completo. |
| RAF2-026 | Resuelto: recuentos y estados documentales actualizados. |

## 8. Aspectos correctos y reutilizables

- DTC-068 evita una descarga pesada y conserva fuente, receta, licencia,
  tamaño, SHA-256, transformación y doble oráculo.
- El plan separa con claridad candidatos, ALT descartados y loci activos; no
  repone un candidato ausente en todas las poblaciones y conserva fijaciones.
- `poisson-count-v1` está limitado a recuentos y tiene frontera, matemática,
  consumo RNG y vectores mucho más precisos que el resto de muestreadores.
- La matriz suministrada exige `alt_af`, referencia identificada desde el
  manifiesto F1 y pareja TSV+JSON íntegra; VCF, MAF, DAF y AC/AN no se aceptan.
- Marimo continúa sin autoridad biológica, usa procesos sin *shell*, digest de
  revisión, destino no-clobber y directorio autocontenido limitado.
- `truth.jsonl` es anterior y autoritativo respecto al VCF, tiene reconstrucción
  independiente y evita usar el propio productor como único oráculo.
- El VCF único, VCF 4.5, BGZF/TBI/CSI, `QUAL=.`, `FILTER=PASS` y ausencia de
  campos inferidos están trazados a herramientas y oráculos independientes.
- Los 122 pasos y seis pausas son mecánicamente contables y todos los pasos
  incluyen los cinco campos exigidos por la convención de planificación.

## 9. Orden residual para autorizar F2

1. Terminar la reparación y revisión independiente de F1.6.
2. Aprobar F1 y fijar su commit base identificable, cerrando RAF2-027.
3. Repetir la validación mecánica y una revisión adversarial del plan corregido;
   no convertir la resolución documental de RAF2-028 a RAF2-040 en aprobación
   automática de F2.

## 10. Condiciones mínimas para cambiar el veredicto

F2 podrá aprobarse cuando F1 esté aprobada y su commit base figure en el
preflight, la validación mecánica del plan corregido continúe verde y una nueva
revisión adversarial no encuentre bloqueantes ni mayores sin destino aprobado.
Las demás condiciones originales de esta sección quedan satisfechas
documentalmente por RAF2-028 a RAF2-040.

## 11. Validaciones realizadas

- Recuento mecánico tras correcciones: 122 pasos, cero pasos sin alguno de los cinco campos y seis
  pausas.
- Contraste de trazabilidad: 88 filas con F2 en el campo de hito, más nueve
  obligaciones parciales cuyo cierre está en otra fase; dos ya figuraban en
  §15.2 y siete estaban omitidas.
- Inspección del estado actual: la CLI aún expone únicamente `validate`, `run`
  y `--version`; el catálogo físico v3 contiene los ámbitos de F1, mientras
  DTC-071 y F2.2.8 definen normativamente la ampliación v4 que F2 implementará
  sin reescribir v1–v3.
- No se ejecutaron pruebas de producto ni se modificó código: F1 está bajo una
  revisión independiente concurrente que ya ha reproducido el defecto descrito
  en RAF2-027, y el objetivo de este documento es evaluar el plan, no certificar
  esa implementación.
- El Markdown y los espacios finales de este archivo se comprobaron con
  `git diff --check --no-index`.
