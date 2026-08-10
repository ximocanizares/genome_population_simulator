# Revisión adversarial del plan de implementación de F2 0.1

## 1. Estado y alcance

- **Documento revisado:**
  [plan-implementacion-f2-0.1.md](./plan-implementacion-f2-0.1.md).
- **Fecha de revisión:** 10 de agosto de 2026.
- **Estado del informe:** correcciones aplicadas; conserva pendientes que
  requieren decisión o consejo.
- **Resultado:** **no aprobar todavía el plan para ejecución**.
- **Alcance:** coherencia del plan con la especificación, la arquitectura, las
  decisiones técnicas, la hoja de ruta y la matriz de trazabilidad vigentes.

Esta revisión distingue los defectos del plan de los heredados de la matriz de
trazabilidad. Se han aplicado al plan y, cuando era inequívoco y no alteraba el
alcance funcional, a la matriz maestra. No se han modificado especificación,
arquitectura, decisiones técnicas ni hoja de ruta: los hallazgos que las
necesitan permanecen pendientes de revisión.

## 2. Veredicto ejecutivo

El plan tiene una estructura de ejecución sólida: 52 pasos pequeños, cinco
pausas obligatorias, fronteras claras entre verdad y VCF, pruebas negativas y
oráculos ajenos al productor. Sin embargo, no es todavía ejecutable ni puede
servir como matriz de cierre de F2.

De los once hallazgos identificados:

- **ocho están resueltos**;
- **tres permanecen pendientes de revisión**;
- los pendientes son RAF2-001, RAF2-004 y RAF2-007 y bloquean la aprobación.

Los principales bloqueos son las decisiones técnicas aún no incorporadas al
registro canónico, una trazabilidad que intenta cerrar requisitos VCF antes de
crear el VCF, cuatro requisitos F2 omitidos del inventario y varias
afirmaciones de cierre que exceden el alcance deliberadamente reducido de F2.

## 3. Hallazgos bloqueantes

### RAF2-001 — Las decisiones que definen F2 siguen pendientes en la fuente canónica

**Estado tras corrección: PENDIENTE DE REVISIÓN.**

**Evidencia.** El plan declara pendientes DTPD-101, DTPD-102 y DTPD-103 y las
convierte en A2, A3 y E1. También reconoce que la longitud exacta de los
*indels* se acordó pero no se trasladó al registro de decisiones. El registro
técnico todavía conserva DTPD-101 a DTPD-103 como pendientes y exige resolver
cada DTPD antes de planificar el hito que bloquea.

**Riesgo.** A1 autoriza una fase cuyo modelo de coordenadas, estado haplotípico,
verdad externa y VCF todavía podría cambiar. Los pasos consumidores parecen
detallados, pero su contrato no es estable. Una aprobación del plan podría
interpretarse erróneamente como aprobación de esas decisiones.

**Corrección necesaria.** Antes de aprobar F2:

1. registrar como decisiones confirmadas las elecciones ya acordadas sobre
   coordenadas, modelo haplotípico y `truth.json`;
2. registrar la decisión VCF ya acordada, con versión, normalización, campos
   `INFO`, compresión, índice, complemento y validadores;
3. registrar por separado el contrato de longitudes exactas de *indels* y la
   concreción de Poisson;
4. actualizar arquitectura y plan con los identificadores confirmados;
5. convertir A2, A3 y E1 en comprobaciones de prerrequisitos o prototipos de
   validación, no en lugares donde la implementación elige el diseño.

**Condición de cierre.** Ninguna DTPD que bloquee F2 figura pendiente y el plan
no contiene un consumidor anterior a su decisión canónica.

**Decisión o consejo necesario.** Deben confirmarse en las fuentes canónicas:

- DTPD-101: base 1 cerrada, base 0 semiabierta o valores etiquetados; el registro
  recomienda evaluar composición y bordes y mantener tipos de frontera;
- DTPD-102: secuencias completas, listas de variantes, segmentos persistentes o
  modelo híbrido, separando el estado interno del fichero de intercambio;
- DTPD-103: versión VCF, normalización, `INFO`, complemento, BGZF/índice y
  validadores;
- Poisson: dependencia Rust pura evaluada conforme a DTC-043 o implementación
  propia versionada y contrastada; se recomienda preferir la opción más simple
  que proporcione un oráculo independiente y consumo aleatorio estable;
- el contrato de longitudes exactas de *indels*, ya descrito en el plan, debe
  recibir un identificador canónico y aprobación explícita.

El plan deja esas elecciones en puertas anteriores a sus consumidores, pero no
puede convertirlas por sí mismo en decisiones confirmadas.

### RAF2-002 — F2.4 intenta cerrar dos requisitos que solo pueden demostrarse con VCF

**Estado tras corrección: RESUELTO.**

**Evidencia.** La matriz asigna REQ-VER-02 y REQ-VER-12 a F2.4, y los pasos
D1–D8 pretenden demostrarlos sobre el fichero estructurado de verdad. Sin
embargo:

- REQ-VER-02 exige que cada genotipo **del VCF** corresponda al individuo del
  archivo;
- REQ-VER-12 exige que los identificadores de muestra **del VCF** correspondan
  a los identificadores estables;
- la pausa F2.4 prohíbe depender del VCF y F2.5 es el primer hito que lo crea.

El problema procede de la matriz maestra y el plan lo reproduce.

**Riesgo.** F2.4 podría aprobar dos pruebas que no ejecutan el comportamiento
normativo. Posteriormente F2.5 no estaría obligado a aportar esa evidencia.

**Corrección necesaria.** Reasignar REQ-VER-02 y REQ-VER-12 a F2.5, o declarar
en F2.4 solo la evidencia parcial de identidad que alimentará las pruebas VCF.
Recalcular los recuentos de F2.4/F2.5, sus ACE y sus pasos de cierre en la matriz
maestra y en el plan.

**Condición de cierre.** Ambos TEST ejecutan un VCF real y comparan sus columnas
y genotipos con los individuos de la verdad canónica.

**Corrección aplicada.** REQ-VER-02 y REQ-VER-12 se han reasignado de F2.4 a
F2.5 en la matriz maestra y en el plan. E5 y E10 aportan la comprobación sobre
un VCF real; F2.4 queda con dos REQ de cierre y F2.5 con veinticinco.

### RAF2-003 — La lista de 58 REQ no es exhaustiva para las filas que involucran F2

**Estado tras corrección: RESUELTO.**

**Evidencia.** La validación del plan cuenta únicamente las filas cuyo destino
es exactamente un solo hito F2. La trazabilidad vigente contiene **62** filas
cuyo campo de hito incluye F2. Las cuatro omitidas son:

- REQ-IND-25, reasignado a F2.3/F6.1;
- REQ-IND-26, reasignado a F2.4/F3.5/F7.7/F8.5;
- REQ-MAN-10, reasignado a F2.3/F5.5/F9.7;
- REQ-MAN-83, reasignado a F2.5/F3.5.

El apartado 13.1 tampoco las reconoce como evidencia adyacente.

**Riesgo.** El creador real de individuos, la coherencia de su ID, la
procedencia efectiva y el primer conjunto real pueden implementarse sin
actualizar las pruebas que la trazabilidad dejó expresamente para F2.

**Corrección necesaria.** Mantener separados:

- los REQ cuyo cierre total corresponde a F2;
- los REQ compartidos para los que F2 debe aportar una evidencia parcial
  identificada y ejecutable.

Añadir las cuatro filas a la sección 13.1, enlazarlas desde C2/C8, D2/D6/D8 y
E8–E12 según corresponda, y hacer que la validación mecánica compruebe ambas
listas.

**Condición de cierre.** La unión de requisitos cerrados y compartidos coincide
exactamente con todas las filas de trazabilidad cuyo hito contiene F2.

**Corrección aplicada.** §13 conserva los 58 cierres exclusivos y §13.1 añade
una tabla separada para REQ-IND-25, REQ-IND-26, REQ-MAN-10 y REQ-MAN-83, con
evidencia concreta y pasos. Las comprobaciones exigen ahora que la unión sea
exactamente las 62 filas cuyo hito contiene F2.

### RAF2-004 — Algunos cierres de F2.2 exceden el alcance implementado

**Estado tras corrección: PENDIENTE DE REVISIÓN.**

**Evidencia.** F2 se limita a una referencia pequeña y a SNP e *indels*. Sin
embargo, F2.2 se propone cerrar requisitos cuyo texto completo incluye
capacidades excluidas:

- REQ-CONF-08 y REQ-CONF-21 incluyen variantes estructurales;
- REQ-IND-02 exige que un individuo pueda corresponder a cualquiera de las
  referencias disponibles;
- REQ-IND-03 exige conservar cuál de esas referencias describe cada muestra;
- F2 declara una sola referencia y aplaza variantes estructurales y referencias
  adicionales.

Además, REQ-IND-11 se incluye entre los 22 REQ, pero no existe un caso explícito
que demuestre que una frecuencia alélica poblacional se rechaza como valor
literal de un individuo.

**Riesgo.** El acta podría marcar requisitos completos después de demostrar
solo el subconjunto de SNP/*indels* o una ejecución monorreferencia.

**Corrección necesaria.** Revisar la matriz maestra requisito por requisito.
Cuando el REQ sea indivisible, asignar el cierre a la primera fase que cubra
todo su texto y conservar F2 como evidencia parcial. Si se pretende cerrar en
F2, ampliar explícitamente sus casos sin contradecir el alcance de la hoja de
ruta. Añadir el caso negativo de REQ-IND-11.

**Condición de cierre.** Ningún REQ se declara superado mediante una prueba que
solo cubra una clase o cardinalidad menor que la exigida por su texto.

**Corrección parcial aplicada.** B14 incluye ya el caso negativo literal de
REQ-IND-11 y prohíbe marcar REQ-CONF-08, REQ-CONF-21, REQ-IND-02 y REQ-IND-03
como cerrados mientras se resuelve este hallazgo.

**Decisión o consejo necesario.** Hay tres opciones: (a) mantener F2 como
evidencia parcial y trasladar el cierre total a la primera fase con variantes
estructurales o varias referencias; (b) dividir los REQ, lo que exige cambiar la
especificación; o (c) ampliar F2, contradiciendo su alcance deliberado actual.
Se recomienda (a), porque conserva los REQ y el alcance pequeño de F2. Requiere
aprobar y modificar la matriz maestra; no se ha hecho de forma unilateral.

### RAF2-005 — La procedencia normativa de varios pasos es semánticamente incorrecta

**Estado tras corrección: RESUELTO.**

**Evidencia.** Las cuatro líneas obligatorias existen mecánicamente, pero varios
enlaces `Procede de` no justifican el cambio del paso. Ejemplos:

- A2/A5 enlazan Poisson, sustituciones individuales o dosis para decidir tipos
  de coordenadas;
- A6 enlaza la regla de deleciones terminales para definir bases y `N`;
- B4/B6 enlazan identidad de referencia del individuo para definir Poisson;
- B7, B10 y B11 enlazan la semántica de frecuencia alélica poblacional;
- B8 enlaza las reglas de ancla de inserción y deleción terminal a un evento
  SNP.

**Riesgo.** Una revisión automatizada da por buena la trazabilidad por presencia
de IDs, aunque el requisito citado no pruebe el comportamiento. También quedan
REQ reales sin un paso que los demuestre directamente.

**Corrección necesaria.** Auditar semánticamente todos los `Procede de`. Cada
REQ debe aparecer en el paso que implementa o prueba su conducta concreta; los
requisitos meramente relacionados deben ir en contexto, no como fuente del
cambio. Repetir después el recuento REQ/ACE.

**Condición de cierre.** Para cada enlace puede explicarse una relación directa
entre el verbo normativo del REQ, el cambio y la comprobación del paso.

**Corrección aplicada.** Se auditaron los campos `Procede de` y se corrigieron
las asociaciones erróneas de A2, A3, A5, A6, B1–B11, C2/C8, D1–D8 y E5–E10.
Entre otros cambios, las reglas de inserciones y deleciones se citan solo en sus
constructores y espacios válidos; Poisson cita sus requisitos de media,
distribución y recuentos; y cada REQ compartido se etiqueta expresamente como
evidencia parcial.

## 4. Hallazgos mayores

### RAF2-006 — BGZF e índice se implementan con fuentes normativas equivocadas

**Estado tras corrección: RESUELTO.**

**Evidencia.** E9 atribuye compresión BGZF, índice y modo texto a REQ-REC-26,
REQ-MAN-27, REQ-MAN-57 y REQ-SAL-21. Esos requisitos tratan, respectivamente,
variantes eliminadas, fase/procedencia, relación VCF–muestras y estados. Los
contratos directos son REQ-SAL-12 y REQ-MAN-58, que la matriz sitúa en fases
posteriores.

**Riesgo.** E9 puede pasar sin demostrar los requisitos que realmente definen
sus artefactos, mientras estos se cerrarían otra vez en F5.5/F9.5.

**Corrección necesaria.** Citar REQ-SAL-12 y REQ-MAN-58 como fuentes directas y
registrar evidencia parcial si su cierre completo permanece en fases
posteriores. Mantener REQ-SAL-21 solo para la parte de estados y REQ-MAN-57 para
la relación con referencia e individuos.

**Corrección aplicada.** E9 cita REQ-SAL-12 y REQ-MAN-58 como fuentes directas
y conserva REQ-SAL-21 a REQ-SAL-23 para estados/publicación. REQ-MAN-57 se
trasladó a E8, REQ-REC-26 a E6 y REQ-MAN-27 a E7. La evidencia de las fuentes
asignadas a fases posteriores se declara parcial.

### RAF2-007 — No se identifica un caso F2 que necesite realmente complemento

**Estado tras corrección: PENDIENTE DE REVISIÓN.**

**Evidencia.** La hoja de ruta y E7 exigen demostrar un caso VCF+complemento,
pero el plan limita F2 a SNP, *indels* inequívocos, GT faseado y ploidía 1–10.
E7 no identifica qué información de uno de esos casos no puede representar la
versión VCF elegida. REQ-VER-24 prohíbe crear una entrada complementaria cuando
el caso es plenamente representable.

**Riesgo.** Para satisfacer la prueba se podría inventar una pérdida o crear un
complemento redundante, incumpliendo precisamente la regla funcional.

**Corrección necesaria.** DTPD-103 debe aportar un ejemplo literal de F2 cuya
pérdida sea real y verificable. Si no existe, corregir la evidencia de cierre de
F2 para demostrar la **ausencia** de complemento en los casos representables y
trasladar el primer caso positivo a F5.5.

**Corrección parcial aplicada.** E7 ya prohíbe inventar una entrada
complementaria y exige la prueba negativa para casos representables. Solo exige
un caso positivo si DTPD-103 demuestra antes una pérdida real.

**Decisión o consejo necesario.** Debe elegirse entre (a) documentar en
DTPD-103 un caso literal de SNP/*indel*/`N` de F2 que VCF no represente sin
pérdida o (b) corregir la evidencia de cierre de la hoja de ruta para demostrar
ausencia de complemento en F2 y trasladar el primer positivo a F5.5. Se
recomienda (b) si la evaluación de VCF no encuentra una pérdida real; es la
única opción compatible con REQ-VER-24 sin fabricar redundancia.

### RAF2-008 — La decisión aleatoria se registra después de diseñar sus consumidores

**Estado tras corrección: RESUELTO.**

**Evidencia.** B1 y B3 definen el TOML y `ParameterKey` para las longitudes
exactas antes de B4, aunque B4 es el paso que registra y aprueba ese contrato.
B4 solo exige estar resuelto antes de B6. Además, no existe una puerta
equivalente a A4/E2 para evaluar una posible dependencia Poisson conforme a
DTC-043.

**Riesgo.** El esquema público y el catálogo pueden fijar accidentalmente una
decisión aún no confirmada, o incorporarse un crate sin la evaluación exigida.

**Corrección necesaria.** Resolver el contrato funcional de longitud y la
concreción técnica de Poisson antes de B1. Añadir una puerta DTC-043, aunque el
resultado sea «implementación propia, sin dependencia nueva».

**Corrección aplicada.** La puerta de decisión es ahora B1; B2–B4, esquema,
parser y catálogo, dependen de ella. B1 evalúa expresamente cualquier
dependencia conforme a DTC-043, incluida la opción de no añadir ninguna.

### RAF2-009 — El orden entre espacio elegible y recuento Poisson es ambiguo

**Estado tras corrección: RESUELTO.**

**Evidencia.** B6 recibe el tamaño elegible para calcular la media, pero B7 es
quien construye el espacio elegible y depende de B6. El flujo normativo exige
calcular oportunidades válidas del ámbito y después sortear el recuento
correspondiente; la aplicación posterior debe trabajar sobre el estado válido
actualizado.

**Riesgo.** Puede usarse la longitud bruta en vez del número de oportunidades
válidas, especialmente con `N`, bordes de *indels* y regiones filtradas.

**Corrección necesaria.** Separar el enumerador de oportunidades de la selección
de posiciones. Hacer que el cálculo del tamaño elegible preceda al sorteo, y
dejar la actualización del espacio después de cada evento en el motor.

**Corrección aplicada.** B6 enumera y cuenta oportunidades válidas sin sortear;
B7 recibe esa cardinalidad y aplica Poisson. B12 mantiene la actualización del
estado vigente entre eventos y las pruebas verifican `media = tasa ×
oportunidades válidas`.

## 5. Correcciones menores

### RAF2-010 — La línea base de F1.4 está desactualizada

**Estado tras corrección: RESUELTO.**

El plan todavía afirma que F1.4 está en desarrollo. F1.4 ya fue aprobada; F1.5
continúa sin cerrar. Debe actualizarse el estado sin eliminar la puerta de
cierre total de F1.

La cabecera y §3.1 reflejan ahora F1.4 aprobada, F1.5 pendiente y la revisión
total de F1 todavía obligatoria.

### RAF2-011 — El intervalo de longitud usa una notación ambigua

**Estado tras corrección: RESUELTO.**

El texto `1..umbral estructural efectivo` puede interpretarse como inclusivo,
aunque el resto del plan rechaza longitudes `>= umbral`. Debe escribirse de
forma inequívoca: `1 ≤ longitud < umbral estructural efectivo`.

La notación se ha sustituido por esa desigualdad explícita.

## 6. Aspectos validados y reutilizables

No requieren rediseño los siguientes elementos:

- 52 pasos abiertos, cada uno con una sola línea `Depende de`, `Procede de`,
  `Cambia` y `Se comprueba con`;
- cinco pausas obligatorias y aprobación separada del cierre de fase;
- dirección correcta verdad canónica → verdad externa → VCF;
- prohibición de usar VCF o el propio productor como oráculo circular;
- conservación de `N`, rechazo de SNP sobre `N` y ausencia de `N` en
  inserciones aleatorias;
- separación de conteo Poisson y longitud exacta de *indel*, sin consumo de RNG
  para la longitud;
- aplicación ordenada de eventos y comprobación contra el estado vigente;
- publicación provisional, SHA-256, no-clobber y validadores independientes
  reutilizados desde F1;
- validación mecánica de formato sin errores de espacios ni casillas cerradas.

## 7. Trabajo restante antes de aprobar

1. Cerrar F1.5 y la revisión final de F1.
2. Resolver RAF2-001 en decisiones técnicas y arquitectura, y propagar los IDs
   confirmados al plan.
3. Resolver RAF2-004 en la matriz maestra sin ampliar silenciosamente F2.
4. Resolver RAF2-007 en DTPD-103 y, si no existe un caso positivo real de F2,
   corregir la evidencia de cierre de la hoja de ruta.
5. Repetir la revisión adversarial de esos tres cambios y aprobar después el
   plan antes de autorizar A1.

## 8. Preguntas reservadas para la revisión posterior

1. ¿Los REQ que incluyen variantes estructurales o varias referencias deben
   quedar como evidencia parcial en F2 o se desea ampliar el alcance inicial?
2. ¿Qué caso concreto de F2 pierde información en VCF 4.5 y exige complemento?
   Si ninguno la pierde, ¿se corrige la evidencia de cierre de la hoja de ruta?
3. ¿La decisión Poisson usará una dependencia evaluada o una implementación
   propia versionada?
Las cuatro filas compartidas ya no requieren pregunta: el plan registra
evidencia parcial en F2 y conserva su cierre total condicionado a las demás
fases de su asignación vigente.

## 9. Comprobaciones mecánicas ejecutadas

- `git diff --no-index --check` sobre plan e informe: sin diagnósticos.
- `git diff --check` sobre el cambio de trazabilidad: sin diagnósticos.
- Casillas de pasos: 52.
- Pausas obligatorias: 5.
- Campos obligatorios por paso: 4/4 en los 52 pasos.
- Identificadores de paso: A1–A10, B1–B14, C1–C8, D1–D8 y E1–E12, sin
  duplicados.
- REQ únicos en la matriz de cierre exclusivo del plan: 58, sin duplicados.
- REQ compartidos en §13.1: 4, sin duplicados.
- Filas de la matriz maestra cuyo hito contiene F2: 62.
- Diferencia entre la unión 58+4 del plan y las filas F2 de la matriz maestra:
  ninguna.
- REQ-VER-02 y REQ-VER-12: reasignados a F2.5 y comprobados sobre VCF real.
- Enlaces relativos consultados: todos sus destinos existen.
- Cercas Markdown: equilibradas en plan e informe.
