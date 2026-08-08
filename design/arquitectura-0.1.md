# Arquitectura del simulador de genomas y poblaciones

- **Versión:** 0.1
- **Estado:** borrador para revisión
- **Fecha:** 7 de agosto de 2026
- **Especificación de origen:**
  [especificaciones-0.1.md](./especificaciones-0.1.md)
- **Hoja de ruta:** [hoja-de-ruta-0.1.md](./hoja-de-ruta-0.1.md)
- **Registro técnico:**
  [decisiones-tecnicas-0.1.md](./decisiones-tecnicas-0.1.md)
- **Lenguaje del núcleo:** Rust, según la especificación y la hoja de ruta 0.1

Este documento es derivado: organiza técnicamente los requisitos normativos sin
ampliarlos ni sustituirlos. La cobertura completa se consulta en
[trazabilidad-requisitos-0.1.md](./trazabilidad-requisitos-0.1.md).
Las capacidades post-v1 se remiten al
[registro de ampliaciones futuras](./registro-ampliaciones-futuras-0.1.md);
las elecciones DTP/DTPD siguen siendo planificación técnica, no funciones.

## 1. Finalidad y estado

Este documento define cómo se divide el sistema y qué contratos deben mantener
sus partes para satisfacer la especificación 0.1. No define archivos de código,
paquetes o bibliotecas de Rust, denominados *crates*, formatos de configuración,
algoritmos concretos ni un plan de
implementación.

La arquitectura cubre toda la versión 1 para evitar contratos incompatibles
entre fases. Sin embargo, su materialización seguirá la secuencia F1–F9 de la
hoja de ruta. Una capacidad futura no autoriza a implementarla antes de su fase.

La especificación y la hoja de ruta confirman Rust para el núcleo. Esta elección
no decide formatos, *crates*, versión del lenguaje, algoritmos, arquitectura ni
plataformas.

## 2. Contexto, alcance y límites

### 2.1. Responsabilidad del sistema

El sistema recibe una configuración, referencias y entradas opcionales; crea
una verdad biológica simulada; puede producir evidencia técnica de
secuenciación y alineamiento; valida la coherencia de los resultados; y conserva
la información necesaria para reconstruir qué ocurrió. Satisface principalmente
las secciones 2, 4, 5 y 6 de la especificación.

Las tres vías biológicas comparten el mismo modelo de referencia, individuo,
haplotipo, variante y verdad:

1. creación directa de individuos;
2. población inicial y evolución por generaciones;
3. familias y poblaciones segregantes definidas mediante un grafo de cruces.

Toda salida asociada a una observación se organiza mediante un evento explícito.
Una ejecución puede declarar cero eventos y, en ese caso, no se crea ninguno de
forma implícita. Todo evento tiene identidad, selección de individuos,
referencias aplicables, salidas y subsemilla. Solo los eventos que producen
FASTQ, BAM o CRAM constituyen su único experimento de secuenciación y requieren
tecnología, perfil, exactamente una tanda y una biblioteca por individuo; solo
los que producen alineamientos requieren además un mapeador. Todas las
bibliotecas se procesan juntas dentro de la tanda. Por tanto, un evento puede
materializar únicamente verdad genotípica o VCF sin parámetros ni entidades
técnicas de secuenciación.

### 2.2. Fronteras externas

Son fronteras del sistema:

- los ficheros estructurados de configuración general y poblacional;
- FASTA de referencia;
- FASTQ suministrados para contaminación;
- los estándares FASTA, FASTQ, SAM/BAM, CRAM y VCF;
- minimap2 como mapeador predeterminado y los mapeadores compatibles que se
  autoricen en el futuro;
- las herramientas externas necesarias para mapear, comprimir, ordenar,
  indexar o validar formatos, cuando la decisión técnica correspondiente las
  seleccione;
- el sistema de ficheros y los recursos de CPU, memoria y almacenamiento.

Cada integración externa se realiza mediante un contrato sustituible. El
dominio biológico no conoce comandos, rutas concretas, cabeceras físicas ni
detalles de procesos externos.

### 2.3. Fuera de alcance arquitectónico

No forman parte del sistema los *SNP callers*, la comparación de sus resultados,
la imputación, el faseado inferido ni la evaluación de herramientas externas.
Tampoco se incluyen las capacidades funcionales excluidas en la sección 3 de la
especificación, entre ellas señales de instrumento, BCL, paneles dirigidos,
Oxford Nanopore, PacBio CLR, aneuploidía individual, reproducción clonal,
conversión génica, calibración empírica automática y todo modelado específico de
elementos transponibles. V1 no define una frontera de entrada para FASTQ y
referencia de calibración ni para perfiles empíricos reutilizables.

Python no forma parte del núcleo. Solo podrá incorporarse más adelante como
herramienta auxiliar de análisis si se aprueba una decisión técnica que defina
su frontera, distribución y reproducibilidad.

## 3. Principios arquitectónicos

1. **La verdad biológica precede a la evidencia técnica.** Lecturas,
   alineamientos, contaminación, errores y calidades consultan la verdad, pero
   nunca la modifican.
2. **Cada dato lleva identidad y procedencia.** Referencia, sistema de
   coordenadas, individuo, haplotipo, evento y conjunto no se deducen de una
   ruta o de una posición implícita.
3. **Las referencias no se mezclan.** Un conjunto de resultados pertenece a
   una única referencia de mapeo y sus identificadores de variante son locales.
   En v1, una unidad reproductiva contiene únicamente individuos basados en la
   misma referencia estructural, aunque puedan observarse contra varias
   referencias de mapeo.
4. **La configuración se resuelve antes de simular.** Se aplican valores
   base y la dominancia parcial de REQ-CONF-116 a REQ-CONF-123, se validan
   máximas incomparables, duplicados, solapamientos y selectores no permitidos y se obtiene
   una configuración efectiva inmutable para la ejecución. Las sustituciones por
   identificador individual solo pueden resolverse contra individuos directos,
   fundadores o parentales expresos que ya existan en esta fase.
5. **Los estados inválidos no cruzan fronteras.** Los identificadores,
   coordenadas, ploidías, proporciones y relaciones se representan mediante
   valores validados o producen errores explícitos.
6. **La aleatoriedad tiene ámbitos independientes.** La evolución, los
   muestreos y cada evento explícito no comparten consumo accidental de una
   única secuencia aleatoria. La reproducibilidad se declara por niveles y
   alcanzar uno no presupone alcanzar el siguiente.
7. **Los estados y la publicación son jerárquicos.** Ejecución, evento, conjunto
   y artefacto distinguen como mínimo `incompleto`, `completo` y `fallido`. Cada
   artefacto se publica atómicamente solo después de validarse y una unidad
   superior solo está completa cuando lo están todas sus unidades obligatorias.
8. **Los recursos se observan antes de optimizar.** La arquitectura expone
   puntos de medida y admite representaciones compactas posteriores, sin fijar
   objetivos cuantitativos no acordados.

## 4. Dirección de dependencias

```text
Entradas y herramientas externas
             │
             ▼
Adaptadores de configuración, formatos y procesos
             │
             ▼
Orquestación de ejecución y de eventos
             │
             ▼
Servicios de dominio biológico y técnico
             │
             ▼
Modelo de dominio e invariantes compartidas

Persistencia y manifiestos implementan contratos definidos hacia el dominio;
el dominio no depende de sus formatos físicos.
```

La dirección permitida es hacia el modelo de dominio. Los componentes de
secuenciación pueden consultar referencias y verdad genotípica. Los componentes
biológicos no pueden depender de FASTQ, BAM, CRAM, VCF, mapeadores ni perfiles
técnicos. La validación puede consultar ambos lados, pero no corregirlos.

## 5. Modelo conceptual del dominio

### 5.1. Identidades

Cada identidad es estable y pertenece a una categoría. Identidades de
categorías distintas no son intercambiables aunque su texto coincida.

| Entidad | Identidad y ciclo de vida | Relaciones obligatorias |
|---|---|---|
| Ejecución | Se crea tras validar y resolver la configuración. No se reutiliza | Versión, entradas, configuración efectiva, semilla raíz, estados y manifiesto general |
| Referencia | Estable durante la ejecución | FASTA exacto, secuencias ordenadas, clasificación cromosoma/*scaffold*, origen y cambios si es derivada |
| Secuencia de referencia | Su identificador procede del primer elemento no vacío de la cabecera FASTA y es único dentro de su referencia | Descripción restante de cabecera, longitud, clasificación, ploidía efectiva y convenio de coordenadas |
| Segmento | Identifica una porción orientada y trazable | Estado actual, origen fundador o evento que lo creó, orientación y correspondencias disponibles |
| Variante o evento genómico | Identidad estable en el historial biológico | Clase, orden, marco de coordenadas previo, haplotipo afectado y estado final |
| Individuo | Único y no reutilizable en la ejecución | Referencia estructural, generación, población, parentales, haplotipos y observaciones |
| Copia haplotípica | Distinta aunque tenga la misma secuencia que otra | Individuo, secuencia, índice de copia, ploidía, segmentos y fundador |
| Población y subpoblación | Identidades no reutilizables | Generación, miembros, estado activa/extinguida y relaciones demográficas |
| Evento biológico | Identifica muestreo, cruce, mutación o cambio demográfico según su categoría | Entradas, orden, resultado y causas de omisión o fallo |
| Evento explícito | Unidad de observación independiente | Individuos, referencias, subsemilla y salidas; tecnología, perfil, tanda, bibliotecas y mapeador solo cuando procedan |
| Experimento de secuenciación | Condición de un evento que genera FASTQ, BAM o CRAM; es el único experimento del evento y no posee otra identidad | Tecnología, perfil, una tanda, una biblioteca por individuo y, si genera alineamientos, mapeador |
| Tanda | Única dentro de un evento de secuenciación | Evento, tecnología, perfil y todas las bibliotecas procesadas conjuntamente |
| Biblioteca | Única para un individuo dentro del evento | Evento, tanda, individuo, tecnología, perfil, objetivo WGS e índices técnicos |
| Grupo de lectura | Único en toda la ejecución y único para el individuo dentro del evento | Evento, individuo, biblioteca, tanda y archivos asociados |
| Instancia de lectura | Observación concreta con identidad técnica única dentro del evento; una pareja constituye una unidad y comparte el nombre exigido por el formato | Biblioteca actual, miembros cuando sea pareja, origen estable y mecanismo de copia o transferencia |
| Conjunto | Local a un evento y una referencia de mapeo | Mismas muestras que los demás conjuntos del evento y resultados de esa referencia |
| Caso de verdad | Local a un conjunto | Diferencia biológica simple o compuesta y uno o varios registros VCF |
| Registro VCF | Local a un conjunto y a un caso de verdad | Referencia, representación estándar y, si pierde detalle, una única entrada complementaria |
| Entrada complementaria | Única para un caso que pierde detalle en VCF | Identificador local estable, registros enlazados y detalle completo no representable |

La sintaxis y el mecanismo de generación de identificadores quedan pendientes.
El contrato exige unicidad, estabilidad, categoría y ausencia de reutilización.

### 5.2. Referencias, coordenadas y segmentos

Una referencia posee una secuencia ordenada de cromosomas o *scaffolds* y una
identidad verificable del FASTA exacto. Una referencia derivada enlaza una única
referencia de origen y la secuencia ordenada de cambios estructurales que la
creó. El catálogo de referencias es un grafo dirigido sin ciclos durante una
ejecución.

En un FASTA, el identificador de secuencia es el primer elemento no vacío que
sigue a `>` y termina en el primer espacio o carácter de separación. Debe ser
único dentro de la referencia. El resto de la cabecera se conserva como
descripción completa y no participa en la unicidad.

Toda secuencia suministrada o generada usa exclusivamente `A`, `C`, `G`, `T` y
`N`, sin distinguir mayúsculas de minúsculas durante la validación. `N` expresa
una base desconocida o no determinada. Su tratamiento operativo debe decidirse
funcionalmente para cada clase de operación; el sistema no puede resolverla al
azar, sustituirla ni atribuirle otro significado de forma implícita.

Toda coordenada contiene, como mínimo, referencia, secuencia, posición o
intervalo y convenio. La configuración utiliza base 1 e intervalos cerrados;
los formatos externos conservan su convenio oficial. La representación interna
se decidirá por separado. Ningún componente puede convertir coordenadas sin
usar el servicio común de coordenadas y registrar la conversión.

Un segmento conserva la secuencia actual o la forma de materializarla y su
procedencia. Puede apuntar a una copia fundadora, a una copia duplicada o a una
inserción sin origen fundador. Una deleción conserva el intervalo ausente. Las
correspondencias entre estados solo existen donde la transformación las hace
significativas; nunca se inventan para cubrir regiones sin homología.

### 5.3. Estado genotípico y verdad

El estado genotípico de un individuo agrupa, para cada secuencia, exactamente el
número de copias que indica su ploidía efectiva. Cada copia contiene una cadena
ordenada de segmentos y eventos suficiente para materializar su secuencia,
conservar fase y dosis y alcanzar un origen fundador o un evento de inserción.

El motor genómico aplica primero los eventos explícitos en su orden declarado.
Después recibe un recuento sorteado por clase y procesa, en orden fijo,
duplicaciones, deleciones, inversiones, translocaciones, repeticiones en tándem,
*indels* y SNP. Para cada evento consulta el espacio válido del estado actual y
lo actualiza tras aceptarlo. La ausencia de destinos suficientes es un error del
motor, no una advertencia ni una autorización para producir menos eventos.

El recuento sorteado de variantes aleatorias procede exclusivamente de Poisson.
Los eventos explícitos y el modo controlado aplicable entregan resultados
concretos y no implementan otra distribución. Las mutaciones de repeticiones en
tándem usan el modelo escalonado simple: cada evento añade o elimina una unidad
completa según el sesgo efectivo entre ambas direcciones.

La **verdad biológica** comprende referencias, haplotipos, genotipos,
parentescos, poblaciones, transmisiones, mutaciones y observaciones. Es la única
fuente autorizada para generar el fichero obligatorio de verdad y los VCF.

La **verdad técnica de lecturas** conserva el origen real de fragmentos,
secuencias original y observada, errores, duplicados, contaminaciones,
quimeras, descartes agregados y categoría sintética cuando corresponda. Esta
verdad se deriva de la biológica y del evento técnico, pero no la sustituye.
También distingue cobertura objetivo, WGS básica anterior a artefactos, final
del FASTQ y alineada observada por referencia de mapeo.

El alcance mínimo de los estados genotípicos completos depende de la vía:

- creación directa: todos los individuos creados;
- evolución: poblaciones inicial y final, y todos los individuos de las
  generaciones o muestreos solicitados; las demás generaciones conservan el
  historial resumido obligatorio;
- familias y poblaciones segregantes: todos los individuos de todos los eventos
  de cruce y poblaciones intermedias.

Un historial resumido nunca se presenta como estado genotípico completo. Un
evento explícito puede solicitar una materialización adicional de verdad o VCF
sin activar el recorrido de secuenciación.

### 5.4. Poblaciones y generaciones

Una generación activa contiene individuos que pueden actuar como progenitores.
La descendencia crea otra generación y sustituye a la anterior. Los individuos
históricos mantienen identidad y linaje, pero no reaparecen como miembros de una
generación posterior.

El estado de una subpoblación distingue N objetivo, intentos reproductivos,
éxitos, fallos por causa, N realizado y Ne realizado o «no evaluable» con causa.
Distingue también subpoblación activa, extinción
planificada y extinción por resultado. Una recolonización crea una identidad
nueva.

### 5.5. Eventos y resultados

Un evento de muestreo crea una observación no destructiva del estado biológico.
Todo evento es explícito y puede producir solo verdad, solo VCF o evidencia
técnica. Cuando genera FASTQ, BAM o CRAM constituye además un experimento de
secuenciación, sin otra identidad, con exactamente una tanda y una biblioteca
por individuo. Todas las bibliotecas se procesan juntas en la tanda. Cada evento
produce cero o más conjuntos, uno por referencia
aplicable. Los FASTQ pertenecen al evento y no a un conjunto concreto, porque
dentro del evento las mismas lecturas pueden mapearse contra varias referencias.

Ejecución, evento, conjunto y artefacto tienen como mínimo los estados
`incompleto`, `completo` y `fallido`. La omisión autorizada de un evento
dependiente se registra como resultado y diagnóstico, sin presentarlo como
completo. El estado de una unidad superior se calcula desde sus obligaciones y
unidades subordinadas; no se deduce de que exista algún fichero.

## 6. Componentes y responsabilidades

### 6.1. Catálogo de componentes

| Componente | Responsabilidad única | Entradas y salidas | Datos que posee o consulta |
|---|---|---|---|
| Admisión de configuración | Leer, combinar, validar y resolver configuración general y poblacional mediante la matriz de ejes y la dominancia normativa | Ficheros y entradas declaradas → configuración efectiva con procedencia o diagnósticos con reglas máximas | Posee modelos de entrada y configuración efectiva; consulta catálogos de capacidades |
| Identidad de ejecución | Crear y comprobar identidades tipadas y el contexto común | Configuración validada → identidad y contexto | Posee el registro de identidades; no posee datos biológicos |
| Aleatoriedad reproducible | Proporcionar flujos independientes por ámbito | Semilla y claves de ámbito → flujo aleatorio auditable | Posee la jerarquía de ámbitos y el registro de semillas; no interpreta decisiones biológicas |
| Catálogo de referencias | Validar, copiar, simular y derivar referencias | FASTA o configuración → referencias disponibles | Posee identidades, secuencias, linaje y recursos exactos de referencia |
| Coordenadas y transformaciones | Validar intervalos y convertir convenios sin ambigüedad | Coordenadas etiquetadas → coordenadas validadas o correspondencias | Posee reglas de conversión; consulta referencias y segmentos |
| Motor de verdad genómica | Crear y transformar haplotipos y eventos | Referencia, cambios y flujos aleatorios → estado faseado | Posee el estado genotípico, eventos, segmentos y correspondencias |
| Creación directa | Traducir la vía directa a operaciones del motor genómico | Configuración efectiva → individuos | Posee solo el estado de progreso; el motor posee la verdad resultante |
| Inicialización poblacional | Crear fundadores, frecuencias, genotipos, LD, HWE y F diploide inicial | Configuración poblacional → población inicial o error previo | Posee estadísticas solicitadas/realizadas y el estado «no aplicable» de F fuera de diploides; entrega verdad al repositorio biológico |
| Reproducción y recombinación | Validar ploidía y estructura canónica, formar bivalentes, gametos y descendencia trazable | Progenitores, mapas/tasas y flujo → descendientes o incompatibilidad reproductiva | Posee eventos de transmisión y comparaciones estructurales; consulta verdad parental |
| Dinámica evolutiva | Ejecutar las ocho etapas generacionales en su orden fijo | Población y fases → población final e historial de etapas aplicadas u omitidas | Posee estado generacional y resúmenes; delega cambios genómicos |
| Selección de viabilidad | Evaluar el SNP bialélico protegido y decidir supervivencia | Individuos elegibles, dosis y función → supervivientes o error | Posee evaluaciones y decisiones; consulta identidad única, dosis y valores genéticos |
| Grafo de familias | Validar y ejecutar el grafo declarativo de cruces | Grafo y parentales → poblaciones intermedias y finales | Posee dependencias y resultados del grafo; consulta selección y reproducción |
| Selección declarativa | Evaluar el catálogo común sobre verdad biológica | Reglas y elegibles → seleccionados y auditoría | No posee genotipos; posee evaluación y resultado de selección |
| Muestreo y resúmenes | Crear observaciones no destructivas y ámbitos estadísticos | Generación y regla → muestra, resumen poblacional y muestral | Posee observaciones y resúmenes; consulta verdad e historial |
| Orquestación de eventos | Validar salidas y parámetros condicionales y aislar cada observación | Observación/configuración → plan de evento y conjuntos | Posee estado del evento y relaciones entre artefactos; solo crea un plan técnico si procede |
| Simulación WGS | Generar fragmentos e instancias de lectura desde haplotipos | Verdad, perfil y objetivo → lecturas con identidad y origen | Posee lecturas previas al demultiplexado y verdad técnica asociada |
| Artefactos y demultiplexado | Convertir conjuntamente la contaminación y aplicar errores, duplicados, adaptadores y reasignaciones | Lecturas y configuración → FASTQ finales y procedencia | Posee transformaciones técnicas e identidades de instancia; consulta fuentes compatibles |
| Alineamiento | Ejecutar el mapeador sobre el FASTQ final | FASTQ y referencia → alineamientos observados | Posee invocación, salida, diagnóstico y metadatos del proceso |
| Alineamiento directo | Crear registros sintéticos coherentes sin FASTQ | Verdad y perfil directo → alineamientos sintéticos | Posee categorías y registros sintéticos; no llama a un mapeador |
| Materialización de verdad | Aplicar el único contrato obligatorio de trazabilidad y escribir verdad, VCF y complemento | Estado faseado, historial y referencia de conjunto → artefactos de verdad y resúmenes | No altera el dominio; materializa el detalle exigido por vía sin ofrecer niveles ni degradarlo por escala |
| Ensamblado de resultados | Organizar recursos comunes, eventos y conjuntos | Artefactos validados → estructura de resultados | Posee catálogo lógico de artefactos y estados de completitud |
| Validación interna | Comparar configuración, verdad y evidencia | Objetivos y resultados → informe, aviso o error | Posee evidencias de comprobación; no corrige resultados |
| Manifiestos y exportación | Registrar procedencia y crear paquetes autocontenidos | Catálogo de artefactos y ejecución → manifiestos/paquete | Enumera individuos materializados y de familias; resume generaciones evolutivas no materializadas y posee enlaces e integridad |
| Estimación y observación de recursos | Estimar y medir tiempo, memoria y almacenamiento | Configuración, trazabilidad y salidas activadas → estimaciones y medidas | Posee mediciones; no decide modelos biológicos ni reduce resultados para ajustarlos al presupuesto |

### 6.2. Contratos comunes de los componentes

Cada componente ofrece una interfaz pequeña expresada mediante conceptos del
dominio. No recibe estructuras físicas de otro componente salvo en las
fronteras de E/S. Toda operación devuelve un resultado válido o un error con:

- categoría y código estable de diagnóstico;
- entidad, requisito y ámbito afectados;
- valor recibido y condición incumplida cuando sea seguro registrarlos;
- causa original y cadena de contexto;
- alcance de la detención: configuración, ejecución, evento, conjunto o rama;
- artefactos provisionales que no deben publicarse como completos.

Los componentes no registran dos veces la misma regla. Admisión decide
precedencias de configuración; coordenadas decide conversiones; identidad
decide unicidad; el motor genómico decide validez estructural; el orquestador
generacional decide el orden; el ensamblador decide completitud.

### 6.3. Dependencias permitidas y prohibidas

- Creación directa, inicialización, evolución y familias pueden depender del
  motor genómico, aleatoriedad, identidad y selección declarativa.
- Reproducción es compartida por evolución y familias; no se duplica por vía.
- Reproducción compara por secuencia la ploidía y el estado estructural canónico
  de todas las copias antes de formar bivalentes. La equivalencia no depende de
  que las historias compartan identificadores, pero no elimina ninguna traza.
- El motor genómico debe consultar la identidad del locus seleccionado para
  rechazar eventos explícitos solapantes y excluirlo de los espacios aleatorios;
  no puede moverlo, copiarlo ni eliminarlo.
- Reproducción y dinámica evolutiva deben rechazar antes de formar descendencia
  cualquier unidad reproductiva que mezcle referencias estructurales. Esta
  validación no restringe las referencias de mapeo de un evento explícito.
- Simulación WGS, artefactos, alineamiento y modo directo pueden consultar
  verdad e identidades, pero no llamar a componentes que cambien genotipos.
- Materialización VCF depende del estado faseado y de la referencia del
  conjunto; no depende de lecturas ni alineamientos.
- Validación depende de vistas de solo lectura de configuración, verdad y
  evidencia. No puede normalizar o reparar archivos después del hecho.
- Manifiestos reciben hechos ya registrados por los responsables. No
  reconstruyen decisiones a partir de nombres de fichero o mensajes de log.
- Ningún componente de dominio depende de crates o herramientas concretas. Los
  adaptadores futuros implementarán los contratos aprobados.

### 6.4. Comprobación aislada

Cada componente admite dobles de prueba para sus dependencias externas. Las
pruebas del dominio usan referencias y poblaciones pequeñas construidas de forma
controlada. Los adaptadores de formato se prueban con casos válidos e inválidos
del estándar elegido. Los procesos externos se prueban con respuestas de éxito,
salida inválida, código de error, interrupción y versión incompatible.

| Componente | Errores propios principales | Comprobación aislada mínima |
|---|---|---|
| Admisión de configuración | Sintaxis, campo desconocido, precedencia conflictiva, distribución o modelo no admitidos, referencia inexistente y solicitud de calibración empírica | Resolver una configuración mínima y rechazar conflictos, calibración o perfiles empíricos y demás alternativas fuera de v1 con diagnósticos que señalen parámetro, valor y ámbitos |
| Identidad de ejecución | Identidad duplicada, categoría incorrecta o reutilización | Crear identidades por categoría y demostrar colisiones y referencias cruzadas inválidas |
| Aleatoriedad reproducible | Ámbito desconocido, clave duplicada o estado no regenerable | Repetir decisiones por ámbito y demostrar que añadir un muestreo no cambia el flujo evolutivo |
| Catálogo de referencias | FASTA inválido, identificador de primer elemento repetido, símbolo distinto de A/C/G/T/N, linaje cíclico o derivación imposible | Separar identificador y descripción, validar el alfabeto, conservar `N` con coordenada y longitud y reconstruir una referencia pequeña desde su origen y cambios |
| Coordenadas y transformaciones | Convenio ausente, intervalo fuera de límites o correspondencia inexistente | Probar extremos, inserciones, inversiones y conversiones de ida y vuelta donde sean definidas |
| Motor de verdad genómica | Evento solapado prohibido, objetivo eliminado, destinos válidos insuficientes, punto de ruptura o ploidía inválidos | Materializar haplotipos conocidos; rechazar SNP o motivo sobre `N`, conservarla dentro de SV e inserciones explícitas y agotar de forma controlada el espacio válido |
| Creación directa | Individuo incompleto o cambio incompatible | Crear individuos pequeños con cambios concretos y aleatorios controlados |
| Inicialización poblacional | Frecuencias, recuentos, F fuera de [0,1], F en ploidía distinta de 2, LD o restricciones incompatibles | Enumerar HWE haploide, diploide y poliploide bialélico/multialélico; comprobar F diploide 0, intermedio y 1; rechazar F en P=1 y P>2; probar distribución explícita poliploide y F realizado «no aplicable» |
| Reproducción y recombinación | Referencias distintas, ploidía impar o desigual, estado estructural no idéntico, mapa incompleto o bivalente imposible | Probar diploides y tetraploides idénticos; SNP, *indels* pequeños y repeticiones distintos; cada clase SV heterocigota; SV fijada con historias distintas; *scaffold* intacto; SV de novo seguida de otro intento |
| Dinámica evolutiva | Orden alternativo solicitado, población activa vacía no autorizada, mezcla de referencias estructurales tras migración o evento demográfico imposible | Ejecutar las ocho etapas, omitir una opcional en su posición, rechazar otro orden, ordenar varios eventos demográficos y comprobar el historial efectivo |
| Selección de viabilidad | Locus ausente, múltiple, ambiguo o situado sobre `N`; dosis sin valor; aptitud negativa, no finita o máximo cero | Probar varias ploidías, rechazo de `N`, fallos antes de supervivencia y recombinación que conserva la identidad única |
| Grafo de familias | Ciclo real, fuente futura, cardinalidad o salida imposible | Validar grafos válidos ramificados y ciclos transitivos mínimos |
| Selección declarativa | Regla desconocida, conjunto ambiguo o elegibles insuficientes | Evaluar Y/O/NO y cantidades exactas/proporcionales sobre un catálogo de verdad fijo |
| Muestreo y resúmenes | Duplicado, estratos solapados, tamaño excesivo o ámbito estadístico confundido | Comparar resumen poblacional y muestral y demostrar observación no destructiva |
| Orquestación de eventos | Evento implícito, parámetro técnico indebido, tecnología múltiple, perfil o mapeador incompatible y referencia ajena | Construir eventos solo de verdad y eventos técnicos independientes, y comprobar requisitos condicionales y aislamiento |
| Simulación WGS | Objetivo ausente/doble, perfil incompatible, identidad repetida o fragmento imposible | Comprobar recuentos, cobertura haploide, aporte por copia y muestreo de `N` con Phred 0 antes de sesgos |
| Artefactos y demultiplexado | Conversión conjunta inviable, exceso de contaminación, colisión de identidad o par huérfano | Resolver una asignación entera conocida; reconstruir origen y diferencia verdadera/observada de cada lectura final, y comprobar agregados de descartes |
| Alineamiento | Herramienta ausente, versión incompatible, código no cero o salida inválida | Sustituir el proceso externo y comprobar comando, contexto y propagación de cada fallo |
| Alineamiento directo | Proporciones, MAPQ, banderas, CIGAR o pares incoherentes | Generar un caso por categoría y validarlo contra la referencia y verdad de origen |
| Materialización de verdad | VCF no conforme, SNP con `REF=N`, alelo no reconstruible, cardinalidad ambigua, enlace roto o pérdida de fase | Validar casos sin/con complemento, BND, dos referencias y cambios con `N`; reconstruir cada caso enlazado sin inventar una base y compararlo con la verdad canónica |
| Ensamblado de resultados | Artefacto activado ausente o mezcla de referencias | Intentar cerrar conjuntos completos e incompletos y verificar estados resultantes |
| Validación interna | Evidencia insuficiente, invariante rota o desviación fuera de política | Inyectar resultados conocidos y comprobar conforme, aviso y error sin corrección automática |
| Manifiestos y exportación | Enlace roto, identidad incoherente o integridad fallida | Comprobar la enumeración de una familia completa y de los individuos evolutivos materializados, junto con los resúmenes no materializados; recorrer enlaces y verificar el paquete |
| Estimación y observación de recursos | Dato no medible, estimación inválida o presupuesto excedido | Comparar estimación y medida y rechazar antes de simular un caso insuficiente sin degradar trazabilidad ni salidas |

## 7. Flujos de datos y ejecución

Esta sección describe la **composición técnica** de los componentes. No redefine
el comportamiento funcional: las condiciones, órdenes, cardinalidades y
errores observables pertenecen a los REQ de la especificación, y los ACE
indican cómo comprobarlos.

### 7.1. Admisión y preparación

Los adaptadores de entrada entregan representaciones sin normalizar al
componente de admisión. Este coordina validadores de sintaxis, identidad,
referencias y reglas; el resolvedor de configuración produce una vista efectiva
inmutable y conserva la procedencia de su resolución. El catálogo de entradas
copiadas y el plan de ejecución solo se crean tras esa frontera.

El estimador de recursos consulta el plan y las salidas activadas sin modificar
la configuración. La identidad de ejecución y el gestor de aleatoriedad se
inyectan después en los motores de dominio. El contrato funcional de esta
frontera está en REQ-CONF-01 a REQ-CONF-123, REQ-FMT-01 a REQ-FMT-07,
REQ-NOF-01 a REQ-NOF-08 y REQ-VAL-01 a REQ-VAL-51; se comprueba
principalmente mediante ACE-01, ACE-23, ACE-35, ACE-36, ACE-64, ACE-65,
ACE-105 y ACE-111.

### 7.2. Referencias y estado biológico

El catálogo de referencias es responsable de identidad, metadatos y acceso a
secuencias. El servicio de coordenadas y transformaciones es la única frontera
que traduce posiciones entre una referencia estructural, sus estados derivados
y las referencias de mapeo. Los motores biológicos reciben estas interfaces y
producen cambios de estado, nunca rutas de fichero ni registros de un formato
externo.

El repositorio de verdad conserva estados faseados, genealogía, procedencia y
capas de cambio. Creación directa, inicialización poblacional, reproducción,
recombinación, dinámica evolutiva y mutación acceden a él mediante operaciones
de dominio explícitas. El orquestador decide qué componente invocar; los
componentes no se llaman entre sí para alterar el orden funcional.

Las fuentes normativas son REQ-REF-01 a REQ-REF-23, REQ-IND-01 a
REQ-IND-59, REQ-DIR-01 a REQ-DIR-05, REQ-EVO-01 a REQ-EVO-73,
REQ-POB-01 a REQ-POB-03, REQ-CRU-01 a REQ-CRU-45, REQ-REP-01 a
REQ-REP-26, REQ-MIG-01 a REQ-MIG-09, REQ-SEL-01 a REQ-SEL-17,
REQ-REC-01 a REQ-REC-29 y REQ-MUT-01 a REQ-MUT-15. Sus recorridos se
aceptan con ACE-02 a ACE-14, ACE-24 a ACE-32, ACE-75 a ACE-91, ACE-95 a
ACE-98, ACE-103 a ACE-107 y ACE-109 a ACE-111 según la asignación individual
de la matriz de trazabilidad.

### 7.3. Observación técnica

La orquestación de eventos recibe una vista de solo lectura de la verdad y crea
un contexto técnico aislado. Simulación WGS, artefactos, demultiplexado y
alineamiento trabajan sobre ese contexto; los adaptadores de procesos externos
solo traducen invocaciones y resultados, y no adquieren autoridad para cambiar
la verdad.

Los FASTQ y demás evidencias nacen en el contexto del evento. Las vistas de
verdad interoperable nacen del repositorio canónico y del contexto de una
referencia, no de inferencias sobre las lecturas. Ambas ramas se reúnen
únicamente al ensamblar y validar el conjunto de resultados.

Esta frontera implementa REQ-TEC-01 a REQ-TEC-10, REQ-RG-01 a REQ-RG-07,
REQ-DAT-01 a REQ-DAT-73, REQ-WGS-01 a REQ-WGS-14, REQ-CONT-01 a
REQ-CONT-16 y REQ-VER-01 a REQ-VER-47. La comprobación funcional se
distribuye entre ACE-15 a ACE-21, ACE-30, ACE-31, ACE-37 a ACE-63,
ACE-66 y ACE-68 a ACE-74.

### 7.4. Cierre y publicación

Cada productor entrega un artefacto provisional acompañado de identidad,
referencia, procedencia y métricas técnicas. El ensamblador forma conjuntos
coherentes; los validadores consumen vistas de solo lectura y emiten resultados
estructurados. El escritor del manifiesto recibe hechos ya producidos y no
reconstruye decisiones desde nombres o mensajes.

La publicación es la última frontera. Solo ella mueve un conjunto desde el
espacio provisional al espacio publicado, mediante una operación recuperable y
atómica desde la perspectiva del consumidor. Los estados y efectos funcionales
de cierre están en REQ-VAL-01 a REQ-VAL-51, REQ-MAN-01 a REQ-MAN-88 y
REQ-SAL-01 a REQ-SAL-26, con ACE-22 a ACE-24, ACE-33, ACE-34, ACE-36,
ACE-46 a ACE-67, ACE-78, ACE-86, ACE-92 a ACE-108.

## 8. Invariantes de frontera

Los invariantes de esta sección protegen la dirección de dependencias. Su
contenido funcional procede de los REQ citados; aquí solo se asignan
responsables y puntos de control.

| Invariante arquitectónico | Responsable y frontera de control | Fuente funcional |
|---|---|---|
| Toda entidad cruza fronteras con identidad tipada, no inferida de un nombre físico | Identidad de ejecución, catálogos y ensamblador de resultados | REQ-CONF-01 a REQ-CONF-14, REQ-TEC-01 a REQ-TEC-10 y REQ-MAN-01 a REQ-MAN-88 |
| Toda coordenada viaja asociada a su referencia y convenio | Catálogo de referencias y servicio de coordenadas | REQ-REF-01 a REQ-REF-23 y REQ-VER-01 a REQ-VER-47 |
| Solo los motores biológicos escriben en el repositorio de verdad | Interfaces de dominio y repositorio de verdad | REQ-IND-01 a REQ-IND-59, REQ-EVO-01 a REQ-EVO-73 y REQ-CRU-01 a REQ-CRU-45 |
| Los componentes técnicos observan instantáneas y no mutan genotipo, genealogía ni trayectoria | Orquestador de eventos y vistas de solo lectura | REQ-TEC-01 a REQ-TEC-10 y REQ-DAT-01 a REQ-DAT-73 |
| La configuración efectiva es inmutable durante cada ámbito de uso y mantiene procedencia | Resolvedor de configuración | REQ-CONF-15 a REQ-CONF-123 y REQ-MAN-01 a REQ-MAN-88 |
| Las decisiones aleatorias consumen flujos identificados por ámbito | Gestor de aleatoriedad y motores consumidores | REQ-CONF-03 a REQ-CONF-08, REQ-MUT-01 a REQ-MUT-15, REQ-EVO-01 a REQ-EVO-73, REQ-REP-01 a REQ-REP-26, REQ-REC-01 a REQ-REC-29 y REQ-DAT-01 a REQ-DAT-73 |
| Un artefacto conserva estado explícito desde su producción provisional hasta su publicación | Productores, validadores, ensamblador y publicador | REQ-VAL-01 a REQ-VAL-51, REQ-MAN-01 a REQ-MAN-88 y REQ-SAL-01 a REQ-SAL-26 |
| Verdad canónica, vistas interoperables y evidencia técnica permanecen como clases distintas | Repositorio de verdad, materializador VCF y productores técnicos | REQ-VER-01 a REQ-VER-47, REQ-DAT-01 a REQ-DAT-73 y REQ-SAL-01 a REQ-SAL-26 |

Ningún adaptador de formato, proceso externo, manifiesto o escritor de salida
puede saltarse estas fronteras. Los casos observables son los ACE asignados a
cada REQ en
[trazabilidad-requisitos-0.1.md](trazabilidad-requisitos-0.1.md).

## 9. Configuración, persistencia y resultados

### 9.1. Clases de información

La arquitectura separa cuatro almacenes lógicos para impedir que una
representación derive autoridad de otra:

| Clase | Propietario | Mutabilidad y uso |
|---|---|---|
| Configuración declarada | Adaptador de entrada | Se conserva para auditoría; no es consultada directamente por los motores |
| Configuración efectiva | Resolvedor de configuración | Vista inmutable que consumen orquestadores y motores |
| Estado interno y verdad | Repositorios de dominio | Evoluciona solo mediante operaciones biológicas y conserva lo necesario para reconstrucción |
| Resultados publicados | Ensamblador, validadores y publicador | Se hacen visibles al exterior después del cierre del conjunto |

Los contenidos exigidos a estas clases están en REQ-CONF-01 a REQ-CONF-123,
REQ-VER-01 a REQ-VER-47, REQ-MAN-01 a REQ-MAN-88 y REQ-SAL-01 a
REQ-SAL-26. Esta separación no decide sus esquemas físicos.

### 9.2. Espacios lógicos

La raíz de ejecución, los recursos compartidos, los eventos y las referencias
de mapeo son espacios lógicos distintos. Las identidades tipadas enlazan esos
espacios; las rutas, extensiones o nombres físicos son responsabilidad de los
adaptadores de persistencia. Así se puede cambiar el empaquetado sin alterar el
dominio ni usar una convención de directorios como clave biológica.

### 9.3. Escritura y recuperación

Los productores pueden escribir por flujo o por fragmentos. Cada artefacto pasa
por estados internos explícitos —provisional, validado, fallido o publicado— y
el manifiesto se construye a partir de esas transiciones. El mecanismo concreto
de temporales, reanudación y limpieza sigue siendo una decisión técnica de
planificación; los efectos observables de completitud pertenecen a REQ-VAL-01 a
REQ-VAL-51, REQ-MAN-01 a REQ-MAN-88 y REQ-SAL-01 a REQ-SAL-26.

## 10. Errores y observabilidad

### 10.1. Estrategia de errores

Los componentes devuelven errores tipados de su propio vocabulario. El
orquestador añade ejecución, etapa y entidad; los adaptadores añaden ruta,
herramienta o posición cuando corresponda. Una capa puede traducir un error al
vocabulario de la capa superior, pero conserva la causa y no repara ni
reclasifica el resultado por su cuenta.

| Origen técnico | Contexto que debe aportar la arquitectura | Fuente de la consecuencia funcional |
|---|---|---|
| Admisión o entrada | Campo, valor, ubicación y reglas implicadas | REQ-CONF-01 a REQ-CONF-123, REQ-FMT-01 a REQ-FMT-07, REQ-REF-01 a REQ-REF-23 y REQ-VAL-01 a REQ-VAL-51 |
| Motor de dominio | Entidad, generación o evento y estado previo identificable | REQ-IND-01 a REQ-IND-59, REQ-EVO-01 a REQ-EVO-73, REQ-CRU-01 a REQ-CRU-45, REQ-REP-01 a REQ-REP-26, REQ-MIG-01 a REQ-MIG-09, REQ-SEL-01 a REQ-SEL-17, REQ-REC-01 a REQ-REC-29 y REQ-MUT-01 a REQ-MUT-15 |
| Productor técnico o proceso externo | Evento, herramienta, versión, invocación y artefacto afectado | REQ-TEC-01 a REQ-TEC-10, REQ-DAT-01 a REQ-DAT-73, REQ-WGS-01 a REQ-WGS-14, REQ-CONT-01 a REQ-CONT-16 y REQ-VAL-01 a REQ-VAL-51 |
| Materialización o publicación | Conjunto, referencia, artefactos y enlaces implicados | REQ-VER-01 a REQ-VER-47, REQ-MAN-01 a REQ-MAN-88, REQ-SAL-01 a REQ-SAL-26 y REQ-VAL-01 a REQ-VAL-51 |
| Defecto de invariante interno | Componente, operación y contexto reproducible | REQ-VAL-01 a REQ-VAL-51 |

La clasificación observable entre rechazo, fallo, aviso, resultado estadístico y
estado no evaluable es exclusivamente la de REQ-VAL-01 a REQ-VAL-51 y sus
ACE; esta tabla no introduce otra política.

### 10.2. Registro observable

Los componentes emiten eventos estructurados a un colector común. El sobre
técnico contiene identidad de ejecución, fase, componente, entidad, gravedad y
correlación causal. La carga útil procede del componente responsable. Los
mensajes para personas se renderizan desde ese registro y no son la única fuente
de diagnóstico.

La observación de tiempo, memoria, almacenamiento y cantidades de dominio se
sitúa en las fronteras que producen o consumen esos recursos. La selección de
métricas, política de datos sensibles y umbrales de rendimiento se mantiene en
las decisiones de planificación, sin modificar los criterios funcionales.

## 11. Estrategia de comprobación

La matriz
[trazabilidad-requisitos-0.1.md](trazabilidad-requisitos-0.1.md) es el índice
obligatorio entre cada REQ, sus ACE, el hito y la evidencia prevista. Esta
sección solo organiza las capas técnicas de comprobación.

| Capa | Propósito técnico | Evidencia |
|---|---|---|
| Dominio unitario | Probar tipos, transformaciones e invariantes sin E/S ni procesos externos | Casos TEST de cada REQ biológico y sus ACE |
| Frontera de componente | Probar traducción de identidades, coordenadas, configuración y errores con dobles de dependencias | Casos TEST de REQ-CONF-01 a REQ-CONF-123, REQ-REF-01 a REQ-REF-23 y REQ-VAL-01 a REQ-VAL-51, además de la familia del componente |
| Adaptador | Comprobar lectura, escritura e invocación contra contratos externos | Casos TEST de REQ-FMT-01 a REQ-FMT-07, REQ-DAT-01 a REQ-DAT-73, REQ-VER-01 a REQ-VER-47, REQ-MAN-01 a REQ-MAN-88 y REQ-SAL-01 a REQ-SAL-26 |
| Propiedad y estadística | Explorar espacios amplios sin confundir variación aleatoria con un invariante | Casos TEST asociados a REQ-MUT-01 a REQ-MUT-15, REQ-EVO-01 a REQ-EVO-73, REQ-REP-01 a REQ-REP-26, REQ-REC-01 a REQ-REC-29 y REQ-DAT-01 a REQ-DAT-73 |
| Recorrido completo | Demostrar que las fronteras componen un paquete trazable | ACE e hitos enlazados por la matriz |

Los validadores independientes y versiones concretas se eligen en la
planificación. La prueba del productor nunca sustituye por sí sola a una
comprobación independiente cuando el hito exige interoperabilidad.

## 12. Recursos y evolución prevista

Las interfaces de secuencia, haplotipo, lecturas y artefactos admiten acceso por
tramos y escritura incremental. La representación interna, el tamaño de bloque,
la concurrencia y la política de caché se eligen con medidas; no forman parte
del contrato funcional. El estimador y el observador de recursos comparten
unidades, pero permanecen separados para poder comparar estimación y uso real.

Los puntos de extensión se sitúan en adaptadores de formatos y herramientas,
estrategias de dominio expresamente aprobadas y productores de artefactos. Una
extensión no puede entrar por una rama genérica y eludir identidad, validación,
trazabilidad o versionado.

Las capacidades funcionales fuera de v1 se registran únicamente en
[registro-ampliaciones-futuras-0.1.md](registro-ampliaciones-futuras-0.1.md);
las decisiones técnicas aún pendientes se mantienen en la sección 10 de
[decisiones-tecnicas-0.1.md](decisiones-tecnicas-0.1.md). Esta arquitectura no
anticipa esos comportamientos con entidades implícitas ni los presenta como
obligaciones vigentes.

## 13. Trazabilidad de requisitos a responsables

| Familia de requisitos | Responsable principal | Colaboradores y comprobación |
|---|---|---|
| REQ-CONF | Admisión de configuración | Identidad, coordenadas, aleatoriedad y validación prueban precedencias y rechazos |
| REQ-TEC, REQ-RG | Orquestación de eventos técnicos | WGS, artefactos, alineamiento y manifiestos prueban aislamiento e identidad técnica |
| REQ-FMT | Ensamblado de resultados | Adaptadores de formato y validación externa prueban orden, índice y referencia |
| REQ-REF | Catálogo de referencias | Coordenadas, integridad y exportación prueban FASTA, linaje e índices necesarios |
| REQ-IND, REQ-DIR | Motor de verdad genómica y creación directa | Configuración y materialización prueban ploidía, fase, eventos y reconstrucción |
| REQ-EVO, REQ-POB | Dinámica evolutiva e inicialización poblacional | Reproducción, selección, muestreo y validación poblacional |
| REQ-CRU | Grafo de familias | Reproducción, selección declarativa, verdad y validación del grafo |
| REQ-REP | Reproducción y dinámica evolutiva | Historial y validación comprueban N, Ne, orden y generaciones |
| REQ-SEL | Dinámica evolutiva: selección de viabilidad | Motor genómico y validación comprueban locus, dosis, función y supervivencia |
| REQ-MIG | Dinámica evolutiva: migración | Selección declarativa e historial prueban asignación conjunta y pertenencia |
| REQ-MUT | Motor de verdad genómica | Evolución, familias y validación prueban orden, tasas y transmisión |
| REQ-REC | Reproducción y recombinación | Coordenadas, motor genómico y salida detallada prueban bivalentes y segmentos |
| REQ-DAT, REQ-WGS | Simulación WGS y alineamiento | Artefactos, modo directo y validación técnica prueban origen y coherencia |
| REQ-CONT | Artefactos y demultiplexado | Procedencia y validación prueban sustitución, transferencias y totales |
| REQ-VER | Materialización de verdad | Motor genómico, referencias y validación VCF prueban fase y aislamiento |
| REQ-VAL | Validación interna | Todos los productores aportan objetivos, hechos y contexto sin autocorregirse |
| REQ-MAN | Manifiestos y exportación | Identidad y ensamblado aportan el catálogo de hechos y artefactos |
| REQ-SAL | Ensamblado de resultados | Validación y manifiestos deciden completitud y exportabilidad |
| REQ-NOF | Estimación y observación de recursos | Todos los componentes exponen medidas sin imponer límites funcionales fijos |

Los criterios ACE-01 a ACE-111 se distribuyen siguiendo esta matriz. Antes de
planificar cada fase, sus criterios se convertirán en casos concretos y se
enlazarán con el componente y la frontera que aportan la evidencia.

## 14. Decisiones abiertas y condición de aprobación

Las decisiones confirmadas, propuestas, pendientes y futuras se registran en
[decisiones-tecnicas-0.1.md](./decisiones-tecnicas-0.1.md). En particular, esta
arquitectura no decide silenciosamente crates, formatos o esquemas, versiones de
estándares, plataformas, algoritmo aleatorio, garantías entre plataformas,
reglas de paralelismo ni objetivos cuantitativos de rendimiento.

Para cerrar F0 deben aprobarse esta división y, como mínimo, resolverse las
decisiones marcadas como bloqueantes de F1. Para cerrar una fase posterior deben
resolverse las decisiones que cambien sus contratos antes de redactar su plan de
implementación.

La arquitectura se considera comprobable cuando:

- cada requisito de la fase tiene un responsable y una evidencia prevista;
- ninguna frontera permite mezclar referencias, identidades o ámbitos;
- verdad y evidencia técnica recorren caminos separados;
- los fallos y resultados parciales tienen un estado no ambiguo;
- un plan de implementación puede dividir el trabajo sin inventar nuevos
  componentes, contratos o decisiones técnicas.
