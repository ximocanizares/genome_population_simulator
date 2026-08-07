# Especificación del simulador de genomas y poblaciones

- **Versión:** 0.1
- **Estado:** borrador actualizado
- **Fecha:** 7 de agosto de 2026
- **Fuentes disponibles:** `README.md`, información confirmada por el usuario y
  la decisión posterior de usar Rust recogida en
  [hoja-de-ruta-0.1.md](./hoja-de-ruta-0.1.md)

## 1. Estado y finalidad del documento

Este documento define el objetivo general y los requisitos confirmados hasta
ahora. También identifica las decisiones que faltan antes de poder cerrar la
especificación y planificar la implementación.

La versión 0.1 sigue siendo un borrador. Los requisitos de esta versión son
obligatorios, pero todavía no se han acordado todos sus detalles ni todos los
casos que permitirán comprobar los resultados.

La sección 4 es la fuente normativa de requisitos funcionales. Las secciones 2,
5 y 6 los resumen y la sección 7 propone criterios de aceptación preliminares.
Si una redacción resumida difiere de un requisito identificado, prevalece el
requisito y debe corregirse el resumen. Esta regla aclara la fuente que debe
consultarse; no autoriza a ignorar contradicciones entre requisitos.

El núcleo del sistema se desarrollará en Rust. Esta elección ya está confirmada
para el producto y no determina por sí sola la arquitectura, los *crates*, la
versión del lenguaje, los formatos, los esquemas ni los algoritmos.

## 2. Objetivo del proyecto

El sistema debe generar una colección coherente de datos genómicos simulados a
partir de un fichero general de configuración. Cuando se use la vía evolutiva o
la vía de familias y poblaciones segregantes, el fichero general referenciará
además un fichero de configuración poblacional propio del escenario. La colección
generada estará formada por una o más referencias, cada una guardada en un único
fichero FASTA, los estados genotípicos completos que sean obligatorios para la
vía utilizada, el historial necesario para reconstruir su procedencia,
copias de las configuraciones utilizadas, un manifiesto general y, cuando haya
eventos, un manifiesto por evento. Cada salida asociada a una observación debe
solicitarse mediante un evento
explícito. Dentro de un evento que active alineamientos, el valor predeterminado
será generar CRAM mediante el modo realista; este modo conserva un FASTQ final
—o una pareja de FASTQ— y un fichero estructurado de procedencia por grupo de
lectura antes de mapearlos. La configuración podrá solicitar un VCF de verdad
sin activar FASTQ, BAM ni CRAM, o usar un modo directo de alineamientos
sintéticos que no requiere FASTQ.

El uso principal es preparar datos para probar programas de detección y
genotipado de SNP. Por ello, cuando un evento solicita alineamientos sin
sustituir los valores predeterminados, simula las lecturas desde los haplotipos,
aplica la contaminación y los errores de demultiplexado, conserva el FASTQ final
resultante y mapea esas mismas lecturas para generar CRAM. Una ejecución sin
eventos no genera por omisión ninguna salida asociada a eventos y conserva el
resultado mínimo obligatorio.

El sistema debe poder representar por separado la contaminación biológica,
producida al mezclar material o lecturas de otras fuentes, y la contaminación
técnica entre muestras procesadas en una misma tanda de secuenciación o conjunto
de bibliotecas multiplexadas. Multiplexar consiste en secuenciar juntas varias
bibliotecas y distinguir después sus muestras mediante índices.

Cada evento que genere FASTQ, BAM o CRAM constituye un experimento de
secuenciación independiente. Usa una única tecnología, un único perfil y,
cuando genera alineamientos, un único mapeador, y produce sus propias lecturas y
resultados. Un evento que solo materialice verdad genotípica o VCF no requiere
tecnología, perfil, biblioteca ni mapeador. Una misma ejecución biológica puede
contener varios eventos con tecnologías o mapeadores distintos, pero no
reutiliza FASTQ entre ellos ni tiene como objetivo comparar mapeadores.

La primera versión simulará exclusivamente secuenciación del genoma completo
(WGS, *whole-genome sequencing*) de DNA genómico. Podrá usar lecturas cortas
tipo Illumina, *single-end* o *paired-end*, y lecturas largas PacBio HiFi, en
experimentos separados. La preparación de biblioteca predeterminada será sin
PCR; también podrá seleccionarse una preparación WGS con amplificación, que
puede producir duplicados de PCR y otros sesgos propios de esa preparación.
Oxford Nanopore y PacBio CLR quedan para futuras ampliaciones con el fin de
reducir el alcance inicial y aplazar perfiles de error más complejos. Esta
decisión de alcance no presupone que una tecnología moderna sea incapaz de
detectar SNP.

Antes de aplicar sesgos o degradaciones, cada copia o haplotipo contribuirá a las
lecturas en proporción a su dosis y a la ploidía. Las diferencias de soporte
entre alelos aparecerán únicamente por el muestreo aleatorio, la cobertura, los
errores, el mapeo y los demás efectos ya admitidos. La primera versión no forzará
un sesgo de biblioteca a favor de un alelo ni la pérdida completa de un alelo.

El simulador será agnóstico del organismo: no incorporará supuestos funcionales
fijos sobre una especie concreta. La primera versión priorizará ejemplos
pequeños, controlados y reproducibles para probar casos concretos. Los escenarios
grandes o complejos deberán seguir siendo posibles por diseño, pero su
optimización será una prioridad secundaria.

La referencia inicial podrá ser simulada por el programa según la configuración
o proceder de un FASTA real suministrado como entrada. Un genoma podrá tener uno
o varios cromosomas o *scaffolds*. Un *scaffold* es una secuencia que organiza
y orienta fragmentos ensamblados y que puede contener regiones cuya secuencia
exacta no se conoce. La cantidad y el tamaño de estas secuencias, y por tanto el
tamaño total del genoma, podrán variar mediante la configuración cuando la
referencia inicial sea simulada. El único FASTA de cada referencia podrá
contener uno o varios de estos cromosomas o *scaffolds*.

Una referencia estructural es una nueva referencia creada a partir de otra
referencia disponible previamente. Solo se diferencia de su referencia de
origen por las duplicaciones, deleciones, inversiones y translocaciones
solicitadas. El sistema debe poder mantener disponibles dos o más referencias
estructurales y conservar el origen y los cambios de cada una.

La verdad conocida es la información exacta, definida durante la simulación,
sobre los genotipos reales de los individuos simulados. El VCF de verdad debe
construirse desde ese estado genotípico y recoger esos genotipos cuando se
solicite como salida; nunca debe inferirlos a partir de FASTQ, BAM o CRAM. Cada
muestra corresponde a un individuo simulado. Su estructura puede corresponder a
cualquiera de las referencias disponibles, pero su BAM o CRAM puede estar
mapeado contra otra referencia disponible.

La ploidía será configurable globalmente y por cromosoma o *scaffold*. Para una
secuencia dada, todos los individuos de la ejecución tendrán la misma ploidía.

Los individuos podrán contener SNP, inserciones o deleciones pequeñas (*indels*),
variantes estructurales y repeticiones en tándem respecto a su referencia.
Podrán crearse directamente, mediante la evolución configurada de una población
inicial o mediante uno o más eventos de autofecundación o cruzamiento sexual.

La coherencia de las salidas no exige que todos los archivos utilicen una misma
referencia. En este documento, un conjunto es la agrupación formada por los BAM
o CRAM de todas las muestras incluidas y, cuando se solicite, el VCF de verdad
correspondiente, expresados respecto a una única referencia de mapeo. La
configuración podrá seleccionar una o varias referencias disponibles. Para cada
referencia seleccionada se generará un conjunto independiente y todos los
conjuntos contendrán las mismas muestras. Los FASTQ del modo realista pertenecen
al evento y, dentro de él, se mapean contra todas las referencias seleccionadas.
Nunca se mezclarán en un conjunto archivos o datos que utilicen referencias
distintas. Los identificadores biológicos de muestra se mantienen entre
conjuntos, pero los registros y los identificadores de variante son locales a
cada referencia y no se relacionan registro a registro.

El manifiesto debe dejar constancia de lo que ocurrió en cada ejecución. Debe
permitir conocer la vía de creación, los ficheros de configuración utilizados,
los parámetros realmente aplicados, el origen y los cambios de las referencias,
y los nombres y las rutas de los archivos generados.

## 3. Fuera del alcance del proyecto

Este proyecto se limita a generar los datos simulados y su verdad conocida. No
ejecutará programas externos de detección y genotipado de SNP, denominados *SNP
callers*. Tampoco comparará sus resultados con la verdad conocida, calculará
medidas de eficacia o fiabilidad ni evaluará herramientas externas.

La ejecución de esas herramientas, la comparación de resultados y su evaluación
pertenecerán a otro proyecto. Ese proyecto podrá utilizar los datos y el VCF de
verdad generados por este simulador, pero no forma parte de esta
especificación.

También quedan fuera de alcance la llamada de variantes y genotipos, la
imputación de genotipos ausentes, el faseado inferido y la corrección de
genotipos a partir de lecturas o alineamientos. Esas operaciones corresponden a
programas externos. La cobertura, la contaminación, los errores, la calidad y
la ambigüedad del mapeo solo afectan a FASTQ, BAM o CRAM y no modifican la verdad
genotípica simulada.

En esta primera versión quedan fuera de alcance las diferencias de ploidía del
mismo cromosoma o *scaffold* entre individuos, incluidos los modelos de
cromosomas sexuales con ploidía individual distinta y las aneuploidías
individuales. Podrán estudiarse como trabajo futuro.

También quedan fuera de alcance los datos crudos del instrumento de
secuenciación, sus señales, los archivos BCL y la generación detallada de
lecturas de índice i5 o i7. La contaminación técnica se representará mediante
el efecto final del demultiplexado, no mediante la simulación completa del
instrumento.

En esta versión también quedan fuera de alcance las lecturas quiméricas formadas
por fragmentos de individuos o bibliotecas diferentes. Los intercambios de
lecturas entre muestras se representan exclusivamente mediante contaminación o
errores de demultiplexado.

Quedan para futuras implementaciones los paneles dirigidos, la secuenciación de
amplicones obtenidos por PCR, la genotipificación por secuenciación (GBS,
*genotyping-by-sequencing*) y los demás protocolos de captura o reducción de
representación. Los duplicados de PCR de una biblioteca WGS amplificada sí
pertenecen al alcance: son un artefacto de su preparación, no una secuenciación
dirigida por PCR.

También se aplazan los modelos explícitos de desequilibrio alélico de biblioteca
y pérdida alélica completa (*allele dropout*). El desequilibrio alélico es un
soporte de lecturas desproporcionado respecto a la dosis; la pérdida alélica
completa es la ausencia de lecturas de un alelo que sí existe en el genotipo.
No deben confundirse con variaciones aleatorias del soporte ni con el
desequilibrio de ligamiento entre loci.

La primera versión usa generaciones discretas y no solapadas. No modela edades,
supervivencia entre generaciones ni reproducción entre cohortes distintas.
Tampoco incluye selección por fecundidad, selección poligénica o epistasis,
elección de pareja por parentesco, distancia genética o genotipo, ni tasas de
mutación SNP dependientes del contexto de bases vecinas.

En este modelo, un cambio demográfico de N no elimina individuos de una
generación que ya existe. Fija el número objetivo de descendientes que formarán
la generación siguiente. Los individuos de la generación actual solo pueden
actuar como progenitores y dejan de ser población activa cuando los sustituye su
descendencia.

La primera versión tampoco modela un componente ambiental del carácter
cuantitativo ni un parámetro de heredabilidad. La heredabilidad describe qué
parte de la variación observada de un carácter se asocia con diferencias
genéticas en una población. Estos modelos quedan para futuras ampliaciones. En
esta versión, el valor del carácter usado para seleccionar coincide con el valor
genético asignado al genotipo o a la combinación de dosis.

La selección cuantitativa de la primera versión es exclusivamente direccional:
puede favorecer valores mayores o valores menores del carácter. Quedan para
futuras ampliaciones la selección estabilizadora, que favorece valores próximos
a un óptimo intermedio, y la selección disruptiva, que favorece ambos extremos
frente a los valores intermedios.

Los elementos transponibles (TE) se introducirán mediante eventos solicitados y
conocidos para estudiar su efecto sobre programas externos. No se modelarán
copias activas o inactivas, tasas de movilización ni transposición espontánea a
lo largo de las generaciones. Las inserciones ya presentes sí se heredarán,
segregarán y recombinarán como las demás variantes.

En reproducción sexual genérica quedan para futuras versiones las ploidías
impares, los multivalentes, la doble reducción y la preferencia de apareamiento
por subgenoma. También se aplazan la interferencia entre entrecruzamientos y los
gametos desequilibrados producidos por recombinación en regiones estructuralmente
incompatibles. Los modelos específicos de herencia polisómica o disómica para
aplicar F también quedan para futuras ampliaciones.

La primera versión solo modela entrecruzamientos como recombinación. La
conversión génica queda para futuras versiones y no debe confundirse con una
translocación no recíproca. También quedan fuera la distorsión de segregación y
el impulso meiótico (*meiotic drive*), que alteran la probabilidad de transmitir
alelos o copias y son distintos de la selección de viabilidad.

La reproducción de la primera versión se limita a autofecundación y cruzamiento
sexual. La reproducción clonal, la apomixis y la producción de dobles haploides
quedan para futuras ampliaciones. Tampoco se modelan cruces recíprocos, herencia
citoplasmática ni efectos maternos o paternos; los parentales son simétricos y no
se asignan sexos.

La composición de variantes estructurales anidadas o de eventos que modifican
una estructura creada previamente en la misma copia haplotípica y su linaje
queda para versiones futuras. Esta exclusión no impide eventos no solapados o
adyacentes en un mismo haplotipo ni solapamientos situados en haplotipos
distintos.

Quedan también para futuras versiones las inserciones antes de la primera base o
después de la última; las deleciones terminales o de secuencias completas; las
inversiones terminales o de secuencias completas; las translocaciones terminales
o de cromosomas o *scaffolds* completos; y los cambios de ploidía o pérdidas
cromosómicas asociados. En una duplicación parcial dentro de la misma secuencia,
también queda fuera insertar la copia dentro de su intervalo fuente. Esta
prohibición no afecta a la excepción confirmada para una fuente que abarque la
secuencia completa. También se aplazan la generación de bloques de LD a partir
de estadísticos objetivo y el cálculo interno de D, D' o r².

El VCF de verdad predeterminado no pretende describir posiciones invariantes ni
emitir bloques gVCF. La representación completa de genotipos y haplotipos
corresponde al fichero estructurado obligatorio de verdad.

## 4. Requisitos confirmados

### 4.1. Ficheros de configuración de entrada

Un fichero estructurado organiza la información en campos y grupos que un
programa puede identificar y leer, en lugar de presentarla como texto libre. No
se fija aquí su formato ni su esquema, es decir, las reglas que definirán sus
campos, tipos y relaciones. Esta decisión se tomará durante la planificación de
la implementación teniendo en cuenta que el núcleo se desarrollará en Rust.

El fichero general contiene las decisiones que afectan a la totalidad de la
ejecución. Un fichero de configuración poblacional describe un escenario
concreto de la vía evolutiva o de la vía de familias y poblaciones segregantes.

- **REQ-CONF-01.** El sistema debe leer los parámetros de cada ejecución desde
  un fichero general de configuración estructurado y legible por el programa.
- **REQ-CONF-02.** El fichero general debe indicar si la referencia inicial será
  simulada o suministrada mediante un FASTA real.
- **REQ-CONF-03.** Cuando la referencia inicial sea suministrada, el fichero
  general debe permitir indicar el archivo FASTA de origen.
- **REQ-CONF-04.** Cuando la referencia inicial sea simulada, el fichero general
  debe permitir configurar el número de cromosomas o *scaffolds*, el tamaño de
  cada uno y la composición general del genoma.
- **REQ-CONF-05.** El sistema debe obtener el tamaño total del genoma a partir de
  los tamaños configurados para sus cromosomas o *scaffolds*.
- **REQ-CONF-06.** El fichero general debe permitir configurar una ploidía global
  predeterminada y sustituciones específicas por cromosoma o *scaffold*. Cada
  valor debe ser un número entero entre 1 y 10, inclusive.
- **REQ-CONF-07.** El fichero general debe seleccionar una de las tres vías de
  creación definidas en la sección 4.3.
- **REQ-CONF-08.** El fichero general debe permitir definir cambios concretos o
  solicitar cambios al azar. Para SNP, *indels* y variantes estructurales por
  separado, debe permitir elegir entre una cantidad media esperada de eventos o
  una densidad por base, kilobase o megabase.
- **REQ-CONF-09.** Cuando se seleccione la vía evolutiva o la vía de familias y
  poblaciones segregantes, el fichero general debe referenciar el fichero de
  configuración poblacional correspondiente.
- **REQ-CONF-10.** Cada escenario de la vía evolutiva o de la vía de familias y
  poblaciones segregantes debe tener su propio fichero de configuración
  poblacional estructurado y legible por el programa.
- **REQ-CONF-11.** El fichero de un escenario evolutivo debe definir cómo se
  obtiene la población final desde la inicial mediante una secuencia ordenada de
  fases evolutivas.
- **REQ-CONF-12.** El fichero de un escenario de familias o poblaciones
  segregantes debe indicar los parentales, los cruces, su orden o secuencia, las
  relaciones necesarias y las generaciones.
- **REQ-CONF-13.** Todos los parámetros configurables que se acuerden en
  versiones posteriores deben proporcionarse mediante el fichero general o el
  fichero poblacional que corresponda.
- **REQ-CONF-14.** La persona usuaria debe poder variar los parámetros sin
  modificar el código del programa.
- **REQ-CONF-15.** El fichero general debe permitir seleccionar una o varias
  referencias de mapeo entre las referencias disponibles.
- **REQ-CONF-16.** Todos los parámetros de generación aleatoria deben
  proporcionarse mediante el fichero general o el fichero de configuración
  poblacional que corresponda.
- **REQ-CONF-17.** Los parámetros de generación aleatoria deben poder definirse
  como valores globales y, cuando se desee, como valores específicos para un
  individuo.
- **REQ-CONF-18.** Un valor específico de un individuo debe sustituir al valor
  global solo para ese individuo. Los demás individuos deben conservar el valor
  global.
- **REQ-CONF-19.** En una población, la configuración debe poder indicar una
  frecuencia alélica como proporción poblacional. Para un individuo concreto no
  debe usar una frecuencia alélica, sino la presencia o ausencia del alelo, su
  genotipo o su dosis alélica, según la ploidía.
- **REQ-CONF-20.** La elección entre cantidad media esperada y densidad, el
  valor y, cuando proceda, la unidad de densidad deben poder definirse como
  valores globales y como sustituciones específicas por individuo, con la
  precedencia establecida en REQ-CONF-18.
- **REQ-CONF-21.** Para SNP, *indels* y variantes estructurales por separado, la
  configuración debe permitir elegir una distribución espacial uniforme o
  definir regiones con una cantidad media esperada o una densidad de cambios
  mayor o menor que el valor general, incluidas regiones con densidad cero.
- **REQ-CONF-22.** La distribución espacial y las regiones deben poder definirse
  como valores globales y como sustituciones específicas por individuo, con la
  precedencia establecida en REQ-CONF-18.
- **REQ-CONF-23.** Toda cantidad configurada directamente debe interpretarse
  como una cantidad media esperada, no como un número exacto de eventos. Cuando
  se configure una densidad, el sistema debe obtener la cantidad media esperada
  multiplicando la densidad por la longitud válida del ámbito al que se aplica.
- **REQ-CONF-24.** El número realizado de eventos debe obtenerse al azar a partir
  de la cantidad media esperada mediante una distribución probabilística y puede
  variar entre ejecuciones.
- **REQ-CONF-25.** La distribución de Poisson debe ser la predeterminada para
  obtener el número realizado. La configuración debe permitir seleccionar otras
  distribuciones probabilísticas y sus parámetros por separado para cada clase
  de variante y para cada región configurada, respetando las sustituciones por
  individuo definidas en REQ-CONF-18.
- **REQ-CONF-26.** Toda ejecución debe usar una semilla aleatoria, es decir, un
  valor que determina las decisiones tomadas al azar. El fichero general debe
  permitir proporcionar una semilla. Si no se proporciona, el sistema debe
  generarla antes de iniciar la simulación.
- **REQ-CONF-27.** Al repetir una ejecución con la misma semilla, las mismas
  configuraciones, las mismas entradas y la misma versión del programa, el
  sistema debe reproducir las mismas decisiones aleatorias y la misma simulación
  biológica. Esta es la garantía mínima de reproducibilidad y corresponde al
  nivel biológico definido en REQ-CONF-111.
- **REQ-CONF-28.** Cada región de generación de variantes debe identificar una
  referencia estructural concreta, uno de sus cromosomas o *scaffolds* y un
  intervalo de coordenadas basado en 1 y cerrado, es decir, con los dos extremos
  incluidos. El formato técnico del intervalo todavía no está definido.
- **REQ-CONF-29.** Los valores espaciales deben seguir, como mínimo, la jerarquía
  de ámbito global del genoma, cromosoma o *scaffold* y región. Para un parámetro
  y una clase de variante, el valor del ámbito más específico debe sustituir al
  del ámbito más general que lo contiene.
- **REQ-CONF-30.** Si dos regiones del mismo nivel se solapan y asignan valores
  incompatibles al mismo parámetro y a la misma clase de variante, el sistema
  debe rechazar la configuración antes de iniciar la simulación. El diagnóstico
  debe identificar las dos regiones y el parámetro en conflicto. Un
  solapamiento con valores compatibles no debe considerarse un conflicto.
- **REQ-CONF-31.** La precedencia espacial de REQ-CONF-29 debe aplicarse junto
  con la precedencia por individuo de REQ-CONF-18. El sistema debe determinar el
  valor efectivo después de considerar ambas jerarquías.
- **REQ-CONF-32.** Cuando una población combine un objetivo de frecuencia
  alélica con presencia, ausencia, genotipos o dosis alélicas fijados para
  individuos concretos, estos valores individuales deben tratarse como
  restricciones obligatorias y tener prioridad sobre el objetivo poblacional.
- **REQ-CONF-33.** Cada fase evolutiva debe indicar, como mínimo, su número de
  generaciones, los valores o cambios de N objetivo y de Ne objetivo, el
  modelo de apareamiento y los parámetros evolutivos que se usarán durante la
  fase. El Ne objetivo puede constar como no configurado; el Ne realizado no es
  una entrada.
- **REQ-CONF-34.** Una fase debe poder sustituir valores generales de la
  configuración durante su propio intervalo. Al terminar la fase, esas
  sustituciones no deben aplicarse a las fases siguientes salvo que su propia
  configuración así lo indique.
- **REQ-CONF-35.** El fichero de configuración poblacional evolutiva debe
  permitir solicitar uno o varios muestreos durante la evolución. Cada muestreo
  debe indicar, como mínimo, la generación y su tamaño, expresado como número de
  individuos o como proporción de la población existente en esa generación.
- **REQ-CONF-36.** Los muestreos evolutivos deben interpretarse como
  observaciones no destructivas. Su presencia o ausencia en la configuración no
  debe cambiar la trayectoria evolutiva cuando el resto de la configuración, las
  entradas, la semilla y la versión del programa sean iguales.
- **REQ-CONF-37.** Cada evento de muestreo debe permitir solicitar expresamente
  la generación de BAM o CRAM y de VCF de verdad para sus individuos. Si no se
  solicita, esas salidas no deben generarse para el evento.
- **REQ-CONF-38.** Los individuos definidos por la persona usuaria deben poder
  aportar su identificador. Antes de simular, el sistema debe rechazar la
  configuración si contiene identificadores de individuo duplicados y señalar
  claramente los duplicados.
- **REQ-CONF-39.** Una población evolutiva debe poder contener uno o varios
  subpoblaciones. Cada subpoblación debe tener un identificador único dentro de la ejecución, y
  ese identificador no debe reutilizarse aunque la subpoblación se extinga o se fusione.
  La configuración debe rechazarse antes de simular si repite un identificador.
- **REQ-CONF-40.** Para la población completa y, cuando corresponda, para cada
  subpoblación, N objetivo y Ne objetivo deben poder configurarse de forma independiente
  y variar mediante valores por generación o mediante fases o tramos
  configurados. Las reglas de obligatoriedad, herencia y expresión de ambos
  valores se definen en la sección 4.3.5.
- **REQ-CONF-41.** La configuración debe permitir definir migración mediante una
  relación dirigida por cada par de subpoblaciones origen y destino. Cada fase
  o intervalo debe elegir de forma excluyente proporciones o números exactos de
  migrantes para todos sus pares; no debe mezclar ambos modos. Cada ruta debe
  seleccionar migrantes mediante el catálogo declarativo común de
  REQ-CONF-108.
- **REQ-CONF-42.** La configuración debe permitir programar, en generaciones o
  fases concretas, eventos de división de una subpoblación, fusión de varias subpoblaciones,
  fundación de una subpoblación nueva mediante individuos seleccionados y extinción de un
  subpoblación. La división y la fundación deben usar el catálogo declarativo
  común de REQ-CONF-108 y realizar la selección sin reemplazo.
  Para una fusión debe permitir indicar dos o más orígenes y una única
  subpoblación nueva. Debe permitir definir N objetivo, Ne objetivo, selección,
  autofecundación, apareamiento, migración y los demás parámetros aplicables a
  las subpoblaciones resultantes. Para la contribución de una fusión debe elegir
  un modo probabilístico con pesos por origen o un modo controlado con recuentos
  parentales exactos por origen. Ambos modos son excluyentes. Una recolonización debe declarar un
  identificador nuevo, distinto del de cualquier subpoblación extinguida.
- **REQ-CONF-43.** La configuración evolutiva debe permitir activar, desactivar y
  parametrizar por fase, generación y subpoblación, cuando corresponda, la selección, la
  deriva genética, la mutación y la migración. También debe poder incluir otros
  procesos evolutivos que el sistema admita.
- **REQ-CONF-44.** El fichero de configuración poblacional debe permitir
  reordenar los procesos opcionales dentro de una generación solo cuando se
  respeten las dependencias biológicas y estructurales obligatorias. Si no se
  indica un orden, debe aplicarse la secuencia predeterminada de REQ-REP-24.
- **REQ-CONF-45.** Antes de simular, el sistema debe rechazar una secuencia que
  contenga procesos desconocidos o combinaciones incompatibles y explicar los
  procesos y el conflicto detectados.
- **REQ-CONF-46.** La configuración debe elegir para toda la ejecución uno de
  dos modos excluyentes de recombinación: una tasa uniforme para todo el genoma,
  expresada en centimorgan por megabase (cM/Mb), o un mapa completo con
  coordenada física y posición genética acumulada en centimorgan (cM).
- **REQ-CONF-47.** Cada mapa debe identificar explícitamente la referencia
  estructural y contener todos sus cromosomas o *scaffolds* recombinantes. Antes
  de simular, el sistema debe rechazar mapas incompletos, identificadores ajenos
  y coordenadas o posiciones genéticas no válidas.
- **REQ-CONF-48.** No debe permitirse combinar una tasa uniforme para unas
  secuencias y mapas para otras dentro de la misma ejecución. Los perfiles
  generales o específicos de grupo deben respetar el único modo elegido y, si
  usan mapas, aportar un mapa completo.
- **REQ-CONF-49.** La configuración debe permitir definir grupos con
  identificadores y miembros, y asignar a cada grupo una tasa o un mapa de
  recombinación específico que sustituya al perfil general para sus individuos.
  Un grupo puede representar sexo, tipo reproductivo u otra clasificación, pero
  el sistema no debe imponerle un significado biológico.
- **REQ-CONF-50.** Cada individuo debe pertenecer a cero o a un único grupo de
  recombinación. Antes de simular, el sistema debe rechazar una configuración que
  asigne un individuo a más de uno e identificar el individuo y todos los grupos
  en conflicto.
- **REQ-CONF-51.** Los cruces complejos deben representarse como un grafo
  dirigido acíclico de eventos. Cada evento debe tener un identificador único,
  sus parentales o fuentes, su generación u orden y el número de descendientes.
- **REQ-CONF-52.** Antes de simular un grafo de cruces, el sistema debe validar
  que los identificadores sean únicos, las fuentes referenciadas existan, el
  orden temporal sea coherente y no haya ciclos. Cualquier error debe detener la
  ejecución e identificar los eventos y relaciones afectados. Reutilizar un
  parental existente en una autofecundación, un retrocruzamiento o varias ramas
  no debe considerarse un ciclo si cada evento produce descendientes nuevos.
- **REQ-CONF-53.** Para usar descendientes como parentales posteriores, la
  configuración debe permitir indicar identificadores concretos; selección
  aleatoria; pertenencia a una población o grupo configurado; y presencia de un
  alelo, genotipo, dosis o haplotipo configurado en uno o varios loci. Estas
  condiciones deben poder combinarse mediante Y, O y NO. La cantidad solicitada
  debe poder expresarse como número exacto o como proporción del conjunto
  elegible.
- **REQ-CONF-54.** Cada población intermedia de un grafo de cruces debe permitir
  solicitar expresamente BAM o CRAM y VCF de verdad para todos sus individuos.
  Si no se solicitan, esos archivos intermedios no deben generarse.
- **REQ-CONF-55.** El fichero general debe permitir modificar el umbral de
  longitud que distingue inserciones y deleciones pequeñas de variantes
  estructurales. Si no se indica, el sistema debe usar 50 pares de bases.
- **REQ-CONF-56.** La configuración debe permitir definir loci concretos de
  repeticiones en tándem y generarlos al azar mediante
  las mismas reglas globales, regionales, individuales y reproducibles de la
  generación de variantes.
- **REQ-CONF-57.** Para cada repetición en tándem, la configuración debe poder indicar
  el motivo repetido, la longitud de la unidad, el número de copias o la longitud
  total y la pureza de la repetición.
- **REQ-CONF-58.** La configuración debe permitir variar las repeticiones en tándem entre
  haplotipos, individuos y poblaciones, incluidos el número de copias y las
  interrupciones del motivo.
- **REQ-CONF-59.** La pureza y la variabilidad deben poder definirse para un locus
  concreto o mediante distribuciones aplicables a conjuntos de loci, con valores
  globales y sustituciones específicas cuando corresponda.
- **REQ-CONF-60.** La longitud del motivo debe ser completamente configurable y
  no debe limitarse a 1–6 pb. El sistema debe rechazar un motivo vacío o una
  longitud declarada que no coincida con la longitud de su secuencia, e
  identificar el locus y el error antes de simular.
- **REQ-CONF-61.** La configuración debe permitir definir para las repeticiones
  en tándem la tasa de mutación, el sesgo entre expansión y contracción y la
  distribución del número de unidades ganadas o perdidas. Debe poder seleccionar
  otros modelos admitidos además del modelo escalonado predeterminado.
- **REQ-CONF-62.** El modelo y sus parámetros deben poder definirse globalmente,
  por locus o región y, en la vía evolutiva, por fase o subpoblación, aplicando las
  precedencias espaciales, individuales y temporales ya establecidas.
- **REQ-CONF-63.** Las mutaciones internas que introducen, eliminan o modifican
  interrupciones y cambian la pureza o el motivo deben estar desactivadas de
  forma predeterminada. Si se activan, la configuración debe permitir definir
  sus tipos, tasas y distribuciones globalmente y por locus, fase o subpoblación según
  las precedencias aplicables.
- **REQ-CONF-64.** Todo evento que solicite alineamientos debe generarlos de
  forma predeterminada en formato CRAM. La configuración del evento debe
  permitir elegir BAM como alternativa o no solicitar alineamientos. Una
  ejecución sin eventos no debe crear un evento implícito ni generar
  alineamientos por omisión. El VCF de verdad debe ser opcional y debe poder
  solicitarse sin FASTQ, BAM ni CRAM. Ninguna de estas salidas, ni los FASTQ
  exigidos por el modo realista, forma parte del mínimo obligatorio de toda
  ejecución.
- **REQ-CONF-65.** Antes de simular, el sistema debe comprobar que cada secuencia
  tiene una ploidía efectiva. Debe rechazar valores no enteros, menores que 1 o
  mayores que 10, o secuencias sin valor efectivo, e identificar la secuencia y
  el problema.
- **REQ-CONF-66.** Cuando la salida BAM o CRAM esté activa, la configuración
  debe permitir elegir una tecnología o perfil de lecturas predefinido y
  sustituir sus parámetros. La primera versión debe incluir perfiles
  predeterminados validados para lecturas tipo Illumina *single-end* y
  *paired-end* y para PacBio HiFi.
- **REQ-CONF-67.** Cuando la salida BAM o CRAM esté activa, la configuración debe
  permitir elegir entre dos modos: el modo realista predeterminado, que simula y
  conserva los FASTQ finales antes de mapearlos, y el modo directo opcional, que
  produce alineamientos sintéticos sin generar ni mapear FASTQ previamente.
- **REQ-CONF-68.** En el modo realista, los FASTQ deben activarse obligatoriamente
  junto con BAM o CRAM. En el modo directo no deben ser necesarios. Una
  configuración incompatible con estas dependencias debe rechazarse antes de
  generar datos.
- **REQ-CONF-69.** Un experimento de mapeo debe seleccionar exactamente un
  mapeador, es decir, un programa que alinea lecturas contra una referencia. El
  predeterminado debe ser minimap2. La configuración puede elegir otro mapeador
  compatible y sus parámetros, pero no varios dentro del mismo experimento.
- **REQ-CONF-70.** Cuando el mapeador sea minimap2, el perfil de tecnología debe
  seleccionar automáticamente el *preset* y los parámetros predeterminados
  apropiados para lecturas cortas tipo Illumina o lecturas largas PacBio HiFi.
  La configuración debe permitir sustituirlos y validar, cuando sea posible,
  las opciones desconocidas o incompatibles, además de advertir de cada
  sustitución.
- **REQ-CONF-71.** Si se seleccionan varias referencias de mapeo dentro de un
  experimento, el único mapeador seleccionado debe aplicarse a todas. La
  primera versión no debe reutilizar FASTQ entre eventos para comparar
  mapeadores; cada evento debe generar sus propias lecturas.
- **REQ-CONF-72.** La cobertura debe configurarse como una cobertura media
  objetivo en «x», de forma global y con sustituciones por individuo,
  cromosoma o *scaffold* y región. Debe admitir valores altos, bajos o cero.
- **REQ-CONF-73.** Las sustituciones de cobertura deben aplicar la precedencia ya
  acordada: región sobre cromosoma o *scaffold* y este sobre el valor global,
  junto con la sustitución por individuo cuando corresponda. Dos valores
  incompatibles en ámbitos solapados del mismo nivel deben rechazarse antes de
  simular, identificando los ámbitos y el conflicto.
- **REQ-CONF-74.** Cada región de cobertura debe identificar la referencia
  estructural del individuo, el cromosoma o *scaffold* y el intervalo sobre el
  que se interpreta.
- **REQ-CONF-75.** Una muestra debe admitir cero, una o varias fuentes de
  contaminación sustitutiva. Para cada fuente, la configuración debe expresar
  la contribución como proporción, cobertura equivalente o número de lecturas o
  parejas que sustituirán datos propios, con una unidad explícita. Cuando se
  combinen unidades distintas, todas las contribuciones deben convertirse de
  forma conjunta a las unidades indivisibles de la biblioteca antes de elegir
  las lecturas: lecturas individuales para *single-end* o lecturas largas y
  parejas completas para *paired-end*.
- **REQ-CONF-76.** Cuando una fuente contaminante se cuantifique por número en
  datos *paired-end*, el valor debe contar pares completos, no lecturas sueltas.
  La configuración de esta operación no debe permitir solicitar pares
  huérfanos.
- **REQ-CONF-77.** En el modo realista, una fuente contaminante puede ser un
  individuo o una población simulados o uno o varios FASTQ externos. En el modo
  directo solo puede ser un individuo o población simulados con verdad conocida.
  La configuración debe identificar sin ambigüedad la fuente, su tecnología,
  perfil y disposición *single-end* o *paired-end*.
- **REQ-CONF-78.** La contaminación técnica debe configurarse por tanda de
  secuenciación o conjunto de bibliotecas multiplexadas, identificando sus
  muestras, bibliotecas e índices. Debe distinguirse expresamente de la
  contaminación biológica definida por fuentes de lecturas.
- **REQ-CONF-79.** La configuración debe admitir índices simples, índices duales
  combinatorios e índices duales únicos (UDI, *Unique Dual Indexes*), en los que
  cada biblioteca usa una combinación exclusiva de dos índices.
- **REQ-CONF-80.** La tasa y procedencia de la asignación errónea de lecturas por
  índices, incluido el *index hopping* o asociación accidental de una lectura
  con índices de otra biblioteca, deben ser configurables y reproducibles.
- **REQ-CONF-81.** La estrategia de índices debe permitir configurar si las
  combinaciones inesperadas reconocibles se filtran o quedan sin asignar durante
  el demultiplexado, además de representar las lecturas que se reasignan a otra
  muestra. El demultiplexado es la separación de las lecturas multiplexadas por
  muestra a partir de sus índices.
- **REQ-CONF-82.** Cada tecnología de secuenciación debe tener un perfil
  predeterminado identificado y versionado que incluya, como mínimo, la tasa y
  los tipos de error habituales que pretende representar. No debe usarse un
  perfil «estándar» implícito o sin versión.
- **REQ-CONF-83.** El perfil de errores debe aplicarse de forma homogénea a todos
  los individuos, cromosomas o *scaffolds* y regiones cuando la configuración no
  proporcione sustituciones.
- **REQ-CONF-84.** La configuración debe permitir sustituir el perfil de errores
  globalmente y por individuo, cromosoma o *scaffold* y región. Debe aplicar la
  misma precedencia por especificidad y el mismo rechazo de conflictos del mismo
  nivel establecidos para la cobertura.
- **REQ-CONF-85.** La configuración debe permitir proporcionar opcionalmente un
  FASTQ real de WGS de DNA genómico y la referencia correspondiente para estimar
  un perfil empírico de simulación WGS.
- **REQ-CONF-86.** La configuración debe distinguir entre estimar un perfil
  empírico y seleccionarlo como perfil efectivo. La calibración no debe sustituir
  al perfil tecnológico predeterminado salvo que la persona usuaria lo seleccione
  expresamente.
- **REQ-CONF-87.** La configuración debe permitir seleccionar como entrada un
  fichero de perfil empírico generado por una ejecución anterior, sin requerir
  de nuevo el FASTQ real del que se estimó.
- **REQ-CONF-88.** La preparación WGS predeterminada debe ser sin PCR. Un perfil
  WGS amplificado debe permitir configurar por separado la tasa de duplicados de
  PCR, originados al representar varias veces una misma molécula amplificada, y
  la de duplicados ópticos, originados al detectar varias veces una misma señal.
- **REQ-CONF-89.** Cada perfil predeterminado debe aportar valores habituales
  para los mecanismos de duplicación que le correspondan, y la configuración
  debe permitir sustituirlos. Un mecanismo no debe aplicarse a una tecnología
  para la que el perfil lo considere no pertinente.
- **REQ-CONF-90.** La configuración debe permitir elegir si los duplicados
  simulados se entregan marcados con la bandera estándar de duplicado en BAM o
  CRAM o se entregan sin marcar. Si se omite esta elección, deben entregarse
  marcados.
- **REQ-CONF-91.** La configuración debe permitir activar de forma opcional la
  presencia de restos de adaptadores sin recortar. Esta degradación debe estar
  desactivada de forma predeterminada y solo aplicarse por solicitud expresa.
- **REQ-CONF-92.** Cuando se active, la configuración debe permitir definir como
  mínimo el perfil o la secuencia del adaptador, la frecuencia o proporción de
  lecturas afectadas y, cuando corresponda, la longitud residual y el extremo de
  la lectura afectado.
- **REQ-CONF-93.** Los perfiles tecnológicos pueden aportar parámetros
  recomendados para restos de adaptadores, pero no deben activar implícitamente
  la degradación.
- **REQ-CONF-94.** Cada perfil tecnológico debe poder definir una curva o
  relación predeterminada entre el contenido GC de la secuencia y la probabilidad
  de obtener lecturas. La configuración debe permitir modificar su intensidad o
  anular este sesgo.
- **REQ-CONF-95.** La configuración efectiva debe determinar de forma definida y
  documentada cómo se combina el sesgo de contenido GC con la cobertura media
  objetivo y con sus sustituciones regionales explícitas.
- **REQ-CONF-96.** La configuración debe permitir definir o generar por regiones
  secuencias repetitivas, duplicadas o altamente homólogas destinadas a producir
  baja mapeabilidad. También debe poder aprovechar las repeticiones y variantes
  estructurales ya presentes en los genomas simulados.
- **REQ-CONF-97.** En el modo directo, la configuración debe permitir definir la
  proporción de lecturas correctamente alineadas, alineadas en una posición
  incorrecta, multimapeadas y no alineadas, junto con una distribución de MAPQ
  para cada categoría. Los valores predeterminados deben proceder del perfil
  tecnológico cuando corresponda.
- **REQ-CONF-98.** Antes de generar alineamientos directos, el sistema debe
  comprobar que las proporciones sean válidas, no negativas y formen una
  distribución completa, que las MAPQ estén dentro de los valores admitidos por
  el formato y que las combinaciones solicitadas de campos y alineamientos sean
  coherentes. Debe rechazar la configuración e identificar el parámetro que
  incumpla estas condiciones.
- **REQ-CONF-99.** La configuración debe permitir generar lecturas quiméricas,
  formadas por fragmentos que no eran contiguos en el genoma o haplotipo de
  origen, con una frecuencia configurable. Los perfiles tecnológicos que
  contemplen este mecanismo deben aportar valores predeterminados. La
  configuración debe poder controlar si se combinan cromosomas o *scaffolds* y
  haplotipos o copias cromosómicas distintos. En esta versión no debe permitir
  combinar fragmentos de individuos o bibliotecas diferentes.
- **REQ-CONF-100.** La configuración debe definir un umbral general para decidir
  qué *scaffolds* reciben un mínimo de un entrecruzamiento por bivalente y
  meiosis. Debe elegir una sola unidad para esa regla: longitud genética en cM o
  longitud física en bases. No debe mezclar ambas unidades dentro de la misma
  regla. Debe permitir sustituir el valor general para *scaffolds* concretos,
  conservando la unidad elegida.
- **REQ-CONF-101.** Una misma ejecución biológica debe poder definir cero, uno o
  varios eventos de observación o muestreo. Cada salida asociada a una
  observación, incluido un VCF, FASTQ, BAM o CRAM, debe pertenecer a un evento
  explícito con identificador único, selección de individuos, referencias
  aplicables, salidas y subsemilla. Cero eventos significa que no se genera
  ninguna salida asociada a eventos; no debe crearse un evento implícito.
- **REQ-CONF-102.** Un evento que genere FASTQ, BAM o CRAM debe seleccionar
  exactamente una tecnología y un perfil tecnológico. Cuando genere BAM o CRAM
  debe seleccionar además exactamente un mapeador conforme al modo configurado.
  Un evento que solo solicite verdad genotípica o VCF no debe requerir
  tecnología, perfil, bibliotecas ni mapeador. Eventos distintos pueden usar
  valores distintos sin compartir lecturas ni resultados técnicos.
- **REQ-CONF-103.** La configuración debe permitir solicitar de forma
  independiente, para cada evento, su VCF de verdad, sus FASTQ y sus BAM o CRAM,
  respetando las dependencias entre el modo realista, el modo directo y las
  salidas solicitadas.
- **REQ-CONF-104.** La configuración evolutiva debe indicar si permite que una
  subpoblación se extinga cuando un resultado de la simulación produce N
  realizado igual a 0. El valor global predeterminado debe ser desactivado y
  debe poder sustituirse por subpoblación y fase evolutiva; el ámbito más
  específico debe prevalecer. Esta opción es distinta de programar un evento
  explícito de extinción.
- **REQ-CONF-105.** Si una generación contiene varios eventos demográficos, la
  configuración debe indicar su orden total de forma explícita y sin
  ambigüedades. No debe inferirse un orden a partir de cómo aparezcan campos no
  ordenados en el fichero.
- **REQ-CONF-106.** Cuando un conjunto de frecuencias o probabilidades deba sumar
  1, el sistema debe comparar la suma con 1 mediante una tolerancia general
  configurable. Si la desviación está dentro de la tolerancia y la suma permite
  normalizar, debe dividir los valores por su suma, registrar los valores antes
  y después y continuar. Si la desviación supera la tolerancia o la suma no
  permite normalizar, debe rechazar la configuración. El valor predeterminado
  exacto de esta tolerancia se fijará durante la implementación.
- **REQ-CONF-107.** La tolerancia general de normalización solo debe corregir
  pequeñas desviaciones de suma. No debe sustituir tolerancias con otra finalidad
  ya acordadas, incluida `0,5/C` para aproximar frecuencias objetivo mediante
  recuentos enteros en el modo controlado.
- **REQ-CONF-108.** La división, la fundación, la migración y la selección de
  descendientes deben usar un único catálogo declarativo de reglas. El catálogo
  debe admitir identificadores concretos, selección aleatoria, pertenencia a una
  población o grupo configurado y presencia de un alelo, genotipo, dosis o
  haplotipo configurado en uno o varios loci. Las condiciones deben poder
  combinarse mediante Y, O y NO, y la cantidad debe expresarse como un número
  exacto o una proporción.
- **REQ-CONF-109.** Una proporción del catálogo común debe interpretarse como la
  probabilidad individual de seleccionar cada elegible y producir un recuento
  aleatorio reproducible. Un número exacto debe seleccionar sin reemplazo y
  fallar si no hay suficientes elegibles. Deben registrarse regla, cantidad o
  proporción solicitada, elegibles, seleccionados y cantidad o proporción
  realizada. Las reglas deben evaluarse sobre la verdad biológica, no sobre
  evidencia técnica posterior.
- **REQ-CONF-110.** Las reglas del catálogo común deben ser declarativas, es
  decir, describir condiciones que el programa conoce y puede validar. La
  primera versión no debe ejecutar scripts, fragmentos de código ni funciones
  personalizadas. Añadir una clase de regla distinta debe requerir ampliar el
  programa y su especificación, no introducir código desde la configuración.

- **REQ-CONF-111.** La reproducibilidad debe declararse y comprobarse en cuatro
  niveles separados: decisiones biológicas; lecturas y decisiones técnicas
  aleatorias; equivalencia semántica de las salidas; e identidad binaria de los
  ficheros. Alcanzar un nivel no debe implicar automáticamente los siguientes.
- **REQ-CONF-112.** El nivel biológico debe garantizar las mismas decisiones y
  la misma verdad biológica con la misma semilla, configuraciones efectivas,
  entradas, perfiles y versión del simulador. Las decisiones de muestreo y de
  cada evento deben usar flujos o subsemillas aislados para que añadir o retirar
  una observación técnica no altere la trayectoria biológica.
- **REQ-CONF-113.** El nivel de lecturas y decisiones técnicas aleatorias debe
  garantizar las mismas secuencias, calidades, procedencias y artefactos
  simulados cuando coincidan la subsemilla del evento, la versión del simulador,
  las configuraciones efectivas, las entradas, los perfiles y las dependencias
  que intervienen en esas decisiones.
- **REQ-CONF-114.** La equivalencia semántica debe evaluarse sobre el contenido
  lógico de cada salida y no sobre sus bytes. Para salidas producidas por una
  herramienta externa solo podrá declararse cuando coincidan su versión,
  parámetros y las demás condiciones que el simulador haya verificado como
  compatibles. El manifiesto debe declarar las condiciones comprobadas y si el
  nivel se alcanzó o no se garantiza.
- **REQ-CONF-115.** La identidad binaria no debe garantizarse de forma general.
  Solo podrá declararse para una clase de salida y unas condiciones concretas si
  se ha verificado expresamente. Las diferencias limitadas a compresión,
  ordenación equivalente o metadatos no biológicos no invalidan por sí mismas
  los niveles anteriores.

- **REQ-TEC-01.** Cada experimento de secuenciación y cada evento que genere
  FASTQ, BAM o CRAM debe seleccionar exactamente una tecnología. No debe mezclar
  lecturas cortas tipo Illumina y lecturas largas PacBio HiFi dentro del mismo
  evento o experimento. Un evento que solo produzca verdad genotípica o VCF no
  constituye un experimento de secuenciación.
- **REQ-TEC-02.** Un experimento puede contener varias bibliotecas, tandas o
  carriles, pero todos deben usar el mismo tipo tecnológico y ser compatibles
  con el perfil de lecturas y el mapeo seleccionados.
- **REQ-TEC-03.** La configuración debe permitir crear experimentos separados
  por tecnología a partir de la misma simulación biológica, incluidos eventos de
  muestreo distintos, con trazabilidad entre ellos y sin mezclar sus archivos ni
  reutilizar sus lecturas.
- **REQ-TEC-04.** El mapeador, el *preset* y sus parámetros deben ser
  compatibles con la única tecnología del experimento y aplicarse de forma
  coherente a todas sus bibliotecas, tandas, carriles y referencias.
- **REQ-TEC-05.** Dos perfiles que requieran modelos de error, calidades o
  parámetros de mapeo distintos deben considerarse tecnologías experimentales
  diferentes aunque pertenezcan al mismo fabricante. Cada tecnología admitida
  debe tener un identificador y una versión explícitos en la configuración.
- **REQ-TEC-06.** Todo experimento de la primera versión debe representar WGS de
  DNA genómico. La configuración debe rechazar paneles dirigidos, amplicones,
  GBS y otros protocolos de captura o reducción de representación.
- **REQ-TEC-07.** Los perfiles WGS deben distinguir preparación sin PCR y
  preparación con amplificación de biblioteca. Si no se indica otra, debe usarse
  la preparación sin PCR.
- **REQ-TEC-08.** Un perfil WGS con amplificación puede introducir duplicados de
  PCR y otros sesgos propios, con parámetros predeterminados del perfil y
  sustituciones configurables. Esto no debe activar ni representar una
  secuenciación dirigida por PCR.
- **REQ-TEC-09.** La primera versión solo debe admitir lecturas cortas tipo
  Illumina, *single-end* o *paired-end*, y lecturas largas PacBio HiFi. Debe
  rechazar Oxford Nanopore y PacBio CLR como tecnologías de simulación o
  calibración e indicar que quedan aplazadas a futuras ampliaciones.
- **REQ-TEC-10.** La configuración efectiva y el manifiesto deben registrar el
  tipo de preparación WGS, el perfil, sus valores predeterminados y todas las
  sustituciones aplicadas.

- **REQ-FMT-01.** La configuración debe permitir elegir CRAM o BAM como formato
  del alineamiento final por individuo. Si se omite la elección, debe usarse
  CRAM; BAM debe ser una alternativa configurable.
- **REQ-FMT-02.** Todos los BAM o CRAM finales deben entregarse ordenados por
  coordenadas de forma predeterminada y acompañados de un índice compatible,
  como CRAI, BAI o CSI según el formato, el tamaño y las necesidades efectivas.
- **REQ-FMT-03.** La configuración efectiva debe registrar el formato, el orden
  y el tipo de índice de los alineamientos finales.
- **REQ-FMT-04.** Antes de considerar completo un resultado, el sistema debe
  validar que cada alineamiento final, su índice y su referencia sean compatibles
  y detener la ejecución e impedir que el conjunto afectado se marque como
  completo si no lo son.
- **REQ-FMT-05.** El VCF multimuesta de verdad debe entregarse de forma
  predeterminada ordenado según el orden de los cromosomas o *scaffolds* de la
  referencia y por coordenada dentro de cada secuencia, comprimido con BGZF e
  indexado mediante TBI o CSI según corresponda.
- **REQ-FMT-06.** La configuración debe permitir solicitar como alternativa un
  VCF de texto sin comprimir. La configuración efectiva y el manifiesto deben
  registrar el orden, la compresión y el tipo de índice aplicables.
- **REQ-FMT-07.** Antes de completar la ejecución, el sistema debe validar la
  coherencia del VCF, su compresión e índice cuando correspondan, y su cabecera y
  referencia. Una salida que no supere esta validación debe detener la ejecución
  e impedir que el conjunto afectado se marque como completo.

### 4.2. Referencias, FASTA y cambios estructurales

- **REQ-REF-01.** La referencia inicial debe poder obtenerse de dos formas:
  simulada por el programa según la configuración o suministrada como un FASTA
  real de entrada.
- **REQ-REF-02.** Un genoma de referencia debe poder contener uno o varios
  cromosomas o *scaffolds*, cuyos tamaños y tamaño total pueden variar según la
  configuración cuando el genoma sea simulado.
- **REQ-REF-03.** Cada genoma de referencia, tanto inicial como derivado, debe
  guardar toda su secuencia en un único fichero FASTA.
- **REQ-REF-04.** El único FASTA de una referencia debe poder contener uno o
  varios cromosomas o *scaffolds*.
- **REQ-REF-05.** El número de índices de una referencia puede ser cero, uno o
  varios.
- **REQ-REF-06.** El sistema debe generar para cada referencia todos los tipos de
  índice que necesiten las operaciones o los formatos utilizados.
- **REQ-REF-07.** El sistema debe conservar si la referencia inicial fue
  simulada o suministrada y, en el segundo caso, cuál fue el archivo FASTA de
  origen.
- **REQ-REF-08.** El sistema debe permitir crear una nueva referencia
  estructural a partir de otra referencia disponible previamente, incluida una
  referencia inicial simulada o suministrada. La referencia de partida se
  denomina referencia de origen.
- **REQ-REF-09.** La nueva referencia estructural debe diferenciarse de su
  referencia de origen únicamente por las duplicaciones, deleciones, inversiones
  y translocaciones solicitadas.
- **REQ-REF-10.** El sistema debe conservar la relación entre cada referencia
  estructural y su referencia de origen, junto con los cambios que las
  distinguen.
- **REQ-REF-11.** El sistema debe permitir que haya dos o más referencias
  estructurales disponibles.
- **REQ-REF-12.** Antes de iniciar la simulación, el sistema debe validar el
  FASTA real suministrado.
- **REQ-REF-13.** El FASTA suministrado debe ser legible, contener al menos una
  secuencia no vacía y tener identificadores de secuencia no vacíos y únicos.
  El identificador debe ser el primer elemento no vacío de la cabecera después
  de `>`, terminado por el primer espacio o carácter de separación. El resto de
  la cabecera debe conservarse por separado como descripción completa y no debe
  participar en la comprobación de unicidad.
- **REQ-REF-14.** Las secuencias del FASTA suministrado deben contener únicamente
  A, C, G, T y N. Se deben aceptar tanto mayúsculas como minúsculas. N es el
  único símbolo admitido para una base desconocida o no determinada.
- **REQ-REF-15.** El FASTA suministrado no debe contener símbolos que representen
  huecos (*gaps*), códigos IUPAC de ambigüedad distintos de N ni ningún otro
  símbolo diferente de los admitidos en REQ-REF-14.
- **REQ-REF-16.** Si el FASTA suministrado incumple cualquiera de estos
  requisitos, el sistema debe detener la ejecución antes de generar resultados
  y explicar claramente los errores detectados.
- **REQ-REF-17.** Una vez validado, el sistema debe copiar el FASTA suministrado
  al directorio de resultados, es decir, a la carpeta que reúne las salidas de
  la ejecución.
- **REQ-REF-18.** La ejecución debe usar la copia guardada en el directorio de
  resultados como referencia inicial y no debe modificar el archivo original.
- **REQ-REF-19.** Una referencia simulada debe poder incluir repeticiones en tándem
  definidas en loci concretos o generadas según la configuración aleatoria y
  regional, y su secuencia FASTA debe reflejar los motivos, copias e
  interrupciones realizados.
- **REQ-REF-20.** Las referencias deben poder contener regiones repetitivas,
  duplicadas o altamente homólogas definidas o generadas para evaluar baja
  mapeabilidad, incluidas las originadas por repeticiones o variantes
  estructurales ya configuradas.
- **REQ-REF-21.** Cada CRAM debe estar asociado inequívocamente al FASTA exacto
  de la referencia contra la que se mapeó. Ese FASTA debe almacenarse entre los
  recursos comunes de la ejecución y enlazarse desde el evento; no debe
  sustituirse por otra referencia. Al exportar el evento debe copiarse al paquete
  autocontenido.
- **REQ-REF-22.** Cada secuencia debe clasificarse explícitamente como cromosoma
  o *scaffold* en la configuración efectiva y en la salida estructurada.
- **REQ-REF-23.** Toda secuencia de DNA suministrada mediante la configuración,
  un catálogo, un motivo, una inserción u otra entrada estructurada, y toda
  secuencia generada por el simulador, debe usar exclusivamente A, C, G, T y N,
  sin distinguir mayúsculas de minúsculas durante la validación. Cualquier otro
  símbolo debe rechazarse antes de usar la secuencia e identificarse en el
  diagnóstico. El tratamiento de N por cada operación queda sujeto a la
  decisión funcional específica registrada como pendiente; aceptar N no
  autoriza a resolverla aleatoriamente de forma implícita.

La copia permite reunir en el directorio de resultados todos los archivos
necesarios para utilizar la simulación. Así se podrán exportar juntos a otros
programas sin depender de la ubicación original del FASTA.

Estos requisitos no seleccionan las herramientas ni los tipos concretos de
índice. La necesidad de cada tipo dependerá de las operaciones y los formatos
que se utilicen durante la implementación.

### 4.3. Individuos, muestras y grupos

En este documento, una muestra es el grupo de datos generado para un individuo
simulado. La ploidía de un cromosoma o *scaffold* es el número de copias de esa
secuencia representadas en un individuo.

#### 4.3.1. Requisitos comunes

Una variante estructural es un cambio que afecta a la estructura de una región
del genoma. Un par de bases (pb) es una unidad de longitud de DNA. Para
inserciones y deleciones, una longitud menor de 50 pb se clasifica de forma
predeterminada como *indel* y una longitud igual o mayor que 50 pb como variante
estructural. Otros tipos estructurales, como las inversiones o duplicaciones,
conservan su clasificación y no dependen únicamente de este umbral.

Una duplicación segmentaria copia una secuencia y la inserta dentro de una
secuencia ya existente. En la primera versión, la fuente puede ser un intervalo
parcial o la secuencia completa de un cromosoma o *scaffold*. Incluso cuando la
fuente es completa, la copia no crea otro cromosoma o *scaffold* independiente
ni aumenta la ploidía.

Las coordenadas escritas directamente en la configuración empiezan en 1 y
describen intervalos con ambos extremos incluidos. Los formatos externos
conservan su convenio oficial; el sistema debe identificarlo y convertirlo de
forma explícita. Una inserción se ancla en la base de referencia inmediatamente
anterior al punto donde se incorpora la secuencia nueva, de forma coherente con
VCF. La secuencia insertada no forma parte del intervalo original de referencia.

Una repetición en tándem es una región en la que un motivo de DNA se repite
consecutivamente. Los motivos de 1–6 pb pueden describirse como microsatélites o
repeticiones cortas en tándem; los motivos mayores permiten representar otros
tipos de repetición en tándem. La unidad repetida es una copia del motivo y el
número de copias determina la longitud principal de la repetición.
La pureza es la proporción del tracto repetido que sigue exactamente el motivo;
la parte restante corresponde a interrupciones. La variabilidad describe las
diferencias entre haplotipos, individuos o poblaciones en el número de copias y
en el patrón de interrupciones.

El modelo escalonado de mutación cambia una repetición en tándem mediante la
ganancia o pérdida de unidades completas del motivo. Es el modelo predeterminado,
pero no se fijan todavía las fórmulas de sus parámetros ni las distribuciones
alternativas admitidas.

- **REQ-IND-01.** El sistema debe crear individuos que contengan SNP, *indels*,
  variantes estructurales y repeticiones en tándem respecto a una referencia.
- **REQ-IND-02.** La estructura de cada individuo debe poder corresponder a
  cualquiera de las referencias disponibles.
- **REQ-IND-03.** El sistema debe conservar qué referencia disponible describe
  la estructura de cada muestra.
- **REQ-IND-04.** Los cambios deben poder crearse al azar configurando, para cada
  clase de variante por separado, una cantidad media esperada de eventos o una
  densidad por base, kilobase o megabase.
- **REQ-IND-05.** Los cambios también deben poder definirse de forma concreta en
  el fichero general de configuración.
- **REQ-IND-06.** El sistema debe poder simular ploidías enteras de 1 a 10
  copias cromosómicas, inclusive.
- **REQ-IND-07.** En una misma ejecución, todos los individuos deben tener la
  misma ploidía efectiva para un cromosoma o *scaffold* dado.
- **REQ-IND-08.** En esta primera aproximación, no se deben permitir cambios de
  ploidía entre individuos de una misma ejecución.
- **REQ-IND-09.** La cantidad media esperada o densidad de SNP, *indels* y
  variantes estructurales debe poder tener un valor global y valores específicos
  para los individuos que se indiquen.
- **REQ-IND-10.** Un valor específico de cantidad media esperada o densidad debe
  sustituir al valor global únicamente para el individuo correspondiente.
- **REQ-IND-11.** La frecuencia alélica debe representar la proporción de copias
  de un alelo en una población y no debe atribuirse literalmente a un solo
  individuo.
- **REQ-IND-12.** Para un individuo concreto, el sistema debe expresar el estado
  de un alelo mediante su presencia o ausencia, su genotipo o su dosis alélica,
  según la ploidía del individuo.
- **REQ-IND-13.** Para cada clase de variante por separado, el sistema debe poder
  distribuir espacialmente los eventos de manera uniforme o aplicar regiones
  configuradas con una cantidad media esperada o densidad mayor o menor que el
  valor general.
- **REQ-IND-14.** El sistema debe admitir regiones configuradas con densidad
  cero, en las que no se generen al azar eventos de la clase de variante
  correspondiente.
- **REQ-IND-15.** Una distribución o región específica de un individuo debe
  sustituir a la configuración global correspondiente solo para ese individuo.
  Los demás individuos deben conservar la configuración global.
- **REQ-IND-16.** Para cada clase de variante, ámbito general o región e
  individuo al que se aplique la generación aleatoria, el sistema debe conservar
  tanto la cantidad media esperada como el número realizado de eventos.
- **REQ-IND-17.** El número realizado debe ser entero y puede diferir de la
  cantidad media esperada y entre ejecuciones que no estén obligadas a repetir
  las mismas decisiones aleatorias.
- **REQ-IND-18.** El número realizado debe obtenerse con la distribución de
  Poisson de forma predeterminada. Cuando la configuración seleccione otra
  distribución probabilística para una clase de variante o una región, el
  sistema debe usar esa distribución y sus parámetros.
- **REQ-IND-19.** Las coordenadas de una región de generación deben interpretarse
  sobre la referencia estructural que corresponde al individuo o a la población
  a los que se aplica.
- **REQ-IND-20.** Si los datos de una muestra se mapean contra una referencia
  distinta de su referencia estructural, la región de generación no debe
  reinterpretarse en las coordenadas de la referencia de mapeo.
- **REQ-IND-21.** Al generar eventos, el sistema debe aplicar el valor espacial
  más específico disponible dentro de la jerarquía genoma, cromosoma o
  *scaffold* y región, junto con la sustitución específica del individuo cuando
  exista.
- **REQ-IND-22.** Cada locus con frecuencias objetivo debe usar exactamente un
  modo. El modo probabilístico, predeterminado, interpreta las frecuencias como
  probabilidades y muestrea la realización de forma estocástica después de
  respetar las restricciones individuales. El modo controlado solo debe
  aplicarse por solicitud explícita y elige recuentos enteros. Los dos modos son
  mutuamente excluyentes.
- **REQ-IND-23.** En el modo controlado y después de respetar las restricciones
  fijas, para un locus bialélico debe elegirse el recuento entero alcanzable más
  próximo al objetivo por C, donde C es el total de copias cromosómicas del
  ámbito. En un locus multialélico, los recuentos de todos los alelos deben
  optimizarse conjuntamente, sumar C y minimizar la desviación global respecto a
  sus objetivos. La tolerancia predeterminada, aplicable solo a este modo, debe
  ser 0,5/C y puede configurarse una más estricta. No alcanzar la tolerancia debe
  ser un error. Todo empate entre soluciones óptimas debe resolverse al azar con
  la semilla, de forma reproducible, para evitar un sesgo sistemático. Una vez
  fijados los recuentos, las copias alélicas deben distribuirse aleatoriamente
  entre individuos, condicionadas al modelo elegido: Hardy-Weinberg, F o una
  distribución explícita. Las restricciones incompatibles deben producir un
  error.
- **REQ-IND-24.** Para cada objetivo poblacional, el sistema debe registrar el
  modo efectivo. En el probabilístico debe conservar frecuencias esperadas y
  realizadas. En el controlado debe conservar C, objetivos, vector de recuentos,
  frecuencias realizadas, desviación global, tolerancia efectiva y restricciones.
  Si hubo empate, debe registrar las soluciones candidatas y la elección aleatoria.
- **REQ-IND-25.** Todo individuo debe tener un identificador único y estable
  dentro de la ejecución. El sistema debe asignarlo cuando el individuo no haya
  sido definido con uno y no debe reutilizarlo aunque el individuo deje de
  formar parte de la población.
- **REQ-IND-26.** El mismo identificador de individuo debe usarse de forma
  coherente en su estado genotípico, sus observaciones, sus relaciones parentales
  o de linaje y todos los archivos y metadatos que le correspondan.
- **REQ-IND-27.** Para cada evento de recombinación, el sistema debe determinar
  el perfil efectivo del individuo: el perfil específico de su grupo cuando sea
  aplicable o, en caso contrario, el perfil general.
- **REQ-IND-28.** El estado genotípico interno de cada individuo debe representar
  por separado cada copia cromosómica o haplotipo según su ploidía y conservar la
  fase hasta el final de la simulación.
- **REQ-IND-29.** Cada SNP, *indel* y variante estructural debe asignarse a uno o
  varios haplotipos concretos del individuo, según su genotipo.
- **REQ-IND-30.** Los haplotipos de un descendiente deben reconstruirse a partir
  de los haplotipos parentales y los eventos de recombinación, conservando el
  origen parental de sus segmentos.
- **REQ-IND-31.** El sistema debe aplicar el umbral efectivo de longitud de forma
  coherente al generar y clasificar inserciones y deleciones y al calcular sus
  cantidades y estadísticas por clase.
- **REQ-IND-32.** Los individuos deben poder diferir de su referencia y entre sí
  en el número de copias, la longitud total, la pureza y las interrupciones de
  las repeticiones en tándem configuradas.
- **REQ-IND-33.** El estado genotípico debe asignar cada alelo de repetición en tándem a
  un haplotipo concreto y conservar su fase, su secuencia realizada y el número
  de copias o longitud total.
- **REQ-IND-34.** La generación aleatoria o regional de repeticiones en tándem debe usar
  la semilla de la ejecución y conservar los parámetros y resultados realizados
  con el mismo nivel de trazabilidad que las demás clases de variación.
- **REQ-IND-35.** Cada evento de mutación de una repetición en tándem debe usar la
  semilla de la ejecución, asignarse a un haplotipo concreto y actualizar su
  estado genotípico con las unidades ganadas o perdidas.
- **REQ-IND-36.** Cuando estén activadas, las mutaciones internas de una
  repetición en tándem deben afectar a un haplotipo concreto y actualizar de
  forma trazable el patrón de interrupciones, la pureza o el motivo realizado.
- **REQ-IND-37.** Para cada cromosoma de un individuo, el estado interno debe
  representar y fasear exactamente el número de haplotipos que indique su
  ploidía efectiva, y conservar la dosis de cada alelo.
- **REQ-IND-38.** Las coordenadas introducidas directamente en la configuración
  deben ser de base 1 y usar intervalos cerrados, con ambos extremos incluidos.
  Esta regla debe aplicarse a regiones, variantes y mapas de recombinación.
- **REQ-IND-39.** Todo formato externo debe mantener su convenio oficial de
  coordenadas. La entrada, la configuración efectiva, la verdad y el manifiesto
  deben identificar el convenio y registrar toda conversión realizada.
- **REQ-IND-40.** Una inserción configurada debe anclarse en la base de
  referencia inmediatamente anterior al punto de inserción. Debe tener una base
  anterior válida y un punto interno; la primera versión debe rechazar
  inserciones antes de la primera base o después de la última.
- **REQ-IND-41.** Toda deleción de la primera versión debe ser interna y dejar al
  menos una base a cada lado. Deben rechazarse las deleciones que toquen el
  inicio o final de un cromosoma o *scaffold* y las que eliminen la secuencia
  completa, tanto si son explícitas como aleatorias.
- **REQ-IND-42.** Toda inversión de la primera versión debe ser interna y dejar
  al menos una base a cada lado. Deben rechazarse las inversiones terminales y
  las que abarquen la secuencia completa, tanto explícitas como aleatorias.
- **REQ-IND-43.** El intervalo fuente de una duplicación puede tocar uno de los
  extremos de la secuencia porque el original no se elimina. También puede
  abarcar la secuencia completa de un cromosoma o *scaffold*. El destino de la
  copia debe ser un punto interno válido de la secuencia receptora.
- **REQ-IND-44.** Una duplicación puede insertar la copia en el mismo cromosoma
  o *scaffold* o en otro del mismo haplotipo o genoma. Debe registrar origen,
  destino, orientación, identificador estable de copia y si es intra- o
  intersecuencia. Aunque la fuente abarque una secuencia completa, el evento
  debe seguir siendo una duplicación segmentaria insertada: no debe crear una
  secuencia independiente ni cambiar la ploidía.
- **REQ-IND-45.** Cuando una duplicación parcial tenga origen y destino en la
  misma secuencia, el destino no debe quedar dentro del intervalo fuente. Puede
  ser adyacente al intervalo, antes o después, o quedar en cualquier otra
  posición externa al intervalo, siempre que sea un punto de inserción interno
  válido. Como única excepción, si la fuente abarca la secuencia completa, la
  copia puede insertarse en un punto interno de esa misma secuencia y, por
  tanto, dentro del intervalo fuente.
- **REQ-IND-46.** La primera versión debe admitir translocaciones entre
  secuencias, tanto recíprocas, que intercambian segmentos, como no recíprocas,
  que trasladan un segmento. Cada translocación debe ser un evento estructural
  compuesto con identificador único y subeventos o puntos de ruptura relacionados.
- **REQ-IND-47.** Todos los puntos de ruptura de una translocación deben ser
  internos y ninguna secuencia puede quedar vacía. La orientación directa o
  invertida debe ser configurable por segmento y conservarse en la verdad, las
  coordenadas y la representación del evento.
- **REQ-IND-48.** Las translocaciones deben poder usarse al construir referencias
  o individuos y deben reflejarse en las lecturas generadas desde ellos. La
  verdad debe conservar sin pérdida el evento compuesto, sus segmentos,
  orientaciones, puntos de ruptura y relaciones.
- **REQ-IND-49.** Las duplicaciones locales deben aplicar un máximo total de
  copias distinto de la ploidía cromosómica y del número de unidades de una
  repetición en tándem. Si no se configura otro valor, el máximo provisional
  debe ser 10: una copia original y hasta nueve derivadas. Este valor debe poder
  aumentarse y revisarse mediante pruebas de rendimiento, porque todavía no se
  ha validado como límite práctico. No debe existir otro máximo absoluto fijo.
- **REQ-IND-50.** El máximo local debe aplicarse de forma acumulativa por
  intervalo fundador o familia de copias y por haplotipo, contando todos los
  eventos que compartan el mismo origen. Los identificadores de origen deben
  impedir que varios eventos eludan el límite. En un poliploide, cada copia
  haplotípica aplica su propio límite.
- **REQ-IND-51.** Antes de ejecutar, el sistema debe estimar el aumento de tamaño
  y recursos causado por las duplicaciones, incluidos memoria, almacenamiento y
  volumen de datos. Configurar un máximo local mayor que 10 debe generar una
  advertencia. Si la estimación supera el presupuesto de recursos configurado,
  el sistema debe rechazar la ejecución; en caso contrario, debe continuar.
  Superar el valor predeterminado no debe confundirse con superar la ploidía
  máxima de 10.
- **REQ-IND-52.** Cada evento de duplicación debe interpretar sus coordenadas de
  origen y destino en el estado de las secuencias inmediatamente anterior a su
  aplicación. Debe conservar la correspondencia entre ese estado, el intervalo
  fuente original, la copia insertada y las coordenadas posteriores cuando
  exista tal correspondencia. Esta regla no elimina las demás restricciones de
  composición y solapamiento de variantes estructurales.

El identificador de individuo distingue a una entidad biológica simulada. Es
distinto del identificador de un evento de muestreo, que distingue una
observación, y del identificador de un conjunto, que agrupa archivos referidos a
la misma referencia de mapeo. La sintaxis de estos identificadores todavía no se
ha definido.

Un haplotipo es la secuencia de variantes que se encuentran en una misma copia
cromosómica. La fase identifica qué alelos pertenecen juntos a cada copia. En
este documento, conservar la fase significa mantener esa asignación explícita en
el estado interno y a través de la herencia y la recombinación.

La cantidad media esperada es el promedio de eventos que se obtendría al repetir
muchas veces la generación aleatoria con los mismos parámetros. No exige que
cada ejecución produzca ese número exacto y puede no ser un número entero. El
número realizado es el número entero de eventos obtenido en una ejecución
concreta.

La densidad relaciona la cantidad media esperada con una longitud de secuencia y
puede expresarse por base, kilobase o megabase. Una kilobase equivale a mil bases
y una megabase, a un millón de bases. La longitud válida de un ámbito es la
longitud en la que se pueden situar eventos de la clase correspondiente después
de aplicar las restricciones pertinentes. La cantidad media esperada derivada
de una densidad se obtiene multiplicando esa densidad por la longitud válida del
ámbito. Cada clase de variante puede usar su propio modo y, si usa densidad, su
propia unidad.

Una distribución probabilística asigna probabilidades a los posibles números
realizados a partir de la cantidad media esperada. La distribución de Poisson es
la opción predeterminada. Podrán seleccionarse otras distribuciones por clase de
variante y por región, pero el catálogo y los parámetros admitidos todavía no
se han definido.

Una distribución espacial uniforme reparte los eventos sin favorecer unas
posiciones válidas sobre otras dentro del ámbito que se determine. Una región
configurable es una parte de la secuencia en la que la cantidad media esperada o
la densidad puede aumentar, disminuir o ser cero respecto al valor general. Cada
clase de variante puede tener su propia distribución espacial y sus propias
regiones. La distribución probabilística del número y la distribución espacial
de las posiciones son decisiones distintas.

La frecuencia alélica es la proporción de copias de un alelo en una población.
No debe confundirse con la cantidad media esperada ni con la densidad de
eventos. La dosis alélica es el número de copias de un alelo que contiene un
individuo y depende de su ploidía.

Todavía deben definirse el formato técnico de las regiones, el catálogo y los
parámetros de las distribuciones probabilísticas adicionales, y cómo se
comprobará estadísticamente que los números realizados concuerdan con las
cantidades medias esperadas.

La ploidía puede variar entre cromosomas o *scaffolds* de 1 a 10 copias, pero la
configuración efectiva de cada secuencia es idéntica para todos los individuos de
la ejecución.

#### 4.3.2. Vía 1: creación directa de individuos

- **REQ-DIR-01.** El sistema debe permitir crear directamente individuos con
  los cambios indicados en el fichero general de configuración.
- **REQ-DIR-02.** La creación directa debe admitir tanto cambios concretos como
  cambios al azar configurados por clase de variante mediante una cantidad
  media esperada o una densidad por base, kilobase o megabase, con valores
  globales y valores específicos por individuo.
- **REQ-DIR-03.** La creación directa debe admitir para cada clase de variante
  una distribución espacial uniforme o regiones con una cantidad media
  esperada o densidad mayor o menor que el valor general, incluidas regiones con
  densidad cero, tanto globales como específicas por individuo.
- **REQ-DIR-04.** En la creación directa o aleatoria de individuos, las
  diferencias de número de copias de las repeticiones en tándem deben poder
  generarse con el modelo escalonado predeterminado o con otro modelo admitido
  seleccionado en la configuración.
- **REQ-DIR-05.** La creación directa debe aplicar mutaciones internas de
  repeticiones en tándem solo cuando se hayan activado expresamente.

#### 4.3.3. Vía 2: población inicial y evolución configurada

El desequilibrio de ligamiento (LD) es la asociación no aleatoria entre alelos
de loci diferentes. El equilibrio de Hardy-Weinberg es una expectativa
probabilística que relaciona las frecuencias alélicas y genotípicas bajo unos
supuestos determinados. Será un modo opcional de inicialización, no una condición
obligatoria. Una población también podrá inicializarse fuera del equilibrio con
distribuciones genotípicas explícitas. Estas condiciones iniciales son distintas
de las desviaciones que aparezcan después por autofecundación, endogamia,
selección, deriva u otros procesos evolutivos.

El coeficiente de endogamia inicial F es una opción de comodidad para aumentar
la homocigosis esperada respecto al equilibrio de Hardy-Weinberg. Para cualquier
ploidía mayor que 1, F=0 conserva las probabilidades genotípicas de ese
equilibrio. Al aumentar F, debe crecer progresivamente la probabilidad de que
las copias homólogas de un individuo porten el mismo alelo. F=1 produce
homocigosis completa en cada locus. En este documento, F es una expectativa
genómica de la población o de una subpoblación, no un valor exacto por
individuo, locus o bloque ni una longitud esperada de tramos de homocigosis
(ROH, *runs of homozygosity*). La primera versión no distingue modelos
específicos de herencia polisómica o disómica para aplicar F. El F realizado se
estima comparando la heterocigosidad observada con la esperada según las
frecuencias alélicas, por subpoblación y generación registrada. En poliploides,
un individuo se considera heterocigoto cuando sus copias homólogas no portan
todas el mismo alelo. El resultado puede diferir del F configurado y ser
negativo cuando hay más heterocigotos de los esperados; ese valor se conserva y
no constituye un error.

La endogamia es el apareamiento entre individuos emparentados. La tasa de
mutación describe con qué frecuencia aparecen cambios genéticos nuevos. Los
cuatro conceptos anteriores son ejemplos, no una lista cerrada de propiedades
o parámetros poblacionales. Esta observación no añade modelos al alcance de v1;
cualquier propiedad nueva requerirá ampliar expresamente la especificación.

La selección es el proceso por el que las diferencias hereditarias pueden
cambiar la contribución de los individuos o genotipos a generaciones futuras.
La deriva genética es el cambio aleatorio de las frecuencias alélicas entre
generaciones debido al muestreo reproductivo, cuya intensidad se relaciona con
N y Ne.

La panmixia es un modelo en el que los apareamientos se producen al azar dentro
de la población. Un cuello de botella es una reducción marcada del tamaño de la
población durante una parte de su historia; en v1 reduce el N objetivo de la
generación siguiente. Una expansión lo aumenta. Ninguno de
los dos eventos elimina retrospectivamente individuos adultos de una generación
ya creada. Una fase evolutiva es un intervalo de una o más generaciones que
comparte el tamaño o cambio demográfico, el modelo de apareamiento y los demás
parámetros configurados para ese intervalo.

El tamaño N objetivo es el número de descendientes que el sistema intenta crear
para una población o subpoblación. N realizado es el tamaño censal obtenido, es
decir, el número de individuos que se crearon realmente. Ambos pueden diferir si
la migración, la selección u otro resultado aleatorio dejan progenitores
insuficientes durante la ejecución. El
tamaño efectivo, representado por Ne, es un parámetro que gobierna la intensidad
de la deriva genética y la endogamia y puede diferir de N. En cada fase y
generación se distinguen N objetivo, N realizado, Ne objetivo y Ne realizado. El
Ne objetivo es el valor estadístico que el esquema reproductivo intenta
aproximar cuando se ha configurado. El Ne realizado es el valor estimado después
de la reproducción. En la primera versión será exclusivamente un Ne
reproductivo: se estimará a partir de los progenitores que dejaron descendencia
y de la desigualdad o varianza de sus contribuciones. No se estimará a partir de
la pérdida de heterocigosidad, el incremento de endogamia ni el desequilibrio de
ligamiento; esos estimadores quedan para futuras ampliaciones. Una población puede
dividirse en subpoblaciones, que son grupos reproductivos internos con parámetros
demográficos o evolutivos propios y que pueden estar conectados por migración.
En este documento, «población» designa siempre el conjunto completo simulado. Si
solo existe un grupo reproductivo, la población contiene una única
subpoblación. La migración es el paso de un individuo de una subpoblación a otra.
N objetivo guía la descendencia y N realizado registra su resultado. Ne
objetivo puede expresarse como valor absoluto o relativo a N, o constar como no
configurado. Ne realizado se calcula siempre.
La migración se aplica mediante rutas dirigidas entre pares de subpoblaciones y
se resuelve una sola vez dentro de cada generación.

Una división de subpoblación reparte los individuos progenitores de una
subpoblación de origen entre dos o más subpoblaciones hijas. Estos individuos se
consideran fundadores de la hija asignada. Después del reparto, cada hija genera
su propia descendencia conforme a su N objetivo y a sus parámetros
reproductivos. El número de progenitores recibidos no determina ese objetivo.

Una fusión reúne todos los progenitores de dos o más subpoblaciones de origen en
una subpoblación nueva. Las subpoblaciones originales dejan de existir después
del evento. La nueva subpoblación conserva la procedencia de cada individuo y
haplotipo, usa sus propios parámetros y genera después su descendencia. Su N
objetivo es independiente del número de progenitores reunidos.

En una fusión, una contribución parental es la participación de un progenitor
de un origen en la creación de descendencia. La ascendencia genómica realizada
es la proporción del genoma descendiente que procede de cada origen, reconstruida
desde los haplotipos transmitidos. No son la misma medida: la recombinación y la
ascendencia previa de los progenitores pueden hacer que difieran.

Una fundación transfiere algunos progenitores de una subpoblación de origen a
una subpoblación nueva sin eliminar la original. Los individuos transferidos
solo pertenecen al conjunto reproductor nuevo durante esa generación. La nueva
subpoblación genera descendencia desde sus fundadores y la original lo hace
desde los individuos que no fueron transferidos. El N objetivo de la nueva
subpoblación es independiente de su número de fundadores.

Una extinción planificada es un evento configurado para una generación
concreta. Una extinción por resultado ocurre cuando la simulación genera N
realizado igual a 0 y la configuración permite esa posibilidad. Ambas terminan
la subpoblación y conservan su historia, pero deben registrarse como causas
distintas.

La recombinación es el intercambio de material entre copias de cromosomas
durante la formación de descendientes. Un mapa de recombinación describe cómo
varía su tasa a lo largo de una referencia estructural.

En los resultados de un muestreo, el ámbito estadístico es el conjunto de
individuos que describe una cifra. El ámbito poblacional comprende la población
completa o una subpoblación. El ámbito muestral comprende solo los individuos
seleccionados para exportar. Una subpoblación es relevante para un evento cuando
forma parte de su ámbito de selección o aporta al menos un individuo a la
muestra. Estos ámbitos deben identificarse siempre de forma explícita.

- **REQ-EVO-01.** A partir de una referencia, el sistema debe permitir crear una
  población inicial con cambios.
- **REQ-EVO-02.** A partir de la población inicial, el sistema debe generar una
  población final simulando explícitamente, generación a generación, las fases,
  propiedades y parámetros evolutivos indicados en el fichero de configuración
  poblacional del escenario. No debe sustituir este proceso por un muestreo
  directo de la población final. Las generaciones deben ser discretas y no
  solapadas: la nueva generación sustituye a la anterior como población activa.
- **REQ-EVO-03.** El fichero de configuración poblacional debe poder representar,
  entre otros ejemplos, LD, equilibrio o desviación de Hardy-Weinberg, endogamia
  y tasa de mutación.
- **REQ-EVO-04.** El sistema debe permitir que un escenario evolutivo incluya
  otras propiedades o parámetros poblacionales además de los ejemplos
  confirmados.
- **REQ-EVO-05.** En un escenario poblacional, las restricciones alélicas
  fijadas para individuos concretos deben conservarse al intentar alcanzar una
  frecuencia alélica poblacional objetivo. La configuración debe rechazarse si
  el objetivo y las restricciones no pueden conciliarse dentro de la tolerancia
  admitida.
- **REQ-EVO-06.** Las fases evolutivas deben ejecutarse en el orden definido. La
  población de salida de cada fase debe ser la población de entrada de la fase
  siguiente.
- **REQ-EVO-07.** Entre generaciones y fases, el sistema debe conservar la
  continuidad de la herencia y aplicar la mutación y los demás procesos
  configurados para cada intervalo.
- **REQ-EVO-08.** La configuración debe poder representar, entre otros
  escenarios, una fase de panmixia de duración configurable seguida de un cuello
  de botella de magnitud y duración configurables. Este ejemplo no limita los
  modelos, el número de fases ni sus valores.
- **REQ-EVO-09.** Cada evento de muestreo debe poder seleccionar individuos por
  identificadores estables, por selección aleatoria o mediante una combinación
  de ambos métodos. Toda selección aleatoria debe ser reproducible y realizarse
  sin reemplazo.
- **REQ-EVO-10.** La ejecución debe conservar la población final y los
  individuos elegidos en cada muestreo solicitado, junto con la generación y la
  población de origen. No debe ser obligatorio conservar poblaciones
  intermedias completas.
- **REQ-EVO-11.** Seleccionar un individuo en un muestreo no debe eliminarlo,
  modificarlo ni excluirlo de la reproducción o de cualquier proceso posterior
  que le corresponda según el modelo. La observación debe conservar su estado en
  la generación muestreada mientras la población continúa evolucionando.
- **REQ-EVO-12.** Las decisiones aleatorias del muestreo deben estar aisladas de
  las decisiones aleatorias que hacen evolucionar la población, de modo que
  solicitar, cambiar o retirar un muestreo no altere la trayectoria evolutiva.
- **REQ-EVO-13.** Por defecto, cada muestreo debe conservar el estado genotípico
  faseado de los individuos seleccionados y, como mínimo, la generación, la fase
  evolutiva o el escenario poblacional, la subpoblación de cada individuo cuando
  corresponda, los identificadores de evento e individuo y la referencia
  estructural que corresponde a cada uno.
- **REQ-EVO-14.** Un muestreo solo debe producir BAM o CRAM y VCF de verdad si el
  evento los solicita expresamente. Cada evento debe constituir una unidad de
  exportación y un experimento tecnológico independiente.
- **REQ-EVO-15.** La creación, reproducción, muestreo y desaparición de
  individuos deben conservar sus identificadores estables y las relaciones con
  sus parentales cuando existan.
- **REQ-EVO-16.** N objetivo y Ne objetivo deben poder variar de manera independiente entre
  generaciones, fases o tramos para representar expansiones, reducciones y
  cuellos de botella. El sistema debe aplicar en cada intervalo los valores que
  resulten de la configuración, sin asumir que N y Ne objetivo sean iguales.
- **REQ-EVO-17.** Cada subpoblación debe poder tener sus propios N objetivo, Ne objetivo,
  modelo de apareamiento y demás parámetros evolutivos, variables por
  generación, fase o tramo. También debe calcular su Ne realizado.
- **REQ-EVO-18.** El sistema debe aplicar la migración direccional configurada
  por pares de origen y destino mediante proporciones o números exactos, según el
  único modo elegido para la fase. Cada migrante debe salir del conjunto
  reproductor de origen y entrar en el de destino antes de la selección y la
  reproducción, sin duplicarse ni migrar a dos destinos en el mismo paso.
- **REQ-EVO-19.** Un individuo migrante debe conservar su identificador y su
  linaje. El sistema debe conservar su ubicación inicial y final dentro de su
  única generación de existencia y en las observaciones registradas. No debe
  presentarlo como si persistiera en la generación posterior.
- **REQ-EVO-20.** El sistema debe ejecutar los eventos configurados de división,
  fusión, fundación y extinción en la generación y fase indicadas.
- **REQ-EVO-21.** Para cada evento que cambie la estructura de subpoblaciones, el sistema
  debe conservar las relaciones de origen y destino entre las subpoblaciones y la
  pertenencia de cada individuo antes y después del evento.
- **REQ-EVO-22.** La selección de fundadores y el reparto de individuos en
  una división deben usar el catálogo declarativo común de REQ-CONF-108. El
  reparto debe realizarse sin reemplazo, ser completo y asignar cada progenitor
  de la subpoblación origen a una sola hija. Los números exactos de las hijas
  deben sumar el número de progenitores del origen. Sus proporciones deben
  resolverse conjuntamente como probabilidades excluyentes y sumar 1 tras
  aplicar la política general de normalización.
  Ningún individuo puede pertenecer a dos hijas ni quedar sin asignar. Después
  del reparto, cada hija debe generar por separado su propia descendencia
  conforme a su N objetivo y sus parámetros efectivos. Cada hija debe poder
  configurar un N objetivo independiente del número de progenitores que recibió.
  La suma de los N objetivo de las hijas no tiene que coincidir con el N de la
  subpoblación origen; la división puede combinarse con una expansión o una
  contracción.
  Debe conservar para cada fundador su subpoblación de origen, la hija asignada,
  el modo o la regla que determinó el reparto y su linaje posterior.
  Antes de ejecutar una división, el sistema debe comprobar que
  existen al menos dos hijas, que los tamaños o proporciones y las reglas o
  identificadores permiten un reparto completo y exclusivo, y que cada hija
  dispone de progenitores compatibles y suficientes para su reproducción
  configurada. Un reparto incompleto, duplicado, ambiguo o no viable debe ser un
  error que detenga la simulación antes de crear descendientes.
- **REQ-EVO-23.** La simulación generación a generación debe aplicar
  explícitamente selección de viabilidad, migración y mutación cuando estén
  activadas. La deriva genética debe emerger del muestreo finito de parentales,
  gametos, copias y descendientes; no debe existir un segundo paso que modifique
  artificialmente las frecuencias alélicas bajo la etiqueta «deriva».
- **REQ-EVO-24.** La intensidad de la deriva debe resultar de N, del esquema que
  intenta aproximar Ne objetivo, de la contribución reproductiva realizada y de
  la estructura poblacional. El sistema debe registrar las frecuencias alélicas
  realizadas por generación para poder observarla.
- **REQ-EVO-25.** En la primera versión, la selección debe ser exclusivamente de
  viabilidad y usar un único locus estable para toda la ejecución. Los valores
  genéticos por genotipo o dosis y la función que transforma el valor del
  carácter en aptitud pueden variar por fase, intervalo y subpoblación.
- **REQ-EVO-26.** En cada generación, el sistema debe aplicar los procesos en el
  orden predeterminado de REQ-REP-24 o en otro orden válido que respete todas sus
  dependencias obligatorias. Debe rechazar antes de simular cualquier orden que
  las infrinja.
- **REQ-EVO-27.** La vía evolutiva debe simular explícitamente la recombinación al
  formar descendientes, usando la semilla de la ejecución y la tasa uniforme o
  el mapa completo y los parámetros efectivos generales o del grupo
  correspondiente.
- **REQ-EVO-28.** Un mapa de recombinación debe interpretarse en las coordenadas
  de su referencia estructural y no debe reinterpretarse sobre una referencia de
  mapeo distinta.
- **REQ-EVO-29.** La evolución debe conservar de generación en generación los
  haplotipos, su fase y el origen de los segmentos heredados o recombinados.
- **REQ-EVO-30.** Cuando la mutación de repeticiones en tándem esté activa, la
  evolución debe aplicar por defecto el modelo escalonado o el modelo alternativo
  configurado, con los parámetros efectivos de la fase, la subpoblación, la región y el
  locus.
- **REQ-EVO-31.** Las mutaciones internas de repeticiones en tándem solo deben
  aplicarse en la evolución cuando estén activadas, usando sus tipos, tasas y
  distribuciones efectivos.
- **REQ-EVO-32.** La población inicial debe contener al menos una subpoblación.
  Cada subpoblación debe poder elegir una distribución inicial distinta. La suma
  de sus tamaños debe coincidir con el N total de la población.
- **REQ-EVO-33.** El equilibrio de Hardy-Weinberg debe ser un modo opcional de
  inicialización por locus. Debe obtener las probabilidades de dosificación
  alélica a partir de las frecuencias alélicas y usar un modelo apropiado para
  todas las ploidías admitidas, no una fórmula limitada a diploides.
- **REQ-EVO-34.** Cuando se use Hardy-Weinberg, los genotipos deben muestrearse
  de forma reproducible. El modo probabilístico no debe forzar recuentos
  exactos. El modo controlado puede fijar recuentos alélicos, pero debe distribuir
  las copias entre individuos aleatoriamente condicionado a Hardy-Weinberg, sin
  forzar recuentos genotípicos exactos. El modelo no debe introducir doble
  reducción ni multivalentes.
- **REQ-EVO-35.** Para inicializar un locus fuera de Hardy-Weinberg, la
  configuración debe poder indicar explícitamente la frecuencia de cada genotipo
  o dosificación alélica válida para su ploidía. Las frecuencias deben someterse
  a la política general de suma, normalización y error de REQ-CONF-106.
- **REQ-EVO-36.** A partir de una distribución genotípica explícita, el sistema
  debe derivar las frecuencias alélicas. Si también se proporcionan frecuencias
  alélicas incompatibles, debe rechazar la configuración antes de simular.
- **REQ-EVO-37.** Para combinar LD inicial y desviación completa de
  Hardy-Weinberg, debe admitirse una distribución avanzada de genotipos
  multilocus con fase. Cada genotipo debe indicar las copias haplotípicas
  completas del individuo conforme a su ploidía y su frecuencia o recuento.
- **REQ-EVO-38.** Cada locus o bloque y subpoblación debe elegir entre el modo
  probabilístico predeterminado y el modo controlado explícito. No debe mezclar
  ambos modos dentro de la misma unidad, aunque loci, bloques o subpoblaciones
  distintos pueden elegir modos distintos. Debe registrarse lo solicitado y lo
  realizado. Un recuento exacto de genotipos, dosis o haplotipos debe ser
  explícito. En Hardy-Weinberg controlado solo se fijan los recuentos alélicos;
  la distribución genotípica sigue siendo aleatoria y condicionada al modelo.
- **REQ-EVO-39.** En el modo probabilístico, los recuentos realizados deben ser
  estocásticos. En el modo controlado, los recuentos genotípicos explícitos deben
  sumar el número de individuos aplicable, mientras que los recuentos alélicos o
  haplotípicos deben sumar C, el número de copias cromosómicas disponibles.
- **REQ-EVO-40.** Para representar LD inicial, cada bloque debe declarar
  explícitamente los haplotipos posibles y sus frecuencias conjuntas. Las
  frecuencias deben someterse a la política general de suma, normalización y
  error de REQ-CONF-106. La
  generación desde D, D', r² u otros estadísticos objetivo no pertenece a la
  primera versión. Las frecuencias haplotípicas deben ser la fuente de verdad y
  las frecuencias alélicas deben derivarse de ellas. Si también se aportan
  frecuencias alélicas, solo deben servir como validación y una discrepancia
  fuera de tolerancia debe ser un error. Por defecto, las frecuencias
  haplotípicas deben usarse como probabilidades. El modo controlado debe admitir
  recuentos exactos de haplotipos por subpoblación cuya suma sea C, el número de
  copias cromosómicas disponibles, y distribuirlos entre individuos de forma
  reproducible condicionada al modelo genotípico o a plantillas explícitas.
- **REQ-EVO-41.** Los bloques de LD de una misma subpoblación no deben
  solaparse. El sistema debe rechazar antes de simular dos bloques que compartan
  cualquier variante. Las variantes fuera de los bloques deben generarse de
  forma independiente respecto a las demás.
- **REQ-EVO-42.** Un bloque de LD puede contener cualquier variante heredable:
  SNP, *indels*, microsatélites u otras repeticiones en tándem, elementos
  transponibles y variantes estructurales. Cada variante debe conservar sus
  reglas de solapamiento, puntos de ruptura y compatibilidad estructural.
- **REQ-EVO-43.** Los bloques de LD deben definir solo la asociación
  haplotípica inicial. Durante la evolución y los cruces, la recombinación puede
  romperlos y crear haplotipos nuevos conforme al mapa o tasa y a las
  restricciones estructurales. El identificador del bloque inicial debe
  conservarse para seguimiento, pero no debe implicar indivisibilidad.
- **REQ-EVO-44.** En cada generación muestreada, el sistema debe calcular desde
  la verdad faseada y conservar las frecuencias de los haplotipos de cada bloque
  inicial. No debe calcular D, D' ni r² en la primera versión; esas medidas deben
  poder obtenerse externamente a partir de las frecuencias guardadas.
- **REQ-EVO-45.** La observación de frecuencias haplotípicas debe servir para
  seguir el efecto de la recombinación y otros procesos biológicos. No debe
  presentarse como evaluación de programas externos ni de *SNP callers*.
- **REQ-EVO-46.** La inicialización debe admitir un coeficiente de endogamia F
  como opción automática para aumentar la homocigosis esperada. Debe poder
  definirse para la población completa y sustituirse por subpoblación; no debe
  poder definirse por locus, región o bloque.
- **REQ-EVO-47.** F debe estar entre 0 y 1, inclusive. F=0 debe corresponder al
  modelo de Hardy-Weinberg. Al aumentar F debe crecer progresivamente la
  probabilidad de que las copias homólogas de un individuo porten el mismo
  alelo. F=1 debe producir homocigosis completa en cada locus. Los valores fuera
  del intervalo deben ser un error.
- **REQ-EVO-48.** El modelo automático de F debe admitirse para todas las
  ploidías mayores que 1 y respetar el significado funcional de REQ-EVO-47. La
  primera versión no debe introducir modelos distintos de F según una herencia
  polisómica o disómica. Esos modelos quedan fuera de alcance. Configurar F para
  ploidía 1 debe ser un error claro.
- **REQ-EVO-49.** Una distribución genotípica explícita por locus o bloque debe
  sustituir localmente al modelo general de F. Ambos métodos solo son
  incompatibles si intentan aplicarse al mismo locus o bloque. F debe seguir
  aplicándose al resto del genoma.
- **REQ-EVO-50.** F debe ser un objetivo esperado de la población o
  subpoblación, no un valor exacto por individuo. La homocigosis realizada debe
  poder variar estocásticamente entre individuos y regiones. Para fijar un
  individuo concreto deben usarse genotipos completos con fase.
- **REQ-EVO-51.** El modo automático de F solo debe ajustar las probabilidades
  genotípicas y el exceso esperado de homocigosis. No debe generar por sí solo
  tramos realistas de homocigosis. Los ROH deben surgir de pedigrí, de
  autofecundación o de evolución con recombinación.
- **REQ-EVO-52.** El sistema debe registrar el F objetivo general o por
  subpoblación, su procedencia, el valor efectivo y el F realizado por
  subpoblación y generación registrada. Debe estimar este último comparando la
  heterocigosidad observada con la esperada según las frecuencias alélicas. En
  poliploides, debe contar como heterocigoto todo individuo cuyas copias
  homólogas no sean todas del mismo alelo. Un valor realizado negativo debe
  registrarse como posible exceso de heterocigotos y no tratarse como error. La
  fórmula técnica exacta se elegirá durante la planificación sin cambiar este
  significado funcional.
- **REQ-EVO-53.** Todos los individuos de la generación inicial deben
  considerarse fundadores. Cada fundador debe tener un identificador estable y
  copias haplotípicas propias con identificadores estables, aunque comparta
  genotipo o secuencia con otros fundadores.
- **REQ-EVO-54.** La verdad y la trazabilidad de los individuos muestreados
  deben llegar al fundador y a la copia haplotípica fundadora concreta, no solo
  a una secuencia equivalente de otro individuo.
- **REQ-EVO-55.** Un individuo no debe aparecer más de una vez dentro del mismo
  evento de muestreo. Puede volver a ser seleccionado por otro evento
  independiente de la misma generación. Las generaciones distintas no comparten
  individuos porque el modelo es discreto y no solapado.
- **REQ-EVO-56.** Si un evento solicita más individuos que los elegibles, debe
  fallar con un error claro antes de generar salidas costosas. No debe reducir el
  tamaño solicitado ni muestrear con reemplazo.
- **REQ-EVO-57.** Un evento puede definir estratos con cantidades propias por
  subpoblación o grupo y con filtros por genotipo u otras reglas configuradas.
  Dentro de cada estrato, la selección aleatoria debe ser sin reemplazo y
  reproducible.
- **REQ-EVO-58.** Los estratos de un mismo evento deben ser no solapados. Si un
  individuo resulta elegible para más de uno, la configuración debe rechazarse
  como ambigua antes del muestreo. También deben validarse los elegibles y el
  tamaño solicitado dentro de cada estrato.
- **REQ-EVO-59.** Para cada identificador explícito, el sistema debe validar que
  exista, pertenezca a la generación del evento y cumpla el filtro o estrato
  aplicable. Los identificadores explícitos son selecciones obligatorias.
- **REQ-EVO-60.** Cuando un evento combine identificadores explícitos y
  selección aleatoria, debe incluir primero los explícitos, retirarlos del
  conjunto elegible y sortear sin reemplazo hasta completar el tamaño. Debe
  fallar si hay más explícitos que el tamaño o no quedan elegibles suficientes.
- **REQ-EVO-61.** Si un individuo aparece en varios eventos de una misma
  generación, debe conservar exactamente el mismo genotipo, fase y ascendencia
  biológica. Cada evento debe materializar de forma independiente su copia de la
  verdad y sus salidas opcionales; no debe deduplicarlas mediante una única
  copia compartida.
- **REQ-EVO-62.** Cada evento debe usar una subsemilla propia para simular de
  forma independiente sus lecturas, degradaciones y alineamientos. Dos eventos
  que incluyan el mismo individuo representan secuenciaciones técnicas distintas
  y deben relacionar el identificador biológico, el evento y la subsemilla.
- **REQ-EVO-63.** Cada evento debe tener un directorio y un manifiesto propios.
  El manifiesto raíz debe enlazar todos los eventos, sus configuraciones
  resueltas, referencias, subsemillas y resultados.
- **REQ-EVO-64.** Cada evento de muestreo debe guardar por separado un resumen
  poblacional y un resumen muestral de la generación observada. El primero debe
  describir la población completa y cada subpoblación relevante. El segundo debe
  calcularse únicamente con los individuos seleccionados para el evento. Esta
  obligación se aplica aunque el evento no solicite BAM, CRAM ni VCF.
- **REQ-EVO-65.** Ambos resúmenes deben incluir, como mínimo, frecuencias
  alélicas, heterocigosidad, F, N y Ne cuando cada medida sea aplicable. Cada
  valor debe indicar su ámbito, generación, población o subpoblación y número de
  individuos usados. El resumen muestral debe denominar N muestral al número de
  individuos seleccionados, sin confundirlo con N objetivo ni N realizado de la
  población. F
  muestral debe estimarse con las frecuencias y la heterocigosidad de la muestra.
  Ne muestral, cuando pueda calcularse, debe ser reproductivo y usar solo los
  progenitores que aportaron descendientes a la muestra y sus contribuciones
  observadas en ella. Toda medida no aplicable debe identificarse como tal y
  explicar el motivo; no debe copiarse un valor poblacional en el ámbito
  muestral.
- **REQ-EVO-66.** Una fusión debe reunir sin duplicación todos los progenitores
  de dos o más subpoblaciones de origen en una subpoblación nueva. Las
  subpoblaciones originales deben dejar de existir después del evento. El
  sistema debe conservar para cada individuo y haplotipo su subpoblación de
  procedencia y su linaje anterior y posterior a la fusión.
- **REQ-EVO-67.** La subpoblación fusionada debe poder configurar sus propios N objetivo,
  Ne objetivo, selección, tasa de autofecundación y modelo reproductivo. Por
  defecto, cualquier individuo elegible debe poder cruzarse con otro elegible
  de cualquier subpoblación de origen. Las restricciones de compatibilidad
  configuradas pueden limitar las parejas admisibles, pero no deben ponderar
  unas procedencias sobre otras. Por defecto, la contribución esperada de cada
  origen debe ser proporcional al número de progenitores elegibles que aportó.
  En este modo probabilístico, la configuración debe poder sustituir esa
  proporción mediante pesos por origen, interpretados como probabilidades. La
  contribución realizada puede diferir de la esperada. Como alternativa
  explícita, el modo controlado debe fijar recuentos exactos de contribuciones
  parentales por origen compatibles con la descendencia solicitada. Los modos
  no deben combinarse en una fusión.
  Su N objetivo debe ser independiente del número de progenitores reunidos y
  debe permitir expansión o contracción durante la primera generación.
  El sistema debe registrar el modo, los pesos o recuentos objetivo, las
  contribuciones parentales realizadas y la contribución genética realmente
  transmitida desde cada origen, conservando la trazabilidad hasta los
  haplotipos transmitidos. El modo controlado debe garantizar los recuentos
  parentales exactos, pero no un porcentaje exacto de bases ancestrales. El
  sistema debe calcular y registrar por separado la ascendencia genómica
  realizada por origen, que puede variar por la recombinación y por la
  ascendencia previa de los progenitores.
- **REQ-EVO-68.** Una fundación debe seleccionar progenitores de una
  subpoblación de origen mediante el catálogo declarativo común de
  REQ-CONF-108 y transferirlos sin reemplazo a una subpoblación nueva. En esa
  generación, cada fundador debe pertenecer y reproducirse solo en el conjunto
  reproductor nuevo y quedar excluido del conjunto de origen.
- **REQ-EVO-69.** Después de una fundación, la subpoblación nueva debe generar
  su descendencia desde los fundadores transferidos y la de origen debe generar
  la suya desde los individuos restantes. La subpoblación de origen debe
  continuar existiendo. El sistema debe conservar la procedencia, el modo o
  regla de selección y el linaje de cada fundador. El N objetivo de la nueva
  subpoblación debe ser independiente del número de fundadores y debe permitir
  expansión o contracción durante la primera generación.
- **REQ-EVO-70.** En una extinción programada, la subpoblación indicada no debe
  producir descendencia y debe dejar de existir desde esa generación. El
  sistema debe conservar su historial y todas las salidas creadas antes de la
  extinción. Cualquier referencia posterior que dependa de esa subpoblación y
  pueda detectarse de antemano debe ser un error de configuración; no debe
  tratarse como una operación vacía.
- **REQ-EVO-71.** Una recolonización posterior a una extinción debe crear una
  subpoblación nueva con un identificador distinto. Puede registrar que ocupa o
  recoloniza el ámbito de la subpoblación extinguida, pero no debe reutilizar su
  identificador ni reactivar su historial. Las dos trayectorias demográficas
  deben permanecer separadas y enlazadas solo mediante la relación explícita de
  recolonización.
- **REQ-EVO-72.** Si una extinción por resultado autorizada hace imposible un
  evento posterior que dependía de la subpoblación, el sistema debe omitir ese
  evento y registrarlo como «no ejecutado por extinción», con la subpoblación y
  la extinción causante. Esta regla se aplica a migraciones, cruces, fusiones,
  fundaciones y muestreos dependientes. La omisión no debe detener la simulación
  de los demás eventos y subpoblaciones. Los eventos deben ser atómicos frente a
  estas dependencias: si requieren varias subpoblaciones y falta al menos una,
  debe omitirse el evento completo, sin reinterpretarlo con las restantes. El
  registro debe conservar el evento, todas las subpoblaciones requeridas y
  cuáles faltaban. La omisión debe propagarse transitivamente a todo evento que
  dependa directa o indirectamente de la subpoblación extinguida o de un
  resultado que no llegó a crearse. Cada omisión debe conservar una cadena de
  causas auditable hasta la extinción original. Las ramas sin esa dependencia
  deben continuar.
- **REQ-EVO-73.** Los eventos demográficos de una misma generación deben
  ejecutarse en el orden explícito configurado. Cada evento debe operar sobre el
  estado completo producido por el anterior. Antes de simular, el sistema debe
  validar sus dependencias y rechazar un orden ausente, ambiguo o imposible,
  identificando los eventos y la dependencia incumplida.
- **REQ-POB-01.** El sistema debe permitir crear poblaciones naturales mediante
  esta vía.

#### 4.3.4. Vía 3: familias y poblaciones segregantes

Una población segregante es una descendencia en la que los cambios genéticos de
los parentales se distribuyen entre los individuos. Un pedigrí describe las
relaciones de parentesco y descendencia. Una población F2 es la segunda
generación filial obtenida al cruzar individuos de la generación F1. Un
retrocruzamiento (*backcross*) cruza un descendiente con uno de sus parentales o
con un individuo genéticamente equivalente.

Las líneas de introgresión incorporan fragmentos de un parental donante en el
fondo genético de otro parental. Las poblaciones de cruce intergeneracional
avanzado con múltiples parentales (*Multiparent Advanced Generation Inter-Cross*,
MAGIC) combinan material genético de varios fundadores mediante cruces
sucesivos.

Estos tipos son ejemplos, no un catálogo cerrado. El sistema podrá ofrecer tipos
de cruce predefinidos y también deberá admitir cruces más complejos definidos por
la persona usuaria.

Un grafo dirigido representa relaciones que avanzan desde una fuente hacia un
resultado. Es acíclico cuando ninguna cadena de relaciones puede volver al punto
del que partió. En un grafo de cruces, cada evento recibe parentales o fuentes y
produce descendientes que pueden alimentar eventos posteriores.

En este contexto, un ciclo es una dependencia imposible: un evento o individuo
necesitaría como antecedente un resultado que depende directa o indirectamente
de él. Una autofecundación no crea una relación del individuo consigo mismo:
usa un parental existente y produce descendientes nuevos. Un retrocruzamiento
puede reutilizar un parental anterior junto con un descendiente posterior y
también produce descendientes nuevos. Ninguno de estos casos constituye un ciclo.

En esta vía, la verdad biológica es el estado genotípico y haplotípico simulado
de cada descendiente antes de generar lecturas o alineamientos. Las reglas para
elegir parentales posteriores se evalúan sobre esa verdad, no sobre evidencias
técnicas que puedan estar incompletas o degradadas.

- **REQ-CRU-01.** El sistema debe permitir crear familias o poblaciones
  segregantes mediante eventos con un parental para autofecundación o dos para
  cruzamiento sexual. Los diseños pueden incorporar más fundadores mediante una
  secuencia de eventos.
- **REQ-CRU-02.** Esta vía debe admitir, entre otros ejemplos, pedigríes,
  poblaciones F2, retrocruzamientos, líneas de introgresión y poblaciones MAGIC.
- **REQ-CRU-03.** La primera versión debe ofrecer plantillas predefinidas para
  F1, F2, retrocruce, autofecundación y pedigrí simple.
- **REQ-CRU-04.** Cuando el cruce solicitado no sea uno de los predefinidos, la
  persona usuaria debe poder definirlo en el fichero de configuración
  poblacional del escenario.
- **REQ-CRU-05.** El fichero de configuración poblacional debe indicar los
  parentales, los cruces, su orden o secuencia, las relaciones necesarias, las
  generaciones y los demás parámetros del escenario.
- **REQ-CRU-06.** La vía de familias y poblaciones segregantes debe simular la
  recombinación al formar cada descendencia mediante la tasa uniforme o los
  mapas generales o específicos de grupo configurados para las referencias
  estructurales correspondientes.
- **REQ-CRU-07.** Cada parental debe proceder de una de estas dos fuentes: la
  referencia mediante su identificador estable a un individuo ya generado y
  disponible en la simulación, o una definición expresa dentro del fichero de
  configuración poblacional.
- **REQ-CRU-08.** Un parental definido en el fichero poblacional debe seguir las
  reglas comunes de los individuos e indicar su referencia estructural, usar la
  ploidía efectiva por secuencia y definir sus cambios por copia haplotípica
  respecto a una referencia compatible. Los loci y regiones no mencionados deben
  conservar el alelo o la secuencia de esa referencia.
- **REQ-CRU-09.** Antes de simular, el sistema debe comprobar que los
  identificadores de parentales referenciados existen, que los parentales
  definidos tienen identificadores únicos y que todos son compatibles con el
  escenario. Si alguna comprobación falla, debe rechazar la configuración e
  identificar el parental y el problema.
- **REQ-CRU-10.** El sistema debe ejecutar los eventos del grafo de cruces en un
  orden compatible con sus generaciones y dependencias.
- **REQ-CRU-11.** La descendencia producida por un evento debe poder utilizarse
  como fuente o parental de uno o varios eventos posteriores, incluidas ramas
  distintas del grafo.
- **REQ-CRU-12.** Cada evento debe producir el número de descendientes indicado y
  asignarles identificadores estables de individuo, conservando su relación con
  el evento y los parentales de origen.
- **REQ-CRU-13.** Antes de ejecutar una selección de parentales posteriores, el
  sistema debe comprobar que existen la población candidata, los identificadores
  y los datos requeridos por la regla. Cuando se solicite un número exacto, si
  hay menos descendientes elegibles debe detener la ejecución e identificar la
  población, la regla, la cantidad solicitada y la cantidad disponible. No debe
  reducir la cantidad ni relajar los criterios de forma silenciosa.
- **REQ-CRU-14.** Una selección aleatoria de descendientes debe usar las
  decisiones aleatorias reproducibles de la simulación. Si una regla produce
  más descendientes elegibles que el número exacto solicitado, debe elegirlos al
  azar sin reemplazo usando la semilla y registrar los candidatos y
  seleccionados.
- **REQ-CRU-15.** El sistema debe conservar obligatoriamente el estado genotípico
  faseado y los metadatos de todos los individuos de todos los eventos y
  poblaciones intermedias, aunque no sean seleccionados ni vuelvan a utilizarse.
- **REQ-CRU-16.** Una población intermedia solo debe producir BAM o CRAM y VCF de
  verdad si su configuración los solicita expresamente. Cuando se soliciten,
  deben incluir a todos los individuos de esa población.
- **REQ-CRU-17.** Cada cruce debe construir los haplotipos de la descendencia a
  partir de los parentales y conservar la fase y el origen de cada segmento a
  través de los eventos posteriores del grafo.
- **REQ-CRU-18.** Cuando el escenario permita mutaciones durante la transmisión,
  la vía de familias debe aplicar a las repeticiones en tándem el modelo y los
  parámetros efectivos configurados.
- **REQ-CRU-19.** Durante la transmisión familiar, las mutaciones internas de
  repeticiones en tándem solo deben aplicarse cuando estén activadas.
- **REQ-CRU-20.** Las líneas endogámicas recombinantes (RIL), las líneas de
  introgresión y las poblaciones MAGIC no deben tener una plantilla predefinida
  en la primera versión, pero deben poder expresarse mediante el sistema general
  de cruces personalizados. Podrán añadirse como plantillas futuras sin cambiar
  ese modelo general.
- **REQ-CRU-21.** La plantilla F2 debe permitir elegir de forma excluyente entre
  autofecundar individuos F1 o cruzar dos individuos F1. El pedigrí y el
  manifiesto deben registrar el modo, los parentales efectivos y los
  descendientes.
- **REQ-CRU-22.** El sistema debe ser aplicable a especies o sistemas autógamos,
  que se reproducen por autofecundación, y alógamos, que cruzan individuos
  distintos, de forma coherente con los cruces personalizados y el modelo
  poblacional.
- **REQ-CRU-23.** Toda inserción de elemento transponible presente en un
  haplotipo parental debe heredarse según la segregación, la ploidía y la
  recombinación, conservando su fase, identificador de evento o copia y linaje.
  No debe tratarse como una anotación que se pierde al cruzar.
- **REQ-CRU-24.** Los cruces, la autofecundación y la evolución con reproducción
  sexual genérica solo deben admitir cromosomas con ploidía par en la primera
  versión. La configuración debe rechazarse antes de simular si cualquiera de
  los cromosomas reproductivos tiene ploidía impar.
- **REQ-CRU-25.** Cada parental debe aportar la mitad de las copias de cada
  cromosoma al descendiente, de modo que se mantenga su ploidía par.
- **REQ-CRU-26.** Los parentales de un cruzamiento deben ser simétricos: cambiar
  su orden no debe cambiar el modelo biológico. La primera versión no debe asignar
  sexos ni modelar cruces recíprocos, herencia citoplasmática o efectos maternos
  o paternos.
- **REQ-CRU-27.** Cada evento reproductivo debe admitir como máximo dos
  parentales. Un parental debe representar autofecundación y dos parentales un
  cruzamiento sexual. Cualquier otra cardinalidad debe rechazarse.
- **REQ-CRU-28.** Un diseño con tres o más fundadores debe expresarse como una
  sucesión de eventos de uno o dos parentales en el grafo dirigido acíclico. Los
  individuos y poblaciones intermedios deben conservarse con la trazabilidad ya
  acordada.
- **REQ-CRU-29.** La primera versión solo debe admitir autofecundación y
  cruzamiento sexual. Debe rechazar reproducción clonal, apomixis y producción de
  dobles haploides como modos reproductivos.
- **REQ-CRU-30.** Las translocaciones pueden estar presentes en referencias,
  individuos y lecturas. Antes de reproducir, las copias homólogas deben ser
  estructuralmente compatibles. Un individuo heterocigoto para una translocación
  no debe poder reproducirse en la primera versión y debe producir un error de
  validación claro.
- **REQ-CRU-31.** Las mutaciones de novo deben estar disponibles en familias y
  poblaciones segregantes mediante el mismo motor usado en evolución, pero
  desactivadas de forma predeterminada. Su activación debe ser explícita y debe
  registrar tasas, eventos realizados y trazabilidad con las mismas reglas.
- **REQ-CRU-32.** El catálogo funcional de selección de descendientes de la
  primera versión debe ser el catálogo declarativo común de REQ-CONF-108.
- **REQ-CRU-33.** Las condiciones de selección deben poder combinarse mediante
  Y, O y NO. La cantidad solicitada debe expresarse como un número exacto o como
  una proporción de los descendientes elegibles.
- **REQ-CRU-34.** Una proporción debe interpretarse como la probabilidad
  individual de seleccionar cada descendiente elegible. El recuento realizado
  debe ser aleatorio y reproducible mediante la semilla. Deben registrarse la
  proporción solicitada, los elegibles, los seleccionados y la proporción
  realizada. Cuando se necesite controlar estrictamente el recuento, la
  configuración debe exigir un número exacto.
- **REQ-CRU-35.** Todas las condiciones deben evaluarse sobre la verdad
  biológica de los descendientes. No deben usar BAM, cobertura, calidad de base
  o de mapeo ni ninguna otra evidencia técnica generada después.
- **REQ-CRU-36.** La insuficiencia de elegibles para un número exacto en un
  cruce planificado debe ser un error estricto. Esta regla no debe confundirse
  con la evolución poblacional, donde un N realizado menor que N objetivo puede
  conservarse con una advertencia cuando lo cause un resultado aleatorio
  admisible.
- **REQ-CRU-37.** Toda variante y todo cambio estructural de un parental definido
  expresamente deben asignarse a una o varias copias haplotípicas concretas. La
  asignación debe ser suficiente para reconstruir la secuencia completa y la
  fase de cada copia a partir de la referencia compatible.
- **REQ-CRU-38.** La primera versión no debe inventar ni sortear la fase de un
  parental definido expresamente. Si la entrada no permite reconstruir sin
  ambigüedad todas sus copias haplotípicas, el sistema debe rechazarla antes de
  simular e identificar el parental y la información que falta.
- **REQ-CRU-39.** La compatibilidad de la referencia de un parental expreso debe
  validarse antes de usarlo en un cruce, incluidas las secuencias necesarias y
  su ploidía efectiva. La omisión de un locus o una región no debe interpretarse
  como dato desconocido ni causar una imputación: debe conservar el estado de
  referencia.
- **REQ-CRU-40.** Reutilizar un parental o un descendiente intermedio en varios
  cruces no debe duplicar su identidad biológica ni alterar su genotipo faseado.
  Cada evento posterior debe formar gametos de manera independiente y producir
  descendientes nuevos con identificadores propios. El grafo debe conservar la
  relación de todos esos eventos con el mismo parental, permitiendo diseños
  ramificados como el uso de una misma F1 en varios retrocruzamientos.
- **REQ-CRU-41.** La validación de ciclos debe analizar dependencias entre
  fuentes, eventos y resultados, no la repetición de un identificador parental
  en eventos distintos. Debe admitir la autofecundación y el retrocruzamiento
  cuando todos sus parentales existen antes del evento y sus descendientes son
  individuos nuevos. Debe rechazar una dependencia en la que un evento o
  individuo requiera, directa o transitivamente, un resultado que depende de él.
- **REQ-CRU-42.** Antes de simular el grafo, el sistema debe comprobar que los
  identificadores de individuos, poblaciones y eventos son únicos dentro de su
  categoría. Toda referencia a uno de ellos debe existir y apuntar a un resultado
  creado antes según el orden de dependencias.
- **REQ-CRU-43.** Cada descendiente debe ser creado exactamente una vez por un
  único evento y recibir un identificador nuevo. Esta regla no impide que un
  individuo ya creado se use después como parental en uno o varios eventos.
- **REQ-CRU-44.** Cada autofecundación debe tener exactamente un parental y cada
  cruzamiento sexual exactamente dos. Antes de ejecutar el evento, el sistema
  debe comprobar la compatibilidad de sus parentales en referencia, ploidía y
  estructura conforme a las reglas reproductivas de la primera versión.
- **REQ-CRU-45.** Toda regla de selección y toda salida solicitada en el grafo
  debe referirse a resultados que puedan producirse según sus dependencias,
  cardinalidades y datos disponibles. Una regla o salida imposible debe
  rechazarse antes de ejecutar el evento e identificar la referencia, el
  requisito incumplido y los resultados afectados.
- **REQ-POB-02.** El sistema debe permitir crear poblaciones segregantes.
- **REQ-POB-03.** El sistema debe permitir crear familias.

#### 4.3.5. Reproducción y dinámica generacional

- **REQ-REP-01.** El sistema de apareamiento debe admitir reproducción mixta
  mediante una tasa de autofecundación entre 0 % y 100 %, inclusive.
- **REQ-REP-02.** La tasa de autofecundación debe poder definirse por población,
  subpoblación, fase o intervalo de generaciones. El valor del ámbito más específico
  debe sustituir al más general durante su intervalo.
- **REQ-REP-03.** La tasa de autofecundación debe interpretarse como la
  probabilidad aplicada a cada evento reproductivo, no como una cuota exacta de
  eventos por generación. Las decisiones deben ser aleatorias y reproducibles.
- **REQ-REP-04.** Los eventos que no resulten en autofecundación deben ser
  alógamos. En la primera versión deben elegir de forma aleatoria uniforme dos
  individuos distintos entre los supervivientes elegibles de la población o
  subpoblación y generación.
- **REQ-REP-05.** La elegibilidad puede depender de selección de viabilidad,
  pertenencia a la población o disponibilidad. La primera versión no debe
  ponderar la elección por parentesco, distancia genética, genotipo ni
  apareamiento selectivo. Se exceptúan los pesos por subpoblación de origen
  configurados para una fusión conforme a REQ-EVO-67.
- **REQ-REP-06.** Un individuo debe poder participar en varios eventos
  reproductivos y producir varios descendientes en una generación.
- **REQ-REP-07.** La configuración debe permitir definir límites o una
  distribución del número de descendientes por parental. Estos parámetros deben
  coordinarse con N y con Ne cuando exista un objetivo para este último.
- **REQ-REP-08.** N debe ser el objetivo configurado del tamaño de cada
  nueva generación y subpoblación. En una transición reproductiva, el N efectivo
  debe fijar cuántos descendientes se intentarán crear para la generación
  siguiente. Un resultado aleatorio puede impedir alcanzarlo conforme a
  REQ-REP-26.
- **REQ-REP-09.** N debe ser obligatorio al crear cada población o subpoblación inicial.
  Si se omite en una generación o fase posterior, debe heredarse el último valor
  efectivo.
- **REQ-REP-10.** N solo debe cambiar mediante un nuevo valor explícito o un
  evento demográfico configurado, como cuello de botella, expansión, fundación o
  fusión. Un cuello de botella debe reducir el N objetivo de la generación
  siguiente y una expansión debe aumentarlo. El historial debe
  distinguir valores heredados y cambios, así como la generación en la que cada
  valor determina la descendencia.
- **REQ-REP-11.** Para cada fase y subpoblación, Ne puede configurarse como valor
  objetivo absoluto o como proporción Ne/N. Ambas formas son mutuamente
  excluyentes. Aunque no se configure un objetivo, el sistema debe representar
  este campo como «sin objetivo» para distinguirlo del Ne realizado.
- **REQ-REP-12.** La proporción Ne/N debe ser válida y convertirse en un Ne
  objetivo usando el N efectivo de la fase o intervalo.
- **REQ-REP-13.** Si Ne se omite, el sistema no debe intentar alcanzar un Ne
  objetivo ni advertir por desviación. Debe calcular el Ne realizado que emerja
  de la reproducción.
- **REQ-REP-14.** Cuando exista, Ne objetivo debe guiar el esquema reproductivo
  para aproximar la desigualdad de contribución necesaria. El sistema debe
  calcular Ne realizado exclusivamente como Ne reproductivo, a partir de los
  progenitores que dejaron descendencia y de la desigualdad o varianza de sus
  contribuciones. No debe usar estimadores basados en pérdida de
  heterocigosidad, incremento de endogamia o desequilibrio de ligamiento. La
  fórmula concreta y su adaptación por ploidía se elegirán durante la
  planificación sin cambiar este significado funcional.
- **REQ-REP-15.** El sistema debe cuantificar y advertir cualquier desviación
  entre Ne objetivo y Ne realizado. Esta desviación no debe invalidar el
  resultado, detener la simulación ni poder elevarse a error. Una configuración
  estructuralmente imposible o incoherente entre N, Ne objetivo y los límites o
  distribuciones de descendencia sí debe rechazarse antes de simular.
- **REQ-REP-16.** El historial y el manifiesto deben registrar por generación y
  subpoblación los parentales, descendientes, distribución familiar, tasa de
  autofecundación objetivo, número esperado y realizado de autofecundaciones, N,
  forma original de Ne, Ne objetivo o su ausencia, Ne realizado y la desviación
  cuantificada con su advertencia cuando exista un objetivo.
- **REQ-REP-17.** Cada generación reproductiva debe crear la siguiente, que
  sustituye a la anterior como población activa. No deben modelarse edades,
  supervivencia individual entre generaciones ni reproducción entre cohortes.
  Un cambio de N no debe eliminar ni añadir retrospectivamente individuos en la
  generación parental ya creada. Sus individuos solo deben actuar como posibles
  progenitores hasta que la descendencia los sustituya.
- **REQ-REP-18.** Los individuos históricos deben conservar identificadores no
  reutilizables y la trazabilidad exigida aunque ya no pertenezcan a la población
  activa.
- **REQ-REP-19.** Deben conservarse por defecto los genotipos completos de la
  población inicial, la final y las generaciones o muestreos solicitados
  expresamente, no los de todas las generaciones intermedias.
- **REQ-REP-20.** Toda generación debe producir un historial estructurado
  resumido de demografía, procesos, N objetivo, N realizado, Ne objetivo, Ne realizado,
  reproducción, eventos, frecuencias alélicas y muestreos.
- **REQ-REP-21.** La configuración debe permitir marcar una generación adicional
  para conservar su estado genotípico completo y solicitar opcionalmente sus BAM
  o CRAM y VCF conforme a las reglas ya acordadas.
- **REQ-REP-22.** La regeneración exacta de decisiones debe requerir la misma
  semilla, configuraciones efectivas, entradas y sus comprobaciones, perfiles,
  versión del simulador y dependencias relevantes. El manifiesto debe conservar
  o identificar estos elementos.
- **REQ-REP-23.** No debe garantizarse reproducción idéntica entre versiones que
  cambien algoritmos o flujos aleatorios, salvo que una versión declare y cumpla
  expresamente una compatibilidad concreta.
- **REQ-REP-24.** El orden predeterminado de cada generación debe ser: eventos
  demográficos y determinación del N de la descendencia siguiente; migración;
  selección de viabilidad;
  elección de parentales y reproducción; segregación y recombinación; mutaciones
  de novo en los haplotipos descendientes; cálculo y registro de N realizado,
  Ne realizado y frecuencias; comparación con N objetivo y Ne objetivo; y muestreo y salidas
  solicitadas. La deriva debe emerger durante la reproducción.
  Cuando haya varios eventos demográficos, su orden interno no debe ser
  predeterminado: debe proceder de la secuencia explícita exigida por
  REQ-CONF-105.
- **REQ-REP-25.** La configuración solo debe poder reordenar procesos opcionales
  si mantiene las dependencias obligatorias. En particular, la mutación de novo
  debe ocurrir después de la segregación y recombinación y antes del registro y
  muestreo de la nueva generación.
- **REQ-REP-26.** Si la migración, la selección u otro resultado aleatorio de la
  ejecución impiden alcanzar N objetivo, el sistema debe conservar la
  subpoblación con el N realizado que haya podido generar. Debe registrar N
  objetivo, N realizado, su desviación y una advertencia. Esta desviación no
  debe detener la simulación ni declarar automáticamente una extinción. Una
  configuración estructuralmente imposible que pueda detectarse antes de
  ejecutar sí debe rechazarse. Si N realizado es 0, debe aplicarse REQ-SEL-09.

#### 4.3.6. Selección de viabilidad

Un carácter cuantitativo es una característica expresada mediante un valor
numérico. En esta versión puede depender del genotipo de un único locus, incluso
si ese locus tiene varios alelos; no convierte el modelo en selección
poligénica. Una función direccional asigna mayor aptitud de forma sistemática
hacia uno de los dos sentidos del carácter: valores mayores o valores menores.

- **REQ-SEL-01.** La primera versión debe usar un único locus seleccionado para
  toda la ejecución. Configurar loci distintos entre fases, intervalos,
  poblaciones o subpoblaciones debe ser un error.
- **REQ-SEL-02.** El locus seleccionado puede ser un SNP, *indel*,
  microsatélite o repetición en tándem, inserción de TE o variante estructural,
  siempre que tenga identificador estable y genotipo o dosis inequívocos. Puede
  ser bialélico o multialélico.
- **REQ-SEL-03.** La configuración debe referenciar el locus por su identificador
  y referencia estructural, no solo mediante una coordenada ambigua, y debe
  admitir ploidías de 1 a 10 en los genotipos que describe.
- **REQ-SEL-04.** Debe asignarse un valor genético numérico a cada genotipo
  completo o combinación de dosis pertinente, incluidos los genotipos
  multialélicos. El valor genético representa la contribución del genotipo y no
  debe confundirse con la aptitud.
- **REQ-SEL-05.** Un conjunto incompleto de valores genéticos para los genotipos
  o combinaciones de dosis debe aportar una regla predeterminada explícita. Si no
  la aporta, el sistema debe rechazar la configuración. La representación
  concreta de esa regla se decidirá con el esquema de configuración.
- **REQ-SEL-06.** Los valores genéticos y la función de aptitud del mismo locus
  pueden variar por fase, intervalo y subpoblación. El valor más específico debe
  prevalecer y el locus no debe cambiar.
- **REQ-SEL-07.** La supervivencia debe decidirse de forma estocástica y
  reproducible para cada individuo. La reproducción posterior debe elegir
  uniformemente entre los supervivientes elegibles.
- **REQ-SEL-08.** Deben registrarse los valores genéticos y del carácter, las
  aptitudes antes y después de normalizar, las decisiones de supervivencia y los
  supervivientes por genotipo y generación.
- **REQ-SEL-09.** Si quedan cero individuos o menos parentales elegibles de los
  necesarios, el sistema no debe declarar automáticamente una extinción ni
  detenerse por una desviación positiva de N realizado respecto a cero. Debe
  aplicar REQ-REP-26. Si N realizado es 0, nunca debe conservar una
  subpoblación activa vacía. Si la configuración permite la extinción por
  resultado en la subpoblación y fase aplicables, debe declararla y registrarla
  como distinta de una extinción planificada. Si no la permite, debe detener
  inmediatamente toda la simulación con un error de ejecución. En ese caso no
  debe omitirse ningún evento ni continuar ninguna otra rama.
- **REQ-SEL-10.** Tras una extinción configurada explícitamente, N deja de ser
  aplicable. Una referencia posterior incompatible y detectable de antemano
  debe ser un error. Tras una
  extinción por resultado autorizada, los eventos dependientes deben omitirse
  conforme a REQ-EVO-72. Una fundación o recolonización posterior válida debe
  crear otra subpoblación con un identificador nuevo; no debe reactivar la
  extinguida.
- **REQ-SEL-11.** El único locus seleccionado debe poder representar un carácter
  cuantitativo mediante los valores genéticos asociados a sus genotipos o
  combinaciones de dosis. El límite de un único locus debe mantenerse.
- **REQ-SEL-12.** El sistema debe mantener por separado el valor genético, el
  valor del carácter y la aptitud relativa de cada individuo. Una función
  configurable debe transformar el valor del carácter en una aptitud relativa
  no negativa.
- **REQ-SEL-13.** Las aptitudes relativas producidas por la función deben
  normalizarse respecto al máximo, cuya probabilidad de supervivencia será 1;
  las demás probabilidades deben ser proporcionales.
- **REQ-SEL-14.** En la primera versión, el valor del carácter debe coincidir
  exactamente con el valor genético asignado al genotipo o combinación de dosis.
  La configuración no debe admitir ruido ambiental, efectos ambientales ni un
  parámetro de heredabilidad.
- **REQ-SEL-15.** La función de aptitud cuantitativa debe ser direccional y la
  configuración debe elegir exactamente un sentido: favorecer valores mayores o
  favorecer valores menores del carácter. Al favorecer valores mayores, un valor
  mayor no debe recibir menor aptitud que uno inferior. Al favorecer valores
  menores, un valor menor no debe recibir menor aptitud que uno superior.
- **REQ-SEL-16.** La primera versión debe rechazar funciones estabilizadoras, con
  un óptimo intermedio favorecido, y funciones disruptivas, que favorezcan ambos
  extremos frente a valores intermedios. No debe aproximarlas silenciosamente
  mediante una función direccional.

#### 4.3.7. Migración

- **REQ-MIG-01.** La migración debe configurarse mediante pares dirigidos de
  subpoblación origen y destino. En cada fase o intervalo, todos los pares deben
  usar exclusivamente proporciones o exclusivamente números exactos de
  migrantes. Mezclar ambos modos debe ser un error.
- **REQ-MIG-02.** En el modo proporcional, los migrantes deben seleccionarse al
  azar y de forma reproducible. La proporción de cada par debe interpretarse
  como la probabilidad de que cada individuo elegible migre por esa ruta, no
  como una cuota que deba cumplirse exactamente. El número realizado puede
  variar entre ejecuciones con semillas distintas. Para exigir una cifra precisa
  debe usarse el modo exacto, que selecciona el número entero configurado. La
  selección aleatoria debe ser la predeterminada, pero cada ruta debe poder usar
  cualquier regla del catálogo declarativo común de REQ-CONF-108.
  Todas las rutas proporcionales salientes de un mismo origen deben resolverse
  conjuntamente como categorías excluyentes: cada tasa es la probabilidad de ir
  a su destino y la probabilidad de permanecer es 1 menos la suma de las tasas
  salientes. Las rutas exactas también deben asignarse conjuntamente y sin
  reemplazo. Un individuo no debe poder migrar a dos destinos en el mismo paso.
- **REQ-MIG-03.** Antes de simular deben validarse los identificadores de origen
  y destino, el modo de la fase, las proporciones o cantidades y la
  disponibilidad de individuos. La suma de las proporciones salientes de cada
  origen debe estar entre 0 y 1, inclusive; una suma superior a 1 debe invalidar
  la configuración. Las cantidades deben ser enteros no negativos y su suma
  saliente no debe superar los individuos disponibles. Las
  reglas e identificadores deben producir conjuntos compatibles y suficientes
  para la asignación conjunta. Si faltan individuos elegibles, el sistema debe
  detenerse con un error claro; no debe reducir cantidades, reutilizar
  individuos ni omitir una ruta.
- **REQ-MIG-04.** La migración debe transferir individuos completos: el migrante
  sale del conjunto reproductor de origen y entra en el de destino, conservando
  identidad y linaje. La transferencia debe completarse antes de aplicar la
  selección de viabilidad y la reproducción.
- **REQ-MIG-05.** Salvo un evento demográfico explícito, la migración debe
  conservar N total, aunque los tamaños de las subpoblaciones puedan desviarse
  temporalmente. Estos tamaños son los de los conjuntos reproductores de esa
  generación y no deben modificar automáticamente el N objetivo de la
  descendencia. Cada subpoblación debe crear la generación siguiente con su N
  configurado. Un cambio demográfico permanente debe exigir una modificación
  explícita de N.
- **REQ-MIG-06.** Deben registrarse los tamaños antes y después de la migración,
  el modo, los valores configurados y realizados por par, las rutas individuales
  y el N objetivo de la descendencia de cada subpoblación. En el modo
  proporcional deben registrarse la tasa solicitada, el número realizado y la
  proporción realizada. En el modo exacto deben registrarse el número solicitado
  y el realizado.
- **REQ-MIG-07.** Un destino extinguido solo puede recibir migrantes si existe un
  evento válido de fundación o recolonización. Una ruta incompatible con una
  extinción planificada debe rechazarse antes de simular. Si la extinción fue un
  resultado aleatorio autorizado, la ruta dependiente debe registrarse como no
  ejecutada y aplicarse REQ-EVO-72.
- **REQ-MIG-08.** Todas las migraciones de una generación deben calcularse
  simultáneamente desde la ubicación inicial de cada individuo. Un migrante solo
  puede completar una ruta en su única generación de existencia; no debe
  encadenar movimientos como A→B→C. En la generación siguiente ese individuo ya
  no existe como población activa. Solo sus descendientes pueden ser candidatos
  a las migraciones configuradas para la nueva generación.

#### 4.3.8. Mutaciones evolutivas

- **REQ-MUT-01.** Las tasas evolutivas de SNP, inserciones, deleciones y
  variantes estructurales deben expresarse por base y generación, por separado,
  con valores globales y sustituciones por cromosoma o *scaffold* y región.
- **REQ-MUT-02.** Cada mutación nueva debe aplicarse a una copia o haplotipo
  concreto ya segregado y recombinado de un individuo de la nueva generación.
  No debe modificar retrospectivamente al parental.
- **REQ-MUT-03.** Cada evento debe conservar generación, individuo, referencia,
  posición, haplotipo, tipo, origen de novo y transmisión posterior.
- **REQ-MUT-04.** El modelo de SNP debe permitir un espectro o matriz válida con
  las probabilidades de cada sustitución, incluida la relación entre
  transiciones y transversiones.
- **REQ-MUT-05.** Si no se proporciona un espectro, debe usarse un modelo
  simétrico documentado entre las tres bases alternativas posibles.
- **REQ-MUT-06.** Las probabilidades del espectro deben validarse, admitir
  sustituciones regionales y registrarse con los eventos esperados y realizados.
  Cada conjunto que deba sumar 1 debe aplicar la política general de
  normalización y error de REQ-CONF-106.
- **REQ-MUT-07.** La primera versión no debe modelar tasas de SNP dependientes de
  CpG, trinucleótidos u otros contextos de bases vecinas.
- **REQ-MUT-08.** Inserciones y deleciones deben tener tasas independientes y
  distribuciones configurables de longitud. El perfil predeterminado debe
  favorecer *indels* cortos; su distribución cuantitativa se decidirá y validará
  durante la planificación.
- **REQ-MUT-09.** Los eventos menores que el umbral efectivo deben clasificarse
  como *indels* y los que lo alcancen o superen como variantes estructurales,
  manteniendo la regla predeterminada de 50 pb.
- **REQ-MUT-10.** Una inserción puede usar una secuencia nueva generada según una
  composición configurada, una secuencia concreta, una copia de una región del
  mismo genoma o una secuencia repetida, incluido un elemento transponible.
- **REQ-MUT-11.** Las repeticiones en tándem deben formar una clase mutacional
  independiente, con tasa por locus y generación o reglas globales o regionales,
  junto con expansión o contracción, sesgo, tamaño de paso, pureza e
  interrupciones opcionales.
- **REQ-MUT-12.** El historial, el manifiesto y la verdad deben conservar por
  clase y ámbito las tasas, el número esperado y realizado de eventos y sus
  transmisiones posteriores.

#### 4.3.9. Elementos transponibles configurados

- **REQ-TE-01.** El catálogo de TE debe poder suministrarse mediante FASTA o ser
  generado por el simulador. Todo catálogo debe relacionar familias o tipos con
  sus secuencias.
- **REQ-TE-02.** Un FASTA externo de TE debe validarse con reglas de secuencia
  compatibles, copiarse a resultados, usarse desde la copia y conservar intacto
  el original. El catálogo simulado también debe guardarse como FASTA.
- **REQ-TE-03.** Los TE deben poder insertarse al construir referencias y como
  variantes configuradas en individuos concretos. La primera versión no debe
  simular transposición espontánea, tasas de movilización ni estados activo o
  inactivo durante la evolución.
- **REQ-TE-04.** Una inserción explícita en un individuo debe indicar individuo,
  haplotipo o copia, familia o secuencia de TE y las demás propiedades del
  evento, sin asignarse por frecuencia poblacional.
- **REQ-TE-05.** Cada inserción individual debe elegir de forma excluyente una
  coordenada exacta o una selección aleatoria reproducible entre los sitios
  permitidos por la familia.
- **REQ-TE-06.** Una coordenada exacta debe validarse contra la referencia
  estructural, el haplotipo, los solapamientos y las reglas de destino. El modo
  aleatorio debe registrar los candidatos, la regla, la semilla y la posición.
- **REQ-TE-07.** Cada familia puede definir orientación directa o inversa,
  truncamiento y duplicación del sitio diana, además del caso simple de copia
  completa sin modificación ni duplicación del sitio diana.
- **REQ-TE-08.** Un evento debe conservar la secuencia fuente y final, el
  fragmento retenido, la orientación y la longitud y secuencia del sitio diana
  original o modificado.
- **REQ-TE-09.** Un evento solicitado puede ser de copiar y pegar, conservando la
  copia donante, o de cortar y pegar, eliminándola. En este último caso,
  eliminación e inserción deben formar un único evento compuesto y trazable.
- **REQ-TE-10.** Cortar y pegar debe actuar en el mismo haplotipo o copia salvo
  que el evento solicitado defina otra regla válida. Debe conservar los
  identificadores, coordenadas y secuencias de las copias donante y nueva.
- **REQ-TE-11.** Cada familia puede definir regiones permitidas, excluidas o
  favorecidas y motivos diana con pesos. Sin reglas, el destino debe elegirse
  uniformemente entre sitios válidos de la referencia o haplotipo efectivo.
- **REQ-TE-12.** Las reglas de destino deben actualizarse tras cambios de
  coordenadas. Si no existe un sitio válido, debe aplicarse la política
  configurada de advertencia o error sin crear una inserción inválida.
- **REQ-TE-13.** Cada copia debe tener identificador estable y conservar familia,
  tipo, secuencia fuente, copia de origen, haplotipo, posición, orientación,
  generación cuando corresponda y descendencia.
- **REQ-TE-14.** La verdad, el VCF y su complemento deben representar sin pérdida
  la inserción y, en eventos compuestos, la eliminación asociada. Deben usar la
  representación VCF estándar cuando sea suficiente.

#### 4.3.10. Recombinación de la primera versión

- **REQ-REC-01.** Toda la ejecución debe usar exclusivamente una tasa uniforme
  para todo el genoma o mapas completos de recombinación, conforme a
  REQ-CONF-46 a REQ-CONF-50.
- **REQ-REC-02.** La primera versión no debe modelar interferencia entre
  entrecruzamientos. Sus posiciones deben generarse de forma independiente según
  la tasa uniforme o el mapa completo y la semilla.
- **REQ-REC-03.** En copias homólogas con inversiones, duplicaciones o deleciones
  heterocigotas, solo deben permitirse entrecruzamientos en intervalos colineales
  y homólogos.
- **REQ-REC-04.** Las regiones estructuralmente incompatibles deben excluirse del
  espacio recombinante y registrarse. No deben generarse gametos desequilibrados
  por pérdida o ganancia de segmentos ni modelarse la fertilidad reducida.
- **REQ-REC-05.** Para cromosomas con ploidía par mayor que dos, cada meiosis debe
  formar bivalentes. Cada bivalente debe recombinar conforme a la tasa o mapa y
  aportar una copia reconstruida al gameto.
- **REQ-REC-06.** Las copias estructuralmente compatibles deben emparejarse al
  azar para formar bivalentes, de manera reproducible mediante la semilla. Los
  emparejamientos y haplotipos transmitidos deben quedar trazados internamente.
- **REQ-REC-07.** La primera versión no debe aplicar preferencia de
  apareamiento por subgenoma, emparejamiento homeólogo preferente, multivalentes
  ni doble reducción.
- **REQ-REC-08.** Cada parental debe aportar la mitad de la ploidía al
  descendiente. El método técnico para generar las posiciones de los
  entrecruzamientos se concretará por separado sin alterar las reglas
  funcionales acordadas.
- **REQ-REC-09.** Sin interferencia, el número de entrecruzamientos por bivalente
  y cromosoma debe obtenerse mediante una distribución de Poisson cuya media se
  derive de su longitud genética en centimorgan o Morgan. Las posiciones deben
  seguir la tasa uniforme o el mapa completo elegido. Si la muestra es cero en
  una secuencia clasificada como cromosoma, debe aplicarse un mínimo de un
  entrecruzamiento por bivalente y meiosis.
- **REQ-REC-10.** En una secuencia clasificada como *scaffold*, el mínimo de un
  entrecruzamiento solo debe aplicarse cuando su longitud supere el umbral
  efectivo. La regla debe usar exclusivamente longitud genética en cM o longitud
  física en bases, conforme a REQ-CONF-100. Cada *scaffold* debe conservar la
  regla efectiva y si procede del valor general o de una sustitución específica.
- **REQ-REC-11.** La salida detallada de recombinación debe generarse solo para
  los individuos incluidos en muestreos solicitados. Para cada uno debe indicar
  los puntos de recombinación y el origen parental de cada segmento. Las demás
  generaciones deben conservar resúmenes, no una exportación individual de todos
  sus segmentos.
- **REQ-REC-12.** Para cada segmento de un individuo muestreado, el origen debe
  poder reconstruirse hasta el individuo fundador y su copia haplotípica
  concreta. El sistema debe conservar internamente durante la evolución el
  historial suficiente para producir esa salida cuando se solicite el muestreo.
- **REQ-REC-13.** Una mutación de novo sobre un segmento heredado no debe
  sustituir su origen fundador. Debe registrarse como una capa separada con su
  identidad estable, el individuo y la generación en que apareció y los demás
  datos de verdad ya exigidos.
- **REQ-REC-14.** Cuando un cambio estructural altere las coordenadas, la
  trazabilidad de un individuo muestreado debe conservar las coordenadas del
  segmento en el haplotipo actual y su correspondencia con las coordenadas de la
  copia haplotípica fundadora cuando exista. Para duplicaciones, inversiones y
  desplazamientos debe conservar también la orientación y la relación entre
  segmentos.
- **REQ-REC-15.** Cada copia producida por una duplicación debe recibir un
  identificador estable propio, aunque varias copias procedan del mismo intervalo
  fundador. Cada copia debe mantener el vínculo a ese intervalo y poder seguirse
  por separado durante la herencia, la recombinación y las mutaciones posteriores.
- **REQ-REC-16.** Una inserción sin correspondencia con un haplotipo fundador,
  sea su secuencia generada o aportada, debe recibir un identificador estable y
  vincularse al evento, individuo y generación en que apareció. Su secuencia y
  descendencia deben rastrearse desde ese origen nuevo y distinguirse
  expresamente de un origen fundador.
- **REQ-REC-17.** La trazabilidad de una deleción heredada debe conservar el
  intervalo fundador ausente y el evento, individuo y generación en que ocurrió
  la pérdida. Debe permitir distinguir una región eliminada de otra que nunca
  formó parte de ese linaje.
- **REQ-REC-18.** La primera versión no debe admitir variantes estructurales
  anidadas ni un evento estructural que modifique una estructura creada por un
  evento previo dentro de la misma copia haplotípica y su linaje. Antes de
  producir resultados, el sistema debe rechazar la configuración, identificar
  los eventos en conflicto y explicar el error; no debe simplificarla. Dos
  haplotipos distintos sí pueden contener variantes estructurales solapadas sobre
  la misma región de referencia.
- **REQ-REC-19.** Una misma copia haplotípica debe poder contener varias
  variantes estructurales no solapadas, incluidas variantes inmediatamente
  adyacentes. Cada evento debe mantener un identificador y una trazabilidad
  propios. La validación debe distinguir la adyacencia permitida del solapamiento
  o anidamiento prohibido.
- **REQ-REC-20.** Los SNP e *indels* pequeños deben poder aparecer dentro de
  segmentos duplicados, invertidos o insertados. Cada copia puede acumular
  cambios propios, asociados a su identificador estable y a sus coordenadas
  actuales, sin perder la procedencia estructural. Estos cambios no deben
  considerarse anidamiento de variantes estructurales. La primera versión debe
  rechazar todo SNP o *indel* que solape un punto de ruptura de una variante
  estructural; un *indel* que atraviese el punto también solapa. Las variantes
  pequeñas inmediatamente adyacentes deben permitirse.
- **REQ-REC-21.** La primera versión debe usar exclusivamente
  entrecruzamientos como mecanismo de recombinación. No debe modelar conversión
  génica ni interpretarla como una translocación no recíproca.
- **REQ-REC-22.** La segregación debe ser mendeliana. Dentro de cada bivalente,
  cada copia estructuralmente compatible debe tener la misma probabilidad de
  transmisión, de forma reproducible mediante la semilla.
- **REQ-REC-23.** La primera versión no debe modelar distorsión de segregación ni
  impulso meiótico. Estas exclusiones no deben alterar ni confundirse con la
  selección de viabilidad aplicada antes de la reproducción.
- **REQ-REC-24.** Los SNP e *indels* pequeños deben poder solaparse dentro de una
  misma copia haplotípica. Cuando sean explícitos, la configuración debe indicar
  su orden y el sistema debe aplicarlos secuencialmente. Las coordenadas de cada
  evento deben interpretarse sobre el estado del haplotipo inmediatamente
  anterior a aplicarlo.
- **REQ-REC-25.** La verdad genotípica actual y las salidas deben derivarse de la
  secuencia haplotípica final. El historial debe marcar cada evento anterior
  como activo, modificado o eliminado según el efecto que conserve después de
  aplicar los eventos posteriores. Debe mantener la correspondencia con la
  referencia original cuando exista y usar identificadores de segmento para las
  secuencias introducidas por inserciones o duplicaciones.
- **REQ-REC-26.** Una variante anterior eliminada por completo no debe aparecer
  en la verdad actual ni en el VCF. En configuración directa, un SNP oculto por
  una deleción no debe emitirse como variante actual. La prohibición de
  REQ-REC-20 para SNP o *indels* que solapen o atraviesen un punto de ruptura
  estructural sigue vigente.
- **REQ-REC-27.** Antes de producir salidas, el sistema debe validar los eventos
  explícitos en su orden. Si uno apunta a una base o segmento eliminado por un
  evento anterior, la configuración es inválida y debe detenerse con un error
  que identifique el evento y el objetivo ausente. No debe ignorarlo ni
  recolocarlo.
- **REQ-REC-28.** Las variantes pequeñas situadas en haplotipos distintos deben
  resolverse de forma independiente. Sus alelos finales pueden coexistir como
  ALT distintos del mismo registro VCF multialélico, y `GT` debe conservar su
  fase y dosis.
- **REQ-REC-29.** Cuando varios eventos pequeños solapados de un mismo haplotipo
  produzcan un único alelo final, el VCF debe representarlo como una sola
  variante o ALT compleja normalizada. El historial estructurado debe conservar
  la cadena ordenada de eventos y el estado activo, modificado o eliminado de
  cada uno.

La validación funcional del grafo queda confirmada. Solo queda aplazada la
sintaxis con la que se representarán técnicamente el grafo de cruces, los
parentales expresos y las combinaciones de condiciones. Las plantillas iniciales
ya están delimitadas y no restringen los cruces personalizados.

### 4.4. Datos por muestra

La referencia de mapeo es la referencia respecto a la que se registran los
alineamientos de una muestra en su BAM o CRAM. Puede ser distinta de la
referencia cuya estructura corresponde a la muestra.

La cabecera es la parte del archivo que contiene metadatos, es decir,
información que describe el propio archivo y su contexto. En los requisitos
siguientes, la cabecera y los metadatos deben ajustarse al estándar del formato,
aunque los campos y las versiones concretas todavía no se han decidido.

Una lectura corta representa un fragmento relativamente breve de DNA. En modo
de extremo único (*single-end*) se observa un extremo; en pares de extremos
(*paired-end*) se observan ambos extremos del mismo fragmento. Una lectura larga
representa fragmentos mucho mayores. Illumina es el ejemplo principal de
lecturas cortas admitidas en la primera versión, y PacBio HiFi es la tecnología
de lecturas largas admitida.

El modo realista simula lecturas desde los haplotipos, aplica la contaminación
y los errores de demultiplexado, guarda el FASTQ final y después mapea
exactamente esas lecturas. Un fichero de procedencia es un registro estructurado
que relaciona una lectura con su fuente real y con el mecanismo que la introdujo
en la muestra. El modo directo genera registros de alineamiento sintéticos sin
crear ni mapear previamente FASTQ; sirve para casos pequeños y controlados y no
debe confundirse con el resultado de un mapeador.
En este documento, la mención confirmada «minimap» se interpreta como minimap2.
Un *preset* es un grupo de opciones predeterminadas del mapeador ajustadas a un
tipo de lecturas.

La puntuación de calidad de una base expresa, en una escala probabilística, la
confianza en que la base observada sea correcta. Una puntuación menor corresponde
a una probabilidad mayor de error.

- **REQ-DAT-01.** Cuando la salida BAM o CRAM esté activa, el sistema debe
  generarla para las muestras incluidas. Debe usar CRAM de forma predeterminada,
  BAM cuando se seleccione esta alternativa y omitirla solo cuando se haya
  desactivado expresamente.
- **REQ-DAT-02.** Los BAM o CRAM de una muestra deben poder generarse mapeados
  contra cada referencia de mapeo seleccionada entre las disponibles.
- **REQ-DAT-03.** La referencia de mapeo debe poder ser una referencia
  estructural distinta de la que corresponde a la muestra.
- **REQ-DAT-04.** Cada BAM o CRAM debe incluir en su propia cabecera o metadatos
  estándar la información necesaria para identificar sin ambigüedad la
  referencia de mapeo, la muestra, el evento, la biblioteca, la tanda o carril y
  cada grupo de lectura, directamente o mediante sus relaciones documentadas.
- **REQ-DAT-05.** El sistema debe permitir introducir errores, contaminación,
  cobertura baja o cero, mala calidad de mapeo y otras degradaciones en los
  datos de lectura o alineamiento según el modo utilizado.
- **REQ-DAT-06.** Cuando estén activos BAM o CRAM, el sistema debe generar un
  conjunto independiente para cada referencia de mapeo seleccionada.
- **REQ-DAT-07.** Si se seleccionan varias referencias de mapeo, todos los
  conjuntos deben contener las mismas muestras.
- **REQ-DAT-08.** Cuando un muestreo evolutivo solicite BAM o CRAM, sus
  individuos deben convertirse en muestras y seguir las mismas referencias de
  mapeo seleccionadas y las mismas reglas de conjuntos independientes que las
  demás muestras. Todos los conjuntos de ese muestreo deben contener los mismos
  individuos muestreados.
- **REQ-DAT-09.** Cada BAM o CRAM y sus metadatos deben usar el identificador
  estable del individuo del que procede la muestra.
- **REQ-DAT-10.** Cuando una población intermedia de cruces solicite BAM o CRAM,
  todos sus individuos deben convertirse en muestras y seguir las referencias
  de mapeo seleccionadas y las reglas de conjuntos independientes. Todos los
  conjuntos de esa población deben contener los mismos individuos.
- **REQ-DAT-11.** Los BAM o CRAM deben reflejar la secuencia realizada de las
  repeticiones en tándem de cada muestra, incluidas sus diferencias respecto a la
  referencia de mapeo.
- **REQ-DAT-12.** El sistema debe simular lecturas cortas tipo Illumina y
  tecnologías similares, tanto *single-end* como, especialmente, *paired-end*.
  Debe permitir configurar longitud, cobertura, errores y calidad y, para
  *paired-end*, orientación y distribución de tamaños de inserto.
- **REQ-DAT-13.** El sistema debe simular lecturas largas PacBio HiFi, con
  distribución configurable de longitudes, cobertura, errores y calidad.
- **REQ-DAT-14.** Las lecturas tipo Illumina *single-end* y *paired-end* y PacBio
  HiFi deben ofrecer perfiles predeterminados validados y permitir ajustar sus
  parámetros. Las distribuciones, tasas y valores cuantitativos exactos de
  longitud, error y calidad se decidirán durante la planificación de la
  implementación y deberán comprobarse con datos adecuados para cada perfil.
- **REQ-DAT-15.** En el modo realista, las lecturas deben generarse desde los
  haplotipos faseados de cada individuo y usar las decisiones aleatorias
  reproducibles de la ejecución.
- **REQ-DAT-16.** Cada FASTQ, BAM o CRAM debe poder relacionarse con la
  tecnología, el perfil y los parámetros efectivos con los que se generaron sus
  lecturas o alineamientos.
- **REQ-DAT-17.** En el modo realista, el sistema debe guardar como salida los
  FASTQ finales después de aplicar la contaminación biológica y los errores de
  demultiplexado, y usar exactamente sus lecturas para generar los BAM o CRAM.
  No debe reconstruir esos FASTQ a partir de los archivos de alineamiento.
- **REQ-DAT-18.** Los FASTQ del modo realista deben conservar la secuencia, la
  calidad, el identificador de cada lectura y su emparejamiento cuando exista.
  El sistema debe conservar los mismos identificadores y pares al mapearlos y no
  debe descartar las lecturas no mapeadas que produzca el mapeador. Cada lectura
  *single-end* o larga, y cada pareja *paired-end* considerada como una unidad,
  debe tener una identidad técnica única dentro del evento.
- **REQ-DAT-19.** En el modo directo, el sistema debe generar BAM o CRAM
  sintéticos sin simular ni mapear FASTQ previamente. Cada archivo y su
  descripción en el manifiesto deben identificar claramente este origen
  sintético y sus parámetros efectivos.
- **REQ-DAT-20.** Los FASTQ, BAM y CRAM de ambos modos deben cumplir los
  estándares de sus formatos y la versión oficial que se elija durante la
  planificación. El modo directo no debe requerir FASTQ.
- **REQ-DAT-21.** Cada experimento de mapeo debe usar exactamente un mapeador. Si
  la configuración no elige otro compatible, debe usar minimap2.
- **REQ-DAT-22.** Si un experimento selecciona varias referencias de mapeo, debe
  aplicar el mismo mapeador y los mismos FASTQ a cada referencia y producir
  conjuntos independientes. El mapeador puede ser distinto entre eventos, pero
  cada evento debe generar sus FASTQ propios y no reutilizarlos para comparar
  mapeadores.
- **REQ-DAT-23.** Cuando se use minimap2, el perfil de lecturas debe seleccionar
  automáticamente el *preset* y los parámetros predeterminados apropiados para
  lecturas cortas tipo Illumina o lecturas largas PacBio HiFi. Las sustituciones
  de parámetros deben prevalecer, generar una advertencia y validarse para
  detectar opciones desconocidas o incompatibles cuando sea posible.
- **REQ-DAT-24.** La cobertura realizada debe ser estocástica y aproximarse a la
  cobertura media objetivo efectiva de cada individuo, cromosoma o *scaffold* y
  región. En el modo realista debe aplicarse al generar los FASTQ; en el modo
  directo, al generar los BAM o CRAM sintéticos. Las regiones se interpretan en
  la referencia estructural del individuo. No debe forzarse una profundidad
  exacta en cada base.
- **REQ-DAT-25.** En el modo realista, las lecturas contaminantes de todas las
  fuentes biológicas configuradas deben sustituir lecturas propias de la muestra
  antes del mapeo. El único mapeador debe procesar el resultado completo contra
  cada referencia seleccionada. En modo directo deben aplicarse reglas
  equivalentes sobre los registros sintéticos y conservarse su procedencia.
- **REQ-DAT-26.** El sistema debe conservar en la verdad y el manifiesto el
  identificador de cada instancia contaminante y el identificador estable de la
  lectura o pareja de origen. Esta información debe
  mantenerse fuera del FASTQ y del BAM o CRAM cuando exponerla aportaría al
  programa evaluado más metadatos de contaminación que los de un archivo real.
- **REQ-DAT-27.** Cuando se contabilice la contaminación realizada en datos
  *paired-end*, cada unidad debe incorporar juntos los dos
  *mates* de la misma fuente y conservar nombres, indicadores del formato,
  orientación y relación de pareja coherentes. Esta operación no debe crear
  pares huérfanos.
- **REQ-DAT-28.** Los FASTQ externos usados como contaminación deben validarse
  antes de generar datos y copiarse al directorio de resultados, exclusivamente
  en modo realista. La simulación debe usar la copia y no modificar el original.
  Si la validación o la compatibilidad con la biblioteca receptora fallan, debe
  detenerse antes de generar datos y explicar los errores. La selección o
  submuestreo de sus lecturas debe ser reproducible con la semilla de la
  ejecución y respetar la cobertura objetivo o el número de secuencias o pares
  configurado. En datos *paired-end* debe preservar el emparejamiento y no crear
  pares huérfanos.
- **REQ-DAT-29.** El sistema debe aplicar y registrar por separado la
  contaminación biológica y la técnica. La contaminación técnica debe ocurrir
  entre muestras pertenecientes a la misma tanda o conjunto de bibliotecas
  multiplexadas configurado.
- **REQ-DAT-30.** La contaminación técnica debe poder producir dos resultados:
  lecturas transferidas de una muestra a otra, que desaparecen de la donante sin
  duplicarse, y lecturas que quedan sin asignar o son filtradas. La procedencia y
  las cantidades realizadas deben ser reproducibles.
- **REQ-DAT-31.** El uso de UDI debe permitir reconocer combinaciones de índices
  inesperadas y simular su filtrado durante el demultiplexado. No debe asumirse
  que los UDI impiden físicamente el *index hopping*.
- **REQ-DAT-32.** El sistema debe simular directamente el efecto final del
  demultiplexado sobre cada muestra: debe mover o reasignar las lecturas
  afectadas a otra muestra, o marcarlas como descartadas o sin asignar según el
  esquema de índices. Debe conservar internamente la procedencia verdadera de
  cada lectura como parte de la verdad conocida.
- **REQ-DAT-33.** En el modo realista, el sistema debe conservar un único FASTQ
  final, o una pareja de FASTQ para datos *paired-end*, por combinación de
  individuo, biblioteca y tanda o carril. No debe guardar además una copia limpia
  anterior a la contaminación o a los errores de demultiplexado.
- **REQ-DAT-34.** Por cada FASTQ final, el sistema debe generar un fichero
  estructurado de procedencia que identifique cada lectura contaminante o, en
  datos *paired-end*, la pareja completa como una unidad, e indique su fuente
  real, la instancia de origen y el mecanismo que la introdujo. Los mecanismos deben distinguir, como
  mínimo, contaminación biológica, cruce entre bibliotecas, reasignación por
  índices y otros mecanismos configurados.
- **REQ-DAT-35.** El fichero de procedencia debe referirse a las lecturas o
  parejas mediante sus identificadores técnicos únicos y enlazar por separado
  su identificador estable de origen. No debe duplicar sus secuencias ni sus
  calidades. Junto con el FASTQ final y el manifiesto, debe permitir reconstruir
  la verdad conocida de la contaminación.
- **REQ-DAT-36.** Las lecturas contaminantes que permanezcan en el FASTQ final
  deben registrarse individualmente por identificador técnico de lectura o de
  pareja y con su origen. Para las lecturas descartadas o sin asignar no deben conservarse
  secuencias ni enumerarse identificadores; el fichero de procedencia debe
  registrar únicamente recuentos agregados, desglosados como mínimo por muestra
  y por causa u origen cuando se conozcan.
- **REQ-DAT-37.** El sistema debe aplicar el perfil de errores efectivo de la
  tecnología y su versión al generar las lecturas. Debe conservar el perfil
  predeterminado y todas las sustituciones globales o específicas realmente
  aplicadas.
- **REQ-DAT-38.** Las puntuaciones de calidad por base del FASTQ deben generarse
  de forma coherente con el modelo de errores efectivo: una base de menor
  calidad debe tener mayor probabilidad de ser errónea conforme a la
  interpretación probabilística de la puntuación.
- **REQ-DAT-39.** Cuando el modelo de la tecnología contemple dependencias de la
  posición dentro de la lectura o del contexto de bases, los errores y las
  puntuaciones de calidad deben reflejar conjuntamente esos patrones.
- **REQ-DAT-40.** El mapeo debe propagar correctamente a BAM o CRAM las
  puntuaciones de calidad del FASTQ final. La verdad conocida debe conservar,
  para cada lectura simulada, la diferencia entre la secuencia original y la
  observada para identificar los errores introducidos.
- **REQ-DAT-41.** A partir del FASTQ real y su referencia correspondiente, el
  sistema debe estimar, según la información disponible, un perfil WGS con
  distribuciones de longitud de lectura, calidad por base y posición, tipos y
  tasas de error y otros parámetros tecnológicos útiles para la simulación. En
  la primera versión, esta calibración debe limitarse a lecturas cortas tipo
  Illumina y lecturas largas PacBio HiFi.
- **REQ-DAT-42.** La estimación debe intentar separar los errores de
  secuenciación de las variantes reales, las regiones ambiguas y los errores de
  mapeo. Debe documentar sus limitaciones y advertir cuando la referencia o los
  datos no permitan una estimación fiable.
- **REQ-DAT-43.** Si la configuración selecciona el perfil empírico estimado, sus
  parámetros y versión deben convertirse en el perfil efectivo reproducible al
  que se aplicarán las sustituciones de ámbito configuradas. Si no lo selecciona,
  el sistema debe mantener el perfil tecnológico predeterminado.
- **REQ-DAT-44.** Cada perfil empírico estimado debe guardarse en los resultados
  como un fichero estructurado reutilizable. Debe contener la versión de su
  esquema, los parámetros estimados, la tecnología declarada o inferida, la
  referencia y los datos de entrada usados, el método y la versión de las
  herramientas, las advertencias y las medidas de calidad de la estimación.
- **REQ-DAT-45.** Una simulación posterior debe poder usar el fichero de perfil
  empírico sin volver a analizar ni necesitar el FASTQ original. El formato
  concreto del artefacto se decidirá durante la planificación de la
  implementación.
- **REQ-DAT-46.** El sistema debe generar por separado duplicados de PCR y
  duplicados ópticos con las tasas efectivas del perfil tecnológico, únicamente
  cuando el mecanismo sea pertinente para ese perfil. Los duplicados de PCR
  deben estar ausentes en la preparación WGS sin PCR predeterminada y poder
  aparecer en un perfil WGS con amplificación.
- **REQ-DAT-47.** En datos *paired-end*, una duplicación debe afectar a la pareja
  completa y conservar una relación de emparejamiento coherente. No debe crear
  pares huérfanos.
- **REQ-DAT-48.** La verdad conocida debe identificar para cada duplicado su tipo
  y la lectura o pareja de origen, y conservar las cantidades esperadas y
  realizadas de cada mecanismo.
- **REQ-DAT-49.** El estado marcado o sin marcar elegido debe aplicarse de forma
  coherente a los BAM o CRAM. El estado predeterminado debe ser marcado y una
  salida sin marcar debe requerir una solicitud expresa. En ambos casos, la
  verdad conocida debe conservar cuáles son los duplicados, su tipo y su lectura
  o pareja de origen.
- **REQ-DAT-50.** Cuando se active la degradación, el sistema debe añadir restos
  de adaptadores sin recortar a las lecturas conforme al perfil o secuencia, la
  frecuencia o proporción, la longitud residual y el extremo efectivos.
- **REQ-DAT-51.** La verdad conocida debe identificar cada lectura afectada y el
  intervalo de secuencia de adaptador añadido.
- **REQ-DAT-52.** Cuando el sesgo de contenido GC esté activo, la probabilidad de
  obtener lecturas de una región debe modificarse conforme a la curva o relación
  y a la intensidad efectivas del perfil tecnológico.
- **REQ-DAT-53.** El sesgo de contenido GC debe combinarse de forma reproducible
  con la cobertura objetivo y las sustituciones regionales sin ocultar ninguna
  excepción explícita. La cobertura realizada debe poder resumirse por
  intervalos adecuados para comprobar el efecto.
- **REQ-DAT-54.** En el modo realista, las lecturas de regiones de baja
  mapeabilidad deben generarse desde sus posiciones y haplotipos verdaderos. El
  mapeador debe determinar los alineamientos primarios o secundarios, las
  posiciones observadas y la calidad de mapeo MAPQ según la referencia usada.
- **REQ-DAT-55.** El sistema no debe sobrescribir artificialmente la MAPQ
  producida por el mapeador para simular baja mapeabilidad en el modo realista.
- **REQ-DAT-56.** La verdad conocida debe conservar para cada lectura su posición
  y haplotipo de origen, de modo que puedan medirse la ambigüedad y los mapeos en
  posiciones incorrectas.
- **REQ-DAT-57.** En el modo directo sintético, cada lectura debe asignarse de
  forma reproducible a una de las categorías configuradas: correctamente
  alineada, alineada en una posición incorrecta, multimapeada o no alineada.
- **REQ-DAT-58.** Los campos, las banderas, el CIGAR, los alineamientos
  secundarios o suplementarios y la MAPQ de cada registro sintético deben ser
  internamente coherentes con su categoría, la referencia y la secuencia de la
  lectura.
- **REQ-DAT-59.** La verdad conocida del modo directo debe conservar la posición
  y el haplotipo de origen y la categoría simulada de cada lectura, con
  independencia de la posición observada en el BAM o CRAM.
- **REQ-DAT-60.** Para cada lectura quimérica, la verdad conocida debe registrar
  todos sus fragmentos, las muestras y haplotipos de procedencia, sus coordenadas
  y orientaciones de origen y los puntos de unión entre fragmentos.
- **REQ-DAT-61.** En el modo realista, las lecturas quiméricas deben incluirse en
  el FASTQ final y el mapeador debe determinar sus alineamientos divididos o
  suplementarios según la referencia usada, sin que el sistema sobrescriba el
  resultado.
- **REQ-DAT-62.** En el modo directo, los CIGAR, banderas, posiciones, MAPQ y
  alineamientos suplementarios de una lectura quimérica deben generarse de forma
  internamente coherente con sus fragmentos y con la categoría de alineamiento
  simulada.
- **REQ-DAT-63.** Todos los fragmentos de una lectura quimérica deben proceder
  del mismo individuo y de la misma biblioteca. La verdad conocida debe conservar
  el origen de cada fragmento. Los intercambios entre muestras no deben
  representarse como quimeras, sino mediante los modelos de contaminación o
  errores de demultiplexado.
- **REQ-DAT-64.** Dentro del mismo individuo y biblioteca, una lectura quimérica
  debe poder unir fragmentos de cromosomas o *scaffolds* diferentes y de
  haplotipos o copias cromosómicas diferentes, según la categoría efectiva del
  perfil o la configuración. Para cada fragmento, la verdad debe registrar el
  cromosoma o *scaffold*, el haplotipo o copia, las coordenadas, la orientación y
  el punto de unión.
- **REQ-DAT-65.** Para lecturas *paired-end*, la configuración debe permitir
  indicar proporciones de parejas con orientación anómala, tamaño de inserto
  inesperado o un miembro no alineado, incluidas las combinaciones admitidas de
  estos estados. Las parejas normales deben seguir la orientación y la
  distribución de inserto del perfil tecnológico efectivo.
- **REQ-DAT-66.** En el modo realista, las parejas deben generarse desde su
  fragmento y estado verdaderos y el mapeador debe determinar el estado
  observado. La verdad conocida debe conservar el fragmento y el estado de
  origen de cada pareja.
- **REQ-DAT-67.** En el modo directo, las banderas, RNEXT y PNEXT, TLEN, CIGAR,
  orientación y estado *proper-pair* deben ser internamente coherentes con la
  anomalía simulada y los alineamientos de ambos miembros.
- **REQ-DAT-68.** Una instancia de lectura debe ser una observación concreta
  presente en una biblioteca. Su identificador técnico debe ser único dentro de
  su evento. Los dos miembros de una pareja deben compartir el nombre de pareja
  exigido por el formato y diferenciarse mediante sus indicadores de miembro;
  ese nombre de pareja no debe identificar otra pareja del mismo evento.
- **REQ-DAT-69.** La contaminación biológica debe crear en la receptora una
  instancia nueva con identificador técnico nuevo y un vínculo al identificador
  estable de la lectura o pareja donante. La instancia original debe permanecer
  en la donante. Una misma unidad de origen no debe copiarse dos veces dentro de
  una misma biblioteca receptora; puede originar una instancia en otra
  receptora compatible. Una duplicación adicional dentro de la receptora solo
  debe aparecer mediante un mecanismo de duplicación configurado y trazable.
- **REQ-DAT-70.** Una reasignación técnica por demultiplexado debe transferir la
  misma instancia y conservar su identificador: debe retirarla de la donante e
  incorporarla a la receptora sin crear otra. El fichero de procedencia debe
  distinguir sin ambigüedad una copia biológica de una transferencia técnica.
- **REQ-DAT-71.** El sistema debe validar antes de combinar FASTQ o producir el
  BAM o CRAM final que no existan colisiones de identificadores técnicos dentro
  del evento. Una colisión debe ser un error; no debe resolverse renombrando una
  lectura sin actualizar de forma atómica FASTQ, BAM o CRAM, verdad,
  procedencia y manifiesto.

- **REQ-WGS-01.** Cada biblioteca debe declarar exactamente uno de dos objetivos
  excluyentes: cobertura media o número de lecturas. Declarar ambos o ninguno
  debe ser un error previo a la generación.
- **REQ-WGS-02.** En datos *paired-end*, un objetivo por número debe expresarse
  en parejas completas. En datos *single-end* o de lecturas largas debe
  expresarse en lecturas individuales.
- **REQ-WGS-03.** El sistema debe derivar y registrar el valor realizado
  equivalente del indicador no elegido cuando sea posible: número de lecturas o
  parejas para un objetivo de cobertura, y cobertura para un objetivo numérico.
- **REQ-WGS-04.** La cobertura en «x» debe definirse como el número de bases de
  lecturas generadas dividido por la longitud haploide de la referencia
  seleccionada, con independencia de la ploidía del individuo.
- **REQ-WGS-05.** Antes de sesgos y degradaciones, el modelo WGS predeterminado
  debe elegir de forma uniforme y reproducible las posiciones iniciales de los
  fragmentos o lecturas entre las copias cromosómicas accesibles.
- **REQ-WGS-06.** Antes de aplicar sesgos, la cobertura local realizada debe ser
  estocástica y aproximadamente compatible con una distribución de Poisson
  alrededor de la media objetivo. No debe imponerse una profundidad exacta por
  base.
- **REQ-WGS-07.** Cuando las secuencias tengan ploidías distintas, el muestreo
  debe repartir las lecturas en proporción a la longitud de cada copia y al
  número de copias. En ausencia de excepciones, cada copia debe recibir
  aproximadamente la misma cobertura y la profundidad total de un cromosoma debe
  ser proporcional a su ploidía.
- **REQ-WGS-08.** El sesgo GC, las regiones con cobertura especial, los límites
  de secuencia y las demás degradaciones deben aplicarse después del muestreo
  WGS básico, de forma documentada y con las precedencias ya acordadas.
- **REQ-WGS-09.** Antes de esos efectos, los haplotipos y alelos deben contribuir
  en proporción a su dosis y ploidía. Por ejemplo, 40× totales en un tetraploide
  uniforme deben producir aproximadamente 10× por copia.
- **REQ-WGS-10.** La primera versión no debe ofrecer parámetros para forzar un
  sesgo de biblioteca por alelo ni la pérdida completa de uno de los alelos.
  Las diferencias realizadas solo deben proceder del muestreo aleatorio, la
  cobertura y los errores, sesgos y resultados de mapeo admitidos.
- **REQ-WGS-11.** El manifiesto y el informe de validación deben registrar la
  cobertura total respecto a la longitud haploide y, cuando pueda medirse, la
  cobertura realizada por cromosoma y copia o haplotipo.
- **REQ-WGS-12.** Solo debe conservarse el BAM o CRAM final que contiene las
  degradaciones configuradas. No debe guardarse simultáneamente una versión
  limpia del mismo experimento.
- **REQ-WGS-13.** Un control sin degradaciones debe ser otro experimento que
  reutilice la misma verdad genotípica, referencias y muestras, y la semilla o
  flujo pertinente, con las degradaciones desactivadas. El manifiesto debe
  enlazar los experimentos relacionados sin mezclar sus resultados.

- **REQ-CONT-01.** La contaminación biológica debe ser sustitutiva: por cada
  lectura o pareja contaminante incorporada a la receptora debe retirarse una
  lectura o pareja propia. La biblioteca donante no debe perder datos por este
  mecanismo.
- **REQ-CONT-02.** La fracción total de contaminación biológica de una biblioteca
  o muestra debe admitir valores entre 0 % y 100 %, inclusive. El 100 % debe
  representar una sustitución completa o intercambio de muestra.
- **REQ-CONT-03.** Con varias fuentes, sus contribuciones deben sumarse y no
  superar el 100 % ni el número de lecturas o parejas disponibles. Un exceso
  debe detener la ejecución como error de configuración. Esta comprobación debe
  realizarse después de convertir conjuntamente todas las unidades conforme a
  REQ-CONT-15 y antes de seleccionar lecturas.
- **REQ-CONT-04.** Toda fuente debe ser compatible con la tecnología, el perfil
  y la disposición de lecturas de la biblioteca receptora. No deben mezclarse
  lecturas cortas y largas, perfiles tecnológicos incompatibles ni disposiciones
  *single-end* y *paired-end* incompatibles.
- **REQ-CONT-05.** Las fuentes predeterminadas deben ser otras muestras o
  bibliotecas compatibles de la misma tanda. Los FASTQ externos solo deben
  admitirse en modo realista y deben declarar o superar la validación de
  compatibilidad.
- **REQ-CONT-06.** La configuración debe poder expresar una matriz dirigida de
  contaminación entre muestras o bibliotecas compatibles de una tanda. Cada
  celda donante→receptora debe indicar su contribución; la matriz no necesita ser
  simétrica y la diagonal representa lecturas propias, no contaminación.
- **REQ-CONT-07.** Las entradas de cada receptora en la matriz, sumadas a sus
  fuentes externas, deben respetar el límite total del 100 %. Antes de generar
  deben validarse identificadores, compatibilidad, unidades y sumas. La diagonal
  debe representar el remanente propio de la receptora después de las
  contribuciones contaminantes, no una fuente adicional que vuelva a sumarse.
- **REQ-CONT-08.** Si se configura una tasa total sin matriz, los donantes deben
  elegirse reproduciblemente entre las demás muestras compatibles de la tanda y
  ponderarse por su número disponible de lecturas o parejas. Una matriz explícita
  debe tener prioridad y desactivar este reparto automático.
- **REQ-CONT-09.** El manifiesto y el informe deben registrar la matriz objetivo
  y realizada o, en el reparto automático, los pesos, donantes seleccionados y
  aportes realizados.
- **REQ-CONT-10.** En errores de índices o demultiplexado, cada lectura o pareja
  reasignada debe eliminarse de la donante e incorporarse a la receptora sin
  duplicarse. El total de la tanda debe conservarse salvo los descartes o datos
  sin asignar.
- **REQ-CONT-11.** Las lecturas o parejas descartadas o sin asignar durante el
  demultiplexado no deben reemplazarse. Deben reducir el número de lecturas y la
  cobertura final de la muestra; en *paired-end* debe descartarse la pareja
  completa conforme a las reglas acordadas.
- **REQ-CONT-12.** El modo directo debe aplicar las mismas reglas funcionales de
  sustitución biológica, matriz o reparto automático, transferencia técnica y
  descarte mediante registros de alineamiento equivalentes, con campos,
  banderas y procedencia coherentes.
- **REQ-CONT-13.** Los objetivos y resultados de contaminación deben ser
  comparables entre modos con la misma configuración y semilla, aunque no se
  exija identidad de las decisiones o archivos. Configurar un FASTQ externo en
  modo directo debe ser un error previo a la generación.
- **REQ-CONT-14.** Por muestra, fuente y mecanismo, el manifiesto y el informe
  deben registrar objetivos, lecturas propias retiradas, contaminantes
  incorporadas, lecturas transferidas, descartadas o sin asignar, totales finales
  y cobertura final realizada.
- **REQ-CONT-15.** Las contribuciones expresadas como proporción, cobertura
  equivalente o número deben convertirse conjuntamente después de conocer el
  número limpio realizado de unidades de la biblioteca receptora. Una
  proporción debe aplicarse a ese número limpio. Un número debe conservar su
  valor entero. Una cobertura equivalente debe convertirse usando las bases que
  aportan las unidades contaminantes y el mismo denominador biológico que la
  cobertura objetivo del evento. En *paired-end*, las bases y el recuento deben
  corresponder siempre a parejas completas.
- **REQ-CONT-16.** La asignación entera conjunta debe elegir, entre las
  soluciones que respeten los límites y la disponibilidad, la que minimice la
  desviación global respecto a todas las contribuciones convertidas. El total
  entero de unidades contaminantes debe ser el entero más próximo a la suma de
  los objetivos convertidos. Los empates deben resolverse de forma reproducible
  con la subsemilla del evento. Deben registrarse por fuente el valor y la unidad
  originales, el objetivo convertido, el recuento entero realizado, la
  desviación y cualquier empate. Si no existe una asignación válida, el evento
  debe rechazarse antes de modificar los FASTQ finales.

Un grupo de lectura reúne lecturas obtenidas para una combinación concreta de
muestra, biblioteca y tanda o carril de secuenciación. En BAM o CRAM se describe
mediante una entrada `@RG` de la cabecera y la etiqueta `RG` de cada registro.

- **REQ-RG-01.** Los BAM o CRAM deben incluir grupos de lectura coherentes que
  identifiquen como mínimo la muestra o individuo, la biblioteca, la tanda o
  carril cuando corresponda y la plataforma o tecnología.
- **REQ-RG-02.** Un mismo individuo debe poder tener varias bibliotecas, tandas o
  carriles, cada uno con un identificador único y parámetros propios
  configurables.
- **REQ-RG-03.** Antes de generar lecturas, el sistema debe rechazar las
  colisiones de identificadores de grupos de lectura, bibliotecas, tandas o
  carriles e identificar los valores en conflicto. El identificador de grupo de
  lectura debe ser único en toda la ejecución.
- **REQ-RG-04.** Los identificadores de muestra, biblioteca, tanda, carril y
  grupo de lectura deben ser consistentes entre FASTQ, BAM o CRAM, manifiesto y
  ficheros de procedencia. El identificador de grupo debe distinguir evento,
  biblioteca, tanda y carril. Los identificadores técnicos de lectura deben
  cumplir además la unicidad y los vínculos de origen de REQ-DAT-68 a
  REQ-DAT-71.
- **REQ-RG-05.** Dentro de cada conjunto definido por una referencia de mapeo,
  el sistema debe generar un único BAM o CRAM final por individuo que combine
  todas sus bibliotecas, tandas y carriles y preserve sus grupos de lectura. No
  debe mezclar referencias ni requerir como producto principal archivos finales
  separados por biblioteca; los intermedios técnicos quedan fuera de las
  salidas principales.
- **REQ-RG-06.** Los FASTQ finales deben entregarse separados por individuo,
  biblioteca y tanda o carril, sin concatenarlos como producto principal. Cada
  archivo o pareja de archivos debe conservar su propio perfil y parámetros y
  corresponder a un único grupo de lectura del BAM o CRAM final combinado.
- **REQ-RG-07.** El campo de muestra `SM` de BAM o CRAM y la columna de muestra
  del VCF deben usar el identificador biológico estable del individuo. El evento
  debe identificarse mediante el directorio, el manifiesto y el grupo de lectura,
  sin modificar `SM`.

La expresión «otras degradaciones» no amplía por sí sola el alcance de v1. Solo
son obligatorias las degradaciones descritas mediante requisitos concretos en
este documento. Añadir otra requerirá una decisión funcional explícita y una
ampliación de la especificación.

### 4.5. VCF de verdad y coherencia de las salidas

- **REQ-VER-01.** Cuando la configuración solicite VCF de verdad, el sistema debe
  generarlo con los genotipos reales simulados de los SNP, *indels*, variantes
  estructurales y repeticiones en tándem de los individuos incluidos.
- **REQ-VER-02.** Cada genotipo del VCF de verdad debe corresponder al individuo
  simulado identificado en el archivo.
- **REQ-VER-03.** Cada VCF de verdad debe estar referido exactamente a la
  referencia de mapeo seleccionada para su conjunto y, cuando existan BAM o CRAM
  correspondientes, al mismo genoma que esos archivos.
- **REQ-VER-04.** Cuando se solicite VCF, el sistema debe generar uno de verdad
  independiente para cada referencia de mapeo seleccionada.
- **REQ-VER-05.** Cada conjunto debe usar exactamente una referencia de mapeo.
  Todos los BAM, CRAM o VCF que existan en él deben estar referidos a esa misma
  referencia. No se deben mezclar datos referidos a genomas distintos.
- **REQ-VER-06.** Cada VCF de verdad debe incluir en su propia cabecera o
  metadatos estándar la información necesaria para identificar sin ambigüedad la
  referencia usada y el conjunto al que pertenece. Debe cumplir la versión
  oficial de VCF que se elija durante la planificación.
- **REQ-VER-07.** Los genotipos de todas las clases de variante deben expresarse
  respecto a la referencia usada por el conjunto correspondiente.
- **REQ-VER-08.** Cuando se generen varios conjuntos, cada VCF de verdad debe
  contener las mismas muestras biológicas, identificadas de forma estable, pero
  sus registros de variante deben considerarse independientes de los registros
  de cualquier otra referencia.
- **REQ-VER-09.** Cada VCF debe construirse directamente contra la secuencia y
  las coordenadas de su referencia de mapeo. No debe obtenerse copiando,
  traduciendo ni transformando registros de un VCF generado contra otra
  referencia.
- **REQ-VER-10.** Cada conjunto por referencia debe tener su propio espacio de
  identificadores de variantes. Ni el VCF, ni su manifiesto ni su fichero
  complementario deben exponer un identificador global de variante o una tabla
  de correspondencias registro a registro entre referencias. El manifiesto
  general puede conservar el linaje entre referencias completas, sin convertir
  ese linaje en enlaces entre variantes.
- **REQ-VER-11.** Cuando un muestreo evolutivo solicite VCF de verdad, debe
  generarse uno independiente para cada referencia de mapeo seleccionada, con
  los genotipos de los mismos individuos muestreados y las reglas de coherencia
  establecidas en REQ-VER-03 a REQ-VER-10.
- **REQ-VER-12.** Los identificadores de muestra del VCF de verdad deben
  corresponder a los identificadores estables de los individuos simulados.
- **REQ-VER-13.** Los genotipos de los descendientes en el VCF de verdad deben
  reflejar los eventos de recombinación realizados y mantener la trazabilidad a
  los haplotipos o copias cromosómicas parentales de los que proceden.
- **REQ-VER-14.** Cuando una población intermedia de cruces solicite VCF de
  verdad, debe generarse uno independiente para cada referencia de mapeo
  seleccionada, con los genotipos de todos sus individuos y las mismas reglas de
  coherencia que los demás conjuntos.
- **REQ-VER-15.** Los VCF de verdad deben expresar genotipos faseados siempre que
  el formato permita representar el caso simulado. Si una limitación de
  representación impide conservar la fase completa dentro del propio registro
  VCF, el sistema debe identificar el caso y conservarla sin pérdida mediante
  el complemento estructurado. No debe convertirla silenciosamente en genotipos
  no faseados ni perderla de la verdad enlazada.
- **REQ-VER-16.** La clasificación de inserciones y deleciones en cada VCF de
  verdad debe usar el mismo umbral efectivo aplicado durante la generación y las
  estadísticas. La versión exacta y las reglas técnicas de normalización todavía
  no están definidas.
- **REQ-VER-17.** Los SNP y *indels* deben representarse mediante alelos
  explícitos REF y ALT, normalizados respecto a la referencia de mapeo del
  conjunto. REF es el alelo de referencia y ALT el alelo alternativo.
- **REQ-VER-18.** Las variantes estructurales deben usar alelos simbólicos cuando
  corresponda y los campos estándar SVTYPE, END y SVLEN para indicar el tipo, la
  posición final y la diferencia de longitud.
- **REQ-VER-19.** Una variante compleja debe conservar información suficiente
  para reconstruir el cambio real simulado. Ninguna limitación de representación
  debe provocar una pérdida silenciosa de la verdad.
- **REQ-VER-20.** Los alelos de repeticiones en tándem deben usar REF y ALT
  explícitos cuando su tamaño sea manejable y una representación simbólica para
  tramos largos o complejos.
- **REQ-VER-21.** El VCF de verdad debe reflejar los cambios de copias y, cuando
  estén activadas, las mutaciones internas de interrupciones, pureza o motivo en
  el haplotipo afectado.
- **REQ-VER-22.** En ambas representaciones, el propio VCF debe exponer sin
  pérdida el motivo, el número de copias de cada alelo o haplotipo, la longitud,
  la pureza y las interrupciones, mantener la fase y permitir reconstruir el
  genotipo real simulado. Esos detalles no deben conservarse únicamente en el
  manifiesto.
- **REQ-VER-23.** La cabecera del VCF debe documentar los campos estándar o
  adicionales empleados para representar las repeticiones en tándem. Los campos
  concretos se elegirán durante la planificación de la implementación.
- **REQ-VER-24.** Si cualquier variante o genotipo complejo no puede
  representarse completamente y sin pérdida en VCF, el sistema debe generar un
  fichero estructurado complementario de verdad. No debe generarlo cuando todos
  los casos sean plenamente representables, salvo que una opción futura lo
  solicite expresamente.
- **REQ-VER-25.** Cada caso conservado en el complemento debe tener un
  identificador estable compartido con su registro VCF. El VCF debe referenciar
  ese identificador y el fichero complementario correspondiente.
- **REQ-VER-26.** El complemento debe conservar toda la información necesaria
  para reconstruir la secuencia o cambio real, el haplotipo y la fase, las
  coordenadas y el genotipo simulado.
- **REQ-VER-27.** El VCF o, cuando corresponda, el fichero de verdad obligatorio
  debe conservar la dosis alélica para la ploidía efectiva de 1 a 10 copias.
- **REQ-VER-28.** Por cada conjunto definido por una referencia de mapeo, el
  sistema debe generar un único VCF de verdad multimuesta con todos los
  individuos del conjunto como columnas. Sus identificadores deben coincidir
  exactamente con los usados en BAM o CRAM, grupos de lectura, genotipos
  completos y manifiesto.
- **REQ-VER-29.** Un VCF de verdad no debe incluir individuos alineados contra
  otra referencia ni mezclar conjuntos. Los VCF individuales no deben formar
  parte de la salida principal.
- **REQ-VER-30.** El VCF multimuesta de verdad debe incluir únicamente loci en
  los que al menos una muestra tenga un alelo no de referencia respecto a la
  referencia de mapeo del conjunto. No debe incluir registros para posiciones
  invariantes ni bloques gVCF de forma predeterminada.
- **REQ-VER-31.** El fichero estructurado obligatorio de genotipos y haplotipos
  debe seguir siendo la fuente completa de verdad, con independencia de los loci
  omitidos del VCF.
- **REQ-VER-32.** El sistema debe construir cada VCF de verdad directamente a
  partir del estado genotípico faseado producido por la simulación. No debe
  llamar, imputar, fasear, corregir ni inferir genotipos a partir de FASTQ, BAM o
  CRAM para generar ese VCF.
- **REQ-VER-33.** En cada locus incluido, el VCF de verdad debe conservar el
  genotipo completo y la fase real de cada individuo siempre que el formato
  pueda representarlos. Una cobertura nula o baja, la contaminación, los
  errores de secuenciación, las calidades bajas, la ausencia de alineamientos o
  un mapeo incorrecto, múltiple o ambiguo no deben convertir un genotipo
  verdadero en ausente ni modificarlo.
- **REQ-VER-34.** Cuando una limitación del formato impida representar sin
  pérdida un genotipo o su fase, debe aplicarse el complemento estructurado de
  verdad definido en REQ-VER-24 a REQ-VER-26. Esta excepción de representación
  no debe confundirse con ausencia de evidencia en las lecturas o
  alineamientos.
- **REQ-VER-35.** Las degradaciones y anomalías configuradas deben afectar solo
  a las lecturas y los alineamientos que correspondan. No deben modificar el
  estado genotípico simulado, el VCF de verdad ni su complemento.
- **REQ-VER-36.** El campo `GT` del VCF de verdad debe contener el genotipo real
  simulado de cada individuo, con la fase y el número de copias que correspondan
  a su ploidía efectiva.
- **REQ-VER-37.** El VCF de verdad no debe incluir campos `FORMAT` derivados de
  observaciones o inferencias, incluidos `DP`, `AD`, `GQ` y `PL`. Estos campos
  describen, respectivamente, profundidad, soporte por alelo, calidad del
  genotipo y verosimilitudes genotípicas. El VCF solo debe añadir a `GT` los
  metadatos biológicos o estructurales necesarios para representar la variante y
  enlazarla con la verdad completa.
- **REQ-VER-38.** La cobertura, el soporte alélico, las calidades y las demás
  métricas observadas en FASTQ, BAM o CRAM deben conservarse en los resultados
  de validación correspondientes y no en el VCF de verdad.
- **REQ-VER-39.** Todos los registros válidos del VCF de verdad deben usar
  `QUAL=.` para indicar que el archivo no asigna una estimación probabilística de
  calidad a la variante verdadera.
- **REQ-VER-40.** Todos los registros válidos del VCF de verdad deben usar
  `FILTER=PASS`. Un registro que no pueda representarse correctamente mediante
  el VCF y, cuando proceda, su complemento enlazado, o que no pueda validarse, no
  debe emitirse como `PASS`: debe tratarse como un error estructural, detener la
  creación del conjunto e impedir que este se marque como completo.
- **REQ-VER-41.** La cabecera de cada VCF de verdad debe documentar que `QUAL=.`
  expresa la ausencia intencionada de una calidad probabilística y que
  `FILTER=PASS` identifica un registro verdadero que ha superado las
  comprobaciones estructurales aplicables.
- **REQ-VER-42.** Las translocaciones deben representarse mediante registros BND
  o alelos simbólicos cuando permitan expresar el evento sin pérdida. Cuando no
  sea suficiente, el VCF debe enlazar un complemento de verdad que conserve el
  evento compuesto, sus puntos de ruptura, segmentos y orientaciones.
- **REQ-VER-43.** El VCF debe incluir toda diferencia respecto a la referencia
  del conjunto cuando al menos una muestra tenga un alelo alternativo, aunque
  ese alelo esté fijado en todas las muestras. Una variante respecto a la
  referencia no necesita ser polimórfica entre las muestras.
- **REQ-VER-44.** Cuando un mismo locus tenga varios alelos alternativos que el
  VCF pueda representar sin pérdida, deben reunirse en un único registro
  multialélico. `GT` debe conservar para cada muestra la ploidía, la dosis y la
  fase de esos alelos después de la normalización.
- **REQ-VER-45.** Si los alelos de un locus no pueden reunirse en un único
  registro normalizado sin pérdida, el sistema debe aplicar la representación
  VCF válida y el complemento de verdad exigido por REQ-VER-24 a REQ-VER-26. En
  variantes pequeñas solapadas, los registros deben derivarse de la secuencia
  haplotípica final y no incluir eventos cuyo efecto haya sido eliminado. Los
  alelos finales de haplotipos diferentes pueden coexistir como ALT distintos.
  Varios eventos solapados dentro del mismo haplotipo que formen un único alelo
  final deben representarse como una sola variante o ALT compleja normalizada.
  Si el VCF no basta, debe usarse el complemento sin pérdida.
- **REQ-VER-46.** Si dos eventos usan los mismos individuos biológicos, la misma
  referencia y la misma verdad genotípica, sus VCF deben contener los mismos
  identificadores locales de variante, coordenadas, REF, ALT y genotipos. Las
  diferencias de subsemilla, tecnología, FASTQ, BAM, CRAM o degradaciones
  técnicas no deben cambiar el VCF de verdad.

Cuando se genera, el VCF de verdad recoge los genotipos definidos por la
simulación, no una estimación basada en la evidencia observada. Su finalidad
dentro de este proyecto es conservar la descripción exacta de lo simulado y
mantener la coherencia entre las salidas. Si se seleccionan varias referencias,
todas las muestras forman parte de todos los conjuntos, pero los archivos de
esos conjuntos deben permanecer separados.

### 4.6. Validación interna de lecturas y alineamientos

En esta sección, el soporte alélico es la cantidad o proporción de lecturas que
contienen un alelo concreto. Una tolerancia estocástica es el intervalo aceptable
de variación causado por el azar de la simulación; no exige que una tasa o una
cobertura realizada coincidan exactamente con su objetivo.

- **REQ-VAL-01.** Siempre que se generen FASTQ o BAM o CRAM, el sistema debe
  comprobar internamente que reflejan de forma estadísticamente coherente la
  verdad biológica simulada y los artefactos configurados.
- **REQ-VAL-02.** La validación debe usar el genoma, los haplotipos, la ploidía,
  los genotipos y la procedencia verdadera de las lecturas definidos por la
  simulación. No debe ejecutar un *SNP caller* ni llamar, imputar, fasear,
  corregir o inferir genotipos a partir de BAM o CRAM.
- **REQ-VAL-03.** La validación debe comprobar que las lecturas proceden de los
  individuos y haplotipos correctos, salvo las procedencias diferentes que la
  configuración introduzca expresamente como contaminación u otro artefacto.
- **REQ-VAL-04.** La validación debe comparar la cobertura objetivo con la
  cobertura realizada en los ámbitos aplicables, incluidas las sustituciones por
  individuo, cromosoma o *scaffold* y región.
- **REQ-VAL-05.** La validación debe comprobar que el soporte alélico observado
  es compatible, dentro de tolerancias estocásticas, con la dosis alélica y la
  ploidía simuladas y con los efectos configurados que puedan alterarlo. No debe
  exigir un balance exacto entre alelos o copias.
- **REQ-VAL-06.** La validación debe comparar los valores objetivo y realizados
  de errores de secuenciación, contaminación biológica y técnica, duplicados y
  las demás degradaciones activadas, usando las tolerancias que correspondan.
- **REQ-VAL-07.** Para datos *paired-end*, la validación debe comprobar la
  coherencia entre ambos miembros de cada pareja, incluidos origen, orientación,
  distancia, estado de alineamiento y anomalías configuradas.
- **REQ-VAL-08.** La validación debe comprobar la coherencia de los
  alineamientos con su referencia, sus posiciones verdaderas cuando proceda,
  sus campos y banderas y los grupos de lectura declarados. Debe distinguir los
  resultados esperados del mapeador de las incoherencias internas del resultado.
- **REQ-VAL-09.** Toda desviación intencionada respecto a la verdad biológica o
  al comportamiento ordinario debe poder explicarse mediante un parámetro,
  perfil, artefacto o decisión aleatoria registrada de la simulación.
- **REQ-VAL-10.** La validación debe producir uno o varios resultados
  estructurados que indiquen las comprobaciones realizadas, los objetivos, los
  valores observados, las tolerancias aplicadas, las desviaciones detectadas y
  su explicación cuando exista. El manifiesto debe enlazar estos resultados con
  los FASTQ, BAM o CRAM, muestras, grupos y conjuntos a los que se refieren.
- **REQ-VAL-11.** Toda inconsistencia estructural o lógica debe tratarse como un
  error que detiene la ejecución e impide marcar como completo el conjunto
  afectado. Esto incluye, como mínimo, archivos inválidos, referencias
  incorrectas, identificadores incompatibles, genotipos imposibles y cabeceras o
  índices incoherentes.
- **REQ-VAL-12.** Una desviación estadística fuera de las tolerancias
  configuradas debe producir una advertencia de forma predeterminada. La
  configuración debe permitir elevar estas desviaciones a errores que detengan
  la ejecución e impidan marcar como completo el conjunto afectado. Se exceptúa
  las desviaciones entre N objetivo y N realizado y entre Ne objetivo y Ne
  realizado: deben advertirse, pero nunca elevarse a error ni detener la
  simulación por sí mismas. El caso N realizado igual a 0 aplica la política de
  REQ-CONF-104 y no debe confundirse con elevar una desviación estadística.
- **REQ-VAL-13.** Para cada comprobación estadística, el resultado estructurado y
  el manifiesto deben registrar la tolerancia, el valor esperado, el valor
  realizado y la decisión aplicada: conforme, advertencia o error. Para N y Ne
  deben registrar la desviación cuantificada y la advertencia no invalidante,
  sin aplicar una política que pueda convertirla en error. N realizado igual a
  0 se resuelve además conforme a REQ-CONF-104.
- **REQ-VAL-14.** La validación WGS debe comprobar el objetivo excluyente de cada
  biblioteca, la cobertura calculada respecto a la longitud haploide, la
  distribución aproximadamente Poisson anterior a sesgos y el reparto realizado
  por cromosoma y copia o haplotipo.
- **REQ-VAL-15.** La validación debe distinguir una variación aleatoria admisible
  del soporte alélico de un modelo intencionado de desequilibrio o pérdida
  alélica, que no existe en la primera versión.
- **REQ-VAL-16.** La validación de contaminación debe comprobar límites,
  compatibilidad, matriz o reparto automático, sustituciones biológicas,
  transferencias técnicas, descartes, totales por tanda y cobertura final.
  También debe comprobar la conversión conjunta de unidades, el carácter óptimo de
  la asignación entera, sus desviaciones y la resolución reproducible de
  empates conforme a REQ-CONT-15 y REQ-CONT-16.
- **REQ-VAL-17.** La validación poblacional debe comprobar N objetivo, N
  realizado, su desviación y su aviso, además de Ne objetivo o su ausencia y Ne
  realizado en cada subpoblación y generación. Debe
  comprobar que Ne realizado procede solo de los progenitores con descendencia
  y de la desigualdad o varianza de sus contribuciones. Cuando haya objetivo,
  debe cuantificar cualquier desviación y emitir una advertencia no invalidante.
  También debe comprobar autofecundaciones esperadas y realizadas, contribución
  familiar, rutas de migración, supervivencia y extinciones con las
  tolerancias y políticas acordadas. Para el locus seleccionado debe validar
  todos los genotipos o combinaciones de dosis bialélicos o multialélicos, sus
  valores genéticos, valores del carácter, aptitudes y probabilidades de
  supervivencia. Debe comprobar que el valor del carácter coincide con el valor
  genético y que no se aplicaron ruido ambiental ni heredabilidad. También debe
  comprobar el sentido direccional configurado y que las aptitudes respetan el
  orden esperado; una función estabilizadora o disruptiva debe rechazarse.
- **REQ-VAL-18.** La validación de mutación y recombinación debe comprobar tasas
  y eventos esperados o realizados, orden recombinación→mutación, mapas
  completos, intervalos recombinantes válidos y emparejamientos bivalentes
  reproducibles. Debe comprobar también el mínimo de un entrecruzamiento en cada
  cromosoma y en los *scaffolds* que superen su umbral efectivo.
- **REQ-VAL-19.** Para cada individuo muestreado, la validación debe comprobar
  que los segmentos reconstruyen sin interrupciones su genoma, que cada segmento
  alcanza una copia haplotípica fundadora y que las mutaciones de novo se
  conservan como capas separadas sin borrar ese origen. Cuando existan cambios
  estructurales, debe comprobar además la correspondencia de coordenadas y la
  orientación entre los segmentos actuales y fundadores que sean relacionables.
- **REQ-VAL-20.** Cuando la población inicial use el modelo de equilibrio de
  Hardy-Weinberg, la validación debe comparar las frecuencias y recuentos
  genotípicos y alélicos esperados y realizados según la ploidía. Debe tratar la
  variación finita como estocástica y no exigir recuentos genotípicos exactos.
  En modo controlado debe comprobar por separado los recuentos alélicos fijados
  y la distribución genotípica aleatoria condicionada.
- **REQ-VAL-21.** La validación de la población inicial debe comprobar si se usó
  asignación independiente o bloques con LD y que los haplotipos resultantes
  conservan la fase. Cuando se configure LD, debe comparar las frecuencias
  conjuntas objetivo y realizadas con las tolerancias que se definan, derivar de
  ellas las frecuencias alélicas y rechazar cualquier frecuencia alélica aportada
  que discrepe fuera de tolerancia. En modo controlado debe comprobar además que
  los recuentos haplotípicos suman C y que su distribución entre individuos
  respeta el modelo genotípico o las plantillas explícitas.
- **REQ-VAL-22.** La validación inicial debe comprobar por subpoblación que no se
  mezclaron los modos probabilístico y controlado dentro del mismo locus o
  bloque, que las frecuencias suman 1, los recuentos genotípicos suman N y los
  recuentos alélicos o haplotípicos suman C. También debe comprobar que las
  frecuencias alélicas derivadas son compatibles con los genotipos, dosis o
  haplotipos aportados.
- **REQ-VAL-23.** La validación debe comprobar que F se aplicó solo en su ámbito
  global o de subpoblación, que no se usó con ploidía 1 ni en ámbitos locales y
  que toda distribución genotípica explícita prevaleció en su locus o bloque.
  Debe comprobar que F=0 conserva las probabilidades de Hardy-Weinberg, que la
  homocigosis crece al aumentar F y que F=1 produce homocigosis completa. Por
  subpoblación y generación registrada, debe comparar la heterocigosidad
  observada y esperada según las frecuencias alélicas para obtener F realizado.
  En poliploides, debe considerar heterocigoto al individuo cuyas copias no sean
  todas del mismo alelo. Un F realizado negativo debe conservarse y no marcarse
  como error.
- **REQ-VAL-24.** Cada evento de muestreo debe validar antes de generar datos que
  sus individuos explícitos y aleatorios existen, pertenecen a la generación y
  estrato correctos, no están duplicados y son suficientes. Una solicitud mayor
  que el conjunto elegible debe detener el evento sin reemplazo ni reducción
  silenciosa del tamaño.
- **REQ-VAL-25.** La validación técnica debe comprobar que cada evento usa una
  sola tecnología y un perfil cuando genera FASTQ, BAM o CRAM, y un solo
  mapeador cuando genera alineamientos. Debe comprobar que sus FASTQ no proceden
  de otro evento y que, dentro del evento, los mismos FASTQ y mapeador se
  aplicaron a todas las referencias seleccionadas. Un evento que solo solicite
  verdad genotípica o VCF debe validarse sin exigirle parámetros de
  secuenciación.
- **REQ-VAL-26.** La validación debe comprobar la unicidad global y coherencia de
  cada grupo de lectura, que distingue evento, biblioteca, tanda y carril, y que
  `SM` y las columnas del VCF conservan el identificador biológico estable.
  También debe comprobar la unicidad de cada instancia de lectura dentro del
  evento y la coherencia de sus vínculos de origen, copias biológicas,
  transferencias técnicas y duplicados configurados.
- **REQ-VAL-27.** La validación estructural debe comprobar el convenio de
  coordenadas, los límites terminales, los puntos de ruptura y orientaciones de
  cada variante estructural, y el límite acumulado de copias locales por copia
  haplotípica. En duplicaciones debe comprobar la excepción exclusiva para una
  fuente que abarque la secuencia completa, así como las coordenadas y
  correspondencias de cada estado anterior y posterior. También debe validar
  secuencialmente que el objetivo de cada evento explícito ordenado sigue
  existiendo antes de aplicarlo.
- **REQ-VAL-28.** Para cada VCF debe comprobarse que se incluyeron los ALT
  fijados, que los ALT representables de un locus se reunieron en un registro
  multialélico y que ploidía, dosis y fase se conservaron. También debe comprobar
  que los espacios de identificadores no crean enlaces entre referencias.
- **REQ-VAL-29.** Cuando varios eventos compartan individuos, referencia y verdad
  genotípica, la validación debe comprobar la identidad de identificadores
  locales de variante, coordenadas, REF, ALT y genotipos entre sus VCF, sin
  exigir identidad de FASTQ, BAM, CRAM o degradaciones.
- **REQ-VAL-30.** La validación de familias debe comprobar la cardinalidad y el
  tipo de cada evento, la simetría de los parentales de un cruzamiento, la
  compatibilidad estructural reproductiva y la trazabilidad de las mutaciones de
  novo cuando estén activadas. Para seleccionar descendientes debe comprobar los
  criterios sobre la verdad biológica, las combinaciones lógicas, el ámbito
  elegible y la cantidad exacta o proporción solicitada. También debe comprobar
  el sorteo sin reemplazo, el uso reproducible de la semilla y la política de
  error cuando falten elegibles para un número exacto. Para parentales expresos
  debe comprobar la referencia, la asignación de cada cambio a copias concretas
  y la reconstrucción completa de la fase sin imputación ni sorteo. En grafos
  ramificados debe comprobar que la reutilización conserva una sola identidad y
  que cada evento produce gametos independientes y descendientes nuevos. Debe
  distinguir una dependencia cíclica imposible de la autofecundación, el
  retrocruzamiento y la reutilización válida de un parental anterior. También
  debe comprobar la unicidad de identificadores dentro de cada categoría, la
  existencia y precedencia de todas las referencias, la creación única de cada
  descendiente, la cardinalidad parental, la compatibilidad de referencia,
  ploidía y estructura y que todas las reglas y salidas solicitadas sean
  producibles.
- **REQ-VAL-31.** Para cada frecuencia alélica objetivo, la validación debe
  comprobar que se usó un solo modo. En el probabilístico debe comparar valores
  esperados y realizados sin aplicar una tolerancia de recuento controlado. En
  el controlado debe comprobar que el vector de recuentos suma C, minimiza la
  desviación global y cumple 0,5/C o la tolerancia más estricta configurada. Debe
  comprobar también que todo empate y la distribución condicionada entre
  individuos se resolvieron y registraron reproduciblemente.
- **REQ-VAL-32.** Para cada evento de muestreo, la validación debe comprobar que
  existen por separado el resumen poblacional y el muestral, que ambos proceden
  de la generación correcta y que cada estadística usa solo los individuos de
  su ámbito declarado. Debe comprobar frecuencias alélicas, heterocigosidad, F,
  N y Ne cuando correspondan, incluidas las reglas diferenciadas de F y Ne
  muestrales de REQ-EVO-65. Una medida muestral etiquetada o presentada como
  poblacional, un valor poblacional copiado como muestral o una omisión sin
  indicar que no es aplicable debe ser un error.
- **REQ-VAL-33.** Para cada cambio de N, la validación debe comprobar que el
  nuevo valor determinó el N objetivo de la generación siguiente y no modificó
  retrospectivamente la generación parental. Debe comprobar que un cuello de
  botella redujo ese objetivo y que una expansión lo aumentó respecto al valor
  anterior aplicable. También debe comprobar N realizado, desviación y aviso.
- **REQ-VAL-34.** Para cada división, la validación debe comprobar que todos los
  progenitores de la subpoblación origen se asignaron una sola vez, que ninguna
  hija comparte fundadores con otra y que cada hija pudo reproducirse conforme
  a sus parámetros. También debe comprobar la procedencia y el modo, regla o
  identificador que determinó cada asignación. Por hija debe comprobar los
  progenitores recibidos, N objetivo y N realizado, sin exigir que la suma de
  objetivos coincida con el N del origen.
- **REQ-VAL-35.** Para cada fusión, la validación debe comprobar que todos los
  progenitores de los orígenes pasaron una sola vez a la nueva subpoblación, que
  los orígenes dejaron de existir y que se conservaron las procedencias
  individuales y haplotípicas. Debe comprobar los parámetros propios de la
  población fusionada y que, sin restricciones explícitas, el origen no limitó
  los cruzamientos. Debe comprobar por separado el tamaño del conjunto
  reproductor reunido, N objetivo y N realizado. También debe comprobar la
  exclusión entre los modos probabilístico y controlado. En el primero debe
  comprobar la contribución esperada desde progenitores elegibles o pesos y la
  variación realizada. En el segundo debe comprobar los recuentos parentales
  exactos y su compatibilidad con la descendencia. En ambos debe comprobar la
  contribución genética realizada y su trazabilidad por haplotipos. Debe validar
  por separado la ascendencia genómica realizada, sin exigir que coincida
  exactamente con la proporción de contribuciones parentales.
- **REQ-VAL-36.** Para cada fundación, la validación debe comprobar que los
  fundadores se transfirieron sin reemplazo, se excluyeron del conjunto
  reproductor de origen y solo se reprodujeron en la nueva subpoblación durante
  esa generación. Debe comprobar que la nueva población descendió de ellos, que
  la original continuó con los individuos restantes y que se conservó la
  procedencia. Debe comprobar por separado el tamaño del conjunto fundador, N
  objetivo y N realizado de la nueva subpoblación.
- **REQ-VAL-37.** Para cada extinción, la validación debe comprobar que no se
  produjo descendencia desde la generación indicada, que se conservaron el
  historial y las salidas previas. Para una extinción planificada debe comprobar
  que no existieron dependencias posteriores inválidas. Toda recolonización debe
  usar un identificador nuevo y mantener separadas las historias de la
  subpoblación extinguida y la nueva. Para N realizado igual a 0 debe comprobar
  el valor global o sustituido de la política en la subpoblación y fase: si
  estaba habilitada debe existir una extinción por resultado; si estaba
  deshabilitada debe existir un error que detuvo inmediatamente toda la
  simulación, sin omisiones ni continuación de otras ramas. Nunca debe quedar una
  subpoblación activa vacía. También debe comprobar que todo evento dependiente
  de una extinción por resultado se omitió de forma atómica, sin usar solo las
  subpoblaciones restantes, y que se registraron todos sus requisitos y los que
  faltaban. Debe comprobar la propagación transitiva, la cadena auditable de
  causas y la continuidad exclusiva de las ramas independientes.
- **REQ-VAL-38.** Para cada fase con migración, la validación debe comprobar que
  se usó exclusivamente el modo proporcional o el de recuento exacto. Debe
  comprobar por origen que las rutas se asignaron conjuntamente, sin reemplazo
  y desde las ubicaciones iniciales de la generación,
  que ningún individuo tuvo dos destinos y que todos los migrantes salieron del
  origen y entraron en el destino antes de la selección y la reproducción. En
  modo proporcional debe comprobar la tasa solicitada, el número realizado y la
  proporción realizada; en modo exacto, la igualdad con el recuento solicitado.
  También debe validar reglas, identificadores y disponibilidad de elegibles.
  Debe rechazar cualquier ruta encadenada o registro que mantenga al mismo
  individuo como migrante en una generación posterior. Debe comprobar además
  los tamaños reproductores antes y después de migrar, el N objetivo aplicado a
  la descendencia y que la migración no modificó ese objetivo por sí sola.
- **REQ-VAL-39.** Para cada generación con varios eventos demográficos, la
  validación debe comprobar que existió un orden explícito único, que cada evento
  recibió el estado producido por el anterior y que se respetaron todas las
  dependencias. Un orden ausente, ambiguo o imposible debe constar como rechazo
  de configuración anterior a la simulación.
- **REQ-VAL-40.** Para cada conjunto de frecuencias o probabilidades que deba
  sumar 1, la validación debe comprobar la suma original, la tolerancia general
  efectiva y la decisión de normalizar o rechazar. Cuando se normalice, debe
  comprobar que los valores finales suman 1 dentro de la precisión aplicable y
  que se conservaron los valores anteriores y posteriores. Debe comprobar por
  separado las tolerancias específicas con otra finalidad, como `0,5/C` en el
  modo controlado.
- **REQ-VAL-41.** La validación debe comprobar que división, fundación,
  migración y selección de descendientes usan exclusivamente el catálogo
  declarativo común. Para cada evento debe comprobar identificadores, selección
  aleatoria, pertenencia, condiciones genéticas, combinaciones Y/O/NO, cantidad
  exacta o proporción, elegibles, seleccionados y resultado realizado según
  proceda. Debe rechazar scripts, código, funciones personalizadas y reglas no
  reconocidas, además de cualquier selección basada en evidencia técnica
  posterior.
- **REQ-VAL-42.** La validación del resultado mínimo debe comprobar por vía el
  alcance de los estados genotípicos completos de REQ-SAL-02 y distinguirlos de
  los historiales resumidos. La ausencia de un estado completo obligatorio o la
  presentación de un resumen como estado completo debe ser un error.
- **REQ-VAL-43.** La validación de estados debe comprobar la jerarquía de
  ejecución, evento, conjunto y artefacto de REQ-SAL-21 a REQ-SAL-26. Ninguna
  unidad debe marcarse como completa si una salida obligatoria subordinada está
  incompleta o fallida, ni debe publicarse o exportarse un artefacto parcial
  como válido.
- **REQ-VAL-44.** Toda afirmación de reproducibilidad debe validarse y registrarse
  en el nivel correspondiente de REQ-CONF-111 a REQ-CONF-115. No debe declararse
  identidad binaria a partir de una comprobación limitada a igualdad biológica
  o semántica.
- **REQ-VAL-45.** La validación de un FASTA debe interpretar como identificador
  el primer elemento de la cabecera, comprobar su unicidad y conservar la
  descripción completa. Debe rechazar cualquier símbolo de secuencia distinto
  de A, C, G, T y N, sin distinguir mayúsculas de minúsculas.

### 4.7. Manifiestos de la ejecución y de sus eventos

Los manifiestos son ficheros estructurados de salida que registran lo que el
sistema utilizó y generó. La raíz contiene el manifiesto general de la ejecución
y cada evento explícito tiene además el suyo, aunque solo materialice verdad o
VCF y no constituya un experimento de secuenciación.

- **REQ-MAN-01.** El sistema debe generar un manifiesto general para cada
  ejecución y un manifiesto adicional para cada evento explícito.
- **REQ-MAN-02.** El manifiesto debe ser un fichero estructurado y legible por el
  programa.
- **REQ-MAN-03.** El manifiesto debe identificar cuál de las tres vías de
  creación se utilizó.
- **REQ-MAN-04.** El manifiesto debe registrar el nombre o la ruta del fichero
  general y de cada fichero de configuración poblacional utilizado, incluida la
  ruta de la copia guardada en el directorio de resultados.
- **REQ-MAN-05.** El manifiesto debe conservar los parámetros realmente
  aplicados, incluida la ploidía global y las sustituciones efectivas por
  cromosoma o *scaffold*, el modo de creación de los cambios y los
  parámetros propios de la vía utilizada. También debe registrar la semilla
  usada, si fue suministrada o generada automáticamente, y la versión del
  programa.
- **REQ-MAN-06.** El manifiesto debe indicar si la referencia inicial fue
  simulada o suministrada. Si fue suministrada, debe registrar la procedencia
  original del FASTA y el nombre o la ruta de la copia utilizada desde el
  directorio de resultados.
- **REQ-MAN-07.** El manifiesto debe registrar cada referencia, su referencia de
  origen cuando sea derivada y los cambios aplicados para obtenerla.
- **REQ-MAN-08.** El manifiesto debe registrar los nombres o las rutas de los
  archivos generados.
- **REQ-MAN-09.** El manifiesto debe registrar las referencias de mapeo
  seleccionadas e identificar el conjunto independiente generado para cada una.
- **REQ-MAN-10.** El manifiesto debe distinguir los valores globales de los
  parámetros de generación aleatoria, los valores específicos definidos para
  individuos concretos y el valor efectivo aplicado a cada individuo afectado
  por una sustitución. También debe distinguir los valores definidos para el
  genoma, el cromosoma o *scaffold* y la región, y el valor efectivo resultante
  de aplicar la precedencia espacial y por individuo. Para cada clase de
  variante debe registrar si se utilizó una cantidad media esperada o una
  densidad y, en el segundo caso, si la densidad se expresó por base, kilobase o
  megabase. También debe registrar si se aplicó una distribución espacial
  uniforme o se configuraron regiones, incluidas las regiones con densidad
  cero.
- **REQ-MAN-11.** Para cada clase de variante, ámbito general o región e
  individuo al que corresponda, el manifiesto debe registrar la cantidad media
  esperada y el número realizado de eventos, junto con la distribución
  probabilística y los parámetros realmente utilizados.
- **REQ-MAN-12.** Para cada región configurada, el manifiesto debe registrar la
  referencia estructural, el cromosoma o *scaffold* y el intervalo indicados.
  Esta procedencia no debe convertirse en una correspondencia entre registros
  de VCF generados contra referencias diferentes.
- **REQ-MAN-13.** Para cada objetivo de frecuencia alélica poblacional, el
  manifiesto debe registrar la frecuencia objetivo, las restricciones de
  presencia, ausencia, genotipo o dosis fijadas para individuos concretos y la
  frecuencia realizada. Debe identificar el modo probabilístico o controlado.
  En el primero debe distinguir valores esperados y realizados; en el segundo,
  C, vector de objetivos y recuentos, desviación global, tolerancia, modelo de
  distribución entre individuos y resultado. Si hubo empate, debe registrar las
  soluciones candidatas, la elección y la semilla o flujo utilizado.
- **REQ-MAN-14.** Para una ejecución evolutiva, el manifiesto debe registrar las
  fases en su orden de ejecución, su duración en generaciones y los parámetros
  efectivos aplicados en cada fase, incluidas sus sustituciones temporales.
- **REQ-MAN-15.** Para cada muestreo evolutivo, el manifiesto debe registrar el
  evento, la generación, el número o la proporción solicitados, los individuos
  seleccionados y la semilla o el flujo aleatorio que corresponda. También debe
  indicar el estado genotípico y los metadatos conservados, si se solicitaron
  salidas adicionales y cuáles fueron los BAM, CRAM o VCF generados.
- **REQ-MAN-16.** Para cada individuo, el manifiesto debe relacionar su
  identificador estable, su generación de origen, sus parentales cuando existan,
  sus observaciones y los archivos generados que le correspondan. Debe distinguir
  los identificadores de individuo, evento de muestreo y conjunto.
- **REQ-MAN-17.** En una ejecución evolutiva, el manifiesto debe registrar las
  trayectorias configuradas de N objetivo y Ne objetivo, así como N objetivo, N
  realizado y su desviación, Ne objetivo o su ausencia, Ne realizado y su
  desviación por generación, fase o tramo. Debe hacerlo para la población
  completa y para cada subpoblación cuando corresponda.
- **REQ-MAN-18.** El manifiesto debe registrar la estructura de subpoblaciones, sus
  identificadores y parámetros, la pertenencia de los individuos por generación
  u observación y las migraciones direccionales configuradas y realizadas. Para
  cada fase debe registrar el modo proporcional o exacto y, por ruta, el método
  de selección, los valores solicitados y realizados y los individuos migrantes.
  Por subpoblación debe registrar los tamaños reproductores antes y después de
  migrar y el N objetivo de su descendencia.
- **REQ-MAN-19.** Para cada división, fusión, fundación o extinción, el manifiesto
  debe registrar el evento configurado, la generación, las subpoblaciones e
  individuos implicados, las relaciones de origen y destino y el resultado
  realizado. Cuando coincidan varios eventos demográficos, debe registrar su
  orden explícito y el estado recibido y producido por cada uno. En una división
  debe conservar el reparto completo de fundadores y su modo o regla y, por
  hija, los progenitores recibidos, N objetivo y N
  realizado. En una fusión debe conservar la procedencia individual y
  haplotípica, el fin de los orígenes, el tamaño del conjunto reproductor, N
  objetivo y N realizado, el modo de contribución, los pesos o recuentos
  parentales objetivo y realizados por origen y la ascendencia genómica
  realizada. En una fundación debe conservar los transferidos, los
  individuos restantes en el origen, el tamaño del conjunto fundador, N objetivo
  y N realizado de la nueva subpoblación. En una extinción debe
  enlazar el historial y las salidas previas y distinguir si fue planificada o
  causada por N realizado igual a 0. En este último caso debe registrar el valor
  efectivo de la política de extinción, su ámbito y si procedía del valor global
  o de una sustitución. Cada evento no ejecutado por una extinción por resultado
  debe conservar su definición completa, todas las subpoblaciones requeridas,
  las ausentes y la cadena de causas hasta la extinción original, incluido
  cualquier resultado intermedio que no llegó a crearse. Una recolonización debe enlazar la
  subpoblación nueva con el ámbito recolonizado sin unir sus identificadores ni
  sus historias.
- **REQ-MAN-20.** El manifiesto debe registrar qué procesos evolutivos estuvieron
  activos o inactivos y los parámetros y valores efectivos aplicados por fase,
  generación y subpoblación cuando corresponda.
- **REQ-MAN-21.** El manifiesto debe registrar el orden efectivo en que se
  aplicaron los procesos dentro de cada escenario o fase e indicar si procedía
  de la configuración o de la secuencia predeterminada.
- **REQ-MAN-22.** El manifiesto debe registrar si se usó una tasa uniforme o un
  mapa de recombinación, la referencia estructural, las tasas, el mapa y los
  parámetros efectivos y los resúmenes generacionales de recombinación. Para cada
  individuo y evento debe registrar el grupo aplicable, cuando exista, y el
  perfil general o específico que se utilizó efectivamente. El manifiesto debe
  conservar también la pertenencia única de cada individuo a su grupo de
  recombinación.
- **REQ-MAN-23.** Para cada parental de una familia o población segregante, el
  manifiesto debe distinguir si se referenció un individuo ya disponible o se
  definió en el fichero poblacional, y conservar su identificador, referencia,
  genotipo faseado por copia y linaje. Para un parental expreso debe conservar
  los cambios declarados, la copia o copias asignadas y que lo no declarado
  mantuvo el estado de referencia.
- **REQ-MAN-24.** Para cada grafo de cruces, el manifiesto debe registrar los
  eventos, su orden realizado, sus parentales o fuentes y sus descendientes,
  usando los identificadores estables de eventos e individuos. Si un individuo
  alimenta varias ramas, debe enlazar todos los eventos con la misma identidad y
  distinguir los gametos y descendientes nuevos producidos en cada uno.
- **REQ-MAN-25.** Para cada selección de descendientes como parentales
  posteriores, el manifiesto debe registrar la población candidata, los
  identificadores o las condiciones aplicadas, su combinación lógica, la
  semilla o flujo aleatorio cuando corresponda y los individuos elegibles y
  seleccionados. Debe distinguir número exacto de proporción. En el segundo
  caso debe registrar las proporciones solicitada y realizada.
- **REQ-MAN-26.** Para cada población intermedia de cruces, el manifiesto debe
  identificar el estado genotípico y los metadatos conservados de todos sus
  individuos, si se solicitaron salidas adicionales y cuáles fueron los BAM,
  CRAM o VCF generados.
- **REQ-MAN-27.** El manifiesto o el estado genotípico trazable debe conservar la
  fase, los haplotipos y el origen de los segmentos heredados o recombinados. El
  manifiesto debe registrar cada limitación de representación en un VCF, su causa
  y el complemento que conserva sin pérdida la fase y la verdad afectadas.
- **REQ-MAN-28.** El manifiesto debe registrar el umbral efectivo que distinguió
  *indels* de variantes estructurales en inserciones y deleciones.
- **REQ-MAN-29.** Para las repeticiones en tándem, el manifiesto debe registrar los loci
  y parámetros configurados, su procedencia concreta o aleatoria, el motivo,
  unidad, copias o longitud y los objetivos de pureza y variabilidad. También
  debe conservar los valores realizados de pureza, interrupciones y variabilidad
  por haplotipo, individuo o población cuando corresponda.
- **REQ-MAN-30.** Para cada mutación de repetición en tándem, el manifiesto debe
  registrar el modelo, los parámetros efectivos, el locus y haplotipo afectados,
  las unidades ganadas o perdidas y, cuando estén activadas, las modificaciones
  internas de interrupciones, pureza o motivo.
- **REQ-MAN-31.** Cuando exista un fichero complementario de verdad, el
  manifiesto debe enlazarlo con el VCF, la referencia de mapeo y el conjunto
  correspondientes y conservar los identificadores compartidos de los casos
  dentro de ese conjunto. Esos identificadores son locales al conjunto y no
  deben enlazar variantes entre referencias.
- **REQ-MAN-32.** El manifiesto debe distinguir las salidas mínimas, las
  activadas de forma predeterminada, las solicitadas expresamente y las
  desactivadas, y registrar todos los archivos generados para satisfacerlas.
- **REQ-MAN-33.** El manifiesto debe registrar las estimaciones previas de
  volumen, almacenamiento, memoria y tiempo, las advertencias emitidas y, cuando
  pueda medirse, el consumo real de la ejecución. Para duplicaciones debe incluir
  el máximo local efectivo, si se elevó sobre el valor predeterminado, el
  presupuesto de recursos configurado y el resultado de compararlo con la
  estimación previa.
- **REQ-MAN-34.** Para cada conjunto y cada FASTQ, el manifiesto debe
  registrar si se usó el modo realista o el directo, la tecnología, el perfil,
  los parámetros efectivos y la semilla o el flujo aleatorio aplicable. En el
  modo directo debe etiquetar claramente el origen sintético.
- **REQ-MAN-35.** Para cada FASTQ final generado, el
  manifiesto debe registrar su procedencia, las muestras y lecturas que contiene
  y los BAM o CRAM para los que se usó. También debe relacionar cada conjunto
  con sus FASTQ y su referencia de mapeo.
- **REQ-MAN-36.** Para cada experimento de mapeo, el manifiesto debe registrar el
  nombre y la versión del único mapeador, el comando efectivo completo y sus
  parámetros, cada referencia y los FASTQ usados. Cuando se use minimap2, debe
  distinguir el *preset*, los valores predeterminados y las sustituciones de la
  configuración.
- **REQ-MAN-37.** Para la cobertura, el manifiesto debe registrar el objetivo y
  la cobertura realizada por ámbito y muestra cuando pueda medirse, incluidos
  los valores globales y las sustituciones efectivas por individuo, cromosoma o
  *scaffold* y región.
- **REQ-MAN-38.** Para cada fuente contaminante, el manifiesto debe registrar su
  identificador y origen, la fracción objetivo, la fracción realizada y los
  recuentos de lecturas o parejas propias retiradas y contaminantes incorporadas,
  además de la semilla o flujo aleatorio de selección. En datos *paired-end*,
  debe registrar las parejas completas como unidad.
- **REQ-MAN-39.** Para una fuente contaminante externa, el manifiesto debe
  identificarla como tal, registrar la ubicación original de cada FASTQ, la
  copia utilizada, el resultado de su validación y los datos necesarios para
  relacionar la fuente con las lecturas seleccionadas y la semilla o el flujo
  aleatorio aplicable. También debe incluir los datos de integridad que se
  definan durante la implementación.
- **REQ-MAN-40.** Para la contaminación técnica, el manifiesto debe registrar la
  tanda, las bibliotecas, las muestras y sus índices, la estrategia de índices,
  las tasas y procedencias configuradas, los errores simulados y las cantidades
  configuradas y realizadas de lecturas reasignadas y sin asignar o filtradas.
  Debe resumir la procedencia verdadera conservada en la verdad conocida.
- **REQ-MAN-41.** El manifiesto debe enlazar cada FASTQ final con su fichero de
  procedencia y registrar que no se conservó un FASTQ limpio adicional. Debe
  resumir las lecturas o parejas contaminantes presentes y los recuentos
  agregados de lecturas descartadas o sin asignar, desglosados al menos por
  muestra y por causa u origen conocido.
- **REQ-MAN-42.** Para los errores de secuenciación, el manifiesto debe registrar
  la tecnología, el nombre y la versión del perfil predeterminado, la tasa y los
  tipos de error representados, todas las sustituciones configuradas y los
  valores efectivos por individuo, cromosoma o *scaffold* y región.
- **REQ-MAN-43.** Cuando se calibre un perfil empírico, el manifiesto debe
  identificar el FASTQ real y su referencia, registrar los parámetros estimados,
  las limitaciones y advertencias de fiabilidad y si el perfil fue seleccionado
  como efectivo o solo estimado.
- **REQ-MAN-44.** El manifiesto debe identificar y enlazar el fichero de perfil
  empírico generado, su versión de esquema, la tecnología, las entradas, el
  método y las versiones de las herramientas, las medidas de calidad y las
  advertencias de la estimación.
- **REQ-MAN-45.** Para los duplicados de PCR y ópticos, el manifiesto debe
  registrar por separado el perfil, la tasa y los valores predeterminados o
  sustituidos realmente usados, las cantidades esperadas y realizadas y la
  relación trazable entre cada duplicado y su lectura o pareja de origen.
- **REQ-MAN-46.** El manifiesto debe registrar si los duplicados se entregaron
  marcados o sin marcar y, cuando se marcaron, la herramienta o el método
  efectivo utilizado.
- **REQ-MAN-47.** Para los restos de adaptadores, el manifiesto debe registrar si
  la degradación estuvo activa, los parámetros recomendados y efectivos, y las
  cantidades configuradas y realizadas de lecturas afectadas.
- **REQ-MAN-48.** Para el sesgo de contenido GC, el manifiesto debe registrar el
  perfil, la curva o relación, la intensidad, si fue anulado y la regla efectiva
  de combinación con la cobertura objetivo y las excepciones regionales. También
  debe registrar la cobertura realizada mediante intervalos o resúmenes
  apropiados.
- **REQ-MAN-49.** El manifiesto debe identificar las regiones configuradas o
  generadas para producir baja mapeabilidad, su relación con repeticiones o
  variantes estructurales y el mapeador y la referencia que determinaron los
  alineamientos y la MAPQ observados.
- **REQ-MAN-50.** En el modo directo, el manifiesto debe registrar las
  proporciones y distribuciones de MAPQ predeterminadas o sustituidas y las
  cantidades realizadas de lecturas correctamente alineadas, mal alineadas,
  multimapeadas y no alineadas.
- **REQ-MAN-51.** Para las lecturas quiméricas, el manifiesto debe registrar el
  perfil, la frecuencia predeterminada o sustituida, la cantidad realizada y la
  ubicación de la verdad conocida que describe sus fragmentos y puntos de unión.
- **REQ-MAN-52.** Para las anomalías *paired-end*, el manifiesto debe registrar
  las tasas configuradas y las cantidades realizadas por tipo y combinación, y
  relacionarlas con la verdad de fragmento y estado de origen.
- **REQ-MAN-53.** El manifiesto debe registrar cada grupo de lectura, su muestra
  o individuo, biblioteca, tanda o carril, plataforma o tecnología, parámetros
  efectivos y archivos asociados mediante los mismos identificadores usados en
  FASTQ, BAM o CRAM y ficheros de procedencia. Debe relacionar todos los grupos
  del individuo con su único BAM o CRAM final dentro de cada conjunto y
  referencia de mapeo. Debe enlazar cada FASTQ final separado con su grupo de
  lectura y con el BAM o CRAM final de destino.
- **REQ-MAN-54.** El manifiesto debe identificar la única tecnología de cada
  experimento, sus bibliotecas, tandas, carriles, perfil y mapeo compatibles, y
  enlazar los experimentos que usen las mismas muestras y verdad
  genotípica sin mezclar sus archivos. Debe registrar el identificador y la
  versión exactos del perfil, incluidos la química y el procesamiento cuando
  formen parte de la definición tecnológica.
- **REQ-MAN-55.** Para cada CRAM, el manifiesto debe registrar el identificador,
  la ubicación y las comprobaciones de integridad de la copia FASTA exacta usada
  para mapearlo. Las cabeceras y el manifiesto no deben permitir asociar el CRAM
  con otra referencia.
- **REQ-MAN-56.** Para cada BAM o CRAM final, el manifiesto debe registrar el
  formato, el orden efectivo, el tipo y la ubicación del índice y el resultado
  de la validación de compatibilidad entre alineamiento, índice y referencia.
- **REQ-MAN-57.** El manifiesto debe relacionar el único VCF multimuesta de cada
  conjunto con su referencia y con todos sus individuos, usando exactamente los
  mismos identificadores de muestra que las demás salidas.
- **REQ-MAN-58.** Para cada VCF de verdad, el manifiesto debe registrar el orden,
  si se entregó comprimido con BGZF o como texto, el tipo y la ubicación del
  índice cuando corresponda y el resultado de validar VCF, compresión, índice,
  cabecera y referencia.
- **REQ-MAN-59.** El manifiesto debe enlazar cada resultado estructurado de la
  validación interna con los FASTQ, BAM o CRAM, muestras, grupos de lectura y
  conjuntos comprobados, y resumir las desviaciones detectadas y sus
  explicaciones registradas. Para cada comprobación estadística debe conservar
  la tolerancia, el valor esperado, el valor realizado y la decisión aplicada.
  También debe identificar los conjuntos que no se marcaron como completos por
  un error estructural, lógico o estadístico elevado a error.
- **REQ-MAN-60.** Para cada biblioteca WGS, el manifiesto debe registrar la
  preparación sin PCR o con amplificación, el objetivo original de cobertura o
  número de lecturas, el indicador equivalente realizado y el perfil y las
  sustituciones aplicados.
- **REQ-MAN-61.** El manifiesto debe registrar la definición de cobertura usada,
  la cobertura total respecto a la longitud haploide y, cuando pueda medirse, la
  cobertura por cromosoma y copia o haplotipo antes y después de sesgos y
  degradaciones.
- **REQ-MAN-62.** Para la contaminación debe registrar la matriz objetivo y
  realizada o los pesos del reparto automático, los límites comprobados y, por
  fuente y mecanismo, lecturas propias retiradas, contaminantes incorporadas,
  transferidas, descartadas, totales finales y cobertura final.
- **REQ-MAN-63.** Cuando un experimento sin degradaciones actúe como control de
  otro, el manifiesto debe enlazar ambos mediante la verdad, referencias,
  muestras biológicas y relación experimental, manteniendo separados sus
  archivos. Cada evento debe conservar su propia subsemilla y generar sus
  propias lecturas.
- **REQ-MAN-64.** El manifiesto debe enlazar el historial resumido de cada
  generación y registrar N objetivo heredado o modificado, N realizado, su
  desviación y aviso, la transición en que el objetivo determinó la descendencia,
  la forma original y el objetivo de Ne o su ausencia, Ne realizado, desviación y advertencia cuando
  haya objetivo, apareamientos, progenitores con descendencia, descendientes,
  distribución familiar y estados genotípicos completos conservados.
- **REQ-MAN-65.** Debe registrar la tasa de autofecundación, sus valores
  esperados y realizados, el único locus seleccionado, sus alelos, genotipos o
  combinaciones de dosis, valores genéticos, valores del carácter, función de
  aptitud, aptitudes antes y después de normalizar,
  supervivientes, extinciones, el modo y los valores de migración solicitados y
  realizados por ruta y las rutas individuales. Debe registrar que el valor del
  carácter coincidió con el valor genético y que no se usó componente ambiental
  ni heredabilidad.
  También debe registrar el sentido direccional elegido y los parámetros
  efectivos de la función.
- **REQ-MAN-66.** Debe registrar las tasas y espectros de mutación, eventos
  esperados y realizados, mutaciones de novo por haplotipo, transmisión y orden
  efectivo de procesos de cada generación.
- **REQ-MAN-67.** Debe enlazar el catálogo FASTA de TE y conservar por copia o
  evento sus identificadores, fuente, familia, mecanismo, secuencias, orientación,
  truncamiento, sitio diana, ubicación, haplotipo y linaje.
- **REQ-MAN-68.** Debe registrar el modo y unidades de recombinación, los mapas
  completos o tasas efectivas y los recuentos agregados de intervalos
  bloqueados, bivalentes formados y entrecruzamientos realizados en los
  resúmenes generacionales. Para cada secuencia debe registrar su clasificación
  como cromosoma o *scaffold*, si se aplicó el mínimo de un entrecruzamiento y,
  en los *scaffolds*, la unidad y el umbral efectivos y su procedencia general o
  específica. Para los individuos muestreados debe enlazar la salida detallada
  de puntos, segmentos, haplotipos transmitidos, copia haplotípica fundadora y
  capas de mutación de novo.
- **REQ-MAN-69.** Debe identificar las entradas, comprobaciones, perfiles,
  herramientas, versiones, parámetros efectivos y dependencias relevantes para
  regenerar decisiones, y declarar por separado los niveles de reproducibilidad
  biológica, técnica aleatoria, semántica y binaria de REQ-CONF-111 a
  REQ-CONF-115. Para cada nivel debe indicar si está garantizado, no garantizado
  o no es aplicable, las condiciones comprobadas y la compatibilidad entre
  versiones cuando exista sin presuponerla.
- **REQ-MAN-70.** Cuando se use el modelo de equilibrio de Hardy-Weinberg o F
  para inicializar una población, debe registrar la ploidía, las frecuencias
  alélicas de entrada y realizadas y las probabilidades, frecuencias y recuentos
  genotípicos esperados y realizados por locus y población o subpoblación. Si se
  usa F, debe registrar su objetivo y el F realizado por subpoblación y
  generación registrada, incluidas las estimaciones negativas. Debe identificar
  el modo y, en el controlado, distinguir los recuentos alélicos fijados de la
  distribución genotípica aleatoria condicionada.
- **REQ-MAN-71.** Para cada población inicial, debe registrar si se usó
  asignación independiente o bloques con LD. Cuando se configure LD, debe
  registrar los bloques, su especificación objetivo, las frecuencias conjuntas
  realizadas, las frecuencias alélicas derivadas y la fase de los haplotipos
  generados. Si se aportaron frecuencias alélicas para validación, debe registrar
  sus discrepancias y la tolerancia aplicada. También debe identificar el modo
  probabilístico o controlado y, en el segundo, C, los recuentos haplotípicos y
  su distribución condicionada entre individuos.
- **REQ-MAN-72.** La raíz de una ejecución con eventos debe contener un
  manifiesto general que enlace la configuración resuelta, las referencias y
  recursos comunes, las subsemillas y el manifiesto propio de cada evento.
- **REQ-MAN-73.** Cada evento debe tener un directorio y un manifiesto propios.
  Este manifiesto debe registrar su identificador, selección de individuos,
  referencias, subsemilla y salidas solicitadas, generadas, omitidas o fallidas.
  Cuando el evento genere FASTQ, BAM o CRAM debe registrar su tecnología, perfil
  y bibliotecas; cuando genere alineamientos debe registrar además su mapeador
  único. Un evento que solo materialice verdad o VCF no debe inventar esos
  campos técnicos como si se hubieran aplicado.
- **REQ-MAN-74.** Dentro de un evento, el manifiesto debe separar los resultados
  por referencia de mapeo y relacionar cada BAM, CRAM o VCF exclusivamente con
  su conjunto. Debe conservar los identificadores biológicos de muestra, sin
  crear correspondencias entre las variantes de conjuntos diferentes.
- **REQ-MAN-75.** Las referencias y sus índices deben almacenarse una sola vez
  entre los recursos comunes de la ejecución. Los manifiestos de evento deben
  referenciarlos mediante su ruta y un identificador verificable, sin crear
  copias redundantes dentro de cada evento.
- **REQ-MAN-76.** El manifiesto general y cada manifiesto de evento deben incluir
  sumas de comprobación de los ficheros que enlazan y el algoritmo utilizado.
  El sistema debe comprobar esas sumas antes y después de exportar un evento.
  El algoritmo concreto se decidirá durante la planificación de la
  implementación.
- **REQ-MAN-77.** Los identificadores de grupo de lectura deben ser únicos en
  toda la ejecución y distinguir evento, biblioteca, tanda y carril. El
  manifiesto debe relacionarlos de forma coherente con FASTQ, BAM o CRAM y
  procedencia, mientras `SM` conserva el identificador biológico del individuo.
- **REQ-MAN-78.** Si el mismo individuo participa en varios eventos, el
  manifiesto general debe enlazar cada observación con el mismo identificador
  biológico y la misma verdad genotípica, fase y ascendencia, y con las
  subsemillas y salidas técnicas independientes de cada evento.
- **REQ-MAN-79.** Para variantes pequeñas explícitas solapadas, el historial
  estructurado debe registrar el orden, el marco de coordenadas anterior a cada
  evento, las relaciones con la referencia o los segmentos identificados y el
  estado final activo, modificado o eliminado. El manifiesto debe enlazar este
  historial con la verdad final y, cuando exista, con el VCF y su complemento.
- **REQ-MAN-80.** El manifiesto de cada muestreo debe enlazar por separado su
  resumen poblacional y su resumen muestral. Para cada uno debe indicar el
  ámbito, la generación, la población o subpoblación, los individuos o el número
  de individuos usados y la procedencia de cada frecuencia alélica,
  heterocigosidad, F, N y Ne. No debe agrupar ambos ámbitos bajo una etiqueta
  común que permita confundirlos.
- **REQ-MAN-81.** Para cada conjunto de frecuencias o probabilidades que deba
  sumar 1, el manifiesto debe registrar la suma y los valores originales, la
  tolerancia general efectiva, la decisión aplicada y, cuando se normalice, la
  suma y los valores finales. Debe identificar por separado cualquier tolerancia
  específica que se haya aplicado con otra finalidad.
- **REQ-MAN-82.** Para cada división, fundación, migración o selección de
  descendientes, el manifiesto debe registrar la regla del catálogo común, sus
  condiciones y combinación lógica, la cantidad exacta o proporción solicitada,
  los elegibles, la semilla cuando corresponda, los seleccionados y la cantidad
  o proporción realizada.
- **REQ-MAN-83.** El manifiesto general y cada manifiesto de evento deben
  registrar los estados incompleto, completo o fallido que correspondan a la
  ejecución, el evento, cada conjunto y cada artefacto. Para todo fallo deben
  enlazar el diagnóstico y distinguir los resultados publicados y validados de
  los temporales, incompletos o fallidos.
- **REQ-MAN-84.** El manifiesto debe indicar qué estados genotípicos completos se
  materializaron conforme a REQ-SAL-02, qué generaciones conservaron solo un
  historial resumido y qué configuraciones o muestreos solicitaron
  materializaciones adicionales. No debe describir un resumen como si fuera un
  estado genotípico completo.
- **REQ-MAN-85.** Para cada conversión conjunta de contaminación, el manifiesto
  debe registrar las unidades originales, los objetivos convertidos, la
  asignación entera, las desviaciones y los empates resueltos con la subsemilla.
  Para cada instancia contaminante o transferida debe conservar el identificador
  técnico único y su vínculo con el origen conforme a REQ-DAT-68 a REQ-DAT-70.
- **REQ-MAN-86.** Para cada secuencia FASTA, el manifiesto o la salida
  estructurada de referencias debe conservar el identificador formado por el
  primer elemento de la cabecera y la descripción completa restante, además del
  resultado de validar el alfabeto A, C, G, T y N.

El formato y el esquema del manifiesto se elegirán durante la planificación de
la implementación teniendo en cuenta que el núcleo se desarrollará en Rust.
Esta decisión no forma parte de la especificación funcional.

### 4.8. Resultado mínimo de una ejecución

- **REQ-SAL-01.** Toda ejecución válida debe producir un único FASTA por cada
  referencia inicial o derivada y los índices que requieran las operaciones que
  se hayan ejecutado.
- **REQ-SAL-02.** Toda ejecución válida debe producir un fichero de verdad con el
  estado genotípico completo, faseado y reconstruible exigido por cada vía. En
  creación directa debe materializar el estado completo de todos los individuos
  creados. En evolución debe materializar las poblaciones inicial y final y
  todos los individuos de las generaciones o muestreos solicitados, mientras el
  resto de generaciones conserva el historial resumido de REQ-REP-20. En
  familias y poblaciones segregantes debe materializar todos los individuos de
  todos los eventos y poblaciones intermedias conforme a REQ-CRU-15.
- **REQ-SAL-03.** Toda ejecución válida debe guardar en el directorio de
  resultados una copia de todos los ficheros de configuración realmente usados.
- **REQ-SAL-04.** Toda ejecución válida debe producir el manifiesto general
  completo y todos los manifiestos de evento que correspondan.
- **REQ-SAL-05.** Una ejecución solo debe terminar con éxito cuando haya
  producido todas las salidas mínimas y todas las salidas activadas, ya sea de
  forma predeterminada o por solicitud expresa, y todos sus eventos y conjuntos
  obligatorios estén completos conforme a REQ-SAL-21 a REQ-SAL-26.
- **REQ-SAL-06.** Si no se solicita VCF, la verdad genotípica debe seguir siendo
  accesible y reconstruible mediante el fichero de verdad obligatorio. Su
  formato y esquema técnicos todavía no están definidos.
- **REQ-SAL-07.** Cuando se use un FASTQ externo como fuente de contaminación,
  su copia en el directorio de resultados debe formar parte de las salidas
  activadas necesarias para que el experimento sea autocontenido y exportable.
- **REQ-SAL-08.** En el modo realista, una ejecución solo debe terminar con
  éxito si existe por cada grupo de lectura su FASTQ final —o pareja de FASTQ— y
  su fichero estructurado de procedencia enlazado desde el manifiesto.
- **REQ-SAL-09.** Cuando se simulen lecturas, la verdad conocida debe conservar
  las diferencias entre la secuencia original y la observada y la procedencia de
  los duplicados, aunque su organización técnica dentro de los ficheros de
  verdad todavía no esté definida.
- **REQ-SAL-10.** Cuando se estime un perfil empírico, su fichero estructurado
  reutilizable debe formar parte de las salidas activadas necesarias para que la
  ejecución termine con éxito.
- **REQ-SAL-11.** Todo BAM o CRAM final debe incluir su índice compatible entre
  las salidas activadas. Una ejecución no debe terminar con éxito si falla la
  validación conjunta del alineamiento, el índice y la referencia.
- **REQ-SAL-12.** Cuando se solicite VCF, la salida activada debe incluir el VCF
  multimuesta comprimido con BGZF y su índice compatible, salvo que se solicite
  expresamente la alternativa de texto sin comprimir. La ejecución no debe
  terminar con éxito si falla la validación correspondiente.
- **REQ-SAL-13.** Siempre que la ejecución genere FASTQ o BAM o CRAM, los
  resultados estructurados de su validación interna deben formar parte de las
  salidas activadas y estar enlazados desde el manifiesto.
- **REQ-SAL-14.** Cuando se utilice un catálogo de TE, su FASTA suministrado o
  generado debe formar parte de las salidas activadas y quedar enlazado desde el
  manifiesto. La copia suministrada debe conservar intacto el original.
- **REQ-SAL-15.** Toda ejecución evolutiva debe producir el historial
  estructurado resumido de todas sus generaciones y los estados genotípicos
  completos exigidos por REQ-REP-19 y REQ-REP-21.
- **REQ-SAL-16.** Cada muestreo evolutivo debe producir para sus individuos la
  salida detallada de recombinación y origen por segmentos exigida en
  REQ-REC-11 a REQ-REC-13. Esta salida no debe extenderse a todos los individuos
  históricos por el mero hecho de generar resúmenes de cada generación.
- **REQ-SAL-17.** Cada evento explícito debe producir su propio
  directorio y manifiesto, separado de los demás eventos y organizado por
  referencia de mapeo cuando haya varias.
- **REQ-SAL-18.** La raíz de la ejecución debe producir un manifiesto general y
  un área común de referencias e índices, enlazados desde los eventos mediante
  rutas e identificadores verificables.
- **REQ-SAL-19.** La exportación independiente de un evento debe crear un
  paquete autocontenido con su manifiesto, sus salidas, la configuración
  necesaria, los metadatos y comprobaciones requeridos y copias de las
  referencias e índices comunes que utilice. Debe validar las sumas de
  comprobación antes y después de crear el paquete.
- **REQ-SAL-20.** Cada evento de muestreo debe incluir entre sus salidas mínimas
  dos resúmenes estructurados y diferenciados: el de la población completa y
  sus subpoblaciones relevantes, y el de la muestra exportada. Ambos deben
  conservarse aunque no se generen salidas técnicas.
- **REQ-SAL-21.** La ejecución, cada evento, cada conjunto y cada artefacto
  generado deben tener un estado explícito. Como mínimo deben distinguirse los
  estados incompleto, completo y fallido. El manifiesto general debe conservar
  el estado de la ejecución y enlazar los estados y diagnósticos de sus unidades
  subordinadas.
- **REQ-SAL-22.** Todo artefacto debe construirse en un área temporal y
  publicarse de forma atómica en los resultados solo después de superar sus
  comprobaciones obligatorias. Un artefacto incompleto o fallido no debe
  presentarse como salida válida ni poder exportarse. Su diagnóstico sí debe
  conservarse cuando sea posible sin confundirlo con el artefacto solicitado.
- **REQ-SAL-23.** Un conjunto solo debe marcarse como completo cuando todos los
  artefactos activados para su referencia estén publicados y validados. El fallo
  de uno debe marcar el conjunto como fallido, sin publicar una versión
  parcial del artefacto afectado.
- **REQ-SAL-24.** Un evento solo debe marcarse como completo cuando sus salidas
  propias y todos sus conjuntos obligatorios estén completos. Un conjunto
  completo que no dependa del resultado fallido de otro debe conservar su estado
  y puede seguir siendo exportable, aunque el evento o la ejecución global
  terminen con estado fallido. La exportación debe declarar el fallo global y no
  presentar el paquete como una ejecución completa.
- **REQ-SAL-25.** La ejecución solo debe marcarse como completa y exitosa cuando
  todas las salidas mínimas y activadas, eventos y conjuntos obligatorios estén
  completos. Un fallo no debe borrar ni invalidar retroactivamente artefactos o
  conjuntos ya publicados y validados que sean independientes, pero el
  manifiesto general debe registrar la ejecución como fallida y la causa.
- **REQ-SAL-26.** Si la validación previa de la configuración o de una entrada
  impide iniciar la simulación, no debe publicarse ningún resultado. Las copias
  de trabajo creadas durante la validación deben permanecer en el área temporal
  o eliminarse de forma segura. El sistema debe emitir un diagnóstico que no se
  presente como manifiesto de una ejecución completa.

### 4.9. Escala y prioridades de uso

- **REQ-NOF-01.** El sistema no debe imponer límites funcionales fijos al número
  o tamaño de cromosomas o *scaffolds*, al tamaño total del genoma, al número de
  individuos o generaciones ni a la complejidad del escenario.
- **REQ-NOF-02.** Los límites prácticos de una ejecución deben depender de la
  memoria, CPU, almacenamiento y tiempo disponibles. El sistema no debe prometer
  capacidad o rendimiento ilimitados.
- **REQ-NOF-03.** La primera versión debe priorizar la creación y comprobación de
  ejemplos pequeños, controlados y reproducibles para casos concretos.
- **REQ-NOF-04.** La arquitectura y los formatos elegidos no deben impedir por
  diseño escenarios mayores o más complejos, aunque su optimización se aplace
  frente a los casos pequeños.
- **REQ-NOF-05.** Antes de ejecutar, el sistema debe estimar de forma aproximada
  el volumen de las salidas y el almacenamiento, la memoria y el tiempo
  necesarios, usando la configuración y los recursos que pueda detectar.
- **REQ-NOF-06.** Si el escenario parece grande para los recursos disponibles,
  el sistema debe mostrar advertencias claras. Una advertencia no debe bloquear
  por sí sola. Debe impedirse el inicio si la configuración es inválida, falta
  inevitablemente un requisito mínimo o una estimación supera un presupuesto de
  recursos configurado, explicando la causa.
- **REQ-NOF-07.** Las estimaciones no deben presentarse como exactas. Los modelos
  y umbrales de advertencia se elegirán y validarán durante la implementación.
- **REQ-NOF-08.** Para reducir el almacenamiento, el sistema no debe conservar
  una segunda copia limpia de los FASTQ ni duplicar las secuencias o calidades en
  el fichero de procedencia. Las estimaciones de espacio deben reflejar esta
  política.

## 5. Flujo funcional confirmado

El flujo funcional describe la relación entre los resultados, no el orden
técnico en que deberá implementarse el sistema.

1. El programa lee el fichero general de configuración y obtiene de él la vía
   de creación, cero, uno o varios eventos explícitos y,
   cuando se haya proporcionado, la semilla aleatoria. Si no se ha proporcionado una semilla,
   genera una antes de iniciar la simulación. Antes de simular, comprueba que
   cada cromosoma o *scaffold* tenga una ploidía efectiva entera de 1 a 10 y
   rechaza cualquier valor o ausencia no admitidos. No crea eventos implícitos:
   una ejecución sin eventos no genera VCF, FASTQ, BAM ni CRAM. Un evento que
   solo solicite verdad genotípica o VCF no necesita parámetros de
   secuenciación. Si un evento genera FASTQ, BAM o CRAM, también rechaza la
   mezcla de tecnologías o los perfiles, bibliotecas o parámetros de mapeo
   incompatibles. Comprueba entonces que el protocolo sea WGS de DNA genómico,
   que cada biblioteca elija exclusivamente cobertura o número de lecturas y
   que declare su preparación sin PCR o con amplificación. Rechaza paneles,
   amplicones, GBS y otros protocolos de captura o reducción de representación.
2. Si se ha seleccionado la vía evolutiva o la vía de familias y poblaciones
   segregantes, el programa obtiene del fichero general la referencia al
   fichero de configuración poblacional propio del escenario y lee este segundo
   fichero. Conserva una copia en resultados de todos los ficheros de
   configuración realmente usados. Con toda la configuración disponible, estima
   el volumen y almacenamiento de las salidas, la memoria y el tiempo necesarios
   según los recursos detectables. Advierte si el escenario parece grande, pero
   solo impide iniciarlo si la configuración es inválida o no puede satisfacerse
   inevitablemente un requisito mínimo, indicando la causa. Si se han indicado
   FASTQ externos como fuentes contaminantes, los valida y los copia al
   directorio de resultados antes de usarlos. La ejecución usa las copias y no
   modifica los originales. Si se han proporcionado un FASTQ real y su
   referencia para calibración, estima el perfil empírico, documenta sus
   limitaciones y advertencias y guarda un fichero estructurado reutilizable. El
   perfil solo sustituye al tecnológico predeterminado si la configuración lo
   selecciona expresamente; también puede seleccionarse un perfil empírico
   generado en una ejecución anterior.
   Si se aporta un catálogo FASTA de elementos transponibles, lo valida, lo
   copia a resultados y usa la copia sin modificar el original. Si el catálogo
   debe simularse, genera y conserva su FASTA enlazado.
3. Si la referencia inicial debe simularse, el programa genera sus cromosomas o
   *scaffolds* con los tamaños y la composición configurados, y obtiene el tamaño
   total a partir de sus tamaños individuales. Cuando se hayan configurado
   repeticiones en tándem concretas o aleatorias, incorpora a la secuencia sus
   motivos, copias e interrupciones realizados. Si la referencia es suministrada,
   el programa valida el FASTA real indicado antes de iniciar la simulación. Solo
   lo utiliza como referencia inicial si es legible, contiene al menos una
   secuencia no vacía, el primer elemento de cada cabecera es un identificador
   no vacío y único, y sus secuencias contienen únicamente A, C, G, T y N, en
   mayúsculas o minúsculas. El resto de cada cabecera se conserva como
   descripción. Si la validación falla, el programa
   detiene la ejecución sin generar resultados y explica claramente los errores
   detectados. Si la validación es correcta, copia el FASTA al directorio de
   resultados y utiliza esa copia, sin modificar el archivo original.
4. El sistema conserva si la referencia inicial fue simulada o suministrada y,
   cuando proceda, la procedencia original del FASTA y el nombre o la ruta de la
   copia utilizada.
5. Toda la secuencia de la referencia inicial se guarda en un único FASTA, que
   puede contener uno o varios cromosomas o *scaffolds*. No se genera ningún
   índice si las operaciones y los formatos posteriores no lo necesitan; en
   caso contrario, se generan todos los tipos necesarios.
6. A partir de una referencia ya disponible, se crea una nueva referencia
   estructural aplicando solo las duplicaciones, deleciones, inversiones y
   translocaciones solicitadas.
7. Toda la secuencia de la referencia derivada se guarda en un único FASTA. Se
   conservan su referencia de origen y los cambios aplicados, y se generan los
   tipos de índice que necesiten las operaciones o los formatos usados. Los
   pasos 6 y 7 pueden repetirse para disponer de dos o más referencias
   estructurales.
8. Los ficheros de configuración indican si los cambios se definirán de forma
   concreta o se crearán al azar. Para la generación aleatoria conservan los
   valores globales y los valores específicos de cada individuo. Un valor
   individual sustituye al global solo para ese individuo. Toda cantidad
   configurada es una cantidad media esperada; cuando se usa una densidad, la
   media se deriva de la densidad y de la longitud válida del ámbito. El sistema
   obtiene al azar un número realizado, que puede diferir de la media, usando la
   distribución de Poisson de forma predeterminada u otra distribución y sus
   parámetros seleccionados para la clase de variante y la región. Para cada
   clase se puede aplicar una distribución espacial uniforme o configurar
   regiones con una cantidad media esperada o densidad mayor o menor que el
   valor general, incluidas regiones con densidad cero. Cada región identifica
   una referencia estructural, un cromosoma o *scaffold* y un intervalo, y se
   interpreta sobre la referencia estructural correspondiente al individuo o la
   población. Los valores espaciales siguen la jerarquía genoma, cromosoma o
   *scaffold* y región: el ámbito más específico sustituye al más general, junto
   con la sustitución por individuo cuando exista. Antes de simular, el sistema
   rechaza los valores incompatibles asignados por regiones del mismo nivel que
   se solapen e identifica ambas regiones y el parámetro en conflicto. En las
   poblaciones, la frecuencia alélica se expresa como una proporción
   poblacional; para un individuo se usa la presencia o ausencia, el genotipo o
   la dosis alélica según su ploidía. Los valores fijados para individuos son
   restricciones obligatorias. El sistema asigna los alelos de los demás
   individuos para aproximarse al objetivo poblacional dentro de la tolerancia
   admitida. Si no existe una asignación válida, rechaza la configuración antes
   de simular, identifica la incompatibilidad y no modifica los valores
   individuales. Para inserciones y deleciones aplica el umbral configurado o,
   si se omitió, 50 pb, tanto al generar como al clasificar y calcular
   estadísticas. Las repeticiones en tándem pueden definirse en loci concretos o
   generarse al azar y por regiones; su motivo, copias o longitud, pureza,
   interrupciones y variabilidad se aplican por haplotipo, individuo o población
   con las sustituciones configuradas. Los cambios de copias usan el modelo
   escalonado de ganancia o pérdida de unidades completas de forma predeterminada
   u otro modelo seleccionado. Las mutaciones internas que cambian
   interrupciones, pureza o motivo solo se aplican si se activan expresamente.
9. El sistema aplica la vía seleccionada:

   - crea directamente individuos con los cambios configurados, incluidas
     inserciones de TE en un haplotipo nombrado y en una coordenada exacta o un
     sitio permitido elegido reproduciblemente;
   - crea una población inicial desde una referencia y simula su evolución
     generación a generación mediante las fases ordenadas del fichero de
     configuración poblacional; la salida de cada fase alimenta a la siguiente.
     Cada locus o bloque y subpoblación inicial elige entre el modo probabilístico
     predeterminado y el modo controlado explícito, sin mezclarlos dentro de la
     misma unidad. Puede usar
     equilibrio de Hardy-Weinberg para cualquier ploidía, frecuencias explícitas
     de genotipo o dosis fuera del equilibrio, o genotipos multilocus
     completamente faseados. Por defecto asigna loci independientemente. El LD
     se define mediante bloques no solapados de haplotipos y frecuencias
     conjuntas, o recuentos controlados que suman las copias disponibles; fuera
     de ellos los loci son independientes. El coeficiente de
     endogamia F puede definirse globalmente o sustituirse por subpoblación, pero
     no por locus, región o bloque. Para ploidías mayores que 1, F progresa desde
     las probabilidades de Hardy-Weinberg en F=0 hasta la homocigosis completa en
     F=1, sin modelos específicos de herencia polisómica o disómica. El F
     realizado se obtiene por subpoblación y generación registrada comparando
     heterocigosidad observada y esperada; puede ser negativo y se conserva sin
     error. Todos los individuos iniciales son fundadores con identificadores
     estables de individuo y copia haplotípica.
     N objetivo es obligatorio al iniciar cada población o subpoblación y se hereda hasta que un
     valor o evento demográfico lo cambia. Ne objetivo puede ser un valor, una
     proporción Ne/N o quedar sin configurar; siempre se calcula Ne realizado
     como Ne reproductivo a partir de los progenitores que dejaron descendencia
     y de la desigualdad de sus contribuciones. Un cambio de N fija el número de
     descendientes de la generación siguiente: un cuello de botella lo reduce y
     una expansión lo aumenta, sin eliminar retrospectivamente progenitores ya
     creados. Cada generación aplica, por defecto, eventos demográficos,
     migración dirigida, selección de viabilidad, elección de parentales,
     reproducción, segregación y recombinación, mutación de novo en la
     descendencia, cálculo de N realizado, Ne realizado y frecuencias,
     comparación con N objetivo y Ne objetivo, y muestreos o salidas. Cualquier desviación de Ne se
     cuantifica y advierte, pero no invalida ni detiene la simulación. La deriva
     surge del muestreo finito y no se aplica como un paso adicional.
     Una división reparte sin reemplazo todos los progenitores del origen entre
     hijas viables. División, fundación, migración y selección de descendientes
     usan el mismo catálogo declarativo de identificadores, azar, pertenencia y
     condiciones genéticas combinables, con cantidad exacta o proporcional. No
     admiten scripts, código ni funciones personalizadas. Cada hija tiene su
     propio N objetivo, independiente de los
     progenitores recibidos, y la suma puede representar expansión o
     contracción. Una fusión reúne todos los orígenes en una población nueva
     con parámetros propios y un N objetivo independiente del conjunto
     reproductor reunido. Su contribución predeterminada por origen es
     probabilística y proporcional a los progenitores elegibles, con pesos
     configurables. El modo controlado alternativo garantiza recuentos
     parentales exactos, no porcentajes exactos de bases ancestrales. Se
     registran por separado las contribuciones parentales y la ascendencia
     genómica realizada con trazabilidad haplotípica. Una fundación
     transfiere sus fundadores a una
     población nueva con N objetivo independiente y deja que el origen continúe
     con los restantes. Ambas pueden expandirse o contraerse en su primera
     descendencia. Una
     extinción no genera descendencia y termina la subpoblación; una
     recolonización posterior usa otro identificador y otra historia. Todos los
     eventos conservan la procedencia individual y haplotípica aplicable. Si
     coinciden varios eventos demográficos, la configuración aporta un orden
     explícito obligatorio y cada evento opera sobre el estado dejado por el
     anterior; las dependencias inválidas se rechazan antes de simular.
     Cada fase usa exclusivamente proporciones o recuentos exactos de migración
     por pares dirigidos. La selección usa el catálogo declarativo común. Las
     rutas salientes se resuelven conjuntamente y sin
     reemplazo antes de selección y reproducción. En modo proporcional, las
     rutas son categorías excluyentes y la permanencia tiene probabilidad uno
     menos la suma de tasas; el recuento realizado es aleatorio y reproducible.
     La selección usa el mismo locus estable durante toda la ejecución. Asigna
     valores genéticos a sus genotipos o dosis y conserva por separado el valor
     del carácter, aunque en esta versión ambos valores coinciden. No añade ruido
     ambiental ni usa heredabilidad. Una función configurable obtiene la aptitud
     relativa y favorece, según se indique, valores mayores o valores menores.
     No admite selección estabilizadora ni disruptiva. Los valores y la función
     pueden variar por fase o subpoblación. La reproducción
     aplica de forma probabilística la autofecundación configurada y, en los
     demás eventos, elige uniformemente dos supervivientes distintos. Genera el
     número de descendientes que permiten los progenitores e intenta alcanzar N
     objetivo. Registra N realizado, su desviación, la distribución de
     descendencia y Ne realizado. No usa para estimar este último la pérdida de
     heterocigosidad, el aumento de endogamia ni el LD. Una desviación de N
     advierte sin detener. Si N realizado es 0, declara una extinción por
     resultado solo cuando la política efectiva lo permite; en otro caso produce
     un error que detiene inmediatamente toda la simulación, sin omisiones ni
     continuación de otras ramas. Una extinción autorizada omite de forma
     atómica y transitiva los eventos dependientes, registra la cadena de causas
     y deja continuar las ramas independientes. Los
     individuos de la generación parental dejan de ser población activa cuando
     los sustituye la descendencia con el N realizado de cada subpoblación.
     La recombinación usa el único modo global elegido. En poliploides pares
     forma al azar bivalentes entre copias estructuralmente compatibles, bloquea
     intervalos no colineales y conserva los haplotipos transmitidos. El número
     de entrecruzamientos sigue Poisson según la longitud genética; cada
     cromosoma recibe al menos uno y cada *scaffold* lo recibe solo si supera su
     umbral efectivo en cM o bases. El umbral es general y admite sustituciones
     por *scaffold* sin mezclar unidades.
     El sistema mantiene internamente el linaje de segmentos hasta la copia
     haplotípica fundadora. Solo los individuos muestreados exportan puntos de
     recombinación, segmentos, coordenadas actuales y fundadoras, orientación,
     copias duplicadas y capas separadas de mutaciones de novo; las demás
     generaciones producen resúmenes.
     Después
     aplica las mutaciones de novo por base y generación a copias descendientes,
     incluidas las reglas separadas para SNP, inserciones, deleciones, variantes
     estructurales y repeticiones en tándem. No moviliza TE espontáneamente; los
     TE ya presentes se heredan como cualquier otra variante.
     Antes de simular rechaza un orden con dependencias incompatibles, una
     ploidía reproductiva impar o un mapa incompleto o inválido. Durante el
     proceso selecciona de forma explícita, aleatoria o combinada y conserva el
     estado de los individuos de los muestreos solicitados, junto con su
     generación y subpoblación o población de origen. La selección aleatoria es
     sin reemplazo dentro del evento y puede estratificarse mediante grupos no
     solapados y filtros de genotipo.
     Los muestreos no eliminan ni modifican individuos y su aleatoriedad se
     mantiene separada para que no cambien la evolución. Por defecto se conserva
     su estado genotípico y los metadatos mínimos de la observación, sin generar
     BAM, CRAM ni VCF. Cada evento conserva por separado un resumen de la
     población completa y las subpoblaciones relevantes y otro de la muestra
     exportada. Ambos identifican inequívocamente su ámbito e incluyen
     frecuencias alélicas, heterocigosidad, F, N y Ne cuando correspondan. Si un
     evento solicita archivos técnicos, constituye un
     experimento técnico independiente con tecnología, perfil, mapeador,
     subsemilla, directorio y manifiesto propios. Sus individuos siguen
     las referencias de mapeo y los conjuntos independientes de los pasos 11 a
     16. Todas las generaciones producen un historial resumido; solo se conserva
     el estado completo de la inicial, la final y las solicitadas; o
   - crea una familia o población segregante desde uno o más parentales según
     los cruces predefinidos o descritos por el usuario. Cada evento admite un
     parental en autofecundación o dos en cruzamiento sexual; los diseños con
     tres o más fundadores se expresan como una secuencia de eventos. El orden de
     los dos parentales no cambia el resultado biológico. Cada parental referencia
     por su identificador estable a un individuo disponible o se define
     expresamente en el fichero poblacional con las reglas comunes; antes de
     simular se comprueban su existencia, unicidad y compatibilidad. El sistema
     valida que el grafo dirigido de eventos no tenga ciclos ni referencias u
     órdenes temporales incorrectos, ejecuta cada evento según sus dependencias y
     número de descendientes y permite que esa descendencia alimente eventos
     posteriores. Los parentales posteriores se eligen mediante el catálogo
     declarativo común, que se valida antes de aplicarse. El sistema conserva
     el estado genotípico y los metadatos de todos los individuos intermedios,
     seleccionados o no. Por defecto no genera para ellos BAM, CRAM ni VCF, pero
     una población intermedia puede solicitarlos y entonces sigue los conjuntos
     y la validación de los pasos 11 a 16. Las plantillas iniciales son F1, F2,
     retrocruce, autofecundación y pedigrí simple. F2 permite autofecundar F1 o
     cruzar individuos F1; RIL, introgresión y MAGIC se expresan mediante cruces
     personalizados. Simula la recombinación con la tasa uniforme o los mapas
     completos generales o específicos de grupo y conserva los
     haplotipos faseados y el origen de los segmentos en toda la descendencia.
     Cuando el escenario admite mutaciones durante la transmisión, aplica el
     modelo de repeticiones en tándem y solo incluye mutaciones internas si están
     activadas.

10. Todos los individuos reciben la misma ploidía efectiva de 1 a 10 copias para
   cada cromosoma o *scaffold*, obtenida del valor global y sus sustituciones.
   Para cada individuo se conserva la referencia que corresponde a su estructura
   y un identificador único y estable, aportado por la persona usuaria o generado
   por el sistema. El identificador no se reutiliza y acompaña al individuo en sus
   relaciones parentales, observaciones y muestras. El sistema guarda el estado
   completo, faseado, con su dosis y reconstruible en el fichero de verdad
   obligatorio para las entidades exigidas por la vía conforme a REQ-SAL-02,
   aunque no se solicite VCF.
11. Para cada evento explícito que solicite alineamientos se generan en CRAM o,
    si se selecciona, BAM. El fichero general selecciona una o varias referencias
    de mapeo entre las disponibles y el modo de obtención. Cada biblioteca WGS
    elige cobertura media o número de lecturas, nunca ambos, y preparación sin
    PCR de forma predeterminada o WGS amplificada de forma opcional. En el modo
    realista se simulan lecturas cortas
    tipo Illumina, *single-end* o *paired-end*, o lecturas largas PacBio HiFi,
    desde los haplotipos faseados con el perfil tecnológico o empírico WGS
    efectivo, los parámetros y la semilla.
    Las posiciones iniciales se eligen uniformemente entre las copias accesibles.
    La cobertura se calcula sobre la longitud haploide y se reparte en proporción
    al número y longitud de las copias, de modo que cada copia recibe una
    cobertura semejante antes de excepciones. La cobertura local es estocástica
    y aproximadamente Poisson; no se fuerza una profundidad por base.
    Los alelos contribuyen según dosis y ploidía, sin sesgo alélico o *allele
    dropout* forzados. Después se aplican la curva de sesgo GC, las regiones con
    cobertura especial, los límites de secuencia y las demás degradaciones con
    su precedencia documentada. Se aplican las tasas y tipos de error
    efectivos, homogéneos salvo las sustituciones configuradas e incluidos los
    patrones por posición o contexto admitidos, y se generan puntuaciones de
    calidad probabilísticamente coherentes. La verdad
    conserva las diferencias entre las secuencias originales y observadas. En la
    preparación WGS sin PCR no se generan duplicados de PCR. Si se elige un
    perfil amplificado, se generan sus duplicados y sesgos configurados. Los
    duplicados ópticos se tratan por separado; en
    *paired-end* se duplica la pareja completa y siempre se conserva la relación
    con la lectura o pareja de origen. Las parejas normales siguen el perfil y
    las proporciones configuradas generan orientación anómala, inserto inesperado
    o un miembro no alineado, conservando fragmento y estado de origen. También
    se generan con la frecuencia efectiva las lecturas quiméricas y se conservan
    sus fragmentos y puntos de unión reales. Si se solicitaron restos de
    adaptadores, se añaden con los parámetros efectivos y se conservan las
    lecturas y los intervalos afectados.
    Antes de elegir lecturas contaminantes, todas las contribuciones expresadas
    como proporción, cobertura equivalente o número se convierten conjuntamente
    a un reparto entero de lecturas o parejas completas. La asignación minimiza
    la desviación global y resuelve reproduciblemente los empates con la
    subsemilla del evento.
    La contaminación biológica sustituye de forma reproducible lecturas propias
    de la receptora por lecturas representativas de fuentes compatibles; no
    reduce los datos de la donante. Puede seguir una matriz dirigida o repartir
    una tasa total entre donantes de la tanda ponderados por su número de
    lecturas. La suma no supera el 100 %. En *paired-end* siempre actúa sobre
    parejas completas. Por separado, el demultiplexado transfiere cada lectura o
    pareja reasignada de la donante a la receptora sin duplicarla. Las lecturas
    descartadas o sin asignar no se reemplazan y reducen la cobertura final.
    Los FASTQ externos compatibles solo se admiten como fuentes en este modo
    realista. Cada copia biológica recibe una identidad técnica nueva y un
    vínculo estable con su lectura de origen; una transferencia técnica conserva
    la identidad de la instancia y la retira de la biblioteca donante. Después
    se crea
    un FASTQ final, o una pareja para datos *paired-end*, por individuo,
    biblioteca y tanda o carril; no se concatenan ni se conserva otro FASTQ
    limpio. Para cada FASTQ o pareja se genera un fichero de procedencia que
    enumera por identificador las lecturas o parejas contaminantes presentes y
    registra solo recuentos
    agregados de las descartadas o sin asignar, sin duplicar secuencias ni
    calidades. No se simulan señales, BCL ni lecturas i5 o i7 detalladas. A
    continuación, un único mapeador procesa el FASTQ final de cada muestra contra
    cada referencia seleccionada. Las lecturas de regiones repetitivas,
    duplicadas o altamente homólogas parten de sus posiciones y haplotipos
    verdaderos; el mapeador determina sus alineamientos y MAPQ sin que el sistema
    los sobrescriba. Por defecto es minimap2; su *preset* y valores
    predeterminados se derivan del perfil de lecturas y las sustituciones
    configuradas prevalecen tras advertirse y validarse cuando sea posible. Las
    calidades del FASTQ se propagan al BAM o CRAM y los duplicados se marcan o no
    según la configuración. Las lecturas se asignan a grupos coherentes de
    muestra, biblioteca, tanda o carril y tecnología. En el modo directo se
    generan BAM o CRAM sintéticos con las mismas reglas funcionales de cobertura,
    sustitución biológica, matriz o reparto automático, transferencia por
    índices y descarte, claramente etiquetados y sin crear ni mapear FASTQ. Solo
    pueden usar fuentes contaminantes simuladas con verdad conocida. Cada lectura
    se asigna a la categoría configurada de
    alineamiento correcto, incorrecto, múltiple o ausente, con MAPQ, campos,
    banderas, CIGAR y alineamientos secundarios o suplementarios coherentes. Las
    lecturas quiméricas y las anomalías *paired-end* mantienen esa coherencia con
    sus fragmentos y estados de origen. La referencia de mapeo puede diferir de
    la estructura de una muestra.
12. Se pueden introducir las demás degradaciones confirmadas en los datos de
    lecturas o alineamientos según el modo utilizado.
13. Si se solicita VCF, para cada referencia de mapeo seleccionada se genera un
    único VCF de verdad multimuesta, independiente de los demás conjuntos, con
    todos sus individuos como columnas y los genotipos simulados de SNP,
    *indels*, variantes estructurales y repeticiones en tándem, expresados
    respecto a la referencia del conjunto. El VCF se construye directamente
    desde el estado genotípico faseado de la simulación y nunca mediante llamada,
    imputación, faseado o corrección a partir de FASTQ, BAM o CRAM. Por ello,
    conserva el genotipo completo y la fase real de cada muestra aunque su
    cobertura sea nula o baja o sus lecturas presenten contaminación, errores,
    calidad baja, ausencia de alineamiento o mapeo incorrecto, múltiple o
    ambiguo. Ninguno de esos casos convierte un genotipo verdadero en ausente.
    `GT` conserva la fase y ploidía reales. No se añaden campos `FORMAT`
    derivados de evidencia o inferencia, como `DP`, `AD`, `GQ` o `PL`; la
    cobertura, el soporte alélico y las calidades se registran en los resultados
    de validación de FASTQ, BAM o CRAM. Todo registro válido usa `QUAL=.` y
    `FILTER=PASS`, y la cabecera documenta que no se estima una calidad
    probabilística. Un registro que no pueda representarse mediante el VCF y su
    complemento cuando corresponda, o que no sea validable, detiene la creación
    del conjunto y no se emite como `PASS`.
    Solo incluye loci donde al menos una muestra tiene un alelo no de referencia,
    aunque ese alelo esté fijado en todas las muestras;
    no añade posiciones invariantes ni
    bloques gVCF de forma predeterminada. El fichero estructurado obligatorio
    conserva la verdad completa. El VCF se refiere al mismo genoma que los BAM o
    CRAM correspondientes cuando existan. Cada VCF se construye directamente
    contra la secuencia y las coordenadas de su referencia, con identificadores
    de variante locales y sin correspondencias registro a registro con los VCF
    de otras referencias. El VCF expresa genotipos faseados cuando el formato permite
    representar el caso. Si el propio registro no admite la fase completa, el
    complemento la conserva sin pérdida y la limitación queda identificada; no
    se convierte silenciosamente en una pérdida de verdad. La clasificación de
    inserciones y deleciones usa el mismo umbral efectivo que la generación. Los
    SNP e *indels* usan REF y ALT explícitos y normalizados. Los ALT
    representables de un mismo locus se reúnen en un único registro multialélico
    que conserva ploidía, dosis y fase. Las variantes estructurales usan alelos
    simbólicos y SVTYPE, END y SVLEN cuando corresponda. Las variantes
    complejas y las repeticiones en tándem conservan información suficiente para
    reconstruir el cambio y el alelo real simulados. Si un caso no cabe sin
    pérdida en VCF, se genera un fichero complementario de verdad; el registro
    VCF lo referencia mediante un identificador estable compartido. Si todos los
    casos son representables, no se genera ese complemento. El VCF se ordena
    según la referencia y por coordenada, se comprime con BGZF y se indexa con
    TBI o CSI de forma predeterminada; el texto sin comprimir solo se entrega si
    se solicita. Antes de completarlo se validan el contenido, la compresión y el
    índice aplicables, la cabecera y la referencia.
14. Cada BAM, CRAM y VCF incluye en su cabecera o metadatos estándar la
    información necesaria para identificar la referencia, la muestra, el evento,
    la biblioteca, la tanda o carril y el grupo de lectura que correspondan.
    `SM` y las columnas de muestra conservan el identificador biológico estable y
    cada identificador `RG` es único en toda la ejecución.
15. Cuando se solicitan estas salidas, los pasos 11 a 14 producen un conjunto
    independiente por cada referencia de mapeo seleccionada. Todos los conjuntos
    contienen las mismas muestras, usan el único mapeador del evento cuando
    corresponde y no mezclan datos de referencias distintas. Los FASTQ del modo
    realista pertenecen exclusivamente al evento y se usan contra todas sus
    referencias seleccionadas; no se reutilizan entre eventos. Dentro de cada
    conjunto se entrega un único BAM o CRAM final por
    individuo, que combina todas sus bibliotecas, tandas y carriles y preserva
    sus grupos de lectura. Solo se conserva este archivo final con las
    degradaciones configuradas. Un control limpio es otro evento enlazado que
    usa la misma verdad, referencias y muestras, pero tiene su propia subsemilla,
    genera lecturas nuevas y mantiene separados sus resultados. Los resultados
    del evento contienen una única tecnología experimental; las mismas muestras
    y verdad pueden alimentar otros eventos técnicos enlazados, que generan
    nuevas lecturas y mantienen separados sus archivos.
    Cada alineamiento final se ordena por coordenadas por defecto, se acompaña de
    un índice compatible y se valida junto con el índice y la referencia. Cuando
    el formato es CRAM, el conjunto enlaza el FASTA exacto almacenado entre los
    recursos comunes. El VCF multimuesta usa los mismos identificadores de individuo y no
    incluye muestras de otros conjuntos.
16. Cuando se generan FASTQ o BAM o CRAM, el sistema los valida contra la verdad
    de origen de las lecturas y los artefactos configurados, sin ejecutar un
    *SNP caller* ni inferir genotipos desde los alineamientos. Comprueba la
    procedencia desde los haplotipos, el objetivo de cada biblioteca, la
    cobertura total y por copia, y el soporte alélico esperado según dosis y
    ploidía sin exigir balance exacto. Comprueba también las tasas realizadas de
    errores, contaminación,
    duplicación y demás degradaciones, la coherencia de las parejas y la
    compatibilidad de alineamientos, referencias y grupos de lectura. Las
    desviaciones intencionadas deben quedar explicadas por la configuración o
    por una decisión aleatoria registrada. El resultado estructurado conserva
    objetivos, valores observados, tolerancias, desviaciones, explicaciones y la
    decisión aplicada. Una incoherencia estructural o lógica detiene la ejecución
    y el conjunto afectado no se marca como completo. Una desviación estadística
    fuera de tolerancia genera una advertencia de forma predeterminada, salvo que
    la configuración la eleve a error y detenga la ejecución. La desviación
    entre Ne objetivo y Ne realizado es una excepción: siempre se cuantifica y
    advierte, pero nunca detiene la simulación. En simulaciones poblacionales
    valida además N objetivo y realizado, Ne objetivo y realizado, reproducción, selección,
    migración, eventos de estructura poblacional, mutación, orden de procesos y
    recombinación. Para cada muestreo comprueba por separado el resumen
    poblacional y el muestral y rechaza cualquier confusión de ámbito.
17. El sistema genera el manifiesto general y los manifiestos de evento con la
    vía utilizada, el fichero general y
    cada fichero de configuración poblacional utilizado, los parámetros
    realmente aplicados, el origen y los cambios de las referencias, y los
    nombres o las rutas de los archivos generados. También registra los valores
    globales de generación aleatoria, si cada clase usó una cantidad media
    esperada o una densidad, la unidad de densidad cuando proceda, la cantidad
    media esperada y el número realizado por clase, ámbito o región e individuo
    cuando corresponda, y la distribución probabilística y sus parámetros
    realmente usados. Registra además la distribución espacial uniforme o las
    regiones configuradas, sus referencias estructurales, cromosomas o
    *scaffolds* e intervalos, los valores de cada nivel espacial, las
    sustituciones individuales y los valores efectivos después de aplicar ambas
    precedencias. Conserva también la semilla usada, si fue
    suministrada o generada automáticamente, la versión del programa, las
    referencias de mapeo seleccionadas y el conjunto generado para cada una.
    Para cada objetivo de frecuencia alélica poblacional, registra la frecuencia
    objetivo, las restricciones individuales y la frecuencia realizada.
    En la vía evolutiva, registra además las fases en su orden, su duración y sus
    parámetros efectivos, y cada muestreo con su generación, tamaño o proporción
    solicitados, individuos seleccionados y semilla o flujo aleatorio aplicable.
    Para cada muestreo indica el estado genotípico y los metadatos conservados,
    enlaza por separado los resúmenes poblacional y muestral e indica las salidas
    adicionales solicitadas y los archivos que se generaron.
    También relaciona cada identificador estable de individuo con su generación
    de origen, sus parentales cuando existan, sus observaciones y sus archivos,
    distinguiéndolo de los identificadores de evento y conjunto.
    Registra las trayectorias de N objetivo y realizado y de Ne objetivo y
    realizado, sus desviaciones y advertencias, la estructura y los parámetros de las
    subpoblaciones, la pertenencia individual, las migraciones configuradas y
    realizadas por par y los eventos de división, fusión, fundación, extinción o
    recolonización con sus subpoblaciones, individuos, procedencias y resultados.
    Conserva además qué procesos evolutivos estuvieron activos y sus parámetros
    y valores efectivos por fase, generación y subpoblación, así como su orden efectivo
    y si procedía de la configuración o de la secuencia predeterminada.
    Registra el modo, las tasas, los mapas y los parámetros efectivos de
    recombinación y los eventos realizados necesarios para seguir el origen de
    los genotipos descendientes, incluidos el grupo y el perfil efectivo usados
    para cada individuo o evento.
    En familias y poblaciones segregantes, distingue los parentales referenciados
    de los definidos en la configuración y conserva su identificador,
    referencia, genotipo faseado por copia y linaje. Para los parentales expresos
    conserva los cambios asignados y el estado de referencia no modificado.
    Registra también los eventos del grafo en el orden realizado y los
    identificadores estables de sus fuentes, parentales y descendientes. Si un
    individuo alimenta varias ramas, mantiene una sola identidad y distingue los
    gametos y descendientes nuevos de cada evento.
    Para cada selección registra la población candidata, los identificadores o
    condiciones y su combinación lógica, la semilla o flujo aplicable y los
    individuos elegibles y elegidos. Distingue el número exacto de la proporción
    probabilística y registra las proporciones solicitada y realizada. Para cada
    población intermedia identifica los estados y metadatos conservados y las
    salidas opcionales solicitadas y generadas.
    Conserva la fase, los haplotipos y el origen de los segmentos y registra toda
    limitación de representación en VCF y el complemento que preserva la verdad.
    Registra también el umbral efectivo usado para clasificar inserciones y
    deleciones. Para las
    repeticiones en tándem registra los objetivos, las distribuciones y los valores
    realizados de motivo, copias o longitud, pureza, interrupciones y
    variabilidad, junto con cada evento de mutación, el modelo y los parámetros
    efectivos y el haplotipo afectado.
    Registra las estimaciones de volumen, almacenamiento, memoria y tiempo, las
    advertencias emitidas y el consumo real cuando pueda medirse.
    Para las lecturas y BAM o CRAM conserva el modo realista o directo, la
    tecnología, el perfil, los parámetros efectivos y la semilla o flujo
    aleatorio aplicable. Para el perfil de errores conserva su nombre, versión,
    tasas, tipos y sustituciones efectivas. Si se estimó un perfil empírico,
    registra sus entradas, parámetros, limitaciones y medidas de calidad, y
    enlaza el artefacto estructurado reutilizable. Para los duplicados conserva
    por tipo las tasas, cantidades esperadas y realizadas, lecturas o parejas de
    origen y si se entregaron marcados; cuando se marcaron, identifica el método
    o la herramienta. Para los restos de adaptadores registra el estado, los
    parámetros y las cantidades realizadas. Para el sesgo GC conserva la curva,
    intensidad, combinación con la cobertura y resúmenes de cobertura realizada.
    Identifica las regiones de baja mapeabilidad y la referencia y el mapeador
    que produjeron sus alineamientos y MAPQ. Relaciona cada FASTQ con su
    procedencia y con los conjuntos que lo utilizaron. Para cada mapeo registra
    el único mapeador, su
    versión, el comando efectivo completo, la referencia y los FASTQ. Si es
    minimap2, distingue el *preset*, los valores predeterminados y las
    sustituciones. También conserva los objetivos y valores realizados de
    cobertura y, para cada contaminación, las fuentes, la matriz o reparto
    automático, las fracciones objetivo y realizadas y los recuentos de lecturas
    propias retiradas y contaminantes incorporadas. En datos *paired-end*
    registra parejas completas. Registra por separado la contaminación técnica
    con su tanda, bibliotecas, muestras, índices, tasas, procedencias, errores y
    cantidades configuradas y realizadas de lecturas transferidas, sin asignar
    o filtradas, incluidos los totales y la cobertura finales.
    Enlaza cada FASTQ final con su fichero de procedencia, indica que no se
    conservó una copia limpia y resume las lecturas o parejas contaminantes y los
    recuentos agregados de descartes o lecturas sin asignar por muestra y causa u
    origen conocido.
    Para cada FASTQ contaminante externo conserva la ubicación original, la
    copia utilizada, su validación y los datos de integridad definidos.
    Cuando exista un complemento de verdad, enlaza sus casos con el VCF, la
    referencia y el conjunto mediante identificadores locales a ese conjunto.
    Enlaza también cada resultado de validación interna con los FASTQ, BAM o
    CRAM, muestras, grupos de lectura y conjuntos comprobados, y resume las
    desviaciones detectadas y sus explicaciones. Para cada comprobación
    estadística conserva la tolerancia, el valor esperado, el valor realizado y
    la decisión aplicada, e identifica los conjuntos que no quedaron completos.
    Para WGS registra preparación, objetivo de biblioteca, cobertura total y por
    copia y experimentos de control enlazados. Para evolución enlaza el historial
    generacional y registra reproducción, selección, migración, mutación y
    recombinación. Para la población inicial registra el modelo de
    Hardy-Weinberg, F objetivo y realizado, valores esperados y realizados y el
    método independiente o por bloques con LD, con sus haplotipos faseados. Para
    contaminación registra matriz o pesos, sustituciones,
    transferencias, descartes, totales y cobertura final. Para TE enlaza el
    catálogo y la verdad completa de cada evento.
    Cuando la referencia inicial sea suministrada, registra además la procedencia
    original del FASTA y el nombre o la ruta de la copia utilizada. El
    manifiesto es estructurado y legible por el programa. La raíz conserva el
    manifiesto general y una sola copia verificable de cada referencia e índice;
    cada evento conserva su directorio y manifiesto propios, con resultados
    separados por referencia y sumas de comprobación. La exportación de un
    evento copia sus recursos necesarios a un paquete autocontenido y valida su
    integridad antes y después.

La ejecución termina con éxito solo cuando existen un FASTA por referencia y
sus índices necesarios, los estados completos y los historiales obligatorios
para la vía, las copias de todas las configuraciones usadas, los manifiestos
completos y todas las salidas activadas. Toda salida de observación requiere un
evento explícito. Dentro de un evento que solicite alineamientos, CRAM es el
formato predeterminado y BAM la alternativa configurable. En el modo realista,
el FASTQ final o pareja de FASTQ y el fichero de procedencia de cada grupo de
lectura son salidas obligatorias dependientes; en el modo directo no son
necesarios. Un evento puede solicitar únicamente verdad genotípica o VCF sin
tecnología, biblioteca ni mapeador. La ejecución, sus eventos, conjuntos y
artefactos solo se publican como completos conforme a la jerarquía y atomicidad
de REQ-SAL-21 a REQ-SAL-26.

Todas las salidas deben describir la misma simulación, aunque no utilicen una
misma referencia. Dentro de cada conjunto, todos los archivos deben usar
exactamente la misma referencia. El flujo termina cuando se han generado esas
salidas y los manifiestos conservan la información necesaria para relacionarlas.

## 6. Entradas y salidas confirmadas

### 6.1. Entradas

| Entrada | Contenido confirmado | Información pendiente |
|---|---|---|
| Fichero general de configuración | Fichero estructurado con referencia, ploidía global de 1 a 10 y sustituciones por secuencia, vía, cambios, semilla, parámetros, umbral entre clases, cero, uno o varios eventos explícitos, referencias aplicables y salidas. Cada evento que genere FASTQ, BAM o CRAM exige WGS de DNA genómico, una tecnología y perfil, preparación sin PCR o amplificada y exactamente un objetivo de cobertura o número por biblioteca; si genera alineamientos exige además un mapeador. Un evento que solo solicite verdad o VCF no requiere esos parámetros. Configura errores, duplicados, contaminación sustitutiva, matriz o tasa total, bibliotecas, tandas y, cuando corresponda, reproducción, selección, migración, mutación y recombinación. Dentro de un evento con alineamientos, CRAM es predeterminado y BAM la alternativa. En modo realista se activan FASTQ; VCF se solicita expresamente por evento. Sin eventos no se generan salidas de observación | Representación técnica de intervalos; distribuciones adicionales; límites prácticos de recursos; formato y esquema generales |
| Fichero de configuración poblacional | Un fichero estructurado y legible por el programa para cada escenario evolutivo o de familias y poblaciones segregantes. En la vía evolutiva contiene fases, N objetivo, Ne objetivo, política de extinción por N=0, subpoblaciones, migraciones, procesos ordenados, eventos demográficos, inicialización, recombinación y muestreos. División, fundación, migración y selección de descendientes comparten el catálogo declarativo confirmado, sin código personalizado. Cada locus o bloque y subpoblación usa modo probabilístico predeterminado o controlado explícito. La inicialización admite equilibrio de Hardy-Weinberg, distribuciones explícitas, bloques de haplotipos con LD y F global o por subpoblación. En la vía de familias define parentales referenciados o expresos y un grafo dirigido acíclico de eventos. Un parental expreso declara cambios faseados por copia respecto a una referencia compatible; lo no declarado permanece como referencia. Los descendientes posteriores pueden seleccionarse por identificador, al azar, por pertenencia o por estado genético en uno o varios loci, con Y, O y NO, usando un número exacto o una proporción probabilística. Un individuo puede alimentar varias ramas sin duplicar su identidad. Cada población intermedia puede solicitar BAM o CRAM y VCF | Tolerancias técnicas; metadatos de muestreo; formato del estado genotípico intermedio; formato y esquema del catálogo común, de parentales expresos, del grafo y de las reglas de selección |
| FASTA real de la referencia inicial | Secuencia real que se usará como referencia inicial cuando se elija el modo suministrado. Debe ser un FASTA legible y contener al menos una secuencia no vacía. El identificador de secuencia es el primer elemento no vacío de la cabecera, debe ser único y el resto se conserva como descripción. Las secuencias solo pueden contener A, C, G, T y N, en mayúsculas o minúsculas; no se admiten huecos, otros códigos IUPAC ni ningún otro símbolo. El sistema lo valida antes de iniciar la simulación y, si detecta errores, detiene la ejecución sin generar resultados y los explica con claridad. Si es válido, lo copia al directorio de resultados, usa la copia y no modifica el original | Queda pendiente definir el tratamiento operativo de N por cada clase de operación |
| FASTQ externos de contaminación | Fuentes permitidas únicamente en modo realista. Deben ser compatibles con tecnología, perfil y disposición de la biblioteca receptora; se validan y copian antes de usarse, se preservan los pares y no se modifican los originales. En modo directo solo se admiten fuentes simuladas con verdad conocida | Validaciones concretas y datos técnicos de integridad |
| FASTQ real WGS y referencia para calibración | Entrada opcional de lecturas cortas tipo Illumina o largas PacBio HiFi usada para estimar un perfil empírico WGS de longitud, calidades, errores y otros parámetros tecnológicos. La referencia permite intentar separar errores de secuenciación, variantes reales, regiones ambiguas y errores de mapeo | Validaciones y método técnico de estimación |
| Perfil empírico reutilizable | Fichero estructurado producido por una calibración anterior y seleccionable como perfil efectivo sin volver a analizar ni necesitar su FASTQ real de origen | Formato y esquema aplazados a la planificación de la implementación |
| Catálogo FASTA de elementos transponibles | Entrada opcional con familias o tipos de TE. Se valida, copia a resultados y usa desde la copia sin modificar el original. Como alternativa, el simulador genera y conserva el catálogo | Propiedades estructuradas exactas del catálogo y validaciones específicas |

El FASTA real solo es una entrada cuando se elige una referencia inicial
suministrada. El fichero de configuración poblacional solo es una entrada en la
vía evolutiva o en la vía de familias y poblaciones segregantes.

### 6.2. Salidas

| Salida | Contenido confirmado | Información pendiente |
|---|---|---|
| FASTA de las referencias | Un único FASTA por cada referencia inicial o derivada; el archivo contiene uno o varios cromosomas o *scaffolds*. Cuando la referencia inicial es suministrada, este archivo es la copia guardada en el directorio de resultados y utilizada por la ejecución; el original no se modifica | Ninguna decisión funcional pendiente sobre su organización, tratamiento o conservación |
| Índices de las referencias | Ningún índice, uno o varios por referencia, según los tipos que necesiten las operaciones o los formatos utilizados | Tipos de índice exigidos por cada operación y herramientas de indexación |
| Estructura de los individuos y grupos | Individuos con identificadores estables y exactamente 1–10 haplotipos por cromosoma o *scaffold* según su ploidía efectiva común entre individuos. El estado asigna todas las variantes a copias concretas y conserva fase, dosis, linaje y origen de segmentos. Las repeticiones en tándem conservan secuencia, motivo, copias, pureza e interrupciones | Representación técnica de intervalos, distribuciones adicionales, tolerancia alélica, modelos y reglas de cruces; estructura interna e identificadores |
| Fichero obligatorio de verdad | Estado genotípico completo, faseado y reconstruible con alcance mínimo específico por vía: todos los individuos creados en la vía directa; población inicial, final y generaciones o muestreos solicitados en evolución, con resúmenes del resto; y todos los individuos de todos los eventos y poblaciones intermedias en familias. Para lecturas conserva secuencia original frente a observada, posiciones y haplotipos de origen, duplicados, adaptadores, quimeras, anomalías de pareja y categorías sintéticas. Existe aunque no se solicite VCF y sigue siendo la fuente completa cuando este omite posiciones invariantes | Formato, esquema y organización de la verdad de lecturas aplazados a la planificación |
| Copias de las configuraciones | Copia en el directorio de resultados de cada fichero general o poblacional realmente utilizado | Organización y nombres técnicos |
| Directorios de eventos | Un directorio independiente por evento explícito, con su manifiesto, selección de individuos, salidas y subsemilla. Dentro del evento, los resultados BAM, CRAM y VCF se separan por referencia. Los eventos que solo generan verdad o VCF no necesitan tecnología, biblioteca ni mapeador | Nombres y estructura técnica de directorios |
| Población final y muestreos evolutivos | En la vía evolutiva se conserva la población final y, para cada muestreo no destructivo, el estado genotípico faseado y al menos generación, fase evolutiva o escenario, subpoblación, identificadores y referencia estructural. Cada individuo muestreado incluye una salida detallada de puntos de recombinación y segmentos trazables hasta la copia haplotípica fundadora, con coordenadas y orientación relacionables, copias duplicadas e inserciones nuevas identificadas por separado, deleciones trazables y capas de mutación de novo. No es obligatorio conservar poblaciones intermedias completas. Por defecto no se generan BAM, CRAM ni VCF, pero cada evento puede solicitarlos | Metadatos adicionales de cada observación; formato y esquema de la traza detallada |
| Resúmenes poblacional y muestral de cada evento | Dos resúmenes estructurados y separados por muestreo y generación. El poblacional describe la población completa y las subpoblaciones relevantes; el muestral usa solo los individuos exportados. Ambos identifican inequívocamente ámbito, población o subpoblación y número de individuos, e incluyen frecuencias alélicas, heterocigosidad, F, N y Ne cuando correspondan. N muestral no se presenta como N objetivo ni N realizado de la población y ninguna estadística muestral se presenta como poblacional | Formato y esquema técnicos; organización de loci en el resumen |
| Poblaciones intermedias de cruces | Estado genotípico faseado, haplotipos, origen de segmentos y metadatos de todos los individuos de todos los eventos y poblaciones intermedias, aunque no sean seleccionados. Por defecto no se generan BAM, CRAM ni VCF; cada población puede solicitarlos para todos sus individuos con las reglas de conjuntos por referencia | Formato del estado genotípico, metadatos adicionales y sintaxis del catálogo común de selección |
| Historial generacional | Resumen estructurado obligatorio de demografía, procesos, N objetivo y realizado, Ne objetivo y reproductivo realizado, F objetivo y realizado, reproducción, selección, migración, mutación, recombinación, frecuencias, extinciones, eventos ejecutados u omitidos y muestreos de cada generación. Los genotipos completos se conservan para población inicial, final y generaciones solicitadas. El resumen no exporta el detalle individual de segmentos de toda la población | Formato, esquema, fórmulas concretas de Ne reproductivo y F realizado y sus adaptaciones por ploidía |
| Catálogo FASTA de elementos transponibles | Copia validada del catálogo suministrado o FASTA generado, enlazado con familias, copias y eventos de TE de la verdad completa | Esquema estructurado de propiedades del catálogo |
| CRAM o BAM de las muestras | CRAM es predeterminado y BAM la alternativa configurable. En cada conjunto por referencia existe un único archivo final por individuo con todos sus grupos y las degradaciones configuradas; no se guarda además una versión limpia. Un control sin degradaciones es otro experimento enlazado con resultados separados. Se entrega ordenado e indexado, con calidades y duplicados marcados salvo solicitud contraria. Cada CRAM se enlaza a la copia FASTA exacta. En modo directo queda etiquetado como sintético | Perfiles y versiones; cabecera y procedimiento de marcado de duplicados |
| Índice del alineamiento final | Índice CRAI, BAI o CSI compatible con el CRAM o BAM final, su tamaño y sus necesidades. Alineamiento, índice y referencia deben superar una validación conjunta | Selección técnica del tipo y herramientas de indexación |
| FASTQ finales de lecturas | Un archivo, o una pareja para *paired-end*, por individuo, biblioteca y tanda o carril cuando BAM o CRAM usan el modo realista. Cada salida conserva su perfil y parámetros y corresponde a un grupo de lectura del BAM/CRAM combinado; no se concatena como producto principal. Contiene errores, calidades, duplicados y contaminación después del demultiplexado. No se guarda otra copia limpia, no se reutiliza entre eventos y no es necesaria en modo directo | Organización técnica y validaciones concretas |
| Fichero de procedencia de lecturas | Fichero estructurado por cada FASTQ final o pareja. Identifica por lectura o pareja las contaminantes presentes, su fuente real y el mecanismo que las introdujo, sin duplicar secuencias ni calidades. Para lecturas descartadas o sin asignar conserva solo recuentos agregados por muestra y causa u origen conocido. Permite reconstruir la verdad de la contaminación junto con el FASTQ y el manifiesto | Formato y esquema aplazados a la planificación de la implementación |
| Copias de FASTQ externos | Copia de cada FASTQ suministrado como fuente contaminante. La simulación usa estas copias y deja intactos los originales, para que el experimento sea autocontenido y exportable | Organización y datos técnicos de integridad |
| Perfil empírico estimado | Artefacto estructurado y reutilizable con versión de esquema, parámetros, tecnología, referencia y datos de entrada, método y versiones de herramientas, advertencias y medidas de calidad. Permite simular después sin el FASTQ real original | Formato y esquema aplazados a la planificación de la implementación |
| Resultado de validación interna | Uno o varios ficheros estructurados que comparan resultados con la verdad y la configuración. Para WGS incluyen objetivo de biblioteca, cobertura total y por copia, soporte alélico sin exigir balance exacto, sesgos y contaminación; para evolución incluyen inicialización Hardy-Weinberg y LD, F objetivo y realizado, N objetivo y realizado, Ne objetivo y realizado, reproducción, selección, migración, mutación, recombinación y linaje de muestreos. Registran objetivos, observaciones, tolerancias y decisiones. Las incoherencias lógicas detienen; las desviaciones estadísticas advierten salvo elevación a error. Las desviaciones de N y Ne acordadas no detienen por sí mismas | Formato, esquema, métricas y tolerancias concretas |
| VCF de verdad multimuesta por referencia | Salida opcional única por conjunto, construida desde el estado genotípico faseado y nunca inferida desde FASTQ o BAM/CRAM. Contiene todos sus individuos y solo loci donde alguna muestra tiene un alelo no de referencia; no contiene posiciones invariantes, bloques gVCF ni VCF individuales como salidas principales. Cada referencia tiene registros e identificadores de variante independientes, sin correspondencias registro a registro con otros VCF. Se ordena por la referencia y coordenada, se comprime con BGZF y se indexa mediante TBI o CSI por defecto; el texto sin comprimir es configurable. `GT` conserva genotipo, fase y ploidía reales aunque falte evidencia. No contiene `DP`, `AD`, `GQ`, `PL` ni otros campos `FORMAT` observados o inferidos; esas métricas van al informe de validación. Todo registro válido usa `QUAL=.` y `FILTER=PASS`, convención documentada en la cabecera. Un registro que no pueda representarse mediante el VCF y su complemento cuando corresponda, o que no sea validable, detiene el conjunto | Metadatos biológicos o estructurales concretos, criterio de tamaño manejable, normalización, cabecera y versión VCF |
| Complemento estructurado de verdad | Fichero generado solo cuando uno o más casos no pueden representarse sin pérdida en VCF. Conserva secuencia o cambio, haplotipo, fase, coordenadas y genotipo completos. Cada caso comparte con el VCF un identificador estable local al conjunto, y el manifiesto enlaza complemento, VCF, referencia y conjunto. No enlaza variantes de referencias distintas | Formato, esquema y campos concretos aplazados a la planificación |
| Manifiestos general y de eventos | La raíz contiene un manifiesto estructurado con vía, semilla, versión, configuraciones resueltas, referencias y recursos comunes, subsemillas y enlaces a eventos. Cada evento tiene directorio y manifiesto propios con selección de individuos, parámetros técnicos cuando procedan, referencias y salidas separadas por conjunto. Los manifiestos conservan parámetros y resultados biológicos y técnicos, los cuatro niveles separados de reproducibilidad, estados y diagnósticos de ejecución, evento, conjunto y artefacto, validaciones, estimaciones, advertencias, consumo y sumas de comprobación. Los identificadores de variante son locales a cada conjunto y no enlazan referencias distintas | Formato, esquema y algoritmo de suma de comprobación aplazados a la planificación de la implementación |
| Paquete exportado de evento | Copia autocontenida de un evento con su manifiesto, configuración, salidas y las referencias e índices comunes que necesita. Las sumas de comprobación se validan antes y después de exportarlo | Formato técnico del paquete y algoritmo de suma de comprobación |

Solo se conserva el alineamiento final de cada experimento. Los controles sin
degradaciones se generan como experimentos relacionados y separados.

## 7. Criterios de aceptación preliminares

Un criterio de aceptación es una comprobación observable que permite decidir si
un requisito se ha cumplido. Los siguientes criterios pueden deducirse de los
requisitos confirmados. Sus casos de prueba y detalles concretos siguen
pendientes.

| Criterio | Evidencia preliminar | Aspectos pendientes |
|---|---|---|
| **ACE-01. Configuración sin cambiar el código** | El programa puede leer ficheros generales con distinto número, tamaño o composición de cromosomas o *scaffolds* y generar las referencias simuladas correspondientes sin modificar el código. Cuando se selecciona la vía evolutiva o la vía de familias y poblaciones segregantes, el fichero general referencia el fichero de configuración poblacional del escenario y el programa usa ambos | Formato, esquema y casos de prueba |
| **ACE-02. Tamaños del genoma simulado** | Cada cromosoma o *scaffold* tiene el tamaño configurado y el tamaño total corresponde a los tamaños de todas sus secuencias | Casos de prueba y comprobación exacta |
| **ACE-03. Modos de referencia inicial** | El sistema puede simular la referencia inicial desde la configuración o utilizar un FASTA real suministrado. Un FASTA legible, con al menos una secuencia no vacía y cuyo primer elemento no vacío de cada cabecera sea un identificador único supera la validación si sus secuencias contienen exclusivamente A, C, G, T y N en mayúsculas o minúsculas. La descripción restante de la cabecera se conserva por separado. El sistema copia el FASTA al directorio de resultados, usa esa copia y deja el original sin cambios. Un archivo ilegible, sin secuencias no vacías, con identificadores vacíos o repetidos, huecos, otros códigos IUPAC o cualquier otro símbolo hace que la ejecución se detenga antes de generar resultados y explique los errores | Tratamiento operativo de N por clase de operación; casos de prueba detallados |
| **ACE-04. Un FASTA por referencia** | Cada referencia inicial o derivada conserva toda su secuencia en un único FASTA, que contiene uno o varios cromosomas o *scaffolds* | Casos de prueba |
| **ACE-05. Índices de las referencias** | Cada referencia tiene cero, uno o varios índices y dispone de todos los tipos que necesitan las operaciones o los formatos que la usan | Tipos exigidos por cada operación y herramientas de indexación |
| **ACE-06. Creación de referencias estructurales** | Se pueden crear dos o más referencias estructurales a partir de referencias disponibles previamente, incluida una referencia inicial suministrada | Casos de prueba, tamaños y reglas exactas |
| **ACE-07. Diferencias entre referencias** | Cada nueva referencia estructural difiere de su referencia de origen solo por las duplicaciones, deleciones, inversiones y translocaciones solicitadas | Casos de prueba y comprobación de cada cambio |
| **ACE-08. Origen de las referencias** | Para cada referencia estructural se pueden conocer su referencia de origen y los cambios que las distinguen | Forma de representar y consultar la relación |
| **ACE-09. Clases de variante** | Se pueden crear individuos con SNP, *indels*, variantes estructurales y repeticiones en tándem respecto a cualquiera de las referencias disponibles | Casos de prueba por clase |
| **ACE-10. Creación de cambios** | Para SNP, *indels* y variantes estructurales por separado, la configuración permite indicar una cantidad media esperada o una densidad por base, kilobase o megabase; la densidad se convierte en una media usando la longitud válida. El número realizado es entero, puede diferir de la media y se obtiene con Poisson de forma predeterminada o con otra distribución seleccionada por clase y región. Los valores espaciales siguen la jerarquía genoma, cromosoma o *scaffold* y región, junto con las sustituciones por individuo. El ámbito más específico prevalece y una configuración con valores incompatibles en regiones solapadas del mismo nivel se rechaza antes de simular identificando ambas regiones y el parámetro. Los valores alélicos fijados tienen prioridad. El modo probabilístico es predeterminado; el controlado, solo explícito, optimiza conjuntamente recuentos enteros, aplica la tolerancia y resuelve empates con la semilla | Formato técnico de coordenadas; catálogo y parámetros de distribuciones adicionales; pruebas estadísticas de las medias |
| **ACE-11. Ploidía** | La configuración admite un valor global y sustituciones por cromosoma o *scaffold*, cada uno entero de 1 a 10. Toda secuencia tiene un valor efectivo idéntico en todos los individuos. Los valores no enteros, fuera de rango o ausentes se rechazan antes de simular. Estado genotípico y VCF representan exactamente sus haplotipos y dosis | Casos de prueba de límites y sustituciones |
| **ACE-12. Creación directa** | Se pueden crear directamente individuos con cambios concretos o aleatorios indicados en la configuración | Casos de prueba |
| **ACE-13. Vía evolutiva** | Cada escenario usa fases ordenadas y generaciones discretas no solapadas. N objetivo guía la descendencia y N realizado registra el tamaño obtenido; una desviación aleatoria advierte sin detener. Ne objetivo es estadístico y opcional; Ne realizado siempre se calcula como Ne reproductivo. El orden predeterminado es demografía, migración, selección de viabilidad, reproducción, segregación y recombinación, mutación de novo, registro y muestreo. La deriva surge del muestreo finito | Algoritmos, fórmulas, tolerancias técnicas y pruebas de cada proceso |
| **ACE-14. Familias y poblaciones segregantes** | Cada escenario representa los cruces como un grafo dirigido acíclico con fuentes, orden y número de descendientes. Un ciclo es una dependencia imposible de un resultado posterior, no reutilizar un parental: autofecundación, retrocruce y ramas múltiples son válidos cuando cada evento produce descendientes nuevos. Antes de simular se validan identificadores, existencia y precedencia de referencias, creación única de descendientes, cardinalidad parental, compatibilidad y producibilidad de reglas y salidas. Existen plantillas F1, F2, retrocruce, autofecundación y pedigrí simple; los demás diseños se expresan como cruces personalizados. Un parental expreso declara sus cambios por copia respecto a una referencia compatible; lo no declarado permanece como referencia y una fase incompleta se rechaza. La descendencia alimenta eventos posteriores por identificador, selección aleatoria, pertenencia o estado genético en uno o varios loci, combinables con Y, O y NO. Un número exacto usa sorteo sin reemplazo si sobran elegibles y falla si faltan; una proporción es una probabilidad individual y produce un recuento aleatorio reproducible. Todo se evalúa sobre la verdad biológica. Un individuo puede alimentar varias ramas con la misma identidad y genotipo, mientras cada evento forma gametos independientes y descendientes nuevos. Se conservan los estados intermedios y el manifiesto reconstruye eventos y linajes | Sintaxis del grafo, de los parentales expresos y de las condiciones; casos de prueba |
| **ACE-15. Referencia de mapeo** | Cuando están activos BAM o CRAM, la configuración selecciona una o varias referencias de mapeo, incluida una distinta de la estructura de una muestra. Para cada referencia se generan los archivos de todas las muestras incluidas en un conjunto independiente | Cardinalidad de archivos dentro de cada conjunto y comprobaciones técnicas |
| **ACE-16. Cabecera de BAM y CRAM** | La cabecera y sus grupos de lectura permiten identificar sin ambigüedad referencia, muestra, evento, biblioteca, tanda o carril y grupo. `SM` conserva el identificador biológico estable y cada identificador `RG` es único en toda la ejecución. El archivo cumple la versión oficial elegida | Versiones exactas de SAM/BAM y CRAM, catálogo concreto de campos y casos de prueba |
| **ACE-17. Degradaciones** | Se pueden generar datos afectados por errores, contaminación, cobertura baja o cero y mala calidad de mapeo. La cobertura media objetivo admite sustituciones por individuo, secuencia y región con la precedencia acordada; antes de sesgos, la realizada es estocástica y aproximadamente Poisson, no exacta por base, y queda registrada cuando se mide | Tolerancias estadísticas; definición, intensidad y comprobación de las demás degradaciones |
| **ACE-18. VCF de verdad** | Cuando se solicita VCF, existe uno independiente por referencia de mapeo con los genotipos y dosis de las mismas muestras biológicas, expresados contra esa referencia. Cada VCF se construye directamente, usa identificadores de variante locales y no expone correspondencias registro a registro con otros VCF | Casos por referencia y comprobación de espacios de identificadores independientes |
| **ACE-19. Cabecera del VCF** | La cabecera o los metadatos estándar de cada VCF permiten identificar sin ambigüedad la referencia, el conjunto y las muestras mediante sus identificadores biológicos estables, y cumplen la versión oficial elegida | Versión exacta de VCF, catálogo concreto de campos y casos de prueba |
| **ACE-20. Separación por referencia** | Cuando se generan BAM, CRAM o VCF, existe un conjunto independiente por cada referencia de mapeo seleccionada. Cada conjunto usa una sola referencia, todos contienen las mismas muestras solicitadas y ninguno mezcla archivos referidos a genomas distintos. En modo realista, dentro del evento todos usan los mismos FASTQ y el único mapeador del evento | Forma de organizar los conjuntos |
| **ACE-21. Coherencia de las salidas** | Las referencias, los individuos y cada conjunto de BAM o CRAM y VCF de verdad se pueden relacionar sin ambigüedad y describen la misma simulación | Identificadores y comprobaciones de consistencia |
| **ACE-22. Manifiestos** | Cada ejecución genera un manifiesto general y uno por evento explícito. Identifican vía, configuraciones, semillas y subsemillas, versión, parámetros efectivos, referencias y recursos comunes, individuos y todos los archivos. Cada evento registra sus salidas y, cuando proceda, modo, tecnología, perfil, contaminación, FASTQ, mapeador único, versión y comando. También conservan evolución, cruces, ploidía, validaciones, los niveles de reproducibilidad declarados y el estado y diagnóstico de cada unidad | Formatos, esquemas y casos de prueba |
| **ACE-23. Niveles de reproducibilidad** | Toda ejecución usa una semilla. Si la configuración no la proporciona, el programa la genera antes de simular y el manifiesto la identifica como generada. Con iguales semilla, configuración efectiva, entradas, perfiles y versión, se garantiza la misma simulación biológica y se aíslan sus flujos de los eventos técnicos. Las mismas lecturas, calidades, procedencias y artefactos se garantizan además cuando coinciden la subsemilla y las dependencias técnicas aplicables. La equivalencia semántica se declara solo bajo condiciones verificadas y la identidad binaria no se promete de forma general. El manifiesto separa los cuatro niveles y nunca infiere uno superior desde otro inferior | Estrategia técnica de subsemillas; condiciones verificables por herramienta, versión o plataforma; casos de prueba |
| **ACE-24. Muestreos evolutivos** | La configuración puede solicitar muestreos mediante generación y número o proporción. Se conservan por defecto el estado genotípico y los metadatos mínimos, sin BAM, CRAM ni VCF. Cada evento guarda por separado un resumen de la población completa y subpoblaciones relevantes y otro calculado solo desde la muestra. Ambos identifican ámbito e incluyen frecuencias alélicas, heterocigosidad, F, N y Ne cuando corresponden; ninguna medida muestral se presenta como poblacional. El muestreo no modifica individuos ni la trayectoria evolutiva | Formato de resúmenes; estrategia técnica de flujos aleatorios; casos de prueba |
| **ACE-25. Identidad individual** | Todo individuo tiene un identificador único, estable y no reutilizable dentro de la ejecución. Los duplicados aportados por la configuración se rechazan. El mismo identificador relaciona genotipo, generación y subpoblación de origen, parentales o linaje, muestreos, BAM o CRAM, VCF y manifiesto, y se distingue de los identificadores de evento y conjunto | Sintaxis técnica de los identificadores y casos de prueba |
| **ACE-26. Estructura poblacional dinámica** | División reparte progenitores entre hijas exclusivas; fusión reúne orígenes; fundación transfiere fundadores y deja continuar al origen. División y fundación usan el catálogo declarativo común. Cada nueva subpoblación tiene N objetivo independiente y registra conjunto fundador, N objetivo y realizado. En fusión, el modo probabilístico usa elegibles o pesos; el controlado garantiza recuentos parentales exactos. Ninguno garantiza porcentaje exacto de bases: se calcula ascendencia genómica realizada con traza haplotípica | Casos de prueba del catálogo común y de cada evento |
| **ACE-27. Recombinación** | La vía evolutiva y la de familias simulan recombinación al formar descendientes. Toda la ejecución usa una tasa uniforme o mapas completos, sin interferencia. Cada individuo pertenece a cero o a un grupo configurable y usa su perfil específico o el general. Solo recombinan intervalos colineales; en poliploides pares se forman bivalentes al azar entre copias compatibles. Los mapas no válidos y las pertenencias múltiples se rechazan | Formato técnico de mapas y posiciones; casos de prueba estructurales |
| **ACE-28. Fase y haplotipos** | El estado interno separa todas las copias cromosómicas según la ploidía, asigna cada variante a haplotipos concretos y conserva fase y origen de segmentos durante evolución, recombinación, cruces, muestreos y poblaciones intermedias. Los VCF emiten genotipos faseados cuando el registro puede representarlos; cualquier limitación usa un complemento enlazado que conserva la fase sin pérdida | Estructura de datos, algoritmo y casos que requieren complemento |
| **ACE-29. Umbral entre *indels* y variantes estructurales** | Sin configuración se usan 50 pb: una inserción o deleción menor de 50 pb es *indel* y una igual o mayor es variante estructural. Otro umbral configurado se aplica de forma coherente a generación, estadísticas, VCF y manifiesto. Otros tipos estructurales conservan su clase | Casos límite y pruebas de coherencia |
| **ACE-30. Representación VCF** | SNP e *indels* usan REF/ALT explícitos y normalizados respecto a la referencia del conjunto. Varios ALT representables del mismo locus forman un único registro multialélico cuyo `GT` conserva ploidía, dosis y fase. Los solapamientos pequeños se representan desde la secuencia haplotípica final. Las variantes estructurales usan alelos simbólicos y SVTYPE, END y SVLEN cuando corresponde. Las variantes complejas conservan información reconstruible y ninguna pérdida de verdad es silenciosa | Versión exacta VCF, normalización técnica y casos complejos |
| **ACE-31. Repeticiones en tándem** | La configuración admite loci concretos o generación aleatoria y regional con motivo no vacío de longitud libre, copias o longitud, pureza e interrupciones. Se mantienen fase y verdad en secuencias y BAM o CRAM. El VCF usa alelos explícitos cuando son manejables o simbólicos para tramos largos o complejos, pero en ambos casos expone motivo, copias por alelo o haplotipo, longitud, pureza e interrupciones y permite reconstruir el genotipo sin depender solo del manifiesto | Campos VCF concretos, criterio técnico de tamaño manejable, límites de recursos y casos de prueba |
| **ACE-32. Mutación de repeticiones en tándem** | El modelo escalonado predeterminado gana o pierde unidades completas y permite configurar tasa, sesgo y distribución del tamaño del cambio, además de modelos alternativos. Sus parámetros pueden variar por los ámbitos acordados. Las mutaciones internas de interrupciones, pureza o motivo están desactivadas por defecto y solo se aplican si se configuran. Cada evento es reproducible, faseado y trazable | Fórmulas y distribuciones concretas, modelos alternativos y casos de prueba |
| **ACE-33. Verdad complementaria sin pérdida** | Si un caso complejo no puede representarse completamente en VCF, se genera un complemento estructurado con identificador compartido y el VCF y el manifiesto lo enlazan. A partir del VCF y el complemento se reconstruyen sin pérdida la secuencia o cambio, coordenadas, haplotipo, fase y genotipo. Si todos los casos son representables, el complemento no se genera | Formato, esquema y casos de reconstrucción |
| **ACE-34. Resultado mínimo y salidas dependientes** | Toda ejecución exitosa produce FASTA e índices necesarios, la verdad completa obligatoria para su vía, los historiales o resúmenes exigidos, copias de configuraciones y manifiesto. Toda salida de observación requiere un evento explícito y cero eventos no crea ninguno implícito. Un evento puede solicitar solo verdad o VCF sin parámetros técnicos. Cuando solicita alineamientos, CRAM es predeterminado y BAM la alternativa; en modo realista también son obligatorios sus FASTQ finales y ficheros de procedencia, mientras el modo directo no los requiere. Una calibración produce su perfil empírico reutilizable. La ejecución falla si falta una salida mínima o activada | Formato del fichero de verdad, organización de copias y casos de fallo |
| **ACE-35. Escala y prioridad de la primera versión** | Se pueden configurar, sin límites funcionales codificados, distintos números y tamaños de secuencias, individuos y generaciones. La validación prioriza escenarios pequeños y reproducibles. Los escenarios mayores están limitados por los recursos disponibles y no se presupone rendimiento ilimitado | Casos pequeños prioritarios, mediciones de recursos y pruebas de escalabilidad sin umbrales fijados |
| **ACE-36. Estimación de recursos** | Antes de ejecutar se muestran estimaciones aproximadas de volumen y almacenamiento de salidas, memoria y tiempo. Un escenario grande provoca una advertencia clara, pero no se bloquea solo por ella. Una configuración inválida o la imposibilidad inevitable de producir el mínimo sí impiden iniciar y explican la causa. El manifiesto conserva estimaciones, advertencias y consumo real medible | Modelos, umbrales y pruebas de calibración de las estimaciones |
| **ACE-37. Tecnologías de lectura** | Con BAM o CRAM activos se pueden representar lecturas cortas tipo Illumina, *single-end* y *paired-end*, y lecturas largas PacBio HiFi. Cada modalidad tiene un perfil predeterminado validado y ajustable; también admite calibración empírica desde FASTQ y referencia. En modo realista las lecturas proceden de haplotipos faseados y resultados y manifiesto registran tecnología y perfil efectivos. Oxford Nanopore y PacBio CLR se rechazan en v1 y quedan para futuras ampliaciones | Distribuciones, parámetros, tasas y valores cuantitativos exactos; datos de validación y casos de prueba |
| **ACE-38. Modos de alineamiento** | En el modo realista predeterminado se guardan los FASTQ finales después de aplicar errores, contaminación y demultiplexado, y los BAM o CRAM contienen esas mismas lecturas, incluidos pares, identificadores y registros no mapeados producidos por el mapeador. En el modo directo no se generan FASTQ previos y los alineamientos estándar quedan claramente etiquetados como sintéticos | Organización y casos de prueba por modo |
| **ACE-39. Un mapeador por experimento** | Un experimento usa exactamente un mapeador, minimap2 salvo selección compatible distinta, contra todas sus referencias. Para minimap2, el perfil elige automáticamente el *preset* y los valores predeterminados, y las sustituciones se advierten y validan. El manifiesto conserva herramienta, versión, comando completo y parámetros efectivos | Versión de minimap2, *presets* y valores exactos, catálogo inicial de alternativas, comprobaciones técnicas de compatibilidad y casos de opciones incompatibles |
| **ACE-40. Independencia de eventos técnicos** | Dos eventos pueden usar mapeadores o tecnologías distintos sobre las mismas muestras biológicas, pero cada uno genera sus propias lecturas y usa un solo mapeador. No se reutilizan FASTQ entre eventos ni se presentan como comparación de mapeadores | Organización, identificadores y casos con varios eventos |
| **ACE-41. Cobertura configurable** | La cobertura media objetivo en «x» puede configurarse globalmente y sustituirse por individuo, cromosoma o *scaffold* y región, incluidos valores altos, bajos o cero. Los conflictos del mismo nivel se rechazan; la realización es aproximadamente Poisson antes de sesgos, se refleja en FASTQ del modo realista o BAM/CRAM del directo y se registra cuando se mide | Tolerancias y casos estadísticos |
| **ACE-42. Contaminación trazable** | Una muestra admite una fracción total de contaminación biológica entre 0 % y 100 %, repartida entre una o varias fuentes mediante matriz dirigida o selección automática. Proporciones, coberturas equivalentes y números se convierten conjuntamente a lecturas o parejas completas y se asignan enteramente minimizando la desviación global; los empates se resuelven con la subsemilla y se registran. Cada unidad contaminante sustituye una propia; la donante no pierde datos. El modo directo aplica reglas funcionalmente equivalentes. Verdad y manifiesto conservan origen, objetivo, conversión, realizado, retiradas y totales finales | Validaciones y casos de prueba |
| **ACE-43. Contaminación *paired-end*** | En datos *paired-end* se seleccionan y sustituyen parejas completas de una misma fuente. Ambos *mates* se incorporan juntos con nombres, indicadores, orientación y relación coherentes, sin crear huérfanos; el manifiesto registra parejas propias retiradas y contaminantes incorporadas | Casos de prueba de emparejamiento |
| **ACE-44. FASTQ contaminante externo** | Los FASTQ externos solo se admiten en modo realista. Se validan y copian al directorio de resultados antes de usarse. La ejecución usa la copia, conserva intacto el original y realiza un submuestreo reproducible que respeta la fracción y el emparejamiento. El manifiesto relaciona ubicación original, copia, validación y lecturas seleccionadas. Configurarlos en modo directo es un error previo | Validaciones exactas y datos técnicos de integridad |
| **ACE-45. Contaminación técnica por índices** | Una tanda multiplexada puede usar índices simples, duales combinatorios o UDI y simular de forma reproducible asignaciones erróneas, incluido *index hopping*. Una lectura o pareja reasignada desaparece de la donante y entra en la receptora sin duplicarse. Las descartadas o sin asignar no se reemplazan y reducen la cobertura. Se conserva la procedencia verdadera, sin generar señales, BCL ni lecturas i5/i7 detalladas | Casos de prueba por estrategia de índices y resultado |
| **ACE-46. FASTQ final, identidad y procedencia** | En modo realista existe un FASTQ final, o pareja *paired-end*, por individuo, biblioteca y tanda o carril después de contaminación y demultiplexado, sin concatenación ni copia limpia adicional. Cada lectura o pareja es una instancia técnica única dentro del evento. Una copia biológica recibe identidad nueva y enlaza su origen; una transferencia por demultiplexado conserva la identidad y abandona la donante; solo un mecanismo de duplicación configurado puede crear copias adicionales. Un fichero estructurado enlazado conserva esas relaciones sin copiar secuencias ni calidades. Las descartadas o sin asignar aparecen solo como recuentos agregados por muestra y causa u origen conocido | Formato, esquema y casos de reconstrucción de la verdad de contaminación |
| **ACE-47. Perfiles de error y calidades** | Cada tecnología usa un perfil de errores identificado y versionado, homogéneo salvo sustituciones globales o por individuo, secuencia o región. Las calidades del FASTQ son probabilísticamente coherentes con los errores y se propagan a BAM/CRAM. La verdad distingue secuencia original y observada | Perfiles exactos, modelos por posición o contexto y casos estadísticos |
| **ACE-48. Calibración empírica reutilizable** | Un FASTQ real y su referencia permiten estimar longitud, calidad y errores, con advertencias sobre variantes, regiones ambiguas, mapeo y fiabilidad. El perfil solo se aplica si se selecciona y se guarda como artefacto estructurado, versionado y reutilizable sin el FASTQ original; el manifiesto lo enlaza | Método de estimación, formato, esquema y medidas concretas de calidad |
| **ACE-49. Duplicados de PCR y ópticos** | Los perfiles pertinentes generan por separado duplicados de PCR y ópticos con tasas sustituibles; en *paired-end* se duplica la pareja completa. La verdad conserva tipo y origen. BAM/CRAM los entrega marcados por defecto y permite una salida sin marcar; el manifiesto registra estado, método y cantidades esperadas y realizadas | Modelos exactos, procedimiento de marcado y casos de prueba |
| **ACE-50. Restos de adaptadores** | Los restos de adaptadores están desactivados por defecto y solo se añaden por solicitud expresa con secuencia o perfil, frecuencia, longitud residual y extremo efectivos. La verdad identifica lectura e intervalo añadido y el manifiesto conserva estado, parámetros y cantidad realizada | Perfiles y casos de prueba concretos |
| **ACE-51. Sesgo de cobertura por GC** | El perfil puede relacionar contenido GC y probabilidad de lectura; la configuración modifica su intensidad o lo anula. El efecto se combina de forma documentada con la cobertura objetivo y las excepciones regionales, y el manifiesto conserva parámetros y resúmenes de cobertura realizada | Curvas, regla de combinación e intervalos de resumen concretos |
| **ACE-52. Baja mapeabilidad realista** | Las lecturas de regiones repetitivas, duplicadas u homólogas parten de posiciones y haplotipos verdaderos. El mapeador determina alineamientos primarios, secundarios, posiciones y MAPQ sin sobrescritura artificial; la verdad permite detectar ambigüedad y posiciones incorrectas | Casos de prueba por tipo de región y referencia |
| **ACE-53. Alineamientos sintéticos coherentes** | En modo directo, las proporciones de lecturas correctas, mal alineadas, multimapeadas y no alineadas y sus distribuciones de MAPQ son configurables. Se rechazan proporciones, MAPQ o combinaciones inválidas; banderas, CIGAR, campos y alineamientos secundarios o suplementarios concuerdan con la categoría y la verdad conserva origen y categoría | Distribuciones predeterminadas y casos de validación por formato |
| **ACE-54. Lecturas quiméricas** | Un perfil puede generar quimeras con frecuencia configurable. Todos los fragmentos pertenecen al mismo individuo y biblioteca, pero pueden proceder de cromosomas, *scaffolds*, haplotipos o copias diferentes. La verdad conserva origen, coordenadas, orientación y uniones; el mapeador resuelve el modo realista y el modo directo produce alineamientos divididos coherentes | Distribuciones de fragmentos y casos de prueba |
| **ACE-55. Anomalías *paired-end*** | Se pueden configurar parejas con orientación anómala, inserto inesperado o un miembro no alineado, incluidas combinaciones. Las normales siguen el perfil; en modo realista decide el mapeador y en modo directo banderas, RNEXT/PNEXT, TLEN, CIGAR, orientación y *proper-pair* son coherentes. La verdad y el manifiesto conservan estado, fragmento y cantidades | Distribuciones y casos de prueba por combinación |
| **ACE-56. Grupos de lectura y archivo por individuo** | Los identificadores únicos de muestra, biblioteca, tanda o carril, plataforma y grupo son consistentes entre FASTQ, `@RG`/`RG` de BAM o CRAM, procedencia y manifiesto. Los FASTQ permanecen separados por grupo. Cada conjunto por referencia contiene un único BAM o CRAM final por individuo con todos sus grupos, y nunca mezcla referencias ni exige productos finales separados por biblioteca | Campos exactos de `@RG`, validaciones y casos con múltiples grupos |
| **ACE-57. Una tecnología por experimento** | Cada experimento o evento que genere FASTQ, BAM o CRAM usa una sola tecnología experimental y mantiene compatibles todos sus perfiles, bibliotecas, tandas, carriles, mapeo y resultados. No mezcla lecturas cortas tipo Illumina y lecturas largas PacBio HiFi ni perfiles incompatibles. Un evento que solo materializa verdad o VCF no necesita tecnología. Eventos separados pueden usar las mismas muestras y verdad, pero producen lecturas independientes con identificadores y versiones trazables | Sintaxis de identificadores y casos de incompatibilidad |
| **ACE-58. CRAM predeterminado y referencia exacta** | Dentro de un evento que solicite alineamientos, cada individuo recibe un CRAM final salvo que se seleccione BAM como alternativa. Cada CRAM está ligado por cabecera y manifiesto a la copia FASTA exacta almacenada entre los recursos comunes, protegida por identificadores y comprobaciones de integridad; no puede asociarse a otra referencia. Una exportación del evento sí copia ese FASTA al paquete autocontenido | Campos, comprobaciones de integridad y casos de referencia incorrecta |
| **ACE-59. Orden, índice y validación** | Cada BAM o CRAM final se entrega ordenado por coordenadas de forma predeterminada y con un índice compatible elegido según formato, tamaño y necesidad. La configuración y el manifiesto registran formato, orden e índice, y el resultado no se considera completo si alineamiento, índice y referencia no superan la validación conjunta | Herramientas, tipos exactos de índice y casos límite de tamaño |
| **ACE-60. VCF multimuesta por referencia** | Cada conjunto solicitado contiene un único VCF de verdad multimuesta con todos sus individuos y exactamente los mismos identificadores usados en BAM/CRAM, grupos de lectura, genotipos y manifiesto. No incluye individuos de otra referencia y no exige VCF individuales como productos principales | Casos de coherencia de identificadores y organización |
| **ACE-61. VCF ordenado, comprimido e indexado** | Por defecto, el VCF sigue el orden de secuencias de la referencia y las coordenadas dentro de cada una, está comprimido con BGZF y tiene un índice TBI o CSI compatible. Puede solicitarse texto sin comprimir. La ejecución falla si VCF, compresión, índice, cabecera y referencia no superan la validación correspondiente | Selección TBI/CSI, herramientas y casos límite |
| **ACE-62. Solo loci diferentes de la referencia en VCF** | El VCF contiene un registro cuando al menos una muestra del conjunto tiene un alelo no de referencia, aunque ese ALT esté fijado en todas las muestras. No emite posiciones idénticas a la referencia ni bloques gVCF por defecto; el fichero estructurado de genotipos y haplotipos conserva la verdad completa | Casos invariantes de referencia, ALT fijado, loci polimórficos, multialélicos y distintas ploidías |
| **ACE-63. VCF independiente de la evidencia observada** | El VCF se construye desde el estado genotípico faseado de la simulación, nunca mediante inferencia desde FASTQ, BAM o CRAM. En cada locus incluido conserva el genotipo completo y la fase real de todas las muestras aunque tengan cobertura cero, contaminación, errores, calidad baja o mapeo ausente, incorrecto, múltiple o ambiguo; esos casos no producen genotipos ausentes. Una limitación de representación se conserva mediante el complemento de verdad y no se confunde con falta de evidencia | Casos con cobertura cero y degradaciones combinadas; comprobación de independencia respecto a FASTQ y BAM/CRAM |
| **ACE-64. Validación interna contra la verdad** | Siempre que existen FASTQ o BAM/CRAM, existe un resultado estructurado enlazado desde el manifiesto. Sin ejecutar un *SNP caller*, compara con la verdad de origen la procedencia desde haplotipos, cobertura, soporte alélico según dosis y ploidía, tasas de errores, contaminación, duplicados y demás degradaciones, coherencia de parejas, y alineamientos, referencias y grupos de lectura. Registra objetivos, valores observados, tolerancias y desviaciones, y toda desviación intencionada queda explicada por la configuración o una decisión aleatoria trazable | Tolerancias, métodos estadísticos y formato y esquema del informe |
| **ACE-65. Consecuencias de la validación** | Un archivo inválido, una referencia incorrecta, identificadores incompatibles, un genotipo imposible o cabeceras o índices incoherentes detienen la ejecución. Una desviación estadística puede elevarse a error, salvo las de N y Ne, que se advierten sin detener. N realizado=0 aplica la política específica: desactivada detiene inmediatamente toda la simulación sin omisiones; activada registra extinción por resultado y permite omisiones dependientes. Informe y manifiesto registran esperado, realizado y decisión | Casos de prueba de cada clase de fallo y sintaxis de la política configurable |
| **ACE-66. Contenido genotípico del VCF verdadero** | En cada locus incluido, `GT` coincide con el genotipo real simulado y conserva fase y ploidía aunque la cobertura sea cero o las lecturas estén degradadas. El VCF no contiene `DP`, `AD`, `GQ`, `PL` ni ningún otro campo `FORMAT` derivado de observaciones o inferencias. La cobertura, el soporte alélico, las calidades y métricas equivalentes aparecen en los resultados de validación de FASTQ/BAM/CRAM, no en el VCF | Casos con distintas ploidías, cobertura cero y degradaciones; inspección de campos `FORMAT` |
| **ACE-67. QUAL y FILTER del VCF verdadero** | Todos los registros válidos usan `QUAL=.` y `FILTER=PASS`. La cabecera explica que no se estima una calidad probabilística y que `PASS` identifica un registro verdadero validado estructuralmente. Si un registro no puede representarse mediante el VCF y su complemento cuando corresponda, o no puede validarse, se produce un error estructural, se detiene la creación del conjunto, este no se marca como completo y el registro no se emite como `PASS` | Casos de cabecera válida, complemento enlazado, registro no representable y fallo de validación |
| **ACE-68. Alcance WGS y preparación** | La primera versión acepta exclusivamente WGS de DNA genómico con lecturas cortas tipo Illumina o largas PacBio HiFi en experimentos separados. Rechaza Oxford Nanopore, PacBio CLR, paneles, amplicones, GBS y captura o reducción de representación. La preparación predeterminada es sin PCR; una preparación WGS amplificada puede generar duplicados de PCR y sesgos configurados y queda identificada en configuración y manifiesto | Casos de rechazo por tecnología y protocolo |
| **ACE-69. Objetivo y muestreo de bibliotecas WGS** | Cada biblioteca declara cobertura o número de lecturas, pero no ambos ni ninguno; en *paired-end* se cuentan parejas. Se deriva el indicador equivalente. La cobertura usa bases de lectura divididas por longitud haploide. Las posiciones se muestrean uniformemente y de forma reproducible y la profundidad previa a sesgos es aproximadamente Poisson, sin profundidad exacta por base | Tolerancias, prueba estadística y modelos futuros alternativos |
| **ACE-70. Ploidía, copias y soporte alélico WGS** | Sin excepciones, las lecturas se reparten por longitud y número de copias, cada copia recibe cobertura semejante y la profundidad cromosómica total es proporcional a su ploidía. Los alelos contribuyen según dosis; la validación admite variación aleatoria y no exige balance exacto. No existen parámetros de sesgo alélico o *allele dropout* forzados | Tolerancias y resúmenes por copia o haplotipo |
| **ACE-71. Un único alineamiento final** | Cada evento conserva solo su BAM/CRAM final con las degradaciones configuradas. Un control sin degradaciones es otro evento que usa la misma verdad, referencias y muestras, pero tiene subsemilla, lecturas y archivos propios; los manifiestos enlazan su relación experimental | Identificadores técnicos de eventos relacionados |
| **ACE-72. Contaminación biológica sustitutiva** | Cada lectura o pareja contaminante sustituye una propia de la receptora sin reducir los datos de la donante. La fracción total admite 0–100 %; varias fuentes no pueden superar el 100 % ni los datos disponibles. El 100 % representa sustitución completa. Exceder límites detiene antes de generar | Unidades exactas de matriz y casos límite |
| **ACE-73. Fuentes y reparto de contaminación** | Las fuentes son compatibles con tecnología, perfil y disposición. Una matriz dirigida puede definir aportes no simétricos; su diagonal es propia. Sin matriz, donantes compatibles de la tanda se ponderan por datos disponibles. La matriz prevalece. Los FASTQ externos validados solo se admiten en modo realista; en directo son error | Validaciones técnicas de FASTQ y casos de matriz |
| **ACE-74. Transferencia y descarte por índices** | Una lectura o pareja reasignada por índices sale de la donante y entra en la receptora sin duplicarse. Las descartadas o sin asignar no se reemplazan y reducen cobertura. El modo directo aplica reglas equivalentes mediante alineamientos sintéticos. Informe y manifiesto conservan objetivos, transferencias, descartes, totales y cobertura final | Casos por estrategia de índices y concordancia funcional entre modos |
| **ACE-75. Plantillas y sistemas de cruce** | Existen plantillas F1, F2, retrocruce, autofecundación y pedigrí simple. F2 permite autofecundar F1 o cruzar F1. RIL, introgresión y MAGIC se expresan con cruces personalizados. El sistema admite especies autógamas y alógamas y registra modo, parentales y descendientes | Sintaxis del sistema general y casos por plantilla |
| **ACE-76. Apareamiento mixto y contribución parental** | La autofecundación es una probabilidad reproducible de 0–100 % por evento, sustituible por población, subpoblación y fase. Los demás eventos eligen uniformemente dos supervivientes distintos, sin ponderación por parentesco, distancia o genotipo. Un parental puede reproducirse varias veces bajo límites o distribuciones configurados | Distribuciones y tolerancias de tasas realizadas |
| **ACE-77. N, Ne y descendencia** | N objetivo es obligatorio al iniciar y fija la descendencia que se intenta crear; un cuello de botella lo reduce y una expansión lo aumenta, sin modificar adultos ya creados. Si un resultado aleatorio impide alcanzarlo, se conservan N realizado, desviación y aviso sin detener. Ne objetivo puede ser absoluto o Ne/N o faltar; Ne realizado usa progenitores con descendencia y desigualdad de contribuciones, nunca heterocigosidad, endogamia o LD. Su desviación tampoco detiene | Fórmulas de Ne reproductivo y F realizado, adaptación por ploidía y casos incompatibles |
| **ACE-78. Generaciones e historial reproducible** | Las generaciones son discretas y no solapadas. Los individuos parentales solo existen como población activa en su generación y los sustituye la descendencia; un cambio de N no los elimina retrospectivamente. Se conservan genotipos completos de la población inicial, final y generaciones o muestreos solicitados, y un historial resumido de todas las demás. La regeneración biológica exige misma semilla, configuración efectiva, entradas, perfiles y versión; los niveles técnicos superiores aplican además las condiciones de ACE-23 | Esquema del historial y política técnica de compatibilidad |
| **ACE-79. Mutaciones evolutivas** | SNP, inserciones, deleciones y SV usan tasas independientes por base y generación; las repeticiones usan su clase y modelo propios. SNP admite espectro validado y modelo simétrico predeterminado, sin contexto vecino. Indels tienen distribuciones de longitud que favorecen tamaños cortos y respetan el umbral. Cada mutación surge después de recombinación en un haplotipo descendiente y conserva origen y transmisión | Distribuciones cuantitativas, tasas predeterminadas y casos estadísticos |
| **ACE-80. TE configurados, no evolución autónoma** | Un catálogo TE validado procede de FASTA o simulación. Se insertan TE al construir referencias o mediante eventos exactos o aleatorios en individuos concretos. No existen tasas de movilización, actividad ni transposición espontánea. Las inserciones presentes se heredan normalmente y se representan en verdad y VCF/complemento | Esquema del catálogo y casos de herencia |
| **ACE-81. Construcción trazable de eventos TE** | Un evento solicitado admite orientación, truncamiento, duplicación de sitio diana, copia completa, copiar-pegar o cortar-pegar y reglas de destino. Conserva secuencias fuente/final, copias, coordenadas, haplotipo, sitio diana y linaje. Sin destino válido advierte o falla según política, pero no crea una inserción inválida | Formato de reglas, pesos y casos complejos VCF/complemento |
| **ACE-82. Selección de viabilidad** | Toda la ejecución usa como máximo un locus estable, bialélico o multialélico, que puede ser SNP, indel, repetición, TE o SV con genotipo/dosis inequívocos para ploidías 1–10. Cada genotipo o combinación de dosis recibe un valor genético y el valor del carácter coincide con él. Una función direccional configurable favorece valores mayores o menores y los transforma en aptitudes relativas no negativas, que se normalizan y producen supervivencia estocástica. Se rechazan funciones estabilizadoras o disruptivas. No se admiten ruido ambiental ni heredabilidad. Los valores pueden variar por fase o subpoblación sin cambiar el locus. No hay selección poligénica, epistasis ni fecundidad directa | Regla para genotipos no enumerados y casos por ploidía, dirección y función de aptitud |
| **ACE-83. Extinción y migración** | La extinción planificada invalida dependencias posteriores detectables. Con N realizado=0, la política global desactivada puede sustituirse por fase y subpoblación: desactivada detiene toda la simulación; activada registra extinción por resultado y omite atómica y transitivamente eventos dependientes con cadenas de causas, mientras continúan ramas independientes. Recolonizar crea otro ID. La migración simultánea usa pares proporcionales o exactos, sin reemplazo ni A→B→C, cambia conjuntos reproductores pero no N objetivo | Reglas técnicas de asignación conjunta y casos de prueba |
| **ACE-84. Deriva y orden generacional** | La deriva emerge del muestreo finito y no existe un paso artificial adicional. El orden predeterminado es demografía, migración, viabilidad, reproducción, segregación/recombinación, mutación de novo, registro y muestreo. Varios eventos demográficos de una generación exigen un orden total explícito; cada uno recibe el estado del anterior. Un orden ausente, ambiguo o imposible se rechaza antes de simular | Casos de orden válido e inválido |
| **ACE-85. Recombinación, bivalentes y estructura** | Toda la ejecución usa tasa uniforme en cM/Mb o mapa físico/genético completo. Sin interferencia, los entrecruzamientos por bivalente siguen Poisson según la longitud genética. Cada cromosoma recibe al menos uno. Un *scaffold* lo recibe si supera un umbral general, sustituible por secuencia, expresado exclusivamente en cM o bases; manifiesto y validación conservan regla, valor y procedencia efectivos. Solo recombinan intervalos colineales. Los poliploides pares forman al azar bivalentes compatibles, sin preferencia de subgenoma, multivalentes ni doble reducción | Método exacto de posiciones, validación de mapas y casos estructurales |
| **ACE-86. Linaje detallado de individuos muestreados** | Solo los individuos incluidos en muestreos exportan los puntos de recombinación y el origen de cada segmento. Cada segmento se reconstruye hasta el individuo y la copia haplotípica fundadores. Las mutaciones de novo forman capas separadas con identidad, individuo y generación de aparición. Los resúmenes generacionales no exportan el detalle de todos los individuos | Formato y esquema de la salida detallada; casos de reconstrucción |
| **ACE-87. Trazabilidad y límites estructurales** | Las coordenadas actuales se relacionan con las fundadoras cuando existe correspondencia, conservando orientación. Cada copia duplicada y cada inserción nueva tiene identificador estable; las deleciones conservan el intervalo fundador ausente. En una copia haplotípica se permiten eventos estructurales no solapados, incluso adyacentes, pero se rechaza solapamiento, anidamiento o modificación de una estructura previa. La única excepción entre fuente y destino de un mismo evento es la duplicación cuya fuente abarca la secuencia completa. Haplotipos distintos pueden solaparse. SNP e indels se permiten en segmentos modificados, salvo si solapan o atraviesan un punto de ruptura; la adyacencia se permite | Formato de relaciones estructurales y casos de validación |
| **ACE-88. Modelos de población inicial** | Cada locus o bloque y subpoblación usa modo probabilístico predeterminado o controlado explícito sin mezclarlos dentro de la unidad. Admite Hardy-Weinberg en ploidías 1–10, frecuencias explícitas de genotipo o dosis y genotipos multilocus faseados. En Hardy-Weinberg controlado se fijan recuentos alélicos, pero los genotipos se distribuyen aleatoriamente condicionados al modelo. Frecuencias y recuentos suman, respectivamente, 1 y el total aplicable; los datos incompatibles se rechazan | Tolerancias numéricas y casos por ploidía y modo |
| **ACE-89. LD y fase iniciales** | Por defecto, la población inicial asigna loci independientemente. El LD usa bloques no solapados con haplotipos y frecuencias conjuntas probabilísticas por defecto o recuentos controlados por subpoblación que suman C. Las frecuencias haplotípicas son la fuente de verdad y de ellas se derivan las alélicas; una frecuencia alélica aportada solo valida y una discrepancia fuera de tolerancia es error. Los recuentos se distribuyen entre individuos condicionados al modelo o plantillas. Los bloques conservan fase e identificador a lo largo de las generaciones, sin generar ni calcular D, D' o r² | Formato técnico de bloques, tolerancias y casos de prueba |
| **ACE-90. Endogamia inicial** | F se configura para la población completa o se sustituye por subpoblación, siempre entre 0 y 1 y nunca para ploidía 1, locus, región o bloque. Para cualquier ploidía mayor que 1 progresa desde Hardy-Weinberg en F=0 hasta homocigosis completa en F=1, sin modelos específicos de herencia polisómica o disómica. No fuerza ROH ni un valor exacto por individuo. Una distribución local explícita prevalece. F realizado compara heterocigosidad observada y esperada por subpoblación y generación; puede ser negativo sin constituir error | Fórmula exacta de F realizado, adaptación por ploidía y tolerancias |
| **ACE-91. Fundadores trazables** | Todos los individuos de la población inicial son fundadores con identificadores estables propios. Cada copia haplotípica también tiene identificador propio aunque su secuencia coincida con otra, y toda traza posterior alcanza una copia fundadora concreta | Formato de identificadores y casos de reconstrucción |
| **ACE-92. Selección de muestreos** | Un evento admite identificadores explícitos, selección aleatoria sin reemplazo o ambas, además de estratos no solapados por subpoblación, grupos y filtros de genotipo. Los explícitos se incorporan primero y se excluyen del sorteo. Duplicados, individuos incompatibles, estratos ambiguos o tamaños mayores que los elegibles detienen el evento sin reemplazo ni reducción silenciosa | Sintaxis de filtros y estratos y casos de prueba |
| **ACE-93. Eventos explícitos independientes** | Una ejecución biológica admite cero, uno o varios eventos explícitos con identificador, selección, referencias, salidas y subsemilla. Solo los eventos que generan FASTQ, BAM o CRAM necesitan una tecnología y perfil; los que generan alineamientos necesitan además un mapeador. Un evento que solo materializa verdad o VCF carece de esos requisitos técnicos. El mismo individuo conserva genotipo, fase y ascendencia entre eventos, pero las lecturas y alineamientos se materializan de nuevo e independientemente | Estrategia técnica de subsemillas y casos multievento |
| **ACE-94. Organización, atomicidad y exportación** | La raíz conserva una sola copia verificable de cada referencia e índice y un manifiesto general. Cada evento tiene directorio y manifiesto propios y separa resultados por referencia. Ejecución, evento, conjunto y artefacto distinguen estados incompleto, completo y fallido. Cada artefacto se publica atómicamente tras validarse; un conjunto solo está completo con todos sus artefactos activados y el éxito global exige todas las unidades obligatorias completas. Los conjuntos completos independientes se conservan y pueden exportarse declarando cualquier fallo global. La exportación produce un paquete autocontenido y valida las sumas antes y después | Esquema, algoritmo de suma y formato del paquete |
| **ACE-95. Coordenadas y límites estructurales** | La configuración usa coordenadas basadas en 1 e intervalos cerrados; cada formato externo usa su convención oficial mediante conversión documentada. Cada duplicación usa el estado inmediatamente anterior para interpretar origen y destino y conserva las correspondencias posteriores. Una fuente completa puede insertar la copia dentro de la misma secuencia; una fuente parcial solo puede usar un destino adyacente o externo a su intervalo. Inserciones, deleciones, inversiones y puntos de ruptura respetan las reglas terminales, y toda translocación conserva subeventos, segmentos y orientaciones sin vaciar una secuencia | Representación técnica y casos límite |
| **ACE-96. Límite local de copias** | El máximo local predeterminado y provisional es 10 copias totales, incluida la original, por familia de copia e intervalo fundador dentro de cada copia haplotípica. Se contabiliza de forma acumulada y separada de ploidía y repeticiones en tándem. Un valor mayor requiere configuración explícita, genera una advertencia y activa una estimación previa de memoria, almacenamiento y volumen de datos. Si supera el presupuesto configurado se rechaza; si no, continúa. No existe otro máximo absoluto fijo y el valor predeterminado debe revisarse mediante pruebas de rendimiento | Pruebas de rendimiento y casos acumulativos y de presupuesto |
| **ACE-97. Reproducción de familias** | Cada evento tiene un parental para autofecundación o dos para cruzamiento sexual; otras cardinalidades se rechazan. El orden de los dos parentales no cambia el resultado biológico y los diseños con tres o más fundadores forman un grafo de eventos. Clonación, apomixis y dobles haploides se rechazan | Casos por tipo de evento y grafo |
| **ACE-98. Segregación y mutación familiar** | La transmisión elige con igual probabilidad una copia compatible de cada bivalente, sin conversión génica, distorsión de segregación ni impulso meiótico. Un heterocigoto para una translocación no puede reproducirse. Las mutaciones de novo están desactivadas por defecto y, al activarse, usan el motor evolutivo y conservan tasa, evento y traza | Casos mendelianos, estructurales y de mutación |
| **ACE-99. Independencia VCF entre referencias** | Cada referencia produce un VCF construido directamente en su secuencia y coordenadas, con un espacio propio de identificadores. Muestras biológicas pueden compartir identificador, pero VCF, complemento y manifiestos no contienen identificadores globales ni correspondencias registro a registro entre variantes de referencias diferentes | Casos con dos referencias y revisión de metadatos |
| **ACE-100. Identidad técnica y biológica** | `SM` de BAM/CRAM y las columnas del VCF usan el identificador biológico estable. El grupo de lectura es único en toda la ejecución, distingue evento, biblioteca, tanda y carril y coincide entre FASTQ, BAM/CRAM, procedencia y manifiesto. Cada instancia de lectura o pareja tiene identidad única dentro del evento y su procedencia distingue copia biológica, transferencia técnica y duplicación configurada sin colisiones | Campos de `@RG`, representación técnica de los vínculos y casos de colisión |
| **ACE-101. VCF estable entre eventos** | Dos eventos con los mismos individuos, referencia y verdad producen VCF con idénticos identificadores locales, coordenadas, REF, ALT y genotipos aunque cambien tecnología, subsemilla, FASTQ, BAM/CRAM o degradaciones | Casos con eventos técnicos diferentes |
| **ACE-102. ALT fijados y multialélicos** | Un ALT fijado en todas las muestras se incluye porque difiere de la referencia. Varios ALT representables en el mismo locus se reúnen en un registro multialélico normalizado y `GT` conserva ploidía, dosis y fase. Los casos no representables usan el complemento sin pérdida | Casos fijados, multialélicos, poliploides y solapados |
| **ACE-103. Variantes pequeñas solapadas** | SNP e *indels* pueden solaparse en una copia haplotípica. Los explícitos exigen un orden, se aplican secuencialmente y usan el estado inmediatamente anterior como marco de coordenadas. La verdad actual y el VCF derivan de la secuencia final; el historial marca cada evento como activo, modificado o eliminado. Un objetivo ya eliminado detiene la configuración; no se ignora ni recoloca. Un SNP eliminado por una deleción no se emite como variante actual. Haplotipos distintos se resuelven independientemente y sus alelos finales pueden coexistir como ALT con fase y dosis. Varios eventos de un haplotipo que formen un alelo final se emiten como una sola variante o ALT compleja normalizada y usan complemento si hace falta | Casos secuenciales, multialélicos, complejos, de eliminación y de error |
| **ACE-104. Modos de frecuencia alélica** | Cada locus o bloque usa exclusivamente modo probabilístico o controlado. El probabilístico es predeterminado, muestrea desde las probabilidades y registra esperado y realizado sin tolerancia controlada. El controlado, solo explícito, elige un vector de recuentos que suma C y minimiza conjuntamente la desviación; exige 0,5/C o una tolerancia más estricta. Después distribuye copias aleatoriamente entre individuos condicionado a Hardy-Weinberg, F o la distribución explícita. En bloques LD controla recuentos haplotípicos por subpoblación. Restricciones incompatibles detienen y los empates se resuelven con la semilla y se registran | Casos por ploidía, multialelismo, bloques, modelo y restricciones |
| **ACE-105. Normalización de probabilidades y frecuencias** | Todo conjunto que deba sumar 1 se compara con una tolerancia general configurable. Una desviación dentro de ella se normaliza y registra con valores y suma antes y después; una desviación mayor o una suma no normalizable produce un error. El valor predeterminado se fija en implementación. Esta tolerancia no sustituye `0,5/C` ni otras tolerancias específicas con otra finalidad | Valor predeterminado, precisión numérica y casos dentro, fuera y no normalizables |
| **ACE-106. Catálogo declarativo común** | División, fundación, migración y selección de descendientes usan el mismo catálogo: identificadores, azar, población o grupo, alelo, genotipo, dosis o haplotipo en uno o varios loci, Y/O/NO y cantidad exacta o proporcional. Las proporciones son probabilidades individuales reproducibles y los números exactos seleccionan sin reemplazo o fallan por insuficiencia. Las reglas usan la verdad biológica. Scripts, código, funciones personalizadas y reglas desconocidas se rechazan; una regla nueva exige ampliar el programa y la especificación | Sintaxis técnica de las reglas y casos compartidos por los cuatro procesos |

Estos criterios son cualitativos. No deben interpretarse como pruebas cerradas
hasta que se resuelvan los aspectos pendientes.

## 8. Decisiones abiertas

La especificación funcional de la versión 1 no está cerrada. La revisión
[revision-adversarial-especificaciones-0.1.md](./revision-adversarial-especificaciones-0.1.md)
mantiene al final las preguntas que requieren valoración. Como mínimo, siguen
abiertas las decisiones sobre:

- el denominador y las etapas de la cobertura;
- el reparto de información entre VCF y complemento;
- el tratamiento operativo de N por cada clase de operación, las colisiones
  aleatorias y las precedencias cruzadas de configuración;
- si las inserciones de elementos transponibles pueden crear referencias
  derivadas;
- la mezcla de referencias estructurales en unidades reproductivas y la
  compatibilidad estructural en reproducción;
- la semántica verificable del tamaño poblacional N, Ne, selección, tasas
  mutacionales, F en poliploides y catálogos de modelos alternativos;
- la calibración empírica;
- la relación entre eventos complejos y registros VCF, y las reglas
  cuantitativas de validación;
- el ámbito real de las sustituciones por individuo, la cardinalidad de las
  entidades técnicas y la política de trazabilidad para escenarios grandes.

Estas decisiones cambian entradas aceptadas, resultados o errores observables y
no deben cerrarse como meros detalles de implementación. Añadir un
comportamiento o una clase de regla no confirmados requerirá revisar
expresamente esta especificación.

## 9. Decisiones técnicas aplazadas a la planificación de la implementación

Esta sección solo delimita elecciones de representación e implementación. No
resuelve ninguna de las decisiones funcionales de la sección 8. Cuando un
párrafo menciona a la vez una cuestión funcional abierta y una elección
técnica, primero debe resolverse la cuestión funcional; después podrá elegirse
la técnica sin cambiar el comportamiento aprobado.

**Formato y esquema de los ficheros de configuración y los manifiestos.** El
formato y el esquema del fichero general, de los ficheros de configuración
poblacional y de los manifiestos se elegirán durante la planificación de la
implementación teniendo en cuenta que el núcleo se desarrollará en Rust. Todos
deberán seguir siendo estructurados y legibles por el programa. Tampoco se fija
todavía una sintaxis para describir el catálogo declarativo común, los cruces
personalizados, los parentales expresos y
las condiciones de selección. Su significado, capacidad expresiva y validación
funcional ya están confirmados. La sintaxis elegida no podrá ejecutar scripts,
código ni funciones personalizadas.

**Integridad y exportación.** Se elegirán el algoritmo concreto de suma de
comprobación, los esquemas exactos de los manifiestos general y de evento y el
formato del paquete autocontenido. La elección deberá conservar sumas para todos
los archivos relevantes y permitir exportar un evento con sus referencias,
índices, configuraciones, metadatos, salidas y verificaciones necesarias.

**Índices exigidos por cada operación.** Durante esa planificación se
identificarán los tipos concretos de índice que exige cada operación o formato
utilizado y se elegirán las herramientas de indexación. Esta decisión no cambia
que cada referencia pueda tener cero, uno o varios índices según sus
necesidades.

**Generación aleatoria y reproducibilidad.** Durante la implementación y la
validación se comparará Poisson con otras distribuciones probabilísticas para
determinar cuáles reproducen mejor las propiedades esperadas en cada escenario.
Todavía deben planificarse las distribuciones candidatas, los criterios de
comparación, los conjuntos de referencia y las pruebas estadísticas; no se
presupone una distribución universalmente mejor. Poisson será el punto de
  partida cuando proceda, pero deberán evaluarse alternativas antes de fijar
  perfiles estables. También se decidirá la estrategia técnica de subsemillas o
  flujos aleatorios independientes y cómo demostrar las garantías separadas de
  REQ-CONF-111 a REQ-CONF-115 entre versiones del programa, dependencias y
  plataformas. Estas elecciones conservarán las medias esperadas configurables,
  los resultados realizados y los niveles de reproducibilidad ya exigidos.

**Tolerancias numéricas.** Durante la implementación se fijarán el valor
predeterminado y la precisión de la tolerancia general usada cuando frecuencias
o probabilidades deban sumar 1. La configuración podrá sustituir ese valor. Las
pruebas deberán cubrir normalización dentro de tolerancia, rechazo fuera de ella
y sumas no normalizables, conservando los valores antes y después. También se
fijarán las tolerancias técnicas específicas aún pendientes. Ninguna sustituirá
`0,5/C` ni otra tolerancia ya acordada con una finalidad distinta.

**Representación técnica de coordenadas.** La configuración usa coordenadas
basadas en 1 e intervalos cerrados. Durante la planificación se elegirá su
representación concreta en los ficheros estructurados y se documentarán las
conversiones a las convenciones oficiales de cada formato externo. Los VCF de
referencias distintas se construirán de forma independiente y no se definirá un
algoritmo de correspondencia entre sus registros.

**Ejecución de los procesos evolutivos.** Se elegirán los algoritmos que apliquen
la secuencia predeterminada confirmada, N objetivo y realizado, Ne objetivo y realizado,
selección, migración, eventos de
subpoblaciones, reproducción y recombinación. Se concretarán la fórmula de Ne
reproductivo, su adaptación por ploidía y las reglas cuantitativas con las que el
esquema reproductivo intentará aproximar Ne objetivo. La fórmula solo podrá usar
los progenitores con descendencia y la desigualdad o varianza de sus
contribuciones. También se concretarán la fórmula de F realizado y su adaptación
por ploidía. Esta deberá comparar la heterocigosidad observada y esperada según
las frecuencias alélicas, y conservar resultados negativos. Estas decisiones
técnicas no deben cambiar los significados funcionales ni las consecuencias de
validación acordadas.

**Estado genotípico faseado.** Se elegirán el formato del fichero de verdad
obligatorio, su esquema, la estructura interna que representa los haplotipos y
el algoritmo que mantiene fase, dosis y origen de segmentos hasta las copias
haplotípicas fundadoras. También se elegirá el formato de la traza detallada de
individuos muestreados y la representación separada de mutaciones de novo,
duplicaciones, inserciones y deleciones. La elección deberá preservar la
reconstrucción completa exigida funcionalmente sin exportar esa traza para todos
los individuos históricos.

**Representación VCF y complemento de verdad.** Se elegirán la versión oficial
exacta de VCF, su catálogo concreto de campos de cabecera, las reglas de
normalización, los campos para repeticiones en tándem, el criterio técnico de
tamaño manejable y el formato y esquema del complemento estructurado. La
validación deberá demostrar reconstrucción sin pérdida a partir de las salidas
enlazadas. El esquema debe respetar la identificación inequívoca de referencia,
conjunto y muestras, usar `GT` con fase y ploidía reales y excluir `DP`, `AD`,
`GQ`, `PL` y cualquier otro campo `FORMAT` derivado de observaciones o
inferencias.

**Repeticiones en tándem.** Se elegirán las fórmulas y distribuciones concretas
del modelo escalonado, los modelos alternativos y los modelos opcionales de
mutaciones internas. Los límites prácticos de longitud y número de loci
dependerán de los recursos, no de un límite funcional de motivo.

**Estimación de recursos.** Se elegirán y calibrarán los modelos que estimen
almacenamiento, memoria y tiempo y los umbrales que generan advertencias. Las
pruebas compararán las estimaciones con consumos reales sin convertirlas en
garantías exactas.

**Perfiles WGS, mapeo y cobertura.** Se concretarán y validarán los perfiles WGS,
parámetros y versiones de Illumina-like *single-end* y *paired-end* y PacBio
HiFi. Se elegirán las distribuciones, tasas y valores cuantitativos exactos de
longitud, error y calidad y se comprobarán con datos adecuados para cada perfil.
También se fijarán y probarán los *presets* y valores que minimap2 aplica
automáticamente y la comprobación del muestreo uniforme y la cobertura
aproximadamente Poisson. La planificación elegirá la versión concreta de
minimap2, el catálogo inicial de mapeadores alternativos y las comprobaciones
técnicas de compatibilidad de herramientas, perfiles, parámetros y referencias.
Se elegirán además las distribuciones, tasas y valores cuantitativos exactos de
duplicados, adaptadores, sesgo GC, quimeras, anomalías *paired-end* y las
categorías y MAPQ sintéticos del modo directo. Estos valores se validarán con
datos adecuados para Illumina-like y PacBio HiFi. La elección cuantitativa no
elimina el requisito de que cada degradación sea simulable, configurable y
  trazable en v1. Tampoco cambia que minimap2 es el mapeador predeterminado,
  sustituible por configuración, que cada evento que genere alineamientos usa
  exactamente un mapeador ni que herramienta, versión y parámetros efectivos
  deben quedar trazados.

**Calibración empírica WGS y verdad de lecturas.** Se elegirán el método de
estimación, las medidas de calidad y los esquemas del perfil WGS reutilizable
y de la verdad que relaciona secuencia original y observada, duplicados,
adaptadores, quimeras, posición real y estado de pareja. También se elegirá el
formato del fichero de procedencia de contaminación, sin duplicar las lecturas.

**Validación interna contra la verdad.** Se elegirán las métricas, métodos
estadísticos, ámbitos de agregación, tolerancias y esquema de los resultados
estructurados que comparan FASTQ y BAM o CRAM con la verdad de origen. También se
planificará cómo aplicar la política confirmada de errores estructurales o
lógicos, advertencias estadísticas predeterminadas y elevación configurable a
error. No se ejecutarán *SNP callers* para realizar esta validación.

**SAM/BAM, CRAM y grupos de lectura.** Se elegirán las versiones oficiales
exactas de SAM/BAM y CRAM, el catálogo concreto de campos de cabecera, el
procedimiento para marcar duplicados y la forma técnica de combinar grupos en el
único archivo final por individuo y referencia. Estas elecciones deben mantener
la identificación inequívoca de referencia, muestra, evento, biblioteca, tanda
o carril y grupo de lectura; `SM` como identificador biológico estable; y el
identificador `RG` único en toda la ejecución. Todos los archivos deberán
cumplir la versión oficial elegida.

**Validación de FASTQ externos y contaminación.** Se elegirán la sintaxis y las
unidades exactas de matrices y tasas, su representación en la configuración,
las tolerancias y validaciones técnicas, las comprobaciones de FASTQ externos y
los datos de integridad que se registrarán. También se definirán las pruebas de
concordancia funcional entre los modos realista y directo. Estas decisiones no
cambian la contaminación sustitutiva, la compatibilidad tecnológica, el uso
principal de fuentes de la misma tanda, la prioridad de una matriz dirigida
sobre la tasa total ponderada por lecturas, el efecto final del demultiplexado ni
la procedencia por identificadores y recuentos ya confirmados.

**Mutación y elementos transponibles.** Se elegirán las tasas predeterminadas,
las distribuciones cuantitativas de longitud de *indels*, el formato del espectro
SNP y el esquema del catálogo y los eventos TE. Estas decisiones no podrán
introducir dependencia de contexto SNP ni transposición espontánea en la primera
versión.

**Recombinación y reproducción poliploide.** Se elegirá el formato de los mapas,
la generación técnica de posiciones, la validación de intervalos colineales y la
traza interna de bivalentes. La implementación deberá respetar el mínimo
confirmado para cromosomas y la regla general con sustituciones por *scaffold*,
usando exclusivamente cM o bases dentro de una misma regla.

Además de confirmar Rust para el núcleo y minimap2 como mapeador predeterminado,
esta especificación no elige formatos, índices, otras herramientas, algoritmos,
*crates*, versión del lenguaje ni arquitectura. La sección 8 enumera las
familias de decisiones funcionales que siguen abiertas para el alcance de v1.

## 10. Próximos pasos de diseño

1. Resolver las decisiones funcionales abiertas de la sección 8 antes de
   planificar la fase que bloquea cada una.
2. Concretar durante la planificación las tolerancias técnicas, la sintaxis del
   catálogo declarativo común y los metadatos adicionales de muestreo. También
   se elegirán las fórmulas de F realizado y Ne reproductivo y su adaptación por
   ploidía, sin reabrir su significado funcional.
3. Definir la representación técnica de intervalos, el catálogo y los parámetros
   de distribuciones alternativas a Poisson y las pruebas y criterios
   estadísticos para compararlas antes de fijar perfiles estables.
4. Elegir durante la planificación la sintaxis del grafo de cruces, de los
   parentales expresos y de las condiciones de selección ya confirmadas.
5. Concretar y validar con datos adecuados los perfiles Illumina-like
   *single-end* y *paired-end* y PacBio HiFi, la versión, *presets* y parámetros
   de minimap2, el catálogo inicial de mapeadores compatibles, la calibración
   empírica y las tolerancias del muestreo uniforme y la cobertura
   aproximadamente Poisson.
6. Definir el algoritmo de sumas de comprobación, los esquemas de los
   manifiestos general y de evento y el formato de exportación autocontenida.
7. Definir las unidades y la representación de la contaminación, las
   validaciones de FASTQ externos y los casos de concordancia entre modos.
8. Completar los perfiles cuantitativos de errores, calidades, duplicados,
   adaptadores, sesgo GC, quimeras, anomalías *paired-end* y demás degradaciones
   WGS admitidas.
9. Definir los campos de cabecera y las versiones de SAM/BAM, CRAM y VCF que
   identificarán referencia, muestra, evento, biblioteca, tanda o carril y grupo
   de lectura.
10. Completar los casos de prueba restantes y tomar las decisiones de
    implementación cuando el comportamiento esperado esté suficientemente
    definido.
