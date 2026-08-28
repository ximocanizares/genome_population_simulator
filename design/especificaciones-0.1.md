# Especificación del simulador de genomas y poblaciones

- **Versión:** 0.2
- **Estado:** borrador F2 simplificado, pendiente de revisión
- **Fecha:** 21 de agosto de 2026
- **Fuentes disponibles:** `README.md`, información confirmada por el usuario y
  la decisión posterior de usar Rust recogida en
  [hoja-de-ruta-0.1.md](./hoja-de-ruta-0.1.md)

## 1. Estado y finalidad del documento

**Simplificación transversal de F2 (21 de agosto de 2026).** F2 adopta el
mismo modelo de confianza local mínimo que cerró F1. La persona usuaria, el
repositorio y el directorio privado de trabajo no se consideran adversarios.
Las entradas externas se validan y los fallos accidentales se diagnostican,
pero F2 no añade defensas frente a manipulación local activa, carreras hostiles,
funciones de entorno inyectadas ni sustituciones deliberadas de ficheros. Ese
endurecimiento solo se evaluará en F9 si el modo de distribución lo necesita.

La aplicación Marimo de F2 es una interfaz local para editar la configuración
que F2 realmente consume, invocar el binario Rust y mostrar diagnósticos y
resúmenes. No crea un protocolo de seguridad, un sistema de expedientes por
paso ni modelos de escenarios futuros. La vía directa usa un único TOML general;
los documentos de escenario se incorporarán en la primera fase que tenga un
consumidor biológico real. Esta decisión prevalece sobre cualquier detalle
técnico anterior que anticipase esas capacidades en F2.

Este documento define el objetivo general y los requisitos funcionales
confirmados para v1. Las elecciones técnicas pendientes se separan del contrato.

La versión 0.1 sigue siendo un borrador documental. Sus requisitos son
obligatorios y las 27 preguntas funcionales de la revisión anterior están
resueltas. La ampliación del 10 de agosto de 2026 confirma la referencia
reducida y la aplicación Marimo; la revisión del 12 de agosto resuelve matrices
y espectros SNV; la revisión del 13 de agosto resuelve mediante DTC-055 el panel
de bloques haplotípicos. Aún deben planificarse las representaciones que siguen
abiertas y convertirse los ACE en pruebas.

**Corrección de alcance F2 posterior al cierre F1.** F1 quedó cerrado en
`e14705a` con cierre documental `efca265` y no contiene vía, ploidía, reglas,
escenarios, eventos ni relaciones de resolución futuras. F2 introduce el
TOML general del recorrido directo y solo los campos que consume. Usa
internamente `scenario_id = "direct"`, literal reservado que no representa
biología ni puede escribir la persona usuaria. Se observa en configuración
efectiva, manifiesto general y contexto RNG. Los escenarios se aplazan a su
primer consumidor real. F2 prueba instalación desde repositorio y arranque
local de Marimo; `.tar.gz`, OCI y plataforma oficial se reservan a F9.
El aprovisionamiento inicial explícito mediante `uv sync --frozen` puede usar
Internet únicamente para descargar las dependencias exactas resueltas por
`uv.lock`. Una vez creado el entorno, el arranque y el uso normal de Marimo y
del simulador no consultan la red ni dependen de servicios remotos. La
instalación distribuible sin descarga previa sigue reservada a F9.

La sección 4 es la fuente normativa de requisitos funcionales. Las secciones 2,
5 y 6 los resumen y la sección 7 propone criterios de aceptación preliminares.
Si una redacción resumida difiere de un requisito identificado, prevalece el
requisito y debe corregirse el resumen. Esta regla aclara la fuente que debe
consultarse; no autoriza a ignorar contradicciones entre requisitos.

La cobertura completa REQ→ACE→fase→prueba se mantiene, sin repetir el contrato,
en [trazabilidad-requisitos-0.1.md](./trazabilidad-requisitos-0.1.md).
Toda capacidad funcional aplazada se identifica y encamina en
[registro-ampliaciones-futuras-0.1.md](./registro-ampliaciones-futuras-0.1.md);
las menciones de futuro de este documento heredan el DTF asociado a su rango de
requisitos.

El núcleo del sistema se desarrollará en Rust. Esta elección ya está confirmada
para el producto y no determina por sí sola la arquitectura, los *crates*, la
versión del lenguaje, los formatos, los esquemas ni los algoritmos.

**Convención documental de planificación.** Todo plan nuevo o ampliación de un
plan posterior al 14 de agosto de 2026 debe identificar cada paso mediante la
coordenada jerárquica completa de su hito en la hoja de ruta. Por ejemplo, los
pasos del hito F6.1 serán `F6.1.1`, `F6.1.2`, etc.; si un paso necesita una
subdivisión verificable, usará `F6.1.1.1`, sin reiniciar la numeración con
letras o números locales ambiguos. Pausas, evidencias, revisiones y matrices
deben citar el mismo identificador. Los planes F1 y F2 ya iniciados conservan
sus identificadores históricos para no romper la trazabilidad existente; esta
excepción no se extiende a planes nuevos.

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
—o una pareja de FASTQ— y un resumen técnico agregado por grupo de lectura
antes de mapearlos. La configuración podrá solicitar un VCF de verdad
sin activar FASTQ, BAM ni CRAM. El modo directo de alineamientos sintéticos
queda aplazado en DTF-051.

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
secuenciación independiente y es el único experimento de ese evento, sin una
identidad adicional. Contiene exactamente una tanda de secuenciación y una
biblioteca por individuo seleccionado. Todas las bibliotecas se procesan juntas
en esa tanda y comparten la única tecnología, el único perfil y, cuando se
generan alineamientos, el único mapeador del evento. Un evento que solo
materialice verdad genotípica o VCF no constituye un experimento de
secuenciación y no tiene tecnología, perfil, biblioteca, tanda ni mapeador. Una
misma ejecución biológica puede contener varios eventos con tecnologías,
perfiles o mapeadores distintos, pero no reutiliza FASTQ entre ellos ni tiene
como objetivo comparar mapeadores.

La primera entrega técnica F4 simulará exclusivamente secuenciación del genoma
completo (WGS, *whole-genome sequencing*) de DNA genómico con el perfil corto
Illumina *paired-end* incorporado. El modo Illumina *single-end* queda aplazado
en DTF-057. La preparación de biblioteca será exclusivamente sin PCR; la
preparación WGS amplificada y sus duplicados o sesgos quedan aplazados en
DTF-056. PacBio HiFi, Oxford Nanopore y PacBio CLR quedan para
futuras ampliaciones con el fin de reducir el alcance inicial y aplazar perfiles
de error más complejos. Esta
decisión de alcance no presupone que una tecnología moderna sea incapaz de
detectar SNP.

F4 materializa inicialmente una sola población seleccionada mediante un único
evento técnico. No acepta varios eventos técnicos ni mezcla varias poblaciones
en la misma ejecución F4. Ambas extensiones quedan aplazadas en DTF-050.

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

Un FASTA real suministrado también podrá originar una referencia reducida para
casos de prueba manejables. El sistema seleccionará regiones aleatorias
reproducibles, conservará su cromosoma de origen y orden relativo y entregará
un mapa completo entre coordenadas derivadas y originales. Esta reducción no es
captura dirigida ni una simulación WGS de representación reducida: transforma
la referencia de entrada antes de crear individuos u observaciones.

Toda capacidad configurable de v1 deberá poder autorizarse y lanzarse desde una
aplicación Marimo además de mediante TOML y el binario de línea de órdenes. La
aplicación será otra forma de crear el mismo conjunto TOML y ejecutar las mismas órdenes;
no definirá reglas biológicas distintas ni sustituirá la configuración efectiva
o los manifiestos como autoridad.

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

V1 tiene un único contrato obligatorio de trazabilidad. No ofrece niveles de
detalle configurables ni reduce automáticamente la verdad, la procedencia o los
metadatos exigidos cuando aumenta el tamaño del escenario. La escala solo puede
cambiar qué materializaciones y salidas solicita la persona usuaria antes de
ejecutar, nunca el contenido mínimo de una salida ya activada.

## 3. Fuera del alcance del proyecto

La ejecución local de v1 no intenta defenderse de una persona con acceso de
escritura simultáneo al directorio privado de trabajo o al mismo sistema de
ficheros. Las entradas proporcionadas por la persona usuaria sí se consideran
no confiables y deben validarse antes de usarse. Una vez copiadas o generadas
dentro del provisional privado, el programa protege los datos frente a errores
de entrada, fallos propios, interrupciones y corrupción accidental mediante
estados, tamaños, SHA-256 y publicación atómica. La defensa adicional frente a
manipulación local activa mediante enlaces físicos o simbólicos, sustitución de
ficheros o carreras TOCTOU pertenece al endurecimiento de despliegue y
publicación de F9; no forma parte de la reproducibilidad científica de F1.

Este proyecto se limita a generar los datos simulados y su verdad conocida. No
ejecutará programas externos de detección y genotipado de SNP, denominados *SNP
callers*. Tampoco comparará sus resultados con la verdad conocida, calculará
medidas de eficacia o fiabilidad ni evaluará herramientas externas.

La ejecución de esas herramientas, la comparación de resultados y su evaluación
pertenecerán a otro proyecto (DTF-024). Ese proyecto podrá utilizar los datos y el VCF de
verdad generados por este simulador, pero no forma parte de esta
especificación.

También quedan fuera de alcance la llamada de variantes y genotipos, la
imputación de genotipos ausentes, el faseado inferido y la corrección de
genotipos a partir de lecturas o alineamientos. Esas operaciones corresponden a
programas externos (DTF-024). La cobertura, la contaminación, los errores, la calidad y
la ambigüedad del mapeo solo afectan a FASTQ, BAM o CRAM y no modifican la verdad
genotípica simulada.

En esta primera versión quedan fuera de alcance las diferencias de ploidía del
mismo cromosoma o *scaffold* entre individuos, incluidos los modelos de
cromosomas sexuales con ploidía individual distinta y las aneuploidías
individuales. Podrán estudiarse como trabajo futuro (DTF-003).

También quedan fuera de alcance los datos crudos del instrumento de
secuenciación, sus señales, los archivos BCL y la generación detallada de
lecturas de índice i5 o i7. La contaminación técnica se representará mediante
el efecto final del demultiplexado, no mediante la simulación completa del
instrumento (DTF-025).

La primera versión no modela carriles de secuenciación ni crea un carril
implícito. Tampoco permite varias tandas dentro de un evento ni varias
bibliotecas del mismo individuo dentro de ese evento. Dividir una tanda en
carriles, repartir una biblioteca entre varios carriles o tandas y representar
réplicas de biblioteca del mismo individuo queda para futuras ampliaciones
(DTF-012).

También quedan para futuras ampliaciones las reglas específicas que apliquen
sustituciones individuales a descendientes todavía no creados y la
preasignación de sus identificadores (DTF-013). En v1, seleccionar un descendiente después
de crearlo solo lo elige para una operación posterior; no cambia los parámetros
que intervinieron en su creación.

En esta versión también quedan fuera de alcance las lecturas quiméricas formadas
por fragmentos de individuos o bibliotecas diferentes. Los intercambios de
lecturas entre muestras se representan exclusivamente mediante contaminación o
errores de demultiplexado (DTF-026).

Quedan para futuras implementaciones los paneles dirigidos, la secuenciación de
amplicones obtenidos por PCR, la genotipificación por secuenciación (GBS,
*genotyping-by-sequencing*) y los demás protocolos de captura o reducción de
representación (DTF-002). Los duplicados de PCR de una biblioteca WGS amplificada sí
pertenecen al alcance: son un artefacto de su preparación, no una secuenciación
dirigida por PCR.

También se aplazan los modelos explícitos de desequilibrio alélico de biblioteca
y pérdida alélica completa (*allele dropout*). El desequilibrio alélico es un
soporte de lecturas desproporcionado respecto a la dosis; la pérdida alélica
completa es la ausencia de lecturas de un alelo que sí existe en el genotipo.
No deben confundirse con variaciones aleatorias del soporte ni con el
desequilibrio de ligamiento entre loci (DTF-027).

La primera versión usa generaciones discretas y no solapadas. No modela edades,
supervivencia entre generaciones ni reproducción entre cohortes distintas.
Tampoco incluye selección por fecundidad, selección poligénica o epistasis,
elección de pareja por parentesco, distancia genética o genotipo, ni tasas de
mutación SNP dependientes del contexto de bases vecinas (DTF-028 y DTF-029).

En este modelo, un cambio demográfico de N no elimina individuos de una
generación que ya existe. Fija el número objetivo de descendientes que formarán
la generación siguiente. Los individuos de la generación actual solo pueden
actuar como progenitores y dejan de ser población activa cuando los sustituye su
descendencia.

La primera versión tampoco modela un componente ambiental del carácter
cuantitativo ni un parámetro de heredabilidad. La heredabilidad describe qué
parte de la variación observada de un carácter se asocia con diferencias
genéticas en una población. Estos modelos quedan para futuras ampliaciones
(DTF-030). En
esta versión, el valor del carácter usado para seleccionar coincide con el valor
genético asignado al genotipo o a la combinación de dosis.

La selección cuantitativa de la primera versión es exclusivamente direccional:
puede favorecer valores mayores o valores menores del carácter. Quedan para
futuras ampliaciones la selección estabilizadora, que favorece valores próximos
a un óptimo intermedio, y la selección disruptiva, que favorece ambos extremos
frente a los valores intermedios (DTF-006).

El único locus de selección de v1 debe ser un SNP bialélico presente exactamente
una vez en cada copia haplotípica aplicable. La selección sobre *indels*,
repeticiones, inserciones, variantes estructurales, loci multialélicos o loci con
número variable de copias queda para futuras ampliaciones (DTF-019).

La primera versión no modela elementos transponibles (TE) como clase funcional.
No admite catálogos, familias, actividad, movilización, mecanismos de copiar o
cortar y pegar, reglas de destino ni metadatos específicos de TE. Un FASTA real
se trata como una secuencia de bases sin anotación TE. Copiar una región
preexistente se representa como duplicación segmentaria genérica y añadir una
secuencia nueva, como inserción genérica. Todo modelado específico de TE queda
para futuras ampliaciones (DTF-007).

También quedan para versiones futuras las distribuciones de recuentos de
variantes aleatorias distintas de Poisson, incluidas opciones como la binomial
negativa o la Poisson inflada con ceros, y los modelos de mutación de repeticiones
en tándem que permitan saltos de varias unidades. V1 no acepta esos parámetros
ni reserva para ellos una semántica implícita (DTF-014 y DTF-015).

Los órdenes generacionales alternativos y los grafos configurables de
dependencias entre procesos quedan para versiones futuras. V1 conserva un orden
único; desactivar una etapa opcional solo la omite y no desplaza las demás
(DTF-022).

La calibración empírica automática desde datos reales queda completamente fuera
de v1. El sistema no acepta FASTQ y referencia para estimar un perfil, no
infiere tecnología ni parámetros desde esos datos y no genera ni reutiliza
ficheros de perfil empírico. Se mantienen los perfiles tecnológicos
incorporados y versionados para lecturas cortas tipo Illumina, así como las
sustituciones explícitas de parámetros ya permitidas. PacBio HiFi queda
aplazado en DTF-049. Una futura
ampliación deberá definir el contrato mínimo de entrada, el esquema del
artefacto, las métricas y umbrales de fiabilidad, la separación entre variantes,
errores, regiones ambiguas y errores de mapeo, y las reglas de reutilización
(DTF-016). La misma ampliación reserva materializar el VCF poblacional de 1000
Genomes u otras fuentes para producir perfiles empíricos exactos; F2 solo
incluye las aproximaciones compactas publicadas de DTC-068.

En reproducción sexual genérica quedan para futuras versiones las ploidías
impares, los multivalentes, la doble reducción y la preferencia de apareamiento
por subgenoma. También se aplazan la interferencia entre entrecruzamientos, la
reproducción de heterocigotos estructurales, los emparejamientos complejos, la
compatibilidad parcial, la recombinación limitada a intervalos colineales entre
estructuras distintas, los gametos desequilibrados y la fertilidad reducida.
Los modelos automáticos de F para ploidías distintas de 2, los
modelos específicos de herencia polisómica o disómica y sus estimadores de F
realizado quedan para futuras ampliaciones (DTF-004, DTF-005, DTF-020 y
DTF-021).

La primera versión solo modela entrecruzamientos como recombinación. La
conversión génica queda para futuras versiones y no debe confundirse con una
translocación no recíproca. También quedan fuera la distorsión de segregación y
el impulso meiótico (*meiotic drive*), que alteran la probabilidad de transmitir
alelos o copias y son distintos de la selección de viabilidad (DTF-005 y
DTF-031).

La reproducción de la primera versión se limita a autofecundación y cruzamiento
sexual. La reproducción clonal, la apomixis y la producción de dobles haploides
quedan para futuras ampliaciones. Tampoco se modelan cruces recíprocos, herencia
citoplasmática ni efectos maternos o paternos; los parentales son simétricos y no
se asignan sexos (DTF-032).

La primera versión tampoco reproduce entre sí individuos cuya estructura se
describa respecto a referencias estructurales distintas. Definir la referencia
de la descendencia, la identidad común de los loci y los mapas de recombinación
entre esas referencias queda para futuras ampliaciones (DTF-011). Esta exclusión no impide
mantener subpoblaciones reproductivamente aisladas basadas en referencias
distintas ni mapear un mismo individuo contra varias referencias disponibles.

La composición de variantes estructurales anidadas o de eventos que modifican
una estructura creada previamente en la misma copia haplotípica y su linaje
queda para versiones futuras (DTF-033). Esta exclusión no impide eventos no solapados o
adyacentes en un mismo haplotipo ni solapamientos situados en haplotipos
distintos.

Quedan también para futuras versiones las inserciones antes de la primera base o
después de la última; las deleciones terminales o de secuencias completas; las
inversiones terminales o de secuencias completas; las translocaciones terminales
o de cromosomas o *scaffolds* completos; y los cambios de ploidía o pérdidas
cromosómicas asociados. En una duplicación parcial dentro de la misma secuencia,
también queda fuera insertar la copia dentro de su intervalo fuente. Esta
prohibición no afecta a la excepción confirmada para una fuente que abarque la
secuencia completa (DTF-034). También se aplazan la generación de bloques de LD a partir
de estadísticos objetivo y el cálculo interno de D, D' o r² (DTF-035).

El VCF de verdad predeterminado no pretende describir posiciones invariantes ni
emitir bloques gVCF (DTF-036). La representación completa de genotipos y haplotipos
corresponde al fichero estructurado obligatorio de verdad.

Los lotes de varias ejecuciones completas quedan fuera de v1 (DTF-037). Una
futura ampliación podrá solicitar X réplicas de la misma configuración, pero
cada réplica será una ejecución normal e independiente, con referencia
aleatoria o derivada, individuos, poblaciones, eventos técnicos y salidas
propios según la vía; también tendrá ID, semilla raíz, subsemillas, directorio y
manifiesto propios. El lote necesitará identidad y resumen, reproducción completa
desde la misma raíz/configuración, estabilidad de prefijo al aumentar X y estado
por réplica, sin presentar un lote incompleto como completo. Esta capacidad solo
podrá definirse después de que existan las vías F3–F8 y deberá incorporarse a
Marimo cuando se implemente; no añade cardinalidad de referencias a F1.6 ni
controles a F2.0.

La conversión de VCF poblacionales a la matriz canónica `alt_af` queda como
utilidad auxiliar futura (DTF-038); F2 no acepta VCF como fuente de frecuencias
ni infiere estas desde genotipos, `AC/AN`, MAF o frecuencia derivada. También se
aplaza la generación aleatoria de loci SNV multialélicos (DTF-039). F2 sí admite
varios `ALT` en una matriz suministrada que supere sus comprobaciones.

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

La incorporación es incremental. F1 solo interpreta del fichero general la
semilla y la referencia inicial o reducida. No acepta todavía vía, ploidía,
reglas, individuos, poblaciones, fases, eventos ni enlaces a escenarios; si se
declaran antes de su fase propietaria, los rechaza como campos fuera de alcance
y no abre ningún fichero enlazado. F2 incorpora la composición documental de
escenarios y los hitos posteriores activan su semántica cuando exista el primer
consumidor real. Los requisitos siguientes conservan el contrato final de v1;
su asignación por fase está en la matriz de trazabilidad.

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
  solicitar cambios al azar. Los SNV heredables de una población inicial deben
  proceder de una matriz de frecuencias suministrada o generada; no se crean por
  separado para cada individuo mediante una densidad. Para *indels*, variantes
  estructurales y mutaciones *de novo* de fases posteriores, debe permitir
  elegir por clase entre una cantidad media esperada de eventos o una densidad
  por base, kilobase o megabase.
- **REQ-CONF-09.** Cuando se seleccione la vía evolutiva o la vía de familias y
  poblaciones segregantes, el fichero general debe referenciar el fichero de
  configuración poblacional correspondiente.
- **REQ-CONF-10.** Cada escenario de la vía evolutiva o de la vía de familias y
  poblaciones segregantes debe tener su propio fichero de configuración
  poblacional estructurado y legible por el programa.
- **REQ-CONF-11.** El fichero de un escenario evolutivo debe definir cómo se
  obtiene la población final desde la inicial mediante una secuencia ordenada de
  fases evolutivas. En la versión 1, el identificador de una fase es único en
  toda la ejecución, también cuando las fases proceden de escenarios distintos;
  una colisión debe rechazarse indicando ambas configuraciones.
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
  individuo declarado y existente antes de iniciar cualquier simulación. En v1,
  estos individuos son los individuos de creación directa, los fundadores y los
  parentales definidos expresamente.
- **REQ-CONF-18.** Durante la validación previa, toda sustitución por individuo
  debe referirse a un identificador existente de una de las categorías admitidas
  por REQ-CONF-17. Un identificador inexistente debe ser un error. El valor
  específico sustituye a los valores de ámbitos más generales solo para ese
  individuo; no se hereda como sustitución individual por sus descendientes. Los
  descendientes futuros deben recibir los valores que correspondan a sus ámbitos
  poblacional, temporal y espacial. Seleccionarlos después de su creación no
  debe modificar retroactivamente los parámetros con los que se generaron.
- **REQ-CONF-19.** En una población, la configuración debe poder indicar una
  frecuencia alélica como proporción poblacional. Para un individuo concreto no
  debe usar una frecuencia alélica, sino la presencia o ausencia del alelo, su
  genotipo o su dosis alélica, según la ploidía.
- **REQ-CONF-20.** La elección entre cantidad media esperada y densidad, el
  valor y, cuando proceda, la unidad de densidad deben poder definirse como
  valores globales y como sustituciones específicas por individuo, con la
  precedencia establecida en REQ-CONF-18.
- **REQ-CONF-21.** Para *indels*, variantes estructurales y mutaciones *de novo*
  por separado, la configuración debe permitir elegir una distribución espacial
  uniforme o definir regiones con una cantidad media esperada o una densidad de
  cambios mayor o menor que el valor general, incluidas regiones con densidad
  cero. En la inicialización desde una matriz de frecuencias SNV, sus filas fijan
  los loci posibles y las regiones solo pueden intervenir al generar esa matriz,
  no al volver a muestrear un individuo.
- **REQ-CONF-22.** La distribución espacial y las regiones deben poder definirse
  como valores globales y como sustituciones específicas por individuo, con la
  precedencia establecida en REQ-CONF-18.
- **REQ-CONF-23.** Toda cantidad configurada directamente debe interpretarse
  como una cantidad media esperada, no como un número exacto de eventos. Cuando
  se configure una densidad, el sistema debe obtener la cantidad media esperada
  multiplicando la densidad por la longitud válida del ámbito al que se aplica.
  Un conjunto exacto de cambios debe declararse mediante eventos explícitos; los
  recuentos alélicos exactos deben usar el modo controlado. Ninguno constituye
  otra distribución probabilística de recuentos. El número exacto de loci de una
  matriz SNV generada es la excepción explícita de REQ-CONF-129 y no fija los
  recuentos alélicos realizados de los individuos.
- **REQ-CONF-24.** El número realizado o sorteado de eventos debe obtenerse al
  azar a partir de la cantidad media esperada mediante una distribución
  probabilística y puede variar entre ejecuciones. El muestreo de alelos desde
  una matriz usa probabilidades por locus y no es un recuento de eventos.
- **REQ-CONF-25.** En v1, la distribución de Poisson debe ser la única
  distribución probabilística admitida para obtener el número realizado o
  sorteado de variantes aleatorias, por clase, región e individuo según las
  sustituciones aplicables. Una configuración que seleccione otra distribución
  debe rechazarse antes de simular e identificar el parámetro y el valor no
  admitido. Esta exclusividad no convierte en Poisson el espectro alélico ni el
  muestreo de genotipos de REQ-IND-60 a REQ-IND-65.

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
- **REQ-CONF-29.** La especificidad espacial debe ordenarse como locus puntual,
  región, cromosoma o *scaffold* y ausencia de selector espacial, pero solo para
  las familias que admitan cada nivel conforme a REQ-CONF-121. Esta jerarquía
  interna no establece prioridad frente a otro eje.
- **REQ-CONF-30.** Dos regiones solapadas para el mismo parámetro y los mismos
  selectores no espaciales deben rechazarse antes de simular, aunque sus valores
  coincidan. No debe suponerse que la región de menor longitud es más específica.
- **REQ-CONF-31.** Cuando coincidan sustituciones de ejes distintos, incluidas
  una espacial y una individual, el valor efectivo debe resolverse mediante la
  relación de dominancia de REQ-CONF-118. No debe combinarse ni priorizarse un
  eje de forma implícita.
- **REQ-CONF-32.** Cuando una población combine un objetivo de frecuencia
  alélica con presencia, ausencia, genotipos o dosis alélicas fijados para
  individuos concretos, estos valores individuales deben tratarse como
  restricciones obligatorias y tener prioridad sobre el objetivo poblacional.
- **REQ-CONF-33.** Cada fase evolutiva debe indicar, como mínimo, su número de
  generaciones, los valores o cambios de N objetivo, el modelo de apareamiento
  y los parámetros evolutivos que se usarán durante la fase. Ne realizado no es
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
  subpoblación, N objetivo debe poder variar mediante valores por generación o
  mediante fases o tramos configurados. V1 debe rechazar antes de simular todo
  Ne objetivo, absoluto o expresado como Ne/N, en cualquier ámbito. El
  diagnóstico debe identificar campo y valor y explicar que el control de Ne
  queda fuera de v1.
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
  subpoblación nueva. Debe permitir definir N objetivo, selección,
  autofecundación, apareamiento, migración y los demás parámetros aplicables a
  las subpoblaciones resultantes. Para la contribución de una fusión debe elegir
  un modo probabilístico con pesos por origen o un modo controlado con recuentos
  parentales exactos por origen. Ambos modos son excluyentes. Una recolonización debe declarar un
  identificador nuevo, distinto del de cualquier subpoblación extinguida.
- **REQ-CONF-43.** La configuración evolutiva debe permitir activar, desactivar y
  parametrizar por fase, generación y subpoblación, cuando corresponda, la selección, la
  deriva genética, la mutación y la migración. También debe poder incluir otros
  procesos evolutivos que el sistema admita.
- **REQ-CONF-44.** V1 no debe permitir reordenar las categorías de procesos de
  una generación. Debe aplicar siempre la secuencia obligatoria de REQ-REP-24.
  Cualquier campo, lista o valor que solicite un orden alternativo debe
  rechazarse antes de simular e identificar el orden solicitado y la regla fija
  de v1.
- **REQ-CONF-45.** Antes de simular, el sistema debe rechazar procesos
  desconocidos, combinaciones incompatibles y toda solicitud de orden
  generacional alternativo. El diagnóstico debe explicar los procesos, valores
  y conflicto detectados. El orden interno explícito de varios eventos
  demográficos conforme a REQ-CONF-105 no constituye una alternativa.
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
  identificar el locus y el error antes de simular. Un motivo que contenga `N`
  tampoco es un locus de repetición modelado válido y debe rechazarse antes de
  simular.
- **REQ-CONF-61.** En v1, toda mutación de una repetición en tándem debe usar el
  modelo escalonado simple: añadir exactamente una unidad completa del motivo o
  eliminar exactamente una. La tasa de mutación y el sesgo o probabilidad entre
  expansión y contracción deben seguir siendo configurables. No debe existir una
  distribución configurable del número de unidades ganadas o perdidas.
- **REQ-CONF-62.** La tasa y el sesgo entre expansión y contracción deben poder
  definirse globalmente, por locus o región y, en la vía evolutiva, por fase o
  subpoblación, aplicando REQ-CONF-116 a REQ-CONF-121. Una configuración que seleccione otro modelo o un salto de
  varias unidades debe rechazarse antes de simular e identificar el parámetro y
  el valor no admitido.
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
- **REQ-CONF-66.** F4 admite exclusivamente el perfil predefinido Illumina
  *paired-end* y sus sustituciones explícitamente permitidas. Una solicitud
  Illumina *single-end* se rechaza e indica DTF-057; PacBio HiFi queda aplazado
  en DTF-049.
- **REQ-CONF-67.** Cuando la salida BAM o CRAM esté activa, la configuración debe
  usar el modo realista: simula y conserva los FASTQ finales antes de mapearlos.
  El modo directo de alineamientos sintéticos queda aplazado en DTF-051.
- **REQ-CONF-68.** En el modo realista, los FASTQ deben activarse obligatoriamente
  junto con BAM o CRAM. Una configuración que solicite alineamientos sin FASTQ
  final debe rechazarse antes de generar datos.
- **REQ-CONF-69.** Un evento con alineamientos debe seleccionar exactamente un
  mapeador, es decir, un programa que alinea lecturas contra una referencia. El
  predeterminado debe ser minimap2. La configuración puede elegir otro mapeador
  compatible y sus parámetros, pero no varios dentro del mismo evento.
- **REQ-CONF-70.** Cuando el mapeador sea minimap2, el perfil de tecnología debe
  seleccionar automáticamente el *preset* y los parámetros predeterminados
  apropiados para lecturas cortas tipo Illumina.
  La configuración debe permitir sustituirlos y validar, cuando sea posible,
  las opciones desconocidas o incompatibles, además de advertir de cada
  sustitución.
- **REQ-CONF-71.** Si se seleccionan varias referencias de mapeo dentro de un
  evento, el único mapeador seleccionado debe aplicarse a todas. La
  primera versión no debe reutilizar FASTQ entre eventos para comparar
  mapeadores; cada evento debe generar sus propias lecturas.
- **REQ-CONF-72.** En F4, la cobertura se configura como una única cobertura
  media objetivo haploide en «x» para todo el evento. Debe admitir valores altos,
  bajos o cero. El denominador de cada biblioteca es la longitud haploide de su
  referencia estructural definida en REQ-WGS-04.
- **REQ-CONF-73.** F4 no admite sustituciones de cobertura por individuo,
  cromosoma, *scaffold* o región; debe rechazarlas antes de generar datos.
- **REQ-CONF-74.** Las regiones de cobertura y sus coordenadas no forman parte
  de F4. Una solicitud de cobertura regional debe rechazarse, sin convertirla en
  una aproximación global.
- **REQ-CONF-75.** La configuración debe expresar un rango inclusivo de
  porcentaje total de contaminación biológica y puede restringir los receptores
  a un subconjunto explícito. Si se omite el subconjunto, todos los individuos de
  la población seleccionada son receptores. Cada receptor sortea
  reproduciblemente su porcentaje total dentro del rango; no se declaran
  porcentajes ni variaciones independientes por donante. El porcentaje sorteado
  se aplica al número de parejas limpias de la receptora y se redondea a la
  pareja entera más próxima, con empate hacia arriba.
- **REQ-CONF-76.** Cuando una fuente contaminante se cuantifique por número en
  datos *paired-end*, el valor debe contar pares completos, no lecturas sueltas.
  La configuración de esta operación no debe permitir solicitar pares
  huérfanos.
- **REQ-CONF-77.** En F4, las fuentes contaminantes son exclusivamente otros
  individuos simulados de la misma población que la receptora. No se admiten
  FASTQ externos, otra población ni una fuente declarada manualmente. Todas las
  fuentes elegibles comparten tecnología, perfil y disposición con la receptora.
- **REQ-CONF-78.** Cada evento que genere FASTQ, BAM o CRAM debe declarar
  exactamente una tanda de secuenciación. La contaminación técnica debe
  configurarse dentro de esa tanda, identificando todas sus muestras,
  bibliotecas e índices. Debe distinguirse expresamente de la contaminación
  biológica definida por fuentes de lecturas.
- **REQ-CONF-79.** F4 exige para cada biblioteca un índice dual único (UDI,
  *Unique Dual Index*) explícito: `i7` e `i5` son cada uno cadenas de exactamente
  ocho bases del alfabeto `ACGT`, en orientación de secuenciación. Cada valor
  `i7` y cada valor `i5` es único dentro de la tanda. La cabecera los escribe
  literalmente como `i7+i5`, sin complementar ni invertir `i5`. F4 rechaza
  índices simples, duales combinatorios, longitudes o alfabetos distintos,
  duplicados y la ausencia de UDI; el marcador `0` queda reservado para perfiles
  futuros que no usan índice y no es válido en F4.
- **REQ-CONF-80.** En F4, la tasa de *index hopping* es la fracción configurable
  de parejas de la tanda que se reasigna reproduciblemente a otra biblioteca,
  elegida uniformemente entre las restantes. Una tasa positiva exige al menos
  dos bibliotecas; de otro modo se rechaza antes de generar. No hay una
  procedencia ni pesos de destino configurables. Cada pareja se selecciona de
  forma independiente mediante un sorteo Bernoulli de esa tasa y el recuento
  realizado se informa, sin imponer el valor esperado. El sorteo se aplica una
  vez a la instantánea completa posterior a la contaminación biológica, con la
  biblioteca pre-salto congelada; todas las transferencias se materializan
  simultáneamente. Una pareja recibida no vuelve a ser candidata y el resultado
  no depende del orden de las bibliotecas.
- **REQ-CONF-81.** F4 no genera combinaciones UDI inesperadas: el *index
  hopping* solo transfiere una pareja al UDI completo ya declarado para otra
  biblioteca de la tanda. Si una combinación inesperada alcanzara la frontera de
  demultiplexado, se descarta defensivamente y se registra solo su recuento
  agregado. La generación de índices i5/i7 y de combinaciones inesperadas queda
  aplazada en DTF-025. El perfil inicial usa UDI y fija la tasa de reasignación a
  cero.
- **REQ-CONF-82.** Cada tecnología de secuenciación debe tener un perfil
  predeterminado identificado y versionado que incluya, como mínimo, la tasa y
  los tipos de error habituales que pretende representar y su longitud máxima
  respaldada. No debe usarse un perfil «estándar» implícito o sin versión. El
  perfil inicial HiSeq X admite una única longitud común a R1 y R2 entre 1 y
  150 pb: una menor recorta sus curvas y una mayor se rechaza hasta seleccionar
  un perfil compatible. Longitudes asimétricas entre mates quedan fuera de F4.
  En F4, el inserto mínimo debe ser igual o mayor que esa longitud común; de
  otro modo la configuración se rechaza porque el perfil no modela adaptadores.
- **REQ-CONF-83.** En F4, el perfil de errores debe aplicarse homogéneamente a
  toda la población seleccionada del evento. Los valores pueden diferir entre
  R1 y R2, pero no por individuo, cromosoma, *scaffold* o región.
- **REQ-CONF-84.** En F4, la configuración puede sustituir globalmente por
  evento el perfil de errores y, para cada mate, su desplazamiento de calidad y
  escala de error. Si se omiten, usa explícitamente `0 Phred` y `1`,
  respectivamente. Debe rechazar selectores por individuo, secuencia o región;
  dichas sustituciones quedan fuera del alcance de F4.
- **REQ-CONF-85.** Antes de iniciar cualquier simulación, la configuración debe
  rechazar todo FASTQ, referencia u otro campo declarado como entrada para
  calibración empírica automática. El diagnóstico debe identificar el campo y
  explicar que esa capacidad queda fuera de v1.
- **REQ-CONF-86.** Antes de iniciar cualquier simulación, la configuración debe
  rechazar toda solicitud de estimar, inferir o generar un perfil empírico,
  incluida la inferencia de tecnología o parámetros desde datos reales. El
  diagnóstico debe identificar el parámetro y explicar que queda fuera de v1.
- **REQ-CONF-87.** Antes de iniciar cualquier simulación, la configuración debe
  rechazar la selección o reutilización de un fichero de perfil empírico como
  perfil efectivo. Este rechazo no debe impedir seleccionar los perfiles
  tecnológicos incorporados para Illumina ni aplicar las
  sustituciones explícitas permitidas a sus parámetros.
- **REQ-CONF-88.** F4 admite exclusivamente preparación WGS sin PCR. Rechaza un
  perfil amplificado y cualquier tasa de duplicados PCR u ópticos; el modelo
  correspondiente queda aplazado en DTF-056.
- **REQ-CONF-89.** El perfil inicial declara explícitamente duplicados PCR no
  aplicables y ópticos a cero. F4 no permite sustituir esos valores.
- **REQ-CONF-90.** El marcado de duplicados en BAM o CRAM queda aplazado en
  DTF-055 y DTF-056; F4 no ofrece esta elección.
- **REQ-CONF-91.** La configuración debe permitir activar o desactivar la
  incorporación de adaptador cuando una lectura sobrepase el inserto, solo en un
  perfil que admita fragmentos más cortos que la lectura. El perfil inicial
  HiSeq X usa insertos de 250–450 pb y lecturas de hasta 150 pb, por lo que no
  produce adaptadores; activar esa capacidad en F4 se rechaza y remite a DTF-054.
- **REQ-CONF-92.** El perfil inicial no admite secuencia, frecuencia, proporción
  ni longitud residual de adaptador. El contrato de un perfil futuro de
  fragmentos cortos corresponde a DTF-054.
- **REQ-CONF-93.** Un perfil futuro puede aportar secuencia o perfil de adaptador
  y solo generarlo por sobrepasar el inserto, nunca como degradación aleatoria
  implícita.
- **REQ-CONF-94.** F4 usa una relación GC neutra y rechaza curvas o intensidades
  no nulas. Un perfil con sesgo GC calibrado queda aplazado en DTF-056.
- **REQ-CONF-95.** En F4 no hay combinación de sesgo GC ni excepciones
  regionales: el muestreo WGS básico conserva su cobertura objetivo común.
- **REQ-CONF-96.** La configuración debe permitir definir o generar por regiones
  secuencias repetitivas, duplicadas o altamente homólogas destinadas a producir
  baja mapeabilidad. También debe poder aprovechar las repeticiones y variantes
  estructurales ya presentes en los genomas simulados.
- **REQ-CONF-97.** La configuración debe rechazar parámetros exclusivos de
  alineamientos sintéticos directos e indicar DTF-051 como destino.
- **REQ-CONF-98.** La configuración debe rechazar antes de generar datos toda
  solicitud de alineamientos directos e identificar el parámetro incompatible.
- **REQ-CONF-99.** La configuración debe permitir generar lecturas quiméricas,
  formadas por fragmentos que no eran contiguos en el genoma o haplotipo de
  origen, con una frecuencia configurable. Los perfiles tecnológicos que
  contemplen este mecanismo deben aportar valores predeterminados. El perfil
  inicial fija frecuencia cero. En F4, esa frecuencia es la fracción de parejas
  limpias que se transforman en quiméricas antes de contaminación y
  demultiplexado: no añade ni retira parejas. Cada quimera une exactamente dos
  segmentos, con punto de unión uniforme entre las bases internas del inserto;
  las fuentes se eligen entre posiciones válidas no contiguas. La primera fuente
  es el prefijo de la plantilla original; para un inserto de longitud `L`, el
  punto `b` se sortea uniformemente en `1..L-1` y la segunda fuente aporta
  exactamente `L-b` bases. Sus candidatas son todos los intervalos válidos de
  esa longitud en las copias del mismo individuo y biblioteca, ordenados por
  copia, secuencia y inicio; se excluye una candidata de la misma copia que
  solape o toque el intervalo del prefijo. Se elige uniformemente una candidata
  restante. Las dos fuentes se leen en orientación directa de su copia y la
  disposición PE normal se aplica después al inserto compuesto. Pueden proceder
  de cromosomas, *scaffolds*, haplotipos o copias distintos, pero nunca de
  individuos o bibliotecas diferentes. Si una pareja seleccionada no tiene
  candidata válida, la ejecución falla antes de publicar salidas. Cada pareja
  limpia se selecciona de forma independiente con un sorteo Bernoulli de esa
  tasa; el recuento realizado se informa y no se fuerza al valor esperado.
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
  exactamente una tecnología, un perfil tecnológico y una tanda. Debe declarar
  exactamente una biblioteca por cada individuo seleccionado y tratar juntas
  todas las bibliotecas dentro de la única tanda. Cuando genere BAM o CRAM debe
  seleccionar además exactamente un mapeador conforme al modo configurado. Un
  evento que solo solicite verdad genotípica o VCF no constituye un experimento
  de secuenciación y no debe declarar tecnología, perfil, bibliotecas, tanda ni
  mapeador. Eventos distintos pueden usar valores distintos sin compartir
  lecturas ni resultados técnicos.
- **REQ-CONF-103.** La configuración debe permitir solicitar de forma
  independiente, para cada evento, su VCF de verdad, sus FASTQ y sus BAM o CRAM,
  respetando las dependencias entre el modo realista y las salidas solicitadas.
- **REQ-CONF-104.** La configuración evolutiva debe indicar si permite que una
  subpoblación se extinga cuando un resultado de la simulación produce N
  realizado igual a 0. El valor global predeterminado debe ser desactivado y
  debe poder sustituirse por subpoblación y fase evolutiva conforme a la
  dominancia de REQ-CONF-118. Esta opción es distinta de programar un evento
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

- **REQ-CONF-116.** Cada parámetro debe tener un valor general o predeterminado,
  o una regla base sin selectores, salvo que su familia exija expresamente otra
  forma. Una sustitución debe contener uno o más selectores admitidos para su
  familia y solo debe aplicarse cuando todos coincidan con el objeto y el estado
  evaluados.
- **REQ-CONF-117.** La especificidad dentro de cada eje debe ordenarse así:
  generación sobre fase sobre ausencia; subpoblación sobre población sobre
  ausencia; locus sobre región sobre secuencia sobre ausencia; evento técnico
  concreto sobre ausencia; e individuo concreto sobre ausencia. En un eje
  nominal ya admitido, como grupo o bloque, una coincidencia concreta debe ser
  más específica que su ausencia. No debe existir prioridad implícita entre
  ejes distintos. Esta ordenación solo se aplica después de que la fase
  propietaria haya validado el contexto con sus datos reales. F2 aportará loci
  y posiciones; F6, población y subpoblación; y F7, fases e intervalos de
  generaciones. F1 solo resolverá los valores y selectores que consuman
  realmente sus recorridos de configuración y referencias; no construirá
  dominios sintéticos ni inferirá relaciones de fases posteriores desde listas
  planas de identificadores.
- **REQ-CONF-118.** Una regla A debe dominar a otra B cuando sea igual o más
  específica en todos los ejes admitidos y aplicables y estrictamente más
  específica en al menos uno. Tras reunir las reglas aplicables se deben retirar
  las dominadas. Una única regla máxima determina el valor efectivo. Varias
  máximas incomparables deben causar un rechazo previo aunque sus valores
  coincidan; el diagnóstico debe identificar parámetro, objeto o ámbito y todas
  las reglas en conflicto. El conflicto solo se resuelve eliminando una regla o
  añadiendo una regla combinada que coincida en los ejes necesarios y domine a
  las anteriores.
- **REQ-CONF-119.** Dos reglas del mismo parámetro con exactamente los mismos
  selectores deben rechazarse antes de simular aunque expresen el mismo valor.
- **REQ-CONF-120.** Un selector no admitido por la familia del parámetro debe
  rechazarse antes de simular e identificar el campo y el ámbito no permitido.
  Una regla base continúa hasta que otra regla aplicable la sustituya; una regla
  de fase solo actúa dentro de esa fase y una de generación solo en esa
  generación. Una sustitución individual se limita a los individuos
  preexistentes de REQ-CONF-17 y REQ-CONF-18: no se hereda a descendientes
  futuros ni modifica retroactivamente su creación.
- **REQ-CONF-121.** Las familias de parámetros deben admitir únicamente los ejes
  de esta matriz normativa; cada celda se limita a los parámetros para los que
  los requisitos de esa familia ya declaran una sustitución:

  | Familia | Ejes de sustitución admitidos en v1 | Herencia o vigencia |
  | --- | --- | --- |
  | Estructura y ploidía | secuencia | toda la ejecución |
  | Creación, frecuencias y mutación | individuo preexistente; población o subpoblación; fase o generación; secuencia, región o locus; grupo o bloque cuando esté definido | base hasta sustitución; fase y generación en su intervalo; individuo sin herencia a descendientes |
  | Demografía, migración, reproducción y selección | población o subpoblación; fase o generación; locus solo para condiciones genéticas ya admitidas; grupo cuando esté definido | base hasta sustitución; fase y generación en su intervalo |
  | Recombinación | grupo y secuencia | toda la ejecución o el grupo declarado |
  | WGS, cobertura, perfiles y errores | evento técnico, individuo, secuencia y región | dentro del evento técnico; una regla sin evento conserva su vigencia declarada |
  | Contaminación y demultiplexado | evento técnico e individuo o muestra representada por su biblioteca | solo dentro de la tanda del evento técnico |
  | Salidas y eventos de observación | evento técnico; individuo como miembro seleccionado, no como modificación biológica | solo dentro del evento técnico |

- **REQ-CONF-122.** Los parámetros biológicos no deben admitir el evento
  técnico como selector. Un evento técnico materializa la verdad ya existente y
  ninguna sustitución técnica puede cambiar genotipo, genealogía, trayectoria
  biológica ni verdad. Los eventos demográficos o reproductivos conservan sus
  campos y ejes biológicos y no deben interpretarse como eventos técnicos.
- **REQ-CONF-123.** La configuración efectiva debe registrar, para cada
  parámetro y objeto evaluado, el valor efectivo, la regla ganadora y la
  procedencia y selectores de esa regla. Ante un rechazo por ambigüedad, el
  diagnóstico debe conservar todas las reglas máximas en conflicto.
  Cuando una fase resuelva una relación jerárquica real, debe conservar además
  la versión y procedencia de los datos que acreditan esa relación. La fase
  propietaria define y valida ese contrato junto con su primer consumidor; no
  se exige a F1 publicar una evidencia sintética independiente ni anticipar el
  orden de parámetros que aún no consume.

- **REQ-CONF-124.** La configuración general debe permitir declarar una
  referencia reducida con identificador propio a partir de un FASTA real
  suministrado, identificando el FASTA de origen y un `region_size_bases`
  entero positivo. La unidad normativa es la base; una interfaz puede mostrar
  el valor convertido a kb sin cambiar ni redondear el valor almacenado.
- **REQ-CONF-125.** Cada referencia reducida debe declarar exactamente uno de
  estos dos objetivos mutuamente excluyentes: `target_bases`, entero positivo,
  o `reduction_percent`, porcentaje finito en el intervalo `(0,100)`. Declarar
  ambos, ninguno o un valor fuera de su dominio debe rechazarse antes de copiar
  regiones. El objetivo retenido debe ser menor que el tamaño elegible; para
  conservar la referencia completa se usa directamente la referencia
  suministrada, no una referencia reducida idéntica.
- **REQ-CONF-126.** La selección de regiones de cada referencia reducida debe
  usar la semilla raíz de la ejecución y una subsemilla propia derivada del
  identificador estable de la referencia. No debe consumir un generador global
  mutable ni cambiar al añadir individuos, eventos técnicos o salidas.
- **REQ-CONF-127.** La configuración efectiva debe registrar el FASTA e ID de
  origen, el ID derivado, `region_size_bases`, la modalidad y valor del objetivo,
  el número exacto de bases derivadas, la semilla y subsemilla aplicables y la
  versión de la regla de selección. Debe conservar los valores declarados sin
  sustituirlos silenciosamente por el resultado realizado.
- **REQ-CONF-128.** La creación aleatoria de individuos desde una población de
  origen debe aceptar exactamente una matriz de frecuencias SNV canónica ya
  preparada o solicitar que el simulador la genere. La matriz debe pertenecer a
  una única referencia estructural y contener una o varias poblaciones. Cada
  individuo debe seleccionar exactamente una de ellas; sus IDs deben coincidir
  literalmente y no se admite mezcla de ancestrías en F2. Varias referencias
  exigen matrices independientes y una población no seleccionada puede
  conservarse como disponible sin intervenir en los cálculos realizados.
- **REQ-CONF-129.** La generación de una matriz de frecuencias SNV debe poder
  configurar exactamente uno de: número exacto de loci candidatos, media
  Poisson o densidad por Mb convertida a media sobre las bases elegibles; además
  debe configurar la distribución espacial, el espectro de sustituciones y un
  espectro de frecuencias `alt_af`. El número exacto se conserva como candidatos
  y la media/densidad registra media y sorteo; un candidato descartado después
  no se reemplaza. En el modo `exact`, las cuotas por fila se calculan primero
  mediante mayores restos. Después, cada población permuta reproduciblemente
  todos los `variant_id` canónicos mediante una subsemilla exclusiva de
  asignación y aplica las cuotas siguiendo esa permutación; está prohibido
  asociar las filas por coordenada u orden canónico de variante. El valor dentro
  de un intervalo se obtiene con otro flujo separado por población y
  `variant_id`; una fila `fixed` asignada no consume RNG para el valor y una fila
  con cuota cero no consume RNG. V1 debe incluir el ejemplo explícito y versionado
  `human_1000g_global_sfs_approx_v1`, aproximación GLOBAL a los recuentos
  publicados finales de 1000 Genomes, con fuente, receta, transformación y
  SHA-256. No constituye perfiles exactos por superpoblación. Es solo un ejemplo
  seleccionable: nunca se carga por defecto. La persona usuaria debe elegir
  expresamente ese perfil o suministrar un perfil válido, que puede ser distinto
  para cada población.
- **REQ-CONF-130.** Al convertir una matriz de frecuencias en genotipos, cada
  población de origen debe elegir entre Hardy-Weinberg, que es el modo
  predeterminado, y un coeficiente de endogamia `F`. El coeficiente `F` solo se
  admite en diploides y debe estar entre 0 y 1, inclusive; en haploides y
  poliploides se mantiene Hardy-Weinberg o una distribución genotípica explícita
  ya admitida, sin inferir una generalización de `F`.

- **REQ-TEC-01.** Cada evento que genere FASTQ, BAM o CRAM constituye también su
  único experimento de secuenciación, sin una identidad de experimento separada,
  y debe seleccionar exactamente una tecnología. No debe mezclar lecturas
  cortas tipo Illumina de perfiles distintos dentro del mismo evento. PacBio
  HiFi queda aplazado en DTF-049. Un evento que solo produzca verdad genotípica o VCF no constituye un
  experimento de secuenciación.
- **REQ-TEC-02.** Cada experimento debe contener exactamente una tanda y una o
  varias bibliotecas: exactamente una por cada individuo seleccionado. Todas las
  bibliotecas deben tratarse conjuntamente en la tanda y ser compatibles con la
  tecnología, el perfil de lecturas y el mapeo seleccionados para el evento. La
  primera versión no debe admitir carriles.
- **REQ-TEC-03.** La configuración debe permitir crear experimentos separados
  por tecnología a partir de la misma simulación biológica, incluidos eventos de
  muestreo distintos, con trazabilidad entre ellos y sin mezclar sus archivos ni
  reutilizar sus lecturas. Otra biblioteca del mismo individuo, otra tanda, otro
  perfil o tecnología y todo control independiente deben declararse mediante
  otro evento.
- **REQ-TEC-04.** El mapeador, el *preset* y sus parámetros deben ser
  compatibles con la única tecnología del experimento y aplicarse de forma
  coherente a todas sus bibliotecas, a su única tanda y a todas sus referencias.
- **REQ-TEC-05.** Dos perfiles que requieran modelos de error, calidades o
  parámetros de mapeo distintos deben considerarse tecnologías experimentales
  diferentes aunque pertenezcan al mismo fabricante. Cada tecnología admitida
  debe tener un identificador y una versión explícitos en la configuración.
- **REQ-TEC-06.** Todo experimento de la primera versión debe representar WGS de
  DNA genómico. La configuración debe rechazar paneles dirigidos, amplicones,
  GBS y otros protocolos de captura o reducción de representación.
- **REQ-TEC-07.** F4 usa exclusivamente preparación WGS sin PCR y la registra
  como tal. Debe rechazar una preparación con amplificación e indicar DTF-056.
- **REQ-TEC-08.** El perfil WGS amplificado, sus duplicados PCR y sus sesgos
  propios quedan aplazados en DTF-056; F4 no acepta sus parámetros ni los
  representa como secuenciación dirigida por PCR.
- **REQ-TEC-09.** F4 solo debe admitir lecturas cortas Illumina *paired-end*.
  Debe rechazar Illumina *single-end* e indicar DTF-057, y rechazar PacBio HiFi,
  Oxford Nanopore y PacBio CLR como tecnologías de simulación con destino
  DTF-049 y DTF-001.
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
- **REQ-REF-05.** En v1, una referencia debe tener cero índices durante una
  ejecución normal y puede tener un único índice FAI cuando una operación lo
  necesite. Ese FAI debe crearse bajo demanda y no debe publicarse como salida
  normal. La coexistencia de varios índices de referencia queda aplazada a
  DTF-048.
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
  único símbolo admitido para una base desconocida o no determinada. Debe
  conservarse exactamente, ocupar su coordenada y contar en la longitud de la
  secuencia; nunca debe resolverse al azar ni sustituirse silenciosamente por
  A, C, G o T.
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
  diagnóstico. Las bases `N` de referencias, haplotipos, secuencias explícitas
  o copias deben conservarse con su posición. Una referencia nueva generada
  aleatoriamente puede declarar una fracción independiente de posiciones `N`;
  entre las posiciones restantes, las frecuencias A, C, G y T deben sumar 1.
- **REQ-REF-24.** El sistema debe poder crear una referencia reducida con ID
  propio a partir de las regiones de los cromosomas o *scaffolds* de un FASTA real
  suministrado. La referencia reducida es una referencia derivada de
  preparación y no debe modificar ni reemplazar la copia autoritativa del FASTA
  de origen.
- **REQ-REF-25.** El FASTA de origen de una reducción debe superar exactamente
  el mismo preflight, copia byte a byte, validación de cabeceras y alfabeto y
  cálculo de SHA-256 que cualquier FASTA suministrado en REQ-REF-12 a
  REQ-REF-18. La selección solo puede leer esa copia validada.
- **REQ-REF-26.** La selección debe ser aleatoria y reproducible, sin reemplazo
  y sin solapamiento entre intervalos de origen. No debe duplicar regiones,
  invertirlas, cambiar su orientación ni crear una correspondencia que no
  exista en el FASTA de origen.
- **REQ-REF-27.** La referencia reducida debe contener exactamente el número de
  bases objetivo. Para `target_bases`, ese número es el entero declarado. Para
  `reduction_percent=p`, `p` expresa el porcentaje eliminado y el tamaño
  retenido debe ser `floor(S × (100 − p) / 100)`, donde `S` es el tamaño total
  elegible. El cálculo debe usar el valor decimal declarado de forma exacta,
  redondear una sola vez hacia abajo y registrarse antes de seleccionar
  regiones. Un resultado menor que una base debe rechazarse.
- **REQ-REF-28.** La selección debe usar intervalos de longitud
  `region_size_bases`. Si el objetivo no es múltiplo de esa longitud, debe usar
  exactamente un intervalo residual más corto para alcanzar el tamaño exacto;
  no debe redondear el objetivo, completar con `N`, repetir bases ni producir
  más de un residual.
- **REQ-REF-29.** Antes de materializar, los intervalos elegidos deben ordenarse
  por el orden de cromosomas o *scaffolds* del FASTA de origen y por coordenada
  creciente dentro de cada secuencia. La secuencia derivada correspondiente debe
  concatenar sus intervalos en ese orden y conservar siempre orientación
  directa.
- **REQ-REF-30.** Un cromosoma o *scaffold* de origen sin intervalos elegidos debe omitirse de
  la referencia reducida. El sistema no debe crear secuencias FASTA vacías y la
  referencia reducida debe contener al menos una secuencia y una base.
- **REQ-REF-31.** Cada secuencia de la referencia reducida debe conservar de
  forma inequívoca la identidad y clase del cromosoma o *scaffold* fuente. La
  secuencia derivada debe reutilizar exactamente su `sequence_id` y su clase;
  el ID de referencia derivada, distinto del ID de origen, separa ambos espacios
  de coordenadas. No se mezclan regiones de secuencias fuente diferentes.
- **REQ-REF-32.** Toda `N` incluida en un intervalo seleccionado debe copiarse
  literalmente, conservar su posición relativa dentro del intervalo y contar
  como una base para el tamaño objetivo. La reducción nunca debe resolver,
  excluir o reemplazar `N` por su contenido.
- **REQ-REF-33.** Un objetivo igual o mayor que el total elegible, un objetivo calculado
  menor que una base, la falta de espacio para completar intervalos no
  solapados o cualquier combinación incompatible debe rechazarse antes de
  materializar la referencia derivada. El rechazo no puede reducir el objetivo
  ni aceptar una referencia parcial.
- **REQ-REF-34.** El resultado válido debe ser una nueva referencia con ID
  propio, un único FASTA canónico, tamaño exacto y SHA-256 sobre sus bytes. El
  FASTA de origen, su copia autoritativa y la referencia reducida deben
  conservar identidades propias y registrar por separado tamaño y suma. La
  igualdad de una suma solo es válida cuando los bytes comparados son realmente
  idénticos y nunca fusiona las identidades de origen y derivada.
- **REQ-REF-35.** La reducción debe producir un mapa estructurado completo con,
  al menos, referencia y cromosoma de origen, intervalo de origen, referencia y
  cromosoma derivados, intervalo derivado, orientación directa y orden de cada
  fragmento. Los intervalos usan coordenadas basadas en 1 y cerradas en el
  contrato de configuración y salida.
- **REQ-REF-36.** El FASTA reducido debe poder reconstruirse exactamente desde
  la copia autoritativa del FASTA de origen, el mapa, los parámetros efectivos,
  la semilla y la versión del algoritmo de selección. El mapa y el manifiesto
  no pueden depender de recordar el orden en que el generador sorteó candidatos.

**Decisión funcional cerrada el 14 de agosto de 2026.** `reduction_percent`
expresa el porcentaje eliminado, usa el cálculo decimal exacto y el redondeo de
REQ-REF-27 y no admite cero. Cromosomas y *scaffolds* forman el espacio elegible
y cada secuencia derivada conserva su ID y clase sin mezclarse con otra.

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

El modelo escalonado simple de mutación cambia una repetición en tándem mediante
la ganancia o pérdida de exactamente una unidad completa del motivo por evento.
Es el único modelo de v1. Puede configurarse la probabilidad relativa de
expansión frente a contracción, pero no el tamaño del salto.

- **REQ-IND-01.** El sistema debe crear individuos que contengan SNP, *indels*,
  variantes estructurales y repeticiones en tándem respecto a una referencia.
- **REQ-IND-02.** La estructura de cada individuo debe poder corresponder a
  cualquiera de las referencias disponibles.
- **REQ-IND-03.** El sistema debe conservar qué referencia disponible describe
  la estructura de cada muestra.
- **REQ-IND-04.** Los cambios deben poder crearse al azar. Los SNV heredables de
  la población inicial usan la matriz de frecuencias de REQ-IND-60 a REQ-IND-65.
  Las demás clases y las mutaciones *de novo* posteriores pueden configurar, por
  clase, una cantidad media esperada de eventos o una densidad por base,
  kilobase o megabase.
- **REQ-IND-05.** Los cambios también deben poder definirse de forma concreta en
  el fichero general de configuración.
- **REQ-IND-06.** El sistema debe poder simular ploidías enteras de 1 a 10
  copias cromosómicas, inclusive.
- **REQ-IND-07.** En una misma ejecución, todos los individuos deben tener la
  misma ploidía efectiva para un cromosoma o *scaffold* dado.
- **REQ-IND-08.** En esta primera aproximación, no se deben permitir cambios de
  ploidía entre individuos de una misma ejecución.
- **REQ-IND-09.** La cantidad media esperada o densidad de *indels*, variantes
  estructurales y mutaciones *de novo* debe poder tener un valor global y valores
  específicos para los individuos que se indiquen. Los SNV iniciales usan una
  frecuencia por población y solo admiten restricciones genotípicas por individuo.
- **REQ-IND-10.** Un valor específico de cantidad media esperada o densidad debe
  sustituir al valor global únicamente para el individuo correspondiente.
- **REQ-IND-11.** La frecuencia alélica debe representar la proporción de copias
  de un alelo en una población y no debe atribuirse literalmente a un solo
  individuo.
- **REQ-IND-12.** Para un individuo concreto, el sistema debe expresar el estado
  de un alelo mediante su presencia o ausencia, su genotipo o su dosis alélica,
  según la ploidía del individuo.
- **REQ-IND-13.** Para cada clase generada por cantidad o densidad, el sistema debe poder
  distribuir espacialmente los eventos de manera uniforme o aplicar regiones
  configuradas con una cantidad media esperada o densidad mayor o menor que el
  valor general.
- **REQ-IND-14.** El sistema debe admitir regiones configuradas con densidad
  cero, en las que no se generen al azar eventos de la clase de variante
  correspondiente.
- **REQ-IND-15.** Una distribución o región específica de un individuo debe
  sustituir a la configuración global correspondiente solo para ese individuo.
  Los demás individuos deben conservar la configuración global.
- **REQ-IND-16.** Para cada clase generada por cantidad o densidad, ámbito general o región e
  individuo al que se aplique la generación aleatoria, el sistema debe conservar
  la cantidad media esperada, el número realizado o sorteado, el número aceptado
  antes de un posible fallo y el número finalmente activo.
- **REQ-IND-17.** El número realizado debe ser entero y puede diferir de la
  cantidad media esperada y entre ejecuciones que no estén obligadas a repetir
  las mismas decisiones aleatorias.
- **REQ-IND-18.** El número realizado de variantes aleatorias generadas desde una
  cantidad media debe obtenerse exclusivamente con la distribución de Poisson.
  Esta regla no se aplica al muestreo de alelos desde una matriz de frecuencias.
  La misma semilla, configuración, entradas y versión deben reproducir el mismo
  recuento.
- **REQ-IND-19.** Las coordenadas de una región de generación deben interpretarse
  sobre la referencia estructural que corresponde al individuo o a la población
  a los que se aplica.
- **REQ-IND-20.** Si los datos de una muestra se mapean contra una referencia
  distinta de su referencia estructural, la región de generación no debe
  reinterpretarse en las coordenadas de la referencia de mapeo.
- **REQ-IND-21.** Al generar eventos, el sistema debe aplicar el valor espacial
  aplicable conforme a la jerarquía interna espacial de REQ-CONF-117 y resolver
  su coincidencia con una sustitución individual mediante REQ-CONF-118.
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
  entre individuos, condicionadas al modelo elegido: Hardy-Weinberg binomial o
  multinomial para ploidías 1–10, F automático solo en diploides o una
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
  formar parte de la población. En la creación directa de F2, la omisión de
  `individual_id` exige un `founder_key` aportado, único y estable dentro de su
  población de origen; el sistema deriva el identificador desde esa clave y las
  causas DTC-039. Dos individuos sin ID ni claves fundadoras distintas se
  rechazan: el orden físico del documento nunca crea una identidad.
- **REQ-IND-26.** El mismo identificador de individuo debe usarse de forma
  coherente en su estado genotípico, sus observaciones, sus relaciones parentales
  o de linaje y todos los archivos y metadatos que le correspondan.
- **REQ-IND-27.** Para cada evento de recombinación, el sistema debe determinar
  el perfil efectivo del individuo: el perfil específico de su grupo cuando sea
  aplicable o, en caso contrario, el perfil general.
- **REQ-IND-28.** El estado genotípico interno de cada individuo debe representar
  por separado cada copia cromosómica o haplotipo según su ploidía y conservar la
  fase hasta el final de la simulación. Toda `N` heredada o incorporada mediante
  una secuencia explícita o copiada debe conservarse sin interpretarse como un
  alelo concreto.
- **REQ-IND-29.** Cada SNP, *indel* y variante estructural debe asignarse a uno o
  varios haplotipos concretos del individuo, según su genotipo.
- **REQ-IND-30.** Los haplotipos de un descendiente deben reconstruirse a partir
  de los haplotipos parentales y los eventos de recombinación, conservando el
  origen parental de sus segmentos. La recombinación debe transportar sin
  cambios las `N` incluidas en un segmento estructuralmente válido.
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
  con el mismo nivel de trazabilidad que las demás clases de variación. Un motivo
  que contenga `N` o un locus cuyo motivo no pueda determinarse sin `N` debe
  excluirse de los candidatos aleatorios.
- **REQ-IND-35.** Cada evento de mutación de una repetición en tándem debe usar la
  semilla de la ejecución, asignarse a un haplotipo concreto y actualizar su
  estado genotípico añadiendo exactamente una unidad completa del motivo o
  eliminando exactamente una, según el sesgo efectivo entre expansión y
  contracción.
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
- **REQ-IND-53.** En la primera versión, todos los individuos de una misma
  unidad reproductiva deben tener la misma referencia estructural. Una unidad
  reproductiva es el conjunto de individuos que pueden elegirse como parentales
  para producir una misma descendencia: el conjunto reproductor de una
  subpoblación en una generación o los parentales de un evento de cruce. Las
  subpoblaciones o ramas que no intercambien individuos ni se reproduzcan entre
  sí pueden usar referencias estructurales distintas. Antes de iniciar la
  simulación o ejecutar un evento demográfico que forme una unidad
  reproductiva, el sistema debe rechazar cualquier mezcla e identificar los
  individuos, las referencias y el evento afectados. Esta regla no limita las
  referencias de mapeo que puedan seleccionarse para observar esos individuos.
- **REQ-IND-54.** Al combinar cambios explícitos y aleatorios, el sistema debe
  aplicar primero todos los eventos explícitos en el orden en que se declararon.
  Cada evento explícito debe recibir el estado producido por el anterior.
- **REQ-IND-55.** Después de los eventos explícitos, el sistema debe sortear por
  separado el número realizado de eventos aleatorios de cada clase y procesar
  las clases en este orden fijo: duplicaciones, deleciones, inversiones,
  translocaciones, repeticiones en tándem, *indels* y SNP. Los elementos
  transponibles no forman una clase de v1.
- **REQ-IND-56.** Cada evento aleatorio debe elegir su ubicación y, cuando
  proceda, su origen exclusivamente entre los candidatos válidos en el estado
  inmediatamente anterior. Tras aceptar un evento, el sistema debe actualizar
  el espacio válido antes de elegir el siguiente. No debe sortear candidatos
  inválidos para reintentarlos ni imponer un límite arbitrario de intentos.
- **REQ-IND-57.** Si el espacio válido actualizado no permite completar el
  número sorteado de una clase, la generación debe fallar sin reducir
  silenciosamente el resultado ni convertir la insuficiencia en advertencia. El
  diagnóstico y el manifiesto deben conservar por clase el número sorteado, el
  aceptado antes del fallo y el finalmente activo, sin publicar la generación
  fallida como verdad completa.
- **REQ-IND-58.** `N` no debe ser origen ni destino de un SNP. Un SNP explícito
  sobre `N` debe rechazarse y la generación aleatoria debe excluir esa posición;
  el modelo SNP solo puede producir A, C, G o T. Una deleción, duplicación,
  inversión, translocación, desplazamiento u otra operación puramente
  estructural puede contener `N` en el interior de su intervalo y debe
  conservarla al copiar o mover el segmento. Una `N` solo invalida un ancla,
  extremo o candidato cuando la operación necesita conocer allí una base
  concreta.
- **REQ-IND-59.** Una inserción aleatoria nueva debe contener solo A, C, G y T.
  Una inserción configurada explícitamente o copiada desde una secuencia
  existente puede contener `N` y debe conservarla exactamente.
- **REQ-IND-60.** Cada fila de una matriz de frecuencias SNV debe identificar un
  alelo `ALT` A, C, G o T distinto del `REF` A, C, G o T de la referencia exacta,
  con `contig`, `position`, `ref`, `alt`, `population_id` y `alt_af`. La única
  semántica aceptada es la frecuencia ALT; no se aceptan ni convierten MAF, DAF,
  `AC/AN` o valores ausentes. `position` usa base 1 cerrada y el nombre de
  secuencia debe coincidir literalmente con el FASTA. Un locus multialélico
  suministrado puede ocupar varias filas, pero no puede repetir un alelo ni
  discrepar en `REF`; se exige una fila por combinación `ALT`×población. Cada
  `alt_af` debe ser un `f64` finito en `[0,1]` y puede usar notación decimal o
  científica. Para cada población, la suma ALT no debe superar 1 fuera de la
  tolerancia general; el resto es la frecuencia de referencia. Un exceso dentro
  de la tolerancia se normaliza y registra; uno mayor se rechaza.
- **REQ-IND-61.** Al generar la matriz, el sistema debe seleccionar loci válidos
  sin reemplazo conforme al número de candidatos y a la distribución espacial,
  obtener `REF` de la referencia, elegir un único `ALT` mediante el espectro de
  sustituciones y asignar frecuencias por población mediante el espectro
  alélico efectivo. F2 genera solo loci bialélicos. `uniform` es el único perfil
  espacial incorporado; `weighted_regions` usa intervalos base 1 cerrados no
  solapados, pesos finitos no negativos y un `outside_weight` explícito. La
  probabilidad por base es proporcional a su peso, peso cero la excluye y debe
  quedar espacio positivo. En modo `exact`, cada población realiza una
  permutación independiente del conjunto completo de `variant_id` canónicos
  antes de asociar las cuotas del espectro. No debe situar SNV sobre `N`.
- **REQ-IND-62.** Para crear individuos desde una matriz de frecuencias, el
  sistema debe muestrear en cada locus los alelos de sus copias conforme a la
  población de origen y al modelo genotípico efectivo. Con Hardy-Weinberg, las
  copias se muestrean de forma independiente; en diploides con `F`, deben usarse
  las fórmulas bialélicas o multialélicas ya fijadas en REQ-EVO-47/48. La
  composición obtenida debe permutarse reproduciblemente entre copias y
  conservar fase aunque la matriz no defina asociación entre loci. Los loci se
  fasean independientemente y esta fase conocida no representa LD. Un genotipo
  explícito incompatible bloquea por defecto; la opción específica
  `--allow-incompatible-genotypes` puede conservarlo como restricción forzada y
  trazada sin modificar la matriz, pero nunca relaja referencia, coordenadas,
  alelos, ploidía, dosis, formato, sumas o publicación segura.
- **REQ-IND-63.** La misma matriz o sus mismos parámetros generadores, la misma
  semilla, los mismos IDs de población e individuo y la misma versión deben
  reproducir matriz, genotipos y fase. Añadir individuos no debe cambiar la
  matriz ni los individuos ya identificados; añadir una población no debe
  cambiar posiciones, alelos ni frecuencias preexistentes y reordenar filas no
  debe cambiar ningún resultado. Deben existir subsemillas separadas para
  selección de loci, `ALT` por variante, frecuencia por población/variante,
  genotipo por población/individuo/variante y fase por individuo/variante. Sus
  nombres, versión y algoritmo consumidor deben registrarse. El sorteo de ALT
  debe identificarse solo mediante datos ya conocidos antes de obtener ALT. El
  modo `exact` debe separar también el
  flujo que asigna filas a variantes del que sortea el valor dentro del
  intervalo por población y variante. No se promete estabilidad de prefijo al
  cambiar el número o la distribución espacial de loci.
- **REQ-IND-64.** En una misma población inicial, el modo de matriz de frecuencias
  SNV y la generación inicial de SNP por cantidad o densidad individual son
  mutuamente excluyentes. Los SNV explícitos pueden actuar como restricciones
  individuales, pero una contradicción con la matriz o el modelo genotípico debe
  rechazarse antes de crear individuos. Las mutaciones *de novo* posteriores
  conservan sus propias tasas y no modifican retroactivamente la matriz inicial.
- **REQ-IND-65.** Para cada población y locus, el sistema debe conservar las
  frecuencias objetivo de la matriz, el modelo Hardy-Weinberg o `F`, el número de
  copias observadas, los recuentos y frecuencias realizados y, cuando proceda,
  los genotipos esperados y realizados. La variación por tamaño finito debe
  informarse y no corregirse modificando genotipos ya sorteados.

  Si una deleción directa elimina el locus en una copia, esa copia debe figurar
  separadamente como ausente por deleción. Los recuentos y frecuencias REF/ALT
  realizadas usan solo las copias que conservan la coordenada; el resumen debe
  declarar tanto el denominador efectivo como el número total de copias.

El espectro de frecuencias alélicas es una tabla versionada de componentes no
solapados. Cada componente es un valor `fixed` de `[0,1]` o un intervalo
`[min_af,max_af)` con peso y muestreo `linear_uniform` o `log_uniform`; se
permite `max_af=1`, pero ningún intervalo continuo produce 1 y
`log_uniform` exige un mínimo positivo. Los extremos exactos `alt_af=0` y
`alt_af=1` solo proceden de componentes `fixed`. Los pesos son
finitos y no negativos, al menos uno es positivo y se normalizan en orden
canónico. El modo predeterminado `stochastic` sortea el componente de cada
candidato; el modo `exact` reparte candidatos por población mediante mayores
restos y resuelve empates por el orden estable de las filas. A continuación,
cada población permuta de forma reproducible e independiente todos los
`variant_id` canónicos con una subsemilla exclusiva de asignación y aplica las
cuotas siguiendo esa permutación. El orden genómico o canónico de variantes no
puede emplearse para esta asociación porque introduciría correlación espacial.
El valor de un intervalo se sortea después mediante otro flujo identificado por
población y `variant_id`; una fila `fixed` participa en la asignación sin
consumir RNG para el valor y una fila con recuento cero no consume RNG. Las
poblaciones pueden asignar a un mismo locus frecuencias distintas, cero o uno.
El perfil incorporado
`human_1000g_global_sfs_approx_v1` usa pesos 64:12:8 para
`[1/5008,0.005)`, `[0.005,0.05)` y `[0.05,1)`. Los dos primeros intervalos usan
muestreo logarítmico y el tercero, uniforme lineal; DTC-049 normaliza los pesos.
`1/5008` representa una copia entre las 5.008 copias autosómicas de las 2.504
muestras publicadas y solo es el mínimo de esta receta aproximada, no un límite
general. El perfil no incluye una masa fija en cero ni frecuencias positivas
inferiores; un perfil suministrado puede incluirlas conforme al contrato general.

El espectro de sustituciones canónico contiene los 12 cambios dirigidos. El
perfil predeterminado `uniform_snv12_v1` da el mismo peso a los tres `ALT` de
cada `REF`. El ejemplo seleccionable `human_harris2017_snv12_approx_v1` resume
de forma descriptiva Harris y Pritchard (2017): normaliza por separado las filas
AFR/AMR/EUR/EAS/SAS de sus 96 contextos, promedia las cinco con el mismo peso,
agrega a seis clases colapsadas por hebra y reparte cada clase por igual entre
la sustitución canónica y su complemento inverso antes de normalizar los 12
cambios. No se presenta como tasa universal de mutación. Una tabla propia
prevalece. El atajo `ti_tv_ratio=R`, finito y no
negativo, asigna peso `2R` a la transición y peso 1 a cada transversión; excluir
todas las transversiones exige una tabla explícita. La configuración efectiva
siempre conserva los 12 pesos resultantes.

Un `ALT` con `alt_af=0` en todas las poblaciones se elimina de la vista canónica
activa y se registra; en un locus multialélico se conservan los demás `ALT` y el
locus solo desaparece si ninguno queda. Un `ALT` fijado, incluso en todas las
poblaciones, se conserva. Los eliminados siguen contando entre los candidatos
generados y no se reemplazan. Una matriz activa vacía es válida, pero `run` debe
exigir la confirmación reproducible `--allow-empty-snv-matrix`; Marimo ofrece la
misma confirmación explícita.

Una matriz suministrada se prepara sin modificar el original. La entrada mínima
contiene las seis columnas de REQ-IND-60 y puede conservar `source_id`; la orden
auxiliar calcula IDs canónicos, valida contra el FASTA, ordena, filtra y publica
atómicamente `population_snvs.tsv` y `population_snvs.metadata.json`. El TSV
canónico tiene exactamente la cabecera `variant_id`, `contig`, `position`, `ref`,
`alt`, `population_id`, `alt_af`, `source_id`. `variant_id` se deriva de
`contig`, `position`, `ref` y `alt`, por lo que es común a todas las poblaciones;
`source_id` es opcional en la entrada y se publica vacío cuando no existe. Las
filas se ordenan por el orden de secuencias del FASTA, posición ascendente,
`ref`, `alt` y `population_id` léxico. El JSON declara versión, la identidad
preparada de referencia (`reference_id` y huella canónica de
ID+secuencias+bases) y, por separado, el SHA-256 físico observado del FASTA,
además de coordenadas y poblaciones,
procedencia y el SHA-256 del TSV; por población puede conservar descripción,
tamaño muestral o número de alelos, fuente, filtros y fecha o versión. En una
matriz empírica esos campos se proporcionan o se declaran explícitamente como
`unknown`, nunca se inventan. Configuración y manifiesto conservan a su vez el
SHA-256 del JSON. `run` solo consume este par canónico y vuelve a comprobar el
FASTA. Una orden independiente y la ruta integrada pueden generar el mismo par
con los mismos bytes; este es reutilizable sin volver a sortear la matriz.

Para F2.2.1, `population_snvs.metadata.json` usa la clase
`population_snvs_metadata` y `schema_version = "1"`. Contiene como campos
obligatorios
`coordinate_system = "base1_closed"`, `reference` (`reference_id`,
`canonical_fingerprint`), `reference_fasta` (`sha256`) y `matrix` (`path =
"population_snvs.tsv"`, `sha256`), además de una lista ordenada de
`populations`. Cada población declara `population_id` y puede declarar
`description`, `sample_size`, `allele_count`, `source`, `filters` y
`data_version`; cuando un dato empírico no se conoce se escribe `"unknown"`.
`provenance` declara si el recurso fue `supplied` o `generated`, sin duplicar
las filas TSV. Los campos adicionales se rechazan hasta que una fase posterior
los necesite.

La configuración general usa una única tabla `[population_snvs]` con
`mode = "supplied"`, `"reuse"` o `"generate"`. El modo `supplied` exige
`input_path`; `reuse` exige `metadata_path`; `generate` no admite rutas y sus
parámetros se incorporan con su primer consumidor en F2.2.5. Cualquier ruta no
pertinente al modo, ambos tipos de ruta o un modo desconocido se rechaza. La
validación de valores `alt_af`, sus sumas, tolerancias y normalización se
incorpora exclusivamente en F2.2.2, para que F2.2.1 no duplique esa autoridad.

En el modo `generate`, `[population_snvs.generate]` declara `count_mode` como
`exact` o `poisson`, `spatial_mode` como `uniform` o `weighted_regions`, y
`substitution_profile = "uniform_snv12_v1"`. Exige exactamente uno de
`locus_count` (solo `exact`), `mean_loci` (solo `poisson`) o
`mean_density_per_base` (solo `poisson`). Cada entrada de
`[[population_snvs.generate.populations]]` declara `population_id`,
`frequency_model = "uniform"`, `min_alt_af` y `max_alt_af`; genera valores en
el intervalo semiabierto `[min_alt_af,max_alt_af)`, sin producir 0 ni 1.
`weighted_regions` añade `outside_weight` y entradas con `contig`, `start`,
`end` y `weight`, en coordenadas base 1 cerradas. F2.2.5 no admite otros
perfiles de frecuencia o sustitución: los espectros fijo, lineal, logarítmico y
los perfiles humanos pertenecen a F2.2.6.

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
número realizado o sorteado es el número entero de eventos obtenido de la
distribución en una ejecución concreta. Debe distinguirse del número aceptado,
que cuenta los eventos situados antes de un posible fallo, y del número
finalmente activo, que cuenta los que permanecen en la verdad después de aplicar
los eventos posteriores.

La densidad relaciona la cantidad media esperada con una longitud de secuencia y
puede expresarse por base, kilobase o megabase. Una kilobase equivale a mil bases
y una megabase, a un millón de bases. La longitud válida de un ámbito es la
longitud en la que se pueden situar eventos de la clase correspondiente después
de aplicar las restricciones pertinentes. La cantidad media esperada derivada
de una densidad se obtiene multiplicando esa densidad por la longitud válida del
ámbito. Cada clase de variante puede usar su propio modo y, si usa densidad, su
propia unidad.

Una distribución probabilística asigna probabilidades a los posibles números
realizados a partir de la cantidad media esperada. En v1, la distribución de
Poisson es la única admitida para el recuento de variantes aleatorias. Los
eventos explícitos y el modo controlado expresan resultados concretos y no son
distribuciones probabilísticas alternativas.

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

Una población de origen de la matriz es un perfil estadístico identificado, no
una población simulada con generaciones, miembros o demografía. En la vía
directa sirve para elegir una columna de frecuencias al crear cada individuo;
no crea por sí misma un escenario evolutivo ni hace que esos individuos formen
una unidad reproductiva.

El espectro de frecuencias alélicas describe qué proporción de los loci usa
cada frecuencia o intervalo de frecuencias. No determina las posiciones ni los
cambios `REF`→`ALT`, y tampoco contiene asociación entre loci. Un espectro
suministrado debe expresar valores o intervalos ordenados, no solapados y con
pesos finitos no negativos que, tras aplicar la política numérica general, sumen
1. V1 solo admite `alt_af`, la frecuencia del alelo alternativo; rechaza
frecuencia derivada, MAF y cualquier entrada que exija inferir o convertir su
semántica.

Los perfiles empíricos dependen de la población, el tamaño de muestra, la región
del genoma y los filtros. Como ejemplo documentado, 1000 Genomes fase 3 observó
aproximadamente 64 millones de variantes autosómicas con frecuencia menor de
0,5 %, 12 millones entre 0,5 % y 5 % y 8 millones por encima de 5 %. ExAC, sobre
exomas de 60.706 personas, informó que el 99 % de sus variantes de alta calidad
tenía frecuencia menor de 1 % y que el 54 % eran *singletons*. Estos datos
justifican perfiles con predominio de frecuencias bajas, pero no constituyen un
perfil humano universal ni autorizan a mezclar estudios. Todo ejemplo incorporado
debe conservar fuente, versión, población, muestra, región, filtros, definición
de frecuencia y transformación aplicada. El SFS incorporado usa solo los tres
recuentos agregados publicados y se etiqueta como aproximación GLOBAL; no
infiere perfiles de superpoblación. El SNV12 incorporado usa [Harris y Pritchard
2017](https://elifesciences.org/articles/24284) y sus [datos fuente
versionados](https://doi.org/10.7554/eLife.24284.004), atribuidos bajo CC BY. El
fichero fuente tiene 3.644 bytes y SHA-256
`970cc71578fc008ef4e0d75c7c460ab874073da7c5b5226c526b387627c97a7f`.
Fuentes primarias adicionales: [1000 Genomes fase
3](https://www.nature.com/articles/nature15393) y
[ExAC](https://www.nature.com/articles/nature19057). La documentación oficial
de IGSR describe las frecuencias de alelos de 1000 Genomes, pero F2 no procesa
su VCF ni deriva perfiles incorporados por superpoblación.

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
  cambios al azar. Los SNV iniciales deben muestrearse desde la matriz y la
  población de origen elegidas; las demás clases admitidas usan una cantidad
  media esperada o una densidad por base, kilobase o megabase, con valores
  globales y valores específicos por individuo.

En F2.3.4, la forma mínima autorada para esos cambios concretos es
`[[direct_events]]`: todos los registros declaran `event_id`, `individual_id`,
`phase`, `kind` y `contig`; un SNP usa `position_1based`, `ref` y `alt`, una
inserción usa `point` y `bases`, y una deleción `start` y `end`. El registro se
resuelve contra una cohorte inicial y su ploidía, mantiene el orden declarado
por individuo y se valida contra matriz, genotipo y fase antes de publicar
artefactos. Esta fase no materializa ni publica individuos; los resultados y
estadísticas realizados siguen siendo F2.3.5.
- **REQ-DIR-03.** La creación directa debe admitir para cada clase generada por cantidad o densidad
  una distribución espacial uniforme o regiones con una cantidad media
  esperada o densidad mayor o menor que el valor general, incluidas regiones con
  densidad cero, tanto globales como específicas por individuo.
- **REQ-DIR-04.** En la creación directa o aleatoria de individuos, las
  diferencias de número de copias de las repeticiones en tándem deben poder
  generarse exclusivamente con el modelo escalonado simple. Cada mutación debe
  añadir o eliminar una unidad completa; la configuración solo puede cambiar el
  sesgo entre ambas direcciones.
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

El equilibrio de Hardy-Weinberg y el coeficiente automático de endogamia F son
opciones distintas. Hardy-Weinberg se admite para ploidías 1–10. En un locus
bialélico con frecuencia alternativa `p` y ploidía `P`, la probabilidad de dosis
`k` es `C(P,k) p^k (1-p)^(P-k)`, donde `C(P,k)` es el coeficiente binomial que
cuenta las formas de distribuir `k` copias alternativas entre `P` copias; en un
locus multialélico se usa la distribución multinomial equivalente sobre las
composiciones válidas. El modelo automático F
solo se admite en diploides. Con dos alelos de frecuencias `p` y `q`, usa
`P(AA)=p²+Fpq`, `P(AB)=2pq(1-F)` y `P(BB)=q²+Fpq`. En un locus multialélico
diploide, para cada alelo `i` usa `P(ii)=p_i²+F p_i(1-p_i)` y para cada par
`i<j`, `P(ij)=2(1-F)p_i p_j`. F=0 coincide con Hardy-Weinberg y F=1 produce
homocigosis completa. En haploides y poliploides solo se admite Hardy-Weinberg
sin F o una distribución explícita. F realizado se calcula únicamente en
diploides; en las demás ploidías figura como «no aplicable».

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

El tamaño N objetivo es el número de intentos reproductivos que deben efectuarse
para una población o subpoblación. Cada intento exitoso crea exactamente un
descendiente; N realizado es el número de descendientes creados. Ambos pueden
diferir si la migración, la selección o el propio sorteo reproductivo dejan sin
opciones válidas uno o más intentos. El tamaño efectivo, representado por Ne,
describe el efecto de la contribución reproductiva y puede diferir de N. V1 no
admite Ne objetivo. Ne realizado es el valor estimado después de la
reproducción y será exclusivamente un Ne
reproductivo: se estimará a partir de los progenitores que dejaron descendencia
y de la desigualdad o varianza de sus contribuciones. No se estimará a partir de
la pérdida de heterocigosidad, el incremento de endogamia ni el desequilibrio de
ligamiento; esos estimadores quedan para futuras ampliaciones. Una población puede
dividirse en subpoblaciones, que son grupos reproductivos internos con parámetros
demográficos o evolutivos propios y que pueden estar conectados por migración.
En este documento, «población» designa siempre el conjunto completo simulado. Si
solo existe un grupo reproductivo, la población contiene una única
subpoblación. La migración es el paso de un individuo de una subpoblación a otra.
N objetivo fija los intentos y N realizado registra sus éxitos. Ne realizado se
calcula siempre que sea evaluable; en caso contrario se registra como «no
evaluable» junto con la causa.
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
  exportación y, cuando genere FASTQ, BAM o CRAM, su único experimento de
  secuenciación independiente.
- **REQ-EVO-15.** La creación, reproducción, muestreo y desaparición de
  individuos deben conservar sus identificadores estables y las relaciones con
  sus parentales cuando existan.
- **REQ-EVO-16.** N objetivo debe poder variar entre generaciones, fases o
  tramos para representar expansiones, reducciones y cuellos de botella. El
  sistema debe aplicar en cada intervalo el valor efectivo de intentos
  reproductivos que resulte de la configuración.
- **REQ-EVO-17.** Cada subpoblación debe poder tener su propio N objetivo,
  modelo de apareamiento y demás parámetros evolutivos, variables por
  generación, fase o tramo. También debe calcular su Ne realizado cuando sea
  evaluable o registrar la causa si no lo es.
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
- **REQ-EVO-24.** La intensidad de la deriva debe resultar de N, de la
  contribución reproductiva realizada y de la estructura poblacional, sin un
  objetivo de Ne ni un paso adicional. El sistema debe registrar las
  frecuencias alélicas realizadas por generación para poder observarla.
- **REQ-EVO-25.** En la primera versión, la selección debe ser exclusivamente de
  viabilidad y usar un único SNP bialélico estable para toda la ejecución. Los
  valores genéticos por dosis del alelo alternativo y la función que transforma
  el valor del carácter en aptitud pueden variar por fase, intervalo y
  subpoblación sin cambiar la identidad del locus.
- **REQ-EVO-26.** En cada generación, el sistema debe aplicar exclusivamente el
  orden obligatorio de REQ-REP-24. Un proceso opcional desactivado o no
  aplicable debe omitirse en su etapa sin mover, adelantar ni reinterpretar las
  demás.
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
  evolución debe aplicar exclusivamente el modelo escalonado simple, con la tasa
  y el sesgo entre expansión y contracción efectivos de la fase, la
  subpoblación, la región y el locus.
- **REQ-EVO-31.** Las mutaciones internas de repeticiones en tándem solo deben
  aplicarse en la evolución cuando estén activadas, usando sus tipos, tasas y
  distribuciones efectivos.
- **REQ-EVO-32.** La población inicial debe contener al menos una subpoblación.
  Cada subpoblación debe poder elegir una distribución inicial distinta. La suma
  de sus tamaños debe coincidir con el N total de la población.
- **REQ-EVO-33.** El equilibrio de Hardy-Weinberg debe ser un modo opcional de
  inicialización por locus para ploidías 1–10. En un locus bialélico con
  frecuencia alternativa `p` y ploidía `P`, la dosis `k` debe seguir
  `C(P,k) p^k (1-p)^(P-k)`. En loci multialélicos debe usar la distribución
  multinomial equivalente sobre todas las composiciones de dosis válidas.
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
  como opción automática exclusivamente en diploides. Debe poder definirse para
  la población completa y sustituirse por subpoblación; no debe poder definirse
  por locus, región o bloque.
- **REQ-EVO-47.** F debe estar entre 0 y 1, inclusive. F=0 debe corresponder al
  modelo de Hardy-Weinberg y F=1 debe producir homocigosis completa. En un locus
  bialélico diploide con frecuencias `p` y `q` deben usarse
  `P(AA)=p²+Fpq`, `P(AB)=2pq(1-F)` y `P(BB)=q²+Fpq`. Los valores de F fuera del
  intervalo deben ser un error.
- **REQ-EVO-48.** En un locus multialélico diploide, el modelo automático F debe
  usar `P(ii)=p_i²+F p_i(1-p_i)` para cada homocigoto y
  `P(ij)=2(1-F)p_i p_j` para cada heterocigoto con `i<j`. Configurar F para
  ploidía 1 o mayor que 2 debe rechazarse antes de simular, identificando ámbito
  y ploidía. V1 no debe inferir una adaptación poliploide.
- **REQ-EVO-49.** Una distribución genotípica explícita por locus o bloque debe
  sustituir localmente al modelo general de F. Ambos métodos solo son
  incompatibles si intentan aplicarse al mismo locus o bloque. F debe seguir
  aplicándose al resto de los loci diploides. En haploides y poliploides debe
  admitirse Hardy-Weinberg sin F o una distribución explícita de dosis o
  genotipos; la distribución explícita prevalece localmente.
- **REQ-EVO-50.** F debe ser un objetivo esperado de la población o
  subpoblación diploide, no un valor exacto por individuo. La homocigosis
  realizada debe poder variar estocásticamente. Para fijar un individuo concreto
  deben usarse genotipos completos con fase.
- **REQ-EVO-51.** El modo automático de F solo debe ajustar las probabilidades
  genotípicas y el exceso esperado de homocigosis. No debe generar por sí solo
  tramos realistas de homocigosis. Los ROH deben surgir de pedigrí, de
  autofecundación o de evolución con recombinación.
- **REQ-EVO-52.** El sistema debe registrar el F objetivo general o por
  subpoblación diploide, su procedencia, el valor efectivo y el F realizado por
  subpoblación y generación registrada. Debe estimar este último únicamente en
  diploides, comparando heterocigosidad observada y esperada según las
  frecuencias alélicas. Un valor negativo debe conservarse como posible exceso
  de heterocigotos. En ploidía 1 o mayor que 2 debe registrar F realizado como
  «no aplicable», nunca como cero ni como una estimación adaptada.
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
  resueltas, referencias y subsemillas.
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
  muestral debe estimarse con las frecuencias y la heterocigosidad de la muestra
  solo en diploides; en ploidía 1 o mayor que 2 debe figurar como «no
  aplicable».
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
  selección, tasa de autofecundación y modelo reproductivo. Por
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
- **REQ-EVO-74.** Una fase posterior debe admitir un panel suministrado formado
  por bloques haplotípicos locales. El paquete canónico debe constar exactamente
  de cinco ficheros enlazados: `panel_blocks.tsv`, `panel_sites.tsv`,
  `panel_haplotype_alleles.tsv`, `panel_frequencies.tsv` y
  `panel.metadata.json`. `panel_blocks.tsv` debe contener una fila por bloque con
  `block_id`, `reference_id`, `sequence_id`, `start` y `end`. Las coordenadas
  deben usar la convención común basada en 1 y cerrada. Los bloques no deben
  solaparse y todos sus sitios deben quedar dentro de su intervalo. El paquete
  debe estar ligado a una única referencia: todos los `reference_id` de
  `panel_blocks.tsv` deben coincidir exactamente con
  `panel.metadata.json.reference.reference_id`; el
  `panel.metadata.json.reference.sha256` debe coincidir con el SHA-256
  recalculado sobre los bytes del FASTA autoritativo; y todo `sequence_id` debe
  existir en ese FASTA. Una mezcla de referencias, un ID o digest divergente o
  una secuencia inexistente debe rechazarse antes de interpretar o crear
  individuos.

  Los cuatro TSV deben usar UTF-8 sin BOM, saltos LF y tabulador como único
  separador. Deben terminar en LF, contener una única cabecera exacta y no
  admitir comillas, líneas vacías, campos vacíos ni columnas adicionales. Sus
  esquemas, en orden, deben ser los siguientes. En la tabla, `\t` representa un
  único carácter de tabulación en la cabecera, no los dos caracteres barra y
  letra t:

  | Fichero | Cabecera exacta |
  |---|---|
  | `panel_blocks.tsv` | `block_id\treference_id\tsequence_id\tstart\tend` |
  | `panel_sites.tsv` | `block_id\tsite_id\tposition\tref\talt` |
  | `panel_haplotype_alleles.tsv` | `block_id\thaplotype_id\tsite_id\tallele_index` |
  | `panel_frequencies.tsv` | `population_id\tblock_id\thaplotype_id\tfrequency` |

  Los identificadores deben conservar la sintaxis ASCII de DTC-039. `start`,
  `end` y `position` deben escribirse como enteros decimales positivos sin signo
  ni ceros iniciales; `allele_index`, como entero decimal no negativo con las
  mismas reglas. `frequency` debe usar la serialización decimal mínima que
  recupera exactamente el mismo `f64`, conforme a DTC-049, sin notación no finita.

  `panel_sites.tsv` debe contener el catálogo de sitios ordenado por bloque, con
  una fila por sitio. `ref` debe contener una base A/C/G/T y `alt`, uno o varios
  ALT A/C/G/T separados por comas, sin espacios, comillas ni alelos repetidos o
  iguales a `ref`. El orden textual de `alt` debe definir `allele_index` 1..N;
  `ref` debe corresponder siempre al índice 0. Cada haplotipo debe contener los alelos faseados de un único
  bloque y su identidad debe ser local a ese bloque: la
  identidad inequívoca es la pareja `block_id` + `haplotype_id`, y reutilizar un
  mismo `haplotype_id` en otro bloque no debe relacionarlos. V1 solo debe admitir
  SNP suministrados, bialélicos o multialélicos: `REF` y cada `ALT` deben ser
  exactamente una base A, C, G o T. Debe rechazar `N`, *indels*, repeticiones,
  variantes estructurales, alelos simbólicos, *breakends* y sitios repetidos o
  solapados.

  `panel_frequencies.tsv` debe contener exactamente una fila por
  (`population_id`, `block_id`, `haplotype_id`), incluidos los valores 0 y 1.
  Para cada población y bloque, las frecuencias deben ser finitas, pertenecer a
  [0,1] y sumar 1 conforme a DTC-049. Una fila ausente o duplicada debe ser un
  error, no debe interpretarse como frecuencia 0 y no debe existir una frecuencia
  global ni herencia de frecuencias entre poblaciones. Las frecuencias alélicas
  marginales deben poder derivarse sin otra fuente de verdad. La tabla
  canónica `panel_haplotype_alleles.tsv` debe ser densa: para cada pareja
  (`block_id`, `haplotype_id`) debe contener exactamente una asignación para cada
  `site_id` del bloque. `allele_index=0` representa `REF` y los valores de 1 a N,
  los `ALT` declarados para el sitio en ese mismo orden. Una asignación ausente,
  duplicada o ajena al catálogo debe ser un error; una ausencia no equivale a
  `REF` y v1 no admite alelos desconocidos ni ausentes en el panel. Dentro de un
  mismo `block_id`, dos `haplotype_id` no deben tener el mismo vector completo y
  ordenado de `allele_index`. Si varios haplotipos son equivalentes, la
  validación debe rechazarlos e identificar el bloque y todos los IDs
  implicados; no debe fusionarlos ni sumar sus frecuencias. La misma secuencia de
  índices puede aparecer en bloques distintos porque representa sitios
  diferentes. Un haplotipo de frecuencia 0 sigue siendo válido si su vector es
  único. Cada bloque debe declarar al menos dos sitios SNP y al menos dos
  vectores HAP distintos. Una población puede asignar frecuencia 1 a un único
  HAP y 0 a todos los demás; el catálogo global sigue siendo válido incluso si
  todas las poblaciones actuales lo hacen, y las filas con 0 deben conservarse.
  Un sitio cuyo mismo ALT esté fijado en todos los HAP debe conservarse. En
  cambio, un sitio con `allele_index=0` en todos los HAP no representa ninguna
  variante respecto a la referencia y debe rechazarse como redundante. La matriz de
  SNP independientes `alt_af` no puede contener ningún locus situado dentro del
  intervalo completo de un bloque, aunque esa posición no sea un `site_id` del
  panel. La validación debe detectar el conflicto antes de crear individuos, sin
  prioridad ni sobrescritura, e identificar como mínimo la referencia,
  secuencia, posición y fila o identidad de matriz conflictiva, además del
  `block_id` y su intervalo. Fuera de los intervalos de bloque debe regir
  DTC-054. `panel.metadata.json` debe declarar la clase y versión del recurso,
  `panel_id`, la referencia y su SHA-256, y la convención de coordenadas basada
  en 1 y cerrada. Para cada uno de los cuatro TSV debe registrar clase, ruta,
  tamaño y SHA-256. También debe registrar la procedencia y la herramienta
  creadora cuando existan. El JSON no debe duplicar las filas de las tablas.

  El orden canónico debe ser único. `panel_blocks.tsv` debe ordenar por
  `reference_id`, orden de `sequence_id` en el FASTA, `start`, `end` y
  `block_id`; `panel_sites.tsv`, por el orden canónico del bloque, `position` y
  `site_id`; `panel_haplotype_alleles.tsv`, por bloque, `haplotype_id` ASCII y
  orden canónico del sitio; y `panel_frequencies.tsv`, por `population_id` ASCII,
  bloque y `haplotype_id` ASCII. Deben ser únicas, respectivamente, las claves
  `block_id`, (`block_id`,`site_id`),
  (`block_id`,`haplotype_id`,`site_id`) y
  (`population_id`,`block_id`,`haplotype_id`). La orden de preparación debe
  rechazar cualquier entrada que no pueda llevarse sin ambigüedad a este orden.
  Una única orden auxiliar del binario Rust, conceptualmente
  `prepare-haplotype-panel`, debe recibir la referencia y las cuatro tablas,
  validar su contrato y relaciones, crear una copia nueva en orden canónico,
  calcular los SHA-256 de la referencia y de las tablas, generar
  `panel.metadata.json` y publicar atómicamente el directorio. No debe modificar
  los originales ni sobrescribir un destino existente. La simulación solo debe
  aceptar el paquete completo: debe recalcular sus sumas y verificar referencia,
  clases, versiones y rutas. Un fichero ausente o modificado o cualquier
  incompatibilidad debe detener la ejecución antes de crear individuos.
- **REQ-EVO-75.** En el modo sin recombinación, cada copia cromosómica fundadora
  debe recibir en cada bloque un haplotipo local completo e intacto conforme a
  sus frecuencias. A efectos del muestreo, cada bloque debe tratarse como un
  locus multialélico cuyas categorías son sus haplotipos locales. Para cada
  individuo de una población y ploidía P, deben asignarse P haplotipos locales a
  sus copias haplotípicas biológicas identificadas. Hardy-Weinberg mediante el
  modelo multinomial debe ser el modo predeterminado para P=1..10. El modelo F
  debe reutilizar exactamente las fórmulas multialélicas aprobadas y solo debe
  admitirse en diploides: F=0 coincide con Hardy-Weinberg y F=1 asigna a ambas
  copias el mismo haplotipo local. F debe rechazarse en P=1 o P>2. Una
  distribución explícita de combinaciones haplotípicas puede prevalecer
  localmente conforme a REQ-EVO-49.

  Cada población y bloque debe elegir uno de los dos modos normativos. El modo
  probabilístico debe ser el predeterminado: sortea cada copia según las
  frecuencias, permite reutilizar un haplotipo en sorteos distintos y produce
  recuentos realizados estocásticos. El modo controlado debe solicitar o derivar
  un vector entero de copias por haplotipo que sume C, construir ese multiconjunto
  y asignarlo reproduciblemente sin alterar sus totales; no significa muestreo
  sin reemplazo del catálogo. Cuando un modo controlado exija además recuentos
  genotípicos finales exactos, la asignación debe ser compatible con HWE, F
  diploide o la distribución explícita aplicable. Una incompatibilidad
  matemática debe rechazarse antes de crear individuos, sin ajustar cantidades
  ni relajar el modelo. Esta restricción no se aplica al modo probabilístico ni
  convierte F en emparejamiento físico, autofecundación o pedigrí. Deben
  registrarse los valores esperados, solicitados y realizados.

  Los bloques deben muestrearse independientemente. Un mismo individuo puede
  combinar haplotipos distintos en bloques diferentes y esa combinación no debe
  expresar un donante común. V1 no debe imponer un LD objetivo entre bloques,
  condicionar el haplotipo de un bloque al elegido en otro ni corregir el LD
  interbloque realizado en una población finita. El modelado posterior de LD
  entre bloques queda en DTF-044. Esta independencia entre bloques no modifica
  la construcción recombinante fundadora dentro de cada bloque, que conserva
  todas las reglas de REQ-EVO-76 y REQ-EVO-77. Cada asignación debe quedar ligada
  de forma reproducible a la copia haplotípica biológica que la recibió. Los SNP situados
  fuera de los intervalos de bloque deben proceder del
  mecanismo vigente de loci independientes de la matriz `alt_af`; esto no
  implica que las bases situadas entre bloques permanezcan iguales a la
  referencia. En modo intacto debe seleccionarse exactamente un HAP por cada
  población×individuo×copia×bloque y copiar conjuntamente todos sus
  `allele_index`, sin cortes ni mosaicos. Antes de aplicar mutaciones posteriores,
  el vector debe ser idéntico a un HAP del panel y un HAP con frecuencia 0 nunca
  debe seleccionarse. Deben registrarse población, bloque, copia biológica y HAP
  fundador. Las mutaciones de novo posteriores pueden modificar ese estado según
  el flujo general, pero deben trazarse en una capa separada y no reescribir la
  identidad ni el origen del HAP fundador.
- **REQ-EVO-76.** En el modo con recombinación, cada individuo y bloque debe
  generarse directamente; no debe existir un depósito prototipo seguido de un
  segundo muestreo ni debe presentarse la operación como una meiosis histórica.
  El productor de una copia final debe seleccionar primero un HAP primario según
  las frecuencias del panel. La configuración debe declarar
  `recombination_probability_per_interbase = r`, finita en (0,1]. Para un bloque
  de longitud L existen L−1 uniones elegibles y cada una debe sortearse
  independientemente en cada ejecución del productor. Por ello, K sigue
  `Binomial(L−1,r)` y la probabilidad de al menos un punto es
  `1−(1−r)^(L−1)`.

  Cuando K=0, el productor debe devolver intacto el HAP primario. Cuando K sea
  mayor o igual que 1, debe sortear un HAP secundario según las frecuencias de la
  población; ambas fuentes pueden coincidir. Los puntos deben ser únicos, estar
  ordenados y caer estrictamente dentro del intervalo del bloque; se deben
  admitir varios cruces. Un punto x debe representar la unión situada después
  de la base x, entre x y x+1; no debe desplazarse hasta un SNP. La salida debe
  comenzar con la fuente primaria en `start` y alternar secundaria y primaria en
  cada punto, siempre en orientación directa y sin inversión. Por ejemplo, para
  `start=100`, `end=300` y puntos 150 y 240, los segmentos deben ser primaria
  100–150, secundaria 151–240 y primaria 241–300. El productor
  debe devolver ese mosaico como una copia final. No se deben corregir ni
  descartar mosaicos para imponer recuentos de fuentes o resultados. DTC-055
  define y cierra para F6.4 un único constructor común de segmentos y su variante
  uniforme: recibe las dos fuentes y los puntos ya ordenados, alterna los
  segmentos con estas convenciones y produce la copia reconstruible. F7.1 debe
  reutilizar ese constructor; DTPD-111 solo debe decidir cómo los mapas y la
  semántica meiótica proporcionan sus puntos, sin bloquear F6.4 ni crear un
  algoritmo paralelo. Los mapas posteriores deben sustituir la probabilidad
  uniforme por probabilidades dependientes de la posición sin crear otro motor. Los
  bloques se procesan independientemente y una fuente elegida en uno no expresa
  continuidad de donante en otro. Los SNP exteriores siguen DTC-054 y el panel
  de entrada no se modifica. Tanto la fuente primaria como la secundaria deben
  seleccionarse exclusivamente del panel original mediante
  `panel_frequencies.tsv`. Un recombinante generado antes en la misma ejecución
  nunca debe convertirse en fuente, las probabilidades no deben actualizarse y
  el resultado no debe depender del orden de procesamiento. La recombinación
  acumulada, la reproducción, las poblaciones evolutivas y el uso de haplotipos
  descendientes como fuentes pertenecen a F7/F8 y quedan fuera de esta
  inicialización. Todos los cruces sorteados deben conservarse en la traza,
  incluidos los que no cambien ningún alelo observado y los pares cuyo efecto
  alélico neto se cancele.

  En diploides, para cada individuo y bloque debe sortearse primero una variable
  de Bernoulli con probabilidad F. Si indica autozigosis —dos copias que proceden
  de la misma realización genética—, se debe ejecutar una sola vez el productor
  y asignar exactamente su mismo vector a ambas copias. En caso contrario, deben
  ejecutarse dos productores independientes. Esto implementa directamente
  `P(ii)=p_i²+F p_i(1−p_i)` y `P(ij)=2(1−F)p_i p_j` para la distribución
  implícita de vectores finales. F=0 equivale a HWE y F=1 produce todos los
  individuos homocigotos, sin representar autofecundación, apareamiento ni
  pedigrí. Cuando P sea distinta de 2, deben ejecutarse P productores
  independientes y F debe rechazarse.

  Para cada población y bloque deben existir C=N×P copias finales. Estas deben
  agruparse por su vector completo y ordenado y su frecuencia realizada debe ser
  `n_h/C`. Un mosaico coincidente con un HAP del panel suma a ese vector sin
  perder su procedencia recombinante; mosaicos nuevos idénticos comparten
  categoría. Las frecuencias del panel son el objetivo para seleccionar las
  fuentes, no frecuencias creadas por un conjunto intermedio.
- **REQ-EVO-77.** Los modos con y sin recombinación deben ser mutuamente
  excluyentes para una misma unidad de creación. La configuración debe separar
  los dos ejes con nombres inequívocos:
  `haplotype_construction = "intact" | "recombinant"` para construir los
  haplotipos y `frequency_sampling = "probabilistic" | "controlled"` para
  tratar sus frecuencias. `haplotype_construction` debe poder declararse como
  valor general y mediante sustituciones por población y bloque compatibles con
  DTC-033; la regla conjunta población+bloque debe prevalecer conforme a la
  precedencia vigente. Poblaciones o bloques distintos pueden usar
  construcciones distintas, pero una misma combinación población×bloque no
  puede mezclar ambas. `intact` debe prohibir cualquier parámetro de
  recombinación. `recombinant` debe exigir
  `recombination_probability_per_interbase` finita y mayor que 0 y menor o igual
  que 1 y solo debe admitir `frequency_sampling = "probabilistic"`. La
  combinación `recombinant` + `controlled` debe rechazarse antes de crear
  individuos e identificar ambos campos y el ámbito población×bloque efectivo.
  `intact` debe admitir HWE, F diploide o una distribución genotípica explícita.
  `recombinant` debe admitir únicamente HWE o F diploide; combinarlo con una
  distribución genotípica explícita debe rechazarse antes de crear individuos e
  identificar población, bloque y campos incompatibles. Condicionar una
  distribución explícita sobre resultados recombinantes queda fuera de v1.
  Que algunas copias obtengan K=0 por azar no cambia el modo efectivo. El
  control exacto de fuentes o resultados recombinantes queda fuera de v1.
  Deben registrar las frecuencias objetivo del panel para las fuentes y los
  recuentos, frecuencias y genotipos realizados sobre las C copias finales por
  población, bloque y vector completo, además de las subsemillas. En el modo
  recombinante deben distinguir el número de ejecuciones del productor con K≥1
  y el número de copias finales derivadas de ellas: una única salida duplicada
  por F cuenta como un productor con cruce y dos copias finales derivadas. Por
  cada productor con cruce deben registrar la fuente primaria, la secundaria,
  los puntos y la procedencia ordenada de todos los segmentos. También deben
  distinguir los vectores ausentes del catálogo de entrada. La
  procedencia intacta y recombinante debe distinguirse aunque el vector final
  coincida. También deben registrarse los genotipos esperados y realizados y el
  F realizado, o «no aplicable» fuera de diploides. No se debe presentar como preservado un LD que la
  recombinación haya modificado.

  La aleatoriedad debe reutilizar exclusivamente el contrato nombrado y local de
  DTC-037. Los ámbitos deben incluir `population_id`, `block_id`, `individual_id`
  y, cuando corresponda, `copy_id`, y separar las decisiones de F/autozigosis,
  fuente primaria, puntos interbase y fuente secundaria. La misma raíz, panel,
  configuración y versión debe reproducir el resultado. Añadir o reordenar
  individuos o bloques no debe alterar los IDs existentes, y cambiar `r` en un
  bloque no debe afectar a otro. Cuando F duplique una salida, debe existir una
  sola ejecución y procedencia del productor y dos copias biológicas enlazadas a
  ella. Deben registrarse la versión de derivación, los ámbitos y las
  subsemillas; no se debe introducir otro generador aleatorio.

  Deben publicarse tres tablas obligatorias. `founder_block_assignments.tsv`
  debe contener una fila por población, individuo, copia y bloque, con
  `final_haplotype_id`, `origin_kind` y `process_id`; una salida duplicada por F
  crea dos asignaciones enlazadas al mismo proceso. Todo productor debe tener un
  `process_id` tipado conforme a DTC-039, incluso cuando K=0. `origin_kind` debe
  admitir exactamente `intact` para K=0 y `recombined` para K≥1, aunque el vector
  final coincida con un HAP del panel. Las dos copias duplicadas por F deben usar
  el mismo `process_id` y `origin_kind`; no debe existir un origen
  `autozygous_copy`.

  `recombination_segments.tsv` debe contener una fila por
  (`process_id`, `segment_index`) con bloque, `start`, `end` y HAP fuente,
  suficiente para reconstruir la salida. Un proceso K=0 debe contener un único
  segmento que cubra el bloque completo con su fuente primaria.
  `realized_haplotype_frequencies.tsv` debe contener una fila por población,
  bloque y HAP final con `copy_count`, `intact_copy_count`,
  `recombined_copy_count` y frecuencia `n_h/C`, agrupando vectores idénticos.
  Debe cumplirse `copy_count = intact_copy_count + recombined_copy_count`; una
  copia duplicada por F debe contarse dos veces en la categoría de su proceso.
  Las tres tablas deben compartir el contrato byte a byte de los TSV del panel:
  UTF-8 sin BOM, LF incluido al final, tabulador, una cabecera exacta y ausencia
  de comillas, líneas o campos vacíos y columnas extra. Sus cabeceras exactas,
  donde `\t` representa un carácter de tabulación, deben ser:

  | Fichero | Cabecera exacta |
  |---|---|
  | `founder_block_assignments.tsv` | `population_id\tindividual_id\tcopy_id\tblock_id\tfinal_haplotype_id\torigin_kind\tprocess_id` |
  | `recombination_segments.tsv` | `process_id\tsegment_index\tblock_id\tstart\tend\tsource_haplotype_id` |
  | `realized_haplotype_frequencies.tsv` | `population_id\tblock_id\thaplotype_id\tcopy_count\tintact_copy_count\trecombined_copy_count\tfrequency` |

  Todos los campos terminados en `_id`, incluido `copy_id`, deben usar su forma
  ASCII canónica conforme a DTC-039. `segment_index` y los tres recuentos deben
  usar decimal no negativo sin signo ni ceros iniciales; `start` y `end`,
  decimal positivo con las mismas reglas;
  `frequency`, la forma decimal mínima *round-trip* de DTC-049. Las filas deben
  ordenarse, respectivamente, por población, individuo, copia y orden canónico
  de bloque; por proceso y `segment_index`; y por población, orden canónico de
  bloque e ID de HAP. Esas mismas tuplas deben ser sus claves únicas. La clase y
  versión de cada tabla pertenecen a su enlace de esquema/manifiesto y no deben
  añadirse como columnas TSV. Cada tabla debe quedar enlazada desde el manifiesto
  mediante clase, versión, ruta, tamaño y SHA-256.
- **REQ-EVO-78.** Cuando un evento produzca resultados biológicos o técnicos,
  su manifiesto propio y el manifiesto raíz deben enlazar todos esos resultados
  y conservar su relación con el evento, el conjunto, la referencia y el
  artefacto correspondientes. Un evento sin una salida activada no debe inventar
  un resultado vacío para satisfacer este requisito.
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
  Debe materializarlos como parte de la verdad de familias; el tamaño del
  escenario no puede reducir este alcance.
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
- **REQ-CRU-23.** SNP, *indels* menores que el umbral estructural y repeticiones
  en tándem pueden diferir entre copias haplotípicas y parentales compatibles.
  Deben heredarse como alelos según la segregación, la ploidía y la
  recombinación, conservando fase y linaje.
- **REQ-CRU-24.** Los cruces, la autofecundación y la evolución con reproducción
  sexual genérica solo deben admitir ploidía par 2, 4, 6, 8 o 10 en cada
  cromosoma o *scaffold* reproductivo. Los parentales de un evento deben tener
  la misma ploidía efectiva en cada secuencia. Una ploidía impar o desigual
  detectable en la configuración o en parentales explícitos debe rechazarse
  antes de simular, identificando secuencia, parentales y valores.
- **REQ-CRU-25.** Cada parental debe aportar la mitad de las copias de cada
  cromosoma o *scaffold* al descendiente, de modo que este conserve la misma
  ploidía par efectiva de sus parentales.
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
- **REQ-CRU-30.** Antes de reproducir, todas las copias haplotípicas de cada
  cromosoma o *scaffold* reproductivo de ambos parentales deben tener exactamente
  el mismo estado estructural canónico: segmentos estructurales, orden,
  orientación, número de copias y puntos de ruptura. Esta regla se aplica a
  deleciones, duplicaciones, inversiones, inserciones clasificadas como
  estructurales y translocaciones. Un estado estructural puede ser compatible
  aunque sus copias procedan de eventos históricos con identificadores distintos,
  siempre que la equivalencia sea inequívoca y se conserve cada historia. Una
  diferencia o heterocigosis estructural dentro de un individuo o entre
  parentales impide usar el individuo o la pareja como opción reproductiva en
  v1; no debe intentarse una recombinación parcial.
- **REQ-CRU-31.** Las mutaciones de novo deben estar disponibles en familias y
  poblaciones segregantes mediante el mismo motor usado en evolución, pero
  desactivadas de forma predeterminada. Su activación debe ser explícita y debe
  registrar tasas, eventos realizados y trazabilidad con las mismas reglas. Una
  variante estructural de novo puede existir en un descendiente y en una
  generación observada. Su portador no debe usarse como parental posterior
  mientras el estado estructural no sea idéntico en todas sus copias y en los
  demás parentales requeridos. V1 no debe rechazar una tasa de SV de novo solo
  porque existan generaciones posteriores ni inventar un mecanismo de fijación;
  debe aplicar la compatibilidad al intento posterior y registrar la causa si
  excluye al portador.
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
  debe comprobar que sus parentales usan la misma referencia estructural y que
  son compatibles en ploidía y estado estructural canónico conforme a
  REQ-CRU-24 y REQ-CRU-30. Un parental explícito incompatible debe rechazarse
  antes de simular; en evolución, cada intento debe considerar únicamente
  individuos o parejas compatibles y registrar la incompatibilidad estructural
  cuando no queden opciones válidas.
- **REQ-CRU-45.** Toda regla de selección y toda salida solicitada en el grafo
  debe referirse a resultados que puedan producirse según sus dependencias,
  cardinalidades y datos disponibles. Una regla o salida imposible debe
  rechazarse antes de ejecutar el evento e identificar la referencia, el
  requisito incumplido y los resultados afectados.
- **REQ-POB-02.** El sistema debe permitir crear poblaciones segregantes.
- **REQ-POB-03.** El sistema debe permitir crear familias.

#### 4.3.5. Reproducción y dinámica generacional

El estado estructural canónico de una copia haplotípica describe sus segmentos
estructurales, orden, orientación, número de copias y puntos de ruptura, sin
depender de los identificadores históricos de los eventos que lo produjeron.
En v1, una secuencia solo es reproductivamente compatible si tiene ploidía par
2, 4, 6, 8 o 10 y todas sus copias haplotípicas en todos los parentales
requeridos presentan exactamente el mismo estado estructural canónico. SNP,
*indels* menores que el umbral estructural y repeticiones en tándem pueden
diferir entre copias y se transmiten como alelos; no cambian esa equivalencia.

- **REQ-REP-01.** El sistema de apareamiento debe admitir reproducción mixta
  mediante una tasa de autofecundación entre 0 % y 100 %, inclusive.
- **REQ-REP-02.** La tasa de autofecundación debe poder definirse por población,
  subpoblación, fase o intervalo de generaciones. El valor efectivo debe
  resolverse conforme a REQ-CONF-116 a REQ-CONF-121.
- **REQ-REP-03.** La tasa de autofecundación debe interpretarse como la
  probabilidad aplicada a cada evento reproductivo, no como una cuota exacta de
  eventos por generación. Las decisiones deben ser aleatorias y reproducibles.
- **REQ-REP-04.** Los eventos que no resulten en autofecundación deben ser
  alógamos. En la primera versión deben elegir de forma aleatoria uniforme dos
  individuos distintos entre los supervivientes elegibles de la población o
  subpoblación y generación.
- **REQ-REP-05.** La elegibilidad puede depender de selección de viabilidad,
  pertenencia a la población, disponibilidad y compatibilidad reproductiva de
  REQ-CRU-24 y REQ-CRU-30. La exclusión por incompatibilidad estructural debe
  identificarse como tal y no presentarse como selección de viabilidad. La primera versión no debe
  ponderar la elección por parentesco, distancia genética, genotipo ni
  apareamiento selectivo. Se exceptúan los pesos por subpoblación de origen
  configurados para una fusión conforme a REQ-EVO-67.
- **REQ-REP-06.** Un individuo debe poder participar en varios eventos
  reproductivos y producir varios descendientes en una generación.
- **REQ-REP-07 (retirado).** Límites o distribuciones configurables del número de
  descendientes por parental; capacidad aplazada fuera de v1. En v1 los
  parentales pueden reutilizarse sin límite de descendientes.
- **REQ-REP-08.** N debe ser el objetivo configurado del tamaño de cada
  nueva generación y subpoblación. En una transición reproductiva, el N efectivo
  debe fijar exactamente cuántos intentos reproductivos se procesan. Cada intento
  debe seleccionar reproduciblemente el tipo de apareamiento y los parentales
  solo entre opciones válidas, incluida la compatibilidad de ploidía y
  estructura, y, si tiene éxito, crear exactamente un descendiente. Un resultado aleatorio puede impedir que todos tengan éxito
  conforme a REQ-REP-26.
- **REQ-REP-09.** N debe ser obligatorio al crear cada población o subpoblación inicial.
  Si se omite en una generación o fase posterior, debe heredarse el último valor
  efectivo.
- **REQ-REP-10.** N solo debe cambiar mediante un nuevo valor explícito o un
  evento demográfico configurado, como cuello de botella, expansión, fundación o
  fusión. Un cuello de botella debe reducir el N objetivo de la generación
  siguiente y una expansión debe aumentarlo. El historial debe
  distinguir valores heredados y cambios, así como la generación en la que cada
  valor determina la descendencia.
- **REQ-REP-11 (retirado).** Configuración de Ne objetivo absoluto; capacidad
  aplazada fuera de v1.
- **REQ-REP-12 (retirado).** Configuración de un cociente Ne/N objetivo;
  capacidad aplazada fuera de v1.
- **REQ-REP-13.** El sistema no debe intentar alcanzar un Ne objetivo ni emitir
  una desviación respecto a él. Debe calcular el Ne realizado que emerja de la
  reproducción cuando sea evaluable o registrar «no evaluable» y la causa.
- **REQ-REP-14.** El sistema debe calcular Ne realizado exclusivamente como Ne
  reproductivo, a partir de los
  progenitores que dejaron descendencia y de la desigualdad o varianza de sus
  contribuciones. No debe usar estimadores basados en pérdida de
  heterocigosidad, incremento de endogamia o desequilibrio de ligamiento. La
  fórmula concreta se elegirá durante la planificación sin cambiar este
  significado funcional.
- **REQ-REP-15 (retirado).** Desviación y advertencia entre Ne objetivo y Ne
  realizado; no existen en v1 porque Ne objetivo no se admite.
- **REQ-REP-16.** El historial y el manifiesto deben registrar por generación y
  subpoblación los parentales, descendientes, tasa de autofecundación objetivo,
  número esperado y realizado de autofecundaciones, N objetivo, número de
  intentos, éxitos y fallos por causa, N realizado, su desviación y aviso, y Ne
  realizado o su estado «no evaluable» con causa. Debe registrar por separado
  las exclusiones y los fallos debidos a ploidía o estructura reproductiva
  incompatible, sin etiquetarlos como selección de viabilidad.
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
  resumido de demografía, procesos, N objetivo, intentos, éxitos, fallos por
  causa, N realizado, Ne realizado o su causa de no evaluación,
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
- **REQ-REP-24.** El orden obligatorio y único de cada generación debe ser:
  (1) eventos demográficos y determinación de N objetivo; (2) migración;
  (3) selección de viabilidad; (4) elección de parentales y reproducción;
  (5) segregación y recombinación; (6) mutaciones de novo en los haplotipos
  descendientes; (7) cálculo y registro de N realizado, Ne realizado y
  frecuencias; y (8) muestreos y salidas solicitadas. La deriva debe emerger
  durante la reproducción. Una etapa opcional desactivada o no aplicable debe
  registrarse como omitida sin cambiar la posición de las demás.
  Cuando haya varios eventos demográficos, su orden interno no debe ser
  implícito: debe proceder de la secuencia explícita exigida por REQ-CONF-105.
  Este orden interno no altera las ocho etapas.
- **REQ-REP-25.** La configuración no debe admitir órdenes alternativos ni un
  grafo configurable de dependencias entre las ocho etapas. En particular, la
  mutación de novo debe ocurrir después de la segregación y recombinación y antes
  del cálculo y registro y de los muestreos y salidas.
- **REQ-REP-26.** Si la migración, la selección u otro resultado aleatorio de la
  ejecución o del propio apareamiento impiden completar N con éxito, el sistema
  debe conservar los descendientes creados. Debe procesar los N intentos o
  registrar inequívocamente que los restantes son imposibles por la misma
  ausencia de opciones válidas. Por intento o causa agregable inequívoca debe
  registrar éxito o fallo y su causa, además de N objetivo, total de intentos,
  N realizado, desviación y advertencia no invalidante. Una imposibilidad
  estructural detectable antes de ejecutar debe rechazarse entonces. Si la
  incompatibilidad solo aparece en el estado realizado de una generación, los
  parentales o parejas incompatibles no deben considerarse opciones válidas y
  la ausencia de alternativas debe registrarse como causa del fallo. Si N
  realizado es 0, debe aplicarse REQ-SEL-09.

#### 4.3.6. Selección de viabilidad

Un carácter cuantitativo es una característica expresada mediante un valor
numérico. En esta versión depende de la dosis alternativa de un único SNP
bialélico; no convierte el modelo en selección poligénica. Una función
direccional asigna mayor aptitud de forma sistemática
hacia uno de los dos sentidos del carácter: valores mayores o valores menores.

- **REQ-SEL-01.** La primera versión debe usar un único locus seleccionado para
  toda la ejecución. Configurar loci distintos entre fases, intervalos,
  poblaciones o subpoblaciones debe ser un error.
- **REQ-SEL-02.** El locus seleccionado debe ser exclusivamente un SNP
  bialélico. V1 debe rechazar como locus de selección un *indel*, una repetición,
  una inserción, una variante estructural, un locus multialélico o un locus con
  número variable de copias.
- **REQ-SEL-03.** La configuración debe identificar el locus mediante referencia
  estructural, cromosoma o *scaffold*, coordenada e identificador estable. Debe
  existir exactamente una vez en cada copia haplotípica aplicable, para ploidías
  de 1 a 10, y su base de referencia debe ser A, C, G o T, nunca `N`.
- **REQ-SEL-04.** El genotipo seleccionado debe expresarse como la dosis del
  alelo alternativo, desde 0 hasta la ploidía aplicable. La configuración debe
  asignar un valor genético numérico inequívoco a cada dosis posible de cada
  ámbito o ploidía que pueda aparecer.
- **REQ-SEL-05.** Debe rechazarse antes de simular cualquier ámbito con una dosis
  posible sin valor genético explícito. V1 no debe completar valores mediante
  una regla implícita ni extrapolar desde otras dosis.
- **REQ-SEL-06.** Los valores genéticos y la función de aptitud del mismo locus
  pueden variar por fase, intervalo y subpoblación. El valor más específico debe
  prevalecer y el locus no debe cambiar.
- **REQ-SEL-07.** La supervivencia debe decidirse de forma estocástica y
  reproducible para cada individuo. La reproducción posterior debe elegir
  uniformemente entre los supervivientes elegibles.
- **REQ-SEL-08.** Deben registrarse la dosis, los valores genético y del
  carácter, las aptitudes antes y después de normalizar, las decisiones de
  supervivencia y los supervivientes por dosis y generación.
- **REQ-SEL-09.** Si quedan cero individuos o menos parentales elegibles de los
  necesarios tras una selección válida, el sistema no debe declarar
  automáticamente una extinción. Debe aplicar REQ-REP-26 y, si N realizado es
  0, la política de extinción por resultado aplicable. Un fallo de identidad,
  dosis o aptitud del locus no debe interpretarse como extinción.
- **REQ-SEL-10.** Tras una extinción configurada explícitamente, N deja de ser
  aplicable. Tras una extinción por resultado autorizada, los eventos
  dependientes deben omitirse conforme a REQ-EVO-72. Una fundación o
  recolonización posterior válida debe crear otra subpoblación con un
  identificador nuevo; no debe reactivar la extinguida.
- **REQ-SEL-11.** El locus debe quedar protegido durante toda la ejecución.
  Ninguna mutación, *indel*, repetición, deleción, duplicación, inversión,
  translocación u otro cambio puede eliminarlo, copiarlo, moverlo, solaparlo o
  sustituir su identidad. La recombinación puede transmitirlo cuando conserva
  su identidad única.
- **REQ-SEL-12.** El sistema debe mantener por separado dosis, valor genético,
  valor del carácter y aptitud de cada individuo. Toda aptitud calculada debe
  ser finita y no negativa.
- **REQ-SEL-13.** Inmediatamente antes de decidir supervivencia, cada individuo
  elegible debe tener locus único, dosis, valor genético, carácter y aptitud
  inequívocos, y al menos uno debe tener aptitud estrictamente positiva. Si el
  máximo es 0 o cualquier locus, dosis o aptitud es ausente, múltiple, ambiguo,
  negativo o no finito, el sistema debe detenerse antes de seleccionar o
  reproducir. El diagnóstico debe identificar generación, población o
  subpoblación, locus y genotipos o dosis afectados. No debe sustituirse por
  supervivencia uniforme ni presentarse como extinción.
- **REQ-SEL-14.** En la primera versión, el valor del carácter debe coincidir
  exactamente con el valor genético asignado a la dosis. La configuración no
  debe admitir ruido ambiental, efectos ambientales ni heredabilidad.
- **REQ-SEL-15.** La función de aptitud debe ser direccional y la configuración
  debe elegir exactamente un sentido: favorecer valores mayores o favorecer
  valores menores del carácter. Al favorecer valores mayores, un valor mayor no
  debe recibir menor aptitud que uno inferior; en el sentido contrario debe
  cumplirse la relación inversa.
- **REQ-SEL-16.** La primera versión debe rechazar funciones estabilizadoras o
  disruptivas. No debe aproximarlas silenciosamente mediante una función
  direccional.
- **REQ-SEL-17.** Todo evento explícito que pueda afectar al locus protegido debe
  rechazarse antes de simular e identificar la operación incompatible. El locus
  y toda operación que pueda volverlo ambiguo deben quedar excluidos de los
  espacios válidos de eventos aleatorios conforme a REQ-IND-54 a REQ-IND-57 y
  REQ-MUT-13 a REQ-MUT-15.

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
- **REQ-MIG-09.** Después de aplicar conjuntamente las rutas migratorias y antes
  de reproducir, todos los individuos del conjunto reproductor de cada
  subpoblación deben usar la misma referencia estructural. Antes de simular, el
  sistema debe rechazar cualquier combinación de rutas que pueda incumplir esta
  regla e identificar las rutas, las subpoblaciones, los individuos o grupos y
  las referencias afectados. Las subpoblaciones basadas en referencias
  estructurales distintas pueden coexistir en la primera versión, pero deben
  permanecer reproductivamente aisladas y no intercambiar migrantes.

#### 4.3.8. Mutaciones evolutivas

- **REQ-MUT-01.** Las tasas evolutivas de SNP, inserciones, deleciones y de cada
  clase de variante estructural deben expresarse por oportunidad válida y
  generación, por separado, con valores globales y sustituciones por cromosoma
  o *scaffold* y región. Para cada clase y ámbito efectivo, después de
  segregación y recombinación y justo antes de procesar esa clase, el número de
  oportunidades debe contarse sobre todas las copias haplotípicas descendientes
  en su estado actual. La media de Poisson debe ser la tasa efectiva multiplicada
  por ese número.
- **REQ-MUT-02.** Cada mutación nueva debe aplicarse a una copia o haplotipo
  concreto ya segregado y recombinado de un individuo de la nueva generación.
  No debe modificar retrospectivamente al parental.
- **REQ-MUT-03.** Cada evento debe conservar generación, individuo, referencia,
  posición, haplotipo, tipo, origen de novo y transmisión posterior.
- **REQ-MUT-04.** El modelo de SNP debe permitir un espectro o matriz válida con
  las probabilidades de cada sustitución, incluida la relación entre
  transiciones y transversiones. El espectro solo debe usar A, C, G y T como
  bases de origen y destino y nunca debe producir `N`.
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
  mismo genoma o una secuencia repetida. Copiar una región preexistente debe
  clasificarse como duplicación segmentaria genérica; ninguna de estas
  operaciones debe atribuir semántica de elemento transponible. La secuencia
  nueva aleatoria solo debe contener A, C, G y T; una secuencia concreta o
  copiada puede conservar `N`.
- **REQ-MUT-11.** Las repeticiones en tándem deben formar una clase mutacional
  independiente, con tasa por locus y generación o reglas globales o regionales,
  junto con un sesgo configurable entre expansión y contracción. Cada evento
  debe cambiar exactamente una unidad completa. La pureza y las interrupciones
  opcionales se rigen por sus requisitos propios y no cambian el tamaño de este
  paso. Un locus cuyo motivo contenga `N` no es elegible ni cuenta como
  oportunidad.
- **REQ-MUT-12.** El historial, el manifiesto y la verdad deben conservar por
  clase y ámbito la tasa efectiva, el número de oportunidades válidas, la media
  de Poisson, el número sorteado, el número aceptado antes de un posible fallo,
  el número finalmente activo y sus transmisiones posteriores.
- **REQ-MUT-13.** Cada copia haplotípica actual debe aportar sus propias
  oportunidades: una región duplicada cuenta una vez por cada copia actual y una
  región eliminada no cuenta. SNP, inserciones, deleciones y cada clase
  estructural deben usar sus propios candidatos válidos. En una clase que
  necesite fuente y destino, varios extremos o un intervalo, solo el ancla, punto
  o intervalo inicial válido determina el número de oportunidades; los demás
  elementos se eligen condicionados entre destinos válidos y no multiplican el
  denominador. Cada locus de repetición actual y elegible de cada copia es una
  oportunidad.
- **REQ-MUT-14.** Las precedencias de ámbitos deben producir espacios efectivos
  disjuntos conforme a REQ-CONF-116 a REQ-CONF-121 para no contar dos veces una oportunidad situada en regiones
  solapadas. El recuento se sortea una sola vez al inicio de cada clase. Después,
  cada evento actualiza el espacio válido conforme a REQ-IND-54 a REQ-IND-57,
  sin volver a sortear el recuento por los cambios ocurridos dentro de esa clase.
- **REQ-MUT-15.** Una oportunidad solo debe contarse si es válida según todas
  las reglas de su clase. Una posición con `N` no cuenta para SNP ni para una
  operación cuyo ancla exija conocer una base concreta. Su presencia dentro de
  un intervalo no invalida una operación puramente estructural que pueda
  conservarla sin interpretarla.

#### 4.3.9. Requisitos de elementos transponibles retirados de v1

Los identificadores siguientes se conservan únicamente para mantener la
trazabilidad con revisiones anteriores. No expresan capacidades funcionales de
v1 y no deben aceptarse como campos de configuración ni producir metadatos o
salidas específicos.

- **REQ-TE-01 (retirado).** Catálogo de TE suministrado o simulado.
- **REQ-TE-02 (retirado).** FASTA externo o generado de TE.
- **REQ-TE-03 (retirado).** Inserción específica de TE en referencias o individuos.
- **REQ-TE-04 (retirado).** Propiedades específicas de una inserción TE.
- **REQ-TE-05 (retirado).** Selección de destino propia de TE.
- **REQ-TE-06 (retirado).** Validación y registro del destino TE.
- **REQ-TE-07 (retirado).** Orientación, truncamiento y sitio diana de TE.
- **REQ-TE-08 (retirado).** Secuencia fuente y final específica de TE.
- **REQ-TE-09 (retirado).** Mecanismos de copiar o cortar y pegar TE.
- **REQ-TE-10 (retirado).** Reglas entre copias donante y nueva de TE.
- **REQ-TE-11 (retirado).** Regiones y motivos diana por familia TE.
- **REQ-TE-12 (retirado).** Actualización y ausencia de destinos TE.
- **REQ-TE-13 (retirado).** Identidad y metadatos propios de copias TE.
- **REQ-TE-14 (retirado).** Representación específica de TE en verdad o VCF.

#### 4.3.10. Recombinación de la primera versión

- **REQ-REC-01.** Toda la ejecución debe usar exclusivamente una tasa uniforme
  para todo el genoma o mapas completos de recombinación, conforme a
  REQ-CONF-46 a REQ-CONF-50.
- **REQ-REC-02.** La primera versión no debe modelar interferencia entre
  entrecruzamientos. Sus posiciones deben generarse de forma independiente según
  la tasa uniforme o el mapa completo y la semilla.
- **REQ-REC-03.** La recombinación solo debe comenzar cuando todas las copias
  haplotípicas reproductivas de ambos parentales tengan el mismo estado
  estructural canónico conforme a REQ-CRU-30. V1 no debe recombinar intervalos
  colineales de copias estructuralmente diferentes ni tratar esa compatibilidad
  parcial como suficiente.
- **REQ-REC-04.** Una diferencia estructural debe invalidar al individuo o a la
  pareja como opción reproductiva completa, no solo excluir una región del
  espacio recombinante. V1 no debe generar gametos desequilibrados ni modelar
  emparejamientos complejos, fertilidad reducida o compatibilidad parcial.
- **REQ-REC-05.** Para cada cromosoma o *scaffold* compatible de ploidía par 2,
  4, 6, 8 o 10, cada meiosis debe formar la mitad de bivalentes que de copias.
  Cada bivalente debe recombinar conforme a la tasa o mapa cuando corresponda y
  aportar una copia reconstruida al gameto.
- **REQ-REC-06.** Como todas las copias de una secuencia reproductiva tienen la
  misma estructura, deben emparejarse al azar para formar bivalentes, de manera
  reproducible mediante la semilla. Los emparejamientos, los estados
  estructurales comparados y los haplotipos transmitidos deben quedar trazados.
- **REQ-REC-07.** La primera versión no debe aplicar preferencia de
  apareamiento por subgenoma, emparejamiento homeólogo preferente, multivalentes
  ni doble reducción.
- **REQ-REC-08.** Cada parental debe aportar la mitad de la ploidía al
  descendiente. Los *scaffolds* deben seguir la misma segregación. Cuando a un
  *scaffold* no le corresponda recombinación por su umbral efectivo, la copia
  elegida de cada bivalente debe transmitirse intacta. El método técnico para
  generar las posiciones se concretará por separado sin alterar estas reglas.
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
  los individuos incluidos en muestreos o materializaciones adicionales
  solicitados expresamente. Para cada uno debe indicar los puntos de
  recombinación y el origen parental de cada segmento. Las demás generaciones
  deben conservar resúmenes, no una exportación individual de todos sus
  segmentos.
- **REQ-REC-12.** Para cada segmento de un individuo incluido en la salida de
  REQ-REC-11, el origen debe
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
  cada copia del mismo estado estructural canónico debe tener la misma probabilidad de
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
(*paired-end*) se observan ambos extremos del mismo fragmento. Illumina es el
ejemplo principal de lecturas cortas admitidas en la primera versión. Las
lecturas largas, incluido PacBio HiFi, quedan aplazadas en DTF-049.

F4 sigue este orden normativo: calcula las parejas de cobertura, sortea el
fragmento biológico, forma una quimera cuando corresponda, realiza calidad y
errores de R1/R2, aplica contaminación biológica y finalmente demultiplexa.
El demultiplexado puede transferir una pareja ya contaminante. El manifiesto conserva solo estadísticas agregadas de
procedencia, nunca un registro por lectura ni pareja. El modo directo de
alineamientos sintéticos queda aplazado en DTF-051.
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
  referencia de mapeo, la muestra, el evento, la biblioteca, la tanda y
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
- **REQ-DAT-12.** F4 debe simular exclusivamente lecturas cortas Illumina
  *paired-end*. Debe permitir configurar longitud, cobertura, errores y calidad,
  y usa orientación `FR`: R1 procede del extremo izquierdo en orientación directa
  y R2 del derecho como complemento inverso, ambas hacia el interior del inserto.
  La distribución de tamaños de inserto pertenece al perfil efectivo.
  Si los dos mates cubren una región común, ambos se emiten completos e
  independientes; F4 no los fusiona ni elimina su solapamiento.
- **REQ-DAT-13.** El sistema debe rechazar en esta versión una solicitud de
  lecturas largas PacBio HiFi e indicar DTF-049 como destino de la ampliación.
- **REQ-DAT-14.** F4 ofrece exclusivamente el perfil predeterminado validado
  Illumina *paired-end* y sus ajustes permitidos. Illumina *single-end* queda
  aplazado en DTF-057; no se selecciona ni se infiere un perfil SE alternativo.
- **REQ-DAT-15.** En el modo realista, las lecturas deben generarse desde los
  haplotipos faseados de cada individuo y usar las decisiones aleatorias
  reproducibles de la ejecución. Cada decisión se deriva de un ámbito RNG local
  identificado por el evento y las identidades estables pertinentes; no puede
  depender del orden físico, de un hilo ni de un flujo RNG compartido.
- **REQ-DAT-16.** Cada FASTQ, BAM o CRAM debe poder relacionarse con la
  tecnología, el perfil y los parámetros efectivos con los que se generaron sus
  lecturas o alineamientos.
- **REQ-DAT-17.** En el modo realista, el sistema debe guardar como salida los
  FASTQ finales después de formar el fragmento, la quimera opcional y sus
  calidades/errores, y de aplicar después contaminación biológica y errores de
  demultiplexado. Una reasignación puede afectar a una pareja ya contaminante.
  F4 no genera BAM ni CRAM; fases posteriores podrán usar esas mismas lecturas.
  No debe reconstruir esos FASTQ a partir de los archivos de alineamiento.
  En F4, cada pareja final se publica como
  `outputs/fastq/<event_id>/<individual_id>_R1.fastq.gz` y
  `outputs/fastq/<event_id>/<individual_id>_R2.fastq.gz`: gzip en el perfil BGZF
  v1, con DEFLATE nivel 6, bloques fijos y `MTIME = 0`. Sus tamaños y SHA-256 se
  calculan sobre los bytes comprimidos finales; la comprobación de formato, IDs,
  secuencias y calidades relee por flujo ambos ficheros descomprimidos antes de
  la publicación atómica. Esta fijación solo hace repetibles los bytes bajo el
  mismo simulador, lock y plataforma admitida; no amplía la garantía general de
  identidad binaria de REQ-REP-23/DTC-038.
- **REQ-DAT-18.** Los FASTQ del modo realista deben conservar la secuencia, la
  calidad, el identificador de cada lectura y su emparejamiento cuando exista.
  El sistema debe conservar los mismos identificadores y pares al mapearlos y no
  debe descartar las lecturas no mapeadas que produzca el mapeador. Cada lectura
  *single-end* o larga, y cada pareja *paired-end* considerada como una unidad,
  debe tener una identidad técnica única dentro del evento.
- **REQ-DAT-19.** El sistema debe rechazar en esta versión el modo directo de
  BAM o CRAM sintéticos e indicar DTF-051 como destino de la ampliación.
- **REQ-DAT-20.** Los FASTQ, BAM y CRAM del modo realista deben cumplir los
  estándares de sus formatos y la versión oficial que se elija durante la
  planificación. Todo BAM o CRAM que se incorpore desde F5.6 debe proceder de
  un FASTQ final.
- **REQ-DAT-21.** Cada evento con alineamientos debe usar exactamente un mapeador. Si
  la configuración no elige otro compatible, debe usar minimap2.
- **REQ-DAT-22.** Si un evento selecciona varias referencias de mapeo, debe
  aplicar el mismo mapeador y los mismos FASTQ a cada referencia y producir
  conjuntos independientes. El mapeador puede ser distinto entre eventos, pero
  cada evento debe generar sus FASTQ propios y no reutilizarlos para comparar
  mapeadores.
- **REQ-DAT-23.** Cuando se use minimap2, el perfil de lecturas debe seleccionar
  automáticamente el *preset* y los parámetros predeterminados apropiados para
  lecturas cortas tipo Illumina. Las sustituciones
  de parámetros deben prevalecer, generar una advertencia y validarse para
  detectar opciones desconocidas o incompatibles cuando sea posible.
- **REQ-DAT-24.** En F4, la cobertura realizada debe ser estocástica y
  aproximarse a la única cobertura media objetivo del evento durante el muestreo
  WGS básico, antes de sesgos y artefactos. En el
  modo realista debe aplicarse al generar las lecturas. Las regiones se interpretan en la
  referencia estructural del individuo. No debe forzarse una profundidad exacta
  en cada base. La cobertura final del FASTQ y la cobertura alineada observada
  contra cada referencia deben medirse después y no constituyen nuevos objetivos.
  Cada fragmento PE se elige directamente entre los inicios que permiten alojar
  el inserto completo en el mismo cromosoma o *scaffold*; no se recorta ni se
  reintenta un extremo inválido.
- **REQ-DAT-25.** En el modo realista, las lecturas contaminantes de todas las
  fuentes biológicas configuradas deben sustituir lecturas propias de la muestra
  antes del demultiplexado. Para cada receptora, las parejas propias retiradas se
  eligen uniformemente y sin reemplazo. La pareja incorporada recibe un
  identificador técnico nuevo de la biblioteca receptora, pero conserva la
  secuencia y calidades ya realizadas para el donante; no vuelve a sortear error
  ni calidad. El mapeo queda aplazado
  en DTF-055; un mapeador posterior debe procesar el resultado completo contra
  cada referencia seleccionada.
- **REQ-DAT-26.** Durante la generación, el sistema puede usar la identidad de
  las parejas contaminantes y de origen para aplicar la sustitución correcta,
  pero no debe conservarla ni publicarla una vez obtenido el FASTQ final. El
  manifiesto informa solo los aportes agregados por donante y receptora conforme
  a REQ-DAT-34–36.
- **REQ-DAT-27.** Cuando se contabilice la contaminación realizada en datos
  *paired-end*, cada unidad debe incorporar juntos los dos
  *mates* de la misma fuente y conservar nombres, indicadores del formato,
  orientación y relación de pareja coherentes. Esta operación no debe crear
  pares huérfanos.
- **REQ-DAT-28.** F4 debe rechazar los FASTQ externos como fuente de
  contaminación e indicar DTF-052. La validación, copia y submuestreo de esas
  entradas pertenecen a la ampliación posterior y no deben crear artefactos ni
  directorios de resultados en este hito.
- **REQ-DAT-29.** El sistema debe aplicar y registrar por separado la
  contaminación biológica y la técnica. La contaminación técnica debe ocurrir
  entre las muestras y bibliotecas procesadas conjuntamente en la única tanda
  del evento.
- **REQ-DAT-30.** La contaminación técnica debe poder producir dos resultados:
  parejas completas transferidas de una muestra a otra, que desaparecen de la
  donante sin duplicarse, y parejas que quedan sin asignar o son filtradas. Los
  recuentos realizados deben ser reproducibles.
- **REQ-DAT-31.** El uso obligatorio de UDI debe permitir reconocer
  combinaciones de índices inesperadas y simular su filtrado defensivo durante
  el demultiplexado. No debe asumirse que los UDI impiden físicamente el *index
  hopping*.
- **REQ-DAT-32.** El sistema debe simular directamente el efecto final del
  demultiplexado sobre cada muestra: debe mover una pareja afectada a otra
  muestra, o descartarla o dejarla sin asignar según el esquema de índices. Solo
  conserva recuentos agregados por mecanismo, origen/destino cuando corresponda
  y causa.
- **REQ-DAT-33.** En el modo realista, el sistema debe conservar un único FASTQ
  final, o una pareja de FASTQ para datos *paired-end*, por individuo y evento.
  Este archivo corresponde a la única biblioteca de ese individuo y a la única
  tanda del evento. No debe guardar además una copia limpia anterior a la
  contaminación o a los errores de demultiplexado.
- **REQ-DAT-34.** Por cada FASTQ final, el manifiesto o resumen de evento debe
  conservar estadísticas agregadas de procedencia: semilla, perfil efectivo,
  cobertura objetivo/real, contaminación objetivo/real y aportes por donante,
  quimeras, duplicados y descartes por causa. El perfil efectivo incluye los
  desplazamientos de calidad y escalas de error de R1 y R2, incluso cuando usan
  sus valores neutros predeterminados.
- **REQ-DAT-35.** F4 no debe generar un fichero de procedencia ni identificar el
  origen de cada lectura o pareja. El resumen no duplica secuencias, calidades ni
  identificadores técnicos individuales.
- **REQ-DAT-36.** Las lecturas contaminantes, transferidas, descartadas o sin
  asignar se informan exclusivamente como recuentos agregados por receptora,
  donante cuando corresponda y causa; no se enumeran identificadores.
- **REQ-DAT-37.** El sistema debe aplicar el perfil de errores efectivo de la
  tecnología y su versión al generar las lecturas. Debe conservar el perfil
  predeterminado y los ajustes globales por mate realmente aplicados. En el
  perfil inicial, cuando ocurre un error, la base observada se
  elige uniformemente entre las otras tres bases A/C/G/T. Una matriz empírica de
  sustituciones queda aplazada en DTF-053. En v1, el modelo de error no debe
  resolver una `N` de la plantilla ni generar nuevas `N` como error o llamada
  ausente. Para longitudes inferiores a 150 pb usa el prefijo del recurso propio
  `resources/f4/gps-hiseqx-pcrfree-pe150-v1.json`, SHA-256
  `ae8d598a95f7a9826b3a8c563a023d745f85031651c1e8c03ee9f54310162d3b`; no
  extrapola ciclos por encima de 150 pb. Es una aproximación sintética explícita
  y ART `HSXn` solo una referencia comparativa externa, sin incorporar ni
  redistribuir sus tablas. Si se configura un
  desplazamiento de calidad o una escala de error, se aplica primero el
  desplazamiento a la calidad del perfil, se deriva su probabilidad de error y
  se aplica después la escala; una combinación que produzca una probabilidad
  fuera del intervalo válido se rechaza, sin recorte silencioso.
- **REQ-DAT-38.** Las puntuaciones de calidad por base del FASTQ deben generarse
  de forma coherente con el modelo de errores efectivo: una base de menor
  calidad debe tener mayor probabilidad de ser errónea conforme a la
  interpretación probabilística de la puntuación. R1 y R2 sortean calidad y
  error de forma independiente con sus curvas respectivas, incluso en posiciones
  del inserto que ambas cubren. Tras cualquier ajuste configurado, la
  probabilidad final determina tanto el sorteo de error como la calidad FASTQ
  final; no se emite una calidad del perfil original incompatible con el error
  efectivamente simulado. F4 codifica FASTQ Sanger (`Phred+33`): redondea la
  calidad final al entero más próximo y exige que esté entre 0 y 93. Un ajuste
  que produzca una calidad fuera de ese intervalo se rechaza antes de generar,
  sin saturarla silenciosamente. La `N` de plantilla sigue la excepción de
  REQ-DAT-72 y usa Phred 0.
- **REQ-DAT-39.** Cuando el modelo de la tecnología contemple dependencias de la
  posición dentro de la lectura o del contexto de bases, los errores y las
  puntuaciones de calidad deben reflejar conjuntamente esos patrones.
- **REQ-DAT-40.** El mapeo, cuando se implemente fuera de F4 en DTF-055, debe
  propagar correctamente a BAM o CRAM las puntuaciones de calidad del FASTQ
  final. F4 no conserva diferencias ni procedencia por lectura o pareja: los
  errores, descartes y reasignaciones se describen únicamente mediante los
  recuentos agregados exigidos por REQ-DAT-34–36.
- **REQ-DAT-41 (retirado).** Estimación de perfiles WGS desde un FASTQ real y su
  referencia; capacidad aplazada fuera de v1.
- **REQ-DAT-42 (retirado).** Separación automática de errores de secuenciación,
  variantes reales, regiones ambiguas y errores de mapeo durante la calibración;
  capacidad aplazada fuera de v1.
- **REQ-DAT-43 (retirado).** Conversión de un perfil empírico estimado en perfil
  efectivo; capacidad aplazada fuera de v1.
- **REQ-DAT-44 (retirado).** Generación de un fichero estructurado de perfil
  empírico; artefacto excluido de v1.
- **REQ-DAT-45 (retirado).** Reutilización de un perfil empírico sin sus FASTQ
  originales; capacidad aplazada fuera de v1.
- **REQ-DAT-46.** F4 no genera duplicados: PCR no aplica al perfil inicial y la
  tasa óptica es cero. Los modelos de duplicación quedan aplazados en DTF-056.
- **REQ-DAT-47.** La duplicación de parejas completas queda aplazada en DTF-056.
- **REQ-DAT-48.** La procedencia individual de duplicados queda aplazada en
  DTF-056.
- **REQ-DAT-49.** El marcado de duplicados en BAM o CRAM queda aplazado en
  DTF-055 y DTF-056.
- **REQ-DAT-50.** F4 no añade restos de adaptador: el perfil inicial no contiene
  insertos menores que la longitud de lectura y rechaza su activación. La
  generación derivada de fragmentos cortos queda en DTF-054.
- **REQ-DAT-51.** La identificación de intervalos de adaptador por lectura queda
  aplazada en DTF-054; F4 solo informa recuento cero en su resumen agregado.
- **REQ-DAT-52.** F4 no aplica sesgo GC: toda pareja elegible conserva la
  probabilidad del muestreo WGS básico. La relación no neutra queda aplazada en
  DTF-056.
- **REQ-DAT-53.** La combinación reproducible de sesgo GC con cobertura y sus
  métricas queda aplazada en DTF-056.
- **REQ-DAT-54.** En F4, las lecturas de regiones de baja mapeabilidad deben
  generarse desde sus posiciones y haplotipos verdaderos, sin publicar esa
  procedencia por lectura. El mapeo y sus alineamientos observados quedan
  aplazados en DTF-055.
- **REQ-DAT-55.** La MAPQ y cualquier evaluación de baja mapeabilidad quedan
  aplazadas en DTF-055; cuando exista mapeo, el sistema no debe sobrescribir la
  MAPQ producida por el mapeador.
- **REQ-DAT-56.** F4 no conserva posición ni haplotipo de origen por lectura.
  La futura evaluación de ambigüedad y mapeos en posiciones incorrectas se
  definirá en DTF-055 sin publicar una verdad técnica individual persistente.
- **REQ-DAT-57.** La clasificación sintética de lecturas para el modo directo
  queda aplazada en DTF-051.
- **REQ-DAT-58.** Los campos de registros sintéticos directos quedan aplazados
  en DTF-051.
- **REQ-DAT-59.** La procedencia específica de categorías de alineamiento
  directo queda aplazada en DTF-051.
- **REQ-DAT-60.** Durante la generación, una quimera debe conservar
  transitoriamente los dos fragmentos, sus coordenadas, orientaciones y punto
  de unión para construir la pareja. F4 no publica ni conserva procedencia por
  lectura; informa solo la tasa efectiva y el recuento agregado de quimeras.
- **REQ-DAT-61.** En F4, las parejas quiméricas se incluyen en el FASTQ final.
  BAM, CRAM y el mapeo de sus alineamientos quedan aplazados en DTF-055.
- **REQ-DAT-62.** Los alineamientos sintéticos y el mapeo real de quimeras quedan
  aplazados en DTF-051 y DTF-055, respectivamente.
- **REQ-DAT-63.** Todos los fragmentos de una lectura quimérica deben proceder
  del mismo individuo y de la misma biblioteca. Los intercambios entre muestras
  no deben representarse como quimeras, sino mediante los modelos de
  contaminación o errores de demultiplexado.
- **REQ-DAT-64.** Dentro del mismo individuo y biblioteca, una lectura quimérica
  une exactamente dos segmentos no contiguos. El primero es un prefijo de la
  plantilla y el segundo un intervalo elegido uniformemente entre los que no
  solapan ni tocan ese prefijo cuando comparten copia; sus longitudes suman el
  inserto y se unen en un único punto interno. Sus fuentes pueden pertenecer a
  cromosomas, *scaffolds*, haplotipos o copias diferentes; se leen en
  orientación directa antes de aplicar la disposición PE habitual. Sus datos de
  origen son transitorios y no forman una salida de F4.
- **REQ-DAT-65.** Las anomalías *paired-end* observables mediante alineamiento
  —orientación anómala, tamaño de inserto observado inesperado o mate no
  alineado— quedan aplazadas en DTF-055. F4 solo emite parejas normales según
  su orientación e inserto efectivos, salvo la quimera definida en REQ-DAT-60–64.
- **REQ-DAT-66.** La generación y validación de estados PE observados por un
  mapeador quedan aplazadas en DTF-055.
- **REQ-DAT-67.** Las parejas sintéticas directas quedan aplazadas en DTF-051;
  no sustituyen el modelo de anomalías PE reservado en DTF-055.
- **REQ-DAT-68.** Una instancia de lectura debe ser una observación concreta
  presente en una biblioteca. Su identificador técnico debe ser único dentro de
  su evento. La cabecera FASTQ es exactamente
  `@GPS:<event_id>:<library_id>:<technical_pair_id> <mate>:N:0:<index>`.
  Los dos miembros de una pareja comparten el prefijo hasta el espacio y usan
  `1` o `2` como `<mate>`. `GPS` declara datos simulados; no se simulan campos
  de instrumento, carril, *tile*, coordenadas, BCL ni control de calidad. En F4,
  `<index>` es siempre el UDI literal `i7+i5`, con ambos componentes de ocho
  bases `ACGT` en orientación de secuenciación; el marcador `0` no es válido.
  Ese nombre de pareja no debe identificar otra pareja del mismo evento.
- **REQ-DAT-69.** La contaminación biológica debe crear en la receptora una
  instancia nueva con identificador técnico nuevo de esa biblioteca y conservar
  la original en la donante. Una misma unidad de
  origen no debe copiarse dos veces dentro de una misma biblioteca receptora;
  la copia conserva secuencia y calidades realizadas y no consume sorteos de
  error o calidad adicionales.
  puede originar una instancia en otra receptora compatible. Una duplicación
  adicional dentro de la receptora solo debe aparecer mediante un mecanismo de
  duplicación configurado. El manifiesto informa los recuentos agregados, no los
  identificadores ni vínculos individuales.
- **REQ-DAT-70.** Una reasignación técnica por demultiplexado debe transferir el
  contenido de la misma pareja: debe retirarla de la donante e incorporarla a la
  receptora sin duplicarla. En la biblioteca final recibe un identificador
  técnico local nuevo; no se publica el identificador previo. El manifiesto
  distingue ambos mecanismos mediante sus recuentos agregados.
- **REQ-DAT-71.** El sistema debe validar antes de combinar FASTQ o producir el
  BAM o CRAM final que no existan colisiones de identificadores técnicos dentro
  del evento. Una colisión debe ser un error; no debe resolverse renombrando una
  lectura sin actualizar de forma atómica FASTQ, BAM o CRAM cuando exista mapeo
  y manifiesto.
- **REQ-DAT-72.** Las regiones que contienen `N` deben seguir siendo elegibles
  para el muestreo WGS. Cada `N` de la plantilla incluida en una lectura debe
  aparecer como `N` en FASTQ con calidad Phred 0 y propagarse sin sustitución al
  BAM o CRAM correspondiente.
- **REQ-DAT-73.** La verdad de lecturas y los informes deben permitir comprobar
  que ninguna `N` de plantilla se resolvió silenciosamente y que el modelo de
  error no creó nuevas `N`. Deben registrar el número o la proporción de `N` en
  las lecturas por biblioteca y ámbito pertinente.

- **REQ-WGS-01.** En F4, el evento debe declarar exactamente una cobertura media
  objetivo común. Un objetivo por número de lecturas o parejas, la ausencia de
  cobertura o más de un objetivo deben ser errores previos a la generación.
- **REQ-WGS-02.** Fuera de F4, si una fase posterior admite un objetivo por
  número, en datos *paired-end* debe expresarlo en parejas completas. En datos
  *single-end* o de lecturas largas debe expresarlo en lecturas individuales.
  F4 rechaza ese objetivo como entrada y solo registra el número de parejas
  derivado de su cobertura.
- **REQ-WGS-03.** El sistema debe derivar y registrar para cada biblioteca el
  número de parejas que corresponde a la cobertura objetivo común y la
  cobertura equivalente realizada.
  En PE, un objetivo de cobertura se convierte al número entero de parejas más
  próximo; un empate se redondea hacia arriba y se registra la cobertura
  equivalente realizada.
- **REQ-WGS-04.** La cobertura en «x» debe definirse como el número de bases de
  lecturas generado para la biblioteca dividido por la longitud haploide de la
  referencia estructural del individuo. Esta longitud es la suma de cada
  cromosoma o *scaffold* de la referencia una sola vez, sin multiplicarla por la
  ploidía. En PE, cada pareja aporta la suma de las longitudes emitidas por R1 y
  R2 —300 bases para PE150— incluso si los mates se solapan. Las posiciones `N`
  cuentan en esta longitud. Las diferencias de
  longitud realizadas en los haplotipos y las
  referencias seleccionadas para el mapeo no deben cambiar el denominador.
- **REQ-WGS-05.** Antes de sesgos y degradaciones, el modelo WGS predeterminado
  debe elegir de forma uniforme y reproducible las posiciones iniciales de los
  fragmentos o lecturas entre las copias cromosómicas accesibles. Una secuencia
  más corta que el inserto mínimo no es accesible para PE, recibe cobertura cero
  y se registra como tal, sin deformar la distribución ni detener la ejecución.
- **REQ-WGS-06.** Antes de aplicar sesgos, la cobertura local realizada debe ser
  estocástica y aproximadamente compatible con una distribución de Poisson
  alrededor de la media objetivo. No debe imponerse una profundidad exacta por
  base.
- **REQ-WGS-07.** Cuando las secuencias tengan ploidías distintas, el muestreo
  debe repartir las lecturas en proporción a la longitud de cada copia y al
  número de copias. En ausencia de excepciones, cada copia debe recibir
  aproximadamente la misma cobertura y la profundidad total de un cromosoma debe
  ser proporcional a su ploidía.
- **REQ-WGS-08.** El sesgo GC, los límites de secuencia y las demás degradaciones
  deben aplicarse después del muestreo WGS básico de forma documentada. F4 no
  admite regiones con cobertura especial.
- **REQ-WGS-09.** Antes de esos efectos, los haplotipos y alelos deben contribuir
  en proporción a su dosis y ploidía. Por ejemplo, 40× totales en un tetraploide
  uniforme deben producir aproximadamente 10× por copia.
- **REQ-WGS-10.** La primera versión no debe ofrecer parámetros para forzar un
  sesgo de biblioteca por alelo ni la pérdida completa de uno de los alelos.
  Las diferencias realizadas solo deben proceder del muestreo aleatorio, la
  cobertura y los errores, sesgos y resultados de mapeo admitidos.
- **REQ-WGS-11.** El manifiesto y el informe de validación deben registrar la
  referencia estructural y la longitud haploide usadas como denominador y
  distinguir, por biblioteca y ámbito aplicable: la cobertura objetivo; la
  cobertura WGS básica realizada antes de sesgos y artefactos; la cobertura
  final del FASTQ después de sesgos, duplicados, contaminación, demultiplexado y
  descartes, separando bases propias y contaminantes; y la cobertura alineada
  observada contra cada referencia de mapeo. Esta última no debe presentarse
  como otro objetivo. Cuando pueda medirse, también debe registrarse la cobertura
  realizada por cromosoma y copia o haplotipo. Deben registrar además el número
  o la proporción de `N` pertinente para interpretar cada medida.
- **REQ-WGS-12.** Solo debe conservarse el BAM o CRAM final que contiene las
  degradaciones configuradas. No debe guardarse simultáneamente una versión
  limpia del mismo experimento.
- **REQ-WGS-13.** F4 rechaza un control sin degradaciones porque exige otro
  evento. Cuando se admita cardinalidad multievento tras DTF-050, ese control
  reutilizará la misma verdad genotípica, referencias y muestras, tendrá su
  propia subsemilla y enlazará la relación experimental sin mezclar resultados.
- **REQ-WGS-14.** Las posiciones `N` deben ser elegibles para el muestreo y
  contar tanto como posiciones cubiertas como bases de lectura al calcular la
  cobertura. No deben excluirse del numerador ni del denominador; los resúmenes
  de `N` de REQ-WGS-11 permiten interpretar la medida.

- **REQ-CONT-01.** La contaminación biológica debe ser sustitutiva: por cada
  lectura o pareja contaminante incorporada a la receptora debe retirarse una
  lectura o pareja propia. La biblioteca donante no debe perder datos por este
  mecanismo. Toda selección de donantes usa las lecturas propias previas a la
  contaminación; no encadena contaminaciones ni depende del orden de receptoras.
- **REQ-CONT-02.** La configuración declara un intervalo de fracción total entre
  0 % y 100 %, inclusive, y opcionalmente un subconjunto de receptoras. Por
  omisión, ese subconjunto es toda la población seleccionada. Cada receptora
  seleccionada sortea de forma independiente y reproducible un valor uniforme
  dentro del intervalo; el 100 % representa sustitución completa.
- **REQ-CONT-03.** Las contribuciones de las fuentes sorteadas para una receptora
  deben sumar exactamente su objetivo realizado y no superar los datos
  disponibles. La conversión a lecturas o parejas completas se realiza antes de
  seleccionarlas; un objetivo imposible detiene la ejecución.
- **REQ-CONT-04.** Toda fuente debe ser compatible con la tecnología, el perfil
  y la disposición de lecturas de la biblioteca receptora. No deben mezclarse
  lecturas cortas y largas, perfiles tecnológicos incompatibles ni disposiciones
  *single-end* y *paired-end* incompatibles.
- **REQ-CONT-05.** Las fuentes son exclusivamente los otros individuos elegibles
  de la misma población de la receptora; la receptora nunca es su propia fuente.
  Un individuo receptor puede ser donante de otra receptora desde sus lecturas
  propias previas a contaminación.
- **REQ-CONT-06.** F4 no admite una matriz dirigida ni pesos declarados de
  donante→receptora. Esa ampliación queda fuera de este hito.
- **REQ-CONT-07.** Antes de generar deben validarse el rango, los receptores y
  que cada receptor con objetivo positivo tenga donantes elegibles con al menos
  tantas parejas limpias distintas en conjunto como su objetivo convertido. Si
  no hay capacidad suficiente, el evento se rechaza antes de escribir FASTQ; no
  se reutiliza una pareja de origen ni se reduce silenciosamente el objetivo.
- **REQ-CONT-08.** Para cada receptora, sus parejas contaminantes se eligen al
  azar y de forma reproducible, sin reemplazo, del conjunto de todas las parejas
  limpias elegibles de los demás individuos de su misma población. Las
  contribuciones emergen proporcionalmente a esa disponibilidad, no de un
  porcentaje ni de pesos configurados por donante.
- **REQ-CONT-09.** El manifiesto y el informe deben registrar por receptora el
  intervalo configurado, objetivo sorteado, número de parejas limpias, conversión
  redondeada a parejas, donantes elegidos y aportes realizados.
- **REQ-CONT-10.** En errores de índices o demultiplexado, cada lectura o pareja
  reasignada debe eliminarse de la donante e incorporarse a la receptora sin
  duplicarse. El total de la tanda debe conservarse salvo los descartes o datos
  sin asignar. En F4, *index hopping* selecciona esa fracción de parejas de la
  tanda y asigna cada una uniformemente a otra biblioteca; la transferencia se
  aplica después de la contaminación biológica sobre una instantánea inmutable y
  simultáneamente. Cada pareja pre-salto se evalúa una sola vez, por lo que no
  existen cadenas de saltos ni una reasignación dependiente del orden.
- **REQ-CONT-11.** Las lecturas o parejas descartadas o sin asignar durante el
  demultiplexado no deben reemplazarse. Deben reducir el número de lecturas y la
  cobertura final de la muestra; en *paired-end* debe descartarse la pareja
  completa conforme a las reglas acordadas.
- **REQ-CONT-12.** La representación de contaminación mediante alineamientos
  sintéticos directos queda aplazada en DTF-051.
- **REQ-CONT-13.** La comparación entre modo directo y realista queda aplazada
  en DTF-051. F4 solo admite contaminación sobre el FASTQ final realista.
- **REQ-CONT-14.** Por muestra, fuente y mecanismo, el manifiesto y el informe
  deben registrar objetivos, lecturas propias retiradas, contaminantes
  incorporadas, lecturas transferidas, descartadas o sin asignar, totales finales
  y cobertura final realizada.
- **REQ-CONT-15.** En F4, la única unidad de contaminación es el porcentaje total
  por receptora. Tras conocer sus parejas limpias, el número objetivo es el
  entero más próximo a `porcentaje × parejas_limpias`, con empate hacia arriba;
  no se admiten cobertura equivalente, número absoluto ni porcentajes por
  donante.
- **REQ-CONT-16.** Para cada receptora, el sistema selecciona exactamente su
  número objetivo de parejas de los donantes elegibles mediante sorteos
  reproducibles sin reemplazo. No resuelve una matriz ni una optimización global;
  el manifiesto informa los aportes realizados por donante.

Un grupo de lectura reúne las lecturas obtenidas para un individuo dentro de un
evento técnico. Como v1 asigna a ese individuo una sola biblioteca y al evento
una sola tanda, existe exactamente un grupo de lectura por individuo y evento.
En BAM o CRAM se describe mediante una entrada `@RG` de la cabecera y la etiqueta
`RG` de cada registro.

- **REQ-RG-01.** Los BAM o CRAM deben incluir grupos de lectura coherentes que
  identifiquen como mínimo el evento, la muestra o individuo, su única
  biblioteca, la única tanda y la plataforma o tecnología.
- **REQ-RG-02.** Cada individuo seleccionado por un evento técnico debe tener
  exactamente una biblioteca, un grupo de lectura y, en el modo realista, un
  FASTQ final o una pareja de FASTQ. Solicitar otra biblioteca del mismo
  individuo debe requerir otro evento.
- **REQ-RG-03.** Antes de generar lecturas, el sistema debe rechazar las
  colisiones de identificadores de grupos de lectura, bibliotecas o tandas e
  identificar los valores en conflicto. El identificador de grupo de lectura
  debe ser único en toda la ejecución.
- **REQ-RG-04.** Los identificadores de muestra, biblioteca, tanda y
  grupo de lectura deben ser consistentes entre FASTQ, BAM o CRAM cuando exista
  mapeo y manifiesto. El identificador de grupo debe distinguir evento,
  biblioteca y tanda. Los identificadores técnicos de lectura deben ser únicos
  dentro de su evento, sin publicar vínculos de origen.
- **REQ-RG-05.** Dentro de cada conjunto definido por una referencia de mapeo,
  el sistema debe generar un único BAM o CRAM final por individuo y preservar
  su único grupo de lectura del evento. No debe mezclar referencias ni requerir
  como producto principal archivos finales separados por biblioteca; los
  intermedios técnicos quedan fuera de las salidas principales.
- **REQ-RG-06.** Los FASTQ finales deben entregarse separados por individuo y
  evento, sin concatenarlos como producto principal. Cada archivo o pareja debe
  conservar su perfil y parámetros y corresponder al único grupo de lectura de
  ese individuo en el BAM o CRAM final.
- **REQ-RG-07.** El campo de muestra `SM` de BAM o CRAM y la columna de muestra
  del VCF deben usar el identificador biológico estable del individuo. El evento
  debe identificarse mediante el directorio, el manifiesto y el grupo de lectura,
  sin modificar `SM`.

La expresión «otras degradaciones» no amplía por sí sola el alcance de v1. Solo
son obligatorias las degradaciones descritas mediante requisitos concretos en
este documento. Añadir otra requerirá una decisión funcional explícita y una
ampliación de la especificación.

### 4.5. VCF de verdad y coherencia de las salidas

Un caso de verdad es una diferencia biológica individual o compuesta respecto a
la referencia de mapeo de un conjunto. Puede necesitar uno o varios registros
VCF para expresarse. Un punto de ruptura (*breakend*, BND) representa uno de los
extremos de una nueva adyacencia entre segmentos; por ello una nueva adyacencia
suele necesitar dos registros BND relacionados.

- **REQ-VER-01.** Cuando la configuración solicite VCF de verdad, el sistema debe
  generarlo con los genotipos reales simulados de los SNP, *indels*, variantes
  estructurales y repeticiones en tándem de los individuos incluidos.
- **REQ-VER-02.** Cada genotipo del VCF de verdad debe corresponder al individuo
  simulado identificado en el archivo.
- **REQ-VER-03.** Cada VCF de verdad debe estar referido exactamente a la
  referencia de mapeo seleccionada para su conjunto y, cuando existan BAM o CRAM
  correspondientes, al mismo genoma que esos archivos.
- **REQ-VER-04.** Cuando se solicite VCF, el sistema debe generar uno de verdad
  independiente para cada conjunto definido por el evento y una referencia de
  mapeo. No debe generar VCF separados por individuo, clase de variante ni caso
  complejo.
- **REQ-VER-05.** Cada conjunto debe usar exactamente una referencia de mapeo.
  Todos los BAM, CRAM o VCF que existan en él deben estar referidos a esa misma
  referencia. No se deben mezclar datos referidos a genomas distintos.
- **REQ-VER-06.** Cada VCF de verdad debe incluir en su propia cabecera o
  metadatos estándar la información necesaria para identificar sin ambigüedad la
  referencia usada y el conjunto al que pertenece. Debe cumplir estrictamente
  VCF 4.5 conforme a DTC-062 y comenzar por `##fileformat=VCFv4.5`.
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
  estadísticas. La normalización técnica debe seguir el algoritmo interno
  determinista y versionado de DTC-064.
- **REQ-VER-17.** Los SNP y *indels* deben representarse mediante alelos
  explícitos REF y ALT, normalizados respecto a la referencia de mapeo del
  conjunto. REF es el alelo de referencia y ALT el alelo alternativo.
- **REQ-VER-18.** Las variantes estructurales deben usar alelos simbólicos cuando
  corresponda y los campos estándar SVTYPE, END y SVLEN para indicar el tipo, la
  posición final y la diferencia de longitud. Toda variante debe usar la
  representación estándar más completa que corresponda y cada campo estándar
  debe conservar su semántica oficial.
- **REQ-VER-19.** Una variante compleja debe conservar información suficiente
  para reconstruir el cambio real simulado. El VCF debe seguir siendo estándar,
  válido e interoperable por sí mismo. Cuando el estándar no pueda conservar
  todo el detalle, la reconstrucción completa debe usar el VCF junto con el
  complemento enlazado; ninguna limitación debe provocar una pérdida silenciosa
  de la verdad.
- **REQ-VER-20.** Los alelos de repeticiones en tándem deben usar REF y ALT
  explícitos cuando su tamaño sea manejable y una representación simbólica para
  tramos largos o complejos.
- **REQ-VER-21.** El VCF de verdad debe reflejar los cambios de copias y, cuando
  estén activadas, las mutaciones internas de interrupciones, pureza o motivo en
  el haplotipo afectado.
- **REQ-VER-22.** En ambas representaciones, el VCF debe exponer mediante campos
  estándar o adicionales válidos todo el motivo, número de copias, longitud,
  pureza, interrupciones y fase que permita expresar sin pérdida la versión
  oficial seleccionada. El detalle restante debe conservarse en la única entrada
  complementaria del caso, nunca únicamente en el manifiesto.
- **REQ-VER-23.** La cabecera del VCF debe documentar los campos estándar o
  adicionales empleados para representar las repeticiones en tándem y cualquier
  enlace con el complemento. Todo campo `INFO` adicional debe declararse en la
  cabecera conforme al estándar. Todo registro usa `GPS_CASE`, declarado con
  `Number=.` y `Type=String`; contiene los `case_id` originales únicos en orden
  canónico de todos los casos que aportan algún ALT al registro. Solo los casos
  con pérdida usan `GPS_LOSS` y `GPS_COMPLEMENT` conforme a DTC-065. Otros
  campos específicos de repeticiones en tándem se concretarán en su fase.
- **REQ-VER-24.** Si cualquier variante o genotipo complejo no puede
  representarse completamente y sin pérdida en VCF, el sistema debe generar una
  única entrada en un fichero estructurado complementario para ese caso. No debe
  crear una entrada para un caso plenamente representable ni generar el fichero
  cuando ningún caso lo necesite.
- **REQ-VER-25.** Un caso de verdad debe corresponder a uno o varios registros
  VCF y, solo cuando exista pérdida de detalle, a exactamente una entrada
  complementaria. Todos los registros del caso que necesiten ese detalle deben
  incluir en `INFO` el código breve y controlado `GPS_LOSS`, nunca texto libre,
  y el identificador local estable `GPS_COMPLEMENT` que enlace la misma entrada.
  Un registro puede enlazar uno o varios casos originales mediante `GPS_CASE`
  sin crear un caso compuesto; un mismo caso aparece una sola vez en la lista,
  aunque aporte varios ALT. `truth.jsonl` debe conservar para cada alelo final
  la relación explícita ALT↔lista canónica de casos fuente. La cardinalidad
  positiva entre varios casos y una entrada complementaria futura se concretará
  en F5.5 únicamente a partir de una pérdida real.
- **REQ-VER-26.** El complemento debe conservar todo el detalle que se perdería
  en VCF y que sea necesario para reconstruir la secuencia o cambio real, el
  haplotipo y la fase, las coordenadas y el genotipo simulado. El manifiesto debe
  enlazar los artefactos, pero no sustituir ese detalle. El fichero estructurado
  obligatorio de verdad sigue siendo la fuente canónica completa.
- **REQ-VER-27.** El VCF o, cuando corresponda, el fichero de verdad obligatorio
  debe conservar la dosis alélica para la ploidía efectiva de 1 a 10 copias.
- **REQ-VER-28.** Por cada conjunto definido por una referencia de mapeo, el
  sistema debe generar, cuando se solicite, exactamente un VCF de verdad
  multimuesta con todos los individuos del conjunto como columnas. Sus
  identificadores deben coincidir exactamente con los usados en BAM o CRAM,
  grupos de lectura, genotipos completos y manifiesto.
- **REQ-VER-29.** Un VCF de verdad no debe incluir individuos alineados contra
  otra referencia ni mezclar conjuntos. No deben generarse VCF individuales.
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
- **REQ-VER-42.** Las translocaciones deben usar la representación estándar más
  completa que corresponda. Una nueva adyacencia representada mediante BND debe
  contener sus dos puntos de ruptura, normalmente como dos registros relacionados
  mediante `MATEID`, `EVENT` u otros campos estándar aplicables. Una translocación
  recíproca puede contener dos nuevas adyacencias y requerir cuatro BND; no debe
  imponerse una cardinalidad universal de dos registros por translocación.
  Cuando el VCF no conserve todo el evento compuesto, todos los registros que
  necesiten ese detalle deben enlazar su única entrada complementaria.
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
- **REQ-VER-47.** El VCF no debe emitir un SNP con `REF=N` ni inventar una base
  concreta para representarlo. Otros cambios deben conservar `N` en `REF` o
  `ALT` solo cuando la versión oficial seleccionada permita una representación
  estándar válida y sin pérdida. Cuando no sea posible, debe aplicarse la única
  entrada complementaria del caso conforme a REQ-VER-24 a REQ-VER-26. La verdad
  estructurada canónica debe conservar siempre cada `N` exactamente.

Cuando se genera, el VCF de verdad recoge los genotipos definidos por la
simulación, no una estimación basada en la evidencia observada. Su finalidad
dentro de este proyecto es conservar la descripción exacta de lo simulado y
mantener la coherencia entre las salidas. Si se seleccionan varias referencias,
todas las muestras forman parte de todos los conjuntos, pero los archivos de
esos conjuntos deben permanecer separados.

### 4.6. Validación interna de lecturas y alineamientos

En esta sección, el soporte alélico es la cantidad o proporción de lecturas que
contienen un alelo concreto. Una comprobación determinista evalúa una regla
exacta que siempre debe cumplirse, como la coherencia de una referencia, una
identidad o una procedencia. Un resultado estadístico describe una magnitud que
puede diferir de su objetivo por el muestreo aleatorio. En v1 estos dos tipos no
comparten una política de aceptación.

- **REQ-VAL-01.** Siempre que se generen FASTQ o BAM o CRAM, el sistema debe
  comprobar internamente sus invariantes deterministas contra la verdad
  biológica simulada y producir por separado informes estadísticos de las
  magnitudes aleatorias y los artefactos configurados.
- **REQ-VAL-02.** La validación debe usar el genoma, los haplotipos, la ploidía,
  los genotipos y la procedencia verdadera de las lecturas definidos por la
  simulación. No debe ejecutar un *SNP caller* ni llamar, imputar, fasear,
  corregir o inferir genotipos a partir de BAM o CRAM.
- **REQ-VAL-03.** La validación debe comprobar que las lecturas proceden de los
  individuos y haplotipos correctos, salvo las procedencias diferentes que la
  configuración introduzca expresamente como contaminación u otro artefacto.
  Antes de destruir las claves técnicas transitorias, un reconstructor
  independiente por flujo debe comparar cada pareja final con la pareja esperada
  desde su clave, verdad biológica y transformación autorizada. Debe comprobar
  bases, calidades, mates, contaminación y transferencia técnica sin publicar ni
  conservar un vínculo individual por lectura.
- **REQ-VAL-04.** La validación debe comparar la cobertura objetivo común con la
  cobertura WGS básica realizada de cada biblioteca. Debe comprobar
  por separado la cobertura final del FASTQ y la cobertura alineada observada
  contra cada referencia, sin compararlas como objetivos equivalentes.
- **REQ-VAL-05.** La validación debe comprobar que el soporte alélico observado
  se informa junto con la dosis alélica, la ploidía simulada y los efectos
  configurados que puedan alterarlo. No debe exigir un balance exacto entre
  alelos o copias ni usar la desviación estadística para aceptar o rechazar el
  resultado.
- **REQ-VAL-06.** La validación debe comparar los valores objetivo y realizados
  de errores de secuenciación, contaminación biológica y técnica, duplicados y
  las demás degradaciones activadas. Estas diferencias se informan y no se
  someten a tolerancias de aceptación en v1.
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
  valores observados y las desviaciones detectadas. Para cada resultado
  estadístico debe registrar como mínimo el ámbito, el objetivo o valor esperado,
  el valor realizado, el tamaño de la muestra o número de unidades observadas,
  las desviaciones absoluta y relativa y los parámetros o la explicación
  pertinentes. El manifiesto debe enlazar estos resultados con los FASTQ, BAM o
  CRAM, muestras, grupos y conjuntos a los que se refieren.
- **REQ-VAL-11.** Toda inconsistencia estructural o lógica debe tratarse como un
  error que detiene la ejecución e impide marcar como completo el conjunto
  afectado. Esto incluye, como mínimo, archivos inválidos, referencias
  incorrectas, identificadores incompatibles, genotipos imposibles, cabeceras o
  índices incoherentes, lecturas o parejas con procedencia o haplotipo imposible
  y artefactos observados sin explicación en la configuración, la verdad o una
  decisión trazada.
- **REQ-VAL-12.** Las diferencias atribuibles al muestreo aleatorio en cobertura,
  soporte alélico, tasas realizadas de errores, contaminación, duplicados y
  demás degradaciones deben registrarse como resultados estadísticos. En v1 no
  se comparan con tolerancias de aceptación, no se clasifican como error y no
  detienen la ejecución; la configuración no puede elevarlas a error. Las
  desviaciones entre N objetivo y N realizado conservan su advertencia explícita
  no invalidante. N realizado igual a 0 aplica además la política de
  REQ-CONF-104. Ne realizado no se compara con un objetivo en v1.
- **REQ-VAL-13.** Para cada comprobación estadística, el resultado estructurado y
  el manifiesto deben registrar los campos mínimos de REQ-VAL-10 sin asignar una
  categoría de aceptación por tolerancia. Si una medida no puede calcularse,
  debe registrarse como «no evaluable» junto con la causa, nunca como cero ni
  como éxito o error. Para N deben registrar la desviación cuantificada y la
  advertencia no invalidante; N realizado igual a 0 se resuelve además conforme
  a REQ-CONF-104. Ne realizado se registra sin objetivo ni advertencia de
  desviación, o como «no evaluable» con causa.
- **REQ-VAL-14.** La validación WGS debe comprobar el objetivo excluyente de cada
  biblioteca, el denominador obtenido al sumar una vez las secuencias de la
  referencia estructural, la distribución aproximadamente Poisson anterior a
  sesgos y el reparto realizado por cromosoma y copia o haplotipo. También debe
  comprobar que haplotipos y referencias de mapeo no alteraron el denominador y
  que las cuatro medidas de REQ-WGS-11 permanecen separadas.
- **REQ-VAL-15.** La validación debe distinguir una variación aleatoria admisible
  del soporte alélico de un modelo intencionado de desequilibrio o pérdida
  alélica, que no existe en la primera versión.
- **REQ-VAL-16.** En F4, la validación de contaminación debe comprobar límites,
  compatibilidad, selección reproducible sin reemplazo entre donantes elegibles
  de la misma población, sustituciones biológicas, transferencias técnicas,
  descartes, totales por tanda y cobertura final. No valida matriz, pesos ni
  optimización global, que quedan en DTF-052.
- **REQ-VAL-17.** La validación poblacional debe comprobar N objetivo, intentos,
  éxitos y fallos por causa, N realizado, su desviación y su aviso, además de Ne
  realizado o su estado «no evaluable» con causa en cada subpoblación y
  generación. Debe
  comprobar que Ne realizado procede solo de los progenitores con descendencia
  y de la desigualdad o varianza de sus contribuciones.
  También debe comprobar autofecundaciones esperadas y realizadas, contribución
  familiar, rutas de migración, supervivencia y extinciones con las políticas
  acordadas. Para el SNP bialélico seleccionado debe validar su identidad única,
  protección y dosis alternativa para ploidías 1–10, la presencia de todos los
  valores genéticos exigidos y la finitud y no negatividad de las aptitudes.
  Antes de seleccionar debe comprobar que existe al menos una aptitud positiva.
  También debe comprobar que el carácter coincide con el valor genético, que no
  se aplicaron ruido ambiental ni heredabilidad y que la función respeta el
  sentido direccional. Cualquier ausencia, duplicación, ambigüedad, valor no
  finito, aptitud negativa o máximo cero debe constar como error previo a la
  supervivencia, nunca como extinción.
- **REQ-VAL-18.** La validación de mutación y recombinación debe comprobar tasas
  efectivas, oportunidades válidas, medias Poisson y eventos sorteados,
  aceptados y activos por clase y ámbito. Debe comprobar que la ploidía y las
  copias actuales cuentan, que ámbitos solapados no duplican oportunidades, que
  destinos o extremos condicionados no multiplican el denominador y que el
  recuento solo se sorteó al inicio de la clase. Debe comprobar además el orden
  recombinación→mutación, mapas completos, ploidías pares iguales, equivalencia
  estructural canónica antes de formar bivalentes y emparejamientos
  reproducibles. Debe comprobar también el mínimo de un entrecruzamiento en cada
  cromosoma y en los *scaffolds* que superen su umbral efectivo, y que los
  *scaffolds* sin recombinación transmitieron intactas las copias elegidas. Para
  `N` debe comprobar las exclusiones de oportunidades específicas de cada clase
  sin excluir intervalos puramente estructurales que puedan conservarla.
- **REQ-VAL-19.** Para cada individuo muestreado, la validación debe comprobar
  que los segmentos reconstruyen sin interrupciones su genoma, que cada segmento
  alcanza una copia haplotípica fundadora y que las mutaciones de novo se
  conservan como capas separadas sin borrar ese origen. Cuando existan cambios
  estructurales, debe comprobar además la correspondencia de coordenadas y la
  orientación entre los segmentos actuales y fundadores que sean relacionables.
- **REQ-VAL-20.** Cuando la población inicial use el modelo de equilibrio de
  Hardy-Weinberg, la validación debe comparar las frecuencias y recuentos
  genotípicos y alélicos esperados y realizados según las fórmulas binomial o
  multinomial de REQ-EVO-33 para ploidías 1–10. Debe tratar la variación finita
  como estocástica y no exigir recuentos genotípicos exactos. En modo controlado
  debe comprobar por separado los recuentos alélicos fijados y la distribución
  genotípica aleatoria condicionada.
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
  global o de subpoblación diploide, que no se usó con ploidía 1 o mayor que 2
  ni en ámbitos locales y que toda distribución explícita prevaleció en su locus
  o bloque. Debe comprobar las fórmulas bialélica o multialélica de
  REQ-EVO-47/48, incluidos F=0, un valor intermedio y F=1. F realizado debe
  calcularse solo en diploides; para las demás ploidías debe figurar como «no
  aplicable». Un F realizado diploide negativo debe conservarse y no marcarse
  como error.
- **REQ-VAL-24.** Cada evento de muestreo debe validar antes de generar datos que
  sus individuos explícitos y aleatorios existen, pertenecen a la generación y
  estrato correctos, no están duplicados y son suficientes. Una solicitud mayor
  que el conjunto elegible debe detener el evento sin reemplazo ni reducción
  silenciosa del tamaño.
- **REQ-VAL-25.** La validación técnica debe comprobar que cada evento usa una
  sola tecnología y un perfil cuando genera FASTQ, BAM o CRAM, y un solo
  mapeador cuando genera alineamientos. Debe comprobar que contiene exactamente
  una tanda y una biblioteca por individuo, que todas sus bibliotecas se tratan
  conjuntamente y que sus FASTQ no proceden de otro evento. Dentro del evento,
  los mismos FASTQ y mapeador deben aplicarse a todas las referencias
  seleccionadas. Un evento que solo solicite verdad genotípica o VCF debe
  rechazarse si declara tecnología, perfil, biblioteca, tanda o mapeador.
- **REQ-VAL-26.** La validación debe comprobar la unicidad global y coherencia de
  cada grupo de lectura, que distingue evento, biblioteca y tanda, y que
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
  que los espacios de identificadores no crean enlaces entre referencias, que
  el archivo cumple la versión oficial seleccionada y que todos los campos
  estándar mantienen su semántica.
- **REQ-VAL-29.** Cuando varios eventos compartan individuos, referencia y verdad
  genotípica, la validación debe comprobar la identidad de identificadores
  locales de variante, coordenadas, REF, ALT y genotipos entre sus VCF, sin
  exigir identidad de FASTQ, BAM, CRAM o degradaciones.
- **REQ-VAL-30.** La validación de familias debe comprobar la cardinalidad y el
  tipo de cada evento, la simetría de los parentales de un cruzamiento, la
  ploidía par e igual por secuencia, la identidad del estado estructural canónico
  de todas las copias reproductivas y la trazabilidad de las mutaciones de novo
  cuando estén activadas. Debe aceptar diferencias en SNP, *indels* pequeños y
  repeticiones, y rechazar heterocigosis o diferencias en cada clase de variante
  estructural. Una SV fijada de forma estructuralmente idéntica debe ser
  compatible aunque conserve historias distintas. Una SV de novo debe quedar
  trazada y, si su portador alimenta otro evento, este debe aplicar de nuevo la
  regla de compatibilidad. Para seleccionar descendientes debe comprobar los
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
  validación debe comprobar que existió un orden interno explícito único, que
  cada evento recibió el estado producido por el anterior y que se respetaron
  todas las dependencias. Para toda generación debe comprobar además que las
  ocho etapas siguieron REQ-REP-24, con cada etapa aplicada u omitida en su lugar.
  Un orden demográfico ausente, ambiguo o imposible o una solicitud de orden de
  procesos alternativo debe constar como rechazo previo a la simulación.
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
  los historiales resumidos. La ausencia de un estado completo obligatorio, la
  presentación de un resumen como estado completo o cualquier reducción
  automática del detalle obligatorio por escala debe ser un error.
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
- **REQ-VAL-46.** Para toda generación aleatoria de cambios, la validación debe
  comprobar que los explícitos se aplicaron primero y en orden, que las clases
  aleatorias siguieron el orden fijo de REQ-IND-55 y que cada evento se eligió
  del espacio válido actualizado. Debe comprobar por clase los números sorteado,
  aceptado antes del fallo y finalmente activo. La insuficiencia de destinos
  debe constar como fallo, nunca como reducción silenciosa, advertencia ni
  agotamiento de un número de reintentos.
- **REQ-VAL-47.** La validación de verdad debe comprobar que cada conjunto
  solicitado contiene exactamente un VCF multimuesta y que referencias distintas
  producen VCF independientes. Para cada caso debe comprobar la cardinalidad
  entre caso, registros y complemento, las declaraciones de todo campo `INFO`
  adicional y la validez de los enlaces locales. Debe cubrir como mínimo un caso
  representable sin complemento, un caso complejo representable mediante un
  registro y una entrada complementaria, una nueva adyacencia con dos BND y una
  entrada, y una translocación recíproca con todos sus BND necesarios y una
  entrada. En cada caso complejo, el par VCF y complemento debe reconstruir todo
  el detalle conservado por la verdad canónica.
- **REQ-VAL-48.** La validación debe comprobar que todo recuento de variantes
  aleatorias se obtuvo mediante Poisson y fue reproducible con la misma semilla,
  configuración, entradas y versión. Para cada mutación de repetición en tándem
  debe comprobar que el cambio fue exactamente una unidad completa y que su
  dirección respetó el sesgo efectivo. Debe constar como rechazo previo a la
  simulación cualquier distribución distinta, modelo alternativo o salto de
  varias unidades, identificando parámetro y valor.
- **REQ-VAL-49.** La admisión debe comprobar que ninguna configuración solicita
  calibración empírica, aporta FASTQ o referencia para calibrarla ni selecciona
  un perfil empírico reutilizable. Cualquiera de esos campos debe rechazarse
  antes de simular con un diagnóstico que identifique el campo y la exclusión de
  v1. Una ejecución válida con un perfil tecnológico incorporado y
  sustituciones permitidas no debe producir ningún artefacto de perfil empírico
  ni enlace a él en sus manifiestos.
- **REQ-VAL-50.** La validación debe comprobar que ninguna `N` se resolvió,
  sustituyó ni creó silenciosamente. Debe cubrir como mínimo: rechazo de SNP y
  locus seleccionado sobre `N`; exclusión de un motivo de repetición con `N`;
  una SV cuyo intervalo contiene `N`; inserciones aleatoria y explícita;
  herencia y recombinación; una lectura con `N` y Phred 0 propagada a BAM o
  CRAM; una ventana compuesta solo por `N` sin sesgo GC; cobertura que la cuenta;
  y representación VCF o complemento sin inventar bases.
- **REQ-VAL-51.** La admisión debe comprobar la resolución de sustituciones de
  REQ-CONF-116 a REQ-CONF-123 antes de simular. Debe cubrir la herencia de la
  regla base, cada jerarquía interna, una regla combinada dominante, conflictos
  entre región y subpoblación y entre individuo y región, duplicados exactos,
  regiones solapadas aunque compartan valor, selectores no permitidos y la
  separación entre parámetros biológicos y técnicos. Cada rechazo debe incluir
  el parámetro, el objeto y las reglas o selectores causantes.
- **REQ-VAL-52.** Antes de derivar una referencia reducida, la validación debe
  comprobar la exclusión entre objetivos, dominios y unidades, tamaño elegible,
  posibilidad de completar los intervalos sin solapamiento, cardinalidad del
  residual y compatibilidad de IDs. Todo fallo debe ocurrir antes de publicar el
  FASTA derivado o su mapa.
- **REQ-VAL-53.** La validación de una referencia reducida debe comprobar con el
  FASTA de origen y el mapa que cada intervalo existe, no se solapa, conserva
  orientación y orden, que la concatenación reconstruye exactamente cada
  cromosoma derivado y que tamaño y SHA-256 coinciden. Debe comprobar por
  separado que toda `N` seleccionada se conservó.
- **REQ-VAL-54.** La aceptación de la aplicación Marimo debe comparar el conjunto
  TOML exportado por la interfaz con un conjunto equivalente escrito a mano y
  demostrar que ambos producen la misma `EffectiveConfig` o los mismos
  diagnósticos, incluidos IDs, orden, enlaces y resolución de rutas relativas
  entre el documento general y sus escenarios. También debe rechazar antes de
  validar o ejecutar cualquier incompatibilidad entre la versión de la
  aplicación y las versiones de esquema o catálogos del binario.
- **REQ-VAL-55.** La validación de la matriz de frecuencias SNV debe comprobar
  el par TSV+JSON y sus SHA-256, referencia y coordenadas, coincidencia literal
  de `REF`, columnas y semántica `alt_af`, combinaciones población×ALT,
  duplicados, sumas multialélicas, filtrado global de ceros, IDs y orden,
  poblaciones, espectros no solapados, regiones no solapadas, reparto exacto o
  estocástico, perfiles, modelo Hardy-Weinberg o `F`, confirmaciones limitadas y
  aislamiento de subsemillas. Para el modo `exact` debe demostrar cantidades
  exactas por fila, estabilidad al reordenar la entrada, ausencia de asociación
  con la coordenada, aislamiento entre poblaciones, consumo cero para filas
  `fixed` o con cuota cero y separación entre el flujo de asignación y el de
  valor. Debe rechazar VCF y semánticas de frecuencia
  alternativas. Sobre corpus y generador independientes debe comparar matriz,
  frecuencias y genotipos esperados y realizados sin exigir que una muestra
  finita reproduzca exactamente las frecuencias objetivo.
- **REQ-VAL-56.** La validación del panel de haplotipos debe reconstruir cada
  haplotipo local de entrada, derivar sus frecuencias alélicas y comprobar sus
  sumas y cobertura completa por población y bloque en
  `panel_frequencies.tsv`, incluidos ceros y fijados, sin valores globales ni
  heredados. Debe verificar que los bloques no se solapen, que cada sitio esté
  dentro de su intervalo, que solo existan SNP con alelos A/C/G/T y que no haya
  sitios repetidos ni solapados. Debe comprobar la cobertura cartesiana exacta
  bloque×haplotipo×sitio de `panel_haplotype_alleles.tsv`, la unicidad de cada
  asignación y que todo `allele_index` pertenezca al catálogo del sitio; no debe
  completar ausencias como `REF`. También debe comparar los vectores completos
  de un mismo bloque y rechazar todos los haplotipos equivalentes, incluidos los
  de frecuencia 0, informando bloque e IDs sin fusionarlos ni sumar
  frecuencias. Antes de crear individuos debe demostrar que
  cada `block.reference_id` coincide con
  `metadata.reference.reference_id`, que `metadata.reference.sha256` coincide
  con el SHA-256 recalculado del FASTA autoritativo y que cada `sequence_id`
  existe en ese FASTA; debe rechazar mezclas, divergencias o secuencias
  inexistentes antes de interpretar el panel. También debe demostrar que
  ningún locus de la matriz `alt_af` cae dentro de un intervalo del panel; cada
  conflicto debe informar referencia, secuencia, posición y fila o identidad de
  matriz, además del `block_id` y su intervalo, sin elegir una fuente ganadora.
  También debe exigir exactamente los cinco ficheros canónicos, comprobar la
  clase y versión del JSON y verificar ruta, tamaño y SHA-256 de los cuatro TSV y
  el enlace a la referencia. El contenido tabular solo debe proceder de los TSV,
  no de una copia en los metadatos. En modo intacto debe detectar cualquier
  combinación intrabloque ausente del panel. En la alternativa recombinante debe
  usar un oráculo independiente que no llame a ninguna función de producción
  examinada. Debe reconstruir los HAP del panel desde referencia, sitios y
  alelos, y cada mosaico desde `recombination_segments.tsv`; recalcular el ID de
  contenido; verificar todas las asignaciones; comprobar que la suma de
  `copy_count` es C=N×P y que cada frecuencia es `copy_count/C`. Debe verificar
  byte a byte las cabeceras, tipos, orden y claves de las tres tablas de salida,
  rechazar columnas extra —incluidas clase o versión— y contrastar esas dos
  propiedades con sus enlaces de esquema/manifiesto. El oráculo debe probar el
  constructor uniforme cerrado por DTC-055 sin depender de mapas ni semántica
  meiótica pendientes en DTPD-111.

  El corpus debe incluir el rechazo de `r=0` en modo recombinante, un cruce
  conocido, cruces múltiples, invisibles y cancelados, resultados repetidos y
  F=0, intermedio y 1; F=1 debe producir homocigosis sin pedigrí ni
  autofecundación. Debe comprobar repetición exacta, aislamiento al añadir o
  reordenar individuos/bloques y paridad Marimo↔TOML. Debe corromper
  deliberadamente cada tabla, hash, referencia, intervalo, alelo, frecuencia y
  enlace y demostrar su detección.

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
  individuos declarados y existentes antes de la simulación y el valor efectivo
  aplicado a cada individuo afectado por una sustitución. También debe distinguir los valores definidos para el
  genoma, el cromosoma o *scaffold* y la región, y el valor efectivo resultante
  de aplicar REQ-CONF-116 a REQ-CONF-123. Para cada clase de
  variante debe registrar si se utilizó una cantidad media esperada o una
  densidad y, en el segundo caso, si la densidad se expresó por base, kilobase o
  megabase. También debe registrar si se aplicó una distribución espacial
  uniforme o se configuraron regiones, incluidas las regiones con densidad
  cero.
- **REQ-MAN-11.** Para cada clase de variante, ámbito general o región e
  individuo al que corresponda, el manifiesto debe registrar la cantidad media
  esperada, el número realizado o sorteado, el número aceptado antes de un
  posible fallo y el número finalmente activo, junto con la distribución
  Poisson y su media efectiva.
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
- **REQ-MAN-16.** El manifiesto debe enumerar individualmente a cada individuo
  cuyo estado completo se materialice y a todos los individuos de familias y
  poblaciones segregantes. Para cada uno debe relacionar identificador estable,
  generación de origen, parentales cuando existan, observaciones y archivos.
  Las generaciones evolutivas no materializadas deben aparecer mediante sus
  resúmenes, no mediante listas individuales completas. Debe distinguir los
  identificadores de individuo, evento de muestreo y conjunto.
- **REQ-MAN-17.** En una ejecución evolutiva, el manifiesto debe registrar la
  trayectoria configurada de N objetivo y, por generación, fase o tramo, N
  objetivo, intentos, éxitos, fallos por causa, N realizado, desviación y aviso.
  Debe registrar Ne realizado o «no evaluable» con causa, sin Ne objetivo ni
  desviación asociada. Debe hacerlo para la población completa y para cada
  subpoblación cuando corresponda.
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
  activos, desactivados o no aplicables y los parámetros y valores efectivos
  aplicados por fase, generación y subpoblación cuando corresponda. Cada proceso
  omitido debe conservar su etapa fija y la causa de omisión.
- **REQ-MAN-21.** El manifiesto debe registrar el orden efectivo en que se
  aplicaron las ocho etapas de cada generación, que siempre debe ser el orden
  fijo de REQ-REP-24. Debe distinguir en cada etapa lo aplicado de lo omitido
  por estar desactivado o no ser aplicable, y registrar por separado el orden
  interno explícito de los eventos demográficos cuando haya varios.
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
  individuos y enumerarlos individualmente, si se solicitaron salidas
  adicionales y cuáles fueron los BAM, CRAM o VCF generados.
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
  registrar el modelo escalonado simple, la tasa y el sesgo efectivos, el locus
  y haplotipo afectados y si añadió o eliminó exactamente una unidad completa.
  Cuando estén activadas, debe registrar por separado las modificaciones internas
  de interrupciones, pureza o motivo.
- **REQ-MAN-31.** Cuando exista un fichero complementario de verdad, el
  manifiesto debe enlazarlo con el VCF, la referencia de mapeo y el conjunto
  correspondientes y conservar los identificadores locales de los casos y sus
  relaciones uno-a-varios con los registros VCF. Debe registrar exactamente una
  entrada por caso con pérdida de detalle, sin usar el manifiesto para almacenar
  el detalle que corresponde al complemento. Esos identificadores son locales
  al conjunto y no deben enlazar variantes entre referencias.
- **REQ-MAN-32.** El manifiesto debe distinguir las salidas mínimas, las
  activadas de forma predeterminada, las solicitadas expresamente y las
  desactivadas, y registrar todos los archivos generados para satisfacerlas.
- **REQ-MAN-33.** El manifiesto debe registrar las estimaciones previas de
  volumen, almacenamiento, memoria y tiempo, las advertencias emitidas y, cuando
  pueda medirse, el consumo real de la ejecución. Para duplicaciones debe incluir
  el máximo local efectivo, si se elevó sobre el valor predeterminado, el
  presupuesto de recursos configurado y el resultado de compararlo con la
  estimación previa. La estimación debe incluir el contrato obligatorio de
  trazabilidad y dejar constancia de cualquier rechazo previo por no poder
  producirlo, sin registrar un nivel reducido alternativo.
- **REQ-MAN-34.** Para cada conjunto y cada FASTQ, el manifiesto debe registrar
  el modo realista, la tecnología, el perfil, los parámetros efectivos y la
  semilla o el flujo aleatorio aplicable. El etiquetado del modo directo queda
  aplazado en DTF-051.
- **REQ-MAN-35.** Para cada FASTQ final generado, el manifiesto debe registrar
  sus hashes y tamaños, muestras, perfil y resumen agregado de procedencia. No
  incluye un fichero ni vínculos por lectura; relaciona cada conjunto con sus
  FASTQ y los recursos técnicos que se activen en fases posteriores.
- **REQ-MAN-36.** Para cada evento con alineamientos, el manifiesto debe registrar el
  nombre y la versión del único mapeador, el comando efectivo completo y sus
  parámetros, cada referencia y los FASTQ usados. Cuando se use minimap2, debe
  distinguir el *preset*, los valores predeterminados y las sustituciones de la
  configuración.
- **REQ-MAN-37.** En F4, el manifiesto debe registrar el objetivo de cobertura
  común, las coberturas WGS básica y final del FASTQ, por biblioteca y muestra,
  y el denominador utilizado. No admite sustituciones por individuo, cromosoma,
  *scaffold* o región. La cobertura alineada observada se añade cuando exista
  mapeo en DTF-055 y no se presenta como objetivo.
- **REQ-MAN-38.** Para cada fuente contaminante, el manifiesto debe registrar su
  identificador y origen, la fracción objetivo, la fracción realizada y los
  recuentos de lecturas o parejas propias retiradas y contaminantes incorporadas,
  además de la semilla o flujo aleatorio de selección. En datos *paired-end*,
  debe registrar las parejas completas como unidad.
- **REQ-MAN-39.** F4 rechaza una fuente contaminante externa y no crea su copia
  ni entrada de manifiesto. Cuando esa capacidad se incorpore desde DTF-052, el
  manifiesto identificará la fuente, ubicación original, copia utilizada,
  resultado de validación, integridad y relación con la selección y su semilla o
  flujo aleatorio.
- **REQ-MAN-40.** Para la contaminación técnica de F4, el manifiesto debe
  registrar la tanda, bibliotecas, muestras y sus UDI `i7+i5`, la tasa de
  *index hopping*, los errores simulados y las cantidades realizadas de parejas
  reasignadas, sin asignar o filtradas. No conserva procedencia individual ni
  declara estrategias de índice distintas de UDI.
- **REQ-MAN-41.** El manifiesto debe registrar para cada FASTQ final sus hashes,
  tamaños y ausencia de FASTQ limpio adicional o fichero de procedencia. Debe
  resumir las parejas contaminantes, reasignadas, descartadas o sin asignar por
  muestra, causa y origen/destino cuando corresponda.
- **REQ-MAN-42.** Para los errores de secuenciación de F4, el manifiesto debe
  registrar la tecnología, perfil y versión, sus tipos de error, y los ajustes
  globales efectivos por mate. No admite ni registra valores por individuo,
  cromosoma, *scaffold* o región.
- **REQ-MAN-43 (retirado).** Registro en el manifiesto de entradas y resultados
  de calibración empírica; capacidad aplazada fuera de v1.
- **REQ-MAN-44 (retirado).** Enlace desde el manifiesto a un fichero de perfil
  empírico; artefacto excluido de v1.
- **REQ-MAN-45.** En F4, el manifiesto declara duplicados PCR no aplicables y
  ópticos a cero. Tasas, cantidades y relaciones de duplicados quedan aplazadas
  en DTF-056.
- **REQ-MAN-46.** El marcado de duplicados queda aplazado en DTF-055 y DTF-056.
- **REQ-MAN-47.** En F4, el manifiesto declara adaptadores desactivados y
  recuento cero; toda solicitud de activarlos se rechaza. El perfil de
  fragmentos cortos y sus parámetros o lecturas afectadas quedan en DTF-054.
- **REQ-MAN-48.** En F4, el manifiesto declara relación GC neutra. Curvas,
  intensidad, combinación con cobertura y métricas de corrección GC quedan
  aplazadas en DTF-056.
- **REQ-MAN-49.** La identificación de regiones para baja mapeabilidad, su
  relación con repeticiones o variantes estructurales y el mapeador, referencia
  y MAPQ observados quedan aplazados en DTF-055.
- **REQ-MAN-50.** El manifiesto de categorías y MAPQ de alineamientos directos
  queda aplazado en DTF-051.
- **REQ-MAN-51.** Para las lecturas quiméricas, el manifiesto debe registrar el
  perfil, la frecuencia predeterminada o sustituida y la cantidad realizada. No
  enlaza una verdad individual de fragmentos o puntos de unión en F4.
- **REQ-MAN-52.** El manifiesto de anomalías *paired-end* observables queda
  aplazado en DTF-055. F4 no publica tasas ni cantidades de esas anomalías.
- **REQ-MAN-53.** El manifiesto debe registrar cada grupo de lectura, su muestra
  o individuo, biblioteca, tanda, plataforma o tecnología, parámetros
  efectivos y archivos asociados mediante los mismos identificadores usados en
  FASTQ y, cuando exista mapeo, BAM o CRAM. Debe relacionar todos los grupos
  del individuo con su único BAM o CRAM final dentro de cada conjunto y
  referencia de mapeo. Debe enlazar cada FASTQ final separado con su grupo de
  lectura y con el BAM o CRAM final de destino cuando corresponda.
- **REQ-MAN-54.** El manifiesto debe identificar la única tecnología de cada
  evento técnico, su única tanda, su única biblioteca por individuo y su perfil
  y mapeo compatibles. Debe enlazar los eventos que usen las mismas muestras y verdad
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
  explicaciones registradas. Para cada resultado estadístico debe conservar el
  ámbito, el objetivo o valor esperado, el valor realizado, el tamaño de muestra
  o número de unidades observadas, las desviaciones absoluta y relativa y la
  explicación o los parámetros pertinentes. Si no se pudo calcular, debe
  registrar «no evaluable» y la causa. También debe identificar los conjuntos
  que no se marcaron como completos por una incoherencia determinista.
- **REQ-MAN-60.** Para cada biblioteca WGS de F4, el manifiesto debe registrar
  la preparación sin PCR, el objetivo original de cobertura, el número de
  parejas derivado y realizado y el perfil y las
  sustituciones aplicados.
- **REQ-MAN-61.** El manifiesto debe registrar la definición de cobertura usada,
  la referencia estructural y la longitud haploide del denominador, la cobertura
  objetivo, la WGS básica anterior a sesgos y artefactos, la final del FASTQ con
  bases propias y contaminantes separadas y la alineada observada por cada
  referencia de mapeo. Cuando pueda medirse, debe desglosar también la cobertura
  por cromosoma y copia o haplotipo. La cobertura alineada no debe etiquetarse
  como objetivo.
- **REQ-MAN-62.** En F4, para la contaminación el manifiesto debe registrar el
  rango, el porcentaje sorteado, los donantes elegibles y los aportes realizados
  por donante, los límites comprobados y, por fuente y mecanismo, parejas propias
  retiradas, contaminantes incorporadas, transferidas, descartadas, totales
  finales y cobertura final. Matriz dirigida, pesos y reparto optimizado quedan
  aplazados en DTF-052.
- **REQ-MAN-63.** F4 no genera ni enlaza eventos de control porque solo admite
  un evento. Tras DTF-050, el manifiesto deberá enlazar un control sin
  degradaciones con su evento relacionado mediante verdad, referencias, muestras
  biológicas y relación experimental, manteniendo separados sus archivos y
  subsemillas.
- **REQ-MAN-64.** El manifiesto debe enlazar el historial resumido de cada
  generación y registrar N objetivo heredado o modificado, intentos, éxitos y
  fallos por causa, N realizado, desviación y aviso, y la transición en que el
  objetivo determinó la descendencia. Debe registrar Ne realizado o «no
  evaluable» con causa, apareamientos, progenitores con descendencia,
  descendientes y estados genotípicos completos conservados.
- **REQ-MAN-65.** Debe registrar la tasa de autofecundación, sus valores
  esperados y realizados, la referencia, secuencia, coordenada e identificador
  del único SNP seleccionado, sus dosis, valores genéticos, valores del carácter, función de
  aptitud, aptitudes antes y después de normalizar,
  supervivientes, extinciones, el modo y los valores de migración solicitados y
  realizados por ruta y las rutas individuales. Debe registrar que el valor del
  carácter coincidió con el valor genético y que no se usó componente ambiental
  ni heredabilidad.
  También debe registrar el sentido direccional elegido y los parámetros
  efectivos de la función.
- **REQ-MAN-66.** Debe registrar las tasas y espectros de mutación, eventos
  realizados, mutaciones de novo por haplotipo, transmisión y orden efectivo de
  procesos de cada generación. Por clase y ámbito debe incluir tasa efectiva,
  oportunidades válidas, media Poisson y recuentos sorteado, aceptado y activo.
  Para `N` debe incluir las oportunidades excluidas por necesitar una base
  concreta y distinguirlas de los intervalos estructurales que la conservaron.
- **REQ-MAN-67 (retirado).** La primera versión no genera ni enlaza catálogos,
  copias, familias, mecanismos ni metadatos específicos de TE.
- **REQ-MAN-68.** Debe registrar el modo y unidades de recombinación, los mapas
  completos o tasas efectivas, los estados estructurales canónicos comparados,
  los individuos o parejas excluidos por incompatibilidad, los bivalentes
  formados y los entrecruzamientos realizados en los resúmenes generacionales.
  Para cada secuencia debe registrar su clasificación
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
  usa F en diploides, debe registrar su objetivo y el F realizado por
  subpoblación y generación, incluidas estimaciones negativas. Para ploidía 1 o
  mayor que 2 debe registrar F realizado como «no aplicable». Debe identificar
  HWE, F diploide o distribución explícita y, en el modo controlado, distinguir
  los recuentos alélicos fijados de la distribución genotípica aleatoria
  condicionada.
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
  Cuando el evento genere FASTQ, BAM o CRAM debe registrar que constituye su
  único experimento de secuenciación, su tecnología, perfil, tanda y la
  biblioteca única de cada individuo; cuando genere alineamientos debe registrar
  además su mapeador único. Un evento que solo materialice verdad o VCF no debe
  contener esos campos técnicos como si se hubieran aplicado.
- **REQ-MAN-74.** Dentro de un evento, el manifiesto debe separar los resultados
  por referencia de mapeo y relacionar cada BAM, CRAM o VCF exclusivamente con
  su conjunto. Debe conservar los identificadores biológicos de muestra, sin
  crear correspondencias entre las variantes de conjuntos diferentes.
- **REQ-MAN-75.** Las referencias y sus índices deben almacenarse una sola vez
  entre los recursos comunes de la ejecución. Los manifiestos de evento deben
  referenciarlos mediante su ruta y un identificador verificable, sin crear
  copias redundantes dentro de cada evento.
- **REQ-MAN-76.** El manifiesto general y cada manifiesto de evento deben incluir
  sumas de comprobación de todos los ficheros que enlazan, incluidas las copias
  utilizadas de configuraciones y entradas, y el algoritmo utilizado.
  El sistema debe comprobar esas sumas antes y después de exportar un evento.
  El algoritmo concreto se decidirá durante la planificación de la
  implementación.
- **REQ-MAN-77.** Los identificadores de grupo de lectura deben ser únicos en
  toda la ejecución y distinguir evento, biblioteca y tanda. El
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
  estado genotípico completo ni declarar un nivel configurable o una degradación
  automática de trazabilidad.
- **REQ-MAN-85.** Para cada receptora de F4, el manifiesto debe registrar el
  número de parejas limpias, porcentaje sorteado, objetivo entero derivado,
  donantes elegibles, aportes realizados y subsemillas aplicables, todo como
  recuentos agregados por receptora, donante y mecanismo. No conserva
  identificadores ni vínculos de origen de parejas individuales. La conversión
  conjunta, matriz y optimización de objetivos quedan en DTF-052.
- **REQ-MAN-86.** Para cada secuencia FASTA, el manifiesto o la salida
  estructurada de referencias debe conservar el identificador formado por el
  primer elemento de la cabecera y la descripción completa restante, además del
  resultado de validar el alfabeto A, C, G, T y N.
- **REQ-MAN-87.** El manifiesto debe registrar, por referencia, biblioteca y
  operación pertinente, el número o la proporción de `N`, sus exclusiones de
  oportunidades, los fragmentos sin bases conocidas y la comprobación de que
  no se resolvió ni se creó ninguna `N` de forma implícita. Debe enlazar los
  casos VCF que necesiten complemento para conservarla sin pérdida.
- **REQ-MAN-88.** El manifiesto y la configuración efectiva deben registrar,
  para cada valor sustituido, el parámetro y objeto, la regla ganadora, sus
  selectores y su procedencia. Una ejecución rechazada por máximas incomparables
  debe conservar esa misma información para todas las reglas máximas en el
  diagnóstico de admisión. Cuando una fase propietaria use relaciones
  jerárquicas, el manifiesto debe enlazar la procedencia real que esa fase haya
  definido; F1 no crea una evidencia de dominio sintética para relaciones que
  todavía no existen.
- **REQ-MAN-89.** Para cada referencia reducida, el manifiesto debe registrar
  los IDs de origen y derivado, rutas y SHA-256 de ambos FASTA, tamaño elegible,
  `region_size_bases`, modalidad y valor del objetivo, número exacto de bases
  derivadas, semilla raíz, subsemilla y versión de la selección.
- **REQ-MAN-90.** El manifiesto debe enlazar el mapa completo de coordenadas de
  cada referencia reducida, registrar su clase, versión, ruta, tamaño y SHA-256
  y permitir comprobar qué cromosomas de origen se incluyeron u omitieron y si
  se utilizó el único fragmento residual permitido.
- **REQ-MAN-91.** La creación desde frecuencias SNV debe registrar si la matriz
  fue suministrada, generada o reutilizada; rutas, clases, versiones, tamaños y
  SHA-256 del TSV, JSON y referencia; poblaciones disponibles y seleccionadas;
  parámetros declarados y efectivos de candidatos, Poisson/densidad, espacio,
  sustituciones y frecuencias; fuente, filtros y versión de perfiles; semilla,
  catálogo/versión de derivación, ámbitos, subsemillas y versiones de cada
  algoritmo consumidor; modelo Hardy-Weinberg o `F`; loci sorteados, descartes
  globales, loci activos y repartos antes/después; frecuencias y genotipos
  esperados y realizados; y toda autorización de genotipo incompatible o matriz
  vacía.
- **REQ-MAN-92.** La creación desde un panel de haplotipos debe registrar su ruta
  y la ruta y SHA-256 de `panel.metadata.json`, `panel_id`, referencia y
  población, bloques e intervalos, identidades locales, frecuencias objetivo y
  realizadas por bloque, modo con o sin recombinación, subsemillas y, cuando
  corresponda, mapas, puntos de cruce y procedencia por segmentos. Cuando se
  combine con una matriz `alt_af`, debe registrar ambos recursos y que se
  comprobó que ningún locus independiente cae dentro de los intervalos de
  bloque. En la alternativa recombinante fundadora debe registrar además la
  configuración efectiva, F, `recombination_probability_per_interbase`, versión
  de derivación, semillas y ámbitos, número de productores con cruce, número
  total de cruces, vectores nuevos y coincidencias con HAP del panel. Debe enlazar
  `founder_block_assignments.tsv`, `recombination_segments.tsv` y
  `realized_haplotype_frequencies.tsv` por clase, versión, ruta, tamaño y SHA-256,
  sin duplicar sus filas en el JSON.

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
  todos los eventos y poblaciones intermedias conforme a REQ-CRU-15. Este es el
  único alcance obligatorio de v1: no debe admitir niveles configurables ni
  reducirlo automáticamente por el tamaño del escenario.
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
  formato es `truth.jsonl` v1 y sus esquemas/validación siguen DTC-061.
- **REQ-SAL-07.** Cuando se use un FASTQ externo como fuente de contaminación,
  su copia en el directorio de resultados debe formar parte de las salidas
  activadas necesarias para que el experimento sea autocontenido y exportable.
- **REQ-SAL-08.** En el modo realista, una ejecución solo debe terminar con
  éxito si existe por cada grupo de lectura su FASTQ final —o pareja de FASTQ—,
  enlazado desde el manifiesto junto con su resumen técnico agregado. Debe
  existir exactamente un grupo y un FASTQ o pareja por individuo y evento; no
  debe existir un fichero de procedencia individual.
- **REQ-SAL-09.** La simulación de lecturas no debe crear una verdad técnica
  individual de diferencias, duplicados o procedencia. El manifiesto conserva
  exclusivamente los parámetros efectivos y los recuentos agregados definidos
  para cada artefacto.
- **REQ-SAL-10 (retirado).** Fichero estructurado reutilizable de perfil
  empírico como salida activada; artefacto excluido de v1.
- **REQ-SAL-11.** Todo BAM o CRAM final debe incluir su índice compatible entre
  las salidas activadas. Una ejecución no debe terminar con éxito si falla la
  validación conjunta del alineamiento, el índice y la referencia.
- **REQ-SAL-12.** Cuando se solicite VCF, la salida activada debe incluir el VCF
  multimuesta comprimido con BGZF y su índice compatible, salvo que se solicite
  expresamente la alternativa de texto sin comprimir. Debe incluir también el
  único fichero complementario del conjunto cuando al menos un caso lo necesite.
  La ejecución no debe terminar con éxito si falta una entrada obligatoria o
  falla la validación correspondiente.
- **REQ-SAL-13.** Siempre que la ejecución genere FASTQ o BAM o CRAM, los
  resultados estructurados de su validación interna deben formar parte de las
  salidas activadas y estar enlazados desde el manifiesto.
- **REQ-SAL-14 (retirado).** La primera versión no admite un catálogo FASTA de
  TE como entrada ni como salida activada.
- **REQ-SAL-15.** Toda ejecución evolutiva debe producir el historial
  estructurado resumido de todas sus generaciones y los estados genotípicos
  completos exigidos por REQ-REP-19 y REQ-REP-21.
- **REQ-SAL-16.** Cada muestreo evolutivo debe producir para sus individuos la
  salida detallada de recombinación y origen por segmentos exigida en
  REQ-REC-11 a REQ-REC-13. La misma obligación se aplica a una materialización
  adicional solicitada expresamente. Esta salida no debe extenderse a todos los
  individuos históricos por el mero hecho de generar resúmenes de cada
  generación.
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
- **REQ-SAL-27.** Toda referencia reducida solicitada debe publicar como unidad
  coherente su FASTA canónico, su mapa estructurado derivada↔origen y sus
  metadatos de integridad y reconstrucción. Si falla cualquiera de ellos, la
  referencia derivada debe quedar fallida y ninguno puede presentarse como
  resultado completo.

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
  necesarios, usando la configuración y los recursos que pueda detectar. La
  estimación debe incluir el único contrato obligatorio de trazabilidad y todas
  las salidas activadas.
- **REQ-NOF-06.** Si el escenario parece grande para los recursos disponibles,
  el sistema debe mostrar advertencias claras. Una advertencia no debe bloquear
  por sí sola. Debe impedirse el inicio si la configuración es inválida, falta
  inevitablemente un requisito mínimo o una estimación supera un presupuesto de
  recursos configurado, explicando la causa. No debe continuar con menos detalle
  ni elegir otro nivel de trazabilidad. Para poder ejecutar, la persona usuaria
  deberá reducir materializaciones, muestreos, individuos, lecturas u otras
  salidas solicitadas hasta que la estimación permita cumplir el mínimo.
- **REQ-NOF-07.** Las estimaciones no deben presentarse como exactas. Los modelos
  y umbrales de advertencia se elegirán y validarán durante la implementación.
- **REQ-NOF-08.** Para reducir el almacenamiento, el sistema no debe conservar
  una segunda copia limpia de los FASTQ ni un fichero de procedencia individual.
  Las estimaciones de espacio deben reflejar esta política.

### 4.10. Aplicación Marimo de configuración y ejecución

- **REQ-UI-01.** Toda capacidad y opción configurable implementada en v1 debe
  poder configurarse y lanzarse desde una aplicación Marimo, además de poder
  expresarse mediante TOML y ejecutarse con el binario. Ninguna opción
  implementada puede quedar disponible únicamente en TOML o línea de órdenes.
- **REQ-UI-02.** La aplicación debe organizar como mínimo las pestañas
  `General/Ejecución`, `Genoma`, `Individuos`, `Poblaciones`, `Evolución`,
  `Familias`, `Eventos técnicos`, `Salidas` y `Revisar y ejecutar`. Puede ocultar
  o deshabilitar pestañas y controles que no correspondan a la vía elegida o a
  una capacidad aún no implementada, indicando la causa sin eliminar su destino
  previsto.
- **REQ-UI-03.** Los escalares, enumeraciones y booleanos deben modificarse con
  controles tipados que conserven su unidad, dominio y ausencia o valor por
  defecto. La aplicación no debe convertir silenciosamente unidades, truncar
  enteros ni redondear el valor normativo mostrado en otra unidad.
- **REQ-UI-04.** Toda lista ordenada, entidad repetible y regla debe usar tablas
  editables que permitan añadir, eliminar y reordenar filas. La identidad y el
  orden normativo deben ser visibles y no pueden depender del orden accidental
  en que Marimo vuelva a ejecutar una celda.
- **REQ-UI-05.** Los selectores y ámbitos de cada regla deben mostrarse de forma
  explícita, incluidos individuo, población o subpoblación, fase o generación,
  referencia, secuencia, región o locus, grupo, bloque y evento técnico cuando
  la familia los admita. La interfaz debe impedir o diagnosticar un selector no
  permitido, pero no resolver conflictos con precedencia propia.
- **REQ-UI-06.** La aplicación debe importar, editar, previsualizar y exportar
  la configuración de la capacidad vigente sin cambiar su significado. En F2
  la unidad visible es un único TOML general de la vía directa. Debe conservar
  IDs, orden y rutas relativas que ya formen parte de ese documento y producir
  una configuración que el binario pueda validar; no se promete conservar
  espacios ni comentarios. El conjunto multidocumento y los TOML de escenario
  se incorporarán en la primera fase que implemente una vía que los consuma, no
  como infraestructura anticipada de F2.
- **REQ-UI-07.** La validación debe invocar explícitamente `validate` del binario
  Rust y asociar cada diagnóstico disponible con su pestaña, tabla, fila y campo.
  El estado validado debe identificar los documentos, el contrato y el conjunto
  de autorizaciones limitadas usados. Cualquier modificación o
  confirmación/retirada posterior de una autorización invalida el resultado
  previo y obliga a validar otra vez antes de habilitar la ejecución.
- **REQ-UI-08.** La aplicación nunca debe ejecutar automáticamente al abrirse,
  importar o modificar un valor. `run` solo se lanza mediante una acción humana
  explícita después de una validación vigente y una revisión visible de la
  configuración y el destino.
- **REQ-UI-09.** La aplicación debe invocar `validate` y `run` como procesos con
  argumentos separados, sin construir ni ejecutar una orden de *shell*. El
  conjunto de autorizaciones de `run` debe ser idéntico al de la validación
  vigente y solo puede contener `--allow-empty-snv-matrix` y
  `--allow-incompatible-genotypes`; una opción general o no catalogada debe
  rechazarse. La aplicación solo debe mostrar o transmitir una de esas
  autorizaciones cuando el binario declare que está efectivamente disponible,
  no solo catalogada para una versión futura. El
  binario, sus esquemas, catálogos, diagnósticos, configuración efectiva y
  manifiestos conservan toda autoridad semántica; Marimo no reimplementa reglas
  biológicas ni decide valores efectivos.
- **REQ-UI-10.** La elección del destino debe impedir selecciones vacías o
  accidentales y mostrar la política no-clobber antes de lanzar. La aplicación
  no puede autorizar sobrescritura, reparación o reutilización de un resultado
  existente que el binario rechazaría.
- **REQ-UI-11.** Durante validación y ejecución, la aplicación debe mostrar
  estado, código de salida, logs y diagnósticos. Al terminar, debe mostrar un
  resumen y enlaces a manifiestos y artefactos; no debe cargar completos en
  memoria FASTA, FASTQ, BAM, CRAM, VCF u otros artefactos genómicos grandes para
  presentar ese resumen.
- **REQ-UI-12.** El contrato aplicación↔binario debe declarar y comprobar las
  versiones admitidas del esquema TOML y la interfaz de línea de órdenes que
  utiliza. La aplicación solo mostrará controles y autorizaciones realmente
  implementados por el binario de la misma versión. Una incompatibilidad debe
  bloquear validación o ejecución con un diagnóstico explícito; F2 no crea un
  protocolo paralelo de capacidades, digests o migraciones.
- **REQ-UI-13.** La aplicación debe ser utilizable localmente en Ubuntu 24.04 y
  en el entorno local de desarrollo sobre un sistema de ficheros Linux. Debe
  ofrecer etiquetas, ayuda, navegación por teclado y estados que no dependan
  solo del color. Cada fase que añada una opción configurable debe comprobar un
  recorrido representativo UI→TOML→binario, sin exigir una auditoría completa
  de navegador o plataforma en cada hito. F2 prueba instalación desde el
  repositorio y arranque local; la acreditación oficial de Ubuntu 24.04,
  accesibilidad de distribución y OCI pertenece a F9.
- **REQ-UI-14.** En F2, Marimo debe permitir elegir una matriz de frecuencias SNV
  suministrada, generada o reutilizada; preparar una entrada mediante la misma
  orden auxiliar del binario y la acción «Vincular y verificar referencia»;
  seleccionar o cargar el espectro alélico; configurar población, candidatos,
  espacio, sustituciones, Hardy-Weinberg o `F`; y previsualizar matriz, sumas,
  descartes y distribuciones sin crear individuos. Las confirmaciones de
  genotipos incompatibles y matriz vacía deben mostrar individuos/loci o el
  resumen afectado y traducirse a los argumentos específicos, nunca a un
  `--force` general.
- **REQ-UI-15.** En la fase posterior correspondiente, Marimo debe permitir
  elegir el directorio del panel en una sección «Bloques haplotípicos en
  población» y ofrecer acciones explícitas para prepararlo mediante la orden
  auxiliar Rust y para validarlo. Debe mostrar un resumen de referencia,
  bloques, sitios, HAP y poblaciones, y permitir mapear cada población del panel
  a una población simulada. Una tabla por población y bloque debe ofrecer una selección excluyente
  de `haplotype_construction` entre `intact` y `recombinant`, separada de
  `frequency_sampling`, y admitir el valor general y las sustituciones por
  población y bloque. Solo debe mostrar y exportar
  `recombination_probability_per_interbase` cuando la construcción efectiva sea
  `recombinant`; el control debe ser tipado. En modo intacto debe mostrar
  `probabilistic` y `controlled`; en recombinante debe desactivar `controlled` y
  explicar la incompatibilidad. F solo debe mostrarse para diploides.

  La previsualización debe incluir frecuencias, `P(K>0)`, homocigosis esperada y
  conflictos con `alt_af`. Para preparar el recurso debe invocar la misma orden
  auxiliar del binario Rust y no reimplementar su validación. Cuando
  también se seleccione una matriz `alt_af`,
  la validación del binario debe mostrar cualquier locus de la matriz situado
  dentro de un intervalo de bloque antes de ejecutar. La
  interfaz debe mostrar que la recombinación puede alterar las asociaciones
  originales y no debe interpretar el panel fuera del binario Rust. Una
  validación vigente debe ser obligatoria antes de ejecutar. Tras la ejecución,
  Marimo debe mostrar enlaces y resúmenes de las tres tablas obligatorias sin
  cargar artefactos grandes completos.

## 5. Mapa conceptual del flujo

Esta sección es un resumen **derivado e informativo** para recorrer el sistema
de principio a fin. El contrato funcional reside exclusivamente en los
requisitos de la sección 4; ante cualquier diferencia, prevalece el requisito
identificado allí. Los criterios ACE citados en la tabla indican cómo se
observará el resultado, sin crear obligaciones adicionales.

| Etapa conceptual | Pregunta que permite seguir el flujo | Fuente normativa | Comprobación principal |
|---|---|---|---|
| Admisión | ¿La petición y sus ficheros forman una ejecución inequívoca? | REQ-CONF-01 a REQ-CONF-130, REQ-FMT-01 a REQ-FMT-07 y REQ-NOF-01 a REQ-NOF-08 | ACE-01, ACE-23, ACE-35, ACE-36, ACE-105, ACE-111, ACE-112 y ACE-116 |
| Referencias | ¿Qué secuencias y sistemas de coordenadas sostienen la verdad? | REQ-REF-01 a REQ-REF-36 | ACE-02 a ACE-08, ACE-95, ACE-110, ACE-112 y ACE-113 |
| Construcción de individuos | ¿Qué genomas faseados existen antes de observarlos? | REQ-IND-01 a REQ-IND-65, REQ-DIR-01 a REQ-DIR-05 y REQ-MUT-01 a REQ-MUT-15 | ACE-09 a ACE-12, ACE-28 a ACE-32, ACE-79 a ACE-81, ACE-95, ACE-96, ACE-103, ACE-107, ACE-109 a ACE-111 y ACE-116/117 |
| Inicialización poblacional | ¿Qué población inicial se deriva de las frecuencias y restricciones declaradas? | REQ-EVO-01 a REQ-EVO-52 y REQ-POB-01 a REQ-POB-03 | ACE-11, ACE-13, ACE-26, ACE-88 a ACE-90, ACE-104 y ACE-105 |
| Evolución | ¿Cómo cambia el estado poblacional entre generaciones? | REQ-EVO-53 a REQ-EVO-78, REQ-REP-01 a REQ-REP-26, REQ-MIG-01 a REQ-MIG-09, REQ-SEL-01 a REQ-SEL-17 y REQ-REC-01 a REQ-REC-29 | ACE-13, ACE-24 a ACE-28, ACE-76 a ACE-87 y ACE-117 |
| Familias y cruces | ¿Qué individuos intermedios y descendientes forman cada familia? | REQ-CRU-01 a REQ-CRU-45, junto con las reglas aplicables de REQ-REP-01 a REQ-REP-26 y REQ-REC-01 a REQ-REC-29 | ACE-14, ACE-27, ACE-75, ACE-85, ACE-87, ACE-97, ACE-98 y ACE-106 |
| Observación técnica | ¿Qué evidencia se genera a partir de una verdad ya fijada? | REQ-TEC-01 a REQ-TEC-10, REQ-RG-01 a REQ-RG-07, REQ-DAT-01 a REQ-DAT-73, REQ-WGS-01 a REQ-WGS-14 y REQ-CONT-01 a REQ-CONT-16 | ACE-15 a ACE-17, ACE-37 a ACE-59, ACE-68 a ACE-74, ACE-93, ACE-94 y ACE-100 |
| Verdad interoperable | ¿Cómo se proyecta la verdad canónica a cada referencia y formato solicitado? | REQ-VER-01 a REQ-VER-47 | ACE-18 a ACE-21, ACE-24, ACE-25, ACE-30, ACE-31, ACE-60 a ACE-63, ACE-66, ACE-100, ACE-101 y ACE-108 |
| Interfaz de autoría | ¿Puede configurarse y lanzarse la misma semántica sin editar TOML? | REQ-UI-01 a REQ-UI-15 y REQ-VAL-54 a REQ-VAL-56 | ACE-114 a ACE-117 |
| Cierre | ¿Cómo se decide que los resultados están completos y trazables? | REQ-VAL-01 a REQ-VAL-56, REQ-MAN-01 a REQ-MAN-92 y REQ-SAL-01 a REQ-SAL-27 | ACE-22 a ACE-24, ACE-33, ACE-34, ACE-36, ACE-46 a ACE-49, ACE-60 a ACE-67, ACE-78, ACE-86, ACE-92 a ACE-117 |

La relación conceptual entre etapas es: primero se admite una petición, luego se
construye o recupera la verdad biológica, después se materializa la evidencia
técnica solicitada y, por último, se valida y publica el conjunto. Esta frase
solo facilita la lectura; los órdenes, condiciones, cardinalidades, errores y
excepciones aplicables son los que define la sección 4.

La correspondencia completa requisito por requisito con criterio, hito y prueba
planificada se mantiene en
[trazabilidad-requisitos-0.1.md](trazabilidad-requisitos-0.1.md).

## 6. Mapa de entradas y salidas

Esta sección también es **derivada e informativa**. Agrupa la interfaz por la
pregunta científica que resuelve, sin volver a definir campos, formatos,
valores, obligatoriedad ni relaciones entre objetos.

### 6.1. Entradas

| Grupo de entrada | Finalidad conceptual | Contrato aplicable |
|---|---|---|
| Petición y configuración | Identificar la ejecución y resolver los parámetros efectivos | REQ-CONF-01 a REQ-CONF-130, REQ-FMT-01 a REQ-FMT-07 y REQ-NOF-01 a REQ-NOF-08 |
| Secuencias y referencias | Establecer las bases y coordenadas sobre las que se construye y compara la verdad | REQ-REF-01 a REQ-REF-36 y las entradas de REQ-DAT-01 a REQ-DAT-73 |
| Individuos, poblaciones y familias | Declarar el material biológico inicial y las relaciones que deben simularse | REQ-IND-01 a REQ-IND-65, REQ-EVO-01 a REQ-EVO-78, REQ-POB-01 a REQ-POB-03 y REQ-CRU-01 a REQ-CRU-45 |
| Procesos biológicos | Especificar creación, mutación, migración, reproducción, selección y recombinación | REQ-DIR-01 a REQ-DIR-05, REQ-MUT-01 a REQ-MUT-15, REQ-MIG-01 a REQ-MIG-09, REQ-REP-01 a REQ-REP-26, REQ-SEL-01 a REQ-SEL-17 y REQ-REC-01 a REQ-REC-29 |
| Eventos y condiciones técnicas | Solicitar observaciones, perfiles, artefactos y tratamientos técnicos admitidos | REQ-TEC-01 a REQ-TEC-10, REQ-RG-01 a REQ-RG-07, REQ-WGS-01 a REQ-WGS-14 y REQ-CONT-01 a REQ-CONT-16 |
| Aplicación Marimo | Importar, editar, revisar, validar y lanzar la misma configuración mediante una interfaz local | REQ-UI-01 a REQ-UI-15 |

La asignación exacta de cada entrada a sus ACE se consulta en la matriz de
trazabilidad; esta agrupación no sustituye esa correspondencia individual.

### 6.2. Salidas

| Grupo de salida | Finalidad conceptual | Contrato aplicable |
|---|---|---|
| Verdad biológica y genealogía | Conservar el estado canónico y la procedencia necesarios para interpretar la simulación | REQ-VER-01 a REQ-VER-47 y los requisitos de procedencia de REQ-IND-01 a REQ-IND-65, REQ-CRU-01 a REQ-CRU-45 y REQ-REC-01 a REQ-REC-29 |
| Evidencia técnica | Entregar los artefactos solicitados sin convertirlos en fuente de verdad | REQ-TEC-01 a REQ-TEC-10, REQ-RG-01 a REQ-RG-07, REQ-DAT-01 a REQ-DAT-73, REQ-WGS-01 a REQ-WGS-14 y REQ-CONT-01 a REQ-CONT-16 |
| Validación e informes | Hacer observables las comprobaciones deterministas y los resultados estadísticos | REQ-VAL-01 a REQ-VAL-56 |
| Manifiestos y paquete publicado | Relacionar configuración efectiva, verdad, evidencia, procedencia y estado de cada artefacto | REQ-MAN-01 a REQ-MAN-92 y REQ-SAL-01 a REQ-SAL-27 |

La matriz de trazabilidad asigna a cada salida sus ACE concretos. Los requisitos
citados determinan qué salida existe en cada caso y qué contenido debe
conservar; esta clasificación no añade una salida ni permite omitir otra.

## 7. Criterios de aceptación preliminares

Un criterio de aceptación es una comprobación observable que permite decidir si
un requisito se ha cumplido. Los siguientes criterios pueden deducirse de los
requisitos confirmados. Sus casos de prueba y detalles concretos siguen
pendientes.

| Criterio | Evidencia preliminar | Aspectos pendientes |
|---|---|---|
| **ACE-01. Configuración sin cambiar el código** | El programa puede leer ficheros generales con distinto número, tamaño o composición de cromosomas o *scaffolds* y generar las referencias simuladas correspondientes sin modificar el código. Cuando se selecciona la vía evolutiva o la vía de familias y poblaciones segregantes, el fichero general referencia el fichero de configuración poblacional del escenario y el programa usa ambos | Formato, esquema y casos de prueba |
| **ACE-02. Tamaños del genoma simulado** | Cada cromosoma o *scaffold* tiene el tamaño configurado y el tamaño total corresponde a los tamaños de todas sus secuencias | Casos de prueba y comprobación exacta |
| **ACE-03. Modos de referencia inicial** | El sistema puede simular la referencia inicial desde la configuración o utilizar un FASTA real suministrado. Un FASTA legible, con al menos una secuencia no vacía y cuyo primer elemento no vacío de cada cabecera sea un identificador único supera la validación si sus secuencias contienen exclusivamente A, C, G, T y N en mayúsculas o minúsculas. `N` se conserva exactamente, ocupa coordenada y cuenta en la longitud; nunca se resuelve. La descripción restante de la cabecera se conserva por separado. El sistema copia el FASTA al directorio de resultados, usa esa copia y deja el original sin cambios. Un archivo ilegible, sin secuencias no vacías, con identificadores vacíos o repetidos, huecos, otros códigos IUPAC o cualquier otro símbolo hace que la ejecución se detenga antes de generar resultados y explique los errores | Casos de alfabeto, integridad y conservación de `N` |
| **ACE-04. Un FASTA por referencia** | Cada referencia inicial o derivada conserva toda su secuencia en un único FASTA, que contiene uno o varios cromosomas o *scaffolds* | Casos de prueba |
| **ACE-05. Índices de las referencias** | Una ejecución normal v1 conserva cero índices de referencia. Una operación indexada puede crear bajo demanda un único FAI, consumirlo y validarlo sin publicarlo como salida normal. Los índices múltiples quedan en DTF-048 | Caso sin índice, caso con FAI y validación independiente |
| **ACE-06. Creación de referencias estructurales** | Se pueden crear dos o más referencias estructurales a partir de referencias disponibles previamente, incluida una referencia inicial suministrada | Casos de prueba, tamaños y reglas exactas |
| **ACE-07. Diferencias entre referencias** | Cada nueva referencia estructural difiere de su referencia de origen solo por las duplicaciones, deleciones, inversiones y translocaciones solicitadas | Casos de prueba y comprobación de cada cambio |
| **ACE-08. Origen de las referencias** | Para cada referencia estructural se pueden conocer su referencia de origen y los cambios que las distinguen | Forma de representar y consultar la relación |
| **ACE-09. Clases de variante** | Se pueden crear individuos con SNP, *indels*, variantes estructurales y repeticiones en tándem respecto a cualquiera de las referencias disponibles | Casos de prueba por clase |
| **ACE-10. Creación de cambios** | Los SNV heredables iniciales proceden de una matriz poblacional suministrada o generada; su aleatoriedad está en las frecuencias y genotipos, no en una densidad independiente por individuo. *Indels*, variantes estructurales y mutaciones *de novo* pueden indicar cantidad media o densidad; la densidad se convierte usando longitud válida y el recuento usa Poisson. Los valores espaciales e individuales se resuelven mediante ACE-111. Una sustitución individual solo admite entidades preexistentes, no se hereda ni se aplica retroactivamente, y un valor alélico fijado es una restricción obligatoria | Matriz SNV y corpus poblacional; coordenadas; pruebas estadísticas de Poisson y de frecuencias; reproducibilidad |
| **ACE-11. Ploidía** | La configuración admite un valor global y sustituciones por cromosoma o *scaffold*, cada uno entero de 1 a 10. Toda secuencia tiene un valor efectivo idéntico en todos los individuos. Los valores no enteros, fuera de rango o ausentes se rechazan antes de simular. Estado genotípico y VCF representan exactamente sus haplotipos y dosis | Casos de prueba de límites y sustituciones |
| **ACE-12. Creación directa** | Se pueden crear directamente individuos con cambios concretos o aleatorios indicados en la configuración | Casos de prueba |
| **ACE-13. Vía evolutiva** | Cada escenario usa fases ordenadas y generaciones discretas no solapadas. N objetivo fija los intentos reproductivos y N realizado registra los descendientes creados; una desviación aleatoria advierte sin detener. Ne objetivo se rechaza y Ne realizado se calcula como Ne reproductivo cuando es evaluable. El orden único es demografía/N objetivo, migración, viabilidad, elección parental/reproducción, segregación/recombinación, mutación de novo, cálculo y registro de N/Ne/frecuencias, y muestreo/salidas. La deriva surge del muestreo finito | Casos de las ocho etapas, etapa opcional omitida y registro efectivo |
| **ACE-14. Familias y poblaciones segregantes** | Cada escenario representa los cruces como un grafo dirigido acíclico con fuentes, orden y número de descendientes. Un ciclo es una dependencia imposible de un resultado posterior, no reutilizar un parental: autofecundación, retrocruce y ramas múltiples son válidos cuando cada evento produce descendientes nuevos. Antes de simular se validan identificadores, existencia y precedencia de referencias, creación única de descendientes, cardinalidad parental, ploidía par e igual, identidad estructural canónica y producibilidad de reglas y salidas. Existen plantillas F1, F2, retrocruce, autofecundación y pedigrí simple; los demás diseños se expresan como cruces personalizados. Un parental expreso declara sus cambios por copia respecto a una referencia compatible; lo no declarado permanece como referencia y una fase incompleta se rechaza. La descendencia alimenta eventos posteriores por identificador, selección aleatoria, pertenencia o estado genético en uno o varios loci, combinables con Y, O y NO. Un número exacto usa sorteo sin reemplazo si sobran elegibles y falla si faltan; una proporción es una probabilidad individual y produce un recuento aleatorio reproducible. Todo se evalúa sobre la verdad biológica. Un individuo puede alimentar varias ramas con la misma identidad y genotipo, mientras cada evento forma gametos independientes y descendientes nuevos. Se materializan el estado faseado y los metadatos de todos los individuos y poblaciones intermedias, y el manifiesto los enumera para reconstruir eventos y linajes | Familia completa, sintaxis del grafo, parentales, condiciones y casos reproductivos |
| **ACE-15. Referencia de mapeo** | Cuando están activos BAM o CRAM, la configuración selecciona una o varias referencias de mapeo, incluida una distinta de la estructura de una muestra. Para cada referencia se generan los archivos de todas las muestras incluidas en un conjunto independiente | Cardinalidad de archivos dentro de cada conjunto y comprobaciones técnicas |
| **ACE-16. Cabecera de BAM y CRAM** | La cabecera y su grupo de lectura permiten identificar sin ambigüedad referencia, muestra, evento, biblioteca, tanda y grupo. `SM` conserva el identificador biológico estable y cada identificador `RG` es único en toda la ejecución. El archivo cumple la versión oficial elegida | Versiones exactas de SAM/BAM y CRAM, catálogo concreto de campos y casos de prueba |
| **ACE-17. Degradaciones** | En F4, un único objetivo de cobertura común puede ser bajo o cero; no admite sustituciones individuales, por secuencia ni regionales. Antes de sesgos, la WGS básica realizada es estocástica y aproximadamente Poisson, no exacta por base. `N` permanece elegible, cuenta en cobertura, aparece con Phred 0 y no la resuelve el modelo de error. El contenido GC excluye `N` y una ventana sin bases conocidas no recibe corrección. Se registran por separado objetivo, WGS básica y FASTQ final con bases propias y contaminantes | Intensidad, métricas y comprobación de las demás degradaciones; casos de `N` |
| **ACE-18. VCF de verdad** | Cuando se solicita VCF, existe uno independiente por referencia de mapeo con los genotipos y dosis de las mismas muestras biológicas, expresados contra esa referencia. Cada VCF se construye directamente, usa identificadores de variante locales y no expone correspondencias registro a registro con otros VCF | Casos por referencia y comprobación de espacios de identificadores independientes |
| **ACE-19. Cabecera del VCF** | La cabecera o los metadatos estándar de cada VCF permiten identificar sin ambigüedad la referencia, el conjunto y las muestras mediante sus identificadores biológicos estables, y cumplen la versión oficial elegida | Versión exacta de VCF, catálogo concreto de campos y casos de prueba |
| **ACE-20. Separación por referencia** | Cuando se generan BAM, CRAM o VCF, existe un conjunto independiente por cada referencia de mapeo seleccionada. Cada conjunto usa una sola referencia, todos contienen las mismas muestras solicitadas y ninguno mezcla archivos referidos a genomas distintos. En modo realista, dentro del evento todos usan los mismos FASTQ y el único mapeador del evento | Forma de organizar los conjuntos |
| **ACE-21. Coherencia de las salidas** | Las referencias, los individuos y cada conjunto de BAM o CRAM y VCF de verdad se pueden relacionar sin ambigüedad y describen la misma simulación | Identificadores y comprobaciones de consistencia |
| **ACE-22. Manifiestos** | Cada ejecución genera un manifiesto general y uno por evento explícito. Identifican vía, configuraciones, semillas y subsemillas, versión, parámetros efectivos, referencias y recursos comunes, individuos y todos los archivos. Cada evento registra sus salidas y, cuando proceda, modo, tecnología, perfil, contaminación, FASTQ, mapeador único, versión y comando. También conservan evolución, cruces, ploidía, validaciones, los niveles de reproducibilidad declarados y el estado y diagnóstico de cada unidad | Formatos, esquemas y casos de prueba |
| **ACE-23. Niveles de reproducibilidad** | Toda ejecución usa una semilla. Si la configuración no la proporciona, el programa la genera antes de simular y el manifiesto la identifica como generada. Con iguales semilla, configuración efectiva, entradas, perfiles y versión, se garantiza la misma simulación biológica y se aíslan sus flujos de los eventos técnicos. Las mismas lecturas, calidades, procedencias y artefactos se garantizan además cuando coinciden la subsemilla y las dependencias técnicas aplicables. La equivalencia semántica se declara solo bajo condiciones verificadas y la identidad binaria no se promete de forma general. El manifiesto separa los cuatro niveles y nunca infiere uno superior desde otro inferior | Estrategia técnica de subsemillas; condiciones verificables por herramienta, versión o plataforma; casos de prueba |
| **ACE-24. Muestreos evolutivos** | La configuración puede solicitar muestreos mediante generación y número o proporción. Se conservan por defecto el estado genotípico y los metadatos mínimos, sin BAM, CRAM ni VCF. Cada evento guarda por separado un resumen de la población completa y subpoblaciones relevantes y otro calculado solo desde la muestra. Ambos identifican ámbito e incluyen frecuencias alélicas, heterocigosidad, F calculado solo en diploides o «no aplicable», N y Ne cuando corresponden; ninguna medida muestral se presenta como poblacional. El muestreo no modifica individuos ni la trayectoria evolutiva | Formato de resúmenes; estrategia técnica de flujos aleatorios; casos de prueba |
| **ACE-25. Identidad individual** | Todo individuo tiene un identificador único, estable y no reutilizable dentro de la ejecución. Los duplicados aportados por la configuración se rechazan. El mismo identificador relaciona genotipo, generación y subpoblación de origen, parentales o linaje, muestreos, BAM o CRAM, VCF y manifiesto, y se distingue de los identificadores de evento y conjunto | Sintaxis técnica de los identificadores y casos de prueba |
| **ACE-26. Estructura poblacional dinámica** | División reparte progenitores entre hijas exclusivas; fusión reúne orígenes; fundación transfiere fundadores y deja continuar al origen. División y fundación usan el catálogo declarativo común. Cada nueva subpoblación tiene N objetivo independiente y registra conjunto fundador, N objetivo y realizado. En fusión, el modo probabilístico usa elegibles o pesos; el controlado garantiza recuentos parentales exactos. Ninguno garantiza porcentaje exacto de bases: se calcula ascendencia genómica realizada con traza haplotípica | Casos de prueba del catálogo común y de cada evento |
| **ACE-27. Recombinación** | La vía evolutiva y la de familias simulan recombinación al formar descendientes. Toda la ejecución usa una tasa uniforme o mapas completos, sin interferencia. Cada individuo pertenece a cero o a un grupo configurable y usa su perfil específico o el general. Solo se forman bivalentes entre copias con estado estructural canónico idéntico; una diferencia estructural invalida al individuo o pareja completos. En ploidías pares 2–10, las copias se emparejan al azar y cada parental aporta la mitad. Los mapas no válidos y las pertenencias múltiples se rechazan | Casos diploide/tetraploide idénticos; variantes pequeñas permitidas; SV heterocigota; ploidía impar/desigual |
| **ACE-28. Fase y haplotipos** | El estado interno separa todas las copias cromosómicas según la ploidía, asigna cada variante a haplotipos concretos y conserva fase y origen de segmentos durante evolución, recombinación, cruces, muestreos y poblaciones intermedias. Los VCF emiten genotipos faseados cuando el registro puede representarlos; cualquier limitación usa la única entrada complementaria del caso para conservar la fase sin pérdida | Estructura de datos, algoritmo y casos que requieren complemento |
| **ACE-29. Umbral entre *indels* y variantes estructurales** | Sin configuración se usan 50 pb: una inserción o deleción menor de 50 pb es *indel* y una igual o mayor es variante estructural. Otro umbral configurado se aplica de forma coherente a generación, estadísticas, VCF y manifiesto. Otros tipos estructurales conservan su clase | Casos límite y pruebas de coherencia |
| **ACE-30. Representación VCF** | SNP e *indels* usan REF/ALT explícitos y normalizados respecto a la referencia del conjunto. Nunca se emite un SNP con `REF=N`. Otros cambios conservan `N` en REF/ALT solo si el estándar elegido lo permite sin pérdida; en otro caso usan complemento sin inventar una base. Varios ALT representables del mismo locus forman un único registro multialélico cuyo `GT` conserva ploidía, dosis y fase. Los solapamientos pequeños se representan desde la secuencia haplotípica final. Las variantes estructurales usan la representación estándar más completa y campos estándar con su semántica oficial. El VCF cumple estrictamente la versión oficial seleccionada y sigue siendo válido e interoperable aunque un caso requiera complemento | Versión exacta VCF, normalización técnica y casos complejos con `N` |
| **ACE-31. Repeticiones en tándem** | La configuración admite loci concretos o generación aleatoria y regional con motivo no vacío de longitud libre, copias o longitud, pureza e interrupciones; un motivo con `N` se rechaza o excluye. Se mantienen fase y verdad en secuencias y BAM o CRAM. El VCF usa alelos explícitos cuando son manejables o simbólicos para tramos largos o complejos y expone cuanto admita sin pérdida el estándar; el resto se conserva en la única entrada complementaria del caso y nunca solo en el manifiesto | Campos VCF concretos, criterio técnico de tamaño manejable, límites de recursos y casos de prueba |
| **ACE-32. Mutación de repeticiones en tándem** | El único modelo de v1 es el escalonado simple: cada mutación añade exactamente una unidad completa o elimina exactamente una. La tasa y el sesgo o probabilidad entre expansión y contracción son configurables en los ámbitos acordados; el tamaño del salto no lo es. Las mutaciones internas de interrupciones, pureza o motivo están desactivadas por defecto y solo se aplican si se configuran. Cada evento es reproducible, faseado y trazable | Casos de expansión, contracción, sesgo, ámbitos y rechazo de saltos múltiples |
| **ACE-33. Verdad complementaria sin pérdida** | Si un caso complejo no puede representarse completamente en VCF, se genera exactamente una entrada complementaria con todo el detalle perdido. Uno o varios registros del caso la enlazan mediante un identificador local estable y una descripción breve en campos `INFO` declarados. Cada registro enlazado pertenece a una sola entrada y caso. A partir del VCF y el complemento se reconstruyen sin pérdida secuencia o cambio, coordenadas, haplotipo, fase y genotipo. Si todos los casos son representables, el complemento no se genera | Formato, esquema y casos de reconstrucción |
| **ACE-34. Resultado mínimo y salidas dependientes** | Toda ejecución exitosa produce FASTA e índices necesarios, la verdad completa obligatoria para su vía, los historiales o resúmenes exigidos, copias de configuraciones y manifiesto. Toda salida de observación requiere un evento explícito y cero eventos no crea ninguno implícito. Un evento puede solicitar solo verdad o VCF sin parámetros técnicos. Cuando solicita alineamientos, CRAM es predeterminado y BAM la alternativa; son obligatorios sus FASTQ finales y sus resúmenes técnicos agregados. El modo directo queda aplazado en DTF-051. V1 no produce perfiles empíricos. La ejecución falla si falta una salida mínima o activada | Formato del fichero de verdad, organización de copias y casos de fallo |
| **ACE-35. Escala y prioridad de la primera versión** | Se pueden configurar, sin límites funcionales codificados, distintos números y tamaños de secuencias, individuos y generaciones. La validación prioriza escenarios pequeños y reproducibles. Los escenarios mayores están limitados por los recursos disponibles y no se presupone rendimiento ilimitado | Casos pequeños prioritarios, mediciones de recursos y pruebas de escalabilidad sin umbrales fijados |
| **ACE-36. Estimación de recursos** | Antes de ejecutar se estiman volumen, almacenamiento, memoria y tiempo incluyendo la trazabilidad obligatoria y todas las salidas. Un escenario grande solo advierte, pero una configuración inválida, un mínimo imposible o un presupuesto configurado superado impiden iniciar y explican la causa. Nunca se degrada el detalle ni se elige otro nivel; la persona usuaria debe reducir materializaciones, muestreos, individuos, lecturas u otras salidas. El manifiesto conserva estimaciones, advertencias y consumo medible | Caso suficiente, preflight insuficiente y comprobación de ausencia de degradación |
| **ACE-37. Tecnologías de lectura** | F4 representa exclusivamente lecturas cortas Illumina *paired-end* con un perfil incorporado, identificado, versionado, validado y ajustable solo mediante las sustituciones permitidas. En modo realista las lecturas proceden de haplotipos faseados y resultado y manifiesto registran tecnología y perfil efectivos. Illumina *single-end* se rechaza hasta DTF-057; PacBio HiFi, Oxford Nanopore y PacBio CLR se rechazan hasta las ampliaciones indicadas | Perfil PE, casos de rechazo y datos de validación |
| **ACE-38. Alineamientos realistas** | Cuando se solicitan, los BAM o CRAM se obtienen exclusivamente al mapear los FASTQ finales después de aplicar errores, contaminación y demultiplexado. Conservan pares, identificadores y registros no mapeados producidos por el mapeador. El modo directo queda aplazado en DTF-051 | Organización y casos de prueba del recorrido realista |
| **ACE-39. Un mapeador por evento técnico** | Un evento con alineamientos usa exactamente un mapeador, minimap2 salvo selección compatible distinta, contra todas sus referencias. Para minimap2, el perfil elige automáticamente el *preset* y los valores predeterminados, y las sustituciones se advierten y validan. El manifiesto conserva herramienta, versión, comando completo y parámetros efectivos | Versión de minimap2, *presets* y valores exactos, catálogo inicial de alternativas, comprobaciones técnicas de compatibilidad y casos de opciones incompatibles |
| **ACE-40. Independencia de eventos técnicos** | Dos eventos pueden usar mapeadores o tecnologías distintos sobre las mismas muestras biológicas, pero cada uno genera sus propias lecturas y usa un solo mapeador. No se reutilizan FASTQ entre eventos ni se presentan como comparación de mapeadores | Organización, identificadores y casos con varios eventos |
| **ACE-41. Cobertura configurable** | F4 declara una única cobertura media objetivo en «x» por evento, incluidos valores altos, bajos o cero. Su denominador es la suma de cada secuencia de la referencia estructural una sola vez, sin ploidía; los haplotipos realizados no lo cambian. Las sustituciones por individuo, secuencia o región se rechazan. La WGS básica es aproximadamente Poisson antes de sesgos. La cobertura final del FASTQ separa bases propias y contaminantes | Métricas y casos estadísticos por las tres medidas F4 |
| **ACE-42. Contaminación trazable** | El rango predeterminado de contaminación es `0–0 %`; al sustituirlo, por defecto cada individuo de la población seleccionada es receptor, aunque la configuración puede restringirlo a un subconjunto. Cada receptora seleccionada sortea su fracción total dentro del rango configurable 0–100 % y la convierte a parejas como el entero más próximo de `fracción × parejas limpias`, con empate hacia arriba. Sus unidades contaminantes se sortean sin reemplazo sobre todas las parejas limpias elegibles de los otros individuos de su misma población; no existe matriz, porcentaje ni pesos por donante. Cada unidad contaminante sustituye una propia; la donante no pierde datos. El manifiesto conserva rango, receptoras, objetivo sorteado, conversión, donantes, realizado, retiradas y totales finales | Validaciones y casos de prueba |
| **ACE-43. Contaminación *paired-end*** | En datos *paired-end* se seleccionan y sustituyen parejas completas de una misma fuente. Ambos *mates* se incorporan juntos con nombres, indicadores, orientación y relación coherentes, sin crear huérfanos; el manifiesto registra parejas propias retiradas y contaminantes incorporadas | Casos de prueba de emparejamiento |
| **ACE-44. FASTQ contaminante externo** | F4 rechaza FASTQ externos como contaminación antes de generar y explica que la capacidad corresponde a DTF-052. No lee, valida, copia ni registra esas entradas | Rechazo previo y ausencia de artefactos |
| **ACE-45. Contaminación técnica por índices** | Tras la contaminación biológica, una tanda multiplexada usa UDI obligatorios y simula de forma reproducible *index hopping*: cada pareja de la instantánea post-contaminación se evalúa exactamente una vez y una fracción se reasigna uniformemente a otra biblioteca. Las transferencias son simultáneas; una pareja reasignada desaparece de la donante y entra en la receptora sin duplicarse, incluso si ya era contaminante, sin cadenas ni dependencia del orden. F4 no genera UDI inesperados; si uno alcanza el demultiplexado se descarta defensivamente. No se generan señales, BCL ni lecturas i5/i7 detalladas | Casos de transferencia simultánea, permutación de bibliotecas y resultado |
| **ACE-46. FASTQ final y resumen técnico** | Existe un FASTQ final, o pareja *paired-end*, por individuo y evento después de contaminación y demultiplexado. Corresponde a la única biblioteca del individuo y a la única tanda, sin concatenación ni copia limpia adicional. El manifiesto conserva semilla, perfil, cobertura, contaminación y aportes por donante, quimeras, duplicados, descartes, hashes y tamaños como recuentos agregados. No conserva origen, diferencias, secuencias, calidades ni identificadores por lectura o pareja | FASTQ final, resumen agregado y ausencia de fichero individual |
| **ACE-47. Perfiles de error y calidades** | F4 usa el perfil identificado y versionado `gps-hiseqx-pcrfree-pe150-v1`, homogéneo para el evento salvo desplazamiento de calidad y escala de error globales por mate. Las calidades FASTQ son probabilísticamente coherentes con los errores; selectores por individuo, secuencia o región se rechazan. La salida conserva perfil, ajustes efectivos y recuentos agregados, no diferencias por lectura | Perfil literal, coherencia calidad-error, rechazos de selectores y casos estadísticos |
| **ACE-48. Calibración empírica fuera de v1** | Una configuración que aporte FASTQ o referencia para calibración, solicite estimar o inferir un perfil o seleccione un perfil empírico reutilizable se rechaza antes de simular. El error identifica el campo y explica que la capacidad queda fuera de v1. Una ejecución válida con un perfil incorporado de Illumina y ajustes explícitos permitidos no genera un fichero de perfil empírico ni lo enlaza desde el manifiesto; PacBio HiFi permanece aplazado en DTF-049 | Casos por clase de campo rechazado y comprobación de ausencia del artefacto |
| **ACE-49. Duplicados de PCR y ópticos** | F4 declara PCR no aplicable y duplicados ópticos a cero; no genera duplicados ni admite sus tasas. Un perfil amplificado y el modelo de duplicados quedan en DTF-056 | Evidencia, modelo y marcado posterior |
| **ACE-50. Restos de adaptadores** | El perfil HiSeq X inicial no produce adaptadores: todos sus insertos son mayores que sus lecturas. F4 rechaza activarlos y conserva recuento cero. Un perfil futuro de fragmentos cortos podrá generarlos únicamente por sobrepasar el inserto, conforme a DTF-054 | Rechazo y recuento cero en F4; perfil futuro |
| **ACE-51. Sesgo de cobertura por GC** | F4 usa GC neutro y rechaza una curva o intensidad no nulas. El perfil con sesgo GC calibrado queda en DTF-056 | Curvas, regla de combinación e intervalos de resumen concretos |
| **ACE-52. Baja mapeabilidad realista** | Las lecturas de regiones repetitivas, duplicadas u homólogas parten de posiciones y haplotipos verdaderos. En F4 se publican solo FASTQ y métricas agregadas; el mapeador, sus alineamientos y la evaluación de ambigüedad o posiciones incorrectas quedan para DTF-055, sin sobrescritura artificial de MAPQ | Casos de prueba por tipo de región y referencia cuando exista mapeo |
| **ACE-53. Alineamientos sintéticos coherentes** | Queda aplazado en DTF-051 el modo directo, incluidas sus categorías, MAPQ, campos y validación propios | Decisión posterior de modo directo |
| **ACE-54. Lecturas quiméricas** | Una tasa configurable transforma parejas limpias en quimeras sin cambiar su número. Para inserto `L`, el primer segmento es el prefijo de longitud `b`, con `b` uniforme en `1..L-1`; el segundo es un intervalo uniforme de longitud `L-b` de la lista canónica del mismo individuo y biblioteca. En la misma copia se excluyen intervalos que solapen o toquen el prefijo; ambos segmentos se leen en orientación directa antes de aplicar PE. Una falta de candidata detiene antes de publicar. F4 publica FASTQ y recuentos agregados, no procedencia individual ni alineamientos | Vector literal de fuentes, uniones, exclusiones y fallo sin candidata |
| **ACE-55. Anomalías *paired-end*** | Aplazadas en DTF-055: orientación anómala, inserto observado inesperado y mate no alineado requieren BAM/CRAM y mapeo. F4 solo publica pares FASTQ normales y las quimeras explícitas de ACE-54 | Modelo, métricas y casos de prueba cuando exista mapeo |
| **ACE-56. Grupo de lectura y archivo por individuo** | Cada individuo tiene exactamente una biblioteca, un grupo de lectura y, en modo realista, un FASTQ o pareja dentro del evento. Los identificadores de muestra, biblioteca, tanda, plataforma y grupo son consistentes entre FASTQ, manifiesto y, desde F5.6, `@RG`/`RG` de BAM o CRAM. F4 no publica procedencia individual; los vínculos transitorios solo permiten validar esas salidas. Cada conjunto futuro por referencia contiene un único BAM o CRAM final por individuo con ese grupo y nunca mezcla referencias | Cabeceras FASTQ, futuros campos `@RG`, validación por reproducción y casos de cardinalidad inválida |
| **ACE-57. Una tecnología y una tanda por experimento** | Cada evento que genere FASTQ, BAM o CRAM constituye su único experimento de secuenciación, sin identidad separada, y contiene una tecnología, un perfil, exactamente una tanda y una biblioteca por individuo. Todas las bibliotecas se procesan conjuntamente y, si existen alineamientos, comparten mapeador. Otra biblioteca del individuo, tanda, perfil, tecnología o control exige otro evento. Los carriles se rechazan en v1. Un evento que solo materializa verdad o VCF no es un experimento y no admite campos técnicos | Casos válidos y rechazos por cardinalidad, carril o campo técnico indebido |
| **ACE-58. CRAM predeterminado y referencia exacta** | Dentro de un evento que solicite alineamientos, cada individuo recibe un CRAM final salvo que se seleccione BAM como alternativa. Cada CRAM está ligado por cabecera y manifiesto a la copia FASTA exacta almacenada entre los recursos comunes, protegida por identificadores y comprobaciones de integridad; no puede asociarse a otra referencia. Una exportación del evento sí copia ese FASTA al paquete autocontenido | Campos, comprobaciones de integridad y casos de referencia incorrecta |
| **ACE-59. Orden, índice y validación** | Cada BAM o CRAM final se entrega ordenado por coordenadas de forma predeterminada y con un índice compatible elegido según formato, tamaño y necesidad. La configuración y el manifiesto registran formato, orden e índice, y el resultado no se considera completo si alineamiento, índice y referencia no superan la validación conjunta | Herramientas, tipos exactos de índice y casos límite de tamaño |
| **ACE-60. VCF multimuesta por referencia** | Cada conjunto solicitado contiene exactamente un VCF de verdad multimuesta con todos sus individuos y los mismos identificadores usados en BAM/CRAM, grupos de lectura, genotipos y manifiesto. No incluye individuos de otra referencia y no genera VCF individuales ni archivos separados por clase o complejidad | Casos de coherencia de identificadores y organización |
| **ACE-61. VCF ordenado, comprimido e indexado** | Por defecto, el VCF sigue el orden de secuencias de la referencia y las coordenadas dentro de cada una, está comprimido con BGZF y tiene un índice TBI o CSI compatible. Puede solicitarse texto sin comprimir. La ejecución falla si VCF, compresión, índice, cabecera y referencia no superan la validación correspondiente | Selección TBI/CSI, herramientas y casos límite |
| **ACE-62. Solo loci diferentes de la referencia en VCF** | El VCF contiene un registro cuando al menos una muestra del conjunto tiene un alelo no de referencia, aunque ese ALT esté fijado en todas las muestras. No emite posiciones idénticas a la referencia ni bloques gVCF por defecto; el fichero estructurado de genotipos y haplotipos conserva la verdad completa | Casos invariantes de referencia, ALT fijado, loci polimórficos, multialélicos y distintas ploidías |
| **ACE-63. VCF independiente de la evidencia observada** | El VCF se construye desde el estado genotípico faseado de la simulación, nunca mediante inferencia desde FASTQ, BAM o CRAM. En cada locus incluido conserva el genotipo completo y la fase real de todas las muestras aunque tengan cobertura cero, contaminación, errores, calidad baja o mapeo ausente, incorrecto, múltiple o ambiguo; esos casos no producen genotipos ausentes. Una limitación de representación se conserva mediante el complemento de verdad y no se confunde con falta de evidencia | Casos con cobertura cero y degradaciones combinadas; comprobación de independencia respecto a FASTQ y BAM/CRAM |
| **ACE-64. Validación interna contra la verdad** | Siempre que existen FASTQ o BAM/CRAM, existe un resultado estructurado enlazado desde el manifiesto. Sin ejecutar un *SNP caller*, comprueba contra la verdad las procedencias, parejas, alineamientos, referencias e identidades, e informa cobertura, soporte alélico, errores, contaminación, duplicados y demás magnitudes aleatorias. Cada resultado estadístico registra ámbito, esperado, realizado, tamaño observado, desviaciones absoluta y relativa y explicación; si no puede calcularse, figura como «no evaluable» con causa | Métodos estadísticos, formato, esquema y casos de medida no evaluable |
| **ACE-65. Consecuencias de la validación** | Un archivo, formato, índice o cabecera inválidos, una referencia o identificadores incompatibles, un genotipo imposible, una procedencia o haplotipo imposible o un artefacto sin explicación trazada detienen la ejecución e impiden completar el ámbito afectado. Una cobertura u otra magnitud aleatoria desviada se informa y permite completar, sin tolerancia de aceptación ni opción de elevarla a error. N desviado conserva un aviso no invalidante. Ne realizado no tiene objetivo; N realizado=0 aplica su política específica | Cobertura desviada que completa, procedencia imposible que falla, medida no evaluable y aviso de N |
| **ACE-66. Contenido genotípico del VCF verdadero** | En cada locus incluido, `GT` coincide con el genotipo real simulado y conserva fase y ploidía aunque la cobertura sea cero o las lecturas estén degradadas. El VCF no contiene `DP`, `AD`, `GQ`, `PL` ni ningún otro campo `FORMAT` derivado de observaciones o inferencias. La cobertura, el soporte alélico, las calidades y métricas equivalentes aparecen en los resultados de validación de FASTQ/BAM/CRAM, no en el VCF | Casos con distintas ploidías, cobertura cero y degradaciones; inspección de campos `FORMAT` |
| **ACE-67. QUAL y FILTER del VCF verdadero** | Todos los registros válidos usan `QUAL=.` y `FILTER=PASS`. La cabecera explica que no se estima una calidad probabilística y que `PASS` identifica un registro verdadero validado estructuralmente. Si un registro no puede representarse mediante el VCF y su complemento cuando corresponda, o no puede validarse, se produce un error estructural, se detiene la creación del conjunto, este no se marca como completo y el registro no se emite como `PASS` | Casos de cabecera válida, complemento enlazado, registro no representable y fallo de validación |
| **ACE-68. Alcance WGS y preparación** | F4 acepta exclusivamente WGS de DNA genómico con Illumina PE y preparación sin PCR. Rechaza preparación amplificada, PacBio HiFi, Oxford Nanopore, PacBio CLR, paneles, amplicones, GBS y captura o reducción de representación. La preparación amplificada, duplicados PCR y sesgos asociados quedan en DTF-056 | Casos de rechazo por tecnología y protocolo |
| **ACE-69. Objetivo y muestreo de bibliotecas WGS** | En F4, cada biblioteca declara exclusivamente cobertura; rechaza un número de lecturas o parejas como objetivo y deriva el número entero de parejas PE. La cobertura usa bases de lectura divididas por la longitud haploide de la referencia estructural, sumando cada secuencia una vez; `N` cuenta tanto en longitud como en bases y posiciones muestreadas y se resume por separado. Las posiciones se muestrean uniformemente y de forma reproducible y la profundidad previa a artefactos es aproximadamente Poisson, sin profundidad exacta por base. Mapeo y cobertura alineada se incorporan desde DTF-055 y no redefinen el objetivo | Prueba de conversión desde cobertura, rechazo numérico, casos con `N` y distribución de muestreo |
| **ACE-70. Ploidía, copias y soporte alélico WGS** | Sin excepciones, las lecturas se reparten por longitud y número de copias, cada copia recibe cobertura semejante y la profundidad cromosómica total es proporcional a su ploidía. Los alelos contribuyen según dosis; la validación admite variación aleatoria y no exige balance exacto. No existen parámetros de sesgo alélico o *allele dropout* forzados | Métricas y resúmenes por copia o haplotipo |
| **ACE-71. Un único alineamiento final** | F4 no genera BAM/CRAM ni un evento de control. Tras DTF-050 y DTF-055, cada evento conservará solo su BAM/CRAM final con las degradaciones configuradas; un control sin degradaciones será otro evento con misma verdad, referencias y muestras, pero subsemilla, lecturas y archivos propios, enlazado por manifiestos | Rechazo F4 y futuros identificadores técnicos de eventos relacionados |
| **ACE-72. Contaminación biológica sustitutiva** | Cada lectura o pareja contaminante sustituye una propia de la receptora sin reducir los datos de la donante. Cada receptora sortea su fracción total dentro del rango 0–100 %; el 100 % representa sustitución completa. Un rango inválido o un receptor sin donante elegible con objetivo positivo detiene antes de generar | Unidades exactas y casos límite |
| **ACE-73. Fuentes y reparto de contaminación** | Las fuentes son otros individuos de la misma población, compatibles por construcción. Para cada unidad, el donante se sortea reproduciblemente entre los elegibles; los aportes por donante son el resultado observado. F4 no admite matrices dirigidas, pesos manuales ni FASTQ externos | Casos de reparto, aislamiento poblacional y ausencia de donantes |
| **ACE-74. Transferencia y descarte por índices** | Una lectura o pareja reasignada por índices sale de la donante y entra en la receptora sin duplicarse. F4 la evalúa una sola vez con su origen pre-salto y materializa todas las transferencias simultáneamente. Las descartadas o sin asignar no se reemplazan y reducen cobertura. Informe y manifiesto conservan objetivos, transferencias, descartes, totales y cobertura final | Casos sin cadena, por permutación de bibliotecas y por descarte |
| **ACE-75. Plantillas y sistemas de cruce** | Existen plantillas F1, F2, retrocruce, autofecundación y pedigrí simple. F2 permite autofecundar F1 o cruzar F1. RIL, introgresión y MAGIC se expresan con cruces personalizados. El sistema admite especies autógamas y alógamas y registra modo, parentales y descendientes | Sintaxis del sistema general y casos por plantilla |
| **ACE-76. Apareamiento mixto y contribución parental** | La autofecundación es una probabilidad reproducible de 0–100 % por evento, sustituible por población, subpoblación y fase. Cada intento elige reproduciblemente el tipo y parentales solo entre opciones con ploidía y estructura reproductivas compatibles; fuera de la autofecundación usa dos supervivientes distintos, sin ponderación por parentesco, distancia o genotipo. La incompatibilidad se registra como causa reproductiva, no como selección de viabilidad. Un parental puede reutilizarse sin límite de descendientes en v1 | Tasas realizadas; opción estructural inválida; ausencia de pareja compatible |
| **ACE-77. N, Ne y descendencia** | N objetivo es obligatorio al iniciar y fija el número de intentos reproductivos; cada éxito crea un descendiente. Se procesan N intentos o se registra inequívocamente la imposibilidad común de los restantes. Se conservan éxitos, fallos y causas, N realizado, desviación y aviso. Una imposibilidad previa se rechaza y N=0 aplica su política. Ne objetivo, absoluto o Ne/N, se rechaza antes de simular; Ne realizado usa progenitores con descendencia y desigualdad de contribuciones, nunca heterocigosidad, endogamia o LD, y si no es evaluable conserva la causa | N alcanzado/parcial, imposibilidad previa, N=0, rechazo de Ne objetivo y Ne realizado/no evaluable |
| **ACE-78. Generaciones e historial reproducible** | Las generaciones son discretas y no solapadas. Los individuos parentales solo existen como población activa en su generación y los sustituye la descendencia; un cambio de N no los elimina retrospectivamente. Se materializan estados completos de población inicial, final y generaciones o muestreos solicitados, y solo resúmenes de las demás. La información interna permite reconstruir el linaje detallado solicitado sin exportar antepasados intermedios. No hay niveles de trazabilidad ni reducción por escala. La regeneración biológica exige misma semilla, configuración efectiva, entradas, perfiles y versión | Evolución materializada y no materializada; esquema del historial |
| **ACE-79. Mutaciones evolutivas** | SNP, inserciones, deleciones y cada SV usan tasas independientes por oportunidad válida y generación. Tras recombinación, cada clase cuenta en todas las copias actuales sus candidatos propios: la ploidía y cada duplicación aumentan oportunidades, una deleción las elimina y un ancla o intervalo inicial cuenta una vez aunque el evento necesite destinos o extremos condicionados. Los ámbitos efectivos no se solapan. La media es tasa×oportunidades y el recuento Poisson se sortea una vez antes de actualizar candidatos evento a evento. `N` se excluye de SNP, motivos de repetición y anclas que necesitan base concreta, pero no de intervalos estructurales que la conservan. Cada locus repetido elegible muta una unidad | Casos de ploidía, duplicación, deleción, región solapada, destino condicionado, locus repetido y exclusiones de `N` |
| **ACE-80. Sin semántica de elementos transponibles en v1** | La configuración rechaza catálogos, familias, actividad, movilización, mecanismos de copiar o cortar y pegar, reglas de destino y metadatos específicos de TE. Un FASTA real se trata como bases sin anotación TE | Casos de rechazo y comprobación de ausencia de campos y salidas TE |
| **ACE-81. Clasificación genérica de copias e inserciones** | Copiar una región preexistente, aunque biológicamente corresponda a un TE, se representa como duplicación segmentaria genérica. Añadir una secuencia nueva se representa como inserción genérica. Verdad, VCF y manifiesto no les atribuyen familia, mecanismo ni metadatos TE | Casos de duplicación e inserción sin anotación TE |
| **ACE-82. Selección de viabilidad** | Toda la ejecución usa exactamente un SNP bialélico seleccionado sobre A, C, G o T, nunca `N`, identificado por referencia estructural, secuencia, coordenada e ID, y presente una vez por copia aplicable. Para ploidías 1–10, cada dosis alternativa posible tiene valor genético explícito. El locus queda protegido: un evento explícito solapante se rechaza y uno aleatorio nunca lo elige; la recombinación puede transmitirlo sin perder identidad. Antes de seleccionar, cada elegible tiene dosis, carácter y aptitud inequívocos; las aptitudes son finitas, no negativas y alguna es positiva. Máximo cero, aptitud negativa/no finita o locus ausente/duplicado detienen antes de supervivencia y no son extinción. Se rechazan otras clases de locus y funciones no direccionales | SNP válido por varias ploidías; `N` rechazada; solapamiento explícito/aleatorio; máximo cero; aptitud inválida; locus ausente/duplicado; recombinación conservadora |
| **ACE-83. Extinción y migración** | La extinción planificada invalida dependencias posteriores detectables. Con N realizado=0, la política global desactivada puede sustituirse por fase y subpoblación: desactivada detiene toda la simulación; activada registra extinción por resultado y omite atómica y transitivamente eventos dependientes con cadenas de causas, mientras continúan ramas independientes. Recolonizar crea otro ID. La migración simultánea usa pares proporcionales o exactos, sin reemplazo ni A→B→C, cambia conjuntos reproductores pero no N objetivo. Antes de simular se rechaza cualquier combinación de rutas que pueda dejar individuos basados en referencias estructurales distintas dentro de una unidad reproductiva; las subpoblaciones basadas en referencias distintas solo pueden coexistir si permanecen reproductivamente aisladas | Reglas técnicas de asignación conjunta y casos con referencias iguales y distintas |
| **ACE-84. Deriva y orden generacional** | La deriva emerge del muestreo finito y no existe un paso artificial adicional. Las ocho etapas de REQ-REP-24 tienen un orden único y no configurable. Una etapa opcional desactivada se omite en su posición. Cualquier campo, lista o valor con otro orden se rechaza antes de simular e identifica lo solicitado y la regla v1. Varios eventos demográficos exigen un orden total explícito dentro de la primera etapa y cada uno recibe el estado del anterior. El historial conserva etapas aplicadas y omitidas y ambos órdenes | Todas las etapas; una opcional omitida; alternativa rechazada; varios eventos demográficos; registro efectivo |
| **ACE-85. Recombinación, bivalentes y estructura** | Toda la ejecución usa tasa uniforme en cM/Mb o mapa físico/genético completo. Sin interferencia, los entrecruzamientos por bivalente siguen Poisson según la longitud genética. Cada cromosoma recibe al menos uno. Un *scaffold* lo recibe si supera un umbral general, sustituible por secuencia, expresado exclusivamente en cM o bases; si no le corresponde recombinación, transmite intacta la copia elegida. Manifiesto y validación conservan regla, valor y procedencia. Las ploidías pares forman al azar bivalentes solo cuando todas las copias y parentales comparten la estructura canónica, sin preferencia de subgenoma, multivalentes, doble reducción ni compatibilidad parcial | Método exacto de posiciones; mapas; *scaffold* intacto; SV fijada idéntica y cada SV heterocigota |
| **ACE-86. Linaje detallado de individuos muestreados** | Solo los individuos incluidos en muestreos o materializaciones adicionales solicitados exportan puntos de recombinación y origen de cada segmento. Cada segmento alcanza al individuo y copia haplotípica fundadores y las mutaciones de novo forman capas separadas. La ejecución conserva internamente lo necesario aunque no exporte antepasados intermedios; los resúmenes no enumeran ese detalle | Evolución muestreada/no muestreada y reconstrucción a través de un antepasado no materializado |
| **ACE-87. Trazabilidad y límites estructurales** | Las coordenadas actuales se relacionan con las fundadoras cuando existe correspondencia, conservando orientación. Cada copia duplicada y cada inserción nueva tiene identificador estable; las deleciones conservan el intervalo fundador ausente. En una copia haplotípica se permiten eventos estructurales no solapados, incluso adyacentes, pero se rechaza solapamiento, anidamiento o modificación de una estructura previa. La única excepción entre fuente y destino de un mismo evento es la duplicación cuya fuente abarca la secuencia completa. Haplotipos distintos pueden solaparse. SNP e indels se permiten en segmentos modificados, salvo si solapan o atraviesan un punto de ruptura; la adyacencia se permite. Historias distintas pueden producir el mismo estado estructural canónico sin perder su trazabilidad | Formato de relaciones y equivalencia canónica; casos de historias distintas |
| **ACE-88. Modelos de población inicial** | Cada locus o bloque y subpoblación usa modo probabilístico predeterminado o controlado explícito sin mezclarlos dentro de la unidad. Admite Hardy-Weinberg en ploidías 1–10: binomial por dosis en loci bialélicos y multinomial por composición en multialélicos. También admite frecuencias explícitas de genotipo o dosis y genotipos multilocus faseados. En Hardy-Weinberg controlado se fijan recuentos alélicos, pero los genotipos se distribuyen aleatoriamente condicionados al modelo, sin exigir recuentos genotípicos exactos salvo configuración explícita. Frecuencias y recuentos suman, respectivamente, 1 y el total aplicable; los datos incompatibles se rechazan | Casos HWE haploide, diploide y poliploide, bialélico y multialélico, por modo |
| **ACE-89. LD y fase iniciales** | Por defecto, la población inicial asigna loci independientemente. El LD usa bloques no solapados con haplotipos y frecuencias conjuntas probabilísticas por defecto o recuentos controlados por subpoblación que suman C. Las frecuencias haplotípicas son la fuente de verdad y de ellas se derivan las alélicas; una frecuencia alélica aportada solo valida y una discrepancia fuera de tolerancia es error. Los recuentos se distribuyen entre individuos condicionados al modelo o plantillas. Los bloques conservan fase e identificador a lo largo de las generaciones, sin generar ni calcular D, D' o r² | Formato técnico de bloques, tolerancias y casos de prueba |
| **ACE-90. Endogamia inicial** | F automático se configura para la población completa o se sustituye por subpoblación exclusivamente diploide, siempre entre 0 y 1 y nunca por locus, región o bloque. Usa las fórmulas bialélica o multialélica de REQ-EVO-47/48: F=0 coincide con Hardy-Weinberg, un valor intermedio modifica las probabilidades según esas fórmulas y F=1 produce homocigosis completa. F en ploidía 1 o mayor que 2 se rechaza antes de simular con ámbito y ploidía. En poliploides se usa HWE sin F o una distribución explícita local, que prevalece. F realizado se calcula solo en diploides, puede ser negativo sin constituir error y figura como «no aplicable» en las demás ploidías | Casos F=0, intermedio y 1; rechazo P=1 y P>2; distribución explícita poliploide; F realizado no aplicable |
| **ACE-91. Fundadores trazables** | Todos los individuos de la población inicial son fundadores con identificadores estables propios. Cada copia haplotípica también tiene identificador propio aunque su secuencia coincida con otra, y toda traza posterior alcanza una copia fundadora concreta | Formato de identificadores y casos de reconstrucción |
| **ACE-92. Selección de muestreos** | Un evento admite identificadores explícitos, selección aleatoria sin reemplazo o ambas, además de estratos no solapados por subpoblación, grupos y filtros de genotipo. Los explícitos se incorporan primero y se excluyen del sorteo. Duplicados, individuos incompatibles, estratos ambiguos o tamaños mayores que los elegibles detienen el evento sin reemplazo ni reducción silenciosa | Sintaxis de filtros y estratos y casos de prueba |
| **ACE-93. Evento explícito inicial** | F4 admite inicialmente un único evento explícito con identificador, selección, una población, FASTQ y subsemilla. Generar FASTQ exige tecnología y perfil; BAM, CRAM y mapeador no forman parte de F4 y se incorporarán desde F5.6 sobre esos FASTQ finales. Un evento que solo materializa verdad o VCF carece de esos requisitos técnicos. Varios eventos o poblaciones se aplazan en DTF-050 | Rechazo de cardinalidad y caso de un evento |
| **ACE-94. Organización, atomicidad y exportación** | La raíz conserva una sola copia verificable de cada referencia e índice y un manifiesto general. Cada evento tiene directorio y manifiesto propios; la raíz enlaza su configuración resuelta, referencias y subsemilla incluso antes de que existan resultados. Cuando el evento produce resultados reales, ambos manifiestos los enlazan y los separan por referencia, sin crear resultados vacíos para salidas no activadas. Ejecución, evento, conjunto y artefacto distinguen estados incompleto, completo y fallido. Cada artefacto se publica atómicamente tras validarse; el éxito exige todas las unidades y la trazabilidad obligatoria completas. No se publica ni exporta una versión degradada por recursos. Los conjuntos completos independientes se conservan y el paquete autocontenido valida sus sumas antes y después | Preflight insuficiente, ausencia de degradación, esquema, sumas y formato |
| **ACE-95. Coordenadas y límites estructurales** | La configuración usa coordenadas basadas en 1 e intervalos cerrados; cada formato externo usa su convención oficial mediante conversión documentada. Cada duplicación usa el estado inmediatamente anterior para interpretar origen y destino y conserva las correspondencias posteriores. Una fuente completa puede insertar la copia dentro de la misma secuencia; una fuente parcial solo puede usar un destino adyacente o externo a su intervalo. Inserciones, deleciones, inversiones y puntos de ruptura respetan las reglas terminales, y toda translocación conserva subeventos, segmentos y orientaciones sin vaciar una secuencia | Representación técnica y casos límite |
| **ACE-96. Límite local de copias** | El máximo local predeterminado y provisional es 10 copias totales, incluida la original, por familia de copia e intervalo fundador dentro de cada copia haplotípica. Se contabiliza de forma acumulada y separada de ploidía y repeticiones en tándem. Un valor mayor requiere configuración explícita, genera una advertencia y activa una estimación previa de memoria, almacenamiento y volumen de datos. Si supera el presupuesto configurado se rechaza; si no, continúa. No existe otro máximo absoluto fijo y el valor predeterminado debe revisarse mediante pruebas de rendimiento | Pruebas de rendimiento y casos acumulativos y de presupuesto |
| **ACE-97. Reproducción de familias** | Cada evento tiene un parental para autofecundación o dos para cruzamiento sexual; otras cardinalidades se rechazan. Todos los parentales usan la misma referencia estructural, ploidía par e igual y estado estructural canónico idéntico en cada secuencia reproductiva. El orden de dos parentales no cambia el modelo y los diseños con tres o más fundadores forman un grafo. Clonación, apomixis, dobles haploides, referencias distintas y parentales estructuralmente heterocigotos se rechazan | Casos por tipo, diploide/tetraploide, ploidía impar/desigual, estructura y referencia |
| **ACE-98. Segregación y mutación familiar** | La transmisión elige con igual probabilidad una copia de cada bivalente estructuralmente idéntico, sin conversión génica, distorsión ni impulso meiótico. SNP, *indels* pequeños y repeticiones distintas se transmiten como alelos. Deleción, duplicación, inversión, inserción estructural o translocación heterocigotas impiden reproducir; el mismo estado estructural fijado en todas las copias y parentales sí se admite aunque tenga historias distintas. Una SV de novo puede existir en la descendencia, pero el portador queda fuera de intentos posteriores mientras sea incompatible; su exclusión y los fallos resultantes se registran conforme a N | Casos por clase SV, SV fijada, SV de novo con reproducción posterior y variantes pequeñas |
| **ACE-99. Independencia VCF entre referencias** | Cada referencia produce un VCF construido directamente en su secuencia y coordenadas, con un espacio propio de identificadores. Muestras biológicas pueden compartir identificador, pero VCF, complemento y manifiestos no contienen identificadores globales ni correspondencias registro a registro entre variantes de referencias diferentes | Casos con dos referencias y revisión de metadatos |
| **ACE-100. Identidad técnica y biológica** | `SM` de BAM/CRAM y las columnas del VCF usan el identificador biológico estable. Cada individuo tiene un grupo de lectura por evento, único en toda la ejecución, que distingue evento, biblioteca y tanda y coincide entre FASTQ, BAM/CRAM cuando existan y manifiesto. Cada instancia de lectura o pareja tiene identidad única dentro del evento. En F4, los vínculos de copia biológica, transferencia técnica o duplicación se conservan solo transitoriamente para generación y validación; el manifiesto no publica identificadores individuales de procedencia | Cabeceras, grupos de lectura, validación por reproducción y casos de colisión |
| **ACE-101. VCF estable entre eventos** | Dos eventos con los mismos individuos, referencia y verdad producen VCF con idénticos identificadores locales, coordenadas, REF, ALT y genotipos aunque cambien tecnología, subsemilla, FASTQ, BAM/CRAM o degradaciones | Casos con eventos técnicos diferentes |
| **ACE-102. ALT fijados y multialélicos** | Un ALT fijado en todas las muestras se incluye porque difiere de la referencia. Varios ALT representables en el mismo locus se reúnen en un registro multialélico normalizado y `GT` conserva ploidía, dosis y fase. Los casos no representables usan el complemento sin pérdida | Casos fijados, multialélicos, poliploides y solapados |
| **ACE-103. Variantes pequeñas solapadas** | SNP e *indels* pueden solaparse en una copia haplotípica. Los explícitos exigen un orden, se aplican secuencialmente y usan el estado inmediatamente anterior como marco de coordenadas. La verdad actual y el VCF derivan de la secuencia final; el historial marca cada evento como activo, modificado o eliminado. Un objetivo ya eliminado detiene la configuración; no se ignora ni recoloca. Un SNP eliminado por una deleción no se emite como variante actual. Haplotipos distintos se resuelven independientemente y sus alelos finales pueden coexistir como ALT con fase y dosis. Varios eventos de un haplotipo que formen un alelo final se emiten como una sola variante o ALT compleja normalizada y usan complemento si hace falta | Casos secuenciales, multialélicos, complejos, de eliminación y de error |
| **ACE-104. Modos de frecuencia alélica** | Cada locus o bloque usa exclusivamente modo probabilístico o controlado. El probabilístico es predeterminado, muestrea desde las probabilidades y registra esperado y realizado sin tolerancia controlada. El controlado, solo explícito, elige un vector de recuentos que suma C y minimiza conjuntamente la desviación; exige 0,5/C o una tolerancia más estricta. Después distribuye copias aleatoriamente entre individuos condicionado a Hardy-Weinberg binomial o multinomial, F automático diploide o la distribución explícita. En bloques LD controla recuentos haplotípicos por subpoblación. Restricciones incompatibles detienen y los empates se resuelven con la semilla y se registran | Casos por ploidía, multialelismo, bloques, modelo y restricciones |
| **ACE-105. Normalización de probabilidades y frecuencias** | Todo conjunto que deba sumar 1 se compara con una tolerancia general configurable. Una desviación dentro de ella se normaliza y registra con valores y suma antes y después; una desviación mayor o una suma no normalizable produce un error. El valor predeterminado se fija en implementación. Esta tolerancia no sustituye `0,5/C` ni otras tolerancias específicas con otra finalidad | Valor predeterminado, precisión numérica y casos dentro, fuera y no normalizables |
| **ACE-106. Catálogo declarativo común** | División, fundación, migración y selección de descendientes usan el mismo catálogo: identificadores, azar, población o grupo, alelo, genotipo, dosis o haplotipo en uno o varios loci, Y/O/NO y cantidad exacta o proporcional. Las proporciones son probabilidades individuales reproducibles y los números exactos seleccionan sin reemplazo o fallan por insuficiencia. Las reglas usan la verdad biológica. Scripts, código, funciones personalizadas y reglas desconocidas se rechazan; una regla nueva exige ampliar el programa y la especificación | Sintaxis técnica de las reglas y casos compartidos por los cuatro procesos |
| **ACE-107. Colisiones en la generación aleatoria** | Los eventos explícitos se aplican primero en su orden declarado. Después, cada clase cuenta su espacio inicial de oportunidades, sortea una vez su recuento Poisson y se procesa en el orden duplicaciones, deleciones, inversiones, translocaciones, repeticiones en tándem, *indels* y SNP. Cada evento se elige solo entre candidatos válidos en el estado anterior y actualiza el espacio para el siguiente, sin reintentar candidatos inválidos, volver a sortear ni aplicar un límite de intentos. Si faltan destinos para completar el recuento, la generación falla sin aceptar menos ni advertir; se registran oportunidades, media y números sorteado, aceptado antes del fallo y finalmente activo | Casos con colisiones, agotamiento progresivo, denominador inicial y diferencia entre los tres recuentos |
| **ACE-108. Cardinalidad entre verdad, VCF y complemento** | Cada conjunto solicitado contiene exactamente un VCF multimuesta estándar. Un registro F2 que agrupa ALT de varios casos usa `GPS_CASE` multivaluado con IDs originales únicos y canónicos; `truth.jsonl` conserva ALT↔casos y no se crea un caso compuesto. Un caso plenamente representable no crea una entrada complementaria. Un caso complejo puede usar un registro VCF y una entrada complementaria. Una nueva adyacencia usa normalmente dos BND relacionados y, si pierde detalle, ambos enlazan una entrada; una translocación recíproca usa los BND necesarios para sus dos adyacencias, que pueden ser cuatro, y una sola entrada si pierde detalle. Todos los campos `INFO` adicionales están declarados en la cabecera. Dos referencias producen VCF y espacios de identificadores separados | Dos ALT/dos casos; un ALT/dos casos; un caso/varios registros; orden/duplicado/ID inexistente; ausencia de complemento F2; caso complejo de un registro, par BND, translocación recíproca y dos referencias; validación con la versión oficial seleccionada |
| **ACE-109. Poisson y modelo escalonado simple** | Todo recuento de variantes aleatorias usa exclusivamente Poisson y se reproduce con la misma semilla, configuración, entradas y versión. Los cambios concretos se expresan mediante eventos explícitos o el modo controlado aplicable, no como otra distribución. Cada mutación de repetición en tándem añade una unidad completa o elimina una, conforme al sesgo configurado. Seleccionar otra distribución de recuentos, otro modelo de repetición o un salto de varias unidades se rechaza antes de simular e identifica parámetro y valor | Casos Poisson reproducibles, eventos explícitos, expansión y contracción unitarias, valores de sesgo y rechazos de alternativas |
| **ACE-110. Conservación de la base desconocida N** | `N` se conserva y nunca se resuelve. No puede ser SNP, locus seleccionado ni motivo de repetición; sí puede permanecer dentro de una SV o secuencia explícita/copiada. Las inserciones aleatorias y errores solo generan A/C/G/T. Se hereda y recombina sin cambios. WGS puede muestrearla, FASTQ la representa con Phred 0, BAM/CRAM la propaga, GC la excluye de su denominador y cobertura la cuenta. VCF nunca emite SNP con `REF=N` y usa complemento cuando el estándar no conserve otro caso sin pérdida. Validación y manifiesto prueban exclusiones, resúmenes y ausencia de resolución | SNP sobre `N`; SV con `N`; inserciones aleatoria y expresa; herencia; lectura Q0; ventana toda `N`; cobertura; VCF/complemento; ausencia de resolución |
| **ACE-111. Sustituciones y precedencia entre ejes** | Cada parámetro parte de su regla base y solo admite los selectores de su familia. Una regla aplica cuando coinciden todos sus selectores; se eliminan las dominadas y una única máxima fija el valor. Una regla combinada domina las reglas parciales. Sin esa regla combinada, región frente a subpoblación e individuo frente a región son máximas incomparables y se rechazan aunque compartan valor; también se rechazan reglas duplicadas, regiones solapadas con iguales selectores no espaciales y selectores no admitidos. Generación domina fase, subpoblación domina población, locus domina región y secuencia, y una coincidencia concreta domina la ausencia dentro de los demás ejes. Un evento técnico nunca sustituye biología. Configuración efectiva y manifiesto conservan valor, regla ganadora, selectores y procedencia | Herencia base; todas las jerarquías internas; regla combinada; los dos conflictos cruzados; duplicado igual; regiones solapadas; selector no permitido; separación biológica/técnica; procedencia |
| **ACE-112. Objetivo y selección de referencia reducida** | Un FASTA real validado origina otra referencia con ID propio. Exactamente uno de `target_bases` o `reduction_percent` fija un objetivo entero positivo y menor que el origen; `region_size_bases` es positivo. El porcentaje expresa bases eliminadas y retiene `floor(S × (100 − p) / 100)` con decimal exacto. El selector directo proporcional usa subsemilla propia, reparte bloques completos por longitud con capacidad comprobada y compone huecos sin enumerar bases; no reemplaza ni solapa, no duplica ni invierte y produce bloques completos más como máximo un residual para alcanzar exactamente el objetivo. Objetivo imposible, espacio insuficiente o parámetros incompatibles fallan antes de materializar | Casos con objetivo exacto/múltiplo/residual, porcentajes próximos a 0 y 100, rechazo de 0/100, objetivo igual/mayor que origen, semilla repetida, cromosomas/*scaffolds*, vector canónico, no solape y espacio insuficiente |
| **ACE-113. Procedencia y reconstrucción de referencia reducida** | Los intervalos se ordenan por secuencia y coordenada, se concatenan en orientación directa, omiten secuencias sin selección y no crean secuencias vacías. Cada cromosoma o *scaffold* conserva `sequence_id` y clase. `N` se copia y cuenta. FASTA canónico, mapa derivada↔origen, SHA-256, configuración efectiva y manifiesto permiten reconstruir exactamente cada base y explican secuencias omitidas y residual | Reconstrucción con lector independiente, mapa con bordes base 1 cerrados y vectores base 0 semiabiertos, `sha256sum`, FASTA mixto con `N` y corrupción deliberada de cada enlace |
| **ACE-114. Equivalencia de autoría Marimo** | Todas las opciones implementadas están en una de las nueve pestañas mínimas. Escalares usan controles tipados; listas y reglas, tablas con alta/baja/reordenación y selectores visibles. Importar, editar, previsualizar y exportar como unidad el TOML general y cero o varios escenarios conserva significado, IDs, orden y enlaces/rutas relativas. El binario Rust valida y resuelve; cualquier cambio de cualquier documento invalida la validación previa y una incompatibilidad, escenario ausente o duplicado se rechaza | Matriz opción→pestaña/control; conjunto sin escenarios y conjuntos multiescenario con rutas relativas/espacios/faltantes/duplicados; comparación UI→conjunto TOML→`EffectiveConfig` o diagnóstico con conjunto manual equivalente |
| **ACE-115. Ejecución segura y accesible desde Marimo** | La aplicación nunca ejecuta automáticamente ni usa shell. Una acción explícita tras validación vigente lanza `run` con argumentos separados y destino no-clobber revisado. Muestra estado, código, logs, diagnósticos y enlaces/resúmenes sin cargar artefactos genómicos grandes. Funciona localmente en Ubuntu/WSL2 Linux y ofrece teclado, foco, etiquetas, ayuda y estados no basados solo en color. Cada capacidad F2–F9 añade controles antes de cerrar | Proceso falso y binario real; cambio posterior a validar; destino existente; logs grandes y artefacto grande; revisión de teclado/foco/color; puerta transversal por hito |
| **ACE-116. Población inicial desde frecuencias SNV** | La ejecución prepara, genera o reutiliza el par canónico TSV+JSON enlazado a una referencia exacta. Solo admite `alt_af`; genera loci bialélicos candidatos sin reemplazo mediante espacio, sustituciones y espectro explícitos y acepta multialelia suministrada. Las poblaciones asignan frecuencias independientemente, permiten 0/1 solo mediante `fixed`, usan siempre intervalos continuos `[min,max)` que no producen esos extremos, eliminan `ALT` globalmente ausentes sin reemplazo y conservan fijados. En modo `exact`, mayores restos fija las cuotas y una permutación reproducible independiente por población asocia las filas a todos los `variant_id`; el valor usa otro flujo, sin consumo para `fixed` ni cuotas cero, y nunca se asigna por coordenada. HWE es predeterminado y `F` solo diploide; cada individuo pertenece a una población y recibe genotipo/fase reproducibles sin LD ni corrección por muestra finita. Los dos perfiles humanos incorporados son aproximaciones descriptivas explícitas, nunca predeterminadas, exactas por superpoblación ni tasas universales; los perfiles suministrados pueden diferir por población. Las dos excepciones de ejecución son específicas, visibles y trazadas | Preparación y hashes; tabla y VCF rechazado; modos exacto con cuotas, permutación y flujo de valor separados, y estocástico; fronteras 0/1; ausencia de correlación espacial; aislamiento poblacional; consumo cero `fixed`/cuota cero; regiones/pesos; 12 sustituciones y Ti/Tv; perfiles uniforme/SFS 1000G aproximado/SNV12 Harris aproximado; receta compacta, fuente CC BY y SHA-256; ceros/fijados; matriz vacía; multialélico; HWE/F; subsemillas y reordenación; flags limitados; errores de suma, duplicado, referencia y `N` |
| **ACE-117. Individuos desde panel de haplotipos** | Un paquete canónico de cinco ficheros contiene una única referencia, bloques locales, SNP A/C/G/T, HAP únicos y frecuencias densas y lo prepara una única orden Rust. Cada bloque tiene al menos dos sitios y HAP; admite población monomórfica y ALT fijado, pero rechaza sitio siempre REF, referencia/secuencia incompatible y loci `alt_af` dentro del intervalo. `haplotype_construction` se resuelve por población×bloque y se separa de `frequency_sampling`; intacto admite probabilístico/controlado y HWE, F diploide o distribución explícita; recombinante exige solo probabilístico, r en (0,1] y HWE o F diploide, y rechaza una distribución explícita. El productor recombinante opera directamente por individuo×bloque con fuentes exclusivas del panel y el constructor común uniforme cerrado por DTC-055; DTPD-111 solo conecta mapas y meiosis en F7.1. En diploides Bernoulli(F) ejecuta uno duplicado o dos independientes; P≠2 ejecuta P sin F. Un punto x es x\|x+1, alterna segmentos *forward* desde la primaria y conserva cruces invisibles. Un vector conocido reutiliza HAP y uno nuevo deriva ID por DTC-039 desde panel+bloque+contenido. RNG DTC-037 separa F, primaria, puntos y secundaria. Todo productor tiene `process_id`; `origin_kind` distingue exactamente K=0 intacto y K≥1 recombinado. C=N×P copias finales se agregan como `n_h/C` y separan recuentos intactos y recombinados. Tres TSV con cabeceras, tipos, orden y claves exactos conservan asignaciones, segmentos y frecuencias; clase/versión están en sus enlaces del manifiesto. Marimo invoca Rust, presenta controles compatibles y exige validación vigente. | Oráculo independiente que reconstruye HAP/mosaicos/IDs/asignaciones y verifica referencia única, C, `n_h/C`, recuentos por origen y contratos TSV; r=0 rechazado, cruces conocidos/múltiples/invisibles, F=0/intermedio/1, combinaciones de modelo admitidas y rechazadas, fuentes solo del panel, repetición/aislamiento, corrupción de tablas y enlaces, tres salidas y paridad Marimo↔TOML |

Estos criterios son cualitativos y deberán convertirse en pruebas durante la
planificación de la implementación.

## 8. Estado de las decisiones funcionales

Las 27 preguntas funcionales de la revisión
[revision-adversarial-especificaciones-0.1.md](./revision-adversarial-especificaciones-0.1.md)
están resueltas para v1. La ampliación posterior confirma la referencia reducida
y la aplicación Marimo; el 14 de agosto de 2026 quedaron también cerrados el
redondeo de `reduction_percent` y la elegibilidad de *scaffolds*. Añadir otro comportamiento,
selector o clase de regla requerirá revisar expresamente esta especificación;
las elecciones restantes de la sección 9 son técnicas y no pueden reabrir la
semántica funcional confirmada.

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
utilizado y se elegirán las herramientas de indexación. En v1, el único índice
de referencia admitido es un FAI creado bajo demanda. La coexistencia de varios
tipos de índice se decidirá posteriormente mediante DTF-048.

**Generación aleatoria y reproducibilidad.** Durante la implementación y la
validación se aplicará `poisson-count-v1` de DTC-067, que resuelve DTPD-120 con
ramas Knuth/PTRS, `lambda<=2^52`, consumo RNG y matemática `libm` exactos,
vectores congelados y oráculo estadístico independiente. No se evaluará otra
distribución como opción de v1 ni un consumidor podrá recortar o volver a
sortear el resultado. También se aplicará la estrategia técnica de subsemillas o
  flujos aleatorios independientes y cómo demostrar las garantías separadas de
  REQ-CONF-111 a REQ-CONF-115 entre versiones del programa, dependencias y
  plataformas. Estas elecciones conservarán las medias esperadas configurables,
  los resultados realizados y los niveles de reproducibilidad ya exigidos.

**DECISIÓN TÉCNICA RESUELTA — coordenadas y selección de intervalos para
referencias reducidas.** DTPD-101 y DTPD-116 quedan resueltas mediante DTC-056 y
DTC-057. El dominio usa coordenadas internas base 0 semiabiertas y tipos
etiquetados, mientras configuración y mapa conservan base 1 cerrada. El
selector directo proporcional canónico reparte bloques completos por longitud,
compone huecos sin recorrido por bases y reserva después el único residual;
materializa en orden de secuencia y coordenada, con vector congelado y
reconstrucción independiente.

**DECISIÓN TÉCNICA RESUELTA — estado haplotípico y fichero obligatorio de
verdad.** DTPD-102 queda resuelta por DTC-061. El estado interno es híbrido y
segmentado, comparte referencias/fuentes y solo materializa secuencias completas
bajo demanda. Su salida obligatoria es `truth.jsonl` v1, UTF-8 sin compresión,
con cabecera inicial, registros tipados de copias/genotipos, orden canónico,
recuentos/digests y longitud/SHA-256 en manifiesto. Esquemas por registro y un
reconstructor semántico independiente deben detectar relaciones
ausentes/duplicadas/desordenadas y exigir igualdad exacta entre genotipos y
haplotipos. Todo VCF deriva de `truth.jsonl`; la compresión de la verdad queda
para una fase posterior sin cambiar este modelo lógico.

**DECISIÓN TÉCNICA RESUELTA — contrato VCF.** DTC-062 fija
VCF 4.5 como única versión de salida, exige la línea exacta
`##fileformat=VCFv4.5`, prioriza las representaciones estándar ampliamente
interoperables y exige validación externa basada en htslib. DTC-063 fija como
salida predeterminada `.vcf.gz` BGZF con exactamente un índice: TBI solo si
todas las secuencias y posiciones caben en su límite normativo, aproximadamente
2^29 bases por secuencia, y CSI si alguna lo supera. La alternativa explícita
es `.vcf` de texto sin compresión ni índice. `EffectiveConfig` y manifiesto
registran modo, tipo y rutas; VCF e índice se validan y publican atómicamente.
DTC-064 fija además una normalización interna determinista y versionada desde
los haplotipos finales: reconstrucción previa de cambios solapados sobre un
intervalo `REF` común, recorte con anclaje, alineamiento a la izquierda y orden
por referencia/coordenada. F2 emite como máximo un registro por `CHROM+POS`,
agrupa todos los alelos compatibles como ALT multialélicos y no los descompone;
`bcftools norm` es solo un oráculo externo. Varios BND relacionados en una
posición quedan como excepción futura para SV/translocaciones. DTC-065 exige
`GPS_CASE` multivaluado, único y canónicamente ordenado en todo registro, y
`truth.jsonl` conserva la relación de cada ALT con sus casos fuente. Solo ante
pérdida aparecen `GPS_LOSS` con código controlado y `GPS_COMPLEMENT`; el único
`vcf-complement.jsonl` versionado del
conjunto contiene una entrada completa por caso perdido y está ausente si no
hay pérdidas. DTC-066 cierra DTPD-103 con aceptación/CI externa mediante
`bcftools` y HTSlib/tabix 1.24, nunca en runtime, más un oráculo semántico
independiente y un corpus positivo/de mutaciones. EBI `vcf-validator` no se usa
como autoridad VCF 4.5.

**DECISIÓN TÉCNICA RESUELTA — recuentos Poisson.** DTC-067 resuelve DTPD-120
con `poisson-count-v1`: `lambda` finita en `0..=2^52`, cero sin consumo RNG,
Knuth para valores menores que 10 y PTRS de Hörmann desde 10, uniformes abiertos
derivados de un `next_u64` ChaCha12 y matemática de `libm` 0.2.16. Un valor o
candidato fuera de los límites se rechaza sin consumir/recortar/resortear; un
recuento superior al espacio falla mediante la regla existente. Configuración
efectiva/manifiesto conservan versión, media, ámbito/subsemilla y sorteo, y los
vectores/oráculos cubren fronteras y estabilidad Linux/WSL2.

**DECISIÓN TÉCNICA RESUELTA — perfiles humanos compactos.** DTC-068 resuelve
DTPD-121 con `human_1000g_global_sfs_approx_v1` y
`human_harris2017_snv12_approx_v1`, ambos explícitos, compactos, aproximados y
nunca predeterminados. El primero normaliza pesos 64:12:8 sobre
`[1/5008,0.005)`, `[0.005,0.05)` y `[0.05,1)`; el límite 1 completa la cobertura
pero nunca se sortea, y el mínimo `1/5008` pertenece a
la receta publicada y no es universal. El segundo transforma las cinco filas
continentales y 96 contextos de Harris y Pritchard a 12 sustituciones mediante
promedio equitativo, colapso por hebra y complemento inverso. Fuente, licencia
CC BY, fichero de 3.644 bytes, SHA-256, transformación y digests de salida se
validan en F2.2.6 con un oráculo independiente. F2 no procesa el VCF Phase 3 ni
usa `bcftools` para estos perfiles; el perfil empírico exacto queda en DTF-016.

**DECISIÓN TÉCNICA RESUELTA — fronteras y selección aleatoria F2.** DTC-075
conserva los intervalos semiabiertos, los extremos `fixed`, la selección SNV
estática y la selección dinámica de indels. DTC-079 permite agrupar las
utilidades y sus pruebas: no se exige un productor, expediente o paso separado
para cada primitiva interna. Los contratos observables y la reproducibilidad
siguen siendo obligatorios.

**DECISIÓN TÉCNICA RESUELTA — aplicación Marimo simplificada.** DTC-079
sustituye para F2 la concreción sobrediseñada de DTC-058. Marimo usa el entorno
fijado por `pyproject.toml`/`uv.lock`, escucha en `127.0.0.1`, edita un
único TOML general, invoca directamente `validate`/`run` y muestra
diagnósticos y resúmenes. Rust conserva la semántica. F2 prueba instalación y
arranque local; F9 conserva distribución, plataforma oficial, accesibilidad de
entrega y endurecimiento.

**DECISIÓN TÉCNICA APLAZADA — importación/exportación multidocumento.** DTC-069
se recuperará con la primera fase que implemente escenarios reales. F2 importa y
exporta únicamente el TOML general de la vía directa, conserva su significado y
no sobrescribe accidentalmente un destino. No crea un grafo de documentos,
copias privadas, preservación opaca ni defensas TOCTOU sin consumidor.
**DECISIÓN TÉCNICA RESUELTA — matrices y espectros SNV.** DTPD-118 fija el par
canónico TSV+JSON, la preparación y generación reutilizables, `alt_af` como única
frecuencia, intervalos no solapados con muestreo fijo/lineal/logarítmico, modos
estocástico y exacto por mayores restos, poblaciones independientes, espacio
uniforme o ponderado y 12 sustituciones dirigidas. Los perfiles humanos
aproximados DTC-068 son ejemplos versionados de selección expresa, no valores
predeterminados ni perfiles exactos por población. La
validación debe aportar vectores y productores/oráculos independientes. DTC-060
liga la preparación a una entrada identificada y verificada del manifiesto F1;
una ruta FASTA por sí sola no puede aportar `reference_id`.

**DECISIÓN TÉCNICA RESUELTA — receta compacta.** DTPD-121 queda cerrada por
DTC-068: F2.2.6 regenera recursos pequeños desde los literales publicados y el
fichero eLife fijado, verifica atribución, tamaño, SHA-256, transformaciones y
digests y no descarga ni procesa el VCF 1000 Genomes. La materialización exacta
del VCF y perfiles empíricos exactos quedan en DTF-016/post-v1 y no bloquean F2.

**DECISIÓN TÉCNICA RESUELTA — longitudes de indels F2.** DTC-059 fija por
separado las longitudes exactas de inserción y deleción aleatorias, ambas con
1 pb por omisión y dominio desde 1 pb hasta justo antes del umbral estructural
efectivo. F2 no sortea longitudes ni admite una distribución para ellas.

**BASE FUNCIONAL APROBADA — panel de haplotipos.** El panel representa
bloques de LD locales, no cromosomas completos ni una identidad de haplotipo que
se prolongue entre bloques. Cada bloque posee intervalo explícito y no solapado;
los SNP exteriores a esos intervalos proceden de la matriz `alt_af` y ninguno de
sus loci puede caer dentro de ellos. Ya se ha fijado que v1 solo admite SNP
suministrados A/C/G/T, que
`panel_haplotype_alleles.tsv` es densa y no admite asignaciones implícitas ni
alelos desconocidos, y que `panel_frequencies.tsv` es densa por población,
bloque y haplotipo. También queda fijado el paquete exacto de cuatro TSV y un
`panel.metadata.json` que los enlaza por ruta, tamaño y SHA-256 sin duplicar sus
tablas, y la orden Rust única que lo prepara sin modificar las entradas. El
runtime solo acepta el paquete íntegro y compatible. La conversión desde VCF
queda fuera de v1 en DTF-041: en esta versión se suministran las cuatro tablas o
se generan externamente. La alternativa recombinante también pertenece a F6.4:
genera directamente cada individuo y bloque mediante productores que usan solo
el panel original; Bernoulli(F) decide uno duplicado o dos independientes en
diploides. La probabilidad uniforme, los puntos x|x+1, la identidad por contenido,
los ámbitos aleatorios, las tres salidas obligatorias y la interfaz Marimo están
fijados.

**DECISIÓN TÉCNICA RESUELTA — paneles de bloques haplotípicos.** DTPD-119 queda
resuelta mediante DTC-055. Los TSV del panel, la semántica de origen/proceso, los
modelos genotípicos y los contratos intacto/recombinante están fijados.
`intact` admite HWE, F diploide o distribución explícita; `recombinant`, solo HWE
o F diploide. Condicionar distribuciones explícitas sobre resultados
recombinantes queda aplazado.
Incorporar *indels*, repeticiones, variantes estructurales,
alelos simbólicos o *breakends* al panel queda en DTF-040 y no presupone un
formato futuro.

**Tolerancias numéricas.** Durante la implementación se fijarán el valor
predeterminado y la precisión de la tolerancia general usada cuando frecuencias
o probabilidades deban sumar 1. La configuración podrá sustituir ese valor. Las
pruebas deberán cubrir normalización dentro de tolerancia, rechazo fuera de ella
y sumas no normalizables, conservando los valores antes y después. También se
fijarán las tolerancias técnicas específicas aún pendientes. Ninguna sustituirá
`0,5/C` ni otra tolerancia ya acordada con una finalidad distinta.

**Representación técnica de coordenadas.** La configuración y el mapa de
reducción usan coordenadas basadas en 1 e intervalos cerrados. DTC-056 fija la
representación interna base 0 semiabierta, los puntos interbase separados y un
único servicio con conversiones comprobadas que reutilizan F1.6, F2 y las fases
posteriores. Cada formato externo conserva su convenio oficial. Los VCF de
referencias distintas se construirán de forma independiente y no se definirá un
algoritmo de correspondencia entre sus registros.

**Ejecución de los procesos evolutivos.** Se elegirán los algoritmos que apliquen
la secuencia obligatoria confirmada, N objetivo, intentos y realizado, y Ne realizado,
selección, migración, eventos de
subpoblaciones, reproducción y recombinación. Se concretarán la fórmula de Ne
reproductivo y las reglas para declarar cuándo no es evaluable. La fórmula solo podrá usar
los progenitores con descendencia y la desigualdad o varianza de sus
contribuciones. También se concretará la fórmula de F realizado exclusivamente
para diploides. Esta deberá comparar la heterocigosidad observada y esperada
según las frecuencias alélicas, y conservar resultados negativos; en ploidía 1
o mayor que 2 el resultado seguirá siendo «no aplicable». Estas decisiones
técnicas no deben cambiar los significados funcionales ni las consecuencias de
validación acordadas.

**Estado genotípico faseado.** DTC-061 fija el estado interno híbrido de cadenas
ordenadas de segmentos/modificaciones y el fichero externo obligatorio
`truth.jsonl` v1, con esquemas por tipo y validación semántica independiente.
Fase, dosis, copia, procedencia y coordenadas deben mantenerse hasta las copias
haplotípicas finales y permitir reconstrucción completa. Las fases posteriores
extenderán los tipos de modificación para mutaciones de novo, duplicaciones y
otras SV sin sustituir este modelo ni exigir materializar todos los genomas
históricos.

**Representación VCF y complemento de verdad.** DTC-062 a DTC-066 resuelven
versión, cabecera, compresión/índice, normalización, campos de enlace,
`vcf-complement.jsonl` y validadores de aceptación de DTPD-103. Las fases de
repeticiones en tándem y SV concretarán únicamente sus campos estándar o
específicos y el criterio técnico de tamaño manejable, sin reabrir el contrato
común. No podrán alterarse la cardinalidad entre casos, registros y entradas ni
la semántica oficial de campos estándar. La validación debe demostrar la
validez independiente del VCF y la reconstrucción sin pérdida a partir de las
salidas enlazadas. El esquema identifica inequívocamente referencia, conjunto y
muestras, usa `GT` con fase y ploidía reales y excluye `DP`, `AD`, `GQ`, `PL` y
cualquier otro campo `FORMAT` derivado de observaciones o inferencias.

**Repeticiones en tándem.** Se elegirá la representación técnica de la tasa y del
sesgo entre expansión y contracción del modelo escalonado simple. Cada mutación
seguirá cambiando exactamente una unidad completa. También se concretarán los
modelos opcionales de mutaciones internas ya admitidos. Los límites prácticos de
longitud y número de loci dependerán de los recursos, no de un límite funcional
de motivo.

**Estimación de recursos.** Se elegirán y calibrarán los modelos que estimen
almacenamiento, memoria y tiempo y los umbrales que generan advertencias. Las
pruebas compararán las estimaciones con consumos reales sin convertirlas en
garantías exactas.

**Perfiles WGS, mapeo y cobertura.** Se concretarán y validarán los perfiles WGS,
parámetros y versiones de Illumina-like *single-end* y *paired-end*. Se
elegirán las distribuciones, tasas y valores cuantitativos exactos de
longitud, error y calidad y se comprobarán con datos adecuados para cada perfil.
También se fijarán y probarán los *presets* y valores que minimap2 aplica
automáticamente y la comprobación del muestreo uniforme y la cobertura
aproximadamente Poisson. La planificación elegirá la versión concreta de
minimap2, el catálogo inicial de mapeadores alternativos y las comprobaciones
técnicas de compatibilidad de herramientas, perfiles, parámetros y referencias.
Se elegirán además las distribuciones, tasas y valores cuantitativos exactos de
duplicados, adaptadores, sesgo GC, quimeras, anomalías *paired-end* y las
categorías y MAPQ sintéticos del modo directo. Estos valores se validarán con
datos adecuados para Illumina-like. La elección cuantitativa no
elimina el requisito de que cada degradación sea simulable, configurable y
  trazable en v1. Tampoco cambia que minimap2 es el mapeador predeterminado,
  sustituible por configuración, que cada evento que genere alineamientos usa
  exactamente un mapeador ni que herramienta, versión y parámetros efectivos
  deben quedar trazados.

**Verdad y procedencia de lecturas.** F4 conserva únicamente el resumen agregado
definido por DTPD-108; no existe un esquema de verdad ni fichero por lectura.
La estimación
automática desde datos reales, el contrato de entrada, el esquema del perfil
empírico, sus métricas de fiabilidad y su reutilización no son decisiones de
planificación de v1, sino trabajo futuro.

**Validación interna contra la verdad.** Se elegirán las métricas, métodos
estadísticos, ámbitos de agregación y esquema de los resultados estructurados
que comparan FASTQ y BAM o CRAM con la verdad de origen. Las comprobaciones
deterministas aplicarán errores; las magnitudes aleatorias se informarán sin
tolerancias de aceptación y sin elevarse a error. Umbrales, intervalos
estadísticos, categorías de aceptación y escalamiento configurable quedan para
una versión futura. No se ejecutarán *SNP callers* para realizar esta validación.

**SAM/BAM, CRAM y grupos de lectura.** Se elegirán las versiones oficiales
exactas de SAM/BAM y CRAM, el catálogo concreto de campos de cabecera, el
procedimiento para marcar duplicados y la forma técnica de combinar grupos en el
único archivo final por individuo y referencia. Estas elecciones deben mantener
la identificación inequívoca de referencia, muestra, evento, biblioteca, tanda
y grupo de lectura; `SM` como identificador biológico estable; y el
identificador `RG` único en toda la ejecución. Todos los archivos deberán
cumplir la versión oficial elegida.

**Contaminación externa y dirigida.** F4 queda limitado al porcentaje total
sorteado por receptora y al reparto aleatorio entre individuos de su misma
población. La sintaxis de matrices, pesos dirigidos, FASTQ externos, sus
validaciones y los datos de integridad se aplazan en DTF-052 sin alterar la
contaminación sustitutiva ni la procedencia por identificadores y recuentos.

**Mutación.** Se elegirán las tasas predeterminadas, las distribuciones
cuantitativas de longitud de *indels* y el formato del espectro SNP. Estas
decisiones no podrán introducir dependencia de contexto SNP ni semántica
específica de elementos transponibles en la primera versión.

**Recombinación y reproducción poliploide.** Se elegirá el formato de los mapas,
la generación técnica de posiciones, la comprobación de equivalencia del estado
estructural canónico y la traza interna de bivalentes. La implementación deberá respetar el mínimo
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
   se elegirán la fórmula diploide de F realizado y la fórmula de Ne
   reproductivo, sin reabrir su significado funcional ni adaptar F a otras
   ploidías.
3. Definir la representación técnica de intervalos y las pruebas estadísticas
   del recuento Poisson de variantes aleatorias y del paso unitario de las
   repeticiones en tándem.
4. Elegir durante la planificación la sintaxis del grafo de cruces, de los
   parentales expresos y de las condiciones de selección ya confirmadas.
5. Concretar y validar con datos adecuados los perfiles Illumina-like
   *single-end* y *paired-end*, la versión, *presets* y parámetros
   de minimap2, el catálogo inicial de mapeadores compatibles y las tolerancias
   del muestreo uniforme y la cobertura
   aproximadamente Poisson.
6. Definir el algoritmo de sumas de comprobación, los esquemas de los
   manifiestos general y de evento y el formato de exportación autocontenida.
7. Implementar y comprobar el rango por receptora y reparto aleatorio de
   contaminación entre individuos de la misma población; las fuentes externas y
   dirigidas quedan en DTF-052.
8. Completar los perfiles cuantitativos de errores, calidades, duplicados,
   adaptadores, sesgo GC, quimeras, anomalías *paired-end* y demás degradaciones
   WGS admitidas.
9. Definir los campos de cabecera y las versiones de SAM/BAM, CRAM y VCF que
   identificarán referencia, muestra, evento, biblioteca, tanda y grupo
   de lectura.
10. Implementar y validar DTC-056/DTC-057 para coordenadas y el selector directo
    proporcional reproducible de referencias reducidas, incluido su vector
    canónico y reconstrucción independiente.
11. Validar mediante F2.0 la instalación bloqueada y el arranque local fijados
   por DTC-079, dejando empaquetado y plataforma oficial para F9,
    versionar su contrato con esquemas, catálogos y binario y no trasladar reglas
    biológicas a la interfaz.
12. Completar los casos de prueba restantes y tomar las decisiones de
    implementación cuando el comportamiento esperado esté suficientemente
    definido.
