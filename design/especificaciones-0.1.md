# Especificación del simulador de genomas y poblaciones

- **Versión:** 0.1
- **Estado:** borrador actualizado
- **Fecha:** 6 de agosto de 2026
- **Fuentes disponibles:** `README.md` e información confirmada por el usuario

## 1. Estado y finalidad del documento

Este documento define el objetivo general y los requisitos confirmados hasta
ahora. También identifica las decisiones que faltan antes de poder cerrar la
especificación y planificar la implementación.

La versión 0.1 sigue siendo un borrador. Los requisitos de esta versión son
obligatorios, pero todavía no se han acordado todos sus detalles ni todos los
casos que permitirán comprobar los resultados.

## 2. Objetivo del proyecto

El sistema debe generar un conjunto coherente de datos genómicos simulados a
partir de un fichero general de configuración. Cuando se use la vía evolutiva o
la vía de familias y poblaciones segregantes, el fichero general referenciará
además un fichero de configuración poblacional propio del escenario. El conjunto
generado estará formado por una o más referencias, cada una guardada en un único
fichero FASTA, individuos organizados en poblaciones o familias, archivos BAM o
CRAM de las muestras, uno o más VCF de verdad y un manifiesto de la ejecución.

La referencia inicial podrá ser simulada por el programa según la configuración
o proceder de un FASTA real suministrado como entrada. Un genoma podrá tener uno
o varios cromosomas o *scaffolds*. Un *scaffold* es una secuencia que organiza
y orienta fragmentos ensamblados y que puede contener regiones cuya secuencia
exacta no se conoce. La cantidad y el tamaño de estas secuencias, y por tanto el
tamaño total del genoma, podrán variar mediante la configuración cuando la
referencia inicial sea simulada. El único FASTA de cada referencia podrá
contener uno o varios de estos cromosomas o *scaffolds*.

Una referencia estructural es una nueva referencia creada a partir de otra
referencia disponible previamente. Solo se diferencia de su referencia de origen
por las duplicaciones, deleciones e inversiones solicitadas. El sistema debe
poder mantener disponibles dos o más referencias estructurales y conservar el
origen y los cambios de cada una.

La verdad conocida es la información exacta, definida durante la simulación,
sobre los genotipos reales de los individuos simulados. El VCF de verdad debe
recoger esos genotipos. Cada muestra corresponde a un individuo simulado. Su
estructura puede corresponder a cualquiera de las referencias disponibles,
pero su BAM o CRAM puede estar mapeado contra otra referencia disponible.

La ploidía será configurable. En una misma ejecución, todos los individuos
tendrán la misma ploidía.

Los individuos podrán contener SNP, inserciones o deleciones pequeñas (*indels*)
y variantes estructurales respecto a su referencia. Podrán crearse directamente,
mediante la evolución configurada de una población inicial o mediante cruces
entre dos o más parentales.

La coherencia de las salidas no exige que todos los archivos utilicen una misma
referencia. Para cada referencia de mapeo utilizada debe existir un conjunto
independiente: sus BAM o CRAM y su VCF de verdad deben estar referidos
exactamente al mismo genoma. Nunca se mezclarán en un conjunto archivos o datos
que utilicen referencias distintas.

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

## 4. Requisitos confirmados

### 4.1. Ficheros de configuración de entrada

Un fichero estructurado organiza la información en campos y grupos que un
programa puede identificar y leer, en lugar de presentarla como texto libre. No
se fija aquí su formato ni su esquema, es decir, las reglas que definirán sus
campos, tipos y relaciones. Esta decisión se tomará durante la planificación de
la implementación, después o al mismo tiempo que la elección del lenguaje, porque
depende de ella.

El fichero general contiene las decisiones que afectan al conjunto de la
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
- **REQ-CONF-06.** El fichero general debe permitir configurar la ploidía común
  de los individuos de la ejecución.
- **REQ-CONF-07.** El fichero general debe seleccionar una de las tres vías de
  creación definidas en la sección 4.3.
- **REQ-CONF-08.** El fichero general debe permitir definir cambios concretos o
  solicitar cambios al azar indicando su tipo y frecuencia.
- **REQ-CONF-09.** Cuando se seleccione la vía evolutiva o la vía de familias y
  poblaciones segregantes, el fichero general debe referenciar el fichero de
  configuración poblacional correspondiente.
- **REQ-CONF-10.** Cada escenario de la vía evolutiva o de la vía de familias y
  poblaciones segregantes debe tener su propio fichero de configuración
  poblacional estructurado y legible por el programa.
- **REQ-CONF-11.** El fichero de un escenario evolutivo debe definir cómo se
  obtiene la población final desde la inicial y qué propiedades y parámetros se
  aplican.
- **REQ-CONF-12.** El fichero de un escenario de familias o poblaciones
  segregantes debe indicar los parentales, los cruces, su orden o secuencia, las
  relaciones necesarias y las generaciones.
- **REQ-CONF-13.** Todos los parámetros configurables que se acuerden en
  versiones posteriores deben proporcionarse mediante el fichero general o el
  fichero poblacional que corresponda.
- **REQ-CONF-14.** La persona usuaria debe poder variar los parámetros sin
  modificar el código del programa.

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
  referencia de origen únicamente por las duplicaciones, deleciones e
  inversiones solicitadas.
- **REQ-REF-10.** El sistema debe conservar la relación entre cada referencia
  estructural y su referencia de origen, junto con los cambios que las
  distinguen.
- **REQ-REF-11.** El sistema debe permitir que haya dos o más referencias
  estructurales disponibles.

Estos requisitos no seleccionan las herramientas ni los tipos concretos de
índice. La necesidad de cada tipo dependerá de las operaciones y los formatos
que se utilicen durante la implementación.

### 4.3. Individuos, muestras y grupos

En este documento, una muestra es el conjunto de datos generado para un
individuo simulado. La ploidía es el número de juegos completos de cromosomas
representados en un individuo.

#### 4.3.1. Requisitos comunes

Una variante estructural es un cambio que afecta a la estructura de una región
del genoma. Esta versión no fija el tamaño que separará una variante estructural
de un *indel* ni su representación detallada.

- **REQ-IND-01.** El sistema debe crear individuos que contengan SNP, *indels*
  y variantes estructurales respecto a una referencia.
- **REQ-IND-02.** La estructura de cada individuo debe poder corresponder a
  cualquiera de las referencias disponibles.
- **REQ-IND-03.** El sistema debe conservar qué referencia disponible describe
  la estructura de cada muestra.
- **REQ-IND-04.** Los cambios deben poder crearse al azar configurando su tipo y
  frecuencia.
- **REQ-IND-05.** Los cambios también deben poder definirse de forma concreta en
  el fichero general de configuración.
- **REQ-IND-06.** El sistema debe poder simular diferentes ploidías.
- **REQ-IND-07.** En una misma ejecución, todos los individuos deben tener la
  ploidía indicada en el fichero general de configuración.
- **REQ-IND-08.** En esta primera aproximación, no se deben permitir cambios de
  ploidía entre individuos de una misma ejecución.

El significado estadístico de la frecuencia usada para crear cambios al azar
todavía no está definido. No debe suponerse si representa una probabilidad por
base, una cantidad esperada, una frecuencia alélica u otra magnitud.

No se han definido los valores de ploidía admitidos ni si la ploidía podrá
variar entre cromosomas. En cualquier caso, todos los individuos de una
ejecución compartirán la misma configuración de ploidía.

#### 4.3.2. Vía 1: creación directa de individuos

- **REQ-DIR-01.** El sistema debe permitir crear directamente individuos con
  los cambios indicados en el fichero general de configuración.
- **REQ-DIR-02.** La creación directa debe admitir tanto cambios concretos como
  cambios al azar configurados por tipo y frecuencia.

#### 4.3.3. Vía 2: población inicial y evolución configurada

El desequilibrio de ligamiento (LD) es la asociación no aleatoria entre alelos
de loci diferentes. El equilibrio de Hardy-Weinberg es una condición de
referencia que relaciona las frecuencias alélicas y genotípicas bajo unos
supuestos determinados. En este documento, el equilibrio o la desviación de
Hardy-Weinberg se trata como una propiedad configurable de la población, no como
una única magnitud numérica.

La endogamia es el apareamiento entre individuos emparentados. La tasa de
mutación describe con qué frecuencia aparecen cambios genéticos nuevos. Los
cuatro conceptos anteriores son ejemplos, no una lista cerrada de propiedades
o parámetros poblacionales. Los modelos y los parámetros exactos todavía no
están definidos.

- **REQ-EVO-01.** A partir de una referencia, el sistema debe permitir crear una
  población inicial con cambios.
- **REQ-EVO-02.** A partir de la población inicial, el sistema debe generar una
  población final aplicando las propiedades y los parámetros evolutivos
  indicados en el fichero de configuración poblacional del escenario.
- **REQ-EVO-03.** El fichero de configuración poblacional debe poder representar,
  entre otros ejemplos, LD, equilibrio o desviación de Hardy-Weinberg, endogamia
  y tasa de mutación.
- **REQ-EVO-04.** El sistema debe permitir que un escenario evolutivo incluya
  otras propiedades o parámetros poblacionales además de los ejemplos
  confirmados.
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

- **REQ-CRU-01.** El sistema debe permitir crear familias o poblaciones
  segregantes a partir de dos o más parentales.
- **REQ-CRU-02.** Esta vía debe admitir, entre otros ejemplos, pedigríes,
  poblaciones F2, retrocruzamientos, líneas de introgresión y poblaciones MAGIC.
- **REQ-CRU-03.** El sistema puede ofrecer tipos de cruce predefinidos.
- **REQ-CRU-04.** Cuando el cruce solicitado no sea uno de los predefinidos, la
  persona usuaria debe poder definirlo en el fichero de configuración
  poblacional del escenario.
- **REQ-CRU-05.** El fichero de configuración poblacional debe indicar los
  parentales, los cruces, su orden o secuencia, las relaciones necesarias, las
  generaciones y los demás parámetros del escenario.
- **REQ-POB-02.** El sistema debe permitir crear poblaciones segregantes.
- **REQ-POB-03.** El sistema debe permitir crear familias.

El contenido mínimo, la semántica, la validación y la capacidad expresiva del
fichero poblacional todavía no están definidos. Tampoco se ha fijado el catálogo
de cruces predefinidos ni un lenguaje para describir cruces complejos.

### 4.4. Datos por muestra

La referencia de mapeo es la referencia respecto a la que se registran los
alineamientos de una muestra en su BAM o CRAM. Puede ser distinta de la
referencia cuya estructura corresponde a la muestra.

La cabecera es la parte del archivo que contiene metadatos, es decir,
información que describe el propio archivo y su contexto. En los requisitos
siguientes, la cabecera y los metadatos deben ajustarse al estándar del formato,
aunque los campos y las versiones concretas todavía no se han decidido.

- **REQ-DAT-01.** El sistema debe generar archivos BAM o CRAM para las muestras.
- **REQ-DAT-02.** Los BAM o CRAM de una muestra deben poder generarse mapeados
  contra cualquiera de las referencias disponibles.
- **REQ-DAT-03.** La referencia de mapeo debe poder ser una referencia
  estructural distinta de la que corresponde a la muestra.
- **REQ-DAT-04.** Cada BAM o CRAM debe incluir en su propia cabecera o metadatos
  estándar la información necesaria para identificar sin ambigüedad la muestra,
  la referencia de mapeo y el conjunto al que pertenece.
- **REQ-DAT-05.** El sistema debe permitir introducir errores, contaminación,
  baja cobertura, mala calidad de mapeo y otras degradaciones en los archivos
  BAM o CRAM.

La expresión «otras degradaciones» todavía no identifica comportamientos
concretos. Deberá detallarse antes de convertirla en requisitos comprobables.
Tampoco se ha fijado cuántos BAM o CRAM se generarán por muestra cuando haya
varias referencias de mapeo posibles.

### 4.5. VCF de verdad y coherencia de las salidas

- **REQ-VER-01.** El sistema debe generar un VCF de verdad que represente los
  genotipos reales simulados de los SNP, *indels* y variantes estructurales de
  los individuos.
- **REQ-VER-02.** Cada genotipo del VCF de verdad debe corresponder al individuo
  simulado identificado en el archivo.
- **REQ-VER-03.** Cada VCF de verdad debe estar referido exactamente al mismo
  genoma usado como referencia de mapeo en los BAM o CRAM correspondientes.
- **REQ-VER-04.** Si se generan BAM o CRAM contra referencias distintas, el
  sistema debe generar un VCF de verdad independiente para cada una de esas
  referencias.
- **REQ-VER-05.** Cada conjunto debe contener únicamente BAM o CRAM y VCF de
  verdad que utilicen la misma referencia. No se deben mezclar datos referidos
  a genomas distintos.
- **REQ-VER-06.** Cada VCF de verdad debe incluir en su propia cabecera o
  metadatos estándar la información necesaria para identificar sin ambigüedad la
  referencia usada y el conjunto al que pertenece.
- **REQ-VER-07.** Los genotipos de todas las clases de variante deben expresarse
  respecto a la referencia usada por el conjunto correspondiente.

El VCF de verdad recoge los genotipos definidos por la simulación. Su finalidad
dentro de este proyecto es conservar la descripción exacta de lo simulado y
mantener la coherencia entre las salidas. Una misma muestra puede formar parte
de conjuntos diferentes si se mapea contra varias referencias, pero los
archivos de esos conjuntos deben permanecer separados.

### 4.6. Manifiesto de la ejecución

El manifiesto es un fichero estructurado de salida que registra lo que el
sistema utilizó y generó en una ejecución.

- **REQ-MAN-01.** El sistema debe generar un manifiesto para cada ejecución.
- **REQ-MAN-02.** El manifiesto debe ser un fichero estructurado y legible por el
  programa.
- **REQ-MAN-03.** El manifiesto debe identificar cuál de las tres vías de
  creación se utilizó.
- **REQ-MAN-04.** El manifiesto debe registrar el nombre o la ruta del fichero
  general y de cada fichero de configuración poblacional utilizado.
- **REQ-MAN-05.** El manifiesto debe conservar los parámetros realmente
  aplicados, incluida la ploidía común, el modo de creación de los cambios y los
  parámetros propios de la vía utilizada.
- **REQ-MAN-06.** El manifiesto debe indicar si la referencia inicial fue
  simulada o suministrada. Si fue suministrada, debe conservar el nombre o la
  ruta del FASTA de origen.
- **REQ-MAN-07.** El manifiesto debe registrar cada referencia, su referencia de
  origen cuando sea derivada y los cambios aplicados para obtenerla.
- **REQ-MAN-08.** El manifiesto debe registrar los nombres o las rutas de los
  archivos generados.

El formato y el esquema del manifiesto se elegirán durante la planificación de
la implementación, después o al mismo tiempo que el lenguaje. Esta decisión no
forma parte de la especificación funcional.

## 5. Flujo funcional confirmado

El flujo funcional describe la relación entre los resultados, no el orden
técnico en que deberá implementarse el sistema.

1. El programa lee el fichero general de configuración y obtiene de él la vía
   de creación seleccionada.
2. Si se ha seleccionado la vía evolutiva o la vía de familias y poblaciones
   segregantes, el programa obtiene del fichero general la referencia al
   fichero de configuración poblacional propio del escenario y lee este segundo
   fichero.
3. Si la referencia inicial debe simularse, el programa genera sus cromosomas o
   *scaffolds* con los tamaños y la composición configurados, y obtiene el tamaño
   total a partir de sus tamaños individuales. Si la referencia es suministrada,
   el programa utiliza el FASTA real indicado como entrada.
4. El sistema conserva si la referencia inicial fue simulada o suministrada y,
   cuando proceda, el archivo FASTA de origen.
5. Toda la secuencia de la referencia inicial se guarda en un único FASTA, que
   puede contener uno o varios cromosomas o *scaffolds*. No se genera ningún
   índice si las operaciones y los formatos posteriores no lo necesitan; en
   caso contrario, se generan todos los tipos necesarios.
6. A partir de una referencia ya disponible, se crea una nueva referencia
   estructural aplicando solo las duplicaciones, deleciones e inversiones
   solicitadas.
7. Toda la secuencia de la referencia derivada se guarda en un único FASTA. Se
   conservan su referencia de origen y los cambios aplicados, y se generan los
   tipos de índice que necesiten las operaciones o los formatos usados. Los
   pasos 6 y 7 pueden repetirse para disponer de dos o más referencias
   estructurales.
8. El fichero general indica si los cambios se definirán de forma concreta o se
   crearán al azar según su tipo y frecuencia.
9. El sistema aplica la vía seleccionada:

   - crea directamente individuos con los cambios configurados;
   - crea una población inicial desde una referencia y genera una población
     final según las propiedades y los parámetros definidos en el fichero de
     configuración poblacional del escenario; o
   - crea una familia o población segregante desde dos o más parentales según
     los cruces predefinidos o descritos por el usuario, su orden o secuencia,
     las relaciones necesarias y el número de generaciones indicados en el
     fichero de configuración poblacional del escenario.

10. Todos los individuos reciben la ploidía común de la ejecución. Para cada
   muestra se conserva la referencia que corresponde a su estructura.
11. Para generar un BAM o CRAM, se elige como referencia de mapeo cualquiera de
   las referencias disponibles. Esta referencia puede ser distinta de la que
   corresponde a la estructura de la muestra.
12. Se pueden introducir degradaciones en los BAM o CRAM, incluidas las
    categorías ya confirmadas.
13. Para la referencia de mapeo elegida se genera un VCF de verdad independiente
    con los genotipos simulados de SNP, *indels* y variantes estructurales,
    expresados respecto al mismo genoma que los BAM o CRAM correspondientes.
14. Cada BAM, CRAM y VCF incluye en su cabecera o metadatos estándar la
    información necesaria para identificar su referencia y su conjunto.
15. Si se usan varias referencias de mapeo, los pasos 11 a 14 producen conjuntos
    independientes. Los datos de referencias distintas no se mezclan.
16. El sistema genera el manifiesto con la vía utilizada, el fichero general y
    cada fichero de configuración poblacional utilizado, los parámetros
    realmente aplicados, el origen y los cambios de las referencias, y los
    nombres o las rutas de los archivos generados. El manifiesto es estructurado
    y legible por el programa.

Todas las salidas deben describir la misma simulación, aunque no utilicen una
misma referencia. Dentro de cada conjunto, todos los archivos deben usar
exactamente la misma referencia. El flujo termina cuando se han generado esas
salidas y el manifiesto conserva la información necesaria para relacionarlas.

## 6. Entradas y salidas confirmadas

### 6.1. Entradas

| Entrada | Contenido confirmado | Información pendiente |
|---|---|---|
| Fichero general de configuración | Fichero estructurado y legible por el programa con las decisiones generales de la ejecución: modo de obtención de la referencia, ploidía común, vía de creación y cambios concretos o aleatorios; para una referencia simulada, también número, tamaño y composición de sus secuencias. Cuando la vía lo requiere, referencia el fichero de configuración poblacional correspondiente | Interpretación de las frecuencias y otros parámetros generales que se acuerden; formato y esquema aplazados a la planificación de la implementación |
| Fichero de configuración poblacional | Un fichero estructurado y legible por el programa para cada escenario evolutivo o de familias y poblaciones segregantes. En la vía evolutiva describe la transformación de la población inicial en la final y sus propiedades y parámetros. En la vía de familias y poblaciones segregantes describe los parentales, los cruces, su orden o secuencia, las relaciones necesarias, las generaciones y los demás parámetros del escenario | Contenido mínimo, semántica, reglas de validación y capacidad expresiva; modelos concretos, cruces predefinidos y personalizados; formato y esquema aplazados a la planificación de la implementación |
| FASTA real de la referencia inicial | Secuencia real que se usará como referencia inicial cuando se elija el modo suministrado | Requisitos detallados, validación, tratamiento del archivo y datos de procedencia |

El FASTA real solo es una entrada cuando se elige una referencia inicial
suministrada. El fichero de configuración poblacional solo es una entrada en la
vía evolutiva o en la vía de familias y poblaciones segregantes.

### 6.2. Salidas

| Salida | Contenido confirmado | Información pendiente |
|---|---|---|
| FASTA de las referencias | Un único FASTA por cada referencia inicial o derivada; el archivo contiene uno o varios cromosomas o *scaffolds* | Ninguna decisión funcional pendiente sobre su organización |
| Índices de las referencias | Ningún índice, uno o varios por referencia, según los tipos que necesiten las operaciones o los formatos utilizados | Tipos de índice exigidos por cada operación y herramientas de indexación |
| Estructura de los individuos y grupos | Individuos con SNP, *indels* y variantes estructurales, con la misma ploidía configurada y creados por la vía directa, la vía evolutiva o la vía de familias y poblaciones segregantes | Valores de ploidía, posible variación entre cromosomas, significado de las frecuencias, modelos evolutivos y reglas de cruzamiento |
| BAM o CRAM de las muestras | Datos de cada muestra mapeados contra una referencia disponible, que puede diferir de la estructura de la muestra, y con posibilidad de degradación; su cabecera o metadatos identifican la muestra, la referencia y el conjunto | Número de archivos por muestra, tecnología de secuenciación, parámetros, campos de la cabecera y versión del formato |
| VCF de verdad por referencia | Genotipos reales simulados de SNP, *indels* y variantes estructurales, expresados contra la misma referencia que los BAM o CRAM de su conjunto | Representación detallada de *indels* y variantes estructurales, campos de la cabecera y versión del formato |
| Manifiesto de la ejecución | Fichero estructurado y legible por el programa con la vía utilizada, el fichero general y los ficheros de configuración poblacional utilizados, los parámetros realmente aplicados, la ploidía común, el origen y los cambios de las referencias, y los nombres o rutas de los archivos generados | Formato y esquema aplazados a la planificación de la implementación |

No se ha decidido si una ejecución producirá BAM, CRAM o permitirá producir
ambos. Tampoco se ha decidido si se conservarán a la vez archivos sin degradar
y archivos degradados. Cuando se seleccionen varias referencias de mapeo, cada
una tendrá su propio conjunto independiente de BAM o CRAM y VCF de verdad.

## 7. Criterios de aceptación preliminares

Un criterio de aceptación es una comprobación observable que permite decidir si
un requisito se ha cumplido. Los siguientes criterios pueden deducirse de los
requisitos confirmados. Sus casos de prueba y detalles concretos siguen
pendientes.

| Criterio | Evidencia preliminar | Aspectos pendientes |
|---|---|---|
| **ACE-01. Configuración sin cambiar el código** | El programa puede leer ficheros generales con distinto número, tamaño o composición de cromosomas o *scaffolds* y generar las referencias simuladas correspondientes sin modificar el código. Cuando se selecciona la vía evolutiva o la vía de familias y poblaciones segregantes, el fichero general referencia el fichero de configuración poblacional del escenario y el programa usa ambos | Formato, esquema y casos de prueba |
| **ACE-02. Tamaños del genoma simulado** | Cada cromosoma o *scaffold* tiene el tamaño configurado y el tamaño total corresponde a los tamaños de todas sus secuencias | Casos de prueba y comprobación exacta |
| **ACE-03. Modos de referencia inicial** | El sistema puede simular la referencia inicial desde la configuración o utilizar un FASTA real suministrado; conserva qué modo se usó y el archivo de origen cuando procede | Requisitos detallados del FASTA suministrado |
| **ACE-04. Un FASTA por referencia** | Cada referencia inicial o derivada conserva toda su secuencia en un único FASTA, que contiene uno o varios cromosomas o *scaffolds* | Casos de prueba |
| **ACE-05. Índices de las referencias** | Cada referencia tiene cero, uno o varios índices y dispone de todos los tipos que necesitan las operaciones o los formatos que la usan | Tipos exigidos por cada operación y herramientas de indexación |
| **ACE-06. Creación de referencias estructurales** | Se pueden crear dos o más referencias estructurales a partir de referencias disponibles previamente, incluida una referencia inicial suministrada | Casos de prueba, tamaños y reglas exactas |
| **ACE-07. Diferencias entre referencias** | Cada nueva referencia estructural difiere de su referencia de origen solo por las duplicaciones, deleciones e inversiones solicitadas | Casos de prueba y comprobación de cada cambio |
| **ACE-08. Origen de las referencias** | Para cada referencia estructural se pueden conocer su referencia de origen y los cambios que las distinguen | Forma de representar y consultar la relación |
| **ACE-09. Clases de variante** | Se pueden crear individuos con SNP, *indels* y variantes estructurales respecto a cualquiera de las referencias disponibles | Límite entre *indels* y variantes estructurales y casos de prueba |
| **ACE-10. Creación de cambios** | La configuración permite definir cambios concretos o solicitar cambios al azar por tipo y frecuencia; el modo y los parámetros utilizados quedan registrados | Significado y control de la frecuencia aleatoria |
| **ACE-11. Ploidía** | Dos ejecuciones pueden usar ploidías configuradas diferentes sin cambiar el código; dentro de cada ejecución, todos los individuos tienen la misma ploidía configurada | Valores admitidos y posible variación entre cromosomas |
| **ACE-12. Creación directa** | Se pueden crear directamente individuos con cambios concretos o aleatorios indicados en la configuración | Casos de prueba |
| **ACE-13. Vía evolutiva** | Cada escenario usa su propio fichero de configuración poblacional para obtener una población final a partir de una inicial. El fichero permite configurar LD, condición de Hardy-Weinberg, endogamia y tasa de mutación como ejemplos no exhaustivos, además de otras propiedades y parámetros que admita el sistema | Contenido mínimo, semántica, validación, modelos, parámetros y comprobaciones de cada propiedad |
| **ACE-14. Familias y poblaciones segregantes** | Cada escenario usa su propio fichero de configuración poblacional para describir dos o más parentales, los cruces, su orden o secuencia, las relaciones necesarias y las generaciones. El sistema puede usar cruces predefinidos o cruces más complejos definidos por el usuario; pedigríes, F2, retrocruzamientos, líneas de introgresión y poblaciones MAGIC son ejemplos no exhaustivos | Contenido mínimo, semántica, validación, capacidad expresiva, cruces predefinidos y personalizados, y casos de prueba |
| **ACE-15. Referencia de mapeo** | Se puede generar un BAM o CRAM de una muestra contra cualquiera de las referencias disponibles, incluida una referencia estructural distinta de la muestra | Cardinalidad de archivos y comprobaciones técnicas |
| **ACE-16. Cabecera de BAM y CRAM** | La cabecera o los metadatos estándar de cada BAM o CRAM permiten identificar sin ambigüedad la muestra, la referencia de mapeo y el conjunto | Campos concretos y versiones de los formatos |
| **ACE-17. Degradaciones** | Se pueden generar BAM o CRAM afectados por errores, contaminación, baja cobertura y mala calidad de mapeo | Definición, intensidad y forma de comprobar cada degradación |
| **ACE-18. VCF de verdad** | Para cada referencia de mapeo usada existe un VCF de verdad independiente con los genotipos simulados de SNP, *indels* y variantes estructurales expresados contra esa misma referencia | Representación detallada de *indels* y variantes estructurales y casos de concordancia |
| **ACE-19. Cabecera del VCF** | La cabecera o los metadatos estándar de cada VCF permiten identificar sin ambigüedad la referencia y el conjunto | Campos concretos y versiones del formato |
| **ACE-20. Separación por referencia** | Ningún conjunto mezcla BAM, CRAM o VCF referidos a genomas distintos; si una muestra se mapea contra varias referencias, sus conjuntos permanecen independientes | Forma de organizar los conjuntos |
| **ACE-21. Coherencia de las salidas** | Las referencias, los individuos y cada conjunto de BAM o CRAM y VCF de verdad se pueden relacionar sin ambigüedad y describen la misma simulación | Identificadores y comprobaciones de consistencia |
| **ACE-22. Manifiesto** | Cada ejecución genera un manifiesto estructurado que identifica la vía utilizada, el fichero general y todos los ficheros de configuración poblacional utilizados, y conserva los parámetros realmente aplicados, la ploidía común, el origen y los cambios de las referencias, y los nombres o rutas de todos los archivos generados | Formato, esquema y casos de prueba |

Estos criterios son cualitativos. No deben interpretarse como pruebas cerradas
hasta que se resuelvan los aspectos pendientes.

## 8. Decisiones abiertas

1. **FASTA real suministrado.** ¿Qué requisitos detallados deberá cumplir y cómo
   se tratará o conservará el archivo de entrada?
2. **Número de conjuntos por muestra.** ¿Se generará un conjunto para una sola
   referencia seleccionada, para varias referencias seleccionadas o para todas
   las referencias disponibles? En cualquier caso, cada referencia tendrá un
   conjunto independiente.
3. **Frecuencia de los cambios al azar.** ¿Qué significa la frecuencia para cada
   tipo de cambio y cómo se controlará que la simulación respeta el valor
   configurado?
4. **Configuración poblacional evolutiva.** ¿Cuál será su contenido mínimo, qué
   significado tendrá cada elemento y qué reglas de validación deberá cumplir?
   También deben concretarse los modelos y parámetros para representar LD,
   equilibrio o desviación de Hardy-Weinberg, endogamia, tasa de mutación y
   cualquier otra propiedad evolutiva que se decida admitir.
5. **Configuración poblacional de familias y poblaciones segregantes.** ¿Cuál
   será su contenido mínimo, qué significado tendrá cada elemento y qué reglas
   de validación deberá cumplir? Debe concretarse la capacidad expresiva
   necesaria para describir parentales, cruces, orden o secuencia, relaciones y
   generaciones, incluidos los cruces más complejos definidos por el usuario.
6. **Cruces predefinidos y personalizados.** ¿Qué tipos de cruce se ofrecerán
   inicialmente como predefinidos y cómo se distinguirán de los definidos por
   el usuario? La respuesta no debe convertir los ejemplos actuales en un
   catálogo cerrado.
7. **Representación de *indels* y variantes estructurales.** ¿Cómo se
   representarán en el VCF de verdad, incluido el límite que distingue ambas
   clases?
8. **Resultado mínimo de una simulación.** ¿Qué conjunto mínimo de datos debe
   producir una ejecución para considerarse completa?
9. **Valores de ploidía.** ¿Qué valores de ploidía debe admitir el sistema?
10. **Ploidía por cromosoma.** ¿La ploidía común de una ejecución se aplicará
   por igual a todos los cromosomas o podrá variar entre cromosomas? En ambos
   casos, todos los individuos de la ejecución deberán compartir la misma
   configuración de ploidía.
11. **Contexto biológico.** ¿Qué organismos, tamaños de genoma y escalas de
   población deben admitirse?
12. **Generación de BAM y CRAM.** ¿Qué tecnologías de secuenciación y qué
    características de los datos deben simularse?
13. **Degradaciones.** ¿Cómo se define y controla cada degradación, y qué otras
    degradaciones deben incluirse?
14. **Formatos BAM y CRAM.** ¿Una ejecución producirá BAM, CRAM o permitirá
    producir ambos?
15. **Datos sin degradar y degradados.** ¿Deben conservarse ambas versiones?
16. **Cabeceras y metadatos.** ¿Qué campos concretos y qué versiones de los
    formatos se usarán para identificar la referencia, la muestra y el conjunto?

## 9. Decisiones técnicas aplazadas a la planificación de la implementación

**Formato y esquema de los ficheros de configuración y del manifiesto.** El
formato y el esquema del fichero general, de los ficheros de configuración
poblacional y del manifiesto se elegirán durante la planificación de la
implementación, después o al mismo tiempo que el lenguaje, porque dependen de
esa elección. Todos deberán seguir siendo estructurados y legibles por el
programa. Tampoco se fija todavía un lenguaje o una sintaxis para describir los
cruces personalizados; su significado y capacidad expresiva son decisiones
funcionales que deben resolverse antes de elegir su representación técnica.

**Índices exigidos por cada operación.** Durante esa planificación se
identificarán los tipos concretos de índice que exige cada operación o formato
utilizado y se elegirán las herramientas de indexación. Esta decisión no cambia
que cada referencia pueda tener cero, uno o varios índices según sus
necesidades.

Esta especificación no elige formatos, índices, lenguaje, herramientas,
algoritmos, bibliotecas ni arquitectura. Las decisiones funcionales que siguen
abiertas se enumeran por separado en la sección 8.

## 10. Próximos pasos de diseño

1. Definir los requisitos detallados del FASTA real suministrado, sin cambiar
   las dos formas ya acordadas de obtener la referencia inicial.
2. Definir el significado y el control de las frecuencias usadas para crear
   cambios al azar.
3. Definir el contenido mínimo, la semántica y las reglas de validación del
   fichero de configuración poblacional evolutiva, incluidos sus modelos,
   propiedades y parámetros.
4. Definir el contenido mínimo, la semántica, las reglas de validación y la
   capacidad expresiva del fichero de configuración poblacional de familias y
   poblaciones segregantes.
5. Decidir qué cruces se ofrecerán inicialmente como predefinidos y qué debe
   poder expresar un cruce personalizado, sin cerrar el catálogo de escenarios
   posibles ni elegir todavía su lenguaje o esquema técnico.
6. Definir la representación detallada de *indels* y variantes estructurales en
   el VCF de verdad.
7. Decidir para cuántas referencias se generarán conjuntos independientes por
   muestra.
8. Definir los valores de ploidía admitidos y decidir si la configuración podrá
   variar entre cromosomas, manteniéndola igual entre individuos.
9. Definir los campos de las cabeceras y las versiones de BAM, CRAM y VCF que
   identificarán cada referencia, muestra y conjunto.
10. Definir el resultado mínimo de una simulación completa.
11. Concretar las demás entradas y propiedades biológicas pendientes.
12. Completar los criterios de aceptación con casos de prueba.
13. Tomar las decisiones de implementación cuando el comportamiento esperado
   esté suficientemente definido.
