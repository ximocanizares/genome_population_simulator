# Especificación, arquitectura mínima y plan de validación externa FASTQ→VCF 0.1

## 1. Estado y autoridad

- **Fecha:** 4 de septiembre de 2026.
- **Estado:** contrato VE-SNP 0.1 cerrado tras la revisión adversarial de §1.2.
  VE.0 está cerrado y VE.1.1 puede empezar; la ejecución externa aún no se ha
  iniciado y VE.2 (haplotipos, LD y fase) sigue aplazada.
- **Finalidad:** definir una validación externa que parta únicamente de los
  FASTQ, la referencia y el VCF de verdad publicados por el simulador.
- **Relación con el simulador:** trabajo paralelo a F5. No modifica F5, no
  presupone que F5 esté terminado y no forma parte de su cierre.
- **Autoridad normativa:** banco externo DTF-024, con coordenada `VE`
  (*validación externa*). La especificación vigente excluye la llamada de
  variantes y la evaluación de herramientas externas del runtime del simulador;
  este plan concreta ese banco sin modificar dicho runtime.
- **Coordenada:** `VE`. No se usa `F9.5`, porque F9 depende de
  todas las fases que entren en la versión y este trabajo se ha solicitado en
  paralelo a F5.

### 1.1. Fuentes inspeccionadas

Este borrador conserva las decisiones de:

- `design/especificaciones-0.1.md`, en especial §§3, 4.4, 4.6, 4.8 y DTF-024;
- `design/arquitectura-0.1.md`, en especial las fronteras entre verdad,
  evidencia técnica, adaptadores externos y resultados;
- `design/hoja-de-ruta-0.1.md`, F1–F5 y F9.5;
- los planes F1–F5, con atención a F1.6, F2.4–F2.5, F3.1/F3.5 y F4;
- los contratos de referencia reducida, `truth.jsonl`, VCF 4.5, FASTQ PE y
  manifiestos presentes en `schemas/`;
- los productores y pruebas actuales de reducción, VCF de verdad y FASTQ PE.

La implementación F5 sin confirmar que aparece en el árbol de trabajo se ha
tratado como trabajo ajeno en curso. No se usa como evidencia de cierre ni se
propone modificarla.

### 1.2. Revisión adversarial del 4 de septiembre de 2026

Una revisión adversarial contrastó este plan contra el código, los esquemas, las
decisiones normativas del proyecto y las herramientas realmente instaladas,
ejecutando FreeBayes, samtools y bcftools locales para no especular. Registró
siete hallazgos graves, once medios y siete menores. Los que exigen reabrir una
decisión se resuelven uno a uno; hasta cerrarlos todos, VE.1.1 no empieza.

| ID | Hallazgo | Estado |
|---|---|---|
| VE-R01 | El comando congelado fusiona SNP cercanos en registros MNP que VE-D06 no sabía comparar y VE-D11 obligaba a contar como falso positivo | **Cerrado** por VE-D13 |
| VE-R02 | VE-D07 prohíbe filtrar por MAPQ, pero el comando aplica el filtro implícito `-m 1` de FreeBayes | **Cerrado** por VE-D07 y §9 |
| VE-R03 | VE-D12 exige identidad binaria como puerta, contra DTC-038 y REQ-CONF-115, y omite sus precondiciones | **Cerrado** por VE-D12 y §8.1 |
| VE-R04 | La geometría de la referencia reducida no está fijada; la que produce el código son trece contigs mosaico con rampas de cobertura en sus extremos | **Cerrado** por VE-D14 y §5.4 |
| VE-R05 | «Igualdad lógica de BAM/VCF» se usa tres veces como criterio y no se define en ninguna parte | **Cerrado** por VE-D16 y §8.2 |
| VE-R06 | El mecanismo de importación A→B no existe como se describe y la referencia derivada reutiliza los nombres de contig del tomate original | **Cerrado** por VE-D15 y §7.1, con una corrección: el simulador ya valida identificador, orden y longitud, de modo que el riesgo se limita al flujo externo |
| VE-R07 | VE-D08 fija un criterio de éxito que la combinación VE-D07 + PE150 + referencia mosaico hace inalcanzable | **Cerrado** por VE-D08 y §6.4.2 |

Los siete graves quedaron cerrados el 5 de septiembre de 2026 mediante VE-D13 a
VE-D16 y las secciones §5.4, §6.4.1, §6.4.2, §7.1, §8.1 y §8.2. Dos se cerraron
corrigiendo el hallazgo, no solo aceptándolo: el simulador ya valida
identificador, orden y longitud de cada secuencia, y dos ejecuciones de
FreeBayes con las mismas entradas resultaron idénticas byte a byte.

Los medios y menores se resolvieron en la misma pasada: se nombró el generador
F2.2 y sus claves en VE-D04, con la banda de recuento segregante esperada; se
corrigió el texto de VE.1.1 que declaraba ausente un FreeBayes instalado; se
declaró la ruta `schema_version = "2"`, `set -o pipefail`, la biblioteca única
por individuo, los denominadores de cobertura y el alcance real de la
exoneración de DTC-066; y se sustituyeron los términos sin definir. VE.0.3
volvió a cerrarse con los ocho ejemplos manuales de §6.4.3.

## 2. Resultado que se quiere demostrar

Un flujo ajeno al simulador debe poder:

1. verificar que referencia, FASTQ y VCF de verdad pertenecen a la misma
   ejecución y a las mismas muestras;
2. alinear con minimap2 los FASTQ Illumina paired-end contra la referencia
   exacta publicada por el simulador;
3. producir alineamientos por muestra que `samtools quickcheck` acepte, que
   `samtools index` indexe y cuyas cabeceras se resuelvan según §7.1;
4. llamar conjuntamente las variantes de las 20 muestras con FreeBayes;
5. comparar el VCF llamado con el VCF de verdad mediante una equivalencia
   semántica acordada, no por identidad byte a byte;
6. evaluar por separado SNP independientes y la conservación observable de la
   señal de haplotipos o desequilibrio de ligamiento (LD). **Fuera de VE 0.1:**
   este punto pertenece a VE.2, que sigue aplazada; el cierre de VE 0.1 no lo
   demuestra ni puede citarse como si lo hiciera.

El **desequilibrio de ligamiento (LD)** es la asociación no aleatoria entre
alelos de loci distintos. Su medida y su interpretación dependen de si se usan
haplotipos faseados o genotipos sin fase; ambas evaluaciones no son
intercambiables.

## 3. Qué significa «recomponer el genoma de referencia»

La expresión admite tres significados distintos. El ensayo debe nombrar cuál
demuestra:

| Significado | Qué observa | Alcance de este borrador |
|---|---|---|
| Cobertura y alineamiento | Las lecturas pueden situarse sobre la referencia y cubrir sus bases según las métricas y denominadores que fija §8.3 | Incluido |
| Consenso por muestra | Se construye una secuencia a partir de alineamientos y genotipos, tomando decisiones ante heterocigosis, falta de cobertura y fase | Excluido de la primera prueba; necesita otro contrato |
| Recuperación de variantes | El VCF llamado recupera alelos y genotipos verdaderos respecto a la misma referencia | Incluido y constituye la demostración principal |

La primera prueba no debe presentarse como ensamblaje *de novo*. Illumina PE de
150 bases, una referencia conocida y un alineador de lecturas permiten un
ensayo de mapeo y llamada de variantes. No demuestran por sí solos que pueda
reconstruirse sin referencia un genoma completo, ni que pueda recuperarse la
fase de largo alcance de cada haplotipo.

Tampoco se exige un consenso en esta versión. En una muestra diploide, un único
consenso obliga a decidir cómo representar dos alelos y qué hacer con regiones
sin cobertura; esa decisión no es necesaria para comprobar la recuperación del
VCF verdadero.

## 4. Decisiones ya confirmadas para el ensayo

| Aspecto | Decisión confirmada |
|---|---|
| Datos de entrada | Cadena de dos ejecuciones del simulador: A publica la referencia reducida; B la importa y publica FASTQ y VCF de verdad contra esa referencia exacta |
| Referencia | 500.000 bases derivadas de la referencia de tomate disponible. Genoma fuente y geometría de la reducción son parámetros de cada ensayo, no constantes del contrato; VE-D14 fija los del primero |
| Cohorte | 20 individuos |
| Cobertura solicitada | 30× por individuo, con el denominador haploide del simulador |
| Lecturas | Illumina paired-end |
| Alineador | minimap2 |
| Detector de variantes | FreeBayes |
| Primera prueba | Sin contaminación, quimeras, errores de demultiplexado ni otros artefactos; se solicita además ausencia de errores de secuenciación |
| Casos biológicos | Un caso de SNP independientes y otro de haplotipos/LD |
| Ejecución | Flujo externo; no incorporar caller ni comparador al runtime del simulador |

La carga nominal es de 300 millones de bases emitidas:
`20 individuos × 500.000 bases × 30`. Con PE150 equivale aproximadamente a
50.000 parejas por individuo y un millón de parejas en total, antes de aplicar
las reglas exactas de redondeo del simulador. Estos valores sirven para
dimensionar el ensayo, no para sustituir los contadores publicados en el
manifiesto.

## 5. Inventario confirmado y límites actuales

### 5.1. Capacidades reutilizables

- F1 puede derivar una referencia de tamaño exacto desde un FASTA real y
  publicar FASTA, mapa de coordenadas, hashes y manifiesto. Existe una
  configuración local de tomate cuyo FASTA fuente tiene 759 MB y SHA-256
  `52d3ebb862265776fb5dbe5dc7f6e57184bef6c0b69d5f845f67773d76f24233`.
- F2 publica `outputs/truth.jsonl` y, cuando se solicita, exactamente un VCF de
  verdad en `outputs/truth.vcf` o `outputs/truth.vcf.gz` con su índice.
- El VCF de verdad es VCF 4.5 multimuesta, conserva genotipos `GT` faseados y
  puede agrupar varios ALT en un registro multialélico.
- F4 publica por individuo
  `outputs/fastq/<event_id>/<individual_id>_R1.fastq.gz` y `_R2.fastq.gz`, más
  tamaños, SHA-256, cobertura y manifiesto de evento.
- Las pruebas actuales demuestran lectura externa de FASTQ, publicación VCF y
  coherencia interna con `truth.jsonl`; no ejecutan minimap2 ni un caller.

### 5.2. Herramientas observadas el 4 de septiembre de 2026

| Herramienta | Estado local observado | Consecuencia |
|---|---|---|
| minimap2 | `2.26-r1175`, disponible | Permite preparar el mapeo, sujeto a aprobar preset y parámetros |
| samtools | `1.21`, disponible | Permite ordenar, indexar y comprobar BAM |
| bcftools | `1.21`, disponible | Es anterior a la versión 1.24 exigida por la puerta externa DTC-066 del simulador |
| bgzip/tabix | HTSlib `1.19`, disponibles | Son anteriores a la versión 1.24 exigida por DTC-066 |
| FreeBayes | `v1.3.10`, disponible | VE.1.1 debe congelar versión, ruta y SHA-256 del ejecutable antes de ejecutar |
| WhatsHap/SHAPEIT4/Beagle | no disponibles | No debe darse por aprobado ningún faseador |
| hap.py/RTG vcfeval/truvari | no disponibles | Debe decidirse el comparador antes de prometer equivalencia compleja |

El simulador no instala ni descarga herramientas. El banco VE usa las
herramientas locales observadas y registra sus identidades en VE.1.1.

### 5.3. Escenario PE neutro disponible

F4 ya admite `technical_event.sequencing_errors = "disabled"` (commit
`07cde2c`). Este modo PE emite exactamente cero sustituciones en ambos mates,
conserva las calidades del perfil y declara explícitamente que no representan
una probabilidad de error realizada. La configuración efectiva y el manifiesto
son coherentes con ese modo; el manifiesto registra oportunidades y cero errores
por mate.

VE-SNP debe usar ese modo y fijar a cero quimeras, contaminación e *index
hopping*. Así queda cerrado el control técnico neutro sin usar
`error_scale = 0` ni un Phred infinito.

### 5.4. Geometría que produce la reducción

`config-general-v2` exige siempre `region_size_bases` y admite `target_bases` o
`reduction_percent`, nunca ambos. Solo `target_bases` da una longitud exacta;
`reduction_percent` la deriva de un porcentaje y no garantiza una cifra redonda.

`region_size_bases` no es un detalle de tamaño: decide la estructura del
artefacto. El número de bloques es `target_bases / region_size_bases`
(`src/reference_reduction.rs`, `preflight_reduced_reference`), y
`allocate_proportional` los reparte entre las secuencias fuente en proporción a
su capacidad. Después, `build_reduced_reference`
(`src/reference_reduction_io.rs`) concatena en **una sola** secuencia derivada
todos los bloques que provienen de la misma secuencia fuente. El resultado no es
un locus contiguo: es un contig por secuencia fuente que recibió bloques, y cada
uno es un mosaico de tramos no adyacentes en el genoma original.

Con el tomate `SL4.00`, sus trece secuencias y los valores del primer ensayo
—`region_size_bases = 10000`, `target_bases = 500000`—, una estimación con
reparto proporcional da cincuenta bloques y trece contigs derivados de 10.000 a
60.000 bases. El reparto exacto lo decide `allocate_proportional`, de modo que
la cifra definitiva es la que publique VE.1.2, no esta estimación.

`spatial_mode="uniform"` pondera por base, así que los mil loci se reparten en
proporción a la longitud de cada contig derivado: alrededor de veinte en uno de
10.000 bases y ciento veinte en uno de 60.000. No se estratifica ni se corrige,
pero el informe publica las métricas por contig (§8.3) para que ese reparto sea
visible, porque los efectos de borde pesan más en los contigs cortos.

Lo que importa para el ensayo es la consecuencia: trece contigs son veintiséis
extremos, y en cada extremo la profundidad decae porque el muestreo solo puede
situar un fragmento completo dentro del contig. Con lecturas de 150 bases eso
afecta a unas 3.900 bases, en torno al 0,78 % de la referencia, donde caerían
unos ocho de los mil SNP. De ahí que §6.4.2 tenga las causas `contig-edge` y
`zero-depth`, y que VE-D14 exija publicar la geometría realizada en lugar de
suponerla.

## 6. Especificación funcional propuesta

### 6.1. Entradas obligatorias

Cada ejecución de validación debe recibir mediante un manifiesto de ensayo:

- ruta, tamaño y SHA-256 de la referencia FASTA derivada;
- identidad de la referencia y, para auditar su procedencia, FASTA fuente,
  configuración de reducción y mapa de coordenadas;
- rutas, tamaños y SHA-256 de los FASTQ R1/R2 de cada individuo;
- ruta, tamaño y SHA-256 del VCF de verdad y su índice cuando exista;
- manifiestos general y de evento y configuración efectiva del simulador;
- identificadores de las 20 muestras y correspondencia exacta
  individuo→R1/R2→columna VCF;
- versión del simulador, semilla y huella del estado del repositorio que produjo
  el conjunto;
- identidad, configuración, semilla y manifiesto de la ejecución A de reducción,
  y los mismos enlaces de la ejecución B que importó su FASTA; los FASTQ y VCF
  solo pertenecen conjuntamente a B;
- ejecutables, versiones y argumentos efectivos del flujo externo.

Una entrada con referencia, muestra, hash o cardinalidad incompatible debe
detenerse antes del alineamiento.

### 6.2. Dos conjuntos de evaluación separados

1. **VE-SNP:** SNP bialélicos independientes. Debe aislar la recuperación de
   alelos y genotipos sin atribuir un fallo a fase, multialelia o *indels*.
2. **VE-LD:** bloques de SNP con haplotipos y LD conocidos. Solo podrá generarse
   desde una salida F5 que haya superado sus propios contratos y evidencias. El
   banco externo podrá prepararse en paralelo, pero no fabricará una sustitución
   local de F5 ni usará sus cambios sin revisar como verdad aprobada.

Ambos conjuntos deben usar la misma referencia derivada y los mismos tamaños de
cohorte y cobertura, salvo que una decisión posterior autorice un control más
pequeño. Deben tener semillas e identidades distintas para no confundir sus
artefactos.

### 6.3. Resultados obligatorios

El flujo debe conservar, al menos:

- índice de minimap2 y su identificación;
- un BAM coordinado e indexado por individuo, con grupo de lectura y `SM`
  coherentes con el identificador biológico;
- tabla canónica biyectiva `individual_id → R1/R2 → RG.ID → RG.SM → RG.LB →
  BAM`, comprobada contra las cabeceras y contra cada alineamiento;
- informes por muestra de lecturas totales, mapeadas, primarias, secundarias,
  suplementarias, correctamente emparejadas y no mapeadas;
- distribución de profundidad por base e intervalo, fracción de bases cubiertas
  a ≥1× y ≥10×, mediana, cuantiles, coeficiente de variación, y número/longitud
  de intervalos con cobertura cero;
- número y proporción de pares no mapeados, además de su identidad y causa
  cuando la herramienta pueda informarla;
- VCF bruto de FreeBayes y su salida normalizada para comparación;
- tablas de verdaderos positivos, falsos negativos y falsos positivos por
  locus, alelo, genotipo y muestra;
- resultados de LD o fase claramente etiquetados por el estimador utilizado;
- tiempos, tiempo de CPU, memoria residente máxima y almacenamiento por
  etapa;
- un manifiesto final con hashes, versiones, comandos, estado y diagnósticos.

Los BAM y el VCF llamado son evidencia externa. Nunca deben sustituir ni
modificar `truth.jsonl` o el VCF de verdad del simulador.

### 6.4. Equivalencia VCF

La recuperación no puede significar igualdad byte a byte. FreeBayes puede usar
otra versión VCF, otros campos `INFO`, `QUAL` y `FILTER`, y no reproducirá
`GPS_CASE`. La comparación debe acordar:

- la misma referencia y el mismo conjunto ordenado de muestras;
- equivalencia de `CHROM`, posición y alelos después de la normalización
  aprobada;
- correspondencia de genotipos por muestra mediante el multiconjunto de alelos
  o dosis, separada de datos ausentes y de la fase;
- tratamiento explícito de registros multialélicos y *indels*;
- ninguna región se excluye por llamabilidad o mapeo ambiguo: VE-D07 lo prohíbe
  y esas condiciones se registran como causa acreditada de §6.4.2;
- campos del caller que son informativos pero no forman parte de la verdad.

El VCF de verdad debe permanecer inmutable. La normalización produce una vista
temporal o derivada identificada; nunca repara el original.

La métrica principal de genotipo ignora el separador de fase: `0|1`, `1|0` y
`0/1` son el mismo heterocigoto porque tienen la misma dosis. `./.` es una no
llamada. La orientación y el separador de fase solo entran en una tabla de fase
separada, y únicamente si el llamado publica fase local explícita y trazable.

Aunque VE-SNP 0.1 solo genere verdad bialélica, el inventario incluye cada ALT
de todos los registros llamados. Un ALT adicional en un locus verdadero, un
SNP sin verdad, un *indel* o un registro multialélico sin correspondencia se
clasifica como falso positivo técnico; no se elimina por estar fuera del alcance
biológico. La correspondencia se establece después de aplicar §6.4.1, no sobre
el registro tal como lo emite el llamador. En el piloto PE neutro, un falso
positivo sin causa acreditada de §6.4.2 detiene el cierre.

### 6.4.1. Descomposición de registros multinucleotídicos

FreeBayes llama por haplotipos locales. Con su valor predeterminado
`--haplotype-length 3`, dos SNP verdaderos separados por tres bases o menos no
se emiten como dos registros bialélicos, sino como un único registro
multinucleotídico (MNP) que abarca el tramo completo. Se ha comprobado
ejecutando FreeBayes v1.3.10 sobre datos construidos: dos SNP a distancia dos
producen `AGC → CGA` en una sola línea, donde la base intermedia no cambia.

Esta fusión es una propiedad real del llamador ante evidencia real: las dos
variantes viajan en las mismas lecturas. Imponer una distancia mínima entre los
loci verdaderos para evitarla falsearía el ensayo, porque eliminaría de la
prueba justo el caso en que dos variantes comparten fragmento. VE-SNP 0.1 no
impone separación mínima entre loci; el comparador descompone.

La descomposición es mecánica y determinista. Para cada ALT de un registro
llamado cuyo `REF` y `ALT` tengan la misma longitud, se recorren sus bases en
paralelo y se emite un SNP derivado por cada posición en que difieren, con la
coordenada `POS + desplazamiento`, la base `REF` y la base `ALT` de esa
posición. Las posiciones en que `REF` y `ALT` coinciden son anclaje del
haplotipo, no variantes, y no generan nada. Un ALT cuya longitud difiera de la
de `REF` no se descompone: es un alelo complejo o con *indel*, queda fuera del
alcance biológico de VE-SNP 0.1, se inventaría íntegro como falso positivo
técnico y se investiga.

El genotipo de cada SNP derivado se obtiene proyectando el `GT` del registro
sobre el alelo del que procede: un índice de alelo en el `GT` aporta la base
alternativa de ese ALT en esa posición, y el índice cero aporta la base de
referencia. Así, los dos casos que produce una verdad faseada quedan cubiertos:

- **En *cis*.** Dos heterocigotos en la misma copia dan un solo haplotipo
  alternativo, es decir un MNP bialélico con `GT` heterocigoto. Descompone a dos
  SNP heterocigotos.
- **En *trans*.** Dos heterocigotos en copias distintas dan dos haplotipos
  alternativos, es decir un MNP multialélico con `GT 1/2`. Descompone también a
  dos SNP heterocigotos, uno por alelo. Es el caso donde la descomposición
  aporta más: sin ella, el registro se descartaría por multialélico pese a
  contener exactamente la verdad.

Un SNP derivado se compara contra la verdad con las mismas reglas de VE-D06 que
un SNP llamado directamente. La clasificación del registro completo es la
conjunción de las de sus SNP derivados: si todos concuerdan y ninguna posición
cambiada carece de verdad, el registro es una concordancia con representación
alternativa; si falta una verdad, esa verdad es un falso negativo; si una
posición cambiada no tiene verdad, es un falso positivo. La tabla de resultados
conserva la representación de origen de cada SNP —directa o derivada de MNP— y
publica el recuento de registros fusionados, para que una tasa de fusión
anómala sea visible aunque las métricas de SNP salgan perfectas.

### 6.4.3. Ejemplos manuales de clasificación

Estos casos fijan la clasificación de cada situación prevista y son la condición
de cierre de VE.0.3. Cada uno admite una sola lectura bajo VE-D06, VE-D08 y
VE-D13; si alguno admitiera dos, el contrato aún no estaría cerrado. La muestra
es `s1` y el contig `c1`.

**E1 — Verdadero positivo.** Verdad `100 A>G`, `s1` `0|1`. Llamada
`c1 100 . A G` con `s1 0/1`. Es TP: VE-D06 compara dosis e ignora el separador
de fase, de modo que `0|1` y `0/1` son el mismo heterocigoto.

**E2 — Falso negativo con causa acreditada.** Verdad `200 C>T`, `s1` `0|1`. No
hay registro llamado en 200. Es FN y permanece en el denominador: VE-D07 impide
retirarlo. Se clasifica `unmappable` solo si las lecturas que solapan el locus
tienen MAPQ 0 en el BAM, dato que hay que exhibir. Sin una causa de la lista
cerrada acreditada por su dato, detiene el cierre.

**E3 — Falso positivo sin causa.** No hay verdad en 300. Llamada
`c1 300 . G A` con `s1 0/1`. Es FP. Sus longitudes `REF` y `ALT` coinciden, así
que no es `complex-allele`, y ninguna otra causa lo explica: detiene el cierre y
se investiga como defecto. En el control PE neutro un FP sin causa es
exactamente lo que el ensayo existe para detectar.

**E4 — Genotipo discordante.** Verdad `400 A>G`, `s1` `1|1`. Llamada
`c1 400 . A G` con `s1 0/1`. El sitio se recupera, así que no es FN ni FP: es
una discordancia de genotipo, que resta en concordancia de GT y no en
sensibilidad de sitio. Las dos magnitudes se publican por separado.

**E5 — Fusión MNP en *cis*.** Verdad `500 A>C` y `502 C>A`, ambos `s1` `0|1` en
la misma copia. Llamada `c1 500 . AGC CGA` con `s1 0/1`. VE-D13 descompone:
posición 500 `A>C`, posición 501 sin cambio —anclaje, no genera nada— y posición
502 `C>A`. Son **dos TP** con representación derivada de MNP. Sin la
descomposición serían dos FN más un FP.

**E6 — Fusión MNP multialélica en *trans*.** Verdad `600 A>C` y `602 C>A`, `s1`
`0|1` y `1|0`, en copias distintas. Llamada `c1 600 . AGC CGC,AGA` con
`s1 1/2`. El `GT` se proyecta por alelo: el alelo 1 aporta la posición 600 y el
2 la 602, cada una heterocigota. Son **dos TP**. Sin la descomposición, VE-D11
lo contaría como registro multialélico sin correspondencia, es decir un FP, pese
a contener exactamente la verdad.

**E7 — Alelo complejo.** Llamada `c1 700 . AT A`. Las longitudes de `REF` y
`ALT` diferen, así que VE-D13 no lo descompone. Se inventaría íntegro como **FP
técnico** con causa `complex-allele`. No se descarta por quedar fuera del
alcance biológico: VE-D11 lo prohíbe expresamente.

**E8 — Par en desequilibrio de ligamiento.** Verdad: dos loci separados 5 kb con
`r²` conocido, ambos recuperados en la llamada. Cada locus se clasifica por
separado como TP. El **par no recibe ninguna métrica**: VE 0.1 no evalúa LD, que
pertenece a VE.2 y sigue aplazada. Se registra «no evaluado en 0.1». Dos TP
individuales no autorizan a afirmar que la señal de LD se conserva.

### 6.4.2. Criterio de cierre por causa verificable

VE-SNP 0.1 no fija un porcentaje de aprobación. Un umbral del 100 % sobre una
referencia real es inalcanzable por propiedades ya conocidas del flujo —fusión
MNP, mapeo ambiguo, extremos de contig, bases `N`—, de modo que garantizaría un
fallo previsible seguido de una renegociación del contrato después de ver los
datos, que es justo lo que VE.0.3 existe para impedir.

En su lugar, el cierre exige que **toda discrepancia quede explicada por una
causa de esta lista cerrada y respaldada por el dato objetivo que la acredita**:

| Causa | Dato que la acredita |
|---|---|
| `mnp` | El registro llamado abarca varias posiciones y su descomposición por VE-D13 recupera la verdad |
| `unmappable` | Las lecturas que solapan el locus tienen MAPQ 0 en el BAM, de modo que el `-m` declarado en §9 las descarta |
| `reference-n` | El FASTA derivado tiene una base `N` en el contexto de alineamiento del locus |
| `contig-edge` | El locus dista de un extremo de contig menos que la longitud de lectura, y su profundidad observada está degradada respecto a la mediana |
| `zero-depth` | La profundidad observada del locus es cero en la muestra discrepante |
| `complex-allele` | El ALT llamado tiene `len(REF) ≠ len(ALT)` y VE-D13 lo deja íntegro como falso positivo técnico |

La lista es cerrada y no tiene categoría de reserva. Clasificar no retira el
locus del denominador de VE-D07: sigue contando como falso negativo, solo se
registra por qué. Y la causa no se atribuye, se acredita: cada clasificación
cita la magnitud observada en el BAM, el VCF o el FASTA que la sostiene, de modo
que un tercero pueda comprobarla sin reejecutar el ensayo.

Detienen el cierre, sin reclasificación posterior:

1. **Una discrepancia sin causa de la lista**, o con una causa que su dato no
   acredita. Se investiga como defecto del simulador o del banco externo.
2. **Un recuento por causa que el revisor de PAUSA-VE.1 no considere explicado**
   por la geometría y los parámetros del ensayo. El informe publica el desglose
   completo por causa, muestra y locus precisamente para hacer visible una
   magnitud anómala: seis fusiones MNP son coherentes con mil loci sobre
   quinientas mil bases; doscientas no lo serían aunque cada una llevara su
   etiqueta correcta.

Sensibilidad, precisión y concordancia se siguen publicando con numeradores y
denominadores completos y con el denominador de VE-D07 —todas las variantes
segregantes realizadas—, pero como medida descriptiva del flujo, no como puerta
de aprobación.

### 6.5. Cobertura, alineamiento y recuperación

El informe debe separar:

- cobertura solicitada por el simulador;
- cobertura de bases presente en FASTQ;
- cobertura alineada total;
- cobertura alineada utilizable bajo el criterio pendiente de calidad de mapeo
  y base;
- soporte de cada alelo y genotipo en el caller.

No se debe concluir que una referencia está «recompuesta» solo porque la tasa
global de mapeo sea alta. Debe comprobarse también la distribución por base o
intervalo y documentar huecos, repeticiones y alineamientos ambiguos.

En VE-SNP 0.1, por tratarse del control PE neutro, se auditan por separado todas
las bases A/C/G/T de la referencia derivada. «No enmascarada» significa aquí
únicamente «distinta de `N`»: el pipeline acepta solo A/C/G/T/N y normaliza a
mayúsculas, de modo que el enmascarado blando del FASTA original no sobrevive a
la reducción y no existe ninguna otra máscara en VE-SNP 0.1. Se auditan: la
fracción con
cobertura ≥1× y ≥10×, los intervalos de profundidad cero y la dispersión de
profundidad por muestra. La sensibilidad usa como denominador todas las
variantes segregantes realizadas desde las posiciones candidatas, fijadas antes
de ver los BAM: profundidad, MAPQ y BQ no pueden retirar ninguna del
denominador. Sí explican un fallo: una ausencia de llamada o una evidencia
técnica insuficiente es un FN, y se clasifica con una causa acreditada de
§6.4.2. El informe no puede sustituir estos valores por una cobertura media
global.

### 6.6. Haplotipos y LD

FreeBayes realiza llamada local basada en haplotipos, pero esto no garantiza que
su `GT` recupere la fase verdadera de largo alcance. Con PE150 y un inserto
acotado, dos variantes lejanas pueden no estar conectadas por ninguna lectura.

La primera evaluación VE-LD puede comprobar la conservación de la señal de LD
a partir de genotipos llamados sin fase, siempre que se nombre el estadístico y
se compare con el mismo estadístico calculado sobre los mismos 20 individuos de
verdad. Esto demuestra que el callset conserva esa señal; no demuestra
recuperación exacta de haplotipos.

Si se exige recuperar haplotipos o fase por bloque, debe aprobarse un faseador,
su versión, entradas y límite de interpretación. Un faseador basado en lecturas
solo podrá unir variantes conectadas por fragmentos. Un faseador poblacional
introduce inferencia estadística y, con 20 individuos y 500 kb, no debe
considerarse un oráculo independiente sin una evaluación específica.

## 7. Arquitectura técnica mínima propuesta

El banco externo queda fuera del runtime Rust del simulador. Se implementa como
una herramienta Python independiente gestionada con `uv`: invoca procesos
locales con argumentos separados, lee BAM/VCF con dependencias versionadas y
publica su propio expediente. No usa contenedores ni incorpora minimap2,
samtools o FreeBayes al simulador. Sus componentes son:

| Componente | Responsabilidad única | Entrada | Salida o error observable |
|---|---|---|---|
| Admisor del conjunto | Resolver manifiestos, hashes, muestra y referencia antes de ejecutar herramientas | Publicación del simulador | Inventario admitido o lista de incompatibilidades |
| Preparador de referencia | Verificar FASTA, crear índices externos y registrar su versión | FASTA derivado | Índices identificados o fallo de referencia |
| Adaptador minimap2 | Alinear cada par R1/R2 y fijar el grupo de lectura de la muestra | FASTA y FASTQ | SAM por flujo y diagnóstico de proceso |
| Adaptador samtools | Ordenar, validar e indexar cada alineamiento sin perder lecturas | SAM | BAM/BAI por muestra e informes |
| Adaptador FreeBayes | Realizar una llamada conjunta de las 20 muestras | FASTA y BAM | VCF bruto y diagnóstico de proceso |
| Normalizador/comparador | Crear vistas comparables sin alterar la verdad, descomponer los MNP conforme a VE-D13 y clasificar discordancias | VCF de verdad y llamado | Tabla de concordancia con procedencia por SNP y métricas acordadas |
| Evaluador LD/fase | Calcular solo el estimador aprobado y declarar si usa fase observada o inferida | Verdad y llamada comparable | Tabla por par/bloque y limitaciones |
| Ensamblador de evidencia | Registrar versiones, comandos, hashes, recursos, estados y resultados | Todas las etapas | Expediente reproducible completo o fallido |

Flujo propuesto:

```text
publicación del simulador
  ├─ referencia FASTA ──> índice minimap2 ───────────────┐
  ├─ FASTQ R1/R2 × 20 ─> minimap2 ─> BAM/BAI × 20 ─────┼─> FreeBayes
  └─ VCF de verdad ──────────────────────────────────────┘       │
                                                                  v
                                      vistas VCF normalizadas y comparación
                                                                  │
                                     métricas SNP ─────── métricas LD/fase
                                                                  │
                                                   expediente reproducible
```

### 7.1. Invariantes de frontera

- Toda etapa usa el mismo FASTA exacto; no se permite mapear contra el genoma
  fuente y comparar contra coordenadas de la referencia derivada. El nombre no
  lo impide: la reducción conserva los identificadores de contig del original
  con longitudes distintas, de modo que ambos ficheros declaran `SL4.0ch01`. Lo
  impide el `checksum` del artefacto publicado por la ejecución A, que el
  admisor resuelve antes de indexar, y la comparación de longitudes `@SQ` de
  VE.1.4.
- Las 20 columnas del VCF de verdad, los `SM` de BAM y las 20 muestras llamadas
  deben formar el mismo conjunto biyectivo.
- El ensayo declara exactamente una biblioteca por individuo. El esquema admite
  varias (`libraries` es un array sin máximo), pero F4 publica un solo par de
  FASTQ por individuo, de modo que con dos bibliotecas las lecturas de dos
  grupos distintos se fundirían en el mismo fichero y la tabla biyectiva de §6.3
  dejaría de existir.
- Cada BAM contiene un único grupo de lectura: `RG.ID` coincide con el
  `read_group_id` publicado por F4, `RG.SM` con el individuo y `RG.LB` con su
  biblioteca. Cada alineamiento porta `RG:Z` resoluble en esa cabecera.
- La lista de BAM dada a FreeBayes es exactamente la tabla canónica, una vez por
  individuo y en orden estable; no admite muestras añadidas, omitidas ni BAM con
  varios `SM`.
- El nombre o la ruta de un fichero no sustituyen su SHA-256 y su identidad en
  manifiesto.
- Un proceso externo con código distinto de cero, salida truncada o versión no
  aprobada falla la etapa y no publica un resultado «completo».
- La normalización de VCF conserva una copia identificada del original y
  registra cada transformación.
- Las métricas de fase y de LD indican si proceden de verdad faseada, llamada
  sin fase o fase inferida.
- Un resumen global nunca oculta una discordancia por muestra, locus o alelo.

## 8. Criterios de aceptación propuestos

VE-D08 ha sustituido el umbral numérico por el criterio de causa verificable de
§6.4.2, de modo que no queda ningún porcentaje pendiente de acordar. Las
condiciones estructurales siguientes lo acompañan:

1. **Admisión:** hashes, referencia y 20 identidades son coherentes; todos los
   artefactos requeridos existen y los manifiestos están completos.
2. **Alineamiento:** cada par FASTQ se consume una sola vez; el total de
   lecturas de entrada se explica como alineado o no alineado; cada BAM es
   válido, coordinado, indexable y conserva la muestra correcta. Cada par no
   mapeado se enumera por muestra y se investiga antes de aprobar el piloto
   neutro.
3. **Cobertura:** se informa la distribución por muestra y región, no solo una
   media global: fracciones ≥1×/≥10×, mediana, cuantiles, dispersión y todos
   los intervalos de cobertura cero. Ningún hueco queda sin clasificar. El
   denominador son las bases A/C/G/T de la referencia derivada, y se publica
   además desglosado por contig, porque la geometría de §5.4 concentra los
   efectos de borde en los contigs cortos. La profundidad se calcula solo con
   alineaciones primarias, con `samtools depth -a` para que las bases sin
   cobertura cuenten como cero en lugar de desaparecer; secundarias y
   suplementarias se informan aparte y no entran en la distribución.
4. **Llamada:** FreeBayes finaliza sobre los 20 BAM y produce un VCF legible,
   normalizable y asociado a la referencia correcta.
5. **SNP:** cada verdad, llamada y discordancia puede trazarse a locus, alelo,
   genotipo y muestra. Sensibilidad, precisión y concordancia no se publican sin
   sus numeradores y denominadores. El cierre lo decide la clasificación
   acreditada de cada discrepancia, no el valor de esas tres medidas.
6. **LD/fase:** el informe separa la conservación estadística de LD de la
   recuperación de fase. No declara haplotipos exactos si no existe evidencia
   que conecte los loci.
7. **Reproducibilidad:** A→B se repite con las mismas configuraciones y
   semillas, y cada artefacto se compara en el nivel que le asigna §8.1. El
   flujo externo se repite sobre esas mismas entradas; tiempos, logs, índices y
   metadatos incidentales se registran por separado.

### 8.1. Nivel de reproducibilidad exigido a cada artefacto

DTC-038 declara por separado cuatro niveles —biología, decisiones y lecturas
técnicas, equivalencia semántica e identidad binaria—, mantiene comparadores
distintos para cada uno y prohíbe deducir una garantía fuerte desde otra más
débil. V1 no garantiza identidad binaria por defecto para ninguna clase: solo
puede declararse por formato y condiciones tras una prueba expresa. REQ-CONF-115
añade que una diferencia limitada a compresión, ordenación equivalente o
metadatos no biológicos no invalida por sí misma los niveles anteriores.

Por eso VE-SNP 0.1 no exige identidad binaria en bloque, que trataría como
fallo de aceptación justo lo que esa regla excluye:

| Artefacto | Nivel exigido | Por qué |
|---|---|---|
| FASTA reducido y mapa de coordenadas | Identidad binaria | Texto sin comprimir y determinista dada la semilla y la configuración |
| `effective-config.json` y manifiestos | Identidad binaria | Sus esquemas no contienen marca de tiempo, duración ni host, de modo que nada incidental los separa; pero copian literalmente la ruta pedida en `configuration_input.original_path`, lo que obliga a la segunda condición de abajo |
| FASTQ `.fastq.gz` | Identidad binaria, bajo las condiciones ya declaradas para el perfil v1: mismo simulador, mismo `Cargo.lock` y misma plataforma admitida | Es la declaración por formato y condiciones que DTC-038 admite. Fuera de ellas, el nivel baja a equivalencia semántica del contenido descomprimido |
| VCF de verdad | Equivalencia semántica | Ningún documento del proyecto declara identidad binaria para su BGZF, y REQ-CONF-115 impide convertir una diferencia de compresión en un fallo |
| Índices `.tbi`, `.csi` y `.fai` | Equivalencia semántica | Son derivados reconstruibles; se registran, pero no deciden la aceptación |

Ningún nivel es exigible si no se cumplen las tres condiciones siguientes, que
el expediente comprueba antes de comparar nada:

1. **`seed` declarada en el TOML de A y de B.** Si falta, `resolve_root_seed` la
   toma del CSPRNG del sistema y registra `SeedOrigin::SystemCsprng`; entonces
   las dos ejecuciones parten de semillas distintas y la repetición no significa
   nada.
2. **La misma cadena de ruta de configuración en ambas invocaciones.** La ruta
   se copia tal cual se pidió, así que invocar una vez con ruta relativa y otra
   con absoluta cambia la configuración efectiva y, con ella, el manifiesto.
3. **El mismo binario, el mismo `Cargo.lock` y la misma plataforma.** Los bytes
   comprimidos dependen de la versión de la biblioteca que los produce.

VE.3.1 es la prueba expresa que DTC-038 exige para poder declarar identidad
binaria en los formatos de la tabla, y el expediente la declara con esas
condiciones. Hoy no existe: la prueba de repetición del repositorio compara el
JSON del manifiesto ya interpretado, no sus bytes, de modo que acredita
equivalencia semántica y no identidad binaria. Cada manifiesto registra además,
como pide DTC-038, si el nivel quedó `garantizado`, `no garantizado` o
`no aplicable`.

### 8.2. Igualdad de BAM y VCF producidos fuera del simulador

Los artefactos del flujo externo no los produce el simulador, así que §8.1 no
los alcanza: caen en el nivel de equivalencia semántica de DTC-038, que ante
herramientas externas exige además sus mismas versiones, parámetros y
condiciones verificadas. VE-D05 ya las congela.

Se midió sobre las herramientas congeladas, repitiendo dos veces el mismo
comando con las mismas entradas. Dos VCF de FreeBayes resultaron **idénticos
byte a byte**. Dos BAM resultaron idénticos en todos sus registros de
alineamiento y en `@HD`, `@SQ` y `@RG`, y difirieron solo en las líneas `@PG`,
cuyo campo `CL:` reproduce la orden ejecutada y con ella el nombre del fichero
de salida. Los índices `.bai` fueron idénticos.

Por tanto la comparación se define en dos pasos, y no al revés:

1. **Primero byte a byte.** Es lo normal cuando la repetición ocurre en la misma
   jornada y con las mismas rutas, y es la comprobación más fuerte disponible.
2. **Solo si falla**, se repite excluyendo las líneas incidentales de la tabla
   siguiente, y el informe publica cuáles causaron la diferencia. La exclusión
   nunca se aplica de entrada: así no puede convertirse en una manta que oculte
   una diferencia real.

| Formato | Decide la igualdad | Incidental, se registra pero no decide |
|---|---|---|
| BAM | Los registros de alineamiento, en el mismo orden, y las cabeceras `@HD`, `@SQ` y `@RG` | Las líneas `@PG`: su `CL:` incorpora la ruta del fichero producido |
| VCF externo | El cuerpo de registros y las cabeceras `##INFO`, `##FORMAT`, `##contig` y la fila `#CHROM` | `##fileDate`, que cambia si la repetición cruza la medianoche; `##commandline` y `##reference`, que incorporan rutas; `##source`, cuya versión ya congela VE-D05 |
| Índice `.bai` | Identidad binaria | — |

Una diferencia en cualquier columna «decide» detiene el cierre. Una diferencia
limitada a la columna incidental no lo detiene, conforme a REQ-CONF-115, pero se
publica con la línea concreta que la produjo: nunca basta con afirmar que las
salidas son «lógicamente iguales».

## 9. Decisiones y dependencias de ejecución

| ID | Decisión o dependencia | Estado acordado | Bloquea |
|---|---|---|---|
| VE-D01 | Propiedad y coordenada | **Aprobado:** banco externo DTF-024 con coordenada `VE`, no F5 ni F9.5 | No bloquea |
| VE-D02 | PE sin errores de secuenciación | **Superado en F4:** `sequencing_errors="disabled"` produce cero sustituciones PE y marca las calidades no probabilísticas; VE fija además quimeras, contaminación e *index hopping* a cero | No bloquea VE-SNP |
| VE-D03 | Contenido exacto de VE-SNP | **Aprobado:** solo SNP bialélicos, sin multialelia ni *indels*, en VE-SNP 0.1 | Configuración y comparador |
| VE-D04 | Densidad y selección de loci | **Aprobado:** 1.000 posiciones candidatas A/C/G/T elegibles, no enmascaradas ni ambiguas, con frecuencia ALT objetivo uniforme en [0,1; 0,4). Los produce el generador F2.2 ya existente, sin código nuevo: `[population_snvs] mode="generate"`, `count_mode="exact"`, `locus_count=1000`, `spatial_mode="uniform"`, y por población `frequency_model="uniform"` con `min_alt_af=0.1` y `max_alt_af=0.4`. El VCF conserva solo las variantes segregantes realizadas: aquellas cuyos cuarenta haplotipos no resultaron todos REF ni todos ALT. Se registran candidatos, frecuencias objetivo, genotipos realizados y recuento final, sin reintentar ni seleccionar para forzar 1.000 SNP. La pérdida esperada por monomorfismo es de 1,08 loci y no supera 5 con probabilidad 0,999, de modo que un recuento en [995, 1000] es normal y uno sensiblemente menor acusa un fallo del generador, no azar. No se impone separación mínima entre loci: la fusión MNP que provoca la proximidad se resuelve en el comparador mediante VE-D13, no alterando la verdad | Generación y criterio de exactitud |
| VE-D05 | Versiones y argumentos | **Aprobado:** congelar minimap2 2.26-r1175, samtools/bcftools 1.21 locales y versión+hash de FreeBayes v1.3.10 antes de ejecutar. HTSlib 1.24 no es requisito del banco externo, pero su cierre no acredita DTC-066: esa puerta exige bcftools y tabix 1.24 y sigue pendiente para el simulador. El expediente lo declara expresamente para que el cierre de VE no se cite como prueba de interoperabilidad que no ha demostrado. Si la herramienta Python del banco embebe su propio HTSlib, su versión se congela también aquí | VE.1.1, comandos reproducibles |
| VE-D06 | Normalización y emparejamiento VCF | **Aprobado para VE-SNP 0.1:** comparar SNP bialélicos por CHROM, posición, REF, ALT y dosis por muestra; `0|1`, `1|0` y `0/1` son equivalentes. Fase se evalúa aparte. Ignorar QUAL/FILTER/INFO/GPS_CASE; normalización compleja se aplaza a multialelia/*indels*. Los registros MNP se descomponen antes de comparar conforme a VE-D13 y sus SNP derivados usan estas mismas reglas | Métricas SNP/indel |
| VE-D07 | Regiones llamables | **Aprobado para VE-SNP 0.1:** todas las variantes segregantes realizadas desde las posiciones candidatas constituyen el denominador de sensibilidad. Ninguna se retira del denominador por profundidad, MAPQ o BQ; esas magnitudes explican un FN, no lo eliminan. Los filtros de entrada del llamador sí existen y se declaran en §9: no se prohíben, se hacen visibles y su efecto se acredita por VE-D08. Una ausencia de llamada es FN y se clasifica por causa. `N`, máscara y repetición se examinan en el ensayo diferido de §13 | Denominadores y criterio de cierre |
| VE-D08 | Criterio de cierre | **Aprobado:** no se exige un porcentaje. El cierre exige que toda discrepancia lleve una causa de la lista cerrada de §6.4.2 y el dato observado que la acredita, y que el informe publique el desglose por causa, muestra y locus para que una magnitud anómala sea visible en PAUSA-VE.1. Sensibilidad, precisión y concordancia se publican siempre con numeradores y denominadores, contados en SNP tras aplicar VE-D13 y no en registros VCF, pero como medida descriptiva y no como puerta. Definición completa en §6.4.2 | Cierre científico |
| VE-D09 | Métrica LD | **Aplazado:** el ensayo de LD poblacional se realizará después de elegir un faseador y su evaluación propia | VE-LD |
| VE-D10 | Recuperación de fase | **Aprobado para VE 0.1:** no usar faseador y no evaluar fase. No es una rama condicional: FreeBayes v1.3.10 declara `##phasing=none` en su cabecera y emite `GT` sin fase, comprobado sobre el ejecutable congelado. El informe registra «fase no evaluable» y no planifica una tabla que no tendría contenido. Fase por bloques queda para después | Afirmación sobre haplotipos |
| VE-D11 | Multialelia e *indels* | **Aprobado:** excluirlos deliberadamente de la verdad VE-SNP 0.1 y crear después casos propios. Todo ALT adicional, *indel* o registro multialélico llamado se inventaría y cuenta como FP técnico, salvo que la descomposición de VE-D13 lo resuelva en SNP con verdad correspondiente; nunca se descarta silenciosamente | Comparación completa |
| VE-D12 | Repeticiones | **Aprobado:** la repetición no usa un criterio propio, sino los cuatro niveles de DTC-038, y exige uno por artefacto en lugar de identidad binaria en bloque. §8.1 fija el nivel de cada uno y las tres condiciones sin las cuales ninguno es exigible. Donde VE.3.1 exija identidad binaria, actúa como la prueba expresa que DTC-038 requiere para poder declararla, y el expediente la declara por formato y condiciones. Réplicas con otras semillas se aplazan | Generalización del resultado |
| VE-D13 | Registros multinucleotídicos | **Aprobado:** no se altera la verdad para evitar la fusión. El comparador descompone cada ALT con `len(REF) = len(ALT)` en un SNP por posición discrepante, proyecta el `GT` por alelo y compara cada SNP derivado con las reglas de VE-D06. Cubre *cis* (MNP bialélico) y *trans* (MNP multialélico `1/2`). Un ALT con longitudes distintas no se descompone: es FP técnico íntegro. Se publica la procedencia de cada SNP y el recuento de registros fusionados. Definición completa en §6.4.1 | Comparador VE.1.6 y métricas VE-D08 |
| VE-D14 | Geometría de la referencia derivada | **Aprobado:** el genoma fuente y la geometría de la reducción son parámetros de cada ensayo, igual que la cohorte o la cobertura, y el banco externo no fija ninguno como constante. Un ensayo declara `region_size_bases` y elige `target_bases` o `reduction_percent`, que el esquema hace mutuamente excluyentes. El primer ensayo usa el tomate `SL4.00`, `region_size_bases = 10000` y `target_bases = 500000`, porque solo `target_bases` da la longitud exacta que exige VE.1.2. La geometría resultante no se predice en el contrato: VE.1.2 publica la realizada y de ella dependen las causas `contig-edge` y `zero-depth` de §6.4.2. Detalle en §5.4 | VE.1.2, cobertura y causas de §6.4.2 |
| VE-D16 | Igualdad de BAM y VCF externos | **Aprobado:** «igualdad lógica» deja de ser un término sin definir. La comparación se intenta primero byte a byte; solo si falla se repite excluyendo las líneas que §8.2 declara incidentales, y el informe publica cuáles causaron la diferencia. Un BAM es igual si coinciden sus registros de alineamiento y sus cabeceras `@HD`, `@SQ` y `@RG`; un VCF externo, si coinciden su cuerpo y sus `##INFO`, `##FORMAT`, `##contig` y la fila `#CHROM`. Se midió sobre las herramientas congeladas: dos VCF de FreeBayes salen idénticos byte a byte, y dos BAM solo difieren en `@PG`. Detalle en §8.2 | VE.3.1 |
| VE-D15 | Enlace A→B e identidad de la referencia | **Aprobado:** el FASTA derivado conserva los nombres de contig del genoma fuente con longitudes distintas, de modo que el nombre nunca acredita identidad. Dentro del simulador la confusión ya es imposible: la admisión rechaza un FASTA cuyas secuencias no coincidan en identificador, orden y longitud con lo declarado. La comprobación que falta es la del banco externo, que no valida nada: su admisor resuelve el FASTA por el `checksum` del artefacto publicado en el manifiesto de A —no por ruta ni por nombre— y rechaza antes de construir el índice si no coincide. VE.1.4 comprueba además que las longitudes de las cabeceras `@SQ` de cada BAM coinciden con las del FASTA derivado, que es lo que delata una confusión que los nombres no delatarían. La configuración de B se genera desde el mapa de A y no se transcribe a mano | Admisor del banco, VE.1.3 y VE.1.4 |

El ensayo se declara con `schema_version = "2"`. Es la única ruta que admite
`population_snvs`, `initial_populations` y `technical_event`; la v3 que
introduce F5 solo acepta `f5_scenario` y no puede expresar VE-SNP. Las dos son
mutuamente excluyentes.

Cada tubería se ejecuta con `set -o pipefail`. Sin él, el estado de salida es
el del último proceso: si minimap2 aborta a mitad, `samtools sort` cierra
ordenadamente lo recibido y devuelve cero, produciendo un BAM válido e
indexable al que solo le faltan lecturas. Además, el recuento de lecturas de
entrada que exige §8.2 se toma del manifiesto de F4 y no del propio BAM, que no
puede acusar su propia falta.

Los argumentos literales de VE-SNP 0.1 son `minimap2 -t 4 -ax sr -R
<RG-canonico> referencia.fa R1.fastq.gz R2.fastq.gz | samtools sort -@ 4 -o
<sample>.bam -`, seguido de `samtools index -@ 4 <sample>.bam` por muestra, y
`freebayes -f referencia.fa -p 2
--genotype-qualities --min-alternate-count 1 --min-alternate-fraction 0
--haplotype-length 3 --min-mapping-quality 1 --min-base-quality 0` sobre
los veinte BAM canónicos. Se conserva el VCF bruto y no se aplica filtrado
posterior que oculte ALT; las alineaciones secundarias/suplementarias se
informan, mientras que cobertura y recuento de pares usan solo primarias.

Los tres últimos argumentos se declaran aunque coincidan con los valores
predeterminados de FreeBayes v1.3.10, porque gobiernan las métricas y no pueden
quedar implícitos en un predeterminado que una versión futura podría cambiar.
`--haplotype-length 3` decide la fusión MNP que VE-D13 descompone.
`--min-mapping-quality 1` descarta las alineaciones de MAPQ 0, que minimap2
asigna a las lecturas que no puede situar de forma única: es un filtro real y
activo, y de él depende la causa `unmappable` de §6.4.2.
`--min-base-quality 0` se declara para dejar constancia de que no hay filtro por
calidad de base.

Se comprobó sobre el ejecutable local que el resto de filtros de entrada de este
comando es inerte y no necesita declararse: `-R`, `-Y`, `-3` y `--min-coverage`
valen cero; `-U` y `-$` no tienen cota; `-z` vale 1,0, de modo que `-Q 10` no
excluye ninguna lectura; y `-4` solo excluiría alineaciones marcadas como
duplicadas, que no existen porque VE-SNP no ejecuta `markdup` y ninguna lectura
lleva el indicador `0x400`. VE.1.1 vuelve a comprobar esta enumeración contra la
versión congelada antes de ejecutar.

## 10. Plan ejecutable

Los identificadores siguientes usan la coordenada aprobada `VE`. VE.0 está
cerrado y VE.1.1 puede empezar. Ninguna ejecución de herramientas se considera
superada hasta que VE.1.1 haya fijado y comprobado sus identidades.

### VE.0 — Cierre del contrato del banco externo

- ✅ **VE.0.1. Aprobar alcance, vocabulario y afirmación científica.**
  - Depende de: ninguno.
  - Procede de: §§2–3 y VE-D01.
  - Cambia: documento propietario de DTF-024; no modifica la especificación del
    simulador.
  - Se comprueba con: acta que distingue alineamiento, consenso, recuperación
    de variantes, ensamblaje y fase.

- ✅ **VE.0.2. Resolver el escenario técnico PE neutro.**
  - Depende de: VE.0.1.
  - Procede de: §5.3 y VE-D02.
  - Cambia: F4 materializa el modo PE `disabled`; VE solo lo configura junto a
    los demás artefactos técnicos neutros.
  - Se comprueba con: una configuración admitida cuyo manifiesto demuestre el
    modo `disabled`, cero sustituciones realizadas, calidades no probabilísticas
    señaladas, y cero contaminación, quimeras e *index hopping*, sin usar F5.

- ✅ **VE.0.3. Fijar contrato de comparación y umbrales.** *(Reabierto el 4 de
  septiembre de 2026 por la revisión adversarial de §1.2 y cerrado de nuevo el 5
  de septiembre.)*
  - Depende de: VE.0.1.
  - Procede de: §§6.4–6.6 y VE-D03–VE-D13.
  - Cambia: tabla normativa de métricas, denominadores, normalización, regiones
    llamables y tratamiento de fase/multialelia/*indels*.
  - Se comprueba con: ejemplos manuales de TP/FN/FP, genotipo discordante,
    registro multialélico, región no llamable y par LD que no admitan dos
    clasificaciones distintas.
  - Motivo de la reapertura: su criterio de cierre exige ejemplos que no
    admitan dos clasificaciones. Un registro MNP las admitía —concordancia o
    falso positivo técnico— según se leyera VE-D06 o VE-D11, y VE-D07 y VE-D08
    se contradecían sobre el denominador.
  - Cerrado con: los siete hallazgos graves resueltos por VE-D13 a VE-D16, la
    pasada de redacción sobre los medios y menores, y los ocho ejemplos
    manuales de §6.4.3, reescritos contra el contrato corregido y con una sola
    clasificación admisible cada uno.

> **PAUSA-VE.0 cerrada de nuevo el 5 de septiembre de 2026.** El contrato de
> comparación queda fijado con VE-D13 a VE-D16 y los ejemplos de §6.4.3. El
> siguiente paso es VE.1.1: congelar las identidades del entorno local antes de
> preparar o ejecutar el ensayo.

### VE.1 — Conjunto pequeño, alineamiento y SNP independientes

- ☐ **VE.1.1. Congelar el entorno externo y sus comprobaciones.**
  - Depende de: PAUSA-VE.0 y VE-D05.
  - Procede de: §§5.2, 6.1 y 7.
  - Cambia: herramienta Python VE gestionada con `uv`, manifiesto de
    dependencias y receta de ejecución; no el runtime del simulador.
  - Se comprueba con: versiones y hashes de ejecutables y dependencias Python,
    ayuda de cada comando y un fixture diminuto válido e inválido. Las cinco
    herramientas están instaladas —FreeBayes en `/usr/local/bin/freebayes`,
    v1.3.10—, de modo que el paso congela identidades y no instala ni descarga
    nada; §13 lo prohíbe. Tampoco usa contenedor.
  - Vuelve a comprobar la enumeración de filtros inertes de §9 contra la versión
    congelada: un valor predeterminado que cambie de versión alteraría las
    métricas sin avisar.

- ☐ **VE.1.2. Derivar y congelar la referencia de tomate de 500.000 bases.**
  - Depende de: VE.1.1.
  - Procede de: F1.6, ACE-112/113 y §4.
  - Cambia: configuración de reducción del ensayo y su expediente externo.
  - Se comprueba con: ejecución A del simulador, reconstrucción independiente
    FASTA fuente+mapa, longitud exacta 500.000, `samtools faidx`, SHA-256 y
    ausencia de mezcla entre coordenadas originales y derivadas.
  - Publica la geometría realizada conforme a VE-D14, leída del FASTA y del mapa
    y no supuesta desde la configuración: número de contigs derivados, longitud
    de cada uno, bloques que lo componen y sus intervalos de origen. De ella
    dependen las causas `contig-edge` y `zero-depth` de §6.4.2, de modo que un
    ensayo con otra geometría cambia el reparto de causas sin cambiar el
    contrato.

- ☐ **VE.1.3. Generar VE-SNP con 20 individuos y 30× PE.**
  - Depende de: VE.0.2–VE.0.3 y VE.1.2.
  - Procede de: §§4, 6.1–6.2 y VE-D03/04.
  - Cambia: ejecución B, cuya referencia principal suministrada es exactamente
    el FASTA publicado por VE.1.2; genera SNP, VCF y FASTQ contra coordenadas
    derivadas, no contra el tomate original ni contra otra referencia.
  - Se comprueba con: manifiestos y hashes; 20 pares R1/R2; las 1.000 posiciones
    candidatas y sus frecuencias objetivo; genotipos realizados y número de SNP
    segregantes finalmente publicado en el VCF; cobertura solicitada 30×;
    enlace A→B, identidad y SHA-256 de referencia iguales entre B, VCF y FASTQ;
    artefactos técnicos realizados a cero según el contrato aprobado; validadores
    FASTQ/VCF externos. No se reintenta ni filtra para alcanzar un número de SNP.

- ☐ **VE.1.4. Alinear e indexar las 20 muestras con minimap2 y samtools.**
  - Depende de: VE.1.1 y VE.1.3.
  - Procede de: §§2, 6.3/6.5 y arquitectura §7.
  - Cambia: índices de referencia, BAM/BAI e informes externos por muestra.
  - Se comprueba con: `samtools quickcheck`, cabeceras `@SQ/@RG`, `SM` exacto,
    tabla biyectiva individuo→FASTQ→RG→BAM, `RG:Z` en cada alineamiento, conteo
    de lecturas explicable, pares no mapeados enumerados e informes de
    profundidad por base/intervalo que localicen huecos y midan homogeneidad.
    Los argumentos concretos se fijan en VE-D05; no se improvisan al ejecutar.
  - Comprueba también, conforme a VE-D15, que cada `@SQ` del BAM coincide con el
    FASTA derivado en identificador **y longitud**. Los identificadores son los
    mismos que los del genoma fuente, así que solo la longitud distingue un BAM
    construido contra la referencia equivocada.

- ☐ **VE.1.5. Llamar conjuntamente VE-SNP con FreeBayes.**
  - Depende de: VE.1.4.
  - Procede de: §§2 y 6.3–6.4.
  - Cambia: VCF bruto y log de FreeBayes.
  - Se comprueba con: proceso exitoso, 20 muestras correctas, referencia y
    contigs compatibles, VCF formalmente legible y ausencia de BAM omitidos o
    duplicados.

- ☐ **VE.1.6. Normalizar y comparar VE-SNP con el VCF de verdad.**
  - Depende de: VE.0.3 y VE.1.5.
  - Procede de: §§6.4–6.5, VE-D06–VE-D08 y VE-D13.
  - Cambia: vistas VCF derivadas, descomponedor de MNP, tabla de discordancias
    con procedencia por SNP e informe de métricas.
  - Se comprueba con: oráculos manuales, numeradores/denominadores, listado
    exhaustivo por muestra/locus y cumplimiento de los umbrales aprobados. No
    compara `GPS_CASE`, `QUAL` o campos del caller salvo decisión expresa.
  - El descomponedor de VE-D13 se valida con fixtures sintéticos deterministas,
    no con los casos que el ensayo produzca por azar: MNP bialélico en *cis*,
    MNP multialélico `1/2` en *trans*, base de anclaje intermedia sin cambio,
    posición cambiada sin verdad correspondiente, verdad ausente dentro del
    tramo fusionado, homocigoto `1/1`, tramo de longitud cuatro con tres
    cambios, ALT con `len(REF) ≠ len(ALT)` que debe quedar íntegro como falso
    positivo, y SNP simple que debe atravesar el módulo sin alterarse. Cada
    fixture fija de antemano la clasificación esperada de cada posición.
  - El cierre aplica §6.4.2: toda discrepancia lleva una causa de la lista
    cerrada y el dato observado que la acredita, y el informe publica el
    desglose por causa, muestra y locus. Una discrepancia sin causa, o con una
    causa que su dato no sostiene, detiene el cierre.

> **PAUSA-VE.1.** Revisar los BAM, el VCF llamado, todas las discordancias y los
> denominadores. No pasar a LD si el caso independiente no explica cada fallo.

### VE.2 — Ampliación aplazada: haplotipos, LD y fase

> **Fuera de VE 0.1.** Esta ampliación solo se abrirá después de VE.3, cuando
> exista una salida F5 aprobada y se haya elegido y evaluado un faseador. No es
> una dependencia de la repetición ni del cierre de VE-SNP.

- ☐ **VE.2.1. Admitir una salida VE-LD de F5 sin acoplarse a su implementación.**
  - Depende de: PAUSA-VE.1 y de una salida F5 aprobada que publique verdad,
    haplotipos, fase, genotipos y LD conforme a sus contratos.
  - Procede de: §§6.2/6.6 y VE-D09–VE-D10.
  - Cambia: adaptador de entrada del banco externo; no F5.
  - Se comprueba con: lector independiente que reconstruya desde la verdad F5
    los haplotipos y el estadístico objetivo antes de ver el VCF llamado.

- ☐ **VE.2.2. Generar, alinear y llamar VE-LD bajo el mismo contrato técnico.**
  - Depende de: VE.2.1 y VE.1.3–VE.1.5.
  - Procede de: §§4 y 6.2/6.6.
  - Cambia: segundo conjunto de resultados, BAM y VCF; mantiene separadas las
    identidades y semillas de VE-SNP.
  - Se comprueba con: las mismas puertas estructurales de VE.1 y comprobación de
    que los bloques/variantes de verdad esperados están presentes antes de
    evaluar la llamada.

- ☐ **VE.2.3. Evaluar genotipos y conservación de LD.**
  - Depende de: VE.0.3, VE.2.2 y VE-D09.
  - Procede de: §§6.4 y 6.6.
  - Cambia: tablas por locus, par, distancia o bloque según el contrato.
  - Se comprueba con: mismo estimador y mismos individuos en verdad y llamada,
    intervalos/pares idénticos, diferencias y tolerancias explícitas y casos
    degenerados documentados.

- ☐ **VE.2.4. Evaluar fase solo si VE-D10 la exige.**
  - Depende de: VE.2.2 y faseador aprobado.
  - Procede de: §6.6 y VE-D10.
  - Cambia: salida de fase inferida y comparación limitada a bloques
    identificables.
  - Se comprueba con: switch error y/o la métrica aprobada, límites de bloque,
    loci no faseables separados y prueba negativa con loci no conectados. No se
    atribuye fase exacta a FreeBayes por el mero hecho de recuperar genotipos.

> **PAUSA-VE.2.** Si esta ampliación se abre, decidir por separado si la
> evidencia demuestra recuperación de genotipos, conservación de señal LD o
> fase por bloque. No fusionar esas conclusiones.

### VE.3 — Repetición y cierre

- ☐ **VE.3.1. Repetir los conjuntos aprobados y conservar recursos observados.**
  - Depende de: PAUSA-VE.1 y VE-D12.
  - Procede de: §§6.3 y 8.
  - Cambia: repetición exacta A→B de VE-SNP 0.1 y repetición del flujo externo,
    con tabla de recursos.
  - Se comprueba con: el nivel que §8.1 asigna a cada artefacto, comprobando
    antes sus tres condiciones —`seed` declarada, misma cadena de ruta de
    configuración y mismo binario, `Cargo.lock` y plataforma—; la comparación de
    BAM y VCF externos que define §8.2; `/usr/bin/time` u observador
    equivalente y tamaños por etapa. Logs, tiempos y metadatos incidentales se
    comparan y registran por separado.
  - Una diferencia se clasifica antes de juzgarla: si afecta solo a compresión,
    ordenación equivalente o metadatos no biológicos, REQ-CONF-115 impide
    tratarla como fallo de los niveles anteriores, y se registra sin detener el
    cierre. Detiene el cierre una diferencia de contenido en cualquier nivel
    exigido, y también una identidad binaria que se pierda en un formato donde
    §8.1 la exige.

- ☐ **VE.3.2. Ensamblar el expediente reproducible y revisar el cierre.**
  - Depende de: VE.3.1.
  - Procede de: §§6.3, 7.1 y 8.
  - Cambia: manifiesto e informe final del banco externo.
  - Se comprueba con: todos los enlaces y hashes resolubles, comandos y versiones
    completos, cero resultado parcial marcado como completo y revisión
    independiente de las afirmaciones científicas.

> **PAUSA-VE.3 / CIERRE VE 0.1.** Aprobar interoperabilidad FASTQ/BAM y
> recuperación SNP. LD y fase pertenecen exclusivamente a VE.2 cuando se abra.

## 11. Recursos y tiempo de planificación

El único volumen calculable antes de ejecutar es el de 300 millones de bases y
aproximadamente un millón de parejas. El tamaño comprimido, el BAM, la memoria
de FreeBayes y el tiempo dependen de la complejidad de la referencia, densidad
de variantes, parámetros y versiones aún abiertas.

Para no presentar una estimación como resultado, VE.1 debe empezar con un
presupuesto operativo provisional y sustituirlo por medidas observadas. Como
reserva de agenda, no como requisito ni predicción, se propone una estación con
4 CPU, 8 GiB de RAM, 10 GiB libres y una ventana de 2 horas para el primer
recorrido.

Esa reserva debe cubrir también la ejecución A, que es la etapa más exigente en
memoria aunque produzca el artefacto más pequeño: la reducción recorre el FASTA
de tomate completo, de 759 MB, para extraer 500.000 bases. El flujo externo, en
cambio, es modesto: un millón de pares, un índice de minimap2 sobre 500 kb y
BAM que suman unos cientos de megabytes. Conviene medir A por separado y no
suponer que el presupuesto del ensayo lo cubre por ser su artefacto pequeño.
VE.1.4–VE.1.6 registrarán uso real por etapa. Si el preflight o las
medidas superan esa reserva, el ensayo se detiene y se revisa el presupuesto; no
reduce cobertura, individuos ni evidencia automáticamente.

## 12. Matriz de evidencia y riesgos

| Hito | Resultado | Evidencia principal | Riesgo protegido |
|---|---|---|---|
| VE.0 | Contrato ejecutable | Decisiones aprobadas y casos manuales inequívocos | Medir algo distinto de lo que se afirma |
| VE.1 | Recuperación externa de SNP | FASTQ→BAM→FreeBayes→comparación por muestra/locus | Mezcla de referencia, muestra, representación VCF o denominador |
| VE.2 (aplazada) | Evaluación separada de LD y fase | Verdad F5 independiente y mismo estimador | Confundir genotipo, LD, haplotipo y fase inferida |
| VE.3 | Cierre reproducible de VE-SNP 0.1 | Manifiesto, hashes, versiones, recursos y revisión | Resultado no repetible o parcial presentado como cierre |

## 13. Trabajo excluido de esta versión

- ensamblaje *de novo* y validación de ensamblados;
- consenso diploide o por haplotipo;
- incorporación de minimap2, FreeBayes o un faseador al simulador;
- calibración de perfiles técnicos;
- contaminación, quimeras, *index hopping*, duplicados, adaptadores u otros
  artefactos;
- variantes estructurales, repetición en tándem y comparación de llamadas
  complejas;
- multialelia e *indels* en la verdad hasta resolver VE-D11; en la llamada se
  inventarían siempre y, cuando la descomposición de VE-D13 los resuelva en SNP
  con verdad correspondiente, dejan de ser falsos positivos;
- VE.2: haplotipos, LD y fase inferida, hasta disponer de salida F5 aprobada y
  faseador evaluado; no bloquea VE 0.1;
- ensayo futuro de regiones con `N`, máscara, repetición o baja llamabilidad:
  usará un contrato separado, métricas estratificadas de mapeo, cobertura, no
  llamada, FP y FN, y no modificará los denominadores de VE-SNP 0.1. VE-SNP 0.1
  no las estudia, pero tampoco las esconde: cuando una de ellas explique una
  discrepancia, §6.4.2 obliga a acreditarla con el dato observado. Lo que queda
  diferido es caracterizarlas, no detectarlas;
- afirmaciones de fase de largo alcance no identificable con PE150;
- comparaciones entre la referencia reducida y coordenadas del tomate original;
- publicación, instalación o descarga de herramientas.
