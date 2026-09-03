# Borrador de especificación, arquitectura mínima y plan de validación externa FASTQ→VCF 0.1

## 1. Estado y autoridad

- **Fecha:** 2 de septiembre de 2026.
- **Estado:** propuesta inicial no aprobada y no ejecutable.
- **Finalidad:** definir una validación externa que parta únicamente de los
  FASTQ, la referencia y el VCF de verdad publicados por el simulador.
- **Relación con el simulador:** trabajo paralelo a F5. No modifica F5, no
  presupone que F5 esté terminado y no forma parte de su cierre.
- **Autoridad normativa:** ninguna. La especificación vigente excluye la
  llamada de variantes y la evaluación de herramientas externas y las remite a
  DTF-024. Este borrador propone el contrato de ese banco externo; no modifica
  la especificación, la arquitectura ni la hoja de ruta del simulador.
- **Coordenada propuesta:** `VE` (*validación externa*). Debe aprobarse o
  sustituirse antes de ejecutar el plan. No se usa `F9.5`, porque F9 depende de
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

## 2. Resultado que se quiere demostrar

Un flujo ajeno al simulador debe poder:

1. verificar que referencia, FASTQ y VCF de verdad pertenecen a la misma
   ejecución y a las mismas muestras;
2. alinear con minimap2 los FASTQ Illumina paired-end contra la referencia
   exacta publicada por el simulador;
3. producir alineamientos por muestra, legibles e indexables con herramientas
   independientes;
4. llamar conjuntamente las variantes de las 20 muestras con FreeBayes;
5. comparar el VCF llamado con el VCF de verdad mediante una equivalencia
   semántica acordada, no por identidad byte a byte;
6. evaluar por separado SNP independientes y la conservación observable de la
   señal de haplotipos o desequilibrio de ligamiento (LD).

El **desequilibrio de ligamiento (LD)** es la asociación no aleatoria entre
alelos de loci distintos. Su medida y su interpretación dependen de si se usan
haplotipos faseados o genotipos sin fase; ambas evaluaciones no son
intercambiables.

## 3. Qué significa «recomponer el genoma de referencia»

La expresión admite tres significados distintos. El ensayo debe nombrar cuál
demuestra:

| Significado | Qué observa | Alcance de este borrador |
|---|---|---|
| Cobertura y alineamiento | Las lecturas pueden situarse sobre la referencia y cubrir sus bases con los criterios acordados | Incluido |
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
| Datos de entrada | FASTQ, referencia y VCF de verdad generados y publicados por el simulador |
| Referencia | 500.000 bases derivadas de la referencia de tomate disponible |
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

### 5.2. Herramientas observadas el 2 de septiembre de 2026

| Herramienta | Estado local observado | Consecuencia |
|---|---|---|
| minimap2 | `2.26-r1175`, disponible | Permite preparar el mapeo, sujeto a aprobar preset y parámetros |
| samtools | `1.21`, disponible | Permite ordenar, indexar y comprobar BAM |
| bcftools | `1.21`, disponible | Es anterior a la versión 1.24 exigida por la puerta externa DTC-066 del simulador |
| bgzip/tabix | HTSlib `1.19`, disponibles | Son anteriores a la versión 1.24 exigida por DTC-066 |
| FreeBayes | no disponible | Bloquea la ejecución, no la redacción del contrato |
| WhatsHap/SHAPEIT4/Beagle | no disponibles | No debe darse por aprobado ningún faseador |
| hap.py/RTG vcfeval/truvari | no disponibles | Debe decidirse el comparador antes de prometer equivalencia compleja |

No se ha instalado ni descargado ninguna herramienta.

### 5.3. Bloqueo de la primera prueba sin errores

La configuración F4 pública permite fijar a cero quimeras, contaminación e
*index hopping*. Sin embargo, el ejecutor aplica actualmente
`ErrorQualitySettings::default()`: escala de error 1 y el perfil incorporado de
calidades. Los campos de escala de error no forman parte del TOML general v2.

Por tanto, el programa actual no ofrece a la vez PE y ausencia de errores de
secuenciación. F3 ofrece un perfil sin errores, pero es single-end. Antes de
ejecutar el ensayo debe aprobarse una de estas salidas, sin resolverla dentro de
este plan:

1. añadir en el documento propietario de F4 un modo PE técnicamente neutro y su
   contrato observable;
2. aceptar para el piloto el perfil PE actual y retirar la condición «sin
   errores de secuenciación»;
3. ejecutar primero un control SE neutro de F3, dejando el control PE solicitado
   para después.

Las opciones 2 y 3 cambian el ensayo pedido. La opción 1 cambia el simulador y
no queda autorizada por este borrador.

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
- informes por muestra de lecturas totales, mapeadas, primarias, secundarias,
  suplementarias, correctamente emparejadas y no mapeadas;
- distribución de profundidad y cobertura según el criterio pendiente;
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
- correspondencia de genotipos por muestra, separando alelos, dosis, datos
  ausentes y fase;
- tratamiento explícito de registros multialélicos y *indels*;
- regiones excluidas por falta de llamabilidad o mapeo ambiguo;
- campos del caller que son informativos pero no forman parte de la verdad.

El VCF de verdad debe permanecer inmutable. La normalización produce una vista
temporal o derivada identificada; nunca repara el original.

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

El banco externo queda fuera del runtime Rust del simulador. Sus componentes
propuestos son:

| Componente | Responsabilidad única | Entrada | Salida o error observable |
|---|---|---|---|
| Admisor del conjunto | Resolver manifiestos, hashes, muestra y referencia antes de ejecutar herramientas | Publicación del simulador | Inventario admitido o lista de incompatibilidades |
| Preparador de referencia | Verificar FASTA, crear índices externos y registrar su versión | FASTA derivado | Índices identificados o fallo de referencia |
| Adaptador minimap2 | Alinear cada par R1/R2 y fijar el grupo de lectura de la muestra | FASTA y FASTQ | SAM por flujo y diagnóstico de proceso |
| Adaptador samtools | Ordenar, validar e indexar cada alineamiento sin perder lecturas | SAM | BAM/BAI por muestra e informes |
| Adaptador FreeBayes | Realizar una llamada conjunta de las 20 muestras | FASTA y BAM | VCF bruto y diagnóstico de proceso |
| Normalizador/comparador | Crear vistas comparables sin alterar la verdad y clasificar discordancias | VCF de verdad y llamado | Tabla de concordancia y métricas acordadas |
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

- Toda etapa usa el mismo FASTA exacto; no se permite mapear contra el tomate
  original y comparar contra coordenadas de la referencia reducida.
- Las 20 columnas del VCF de verdad, los `SM` de BAM y las 20 muestras llamadas
  deben formar el mismo conjunto biyectivo.
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

Los umbrales numéricos permanecen pendientes. Sí pueden aprobarse ya las
condiciones estructurales siguientes:

1. **Admisión:** hashes, referencia y 20 identidades son coherentes; todos los
   artefactos requeridos existen y los manifiestos están completos.
2. **Alineamiento:** cada par FASTQ se consume una sola vez; el total de
   lecturas de entrada se explica como alineado o no alineado; cada BAM es
   válido, coordinado, indexable y conserva la muestra correcta.
3. **Cobertura:** se informa la distribución acordada por muestra y región, no
   solo una media global.
4. **Llamada:** FreeBayes finaliza sobre los 20 BAM y produce un VCF legible,
   normalizable y asociado a la referencia correcta.
5. **SNP:** cada verdad, llamada y discordancia puede trazarse a locus, alelo,
   genotipo y muestra. Sensibilidad, precisión y concordancia no se publican sin
   sus numeradores y denominadores.
6. **LD/fase:** el informe separa la conservación estadística de LD de la
   recuperación de fase. No declara haplotipos exactos si no existe evidencia
   que conecte los loci.
7. **Reproducibilidad:** una repetición con las mismas entradas y versiones
   reproduce los resultados lógicos definidos; diferencias de metadatos o
   compresión se registran por separado.

## 9. Decisiones abiertas y recomendación inicial

| ID | Decisión pendiente | Opciones o recomendación no aprobada | Bloquea |
|---|---|---|---|
| VE-D01 | Propiedad y coordenada | Recomendada: banco externo DTF-024 con coordenada `VE`, no F5 ni F9.5 | Numeración ejecutable y documentos de origen |
| VE-D02 | PE sin errores de secuenciación | Resolver la incompatibilidad de §5.3. Recomendada: definir primero el contrato PE neutro; no fingir neutralidad con el perfil actual | Generación de VE-SNP y VE-LD |
| VE-D03 | Contenido exacto de VE-SNP | Recomendada: solo SNP bialélicos, sin multialelia ni *indels*, en la primera prueba | Configuración y comparador |
| VE-D04 | Densidad y selección de loci | Fijar número/densidad, frecuencias y si se restringen a regiones únicas A/C/G/T | Generación y criterio de exactitud |
| VE-D05 | Versiones y argumentos | Aprobar FreeBayes, minimap2, samtools y bcftools; decidir si se homologa todo con HTSlib 1.24 | Comandos reproducibles |
| VE-D06 | Normalización y emparejamiento VCF | Decidir recorte/alineamiento izquierdo, descomposición multialélica, comparación por alelo y herramienta | Métricas SNP/indel |
| VE-D07 | Regiones llamables | Fijar umbrales de profundidad, calidad de mapeo/base, repetición y tratamiento de `N` | Denominadores y umbrales |
| VE-D08 | Umbrales de aprobación | Fijar sensibilidad, precisión, concordancia genotípica, tasa de no llamada y cobertura | Cierre científico |
| VE-D09 | Métrica LD | Elegir estadístico, pares/bins y tolerancia. Recomendada para el piloto: comparar el mismo estadístico sin fase en verdad y llamada y etiquetarlo como conservación de señal, no fase | VE-LD |
| VE-D10 | Recuperación de fase | Decidir si se exige. Si se exige, evaluar faseador y limitar la afirmación a bloques conectados/identificables | Afirmación sobre haplotipos |
| VE-D11 | Multialelia e *indels* | Recomendada: excluirlos de VE-SNP 0.1 y crear después casos propios; no descartarlos silenciosamente del VCF | Comparación completa |
| VE-D12 | Repeticiones | Decidir si basta una semilla fija para el piloto o se exigirán réplicas | Generalización del resultado |

## 10. Plan propuesto

Los identificadores siguientes usan la coordenada provisional `VE`. No se
marcará ninguna casilla como completada ni se ejecutará un paso hasta aprobar
VE-D01 y las decisiones que lo bloqueen.

### VE.0 — Cierre del contrato del banco externo

- ☐ **VE.0.1. Aprobar alcance, vocabulario y afirmación científica.**
  - Depende de: ninguno.
  - Procede de: §§2–3 y VE-D01.
  - Cambia: documento propietario de DTF-024; no modifica la especificación del
    simulador.
  - Se comprueba con: acta que distingue alineamiento, consenso, recuperación
    de variantes, ensamblaje y fase.

- ☐ **VE.0.2. Resolver el escenario técnico PE neutro.**
  - Depende de: VE.0.1.
  - Procede de: §5.3 y VE-D02.
  - Cambia: solo el documento propietario que se apruebe; si exige ampliar F4,
    esa ampliación necesitará especificación, arquitectura y plan propios.
  - Se comprueba con: una configuración admitida cuyo manifiesto demuestre cero
    contaminación, quimeras, *index hopping* y errores realizados, sin usar F5.

- ☐ **VE.0.3. Fijar contrato de comparación y umbrales.**
  - Depende de: VE.0.1.
  - Procede de: §§6.4–6.6 y VE-D03–VE-D12.
  - Cambia: tabla normativa de métricas, denominadores, normalización, regiones
    llamables y tratamiento de fase/multialelia/*indels*.
  - Se comprueba con: ejemplos manuales de TP/FN/FP, genotipo discordante,
    registro multialélico, región no llamable y par LD que no admitan dos
    clasificaciones distintas.

> **PAUSA-VE.0.** No preparar scripts ni instalar herramientas hasta aprobar el
> alcance, resolver PE neutro y cerrar la comparación que determinará éxito o
> fallo.

### VE.1 — Conjunto pequeño, alineamiento y SNP independientes

- ☐ **VE.1.1. Congelar el entorno externo y sus comprobaciones.**
  - Depende de: PAUSA-VE.0 y VE-D05.
  - Procede de: §§5.2, 6.1 y 7.
  - Cambia: manifiesto de herramientas y receta de ejecución; no el runtime del
    simulador.
  - Se comprueba con: versiones exactas, hashes de ejecutables o imagen, ayuda
    de cada comando y un fixture diminuto válido e inválido. Debe incluir
    FreeBayes, hoy ausente.

- ☐ **VE.1.2. Derivar y congelar la referencia de tomate de 500.000 bases.**
  - Depende de: VE.1.1.
  - Procede de: F1.6, ACE-112/113 y §4.
  - Cambia: configuración de reducción del ensayo y su expediente externo.
  - Se comprueba con: ejecución del simulador, reconstrucción independiente
    FASTA fuente+mapa, longitud exacta 500.000, `samtools faidx`, SHA-256 y
    ausencia de mezcla entre coordenadas originales y derivadas.

- ☐ **VE.1.3. Generar VE-SNP con 20 individuos y 30× PE.**
  - Depende de: VE.0.2–VE.0.3 y VE.1.2.
  - Procede de: §§4, 6.1–6.2 y VE-D03/04.
  - Cambia: configuración de ensayo y directorio de resultados del simulador;
    no F5.
  - Se comprueba con: manifiestos y hashes; 20 pares R1/R2; VCF con las mismas
    20 muestras; cobertura solicitada 30×; artefactos técnicos realizados a
    cero según el contrato aprobado; validadores FASTQ/VCF externos.

- ☐ **VE.1.4. Alinear e indexar las 20 muestras con minimap2 y samtools.**
  - Depende de: VE.1.1 y VE.1.3.
  - Procede de: §§2, 6.3/6.5 y arquitectura §7.
  - Cambia: índices de referencia, BAM/BAI e informes externos por muestra.
  - Se comprueba con: `samtools quickcheck`, cabeceras `@SQ/@RG`, `SM` exacto,
    conteo de lecturas explicable e informes de alineamiento y profundidad. Los
    argumentos concretos se fijan en VE-D05; no se improvisan al ejecutar.

- ☐ **VE.1.5. Llamar conjuntamente VE-SNP con FreeBayes.**
  - Depende de: VE.1.4.
  - Procede de: §§2 y 6.3–6.4.
  - Cambia: VCF bruto y log de FreeBayes.
  - Se comprueba con: proceso exitoso, 20 muestras correctas, referencia y
    contigs compatibles, VCF formalmente legible y ausencia de BAM omitidos o
    duplicados.

- ☐ **VE.1.6. Normalizar y comparar VE-SNP con el VCF de verdad.**
  - Depende de: VE.0.3 y VE.1.5.
  - Procede de: §§6.4–6.5 y VE-D06–VE-D08.
  - Cambia: vistas VCF derivadas, tabla de discordancias e informe de métricas.
  - Se comprueba con: oráculos manuales, numeradores/denominadores, listado
    exhaustivo por muestra/locus y cumplimiento de los umbrales aprobados. No
    compara `GPS_CASE`, `QUAL` o campos del caller salvo decisión expresa.

> **PAUSA-VE.1.** Revisar los BAM, el VCF llamado, todas las discordancias y los
> denominadores. No pasar a LD si el caso independiente no explica cada fallo.

### VE.2 — Haplotipos y LD

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

> **PAUSA-VE.2.** Decidir si la evidencia demuestra recuperación de genotipos,
> conservación de señal LD o fase por bloque. No fusionar esas conclusiones.

### VE.3 — Repetición y cierre

- ☐ **VE.3.1. Repetir los conjuntos aprobados y conservar recursos observados.**
  - Depende de: PAUSA-VE.2 y VE-D12.
  - Procede de: §§6.3 y 8.
  - Cambia: repetición o réplicas aprobadas y tabla de recursos.
  - Se comprueba con: igualdad lógica definida, hashes cuando procedan,
    `/usr/bin/time` u observador equivalente y tamaños por etapa.

- ☐ **VE.3.2. Ensamblar el expediente reproducible y revisar el cierre.**
  - Depende de: VE.3.1.
  - Procede de: §§6.3, 7.1 y 8.
  - Cambia: manifiesto e informe final del banco externo.
  - Se comprueba con: todos los enlaces y hashes resolubles, comandos y versiones
    completos, cero resultado parcial marcado como completo y revisión
    independiente de las afirmaciones científicas.

> **PAUSA-VE.3 / CIERRE.** Aprobar por separado: interoperabilidad FASTQ/BAM,
> recuperación SNP, conservación de LD y recuperación de fase si se ejecutó.

## 11. Recursos y tiempo de planificación

El único volumen calculable antes de ejecutar es el de 300 millones de bases y
aproximadamente un millón de parejas. El tamaño comprimido, el BAM, la memoria
de FreeBayes y el tiempo dependen de la complejidad de la referencia, densidad
de variantes, parámetros y versiones aún abiertas.

Para no presentar una estimación como resultado, VE.1 debe empezar con un
presupuesto operativo provisional y sustituirlo por medidas observadas. Como
reserva de agenda, no como requisito ni predicción, se propone una estación con
4 CPU, 8 GiB de RAM, 10 GiB libres y una ventana de 2 horas para el primer
recorrido. VE.1.4–VE.1.6 registrarán uso real por etapa. Si el preflight o las
medidas superan esa reserva, el ensayo se detiene y se revisa el presupuesto; no
reduce cobertura, individuos ni evidencia automáticamente.

## 12. Matriz de evidencia y riesgos

| Hito | Resultado | Evidencia principal | Riesgo protegido |
|---|---|---|---|
| VE.0 | Contrato ejecutable | Decisiones aprobadas y casos manuales inequívocos | Medir algo distinto de lo que se afirma |
| VE.1 | Recuperación externa de SNP | FASTQ→BAM→FreeBayes→comparación por muestra/locus | Mezcla de referencia, muestra, representación VCF o denominador |
| VE.2 | Evaluación separada de LD y fase | Verdad F5 independiente y mismo estimador | Confundir genotipo, LD, haplotipo y fase inferida |
| VE.3 | Expediente reproducible | Manifiesto, hashes, versiones, recursos y revisión | Resultado no repetible o parcial presentado como cierre |

## 13. Trabajo excluido de esta versión

- ensamblaje *de novo* y validación de ensamblados;
- consenso diploide o por haplotipo;
- incorporación de minimap2, FreeBayes o un faseador al simulador;
- calibración de perfiles técnicos;
- contaminación, quimeras, *index hopping*, duplicados, adaptadores u otros
  artefactos;
- variantes estructurales, repetición en tándem y comparación de llamadas
  complejas;
- multialelia e *indels* hasta resolver VE-D11;
- afirmaciones de fase de largo alcance no identificable con PE150;
- comparaciones entre la referencia reducida y coordenadas del tomate original;
- publicación, instalación o descarga de herramientas.
