# Punto de reanudación del proyecto

- **Fecha:** 14 de agosto de 2026
- **Último cierre documental:** `docs: resolve haplotype panel design`.
- **Estado:** DTPD-118 y DTPD-119 están resueltas y publicadas. DTPD-119 queda
  cerrada mediante DTC-055 después de una revisión independiente y dos rondas
  de correcciones.
- **Copia de la conversación:**
  [chat-dtpd-119-2026-08-12.md](../docs/chat-dtpd-119-2026-08-12.md) y
  [chat-generacion-automatica-bloques-ld-2026-08-13.md](../docs/chat-generacion-automatica-bloques-ld-2026-08-13.md).

## DTPD-119: decisiones ya aprobadas

- El panel contiene bloques de LD locales independientes, no haplotipos
  cromosómicos completos. La identidad es `(block_id, haplotype_id)`.
- La población inicial no impone LD ni continuidad de donante entre bloques.
  Esa dependencia cromosómica queda en DTF-044. La recombinación fundadora
  dentro de cada bloque permanece disponible exactamente como define DTC-055.
- Cada bloque declara referencia, secuencia y límites explícitos; los bloques no
  se solapan y todos sus sitios quedan dentro del intervalo.
- Dentro del bloque se usan haplotipos locales. Fuera de los bloques, los SNP se
  muestrean independientemente desde la matriz `alt_af` de DTC-054.
- Un locus `alt_af` dentro del intervalo de un bloque es un error, aunque no sea
  un sitio del panel; no existe prioridad ni sobrescritura.
- La tabla de alelos es densa: cada haplotipo declara un `allele_index` para cada
  sitio, con `0=REF` y `1..N=ALT`. Ausencias y duplicados son errores.
- Esta versión del panel solo admite SNP A/C/G/T suministrados, bialélicos o
  multialélicos. Las demás clases quedan para una ampliación posterior.
- Las frecuencias son densas por población, bloque y haplotipo, admiten 0 y 1 y
  suman 1 conforme a DTC-049. No hay valores globales ni herencia entre
  poblaciones.
- El paquete canónico contiene exactamente cuatro TSV y
  `panel.metadata.json`; el JSON enlaza tablas y referencia mediante SHA-256 sin
  duplicar filas.
- El paquete solo admite una referencia: todos los bloques comparten el
  `reference_id` del metadato; el SHA-256 recalculado del FASTA coincide con el
  digest declarado; y cada `sequence_id` existe en ese FASTA.
- La orden Rust `prepare-haplotype-panel` prepara y publica una copia canónica
  atómica sin modificar entradas ni sobrescribir. Runtime recalcula y verifica;
  Marimo invoca la misma orden. VCF queda en DTF-041.
- Los vectores HAP deben ser únicos dentro de cada bloque, incluso con frecuencia
  0; no se fusionan IDs ni frecuencias.
- Cada bloque se trata como locus multialélico. Reutiliza HWE para P=1..10 y F
  multialélico solo diploide; la distribución explícita local se admite
  únicamente con construcción intacta.
- Cada bloque exige al menos dos sitios y dos HAP distintos. Se admiten
  poblaciones con un único HAP de frecuencia 1 y sitios fijados en ALT; un sitio
  `REF` en todos los HAP se rechaza como redundante.
- `haplotype_construction` intacto/recombinante se configura por
  población×bloque y se separa de `frequency_sampling`
  probabilístico/controlado. Intacto prohíbe `r`; recombinante exige `r` en
  (0,1] y solo admite muestreo probabilístico; el control recombinante queda en
  DTF-042.
- El modo probabilístico sortea copias con reemplazo. El controlado asigna un
  multiconjunto exacto que suma C y rechaza incompatibilidades sin ajustes.
- El modo intacto asigna exactamente un HAP completo por copia y bloque, conserva
  el HAP fundador y separa las mutaciones de novo posteriores.
- La alternativa recombinante inicializa directamente cada individuo×bloque y no
  representa una meiosis histórica ni usa un depósito prototipo intermedio.
- `recombination_probability_per_interbase=r` es finita en (0,1]. Se sortean de
  forma independiente las L−1 uniones: K~Binomial(L−1,r). K=0 conserva la copia;
  K≥1 usa su HAP inicial como primaria, sortea una secundaria por frecuencias y
  crea el mosaico con puntos únicos, ordenados e internos mediante el constructor
  común y la variante uniforme cerrados por DTC-055. DTPD-111 solo conectará
  mapas y semántica meiótica en F7.1 y no bloquea F6.4.
- Un punto x representa la unión x|x+1. La salida comienza con la primaria,
  alterna segmentos *forward* y no ajusta a SNP; también conserva en traza cruces
  invisibles o cuyo efecto alélico se cancela.
- En diploides, Bernoulli(F) ejecuta un productor y duplica su salida o ejecuta
  dos independientes; con P distinta de 2 ejecuta P y rechaza F. F no es
  autofecundación. Las C copias finales producen frecuencias realizadas `n_h/C`.
- Primaria y secundaria proceden siempre del panel original; no hay
  recombinación acumulada ni dependencia del orden. Un vector conocido reutiliza
  su HAP y uno nuevo deriva ID por DTC-039 desde panel+bloque+vector.
- El RNG DTC-037 separa F, primaria, puntos y secundaria por población, bloque,
  individuo y copia. Se distinguen productores con cruce y copias finales
  derivadas y se conserva su procedencia completa.
- Son obligatorias `founder_block_assignments.tsv`,
  `recombination_segments.tsv` y `realized_haplotype_frequencies.tsv`, enlazadas
  por el manifiesto. Marimo prepara/valida mediante Rust, previsualiza el resumen
  y enlaza estas salidas sin cargarlas completas.
- Las tres salidas tienen cabeceras, claves, orden, tipos y serialización TSV
  canónicos exactos; clase y versión residen en sus enlaces del manifiesto, no
  en columnas adicionales.
- Los cuatro TSV del panel tienen cabeceras, serialización UTF-8/LF/tabulador,
  tipos, claves y orden canónico exactos; `alt` separa varios ALT con comas sin
  espacios y su orden define `allele_index`.
- Todo productor posee `process_id`. `origin_kind` es `intact` para K=0 o
  `recombined` para K≥1. La frecuencia final separa recuentos totales, intactos y
  recombinados.
- `intact` admite HWE, F diploide o distribución genotípica explícita;
  `recombinant` solo HWE o F diploide. La distribución explícita recombinante
  queda en DTF-043.

## Punto exacto de reanudación

Continuar la decisión separada sobre la generación automática de bloques de LD.
La geometría, los perfiles, las frecuencias marginales, el constructor
haplotípico, la puntuación, el calendario de recocido y el presupuesto de
movimientos ya están acordados. También quedan cerrados el preflight conservador
de imposibilidad de LD, Hyndman–Fan Type 7, la representabilidad exacta de
cobertura, el aislamiento poblacional y los seis ámbitos físicos RNG con su
futuro catálogo v3. No queda ninguna cuestión abierta dentro de DTPD-119 ni una
elección funcional o técnica sustantiva en el preplan de generación automática.
El nombre físico `human_ld_ragsdale_approx_v1` queda ratificado como
`profile_id` definitivo.

El próximo trabajo no es otra decisión funcional: hay que trasladar el borrador
a especificaciones y trazabilidad,
redactar la arquitectura y la decisión técnica y, solo después, convertir el
preplan F6.4 en plan ejecutable.

Los acuerdos y propuestas de esta conversación se conservan en
[borrador-generacion-automatica-bloques-0.1.md](borrador-generacion-automatica-bloques-0.1.md).

## Siguiente conversación

Discutir, como idea separada y todavía sin requisitos, la generación automática
de bloques, SNP, haplotipos y frecuencias, produciendo el mismo paquete canónico
de cinco ficheros. Ya se acordó usar `r²` haplotípico faseado, dependiente de la
distancia física, mediante perfiles por población. El programa incluirá el
perfil sintético aproximado `human_ld_ragsdale_approx_v1` como ejemplo
seleccionable pero nunca predeterminado; un perfil empírico derivado de 1000
Genomas podrá añadirse en el futuro. Cada población podrá cargar su propio
perfil. En esta versión los bloques se mantendrán independientes entre sí
y DTF-044 reserva su LD posterior; la recombinación fundadora intrabloque de
DTC-055 no cambia. La colocación inicial será uniforme sobre inicios elegibles;
DTF-045 reserva regiones pericentroméricas con baja recombinación y bloques
mayores. El constructor haplotípico aprobado usa intercambios que conservan las
marginales y recocido simulado para ajustar los cuantiles agregados de `r²`;
DTF-046 reserva el signo del LD. Los resultados se resumen en cuatro TSV de
geometría, frecuencias realizadas, perfil de LD e intentos, enlazados por el
manifiesto; DTF-047 reserva una tabla futura con el `r²` de cada pareja.
La representabilidad de cobertura ya queda separada de la colocación: un grafo
exacto de residuos módulo la menor longitud admisible calcula la menor suma por
residuo y reconstruye un testigo exacto. El testigo no sustituye el muestreo
ponderado; pasar la prueba no garantiza que el layout fragmentado pueda colocar
la lista. El agotamiento de layout y la falta de recursos del preflight se
diagnostican por separado y nunca como suma imposible. Para `block_count`, solo
superar `sum(floor(segment_length/min_allowed_length))` demuestra imposibilidad.
El perfil
`human_ld_ragsdale_approx_v1` ya queda definido como aproximación
sintética: usa la tendencia `r²_÷` de Ragsdale y Gravel, una distribución beta
por población e intervalo y el valor global/versionado `kappa = 8.0`; nunca es
predeterminado ni se presenta como cuantiles empíricos. El preflight enumera el
dominio discreto de `r²` y deriva cotas Type 7 conservadoras; fixtures exhaustivos
pequeños protegen frente a falsos imposibles. Los ámbitos físicos aprobados son
`ld_geometry`, `ld_site_alleles`, `ld_population_frequencies`,
`ld_haplotype_initialization`, `ld_annealing_proposals` y
`ld_annealing_acceptance`, todos con `derivation_version = 2`. La geometría es
común e inmutable; frecuencias, cuantización, inicialización, cálculo de
cuantiles Type 7, score y recocido son independientes por población y usan un
`population_attempt_index`. Cambiar, añadir o reordenar una población no altera
los resultados no relacionados. Los catálogos v1 y v2 no se modifican ni
reinterpretan. El futuro `catalogs/subseed-scopes-v3.tsv` tendrá
`catalog_version = 3`, `derivation_version = 2`, los ámbitos anteriores
necesarios y los seis ámbitos LD con vectores normativos propios. Una ejecución
automática seleccionará coherentemente catálogo 3 y derivación 2, sin mezcla
silenciosa, y lo registrará conforme a DTC-037/DTC-038. El fichero real no se
creará antes de cerrar las puertas normativas y autorizar la implementación.

La puntuación, el calendario y el presupuesto del recocido ya están cerrados.
Cada población usa
`population_haplotype_cells=sum(C*snp_count)` sobre sus bloques y
`max_annealing_moves=min(10_000_000,max(100_000,100*population_haplotype_cells))`.
El límite cuenta propuestas y no es configurable. Cada población dispone de
`max_haplotype_attempts=1000`; al fallar reintenta solo ella. Si una agota el
límite, falla el panel completo sin resortear las demás y sin publicación
parcial. La unión densa de HAP se forma únicamente cuando todas han pasado.

La orden futura deberá producir exactamente el paquete canónico de DTC-055 y
pasar el mismo validador que un panel suministrado. Esta capacidad no forma parte
de DTC-055 y, aunque su encaje preliminar es F6.4, todavía no tiene REQ, ACE,
DTC/DTPD ni plan ejecutable.
