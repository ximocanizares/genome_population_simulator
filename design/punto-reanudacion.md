# Punto de reanudación del proyecto

- **Fecha:** 13 de agosto de 2026
- **Último cierre documental:** `docs: resolve haplotype panel design`.
- **Estado:** DTPD-118 y DTPD-119 están resueltas y publicadas. DTPD-119 queda
  cerrada mediante DTC-055 después de una revisión independiente y dos rondas
  de correcciones.
- **Copia de la conversación:**
  [chat-dtpd-119-2026-08-12.md](../docs/chat-dtpd-119-2026-08-12.md).

## DTPD-119: decisiones ya aprobadas

- El panel contiene bloques de LD locales independientes, no haplotipos
  cromosómicos completos. La identidad es `(block_id, haplotype_id)`.
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

Iniciar una decisión separada sobre la generación automática, por el propio
programa, de bloques, SNP, haplotipos y frecuencias. No queda ninguna cuestión
funcional o técnica abierta dentro de DTPD-119.

## Siguiente conversación

Discutir, como idea separada y todavía sin requisitos, la generación automática
de bloques, SNP, haplotipos y frecuencias a partir de un valor o rango de LD,
produciendo el mismo paquete canónico de cinco ficheros. La propuesta inicial es
usar `r²` entre SNP consecutivos dentro de cada bloque, pero no está aprobada.
Antes de especificar o planificar hay que decidir:

1. el estadístico de LD y las parejas de SNP a las que se aplica;
2. si se declara un valor o un intervalo objetivo;
3. el número/densidad de SNP y la distribución del tamaño de los bloques;
4. el espectro de frecuencias alélicas y el modelo que construye haplotipos
   compatibles;
5. la tolerancia, los casos imposibles y la política de intentos/fallo;
6. la semilla, los ámbitos y la versión del algoritmo;
7. los valores de LD solicitados y realizados que deben registrarse y validar.

La orden futura deberá producir exactamente el paquete canónico de DTC-055 y
pasar el mismo validador que un panel suministrado. Esta capacidad no forma parte
de DTC-055 y todavía no tiene REQ, ACE, DTC, DTPD ni asignación de fase.
