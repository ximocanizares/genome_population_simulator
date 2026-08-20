# Punto de reanudación del proyecto

> **Punto vigente (20 de agosto de 2026):** F1.1–F1.6 del
> [plan de reimplementación F1](plan-reimplementacion-f1-0.2.md) están
> completadas y todas sus pausas están aprobadas. H-F16-R1/R2 quedaron aplicados
> 2/2; las dos revisiones independientes posteriores fueron consecutivas y
> limpias. F1 está cerrada sobre el commit publicado `e14705a`; F2 todavía no
> se ha iniciado. La documentación F2 fue corregida después de DTC-078 y queda
> pendiente de nueva revisión adversarial y aprobación; la revisión 0.4 no la
> aprueba.
> Las indicaciones históricas
> posteriores de este documento no sustituyen este punto.

- **Fecha:** 20 de agosto de 2026
- **Implementación F1.5 histórica del checkpoint anterior:** `b5dc220 feat:
  complete F1.5 reference evidence`; no acredita el plan 0.2 vigente.
- **Base publicada de cierre F1:** `e14705a feat: complete minimal F1.6
  reference reductions`, publicada en `origin/main`.
- **Checkpoint de seguridad anterior:** `66ad8b2 checkpoint: preserve F1 and
  F2 work before simplification`.
- **Estado actual de ejecución:** la auditoría sistemática bloqueó el alcance
  anterior y recomendó simplificar en vez de seguir parcheando. DTC-078 fija el
  modelo de confianza local mínimo; DTC-076/DTC-077 quedan superadas. El plan
  vigente es `plan-reimplementacion-f1-0.2.md`. F1.1–F1.3 están
  completados. F1.3 conserva raíz de 256 bits y derivación BLAKE3/ChaCha12,
  activa solo `reference_generation` y `reference_reduction`, y mantiene
  intactos los catálogos v1–v3. La revisión independiente aprobó PAUSA-F1.3
  después de corregir destinos y símbolos de trazabilidad. F1.4 conserva solo
  manifiesto general `complete`/`failed`, catálogo único y enlaces por ID,
  wrapper privado 0700 con `INCOMPLETE` exterior al payload publicable y rename no-clobber; retiró unidades
  ficticias, capabilities, inventario de inodos/hardlinks/TOCTOU y el crate
  `jsonschema`. La primera ronda de seis hallazgos quedó aplicada 6/6, y la
  segunda revisión abrió R1–R3 y los tres están aplicados: cruce efectivo
  completo, fallback `failed` sin artefacto efectivo inválido y alineación de
  todos los diagnósticos activos con el catálogo v2. La tercera revisión abrió
  A1 (autoridad física FASTA/mapa) y A2 (diagnósticos `failed` contra catálogo
  v2); ambos están aplicados 2/2 con regresiones y gates completos. La cuarta
  revisión abrió Q1–Q4: relación exacta entre objetivo declarado y bases,
  fallback `failed` conservador, clasificación semántica exacta de salidas y
  paridad Rust/Python sobre todo `diagnostics-v2`. Los cuatro están aplicados;
  la suite completa ejecutó 138 pruebas correctamente y dejó un helper de
  subproceso ignorado en el padre. La quinta revisión detectó dos defectos:
  autoridad TOML incompleta y categorías solapables. Ambos se corrigieron con
  reconstrucción exhaustiva desde TOML y disjunción global; la suite actual
  ejecuta 140 pruebas y deja un helper ignorado. La sexta revisión detectó la
  falta de ligadura de `execution_id`, una divergencia del oráculo Python y la
  ruta reservada `INCOMPLETE`; las tres están corregidas. La séptima revisión
  fue funcionalmente limpia y detectó solo que TEST-REQ-SAL-03 no apuntaba a
  una prueba de la copia TOML y que REQ-MAN-32 conservaba un símbolo antiguo.
  Ambos enlaces están corregidos; la prueba nueva recorre el CLI real y cruza
  bytes, manifiesto, clase, ruta, tamaño y `sha256sum`. La suite actual ejecuta
  143 pruebas y deja un helper ignorado. F1.5 reutiliza el dominio, adaptador,
  generador y FAI ya correctos; completa procedencia FASTA, formato canónico,
  cabecera/descripción reversibles y alfabeto en la salida efectiva, ligados al
  FASTA físico y al catálogo mínimo. A y B se midieron tres veces; C queda
  programado. La revisión independiente de F1.5 quedó aplicada 3/3: ejecución
  explícita con raíz declarada y CSPRNG, testigo compartido para las cláusulas
  de referencia derivada y RSS corregido a 393832 KiB (aproximadamente 385
  MiB). PAUSA-F1.5 está aprobada. Sobre esa base se completó F1.6; los gates
  posteriores a H-F16-R1/R2 listan 153 pruebas: 152 superadas y un helper
  ignorado. Dos revisiones independientes consecutivas posteriores terminaron
  limpias. PAUSA-F1.6 está aprobada, F1 está cerrada y F2 no ha comenzado.
- **Estado de los agentes al pausar:** implementador y revisor F1.3 terminaron;
  la revisión adversarial F2 0.4 es histórica, anterior a DTC-078 y no aprueba
  el plan F2 corregido. F2 requiere una nueva revisión adversarial y aprobación
  explícita antes de iniciarse.
- **Estado de diseño:** DTPD-101, DTPD-116, DTPD-118, DTPD-119 y DTPD-121 están resueltas.
  DTC-056 fija las coordenadas comunes y DTC-057,
  `reduced_reference_selection_v1`; DTPD-119 queda cerrada mediante DTC-055
  después de una revisión independiente y dos rondas de correcciones. DTC-068
  sustituye la materialización pesada 1000 Genomes de F2 por los perfiles
  aproximados compactos `human_1000g_global_sfs_approx_v1` y
  `human_harris2017_snv12_approx_v1`; el VCF y perfiles exactos quedan en
  DTF-016/post-v1 y no bloquean F2.2. DTC-070 resuelve RAF2-029: el modo
  `exact` usa cuotas por mayores restos, permutación independiente por población
  y flujo separado para el valor; DTC-071 publica el catálogo RNG v4 sin
  circularidad y resuelve RAF2-028. DTC-072 fija los consumidores RNG no Poisson,
  DTC-073 hace efectivas las capacidades solo con su semántica y DTC-074 exige
  evidencia de plataforma detectada. DTC-075 resuelve RAF2-041 a RAF2-044:
  intervalos semiabiertos sin extremos continuos, selector indel dinámico,
  propietarios RNG únicos y `matrix_candidate_count` explícito. RAF2-030 queda resuelto mediante DTC-065/066: `GPS_CASE` es
  multivaluado y canónico, `truth.jsonl` conserva ALT↔casos y F2 no genera
  complemento; la primera cardinalidad positiva se define en F5.5. La revisión
  adversarial F2 0.3 queda corregida: RAF2-028 a RAF2-044 están
  resueltos documentalmente; el plan corregido tiene 130 pasos con cinco
  campos, seis pausas explícitas, 78 cierres, catorce compartidos y nueve
  parciales directos; requiere nueva revisión adversarial;
  REQ-REF-14 pasó justificadamente a F1.5 y F2.1 conserva su invariante,
  mientras REQ-CONF-38 pasó a F2.3 con individuos reales.
  RAF2-027 queda resuelto por el cierre y publicación de F1 en `e14705a`. La
  revisión adversarial F2 0.4 fue anterior a DTC-078 y queda sustituida como
  aprobación por la corrección vigente; iniciar F2 requiere una nueva revisión
  adversarial y aprobación explícita del usuario.
- **Copia de la conversación:**
  [chat-reimplementacion-f1-2026-08-15.md](../docs/chat-reimplementacion-f1-2026-08-15.md),
  [chat-dtpd-119-2026-08-12.md](../docs/chat-dtpd-119-2026-08-12.md) y
  [chat-generacion-automatica-bloques-ld-2026-08-13.md](../docs/chat-generacion-automatica-bloques-ld-2026-08-13.md).

## Cierre vigente de F1

F1.6.1–F1.6.5 están completadas y publicadas en `e14705a`. La selección publica sus
parámetros reconstruibles; FASTA y mapa se contrastan mediante reconstructor,
`cmp`, `sha256sum`, `samtools` y validación JSON independiente; y el expediente
vigente contiene exactamente 54 REQ distribuidos 4/1/6/22/21. Los dos defectos
encontrados durante la implementación tienen regresiones rojo-verde y están
corregidos. Todos los gates propios de F1.6.5 pasan: 152 pruebas superadas, un
helper ignorado, doctest, fmt, Clippy estricto, doc, release, metadata, audit y
`git diff --check`.

H-F16-R1/R2 están aplicados 2/2. Las dos revisiones independientes ejecutadas
después de esas correcciones terminaron limpias. PAUSA-F1.6 y F1 están
aprobadas; no iniciar F2 sin la autorización expresa del usuario. Conservar la superficie TOML mínima, los dos únicos
ámbitos RNG activos y la ausencia de `ResolutionDomain`, `resolution-evidence`,
eventos/manifiestos sintéticos, normalización poblacional y hardening frente a
atacante local.

Las rondas H-F16 descritas más abajo son historia del checkpoint `66ad8b2`; no
constituyen el plan vigente ni obligan a conservar código trasladado.

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

F1.6 conserva H-F16-R09-R2/R3 y cierra H-F16-R11–R16 tras la decimocuarta ronda.
La matriz ejecutable recorre una/dos reducciones, índices 0/1 y las tres fronteras;
los corpus vuelven a ejecutar positivos completo/fallido, mezclas contextuales,
contradicciones entre colecciones y Python frente a Rust. El catálogo trazable
comprueba los 22 nombres contra sus fuentes y contra el listado ejecutado.

H-F16-R13 añade un corpus real con cero, una y dos reducciones y cierres
completo/fallido. Referencias, recurso autoritativo, configuración efectiva y
manifiestos de evento se cruzan íntegramente contra el catálogo; 29 mutaciones
individuales se contrastan con un oráculo semántico independiente ampliado.

H-F16-R14 cierra el resto del grafo dinámico interno: configuraciones, unidad de
ejecución, biyección eventos/unidades, categorías de salida, propietarios,
clases y unicidad global de rutas. Su corpus catalogado recorre ocho positivos y
24 contradicciones; el expediente incluye un inventario que separa qué igualdad
pertenece al esquema, al grafo del manifiesto o a la relectura física de otros
ficheros.

H-F16-R15 añade el contrato autoritativo de diagnósticos y completa la semántica
de catálogos de secuencias y reducciones: siete positivos y 27 mutaciones
individuales se recorren desde un catálogo literal y una tabla declarativa de
oráculos. Se permiten códigos compartidos por unidades fallidas, pero no códigos
ausentes, duplicados ni huérfanos; particiones, intervalos, sumas, objetivo,
región y residual deben reconstruir exactamente la referencia derivada.

H-F16-R16 sustituye el `f64` del oráculo porcentual por aritmética decimal larga
independiente y cierra la relectura semántica del paquete. Raíz y subsemillas se
derivan con catálogo v3; orden y cabeceras/descripciones son reversibles;
normalizaciones reconstruyen Neumaier; y cada evento coincide entre
configuración efectiva, manifiesto general y manifiesto de evento. El corpus
literal contiene seis positivos y 20 mutaciones correlacionadas.

H-F16-R17 valida la gramática decimal TOML antes de retirar separadores, relee
el FASTA autoritativo mediante `openat`/`O_NOFOLLOW` y cruza sus bytes con ambos
documentos JSON. También reconstruye la configuración efectiva desde las copias
TOML y cierra identidad, orden, dominio, regla ganadora, selectores, procedencia
y colecciones exclusivas. Su corpus literal contiene seis positivos y 34
negativos; sus validadores externos, incluida OCI con Docker 29.6.1, superaron
la validación.

H-F16-R18 conserva el contexto completo de cada `ResolutionRequest` en el valor
efectivo, valida que `object_id` pertenezca a una entidad declarada compatible y
vuelve a resolver la regla dominante. El oráculo independiente parsea reglas e
identidades desde TOML y calcula aplicabilidad/dominancia sin llamar a
producción. Dos mutaciones con tamaño/SHA recalculados cierran la copia de regla
base y el objeto ajeno; las pruebas adicionales cubren regla combinada, orden y
máximas incomparables.

H-F16-R19–R21 contrastan cada dimensión presente con los ejes permitidos del
parámetro, permiten jerarquías completas sin mezclar ejes no relacionados y
exigen que posición/locus lleven referencia y secuencia declaradas y queden en
rango. El oráculo independiente usa una tabla parámetro→ejes propia y reconstruye
longitudes desde el TOML. Dos mutaciones correlacionadas recalculan tamaño y
SHA-256; la suite suma 189 pruebas y los validadores externos, incluida OCI con
Docker 29.6.1, están verdes.

H-F16-R22–R24 añaden la matriz ejecutable de ocho parámetros por tres niveles
espaciales, exigen el objeto principal correspondiente y sustituyen la
inferencia desde listas planas por `ResolutionDomain`. F1 conserva el contrato
tipado y un snapshot versionado, pero no introduce poblaciones, fases
biológicas ni loci reales: los proveedores pertenecen a F2 y fases posteriores.
El oráculo relee por separado todas las copias de escenarios y agrega sus fases.
Ocho mutaciones correlacionadas y positivos completo/fallido con cero/una/dos
reducciones pasan producción y oráculo independiente según corresponda. La
suite suma 193 pruebas; puertas Cargo, validadores externos y OCI con Docker
Engine 29.6.1, commit `8ec5ab3`, están verdes.

H-F16-R25–R27 publican `resolution-evidence` v1 como fuente autoritativa
independiente de dominio y solicitudes, enlazada por identidad física y releída
sin seguir enlaces; la configuración efectiva queda como derivada. Los IDs de
fase son globalmente únicos entre escenarios y el diagnóstico conserva ambas
fuentes. El orden de parámetros usa una tabla literal
`parameter_order_version=1`, con implementación separada en el oráculo. Tres
pruebas nuevas y doce mutaciones del corpus cubren TOCTOU, enlaces, cambios
correlacionados, completo/fallido, colisión de fases y permutaciones. La suite
suma 196 pruebas; puertas Cargo, cuatro esquemas externos, cinco artefactos,
`samtools faidx` y OCI con Docker Engine 29.6.1 están verdes.

La puerta física añade una mutación que mueve un directorio intermedio real y
lo sustituye por un *symlink*: falla antes del publicador, no crea destino y
conserva `INCOMPLETE`. La canonicalización alcanza directamente `a//b`, `.`,
`..`, separador final y barra inversa; ya no se atribuye al alias de manifiesto
un rechazo que en realidad se producía primero por SHA-256.

**Primera acción al reanudar:** solicitar o confirmar la autorización para
iniciar F2 sobre la base inmutable `e14705a`. F1 no requiere nuevas correcciones.

El cierre de F1 está publicado en `origin/main` mediante `e14705a`. La
validación JSON externa no es una dependencia runtime de Rust. Al retomar,
partir de ese commit y preservar su expediente de 54 requisitos.

Queda una incoherencia futura no bloqueante para esta pausa: el preplan LD aún
reserva `subseed-scopes-v3.tsv`, pero F1.6 ya publica v3 con el ámbito
`reference_reduction`. Antes de convertir el preplan LD en plan ejecutable se
debe asignar la siguiente versión disponible o integrar formalmente los nuevos
ámbitos sin reescribir v3.

El preplan de generación automática de bloques de LD no tiene elecciones
funcionales o técnicas materiales abiertas y conserva el `profile_id`
definitivo `human_ld_ragsdale_approx_v1`. Su próximo trabajo, cuando vuelva a
ser prioritario, será trasladar el borrador a especificaciones y trazabilidad,
redactar arquitectura y decisión técnica y solo después crear el plan ejecutable
de F6.4.

Los acuerdos y propuestas de esta conversación se conservan en
[borrador-generacion-automatica-bloques-0.1.md](borrador-generacion-automatica-bloques-0.1.md).

## Estado conservado del preplan LD

La idea separada, todavía sin requisitos, generará automáticamente bloques, SNP,
haplotipos y frecuencias y producirá el mismo paquete canónico de cinco ficheros.
Ya se acordó usar `r²` haplotípico faseado, dependiente de la
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
