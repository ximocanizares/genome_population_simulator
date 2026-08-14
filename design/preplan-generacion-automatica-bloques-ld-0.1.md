# Preplan de generación automática de bloques de LD 0.1

## 1. Cabecera y estado

- **Fecha:** 13 de agosto de 2026.
- **Estado:** preplan bloqueado; no ejecutable.
- **Encaje candidato:** F6.4, antes de admitir el panel y crear la población
  fundadora.
- **Fuente de conversación:**
  [borrador-generacion-automatica-bloques-0.1.md](./borrador-generacion-automatica-bloques-0.1.md).
- **Contrato de salida ya existente:** paquete canónico de cinco ficheros de
  DTC-055.

Este documento conserva el encaje de planificación y las decisiones ya
aprobadas sin convertir el borrador de conversación en una especificación. No
contiene pasos autorizables de implementación: la puerta de entrada de
`planificar-implementaciones` permanece cerrada hasta que existan requisitos,
criterios de aceptación, arquitectura y decisiones técnicas vigentes para el
generador.

## 2. Encaje en la hoja de ruta

La capacidad pertenece a F6.4 porque produce una fuente poblacional de bloques
haplotípicos y su primer consumidor es la creación de la población fundadora.
No pertenece a F2: esa fase representa loci marginalmente independientes y
excluye expresamente paneles de haplotipos y LD.

El orden lógico provisional, sin asignar todavía componentes o archivos de
código, es:

```text
configuración y perfiles aprobados
  → generación automática del panel
  → paquete canónico DTC-055
  → mismo validador que un panel suministrado
  → inicialización fundadora intacta o recombinante DTC-055
```

La generación puede ofrecerse como preparación reutilizable antes de ejecutar
la población, pero F6.4 no podrá cerrarse hasta demostrar su integración y la
paridad Marimo exigida por DTC-051. El producto no debe mantener dos contratos
de panel ni dos validadores según el origen suministrado o generado.

El LD y la continuidad haplotípica entre bloques siguen fuera de esta capacidad
y conservan su destino en DTF-044. La recombinación fundadora dentro de cada
bloque permanece exactamente como la define DTC-055.

## 3. Decisiones funcionales ya cerradas para trasladar a fuentes normativas

Estas decisiones pueden condicionar la futura secuencia de trabajo, pero aún no
son pasos implementables porque solo constan en el borrador de conversación:

1. El generador usa `r²` haplotípico faseado y dependiente de la distancia
   física. La generación aleatoria se limita a SNP bialélicos.
2. Cada población puede elegir explícitamente un perfil de `r²` incorporado o
   suministrado. Los perfiles incorporados son ejemplos versionados, nunca
   predeterminados ni inferidos desde el nombre de la población; pueden ser
   aproximaciones humanas sintéticas si se declaran y documentan como tales.
3. Los bloques se muestrean independientemente en la población inicial, sin LD,
   continuidad de donante ni transiciones condicionadas entre bloques.
4. La cantidad se expresa mediante exactamente uno de `block_count`,
   `covered_bases` o `coverage_fraction`, sin modalidad predeterminada. Una
   solicitud imposible falla sin reducir silenciosamente el objetivo.
5. Las longitudes usan exclusivamente una distribución tabulada global para el
   panel; el usuario suministra la tabla o selecciona expresamente un perfil
   incorporado derivado de 1000 Genomas, nunca predeterminado.
6. Cada intervalo de longitud es inclusivo, tiene límites enteros positivos y
   `length_min_bp <= length_max_bp`. Los intervalos están ordenados y no se
   solapan; pueden dejar huecos con probabilidad cero.
7. Los pesos son finitos, estrictamente positivos y suman 1 conforme a DTC-049.
   Primero se sortea el intervalo por peso y después una longitud entera
   uniforme entre ambos límites incluidos.
8. Una tabla vacía, desordenada, solapada, con límites o pesos inválidos o
   incompatible con el espacio elegible se rechaza antes de generar bloques,
   sin corrección implícita.
9. La densidad y la posición relativa de los SNP dentro de cada bloque se
   obtienen de una tabla global de distancias físicas entre SNP consecutivos,
   no de un número de SNP por kilobase. La tabla se suministra o se selecciona
   expresamente de un perfil incorporado, nunca predeterminado, y aplica las
   mismas reglas de intervalos, pesos y suma que la distribución de longitudes.
   El perfil incorporado de distancias entre SNP procede de los mismos datos y
   filtros de 1000 Genomas que sus componentes compatibles; esto no convierte
   el perfil aproximado de `r²` en cuantiles empíricos. El número de sitios es
   un resultado de la longitud del bloque y de los saltos sorteados.
10. Dentro de cada bloque se sortean saltos consecutivos hasta que el siguiente
    deja de caber. Ese último salto se descarta sin truncarlo. La suma de los
    saltos aceptados define el tramo comprendido entre el primer y el último
    SNP; el tramo completo se desplaza uniformemente entre todas las posiciones
    donde cabe dentro del bloque, incluidos ambos extremos. Cada bloque debe
    contener al menos dos SNP. Si no cabe ningún salto, se regenera el bloque
    completo sin acortarlo, reducir sitios ni relajar la distribución, hasta el
    límite `max_block_attempts` aprobado más abajo.
11. Los bloques se colocan uniformemente entre todos los inicios elegibles. El
    intervalo completo debe caber en una única secuencia y no atravesar la
    máscara de exclusión, bases `N`, loci reservados ni bloques ya aceptados.
    Por tanto, el reparto entre secuencias es proporcional al número de
    colocaciones completas que admite cada una, no igual por `sequence_id`. V1
    admite máscara de exclusión, pero no pesos espaciales ni perfiles
    cromosómicos. Las regiones pericentroméricas con baja recombinación y
    bloques mayores quedan fuera de esta versión en DTF-045.
12. La separación mínima entre bloques es cero: un bloque puede terminar en
    `x` y el siguiente comenzar en `x+1`, pero nunca comparten una base. V1 no
    añade una distribución de espacios ni rechaza la adyacencia. La distancia
    interbloque realizada puede registrarse como medida descriptiva, pero no
    condiciona el sorteo ni representa LD. Los modelos espaciales posteriores
    permanecen en DTF-044 y DTF-045.
13. La disposición empieza después de generar la lista completa de longitudes.
    Se colocan de mayor a menor y cada longitud se sortea uniformemente entre
    sus inicios todavía elegibles. Si alguna no cabe, se descarta toda la
    disposición y se repite con la misma lista, sin sustituir ni reducir el
    bloque fallido. Los `block_id` se asignan solo tras lograr una disposición
    válida y según el orden canónico `sequence_id`, inicio y final. Este
    algoritmo debe ser versionado, reproducible y registrar intentos, pero no
    promete uniformidad entre todas las disposiciones globales posibles. La
    recolocación se limita mediante `max_layout_attempts`.
14. `covered_bases` exige un objetivo entero exacto. Para
    `coverage_fraction` se calcula
    `target_bases = floor(eligible_bases * coverage_fraction)`, por lo que la
    cobertura nunca excede la fracción solicitada. Se sortean longitudes que
    quepan en el resto; si una lo supera, se repite únicamente ese sorteo. La
    última longitud debe igualar exactamente el resto y pertenecer a un
    intervalo admitido; en caso contrario se descarta la lista completa. Nunca
    se trunca un bloque ni se acepta una cobertura aproximada. Un objetivo menor
    que la longitud mínima o no expresable como suma de longitudes válidas se
    rechaza antes de colocar. La representabilidad numérica se separa de la
    colocación física. Sea `S` el conjunto de longitudes enteras admitidas que
    caben al menos una vez en un segmento elegible y `m = min(S)`. El preflight
    construye el grafo exacto de residuos `0..m-1`, con transición
    `r -> (r+l) mod m` y coste `l` para cada `l` de `S`; un camino mínimo desde
    cero obtiene `d[r]` y sus predecesores. Un objetivo `T` es representable si y
    solo si `d[T mod m] <= T`. La diferencia es un múltiplo de `m`, que pertenece
    a `S`, y permite reconstruir un testigo exacto. El cálculo usa aritmética
    entera comprobada. El testigo solo demuestra la suma: no sustituye el
    muestreo ponderado ni obliga a usar esa lista. La regeneración de listas se
    limita mediante `max_length_set_attempts`.
15. La política de agotamiento usa tres enteros positivos, configurables y
    versionados: `max_block_attempts = 100`,
    `max_length_set_attempts = 1000` y `max_layout_attempts = 1000`. Antes de
    consumir intentos, el preflight rechaza cobertura mayor que el espacio
    elegible, objetivo menor que la longitud mínima, ausencia de una distancia
    SNP que quepa en alguna longitud admitida y cobertura demostrada no
    representable por el grafo de residuos. Para `block_count`, el máximo
    demostrable es
    `sum(floor(segment_length / min_allowed_length))` sobre los segmentos
    elegibles: excederlo es imposible, pero quedar por debajo no garantiza una
    disposición sorteada concreta. Agotar cualquier límite falla sin relajar
    bloques, cobertura, SNP ni tolerancias.
    Diagnóstico y manifiesto de fallo conservan contador, límite, intentos,
    semilla, objetivo y último rechazo; el usuario puede elevar expresamente el
    límite y repetir. Pasar el preflight numérico no demuestra colocabilidad
    conjunta: la fragmentación o las máscaras pueden agotar layout y se informan
    como tal, no como suma imposible. Si el cálculo no puede completarse por
    límites de recursos, se emite un diagnóstico de recursos y nunca se declara
    imposible.
16. El `r²` se evalúa de forma agregada por población e intervalo de distancia.
    Se incluyen todas las parejas de SNP del mismo bloque y se agrupan las de
    todos los bloques; nunca se incluyen parejas interbloque. Los cuantiles
    realizados se comparan con los solicitados, sin exigir un rango a cada
    pareja ni corregir parejas individuales. Por población e intervalo se
    registran cuantiles solicitados y realizados. Cada fila objetivo contiene
    obligatoriamente `r2_p10`, `r2_p25`, `r2_p50`, `r2_p75` y `r2_p90`; son
    `f64` finitos en `[0,1]`, no decrecientes, y `p50` es la mediana. V1 no añade
    `p0`, `p5`, `p95` ni `p100`. `pair_count` no es un objetivo: en perfiles
    incorporados es evidencia de derivación y en la generación es un recuento
    realizado. Cada fila incluye además `r2_max_abs_error`, un `f64` finito en
    `[0,1]` obligatorio y común a sus cinco cuantiles. Cada cuantil se acepta si
    y solo si `abs(realizado - solicitado) <= r2_max_abs_error`; no existe
    tolerancia relativa o implícita. Un perfil incorporado la deriva y documenta
    y uno suministrado debe declararla. Los intervalos usan límites enteros
    positivos e inclusivos, están ordenados y forman una partición contigua: el
    primero comienza en 1 y cada mínimo posterior es el máximo anterior más 1,
    sin huecos ni solapamientos. El máximo del último intervalo alcanza al
    menos `max_block_length - 1`, aunque puede excederlo. Todos los perfiles
    poblacionales de un mismo panel comparten exactamente esos límites; solo
    pueden diferir en mínimos, tolerancias y cuantiles. Una fila completamente
    fuera de las distancias posibles se informa como no aplicable; cualquier
    fila que interseque, aunque sea parcialmente, el soporte posible permanece
    activa. Todo cuantil empírico usa Hyndman–Fan Type 7: para el vector ordenado
    `x[0..n-1]`, si `n = 1`, `Q = x[0]`; en otro caso,
    `t = (n - 1) * q`, `lo = floor(t)`, `hi = ceil(t)`, `gamma = t - lo` y
    `Q = (1 - gamma) * x[lo] + gamma * x[hi]`. La misma fórmula y los valores
    `q = 0.10`, `0.25`, `0.50`, `0.75` y `0.90` rigen generador, preflight,
    oráculo e informes.
17. Un SNP con `alt_af` 0 o 1 se conserva, pero está fijado en esa
    población. Todas sus parejas se excluyen allí porque el denominador de `r²`
    es cero; nunca se convierten en `r² = 0`. Una misma pareja puede ser
    evaluable en una población y no en otra. Por población e intervalo se
    registran las parejas intrabloque totales, evaluables y excluidas por
    contener uno o dos SNP fijados.
18. Cada fila activa del perfil declara `min_evaluable_pairs`, un entero
    estrictamente positivo específico del intervalo. Un perfil incorporado
    deriva y documenta el mínimo con su receta; un perfil suministrado debe
    declararlo. Las filas completamente fuera del soporte posible se informan
    como no aplicables y no tienen que cumplir ese mínimo; una intersección
    parcial basta para mantener activa la fila. Si una población queda por
    debajo, sus cuantiles de ese intervalo son inválidos y se repite únicamente
    la construcción haplotípica de esa población. Nunca se fusionan intervalos,
    se reduce el mínimo ni se incluyen parejas con SNP fijados.
    `max_haplotype_attempts = 1000` es un entero positivo, configurable y
    versionado y se aplica por población. La geometría común permanece inmutable,
    mientras frecuencias, cuantización, inicialización, cálculo de cuantiles
    Type 7, puntuación y recocido tienen un `population_attempt_index` propio. Un incumplimiento en
    una población o intervalo, incluido un único cuantil fuera de tolerancia,
    rechaza solo su intento y permite conservar internamente las demás
    poblaciones válidas. Agotar el límite de una población falla el panel
    completo, sin resortear las otras ni relajar objetivos, y registra población,
    intento y último incumplimiento. No se publica un panel parcial: todas las
    poblaciones deben superar el validador común antes de la publicación atómica
    y de formar la unión densa final de HAP por bloque.
19. Las frecuencias alélicas marginales de los SNP intrabloque reutilizan los
    espectros de DTC-054. Cada SNP y población obtiene de forma independiente
    su `alt_af` desde el espectro efectivo de esa población; 0 y 1 son valores
    válidos. Estas frecuencias son objetivos marginales que deben conservarse al
    ajustar el LD. Un intento incompatible regenera frecuencias y haplotipos
    solo de esa población, pero mantiene la geometría de bloques, posiciones SNP
    y resultados de las demás. Los perfiles incorporados agrupan espectros
    alélicos y LD cuya compatibilidad y procedencia están documentadas; un usuario puede
    mezclar los suyos, asumiendo el riesgo de incompatibilidad. No se reutiliza
    una tabla `population_snvs.tsv` concreta: esa matriz corresponde a loci
    independientes externos a los bloques y aquí solo se reutiliza el modelo de
    espectros.
20. Cada población declara `panel_haplotype_copies`, un entero estrictamente
    positivo. Para cada bloque se construyen exactamente `C` copias
    haplotípicas, donde `C` es ese valor; las copias idénticas se deduplican y
    su frecuencia es `copy_count / C`. El número de haplotipos distintos es un
    resultado, no un parámetro. El catálogo canónico de cada bloque es la unión
    de los vectores observados en todas las poblaciones, con frecuencia 0
    explícita donde un vector no aparece, conforme al contrato denso de
    DTC-055. Un perfil incorporado conserva el tamaño haploide usado en su
    derivación o parametrización. Cada `alt_af` solicitada se cuantiza mediante
    redondeo estocástico reproducible: para `x = alt_af * C`, el recuento es
    `floor(x)` más una Bernoulli cuya probabilidad es la parte fraccionaria de
    `x`. Los enteros se conservan exactamente, por lo que 0 produce 0 y 1
    produce `C`. El sorteo usa una subsemilla propia por población, bloque y
    SNP. La frecuencia efectiva `count / C`, registrada junto con la solicitada
    y el recuento, es la marginal exacta que debe conservar el constructor. Una
    variante rara puede desaparecer y no se rescata.
21. El constructor representa cada población y bloque mediante una matriz
    binaria de `C` filas haplotípicas y una columna por SNP. Cada columna se
    inicializa reproduciblemente con su recuento ALT exacto. Un movimiento
    intercambia 0 y 1 entre dos filas de una columna y, por construcción,
    conserva todas las marginales. Para cada población, intervalo activo y
    cuantil `q` en p10, p25, p50, p75 y p90, la función objetivo calcula
    `delta = abs(realized - target)`,
    `excess = max(0, delta - r2_max_abs_error)` y el término `excess^2`. El
    `score` de cada población es la suma `f64` de esos términos sobre todos sus
    intervalos activos, en orden canónico estable y mediante suma
    compensada de Neumaier, sin reducción paralela; los intervalos no aplicables
    no participan. Si el score inicial es 0, el intento termina con éxito sin
    movimientos. En otro caso, `T0 = max(score_inicial, 1e-12)`. Cada movimiento
    propuesto usa la temperatura vigente: una mejora o empate con
    `delta_score <= 0` se acepta siempre; un empeoramiento se acepta si y solo si
    `U < exp(-delta_score / T)`, donde `U` es uniforme en `[0,1)` y procede del
    ámbito RNG del constructor. Tras cada bloque exacto de 1000 movimientos
    propuestos, aceptados o no, se aplica `T *= 0.95`; la temperatura no se
    recalcula desde scores posteriores. Los mínimos evaluables
    se comprueban antes de optimizar; el intento se acepta si cumple esa
    población y, al final, se deduplican sus filas. Como `r²` no conserva
    el signo de la asociación, acoplamiento y repulsión equivalentes son válidos
    y la semilla decide el resultado; el control futuro mediante `r` firmado se
    reserva a DTF-046. La aceptación final exige `score == 0.0`, todos los
    mínimos y la comprobación individual de cada celda; el score guía la
    búsqueda, pero no sustituye esa validación. Se registran el mejor score y
    sus mejores cuantiles. El proceso evalúa y se detiene inmediatamente cuando
    score, mínimos y celdas son válidos. Se registran `T0`, factor 0.95, longitud
    de etapa 1000, movimientos propuestos, aceptados y empeoramientos aceptados,
    además de la temperatura final; score, temperatura y mejor resultado son
    poblacionales y el calendario está versionado. El tamaño del problema se
    define para cada población como `population_haplotype_cells`, la suma en
    orden canónico de `panel_haplotype_copies * snp_count` sobre sus bloques. Con
    aritmética entera comprobada y sin overflow,
    `max_annealing_moves = min(10_000_000, max(100_000,
    100 * population_haplotype_cells))`. La fórmula es determinista, versionada
    y no configurable en esta primera versión, y cuenta movimientos propuestos. Al
    agotarla sin un estado válido se rechaza el intento de esa población y el
    siguiente usa su propio `population_attempt_index` y subsemilla, hasta
    `max_haplotype_attempts = 1000` para esa población. Agotar sus intentos
    produce `sin solución encontrada dentro del límite`, nunca imposibilidad
    matemática. Se registran las celdas, fórmula y
    versión, límite, propuestas y mejor score y cuantiles.
22. Antes del recocido, el preflight examina cada pareja bialélica no fijada.
    Para `C` copias y recuentos ALT `a` y `b`, enumera el recuento conjunto `h`
    sobre todos los enteros de `[max(0, a + b - C), min(a, b)]` y calcula su
    dominio exacto de `r²`, con `p_A = a/C`, `p_B = b/C`, `p_AB = h/C` y
    `r² = (p_AB - p_A*p_B)^2 / (p_A*(1-p_A)*p_B*(1-p_B))`. Por población e
    intervalo construye los vectores de mínimos y máximos por pareja, los ordena
    por separado y aplica Type 7 para obtener, por cuantil, un intervalo relajado
    conservador de valores alcanzables, ignorando compatibilidades multilocus.
    Solo falla como `matemáticamente imposible` si el intervalo solicitado
    `[target - r2_max_abs_error, target + r2_max_abs_error]` no interseca el
    alcanzable de algún cuantil o si las marginales determinan menos parejas
    evaluables que `min_evaluable_pairs`. Superar esta prueba no demuestra
    factibilidad. Una incompatibilidad multilocus puede agotar
    el `max_haplotype_attempts` de esa población como `sin solución encontrada
    dentro del límite`, nunca como prueba matemática. Ningún fallo admite
    `--force`. El diagnóstico
    conserva población, intervalo, cuantil, rango alcanzable, rango solicitado
    con tolerancia, `C`, marginales relevantes y motivo.
23. El azar usa seis ámbitos físicos exactos con `derivation_version = 2`:
    `ld_geometry` para longitudes, disposición, saltos SNP y desplazamiento;
    `ld_site_alleles` para ALT; `ld_population_frequencies` para sorteo de
    `alt_af` y redondeo; `ld_haplotype_initialization` para repartir recuentos en
    C filas; `ld_annealing_proposals` para elegir columnas, filas e intercambios;
    y `ld_annealing_acceptance` solo para el uniforme de movimientos
    empeorantes. La separación entre propuestas y aceptación evita que consumir
    o no un uniforme desplace propuestas posteriores.

    Cada derivación incorpora únicamente los IDs estables aplicables entre
    `scenario_id`, `reference_id`, `population_id`, `block_id`, `site_id` y
    `population_attempt_index`; prohíbe ruta, reloj, hilo, *scheduler*, orden de
    iteración y posición accidental. Los ámbitos poblacionales incorporan
    `population_id` y su índice de intento. Cada reintento conserva los sorteos
    anteriores. Cambiar, añadir o reordenar una población no altera la geometría
    ni las frecuencias, haplotipos o recocido de las demás; reordenar otras
    entidades identificadas o salidas tampoco cambia elementos no relacionados.
    Cambiar cantidad, cobertura o distribución sí reconstruye la geometría.

    Los catálogos v1 y v2 permanecen intactos, legibles y reproducibles; ninguna
    semilla antigua se reinterpreta. El futuro catálogo exacto será
    `catalogs/subseed-scopes-v3.tsv`, con `catalog_version = 3` y
    `derivation_version = 2`. Incluirá los ámbitos anteriores necesarios y los
    seis nuevos de LD, con vectores normativos propios. La generación automática
    seleccionará explícitamente catálogo 3 y derivación 2 sin mezclar versiones
    silenciosamente. Configuración efectiva y manifiesto registrarán catálogo,
    derivación, contexto, ámbito, componentes y subsemillas conforme a
    DTC-037/DTC-038. El catálogo real no se crea en este preplan; su publicación
    requiere cerrar las puertas normativas y autorizar la implementación.

24. La salida es exactamente el paquete canónico de DTC-055 y atraviesa el mismo
    validador que un panel suministrado.
25. Todo panel generado se valida de nuevo leyendo desde cero sus cinco
    artefactos publicados, sin matrices, cachés ni puntuaciones del optimizador.
    El oráculo aplica DTC-055 y reconstruye recuentos y frecuencias ALT, parejas
    intrabloque, exclusiones por SNP fijados, `r²`, cuantiles, mínimos y
    tolerancias. En CI y aceptación, fixtures faseados pequeños con `r²` 0, 1 e
    intermedio, SNP fijados y varias poblaciones se contrastan además con
    PLINK 2. PLINK 2 es solo un validador externo de pruebas, no una dependencia
    de ejecución. Para `human_ld_ragsdale_approx_v1`, fixtures independientes
    comprueban además los cuantiles Beta con `kappa = 8.0` para medias
    intermedias y las distribuciones degeneradas con `mu` 0 y 1. Si la
    recomputación desde los artefactos no reproduce las marginales y el perfil
    de LD, el panel es inválido. Fixtures pequeños enumeran exhaustivamente las
    realizaciones compatibles y comprueban Type 7 y que el preflight no produzca
    falsos `matemáticamente imposible`. Un oráculo independiente por enumeración
    o programación dinámica exhaustiva comprueba en fixtures pequeños el grafo
    de residuos, intervalos degenerados y con huecos, la recomposición exacta del
    testigo y la frontera entre suma imposible, layout agotado y falta de
    recursos del preflight. Otros oráculos cambian el perfil, reordenan o añaden
    una población y exigen bytes idénticos en artefactos poblacionales no
    relacionados cuando estén físicamente separados o, en su defecto, filas y
    proyecciones canónicas idénticas. Un fallo poblacional forzado no publica un
    panel parcial ni resortea poblaciones ya válidas.
26. La generación publica cuatro informes TSV resumidos, externos al directorio
    del paquete canónico DTC-055 y serializados en orden canónico:
    `block_generation_summary.tsv`, con una fila por bloque y los campos
    `block_id`, `sequence_id`, `start`, `end`, `length`, `snp_count`,
    `first_snp`, `last_snp` y `layout_attempt`;
    `allele_frequency_realization.tsv`, con una fila por población, bloque y SNP
    y los campos `population_id`, `block_id`, `site_id`, `requested_alt_af`,
    `alt_count`, `panel_haplotype_copies`, `effective_alt_af` y `fixed`;
    `ld_profile_realization.tsv`, con una fila por población e intervalo y
    `population_id`, límites de distancia, estado `evaluated|not_applicable`,
    parejas totales y evaluables, recuentos separados de exclusión por uno y por
    dos SNP fijados, `min_evaluable_pairs`, `r2_max_abs_error`, los cinco
    cuantiles objetivo y realizados y `passed`; y
    `generation_attempts.tsv`, con `stage`, `population_id`, `block_id`,
    `attempt_limit`, `attempts_used`, `final_status` y `last_rejection`. El
    manifiesto de la ejecución generadora los enlaza por clase, ruta, tamaño y
    SHA-256. Su contenido permite reconstruir la aceptación o el rechazo. Los
    cinco ficheros DTC-055 no cambian y no se publica por defecto la tabla de
    todas las parejas de `r²`; esa ampliación se reserva a DTF-047.
27. V1 puede incorporar, bajo el identificador definitivo
    `human_ld_ragsdale_approx_v1`, una aproximación humana transparente
    inspirada en las curvas de caída de LD de Ragsdale y Gravel (2020) para
    cinco poblaciones de 1000 Genomas. Usa la curva `r_div^2` del panel 2C como
    tendencia media `mu` por distancia, no Rogers-Huff ni `sigma_D^2`. Por
    población e intervalo genera una distribución Beta acotada en `[0,1]`: si
    `0 < mu < 1`, usa `alpha = mu * 8.0` y
    `beta = (1 - mu) * 8.0`, y calcula `p10`, `p25`, `p50`, `p75` y `p90`
    mediante la inversa de la CDF en 0.10, 0.25, 0.50, 0.75 y 0.90. Para
    `mu = 0` o `mu = 1`, la distribución es degenerada y los cinco cuantiles
    valen `mu`. `kappa = 8.0` es un único `f64` finito y positivo, común a todas
    las poblaciones e intervalos y parte de la versión del builtin; no es
    configurable, inferido ni ajustado durante la ejecución. Un perfil
    suministrado puede declarar directamente sus cuantiles y no está obligado a
    usar Beta.
    El recurso conserva fuente y cita, método de digitalización o tabla
    reconstruida, ajuste, intervalos, parámetros y hashes. Se declara siempre
    aproximado y sintético, nunca como cuantiles empíricos, y requiere selección
    explícita. Esta vía no exige procesar de inicio los 2504 genomas y no impide
    añadir después un perfil empírico exacto.
28. Todo perfil reutilizable del generador es un paquete versionado con
    exactamente cuatro ficheros obligatorios y ningún archivo adicional:
    `block_lengths.tsv`, para intervalos y pesos globales de longitud;
    `snp_gaps.tsv`, para intervalos y pesos globales de saltos;
    `ld_profile.tsv`, para identidades poblacionales del recurso, intervalos y
    cuantiles de LD; y `profile.metadata.json`. El JSON declara clase, versión,
    `profile_id`, carácter `builtin|supplied` y `approximate|empirical`,
    procedencia, receta y parámetros, y clase, ruta, tamaño y SHA-256 de cada
    TSV, sin duplicar sus tablas. Las rutas son relativas al paquete. La máscara
    de exclusión queda fuera porque pertenece a la referencia concreta de la
    ejecución y se declara y verifica por SHA-256 de forma separada. Los
    paquetes incorporados y suministrados atraviesan el mismo validador.
29. `block_lengths.tsv` usa la cabecera exacta
    `length_min_bp<TAB>length_max_bp<TAB>weight`. Se serializa en UTF-8 sin BOM,
    con tabuladores, finales de línea LF y LF final; no admite comillas, valores
    vacíos, columnas adicionales ni comentarios. `length_min_bp` y
    `length_max_bp` son enteros decimales positivos canónicos y cumplen
    `length_min_bp <= length_max_bp`; `weight` es un `f64` decimal finito,
    estrictamente positivo y con representación de ida y vuelta exacta. Las
    filas se ordenan ascendentemente por `(length_min_bp, length_max_bp)`, con
    una sola fila por intervalo, sin duplicados ni solapamientos. Los intervalos
    son inclusivos y pueden dejar huecos. La suma de pesos se valida conforme a
    DTC-049; un incumplimiento se rechaza sin normalizar ni corregir la tabla.
30. `snp_gaps.tsv` usa la cabecera exacta
    `gap_min_bp<TAB>gap_max_bp<TAB>weight` y las mismas reglas físicas y
    semánticas: UTF-8 sin BOM, tabuladores, LF y LF final; sin comillas, vacíos,
    columnas adicionales ni comentarios; límites enteros decimales positivos
    canónicos con `gap_min_bp <= gap_max_bp`; y `weight` como `f64` decimal
    finito, estrictamente positivo y con ida y vuelta exacta. Las filas se
    ordenan ascendentemente por `(gap_min_bp, gap_max_bp)` y contienen intervalos
    inclusivos únicos, sin duplicados ni solapamientos y con huecos permitidos.
    La suma se valida conforme a DTC-049 y se rechaza sin normalización ni
    corrección si incumple. El ejemplo con intervalos 1–10, 11–100, 101–500,
    501–1000 y 1001–5000 y pesos 0.05, 0.25, 0.40, 0.20 y 0.10 es ilustrativo,
    no un valor predeterminado.
31. `ld_profile.tsv` usa la cabecera exacta
    `population_profile_id<TAB>distance_min_bp<TAB>distance_max_bp<TAB>`
    `min_evaluable_pairs<TAB>r2_max_abs_error<TAB>r2_p10<TAB>r2_p25<TAB>`
    `r2_p50<TAB>r2_p75<TAB>r2_p90`. Aplica la misma serialización canónica:
    UTF-8 sin BOM, tabuladores, LF y LF final, sin comillas, vacíos, columnas
    adicionales ni comentarios. `population_profile_id` cumple DTC-039 y es
    una identidad única dentro del paquete; los límites son enteros decimales
    positivos e inclusivos; `min_evaluable_pairs` es un entero estrictamente
    positivo; y el error y los cinco cuantiles son `f64` decimales finitos en
    `[0,1]`, con ida y vuelta exacta. Los cuantiles son no decrecientes. La clave
    `(population_profile_id, distance_min_bp, distance_max_bp)` es única y las
    filas se ordenan por `population_profile_id` ASCII y después por ambos
    límites. Cada población contiene exactamente una fila por intervalo; todas
    comparten fronteras idénticas y contiguas desde 1, y el último máximo
    alcanza al menos la longitud máxima de bloque menos 1. El mínimo permanece
    declarado aunque una fila resulte no aplicable en una ejecución concreta,
    caso en el que no se exige. Los valores CEU/YRI de ejemplo son ilustrativos.
32. `profile.metadata.json` es un objeto cerrado, sin campos desconocidos ni
    tablas duplicadas, serializado en UTF-8 sin BOM y con LF final. Declara
    obligatoriamente `profile_kind = "ld_generation_profile"`,
    `schema_version = 1`, un `profile_id` conforme a DTC-039,
    `package_origin` en `builtin|supplied`, `derivation_kind` en
    `approximate|empirical`, `statistic = "phased_haplotype_r2"`,
    `variant_class = "biallelic_snp"`, `population_profile_ids`,
    `quantile_model`, `provenance` y `artifacts`. `population_profile_ids` es un
    array no vacío de IDs DTC-039 únicos, ordenados por ASCII, y coincide
    exactamente con el conjunto de `ld_profile.tsv`. `quantile_model` es
    discriminado: el builtin de Ragsdale usa
    `kind = "beta_from_mean_curve"` y `kappa = 8.0`; los perfiles con cuantiles
    directos usan `kind = "direct_quantiles"`. Una derivación `approximate`
    exige `citation`, `source_detail` y `derivation_recipe` no vacíos, además de
    su modelo. Una derivación `empirical` exige procedencia y receta explícitas,
    con `source_dataset`, `source_release`, `reference_build`, `samples`,
    `regions` y `filters` no vacíos; nada se infiere. `artifacts` contiene
    exactamente tres elementos, uno por clase `block_lengths`, `snp_gaps` y
    `ld_profile`; sus rutas relativas simples son exactamente los nombres
    canónicos, `size` es un entero mayor o igual que 0 y `sha256` contiene 64
    dígitos hexadecimales minúsculos y se recalcula al validar.
33. La máscara de usuario es un fichero opcional `excluded_regions.tsv`, externo
    al paquete de perfil y ligado a la referencia de la ejecución. Usa la
    cabecera exacta `region_id<TAB>sequence_id<TAB>start<TAB>end` y la
    serialización canónica UTF-8 sin BOM, tabuladores, LF y LF final, sin
    comillas, vacíos, columnas adicionales ni comentarios. `region_id` cumple
    DTC-039 y es único; `sequence_id` existe en el FASTA; `start` y `end` son
    enteros decimales positivos en base 1 y forman un intervalo cerrado con
    `start <= end`, enteramente contenido en la secuencia. Las filas se ordenan
    por `sequence_id` ASCII, `start`, `end` y `region_id`. Las regiones no se
    solapan, aunque pueden ser adyacentes. Su SHA-256 se calcula y registra junto
    con `reference_id` y `reference_sha256` en configuración efectiva y
    manifiesto. Si falta el fichero no existe máscara de usuario; siguen
    excluyéndose automáticamente `N`, loci reservados y bloques ya aceptados.
    Una máscara nunca puede reactivar esos espacios.
34. La configuración común reside exactamente en
    `[initial_population.ld_blocks]`, con `source = "generated"`; exactamente
    uno de `profile = "builtin:<profile_id>"` o `profile_path` relativo al
    paquete; exactamente uno de `block_count`, `covered_bases` o
    `coverage_fraction`; y `excluded_regions_path` opcional. Cada población
    simulada alcanzada declara obligatoriamente un `population_profile_id`
    existente en `ld_profile.tsv` bajo
    `[initial_population.ld_blocks.populations.<population_id>]`. `profile` y
    `profile_path` son obligatorios y mutuamente excluyentes. Las rutas se
    resuelven respecto del documento propietario conforme a DTC-101 cuando se
    cierre. Las sustituciones siguen DTC-033 y solo pueden cambiar la asignación
    `population_profile_id` de un ámbito poblacional específico; no modifican
    parcialmente `block_lengths.tsv`, `snp_gaps.tsv`, `ld_profile.tsv` ni los
    parámetros de un builtin. Para cambiar tablas o builtin se aporta un paquete
    supplied completo. La geometría común se resuelve una vez y los remapeos no
    cambian bloques ni posiciones. Un `source` distinto o desconocido, una
    población sin mapeo o campos de modos incompatibles se rechazan antes de
    generar.

Los valores de los ejemplos TSV, JSON y TOML del borrador siguen siendo
ilustrativos. Este preplan conserva los nombres funcionales ya cerrados, pero no
aprueba todavía el esquema físico completo ni otros campos provisionales.

## 4. Puertas normativas y de diseño

### 4.1. Contradicción que bloquea la planificación ejecutable

La especificación vigente define en REQ-EVO-74 un «panel suministrado» y su
sección 9 afirma que en v1 se suministran las cuatro tablas o se generan
externamente. La generación integrada cambia entradas, comportamiento,
resultados y errores observables. Por ello no puede incorporarse únicamente
mediante un plan o una decisión técnica.

Antes de redactar pasos ejecutables deben completarse, en este orden:

1. **PENDIENTE DE REVISIÓN — especificación:** añadir REQ propios para la
   generación, ampliar los requisitos de configuración, validación, manifiesto
   y Marimo que correspondan y definir ACE independientes. Debe retirarse o
   acotarse la exclusión de generación interna sin debilitar el contrato de
   panel suministrado.
2. **PENDIENTE DE REVISIÓN — trazabilidad:** asignar cada REQ/ACE a F6.4 y a
   una evidencia prevista, sin presentar DTC-055 como fuente suficiente de la
   generación.
3. **PENDIENTE DE REVISIÓN — arquitectura:** asignar responsabilidades para
   perfiles, generación, validación, publicación, aleatoriedad y Marimo;
   definir contratos, flujo, invariantes y errores sin duplicar el productor y
   validador canónicos, incorporando el método de residuos acordado, su testigo y
   la separación entre representabilidad, colocación y límites de recursos.
4. **PENDIENTE DE REVISIÓN — decisión técnica:** abrir y resolver una DTPD para
   formatos de perfiles, algoritmos de construcción y cálculo de cuantiles,
   derivación y versionado del perfil aproximado, RNG y oráculos. DTC-055 solo
   fija el recurso de salida y su consumo.
5. **PENDIENTE DE REVISIÓN — hoja de ruta:** confirmar F6.4 y dividir, si la
   arquitectura lo justifica, preparación automática, admisión común e
   inicialización fundadora como resultados revisables.
6. **PENDIENTE DE REVISIÓN — plan F6:** redactar el plan detallado únicamente
   cuando sus fases previas y todas las puertas anteriores estén cerradas.

## 5. Estado de las decisiones funcionales

`human_ld_ragsdale_approx_v1` queda ratificado como `profile_id` físico
definitivo para TOML, metadatos, manifiestos y ejemplos. No quedan decisiones
funcionales o técnicas materiales abiertas.

Siguen pendientes, por separado, las puertas normativas y de diseño de §4:
estas decisiones continúan en un borrador y no autorizan un plan ejecutable ni
sustituyen REQ, ACE, arquitectura, trazabilidad o DTC.

## 6. Mapa preliminar de dependencias y riesgos

| Resultado futuro | Depende de | Riesgo silencioso que deberá protegerse |
|---|---|---|
| Perfiles admitidos e incorporados | Traslado normativo del contrato cerrado; DTC-040/DTC-049 | Aplicar un perfil implícito, mezclar población, referencia, filtros o versión, o presentar una aproximación sintética como cuantiles empíricos |
| Paquete de perfil reutilizable | Contratos de los cuatro ficheros, rutas relativas, hashes y validador común aprobados | Admitir archivos extra o campos desconocidos, duplicar tablas en JSON, usar rutas no canónicas, no recalcular hashes o validar distinto un builtin |
| Procedencia del perfil | Metadatos cerrados para `approximate|empirical` y modelos Beta/directo | Inferir procedencia, omitir muestras/filtros empíricos, aceptar poblaciones divergentes del TSV o etiquetar como empírico un perfil aproximado |
| Configuración TOML | Sección, exclusividades, mapeos y precedencia poblacional aprobados; resolución de rutas remite a DTC-101 | Combinar builtin y ruta, omitir una población, parchear tablas parcialmente, aceptar otro source o reconstruir geometría al remapear una población |
| Perfil humano aproximado | Tendencia `r_div^2`, Beta, `kappa = 8.0`, fuente y trazabilidad aprobadas; `profile_id` físico por ratificar | Usar otra estadística del artículo, configurar o ajustar `kappa` en ejecución, cambiar el nombre solo en parte de los contratos, ocultar digitalización/ajuste, omitir hashes o afirmar reproducción exacta de 1000 Genomas |
| Cantidad y longitudes de bloques | Decisiones, límites y contrato físico de `block_lengths.tsv` cerrados en §3 | Cobertura o recuento reducido, sesgo por intervalos de distinto ancho, borde inclusivo incorrecto, orden no canónico o normalización silenciosa |
| Cobertura solicitada | Objetivo exacto, grafo de residuos, testigo y separación de la colocación aprobados | Exceder la fracción, truncar el último bloque, aceptar aproximación, desbordar la aritmética, convertir falta de recursos o layout agotado en falso imposible, o usar el testigo en lugar del muestreo ponderado |
| Colocación entre secuencias | Inicios elegibles, máscara canónica, reparto, algoritmo y límite aprobados | Atravesar `N`, máscara o loci reservados, repartir por igual entre secuencias, sustituir solo el bloque fallido o introducir pesos espaciales de DTF-045 |
| Máscara de exclusión | `excluded_regions.tsv`, coordenadas, orden, hash y ausencia aprobados | Interpretar coordenadas como base 0, admitir solapes o secuencias inexistentes, omitir el hash o usar la máscara para reactivar espacios excluidos |
| Identidad de bloques | Disposición válida y orden canónico aprobados; contrato físico de IDs remite a DTC-039 | Asignar identidad antes de colocar, conservar prioridad por longitud o derivar el ID del orden aleatorio |
| Separación interbloque | Adyacencia y no solapamiento aprobados | Rechazar bloques adyacentes, introducir una distribución de espacios o interpretar la distancia realizada como LD |
| Distancias y número realizado de sitios | Contrato de `snp_gaps.tsv`, descarte, desplazamiento, mínimo y límite aprobados | Truncar el último salto, sesgar el desplazamiento, conservar un bloque con menos de dos SNP, normalizar pesos o tratar el ejemplo como predeterminado |
| Agotamiento | Tres contadores estructurales y uno haplotípico, límites y fallo aprobados | Reiniciar el contador incorrecto, relajar objetivos, perder semilla, objetivo y último rechazo, o confundir layout agotado y límite de recursos con imposibilidad numérica |
| Reproducibilidad RNG | Seis ámbitos físicos de derivación 2, intentos por población y futuro catálogo v3 aprobados | Acoplar propuestas con aceptación, usar componentes accidentales, mezclar versiones, sobrescribir catálogos publicados, reinterpretar semillas o permitir que una población altere resultados no relacionados |
| Frecuencias marginales intrabloque | Espectros DTC-054, redondeo estocástico y movimientos conservadores aprobados | Reutilizar una matriz concreta de loci independientes, correlacionar poblaciones, alterar `count/C` al ajustar LD, rescatar variantes raras o regenerar la geometría |
| Sitios y haplotipos por población | Perfil de LD, tabla de saltos, frecuencias marginales y constructor por intercambios aprobados | `r²` calculado sobre pares distintos, marginales incompatibles, fase incorrecta o aceptar solo algunas poblaciones |
| Tamaño y catálogo haplotípicos | `panel_haplotype_copies`, cuantización, deduplicación, frecuencia por recuento y unión densa aprobados | Construir un número distinto de copias por bloque, tratar el número de HAP distintos como entrada, perder frecuencias 0 o emitir catálogos distintos por población |
| Puntuación haplotípica | Exceso cuadrático, orden, Neumaier, score poblacional y validación por celda aprobados | Sumar error dentro de tolerancia, mezclar poblaciones, usar reducción paralela, incluir bins no aplicables, aceptar solo por redondeo del score u omitir los mejores cuantiles |
| Recocido haplotípico | Intercambios, aceptación, calendario y límite por celdas aprobados | Enfriar por aceptados, contar movimientos aceptados, desbordar la fórmula, configurar el límite, no renovar subsemilla o confundir agotamiento con imposibilidad |
| Viabilidad haplotípica | Dominios discretos, cotas Type 7, clasificaciones, diagnóstico y fallo sin `--force` aprobados | Declarar imposible un caso no demostrado, introducir compatibilidad multilocus en la cota relajada, confundir agotamiento con prueba matemática u omitir rangos del diagnóstico |
| Cuantiles empíricos | Hyndman–Fan Type 7 común y fixtures exhaustivos aprobados | Usar otra interpolación entre generador y oráculo, desordenar vectores de cotas o aceptar falsos imposibles |
| Evaluación agregada de LD | Todas las parejas intrabloque, población×intervalo, cinco cuantiles, tolerancia, límites y salidas solicitadas/realizadas aprobados | Incluir parejas interbloque, validar cada pareja, ponderar bloques en lugar de parejas, omitir pares o corregir SNP individuales |
| Partición de distancias de LD | Límites inclusivos, contiguos, soporte, fronteras comunes y serialización aprobados | Admitir huecos o solapamientos, usar fronteras distintas por población, omitir una fila parcialmente aplicable o exigir mínimo a una fila fuera del soporte |
| Perfil objetivo de `r²` | Cinco cuantiles, tolerancia y contrato físico de `ld_profile.tsv` aprobados | Aceptar cuantiles desordenados, fronteras distintas por población, omitir bins, interpretar `pair_count` como objetivo o usar tolerancia relativa/implícita |
| SNP fijados en LD | Exclusión y recuentos por población×intervalo aprobados | Convertir denominador cero en `r²=0`, excluir globalmente una pareja evaluable en otra población o perder el motivo de exclusión |
| Soporte mínimo de LD | `min_evaluable_pairs`, regeneración poblacional y límite por población aprobados | Fusionar intervalos, rebajar el mínimo, incluir fijados, conservar cuantiles con soporte insuficiente o resortear poblaciones no causantes |
| Paquete generado | DTC-055, DTC-039, DTC-040 y publicación común | Diferencias entre panel suministrado y generado, IDs/hashes/orden divergentes |
| Oráculo desde artefactos | Cinco artefactos publicados, DTC-055 y contraste de fixtures con PLINK 2 aprobados | Validar cachés del optimizador, ocultar una serialización incorrecta, depender de PLINK en ejecución o aceptar marginales/LD no reproducibles |
| Informes de generación | Cuatro TSV resumidos y enlaces de manifiesto aprobados; permanecen fuera del paquete DTC-055 | Mezclar evidencias con el paquete canónico, omitir un rechazo, impedir reconstruir la decisión, publicar todos los pares por defecto o no verificar el SHA-256 |
| Inicialización fundadora | Panel validado y DTC-055 | Usar el generador para alterar intacto/recombinante, F o frecuencias durante el consumo |
| Autoría Marimo | DTC-051 y contrato Rust futuro | Duplicar semántica o dejar opciones disponibles solo en TOML |

## 7. Condición para sustituir este preplan

Este documento debe sustituirse por un plan ejecutable de F6.4 cuando:

- todas las decisiones funcionales de §3 estén trasladadas a fuentes normativas
  o técnicas adecuadas;
- la especificación ya no contradiga la generación integrada;
- arquitectura y trazabilidad identifiquen responsables y evidencias;
- cada paso pueda citar REQ/ACE y arquitectura concretos;
- los riesgos de LD, fase, frecuencias, coordenadas y reproducibilidad tengan
  oráculos independientes definidos;
- cada hito termine en una pausa de revisión y Marimo forme parte del cierre.
