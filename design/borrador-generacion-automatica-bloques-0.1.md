# Borrador de generación automática de bloques de LD 0.1

- **Estado:** borrador de conversación; no normativo.
- **Fecha:** 13 de agosto de 2026.
- **Propósito:** conservar las decisiones y propuestas discutidas antes de
  incorporarlas a las especificaciones, arquitectura y planes.
- **Contrato de salida previsto:** el paquete canónico de cinco ficheros fijado
  por DTC-055, sometido al mismo validador que un panel suministrado.

## 1. Acuerdos alcanzados

1. La generación automática utilizará el `r²` haplotípico faseado como medida
   de LD. La generación aleatoria de esta versión se limitará a SNP bialélicos;
   los multialélicos aleatorios continúan aplazados.
2. El `r²` objetivo no será constante por definición. Su distribución dependerá
   de la distancia física entre los SNP.
3. Cada población podrá utilizar un perfil de LD diferente. El usuario podrá
   seleccionar un perfil incorporado o aportar uno propio.
4. El primer perfil humano incorporado podrá ser una aproximación sintética
   versionada inspirada en curvas publicadas de 1000 Genomas. Será un ejemplo
   seleccionable y nunca se aplicará por defecto ni se inferirá a partir del
   nombre de una población simulada. Su identidad y sus metadatos indicarán
   `approx` o `synthetic`; no podrá presentarse como una reconstrucción de
   cuantiles empíricos.
5. Un mismo perfil podrá reutilizarse en varias poblaciones únicamente mediante
   una selección explícita.
6. En la población inicial, los bloques se muestrearán independientemente. No
   habrá donante común, transición condicionada ni LD objetivo entre bloques.
   El LD y la continuidad haplotípica interbloque quedan en DTF-044.
7. La independencia entre bloques no modifica la recombinación fundadora dentro
   de cada bloque aprobada en DTC-055. Durante cruces y evolución posteriores,
   la recombinación operará sobre las copias cromosómicas ya construidas y podrá
   romper los bloques conforme a la tasa o al mapa aplicable.
8. La cantidad de bloques se declarará mediante exactamente una de tres
   modalidades mutuamente excluyentes: `block_count`, `covered_bases` o
   `coverage_fraction`. Ninguna será predeterminada. `block_count` exige un
   número exacto y deja que la cobertura resulte de las longitudes generadas;
   `covered_bases` y `coverage_fraction` fijan una cobertura objetivo y dejan
   que el número de bloques sea un resultado. La fracción se calcula sobre las
   bases elegibles. Solo son elegibles los intervalos en los que cabe un bloque
   completo sin salir de la secuencia, atravesar una región excluida, solaparse
   con otro bloque ni entrar en conflicto con loci reservados. Una solicitud
   imposible debe rechazarse sin reducir silenciosamente cantidad o cobertura.
9. Las longitudes se generarán exclusivamente desde una distribución tabulada
   por intervalos y pesos. No habrá modos alternativos de longitud fija ni
   intervalo uniforme global. La distribución será única para el panel porque
   los límites de los bloques son comunes a todas sus poblaciones; no debe
   confundirse con los perfiles de `r²`, que sí pueden variar por población.
   El usuario podrá suministrar la tabla o seleccionar un perfil versionado
   incorporado. El programa incluirá al menos un perfil derivado de datos de
   1000 Genomas como ejemplo seleccionable, nunca predeterminado.
10. La tabla de longitudes debe contener intervalos inclusivos con límites
    enteros positivos y `length_min_bp <= length_max_bp`, ordenados de menor a
    mayor y sin solapamientos. Se permiten huecos, que representan longitudes
    con probabilidad cero. Todos los pesos deben ser finitos y estrictamente
    positivos y su suma debe ser 1 conforme a DTC-049. Primero se sortea un
    intervalo por sus pesos y después una longitud entera uniforme entre ambos
    límites, incluidos. Una tabla vacía, desordenada, solapada, con límites o
    pesos inválidos o incompatible con el espacio elegible debe rechazarse antes
    de generar bloques, sin normalización ni corrección implícita.
11. La densidad y la posición relativa de los SNP dentro de un bloque se
    modelarán mediante una distribución tabulada de distancias físicas entre
    SNP consecutivos, no mediante un número fijo de SNP por kilobase. La tabla
    será única para el panel, suministrada o seleccionada de un perfil
    incorporado, y nunca tendrá un valor predeterminado. El perfil humano
    aproximado agrupará tablas de saltos, longitudes, frecuencias y LD bajo una
    receta sintética única, versionada y documentada. El número de sitios de
    cada bloque será un resultado de su longitud y de las distancias sorteadas.
12. Para cada bloque se sortearán distancias consecutivas hasta que el siguiente
    salto ya no quepa. Ese último salto se descartará sin acortarlo. La suma de
    los saltos aceptados define el tramo entre el primer y el último SNP; el
    tramo completo se desplazará después uniformemente entre todas las
    posiciones donde cabe dentro del bloque, incluidos ambos extremos. Así el
    primer SNP no queda fijado al inicio. Todo bloque debe contener al menos dos
    SNP conforme a DTC-055. Si no cabe ninguna distancia válida, se repetirá la
    generación completa de ese bloque, hasta el límite configurable
    `max_block_attempts = 100`. Este límite no podrá sustituirse por acortar un
    salto, reducir sitios o relajar la tabla.
13. Los bloques se colocarán uniformemente entre todos los inicios elegibles de
    la referencia. Un inicio solo será elegible si el intervalo completo cabe en
    una secuencia y no atraviesa regiones excluidas, bases `N`, loci reservados
    ni bloques ya aceptados. Por ello, una secuencia recibirá bloques en
    proporción al número de colocaciones completas que admita, no un número
    igual por `sequence_id`. Esta versión admite una máscara de exclusión, pero
    no pesos espaciales ni perfiles cromosómicos de colocación.
14. La creación explícita de regiones pericentroméricas con menor recombinación
    y bloques mayores queda aplazada en DTF-045. Ese modelo futuro deberá ligar
    anotación espacial, distribución de longitudes y mapa de recombinación sin
    modificar silenciosamente la colocación uniforme de esta versión.
15. La separación mínima entre bloques será cero. Dos bloques pueden ser
    adyacentes —uno termina en `x` y otro comienza en `x+1`—, pero nunca
    compartir una base. No se añadirá una distribución de espacios ni se
    rechazará la adyacencia. Las distancias interbloque realizadas podrán
    registrarse, pero no condicionarán el sorteo ni representarán LD entre
    bloques. Los perfiles de separación o dominios espaciales quedan en
    DTF-044/DTF-045.
16. La disposición se construirá después de generar la lista completa de
    longitudes. Los bloques se colocarán de mayor a menor longitud; para cada
    uno se sorteará uniformemente entre todos sus inicios todavía elegibles. Si
    una longitud no tiene ningún inicio disponible, se descartará la disposición
    completa y se repetirá con la misma lista de longitudes, sin sustituir ni
    reducir solo el bloque fallido. Los `block_id` se asignarán después de una
    colocación válida, en el orden canónico `sequence_id`, inicio y final, por lo
    que la prioridad de colocación no otorga identidad preferente a los bloques
    largos. El algoritmo no promete uniformidad entre todas las disposiciones
    globales posibles; sí debe ser versionado, reproducible y registrar sus
    intentos. La disposición se repetirá hasta el límite configurable
    `max_layout_attempts = 1000` de la política general de agotamiento.
17. `covered_bases` exigirá una cobertura entera exacta. Para
    `coverage_fraction`, primero se calculará
    `target_bases = floor(eligible_bases * coverage_fraction)`, de modo que la
    cobertura realizada nunca supere la fracción solicitada. Se sortearán
    longitudes mientras quepan en el resto; una longitud mayor que el resto se
    rechazará y se repetirá únicamente ese sorteo. La última longitud deberá ser
    exactamente igual al resto y pertenecer a algún intervalo admitido. Si el
    resto no pertenece a ninguno, se descartará la lista completa y se volverá
    a generar. Nunca se recortará un bloque ni se aceptará cobertura aproximada.
    Un objetivo menor que la longitud mínima o que no pueda expresarse como suma
    de longitudes admitidas deberá rechazarse antes de colocar bloques. La
    representabilidad numérica se comprobará por separado de la colocación
    física. Sea `S` el conjunto de longitudes enteras admitidas que caben al
    menos una vez en algún segmento elegible y sea `m = min(S)`. Se construirá el
    grafo exacto de residuos `0..m-1`: cada `l` de `S` aporta desde cada residuo
    `r` la transición `r -> (r+l) mod m` con coste `l`. Un camino mínimo desde
    cero obtendrá la menor suma alcanzable `d[r]` y conservará predecesores. El
    objetivo `T` es representable si y solo si `d[T mod m] <= T`; como `m`
    pertenece a `S`, la diferencia es un múltiplo de `m` y permite reconstruir
    una combinación testigo exacta. Todo el cálculo usará aritmética entera
    comprobada.

    La combinación testigo será evidencia del preflight: no sustituirá el
    muestreo ponderado de longitudes ni obligará al generador a usarla. Superar
    esta prueba tampoco demostrará que los bloques puedan colocarse
    conjuntamente en el espacio fragmentado. Si máscaras o fragmentación hacen
    agotar la disposición, se informará agotamiento de layout, nunca suma
    imposible. Si el cálculo exacto no puede ejecutarse por límites de recursos,
    se emitirá un diagnóstico de recursos y no se clasificará el objetivo como
    imposible. El máximo de intentos se rige por la política general de
    agotamiento siguiente.
18. La política de agotamiento tendrá tres límites enteros positivos,
    configurables y versionados: `max_block_attempts = 100` para regenerar la
    estructura de SNP de un bloque que no alcance dos sitios;
    `max_length_set_attempts = 1000` para regenerar una lista que no complete la
    cobertura exacta; y `max_layout_attempts = 1000` para recolocar la misma
    lista cuando el espacio se haya fragmentado. Antes de consumir intentos, el
    sistema rechazará los casos demostrablemente imposibles: cobertura superior
    al espacio elegible, objetivo menor que la longitud mínima, ausencia de una
    distancia SNP que quepa en alguna longitud admitida o cobertura demostrada
    no expresable mediante el grafo de residuos del punto 17. Para `block_count`,
    el máximo demostrable será
    `sum(floor(segment_length / min_allowed_length))` sobre los segmentos
    elegibles; superarlo será imposible, pero no superarlo no garantizará que una
    lista sorteada concreta pueda disponerse. Agotar un límite detendrá la
    generación sin reducir
    bloques, cobertura, SNP ni tolerancias. El diagnóstico y el manifiesto de
    fallo identificarán el contador, límite e intentos consumidos, semilla,
    objetivo y último motivo de rechazo. El usuario podrá aumentar
    explícitamente el límite y repetir. Solo una demostración exacta permitirá
    declarar imposibilidad; el agotamiento de colocación y la falta de recursos
    del preflight conservarán clasificaciones distintas.
19. El perfil de `r²` se aplicará de forma agregada por población e intervalo de
    distancia, no como un rango obligatorio para cada pareja. Para cada
    población se calculará el `r²` haplotípico faseado de todas las parejas de
    SNP pertenecientes al mismo bloque, se agruparán las parejas por distancia
    y se calcularán sobre el conjunto de todos los bloques los cuantiles
    realizados. Esos cuantiles se compararán con los solicitados. No se incluirán
    parejas de bloques diferentes ni se corregirá una pareja individual para
    acercarla a la mediana. La salida conservará, por población e intervalo, el
    número de parejas y los cuantiles solicitados y realizados. Los intervalos
    del perfil deberán cubrir todas las distancias intrabloque que el generador
    pueda producir. Los cinco cuantiles, la tolerancia, el mínimo evaluable y el
    tratamiento de SNP fijados quedan cerrados en los puntos 20–25 y el cálculo
    empírico exacto de cuantiles, en el punto 31.
20. Un SNP con frecuencia ALT 0 o 1 en una población se conservará, pero se
    considerará fijado y todas sus parejas se excluirán de la evaluación de
    `r²` de esa población porque el estadístico no está definido. Nunca se
    sustituirá por `r² = 0`. La misma pareja podrá ser evaluable en una población
    y excluida en otra. Por población e intervalo se registrarán al menos las
    parejas intrabloque totales, las evaluables y las excluidas por contener uno
    o dos SNP fijados. Si quedan menos parejas evaluables que el mínimo exigido,
    se aplicará la política de repetición o fallo que se cierre para ese mínimo.
21. Cada fila del perfil de `r²` declarará `min_evaluable_pairs`, un entero
    estrictamente positivo específico del intervalo. Los perfiles incorporados
    incluirán mínimos derivados y documentados junto con su receta; los perfiles
    suministrados deberán declararlos explícitamente. Si una población queda por
    debajo, sus cuantiles de ese intervalo no se considerarán válidos y se
    repetirá únicamente la construcción haplotípica de esa población conforme a
    su límite propio.
    Nunca se fusionarán intervalos, se rebajará el mínimo ni se incluirán parejas
    con SNP fijados para alcanzarlo. El agotamiento deberá identificar población,
    intervalo, mínimo, parejas evaluables y exclusiones por fijación.
22. La geometría de bloques y SNP será común e inmutable durante la construcción
    poblacional. Frecuencias, cuantización, inicialización, cálculo de cuantiles
    Type 7, puntuación y recocido se generarán de forma independiente por
    población, cada una con su propio `population_attempt_index`. El límite
    configurable y versionado
    `max_haplotype_attempts = 1000` se aplicará por población. Un fallo repetirá
    solo la población afectada y permitirá conservar internamente los resultados
    válidos de las demás. Agotar el límite de una población hará fallar el panel
    completo sin resortear las otras y sin relajar cuantiles, frecuencias ni
    mínimos. Informes y manifiesto identificarán población, intento y último
    incumplimiento.

    No se publicará ningún panel parcial. La publicación será atómica y solo se
    realizará después de que todas las poblaciones superen el mismo validador
    común. Solo entonces se construirá, por bloque, el catálogo final de HAP como
    unión densa de sus resultados poblacionales.
23. Cada fila del perfil contendrá obligatoriamente los cinco cuantiles
    `r2_p10`, `r2_p25`, `r2_p50`, `r2_p75` y `r2_p90`. Deben ser `f64` finitos
    en `[0,1]` y cumplir
    `r2_p10 <= r2_p25 <= r2_p50 <= r2_p75 <= r2_p90`; `p50` es la mediana. Esta
    versión no añade mínimos, máximos, `p5` ni `p95`. `pair_count` no es un
    objetivo de entrada: en un perfil empírico se conserva como evidencia de
    derivación y en un panel generado se informa como recuento realizado. Todo
    cuantil empírico se calcula mediante Hyndman–Fan Type 7 conforme a la fórmula
    reproducible del punto 31.
24. Cada fila declarará además `r2_max_abs_error`, un `f64` finito en `[0,1]`
    que se aplicará a sus cinco cuantiles. Para aceptar una población e
    intervalo deberá cumplirse
    `abs(r2_realizado - r2_solicitado) <= r2_max_abs_error` para cada cuantil.
    La tolerancia es una diferencia absoluta de `r²`, no un porcentaje relativo.
    Los perfiles incorporados incluirán valores derivados mediante una receta
    documentada y los suministrados deberán declararlos; no habrá tolerancia
    implícita. El incumplimiento de un solo cuantil rechazará el intento
    haplotípico de esa población.
25. Los intervalos de distancia usarán límites enteros positivos e inclusivos,
    estarán ordenados y formarán una partición contigua: la primera fila
    comenzará en 1 y cada `distance_min_bp` posterior será el
    `distance_max_bp` anterior más 1. No se admitirán huecos ni solapamientos.
    La última fila deberá alcanzar al menos `longitud_máxima_admitida - 1`,
    aunque podrá extenderse más para reutilizar el perfil. Todos los perfiles
    poblacionales de un mismo panel deberán compartir exactamente estos límites;
    solo podrán diferir en mínimos, tolerancias y cuantiles. Una fila situada
    completamente fuera de las distancias posibles de la configuración se
    marcará `no aplicable`, se informará y no tendrá que cumplir
    `min_evaluable_pairs`; cualquier fila que sí interseque el soporte posible
    permanecerá activa.
26. Las frecuencias alélicas marginales de los SNP intrabloque se generarán
    reutilizando los espectros de DTC-054. Cada SNP y población obtendrá su
    `alt_af` desde el espectro efectivo de esa población, independientemente de
    las demás; se admitirán 0 y 1. Estas frecuencias serán objetivos marginales
    que deberá conservar la distribución haplotípica mientras reproduce el
    perfil agregado de `r²`. Un intento incompatible regenerará las frecuencias
    y haplotipos únicamente de esa población sin cambiar bloques, posiciones SNP
    ni resultados de las demás. Los perfiles incorporados agruparán espectros
    alélicos y LD derivados de datos y filtros compatibles; el usuario podrá
    combinar los suyos, asumiendo que el intento puede resultar imposible. No se
    reutilizará
    `population_snvs.tsv`: esa matriz contiene loci independientes externos a
    los bloques; solo se reutiliza su modelo de espectros.
27. Cada población declarará un `panel_haplotype_copies` entero estrictamente
    positivo. Para cada bloque se construirán exactamente ese número de copias
    haplotípicas; las copias idénticas se agruparán y su frecuencia será
    `copy_count / panel_haplotype_copies`. El número de HAP distintos será un
    resultado, no otro parámetro. El catálogo canónico de un bloque será la
    unión de los vectores observados en todas las poblaciones, con frecuencia 0
    explícita donde un vector no aparezca, conforme al contrato denso de
    DTC-055. Todo perfil incorporado conservará el tamaño haploide usado en su
    derivación o construcción. La cuantización de `alt_af` a recuentos enteros
    queda fijada en el punto 28.
28. Cada `alt_af` solicitada se convertirá a un recuento entero mediante redondeo
    estocástico reproducible. Para `C = panel_haplotype_copies`, se calculará
    `x = alt_af * C` y el recuento será `floor(x)` más una Bernoulli cuya
    probabilidad sea la parte fraccionaria de `x`. Un `x` entero se conservará;
    `alt_af=0` producirá 0 y `alt_af=1`, C. El sorteo tendrá una subsemilla
    separada por población, bloque y SNP. La frecuencia efectiva `count/C` será
    la marginal exacta que conservará el constructor y se registrará junto con
    la frecuencia solicitada y el recuento. Una variante rara podrá quedar en 0
    y no se rescatará; aumentar C será la vía explícita para obtener mayor
    resolución.
29. El constructor representará cada población×bloque como una matriz binaria
    de C filas haplotípicas y una columna por SNP. Cada columna se inicializará
    con exactamente su recuento ALT efectivo, repartido reproduciblemente entre
    las filas. Un movimiento elegirá una columna y dos filas con valores
    opuestos e intercambiará 0 y 1; por construcción, ningún movimiento podrá
    alterar las marginales. Para cada población, intervalo activo y cuantil se
    calculará `delta = abs(realizado - objetivo)` y
    `exceso = max(0, delta - r2_max_abs_error)`. Cada población tendrá su propia
    puntuación, suma de `exceso²` para los cinco cuantiles y todos sus intervalos
    activos. Se recorrerán en orden canónico y se sumarán con Neumaier, sin
    reducción paralela; los intervalos `not_applicable` no participarán.

    La puntuación guiará la búsqueda, pero no sustituirá la comprobación por
    celda ni los mínimos evaluables. El resultado solo podrá aceptarse cuando la
    puntuación sea exactamente 0 y todos los mínimos sean válidos. Se aceptarán
    siempre las mejoras y los empates. Si la puntuación inicial ya es 0, el
    intento terminará sin proponer movimientos. En otro caso usará
    `T0 = max(puntuación_inicial, 1e-12)`. Un empeoramiento de magnitud
    `delta_score > 0` se aceptará si un uniforme `U` en `[0,1)` cumple
    `U < exp(-delta_score / T)`. Después de cada 1.000 movimientos propuestos,
    aceptados o no, la temperatura se multiplicará por 0,95. No se recalculará
    a partir de puntuaciones posteriores.

    El proceso comprobará y terminará inmediatamente cuando la puntuación sea
    0 y todos los mínimos y celdas sean válidos; solo entonces se
    agruparán filas idénticas y se publicarán sus frecuencias. Los mínimos de
    parejas evaluables se comprobarán antes de optimizar porque dependen de las
    marginales y posiciones, no de sus intercambios. El ajuste y la aceptación
    se evaluarán sobre todos los bloques de una sola población. La temperatura,
    la mejor puntuación y los mejores cuantiles pertenecerán también a esa
    población y su intento.
30. El objetivo usa exclusivamente `r²`, por lo que no impondrá el signo de la
    asociación. Fases en acoplamiento o repulsión con el mismo `r²` serán
    igualmente válidas y la semilla determinará el resultado. El control futuro
    mediante `r` firmado queda en DTF-046. El calendario y el límite de
    movimientos quedan fijados y versionados en este mismo punto, sin reabrir la
    puntuación, las marginales ni los criterios de aceptación ya acordados. El
    manifiesto registrará `T0`, factor, longitud de etapa, temperatura final y
    movimientos propuestos, aceptados y aceptados pese a empeorar.

    El presupuesto se calculará con aritmética entera comprobada. Primero se
    obtendrá para cada población
    `population_haplotype_cells = sum(panel_haplotype_copies * snp_count)` sobre
    todos sus bloques en orden canónico. Después se fijará para esa población
    `max_annealing_moves = min(10_000_000, max(100_000,
    100 * population_haplotype_cells))`. La fórmula será determinista,
    versionada y no configurable en esta versión; contará movimientos
    propuestos, no aceptados.

    Si se agota el presupuesto sin alcanzar el estado válido, se descartará el
    intento de esa población y se iniciará el siguiente con su
    `population_attempt_index` y subsemilla, hasta
    `max_haplotype_attempts = 1000` para esa población. Agotar sus intentos
    producirá `sin solución encontrada dentro del límite`, no `matemáticamente
    imposible`.
    El manifiesto registrará las celdas, la versión de la fórmula, el límite,
    los movimientos propuestos y la mejor puntuación y cuantiles alcanzados.
31. Todos los cuantiles empíricos usarán Hyndman–Fan Type 7. Para un vector no
    vacío ordenado `x[0..n-1]` y un cuantil `q`, si `n = 1` se devolverá `x[0]`;
    en otro caso se calculará `t = (n - 1) * q`, `lo = floor(t)`,
    `hi = ceil(t)`, `gamma = t - lo` y
    `Q = (1 - gamma) * x[lo] + gamma * x[hi]`. Se aplicará exactamente a
    `q = 0.10`, `0.25`, `0.50`, `0.75` y `0.90` en el generador, preflight,
    oráculo e informes.

    Antes de optimizar, el preflight evaluará cada pareja bialélica no fijada.
    Con `C` copias y recuentos ALT `a` y `b`, el recuento conjunto `h` recorrerá
    todos los enteros de
    `[max(0, a + b - C), min(a, b)]`. Para cada `h` se enumerará el `r²`
    permitido, usando `p_A = a/C`, `p_B = b/C`, `p_AB = h/C` y
    `r² = (p_AB - p_A*p_B)^2 / (p_A*(1-p_A)*p_B*(1-p_B))`. Para cada población
    e intervalo se formarán los
    vectores de mínimos y máximos de esos dominios por pareja, se ordenarán por
    separado y se les aplicará Type 7. Así se obtiene, para cada cuantil, un
    intervalo relajado conservador de valores alcanzables que ignora
    compatibilidades multilocus.

    Solo se declarará `matemáticamente imposible` si el intervalo solicitado
    `[target - r2_max_abs_error, target + r2_max_abs_error]` no interseca el
    intervalo alcanzable de algún cuantil, o si los marginales determinan menos
    parejas evaluables que `min_evaluable_pairs`. Superar el preflight no
    demuestra factibilidad: se ejecutará el recocido de esa población y agotar
    sus intentos será `sin solución encontrada dentro del límite`, no una prueba
    matemática. Ambos fallos detendrán la generación sin `--force`. El
    diagnóstico conservará población, intervalo, cuantil, rango alcanzable,
    rango solicitado con
    tolerancia, `C`, marginales relevantes y motivo.
32. El azar usará seis ámbitos físicos exactos con `derivation_version = 2`:
    `ld_geometry` para longitudes, disposición, saltos entre SNP y desplazamiento
    intrabloque; `ld_site_alleles` para elegir ALT; `ld_population_frequencies`
    para sortear `alt_af` y redondearlo; `ld_haplotype_initialization` para
    repartir los recuentos en las C filas; `ld_annealing_proposals` para elegir
    columnas, filas e intercambios; y `ld_annealing_acceptance` para consumir
    exclusivamente el uniforme de cada movimiento empeorante. Separar propuestas
    y aceptación impedirá que consumir o no ese uniforme desplace propuestas
    posteriores.

    Cada solicitud de subsemilla usará únicamente los identificadores estables
    que le sean aplicables entre `scenario_id`, `reference_id`, `population_id`,
    `block_id`, `site_id` y `population_attempt_index`. Quedan prohibidos ruta,
    reloj, hilo, *scheduler*, orden de iteración y posición accidental. Cada
    reintento incluirá su índice poblacional sin modificar los anteriores. Los ámbitos de
    frecuencias, inicialización, propuestas y aceptación incorporarán
    `population_id` y el `population_attempt_index` aplicable. Añadir, cambiar o
    reordenar una población no alterará la geometría común ni las frecuencias,
    haplotipos o recocido de las demás. Reordenar bloques ya identificados,
    secuencias o salidas tampoco cambiará identidades no relacionadas; cambiar
    cantidad, cobertura o distribución de bloques sí podrá reconstruir toda la
    geometría.

    Los catálogos v1 y v2 permanecerán intactos, legibles y reproducibles; no se
    reinterpretarán sus semillas. La futura materialización exacta será
    `catalogs/subseed-scopes-v3.tsv`, con `catalog_version = 3` y
    `derivation_version = 2`. Contendrá los ámbitos anteriores necesarios y los
    seis ámbitos LD nuevos, con vectores normativos propios. Una ejecución con
    generación automática seleccionará explícita y coherentemente catálogo 3 y
    derivación 2; no mezclará silenciosamente versiones. Configuración efectiva
    y manifiesto registrarán catálogo, derivación, contexto, ámbito, componentes
    y subsemillas conforme a DTC-037/DTC-038. El fichero real no se creará hasta
    superar las puertas normativas y autorizar su implementación.

33. Todo panel generado deberá superar una validación independiente que vuelva
    a leer desde cero los cinco ficheros publicados y no use matrices, cachés ni
    puntuaciones del optimizador. Aplicará el validador DTC-055, reconstruirá
    recuentos y frecuencias ALT desde HAP y frecuencias, recalculará parejas
    intrabloque, exclusiones por fijación, `r²`, cuantiles y mínimos, y verificará
    las tolerancias solicitadas. En pruebas y aceptación, fixtures pequeños con
    `r²` 0, 1 e intermedio, SNP fijados y varias poblaciones se expandirán a
    haplotipos faseados y se contrastarán también con PLINK 2. PLINK será solo
    un validador externo de CI/aceptación, no una dependencia runtime. El panel
    no se considerará válido si la recomputación desde sus artefactos no
    reproduce marginales y perfil de LD. Fixtures pequeños enumerarán además de
    forma exhaustiva todas las realizaciones compatibles y comprobarán tanto la
    fórmula Type 7 como que el preflight nunca clasifique `matemáticamente
    imposible` una realización factible.
    Un oráculo independiente enumerará o resolverá por programación dinámica
    fixtures pequeños de cobertura. Incluirá intervalos degenerados y con
    huecos, recompondrá el testigo hasta la suma exacta y protegerá la frontera
    entre suma imposible, layout agotado y preflight no ejecutable por recursos.
    Los oráculos de aislamiento comprobarán que cambiar el perfil, reordenar o
    añadir una población conserva los bytes de los artefactos poblacionales no
    relacionados cuando el contrato físico los separe y, en otro caso, conserva
    exactamente sus filas o proyecciones canónicas. También forzarán el fallo de
    una población y verificarán que no se publica ningún panel parcial ni se
    resortean las poblaciones ya válidas.

## 2. Representación propuesta del perfil

Un perfil relacionará intervalos de distancia física con la distribución
observada de `r²`. Guardar solo el mínimo y el máximo no basta: en un conjunto
empírico ambos extremos pueden aproximarse a 0 y 1 en muchos intervalos. Se
conservan cinco cuantiles obligatorios. Los recuentos de parejas se registran
como evidencia de derivación o resultado, no como objetivo.

`ld_profile.tsv` usará esta cabecera exacta:

```tsv
population_profile_id	distance_min_bp	distance_max_bp	min_evaluable_pairs	r2_max_abs_error	r2_p10	r2_p25	r2_p50	r2_p75	r2_p90
```

Ejemplo exclusivamente ilustrativo; no contiene valores calculados de
1000 Genomas:

```tsv
population_profile_id	distance_min_bp	distance_max_bp	min_evaluable_pairs	r2_max_abs_error	r2_p10	r2_p25	r2_p50	r2_p75	r2_p90
CEU	1	1000	1000	0.05	0.08	0.18	0.42	0.72	0.91
CEU	1001	5000	500	0.05	0.03	0.09	0.25	0.51	0.76
YRI	1	1000	1000	0.05	0.05	0.13	0.34	0.63	0.86
YRI	1001	5000	500	0.05	0.02	0.07	0.20	0.44	0.70
```

La evaluación agregará todas las parejas intrabloque de la población en el
intervalo correspondiente. Los valores individuales pueden quedar fuera de los
cuantiles centrales. Los intervalos inclusivos y contiguos y la tolerancia
absoluta por fila ya están acordados. Todas las poblaciones del recurso
declararán exactamente los mismos intervalos; cada una podrá declarar su propio
mínimo, tolerancia y cuantiles.

`population_profile_id` cumplirá DTC-039. Los límites y
`min_evaluable_pairs` serán enteros decimales positivos canónicos. La tolerancia
y los cinco cuantiles serán `f64` finitos en `[0,1]`, con serialización de ida y
vuelta exacta; los cuantiles serán no decrecientes. La clave será
`(population_profile_id, distance_min_bp, distance_max_bp)` y las filas se
ordenarán por esos campos. Los intervalos comenzarán en 1, serán contiguos y el
último cubrirá como mínimo `max_block_length - 1`.

El fichero seguirá el mismo contrato TSV canónico de los otros dos: UTF-8 sin
BOM, tabuladores y LF final, sin comillas, vacíos, comentarios ni columnas
adicionales.

### 2.1. Perfil humano aproximado basado en Ragsdale y Gravel

La primera versión podrá incluir el perfil identificado de forma definitiva como
`human_ld_ragsdale_approx_v1`. Se inspirará en las curvas de caída de LD de
Ragsdale y Gravel (2020) para cinco poblaciones de 1000 Genomas. Se usará la
curva de su estimador `r²_÷` como tendencia media por distancia; no se usará
la curva Rogers–Huff ni `sigma_D²` como si fueran la misma magnitud.

El artículo publica curvas agregadas, no los cuantiles por pareja exigidos por
este borrador. Para cada población e intervalo, la curva digitalizada o ajustada
aportará una media `mu` en `[0,1]`. Si `0 < mu < 1`, se construirá una
distribución beta con `alpha = mu * kappa` y
`beta = (1 - mu) * kappa`, donde `kappa` es finito y estrictamente positivo.
Los cuantiles `p10`, `p25`, `p50`, `p75` y `p90` serán los inversos de su
función de distribución acumulada en 0,10; 0,25; 0,50; 0,75 y 0,90. Si `mu` es
exactamente 0 o 1, la distribución será degenerada y los cinco cuantiles serán
iguales a `mu`.

El perfil incorporado usará el único valor `kappa = 8.0` para todas las
poblaciones e intervalos. Este valor forma parte de la versión del perfil: no
será configurable, inferido ni ajustado durante una simulación. Un perfil
suministrado podrá declarar directamente sus cinco cuantiles y no estará
obligado a generarlos mediante esta distribución. El perfil conservará la cita,
la tabla digitalizada o reconstruida, el procedimiento de ajuste, los
intervalos de distancia, los parámetros y los SHA-256.

Este perfil pretende producir LD humano plausible, no reproducir exactamente
1000 Genomas. Un perfil empírico recalculado desde VCF faseados podrá añadirse
después como otro recurso y nunca reemplazará silenciosamente al aproximado.

## 3. Metadatos propuestos

Cada perfil incorporado o suministrado será un paquete versionado con
exactamente estos cuatro ficheros obligatorios:

1. `block_lengths.tsv`, con los intervalos y pesos globales de longitud;
2. `snp_gaps.tsv`, con los intervalos y pesos globales de distancia entre SNP
   consecutivos;
3. `ld_profile.tsv`, con las identidades poblacionales del recurso, los
   intervalos de distancia y los cuantiles de `r²`;
4. `profile.metadata.json`, con clase, versión, `profile_id`, naturaleza
   incorporada o suministrada, carácter aproximado o empírico, procedencia,
   receta, parámetros y, para cada TSV, clase, ruta relativa, tamaño y SHA-256.

El JSON no duplicará las tablas y el paquete no admitirá archivos adicionales.
Los perfiles incorporados y suministrados atravesarán el mismo validador. La
máscara de exclusión no formará parte del paquete porque depende de la
referencia concreta de la ejecución; se declarará, identificará y verificará
por SHA-256 como una entrada separada.

Cada perfil deberá permitir reconstruir su procedencia. El esquema base de
`profile.metadata.json` queda representado por este ejemplo:

```json
{
  "profile_kind": "ld_generation_profile",
  "schema_version": "1",
  "profile_id": "human_ld_ragsdale_approx_v1",
  "package_origin": "builtin",
  "derivation_kind": "approximate",
  "statistic": "phased_haplotype_r2",
  "variant_class": "biallelic_snp",
  "population_profile_ids": ["CEU", "CHB", "MXL", "PUR", "YRI"],
  "quantile_model": {
    "kind": "beta_from_mean_curve",
    "kappa": 8.0
  },
  "provenance": {
    "citation": "Ragsdale and Gravel (2020)",
    "source_detail": "Figure 2C",
    "derivation_recipe": "..."
  },
  "artifacts": [
    { "class": "block_lengths", "path": "block_lengths.tsv", "size": 0, "sha256": "..." },
    { "class": "snp_gaps", "path": "snp_gaps.tsv", "size": 0, "sha256": "..." },
    { "class": "ld_profile", "path": "ld_profile.tsv", "size": 0, "sha256": "..." }
  ]
}
```

`profile_kind`, `schema_version`, `statistic` y `variant_class` tendrán
exactamente los valores mostrados. `profile_id` y cada identidad poblacional
cumplirán DTC-039. `package_origin` admitirá `builtin` o `supplied` y
`derivation_kind`, `approximate` o `empirical`. La lista poblacional será no
vacía, única, ordenada por ASCII e idéntica al conjunto de
`population_profile_id` de `ld_profile.tsv`.

`quantile_model.kind` distinguirá `beta_from_mean_curve` de
`direct_quantiles`. El primero exige un `kappa` finito y positivo; el perfil de
Ragsdale y Gravel exige exactamente 8.0. Un perfil aproximado exigirá en
`provenance` `citation`, `source_detail` y `derivation_recipe` no vacíos. Uno
empírico exigirá, además de receta y procedencia, conjunto y versión de origen,
referencia, muestras, regiones y filtros explícitos.

`artifacts` contendrá exactamente las tres clases mostradas, una vez cada una,
con sus nombres canónicos como rutas relativas simples. `size` será un entero
no negativo y `sha256`, 64 caracteres hexadecimales minúsculos verificados
contra los bytes. El objeto será cerrado, sin campos desconocidos ni tablas
duplicadas. Se serializará en UTF-8 sin BOM y con LF final.

El perfil aproximado de Ragsdale y Gravel es un producto sintético versionado
por el proyecto, no un fichero oficial de 1000 Genomas ni del artículo.

### 3.1. Máscara externa de exclusión

La máscara opcional se suministrará por separado como
`excluded_regions.tsv`, con esta cabecera exacta:

```tsv
region_id	sequence_id	start	end
MASK_001	chr1	100000	120000
MASK_002	chr2	450000	451500
```

El fichero seguirá el contrato TSV canónico del perfil. `region_id` cumplirá
DTC-039 y será único; `sequence_id` existirá en el FASTA. `start` y `end`
serán enteros decimales positivos en coordenadas base 1 cerradas, con
`start <= end`, y el intervalo completo pertenecerá a la secuencia. Las filas
se ordenarán por `sequence_id` ASCII, `start`, `end` y `region_id`. Los
intervalos no podrán solaparse, aunque sí ser adyacentes.

El programa calculará el SHA-256 de la máscara y lo registrará junto con
`reference_id` y `reference_sha256` en la configuración efectiva y el
manifiesto. La ausencia del fichero significa que no hay máscara aportada por
el usuario; no desactiva la exclusión automática de bases `N`, loci reservados o
bloques ya aceptados. Una máscara nunca puede volver elegible uno de esos
espacios.

## 4. Configuración prevista

La configuración separará la geometría común de la asignación del perfil a
cada población:

```toml
[initial_population.ld_blocks]
source = "generated"
profile = "builtin:human_ld_ragsdale_approx_v1"
block_count = 100
excluded_regions_path = "inputs/excluded_regions.tsv"

[initial_population.ld_blocks.populations.POP_A]
population_profile_id = "CEU"

[initial_population.ld_blocks.populations.POP_B]
population_profile_id = "YRI"
```

Un paquete suministrado se seleccionará en la misma sección:

```toml
[initial_population.ld_blocks]
source = "generated"
profile_path = "profiles/custom-human-ld"
covered_bases = 5_000_000
```

`profile` y `profile_path` serán obligatorios y mutuamente excluyentes. También
se declarará exactamente uno de `block_count`, `covered_bases` o
`coverage_fraction`. `excluded_regions_path` será opcional. Las rutas serán
relativas al documento propietario conforme al contrato común que cierre
DTPD-101.

Cada población simulada alcanzada declarará un `population_profile_id` existente
en `ld_profile.tsv`. La geometría común se resolverá una sola vez y esas
asignaciones no cambiarán bloques ni posiciones. Las sustituciones seguirán
DTC-033: podrán cambiar la asignación de una población, pero nunca modificar
parcialmente las tablas globales ni los parámetros de un perfil incorporado.
Para cambiar cualquiera de ellos se suministrará un paquete propio completo.

Un `source` distinto de `generated`, una identidad poblacional ausente o una
combinación de campos incompatible se rechazará antes de generar. Marimo
expondrá las mismas alternativas y validaciones.

`block_lengths.tsv` usará esta cabecera exacta:

```tsv
length_min_bp	length_max_bp	weight
1000	4999	0.20
5000	9999	0.50
10000	25000	0.30
```

Primero se seleccionará un intervalo según sus pesos y después se sorteará una
longitud entera dentro del intervalo. Los valores son ilustrativos. La tabla
será global para el panel; cada población conservará por separado su perfil de
`r²`. La tabla se valida completa antes de generar bloques: límites positivos e
inclusivos, intervalos ordenados y sin solapamientos, huecos permitidos, pesos
finitos y positivos que suman 1 conforme a DTC-049. `length_min_bp` y
`length_max_bp` son enteros decimales positivos canónicos y `weight` es un
`f64` decimal finito con serialización de ida y vuelta exacta. Las filas se
ordenan por `length_min_bp` y después por `length_max_bp`; no se admiten
intervalos duplicados.

El fichero será UTF-8 sin BOM, separado por tabuladores y terminado por LF.
No admitirá comillas, campos vacíos, comentarios ni columnas adicionales.

`snp_gaps.tsv` usará esta cabecera exacta:

```tsv
gap_min_bp	gap_max_bp	weight
1	10	0.05
11	100	0.25
101	500	0.40
501	1000	0.20
1001	5000	0.10
```

Los valores son ilustrativos. La distribución será global para el panel,
suministrada o incorporada y nunca predeterminada. El fichero usará exactamente
las mismas reglas físicas y semánticas que `block_lengths.tsv`:
`gap_min_bp`/`gap_max_bp` serán enteros decimales positivos canónicos,
`weight` será un `f64` finito positivo con ida y vuelta exacta, y las filas se
ordenarán por ambos límites. Los intervalos serán inclusivos, únicos y no
solapados; se admitirán huecos. La suma seguirá DTC-049 y una tabla inválida se
rechazará sin normalización ni corrección silenciosa.

Se sortean saltos hasta que el siguiente no cabe, se descarta ese salto sin
truncarlo y se desplaza uniformemente el tramo aceptado dentro del bloque. Se
exige un mínimo de dos sitios; la regeneración se limita mediante
`max_block_attempts = 100`.

## 5. Resultados y trazabilidad previstos

El generador producirá el paquete canónico de cinco ficheros fijado por
DTC-055 y cuatro informes TSV resumidos. Los informes no formarán parte del
directorio canónico del panel: el manifiesto de la ejecución generadora los
enlazará mediante clase, ruta, tamaño y SHA-256.

`block_generation_summary.tsv` contendrá una fila por bloque y esta cabecera
exacta:

```tsv
block_id	sequence_id	start	end	length	snp_count	first_snp	last_snp	layout_attempt
```

`allele_frequency_realization.tsv` contendrá una fila por población, bloque y
SNP y esta cabecera exacta:

```tsv
population_id	block_id	site_id	requested_alt_af	alt_count	panel_haplotype_copies	effective_alt_af	fixed
```

`ld_profile_realization.tsv` contendrá una fila por población e intervalo de
distancia y esta cabecera exacta:

```tsv
population_id	distance_min_bp	distance_max_bp	status	total_pairs	evaluable_pairs	excluded_one_fixed_snp_pairs	excluded_two_fixed_snps_pairs	min_evaluable_pairs	r2_max_abs_error	target_p10	target_p25	target_p50	target_p75	target_p90	realized_p10	realized_p25	realized_p50	realized_p75	realized_p90	passed
```

`status` distinguirá `evaluated` de `not_applicable`. La fila conservará por
separado las parejas excluidas por uno y por dos SNP fijados.

`generation_attempts.tsv` contendrá una fila por etapa e identidad afectada y
esta cabecera exacta:

```tsv
stage	population_id	block_id	attempt_limit	attempts_used	final_status	last_rejection
```

Las filas de los cuatro informes seguirán un orden canónico. Su contenido y
los cinco artefactos del panel deberán bastar para reconstruir por qué la
realización se aceptó o rechazó. No se publicará por defecto una tabla con el
`r²` de cada pareja, porque su tamaño puede crecer de forma cuadrática; ese
diagnóstico detallado queda aplazado en DTF-047.

## 6. Decisiones pendientes

No quedan elecciones funcionales o técnicas materiales abiertas en este
borrador. Esto no lo convierte en contrato: siguen pendientes su traslado y
revisión en especificaciones, trazabilidad, arquitectura y decisión técnica, y
solo después podrá redactarse un plan ejecutable.
