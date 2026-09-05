# Plan de implementación F5 0.1 — población fundadora y LD sintético

## 1. Estado, autoridad y alcance

- **Fecha:** 1 de septiembre de 2026.
- **Estado:** ejecutable de forma secuencial después del cierre verificable de
  F4. No autoriza *commit*, publicación externa ni inicio de F6/F7.
- **Fuentes normativas:** REQ-EVO-32–34/38–41/43/44/46–53/74–78 y los REQ de
  configuración, individuo, validación, manifiesto e interfaz asignados a F5 en
  la matriz; ACE-89/90/104/117, DTC-089, arquitectura §15 y hoja de ruta §F5.

F5 genera internamente bloques y haplotipos SNP A/C/G/T, ajusta marginales y
media `r²(d)` y asigna fundadores completos e intactos. No admite panel,
haplotipos externos, `alt_af` intrabloque, *indels*, repeticiones, SV, BAM/CRAM,
mosaicos, recombinación ni genealogía. F6 añade variación compleja; F7.1 diseña
la primera meiosis. F5.3 materializa identidad y pertenencia poblacional antes
de artefactos; F5.5 solo audita que el paquete resultante no las modifique.

La jerarquía única es F5.0–F5.6. Cada hito depende únicamente del anterior; no
existe dependencia de revisión futura ni ciclo F5→F6/F7→F5. Los pasos de
esquema crean versiones nuevas sin reescribir las ya publicadas. No queda
ninguna decisión funcional o técnica F5 pendiente: DTPD-115 pertenece a F9 y no
es una puerta de admisión ni de ejecución de este plan.

## 2. Criterios de orden y puertas comunes

Primero se cierran los contratos que consumen todos los componentes. Después se
construyen layout y marginales; luego H×S/LD; después se materializan identidad
y pertenencia fundadora, y solo entonces se construye y publica el paquete
atómico. Ningún paso posterior puede modificar identidad, pertenencia, copias o
asignaciones ya incluidas en ese paquete. Los riesgos de error silencioso
—decimal de cobertura, referencia, fase, D/`r²`, reducción numérica, RNG e IDs—
se contrastan con lectores u oráculos que no comparten el productor.

Cada paso ejecuta las pruebas específicas indicadas y, antes de su pausa:
`cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`
y el subconjunto `cargo test` pertinente. El cierre ejecuta la suite completa.

## 3. F5.0 — contratos y preflight

- ☐ **F5.0.1.** Crear y cerrar los contratos externos completos general v3,
  escenario F5 v1, efectivo v5, manifiesto general v5, perfil F5 v1, los siete
  TSV, fallo v1 y `truth.jsonl` v2; incorporar `f5_scenario` v1 como el único
  documento enlazado y resolver cada ruta desde el documento que la declara,
  sin abrir workspace ni modificar general v2/efectivo v4/manifiesto v4 ni
  `truth.jsonl` v1.
  - Depende de: cierre F4.
  - Procede de: REQ-EVO-74/77, REQ-UI-15 y DTC-089 «Escenario TOML».
  - Cambia: enlace general→escenario, adaptador TOML, `DeclaredConfig`,
    `EffectiveConfig` y todos los esquemas F5 definitivos; el general posee semilla/ploidía/referencia,
    el escenario perfil/máscara/política, y efectivo/manifiesto comprueban sus
    copias iguales con versiones general 3, escenario 1, efectivo 5/F5 1 y
    manifiesto 5/F5 1, verdad 2 y
    `truth_projection="f5-founder-truth-projection-v2"`. Ningún hito posterior amplía
    estos esquemas.
  - Se comprueba con: pareja mínima válida; rechazo de cero, lista o varios
    escenarios; rechazo de v2/v4 como F5; inmutabilidad byte a byte de los tres
    esquemas heredados; campos ausentes/desconocidos;
    one-of de objetivo, máscara/hash, suma de fundadores, IDs, F/P y límites;
    CLI relativa, escenario relativo al general, recurso relativo al escenario,
    duplicación prohibida; `seed` exacta `^[0-9a-f]{64}$` e igualdad entre
    general resuelto, raíz efectiva, `f5.seed` efectivo y `f5.seed` manifiesto;
    discordancia de cualquier otra copia efectiva/manifiesto; y verdad v2 con
    población/copia/locus estructurados. Los esquemas legacy se comparan byte a byte.
- ☐ **F5.0.2.** Implementar el lector único de `profile.metadata.json` y los
  cinco TSV, con sus cabeceras, dominios, orden, claves, normalización y hashes.
  - Depende de: F5.0.1.
  - Procede de: REQ-EVO-74/75, ACE-117 y DTC-089 «Paquete de perfiles».
  - Cambia: recursos F5 y validación de entrada; no genera datos biológicos.
  - Se comprueba con: corpus válido y una mutación aislada de BOM, LF, cabecera,
    fila, intervalo, peso, población, versión, tamaño y SHA-256.
- ☐ **F5.0.3.** Publicar los catálogos inmutables `f5-rng-catalog-v1` e ID F5,
  incluidos material tipado, prefijos y vectores por ámbito.
  - Depende de: F5.0.2.
  - Procede de: REQ-EVO-77, REQ-VAL-56 y DTC-037/039/089.
  - Cambia: catálogos RNG/ID y fixtures, incluidas `FCP_` por secuencia/fase y
    `FHA_` como asignación de bloque distinta; no cambia los catálogos previos.
  - Se comprueba con: derivador independiente BLAKE3/ChaCha12, enteros con
    rechazo, Fisher–Yates, Bernoulli 0/1 sin consumo, colisión dirigida y bytes
    canónicos exactos de layout (`u64` big-endian) y matriz (cabecera H/S,
    orden de filas, bits MSB y relleno cero).
- ☐ **F5.0.4.** Implementar el preflight puro
  `f5-resource-accounting-v7` y su ledger acumulativo, que resuelve racional de
  cobertura, capacidad, dimensiones y cotas de operaciones, pico de memoria y
  bytes publicados con símbolos dimensionalmente distintos. Este núcleo no
  abre ficheros: recibe observaciones de las fases físicamente acotadas de
  F5.1.2 y no acepta tamaños obtenidos por `metadata` como autoridad.
  - Depende de: F5.0.3.
  - Procede de: REQ-EVO-74/76/77, ACE-117 y arquitectura §15.
  - Cambia: admisión, escáner secuencial referencia+máscara, lista acotada de
    tramos elegibles y catálogo diagnóstico F5.
  - Se comprueba con: decimal extremo, T=0, soporte discontinuo, T>E,
    layout imposible, `S_bound` previo al RNG, H×S excesivo, `bin_count` por
    población, `reference_bases_processed`, bytes físicos FASTA/máscara y de
    TOML/metadata/TSV de perfil observados por el lector,
    `tract_count`, pico de IDs FASTA
    retenidos, filas de
    perfil/normalización/TSV/JSON, aportaciones ajenas F5, cada término y
    overflow de las fórmulas,
    fronteras inclusivas 2^28 unidades/2^30 bytes de trabajo/2^31 bytes de
    salida —incluidas copias/configuración, salidas activadas no F5, manifiesto
    y `truth.jsonl`—, límites menores declarados,
    cero consulta de RAM/CPU o perfil externo, cero consumo RNG y cero destino.
    Dos FASTA con la misma geometría y bases pero descripciones mínima y de
    80 MiB aumentan la cota exactamente por sus bytes recorridos en preflight y
    recuperación REF; el límite igual se admite y el límite menos uno se
    rechaza antes de RNG. La misma prueba se aplica a una máscara con igual
    geometría y distinto tamaño de `region_id`.
    Un fixture usa el mismo escenario con `block_count=1`, el mismo FASTA de
    una secuencia y dos máscaras con igual E: una deja un único tramo compacto
    y la otra 4096 tramos de longitud `Lmin`, separados por una base
    enmascarada. Debe reconstruir
    `reference_bases_processed=4096*Lmin+4095` y, manteniendo iguales todos los
    demás valores, obtener una diferencia exacta de operaciones
    `4095*(3+max_layout_attempts+5*max_layout_attempts*block_count)` —emisión,
    capacidad/request, copias y las dos consultas directas—,
    `512*4095` bytes de trabajo y cero bytes de salida entre ambas cuentas. Con
    memoria y salida holgadas, la máscara fragmentada se acepta con el límite de
    operaciones igual a su cota y se rechaza con cota−1; en casos separados,
    con operaciones y salida holgadas, se acepta con el límite de memoria igual
    a su cota y se rechaza con cota−1. Cada rechazo usa `GPS-RES-101` antes de
    RNG y destino. El fixture reconstruye la fórmula con unidades homogéneas y
    comprueba igualdad escalar exacta del `resource_preflight` entre efectivo y
    manifiesto; el sobre fijo de salida no depende de la fragmentación. El
    fixture fuerza que `final_formula_peak` sea el máximo para que la diferencia
    de memoria sea atribuible a la fragmentación, no a otro pico de fase.

> **Revisión F5.0.** Aprobar TOML, perfiles, vectores, diagnósticos y pureza del
> preflight. No continuar sin reconstrucción independiente de cada contrato.

## 4. F5.1 — geometría y marginales

- ☐ **F5.1.1.** Construir longitudes exactas condicionadas para
  `covered_bases`/`coverage_fraction` y estocásticas para `block_count`.
  - Depende de: F5.0.4.
  - Procede de: REQ-EVO-74, ACE-117 y DTC-089 «Layout».
  - Cambia: geometría F5.
  - Se comprueba con: suma T exacta, residual final, pesos restringidos,
    overflow y un oráculo entero que prueba la factibilidad previa.
- ☐ **F5.1.2.** Integrar la admisión F5 físicamente acotada y, solo tras
  sellarla, colocar bloques completos e insertar SNP mediante gaps con intentos
  transaccionales y orden/IDs canónicos. Este hito sustituye las correcciones
  incrementales v5/v6; se ejecuta en los subhitos siguientes y no incorpora
  crates nuevas.
  - Depende de: F5.1.1.
  - Procede de: REQ-EVO-74/76/77, REQ-VAL-54/56, REQ-MAN-70/92,
    REQ-SAL-02, DTC-040/049/089 y arquitectura §15.
  - Cambia: frontera CLI F5, carga de perfil, preflight, binding físico,
    geometría e identidades `BLK_`/`SITE_`; no cambia RNG, biología, formatos de
    perfil ni algoritmos de layout.

  - ✅ **F5.1.2.1.** Sustituir la sonda CLI por el dispatcher completo acotado y
    construir el bootstrap general+escenario.
    - Depende de: F5.1.1 y núcleo F5.0.4.
    - Procede de: arquitectura §15 y DTC-089 «Admisión físicamente acotada».
    - Cambia: una sola apertura del general, límite inclusivo 65536+byte de
      exceso, SHA-256/tamaño sobre el mismo handle, clasificación tras EOF y la
      misma lectura acotada del escenario. Tras parsearlo crea el ledger y
      contrasta el coste del bootstrap antes de cualquier otra apertura.
    - Corrección de revisión v7: la misma frontera física precede también a los
      cargadores públicos que todavía devuelven el general v2; hasta EOF no se
      atribuye el documento a F5. El dispatcher conserva la ruta de autoría v2
      solicitada para `source_path` y `configuration_input.original_path`, y
      usa su forma normalizada solo para la apertura física y los enlaces F5.
    - Se comprueba con: F5 cuyos discriminadores aparecen después de 4096
      bytes; general y escenario de 65536/65537 bytes; TOML inválido al final;
      espías de apertura que prueban cero perfil/referencia/máscara/destino/RNG
      al fallar; v3 exacto, 65537 y muy grande por cada cargador público
      pertinente, sin parser ni bytes posteriores al byte de exceso; diagnóstico
      neutral para v2/v3/inclasificable grande; ruta relativa v2 publicada; y
      fixture de anchura fija que alcanza el punto fijo de
      `max_working_memory_bytes`, admite límite exacto y rechaza límite−1.

  - ✅ **F5.1.2.2.** Implementar `BoundedInput`,
    `PhysicalInputIdentity` y las transiciones del ledger sin `fs::read` ni
    longitud de `metadata` como autorización.
    - Depende de: F5.1.2.1.
    - Procede de: DTC-040, REQ-VAL-54/56 y DTC-089.
    - Cambia: contador/hash sobre los mismos bytes, comprobación prospectiva de
      cada crecimiento, lectura `límite+1`, aritmética comprobada y diagnóstico
      con algoritmo, fase, rol, dimensión, acumulado y límite.
    - Se comprueba con: lector que entrega fragmentos de 1 byte, EOF exacto,
      byte de exceso, error intermedio, overflow dirigido y mutación entre una
      consulta de atributos y la lectura. El oráculo independiente recalcula
      tamaño y SHA desde el fichero y demuestra que los atributos no deciden la
      admisión.

  - ✅ **F5.1.2.3.** Leer el paquete de perfil bajo el presupuesto restante.
    - Depende de: F5.1.2.2.
    - Procede de: REQ-EVO-74/75, DTC-049 y DTC-089 «Paquete de perfiles».
    - Cambia: metadata bajo límite derivado del ledger; reserva previa del coste
      total conocido de los cinco enlaces; TSV secuenciales bajo su tamaño
      declarado+1; liberación del búfer antes de la tabla siguiente; retención
      exclusiva de tipos de dominio y auditorías.
    - Se comprueba con: metadata o cada TSV mayor que el restante, tamaño/SHA
      declarado incorrecto, padding semánticamente neutro, presupuesto exacto y
      exacto−1, pico por tabla y orden de aperturas. Una prueba independiente
      registra handles y demuestra que ninguna tabla posterior se abre después
      del primer rechazo.
    - Estado reabierto (5 de septiembre de 2026): la revisión integral invalidó
      la evidencia de cierre del día 4 mediante A1–A7. Las correcciones y sus
      regresiones quedan registradas en
      `evidence/f5.1.2/aplicacion-correcciones-integrales-2026-09-05.md`.
      Una revisión independiente posterior detectó R1–R3 en versión externa,
      reserva prospectiva de TSV y dimensión diagnóstica; las correcciones
      locales también constan allí. Una revisión de cierre posterior rechazó
      todavía R1 porque la proyección escalar v7 había heredado el campo
      interno `input_parser_peak_bytes`: se retiró del productor y del esquema
      efectivo —y, por referencia, del manifiesto— y se añadió un oráculo de
      claves literal e independiente para ambas rutas. Esta proyección de .3
      no se presenta como el contrato final: `input_identities` y `phases`
      siguen reservados a .5. La revisión independiente de cierre del 5 de
      septiembre de 2026 aprobó la pausa A: verificó los 27 campos escalares
      publicables en .3, la reserva prospectiva de TSV, la dimensión
      diagnóstica y las rutas públicas. H10 sigue reservado a F5.1.2.6 y no se
      presenta como gate verde.

  **✅ Pausa de revisión F5.1.2-A.** Aprobar dispatcher, bootstrap, lector, perfil
  y ledger. Exigir ausencia mecánica de `fs::read` y de sondas parciales en la
  ruta F5, pruebas exacto/exacto−1 y cero apertura posterior. No continuar con
  referencia hasta corregir todos los hallazgos.

  - ☐ **F5.1.2.4.** Integrar referencia y máscara en el mismo ledger.
    - Depende de: revisión F5.1.2-A aprobada.
    - Procede de: REQ-EVO-74/76, DTC-040 y DTC-089 «Layout».
    - Cambia: primer FASTA y única máscara en streaming; reserva de dos visitas
      por byte FASTA antes de leerlo; comprobación previa a retener IDs, filas y
      tramos; receipt de tamaño/SHA observado en el mismo handle.
    - Se comprueba con: descripciones FASTA y `region_id` largos, línea de
      máscara excesiva, referencia/máscara vacía o truncada, límite exacto y
      exacto−1 en operaciones y memoria, SHA incorrecto y espías que prueban
      cero RNG/destino. El oráculo recuenta bytes, bases y tramos sin reutilizar
      estructuras productivas.

  - ☐ **F5.1.2.5.** Sellar `F5AdmissionReceipt`, exigir
    `SealedF5Admission` en layout y acotar la recuperación REF.
    - Depende de: F5.1.2.4.
    - Procede de: REQ-VAL-54/56, REQ-MAN-70/92, REQ-SAL-02 y arquitectura §15.
    - Cambia: binding conjunto de identidad semántica y de rol+ruta
      normalizada+tamaño+SHA de general, escenario, metadata, cinco TSV,
      referencia y máscara opcional; efectivo/manifiesto v5 copian algoritmo
      v7, identidades y fases. La segunda pasada FASTA usa tamaño recibido+1 y
      exige tamaño/SHA idénticos antes de aceptar REF.
    - Se comprueba con: cruce de dos TOML semánticamente iguales con distinto
      comentario/padding; intercambio de TSV; modificación, truncado y
      crecimiento del FASTA entre pasadas; ausencia/presencia de máscara; y
      desigualdad dirigida entre efectivo y manifiesto. Ningún receipt v5/v6 se
      acepta ni se migra.

  **Pausa de revisión F5.1.2-B.** Aprobar el modelo local de cambios
  accidentales, el orden de fases, las identidades físicas, la reserva de la
  segunda pasada y el binding sellado. No continuar con layout si un consumidor
  todavía acepta escalares o receipts sueltos.

  - ☐ **F5.1.2.6.** Revalidar el layout productivo y cerrar la sustitución v7.
    - Depende de: revisión F5.1.2-B aprobada.
    - Procede de: REQ-EVO-74/77 y DTC-089 «Layout».
    - Cambia: integración final; elimina las rutas v5/v6 y sus pruebas que
      afirmen admisión por metadata o por sonda, sin cambiar el algoritmo
      biológico.
    - Se comprueba con: dos contigs, adyacencia, N, máscara, fragmentación,
      menos de dos SNP, reintento, agotamiento e invariancia ante planificación.
      El fixture instrumental recorre carga, admisión, binding, streaming FASTA
      y layout; fuerza ganador interior, último tramo, dos residuos y reintento;
      desglosa capacidad/request, copias, reservas, total/prefijo, crecimiento,
      división, `remove`/`insert`, reintentos y picos de fase. Un segundo oráculo
      reconstruye `resource_preflight` v7 desde entradas y salidas, sin tipos ni
      funciones productivas.
    - Gates: `cargo fmt --check`, Clippy estricto, suite F5.1.2 dirigida, suite
      completa, validación JSON Schema de efectivo/manifiesto y búsquedas
      negativas de `f5-resource-accounting-v5`, v6, sonda 4096 y lectura
      completa no acotada en la ruta F5.

  **Pausa de cierre F5.1.2-C.** Cerrar solo con todos los gates verdes,
  reconstrucción independiente v7, evidencia de fronteras exactas y
  confirmación de cero RNG, destino o publicación antes de `sealed`.
- ☐ **F5.1.3.** Generar `p_i` por bloque×subpoblación desde el espectro, aislando
  fila y valor y conservando 0/1 sin consumo indebido.
  - Depende de: F5.1.2.
  - Procede de: REQ-EVO-40/75, ACE-89/117 y DTC-089.
  - Cambia: generador marginal y ámbitos `marginal_spectrum`.
  - Se comprueba con: extremos, intervalo degenerado, uniformes de 53 bits,
    subpoblación añadida/reordenada y perfil sin selección implícita.

> **Revisión F5.1.** Aprobar layout reconstruible, suma exacta e independencia
> de marginales antes de construir haplotipos.

## 5. F5.2 — tabla haplotípica y LD

- ☐ **F5.2.1.** Crear cada columna mediante permutación y prefijo ponderado, sin
  enumerar subconjuntos, y verificar la cota única `max(q)`.
  - Depende de: F5.1.3.
  - Procede de: REQ-EVO-40/75/76, ACE-117 y DTC-089 «Layout y marginales».
  - Cambia: Constructor F5 y matriz bit-packed.
  - Se comprueba con: pesos uniformes/no uniformes, p 0/1, empate, H>`2^S`,
    duplicados eliminados por swaps o intento rechazado, H=4096 y prueba
    instrumental de ausencia de `2^H`.
- ☐ **F5.2.2.** Calcular p, p_ij, D, `r²`, bins y signos con orden y Neumaier
  normativos.
  - Depende de: F5.2.1.
  - Procede de: REQ-EVO-43/44/75, ACE-89/117 y arquitectura §15.
  - Cambia: Evaluador LD.
  - Se comprueba con: oráculo racional pequeño, D positivo/negativo/cero,
    r² 0/1/intermedio, fijados, bin vacío y ausencia de pares interbloque.
- ☐ **F5.2.3.** Ejecutar propuestas 2×2 y seleccionar el mejor intento con la
  tupla, consumo RNG, terminación y límites exactos de DTC-089.
  - Depende de: F5.2.2.
  - Procede de: REQ-EVO-76/77, ACE-117 y DTC-089 «Búsqueda».
  - Cambia: búsqueda local y diagnóstico de intentos; el coste medido y
    declarado es `O(H*S² + Q*S)` por bloque×subpoblación×intento, no una cota
    que omita la evaluación ponderada de todos los pares.
  - Se comprueba con: propuesta no candidata, mejora, empate rechazado, éxito
    temprano, agotamiento, desempate por índice, complejidad medida y ausencia
    de relajación o falsa declaración de imposibilidad.

> **Revisión F5.2.** Revisión independiente del oráculo LD, consumo RNG,
> complejidad y mejor intento. No continuar con un resultado solo plausible.

## 6. F5.3 — fundadores intactos y F

- ☐ **F5.3.1.** Resolver para cada bloque×subpoblación el modo solicitado y
  efectivo —predeterminado de subpoblación o sustitución única por índice— y la
  tolerancia controlada declarada/efectiva; implementar muestreo probabilístico
  y cuotas por mayores restos con factibilidad previa y permutación de copias.
  - Depende de: F5.2.3.
  - Procede de: REQ-EVO-40/76, ACE-89/117 y DTC-089 «Fundadores».
  - Cambia: Inicializador fundador.
  - Se comprueba con: matriz de dos bloques×dos subpoblaciones con modos
    distintos, sustitución ausente/duplicada/fuera de rango, `null` en
    probabilístico, suma C, cota inclusiva predeterminada 0,5/C y una más
    estricta, empate reproducible por RNG, cuota imposible y F>0 incompatible.
- ☐ **F5.3.2.** Aplicar F diploide y publicar F realizado; rechazar F en otras
  ploidías y registrar Ne no aplicable en F5.
  - Depende de: F5.3.1.
  - Procede de: REQ-EVO-46–52/76, DTPD-109 resuelta y DTC-089.
  - Cambia: emparejamiento y estadísticas fundadoras.
  - Se comprueba con: fórmulas bialélica/multialélica, F=0/intermedio/1,
    H_expected=0, P=1/P>2 y vector de Bernoulli.
- ☐ **F5.3.3.** Materializar antes de toda salida las identidades de fundador y
  su pertenencia única: `IndividualId`, `FCP_` por secuencia/copia y `FHA_` por
  bloque como asignación separada.
  - Depende de: F5.3.2.
  - Procede de: REQ-EVO-32/43/53/76, REQ-POB-01, ACE-89/91/117 y arquitectura §15.
  - Cambia: repositorio de verdad, pertenencia y enlaces `FCP_`/`FHA_`→`HAP_`.
  - Se comprueba con: reconstructor de cada copia y asignación, fase/ploidía,
    identidad a través de bloques/secuencias, pertenencia ausente/duplicada e
    IDs aportados/derivados; búsqueda negativa de `HaplotypeProcess`, cruces,
    mosaicos o mapas.

> **Revisión F5.3.** Aprobar fórmulas, cuotas y verdad fundadora antes de
> persistir el resultado completo.

## 7. F5.4 — artefactos, fallo y oráculo

- ☐ **F5.4.1.** Escribir conforme a los contratos ya cerrados en F5.0
  exactamente los siete TSV de éxito y la proyección fundadora obligatoria
  `outputs/truth.jsonl` v2 mediante `f5-founder-truth-projection-v2`, con sus
  IDs, clases, claves y orden.
  - Depende de: F5.3.3.
  - Procede de: REQ-EVO-77, REQ-SAL-02, REQ-MAN-92, REQ-VAL-56 y DTC-089 «Formatos».
  - Cambia: escritores F5; no cambia ni extiende ningún esquema.
  - Se comprueba con: lector distinto y corrupción aislada de cada campo,
    clave, orden, `NA`, decimal, alelos, tamaño y SHA-256, incluido que
    `alleles`, `matrix_bytes` y `placed_layout_bytes` respetan su codificación
    canónica; desde `truth.jsonl` reconstruye por sí solo toda tupla FCP, fase,
    SNP, alelo, dosis y genotipo. Cada copia aporta `founder_membership`, cada
    genotipo `f5_locus`, y falla ante población, versión/path/hash alterados o v1.
- ☐ **F5.4.2.** Implementar el oráculo independiente desde los bytes cerrados de
  los siete TSV y `truth.jsonl`, incluido el recuento MAN-70 que no comparte
  contadores, asignaciones ni acumuladores del productor.
  - Depende de: F5.4.1.
  - Procede de: REQ-MAN-70, REQ-VAL-40/56, ACE-117 y DTC-049/089.
  - Cambia: validador de aceptación y reconstructor, no productor.
  - Se comprueba con: mutación aislada de layout, H×S/q, marginal, D/`r²`, signo,
    intento, cuota/tolerancia/modo, F, ID, asignación, RNG, truth, auditoría
    DTC-049, cota de recurso y cada vector esperado/realizado/recontado MAN-70.
    El recuento agrupa `site_id`×`population_id` leyendo solo `truth.jsonl` v2,
    sin abrir los TSV ni consultar el manifiesto;
    cada mutación debe fallar por la causa correspondiente.
- ☐ **F5.4.3.** Instanciar la configuración efectiva y el manifiesto conforme
  a los esquemas completos de F5.0, con las claves y siete enlaces TSV exactos
  —ID, clase, versión, ruta, tamaño y hash—, el enlace
  `truth` separado, las auditorías físicas DTC-049, las cotas/límites de recursos
  y los resúmenes MAN-70 producidos por F5.4.2; implementar
  `f5-failure.json` sin artefactos biológicos y la publicación atómica.
  - Depende de: F5.4.2.
  - Procede de: REQ-MAN-70/92, REQ-VAL-40/56, ACE-117 y DTC-044/045/049/089.
  - Cambia: serializadores, validación y publicación atómica; no los esquemas.
  - Se comprueba con: éxito 7/7+truth, fallo 0/7+sin truth+diagnóstico, rechazo
    preflight sin destino, interrupción `INCOMPLETE`, JSON Schema independiente,
    igualdad exacta de las tres copias `seed`, cuatro normalizaciones completas,
    `reference_bases_processed`, `tract_count`, pico de IDs FASTA retenidos,
    `S_bound`, bins, filas,
    aportaciones ajenas y tres cotas/límites, con igualdad exacta entre efectivo
    y manifiesto, y resumen MAN-70 por dosis con `recount_source` e igualdad
    obligatoria.

> **Revisión F5.4.** No continuar sin demostrar publicación atómica y oráculo
> independiente sobre éxito y fallo.

## 8. F5.5 — auditoría de identidad y publicación cerrada

- ☐ **F5.5.1.** Auditar que el paquete publicado solo lee identidad, copias,
  pertenencia y asignaciones ya cerradas en F5.3; ningún paso F5.5 puede crear,
  sustituir o modificar esos datos.
  - Depende de: F5.4.3.
  - Procede de: REQ-EVO-32/53/74/77, REQ-POB-01 y arquitectura §15.
  - Cambia: evidencia y validación, no el estado biológico ni el paquete atómico.
  - Se comprueba con: mutación posterior prohibida, pertenencia ausente/duplicada,
    IDs aportados/derivados, igualdad FCP/FHA y ausencia de migración, selección
    o genealogía.

> **Revisión F5.5.** Confirmar que identidad y pertenencia precedieron a la
> publicación y que el hito no modifica el paquete atómico.

## 9. F5.6 — Marimo y cierre

- ☐ **F5.6.1.** Añadir los controles F5 a Marimo sin panel ni recombinación y
  ejecutar el mismo `validate`/`run` Rust.
  - Depende de: F5.5.1.
  - Procede de: REQ-UI-15, ACE-114/115/117 y DTC-051/089.
  - Cambia: pestaña poblacional y adaptador de proceso.
  - Se comprueba con: matriz campo→control, importación/exportación de exactamente
    un general v3+escenario v1, rechazo de cero/varios/v2-v4, edición de
    sustituciones modo/tolerancia y presupuestos, invalidación tras editar y
    paridad de efectivo, diagnóstico y hashes con TOML manual.
- ☐ **F5.6.2.** Incorporar el perfil opcional
  `human_ld_ragsdale_approx_v1` mediante su recurso versionado y receta de medias,
  sin beta/cuantiles activos ni selección predeterminada.
  - Depende de: F5.6.1.
  - Procede de: DTC-089 «Paquete de perfiles» y decisión funcional preservada.
  - Cambia: recurso incorporado, licencia/cita y fixtures.
  - Se comprueba con: SHA-256 de fuente/recurso, reproducción de cada `mu`,
    ausencia de `kappa` en cálculo y selección exclusivamente explícita.
- ☐ **F5.6.3.** Ejecutar todos los gates, guardar evidencia y auditar las 46
  filas exactas de la matriz REQ→ACE→hito→prueba con hito F5, incluidas
  `REQ-VAL-54`, `REQ-MAN-70` y `REQ-SAL-02`, sin omitir filas compartidas ni
  incorporar filas ajenas al alcance.
  - Depende de: F5.6.2.
  - Procede de: todos los REQ/ACE F5 y reglas de avance de la hoja de ruta.
  - Cambia: evidencia de cierre y estado de trazabilidad, no contrato biológico.
  - Se comprueba con: recuento mecánico exacto 46, unicidad de los 46 REQ,
    presencia expresa de `REQ-VAL-54`/`REQ-MAN-70`/`REQ-SAL-02`,
    `cargo fmt --check`, Clippy estricto, `cargo test`,
    lectores JSON/TSV/JSONL, oráculo, vectores, paridad y búsquedas negativas de panel,
    SV, repeticiones, BAM/CRAM y recombinación bajo F5.

> **Revisión F5.6.** No iniciar F6 ni F7 hasta aprobar las evidencias y confirmar
> que no quedan decisiones funcionales o técnicas F5 abiertas.

## 10. Matriz de cierre

| Hito | Resultado | Evidencia principal | Riesgo protegido |
|---|---|---|---|
| F5.0 | contratos v3/v1/v5 y núcleo contable v7 | cardinalidad/versiones, fórmulas/fronteras y cero publicación | entrada/identidad/recurso ambiguo |
| F5.1 | admisión sellada, layout y marginales | lectores acotados, receipt físico, oráculos de cuenta/geometría y aislamiento | lectura no acotada, cruce de preflight o coordenada/objetivo falso |
| F5.2 | H×S y LD acotados | oráculo D/`r²` y complejidad | LD plausible pero incorrecto |
| F5.3 | fundadores, copias y pertenencia cerrados | fórmulas y reconstrucción | fase/F/cuota o identidad errónea |
| F5.4 | siete TSV, verdad fundadora v2 y fallo inequívoco | reconstrucción `truth` por población/locus, MAN-70 solo desde verdad, DTC-049 y JSON Schema | procedencia parcial o autorreferente |
| F5.5 | paquete inmutable auditado | mutación posterior prohibida | historia implícita |
| F5.6 | interfaz y cierre | paridad, gates y auditoría | contrato paralelo o alcance filtrado |
