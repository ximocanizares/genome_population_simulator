# Plan de implementación F5 0.1 — población fundadora y LD sintético

## 1. Estado, autoridad y alcance

- **Fecha:** 1 de septiembre de 2026.
- **Estado:** ejecutable de forma secuencial después del cierre verificable de
  F4. No autoriza publicación externa ni inicio de F6/F7. La cláusula original
  decía además «no autoriza *commit*», que ningún otro plan del proyecto
  contiene y que el historial contradice —hay commits de F5 en `main` desde su
  aprobación—: se retira por inaplicada en lugar de mantenerla como norma
  muerta.
- **Fuentes normativas:** REQ-EVO-32–34/38–41/43/44/46–53/74–78 y los REQ de
  configuración, individuo, validación, manifiesto e interfaz asignados a F5 en
  la matriz; ACE-89/90/91/104/114/115/117, DTC-089, arquitectura §15 y hoja de
  ruta §F5. La matriz asigna a F5 también REQ-EVO-01 y REQ-EVO-03, fuera del
  rango citado arriba.

F5 genera internamente bloques y haplotipos SNP A/C/G/T, ajusta marginales y
media `r²(d)` y asigna fundadores completos e intactos. No admite panel,
haplotipos externos, `alt_af` intrabloque, *indels*, repeticiones, SV, BAM/CRAM,
mosaicos, recombinación ni genealogía. F6 añade variación compleja; F7.1 diseña
la primera meiosis. F5.3 materializa identidad y pertenencia poblacional antes
de artefactos; F5.5 solo audita que el paquete resultante no las modifique.

La jerarquía única es F5.0–F5.6. Cada hito depende únicamente del anterior; no
existe dependencia de revisión futura ni ciclo F5→F6/F7→F5. Los pasos de
esquema crean versiones nuevas sin reescribir las publicadas **externamente**;
ningún esquema F5 lo está todavía, de modo que el contrato de
`resource_preflight` puede completarse dentro de v5 sin crear una versión
nueva. Esa distinción es necesaria: DTC-089 define `resource_preflight` como
objeto exacto con `input_identities` y `phases`, y exige actualizar ambos
esquemas y sus lectores antes de reanudar la implementación. No queda
ninguna decisión funcional o técnica F5 pendiente: DTPD-115 pertenece a F9 y no
es una puerta de admisión ni de ejecución de este plan.

### 1.1. Revisión adversarial del 8 de septiembre de 2026

Una revisión adversarial contrastó este plan contra DTC-089, los REQ, la
arquitectura, la matriz de trazabilidad, el código F5 ya implementado, el
historial de Git y la ejecución real: la batería completa, el test H10 aislado y
un experimento con el binario sobre el fixture de fragmentación. Registró cuatro
hallazgos graves, seis medios y cinco menores, y retiró siete sospechas al
comprobarlas.

| ID | Hallazgo | Estado |
|---|---|---|
| F5-R01 | El criterio «límite exacto admite / exacto−1 rechaza» omite la construcción de anchura fija que DTC-089 exige, de modo que seguirlo al pie de la letra produce fixtures imposibles; H10 falla exactamente por eso | **Cerrado** en F5.0.4, F5.1.2.4 y F5.1.2.6 |
| F5-R02 | El esquema efectivo v5 ya publicado no tiene `input_identities` ni `phases`, y el plan promete a la vez no reabrir esquemas y copiarlos en .5 | **Cerrado** el 9 de septiembre de 2026: DTC-089 fija la actualización atómica en F5.1.2.5, junto al sellado que produce esos objetos |
| F5-R03 | Subhitos cerrados con sus dependencias abiertas, estado del plan divergente de Git y verificación circular de F5.0.4 | **Cerrado** en F5.0.4 y en el estado de los hitos |
| F5-R04 | `phase_id`, componente del identificador `FCP_`, no está definido en ningún documento ni tiene constructor F5 | **Cerrado** el 9 de septiembre de 2026: `FCP_` se deriva de `copy_index_u32` |

Los medios y menores se corrigieron en la misma pasada: el fixture de máscaras
comparables, la contradicción entre la cláusula de *commit* y el historial, el
techo implícito que la cota de 65536 bytes impone al número de fundadores, la
afirmación inexacta sobre `input_parser_peak_bytes`, el momento en que se valida
el SHA del FASTA, el desvío real de H10, los comandos de las puertas, las
referencias ACE/REQ omitidas en la cabecera y el volumen de registro de estado
incrustado en el plan.

Lo que la revisión **retiró** tras comprobarlo también importa, porque acota
dónde no hay problema: las 46 filas de la matriz de trazabilidad son exactas y
únicas e incluyen las tres que el plan nombra; el fixture de anchura fija de
memoria existe y hace punto fijo; los catálogos RNG e ID de F5.0.3 están
completos; no se han incorporado crates nuevas; las fórmulas del fixture de
fragmentación se derivan exactamente de la fórmula v7; no queda `fs::read` en la
ruta F5 fuera de pruebas; y **el modelo de presupuesto no tiene doble cuenta ni
divergencia entre la reserva prospectiva y la fórmula final**.

## 2. Criterios de orden y puertas comunes

Primero se cierran los contratos que consumen todos los componentes. Después se
construyen layout y marginales; luego H×S/LD; después se materializan identidad
y pertenencia fundadora, y solo entonces se construye y publica el paquete
atómico. Ningún paso posterior puede modificar identidad, pertenencia, copias o
asignaciones ya incluidas en ese paquete. Los riesgos de error silencioso
—decimal de cobertura, referencia, fase, D/`r²`, reducción numérica, RNG e IDs—
se contrastan con lectores u oráculos que no comparten el productor.

Cada paso ejecuta las pruebas específicas indicadas y, antes de su pausa:
las puertas exactas de `CONTRIBUTING.md` —`cargo fmt --all -- --check`,
`cargo clippy --workspace --all-targets --locked -- -D warnings` y
`cargo test --workspace --all-targets --locked`, más `cargo test --workspace
--doc --locked`— restringiendo la última al subconjunto pertinente. El cierre
ejecuta la suite completa. No se usa `--all-features`: el proyecto no declara
ninguna.

**Techo implícito de la cota de documento (F5-R08).** DTC-089 fija un límite
físico inclusivo de 65536 bytes para el general y el escenario, y F5 exige que
`founder_ids`/`founder_keys` tengan exactamente `individual_count` elementos.
Ambas cosas juntas ponen un techo al tamaño de la población que ningún REQ
declara: con claves cortas, del orden de cuatro o cinco mil fundadores, y menos
si el escenario lleva algo más. La misma cota alcanza ya al general v2 y a la
lista `individual_ids` de F4. El techo se declara aquí como límite conocido, y
F5.0.1 lo comprueba con un documento legítimo próximo al límite; ampliarlo
exigiría reabrir DTC-089 y no pertenece a F5.

## 3. F5.0 — contratos y preflight

**Estado de los hitos de F5.0 y F5.1.1 (8 de septiembre de 2026).** Sus casillas
figuraban sin marcar mientras F5.1.2.1–.3, que dependen de ellos, figuraban
cerrados. Git registra `307995f` «cerrar contratos y preflight de F5» y
`c73a14e` «cerrar planificador de longitudes F5.1.1», pero ninguno marcó una
casilla, de modo que el plan nunca reflejó ese cierre. La casilla de cada hito
se conservó abierta a propósito hasta que la revisión independiente de F5.1.2
confirmara que la sustitución v7 no los reabrió; lo que no puede repetirse es
cerrar un subhito dependiente sin dejar constancia del estado de aquello de lo
que depende.

Las cinco casillas se cerraron el 10 de septiembre de 2026, al cerrarse F5.1.2.
La condición se cumplió: las revisiones independientes de F5.1.2.4, .5 y .6
recorrieron el preflight puro, el ledger, los catálogos y los esquemas —que es
el contenido de F5.0— sin hallar en ellos ningún defecto, y el contrato de
`resource_preflight` que F5.0.1 dejó incompleto se cerró en F5.1.2.5. El vector
`FCP_` del catálogo de F5.0.3 se regeneró el 9 de septiembre al derivar el
identificador del índice de copia.

- ✅ **F5.0.1.** Crear y cerrar los contratos externos completos general v3,
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
    `truth_projection="f5-founder-truth-projection-v2"`. Ningún hito posterior
    amplía estos esquemas, con una excepción declarada: el objeto
    `resource_preflight` se cierra en F5.1.2.5 con `input_identities` y
    `phases`, que DTC-089 exige y que la proyección escalar de F5.1.2.3 dejó
    fuera. Hasta entonces el esquema efectivo v5 publica solo los 27 campos
    escalares y su `additionalProperties: false` los rechazaría, así que
    completarlo es trabajo de .5 y no una versión nueva. DTC-089 lo fija así
    desde el 9 de septiembre de 2026: la actualización atómica de esquemas,
    lectores y productor ocurre en F5.1.2.5, que es cuando el sellado produce
    las identidades y las fases.
  - Se comprueba con: pareja mínima válida; rechazo de cero, lista o varios
    escenarios; rechazo de v2/v4 como F5; inmutabilidad byte a byte de los tres
    esquemas heredados; campos ausentes/desconocidos;
    one-of de objetivo, máscara/hash, suma de fundadores, IDs, F/P y límites;
    CLI relativa, escenario relativo al general, recurso relativo al escenario,
    duplicación prohibida; `seed` exacta `^[0-9a-f]{64}$` e igualdad entre
    general resuelto, raíz efectiva, `f5.seed` efectivo y `f5.seed` manifiesto;
    discordancia de cualquier otra copia efectiva/manifiesto; y verdad v2 con
    población/copia/locus estructurados. Los esquemas legacy se comparan byte a byte.
- ✅ **F5.0.2.** Implementar el lector único de `profile.metadata.json` y los
  cinco TSV, con sus cabeceras, dominios, orden, claves, normalización y hashes.
  - Depende de: F5.0.1.
  - Procede de: REQ-EVO-74/75, ACE-117 y DTC-089 «Paquete de perfiles».
  - Cambia: recursos F5 y validación de entrada; no genera datos biológicos.
  - Se comprueba con: corpus válido y una mutación aislada de BOM, LF, cabecera,
    fila, intervalo, peso, población, versión, tamaño y SHA-256.
- ✅ **F5.0.3.** Publicar los catálogos inmutables `f5-rng-catalog-v1` e ID F5,
  incluidos material tipado, prefijos y vectores por ámbito.
  - Depende de: F5.0.2.
  - Procede de: REQ-EVO-77, REQ-VAL-56 y DTC-037/039/089.
  - Cambia: catálogos RNG/ID y fixtures, incluidas `FCP_` por secuencia/fase y
    `FHA_` como asignación de bloque distinta; no cambia los catálogos previos.
  - Se comprueba con: derivador independiente BLAKE3/ChaCha12, enteros con
    rechazo, Fisher–Yates, Bernoulli 0/1 sin consumo, colisión dirigida y bytes
    canónicos exactos de layout (`u64` big-endian) y matriz (cabecera H/S,
    orden de filas, bits MSB y relleno cero).
- ✅ **F5.0.4.** Implementar el preflight puro
  `f5-resource-accounting-v8` y su ledger acumulativo, que resuelve racional de
  cobertura, capacidad, dimensiones y cotas de operaciones, pico de memoria y
  bytes publicados con símbolos dimensionalmente distintos. Este núcleo no
  abre ficheros: recibe observaciones de las fases físicamente acotadas de
  F5.1.2 y no acepta tamaños obtenidos por `metadata` como autoridad.
  - Depende de: F5.0.3.
  - Procede de: REQ-EVO-74/76/77, ACE-117 y arquitectura §15.
  - Cambia: admisión, escáner secuencial referencia+máscara, lista acotada de
    tramos elegibles y catálogo diagnóstico F5.
  - El hito se cierra con las fórmulas y el ledger, que no abren ficheros y son
    comprobables con observaciones inyectadas. Los fixtures que exigen recorrer
    físicamente un FASTA o una máscara —el de 80 MiB de descripción, el de
    `region_id` largo y el de 4096 tramos— se verifican en F5.1.2.4, que es
    donde existe el escáner. Sin esa separación la comprobación sería circular:
    F5.1.2.1 depende de este núcleo y este núcleo se verificaría con el lector
    de F5.1.2.4.
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
    Toda frontera «límite exacto admite, límite menos uno rechaza» se construye
    con el literal de anchura fija y punto fijo de DTC-089, y en las tres
    dimensiones, no solo en memoria de trabajo: los bytes del escenario entran
    en la cota de operaciones como `2*scenario_toml_bytes`, de modo que escribir
    el límite después de calcular la cota la cambia y hace el fixture
    irrealizable. El documento base ya contiene la línea del límite con su
    anchura definitiva; el punto fijo solo sustituye sus dígitos.
    Dos FASTA con la misma geometría y bases pero descripciones mínima y de
    80 MiB aumentan la cota exactamente por sus bytes recorridos en preflight y
    recuperación REF; el límite igual se admite y el límite menos uno se
    rechaza antes de RNG. La misma prueba se aplica a una máscara con igual
    geometría y distinto tamaño de `region_id`.
    Un fixture usa el mismo escenario con `block_count=1`, el mismo FASTA de
    una secuencia y dos máscaras con igual E: una deja un único tramo compacto
    y la otra 4096 tramos de longitud `Lmin`, separados por una base
    enmascarada. Las dos máscaras se igualan en bytes rellenando `region_id`,
    porque 4095 filas frente a una cambian `mask_bytes_processed` y sin ese
    relleno «iguales todos los demás valores» es falso. Debe reconstruir
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

- ✅ **F5.1.1.** Construir longitudes exactas condicionadas para
  `covered_bases`/`coverage_fraction` y estocásticas para `block_count`.
  - Depende de: F5.0.4.
  - Procede de: REQ-EVO-74, ACE-117 y DTC-089 «Layout».
  - Cambia: geometría F5.
  - Se comprueba con: suma T exacta, residual final, pesos restringidos,
    overflow y un oráculo entero que prueba la factibilidad previa.
- ✅ **F5.1.2.** Integrar la admisión F5 físicamente acotada y, solo tras
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
      interno `input_parser_peak_bytes`: se retiró de la proyección publicada
              y del esquema
      efectivo —y, por referencia, del manifiesto— y se añadió un oráculo de
      claves literal e independiente para ambas rutas. Esta proyección de .3
      no se presenta como el contrato final: `input_identities` y `phases`
      siguen reservados a .5. La revisión independiente de cierre del 5 de
      septiembre de 2026 aprobó la pausa A: verificó los 27 campos escalares
      publicables en .3, la reserva prospectiva de TSV, la dimensión
      diagnóstica y las rutas públicas. H10 sigue reservado a F5.1.2.6 y no se
      presenta como gate verde.

  > Los bloques de estado fechados de F5.1.2.3 y F5.1.2.4 se conservan porque
  > documentan reaperturas que cambiaron el alcance, no avances rutinarios. El
  > detalle contable vive en `evidence/f5.1.2/`; el plan no lo duplica.

  **✅ Pausa de revisión F5.1.2-A.** Aprobar dispatcher, bootstrap, lector, perfil
  y ledger. Exigir ausencia mecánica de `fs::read` y de sondas parciales en la
  ruta F5, pruebas exacto/exacto−1 y cero apertura posterior. No continuar con
  referencia hasta corregir todos los hallazgos.
  Esta pausa se aprobó con el target `--lib` en rojo por H10, contra la puerta
  común de §2. Se acepta como excepción declarada y acotada, no como práctica:
  H10 vive en `src/f5_layout.rs`, pertenece a este mismo hito y su causa quedó
  establecida (F5-R01). Ninguna pausa posterior se aprueba con un target en
  rojo; F5.1.2.6 lo deja verde reescribiendo el fixture.

  - ✅ **F5.1.2.4.** Integrar referencia y máscara en el mismo ledger.
    - Depende de: revisión F5.1.2-A aprobada.
    - Procede de: REQ-EVO-74/76, DTC-040 y DTC-089 «Layout».
    - Cambia: primer FASTA y única máscara en streaming; reserva de dos visitas
      por byte FASTA antes de leerlo; comprobación previa a retener IDs, filas y
      tramos; receipt de tamaño/SHA observado en el mismo handle.
    - El SHA del FASTA se observa **dentro** del escaneo y bajo el ledger, no
      antes del preflight: donde DTC-089 dice «antes del preflight de recursos»
      describe el flujo anterior a v7, superado por el propio DTC-089 al
      introducir el ledger. El escaneo es el único lector físico.
    - Se comprueba con: descripciones FASTA y `region_id` largos, línea de
      máscara excesiva, referencia/máscara vacía o truncada, límite exacto y
      exacto−1 en operaciones y memoria —siempre con el literal de anchura fija
      y punto fijo de F5.0.4, porque los bytes del escenario entran en la cota
      de operaciones—, SHA incorrecto y espías que prueban cero RNG/destino. El oráculo recuenta bytes, bases y tramos sin reutilizar
      estructuras productivas.
    - Estado en corrección (5 de septiembre de 2026): C1–C6 están aplicados
      localmente con regresiones focales y gates de formato, Clippy y
      documentación correctos. La suite completa conserva únicamente H10
      (`4961/4952`, +9), aplazado expresamente a F5.1.2.6 como C7. La
      contabilidad completa consta en
      `evidence/f5.1.2/aplicacion-correcciones-f5.1.2.4-2026-09-05.md`.
      F5.1.2.4 permanece abierto hasta la revisión independiente; no existe
      todavía receipt sellado ni segunda visita física FASTA.
    - Corrección residual (5 de septiembre de 2026): una revisión
      independiente posterior reabrió C3/C5 mediante L1–L3: `other` podía
      sustituirse después del escaneo, `covered_bases` reservaba con B=0 y la
      geometría se clonaba antes de admitir su coexistencia. Los tres están
      corregidos localmente con binding inicial, reserva conocida temprana y
      backing único compartido. El checkbox continúa abierto hasta una nueva
      revisión independiente; no se ha avanzado .5/.6 ni H10.
    - Corrección residual final (5 de septiembre de 2026): C8 elimina de la
      magnitud publicada las sondas EOF no observadas y las respalda con
      trabajo futuro ya cargado, sin permitir lecturas físicas sin reserva;
      un oráculo literal independiente, con y sin máscara, acepta la fórmula
      exacta de DTC-089 y rechaza exacta−1. C9 clasifica un SHA de máscara
      distinto como `GPS-INPUT-002`, con esperado/observado, mientras la
      sintaxis malformada conserva `GPS-FMT-101`; ambas rutas rechazan antes de
      entropía o destino. F5.1.2.4 sigue sin marcar hasta otra revisión
      independiente. No se ha modificado la visibilidad pública provisional,
      que podrá reevaluarse al sustituirla por el sellado interno de .5.
    - Causa de H10, establecida el 8 de septiembre de 2026 (F5-R01): **no es un
      defecto del modelo de presupuesto**. El fixture calcula la cota sin la
      línea del límite y después añade `max_operation_units = 4952`, cuyos
      27 bytes entran dos veces en la cota; el exceso observado es exactamente
      2×27 = 54, y declarar 5006 = 4952 + 54 admite el mismo fixture. Se
      comprobó con el binario sobre el fixture, variando solo el literal a
      anchura constante: sin línea admite, 4952 y 5005 rechazan con acumulado
      5006, y 5006 admite. La reserva prospectiva del escaneo coincide con la
      fórmula final, de modo que no hay doble cuenta ni divergencia.
      Por tanto la aserción `preflight.operation_units() ==
      baseline.operation_units()` del fixture es insatisfacible por
      construcción en v7, y H10 no puede cerrarse ajustando una constante:
      exige reescribir el fixture con la línea de anchura fija ya presente en
      el documento base, como hace
      `f5_bootstrap_dispatch::fixed_width_bootstrap_memory_limit_…` para
      memoria.

  - ✅ **F5.1.2.5.** Sellar `F5AdmissionReceipt`, exigir
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

    - Revisión independiente de cierre (9 de septiembre de 2026). Verificó C1–C9
      y L1–L3 uno a uno contra el código, no contra la evidencia: todos
      aplicados, ninguno reaparecido, ningún hallazgo grave. Confirmó además que
      las cuatro fronteras exacto/exacto−1 de operaciones se construyen con
      punto fijo o marcador de anchura fija, y que las dos que no lo hacían
      acreditaban la frontera real por coincidencia; se les añadió la guarda.
      El criterio «referencia o máscara vacía o truncada» carecía de regresión
      en esta ruta —la cobertura de F1 ejercita otro lector— y se añadieron los
      casos de referencia vacía, con solo cabecera y truncada, y de máscara
      vacía y con la última fila sin LF. Se completaron también la frontera
      inclusiva exacta de la fila de máscara en 65536 bytes y la comprobación de
      destino, que no pedía ninguno y por tanto no acreditaba nada.
    - Deuda conocida y acotada: `EligibleScan::ensure_memory_bound` es
      inalcanzable, porque el ledger compara el mismo valor más la retención de
      perfil y rechaza antes, y su diagnóstico no llevaría `dimension`, `phase`,
      `role`, `accumulated` ni `limit`. No se retira porque en la ruta de IDs
      adelanta la comprobación antes de `SequenceId::parse`, que es la
      «comprobación previa a retener» que este hito exige; delegarla en el
      ledger exigiría un método de consulta sin cobro que hoy no existe. Se
      revisará al sellar el receipt en F5.1.2.5.
    - Excepción declarada al bullet «reserva de dos visitas por byte FASTA antes
      de leerlo»: el escrow de C8 permite entregar como máximo un byte de FASTA
      y uno de máscara contra trabajo futuro ya cargado. Es de un solo uso, la
      lectura siguiente rechaza, y el ledger exige que el escrow no supere los
      totales. Sin esa excepción la magnitud publicada tendría que incluir las
      sondas EOF que C8 retiró.

  **✅ Pausa de revisión F5.1.2-B.** Aprobar el modelo local de cambios
  accidentales, el orden de fases, las identidades físicas, la reserva de la
  segunda pasada y el binding sellado. No continuar con layout si un consumidor
  todavía acepta escalares o receipts sueltos.
  Aprobada el 10 de septiembre de 2026. `SealedF5Admission` es el único tipo
  público y ningún consumidor acepta ya escalares ni receipts sueltos: no hay
  constructor desde documento serializado, de modo que un receipt v5 o v6 no
  puede migrarse. La revisión de código del hito no halló ningún defecto grave;
  sus cinco hallazgos menores se corrigieron antes de cerrar, incluido el único
  con consecuencia real: `F5BlockLengthPlan::sample` era público y, como
  `from_preflight` solo exige la mitad semántica del binding, existía una cadena
  pública capaz de abrir un ámbito RNG con una configuración físicamente
  distinta de la admitida.
  - Queda declarado para F5.1.2.6: no existe todavía productor de
    `effective-config.json` ni de `manifest.json` v5 para F5, de modo que la
    igualdad exigida entre ambos está acreditada por pruebas pero no ejercida
    por código productivo. El comparador pertenece a ese hito.

  - ✅ **F5.1.2.6.** Revalidar el layout productivo y cerrar la sustitución v7.
    - Depende de: revisión F5.1.2-B aprobada.
    - Procede de: REQ-EVO-74/77 y DTC-089 «Layout».
    - Cambia: integración final; elimina las rutas v5/v6 y sus pruebas que
      afirmen admisión por metadata o por sonda, sin cambiar el algoritmo
      biológico. Reescribe el fixture instrumental de H10 con la línea del
      límite de anchura fija presente en el documento base, que es lo que ese
      test necesita para ser satisfacible; su aserción actual de igualdad entre
      la cota con y sin límite declarado no puede cumplirse en v7.
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
    - Estado implementado (10 de septiembre de 2026), a la espera de la pausa
      F5.1.2-C. H10 quedó verde reescribiendo el fixture con la línea del
      límite siempre presente y anchura fija: la sonda `9999` mide la cota
      sobre el mismo documento, la cota resulta `5006 = 4952 + 2*27` y la
      variante exacta−1 conserva la anchura. La igualdad insatisfacible se
      sustituyó por el punto fijo `operation_units == max_operation_units`, y
      el rechazo exacto−1 se comprueba donde v7 lo emite: dentro de la carga,
      no al recuperar la admisión sellada. El fixture añade además el desglose
      de picos de fase sobre el receipt sellado. Un segundo oráculo
      independiente reconstruye la proyección v7 completa —veintisiete
      escalares, `algorithm` incluido, más nueve identidades, diez con máscara,
      y nueve fases— con y sin máscara. Una prueba de mutación lo acredita:
      perturbar un solo campo del desglose por fase hace fallar exactamente ese
      oráculo y ningún otro test del repositorio.
      El comparador que .5 dejó pendiente existe ya como código productivo:
      `SealedF5Admission::ensure_contract_copies_match` rechaza con
      `GPS-FMT-006` cualquier copia de `/f5/resource_preflight` que difiera de
      la admisión sellada, citando el puntero JSON divergente. Queda para
      F5.4.3 el productor de efectivo/manifiesto que lo invoque: hoy el
      comparador no tiene llamador productivo porque ese productor no existe.

  **✅ Pausa de cierre F5.1.2-C.** Cerrar solo con todos los gates verdes,
  reconstrucción independiente v7, evidencia de fronteras exactas y
  confirmación de cero RNG, destino o publicación antes de `sealed`.
  Aprobada el 10 de septiembre de 2026, con las cinco puertas en verde y la
  batería completa sin ningún fallo: 718 pasan, 0 fallan. La reconstrucción
  independiente v7 existe y es discriminante —una mutación de un solo campo del
  desglose por fase hace fallar ese oráculo y ningún otro test—; las fronteras
  exactas están acreditadas con literal de anchura fija en las tres dimensiones;
  y ninguna vía pública abre RNG, destino o publicación antes del sellado, tras
  cerrar en F5.1.2.5 la única que quedaba.
  - H10 estuvo rojo desde el 5 de septiembre. No solo fallaba: abortaba en la
    admisión, de modo que sus aserciones de layout nunca llegaban a ejecutarse.
    Conviene recordarlo cuando un plan aplace un test a un hito posterior: lo
    que se aplaza no es solo el fallo, sino toda la cobertura que ese test
    afirmaba dar.
- ✅ **F5.1.3.** Generar `p_i` por bloque×subpoblación desde el espectro, aislando
  fila y valor y conservando 0/1 sin consumo indebido.
  - Depende de: F5.1.2.
  - Procede de: REQ-EVO-40/75, ACE-89/117 y DTC-089.
  - Cambia: generador marginal y ámbitos `marginal_spectrum`.
  - Se comprueba con: extremos, intervalo degenerado, uniformes de 53 bits,
    subpoblación añadida/reordenada y perfil sin selección implícita.

> **✅ Revisión F5.1.** Aprobar layout reconstruible, suma exacta e independencia
> de marginales antes de construir haplotipos.
> Aprobada el 10 de septiembre de 2026. El layout es reconstruible y su
> pertenencia se comprueba: la frontera de marginales rederiva cada `block_id` y
> `site_id` y exige además el recuento sellado, porque dos ejecuciones que solo
> difieren en el número de bloques comparten el identificador del bloque cero.
> Las sumas usan el acumulador compensado en el orden escrito, y ahora hay una
> prueba con pesos no diádicos que lo distingue de una suma ingenua. Los
> marginales son independientes por subpoblación: cada sitio abre su propia
> subsemilla, reordenar las subpoblaciones no cambia un valor y añadir una no
> altera las anteriores **con la misma identidad de ejecución**. Ese matiz
> importa: `block_id` y `site_id` incorporan `execution_id`, de modo que si la
> política de F5 lo deriva de los bytes del escenario —como hacen F1 a F4 con
> el general—, añadir una subpoblación cambiaría los identificadores y con
> ellos todos los marginales, sin alterar la geometría. La prueba que acredita
> la invariancia fija la ejecución, y el material del `execution_id` de F5 se
> decide al encadenar la cadena productiva.
> - Queda declarado para F5.2.1: `build_f5_marginals` no tiene aún llamador
>   productivo, y es ese hito quien debe encadenar layout, marginales y
>   haplotipos.

## 5. F5.2 — tabla haplotípica y LD

- ✅ **F5.2.1.** Crear cada columna mediante permutación y prefijo ponderado, sin
  enumerar subconjuntos, y verificar la cota única `max(q)`.
  - Depende de: F5.1.3.
  - Procede de: REQ-EVO-40/75/76, ACE-117 y DTC-089 «Layout y marginales».
  - Cambia: Constructor F5 y matriz bit-packed.
  - Se comprueba con: pesos uniformes/no uniformes, p 0/1, empate, H>`2^S`,
    duplicados eliminados por swaps o intento rechazado, H=4096 y prueba
    instrumental de ausencia de `2^H`.
  - Estado implementado (10 de septiembre de 2026), a la espera de revisión.
    `src/f5_haplotypes.rs` encadena por primera vez layout, marginales y
    haplotipos: `build_f5_haplotype_attempt` llama a `build_f5_marginals`,
    comprueba `H<=2^V` en todos los bloques antes de abrir ningún intento y
    construye cada columna con permutación y prefijo compensado. La matriz
    bit-packed almacena exactamente los bytes canónicos de `matrix_bytes`. Los
    ámbitos `initial_column` y `polarity` reproducen su vector publicado y el
    vector se cierra hasta el valor producido. Las filas duplicadas se cuentan
    y se publican sin fusionarlas: DTC-089 encarga su eliminación a los swaps
    de F5.2.3, de modo que este hito no rechaza el intento por sí mismo. La
    ausencia de `2^H` la acreditan dos hechos independientes del algoritmo: la
    prueba con H=4096 **termina** sobre la ruta productiva —18,6 ms en debug,
    1,8 ms en release—, cosa que una enumeración de `2^4096` no haría nunca; y
    el espía de palabras mide H−1 enteros exactos por permutación. El contador
    de operaciones por bloque (469, 947, 1898, 3802 y 7610 para H=16..256 con
    S=12, y 121 866 con H=4096) confirma que el coste sigue al modelo lineal,
    pero **no demuestra** la ausencia de enumeración: solo se incrementa dentro
    de bucles acotados por H, así que no podría crecer exponencialmente aunque
    otra parte del código enumerase subconjuntos.
- ✅ **F5.2.2.** Calcular p, p_ij, D, `r²`, bins y signos con orden y Neumaier
  normativos.
  - Depende de: F5.2.1.
  - Procede de: REQ-EVO-43/44/75, ACE-89/117 y arquitectura §15.
  - Cambia: Evaluador LD.
  - Se comprueba con: oráculo racional pequeño, D positivo/negativo/cero,
    r² 0/1/intermedio, fijados, bin vacío y ausencia de pares interbloque.
  - Estado implementado y corregido tras revisión (11 de septiembre de 2026).
    `src/f5_ld.rs` evalúa un intento completo contra el layout y el perfil LD:
    recalcula `p_hat_i` desde la matriz —y exige que coincida bit a bit con el
    que publicó la columna—, suma `p_hat_ij` en orden canónico `h` con el
    `NeumaierSum` reutilizado de `f5_marginals`, y escribe `D` y `r²` en el
    orden literal de DTC-089 sin `mul_add`. La agrupación del denominador es la
    de los corchetes del contrato, `[p_i(1-p_i)] * [p_j(1-p_j)]`, y queda
    fijada: la alternativa `((a*b)*c)*d` cambia el último bit en una tercera
    parte de las parejas. `not_evaluable` es de pareja y `not_applicable` de
    bin; ninguno se colapsa en un cero. El oráculo racional es exhaustivo sobre
    las 256 matrices `4x2`, con pesos uniformes `1/4` y también no uniformes
    `4/8, 2/8, 1/8, 1/8`, donde `p`, `p_ij`, `D` y `r²` son cocientes exactos
    de enteros de 128 bits.
  - Dos decisiones locales de aritmética, ambas cerradas en revisión:
    - **La evaluabilidad se decide por factor, no por el producto.**
      `normalize_slice` solo exige `|suma-1| <= 1e-9` *después* de dividir, así
      que un perfil legal puede dejar pesos cuya suma compensada valga
      `1.0000000000000002`; una columna monomórfica realiza entonces esa
      marginal y su factor `p(1-p)` es negativo. Con dos columnas así el
      producto de dos factores negativos es positivo y declararía evaluable una
      pareja sin variación, con `r²` cercano a 1 y el signo incumplido. Exigir
      `factor > 0` en cada lado cubre a la vez el cero del contrato y ese caso.
      Una marginal fuera de `[0,1]` **no aborta**: DTC-089 no contempla abortar
      ahí, y la primera implementación lo hacía.
    - **`r²` se recorta a 1 donde nace.** En el acoplamiento perfecto
      `D = fl(p - fl(p*p))` y el factor es `fl(p * fl(1-p))`: dos redondeos
      distintos del mismo real, cuyo cociente supera 1 en una fracción
      apreciable de los pesos decimales —máximo observado `1.0000000000442435`
      con `w = 0,999999`—. Es la reconciliación entre la fórmula literal de
      DTC-089 y el dominio `[0,1]` que `f5-ld-evaluation.tsv` exige a medias y
      errores. Un solo `min` en el origen mantiene dentro del dominio `r²`, la
      media del bin y su error, sin recortarlos por separado.
  - Queda declarado para F5.2.3, y son condiciones de su diseño:
    - `evaluate_f5_linkage_disequilibrium` **no tiene aún llamador productivo**,
      igual que le ocurrió a `build_f5_marginals` antes de F5.2.1: es F5.2.3
      quien debe encadenar layout, marginales, haplotipos y evaluación.
    - El valor normativo de `p_hat_i` es la **suma canónica completa** sobre las
      H filas, y el evaluador exige que la columna publicada coincida con ella
      bit a bit. Las actualizaciones incrementales `O(S)` de las propuestas 2×2
      sirven solo para **puntuar**; antes de volver a evaluar hay que recomputar
      canónicamente la columna publicada. Si no, cada propuesta aceptada
      abortará con `GPS-INTERNAL-101`.
    - Un bloque sin ninguna pareja evaluable publica `bins_outside=0`,
      `signs_failed=0`, `max_r2_error=0` y `sum_r2_error=0`: el mejor valor
      posible en cinco de los siete campos de la tupla, de modo que un intento
      degenerado ganaría la selección. F5.2.3 **debe consultar `is_valid()`
      antes de la tupla**, tanto para terminar como para seleccionar.
    - El veredicto de intento inválido se publica en `is_valid()` y no se
      convierte aquí en `GPS-DOMAIN-102`/`104`: esa decisión pertenece a la
      búsqueda, que es quien conoce todos los intentos.
- ✅ **F5.2.3.** Ejecutar propuestas 2×2 y seleccionar el mejor intento con la
  tupla, consumo RNG, terminación y límites exactos de DTC-089.
  - Depende de: F5.2.2.
  - Procede de: REQ-EVO-76/77, ACE-117 y DTC-089 «Búsqueda».
  - Cambia: búsqueda local y diagnóstico de intentos; el coste medido y
    declarado es `O(H*S² + Q*S)` por bloque×subpoblación×intento, no una cota
    que omita la evaluación ponderada de todos los pares.
  - Se comprueba con: propuesta no candidata, mejora, empate rechazado, éxito
    temprano, agotamiento, desempate por índice, complejidad medida y ausencia
    de relajación o falsa declaración de imposibilidad.
  - Estado implementado y corregido tras revisión adversarial (11 de septiembre
    de 2026). `src/f5_search.rs` es el primer llamador productivo de la cadena
    completa: `search_f5_haplotypes` recibe configuración, admisión sellada,
    identidad de ejecución y layout, deriva los marginales una sola vez,
    **resuelve cada bloque×subpoblación por separado** —construyendo cada
    intento con `build_f5_block_attempt`, refinándolo con las propuestas 2×2 y
    evaluándolo con `evaluate_f5_block_linkage_disequilibrium`— y compone el
    resultado con el ganador de cada unidad. Publica una fila por
    bloque×subpoblación×intento con los campos de `f5-attempts.tsv` y, al
    agotar, un `GPS-DOMAIN-104` que conserva el historial y el mejor intento.
  - Decisiones locales registradas:
    - **El alcance de la tupla y el de los intentos es el bloque×subpoblación.**
      Los siete campos son magnitudes de una matriz, y los intentos se cuentan
      igual: DTC-089 escribe el coste «por bloque×subpoblación×intento»,
      `f5-attempts.tsv` se indexa por `(population_id, block_id,
      attempt_index)` con «una sola fila `true` por bloque×población», el
      contexto de `GPS-DOMAIN-104` nombra población y bloque, REQ-EVO-75 dice
      que «solo falla el bloque×subpoblación sin ninguna pareja evaluable» y
      `attempt_summary` publica `haplotype_attempts_total` y `max_attempt_index`
      **como escalares distintos**, cosa que solo tiene sentido si cada unidad
      lleva su cuenta. La primera versión hacía el intento global y todo-o-nada:
      con `U` unidades que terminan de forma independiente con probabilidad `p`,
      el intento global triunfa con `p^U`, de modo que un escenario resoluble
      unidad a unidad se declaraba agotado —16 unidades que terminaban por su
      cuenta como tarde en el intento 3 necesitaban ~320 intentos globales—. Eso
      es justo lo que DTC-089 prohíbe al decir que el agotamiento «no afirma
      imposibilidad». La independencia entre unidades está acreditada: DTC-089
      prohíbe parejas LD interbloque, el layout las sitúa en intervalos
      disjuntos, los tres ámbitos RNG llevan ya población, bloque e intento, y
      `f5-founders.tsv` asigna un haplotipo por bloque por separado. Ningún
      documento exige coherencia entre bloques en esta fase.
    - **Dos aritméticas, una sola verdad, y una guarda que atribuye cada
      divergencia.** DTC-089 fija la actualización de una propuesta en `O(S)`
      —«las actualizaciones 2×2 aceptadas mantienen las estadísticas en
      `O(S)`»— y la fórmula de admisión contabiliza ese coste exacto, de modo
      que las marginales y las conjuntas se mueven por incremento y **la
      puntuación del bucle es determinista pero no canónica**. Solo los valores
      publicados lo son: al cerrar cada intento se recomputan canónicamente las
      columnas (`refresh_realized_columns`) y se reevalúa la unidad.
      La igualdad bit a bit entre las dos aritméticas no puede exigirse
      —`m - q_a + q_b` no es la suma de Neumaier sobre el nuevo conjunto, y
      hacer que lo fuera costaría `O(H)` por propuesta—, pero sí puede exigirse
      que ninguna divergencia quede **sin atribuir**. Las tres reglas discretas
      del contrato son discontinuas en un punto exacto —la de signos en
      `D == 0.0`, la de evaluabilidad en `factor == 0.0`, la de bins en
      `error == max_abs_error`—, así que `ensure_score_did_not_diverge`
      contrasta las dos aritméticas **pareja a pareja** y **bin a bin** y solo
      exime una divergencia si uno de sus dos valores está justo en la
      frontera. Cualquier otra, y cualquier diferencia en `filas_duplicadas`
      —que no cuelga de ninguna frontera de coma flotante—, es
      `GPS-INTERNAL-101`. Los veredictos canónicos que la guarda compara se
      recogen **durante** la evaluación canónica que ya se hace de todos modos,
      de modo que no cuestan ninguna visita de fila adicional.
      La medida que motivó la guarda queda registrada: con pesos
      `0,3/0,2/0,3/0,2`, tres marginales incrementales quedaban en
      `0,5000000000000001` frente a `0,5` y seis de quince parejas cambiaban de
      veredicto de signo, todas sobre una `D` canónicamente nula. Con la
      atribución por pareja esas seis quedan explicadas y ninguna otra
      divergencia pasa desapercibida.
    - **La validez y la terminación son claves anteriores a la tupla.** Un
      bloque sin parejas evaluables tiene el mejor valor posible en cuatro de
      los siete campos, así que la tupla sola premiaría un intento degenerado.
      `block_has_terminated` consulta `is_valid()` primero y `order_key_is_better`
      antepone terminación y validez, con el menor `attempt_index` como último
      desempate.
    - **`matrix_bytes` se compara byte a byte, dentro y fuera del bucle.**
      Dentro, la regla `O(1)` del bit `(a,i)`: con `a<b` e `i<j` es el primer
      bit que el intercambio cambia, de modo que la codificación nueva es menor
      exactamente cuando ese bit valía 1. Fuera, entre intentos de una unidad,
      `F5AlleleMatrix::canonical_cmp` compara la misma codificación sin
      materializarla, en `O(H*S/8)` por intento. La primera versión comparaba el
      SHA-256, que invierte el orden lexicográfico en la mitad de las parejas e
      **incumplía** el «mínimo lexicográfico» del contrato; con el bucle dentro
      del bloque solo hay que conservar en línea la matriz del mejor candidato
      de esa unidad —la que ya hace falta para componer la salida—, así que
      cumplirlo sale gratis. El digest se sigue publicando porque el esquema de
      fallo lo exige, pero no decide nada.
    - **`selected_best` significa una sola cosa.** Con el bucle por unidad, la
      terminación es la primera clave de orden y el bucle para en el primer
      intento que termina, de modo que el mínimo de la clave y la matriz
      publicada son siempre el mismo intento. La ambigüedad que tenía la
      primera versión —una cosa al tener éxito y otra al agotarse— se disuelve.
    - **`GPS-RES-101` no se duplica aquí.** La admisión sellada, que la
      frontera exige y verifica, ya rechaza `H*site_bound_per_block >
      max_haplotype_cells`, cota que domina el producto exacto del layout.
    - **`GPS-DOMAIN-102` es una demostración y se comprueba antes que la
      cardinalidad.** Se emite sobre los marginales cuando un bloque tiene menos
      de dos columnas no fijadas, y basta uno para detener el escenario
      (REQ-EVO-75). Va **antes** de `H <= 2^V` a propósito: un bloque sin
      columnas variables incumple las dos cosas, y de los dos diagnósticos el
      que describe la causa es que no hay variación que medir, no que falte
      cardinalidad.
    - **`GPS-DOMAIN-103` no lo emite ninguna ruta de este hito**: la cuota
      controlada pertenece a F5.3 y su diagnóstico nacerá allí. La ausencia es
      deliberada, no un olvido.
    - **El mapa de duplicados usa FNV-1a fijado y contrasta bytes**, y los
      decrementos de los contadores usan `checked_sub`: un desajuste futuro
      envolvería en release y dejaría `evaluable_pairs` gigante con
      `is_valid()` siempre cierto.
    - **La tolerancia de bin se escribe una sola vez**, en `f5_bin_abs_error` y
      `f5_bin_is_outside` de `src/f5_ld.rs`, con una prueba que fija su carácter
      inclusivo en el borde exacto y un ULP por encima.
  - Coste medido, por bloque×subpoblación×intento y con `B` bins del perfil:
    la siembra canónica cuesta `H*S + H*S*(S-1)/2` visitas de fila; cada
    propuesta candidata cuesta dos marginales movidas en `O(1)`, `2S-3` parejas
    actualizadas, **solo los bins que esas parejas tocan**, `S` marginales
    recorridas y dos filas rehasheadas; una propuesta no candidata cuesta
    cuatro enteros acotados y nada más; el cierre cuesta `O(H*S)` para
    recomputar las columnas, `O(H*S²)` para la evaluación LD —que de paso
    recoge los veredictos de la guarda—, `O(S²+B)` para la atribución y
    `O(H*S/8)` para el desempate lexicográfico. El total es
    `O(H*S² + Q*S + A*B)`, con `A` los intentos: el término `A*B` es la siembra
    de los agregados de bin, que se hace una vez por intento y **no** dentro del
    bucle. El bucle cumple la cota `O(Q*S)` de DTC-089.
    La corrección que lo hizo cierto llegó tras la revisión de cierre: el bucle
    respaldaba los dos vectores de bin **enteros** por propuesta y rehacía los
    tres agregados sobre los `B` bins, de modo que el coste crecía como
    `B/S_bound` sin cota. Con un perfil legal de resolución por par de bases
    —REQ-EVO-75 exige bins contiguos de 1 a `Lmax-1`— y un `gap_min` grande
    frente a `Lmax`, un bloque de dos sitios con 255 bins llegó a costar 62
    veces el presupuesto sellado de la admisión. Ahora `bins_outside`,
    `sum_r2_error` y `max_r2_error` se mantienen por incrementos sobre los bins
    tocados, y solo esos se respaldan. `max_r2_error` es un máximo y no admite
    mantenimiento por incrementos cuando el bin que lo ostentaba baja: se
    resuelve con el recuento de bins que alcanzan el máximo, y cuando cae a
    cero —el único caso que obliga a rehacerlo— se rehace sobre los bins
    **activos**, los que tienen alguna pareja asignada, cuyo conjunto es fijo
    porque las distancias no cambian y cuyo tamaño está acotado por
    `S*(S-1)/2`. Lo cierto, y es lo que importa, es que **el bucle de propuestas
    no depende de `B`**: cuando `B <= S*(S-1)/2` los bins activos pueden llegar
    a ser todos, y la siembra sí los recorre una vez por intento. La batería lo
    fija por dos vías: dos fixtures de la misma geometría con uno o `L-1` bins,
    cuyo recuento de visitas elementales del bucle debe ser **idéntico**; y la
    identidad exacta `max_rebuild_visits = max_rebuilds * bins_activos`, con los
    bins activos calculados por la propia prueba desde el layout y el perfil.
    La medida no es una estimación: `pair_updates` se cuenta dentro de
    `refresh_pair`, `loop_operations` cuenta las visitas elementales del bucle
    y `row_visits` mide las **tres** sumas canónicas del intento —siembra,
    recomputación de columnas y evaluación LD de cierre— en un único medidor.
    La batería exige las identidades exactas
    `row_visits = 2*H*(S + S*(S-1)/2) + H*S + H`, que se simplifica a
    `H*(S+1)²` y es **invariante frente al presupuesto de propuestas**, y
    `pair_updates = candidate_proposals*(2S-3)`, más la cota declarada
    `loop_operations <= Q*(24 + 8*S_bound)`, que escribe la **forma** del coste
    —una constante por propuesta más un término lineal en `S`— en lugar de un
    múltiplo del presupuesto grueso de v7. Los cuatro bucles que reajustan bins
    cargan su unidad **dentro** de la iteración, no en una línea escrita junto
    al bucle: escrita fuera sería otra identidad a mano y un recorrido que
    cambiase de longitud no se vería en la cifra, que es exactamente el defecto
    que la revisión de cierre encontró y que esta corrección cierra.
    Lo que `row_visits` **no** acredita: mide las sumas canónicas que pasan por
    `site_marginal`, `joint_marginal` y `refresh_realized_columns` con el
    medidor que se les pasa. Un recorrido de las H filas escrito a mano, o uno
    que declare su propio medidor, queda fuera de la cifra y la batería no lo
    vería. La afirmación anterior —«cualquier suma canónica lo incrementa venga
    de donde venga»— era falsa y se ha retirado: el medidor vive ahora en el
    estado de la búsqueda, de modo que una suma añadida dentro de un método
    tiene a mano el correcto y sí rompe la identidad, pero eso es una
    probabilidad, no una garantía.
    La memoria añadida son 32 bytes por pareja intrabloque más los 32 del
    veredicto canónico, dentro de los `192*S*(S-1)/2` que
    `f5-resource-accounting-v8` ya presupuesta.
  - Las cotas de dimensión las impone cada cual en su sitio, y la búsqueda no
    las duplica: `H` fuera de `2..=4096` lo rechaza el lector del perfil con
    `GPS-FMT-101`, y `site_bound_per_block > 4096` o
    `H*site_bound_per_block > max_haplotype_cells` los rechaza la admisión
    sellada con `GPS-RES-101`. Esa cota domina el producto exacto del layout,
    de modo que una comprobación aquí sería una rama inalcanzable.
  - M-2 de la revisión queda **cerrado sin cambiar el código**: los agregados
    de bin se mantienen con `+=`/`-=` planos sobre un total sembrado con
    Neumaier, y DTC-089 avala explícitamente esas actualizaciones incrementales
    `O(S)`. La revisión no consiguió provocar un cambio de terminación con cinco
    semillas y tres perfiles, y la búsqueda de una deriva que moviese un campo
    discreto quedó escrita en
    `the_incremental_score_never_diverges_from_the_canonical_one` con resultado
    negativo. La aritmética de Neumaier gobierna los valores publicados y el
    veredicto; la puntuación del bucle es determinista pero no canónica, y la
    guarda por pareja es quien vigila que la deriva no cruce una frontera
    discreta sin atribución.
  - Alcance de archivos y cambios ajenos al módulo, todos aditivos: dos
    constructores nuevos en `src/diagnostics.rs` para `GPS-DOMAIN-102` y
    `GPS-DOMAIN-104`, que el catálogo ya declaraba y ninguna ruta construía;
    `f5_uniform_index_with_source` pasa a `pub(crate)` en `src/randomness.rs`;
    los tres límites predeterminados de `src/f5_preflight.rs` pasan a `pub` con
    prefijo `F5_`; `src/f5_haplotypes.rs` suma `build_f5_block_attempt`,
    `ensure_f5_haplotype_cardinality`, `canonical_cmp`, el intercambio 2×2 y la
    recomputación canónica de columnas; `src/f5_ld.rs` suma
    `evaluate_f5_block_linkage_disequilibrium`, el contador de visitas de fila y
    la regla compartida de tolerancia de bin; `src/f5_profile.rs` suma un
    constructor `LdBin::for_test` bajo `cfg(test)`. Ningún cambio deshace
    trabajo previo.
  - Discriminación de la batería, medida contra doce variantes defectuosas:
    intento global todo-o-nada con 16 unidades (2 pruebas), recorrido canónico
    `O(H*S)` dentro de la propuesta usando el medidor del estado (3),
    `p_hat_j` obsoleta (4, con `GPS-INTERNAL-101` desde la guarda), `p_hat_i`
    obsoleta (4), signo invertido en el incremento de la conjunta (14), pareja
    `(i,j)` sin rehacer (6), séptimo campo comparado por digest (2 unitarias),
    ámbito RNG con índice constante en el ayudante (1 unitaria + 1 de
    integración) y en el llamador (1), recorrido de los `B` bins por propuesta
    (1), guarda de divergencia devuelta a `Ok(())` (1 unitaria) y contraste bit
    a bit de la columna publicada desactivado (1 unitaria).
    Dos variantes **no** discriminan y quedan declaradas: un recorrido canónico
    que declare su propio medidor local, que ninguna cifra publicada refleja; y
    retirar de `build_block` la comprobación `abs_error > marginal_max_abs_error`,
    que es una invariante interna cuya violación es inalcanzable por
    construcción —el prefijo elegido y su vecino difieren en el peso de una
    fila— y que por tanto ninguna prueba puede exigir. En su lugar, la batería
    contrasta el `max_abs_error` **publicado** contra el máximo
    `|p_hat_i - p_i|` de un oráculo entero, que es la parte que sí podría estar
    mal.
    Los cuatro campos instrumentales del registro (`candidate_proposals`,
    `row_visits`, `pair_updates`, `loop_operations`) son `pub` porque la
    batería vive en un crate aparte; no tienen destino en ningún artefacto y su
    documentación lo dice.
  - **Condición normativa pendiente, ajena a F5.2.** El término `Q*S_bound` de
    `f5-resource-accounting-v8` es **insatisfacible por construcción**: con
    `S_bound = 2` concede dos unidades por propuesta y DTC-089 exige cuatro
    sorteos por propuesta, de modo que ninguna implementación puede cumplirlo.
    Eso explica el factor residual que se mide incluso con el bucle óptimo. Se
    resuelve de una de dos maneras —definiendo «unidad de operación» a una
    granularidad más gruesa que la visita elemental, cosa que hoy **ningún
    documento hace**, o pasando el término a `Q*(c + S_bound)`— y pertenece al
    hito de la admisión, no a este. Debe resolverse **antes de que exista
    contabilidad de recursos en tiempo de ejecución**; hasta entonces no afecta
    a ninguna salida, porque la contabilidad es una cota de admisión y no un
    medidor de ejecución.
  - Residual declarado y no corregido: `settle_all_bins` recorre los `B` bins
    una vez por intento y la siembra reserva cinco vectores de tamaño `B`, lo
    que deja un término `A*bloques*N_pop*B` que v7 no presupuesta. Solo muerde
    con `Q` diminuto, y se cerraría sembrando desde los bins que las parejas
    nombran en lugar de desde los declarados. Queda para quien retome el
    módulo.
  - Limitación declarada de los oráculos: `D` y `r²` solo se contrastan con
    pesos **diádicos** —uniformes o `1/2, 1/4, 1/8, 1/8`—, donde todo el
    álgebra de DTC-089 cabe en `u128` sin perder un bit. Con pesos no diádicos
    el oráculo racional desbordaría, de modo que esos fixtures comprueban el
    desenlace y la guarda de divergencia, no los valores.

> **Revisión F5.2.** Revisión independiente del oráculo LD, consumo RNG,
> complejidad y mejor intento. No continuar con un resultado solo plausible.

## 6. F5.3 — fundadores intactos y F

- ✅ **F5.3.1.** Resolver para cada bloque×subpoblación el modo solicitado y
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
  - Estado implementado (11 de septiembre de 2026) y corregido tras la revisión
    adversarial. `src/f5_founders.rs` resuelve, para cada bloque×subpoblación
    de la búsqueda ganadora, el modo solicitado y efectivo, su origen y sus dos
    tolerancias, y asigna una fila haplotípica a cada copia: en probabilístico,
    un U53 por copia en su propio ámbito `founder_sample`; en controlado, el
    reparto por mayores restos con factibilidad previa `O(H)`, desempate único
    del grupo fronterizo con `controlled_quota_tie` y permutación Fisher–Yates
    de las `C` identidades de copia con `founder_controlled_permutation`. Es la
    primera ruta que emite `GPS-DOMAIN-103`, que hasta hoy no tenía ningún
    emisor.
  - Decisiones locales registradas:
    - **`requested_mode` es el modo que el documento pide para esa unidad, no
      el predeterminado de la subpoblación.** `source` es quien distingue de
      dónde vino —`subpopulation_default` o `block_override`— y hace
      innecesario reservar `requested_mode` para el predeterminado; con esta
      lectura una sustitución que repite el modo de su subpoblación sigue
      siendo distinguible. `effective_mode` coincide **siempre** con el
      solicitado: DTC-089 prohíbe corregir o relajar una cuota, de modo que una
      unidad inviable se rechaza y no se degrada a probabilístico. Es la misma
      relación que F4 ya exige entre `requested_mode` y `effective_mode` de
      `sequencing_errors`, donde el manifiesto y el efectivo se rechazan si
      difieren.
    - **Las filas se publican por `population_id`+`block_id`, no por
      `block_index`.** Es el orden literal de la lista `sampling` del efectivo
      y del manifiesto. No es el mismo orden: `block_id` es un identificador
      derivado con BLAKE3, así que con dos o más bloques ambos órdenes divergen
      —el fixture de `rows_are_ordered_by_population_and_block_id` lo
      demuestra—. Las unidades asignadas se publican en ese mismo orden para
      que no convivan dos órdenes canónicos en el mismo resultado; quien
      necesite el orden del layout ordena por `block_index`, que cada fila
      conserva.
    - **`GPS-DOMAIN-103` cubre la cuota, no la combinación `F>0`×controlado.**
      El catálogo dice literalmente que la cuota «no puede satisfacerse
      exactamente dentro de su tolerancia», y con `F>0` no hay ninguna cuota
      que no quepa: hay una configuración contradictoria. Por eso `F>0` con
      controlado es `GPS-CONF-101`, que es además el código con el que la
      admisión ya la rechaza hoy —`adapters.rs` comprueba el predeterminado de
      la subpoblación y cada sustitución controlada—, y `GPS-DOMAIN-103` queda
      para la cota incumplida y para el residual que ninguna cuota entera puede
      sumar. REQ-EVO-76 agrupa los dos rechazos en la misma frase; lo que
      comparten es el **momento**, no el código, y aquí los dos ocurren antes
      de abrir ningún ámbito.
    - **La cota se evalúa en copias, no en `n_h/C`.** `|n_h/C - q_h| <= t` se
      comprueba multiplicada por `C`: `|n_h - C*q_h| <= t*C`, con el
      predeterminado expresado exactamente como `0.5` copias y una tolerancia
      declarada convertida con `t*C`. Es la misma desigualdad sin el cociente
      `n_h/C`, que introduce un redondeo ajeno al contrato: `C*q_h` y su
      distancia a un entero contiguo son exactos en binary64 con pesos
      diádicos, mientras que `n_h/C` no lo es salvo con `C` potencia de dos.
      - La primera versión comparaba `|n_h/C - q_h|` y **rechazaba cuotas que
        el contrato admite**. El caso mínimo son tres fundadores haploides y
        dos haplotipos equiprobables: el error exacto y la cota valen ambos
        `1/6`, y en binary64 daban `0.16666666666666669` frente a
        `0.16666666666666666`. La revisión midió el alcance en racionales sobre
        `C = 1..64` —29 de 64 valores con pesos `1/2,1/4,1/8,1/8`, 22 con
        `1/2,1/2` y 13 con `3/4,1/4`—, de modo que no era un caso de esquina.
      - No es un cambio de contrato, sino su aplicación: DTC-089 tipa como
        binary64 los **restos** y prohíbe la tolerancia en su **igualdad**, no
        en la evaluación de una desigualdad que enuncia sobre `n_h/C` y `q_h`;
        `0,5/C` es exactamente el umbral que hace admisible el redondeo al
        entero más próximo, de modo que el resto `0,5` —el que genera el empate
        fronterizo que el contrato regula con tanto detalle— cae **sobre** la
        cota por diseño; y el proyecto ya evalúa condiciones así en aritmética
        exacta. La tolerancia publicada sigue siendo la del contrato, por
        haplotipo; la de copias es interna y el diagnóstico publica ambas.
    - **El residuo de redondeo del paseo de prefijos lo absorbe el último rango
      con peso positivo.** La mecánica es la misma con la que DTC-089 elige
      banda del espectro —recorrer `q` en orden y quedarse con el primer
      prefijo compensado estrictamente mayor que el U53—, pero `choose_band`
      devuelve la última fila cuando ninguna cruza, y aquí eso podría dar una
      copia a un `q_h = 0`. Un haplotipo de peso nulo no puede recibir copias
      sin romper la frecuencia realizada. En el modo controlado la protección
      es doble: un resto nulo nunca entra en el residual —los restos suman
      exactamente el residual— y, si entrara, `|1 - 0| = 1 > 0.5` copias lo
      rechazaría.
    - **`H = 1` no consume RNG.** Es la regla general del catálogo —«una
      elección con una sola alternativa no consume RNG»—, aplicada igual que en
      `choose_band`. No es alcanzable por la cadena completa: la admisión del
      perfil exige `H` en `2..=4096` y rechaza antes un perfil de un solo
      haplotipo con `GPS-FMT-101`. La frontera del muestreo no puede depender
      de esa comprobación ajena, y ese rango fija el techo real de 4096
      haplotipos.
    - **Un grupo solo es «frontera empatada» si el residual lo parte.** Con un
      `select_count` estrictamente entre cero y `g` se abre
      `controlled_quota_tie` y se publica el objeto de `quota_ties`; con
      `select_count == g` o `0` el grupo está incluido o excluido por completo,
      no consume RNG y no publica nada. Un grupo de tamaño 1 nunca cruza nada,
      de modo que tampoco aparece: llamar «empate» a un único rango obligado
      llenaría `quota_ties` de filas que no representan ninguna elección.
    - **Nada del desempate ocurre antes de la factibilidad, y se puede
      comprobar.** `plan_controlled_quota` no recibe ninguna fuente aleatoria:
      no puede consumir una palabra aunque alguien reordene su cuerpo. Y
      `select_tie_prefix` recibe el ámbito como un cierre que **solo invoca
      cuando hay frontera**, de modo que «un grupo incluido o excluido por
      completo no consume RNG» se prueba por la ausencia de la llamada, no por
      la lectura del código.
    - **El censo fundador se resuelve en este hito y se publica ordenado por
      `individual_id`.** El ámbito `founder_sample` lleva `individual_id` entre
      sus componentes, así que sin identidades no hay muestreo posible.
      `resolve_f5_founder_roster` usa el `founder_id` aportado tal cual y
      deriva el de cada `founder_key` con el contrato población×clave de
      REQ-IND-25, el mismo `IdentityRegistry` que ya usa F1/F2. El orden por ID
      es el que DTC-089 exige a la lista `founders` del efectivo y no cambia
      ningún sorteo, porque las identidades de copia se ordenan igualmente
      antes de permutarse. F5.3.3 debe **reutilizar** este censo para
      materializar `FCP_`/`FHA_` y la pertenencia.
    - **La frontera con F5.3.2 queda abierta y vacía.** El muestreo recorre las
      copias de cada individuo juntas y en orden de `copy_index`; F5.3.2 solo
      tiene que insertar, para `copy_index == 1` con `P = 2`, el `founder_f` y
      la reutilización de la fila de la copia 0. Este hito no publica
      `F_realized` ni calcula `P(ii)`/`P(ij)`, y una búsqueda negativa sobre el
      propio módulo fija que no abre el ámbito `founder_f`.
  - Guardas declaradas sin prueba propia, por ser inalcanzables desde la vía
    productiva: la incompatibilidad `F>0`×controlado, la sustitución duplicada,
    la sustitución que no corresponde a ninguna unidad y la comprobación de que
    las cuotas suman `C` las garantizan antes la admisión —o el preflight, que
    resuelve `B` y comprueba `block_index` en `0..B-1`— o la propia
    construcción del reparto. Se conservan porque `assign_f5_founders` es una
    frontera pública, pero **borrarlas no haría fallar ninguna prueba** y así
    queda dicho; la revisión lo verificó mutación a mutación. La propiedad que
    protege la última —la cuota suma `C` y las copias la realizan— sí la fija
    `controlled_quotas_sum_c_and_match_the_integer_oracle`.
  - Coste, por bloque×subpoblación: `O(C + H)` en los dos modos, medido por un
    contador que se incrementa **dentro** de cada bucle y se publica en
    `F5FounderUnit::operations`. La primera versión crecía como `C log C` y
    declaraba el término menor: ordenaba los `H` restos —seis comparaciones con
    los perfiles reales— y, sobre todo, **ordenaba las `C` asignaciones** al
    final de cada unidad, cosa que el plan no mencionaba. La revisión lo midió
    por la vía productiva: 25 visitas con `C=2` frente a un presupuesto de 8 y
    14 344 frente a 2 048 con `C=512`, con la razón creciendo como `log2(C)`.
    Las tres correcciones no tocan el RNG:
    - la permutación controlada permuta los **índices** `0..C` en vez de las
      identidades de copia. Fisher–Yates intercambia posiciones, de modo que la
      secuencia RNG y la permutación resultante son idénticas, pero el
      resultado se escribe ya en el orden canónico y no hay que reordenar `C`
      asignaciones;
    - el umbral de los mayores restos sale de una selección parcial
      (`select_nth_unstable_by`) en `O(H)` y el grupo empatado de una pasada
      lineal: ordenar los `H` restos no aportaba nada, porque dentro del grupo
      el orden lo decide el desempate y fuera de él solo importa el lado;
    - el paseo de prefijos del muestreo probabilístico se compone **una vez por
      unidad** —el vector de prefijos compensados no depende de la copia— y
      cada copia lo consulta por partición binaria, de modo que el coste baja
      de `C*H` visitas a `C*log2(H)`. La tabla produce exactamente el mismo
      rango que el recorrido secuencial, y una prueba lo contrasta sobre una
      rejilla de 512 umbrales y siete perfiles.
  - La identidad de los resultados se comprobó comparando todas las
    asignaciones, cuotas y empates de seis escenarios —controlado con 2, 6, 9
    y 16 copias y probabilístico con 8 y 9— antes y después de las tres
    correcciones: son idénticos. Con ellas el término `B*C_total` de v7 cubre
    el trabajo real y no queda ninguna condición de contabilidad declarada. El
    contador mide el trabajo que el módulo ejecuta en sus bucles; una
    ordenación reintroducida **sin contabilizar** no la vería, y por eso la
    cota se declara también aquí.
  - Limitación declarada de los oráculos: el contraste entero exige pesos
    **diádicos**, donde `C*q_h`, sus pisos y sus restos son exactos en
    binary64. Con la cota medida en copias, la estructura del reparto —pisos,
    estados, grupo fronterizo y prefijo— **y el veredicto** de factibilidad se
    contrastan para todo `C` de la barrida, sin excepción. La primera versión
    dejaba los `C` no potencia de dos dentro de un `if let Ok(...)`: 156 de los
    192 casos no comprobaban ningún rechazo —justo donde vivía el defecto de la
    cota— y el recuento de suficiencia solo medía los 36 restantes. Ahora la
    barrida exige comparar los 192 e incluir rechazos.
  - Pruebas y variante defectuosa que cada una atrapa. `tests/f5_founders.rs`
    recorre la cadena completa con 29 pruebas y `src/f5_founders.rs` añade
    catorce unitarias con oráculo entero, fuente espía y oráculo del paseo
    secuencial; la línea base pasa de 824 a 867. Cada una se verificó **contra
    su variante defectuosa**, aplicada una a una sobre el productor:
    - cota otra vez en unidades de `n_h/C` → caen las tres regresiones del
      bloqueante —el caso haploide mínimo, el reparto de seis copias y la
      barrida de cardinalidades—, el oráculo entero y la prueba de los dos
      estados (5);
    - restos ordenados de menor a mayor → el oráculo entero, la suma `C`, la
      cota inclusiva, el empate y el orden de filas (9);
    - comprobar solo el piso del grupo fronterizo → la cuota imposible, la de
      los dos estados y el oráculo (4);
    - cota exclusiva (`>=` en lugar de `>`) → 9 pruebas, la inclusiva incluida;
    - umbral de la selección parcial con `>=` → 14 pruebas, entre ellas el
      oráculo y la cuota imposible;
    - prefijo empatado sin reiniciar al incluir el grupo entero → el oráculo
      entero;
    - sufijo en vez de prefijo del grupo empatado → el empate reproducible;
    - candidatos empatados ordenados de mayor a menor → 5, el empate incluido;
    - ámbito del desempate abierto siempre → la prueba de que un grupo sin
      frontera no lo abre;
    - desempate abierto con el ámbito de la permutación → el empate;
    - no permutar los índices de copia → la permutación reproducible;
    - asientos ordenados por copia y luego por individuo → la permutación y el
      orden canónico;
    - partición binaria con `<` en vez de `<=` → el oráculo del paseo
      secuencial y la frontera exacta del U53;
    - residuo de redondeo al último rango en vez del último positivo → el peso
      nulo inalcanzable;
    - modo solicitado tomado siempre del predeterminado → los dos modos y el
      último índice legal;
    - tolerancia declarada ignorada → la cota más estricta;
    - filas ordenadas por `block_index` → el orden por `block_id`;
    - censo sin ordenar por ID → las claves fundadoras en controlado;
    - muestreo con `copy_index` fijado a 0 → el ámbito por copia;
    - ordenación de las `C` asignaciones contabilizada → la cota lineal del
      coste.
  - Hueco declarado y no cubierto: **ningún fixture alcanza un grupo fronterizo
    de más de dos rangos por la vía productiva**. Exige tres o más pesos con el
    mismo resto y `C` que parta el grupo, y la búsqueda haplotípica no converge
    con esos perfiles: se probaron 216 combinaciones de pesos uniformes de 4 y
    8 rangos, cuatro longitudes de bloque, seis objetivos `r²` y tres semillas,
    sin un solo caso terminado. Las unitarias sí lo cubren —un grupo de cuatro
    que abarca el residual entero, y grupos de hasta ocho en la barrida del
    oráculo—, de modo que el hueco es del fixture y no del algoritmo.
  - **Correcciones tras la revisión independiente de cierre F5.3 (12 de
    septiembre de 2026).** Se aplican sobre F5.3.1 sin reabrir ninguna de sus
    decisiones.
    - **La semántica de la cota controlada queda declarada en DTC-089.** La
      revisión observó que la desigualdad literal del contrato,
      `|n_h/C-q_h|<=t`, no es la que el código evalúa, y recomendó evaluarla de
      forma exacta en enteros sobre los pesos binary64. El usuario resolvió el
      12 de septiembre de 2026 **mantener la evaluación en copias y escribirla
      en el contrato**: la cuota ideal normativa es `fl(C*q_h)` y la cota es
      `|n_h - fl(C*q_h)| <= t*C`, con el predeterminado aportando `0,5` exacto.
      La raíz es que `0,5/C` es exactamente el umbral de redondeo al entero más
      próximo, y ese umbral solo es satisfacible en aritmética exacta: con `q`
      normalizado en binary64, un ideal semientero cae un ulp fuera.
      - Alternativa descartada 1, **evaluación exacta estricta sobre los pesos
        binary64**: rechaza configuraciones ordinarias. Verificado con
        `q={0,7, 0,3}` y cinco fundadores haploides —`fl(5*0,7)=3,5` y
        `fl(5*0,3)=1,5`, ambos semienteros exactos—: en decimales exactos la
        cuota `[3,2]` cumple la cota con igualdad, `|3/5-7/10|=|2/5-3/10|=1/10`,
        y evaluando sobre los pesos binary64 **ninguna** de las seis cuotas
        enteras la cumple. Es el defecto de F5.3.1 con otro disfraz: allí se
        rechazaba la mitad de los censos legales por un ulp; aquí se
        rechazarían las configuraciones cuyo ideal cae en un semientero.
      - Alternativa descartada 2, **ensanchar la tolerancia unos ulp**:
        obligaría a elegir y justificar un margen y cambiaría el valor
        publicado en `controlled_max_abs_error_effective` del efectivo y del
        manifiesto.
      - Consecuencia declarada en DTC-089 y no suavizada: un validador
        independiente que compruebe la desigualdad literal sobre los pesos
        binary64 puede encontrar cuotas publicadas que la exceden en un ulp.
        La desigualdad normativa es la de copias y quien valide debe usar esa.
      - Prueba:
        `a_half_integer_ideal_is_accepted_because_the_bound_is_evaluated_in_copies`
        fija el caso, comprueba que la cuota publicada `[3,2]` cumple la cota
        con **igualdad** y que las dos alternativas descartadas la rechazarían.
        Discriminación por los dos lados: evaluar `|n/C-q|` en binary64 mata 7
        pruebas de integración y 2 unitarias; evaluar un ulp por debajo de
        `0,5` —que es el efecto de la evaluación exacta en este dato— mata 13.
    - **El `0,5` exacto en copias ya tiene prueba.** Era la única defensa del
      cambio que motivó la revisión de F5.3.1 y no existía: sustituirlo por
      `default_tolerance * C` sobrevivía a la batería entera, porque los dos
      valores solo divergen para `C` como 49, 98 o 103 y **ninguna prueba usaba
      `C > 16`**. `the_default_copy_tolerance_is_exactly_one_half_and_not_a_rescaled_one`
      lo fija con `C=49`, `q={1/2,1/2}` y controlado: el exacto acepta la cuota
      `[25,24]` con error exactamente medio y el reescalado —`0,49999999999999994`—
      la rechaza con `GPS-DOMAIN-103`.
    - **El oráculo entero de cuotas ya admite pesos no diádicos.** Abortaba con
      «el oráculo exige pesos diádicos», de modo que la clase de entrada que
      revela la semántica de la cota no podía llegar a la comprobación. Ahora
      convierte el ideal normativo `fl(C*q_h)` a su racional exacto —todo
      binary64 es un racional diádico— y opera en `u128`.
      `the_integer_quota_oracle_accepts_non_dyadic_profiles` barre cinco
      perfiles de dos rangos con decimales no diádicos por `C ∈ {5, 7, 10}`:
      quince combinaciones, treinta unidades bloque×subpoblación, incluidos los
      ideales semienteros, los que caen un ulp por debajo de un entero
      —`fl(7*0,7)=4,8999999999999995`— y los que dejan un rango con cero copias.
      Al extenderlo se corrigió una propiedad que era un artefacto de los pesos
      diádicos: el residual se cuenta desde los **pisos**, no desde la suma de
      restos, porque `Σ fl(C*q_h)` no vale exactamente `C`.
      Los perfiles de tres rangos quedan fuera del barrido porque la geometría
      mínima del fixture no le da a la búsqueda LD ninguna matriz que cumpla a
      la vez bins y signos; el perfil diádico de `H=4` cubre el caso
      multirrango.
    - **El modo controlado ya está probado fuera de diploides.**
      `the_controlled_mode_works_outside_diploids` fija `P=1` y `P=3` con el
      perfil diádico: cuota exacta, `copy_index` recorriendo `0..P-1`,
      identidades materializadas y `F_realized` no aplicable.
    - **Precedencia de diagnóstico corregida.** Con `P=3`, `inbreeding_f=0,5` y
      controlado salía «la unidad controlada exige F efectivo igual a cero» en
      vez del de ploidía, porque `resolve_f5_sampling` corría antes que
      `resolve_founder_inbreeding` y releía el literal. Ahora el F se resuelve
      y se valida contra la ploidía **antes**, y llega a `resolve_f5_sampling`
      ya tipado en `PopulationFounders`, que además sustituye la búsqueda lineal
      del censo por bloque. La combinación sigue siendo inalcanzable desde la
      vía productiva —`adapters.rs` la rechaza en la admisión—, así que el
      cambio no tiene prueba de integración posible y así queda dicho; lo que sí
      se prueba es la lectura de la que depende el rechazo,
      `FounderInbreeding::effective_value`, cuyas dos direcciones enmascaradas
      por la admisión mata la unitaria correspondiente.
    - **Son dos, no una, las guardas duplicadas en `adapters.rs`.** La entrada
      anterior nombraba solo la incompatibilidad `F>0`×controlado; el tope
      `t <= 0.5/C` de la tolerancia declarada también lo comprueba antes
      `adapters.rs` para cada sustitución controlada, y neutralizar cualquiera
      de las dos sobrevive a la batería.
    - **`HaplotypePrefix::single` pasa a llamarse `single_rank`.** La regla del
      catálogo exime del sorteo a «una elección con una sola alternativa», y el
      campo fijaba «un solo haplotipo». Coinciden porque el perfil exige
      `weight > 0` en cada fila y `2 <= H <= 4096`, pero el nombre y el
      comentario decían la regla y comprobaban la implementación. No se cambia
      el comportamiento: derivarlo del número de pesos **positivos** alteraría
      el consumo de RNG en un caso que el perfil no admite.
    - Anotado para la revisión de catálogos, sin cambiar nada: `GPS-DOMAIN-103`
      no puede llevar `block_id` estructurado porque el catálogo de
      diagnósticos fija su contexto como `scenario_id,population_id,detail`, y
      está congelado por huella. F5.3.1 lo mete en `detail`, que es lo correcto.
    - Anotado, sin cambiarlo: `adapters.rs` rechaza una subpoblación con
      `controlled` y `F>0` aunque **todas** sus unidades estén sustituidas a
      `probabilistic`. DTC-089 impone la restricción por unidad, no por
      predeterminado; es un sobre-rechazo teórico y la guarda por unidad de
      `resolve_f5_sampling` es la correcta.
- ✅ **F5.3.2.** Aplicar F diploide y publicar F realizado; rechazar F en otras
  ploidías y registrar Ne no aplicable en F5.
  - Depende de: F5.3.1.
  - Procede de: REQ-EVO-46–52/76, DTPD-109 resuelta y DTC-089.
  - Cambia: emparejamiento y estadísticas fundadoras.
  - Se comprueba con: fórmulas bialélica/multialélica, F=0/intermedio/1,
    H_expected=0, P=1/P>2 y vector de Bernoulli.
  - Implementado en `src/f5_founders.rs` sobre la frontera que dejó F5.3.1:
    `sample_probabilistic_unit` recorre las copias juntas y en orden de
    `copy_index`, y la copia 1 de un diploide consulta `founder_f` **antes** de
    abrir su `founder_sample`. El resultado publica además `f_realized`, una
    fila por población×bloque en el orden de `sampling`.
    Decisiones locales:
    - **El orden de consumo entra como cierres, no como comentario.** DTC-089
      fija que `founder_f` se consume antes que el `founder_sample` de la copia
      1 y que, al reutilizar, ese segundo flujo no se abre; ninguna salida
      muestra ese orden, porque los dos ámbitos derivan subsemillas
      independientes. `draw_copy_rank` recibe los dos ámbitos como cierres
      —mismo patrón que `select_tie_prefix` en F5.3.1— y las unitarias
      observan cuál se abre, en qué orden y cuál no llega a invocarse, pasando
      un cierre que entra en pánico si se invoca.
    - **`F=0` y `F=1` no derivan siquiera la subsemilla.** El contrato regula
      el **consumo**, no la apertura: `f5_bernoulli` ya devuelve sin gastar
      palabra en los extremos, y el crate tiene precedente de derivar y no
      consumir —`F5SpectrumProfile::sample_site` y el propio `founder_sample`
      de un perfil con un solo haplotipo—. No derivarla es entonces una
      elección local, no una exigencia: ahorra una derivación BLAKE3 por
      individuo y es inobservable. `FounderInbreeding` distingue los tres
      regímenes —`Independent`, `Bernoulli(F)` y `AlwaysReuse`— y solo el del
      medio invoca el cierre.
    - **La guarda de ploidía recibe el F declarado, no un cero.** REQ-EVO-48
      rechaza *configurar* F con `P=1` o `P>2`, no solo aplicarlo, y un `0.0`
      efectivo no distingue «sin declarar» de «declarado 0». Por eso
      `effective_inbreeding_f` devuelve `Option<f64>` y
      `resolve_founder_inbreeding(ploidy, declarado)` rechaza cualquier
      declaración fuera de `P=2`, citando ámbito y ploidía. La admisión ya lo
      rechaza antes (`GPS-CONF-101`, campos `inbreeding_f,ploidy`), de modo que
      la guarda del módulo es una segunda frontera; a diferencia de las
      declaradas en F5.3.1, **esta sí tiene prueba propia**, porque la función
      es pura y las unitarias la llaman directamente.
    - **`H_expected = 0` es `not_applicable`, nunca cero; el resto de la rama
      es defensa inalcanzable.** Con una sola categoría realizada `q` vale 1
      exacto y la resta da cero **exacto**: ese caso, el del contrato, no
      depende del redondeo y sí tiene prueba. Las otras dos formas de no ser
      positivo —negativo o `NaN`— son inalcanzables: con más de una categoría
      el mínimo exacto de `C² - Σ n_h²` es `2C - 2`, de modo que
      `H_expected >= 2/C - 2/C²` domina cualquier error de redondeo. Restringir
      la rama a `== 0.0` no hace fallar ninguna prueba, y así queda dicho: se
      conserva para no dividir nunca sin denominador, no como lectura del
      contrato.
    - **`Σ q_realized²` se suma en orden de rango, sin compensación.** DTC-049
      gobierna las normalizaciones del perfil, no esta fórmula, y DTC-089 exige
      evaluar binary64 «en el orden escrito». Con una sola categoría realizada
      la suma vale exactamente 1 y `H_expected` exactamente 0, que es el caso
      del contrato.
    - **`target_f` vale `Some(F efectivo)` en `P=2` aunque nadie lo declare.**
      El contrato solo declara `null` para `P≠2`; en diploides el objetivo
      efectivo es 0,0 cuando no se declara, igual que el `inbreeding_f` del
      efectivo.
    - **La precedencia `[f5]` → subpoblación es ahora gobernadora del
      muestreo, y por eso tiene cobertura propia.** `effective_inbreeding_f`
      existía desde F5.3.1, donde solo alimentaba un rechazo; aquí decide los
      genotipos fundadores. El fixture de la batería aprende a escribir
      `[f5] inbreeding_f` y tres pruebas de cadena fijan que el F de población
      gobierna a quien no declara el suyo, que el de la subpoblación
      **sustituye** al de población —no al revés— y que cada subpoblación
      recibe el suyo y no el de otra. Sin ellas, las tres variantes
      correspondientes pasaban la batería entera, y la de la precedencia
      invertida produce homocigosis completa donde el contrato pide HWE.
    - **`F_realized` se publica también en modo controlado.** El contrato lo
      pide «por población×bloque diploide» sin restringirlo a un modo; con
      cuota exacta la frecuencia realizada coincide con el objetivo `q`, y la
      prueba lo comprueba con el mismo oráculo racional.
    - **Ne no aplicable se registra por ausencia, no con un tipo nuevo.** La
      primera versión publicaba `F5ReproductiveNe::NotApplicable`. Se ha
      **retirado**: el objeto `f5` del manifiesto enumera sus claves de forma
      cerrada y no incluye Ne, `manifest-general-v5.schema.json` tampoco y es
      `additionalProperties:false`, de modo que el tipo no tenía ni tiene
      consumidor, y su aserción era una tautología de tipos que ninguna edición
      del productor podía romper. Lo que queda es lo que sí se puede demostrar:
      F5 no calcula ningún Ne, y una búsqueda negativa de `k_bar`, `V_k`, `k_j`
      y `N_p` sobre el módulo impide que la fórmula de F7.3 aparezca.
      Observación declarada, no bloqueante: DTC-089 dice que F5 «registra Ne
      reproductivo como `not_applicable`», pero ningún artefacto F5 publicado
      tiene dónde escribirlo; si F5.4 necesita ese registro, hará falta abrir
      un campo, y eso es decisión de contrato, no de este hito.
    - **La búsqueda negativa de F5.3.1 se vuelve positiva.** Donde aquel hito
      exigía que el módulo no nombrara `founder_f`, este exige que lo abra
      **exactamente una vez**: un segundo consumidor no pasaría por las
      unitarias que fijan el orden.
    - **`realized_diploid_f` es total.** Un censo sin copias devuelve no
      aplicable —`C=0` no tiene frecuencia realizada— y un emparejamiento roto
      o un rango fuera de la matriz son invariante interna, no cálculo. Sin la
      primera comprobación, dos copias 0 de individuos distintos contarían como
      un heterocigoto inexistente.
  - Guarda declarada sin prueba propia, por ser inalcanzable desde la vía
    productiva: la que exige que cada fila de muestreo tenga F resuelto —el
    mismo censo alimenta los asientos y el mapa de F—. **Borrarla no haría
    fallar ninguna prueba** y así queda dicho. Las otras tres que lo parecían
    sí la tienen, y las tres se comprobaron mutación a mutación: la identidad
    del individuo antes de reutilizar —una unitaria llama a
    `sample_probabilistic_unit` con una pareja ajena y `F=1`—, la condición
    `copy_index == 1` —otra lo hace con tres asientos del mismo individuo y
    `P=3`, donde sin ella la reutilización se propagaría a la copia 2— y el
    censo vacío, cuya unitaria fija también el coste: sin la guarda, `0/0`
    daría `NaN` y el mismo `None`, y lo único observable es que se recorrerían
    los `H` rangos.
  - Coste: la copia 1 diploide con `0<F<1` añade **una** visita, contabilizada
    dentro de `draw_copy_rank` y devuelta al bucle que la suma, no junto a él.
    El recuento de `F_realized` es `O(C + H)` —una visita por individuo y otra
    por rango, contadas **dentro** de sus dos bucles— y se suma al
    `F5FounderUnit::operations` de su unidad. La cota lineal de F5.3.1 sigue
    vigente y su prueba pasa sin cambiarla; una cota inferior nueva, escrita
    solo con las partes deterministas del modo controlado, falla si el recuento
    se ejecuta fuera de la contabilidad.
  - Pruebas y variante defectuosa que cada una atrapa. `tests/f5_founders.rs`
    añade catorce pruebas de la cadena completa y `src/f5_founders.rs` doce
    unitarias; la línea base pasa de 867 a 893. El oráculo de `F_realized` es
    racional en `i128` —`1 - het*C²/(N*(C²-Σn_h²))`, una expresión distinta de
    la del productor— y el de la Bernoulli reconstruye `bits53 < floor(F*2^53)`
    desde el ámbito publicado. Cada prueba se verificó **contra su variante
    defectuosa**, aplicada una a una sobre el productor:
    - `F=0` reutiliza la fila de la copia 0 → el HWE copia a copia, el ámbito
      por copia de F5.3.1 y los extremos que no abren `founder_f` (3);
    - `F=1` no reutiliza → la homocigosis completa, `H_expected=0` y los
      extremos (3);
    - `F=0` o `F=1` resueltos como Bernoulli, que abriría el ámbito → los tres
      regímenes de consumo, y con `F=0` también el HWE copia a copia (2);
    - `founder_sample` abierto antes que `founder_f` → el orden de apertura y
      los extremos (2);
    - reutilizar abriendo igualmente el segundo flujo → los mismos dos;
    - Bernoulli invertida → la reutilización intermedia, el vector de Bernoulli
      y el orden (3);
    - Bernoulli con `<=` en vez de `<` —mutando `randomness.rs`— → el vector de
      Bernoulli y el orden (2);
    - sin la guarda de ploidía → la unitaria de la declaración fuera de `P=2`;
    - `H_expected=0` informado como cero → las dos pruebas de no aplicable, la
      unitaria y la de la cadena;
    - `H_observed` dividido por `C` en vez de por `N` → los dos oráculos de
      `F_realized`;
    - `F_realized` en valor absoluto, que borraría el exceso de heterocigotos →
      los dos oráculos;
    - frecuencia **objetivo** del perfil en lugar de la realizada —sustituyendo
      `n_h/C` por el `q_h` diádico— → los dos oráculos de `F_realized`, el
      reparto mixto y las dos pruebas de no aplicable (5);
    - sin la guarda de emparejamiento → la unitaria del censo roto;
    - `target_f` publicado también fuera de `P=2` → la prueba poliploide;
    - `F_realized` calculado en toda ploidía → la poliploide, el caso haploide
      mínimo y las cuotas controladas (3);
    - reutilizar sin comprobar el individuo → la unitaria de la pareja ajena;
    - recuento de `F_realized` sin contabilizar dentro de su bucle → el oráculo
      unitario y la cota inferior del coste (2);
    - recuento contabilizado pero no sumado al coste de la unidad → la cota
      inferior;
    - sin el respaldo al F de población (`.or(…)`) → la prueba del F de `[f5]`;
    - precedencia invertida, el F de población ganando al de la subpoblación →
      la prueba de la sustitución;
    - F tomado de otra subpoblación (`values().next()`) → la prueba de las dos
      subpoblaciones;
    - F efectivo a la mitad → la reutilización intermedia **y** la prueba
      distribucional;
    - sin la guarda del censo vacío → su unitaria;
    - sin la condición `copy_index == 1` → la unitaria de la copia 2 con `P=3`.
  - Las demás pruebas fijan el **mecanismo** —la copia 1 reutiliza si y solo si
    su Bernoulli lo dice—; una prueba distribucional añadida fija además las
    **fórmulas** de REQ-EVO-47/48, que es lo que el plan enuncia como criterio
    y lo que REQ-VAL-23 manda comprobar: con `F=1/2`, `q=(1/2,1/4,1/8,1/8)` y
    8 000 individuos diploides contrasta las diez categorías `i<=j` contra
    `q_i²+F·q_i(1-q_i)` y `2(1-F)·q_i·q_j`, con desviación máxima observada
    0,0045 frente a una cota de 0,02. Su sensibilidad está medida: mata un F
    efectivo a la mitad, pero **no** uno un 10 % menor, que sí mata la prueba
    de mecanismo. Las dos son necesarias y ninguna sustituye a la otra.
  - Bordes cubiertos: `F` exactamente 0 y 1, `F=1/2` con reutilización y
    muestreo independiente en la misma unidad, un homocigoto **por muestreo
    independiente** —que no es una reutilización—, `H_expected=0` por una sola
    categoría realizada, `C` impar en `P=1` y `P=3` —imposible en diploides—,
    censo vacío, censo mal emparejado, rango fuera de la matriz, cuota
    controlada **no exacta** —`C=6` reparte `[3,1,1,1]` y la frecuencia
    realizada deja de ser el objetivo, que es donde la variante «frecuencia
    objetivo» se distingue en controlado— y reparto **mixto**, con una
    sustitución que deja un bloque controlado y otro probabilístico en la misma
    subpoblación.
  - Deuda nombrada para F5.4: REQ-EVO-52 exige registrar también la
    **procedencia** del F (`inbreeding_f_source`, `population` o
    `subpopulation`). `effective_inbreeding_f` la conoce y la descarta, de modo
    que F5.4 tendrá que repetir la resolución de precedencia al escribir el
    efectivo y el manifiesto. No se devuelve aquí por el mismo motivo por el
    que se retiró el tipo de Ne —superficie pública sin consumidor—, y la
    lógica que tendría que replicar está ahora cubierta por las tres pruebas de
    precedencia.
- ✅ **F5.3.3.** Materializar antes de toda salida las identidades de fundador y
  su pertenencia única: `IndividualId`, `FCP_` por secuencia/copia y `FHA_` por
  bloque como asignación separada.
  - **Resuelto el 9 de septiembre de 2026 (F5-R04).** `FCP_` se deriva de
    `copy_index_u32`, no de `phase_id`. La fila anterior contradecía a la propia
    DTC-089, que en la proyección de verdad fija la identidad de una copia como
    `(reference_id, individual_id, sequence_id, copy_index)` y da esa columna al
    TSV de fundadores, mientras `PhaseId` es en el código un identificador
    aportado por la usuaria y sin constructor F5. Biológicamente el nombre
    también sobraba: en F5 no hay meiosis, de modo que las copias de un
    individuo no tienen origen parental que las distinga; cuando F7 introduzca
    la fase, será una propiedad de la copia y no su identidad. Con el cambio, el
    oráculo de F5.4.2 puede rederivar `FCP_` leyendo solo `truth.jsonl`.
    Se regeneraron el vector de `id-catalog-v4` y su huella de contrato; F5.0.3
    sigue abierto y ningún artefacto F5 se ha publicado fuera.
  - Depende de: F5.3.2.
  - Procede de: REQ-EVO-32/43/53/76, REQ-POB-01, ACE-89/91/117 y arquitectura §15.
  - Cambia: repositorio de verdad, pertenencia y enlaces `FCP_`/`FHA_`→`HAP_`.
  - Se comprueba con: reconstructor de cada copia y asignación, fase/ploidía,
    identidad a través de bloques/secuencias, pertenencia ausente/duplicada e
    IDs aportados/derivados; búsqueda negativa de `HaplotypeProcess`, cruces,
    mosaicos o mapas.
  - Estado implementado (12 de septiembre de 2026), pendiente de la revisión
    adversarial de la pausa «Revisión F5.3». `src/f5_founders.rs` añade
    `materialize_f5_founder_identities`, que consume la asignación de F5.3.1 y
    F5.3.2 —no decide ninguna, no abre ningún ámbito RNG— y publica un
    `F5FounderIdentities` con tres tablas y una consulta de pertenencia: las
    copias `FCP_` por subpoblación, individuo, secuencia y `copy_index`; las
    filas `HAP_` de cada bloque×subpoblación con sus alelos; y las asignaciones
    `FHA_` que enlazan copia y fila en cada bloque. `copy()` resuelve la
    pertenencia de un `FCP_` publicado, que es lo que `founder_membership`
    necesita para no inferir la subpoblación del identificador.
  - Decisiones locales registradas:
    - **Vive en `src/f5_founders.rs`, no en un módulo nuevo.** Es el
      «Inicializador fundador» de arquitectura §15, que posee a la vez `q`/`n_h`,
      F, `IndividualId`, `FCP_` y `FHA_`; separar la materialización habría
      dejado la guarda de pertenencia sin poder probarse, porque los censos y
      unidades que necesita tienen campos privados y solo un módulo hijo puede
      construirlos. La búsqueda negativa de F5.3.1/F5.3.2 sigue apuntando al
      mismo fichero y ahora cubre también lo que F5.3.3 prohíbe.
    - **Las copias existen por cada secuencia que el layout puebla, no por cada
      secuencia de la referencia. Decidido y cerrado; escrito en DTC-089.**
      DTC-089 identifica la copia «por secuencia … a través de todos los bloques
      de esa secuencia», `f5-founders.tsv` no tiene ninguna fila sin bloque y la
      proyección de verdad no tendría ningún genotipo que emitir para una
      secuencia vacía. La revisión adversarial dejó abierta la semántica: una
      secuencia **retenida** por la referencia y sin bloques no recibe copia, y
      eso podía ser deliberado o un defecto. El usuario lo cerró el 12 de
      septiembre de 2026 manteniendo esta conducta. Las otras dos opciones se
      descartaron por su coste: materializar también las copias vacías exigiría
      que la admisión sellara los **identificadores** de secuencia y no solo su
      recuento —reabrir el sellado—, y rechazar la divergencia tumbaría
      configuraciones legítimas, porque con muchas secuencias y pocos bloques lo
      normal es que la mayoría no reciba ninguno. Es, por tanto, la única opción
      que no reabre un contrato cerrado, y **su coste está declarado en
      DTC-089**: en una ejecución que retiene 25 secuencias y puebla 24, ningún
      fundador tiene copia de la secuencia 25 en `truth.jsonl`. La afirmación
      anterior de esta entrada —«la admisión solo conserva el recuento»— describía
      una limitación en vez de la decisión, y además sugería que el recuento no
      bastaba para **detectar** la divergencia, cuando sí basta: lo que no
      permite es nombrar la secuencia ausente.
      La prueba `a_retained_sequence_without_blocks_publishes_no_copy` fija la
      semántica: con `slack: 3`, dos secuencias, dos bloques y una semilla que
      los coloca ambos en `chr1`, `chr2` queda retenida y sin copias, y eso no
      es un error. Con el `slack: 1` del fixture por defecto el caso es
      estructuralmente imposible, que es por lo que ninguna prueba anterior lo
      tocaba.
    - **Cada tabla se publica en la clave literal que DTC-089 da a su fichero.**
      «Población, individuo, secuencia, copia, bloque» nombra columnas, y el
      orden es el de sus valores; donde DTC-089 quiere el orden FASTA lo dice
      expresamente, como en `f5-block-layout.tsv` («orden FASTA, inicio, fin,
      ID»). Por eso `sequence_id` y `block_id` ordenan lexicográficamente y no
      por ordinal de secuencia ni por `block_index`. Es la misma lectura que
      F5.3.1 registró para la lista `sampling`. El orden contractual de los
      artefactos sigue siendo de F5.4.1; aquí solo se entrega ya ordenado para
      que no convivan dos órdenes canónicos.
    - **`HAP_` se materializa aquí, para todos los rangos de cada
      bloque×subpoblación.** El plan asigna a F5.3.3 los «enlaces
      `FCP_`/`FHA_`→`HAP_`», y `FHA_` toma `haplotype_id` como componente: sin
      `HAP_` no hay `FHA_`. Se derivan los `H` rangos, no solo los asignados,
      porque `f5-haplotypes.tsv` publica la fila aunque ninguna copia la reciba
      y porque la identidad debe estar materializada antes de toda salida.
      `alleles` es la fila de la matriz ganadora tal cual: la polaridad de
      DTC-089 es una propiedad publicada del sitio en `f5-frequencies.tsv` y no
      reescribe esos bits.
    - **El `block_index` se conserva en las filas publicadas.** Es componente
      DTC-039 de `BLK_`, `HAP_` y `FHA_`, y ninguno de los siete TSV lo escribe;
      un lector externo lo recupera del ordinal de fila de `f5-block-layout.tsv`,
      que está en orden FASTA/inicio/fin. Conservarlo aquí evita que F5.4 tenga
      que reconstruirlo para derivar nada.
    - **La unicidad global del `IndividualId` se cierra en el censo, no en la
      admisión.** La primera redacción de esta entrada afirmaba que ninguna
      configuración admitida podía romper la pertenencia; era falso para la
      mitad del dominio y la revisión adversarial lo reprodujo: `adapters.rs`
      cruza los `founder_ids` **literales** entre sí, pero los `IndividualId`
      derivados de `founder_keys` viven en otro registro y nunca se contrastan
      contra los literales, y `validate_provided_id` acepta el texto `IND_…`.
      Un `founder_id` literal igual a un `IND_` derivado de otra subpoblación
      pasaba la admisión. El cruce se añade en `resolve_f5_founder_roster`, que
      es el primer punto donde ambos conjuntos existen: la admisión trabaja
      sobre el documento declarado, donde el identificador derivado todavía no
      se ha construido, y derivarlo allí obligaría a la máquina de admisión a
      abrir la derivación de identidades con la semilla raíz, cosa que hoy no
      hace. `conflicting_membership` se conserva como red del repositorio
      —una pertenencia duplicada **no** produciría colisión de `FCP_`, porque
      `population_id` es componente del material— y se prueba de forma unitaria
      sobre censos construidos a mano.
    - **El layout se vincula por identidad y por geometría.** La secuencia de
      cada copia viene del layout, de modo que un layout ajeno cambiaría
      silenciosamente los `FCP_`. Rederivar cada `block_id` —como F5.1.3 ya hace
      sobre sus marginales— vincula la ejecución pero **no** la geometría:
      `BLK_` y `SITE_` se derivan de índices, así que dos layouts de la misma
      ejecución con las mismas dimensiones y distintas coordenadas comparten
      identificadores. La primera redacción declaró ese hueco irresoluble «sin
      cambiar el material de una clase publicada»; era falso, y la revisión
      adversarial lo explotó con dos referencias de la misma ejecución, una de
      una secuencia y otra de dos, publicando `FCP_` sobre un `chr2` inexistente.
      Se cierra sin tocar ningún material publicado: cada bloque debe caber
      dentro de un **tramo elegible de esta admisión sellada**, con su mismo
      ordinal de secuencia y sus mismas coordenadas. La contención es exacta
      porque el layout coloca cada bloque dentro de esos tramos y solo los parte
      al consumirlos, y cierra las dos direcciones —secuencia que no existe y
      coordenadas fuera del tramo propio—. Residuo declarado: la admisión
      sellada conserva ordinales de secuencia, no identificadores, de modo que
      un layout cuyo `sequence_id` difiera pero cuyo ordinal y coordenadas
      coincidan seguiría pasando; publicar esos identificadores es de F5.4.
      Seguimiento nombrado: la contención no comprueba `Σ length_bp` contra el
      objetivo de DTC-089 porque `sealed.target_bases()` es `None` en toda la
      familia de fixtures —el objetivo se declara por `block_count`—, de modo
      que la guarda sería inalcanzable y no probable. Le corresponde a quien
      añada al fixture el objetivo `covered_bases`.
    - **La comparación `block_id` del layout contra el de la unidad se conserva
      y se prueba.** No forma parte del vínculo del layout —ambos se derivan del
      mismo índice y contexto—, pero sí vincula la **asignación**: una
      construida sobre la búsqueda de otra ejecución nombra otros `BLK_`. Emite
      `GPS-CONF-101` sobre `assignment.identity`, no un interno.
    - **La matriz ganadora debe tener los sitios de su bloque.** `alleles` es
      componente del `HAP_`, así que una matriz de otra geometría produciría
      identificadores internamente consistentes sobre una cadena que no
      describe el bloque. Se comprueba `site_count` contra los sitios colocados
      antes de derivar ningún `HAP_`.
    - **No se publica ningún contador de coste.** F5.3.1 y F5.3.2 publican
      `operations` porque el plan les pide complejidad medida; F5.3.3 no tiene
      ningún requisito de coste y `f5-resource-accounting-v8` está cerrado, de
      modo que añadir un contador sería superficie pública sin consumidor, que
      es justo el motivo por el que F5.3.2 retiró el tipo de Ne.
    - **La deuda de `inbreeding_f_source` sigue en F5.4.** F5.3.3 no escribe
      `effective-config.json#/f5/subpopulations` ni el manifiesto, que son los
      únicos consumidores de la procedencia del F; no le da consumidor y se deja
      donde F5.3.2 la dejó.
  - Cómo se comprobó, y contra qué falla cada prueba (32 de integración en
    `tests/f5_founders.rs`, 4 unitarias en `src/f5_founders.rs` y una de
    coherencia cruzada en `tests/normative_catalogs.rs`; batería completa
    930 pasan / 0 fallan / 12 ignoradas, frente a las 893 previas):
    - **Reconstructor.** `assert_every_identity_is_rederivable` rederive cada
      `FCP_`, `HAP_` y `FHA_` desde lo publicado —población, individuo,
      secuencia, copia, `block_index`, alelos e identificadores enlazados— más
      escenario, referencia, ejecución y semilla efectiva, con una
      implementación propia de DTC-037, BLAKE3 y Base32 que no llama a ninguna
      función del productor. Antes de juzgarlo,
      `the_independent_id_oracle_reproduces_the_frozen_catalog_vectors`
      reproduce material, digest e identificador de las **cinco** filas F5 de
      `id-catalog-v4`: un oráculo equivocado se detecta contra el vector
      congelado, no contra el productor.
    - Sensibilidad medida, aplicando cada defecto y ejecutando la batería:
      - `FCP_` derivado sin `copy_index` → 9 pruebas;
      - copias materializadas solo de la primera secuencia → 4, y la guarda de
        pertenencia ausente emite su `GPS-INTERNAL-101`;
      - `FHA_` derivado con `block_index` constante → 6;
      - `FHA_` enlazando siempre la fila 0 → 7;
      - `conflicting_membership` neutralizado → las 2 unitarias de pertenencia;
      - sin el orden canónico de las asignaciones → la prueba de orden;
      - alelos leídos siempre de la fila 0 → la prueba de filas haplotípicas;
      - sin rederivar el `block_id` del layout → la prueba del layout ajeno;
      - sin el vínculo con la admisión sellada → la suya;
      - productor invirtiendo `individual_id` y `sequence_id` en el material de
        `FCP_` → 6, incluido el reconstructor;
      - oráculo con versión de derivación 2 → 7, incluida la del catálogo;
      - pertenencia apuntando siempre a la primera copia → 7;
      - `copy()` devolviendo la primera copia ante un `FCP_` desconocido → la de
        pertenencia;
      - conjunto de secuencias sin orden (`HashSet`) → 2, por la aserción de
        determinismo;
      - «mosaico» escrito en el módulo → la búsqueda negativa;
      - sin la contención en tramos elegibles → las dos del layout ajeno;
      - sin el control de recuento de bloques → la suya;
      - sin la comparación `block_id` layout↔unidad → la de la asignación de
        otra ejecución;
      - sin la comprobación de sitios de la matriz → la suya;
      - sin el cruce literal↔derivado del censo → la del `founder_id` literal
        igual a un `IND_` derivado;
      - asignaciones ordenadas por `block_index` en vez de `block_id` → la de
        orden, que antes de corregir su fixture sobrevivía;
      - vector `FHA_` del catálogo con el `FCP_` antiguo → la coherencia cruzada
        de `normative_catalogs`, además de la huella y la del superconjunto;
      - copias publicadas también para las secuencias retenidas sin bloques → 4,
        empezando por la que fija la semántica de DTC-089.
    - Bordes cubiertos: un individuo con **dos** secuencias, la misma copia a
      través de **dos** bloques de una secuencia, `P=1` y `P=2`, `founder_ids`
      aportados y `founder_keys` derivados en la misma ejecución, la misma
      `founder_key` en dos subpoblaciones —que son dos fundadores, no uno—, el
      `individual_id` literal compartido entre dos subpoblaciones —rechazado ya
      en la admisión— y el literal igual a un derivado —rechazado en el censo—,
      un `FCP_` desconocido sin pertenencia, dos subpoblaciones declaradas en
      orden inverso al canónico, `P=1`, `P=2` y `P=3`, el modo controlado,
      `inbreeding_f>0` con reutilización de fila, diez secuencias donde el
      orden lexicográfico y el FASTA divergen, una secuencia retenida y sin
      bloques, y la materialización repetida, que da el mismo repositorio.
  - Correcciones tras la revisión adversarial Opus 5 alto (12 de septiembre de
    2026). La revisión confirmó el núcleo —los cinco materiales DTC-037 con sus
    componentes exactos y en orden, `FCP_` desde `copy_index_u32`, `FHA_` como
    asignación separada, `copy()` total e inyectiva, cero ámbitos RNG, `HAP_`
    para los `H` rangos, idempotencia— y señaló cinco graves, cuatro de ellos
    aplicados aquí:
    - **Layout ajeno de la misma ejecución (grave, aplicado).** Ver la decisión
      de vinculación del layout, reescrita arriba. Pruebas: layout de otra
      referencia con una secuencia que la propia no tiene; layout con
      coordenadas fuera de los tramos propios; y layout con otro recuento de
      bloques, con el fixture ajustado (`slack: 2`) para que sus tres bloques
      quepan dentro del tramo propio y el recuento quede aislado de la
      contención. Las tres mueren si se borra su guarda.
    - **Unicidad global del `IndividualId` (grave, aplicado).** Ver la decisión
      reescrita arriba. Prueba: se deriva primero el `IND_` de una `founder_key`
      de `P2` y después se declara ese texto como `founder_id` literal de `P1`;
      la admisión lo acepta y `resolve_f5_founder_roster` lo rechaza con
      `GPS-CONF-101`.
    - **La prueba de orden no ejercitaba `block_id` (grave, aplicado).** Su
      fixture tenía un bloque por secuencia, así que la tupla anterior ya era
      única y la quinta clave nunca desempataba: ordenar por `block_index` no
      mataba ninguna prueba. Ahora usa una secuencia con dos bloques, comprueba
      que cada copia recibe los dos y que los dos órdenes divergen de verdad.
    - **Vector `FHA_` del catálogo incoherente (grave preexistente, aplicado).**
      `id-catalog-v4` incrustaba `founder_copy_id=FCP_7TEWR4IDBCQZGCB2HTDB2YMQHY`,
      de la era `phase_id`, mientras su fila `f5_founder_copy` deriva
      `FCP_3MYOBVN22OO7JPHBBXMUFWFJQQ`: F5-R04 regeneró una fila y no la otra.
      La fila seguía siendo internamente consistente, de modo que ni el
      reconstructor ni `tests/f5_identity.rs` podían verlo. Se regeneró con el
      precedente explícito de F5-R04 —F5.0.3 sigue abierto y ningún artefacto F5
      se ha publicado fuera—: componente, material, digest e identificador. El
      nuevo identificador es `FHA_3H6PQOV7XJRT7HA3RGRDRDONFY`, que es
      exactamente el que `tests/f5_identity.rs` ya exigía al productor, así que
      catálogo y código quedan de acuerdo. Se actualizó la huella de contrato de
      `tests/normative_catalogs.rs` —único cambio de esa prueba junto a la nueva
      comprobación— y se añadió la **coherencia cruzada**: todo identificador
      que un vector nombre como componente debe derivarlo otra fila del mismo
      catálogo. Restaurar el `FCP_` antiguo la mata.
    - **Índice de la búsqueda (observación, aplicado).** El barrido lineal de
      `search.block(...)` por unidad recorría `B*P` filas `B*P` veces; se
      sustituye por un `BTreeMap` construido una vez, sin superficie pública
      nueva. El revisor midió y no observó el término cuadrático a las escalas
      del fixture; el cambio es barato y se aplica igualmente. Los `B*P*H`
      vectores de alelos siguen siendo `String`: es la forma que DTC-039 hace
      componente del `HAP_` y la que `f5-haplotypes.tsv` publica.
    - **H-2 (semántica con consecuencia biológica): resuelto por el usuario el
      12 de septiembre de 2026.** Se mantiene la conducta —copias solo para las
      secuencias que el layout puebla— y se documenta en DTC-089, no solo aquí.
      Ver la decisión reescrita arriba.
  - Huecos de cobertura cerrados por la revisión: el orden lexicográfico de
    `sequence_id` frente al orden FASTA —con diez secuencias `chr1…chr10` los
    dos divergen y ningún fixture anterior los distinguía, de modo que la
    decisión de orden no estaba probada—; el modo **controlado**, que las
    catorce pruebas iniciales no tocaban y donde la permutación de F5.3.1
    reparte los rangos; `inbreeding_f > 0`, donde la copia 1 reutiliza la fila de
    la copia 0 y debe conservar aun así su propio `FCP_` y su propio `FHA_`; y
    la ploidía 3, añadida al recorrido de ploidías.
  - **B2 — el repositorio deja de retener lo que crece con `B` (12 de
    septiembre de 2026).** La revisión de cierre midió que lo retenido excedía
    entre 7 y 8 veces el término `128*C_total` que `final_formula_peak`
    presupuesta, con `B*C_total*136 B` de asignaciones que **ninguna fórmula
    contabiliza**. El usuario eligió la forma perezosa y, después, contabilizar
    lo que falte.
    - **Lo que se conserva**: las copias `FCP_` con su índice de pertenencia
      —lo que F5.4 necesita en acceso aleatorio para `founder_membership`— y,
      por bloque×subpoblación, los `H` identificadores `HAP_` a 16 bytes. Se
      deja de conservar `alleles`, recomputable de la matriz ganadora que
      F5.2.3 ya retiene; `f5_haplotype_alleles` es ahora pública para que quien
      escriba el TSV lo haga desde la misma matriz de la que salió el material
      del identificador.
    - **Lo que se produce a demanda**: `for_each_assignment` deriva cada `FHA_`
      y lo entrega sin retenerlo, recorriendo copia × bloque, que es la clave de
      `f5-founders.tsv`. El otro orden publicado, el de `f5-haplotypes.tsv`, lo
      sirve `units()` sin recorrer nada: las unidades ya están por
      `population_id`+`block_id` y los rangos son el índice del vector. Ninguna
      de las dos claves obliga a una segunda pasada ni a ordenar `B*C_total`
      filas; desaparece el `sort_by` que había.
    - **«Materializar antes de toda salida» se sostiene así: la pasada previa es
      la derivación, y lo que prueba que ocurrió es el registro.** El registro
      de identidades vive dentro del repositorio y detecta toda colisión
      conforme a DTC-039. La garantía de que nada se publica a medias no la da
      la retención, sino DTC-089: el éxito publica los siete TSV y la
      proyección de verdad de forma atómica o falla entero, de modo que una
      colisión a mitad del recorrido aborta antes de publicar.
    - **El recorrido es único.** Derivar un `FHA_` dos veces es una colisión, de
      modo que `for_each_assignment` se consume una sola vez y lo dice con un
      diagnóstico propio en lugar de un `GPS-DOMAIN-002` confuso. Antes de
      recorrer comprueba que la asignación recibida es la que materializó el
      repositorio, posición a posición por subpoblación y bloque.
    - **Medición, en el mismo fixture, antes y frente al presupuesto**
      —dos subpoblaciones de dos individuos, `P=2`, `H=2`, `S=3`,
      `128*C_total = 1024 B`—:

      | caso | antes | después | factor antes | factor después |
      |---|---|---|---|---|
      | 1 secuencia, `B=2` | 3 672 B | 1 248 B | 3,6× | 1,2× |
      | 2 secuencias, `B=2` | 4 440 B | 2 016 B | 4,3× | 2,0× |
      | 1 secuencia, `B=8` | 12 384 B | 2 688 B | 12,1× | 2,6× |

      Proyectado a `B=1000`, una subpoblación, `H=4`, `C_total=2000` y una
      secuencia: las tablas pasan de ≈272 MB a ≈344 KB, de ~1 063× a ~1,34×.
    - **La razón no baja de 1, y el número es ese.** `128*C_total` concede 128
      bytes por copia del **individuo**, y el repositorio tiene una copia por
      individuo **y secuencia** —72 B— más su entrada de pertenencia —24 B—:
      con una sola secuencia eso ya consume el 75 % del presupuesto, y las
      unidades, que la fórmula no contempla en absoluto, lo desbordan. Es
      exactamente lo que el paso 2 tiene que contabilizar.
    - **El término dominante que queda no son las tablas, sino el registro de
      colisiones**, que conserva una entrada por identidad derivada:
      `C_total*secuencias + B*subpoblaciones*H + B*C_total`. A escala de
      producción son ~2·10⁶ entradas, del orden de 180 MB, y ninguna fórmula lo
      contabiliza. `F5FounderIdentities::registered_identities` lo publica para
      que el paso 2 pueda presupuestarlo y para que la prueba pueda medirlo.
    - **Prueba de la relación que faltaba**:
      `the_retained_repository_does_not_grow_with_the_block_count` compara dos
      cadenas idénticas salvo `B` —2 y 8— y fija que las copias no cambian, que
      lo retenido solo crece por los `HAP_` de cada unidad y **exactamente** por
      su tamaño, que lo retenido no excede el pico sellado, que las asignaciones
      emitidas sí crecen como `B*C_total` y que el registro crece con ellas. Un
      alambre de tropiezo sobre `size_of::<F5FounderIdentities>()` detecta una
      tabla retenida nueva, que `retained_bytes` no vería porque solo recorre la
      superficie pública.
    - **Discriminación verificada**: un identificador retenido de más por unidad
      mata 18 pruebas, incluida la de la relación; una tabla retenida nueva mata
      el alambre de tropiezo —256 frente a 232 bytes—; sin la guarda del
      recorrido único muere su prueba; sin la de la asignación propia, la suya.
    - **Ninguna prueba se perdió**: las que recorrían `assignments()` acumulan
      ahora dentro del cierre, recogidas una vez en un tipo `Materialized` de la
      batería, y la idempotencia compara copias, unidades, pertenencia y los dos
      recorridos en lugar de un `PartialEq` del repositorio.
    - **Paso 2: la contabilidad pasa a `f5-resource-accounting-v8` (12 de
      septiembre de 2026).** El paso 1 se detuvo porque
      `schemas/effective-config-v5.schema.json` fija
      `"algorithm": {"const": ...}` y `schemas/manifest-general-v5.schema.json`
      lo hereda por `$ref`, de modo que cambiar el nombre reabre dos esquemas.
      El usuario lo autorizó con el fundamento del §1 de este plan —«los pasos
      de esquema crean versiones nuevas sin reescribir las publicadas
      **externamente**; **ningún esquema F5 lo está todavía**»— y con el
      precedente de DTC-089 l. 1210, donde v7 «sustituye por completo» los
      parches v5 y v6. **No se crea una versión nueva de esquema**: efectivo y
      manifiesto siguen siendo v5 y lo único que cambia es el valor del `const`
      y sus lectores.
      - **Los cuatro términos que v7 no acotaba**, todos medidos sobre la
        cadena productiva:
        - `192*X` sustituye a `128*C_total`. La razón que no bajaba de 1 en el
          paso 1 era un desajuste de **recuento**, no un exceso del código: el
          repositorio conserva una copia por individuo **y secuencia**, que es
          `X=C_total*sequence_count`, ya definido por v7 para otros términos.
          Cada entrada mide 72 B de copia, 24 B de pertenencia y 48 B del texto
          `FCP_`: 144 B, que `128` tampoco cubría con una sola secuencia;
        - `64*B*C_total` presupuesta los asientos de
          `F5FounderUnit::assignments` de F5.3.1 —32 B exactos— que ninguna
          fórmula contenía;
        - `128*(X+B*N_pop*H+B*C_total)` presupuesta el registro de colisiones,
          el término dominante. DTC-039 exige que toda colisión se detecte, de
          modo que el registro conserva una entrada por identidad derivada;
        - `(160+S_bound)*H` sustituye a `64*H`: una fila haplotípica retenida
          cuesta 16 B de identificador, 48 B del texto `HAP_` y la parte
          proporcional de su unidad —108 B—, y `S_bound` cubre además su fila
          de la matriz ganadora, que se materializa por filas de `ceil(S/8)` y
          no en el empaquetado exacto de `ceil(H*S_bound/8)`.
      - **Medición en el mismo fixture del paso 1** —dos subpoblaciones de dos
        individuos, `P=2`, `H=2`, `S_bound=3`, `C_total=8`—, contando ahora
        todo lo que la fase conserva viva: tablas, textos de identidad,
        registro y asientos.

        | caso | vivo | v7 concedía | razón v7 | v8 concede | razón v8 |
        |---|---|---|---|---|---|
        | 1 secuencia, `B=2` | 4 544 B | 1 024 B | 4,44× | 7 960 B | **0,571** |
        | 2 secuencias, `B=2` | 6 192 B | 1 024 B | 6,05× | 10 520 B | **0,589** |
        | 1 secuencia, `B=8` | 13 232 B | 1 024 B | 12,92× | 24 160 B | **0,548** |

        La razón baja de 1 con un margen del 41 al 45 %, y ningún término
        individual queda corto: copias 1 136/1 536, filas haplotípicas
        848/1 304, asientos 512/1 024 y registro 2 048/4 096 en el primer caso.
        El pico sellado pasa de 80 900 a 87 324 B con `B=2` y de 110 864 a
        131 952 B con `B=8`.
      - **Ningún escenario legítimo se rechaza.** Cada término nuevo tiene un
        gemelo de salida que v7 ya imponía y que muerde antes: `X` está acotado
        por `truth_bound≥4096*X` contra el tope duro de salida —`X≤524 288`—
        mientras la memoria admitiría `X≤3,3·10⁶`; `B*C_total` aparece en
        `tsv_rows_bound` a 4 096 B por fila frente a los 192 B que v8 le suma;
        y `B*N_pop*H` igual. La memoria nunca pasa a ser la cota que decide.
      - **Punto fijo de los fixtures de anchura fija.** Ninguno de los términos
        nuevos depende de los bytes de los documentos, de modo que
        `final_formula_peak` sigue siendo el máximo y sigue siendo insensible a
        la anchura del literal que la prueba escribe después. Las dos pruebas
        que reescriben el límite —`profile_padding_changes_its_phase...` y
        `layout_streams_large_fasta_...`— vuelven a calcular desde los bytes
        finales y comprueban la igualdad: el punto fijo se alcanza en una
        iteración, como antes, y no hubo que reabrirlo.
      - **Sitios donde cambió el literal**: `src/f5_admission.rs`
        —`F5_RESOURCE_ACCOUNTING_ALGORITHM`, único lector—,
        `schemas/effective-config-v5.schema.json` —el `const`, que
        `manifest-general-v5` hereda por `$ref`—, las cadenas de
        `tests/f5_0_1_external_schemas.rs`, `tests/f5_resource_preflight.rs`,
        `tests/f5_founders.rs` y `tests/f5_search.rs`, y la prosa de
        `src/f5_preflight.rs`, `src/f5_search.rs`, `src/f5_layout.rs`,
        `src/application.rs`, DTC-089 y la arquitectura. `catalogs/` no publica
        ningún vector de esta contabilidad, así que no hay huella que regenerar.
        `evidence/f5.1.2/` conserva «v7» porque son actas fechadas de aquel
        hito, no contrato vigente.
      - **Discriminación verificada**, mutando el productor y midiendo con
        `--no-fail-fast`: devolver `192*X` a `128*C_total` mata la prueba nueva
        del diferencial de ploidía —2 560 frente a 2 816— y las tres del
        oráculo literal; suprimir `64*B*C_total` o el término del registro, lo
        mismo; devolver `(160+S_bound)*H` a `64*H` mata las tres del oráculo
        literal pero **no** la batería fundadora, porque el margen del registro
        tapa el déficit en el agregado. Queda dicho: ese coeficiente lo fija el
        oráculo de `tests/f5_resource_preflight.rs`, no la relación.
      - **La prueba de la relación, ampliada al registro.**
        `the_retained_repository_does_not_grow_with_the_block_count` compara
        ahora **término a término** —copias, filas haplotípicas, asientos y
        registro— contra los cuatro términos de DTC-089 escritos en la batería,
        no leídos del código: agregarlos antes de comparar dejaría que un
        término sobrado tapase a otro corto. La prueba nueva
        `the_sealed_budget_grows_exactly_with_the_founder_structures_it_retains`
        usa la ploidía, que es la única dimensión que mueve `C_total` sin tocar
        bloques, tramos, IDs, filas de perfil, `H`, `S_bound` ni bins, y exige
        que el pico sellado crezca **exactamente** por los tres términos
        fundadores.
      - **Observación fuera de alcance, no aplicada.** `F5HaplotypeSearch`
        retiene un `F5AttemptRecord` por unidad×intento —hasta `B*N_pop*A`, del
        orden de 264 B cada uno— que sigue vivo mientras se asignan fundadores
        y que ninguna fórmula de memoria contabiliza. No se toca: es una
        decisión de F5.2, cerrada tras su revisión adversarial, y presupuestarla
        con `A=100` por omisión multiplicaría el pico de los fixtures actuales
        muy por encima de sus límites declarados, que es la condición de parada
        de este encargo. Queda anotado para quien reabra F5.2.

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
    las tres puertas de `CONTRIBUTING.md` con `--locked`,
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
