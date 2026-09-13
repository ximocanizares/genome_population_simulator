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

- ✅ **F5.4.1.** Escribir conforme a los contratos ya cerrados en F5.0
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
  - **Implementado el 12 de septiembre de 2026.** Módulo nuevo
    `src/f5_artifacts.rs` y batería `tests/f5_artifacts.rs`; `src/lib.rs` lo
    registra y `src/truth_projection.rs` eleva `case_id` a `pub(crate)` para
    que la proyección v2 use el mismo identificador de caso que v1 en lugar de
    reimplementarlo. Ningún esquema cambia y ningún hito anterior se reabre.
    - **Frontera.** `F5SuccessRequest::new` enlaza configuración, admisión
      sellada, ejecución, layout, búsqueda, asignación y longitudes de
      secuencia, y rechaza antes de abrir ningún fichero una admisión que no
      sea la de esa configuración. `write_f5_success_artifacts` devuelve los
      siete `F5ArtifactLink` en el orden de la tabla —que es el lexicográfico
      de `artifact_id`— y el `F5TruthLink` aparte, cada uno con su objeto
      cerrado en `to_json`. F5.4.3 los copiará al manifiesto; este hito no
      escribe manifiesto, efectivo ni fallo.
    - **Decisión 0 — cada término de una clave nombra una columna de su
      propio fichero.** Es la regla completa de la que cuelga la decisión 1, y
      cierra sola las siete claves de DTC-089 sin dejar ningún término
      ambiguo:

      | Fichero | Clave de DTC-089 | Columnas que la nombran |
      |---|---|---|
      | `f5-block-layout.tsv` | FASTA, inicio, fin, ID | `sequence_id` (por ordinal FASTA), `start_1based`, `end_1based`, `block_id`; la clave está escrita entera porque «FASTA» no es columna |
      | `f5-sites.tsv` | bloque, índice | `block_id`, `site_index` |
      | `f5-haplotypes.tsv` | población, bloque, rango | `population_id`, `block_id`, `haplotype_rank` |
      | `f5-frequencies.tsv` | población, bloque, sitio | `population_id`, `block_id`, **`site_id`** |
      | `f5-ld-evaluation.tsv` | población, bloque, inicio de bin | `population_id`, `block_id`, `distance_start_bp` |
      | `f5-founders.tsv` | población, individuo, secuencia, copia, bloque | `population_id`, `individual_id`, `sequence_id`, `copy_index`, `block_id` |
      | `f5-attempts.tsv` | población, bloque, intento | `population_id`, `block_id`, `attempt_index` |

      La consecuencia que faltaba: **`f5-frequencies.tsv` no publica
      `site_index`**, de modo que su «sitio» solo puede ser `site_id`. El
      escritor lo emitía en el orden posicional de `block.columns()` —el del
      `site_index`— y ahora ordena por `site_id` dentro de cada bloque, con un
      vector de hasta `S_bound` índices por bloque. Los genotipos de
      `truth.jsonl`, cuya clave DTC-089 escribe «referencia, población,
      bloque, sitio e individuo», siguen la misma lectura: el único campo de
      sitio de un genotipo es el `site_id` de `f5_locus`. El oráculo de la
      batería ya no deriva el orden esperado de `f5-sites.tsv` —eso consagraba
      al productor—: exige la terna literal `(population_id, block_id,
      site_id)` estrictamente creciente.
    - **Decisión 1 — «bloque» como clave de orden es `block_id`.** Es la única
      lectura compatible con F5.3, que el encargo prohíbe reabrir:
      `resolve_f5_sampling` ordena por `population_id`+`block_id`, DTC-089 pide
      esa misma pareja para la lista `sampling` del efectivo, y F5.3.3 declara
      que su recorrido de asignaciones **es** el orden de `f5-founders.tsv` y
      que sus unidades **son** el de `f5-haplotypes.tsv`. Los seis ficheros
      cuya clave dice «bloque» usan por tanto `block_id`. El séptimo,
      `f5-block-layout.tsv`, no dice «bloque»: su clave está escrita entera
      —«FASTA, inicio, fin, ID»— y es la geometría, es decir `block_index`
      ascendente. La distinción no es cosmética: `block_index` es componente
      del material del `BLK_` y ninguna columna lo publica, de modo que un
      lector independiente solo puede recuperarlo del ordinal de fila del
      layout. La batería lo usa como oráculo de orden: rederiva los `BLK_`
      desde ese ordinal y falla si el layout se desordena.
    - **Decisión 2 — las longitudes de secuencia son una entrada explícita.**
      `truth.jsonl` exige `sequence_length` en cada fuente de segmento y
      ninguna pieza F5 previa la conserva: el layout retiene bloques y la
      admisión tramos elegibles, cuyo último extremo no es la longitud.
      `F5ReferenceSequences` la aporta, el escritor exige una entrada por
      secuencia poblada y comprueba que cada bloque cabe en ella. Es una
      entrada por secuencia retenida, exactamente lo que
      `256*retained_sequence_ids` de `final_formula_peak` ya concede. **F5.4.3
      debe cablearla desde el escaneo de admisión**; mientras tanto solo la
      batería la construye.
    - **Decisión 3 — forma de la copia fundadora.** Un fundador F5 es la
      referencia con los SNP de sus bloques aplicados, de modo que cada
      `haplotype_copy` publica **un** segmento por secuencia completa,
      `[0, sequence_length)`, cuyas modificaciones son los sitios donde su fila
      haplotípica lleva el alelo alternativo; `current_pieces` repite esa vista
      y `history` los mismos hitos. Un bit a 1 de la matriz es el alelo
      alternativo con independencia de la polaridad: DTC-089 usa la polaridad
      solo en la regla de signos de LD.
      **Corrección de la revisión F5.4.1:** `truth_bound` acota **bytes de
      salida**, no memoria de trabajo —solo entra en `output_bytes`,
      `src/f5_preflight.rs`— y `final_formula_peak` no reserva ni un byte por
      copia serializada. Materializar el registro con tres
      `collect::<Vec<Value>>()` vivos a la vez costaba 2,3× el pico sellado
      con `B=8` y 6,2× con `B=128`. El registro `haplotype_copy` se serializa
      ahora **por flujo**: se escriben sus claves en el orden canónico de
      `serde_json::Map` y las tres vistas recorren el búfer de `CopySnp` tres
      veces sin construir ningún `Value`. Los bytes publicados son idénticos,
      byte a byte, a los de la versión materializada.
    - **Decisión 4 — `reference_content_sha256` es el SHA-256 del FASTA
      sellado.** F5 no transforma la referencia, de modo que la huella de su
      contenido y la del artefacto declarado coinciden por construcción; no se
      inventa una huella semántica que nadie calcula.
    - **Decisión 5 — dos `source_digests` constantes.** La referencia
      —`kind="reference"`— y el perfil —`kind="matrix"`, con el SHA-256 de su
      metadato—. La procedencia fina viaja **por SNP** en
      `provenance.matrix`, que nombra el `HAP_` del que salió el alelo, la
      subpoblación y el `SITE_`. Así la cabecera no crece con los bloques y la
      traza FHA→HAP sigue siendo recuperable desde `truth.jsonl`.
    - **Decisión 6 — cada genotipo publica siempre los dos alelos del sitio**,
      también cuando la dosis de uno es cero, para que el recuento MAN-70 por
      dosis `0..ploidy` lea un índice de alelo estable en vez de deducirlo de
      los alelos que resultaron presentes.
    - **Decisión 7 — los SNP de una copia se ordenan por coordenada antes de
      serializarse.** Los bloques llegan en la clave publicada, que es un
      digest, y la secuencia de la copia es física. Se ordena el búfer de
      **una** copia, acotado por los sitios de su secuencia; ninguna tabla de
      filas se ordena.
    - **Streaming, y qué se retiene.** Las ocho salidas se escriben por flujo.
      `f5-founders.tsv` y `truth.jsonl` se abren a la vez y comparten la
      **única** pasada de `for_each_assignment`, porque comparten su clave
      —población, individuo, secuencia y fase—; derivar dos veces un `FHA_`
      es una colisión y la prueba lo fija. Lo único que el escritor retiene son
      índices de posición: `B` para el layout, `B` ordenados por `block_id`
      para los sitios, `B*N_pop` para las unidades y otros tantos para los
      tramos del historial de intentos. `final_formula_peak` ya concede
      `256*4*B` y `B*N_pop*(…+256*S_bound+…)`, muy por encima. Ninguna tabla de
      filas vive en memoria: `f5-attempts.tsv` publica el tramo contiguo de
      cada unidad sin copiarlo ni reordenarlo, y los genotipos se emiten unidad
      a unidad desde la matriz ganadora.
    - **La retención de `F5AttemptRecord` sigue sin contabilizar, y este hito
      no la agrava ni puede resolverla.** Escribir `f5-attempts.tsv`
      incrementalmente **no** libera nada: `search_f5_haplotypes` devuelve el
      historial completo y el escritor es el último paso de la cadena, de modo
      que los `B*N_pop*A` registros ya han vivido durante la asignación y la
      materialización antes de que exista un `ExecutionWorkspace`. Eliminar la
      retención exige que la búsqueda acepte un sumidero y escriba cada intento
      según lo ejecuta, lo que reabre la frontera pública de F5.2 y obliga a
      abrir el destino antes de saber si la ejecución tendrá éxito —y un fallo
      no puede publicar ninguno de los siete—. Queda anotado, como en F5.3.3,
      para quien reabra F5.2.
    - **Batería.** 33 pruebas en `tests/f5_artifacts.rs` —28 originales más
      las cinco de la revisión F5.4.1: fidelidad de `truth.jsonl`, ejecución
      ajena, declaración de secuencias más ancha que la sellada, validación
      del documento publicado contra su esquema externo y pico del escritor
      medido con asignador contador—. El lector es
      independiente: reimplementa DTC-037, BLAKE3 y Base32, se valida antes
      contra los cinco vectores congelados de `id-catalog-v4` —material
      canónico, digest e identificador—, y desde ahí rederiva los `BLK_`,
      `SITE_`, `HAP_`, `FCP_` y `FHA_` de los ficheros. `reconstruct_from_truth`
      no abre ningún TSV: recupera identidad, pertenencia, fase, alelos, dosis y
      recuentos por locus×subpoblación desde `truth.jsonl` v2 y exige que las
      dosis salgan de los portadores; `check_package` cruza después esa
      reconstrucción con los siete ficheros. El recuento por
      `site_id`×`population_id` se comprueba aparte: su histograma de dosis
      suma exactamente los individuos declarados en la cabecera v2, ninguna
      dosis excede la ploidía y el total de alelos alternativos coincide con la
      suma por individuo. Cobertura de ámbitos: dos
      subpoblaciones con censos distintos, dos secuencias pobladas y una
      estéril, ploidías 1, 2 y 3, modos probabilístico y controlado con
      sustitución por unidad, `F` de población 0 y 1, identificadores
      fundadores aportados y derivados, y bins LD `evaluated` y
      `not_applicable` en el mismo fichero.
    - **Discriminación verificada**, mutando el productor con `--no-fail-fast`.
      Primera ronda, 15 mutaciones: cabecera v1, copia sin
      `founder_membership`, `f5_locus` de otra subpoblación, sitios en orden de
      layout, genotipos por individuo antes que por sitio, fase constante,
      alelos invertidos, intentos en el orden de la búsqueda, recuento de
      copias desviado, dosis intercambiadas, SNP de la copia sin ordenar,
      enlace de verdad con versión 1, `NA` en un bin evaluado y haplotipos en
      el orden de la búsqueda matan entre 11 y 13 pruebas cada una. **La
      decimoquinta no mató ninguna**: fijar `polarity` a `positive` producía un
      paquete internamente coherente, porque la polaridad no es rederivable
      desde los propios ficheros. El agujero se cerró con cinco pruebas de
      fidelidad que contrastan contra el valor que F5.2 publicó las columnas
      cuya verdad no vive en los bytes —polaridad, las tres frecuencias, los
      diez campos LD, los nueve contadores de intento y `q`—, más la exigencia
      de que la polaridad publicada no sea constante. Segunda ronda, 10
      mutaciones de sustitución de columna —polaridad fijada, realizada por
      objetivo, `proposals` por `accepted`, objetivo LD por su tolerancia, `q`
      constante, `length_bp` desplazado, `abs_error` por el objetivo,
      `site_index` fijado, `block_assignment_id` por el `FCP_` y genotipo sin
      portadores—: fijar la
      polaridad y sustituir `proposals` por `accepted` matan exactamente la
      prueba de fidelidad que les corresponde, y las siete restantes matan 16
      pruebas cada una. **La `q` constante tampoco mataba ninguna al
      principio**, porque el fixture ancho usa pesos uniformes `0,5/0,5` y hacía
      degenerada la columna; con el perfil diádico `1/2, 1/4, 1/8, 1/8` la mata
      su prueba de fidelidad. Las dos lagunas son la misma lección de F5.3.2:
      el fixture tiene que distinguir la clave que decide.
    - **Corrupción aislada.** `every_isolated_corruption_is_detected_by_an_independent_reader`
      aplica 52 mutaciones a los bytes ya publicados —BOM, CRLF, LF final
      ausente, línea vacía, comentario, columna de más, cabecera alterada,
      celda entrecomillada, desórdenes de las seis claves, identificadores
      dañados de las cinco clases, `NA` fuera de su único caso, decimal no
      mínimo, alelo volteado, `selected_best` ausente, duplicado o no booleano,
      `truth` publicado como v1, versión y clase de proyección ajenas,
      proyección ausente, censo y recuentos de cabecera alterados, huella no
      hexadecimal, pertenencia y locus de otra subpoblación, dosis y portador
      alterados, fase distinta del índice de copia, copias desordenadas, alelo
      alternativo cambiado y sitio del locus cambiado— y exige que **cada una**
      sea rechazada. `a_missing_artifact_is_detected` cubre la ausencia del
      fichero y `every_link_declares_the_literal_path_size_and_digest_of_its_bytes`,
      la ruta literal, el tamaño `u64` y el SHA-256 de sus bytes.
    - **Revisión adversarial F5.4.1, correcciones aplicadas.** Nueve
      hallazgos; ocho aplicados y uno pendiente de una decisión del usuario.
      - **H1 — el escritor no cabía en su propio presupuesto (grave).**
        Corregido con la serialización por flujo de la decisión 3. Medido con
        asignador contador sobre `write_f5_success_artifacts`: antes
        233 868 / 847 360 / 3 502 772 B con `B=8/32/128` contra 100 792 /
        193 120 / 562 432 B sellados —2,32×, 4,39× y 6,23×—; después
        31 560 / 40 130 / 96 637 B, es decir 0,31, 0,21 y 0,17 del sellado, y
        el crecimiento total de `B=8` a `B=128` baja de 3 268 904 B a
        65 077 B, que es lo que cuestan los índices de posición.
        La prueba `the_writer_peak_fits_the_sealed_budget_and_does_not_grow_with_the_block_count`
        fija la relación con dos valores de `B`, como el gemelo de F5.3
        `the_retained_repository_does_not_grow_with_the_block_count`: exige
        `pico <= sealed.working_memory_bytes()` y acota el crecimiento por
        bloque a un kibibyte, un orden de magnitud por debajo de los ~25 KiB
        por bloque que costaba materializar los registros. Con el escritor
        anterior falla en los tres casos.
      - **Decisión 8 — `#![forbid(unsafe_code)]` vive en los raíces de
        `src/`, no en `Cargo.toml`.** Medir el pico exige instalar un
        `#[global_allocator]` contador, y `forbid` no puede levantarse en el
        objetivo donde se aplica: desde `[lints.rust]` de `Cargo.toml`
        alcanzaba también a `tests/` e impedía la única medición honesta de
        la relación que H1 destapó. El atributo pasa a `src/lib.rs` y a
        `src/main.rs` —dos raíces de crate distintos, de modo que la
        prohibición se repite en ambos— con una nota que explica por qué está
        ahí, para que nadie lo devuelva a `Cargo.toml` creyendo que se relajó
        por descuido. **No se relajó ningún contrato:** ningún documento
        normativo —`CONTRIBUTING.md`, `AGENTS.md` ni los de `design/`— fijaba
        la ubicación, y la garantía sobre lo publicado queda idéntica: ni un
        `unsafe` en `src/`. Comprobado introduciendo un bloque `unsafe` en
        `src/f5_artifacts.rs` y otro en `src/main.rs`: ambos son error de
        compilación —«usage of an `unsafe` block», con el `forbid` citado
        como origen— y se retiraron tras la comprobación. Lo único que queda
        fuera son los objetivos de prueba, donde el asignador contador
        declara su `unsafe impl GlobalAlloc` con la justificación de
        seguridad escrita: cada método delega en `System` y solo añade un
        contador de enteros.
      - **H2 — «sitio» leído de dos maneras (grave).** Aplicado: decisión 0,
        `write_frequencies` y `stream_genotypes` ordenan por `site_id`, y
        `check_package` compara contra la clave literal en vez de derivarla de
        `f5-sites.tsv`.
      - **H3 — `execution_id` ajeno (grave).** Aplicado. `F5FounderIdentities`
        publica `execution_id()` y el escritor rechaza con `GPS-CONF-101`
        —campo `execution.identity`— una solicitud cuya ejecución no sea la
        que derivó los identificadores. Se conserva el parámetro, que es lo
        que permite probar el rechazo.
      - **H4 — seis mutaciones supervivientes (grave).** Aplicado. El lector
        deja de comprobar solo la forma: exige que las fuentes de cada copia
        sean exactamente las de `source_digests`, que la procedencia de cada
        SNP nombre la matriz declarada y la subpoblación de su copia, que las
        tres vistas repitan la misma procedencia y que ningún `case_id` se
        repita. La prueba nueva
        `the_published_truth_carries_the_digests_the_trace_and_the_cases_it_received`
        contrasta contra la cadena las dos huellas de la cabecera, los dos
        `source_digests`, el `HAP_` que `f5-founders.tsv` asignó a cada copia
        en el bloque del sitio y el `case_id` de cada genotipo con un oráculo
        propio de DTC-061. Las seis mutaciones pasan a matar: 1 prueba las de
        `matrix_source_id` y `reference_content_sha256`, y 20 las otras
        cuatro.
      - **H5 — `F5ReferenceSequences` sin atadura (menor).** Aplicado: la
        declaración no puede nombrar más secuencias que
        `retained_sequence_ids` ni sumar más bases que
        `reference_bases_processed`, de modo que un `chr1 = 1e9` sobre una
        referencia de seis bases se rechaza.
      - **H6 — el `truth.jsonl` publicado no se validaba (menor).** Aplicado:
        cada línea publicada pasa por `truth-v2.schema.json` con ploidías 1, 2
        y 3, y la prueba fija que el `case_id` mide exactamente 64 caracteres,
        el máximo de `$defs/id`.
      - **H7, H8 y H9 (menores).** Aplicados: `mean_r2` y `abs_error` faltan
        si y solo si el bin es `not_applicable`; el layout afirma que el
        ordinal de fila es el `block_index` del que cuelga el `BLK_`; e
        `insert` de `F5ReferenceSequences` es `#[must_use]`.
      - **Huecos de cobertura cerrados.** El escritor comprueba además que
        `attempt_index` y `proposals` caben en los límites efectivos, que
        `copy_index` pertenece a `0..ploidy-1` y que los pesos haplotípicos
        suman 1 por población y bloque.
      - **Oráculos circulares.** Se corrigió el de H2 y el de la polaridad,
        que ahora ordena su `expected` por la clave literal en vez de heredar
        el orden de emisión del productor. El tercero
        —`the_published_founder_rows_are_the_single_assignment_traversal`,
        que compara contra `for_each_assignment` de una segunda cadena— **se
        conserva a propósito**: su objeto no es la verdad de las filas, que
        `check_package` rederiva con el oráculo independiente de DTC-037, sino
        la invariante de F5.3.3 de que el escritor publica esa pasada **sin
        reordenarla**; expresarla sin nombrar la pasada la vaciaría de
        contenido.

- ✅ **F5.4.2.** Implementar el oráculo independiente desde los bytes cerrados de
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
  - **Implementado el 12 de septiembre de 2026; revisión adversarial aplicada
    el mismo día.** Módulo nuevo `src/f5_oracle.rs` y batería
    `tests/f5_oracle.rs`; `src/lib.rs` lo registra. No se tocó
    `src/f5_artifacts.rs` ni ningún otro productor: el hito no destapó ningún
    defecto del escritor. La revisión sí destapó **un bloqueante y cinco
    graves en el oráculo**, todos aplicados: el alelo alternativo y el objetivo
    marginal se creían en vez de rederivarse, un ayudante de la batería mutaba
    la clave en vez del valor, una guarda del cierre de unidad marginal estaba
    muerta, una decimotercera variante defectuosa sobrevivía sin declararse y
    el recuento de filas JSON no contaba lo que decía contar.
    - **Frontera e independencia.** `audit_f5_package(payload_root,
      declaration)` recibe la ruta del paquete y una `F5PackageDeclaration`
      con las cifras que F5.4.3 copiará al manifiesto, y contrasta las
      segundas contra los primeros. **No recibe** `F5PhysicalLayout`,
      `F5HaplotypeSearch`, `F5FounderAssignment`, `F5FounderIdentities` ni
      `F5SuccessArtifacts`: ni un contador, ni una asignación, ni un
      acumulador del productor entra en el módulo. `recount_f5_man70(path)`
      es una función aparte que abre **solo** `truth.jsonl`; una prueba la
      ejecuta en un directorio donde ese fichero es lo único que existe.
    - **Decisión 0 — promoción, no duplicación ni sustitución.** La batería
      de F5.4.1 ya traía un lector propio y un `check_package` con 52
      mutaciones. Se **promueve** su disciplina a código productivo y se
      extiende al alcance entero del hito —RNG, D/`r²`, intentos, cuota,
      F, DTC-049, cotas y MAN-70—, con errores tipados `Diagnostic` y una
      causa enumerada en lugar de `String`. No se sustituye la de F5.4.1
      porque `tests/f5_artifacts.rs` está fuera del alcance de este encargo y
      porque es la evidencia propia del **escritor**: el oráculo juzga el
      paquete, aquella juzga quien lo escribe. La superposición está acotada
      —lector de bytes y rederivación de IDs— y F5.6.3 puede retirarla
      haciendo que aquella batería llame a `audit_f5_package`.
    - **Decisión 1 — la causa forma parte del contrato.** El plan exige que
      cada mutación falle «por la causa correspondiente», y un `String` no
      permite comprobarlo. `F5OracleRejection` lleva una `F5OracleCause` de
      21 variantes junto al diagnóstico público (`GPS-FMT-101`,
      `GPS-RES-101` o `GPS-INTERNAL-101`, con la causa también en su
      contexto). Las dos tablas de mutación contrastan la causa, no solo el
      fallo.
    - **Decisión 2 — el oráculo se valida antes de usarse.** `IdOracle` y
      `SubseedOracle` reimplementan el material canónico, BLAKE3 `derive_key`,
      Base32 y el ChaCha12 local, y sus constructores reproducen **primero**
      los cinco vectores F5 de `catalogs/id-catalog-v4.tsv` y todas las
      subsemillas y palabras de `catalogs/f5-rng-catalog-v1.tsv`. Si no los
      reproducen, el oráculo no juzga a nadie: devuelve `GPS-INTERNAL-101`.
    - **Decisión 3 — se rederivan las tres decisiones RNG cuyos componentes
      los bytes publican.** La primera versión decía que la polaridad era «la
      única», y **era falsa**: la revisión adversarial la destapó y con ella
      dos huecos, uno bloqueante.
      - `site_alt(scenario_id, block_id, site_id)` fija el **alelo
        alternativo** con `f5_uniform_index(rng, 3)` entre las tres bases
        distintas de REF, en orden `ACGT`. Sus tres componentes están
        publicados: `block_id` y `site_id` en `f5-sites.tsv` y `scenario_id`
        en la declaración.
      - `marginal_spectrum(scenario_id, population_id, block_id, site_id)`
        fija el **objetivo marginal**, con la mecánica cerrada de DTC-089
        —fila y valor, cada uno con su propio U53 y solo cuando tiene más de
        una alternativa—. Sus cuatro componentes están en la propia fila de
        `f5-frequencies.tsv`.
      - `polarity(scenario_id, population_id, block_id, attempt_index,
        site_index)` fija la **polaridad** con una Bernoulli(0,5). Su
        `attempt_index` es el único componente que no vive en una columna: lo
        da la fila `selected_best` de `f5-attempts.tsv`, que es el intento del
        que salió la matriz ganadora, y eso convierte `selected_best` en un
        dato necesario y no decorativo.

      Las demás decisiones del catálogo —longitudes, posiciones y huecos del
      layout, columnas iniciales, propuestas de intercambio, muestreo
      fundador, F, desempate de cuota y permutación controlada— no publican
      alguno de sus componentes y quedan declaradas fuera.
    - **Decisión 3 bis — las bandas del espectro entran en la declaración.**
      Rederivar el objetivo marginal exige `alt_af_min`, `alt_af_max` y el
      peso efectivo de cada banda. El efectivo publica las dos primeras cifras
      de cada fila —su auditoría DTC-049 lleva `alt_af_min` como clave y el
      peso como valor efectivo— y el oráculo **las cruza**; `alt_af_max` no
      aparece en ningún documento publicado y viene del perfil, atado por
      `profile_metadata_sha256`, que el oráculo contrasta contra
      `source_digests` de la verdad. Es la misma frontera que ya tenían
      `ld_bins` y `haplotype_weights`, y queda declarada como tal en el
      módulo. Se eligió rederivar en vez de declarar la frontera porque
      `input_alt_af` de los resúmenes MAN-70 es precisamente ese valor: creerlo
      habría dejado el vector esperado de REQ-MAN-70 sin ninguna atadura.
    - **Decisión 4 — `expected_*` de MAN-70 usa la frecuencia realizada.**
      `realized_alt_af` es `Σ_h q_h x_hi`, la probabilidad exacta de que una
      copia muestreada probabilísticamente lleve el alelo alternativo. El
      modelo es binomial en `hwe`, la mezcla de identidad por descendencia en
      `diploid_f` —`p²+Fp(1-p)`, `2p(1-p)(1-F)`, `(1-p)²+Fp(1-p)`, que es
      exactamente lo que `sample_probabilistic_unit` realiza— y la
      hipergeométrica inducida por permutar sin reemplazo las `C` copias en
      `controlled_hwe_conditioned_counts`. En este último el recuento alélico
      fijado `K` **no** es `realized_alt_af*C`: la cuota reparte `n_h` enteros,
      de modo que `K` es el número de copias de la unidad cuya fila lleva el
      alelo alternativo, que `f5-founders.tsv` y `f5-haplotypes.tsv` publican.
    - **Decisión 5 — `json_rows_bound` no acota las líneas de
      `truth.jsonl`.** Es una interpretación que el propio oráculo destapó al
      medirse: con censo 4 y ocho bloques la proyección publica 105 líneas y
      el bound vale 104. Sus términos son los del efectivo y el manifiesto
      —censo, muestreo, F realizado, resúmenes MAN-70 y filas de
      normalización— y su término de locus es `B*N_pop*S_bound`, **sin** el
      censo; las líneas de la proyección crecen con `individual_count` y viven
      en `output_bytes` a través de `truth_bound`. El oráculo acota por tanto
      las filas JSON que el manifiesto publicará, y los bytes de la proyección
      por `output_bytes`. No se cambió ninguna fórmula: se corrigió la lectura
      antes de que F5.4.3 la heredara.

      La revisión añadió que la cifra no contaba lo que decía contar: faltaban
      `founder_count` y los nueve objetos fijos. Con ellos son **64** filas en
      el fixture ancho —3 fundadores, 2 subpoblaciones, 8 de muestreo, 8 de F
      realizado, 24 resúmenes MAN-70, 2×5 de normalización y 9 fijas—, y la
      prueba fija ese entero exacto en vez de comprobar solo `<= bound`. Los
      `2*B*N_pop*H` rangos de `quota_ties` **no son observables desde los
      bytes** —ningún TSV los enumera—, de modo que la cifra es una **cota
      inferior** y el accesor se llama `observable_json_rows` para decirlo.
    - **Una guarda muerta y su riesgo, cerrados.** El cierre de unidad de
      `f5-frequencies.tsv` comparaba `realized_alt_af.len()` con el número de
      sitios del bloque, y ese vector se dimensiona al visitar la unidad: la
      comparación nunca podía ser cierta. Borrar la última fila por `site_id`
      de una unidad moría de rebote, por `max_marginal_error` o por el
      `pair_count` de un bin, con causa `Attempt` en lugar de `Marginal`. La
      sustituye un cierre real —los sitios publicados contra los del bloque, al
      cambiar de unidad y al terminar el fichero—. Eso cierra **por
      construcción** el riesgo que la revisión marcó sin reproducir: un sitio
      omitido cuya marginal realizada y objetivo fuesen cero no habría movido
      ni `max_marginal_error` ni ningún bin, y se habría publicado un resumen
      MAN-70 con `recount_matches = true` para un sitio polimórfico. Ahora la
      omisión se detecta por cardinalidad, sin depender de ningún valor, de
      modo que el fixture que lo exhibiría ya no puede construirse.
    - **Guardas inalcanzables, declaradas.** Cuatro defensas se conservan sin
      poder alcanzarse desde `audit_f5_package` y se anotan como tales en el
      código, igual que hace el resto del proyecto: la ploidía del recuento
      —que la cabecera ya contrastó—, el reparto total de la cuota —que el
      censo ya implica—, el recuento de bloques de una copia y la suma de
      `allele_dosages` contra la ploidía. Una quinta, `position <
      sequence_length`, no es alcanzable por la auditoría completa —el cruce
      con `f5-sites.tsv` mata antes cualquier posición movida— pero **sí por la
      frontera pública del recuento**, que no abre ningún TSV, y ahí se prueba.
    - **Fronteras declaradas, no supuestas.** Tres cifras se creen y el módulo
      lo dice: el **orden FASTA** del layout se aproxima por orden de primera
      aparición, porque ningún documento publica el diccionario ordenado de
      secuencias —la admisión sella su recuento, no su lista—; la **auditoría
      DTC-049**, el espectro, los bins LD y los pesos se rehacen desde los
      valores que la propia declaración trae, atados al perfil solo por
      `profile_metadata_sha256`, que el oráculo sí cruza contra
      `source_digests`; y `recompute_realized_f` es una **transcripción
      literal** de la fórmula de DTC-089 —el orden binary64 no compensado es
      parte del contrato—, con entrada independiente y tres casos racionales
      exactos, incluida la estimación negativa que REQ-EVO-52 conserva. Se
      declara además que el módulo comparte con el productor los literales de
      `F5_TSV_CONTRACTS` y `F5_TRUTH_*`: son el contrato que se juzga, no
      estado, pero son la única dependencia textual entre los dos módulos.
    - **Decisión 6 — dos bloques de una secuencia no pueden solaparse.**
      DTC-089 escribe la clave del layout como «FASTA, inicio, fin» y el orden
      estricto no prohíbe por sí solo que un bloque termine dentro del
      siguiente. El oráculo lo prohíbe porque de ello depende que `case_id`
      identifique un caso: dos loci de un mismo individuo en la misma
      coordenada darían el mismo `case_id`. Con el solape prohibido, la
      unicidad de `case_id` se comprueba con una entrada por **locus** y no
      por genotipo, que es lo que impedía que el recuento creciera con el
      paquete.
    - **Memoria.** El oráculo no materializa el paquete: cruza
      `f5-founders.tsv` con las copias de `truth.jsonl` en paso acompasado
      —ambas llegan en la misma clave— y retiene lo mismo que el productor:
      geometría `B` y `B*S`, matrices ganadoras en bits `B*N_pop*H*S`, un
      rango por asiento `B*C_total` y un vector de dosis por locus. El
      recuento MAN-70, que no puede compartir nada con lo anterior, retiene
      una entrada por copia. La única concesión sobre `working_memory_bytes`
      es el registro `haplotype_copy` que el lector tiene **en vuelo**: el
      escritor lo emite por flujo y la fórmula lo contabiliza en
      `output_bytes`, pero interpretarlo obliga a materializarlo una vez.
      **El exceso crece linealmente con `B`**, porque una copia cubre los
      bloques de su secuencia; lo que lo salva es que el registro más largo
      crece también linealmente, de modo que **la razón entre ambos
      converge**: medida, vale 4,56 con `B=8`, 7,67 con `B=32` y 8,48 con
      `B=128`, que es ya la asíntota del árbol `serde_json::Value` con
      `arbitrary_precision`. Decir «acotado por el registro más largo, no por
      el paquete» sin más sería engañoso, y la primera versión lo decía; lo
      acotado es la razón. No es el defecto de F5.4.1 —allí el pico crecía con
      el **número** de registros, no con su tamaño—.
      `the_oracle_peak_keeps_a_bounded_ratio_to_the_record_it_holds` mide los
      tres tamaños, exige la razón por debajo de diez y comprueba además que
      multiplicar por dieciséis los genotipos no duplica el pico.
    - **Condición para F5.4.3 y F5.5.** `f5-resource-accounting-v8` **no
      reserva ni un byte para el oráculo**, cuyo pico es del orden de 4,2 veces
      `working_memory_bytes` en el fixture medido. Ejecutarlo dentro del flujo
      de una ejecución exigirá antes un término propio en la fórmula o un
      lector por flujo del registro `haplotype_copy`. Mientras sea una
      herramienta de aceptación posterior, la cuestión no se plantea.
    - **Cobertura.** 18 pruebas de integración y 9 unitarias.
      `every_isolated_byte_mutation_dies_of_its_own_cause` recorre **74**
      corrupciones de bytes —contrato de bytes, layout, sitios, H×S/q,
      marginal, polaridad, D/`r²`, intento, asignación, identidades y
      `truth.jsonl`— y `every_isolated_declaration_mutation_dies_of_its_own_cause`
      **31** de la declaración —semilla, catálogos, DTC-049, modo, cuota,
      tolerancia, F y las siete cotas—. Cada una contrasta la causa exacta.
      Las 21 causas del enumerado quedan cubiertas. Los ámbitos: ploidías 1,
      2 y 3, ambos modos, `F=0` y `F=1`, IDs aportados y derivados, secuencia
      estéril, dos subpoblaciones, bins `evaluated` y `not_applicable` y un
      perfil diádico de cuatro pesos distintos.
    - **Mutación de bytes y enlace.** Una corrupción cambia siempre el
      SHA-256, de modo que juzgarla con el enlace intacto la mataría siempre
      por `Link` y taparía su causa. La batería rehace el enlace con los bytes
      mutados —el caso del defecto que reescribe artefacto y declaración a la
      vez— y comprueba la mitad del enlace por separado, con cuatro casos
      dedicados.
    - **Verificación de discriminación.** **Dieciocho** variantes defectuosas
      del oráculo se introdujeron y se retiraron: ALT creído, objetivo
      marginal creído, marginal realizada creída, media y error LD creídos,
      cierre de unidad marginal retirado, cruce copia/TSV de SNP retirado,
      `case_id` creído, portador repetido admitido, igualdad
      realizado/recontado retirada, tolerancia de bin LD retirada, cota de
      secuencia del recuento retirada, polaridad no rederivada, identidades
      creídas, cuota controlada no comprobada, F realizado creído, DTC-049 no
      auditado, contadores del mejor intento creídos y dosis recontada desde
      `allele_dosages`. **Diecisiete fallan** la batería; varias cambiando la
      causa de una mutación concreta, lo que es la prueba más fuerte de que la
      causa forma parte del contrato.
    - **La superviviente y su atribución corregida.** Retirar **solo** la
      comparación entre realizado y recontado no mata ninguna mutación. La
      primera versión de esta entrada dijo que la igualdad era «redundante por
      construcción» y que su papel era disparar si se debilitaba cualquier
      eslabón; **la revisión demostró que ambas afirmaciones estaban mal
      atribuidas**. Lo cierto, y reproducido en las dos direcciones: la
      igualdad **sí tiene dientes**, y los tiene frente a la guarda de
      portador repetido. Con esa guarda en pie, «portador repetido con sus
      dosis ajustadas» muere por la verdad; retirando solo la guarda, muere
      **por MAN-70**, es decir por la igualdad; retirando las dos, **pasa
      inadvertida**. Lo que sí es cierto es que ninguna mutación la aísla
      mientras el resto de la cadena está en pie.
    - **Cuatro mutaciones que cambiaron de causa al escribirlas.** «frecuencia
      realizada movida» y «media r² movida» morían por el error asociado, no
      por la recomputación; se añadieron sus gemelas coherentes —que mueven
      también `abs_error` y se quedan dentro de la tolerancia— y solo rehacer
      `Σ_h q_h x_hi` y `D`/`r²` las mata. «objetivo marginal movido» no movía
      el error, y su gemela lo mueve dejando intacto el máximo de su unidad.
      Y la sustitución del ALT solo se vuelve concluyente cuando cambia sus
      **cinco** vistas y el `case_id`: con una sola, la coherencia interna la
      mata antes de llegar al sorteo.
    - **Un falso negativo en la propia batería.** `replace_after` buscaba el
      literal a sustituir **desde el inicio del ancla**, de modo que «case_id
      ajeno» y «huella de la referencia alterada» acertaban una letra de la
      clave —`"cbse_id"`, `"reference_content_shz256"`— y morían por
      `exact_object` sin tocar el valor. Es exactamente el patrón
      `false && A || B` contra el que el encargo advertía. Corregido con
      `+ anchor.len()`; sin la corrección, desactivar la guarda de `case_id`
      dejaba la batería verde y el oráculo DTC-061 sin probar.
- ✅ **F5.4.3.** Instanciar la configuración efectiva y el manifiesto conforme
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
  - **Módulo y frontera.** `src/f5_publication.rs`, con
    `publish_f5_success`/`publish_f5_failure` y la puerta física pública
    `verify_published_f5_package`. `src/publication.rs` gana tres piezas
    genéricas: `payload_root`, `write_root_document` y `write_f5_manifest`, más
    `discard_all_artifacts`. v5 **no** es una ampliación de v4 —su raíz no
    tiene `artifacts`, ni `references`, ni `effective_config_artifact_id`—, de
    modo que la familia del manifiesto elige la puerta física en el cierre y
    otra vez antes del renombrado; la vía v4 queda intacta.
  - **`effective-config.json` es documento raíz, no artefacto.** La lista
    `artifacts` del manifiesto v5 es cerrada a los siete TSV —`prefixItems`
    con `items:false`— y no tiene campo donde enlazar el efectivo. Catalogarlo
    lo dejaría sin declarar, que es justo lo que v4 evitaba con
    `effective_config_artifact_id`. Viaja junto a `manifest.json` con la misma
    disciplina de temporal y renombrado. Consecuencia deliberada: un fallo
    controlado **no** lo publica, y la puerta física comprueba su ausencia en
    lugar de confiar en la purga del catálogo, porque la purga solo alcanza a
    los artefactos.
  - **Rutas publicables.** `$defs/path` prohíbe la ruta absoluta, de modo que
    la ruta física no es publicable: vive, completa y con su SHA-256 observado,
    en `resource_preflight.input_identities[].normalized_path`, cuyo esquema no
    tiene patrón. Las rutas del efectivo son **relativas a la raíz de
    autoría** —el directorio del general—, y `source_path` es el nombre del
    general. Una ruta de autoría con un carácter fuera de `[A-Za-z0-9._-]` no
    es representable en el contrato publicado, y **la admisión la rechaza
    antes** con `GPS-CONF-101`: los cinco enlaces declarados por
    `require_relative_path` y el nombre del propio general por
    `require_publishable_general_name`. La rama equivalente de
    `publishable_path` es profundidad defensiva y hoy es inalcanzable. Es una
    consecuencia del esquema, no una restricción añadida aquí.
  - **Los resúmenes MAN-70 son entrada, no cálculo.** `recounted_*` solo puede
    nacer de releer `outputs/truth.jsonl`, y REQ-MAN-70 prohíbe reutilizar
    acumuladores del productor. `publish_f5_success` recibe el recuento como
    una función del llamador, a la que entrega la raíz de carga útil ya escrita
    y la declaración que este hito construye. El productor **no puede**
    fabricar un `F5InitializationSummary`: sus campos son privados de F5.4.2,
    de modo que el tipo impide la autocertificación. Si la vía productiva
    comparase realizado y recontado producidos por el mismo lector, la
    igualdad sería vacua; aquí realizado y recontado vienen ambos del oráculo
    y la igualdad que el manifiesto publica es la que ese lector ya exigió
    sobre los bytes.
  - **El oráculo queda fuera del flujo, y se dice.** `f5-resource-accounting-v8`
    no reserva ni un byte para él y su pico medido en F5.4.2 es ≈4,2 veces
    `working_memory_bytes`. F5.4.3 **no lo ejecuta**: la frontera lo recibe.
    Darle un término propio exigiría v9, y el literal v8 está en los dos
    esquemas que este hito no puede tocar; declararlo como aportación ajena
    exigiría sellar la admisión con `F5OtherResources` distinto de
    `default()`, que las cinco fronteras F5 anteriores fijan literalmente. La
    condición sigue abierta y se vuelve bloqueante en el primer `run`
    productivo, es decir en F5.6.
  - **Deuda 1 cerrada: `inbreeding_f_source`.** `effective_inbreeding_f` ya no
    descarta la procedencia: `resolve_f5_founder_roster` la resuelve y el
    censo la publica junto al `F` ya resuelto contra la ploidía. La regla,
    declarada: la fuente es `subpopulation` **si y solo si** la subpoblación
    declara su propio `inbreeding_f`; en cualquier otro caso —el escenario lo
    declara para todas, o nadie lo declara y el efectivo vale `0` en
    diploides— es `population`, porque el valor efectivo procede entonces del
    ámbito poblacional. El censo absorbe también el modo predeterminado de la
    subpoblación, de modo que el literal `sampling` se interpreta **una** vez
    y `resolve_f5_sampling` lo consume en lugar de repetir el análisis.
  - **Deuda 2 cerrada: las longitudes de secuencia.** El escaneo de admisión es
    el único punto de la ejecución donde la longitud existe —nace de contar
    bases y muere al cambiar de cabecera; el tramo elegible más alto no la da,
    porque una cola de `N` la deja corta—. `EligibleScan::record_retained_id`
    la retiene junto al identificador, bajo el mismo cargo de 256 bytes por
    secuencia que v8 ya concede a `retained_sequence_ids`: **no añade ningún
    término a la fórmula**. `F5ReferenceSequences::from_sealed` la cablea, y
    las guardas de F5.4.1 —recuento, cobertura de las secuencias pobladas y
    suma de bases— pasan a ser redundantes por construcción y valen como
    comprobación cruzada. **Medido por la revisión adversarial**: 96 B por
    secuencia frente a los 256 del presupuesto, con 400 secuencias de
    identificador máximo. Con un matiz que no estaba dicho y ahora sí:
    `from_sealed` construye un **segundo** catálogo que coexiste con el vector
    sellado y cuesta 127 B por secuencia, de modo que el total real es
    **223 B por secuencia, un 87 % del presupuesto**. Cabe, y queda declarado.
  - **Deuda 3 cerrada corrigiendo el contrato: Ne reproductivo (2026-09-12).**
    DTC-089 decía que «F5 no reproduce y registra Ne reproductivo como
    `not_applicable`». Esa frase estaba **sola contra tres fuentes**, una de
    ellas la propia decisión que la contenía: la enumeración del objeto `f5`
    del manifiesto —escrita en ese mismo DTC-089— no contiene ningún campo
    donde escribirlo; `manifest-general-v5.schema.json` es
    `additionalProperties:false` y tampoco lo tiene; y REQ-MAN-92, que enumera
    lo que F5 debe registrar, no lo pide. **Se corrige DTC-089**, no el
    esquema: la frase queda en que F5 no reproduce y que el Ne reproductivo no
    le aplica —no lo calcula ni lo registra—, con la razón escrita en el
    propio contrato para que nadie vuelva a tropezar con ella, y remitiendo la
    fórmula por contribuciones a F7.3, que es la primera fase con parentales.
    **La enmienda es decisión expresa del usuario**, tomada tras plantearle las
    tres opciones —corregir la frase, añadir el campo al manifiesto o aplazarlo
    a la Revisión F5.4— y no una resolución del implementador.
    Dos alternativas descartadas: **añadir el campo** cumpliría la letra de
    una frase contradiciendo las otras tres, reabriría un esquema que este
    hito tiene prohibido tocar y añadiría superficie que ningún consumidor
    lee, que es exactamente por lo que `F5ReproductiveNe` se retiró en F5.3.2;
    **dejarla abierta** arrastraría la contradicción a F5.5 y F5.6. La
    búsqueda negativa `f5_never_computes_a_reproductive_ne` de F5.3.2 sigue
    pasando y ni `src/`, ni `schemas/`, ni `catalogs/` mencionan
    `reproductive_ne`.
  - **Contadores de un fallo controlado.** `layout_attempts_used` de un
    agotamiento de layout es `max_layout_attempts`, porque `F5BestLayout` solo
    existe por la vía de agotamiento y esa vía recorre el bucle entero; en
    cualquier otra causa de layout es cero. En un agotamiento de búsqueda, los
    cuatro contadores describen **la unidad que el diagnóstico nombra**: los
    agregados de otra unidad no explicarían por qué esta se agotó.
  - **Cobertura.** 30 pruebas de integración y 6 unitarias. Ámbitos: ploidías
    1, 2 y 3, ambos modos con sus dos orígenes, `F` en `0` y `1`, las dos
    procedencias de `inbreeding_f` y su `null` fuera de `P=2`, IDs aportados y
    derivados, secuencia estéril, dos subpoblaciones, perfil diádico de cuatro
    pesos y bins `evaluated`/`not_applicable`. Los tres documentos que el hito
    instancia validan contra sus esquemas con el validador Python externo.
  - **Verificación de discriminación.** Catorce variantes defectuosas se
    introdujeron y se retiraron; **trece** tumban la batería. Fuente de F
    siempre `population`, longitud de secuencia desplazada en uno, copia
    doctorada de `resource_preflight` en el manifiesto, igualdad
    efectivo/manifiesto retirada, fallo que no purga los artefactos,
    `verify_link` que cree el tamaño y el hash declarados, `source_path`
    absoluto, normalizaciones resumidas a su primer valor, octavo enlace en
    `artifacts`, revalidación previa al renombrado retirada, comparación
    catálogo/enlaces retirada, `inbreeding_f` nunca `null` y comprobación de
    artefactos biológicos supervivientes retirada. La batería de divergencias
    recorre **veintiuna** mutaciones de los dos documentos publicados y exige
    que el rechazo cite el término tocado, no solo que falle.
  - **La superviviente.** Retirar la guarda `!summary.recount_matches()` del
    constructor del manifiesto no mata ninguna mutación, y es exacto: el
    oráculo rechaza antes de devolver un resumen desigual, de modo que hoy la
    guarda es inalcanzable por construcción del tipo. Lo que sí tiene dientes
    es la comprobación equivalente de la puerta física sobre el documento
    publicado —`recount_source` y `recount_matches`—, que dos mutaciones matan.
    Se conserva la guarda porque el llamador del recuento es una función
    ajena y el contrato no puede descansar en que ninguna futura la relaje.
  - **Revisión adversarial aplicada (2026-09-12).** Un bloqueante, dos graves,
    un menor y varias notas. El revisor retiró ocho sospechas tras
    comprobarlas: catalogar el efectivo es **imposible** por
    `prefixItems`+`items:false`; las 21 mutaciones de divergencia no esconden
    ningún falso negativo —la prueba afirma primero que el paquete íntegro
    pasa—; la igualdad efectivo/manifiesto no es circular porque relee del
    disco; y la superviviente declarada es inalcanzable por construcción del
    tipo, exactamente donde se dijo.
  - **B-1: la puerta física picaba hasta 13,3 veces el presupuesto.** Ninguna
    de las tres pruebas que vigilaban el pico —escritor, repositorio y
    oráculo— cubría esta fase, de modo que el defecto era invisible; y la
    versión anterior de esta entrada afirmaba que el oráculo era la única
    pieza fuera de la fórmula, lo que era **falso**:
    `verify_published_f5_package` sí se ejecuta, y dos veces. Corregido en
    cuatro pasos, midiendo con el arnés `peak_bytes` después de cada uno:
    prestar las dos listas de subpoblaciones en vez de clonarlas
    (9,57 → 9,56×); **tipar la puerta** con `serde` para que
    `initialization_summaries` —el término que crece con `B*N_pop*S` y que,
    con `arbitrary_precision`, cuesta unas siete veces su texto— se **salte**
    en vez de materializarse (9,56 → **0,87×**); emitir esos mismos resúmenes
    por flujo al serializar y escribir el manifiesto **directamente al
    fichero**; y leer los dos documentos por flujo en la puerta. Resultado:
    la puerta pasa de **9,57× a 0,87×** con `B=32` y de **12,59× a 0,71×** con
    `B=128`.
    `the_verification_peak_fits_the_sealed_budget_and_does_not_grow_with_the_block_count`
    es la prueba que faltaba y mata la regresión con 9,29×.
  - **Lo que sigue sin caber, con su cifra.** La publicación **completa** pasa
    de **19,08× a 1,63×** con `B=32` —314 412 frente a 193 120— y a **1,36×**
    con `B=128` —763 759 frente a 562 432—. Lo que sobra es el objeto `f5` del
    manifiesto más lo que la puerta obliga a materializar para compararlo con
    el efectivo, ambos `O(B*N_pop)`: `sampling`, `f_realized` y `quota_ties`,
    que DTC-089 obliga a comparar término a término entre los dos documentos.
    `f5-resource-accounting-v8` **nunca modeló la construcción de los dos
    documentos**: su `final_formula_peak` presupuesta el escritor por flujo de
    F5.4.1, que cumple con 0,2-0,4×. La desviación queda acotada por
    `the_whole_publication_stays_within_twice_the_sealed_budget_without_the_oracle`,
    declarada aquí con su número y **elevada como bloqueante de F5.6**, donde
    coincide con la condición del oráculo.
  - **G-1: se publicaba un `code` que el propio esquema rechaza.**
    `f5-failure-v1.schema.json` cierra `code` a diez literales y el catálogo
    público de diagnósticos es más ancho: `GPS-INTERNAL-001`, `GPS-IO-001` y
    `GPS-FMT-006` son las causas **naturales** de la etapa `publication` y
    ninguna pertenece al enumerado. `from_diagnostic` tomaba el código sin
    filtrar y publicaba atómicamente, con destino final creado, un manifiesto
    que su propio esquema rechaza. Corregido con el tipo cerrado
    `F5FailureCode`: una causa ajena se publica como `GPS-INTERNAL-101` —que
    sí pertenece— **conservando el código original como código del
    diagnóstico**, y anotando bajo `published_failure_code` cuál se publicó en
    su lugar. La puerta comprueba además `code` y
    `stage` contra sus enumerados, y una prueba de esquema externo recorre las
    cinco etapas y la causa ajena.
  - **G-2: un general con un espacio en el nombre reventaba al final.**
    `require_relative_path` impone el alfabeto de `$defs/path` a los cinco
    enlaces declarados, pero el nombre del **propio** general lo elige el
    usuario en la línea de órdenes y nadie lo validaba: con `mi general.toml`
    se ejecutaban admisión, layout, búsqueda, fundadores, identidades, los
    ocho artefactos y el recuento MAN-70 **completos** antes del rechazo, con
    un diagnóstico que culpaba al esquema del documento producido. La frontera
    correcta es la admisión: `require_publishable_general_name` lo rechaza con
    `GPS-CONF-101` y `fields=general_path` antes de abrir el fichero, y el
    detalle dice qué renombrar.
  - **M-1 y las notas.** El invariante de `root_documents` solo se cumplía en
    una dirección: un artefacto podía reclamar `effective-config.json` ya
    escrito, sobrescribirlo y catalogarse; `check_new_artifact` lo consulta
    ahora y la prueba cubre las dos direcciones. `layout_attempts_used` se
    extrae a una función pura y se prueba, porque `F5BestLayout` no tiene
    constructor público y un fixture que agote el layout no sirve —el propio
    fixture redimensiona el FASTA con `block_count`, que es el falso negativo
    ya documentado—. `recount_matches` se emite desde el resumen y no como
    literal. `RECOUNT_SOURCE` se importa de `f5_oracle` en vez de duplicarse.
    `F5FailureStage::Preflight` se documenta como variante **sin productor
    posible**: DTC-089 prohíbe que un rechazo de preflight abra destino.
    `F5FounderRoster` pierde `Eq` —cambio de API pública— porque conserva un
    `f64`; `PartialEq` basta para todo su uso.
  - **Un hueco que la corrección cerró de paso.** Al tipar la puerta con
    `deny_unknown_fields` y campos obligatorios, el objeto `/f5` de los dos
    documentos queda **cerrado dentro del productor** y no solo en el
    validador externo: sus veinticinco y dieciocho claves exactas son ahora
    las de un tipo.
  - **Revisión F5.4 (cierre independiente) aplicada (2026-09-12).** Devolvió
    **NO CERRADO** por un motivo exacto: su criterio es «publicación atómica y
    oráculo independiente sobre éxito **y fallo**», y la mitad «fallo» no
    estaba demostrada **porque no estaba implementada**. Dos bloqueantes y
    cuatro menores, todos aplicados.
  - **H1: el paquete de fallo no tenía oráculo.** `audit_f5_package` exige los
    siete TSV y `truth.jsonl` y no tiene modo de fallo; y la puerta física
    leía `code`, `stage` y el enlace **del manifiesto** sin abrir nunca
    `outputs/f5-failure.json`. Reproducido: reescribiendo el diagnóstico con
    otra causa, otra etapa, `layout_attempts_used=999` y un artefacto
    biológico, y **rehaciendo el enlace**, la puerta devolvía `Ok(())`. Dos
    piezas: `verify_failure_diagnostic` promueve al productor el bucle de los
    **nueve términos comunes** que ya estaba escrito en la batería —igual que
    F5.4.2 promovió el lector de F5.4.1—, y `audit_f5_failure_package` es el
    gemelo del oráculo para el otro estado: doce claves, enumerados cerrados
    de código y etapa, `artifacts=[]`, los nueve términos contra la
    declaración, los contadores dentro de los límites efectivos, el árbol
    **enumerado** y el enlace recomputado.
  - **H2: un fichero que nadie declara llegaba al destino.** Ninguna de las dos
    ramas recorría el árbol: la de éxito no miraba y la de fallo comprobaba
    `exists()` sobre nueve rutas literales. Reproducido en los **dos** estados
    con un documento raíz `orphan.json`. `verify_payload_inventory` enumera la
    carga útil y exige igualdad exacta con el conjunto declarado: los siete
    TSV, la verdad y los dos documentos raíz en éxito; el manifiesto y el
    diagnóstico en fallo. No es hipotético: DTC-035 exige conservar la
    configuración original como copia TOML y en v5 esa copia **solo puede ser**
    un documento raíz, porque `artifacts` es `prefixItems`+`items:false`.
  - **H3: la etapa `publication` no tenía productor posible, y ahora lo tiene.**
    DTC-045 lo resuelve sin ambigüedad: «Tras comenzar, un fallo controlado …
    publica atómicamente el directorio final como fallido». Un recuento MAN-70
    rechazado es un fallo controlado posterior al comienzo, y con la firma
    anterior el `ExecutionWorkspace` entraba por valor: al devolver `Err` ya
    estaba destruido y el llamador se quedaba sin con qué publicar, dejando
    solo el provisional `INCOMPLETE` que DTC-045 reserva a la **interrupción no
    controlable**. `publish_f5_success` devuelve ahora
    `Result<PathBuf, Box<F5PublicationFailure>>`, que **entrega el provisional
    vivo**. Se eligió el tipo suma y no `&mut ExecutionWorkspace` porque hace
    explícita la entrega en la firma: con `&mut`, el llamador se quedaría un
    provisional a medio escribir sin ninguna señal de que le toca cerrarlo.
    **Qué no puede publicar un fallo, escrito y no supuesto**: `workspace` es
    `None` cuando el renombrado ya había empezado —destino existente o sistema
    de ficheros sin garantía no-clobber— y cuando publicar el fallo vuelve a
    fallar; un `ENOSPC` real no deja escribir ni el diagnóstico. Esos casos
    **son** la interrupción de DTC-045. Esto da sentido retroactivo a G-1: el
    tipo `F5FailureCode` protegía una etapa que ningún llamador podía producir.
  - **H4: tres términos sin juez.** `eligible_bases`, `attempt_summary` y
    `quota_ties` no estaban en las dieciocho claves del efectivo ni en la
    declaración. Ahora: `eligible_bases` queda acotado por arriba en la puerta
    —no puede exceder `resource_preflight.reference_bases_processed`, que el
    propio manifiesto publica— y por abajo en el oráculo —no puede ser menor
    que las bases que el layout colocó, que el oráculo suma de
    `f5-block-layout.tsv`—; `attempt_summary` viaja a la declaración y el
    oráculo contrasta sus dos términos derivables contra las filas de
    `f5-attempts.tsv`, dejando `layout_attempts_used` declarado porque ningún
    TSV lo publica. `quota_ties` queda **declarado no observable**, como en
    F5.4.2, y ahora también en `f5_publication.rs`.
  - **H5: `layout_attempts_used` publicaba un dato falso.** La regla decía que
    cualquier causa distinta del agotamiento aborta antes de dejar un mejor
    intento, y es **falsa**: `consider_best` corre **antes** de
    `materialize_attempt` en la misma iteración. Lo que descartaba el candidato
    era `impl From<Diagnostic> for F5LayoutFailure`, que fijaba
    `best_layout: None`, de modo que un fallo en el intento 3 de 10 publicaba
    `layout_attempts_used = 0`. `F5LayoutFailure` conserva ahora
    `attempts_used` y el informe lo lee en vez de deducirlo. La prueba tautológica —cuatro `assert_eq!` sobre una
    función pura de dos argumentos, séptima aparición del patrón— se sustituye
    por un agotamiento **real**: hueco mayor que el bloque, ningún intento
    coloca dos sitios, `GPS-DOMAIN-106` y los ocho intentos publicados.
  - **H6: `validate_rng_catalog` solo exigía que el catálogo no estuviera
    vacío** mientras el reparto 3+9 de ese módulo se sostiene sobre doce filas.
    Cierra con `seen != F5_RNG_CATALOG_DECISIONS`, y una prueba contrasta la
    constante contra el fichero publicado.
  - **La costura del oráculo, cerrada.** `publish_f5_success` entregaba al
    oráculo una declaración construida desde la memoria del productor y luego
    armaba el manifiesto **por segunda vez** desde las mismas estructuras: el
    oráculo aprobaba una construcción paralela a la publicada. Ahora la
    declaración se construye **una vez** y el objeto `f5` se **deriva de
    ella**. Los tres únicos términos que no puede transportar —
    `inbreeding_f_source` y el `block_index` y `source` de cada fila de
    muestreo— siguen viniendo del censo y de la resolución, porque el contrato
    de `F5PackageDeclaration` es de F5.4.2 y este hito no lo reabre; queda
    dicho término a término en el código.
  - **La demostración de atomicidad es de frontera de biblioteca.** Todo lo que
    la Revisión F5.4 dio por demostrado —`RENAME_NOREPLACE` sin fallback,
    `OpenArtifact::Drop`, `INCOMPLETE` en el envoltorio y nunca en el árbol
    renombrable, la revalidación previa al renombrado— se prueba **con
    fixtures sobre la frontera de biblioteca**: `application.rs` y `main.rs` no
    invocan ninguna de las cuatro fronteras F5, y las únicas pruebas de
    interrupción y concurrencia **entre procesos** recorren F1. **F5.6 debe
    repetirla sobre el proceso**, con la interrupción y la concurrencia que
    DTC-045 exige probar en WSL2/ext4 y Ubuntu nativo.
  - **Dictamen sobre las dos piezas fuera de v8.** La revisión independiente
    concluye —y se registra así— que **es v8 el que está incompleto**, no las
    fases las que deben caber: los dos excesos los crea el propio contrato,
    porque `final_formula_peak` presupuestó el escritor por flujo de F5.4.1 y
    **nunca modeló** ni la construcción de los dos documentos ni un lector del
    paquete. La recomendación para cuando F5.6 la tome es
    **`f5-resource-accounting-v9` con dos términos nuevos**, uno por cada
    pieza, con sus órdenes ya medidos: la construcción y comparación de los dos
    documentos crece con **`O(B*N_pop)`** —`sampling`, `f_realized` y
    `quota_ties`, que el contrato obliga a comparar término a término— y el
    lector del paquete está acotado por **el registro más largo**, no por el
    paquete, con la razón convergente que F5.4.2 midió. Y con su coste dicho:
    v9 **reabre los dos esquemas**, porque el literal
    `"f5-resource-accounting-v8"` es un `const` en
    `effective-config-v5.schema.json` y, por referencia, en
    `manifest-general-v5.schema.json`. Por eso no se hizo aquí.
  - **Comprobación de cierre: cuatro flecos (2026-09-13).** Da **F5.4 cerrado**
    tras reproducir sus **propios** ataques —el documento de fallo reescrito
    con otra causa, otra etapa, otros contadores y un artefacto biológico, con
    el enlace rehecho, muere dos veces: en la puerta con `GPS-FMT-006` y en el
    oráculo con causa `Link`; el `orphan.json` cae en los dos estados citando
    la ruta; y un fichero oculto en la raíz, otro dentro de `outputs/` y un
    subdirectorio inesperado caen los tres—. Los picos, además, **mejoraron**
    pese a recorrer el árbol y parsear un documento más: la puerta mide 0,79×
    con `B=32` y 0,62× con `B=128`, frente a 0,87× y 0,71×.
  - **N1: la afirmación anterior de H5 era falsa y queda corregida.** Se
    escribió que el índice se conservaba «en todas sus vías, `From<Diagnostic>`
    incluida», y **no era así**: esa conversión seguía fijando cero y seguía
    siendo alcanzable desde **diez `?` del cuerpo del bucle**. La conversión
    implícita **se ha retirado**: el cuerpo vive ahora en un cierre que
    devuelve `Diagnostic` y una sola conversión le añade lo consumido, de modo
    que ninguna vía del bucle puede volver a perderlo; los tres únicos sitios
    anteriores al primer intento usan `before_any_attempt`, que declara el cero
    como verdad y no como comodín. Dos pruebas nuevas cubren las dos mitades.
    **Un matiz medido que corrige al revisor**: el error de dominio «residual
    final fuera del soporte factible» **no** es alcanzable desde una
    configuración admitida —cinco fixtures con `covered_bases` no
    representable mueren antes con `GPS-DOMAIN-101`, y la aritmética de
    intervalos del reparto condicionado mantiene el residuo factible por
    construcción—, de modo que el cero en producción era alcanzable por las
    guardas de desbordamiento, no por ese camino. La corrección es la misma y
    es estructural, que es lo que importa.
  - **N2: la comparación catálogo/manifiesto se quedó sin prueba.** Efecto
    colateral de H2: el inventario adelantó a la única que la ejercitaba. No es
    código muerto —el inventario compara **rutas** y el catálogo compara
    **(id, ruta)**—, de modo que se añade la vigésimo segunda divergencia: un
    `artifact_id` ajeno sobre la **misma ruta**, que solo el catálogo ve.
  - **N3: la frontera de «no se puede publicar el fallo» era más ancha que su
    docstring.** `ExecutionWorkspace::publish` toma `self` por valor, así que
    **cualquiera** de sus comprobaciones previas al renombrado —catálogo,
    marcador `INCOMPLETE` y la **segunda pasada** de la puerta— devuelve
    también `None`. Ese último caso es justamente un `GPS-FMT-006` de etapa
    `publication`, y no poder publicarlo es defendible: significa que el
    provisional se alteró entre dos revalidaciones y ya no se sabe qué se
    declararía. Queda escrito.
  - **N4: el inventario seguía enlaces simbólicos.** Decidía con `is_dir()`,
    que los sigue, y la recursión no tiene tope: un enlace a un ancestro la
    colgaría. Tipa ahora con `entry.file_type()` y rechaza lo que no sea
    fichero regular o directorio, en la puerta y en el oráculo. Lo que cambia
    el veredicto es el tipado; la guarda cambia la **razón**, y así se dice.
  - **Un falso negativo evitado.** La primera versión de la prueba de rechazo
    de preflight usaba un `block_count` mayor esperando `GPS-DOMAIN-101`; el
    fixture **redimensiona el FASTA con el número de bloques**, de modo que el
    objetivo siempre cabía y la admisión no rechazaba nada. La prueba habría
    pasado solo cuando el `block_count` era tan grande que excedía el
    presupuesto, midiendo una causa distinta de la anunciada. Se sustituyó por
    dos rechazos de admisión explícitos y de causas distintas —presupuesto
    excedido e identidad de perfil que no coincide—, ambos con destino y
    provisional inexistentes.

> **Revisión F5.4.** No continuar sin demostrar publicación atómica y oráculo
> independiente sobre éxito y fallo.

## 8. F5.5 — auditoría de identidad y publicación cerrada

- ✅ **F5.5.1.** Auditar que el paquete publicado solo lee identidad, copias,
  pertenencia y asignaciones ya cerradas en F5.3; ningún paso F5.5 puede crear,
  sustituir o modificar esos datos.
  - Depende de: F5.4.3.
  - Procede de: REQ-EVO-32/53/74/77, REQ-POB-01 y arquitectura §15.
  - Cambia: evidencia y validación, no el estado biológico ni el paquete atómico.
  - Se comprueba con: mutación posterior prohibida, pertenencia ausente/duplicada,
    IDs aportados/derivados, igualdad FCP/FHA y ausencia de migración, selección
    o genealogía.
  - **Interpretación elegida: las dos mitades, separadas por su naturaleza, y
    ninguna publica nada.** El hito pide dos cosas distintas y solo una de
    ellas puede salir de los bytes. (a) Que el paquete **coincida** con las
    identidades que F5.3 cerró y que ninguna se creara después es una
    propiedad **temporal**: un fichero no atestigua en qué orden se derivó lo
    que contiene, de modo que la precedencia solo se demuestra relacionando el
    repositorio cerrado con los bytes publicados. (b) Que **ninguna ruta de
    F5.5 pueda escribir** en esos datos es una propiedad **del tipo y del
    recorrido**: la auditoría recibe `&F5FounderIdentities` y derivar exige
    `&mut self`, de modo que lo impide el comprobador de préstamos, no una
    promesa. Se implementan las dos. Lo que **no** se hace es repetir F5.4.2:
    el oráculo ya rederiva cada `FCP_`, `HAP_` y `FHA_` desde los bytes de
    forma independiente, y volver a validarlos aquí sería el oráculo circular
    que esta base de código ya ha pagado siete veces.
  - **Dónde vive la evidencia: en una comprobación que no publica nada.** No
    hay octavo artefacto —`artifacts` es un `prefixItems` de siete con
    `items:false`—, no se toca el manifiesto y no se reabre ningún esquema. La
    auditoría devuelve un `F5ClosureReport` en memoria y la batería lo fija.
    El paquete atómico queda intacto: una prueba mide el árbol publicado byte
    a byte antes y después de auditarlo, y valida el instrumento aparte
    cambiando un byte sobre una copia.
  - **Módulo y frontera.** `src/f5_closure.rs`, con `seal_and_publish_f5_success`,
    `seal_and_publish_f5_failure`, `F5ClosureWitness`, `audit_f5_identity_closure`
    y `audit_f5_failure_closure`, y la batería `tests/f5_closure.rs`. Ningún
    otro módulo cambia salvo el registro en `src/lib.rs`: F5.4.3 se consume,
    no se modifica.
  - **Corrección de la Revisión F5.5: un sello suelto no atestiguaba nada, y
    ahora hay testigo.** La primera versión guardaba cuatro cifras del
    repositorio sin consumir y las comparaba después de publicar. La revisión
    adversarial demostró que eso no prueba precedencia: las cuatro cifras son
    una **función determinista del repositorio sin consumir**, de modo que
    rematerializar un repositorio equivalente **después** de que el paquete
    existiera, sellarlo entonces y consumir su recorrido con un visitante vacío
    avalaba el paquete viejo. La propiedad temporal era cierta por construcción
    —`for_each_assignment` es el único `&mut self` y la publicación lo
    consume—, pero **la auditoría, que es el producto del hito, no la
    atestiguaba**. Ahora la auditoría exige un `F5ClosureWitness`, que solo
    nace de `seal_and_publish_f5_success` o `seal_and_publish_f5_failure`
    —sellan y **luego** publican— y que guarda el destino publicado. Un
    repositorio reconstruido a posteriori no puede acreditar un paquete ya
    publicado: obtener un testigo obliga a publicar, publicar obliga a un
    destino nuevo —`ExecutionWorkspace` rechaza uno existente— y el testigo
    nombra ese destino nuevo. El constructor del sello pasa a ser privado del
    módulo.
  - **Adaptación declarada: no se reabre la frontera de F5.4.3.** La revisión
    pedía que `publish_f5_success` recibiera y devolviera el testigo. Eso
    cambiaría su firma pública y la de todos sus llamadores, que son de un hito
    ya cerrado. Se obtiene la **misma** semántica —sellar, publicar, devolver
    el testigo que ata las dos cosas— con dos envoltorios en `f5_closure` que
    delegan en `publish_f5_success` y `publish_f5_failure` sin tocarlos ni
    cambiar sus errores. La garantía es idéntica porque el testigo no tiene
    constructor público: solo esos dos envoltorios lo crean, y ambos publican.
  - **Ley de conservación, en sus dos mitades.** En el éxito, entre el sello y
    la auditoría el registro solo puede crecer en **exactamente una `FHA_` por
    fila de `f5-founders.tsv`**: las que deriva la pasada única de
    `for_each_assignment`. Ni una más —sería una identidad creada después de
    F5.3— ni una menos —sería una fila publicada sin derivar—. No es una
    comparación de valores, de modo que no puede ser circular.
  - **Corrección de la Revisión F5.5: la mitad «fallo» no aplicaba su ley.** La
    primera versión solo comprobaba que el registro no **encogiese**, y el
    registro nunca encoge: derivar el recorrido completo después de publicar un
    fallo pasaba la auditoría, y la única aserción que sostenía esa mitad
    miraba el repositorio, no el veredicto. Ahora la mitad «fallo» exige
    **igualdad estricta** —publicar un fallo no deriva ninguna identidad— con
    causa `Conservation`, y tiene prueba que muere si se desactiva. La rama del
    registro encogido, que era inalcanzable, se retiró.
  - **Qué se comprueba contra los bytes.** Cada fila de `f5-founders.tsv`
    resuelve su `FCP_` en el repositorio y su tupla —pertenencia, individuo,
    secuencia e índice de copia— tiene que ser la cerrada; su `HAP_` tiene que
    ser uno de los que F5.3 cerró para ese bloque×subpoblación; las tres clases
    de identificador ocupan su propio lugar, que es la mitad de «igualdad
    FCP/FHA» que un paquete puede romper. El conjunto de copias publicadas es
    **exactamente** el cerrado, y `outputs/truth.jsonl` —único artefacto con
    `founder_membership`— publica esa misma pertenencia, completa y sin
    repeticiones.
  - **IDs aportados y derivados.** El censo cerrado dice cuál de las dos vías
    resolvió cada `IndividualId`, y el informe publica los dos recuentos. La
    batería los contrasta contra la declaración del fixture y, aparte, contra
    el prefijo de clase `IND_` de los derivados.
  - **Ausencia de migración, selección y genealogía.** Se comprueba en dos
    planos. Sobre bytes, la cabecera exacta de `f5-founders.tsv`: el fichero
    que sostiene la identidad fundadora no admite una columna de parental, de
    origen migratorio ni de eficacia biológica porque no admite ninguna columna
    más. Sobre el texto del módulo, una **guarda léxica declarada como tal**,
    con control positivo, igual que `f5_never_computes_a_reproductive_ne`.
  - **Memoria: no hay tercera pieza fuera de `f5-resource-accounting-v8`.** La
    auditoría sostiene **un registro cada vez** y no retiene nada proporcional
    a `B*C_total`, al censo ni al número de identidades. Los registros de
    verdad se deserializan contra una estructura de cuatro campos y no contra
    un árbol JSON —lo que descarta `segments`, `history` y `current_pieces`
    conforme se leen—; la primera versión, que sí materializaba el árbol,
    picaba 132 566 bytes frente a los 100 792 sellados y **no cabía**. La
    versión final pica 23 120 bytes con `B=8` y 70 672 con `B=32`, frente a
    100 792 y 193 120 sellados, de modo que cabe dentro del presupuesto de la
    admisión y no necesita ninguna declaración fuera de la fórmula. La prueba
    acota además el pico contra el mayor registro publicado.
  - **Cada prueba discrimina, y se midió dos veces.** En la primera ronda se
    desactivaron doce comprobaciones; la Revisión F5.5 midió veintinueve y
    encontró **diecisiete sin prueba**, cuatro de ellas nombradas por este
    mismo plan como la evidencia del hito. La segunda ronda cubre **treinta
    mutaciones con una sola superviviente**, la precondición declarada abajo.
    La tabla completa vive en el encabezado de `tests/f5_closure.rs`. Se
    añadieron pruebas discriminantes para la cabecera exacta, las clases `FCP_`
    y `HAP_`, la tupla `sequence_id`/`copy_index`, la pertenencia ausente en el
    registro, la unidad no cerrada, la copia no cerrada en la verdad, la
    identidad de la verdad, la copia repetida, el avance de bloques, el LF
    final, las columnas de más y de menos, el documento vacío, el documento sin
    cabecera y las dos ramas del sello.
  - **Lo que el estudio corrigió de la propia batería.** En la primera ronda,
    **cuatro pruebas pasaban por una comprobación distinta de la que
    anunciaban** —la de pertenencia por el orden, la de renombrado por la
    búsqueda en el censo—: el orden se juzga ahora **después** de las
    comparaciones de identidad y toda prueba fija también el campo y un
    fragmento del detalle. En total se retiraron **cuatro aserciones que no
    podían fallar**: la igualdad de individuos contra la suma de
    `individual_count`, los dos ceros del informe de fallo —sustituidos por el
    contraste contra un paquete de éxito del mismo fixture—, la tautología
    `aportados + derivados == individuos` y el recuento de columnas implicado
    por la igualdad literal de la cabecera.
  - **Las corrupciones se aplican sobre el destino publicado, no sobre una
    copia.** Es consecuencia del testigo —acredita un destino concreto— y es
    además la lectura correcta del hito: alterar el paquete ya publicado *es*
    la mutación posterior que F5.5 prohíbe.
  - **Una precondición inalcanzable, declarada.** `check_repository_order`
    comprueba el orden canónico de unidades antes de su búsqueda binaria. Desactivarla no mata ninguna prueba porque ninguna configuración
    admitida la alcanza. Se conserva porque lo que protege es un **falso
    negativo de la propia auditoría**: una búsqueda binaria sobre datos
    desordenados devolvería «ausente» sobre identidades presentes y la clausura
    rechazaría un paquete correcto citando una causa falsa. Queda escrito en
    vez de disimulado, igual que las dos ramas defensivas del testigo.
  - **Sin llamador en producción, como el oráculo.** F5 todavía no tiene
    orquestador —es F5.6—, de modo que `audit_f5_identity_closure` se ejercita
    desde la batería exactamente igual que `audit_f5_package`.
  - **Ámbitos cubiertos.** Ploidías 1/2/3, modos probabilístico y controlado,
    `F` en 0.0 y 1.0, identidades aportadas y derivadas, dos subpoblaciones,
    dos secuencias pobladas más una estéril, y paquete de éxito y de fallo.
  - **Corrección de la Revisión F5.5 (cierre): el testigo sellaba el
    repositorio pero no el censo.** Era la pieza que **ningún auditor
    cubría** —el censo no está en el paquete, de modo que el oráculo de F5.4.2
    no puede verlo— y la clausura lo recibía como cuarto parámetro suelto.
    Publicando con `founder_keys` y auditando con una segunda cadena que
    declarase `founder_ids` literales iguales a los `IND_` derivados, los dos
    recuentos que este plan nombra como evidencia del hito **se invertían** y
    la auditoría devolvía `Ok`. La clase fundadora de cada subpoblación se
    congela ahora en el testigo al sellar, tomada de la configuración que
    publica y solo si `matches_admission` la admite —`SealedF5Admission` cubre
    `subpopulation.founder_kind` en su material canónico—, y **el censo deja de
    ser un parámetro de la auditoría**: la sonda ya no es representable, que es
    más fuerte que rechazarla. Con el censo desaparece también la mitad de
    `check_repository_order` que validaba su orden, que protegía un objeto que
    no era el que publicó. `F5PublicationRequest` gana dos accesores de solo
    lectura —`config()` y `sealed()`— y nada más de F5.4.3 cambia.
  - **Corrección de la Revisión F5.5 (cierre): la clausura de fallo acreditaba
    un paquete inexistente.** No abría ni un byte: el sello miraba memoria, la
    acreditación caía a comparación literal cuando ninguna ruta canonicalizaba
    y el bucle de artefactos biológicos estaba construido sobre una **ausencia**
    que un destino borrado satisface trivialmente. Con `rm -rf` del destino, y
    también borrando solo el documento de fallo, devolvía `Ok`. Ahora la
    canonicalización es **obligatoria** en las dos mitades —una ruta que no
    resuelve es `Missing`, que es lo que es— y la mitad «fallo» abre
    `outputs/f5-failure.json` sin leerlo: una afirmación positiva sobre el
    disco. Juzgar su contenido sigue siendo de `audit_f5_failure_package`.
  - **`F5ClosureCause::Missing` pasa a tener pruebas.** No tenía ninguna
    aparición en la batería y es alcanzable por la mutación posterior más
    canónica que existe: borrar un artefacto publicado. Tres pruebas nuevas —el
    TSV fundador, la proyección de verdad y el destino entero— la fijan.
  - **Dos afirmaciones rebajadas a lo que es cierto.** El contrato de LF **no**
    estaba unificado: `walk_truth_membership` rompe en el primer genotipo y la
    última línea del documento es siempre un genotipo, de modo que el LF final
    de `truth.jsonl` nunca se mira. La afirmación pasa a ser «los dos recorridos
    exigen LF en cada línea que leen»; cerrar el recorrido costaría una pasada
    entera por un contrato que la puerta física y el oráculo ya sostienen. Y la
    cota diferencial del pico deja de apoyarse en un `+8192` que dejaba pasar
    8 KiB sin explicar: con el par 8→16 copias, que no atraviesa ninguna
    duplicación de búfer, se exige la diferencia del mayor registro **sin
    holgura** —y medido, el pico ni sube: baja 178 bytes mientras el censo se
    dobla—.
  - **Guarda léxica de «no escribe», ampliada.** Enumeraba seis vías; enumera
    ahora trece, con `File::options`, `fs::copy`, `fs::hard_link`,
    `fs::set_permissions`, `::symlink(` y `Command`. Lo que sostiene la
    propiedad de verdad sigue siendo la huella byte a byte, que sí es completa.
  - **Cuatro notas menores de la Revisión F5.5, aplicadas.** El contrato de
    bytes se **unifica**: los dos recorridos exigen LF final, como el oráculo,
    con una sola función. La primera línea de `truth.jsonl` tiene que **ser**
    la cabecera v2 y no solo existir: si faltara, la primera copia se
    consumiría en su lugar y el rechazo saldría con una causa falsa, que es el
    defecto ya corregido en el otro recorrido. El `allow(too_many_lines)` lleva
    justificación y condición de retirada. Y la prueba de pico mueve ahora
    **las dos** magnitudes por separado —bloques con censo fijo y censo con
    bloques fijos—: multiplicar por ocho las copias mueve el pico solo lo que
    se movió el mayor registro publicado.
  - **Nota fuera de alcance.** `src/f5_founders.rs:2285` documenta
    `GPS-CONF-101` para «el recorrido ya se consumió» y el código devuelve
    `GPS-INTERNAL-101`. La prueba de este hito fija el código real; el docstring
    de F5.3 es el que no corresponde. No se toca por la prohibición de revertir
    cambios ajenos: queda anotado para quien reabra F5.3.
  - **Obligación de proceso, escrita en F5.6.** La garantía «para auditar hay
    que haber publicado» basta para la precedencia y no como obligación de
    proceso: F5.6.1 y F5.6.3 recogen ahora que el orquestador publique
    **exclusivamente** por los envoltorios que sellan y componga los dos
    auditores sobre el destino ya renombrado.
  - **Puertas.** `fmt`, `clippy -D warnings`, `test --workspace --all-targets`
    (1076 pasan, 0 fallan, 12 ignoradas, frente a 1070/0/12 al cerrar el hito y
    1029/0/12 antes de empezarlo), `test --doc`, `build --release` y
    `doc --no-deps`, todas en verde.

> **Revisión F5.5.** Confirmar que identidad y pertenencia precedieron a la
> publicación y que el hito no modifica el paquete atómico.

## 9. F5.6 — Marimo y cierre

- ✅ **F5.6.1.** Añadir los controles F5 a Marimo sin panel ni recombinación y
  ejecutar el mismo `validate`/`run` Rust.
  - Depende de: F5.5.1.
  - Procede de: REQ-UI-15, ACE-114/115/117 y DTC-051/089.
  - Cambia: pestaña poblacional y adaptador de proceso.
  - Se comprueba con: matriz campo→control, importación/exportación de exactamente
    un general v3+escenario v1, rechazo de cero/varios/v2-v4, edición de
    sustituciones modo/tolerancia y presupuestos, invalidación tras editar y
    paridad de efectivo, diagnóstico y hashes con TOML manual.
  - **Obligación de proceso heredada de la Revisión F5.5.** El adaptador de
    proceso debe publicar **exclusivamente** por `seal_and_publish_f5_success` y
    `seal_and_publish_f5_failure`, nunca por `publish_f5_success` ni
    `publish_f5_failure` directamente. La garantía de F5.5.1 es «para auditar
    hay que haber publicado», que basta para la precedencia y **no** basta como
    obligación de proceso: un orquestador que llame a las funciones de F5.4.3
    —que siguen siendo `pub`, y así se quedan: bajarlas a `pub(crate)`
    reabriría la frontera de un hito cerrado— obtiene un destino **sin
    testigo** y deja la clausura inaplicable sin ninguna señal.
  - **Implementado (2026-09-13), pendiente de revisión adversarial.** La casilla
    sigue en ☐ hasta que la revisión la apruebe.
  - **Decisión material: `run` F5 queda cableado.** Antes de este hito
    `src/main.rs` devolvía «run F5 todavía no está habilitado» y ninguna de las
    cuatro fronteras F5 tenía llamador productivo. La comprobación que el hito
    pide —«paridad de **efectivo**, diagnóstico y **hashes** con TOML manual»—
    no es formulable sin `effective-config.json`, `manifest.json` ni los ocho
    artefactos, de modo que encadenar la ejecución es parte del hito y no un
    añadido. `application.rs` gana `execute_validated_f5_config`
    —admisión → layout → búsqueda → asignación → identidades → publicación— y
    `main.rs` enruta `(Run, F5)` a ella.
  - **Material del `execution_id` de F5, que la Revisión F5.1 dejó abierto
    «al encadenar la cadena productiva».** Es la concatenación de los SHA-256
    **observados** de los dos documentos de autoría, en orden de rol (general y
    después escenario), pasada por `configuration_digest` y por
    `derive_execution_id`, que es el mismo camino que F1–F4 usan con su único
    documento. No puede ser solo el general: el general enlaza el escenario
    **por ruta y sin hash**, así que dos escenarios distintos compartirían
    identidad de ejecución y, con ella, `block_id` y `site_id`. La prueba
    `the_execution_identity_covers_both_authored_documents` lo fija y su
    mutación —suprimir el término del escenario— produce exactamente el
    layout **idéntico** que la Revisión F5.1 anticipó.
  - **DTC-045: las cinco etapas publican `failed`.** La primera redacción abría
    el provisional **después** de las cuatro etapas en memoria y dejaba sin
    destino los fallos de layout, búsqueda y fundadores. La revisión adversarial
    lo rechazó con razón: definir «comenzar» como «cuando este código abre el
    provisional» hace la regla **autocumplida** y dependiente del orden de las
    líneas, la Revisión F5.4 ya resolvió el caso análogo —un recuento MAN-70
    rechazado **sí** publica `failed`— y esas etapas ya han abierto ámbitos RNG
    derivados del `execution_id`: la simulación ha comenzado. Ahora el
    provisional se crea **antes** del layout y `F5FailureStage` tiene productor
    en cuatro de sus cinco variantes: `Layout`, `Haplotypes`, `Founders` y
    `Publication`. **`Preflight` queda sin productor por prohibición del
    contrato, no por omisión**: DTC-089 impide que un rechazo de preflight abra
    destino, de modo que darle uno sería romper el contrato y no completar una
    cobertura. Queda escrito para que nadie la «arregle» más adelante.
    Para ello `f5_founders.rs` gana `empty_f5_founder_identities`, un
    repositorio **vacío** ligado a configuración y ejecución: es **aditivo**, no
    deriva ninguna identidad, no toca el registro y no revierte nada de F5.3.
    La ley de conservación del fallo —publicar un fallo no deriva ninguna
    identidad— se cumple por construcción.
  - **Los dos auditores, compuestos.** `audit_f5_package` entra como el
    recuento MAN-70 que `publish_f5_success` exige, sobre la raíz de carga útil
    que `publish` renombra **sin tocar** al destino; `audit_f5_identity_closure`
    corre después sobre el destino ya renombrado, que es el que su testigo
    acredita. En la mitad «fallo» el gemelo es `audit_f5_failure_closure`.
    `audit_f5_failure_package` necesita una `F5FailureDeclaration` que el
    orquestador no conserva: queda para F5.6.3, que ya la pide.
  - **Obligación de proceso, comprobada con control positivo.**
    `the_orchestrator_publishes_only_through_the_sealing_wrappers` busca
    `publish_f5_success(` y `publish_f5_failure(` **no precedidos por
    `seal_and_`** —una búsqueda ingenua daría siempre positivo, porque el
    nombre del envoltorio contiene el del envuelto— en `application.rs` y
    `main.rs`, con dos controles positivos: que el orquestador sí nombra los
    envoltorios y que la **misma** búsqueda encuentra la llamada desnuda donde
    vive, en `f5_closure.rs`.
  - **Componente Python.** `ui/f5_editor.py` (nuevo) modela el par: tipos,
    enumerados, las tres expresiones regulares de los esquemas y una tabla de
    cotas enteras que una prueba contrasta contra `minimum`/`maximum` de los
    dos ficheros JSON Schema, para no duplicar reglas sin juez (DTC-051).
    **El renderizado se escribe a mano y no con `tomlkit`**: los bytes
    exportados son parte del contrato —el nombre del general se publica en
    `effective-config.json#/source_path` y la identidad de ejecución sale de
    los SHA-256 de los dos documentos—, de modo que el orden de claves, los
    espacios y los literales decimales tienen que quedar bajo una sola función
    que una prueba pueda fijar contra un TOML manual. `ui/process_adapter.py`
    gana `write_f5_pair`, `validate_f5_pair`, `run_f5_pair`, la sesión
    `ValidatedF5Session` —que conserva los **bytes** de los dos documentos y el
    nombre del general, no un resumen— y `f5_package_summary`, que resume el
    paquete leyendo solo el manifiesto acotado. Los dos documentos se escriben
    en el directorio de autoría y no en un temporal anónimo: sus enlaces son
    rutas relativas que Rust resuelve desde ahí, y el nombre del general se
    publica.
  - **Pestaña poblacional.** Las nueve pestañas mínimas de ACE-114 no cambian:
    el escenario F5 entra en «Poblaciones» con sus veinticinco controles, la
    importación de pareja, los dos documentos exportables y los botones de
    `validate`/`run` propios. `F` solo se edita con ploidía 2 y se retira del
    documento cuando la ploidía deja de serlo; `controlled` se declara
    desactivado cuando `F` es mayor que cero. No hay panel, haplotipos
    suministrados, construcción recombinante ni mapas, y una prueba lo busca
    con control positivo.
  - **Dos pruebas ajenas actualizadas, y por qué no es revertir.**
    `tests_ui/test_ui_contract.py` listaba `scenario_id` como «opción futura
    ausente»: F5.6.1 es el hito que la incorpora, de modo que pasa a exigir su
    control y conserva prohibidos panel, construcción recombinante y mapas.
    `tests/cli_config_diagnostics.rs` fijaba que `run` F5 devolvía
    `GPS-CONF-101` con el texto «la publicación pertenece a F5.4»; ahora fija
    que publica los diez ficheros, que `source_path` es el **nombre** del
    documento aunque se invoque por una ruta relativa desde un directorio
    ajeno, y que un segundo `run` al mismo destino se rechaza con `GPS-IO-003`
    sin sobrescribir. Su perfil de fixture pasa de `target_mean_r2=0.5` a `1`
    para que el escenario sea **factible**: la admisión que esas pruebas miden
    no depende de ese valor. El espectro **no** se toca —la primera redacción
    afirmaba que también hacía falta y era falso: comprobado aislando el
    cambio, con `[0,1)` y `target_mean_r2 = 1` las siete pruebas pasan—.
  - **Medida del primer `run` productivo** (binario `--release`, WSL2/ext4):
    `B=2,N_pop=1` 7,4 MiB RSS y <10 ms; `B=8,N_pop=2` 7,5 MiB y 20 ms;
    `B=16,N_pop=2` 8,1 MiB y 40 ms, con paquetes de 108 KiB, 372 KiB y 668 KiB.
    Las **dos piezas fuera de `f5-resource-accounting-v8`** —la construcción y
    comparación de los dos documentos, `O(B*N_pop)`, y el lector del paquete,
    acotado por el registro más largo— caben con holgura dentro del presupuesto
    de memoria de trabajo sellado (1 GiB por omisión) en los tres tamaños. La
    deuda de v9 **no se toca**: reabriría `effective-config-v5.schema.json` y
    `manifest-general-v5.schema.json`, que este hito tiene prohibido. Se
    observa además que el libro mayor v8 rechaza un escenario sobredimensionado
    **antes** de abrir la referencia (`GPS-RES-101`, fase
    `reference_and_mask_scan`), que es la guarda que se vuelve exigible al
    haber `run`.
  - **Deuda heredada que sigue abierta.** La demostración de atomicidad es de
    **frontera de biblioteca**: este hito añade la primera prueba de proceso
    —destino existente rechazado sin sobrescribir— pero **no** cubre
    interrupción ni concurrencia entre procesos, que DTC-045 exige probar y que
    hoy solo recorren F1. Queda para F5.6.3, que ahora sí puede hacerlo porque
    existe un `run` F5.
  - **Revisión adversarial F5.6.1 aplicada (2026-09-13).** Devolvió
    **integrable con condiciones**, con un bloqueante de pérdida de datos.
    - **B1 — el par F5 machacaba un general v2 ajeno.** `_replaceable_f5_document`
      aceptaba por `config_kind = "general"`, que es **también** el de los
      generales v2 que las otras ocho pestañas editan; con el destino por
      omisión `<repo>/general.toml`, la primera validación F5 lo destruía. Ahora
      la cabecera admisible es `config_kind` **y** `schema_version`, las dos
      primeras líneas exactas que emite el renderizador, en
      `F5_DOCUMENT_HEADERS`. La prueba colocaba el documento ajeno en la ruta
      del **escenario**, donde el rol no coincidía y no cazaba nada; ahora cubre
      el general v2 en la ruta del general.
    - **G2 — importar y exportar no conservaba el TOML.** El renderizador impone
      un orden de claves único, de modo que un escenario legal en otro orden no
      volvía igual y, como el `execution_id` sale de los SHA-256 observados, la
      misma autoría dejaba de ser la misma ejecución. `F5Pair` conserva ahora
      los **bytes importados** y `f5_pair_documents` los reemite mientras sigan
      describiendo el formulario vigente —lo comprueba **reparseando el
      origen**, no lo supone— y cede al renderizador en cuanto un control
      cambia. Dos órdenes legales distintos quedan fijados en la batería.
    - **G3 — la rama que publica el `failed` no tenía prueba de comportamiento.**
      Ahora hay dos: una inyecta un recuento MAN-70 que rechaza —por
      `execute_validated_f5_config_with_recount`, instrumental y `doc(hidden)`,
      que existe porque esa rama no se alcanza de otro modo sin corromper el
      provisional desde fuera del proceso— y otra agota la búsqueda haplotípica
      con un objetivo de LD inviable. Las dos comprueban el destino publicado,
      `state="failed"`, la etapa, la ausencia de los ocho ficheros biológicos y
      que el diagnóstico devuelto es el del informe.
    - **G4 — el material del `execution_id` no tenía vector ni oráculo.**
      Intercambiar sus dos líneas dejaba 1080/0. Ahora hay un vector congelado
      con dos SHA-256 literales y, sobre todo, una **rederivación independiente**
      desde semilla y digest con las mismas primitivas públicas que
      `src/results.rs` usa para F1–F4, más la comprobación de que el material es
      exactamente `general ‖ escenario` y de que el orden contrario da otra
      identidad.
    - **G5 — DTC-045.** Ver la entrada propia de más arriba.
    - **M6 — `run` repetía entera la admisión física.** El dispatcher
      materializaba, descartaba y `run` volvía a admitir: medido con `strace`,
      `ref.fa` se abría **tres** veces y el metadata y `ld_profile.tsv` **dos**.
      `CliValidatedConfig::F5` lleva ahora el `EffectiveF5Config` ya admitido,
      como la rama `Legacy`. Medido de nuevo: `validate` abre cada entrada
      **una** vez y `run` abre `ref.fa` **dos** —la admisión y la lectura
      productiva de bases, que no es una segunda admisión— y el metadata y el
      perfil **una**. `validate_f5_admission_after_bootstrap` desaparece por
      falta de llamador.
    - **M7a — la regla del nombre publicable estaba duplicada sin juez.** Vivía
      en el orquestador Rust y en `_PUBLISHABLE_NAME` del adaptador, y divergía:
      `validate` la aceptaba y `run` la rechazaba. Ahora vive **una sola vez**
      en la admisión del dispatcher, de modo que `validate` y `run` la aplican
      igual, y el adaptador conserva solo una guarda de **contención** —no
      escribir fuera del directorio de autoría—, que es otra regla. Hay prueba
      Rust y prueba Python.
    - **M7b — `controlled` se anunciaba pero no se desactivaba.** REQ-UI-15 dice
      «desactivarse». `render_f5_scenario_toml` rechaza ahora la combinación, de
      modo que el documento no llega a existir; Rust ya la rechaza por su cuenta
      en la subpoblación y por unidad, así que la interfaz no inventa regla.
    - **M7c — `run` F5 aceptaba y descartaba en silencio** las dos
      autorizaciones estrechas de la matriz SNV de F2. Ahora las rechaza con
      diagnóstico, en `validate` y en `run`.
    - **M7d — `f5_package_summary` podía levantar `ValueError`** por su
      `zip(..., strict=True)` ante un manifiesto con un enlace de más o de
      menos. Comprueba el recuento antes.
    - **M7e — con cero subpoblaciones el renderizador emitía un escenario que su
      propio parser rechazaba**, y la interfaz lo mostraba en verde. Ahora lo
      exige también al renderizar.
    - **Notas aplicadas.** El control positivo de la guarda de auditores
      distingue un nombre presente y cercano de uno ausente y cercano; las
      aserciones de clausura pasan de «mayor que cero» a la fórmula del fixture
      —un individuo, dos copias, aportado y no derivado—;
      `execute_validated_f5_config_with_entropy` desaparece, que era `pub` sin
      llamador; las tres expresiones regulares de los esquemas se contrastan
      ahora contra sus `pattern`, igual que `F5_INTEGER_BOUNDS`; y el comentario
      «Once ficheros» pasa a «Diez», que es lo que la aserción exige.
  - **Puertas.** `fmt`, `clippy -D warnings`, `test --workspace --all-targets`
    (1084 pasan, 0 fallan, 12 ignoradas, frente a 1076/0/12 al empezar y
    1080/0/12 antes de aplicar la revisión), `test --doc` (2/0/0),
    `build --release --target x86_64-unknown-linux-gnu` y `doc --no-deps` sin
    avisos, todas en verde. **Lado Python:** `marimo check --strict
    ui/gps_app.py` pasa de **rojo** (1 crítico `multiple-definitions` sobre
    `selected_technical_event`, preexistente) a **verde**; `pytest tests` 7
    pasan y 1 se omite, sin cambio; `pytest tests_ui` pasa de 65 a 94 pruebas y
    de 80 a 180 subpruebas, conservando **los dos fallos preexistentes**
    —`test_f2_5_vcf_ui` y
    `test_ui_technical_event_v2_round_trips_to_rust_validation`— que ya fallaban
    en `ed30562` y quedan fuera del alcance de este hito.
- ☐ **F5.6.2 — aplazada como ampliación futura (13 de septiembre de 2026).**
  Incorporar el perfil opcional `human_ld_ragsdale_approx_v1` mediante su
  recurso versionado y receta de medias, sin beta/cuantiles activos ni
  selección predeterminada.
  - Depende de: F5.6.1.
  - Procede de: DTC-089 «Paquete de perfiles» y decisión funcional preservada.
  - Cambia: **nada**. El hito no se implementa: no hay recurso incorporado, ni
    licencia o cita nuevas, ni fixtures. No existe `resources/f5/`.
  - Se comprobaría con: SHA-256 de fuente/recurso, reproducción de cada `mu`,
    ausencia de `kappa` en cálculo y selección exclusivamente explícita. Esas
    cuatro comprobaciones siguen siendo el encargo de quien lo retome.
  - **Decisión del usuario, 13 de septiembre de 2026: aplazamiento.** El perfil
    incorporado nunca fue obligatorio. El borrador retirado que lo propuso decía
    que la primera versión «**podrá** incluir» el perfil y nunca «deberá»
    (`git show 81e6589:design/borrador-generacion-automatica-bloques-0.1.md`,
    §2.1). F5 se cierra sin él: un perfil suministrado por la persona usuaria
    cubre el mismo contrato y atraviesa el mismo validador. El aplazamiento
    queda registrado como **DTF-058** en
    `registro-ampliaciones-futuras-0.1.md`, que es el registro único de
    capacidades funcionales excluidas de v1.
  - **Motivo que decide: la fuente no cubre el rango de distancias que el perfil
    necesita.** El eje x de la figura está en **centimorgans**, de 0,03 a 1 cM;
    `ld_profile.tsv` trabaja en **pares de base** y sus bins deben cubrir todas
    las distancias enteras de 1 a `Lmax-1`. Con la equivalencia humana habitual
    de 1 cM ≈ 1 Mb —una aproximación, no un dato del artículo— la curva cubre de
    unos 30 kb a 1 Mb, con dos consecuencias. Primera: convertir centimorgans a
    pares de base exige una **tasa de recombinación supuesta** que ningún
    contrato fija, que varía por región y población y que es justo lo que la
    figura promedia. Segunda: de 1 pb a unos 30 kb **no hay ningún dato**, y ese
    tramo es el grueso del rango exigido, de modo que habría que extrapolar
    varios órdenes de magnitud por debajo del primer punto medido y publicarlo
    bajo una cita ajena, que es precisamente lo que DTC-089 evita al exigir
    «bytes de digitalización» y no una curva ajustada. **No es un problema de
    licencia ni de precisión: la premisa del hito era errónea.** Quien lo retome
    con esa misma figura se encontrará el mismo muro, y lo descubriría después
    de digitalizar.
  - **Dos hallazgos de la fuente, comprobados sobre la figura del preprint CC BY
    y conservados para que nadie los «corrija» al revés.** El panel correcto es
    el **2C**, como decía el borrador: su eje es `r²_÷`, el estimador que
    DTC-089 manda usar, mientras que el **2B** es `r²_RH`, el de Rogers-Huff,
    que el propio borrador prohíbe confundir con él; un resumen de búsqueda web
    afirmó lo contrario y la imagen lo desmiente. Y las poblaciones publicadas
    son **YRI, CDX, CEU, MXL y PUR**: el borrador retirado declaraba
    `["CEU", "CHB", "MXL", "PUR", "YRI"]`, es decir **CDX, no CHB**, que son
    poblaciones distintas —Dai de Xishuangbanna frente a Han de Pekín—. Quien
    retome el hito debe usar CDX o no podrá llamar al perfil una aproximación a
    esa figura.
  - **Licencia: solo el preprint es digitalizable.** El artículo de revista
    —*Molecular Biology and Evolution* 37(3):923-932 (2020), Oxford University
    Press— es «all rights reserved». El preprint en bioRxiv
    ([10.1101/557488](https://www.biorxiv.org/content/10.1101/557488v4)), mismo
    título y mismos autores, es CC BY 4.0. Cualquier digitalización futura debe
    partir de esta segunda versión y citarla como tal.
  - **Cuatro decisiones materiales que siguen abiertas, además de los bytes.**
    Primera: el paquete son **cinco** TSV y DTC-089 fija la procedencia solo de
    las medias de LD; `block_lengths.tsv`, `snp_gaps.tsv`, `snp_spectrum.tsv` y
    `haplotype_weights.tsv` del perfil incorporado no tienen fuente normativa
    —el precedente F2.2.6 resuelve un caso análogo con literales propios
    (64:12:8) y no con una digitalización—. Segunda: `src/f5_profile.rs` exige
    que los bins de cada población terminen **exactamente** en `Lmax-1`, con
    `Lmax` el máximo de `block_lengths.tsv` del mismo paquete, así que elegir
    las longitudes de bloque y elegir hasta dónde llega la curva son una sola
    decisión. Tercera: DTC-089 exige `population_profile_id` **igual** al
    `population_id` de la subpoblación, sin remapeo, de modo que publicar
    YRI/CDX/CEU/MXL/PUR obliga a quien use el perfil a nombrar así sus
    subpoblaciones. Cuarta: `source_citation` y `source_license` son campos
    obligatorios del metadata y hay que rellenarlos con las condiciones reales
    de la versión digitalizada.
  - **Tres verificaciones del estado actual que el aplazamiento no pierde.**
    La **selección exclusivamente explícita ya se cumple estructuralmente**:
    `schemas/config-f5-scenario-v1.schema.json` exige siempre
    `profile_metadata_path` y `profile_metadata_sha256`, no hay predeterminado y
    no existe ninguna resolución tipo `builtin:` en `src/` —la sintaxis
    `profile = "builtin:human_ld_ragsdale_approx_v1"` murió con el borrador
    retirado—; el único valor con forma de predeterminado es el marcador de
    formulario `"profiles/profile.metadata.json"` de `ui/f5_editor.py`, que no
    debe repuntarse nunca a un paquete incorporado. Un perfil incorporado
    **no necesitaría vía de carga nueva**: entra por la misma frontera sellada
    que uno suministrado, con sus **seis** identidades físicas —metadata y cinco
    TSV— en `input_identities`, y `provenance_kind:"built_in"` ya está admitido
    por el lector y por el esquema. Y **`kappa` no aparece hoy en ninguna ruta
    de cálculo**: sus únicas apariciones en el árbol son prosa de este plan y de
    DTC-089, ninguna en `src/`, `tests/` ni `ui/`.
- ✅ **F5.6.3.** Ejecutar todos los gates, guardar evidencia y auditar las 46
  filas exactas de la matriz REQ→ACE→hito→prueba con hito F5, incluidas
  `REQ-VAL-54`, `REQ-MAN-70` y `REQ-SAL-02`, sin omitir filas compartidas ni
  incorporar filas ajenas al alcance.
  - Depende de: F5.6.1. F5.6.2 queda aplazada como ampliación futura y no
    bloquea el cierre.
  - Procede de: todos los REQ/ACE F5 y reglas de avance de la hoja de ruta.
  - Cambia: evidencia de cierre y estado de trazabilidad, no contrato biológico.
  - Se comprueba con: recuento mecánico exacto 46, unicidad de los 46 REQ,
    presencia expresa de `REQ-VAL-54`/`REQ-MAN-70`/`REQ-SAL-02`,
    las tres puertas de `CONTRIBUTING.md` con `--locked`,
    lectores JSON/TSV/JSONL, oráculo, vectores, paridad y búsquedas negativas de panel,
    SV, repeticiones, BAM/CRAM y recombinación bajo F5, **los dos auditores
    compuestos sobre el destino ya renombrado** —`audit_f5_package` y
    `audit_f5_identity_closure`, éxito y fallo— y una búsqueda negativa que
    demuestre que ninguna vía de publicación de F5.6 llama a
    `publish_f5_success` ni a `publish_f5_failure` fuera de los envoltorios que
    sellan.
  - **Ejecutado (13 de septiembre de 2026), pendiente de revisión.** La casilla
    sigue en ☐ hasta que la Revisión F5.6 la apruebe. El acta está en
    `evidence/f5.6/cierre.md` y las comprobaciones nuevas en
    `tests/f5_6_closure.rs`.
  - **Recuento mecánico, y por qué lo es.** 46 filas, 46 REQ distintos y los
    tres nombrados. El recuento sale de **contar la tabla** término a término
    —el campo Hito admite varios hitos separados por `/`, así que `F2.0/F5.6`
    cuenta y `F50.1` no— y se contrasta contra tres verdades que se mantienen
    por separado: el literal del contrato, la prosa de la §1 de la matriz y la
    prosa de esta misma entrada. No hay ninguna lista de 46 identificadores
    escrita a mano. `the_projection_discriminates_a_changed_matrix` demuestra
    sobre una tabla sintética que añadir, quitar o reasignar una fila mueve el
    recuento, y cuatro mutaciones sobre la matriz real lo confirman.
  - **Las 46 filas estaban todas en «Planificada».** La matriz describía la
    evidencia *prevista* y no la observada, de modo que no servía como gate de
    cierre. Las 46 pasan ahora a declarar lo que F5 cerró, con citas que una
    prueba comprueba que existen y que son pruebas. Las seis filas compartidas
    con F6, F7 y F8 conservan expresamente su parte pendiente.
  - **Hueco encontrado y cerrado: la cota de ploidía de F5.** REQ-CONF-65 exige
    rechazar una ploidía fuera de `1..=10` identificando el problema. La guarda
    existe en `src/adapters.rs` y **no tenía prueba**: la batería de F5.0.1 solo
    sustituía `ploidy = 2` por `ploidy = 1`, que es legal, y la construcción del
    efectivo se limita después a una conversión a `u8`.
    `the_f5_admission_rejects_a_ploidy_outside_one_to_ten` la cubre y muere
    contra esa misma variante defectuosa.
  - **Defecto material encontrado y corregido: la cota de propuestas del
    auditor de fallo.** `audit_f5_failure_package` rechazaba con `GPS-RES-101`
    el paquete que la propia cadena productiva publica al agotarse la búsqueda,
    porque comparaba los tres contadores contra `limits` **posición a
    posición** y a `swap_proposals` le aplicaba `max_swap_proposals`. El
    contrato dice otra cosa en dos sitios: DTC-089 l. 1196 fija
    `max_swap_proposals` como presupuesto de **un intento** —«Un intento
    termina … o al consumir `max_swap_proposals`»— y l. 1847-1848 enumera
    `swap_proposals` junto a `haplotype_attempts_used`, que es **acumulado**.
    El productor hacía lo correcto. Reproducido con el binario: `32 * 256 =
    8192` frente a un límite de `256`. La cota correcta es el producto,
    `swap_proposals <= max_swap_proposals * haplotype_attempts_used`, con
    `checked_mul` —techo contractual `10^7 * 10^4 = 10^11`— y `swap_proposals`
    **fuera** del bucle posicional, que es lo que hacía pasar por simetría una
    coincidencia de índices. El porqué queda escrito en el código citando las
    dos líneas de DTC-089. La cota es **ajustada**: el productor acumula los
    dos contadores en el mismo recorrido. El defecto era de **F5.4.2** y solo
    podía aparecer al componer el auditor sobre un agotamiento **real**, que es
    por lo que ninguna revisión por hito lo vio. Discriminación comprobada por
    los dos lados: con la cota vieja la prueba falla al aceptar, y con una cota
    `u64::MAX` falla al rechazar.
  - **Los cuatro auditores quedan compuestos en producción.** Resuelta la cota,
    `audit_f5_failure_package` ya puede componerse:
    `audit_published_f5_failure` reconstruye la `F5FailureDeclaration` desde el
    **manifiesto publicado** —no desde el informe en memoria, que haría la
    auditoría autocertificante— y corre sobre el destino ya renombrado, después
    de `audit_f5_failure_closure`. La composición se demuestra por
    **comportamiento**: un agotamiento real devuelve `GPS-DOMAIN-104` y no el
    diagnóstico de ningún auditor, y devolviendo la cota vieja el orquestador
    pasa a devolver `GPS-RES-101`, lo que solo puede ocurrir si el auditor está
    compuesto y se alcanza.
  - **Atomicidad de proceso, cerrada en WSL2/ext4.** Interrupción real por
    `SIGKILL` y concurrencia **entre procesos** sobre el `run` F5, con el solape
    comprobado y no supuesto. Tres mutaciones de `src/publication.rs` —retirar
    `RENAME_NOREPLACE`, abrir los permisos del provisional y renombrar el
    marcador `INCOMPLETE`— matan las pruebas. La primera mutación de permisos
    **no mutó** —una llamada posterior a `set_permissions` la revertía— y quedó
    registrada. **Pendiente la mitad que DTC-045 pide en Ubuntu nativo.**
  - **Deudas recontadas.** Siguen exactas: el dictamen de que **es v8 el
    incompleto** y la recomendación de `f5-resource-accounting-v9` con dos
    términos nuevos que reabriría los dos esquemas; la desviación de hasta 2x
    del presupuesto sellado, misma raíz; el docstring de F5.3 que anuncia
    `GPS-CONF-101` donde el código devuelve `GPS-INTERNAL-101`; y DTC-066. La
    condición «bloqueante en el primer `run` productivo» **ya no es teórica**:
    ese `run` existe desde F5.6.1 y compone `audit_f5_package`, de modo que la
    cadena productiva ejecuta una pieza que su contrato de recursos no modela.
    Caber con holgura —7,4 a 8,1 MiB frente a 1 GiB— no es estar declarado. La
    superposición entre el oráculo de F5.4.2 y la batería del escritor de
    F5.4.1 **no se retira**: el plan decía «puede», no hay defecto y
    `tests/f5_artifacts.rs` es evidencia propia del escritor. DTF-058 sigue
    aplazada y no bloquea.
  - **Puertas.** `fmt`, `clippy -D warnings`, `test --workspace --all-targets
    --locked --no-fail-fast` (**1098 pasan, 0 fallan, 12 ignoradas**, frente a
    1085/0/12 al empezar; las ignoradas **no cambian**),
    `test --doc --locked` (2/0/0), `build --release --locked --target
    x86_64-unknown-linux-gnu` y `doc --no-deps`, todas en verde y sin avisos.
    **Lado Python:** `marimo check --strict` exit 0; `pytest tests` 7 pasan y 1
    se omite; `pytest tests_ui` 94 pasan y **2 fallan**, `test_f2_5_vcf_ui`
    (F2.5) y `test_ui_technical_event_v2_round_trips_to_rust_validation` (F3,
    desplazada por F4.1). Los dos son preexistentes, ajenos a F5 y quedan
    declarados con su causa en el acta en vez de escondidos.

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
