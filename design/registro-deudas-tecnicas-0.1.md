# Registro de deudas técnicas 0.1

## 1. Finalidad, y en qué se distingue de los otros registros

Registro único de **deficiencias conocidas del código y de las pruebas que ya
existen**. Cada entrada describe algo que **está mal o sin cubrir hoy**, no algo
que falte por construir.

Dónde va cada cosa, para que nadie vuelva a dudarlo:

| Si es… | Va a |
|---|---|
| Una capacidad funcional excluida de la versión 1 | `registro-ampliaciones-futuras-0.1.md` |
| Una decisión de representación o implementación | `decisiones-tecnicas-0.1.md` |
| **Una deficiencia del código o de las pruebas actuales** | **este documento** |
| Un paso pendiente de una fase concreta | el plan de esa fase |

Hasta ahora las deudas vivían **dispersas dentro de los planes de fase**, lo que
hacía que una deuda transversal no tuviera sitio y que nadie supiera cuántas hay.
Este registro las reúne.

Una entrada se cierra cuando la deficiencia deja de existir, **no** cuando se
documenta mejor.

---

## 2. Registro

| ID | Deuda | Gravedad | Detectada | Estado |
|---|---|---|---|---|
| DT-001 | Dos tercios de las comprobaciones de los verificadores independientes no las ejercita ninguna prueba | **Grave** | 21/09/2026 | Abierta, muy reducida (29/09/2026) |
| DT-002 | Dos pruebas ignoradas están vaciadas, no solo apagadas | Media | 21/09/2026 | Abierta |
| DT-003 | Ocho pruebas ignoradas corresponden a contratos que ya no existen | Baja | 20/09/2026 | Abierta |
| DT-004 | La puerta de Python no corre desde el 25/08/2026 | Media | 19/09/2026 | Abierta |
| DT-005 | La marca del lector de verdad que distingue un tramo sustituido por un tracto no la ejercita ninguna prueba | Baja | 22/09/2026 | Abierta |
| DT-006 | El sorteo del inserto por rechazo es inutilizable cuando la copia más larga acerca el máximo al mínimo del perfil | Media | 24/09/2026 | Abierta, con guarda de admisión (25/09/2026) |
| DT-007 | Ningún diagnóstico público nuevo puede publicarse en un manifiesto del recorrido general | Baja | 25/09/2026 | Abierta |
| DT-008 | El orden entre el paso de linaje y la decisión de base no lo separaba ninguna prueba | Baja | 26/09/2026 | **Cerrada (26/09/2026)** |
| DT-009 | Las dos fuentes de la calibración discrepan en cuánto multiplica la PCR el tartamudeo de los homopolímeros | Media | 27/09/2026 | Abierta |
| DT-010 | El recuento de Poisson devolvía error donde el método pide repetir | Media | 27/09/2026 | **Cerrada (27/09/2026; el oráculo, el 30/09/2026)** |
| DT-011 | Dos rechazos del generador de alelos no se alcanzan desde una configuración | Baja | 27/09/2026 | Abierta |
| DT-012 | El barrido de opciones comprueba que el control existe, no que atienda a la opción | Baja | 27/09/2026 | Abierta |
| DT-013 | El generador de alelos no se puede usar sobre una referencia que construye el propio programa | Media | 29/09/2026 | **Cerrada (01/10/2026)** |
| DT-014 | El esquema de la verdad de F5 admite todavía la copia entera del haplotipo en un alelo que F5 nunca publica | Baja | 30/09/2026 | Abierta |
| DT-015 | El validador de la verdad por flujo rechaza toda verdad con alelos de repetición | Baja | 30/09/2026 | Abierta |
| DT-016 | El tracto de un microsatélite compuesto no se coloca en una referencia simulada | Media | 01/10/2026 | **Cerrada (01/10/2026)** |
| DT-017 | El juez de los recursos publicados cita el esquema efectivo v4 en sus rechazos aunque el documento sea v6 | Baja | 01/10/2026 | Abierta, corregida para la regeneración de la referencia simulada |
| DT-018 | Una copia con cambios de base aplicados fuera de orden de posición aborta la ejecución al reconstruir la verdad | **Grave** | 02/10/2026; redefinida el 03/10/2026 | **Cerrada (03/10/2026)** |
| DT-019 | Sin VCF, nadie juzga los genotipos de la verdad al publicar | **Grave** | 02/10/2026 | **Cerrada (03/10/2026)** |
| DT-020 | Sin VCF, quitar un caso indel o de repetición de la verdad se publica | Media | 03/10/2026 | **Cerrada (03/10/2026)** |
| DT-021 | El juez de la verdad no ancla los casos todo-REF, deja pasar un SNP quitado y clasifica mal una corrupción estructural | Media | 03/10/2026 | **Cerrada (04/10/2026)** |
| DT-022 | El juez de los recursos publicados resuelve la ruta de una matriz reutilizada respecto al directorio de trabajo, no al paquete | Media | 04/10/2026 | **Cerrada (04/10/2026)** |
| DT-023 | Resolver una referencia desde un manifiesto juzga la verdad entera de aquel paquete | Baja | 04/10/2026 | **Cerrada por aceptación (04/10/2026)** |

---

## DT-001. Los verificadores independientes están casi sin probar

**Gravedad: grave.** Sigue siendo la deuda más importante del proyecto, pero el
23 de septiembre de 2026 se cerró su mitad más peligrosa y se puso una puerta
que impide que vuelva a crecer sin que nadie se entere.

### Qué pasa

Los tres módulos que actúan como **jueces independientes** —los que permiten
afirmar que la verdad publicada es correcta— tenían la mayoría de sus
comprobaciones **sin ejercitar por ninguna prueba**. Estas son las dos
mediciones, hechas con el mismo método y con la misma regla de recuento:

| Verificador | Comprobaciones | Nunca se disparaban (21/09) | Comprobaciones hoy | Nunca se disparan hoy |
|---|---|---|---|---|
| Oráculo de aceptación del paquete (`src/f5_oracle.rs`) | 365 | 252 | 365 | **52** |
| Oráculo semántico del VCF (`src/vcf_semantic_validation.rs`) | 77 | 39 | 98 | **22** |
| Auditoría de identidad y publicación (`src/f5_closure.rs`) | 42 | 11 | 42 | **5** |

**Desde el 29/09/2026 la puerta vigila también `src/results.rs` y
`src/event_manifest.rs`**, que comprueban la configuración efectiva y el
manifiesto de evento. El 27/09/2026 una revisión adversarial había encontrado
ahí ciento treinta líneas de guarda **sin ejercitar como rechazo por ninguna
prueba** —sustituirlas enteras por una línea dejaba pasar la batería completa—
y la puerta no podía enterarse porque no miraba esos ficheros. Se cubrieron con
una prueba de mutación; ahora la puerta los mide. Primera medición, sin cubrir
nada a propósito:

| Módulo vigilado | Comprobaciones | Nunca se disparan |
|---|---|---|
| Configuración efectiva y resultados (`src/results.rs`) | 364 | **311** |
| Manifiesto de evento (`src/event_manifest.rs`) | 147 | **117** |

**Es decir: de 511 comprobaciones que cuenta la puerta en estos dos módulos,
428 no las ejercita ninguna prueba** (el 84 %). El techo de la línea base es
esa cifra y el objetivo de esta entrega es solo que no crezca; cubrirlas es un
trabajo aparte, pendiente de decisión. Las cifras de los tres jueces no se
movieron (365/52, 98/22, 42/5).

**El validador del manifiesto general entró el 29/09/2026, y confirmó el
pronóstico.** Era la primera de las dos familias que quedaban fuera del
denominador: otro validador del mismo fichero, con sus propias constructoras,
que nadie vigilaba sin que se hubiera decidido. Al anotar `manifest_error` y
`manifest_error_for_schema` y meterlas en la puerta, `src/results.rs` pasa de
**364 a 499 comprobaciones** y de **311 a 424 sin ejercitar**: de las 135 que
estaban en la sombra, **113 no las toca ninguna prueba** y 22 sí. La cifra no
empeora porque el código empeore, sino porque antes no se veía — el mismo
efecto que cuando se enderezaron los trece rechazos construidos a mano del
verificador del VCF y las guardas sin disparar subieron de 17 a 37.

**Y 499 y 147 siguen siendo un mínimo, no el total**, porque queda la segunda
familia, que es estructural:

- los validadores auxiliares que devuelven `Err("texto")` en vez de un
  diagnóstico (catorce funciones que devuelven `Result<_, &'static str>`,
  con setenta y seis `Err("texto")` a pelo, y además guardas con `ok_or("texto")`):
  su texto lo convierte, en el llamante, **una sola llamada** a `effective_error` o
  `repeat_section_error` así que decenas de guardas cuentan
  como un único punto, que basta que se dispare una vez para darlas por
  cubiertas.

En `src/event_manifest.rs` no hay este segundo problema (no tiene ningún
`Err("...")` a pelo), y su única construcción de diagnóstico fuera de la
constructora es `retag_event_manifest_schema`, que reetiqueta un error ya
rechazado. Una prueba léxica (`tests/guard_constructors.rs`) prohíbe construir
un `GPS-FMT-006` fuera de las constructoras de ambos módulos; **no cubre** las
dos familias anteriores.

Tres advertencias sobre la comparación, para que nadie la lea mal:

- **El verificador del VCF creció dos veces.** Tenía 58 comprobaciones el 21 de
  septiembre y 77 el 23: los microsatélites le añadieron veinte, de las que
  diecisiete ya nacían ejercitadas y tres no. Y después pasó a 98, que es lo
  que **siempre** tuvo: trece rechazos se construían a mano en vez de pasar por
  una constructora, y con ellos otros ocho quedaban sin contar. La revisión
  adversarial lo encontró. Al enderezarlos, las guardas sin disparar saltaron
  de 17 a 37 —de modo que veinte estaban muertas y nadie podía saberlo— y
  bajaron a 25 al cubrirlas.
- **El oráculo grande mide 365 y no 366.** La cifra anterior contaba también
  la llamada que hace el envoltorio `truth_reject` a `reject`, que no es un
  punto de rechazo sino un renombre. El resto de las cifras coincide
  exactamente, de modo que las dos mediciones son comparables.
- **Un rechazo fuera de la constructora es invisible para la puerta**: no se
  cuenta en el denominador y no se anota al dispararse. Por eso hay ahora una
  prueba léxica que prohíbe construir el error del verificador del VCF fuera
  de sus cinco constructoras, y una que exige que cada constructora lleve la
  anotación que la ata al barrido.

En la auditoría de identidad, **todas** las comprobaciones que quedan sin
disparar están documentadas en su propio código como inalcanzables, con su
porqué. En el verificador del VCF, tras la medición del 29/09/2026, también
todas las que quedan (22). Lo que queda del oráculo grande se detalla más abajo.

### La puerta

`scripts/check-guard-coverage.sh <directorio-nuevo>` corre la batería con el
barrido encendido y compara el resultado con `scripts/guard-coverage-baseline.tsv`,
que declara por juez **cuántas comprobaciones tiene** y **cuántas pueden quedar
sin disparar**. Las dos columnas se comprueban: sin la primera, borrar una
guarda bajaría las dos cifras a la vez y la puerta lo tomaría por una mejora.
Si el número de comprobaciones baja, la puerta falla; si sube, avisa. Las piezas son tres: la anotación
`#[track_caller]` en las constructoras de error de cada juez,
`src/guard_trace.rs`, que anota la línea de cada rechazo cuando la variable de
entorno `GPS_GUARD_TRACE_DIR` nombra un directorio, y `scripts/guard-sites.awk`,
que enumera los puntos de rechazo del texto fuente.

**Histórico: hasta el 29/09/2026 la puerta solo miraba tres módulos.**
Se anotó aquí el 26 de septiembre de 2026 porque la entrega 4 de microsatélites
añadió dos más que **no** están en la línea base ni los mide nadie: las guardas
de `validate_lineage_counts` en `src/event_manifest.rs` —que comprueban que los
tres recuentos de linaje nombran los mismos loci, que los heredados no superan
a los resbalones y que los pasos escritos no superan a las parejas— y las de
`validate_effective_clonal_stutter` en `src/results.rs` —que comprueban que la
probabilidad por ciclo publicada sale de las dos fracciones publicadas y del
número de ciclos, y que la marca de tosquedad está—.

**Medido, no supuesto.** La revisión independiente de cierre las neutralizó una
a una y ejecutó la batería completa: **1449 pasan, 0 fallan, 11 ignoradas** con
tres guardas de los recuentos de linaje y dos de la configuración efectiva
desactivadas. Ninguna prueba construye un manifiesto ni una configuración
efectiva malformados por ese lado.

**La excepción, que conviene no exagerar la desnudez.** De las guardas de los
recuentos de linaje, la que exige que no haya más pasos escritos que parejas
**sí se dispara**: se comprobó mutando el contador para que sumara dos en vez
de uno, y la publicación se rechazó a sí misma con el mensaje exacto. Y parte
de lo que las demás comprueban —la forma de las listas, los dominios de cada
campo, que no haya claves desconocidas— lo cubre además el validador externo de
Python contra los esquemas publicados, que sí corre en la batería. Lo que queda
de verdad sin cubrir son las relaciones entre campos que un esquema JSON no
puede expresar.

**Qué mide esa puerta y qué no**, dicho sin adornos: mide si una guarda **se
dispara**, no si una prueba **fallaría** en caso de que la guarda estuviera mal
escrita. Son cosas distintas y la segunda es la que importa. Una prueba que
provoca un rechazo sin exigir **cuál** rechazo satisface la puerta igual que
una que exige el motivo exacto. Por eso la puerta es un **suelo** y no un
techo: impide que la cobertura caiga sin aviso, y nada más. Lo que demuestra
que una guarda está probada sigue siendo mutarla y ver morir la prueba.

### Lo que se cubrió, y cómo se demostró

Las guardas nuevas se cubrieron con corpus de corrupciones que exigen el
**motivo exacto** del rechazo, no solo que haya rechazo, porque casi todas
comparten causa tipada. Cada una se verificó **mutando la guarda una sola vez y
comprobando que la prueba muere**, con la huella del fichero antes y después y
la copia de seguridad fuera de git. **132 mutaciones en total y ninguna
sobrevivió**: 127 sobre comprobaciones concretas —126 distintas, porque la que
ata la columna de referencia al genoma real se comprobó de dos maneras— y cinco
sobre las dos guardas léxicas nuevas, que también tienen que poder fallar.

Las peor cubiertas eran las que vigilan la verdad publicada, y son las que más
bajaron:

| Función | Puntos de rechazo | Muertos antes | Muertos ahora |
|---|---|---|---|
| `recount_f5_man70`, recuento independiente desde la verdad | 31 | 30 | 8 |
| `read_truth_copy`, validación de cada copia haplotípica | 30 | 26 | 5 |
| `recount_genotype`, validación de cada genotipo | 29 | 24 | 5 |
| `read_truth_header`, cabecera de la verdad | 17 | 13 | 3 |

De los 21 que quedan en esas cuatro funciones, **18 están documentados como
inalcanzables**, cada uno con su porqué escrito en el módulo —unos junto a la
línea y otros en la documentación de la función que los contiene, cuando una
misma razón cubre a varios—: once porque
la comprobación de claves exactas que los precede ya exigió esa clave, uno
porque el lector de genotipos acaba de resolver esa misma subpoblación, dos
porque son conversiones que no pueden fallar en 64 bits, uno porque una línea
sin salto solo puede ser la última del fichero, uno porque el primer genotipo
siempre abre locus, uno porque las dos comprobaciones de orden anteriores lo
implican, y uno que ya estaba declarado inalcanzable desde antes. Y **tres
siguen sin cubrir**: que un sitio aparezca dos veces dentro de la misma copia, que una
copia aparezca dos veces, y que una subpoblación publique dos loci en la misma
posición. Las tres exigen un documento saturado que hay que construir a mano.

Y en el verificador del VCF se cubrió **la comprobación que ata la columna de
referencia al genoma real**, que la medición anterior señalaba por su nombre.
No se alcanza corrompiendo el VCF —la comparación contra la derivación
independiente la caza antes—, de modo que se ejercita con una prueba unitaria
en el propio módulo que la llama directamente, y su documentación dice quién
caza antes y por qué esa prueba tiene que ser unitaria.

### Lo que queda

**Medición del 29/09/2026 (tras la revisión adversarial): 52 comprobaciones del oráculo grande (antes 168) y
22 del verificador del VCF (antes 25), con el mismo número de comprobaciones en
los dos, de modo que la mejora es real y no borra guardas.** Se añadieron 28
pruebas (1523 pasan, 0 fallan, 11 ignoradas). Cada prueba nueva exige la causa
tipada **y** un fragmento del detalle que solo escribe su guarda, y se verificó
**mutando la guarda sobre una copia** con huella SHA antes y después: 109
comprobaciones del oráculo (ninguna sobrevivió salvo la que se dice abajo), 20
mutaciones adicionales sobre disyunciones sueltas de condiciones compuestas y
las 3 del VCF. Al comprobarlo aparecieron dos huecos que se cerraron: una
disyunción del orden del layout que otra cubría por solapamiento, y las tres
comprobaciones de brazos de `match` que el mutador automático no sabía aplicar.

**Lo que sigue sin dispararse, y por qué (52 + 22):**

- **Oráculo grande, 52.** 51 están documentadas en su código como
  inalcanzables, cada una con su porqué: conversiones a `u32`/`usize` que no
  fallan, claves que un `exact_object` anterior ya exigió, implicaciones entre
  comprobaciones que corren antes, errores de entrada y salida que una prueba no
  puede provocar de forma determinista, y la construcción de una auditoría
  DTC-049 que ningún llamador externo puede fabricar. Una, la de
  `TruthCopies::next_copy` sobre un fichero de copias sin salto de línea final, la
  **dispara** una prueba pero está **sombreada**: el recuento MAN-70 rechaza el
  mismo defecto con el mismo texto después, y mutada sobrevive.
- **Corregido por la revisión.** La guarda «una copia aparece dos veces» se
  había declarado inalcanzable con un razonamiento falso: la clave del orden
  tiene cinco campos y la identidad del mapa cuatro, sin subpoblación; dos
  subpoblaciones que compartan identificador de individuo la alcanzan. Ahora
  tiene prueba. También se cubrieron el bloque que empieza en la coordenada cero
  y el orden estricto de copias (`>=` frente a `>`), y `validate_id_catalog` se
  parametrizó (`validate_id_catalog_text`), con lo que sus seis rechazos tienen
  prueba.
- **Redundantes, dichas en el código:** `locus.block_id != placed.block_id` (un
  genotipo de otro bloque abre grupo nuevo) y `!value.is_finite()` en el F
  efectivo (`contains` ya es falso para NaN e infinitos).
- **Aviso:** la mutación «catálogo RNG ajeno» fija el literal de la versión 2 y
  solo discrimina con fixtures sin distribución de recuentos.
- **Verificador del VCF, 22.** Las 21 que ya estaban documentadas como
  inalcanzables, más la construcción final de una secuencia FASTA, que se
  demostró inalcanzable: `ReferenceSequence::new` solo rechaza que la identidad
  difiera del texto de la cabecera, y la identidad sale de ese mismo texto. Se
  ejercitaron las otras tres que llevaban «no se ha demostrado»: la identidad de
  secuencia que el analizador de cabeceras acepta y el tipo de identidad rechaza
  (`chr#1`), y las dos de `expected_gt` con una prueba unitaria sobre un registro
  armado a mano. De estas dos, el alelo ausente es inalcanzable desde
  `expected_records`, y la muestra sin copias en el contig **sigue sin
  demostrarse inalcanzable de extremo a extremo**.
- **Disyunciones no aisladas**, dicho para que nadie las crea cubiertas: en
  la comprobación del locus recontado contra `f5-sites.tsv` (6058) no se aisló la
  del sitio, la de la secuencia ni la del censo (esta última es inalcanzable: la
  cabecera ya se cruzó con el censo declarado); en el espectro, la de que el
  mínimo esté en [0,1] la implica la clave de la auditoría, y la de que el valor
  sea finito la implica la pertenencia al intervalo.

Sigue sin medirse lo mismo que entonces: las fases 3 y 4 al completo, la
admisión del productor, el contrato de la primera versión de la verdad, y
`src/results.rs` y `src/event_manifest.rs`.

---

## DT-002. Dos pruebas ignoradas están vaciadas, no solo apagadas

**Gravedad: media.** En `src/execution.rs`, dos pruebas marcadas como ignoradas
llevan además `#[allow(unreachable_code)]` y un `return;` **a media función**.
Todo lo que viene detrás —las aserciones sobre el manifiesto fallido, el estado
del evento, las referencias descartadas— es código muerto **incluso si alguien
las reactivara**.

No son pruebas apagadas: son **pruebas vaciadas que conservan el nombre de lo
que decían comprobar**. Quien las encienda creerá haber recuperado cobertura que
no existe.

O se restauran enteras, o se borran. Conservarlas así es peor que no tenerlas.

---

## DT-003. Ocho pruebas ignoradas son de contratos que ya no existen

**Gravedad: baja.** De las once pruebas ignoradas del proyecto:

- **ocho** corresponden a contratos que la fase 4 dejó obsoletos al pasar de
  lecturas sueltas a lecturas en pareja: una en `src/publication.rs`, dos en
  `src/execution.rs` —las de DT-002—, una en `tests/f3_event_assembly.rs` y
  cuatro en `tests/f3_event_admission.rs`;
- **dos** son experimentos de rendimiento que se ejecutan a mano, en
  `tests/f3_wgs_sampling.rs` y `tests/f4_demultiplexing.rs`;
- **una** es un ayudante de subproceso, en `tests/publication_process.rs`.

Las tres últimas están bien así. Las ocho primeras **deberían borrarse**: una
prueba apagada porque su contrato desapareció no es inofensiva, porque nadie
sabe si lo es sin volver a averiguarlo — y en dos de los ocho casos la respuesta
resulta ser que **ni siquiera dirían nada si se encendieran**.

---

## DT-004. La puerta de Python no corre desde el 25 de agosto de 2026

**Gravedad: media.** `tests/test_process_adapter.py` importa `pytest`, que no
está declarado en `pyproject.toml`. Bajo la puerta canónica —`unittest`— el
módulo aborta al importar y se cuenta como error; solo corre invocándolo con la
dependencia añadida a mano.

Viene del commit `6139fe9`. Registrada también en el plan de la fase 5; se
recoge aquí porque es transversal. Quien la cierre debe decidir entre declarar
`pytest` como dependencia de desarrollo o reescribir el módulo con `unittest`.

---

## 3. Lección transversal

Estas deudas comparten una forma: **algo que parece cubierto y no lo
está**. Una guarda que existe pero no se ejercita, una prueba que conserva su
nombre pero no su cuerpo, una prueba apagada que nadie sabe si importa, una
puerta que se cita en los informes y no se ejecuta.

De ahí la regla que el proyecto ya aplica a las regresiones y que conviene
extender a la verificación entera: **no basta con que exista la comprobación;
hay que demostrar que puede fallar**.

---

## DT-005. La marca de «este tramo lo ocupa un tracto» no se ejercita

**Gravedad: baja.** Es un caso concreto de DT-001, y se registra aparte porque
nace hoy y con una causa distinta: no es una guarda olvidada, es una guarda que
**otra decisión dejó inalcanzable desde dentro**.

**Qué es.** El lector independiente de la verdad anota, por cada tramo de
referencia que una copia conserva, si lo ocupan las bases del genoma o un
**tracto de repetición realizado**. La marca sirve para que un caso de cambio de
base que cayera dentro de un tracto sustituido **no lea de ahí la base de la
referencia**, que no es la que la copia tiene.

**Por qué no se ejercita.** Desde que la decisión S-32 rechaza **antes de
simular** cualquier cambio pequeño dentro del intervalo de un locus de
repetición, este programa **no puede producir** un documento donde eso ocurra.
Ninguna prueba la dispara porque ninguna ejecución la alcanza.

**Por qué no se retira.** El lector **no produce la verdad: la verifica**, y su
trabajo es no fiarse de quien la escribió —esta versión, otra anterior u otro
productor—. Sin la marca, un documento así daría por bueno un dato falso en
silencio, que es exactamente la clase de fallo que DT-001 describe.

**Cómo se cerraría.** Construyendo a mano un documento de verdad que combine un
fragmento de tracto sustituido con un caso de cambio de base dentro de él, y
exigiendo que el lector lo rechace. Es la misma técnica que DT-001 propone para
las demás: **demostrar que la comprobación puede fallar**, en vez de confiar en
que existe.

---

## DT-006. El sorteo del inserto por rechazo se vuelve inutilizable en el extremo corto

**Gravedad: media.** No da un resultado incorrecto: da el correcto tardando
tanto que parece un cuelgue, y un fixture escrito de buena fe cae dentro sin
avisar.

### Qué pasa

El perfil técnico vigente sortea la longitud del inserto de una **normal de
media 350 y desviación 20**, truncada al intervalo `[250, 450]` **por
rechazo**: se generan pares de uniformes y se descarta el valor hasta que cae
dentro. El muestreador, además, **condiciona el máximo a la copia más larga**
del individuo, para que no salga un inserto que ninguna copia pueda alojar.

Las dos cosas juntas tienen una consecuencia que no está escrita en ningún
sitio: **cuando la copia más larga se acerca al mínimo del perfil, el intervalo
admitido queda entero en la cola de la normal y la aceptación se desploma**. El
mínimo, 250, está a **cinco desviaciones** de la media.

| Copia más larga | Intervalo admitido | Sorteos por pareja, del orden de |
|---|---|---|
| 450 pb o más | `[250, 450]` | 1 |
| 350 pb | `[250, 350]` | 2 |
| 280 pb | `[250, 280]` | 4 000 |
| 250 pb | `[250, 250]` | 13 000 000 |

### Cómo se detectó

Escribiendo el fixture saturado del tartamudeo (F4.8.4), que el plan pedía con
una copia de **exactamente 250 pb** para forzar el inserto mínimo. Una ejecución
de **dos parejas** tardó **111 segundos** medidos. Con doscientas parejas el
proceso llevaba veinticuatro minutos girando cuando se interrumpió.

### Qué se hizo y qué queda

El fixture se rehízo con una copia de 280 pb y el tracto colocado de modo que
los dos mates lo cubran entero; la propiedad que la prueba necesita se conserva
y el coste baja a milisegundos. La aritmética está escrita en
`tests/stutter_pipeline.rs`.

**Lo que no se ha hecho** es tocar el muestreador, y es una decisión del
responsable: cambiar el rechazo por un método de inversión cambiaría los
vectores congelados de F4 y todas las salidas ya publicadas. Queda registrado
como ampliación futura DTF-060.

### La guarda de admisión (25/09/2026, corregida el mismo día)

La primera versión de la guarda predecía el coste **por la forma de la
distribución**: comparaba la copia más larga admisible contra la media y la
desviación del perfil, en desviaciones típicas, y rechazaba por debajo de un
umbral fijo. Una revisión adversarial encontró que eso **mide la magnitud
equivocada**. Lo que cuesta tiempo de verdad es el número de parejas pedidas
dividido por la probabilidad de acertar en cada sorteo, no la probabilidad
sola:

- una copia de 260 pb con **una sola pareja** pedida se resolvía en 81 ms sin
  guarda, y la primera versión la rechazaba igualmente, prohibiendo una
  ejecución que acababa en menos de un décimo de segundo;
- un panel de 60 tramos de 271 pb a 30x —muchas parejas, cada una barata—
  tardaba 78 segundos, cinco veces más que el mismo panel con tramos de 350
  pb, y la primera versión lo admitía sin más porque solo miraba la forma, no
  el volumen pedido.

Un umbral por forma, además, es en la práctica un límite de tamaño de
secuencia codificado, y la especificación exige que ese tamaño no tenga
límites codificados y que un escenario grande y legítimo solo pueda advertir,
nunca rechazar.

**El criterio vigente cuenta intentos de verdad**, en vez de predecirlos.
`f4_truncated_normal`, en `src/randomness.rs`, cuenta sus propios intentos de
rechazo y aborta ese único sorteo con
`RandomSamplingError::TruncatedNormalAttemptsExhausted` si agota
`F4_TRUNCATED_NORMAL_MAX_ATTEMPTS` (un millón) sin acertar en el intervalo
admitido; no consume ni un número de la cinta de azar de más al abortar, y
para cualquier sorteo que sí termina, el resultado y el consumo son idénticos
a los de antes de que el contador existiera —comprobado byte a byte contra un
binario del commit anterior, compilado aparte, sobre un escenario congelado
que completa sin agotar la cota—. En `src/wgs_sampling.rs`,
`WgsPairedEndSampler::insert_sampling_error` traduce ese fallo concreto —y
solo ese— al mismo diagnóstico público de siempre,
`WgsSamplingError::ImpracticalInsertSampling`, con el individuo, su copia más
larga y el mínimo del perfil.

Medido: una copia de 271 pb necesita del orden de **23 000** intentos por
sorteo, el peor caso legítimo conocido; 280 pb, unos 4 000; 250 pb, el propio
mínimo del perfil, unos **13 000 000**. Un millón deja unas cuarenta veces de
margen sobre el peor caso legítimo y sigue cazando el desbocamiento de trece
millones muy por debajo de su coste real. No hay cota agregada sobre el total
de intentos de una biblioteca entera: cada sorteo se acota por separado, así
que un panel grande y proporcionado a lo que pidió la persona usuaria —muchas
parejas, cada una barata— no se rechaza por volumen; solo se rechaza el
sorteo concreto que de verdad se desboca.

Las pruebas viven en `tests/f4_insert_sampling_guard.rs` (nivel de biblioteca,
con dos de la revisión adversarial que comprueban que la guarda mira la copia
haplotípica realizada y no la secuencia de referencia) y, para el contador en
sí, en `src/randomness.rs` (dos pruebas unitarias con una fuente de azar
controlada: una agota la cota exactamente y sin un intento de más, la otra
demuestra que el último intento permitido todavía puede acertar).

**Lo que queda:** DTF-060, el método de muestreo sin rechazo, sigue aplazado
y no es de esta guarda. Y el diagnóstico reutiliza `GPS-DOMAIN-003` en vez de
un código propio, por la limitación que registra DT-007.

---

## DT-008. El orden entre las dos partes del tartamudeo no lo separaba ninguna prueba

**Cerrada el 26 de septiembre de 2026.** Se abrió ese mismo día con un
razonamiento **equivocado**, y la revisión adversarial lo refutó ejecutando. Se
deja escrita entera porque el error es instructivo.

### Qué se afirmó, y por qué era falso

Se dijo que ningún fixture podía separar el orden entre el paso neto heredado
del linaje y la decisión de tartamudeo de base, «por una razón estructural: las
dos sortean su paso del mismo reparto de pesos y lo filtran con la misma guarda
de viabilidad, de modo que intercambiarlas deja el mismo par de pasos aplicados
sobre el mismo par de tiradas».

**La guarda no es fija: depende del estado que deja el primer paso.** La parte
de base calcula sus candidatos con la tirada que **encuentra**, es decir con la
que el paso de linaje le deja. Si el primero deja la tirada corta, el segundo
pierde alternativas. Las dos partes **dejan de conmutar en cuanto la guarda
muerde**, y ese régimen es el de los tractos cortos con los pesos cargados a la
contracción.

El error de la medición que llevó a la afirmación fue elegir un fixture con los
dos pasos forzados **al mismo valor**. Ahí, y solo ahí, los dos órdenes
escriben los mismos bytes.

### Cómo se cierra

Con la prueba `the_lineage_step_goes_before_the_base_decision`, en
`tests/clonal_lineage.rs`: dinucleótido de **cuatro** unidades, un ciclo con
probabilidad de resbalón 1, tasa de base 1 y pesos `−3 → 0,9` y `+1 → 0,1`.

- Orden vigente: el linaje contrae de cuatro a una unidad, y la parte de base
  ya solo puede expandir. **Ninguna lectura muestra una unidad.**
- Orden invertido: la base contrae primero, el paso neto del linaje ya no cabe
  y se descarta. **237 de 300 lecturas muestran una unidad.**

La aserción es «cero lecturas con una unidad», y está verificada por mutación.

### Lo que queda anotado, por si vuelve a aparecer

Hay un segundo camino por el que el orden se observa y que la primera versión
también descartó mal: `slip_probability` recibe la tirada **posterior** al paso
de linaje, de modo que con
`technical_event.stutter.run_length_coefficient` distinto de cero —un campo que
el autor del escenario puede declarar— la propia probabilidad de la parte de
base depende del orden. Hoy ninguna prueba lo ejercita, porque el coeficiente
arranca en cero y no está calibrado (pendiente P-1).

---

## DT-007. Ningún diagnóstico público nuevo puede publicarse en un manifiesto del recorrido general

**Gravedad: baja.** No es un defecto que afecte a ningún resultado publicado
hoy: es un techo que hay que conocer antes de intentar cerrar el siguiente.

### Qué pasa

Todo diagnóstico público lleva un código estable de un catálogo versionado
(DTC-044), publicado en `catalogs/diagnostics-v4.tsv`. Pero ese fichero no es
el que de verdad rige cuando un fallo controlado publica un **manifiesto del
recorrido general** —el general y el de cada evento técnico, los documentos
de F1 a F4—: esos tres puntos de anclaje están **fijados en el código** a la
versión histórica **`"3"`** del catálogo, la que existía antes de que F5
añadiera sus propios códigos:

- `src/results.rs`, dos veces, dentro de `validate_diagnostics(diagnostics,
  "3", path)`: una en `parse_general_manifest_bytes_unretagged`, la lectura
  vigente y alcanzable del manifiesto general; la otra en
  `parse_general_manifest_v1_object`, marcada `#[allow(dead_code)]` y sin
  ninguna llamada hoy —el contrato v1 archivado que describe DTC-036—, así
  que hoy solo la primera se ejercita de verdad;
- `src/event_manifest.rs`, la constante `EVENT_MANIFEST_DIAGNOSTIC_CATALOG_VERSION
  = "3"` — valida el manifiesto de cada evento técnico.

Esto es coherente con DTC-036: el manifiesto general conserva su versión de
esquema `"2"` desde F2, cerrada, y por eso solo puede citar los códigos que
existían cuando se cerró. F5 introdujo su propia familia de códigos
(`GPS-DOMAIN-101` en adelante, `GPS-CONF-101`, `GPS-RES-101`) sin tocar esa
frontera, porque sus fallos se publican en un documento propio —el de fallo
F5, con su propio esquema— y no en el manifiesto general. Hay una prueba que
demuestra justo esto por el lado negativo: `tests/f3_event_manifest.rs`, el
caso `f4-event-f5-diagnostic`, construye un manifiesto de evento con el código
F5 `GPS-CONF-101` y comprueba que **tanto el validador Rust interno como el
externo lo rechazan**.

**La consecuencia práctica:** ningún código nuevo puede añadirse al catálogo y
usarse en un fallo del recorrido general —F1 a F4, y ahora F6— sin que la
publicación del propio fallo controlado se rompa contra su propio manifiesto.
Se descubrió así, no por inspección: la primera versión de la guarda de
DT-006 iba a abrir `GPS-DOMAIN-107`, y al escribir la prueba de extremo a
extremo, el intento de publicar el fallo controlado se rechazaba a sí mismo
con `GPS-FMT-006` («el diagnóstico GPS-DOMAIN-107 no coincide con
diagnostics-v3.tsv»). La guarda de DT-006 reutiliza `GPS-DOMAIN-003`, un
código ya presente en la versión histórica, para esquivarlo.

### Qué se hizo y qué queda

Nada, todavía: es un techo conocido, no un arreglo. Cualquier defecto nuevo
del recorrido general que necesite nombrar algo que el catálogo histórico no
distingue tendrá que reutilizar un código existente y llevar el matiz en
`detail`, exactamente como hizo DT-006, o abrir una versión nueva del esquema
del manifiesto general y del de evento técnico —una decisión de arquitectura,
no de implementación, porque cambia un contrato que DTC-036 declaró cerrado—.

---

## DT-009. Las dos fuentes discrepan en cuánto multiplica la PCR el tartamudeo de los homopolímeros

**Detectada el 27 de septiembre de 2026**, al ir a las fuentes originales para
decidir la longitud de tracto de referencia de la dependencia con la longitud
del tracto. No se ha cambiado nada por ella: queda registrada para decidirse
aparte, porque toca un modelo ya commiteado.

### Las dos cifras, y quién dice cada una

| Motivo | Gymrek (2016), tabla 1 | Fungtammasan y col. (2015), texto del artículo |
|---|---|---|
| Mononucleótido | **6,30×** (0,17 frente a 0,027) | **~2,5×** («la diferencia aumenta hasta ~11 pb, tras lo cual se estabiliza en una diferencia de ~2,5 veces») |
| Dinucleótido | **3,45×** (0,038 frente a 0,011) | **~3,6×** («aumenta hasta ~16 pb, estabilizándose en ~3,6 veces») |

En dinucleótidos concuerdan bien. En mononucleótidos difieren por un factor
**2,5**, y no es una lectura a ojo de una figura: las dos cifras están escritas
con todas las letras en el texto de sus artículos.

### Por qué importa

Porque el tartamudeo clonal —DTC-092, el modelo que reparte entre los ciclos de
la PCR el incremento que la amplificación añade— **despeja su probabilidad por
ciclo exactamente de esa diferencia**, con
`r_pcr = 1 − (1−r_con)/(1−r_sin)`. Si el multiplicador real de los
homopolímeros fuera 2,5 y no 6,30, la probabilidad por ciclo que el simulador
usa hoy para los homopolímeros estaría sobrestimada.

### Qué no se ha hecho, y por qué

**No se ha tocado nada.** El proyecto adopta las cifras de Gymrek porque son las
que la calibración declara como fuente de la tabla por longitud de motivo, y
cambiarlas sería rehacer un modelo cerrado a partir de una frase de otro
artículo, medida sobre otra muestra y con otro programa de genotipado. La
elección corresponde al responsable, no a quien encontró la discrepancia.

### Qué haría falta para cerrarla

Decidir cuál de las dos fuentes manda en el multiplicador de los homopolímeros,
o declarar que la discrepancia entra dentro de la tosquedad ya declarada de las
tasas por ciclo. Si se decidiera por Fungtammasan, habría que rehacer la
columna sin PCR de los mononucleótidos o la de con PCR, y con ella las seis
probabilidades por ciclo publicadas.

---

## DT-010. El recuento de Poisson devolvía error donde el método pide repetir

**Encontrada** el 27 de septiembre de 2026 por la revisión de cierre del
generador de alelos. **Corregida el mismo día.** Se registra igualmente, porque
tuvo alcance más allá de quien la destapó y porque conviene que conste cuánto
tiempo vivió.

### Qué pasaba

`poisson-count-v1` usa, para medias de diez o más, el método de rechazo con
transformación. Cuando el candidato que propone cae **fuera del dominio**, el
método canónico **lo rechaza y repite**. Esta implementación **devolvía un
error interno**. El resultado era que un parámetro perfectamente legal abortaba
la ejecución con `NumericFailure`, y lo hacía **según dónde cayera la secuencia
de azar**: el mismo valor aceptaba sobre un tracto y fallaba sobre otro.

La revisión lo reprodujo con una media de veinte pasos sobre un tracto de tres
unidades.

### Desde cuándo, y a quién afecta

El cuerpo era **idéntico byte a byte al del commit anterior** al generador de
alelos, de modo que la deficiencia es anterior a ese trabajo. Lo que el
generador añadió fue **exponerla a un parámetro que escribe el usuario** y cuyo
esquema no ponía techo; los demás consumidores del recuento —la composición
genotípica y el recuento estocástico por biblioteca del recorrido técnico—
alimentan la media desde valores acotados por otras guardas, y por eso nadie se
había tropezado con ella.

### Qué se hizo

Un candidato **negativo** se rechaza y se repite, que es lo que el método pide.
Un candidato **no numérico** sigue siendo error: no se puede repetir lo que no
tiene valor, y reintentarlo sería un bucle sin salida. Y la media de pasos del
generador tiene ahora **techo declarado** en el esquema y en el código, mil
pasos por linaje, holgadamente por encima de cualquier tracto que este
simulador admita.

**Ninguna ejecución que antes terminara cambia de resultado.** La rama nueva
solo se toma donde la anterior devolvía error, de modo que un recorrido que
acabó consumió las mismas palabras y obtuvo lo mismo. Las huellas congeladas del
recorrido técnico no se mueven, y se comprobó.

### La corrección dejó atrás al juez (30 de septiembre de 2026)

El arreglo del 27 de septiembre tocó solo el productor. **El oráculo de
aceptación de F5**, que reimplementa el mismo muestreador para rederivar desde
la semilla lo que se publicó, **conservó la regla vieja**: ante un candidato
negativo devolvía «sin valor», y de ahí salía el rechazo «la lambda declarada no
admite el muestreador de recuentos». La revisión adversarial del 30 de
septiembre lo descubrió y lo reprodujo de dos maneras. Con la misma cinta
aleatoria, el productor daba un recuento donde el oráculo daba nada en 8 de 256
cintas con media 12 y en 7 de 256 con media 20 (cintas derivadas de las
semillas por otro procedimiento que el de la prueba de regresión, que da 7 de
256 en las dos medias). Y de extremo a extremo, 3 de 61
semillas de un escenario F5 correcto, con el recuento de haplotipos de media
12, acababan en `GPS-INTERNAL-101` y un paquete en estado `failed`: el productor
publicaba bien y el juez rechazaba lo que acababa de publicarse.

**Corregido el mismo 30 de septiembre.** El oráculo aplica ahora la misma
regla que el productor: un candidato no numérico o por encima de `2^53` se
rechaza, y uno negativo se descarta y se vuelve a tirar. Sigue siendo una
reimplementación propia, sin llamar al productor. La prueba de regresión da la
misma cinta, derivada de 256 semillas, a los dos muestreadores con medias 12 y
20 y exige que coincidan en el recuento y en las palabras consumidas; y **cuenta
las cintas con candidato negativo** (7 de 256 en cada media, y son las mismas
siete semillas en las dos: no son dos comprobaciones independientes) para que no
sea un vector donde el efecto pasa inadvertido. Se comprobó que falla, en la
semilla 69 de la media 12, con la regla vieja.

Una segunda prueba cubre el límite superior, con la media máxima admitida, `2^52`,
que es donde un candidato por encima de `2^53` se alcanza de verdad (del orden de
`2^-30` por vuelta): las semillas 0 a 7 coinciden en productor y oráculo, y la
semilla 1 392 399 649 hace que el productor dé `NumericFailure` y el oráculo
nada. Mata dos mutantes que antes sobrevivían a la batería entera: el oráculo
sin cota superior y el oráculo con la cota bajada a `2^52`.

La lección es la habitual de este proyecto: una corrección en un productor que
tiene un juez independiente se hace **en los dos**, y se busca el juez antes de
cerrar la deuda.

### Qué queda abierto

Nada de esta deficiencia. Queda anotado, eso sí, que **el techo de la media es
una guarda de la frontera y no del método**: `poisson-count-v1` sigue admitiendo
medias hasta `2^52` para sus demás consumidores.

Una cautela sobre el alcance de ese límite: con medias pequeñas el candidato por
encima de `2^53` es inalcanzable en la práctica (del orden de `2^-50` por vuelta
para medias entre 270 y mil), y solo se vuelve probable con medias enormes. El
perfil F5 admite hasta `2^52`, así que el rechazo es alcanzable desde un
escenario, aunque con una media que ningún uso razonable escribe.

---

## DT-011. Dos rechazos del generador de alelos no se alcanzan desde una configuración

**Anotada** el 27 de septiembre de 2026 por la revisión de cierre. No es un
defecto: es una laguna de cobertura **declarada**, con la forma que este
proyecto ya usó para otras.

Los dos rechazos son:

- **El tramo de la referencia sin ni una copia completa del motivo.** Con
  referencia suministrada lo caza antes la guarda de admisión de los loci, que
  exige una tirada mínima —tres unidades por omisión—; con referencia simulada
  esa guarda no corre, pero el tracto lo pone el propio simulador. *Corregido
  el 01/10/2026:* esto último era falso hasta ese día —el tracto no se ponía y
  el rechazo **sí** se alcanzaba desde una configuración, que es DT-013—. Desde
  la corrección es cierto para los loci de un motivo, los únicos que el
  generador admite: el tramo declarado mide al menos un motivo y se rellena con
  él. Y desde el mismo día la guarda de admisión **también corre** con
  referencia simulada, sobre el tracto colocado (corrección de S-06), de modo
  que este rechazo vuelve a estar detrás de ella en los dos tipos de
  referencia.
- **Las frecuencias generadas degeneradas.** Son proporciones de linajes con
  denominador positivo, de modo que su suma solo podría no ser finita si el
  recuento de eventos desbordara, y eso lo caza antes su propio rechazo.

**Los dos se conservan y ninguno se borra**: son el contrato del generador y no
pueden depender de que otra guarda siga existiendo donde está hoy. Los dos se
ejercitan por las pruebas del módulo, que llaman al generador directamente; lo
que **no** existe es un escenario de configuración que los alcance, y así se
dice en la documentación de cada uno en vez de dejar creer que discriminan algo
que nadie puede provocar desde fuera.

### Qué haría falta para cerrarla

Decidir si la guarda de admisión de los loci debe aplicarse también con
referencia simulada. Es una decisión del responsable —la decisión S-06 la
excluye a propósito— y no se toma aquí.

---

## DT-012. El barrido de opciones comprueba que el control existe, no que atienda a la opción

**Anotada** el 27 de septiembre de 2026 por la revisión de cierre, que la
demostró: añadió una opción al esquema del documento general, la apuntó en la
tabla a un **control ajeno que ya existía**, y el barrido pasó.

### Qué mide de verdad

El barrido recorre las opciones declarables que el trabajo de microsatélites y
amplificación añadió al esquema, y exige que **la tabla de controles cubra el
esquema entero** y que el nombre de cada control aparezca en el cuaderno o en su
panel. Eso impide lo que existe para impedir: **añadir una opción y olvidarse de
darle control**, porque una opción sin fila en la tabla rompe la prueba.

Lo que **no** demuestra es que el control apuntado **atienda** a esa opción. Es
un suelo contra el olvido, no una demostración de que la interfaz hace lo que su
fila dice.

### Qué haría falta para cerrarla

Ejecutar el cuaderno y observar el formulario que produce, opción a opción, en
vez de buscar nombres en el texto de los ficheros. Es trabajo de otra clase —una
prueba de la interfaz viva, no de su fuente— y no estaba en el alcance de esta
entrega.

---

## DT-013. El generador de alelos no se puede usar sobre una referencia que construye el propio programa

**Cerrada el 1 de octubre de 2026.** Se conserva la descripción, corregida: la
original contenía una afirmación falsa que escondía la causa.

**Gravedad: media**, y era más grave de lo que decía. La entrada original
afirmaba que «no hay ningún resultado incorrecto: el programa rechaza y lo
dice». No lo era para los alelos declarados a mano, que publicaban resultados
incorrectos sin avisar, ni siempre para el generador, que a veces aceptaba
(las dos cosas, más abajo).

### Qué pasaba

Un microsatélite se puede declarar sobre los dos tipos de referencia, la que
aporta el usuario y la que construye el propio programa. Con el **generador de
alelos por deslizamiento**, la entrada decía que la referencia construida por
el programa **se rechazaba siempre**:

```
detail="locus STR-1: el tramo de la referencia no contiene ni una copia
completa del motivo, de modo que no hay alelo ancestral del que partir"
```

El generador **lee el alelo de partida de las bases de la referencia**, y una
referencia construida por el programa llevaba bases sorteadas al azar en la
ventana declarada.

**Y tampoco se rechazaba siempre**, como decía la entrada: el rechazo solo
saltaba si en esas bases al azar no aparecía ni una copia del motivo. Con un
motivo corto casi siempre aparece alguna. El 01/10/2026, con `CA` sobre treinta
bases al azar que contenían algún `CA` suelto, el programa de antes de la
corrección **aceptó** el escenario y publicó como alelo ancestral el tramo al
azar, con pureza 0,27, sin avisar de nada.

### Cómo se reprodujo

El 29 de septiembre de 2026, al preparar la primera simulación de uso real —un
cromosoma de 2000 bases con `(AGAT)₁₂`, seis individuos diploides—. El mismo
escenario, cambiando **solo** el bloque del generador por dos alelos
declarados, validaba con código de salida 0; con el generador, se rechazaba con
el mensaje de arriba. Se resolvió entonces construyendo el FASTA a mano con el
tracto dentro.

### La afirmación falsa, corregida

La entrada decía: «cuando los alelos se declaran, el programa sí sabe
materializar el tracto sobre una referencia construida por él». **No lo
sabía.** Lo que hacía era **sobrescribir el intervalo en las copias** de los
individuos con el alelo sorteado; la **referencia** seguía llevando bases al
azar en ese tramo. Las consecuencias, que nadie había visto porque no había
rechazo que las delatara: el FASTA publicado no tenía el microsatélite, **todas
las copias** recibían un alelo alternativo —ningún alelo declarado coincidía
con las bases al azar— y el alelo de la referencia no lo podía tener nadie.

### La causa de fondo

La preparación de la referencia (`prepare_reference`, en `src/execution.rs`)
llamaba al generador de referencia **sin regiones**. El generador sabía colocar
regiones declaradas antes del relleno desde el 23 de septiembre de 2026, pero
esa capacidad no tenía ningún llamador fuera de sus propias pruebas. Lo
encontró la revisión adversarial del 30 de septiembre de 2026.

### Cómo se cerró

Por decisión del responsable del 1 de octubre de 2026: en un locus de **un solo
motivo**, la preparación coloca en la referencia el motivo repetido desde su
primera base sobre el tramo declarado, de `start` a `end` inclusive, **antes**
del relleno al azar; si el tramo no es múltiplo del motivo, la última copia
queda parcial. No hace falta campo nuevo: el tramo ya dice cuánto mide el
tracto. El juez que regenera la referencia simulada desde la configuración
efectiva coloca los mismos tractos con su propio cálculo, de modo que una
ejecución que dejara de pasarlos se rechaza al publicar.

Con eso el generador de alelos parte del tracto colocado, y la combinación
funciona para **todo locus que el generador admite**: el generador solo admite
loci de un motivo, con cualquier tipo de referencia. Lo prueban
`tests/simulated_reference_repeat_tract.rs` y, para la reproducibilidad, una
huella congelada de una referencia simulada **sin** microsatélites tomada con
el commit `7dd1fbe`, anterior a la corrección.

### Lo que no cerraba

El tracto de un locus **compuesto** sobre una referencia simulada. Se anotó
como DT-016 y se cerró el mismo día.

---

## DT-014. El esquema de la verdad de F5 admite todavía la copia entera del haplotipo

**Gravedad: baja.** Es una holgura del esquema publicado, no un resultado
incorrecto: ningún fichero real la aprovecha y el único lector la rechaza.

### Qué pasa

Hasta el 30 de septiembre de 2026, cada alelo de inserción o de deleción de la
verdad estructurada llevaba en `final_haplotype_bases` **la secuencia entera de
la copia portadora**. Eso hacía crecer `truth.jsonl` como el número de *indels*
por la longitud del genoma, y ese día se quitó del recorrido general: las
versiones `"1"` y `"3"` pasaron a la `"4"` y la `"5"`, sin el campo.

La verdad de F5, la de los fundadores sintéticos, es otra versión, la `"2"`, y
**se dejó como estaba**. Su esquema, `truth-v2.schema.json`, fija que todo caso
genotípico es un SNP, pero su definición de alelo sigue enumerando también los
alelos de inserción y de deleción, con el campo dentro. Un documento que metiera
un alelo de inserción con la copia entera dentro de un caso SNP **pasaría ese
esquema**.

### Por qué es inofensiva

- **Nadie la escribe.** El productor de la verdad de F5
  (`src/f5_artifacts.rs`) solo emite casos SNP con alelos de referencia y SNP.
- **El único lector la rechaza.** El oráculo de F5 (`src/f5_oracle.rs`)
  rechaza un caso que no sea SNP («F5 solo publica casos SNP») y todo alelo que
  no sea de referencia o SNP («alelo ... ajeno a F5»). Lo demuestra, fallando
  por su motivo, la prueba del oráculo que corrompe cada campo de un genotipo
  publicado.

### Por qué no se toca ahora

Cambiar el esquema obliga a cambiar la versión, y la versión de la verdad de F5
está fijada en dos contratos más: el manifiesto de F5
(`manifest-general-v5.schema.json`) la exige como constante `"2"`, y DTC-089,
la decisión que fija el F5 sintético reproducible, la escribe literalmente y
prohíbe ampliar después sus esquemas. Quitar una rama muerta costaría una
versión nueva del manifiesto de F5, con su lector y su oráculo, sin ahorrar un
solo byte en ningún fichero.

### Qué haría falta para cerrarla

Aprovechar la próxima vez que haya que cambiar la verdad o el manifiesto de F5
por otro motivo para quitar de la definición de alelo de su esquema las ramas
de inserción y de deleción.

---

## DT-015. El validador de la verdad por flujo rechaza toda verdad con alelos de repetición

**Gravedad: baja.** Hoy no afecta a ninguna ejecución: esa vía solo la usan
las pruebas. Pero es una función pública que dice validar cualquier verdad y no
lo hace.

### Qué pasa

El validador de la verdad tiene dos formas de leer `truth.jsonl`: una que carga
el documento entero y otra **por flujo**, que lee copia a copia y las libera
(`validate_truth_jsonl_streaming`, en `src/truth_validation.rs`). La segunda
comprueba la cabecera **como si fuera un documento sin registros**, antes de
haber leído ninguno. Y la regla de versiones dice que un documento sin alelos
de repetición tiene que declarar la versión general. Así que toda verdad que
declare la versión con repeticiones —la `"5"`, y antes la `"3"`— se rechaza
por versión no admitida, aunque sea correcta.

### Cómo se comprobó

El 30 de septiembre de 2026, sobre una verdad con alelos de repetición
publicada por el propio programa: la lectura completa la acepta y la lectura
por flujo la rechaza con `Contract(UnsupportedVersion)`. Se reprodujo igual
sobre el código anterior a la versión 5, con la verdad `"3"`: la deuda no la
introdujo el cambio de versiones.

### Por qué no se arregla ahora

Ninguna ejecución la usa, y arreglarla obliga a decidir cómo se comprueba por
flujo una regla que depende de registros que todavía no se han leído: aplazar
la comprobación de versión al final de la lectura, o exigirla por registro. Es
un cambio del validador que merece su propio encargo.

### Qué haría falta para cerrarla

Que la lectura por flujo compruebe la cabecera sin la regla de promoción,
acumule si ha visto algún alelo de repetición y aplique la regla al terminar,
con una prueba que lea por flujo una verdad `"5"` publicada.

---

## DT-016. El tracto de un microsatélite compuesto no se coloca en una referencia simulada

**Cerrada el 1 de octubre de 2026, el mismo día en que se anotó.**

**Gravedad: media.** Era el mismo defecto que DT-013 describía, en la parte que
su primer cierre no cubría: un resultado incorrecto que no avisaba.

### Qué pasaba

Con referencia simulada, el tracto de un locus de **varios motivos**
(`motifs`, como `(TCTA)n(TCTG)n`) seguía siendo **relleno al azar** en el FASTA
publicado. Ningún alelo declarado coincidía con esas bases: todas las copias
recibían un alelo alternativo y el de la referencia no lo tenía nadie.

### Por qué no bastaba el tramo

Con un solo motivo el tramo basta para construir el tracto. Con varios, no: hay
que saber cuántas unidades lleva **cada bloque** en la referencia.

### Cómo se cerró

Por decisión del responsable del 1 de octubre de 2026, que eligió implementar
el campo entero y no un rechazo provisional:

- **`reference_units`**, la lista de unidades de cada bloque en la referencia,
  en el orden de los motivos. Obligatoria en un compuesto sobre referencia
  simulada; se rechaza en un locus de un solo motivo y con referencia
  suministrada. Cada valor es un entero de al menos una unidad, como el número
  de unidades por bloque de un alelo compuesto.
- La suma de unidades por longitud de cada motivo tiene que ser exactamente la
  longitud del tramo, o la configuración se rechaza antes de simular nombrando
  el locus, la suma y el tramo.
- La preparación de la referencia coloca cada bloque con su motivo repetido, en
  orden, antes del relleno al azar.
- La configuración efectiva v6 publica el campo en el locus, y el juez que
  regenera la referencia simulada **lee de ahí** la geometría de cada tracto
  —secuencia, tramo, motivos y unidades por bloque— y lo coloca con su propio
  cálculo. Un productor que coloque mal los bloques, que ignore el campo o que
  lo publique distinto de como lo usó, se rechaza al publicar.
- El esquema del documento general v4 y el de la efectiva v6 lo declaran, y el
  editor de loci de la interfaz lo ofrece y lo conserva en la ida y la vuelta.

Cabe en las versiones abiertas sin versión nueva por la decisión S-26 del plan
de repeticiones. Lo prueban `tests/simulated_reference_repeat_tract.rs` y
`tests_ui/test_repeat_loci_ui.py`.

### Lo que el esquema no puede comprobar

Que la lista tenga tantos valores como motivos y que la suma cuadre con el
tramo: comparar propiedades hermanas entre sí no es expresable en JSON Schema.
Lo dice la descripción del campo y lo comprueban el lector de Rust y sus
pruebas.

---

## DT-017. El juez de los recursos publicados cita el esquema efectivo v4 aunque el documento sea v6

**Gravedad: baja.** No acepta nada que deba rechazar: lo que falla es el
**nombre del esquema** que el diagnóstico manda a mirar.

### Qué pasa

La función que construye los rechazos de la configuración efectiva en
`src/results.rs` (`effective_error`) escribe fijo `effective-config-v4.schema.json`.
El lector cerrado de la cabecera del documento ya corrige el nombre según la
versión declarada, pero el juez que relee los recursos publicados
(`validate_published_resources`) tiene cientos de rechazos que pasan por esa
función, y con un documento v6 todos mandan a mirar el esquema v4.

### Qué se hizo y qué queda

Lo encontró la revisión adversarial del 1 de octubre de 2026 en el rechazo de
«FASTA físico … no coincide base a base». Ese rechazo, y todos los de la
regeneración de la referencia simulada, se reetiquetan ya con el esquema de la
versión declarada; lo comprueba
`the_reference_judge_reads_the_blocks_from_the_effective_config`. El resto de
rechazos del juez sigue citando v4.

### Qué haría falta para cerrarla

Que el juez de recursos publicados reetiquete sus rechazos de la configuración
efectiva con el esquema de la versión declarada, sin tocar los del manifiesto,
que comparten código de diagnóstico y no deben reetiquetarse con un esquema
efectivo. No es una línea: hay que separar los dos tipos de rechazo antes de
reetiquetar.

---

## DT-018. Una copia con cambios de base aplicados fuera de orden de posición aborta la ejecución

**Redefinida el 3 de octubre de 2026.** Se registró como «un cambio de base
directo y otro de la matriz en la misma copia», y estaba mal acotada: no hace
falta la matriz. La revisión adversarial de ese día encontró que **dos cambios
de base directos** en la misma copia, declarados en orden descendente de
posición —el de la base 180 y después el de la 150, en la fase 0 de `ind-1`—,
abortan igual, con la versión base y con la corregida; los mismos dos en orden
ascendente publican. La condición es que los cambios de base de una copia
entren en su historia en un orden distinto del de su posición, vengan de la
matriz o de eventos directos. Los dos ejemplos: un directo seguido de uno de la
matriz en una posición más adelante, y dos directos en orden descendente.

**Gravedad: grave.** No publica nada falso, pero **impide simular** casi
cualquier escenario que combine un cambio de base declarado a mano con la
matriz de variantes de población, o varios cambios declarados a mano sin
ordenar: la ejecución termina con un error interno en vez de con un paquete.

### Qué pasa

En el primer ejemplo, una copia recibe un cambio de base directo —un evento
`snp` de `[[direct_events]]`— y, **más adelante en la misma copia**, un cambio
de base de la matriz de población; en el segundo, dos directos con la
posición en orden descendente. La ejecución aborta con `GPS-INTERNAL-001`, en la
operación de reconstruir la verdad de forma independiente, porque la historia
de la copia «no se reproduce» (`HistoryDoesNotReplay`).

La causa está en el reproductor de historias del validador de la verdad
(`src/truth_validation.rs`). Al reproducir los eventos de una copia, los
cambios de base se añaden a su lista **en el orden de la historia**
(`truth_validation.rs`, en `replay_event`, hacia la línea 773 antes de la corrección), y en la historia el de la matriz va antes que el
directo —o el de la base 180 antes que el de la 150—; después, la lista reconstruida se compara **ordenada tal cual** con
las piezas vigentes que publica la verdad (`truth_validation.rs`, en `reconstruct_copy`, hacia la línea 1149), que van
por posición. Las dos listas tienen los mismos cambios en distinto orden, y la
comparación las da por distintas.

### Cómo se reproduce

La segunda revisión adversarial del 2 de octubre de 2026 lo encontró con un
escenario mínimo: referencia de 360 bases, dos individuos diploides, la matriz
con dos cambios de base en las posiciones 21 y 200, y un cambio de base directo
en la posición 180 de la fase 1 de `ind-1`. Aborta con la versión actual y con
la anterior a la fusión de variantes solapadas del VCF, de modo que **no lo
introdujo esa fusión**: es previo.

### Cómo se cerró

**Cerrada el 3 de octubre de 2026.** El diagnóstico se confirmó ejecutando: la
verdad publica los cambios de base de cada tramo ordenados por posición (el
productor los inserta en su sitio, y la lectura estructural rechaza piezas
fuera de orden), mientras que el juez los añadía al final de la lista en el
orden del historial. Las dos listas tenían los mismos cambios, con las mismas
bases y posiciones, en distinto orden; no había ningún otro estado distinto.

La corrección está solo en el juez (`src/truth_validation.rs`): al reproducir
un cambio de base, lo inserta **en su posición** dentro de la lista del tramo,
igual que el productor. La comparación final no se tocó y sigue siendo exacta:
un cambio de más, de menos, con otra base o en otra posición se sigue
rechazando; dos cambios en la misma posición los rechaza la lectura estructural
(`ModificationsNotOrdered`) antes de llegar al reproductor; la guarda de
posición repetida del propio reproductor es defensiva y ninguna prueba la
ejercita.

Pruebas: `tests/snp_order_replay.rs`, con los dos escenarios del revisor
(`tests/fixtures/snp-order/`: dos directos en orden descendente, y un directo
seguido de uno de la matriz más adelante) y una verdad escrita a mano con
historial descendente, y otra con tres cambios en orden no monótono (4, 1, 6),
que distingue insertar en su sitio de insertar al principio o junto al último. Las pruebas de rechazo fallan si se relaja la
comparación (ignorar la base, ignorar las posiciones o aceptar longitudes
distintas). Una exploración de 2000 escenarios al azar dio 0 abortos al
reconstruir la verdad, frente a 106 con la versión anterior.

Nota: el comprobador de la exploración señalaba en 39 escenarios un caso sin
registro «propio» con un registro que toca su locus (`AUSENTE-CON-ALT`). La
revisión adversarial lo identificó como un falso positivo del comprobador: son
casos de efecto nulo adyacentes a otro registro, y el VCF es exacto (ver
REQ-VER-25).

## DT-019. Sin VCF, nadie juzga los genotipos de la verdad al publicar

**Estado: cerrada el 3 de octubre de 2026.**

**Gravedad: grave.** No es un defecto que publique algo falso por sí mismo,
sino la falta de la red que lo pararía, y falta **en la configuración por
omisión**: el VCF solo se escribe si se pide `[vcf_output]`, y sin él los casos
genotípicos de `truth.jsonl` —la fuente canónica de la verdad— salen sin que
ningún juez los haya contrastado. Ya ha dejado pasar un defecto real.

### Qué pasa

Al publicar, el productor de la verdad solo contrasta las **copias**: compara
la huella SHA-256 de cada copia que materializa con la de la reconstrucción
independiente (`verify_copies_against_independent_reconstruction`, en
`src/truth_projection.rs`). Los **casos genotípicos** —locus, alelos,
portadores y dosis— no los mira nadie. El juez semántico de la verdad
(`validate_and_reconstruct_truth` y su variante para el VCF, en
`src/truth_validation.rs`), que sí los contrasta con las copias reconstruidas,
solo se ejecuta dentro de la publicación del VCF y en las pruebas.

### Cómo se detectó

Al corregir, el 2 de octubre de 2026, el locus de los casos indel que salía en
coordenadas de la copia y no de la referencia. Hasta ese día, una ejecución sin
VCF publicaba ese locus falso sin error; con VCF abortaba, pero solo porque el
agrupador no encontraba el caso, no porque nadie juzgara la verdad. Tras la
corrección, dos mutantes de verificación lo confirmaron sobre una copia del
árbol: con el productor devolviendo el locus viejo (M1) y con el juez anclando
las inserciones a la derecha (M6), **las ejecuciones sin VCF publicaron igual**;
solo una prueba que compara lo publicado con una tabla escrita a mano
(`tests/indel_case_reference_locus.rs`) lo vio.

### Cómo se cerró

Al publicar, el productor ejecuta ahora el juez semántico de la verdad sobre el
documento ya proyectado, **con o sin VCF** y antes de escribirlo. Si rechaza,
la ejecución falla con un defecto interno (`GPS-INTERNAL-001`) que lleva el
motivo del juez.

- `src/truth_validation.rs`: nueva `validate_truth_with_copy_inspection`, **la
  única implementación del juez**: `validate_truth_for_vcf` es ahora esa misma
  función con una inspección vacía, de modo que no pueden divergir. Entrega a
  quien la llama la huella de cada copia reconstruida y suelta las bases
  **antes** de la inspección. Se retiró `reconstructed_copy_sha256`, que dejó de
  usarse.
- `src/truth_projection.rs`: el contraste de huellas
  (`verify_copies_against_independent_reconstruction`) pasó a llamarse
  `verify_publication_against_simulation` y hace **los dos contrastes con una
  sola reconstrucción por copia**: la huella de cada copia y, a continuación,
  los casos genotípicos. Hay una variante de error nueva,
  `GenotypeJudgement`. Y un gancho solo de pruebas (`test_tamper`) que altera
  el documento justo antes de verificarlo.
- `src/execution.rs`: esa variante se traduce a defecto interno, como ya se
  hacía con las de las copias.

Con VCF, el juez ya se ejecutaba varias veces dentro de la publicación del
VCF (proyección y agrupación lo repiten a propósito, para que cada frontera
pública no acepte un documento sin reconstruir las copias). Eso no se tocó. La
ejecución nueva, en la publicación de la verdad, es **una más** también con
VCF: no se pudo dejar una sola porque esas fronteras reciben solo el documento
y no tienen cómo saber que otro ya lo juzgó. A cambio, la comparación de
huellas y el juez comparten reconstrucción, de modo que la verdad no reconstruye
cada copia dos veces.

### Coste medido

Binario release, genoma simulado de 40 Mb, 4 individuos diploides, 1000 loci
de matriz y 30 inserciones directas, sin VCF y sin lecturas; tres ejecuciones
de cada binario, alternadas. Tiempo: antes 48,6 / 43,7 / 45,2 s, después
44,7 / 51,9 / 44,1 s (medianas 45,2 y 44,7 s: dentro del ruido). Las verdades
son idénticas byte a byte. Con 20 000 loci de matriz el escenario tarda más de
diez minutos; no se midió con ambos binarios.

**Memoria.** La primera versión de la corrección tenía una regresión que esta
medición de 40 Mb no veía: con un cromosoma de 100 Mb y un haploide, el máximo
subía de 204 MB a 300 MB, con VCF y sin él, porque la inspección materializaba
la copia del productor mientras las bases reconstruidas seguían vivas. Se
corrigió calculando la huella y soltando las bases antes de la inspección:
medido con `/usr/bin/time -v` sobre ese escenario, 204,5 MB antes, 300,6 MB con
la primera versión y 204,5 MB con la final, con tiempos de 11,4 / 11,4 / 11,6 s.

### Pruebas

En `src/execution.rs`, de extremo a extremo con el gancho de alteración:
`shifted_case_locus_without_vcf_is_rejected_as_internal_defect` (locus
desplazado sin VCF), `changed_case_dosage_without_vcf_is_rejected_as_internal_defect`
(dosis y portadores intercambiados de forma coherente: la estructura del
documento sigue siendo válida, así que solo el juez lo ve) y
`shifted_case_locus_with_vcf_is_rejected_before_any_vcf` (con VCF). Cada una
exige que el rechazo sea el del juez (`GenotypeJudgement`). Se verificaron
mutando sobre una copia del árbol: quitar el juez, ejecutarlo ignorando su
error, ejecutarlo condicionado a algo externo (el equivalente de «solo con
VCF» en este diseño, donde el juez es incondicional) y quitar la llamada a la
verificación en `project_truth` matan las tres pruebas.

Una exploración de 4000 escenarios al azar sin VCF (dos generadores de 2000) y
600 con VCF dio 0 escenarios que antes publicaran y ahora aborten, y verdad y
VCF idénticos byte a byte frente a la base `8d9da4d`.

### Lo que esta corrección no cubre

Sigue habiendo cosas que el juez no ve y que no son de esta deuda: que **falte**
un caso entero (DT-020) y tres limitaciones anteriores del juez (DT-021).

## DT-020. Sin VCF, quitar un caso indel o de repetición de la verdad se publica

**Estado: cerrada el 3 de octubre de 2026**, junto con DT-021; reabierta y cerrada de nuevo el 4 de octubre tras dos revisiones adversariales.

**Gravedad: media.** Aparece al revisar DT-019. El juez semántico contrasta cada
caso **que está** con las copias, pero no exige que estén **todos**.

### Qué pasa

Si un productor defectuoso omite un caso indel o de repetición (y baja el
recuento de casos de la cabecera para que el contrato estructural no proteste),
una ejecución sin VCF publica la verdad sin el caso. Con VCF, el agrupador lo
rechaza con `AlternateWithoutCase` (`src/vcf_grouping.rs`), porque encuentra un
alelo alternativo en las copias sin caso que lo explique.

### Qué haría falta para cerrarla

No basta con exigir un caso por evento: en un evento insertado y borrado después
(*inserted-then-deleted*) el primero no tiene caso, a propósito. Cerrarla exige
una comprobación de **completitud basada en la secuencia**, como la del
agrupador: toda diferencia de una copia con la referencia debe quedar explicada
por algún caso, salvo las ventanas en que la copia vuelve a ser la referencia.

### Cómo se cerró

Ver «Cómo se cerraron DT-020 y DT-021», al final de DT-021: las dos se
cerraron con la misma comprobación del conjunto de casos.

## DT-021. Tres limitaciones del juez de la verdad anteriores a DT-019

**Estado: cerrada el 3 de octubre de 2026**, junto con DT-020; reabierta y cerrada de nuevo el 4 de octubre tras dos revisiones adversariales.

**Gravedad: media.** Ya estaban antes de DT-019 y valen también con VCF.

- Un caso todo-REF (dosis `[2,0]`) no está anclado a ninguna posición: se puede
  desplazar, fabricar o duplicar con otro identificador de caso y se publica.
- Quitar un caso SNP de dosis `[0,4]` se publica, incluso con VCF, si el caso de
  otro individuo cubre ese mismo alelo alternativo.
- Una corrupción estructural del productor (un portador de otro individuo, un
  identificador fuera de léxico) sale como `GPS-FMT-006`, error de formato, y no
  como defecto interno.

### Cómo se cerraron DT-020 y DT-021

**El contrato.** Qué casos tiene que publicar la verdad estaba escrito a
trozos: DTC-061 fija que la matriz conserva también los loci todo-REF;
`schemas/truth-v4.md`, que un caso de efecto nulo sigue en el fichero porque el
evento ocurrió, y que un locus de repetición sin variantes no publica nada; y
`schemas/truth-v1.md`, que una inserción solo es válida si conserva algún
fragmento. El resto —un caso por cambio de base vigente, uno por deleción del
historial, uno por inserción con bases supervivientes, uno por tracto de
repetición realizado— solo lo fijaba el productor (`src/truth_projection.rs`).
No había discrepancia entre los dos; se ha escrito entero en la sección nueva
«El conjunto de casos» de `schemas/truth-v4.md`, y DTC-061 lo resume.

**El juez** (`src/truth_validation.rs`) juzga ahora, individuo a individuo, el
**conjunto** de casos, en las dos direcciones y con su propia réplica de las
copias, sin código del productor ni del oráculo del VCF:

- **del cambio al caso** (DT-020 y la segunda limitación de DT-021): cada
  cambio de base vigente, cada fragmento de bases insertadas, cada tracto de
  repetición y cada deleción del historial tienen su caso; las bases
  insertadas tienen que proceder de una inserción del historial, y el estado
  heredado de la copia —sus segmentos, antes del historial— tiene que
  conservar la referencia entera, de modo que toda base que falte la quitó
  una deleción del historial (en la primera versión se pedía que el hueco
  cayera en el locus de alguna deleción, y eso dejaba esconder una pérdida
  heredada; ver «Segunda revisión»); y, cuando se conoce la matriz, cada fila
  tiene el caso de cada individuo de su población;
- **del caso a su fuente** (primera limitación de DT-021): cada caso es una de
  esas fuentes, y dos casos con la misma —posición SNP, evento o intervalo de
  repetición— son un duplicado, aunque lleven `case_id` distintos.

Se trabaja sobre los fragmentos que el juez ya replicaba y no sobre las bases,
de modo que no cambia lo que se retiene en memoria por copia (DT-019): las
bases se siguen soltando antes de la inspección.

No hacía falta la excepción de «las ventanas en que la copia vuelve a ser la
referencia» que preveía la entrada: el juez compara **fragmentos**, no
secuencia alineada, y en la inserción borrada entera la copia no conserva
ningún fragmento de ella, así que no hay cambio que explicar. Lo que sí queda
es la deleción, cuyo caso se exige por estar en el historial.

**`truth.jsonl` no lleva las filas de la matriz**, solo su huella; las filas
van en su propio artefacto del paquete (`population_snvs.tsv` y sus
metadatos, ligados por SHA-256 en el manifiesto y en `group_digests`), y la
población de cada individuo está en la configuración efectiva. Con ellas
(`MatrixCaseLoci`) el juez exige los casos todo-REF y ata a su fila un caso
SNP que no es cambio de ninguna copia, en dos momentos que corren siempre,
con VCF o sin él: al proyectar la verdad
(`verify_publication_against_simulation`, en `src/truth_projection.rs`), con
la matriz en memoria, y —desde la segunda revisión— en el juez de los
recursos publicados (`validate_physical_truth`, en `src/results.rs`), sobre
los bytes del paquete. Los lectores de solo `truth.jsonl` (el lector completo,
`validate_and_reconstruct_truth`, que usa en las pruebas el oráculo del VCF
`validate_published_vcf_semantics`, y el lector por flujo, que además rechaza
toda verdad v5 por DT-015) comprueban el conjunto **sin** la matriz: un caso
SNP sin cambio en ninguna copia lo aceptan sin poder comprobarlo. Las
fronteras internas del VCF (`validate_truth_for_vcf`) no juzgan el conjunto:
solo reciben el documento y corren después de la publicación de la verdad.

*Corregido en la segunda revisión*: la primera versión de esta entrada decía
que «las filas de la matriz no se publican», y era falso; por esa premisa el
juez de los recursos publicados se quedó sin juzgar el conjunto.

**La clasificación** (tercera limitación): en `src/execution.rs`, un rechazo
del contrato estructural sobre el documento recién proyectado, un `case_id`
duplicado o un alelo de repetición fuera de su catálogo son ahora defecto
interno (`GPS-INTERNAL-001`, operación `truth_projection`). El resto de
errores de la proyección sigue saliendo como `GPS-FMT-006`, y la coincidencia
es ahora exhaustiva, de modo que una variante nueva obliga a decidir. La
lectura de una verdad suministrada no pasa por ese camino.

### Pruebas

En `src/execution.rs`, de extremo a extremo con el gancho de alteración del
productor, **cada una con VCF y sin él**, y todas exigen el error concreto del
juez del conjunto (`MissingCase` o `GenotypeMismatch` dentro de
`GenotypeJudgement`), no un rechazo cualquiera:

- quitar el caso de una inserción, de una deleción, de una repetición y de un
  SNP `[0,4]` cubierto por otro individuo;
- desplazar, fabricar y duplicar con otro `case_id` un caso todo-REF (el
  escenario añade una fila de matriz de frecuencia 0,01; una fila de
  frecuencia 0 la retira la preparación);
- quitar el caso de la deleción que borra entera una inserción
  (`inserted-then-deleted`), que ninguna diferencia de secuencia delata;
- **dos escenarios saturados**: llevar una inserción, o una deleción, del
  historial a los segmentos y quitar su caso. La copia reconstruida es la
  misma y el historial ya no tiene el evento, así que solo la comprobación de
  bases insertadas sin evento, o la del estado heredado, lo ve;
- un portador de otro individuo, que sale como defecto interno;
- controles: el escenario sin alterar y `inserted-then-deleted` publican, con
  VCF y sin él.

En `tests/truth_validation.rs`, los dos lectores (completo y por flujo) sobre
verdades escritas a mano: una deleción sin caso, un cambio de base sin caso,
un caso duplicado en el primer individuo y en el último, un caso *indel*
inventado y, como control, dos individuos completos. Dos pruebas antiguas
que construían a mano una verdad con una deleción en el historial y sin su
caso se completaron con él.

Las once pruebas de extremo a extremo de rechazo fallan contra el código de
`9c1fb5f`, con VCF y sin él. Sin VCF, diez publican y la del portador de otro
individuo sale como `GPS-FMT-006`. Con VCF, cinco publican (el SNP `[0,4]`,
los tres ataques al caso todo-REF y la deleción de efecto nulo), cinco
abortan por el agrupador del VCF (`AlternateWithoutCase`) y no por el juez
—la prueba lo distingue— y la del portador sale como `GPS-FMT-006`.

**Mutantes**, repetidos sobre el árbol final de la tercera revisión, en una
copia con su propio directorio de compilación y comprobando la huella SHA-256
del fichero antes, mutado y restaurado: treinta y ocho. **Mueren treinta y
siete.** Quitan o aflojan, uno a uno: la detección de duplicados; la
exigencia de caso para un cambio de base (dos formas), un tracto, unas bases
insertadas o una deleción del historial; el rechazo de bases insertadas sin
evento; la exigencia de las filas de la matriz; el rechazo de un caso sin
fuente, en general y solo para las repeticiones; la matriz o el conjunto
entero al proyectar; el conjunto en el lector completo; cada uno de los tres
momentos en que el lector por flujo juzga a un individuo, y su rango de
copias; la clasificación como defecto interno; el hueco inicial, el final y
la deleción de segmento del estado heredado; los alelos exactos de un caso
SNP, con y sin la fila, sin ordenar los observados o tomándolos de otras
secuencias; el orden de los alelos de repetición y sus tractos de otras
secuencias; la población de las copias, quitada o mirada solo en la
primera; el identificador de caso, quitado y con la posición en base 0; el
juez de la verdad en los recursos publicados y su matriz; la copia que falta,
la que sobra y la ploidía en ese juez; y uno común a las dos construcciones
de las filas: la posición sin pasar a base 0.

**Sobrevive uno, y se sabe por qué**: no añadir al catálogo de un caso SNP de
una fila los observados que la fila no tiene. Esa rama no se alcanza desde
una configuración (`ExplicitSnpAlleleOutsideMatrix`), y así lo dice su
comentario. En rondas anteriores sobrevivieron otros, y se corrigieron: una
mutación del juez de los recursos publicados que conservaba el `?` de la
llamada —no mutaba nada—; el rango del lector por flujo, que la prueba de dos
individuos no veía porque el segundo era un clon del primero; y, en la
tercera ronda, los tractos de otras secuencias, que necesitaron una prueba
con dos loci de repetición en las mismas coordenadas de dos secuencias.

Lo que **no** se ejercita, dicho también en su comentario: el rechazo de un
caso *indel* sin evento estructural, que el contrato estructural rechaza
antes. *Corregido en la tercera revisión*: aquí se decía lo mismo del
rechazo «el individuo no tiene población conocida», y era falso para el juez
de los recursos publicados, donde un paquete con un individuo renombrado
llega a él; ahora lo ejercita una prueba.

### Segunda revisión adversarial (3 de octubre de 2026)

La revisión del cierre confirmó que la publicación funcionaba —ningún aborto
nuevo, salidas idénticas, ninguno de sus 48 ataques en publicación pasaba—,
pero encontró el cierre incompleto. Lo que se corrigió:

- **La premisa de que la matriz no se publica era falsa.** La matriz generada
  o suministrada se publica en el paquete (`population_snvs.tsv` y sus
  metadatos, `publish_population_snvs` en `src/execution.rs`), y la
  reutilizada la lee de su ruta el juez de los recursos publicados. Ese juez
  (`validate_published_resources`) solo miraba la cabecera de la verdad, y los
  lectores de un paquete ya escrito aceptaban un caso todo-REF inventado,
  desplazado u omitido, un portador cambiado con su dosis y un caso *indel*
  quitado. Ahora ese juez, en `validate_physical_truth`
  (`src/results.rs`), lee el `truth.jsonl` del paquete y le pasa el juez
  completo con la matriz publicada y la población de cada individuo según la
  configuración efectiva (`MatrixCaseLoci::from_matrix`, la misma
  construcción que usa la proyección). Corre al publicar, desde
  `src/publication.rs`, siempre que haya poblaciones iniciales; sin ellas la
  verdad no lleva copias ni casos. Sus rechazos salen como el resto de los de
  ese juez, `GPS-FMT-006`. Es una comprobación más del juez que vigila la
  puerta de cobertura: `src/results.rs` pasa de 538 a 539 en
  `scripts/guard-coverage-baseline.tsv`, con el techo de sin disparar en 425.
- **Un caso todo-REF estaba atado a su posición, pero no a sus alelos.**
  `MatrixCaseLoci` guarda ahora los alternativos de cada fila, y el juez exige
  el catálogo exacto de cada caso SNP en el orden del productor: la
  referencia, los alternativos de la fila en orden de base y después los
  observados que la fila no tiene; sin fila, solo los observados, en orden de
  base. Cambiar el ALT de un caso todo-REF o añadir uno que nadie lleva se
  rechaza. La rama de «observados que la fila no tiene» no se alcanza desde
  una configuración: un cambio de base directo en una fila con otro
  alternativo se rechaza antes de simular (`ExplicitSnpAlleleOutsideMatrix`,
  comprobado ejecutándolo).
- **Una deleción heredada podía esconderse dentro de otra del historial.** La
  comprobación de huecos aceptaba un hueco si caía en el locus de alguna
  deleción del historial, y ese locus va de la primera a la última base que
  la deleción elimina, aunque entre medias haya un hueco anterior. Se
  sustituyó por una más simple y exacta: el estado heredado de la copia —sus
  segmentos, antes del historial— tiene que conservar la referencia entera,
  sin deleciones de segmento ni tramos que ningún segmento cubra
  (`inherited_loss`). Con eso, toda base que falte la quitó una deleción del
  historial, y esa deleción tiene su caso. El productor no escribe nunca un
  cambio estructural en los segmentos.
- **Población e identificador de caso.** Con la población conocida, el juez
  exige que cada copia declare `provenance.origin =
  "initial_population:<población>"` de su individuo y recalcula el `case_id`
  con su propia implementación —unas cuarenta líneas, sin código del
  productor—: intercambiar entre individuos dos casos todo-REF que solo se
  distinguen por el identificador ya no pasa. El formato del `case_id` y el
  del origen, que no estaban escritos, quedan fijados en
  `schemas/truth-v4.md`.
- **Mutantes que sobrevivían a la batería entera**: sin el hueco final, sin
  el hueco inicial (ahora en `inherited_loss`) y con un caso de repetición
  inventado. Tienen prueba propia. Los de `DuplicateCaseId` y
  `RepeatAlleleOutsideItsCase` devueltos a error de formato son
  inalcanzables, y así se dice ahora en el comentario de la clasificación.
- **Guardas sin ejercitar, dichas en su comentario** y no solo aquí: «un
  caso indel nombra su evento estructural». (Se decía también de «el
  individuo no tiene población conocida»; era falso, ver la tercera
  revisión.)
- **El lector por flujo** filtraba todas las copias para cada individuo, con
  coste cuadrático; ahora toma las del individuo por rango. Recuérdese que
  rechaza toda verdad v5 (DT-015): su juicio del conjunto solo vale para la
  v4.

Pruebas nuevas, además de las de la primera ronda:

- en `src/execution.rs`, con VCF y sin él: el ALT de un caso todo-REF
  cambiado; un ALT de más en el caso de un cambio de base directo y en el de
  una fila; la deleción heredada escondida dentro de otra del historial; las
  copias de un individuo con otra población; dos casos todo-REF
  intercambiados entre individuos; y un caso de repetición inventado. El
  escenario gana un cambio de base directo fuera de la matriz;
- en `tests/truth_validation.rs`, los dos lectores con un segmento recortado
  al principio y al final;
- en `tests/published_resource_validation.rs`, el juez de los recursos
  publicados sobre un paquete real alterado: caso todo-REF inventado,
  desplazado y omitido, portador cambiado con su dosis y caso *indel*
  quitado. Exige el motivo exacto del juez de la verdad.

Los ataques de la revisión se repitieron: los 48 de la primera tanda más la
deleción escondida y el caso de repetición inventado, cada uno con VCF y sin
él (100 ejecuciones), se rechazan todos; solo publica el escenario sin
alterar. Los 26 paquetes alterados de la revisión los rechaza el juez de los
recursos publicados con el motivo del juez de la verdad, y acepta el paquete
sin alterar.

### Tercera revisión adversarial (4 de octubre de 2026)

Sin ningún rechazo falso en más de 9000 escenarios, la revisión encontró dos
defectos graves y algunos huecos. Lo que se corrigió:

- **El juez de los recursos publicados no comprobaba qué copias hay.** El
  juez de la verdad recorre los individuos y las copias que el documento
  trae; al proyectar, el contraste con lo simulado exige que estén todas
  (`CopySetMismatch`), pero en el paquete no hay estado simulado, y se
  aceptaba un paquete sin un individuo entero o sin una copia de un
  triploide con sus portadores rehechos. `validate_physical_truth` exige
  ahora que las identidades de copia sean exactamente cada individuo de la
  configuración efectiva por cada secuencia de la referencia por cada copia
  de la ploidía de su población, con un rechazo propio para la que falta y
  otro para la que sobra. Son dos comprobaciones más que vigila la puerta de
  cobertura: `src/results.rs` pasa de 539 a 541, con el techo en 425.
- **La memoria subía un 41 %, no un 27 %.** La medida anterior se hizo con
  un binario de las 21:53, anterior a los últimos cambios de las 22:12; no
  era el del árbol final. La causa real: la verdad en memoria seguía viva
  mientras el juez de los recursos publicados releía la del paquete. Se
  suelta ahora (`drop(truth)` en `src/execution.rs`, junto a la liberación de
  la referencia). La afirmación de que «solo lo arreglaría juzgar por flujo»
  era también falsa. Las cifras medidas con el binario final están en «Coste».
- **Una guarda que se daba por inalcanzable se alcanza**: «el individuo no
  tiene población conocida», con un paquete en que un individuo se renombra.
  Su comentario lo dice ahora y la prueba del juez de los recursos publicados
  lo ejercita.
- **Tres mutantes sobrevivían a la batería entera**: ordenar los observados
  de un caso SNP sin fila, filtrar los observados por secuencia y comprobar el
  origen de todas las copias y no solo de la primera. Pruebas nuevas: dos
  copias de `ind-2` con cambios directos distintos en la misma posición, la
  mayor en la primera copia; dos secuencias con un cambio en el mismo
  desplazamiento en el paquete del juez de los recursos publicados; y el
  origen cambiado solo en la última copia.
- **El orden de los alelos de un caso de repetición** queda fijado en
  `schemas/truth-v4.md` —la referencia y después cada tracto distinto en el
  orden en que aparece de la primera copia a la última— y el juez lo exige
  cuando conoce la población. Antes, permutar los alelos y rehacer el
  `case_id` se aceptaba: era equivalente, pero el identificador dejaba de ser
  canónico.

Pruebas nuevas: en `tests/published_resource_validation.rs`, el paquete gana
una segunda secuencia y dos cambios de base directos, y las alteraciones
pasan de cinco a nueve —un individuo quitado, una copia quitada con sus
portadores rehechos, una copia de más y un individuo renombrado—, y una
prueba nueva publica dos loci de repetición en las mismas coordenadas de dos
secuencias; en
`src/execution.rs`, el origen cambiado solo en la última copia y los alelos
de repetición permutados, con VCF y sin él.

Dos deficiencias anteriores que encontró la revisión se registran aparte,
como DT-022 y DT-023.

### Exploración y coste

**Todo lo que sigue se midió con el binario del árbol final**, SHA-256
`25047a1db5718d9af1f9fb8edf3d9b8199c6b8c0da69f950b8ac2e70faa4049d`, frente al
de `9c1fb5f`, `1af436eaa4d6addb4d96a84fa42df915052a48fdc4adeea2092656754bae8bb2`.
Las cifras de las dos primeras rondas se hicieron con binarios intermedios y
se sustituyen; las de memoria de la segunda eran además falsas (ver la
tercera revisión).

**Exploración.** Cinco generadores: los dos de la revisión de DT-018 (SNP
directos densos y desordenados con *indels* y la matriz, ploidía 1 a 4), el de
la primera ronda (repeticiones declaradas, una o dos poblaciones, filas de
frecuencia 0,01, eventos alrededor del tracto), el del primer revisor de esta
deuda (inserciones anidadas, deleciones que borran bases insertadas o cubren
filas, SNP directos en filas, dos secuencias, hasta tres poblaciones, VCF de
texto, BGZF o ninguno) y el del segundo (matriz generada, referencia simulada
con varios cromosomas, tracto insertado en una referencia suministrada,
fundadores por clave, identificadores con puntos, guiones y subrayados,
alelos de repetición por deslizamiento). En total 6900 ejecuciones (2000,
2000, 1500, 700 y 700), de las que 5098 publican con la base. **Ninguna que
publicaba con la base aborta ahora, y ninguna publica distinto**: verdad y VCF
idénticos byte a byte. Las catorce rutas de publicación de la revisión de
DT-019 —referencia con tracto insertado y lecturas, loci con *indels* y
lecturas, referencia reducida, referencia existente, los tres casos de efecto
nulo; con VCF y sin él— publican con todos sus ficheros idénticos.

**Coste**, `/usr/bin/time -v`, ejecuciones alternadas de cada binario:

| Escenario | Tiempo base | Tiempo nuevo | Memoria máxima base | nueva |
|---|---|---|---|---|
| 100 Mb, un haploide, 10 filas, sin VCF (mediana de 3) | 12,8 s | 13,6 s | 204,0 MB | 204,7 MB |
| 2 Mb, 60 diploides, 1000 filas, 40 *indels* (60 040 casos), sin VCF (mediana de 3) | 31,2 s | 34,7 s | 307,7 MB | 313,9 MB |
| 2 Mb, 120 diploides, 2500 filas, 40 *indels* (300 040 casos), sin VCF (2 de cada) | 310,4 y 310,9 s | 320,6 y 324,1 s | 1303,8 MB | 1273,3 y 1273,6 MB |

El tiempo sube entre un 3 y un 11 %: el juez de los recursos publicados lee y
juzga la verdad entera una segunda vez. La memoria máxima sube un 2 % en el
segundo escenario y **baja** un 2 % en el tercero, porque la verdad en memoria
se suelta antes de esa segunda lectura. Las salidas de todas las parejas son
idénticas byte a byte.

## DT-022. El juez de los recursos publicados resuelve la ruta de una matriz reutilizada respecto al directorio de trabajo

**Gravedad: media. Cerrada el 04/10/2026.** La encontró la tercera revisión
adversarial de DT-020 y DT-021, el 4 de octubre de 2026; es anterior a ellas.

### Qué pasaba

Cuando la matriz de SNV se reutiliza (`population_snvs.mode = "reuse"`), el
juez de los recursos publicados reconstruía la configuración desde el TOML
copiado en el paquete, y la ruta de los metadatos de la matriz
(`metadata_path`) se resolvía respecto a la carpeta del TOML **original**
(`resolve_link` en `src/application.rs`, con la ruta de
`configuration_input.original_path` del manifiesto). `original_path` se copia
tal como se pidió, así que con un TOML relativo la ruta acababa resuelta
contra el directorio de trabajo del proceso: **el mismo paquete se aceptaba o
se rechazaba según desde dónde se lanzara el juez**. Reproducido ejecutando el
juez desde `/` y desde el directorio de la ejecución. No afectaba a la
publicación, que corre con el directorio de la ejecución, sino a quien
juzgase un paquete ya publicado: `resolve_manifest_reference` y las pruebas.
Ninguna orden de la CLI llamaba a esas funciones.

La matriz reutilizada no viaja en el paquete **a propósito**: el esquema del
manifiesto fuerza `metadata_artifact_id` nulo en modo `reuse` y el manifiesto
solo registra `metadata_sha256` y la auditoría de normalización; los
metadatos llevan dentro la huella del TSV. Copiarla o registrar una ruta
absoluta habría exigido cambiar el formato publicado, y `original_path` es
literal por una condición documentada de identidad binaria de los manifiestos.

### Cómo se cerró

Sin cambiar ningún formato publicado: **quien juzga dice dónde está la
matriz**. `validate_published_resources` y `resolve_manifest_reference`
(`src/results.rs`) reciben `reused_snv_metadata: Option<&Path>`, la ruta de
`population_snvs.metadata.json` con el TSV al lado. Es **obligatorio si y solo
si el paquete es de `mode = "reuse"`**: sin él, el juez rechaza diciendo qué
falta; sin reutilización se ignora. El juez lee la matriz de esa ubicación y
nunca de la ruta del TOML copiado, y la acepta solo si coinciden las huellas:
`metadata_sha256` del manifiesto con los bytes de los metadatos y la del TSV
dentro de los metadatos. `ExecutionWorkspace` guarda la ruta que la ejecución
acaba de usar (`set_reused_snv_metadata`, llamado desde `src/execution.rs`) y
la pasa al juez de publicación. Todas las pruebas que llamaban al juez pasan
`None` salvo la de F2.2, que pasa la ruta de la matriz reutilizada.

Pruebas nuevas, en `tests/published_reused_matrix_location.rs`: el juez
acepta un paquete con TOML y `metadata_path` relativos cuando la matriz se ha
movido a un sitio que ninguna ruta del TOML alcanza, y sin la ubicación lo
rechaza (control negativo). **Lo que discrimina es mover la matriz, no cambiar
de directorio de trabajo**: una versión anterior lanzaba el juez en un
subproceso y esa parte no aportaba nada, así que se quitó. Además: el error
claro cuando falta la ubicación, con el código de entrada y no el de esquema
(con una ruta absoluta en el TOML, para que el juez antiguo la hubiera
encontrado sin que se la dijeran); el nombre del fichero cuando falta el TSV
junto a los metadatos; y el rechazo de otra matriz con los mismos nombres de
fichero, de unos metadatos con otros bytes (solo la huella de los metadatos los
separa), de un par híbrido y de una auditoría de normalización falsificada a
la vez en la configuración efectiva y en el manifiesto.
Mutaciones sobre copia, todas muertas salvo una documentada: el juez que
ignora el parámetro y usa la ruta del TOML, el que no exige el dato, el que
no compara la huella de los metadatos, el que no compara la auditoría
enlazada. **La comparación de `matrix.sha256` de los metadatos con el TSV no la
discrimina ninguna prueba** (mutante equivalente): las comprobaciones
posteriores de bytes canónicos y de auditoría rechazarían igualmente un TSV
retocado o el de otra matriz. Es una guarda anterior, no tocada aquí.

Dos precisiones. Una ruta relativa para la ubicación se resuelve contra el
directorio de trabajo del proceso que llama: conviene pasar una absoluta, y
está dicho en la documentación de las dos funciones. Y el error de ubicación
ausente usa `GPS-INPUT-001`, el código más cercano del catálogo, porque este
no tiene uno de «argumento ausente» y añadirlo cambiaría un catálogo
congelado.

## DT-023. Resolver una referencia desde un manifiesto juzga la verdad entera de aquel paquete

**Gravedad: baja. Cerrada por aceptación el 04/10/2026**, con una condición
para reabrirla.

### Qué pasa

`resolve_manifest_reference` (`src/results.rs`) llama a
`validate_published_resources` sobre el paquete del que toma la referencia, así
que resolver una referencia cuesta lo que juzgar la verdad entera de aquel
paquete. No es un fallo de corrección. **El coste no es nuevo**: la base
(`9c1fb5f`) ya leía la verdad entera en memoria, solo para su cabecera;
`99d70ac` añadió juzgarla.

### Medidas

`/usr/bin/time -v` sobre paquetes publicados con el simulador del árbol de
`99d70ac`, sondas de `resolve_manifest_reference` compiladas de `9c1fb5f` y
de `99d70ac`, alternadas:

| Paquete | Verdad | Tiempo base | Tiempo actual | Memoria máxima base | actual |
|---|---|---|---|---|---|
| 2 Mb, 60 diploides, 1000 filas (3 repeticiones) | 65 MB | 7,2 a 8,0 s | 8,3 a 8,5 s | 106,8 MB | 177,7 MB |
| 2 Mb, 120 diploides, 2500 filas (2 repeticiones) | 329 MB | 1:46 y 2:04 | 1:39 y 2:26 | 421 MB | 780 MB |

En el primero el tiempo sube un 16 % y la memoria un 66 %; en el segundo la
memoria sube un 85 % y el tiempo es ruido (la dispersión entre repeticiones
supera la diferencia). Resolver sin juzgar costaría milisegundos, pero perdería
la garantía que DTC-060 busca: que la referencia, la verdad y las relaciones
del paquete son coherentes antes de entregarla.

### Decisión

Aceptarla. Ninguna orden de la CLI usa la función hoy: solo la usan las
pruebas, con paquetes pequeños. **Reconsiderarla si una orden de la CLI llega
a usarla con paquetes grandes**; la alternativa entonces sería un juez ligero,
solo de la referencia y del manifiesto, lo que exigiría decidir qué
comprobaciones bastan. Documentado también en la documentación de la función.
