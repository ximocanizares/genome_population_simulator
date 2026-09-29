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
| DT-001 | Dos tercios de las comprobaciones de los verificadores independientes no las ejercita ninguna prueba | **Grave** | 21/09/2026 | Abierta, muy reducida (23/09/2026) |
| DT-002 | Dos pruebas ignoradas están vaciadas, no solo apagadas | Media | 21/09/2026 | Abierta |
| DT-003 | Ocho pruebas ignoradas corresponden a contratos que ya no existen | Baja | 20/09/2026 | Abierta |
| DT-004 | La puerta de Python no corre desde el 25/08/2026 | Media | 19/09/2026 | Abierta |
| DT-005 | La marca del lector de verdad que distingue un tramo sustituido por un tracto no la ejercita ninguna prueba | Baja | 22/09/2026 | Abierta |
| DT-006 | El sorteo del inserto por rechazo es inutilizable cuando la copia más larga acerca el máximo al mínimo del perfil | Media | 24/09/2026 | Abierta, con guarda de admisión (25/09/2026) |
| DT-007 | Ningún diagnóstico público nuevo puede publicarse en un manifiesto del recorrido general | Baja | 25/09/2026 | Abierta |
| DT-008 | El orden entre el paso de linaje y la decisión de base no lo separaba ninguna prueba | Baja | 26/09/2026 | **Cerrada (26/09/2026)** |
| DT-009 | Las dos fuentes de la calibración discrepan en cuánto multiplica la PCR el tartamudeo de los homopolímeros | Media | 27/09/2026 | Abierta |

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
| Oráculo de aceptación del paquete (`src/f5_oracle.rs`) | 365 | 252 | 365 | **168** |
| Oráculo semántico del VCF (`src/vcf_semantic_validation.rs`) | 77 | 39 | 98 | **25** |
| Auditoría de identidad y publicación (`src/f5_closure.rs`) | 42 | 11 | 42 | **5** |

**La puerta solo mira estos tres módulos**, y conviene tenerlo presente al
leer las cifras: no mira `src/results.rs` ni `src/event_manifest.rs`, que son
los que comprueban la configuración efectiva y el manifiesto de evento. El
27/09/2026 una revisión adversarial encontró ahí ciento treinta líneas de
guarda **sin ejercitar como rechazo por ninguna prueba** —sustituirlas enteras
por una línea dejaba pasar la batería completa—, y esta puerta no podía
enterarse porque no mira ese fichero. Se cubrieron con una prueba de mutación;
ampliar la puerta a esos dos módulos sigue sin decidir.

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
porqué. En el verificador del VCF lo están veintiuna de las veinticinco; las
**cuatro** restantes llevan escrito que **no se ha demostrado que sean
inalcanzables** y que quedan pendientes, que es más honrado que suponerlo: las
dos del cálculo del genotipo esperado y dos del análisis de una secuencia
FASTA.

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

**La puerta solo mira tres módulos, y eso deja verificadores fuera del radar.**
Se anota aquí el 26 de septiembre de 2026 porque la entrega 4 de microsatélites
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

**168 comprobaciones del oráculo grande y 4 del verificador del VCF.** Las
cuatro del VCF son las únicas de ese módulo de las que no se ha demostrado
nada: la muestra sin copias en el contig, el alelo alternativo ausente del
registro esperado, la identidad de secuencia que el analizador de cabeceras
acepta y el constructor rechaza, y el fallo de construcción de una secuencia
cuyas tres piezas son válidas.

De las del oráculo grande, 21 son autocomprobaciones
contra catálogos compilados —inalcanzables a propósito, y las veintiuna
muertas— y las 147 restantes se
reparten sobre todo entre la evaluación del desequilibrio de ligamiento (13),
la auditoría del paquete fallido (11), la lectura del reparto físico (10), el
recorrido de fundadores y copias, la lectura de haplotipos y la auditoría del
muestreo.

Sigue sin medirse lo mismo que entonces: las fases 3 y 4 al completo, la
admisión del productor y el contrato de la primera versión de la verdad.

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

### Qué queda abierto

Nada de esta deficiencia. Queda anotado, eso sí, que **el techo de la media es
una guarda de la frontera y no del método**: `poisson-count-v1` sigue admitiendo
medias hasta `2^52` para sus demás consumidores.

---

## DT-011. Dos rechazos del generador de alelos no se alcanzan desde una configuración

**Anotada** el 27 de septiembre de 2026 por la revisión de cierre. No es un
defecto: es una laguna de cobertura **declarada**, con la forma que este
proyecto ya usó para otras.

Los dos rechazos son:

- **El tramo de la referencia sin ni una copia completa del motivo.** Con
  referencia suministrada lo caza antes la guarda de admisión de los loci, que
  exige una tirada mínima —tres unidades por omisión—; con referencia simulada
  esa guarda no corre, pero el tracto lo pone el propio simulador.
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
