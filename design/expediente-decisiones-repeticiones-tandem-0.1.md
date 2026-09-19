# Expediente de decisiones — repeticiones en tándem anticipadas y tartamudeo de PCR 0.1

## 1. Qué es este documento

Este expediente prepara **decisiones**, no las toma. Responde a una petición
concreta: tener pronto microsatélites que varíen entre individuos y ver el
tartamudeo de la PCR en las lecturas, **sin implementar antes toda la fase de
variación compleja**.

Contiene tres cosas:

- **Parte A.** El conjunto más pequeño de trabajo que hace falta, qué queda
  fuera y qué se rompería si se recortara más.
- **Parte B.** Una ficha por cada decisión abierta que ese corte necesita.
- **Parte C.** El borrador del requisito del tartamudeo, que hoy no existe en
  ningún documento del proyecto.

**Alcance de este fichero.** No modifica ningún otro documento. Todas las
correcciones que propone a la especificación, a las decisiones técnicas, al
plan de la fase de variación compleja o a la hoja de ruta están escritas aquí
como propuestas, con su ubicación exacta, para que se apliquen —o no— después
de decidir.

**Cómo leerlo.** Las fichas de la parte B están numeradas `D-01`…`D-18`. Cada
una dice qué hay que decidir, qué opciones hay, qué se recomienda y qué se
bloquea si no se decide. Si solo hay tiempo para leer tres, léanse
[D-12](#d-12), [D-02](#d-02) y [D-01](#d-01).

---

## 2. Resumen: lo esencial en una página

**El corte mínimo es mucho más pequeño de lo que el plan de la fase de
variación compleja sugiere**, y la razón es que tres piezas que ese plan da por
necesarias ya existen en el código:

| Lo que el plan supone | Lo que hay en el código | Consecuencia |
|---|---|---|
| Representar alelos de repetición por haplotipo depende del primer paso del hito de variantes estructurales, porque allí vive la maquinaria que permite que cada copia lleve secuencia distinta | `src/segments.rs` ya tiene esa maquinaria: `SegmentChain` con modificaciones tipadas `Snp`, `Insertion` y `Deletion` por copia haplotípica, y `materialize_copy_interval` ya las materializa | **No hace falta el hito de variantes estructurales.** Un alelo de repetición que difiere en unidades completas es, mecánicamente, una inserción o una deleción anclada en el tracto |
| Incluir repeticiones en una referencia simulada depende del paso que deriva una referencia estructural | `src/reference_generation.rs` ya genera referencias simuladas y `src/reference_fasta.rs` ya escribe el FASTA canónico de 80 bases por línea | **No hace falta derivar referencias estructurales.** Basta con una referencia suministrada que ya contenga el tracto, o con una extensión pequeña del generador que ya existe |
| El VCF de repeticiones necesita el fichero complementario y el catálogo de códigos de pérdida | El proyecto emite `VCFv4.5` (`src/vcf_contract.rs`), y esa versión del estándar ya define los campos de repetición en tándem | **No hace falta el complemento**, siempre que se usen esos campos: si el estándar lo expresa, no hay pérdida que complementar |

**El tartamudeo choca con un requisito cerrado.** La fase de realismo técnico
declara, de forma expresa y en cinco sitios, que su preparación de biblioteca
es **sin PCR**: REQ-TEC-07, REQ-TEC-08, REQ-CONF-88, REQ-DAT-46 y ACE-49, con
la preparación amplificada aplazada en DTF-056. Un artefacto de PCR no cabe hoy
en ese contrato. Es la decisión que más bloquea y está en la ficha
[D-12](#d-12).

**Cuántas decisiones hay que cerrar.** Dieciocho fichas. Quince son decisiones
abiertas de verdad; tres son confirmaciones de consecuencias que otras
decisiones ya tomadas imponen. De las quince, cuatro proceden del plan de la
fase de variación compleja y once son nuevas, descubiertas al leer el código.

**Qué recomiendo hacer primero.** Cerrar [D-12](#d-12), [D-02](#d-02) y
[D-01](#d-01), en ese orden, y después empezar por el modelo de dominio del
locus, que no depende de nada.

---

# Parte A — El corte mínimo viable

## 3. Qué pide el usuario, traducido a piezas

Dos cosas:

1. **Microsatélites en el genoma que varíen entre individuos.** Es decir: un
   locus de repetición en coordenadas conocidas de la referencia, y alelos
   —números de copias— que difieren entre individuos y entre las dos copias
   haplotípicas de un mismo individuo.
2. **Que la PCR cambie el número de repeticiones observado.** El artefacto
   técnico: la verdad del individuo no cambia, cambian las lecturas.

Nada de esto exige translocaciones, inversiones, duplicaciones, ploidía
generalizada, referencias de mapeo múltiples ni el complemento del formato VCF.

## 4. Qué entra en el corte mínimo

Siete piezas. Las tres primeras son el núcleo; sin ellas no hay microsatélite.

### Pieza 1 — Modelo de dominio del locus de repetición

Un locus de repetición es: secuencia, intervalo en la referencia, motivo
(A/C/G/T, longitud libre, sin `N`), número de copias de referencia, y patrón
de interrupciones realizado.

Valida antes de simular: motivo vacío rechazado, motivo con `N` rechazado,
longitud declarada que no coincide con la del motivo rechazada, cada error
identificando el locus (REQ-CONF-60, REQ-IND-34, REQ-REF-23).

**No depende de nada.** Es código puro sobre secuencias en memoria. Se puede
empezar hoy, en cuanto se cierre [D-06](#d-06).

Corresponde al paso F6.4.1 del plan de la fase de variación compleja, sin
cambios.

### Pieza 2 — El alelo por copia haplotípica

Cada copia haplotípica de cada individuo lleva su número de unidades en cada
locus. La diferencia respecto a la referencia se representa como una
modificación pequeña tipada sobre la cadena de segmentos: una inserción de
`k` unidades completas, o una deleción de `k` unidades completas.

Esto es lo que hace que los individuos **varíen**. Es la pieza que el usuario
más quiere y, como dice la tabla del resumen, **no necesita el hito de
variantes estructurales**: `src/segments.rs` ya sostiene modificaciones
tipadas por copia, y `src/materialization.rs` ya las materializa.

Lo que sí hace falta decidir es si el alelo se marca como «repetición» o viaja
como un *indel* corriente; ficha [D-02](#d-02).

Corresponde al paso F6.4.3 del plan, **con su dependencia del hito estructural
eliminada**.

### Pieza 3 — La configuración

Declarar los loci y declarar los alelos. Dos cosas distintas:

- el **catálogo de loci**, que es propiedad de la referencia y se declara una
  sola vez;
- los **alelos por individuo y copia**, que son propiedad de cada individuo.

La vía directa ya admite cambios explícitos por individuo ya faseado
(`direct_events` en el esquema general v2, con `kind` en `snp`, `insertion`,
`deletion`). Añadir una clase de repetición es la extensión más pequeña
posible. Fichas [D-01](#d-01) y [D-04](#d-04).

### Pieza 4 — La verdad y el VCF

La verdad publicada (`truth.jsonl`) y el VCF de verdad deben mostrar el alelo
de cada copia con su fase. El aparato ya existe: `GenotypeCaseKind::Indel`,
`TruthAllele::Insertion` y `TruthAllele::Deletion` ya viven en
`src/truth_schema.rs`, y la normalización de la proyección VCF ya agrupa los
alelos de un mismo sitio en un único registro multialélico con genotipo
faseado (DTC-064).

Lo que falta es decidir si eso basta o si el número de copias y el motivo
tienen que aparecer explícitamente. Fichas [D-08](#d-08) y [D-09](#d-09).

### Pieza 5 — Configuración efectiva y manifiesto

REQ-MAN-29 ya exige registrar los loci y parámetros configurados, su
procedencia concreta o aleatoria, el motivo, la unidad, las copias o la
longitud, y los valores realizados de pureza e interrupciones. Es trabajo
mecánico una vez cerrado el contrato de configuración.

### Pieza 6 — Control en la interfaz Marimo

La hoja de ruta, sección 7, dice que una fase no se cierra si una opción
configurable carece de control Marimo y de equivalencia interfaz → TOML →
configuración efectiva. Si este corte se declara hito, arrastra su control.
Ficha [D-17](#d-17).

### Pieza 7 — El tartamudeo

Un paso nuevo en el recorrido de lecturas, entre la formación de la quimera y
la transformación en par observado. Toda la parte C del expediente.

## 5. Qué queda fuera, con destino

| Queda fuera | Dónde estaba | Por qué no hace falta |
|---|---|---|
| Duplicaciones, deleciones estructurales, inversiones y translocaciones | Hitos de referencias estructurales y de variantes estructurales del plan de la fase de variación compleja | Un alelo de microsatélite no es una variante estructural mientras su diferencia quede por debajo del umbral de 50 pb |
| Referencias estructurales derivadas, su mapa, su linaje y su comprobador | Hito de referencias estructurales | El tracto puede vivir en la referencia suministrada o en la simulada que ya se genera hoy |
| Ploidía generalizada por secuencia | Hito de ploidía generalizada | Diploide global basta para que un individuo tenga dos alelos distintos |
| El modelo escalonado simple de mutación | Paso F6.4.4 | Sirve para que las repeticiones **muten** entre generaciones. El usuario pide que **varíen entre individuos**, que es otra cosa: se consigue declarando alelos distintos |
| Mutaciones internas de pureza, interrupciones o motivo | Paso F6.4.5 | REQ-CONF-63 las exige desactivadas de forma predeterminada. Desactivadas no cuestan nada |
| Generación aleatoria de loci por regiones | Paso F6.4.6 | Con loci declarados explícitamente no se consume aleatoriedad ninguna |
| Fichero complementario del VCF y catálogo de códigos de pérdida | Hito de variantes complejas y complemento VCF | Solo hace falta si el VCF pierde detalle. Con los campos de repetición del estándar 4.5 no lo pierde. Ficha [D-09](#d-09) |
| Conjuntos y matrices por referencia de mapeo, BAM, CRAM y mapeo real | Hito de conjuntos y matrices por referencia | El usuario quiere ver el tartamudeo en las lecturas; la fase de realismo técnico termina en FASTQ |
| Fundadores de la fase de población inicial con alelos de repetición | — | El corte usa la vía directa. Ficha [D-18](#d-18) |

## 6. Qué se rompería si se recortara más

Cuatro recortes parecen tentadores y no funcionan. Conviene dejarlo dicho para
no descubrirlo a mitad de la implementación.

**Recorte 1: «que el alelo sea un *indel* corriente y ya está».**

Funciona para la secuencia, pero rompe dos cosas. Primera, REQ-IND-33 exige
que el estado genotípico conserve, por haplotipo, la fase, la secuencia
realizada **y el número de copias o longitud total**; un *indel* anónimo no
conserva el número de copias. Segunda, y más grave: **el tartamudeo necesita
saber dónde están los tractos**. Si el alelo es un *indel* anónimo y el locus
no está declarado, el motor de lecturas no tiene forma de saber que esas bases
son un microsatélite. El tartamudeo quedaría sin definir.

**Recorte 2: «que el tartamudeo se aplique lectura a lectura».**

Rompe la coherencia del par. Las dos lecturas de un par salen de la misma
molécula. Si ambas cubren el tracto —cosa que ocurre cuando el inserto mide
menos del doble de la lectura, y el perfil vigente sortea insertos entre 250 y
450 pb con lecturas de 150, así que ocurre en toda la cola corta de la
distribución— deben mostrar **el mismo** recuento. Un
tartamudeo por lectura produce pares incoherentes que en datos reales no
existen, y que además desconcertarían a cualquier herramienta de genotipado de
microsatélites. Está en el requisito propuesto REQ-DAT-75 de la parte C.

**Recorte 3: «no tocar la configuración, declarar los loci en un fichero
aparte y punto».**

El proyecto no permite que un valor que cambia el comportamiento observable
quede fuera de la configuración efectiva y del manifiesto. Un locus no
declarado en la configuración efectiva hace la ejecución irreproducible y
rompe REQ-MAN-29. Sí es posible que el **contenido** viva en un fichero
externo, como ya hace la matriz de frecuencias de SNV; ficha [D-04](#d-04).

**Recorte 4: «meter el tartamudeo dentro del muestreo del fragmento».**

Esto sí lo he comprobado ejecutando la lectura del código, y **falla**. En
`src/chimera.rs`, la función que forma la quimera vuelve a materializar el
prefijo del fragmento desde la copia haplotípica y lo compara **byte a byte**
con el prefijo del fragmento que recibe; si no coinciden devuelve
`ChimeraError::PairSourceMismatch`. Si el tartamudeo hubiera modificado ya el
fragmento, esa comprobación fallaría cada vez que una molécula con tartamudeo
resultara elegida para formar quimera. El sitio correcto es **después** de la
quimera; ficha [D-13](#d-13).

## 7. Conclusión de la parte A

El corte es posible y es pequeño. **No arrastra ninguno de los cinco hitos
grandes** de la fase de variación compleja. Arrastra:

- un paso nuevo y aislado (modelo de dominio del locus);
- una extensión pequeña de una estructura que ya existe (la cadena de
  segmentos, sin clases nuevas de segmento);
- una versión nueva del esquema general y de la configuración efectiva;
- una extensión de la cabecera y la proyección del VCF;
- un paso nuevo en el recorrido de lecturas, que **reabre el contrato de una
  fase cerrada**.

Ese último punto es el precio real del corte, y no se puede evitar: el
tartamudeo vive en el motor de lecturas, y el motor de lecturas está cerrado.

---

# Parte B — Fichas de decisión

Dieciocho fichas. Las marcadas **(confirmación)** no son decisiones abiertas:
son consecuencias que otras decisiones ya tomadas imponen, y se listan para que
nadie las descubra tarde.

---

## Grupo 1 — El corte y su contrato

<a id="d-01"></a>
### D-01. Dónde se declaran los loci y los alelos, y qué versión de esquema nace

Procede de DP-F6-01 del plan de la fase de variación compleja, recortada a lo
que este corte necesita.

**Qué hay que decidir.** El fichero general de configuración tiene hoy dos
versiones vivas: la versión 2, que es la vía directa —individuos declarados
uno a uno con sus cambios explícitos—, y la versión 3, que es la vía de la
población inicial sintética y exige un escenario enlazado. Ninguna de las dos
conoce las repeticiones. Los esquemas publicados no se reabren: hay que crear
una versión nueva. Falta decidir cuál y con qué forma.

**Opciones.**

| Opción | Qué implica | Consecuencia práctica |
|---|---|---|
| (a) Esquema general **v4** que extiende la vía directa: una tabla nueva `tandem_repeat_loci` y una clase nueva `tandem_repeat` dentro de `direct_events` | Un solo fichero, una sola versión nueva, sin escenario enlazado | La más pequeña. El usuario escribe un TOML y ejecuta. Choca con el techo de 65536 bytes si hay muchos individuos (ver [D-04](#d-04)) |
| (b) Un escenario de repeticiones propio, enlazado desde el general, al estilo del escenario de la población inicial | Dos ficheros, dos esquemas nuevos | Más ceremonia. Gana sitio frente al techo de tamaño, pero el techo se aplica también al escenario |
| (c) Esperar a que se cierre el contrato completo de la fase de variación compleja | Nada nuevo ahora | Bloquea el encargo entero. Es la opción de no hacer el corte |

**Recomendación: (a).** Es la extensión más pequeña que existe y no crea una
segunda vía de configuración que después haya que reconciliar. La forma
concreta que propongo:

- `[[tandem_repeat_loci]]` — catálogo del locus: `locus_id`, `contig`,
  `start_1based`, `end_1based`, `motif`, y las interrupciones si se declaran
  (ver [D-06](#d-06)). Es propiedad de la referencia, no del individuo.
- En `direct_events`, `kind = "tandem_repeat"` con `locus_id`, `phase` y
  `units`. Reutiliza el selector de individuo y fase que ya existe.

**Qué se bloquea si no se decide.** Todo lo que escribe configuración,
configuración efectiva, manifiesto o esquema. Es decir: todo menos el modelo de
dominio del locus, que es puro.

---

<a id="d-02"></a>
### D-02. Cómo se representa el alelo dentro del estado genotípico

Decisión nueva. **Es la que determina si el corte mínimo es pequeño o es toda
la fase de variación compleja.**

**Qué hay que decidir.** Un individuo con 7 unidades donde la referencia tiene
5 lleva 2 unidades de más en esa copia haplotípica. ¿Cómo se guarda eso?

**Opciones.**

| Opción | Qué implica | Qué arrastra |
|---|---|---|
| (a) **Modificación pequeña tipada**: una variante nueva `SmallModification::TandemRepeat` junto a `Snp`, `Insertion` y `Deletion`, que lleva el `locus_id`, el número de unidades y el ancla, y que se materializa como la inserción o la deleción equivalente | Cambio contenido en `src/segments.rs` y `src/materialization.rs` | **Nada.** La cadena de segmentos ya sostiene modificaciones tipadas por copia |
| (b) **Inserción o deleción corriente** más una tabla lateral que dice «este *indel* es el locus X» | Cambio aún menor | Pierde el número de copias en el estado, contra REQ-IND-33, y deja el vínculo locus↔*indel* fuera de la estructura que la verdad proyecta |
| (c) **Clase nueva de segmento estructural**, como propone el primer paso del hito de variantes estructurales | Cambio profundo en la cadena de segmentos, con identidad propia por copia, orientación y vínculo al intervalo fundador | **Arrastra el hito de variantes estructurales entero**, y con él la verdad versión 3 y el reconstructor independiente |

**Recomendación: (a).** La opción (c) es lo que el plan de la fase de variación
compleja da por supuesto, y es innecesaria: una repetición en tándem no cambia
la topología de la copia, solo su longitud en un punto. La opción (b) ahorra
media hora y pierde un requisito.

Conviene decir en voz alta la consecuencia de (a): el alelo de repetición es
**una inserción o una deleción con etiqueta**. Materializa igual, se proyecta
igual, se valida igual. Lo único que añade es saber que esas bases son un
tracto y cuántas unidades tiene, que es justo lo que hacen falta para la verdad
y para el tartamudeo.

**Qué se bloquea si no se decide.** La pieza 2 entera, que es el corazón del
encargo. Y, por dependencia, la verdad, el VCF y el tartamudeo.

---

<a id="d-03"></a>
### D-03. Dónde vive el tracto en la referencia

Decisión nueva.

**Qué hay que decidir.** Para que un individuo tenga 7 unidades y otro 5, la
referencia tiene que tener un número de unidades en ese sitio. ¿De dónde sale?

**Opciones.**

| Opción | Qué implica |
|---|---|
| (a) **Referencia suministrada que ya lo contiene.** El usuario da un FASTA y declara el locus por coordenadas. El sistema comprueba que el intervalo declarado es efectivamente el motivo repetido un número entero de veces (con sus interrupciones declaradas) y rechaza si no | Cero código nuevo de referencia. Una comprobación de admisión |
| (b) **Referencia simulada con el tracto sembrado.** Extensión de `src/reference_generation.rs`, que ya existe y ya genera secuencias a partir de frecuencias de base | Código nuevo pequeño y bien acotado. No exige derivar referencias estructurales |
| (c) Las dos | Suma de las dos |

**Recomendación: (c), con (a) primero.** La opción (a) es gratis y desbloquea
el trabajo inmediatamente; permite probar el corte con un FASTA de juguete
escrito a mano. La (b) es cómoda para generar escenarios sin buscar un genoma
real y cuesta poco, porque el escritor FASTA canónico —80 bases por línea,
mayúsculas, salto de línea `LF`— ya existe en `src/reference_fasta.rs`.

Lo que **no** hace falta en ninguna de las dos es el paso que deriva una
referencia estructural a partir de otra, que es de lo que el plan las hacía
depender.

**Qué se bloquea si no se decide.** No bloquea el modelo de dominio ni el
alelo. Bloquea poder ejecutar un caso completo de punta a punta.

---

<a id="d-04"></a>
### D-04. Cómo se asignan los alelos a los individuos

Decisión nueva, con una trampa de tamaño detrás.

**Qué hay que decidir.** Declarar el alelo de cada copia de cada individuo en
cada locus, ¿a mano o por frecuencias?

**La trampa.** DTC-089, el contrato que fija cuánto puede ocupar un fichero de
configuración, impone un techo duro de 65536 bytes al fichero general. Una
declaración explícita ocupa del orden de 70 a 90 bytes por entrada. Con 100
individuos diploides y 5 loci son 1000 entradas: entre 70 y 90 kB. **El techo
se rompe antes de llegar a un estudio de tamaño normal.**

**Opciones.**

| Opción | Qué implica | Techo |
|---|---|---|
| (a) **Explícito en el general.** Una entrada por individuo, copia y locus | Nada nuevo más allá de [D-01](#d-01) | Cómodo hasta unas diez muestras; imposible a partir de unas cincuenta |
| (b) **Fichero externo de alelos**, al estilo de la matriz de frecuencias de SNV, que ya admite los modos `supplied`, `reuse` y `generate` con ruta a un fichero | Un esquema y un lector más | Sin techo práctico |
| (c) **Frecuencias alélicas por locus y población**, y el sistema sortea el genotipo de cada individuo | Consume aleatoriedad: hace falta un ámbito de semilla y un algoritmo normativo. `categorical-f64-v1` y `hwe-composition-v1` ya existen en la familia de consumidores de DTC-072 | Sin techo. Es además lo que un genetista querría para simular una población |

**Recomendación: (a) ahora, (c) inmediatamente después.** La opción (a) no
sirve para un estudio, pero desbloquea hoy mismo todo el trabajo de las piezas
1 a 4 y del tartamudeo, con un caso de dos individuos escrito a mano. La (c) es
la que el usuario acabará usando, y conviene no cerrar la forma del TOML de (a)
de manera que impida añadir (c) después.

La (b) solo hace falta si se quieren alelos individuales en masa, que es un uso
menos frecuente que las frecuencias.

**Qué se bloquea si no se decide.** Nada del núcleo. Pero si se elige (a) y se
publica, y luego hace falta (c), sale otra versión de esquema.

---

<a id="d-05"></a>
### D-05. Qué pasa cuando el alelo se separa más de 50 pb de la referencia

Decisión nueva.

**Qué hay que decidir.** REQ-CONF-55 y REQ-IND-31 fijan un umbral
configurable, 50 pb por omisión, que separa las inserciones y deleciones
pequeñas de las variantes estructurales, y exigen aplicarlo **igual al generar,
al clasificar y al contar**. Un alelo de repetición que difiere en 60 pb de la
referencia cae, por esa regla, del lado estructural. ¿Qué se hace?

**Opciones.**

| Opción | Consecuencia |
|---|---|
| (a) **Rechazar antes de simular** todo alelo cuya diferencia con la referencia alcance el umbral, identificando el locus, el individuo y la copia | El corte queda íntegramente dentro del territorio de los *indels* pequeños. Para microsatélites de 1 a 6 pb con rangos alélicos normales sobra: un tetranucleótido admite ±12 unidades |
| (b) **Admitirlo y clasificarlo como variante estructural** | Arrastra la representación simbólica en VCF, y con ella el hito de variantes complejas |
| (c) **Declarar que las repeticiones no se someten al umbral** | Contradice REQ-IND-31, que exige aplicar el umbral de forma coherente. Habría que corregir la especificación |

**Recomendación: (a).** Es la restricción que mantiene el corte pequeño, es
honesta —dice lo que no hace y falla pronto— y no es limitante para
microsatélites, que es lo que el usuario ha pedido. Las repeticiones largas
quedan con destino: el hito de variantes complejas de la fase de variación
compleja.

**Qué se bloquea si no se decide.** La validación de admisión, y con ella el
primer caso que se ejecute con un alelo grande dará un resultado no acordado.

---

## Grupo 2 — Cómo se realiza la secuencia

<a id="d-06"></a>
### D-06. Pureza e interrupciones: qué se declara y qué se calcula

Procede de DP-F6-07. **El usuario ha pedido atención especial a esta ficha.**

**El hueco exacto.** La especificación define la pureza en su glosario: «la
proporción del tracto repetido que sigue exactamente el motivo; la parte
restante corresponde a interrupciones». REQ-CONF-57 permite indicarla y
REQ-CONF-59 permite definirla por locus o por distribuciones. **Nadie ha dicho
dónde caen las interrupciones ni qué bases toman.** Son dos huecos distintos y
hay un tercero, que nadie ha planteado: qué le pasa al patrón de
interrupciones cuando el alelo cambia de longitud.

Un ejemplo concreto del problema. Motivo `AC`, 10 unidades, pureza 0,9. El
tracto mide 20 bases; el 10 % son 2 bases interrumpidas. ¿Cuáles? ¿La 7 y la
13? ¿Las dos seguidas? ¿Y qué bases llevan: `G`, `T`, una al azar? Dos
implementaciones razonables producen dos FASTA distintos, y ninguna de las dos
es más correcta que la otra según lo escrito hoy.

**Opciones.**

| Opción | Qué implica | Coste |
|---|---|---|
| (a) **Interrupciones declaradas, pureza calculada.** El usuario declara la lista de interrupciones como pares (desplazamiento dentro del tracto, base), y el sistema calcula y publica la pureza realizada | No hay que inventar ninguna regla de colocación ni ninguna elección de base. Cero aleatoriedad. El FASTA es exactamente lo que el usuario escribió | Mínimo |
| (b) **Pureza declarada, colocación determinista y versionada.** El usuario da un número; un algoritmo reparte las interrupciones de forma regular (por ejemplo, repartidas lo más uniformemente posible sobre las unidades, y siempre en el mismo desplazamiento dentro de la unidad) y sustituye por una base fija derivada del motivo | Hay que fijar un algoritmo y versionarlo. Sin aleatoriedad, reproducible | Medio |
| (c) **Pureza declarada, colocación sorteada.** Posiciones uniformes sobre el tracto y sustitución uniforme entre las tres bases distintas de la del motivo en ese desplazamiento | Hay que fijar el algoritmo, abrir un ámbito de semilla propio y congelar vectores de prueba | Alto |

**Recomendación: (a) para el corte mínimo, y (b) después si hace falta.**

La razón es la que el propio proyecto repite: la pregunta «¿no sería más
fácil…?» es una corrección, no una duda. La pureza como **entrada** obliga a
inventar una regla de colocación que nadie ha pedido y que no cambia nada de
lo que el usuario quiere ver. La pureza como **salida** cumple la definición
del glosario, cumple REQ-MAN-29 —que pide registrar «los valores realizados de
pureza e interrupciones»— y elimina de un plumazo un algoritmo, un ámbito de
semilla y un juego de vectores congelados.

La opción (c) además tiene un defecto silencioso propio: si las interrupciones
se sortean, un cambio de semilla mueve la secuencia de la referencia y, con
ella, todo lo demás.

**Redacción propuesta para el glosario** (sección 4.3 de la especificación,
párrafo de la repetición en tándem), como añadido a lo que ya dice:

> Las interrupciones de un tracto se declaran de forma explícita como pares de
> desplazamiento dentro del tracto y base sustituta. La pureza es un valor
> **realizado** que el sistema calcula y publica, no un parámetro de entrada de
> la primera versión. Una configuración que declare una pureza objetivo debe
> rechazarse antes de simular, identificando el locus y el parámetro, con
> destino a la ampliación que fije la regla de colocación.

**Qué se bloquea si no se decide.** El modelo de dominio del locus, que es el
primer paso y el único que no depende de nada más. Es decir: **sin esta ficha
no se puede empezar**.

---

<a id="d-07"></a>
### D-07. Por qué extremo se ganan o se pierden unidades

Decisión nueva, pequeña y con consecuencias grandes.

**Qué hay que decidir.** Un alelo de 7 unidades donde la referencia tiene 5
lleva 2 unidades de más. ¿Dónde se insertan? ¿Y de dónde se quitan cuando son
menos?

**Por qué importa.** Si el patrón de interrupciones no es uniforme —y en un
microsatélite real casi nunca lo es— el sitio donde se añade o se quita cambia
la secuencia resultante. Peor: si cada alelo se **regenera** desde cero a
partir del motivo y las interrupciones, el alelo de 7 unidades y el de 5 pueden
diferir en sitios que no son el extremo, y entonces la diferencia deja de ser
una inserción o una deleción limpia. Eso rompería la representación elegida en
[D-02](#d-02) y complicaría el VCF sin necesidad.

**Opciones.**

| Opción | Consecuencia |
|---|---|
| (a) **Extremo derecho del tracto, unidades completas.** Las unidades se añaden o se quitan justo antes del final del tracto; el resto del tracto queda byte a byte igual al de la referencia | La diferencia es siempre una inserción o una deleción limpia de `k·(longitud del motivo)` bases. El patrón de interrupciones existente no se toca |
| (b) **Extremo izquierdo** | Equivalente, con el ancla en el otro lado |
| (c) **Regenerar el tracto entero desde el motivo y las interrupciones** | La diferencia deja de ser un *indel* limpio en cuanto hay una interrupción. No recomendable |

**Recomendación: (a).** Y decirlo en el requisito, no dejarlo a la
implementación. Es exactamente el tipo de detalle que dos implementaciones
resuelven distinto sin que ninguna prueba lo note.

Nota práctica: la normalización del VCF alinea los *indels* a la izquierda
(DTC-064), así que la posición publicada en el VCF será el ancla izquierda del
tracto **independientemente** de por dónde se haya insertado internamente. Eso
está bien y es lo estándar; ver [D-10](#d-10).

**Qué se bloquea si no se decide.** La materialización del alelo, y con ella
cualquier prueba que compare la secuencia resultante con una escrita a mano.

---

## Grupo 3 — Cómo se publica

<a id="d-08"></a>
### D-08. Forma de la verdad publicada

Procede de DP-F6-03, recortada.

**Qué hay que decidir.** El fichero de verdad `truth.jsonl` publica hoy casos
de tipo SNP y de tipo *indel*, con alelos `Reference`, `Snp`, `Insertion` y
`Deletion`. ¿Un alelo de repetición es un caso *indel* o un caso nuevo?

**Opciones.**

| Opción | Consecuencia |
|---|---|
| (a) **Caso *indel* con metadatos de repetición**: mismos tipos de registro, más los campos `locus_id`, `motif` y `units` en el alelo. Exige una versión nueva del esquema de verdad —llamémosla v3— porque añade campos, pero **no** añade tipos de registro ni toca los reconstructores existentes | Extensión aditiva. El validador de verdad que ya existe sigue funcionando sobre lo que ya sabía |
| (b) **Tipo de caso propio**, `TandemRepeat`, con su alelo propio | Más limpio conceptualmente, más caro: obliga a extender el validador de verdad y el oráculo semántico del VCF con una rama nueva |
| (c) **No tocar la verdad**: publicar el alelo como *indel* anónimo y poner el número de copias solo en el manifiesto | Contradice REQ-VER-22, que dice expresamente que el detalle restante «nunca» debe quedar únicamente en el manifiesto |

**Recomendación: (a).** La verdad ya distingue el caso por su contenido; lo que
falta es el número de unidades y el motivo, que son dos campos. Y crear una
versión nueva del esquema de verdad no es un problema: el plan de la fase de
variación compleja ya prevé que las fases nuevas creen versiones nuevas y dejen
intactas las heredadas. Lo que hay que garantizar es que `truth-v1` y
`truth-v2` no cambien ni un byte, y que una ejecución sin loci de repetición
publique exactamente lo que publica hoy.

**Qué se bloquea si no se decide.** La proyección de la verdad, el VCF y la
validación de cierre.

---

<a id="d-09"></a>
### D-09. Campos del VCF y qué es un «tamaño manejable»

Procede de DP-F6-06.

**Qué hay que decidir.** REQ-VER-20 exige alelos explícitos en `REF` y `ALT`
«cuando su tamaño sea manejable» y una representación simbólica para tramos
largos. REQ-VER-22 exige exponer «todo el motivo, número de copias, longitud,
pureza, interrupciones y fase» mediante campos estándar o adicionales válidos,
y conservar el detalle restante en la entrada complementaria del caso. Nadie ha
decidido qué campos ni dónde está el umbral.

**Por qué esta ficha decide si entra o no el fichero complementario.** El
proyecto emite `VCFv4.5`. La versión 4.4 del estándar introdujo campos
específicos de repetición en tándem —`RU`, `RN`, `RUS`, `RUL`, `RUC`, `RB` y
compañía— y la 4.5 los conserva. **Si se usan, no hay pérdida de detalle, y sin
pérdida no hay entrada complementaria, ni código de pérdida, ni catálogo, ni
esquema del complemento.** Es decir: esta ficha es la que mantiene fuera del
corte el hito entero de variantes complejas y complemento VCF.

**Opciones.**

| Opción | Consecuencia |
|---|---|
| (a) **Campos de repetición del estándar 4.5**, declarados en la cabecera, con alelos explícitos en `REF` y `ALT` | Sin pérdida, sin complemento. Interoperable con herramientas que ya entienden esos campos |
| (b) **Campos propios `INFO/GPS_*`** declarados en la cabecera | Válido como VCF, pero reinventa lo que el estándar ya tiene, y obliga a justificar por qué |
| (c) **Alelos explícitos y nada más**, con el motivo y las copias solo en el manifiesto | Contra REQ-VER-22 |

**Recomendación: (a).** Y, sobre el umbral de «tamaño manejable»: fijarlo en el
mismo umbral que separa *indel* de variante estructural, 50 pb por omisión y
configurable. Así hay **un solo umbral** en todo el sistema en vez de dos que
se puedan desincronizar, y encaja con [D-05](#d-05): lo que este corte admite es
exactamente lo que se publica explícito. Los tramos que exigen representación
simbólica quedan, con destino, en el hito de variantes complejas.

Conviene confirmar sobre el propio estándar 4.5 los nombres, números y tipos
exactos de esos campos antes de escribir la cabecera. El proyecto ya ha pagado
una vez el precio de dar por buena una cifra no comprobada.

**Qué se bloquea si no se decide.** La cabecera del VCF, la proyección y la
puerta de validación externa.

---

<a id="d-10"></a>
### D-10. Dónde cae la posición del alelo en el VCF **(confirmación)**

No es una decisión abierta: DTC-064, que fija la normalización de la proyección
VCF, ya decide «alineamiento a la izquierda contra la referencia». Se anota
para que no sorprenda.

**Consecuencia.** Un alelo de repetición publicado como inserción o deleción
quedará anclado en el **extremo izquierdo del tracto**, no donde se insertó
internamente ([D-07](#d-07)) ni en el inicio declarado del locus si este no
coincide. Es el comportamiento estándar y el que produce `bcftools norm`, así
que la puerta externa lo aceptará.

**Segunda consecuencia.** DTC-064 fija además, para la fase de variantes
pequeñas, **como máximo un registro por combinación de cromosoma y posición**,
con todos los alelos compatibles integrados como varios `ALT` de un único
registro multialélico. Eso es exactamente lo que se quiere: un diploide con 5 y
7 unidades produce **un** registro con dos `ALT` y un genotipo faseado. No hace
falta tocar nada.

---

## Grupo 4 — Aleatoriedad e identidad

<a id="d-11"></a>
### D-11. Ámbitos de semilla e identificadores

Procede de DP-F6-04, muy recortada.

**Qué hay que decidir.** El proyecto exige que toda decisión aleatoria se
derive de un ámbito de semilla propio, identificado por el evento y las
identidades estables pertinentes, y que nunca dependa del orden físico ni de un
flujo compartido (REQ-DAT-15, DTC-037, DTC-081).

**Lo bueno: el corte mínimo casi no usa aleatoriedad.** Con loci declarados y
alelos declarados ([D-04](#d-04) opción (a)) y con interrupciones declaradas
([D-06](#d-06) opción (a)), **no se consume ni un solo número aleatorio** por
las repeticiones. Solo el tartamudeo necesita ámbitos nuevos.

**Qué hay que decidir, entonces.**

1. Si el locus necesita una clase de identificador automático propia, o basta
   con un identificador declarado por el usuario y validado por unicidad.
2. Qué ámbitos abre el tartamudeo y en qué versión del catálogo técnico.

**Recomendación.**

- Identificador del locus **declarado por el usuario** y validado por unicidad
  dentro de la referencia, como ya se hace con los identificadores de
  individuo, población y secuencia. Sin clase de identificador automático
  nueva. Es más simple y el locus es una entidad que el usuario nombra.
- Para el tartamudeo, **un catálogo técnico nuevo, `f4-technical-v4`**, que
  extiende el vigente `f4-technical-v3` (DTC-086) con una etapa `stutter`. Los
  componentes, siguiendo exactamente el patrón de las etapas existentes:
  `[event_id, library_id, stutter, template_ordinal, locus_id, decision]`, con
  `decision` en `bernoulli` y `step`.

Dos detalles que no son adorno:

- El ámbito se indexa por `locus_id`, **no** por la posición del tracto dentro
  de la molécula. Así, añadir un locus en otro cromosoma no desplaza las
  decisiones de los que ya había.
- Cada sorteo arranca de su clave con un generador nuevo, como ya exige
  DTC-081: «no existe consumo acumulado entre decisiones». Esto tiene una
  consecuencia excelente y es la base de la prueba de compatibilidad de
  [D-16](#d-16): **añadir el tartamudeo no puede desplazar ninguna decisión
  anterior**.

**Qué se bloquea si no se decide.** El tartamudeo entero, y la reproducibilidad
de cualquier ejecución que lo use.

---

## Grupo 5 — El tartamudeo

<a id="d-12"></a>
### D-12. ¿Cabe un artefacto de PCR en una fase que declara no usar PCR?

Decisión nueva. **Es la que más bloquea, y es de política del proyecto, no
técnica.**

**El hecho.** La fase de realismo técnico de secuenciación declara de forma
expresa, y en cinco sitios, que su preparación de biblioteca es **sin PCR**:

- **REQ-TEC-07**: «F4 usa exclusivamente preparación WGS sin PCR y la registra
  como tal. Debe rechazar una preparación con amplificación e indicar DTF-056.»
- **REQ-TEC-08**: el perfil amplificado, sus duplicados y sus sesgos quedan
  aplazados.
- **REQ-CONF-88**: rechaza un perfil amplificado y cualquier tasa de duplicados.
- **REQ-DAT-46**: «F4 no genera duplicados: PCR no aplica al perfil inicial.»
- **ACE-49** y **ACE-68**, que lo repiten como criterio de aceptación.
- **DTF-056** registra la preparación amplificada como ampliación futura y dice
  que el perfil vigente es «HiSeq X sin PCR».

**Por qué importa.** El usuario ha confirmado que quiere el **artefacto
técnico**, el tartamudeo de la PCR, no el deslizamiento biológico entre
generaciones. Pero la fase que produce las lecturas dice que no hay PCR.

**El matiz biológico, que es real.** Una biblioteca «sin PCR» evita la
amplificación en la preparación, pero **no** evita toda polimerasa: la
generación de agrupaciones en la celda de flujo es amplificación en puente, y
la propia secuenciación por síntesis es copia con polimerasa. El deslizamiento
de la polimerasa sobre un microsatélite ocurre en ambas, y por eso los
microsatélites muestran tartamudeo también en datos de genoma completo sin PCR,
aunque menos marcado que en amplicones. Es decir: la petición del usuario es
biológicamente legítima incluso dentro del perfil vigente. Lo que hace falta es
decidir **cómo se nombra**, porque el nombre determina si contradice o no un
requisito cerrado.

**Opciones.**

| Opción | Qué implica | Coste documental |
|---|---|---|
| (a) **Nombrarlo por su mecanismo, no por la PCR**: «deslizamiento de la polimerasa durante la amplificación en agrupación y la síntesis», un artefacto del perfil vigente sin PCR | No contradice REQ-TEC-07 ni REQ-CONF-88, porque esos hablan de la **preparación** de biblioteca. Hay que decirlo explícitamente en el requisito nuevo para que nadie lo lea como una grieta | Bajo: un requisito nuevo y una frase aclaratoria |
| (b) **Reabrir DTF-056** y adelantar la preparación amplificada | Trae consigo duplicados de PCR y sesgos propios, que DTF-056 aplaza **por falta de calibración empírica** | Alto. Arrastra un modelo de duplicados que nadie ha calibrado |
| (c) **Declarar el tartamudeo como un modo opcional** que solo se activa con una preparación amplificada, y dejarlo inerte hasta que exista | Coherente, pero no entrega lo que el usuario pide | Medio, y sin resultado |

**Recomendación: (a).** Es la única que entrega el resultado sin abrir un
modelo sin calibrar. Exige, eso sí, ser honesto en la redacción: el requisito
debe decir que este artefacto **no** representa una preparación amplificada,
que no introduce duplicados y que el perfil sigue registrándose como sin PCR.
Y debe decir que, cuando DTF-056 se desarrolle, la tasa de tartamudeo del
perfil amplificado será distinta y mayor.

Sobre el nombre en el chat y en los documentos: el usuario lo llama
tartamudeo, que es el término que usa todo el mundo. Propongo conservar
«tartamudeo» como nombre y añadir en su primera aparición la explicación del
mecanismo, que es lo que la guía de redacción del proyecto pide.

**Qué se bloquea si no se decide.** Toda la parte C. No se puede redactar el
requisito sin saber si contradice o no un requisito cerrado.

---

<a id="d-13"></a>
### D-13. En qué punto exacto del recorrido se aplica

Decisión nueva, con la lectura del código ya hecha.

**Lo que dice el código.** El recorrido de una pareja, hoy:

1. `WgsPairedEndSampler::sample_next_pair` (`src/wgs_sampling.rs`) sortea la
   longitud de inserto, elige la copia haplotípica, elige el inicio y
   **materializa la secuencia del fragmento** con `materialize_copy_interval`.
2. Construye `PairedEndFragment`, que es **la molécula**: una secuencia, su
   copia de origen y su inicio.
3. `emit_paired_end` (`src/paired_end.rs`) corta R1 del extremo izquierdo y R2
   como complemento inverso del extremo derecho. **Aquí deja de haber una
   molécula y pasa a haber dos lecturas.**
4. `chimera::transform_paired_end` (`src/chimera.rs`) decide si la pareja se
   vuelve quimérica y, si sí, **reconstruye un fragmento nuevo** con las dos
   fuentes unidas y vuelve a llamar a `emit_paired_end`.
5. `error_quality::transform_paired_end` (`src/error_quality.rs`) transforma
   cada lectura en su par observado, con errores y calidades por ciclo.

**La lectura del encargo era «entre el muestreo del fragmento y la
transformación en par observado». Es correcta en la dirección y necesita una
precisión: hay dos sitios posibles dentro de ese hueco, y uno de los dos
falla.**

| Sitio | Qué pasa |
|---|---|
| Antes del paso 2, sobre la secuencia recién materializada | **Falla.** En el paso 4, la formación de la quimera vuelve a materializar el prefijo desde la copia haplotípica y lo compara byte a byte con el prefijo del fragmento que recibe; si no coinciden devuelve `ChimeraError::PairSourceMismatch`. Un fragmento ya modificado por el tartamudeo rompería esa comprobación |
| **Después del paso 4 y antes del paso 5** | **Funciona.** La pareja todavía lleva su fragmento completo (`PairedEndPair::fragment()`), así que la molécula sigue disponible. El tartamudeo reescribe la secuencia del fragmento y vuelve a emitir R1 y R2 con `emit_paired_end`, exactamente como ya hace la quimera |

**Recomendación: después de la quimera, antes de calidad y errores.** Un
`stutter::transform_paired_end` con la misma forma que los dos que ya existen.
Esto además es defendible biológicamente: quimera y tartamudeo son ambos
artefactos de la copia con polimerasa, y aplicar el tartamudeo sobre la
molécula quimérica ya formada es tan razonable como al revés.

**Un detalle que hay que decir en el requisito.** Tras el tartamudeo, el
intervalo `start`–`end` del fragmento deja de describir un intervalo contiguo
de la copia haplotípica, porque su longitud ha cambiado. No es un problema
nuevo: la quimera ya deja ese intervalo en esa situación, conservando el
`start` original con una secuencia unida. Y no rompe nada aguas abajo, porque
la transformación de calidad y errores solo usa las secuencias de las lecturas
y las identidades, no las coordenadas del fragmento, y porque REQ-DAT-56 dice
que la fase de realismo técnico no conserva posición ni haplotipo de origen por
lectura.

**Qué se bloquea si no se decide.** La implementación entera del tartamudeo.

---

<a id="d-14"></a>
### D-14. El modelo de probabilidad

Decisión nueva.

**Qué hay que decidir.** Tres hechos de los datos reales que el modelo debe
reproducir: el tartamudeo **crece con la longitud del tracto**, es **más
frecuente cuanto más corto es el motivo** —los mononucleótidos son los peores—
y las **contracciones son más comunes que las expansiones**.

**Propuesta de modelo, `stutter-v1`.** Dos decisiones por molécula y tracto:

1. **¿Tartamudea?** Un ensayo de Bernoulli con probabilidad
   `p = min(p_max, max(0, a + b·(u − u₀)))`, donde `u` es el número de unidades
   completas del tracto **contenidas en esa molécula**, y `a`, `b`, `u₀` y
   `p_max` se configuran **por longitud de motivo**. La dependencia lineal de la
   longitud del tracto es la aproximación habitual y es honesta; una logística
   no añade nada que el usuario haya pedido.
2. **¿Cuánto y en qué dirección?** Un sorteo categórico sobre un conjunto
   declarado de pasos en unidades enteras con sus pesos, por ejemplo
   `{−2, −1, +1}` con pesos configurables. Ahí es donde se expresa que las
   contracciones son más comunes.

**Cómo encaja con lo que el proyecto ya hace.**

- El proyecto exige que los **recuentos** aleatorios salgan exclusivamente del
  algoritmo normativo `poisson-count-v1` (DTC-067). Aquí no hay ningún recuento
  de sucesos sobre oportunidades: hay un ensayo binario y una elección entre
  alternativas. Los consumidores normativos que corresponden ya existen:
  `bernoulli-v1` y `uniform-index-v1` en las primitivas de la fase de realismo
  técnico (DTC-088), y `categorical-f64-v1` en la familia general (DTC-072).
  **Conviene confirmarlo expresamente en la decisión técnica nueva**, para que
  nadie lea la regla de DTC-067 como una prohibición.
- **Sin reintentos.** El proyecto rechaza los reintentos en la selección
  aleatoria: se filtran los candidatos válidos **antes** de sortear, como hace
  `dynamic-indel-selection-v1`. Aquí significa: un paso solo entra en el sorteo
  categórico si el resultado deja al menos una unidad completa en el tracto y
  deja la molécula con longitud suficiente para las dos lecturas. Si ningún
  paso es viable, no hay tartamudeo y **no se consume aleatoriedad**.
- **Cero consumo en los extremos.** `bernoulli-v1` no sortea para probabilidad
  0 ni 1 (DTC-088). Con la tasa a cero, el tartamudeo no consume ni un número.
  Esa es la base de la prueba de compatibilidad de [D-16](#d-16).
- **Varios tractos en una molécula:** se procesan en orden canónico de inicio
  del tracto sobre la molécula, cada uno con su ámbito propio indexado por
  `locus_id`.

**Sobre los valores predeterminados: no inventarlos.** El proyecto tiene la
regla de no afirmar cifras sin comprobarlas, y ya pagó por incumplirla. No
existe en el repositorio ninguna calibración de tartamudeo. Propongo que los
valores predeterminados sean **cero** —tartamudeo desactivado— y que cualquier
valor que el usuario ponga se registre en el manifiesto como **no calibrado**,
hasta que exista una calibración aprobada. Si el usuario quiere un punto de
partida orientativo para sus pruebas, que lo elija él y quede constancia de que
es suyo.

**Qué se bloquea si no se decide.** El requisito no puede redactarse sin el
modelo, y la decisión técnica que fija los algoritmos tampoco.

---

<a id="d-15"></a>
### D-15. Qué pasa con la longitud del fragmento y el tamaño de inserto

Decisión nueva.

**Qué hay que decidir.** Si la molécula gana o pierde `k` unidades, cambia de
longitud en `k·(longitud del motivo)` bases. ¿Se deja cambiar, o se compensa?

**Opciones.**

| Opción | Consecuencia |
|---|---|
| (a) **La molécula cambia de longitud.** El inserto observado se mueve en `±k·(longitud del motivo)` | Es lo que ocurre físicamente: la molécula ya estaba seleccionada por tamaño cuando la polimerasa la copió, así que el amplicón sale más corto o más largo. Ensancha ligerísimamente la distribución de insertos, en la dirección correcta |
| (b) **Se compensa** volviendo a materializar más o menos secuencia flanqueante para conservar la longitud sorteada | Conserva exacta la distribución de insertos del perfil, pero es falso: ninguna polimerasa añade bases al otro extremo para compensar |
| (c) Se descarta la molécula cuando el tartamudeo la cambiaría de longitud | Sesga la cobertura del locus, que es el peor efecto posible |

**Recomendación: (a).** Con dos guardas explícitas:

1. Si la contracción dejara la molécula **más corta que la longitud de
   lectura**, ese paso no es viable y queda fuera del sorteo categórico
   ([D-14](#d-14)), porque `emit_paired_end` exige que el inserto sea al menos
   tan largo como la lectura. Con insertos de 250 a 450 pb y lecturas de 150,
   este caso es prácticamente imposible, pero la guarda tiene que estar: el
   error que produciría sería un fallo duro de la ejecución.
2. Las dos lecturas siguen midiendo exactamente la longitud del perfil. Lo que
   cambia es el trozo no observado del centro, salvo cuando las lecturas se
   solapan, que es justo el caso en que ambas ven el tracto.

**Qué se bloquea si no se decide.** La implementación y la prueba del tamaño de
inserto.

---

<a id="d-16"></a>
### D-16. Valores predeterminados y compatibilidad con lo ya publicado

Decisión nueva, y es la que hace el cambio aceptable en una fase cerrada.

**Qué hay que decidir.** La fase de realismo técnico está cerrada con sus
evidencias aprobadas. Añadirle un paso puede cambiar sus resultados. ¿Qué
garantía se da?

**Recomendación: la misma garantía que la fase de variación compleja se exige a
sí misma respecto de la población inicial.** Textualmente, adaptada:

> Con el tartamudeo desactivado —el valor predeterminado— una ejecución publica
> **exactamente los mismos bytes** que publicaba antes de que el tartamudeo
> existiera: los FASTQ finales, sus resúmenes, la verdad, el VCF y los
> manifiestos.

Esto es alcanzable, y no por casualidad: DTC-081 exige que cada sorteo arranque
de su propia clave y que «no existe consumo acumulado entre decisiones». Un
ámbito nuevo no puede desplazar los flujos existentes. Y con la tasa a cero,
`bernoulli-v1` no consume nada.

Es además la prueba más barata y más discriminante de todo el corte: un
`sha256sum` de cada artefacto antes y después.

**Qué se bloquea si no se decide.** La posibilidad de justificar el cambio en
una fase cerrada sin volver a ejecutar toda su evidencia de cierre.

---

## Grupo 6 — Cierre y gobierno

<a id="d-17"></a>
### D-17. ¿La interfaz Marimo es condición de cierre de este corte?

**Qué hay que decidir.** La hoja de ruta, sección 7, dice que una fase no se
cierra si una opción configurable carece de control en la interfaz y de
equivalencia interfaz → TOML → configuración efectiva. Este corte introduce
opciones configurables: los loci, los alelos y los parámetros del tartamudeo.

**Opciones.** (a) Sí, con su control y su prueba de ida y vuelta. (b) No,
porque el corte es un adelanto y no un cierre de fase; el control llega cuando
se cierre el hito de repeticiones de la fase de variación compleja.

**Recomendación: (a) para el tartamudeo, (b) diferible para los loci** si el
corte se ejecuta como adelanto y no como cierre de hito. La razón es
asimétrica: los parámetros del tartamudeo son los que el usuario va a querer
tocar una y otra vez para ver el efecto, y tocarlos a mano en un TOML es
exactamente la fricción que la interfaz existe para evitar. El catálogo de loci
se escribe una vez.

**Qué se bloquea si no se decide.** La condición de cierre del trabajo.

---

<a id="d-18"></a>
### D-18. ¿Esto es un hito propio o adelanta el hito de repeticiones?

**Qué hay que decidir.** El trabajo cae en dos fases distintas: las
repeticiones pertenecen a la fase de variación compleja, y el tartamudeo a la
fase de realismo técnico, que está cerrada. ¿Cómo se registra?

**Opciones.**

| Opción | Consecuencia |
|---|---|
| (a) **Un hito nuevo y propio**, transversal, con su plan y su pausa | Claro, y deja el plan de la fase de variación compleja intacto. Hay que actualizar la hoja de ruta |
| (b) **Adelantar los pasos F6.4.1 y F6.4.3** del plan vigente y añadir un paso nuevo a la fase de realismo técnico | Reutiliza numeración existente. Obliga a declarar que F6.4.3 pierde su dependencia del hito de variantes estructurales, que es precisamente el hallazgo de [D-02](#d-02) |
| (c) Esperar a la fase de variación compleja | No hay corte |

**Recomendación: (b), con la corrección de dependencia escrita.** Y, para el
tartamudeo, un paso nuevo en la fase de realismo técnico —numerado a
continuación de los existentes— que se ejecuta con la evidencia de cierre de
esa fase repetida. Reabrir una fase cerrada exige, según la propia hoja de
ruta, dejar constancia; ese es el sitio.

**Qué se bloquea si no se decide.** El registro del trabajo y su trazabilidad.

---

# Parte C — Borrador del requisito del tartamudeo

## 8. Antes del texto: qué hay que entender

**Hoy no existe nada.** He buscado en todo el repositorio, fuera del directorio
de compilación, las palabras «stutter», «tartamudeo», «slippage» y
«deslizamiento», en código, documentos, esquemas, configuraciones y la
interfaz: **cero apariciones**. No hay requisito funcional, ni decisión
técnica, ni ampliación registrada, ni código.

**Esto modifica una fase cerrada.** La fase de realismo técnico de
secuenciación está cerrada. Su objetivo declarado en la hoja de ruta es
«incorporar pronto los perfiles, degradaciones y eventos técnicos de WGS
confirmados, **sin alterar la verdad genotípica**», que es exactamente la
naturaleza del tartamudeo; encaja perfectamente en su propósito, pero no en su
contrato publicado. El cambio exige, como mínimo:

1. una decisión técnica nueva que fije el algoritmo y el catálogo de ámbitos
   —propongo numerarla a continuación de la última existente—;
2. los requisitos funcionales que siguen;
3. la modificación de la **frase de orden normativo** de la sección 4.4 de la
   especificación, que hoy dice: «F4 sigue este orden normativo: calcula las
   parejas de cobertura, sortea el fragmento biológico, forma una quimera
   cuando corresponda, realiza calidad y errores de R1/R2, aplica contaminación
   biológica y finalmente demultiplexa»;
4. la modificación de REQ-DAT-17, que enumera el mismo orden al describir los
   FASTQ finales;
5. una frase aclaratoria en REQ-TEC-07 o REQ-CONF-88 que deje dicho que este
   artefacto no es una preparación amplificada (ver [D-12](#d-12));
6. la repetición de la evidencia de cierre de la fase de realismo técnico.

**No es un añadido silencioso, y no debe presentarse como tal.**

## 9. Numeración propuesta

La familia REQ-DAT llega hoy a 73 y es donde viven las quimeras (REQ-DAT-60 a
64), que son el artefacto más parecido. La familia REQ-CONF llega a 130,
REQ-MAN a 92, REQ-VAL a 56 y los criterios de aceptación a 117.

Propongo, todos **como propuesta sujeta a aprobación**: REQ-DAT-74 a
REQ-DAT-80, REQ-CONF-131 a REQ-CONF-133, REQ-MAN-93, REQ-VAL-57 y ACE-118.

## 10. Texto propuesto

### 10.1. Párrafo explicativo, para la sección 4.4

Este párrafo va antes de los requisitos, junto a los que ya explican la lectura
corta y la calidad de una base.

> El **tartamudeo** es el cambio del número de unidades observado en una
> repetición en tándem, provocado por el deslizamiento de la polimerasa al
> copiar el tracto. Ocurre durante la amplificación en agrupación y la síntesis
> de la secuenciación, y también durante la amplificación por PCR cuando la
> hay. Es un artefacto técnico: la molécula copiada muestra un número de
> unidades distinto del que tiene el individuo, y la verdad genotípica no
> cambia. En datos reales, la frecuencia del tartamudeo crece con la longitud
> del tracto, es mayor cuanto más corto es el motivo y produce contracciones
> con más frecuencia que expansiones.

### 10.2. Requisitos propuestos

> - **REQ-DAT-74.** F4 debe poder aplicar tartamudeo a las repeticiones en
>   tándem **declaradas**. El tartamudeo cambia el número de unidades del
>   tracto contenido en una molécula, en unidades completas del motivo, y no
>   modifica en ningún caso el estado genotípico del individuo, la verdad
>   publicada ni el VCF de verdad. F4 no detecta tractos no declarados: una
>   región repetitiva que no corresponda a un locus de repetición en tándem
>   configurado no produce tartamudeo.
>
> - **REQ-DAT-75.** El tartamudeo es una **propiedad de la molécula, no de la
>   lectura**. Para cada molécula y cada tracto declarado contenido en ella, el
>   sistema debe tomar una única decisión y aplicarla a la molécula completa
>   antes de volver a formar sus dos lecturas. Cuando las dos lecturas de una
>   pareja cubran el mismo tracto, deben mostrar **el mismo** número de
>   unidades. El sistema no debe aplicar decisiones independientes a R1 y a R2.
>
> - **REQ-DAT-76.** El orden normativo de F4 pasa a ser: calcular las parejas
>   de cobertura, sortear el fragmento biológico, formar una quimera cuando
>   corresponda, **aplicar el tartamudeo cuando corresponda**, realizar calidad
>   y errores de R1/R2, aplicar contaminación biológica y demultiplexar. Cuando
>   una molécula sea quimérica, el tartamudeo se aplica a la molécula quimérica
>   ya formada.
>
> - **REQ-DAT-77.** La decisión de tartamudeo debe ser reproducible y derivarse
>   de un ámbito de semilla propio, identificado por el evento, la biblioteca,
>   la plantilla y el identificador del locus, y nunca del orden físico de
>   generación ni de un flujo aleatorio compartido. Debe constar de dos
>   decisiones: si la molécula tartamudea y, en caso afirmativo, cuántas
>   unidades gana o pierde. La probabilidad de tartamudear debe crecer con el
>   número de unidades contenidas en la molécula y depender de la longitud del
>   motivo; la elección entre expansión y contracción debe ser configurable. Los
>   pasos que no dejarían al menos una unidad completa en el tracto, o que
>   dejarían la molécula más corta que la longitud de lectura, deben excluirse
>   **antes** de sortear; el sistema no debe reintentar. Si ningún paso es
>   viable, la molécula no tartamudea.
>
> - **REQ-DAT-78.** Cuando una molécula gana o pierde unidades, **su longitud
>   cambia** en el número de unidades multiplicado por la longitud del motivo, y
>   con ella el tamaño de inserto observado de la pareja. El sistema no debe
>   compensar ese cambio materializando secuencia flanqueante adicional ni
>   recortando la molécula para conservar la longitud sorteada. Las dos lecturas
>   conservan la longitud del perfil.
>
> - **REQ-DAT-79.** Una molécula que contiene solo parte de un tracto puede
>   tartamudear si contiene al menos el número mínimo configurable de unidades
>   completas; por debajo de ese mínimo no tartamudea y no consume ninguna
>   decisión aleatoria. Una lectura que cubre parcialmente el tracto, o que lo
>   toca por un extremo, muestra la parte del tracto ya modificado que le
>   corresponda: no se aplica ninguna corrección por lectura.
>
> - **REQ-DAT-80.** El tartamudeo **no altera la verdad**. Dos ejecuciones
>   idénticas salvo en la tasa de tartamudeo deben publicar la verdad, el VCF de
>   verdad y los manifiestos de verdad byte a byte iguales, y diferir únicamente
>   en los FASTQ finales y en los recuentos técnicos agregados. Con el
>   tartamudeo desactivado, todos los artefactos publicados deben ser byte a
>   byte iguales a los de una ejecución equivalente anterior a la existencia del
>   tartamudeo.
>
> - **REQ-CONF-131.** La configuración del evento técnico debe permitir activar
>   el tartamudeo y definir, por longitud de motivo, los parámetros de su
>   probabilidad y el número mínimo de unidades contenidas que lo hacen
>   posible; y definir el conjunto de pasos admitidos en unidades enteras con
>   sus pesos relativos. Debe permitir sustituciones por locus.
>
> - **REQ-CONF-132.** El tartamudeo debe estar **desactivado de forma
>   predeterminada**. Una configuración que seleccione un modelo distinto del
>   admitido, un paso expresado en bases en vez de unidades o un parámetro fuera
>   de su dominio debe rechazarse antes de simular, identificando el parámetro y
>   el valor.
>
> - **REQ-CONF-133.** El tartamudeo no constituye una preparación de biblioteca
>   amplificada: no introduce duplicados de PCR ni sesgos de amplificación, y el
>   perfil de preparación sigue registrándose como sin PCR conforme a
>   REQ-TEC-07. La preparación amplificada y su tasa de tartamudeo propia siguen
>   aplazadas en DTF-056.
>
> - **REQ-MAN-93.** El manifiesto del evento debe registrar la activación del
>   tartamudeo, los parámetros efectivos por longitud de motivo y por locus, la
>   versión del algoritmo, la versión del catálogo de ámbitos aleatorios y, como
>   **recuentos agregados**, el número de moléculas evaluadas, el número que
>   tartamudeó y la distribución de los pasos aplicados por locus. No debe
>   registrar procedencia por lectura ni por pareja, conforme a REQ-DAT-56. Los
>   parámetros deben registrarse indicando que no proceden de una calibración
>   empírica aprobada.
>
> - **REQ-VAL-57.** La validación interna debe poder reproducir las decisiones
>   de tartamudeo desde la semilla de la ejecución, sin reutilizar el código
>   productor, y comprobar sobre los FASTQ finales que las dos lecturas de cada
>   pareja que cubren un mismo tracto muestran el mismo número de unidades.

### 10.3. Criterio de aceptación propuesto

> | **ACE-118. Tartamudeo de repeticiones en tándem** | El tartamudeo se aplica por molécula y tracto declarado, nunca por lectura: las dos lecturas de una pareja que cubren el mismo tracto muestran el mismo número de unidades. Está desactivado de forma predeterminada y, desactivado, la ejecución publica los mismos bytes que sin la capacidad. Su probabilidad crece con el número de unidades contenidas y depende de la longitud del motivo; las contracciones y expansiones se ponderan por configuración. La molécula cambia de longitud y con ella el inserto observado. La verdad, el VCF y sus manifiestos no cambian con la tasa. Cada decisión es reproducible desde un ámbito de semilla propio | Vectores congelados, prueba de coherencia del par con tracto en el solapamiento, identidad byte a byte con tasa cero y comparación de verdad entre tasas |

## 11. Pruebas, y contra qué variante defectuosa tiene que fallar cada una

El defecto recurrente número uno de este proyecto son pruebas que no pueden
fallar por el defecto que dicen cubrir. Por eso cada comprobación lleva aquí su
variante defectuosa concreta.

| # | Qué comprueba | Contra qué variante defectuosa debe fallar | Cuidado |
|---|---|---|---|
| 1 | **Coherencia del par.** Fixture con el tracto entero dentro del solapamiento de R1 y R2. Con el perfil vigente —inserto entre 250 y 450 pb, media 350, lectura 150— el solapamiento existe solo en el extremo corto de la distribución: con inserto 250 las dos lecturas comparten las 50 bases centrales. El fixture usa una copia haplotípica lo bastante corta para forzar el inserto mínimo y coloca un tracto de 30 pb dentro de esas 50 bases. Cuenta unidades en R1 y en el complemento inverso de R2 y exige que sean iguales | Un tartamudeo aplicado **por lectura** | Un fixture en que solo R1 cubre el tracto **no discrimina**: la aserción se cumple trivialmente, y con la media de 350 la mayoría de las parejas están en ese caso. Hace falta el fixture **saturado**, con inserto en el mínimo y el tracto íntegramente dentro del solapamiento. Decirlo en la *docstring* |
| 2 | **Tasa cero, mismos bytes.** `sha256sum` de los dos FASTQ, la verdad, el VCF y los manifiestos, con la capacidad presente y desactivada, contra los de la ejecución equivalente sin la capacidad | Una implementación que consuma aleatoriedad con la tasa a cero, o que desplace el ordinal de plantilla | Exige guardar los hashes de referencia **fuera de git** antes de tocar el código |
| 3 | **La verdad no cambia.** Misma ejecución con tasa 0 y con tasa alta: verdad, VCF y manifiestos de verdad con el mismo `sha256`; FASTQ distintos | Una implementación que aplique el tartamudeo a la **copia haplotípica** en vez de a la molécula, es decir, que cambie el genotipo | Si los FASTQ salieran iguales, la prueba pasaría por el motivo equivocado: hay que exigir además que **difieran** |
| 4 | **Dirección.** Peso 1 a la contracción y 0 a la expansión: todas las moléculas que tartamudean salen más cortas | Un signo invertido | El caso con peso 1 no consume aleatoriedad en el categórico; conviene un segundo caso con pesos repartidos |
| 5 | **Longitud del tracto.** Dos loci con el mismo motivo, uno de 6 unidades y otro de 30: la probabilidad calculada por la función del modelo es estrictamente mayor para el largo, y sobre un número fijo de moléculas con semilla congelada la fracción observada también | Un modelo con probabilidad constante que ignore la longitud | Probar **la función de probabilidad** directamente, además de la fracción muestral; una prueba solo muestral puede no discriminar con pocas moléculas |
| 6 | **Longitud del motivo.** Un mononucleótido y un tetranucleótido con la **misma longitud en pares de base**: probabilidad mayor para el mononucleótido | Un modelo que dependa solo de los pares de base del tracto y no de la longitud del motivo | El fixture tiene que igualar la longitud en bases, si no la diferencia podría venir del número de unidades |
| 7 | **Inserto.** Una molécula que pierde una unidad de un motivo de 4 pb produce un fragmento 4 pb más corto | Una implementación que recorte o rellene para conservar la longitud sorteada | Hay que observarlo sobre la molécula, no sobre las lecturas, que siguen midiendo lo mismo |
| 8 | **Cobertura parcial.** Una molécula que contiene menos unidades completas que el mínimo no tartamudea y **no consume aleatoriedad**, contado con un generador instrumentado | Una implementación que tartamudee sobre un resto de una unidad, o que sortee y descarte | Contar el consumo, no solo el resultado: una implementación que sortee y descarte da el mismo resultado y desplaza el flujo |
| 9 | **Quimera.** Corpus con tasa de quimera y tasa de tartamudeo altas a la vez: ninguna ejecución produce `ChimeraError::PairSourceMismatch` | Una implementación que modifique el fragmento **antes** de la quimera | Sin tasa de quimera alta la prueba no ejercita el camino |
| 10 | **Reproducibilidad.** Misma semilla, mismo `sha256` de los FASTQ; otra semilla, distinto | Un ámbito de semilla que dependa del orden físico o de un flujo compartido | Comprobar también que **añadir un locus en otro cromosoma no cambia** las decisiones de los loci previos: discrimina un ámbito indexado por posición en vez de por identificador del locus |
| 11 | **Validador independiente.** El validador reproduce las decisiones desde la semilla y regenera las lecturas sin llamar al productor | Un validador que lea los bytes publicados por el productor en vez de recomponerlos | Es la comprobación que el proyecto exige con `./scripts/check-independent-validators.sh` |

Para todas las que exijan mutar el código: copia de seguridad **fuera de git**,
comprobación por `SHA` de que la mutación mutó y de que la restauración devolvió
el fichero byte a byte. Nunca `git checkout --`.

---

## 12. Resumen de las modificaciones propuestas a otros documentos

Ninguna está aplicada. Se enumeran con su ubicación exacta.

| Documento | Dónde | Qué |
|---|---|---|
| `design/especificaciones-0.1.md` | §4.3, párrafo del glosario de la repetición en tándem | Añadir la regla de interrupciones declaradas y pureza realizada ([D-06](#d-06)) |
| `design/especificaciones-0.1.md` | §4.4, frase «F4 sigue este orden normativo…» | Insertar el tartamudeo entre la quimera y la calidad |
| `design/especificaciones-0.1.md` | REQ-DAT-17 | Mismo cambio de orden |
| `design/especificaciones-0.1.md` | §4.4, tras REQ-DAT-73 | Añadir REQ-DAT-74 a REQ-DAT-80 |
| `design/especificaciones-0.1.md` | §4.1, tras REQ-CONF-130 | Añadir REQ-CONF-131 a REQ-CONF-133 |
| `design/especificaciones-0.1.md` | REQ-TEC-07 o REQ-CONF-88 | Frase aclaratoria de [D-12](#d-12) |
| `design/especificaciones-0.1.md` | §4.7 y §4.6 | Añadir REQ-MAN-93 y REQ-VAL-57 |
| `design/especificaciones-0.1.md` | §7 | Añadir ACE-118 |
| `design/decisiones-tecnicas-0.1.md` | Tabla de decisiones | Decisión nueva con el algoritmo `stutter-v1`, el catálogo `f4-technical-v4` y la confirmación de que DTC-067 no prohíbe usar `bernoulli-v1` y `categorical-f64-v1` aquí ([D-11](#d-11), [D-14](#d-14)) |
| `design/decisiones-tecnicas-0.1.md` | Tabla de decisiones | Decisión nueva con la realización de interrupciones y pureza ([D-06](#d-06)) |
| `design/plan-implementacion-f6-0.1.md` | §8, paso F6.4.3 | Corregir la dependencia: ya no depende del primer paso del hito de variantes estructurales ([D-02](#d-02)) |
| `design/plan-implementacion-f6-0.1.md` | §8, paso F6.4.2 | Corregir la dependencia: no depende de derivar una referencia estructural ([D-03](#d-03)) |
| `design/plan-implementacion-f6-0.1.md` | §3, tabla de decisiones pendientes | Recortar DP-F6-03, DP-F6-04 y DP-F6-06 a lo que este corte necesita |
| `design/hoja-de-ruta-0.1.md` | §F4 y §F6 | Registrar el adelanto y la reapertura del contrato de la fase de realismo técnico ([D-18](#d-18)) |
| `design/registro-ampliaciones-futuras-0.1.md` | DTF-056 | Anotar que la tasa de tartamudeo del perfil amplificado pertenece a esa ampliación |
| `design/trazabilidad-requisitos-0.1.md` | Matriz | Asignar los requisitos nuevos |

---

## 13. Orden de trabajo recomendado

1. **Decidir [D-12](#d-12)**, el contrato de la fase cerrada. Si la respuesta es
   que no cabe, el resto de la parte C se cae y hay que replantear el encargo.
2. **Decidir [D-02](#d-02) y [D-06](#d-06).** Con esas dos se puede empezar el
   modelo de dominio del locus, que no depende de nada más.
3. **Implementar el modelo de dominio del locus.** Es puro, aislado y
   verificable contra secuencias escritas a mano.
4. **Decidir [D-01](#d-01), [D-03](#d-03), [D-04](#d-04), [D-05](#d-05) y
   [D-07](#d-07)**, que son el contrato de entrada.
5. **Implementar el alelo por copia haplotípica** y probar la materialización.
6. **Decidir [D-08](#d-08) y [D-09](#d-09)** y publicar verdad y VCF.
7. **Decidir [D-11](#d-11), [D-13](#d-13), [D-14](#d-14), [D-15](#d-15) y
   [D-16](#d-16)** e implementar el tartamudeo.
8. **Decidir [D-17](#d-17) y [D-18](#d-18)** y cerrar.

El paso 3 puede empezar hoy mismo en cuanto estén [D-02](#d-02) y
[D-06](#d-06), y el paso 7 es independiente del 5 y del 6 salvo por el
identificador del locus: si hace falta, se puede trabajar en paralelo.
