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

**Estado «Aceptada».** La deficiencia sigue existiendo y **sigue abierta**: es
conocida, el responsable ha decidido convivir con ella y no hay trabajo
previsto. Se distingue de «Cerrada» porque nada la ha arreglado, y de «Abierta»
porque no se espera arreglarla. Cada entrada aceptada dice **qué** se paga por
convivir con ella y **qué condición** obligaría a reconsiderarla.

---

## 2. Registro

| ID | Deuda | Gravedad | Detectada | Estado |
|---|---|---|---|---|
| DT-001 | Dos tercios de las comprobaciones de los verificadores independientes no las ejercita ninguna prueba | **Grave** | 21/09/2026 | Abierta, muy reducida (29/09/2026); punto ciego estructural cerrado (04/10/2026); nueve jueces más vigilados (05/10/2026); tanda 1 de cobertura, jueces pequeños (08/10/2026); tanda 2, jueces medianos (08/10/2026); tanda 3, matriz de SNV (08/10/2026); tanda 4, manifiesto de evento (08/10/2026); tanda 5, verdad (09/10/2026); tanda 6, primera parte de `results.rs`: microsatélites, tartamudeo y PCR (10/10/2026); tanda 7, segunda parte de `results.rs`: la configuración efectiva (11/10/2026); tanda 8, tercera parte de `results.rs`: estadísticas iniciales de las poblaciones y auditoría de la matriz (12/10/2026) |
| DT-002 | Dos pruebas ignoradas están vaciadas, no solo apagadas | Media | 21/09/2026 | **Cerrada (06/10/2026)**, borradas; los dos huecos que dejaban, cubiertos con pruebas nuevas de la vía de parejas |
| DT-003 | Ocho pruebas ignoradas corresponden a contratos que ya no existen | Baja | 20/09/2026 | Abierta |
| DT-004 | La puerta de Python no corre desde el 25/08/2026 | Media | 19/09/2026 | **Cerrada (06/10/2026)** |
| DT-005 | La marca del lector de verdad que distingue un tramo sustituido por un tracto no la ejercita ninguna prueba | Baja | 22/09/2026 | Abierta |
| DT-006 | El sorteo del inserto por rechazo es inutilizable cuando la copia más larga acerca el máximo al mínimo del perfil | Media | 24/09/2026 | **Aceptada (06/10/2026)**, con guarda de admisión (25/09/2026) |
| DT-007 | Ningún diagnóstico público nuevo puede publicarse en un manifiesto del recorrido general | Baja | 25/09/2026 | Abierta |
| DT-008 | El orden entre el paso de linaje y la decisión de base no lo separaba ninguna prueba | Baja | 26/09/2026 | **Cerrada (26/09/2026)** |
| DT-009 | Las dos fuentes de la calibración discrepan en cuánto multiplica la PCR el tartamudeo de los homopolímeros | Media | 27/09/2026 | **Cerrada por decisión (06/10/2026)**; implementada el 08/10/2026 (DTC-099) |
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
| DT-023 | Resolver una referencia desde un manifiesto juzga la verdad entera de aquel paquete | Baja | 04/10/2026 | **Aceptada (04/10/2026)** |
| DT-024 | `validate_published_resources` no termina con `lineages` = 18446744073709551615 en la configuración efectiva | Media | 04/10/2026 | **Cerrada (06/10/2026)** |
| DT-025 | `validate_published_resources` aborta con `capacity overflow` si `reference_units[1]` = 18446744073709551615 en la efectiva | Media | 04/10/2026 | **Cerrada (06/10/2026)** |
| DT-026 | El juez de cada mate del FASTQ final no leía ni un registro: toda la validación de registros era código muerto en producción | **Grave** | 05/10/2026 | **Cerrada (05/10/2026)** |
| DT-027 | La verdad de un paquete sin poblaciones no se juzga: 31 alteraciones de la cabecera pasan el juez de los recursos publicados | Media | 05/10/2026 | **Cerrada (06/10/2026)** |
| DT-028 | Sospecha: un alelo de referencia con una clave desconocida se acepta, el alelo SNP las rechaza | Baja | 05/10/2026 | **Cerrada (06/10/2026)**, reproducida |
| DT-029 | El juez del FASTQ final no comprueba siete propiedades del contrato: identidad del fichero, unicidad de la pareja, formato de los identificadores, longitud, `N` con Phred 0, R1 frente a R2 y resúmenes del manifiesto | Media | 05/10/2026 | **Cerrada (05/10/2026)** |
| DT-030 | Un bloque BGZF que declara más de 64 KiB descomprimidos hace abortar al juez con un pánico en vez de rechazar el fichero | Media | 05/10/2026 | **Cerrada (06/10/2026)**, sin cambiar dependencias |
| DT-031 | Una referencia corta con tasa de quimeras alta aborta como «invariante imposible» en vez de rechazarse como configuración | Baja | 05/10/2026 | **Cerrada (06/10/2026)** |
| DT-032 | El lector de la verdad acepta como lista JSON objetos que el esquema exige objeto (`"counts":[6,8]`) | Baja | 06/10/2026 | **Cerrada (06/10/2026)** |
| DT-033 | La verificación del paquete F5 publicado acepta manifiestos con enlaces escritos como lista JSON, que el esquema rechaza | Baja | 06/10/2026 | **Cerrada (06/10/2026)** |
| DT-034 | Lo publicado atribuye las cifras de calibración del tartamudeo solo a Raz et al. (2019), cuando son de Gymrek (2016), cuya licencia exige atribuirlas | Baja | 06/10/2026 | Abierta |
| DT-035 | La puerta del paquete F5 acepta manifiestos que el esquema rechaza en campos que no lee con forma tipada (`attempt_summary`, `quota_ties`, `limits`, `profile`) | Baja | 06/10/2026 | Abierta |
| DT-036 | Un ensamblado de evento fallido no deja el evento enlazado con su estado (REQ-SAL-21); el campo `link` nunca lleva valor | Media | 06/10/2026 | Abierta |
| DT-037 | `tests/test_population_snvs_editor.py` está fuera de la puerta de Python y falla 1 de 7 | Baja | 06/10/2026 | Abierta |
| DT-038 | Ante un fallo posterior al evento, el VCF ya publicado y validado se descarta del paquete fallido aunque es independiente | Media | 06/10/2026 | Abierta |
| DT-039 | Con evento, si se altera el FASTA de la referencia suministrada antes del cierre, la rama de fallo controlado no consigue publicar y se pierde la causa original | Media | 06/10/2026 | Abierta |
| DT-040 | Con tartamudeo clonal, la configuración efectiva publica `technical_event.rng_catalog = "f4-technical-v5"` y su esquema v6 exige `f4-technical-v6`: el validador externo la rechaza | Media | 08/10/2026 | Abierta; decidido el 08/10/2026 que manda el esquema, en un cambio aparte junto con DT-041 |
| DT-041 | Toda ejecución amplificada, con microsatélites o sin ellos, publica un manifiesto general que no pasa su propio esquema: el perfil está fijado al de sin PCR y los catálogos a v3 y v4 | Media | 08/10/2026 | Abierta; decidido que manda el contrato, ampliándolo; junto con DT-040 |
| DT-042 | En el cuaderno, con un evento técnico activo, la celda que reconstruye el formulario falla con `NameError` sobre `TechnicalEventForm`: ningún control del panel de PCR llega al TOML | Media | 08/10/2026 | Abierta |
| DT-043 | El juez acepta en el manifiesto de evento los recuentos realizados de la amplificación (familias, parejas, histograma, linaje) alterados, con las huellas rehechas | Media | 08/10/2026 | Abierta, reducida (10/10/2026): los parámetros de la amplificación ya se atan a la efectiva; quedan los recuentos realizados |
| DT-044 | REQ-CONF-134 a 138, REQ-DAT-81 a 85 y REQ-MAN-94 no tienen fila en la matriz de trazabilidad | Baja | 08/10/2026 | Abierta |
| DT-046 | Referencias cruzadas internas del manifiesto de evento no atadas por su propio juez (`id_catalog`, bibliotecas repetidas, receptora y donantes ajenos, UDI repetidos, esquema v2 en los hechos del escritor) | Baja | 08/10/2026 | Abierta |
| DT-045 | El juez de la matriz de SNV tiene catorce comprobaciones de código muerto y un conjunto redundante, porque repite lo que `parse_metadata_json` ya exige | Baja | 08/10/2026 | Abierta; encargo aparte |
| DT-047 | El contrato de la verdad admite un cambio de base de segmento cubierto por una deleción del mismo segmento y la réplica del juez lo rechaza como «no se reproduce» | Baja | 09/10/2026 | **Cerrada (09/10/2026)**: el error estaba en el contrato |
| DT-048 | El juez no exige que el estado heredado de una copia recorra la referencia una vez y en orden creciente (colinealidad): acepta rotaciones y duplicaciones | Media | 09/10/2026 | **Cerrada (09/10/2026)**: el juez exige que los segmentos recorran la referencia una vez y en orden |
| DT-049 | El lector por flujo de la verdad no aplica a los casos genotípicos el paso léxico del contrato | Media-baja | 09/10/2026 | Abierta (junto a DT-015) |
| DT-050 | En memoria, el juez de la verdad aceptaba una copia con un segmento de otra secuencia en medio | Baja | 09/10/2026 | **Cerrada (09/10/2026)**, en el mismo cambio que DT-048 |
| DT-051 | El juez de los recursos publicados no ata los loci de repetición publicados (manifiesto y configuración efectiva) al TOML copiado: acepta alterar secuencias, frecuencias, pureza, intervalos, topes y mínimos | Media | 10/10/2026 | **Cerrada (10/10/2026)**, en un cambio propio tras la tanda 6 |
| DT-052 | La contabilidad de bases de tracto no la juzga nadie al publicar: su juez no está cableado y, aunque lo estuviera, no ata las cotas al TOML | Media | 10/10/2026 | **Cerrada (10/10/2026)**, en un cambio propio tras DT-051 |
| DT-053 | Las copias realizadas de cada alelo no se atan al sorteo rehecho con la semilla | Baja | 10/10/2026 | Abierta; decidido no arreglarla por ahora |
| DT-054 | Nadie compara el `execution_id` de la cabecera de la verdad y del VCF con el del manifiesto general: un TOML copiado cambiado de forma coherente pasaba con la verdad de la ejecución original | Media | 11/10/2026 | **Cerrada (11/10/2026)**, en el mismo cambio que DT-052 |
| DT-055 | El juez de la pareja efectiva/manifiesto no comprueba dos desigualdades de las estadísticas iniciales que la arquitectura fija (los genotipos realizados cuentan solo observaciones completas); las caza después el juez completo | Baja | 12/10/2026 | Abierta |

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

**Y 499 y 147 eran un mínimo, no el total**, porque quedaba la segunda
familia, que es estructural, y que **se cerró el 04/10/2026**:

- los validadores auxiliares que devuelven un motivo como texto en vez de un
  diagnóstico (quince funciones de `src/results.rs`: catorce con
  `Result<_, &'static str>` y una, el rehacedor del abanico de deslizamiento,
  con `Result<_, String>`, con ciento setenta rechazos entre `Err("texto")` y
  `ok_or("texto")`): su texto lo convertía, en el llamante, **una sola
  llamada** a `effective_error`, `manifest_error` o `repeat_section_error`, así
  que decenas de guardas contaban como un único punto, y bastaba con que se
  disparase uno para darlas por cubiertas. Había además un cierre `reject` y
  `repeat_section_error` sin `#[track_caller]`, con el mismo efecto sobre
  veinticinco rechazos del juez del abanico y nueve llamadores de
  `repeat_section_error` (ocho en `validate_published_tract_output_accounting`
  y uno en el lector del manifiesto).

**Cómo se cerró.** Cada rechazo de texto se escribe ahora
`Err(reason("..."))` (o `ok_or_else(|| reason("..."))`, nunca `ok_or(reason(..))`,
que anotaría siempre). `reason` es una constructora genérica con
`#[track_caller]` que anota la línea y devuelve el texto tal cual, de modo que
**las firmas, los llamantes, los mensajes, los códigos y el orden de las
comprobaciones no cambian**: se comprobó que la secuencia de literales de
`src/results.rs` es idéntica salvo cuatro nombres de esquema nuevos que
repiten los que ya usaban los envoltorios. Los envoltorios que reciben texto ya
anotado construyen el diagnóstico sin anotar (`manifest_diagnostic` y
`effective_diagnostic`) para no contar como rechazo una línea que no decide
nada. Dos pruebas léxicas nuevas en `tests/guard_constructors.rs` impiden que
vuelva a escribirse un `Err("literal")` sin `reason` y comprueban que su
instrumento lo ve.

**Cifras tras el cierre**, medidas con la puerta sobre el árbol final:

| Módulo vigilado | Comprobaciones antes | Sin disparar antes | Comprobaciones hoy | Sin disparar hoy |
|---|---|---|---|---|
| `src/results.rs` | 541 | 423 | **739** | **586** |

Los otros cuatro módulos no se mueven. **Suben las dos cifras y no es una
regresión**: es lo que antes estaba oculto. Ejemplo medido: el validador de la
sección `repeat_loci` del manifiesto tenía 24 rechazos y **ningún** punto
propio —todos eran una línea que disparaba cualquier otra prueba—; ahora tiene
24, y **17 no los dispara ninguna prueba**. El techo de la línea base pasa a la
cifra medida, 586; cubrirlas sigue siendo el trabajo de esta deuda.

**Corregido el mismo día tras la revisión adversarial: 773 y 619.** La
revisión encontró dos puntos ciegos más del mismo tipo en `src/results.rs`:
`physical_error` agrupaba 33 llamadores en una sola llamada a `effective_error`
(ahora lleva `#[track_caller]` y está vigilado), y dos `GPS-IO-002` de
identidad de artefacto se construían a mano (ahora pasan por
`artifact_identity_error`, vigilada). Las comprobaciones pasan de 739 a 773 y
las sin disparar de 586 a 619: de las 34 netas, 33 no las ejercita ninguna
prueba. Además se cerraron tres vías por las que un rechazo nuevo podía
esquivar la puerta, todas con una prueba léxica que se verificó **muriendo**
contra un mutante: construir un rechazo con `effective_diagnostic` o
`manifest_diagnostic` fuera de sus envoltorios, llamar a `*_from_reason` con un
literal, y escribir `ok_or(reason(..))`, que anota siempre.

**Lo que sigue sin verse.** Los ayudantes que reciben el motivo como argumento
y hacen un único rechazo para muchos llamantes (`required_non_empty_text` y
similares) se cuentan una vez: son comprobaciones genéricas, no texto suelto.

**Nueve módulos más entran en la puerta (05/10/2026).** Los que juzgan lo
publicado o lo leído y que este registro listaba sin vigilar pasan a la puerta,
por decisión del responsable. Cifras **medidas**, no estimadas (la estimación
léxica anterior quedaba corta en casi todos), con el mismo número de pruebas y
los mismos mensajes, variantes y códigos:

| Módulo | Qué juzga | Comprobaciones | Sin disparar |
|---|---|---|---|
| `src/truth_validation.rs` | el juez semántico de la verdad | 118 | 95 |
| `src/truth_schema.rs` | el contrato de la verdad | 69 | 46 |
| `src/fastq_validation.rs` | el FASTQ final y su repetición | 71 | 57 |
| `src/population_snvs.rs` | la matriz de SNV y la auditoría de normalización | 143 | 123 |
| `src/f5_publication.rs` | `verify_published_f5_package` y los documentos que publica | 52 | 36 |
| `src/reference_fasta.rs` | el FASTA suministrado | 41 | 30 |
| `src/reference_reduction_io.rs` | el mapa de reducción | 55 | 45 |
| `src/integrity.rs` | la identidad de un artefacto | 1 | 0 |
| `src/technical_profile.rs` | el recurso del perfil técnico | 15 | 14 |
| **Total** | | **565** | **446** |

**Cómo se hizo.** Una constructora genérica nueva, `rejected` en
`src/guard_trace.rs` (`#[track_caller]`, anota la línea y devuelve el error
tal cual), envuelve cada rechazo: `Err(rejected(Variante { .. }))`,
`.ok_or_else(|| rejected(..))` y `.map_err(|_| rejected(..))`; nunca
`ok_or(rejected(..))`, que anotaría siempre. Los módulos que ya tenían
funciones constructoras (`mismatch`, `audit_error`, `invalid`, `missing`,
`document_error`, `manifest_error`, `empty_header`, `invalid_symbol`,
`map_error`, `configuration_mismatch` y otras) las anotan ellas mismas con
`#[track_caller]` y `record`, y se vigilan por nombre en el array `judges`.
Dos clausuras que escondían varios rechazos tras una sola línea
(`out_of_range` en `reference_reduction_io.rs` y `mismatch` en
`f5_publication.rs`) pasaron a ser funciones; en `fastq_validation.rs` y
`truth_validation.rs` las clausuras `structural` y `missing` siguen siendo
clausuras (no admiten `#[track_caller]`) y cada uso que decide un rechazo va
envuelto en `rejected(..)`; la excepción es `.map_err(structural)` de
`parse_final_header` en `fastq_validation.rs`, que reetiqueta un texto ya
anotado. Se comprobó que la secuencia de literales de cada módulo es la
misma, salvo la reordenación de los dos textos de las funciones nuevas. La
prueba `no_rejection_of_the_added_judges_skips_its_constructor`, en
`tests/guard_constructors.rs`, impide que un rechazo nuevo se salte la
constructora, y se verificó **muriendo** contra tres mutantes.

**Lo que queda fuera de estos módulos, y por qué.**

- `src/f5_artifacts.rs` **no entra**: su propia cabecera dice «este módulo no
  decide nada: escribe lo que F5.1, F5.2 y F5.3 ya cerraron». Sus `internal` y
  `configuration` son invariantes internas del escritor, no un juicio sobre lo
  publicado. Lo que sí juzga esos ficheros ya lo cubren `src/f5_oracle.rs`,
  `src/f5_closure.rs` y `src/results.rs`.
- Los fallos de **entrada y salida** (`input_error`, `io_error`, `Diagnostic::io`
  y similares) y las invariantes **internas** (`internal_error`, `size_error`)
  no se cuentan en `reference_fasta.rs`, `reference_reduction_io.rs`,
  `integrity.rs` ni `fastq_validation.rs`: no son un veredicto sobre el
  contenido que se juzga. **Una excepción en `fastq_validation.rs`:**
  `validate_final_artifact` informa de un fallo de `identify_file` como
  `Artifact` (campo `sha256`), y ese rechazo sí se cuenta. En los módulos de enumeración de errores sí se
  cuenta toda variante (también los desbordes de contador), por sencillez; son
  guardas casi inalcanzables y engordan la cifra de sin disparar.
- La validación de la configuración del **generador** de matrices de
  `population_snvs.rs` (`PopulationSnvsGenerationError`): juzga la entrada del
  usuario antes de generar, no lo publicado. La prueba léxica la exime por su
  nombre.
- Los ayudantes que reciben el motivo como argumento y hacen un único
  rechazo para muchos llamantes se cuentan **una vez**: `require` en el juez de
  casos de `truth_validation.rs` (seis llamantes); las clausuras `text` y
  `number` de la auditoría de normalización de `population_snvs.rs`; y
  `required_object`, `required_string` y `required_literal`, de la lectura de
  los metadatos de la misma matriz.
- `src/vcf_contract.rs` (unas 20 variantes de `VcfContractError`) es sobre
  todo un escritor, no un juez, y sigue fuera.
- Los rechazos de fallo de fichero o de otro módulo que estos jueces
  propagan con `?` o `map_err(Variante::Envoltorio)` no son puntos propios:
  los cuenta el módulo que los origina, si está vigilado.

**Regla de los envoltorios, y lo que la prueba léxica no ve.** Una función que
solo reetiqueta un error que ya anotó su origen (`parse_final_header`,
`parse_metadata_json`, `parse_canonical_tsv`, `reconstruct_copy`,
`validate_genotype_semantics`, `validate_case_set`,
`validate_reference_digest`) **no** lleva `rejected(..)`: contaría como punto un
rechazo que no decide nada. Siete lo llevaban al entrar y se quitaron tras la
revisión adversarial, como `reference_content_sha256`, que usan también los
productores y cuyo fallo el juez de recursos publicados se traga a propósito;
sus llamadas que sí rechazan anotan ellas. La prueba léxica admite estos
reetiquetados por **prefijo con recuento exacto** (cada exención declara
cuántos usos tiene y la prueba falla si cambia). Lo que **no** ve, y solo caza
la puerta por recuento: un rechazo escrito con una variable intermedia,
`<D::Error as serde::de::Error>::custom(..)`, y un rechazo nuevo escrito con un
prefijo exento si a la vez se sube el recuento declarado. Un `rejected(..)`
puesto en un camino que se ejecuta siempre, sin rama que rechace, tampoco lo
ve nadie: la puerta lo cuenta como disparado y solo deja un AVISO por
comprobaciones nuevas.

**Estas cifras no mejoran nada por sí solas.** Es lo mismo que pasó al meter
el manifiesto general y `results.rs`: se hace visible lo que no se veía. Cubrir
las 446 sin disparar es el trabajo de esta deuda; el techo de cada módulo es la
cifra medida.

La lista autoritativa de lo que la puerta vigila es el array `judges` de
`scripts/check-guard-coverage.sh`.

En `src/event_manifest.rs` no hay este segundo problema (no tiene ningún
`Err("...")` a pelo), y su única construcción de diagnóstico fuera de la
constructora es `retag_event_manifest_schema`, que reetiqueta un error ya
rechazado. Una prueba léxica (`tests/guard_constructors.rs`) prohíbe construir
un `GPS-FMT-006` fuera de las constructoras de ambos módulos y, desde el
04/10/2026, limita dónde pueden llamarse los constructores sin anotar
(`effective_diagnostic`, `manifest_diagnostic` y los `*_from_reason`); lo que
**no** puede ver es un rechazo nuevo escrito dentro de una de las funciones
permitidas.

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

### Cobertura por tandas (desde el 08/10/2026)

Decisión del responsable (08/10/2026): reducir las sin disparar **por jueces,
empezando por los pequeños**. Cada tanda cubre con pruebas que exigen el motivo
exacto del rechazo, demuestra con razonamiento **y** un intento ejecutado lo
que es inalcanzable (anotado junto a la guarda: «Inalcanzable porque…») y
verifica la discriminación **mutando cada guarda en una copia**, con huella
SHA antes, mutada y restaurada. Una línea por tanda:

| Tanda | Fecha | Jueces | Antes (comprobaciones / sin disparar) | Después |
|---|---|---|---|---|
| 1 | 08/10/2026 | `src/technical_profile.rs` | 15 / 14 | 15 / **5** |
| 1 | 08/10/2026 | `src/reference_fasta.rs` | 41 / 30 | 41 / **7** |
| 1 | 08/10/2026 | `src/reference_reduction_io.rs` | 55 / 45 | 55 / **12** |
| 2 | 08/10/2026 | `src/fastq_validation.rs` | 84 / 36 | 84 / **26** |
| 2 | 08/10/2026 | `src/f5_publication.rs` | 53 / 36 | 53 / **10** |
| 3 | 08/10/2026 | `src/population_snvs.rs` | 143 / 123 | 143 / **25** |
| 4 | 08/10/2026 | `src/event_manifest.rs` | 147 / 117 | 147 / **8** |
| 5 | 09/10/2026 | `src/truth_validation.rs` | 118 / 93 | 118 / **32** |
| 5 | 09/10/2026 | `src/truth_schema.rs` | 69 / 44 | 69 / **8** |

**Tanda 1 (08/10/2026): 65 sin disparar menos (de 89 a 24; ninguna guarda
borrada, el número de comprobaciones no cambia).** Las 24 que quedan (5 + 7 + 12) están todas anotadas en el código
como inalcanzables, cada una con su porqué:

- *Perfil técnico (5).* La huella y el JSON ilegible de las dos
  `verify_*_profile_resource` comparan los bytes incrustados con una constante
  (cuatro); `required` sigue siempre a `exact_fields` sobre el mismo objeto
  (una). Las dos premisas las ejecutan pruebas, que no discriminan esas guardas.
- *FASTA suministrado (7).* El ordinal sin identificador del recorrido por
  flujo; `Reference::new` tras el bucle del lector que materializa; el «ID no
  declarado» de `finish_summary` (el lector lo rechaza antes con el mismo
  catálogo); la conversión de bases y `from_owned` de `finish_record`; y el
  «primer carácter no blanco» y la lectura de la definición de
  `parse_header_with_noodles`.
- *Mapa de reducción (12).* Cinco conversiones entre `usize` y `u64` que no
  fallan en 64 bits (cuatro `usize::try_from` y la de `segments.len()`); tres
  sumas que no desbordan (el recuento de segmentos, el orden de
  materialización y la coordenada derivada); la región ilegible del FAI; que
  los órdenes de selección formen `1..=n`, que ya implican su unicidad y su
  rango; la secuencia derivada ausente de `rebuilt`, que `rebuilt_order` ya
  garantizó; y el objetivo exacto, que implican las longitudes ya comprobadas.

Cómo se cubrieron: tres módulos de pruebas unitarias dentro de cada fichero
(las guardas son privadas o `pub(crate)` y el recurso del perfil está
incrustado), cada una **con el recurso o el documento real alterado en un solo
punto** y el motivo exacto. En el mapa de reducción el documento se escribe a
mano para un caso conocido y se compara `invariant` y `detail`; en el perfil
técnico, el JSON incrustado con un puntero; en el FASTA, ficheros mínimos para
el recorrido por flujo, el lector de resumen y el índice de aceptación, y un
resumen construido a mano para la construcción de una referencia sin bases.

**Verificación por mutación.** Sobre una copia fuera del repositorio se mutó
cada guarda cubierta (condición anulada con `if false && (..)`, motivo o código
cambiado, o constructora sustituida) y, además, cada mitad de las
disyunciones que una prueba tabular comparte con otras. Resultado: **88
ejecuciones y ninguna sobrevive en ese conjunto** (la ronda 2 amplió el
conjunto). Dos supervivientes iniciales eran guardas
que esta tanda no cubre (la longitud declarada de `finish_record` y la clase
derivada del juez del mapa): la primera la mata la prueba existente de
`tests/fasta_validation.rs`, que el filtro de la mutación no incluía; la
segunda no tenía prueba que exigiera su motivo y se añadió. Una mutación
(anular `selection_order == 0 || !insert`) dejó de insertar y mató también la
prueba del mapa válido; se repitió mitad por mitad con la inserción intacta.

**Ronda 2 tras la revisión adversarial (08/10/2026).** El revisor confirmó
las cifras y que las 24 inalcanzables lo son (fuzzers de 120 000 FASTA,
150 000 mapas y 100 000 perfiles: ninguna se disparó), pero encontró que «ninguna
mutación sobrevive» no valía fuera del conjunto probado. Se añadieron casos y
se mataron con mutación (19 ejecuciones más, una repetida): la contigüidad de
la proyección derivada (inicio y longitud, con un caso «desplazada y
compensada»), cinco guardas del juez del mapa que solo disparaba una prueba de
aceptación con `.is_err()` (orden de materialización, etiqueta y longitud
completo/residual, proyección, orden del catálogo derivado y campos ausentes o
desconocidos), los cuatro campos de un segmento de calidad, y los literales de
orientación y de ciclos máximos de PCR. Una mitad de condición,
la cobertura derivada de `validate_reduction_map_document`, es **equivalente**
(quitarla sobrevive porque lo impone la contigüidad) y queda anotada en el
código sin quitarla.

**Alcance de «0 pendientes».** Se refiere a las guardas que cuenta la puerta,
no a cada llamante de un ayudante genérico. Los ayudantes `required_text` y
`required_u64` del perfil técnico y `require_u64`, `require_exact_fields` y
`require_interval` del mapa de reducción cuentan como **una** guarda para
decenas de llamantes: la puerta las da por disparadas con un solo uso, y
borrar la comprobación de un llamante concreto no lo detecta salvo que su
literal tenga caso propio (se añadieron la orientación y los ciclos máximos de
PCR; el resto de literales no tienen caso individual).

**Lo que esta tanda no hace.** La puerta sigue midiendo el disparo, no la
discriminación; las pruebas de esta tanda sí están mutadas, y las anteriores a
ella no por ello lo están. Quedan los otros once jueces (la siguiente tanda son
los medianos).

**Tanda 2 (08/10/2026): 36 sin disparar menos (de 72 a 36; ninguna guarda
borrada, el número de comprobaciones no cambia).** Las 36 que quedan (26 + 10)
están todas anotadas en el código como inalcanzables, cada una con su porqué:

- *Juez del FASTQ final (36 → 26).* Cubiertas con prueba y motivo exacto: el
  resumen de cobertura (objetivo no finito o negativo, denominador cero, más `N`
  que bases y desbordamiento de la suma, cada mitad con su caso), el fichero
  más corto que el marcador de fin BGZF con su límite exacto (27 bytes se
  rechaza, 28 se acepta), y las seis negativas del replay cuando las entradas no
  son coherentes entre sí (plan de otra referencia, individuo de otra
  referencia, secuencias ausentes, individuo sin biblioteca, ninguna secuencia
  de la longitud de una lectura y media de Poisson no finita). Las 26 que quedan:
  dieciséis desbordamientos de contadores, sumas y conversiones que no fallan en
  64 bits (registros, ordinal, longitud de línea, cabecera, bases, `N`, tamaño
  acumulado, marcador de fin, longitudes de las secuencias y de las copias, y
  las sumas parciales y el inicio de la reproducción); las cuatro comprobaciones
  del lector de cabeceras del juez de la pareja (registro truncado, no UTF-8,
  sin salto de línea y sin `@`), que ya rechazó antes el juez de cada mate sobre
  el mismo fichero; «la cabecera no termina en salto de línea» y «la cabecera no
  tiene identidad»; la base ajena a `ACGTN` del replay, que la igualdad con la
  secuencia reproducida rechaza antes; la identidad del fichero que el juez de la
  pareja acaba de abrir; el candidato ausente de la reproducción, que implica la
  suma de longitudes; y la biblioteca ausente de `validate_declared_artifact`,
  que ya rechazó `ReplayedWgsPlan::new`.
- *Verificación del paquete F5 (36 → 10).* Cubiertas con prueba y motivo exacto
  (en el propio módulo, con un paquete de éxito y otro de fallo sintéticos
  alterados en un solo término, y una prueba de integración): cada término de la
  raíz, de la lista de enlaces, de la verdad, de los resúmenes del recuento, de las
  claves del recuento haplotípico, de la cota de bases elegibles, del efectivo, de
  las subpoblaciones, del fallo y de su diagnóstico; las rutas publicables; el
  fichero con nombre no UTF-8 y el directorio ilegible del inventario; el
  manifiesto inexistente; y que los límites sellados coincidan con los
  del escenario. Las 10 que quedan: `general_path` sin padre ni nombre (dos), la
  identidad del metadata del perfil, las cuatro auditorías de normalización
  (número y orden), el índice de bloque fuera de `u32`, el recuento desigual
  que el oráculo nunca devuelve, y tres fallos de E/S del recorrido del
  inventario (entrada del iterador, tipo del fichero y prefijo de la ruta).

**Pendientes: ninguno.** Ninguna guarda exigía cambiar código de producción.

**Defectos encontrados: ninguno.** Ninguna prueba reveló una guarda mal escrita.

**Verificación por mutación.** Sobre una copia fuera del repositorio se mutó
cada guarda cubierta, **cada mitad de sus disyunciones y cada motivo** (la
lección de la revisión de la tanda 1), con huella SHA antes, mutada y
restaurada: 117 ejecuciones en el paquete F5 y 113 en el FASTQ (una no compiló
y se rehízo). Además se mutó una muestra de guardas que ya estaban cubiertas por
pruebas que solo pedían error o que nadie ejercitaba como llamante. En la primera
ronda sobrevivieron 12 mutaciones reales del paquete F5 —la mitad del enlace por el
SHA-256, el tipo del fichero del inventario, las cinco comprobaciones de
identidad de la verdad y las tres de cada enlace por posición, el tamaño de la
verdad, un motivo y el motivo de la discrepancia entre diagnóstico y
manifiesto— y 5 del FASTQ —el alfabeto Base32 y la longitud exacta del UDI, la
unicidad solo en R1, el ordinal de la pareja repetida y el artefacto del mate
R1—: se añadieron casos y todas murieron. **Un ayudante que cuenta como una sola
guarda para varios llamantes** —`validate_final_artifact`, que sirve a R1 y a
R2— tenía el de R1 sin cubrir; se añadió su caso. Sobreviven cuatro, todas
**equivalentes o redundantes** y anotadas en el código: la mitad `!text.is_empty()`
de `publishable_relative` (un texto vacío ya lo rechaza el componente vacío); las
dos primeras mitades de `left.is_none() || right.is_none() || left != right` del
juez de la pareja (`left != right` ya es cierto con un solo `None`); `line.is_empty()`
de `read_line` (la última byte de una línea vacía ya no es un salto de línea); y
el rebobinado explícito del juez de un mate tras comprobar el marcador de fin,
redundante porque `reject_oversized_blocks` deja el fichero al principio (no es
una guarda de la puerta). La revisión adversarial encontró que la biblioteca
ausente del plan **no** era equivalente: con objetivo de cobertura, sin esa guarda
rechaza `derive_library_count_subseed` con otro error; se añadió su caso. Tampoco
se distingue `> 256` de `>= 256` en `read_final_line`: ninguna cabecera válida
llega a 256 bytes (máximo 155).

**Alcance de «0 pendientes» en esta tanda.** Igual que en la anterior se refiere
a las guardas que cuenta la puerta. En el paquete F5, `verify_link`,
`failure_mismatch` y `gate::Strict` cuentan como **una** guarda para varios
llamantes; se probó cada enlace del manifiesto, de la verdad y del fallo, cada
término del diagnóstico y las listas ya cubiertas por la prueba de la tanda
anterior. En el FASTQ, `validate_final_artifact` y `read_final_line` sirven a
varios llamantes; se probaron R1 y R2 del primero.

**Tanda 3 (08/10/2026): 98 sin disparar menos (de 123 a 25) en la matriz de
SNV, `src/population_snvs.rs`; ninguna guarda borrada, las 143 comprobaciones
siguen siendo 143.** Las 25 que quedan están todas anotadas en el código como
inalcanzables, cada una con su porqué junto a la guarda:

- *Cubiertas con prueba y motivo exacto (98).* Los analizadores de forma del TSV
  canónico y del suministrado (UTF-8, VCF, cabecera, número de columnas con su
  línea) y el de los metadatos (cada campo ausente, de tipo erróneo, de más o con
  literal distinto, y los tres ayudantes `required_object`, `required_string` y
  `required_literal` y el `closed_object` **por llamante**); la preparación de
  filas (contig, posición, alelos, población, REF distinto del de la referencia y
  posición fuera del contig, cada uno con su fila y valor); la normalización de
  frecuencias (tolerancia fuera de rango, frecuencia ilegible o fuera de [0, 1],
  exceso de suma con sus números exactos y el límite inclusivo); la autorización
  de la matriz vacía; el juez de la pareja reutilizada (referencia, huella, SHA
  del FASTA y de la matriz, orden canónico, poblaciones, referencia sin bases,
  metadatos o matriz que no se pueden leer); la restauración de las frecuencias
  desde la auditoría publicada; el juez de la auditoría (cada campo del grupo y
  del valor, identidad, orden, sumas, las tres decisiones y los pares descartados);
  y el lector del FASTA autoritativo. Las pruebas viven en el módulo
  `guard_tests` del propio fichero, que parte de una matriz y unos metadatos
  reales preparados con el código de producción y alterados **en un solo punto**,
  y en una prueba de integración de `tests/population_snvs_preparation.rs`.
- *Inalcanzables (25).* Trece comprobaciones de `validate_reusable_population_snvs_bytes`
  que repiten una forma de los metadatos que `parse_metadata_json` ya exigió
  sobre esos mismos bytes (objeto en la raíz, `reference`, `reference_fasta`,
  `matrix` y su `sha256`, `normalization_audit`, `populations` y cada
  `population_id`, `provenance`, y la lista, los objetos y los dos identificadores
  de los pares descartados); la lista de la restauración (la misma razón, y la
  comprobación gemela del juez público sí está cubierta); en la preparación, el
  recorrido de bases de `prepare_population_snvs_from_reference` (ya lo
  rechazó `PreparedReferenceIdentity::from_reference`), la población del segundo
  recorrido de filas (la rechazó el primero), la posición textual de la
  auditoría y la clave de identidad de cada valor, que proceden de filas ya
  analizadas, y el `find` de la cobertura de población, que usa el mismo
  predicado que el `any` anterior; tres derivaciones de `derive_matrix_variant_id`
  y la conversión `usize::try_from(position - 1)` del juez de la auditoría, que
  no fallan (componentes de 64 bytes como máximo, 64 bits); y dos de la
  normalización (el efectivo fuera de [0, 1] y la suma normalizada que se desvía
  de 1 más que la tolerancia), que acota el redondeo (un ulp) y cuyo intento
  ejecutado fue una búsqueda al azar de 1,19 millones de grupos sin ningún caso.
  La prueba `the_shape_guards_of_the_reuse_judge_are_shadowed_by_the_metadata_parser`
  ejecuta el intento sobre las de forma de los metadatos: cada metadato roto de esa forma lo
  rechaza antes el analizador con su propio texto (y por eso **no discrimina**
  las guardas inalcanzables que ilustra). Un metadato o una matriz publicados y
  manipulados **no** pueden llegar a ninguna de las 25: se comprobó una a una.

**Pendientes: ninguno.** Ninguna guarda exigía cambiar código de producción
(solo comentarios).

**Defectos encontrados: ninguno.** Ninguna prueba reveló una guarda mal
escrita.

**Verificación por mutación.** Sobre una copia fuera del repositorio, con huella
SHA antes, mutada y restaurada de cada una (el arnés se niega a trabajar bajo
`/home/ximo/devel` y comprueba que la restauración devuelve la huella original),
**249 ejecuciones**: cada guarda cubierta, cada mitad de sus disyunciones, cada
motivo (todos los `audit_error`, `missing` e `invalid` con su literal), cada
llamante de los ayudantes y de las clausuras `text` y `number`, y las
comparaciones por bits (`-0.0` frente a `0.0`) de las sumas, la tolerancia y, tras
la ronda 2, la frecuencia de referencia. 235 mueren, 1 no compiló y se
rehízo, 1 es el control (un cambio sin efecto, que **debe sobrevivir** y
sobrevive) y **12 sobreviven, todas demostradas**:

- *Inalcanzables (3 de los mutantes de motivo):* la lista de la restauración, la
  derivación del identificador de la auditoría y la conversión de la posición.
- *Mitades que no discriminan, anotadas en el código y sin quitar (9):*
  `is_finite()` de la tolerancia, de `reference_frequency` y `effective` en la
  restauración, de `original` en el juez de la auditoría y del lector `number`
  (en la tolerancia de `normalize_alt_frequencies` la mitad es
  redundante porque `NaN` y el infinito tampoco están en el intervalo; en los
  campos del JSON, `1e999` se lee como número pero `as_f64()` da `None`, no un
  infinito: lo demuestra `an_out_of_range_json_number_is_not_a_finite_float_for_the_audit_judge`);
  `!value.is_empty()` de `variant_id` (`VariantId::parse("")` ya falla con el mismo
  motivo); la mitad `effective_delta > tolerance` de la decisión `normalized`
  (inalcanzable por la cota de redondeo); y las dos mitades que se solapan del
  orden de grupos, `last >= key` y `!seen.insert(..)`: cada una por separado es
  redundante (la secuencia aceptada es estrictamente creciente, así que una clave
  repetida siempre la caza la otra), pero **debilitar las dos a la vez sí
  mata la prueba** (mutación `J-order-and-seen-both-weakened`).

**Lecciones aplicadas de las tandas 1 y 2 y lo que mostró la ronda 1.** Un primer
conjunto de pruebas dejó vivas 13 mutaciones reales: cuatro lecturas de texto de
la cabecera del grupo (`population_id`, `contig`, `position`, `ref`), cuya
caída a «vacío» daba el mismo motivo por otra guarda con un solo fallo (hizo falta
un caso con **dos fallos a la vez**, que distingue qué lectura falló primero); las
cuatro sumas y la tolerancia, que se comparaban por bits sin caso que distinguiera
el `-0.0`; y las guardas de duplicado, cobertura y exceso de suma que solo
disparaban pruebas existentes sin exigir el motivo. Todas se cubrieron y mueren.

**Ronda 2 tras la revisión adversarial (08/10/2026).** La revisión aprobó las 25
inalcanzables y los 12 supervivientes y halló cuatro cosas, ya resueltas:

- *`-0.0` en la rama `reference_remainder`.* Ninguna prueba distinguía `-0.0` de
  `0.0` en la frecuencia de referencia publicada (la afirmación anterior de que
  «todas las comparaciones por bits» estaban cubiertas era falsa para esta). Caso
  nuevo en `the_audit_judge_rejects_a_decision_that_the_values_do_not_support`
  (P2 chr2:2 con REF `-0.0`, por la vía directa y por la de bytes): la
  comparación numérica lo deja pasar y esa mutación **muere**.
- *Cambio de producción de una línea (decidido por el responsable).* En esa misma
  rama el juez comparaba `originals != effectives` numéricamente: unos metadatos
  manipulados con `original` `-0.0` y `effective` `0.0` pasaban toda la vía de
  reutilización. Ahora se compara por bits, como el resto de ramas; el caso
  (P2 chr1:2) está en la misma prueba y la versión numérica lo hace fallar. No
  cambia el número de comprobaciones (143).
- *Campos admitidos de cada población.* La prueba
  `a_population_with_every_admitted_field_is_accepted` acepta una población con
  los seis campos que admite la especificación; quitar cualquiera de la lista
  cerrada la hace fallar (se mutaron `data_version`, `filters`, `description` y
  `sample_size`).
- *Dos mitades `is_finite()` que no se habían mutado*: la de la frecuencia en
  `normalize_alt_frequencies` y la de `effective` en el juez de la auditoría.
  Sobreviven y son equivalentes (el intervalo cerrado ya excluye `NaN` y el
  infinito); la primera no estaba anotada y ahora lo está. Con ellas, **siete**
  mitades `is_finite()` equivalentes, y 14 supervivientes demostrados en total.

**Alcance de «0 pendientes» en esta tanda.** Se refiere a las guardas que cuenta
la puerta. `required_object`, `required_string`, `required_literal`,
`closed_object` y las clausuras `text` y `number` de la auditoría cuentan como
**una** guarda para varios llamantes: se probó cada llamante (tres, cinco, cuatro
y seis de los ayudantes de metadatos; siete campos de texto y seis numéricos de la
auditoría). No se vigila aquí la configuración del **generador** de matrices
(fuera de la puerta, ver más arriba).

**Tanda 4 (08/10/2026): 109 sin disparar menos (de 117 a 8) en el manifiesto de
evento, `src/event_manifest.rs`; ninguna guarda borrada, las 147 comprobaciones
siguen siendo 147.** Las 8 que quedan están todas anotadas junto a la guarda
como inalcanzables, con su porqué.

- *Cubiertas con prueba y motivo exacto (109).* `tests/event_manifest_guards.rs`
  parte de cuatro manifiestos de evento **publicados de verdad** por el programa
  (dos bibliotecas con contaminación biológica y salto de índice; con
  amplificación; con amplificación y microsatélites; con tartamudeo sin
  amplificar), altera cada uno **en un solo punto** y exige el motivo exacto, con
  control positivo, comprobación de que la alteración alteró y casos de «lado
  aceptado» (límites inclusivos, `-0.0` frente a `0.0` en las coberturas
  comparadas por bits). Tablas por contexto para los ayudantes compartidos
  (`required_object`, `required_array`, `validate_closed`: una clave de más, de
  menos y con otro nombre, **por llamante**; `parse_id`, `require_text` y
  `validate_hex`, por llamante; `locus_slip_counts` y
  `validate_demultiplexing_aggregates`, por lista), y los hechos del escritor
  (`validate_writer_facts`) armando un `EventManifest` incoherente.
- *Inalcanzables (8).* (1) la sección `artifacts` del literal `json!` que
  acaba de construir `event_manifest_document`; (2) «falta resumen FASTQ de una
  biblioteca completa» y (3) «un evento no completo no publica resúmenes FASTQ»,
  que otras comprobaciones anteriores ganan siempre; (4, 5) `cobertura propia` y
  `cobertura contaminante inválida` del cruce resumen-cobertura, que
  `validate_coverage` rechaza antes; (6) «descarte de demultiplexado con
  destino», que `validate_closed` rechaza antes; (7) «porcentaje de
  contaminación biológica ausente», que el control del porcentaje rechaza antes;
  (8) el `falta {campo}` de `required_value`, porque todo llamante lo precede un
  `validate_closed` sobre el mismo objeto. La prueba
  `unreachable_guards_are_shadowed_by_other_checks` ejecuta el intento y **no
  discrimina** esas guardas.

**Pendientes: ninguno.** Ninguna guarda exigía cambiar código de producción
(solo comentarios, más la comparación por bits de la ronda 2).

**Defectos encontrados: holguras menores, registradas en DT-046.** No son
guardas mal escritas en sentido estricto pero dejan pasar documentos que no
deberían; no se arreglan aquí.

**Verificación por mutación.** Sobre copias fuera del repositorio (el arnés se
niega a trabajar bajo `/home/ximo/devel`), **un mutante por ejecución**, con
huella SHA antes, mutada y restaurada de cada una: **540 ejecuciones**: cada
motivo (147), cada guarda entera y cada mitad de cada disyunción o conjunción
(215), cada llamante de `require_text`, `validate_closed`, `parse_id` y
`validate_hex` (70), comparaciones y límites a mano (91, incluidos `>=` frente a
`>`, rangos cerrados, aritmética que desborda y bits frente a valor), los diez
pares de las cinco igualdades de la contaminación y siete variantes que
conservan el efecto lateral de un `insert`. **493 mueren y 47 sobreviven
(tras añadir el caso del `-0.0` propio); todos demostrados**:

- *Inalcanzables (11 mutantes de las 8 guardas).*
- *Mitades equivalentes anotadas en el código (36).* Las de `is_finite()` (JSON no
  trae infinitos: cuatro lugares) y las de la tasa de quimera (`!is_finite`,
  `!(0..=1)`) y de los individuos repetidos del escritor, que la configuración ya
  excluye; `valid_manifest_id` en origen y destino (misma sintaxis que
  `LibraryId`); `>=` frente a `>` y `insert().is_some()` en los órdenes de
  individuos, coberturas, resúmenes, diagnósticos y salidas FASTQ (cada una la
  caza la otra; en las salidas, `ids.insert` caza antes la tupla repetida);
  `!paths.insert` de las salidas; las mitades de `artifact_mates` en el escritor
  y en el validador (**mutuamente** redundantes: quitar las dos a la vez muere);
  las cinco igualdades de conservación de la contaminación (**cada una se deduce
  de las otras cuatro**: ninguna sola muere, los diez pares sí); la suma de
  familias del histograma (el desborde de parejas llega a la vez y da el mismo
  motivo); las dos mitades del estado del evento que implica `validate_outputs`;
  la rama nula de `reference_tract_bases`; el índice del arreglo de errores; y
  `unwrap_or(0)` de los heredados.
- *Ronda 2 (09/10/2026), revisión adversarial.* Tres mutantes no equivalentes
  sobrevivían y se cubrieron: `len() != 64` frente a `< 64` en `validate_hex` y
  en `validate_subseed_hex` (casos de 65 caracteres: subsemilla del evento,
  huella de la referencia y las tres subsemillas de la contaminación) y
  `expects_lineage` sin `expects_pcr &&` (caso con `"stutter": {}` en
  `artifacts.pcr` de un evento sin amplificar, que exige «campos inválidos en
  artifacts.pcr»). Además, por decisión del responsable, **cambio de producción
  de una línea**: la fracción realizada de duplicados se compara ahora por bits
  (`to_bits()`), como el cruce resumen-cobertura, y un caso con `-0.0` falla con
  la comparación numérica. Los mecanismos de `is_finite()` quedan corregidos
  (`as_f64()` filtra lo no finito, también con `arbitrary_precision`). Las cifras
  de la puerta no cambian (147 / 8). Las siete ejecuciones de esta ronda: los
  cuatro mutantes nuevos mueren y los tres de `>=` frente a `>` del orden
  sobreviven (ya contados).

**Lecciones aplicadas.** Tres rechazos del escritor repetían el texto del
validador final (ajustes de calidad, mate ausente y cobertura de un evento no
completo): quitar la guarda del escritor daba el mismo motivo por otra vía y
**ninguna prueba lo notaba**; se añadieron casos con una incoherencia anterior
que el validador ve antes. Y un primer arnés por lotes daba **muertes falsas**
(un mutante hacía alcanzable la guarda inalcanzable de otro: el porcentaje
ausente): por eso las 540 ejecuciones son una a una.

**Alcance de «0 pendientes» en esta tanda.** Se refiere a las guardas que
cuenta la puerta. DT-043 (el juez de recursos publicados no ata la sección de
amplificación) y DT-036 (evento fallido sin estado en el manifiesto general) no
se tocan: los documentos de evento fallido que aquí se juzgan se derivan del
completo alterándolo, porque el programa no publica hoy ese documento por la vía
del ensamblado.

**Tanda 5 (09/10/2026): 61 sin disparar menos (de 93 a 32) en el juez semántico de
la verdad, `src/truth_validation.rs`, y 36 menos (de 44 a 8) en su contrato,
`src/truth_schema.rs`; ninguna guarda borrada, las comprobaciones siguen siendo
118 y 69.** Las 40 que quedan están todas anotadas junto a la guarda como
inalcanzables, con su porqué real.

| Módulo | Comprobaciones | Sin disparar antes | Sin disparar hoy |
|---|---|---|---|
| `src/truth_validation.rs` | 118 | 93 | **32** |
| `src/truth_schema.rs` | 69 | 44 | **8** |

- *Cubiertas con prueba y motivo exacto (97).* `tests/truth_validation_guards.rs`
  (84 pruebas) y `tests/truth_schema_guards.rs` (45) parten de una verdad **publicada
  de verdad** por el programa y la alteran en un solo punto. El escenario
  compartido (`tests/truth_guard_scenario/`) lleva dos secuencias (600 y 400
  bases), tres individuos de dos poblaciones, una matriz con un locus
  multialélico, cambios de base directos, inserciones y deleciones directas (una
  inserción en el punto uno, una deleción desde la base uno, una deleción que
  se come bases insertadas y de referencia, otra que solo se come insertadas, una
  inserción que sigue a otra), un tracto simple `(CA)50` y uno **compuesto**
  `(CAG)4(CA)9`. Se publica en versión 4 (sin repeticiones: la leen los dos
  lectores) y en versión 5 (con ellas: solo el lector completo, ver DT-015). Cada
  prueba exige el error exacto (variante, identidad, campo y valor); las del juez
  semántico, con **los dos lectores** cuando los dos pueden verlo (las tablas
  léxicas del contrato se probaron solo sobre el documento en memoria: el lector por
  flujo no aplica ese paso a los casos, ver DT-049), con control positivo y comprobación de
  que la alteración alteró. Las comprobaciones que el contrato adelanta
  (rechazos que el lector por bytes nunca ve porque `validate_structure` corre
  antes) se ejercitan con el documento **en memoria**, que es lo que recibe
  `validate_and_reconstruct_truth`; los documentos que armó a mano la prueba van
  con el fragmento vigente que resultaría de **ignorar la guarda**, para que sin
  ella no los pare un rechazo posterior con el mismo error. Dos pruebas en línea
  (`src/truth_validation.rs`) cubren las comprobaciones `pub(crate)`: una vista
  validada que se reutiliza contra otra referencia, y los jueces ante una
  referencia de otro identificador. Tablas por llamante para los ayudantes
  léxicos del contrato (`valid_id`: cada uno de sus llamantes, con el límite de 64
  y el de 65;
  `valid_sha`; `valid_base`; `valid_bases`; `valid_repeat_blocks`, los tres sitios
  donde se publica un bloque; `validate_modification_provenance`, ocho sitios).
- *Inalcanzables (40; 32 + 8).* Ningún `truth.jsonl` publicado y manipulado, ni
  documento en memoria, llega a ellas. Por clases:
  - *Conversiones de ancho* (`usize` y `u64` tienen el mismo ancho en x86_64):
    doce en el juez (réplica, materialización, locus de repetición, bloques,
    entero fuera de `usize`) y dos en el contrato. En un objetivo de 32 bits
    serían alcanzables.
  - *Desbordes de contador* que exigen 2^64 registros (copias y casos, lector
    completo y por flujo; 2 en el contrato) o un alelo con más de diez
    portadores (el paso léxico limita el índice de copia a 0..9).
  - *Sombreadas por otra comprobación anterior*: `BasesNotMaterialized` en
    `bases_of`, la réplica y la materialización (la huella de la referencia
    falla antes: `a_reference_without_bases_is_rejected_before_any_copy_is_rebuilt`);
    `DuplicateHeader` (la cabecera solo se acepta en la línea 1); dos
    modificaciones de segmento con el mismo evento (las rechaza
    `validate_structural_event_ids` antes); la longitud de un fragmento y la suma de longitudes dentro de
    `piece_boundary_index`, `split_replay_at` y el bucle de la deleción (las
    validó `replay_length` en la misma rama), el rechazo final de ambas funciones
    (el punto siempre cae en un fragmento) y el `else` de `let ... Source`;
    en el contrato, las bases vacías de un fragmento vigente o de un tracto (el
    paso léxico), el punto cero del evento fuente de una inserción y los
    extremos del evento fuente de una deleción (el alelo ya los exigió).
  La prueba `unreachable_guards_are_shadowed_by_other_checks` de cada fichero
  ejecuta el intento donde cabe y **no discrimina** esas guardas.
- *Pendientes: ninguno.* Ninguna guarda exigía cambiar código de producción
  (solo comentarios).

**Defectos encontrados.** **DT-047** (baja, **cerrada en la ronda 2**): el contrato
admitía un cambio de base de segmento cubierto por una deleción del mismo
segmento y la réplica lo rechazaba. El error estaba en el contrato
(`Segment::new` rechaza todo solape y `schemas/truth-v1.md` exige modificaciones
«ordenadas y no solapadas»): se quitó la excepción `deletion_covers_snp`; ahora
el contrato rechaza ese documento, en los dos órdenes y con su motivo exacto
(`ModificationsNotOrdered` y `OverlappingModifications`), y la guarda de la
réplica solo se alcanza con un documento en memoria. **DT-048** (media, nueva): el
juez no exige que el estado heredado recorra la referencia una vez y en orden
(rotaciones y duplicaciones se aceptan); se quitaron las dos pruebas que lo
consagraban. **DT-049** (media-baja, nueva): el lector por flujo no aplica a los
casos el paso léxico del contrato. Y se **reconfirma DT-015**: el lector por flujo
rechaza toda verdad de versión 5, de modo que ninguna prueba de esta tanda con
repeticiones lo usa.

**Verificación por mutación.** Sobre cinco copias fuera del repositorio (el arnés
se niega a trabajar bajo `/home/ximo/devel`, cada una con su directorio de
compilación y sus copias de seguridad fuera de git), **un mutante por ejecución**,
con la huella SHA antes, mutada y restaurada de cada una (todas coincidieron
entre «antes» y «restaurada», con `touch` tras restaurar). La pasada final sobre
el árbol de la ronda 2 fue de **746 ejecuciones** (la de la ronda 1 había sido de 752): cada motivo de los 187 rechazos,
cada guarda entera, cada mitad de cada disyunción o conjunción (también las
anidadas), cada operador de comparación y cada literal entero de las condiciones
(247 + 296 mutantes generados) y 186 escritos a mano para las **derivaciones**
que alimentan guardas (ancla de la inserción y su `.rev()`, `min`/`max` del
intervalo de una deleción, reparto de cambios al partir un fragmento, límites de
`<`/`<=` en la búsqueda del fragmento, filtros de copias por referencia,
individuo y secuencia, tolerancia del número de unidades, cada letra de los
alfabetos, `contained`, `overlaps`, `deletion_covers_snp`, `sort_key`). Resultado:
**595 mueren** (588 con las pruebas de integración, 7 con las en línea), **29 no
compilan** (se reescribieron conservando el efecto lateral, como en la ronda 1) y
**122 sobreviven**; tras añadir una prueba (tracto al final de su fragmento) mueren
dos más y quedan **120**, demostrados:
  - 54 son mutantes de los sitios inalcanzables de arriba.
  - 10 los mueren con la **batería entera** y no con estas pruebas (la rama de la
    matriz de `validate_individual_case_set`: población de las copias, alelos
    exactos de un caso SNP y de repetición, identificador del caso, filtros por
    secuencia de `expected_snp_alleles` y `expected_repeat_alleles`; las
    ejercita `tests/published_resource_validation.rs` y las pruebas de
    `execution`), y 1 con la prueba en línea que faltaba en la lista del arnés.
  - 55 son **equivalentes**, con su motivo anotado junto al código: solo
    `summaries.contains_key` y `seen_genotypes.insert` del orden del lector por
    flujo (`previous >= identidad` y `previous >= clave` **no** lo son: anular
    cada una sola muere por `two_copies_in_decreasing_order_are_out_of_order` y
    `two_genotypes_in_decreasing_order_are_out_of_order`); `population_id.is_empty()`
    y `variant_id.is_empty()` de la procedencia (los caza antes `valid_id`) y
    `historical_position == position` de la réplica inicial (lo caza antes
    `validate_structural_event_ids`); las tres igualdades de
    `segment_snp_matches_history` y la clase y el tipo del evento, que repiten
    el filtro de la réplica o la guarda previa con el **mismo error**;
    `*point == 0` y `start >= end` de la réplica (el siguiente rechazo da el
    mismo error); `<` frente a `<=` en `split_replay_at`, porque la frontera se
    ve desde los dos lados; `start >= end` de `validated_source`; el
    `motif_length == 0` de los bloques; `source_start >= 1`/`source_end >= 2`
    en el alelo de deleción; las bases vacías de `InvalidInsertion` y el
    `start >= end` de `InvalidDeletion` (el paso léxico y `contained` los
    rechazan antes); y las comparaciones de `overlaps`, `contained`, `sort_key` y
    la etiqueta de deleción de un punto que ninguna pareja ordenada distingue.
    Los dos que la revisión dio por equivalentes sin serlo
    (`other.start < self.start` y `other.start <= self.end` de `overlaps`)
    **mueren** con dos filas aceptadas nuevas (inserción estrictamente antes de
    una deleción, y deleción seguida de una inserción tras su final). La comparación `>` frente
    a `>=` de la tolerancia del número de unidades es equivalente **en la
    práctica** (no hay entrada de coma flotante en la que la diferencia sea
    exactamente la tolerancia).

**Lecciones aplicadas.** Una primera encuesta con un solo `chr1` dejó vivas las
mitades de `identity.sequence_id` de los tres filtros de copias: **se añadió una
segunda secuencia** al escenario. Otra, con tractos simples, dejó sin probar el
reparto de bloques: se añadió el tracto compuesto. Los mutantes de `.rev()`
del ancla, `min`/`max` del intervalo y el reparto de cambios al partir
exigieron **eventos directos nuevos** en el escenario (inserción tras otra
inserción, deleción que cruza una inserción, deleción solo de bases
insertadas, cambio de base justo en la frontera de una inserción). Cuatro
documentos armados a mano pasaban por la guarda y por el mutante con el mismo
error (el fragmento vigente que armé no era el que resultaba): se corrigieron
hasta que el mutante daba **otro** error. Y la comparación de ese día con los
supervivientes de la batería entera evitó dar por sin cubrir 10 guardas que
cubren pruebas de otros módulos.

Tras la pasada final se amplió una prueba (`an_inherited_base_change_covered_by_a_segment_deletion_does_not_replay`,
con el orden inverso): solo añade una aserción, así que no puede restar
discriminación.


**Tanda 6, primera parte (10/10/2026): 120 sin disparar menos (de 618 a 498) en
`src/results.rs`, con microsatélites, tartamudeo y PCR; ninguna guarda borrada,
las comprobaciones siguen en 779.** `results.rs` es el juez de los recursos
publicados, de la configuración efectiva y del manifiesto general; es el
módulo con más guardas sin disparar y se cubre por temas en varias tandas. Esta
toca diez funciones y sus ayudantes privados exclusivos: 130 guardas sin
disparar al empezar (113 en las diez funciones y 17 en los ayudantes); quedan
**12**, todas anotadas junto a la guarda como inalcanzables, con su porqué.

| Función | Sin disparar antes | Hoy |
|---|---|---|
| `validate_repeat_loci_section` (sección de loci del manifiesto) | 17 | 4 |
| `validate_published_slippage_fans` | 7 | 2 |
| `published_ancestral_tract` | 13 | 2 |
| `compare_published_fan_against_oracle` | 6 | 0 |
| `validate_repeat_slippage` (bloque del generador) | 7 | 0 |
| `validate_repeat_allele` | 21 | 1 |
| `validate_effective_stutter` | 8 | 0 |
| `validate_effective_tract_length_dependence` | 11 | 0 |
| `validate_effective_clonal_stutter` | 12 | 0 |
| `validate_effective_pcr` | 11 | 0 |
| ayudantes: `recompute_slippage_fan`, `validate_repeat_tract_output`, `validate_repeat_admission`, `repeat_motif_list`, `repeat_u64_list` | 17 | 3 |
| **Total del alcance** | **130** | **12** |

(La tabla cuenta guardas con la numeración de la medición del 10/10/2026; las 618
y 498 son las de la puerta sobre todo el módulo.)

- *Cubiertas con prueba y motivo exacto (118).* `tests/results_repeat_guards.rs`
  (24 pruebas) parte de paquetes **publicados de verdad** por el programa y los
  altera en un solo punto: `mixed` (un locus declarado con un alelo interrumpido
  y otro puro, y uno generado por deslizamiento), `edges` (cuatro loci generados
  cuyos linajes llegan a los bordes: tirada agotada, tope, expansión que ya no
  cabe y un tracto interrumpido), `inserting` (inserción de un tracto en la
  referencia y un locus generado sobre la final), la ploidía 4, y el escenario
  amplificado con microsatélites (tartamudeo de base y clonal, con y sin
  multiplicador de los homopolímeros) y sin ellos. El lector del manifiesto
  juzga la sección de loci (`validate_general_manifest_bytes`); el juez completo
  (`validate_published_resources`), el abanico; y la estructura de las
  secciones efectivas, `validate_effective_config_against_manifest_bytes` con la
  efectiva alterada y el manifiesto intacto (corre antes de todo lo que ata la
  efectiva al TOML). Cada tabla lleva control positivo, comprueba que **la
  alteración alteró** y fija el lado aceptado de los límites (inicio uno, final
  igual al inicio, extremos cerrados de pureza y sesgo, cuatro clases de
  eficiencia, tres nudos, mínimo de ciclos de 32 bits justo) y de las listas
  cerradas; las frecuencias y las fracciones se comparan al bit con el decimal
  siguiente, y un cero negativo se rechaza donde la frecuencia realizada es
  cero. Ayudantes compartidos con casos por llamante (`validate_repeat_admission`,
  dos; `repeat_motif_list` y `repeat_u64_list`, locus y alelo).
- *Inalcanzables (12).* Ningún paquete publicado y manipulado llega a ellas:
  - *Sombreadas por una comprobación anterior (5):* los `ok_or_else` de los
    mínimos de omisión, de la contabilidad de tracto, de los mínimos del locus y
    de la frecuencia realizada del alelo (el cierre exige las claves exactas, de
    modo que ninguna falta), y el origen del alelo (lo exige antes
    `validate_repeat_allele`). Intento ejecutado:
    `unreachable_guards_are_shadowed_by_other_checks`.
  - *Sombreada por la lectura de la referencia (1):* el identificador de la
    referencia (lo exigen `read_physical_reference` y `judge_inserted_reference`).
    Intento ejecutado: `unreachable_fan_guards_are_shadowed_by_other_checks`.
    **Las bases no materializadas sí se alcanzan** (corrección de la revisión
    adversarial): con un paquete rehecho entero sin poblaciones ni matriz (copia
    del TOML, su huella y el identificador de ejecución, efectiva, manifiesto, dos
    roles y dos artefactos) el juez completo rechaza con ese motivo exacto; lo
    cubre `a_package_forged_without_populations_reaches_the_unmaterialized_reference`
    y su mutante muere.
  - *Conversiones de ancho (2):* `u64` a `usize`, iguales en x86_64.
  - *Imposibles en la práctica (4):* el recuento de Poisson fallido (la media está
    acotada a [0, 1000] y el error exige 2^64 iteraciones o un sorteo a menos de
    2^-52 del extremo, sobre una semilla atada al TOML), el desborde de la suma de
    pasos y el paso de tracto fallido (el filtro de viabilidad exige una tirada de
    al menos dos unidades); y, además, la derivación de la subsemilla, que solo falla con
    un identificador de locus de más de 4 GiB (**no es imposible en lógica**, no
    cabe en una prueba).
- *Mitades equivalentes demostradas.* Los dos atajos de sesgo 1 y sesgo 0 de
  `recompute_slippage_fan` (`f4_bernoulli` ya devuelve sin consumir azar con
  probabilidad 1 o 0: quitar cada uno, o los dos, da el mismo recorrido; anotados
  junto al código y conservados). Las seis mitades `!x.is_finite()` de los
  decimales del generador, de los mínimos, del alelo y del recuento de copias por
  bloque: con `serde_json` y `arbitrary_precision`, `as_f64` devuelve `None` para
  un número que no cabe en un `f64` finito, así que un decimal no finito no llega
  a la comprobación; lo que rechaza `1e999` es la lectura, y tiene caso en cada
  sitio. Además: la rama nula de `reference_tract_bases` (`null` da `None` con
  `as_u64`), `>` frente a `>=` de la tolerancia de las copias por bloque (no hay
  decimal de entrada con diferencia exactamente igual a la tolerancia) y
  `stopped = true` del recorrido (el filtro de viabilidad es idempotente sobre un
  estado que no cambia).
- *Pendientes: ninguno.* Ninguna guarda exigía cambiar código de producción (solo
  comentarios).

**Defecto encontrado: DT-051** (media; sin arreglar en esta tanda, cerrada el
10/10/2026 en un cambio propio). Al preguntarse qué cambiaría
alguien para que el juez completo aceptara un manifiesto falso, se ejecutó y el
juez **acepta** alterar la sección de loci del manifiesto (secuencias, frecuencias,
pureza, unidades, intervalos, topes, mínimos, procedencia de la referencia,
parámetros del generador) y, en la efectiva, los loci declarados y los campos de los
generados que su juez no lee. Ver DT-051.

**Verificación por mutación.** Sobre cuatro copias fuera del repositorio (el arnés
se niega a trabajar bajo `/home/ximo/devel`, cada una con su directorio de
compilación y sus copias de seguridad fuera de git), **un mutante por ejecución**,
con la huella SHA antes (`1ffd7ba82c…`), mutada y restaurada de cada una (todas
coincidieron entre «antes» y «restaurada», con `touch` tras restaurar y comprobación
final contra el árbol real). La pasada final fue de **623 ejecuciones** sobre el árbol
final: el motivo de cada rechazo, cada guarda entera, cada mitad de cada
disyunción o conjunción (también las anidadas), cada operador de comparación y cada
literal entero de las condiciones, los cierres `filter`/`is_some_and`/`is_none_or`
y 184 mutantes a mano de derivaciones y límites (catálogo derivado de la presencia
del tartamudeo clonal, iteraciones de la inversión, longitud de la tabla de factores,
suelo de la tolerancia, aritmética comprobada frente a la que da la vuelta,
viabilidad del paso, sentido del sesgo, orden canónico de los alelos, ploidía de los
linajes, identificador de la referencia final). Resultado: **588 mueren**, **8 no
compilan** (se reescribieron a mano conservando el efecto y mueren) y **27
sobreviven**, todos demostrados: 12 son los rechazos de los sitios inalcanzables, 6
las mitades `is_finite`, 5 mutan texto de un comentario que quedó dentro de la
condición, y 4 equivalentes (`>=` de la tolerancia, la rama nula, `stopped = true` y
una mutación mía que era la identidad). La primera pasada (606 mutantes) dejó 74
vivos; la mayoría eran guardas que solo cubrían pruebas de otros ficheros, o a las
que les faltaba un lado del límite, y se añadieron casos (linajes distintos de
las copias, tope de media por los dos lados, inicio uno, lista de loci vacía,
procedencia de la referencia, listas de bloque más largas y más cortas que el
recuento, bases de bloque por encima y por debajo, tolerancia con suelo uno,
marcas de tosquedad que no son booleanos, fracción y `m` negativas con fracción
cero) y los tres escenarios que alcanzan los bordes del recorrido, la inserción y la
ploidía 4.

**Ronda 2 (revisión adversarial).** Se corrigió la inalcanzable de las bases no
materializadas (ver arriba) y se anotaron los dos atajos de sesgo como
equivalentes. Cinco mutantes sobre el árbol corregido, huella `e3e9182a69…` antes y
restaurada: el motivo de las bases no materializadas y su guarda neutralizada
**mueren** con la prueba nueva; los tres de los atajos de sesgo (los dos, solo el
de 1, solo el de 0) **sobreviven**, equivalentes. La revisión confirmó además que
los cinco supervivientes «de comentario» eran ruido del generador (repetían las
mitades reales de la condición, que mueren).

**Tanda 7, segunda parte de `results.rs` (11/10/2026): 104 sin disparar menos (de
497 a 393) y 97 de ellas en la configuración efectiva; ninguna guarda borrada, las
comprobaciones siguen en 787.** Las otras 7 son guardas de ayudantes compartidos y
del lector del manifiesto que alcanzan los intentos de llegar a las inalcanzables.
`results.rs` es el juez de los recursos publicados,
de la configuración efectiva y del manifiesto general. Esta tanda toca diez funciones
de la configuración efectiva y sus ayudantes privados exclusivos: 161 guardas sin
disparar al empezar (123 en código vivo y 38 en las dos funciones de la versión 1, que
eran **código muerto**, ver abajo); quedan **26 de las 123 vivas**, todas anotadas
junto a la guarda como inalcanzables con su porqué real. Las de la versión 1 desaparecen
con el código (borrado el 11/10/2026).

| Función | Sin disparar antes | Hoy |
|---|---|---|
| `validate_effective_config_v3_against_manifest_bytes` (la efectiva frente al manifiesto) | 15 | 7 |
| `validate_effective_technical_event_link` (enlace del evento técnico) | 7 | 3 |
| `validate_effective_targets_against_configuration` (la efectiva atada al TOML copiado) | 21 | 11 |
| `validate_effective_initial_populations` | 10 | 0 |
| `validate_effective_direct_events` | 9 | 0 |
| `validate_effective_technical_event` (sin amplificación) | 25 | 1 |
| `validate_effective_reference` | 19 | 2 |
| `validate_effective_derived_sequences` | 9 | 2 |
| ayudantes: `validate_effective_reduction`, `validate_reversible_text` | 8 | 0 |
| **Total vivo** | **123** | **26** |
| `validate_effective_config_v1_against_manifest_bytes`, `validate_effective_targets_against_configuration_v1` (**código muerto, borrado**) | 38 | — |

- *Las funciones de la versión 1 eran código muerto y se borraron (11/10/2026, decisión
  del responsable).* Llevaban `#[allow(dead_code)]` y **ninguna ruta de producción ni de
  prueba las llamaba**: `validate_effective_config_v1_against_manifest_bytes` (16 guardas),
  `validate_effective_targets_against_configuration_v1` (22; solo la llamaba
  `validate_published_resources_v1`), `validate_published_resources_v1` (21),
  `parse_general_manifest_v1_object` (7) y `retag_v1_pair_error`: 66 comprobaciones, todas
  sin disparar, y 685 líneas. Los esquemas v1 de `schemas/` son ficheros congelados y no se
  tocaron. DTC-036 se actualizó a las versiones que el programa lee de verdad.
- *Cubiertas con prueba y motivo exacto (97).* `tests/truth_validation_guards.rs`
  (84 pruebas) y `tests/truth_schema_guards.rs` (45) parten de una verdad **publicada
  de verdad** por el programa y la alteran en un solo punto. El escenario
  compartido (`tests/truth_guard_scenario/`) lleva dos secuencias (600 y 400
  bases), tres individuos de dos poblaciones, una matriz con un locus
  multialélico, cambios de base directos, inserciones y deleciones directas (una
  inserción en el punto uno, una deleción desde la base uno, una deleción que
  se come bases insertadas y de referencia, otra que solo se come insertadas, una
  inserción que sigue a otra), un tracto simple `(CA)50` y uno **compuesto**
  `(CAG)4(CA)9`. Se publica en versión 4 (sin repeticiones: la leen los dos
  lectores) y en versión 5 (con ellas: solo el lector completo, ver DT-015). Cada
  prueba exige el error exacto (variante, identidad, campo y valor); las del juez
  semántico, con **los dos lectores** cuando los dos pueden verlo (las tablas
  léxicas del contrato se probaron solo sobre el documento en memoria: el lector por
  flujo no aplica ese paso a los casos, ver DT-049), con control positivo y comprobación de
  que la alteración alteró. Las comprobaciones que el contrato adelanta
  (rechazos que el lector por bytes nunca ve porque `validate_structure` corre
  antes) se ejercitan con el documento **en memoria**, que es lo que recibe
  `validate_and_reconstruct_truth`; los documentos que armó a mano la prueba van
  con el fragmento vigente que resultaría de **ignorar la guarda**, para que sin
  ella no los pare un rechazo posterior con el mismo error. Dos pruebas en línea
  (`src/truth_validation.rs`) cubren las comprobaciones `pub(crate)`: una vista
  validada que se reutiliza contra otra referencia, y los jueces ante una
  referencia de otro identificador. Tablas por llamante para los ayudantes
  léxicos del contrato (`valid_id`: cada uno de sus llamantes, con el límite de 64
  y el de 65;
  `valid_sha`; `valid_base`; `valid_bases`; `valid_repeat_blocks`, los tres sitios
  donde se publica un bloque; `validate_modification_provenance`, ocho sitios).
- *Inalcanzables (40; 32 + 8).* Ningún `truth.jsonl` publicado y manipulado, ni
  documento en memoria, llega a ellas. Por clases:
  - *Conversiones de ancho* (`usize` y `u64` tienen el mismo ancho en x86_64):
    doce en el juez (réplica, materialización, locus de repetición, bloques,
    entero fuera de `usize`) y dos en el contrato. En un objetivo de 32 bits
    serían alcanzables.
  - *Desbordes de contador* que exigen 2^64 registros (copias y casos, lector
    completo y por flujo; 2 en el contrato) o un alelo con más de diez
    portadores (el paso léxico limita el índice de copia a 0..9).
  - *Sombreadas por otra comprobación anterior*: `BasesNotMaterialized` en
    `bases_of`, la réplica y la materialización (la huella de la referencia
    falla antes: `a_reference_without_bases_is_rejected_before_any_copy_is_rebuilt`);
    `DuplicateHeader` (la cabecera solo se acepta en la línea 1); dos
    modificaciones de segmento con el mismo evento (las rechaza
    `validate_structural_event_ids` antes); la longitud de un fragmento y la suma de longitudes dentro de
    `piece_boundary_index`, `split_replay_at` y el bucle de la deleción (las
    validó `replay_length` en la misma rama), el rechazo final de ambas funciones
    (el punto siempre cae en un fragmento) y el `else` de `let ... Source`;
    en el contrato, las bases vacías de un fragmento vigente o de un tracto (el
    paso léxico), el punto cero del evento fuente de una inserción y los
    extremos del evento fuente de una deleción (el alelo ya los exigió).
  La prueba `unreachable_guards_are_shadowed_by_other_checks` de cada fichero
  ejecuta el intento donde cabe y **no discrimina** esas guardas.
- *Pendientes: ninguno.* Ninguna guarda exigía cambiar código de producción
  (solo comentarios).

**Defectos encontrados.** **DT-047** (baja, **cerrada en la ronda 2**): el contrato
admitía un cambio de base de segmento cubierto por una deleción del mismo
segmento y la réplica lo rechazaba. El error estaba en el contrato
(`Segment::new` rechaza todo solape y `schemas/truth-v1.md` exige modificaciones
«ordenadas y no solapadas»): se quitó la excepción `deletion_covers_snp`; ahora
el contrato rechaza ese documento, en los dos órdenes y con su motivo exacto
(`ModificationsNotOrdered` y `OverlappingModifications`), y la guarda de la
réplica solo se alcanza con un documento en memoria. **DT-048** (media, nueva): el
juez no exige que el estado heredado recorra la referencia una vez y en orden
(rotaciones y duplicaciones se aceptan); se quitaron las dos pruebas que lo
consagraban. **DT-049** (media-baja, nueva): el lector por flujo no aplica a los
casos el paso léxico del contrato. Y se **reconfirma DT-015**: el lector por flujo
rechaza toda verdad de versión 5, de modo que ninguna prueba de esta tanda con
repeticiones lo usa.

**Verificación por mutación.** Sobre cinco copias fuera del repositorio (el arnés
se niega a trabajar bajo `/home/ximo/devel`, cada una con su directorio de
compilación y sus copias de seguridad fuera de git), **un mutante por ejecución**,
con la huella SHA antes, mutada y restaurada de cada una (todas coincidieron
entre «antes» y «restaurada», con `touch` tras restaurar). La pasada final sobre
el árbol de la ronda 2 fue de **746 ejecuciones** (la de la ronda 1 había sido de 752): cada motivo de los 187 rechazos,
cada guarda entera, cada mitad de cada disyunción o conjunción (también las
anidadas), cada operador de comparación y cada literal entero de las condiciones
(247 + 296 mutantes generados) y 186 escritos a mano para las **derivaciones**
que alimentan guardas (ancla de la inserción y su `.rev()`, `min`/`max` del
intervalo de una deleción, reparto de cambios al partir un fragmento, límites de
`<`/`<=` en la búsqueda del fragmento, filtros de copias por referencia,
individuo y secuencia, tolerancia del número de unidades, cada letra de los
alfabetos, `contained`, `overlaps`, `deletion_covers_snp`, `sort_key`). Resultado:
**595 mueren** (588 con las pruebas de integración, 7 con las en línea), **29 no
compilan** (se reescribieron conservando el efecto lateral, como en la ronda 1) y
**122 sobreviven**; tras añadir una prueba (tracto al final de su fragmento) mueren
dos más y quedan **120**, demostrados:
  - 54 son mutantes de los sitios inalcanzables de arriba.
  - 10 los mueren con la **batería entera** y no con estas pruebas (la rama de la
    matriz de `validate_individual_case_set`: población de las copias, alelos
    exactos de un caso SNP y de repetición, identificador del caso, filtros por
    secuencia de `expected_snp_alleles` y `expected_repeat_alleles`; las
    ejercita `tests/published_resource_validation.rs` y las pruebas de
    `execution`), y 1 con la prueba en línea que faltaba en la lista del arnés.
  - 55 son **equivalentes**, con su motivo anotado junto al código: solo
    `summaries.contains_key` y `seen_genotypes.insert` del orden del lector por
    flujo (`previous >= identidad` y `previous >= clave` **no** lo son: anular
    cada una sola muere por `two_copies_in_decreasing_order_are_out_of_order` y
    `two_genotypes_in_decreasing_order_are_out_of_order`); `population_id.is_empty()`
    y `variant_id.is_empty()` de la procedencia (los caza antes `valid_id`) y
    `historical_position == position` de la réplica inicial (lo caza antes
    `validate_structural_event_ids`); las tres igualdades de
    `segment_snp_matches_history` y la clase y el tipo del evento, que repiten
    el filtro de la réplica o la guarda previa con el **mismo error**;
    `*point == 0` y `start >= end` de la réplica (el siguiente rechazo da el
    mismo error); `<` frente a `<=` en `split_replay_at`, porque la frontera se
    ve desde los dos lados; `start >= end` de `validated_source`; el
    `motif_length == 0` de los bloques; `source_start >= 1`/`source_end >= 2`
    en el alelo de deleción; las bases vacías de `InvalidInsertion` y el
    `start >= end` de `InvalidDeletion` (el paso léxico y `contained` los
    rechazan antes); y las comparaciones de `overlaps`, `contained`, `sort_key` y
    la etiqueta de deleción de un punto que ninguna pareja ordenada distingue.
    Los dos que la revisión dio por equivalentes sin serlo
    (`other.start < self.start` y `other.start <= self.end` de `overlaps`)
    **mueren** con dos filas aceptadas nuevas (inserción estrictamente antes de
    una deleción, y deleción seguida de una inserción tras su final). La comparación `>` frente
    a `>=` de la tolerancia del número de unidades es equivalente **en la
    práctica** (no hay entrada de coma flotante en la que la diferencia sea
    exactamente la tolerancia).

**Lecciones aplicadas.** Una primera encuesta con un solo `chr1` dejó vivas las
mitades de `identity.sequence_id` de los tres filtros de copias: **se añadió una
segunda secuencia** al escenario. Otra, con tractos simples, dejó sin probar el
reparto de bloques: se añadió el tracto compuesto. Los mutantes de `.rev()`
del ancla, `min`/`max` del intervalo y el reparto de cambios al partir
exigieron **eventos directos nuevos** en el escenario (inserción tras otra
inserción, deleción que cruza una inserción, deleción solo de bases
insertadas, cambio de base justo en la frontera de una inserción). Cuatro
documentos armados a mano pasaban por la guarda y por el mutante con el mismo
error (el fragmento vigente que armé no era el que resultaba): se corrigieron
hasta que el mutante daba **otro** error. Y la comparación de ese día con los
supervivientes de la batería entera evitó dar por sin cubrir 10 guardas que
cubren pruebas de otros módulos.

Tras la pasada final se amplió una prueba (`an_inherited_base_change_covered_by_a_segment_deletion_does_not_replay`,
con el orden inverso): solo añade una aserción, así que no puede restar
discriminación.


**Tanda 6, primera parte (10/10/2026): 120 sin disparar menos (de 618 a 498) en
`src/results.rs`, con microsatélites, tartamudeo y PCR; ninguna guarda borrada,
las comprobaciones siguen en 779.** `results.rs` es el juez de los recursos
publicados, de la configuración efectiva y del manifiesto general; es el
módulo con más guardas sin disparar y se cubre por temas en varias tandas. Esta
toca diez funciones y sus ayudantes privados exclusivos: 130 guardas sin
disparar al empezar (113 en las diez funciones y 17 en los ayudantes); quedan
**12**, todas anotadas junto a la guarda como inalcanzables, con su porqué.

| Función | Sin disparar antes | Hoy |
|---|---|---|
| `validate_repeat_loci_section` (sección de loci del manifiesto) | 17 | 4 |
| `validate_published_slippage_fans` | 7 | 2 |
| `published_ancestral_tract` | 13 | 2 |
| `compare_published_fan_against_oracle` | 6 | 0 |
| `validate_repeat_slippage` (bloque del generador) | 7 | 0 |
| `validate_repeat_allele` | 21 | 1 |
| `validate_effective_stutter` | 8 | 0 |
| `validate_effective_tract_length_dependence` | 11 | 0 |
| `validate_effective_clonal_stutter` | 12 | 0 |
| `validate_effective_pcr` | 11 | 0 |
| ayudantes: `recompute_slippage_fan`, `validate_repeat_tract_output`, `validate_repeat_admission`, `repeat_motif_list`, `repeat_u64_list` | 17 | 3 |
| **Total del alcance** | **130** | **12** |

(La tabla cuenta guardas con la numeración de la medición del 10/10/2026; las 618
y 498 son las de la puerta sobre todo el módulo.)

- *Cubiertas con prueba y motivo exacto (118).* `tests/results_repeat_guards.rs`
  (24 pruebas) parte de paquetes **publicados de verdad** por el programa y los
  altera en un solo punto: `mixed` (un locus declarado con un alelo interrumpido
  y otro puro, y uno generado por deslizamiento), `edges` (cuatro loci generados
  cuyos linajes llegan a los bordes: tirada agotada, tope, expansión que ya no
  cabe y un tracto interrumpido), `inserting` (inserción de un tracto en la
  referencia y un locus generado sobre la final), la ploidía 4, y el escenario
  amplificado con microsatélites (tartamudeo de base y clonal, con y sin
  multiplicador de los homopolímeros) y sin ellos. El lector del manifiesto
  juzga la sección de loci (`validate_general_manifest_bytes`); el juez completo
  (`validate_published_resources`), el abanico; y la estructura de las
  secciones efectivas, `validate_effective_config_against_manifest_bytes` con la
  efectiva alterada y el manifiesto intacto (corre antes de todo lo que ata la
  efectiva al TOML). Cada tabla lleva control positivo, comprueba que **la
  alteración alteró** y fija el lado aceptado de los límites (inicio uno, final
  igual al inicio, extremos cerrados de pureza y sesgo, cuatro clases de
  eficiencia, tres nudos, mínimo de ciclos de 32 bits justo) y de las listas
  cerradas; las frecuencias y las fracciones se comparan al bit con el decimal
  siguiente, y un cero negativo se rechaza donde la frecuencia realizada es
  cero. Ayudantes compartidos con casos por llamante (`validate_repeat_admission`,
  dos; `repeat_motif_list` y `repeat_u64_list`, locus y alelo).
- *Inalcanzables (12).* Ningún paquete publicado y manipulado llega a ellas:
  - *Sombreadas por una comprobación anterior (5):* los `ok_or_else` de los
    mínimos de omisión, de la contabilidad de tracto, de los mínimos del locus y
    de la frecuencia realizada del alelo (el cierre exige las claves exactas, de
    modo que ninguna falta), y el origen del alelo (lo exige antes
    `validate_repeat_allele`). Intento ejecutado:
    `unreachable_guards_are_shadowed_by_other_checks`.
  - *Sombreada por la lectura de la referencia (1):* el identificador de la
    referencia (lo exigen `read_physical_reference` y `judge_inserted_reference`).
    Intento ejecutado: `unreachable_fan_guards_are_shadowed_by_other_checks`.
    **Las bases no materializadas sí se alcanzan** (corrección de la revisión
    adversarial): con un paquete rehecho entero sin poblaciones ni matriz (copia
    del TOML, su huella y el identificador de ejecución, efectiva, manifiesto, dos
    roles y dos artefactos) el juez completo rechaza con ese motivo exacto; lo
    cubre `a_package_forged_without_populations_reaches_the_unmaterialized_reference`
    y su mutante muere.
  - *Conversiones de ancho (2):* `u64` a `usize`, iguales en x86_64.
  - *Imposibles en la práctica (4):* el recuento de Poisson fallido (la media está
    acotada a [0, 1000] y el error exige 2^64 iteraciones o un sorteo a menos de
    2^-52 del extremo, sobre una semilla atada al TOML), el desborde de la suma de
    pasos y el paso de tracto fallido (el filtro de viabilidad exige una tirada de
    al menos dos unidades); y, además, la derivación de la subsemilla, que solo falla con
    un identificador de locus de más de 4 GiB (**no es imposible en lógica**, no
    cabe en una prueba).
- *Mitades equivalentes demostradas.* Los dos atajos de sesgo 1 y sesgo 0 de
  `recompute_slippage_fan` (`f4_bernoulli` ya devuelve sin consumir azar con
  probabilidad 1 o 0: quitar cada uno, o los dos, da el mismo recorrido; anotados
  junto al código y conservados). Las seis mitades `!x.is_finite()` de los
  decimales del generador, de los mínimos, del alelo y del recuento de copias por
  bloque: con `serde_json` y `arbitrary_precision`, `as_f64` devuelve `None` para
  un número que no cabe en un `f64` finito, así que un decimal no finito no llega
  a la comprobación; lo que rechaza `1e999` es la lectura, y tiene caso en cada
  sitio. Además: la rama nula de `reference_tract_bases` (`null` da `None` con
  `as_u64`), `>` frente a `>=` de la tolerancia de las copias por bloque (no hay
  decimal de entrada con diferencia exactamente igual a la tolerancia) y
  `stopped = true` del recorrido (el filtro de viabilidad es idempotente sobre un
  estado que no cambia).
- *Pendientes: ninguno.* Ninguna guarda exigía cambiar código de producción (solo
  comentarios).

**Defecto encontrado: DT-051** (media; sin arreglar en esta tanda, cerrada el
10/10/2026 en un cambio propio). Al preguntarse qué cambiaría
alguien para que el juez completo aceptara un manifiesto falso, se ejecutó y el
juez **acepta** alterar la sección de loci del manifiesto (secuencias, frecuencias,
pureza, unidades, intervalos, topes, mínimos, procedencia de la referencia,
parámetros del generador) y, en la efectiva, los loci declarados y los campos de los
generados que su juez no lee. Ver DT-051.

**Verificación por mutación.** Sobre cuatro copias fuera del repositorio (el arnés
se niega a trabajar bajo `/home/ximo/devel`, cada una con su directorio de
compilación y sus copias de seguridad fuera de git), **un mutante por ejecución**,
con la huella SHA antes (`1ffd7ba82c…`), mutada y restaurada de cada una (todas
coincidieron entre «antes» y «restaurada», con `touch` tras restaurar y comprobación
final contra el árbol real). La pasada final fue de **623 ejecuciones** sobre el árbol
final: el motivo de cada rechazo, cada guarda entera, cada mitad de cada
disyunción o conjunción (también las anidadas), cada operador de comparación y cada
literal entero de las condiciones, los cierres `filter`/`is_some_and`/`is_none_or`
y 184 mutantes a mano de derivaciones y límites (catálogo derivado de la presencia
del tartamudeo clonal, iteraciones de la inversión, longitud de la tabla de factores,
suelo de la tolerancia, aritmética comprobada frente a la que da la vuelta,
viabilidad del paso, sentido del sesgo, orden canónico de los alelos, ploidía de los
linajes, identificador de la referencia final). Resultado: **588 mueren**, **8 no
compilan** (se reescribieron a mano conservando el efecto y mueren) y **27
sobreviven**, todos demostrados: 12 son los rechazos de los sitios inalcanzables, 6
las mitades `is_finite`, 5 mutan texto de un comentario que quedó dentro de la
condición, y 4 equivalentes (`>=` de la tolerancia, la rama nula, `stopped = true` y
una mutación mía que era la identidad). La primera pasada (606 mutantes) dejó 74
vivos; la mayoría eran guardas que solo cubrían pruebas de otros ficheros, o a las
que les faltaba un lado del límite, y se añadieron casos (linajes distintos de
las copias, tope de media por los dos lados, inicio uno, lista de loci vacía,
procedencia de la referencia, listas de bloque más largas y más cortas que el
recuento, bases de bloque por encima y por debajo, tolerancia con suelo uno,
marcas de tosquedad que no son booleanos, fracción y `m` negativas con fracción
cero) y los tres escenarios que alcanzan los bordes del recorrido, la inserción y la
ploidía 4.

**Ronda 2 (revisión adversarial).** Se corrigió la inalcanzable de las bases no
materializadas (ver arriba) y se anotaron los dos atajos de sesgo como
equivalentes. Cinco mutantes sobre el árbol corregido, huella `e3e9182a69…` antes y
restaurada: el motivo de las bases no materializadas y su guarda neutralizada
**mueren** con la prueba nueva; los tres de los atajos de sesgo (los dos, solo el
de 1, solo el de 0) **sobreviven**, equivalentes. La revisión confirmó además que
los cinco supervivientes «de comentario» eran ruido del generador (repetían las
mitades reales de la condición, que mueren).

**Tanda 7, segunda parte de `results.rs` (11/10/2026): 104 sin disparar menos (de
497 a 393) y 97 de ellas en la configuración efectiva; ninguna guarda borrada, las
comprobaciones siguen en 787.** Las otras 7 son guardas de ayudantes compartidos y
del lector del manifiesto que alcanzan los intentos de llegar a las inalcanzables.
`results.rs` es el juez de los recursos publicados,
de la configuración efectiva y del manifiesto general. Esta tanda toca diez funciones
de la configuración efectiva y sus ayudantes privados exclusivos: 161 guardas sin
disparar al empezar (123 en código vivo y 38 en las dos funciones de la versión 1, que
eran **código muerto**, ver abajo); quedan **26 de las 123 vivas**, todas anotadas
junto a la guarda como inalcanzables con su porqué real. Las de la versión 1 desaparecen
con el código (borrado el 11/10/2026).

| Función | Sin disparar antes | Hoy |
|---|---|---|
| `validate_effective_config_v3_against_manifest_bytes` (la efectiva frente al manifiesto) | 15 | 7 |
| `validate_effective_technical_event_link` (enlace del evento técnico) | 7 | 3 |
| `validate_effective_targets_against_configuration` (la efectiva atada al TOML copiado) | 21 | 11 |
| `validate_effective_initial_populations` | 10 | 0 |
| `validate_effective_direct_events` | 9 | 0 |
| `validate_effective_technical_event` (sin amplificación) | 25 | 1 |
| `validate_effective_reference` | 19 | 2 |
| `validate_effective_derived_sequences` | 9 | 2 |
| ayudantes: `validate_effective_reduction`, `validate_reversible_text` | 8 | 0 |
| **Total vivo** | **123** | **26** |
| `validate_effective_config_v1_against_manifest_bytes`, `validate_effective_targets_against_configuration_v1` (**código muerto, borrado**) | 38 | — |

- *Las dos funciones de la versión 1 son código muerto.* Las dos llevan
  `#[allow(dead_code)]` y **ninguna ruta de producción ni de prueba las llama**:
  `validate_effective_targets_against_configuration_v1` solo la llama
  `validate_published_resources_v1`, que tampoco tiene llamantes, y
  `validate_effective_config_v1_against_manifest_bytes` no tiene ninguno (el único lector de
  una efectiva es `validate_effective_config_against_manifest_bytes`, que solo admite las versiones 4 y 6).
  Lo mismo vale para `validate_published_resources_v1` (21 guardas) y
  `parse_general_manifest_v1_object` (7), fuera del alcance de esta tanda. **No se
  probaron ni se borraron**: decide el responsable si se borran o se anotan. Las 38
  guardas siguen contando en la cifra de la puerta mientras existan.
- *Cubiertas con prueba y motivo exacto (97).* `tests/results_effective_guards.rs`
  (13 pruebas) parte de paquetes **publicados de verdad** por el programa y altera la
  efectiva, o el manifiesto, **en un solo punto**, con el motivo exacto del rechazo: una
  referencia suministrada con dos secuencias, una reducción, dos cohortes (una con
  `founder_keys` y endogamia), tres eventos directos, un evento técnico con
  contaminación biológica y VCF de texto (ploidía global 3 y tolerancia propia, para
  que los valores esperados no coincidan con los de por defecto); una referencia simulada
  con composición y bases desconocidas, con y sin semilla declarada; y una referencia final
  por inserción. El juez de la efectiva frente al manifiesto
  (`validate_effective_config_against_manifest_bytes`) recibe los dos documentos
  alterados; el juez completo, además, lee el TOML copiado. Cada tabla lleva
  control positivo, comprueba que **la alteración alteró** y fija el lado aceptado de
  los límites (ploidía 10, fase 9, longitud de lectura 1 y 150, contaminación de 0 a 0
  y de 100 a 100, tasas 0 y 1, cobertura 0, endogamia 0 y 1, composición con extremos
  cerrados, bases desconocidas iguales a la longitud) y de las listas cerradas. Varias
  guardas **solo se alcanzan con un manifiesto válido distinto**, que las pruebas
  construyen: sin eventos, con la lista de eventos vacía, con dos eventos, de una
  ejecución fallida con el evento fallido, de una fallida sin enlace a la efectiva,
  sin la reducción, sin el VCF y con una reducción que reutiliza el identificador de la
  referencia (el rol de cada FASTA se reajusta para que el lector del manifiesto la
  acepte). Las del TOML copiado se alcanzan **alterando el fichero copiado** (y
  restaurándolo): sintaxis, campo desconocido y objetivo de reducción que no cabe, con
  el código y el significado del catálogo de diagnósticos.
- *Inalcanzables (26).* Ningún paquete publicado y manipulado llega a ellas.
  - *Sombreadas por el lector del manifiesto (11):* manifiesto que no es JSON ni
    objeto (2), versión no admitida (1), una reducción que no es objeto (1), un
    enlace de evento que no es objeto, sin `manifest_artifact_id` o a un artefacto que no
    es manifiesto de evento (3), y el JSON del manifiesto releído por el juez completo,
    el `artifact_id` de la copia TOML, su ausencia en el catálogo y su `original_path` (4);
    `parse_general_manifest_bytes` los rechaza antes, con otro texto.
  - *Sombreadas por la lectura de la efectiva (2):* JSON ilegible y raíz que no es objeto
    (`effective_schema_version`, la única entrada).
  - *Sombreadas por la validación estructural de la efectiva (8):* una reducción o una
    secuencia que no son objeto (3), la raíz de la semilla (1: la valida
    `validate_effective_seed`), y la cabecera y la descripción ausentes (4: la lista
    cerrada de campos exige las claves).
  - *Sombreada por la comparación de origen (1):* el formato FASTA canónico frente
    al origen del TOML; el origen se compara antes.
  - *Imposibles en la práctica (4):* el resumen del TOML solo falla con más de 2^64 bytes;
    la derivación del identificador de ejecución, sobre un registro recién creado; la
    conversión de error de los cambios de inserción (solo corre con plan de inserción, y
    solo falla sin él); y la del plan de identidades técnicas (colisión de identificadores
    derivados con individuos distintos), que **no es imposible en lógica** pero sí en una
    prueba. Intento ejecutado de las sombreadas:
    `unreachable_effective_guards_are_shadowed_by_other_checks`.
- *Mitades equivalentes demostradas (12 mutantes que sobreviven a la batería
  entera).* Ninguna entrada las distingue, y se buscó una para cada una:
  - el filtro de vacío del identificador de individuo inicial (`IndividualId::parse("")`
    ya falla, con el mismo motivo);
  - siete mitades `is_finite` (la endogamia, el mínimo y el máximo del rango de
    contaminación, la tasa de quimeras, la de index hopping, la cobertura y la N de la
    composición): con `serde_json` y `arbitrary_precision`, `as_f64` da `None` para
    `1e999`, así que un decimal no finito no llega a la comprobación; la premisa y cada
    sitio, con el mismo motivo, los fija `unreachable_effective_guards_are_shadowed_by_other_checks`;
  - el techo del mínimo y el suelo del máximo del rango de contaminación (el orden
    `mínimo <= máximo` con `0 <= mínimo` y `máximo <= 100` ya los implica);
  - la mitad de pertenencia del individuo en la comprobación de la biblioteca (el plan
    de identidades solo trae las de los individuos del evento, y un ajeno se rechaza
    antes con «library técnica no derivable»; anotada junto a la guarda);
  - `!sum.is_finite()` de la suma de la composición (las cuatro frecuencias ya están en
    [0, 1]; anotada junto a la guarda).
  La batería entera la corrieron los doce de uno en uno sobre copias; uno
  (`re-n-finite`) «murió» en una ejecución con tres baterías simultáneas por una
  prueba ajena (`current_schemas_close_copy_and_genotype_records_to_unknown_properties`)
  que **no falla** al repetirla sola ni con el mutante aplicado: las tres copias
  compartían ficheros temporales. Repetido solo, sobrevive.
- *Pendientes: ninguno.* Ninguna guarda exigía cambiar código de producción (solo
  comentarios).

**Verificación por mutación.** Sobre tres copias fuera del repositorio (el arnés se
niega a trabajar bajo `/home/ximo/devel`, cada copia con su directorio de compilación y
su copia de seguridad fuera de git), **un mutante por ejecución**, con la huella SHA antes
(`72b840deb3…`), mutada y restaurada de cada uno (coincidieron todas entre «antes» y
«restaurada», con `touch` tras restaurar y comprobación final de que las copias son
idénticas al árbol real). La pasada final fue de **308 ejecuciones** sobre el árbol final: el
motivo de cada rechazo, cada guarda entera, cada mitad de cada disyunción (también la
mitad anulada con `false`), cada operador de comparación y cada extremo cerrado o abierto de
los rangos, cada letra de las listas cerradas (las bases en minúscula, una por una),
y derivaciones (la ploidía global esperada, el origen de la semilla, la tolerancia
comparada por bits, la suma desbordada de longitudes) más 32 motivos de guardas que ya
disparaban otras pruebas. Resultado: **294 mueren**, **2 no compilan** (se reescribieron
y mueren) y **12 sobreviven**, todos equivalentes demostrados (arriba). La primera pasada
dejó 19 vivos: 12 equivalentes y 7 que **no lo eran** (los lados que sobran: más
`founder_keys`, más bibliotecas, más secuencias y más reducciones que el TOML; la semilla del
sistema frente a la declarada; y cada base en minúscula por separado en la inserción y
en la UDI), que se cubrieron con casos nuevos. Una segunda pasada con los dos mutantes
que no compilaban, reescritos, dejó uno más vivo: las bibliotecas del mismo individuo,
cuyo rechazo se escondía tras el `library_id`; se añadió el caso con todo lo demás
correcto. **Veinticuatro de las 32 guardas que ya disparaban otras pruebas tenían un
motivo que ninguna prueba del fichero nuevo afirmaba**, y las tres primeras que se
midieron contra la batería entera sobrevivieron: tampoco lo afirmaba ninguna otra. Se
fijaron con `the_guards_other_tests_only_triggered_have_their_exact_motive` y las 32
mueren.

**Ronda 2 (revisión adversarial).** Se corrigieron dos huecos. (1) La mitad de la
**raíz de la semilla** en la comparación con el TOML **no es equivalente** (la pasada de
arriba la daba por tal): si el TOML declara la semilla, la reconstrucción usa la del TOML,
y el `execution_id` se puede rehacer (`derive_execution_id` y `configuration_digest` son
públicas). Con la raíz cambiada y el identificador rehecho, solo esa mitad rechaza («la
semilla efectiva no coincide con la declaración TOML»), y es la única defensa: con ella
anulada la batería entera sobrevivía y una falsificación completa se aceptaba. Caso nuevo
en dos paquetes (referencia simulada y referencia final por inserción) y comentario
corregido. (2) La mitad `id_catalog` de «catálogo técnico F4 inválido» no la probaba nadie;
caso nuevo (`id_catalog.version = "99"`). Resultado de los dos mutantes (mitad anulada), uno
por ejecución contra la batería entera sobre una copia: **los dos mueren** (huella `10aff0c2a0…` antes y restaurada; mutadas `daf09fd8af…` y `20dd5ec62d…`; los mata `the_effective_config_must_be_the_one_the_toml_resolves` y `the_technical_event_is_judged`). La puerta queda en 721 comprobaciones y 327 sin disparar tras borrar la versión 1. Los equivalentes siguen
siendo 12: los 11 de arriba y `!sum.is_finite()`.

**Alcance de «0 pendientes» en esta tanda.** Se refiere a las guardas del alcance que
cuenta la puerta. DT-040 (catálogo de azar con
tartamudeo clonal), DT-041 (esquema del manifiesto sin perfil con PCR), DT-043, DT-046
y DT-053 no se tocan y ninguna prueba nueva choca con ellas: el evento técnico que se
prueba es el sin amplificación y sin tartamudeo, con el catálogo v3.

**Tanda 8, tercera parte de `results.rs` (12/10/2026): 75 sin disparar menos (de 327
a 252); ninguna guarda borrada, las comprobaciones siguen en 721.** `results.rs` es el
juez de los recursos publicados, de la configuración efectiva y del manifiesto general.
Esta tanda toca las estadísticas iniciales de las poblaciones (lo que el manifiesto
publica sobre cuántas copias de cada alelo salieron y cuántos genotipos de cada clase),
la auditoría de la normalización de frecuencias y la parte del juez físico que relee la
matriz de SNV: 84 guardas sin disparar al empezar; quedan **9**, todas anotadas junto a
la guarda como inalcanzables con su porqué real.

| Función | Sin disparar antes | Hoy |
|---|---|---|
| `validate_initial_population_statistics` (la raíz) | 7 | 1 |
| `validate_initial_population_statistics_provenance` (procedencia y restricciones SNV explícitas) | 7 | 0 |
| `validate_initial_population_statistics_population` | 6 | 0 |
| `validate_initial_population_statistics_locus` | 38 | 0 |
| `validate_serialized_normalization_audit` | 16 | 0 |
| `validate_physical_population_snvs` (el juez físico de la matriz) | 10 | 8 |
| **Total** | **84** | **9** |

- *Cubiertas con prueba y motivo exacto (75).* `tests/results_population_snvs_guards.rs`
  (13 pruebas) parte de un paquete **publicado de verdad** por el programa (referencia
  simulada de dos secuencias, matriz generada de tres loci, tres poblaciones con
  ploidías 2, 2 y 3 y los dos modelos genotípicos) y altera **un solo punto del
  manifiesto**; el juez lee el manifiesto antes que nada, de modo que el rechazo sale del
  lector del manifiesto y, si la guarda no rechazara, el caso cae después en la
  comparación efectiva/manifiesto con **otro** texto (esto es lo que discrimina la
  auditoría, cuyas guardas comparten un único motivo). Los casos «aceptado» se aplican a
  los dos documentos por igual. Cada tabla lleva control positivo y comprueba que la
  alteración alteró. Para las ramas que el generador no produce, la base añade un locus
  **sintético** (tres alelos; o sin una sola copia observada, para la frecuencia
  realizada +0.0/-0.0 comparada por bits) o una población sintética (ploidías 1 y 10);
  lo dice el nombre de cada función. Lado aceptado fijado: ploidías 1 y 10, fase 9,
  coeficiente de endogamia 0 y 1, frecuencias y probabilidades 0 y 1, tolerancias 1e-12
  (aceptado a 5e-13, rechazado a 1.5e-12), tolerancia de la auditoría 0 y 1e-6, posición 1,
  observaciones genotípicas iguales a los individuos, dos copias ausentes por deleción, la
  clave de unicidad de restricciones, loci y grupos difiriendo en una sola componente.
  El juez físico se alcanza con el juez completo: una matriz presentada como
  reutilizada que el TOML genera, y otra que el TOML ni declara (con el TOML copiado
  alterado y el identificador de ejecución rehecho: la derivación es pública), además de
  tres rechazos que antes solo disparaban otras pruebas sin afirmar su texto (el
  enlace de los documentos frente al par publicado, por auditoría y por origen, y el par
  con otro TSV y su identidad rehecha en el catálogo).
- *Inalcanzables (9).* Ningún paquete publicado y manipulado llega a ellas; la prueba
  `unreachable_physical_snv_guards_are_shadowed_by_other_checks` ejecuta el intento y
  **no las discrimina** (otra comprobación rechaza antes, con su texto).
  - *Sombreadas por la lectura de la referencia y de la efectiva (4 + 1 de la
    tolerancia):* el artefacto FASTA ausente de la referencia o del catálogo, en la rama
    de reutilización y en la otra (dos guardas por rama), con y sin referencia final por
    inserción (`read_physical_reference`, `judge_inserted_reference` y la validación de la
    referencia efectiva lo rechazan antes), y la tolerancia de normalización (la efectiva
    ya se ató al TOML por bits).
  - *Sombreadas por el lector del manifiesto (3):* metadatos sin identificador o fuera del
    catálogo y TSV ausente del catálogo.
  - *Sombreada por la validación de la población (1):* el identificador de una población
    estadística (lo exige antes `validate_initial_population_statistics_population`).
  - Quedan también inalcanzables **por llamante** las del ayudante compartido
    `validate_manifest_population_snvs` cuando lo llama la efectiva: llega ya igualada
    al manifiesto, que el mismo ayudante rechazó primero.
- *Pendientes: ninguno.* Ninguna guarda exigía cambiar código de producción (solo
  comentarios).

**Defecto encontrado y cerrado en la ronda 2 (12/10/2026, baja).** El juez de la pareja
efectiva/manifiesto aceptaba una posición de la auditoría escrita `"08"` o `"+8"`, que el
esquema del manifiesto prohíbe (`^[1-9][0-9]*$`): la lectura usaba `str::parse::<u64>`. El
mismo patrón laxo estaba en el lector del TSV de la matriz suministrada
(`prepare_population_snvs`): con `"08"` o `"+8"` **`validate` pasaba** y `run` fallaba después con
un diagnóstico engañoso («efectivo activo ausente»). Por decisión del responsable se arregló
en los dos sitios con una función común, `parse_canonical_position` (`src/population_snvs.rs`):
exige un primer carácter de 1 a 9 antes de leer el número, de modo que rechaza el cero
inicial, el signo y el 0. En el TSV el rechazo es el `InvalidPosition` (GPS-CONF-008) que ya
usaban las posiciones inválidas, y sale en `validate`. Pruebas con motivo exacto para `"08"`,
`"+8"` y `"0"`: la auditoría en `the_normalization_audit_groups_are_judged`; el TSV en
`a_supplied_matrix_position_must_be_canonical_when_validating` (`validate` real, con control
positivo) y en dos filas nuevas de la tabla de `prepare_population_snvs`. **La lectura de la auditoría en
los metadatos (`validate_normalization_audit_document`, la vía de la reutilización) se
arregló también**, con la misma función y su `audit_error("position")`; prueba con motivo
exacto para `"08"`, `"+8"` y `"0"` en `the_audit_judge_rejects_an_identity_that_does_not_fit_the_matrix_or_the_reference`
(el mutante que vuelve al `parse` laxo muere). **Sin tocar, y no son validadores:** las dos
relecturas de un texto que el bucle anterior ya validó con `parse_canonical_position`
(`prepare_population_snvs`, y el orden de grupos, que usa `unwrap_or(u64::MAX)` solo para
ordenar), y el `parse` de la posición de un registro VCF en `src/results.rs`
(`unwrap_or(0)`), que mide bases de un VCF que el propio programa escribió y no juzga
la forma del texto.

**Mitades equivalentes demostradas (32 mutantes que sobreviven, de los que 9 son los de
las inalcanzables).**
  - 13 son mitades `is_finite` (un decimal no finito, `1e999`, lo lee `as_f64` como
    `None` y lo rechaza antes la lectura con su propio motivo; la premisa y cada sitio los
    fija `non_finite_decimals_never_reach_the_finite_halves`).
  - 2 son `>` frente a `>=` en las tolerancias de 1e-12 de las sumas de frecuencias y de
    probabilidades: la diferencia de dos `f64` con 1.0 es múltiplo de 2^-52 (o 2^-53
    por debajo) y 1e-12 no lo es, así que nunca vale exactamente 1e-12.
  - 6 son las mitades de la lista cerrada de claves (longitud y clave desconocida) de
    la restricción explícita, del grupo de la auditoría y de su valor: todas las claves
    son obligatorias y se exigen después con el mismo motivo. **En la raíz, la procedencia
    y la población sí se distinguen** (quitar o renombrar la procedencia, los algoritmos,
    el modelo o los loci cambia el motivo del rechazo posterior): casos añadidos tras la
    primera pasada, que dejó vivos seis mutantes.
  - 1 es `index <= 0` frente a `index == 0` (un índice `usize` no es negativo).
  - 1 es la mitad `metadata_sha256` del enlace del par físico: el lector del manifiesto
    ya exige que sea el resumen del artefacto del catálogo y los bytes se acaban de
    verificar contra él; anotada junto a la guarda.

**Verificación por mutación.** Sobre cinco copias fuera del repositorio (el arnés se
niega a trabajar bajo `/home/ximo/devel`, cada copia con su directorio de compilación y
su copia de seguridad fuera de git), **un mutante por ejecución**, con la huella SHA antes
(`6138a3f136…`), mutada y restaurada de cada uno (coincidieron todas entre «antes» y
«restaurada», con `touch` tras restaurar y comprobación final de que las cinco copias
son idénticas). Primera pasada contra el binario de pruebas nuevo: 422 mutantes distintos
(335 generados: cada motivo, cada guarda entera, cada mitad de cada disyunción, cada
operador y cada literal entero; 78 escritos a mano para las derivaciones que alimentan
guardas —extremos cerrados y abiertos de los rangos, cada ploidía límite, las
tolerancias, los desbordes con `wrapping_*`, las componentes de las claves de unicidad,
el orden de grupos y valores, el origen y la existencia de la fuente física—; 9 más
reescritos para los que no compilaban). Resultado: **376 mueren**, **14 no compilan**
(los 9 que importaban se reescribieron y mueren; los otros 5 son mitades `is_finite`
mal formadas) y **32 sobreviven**, todos equivalentes o inalcanzables (arriba). La
primera pasada dejó vivos 7 que **no lo eran** (las dos mitades de la lista de claves de
la raíz, de la procedencia y de la población; la posición 1 de la auditoría, que un
literal `> 1` rechazaba; y la suma desbordada `observadas + ausentes`), que se cubrieron
con casos nuevos y mueren. Los 14 mutantes del juez físico que sobrevivían al binario
nuevo se corrieron contra la **batería entera**, de uno en uno sobre copias con su propio
directorio temporal: 9 mueren por pruebas ajenas, 3 mueren con casos nuevos (los textos
del par físico y el origen del enlace, que nadie afirmaba) y 2 sobreviven (la tolerancia
inalcanzable y la mitad `metadata_sha256`).

**Ronda 2 (revisión adversarial, 12/10/2026).** (1) Se quitó el caso «posición con un cero
a la izquierda» del juez físico: exigía el rechazo del enlace, que solo salía porque el
juez de la pareja **aceptaba** el `"08"` (ver arriba el defecto, ya cerrado); la mitad de la
auditoría que discriminaba la sigue discriminando «otra decisión de normalización». (2) Se
cerró el `"08"` y el `"+8"` en la auditoría y en el TSV suministrado. (3) **Los documentos
sintéticos y varios casos «aceptados» no cumplían el contrato**: la arquitectura fija que los
genotipos realizados cuentan solo observaciones completas, y el locus de tres alelos, la
población sintética, «dos copias ausentes por deleción» y «tantas observaciones como
individuos» las incumplían aunque el docstring dijera lo contrario. Ahora la población
sintética saca los recuentos de los alelos de sus genotipos, el locus de tres alelos cuadra
(A = 2, G = 3, C = 1 con un individuo A/G, otro A/C y otro G/G), los casos aceptados fijan
los genotipos que los cumplen y `assert_architecture_inequalities` lo comprueba sobre la
base y sobre cada caso aceptado. El caso «tantas observaciones como individuos» se quitó: el
paquete ya trae el límite y lo cubre el control positivo (la prueba lo afirma). El juez de
la pareja **no** comprueba esas dos desigualdades; por decisión del responsable no se
añaden y se registran como DT-055. (4) Los límites, que se probaban lejos del borde, se
prueban **justo fuera**: `next_up(1.0)`, `next_up(1e-6)`, el negativo más pequeño
representable (`-5e-324`) y la primera desviación representable por encima de 1e-12 (4504
epsilon, frente a 4503 que se acepta), en la tolerancia y la frecuencia de referencia de la
auditoría, en los valores de la auditoría, en las frecuencias objetivo y realizada del
locus, en la probabilidad y el recuento esperados, en las dos sumas y en el coeficiente de
endogamia. (5) `publish_files` crea ahora la guarda que borra el directorio antes de
crearlo. (6) Los comentarios de las mitades `is_finite` se movieron al brazo `f_diploid` y se
corrigieron: NaN y los infinitos nunca están en `0.0..=1.0`; en el filtro del recuento
esperado (que sí aceptaría `+inf`) solo vale la premisa de la lectura. (7) Registrada en
DT-045 la simplificación de leer la tolerancia de la configuración reconstruida.

**Verificación por mutación de la ronda 2** (tres copias fuera del repositorio, un mutante
por ejecución, huella antes `3583e33589…` en `results.rs` y `4b138c2fe4…` en
`population_snvs.rs`, igual restaurada; `touch` tras restaurar): 30 mutantes de límites
(cada rango `0.0..=1.0` por arriba y por abajo, `0.0..=1e-6`, `>= 0.0`, las dos tolerancias
de 1e-12 en 1.4e-12 y 1.9e-12; los supervivientes son solo los de las líneas de
comentario, que el generador también toca) más la vuelta al `parse` laxo en los dos sitios
y cuatro del ayudante (`b'0'..=b'9'`, `b'1'..=b'8'`, admitir `+`, aceptarlo todo): **mueren
todos los que tocan código**. Se repitieron los 77 escritos a mano y 5 reescritos de la
primera ronda contra las pruebas nuevas: mueren todos salvo los equivalentes ya demostrados.

**Alcance de «0 pendientes» en esta tanda.** Se refiere a las guardas del alcance que
cuenta la puerta. DT-045 (las catorce guardas muertas de la reutilización de la matriz en
`src/population_snvs.rs`), DT-049 y DT-053 no se tocan y ninguna prueba nueva choca con
ellas: la reutilización de la matriz se prueba solo hasta la presentación como
reutilizada.

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

### Cómo se cerró (06/10/2026)

**Se borraron las dos**, por decisión del responsable, y con ellas el único
lector del ayudante de prueba `FORCE_FASTQ_REVALIDATION_FAILURE`, que quedó
huérfano y se borró también. Las pruebas ignoradas del proyecto pasan de once a
nueve.

Antes de borrarlas se comprobó, **mutando una copia** del árbol, qué parte de lo
que decían cubrir tiene hoy una prueba viva que pueda fallar:

| Propiedad de las pruebas vaciadas | Cobertura viva | Mutación en la copia |
|---|---|---|
| Un fallo conjunto de R1 y R2 al publicar el FASTQ descarta los dos mates y no deja nada en el catálogo | `fastq_output::tests::paired_publication_rolls_back_on_injected_r2_or_joint_validation_failure` | quitar el descarte de R1: la prueba cae; sustituir la inyección del juez conjunto por `Ok(())`: la prueba cae |
| Con el evento ya publicado, el paquete fallido conserva la referencia, su mapa y el manifiesto del evento | `execution::tests::a_failed_run_with_reads_keeps_the_final_reference_with_its_map` | no retener el manifiesto del evento, no retener los FASTQ o conservar todo el catálogo: la prueba cae las tres veces |
| Las reducciones de un paquete fallido quedan en estado `failed` sin identificadores de artefacto, y sus ficheros no se publican | `execution::tests::physical_resource_failure_marks_reduction_failed_without_partial_artifacts` | quitar la línea que las marca `failed`: la prueba cae |

**Dos propiedades no tenían equivalente en la vía actual** (mutadas, la
batería entera seguía en verde), y se cubrieron el mismo día con pruebas nuevas
de la vía de parejas, no restaurando las de lecturas sueltas:

1. **Fallo con referencias reducidas ya completas.** Al cubrirlo se vio que el
   código **descartaba** la reducción completa y validada, contra REQ-SAL-25 (un
   fallo no invalida retroactivamente lo ya publicado y validado que sea
   independiente) y DTC-045 (el fallo controlado conserva y declara
   exportables esos artefactos); REQ-SAL-27 solo obliga a marcar fallida la
   reducción cuando falla **ella misma**. Por decisión del responsable se
   **corrigió el código con una regla única**: toda reducción completa cuyos
   bytes en disco siguen siendo los que se validaron al escribirlos (FASTA y
   mapa, comparados con el registro de cuando se catalogaron) se conserva
   entera —FASTA, mapa, enlaces y estado `complete`—, haya evento o no y falle
   lo que falle; la que no pasa esa comprobación, o la que causó el fallo, queda
   `failed`, sin enlaces y sin ficheros (REQ-SAL-27). Sin evento se conserva
   además la referencia suministrada de la que deriva. Una primera versión de la
   corrección distinguía «con evento», y una revisión adversarial encontró que
   publicaba como completa una reducción que el juez final rechazaba (FASTA con
   una base cambiada, mapa con un campo ajeno), porque el manifiesto fallido no
   lleva configuración efectiva y ese juez ya no volvía a correr. Pruebas, todas
   en `src/execution.rs`:
   - `a_failure_after_the_event_keeps_the_completed_reduced_references` y
     `a_failure_without_event_keeps_the_completed_reduced_references_and_their_source`
     (dos reducciones, ambas deben conservarse);
   - `a_reduction_altered_after_validation_is_discarded_and_the_other_is_kept`
     (con y sin evento, alterando el FASTA o el mapa: solo la alterada queda
     fallida);
   - `physical_resource_failure_marks_reduction_failed_without_partial_artifacts`
     (la reducción misma falla) sigue pasando.
   El manifiesto fallido sin evento y con reducción completa lo aceptan el juez
   del manifiesto y el esquema `manifest-general-v4` (jsonschema 4.10.3), con y
   sin evento.

   **Limitación conocida (no es deuda nueva).** La comprobación es que los bytes
   en disco son los que se validaron al escribirlos; **no** equivale a repetir el
   juicio final. Al escribir, `validate_reduction_map_document` compara el mapa
   con la derivada que construye el propio productor; el juez final
   (`validate_published_resources`) además reconstruye la derivada desde la
   subsemilla y la selección y la compara base a base. Si el productor tuviera un
   defecto en la selección, el juez final rechazaría la reducción con los bytes
   intactos, la comprobación de huellas pasaría y se conservaría como completa
   aunque fuera la causa del fallo. Exige un defecto del productor, no una
   manipulación, y queda fuera de lo decidido.

   Mutaciones (SHA de `src/execution.rs` `68c26a929762` antes y
   restaurado; entre paréntesis el mutado y lo que cae):
   - conservar sin revalidar el disco (`b97471fc6588`): la prueba de la
     reducción alterada y la de fallo físico;
   - nunca conservar (`76edda1e4481`): las cuatro pruebas de conservación;
   - conservar solo la reducción de índice 0 (`0b2d3bc18d5c`): las cuatro;
   - no comprobar el FASTA (`e2e7d84552d5`) o no comprobar el mapa
     (`c40ab2c2150a`): la de la reducción alterada (y, en el primer caso, la de
     fallo físico);
   - sin evento, no conservar la suministrada (`21a831263125`): tres;
   - sacar la reducción de `references` (`0a1bce0be540`) o no retener su mapa
     (`014e13dbdcba`): cuatro, pero **por el juez del catálogo** en la primera
     aserción (`GPS-FMT-006`), no por las aserciones de lista o de disco, que
     quedan implicadas por él (dicho en el docstring de la prueba);
   - sustituir la causa registrada (`888b53f0aa13`): cae la prueba de fallo tras
     el evento, que usa un código propio.
   La mutación original de esta deuda (la primera versión, que no excluía la
   derivada) caía **primero por el juez del manifiesto**, no por una aserción
   propia; por eso se rehízo la prueba.
2. **Fallo del ensamblado del evento, incluida la validación del FASTQ.**
   `execution::tests::a_failed_event_assembly_keeps_the_completed_reductions_but_no_fastq_nor_event`
   (REQ-SAL-22, REQ-SAL-23 y REQ-SAL-25): con dos individuos, la validación
   conjunta falla en la **segunda** biblioteca, cuando la primera ya cerró sus
   FASTQ; el paquete fallido no conserva ningún FASTQ ni fichero de evento ni
   artefactos parciales (el catálogo es exactamente la configuración, la
   referencia suministrada y las dos reducciones con sus mapas), y registra la
   causa. **No exige `events: null`**: que el evento fallido no quede enlazado
   contradice REQ-SAL-21, y está en DT-036. El gancho de prueba de
   `src/fastq_output.rs` (solo `#[cfg(test)]`) pasó de booleano a cuenta atrás
   para poder fallar la n-ésima biblioteca. Mutaciones: la guarda
   `if events.is_some()` a `if true` (`68c26a929762` a `47c507f623ff`), cae; un
   gancho que nunca dispara (`fastq_output.rs` `2df62d66de71` a `67ffdb7d82e7`),
   cae; `retain_only_artifacts` sin borrar los mapas de reducción
   (`publication.rs` `816eba15b479` a `aebc945a9785`), cae (con la versión
   anterior de la prueba; no se repitió sobre la actual).

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

**Actualización (06/10/2026):** las dos de `src/execution.rs` se borraron al
cerrar DT-002. Quedan nueve ignoradas en total, de las cuales **seis** son de
contratos obsoletos (una en `src/publication.rs`, una en
`tests/f3_event_assembly.rs`, cuatro en `tests/f3_event_admission.rs`).

---

## DT-004. La puerta de Python no corre desde el 25 de agosto de 2026

**Gravedad: media.** `tests/test_process_adapter.py` importa `pytest`, que no
está declarado en `pyproject.toml`. La puerta canónica
(`unittest discover -s tests_ui`) **nunca lo recogía**, porque el fichero estaba
en `tests/`, fuera de lo que recorre; solo corría invocándolo a mano, y entonces
abortaba al importar `pytest` si no estaba instalado.

Viene del commit `6139fe9`. Registrada también en el plan de la fase 5; se
recoge aquí porque es transversal. Quien la cierre debe decidir entre declarar
`pytest` como dependencia de desarrollo o reescribir el módulo con `unittest`.

### Cómo se cerró (06/10/2026)

**Decisión del responsable: reescribirlo con `unittest`, sin añadir
dependencias** (`pyproject.toml` y `uv.lock` no se tocan). El módulo tenía una
sola prueba, la que exige que el adaptador rechace un binario que declara la
versión 0.1.0 antes de ejecutar ninguna acción del formulario. Se pasó, con sus
mismas aserciones, a `tests_ui/test_process_adapter.py` como
`test_e147_f1_version_is_rejected_before_any_form_action` —`tmp_path` es ahora
`tempfile.TemporaryDirectory`, `pytest.raises(match=...)` es
`assertRaisesRegex`— y se borró `tests/test_process_adapter.py`. Con eso la
descubre la puerta canónica, `uv run --frozen python -m unittest discover -s
tests_ui`, **sin ningún comando nuevo**: el fichero estaba en `tests/`, fuera
de lo que esa puerta recorre, y además importaba `pytest`, así que ni
`discover` ni `unittest` podían ejecutarlo.

Comprobado mutando el adaptador en una copia: si la comparación de versión
acepta cualquiera, la prueba nueva cae (y con ella las dos de
`test_version_check_accepts_only_the_supported_literal`); si solo compara el
nombre y no la versión, vuelve a caer.

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

### Estado: aceptada (06/10/2026)

**Aceptada, sigue abierta**, por decisión del responsable: el muestreador se
deja como está y no hay trabajo previsto. El tope de intentos convierte lo que
parecía un cuelgue en un diagnóstico claro, y el resultado nunca es incorrecto.
**Lo que se paga** es doble: en los casos extremos, una ejecución rechazada; y,
antes de llegar a ellos, un **régimen lento pero admitido** —del orden de
23 000 intentos por sorteo con una copia de 271 pb, y hasta un millón admitidos
por sorteo—, en el que la ejecución termina, pero tarda. Arreglarlo de fondo
cambiaría cada byte publicado y exigiría rehacer toda la evidencia de
congelación, un coste que hoy no compensa. **Se reconsidera si** una fase puede
permitirse reescribir la línea base de azar de la fase de realismo técnico; el
método sin rechazo sigue registrado como ampliación futura DTF-060. No se tocó
código.

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

### Cómo se cerró (06/10/2026)

El responsable decidió que **manda Gymrek por omisión y el multiplicador de los
homopolímeros es configurable** (DTC-099, REQ-CONF-141). Se ajusta sobre la
columna sin PCR, que no cambia, de modo que ninguna ejecución sin PCR se mueve
y una amplificada sin el campo publica los mismos bytes. Quien quiera la cifra
de Fungtammasan escribe `homopolymer_pcr_fold = 2.5`; su propia figura da 2,82
a 15 pb, el ancla, y el usuario elige.

**Implementado el 08/10/2026.** La igualdad byte a byte sin el campo se
comprobó de dos maneras: cinco escenarios ejecutados con el binario del commit
de partida y con el nuevo dan los mismos 55 ficheros
(`evidence/dt-009-multiplicador-congelacion/`), y una huella congelada medida
en el commit de partida antes de tocar el código vive ahora en
`tests/clonal_homopolymer_fold.rs`. Al implementarlo se encontró además que el
juez no comparaba el bloque del tartamudeo clonal del manifiesto de evento con
el de la configuración efectiva; se añadió esa igualdad (REQ-VAL-61). Con ella
y con la comparación con el TOML copiado que ya existía, el responsable aprobó
el 08/10/2026 que el juez **no** lleve las guardas campo a campo que proponía
el estudio previo (DTC-099). La igualdad cubre el bloque del tartamudeo clonal
entero; el resto de la sección de amplificación del manifiesto de evento sigue
sin juez propio (DT-043).

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

### Ampliación del 04/10/2026

El juez del abanico de deslizamiento (`validate_published_slippage_fans`,
`reject`) también cita `effective-config-v4.schema.json` sobre una
configuración efectiva v6: el corpus de la revisión adversarial del 04/10/2026
produjo **630 diagnósticos** así, idénticos en la base `aef5ed4` y en el árbol
corregido (no es una regresión). **Ninguna prueba fija el nombre del esquema de
ese juez.** Encaja aquí: es el mismo mecanismo, un rechazo de la configuración
efectiva que no se reetiqueta.

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

**Gravedad: baja. Aceptada el 04/10/2026**, con una condición para
reconsiderarla.

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

---

## DT-024. El juez del abanico no termina con un número de linajes enorme

**Gravedad: media.** Un juez que cuelga en vez de rechazar no protege a nadie.

Lo encontró el corpus de corrupciones de la revisión adversarial del
04/10/2026. Con `lineages = 18446744073709551615` en la configuración efectiva
de un paquete publicado, `validate_published_resources` no termina: con
`timeout 30` sale con el código 124. **Igual en la base `aef5ed4`**: no es una
regresión. No hay ningún límite sobre `lineages` antes de `recompute_slippage_fan`,
que recorre un linaje tras otro. Haría falta rechazar, antes de rehacer el
abanico, un número de linajes mayor que lo que el generador admite o que el
presupuesto de recursos declarado. No se ha corregido: es una deuda registrada.

### Cómo se cerró (06/10/2026)

**Corrección.** Antes de rehacer el abanico, `validate_published_slippage_fans`
(`src/results.rs`) exige dos cotas, cada una con su rechazo propio:

- **Los linajes han de ser las copias haplotípicas de las poblaciones
  iniciales publicadas**: individuos por ploidía, sumados con aritmética
  comprobada, que es el recuento que usa el productor
  (`haplotype_copy_count` en `src/application.rs`). No es una cota nueva: el
  generador no admite otro número, y el documento publica además la lista de
  los individuos, de modo que el valor queda acotado por lo publicado.
- **La media de pasos no pasa de `SLIPPAGE_MAX_MEAN_EVENTS`** (1000), el tope
  que el productor ya impone. Es el otro número del documento que gobierna el
  coste del recorrido (un paso por linaje y por evento); con el máximo del
  muestreador de Poisson (2^52) el juez tampoco habría terminado.

En el mismo juez no queda ningún otro bucle ni reserva gobernados por un número
leído sin acotar: `max_units` solo se compara, y los alelos y los motivos
recorren listas del propio documento.

**Prueba.** `the_fan_judge_rejects_unsupportable_lineages_and_mean_events_without_hanging`
(`tests/repeat_slippage.rs`). Va primero el caso de dos millones de linajes,
que sin la cota **termina** (37 s en depuración) y rechaza por otro motivo, de
modo que la prueba muere por el motivo y no cuelga la batería; después, el de
`u64::MAX`. Con la media pasa lo mismo, con 10 000.

---

## DT-025. El juez de recursos publicados aborta con `capacity overflow`

**Gravedad: media.** Aborta el proceso en vez de devolver un rechazo.

Lo encontró el mismo corpus. Con `reference_units[1] = 18446744073709551615`
en la configuración efectiva, `validate_published_resources` aborta con el
pánico `capacity overflow` en vez de rechazar el documento. **Igual en la base
`aef5ed4`**: no es una regresión. Falta validar el valor antes de reservar con
él. No se ha corregido: es una deuda registrada.

### Cómo se cerró (06/10/2026)

**Corrección**, en `published_simulated_tracts` (`src/results.rs`), el lector de
los tractos de la referencia simulada que publica la configuración efectiva:

- **Cada bloque de `reference_units`** se calcula con `checked_mul` y se compara
  con el hueco que queda en el tramo (`length - bases.len()`, que no desborda
  porque ningún bloque aceptado pasa de ahí), **antes** de reservar. Antes
  `units.saturating_mul(..)` y una suma sin comprobar dejaban pasar un valor
  imposible hasta `repeat_motif_over`, que reserva tantos bytes como le piden.
- **El final del tramo** ha de caber en la secuencia de la referencia publicada.
  Es la misma reserva por la otra vía: un `end` enorme en un locus de un solo
  motivo pedía `end - start + 1` bytes. No estaba en el encargo; se encontró
  revisando las demás reservas del juez.

- **El juez de la referencia final por inserción**
  (`judge_inserted_reference`) acotaba la suma de `reference_units` solo con
  `max_units`, que se lee del mismo documento y no acota nada: con
  `reference_units = [6148914691236517205]` y `max_units` igual reservaba 8 GiB
  («memory allocation of 8589934592 bytes failed», salida 134), y es alcanzable
  desde un paquete ajeno al resolver una referencia desde su manifiesto. Lo
  encontró la revisión adversarial. Ahora, antes del bucle que construye el
  tracto, calcula Σ unidades × longitud del motivo con `checked_mul` y
  `checked_add`, **acumulada entre loci**, y rechaza si pasa del tamaño en
  bytes del FASTA final publicado, que contiene todo lo insertado. La cota
  sale, pues, de un fichero y no del documento.

**Reservas del juez revisadas** (todo `src/results.rs` fuera de pruebas: cada
`with_capacity`, `vec![..; n]`, `repeat`, `resize` y bucle `0..n`): el resto se
dimensiona con la longitud de listas del propio documento, con datos físicos ya
leídos o con las cotas de arriba; la única otra reserva con un número del
documento, el largo de las secuencias que pasa al regenerador de la referencia
simulada, lo ata antes a la declaración TOML copiada en el paquete (se comprobó con
`length` y `total_genome_size` de un terabase: «la referencia efectiva no
coincide con la declaración TOML»).

**Pruebas**, cuatro casos nuevos en
`every_corruption_of_a_published_tract_is_rejected_by_its_own_branch`
(`tests/simulated_reference_repeat_tract.rs`): el de `u64::MAX` del hallazgo y
tres elegidos para que, sin su cota, el juez **termine** y rechace por otro
motivo en vez de pedir terabytes —un bloque de diez millones de unidades, un
producto que da la vuelta a `u64` y cae justo en el hueco, y un final de veinte
millones en una secuencia de mil bases—; una prueba aparte para que el final
tenga que caber en **su** secuencia y no en otra más larga; y dos casos en
`every_corruption_of_an_inserted_publication_is_rejected_by_its_own_branch`
(`tests/reference_tract_insertion_run.rs`) para el juez de la inserción: tres
millones de unidades, que sin la cota terminan y rechazan por otro motivo, y
`u64::MAX / 3`.

---

## DT-026. El juez de cada mate del FASTQ final no leía ningún registro

**Gravedad: grave. Estado: cerrada el 5 de octubre de 2026, el mismo día en que
se anotó.**

El juez que debía decidir si cada mate del FASTQ final publicado era correcto
daba por bueno cualquier contenido.

Lo encontró el corpus de la revisión adversarial del 05/10/2026, al comparar la
ampliación de la puerta con la base `a4beae1`, y era anterior a ese cambio. En
`validate_final_fastq_mate` (`src/fastq_validation.rs`, entonces líneas
234-242), `file.try_clone()` **compartía la posición de lectura** con el
fichero original; `validate_canonical_bgzf_eof` la dejaba al final del fichero,
y el lector BGZF que venía después leía cero líneas y devolvía `Ok`. Por tanto
**toda `validate_final_record`** —la estructura de cuatro líneas, el alfabeto
de las bases y el de las calidades— **era código muerto en producción**.

*Reproducción, sobre el código de entonces.* `printf 'hello\nworld\n+\nxx\n' |
bgzip -c > a.fq.gz` y `validate_final_fastq_mate` sobre ese fichero devolvía
`Ok`. De los 268 FASTQ corruptos del corpus pasaban 260, y el mutante
«`validate_final_record` devuelve `Ok`» sobrevivía a la batería entera.

### Cómo se cerró

**Corrección.** Se abre el fichero una sola vez, se comprueba el marcador de
fin sobre ese descriptor y se **rebobina** al principio antes de leerlo, igual
que ya hacía el juez del FASTQ WGS de F3 (`validate_wgs_fastq`). Desaparece el
`try_clone`; no queda ninguno en `src/`.

**Una sola pasada del juez de cada mate.** Encendido el juez, cada mate se
juzgaba dos veces enteras al publicar: una en el escritor
(`src/fastq_output.rs`, al catalogar cada mate) y otra al empezar el juez de la
pareja (`validate_final_paired_fastq`). Se quita la del escritor. Queda la del
juez de la pareja, que es la única ruta por la que pasa toda publicación —el
escritor por flujo y el de colección acaban en el mismo cierre— y que, si algo
falla, descarta los dos artefactos. Ningún contrato se rompe: la comprobación
sigue ocurriendo antes de la publicación atómica (REQ-DAT-17). El juez de la
pareja sigue releyendo después las cabeceras de los dos mates para
compararlas; esa segunda lectura no se ha fundido con la primera.

**Lo más importante era que el juez encendido no rechazara nada válido**,
porque su código no se había ejecutado nunca fuera de las pruebas. Se comprobó
de tres formas:

- **la batería entera**, que publica FASTQ finales en decenas de escenarios y
  ahora los pasa todos por el juez de registros: con el rebobinado solo, 1670
  pasan, 0 fallan y 11 ignoradas, las mismas cifras que sin él; con todo el
  cambio, 1690, 0 y 11;
- **una exploración de 2240 escenarios con lecturas** (240 + 1000 + 1000,
  generados al azar con tres semillas fijas), ejecutados con el binario de la
  base `2eaced4` (SHA-256 `9e0c07e1…802bf0`) y con el del árbol final
  (`e5a38618…f9e684`). Combinan referencia simulada y suministrada —esta con
  minúsculas y tiradas de `N`—, microsatélites de tracto existente, insertados
  en la referencia, colocados en la simulada y con deslizamiento germinal,
  referencia reducida, tartamudeo técnico, PCR con y sin tartamudeo clonal,
  quimeras, contaminación, saltos de índice, errores de secuenciación
  encendidos y apagados, y longitudes de lectura de 1 a 150. **Ninguno cambia
  de resultado**: los 2201 que publican en la base publican con el juez
  encendido, con los FASTQ **idénticos byte a byte**, y los 39 que la
  configuración rechaza se rechazan igual en los dos. Se midió dos veces: con
  el rebobinado solo (binario `79125208…cc90f9`) y con el árbol final. En las
  dos tandas de mil, 1613 escenarios publican al menos una pareja (5 485 685
  parejas); el resto pide cobertura cero. La revisión adversarial añadió 804
  escenarios propios, con lecturas extremas, sin ningún cambio de resultado,
  y 46 casos BGZF límite, cada uno con el veredicto que el juez debe dar hoy
  (dos son el pánico de DT-030; nueve, huecos de DT-029; y uno, una calidad
  Q93 por encima del perfil, la admite el contrato, que con un desplazamiento
  de calidad permite cualquier Phred de 0 a 93);
- **dos jueces externos en Python** sobre las dos tandas de mil: uno con la
  cabecera de REQ-DAT-68, el alfabeto `ACGTN`, el separador y las calidades
  Sanger (8068 FASTQ, 10 971 370 registros), y el juez estricto de la
  revisión, que mira además las siete propiedades de DT-029 (los 2000
  paquetes en orden: 5 485 685 parejas y 1 811 988 `N`).
  Ningún registro incumple el contrato. El escritor era correcto; lo que
  faltaba era el juez.

**Pruebas.** En `tests/final_fastq_record_validation.rs`, diecinueve:

- una por cada clase de rechazo del juez del mate y del analizador de su
  cabecera, que exige el **motivo exacto** —varias guardas se solapan y sin el
  motivo una guarda borrada quedaría tapada por la siguiente—, incluido el
  caso de la reproducción y el de unos registros válidos **sin el marcador de
  fin**, que es la única comprobación del marcador en este camino;
- un control positivo con los dos extremos de las calidades (`!` y `~`), y un
  FASTQ publicado de verdad, con `N` de plantilla y errores de secuenciación,
  que se acepta, mientras que la misma copia con un byte de calidad cambiado
  se rechaza;
- cinco del juez de la pareja: un control positivo, dos mates que nombran
  parejas distintas, un artefacto con el SHA-256 o el tamaño alterados, y una
  pareja cuyas cabeceras y artefactos coinciden pero con un registro inválido
  en R1, que solo el juicio de los mates puede rechazar.

Y en `tests/f4_final_fastq_provenance.rs`, una más: el escritor rechaza y
descarta una pareja cuyos registros incumplen el contrato, sin dejar ningún
mate ni ninguna entrada en el catálogo. Es la que impide quitar la llamada del
escritor al juez de la pareja, que desde este cambio es la única pasada.

Contra el código de la base `2eaced4`, extraída con `git archive`, fallaban
trece de las catorce pruebas de la primera ronda —todas menos el control
positivo—.

**Mutantes**, sobre una copia del árbol, con SHA-256 antes, mutado y
restaurado: veintinueve. Mueren veintiséis:

- quitar el rebobinado, `validate_final_record` devolviendo `Ok`, cada guarda
  de estructura y de alfabeto, sacar la `N` del alfabeto, mover un extremo del
  intervalo de calidades y cada guarda del analizador de la cabecera;
- de la revisión: quitar la comprobación del marcador de fin, que el juez de la
  pareja no juzgue los mates, que no compare las cabeceras de R1 y R2, que no
  compare el SHA-256 publicado; y además que no compare el tamaño y que el
  escritor no llame al juez de la pareja.

**Sobreviven tres, y son inalcanzables desde un fichero**: las comprobaciones
de que la cabecera, la secuencia y las calidades terminan en salto de línea.
El lector de líneas rechaza antes toda línea final sin salto y el bucle
rechaza antes todo registro con una línea vacía; se dice así en la cabecera
del fichero de pruebas.

**Corpus de corrupciones.** El de la revisión ya no estaba disponible al
cerrar la primera ronda, así que se rehízo: 334 corrupciones de un solo cambio
sobre dos mates publicados de verdad, en el primer registro, en uno
intermedio y en el último, más los dos mates intactos. Con la base se aceptan
los 336. Con el juez encendido se rechazan 328 y se aceptan los dos controles
y seis corrupciones, las que cambian `TP_` por `XX_` (DT-029).

**Puerta de cobertura.** `src/fastq_validation.rs` pasa de 57 a **36**
comprobaciones sin disparar, con las mismas 71. El techo baja a 36.

**Lo que el juez sigue sin mirar** no vive aquí: son siete propiedades del
contrato, reunidas en DT-029, cerrada después. La revisión encontró además que un
bloque BGZF con un tamaño declarado imposible hace abortar al juez en vez de
rechazar el fichero (DT-030), y que una referencia corta con muchas quimeras
aborta como defecto interno (DT-031). Ninguna de las tres es parte de este
defecto.

---

## DT-027. La verdad de un paquete sin poblaciones no se juzga

**Gravedad: media.**

También anterior al cambio. En `src/results.rs` (líneas 4171-4178), un paquete
sin poblaciones no pasa por el juez completo de la verdad. El comentario dice
que `validate_published_reference_names` juzga la cabecera, y esa función solo
mira la referencia. Resultado medido por el corpus: **31 alteraciones de la
cabecera de `truth.jsonl`** que el lector rechaza **pasan** el juez de los
recursos publicados. No se ha corregido.


### Cómo se cerró (06/10/2026)

**Corrección, en dos partes.**

1. `validate_physical_truth` (`src/results.rs`) ya no se salta la verdad
   cuando no hay poblaciones: el lector de la verdad (`read_truth_jsonl`)
   la lee siempre y aplica el contrato de la cabecera. El juez semántico
   completo sigue corriendo solo cuando hay poblaciones, porque solo entonces
   hay copias y casos que juzgar. Se midió con el corpus de la revisión
   (105 alteraciones de la cabecera, con referencia simulada y con
   suministrada) que correrlo también con bases y sin poblaciones da
   exactamente los mismos veredictos: no aporta nada medible, y por eso no se
   deja.
2. **La huella del contenido de la referencia ya no se omite sin bases.**
   `validate_published_reference_names` solo la comparaba si podía
   calcularla sobre las bases materializadas, y con una referencia suministrada
   leída solo con su ficha la omitía en silencio, contra DTC-061 («los
   digests discordantes fallan»). Ahora usa la huella que el lector del FASTA
   publicado calcula en la misma pasada que lo valida, sin retenerlo en
   memoria (`ValidatedFastaSummary::content_sha256`), y si no hay ninguna de
   las dos, o la cabecera no trae la misma, rechaza. Esa pasada la hace el
   propio juez sobre el FASTA publicado, así que sigue siendo independiente de
   quien escribió la verdad. Los dos cálculos dan la misma huella:
   `reference_fasta_summary_matches_full_read` y
   `summary_content_digest_equals_the_materialized_one_on_wrapped_sequences`
   (`tests/fasta_validation.rs`) lo comprueban con minúsculas, `N`, CRLF y
   secuencias en varias líneas de ancho desigual.

**Prueba.** `truth_header_of_a_package_without_populations_is_judged`
(`tests/published_resource_validation.rs`): trece alteraciones, con la
referencia simulada y con la suministrada (esta última sin bases). Nueve
alteraciones de la cabecera y del documento (tipo de documento, versión,
versión promovida, convención de coordenadas, clave desconocida, recuento de
copias, huella de grupo mal formada, registro de tipo desconocido, segunda
cabecera) pasaban el juez antes; dos casos de huella del contenido (a ceros y
con un dígito cambiado) se rechazaban con bases y se aceptaban sin ellas; dos
controles ya los cazaba la comparación de nombres. Se verificó fallando contra
tres variantes: restaurar el filtro por poblaciones (mueren las 18 alteraciones
de la cabecera en las dos referencias), saltarse el lector solo sin bases
(mueren las 9 con la suministrada) y volver a omitir la huella sin bases
(mueren las 2 de huella con la suministrada). Lo que la prueba no discrimina:
el juez semántico con bases y sin poblaciones, que por lo dicho no cambia
ningún veredicto.

**Lo que sigue abierto, y no es de esta deuda.** En el corpus de la revisión
pasan el juez, también con poblaciones, un `execution_id` cambiado en la
cabecera y un `source_digests` de más; no hay norma que obligue a que el
`execution_id` de la verdad coincida con el del manifiesto. Queda anotado en
DT-032.

---

## DT-028. Sospecha: un alelo de referencia admite claves desconocidas

**Gravedad: baja, sin confirmar del todo.**

Un alelo `{"kind":"reference","zz_extra":1}` se acepta, mientras que el alelo
`snp` rechaza las claves desconocidas. Lo anotó la revisión adversarial del
05/10/2026 como sospecha; no se ha investigado si el contrato de la verdad
quiere ese rechazo también para el alelo de referencia.


### Cómo se cerró (06/10/2026)

**Reproducida.** El esquema de la verdad (`schemas/truth-v5.schema.json`, y los
anteriores) cierra el alelo de referencia con `additionalProperties: false`,
igual que el resto, así que el contrato es inequívoco. El lector aceptaba
`{"kind":"reference","zz_extra":1}` porque una variante sin campos de un enum
etiquetado por dentro no aplica `deny_unknown_fields` de serde. Es la única
variante de ese tipo en el lector de la verdad.

**Corrección.** `WireAllele::Reference` pasa a variante con llaves vacías
(`Reference {}`) en `src/truth_validation.rs`.

**Prueba.** `schema_reader_rejects_unknown_keys_in_every_allele_kind`
(`tests/truth_validation.rs`): una clave de más en un alelo de referencia, uno
SNP y uno de inserción, cada uno con su control (el fichero intacto se lee).

---

## DT-029. El juez del FASTQ final no comprobaba siete propiedades del contrato

**Gravedad: media. Estado: cerrada el 5 de octubre de 2026.**

Al cerrar DT-026, el juez del FASTQ final —el de cada mate y el de la pareja,
en `src/fastq_validation.rs`— exigía la forma de cada registro, pero no siete
propiedades del contrato. **Nada se publicaba mal por ellas**: el juez estricto
de la revisión adversarial del 05/10/2026 (`judge_full.py`) no encontró ningún
incumplimiento en 16,4 millones de parejas. Si el escritor las hubiera
incumplido, el juez no lo habría visto.

**Cómo se cerró.** El juez ya no deduce lo esperado del fichero: quien publica
se lo dice. `validate_final_fastq_mate` recibe una `FinalFastqIdentity` (evento,
biblioteca, UDI `i7+i5` y longitud de lectura) y `validate_final_paired_fastq`
recibe además unos `FinalFastqCounts` (parejas, bases y `N`). El escritor
(`FinalPairedFastqStream`) fija la identidad al abrir y los recuentos al cerrar.
Cada rechazo tiene su motivo propio:

1. **Identidad del fichero.** Cada registro se compara con el evento, la
   biblioteca y el UDI esperados. **Salto de índice:** una pareja que salta se
   escribe en el fichero de su biblioteca *final* y lleva el UDI de esa
   biblioteca, porque el demultiplexador asigna a la vez biblioteca y UDI
   finales (`src/demultiplexing.rs`, `accept`). Así que un solo UDI esperado vale
   para todo el fichero y un salto legítimo no se rechaza; lo comprueba el
   control positivo de `tests/f3_end_to_end.rs`, con dos bibliotecas, salto a
   tasa uno y contaminación.
2. **Pareja repetida** dentro de un fichero, comprobada solo en R1: la
   comparación de cabeceras de R1 con R2 ya ata cada registro de R2 al de R1.
   Se guardan los 128 bits de cada identificador en un `Vec<u128>` (16 bytes por
   pareja) y se ordenan al final. Un primer intento con un `HashSet` de los 26
   símbolos en cada mate costaba memoria: medido con una ejecución de 5 Mb de
   referencia a 30x (alrededor de medio millón de parejas), el pico pasó de
   510 264 kB (base) a 551 716 kB (+8,1 %); con el `Vec<u128>` en R1 es de
   510 480 kB (+0,04 %), con los FASTQ idénticos byte a byte.
   **Límite que conviene decir:** REQ-DAT-68 pide la unicidad por evento y esta
   comprobación es por fichero. La de todo el evento la impone DTC-039 al
   derivar los identificadores (`TechnicalIdentityRegistry` detecta la
   colisión), sin comprobación independiente posterior.
3. **Forma de los identificadores:** `LIB_` y 26 símbolos Base32 para la
   biblioteca (las de un FASTQ publicado son siempre automáticas,
   `derive_library`); para la pareja, `TP_` y los 128 bits canónicos, es decir,
   con los 2 bits sobrantes del último símbolo a cero, reutilizando
   `decode_derived_id` de `src/identity.rs` mediante
   `decode_technical_pair_id` (DTC-087). Un `TP_` con 26 símbolos Base32 pero
   no canónico se rechaza.
4. **Longitud de lectura** igual a la del evento, y la propia longitud esperada
   dentro de 1 a 150.
5. **`N` con Phred 0.** Solo en ese sentido: una `N` exige `!`; el contrato no
   prohíbe un Phred 0 sobre una base `ACGT`.
6. **R1 frente a R2 en longitud.** Queda cubierto por el punto 4, porque los dos
   mates deben medir lo mismo que el perfil. Una comparación directa entre R1 y
   R2 no podría dispararse sola, así que no se escribió; la prueba
   `mates_of_different_lengths_are_rejected` lo ejercita por el motivo de la
   longitud.
7. **Recuentos de parejas, bases y `N`.** Cada mate devuelve lo que recuenta al
   leerse, sin una pasada más, y el juez de la pareja lo compara con los
   contadores que el escritor acumuló, que son los que el manifiesto del evento
   publica. **Límite que conviene decir:** la comparación ocurre al publicar,
   contra el resumen que se va a escribir; no hay una relectura posterior de un
   manifiesto ya publicado contra los ficheros, porque `src/event_manifest.rs`
   solo recibe bytes de manifiesto, sin acceso a los ficheros. No cambia ningún
   formato publicado.

**Pruebas.** Nueve nuevas en `tests/final_fastq_record_validation.rs` (una por
hueco, más la longitud esperada fuera de 1 a 150 y el recuento por mate) que exigen el motivo exacto, con casos de un `TP_` de 27 símbolos, uno no canónico, un i7 igual con otro i5 y un repetido no contiguo;
el control positivo de `tests/f3_end_to_end.rs` ahora pasa la identidad y los
recuentos del manifiesto; y el de un FASTQ publicado de verdad toma la
biblioteca y los recuentos del manifiesto y los compara con los ficheros, en un
escenario que sí tiene `N`. Los casos `undecl_mixed_*`, `undecl_duplicate_pair_id`
y `hole_*` de la revisión, regenerados con identificadores canónicos, se
rechazan todos con su motivo (la base los aceptaba).

La batería pasa de 1690 a 1699 pruebas, sin fallos y con 11 ignoradas.
Exploración: 140 escenarios (extremos de longitud 1 y 150, varias bibliotecas,
salto de índice, contaminación, quimeras, muchas `N`, PCR y tartamudeo) con la
base `8ac6e02` y con el árbol final: ninguno aborta y los 818 FASTQ son
byte a byte idénticos; el juez externo estricto da OK en los 140.
La revisión adversarial repitió la comparación con 2200 escenarios, sin abortos
nuevos, y tras sus correcciones se repitió en 320 más (1712 FASTQ y los
manifiestos de evento, idénticos byte a byte frente a la base; juez externo OK
en los 320).

---

## DT-030. Un bloque BGZF con un tamaño imposible hace abortar al juez

**Gravedad: media.**

`noodles-bgzf` 0.47 no acota el tamaño descomprimido que declara el pie de cada
bloque. Un bloque que declara más de 64 KiB hace que el lector **entre en
pánico** —«range end index 65537 out of range for slice of length 65536», en
`noodles-bgzf-0.47.0/src/io/block/data.rs:51`— y el programa termina con
código 101 en vez de rechazar el fichero con un diagnóstico. Lo encontró la
revisión adversarial del 05/10/2026 (casos `panic_isize_65537` y
`rej_truncated_compressed_block_plus_eof`) y se reprodujo al cerrar DT-026
con el juez de cada mate ya corregido.

- **Afectados confirmados**: el juez de cada mate del FASTQ final y la
  relectura de cabeceras del juez de la pareja
  (en el momento de abrir este asunto, `src/fastq_validation.rs:436`; hoy está
  en las líneas 384 a 389), que ya podía abortar así antes de DT-026.
- **Afectados, reproducidos el 06/10/2026 con un VCF manipulado**: los dos
  lectores de VCF comprimidos del juez de los recursos publicados
  (`src/results.rs:659` y `:2689`) y el juez semántico del VCF
  (`src/vcf_semantic_validation.rs`), que usan el mismo lector.

El fichero lo produce el propio simulador, así que el riesgo real es un
fichero corrompido entre la escritura y la relectura, no una entrada del
usuario. Aun así, un juez no debería abortar.

**Estado: cerrada el 06/10/2026, sin tocar `Cargo.toml` ni `Cargo.lock`.**

### Las dos salidas, comparadas

**(a) Comprobar el pie de cada bloque antes de entregarlo al lector. Es la que
se aplicó.** El formato es sencillo: cabecera fija de 18 bytes con el tamaño del
bloque (`BSIZE`) en los dos últimos, y `ISIZE` en los cuatro últimos del bloque.
`src/bgzf_blocks.rs` recorre el fichero con el mismo encuadre que el lector,
lee 18 bytes de cabecera y 4 de pie de cada bloque, salta el resto con un
desplazamiento sin descomprimir nada, y rechaza un `ISIZE` mayor de 65 536.
Coste: dos lecturas pequeñas y un salto por bloque, y los bloques son de 64 KiB
**descomprimidos**: 784 MB descomprimidos son 12 026 bloques. Medido por la
revisión, la guarda tarda de 26 a 33 ms, frente a 0,56 a 0,84 s de leer el
fichero entero. Una cabecera que no es BGZF o un bloque cortado
**no** los rechaza la guarda: detiene el recorrido y deja que el lector los
diagnostique, como ya hacía. No cambia ninguna salida publicada (solo lee).

**(b) Subir el grafo.** `Cargo.lock` ya lleva **dos** `noodles-bgzf`: la 0.47
(la fija `Cargo.toml` para `noodles-csi` 0.56, `noodles-tabix` 0.62 y
`noodles-vcf` 0.88) y la 0.51 (de `noodles-fasta` 0.66). Que la 0.51 sea la
única exigiría subir a `noodles-csi` 0.60 o 0.61, `noodles-tabix` 0.66 o 0.67 y
`noodles-vcf` 0.92 o 0.93, que piden la 0.51 (las 0.62, 0.68 y 0.94 piden la
0.52). DTC-043 pide documentar propósito, mantenimiento, características y
scripts de construcción de cada dependencia nueva, y DTC-048 exige que toda
actualización sea «deliberada y aislada» y vaya con pruebas, validadores
independientes y escenarios de rendimiento. El riesgo está en la salida: cuatro
versiones menores de crates 0.x pueden cambiar la API y, sobre todo, los bytes
del VCF comprimido y de su índice, que deben salir idénticos. Haría falta
comparar el VCF `.gz`, el `.tbi` y el `.csi` byte a byte en un conjunto amplio
de escenarios.

**Elección.** (a): es pequeña, no cambia ninguna salida y no condiciona la
decisión de (b), que sigue siendo razonable por otros motivos (quitar una
versión duplicada) pero ya no es urgente.

### Qué cubre

La guarda se llama **antes** de cada lector BGZF de un juez de ficheros
publicados:

- cada mate del FASTQ final (`validate_final_fastq_mate`);
- el recuento de bases de tracto del VCF comprimido (`recount_vcf_tract_bases`)
  y la comprobación de que el VCF nombra la referencia
  (`validate_published_reference_names`), en `src/results.rs`;
- el juez semántico del VCF (`validate_published_vcf_semantics`).

**Dónde no está, y por qué.**

- La relectura de cabeceras de `validate_final_paired_fastq` no la lleva:
  `validate_final_fastq_mate` acaba de recorrer y leer enteros esos mismos dos
  ficheros, y un fichero que pasara una vez y no la otra solo existe si alguien
  lo cambia entre las dos aperturas.
- El juez del FASTQ de F3 (`validate_wgs_fastq`) la tiene sin llevar, pero no
  lo llama ninguna ruta de publicación (solo sus pruebas); queda anotado.
- Tampoco la llevan los lectores del **productor** que releen lo que él mismo
  acaba de escribir (`src/fastq_output.rs`, `src/vcf_publication.rs` y
  `src/vcf_contract.rs`).
- **Los índices `.tbi` y `.csi` también son BGZF** (`noodles-tabix` y
  `noodles-csi` entran por la 0.47), y uno manipulado **sí** hace entrar en
  pánico a `read_index()` (reproducido por la revisión). Pero **ningún juez de
  lo publicado lee esos índices**: solo los lee el productor
  (`src/vcf_publication.rs`), justo después de escribirlos, de modo que no
  hay entrada ajena que lo alcance. No es una sospecha, es un hecho acotado:
  si algún día un juez de lo publicado los lee, hará falta la misma guarda.

### Pruebas

Un FASTQ con un `ISIZE` de 65 537 y de `u32::MAX`, el caso de un bloque
truncado seguido del marcador de fin, el de un tercer bloque corrupto en un
fichero de varios bloques (`tests/final_fastq_record_validation.rs`), un VCF
comprimido con un bloque imposible en los dos lectores de `src/results.rs`
(`tests/repeat_tract_budget.rs`) y en el juez semántico
(`tests/vcf_publication.rs`), más cuatro pruebas de la propia guarda
(`src/bgzf_blocks.rs`, una de ellas para un `BSIZE` menor que cabecera más pie). Cada una tiene su control positivo: el fichero válido
se acepta.

---

## DT-031. Una referencia corta con muchas quimeras aborta como defecto interno

**Gravedad: baja.**

En 3 de los 804 escenarios válidos de la revisión adversarial del 05/10/2026
—referencias de 1200 bases o menos con una tasa de quimeras alta—, la
ejecución termina con «Una invariante interna considerada imposible se ha
roto» (`GPS-INTERNAL-001`, `Chimera(NoSecondSource)`, salida 70). La base
`2eaced4` hace lo mismo, así que no es efecto de DT-026. Se reprodujo al
cerrarla con el caso `x0028` de la revisión: referencia simulada de 600 bases,
lecturas de 2 bases, cobertura 0,5 y `chimera_rate = 1`.

DTC-084 prevé el fallo: «una candidata ausente para una pareja seleccionada
falla antes de publicar, sin reintento ni reducción silenciosa de la tasa».
Lo que no prevé es que se presente como un **defecto interno** del programa:
para quien lo ejecuta, es una configuración que no cabe en su referencia y
debería recibir un diagnóstico de configuración que lo diga.

### Cómo se cerró (06/10/2026)

**Reproducida** con una sola copia (`ploidy = 1`). El caso del registro con
dos copias publica, y el fallo con una sola **depende de la semilla**: con la
configuración de la prueba, 19 de 20 semillas publican y `validate` la acepta.
Que con dos copias la otra copia ofrezca siempre un segundo origen no está
demostrado: no se ha comprobado con una copia que lleve deleciones grandes.

**Corrección.** `wgs_paired_sampler_error` (`src/execution.rs`) traduce
`Chimera(NoSecondSource)` a `GPS-DOMAIN-003` (salida 1), el mismo código del
precedente DT-006, con el individuo como `entity`, `segundo_origen_de_quimera`
como `invariant` y, en `detail`, la dependencia de la semilla y la sugerencia
de una referencia más larga o una tasa de quimeras menor. Se descarta una
salida de configuración (`GPS-CONF-008`, salida 2) porque el fallo depende de
la semilla y DTC-044 reserva la salida 2 a la configuración inválida. Sin
código nuevo, sin reintento, sin bajar la tasa y sin tocar DTC-084. Además,
`select_second_source` (`src/chimera.rs`) devuelve ahora su propia variante,
`SelectedSourceOutOfRange`, cuando el índice sorteado se sale de las candidatas
contadas, que es una invariante interna y no la falta real de candidatas: la
traducción solo atrapa esta última.

**Pruebas.** `a_chimera_rate_that_does_not_fit_a_short_reference_is_a_domain_diagnostic`
(`tests/f4_chimeras_pairs.rs`), de extremo a extremo: código, categoría,
salida 1, individuo, invariante y sugerencia, con control de tasa 0 que
publica. **Usa la semilla de las pruebas porque con ella el sorteo falla**; no
demuestra que 600 pb con una copia y tasa 1 nunca quepan. Y
`a_selected_index_outside_the_counted_candidates_is_not_a_missing_source`
(`src/chimera.rs`), que ejercita la rama de la variante nueva.

---

## DT-032. El lector de la verdad acepta objetos escritos como lista JSON

**Gravedad: baja.** Nada se publica mal por esto: el escritor siempre emite
objetos. Es un hueco del lector, anterior al cambio de DT-027 a DT-031, que
encontró la revisión adversarial del 06/10/2026.

`read_truth_jsonl` (`src/truth_validation.rs`) deserializa con `serde` unas
estructuras que, por defecto, **aceptan también una secuencia de valores en el
orden de los campos**. `deny_unknown_fields` no lo impide. Así, el
`truth.jsonl` con `"counts":[6,8]` en lugar de `"counts":{"copies":6,"genotypes":8}`
se lee, y el esquema JSON de la verdad (`schemas/truth-v4.schema.json` y
`truth-v5.schema.json`, validados con `jsonschema` 4.10.3) lo rechaza.

**Reproducción medida el 06/10/2026**, sobre un paquete real con poblaciones
(dos individuos diploides, una fila de matriz): se sustituye, una a una, cada
estructura del documento por la lista de sus valores, 85 alteraciones. El juez
de los recursos publicados **acepta 40** y `jsonschema` rechaza las 85. Entre
ellas, `counts` de la cabecera, la fuente de cada huella (`source_digests[*].source`
y `group_digests[0].source`), la procedencia de cada copia y de cada
modificación, la fuente de cada segmento y de cada pieza, la procedencia
de cada SNP y cada alelo y cada portador de un genotipo. El revisor había anotado
tres casos con `counts`, la fuente de la huella y la procedencia de una copia;
el barrido los amplía.

**Qué haría falta.** Un deserializador de estructuras que rechace secuencias
(una pasada por un `serde_json::Value` o `deserialize_map` explícito), o una
comprobación del esquema en el lector.

### Cómo se cerró (06/10/2026)

**Solución.** Antes de deserializar cada línea (`parse_wire_text`, que usan el
lector completo y el de flujo), `ObjectShape` recorre el JSON una vez, sin
construir nada, y rechaza toda lista que no aparezca bajo uno de los trece
nombres que el esquema declara como lista (`TRUTH_LIST_FIELDS`: `alleles`,
`allele_dosages`, `blocks`, `carriers`, `current_pieces`, `group_digests`,
`history`, `modifications`, `piece_links`, `segments`, `snps`,
`source_digests`, `sources`), y toda lista cuyos elementos sean a su vez listas.
Cubre de golpe todos los tipos, los presentes y los futuros, sin tocar ninguno.
Falla cerrado, con una condición: un tipo nuevo queda protegido **mientras sus
nombres de campo no coincidan con uno de la lista blanca** (un objeto que se
llamase `sources` quedaría abierto), y la prueba de nombres lo exige tras la
revisión; un campo de lista nuevo que se olvide aquí rechaza una verdad válida
y la prueba lo delata. Se descartaron
un envoltorio por campo o por tipo (obliga a no olvidarlo en cada tipo nuevo y
no alcanza los structs dentro de enums etiquetados por dentro, que `serde`
vuelve a leer desde un almacén que sí acepta secuencias) y un
`serde_json::Value` previo (doble asignación por línea y, por sí solo, tampoco
rechaza secuencias). No cambia el esquema ni lo que escribe el productor.

**Medido el 06/10/2026.** El corpus de 85 de la revisión no se pudo
reconstruir (sus generadores no estaban); se sustituyó por un barrido
sistemático **mayor**: cada objeto de cada línea de cuatro verdades (la
fixture con matriz, una inserción, una deleción y una verdad real publicada
con cuatro individuos, matriz y un locus compuesto de tres bloques) se
sustituye por la lista de sus valores: **350 alteraciones, 33 formas de objeto
distintas**. Antes, el lector aceptaba **139** y rechazaba 211; ahora rechaza
las **350**. Frente a `jsonschema` 4.10.3 (esquema de la versión de cada
documento), antes había 139 casos aceptados por el lector y rechazados por el
esquema, ahora **ninguno**; los cuatro controles sin alterar los aceptan los
dos. El barrido no puede reproducir que «40 de 85» del informe fuesen las
mismas 40: la cifra comparable es la de este corpus.

**Coste medido por el revisor:** unos 1 ms por MB de verdad (145 a 167 ms en
21 MB, +15 % en el parseo; el juez completo no se midió). **Efecto en la puerta
de cobertura:** como el lector por flujo y el completo comparten ahora la misma
línea de rechazo, la puerta ya no distingue cuál de los dos la dispara.

**Pruebas.** `every_struct_written_as_a_list_is_rejected_by_the_truth_readers`
(`tests/truth_validation.rs`), el barrido, con la fixture nueva
`tests/fixtures/truth-v5/repeat-population.jsonl`; y
`list_field_names_are_exactly_the_schema_array_properties`
(`src/truth_validation.rs`), que fija la lista contra las propiedades de lista
del esquema v5 (y la v4 como subconjunto). Seis mutaciones, todas muertas
(sin la comprobación, `counts` permitido como lista, elementos de lista
aceptados, el lector de flujo sin la comprobación, la lista aceptada siempre y
`blocks` fuera de la lista).

**Corrección tras la revisión adversarial.** (1) La lista de valores del
barrido sale en orden **alfabético** de claves, no en el de los campos de
Rust, de modo que en muchos casos el propio tipo ya rechazaba; el comentario
decía lo contrario. Generarla en el orden de Rust exigiría conocer cada tipo
campo a campo, que es justo lo que esta solución evita, y no se hizo. En su
lugar el barrido exige que el error lleve el texto del rechazo nuevo, y así
solo cuenta un caso si actúa la comprobación nueva: añadir `provenance`,
`identity`, `copy`, `locus` o `source_event` a la lista blanca ahora lo mata el
barrido solo. `units` (un escalar, no un objeto) sigue sin matarlo el barrido,
que solo sustituye objetos; lo mata la prueba de nombres. (2) La prueba de
nombres separaba solo las listas; ahora reúne también las propiedades que no
lo son y exige que ninguna lleve un nombre de la lista blanca, y reconoce una
lista por `type` (cadena o lista de tipos), `$ref` y `oneOf`/`anyOf`/`allOf`.
Mutando el esquema en una copia, mueren un objeto llamado `sources`, una lista
`["array","null"]`, una lista por `oneOf` y el control `"type":"array"`;
reescribir el esquema sin cambiarlo no la rompe. No reconoce otras formas de
declarar una lista (`items` sin `type`, por ejemplo).

**Otros lectores con el mismo hueco.** DT-033, ya cerrada, recogió el de la
verificación del paquete F5. La auditoría de identidad
(`audit_f5_identity_closure`, `src/f5_closure.rs`) acepta la `identity` de una
copia escrita como lista y pierde su `reference_id`, pero el oráculo F5 la
rechaza (por huella y, recalculadas, con «copy.identity no es un objeto»), de
modo que queda cubierto de extremo a extremo. `src/results.rs` y
`src/event_manifest.rs` leen `serde_json::Value` a mano y no tienen el hueco.

**Otras dos observaciones de la misma revisión, que no son este defecto y no
se corrigen aquí.** Un `execution_id` cambiado en la cabecera de la verdad y
un `source_digests` de más pasan el juez, con y sin poblaciones: no hay norma
que fije que el `execution_id` de la verdad coincida con el del manifiesto, ni
que prohíba entradas de fuente que ningún registro cita.

---

## DT-033. La verificación del paquete F5 publicado acepta enlaces escritos como lista JSON

**Gravedad: baja.** El productor no escribe así; es un hueco del juez, hallado
por la revisión adversarial de DT-032 y reproducido sobre `ed50bbf`.

`verify_published_f5_package` (`src/f5_publication.rs`) lee el manifiesto con
las estructuras `gate::Link` y `gate::Root`, que usan `derive(Deserialize)`
con `deny_unknown_fields`. Por defecto `serde` acepta también una secuencia con
los valores en el orden de los campos. El juez **acepta cuatro manifiestos que
el esquema `manifest-general-v5` rechaza**: el enlace de la verdad escrito como
lista, un enlace de artefacto como lista, todos los enlaces como lista y el
manifiesto entero como lista. Ningún otro juez de Rust lee ese manifiesto, así
que el hueco existe de extremo a extremo. Es el mismo mecanismo de DT-032.

**Cómo se cerró (06/10/2026).** No hizo falta la pasada previa de forma de
DT-032. Aquí los objetos que se leen con forma tipada son pocos y están todos
en `mod gate`, y los campos que sí son lista (`artifacts`,
`initialization_summaries`) ya están tipados como `Vec`, de modo que no hay
lista blanca que mantener ni que fijar contra el esquema. Se añadió un
envoltorio `gate::Strict<T>` que lee el valor con `deserialize_any`, acepta
solo un objeto y delega en `T` con `MapAccessDeserializer`: el flujo y
`deny_unknown_fields` quedan intactos, no se materializa nada y la lista se
rechaza en `visit_seq` con `rejected(..)` («an object, not a list of its
values»). Se envuelven la raíz del manifiesto, `/f5` de éxito y de fallo, la
verdad, cada enlace y cada resumen de inicialización, `f5_failure`, la
configuración efectiva (raíz y `/f5`) y el diagnóstico `f5-failure.json`: los
doce tipos que `gate` deserializa con `derive`. Se usa `deserialize_any` y no
`deserialize_map` porque `serde_json` rechaza entonces la lista por su cuenta,
sin pasar por `rejected`, y la puerta de cobertura no vería la comprobación.
Sin cambios en el esquema ni en el productor.

Prueba: `the_gate_rejects_every_object_it_reads_when_written_as_a_json_list`
(`tests/f5_publication.rs`) sustituye por la lista de sus valores cada objeto
anterior (raíz, `/f5`, verdad, los siete enlaces, los resúmenes y la raíz y
`/f5` del efectivo, en el éxito; raíz, `/f5`, `f5_failure`, su enlace y el
diagnóstico, con el enlace rehecho, en el fallo) y exige el texto del rechazo
nuevo, no uno cualquiera; y contrasta cada caso con el esquema
`manifest-general-v5` (o `effective-config-v5`, `f5-failure-v1`) con el
`python3` del entorno virtual: todos rechazados por los dos y los controles
intactos aceptados por los dos. La lista sale en orden alfabético de claves,
así que contra la puerta sin arreglar muchos casos mueren ya por tipo; lo que
discrimina es el texto. La aceptación real se reprodujo aparte con la verdad y
un enlace en orden de campos: la puerta anterior devolvía `Ok` y el esquema
lo rechazaba. Mutaciones sobre copia (nueve, todas muertas): lista aceptada en
`visit_seq` y el envoltorio retirado, uno a uno, de la raíz, `/f5` de éxito,
`f5_failure`, la verdad, los enlaces, los resúmenes, el efectivo y el
diagnóstico.

Cobertura: `src/f5_publication.rs` pasa de 52 / 36 a 53 / 36 (la comprobación
nueva la dispara la prueba, el techo no cambia); el juez incorpora `rejected` a
su entrada del array de `scripts/check-guard-coverage.sh`; sin ello la puerta
no contaba la comprobación nueva. `tests/guard_constructors.rs` no cambia: ya
trata `rejected(` de forma fija.

**Fuera de alcance, sin cambios.** Los campos `Value` del manifiesto
(`limits`, `profile`, `sampling`, etc.) no se leen con forma tipada: la puerta
solo los compara con la copia efectiva (y, en el éxito, lee `reference_bases_processed`
de `resource_preflight`); no llama al oráculo, que en producción recibe su
declaración desde memoria y no desde el manifiesto. Además salta sin leer
`attempt_summary`, `quota_ties` y `f_realized`. Ver DT-035.
`audit_f5_identity_closure` sigue como en DT-032 y no se tocó.

Los envoltorios `Strict` de `SuccessDocument` y `FailedDocument` son defensa
redundante que ninguna prueba alcanza: la puerta lee antes el mismo fichero
como raíz envuelta y la lista se rechaza allí. Se dejan para que la lectura no
dependa del orden de las lecturas.

---

## DT-034. La atribución de la calibración del tartamudeo omite a Gymrek (2016)

**Gravedad: baja.** Decidido por el usuario: se corrige en un cambio aparte.

Lo publicado atribuye las cifras de la calibración del tartamudeo solo a «Raz
et al. (2019), Nucleic Acids Research 47(5):2436-2445, CC BY 4.0» (constante en
`src/stutter.rs`, hacia la línea 153). Las cifras de la tabla por longitud de
motivo son de Gymrek (2016), tabla 1, CC BY 4.0, cuya licencia exige
atribuirlas; Raz aporta la forma del modelo. `design/calibracion-tartamudeo-0.1.md`
(líneas ~51 y 223-224) lo dice bien. Gymrek no aparece en el código ni en los
esquemas.

**Qué haría falta.** Corregir la constante y lo que la cite. Cambia los bytes de
toda ejecución con tartamudeo, porque el texto de la fuente va en la salida; por
eso va aparte del multiplicador configurable de DT-009. No se ha corregido.

---

## DT-035. La puerta del paquete F5 acepta formas que el esquema rechaza en campos no tipados

**Gravedad: baja.** El productor no escribe así, pero la puerta también juzga
paquetes ya publicados. Hallado por la revisión adversarial de DT-033 y
reproducido con la puerta arreglada y con la antigua.

La puerta acepta cuatro manifiestos que `manifest-general-v5` rechaza, en campos
que no lee con forma tipada: `attempt_summary` escrito como lista, `quota_ties`
como texto, `limits` como lista en el manifiesto y en la configuración efectiva
a la vez, y lo mismo con `profile`. Los campos `Value` solo se comparan entre
las dos copias, de modo que dos copias igual de malformadas coinciden, y la
puerta salta sin leer `attempt_summary`, `quota_ties` y `f_realized`.

**Qué haría falta.** Validar la forma de esos campos contra el esquema o con
comprobaciones propias. No se ha corregido.

---

## DT-036. Un ensamblado de evento fallido no deja el evento enlazado con su estado

**Gravedad: media.** Defecto contra REQ-SAL-21, que exige un estado explícito
para cada evento y que el manifiesto general enlace los estados de sus unidades.
Hallado el 06/10/2026 al cerrar DT-002. **Sin arreglar**, a resolver en un
trabajo aparte.

`TechnicalEventAssemblyFailure` (en `src/execution.rs`) tiene un campo
`link: Option<ManifestEventLink>`, y la ejecución lo lee tras un fallo
(`event_links = failure.link.map(|link| vec![link])`) para decidir si el paquete
fallido enlaza el evento. **Todas** las construcciones del fallo escriben
`link: None`, así que ese camino no se recorre nunca: un fallo del ensamblado
deja el manifiesto con `events: null`, sin el evento ni su manifiesto fallido.
La prueba ignorada
`tests/f3_event_assembly.rs::controlled_failed_event_keeps_only_closed_fastq_and_its_failed_manifest`
codificaba justo lo contrario —un manifiesto de evento fallido que conserva solo
los FASTQ cerrados—, y la de lecturas sueltas de `src/execution.rs` borrada en
DT-002 afirmaba un evento con estado `failed`.

Lo que sí cubre una prueba hoy
(`a_failed_event_assembly_leaves_no_fastq_nor_derived_reference`) es lo demás:
ni FASTQ, ni referencia derivada, ni artefactos parciales, y la causa
registrada. **Esa prueba deliberadamente no exige `events: null`**, para no
fijar el defecto. Arreglarlo es una decisión de contrato: o el campo `link` se
rellena con un evento fallido y su manifiesto, o se quita y el evento fallido se
declara de otro modo.

---

## DT-037. `tests/test_population_snvs_editor.py` está fuera de la puerta de Python y falla una prueba de siete

**Gravedad: baja.** Mismo tipo de problema que DT-004. La puerta canónica
(`unittest discover -s tests_ui`) no recorre `tests/`, así que este módulo no se
ejecuta nunca desde ella. Lanzado a mano (`python -m unittest
tests.test_population_snvs_editor`) falla 1 de 7:
`test_draft_schema_rejects_generate_with_supplied_or_reuse`, con
`'2' was expected`. **Sin arreglar.**

---

## DT-038. Un fallo posterior al evento descarta el VCF ya publicado y validado

**Gravedad: media.** Contra REQ-SAL-25 y DTC-045: un fallo no debe invalidar
retroactivamente artefactos ya publicados y validados que sean independientes, y
el VCF, una vez escrito y validado, lo es. `publish_controlled_failure`
(`src/execution.rs`) retiene solo la copia de la configuración, las
referencias, los FASTQ, los manifiestos de evento, el mapa de la inserción y,
desde el cierre de DT-002, las reducciones completas; el VCF y su índice se
borran del paquete fallido. Ya ocurría antes de DT-002. Reproducido por el
revisor de ese cierre: con VCF y evento, la reducción se conserva y el VCF se
descarta. **Sin arreglar**: exigiría decidir también qué hacer con la verdad y
la matriz que el VCF cita, que tampoco se retienen.

---

## DT-039. Con evento, alterar la referencia suministrada antes del cierre impide publicar el paquete fallido

**Gravedad: media.** Con un evento ya publicado, si el FASTA de la referencia
suministrada se altera antes del cierre, la propia rama de fallo controlado no
consigue publicar: falla con `GPS-FMT-006` «artefacto de evento no coincide con
el catálogo físico general», **no queda ningún destino y se pierde la causa
original**. Reproducido por el revisor del cierre de DT-002 con el gancho previo
al manifiesto. Un fallo controlado debería publicar siempre un paquete fallido
coherente o, al menos, conservar y devolver la causa original. **Sin arreglar.**

---

## DT-040. Con tartamudeo clonal, la configuración efectiva contradice su esquema en el catálogo de azar

**Gravedad: media.** Encontrada el 08/10/2026 al pasar por primera vez por el
validador externo de esquemas la configuración efectiva de una ejecución
**amplificada con loci de repetición**. Ninguna prueba lo hacía: las que
validan la efectiva por fuera usan escenarios amplificados sin loci.

**Qué pasa.** Con tartamudeo clonal, la efectiva publica
`technical_event.rng_catalog = "f4-technical-v5"` y, dentro,
`technical_event.artifacts.pcr.rng_catalog = "f4-technical-v6"`. El esquema
`effective-config-v6` exige `f4-technical-v6` en el primero cuando existe el
bloque `pcr.stutter` (su `allOf` con el comentario «El tartamudeo clonal
estrena el v6»), y el validador externo rechaza el documento. El productor
(`technical_event_json` en `src/results.rs`), el juez de la efectiva
(`validate_effective_technical_event`) y el del manifiesto general (que solo
admite v3, v4 y v5 en el inventario técnico) coinciden entre sí en el v5; el
manifiesto de evento publica v6. El código y el esquema se contradicen desde el
commit del tartamudeo clonal, `a40ce85` (26/09/2026); se reprodujo con el
binario del commit `fc2b205`.

**Decisión del responsable (08/10/2026): manda el contrato, es decir, el
esquema.** Se corregirá en un cambio **aparte**, junto con DT-041: el productor
publicará `f4-technical-v6` en `technical_event.rng_catalog` de la efectiva
cuando haya tartamudeo clonal —DTC-091: una ejecución publica el catálogo más
alto cuya etapa abre; DTC-092: solo la amplificada **sin** loci sigue en v5—,
el juez de la efectiva (`validate_effective_technical_event`) y el del
manifiesto general aceptarán el v6, y se ampliará el esquema del manifiesto
general. Mueve los bytes de toda ejecución amplificada con microsatélites, de
modo que hará falta una huella nueva decidida por escrito.

**Sin arreglar todavía**: el multiplicador de DT-009 no podía tocarlo sin romper
la igualdad byte a byte que exigía. `tests/clonal_homopolymer_fold.rs` comprueba
que este es el único motivo por el que la efectiva no pasa su esquema, y fallará
cuando se corrija para que se quite el rodeo.

---

## DT-041. El manifiesto general de toda ejecución amplificada incumple su propio esquema

**Gravedad: media.** Hallada el 08/10/2026 por la revisión adversarial de
DT-009 y comprobada después contra el binario de esa entrega: el manifiesto
general de los escenarios `amplificado-loci` y `amplificado-sin-loci` de
`evidence/dt-009-multiplicador-congelacion/` lo rechaza el validador externo con
`'gps-hiseqx-pcrfree-pe150-v1' was expected`.

**Qué pasa.** Los esquemas del manifiesto general v4 y v6 fijan, en
`technical_events[].profile.id`, el perfil sin PCR
(`gps-hiseqx-pcrfree-pe150-v1`), y en `rng_catalog` solo admiten `f4-technical-v3`
y `f4-technical-v4` (el v4 del esquema, solo el v3). Una ejecución amplificada
publica el perfil `gps-hiseqx-pcr-pe150-v1` y el catálogo v5 o v6, de modo que
**toda** ejecución amplificada, con microsatélites o sin ellos, publica un
manifiesto general que no pasa su esquema. El juez de Rust lo acepta (admite v3,
v4 y v5 y los dos perfiles). Ninguna prueba valida por fuera el manifiesto
general de una ejecución amplificada.

**Decisión del responsable (08/10/2026): manda el contrato, ampliándolo** para
admitir el perfil con PCR y los catálogos que la amplificación abre. Se resuelve
**junto con DT-040**, en el mismo cambio aparte. **Sin arreglar.**

---

## DT-042. En el cuaderno, la celda que reconstruye el formulario falla con un evento técnico activo

**Gravedad: media.** Hallada el 08/10/2026 por la revisión adversarial de
DT-009, que la reprodujo ejecutando el cuaderno. Con un evento técnico activo,
la celda que reconstruye el formulario (`ui/gps_app.py`, hacia la línea 1046)
falla con `NameError: name 'TechnicalEventForm' is not defined`: la celda lo
declara como parámetro, pero **ninguna celda lo define**, porque la de
importación no lo importa de `config_editor` (comprobado leyendo el fichero).
Viene de F4.1, commit `e341b57` (29/08/2026). Consecuencia: **ningún control
del panel de PCR llega al TOML desde el cuaderno en ejecución**, y tampoco los
del resto del evento técnico. `marimo check --strict` no lo detecta.

**No explica los dos fallos conocidos de `tests_ui`**, como se sospechaba: se
comprobó el 08/10/2026 que los dos fallan por un diagnóstico de Rust
(`GPS-CONF-008`, «la matriz SNV falta después de la validación efectiva» y
«truth.jsonl no contiene copias ni genotipos cuando initial_populations está
vacío»), no por este `NameError`. Por la misma razón, la prueba de la interfaz
del multiplicador (`test_the_control_value_reaches_the_form`) es estructural y
lo dice en su docstring. **Sin arreglar.**

---

## DT-043. El juez no ata los recuentos realizados de la amplificación del manifiesto de evento

**Gravedad: media.** Hallada el 08/10/2026 por la revisión adversarial de
DT-009, que la reprodujo con las huellas rehechas en el manifiesto general: el
juez de los recursos publicados **aceptaba** un manifiesto de evento con
`artifacts.pcr.cycles` cambiado a 12 —incoherente con los `cycles` del bloque
del tartamudeo clonal, que siguen en 10— y con `lineage_applied_steps_by_locus`
aumentado en uno. Desde DT-009 el **bloque** del tartamudeo clonal tiene que
ser idéntico al de la configuración efectiva (REQ-VAL-61).

**Cerrada en parte el 10/10/2026** (decisión del responsable: arreglarlo ya, en
un cambio propio, junto con DT-052). Los **parámetros** de la amplificación ya
se atan: el bloque `technical_event.artifacts.pcr` de la configuración efectiva
está atado exactamente a lo que resuelve el TOML copiado, sus claves (catorce,
quince con el tartamudeo clonal) son un subconjunto exacto de las del bloque del
manifiesto de evento —los dos salen de la misma función del productor—, y el
juez completo exige ahora que cada una aparezca **idéntica** en el manifiesto de
evento (`validate_physical_technical_event`, al final). El `cycles` = 12 se
rechaza. El catálogo de azar del bloque (`artifacts.pcr.rng_catalog`) entra en
la comparación sin decidir nada sobre DT-040: los dos documentos publican el
mismo, y DT-040 trata del catálogo de la **raíz** del evento en la efectiva, que
no forma parte del bloque. Pruebas, en `tests/results_repeat_guards.rs`:
`the_amplification_of_the_event_manifest_is_the_effective_one` (cada una de las
once claves alcanzables, alterada una a una con las huellas rehechas, sobre los
tres paquetes amplificados; todas aceptadas por el código anterior) y
`real_packages_pass_the_tract_and_amplification_comparisons` (lado aceptado).
Cinco claves no las alcanza esta comparación porque otra guarda las caza antes:
`applicable`, `algorithm`, `rng_catalog` y el método de la inversión los fija
como literales el juez propio del manifiesto de evento, y el sub-bloque
`stutter` lo compara antes la igualdad de REQ-VAL-61.

**Lo que queda abierto.** Los **recuentos realizados** del bloque del manifiesto
de evento —`families`, `pairs`, `duplicate_pairs`,
`realized_duplicate_fraction`, `family_size_histogram`,
`families_truncated_at_end`, `gc_attempts_total`,
`gc_attempts_exhausted_families`, `unique_molecule_coverage_x`,
`resource_bound` y los tres recuentos de linaje— no están en la efectiva y
nada los ata a las lecturas publicadas ni a una amplificación rehecha: solo los
comprueban su esquema y las relaciones internas que `validate_pcr_manifest`
mira (familias frente a parejas, histograma frente a totales, fracción
realizada frente a recuentos). Un `lineage_applied_steps_by_locus` aumentado
en uno con las huellas rehechas **sigue aceptándose** (comprobado ejecutando el
10/10/2026). Atarlos exigiría rehacer la amplificación con la semilla o recontar
sobre los FASTQ. El camino ya está pedido en otro sitio: REQ-VAL-58 exige, en la
validación de cierre, contar las parejas duplicadas sobre los FASTQ finales y
rederivar desde la semilla las familias y sus tamaños; llevar ese oráculo al
juez de publicación cerraría esta deuda (según la revisión adversarial de
DT-052, REQ-VAL-58 no tiene fila en la matriz de trazabilidad; no lo he
comprobado aparte). La regla de lo que sí se ata quedó en REQ-VAL-63.
**Sin arreglar.**

---

## DT-044. Once requisitos de la amplificación no tienen fila en la matriz de trazabilidad

**Gravedad: baja.** Hallada el 08/10/2026 por la revisión adversarial de
DT-009 y comprobada con `grep`: REQ-CONF-134 a 138, REQ-DAT-81 a 85 y REQ-MAN-94
no tienen fila en `design/trazabilidad-requisitos-0.1.md`, aunque la
especificación ya los presenta como implementados. **Sin arreglar.**

## DT-045. El juez de la matriz de SNV repite comprobaciones que ya hizo el analizador de metadatos

**Gravedad: baja.** Hallada en la tanda 3 de DT-001 (08/10/2026) y confirmada por
la revisión adversarial. Trece comprobaciones de
`validate_reusable_population_snvs_bytes` y la de la lista de
`restore_published_allele_frequencies` son **código muerto**: repiten sobre el
mismo `Value` lo que `parse_metadata_json` ya exigió, de modo que ningún dato
publicado puede alcanzarlas (la prueba
`the_shape_guards_of_the_reuse_judge_are_shadowed_by_the_metadata_parser` lo
ejecuta). Además, el conjunto `seen` del orden de grupos de
`validate_normalization_audit_document` es redundante: la secuencia aceptada es
estrictamente creciente.

**Simplificación propuesta** (decidida como encargo aparte, **sin hacer**): que
`parse_metadata_json` entregue tipados los campos que comprueba y que la
restauración reciba la lista ya extraída. El total del juez bajaría a unas 129
comprobaciones con 11 sin disparar. Las otras inalcanzables de la tanda 3 (cotas
de redondeo, conversiones, invariantes entre recorridos) no se eliminan así.

**Otra simplificación propuesta, de la tanda 8 de DT-001 (12/10/2026, sin hacer):**
el juez físico de la matriz de SNV de `src/results.rs`
(`validate_physical_population_snvs`) lee la tolerancia de normalización de la
efectiva publicada, y esa lectura tiene una guarda inalcanzable (la efectiva ya se
ató al TOML por bits). Podría leerla de la configuración ya reconstruida
(`physical_config.normalization().tolerance()`), como hace la rama de
reutilización, y esa inalcanzable desaparecería.

## DT-046. Referencias cruzadas internas del manifiesto de evento no atadas por su propio juez

**Gravedad: baja.** Halladas en la tanda 4 de DT-001 (08/10/2026) al cubrir
`src/event_manifest.rs`, ampliadas en la ronda 2 (09/10/2026). **Sin arreglar**
(decisión del agente raíz: registrar).

1. `validate_f4_id_catalog` compara `values.iter().filter_map(Value::as_str)` con
   la lista esperada: un elemento que no es texto se descarta, así que
   `"v1": ["execution", ..., 7]` **se acepta** (reproducido). Valoración del
   revisor: el juez lo acepta y `jsonschema` 4.10.3 lo **rechaza**.
2. `validate_selected_individuals` exige individuos y grupos de lectura únicos,
   pero **no** bibliotecas: dos individuos con el mismo `library_id` se aceptan
   (reproducido). Valoración del revisor: **los dos** (juez y `jsonschema`) lo
   aceptan.
3. Los diagnósticos de `validate_writer_facts` llevan el nombre provisional del
   esquema v2 aunque el documento sea v3, porque `event_manifest_document` los
   devuelve sin el retiquetado que sí aplica `validate_event_manifest_bytes`
   (mismo defecto que DT-017 en el otro juez; etiqueta reproducida por el
   revisor).
4. El juez aislado también acepta, y `jsonschema` 4.10.3 tampoco rechaza: una
   receptora de contaminación que no está entre los seleccionados, un donante
   ajeno, una biblioteca receptora que es la de otro individuo y dos bibliotecas
   con el mismo par de UDI. **No verificado por ejecución:** por lectura, el
   juez del paquete completo (`src/results.rs`, hacia las líneas 3281 y 3399)
   los cazaría al comparar `selected_individuals` con el plan rederivado.

## DT-047. El contrato admitía un cambio de base de segmento cubierto por una deleción del mismo segmento y la réplica lo rechazaba

**Gravedad: baja. Cerrada (09/10/2026, ronda 2 de la tanda 5 de DT-001).** El
contrato (`validate_structure`) admitía, en la lista de modificaciones de un
segmento, un cambio de base y una deleción que lo cubre, en los dos órdenes
(`deletion_covers_snp`), y la réplica del juez lo rechazaba como «el historial
no se reproduce». **El error estaba en el contrato**: `Segment::new`
(`src/segments.rs`) rechaza cualquier solape y `schemas/truth-v1.md` exige
modificaciones «ordenadas y no solapadas». Se eliminó la excepción; ahora el
contrato rechaza el documento con `ModificationsNotOrdered` (el cambio y después
la deleción) u `OverlappingModifications` (la deleción y después el cambio, o una
inserción dentro o al final de la deleción). Las pruebas
`modifications_of_a_segment_are_ordered_and_never_overlap` y
`an_inherited_base_change_covered_by_a_segment_deletion_is_rejected` exigen el
motivo exacto, y los mutantes que anulan el orden o el solape mueren. La
guarda de la réplica (`push_initial_source` con un inicio por detrás del cursor)
solo se alcanza ya con un documento en memoria que no pasó el contrato
(`the_replay_rejects_a_covered_change_in_a_document_that_skipped_the_contract`).
Cambio de producción: una condición menos en `validate_copy`; el número de
comprobaciones no cambia (69).

## DT-048. El juez no exige que el estado heredado de una copia recorra la referencia una vez y en orden creciente

**Gravedad: media.** Hallada por la revisión adversarial de la tanda 5 de DT-001
(09/10/2026) y confirmada por el responsable: **el estado heredado debe cubrir la
referencia exactamente una vez y en orden creciente**. **Cerrada el 09/10/2026**,
en un cambio aparte justo después de la tanda 5 (ver «Cómo se cerró»).

`inherited_loss` (`src/truth_validation.rs`) solo busca bases que **faltan**: no
detecta duplicaciones ni reordenaciones. Reproducido sobre la verdad publicada por
el programa:

1. La copia 1 de I3 con los segmentos `[300, 600)` y después `[0, 300)` (el
   cromosoma rotado) se acepta con los dos lectores.
2. La copia 1 de I1 en `chr1` con los segmentos `[0, 424)`, `[420, 424)` y
   `[424, 600)` más la deleción del historial `del-1` (420..424) se acepta: la
   copia es la referencia más su cambio en 49 y, aun así, el caso de deleción la
   declara portadora.

*Referencia normativa:* `design/especificaciones-0.1.md` (línea ~367, «intervalos
colineales») y `schemas/truth-v4.md` (exige «conservar la referencia entera»).

*Corrección propuesta* (se hizo con una adaptación, ver abajo): al formar la réplica inicial, exigir segmentos
**ordenados, sin solaparse y sin huecos**, y convertir en rechazos las dos pruebas
que se quitaron de la tanda 5 (`touching_source_fragments_in_reverse_order_are_accepted`
y `a_repeat_tract_touching_the_previous_fragment_is_accepted`, que daban por
válida la rotación).

*Mutante relacionado:* `cursor = cursor.max(end)` por `cursor = end` en
`inherited_loss` **sobrevive y no es equivalente**: solo lo distinguen documentos
con segmentos solapados.

### Cómo se cerró (09/10/2026)

**La regla.** Los segmentos de cada copia (`segments[*].source`), en el orden en
que se publican, recorren la secuencia de la referencia de la copia
**exactamente una vez y en orden creciente**: cada segmento es de esa secuencia y
no está vacío ni invertido (añadido en la ronda 2), el primero empieza en 0, cada uno
empieza donde acaba el anterior y el último acaba en la longitud de la secuencia
según la referencia. Las modificaciones de segmento no cuentan para el recorrido:
un cambio de base y una inserción heredados no consumen referencia, y una deleción
de segmento la sigue rechazando el juez del conjunto de casos. Los fragmentos
vigentes no necesitan regla propia: salen de replicar el historial sobre los
segmentos (inserciones entre fragmentos, deleciones que los recortan, tractos que
sustituyen un tramo de un fragmento), lo que nunca reordena ni duplica, y el juez
exige que coincidan con esa réplica.

**Dónde.** `require_collinear_segments` (`src/truth_validation.rs`), con el rechazo
nuevo `SegmentsNotCollinear { identity, reached, next_start }` —hasta dónde había
llegado el recorrido y dónde empieza el segmento siguiente (tras el último, la
longitud de la secuencia)—. Corre en `reconstruct_copy`, que es común a los dos
lectores (completo y por flujo), al documento en memoria, a las fronteras del VCF y
a la publicación. Sustituye a la búsqueda de huecos de `inherited_loss`, que ordenaba
los segmentos antes de mirarlos y por eso no veía rotaciones ni desórdenes; queda en
`inherited_loss` solo la deleción de segmento. La regla es del juez y no del
contrato porque necesita la longitud de la secuencia **de la referencia**.

**Adaptación respecto a la propuesta.** Se exige **después** de replicar y
materializar la copia, no al formar la réplica inicial. Exigida al principio, la
batería mostró once pruebas de guardas que dejaban de alcanzar su guarda (dos
fragmentos que se solapan, longitudes que desbordan, un fragmento más allá del
final, un cambio en la última posición representable, dos cambios en una posición,
inserciones y deleciones sin base de referencia delante): todas usan segmentos no
colineales para llegar a ellas, y habrían pasado a inalcanzables. Exigida al final,
esas guardas siguen parando lo suyo con su motivo, y la nueva para lo que ninguna
otra ve: el hueco, la rotación, el desorden y el tramo duplicado o solapado que una
deleción del historial vuelve a quitar. Un solape que llega a los fragmentos
vigentes sigue rechazándose antes como `DuplicateCurrentSourcePosition`.

**Ninguna verdad del programa la incumple.** Comprobado antes de tocar el juez, con
un comprobador independiente en Python sobre una copia del árbol instrumentada para
volcar cada verdad que el programa proyecta: **409 verdades (2857 copias)
publicadas por la batería entera** (1972 pruebas en verde), con inserciones y
deleciones (102 y 113 verdades), fragmentos insertados (96), tractos de repetición
(193; 80 con un *indel* en la misma copia), referencias simuladas con
microsatélites plantados y suministradas, varias secuencias (31), varias
poblaciones (36), ploidías 1, 2, 3, 4 y 10 y una secuencia de más de 100 kb; **35
verdades (1227 copias)** publicadas con el binario de *release* desde cada
configuración de los *fixtures* que publica (`vcf-fusion`, `snp-order`,
`vcf-null-effect`, `reference_reduction`); y las tres verdades ya publicadas en
`results/` que tienen copias (1300). Todas las copias tienen un único segmento
`[0, longitud)`: el productor solo crea cadenas intactas (`SegmentChain::intact`).
Cero incumplimientos.

**Pruebas** (`tests/truth_validation_guards.rs`, sección «El estado heredado recorre
la referencia una vez y en orden»), todas sobre la verdad publicada de verdad, con
el motivo exacto y los dos lectores cuando los dos pueden verlo:
`touching_source_fragments_in_reverse_order_are_rejected` y
`a_repeat_tract_touching_the_previous_fragment_is_rejected` (las dos de la tanda 5,
ahora rechazos), `a_duplicated_stretch_that_the_history_deletes_again_is_rejected`
(el caso 2), `a_partial_overlap_hidden_by_a_deletion_is_rejected` (segmentos
`[0, 424)` y `[422, 600)` con `del-1`: los fragmentos vigentes ya no se solapan),
`a_gap_between_segments_is_rejected`,
`segments_out_of_order_on_another_sequence_are_rejected` (`chr2`, tres segmentos
fuera de orden), `a_copy_without_segments_does_not_walk_the_reference` (en
memoria; también la secuencia ausente) y, del lado aceptado,
`collinear_segments_cut_anywhere_are_accepted` (las seis copias de `chr1` partidas
en tres segmentos, con inserciones, deleciones y, en la versión 5, el tracto). En
`tests/truth_validation.rs`, `both_readers_reject_a_segment_trimmed_at_either_end`
pasa del rechazo del conjunto de casos al nuevo. Contra el comportamiento anterior
(el juez de antes reconstruido en una copia), las ocho pruebas de rechazo fallan: la
rotación, el tracto rotado, el tramo duplicado, el solape escondido y el desorden en
`chr2` se **aceptaban**; el hueco y los recortes los paraba el conjunto de casos con
otro motivo.

**Mutaciones** (copia fuera del repositorio, un mutante por ejecución, huella SHA
antes, mutada y restaurada): quitar la llamada, `!=` por `<` o por `>`, quitar el
centinela de la longitud, empezar el recorrido en 1, tomar longitud 0 si falta la
secuencia, intercambiar los dos campos del error y que `inherited_loss` no vea la
deleción de segmento: **todos mueren**. `reached = end` por
`reached = reached.max(end)` sobrevivía y se anotó como equivalente porque un
segmento invertido lo rechazaría la réplica; **la revisión mostró que no siempre**
(ver la ronda 2, donde se endureció la regla). El mutante antiguo de `cursor`
desaparece con el código que lo contenía.

**Ronda 2 (09/10/2026, tras la revisión adversarial).** El revisor lo aprobó con un
fuzzer propio (cadenas colineales y perturbadas en `chr1` y `chr2`, versiones 4 y
5, los dos lectores): cero no colineales aceptados, frente a unos 190 (v4) y 40
(v5) con el juez anterior; 3000 documentos basura sin pánicos; ningún colineal
válido rechazado. Sus hallazgos y lo hecho:

- **La regla no dependía solo de sí misma.** Un segmento invertido en memoria
  (`[200, 100)` entre `[0, 200)` y `[100, 400)`, con deleciones de segmento que
  reponen el cursor) superaba la réplica; lo paraba el conjunto de casos, y las
  fronteras del VCF, que no lo juzgan, lo aceptaban. Ahora cada segmento tiene que
  **avanzar** el recorrido: se rechaza también un segmento invertido o **vacío**
  (`end <= start`). El vacío se rechaza porque no recorre nada y el contrato ya lo
  prohíbe por bytes (`InvalidSegmentInterval`); admitirlo solo en memoria no
  aportaría nada. Con eso, `reached = reached.max(end)` sí es equivalente por la
  propia condición (`end > start == reached`), y lo dice el comentario.
- **DT-050 se cierra de paso**: la misma condición exige que cada segmento sea de
  la secuencia de la copia.
- **Huecos de prueba**: un segmento que vuelve a empezar en 0 (`[0, 50)` y
  `[0, 400)` en `chr2`, con una deleción del historial que quita el duplicado) y
  una deleción de segmento detrás de otros cambios del mismo segmento (copia 1 de
  I3, 550..560), que distinguen una regla que tolerara reiniciar en 0 y un
  `inherited_loss` que solo mirara la primera modificación. Al mutar se vio un
  tercero: nada distinguía un recorrido que acaba **más allá** de la longitud
  (`!=` por `<` en la comprobación final sobrevivía); lo cubre
  `a_walk_past_the_end_of_the_sequence_is_rejected` (en memoria, el último
  segmento acaba en 610 y una deleción de segmento quita lo que sobra).
- **Las guardas «posición de fuente duplicada» (dos) y «un cambio en dos
  fragmentos»** son inalcanzables con documentos colineales y redundantes para la
  seguridad: quitarlas no dejaría pasar nada, solo cambiaría el motivo. Se
  mantienen por el diagnóstico, con una nota junto a cada una.
- **Textos normativos.** La regla se escribió en `schemas/truth-v4.md` como regla
  de la réplica, en su propio párrafo, y en DTC-061: rige en las versiones `"4"` y
  `"5"`, y una versión futura con variantes estructurales (duplicaciones,
  inversiones, translocaciones; REQ-IND-46 a 49 y REQ-REC-20) tendrá que
  sustituirla. **No se sube la versión del contrato** (decisión del responsable):
  el productor nunca escribió documentos no colineales, mismo criterio que DT-047.

Pruebas nuevas de la ronda 2: `an_inverted_segment_that_the_replay_lets_through_is_rejected`,
`an_empty_segment_is_rejected`, `a_segment_of_another_sequence_is_rejected` (las
tres en memoria), `a_segment_restarting_at_zero_is_rejected`,
`a_segment_deletion_after_other_inherited_changes_has_no_case` (las dos por bytes,
con los dos lectores) y `a_walk_past_the_end_of_the_sequence_is_rejected`. Mutantes
de la ronda 2, uno por ejecución en una copia: quitar la mitad «vacío o
invertido», `<=` por `<`, quitar la mitad de la secuencia, tolerar reiniciar en 0,
`!=` por `<` o por `>` en el inicio y en la comprobación final, quitar la
comprobación final, empezar en 1, longitud 0 si falta la secuencia, intercambiar
los campos del error, `inherited_loss` solo con la primera modificación y quitar la
llamada: **todos mueren**; `reached.max(end)` sobrevive y es equivalente.

**Cobertura.** `src/truth_validation.rs` pasa de 118 a 121 comprobaciones: el
rechazo nuevo en dos sitios (un segmento que no avanza donde debe y el recorrido
que no acaba en la longitud) y el de la secuencia ausente que necesita para
medirlo. Las tres se disparan y las sin disparar siguen en 32.

## DT-049. El lector por flujo de la verdad no aplica a los casos genotípicos el paso léxico del contrato

**Gravedad: media-baja.** Hallada por la revisión adversarial de la tanda 5 de
DT-001 (09/10/2026). **Sin arreglar** (decisión del responsable: registrar, junto a
DT-015). Solo afecta a la vía por flujo, que hoy solo usan las pruebas.

`validate_truth_jsonl_streaming` comprueba cada caso genotípico con
`validate_genotype` pero **no** con el paso léxico de `validate_structure`
(`validate_lexical_document`). Reproducido: acepta un `case_id` con un espacio o de
más de 64 caracteres, y un alelo de deleción cuyo evento fuente empieza en 421 en
vez de 420, que el lector completo rechaza. Las tablas léxicas de
`tests/truth_schema_guards.rs` se prueban por eso sobre el documento en memoria, no
con los dos lectores.

## DT-050. En memoria, el juez de la verdad aceptaba una copia con un segmento de otra secuencia en medio

**Gravedad: baja. Cerrada (09/10/2026), en el mismo cambio que DT-048.** Hallada
por la revisión adversarial de DT-048. Un documento **en memoria**, que no pasa
por el contrato, con una copia de `chr2` cuyos segmentos eran `chr2[0, 100)`,
`chr1[100, 300)` y `chr2[300, 400)` lo aceptaba `validate_and_reconstruct_truth`
(pública), y la copia reconstruida no era `chr2`: ni el recorrido de los segmentos
ni la detección de fragmentos duplicados miraban la secuencia de cada segmento.
Por bytes lo rechaza el contrato (`SegmentSourceMismatch`). Es de la misma familia
que las demás vías en memoria que saltan el contrato.

*Cómo se cerró.* La regla de DT-048 (`require_collinear_segments`) exige además que
cada segmento sea de la secuencia de la copia; si no, `SegmentsNotCollinear` en ese
segmento. Lo prueba `a_segment_of_another_sequence_is_rejected`, que falla al
quitar esa condición.

## DT-051. El juez de los recursos publicados no ata los loci de repetición publicados al TOML copiado

**Gravedad: media. Cerrada (10/10/2026)**; ver «Cómo se cerró» al final. Hallada el
10/10/2026 en la tanda 6 de DT-001, **reproducida ejecutando** sobre un paquete
publicado de verdad (un locus declarado con un alelo interrumpido y uno generado,
cinco individuos). El encargo de la tanda era registrarla; ninguna prueba de la
tanda daba por válido ese comportamiento.

`validate_published_resources` **acepta** estas alteraciones:

- en el **manifiesto general**, sección `repeat_loci`: la pureza, las unidades, la
  secuencia, la frecuencia declarada y las copias sorteadas de un alelo declarado;
  el tope y los mínimos de admisión de un locus y los de omisión; la secuencia de
  su identificador, su inicio, sus motivos; la procedencia de la referencia; la
  frecuencia de un alelo generado y la media de pasos de su bloque del generador.
  Solo lo comprueba la coherencia interna de la sección (forma, aritmética,
  identificador único) y el identificador de la referencia. (El texto original
  decía que también la contabilidad de tracto; **no era así**: su juez,
  `validate_published_tract_output_accounting`, no lo llama nadie al publicar.
  Corregido el 10/10/2026; ver DT-052);
- en la **configuración efectiva**, `repeat_loci`: para un locus declarado, cualquiera
  de sus campos (secuencia del alelo, frecuencia, pureza, inicio, tope, mínimos); para
  uno generado, la pureza y la frecuencia efectiva (su juez, `REQ-VAL-60`, solo
  rehace la secuencia, la frecuencia declarada, la contabilidad, los linajes, la
  media, el sesgo y el intervalo).

Reproducido además por la revisión adversarial:

- **(a) Trasplante del abanico.** El mismo escenario publicado con
  `mean_events = 1.0`, llevando su locus generado (parámetros y abanico) a la
  efectiva y al manifiesto del paquete original (3.0), lo aceptan el juez completo
  y `validate_published_tract_output_accounting`: los documentos afirman un abanico
  (47, 48, 50, 51, 52) que no es el de la verdad. El juez del abanico (REQ-VAL-60)
  lo rehace desde los parámetros **publicados**, que nada ata al TOML; el comentario
  de `results.rs` junto a `validate_published_slippage_fans` no se cumple si cambian
  también los parámetros.
- **(b)** `realized_copies` y `realized_frequency` del manifiesto no se comparan con
  la verdad: intercambiarlos entre alelos conservando la suma se acepta.
- **(c)** Manifiesto y efectiva no se comparan entre sí: quitar un alelo generado del
  manifiesto se acepta.

Es decir, los documentos pueden mentir sobre los microsatélites que describen sin
que el juez lo note; las lecturas, la verdad y el VCF se juzgan contra el TOML y no
se ven afectados. Se parece a DT-043 (la sección de amplificación del manifiesto de
evento).

**Corrección propuesta** (decisión del responsable: se hará en un cambio **aparte**
justo después de la tanda 6; aquí solo se registra). Tras
`validate_published_slippage_fans` y `with_generated_repeat_alleles`:

1. exigir `effective["repeat_loci"] == repeat_loci_json(&physical_config)`;
2. exigir que cada locus del manifiesto (sin `drawn_copies` ni `realized_*`) coincida
   con el de la efectiva, y que `reference_id`, `reference_source` y
   `admission_defaults` coincidan;
3. como segundo paso, contar las copias realizadas desde la verdad y compararlas con
   `drawn_copies` y `realized_*`.

La comparación va **después** del juez del abanico para conservar su independencia
(REQ-VAL-60): el abanico se rehace de lo publicado, y solo después se ata lo
publicado al TOML.

### Cómo se cerró (10/10/2026)

Decisión del responsable del 10/10/2026: arreglarlo ya, en un cambio propio, con
el precedente de REQ-VAL-61 (el juez reconstruye el bloque del tartamudeo clonal
desde el TOML y exige que los documentos lo repitan idéntico). Los tres pasos
propuestos se hicieron, en una función nueva de `src/results.rs`,
`validate_published_repeat_loci`, que el juez completo llama **después** del juez
del abanico, del de la referencia final y del de la referencia simulada, de modo
que cada uno conserva su motivo y el abanico se sigue rehaciendo desde lo
publicado (REQ-VAL-60). Con la configuración reconstruida desde el TOML copiado y
los alelos generados desde la referencia publicada, exige:

1. **Configuración efectiva.** `repeat_loci` idéntico al que publica esa
   configuración (`repeat_loci_json`), y `repeat_admission` con los mínimos del
   TOML. Cierra las alteraciones de los loci declarados y generados y el punto
   (a): el trasplante del abanico con sus parámetros ya no cuadra con el TOML.
2. **Manifiesto general.** La sección de loci idéntica a la que publica esa
   configuración (`repeat_loci_manifest_json`), salvo lo que solo sale del
   sorteo —copias sorteadas, copias y frecuencia realizadas— y la contabilidad de
   bases de tracto, que se ponen a cero en los dos lados antes de comparar. Ata
   `reference_source`, `admission_defaults`, cada campo de cada locus y alelo, y
   el número de loci y alelos (punto (c)).
3. **Copias realizadas (punto (b)).** Se cuentan **desde la verdad ya juzgada**
   (`truth_repeat_tally`, que devuelve ahora `validate_physical_truth`): por
   locus, cuántas copias llevan cada tracto según los casos genotípicos de
   repetición; una copia sin caso lleva el tracto de la referencia, porque el
   productor no aplica nada cuando el alelo sorteado es ese tracto. La verdad no
   publica qué índice de alelo le tocó a cada copia, así que dos alelos con las
   mismas bases y los mismos bloques son el mismo para ella: el juez compara la
   suma de las copias realizadas de cada grupo de alelos indistinguibles. Las
   copias sorteadas y la frecuencia realizada quedan atadas con eso, porque el
   lector del manifiesto ya exige que la frecuencia sea copias entre sorteadas,
   bit a bit, y que las copias realizadas sumen las sorteadas.

El comentario junto a `validate_published_slippage_fans` se corrigió: el juez del
abanico caza otros alelos **con los mismos parámetros**; lo que ata los parámetros
al TOML es la comparación nueva.

**Pruebas** (`tests/results_repeat_guards.rs`, sección DT-051), todas con motivo
exacto y control positivo, y todas **fallan contra el código anterior**
(ejecutado sobre una copia del commit anterior; la sonda de la revisión de la
tanda 6 dio «aceptado» en cada una de las alteraciones):

- `the_effective_repeat_loci_are_the_ones_the_toml_produces`: las dieciocho
  alteraciones de la efectiva (diecisiete de loci y la de `repeat_admission`);
- `the_repeat_section_of_the_manifest_is_the_one_the_toml_produces`: doce
  alteraciones del manifiesto, incluida la de quitar un alelo generado sin copias;
- `the_realized_copies_are_the_ones_the_truth_carries`: el intercambio de copias
  entre los alelos declarados, entre dos generados, y la copia del alelo igual a
  la referencia dada a otro (discrimina que una copia sin caso cuente como ese
  alelo);
- `a_fan_transplanted_with_its_parameters_is_rejected`: el trasplante en la
  efectiva, en el manifiesto y en los dos;
- `alleles_the_truth_cannot_tell_apart_are_compared_by_their_sum`: lado
  aceptado de dos alelos gemelos (repartir sus copias de otra forma se acepta) y
  lado rechazado (dar una copia de un gemelo al alelo distinto).

El lado aceptado lo da además toda la batería: cada paquete que el programa
publica pasa por el juez completo al publicar (loci declarados, generados,
compuestos, tractos insertados, ploidía 4, amplificado con y sin multiplicador
de los homopolímeros, referencia simulada), y ninguno cambió de resultado.

**Verificación por mutación** (copia fuera del repositorio, con su propio
directorio de compilación, un mutante por ejecución, huella antes, mutada y
restaurada, y `touch` tras restaurar): diecinueve mutantes. Murieron dieciocho:
anular cada una de las cuatro comprobaciones nuevas, invertir la condición de
versión, no tratar como referencia el alelo igual al tracto de la referencia,
ignorar los bloques del alelo o los de la verdad, quitar cada una de las cuatro
neutralizaciones del sorteo, no quitar los grupos sin copias, sumar mal un grupo
de gemelos (solo lo caza la prueba de los gemelos), no reconocer ninguna copia
en un caso, leer un alelo de repetición como referencia, no reconocer el caso por
su intervalo, y no filtrar las copias por secuencia (lo caza una publicación con
un locus en otro cromosoma, `repeat_slippage.rs`). Sobrevivió a esa pasada
quitar el inicio de la identificación del caso por intervalo, y se declaró
equivalente («como dos loci no se solapan, el final basta»). **No lo era**: la
segunda revisión adversarial lo mató con cinco pruebas de
`tests/vcf_overlap_fusion.rs`, que ejecutan la batería entera y no solo el
fichero de la sección. Una inserción directa en el borde derecho del tracto
publica en la verdad un caso de intervalo vacío `[final, final)` que comparte el
final con el locus y va después del caso de repetición; sin el inicio, se tomaría
por el caso del locus. Hace falta el intervalo entero. Quitar el final sobrevive
solo porque la verdad ordena la inserción del borde izquierdo `[inicio, inicio)`
**antes** del caso del locus y este la sobrescribe: es frágil si cambia ese orden,
y está anotado junto al código.

La regla queda fijada en REQ-VAL-62 (`design/especificaciones-0.1.md`), que
remite a DT-052 para la parte de la contabilidad de bases de tracto.

La puerta de cobertura de los jueces: `src/results.rs` pasa de 779 a 783
comprobaciones y sigue con 498 sin disparar; las cuatro nuevas se disparan.

La segunda revisión halló además un hueco de prueba: comparar en la efectiva solo
la pureza de los mínimos de admisión, y no la tirada mínima, sobrevivía a toda la
batería. Se añadió el caso «mínimos de admisión: tirada» a
`the_effective_repeat_loci_are_the_ones_the_toml_produces`, y ese mutante muere.

**Limitación, por diseño.** Las copias realizadas se atan a la verdad, y la verdad
no distingue dos alelos con las mismas bases y los mismos bloques. Con dos `(GT)24`
declarados a 0,85 y 0,05, pasar al de 0,05 las catorce copias del de 0,85 **se
acepta** (lo reprodujo la segunda revisión). Lo cerraría atar las copias al
sorteo rehecho con la semilla: DT-053.

**No se cubre aquí** (fuera del alcance de DT-051): la contabilidad de bases de
tracto (DT-052) y las copias realizadas frente al sorteo rehecho (DT-053).

## DT-052. La contabilidad de bases de tracto no la juzga nadie al publicar

**Gravedad: media. Cerrada (10/10/2026)**; ver «Cómo se cerró» al final. Hallada
el 10/10/2026 por la segunda revisión adversarial de DT-051, que la **reprodujo
ejecutando**. Decisión del responsable del 10/10/2026: se arregla en un cambio
propio justo después de DT-051.

`validate_published_tract_output_accounting` (`src/results.rs`), el juez que
recuenta las bases de tracto releyendo la verdad y el VCF, **no lo llama nadie en
`src/`** desde que entró (commit 12fd93c, 23/09/2026); solo lo llaman pruebas. El
juez completo (`validate_published_resources`), que es el que corre al publicar,
**acepta**:

- en el manifiesto general, la cota de bases de tracto de la verdad aumentada en
  1000 (con su total), las bases realizadas de la verdad disminuidas en 2 y, con
  VCF, las del VCF disminuidas en 1;
- en la configuración efectiva, la cota de la verdad aumentada en 1000, una cota
  por locus aumentada en 5 y `declared_max_tract_output_bases = 7` con 3814 bases
  realizadas.

Aunque se cablease tal como está, esa función solo ata el manifiesto y la efectiva
entre sí y con lo recontado, **no al TOML**.

**Corrección propuesta** (la probó la revisión en una copia, con la batería en
verde): llamar a esa función desde el juez completo; comparar la cota de la verdad
del manifiesto con la que resuelve el TOML (poniendo a cero solo lo medido), y la
cota de la efectiva (`repeat_output_budget`) con la que resuelve el TOML salvo las
partes del VCF, que se miden al publicar (DTC-097). REQ-VAL-62 remite a esta deuda
hasta que se cierre.

### Cómo se cerró (10/10/2026)

Los tres pasos propuestos, en `src/results.rs`, dentro del juez completo y
**después** de `validate_published_repeat_loci` (DT-051), que a su vez va después
del juez del abanico, del de la referencia final y del de la referencia simulada:
ninguno de ellos cambia de motivo.

1. **El juez que recuenta, cableado.** Su cuerpo pasa a
   `judge_tract_output_accounting`, que hace las mismas comprobaciones y
   devuelve además lo recontado; `validate_published_tract_output_accounting`
   queda como envoltorio público con la misma firma. El juez completo lo llama
   con la verdad y el VCF del catálogo (texto o comprimido); sin verdad, lo
   medido es cero. **Relee la verdad y el VCF a propósito**: el recuento va
   sobre los bytes, por un camino que no comparte nada con el lector tipado
   que acaba de juzgar la verdad; contar sobre el documento ya leído sería
   usar la estructura del productor (`realized_truth_tract_bases` recorre ese
   mismo tipo) y la comparación perdería su independencia. **Lee línea a
   línea**, la verdad con un lector con búfer y el VCF también sobre el lector
   BGZF (`published_vcf_lines`). La primera versión leía la verdad entera en una
   cadena y el VCF entero en memoria, y la revisión adversarial midió que eso
   subía el pico del juez un 38 %; el texto de entonces decía que las dos
   lecturas no coincidían en memoria, y **era falso**. Medido de nuevo el
   11/10/2026 con compilación optimizada, sobre los escenarios de la revisión
   (cuatro loci, VCF comprimido), ejecutando solo el juez completo sobre el
   paquete ya publicado, dos repeticiones:

   | Paquete | Antes del cambio (36d78eb) | Primera versión | Línea a línea |
   |---|---|---|---|
   | 400 individuos, verdad de 121 MB | 267,7 MB; 7,3 s | 369,3 MB; 8,1–8,2 s | 267,6 MB; 8,1 s |
   | 1200 individuos, verdad de 365 MB | 784,4 MB; 39,5–42,0 s | 1090,3 MB; 42,2–44,4 s | 784,4 MB; 42,9–43,3 s |

   El pico vuelve exactamente al de antes; el tiempo sube entre un 3 % y un 10 %
   (incluye también la comparación de DT-054, que lee la cabecera del VCF). En
   la ejecución entera no se nota en el pico, que lo fija la simulación (488 MB
   con 400 individuos, en las tres versiones; 1,41 GB con 1200, en la primera
   versión y en la de línea a línea).
2. **La cota del manifiesto, atada al TOML** (`validate_published_tract_bounds`).
   La contabilidad de tracto del manifiesto tiene que ser **exactamente** la que
   el productor escribiría con la configuración reconstruida del TOML copiado y
   las cifras **recontadas**: la cota de la verdad sale del TOML; la parte VCF y
   las bases realizadas, de lo recontado; los totales, de sumar. Como el juez
   que recuenta ya exige que el manifiesto coincida con la efectiva y con lo
   recontado, aquí solo llega una cota de la verdad inflada **a la vez** en los
   dos documentos.
3. **La cota de la efectiva, atada al TOML**, en la misma función:
   `repeat_output_budget` tiene que ser **exactamente** el que publica la
   configuración reconstruida, completado con lo que el VCF lleva, en total y
   por locus (`publish_measured_vcf_tract_bases` con lo recontado, DTC-097). No
   se neutraliza nada: la parte VCF y el total también se comparan, con lo
   medido. Eso cierra además un hueco que la propuesta no veía: **sin loci**, el
   manifiesto no publica sección y el juez que recuenta no mira nada, de modo
   que unas bases de VCF inventadas en la efectiva de un documento v4 sin loci
   se aceptaban; ahora tiene que ser cero entera.

**Pruebas** (`tests/results_repeat_guards.rs`, sección DT-052), todas con motivo
exacto y control positivo; ejecutadas contra una copia del commit anterior,
**cada uno** de los dieciséis casos se aceptaba:

- `the_tract_accounting_is_judged_when_publishing_and_bound_to_the_toml`: con y
  sin VCF, la cota de la verdad del manifiesto +1000 con su total y las bases
  realizadas de la verdad −2 (y, con VCF, las del VCF −1), que rechaza el juez
  que recuenta; en la efectiva, la cota de la verdad +1000, una cota por locus
  +5 y el tope declarado a 7, que rechaza la comparación con el TOML, y la cota
  de la verdad +1000 con su total, que rechaza el juez que recuenta; y la cota
  de la verdad +1000 con su total **en los dos documentos**, que rechaza la
  comparación del manifiesto con el TOML. Cada caso comprueba antes que el
  lector del manifiesto lo acepta (la aritmética sigue cuadrando), para que el
  rechazo sea del juez nuevo;
- `without_loci_the_effective_tract_bound_is_all_zero`: el hueco sin loci;
- `real_packages_pass_the_tract_and_amplification_comparisons`: lado aceptado
  con paquetes publicados de verdad —loci sin VCF, con VCF de texto, tracto
  insertado con VCF comprimido, amplificado con loci, con multiplicador de los
  homopolímeros, con loci y VCF, y sin loci—, comprobando que el VCF lleva
  bases de tracto donde se espera;
- añadidos tras la segunda revisión adversarial (11/10/2026): el caso «copias
  de la verdad por copia de haplotipo» en la primera prueba (sin él, no
  comparar ese campo sobrevivía a la batería entera), y
  `with_a_vcf_the_loci_are_bound_to_the_toml_before_the_tract_recount`, que con
  VCF exige el motivo de la comparación de los loci con el TOML para una
  secuencia de locus cambiada y para los identificadores de los dos loci
  intercambiados: si el recuento se adelantara a esa comparación, los dos se
  rechazarían por un síntoma del recuento y no por su causa, y las pruebas de
  DT-051, sin VCF, no lo veían.

El lado aceptado lo da además toda la batería: cada paquete que el programa
publica pasa por el juez completo, y ninguno cambió de resultado (2015 pruebas
antes del cambio, todas en verde con él).

**Verificación por mutación** (copia fuera del repositorio, con su propio
directorio de compilación, un mutante por ejecución contra la batería entera,
huella antes, mutada y restaurada, y `touch` tras restaurar). Murieron los nueve
(dos primeras redacciones del de «no recontar» no compilaban y no cuentan):
anular la comparación del manifiesto (la mata solo el caso coherente en los dos
documentos) o la de la efectiva (la matan los casos de la efectiva y el de sin
loci); no llamar al juez que recuenta; rehacer la efectiva sin lo medido o el
manifiesto sin lo recontado; no pasarle el VCF; no devolver las bases del VCF por
locus; invertir la condición de versión; y, para DT-043, anular la comparación de
la amplificación (la mata solo `the_amplification_of_the_event_manifest_is_the_effective_one`).
Las que rehacen sin lo medido, no pasan el VCF o no recuentan las matan además el
lado aceptado y decenas de pruebas que publican paquetes con loci.

Tras la segunda revisión adversarial (11/10/2026), dos mutantes más que entonces
sobrevivían a la batería entera, ahora muertos, cada uno solo por su prueba nueva
(huella antes y restaurada `5a2b3d5a…7a18`): no comparar las copias de la verdad
por copia de haplotipo en la efectiva (mutada `aff30acc…cc87`; lo mata
`the_tract_accounting_is_judged_when_publishing_and_bound_to_the_toml`) y
adelantar el recuento a la comparación de los loci con el TOML (mutada
`555fbc51…b3b1`; lo mata
`with_a_vcf_the_loci_are_bound_to_the_toml_before_the_tract_recount`).

La puerta de cobertura de los jueces: `src/results.rs` pasa de 783 a 786
comprobaciones y sigue con 498 sin disparar; las tres nuevas se disparan.
Con la segunda ronda (11/10/2026), 787 comprobaciones y 497 sin disparar: entran
las dos de DT-054, que se disparan, y la relectura del VCF por flujo pasa de
cuatro puntos de rechazo a tres con los mismos motivos (de ellos solo se dispara,
antes y ahora, el del bloque BGZF que declara más de 64 KiB). Las mutaciones de
la segunda ronda se hicieron antes de unir en un solo punto los tres motivos de
la lectura de una línea del VCF, que no tocan ninguna de las comprobaciones
mutadas; las medidas de memoria y tiempo, también.

La regla queda en REQ-VAL-62, que ya no remite a esta deuda.

## DT-053. Las copias realizadas de cada alelo no se atan al sorteo rehecho con la semilla

**Gravedad: baja.** Hallada el 10/10/2026 por la segunda revisión adversarial de
DT-051, que la **reprodujo ejecutando**. **Sin arreglar** (decisión del
responsable: no arreglarla por ahora).

Desde DT-051 las copias realizadas del manifiesto se atan a lo que la verdad
lleva, no al sorteo. Con un TOML copiado falsificado de forma coherente
—frecuencias 0,6/0,4 cambiadas a 0,9/0,1, con su huella, el identificador de
ejecución y los documentos rehechos— el juez acepta las copias 4/6 que el paquete
publica, cuando el sorteo verdadero con esas frecuencias da 10/0. Tampoco
distingue el reparto entre dos alelos gemelos (la limitación anotada en DT-051).

**Corrección posible:** rehacer el sorteo con la semilla y comparar las copias por
índice de alelo. El juez ya lo rehace para las estadísticas de la población
inicial y lo descarta; aprovecharlo cerraría también el caso de los gemelos.

## DT-054. Nadie compara el `execution_id` de la verdad y del VCF con el del manifiesto general

**Gravedad: media. Cerrada (11/10/2026)**; ver «Cómo se cerró» al final. Hallada
el 11/10/2026 por la revisión adversarial de DT-052, que la **reprodujo
ejecutando**. Decisión del responsable: arreglarla en el mismo cambio.

El juez completo ata el `execution_id` del manifiesto general a la semilla y a
los bytes del TOML copiado, pero nadie comparaba con él el que llevan las
cabeceras de la verdad (`execution_id`) y del VCF (`##gps_execution_id=`).
Reproducción: en el TOML copiado se añade `max_tract_output_bases = 999999`, se
rehace su huella en el manifiesto, se pone el mismo tope en la efectiva y se
sustituye el `execution_id` del manifiesto por el que sale del TOML nuevo. El
juez lo **aceptaba**, aunque la cabecera de la verdad seguía con el
identificador de la ejecución original: el paquete afirmaba una configuración
que no es la que produjo su verdad. El manifiesto de evento no publica
`execution_id` (comprobado en el productor y en su esquema v3), de modo que no
hay nada que comparar allí.

### Cómo se cerró (11/10/2026)

En `src/results.rs`, función nueva `validate_published_execution_ids`, llamada
por el juez completo después de todo lo de DT-051 y DT-052:

- **verdad**: el `execution_id` de la cabecera, que ya lee y juzga el juez de la
  verdad (`validate_physical_truth` lo devuelve ahora junto al recuento de
  copias), tiene que ser el del manifiesto general;
- **VCF**, si lo hay: la línea `##gps_execution_id=` de la cabecera, leída por
  flujo y solo hasta la primera línea que no es de metadatos, tiene que existir
  y llevar el mismo identificador.

Los dos rechazos pasan por `manifest_error`, con `#[track_caller]`.

Antes de exigirlo se pasó la batería entera: **todos los paquetes que publica el
programa lo cumplen**. Fallaron solo las dos pruebas de
`tests/f4_physical_technical_manifest.rs`, cuya fixture **no la publica el
programa**: publica sin evento y después cambia a mano el TOML copiado y rehace
el `execution_id` del manifiesto, que es justo el patrón que esta deuda cierra.
La fixture rehace ahora también la cabecera de la verdad, con su huella, como
tendría que hacer quien quisiera colarlo de forma coherente.

**Prueba**: `the_truth_and_the_vcf_carry_the_execution_id_of_the_manifest`
(`tests/results_repeat_guards.rs`), con la reproducción de la revisión y el
motivo exacto: sin VCF la rechaza la comparación de la verdad; con VCF y la
verdad sin rehacer, también la de la verdad, que va antes; con la cabecera de la
verdad rehecha (y su huella), la del VCF. Falla contra el código anterior
(ejecutada sobre una copia de 36d78eb: el juez la aceptaba).

**Verificación por mutación** (copia fuera del repositorio, batería entera, huella
antes y restaurada `5a2b3d5a…7a18`): anular la comparación de la verdad (mutada
`76cc2f52…549e`) y anular la del VCF (mutada `e0fd338b…bded`) mueren los dos, y
solo por esta prueba.

## DT-055. El juez de la pareja no comprueba que los genotipos realizados cuenten solo observaciones completas

**Gravedad: baja.** Hallada por la revisión adversarial de la tanda 8 de DT-001
(12/10/2026). La arquitectura (`design/arquitectura-0.1.md`, F2.3.5) fija que las
estadísticas iniciales cuentan en los genotipos realizados **solo observaciones
completas**. De ahí salen dos desigualdades por locus que el juez de la pareja
efectiva/manifiesto (`validate_initial_population_statistics_locus` en
`src/results.rs`) **no comprueba**:

1. por alelo, la suma de dosis por recuento realizado de cada genotipo no supera el
   recuento realizado de ese alelo;
2. las observaciones genotípicas completas por la ploidía no superan las copias
   observadas.

Las caza después el juez completo, que rehace las estadísticas desde el TOML, la
semilla y la matriz y las compara con las publicadas (`validate_physical_initial_population_statistics`),
de modo que un paquete publicado de verdad nunca las incumple y un paquete
manipulado no se acepta. Pero un documento que solo vea el juez de la pareja sí
puede incumplirlas. Decisión del responsable: **registrarlo, no añadir las
comprobaciones**. La prueba `tests/results_population_snvs_guards.rs` no consagra
ningún documento que las incumpla: la base y cada caso «aceptado» pasan por
`assert_architecture_inequalities`.
