# Acta de cierre de F5 — F5.6.3

**Fecha:** 13 de septiembre de 2026.
**Alcance:** ejecutar todas las puertas, guardar su salida real, auditar las 46
filas de la matriz de trazabilidad cuyo hito proyecta sobre F5 y recontar las
deudas que F5 deja declaradas.
**Entorno:** WSL2 sobre ext4, Ubuntu 24.04, Rust 1.97.1, edición 2024, target
`x86_64-unknown-linux-gnu`. Compilación dentro del repositorio.

> Este acta no cierra el hito por sí sola: la casilla de F5.6.3 sigue en ☐ hasta
> que la Revisión F5.6 la apruebe. La auditoría destapó **un defecto material**,
> descrito y resuelto en la §5.

## 1. Recuento mecánico de la proyección F5

La proyección son **46 filas**, con **46 REQ distintos**, e incluye
expresamente `REQ-VAL-54`, `REQ-MAN-70` y `REQ-SAL-02`.

El recuento es mecánico y vive en
`tests/f5_6_closure.rs::the_f5_projection_of_the_matrix_is_exactly_forty_six_unique_requirements`,
que relee `design/trazabilidad-requisitos-0.1.md`, parte cada fila por sus
cuatro columnas y proyecta el campo Hito **término a término** —el campo admite
varios hitos separados por `/`, de modo que `F2.0/F5.6` cuenta y `F50.1` no—.
No hay ninguna lista de 46 identificadores escrita a mano que alguien tenga que
mantener en paralelo: la lista **es** la tabla.

Tres verdades independientes sostienen la cifra y ninguna se compara consigo
misma:

1. el literal del contrato, `DECLARED_F5_ROWS = 46`, que es lo que el plan pide;
2. la prosa de la §1 de la propia matriz, que se mantiene a mano;
3. la prosa de la entrada F5.6.3 del plan, que también se mantiene a mano.

Si la tabla cambia sin que cambien las dos prosas —o al revés— la prueba muere.
Comprobado con cuatro mutaciones verificadas: borrar la fila de `REQ-POB-01`
(46→45, falla), cambiar la prosa de la matriz de 46 a 47 (falla), citar una
prueba que no existe (falla) y devolver una fila a «Planificada» (falla).

El proyector se prueba aparte, en
`tests/f5_6_closure.rs::the_projection_discriminates_a_changed_matrix`: fija el
criterio término a término sobre siete aceptaciones y siete rechazos, y sobre
una tabla sintética demuestra que **añadir, quitar o reasignar una fila F5 mueve
el recuento**. Sin ella, un proyector por subcadena contaría `F50.1` y uno por
campo completo perdería `F2.0/F5.6`, y el recuento de 46 sería una coincidencia.

## 2. Auditoría fila a fila

Las 46 filas estaban **todas** en «Planificada» —salvo `REQ-VAL-54`, en «Parcial
F2.0 completada»— pese a que F5.0 a F5.6.1 están implementadas. La matriz
describía la evidencia *prevista*, no la observada, y por tanto no servía como
gate de cierre de nada.

Este hito sustituye la columna de evidencia de las 46 filas por lo que F5 cerró
realmente, con citas a pruebas concretas.
`tests/f5_6_closure.rs::every_f5_row_states_what_f5_did_and_cites_tests_that_exist`
exige, fila a fila, que la evidencia empiece por «Implementada» o «Parcial»,
nombre F5, no empiece por «Planificada» y cite al menos una prueba; y que
**cada** símbolo citado exista y sea una prueba de verdad —`#[test]` en Rust,
`def test_…` en Python—. El lector de citas se prueba a su vez en
`tests/f5_6_closure.rs::the_citation_reader_discriminates_tests_from_helpers`,
que distingue una prueba de una función auxiliar y de una prueba anterior: sin
él, una cita mal escrita pasaría por no encontrarse nada.

Las filas compartidas con hitos posteriores conservan expresamente su parte
pendiente y siguen en «Parcial»: `REQ-EVO-03` (F7.6), `REQ-EVO-43` (F7.1),
`REQ-EVO-44` (F7.7), `REQ-EVO-49` (F6.3), `REQ-EVO-78` (F6.6) y `REQ-SAL-02`
(F8.5). Ninguna asignación de REQ a hitos cambia y ninguna fila ajena entra: el
total de filas `REQ` de la matriz sigue siendo 844.

### Dónde quedan los lectores, el oráculo, los vectores y la paridad

El gate de cierre pide nombrarlos. No se reimplementan aquí: viven en las
baterías de F5.1 a F5.6.1 y ahora están **citados desde la matriz**, de modo que
la prueba de citas comprueba que existen y que son pruebas.

| Pieza | Dónde vive | Fila que la cita |
|---|---|---|
| Lector JSON acotado del par de autoría | `tests/f5_bounded_input.rs::fragmented_reader_reaches_exact_eof_and_hashes_the_same_literal_bytes`, `tests/f5_bootstrap_dispatch.rs::dispatcher_classifies_discriminators_after_4096_bytes_and_preserves_non_f5` | REQ-EVO-74 |
| Lector TSV independiente de los siete artefactos | `tests/f5_artifacts.rs::every_isolated_corruption_is_detected_by_an_independent_reader` | REQ-VAL-56 |
| Lector JSONL de `truth.jsonl` v2 | `tests/f5_oracle.rs::the_man70_recount_needs_only_the_founder_truth`, `tests/f5_artifacts.rs::every_published_truth_line_validates_against_the_external_schema` | REQ-MAN-70, REQ-SAL-02 |
| Oráculo independiente del constructor | `tests/f5_oracle.rs::{the_oracle_accepts_a_whole_published_package_and_summarises_it, every_isolated_byte_mutation_dies_of_its_own_cause}` | REQ-VAL-56 |
| Vectores congelados RNG e ID | `tests/f5_rng_catalog.rs::f5_rng_catalog_v1_has_literal_subseed_and_chacha12_vectors`, `tests/f5_founders.rs::the_independent_id_oracle_reproduces_the_frozen_catalog_vectors` | REQ-EVO-77 |
| Paridad Marimo ↔ TOML manual ↔ Rust | `tests_ui/test_f5_6_1_scenario_ui.py::{test_the_interface_publishes_the_same_effective_manifest_and_digests, test_a_rejected_pair_produces_the_same_diagnostic_through_both_ways}` | REQ-VAL-54, REQ-UI-15 |

### Hueco encontrado y cerrado: la cota de ploidía de F5

`REQ-CONF-65` exige rechazar una ploidía no entera, menor que 1 o mayor que 10
identificando el problema. La guarda existe —`src/adapters.rs`, en la admisión
del general v3— pero **no tenía ninguna prueba**: la batería de F5.0.1 solo
sustituía `ploidy = 2` por `ploidy = 1`, que es un valor legal, y la
construcción del efectivo F5 se limita después a una conversión a `u8`, que
admitiría hasta 255.

`tests/f5_6_closure.rs::the_f5_admission_rejects_a_ploidy_outside_one_to_ten`
cubre el hueco: acepta los dos extremos legales —1 y 10, la cota es inclusiva
por los dos lados— y rechaza 0, 11, 250 y 4096 exigiendo `GPS-CONF-101`, el
campo `ploidy` y el detalle `debe pertenecer a 1..=10`. **Discrimina**:
sustituida la guarda por la conversión a `u8` que la seguía, la prueba falla.

## 3. Búsquedas negativas de alcance, con control positivo

`tests/f5_6_closure.rs::no_f5_source_or_schema_names_a_panel_sv_repeat_alignment_or_recombination`
recorre los quince módulos `src/f5_*.rs` y los seis esquemas F5/v5 buscando doce
términos: `panel`, `structural_variant`, `tandem_repeat`, `microsatellite`,
`bam`, `cram`, `recombination`, `recombinant`, `crossover`, `meiosis`, `mosaic` y
`linkage_map`. Ninguno aparece.

La guarda es léxica y lleva tres controles positivos, porque una búsqueda mal
escrita pasaría siempre: que **encuentra** el término cuando está en el código,
que **no** lo encuentra cuando solo está en la prosa que lo prohíbe —que es el
modo real en que los módulos F5 nombran lo que no hacen, y por eso la búsqueda
retira los comentarios antes de mirar— y que no lo encuentra en un comentario
al final de una línea de código. **Discrimina**: plantado `fn panel_stub()` al
final de `src/f5_ld.rs`, la prueba falla.

La mitad de interfaz la cubre
`tests_ui/test_f5_6_1_scenario_ui.py::test_no_panel_supplied_haplotypes_recombination_or_map_control_exists`,
y la de los módulos de fundadores y clausura,
`tests/f5_founders.rs::the_module_never_names_a_haplotype_process_a_cross_a_mosaic_or_a_map`
y `tests/f5_closure.rs::the_closure_module_names_no_migration_selection_or_genealogy`.

### Ninguna vía de publicación de F5.6 llama a los publicadores desnudos

`tests/f5_6_closure.rs::no_f5_6_publication_path_calls_the_bare_publishers`
extiende la obligación de proceso de F5.6.1 —que cubría `application.rs` y
`main.rs`— a **todos** los módulos del crate menos `f5_closure.rs`, que es el
que envuelve. La búsqueda descarta lo precedido por `seal_and_` —el nombre del
envoltorio contiene el del envuelto— y lo precedido por `fn `, que es la
declaración de `src/f5_publication.rs` y no una llamada. Tres controles
positivos: el orquestador sí nombra los dos envoltorios, la misma búsqueda sí
encuentra las llamadas desnudas donde viven, y distingue envoltorio, declaración
y llamada sobre tres cadenas literales. **Discrimina**: plantado
`publish_f5_failure(workspace)` en `src/application.rs`, la prueba falla.

## 4. Los dos auditores compuestos sobre el destino ya renombrado

La mitad «éxito» la fijó F5.6.1 y sigue verde:
`tests/f5_cli_run.rs::the_productive_chain_publishes_and_both_auditors_count_the_package`
—`audit_f5_package` como recuento MAN-70 sobre la raíz de carga útil y
`audit_f5_identity_closure` sobre el destino ya renombrado— y
`tests/f5_cli_run.rs::the_orchestrator_composes_both_audits`.

La mitad «fallo» la cierra este hito en la parte que F5.6.1 dejó abierta:
`tests/f5_6_closure.rs::the_failure_package_auditor_accepts_the_renamed_destination_and_dies_of_each_mutation`
publica un fallo con la cadena productiva real —una cuota controlada que su
tolerancia declarada de `0` no sostiene, `GPS-DOMAIN-103`, etapa `founders`, sin
ninguna instrumentación—, comprueba que el destino está **renombrado** y que no
queda provisional vivo, reconstruye la `F5FailureDeclaration` desde el
**manifiesto** y se la da a `audit_f5_failure_package`, que relee
`outputs/f5-failure.json`, que es otro documento. No es circular: el auditor
contrasta dos documentos independientes y reabre los bytes del enlace.
**Discrimina** con tres mutaciones dentro de la propia prueba: una etapa que el
documento no sostiene, un enlace que no describe sus bytes y un artefacto
biológico intruso en el árbol publicado; retirado el intruso, vuelve a aceptar.

El gemelo `audit_f5_failure_closure` ya lo componía el orquestador desde
F5.6.1, y `audit_f5_failure_package` queda compuesto **ahora**: la §5 explica
por qué no podía estarlo antes.

## 5. Defecto encontrado y corregido: la cota de propuestas del auditor de fallo

**`audit_f5_failure_package` rechazaba el paquete que la propia cadena
productiva publica cuando la búsqueda haplotípica se agota.** El defecto era de
**F5.4.2** y solo podía aparecer al componer ese auditor sobre un agotamiento
**real**: las baterías de F5.4 construían sus declaraciones con contadores
pequeños, de modo que ninguna revisión por hito lo vio. Es exactamente la clase
de defecto para la que existe un hito de cierre.

### Qué estaba mal

`src/f5_oracle.rs` comparaba los tres contadores contra `limits` **posición a
posición**, de modo que a `swap_proposals` le aplicaba `max_swap_proposals`. El
contrato dice otra cosa, y lo dice en dos sitios:

- **DTC-089 l. 1196:** `max_swap_proposals` es el presupuesto de **un intento**
  —«Un intento termina … o al consumir `max_swap_proposals`»—.
- **DTC-089 l. 1847-1848:** enumera los cuatro contadores juntos, y
  `haplotype_attempts_used` es **acumulado** sobre los intentos. Leer el tercero
  como «de un intento» rompe la coherencia del propio grupo.

El productor hacía lo correcto. Reproducido con el binario: un `run` con
`max_swap_proposals = 256` y `max_haplotype_attempts = 32` publica
`counters.swap_proposals = 8192`, y el auditor lo rechazaba con `GPS-RES-101`.
Que los tres contadores compartieran bucle posicional es lo que hizo invisible
el desajuste: la correspondencia índice-a-índice parecía una simetría y era una
coincidencia.

### La cota correcta, literalmente

    swap_proposals <= max_swap_proposals * haplotype_attempts_used

`swap_proposals` sale del bucle posicional y se comprueba aparte; los otros dos
sí son contador contra su límite homónimo. El producto usa `checked_mul` —el
techo contractual es `10^7 * 10^4 = 10^11` y cabe en `u64`, pero un
desbordamiento silencioso aquí abriría la cota en lugar de cerrarla— y el
porqué queda escrito en el código citando las dos líneas de DTC-089, para que
nadie lo «simplifique» de vuelta.

La cota es **ajustada, no holgada**, y por una razón estructural: el productor
acumula los dos contadores en el **mismo recorrido** —un intento sumado y sus
propuestas sumadas por cada registro de la unidad que el diagnóstico nombra—,
de modo que `swap_proposals` es la suma de `haplotype_attempts_used` sumandos y
ninguno puede pasar de `max_swap_proposals`. Un fallo que no nace de la búsqueda
deja los dos contadores en cero y la cota sigue valiendo.

### Discriminación, por los dos lados

`tests/f5_6_closure.rs::the_failure_package_auditor_accepts_an_exhausted_search_package`
exige las dos mitades, y se ha comprobado que muere con cada una de las dos
formas de equivocar la cota:

| Cota en `src/f5_oracle.rs` | Qué hace la prueba |
|---|---|
| `max_swap_proposals` (la vieja, sin producto) | **falla** al aceptar: `swap_proposals=8192 excede max_swap_proposals=256 por haplotype_attempts_used=32` |
| `u64::MAX` (demasiado laxa, no rechaza nada) | **falla** al rechazar: «una propuesta por encima de intentos × presupuesto tiene que rechazarse» |
| `max_swap_proposals * haplotype_attempts_used` | pasa |

La mitad que cierra se ejerce mutando los **bytes del paquete**, no la
declaración: el auditor compara antes declaración y documento término a
término, así que tocar solo la declaración moriría de esa otra causa y no diría
nada de la cota. La prueba reescribe `swap_proposals` en `f5-failure.json` **y**
en `manifest.json`, rehace tamaño y SHA-256 del enlace para que el paquete siga
siendo coherente consigo mismo, y exige que el rechazo sea `GPS-RES-101` con un
detalle que nombre `swap_proposals` y `haplotype_attempts_used`.

### Y el auditor queda compuesto en el orquestador

Era la mitad que F5.6.1 no pudo componer. Con la cota corregida ya es posible.
El impedimento que quedaba —la `F5FailureDeclaration` que el orquestador no
conserva— se resuelve **sin** reutilizar el informe en memoria, que haría la
auditoría autocertificante: `audit_published_f5_failure` la reconstruye desde el
**manifiesto publicado**, de modo que el oráculo contrasta dos documentos del
paquete y reabre los bytes del enlace. Corre sobre el destino **ya renombrado**,
después de `audit_f5_failure_closure`.

La composición se demuestra por **comportamiento**, no solo por el texto de la
fuente:
`tests/f5_6_closure.rs::the_orchestrator_composes_the_failure_package_auditor`
exige que un agotamiento real devuelva `GPS-DOMAIN-104` —el diagnóstico del
informe— y no el de ningún auditor. Discrimina porque el agotamiento es
justamente el caso que la cota vieja rechazaba: devuelta la cota vieja, el
orquestador pasa a devolver `GPS-RES-101`, lo que solo puede ocurrir si el
auditor está compuesto y se alcanza. Comprobado.

Con esto, los **cuatro** auditores están compuestos en producción:
`audit_f5_package` y `audit_f5_identity_closure` en la mitad «éxito»,
`audit_f5_failure_closure` y `audit_f5_failure_package` en la de «fallo».

## 6. Atomicidad en frontera de proceso — deuda heredada, cerrada aquí

F5.4 y F5.5 demostraron `RENAME_NOREPLACE` sin *fallback*, `OpenArtifact::Drop`,
el marcador `INCOMPLETE` en el envoltorio y la revalidación previa al
renombrado **sobre la frontera de biblioteca**, con fixtures dentro del proceso
de prueba. Las únicas pruebas de interrupción y concurrencia entre procesos
recorrían F1. F5.6.1 cableó el primer `run` F5 productivo, de modo que ahora sí
es posible repetirla sobre el proceso. Dos pruebas nuevas la cierran:

- `tests/f5_6_closure.rs::an_interrupted_f5_run_leaves_its_provisional_and_creates_no_destination`
  lanza un `run` F5 real, espera a que aparezca su provisional marcado,
  **comprueba que el proceso sigue vivo** y lo mata con `SIGKILL`, que no
  ejecuta ningún código de limpieza. Lo que queda es el provisional con su
  marcador `INCOMPLETE` y permisos `0700`, ningún destino y ningún byte
  biológico fuera del árbol marcado.
- `tests/f5_6_closure.rs::two_concurrent_f5_runs_have_one_winner_and_never_clobber_the_destination`
  comprueba primero que un destino que ya existe se rechaza con `GPS-IO-003` sin
  tocar su contenido, y después lanza **dos** `run` F5 solapados sobre el mismo
  destino libre: exige que los dos provisionales coexistan antes de que ninguno
  renombre —el solape es real y se comprueba, no se supone—, que exactamente uno
  publique, que el perdedor muera con `GPS-IO-003` conservando su provisional
  marcado y que el paquete del ganador esté completo.

Las dos usan un par de autoría dimensionado para que el `run` dure segundos, de
modo que la ventana no dependa de la carga de la máquina.

**Discriminan**, con la mutación verificada en cada caso:

| Mutación en `src/publication.rs` | Prueba que muere |
|---|---|
| `RenameFlags::NOREPLACE` → `RenameFlags::empty()` | la carrera: los dos procesos publican |
| permisos efectivos del provisional `0700` → `0755` | la interrupción: `493 != 448` |
| el marcador `INCOMPLETE` pasa a llamarse `PLACEHOLDER` | la interrupción: no hay provisional marcado |

La primera mutación de permisos que se intentó —cambiar solo el `builder.mode`—
**no mutó**: una llamada posterior a `set_permissions` la revertía y la prueba
seguía en verde. Hizo falta mutar las dos escrituras. Queda escrito porque es
exactamente el falso superviviente contra el que el encargo avisa.

**Lo que no cubre:** DTC-045 pide la comprobación en WSL2/ext4 **y en Ubuntu
nativo**. Aquí solo se ha ejecutado en WSL2/ext4. La mitad nativa sigue
pendiente y no la puede dar por hecha este acta.

## 7. Deudas recontadas

| Deuda | Origen | Estado a 13/09/2026 |
|---|---|---|
| **Dos piezas fuera de `f5-resource-accounting-v8`**: el oráculo del paquete y la construcción/comparación de los dos documentos, más el lector del paquete | Revisión F5.4 / F5.4.3 | **Abierta, y ya no es teórica.** El dictamen de la revisión independiente sigue siendo exacto: **es v8 el incompleto**, porque `final_formula_peak` presupuestó el escritor por flujo de F5.4.1 y nunca modeló ni la construcción de los dos documentos —`O(B*N_pop)`— ni un lector del paquete —acotado por el registro más largo—. La recomendación sigue siendo **v9 con dos términos nuevos**, y su coste sigue siendo que **reabre los dos esquemas** que llevan el literal `"f5-resource-accounting-v8"` como `const`. El plan la declaró «bloqueante en el primer `run` productivo»: ese `run` existe desde F5.6.1 y **compone `audit_f5_package`**, de modo que la cadena productiva ejecuta hoy una pieza que su propio contrato de recursos no modela. F5.6.1 midió que las tres piezas caben con holgura en el presupuesto sellado de 1 GiB en los tres tamaños probados (7,4–8,1 MiB RSS). Caber no es estar declarado. |
| **Desviación de hasta 2× del presupuesto sellado en la publicación completa** | Revisión F5.4, «elevada como bloqueante de F5.6» | **Abierta, misma raíz que la anterior.** Acotada y declarada por `tests/f5_publication.rs::the_whole_publication_stays_within_twice_the_sealed_budget_without_the_oracle`, y medida en el primer `run` productivo. Se resuelve con v9 o no se resuelve. |
| **Atomicidad demostrada en frontera de biblioteca, no de proceso** | Revisión F5.4 / F5.5, con destino F5.6 | **Cerrada en WSL2/ext4** por las dos pruebas de la §6. **Pendiente en Ubuntu nativo.** |
| **`audit_f5_failure_package` sin composición en el orquestador** | F5.6.1, con destino F5.6.3 | **Cerrada.** Compuesta en producción sobre el destino ya renombrado, con la declaración reconstruida desde el manifiesto publicado. Antes era imposible por el defecto de la §5. Los cuatro auditores están ahora compuestos. |
| **DTF-058 — perfil humano de LD `human_ld_ragsdale_approx_v1`** | F5.6.2, aplazada el 13/09/2026 | **Aplazada, y no bloquea el cierre.** Registrada en `design/registro-ampliaciones-futuras-0.1.md`. El motivo que decide sigue siendo exacto: la fuente no cubre el rango de distancias que el perfil necesita. No existe `resources/f5/`, no hay resolución `builtin:` en `src/` y `kappa` no aparece en ninguna ruta de cálculo. |
| **Superposición acotada entre el oráculo de F5.4.2 y la batería del escritor de F5.4.1** | F5.4.2, «F5.6.3 puede retirarla» | **No retirada, deliberadamente.** El plan dice «puede», no «debe». Retirarla exigiría reescribir `tests/f5_artifacts.rs`, que es la evidencia propia del **escritor** —el oráculo juzga el paquete, aquella juzga a quien lo escribe— y está fuera del alcance de este encargo. No hay defecto: hay solapamiento acotado. |
| **Docstring de F5.3 que anuncia `GPS-CONF-101` donde el código devuelve `GPS-INTERNAL-101`** | F5.5, anotado para quien reabra F5.3 | **Abierta y anotada.** No se toca por la prohibición de revertir cambios ajenos. |
| **DTC-066: `bcftools`/`tabix` 1.24 externo** | F2.5 | **Abierta y ajena a F5.** Sostiene una de las 13 pruebas ignoradas. |

## 8. Puertas, con su salida real

### Rust — las seis, todas con `--locked` donde el contrato lo pide

| Puerta | Resultado |
|---|---|
| `cargo fmt --all -- --check` | exit 0, sin diferencias |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | exit 0, **0 avisos** |
| `cargo test --workspace --all-targets --locked --no-fail-fast` | exit 0 — **1098 pasan, 0 fallan, 12 ignoradas**, 0 avisos |
| `cargo test --workspace --doc --locked` | exit 0 — 2 pasan, 0 fallan, 0 ignoradas |
| `cargo build --release --locked --target x86_64-unknown-linux-gnu` | exit 0, 0 avisos, 1 m 39 s |
| `cargo doc --no-deps` | exit 0, **0 avisos** |

Línea base al empezar el hito, medida en el mismo árbol: **1085 / 0 / 12**, 0
avisos. El delta es **+13 pruebas que pasan y ninguna ignorada nueva**: las
trece que añade `tests/f5_6_closure.rs`, todas verdes. El recuento de ignoradas
vuelve a **12** porque el defecto de la §5 se corrigió en vez de aplazarse.

### Python — dir a dir

| Puerta | Resultado |
|---|---|
| `MARIMO_SKIP_UPDATE_CHECK=1 uv run --frozen marimo check --strict ui/gps_app.py` | exit 0, sin salida |
| `MARIMO_SKIP_UPDATE_CHECK=1 uv run --frozen --with pytest python -m pytest tests` | exit 0 — **7 pasan, 1 se omite** |
| `MARIMO_SKIP_UPDATE_CHECK=1 uv run --frozen --with pytest python -m pytest tests_ui` | exit 1 — **94 pasan, 2 fallan** |

### Los dos fallos de `tests_ui`, con su causa y su fase

Son **preexistentes y ajenos a F5**: fallaban antes de F5.6.1 y este hito no los
toca. No se esconden.

1. `tests_ui/test_f2_5_vcf_ui.py::F25VcfUiRouteTests::test_ui_toml_routes_configuration_to_truth_and_optional_vcf`
   — **fase F2.5**. El TOML que la interfaz exporta con `vcf_output_mode="bgzf"`
   ya no valida: el binario responde `GPS-CONF-008`,
   `component = "vcf_output,initial_populations"`, detalle «truth.jsonl no
   contiene copias ni genotipos cuando `initial_populations` está vacío». El
   fixture de la prueba no declara población inicial y una regla posterior a
   F2.5 lo exige.
2. `tests_ui/test_ui_contract.py::UiContractTests::test_ui_technical_event_v2_round_trips_to_rust_validation`
   — **fase F3, desplazada por F4.1**. La prueba espera el texto «F4.1 rechaza
   la ejecución» y el binario devuelve `GPS-CONF-008`,
   `component = "direct_events,population_snvs"`, detalle «la matriz SNV falta
   después de la validación efectiva». El rechazo sigue ocurriendo; el mensaje
   que la prueba fija es el de otra redacción.

Ninguno de los dos es un fallo de F5 ni una regresión de este hito, y ninguno se
corrige aquí: están fuera del alcance del encargo.

## 9. Lo que este hito deja abierto

1. **La mitad nativa de la atomicidad de proceso** (Ubuntu nativo), que DTC-045
   pide y aquí solo se ha ejecutado en WSL2/ext4. **No la tenemos.** Destino: la
   Revisión F5.6 decide si la exige antes de F6 o la traslada al primer hito que
   disponga del entorno.
2. **`f5-resource-accounting-v9`**, con el dictamen y la recomendación intactos
   y ahora con consumidor productivo real. La composición del cuarto auditor
   añade a la rama de fallo una relectura del manifiesto publicado, acotada por
   ese documento, que v8 tampoco modela: es la misma deuda, no una nueva.
3. **Los dos fallos preexistentes de `tests_ui`**, de F2.5 y F3.
4. **DTF-058**, aplazada por decisión de la persona usuaria.

El defecto de la §5 **queda resuelto**, no aplazado. Nada de lo anterior es una
regresión introducida por F5.6.3 y nada estaba oculto.
