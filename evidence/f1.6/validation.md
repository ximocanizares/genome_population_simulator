# Validación de F1.6

- Fecha: 15 de agosto de 2026.
- Plataforma: Ubuntu 24.04 sobre WSL2, x86-64, ext4.
- Estado de la ronda: H-F16-R25–R27 están 3/3 aplicados con adaptación. Sus
  regresiones, las 196 pruebas y todas las puertas, incluidos validadores
  externos y OCI, son verdes. Se restauran F14, MAN-90, SAL-27, F1.6 22/22 y
  F1 83/83 para una nueva revisión independiente; la pausa sigue abierta.

## Pruebas específicas

- `reference_reduction_config`: 4 casos de autoría, literal/ubicación,
  exclusividad, dominios, objetivo exacto y subsemilla.
- Módulo `coordinates`: 4 casos de bordes, ida/vuelta, espacios incompatibles y
  rechazo sin *clipping*.
- Módulo `exact_decimal`: 6 casos; equivalencias/exponentes, redondeo único,
  frontera exacta que diverge de `f64`, precisión superior a `u128` y dominios
  abiertos, más gramática TOML exacta antes de retirar separadores.
- Módulo `reference_reduction`: 10 casos; suma/residual, múltiplo sin residual, preflight sin consumo
  RNG, orden de materialización, dos oráculos exhaustivos de huecos, prueba
  estadística, bordes y dos vectores congelados.
- `f1_6_acceptance`: 5 recorridos públicos; bases y porcentaje,
  misma/diferente semilla, FASTA mixto con `N`, reconstrucción ajena, corrupción
  del mapa, disposición inviable, omisión observable, aislamiento causal y
  categorías con cero/una/dos reducciones. El aislamiento añade una regla
  efectiva `output_enabled` y compara subsemilla, mapa y FASTA.
- El recorrido H-F16-R06 conserva los positivos `positive-complete` y
  `positive-failed`. Python/jsonschema acepta como estructuralmente válidas diez
  mutaciones semánticas: `wrong-fasta-digest`, `wrong-map-path`, `wrong-source`,
  `wrong-derived-origin`, `orphan-reduced-reference`,
  `reused-reference-fasta`, `incompatible-origin-kind`,
  `crossed-derived-fasta-id`, `crossed-map-id` y
  `crossed-derived-reference-id`.
- El validador Rust y un oráculo semántico independiente que reinterpreta los
  bytes releídos rechazan las diez mutaciones semánticas. El oráculo de prueba
  no llama al productor, al validador Rust, ni a Python/jsonschema; cruza por su
  cuenta catálogo, referencias, reducciones, unidades, origen y propiedad de
  FASTA. Ambos oráculos aceptan también los dos positivos.
- H-F16-R13 añade las regresiones
  `execution::tests::global_manifest_links_must_match_catalog_exactly_without_reductions`
  y
  `execution::tests::global_manifest_links_match_catalog_for_complete_failed_zero_one_and_two_reductions`.
  El segundo recorrido ejecuta seis positivos: `positive-zero-complete`,
  `positive-zero-failed`, `positive-one-complete`, `positive-one-failed`,
  `positive-two-complete` y `positive-two-failed`.
- El mismo corpus ejecuta 29 mutaciones individuales. Sobre los cinco contextos
  cero-completo, uno-completo, uno-fallido, dos-completo y dos-fallido prueba
  `zero-complete-reference-digest`, `zero-complete-authoritative-resource`,
  `zero-complete-effective-config-digest`, `one-complete-reference-digest`,
  `one-complete-authoritative-resource`,
  `one-complete-effective-config-digest`, `one-failed-reference-digest`,
  `one-failed-authoritative-resource`, `one-failed-effective-config-digest`,
  `two-complete-reference-digest`, `two-complete-authoritative-resource`,
  `two-complete-effective-config-digest`, `two-failed-reference-digest`,
  `two-failed-authoritative-resource` y
  `two-failed-effective-config-digest`. Completa la matriz con
  `zero-failed-effective-config-digest`, `zero-reference-id`,
  `zero-reference-kind`, `zero-reference-path`, `zero-reference-size`,
  `zero-effective-config-id`, `zero-effective-config-kind`,
  `zero-effective-config-path`, `zero-effective-config-size`,
  `zero-event-manifest-id`, `zero-event-manifest-kind`,
  `zero-event-manifest-path`, `zero-event-manifest-size` y
  `zero-complete-event-digest`.
- Cada mutación conserva forma Draft 2020-12 válida y es rechazada por Rust y
  por el oráculo semántico independiente. El catálogo constante del propio test
  se compara con los nombres realmente recorridos, de modo que el expediente no
  puede acreditar casos omitidos por una lista meramente documental.
- H-F16-R14 añade
  `execution::tests::manifest_graph_rejects_every_dynamic_edge_contradiction`.
  Sus ocho positivos son `graph-positive-zero-complete`,
  `graph-positive-zero-failed`, `graph-positive-one-complete`,
  `graph-positive-one-failed`, `graph-positive-two-complete`,
  `graph-positive-two-failed`, `graph-positive-multiple-complete` y
  `graph-positive-multiple-failed`: cubren cierres completo/fallido, cero/una/dos
  reducciones y una ejecución con tres configuraciones y dos eventos.
- Las 24 contradicciones iniciales de R14, estructuralmente válidas, son
  `graph-configuration-copy-digest`, `graph-configuration-copy-id`,
  `graph-configuration-copy-class`, `graph-configuration-copy-path`,
  `graph-configuration-copy-size`, `graph-configuration-copy-reused`,
  `graph-configuration-input-duplicate`,
  `graph-configuration-kind-incompatible`,
  `graph-general-configuration-missing`, `graph-event-unit-crossed`,
  `graph-event-unit-orphan`, `graph-event-unit-missing`,
  `graph-execution-unit-crossed`, `graph-output-artifact-missing`,
  `graph-output-wrong-category`, `graph-artifact-path-duplicate`,
  `graph-orphan-event-artifact`, `graph-failed-configuration-copy-digest`,
  `graph-failed-execution-unit-crossed`,
  `graph-one-complete-output-missing`, `graph-one-failed-output-missing`,
  `graph-two-complete-output-missing`, `graph-two-failed-output-missing` y
  `graph-two-complete-artifact-path-duplicate`. El catálogo constante se compara
  al final; R25 añade las cinco aristas de `resolution-evidence`, por lo que el
  catálogo actual contiene 37 nombres efectivamente recorridos.
- H-F16-R15 añade
  `execution::tests::manifest_dynamic_graph_rejects_diagnostic_and_reference_inconsistencies`.
  Sus siete positivos reales son `dynamic-positive-zero-complete`,
  `dynamic-positive-zero-failed`, `dynamic-positive-one-complete`,
  `dynamic-positive-one-failed`, `dynamic-positive-two-complete`,
  `dynamic-positive-two-failed` y `dynamic-positive-percent-complete`; cubren
  completo/fallido con cero, una y dos reducciones, diagnóstico compartido por
  ejecución y reducción fallidas, y el cálculo porcentual.
- Sus 27 contradicciones individuales, aceptadas estructuralmente por Python y
  rechazadas por Rust y por el oráculo declarativo independiente, son
  `dynamic-diagnostic-execution-absent`,
  `dynamic-diagnostic-reduction-absent`, `dynamic-diagnostic-invented`,
  `dynamic-diagnostic-category`, `dynamic-diagnostic-meaning`,
  `dynamic-diagnostic-required-context`, `dynamic-diagnostic-duplicate`,
  `dynamic-diagnostic-orphan`, `dynamic-reference-sequence-duplicate`,
  `dynamic-reference-total`, `dynamic-reference-unknown-count`,
  `dynamic-reduction-partition-overlap`,
  `dynamic-reduction-partition-missing`,
  `dynamic-reduction-selection-order`,
  `dynamic-reduction-sequence-missing`, `dynamic-reduction-out-of-range`,
  `dynamic-reduction-overlap`, `dynamic-reduction-global-sum`,
  `dynamic-reduction-per-sequence-sum`, `dynamic-reduction-eligible`,
  `dynamic-reduction-declared-target`, `dynamic-reduction-derived-target`,
  `dynamic-reduction-region`, `dynamic-reduction-residual`,
  `dynamic-reduction-derived-class`, `dynamic-reduction-derived-ploidy` y
  `dynamic-reduction-percent-target`. El propio test compara este conjunto con
  `MANIFEST_DYNAMIC_CORPUS_CASES`; el expediente verifica además que cada nombre
  está en fuente, listado ejecutado y documentación.
- H-F16-R16 añade
  `execution::tests::manifest_bundle_closes_exact_percent_seed_header_normalization_and_event_edges`.
  Sus seis positivos son
  `exhaustive-positive-percent-above-integer-boundary`,
  `exhaustive-positive-percent-exponent-notation`,
  `exhaustive-positive-percent-below-integer-boundary`,
  `exhaustive-positive-non-utf8-header`,
  `exhaustive-positive-multiple-events-normalizations` y
  `exhaustive-positive-two-reductions`. Los tres primeros son
  ejecuciones reales, no manifiestos retocados: cubren S=20 a ambos lados de la
  frontera de 30 %, incluida p=`30.0000000000000000001` con target=13, y
  notación exponencial. El cuarto conserva descripción no UTF-8 reversible; el
  quinto aporta varias configuraciones, dos eventos y tres normalizaciones; el
  sexto acredita dos reducciones en el mismo expediente.
- Sus 20 contradicciones por campo/arista son
  `exhaustive-reduction-root-derived`,
  `exhaustive-reduction-subseed-derived`,
  `exhaustive-derived-reference-order`, `exhaustive-header-token`,
  `exhaustive-header-description`, `exhaustive-header-description-non-utf8`,
  `exhaustive-normalization-sum-before`,
  `exhaustive-normalization-delta-before`, `exhaustive-normalization-before`,
  `exhaustive-normalization-after`, `exhaustive-normalization-sum-after`,
  `exhaustive-normalization-delta-after`,
  `exhaustive-normalization-decision`, `exhaustive-normalization-order`,
  `exhaustive-normalization-duplicate-value`,
  `exhaustive-normalization-duplicate-parameter`,
  `exhaustive-event-general-subseed`, `exhaustive-event-effective-subseed`,
  `exhaustive-event-manifest-subseed` y
  `exhaustive-event-correlated-derived-subseed`. Las mutaciones correlacionadas
  alteran también la configuración efectiva o el manifiesto de evento para que
  una mera desigualdad entre copias no sustituya al oráculo semántico.
- H-F16-R17 añade
  `execution::tests::manifest_bundle_closes_toml_fasta_and_effective_value_edges`.
  Sus seis positivos reales son `r17-positive-percent-separators`,
  `r17-positive-percent-plus-exponent`,
  `r17-positive-percent-high-precision`, `r17-positive-fasta-non-utf8`,
  `r17-positive-effective-multiple-scopes` y
  `r17-positive-two-reductions`. Cubren separadores TOML válidos, signo y
  exponente, precisión superior a `f64`, cabecera no UTF-8, ámbitos global,
  individual y de secuencia y varias reducciones.
- Las 34 contradicciones R17, estructuralmente válidas salvo la mutación física,
  son `r17-percent-leading-separator`, `r17-percent-trailing-separator`,
  `r17-percent-double-separator`, `r17-percent-separator-before-dot`,
  `r17-percent-separator-after-dot`, `r17-percent-missing-integer`,
  `r17-percent-missing-fraction`, `r17-percent-leading-zero-integer`,
  `r17-percent-leading-zero-float`, `r17-percent-signed-leading-zero`,
  `r17-percent-empty-exponent`, `r17-percent-empty-signed-exponent`,
  `r17-percent-leading-exponent-separator`,
  `r17-percent-trailing-exponent-separator`,
  `r17-fasta-correlated-id-header`, `r17-fasta-correlated-description`,
  `r17-fasta-correlated-length`, `r17-fasta-correlated-unknown-count`,
  `r17-fasta-correlated-order`, `r17-fasta-intermediate-symlink`,
  `r17-effective-duplicate-identical`,
  `r17-effective-duplicate-contradictory`, `r17-effective-reordered`,
  `r17-effective-parameter-cross`, `r17-effective-value-cross`,
  `r17-effective-winner-cross`, `r17-effective-selectors-cross`,
  `r17-effective-provenance-path`, `r17-effective-provenance-location`,
  `r17-effective-automatic-missing`, `r17-effective-correlated-ploidy`,
  `r17-effective-source-path`, `r17-effective-controlled-duplicate` y
  `r17-effective-controlled-count`. Las mutaciones FASTA y ploidía alteran de
  forma correlacionada manifiesto y configuración efectiva, pero dejan intactos
  FASTA o TOML; por ello no se rechazan mediante una mera desigualdad JSON.
- H-F16-R18 añade `r18-effective-dominated-base-copy` y
  `r18-effective-foreign-object`. Ambas mutaciones reescriben el documento
  enlazado y recalculan tamaño y SHA-256 en cada aparición del artefacto, por lo
  que alcanzan la semántica: la primera conserva el contexto individual pero
  sustituye la ganadora dominante `variant-ind-1` por `variant-base`; la segunda
  usa `object_id=zzz-foreign`. Producción y oráculo independiente rechazan las
  dos después de que Python/jsonschema acepte su estructura.
- H-F16-R19 y H-F16-R21 añaden
  `r19-biological-parameter-technical-context` y
  `r21-orphan-position-with-individual`. La primera reemplaza un valor
  biológico individual por `variant-base` aplicado a `event-1` con contexto
  técnico; la segunda conserva un individuo válido pero añade una posición sin
  referencia ni secuencia. En ambas se reescribe la configuración efectiva y
  se recalculan tamaño y SHA-256 en todos sus enlaces. Python/jsonschema acepta
  la forma y producción y oráculo independiente rechazan la semántica.
- El mismo paquete contiene como testigo positivo `coverage-base` resuelto para
  `event-1` en su eje técnico permitido. La batería global conserva además
  positivos completo/fallido con cero, una y dos reducciones, de modo que los
  cruces nuevos no rompen ninguno de esos cierres.
- `effective_config` ejecuta 15 casos. Entre los cinco previos,
  `request_context_drives_combined_dominance_independently_of_rule_order`
  recorre ámbitos individual+secuencia y dos órdenes de reglas, y
  `foreign_object_and_incomparable_maxima_are_rejected_before_snapshot`
  rechaza el objeto ajeno y dos máximas incomparables antes de crear la
  instantánea. Los cuatro casos R19–R21 comprueban: parámetro biológico en eje
  técnico prohibido; dominancia base→población→subpoblación;
  base→fase→generación; secuencia→región→locus con posición; y rechazo de
  posición huérfana, locus sin posición o coordenada fuera de rango.
- H-F16-R22 añade la matriz
  `every_parameter_is_checked_at_sequence_region_and_locus_resolution_levels`:
  cruza los ocho parámetros con Sequence/Region/Locus y compara literalmente
  los niveles admitidos. `spatial_level_and_primary_object_follow_parameter_metadata`
  conserva además los negativos `coverage+locus`, `recombination+region` y
  referencia usada indebidamente como objeto espacial principal.
- H-F16-R23 añade
  `hierarchical_contexts_are_not_inferred_from_flat_identifier_lists` y
  `hierarchical_contexts_must_match_parent_interval_and_coordinate_evidence`.
  Las listas planas no prueban relaciones; un `ResolutionDomain` sintético
  explícito demuestra los tres recorridos positivos de dominancia y los
  negativos de padre, intervalo o coordenada discrepantes.
- El corpus correlacionado R22 ejecuta `r22-coverage-locus`,
  `r22-recombination-region`, `r22-ploidy-locus` y `r22-reference-object`.
  Reescribe la configuración efectiva, recalcula tamaño/SHA-256 en cada enlace
  y llega a la semántica tras ser aceptado por Python/jsonschema.
- El paquete con dos escenarios contiene un valor efectivo cuya fase
  `phase-1` solo existe en la copia del escenario y cuya generación 5 domina a
  la fase mediante un dominio sintético versionado. Producción y oráculo
  aceptan los positivos completo/fallido. Los casos correlacionados
  `r23-complete-missing-phase-domain-proof` y
  `r23-failed-missing-phase-domain-proof` eliminan la prueba del intervalo; los
  casos `r24-complete-scenario-copy-id` y `r24-failed-scenario-copy-id` alteran
  el ID dentro del TOML copiado y recalculan su identidad física. Ambos
  validadores semánticos rechazan los cuatro.
- H-F16-R25 publica `inputs/resolution/evidence.json` con clase
  `resolution_evidence`, esquema y orden versionados, dominio y solicitudes
  admitidos antes de derivar `effective-config.json`. Producción lo relee con
  `O_NOFOLLOW`, valida tamaño/SHA-256 y reconstruye el resultado; el oráculo lo
  lee y resuelve por separado. Los casos correlacionados
  `r25-complete-correlated-domain-provider-interval`,
  `r25-failed-correlated-domain-provider-interval` y
  `r25-correlated-global-output-to-event` alteran solo la salida derivada y
  recalculan sus enlaces, por lo que la evidencia independiente los rechaza.
  `execution::tests::resolution_evidence_toctou_is_rejected_before_publication`
  sustituye físicamente sus bytes antes de la revalidación final y demuestra
  rechazo previo al publicador, destino ausente e `INCOMPLETE` conservado.
  `graph-resolution-evidence-digest`, `graph-resolution-evidence-id`,
  `graph-resolution-evidence-kind`, `graph-resolution-evidence-path` y
  `graph-resolution-evidence-size` cubren por separado el enlace del manifiesto
  contra el catálogo.
- H-F16-R26 añade
  `effective_config::phase_ids_are_globally_unique_across_scenarios_independently_of_order`:
  una misma fase en dos escenarios se rechaza en ambos órdenes y el diagnóstico
  conserva ID y ruta de las dos fuentes. Los casos correlacionados
  `r26-complete-duplicate-phase-across-scenarios` y
  `r26-failed-duplicate-phase-across-scenarios` reescriben la copia TOML,
  recalculan tamaño/SHA-256 y son rechazados por producción y oráculo.
- H-F16-R27 añade
  `effective_config::parameter_order_v1_is_literal_and_independent_of_rule_and_request_permutations`:
  la tabla literal v1 ordena `variant_rate` antes que `output_enabled` para el
  mismo individuo bajo permutaciones de reglas y solicitudes. El productor no
  usa discriminantes ni texto y el oráculo usa una tabla separada. Los casos
  `r27-complete-parameter-order` y `r27-failed-parameter-order` intercambian
  ambos valores después de recalcular la identidad del documento enlazado y son
  rechazados semánticamente.
- El recorrido de esquema contextual usa Python/jsonschema como oráculo externo:
  acepta positivos completo/fallido y rechaza seis mezclas contextuales no
  permitidas.
- `f1_6_traceability`: 4 casos; el catálogo explícito enlaza cada uno de los 22
  requisitos con funciones presentes en su fuente, citadas literalmente en la
  matriz y conservadas en `executed-test-list.txt`; además comprueba el balance
  automático 83 y los enlaces DTPD-101/DTPD-116→DTC-056/DTC-057.
- `results_publication`: 14 recorridos públicos; incluye capacidad opaca, puerta
  Draft 2020-12, cruce exacto del catálogo y los ficheros del workspace, rechazo
  de workspace vacío/ajeno aunque sea idéntico, carrera no-clobber, rechazo y
  SIGKILL.
- Módulos unitarios internos: 16 casos de publicación y ejecución. La matriz de
  nueve fallos recorre una/dos reducciones, índices 0/1 y FASTA/mapa/metadatos;
  los corpus ejecutan positivos completo/fallido, mezclas contextuales y
  contradicciones entre colecciones. Publicación cubre extra, enlace simbólico
  de fichero, contenido/TOCTOU, enlace físico externo, padre simbólico real y
  canonicalización directa, además del positivo que limita `st_nlink == 1` a
  ficheros regulares.

## Oráculos independientes

- Python 3.12.3 con jsonschema 4.10.3 acepta configuración efectiva, manifiesto
  y mapa válidos; el contrato del mapa es Draft 2020-12 cerrado.
- Un reconstructor de prueba que no llama al productor lee el FASTA fuente y el
  JSON, usa intervalos origen y derivados, clases, ambos órdenes, orientación y
  parámetros efectivos, y emite los bytes del FASTA canónico completo.
- Un oráculo semántico separado deserializa de nuevo los bytes del manifiesto y
  reconstruye mapas independientes de artefactos, referencias, reducciones y
  unidades. Comprueba propiedad única de FASTA, orígenes, estados, recurso
  autoritativo, configuración efectiva, eventos e identidades cruzadas sin
  llamar a ninguna función del validador de producción ni usar Python como
  sustituto de semántica dinámica.
- El oráculo de grafo parte de raíces y aristas, no del orden del validador Rust:
  construye por separado catálogos ID/ruta, propietarios y alcanzabilidad;
  comprueba la única configuración general, escenarios canónicos, unidad raíz,
  biyección evento↔unidad, referencias/reducciones y la partición exacta de
  salidas. Rechaza también todo artefacto de clase desconocida o no alcanzable.
- H-F16-R15 añade la tabla declarativa `INDEPENDENT_DYNAMIC_INVARIANTS`. Sus dos
  reglas construyen contratos diagnostics-v2 y catálogos de
  referencia/reducción desde los bytes, sin llamar a producción. Python queda
  limitado a demostrar que cada mutación conserva la forma Draft 2020-12.
- H-F16-R17 añade tres recorridos independientes: una máquina de estados para la
  gramática decimal TOML, un parser FASTA directo sobre bytes y un catálogo de
  reglas/colecciones leído desde la copia TOML sin llamar al resolvedor de
  producción. Python continúa limitado a estructura JSON.
- H-F16-R18 conserva `resolution_context` en cada valor efectivo. Producción
  reutiliza el `RuleSet` canónico y el catálogo de identidades declarado; el
  oráculo parsea por separado reglas e identidades desde TOML, calcula
  aplicabilidad y rangos de especificidad y obtiene el máximo sin llamar al
  resolvedor de producción. La excepción automática queda limitada a la regla
  global realmente sustituida por una solicitud explícita del mismo parámetro.
- H-F16-R19–R21 añaden al oráculo una tabla local parámetro→ejes permitidos,
  separada del catálogo de producción. También reconstruyen jerarquías de
  contexto y longitudes de secuencia desde el TOML: cada hijo exige su padre,
  referencia y secuencia son inseparables y posición/locus requieren contexto
  espacial y límites válidos. En F1 las entidades población/subpoblación no
  contienen todavía una arista genealógica explícita; por tanto la coherencia
  comprobable aquí es que ambos niveles estén presentes y declarados, sin
  inventar una relación no representada por el formato. H-F16-R23 sustituye
  esa adaptación provisional: desde esta ronda la jerarquía solo se acepta con
  prueba exacta en `ResolutionDomain`; la mera presencia ya no basta.
- H-F16-R22 añade una tabla independiente parámetro→niveles espaciales y deriva
  el nivel más específico directamente del contexto. H-F16-R23 reconstruye el
  snapshot del dominio, valida identidades, intervalos, loci, orden y proveedor
  y lo usa para cada relación. H-F16-R24 relee todas las copias de escenario
  catalogadas, exige ruta e ID canónicos y agrega sus fases/eventos sin llamar a
  producción.

### Inventario exhaustivo de repeticiones dinámicas dentro del manifiesto

| Repetición/arista | Puerta semántica | Oráculo independiente |
|---|---|---|
| `configuration_inputs[*].copy` ↔ `artifacts` | Identidad completa, clase/ID/ruta por `kind`, una general, escenarios y propietario únicos | Catálogo y raíces de configuración reconstruidos |
| `references[*].fasta` y `authoritative_resource` ↔ `artifacts` | Identidad completa, ruta autoritativa y propietario FASTA | Mapa de referencias/FASTA independiente |
| `references[*].sequences[*]` y `total_genome_size` ↔ FASTA autoritativo | IDs únicos, suma sin desbordamiento y relectura por dirfd/`O_NOFOLLOW`; orden, cabecera/descripción reversible, ID, longitud y N deben coincidir con los bytes | Parser FASTA independiente directo, incluido no UTF-8 |
| `effective_config` ↔ `artifacts` | ID/clase canónicos e identidad completa | Raíz opcional y alcanzabilidad |
| `events[*].manifest` ↔ `artifacts`; eventos ↔ unidades | Identidad completa, propietario único y biyección 1:1 | Conjuntos de eventos, manifiestos y unidades reconstruidos |
| `execution_id`/`state` ↔ unidad `execution` | Una unidad y coincidencia exacta | Raíz de ejecución independiente |
| `reference_reductions[*]`, referencias, unidades y artefactos | Fuente, derivada, mapa, estados, origen y propiedad | Recorrido independiente conservado de R11–R13 |
| `included_sequence_ids`/`omitted_sequence_ids` | Conjuntos únicos, disjuntos y partición exacta de la fuente; `included` coincide con selección y catálogo derivado | Álgebra de conjuntos reconstruida |
| `selection_order[*]`, región, residual y objetivo | Orden 1..N, secuencia/rango válidos, no solape, longitudes full/residual y sumas globales/por secuencia exactas; reutiliza `preflight_reduced_reference` | Intervalos agrupados/ordenados y ecuaciones recalculadas |
| `target.mode/declared_value/derived_bases` | Entero exacto o `ExactDecimal`; gramática TOML validada antes de retirar `_`; coincide con objetivo, suma y tamaño derivado | Máquina de estados decimal y aritmética larga independientes |
| `outputs.minimum/default_enabled/explicitly_requested` ↔ `artifacts` | Partición exacta; `manifest-general` es la única salida raíz no catalogada | Conjuntos esperados recalculados por clase |
| `artifacts[*].id/path/kind` | IDs y rutas globalmente únicos; clases v1 cerradas; todo artefacto tiene propietario | Catálogo ID/ruta y conjunto de alcanzabilidad |
| `units[*].state/diagnostic_code` ↔ `diagnostics[*]` | Toda unidad fallida resuelve código; completas usan `null`; compartir se permite; códigos extra/huérfanos no | Conjuntos de códigos embebidos/referidos iguales |
| `diagnostics[*].code/category/message/context` ↔ `diagnostics-v2.tsv` | Código único y existente, prefijo/categoría, significado y todos los contextos requeridos no vacíos | Catálogo TSV parseado directamente de forma independiente |
| `randomness.root_seed` y subsemilla de cada reducción | Raíz idéntica a la configuración efectiva; derivación DTC-037 con ámbito `reference_reduction` y catálogo v3 | Material canónico y BLAKE3 reconstruidos localmente sin llamar a `derive_subseed` |
| Orden fuente→derivada y `header`/`description` | Orden derivado igual al filtrado del fuente; token FASTA igual a `sequence.id`; descripción reconstruida byte a byte, incluido `bytes-hex` | Parser reversible independiente sobre bytes |
| `normalizations[*]` y `values[*]` | Parámetros e identidades únicos y ordenados; Neumaier, sumas, deltas, `before`/`after`, tolerancia y decisión exactos; igualdad efectiva↔general | Implementación Neumaier separada y comparación por bits |
| `event_ids`, `event_subseeds`, `events[*].subseed` y manifiesto de evento | Derivación v3 `technical_event`, identidad/orden y cruce configuración efectiva↔general↔evento; referencia y unidad exactas | Material BLAKE3 separado y relectura de los tres documentos |
| Campos exclusivos de configuración efectiva | `source_path`, ruta, IDs, aleatoriedad, normalizaciones y conjuntos controlados se reconstruyen desde las copias TOML; `effective_values` conserva el contexto, exige objeto declarado de clase compatible, identidad/orden únicos, dominio y máxima dominante con regla/valor/selectores/procedencia exactos y ploidía cruzada | Catálogo independiente de reglas e identidades leído desde TOML; aplicabilidad y dominancia recalculadas por rangos separados |

El resto de repeticiones de nivel superior también queda clasificado: Draft
2020-12 fija cardinalidad y pertenencia exactas de `dependencies`, `checks` y
`reproducibility`; obliga a vaciar `tools`, `profiles` y
`outputs.default_enabled`, y fija literalmente `outputs.disabled`.
`configuration_inputs` se cruza con `scenario_ids`; `event_ids` se cruza con
subsemillas y eventos. Aunque `individual_ids`, `controlled_frequency_sets` y
`effective_values` solo existen en la configuración efectiva, R17 los deduce de
las copias TOML autoritativas: ya no se confía únicamente en el enlace de
tamaño/SHA. Las solicitudes explícitas pueden sustituir el valor base global
del mismo parámetro; toda entrada conservada debe corresponder exactamente a su
regla declarada y los valores automáticos no sustituidos son obligatorios.
`normalizations[*].values` no contiene IDs de artefacto, pero sus identidades,
orden y ecuaciones internas se validan y la colección completa se cruza con la
configuración efectiva; `specific_tolerances` está vacío. `references[*].indexes`
también está vacío en v1. Los eventos generales solo admiten unidades completas
con `diagnostic_code=null`; la misma política genérica cubre cualquier unidad
fallida que el esquema llegue a admitir. Estados y
`validated` se comprueban por el esquema y, cuando forman parte de una arista,
por identidad íntegra en el grafo. Las repeticiones internas de la configuración
efectiva y de los manifiestos de evento pertenecen a sus propios documentos y a
su relectura; no se presentan como deducibles de un único `manifest.json`.

Los escalares restantes también quedan delimitados: versión del programa,
rutas originales de configuración, alfabeto y formato canónico están cerrados
por esquema o pertenecen al dominio del fichero físico. `route`, la referencia
fuente, escenarios, reducciones y normalizaciones se cruzan con la configuración
efectiva. Las semillas sí son derivables después de releer esa configuración: la
puerta de producción usa el catálogo v3 canónico y el oráculo reconstruye por
separado material, longitudes, BLAKE3 y codificación hexadecimal.

Las igualdades que requieren leer configuración efectiva y manifiestos de evento
pertenecen a la relectura semántica del paquete. La identidad física de todos los
artefactos se vuelve a cerrar después mediante el recorrido no-follow y los
validadores externos; no se confunden ambos oráculos.
- El validador relee los bytes del mapa escritos y rechaza propiedades
  desconocidas en raíz/segmento/intervalo, JSON truncado, etiquetas
  full/residual intercambiadas, hueco/solape, orden alterado y coordenadas fuera
  de rango; también comprueba orientación, cobertura, bases y orden/clase
  directos del catálogo derivado mediante mutaciones negativas.
- Los recorridos fallidos posteriores a metadatos releen configuración efectiva
  y manifiesto, cruzan cada ID/estado/validación y validan ambos esquemas. La
  segunda reducción puede fallar sin degradar la primera ya completa.
- La publicación solo admite `manifest_general` por la vía estructural y
  semántica dedicada, vuelve a leer `manifest.json` justo antes del cambio
  atómico y rechaza una mutación posterior a su primera validación. Los hooks y
  publicadores inyectables son internos; la API pública solo ofrece
  `ExecutionWorkspace::publish()`.
- La capacidad solo nace desde el workspace, conserva token privado y snapshot
  del catálogo y no autoriza otra instancia idéntica. Antes de escribir y de
  nuevo tras el hook pre-rename, el árbol se recorre desde dirfd con
  `openat`/`statat`/`fstat` y `O_NOFOLLOW`: no admite symlinks, especiales,
  extras, ausentes, directorios ajenos, rutas no canónicas ni dos rutas para el
  mismo dispositivo/inodo; cada fichero se relee desde su fd y contrasta
  tamaño/SHA-256. La última comprobación precede inmediatamente al rename, sin
  prometer imposibilidad absoluta frente a una carrera externa posterior.
- `manifest-general-runtime-validation-v1` compila con `jsonschema` 0.49.9
  exactamente el esquema raíz y el recurso común Draft 2020-12 publicados y
  embebidos. El corpus diferencial restaurado ejecuta Python/jsonschema y Rust
  sobre ocho casos, incluidos los positivos completo/fallido y seis mezclas
  contextuales; otro corpus comprueba diez contradicciones entre colecciones que
  el esquema admite y que tanto la semántica Rust como el oráculo semántico
  independiente rechazan. La escritura recibe una capacidad opaca;
  la publicación contrasta catálogo, bytes, tamaño y SHA-256 antes de escribir
  y otra vez inmediatamente antes del cambio atómico.
- GNU `sha256sum` 9.4 coincide con el digest del manifiesto.
- `samtools faidx` 1.21 indexa una copia de la derivada fuera del resultado.
- `scripts/check-independent-validators.sh
  /tmp/gps-f16-validators-r15-20260815` volvió a superar las
  categorías FASTA, SHA-256, tar y JSON. El intento OCI dentro del aislamiento
  terminó en código 1 por falta de permisos sobre `/var/run/docker.sock`; no se
  contabilizó entonces como puerta superada ni como ausencia del daemon.
- Una ejecución actual en
  `/tmp/gps-f16-independent-r15-20260815` superó por separado Python/jsonschema
  para manifiesto, configuración efectiva y mapa, `samtools faidx` para una
  copia externa del FASTA reducido y `sha256sum` para los cinco artefactos del
  catálogo; tamaños y digest coinciden exactamente con `manifest.json`.
- La ejecución R16 en
  `/tmp/gps-f16-independent-r16-HtlQprzX/result` superó Python/jsonschema para
  manifiesto general, configuración efectiva y manifiesto de evento; Python
  `hashlib` contrastó tamaño/SHA de sus cinco artefactos, GNU `sha256sum`
  calculó los mismos cinco ficheros y `samtools faidx` indexó el FASTA (el
  índice se conservó fuera del resultado).
- `scripts/check-independent-validators.sh
  /tmp/gps-f16-validators-r16-oci-20260815-final`, ejecutado fuera del
  aislamiento del socket, superó las cinco categorías: FASTA, SHA-256, tar,
  JSON Schema y OCI. La prueba OCI obtuvo el rechazo esperado al inspeccionar
  una imagen inexistente fijada por digest (`status=1`, «no such image») y no un
  error de infraestructura. Servidor Docker 29.6.1, Linux/amd64; cliente Docker
  29.6.1, build `8900f1d`.
- La ejecución R17 en `/tmp/gps-f16-independent-r17-20260815` superó
  Python/jsonschema para manifiesto general, configuración efectiva y
  manifiesto de evento. Python `hashlib` contrastó tamaño y SHA-256 de los cinco
  artefactos, GNU `sha256sum` produjo los mismos cinco digests y `samtools
  faidx` indexó una copia externa del FASTA.
- `scripts/check-independent-validators.sh
  /tmp/gps-f16-validators-r17-20260815-final`, ejecutado fuera del aislamiento
  del socket, superó FASTA, SHA-256, tar, JSON Schema y OCI. La prueba OCI
  obtuvo el rechazo esperado de la imagen inexistente (`status=1`, «no such
  image»), con servidor Docker 29.6.1 Linux/amd64 y cliente Docker 29.6.1 build
  `8900f1d`.
- La ejecución R18 en `/tmp/gps-f16-independent-r18-20260815` superó
  Python/jsonschema para manifiesto general, configuración efectiva con
  `resolution_context` y manifiesto de evento. Python `hashlib` contrastó
  tamaño/SHA-256 de sus cinco artefactos, GNU `sha256sum` produjo los mismos
  cinco digests y `samtools faidx` indexó una copia externa del FASTA.
- `scripts/check-independent-validators.sh
  /tmp/gps-f16-validators-r18-20260815-final`, ejecutado fuera del aislamiento
  del socket, superó FASTA, SHA-256, tar, JSON Schema y OCI. La prueba OCI
  obtuvo el rechazo esperado de la imagen inexistente (`status=1`) con Docker
  29.6.1 build `8900f1d`, no un fallo de infraestructura.
- La ejecución R19–R21 en
  `/tmp/gps-f16-independent-r19-r21-20260815` superó Python/jsonschema para
  manifiesto general, configuración efectiva y manifiesto de evento. Python
  `hashlib` contrastó tamaño/SHA-256 de los cuatro artefactos catalogados; GNU
  `sha256sum` produjo los digests de los cinco ficheros del paquete, incluido
  el manifiesto general, y `samtools faidx` indexó una copia externa del FASTA.
- `scripts/check-independent-validators.sh
  /tmp/gps-f16-validators-r19-r21-20260815-final`, ejecutado fuera del
  aislamiento del socket, superó FASTA, SHA-256, tar, JSON Schema y OCI. La
  prueba OCI obtuvo el rechazo semántico esperado de la imagen inexistente
  (`status=1`), con Docker 29.6.1 build `8900f1d`.
- La ejecución R22–R24 en
  `/tmp/gps-f16-independent-r22-r24-20260815` superó Python/jsonschema para
  manifiesto general, configuración efectiva con `resolution_domain` y
  manifiesto de evento. Python `hashlib` contrastó tamaño/SHA-256 de los cuatro
  artefactos catalogados; GNU `sha256sum` produjo los digests de los cinco
  ficheros del paquete y `samtools faidx` indexó una copia externa del FASTA.
- `scripts/check-independent-validators.sh
  /tmp/gps-f16-validators-r22-r24-20260815-final`, ejecutado fuera del
  aislamiento del socket, superó FASTA, SHA-256, tar, JSON Schema y OCI. La
  prueba OCI obtuvo el rechazo semántico esperado de la imagen inexistente
  (`status=1`), con Docker Engine 29.6.1, commit `8ec5ab3`.
- La ejecución R25–R27 en
  `/tmp/gps-f16-independent-r25-r27-20260815` superó Python/jsonschema 4.10.3
  para manifiesto general, configuración efectiva, manifiesto de evento y la
  evidencia de resolución v1. Python `hashlib` contrastó tamaño/SHA-256 de los
  cinco artefactos catalogados; GNU `sha256sum` produjo los seis digests del
  paquete, incluido el manifiesto general, y `samtools faidx` 1.21 indexó y
  consultó una copia externa del FASTA.
- `scripts/check-independent-validators.sh
  /tmp/gps-f16-validators-r25-r27-20260815-final`, ejecutado fuera del
  aislamiento, superó FASTA, SHA-256, tar, JSON Schema y OCI. La categoría OCI
  obtuvo el rechazo semántico esperado de la imagen inexistente (`status=1`),
  no un fallo de infraestructura. Cliente Docker 29.6.1 build `8900f1d` y
  servidor Docker Engine 29.6.1 commit `8ec5ab3`, Linux/amd64.

## Puertas completas

- `cargo fmt --all -- --check`: código 0.
- `cargo clippy --workspace --all-targets --all-features --locked -- -D warnings`:
  código 0.
- `cargo test --workspace --all-targets --all-features --locked`: 196 pruebas
  enumeradas, todas superadas. El listado F1.6 relevante y el total se conservan
  en `executed-test-list.txt`.
- `cargo test --workspace --doc --all-features --locked`: 2 `compile_fail`
  superadas, incluida la imposibilidad de importar `AtomicPublisher`.
- `cargo doc --workspace --no-deps --all-features --locked`: código 0.
- `cargo build --release --workspace --all-features --locked`: código 0.
  Binario producido: 9.963.992 bytes; no existía una línea base equivalente
  previa a la dependencia Draft 2020-12 y no se afirma una mejora de tamaño.
- `cargo metadata --locked --format-version 1 --no-deps` y `cargo tree --locked -e
  features`: código 0; `jsonschema` 0.49.9 queda fijado sin sus *features*
  predeterminadas y detrás del adaptador privado.
- `cargo audit`: 149 dependencias examinadas contra 1.216 avisos RustSec,
  código 0.
- `git diff --check`: código 0.

## Resultado de revisión

Los cinco hallazgos de revisión propia, los cuatro graves de la primera revisión
independiente, los graves de seguimiento y M1/M2 se corrigieron y tienen
regresión en `review.md`. H-F16-R25–R27 quedan 3/3 aplicados con adaptación y
H-F16-R05–R27 están aplicados; no quedan defectos bloqueantes o graves
conocidos. F14, MAN-90, SAL-27, 22/22 y 83/83 se restauran para una
nueva revisión independiente, sin aprobar automáticamente la pausa. La
revisión de rendimiento conserva
explícito el peor caso del recorrido de huecos sin prometer una escala no
medida. La nueva puerta física recorre y resume cada artefacto antes del
manifiesto y la publicación; no existía una línea base específica previa y no
se reconstruye una comparación retrospectiva. Es el coste aceptado de cerrar la
integridad y la ventana controlada. No se ha iniciado F2 ni se ha aprobado
automáticamente la pausa.
