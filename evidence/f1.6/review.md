# Revisión vigente de F1.6

La evidencia histórica con símbolos retirados se sustituye por este registro
revalidable. La autoridad de requisitos y de símbolos es la matriz F1.6 y la
lista ejecutada v2.

| Hallazgo | Gravedad | Estado | Regresión actual |
|---|---|---|---|
| El selector previo permitía bifurcar el contrato por algoritmo configurable. | Grave | Aplicado | `reference_reduction_config::selection_algorithm_is_rejected_at_the_toml_boundary`, `effective_config_schema::effective_reduction_rejects_removed_selection_algorithm` y `f1_6_acceptance::reduction_publishes_a_reproducible_map_without_algorithm_metadata`. |
| La disposición podía recorrer una vez por cada base de holgura. | Grave | Aplicado con adaptación | `reference_reduction::tests::direct_layout_handles_hundred_megabases_with_ten_blocks_without_per_base_work`; el benchmark compara 100 MB con 10, 100 y 1.000 bloques. |
| El accessor denominado `segments_in_selection_order` exponía el orden de inserción por secuencia. | Grave | Aplicado | `reference_reduction::tests::canonical_selector_reproduces_the_frozen_direct_vector` exige el orden 1, 2, residual 3. |
| La subsemilla activa no tenía fila de catálogo coherente con el selector canónico. | Grave | Aplicado | `subseed_derivation::published_catalogs_are_immutable_and_v6_has_an_independent_oracle` y `f1_runtime_scopes_remain_on_derivation_v1_and_match_active_catalog_v6`. |
| La lista F1.6 incluía símbolos ya retirados. | Grave | Aplicado | `f1_6_traceability::simplified_f1_expedient_closes_exactly_54_requirements_with_real_symbols` compara exactamente sus símbolos vivos con la matriz vigente. |
| Una prueba de fixture histórica permanecía ignorada. | Menor | Eliminado | Se retiró `archived_e14705a_reduction_tree_revalidates_physically_and_resolves_its_derived_reference`; el fixture que aún se usa cubre solo el rechazo de manifiesto v1, no el selector retirado. |

Puertas requeridas: formato, Clippy estricto, pruebas de todos los targets y
features, documentación y `git diff --check`.
