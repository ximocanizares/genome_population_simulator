# Inventario de conservación para F1 0.2

- **Base:** `66ad8b2`.
- **Decisión:** DTC-078.
- **Auditoría:** `design/auditoria-sistematica-cierre-f1-0.1.md`.
- **Estado:** ejecutado hasta PAUSA-F1.4.

## Componentes

| Ruta o responsabilidad | Clase | Destino en F1 0.2 |
| --- | --- | --- |
| `src/lib.rs` y `src/domain.rs` | A/D | Conservar la frontera del crate; retirar exportaciones o anticipaciones sin consumidor al simplificar sus módulos |
| `src/main.rs` y `src/diagnostics.rs` | A | Conservar CLI y diagnósticos; reducir comandos/campos al recorrido F1 real |
| `configuration.rs`, `adapters.rs`, `application.rs` | A/D | Conservar declared/validated/effective; retirar semántica futura y autoridades duplicadas |
| `rules.rs` | A/B/D | Conservar solo reglas con consumidor F1; trasladar ejes y jerarquías futuros a sus fases y la auditoría conjunta a F9.7 |
| `numeric.rs` | B | Trasladar la normalización DTC-049 a F6.2; conservar temporalmente hasta retirar sus consumidores del checkpoint |
| `randomness.rs`, BLAKE3, ChaCha12 y raíz | A | Conservar; runtime F1 consume generación y reducción |
| catálogos de subsemillas v1–v3 | A/B | Conservar bytes publicados; F1 no ejecuta ámbitos sin consumidor |
| `identity.rs` | A | Conservar IDs producidos por configuración, ejecución, referencia y artefacto |
| `reference.rs`, `reference_fasta.rs`, `reference_generation.rs` | A | Conservar generación, validación y copia exacta FASTA |
| `exact_decimal.rs`, `coordinates.rs`, `reference_reduction*.rs` | A | Conservar para F1.6 y desacoplar del grafo retirado |
| `integrity.rs` | A | Conservar tamaño y SHA-256 por fichero |
| `execution.rs` | A/B/D | Conservar la orquestación de referencia y publicación; retirar eventos, jerarquías y validaciones duplicadas sin consumidor F1 |
| `publication.rs` | A/C/D | Conservar provisional, INCOMPLETE y rename no-clobber; retirar hardlink/inodo/TOCTOU/capabilities de atacante activo |
| `results.rs`, esquemas de resultados | A/D | Sustituir por manifiesto general y catálogo simple; retirar enlaces/unidades duplicados |
| `manifest_schema.rs` y `jsonschema` runtime | C/D | Retirados en F1.4.5; los esquemas se validan mediante CPython/jsonschema externo en pruebas |
| `ResolutionDomain`, `resolution-evidence-v1.schema.json` | D/B | Eliminar de F1; reintroducir solo con una fuente real en la fase propietaria |
| matriz completa de reglas/ejes y normalización poblacional | B | Trasladar a F2/F6/F7; auditoría completa en F9.7 |
| eventos técnicos y manifiestos de evento | B | Introducir en F3/F4 con sus primeros consumidores |
| tests de symlink/hardlink/inodo/TOCTOU hostil | C | Trasladar al hardening F9/CI, no a la puerta F1 |
| tests 0/1/N eventos × 0/1/2 reducciones | D | Sustituir por casos F1 reales sin producto cartesiano ficticio |

## Artefactos no Rust

| Rutas | Clase | Destino en F1 0.2 |
| --- | --- | --- |
| `catalogs/diagnostics-v*.tsv`, `id-classes-v1.tsv` | A | Conservar los contratos publicados usados por CLI, IDs y manifiesto mínimo |
| `catalogs/parameter-families-v1.tsv` | A/B | Conservar la compatibilidad publicada; F1 solo consume las familias que permanezcan tras F1.2 |
| `catalogs/subseed-scopes-v1.tsv` a `v3.tsv` | A/B | Conservar bytes y vectores; activar solo generación y reducción en F1.3 |
| `schemas/config-general-v1.md`, `config-scenario-v1.md` | A/B | Mantener el esquema general F1 y conservar el esquema de escenario publicado para F2; el runtime F1 no carga ni compone escenarios |
| `schemas/diagnostics-v1.schema.json` | A | Conservar y probar externamente |
| `schemas/effective-config-v1.schema.json` | A/D | Simplificar a consumidores F1 en F1.2 |
| `schemas/manifest-general-v1.schema.json`, `results-common-v1.schema.json` | A/D | Simplificar al manifiesto general y catálogo únicos en F1.4 |
| `schemas/manifest-event-v1.schema.json` | B | Trasladar a F3/F4 con el primer evento real |
| `schemas/reference-reduction-map-v1.schema.json` | A | Conservar para F1.6 |
| `schemas/resolution-evidence-v1.schema.json` | D/B | Eliminar de F1.2; solo reintroducir con una fuente real propietaria |
| `scripts/check-independent-validators.sh` | A/C | Conservar comprobaciones de formato e integridad; retirar casos exclusivos de hardening activo |
| `scripts/run-f1-reference-benchmarks.sh` | A | Conservar como observación no normativa de F1.5 |
| `scripts/package-f1-distribution.sh` | B/C | Trasladar empaquetado y OCI a F9; no usarlo para cerrar el runtime mínimo |

## Familias de pruebas

| Familia | Clase | Destino en F1 0.2 |
| --- | --- | --- |
| admisión TOML, configuración declarada y diagnósticos CLI | A | Reescribir contra la superficie F1.2 y conservar negativos de formato |
| `rule_*`, catálogo de parámetros y configuración efectiva | A/B/D | Mantener solo reglas consumidas en F1; trasladar matrices jerárquicas y combinatorias |
| semilla y derivación de subsemillas | A/B | Conservar vectores publicados; ejecutar solo los dos ámbitos F1 |
| identidad e integridad | A | Conservar vectores de IDs realmente producidos, tamaño y SHA-256 |
| FASTA, referencia simulada/suministrada e índice | A | Conservar y repetir en F1.5 |
| decimal, coordenadas y reducción | A | Conservar y desacoplar en F1.6 |
| normalización numérica poblacional | B | Trasladar a F6.2 |
| publicación, manifiestos y trazabilidad F1.4/F1.6 | A/C/D | Reescribir contra el paquete mínimo; retirar evento ficticio y threat model activo |
| esquemas de autoría y catálogos normativos | A/B | Mantener compatibilidad publicada y ajustar las afirmaciones al alcance vigente |

## Dependencias

| Crate | Consumidor F1 0.2 | Decisión |
| --- | --- | --- |
| `blake3` | subsemillas e IDs | conservar |
| `getrandom` | semilla raíz generada | conservar |
| `jsonschema` | sin consumidor runtime | retirado de `Cargo.toml` y `Cargo.lock` en F1.4.5; CPython/jsonschema permanece como herramienta externa de pruebas |
| `noodles-fasta` | FASTA e índice de aceptación | conservar |
| `rand_chacha` | generación/reducción reproducibles | conservar |
| `rustix` | rename no-clobber y operaciones Linux mínimas | conservar con features mínimas |
| `serde_json` | configuración efectiva, mapa y manifiesto | conservar |
| `sha2` | integridad publicada | conservar |
| `toml` | configuración de autoría | conservar |

## Pruebas y evidencia

Se conservan pruebas unitarias de TOML, semilla, FASTA, SHA, decimal,
coordenadas y reducción. Las pruebas de configuración efectiva y publicación se
reescriben contra el contrato reducido. La evidencia histórica
`evidence/f1.6/` permanece asociada al checkpoint y no se reutiliza para cerrar
F1 0.2 salvo que el caso siga existiendo y vuelva a ejecutarse.
