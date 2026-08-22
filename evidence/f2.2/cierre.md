# Cierre F2.2 — matriz SNV y variantes pequeñas

## Alcance comprobado

F2.2 queda cerrada como frontera de matriz poblacional común y de núcleo
transaccional de eventos pequeños. No incorpora otro ejecutor ni configura
destinos de individuo, copia o fase: esos nombres aún no existen en TOML y su
primer consumidor observable es F2.3.

`tests/f2_2_integration.rs` cubre el recorrido completo ya disponible:

- generación reproducible de una matriz canónica con una semilla declarada;
- oráculo independiente que contrasta cada REF publicado con los bytes FASTA,
  la cardinalidad, los alelos y el intervalo de frecuencia;
- reutilización que valida el recurso sin volver a publicarlo ni modificarlo;
- auditoría de normalización por población×locus (originales, sumas, decisión
  y valores efectivos) conservada desde metadata hasta configuración efectiva
  y manifiesto, y rechazada si deja de coincidir con el par canónico; los
  ceros ausentes sólo se admiten para pares canónicos
  `discarded_variant_pairs` de población×variante;
- preparación de una matriz suministrada sin alterar su TSV de origen;
- procedencia congelada (`variant_id`, población y `source_id`) de la fila
  canónica al convertirla en SNP explícito, seguida de un bloque de indels
  aleatorios repetible y de bytes finales contrastados con una verdad literal;
- validación física de los artefactos publicados mediante el validador de
  recursos, no sólo por sus documentos JSON correlacionados.

Los oráculos especializados permanecen en las pruebas propietarias:
normalización (`tests/population_snvs_preparation.rs`), generación/espectros
(`tests/population_snvs_generation.rs`), publicación/reutilización
(`tests/population_snvs_publication.rs`), RNG/Poisson (`tests/f2_rng.rs`), SNP
(`tests/snp_application.rs`), indels (`tests/indel_application.rs`) y
composición (`tests/small_event_composition.rs`). La paridad Marimo↔TOML está
en `tests/test_population_snvs_editor.py`, `tests_ui/test_ui_contract.py` y
`tests_ui/test_process_adapter.py`.

## Validación ejecutada

El cierre se comprobó con `cargo fmt --all -- --check`,
`cargo clippy --workspace --all-targets --all-features --locked -- -D warnings`,
`cargo test --workspace --all-targets --all-features --locked`,
`cargo test --workspace --doc --all-features --locked`,
`cargo doc --workspace --no-deps --all-features --locked` y `git diff --check`.
La suite UI focal también pasó con
`uv run python -m unittest tests/test_population_snvs_editor.py tests_ui/test_ui_contract.py tests_ui/test_process_adapter.py`
(23 pruebas, una omitida por una condición declarada del entorno).

## Límites que se conservan

- Ninguna frecuencia poblacional se interpreta como alelo de un individuo.
- Las filas suministradas/generadas/reutilizadas describen la matriz común; no
  añaden individuos, VCF, red, servidor ni publicación adicional.
- La matriz vacía requiere la autorización explícita ya expuesta por Marimo y
  la CLI.
- F2.3 debe crear identidades de individuo/copia/fase antes de publicar en
  TOML controles de SNP o indels dirigidos.
