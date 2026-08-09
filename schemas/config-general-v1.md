# Esquema de autoría TOML: configuración general v1

Todo documento general declara, antes de cualquier contenido semántico:

```toml
config_kind = "general"
schema_version = "1"
```

Las claves raíz admitidas en F1.2 son `config_kind`, `schema_version`, `route`,
`normalization_tolerance`, `reference`, `scenarios`, `events`, `individuals`,
`populations`, `subpopulations`, `phases`, `groups`, `blocks`, `rules` y
`frequency_sets`. No se admiten claves desconocidas ni duplicadas.

- `route` es exactamente `direct`, `evolution` o `families`.
- `reference` es una tabla. `kind = "simulated"` exige `id`, `default_ploidy`,
  `sequences` y `composition`; `kind = "supplied"` exige `id`, `fasta` y
  `sequences`. Cada secuencia es una tabla con `id`, `class`, `length` y
  `ploidy` opcional. `class` es `chromosome` o `scaffold`.
- `scenarios` es una lista de tablas con `id`, `path` y `route`. La vía directa
  exige cero escenarios; las vías evolutiva y de familias exigen al menos uno
  de su misma vía.
- `events`, `individuals`, `populations`, `subpopulations`, `phases`, `groups`
  y `blocks` son listas de tablas con un `id` explícito. No se crea ningún
  evento implícito.
- `rules` es una lista de tablas con `id`, `parameter`, `value` y cero o más
  selectores. Los selectores admitidos por familia están en el catálogo de
  F1.2; una regla sin selectores es la regla base. Los nombres de selector son
  `individual`, `population` o `subpopulation`, `phase` o `generation`,
  `reference` con `sequence` y opcionalmente una región o `locus`, `group`,
  `block` y `technical_event`. Los dos niveles de un mismo eje son excluyentes.
  El tipo y dominio de `value` dependen de `parameter`: `ploidy` es un entero
  entre 1 y 10 inclusive; `selfing_rate`, `error_rate` y
  `contamination_rate` son `f64` finitos en `[0,1]`; `variant_rate`,
  `recombination_rate` y `coverage` son `f64` finitos no negativos; y
  `output_enabled` es booleano. Un entero decimal es admisible como valor real
  cuando el parámetro exige `f64`.
- `frequency_sets` es una lista de tablas identificadas. El modo
  `probabilistic` usa `values = [{ id, category, value }, ...]`; el modo
  `controlled` usa `counts = [{ id, category, count }, ...]`. Las formas no se
  pueden combinar. `tolerance` solo pertenece al modo probabilístico.

`normalization_tolerance` y cada `frequency_sets.tolerance` son números finitos
en `[0, 1e-6]`; la tolerancia local sustituye a la general y, si ambas faltan,
se usa `1e-9`. Los valores probabilísticos son finitos en `[0,1]` y se admiten
solo si su suma Neumaier queda dentro del borde inclusivo de la tolerancia.

Las coordenadas de región son `region_start` y `region_end`, basadas en 1 y
cerradas, junto con `reference` y `sequence`. Las listas con orden normativo
conservan su orden; ninguna otra posición física concede identidad o precedencia.
