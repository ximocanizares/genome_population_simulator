# Esquema de autoría TOML: configuración general v2

La vista estructural activa para la interfaz local se publica en
`config-general-v2.schema.json`. El binario Rust conserva la admisión y la
validación semántica final.

```toml
config_kind = "general"
schema_version = "2"
```

La raíz es cerrada. `reference` es obligatorio; `seed`, `ploidy`,
`normalization_tolerance`, `reduced_references`, `population_snvs`,
`initial_populations`, `direct_events`, `technical_event` y `vcf_output` son
opcionales. Los campos desconocidos y las claves duplicadas se rechazan.

## Evento técnico F4

`[technical_event]` es el único evento técnico paired-end. Exige `event_id`,
`individual_ids`, `target_coverage_wgs` finita y mayor o igual que cero, y al
menos una entrada `[[technical_event.libraries]]`. Cada biblioteca declara
`individual_id`, `i7` e `i5`; ambos componentes UDI son exactamente ocho bases
mayúsculas A/C/G/T. `read_length_bases`, si aparece, está entre 1 y 150; al
omitirse se resuelve a 150. `chimera_rate`, si aparece, es una fracción finita
en `[0,1]`; al omitirse se resuelve a `0`. `index_hopping_rate`, si aparece,
es una fracción finita en `[0,1]`, se resuelve a `0` al omitirse y se aplica
después de la contaminación biológica. Los campos planos
`contamination_min_percent` y `contamination_max_percent` son porcentajes
finitos en `[0,100]`, se resuelven ambos a `0` al omitirse y deben respetar
`min <= max`. `contamination_recipient_ids` es una lista opcional, no vacía y
sin repetidos de individuos seleccionados; al omitirse se selecciona toda la
población del evento. No existe una tabla anidada de contaminación.
`read_count` no forma parte de v2 y se rechaza.
`library_id` se deriva de forma determinista por Rust y no se autoriza en TOML.
`sequencing_errors` admite `"enabled"` (predeterminado) o `"disabled"`.
El modo `disabled` es el control PE neutro de VE-D02: no representa una
probabilidad cero ni usa `error_scale = 0`; emite exactamente cero
sustituciones en R1 y R2, conserva los Phred del perfil y publica en la
configuración efectiva y el manifiesto que esas calidades no representan una
probabilidad de error realizada.

```toml
[technical_event]
event_id = "reads-main"
individual_ids = ["sample-A"]
target_coverage_wgs = 30
sequencing_errors = "disabled" # control PE neutro; omitir para el modo actual
read_length_bases = 150
chimera_rate = 0
index_hopping_rate = 0
contamination_min_percent = 0
contamination_max_percent = 0
# contamination_recipient_ids = ["sample-A"]

[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
```
