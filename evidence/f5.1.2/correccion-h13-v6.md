# Corrección H13 — admisión de entradas F5/v6

## Hallazgo y decisión

| ID | Gravedad | Decisión inicial | Estado final | Archivos | Pruebas | Seguimiento |
|---|---|---|---|---|---|---|
| H13 | Grave | Aplicar con adaptación | Aplicado con adaptación | `application.rs`, `f5_profile.rs`, `f5_preflight.rs`, esquema efectivo y pruebas F5 | CLI: límite físico TOML exacto/−1; metadata `source_citation` con límite de memoria exacto/−1 e invariancia de layout | Ninguno |

La adaptación evita usar `F5OtherResources`, que por contrato representa solo
recursos ajenos a F5 y que el CLI entrega vacío. `F5InputDocuments` obtiene
por `metadata` las longitudes del general, escenario, metadata y los cinco TSV
antes de abrir el metadata o materializar un TSV. El gate de admisión rechaza
`2*sum(documentos)` si no cabe en `max_working_memory_bytes`; el factor dos
cubre bytes crudos retenidos y parser/búfer. F5 limita además cada TOML a
65536 bytes antes de materializarlo.

## Fórmula y compatibilidad

`f5-resource-accounting-v6` añade:

```text
input_document_bytes = general_toml_bytes + scenario_toml_bytes
  + profile_metadata_bytes + profile_tsv_bytes
input_parser_peak_bytes = 2 * input_document_bytes
operation_units += input_document_bytes
working_memory_bytes += input_parser_peak_bytes
```

Los cinco escalares nuevos son `u64` cerrados en
`effective-config-v5.schema.json`; el manifiesto reutiliza esa definición. V6
no cambia RNG, identidad semántica, geometría ni layout: únicamente el cálculo
de recursos y su auditoría.

## Gates

- Un comentario TOML de 65536 bytes es admitido y uno de 65537 recibe
  `GPS-RES-101` por CLI antes de materializarlo.
- Un `source_citation` mayor incrementa exactamente en dos veces su delta la
  reserva de parser; el límite igual se admite y el límite menos uno se rechaza
  por CLI antes de abrir el metadata/TSV.
- Los layouts de paquetes admitidos que solo difieren en comentario TOML o
  cita son idénticos con la misma semilla.

Riesgo residual: la ruta F5 exige sus discriminadores de cabecera dentro de la
sonda CLI acotada; el parser TOML completo sigue siendo la autoridad de
sintaxis y semántica.
