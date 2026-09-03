# Contratos TSV F5 v1

Todos los TSV F5 v1 son UTF-8 sin BOM, separador TAB, LF final obligatorio,
cabecera única, sin comillas, comentarios, líneas vacías ni columnas extra. Los
enteros son decimales canónicos sin signo y los `f64` usan el literal mínimo de
ida y vuelta. Cada fichero se ordena por la clave indicada.

Cada enlace de manifiesto es un objeto cerrado `artifact_id`, `artifact_kind`,
`schema_version="1"`, `path`, `size_bytes:u64` y `sha256` hexadecimal en
minúsculas. La ruta es literalmente `outputs/<fichero>` y no admite `..`.
Cada fila contiene exactamente las columnas de su cabecera —ni una menos ni una
más—; el layout tiene al menos un bloque, sitios/haplotipos/frecuencias y
fundadores tienen una fila por entidad materializada, y `attempts` una fila por
intento ejecutado. Un éxito contiene exactamente los siete enlaces, una vez y
en el orden de la tabla; un fallo no contiene ninguno de ellos.

| Fichero | artifact_id / kind | Ruta literal | Cabecera exacta | Clave/orden |
|---|---|---|---|
| `f5-block-layout.tsv` | `f5-01-block-layout-v1` / `f5_block_layout` | `outputs/f5-block-layout.tsv` | `block_id\tsequence_id\tstart_1based\tend_1based\tlength_bp` | FASTA, inicio, fin, ID |
| `f5-sites.tsv` | `f5-02-sites-v1` / `f5_sites` | `outputs/f5-sites.tsv` | `block_id\tsite_index\tsite_id\tsequence_id\tposition_1based\tref\talt` | bloque, índice |
| `f5-haplotypes.tsv` | `f5-03-haplotypes-v1` / `f5_haplotypes` | `outputs/f5-haplotypes.tsv` | `population_id\tblock_id\thaplotype_id\thaplotype_rank\tq\talleles` | población, bloque, rango |
| `f5-frequencies.tsv` | `f5-04-frequencies-v1` / `f5_frequencies` | `outputs/f5-frequencies.tsv` | `population_id\tblock_id\tsite_id\tpolarity\ttarget_alt_af\trealized_alt_af\tabs_error` | población, bloque, sitio |
| `f5-ld-evaluation.tsv` | `f5-05-ld-evaluation-v1` / `f5_ld_evaluation` | `outputs/f5-ld-evaluation.tsv` | `population_id\tblock_id\tdistance_start_bp\tdistance_end_bp\tstatus\tpair_count\tmean_r2\ttarget_mean_r2\tmax_abs_error\tabs_error` | población, bloque, inicio de bin |
| `f5-founders.tsv` | `f5-06-founders-v1` / `f5_founders` | `outputs/f5-founders.tsv` | `population_id\tindividual_id\tsequence_id\tcopy_index\tblock_id\tfounder_copy_id\tblock_assignment_id\thaplotype_id` | población, individuo, secuencia, copia, bloque |
| `f5-attempts.tsv` | `f5-07-attempts-v1` / `f5_attempts` | `outputs/f5-attempts.tsv` | `population_id\tblock_id\tattempt_index\tproposals\taccepted\tduplicate_rows\tbins_outside\tsigns_failed\tmax_r2_error\tsum_r2_error\tmax_marginal_error\tselected_best` | población, bloque, intento |

`alleles` tiene exactamente S caracteres `0`/`1`; `polarity` es `positive` o
`negative`; `status` es `evaluated` o `not_applicable`, único caso en que
`mean_r2` y `abs_error` son `NA`; `selected_best` es `true`/`false` y hay uno
solo `true` por bloque×población.
