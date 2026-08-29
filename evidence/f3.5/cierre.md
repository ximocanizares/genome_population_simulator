# Cierre F3.5.4

La prueba `tests/f3_end_to_end.rs` publica un caso de 150 bases con FASTA
suministrado, un único individuo diploide y un SNP directo `A>G` en fase 0 y
posición 75. La matriz suministrada fija un `A>T` en la posición 25, distinta:
así no puede ocultar el efecto del evento directo. Reconstruye `truth.jsonl`
contra la referencia publicada y comprueba una copia con ambos cambios y otra
solo con el cambio de matriz, sin consumir FASTQ. Después valida el FASTQ mediante la
reproducción independiente de su origen y los esquemas JSON Schema externos de
los manifiestos general v3 y de evento v1.

El mismo caso incluye una `N` manual en la referencia: queda elegible, aparece
literalmente en el FASTQ con Phred 0 y el validador independiente la contabiliza.
También calcula por separado la cobertura WGS básica (un fragmento de 150 bases
sobre 150 bases haploides) y la cobertura del FASTQ final publicada. Ambas dan
1x en F3 porque su contrato mínimo no introduce descartes ni lecturas ajenas;
la igualdad no se usa como una sustitución de la cobertura alineada externa.
El manifiesto del evento se contrasta además con el evento, referencia,
individuo y artefacto publicados, además de sus validadores de esquema.

Al corromper el FASTQ después de la publicación, la verdad continúa siendo
reconstruible y la validación independiente del FASTQ falla. Esto conserva la
dirección `configuración → referencia → individuo → verdad → FASTQ`: la
evidencia no crea ni repara la verdad biológica.

Comprobación focal:

```text
cargo test --test f3_end_to_end --locked
```

El caso no añade perfiles, pares, mapeadores, alineamientos, escenarios,
persistencia, empaquetado ni hardening. Es entrada FASTQ válida para un
mapeador o detector externo, que permanece fuera de F3.
