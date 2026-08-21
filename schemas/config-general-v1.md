# Esquema de autoría TOML: configuración general v1 en F2.2.1

La vista estructural versionada que puede consumir la interfaz local se publica
en `config-general-v1.schema.json`. Describe forma, tipos, dominios básicos y
anotaciones de autoría, pero no sustituye la admisión ni la validación semántica
del binario Rust.

F2.2.1 conserva un único documento general, con cabecera obligatoria:

```toml
config_kind = "general"
schema_version = "1"
```

La raíz es cerrada y contiene esas claves, `seed`, `ploidy` y
`normalization_tolerance` opcionales, `reference` obligatorio,
`reduced_references` opcional y `population_snvs`
opcional. `ploidy` es una
ploidía global uniforme entera entre 1 y 10; al omitirse conserva el valor
efectivo diploide (`2`) para mantener compatibles los TOML F1. Rust valida el
rango y ninguna secuencia ni individuo puede sustituirla localmente. Los campos
desconocidos y las claves duplicadas se rechazan con su ubicación. En
particular, F2.2.2 no admite escenarios, eventos, individuos configurables,
poblaciones simuladas, fases configurables, reglas ni parámetros de generación
o frecuencias. La única excepción estructural es `[population_snvs]`: declara
el origen de una matriz común y F2.2.2 valida sus sumas ALT multialélicas.

`seed`, cuando aparece, contiene 256 bits como 64 hexadecimales minúsculos.

`normalization_tolerance`, cuando aparece, conserva su literal decimal TOML;
Rust lo interpreta como `f64` finito en `[0,1e-6]`. Si se omite, su valor
efectivo es `1e-9`. La interfaz solo conserva el literal y no replica esa
validación semántica.

Una referencia simulada tiene esta forma:

```toml
[reference]
kind = "simulated"
id = "ref-main"

[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 1000
```

`composition` es una tabla cerrada que exige A, C, G y T exactamente una vez.
Cada valor es un `f64` finito en `[0,1]`. La suma determinista de Neumaier debe
diferir de 1 como máximo `1e-9`; los cuatro valores efectivos se dividen por
esa suma antes de generar. Esta normalización es específica de la referencia y
no introduce un modelo de frecuencias poblacionales.

Una referencia suministrada cambia `kind` por `supplied` y exige `fasta`; no
admite `composition`. Cada secuencia exige `id`, `class` (`chromosome` o
`scaffold`) y `length`. En una referencia suministrada, `length` es una
aserción que se contrasta posteriormente con el FASTA, no otra autoridad.

Cada elemento de `reduced_references` exige `id`, `source_reference_id`,
`region_size_bases` positivo y exactamente uno entre `target_bases` positivo o
`reduction_percent`. Esta operación solo se admite cuando la referencia inicial
es `kind = "supplied"`; una referencia simulada no puede ser fuente de una
reducción en F1. La referencia fuente es la referencia inicial del mismo
documento y el ID derivado debe ser distinto.

La tabla opcional `[population_snvs]` fija una fuente excluyente de matriz SNV:
`mode = "supplied"` exige solo `input_path`; `mode = "reuse"` exige solo
`metadata_path`; `mode = "generate"` no admite rutas. Los parámetros de
generación no se admiten aún: se incorporan en F2.2.5. La tabla canónica
`population_snvs.tsv` usa exactamente las ocho columnas ordenadas
`variant_id`, `contig`, `position`, `ref`, `alt`, `population_id`, `alt_af` y
`source_id`; VCF no es una entrada admitida. Su JSON asociado declara la clase
`population_snvs_metadata`, versión `1`, coordenadas `base1_closed`, referencia,
enlace del TSV, poblaciones y procedencia `supplied` o `generated`. F2.2.2
interpreta cada `alt_af` como `f64` finito en `[0,1]`, suma por
población×locus en orden canónico con Neumaier y aplica la tolerancia efectiva:
si la suma ALT no supera 1 conserva el resto REF; un exceso dentro del borde
inclusivo se normaliza y registra; otro exceso se rechaza. La validación de
FASTA, hashes e IDs canónicos sigue perteneciendo a los pasos posteriores.

El borrador de escenarios permanece en `config-scenario-v1.md` como artefacto
documental futuro; el binario F1 no lo carga ni lo interpreta y F2 tampoco lo
incorpora antes de que exista un consumidor biológico real.
