# Esquema de autoría TOML: configuración general v1 en F2.3.4

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
`reduced_references` opcional, `population_snvs`
opcional, `initial_populations` opcional y `direct_events` opcional. `ploidy` es una
ploidía global uniforme entera entre 1 y 10; al omitirse conserva el valor
efectivo diploide (`2`) para mantener compatibles los TOML F1. Rust valida el
rango. Las cohortes directas de F2.3.1 declaran su propia ploidía efectiva; no
existen sustituciones por secuencia ni por individuo. Los campos
desconocidos y las claves duplicadas se rechazan con su ubicación. En
particular, no admite escenarios, poblaciones evolutivas ni reglas. Además de las cohortes directas,
`[population_snvs]` declara
el origen de una matriz común y F2.2.2 valida sus sumas ALT multialélicas.

F2.3.4 declara destinos explícitos de copia con `[[direct_events]]`, resuelve
individuo y fase contra las cohortes y preflighta SNP e indels antes de abrir la
publicación. No materializa ni publica individuos y Marimo no recibe controles
nuevos. `--allow-incompatible-genotypes` solo conserva y registra una
restricción SNP de probabilidad cero; no es `--force` ni relaja referencia,
alelos, coordenadas o ploidía. `--allow-empty-snv-matrix` conserva su alcance
independiente.

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
`metadata_path`; `mode = "generate"` no admite rutas y exige
`[population_snvs.generate]`. Este bloque usa `count_mode = "exact"` con
`locus_count`, o `count_mode = "poisson"` con exactamente uno de `mean_loci` y
`mean_density_per_base`; fija `spatial_mode = "uniform"` o
`"weighted_regions"`, y selecciona explícitamente `substitution_profile`:
`uniform_snv12_v1` (predeterminado) o
`human_harris2017_snv12_approx_v1`. El segundo es una aproximación versionada
de los doce cambios SNV dirigida por el recurso Harris 2017; no constituye una
tasa universal.

Cada población declara `population_id` y uno de estos modelos de frecuencia:
`frequency_model = "uniform"` con `min_alt_af` y `max_alt_af` (el intervalo
heredado semiabierto `[min_alt_af,max_alt_af)`); `frequency_model =
"spectrum"` con una lista no vacía `spectrum`; o el perfil explícito
`human_1000g_global_sfs_approx_v1`. Este último es una aproximación global SFS
versionada, nunca una elección implícita. Los modelos `spectrum` y
`human_1000g_global_sfs_approx_v1` admiten `frequency_mode = "stochastic"`
(predeterminado) o `"exact"`: el primero sortea la componente por VariantId y
el segundo fija sus cuotas mediante mayores restos, con una permutación
determinista por VariantId. El modelo `uniform` debe omitir `frequency_mode`.

Cada entrada de `spectrum` declara un `weight` finito no negativo (al menos
uno debe ser positivo) y es `kind = "fixed"` con `fixed_alt_af` en `[0,1]`, o
`kind = "interval"` con `min_alt_af`, `max_alt_af` y `sampling`. Los
intervalos son semiabiertos, no se solapan y usan `sampling =
"linear_uniform"` o `"log_uniform"`; el muestreo logarítmico exige mínimo
estrictamente positivo. Los valores `0` y `1` solo se expresan con componentes
fijos. Las regiones ponderadas incluyen `outside_weight`, `contig`, `start`,
`end` y `weight` en base 1 cerrada. La tabla canónica
`population_snvs.tsv` usa exactamente las ocho columnas ordenadas
`variant_id`, `contig`, `position`, `ref`, `alt`, `population_id`, `alt_af` y
`source_id`; VCF no es una entrada admitida. Su JSON asociado declara la clase
`population_snvs_metadata`, versión `1`, coordenadas `base1_closed`, referencia,
enlace del TSV, poblaciones, procedencia `supplied` o `generated` y una
auditoría por `población×locus`. La auditoría conserva los ALT originales y su
suma, la decisión (`reference_remainder` o `normalized`), los valores
efectivos, su suma y la frecuencia REF resultante. La misma auditoría se enlaza
sin cambios en la configuración efectiva y el manifiesto. La metadata además
enumera de forma canónica `discarded_variant_pairs`: un ALT efectivo cero que
ya no está en el TSV sólo es admisible si su par `population_id`×`variant_id`
pertenece exactamente a esa lista. Cada variante descartada debe aparecer una
vez para cada población declarada, y su `variant_id` se deriva de la
referencia, contig, posición, REF y ALT. Al
reutilizar, Rust contrasta identidad, orden, alelos, TSV, tolerancia y metadata.
F2.2.2
interpreta cada `alt_af` como `f64` finito en `[0,1]`, suma por
población×locus en orden canónico con Neumaier y aplica la tolerancia efectiva:
si la suma ALT no supera 1 conserva el resto REF; un exceso dentro del borde
inclusivo se normaliza y registra; otro exceso se rechaza. La validación de
FASTA, hashes e IDs canónicos sigue perteneciendo a los pasos posteriores.

## Cohortes iniciales directas (F2.3.1)

Cada `[[initial_populations]]` es una cohorte de creación directa, no una
población evolutiva ni reproductiva. Selecciona una `population_id` de la matriz
canónica cuando esta se consuma en F2.3.2, declara `individual_count` positivo,
la `ploidy` efectiva uniforme de la cohorte entre 1 y 10 y `genotype_model`.
`hwe` admite toda ploidía; `f_diploid` exige `ploidy = 2` e
`inbreeding_coefficient` finito en `[0,1]`. HWE no admite ese coeficiente. Esta
fase solo valida y conserva el modelo: no sortea genotipos ni fase.

`individual_ids` es una lista opcional no vacía de IDs aportados. Si existe, debe
contener exactamente `individual_count` IDs válidos y cada uno debe ser único en
toda la ejecución. JSON Schema impone que la lista no esté vacía y que sus
elementos no se repitan, pero la igualdad entre su longitud y
`individual_count` es una semántica exclusiva de la validación Rust. Si se
omite, Rust deriva un ID estable por individuo con DTC-039,
usando como ámbito `population_id` y una clave fundadora canónica de ordinal
uno-basado. Ese ordinal pertenece a la cohorte, no al orden físico del TOML, de
modo que reordenar las cohortes no cambia los IDs derivados y aumentar N conserva
los anteriores. Las IDs aportadas y derivadas nunca se mezclan en una misma
cohorte.

```toml
[[initial_populations]]
population_id = "P1"
individual_count = 2
ploidy = 2
genotype_model = "f_diploid"
inbreeding_coefficient = 0.25

[[initial_populations]]
population_id = "P2"
individual_count = 2
ploidy = 4
genotype_model = "hwe"
individual_ids = ["sample-A", "sample-B"]
```

Marimo no muestra aún estos campos: sus controles se incorporan exclusivamente
en F2.3.6. La configuración efectiva publica las cohortes, sus IDs ya resueltos,
ploidía y modelo para que F2.3.2 sea el primer consumidor biológico.

## Cambios explícitos por individuo (F2.3.4)

Cada `[[direct_events]]` exige `event_id`, `individual_id`, `phase`, `kind` y
`contig`. `kind = "snp"` exige exactamente `position_1based`, `ref` y `alt`;
`kind = "insertion"`, `point` y `bases`; y `kind = "deletion"`, `start` y
`end`. SNP usa posición base 1 sobre la referencia; los indels usan las
coordenadas del estado vigente del compositor y se aplican en el orden declarado
para ese individuo tras la proyección de SNP de matriz. `individual_id` debe
pertenecer a una cohorte y `phase` debe existir para su ploidía. Se exige una
matriz SNV para recorrer genotipo y fase, cuya identidad y huella canónica se
comprueban antes de los cambios. Tanto `validate` como `run` ejecutan este
preflight de solo lectura; solo `run` abre la publicación.

```toml
[[direct_events]]
event_id = "sample-A-snp"
individual_id = "sample-A"
phase = 0
kind = "snp"
contig = "chr1"
position_1based = 42
ref = "A"
alt = "G"

[[direct_events]]
event_id = "sample-A-insertion"
individual_id = "sample-A"
phase = 0
kind = "insertion"
contig = "chr1"
point = 43
bases = "AC"
```

El borrador de escenarios permanece en `config-scenario-v1.md` como artefacto
documental futuro; el binario F1 no lo carga ni lo interpreta y F2 tampoco lo
incorpora antes de que exista un consumidor biológico real.
