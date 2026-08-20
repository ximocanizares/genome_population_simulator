# Esquema de autoría TOML: configuración general v1 en F1

F1 admite un único documento general, con cabecera obligatoria:

```toml
config_kind = "general"
schema_version = "1"
```

La raíz es cerrada y contiene exactamente esas claves, `seed` opcional,
`reference` obligatorio y `reduced_references` opcional. Los campos
desconocidos y las claves duplicadas se rechazan con su ubicación. En
particular, F1 no admite rutas biológicas, ploidía, escenarios, eventos,
individuos, poblaciones, fases, reglas ni conjuntos de frecuencias.

`seed`, cuando aparece, contiene 256 bits como 64 hexadecimales minúsculos.

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

El esquema de escenarios publicado permanece en `config-scenario-v1.md` como
artefacto documental para F2; el binario F1 no lo carga ni lo interpreta.
