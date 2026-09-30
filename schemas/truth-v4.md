# Contrato `truth.jsonl` v4 y v5

Estas son las dos versiones vigentes de la verdad estructurada del recorrido
general. Sustituyen el 30 de septiembre de 2026 a la versión `"1"`
(`truth-v1.schema.json`, descrita en `truth-v1.md`) y a la `"3"`
(`truth-v3.schema.json`), que quedan como registro histórico, sin lector
activo. La verdad de F5 sigue siendo la `"2"` (`truth-v2.schema.json`) y no
cambia.

| Versión | Esquema | Cuándo se publica | Sustituye a |
|---|---|---|---|
| `"4"` | `truth-v4.schema.json` | la verdad no publica ningún alelo de repetición | `"1"` |
| `"5"` | `truth-v5.schema.json` | la verdad publica al menos un alelo de repetición | `"3"` |

La versión se promueve **por lo que se publica, no por lo que se declara**,
igual que entre la `"1"` y la `"3"`: un escenario que declara un locus de
repetición cuyo único alelo es el de la referencia no publica ninguna variante
de repetición y sigue siendo `"4"`. El validador lo exige en las dos
direcciones. Una versión ausente, anterior o desconocida se rechaza, sin
intentar reinterpretarla ni migrarla.

## Qué cambia respecto a las versiones 1 y 3

**Una sola cosa: los alelos de inserción y de deleción ya no llevan
`final_haplotype_bases`.** Ese campo copiaba en cada registro de genotipo la
secuencia entera de la copia portadora, de modo que el fichero crecía como el
número de *indels* multiplicado por la longitud del genoma: con un cromosoma de
500 kb y diez deleciones, `truth.jsonl` medía 5.018.930 bytes y el 99,6 % era
ese campo; con un cromosoma humano, cada *indel* de cada individuo habría
añadido unos 250 MB. Contradecía REQ-NOF-01 y REQ-NOF-04.

El campo era además **redundante** como dato: cada copia se publica ya en su
registro `haplotype_copy`, y la secuencia final se reconstruye desde él y la
referencia. No se sustituye por nada en el fichero, ni por una huella por copia
ni por una ventana local. Los esquemas nuevos lo rechazan como propiedad
desconocida y el lector de la verdad también.

Lo que **no** era redundante es el contraste que hacía de rebote. La secuencia
del campo la materializaba el productor desde su estado interno, y el validador
la comparaba con la que él reconstruía desde la verdad: era el único punto en
que se contrastaban las dos. Sin él, una proyección que publicara bases
distintas de las simuladas, de forma coherente en todos sus registros, pasaría
todas las demás comprobaciones —las lecturas saldrían de una secuencia y la
verdad y el VCF describirían otra—. Ese contraste **se conserva en memoria**:
al publicar, el productor compara, para **cada** copia y no solo las que llevan
*indels*, la huella SHA-256 de su materialización con la de la reconstrucción
del validador independiente, y la ejecución falla con un defecto interno si
difieren. No se escribe nada en `truth.jsonl`.

Todo lo demás es exactamente lo que describe `truth-v1.md` para la versión
`"1"` —cabecera, copias, fragmentos vigentes, historial, procedencias, orden
canónico y relaciones del conjunto— y, para la `"5"`, lo que añadió la `"3"`:
el fragmento vigente y el evento de historial `repeat_allele`, el caso
genotípico `repeat` y su alelo alternativo `repeat`, con su secuencia completa
y su estructura en bloques.

## El alelo estructural

Un caso `indel` usa `reference` y una inserción o deleción tipada. El alelo
estructural enlaza su **evento fuente** y sus **`piece_links`**: los fragmentos
vigentes de la copia portadora que todavía proceden de esa inserción, en
offsets de la secuencia final de la copia. Una deleción no enlaza ningún
fragmento. Una inserción parcialmente eliminada es válida si conserva enlaces
no vacíos.

El validador semántico reconstruye cada copia **sin usar código del
productor**: reproduce sus segmentos y su historial sobre la referencia y exige
que el resultado coincida con sus `current_pieces`. Después exige, para cada
copia del individuo y la secuencia del caso, que sea portadora del alternativo
si y solo si:

- su historial contiene el evento con el mismo identificador y los mismos
  datos (punto y bases de la inserción; inicio y fin de la deleción);
- el evento fuente del alelo coincide con ese evento;
- los `piece_links` del alelo son exactamente los fragmentos insertados con esa
  procedencia en la copia reconstruida, con sus offsets; para una deleción, la
  lista vacía.

La cota de los `piece_links` —que no se salgan de la copia— la implica esa
igualdad con la reconstrucción; antes se comprobaba además contra la longitud
de `final_haplotype_bases`, que ya no existe.
