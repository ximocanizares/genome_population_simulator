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

## Orden de los cambios de base de un fragmento

En cada fragmento vigente de tipo `source`, los cambios de base (`snps`) se
publican en orden **estrictamente creciente de posición**, sin posiciones
repetidas, con independencia del orden en que entraron en el historial de la
copia. El validador semántico replica el historial insertando cada cambio en su
sitio y compara la lista resultante con la publicada, igual en bases, posiciones
y longitud.

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

## El locus de un caso *indel*

**Desde el 2 de octubre de 2026** el `locus` de un caso `indel` se publica en
coordenadas de la **referencia**, como manda `coordinate_convention`, y no en
las coordenadas **vigentes** de la copia con que se declaró el evento. Las dos
coinciden mientras ningún *indel* anterior de la misma copia las desplace; en
cuanto uno lo hace, dejan de coincidir. Hasta esa fecha el locus copiaba las
vigentes: una inserción `AAT` declarada en el punto vigente 175, detrás de una
inserción de dos bases y una deleción de una, quedaba en la base 173 de la
referencia y su caso decía 175, de modo que el VCF no encontraba el caso de
ese ALT y la ejecución abortaba. No cambia la versión: el esquema no fijaba
esa igualdad, la imponía el código, y lo que se publicaba contradecía el
convenio declarado en la cabecera.

El evento fuente del alelo (`source_event`), los `start`/`end` de la deleción
y el historial de la copia **siguen en coordenadas vigentes**: son los datos
con que el evento se aplicó y con los que se replica.

La regla de traducción es una sola y es la misma con la que el VCF ancla un
fragmento insertado:

- **Inserción.** El locus es el punto vacío `[a, a)`, donde `a` es el final de
  la última base de referencia que precede al punto de inserción en la copia
  tal como estaba al aplicarla. Una inserción dentro de otra hereda así el
  ancla de la exterior, y una inserción justo en la junta que deja una
  deleción se ancla a la izquierda de la junta, no a la derecha.
- **Deleción.** El locus va de la primera a la última base de referencia que
  la deleción elimina. Si solo elimina bases insertadas, no ocupa ninguna base
  de referencia y su locus es el punto vacío del ancla de su inicio.

Dos consecuencias de la regla que un consumidor tiene que conocer:

- **El locus de una deleción puede incluir bases que ya había borrado un
  evento anterior.** Es el intervalo entre la primera y la última base de
  referencia que elimina, y entre ellas puede haber un hueco previo: una
  deleción que se lleva las bases de los dos lados de una deleción anterior
  de cinco bases tiene un locus que abarca también esas cinco.
- **El locus de una inserción no tiene por qué caer dentro del `REF` de su
  registro VCF.** El VCF normaliza cada alelo —lo alinea a la izquierda y
  recorta lo común—, de modo que el locus de una inserción puede quedar en el
  borde derecho del `REF` o fuera de él, por ejemplo cuando una deleción
  posterior más a la izquierda desplaza el registro. Para asociar un caso con
  su registro hay que usar `GPS_CASE`, no `POS` ni `REF`.

La estructura de una línea solo puede comprobar la forma —el locus de una
inserción no tiene anchura—. La igualdad exacta la exige el validador
semántico: al replicar el historial de cada copia calcula, con su propia
implementación de la regla, el locus de cada inserción y deleción, y rechaza
un caso cuya copia portadora asigne a su evento un locus distinto del
publicado.

## Un caso sin registro VCF

**Desde el 3 de octubre de 2026** un caso de `truth.jsonl` puede no aportar
ningún ALT al VCF. Ocurre cuando cada copia portadora es igual a la
referencia alrededor del locus: una inserción y una deleción de la misma
copia que se anulan, una inserción que una deleción posterior borra entera, o
un alelo de repetición más una deleción del flanco que juntos reproducen la
referencia. La verdad y el VCF difieren ahí **a propósito**:

- la verdad registra **eventos**: el caso sigue en el fichero con sus
  portadores y su dosis, porque la copia es portadora del alternativo si y
  solo si su historial contiene el evento, y lo contiene;
- el VCF registra **secuencia**: en ese sitio la copia es igual a la
  referencia y no hay ALT que publicar por ese caso.

Lo habitual es que entonces el caso no aparezca en ningún registro. **Sí
aparece** en el `GPS_CASE` de un registro cuando su copia lleva en él otro ALT
cuyo tramo compuesto —o fundido con el de otras copias— contiene el locus del
caso: la relación ALT↔casos se establece por tramo, y el caso queda listado
con ese ALT aunque su efecto propio sea nulo. Un ejemplo: una copia con una
inserción y una deleción que se anulan y, unas bases más allá, otra deleción,
todo dentro de un registro fundido con la deleción grande de otra copia; los
tres casos de esa copia figuran en el `GPS_CASE` del registro.

No cambia la versión ni ningún campo: nada de lo que se escribe es distinto,
solo deja de ser obligatorio que cada caso tenga un registro (excepción de
REQ-VER-25 y nota de DTC-065). El agrupador del productor y el oráculo VCF
independiente aceptan un caso sin ALT **solo** si comprueban que su copia
portadora es la referencia alrededor del locus, de una de dos maneras:

- una ventana que el VCF compuso, y en la que esa copia queda igual a la
  referencia, contiene el locus entero;
- ningún cambio de la copia —hueco, punto de inserción, cambio de base o
  tracto sustituido, leídos de sus piezas finales— toca el locus, extremos
  incluidos. Es la que cubre una inserción borrada entera en un sitio que
  ninguna otra copia cambia: no queda cambio alguno y no se compone ventana.

Un caso cuya copia difiere y no tiene ALT sigue siendo un defecto interno.
