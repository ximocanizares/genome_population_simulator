# Trabajo pendiente

## Mejoras de rendimiento resueltas

- [x] Reducir a la mitad el pico de memoria de una ejecución con referencia
  suministrada. Una ejecución cargaba el genoma tres veces —preparación,
  recarga autoritativa y validación física de la publicación— y dos de esas
  cargas se solapaban, de modo que el pico doblaba el tamaño del genoma de
  partida. Liberar la referencia en los dos puntos donde deja de usarse baja el
  pico de 1.631 MB a 867 MB con el tomate de 759 MB, y el tiempo de CPU un 11 %.
  Los siete artefactos publicados quedan idénticos byte a byte.
  - Hace falta corregir los dos solapamientos a la vez: al eliminar solo uno, el
    otro sostiene el mismo pico y la medición no acusa ninguna mejora.
- [x] Eliminar las validaciones físicas redundantes durante la publicación:
  una única validación completa se realiza tras escribir y releer
  `manifest.json`; antes del renombrado atómico se conservan las comprobaciones
  de catálogo, tamaño, SHA-256 y manifiesto.

## Configuración de cohortes grandes

- [ ] Añadir una forma compacta de declarar fundadores consecutivos mediante
  una raíz y numeración determinista (por ejemplo, `founder_id_prefix = "tom"`
  para producir `tom-0001`, `tom-0002`, … hasta `individual_count`).
  - Debe conservar reproducibilidad, unicidad dentro de la población y la
    derivación estable de IDs técnicos.
  - `individual_ids` y `founder_keys` explícitos deben seguir disponibles para
    nombres no consecutivos o asignados externamente.
  - La sintaxis exacta, los conflictos entre modalidades y la interfaz Marimo
    requieren especificación antes de implementarse.

## Salidas FASTQ

- [x] Publicar los FASTQ finales comprimidos como `.fastq.gz`.
  - El perfil v1 usa gzip/BGZF de `noodles-bgzf 0.47.0`, DEFLATE nivel 6,
    bloques BGZF y `MTIME = 0`; no añade dependencias ni una opción pública de
    compresión. La repetibilidad binaria se limita al mismo simulador, lock y
    plataforma admitida.
  - La escritura es por flujo; se relee y valida el FASTQ descomprimido antes
    de publicar atómicamente. Tamaño y SHA-256 se calculan y verifican sobre
    los bytes gzip finales. Las pruebas incluyen decodificación con `gzip -cd`.

## Controles de error de secuenciación PE

- [ ] Exponer en el TOML del evento técnico F4 escalas de error independientes
  para R1 y R2.
  - Conservarán `1` como perfil actual. El valor `0` no se usará para
    representar ausencia de errores: una probabilidad cero no tiene un Phred
    finito compatible con Q0–Q93.
- [x] Añadir `sequencing_errors = "disabled"` como modo PE neutro. Produce
  exactamente cero sustituciones en ambos mates sin usar `error_scale = 0` y
  registra de forma visible que las calidades FASTQ conservadas no representan
  una probabilidad de error realizada.
  - La configuración efectiva registra los modos solicitado/efectivo y la
    marca de calidad. El manifiesto del evento registra además, por mate,
    `quality_shift = 0`, `error_scale = 1`, `opportunities` y `errors`; estos
    contadores agregados no forman parte de la configuración efectiva y el modo
    `disabled` exige `errors = 0`.
  - Este control deja VE-D02 listo para un ensayo PE neutro que fije también a
    cero contaminación, quimeras e *index hopping*.

## Lista de permisos de Claude Code

El 4 de septiembre de 2026 se añadió `.claude/settings.json` con una lista de
permisos que evita las confirmaciones repetidas de las puertas Cargo. Quedan
dos concesiones que se aceptaron a sabiendas y un comportamiento del emparejado
que conviene tener presente.

- [ ] `Bash(cargo test *)` y `Bash(cargo build *)` están permitidos sin
  confirmación aunque no sean estrictamente de solo lectura: escriben en
  `target/` y los tests crean directorios temporales. Se admitieron porque son
  puertas obligatorias de `CONTRIBUTING.md` y sin ellas cada verificación
  interrumpe con una pregunta.
  - Revisar si conviene acotarlos a las formas exactas de `CONTRIBUTING.md`
    —`cargo test --workspace --all-targets --locked` y el `cargo build` de
    entrega— en lugar del comodín.
  - Relacionado: quedan residuos `tmp-*/` en la raíz del árbol de trabajo, como
    `tmp-genome-TDvk2K` del 10 de agosto de 2026. `.gitignore` los excluye con
    `tmp*/`, así que no contaminan el repositorio, pero ninguna prueba los
    limpia al terminar. Decidir si el destino temporal debe salir de la raíz o
    si su limpieza pertenece a la propia prueba.
- [ ] `Bash(cargo fmt --all -- --check)` se permitió como forma exacta y no como
  comodín, porque `cargo fmt` sin `--check` reescribe el código fuente. No
  ampliarlo a `Bash(cargo fmt *)`: concedería reescritura de fuentes sin
  confirmación.
  - `.claude/` no está bajo seguimiento de Git. Antes de versionarlo, decidir si
    la lista pertenece al repositorio compartido o a `settings.local.json`.
- [ ] Los patrones de la lista no cubren un comando compuesto: `samtools
  --version` coincide, pero `samtools --version 2>&1 | head -2` no, porque la
  tubería forma otra orden. La comprobación del 4 de septiembre de 2026 lo
  confirmó ejecutando ambas formas.
  - No ampliar los patrones con comodines para abarcarlas: una regla lo bastante
    ancha para admitir cualquier cosa después de la tubería concedería
    precisamente lo que no se ha revisado.
  - Registrado solo como comportamiento conocido; no requiere cambio.

## Rendimiento de la generación de lecturas (F3/F4)

Medido el 5 de septiembre de 2026 con el binario de `2da9ea5`, referencia
derivada de tomate de 500.000 bases en trece contigs, sobre cuatro núcleos.

| Individuos | Cobertura | Loci | Tiempo |
|---:|---:|---:|---:|
| 2 | 1x | 1000 | 25–29 s |
| 4 | 1x | 1000 | 137 s |
| 2 | 2x | 1000 | 102 s |
| 2 | 2x | 10 | 112 s |
| 2 | 4x | 1000 | 348 s |

- [x] **Resuelto.** El coste escalaba de forma casi cuadrática porque, al
  publicar un evento técnico, la segunda pasada reiniciaba el muestreador por
  cada pareja y recorría su serie desde el principio. Reproducirlas en una sola
  pasada dejó el escalado lineal: el ensayo de 20 individuos a 30x pasó de
  proyectar semanas a tardar 76 segundos. La corrección alcanza también las
  parejas que aporta una donante, que solo se ejercitan con contaminación
  activa y que la primera versión dejó sin tratar.
  - Medir solo el escenario neutro ocultó la mitad del problema: con los
    artefactos técnicos desactivados no se recorre la vía de la contaminación.
    Conviene medir cada cambio de rendimiento con al menos dos configuraciones,
    una neutra y otra con errores, quimeras, *index hopping* y contaminación.

## Memoria del demultiplexado y la escritura FASTQ

Medido el 8 de septiembre de 2026 con la referencia reducida de tomate de
500.000 bases, sobre cuatro núcleos.

| Ensayo | Antes | Después |
|---|---:|---:|
| 4 individuos, 16x, PE neutro | 202 MB | 100 MB |
| 4 individuos, 16x, con artefactos | 206 MB | 103 MB |
| 20 individuos, 30x, PE neutro | 1.615 MB | 468 MB |

- [x] **Resuelto.** La ejecución acumulaba todas las parejas de todas las
  bibliotecas antes de demultiplexar, de modo que el pico crecía con la
  cobertura: cada pareja retenida costaba unos 1,5 kB y el ensayo de 20
  individuos a 30x llegaba a 1.615 MB. Ahora los FASTQ R1 y R2 de todas las
  bibliotecas se mantienen abiertos a la vez y cada pareja se escribe en cuanto
  se conoce su biblioteca final, sin conservarla. El pico bajó un 71 % y el
  tiempo de CPU un 20 % (81,8 a 65,8 segundos). Los artefactos publicados
  quedan idénticos byte a byte, comprobado sobre los 173 MB y las veinte
  bibliotecas del ensayo grande, y también con errores de secuenciación,
  quimeras, *index hopping* y contaminación activos.
  - Lo que lo hacía posible ya estaba en el diseño: el destino de una pareja
    tras el *index hopping* depende solo de su biblioteca e instancia
    anteriores al salto, nunca de las demás parejas. Por eso el demultiplexado
    puede resolverse pareja a pareja sin cambiar ningún sorteo.
  - `ExecutionWorkspace` necesitó un artefacto en flujo con vida propia
    (`begin_streamed_artifact` / `finish_streamed_artifact`), porque la vía
    anterior abría, validaba y publicaba en una sola llamada y no admitía dos
    salidas simultáneas. Las garantías se conservan: el validador sigue
    recibiendo el temporal, nada entra en el catálogo antes de renombrarse y
    abandonar un flujo elimina su temporal.
  - Al variar solo la cobertura, el pico dejó de moverse: 96,2 / 96,3 / 99,7 MB
    para 4 individuos a 4x, 8x y 16x, donde antes eran 96,4 / 119,0 / 202,1 MB.

### Las claves PH que la primera pasada conserva

Medido el 8 de septiembre de 2026 con la misma referencia reducida de 500.000
bases, sobre cuatro núcleos. La línea de partida es `ab33d05`.

| Ensayo | Antes | Después | Δ |
|---|---:|---:|---:|
| 4 individuos, 16x, PE neutro | 99,7 MB | 91,8 MB | −8 % |
| 4 individuos, 16x, con artefactos | 102,8 MB | 91,7 MB | −11 % |
| 20 individuos, 8x, PE neutro | 190,3 MB | 161,0 MB | −15 % |
| 20 individuos, 15x, PE neutro | 279,2 MB | 229,6 MB | −18 % |
| 20 individuos, 30x, PE neutro | 467,7 MB | 374,7 MB | −20 % |

- [x] **Resuelto.** Una clave retenida ya no cuesta dos asignaciones de texto.
  El texto de todo identificador tipado vive ahora en una asignación compartida
  e inmutable (`Arc<str>` dentro de `typed_id!`), de modo que clonar una
  identidad es incrementar un contador: el individuo, que se repite en todas
  las claves de su biblioteca, se sostiene una sola vez, y el registro de
  colisiones técnicas de F4 anota el mismo texto que la clave ya sostiene en
  lugar de duplicarlo. El orden, la igualdad, el hash y el `Debug` siguen
  siendo los del texto, así que ningún sorteo ni identificador se mueve: los
  artefactos publicados quedan idénticos byte a byte (`diff -r`) en las cinco
  configuraciones, incluida la que activa errores de secuenciación, quimeras,
  *index hopping* y contaminación. La pendiente bajó de 12,6 a 9,7 MB por
  unidad de cobertura sobre un término fijo de unos 84 MB, y el tiempo mejoró
  también: el ensayo de 30x pasó de 69,8 a 60,5 segundos y las llamadas al
  asignador cayeron un 27 % (118,6 a 86,2 millones en el ensayo de 8x).
  - Reparto del pico según `heaptrack` con símbolos, ensayo de 20 individuos a
    8x, en megabytes de montículo pedido:

    | Consumidor | Antes | Después |
    |---|---:|---:|
    | Proyección de verdad (no depende de la cobertura) | 33,7 | 33,7 |
    | Vector de claves limpias de las 20 bibliotecas | 15,0 | 10,0 |
    | Texto de las instancias `PH_` (clave, registro y derivación) | 15,5 | 12,8 |
    | Texto del individuo repetido en cada clave | 1,6 | 0,0 |
    | Nodos del registro de colisiones de instancias `PH_` | 12,9 | 9,9 |
    | Registro de colisiones de parejas técnicas `TP_` | 19,3 | 21,0 |
    | Parejas regeneradas de la biblioteca en curso | 11,9 | 10,9 |
    | **Pico de montículo** | **145,2** | **129,5** |
    | **RSS máximo** | **190,3** | **161,0** |

    El RSS baja más que el montículo pedido porque desaparecen unas 530.000
    asignaciones vivas, y cada una costaba además su cabecera y su redondeo de
    bloque.

  - La parte variable de `PreHopInstanceId` **no** es un ordinal derivable en
    la práctica: es base32 de los 16 primeros bytes de un BLAKE3 sobre
    (semilla, evento, biblioteca, origen, plantilla). Se deriva del ordinal de
    plantilla, sí, pero el orden de las claves es el del texto y ese orden
    alimenta los sorteos de donante y de retirada; el alfabeto base32 no
    conserva el orden de los bytes (`2` vale 50 y `A` vale 65), de modo que
    ordenar por ordinal o por el dígito crudo daría otra permutación y otros
    FASTQ. Por eso se descartó guardar solo el ordinal, y por eso una
    representación compacta de longitud fija tampoco compensa: el texto tiene
    que seguir estando para comparar.
  - No hizo falta la vía invasiva —formar el plan de contaminación sin tener a
    la vez las claves de todas las bibliotecas—. La representación compartida
    bastó para bajar un 20 % sin tocar ninguna decisión técnica.
  - Lo que no funcionó: `shrink_to_fit()` sobre el inventario de cada
    biblioteca. La holgura del crecimiento duplicado existe (unos 1,5 MB de los
    10 del vector a 8x), pero no se traduce en RSS, porque el asignador
    reutiliza los búferes que liberó la biblioteca anterior. El pico de 30x
    solo pasó de 374,7 a 373,6 MB, dentro del ruido, así que se descartó por no
    pagar una copia por biblioteca a cambio de nada.
  - Coste conocido del cambio: para un texto que nadie más sostiene, `Arc<str>`
    gasta 16 bytes de cabecera de contadores en el montículo, así que el
    registro de parejas técnicas subió unos 8 bytes por pareja (19,3 a
    21,0 MB). Lo compensa con mucho el ahorro de las `PH_`, que antes se
    copiaban tres veces y ahora las sostienen a la vez la clave, el registro y
    el conjunto de la regeneración.

- [ ] El pico sigue creciendo con la cobertura, ahora a 9,7 MB por unidad, pero
  las claves PH ya no son el consumidor principal. Con el reparto de arriba, el
  siguiente candidato es el registro de colisiones de F4.5: conserva el texto
  de cada `TechnicalPairId` de todo el evento (21,0 MB a 8x, del orden de
  80 MB a 30x) y ninguna otra estructura lo sostiene, de modo que no hay nada
  que compartir. Tocarlo es tocar una garantía de identidad (DTC-039), así que
  conviene tratarlo como trabajo propio y con su evidencia.
  - El otro término que crece es `RegeneratedPairs`: retiene las parejas
    transformadas de una biblioteca entera, unos 800 bytes por pareja, porque
    las claves llegan en orden de identificador y no en el de generación.
    Emitirlas en el orden de la serie exigiría cambiar el orden de escritura de
    los FASTQ, que es parte de la identidad binaria.
  - El plan de contaminación, además, clona la lista completa de candidatas de
    todas las demás bibliotecas por cada receptora. Con el texto compartido ya
    no asigna cadenas, pero sigue copiando 32 bytes por clave y por receptora,
    lo que es cuadrático en número de individuos. Ordenar referencias en lugar
    de clones lo evitaría sin cambiar la permutación, porque el criterio de
    orden es el mismo.

## Memoria de la referencia

Medido el 5 de septiembre de 2026 sobre el tomate de 759 MB y sobre genomas
simulados, con el binario de `d950aca`.

| Ruta | Tiempo | Memoria | ¿Crece con el genoma? |
|---|---:|---:|---|
| Suministrada sin reducir ni matriz SNV | 12 s | 5,6 MB | no |
| Suministrada reducida | 16 s | 8,3 MB | no |
| Generación al azar de 50 Mb | 4,6 s | 152 MB | sí, unas 3 veces |
| Suministrada con matriz SNV | 76 s | 1.544 MB | sí, unas 2 veces |

Las dos primeras quedaron resueltas: validar sin conservar las bases y pedir los
tramos por índice FAI deja su memoria independiente del genoma de partida.

Las dos últimas se redujeron a la mitad el 8 de septiembre de 2026 al dejar de
copiar bases que ya estaban en memoria: la matriz SNV las duplicaba dos veces y
la generación al azar una. La ruta con matriz SNV pasó de 1.544 a 780 MB con el
tomate, y generar 50 Mb de 152 a 103 MB.

- [ ] Cada una conserva todavía una copia del genoma, y esa sí la necesita su
  comprobación: la matriz se valida contra las bases, y una referencia simulada
  se compara con la que la semilla regenera. Bajar de ahí exige comparar
  secuencia a secuencia en flujo, en lugar de tener las dos versiones completas
  a la vez.
  - No es una optimización más: cambia cómo se comprueba que lo publicado
    coincide con la semilla, que es una garantía de reproducibilidad. Conviene
    tratarlo como trabajo propio, con su evidencia.
  - El pico quedaría en una secuencia en lugar de un genoma, de modo que un
    cromosoma humano grande costaría unos 250 MB en vez de varios gigabytes.
  - `heaptrack` es la única herramienta de perfilado disponible en el entorno:
    no hay `perf`, `valgrind` ni `gdb`. Basta para un cuello de memoria; para
    uno de CPU haría falta instalar algo.
