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

- [ ] Las dos últimas siguen creciendo con el genoma. El pico son dos lecturas
  del FASTA que conviven dentro de la validación física de la publicación, cada
  una con el genoma entero; con 50 Mb son exactamente dos bloques de 50 MB.
  - Ambas arrancan en `validate_published_resources`, pero la segunda se pierde
    en un `collect` de la biblioteca estándar y localizarla exige más perfilado.
  - Un genoma humano proyecta del orden de 6 GB en la ruta con matriz SNV y
    unos 9 GB generando al azar.
  - `heaptrack` es la única herramienta de perfilado disponible en el entorno:
    no hay `perf`, `valgrind` ni `gdb`. Basta para un cuello de memoria; para
    uno de CPU haría falta instalar algo.
