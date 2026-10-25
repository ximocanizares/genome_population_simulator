# Simulaciones con un agente: guía de uso

El simulador genera, a partir de una descripción en un fichero TOML, una
población de individuos diploides con su verdad genotípica, un VCF y, si se
pide, lecturas Illumina pareadas. Esta guía explica cómo instalarlo y cómo
manejarlo a través de un agente de IA (Claude Code, por ejemplo): tú pides la
simulación en lenguaje natural y el agente escribe la configuración, la valida,
la ejecuta y analiza el resultado.

Está pensada para estudiar la **endogamia (F)** y los **microsatélites (SSR)**,
con y sin PCR.

## Qué hace y qué no hace todavía

Hace:

- Genera individuos diploides en equilibrio de Hardy-Weinberg o con un
  coeficiente de endogamia F, sobre una referencia FASTA que aportas (o que
  simula el propio programa).
- Coloca SNV con frecuencias que declaras, y microsatélites cuyos alelos
  declaras o se generan por deslizamiento.
- Escribe la verdad (cada copia y cada genotipo) y un VCF.
- Simula lecturas Illumina pareadas de hasta 150 pb con errores, con o sin PCR
  (duplicados, familias de PCR) y con tartamudeo en los microsatélites.

No hace: BAM o CRAM, variantes estructurales, evolución por generaciones,
cruces, lecturas largas ni genomas grandes. Las limitaciones que afectan a F y a
los microsatélites están en la sección final.

## Instalación (Linux o WSL)

Necesitas `git`, `curl` y un compilador de C. El resto se instala así.

1. **Rust.** El repositorio fija la versión 1.97.1 en `rust-toolchain.toml`
   (con `rustfmt` y `clippy`); `rustup` la instala sola al compilar.

   ```bash
   curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
   ```

2. **Compilar el binario optimizado**, desde el repositorio del simulador:

   ```bash
   cargo build --release --locked --target x86_64-unknown-linux-gnu
   ```

   El binario queda en
   `target/x86_64-unknown-linux-gnu/release/genome-population-simulator`.
   Usa siempre este, no el de depuración: 100 individuos a 30X tardan unos 3 a
   5 s, frente a 137 s.

3. **Cuaderno marimo (opcional).** Necesita [`uv`](https://docs.astral.sh/uv/)
   y Python 3.12:

   ```bash
   uv sync --all-groups
   uv run --frozen python scripts/run_marimo.py      # http://127.0.0.1:2718
   ```

   `run_marimo.py` ejecuta `marimo run ui/gps_app.py` en el puerto 2718 (`--port`
   para cambiarlo, `--no-open` para no abrir el navegador). Abrir el cuaderno no
   lanza ninguna simulación. Con el agente no hace falta, pero sirve para
   revisar configuraciones a mano.

## Trabajar con el agente

**No trabajes dentro del repositorio del simulador.** Su `CLAUDE.md` es política
de *desarrollo* (ciclos de implementación y revisión, commits…) y confunde a un
agente que solo quiere usar el programa. Haz esto:

1. Crea un directorio de trabajo **fuera** del repositorio, por ejemplo
   `~/simulaciones`.
2. Copia en él la plantilla:
   `cp <repositorio>/uso-con-agentes/CLAUDE.md ~/simulaciones/CLAUDE.md`
3. Arranca el agente desde ese directorio (`cd ~/simulaciones && claude`). Le
   dirás dónde está el repositorio y el binario; la plantilla le prohíbe
   modificarlo.

La plantilla enseña al agente el ciclo (escribir el TOML, `validate`, corregir
según el diagnóstico, `run`, analizar), dónde está la referencia de campos
(`schemas/config-general-v4.schema.json`), las trampas conocidas y cómo leer los
resultados.

### Qué pedirle

Pide en términos biológicos y di lo que quieres medir:

- «Simula 50 individuos con F = 0,3 y un microsatélite (AG)15 a 30X con PCR.»
- «Lo mismo sin PCR, con la misma semilla, para comparar.»
- «Calcula Ho, He y F_IS por locus y dime si F_IS sale cerca de 0,3.»
- «Cuenta las unidades del microsatélite en las lecturas del individuo I001 y
  compáralas con su genotipo real.»

Pídele una **semilla fija** si quieres repetir o comparar: sin ella, cada
ejecución sale distinta. Los tres ejemplos de `ejemplos/` sirven de punto de
partida; el agente puede copiarlos.

### Cómo saber que lo que devuelve tiene sentido

- **Pide el TOML.** Revisa `individual_count`, `inbreeding_coefficient`, el
  motivo y las posiciones del microsatélite, la cobertura y la tabla de PCR. El
  agente debe guardar una copia en `sim-.../general.toml`.
- **Declarado frente a realizado.** `manifest.json` y
  `events/<evento>/manifest.json` dicen lo que se hizo de verdad: fracción de
  duplicados, familias, deslizamientos aplicados, frecuencias alélicas
  realizadas. Deben cuadrar con lo pedido, salvo ruido de muestreo.
- **Con muestras pequeñas, el ruido es grande.** Con 20 individuos, una F_IS
  puede variar mucho de un SNV a otro (en `f-ssr`, con 200 individuos y F = 0,5,
  entre 0,47 y 0,64). Pide varios loci o más individuos antes de concluir.
- **Pregúntale qué no puede concluir.** La plantilla le pide avisar de las
  limitaciones; si no lo hace, pregúntalo.
- **Ejecuciones grandes.** El agente debe avisarte antes de lanzarlas. Orden de
  magnitud: 100 individuos, 2000 pb y 30X, unos 3 a 5 s y 41 MB.

## Qué contiene una carpeta de resultados

| Fichero | Qué es |
|---|---|
| `outputs/truth.jsonl` | La verdad: qué lleva cada copia de cada individuo, y cada genotipo. Coordenadas desde cero. |
| `outputs/truth.vcf` | La misma verdad como VCF (coordenadas desde uno). Un microsatélite es un registro con varios alelos y las unidades de repetición en `RUC`. |
| `outputs/fastq/<evento>/<individuo>_R1/_R2.fastq.gz` | Lecturas pareadas, solo si hay evento técnico. |
| `effective-config.json` | La configuración como se aplicó, con los alelos generados, sus frecuencias, la calibración del tartamudeo y la semilla usada. |
| `manifest.json` | Resumen de la ejecución con el SHA-256 de cada fichero y los recuentos realizados por locus. |
| `events/<evento>/manifest.json` | Lo realizado en las lecturas: duplicados, familias de PCR, deslizamientos y tartamudeos, cobertura. |
| `resources/references/<id>.fa` | La referencia usada. |
| `inputs/config/general.toml` | Copia de la configuración. |
| `population_snvs.tsv`, `population_snvs.metadata.json` | La matriz de SNV usada y su descripción. |

## Limitaciones que debes conocer

**Cómo funciona F.** F se aplica **locus a locus y de forma independiente**. Cada
locus tiene el déficit de heterocigotos correcto, pero todos los individuos
tienen la misma endogamia esperada, no hay correlación entre loci y no hay
tramos de homocigosis. Sirve para estudiar F a escala de población (F_IS,
frecuencias genotípicas), no para individuos más o menos consanguíneos ni para
estimar la F individual. La endogamia por genealogía llegará con los cruces de
una fase futura.

**F también se aplica a los microsatélites,** con la misma fórmula multialélica
(desde el 25/10/2026; antes los ignoraba). En `ejemplos/f-ssr`, con F = 0,5, el
microsatélite da F_IS = 0,49. Si sale cerca de 0, el binario es anterior al
arreglo: recompílalo.

**El tartamudeo heredado por familias de PCR es raro por diseño.** Las copias de
una familia divergen en los primeros ciclos. Para verlo hacen falta homopolímeros,
pocos ciclos o mucha cobertura. No verlo con un dinucleótido a 30X no es un fallo.

**La cobertura pedida no es la muestra útil sobre un microsatélite.** Solo
cuentan las lecturas que lo cruzan con flanco. En `ssr-pcr`, un individuo tiene
400 lecturas y solo 28 cruzan el tracto (17 sin PCR). Con PCR, además, muchas son
duplicados.

**Lo que no hay:** BAM o CRAM; variantes estructurales; evolución por
generaciones; cruces; lecturas largas (solo Illumina pareadas de hasta 150 pb);
escala para genomas grandes (no probada).

Otras cosas a tener presentes:

- Los alelos generados por deslizamiento salen de linajes independientes que
  parten del alelo de la referencia (una genealogía en estrella, sin
  coalescencia): sus frecuencias no reflejan parentesco ni estructura entre
  alelos. Después, cada copia sortea su alelo de ese abanico.
- Las calibraciones del tartamudeo y de la PCR son aproximadas; el manifiesto las
  marca (`coarse`) y cita sus fuentes.
- Los ejemplos usan una referencia FASTA aportada, con el tracto del
  microsatélite ya escrito entre `start` y `end`; los SNV deben quedar fuera de
  él. El programa también puede simular la referencia (`kind = "simulated"`),
  pero esa vía no se ha probado en esta guía.

## Ejemplos

En `ejemplos/`, cada uno con su `general.toml`, su referencia, su matriz de SNV y
un `README.md`. Todos se han validado y ejecutado con el binario optimizado:

| Carpeta | Qué es |
|---|---|
| `ssr-pcr` | 20 individuos, microsatélite (CA)24, 30X con PCR y tartamudeo. |
| `ssr-sin-pcr` | El mismo, sin PCR, con la misma semilla. |
| `f-ssr` | F = 0,5, cuatro SNV y un microsatélite; solo verdad y VCF. |

Desde la carpeta del ejemplo: `gps validate general.toml` y
`gps run general.toml --output <carpeta-nueva>`, con `gps` el binario optimizado.
