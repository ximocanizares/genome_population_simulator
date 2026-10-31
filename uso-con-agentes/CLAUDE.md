# Instrucciones para el agente: simulaciones con genome-population-simulator

Eres el agente de una persona que **usa** el simulador de poblaciones genómicas
para generar datos sintéticos (verdad genotípica, VCF y lecturas Illumina). Tu
trabajo es traducir lo que pide en una configuración TOML, validarla, ejecutarla
y analizar el resultado. Responde en castellano.

Sustituye `SIMULADOR` por la ruta del repositorio del simulador. Esa ruta, y la
del binario, se las pides a la persona la primera vez y las anotas aquí:

- Repositorio: `SIMULADOR = ...`
- Binario: `$SIMULADOR/target/x86_64-unknown-linux-gnu/release/genome-population-simulator`
  (alias `gps` en lo que sigue). Si no existe, dile a la persona que lo compile
  con `cargo build --release --locked --target x86_64-unknown-linux-gnu` en el
  repositorio; no lo compiles tú.

## Reglas inviolables

1. **No modifiques el repositorio del simulador.** Solo lees (esquemas, ejemplos)
   y ejecutas el binario. Nada de editar, compilar, `git`, `cargo`, ni borrar
   ficheros allí.
2. **No trabajes dentro del repositorio.** Su `CLAUDE.md` es política de
   desarrollo y no te corresponde. Tu directorio de trabajo es este, fuera de él.
3. **Una simulación = una carpeta** de este directorio con su `general.toml`,
   sus entradas (FASTA, matriz de SNV) y su salida. No sobrescribas la de otra.
4. **Avisa a la persona antes de lanzar una ejecución grande** (ver «Rendimiento»)
   y espera su confirmación.
5. **No inventes campos.** Si dudas de uno, consulta el esquema (abajo). Una
   afirmación sobre el programa que no hayas comprobado ejecutándolo, di que no
   la has comprobado.

## El ciclo de trabajo

1. Entiende la petición y fija con la persona lo que falte: individuos, F,
   locus y motivo, cobertura, PCR sí o no, semilla. Si pide reproducibilidad,
   fija `seed` (64 caracteres hexadecimales); sin `seed` cada ejecución sale
   distinta (comprobado).
2. Crea `sim-NOMBRE/` con `general.toml`, `reference.fa` y la matriz de SNV.
   Parte de un ejemplo del repositorio, que no se toca, se copia:
   `$SIMULADOR/uso-con-agentes/ejemplos/` (`ssr-pcr`, `ssr-sin-pcr`,
   `f-ssr`).
3. Valida: `gps validate sim-NOMBRE/general.toml`. Sin salida y código 0 es
   válido. Con error, sale `GPS-xxx-nnn` y un `detail` que dice qué campo falla;
   `--diagnostics json` da lo mismo en JSON. Corrige y repite hasta que pase.
   Un error de validación sale con código 2.
4. Ejecuta: `gps run sim-NOMBRE/general.toml --output sim-NOMBRE/salida`. El
   directorio `salida` **no debe existir** (si existe, falla con `GPS-IO-003`).
   Para repetir, usa otro nombre o bórralo tú, solo si es tuyo.
5. Analiza (secciones siguientes) y **compara con lo que declaraste**: tasas
   efectivas, recuentos realizados, frecuencias. Informa a la persona de lo
   que ves y de lo que no se puede concluir.

## Referencia de la configuración

- **Esquema, con descripciones de cada campo:**
  `$SIMULADOR/schemas/config-general-v4.schema.json`. Cada propiedad lleva su
  `description` y anotaciones `x-gps-tab`, `x-gps-control` y `x-gps-unit`
  (unidad: `fraction`, `percent`, `bases`). Léelo con Python en vez de volcarlo
  entero: tiene decenas de miles de caracteres. Las definiciones útiles están
  en `$defs`: `repeat_locus`, `repeat_allele`, `repeat_slippage`,
  `initial_population`, `technical_event` (con `pcr` y `stutter`),
  `population_snvs`, `supplied_reference`.

```bash
python3 - <<'E'
import json
s = json.load(open("SIMULADOR/schemas/config-general-v4.schema.json"))
d = s["$defs"]["repeat_locus"]
for k, v in d["properties"].items():
    print(k, "-", v.get("description", "")[:300])
E
```
- **Especificaciones** (`$SIMULADOR/design/especificaciones-0.1.md`), para el
  porqué: 4.2 referencias y microsatélites, 4.3 individuos y F, 4.4 datos por
  muestra (PCR, tartamudeo, lecturas), 4.5 VCF de verdad, 4.7 manifiestos. Son
  largas; busca por palabra clave.
- Los ejemplos validados están en `$SIMULADOR/uso-con-agentes/ejemplos/`.

### Estructura mínima de una configuración

```toml
config_kind = "general"
schema_version = "4"          # obligatorio con microsatélites
seed = "<64 hex>"             # opcional

[population_snvs]             # matriz TSV: contig position ref alt population_id alt_af
mode = "supplied"
input_path = "population-snvs.tsv"

[reference]
kind = "supplied"
id = "ref-ssr"
fasta = "reference.fa"
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 2000

[[reference.repeat_loci]]     # el tramo start..end de la referencia ES el tracto
id = "STR-CA"
sequence = "chr1"
start = 1001                  # base uno, inclusiva
end = 1048
motif = "CA"
max_units = 60
[reference.repeat_loci.slippage]   # alelos generados; excluyente con [[...alleles]]
model = "stepwise"
mean_events = 3.0
expansion_bias = 0.5

[reference.repeat_detection]       # OBLIGATORIA con microsatélites: activa la detección
[reference.repeat_detection.slippage]   # deslizamiento común; model y mean_events obligatorios
model = "stepwise"
mean_events = 1.0
# opcionales: min_copies = [8, 6, 6, 6, 5, 4] (motivos de 1 a 6 pb), max_units = 200,
# y [[reference.repeat_detection.regions]] (sequence, start, end) para limitar la detección

[[initial_populations]]
population_id = "P1"
individual_count = 20
individual_ids = ["I001", ...]     # exactamente individual_count
ploidy = 2
genotype_model = "hwe"             # o "f_diploid" + inbreeding_coefficient

[vcf_output]
mode = "text"

[technical_event]                  # lecturas Illumina; omítelo si solo quieres verdad y VCF
event_id = "reads-ssr"
individual_ids = [...]
target_coverage_wgs = 30
library_preparation = "pcr"        # o "pcr_free"
[technical_event.pcr]              # solo con "pcr"; los tres campos son obligatorios
duplicate_fraction = 0.7
cycles = 10
efficiency = 0.9
[technical_event.stutter]
enabled = true                     # declarar un microsatélite NO enciende el tartamudeo
# cuatro [[technical_event.stutter.steps]] (units -3, -2, -1, 1, cada uno con weight)
# y una [[technical_event.libraries]] por individuo (individual_id, i7, i5)
```

**Con la sección de microsatélites el simulador barre la referencia y detecta
todos sus microsatélites de 1 a 6 pb**; cada uno varía con el deslizamiento
común. Un locus declarado solo ajusta uno detectado: tiene que solaparse con un
microsatélite detectado de su motivo (o una rotación), y entonces mandan sus
fronteras y sus alelos. Si la persona quiere solo los loci que declara, limita
la detección con regiones que los cubran. En los tres ejemplos el único
microsatélite que se detecta es el declarado (comprobado).

Si pide alelos explícitos en vez de deslizamiento, usa `[[reference.repeat_loci.alleles]]`
con `sequence` (o `units`) y `frequency` (suman 1); no los mezcles con `slippage`.
Con `ploidy = 2` y `genotype_model = "f_diploid"`, F va en
`inbreeding_coefficient` (0 a 1).

## Trampas conocidas (todas comprobadas ejecutando)

| Trampa | Qué ocurre |
|---|---|
| Microsatélite sin `[reference.repeat_detection]` | `los loci de repetición, los tractos insertados y los plantados solo existen con la sección de microsatélites: declare [reference.repeat_detection] ...`. Añade la sección con su `slippage`. |
| Locus declarado donde la referencia no tiene ese microsatélite | `locus ...: no coincide con ningún microsatélite detectado en la referencia final (...)`, con el detectado más cercano y el criterio que no alcanzó el tramo. Comprueba posiciones y motivo. |
| Microsatélite con `schema_version = "3"` (u otra) | Error `GPS-CONF-101` que dice «campo desconocido: initial_populations». El mensaje confunde: la causa es la versión. Pon `"4"`. |
| Tartamudeo sin los cuatro pasos | `hay que declarar el peso de los cuatro pasos admitidos [-3, -2, -1, 1]; faltan -3`. Declara los cuatro (pueden tener peso 0). |
| Declarar el microsatélite y esperar tartamudeo | No se enciende solo: hace falta `[technical_event.stutter] enabled = true` (especificaciones, 4.4). |
| `f_diploid` con ploidía distinta de 2 | `f_diploid solo se admite con ploidía 2` (probado con 1 y 3). |
| `f_diploid` sin `inbreeding_coefficient` | `f_diploid exige inbreeding_coefficient`. |
| `hwe` con `inbreeding_coefficient` | `hwe no admite inbreeding_coefficient`. |
| `individual_ids` de otra longitud que `individual_count` | Rechazado. |
| Individuo del evento técnico que no está en la población | `individuo no declarado en initial_populations`. |
| Bibliotecas: faltan todas, o hay i7 o i5 repetidos | `falta un campo obligatorio` o `cada componente UDI i7 e i5 debe ser único en la tanda`. Hace falta una biblioteca por individuo, con i7 e i5 distintos en toda la tanda (los ejemplos usan 8 bases). **Genéralas con un script** (abajo). |
| `library_preparation = "pcr"` sin tabla `[technical_event.pcr]` | Rechazado: los tres campos son obligatorios, sin valores por omisión. Y al revés: una tabla `[technical_event.pcr]` con `library_preparation = "pcr_free"` también se rechaza. |
| Directorio de salida existente | `GPS-IO-003`. El directorio de `--output` no debe existir. |
| `--output` omitido en `run` | `run exige --output <ruta>`. |
| Rutas del TOML (`fasta`, `input_path`) | Son **relativas al directorio del TOML**, no al directorio actual (comprobado lanzando desde otro directorio y con un decoy en el actual). `--output` sí es relativo al directorio actual. |
| SNV de la matriz dentro de un microsatélite (incluidos sus extremos) | `la variante ... cae dentro del tracto declarado`. Mueve el SNV fuera del tramo `start..end`. |
| SNV cuya base `ref` no coincide con el FASTA | `ReferenceMismatch`. La columna `ref` debe ser la base de la referencia en esa posición (base uno). |
| Matriz de SNV vacía (solo la cabecera) | `EmptyActiveMatrixRequiresAuthorization`: hay que añadir `--allow-empty-snv-matrix`. |
| Alelos explícitos: qué se rechaza | Pureza menor de 0,8 respecto al motivo (por ejemplo, otro motivo) y unidades por encima de `max_units`. **No** se rechaza un alelo más largo que el tracto de la referencia (probado con 60 pb sobre un tracto de 48). Los umbrales se cambian con `min_purity` y `min_longest_run_units`. |
| `slippage` y `alleles` juntos | Excluyentes: o uno u otro. |
| `[reference.repeat_loci]` fuera de la secuencia | `el intervalo ... no cabe en chr1`. |
| Campo mal escrito | `GPS-CONF-006 campo desconocido` con línea y columna. |
| **F y microsatélites** | F se aplica también a los microsatélites (comprobado: con F = 0,5, `ejemplos/f-ssr` da F_IS 0,49 en el microsatélite). Si sale cerca de 0, el binario es anterior al 25/10/2026: pide a la persona que lo recompile. |

El programa solo admite una referencia suministrada (FASTA) o simulada
(`kind = "simulated"`); los ejemplos usan la suministrada, que es la que se ha
comprobado. Para un microsatélite en una referencia suministrada, el tracto
debe estar ya en el FASTA entre `start` y `end` (o declararse con
`insert_after` y `reference_units`, no comprobado aquí).

### Generar individuos y bibliotecas

```python
import random
def bloques(n, semilla=1, prefijo="I"):
    """Devuelve (ids, texto TOML de las bibliotecas) con i7 e i5 únicos."""
    rng, vistos, libs = random.Random(semilla), set(), []
    def nuevo():
        while True:
            s = "".join(rng.choice("ACGT") for _ in range(8))
            if s not in vistos:
                vistos.add(s); return s
    ids = [f"{prefijo}{i:03d}" for i in range(1, n + 1)]
    for i in ids:
        libs.append(f'[[technical_event.libraries]]\nindividual_id = "{i}"\ni7 = "{nuevo()}"\ni5 = "{nuevo()}"\n')
    return ids, "\n".join(libs)
lista = lambda ids: "[" + ", ".join(f'"{x}"' for x in ids) + "]"
```

Úsalo para escribir `individual_ids` en `[[initial_populations]]` y en
`[technical_event]`, y añade el texto de las bibliotecas al final del TOML.

## Rendimiento orientativo (medido)

Con el binario de `release`: 100 individuos, 2000 pb y 30X, con PCR y
tartamudeo, tardan unos 3 a 5 s y 41 MB (4,85 s medidos aquí). Con el binario de
depuración, 137 s. 20 individuos con el mismo escenario, 1 s y 17 MB. Usa
siempre el de `release`. El coste crece con individuos × longitud × cobertura:
un genoma grande o cientos de individuos a 30X **no se han probado**; avisa a
la persona antes de lanzar algo así y empieza por una versión reducida.

## Qué hay en el directorio de salida

| Fichero | Contenido |
|---|---|
| `outputs/truth.jsonl` | Verdad por copia haplotípica y por genotipo. Coordenadas **base cero, semiabiertas**. |
| `outputs/truth.vcf` | VCF de la verdad (si hay `[vcf_output]`), coordenadas base uno. Los SNV son biálelicos; un microsatélite es un registro multialélico con `RUC` (unidades por ALT), `RUS`, `RUL` y `RB`. |
| `effective-config.json` | Configuración completa tal como se aplicó: alelos generados con su frecuencia efectiva, parámetros del tartamudeo y de la PCR (con la fuente de calibración), semilla de la ejecución (`seed.root`). |
| `manifest.json` | Resumen: artefactos con su SHA-256, estado, recuentos realizados por locus (`repeat_loci.loci[].alleles[].realized_copies`, `slippage.applied_events`) y frecuencias genotípicas realizadas de la matriz de SNV. |
| `events/<evento>/manifest.json` | Lo realizado en las lecturas: `artifacts.pcr` (`duplicate_pairs`, `families`, `realized_duplicate_fraction`, `family_size_histogram`, `lineage_slips_by_locus`, `lineage_inherited_slips_by_locus`), `artifacts.stutter` (`evaluated_tracts`, `stuttered_tracts`, `steps_by_locus`), `coverage`, `final_fastq_summary` (parejas y bases por biblioteca). |
| `outputs/fastq/<evento>/<id>_R1.fastq.gz`, `_R2` | Lecturas pareadas por individuo. |
| `resources/references/<id>.fa` | Referencia final usada. Úsala para los flancos. |
| `inputs/config/general.toml` | Copia de la configuración. |

Antes de interpretar nada, mira `manifest.json` (`state` debe ser `complete`) y
`diagnostics` (vacío). Compara lo **declarado** con lo **realizado**: la fracción
de duplicados realizada (`realized_duplicate_fraction`) cae cerca de la
declarada pero no es igual; las frecuencias alélicas realizadas varían por
muestreo con pocos individuos.

### Genotipos, Ho y He desde el VCF

```python
import sys
from collections import Counter

def leer_vcf(ruta):
    muestras, filas = [], []
    for linea in open(ruta):
        if linea.startswith("##"):
            continue
        f = linea.rstrip("\n").split("\t")
        if linea.startswith("#CHROM"):
            muestras = f[9:]
            continue
        info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
        alelos = [f[3]] + f[4].split(",")
        unidades = None   # solo en microsatélites: REF por longitud, ALT por RUC
        if "RUC" in info:
            motivo = info["RUS"].split(",")[0]
            unidades = [len(f[3]) / len(motivo)] + [float(x) for x in info["RUC"].split(",")]
        gts = {s: tuple(int(a) for a in g.split(":")[0].replace("/", "|").split("|"))
               for s, g in zip(muestras, f[9:])}
        filas.append(dict(chrom=f[0], pos=int(f[1]), alelos=alelos, unidades=unidades, gts=gts))
    return filas

def ho_he(fila):
    gts = list(fila["gts"].values())
    ho = sum(len(set(g)) > 1 for g in gts) / len(gts)
    c = Counter(a for g in gts for a in g)
    n = sum(c.values())
    he = 1 - sum((k / n) ** 2 for k in c.values())
    return ho, he

for fila in leer_vcf(sys.argv[1]):
    ho, he = ho_he(fila)
    fis = 1 - ho / he if he > 0 else float("nan")
    print(f'{fila["chrom"]}:{fila["pos"]} {"SSR" if fila["unidades"] else "SNV"} '
          f'alelos={len(fila["alelos"])} Ho={ho:.3f} He={he:.3f} F_IS={fis:.3f}')
```

`F_IS = 1 - Ho/He`. Con F de la configuración, espera F_IS cercana a F en cada
SNV, con ruido de muestreo (probado: F = 0,5, 200 individuos, SNV entre 0,47 y
0,64). Para el genotipo en unidades de un individuo:
`[fila["unidades"][a] for a in fila["gts"]["I001"]]`.

Desde `truth.jsonl` se puede ver lo mismo por copia: cada línea
`record_type = "haplotype_copy"` trae `identity.individual_id`, `identity.copy_index`
y `current_pieces` (la pieza con `kind = "repeat_allele"` lleva `locus_id` y
`blocks[].units`; una copia sin esa pieza lleva el alelo de la referencia); y cada `record_type = "genotype"` con `case_kind = "repeat"`
trae `alleles`, `allele_dosages` y `carriers`.

### Contar unidades del microsatélite en las lecturas

```python
import gzip, sys
from collections import Counter

def rc(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))

def cargar_fa(ruta, nombre):
    seq, dentro = [], False
    for l in open(ruta):
        if l.startswith(">"):
            dentro = l[1:].split()[0] == nombre
        elif dentro:
            seq.append(l.strip())
    return "".join(seq)

# uso: python3 cuenta.py referencia.fa chr1 1001 1048 CA R1.fastq.gz R2.fastq.gz
fa, nombre, ini, fin, motivo, *fqs = sys.argv[1:]
ini, fin = int(ini), int(fin)
ref = cargar_fa(fa, nombre)
izq, der = ref[ini - 13:ini - 1], ref[fin:fin + 12]      # flancos de 12 pb
cuenta, total = Counter(), 0
for fq in fqs:
    with gzip.open(fq, "rt") as h:
        for i, l in enumerate(h):
            if i % 4 != 1:
                continue
            total += 1
            for s in (l.strip(), rc(l.strip())):
                a = s.find(izq)
                b = s.find(der, a + len(izq)) if a >= 0 else -1
                if a >= 0 and b >= 0:
                    cuenta[(b - a - len(izq)) / len(motivo)] += 1
                    break
print(f"lecturas={total} cruzan={sum(cuenta.values())}")
for u in sorted(cuenta):
    print(f"{u:g}\t{cuenta[u]}")
```

Solo cuenta las lecturas que contienen entero el flanco izquierdo y el derecho
(los flancos con una variante o un error no casan: es una cota inferior). Para
un individuo, los valores esperados son los de su genotipo en el VCF más algún
paso de tartamudeo (-1, -2, -3 o +1 unidades).

**Muestra útil.** Pedir 30X no significa 30 lecturas sobre el microsatélite:
solo cuentan las que lo cruzan con flanco. Medido: 20 individuos, 30X, 2000 pb,
tracto de 48 pb. Un individuo tiene 400 lecturas y solo 28 cruzan el tracto con
PCR y 17 sin PCR. Con PCR, además, muchas son duplicados (≈70 % declarado). Si
la persona quiere llamar genotipos sobre un microsatélite, calcula con ella la
cobertura efectiva y propón subir la cobertura.

**Tartamudeo.** Hay dos efectos distintos:
- Sin PCR (`pcr_free`): errores de la polimerasa independientes por lectura
  (`artifacts.stutter`). Se ven como lecturas sueltas con ±1 a ±3 unidades.
- Con PCR: el deslizamiento ocurre en un ciclo y lo **heredan todas las copias
  de esa familia** (`artifacts.pcr.lineage_slips_by_locus`). Por diseño es
  raro: las copias de una familia divergen en los primeros ciclos. Para verlo
  hacen falta homopolímeros (motivo de 1 base), pocos ciclos o mucha cobertura.
  No lo des por fallo si no aparece con un dinucleótido a 30X.

## Limitaciones: avísalas a la persona cuando afecten a su pregunta

- **F se aplica locus a locus e independientemente.** Cada locus tiene el
  déficit de heterocigotos correcto, pero todos los individuos tienen la misma
  endogamia esperada y no hay correlación entre loci ni tramos de homocigosis.
  Sirve para F a escala de población (F_IS, frecuencias genotípicas), no para
  individuos más o menos consanguíneos ni para estimar F individual. La
  endogamia por genealogía llegará con los cruces de una fase futura.
- **No hay:** BAM o CRAM, variantes estructurales, evolución por generaciones,
  cruces, lecturas largas (solo Illumina pareadas de hasta 150 pb), ni escala
  para genomas grandes.
- Los alelos generados por deslizamiento salen de linajes independientes que
  parten del alelo de la referencia (genealogía en estrella, sin coalescencia);
  cada copia sortea después su alelo de ese abanico. Las frecuencias no reflejan
  parentesco ni estructura entre alelos.
- Las calibraciones del tartamudeo son aproximadas (`coarse: true` en el
  manifiesto; fuentes citadas allí). No las presentes como exactas.

## Al terminar

Entrega a la persona: la ruta de la carpeta, el `general.toml` usado, un
resumen de lo realizado frente a lo declarado, las cifras pedidas (Ho, He, F_IS,
recuentos de unidades) y las limitaciones que condicionen su lectura.
