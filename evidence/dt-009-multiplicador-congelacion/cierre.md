# Multiplicador de la PCR en los homopolímeros: lo de antes no se mueve (08/10/2026)

Evidencia de DTC-099, la decisión que hace configurable cuántas veces multiplica
la PCR el tartamudeo de los homopolímeros, y de su requisito REQ-CONF-141: **sin
el campo nuevo, todo lo publicado tiene que ser byte a byte lo de antes**.

## Qué se comparó, y contra qué

Se extrajo el commit de partida, `fc2b205`, con `git archive`, se compiló
**aparte**, con su propio directorio de compilación, y se ejecutaron cinco
escenarios con ese binario **antes de tocar una línea del código**. Después se
ejecutaron los mismos cinco con el binario final. Los dos binarios se lanzaron
desde la raíz del repositorio, con la ruta de cada escenario dentro de
`inputs/`, porque la configuración efectiva registra la ruta relativa a la
invocación (ver la corrección de `evidence/stutter-f4-8-congelacion/cierre.md`).

| Escenario | Qué ejercita |
|---|---|
| `amplificado-loci` | preparación con PCR, un homopolímero de 15 unidades con dos alelos y un dinucleótido de 17, tartamudeo de base, duplicados al 30 %: el tartamudeo clonal entero |
| `amplificado-sustitucion` | lo mismo con la probabilidad por ciclo del dinucleótido sustituida |
| `sin-pcr-loci` | los mismos loci con tartamudeo de base y sin PCR |
| `amplificado-sin-loci` | preparación con PCR sin microsatélites: el bloque clonal no aparece |
| `sin-evento` | sin evento técnico |

## Resultado

**Los 55 ficheros publicados coinciden byte a byte.** Las huellas están en
`linea-base.sha256` (binario de partida) y `con-el-multiplicador.sha256`
(binario final); los dos ficheros son idénticos. Se ejecutó además el binario de
partida dos veces en directorios de salida distintos, con el mismo resultado: la
comparación no depende de dónde se escriba.

`multiplicador-2.5.sha256` es el control positivo: `amplificado-multiplicador`
es `amplificado-loci` con `homopolymer_pcr_fold = 2.5`. Cambian la copia del
TOML, la configuración efectiva, los dos manifiestos, la cabecera de la verdad
—que lleva el identificador de ejecución y la huella de la configuración— y los
R1 de las dos muestras, que son los que leen el homopolímero en este montaje.
Los registros de la verdad, la referencia, la matriz y los R2 no cambian.

## La huella dentro de la batería

`tests/clonal_homopolymer_fold.rs` congela, con los mismos ficheros de
`inputs/amplificado-loci/`, las huellas de los cuatro FASTQ, de la verdad, del
manifiesto de evento y del bloque del tartamudeo clonal de la efectiva. Se
midieron con una sonda ejecutada sobre la copia del commit de partida, no
después del cambio, y coinciden con las de `linea-base.sha256`.

## Cómo reproducirlo

```text
git archive fc2b205 | tar -x -C <copia>
cd <copia> && CARGO_TARGET_DIR=<copia-target> cargo build --locked --bin genome-population-simulator
# desde la raíz del repositorio:
<copia-target>/debug/genome-population-simulator run evidence/dt-009-multiplicador-congelacion/inputs/<escenario>/general.toml --output <salida-base>/<escenario>
cargo build --locked --bin genome-population-simulator
./target/debug/genome-population-simulator run evidence/dt-009-multiplicador-congelacion/inputs/<escenario>/general.toml --output <salida-actual>/<escenario>
# y en cada salida: find . -type f | sort | xargs sha256sum
```

## Lo que esta evidencia no dice

No dice que la configuración efectiva de una ejecución amplificada con
microsatélites pase su esquema: no lo pasa, ni antes ni después, por DT-040 —su
catálogo de azar dice `f4-technical-v5` y el esquema exige `f4-technical-v6`—.
Es anterior a esta entrega y se registró al encontrarlo.

*Nota del 12/10/2026.* DT-040 y DT-041 se cerraron ese día: desde entonces la
efectiva de `amplificado-loci` publica `f4-technical-v6` y pasa su esquema, y
los manifiestos generales de las dos ejecuciones amplificadas pasan el suyo.
Con el binario corregido, en los tres escenarios amplificados con
microsatélites (`amplificado-loci`, `amplificado-multiplicador` y
`amplificado-sustitucion`) cambian solo `effective-config.json` y
`manifest.json`; los demás escenarios no cambian ni un byte. Las huellas de
antes y de después están en el cierre de DT-040 del registro de deudas.
Los ficheros `.sha256` de aquí siguen siendo los de la entrega de DT-009.

*Nota del 13/10/2026.* Con el cierre de DT-034 el texto de la fuente de la
calibración del tartamudeo atribuye las cifras a Gymrek (2016) además de Raz
(2019). Cambian de bytes, en los escenarios con microsatélites
(`sin-pcr-loci`, `amplificado-loci`, `amplificado-multiplicador` y
`amplificado-sustitucion`), `effective-config.json`, `events/reads-main/manifest.json`
y, por arrastre, `manifest.json`; los FASTQ, la verdad y los escenarios sin
microsatélites no cambian. Por eso `linea-base.sha256`, `con-el-multiplicador.sha256`
y `multiplicador-2.5.sha256` guardan huellas anteriores a ese texto y hoy no se
reproducirían tal cual en esos tres ficheros. Las huellas de antes y de después
están en el cierre de DT-034 del registro de deudas.
