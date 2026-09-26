# Congelación acotada de la fase de realismo técnico — tartamudeo (24/09/2026)

Esta es la evidencia que el punto **10.4** del plan de microsatélites y
amplificación exige repetir al cerrar la entrega 2: la fase de realismo técnico
estaba cerrada con su evidencia aprobada, se le ha añadido un paso, y hay que
demostrar que **un escenario sin loci declarados publica exactamente los mismos
bytes que antes**.

## Qué se comparó, y contra qué

Se extrajo el commit **anterior** al del tartamudeo (`c4a55f3~1`) con
`git archive`, se compiló **aparte**, con su propio directorio de compilación,
y se ejecutó el mismo escenario con los dos binarios. Compilarlo aparte no es
un adorno: con el mismo directorio, el binario de una de las dos versiones
sobreviviría a la otra y la comparación no significaría nada.

El escenario está íntegro en `inputs/`. Es deliberadamente **poco favorable**:
preparación sin PCR, **ningún locus de repetición declarado**, y los tres
artefactos técnicos que sí consumen aleatoriedad activos a la vez —quimera al
5 %, salto de índice al 2 % y contaminación biológica del 1 al 3 %— más dos
individuos diploides, una matriz SNV suministrada y vista VCF comprimida con su
índice. Un escenario neutro habría dado la misma respuesta sin ejercitar nada.

## Resultado

**Los catorce artefactos publicados coinciden byte a byte.** Las huellas están
en `linea-base.sha256` y `con-la-capacidad.sha256`; los dos ficheros son
idénticos, línea a línea.

| Artefacto | Qué demostraría un cambio |
|---|---|
| Los cuatro FASTQ finales | que la etapa nueva consume aleatoriedad o desplaza el ordinal de plantilla |
| `outputs/truth.jsonl` | que el tartamudeo toca el genotipo, que es lo que nunca debe hacer |
| `outputs/truth.vcf.gz` y su índice `.tbi` | lo mismo, en la vista publicada |
| `effective-config.json` | que la configuración efectiva publica algo del tartamudeo sin haberlo |
| `manifest.json` y `events/reads-main/manifest.json` | que la versión del documento o el catálogo de ámbitos subieron sin motivo |
| `resources/references/ref-main.fa`, la matriz y sus metadatos | que algo aguas arriba se movió |

La razón de que esto sea alcanzable, y no una casualidad, está en DTC-081: cada
sorteo arranca de su propia clave y «no existe consumo acumulado entre
decisiones». Una etapa nueva no puede desplazar un flujo existente. Y el
catálogo `f4-technical-v4` **solo** presta su nombre como ámbito a las dos
decisiones del tartamudeo: las trece heredadas siguen derivándose bajo el v3.

## Cómo reproducirlo

**Desde la raíz del repositorio, y no desde otro sitio.** La configuración
efectiva registra las rutas **relativas a la invocación**, de modo que ejecutar
lo mismo desde otro directorio cambia dos artefactos —la configuración efectiva
y el manifiesto general, que resume su huella— sin que nada del simulador haya
cambiado. Siguiendo la receta al pie de la letra desde otro sitio, las huellas
archivadas no salen.

> **Corrección del 25/09/2026, al rehacer esta evidencia para la entrega 3.**
> La receta de abajo **no reproduce las huellas archivadas**: da doce iguales y
> dos distintas, `effective-config.json` y `manifest.json`. No es un defecto del
> simulador —cuatro commits distintos dan las catorce iguales entre sí—, sino de
> la receta: la ejecución archivada **no** copió el escenario a `inputs/`, sino
> que lo invocó por su ruta real dentro del directorio de evidencia,
> `evidence/stutter-f4-8-congelacion/inputs/general.toml`. Con esa ruta, las
> catorce salen. Es justo lo que el párrafo de arriba advierte y esta receta no
> aplicaba. La receta correcta está en
> `evidence/pcr-f4-9-congelacion/cierre.md`.

```text
git archive c4a55f3~1 | tar -x -C <copia>
cd <copia> && CARGO_TARGET_DIR=<copia-target> cargo build --locked --bin genome-population-simulator
<copia-target>/debug/genome-population-simulator run inputs/general.toml --output <salida-base>
cargo build --locked --bin genome-population-simulator
./target/debug/genome-population-simulator run inputs/general.toml --output <salida-actual>
diff <(cd <salida-base> && find . -type f | sort | xargs sha256sum) \
     <(cd <salida-actual> && find . -type f | sort | xargs sha256sum)
```

## Lo que esta evidencia **no** dice

No dice que un escenario **con** loci declarados publique los mismos bytes, y no
debe decirlo: declarar un locus sube el documento general a la versión 4 y con
él la configuración efectiva y el manifiesto general, de modo que exigir su
identidad sería pedir lo imposible. Esa acotación es la contradicción **C-1** y
está recogida en REQ-DAT-80. Lo que sí se exige con un locus declarado cuyo
único alelo es el de la referencia —FASTQ distintos, verdad y VCF idénticos— lo
comprueba `tests/stutter_pipeline.rs`, no esta evidencia.
