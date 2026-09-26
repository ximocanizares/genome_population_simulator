# Congelación acotada de la fase de realismo técnico — preparación amplificada (25/09/2026)

Esta es la evidencia que el punto **10.4** del plan de microsatélites y
amplificación exige repetir al cerrar la **entrega 3**: la fase de realismo
técnico estaba cerrada con su evidencia aprobada, se le ha añadido la
preparación amplificada por PCR, y hay que demostrar que **un escenario que no
la declara publica exactamente los mismos bytes que antes**.

## Resultado

**Los catorce artefactos publicados coinciden byte a byte** con los de la
entrega 2. `linea-base.sha256` y `con-la-amplificacion.sha256` son idénticos,
línea a línea.

El escenario es el mismo de la entrega 2, y sigue siendo deliberadamente **poco
favorable**: preparación sin PCR, **ningún locus de repetición declarado** y los
tres artefactos técnicos que sí consumen aleatoriedad activos a la vez —quimera
al 5 %, salto de índice al 2 % y contaminación biológica del 1 al 3 %— más dos
individuos diploides, una matriz SNV suministrada y vista VCF comprimida con su
índice.

## Cuatro commits, las mismas catorce huellas

Para que esto no dependa de una sola medición, el escenario se ejecutó con
**cuatro** versiones distintas del programa, cada una compilada en su propio
directorio:

| Commit | Qué es | `effective-config.json` | `manifest.json` |
|---|---|---|---|
| `c4a55f3` | el tartamudeo de base | `80af2982…` | `afe31f90…` |
| `ca35e6a` | cierre de la entrega 2 | `80af2982…` | `afe31f90…` |
| `855f0d7` | corte del sorteo del inserto | `80af2982…` | `afe31f90…` |
| árbol de la entrega 3 | esta entrega | `80af2982…` | `afe31f90…` |

Las cuatro coinciden entre sí y con las catorce huellas archivadas.

## La receta de la entrega 2 estaba mal escrita, y aquí queda corregida

Al rehacer esta evidencia apareció una discrepancia que conviene dejar
explicada, porque cuesta media hora descubrirla y parece un defecto del
programa sin serlo.

Siguiendo **al pie de la letra** la receta que archivó la entrega 2 —copiar el
escenario a `inputs/` en la raíz del repositorio y ejecutar
`run inputs/general.toml`— salen **doce** huellas iguales y **dos distintas**:
`effective-config.json` y `manifest.json`, que resume la huella del anterior.

No es un defecto del simulador ni de ningún commit. La causa está en la propia
`cierre.md` de la entrega 2, que la advierte sin aplicarla a su propia receta:
**la configuración efectiva registra la ruta del fichero de configuración tal
como se escribió en la invocación**. La ejecución archivada no usó
`inputs/general.toml`: usó la ruta real del escenario dentro del directorio de
evidencia. Comprobado, las tres invocaciones dan tres huellas distintas y solo
una es la archivada:

```text
evidence/stutter-f4-8-congelacion/inputs/general.toml  ->  93d7fb14…   (la archivada)
inputs/general.toml                                    ->  c717c96a…
/home/…/genome_population_simulator/inputs/general.toml ->  b4c2f731…
```

**La evidencia de la entrega 2 reproduce perfectamente**; lo que no reproducía
era su receta. La de abajo es la correcta y no necesita copiar nada.

## Cómo reproducirlo

**El escenario vive en `evidence/stutter-f4-8-congelacion/inputs/`, y esta
carpeta no tiene copia propia a propósito.** Es el mismo escenario de la entrega
2 —tiene que serlo, o no habría nada que comparar—, y **la ruta con la que se
invoca forma parte del resultado**: duplicarlo aquí y llamarlo por su ruta nueva
cambiaría `effective-config.json` y `manifest.json` sin que nada del simulador
hubiera cambiado.

> Lo comprobé en mis propias carnes al escribir esta evidencia: copié el
> escenario a esta carpeta, puse la ruta nueva en la receta, y las dos huellas
> se movieron. Es **exactamente** el defecto que este documento describe tres
> párrafos más arriba, cometido por quien lo estaba describiendo. Por eso la
> copia ya no está y la receta apunta al original.

Desde la raíz del repositorio, y con la **misma** ruta en las dos invocaciones:

```text
git archive 855f0d7 | tar -x -C <copia>
cd <copia> && CARGO_TARGET_DIR=<copia-target> cargo build --locked --bin genome-population-simulator
<copia-target>/debug/genome-population-simulator \
    run evidence/stutter-f4-8-congelacion/inputs/general.toml --output <salida-base>

cd <repo> && cargo build --locked --bin genome-population-simulator
./target/debug/genome-population-simulator \
    run evidence/stutter-f4-8-congelacion/inputs/general.toml --output <salida-actual>

diff <(cd <salida-base>   && find . -type f | sort | xargs sha256sum) \
     <(cd <salida-actual> && find . -type f | sort | xargs sha256sum)
```

Con esa invocación, `con-la-amplificacion.sha256` y `linea-base.sha256` salen
idénticos, línea a línea.

## Lo que esta evidencia **no** dice

No dice que un escenario **con** preparación amplificada publique los mismos
bytes: no debe, y ahí está el asunto. Lo que sí se exige de la preparación
amplificada con todos sus parámetros neutros —fracción de duplicados cero y
curva neutra— es que sus **FASTQ** sean byte a byte los de la preparación sin
PCR, y eso lo comprueba `tests/pcr_pipeline.rs`, con los errores de
secuenciación encendidos y apagados.

Tampoco dice nada de un escenario con loci de repetición declarados: esa
acotación es la contradicción **C-1** y sigue igual.

## El ancla que esta evidencia no puede dar, y quién la da

Esta evidencia compara **dos binarios** entre sí. Es necesaria y no es
suficiente: un cambio que moviera los dos por igual le sería invisible. El ancla
**absoluta** —huellas congeladas escritas a mano, del fichero comprimido, del
contenido descomprimido y de la primera pareja que el muestreador produce por su
camino real— vive en `tests/pcr_frozen_identity.rs`, corre en cada batería y se
comprobó que también pasa contra el commit anterior a esta entrega.
