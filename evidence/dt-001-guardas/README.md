# Evidencia de DT-001: las guardas de los verificadores independientes

Qué hay aquí y para qué sirve cada cosa. Todo se midió el 23 de septiembre de
2026 sobre el árbol de trabajo, con `bcftools` y `tabix` 1.24 en `/usr/local` y
las dependencias de desarrollo de Python instaladas.

## El barrido

`scripts/check-guard-coverage.sh` corre la batería con el barrido encendido y
cuenta, por cada juez independiente, cuántos de sus puntos de rechazo llegan a
dispararse alguna vez. Tres ejecuciones, en orden:

| Fichero | Cuándo | Oráculo grande | Oráculo del VCF | Auditoría de identidad |
|---|---|---|---|---|
| `barrido-0-*` | antes de tocar nada | 252 de 365 sin disparar | 39 de 77 | 11 de 42 |
| `barrido-1-*` | tras cubrir el oráculo del VCF | 252 de 365 | 18 de 77 | 11 de 42 |
| `barrido-3-*` | al terminar | **168 de 365** | **25 de 98** | **5 de 42** |

El denominador del verificador del VCF cambia de 77 a 98 entre la segunda
medición y la tercera, y no porque se añadieran comprobaciones: la revisión
adversarial encontró que trece rechazos se construían a mano en vez de pasar
por una constructora, de modo que **el barrido no los veía**. Al enderezarlos,
el denominador pasó a ser el verdadero y salieron a la luz veinte guardas
muertas que nadie podía contar.

De cada ejecución se guarda la tabla (`-tabla.md`), la lista exacta de las
guardas que no se dispararon (`-sin-disparar.txt`) y el recuento de la batería
(`-bateria.txt`).

**Lo que este barrido no mide.** Mide si una guarda **se dispara**, no si una
prueba **fallaría** en caso de que la guarda estuviera mal escrita. Una prueba
que provoca un rechazo sin exigir **cuál** rechazo lo satisface igual que una
que exige el motivo exacto. Por eso hace falta lo siguiente.

## El estudio de mutación

`mutaciones.md` es la tabla de los resultados de las **132** mutaciones y
`mutaciones.json` su lista, cada una con el texto exacto que sustituye y la
línea del árbol final a la que apunta. `arnes-de-mutacion.py` es el arnés que
las aplica **de una en una**, y toma el árbol y el directorio de la copia de
seguridad como argumentos: con la ruta del repositorio escrita a fuego, quien
lo ejecutara mañana mutaría el árbol real sin saberlo.

1. huella SHA-256 del fichero antes;
2. sustitución única —si el texto no aparece exactamente una vez, aborta—;
3. huella después, que tiene que diferir: así se demuestra que la mutación mutó;
4. ejecución de la prueba que dice cubrir esa guarda, que tiene que **fallar**;
5. restauración desde una copia de seguridad **fuera de git**, y comprobación de
   que la huella vuelve a ser la de antes.

Dos trampas que este arnés aprendió a evitar, y que conviene no repetir:

- **Una mutación que no compila no demuestra nada.** El fallo es del
  compilador, no de la prueba. El arnés lo detecta y lo marca como nulo.
- **Restaurar conservando la fecha del fichero deja el binario mutado.** `cargo`
  compara fechas y cree que su compilación sigue al día, de modo que el árbol
  queda restaurado y la biblioteca compilada no. Ocurrió una vez durante este
  trabajo y obligó a revalidar tres corpus enteros. El arnés toca la fecha
  después de restaurar.

## Cómo reproducirlo

```bash
./scripts/check-guard-coverage.sh <directorio-nuevo>
python3 evidence/dt-001-guardas/arnes-de-mutacion.py \
    evidence/dt-001-guardas/mutaciones.json . /tmp/copia-de-guardas
```

El segundo argumento es el árbol que se va a mutar y el tercero el directorio
de la copia de seguridad, que **tiene que estar fuera de ese árbol**: el arnés
lo comprueba y se niega si no lo está.
