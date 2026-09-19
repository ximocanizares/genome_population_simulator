# Calibración del tartamudeo en repeticiones en tándem

Propuesta de valores por defecto para el artefacto de deslizamiento de la
polimerasa («tartamudeo», *stutter*) en repeticiones en tándem, con sus fuentes
publicadas y la situación legal de cada una.

**Estado: propuesta.** Ninguna cifra de este documento es normativa todavía. Se
redacta para que el responsable del proyecto decida si las adopta, con qué
nombre y bajo qué requisito funcional.

---

## 1. Advertencia de alcance: estos valores son humanos

**Todas las cifras de este documento proceden de secuenciación de genoma humano
completo en Illumina.** No son universales y no deben publicarse como si lo
fueran. Si se adoptan, se adoptan con esa etiqueta explícita, igual que el
proyecto ya hace con sus dos perfiles humanos aproximados, que declaran su
procedencia en el nombre.

Matiz a favor de trasladarlas a otras especies: el tartamudeo **no es una
propiedad del organismo**, sino del tracto repetido y de la polimerasa que lo
copia. Por eso traslada mejor de lo que trasladaría, por ejemplo, una tasa de
mutación. Aun así, la composición de microsatélites de un genoma vegetal difiere
de la del humano, de modo que la calibración debe seguir siendo **configurable**
y estos valores solo el punto de partida.

---

## 2. El protocolo que importa aquí es el que no lleva PCR

El proyecto admite **exclusivamente preparación de biblioteca sin PCR**, y el
perfil amplificado está aplazado. Por eso la columna que se propone adoptar es
la del protocolo sin PCR. La columna con PCR se conserva solo para documentar la
diferencia y por si algún día se reabre el perfil amplificado.

Que el tartamudeo exista sin PCR no es una excepción: el deslizamiento de la
polimerasa ocurre también en la amplificación en agrupación sobre la celda y en
la propia síntesis durante la secuenciación. Es menor, pero está.

---

## 3. Probabilidad de tartamudeo por longitud del motivo

Fracción de lecturas afectadas por tartamudeo, según la longitud de la unidad de
repetición. Fuente: Gymrek (2016), tabla 1.

| Motivo (pb) | Con PCR | **Sin PCR — valor propuesto** | Reducción |
|---|---|---|---|
| 1 (homopolímero) | 0,17 | **0,027** | 6,30× |
| 2 | 0,038 | **0,011** | 3,45× |
| 3 | 0,011 | **0,0029** | 3,79× |
| 4 | 0,0069 | **0,0016** | 4,31× |
| 5 | 0,0074 | **0,0015** | 4,93× |
| 6 | 0,014 | **0,0010** | 14,0× |
| Media | 0,051 | **0,011** | 4,63× |

**La tasa decrece con la longitud del motivo**, y los homopolímeros son con
diferencia los peores.

**No adoptar el 14× de los hexanucleótidos como dato.** El propio autor advierte
que esa diferencia «probablemente se debe a ruido derivado de la escasa cantidad
de loci analizados». La probabilidad sin PCR de 0,0010 se propone igualmente,
pero marcada como la menos fiable de la tabla.

---

## 4. Tamaño y dirección del error

Dos hechos que la misma fuente establece y que el modelo debe respetar:

- El error viene **en múltiplos enteros de la unidad de repetición**.
- Hay un **sesgo fuerte hacia perder repeticiones** en vez de ganarlas.

Y un hallazgo que simplifica el diseño: **la distribución del tamaño del error
es parecida en los dos protocolos**, con un sesgo hacia la contracción
ligeramente menor sin PCR. Es decir, entre biblioteca amplificada y sin
amplificar cambia sobre todo **la frecuencia**, no la forma. Un modelo que
comparta la distribución de tamaños y varíe solo la tasa cubre ambos casos.

**Pendiente de fuente abierta:** la tasa también **crece con la longitud del
tracto** hasta estabilizarse. Fungtammasan y colaboradores (2015) lo describen,
pero el artículo está tras muro de pago y **no se ha podido extraer la curva
ajustada**. Mientras no se consiga, la dependencia con la longitud del tracto
queda como parámetro configurable sin valor por defecto derivado de fuente, y no
debe inventarse uno.

---

## 5. Forma del modelo

Se propone la forma que usa HipSTR, herramienta ampliamente empleada para
genotipar repeticiones a partir de secuenciación Illumina. Son tres parámetros:

1. probabilidad de que el tartamudeo **añada** repeticiones;
2. probabilidad de que las **quite**;
3. un parámetro **geométrico** que gobierna cuántas unidades se ganan o pierden.

Las dos primeras suman la probabilidad total de la tabla de la sección 3, y el
sesgo entre ambas recoge la preferencia por la contracción. La tercera produce
la caída rápida al alejarse del alelo verdadero.

---

## 6. Situación legal de cada fuente

| Fuente | Qué se toma | Licencia | Obligación |
|---|---|---|---|
| Gymrek (2016), preprint en bioRxiv | **Las cifras** de la tabla 1 | **CC BY 4.0** | Atribuir |
| Fungtammasan y col. (2015), *Genome Research* | Solo la afirmación cualitativa sobre la longitud del tracto | Artículo de pago | Citar; no se reproduce contenido |
| HipSTR | **Solo la forma del modelo**, tres parámetros | Código bajo **GPL v2** | Citar el método |

**Sobre HipSTR, para que conste.** Se comprobó su repositorio: **no distribuye
ningún fichero de parámetros precalculados**; los aprende de los datos de cada
usuario mediante un algoritmo iterativo. Por tanto **no se toma dato alguno de
esa herramienta**. Lo que se reutiliza es la descripción publicada de su modelo,
que es un método y se cita, no se licencia. Su código está bajo GPL versión 2,
licencia que solo obligaría si se copiara código, cosa que **no se hace ni debe
hacerse**: la implementación es propia, escrita a partir de la descripción del
método.

---

## 7. Referencias

- **Gymrek M. (2016).** *PCR-free library preparation greatly reduces stutter
  noise at short tandem repeats.* bioRxiv 043448.
  <https://doi.org/10.1101/043448> — **CC BY 4.0**.
  Compara dos conjuntos de genoma completo a alta cobertura, con y sin el
  protocolo sin PCR de Illumina. Es la fuente de toda la sección 3 y de los
  hechos de la sección 4.

- **Fungtammasan A. y col. (2015).** *Accurate typing of short tandem repeats
  from genome-wide sequencing data and its applications.* Genome Research.
  PMID 25823460. <https://pubmed.ncbi.nlm.nih.gov/25823460/>
  Corrobora de forma independiente la ventaja del protocolo sin PCR («hasta
  nueve veces menos errores») y describe el crecimiento de la tasa con la
  longitud del tracto.

- **HipSTR.** <https://hipstr-tool.github.io/HipSTR/> y
  <https://github.com/gymrek-lab/HipSTR> — código bajo GPL v2.
  Origen de la **forma** del modelo de tres parámetros de la sección 5.

---

## 8. Lo que este documento no resuelve

- **Dónde ocurre el deslizamiento en el recorrido**, y por tanto si las dos
  lecturas de un par deben coincidir. Las fuentes no lo zanjan: es una decisión
  de diseño del simulador y vive en el expediente de decisiones de repeticiones
  en tándem.
- **La dependencia con la longitud del tracto**, sin valor por defecto derivado
  de fuente abierta, como se explica en la sección 4.
- **Si estos valores valen para tomate**, más allá del argumento de la sección 1
  de que el artefacto depende del tracto y no del organismo.
