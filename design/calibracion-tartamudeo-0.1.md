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

## 2. Las dos preparaciones, y por qué ahora hacen falta las dos

Cuando se escribió este documento, el proyecto admitía **exclusivamente
preparación sin PCR** y el perfil amplificado estaba aplazado, de modo que la
columna con PCR era solo informativa. **El 19 de septiembre de 2026 el
responsable decidió simular también bibliotecas amplificadas**, así que ahora
**las dos columnas son de uso**, y el diseño correspondiente vive en
`design/diseno-bibliotecas-pcr-0.1.md`.

Que el tartamudeo exista sin PCR no es una excepción: el deslizamiento de la
polimerasa ocurre también en la amplificación en puente sobre la celda y en la
propia síntesis durante la secuenciación, y a ello se suma el mosaicismo
somático de las células de partida. Es menor que con amplificación, pero está —
y la sección 3.1 detalla qué contiene exactamente cada columna, porque de ahí
depende no contar dos veces el mismo fenómeno.

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

### 3.1. Las dos columnas son totales, no componentes

Es el malentendido más fácil de cometer con esta tabla, y llevaría a contar dos
veces el mismo fenómeno.

Los dos conjuntos que compara la fuente pasaron por las **mismas** células, la
**misma** amplificación en puente sobre la celda y la **misma** síntesis
durante la secuenciación. Lo único que los distingue es el paso de PCR de
preparación. Por tanto:

- la **columna sin PCR** ya contiene, sumados y sin posibilidad de separarlos,
  el **mosaicismo somático** de las células de partida, el deslizamiento de la
  **amplificación en puente** y el de la **síntesis**;
- la **columna con PCR** contiene todo lo anterior **más** la PCR de
  biblioteca.

**Decisión del proyecto (19 de septiembre de 2026):** el mosaicismo somático se
queda **dentro** de la tasa sin PCR y **no se modela aparte**. Modelarlo por
separado obligaría a descontarlo de esta cifra, y la medición no permite
hacerlo. Nótese además que el mosaicismo somático no es un artefacto: son
células que de verdad llevan longitudes distintas. Si algún día se modelara,
correspondería a la **verdad** del individuo y no a las lecturas, de modo que
no cabe en una fase cuyo contrato prohíbe que los eventos técnicos alteren la
verdad.

### 3.2. Cómo se componen, sin contar dos veces

Para una biblioteca amplificada **no** se suma nada a la columna con PCR: ya
está todo dentro. Lo que hay que hacer es **despejar** el incremento propio de
la PCR de biblioteca, suponiendo que los dos procesos actúan de forma
independiente:

```
r_pcr = 1 − (1 − r_con) / (1 − r_sin)
```

y aplicar después las dos partes, la de base y la de PCR. La composición
devuelve exactamente la fracción medida: `1 − (1−r_sin)(1−r_pcr) = r_con`.

| Motivo (pb) | Con PCR (medido) | Sin PCR (medido) | Solo la PCR (despejado) | Recompuesto |
|---|---|---|---|---|
| 1 | 0,1700 | 0,0270 | 0,1470 | **0,1700** |
| 2 | 0,0380 | 0,0110 | 0,0273 | **0,0380** |
| 3 | 0,0110 | 0,0029 | 0,0081 | **0,0110** |
| 4 | 0,0069 | 0,0016 | 0,0053 | **0,0069** |
| 5 | 0,0074 | 0,0015 | 0,0059 | **0,0074** |
| 6 | 0,0140 | 0,0010 | 0,0130 | **0,0140** |

**Consecuencia de orden.** Una biblioteca amplificada a la que le falte la
parte de base se queda solo en la columna despejada, y **subestima** el
tartamudeo: un 14 % menos en homopolímeros, 28 % en dinucleótidos, 26 % en
trinucleótidos, 23 % en tetranucleótidos y 20 % en pentanucleótidos. La
preparación con PCR **no puede implantarse sin el tartamudeo sin PCR**; no es
una preferencia, es que la aritmética no cierra.

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
