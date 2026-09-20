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
| Media | 0,051 | **0,011** | 4,64× |

**La tasa decrece con la longitud del motivo**, y los homopolímeros son con
diferencia los peores.

### 3.1. Las dos columnas son totales, no componentes

Es el malentendido más fácil de cometer con esta tabla, y llevaría a contar dos
veces el mismo fenómeno.

Cada uno de los dos conjuntos que compara la fuente pasó por una celda de flujo
con su **amplificación en puente** y por la **síntesis** de la secuenciación. La
diferencia **de protocolo** entre ellos es el paso de PCR de preparación,
aunque, como advierte el recuadro de abajo, no es la única diferencia entre los
dos conjuntos. Por tanto:

- la **columna sin PCR** ya contiene, sumados y sin posibilidad de separarlos,
  el **mosaicismo somático** de las células de partida, el deslizamiento de la
  **amplificación en puente** y el de la **síntesis**;
- la **columna con PCR** contiene todo lo anterior **más** la PCR de
  biblioteca.

> ### Corrección del 20 de septiembre de 2026: no son las mismas células
>
> Hasta esta fecha, este documento afirmaba que «los dos conjuntos pasaron por
> las **mismas** células, la **misma** amplificación en puente y la **misma**
> síntesis», y que «lo único que los distingue es el paso de PCR». **Es falso**,
> y conviene decirlo porque de ahí se derivan las tasas por ciclo del modelo con
> PCR.
>
> Son **dos individuos distintos**, de proyectos distintos: uno de Papuasia
> Nueva Guinea secuenciado **con** PCR a 25,9× y otro de Pakistán **sin** PCR a
> 37,6×, y el primero además en una corrida agrupada con otras muestras.
>
> **Qué sigue en pie.** Que cada columna es un **total** y no un componente,
> que es lo que obliga a **restar** en vez de sumar, y con ello toda la
> aritmética de la sección 3.2 y la consecuencia de orden que cierra esta
> sección.
>
> **Qué se debilita.** Atribuir la diferencia entre columnas **solo** a la PCR.
> La resta puede llevar dentro efectos de individuo y de lote además del efecto
> de la amplificación. En consecuencia, **las tasas por ciclo derivadas de esta
> tabla son toscas**, y deben publicarse marcadas como tales en el manifiesto
> de cualquier ejecución que las use. No se propone corregirlas: no hay con qué,
> y una corrección inventada sería peor que una cifra tosca declarada como tal.

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

**Dependencia con la longitud del tracto: pendiente de digitalizar, no de
acceso.** La tasa también **crece con la longitud del tracto** hasta
estabilizarse, y Fungtammasan y colaboradores (2015) lo describen.

> **Corrección del 20 de septiembre de 2026.** Este documento afirmaba que ese
> artículo está «tras muro de pago» y que «no se ha podido extraer la curva
> ajustada». **Es falso**: el artículo está en **acceso abierto** en Europe PMC
> con el identificador **PMC4417121** y licencia **CC BY-NC**, y la curva **sí
> se puede recuperar**, de su figura 1A. Lo único que falta es el trabajo de
> digitalizarla, y la licencia no lo impide mientras el uso no sea comercial.

Mientras no se haga, la dependencia con la longitud del tracto queda como
parámetro configurable **sin valor por defecto derivado de fuente**, y no debe
inventarse uno. Nótese además que el predictor que la literatura forense
considera mejor no es la longitud total del tracto sino la **tirada
ininterrumpida más larga** del motivo.

---

## 5. Forma del modelo

**Actualización del 20 de septiembre de 2026.** Existe un modelo **publicado y
en licencia CC BY 4.0** que infiere la forma midiendo el ruido directamente, y
es el que el proyecto adopta: **Raz O., Biezuner T., Spiro A. y col. (2019)**,
*Short tandem repeat stutter model inferred from direct measurement of in vitro
stutter noise*, Nucleic Acids Research 47(5):2436–2445,
doi:10.1093/nar/gky1318. Es una cadena de Markov en la que **cada ciclo de PCR
es un paso**, con probabilidades de transición **linealmente dependientes de la
longitud actual del tracto**, contracciones de hasta tres unidades por ciclo y
expansiones de una sola. **No publican los valores de los parámetros**, solo el
método de ajuste: se toma la forma, no los números. Y sus predictores incluyen
el **tipo** de motivo, no solo su longitud, lo que es un desajuste conocido con
la tabla de la sección 3, que va por longitud.

Ojo a una confusión fácil: la dependencia **lineal** es de la probabilidad **por
ciclo**; la curva con meseta que describe la sección 4 es la **fracción
observada de extremo a extremo**, que sale de acumular la primera sobre muchos
ciclos. Las dos cosas son ciertas y no se contradicen.

La descripción que sigue, la de HipSTR, se conserva porque es la forma en que
la mayoría de las herramientas de genotipado expresan el mismo fenómeno y
ayuda a leer sus salidas. Son tres parámetros:

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
| Raz y col. (2019), *Nucleic Acids Research* | **La forma del modelo**: cadena de Markov por ciclo, dependencia lineal, asimetría del paso | **CC BY 4.0** | Atribuir |
| Fungtammasan y col. (2015), *Genome Research* | La afirmación sobre la longitud del tracto y, en el futuro, la curva de su figura 1A | **Acceso abierto**, PMC4417121, bajo **CC BY-NC** (corregido el 20/09/2026: no es de pago, y la licencia sí se conoce) | Atribuir. La cláusula **no comercial** permite reproducir la curva en este proyecto mientras su uso no lo sea; si alguna vez lo fuera, hay que revisarlo |
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

- **Raz O., Biezuner T., Spiro A. y col. (2019).** *Short tandem repeat stutter
  model inferred from direct measurement of in vitro stutter noise.* Nucleic
  Acids Research 47(5):2436–2445. <https://doi.org/10.1093/nar/gky1318>
  — **CC BY 4.0**. Origen de la **forma** del modelo que el proyecto adopta.

- **Fungtammasan A. y col. (2015).** *Accurate typing of short tandem repeats
  from genome-wide sequencing data and its applications.* Genome Research.
  PMID 25823460. **Acceso abierto: PMC4417121, licencia CC BY-NC.**
  <https://pubmed.ncbi.nlm.nih.gov/25823460/>
  Corrobora de forma independiente la ventaja del protocolo sin PCR («hasta
  nueve veces menos errores») y describe el crecimiento de la tasa con la
  longitud del tracto, con la curva en su figura 1A.

- **HipSTR.** <https://hipstr-tool.github.io/HipSTR/> y
  <https://github.com/gymrek-lab/HipSTR> — código bajo GPL v2.
  Origen de la **forma** del modelo de tres parámetros de la sección 5.

---

## 8. Lo que este documento no resuelve

- **Dónde ocurre el deslizamiento en el recorrido**, y por tanto si las dos
  lecturas de un par deben coincidir. Las fuentes no lo zanjan: es una decisión
  de diseño del simulador y vive en el expediente de decisiones de repeticiones
  en tándem.
- **La dependencia con la longitud del tracto**, sin valor por defecto: la
  fuente **sí** es accesible y su licencia **sí** se conoce, CC BY-NC; lo que
  falta es digitalizar la curva de su figura 1A, como se explica en la
  sección 4.
- **Si estos valores valen para tomate**, más allá del argumento de la sección 1
  de que el artefacto depende del tracto y no del organismo.

---

## 9. La dependencia con la longitud del tracto, digitalizada (20/09/2026)

Cierra parcialmente lo que estaba pendiente. La fuente es **Fungtammasan y
col. (2015)**, figura 1A, panel superior izquierdo: proporción de lecturas
erróneas en **escala logarítmica** frente a la longitud del tracto en pares de
base, con curvas separadas para mononucleótidos y dinucleótidos, con y sin PCR.

**Valores leídos de la figura**, como `log₁₀` de la proporción. Mononucleótidos:

| Tracto (pb) | Sin PCR | Con PCR |
|---|---|---|
| 5 | −5,4 | −5,3 |
| 7 | −4,1 | −3,8 |
| 9 | −3,0 | −2,7 |
| 11 | −2,2 | −1,8 |
| 13 | −1,65 | −1,25 |
| 15 | −1,3 | −0,85 |

Dinucleótidos, que abarcan tractos más largos:

| Tracto (pb) | Sin PCR | Con PCR |
|---|---|---|
| 10 | −3,6 | −3,4 |
| 18 | −2,3 | −2,0 |
| 26 | −1,6 | −1,35 |
| 34 | −1,25 | −0,9 |

**La forma, que es lo más sólido de todo esto:** en escala logarítmica la curva
**sube casi recta hasta unos 10–12 pb y después se aplana**. Eso es crecimiento
exponencial con meseta, y es coherente con una probabilidad **lineal por
ciclo**: al acumular ciclos, una tasa lineal produce una fracción observada que
satura. Las dos literaturas no se contradicen.

**Contraste con la tabla de la sección 3, que es tranquilizador.** Gymrek mide
0,17 con PCR en homopolímeros, es decir `log₁₀ = −0,77`; la figura da −0,85 a
15 pb. Y 0,027 sin PCR es `−1,57`, frente a −1,3 leído. Mismo orden, con las
diferencias esperables entre muestras y profundidades distintas.

### Tres advertencias, y la tercera impide cerrar el pendiente del todo

1. **Estos valores se leyeron a ojo de una figura**, en escala logarítmica y a
   partir de una imagen de 680 píxeles de ancho. La incertidumbre es de unas
   **±0,15 unidades logarítmicas**, que en lineal es del orden del ±40 %. Sirven
   para la forma y el orden de magnitud, **no como calibración fina**.
2. **Los números exactos existen y son accesibles.** Los autores publican los
   ficheros de tasas de error en un historial de Galaxy enlazado desde el
   artículo. No se pudieron descargar automáticamente porque la página es una
   aplicación interactiva, pero **quien quiera precisión debe bajarlos de ahí**
   en vez de usar esta tabla.
3. **El artículo no dice cuántos ciclos de PCR usó.** Y esto es lo decisivo:
   estas cifras son **de extremo a extremo**, y nuestro modelo necesita la
   **pendiente por ciclo**. Convertir unas en otra exige conocer el número de
   ciclos, que no está. Por tanto el coeficiente **sigue arrancando en cero** y
   el pendiente sigue abierto en su parte cuantitativa, aunque ya no por falta
   de acceso.

**Condiciones de la medición**, que conviene tener presentes: los datos con PCR
son de un varón humano a 150× y los de sin PCR de la muestra NA12882 a 245×,
ambos en el mismo tipo de instrumento con lecturas de 100 bases en pareja. Son
**muestras y profundidades distintas**, la misma salvedad que afecta a la tabla
de la sección 3.

---

## 10. La fracción de duplicados: hay teoría, y valida la fórmula (20/09/2026)

Se buscó una fuente mejor para el valor de duplicados y apareció algo más útil
que un número: **una teoría con anclajes publicados**.

**Rochette N. C., Rivera-Colón A. G., Walsh J., Sanger T. J.,
Campbell-Staton S. C. y Catchen J. M. (2022–2023).** *On the causes,
consequences, and avoidance of PCR duplicates: towards a theory of library
complexity.* bioRxiv doi:10.1101/2022.10.10.511638; publicado en *Molecular
Ecology Resources* (2023), PMID 37062860.

**Su hallazgo central desmonta la intuición habitual.** La tasa de duplicados
la determina **la razón entre profundidad y complejidad** —parejas de lecturas
secuenciadas dividido por moléculas distintas de partida—, y **solo esa razón
importa**, no los valores absolutos. El **número de ciclos de PCR juega un
papel menor**: el ruido de amplificación solo influye apreciablemente cuando la
razón está entre un décimo y uno.

**Y eso valida el parámetro interno de este diseño**, porque las parejas por
molécula **son** exactamente esa razón. No es una variable inventada por
conveniencia del modelo: es la magnitud que la teoría señala como la que manda.

**Anclajes publicados frente a nuestra fórmula.** Los autores dan dos, «con
independencia del ruido de la PCR»:

| Razón | Publicado | Nuestra fórmula `1 − (1 − e^(−m))/m` |
|---|---|---|
| 0,1 | **5–15 %** | **4,8 %** |
| 1,0 | **35–60 %** | **36,8 %** |

Nuestra fórmula cae **en el borde inferior de los dos rangos**, que es
exactamente donde debe caer un modelo **sin ruido de amplificación**: el ruido
real solo puede añadir duplicados, nunca quitarlos. La coincidencia del extremo
inferior en los dos anclajes es una comprobación independiente de que la
fórmula es la correcta y de que su sesgo es el esperado y va en la dirección
segura.

**Ejemplo trabajado, ya con respaldo**, para la documentación del parámetro:

| Fracción declarada | Razón implicada |
|---|---|
| 1 % | 0,02 |
| 2,5 % | 0,05 |
| 9,4 % | 0,2 |
| 21,3 % | 0,5 |

**Sigue sin fijarse un valor por omisión**, y ahora por una razón mejor
fundada: la teoría demuestra que la tasa depende de la **complejidad de la
biblioteca y la profundidad de cada experimento**, no de una constante del
protocolo. Un valor por omisión afirmaría algo que la propia teoría niega. El
número correcto es el del informe de calidad del laboratorio.

**Nota de licencia:** el programa que acompaña al artículo está bajo licencia
GPL versión 3; de aquí no se toma código, solo los anclajes numéricos y el
resultado teórico, que se citan.
