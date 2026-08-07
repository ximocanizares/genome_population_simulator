# Memoria, asignaciones y disposición de datos

## Asignaciones

- Localizar `Vec::new`, `String::new`, `format!`, `Box::new`, `collect`, clones y
  mapas creados dentro de bucles frecuentes.
- Reutilizar búferes con `clear` cuando su contenido no escape y conservar la
  capacidad sea seguro.
- Preasignar con una estimación o cota razonable; evitar reservas gigantes por
  un máximo teórico.
- Encadenar iteradores cuando una colección intermedia solo se vuelve a
  recorrer.
- Escribir en un búfer reutilizable en lugar de crear cadenas por registro.
- Aplicar `$auditar-clonaciones-rust` a copias de secuencias, haplotipos,
  segmentos y registros.

## Representación del simulador

- Evitar almacenar una secuencia completa por individuo. Medir representaciones
  basadas en segmentos fundadores, cambios o árboles de procedencia.
- Usar identificadores compactos para referencias, individuos, haplotipos,
  loci, eventos y grupos; resolver nombres en las fronteras.
- Internar nombres repetidos de cromosomas, muestras y referencias cuando una
  medida muestre muchas copias.
- Separar metadatos fríos de campos recorridos en cada generación.
- Considerar vectores paralelos o una estructura compacta cuando un bucle use
  solo pocos campos de una estructura grande; medir antes la pérdida de claridad.
- Representar genotipos y dosis de forma compacta sin perder ploidía, fase ni
  valores ausentes.

## Localidad y caché

- Recorrer datos contiguos en el orden de uso.
- Evitar grafos de objetos pequeños asignados por separado si se recorren
  secuencialmente.
- Agrupar operaciones por cromosoma, población o bloque cuando mejore localidad
  sin cambiar el orden biológico contractual.
- Revisar padding, tamaño de enumeraciones y alineamiento solo después de medir
  memoria o fallos de caché.

## Técnicas condicionadas

- Usar `SmallVec` solo si la distribución real muestra colecciones pequeñas con
  frecuencia alta.
- Usar arenas si muchos objetos comparten exactamente una fase de vida y pueden
  liberarse juntos.
- Usar mapas especializados o hashers alternativos solo tras medir colisiones,
  seguridad y distribución de claves.
- Cambiar el asignador global únicamente con una comparación completa.

## Comprobación

Para cada cambio medir memoria máxima, asignaciones y tiempo. Verificar además
que siguen intactos fase, linaje, orden, referencia y reproducibilidad.
