# Concurrencia, contención y escalado

## Trabajo y granularidad

- Paralelizar unidades independientes con suficiente trabajo: eventos,
  referencias, cromosomas, poblaciones o bloques, según el contrato.
- Evitar crear una tarea por base, variante o lectura si la coordinación domina
  el cálculo.
- Medir 1, 2, 4 y más hilos hasta el límite de la máquina; informar velocidad y
  eficiencia, no solo el mejor tiempo.
- Mantener una ruta secuencial para entradas pequeñas cuando el umbral se
  justifique mediante medida.

## Estado compartido

- No adquirir un bloqueo dentro de un bucle crítico si puede acumularse trabajo
  local y combinarse después.
- Usar reducción por hilo para contadores y estadísticas frecuentes.
- No sustituir `Mutex` por `RwLock`, un mapa concurrente u otra biblioteca sin
  medir el patrón real de lecturas, escrituras y espera.
- Revisar falsa compartición en contadores o estados por hilo cuando el perfil
  muestre tráfico de caché.
- Justificar órdenes atómicos y separar defectos de corrección de mejoras de
  rendimiento.

## Canales y tuberías

- Usar canales acotados para aplicar contrapresión.
- Medir tamaño de lote, capacidad y equilibrio entre productor y consumidor.
- Evitar bloquear un hilo de Rayon esperando I/O o un canal lleno.
- Separar cálculo, compresión y escritura solo si la tubería mejora el escenario
  completo.

## Exceso de paralelismo

- Contar todos los grupos de hilos vivos: Rayon, minimap2, HTSlib, compresión y
  cualquier proceso externo.
- Evitar que cada grupo asuma que posee todos los núcleos.
- Detectar paralelismo anidado, tareas bloqueantes dentro de iteradores paralelos
  y dos grupos de CPU activos a la vez.
- Confirmar con número de hilos, cambios de contexto, espera y escalado; no
  rediseñar la topología solo por reconocer un patrón.
- Probar primero limitar o repartir hilos; después considerar unificar grupos o
  serializar la fase más barata.

## Determinismo

- No compartir un único generador aleatorio protegido por bloqueo.
- Derivar flujos por identidad lógica, no por el orden en que un hilo recibe la
  tarea.
- Comparar resultados con distinto número de hilos usando la misma semilla.
- Mantener orden estable en archivos cuando el formato o la especificación lo
  exijan; ordenar en una frontera definida si el cálculo paralelo no lo conserva.

## Medición

Registrar tiempo, utilización, hilos, espera, memoria y salida para cada
configuración. Una mejora paralela no es válida si cambia la simulación, pierde
errores o hace depender el resultado del planificador.
