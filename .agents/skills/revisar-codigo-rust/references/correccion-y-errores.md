# Corrección, errores y concurrencia

## Entradas y valores predeterminados

- Validar rangos, unidades, exclusiones mutuas y dependencias antes de simular.
- Aplicar los valores predeterminados en un único lugar y registrar el valor
  efectivo, no solo el solicitado.
- No usar `Default` para decisiones biológicamente significativas sin que el
  contrato lo declare.
- Rechazar claves desconocidas o estados incompatibles cuando aceptarlos pueda
  ocultar un error de configuración.
- Comprobar sumas, conversiones numéricas, desbordamientos y pérdidas de
  precisión.

## Tratamiento de errores

- No ignorar `Result`, códigos de salida, escritura parcial ni errores al cerrar
  o sincronizar archivos.
- Evitar `unwrap`, `expect`, `panic!`, índices directos y `unreachable!` en rutas
  alcanzables desde entradas del usuario.
- Usar errores tipados cuando el llamador deba distinguir causas y añadir
  contexto en la frontera de la aplicación.
- Conservar la causa original al envolver errores.
- No sustituir un error por una colección vacía, valor cero o advertencia si el
  contrato exige detenerse.
- Limpiar o marcar como incompletas las salidas parciales; no dejarlas con
  apariencia válida.
- Comprobar explícitamente los procesos externos: inicio, salida, señales,
  `stderr`, archivos esperados y compatibilidad de versiones.

## Índices, tamaños y aritmética

- Revisar todas las conversiones entre enteros con signo, sin signo y tamaños
  de plataforma.
- Evitar restar uno antes de validar posiciones basadas en uno.
- Comprobar suma y multiplicación de cobertura, ploidía, longitud y número de
  copias.
- No truncar contadores de lecturas, individuos, variantes o bytes.
- Tratar EOF inesperado, registros truncados y campos ausentes como errores
  definidos.

## `unsafe`

- Exigir una razón por la que no basta código seguro.
- Exigir un comentario `SAFETY` que enumere precondiciones y por qué se cumplen.
- Documentar obligaciones del llamador en funciones o rasgos inseguros.
- Revisar validez, alineamiento, inicialización, aliasing, límites y vida de cada
  puntero o vista.
- Exigir pruebas de regresión para la invariante que sostiene el bloque.

## Concurrencia

- No mantener guardas de bloqueo durante operaciones largas, I/O o esperas.
- Definir un orden total si se adquieren varios bloqueos.
- Justificar cada orden de memoria atómico; usar órdenes débiles solo con una
  relación de sincronización demostrada.
- Usar canales acotados para productores capaces de superar al consumidor.
- Recoger y propagar errores y pánicos de hilos o tareas.
- Revisar que la cancelación o el fallo no deje manifiestos ni archivos marcados
  como completos.
- Confirmar que el paralelismo no cambia la aleatoriedad, el orden contractual
  ni el resultado con la misma semilla.
