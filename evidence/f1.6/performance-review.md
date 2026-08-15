# Revisión de rendimiento F1.6

Fecha: 14 de agosto de 2026.

- El selector representa huecos y clases de candidatos mediante rangos
  compactos; no crea un vector por base del FASTA real.
- La materialización reserva el número exacto de segmentos comprobado y copia
  únicamente las bases del objetivo derivado.
- La validación de publicación relee el mapa serializado y reconstruye las bases
  una vez; añade tiempo y memoria lineales en el tamaño del objetivo para evitar
  acreditar bytes distintos de los escritos. No duplica el espacio de
  candidatos ni altera el algoritmo de selección.
- La elección uniforme usa rechazo entero y evita el sesgo de `%` directo.
- El coste depende del número de bloques seleccionados y huecos, no del número
  de bases candidatas enumeradas. La búsqueda segura vuelve a recorrer los
  huecos en cada bloque y, por tanto, puede crecer cuadráticamente con el número
  de segmentos en el peor caso. F1.6 no fija un umbral humano de rendimiento;
  esta observación se conserva para medir antes de prometer escalabilidad a
  cientos de miles de bloques.
- No se aplicó una optimización especulativa que cambiase el espacio muestral.
  Los oráculos exhaustivos pequeños, el vector congelado y la prueba estadística
  tienen prioridad sobre una aceleración sin equivalencia demostrada.

Resultado: sin regresión de memoria proporcional al tamaño del genoma ni
hallazgo de rendimiento que invalide el contrato funcional de F1.6. La posible
optimización del recorrido de huecos queda como trabajo medible, no como promesa
de rendimiento cerrada.
