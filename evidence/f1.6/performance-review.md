# Revisión de rendimiento F1.6

> **Expediente histórico no vigente.** Las observaciones pueden orientar la
> reauditoría, pero no acreditan por sí solas el runtime simplificado actual.

Fecha: 14 de agosto de 2026.

- El selector directo proporcional compone los huecos mediante separadores
  sorteados; no crea un vector ni recorre una vez por cada base del FASTA real.
- La materialización reserva el número exacto de segmentos comprobado y copia
  únicamente las bases del objetivo derivado.
- La validación de publicación relee el mapa serializado y reconstruye las bases
  una vez; añade tiempo y memoria lineales en el tamaño del objetivo para evitar
  acreditar bytes distintos de los escritos. No duplica el espacio de
  candidatos ni altera el algoritmo de selección.
- La elección uniforme usa rechazo entero y evita el sesgo de `%` directo.
- El coste del muestreo de huecos es O(bloques log bloques), no depende de las
  bases de holgura. El benchmark compara 100 MB con 10, 100 y 1.000 bloques y
  exige hashes estables entre repeticiones.
- El vector canónico y la reconstrucción independiente fijan el contrato; el
  selector no promete compatibilidad con el algoritmo transitorio anterior.

Resultado: sin regresión de memoria proporcional al tamaño del genoma ni
hallazgo de rendimiento que invalide el contrato funcional de F1.6. La posible
optimización del recorrido de huecos queda como trabajo medible, no como promesa
de rendimiento cerrada.
