# Diseño, Rust idiomático y legibilidad

## Tipos e invariantes

- Hacer irrepresentables los estados inválidos mediante tipos nuevos,
  enumeraciones y constructores validados.
- En parámetros o estructuras con varios booleanos relacionados, enumerar las
  combinaciones. Sustituirlos si alguna combinación carece de sentido.
- No transportar cantidades de dominio como enteros desnudos cuando puedan
  confundirse unidades o identidades.
- Mantener privados los campos cuya escritura directa rompa invariantes.
- Marcar con `must_use` los resultados cuya omisión pierda una obligación.

## Propiedad

- Preferir préstamos y cortes a copias cuando el propietario sea claro.
- No añadir `'static` para silenciar el compilador.
- Usar `Arc` únicamente cuando exista propiedad compartida real; no ocultar una
  jerarquía de propiedad mal definida.
- Revisar cada clonación grande en un bucle mediante
  `$auditar-clonaciones-rust`.
- Aceptar una copia pequeña y fría si la alternativa propaga vidas útiles o
  reduce la claridad sin beneficio medido.

## Control y patrones

- Usar correspondencias exhaustivas para enumeraciones internas.
- Evitar comodines que absorban silenciosamente variantes futuras.
- Preferir retornos tempranos y `let ... else` a anidamientos profundos.
- Limitar la mutabilidad al tramo donde se construye o modifica el valor.
- No recopilar iteradores en una colección si el siguiente paso puede consumir
  el flujo directamente.

## Módulos y dependencias

- Mantener una dirección clara: los módulos de datos y formatos no deben
  depender de fases concretas del flujo.
- Colocar los tipos usados por varias fases en una frontera compartida.
- Evitar rasgos con una sola implementación salvo que definan una frontera
  pública o una sustitución prevista.
- Preferir rutas absolutas desde la raíz del crate para referencias lejanas.
- Evitar módulos genéricos llamados `utils`, `helpers`, `manager` o `common` sin
  una responsabilidad de dominio.
- Mantener módulos hoja como archivos planos y evitar `foo/foo.rs`.

## Nombres

- Nombrar funciones con verbos y tipos con sustantivos del dominio.
- Evitar `data`, `item`, `value`, `thing`, `tmp`, `handler` y abreviaturas no
  acordadas.
- Nombrar booleanos como predicados: `is_`, `has_`, `should_`.
- Usar el mismo término para el mismo concepto en todo el crate.
- Extraer literales con significado, unidades o necesidad de sincronización a
  constantes documentadas.

## Seguridad ante refactorizaciones

- Evitar `..Default::default()` cuando un campo futuro requiera una decisión.
- Desestructurar exhaustivamente en implementaciones manuales de igualdad,
  hash, clonación o depuración que dependan del conjunto de campos.
- Usar patrones de cortes para formas fijas en lugar de comprobar longitud y
  después indexar.
- Revisar comodines y `..` que permitirían compilar tras añadir una variante o
  campo sin adaptar el comportamiento.

## Olores

- Examinar funciones largas o con demasiadas decisiones, estructuras cuyos
  métodos usan grupos disjuntos de campos y lógica repetida en tres lugares.
- Reemplazar parámetros booleanos posicionales por enumeraciones o opciones con
  nombre.
- Reemplazar cadenas que representan conjuntos cerrados por enumeraciones.
- Eliminar código muerto y bloques comentados.
- Exigir justificación y condición de retirada para cada `allow`, `TODO` o
  `FIXME`.
- No usar `Deref` para simular herencia de métodos en una envoltura normal.
