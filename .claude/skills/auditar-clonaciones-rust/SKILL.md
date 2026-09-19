---
name: auditar-clonaciones-rust
description: Auditar código Rust para localizar clonaciones, conversiones y asignaciones innecesarias sin promover refactorizaciones especulativas. Usar ante `.clone()`, `.to_owned()`, `.to_string()`, `.to_vec()`, `Arc::clone`, recopilaciones intermedias o parámetros propietarios, especialmente en recorridos de poblaciones, haplotipos, secuencias y lecturas.
---

# Auditar clonaciones Rust

## Objetivo

Separar las copias realmente costosas de las que son correctas, infrecuentes o
más sencillas que su alternativa. Evaluar conjuntamente coste probable,
frecuencia y complejidad de eliminar la copia.

Esta auditoría no sustituye una revisión completa ni una medición de
rendimiento.

## Preparación

1. Delimitar los archivos y recorridos incluidos.
2. Identificar si el código se ejecuta por base, variante, segmento, individuo,
   generación, lectura o solo una vez por ejecución.
3. Ejecutar, si el proyecto lo admite:

```text
cargo clippy --all-targets --all-features -- \
  -W clippy::redundant_clone \
  -W clippy::clone_on_copy \
  -W clippy::needless_pass_by_value
```

4. Buscar `.clone()`, `.cloned()`, `.to_owned()`, `.to_string()`, `.to_vec()`,
   `Arc::clone`, `Rc::clone`, `format!` y `collect` intermedios.
5. Leer el contexto y los usos posteriores antes de clasificar cada caso.

## Clasificación

### Nivel A: corregir

Usar cuando la copia es claramente innecesaria, la corrección es local y no
cambia la interfaz. Ejemplos:

- clonar un valor que puede moverse y no vuelve a usarse;
- clonar un tipo `Copy`;
- crear un `String` para pasarlo inmediatamente como `&str`;
- recopilar clones cuando basta mover o mantener un iterador;
- clonar una clave justo antes de que deje de existir.

Proponer un cambio concreto y pequeño.

### Nivel B: corregir si es local

Usar cuando la copia puede repetirse o ser grande y se elimina relajando una
firma privada, separando préstamos de campos o reutilizando un búfer. La
corrección no debe atravesar una interfaz pública ni extender vidas útiles por
varios módulos.

En un recorrido interno crítico, una clonación de secuencia, vector de segmentos
o estado haplotípico suele pertenecer aquí si la alternativa es sencilla.

### Nivel C: medir antes

Usar cuando eliminar la copia exige:

- cambiar propiedad o vidas útiles entre módulos;
- introducir `Cow`, `Arc`, arenas o nuevas abstracciones;
- modificar una interfaz pública;
- cambiar propiedades `Send` o `Sync`;
- reestructurar un algoritmo o una representación genómica.

No aplicar el cambio sin una prueba de rendimiento representativa. Proponer la
medida y el criterio de aceptación sin inventar porcentajes.

### Nivel D: conservar

No informar individualmente, salvo que el usuario pida un inventario completo:

- clonaciones en pruebas, diagnósticos o errores;
- configuración y preparación ejecutadas una sola vez;
- `Arc::clone` o `Rc::clone` fuera de un bucle crítico demostrado;
- copias pequeñas cuyo reemplazo empeora la claridad;
- cambios que exigirían `unsafe` o romper una interfaz sin evidencia.

Puede resumirse cuántos casos de este tipo se descartaron.

## Patrones importantes

### Valor que ya puede moverse

Si todas las ramas consumen el mismo valor y este no se usa después, moverlo en
cada rama evita la copia.

### Conversión seguida de préstamo

Eliminar `.to_string()`, `.to_owned()` o `format!` cuando el receptor ya acepta
un préstamo adecuado.

### Recopilación intermedia

Evitar `collect::<Vec<_>>()` cuando el siguiente paso solo recorre el vector y
puede consumir el iterador. Si se necesita propiedad, preferir `.cloned()` o
`.copied()` por claridad, sin atribuirles una mejora que no existe.

### Clonación para evitar un conflicto de préstamos

Comprobar si los campos son disjuntos y pueden prestarse por separado:

```rust
let Estado { segmentos, historial, .. } = estado;
procesar(segmentos, historial);
```

Clasificar como B si el cambio es local; como C si propaga vidas útiles o cambia
la arquitectura.

### Contadores de referencias

`Arc::clone` y `Rc::clone` incrementan un contador, no copian profundamente el
valor. Solo tratarlos como hallazgo si una medición o la estructura del bucle
demuestra una frecuencia relevante.

## Medición para el nivel C

Usar una prueba existente o añadir una que reproduzca el tamaño y el patrón de
acceso reales. Ejecutarla en modo optimizado y comparar antes y después con la
misma entrada, semilla y número de hilos. Registrar tiempo, memoria o
asignaciones según el mecanismo propuesto.

Aceptar el cambio solo si la mejora es reproducible y compensa la complejidad.
Revertirlo si no se distingue del ruido o empeora otro recurso importante.

## Formato de salida

Para cada hallazgo A o B:

```markdown
- `ruta:línea` — **Nivel A|B** — título
  - Coste: qué se copia y con qué frecuencia.
  - Evidencia: uso posterior o contexto observado.
  - Corrección: cambio mínimo propuesto.
```

Para C, sustituir la corrección por un experimento concreto y terminar con
“No aplicar antes de medir”. No mostrar D como hallazgos.

## Límites

- No estimar aceleraciones sin medir.
- No eliminar una clonación en una frontera `unsafe` sin demostrar su contrato.
- No asumir que una implementación de `Clone` carece de efectos; leerla si es
  propia o está documentada de forma especial.
- Si una función contiene muchas clonaciones relacionadas, proponer revisar su
  propiedad global en lugar de emitir una lista mecánica.
- No modificar código durante una auditoría salvo petición expresa.
