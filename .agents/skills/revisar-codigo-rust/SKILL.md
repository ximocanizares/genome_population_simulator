---
name: revisar-codigo-rust
description: Revisar código Rust del simulador para detectar defectos de corrección, errores genómicos silenciosos, problemas de diseño, concurrencia, fiabilidad, pruebas y dependencias. Usar al revisar un diff, commit, módulo, funcionalidad o repositorio Rust; producir hallazgos priorizados y verificables sin modificar el código salvo petición expresa.
---

# Revisar código Rust

## Objetivo

Determinar si el cambio cumple su intención y conserva las invariantes del
simulador. Informar primero de defectos reales y riesgos demostrables; evitar
listas mecánicas de preferencias estilísticas.

Una revisión es una tarea de lectura. No aplicar correcciones, crear commits ni
editar informes salvo que el usuario lo solicite.

## Referencias

Leer siempre [formato-hallazgos.md](references/formato-hallazgos.md).

Seleccionar después las listas pertinentes:

- Corrección, errores, valores predeterminados y concurrencia:
  [correccion-y-errores.md](references/correccion-y-errores.md).
- Tipos, Rust idiomático, módulos, nombres y refactorización:
  [diseno-y-legibilidad.md](references/diseno-y-legibilidad.md).
- Pruebas, dependencias y herramientas:
  [pruebas-y-herramientas.md](references/pruebas-y-herramientas.md).
- Referencias, coordenadas, haplotipos, semillas y formatos del simulador:
  [dominio-genomico.md](references/dominio-genomico.md).

Leer una referencia completa cuando su categoría aparezca en el alcance. En una
revisión general, leerlas todas.

## Principios

- Priorizar corrección, seguridad, pérdida de datos y resultados biológicamente
  incorrectos.
- Revisar contra la especificación y la arquitectura, no contra preferencias
  personales.
- Citar únicamente rutas, líneas, comandos y comportamientos comprobados.
- No fabricar resultados de compilación, pruebas ni perfiles.
- Explicar el mecanismo del fallo y aportar una corrección concreta.
- No informar el mismo defecto varias veces; elegir la categoría principal.
- Reconocer lo que está bien cuando protege una invariante importante.
- Mantenerse dentro del alcance y separar observaciones externas.

## Procedimiento

### 1. Establecer el alcance

Identificar archivos, diff o commits incluidos, código de producción, pruebas,
configuración y documentos de diseño aplicables. Anotar exclusiones y cambios
ajenos del directorio de trabajo.

### 2. Determinar la intención

Resumir en un párrafo qué pretende hacer el cambio, qué contrato ofrece y qué
invariantes debe conservar. Si la intención no puede deducirse de la petición,
la especificación, las pruebas o el código, plantear una pregunta concreta antes
de emitir juicios sobre diseño.

### 3. Inventariar

Localizar:

- interfaces públicas y tipos modificados;
- rutas de entrada, transformación y salida;
- `unsafe`, concurrencia, aleatoriedad y procesos externos;
- parsers, serializadores y conversiones de coordenadas;
- asignaciones, clonaciones y materializaciones en bucles;
- pruebas existentes y huecos de cobertura;
- dependencias o características nuevas.

### 4. Ejecutar comprobaciones

Usar los comandos definidos por el proyecto. Cuando estén disponibles:

```text
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-targets --all-features
cargo doc --no-deps
```

Añadir pruebas específicas, `cargo audit` y combinaciones de características si
forman parte del alcance. Registrar los resultados reales y no repetir comandos
que ya dispongan de una salida fiable aportada por el usuario.

### 5. Realizar pasadas temáticas

Revisar por categorías para reducir omisiones. Si se permite usar agentes,
pueden repartirse categorías independientes, pero cada uno debe recibir el
mismo alcance, intención y formato. La revisión principal debe verificar y
deduplicar sus hallazgos. Si no se permite, realizar las pasadas de forma
secuencial.

### 6. Desafiar las pruebas

Para cada función no trivial modificada, buscar al menos una clase de entrada no
cubierta y explicar qué defecto revelaría. No exigir pruebas irrelevantes a una
función trivial ni confundir cobertura de líneas con corrección.

### 7. Sintetizar

Ordenar hallazgos por gravedad y después por confianza. Unificar causas comunes,
separar observaciones fuera de alcance y ofrecer comandos concretos para
reproducir o verificar.

## Formato del informe

1. Alcance y documentos de autoridad.
2. Veredicto breve.
3. Estado de compilación y pruebas.
4. Preguntas o supuestos.
5. Tres prioridades principales.
6. Hallazgos ordenados por gravedad.
7. Observaciones fuera de alcance.
8. Pruebas que faltan.
9. Aspectos correctamente resueltos.
10. Comandos para volver a verificar.

Guardar el informe únicamente cuando se solicite. Usar la ubicación acordada
por el proyecto y enlaces relativos válidos. En caso contrario, entregarlo en la
respuesta.

## Cierre

Antes de terminar:

- comprobar que cada hallazgo cita código leído;
- eliminar hallazgos puramente hipotéticos sin mecanismo ni prueba;
- confirmar que la gravedad coincide con el impacto;
- diferenciar fallos preexistentes de los introducidos por el cambio;
- verificar que no se mezclan recomendaciones de rendimiento sin medición con
  defectos de corrección;
- indicar expresamente si no hay hallazgos.
