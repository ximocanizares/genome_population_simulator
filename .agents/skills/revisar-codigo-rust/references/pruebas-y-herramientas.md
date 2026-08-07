# Pruebas, dependencias y herramientas

## Cobertura de comportamiento

- Exigir pruebas para interfaces públicas y lógica privada no trivial.
- Cubrir recorrido normal, cada error estable, vacío, uno, máximos, límites y
  entrada mal formada.
- Añadir pruebas de regresión para defectos corregidos e invariantes
  documentadas.
- Usar nombres que expresen condición y resultado esperado.
- Evitar pruebas dependientes del tiempo real, red, orden del sistema o sueño.
- Fijar la semilla y mostrarla al fallar una prueba aleatoria.

## Propiedades y generación de casos

- Usar pruebas basadas en propiedades o fuzzing para parsers, serializadores,
  conversiones de coordenadas, normalización y operaciones con leyes de ida y
  vuelta.
- Definir oráculos independientes; no comparar dos caminos que comparten la
  misma función defectuosa.
- Incluir entradas mínimas capaces de revelar desbordamientos, solapamientos,
  pérdidas de fase y mezcla de referencias.
- Probar concurrencia cuando haya estado compartido, canales, atomics o código
  `unsafe`.

## Desafío de pruebas

Para cada función modificada no trivial, indicar una clase de entrada que no
esté cubierta, el defecto concreto que podría revelar y una prueba propuesta.
No inflar la lista con casos equivalentes.

## Dependencias

- Justificar qué aporta cada dependencia nueva y por qué no basta la biblioteca
  estándar o una dependencia ya presente.
- Revisar mantenimiento, licencia, características, dependencias transitivas y
  compatibilidad con la versión mínima de Rust acordada.
- Mantener privadas las estructuras de bibliotecas externas cuando exponerlas
  convierta su sustitución en un cambio público.
- Probar las combinaciones de características realmente admitidas.
- Confirmar las dependencias nativas y la distribución en los sistemas objetivo,
  especialmente para HTSlib, compresión y CRAM.

## Herramientas y política

- Ejecutar formato, Clippy, pruebas y documentación con las características
  definidas por el proyecto.
- Usar `-D warnings` en integración continua, sin depender de
  `#![deny(warnings)]` frente a cambios de compilador.
- No permitir `unsafe` si el crate pretende ser seguro sin declararlo y
  auditarlo expresamente.
- Revisar que archivos de bloqueo, versión mínima y perfiles de compilación
  correspondan al tipo de aplicación.
- Usar registros estructurados para diagnósticos; reservar la salida estándar
  para el resultado principal de la interfaz de línea de órdenes.

## Resultado de las herramientas

Registrar comando, versión relevante, código de salida y resumen. Si una
herramienta no está instalada o un fallo es preexistente, indicarlo sin
presentarlo como validación superada.
