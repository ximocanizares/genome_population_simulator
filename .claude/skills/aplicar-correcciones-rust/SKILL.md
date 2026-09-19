---
name: aplicar-correcciones-rust
description: Verificar y aplicar de forma conservadora hallazgos procedentes de una revisión de código o rendimiento Rust. Usar cuando exista un informe de revisión y el usuario pida corregirlo, contabilizando cada hallazgo, añadiendo pruebas de regresión y validando los cambios sin copiar ciegamente las soluciones sugeridas.
---

# Aplicar correcciones Rust

## Objetivo

Convertir hallazgos revisados en cambios pequeños, demostrables y trazables. Un
informe es evidencia y orientación, no un parche autorizado ni necesariamente
actualizado.

## Principios

- Verificar cada hallazgo contra el código actual.
- Priorizar corrección y seguridad sobre estilo o velocidad.
- No resolver silenciosamente una ambigüedad de comportamiento.
- Contabilizar todos los hallazgos, incluidos los no aplicados.
- Limitar el diff a la corrección y sus pruebas.
- No afirmar una validación que no se haya ejecutado.
- No cambiar una interfaz pública sin autorización expresa.
- No crear commits ni publicar cambios salvo autorización del usuario.

## Preparación

1. Leer alcance, estado, preguntas y prioridades del informe completo.
2. Extraer para cada hallazgo identificador, gravedad, confianza, supuestos,
   archivos, comportamiento, propuesta, impacto público y prueba posible.
3. Localizar el código por significado, no solo por el número de línea.
4. Confirmar si el defecto sigue presente.
5. Registrar una decisión inicial para cada hallazgo antes de editar.
6. Guardar una línea base de rendimiento antes de cambiar un recorrido medido,
   si existe una prueba adecuada. Si no se guardó a tiempo, declararlo y no
   reconstruir una línea base artificial.

Detener la aplicación si el informe no puede interpretarse, faltan los archivos
o el código ha cambiado tanto que los hallazgos no pueden relacionarse de forma
segura.

## Decisiones iniciales

- **Aplicar:** el problema existe, la solución es clara, no inventa política ni
  cambia una interfaz pública y puede validarse.
- **Preguntar:** el defecto es creíble, pero dos decisiones del usuario producen
  comportamientos distintos.
- **Aplazar:** requiere rediseño, un cambio público o información demasiado
  amplia para una pregunta concreta.
- **Disputar:** la evidencia del código contradice el hallazgo.
- **Ya corregido:** el código actual satisface la regla y se ha comprobado.

La gravedad ayuda a ordenar, pero no convierte una propuesta ambigua en
aplicable.

## Estados finales

Usar exclusivamente:

- **Aplicado**;
- **Aplicado con adaptación**;
- **Ya corregido**;
- **Aplazado**;
- **Disputado**;
- **Validación fallida**;
- **Bloqueado por cambio de contexto**;
- **Sustituido por otro hallazgo**;
- **Pendiente de respuesta**.

Todo estado distinto de “Aplicado” debe incluir una razón y un destino.

## Flujo por hallazgo

Procesar primero bloqueantes y graves. Trabajar de uno en uno, salvo que dos
hallazgos sean inseparables y se explique la agrupación.

1. Volver a confirmar el defecto justo antes de editar.
2. Para errores de corrección, límites o invariantes, añadir primero una prueba
   que falle o demuestre el problema cuando sea práctico.
3. Implementar el cambio mínimo. Adaptar la sugerencia a la API real; no pegar
   código del informe sin verificarlo.
4. Ejecutar la prueba específica.
5. Ejecutar las pruebas del módulo afectado y la validación general.
6. Actualizar inmediatamente el registro del hallazgo.

Si la validación falla por el cambio, corregirlo o revertir completamente el
intento y marcar “Validación fallida”. No dejar una solución parcial.

## Calidad de las pruebas

- Reproducir el defecto real, no un indicador indirecto.
- Comprobar el valor o la variante de error exactos cuando el contrato sea
  estable.
- Cubrir el límite que reveló el problema.
- Usar referencia, coordenadas, ploidía, fase y semilla explícitas en casos
  genómicos.
- Conservar una verdad independiente del código corregido.

## Validación

Ejecutar lo exigido por el proyecto. Como mínimo, cuando esté disponible:

```text
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-targets --all-features
cargo doc --no-deps
```

Usar `cargo audit` si forma parte del entorno. Registrar comandos, códigos de
salida y resultados concisos. Distinguir fallos preexistentes de regresiones
nuevas mediante evidencia.

Si el cambio afecta a una ruta con prueba de rendimiento, comparar contra la
línea base con la misma entrada, semilla, número de hilos y compilación
optimizada. Una corrección de seguridad o exactitud no se revierte
automáticamente por ser más lenta; registrar el compromiso y elevarlo.

## Informe de aplicación

Entregar una tabla con una fila por hallazgo:

| ID | Gravedad | Decisión inicial | Estado final | Archivos | Pruebas | Seguimiento |
|---|---|---|---|---|---|---|

Añadir:

1. resumen de resultados;
2. preguntas y respuestas;
3. explicación por hallazgo;
4. comandos y resultados reales;
5. hallazgos graves pendientes;
6. comparación de rendimiento si correspondía.

Guardar el informe solo si el usuario o el plan lo solicita, en la ubicación
definida por el proyecto. Una ejecución sin archivo guardado sigue necesitando
el mismo recuento completo en la respuesta.

La aplicación está completa cuando ningún hallazgo ha desaparecido, todas las
correcciones tienen pruebas y cada elemento pendiente conserva una razón y una
acción posterior.
