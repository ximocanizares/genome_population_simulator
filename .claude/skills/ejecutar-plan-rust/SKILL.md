---
name: ejecutar-plan-rust
description: Ejecutar paso a paso un plan escrito de implementación del simulador en Rust, integrando implementación, revisión, correcciones y validación. Usar cuando el usuario pida implementar un plan, un hito o una secuencia de pasos ya aprobados, conservando su trazabilidad y deteniéndose ante contradicciones de diseño o en las pausas de revisión.
---

# Ejecutar plan Rust

## Objetivo

Llevar un plan aprobado hasta código Rust validado, un paso cada vez. Coordinar
las habilidades del proyecto sin duplicar sus listas de comprobación:

1. `$implementar-funcionalidad-rust`;
2. `$revisar-codigo-rust`;
3. `$aplicar-correcciones-rust`;
4. `$revisar-rendimiento-rust` cuando el paso tenga un requisito medible o
   modifique un recorrido crítico.

## Condiciones previas

Detener el inicio si no existe:

- un plan escrito con pasos, dependencias y comprobaciones;
- una especificación que defina el comportamiento;
- una arquitectura que resuelva componentes y contratos;
- un proyecto Rust inicializado para los pasos que requieran compilar.

Leer completamente esos documentos y comprobar que sus requisitos previos
coinciden con el repositorio real. El código prevalece para describir el estado
actual; la especificación prevalece para el comportamiento acordado.

## Regla de divergencia

Aceptar y registrar una adaptación local que conserve el propósito: un nombre
mejor, una función auxiliar distinta o el ajuste a la firma real de un componente
reutilizado.

Detenerse y pedir una decisión cuando:

- el plan contradiga la especificación;
- falte una regla que cambie el comportamiento observable;
- el cambio alcance componentes o pasos no previstos;
- sea necesario modificar un contrato público o una invariante compartida;
- una prueba no pueda corregirse sin cambiar el diseño;
- se alcance una pausa de revisión prevista.

No alterar silenciosamente el orden, fusionar pasos independientes, eliminar
alcance ni inventar requisitos.

## Ciclo por paso

### 1. Preparar

- Confirmar identificador, resultado, dependencias, origen y comprobación.
- Revisar el estado del directorio de trabajo y preservar cambios ajenos.
- Identificar riesgos de referencia, coordenadas, fase, ploidía, semillas,
  memoria y concurrencia.

### 2. Implementar

Aplicar `$implementar-funcionalidad-rust` solo al paso actual. Añadir tipos,
código y pruebas necesarios. Registrar cualquier adaptación local.

### 3. Validar

Ejecutar primero las pruebas específicas y después las comprobaciones generales
disponibles:

```text
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-targets --all-features
```

No continuar con una validación roja ni atribuir al cambio fallos preexistentes
sin comprobarlo.

### 4. Revisar

Aplicar `$revisar-codigo-rust` al diff del paso. No filtrar previamente los
hallazgos. Si el cambio toca un recorrido crítico, aplicar también
`$revisar-rendimiento-rust` con una medida representativa.

### 5. Corregir

Aplicar `$aplicar-correcciones-rust` al informe. Resolver hallazgos claros y
verificables; pedir decisión para los que cambien diseño o comportamiento.
Repetir la validación tras las correcciones.

### 6. Cerrar

- Marcar el paso `✅` únicamente cuando código, pruebas y revisión estén
  completos.
- Resumir archivos, pruebas, comandos y divergencias.
- Crear un commit solo si el usuario lo ha autorizado. Si no, dejar el paso
  validado y listo para confirmar.
- Empezar el paso siguiente únicamente después de cerrar el actual.

## Granularidad y pausas

Respetar la unidad mínima definida por el plan. Fusionar dos pasos solo si son
inseparables para compilar o probar, explicándolo antes o en el informe. No
convertir un hito completo en un único cambio difícil de revisar.

Las pausas de revisión son obligatorias. Al terminar el último paso del hito,
entregar:

- pasos completados;
- resultados de validación;
- commits, si estaban autorizados;
- decisiones locales registradas;
- riesgos o tareas pendientes;
- siguiente hito previsto.

Esperar la revisión del usuario antes de continuar.

## Finalización del plan

Al completar el último paso:

1. ejecutar la validación completa en el entorno objetivo disponible;
2. comprobar que todas las casillas y evidencias coinciden con el estado real;
3. revisar que no queden hallazgos críticos sin destino;
4. resumir los pasos, commits autorizados, pruebas y trabajo aplazado.

El plan no está completo si un paso solo “parece funcionar”, si una validación
no se ejecutó o si una divergencia de diseño quedó oculta.
