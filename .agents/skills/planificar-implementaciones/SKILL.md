---
name: planificar-implementaciones
description: Convertir una especificación ya acordada de un proyecto de análisis genético en un plan de implementación ordenado, dividido en cambios pequeños y verificables. Usar al crear, revisar o actualizar el orden de desarrollo, los hitos, las dependencias, las pausas de revisión y las pruebas necesarias antes de empezar a programar o durante el seguimiento de la implementación.
---

# Planificar implementaciones

## Objetivo

Convertir un diseño ya acordado en un orden de trabajo que pueda ejecutarse paso
a paso. El plan debe indicar qué se construirá, en qué orden, de qué decisión de
la especificación procede y cómo se comprobará.

Redactar el plan en castellano claro para genetistas. No dar por supuestos
conocimientos avanzados de programación, arquitectura de software, pruebas o
estadística.

## Límite principal

No tomar decisiones de diseño dentro del plan de implementación. Si falta una
decisión que cambia lo que se debe construir, detener esa parte del plan y
devolver la cuestión a la especificación.

No presentar como decidido un supuesto elegido solo para poder continuar.
Puede avanzarse con partes independientes, pero el plan completo no estará
listo mientras queden decisiones necesarias sin resolver.

## Antes de escribir el plan

1. Leer por completo la especificación y los documentos que explican cómo se
   organizan las partes del sistema, es decir, su arquitectura técnica.
2. Comprobar que el alcance, el comportamiento esperado y las decisiones
   técnicas necesarias están cerrados.
3. Identificar las secciones concretas que justifican cada parte del trabajo.
4. Localizar componentes existentes que deban reutilizarse.
5. Anotar las decisiones abiertas y resolverlas en la especificación antes de
   planificar el trabajo que dependa de ellas.
6. Respetar las herramientas y versiones ya acordadas. En este proyecto, los
   pasos de Python deben usar `uv` salvo que la especificación indique un cambio.

## Estructura del plan

Mantener esta estructura, salvo que el proyecto ya tenga una plantilla más
específica:

1. **Cabecera y estado.** Indicar título, fecha, estado del plan y documentos de
   diseño de los que procede. Aclarar que el plan ordena la implementación y no
   sustituye a la especificación.
2. **Alcance incluido y excluido.** Enumerar qué se construirá y qué no. Asociar
   cada elemento excluido con un plan posterior o indicar expresamente que no
   está previsto.
3. **Criterios para decidir el orden.** Explicar brevemente por qué los hitos
   siguen esa secuencia.
4. **Requisitos previos.** Indicar qué debe existir y poder comprobarse antes
   del primer paso.
5. **Hitos y pasos.** Dividir el trabajo en grupos de pasos pequeños que puedan
   implementarse, probarse y confirmarse por separado.
6. **Pausas de revisión.** Añadir una revisión obligatoria al terminar cada
   hito, antes de continuar con el siguiente.
7. **Resumen de comprobaciones.** Incluir una tabla que relacione cada hito con
   las pruebas o evidencias que demostrarán que funciona.
8. **Trabajo posterior.** Indicar qué queda para otros planes para que ninguna
   tarea desaparezca sin explicación.

## Cómo escribir cada paso

Usar una casilla por cada cambio que pueda quedar registrado como una unidad.
Esta unidad se denomina *commit* en el sistema de control de versiones:

```markdown
- ☐ A1. Resultado concreto que debe construirse.
  - Depende de: pasos anteriores necesarios, o «ninguno».
  - Procede de: sección concreta de la especificación o de la arquitectura.
  - Se comprueba con: prueba o evidencia esperada.
```

Aplicar estas reglas:

- Describir un resultado observable, no una actividad imprecisa como
  «trabajar en el adaptador».
- Hacer cada paso lo bastante pequeño para implementarlo, revisarlo, probarlo
  y registrarlo sin mezclar cambios independientes.
- Indicar siempre sus dependencias. No confiar solo en el orden de la lista.
- Citar siempre la sección de diseño que lo justifica. Un paso sin origen suele
  señalar una decisión de diseño ausente.
- Indicar cómo se comprobará. Evitar pasos cuya finalización dependa únicamente
  de una opinión.
- Actualizar la casilla de `☐` a `✅` solo cuando el cambio y sus pruebas estén
  terminados.

## Criterios para ordenar el trabajo

Elegir solo los criterios que se apliquen y explicarlos en el plan:

- Definir primero las estructuras de datos y la forma en que se comunican las
  partes del sistema, es decir, sus tipos e interfaces, cuando las necesiten
  los pasos posteriores.
- Probar primero la lógica principal de forma aislada y conectar después los
  archivos, procesos externos y demás entradas y salidas.
- Construir primero la versión más sencilla cuando pueda servir como referencia
  para comprobar versiones posteriores.
- Reutilizar funciones, reglas y constantes existentes antes de proponer una
  segunda implementación.
- Ejecutar antes las comprobaciones rápidas cuando no cambien el significado
  del resultado.
- Avanzar de forma incremental: completar y revisar un hito antes de empezar el
  siguiente.

En proyectos de análisis genético, separar en su propio paso y en su propio
registro cualquier cambio que pueda producir un resultado incorrecto sin
provocar un error visible. Algunos ejemplos habituales son la conversión de
coordenadas, la interpretación de genotipos, la asociación de muestras y la
normalización de variantes. Indicar qué prueba de referencia protege ese
cambio.

## Comprobación de los hitos

Definir pruebas que puedan demostrar el resultado esperado. Preferir, cuando
exista, una verdad conocida, un conjunto de referencia o una implementación ya
validada. Una prueba que solo compara una parte del sistema consigo misma no es
suficiente para demostrar que el resultado biológico es correcto.

Combinar según corresponda:

- pruebas pequeñas para reglas aisladas;
- pruebas de integración para comprobar el recorrido completo;
- casos sintéticos con verdad conocida;
- comparación con una herramienta o resultado de referencia;
- repetición de casos guardados que antes revelaron un fallo;
- comprobación de la interfaz de línea de órdenes y de sus diagnósticos.

La tabla final debe tener, como mínimo, estas columnas:

| Hito | Qué queda construido | Cómo se demuestra |
|---|---|---|
| A | Resultado del hito | Pruebas o evidencia concreta |

## Pausas de revisión

Cerrar cada hito con un bloque como este:

```markdown
> Revisión del hito A: comprobar los cambios y sus pruebas. No continuar hasta
> confirmar el resultado.
```

La pausa permite corregir pronto un error de enfoque y evita que los hitos
posteriores se apoyen en una base no revisada.

## Lenguaje y términos técnicos

- Escribir frases cortas, directas y en voz activa.
- Explicar un término técnico la primera vez que aparezca si puede resultar
  confuso para un genetista.
- Desarrollar las siglas no habituales en su primera aparición.
- Explicar primero para qué sirve un elemento y después dar su nombre técnico.
- Mantener los nombres técnicos necesarios para que el plan sea preciso.
- Distinguir claramente un requisito obligatorio de una recomendación o una
  posibilidad futura.

## Revisión final

Antes de dar el plan por terminado, comprobar:

- que la especificación y la arquitectura resuelven las decisiones necesarias;
- que el alcance incluido y excluido está claro;
- que ningún trabajo excluido se ha perdido sin indicar su destino;
- que cada paso tiene casilla, dependencia, origen y comprobación;
- que los pasos están ordenados por sus dependencias;
- que cada paso es una unidad pequeña que puede revisarse y registrarse;
- que los cambios con riesgo de errores silenciosos están aislados y protegidos
  por una prueba de referencia;
- que cada hito termina con una pausa de revisión;
- que la tabla final explica cómo se demostrará cada hito;
- que el plan puede ejecutarse sin tomar nuevas decisiones de diseño;
- que las herramientas, sus versiones y los comandos respetan lo acordado en
  la especificación, incluido el uso de `uv` en este proyecto;
- que el texto está en castellano claro y resulta comprensible para un
  genetista.

El plan está listo cuando otra persona puede ejecutarlo de principio a fin sin
tener que decidir durante la implementación qué debía hacer el sistema.
