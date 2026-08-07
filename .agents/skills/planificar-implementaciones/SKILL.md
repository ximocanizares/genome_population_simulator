---
name: planificar-implementaciones
description: Convertir una especificación y una arquitectura ya acordadas de un proyecto de análisis genético en un plan de implementación ordenado, trazable y verificable. Usar al crear, revisar o actualizar hitos, dependencias, cambios pequeños, pausas de revisión y pruebas, respetando el lenguaje y las herramientas que se hayan elegido sin introducir nuevas decisiones de diseño.
---

# Planificar implementaciones

## Objetivo

Convertir el diseño acordado en una secuencia de cambios pequeños que otra
persona pueda ejecutar sin decidir durante la implementación qué debía hacer el
sistema.

Redactar el plan en castellano claro para genetistas y desarrolladores. Explicar
los términos de software que condicionen el orden o la comprobación del trabajo.

## Puerta de entrada

No empezar el plan completo hasta disponer de:

- una especificación vigente, sin decisiones funcionales críticas abiertas;
- una arquitectura que defina componentes, responsabilidades, contratos,
  flujos de datos e invariantes;
- el lenguaje, las herramientas y las versiones necesarias para escribir
  comandos concretos, si el plan llega a ese nivel.

No resolver decisiones de diseño dentro del plan. Si falta una decisión que
cambia lo que debe construirse, devolverla al documento correspondiente y
continuar solo con las partes independientes.

## Preparación

1. Leer completamente la especificación, la arquitectura y las decisiones que
   las complementen.
2. Inspeccionar el repositorio para localizar componentes reutilizables,
   convenciones, pruebas y comandos existentes.
3. Crear un mapa de dependencias entre capacidades, componentes y formatos.
4. Señalar los cambios capaces de producir resultados incorrectos sin un fallo
   visible, por ejemplo coordenadas, fase, genotipos, muestras o referencias.
5. Confirmar que cada elemento aplazado tiene un destino explícito.

No suponer ningún lenguaje, gestor de dependencias, compilador ni herramienta
mientras no se haya elegido en la arquitectura o en una decisión posterior
aprobada.

## Estructura del plan

Mantener, salvo que exista una plantilla acordada:

1. cabecera, fecha, estado y documentos de origen;
2. alcance incluido, excluido y trabajo posterior;
3. criterios usados para ordenar el trabajo;
4. requisitos previos y decisiones bloqueantes;
5. hitos con pasos pequeños y dependencias explícitas;
6. pausa de revisión al cerrar cada hito;
7. matriz final de hitos y evidencias de validación.

## Formato de cada paso

Usar una casilla por cambio independiente:

```markdown
- ☐ A1. Resultado concreto y observable.
  - Depende de: A0 o «ninguno».
  - Procede de: especificación §X y arquitectura §Y.
  - Cambia: rutas, componentes o interfaces previstos.
  - Se comprueba con: comando, prueba y evidencia esperada.
```

Aplicar estas reglas:

- describir un resultado, no una actividad imprecisa;
- limitar el paso a una unidad que pueda implementarse, revisar y validar por
  separado;
- citar siempre el origen; un paso sin origen suele ocultar diseño nuevo;
- declarar dependencias aunque el orden de la lista parezca evidente;
- indicar archivos o componentes concretos cuando la arquitectura permita
  conocerlos;
- usar comprobaciones reales del proyecto, no marcadores genéricos;
- cambiar `☐` por `✅` solo después de completar el cambio y sus pruebas.

Un paso puede corresponder a un *commit*, pero el plan no autoriza por sí solo a
crear commits ni a publicar cambios. Respetar siempre la autorización del
usuario y el flujo de trabajo del repositorio.

## Orden de trabajo

Elegir y explicar solo los criterios pertinentes:

- establecer primero contratos, formatos y estructuras que necesiten los pasos
  posteriores;
- probar la lógica central de forma aislada antes de conectar archivos,
  procesos externos o interfaces;
- construir una versión sencilla de referencia antes de optimizaciones o
  variantes complejas;
- reutilizar componentes existentes antes de crear implementaciones paralelas;
- ejecutar pronto las comprobaciones rápidas que detecten errores estructurales;
- cerrar y revisar un hito antes de apoyar otro sobre él.

Aislar en pasos propios los cambios con riesgo de error silencioso y protegerlos
con una verdad conocida, un caso sintético independiente o una referencia
externa adecuada. No aceptar como única prueba que una parte del sistema se
compare consigo misma.

## Validación

Combinar según corresponda:

- pruebas unitarias de reglas aisladas;
- pruebas de integración de recorridos completos;
- casos sintéticos con verdad conocida y semillas reproducibles;
- pruebas de formatos e interoperabilidad con herramientas externas;
- pruebas de regresión para fallos ya encontrados;
- comprobaciones de errores, diagnósticos y casos límite;
- medidas de rendimiento solo cuando exista un requisito medible.

Definir para cada comprobación el comando, el resultado esperado y la evidencia
que se conservará. Usar únicamente el gestor de dependencias, compilador,
entorno y comandos elegidos por el proyecto.

Cerrar cada hito con una pausa explícita:

```markdown
> Revisión del hito A: comprobar los cambios y sus evidencias. No continuar
> hasta confirmar el resultado.
```

La matriz final debe incluir al menos:

| Hito | Resultado construido | Evidencia | Riesgo protegido |
|---|---|---|---|
| A | Resultado del hito | Pruebas o comprobación concreta | Riesgo principal |

## Comprobación final

- Confirmar que el plan no decide comportamiento ni arquitectura nuevos.
- Comprobar que cada paso tiene dependencia, origen, cambio y validación.
- Verificar que los pasos respetan el orden técnico real.
- Confirmar que los riesgos silenciosos tienen pruebas independientes.
- Revisar que cada hito termina con una pausa y una evidencia observable.
- Verificar que el trabajo excluido conserva un destino.
- Confirmar que herramientas y comandos coinciden con las decisiones vigentes.
- Eliminar pasos demasiado grandes, duplicados o imposibles de verificar.
- Comprobar que el texto está en castellano y es comprensible para sus dos
  públicos.

El plan está listo cuando puede ejecutarse de principio a fin sin inventar
decisiones y cada avance puede demostrarse antes de continuar.
