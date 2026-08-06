# Especificación del simulador de genomas y poblaciones

- **Versión:** 0.1
- **Estado:** borrador actualizado
- **Fecha:** 6 de agosto de 2026
- **Fuentes disponibles:** `README.md` e información confirmada por el usuario

## 1. Estado y finalidad del documento

Este documento define el objetivo general y los requisitos confirmados hasta
ahora. También identifica las decisiones que faltan antes de poder cerrar la
especificación y planificar la implementación.

La versión 0.1 sigue siendo un borrador. Los requisitos de esta versión son
obligatorios, pero todavía no se han acordado todos sus detalles ni las medidas
cuantitativas que permitirán aceptar los resultados.

## 2. Objetivo del proyecto

El sistema debe generar datos genómicos simulados con una verdad conocida. La
verdad conocida es la información exacta, definida durante la simulación, sobre
los genotipos reales de los individuos simulados. Esta información servirá como
referencia para comprobar los resultados de otras herramientas.

El objetivo es evaluar de manera rigurosa la eficacia y la fiabilidad de
programas externos que detectan SNP y determinan sus genotipos a partir de
datos de secuenciación. Estos programas se denominan *SNP callers* en este
documento.

La evaluación comparará los resultados de los *SNP callers* con la verdad
conocida. Las medidas concretas de eficacia y fiabilidad todavía no están
definidas.

## 3. Requisitos confirmados

### 3.1. Genoma de referencia

- **REQ-REF-01.** El sistema debe permitir configurar la composición y la
  estructura del genoma de referencia.
- **REQ-REF-02.** La estructura configurable debe permitir incluir
  duplicaciones, deleciones e inversiones.

Estos requisitos no determinan todavía respecto a qué secuencia de partida se
definen las duplicaciones, deleciones e inversiones.

### 3.2. Individuos y grupos

- **REQ-IND-01.** El sistema debe crear individuos que contengan cambios
  respecto al genoma de referencia.
- **REQ-POB-01.** El sistema debe permitir crear poblaciones naturales.
- **REQ-POB-02.** El sistema debe permitir crear poblaciones segregantes, es
  decir, poblaciones en las que los cambios genéticos se distribuyen entre los
  descendientes.
- **REQ-POB-03.** El sistema debe permitir crear familias.

No se han definido la ploidía, las clases exactas de cambios de los individuos
ni los modelos que describirán cada tipo de grupo.

### 3.3. Datos por individuo

- **REQ-DAT-01.** El sistema debe generar un archivo BAM o CRAM por individuo
  respecto al genoma de referencia.
- **REQ-DAT-02.** El sistema debe permitir introducir errores, contaminación,
  baja cobertura, mala calidad de mapeo y otras degradaciones en los archivos
  BAM o CRAM.

La expresión «otras degradaciones» todavía no identifica comportamientos
concretos. Deberá detallarse antes de convertirla en requisitos comprobables.

### 3.4. VCF de verdad

- **REQ-VER-01.** El sistema debe generar un VCF de verdad que represente el
  genotipo real simulado de los individuos.

El VCF de verdad es el archivo que recoge los genotipos definidos por la
simulación. Actúa como referencia para comparar el resultado de un *SNP caller*
con lo que realmente se simuló.

## 4. Flujo funcional confirmado

El flujo funcional describe la relación entre los resultados, no el orden
técnico en que deberá implementarse el sistema.

1. Se configura la composición y la estructura del genoma de referencia.
2. Se crean individuos con cambios respecto a esa referencia.
3. Los individuos se organizan como poblaciones naturales, poblaciones
   segregantes o familias, según el escenario que se quiera simular.
4. Se genera un BAM o CRAM por individuo respecto a la referencia.
5. Se pueden introducir degradaciones en los BAM o CRAM, incluidas las
   categorías ya confirmadas.
6. Se genera el VCF de verdad con los genotipos reales simulados.
7. Un *SNP caller* externo analiza los datos y sus resultados se comparan con el
   VCF de verdad.

Todavía no se ha decidido si el propio sistema ejecutará los *SNP callers* y
realizará la comparación, o si preparará los datos para que esas tareas se hagan
por separado.

## 5. Salidas confirmadas

| Salida | Contenido confirmado | Información pendiente |
|---|---|---|
| Genoma de referencia | Genoma configurable en composición y estructura | Formato, organismo, secuencia de partida y reglas de configuración |
| Estructura de los individuos y grupos | Individuos con cambios respecto a la referencia, organizados en los tipos de grupo confirmados | Forma de representación y modelos biológicos |
| BAM o CRAM por individuo | Datos de cada individuo respecto a la referencia, con posibilidad de degradación | Tecnología de secuenciación, parámetros y características técnicas |
| VCF de verdad | Genotipo real simulado de los individuos | Clases de variantes incluidas y reglas detalladas de representación |

No se ha decidido si una ejecución producirá BAM, CRAM o permitirá producir
ambos. Tampoco se ha decidido si se conservarán a la vez archivos sin degradar
y archivos degradados.

## 6. Criterios de aceptación preliminares

Un criterio de aceptación es una comprobación observable que permite decidir si
un requisito se ha cumplido. Los siguientes criterios pueden deducirse de los
requisitos confirmados. Las cantidades, tolerancias y medidas concretas siguen
pendientes.

| Criterio | Evidencia preliminar | Aspectos pendientes |
|---|---|---|
| **ACE-01. Configuración de la referencia** | Se puede generar una referencia con la composición y estructura solicitadas, incluidas duplicaciones, deleciones e inversiones | Casos de prueba, tamaños y reglas exactas |
| **ACE-02. Creación de individuos** | Los individuos contienen los cambios simulados respecto a la referencia | Clases de cambios, ploidía y número de individuos |
| **ACE-03. Tipos de grupos** | Se puede generar al menos un caso comprobable de población natural, población segregante y familia | Modelos, relaciones y parámetros de cada tipo |
| **ACE-04. Datos por individuo** | Cada individuo tiene el BAM o CRAM correspondiente respecto a la referencia y puede asociarse sin ambigüedad con él | Contenido requerido y comprobaciones técnicas del archivo |
| **ACE-05. Degradaciones** | Se pueden generar BAM o CRAM afectados por errores, contaminación, baja cobertura y mala calidad de mapeo | Definición, intensidad y medida de cada degradación |
| **ACE-06. VCF de verdad** | El VCF de verdad representa los genotipos que el sistema ha simulado para los individuos | Representación de cada clase de variante y comprobaciones de concordancia |
| **ACE-07. Uso para evaluación** | Los resultados de un *SNP caller* externo pueden compararse con el VCF de verdad | Procedimiento de comparación, normalización, medidas y umbrales |

Estos criterios son cualitativos. No deben interpretarse como pruebas cerradas
hasta que se resuelvan los aspectos pendientes.

## 7. Decisiones abiertas

1. **Tipos de variantes evaluados.** ¿El objetivo se limita a evaluar SNP o
   también debe evaluar otros tipos de variante?
2. **Referencia frente a individuos.** ¿Las duplicaciones, deleciones e
   inversiones son solo elementos configurables del genoma de referencia,
   también pueden ser cambios de los individuos, o cumplen ambas funciones?
3. **Base de las alteraciones de la referencia.** ¿Respecto a qué secuencia se
   definen las duplicaciones, deleciones e inversiones de la referencia?
4. **Resultado mínimo útil.** ¿Qué conjunto mínimo de datos y comprobaciones
   debe producir el sistema para que permita una evaluación útil?
5. **Contexto biológico.** ¿Qué organismos, ploidías, tamaños de genoma y
   escalas de población deben admitirse?
6. **Cambios de los individuos.** ¿Qué clases de cambios respecto a la
   referencia deben poder simularse?
7. **Modelos de grupos.** ¿Qué propiedades deben representar las poblaciones
   naturales, las poblaciones segregantes y las familias?
8. **Generación de BAM y CRAM.** ¿Qué tecnologías de secuenciación y qué
   características de los datos deben simularse?
9. **Degradaciones.** ¿Cómo se define y controla cada degradación, y qué otras
   degradaciones deben incluirse?
10. **Salidas conservadas.** ¿Deben producirse BAM y CRAM, o uno de ellos por
    ejecución? ¿Deben conservarse versiones sin degradar y degradadas?
11. **Alcance de la evaluación.** ¿El sistema debe ejecutar los *SNP callers*,
    comparar sus resultados y elaborar medidas, o solo generar los datos y la
    verdad conocida?
12. **Medidas de eficacia y fiabilidad.** ¿Qué medidas y umbrales se usarán
    para evaluar un *SNP caller*?
13. **Representación en el VCF de verdad.** ¿Cómo se representarán los cambios
    de los individuos y cómo se relacionarán con las alteraciones estructurales
    de la referencia?
14. **Uso del sistema.** ¿Qué datos y parámetros aportará la persona usuaria y
    cómo solicitará una simulación?

## 8. Aspectos no decididos de la implementación

Esta versión no selecciona herramientas, algoritmos, lenguajes, bibliotecas ni
arquitectura. Tampoco fija un organismo, una ploidía, una tecnología de
secuenciación, las clases exactas de variantes de los individuos ni medidas de
evaluación.

Estos aspectos permanecen abiertos. No se consideran excluidos del proyecto.

## 9. Próximos pasos de diseño

1. Resolver si la evaluación cubrirá solo SNP o también otros tipos de variante.
2. Separar con precisión las alteraciones de la referencia y los cambios de los
   individuos.
3. Definir el resultado mínimo útil y el alcance de la evaluación que realizará
   el propio sistema.
4. Concretar los modelos biológicos, las entradas y el contenido de cada salida.
5. Completar los criterios de aceptación con medidas y casos de prueba.
6. Tomar las decisiones de implementación cuando el comportamiento esperado
   esté suficientemente definido.
