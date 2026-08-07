---
name: revisar-rendimiento-rust
description: Revisar y medir el rendimiento de código Rust del simulador genómico, con atención a CPU, memoria, asignaciones, disposición de datos, entrada/salida y paralelismo. Usar ante lentitud, consumo elevado, falta de escalabilidad, recorridos críticos o cambios de optimización; exigir perfiles o pruebas reproducibles antes de afirmar mejoras.
---

# Revisar rendimiento Rust

## Objetivo

Localizar costes que importen en escenarios representativos y proponer
experimentos capaces de confirmar o refutar cada mejora. No tratar una forma de
código como rápida por intuición ni optimizar rutas frías.

Una revisión no autoriza modificaciones. Aplicar cambios solo si el usuario lo
pide después de conocer la evidencia y el coste de complejidad.

## Referencias

Leer siempre:

- [formato-hallazgos.md](references/formato-hallazgos.md);
- [metodologia.md](references/metodologia.md).

Leer según el problema:

- asignaciones, copias, estructuras y caché:
  [memoria-y-datos.md](references/memoria-y-datos.md);
- bucles, ramas, despacho, compresión y archivos:
  [bucles-e-io.md](references/bucles-e-io.md);
- hilos, bloqueos, canales y exceso de paralelismo:
  [concurrencia.md](references/concurrencia.md).

En una revisión general de rendimiento, leerlas todas.

## Principios

- Medir antes de prescribir.
- Comparar siempre la misma entrada, semilla, versión, compilación y número de
  hilos.
- Evaluar tiempo, memoria, asignaciones e I/O; mejorar uno puede empeorar otro.
- Diferenciar ruta crítica, patrón probable y posibilidad especulativa.
- No estimar porcentajes ni multiplicadores sin resultados reales.
- Incluir el coste de complejidad y mantenimiento de cada propuesta.
- Conservar exactamente la salida biológica y técnica, salvo que el experimento
  estudie expresamente otro compromiso.

## Procedimiento

### 1. Definir el escenario

Especificar tamaño de referencia, cromosomas, ploidía, individuos, generaciones,
variantes, cobertura, tecnología, eventos y hilos. Identificar qué métrica
importa: individuos por segundo, generaciones por segundo, bases o lecturas por
segundo, tiempo total, memoria máxima, archivos producidos o latencia de una
operación.

No usar solo ejemplos diminutos si el problema aparece con poblaciones grandes.
Mantener también un caso pequeño para detectar costes fijos.

### 2. Inventariar medidas existentes

Buscar pruebas de rendimiento, perfiles, registros de memoria y resultados
anteriores. Comprobar que se ejecutaron en modo optimizado y que el trabajo no
fue eliminado por el compilador.

### 3. Obtener una línea base

Ejecutar varias repeticiones, conservar datos brutos y registrar entorno,
compilador, características, entrada, semilla e hilos. Usar Criterion u otra
herramienta ya elegida para micropruebas; usar ejecuciones completas y un
perfilador para recorridos de extremo a extremo.

Elegir las herramientas disponibles en el entorno real. No suponer sistema
operativo, contenedor ni contadores de hardware. Si una herramienta no está
disponible, indicarlo y usar una alternativa menos detallada sin inventar datos.

### 4. Perfilar

Identificar primero dónde se consume el tiempo o la memoria. Relacionar las
muestras con funciones y frecuencia de llamada. Comprobar si el límite es CPU,
memoria, caché, asignaciones, compresión, disco, procesos externos o espera entre
hilos.

### 5. Realizar pasadas temáticas

Aplicar las referencias pertinentes. Si se permite usar agentes, repartir
categorías independientes y entregar a cada uno la misma línea base. La revisión
principal debe verificar y deduplicar los hallazgos. Sin autorización, realizar
las pasadas secuencialmente.

### 6. Proponer experimentos

Cada propuesta debe incluir:

- mecanismo de coste observado;
- cambio mínimo que lo eliminaría;
- medida que confirma o refuta la hipótesis;
- umbral acordado o criterio estadístico de la herramienta;
- complejidad y riesgo añadidos;
- comprobación de identidad de resultados.

No recomendar una reescritura amplia cuando un experimento local pueda probar
primero la causa.

### 7. Volver a medir

Comparar antes y después con idénticas condiciones. Repetir, conservar el ruido
y comprobar que la mejora no procede de cambiar la cantidad de trabajo. Medir
escalado con varios números de hilos cuando el paralelismo forme parte de la
hipótesis.

## Formato del informe

1. Alcance, entorno y escenario.
2. Veredicto y limitación principal.
3. Línea base y calidad de la evidencia.
4. Plan de medición pendiente.
5. Hallazgos ordenados por evidencia.
6. Observaciones fuera de alcance.
7. Aspectos ya eficientes.
8. Comandos y datos necesarios para repetir.

Guardar el informe solo cuando se solicite. No presentar un patrón observado en
el código como una mejora demostrada: etiquetarlo como hipótesis hasta medirlo.

## Cierre

- Todos los hallazgos citan código y evidencia leídos.
- Las rutas frías no aparecen como prioridades.
- Cada cambio propuesto tiene experimento y comprobación de identidad.
- Se distinguen tiempo, memoria, asignaciones e I/O.
- Se registra el efecto del número de hilos y de los procesos externos.
- Las conclusiones no exceden la calidad de la medida.
