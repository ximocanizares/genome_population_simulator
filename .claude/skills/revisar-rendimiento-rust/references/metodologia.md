# Metodología y configuración de compilación

## Línea base

- Definir una entrada versionada y una semilla fija.
- Registrar CPU, memoria, sistema, compilador, perfil, características y número
  de hilos.
- Compilar en modo optimizado. No extraer conclusiones de rendimiento de una
  compilación de desarrollo.
- Ejecutar calentamiento y varias repeticiones; conservar resultados brutos y
  dispersión.
- Usar tamaños pequeño, representativo y grande para distinguir coste fijo,
  complejidad y límite de memoria.

## Medidas

- Medir tiempo total y, cuando proceda, CPU, memoria máxima, asignaciones, bytes
  leídos y escritos, compresión, espera y cambios de contexto.
- Para el simulador, informar unidades de dominio: individuos, generaciones,
  segmentos, bases, lecturas y registros por segundo.
- Medir por separado procesos externos como minimap2 para no atribuir su coste
  al núcleo Rust.
- Comprobar siempre que todas las ejecuciones realizan el mismo trabajo y
  producen resultados equivalentes.

## Herramientas

- Usar Criterion para funciones pequeñas y estables.
- Usar mediciones de extremo a extremo para escenarios completos.
- Elegir perfilador de CPU, memoria o asignaciones según el sistema disponible.
- No asumir que existen contadores de hardware, permisos de perfilado o un
  contenedor concreto.
- No instalar herramientas ni cambiar el sistema sin autorización.

## Experimentos

- Cambiar una causa cada vez.
- Comparar contra una línea base guardada antes de editar.
- Definir previamente qué resultado hace que la complejidad merezca la pena.
- Repetir la comparación con varios números de hilos si afecta al paralelismo.
- Revertir una optimización que no mejora de forma reproducible o perjudica un
  recurso más importante.

## Compilación

- Revisar perfiles de Cargo, optimización, información de depuración, LTO,
  unidades de generación y características solo mediante experimentos.
- Tratar opciones específicas de CPU como una decisión de distribución: pueden
  impedir ejecutar el binario en otras máquinas.
- Comparar asignadores, bibliotecas de compresión o características opcionales
  mediante pruebas A/B; no asumir un ganador universal.
- Registrar todas las opciones necesarias para repetir el resultado.

## Salida

Separar hechos medidos, inferencias del perfil e hipótesis pendientes. Una
medida sin el comando, la entrada y el entorno no constituye una línea base
reproducible.
