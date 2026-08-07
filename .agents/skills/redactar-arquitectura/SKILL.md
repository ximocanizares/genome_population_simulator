---
name: redactar-arquitectura
description: Convertir una especificación funcional ya acordada de un proyecto de análisis genético en una arquitectura técnica clara, trazable e independiente del lenguaje. Usar al definir o revisar componentes, responsabilidades, contratos, flujos de datos, invariantes, errores, reutilización y límites del sistema antes de planificar la implementación.
---

# Redactar arquitectura

## Objetivo

Explicar cómo se dividirá el sistema y cómo colaborarán sus partes para cumplir
la especificación, sin convertir el documento en código ni en un plan de tareas.

Escribir en castellano. Hacer comprensible el significado para un genetista y
dar a quien implemente contratos suficientemente precisos.

## Puerta de entrada

Trabajar a partir de una especificación vigente. Antes de diseñar:

1. leerla completa y localizar objetivos, no objetivos, invariantes, entradas,
   salidas y decisiones pendientes;
2. inspeccionar el repositorio para conocer componentes y convenciones reales;
3. distinguir lo que se reutiliza, modifica, crea o aplaza;
4. confirmar que no quedan preguntas funcionales que cambien la arquitectura.

Si falta una decisión funcional, devolverla a la especificación. No disfrazar
una preferencia técnica de requisito acordado.

Mientras no se haya elegido lenguaje, describir contratos, datos y dependencias
sin sintaxis, bibliotecas, tipos o herramientas propios de un lenguaje.

## Trazabilidad

Relacionar cada decisión arquitectónica con una sección concreta de la
especificación. Cuando dependa del código existente, citar rutas, módulos,
interfaces o comandos comprobados.

Para cada decisión importante indicar:

- requisito que satisface;
- opciones consideradas si hubo una elección real;
- motivo confirmado de la elección;
- consecuencias, restricciones y trabajo futuro.

No fabricar razones ni comparar alternativas que no se hayan evaluado.

## Contenido mínimo

### 1. Contexto y límites

Resumir qué problema resuelve esta arquitectura, qué queda dentro, qué queda
fuera y qué sistemas o formatos externos forman sus fronteras.

### 2. Modelo conceptual

Definir las entidades principales, su identidad, sus relaciones y su ciclo de
vida. En este proyecto pueden incluir referencias, individuos, poblaciones,
haplotipos, eventos, lecturas, alineamientos, genotipos y conjuntos de
resultados, pero usar solo las entidades aprobadas por la especificación.

### 3. Componentes y responsabilidades

Para cada componente describir:

- responsabilidad única;
- entradas y salidas;
- datos que posee y datos que solo consulta;
- contratos que ofrece y consume;
- errores que produce;
- dependencias permitidas y prohibidas;
- forma de comprobarlo de manera aislada.

Evitar componentes genéricos como “utilidades” o “gestor” si no expresan una
responsabilidad del dominio.

### 4. Flujos de datos y ejecución

Seguir cada recorrido importante desde la entrada hasta los resultados. Indicar
orden, transformaciones, estados persistidos, metadatos y puntos de fallo.

Usar un diagrama solo cuando aclare relaciones o secuencias que la prosa no
muestre bien. Acompañarlo con una explicación textual de sus reglas.

### 5. Contratos e invariantes

Definir qué debe cumplirse siempre en las fronteras entre componentes. Incluir,
cuando corresponda:

- identidad y versión de referencias y muestras;
- sistemas de coordenadas;
- ploidía, fase y procedencia de haplotipos;
- compatibilidad entre archivos y metadatos;
- semillas y reproducibilidad;
- reglas de ordenación, unicidad e integridad;
- comportamiento ante datos incompletos o incompatibles.

Dar una representación imposible o un error explícito a los estados inválidos;
no confiar únicamente en convenciones narrativas.

### 6. Configuración, persistencia y resultados

Separar configuración de entrada, estado interno, manifiestos y resultados.
Definir quién valida cada dato, cuándo queda fijado y cómo puede rastrearse hasta
la ejecución que lo produjo. No decidir todavía el formato físico si la
especificación lo aplaza hasta elegir lenguaje.

### 7. Errores y observabilidad

Clasificar fallos de configuración, incompatibilidades, límites de recursos y
fallos externos. Indicar qué se detiene, qué puede continuar, qué aviso recibe
el usuario y qué contexto debe conservarse para diagnosticarlo.

### 8. Estrategia de comprobación

Definir las fronteras que permiten pruebas aisladas y las verdades independientes
que validarán los recorridos completos. Prestar especial atención a errores que
podrían producir resultados biológicamente plausibles pero incorrectos.

### 9. Evolución prevista

Registrar extensiones futuras y señalar qué contrato actual evita bloquearlas,
sin construir abstracciones para casos que todavía no están definidos.

## Reglas de diseño

- Mantener una dirección clara de dependencias.
- Colocar los conceptos compartidos en fronteras neutrales, no dentro de una
  fase concreta del flujo.
- Separar lógica biológica, orquestación y entrada/salida cuando tengan ritmos
  de cambio o pruebas distintos.
- Hacer explícitas identidad, unidades, coordenadas y procedencia de los datos.
- Evitar duplicar reglas en varios componentes; asignar un único responsable.
- Preferir interfaces pequeñas basadas en necesidades observables.
- No optimizar sin un requisito y una medida; sí dejar visibles los límites de
  recursos y los puntos donde medir.
- No seleccionar dependencias o herramientas por costumbre. Registrar esas
  elecciones cuando se conozcan sus requisitos.

## Comprobación final

- Cada requisito arquitectónicamente relevante tiene un responsable.
- Cada componente tiene límites, contratos y errores claros.
- Los flujos principales y excepcionales llegan a un resultado definido.
- Las invariantes se conservan en todas las fronteras.
- No se mezclan datos asociados a referencias incompatibles.
- Las decisiones se apoyan en la especificación o en evidencia del repositorio.
- No quedan elecciones funcionales ocultas ni detalles prematuros del lenguaje.
- La arquitectura puede convertirse en un plan sin inventar diseño nuevo.

La arquitectura está lista cuando explica con precisión cómo se satisfará la
especificación y dónde se comprobará cada responsabilidad.
