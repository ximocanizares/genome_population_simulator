---
name: redactar-especificaciones
description: Redactar y revisar especificaciones de proyectos de análisis genético en castellano claro, trazables y comprensibles para genetistas. Usar al crear, ampliar, corregir o simplificar requisitos y decisiones funcionales, separar lo decidido de lo pendiente y mantener explícitos el alcance, las restricciones y el trabajo futuro sin elegir todavía la arquitectura ni el lenguaje de implementación.
---

# Redactar especificaciones

## Objetivo

Convertir las decisiones del proyecto en un contrato funcional que puedan leer
dos públicos: la persona experta en genética, que debe validar el significado
biológico, y quien vaya a diseñar o implementar el sistema, que necesita
requisitos inequívocos.

Escribir siempre en castellano. Mantener los nombres oficiales de formatos,
herramientas y conceptos cuando traducirlos reste precisión.

## Límites

- Describir qué debe hacer el sistema, para qué y bajo qué restricciones.
- No elegir arquitectura, lenguaje, bibliotecas ni estructura interna salvo que
  ya estén acordados y formen parte del requisito.
- No convertir la especificación en un plan de tareas.
- No inventar motivaciones, ventajas, cifras ni decisiones para rellenar huecos.
- No atribuir al código o a una herramienta un comportamiento que no se haya
  comprobado en una fuente disponible.

Si falta una decisión que cambia el comportamiento observable, marcarla como
pendiente y plantearla de forma explícita. Avanzar con las partes que no dependan
de ella.

## Fuentes y trazabilidad

Antes de redactar:

1. Leer la versión vigente de la especificación y los documentos relacionados.
2. Recoger las decisiones expresas del usuario y conservar su significado.
3. Inspeccionar el repositorio cuando se hagan afirmaciones sobre lo que ya
   existe; citar rutas, secciones o interfaces concretas.
4. Distinguir en las notas de trabajo entre hecho confirmado, propuesta,
   decisión pendiente y ampliación futura.
5. Resolver contradicciones con la decisión más reciente y dejar constancia si
   el cambio afecta a otras secciones.

No copiar literalmente las frases del usuario cuando una redacción más clara
mantenga exactamente el mismo requisito. Sí conservar los términos cuya forma
sea parte de un formato, una interfaz o una decisión.

## Preguntas que debe responder cada requisito

Comprobar, cuando sean pertinentes, estas cinco cuestiones:

1. **Qué** comportamiento o resultado se exige.
2. **Cómo se observa**, mediante entradas, salidas, estados o errores visibles.
3. **Por qué importa**, solo cuando la razón esté acordada o sea necesaria para
   interpretar correctamente la regla.
4. **Qué restricciones** y casos límite debe respetar.
5. **Qué queda fuera** de esta versión y dónde se registra para el futuro.

Evitar justificar una decisión con beneficios genéricos. Si la razón no consta,
escribir el requisito sin inventarla.

## Organización recomendada

Adaptar la estructura al documento, manteniendo como mínimo:

1. objetivo y contexto;
2. alcance incluido;
3. alcance excluido y ampliaciones futuras;
4. términos y modelo conceptual;
5. requisitos funcionales, agrupados por capacidad;
6. reglas e invariantes que nunca deben romperse;
7. entradas, salidas y metadatos observables;
8. errores y situaciones excepcionales;
9. configuración que debe poder expresar el usuario;
10. criterios de aceptación o ejemplos que permitan validar el significado;
11. decisiones pendientes.

Dar un destino explícito a todo requisito aplazado. No usar un cajón genérico de
“fuera de alcance” si ya se conoce la fase o documento que deberá recuperarlo.

## Reglas de redacción

- Pensar primero en un lector genetista. Dar por conocidos los conceptos
  habituales de genética y análisis de variantes, pero no la arquitectura de
  software, las pruebas ni la estadística avanzada.
- Usar frases cortas, voz activa y una idea principal por párrafo.
- Explicar primero para qué sirve un concepto técnico y después nombrarlo.
- Definir términos ambiguos y desarrollar las siglas no habituales en su
  primera aparición.
- Usar siempre el mismo nombre para el mismo concepto.
- Separar requisitos obligatorios, opciones configurables, valores
  predeterminados y posibilidades futuras.
- Preferir prosa para explicar decisiones y listas o tablas para conjuntos de
  requisitos comparables.
- Incluir ejemplos solo si aclaran una regla, una coordenada, una interacción o
  un caso límite. No introducir mediante el ejemplo un requisito nuevo.
- Mantener fórmulas, pseudocódigo y detalles estadísticos solo cuando ayuden a
  validar una decisión.
- Ser conciso: eliminar repetición, antecedentes innecesarios y secciones que
  no cambien la comprensión o la implementación posterior.

Crear un glosario breve únicamente cuando varios términos nuevos necesiten una
referencia común.

## Revisión de coherencia

Al modificar una decisión:

1. buscar todas sus apariciones y reglas dependientes;
2. actualizar ejemplos, valores predeterminados y casos límite afectados;
3. comprobar que no se mezclan entidades incompatibles, como datos asociados a
   referencias genómicas distintas;
4. conservar identificadores o numeración estables cuando sean usados por otros
   documentos;
5. anotar el cambio en el historial del documento si existe.

## Comprobación final

Antes de entregar:

- verificar que objetivos, no objetivos y versión prevista están claros;
- comprobar que cada obligación puede observarse o validarse;
- confirmar que lo decidido, lo configurable, lo pendiente y lo futuro se
  distinguen sin ambigüedad;
- revisar que no se han añadido decisiones de diseño o implementación;
- comprobar la coherencia entre secciones, ejemplos y términos;
- confirmar que las afirmaciones sobre archivos o código están fundamentadas;
- leer el documento desde el punto de vista de un genetista y de quien tendrá
  que diseñar el sistema;
- eliminar frases que no aporten una decisión, una restricción o una
  explicación necesaria.

La especificación está lista cuando ambos públicos pueden acordar qué debe hacer
el sistema sin tener que adivinar requisitos ni decisiones de implementación.
