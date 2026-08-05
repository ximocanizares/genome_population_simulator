---
name: redactar-especificaciones
description: Redactar y revisar especificaciones de proyectos de análisis genético en castellano claro para genetistas que no tienen por qué ser informáticos ni estadísticos. Usar al crear, ampliar, corregir o simplificar requisitos, diseños, planes técnicos y decisiones del proyecto, y al explicar conceptos de software, testing, formatos genómicos o estadística.
---

# Redactar especificaciones

## Objetivo

Redactar especificaciones de proyectos de análisis genético que sean
técnicamente correctas y que un genetista pueda entender sin necesitar
conocimientos avanzados de informática o estadística.

## Reglas obligatorias

1. Escribir en castellano.
2. Usar un lenguaje claro, directo y preciso.
3. Pensar primero en un lector genetista. Dar por conocidos los conceptos
   habituales de genética y variant calling, pero no asumir experiencia en
   arquitectura de software, programación, testing o estadística avanzada.
4. Definir un término la primera vez que aparezca si puede resultar técnico,
   ambiguo o confuso para ese lector.
5. Explicar primero qué significa o para qué sirve un concepto y después, si es
   útil, indicar su nombre técnico.
6. Desarrollar las siglas en su primera aparición, salvo las de uso habitual
   para la audiencia, como DNA, SNP, BAM o VCF.
7. Mantener la precisión técnica. Simplificar la explicación, no el requisito ni
   su significado.

## Forma de redactar

- Preferir frases cortas y voz activa.
- Presentar una idea principal por párrafo.
- Usar títulos que indiquen claramente el contenido de la sección.
- Usar listas para requisitos o pasos, sin convertir todo el documento en una
  sucesión de viñetas.
- Evitar anglicismos cuando exista una expresión castellana clara. Mantener el
  término inglés cuando sea el nombre reconocido de una herramienta, formato o
  concepto, y definirlo si hace falta.
- Evitar palabras abstractas o formales cuando pueda decirse lo mismo de manera
  directa.
- Separar la explicación funcional de los detalles de implementación.
- Indicar por qué importa una regla cuando su finalidad no resulte evidente.
- No añadir fórmulas, pseudocódigo o detalles estadísticos si no ayudan a tomar
  una decisión o comprender un requisito.

## Cómo definir términos técnicos

Al introducir un término técnico:

1. Explicar su función en lenguaje común.
2. Dar el nombre técnico entre paréntesis o en la misma frase.
3. Aclarar cómo se usa en el documento si puede tener varios significados.
4. No repetir la definición en cada aparición posterior.

Si una especificación contiene muchos términos nuevos, crear un glosario breve.
No crear un glosario para términos que ya estén bien explicados en el texto.

## Proceso de trabajo

1. Identificar el objetivo de la especificación y las decisiones que debe
   permitir tomar.
2. Identificar qué conocimientos puede darse por supuestos para un genetista.
3. Localizar términos de informática, testing o estadística que necesiten una
   definición.
4. Organizar el contenido desde la visión general hasta los detalles.
5. Redactar con las reglas anteriores.
6. Revisar la especificación completa desde el punto de vista del lector
   objetivo.

## Revisión final

Antes de entregar una especificación, comprobar:

- que está escrita en castellano;
- que un genetista puede entender el objetivo de cada sección;
- que no se presuponen conocimientos de informática o estadística avanzada;
- que los términos técnicos o confusos están definidos;
- que las siglas no habituales están desarrolladas;
- que no hay frases innecesariamente largas o formales;
- que la simplificación no ha cambiado el significado técnico;
- que los requisitos obligatorios se distinguen de recomendaciones y opciones;
- que la especificación mantiene una terminología coherente.

## Evolución de esta guía

Añadir a esta skill las nuevas recomendaciones de redacción acordadas durante el
proyecto. Integrarlas en la sección adecuada y evitar reglas duplicadas o
contradictorias.
