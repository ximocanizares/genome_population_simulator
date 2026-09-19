# Formato y gravedad de los hallazgos

## Gravedad

- **Bloqueante:** puede corromper datos, producir resultados biológicos
  incorrectos, romper seguridad de memoria o impedir el uso principal. Exige
  corrección antes de continuar.
- **Grave:** defecto de corrección, pérdida de información, contrato público
  roto o ausencia de una prueba esencial. Debe resolverse antes de publicar.
- **Menor:** problema real de mantenibilidad, diagnóstico, robustez o cobertura
  con impacto limitado.
- **Nota:** mejora pequeña y concreta. Agrupar las notas relacionadas y evitar
  ruido estilístico.

No elevar gravedad por el tamaño del diff. Basarla en impacto y alcance.

## Confianza

- **Alta:** el código y una prueba, contrato o ejecución demuestran el problema.
- **Media:** el mecanismo es claro, pero falta una condición del entorno o un
  caso ejecutable.
- **Baja:** existe una posibilidad sin evidencia suficiente. No presentarla como
  hallazgo; convertirla en pregunta o comprobación propuesta.

Un hallazgo bloqueante necesita confianza alta o una explicación explícita de
la incertidumbre que por sí misma impide validar el sistema.

## Formato

```markdown
- `ruta/al/archivo.rs:LÍNEA` — **[Gravedad]** Título
  - Confianza: Alta|Media.
  - Contrato: requisito, invariante o comportamiento esperado.
  - Evidencia: código, prueba o salida observada.
  - Mecanismo: cómo aparece el fallo y qué datos afecta.
  - Corrección: cambio mínimo o decisión necesaria.
  - Prueba: caso que fallaría antes y pasaría después.
```

Usar una ubicación primaria por hallazgo. Añadir otras ubicaciones en la
explicación si comparten la misma causa.

## Reglas de integridad

- Leer la línea y su contexto antes de citarla.
- No inventar comandos, salidas, perfiles, firmas ni versiones.
- No afirmar que una prueba pasa si no se ejecutó.
- No emitir un hallazgo por una preferencia sin impacto demostrable.
- No duplicar una causa en varias categorías.
- Distinguir defecto introducido, defecto preexistente y observación fuera de
  alcance.
- Si no hay hallazgos, escribirlo de forma explícita y resumir qué se comprobó.

## Respuesta del autor

Permitir que cada hallazgo termine como: aceptado, aceptado con adaptación,
aplazado con destino, disputado con evidencia o ya corregido. No considerar
“ignorado” un estado válido.
