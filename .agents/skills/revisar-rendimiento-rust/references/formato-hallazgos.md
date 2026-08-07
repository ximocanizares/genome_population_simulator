# Formato y evidencia de rendimiento

## Niveles

- **Ruta crítica:** un perfil o prueba reproducible identifica el lugar y el
  mecanismo explica una parte relevante del coste. Requiere confianza alta.
- **Probable:** el patrón y la frecuencia observada hacen plausible la mejora,
  pero todavía falta medir el cambio. Proponer un experimento antes de aplicar.
- **Especulativo:** el patrón existe, pero no se sabe si importa. No modificar;
  registrar solo cuando ayude a orientar una medición.
- **Nota:** enseñanza o posibilidad en código frío. No presentarla como trabajo
  prioritario.

Una coincidencia de patrón nunca alcanza por sí sola “Ruta crítica”.

## Formato

```markdown
- `ruta/al/archivo.rs:LÍNEA` — **[Nivel]** Título
  - Confianza: Alta|Media|Baja.
  - Evidencia: salida citada o “solo inspección del patrón”.
  - Frecuencia: unidad y número observado o desconocido.
  - Mecanismo: asignación, caché, bloqueo, syscall, compresión, rama, etc.
  - Experimento: comando, entrada, métrica y criterio de aceptación.
  - Coste de complejidad: tipos, dependencia, memoria, mantenimiento o riesgo.
  - Cambio candidato: modificación mínima que se probaría.
```

Agrupar repeticiones del mismo patrón en un único hallazgo con ubicaciones
representativas.

## Integridad

- Citar solo perfiles y pruebas realmente ejecutados.
- No calcular aceleraciones a partir de intuición.
- No llamar caliente a código de inicio, configuración o error sin una medida.
- No confundir menos instrucciones con menos tiempo total.
- No ocultar regresiones de memoria, exactitud, determinismo o legibilidad.
- Mantener la salida idéntica y comprobarla con un oráculo adecuado.
