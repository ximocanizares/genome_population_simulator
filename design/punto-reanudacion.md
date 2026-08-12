# Punto de reanudación del proyecto

- **Fecha:** 10 de agosto de 2026
- **Último commit publicado:** `9df8b7c docs: plan reduced references and Marimo UI`
- **Estado posterior al commit:** existen cambios documentales sin confirmar en
  especificaciones, arquitectura, decisiones técnicas, hoja de ruta, plan F2 y
  trazabilidad. El directorio `tmp-genome-TDvk2K/` es intencionado y no debe
  añadirse ni eliminarse sin una petición expresa.

## Decisión que se estaba desarrollando

F2 debe admitir o generar una matriz de frecuencias SNV común a una población de
origen estadística. La matriz generada usa número exacto o media/densidad de loci,
distribución espacial, espectro de sustituciones y un espectro de frecuencias
alélicas explícito. Debe existir al menos un ejemplo empírico público y versionado
y también se puede suministrar un espectro propio; ningún perfil se aplica
implícitamente.

Los individuos se muestrean desde esa matriz. Hardy-Weinberg es el modo
predeterminado y el coeficiente de endogamia `F` solo se admite en diploides. La
matriz contiene frecuencias marginales y no crea LD. F6.4 admitirá un panel de
haplotipos y frecuencias en modo intacto; F7.1 añadirá la creación recombinante
reutilizando el motor común y conservando la procedencia por segmentos.

## Estado comprobado

- 843 requisitos y 843 filas de trazabilidad, sin ausencias ni duplicados.
- ACE-01 a ACE-117, DTC-001 a DTC-053 y DTPD-101 a DTPD-119.
- Plan F2: 77 pasos abiertos, 6 pausas, 84 requisitos exclusivos y 4 compartidos.
- Los 12 requisitos nuevos de F2 se reparten en 4 para F2.2 y 8 para F2.3.
- `git diff --check` estaba limpio.
- No se hizo commit de esta ampliación de matrices y espectros.

## Siguiente conversación

Revisar DTPD-118: decidir el formato exacto del espectro suministrado —valores
discretos, intervalos o ambos— y la regla para seleccionar una frecuencia dentro
de un intervalo. Después se revisará DTPD-119, el formato y muestreo del panel de
haplotipos. RAF2-001, RAF2-004 y RAF2-007 siguen pendientes de revisión.
