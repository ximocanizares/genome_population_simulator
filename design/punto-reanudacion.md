# Punto de reanudación del proyecto

- **Fecha:** 12 de agosto de 2026
- **Último commit publicado al iniciar la sesión:** `2eac575 docs: plan population SNV matrices`
- **Estado de esta revisión:** DTPD-118 se ha discutido y resuelto; sus cambios
  documentales deben validarse, confirmarse y publicarse antes de continuar con
  DTPD-119. Los directorios `tmp*/` están ignorados por Git.

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
- ACE-01 a ACE-117, DTC-001 a DTC-054, DTF-001 a DTF-039 y DTPD-101 a
  DTPD-119; DTPD-118 ya está resuelta.
- El plan F2 conserva seis pausas; su recuento de pasos se actualiza para
  materializar DTC-054 antes de este commit.
- Los 12 requisitos nuevos de F2 se reparten en 4 para F2.2 y 8 para F2.3.
- `git diff --check` estaba limpio.
- No se hizo commit de esta ampliación de matrices y espectros.

## Siguiente conversación

Revisar DTPD-119: formato canónico y muestreo del panel de haplotipos. La primera
propuesta aún no aprobada es un paquete con catálogo de sitios, alelos por
haplotipo, frecuencias por población y metadatos JSON. Después deberán cerrarse
el muestreo intacto y el modo recombinante. RAF2-001, RAF2-004 y RAF2-007 siguen
pendientes de revisión.
