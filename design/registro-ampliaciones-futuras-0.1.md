# Registro de ampliaciones funcionales futuras 0.1

## 1. Finalidad y uso

Este es el registro único de capacidades funcionales excluidas de v1. Conserva
los identificadores DTF del registro técnico, pero añade fuente, destino y
dependencias recuperables. No autoriza ninguna capacidad: incorporarla exige
una decisión funcional posterior y la actualización de los REQ normativos.

Las elecciones de representación o implementación que no cambian entradas,
resultados ni errores pertenecen a DTP/DTPD en
[decisiones-tecnicas-0.1.md](./decisiones-tecnicas-0.1.md), no a este registro.

## 2. Registro canónico

| ID | Ampliación excluida de v1 | Fuente normativa o sección | Motivo y estado | Destino | Dependencias mínimas |
| --- | --- | --- | --- | --- | --- |
| DTF-001 | Oxford Nanopore y PacBio CLR | REQ-CONF-16; REQ-TEC-04 | Fuera de alcance; abierta post-v1 | Post-v1, decisión de tecnologías | Perfiles, errores, calidad y validación |
| DTF-002 | Paneles, amplicones, GBS y otros protocolos | REQ-TEC-08; REQ-WGS-12 | v1 solo WGS; abierta post-v1 | Post-v1, decisión de protocolos | Objetivos, regiones capturadas, cobertura y salidas |
| DTF-003 | Ploidía variable por individuo, cromosomas sexuales y aneuploidía | REQ-IND-07/08; REQ-CRU-23 | v1 exige ploidía común por secuencia; abierta post-v1 | Post-v1, decisión de ploidía variable | Genotipo, reproducción, dosis, VCF y validación |
| DTF-004 | Reproducción con ploidías impares, multivalentes y doble reducción | REQ-CRU-23..26; REQ-REC-03..08 | Fuera de v1; abierta post-v1 | Post-v1, decisión de meiosis poliploide | DTF-003, emparejamiento y segregación |
| DTF-005 | Conversión génica e interferencia | REQ-REC-01/02/09 | v1 usa entrecruzamientos sin interferencia; abierta post-v1 | Post-v1, decisión de recombinación avanzada | Mapas, eventos de transmisión y procedencia |
| DTF-006 | Selección poligénica, ambiental, estabilizadora o disruptiva | REQ-SEL-01..17 | v1 limita selección a un SNP y modelo direccional; abierta post-v1 | Post-v1, decisión de selección avanzada | DTF-019, caracteres, aptitud y diagnóstico |
| DTF-007 | Modelado específico de elementos transponibles | REQ-TE-01..14 retirados | Clase retirada de v1; abierta post-v1 | Post-v1, decisión completa de TE | Catálogo, mecanismos, destinos, identidad, verdad y metadatos |
| DTF-008 | Correspondencia de variantes entre referencias | REQ-VER-14/15; REQ-MAN-12 | Prohibida en v1; abierta post-v1 | Post-v1, decisión de comparación entre referencias | Coordenadas, normalización y linaje de referencias |
| DTF-009 | Reproducibilidad entre versiones | REQ-CONF-111..115 | No garantizada en v1; abierta post-v1 | Post-v1, política de compatibilidad | Versiones, algoritmos, perfiles y dependencias |
| DTF-010 | Reanudación distribuida y ejecución remota | Sección 9 de la especificación; DTP-007 | Sin contrato funcional en v1 | Post-v1, decisión de ejecución distribuida | Atomicidad, persistencia, seguridad y recursos |
| DTF-011 | Reproducción entre referencias estructurales distintas | REQ-IND-53; REQ-CRU-44; REQ-MIG-09 | Rechazada en v1; abierta post-v1 | Post-v1, decisión de compatibilidad entre referencias | DTF-021, identidad de loci, recombinación y descendencia |
| DTF-012 | Carriles y cardinalidades técnicas mayores | REQ-TEC-01..04; REQ-RG-01..07 | v1 no modela carriles; abierta post-v1 | Post-v1, decisión de organización técnica | Bibliotecas, tandas, FASTQ, RG, contaminación y manifiestos |
| DTF-013 | Sustituciones para descendientes futuros | REQ-CONF-17/18/120 | Rechazadas en v1; abierta post-v1 | Post-v1, decisión de identidad futura | Preasignación, momento de aplicación, herencia y precedencia |
| DTF-014 | Distribuciones adicionales de recuentos aleatorios | REQ-CONF-25; REQ-IND-18; REQ-VAL-48 | v1 solo Poisson; abierta post-v1 | Post-v1, decisión de distribuciones | Parámetros, ámbitos, validación y reproducibilidad |
| DTF-015 | Repeticiones con saltos de varias unidades | REQ-CONF-61/62; REQ-MUT-11 | v1 solo paso unitario; abierta post-v1 | Post-v1, decisión de modelos de repetición | Distribución del salto, límites, verdad y manifiesto |
| DTF-016 | Calibración empírica y reutilización de perfiles | REQ-CONF-85..87; REQ-DAT-41..45 retirados | Rechazada en v1; abierta post-v1 | Post-v1, decisión de calibración | Entradas, inferencia, fiabilidad, esquema y reutilización |
| DTF-017 | Umbrales y consecuencias de validación estadística | REQ-VAL-10..13 | v1 solo informa desviaciones; abierta post-v1 | Post-v1, decisión de aceptación estadística | Métricas, tamaños muestrales, intervalos y severidad |
| DTF-018 | Ne objetivo y control de contribución parental | REQ-REP-07 retirado; REQ-REP-11..16 | Rechazado en v1; abierta post-v1 | Post-v1, decisión de Ne objetivo | Ploidía, N, contribuciones, fórmula y optimización |
| DTF-019 | Selección sobre variantes distintas del SNP bialélico único | REQ-SEL-01..17 | Fuera de v1; abierta post-v1 | Post-v1, decisión de locus seleccionado general | Identidad, dosis, protección y número de copias |
| DTF-020 | F automático o realizado fuera de diploides | REQ-EVO-33/34/46..52 | Fuera de v1; abierta post-v1 | Post-v1, decisión de endogamia poliploide | Herencia polisómica/disómica y estimadores |
| DTF-021 | Reproducción con diferencias estructurales | REQ-CRU-23..31; REQ-REC-03..08 | Incompatible en v1; abierta post-v1 | Post-v1, decisión de compatibilidad estructural | Emparejamiento, gametos, fertilidad y colinealidad |
| DTF-022 | Órdenes generacionales alternativos | REQ-REP-24/25; REQ-CONF-44/45 | Rechazados en v1; abierta post-v1 | Post-v1, decisión de orden configurable | Catálogo o grafo, validación y trazabilidad |
| DTF-023 | Niveles configurables o degradación de trazabilidad | REQ-SAL-02/16; REQ-NOF-05/06 | Prohibidos en v1; abierta post-v1 | Post-v1, decisión de niveles de trazabilidad | Estados completos, resúmenes, presupuesto y auditoría |
| DTF-024 | Llamada, imputación, faseado inferido y evaluación de herramientas externas | Especificación §3, párrafos 1–3 | Pertenece a otro proyecto, no a v1 | Otro proyecto o futura decisión expresa de alcance | VCF de verdad, métricas y contratos de herramientas externas |
| DTF-025 | Señales crudas, BCL y lecturas detalladas de índices | Especificación §3, párrafo sobre datos del instrumento | Fuera de v1; abierta post-v1 | Post-v1, decisión de simulación del instrumento | Tecnología, señales, índices y demultiplexado |
| DTF-026 | Quimeras entre individuos o bibliotecas | Especificación §3; REQ-CONF-99 | Prohibidas en v1; abierta post-v1 | Post-v1, decisión de quimeras entre fuentes | Identidad, procedencia, contaminación y bibliotecas |
| DTF-027 | Desequilibrio alélico de biblioteca y pérdida alélica completa | Especificación §3; REQ-WGS-10 | Fuera de v1; abierta post-v1 | Post-v1, decisión de sesgo alélico | Dosis, muestreo WGS, cobertura y verdad de lecturas |
| DTF-028 | Generaciones solapadas, edades, supervivencia y cohortes | Especificación §3; REQ-REP-17..20 | v1 usa generaciones discretas no solapadas | Post-v1, decisión de ciclo vital | Identidad temporal, reproducción, N y genealogía |
| DTF-029 | Selección por fecundidad, elección de pareja y mutación dependiente del contexto | Especificación §3; REQ-SEL-01..17; REQ-MUT-01..15 | Modelos excluidos de v1 | Post-v1, decisiones separadas de evolución avanzada | Aptitud, apareamiento, contexto de secuencia y tasas |
| DTF-030 | Componente ambiental y heredabilidad del carácter | Especificación §3; REQ-SEL-06..10 | v1 iguala carácter y valor genético | Post-v1, decisión de carácter cuantitativo | Ambiente, varianza fenotípica, heredabilidad y selección |
| DTF-031 | Distorsión de segregación e impulso meiótico | Especificación §3; REQ-REC-01..29 | Fuera de v1; abierta post-v1 | Post-v1, decisión de transmisión no mendeliana | Segregación, recombinación, probabilidades y verdad |
| DTF-032 | Reproducción clonal, apomixis, dobles haploides y herencia parental asimétrica | Especificación §3; REQ-CRU-01..45; REQ-REP-01..26 | Fuera de v1; abierta post-v1 | Post-v1, decisión de modos reproductivos | Sexo parental, citoplasma, ploidía y genealogía |
| DTF-033 | Composición de variantes estructurales anidadas | Especificación §3; REQ-IND-49..52 | Excluida de v1 | Post-v1, decisión de composición estructural | Coordenadas, estado intermedio, linaje y VCF/complemento |
| DTF-034 | Variantes terminales, secuencias completas y copia dentro de su fuente | Especificación §3; REQ-IND-49..52 | Casos límite excluidos de v1 | Post-v1, decisión de límites estructurales | Coordenadas, cromosomas completos, ploidía y linaje |
| DTF-035 | Inferencia o generación de LD desde D, D' o r² | Especificación §3; REQ-EVO-40..45 | v1 recibe haplotipos y frecuencias explícitos | Post-v1, decisión de estadísticos de LD | Modelo poblacional, estimadores, bloques y tolerancias |
| DTF-036 | Posiciones invariantes y bloques gVCF | Especificación §3; REQ-VER-01..18 | VCF de verdad v1 solo representa diferencias | Post-v1, decisión de vista gVCF | Estándar seleccionado, bloques, genotipos y tamaño de salida |
| DTF-037 | Lotes reproducibles de ejecuciones completas | Especificación §3, párrafo de lotes | Excluidos de F1.6/F2.0 y de v1; cada réplica futura será una ejecución completa e independiente, no otra referencia dentro de la misma ejecución | Post-v1, después de F9.7 y de una decisión funcional específica de lotes | Vías F3–F8 y consolidación F9; identidad/raíz/semillas, estabilidad de prefijo, directorios/manifiestos y estados por réplica, identidad/resumen del lote, preflight agregado y paridad Marimo |

## 3. Cobertura y mantenimiento

- Toda frase normativa «fuera de v1», «futura», «aplazada» o equivalente debe
  estar cubierta por la fuente de una fila DTF. La fila enlaza todas las
  menciones del mismo alcance aunque el texto histórico use otra formulación.
- Toda elección pendiente puramente técnica debe usar DTP o DTPD y no DTF.
- Una ampliación nueva exige un DTF estable en el mismo cambio; no se admiten
  cajones genéricos sin destino y dependencias.
- Cerrar un DTF exige una nueva especificación funcional; nunca se cierra por
  una decisión de implementación.
- La revisión mecánica debe buscar las construcciones de futuro en `design/*.md`
  y comprobar que pertenecen a una fuente DTF o a una decisión DTP/DTPD.

## 4. Auditoría de expresiones diferidas

La auditoría del 8 de agosto de 2026 buscó, sin distinguir mayúsculas, las
formas `futuro`, `futura(s)`, `pendiente(s)`, `aplazado/a(s)`, `fuera de v1` y
`se decidirá`. Se contó una aparición por línea y se revisó su sección y su
significado, no solo la palabra encontrada. El recuento excluye esta propia
sección de auditoría para evitar que la tabla y su explicación se autocuenten.

| Documento | DTF | DTP | DTPD | Histórica o diagnóstica | Total |
| --- | ---: | ---: | ---: | ---: | ---: |
| especificaciones-0.1.md | 113 | 0 | 24 | 0 | 137 |
| arquitectura-0.1.md | 13 | 0 | 0 | 29 | 42 |
| decisiones-tecnicas-0.1.md | 2 | 4 | 12 | 15 | 33 |
| hoja-de-ruta-0.1.md | 0 | 12 | 0 | 0 | 12 |
| registro-ampliaciones-futuras-0.1.md | 14 | 0 | 0 | 0 | 14 |
| revision-adversarial-especificaciones-0.1.md | 0 | 0 | 0 | 39 | 39 |
| trazabilidad-requisitos-0.1.md | 0 | 0 | 0 | 14 | 14 |
| **Total** | **142** | **16** | **36** | **97** | **291** |

Las menciones históricas o diagnósticas describen el estado de una auditoría,
una prueba planificada o la ausencia de una decisión; no crean alcance futuro.
Las exclusiones funcionales normativas se enlazan mediante el DTF explícito de
la sección 3 o mediante los REQ estables citados en la fila DTF. La auditoría
detectó trece grupos que no estaban en DTF-001..023 y añadió DTF-024..036. Tras
esa corrección quedaron **0 exclusiones funcionales futuras huérfanas** en el
corte auditado. DTF-037 se añadió el 10 de agosto de 2026 por una petición
posterior y queda fuera de los recuentos históricos de la tabla.
