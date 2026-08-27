# Plan de implementación de F3 0.1

## 1. Cabecera y estado

- **Título:** generación mínima de lecturas desde haplotipos.
- **Versión:** 0.1.
- **Fecha:** 27 de agosto de 2026.
- **Estado:** F3.1.1–F3.1.2 completados y revisados; F3.1.3 es el siguiente
  paso.
- **Ámbito:** F3.1 y F3.5 simplificados de la hoja de ruta vigente.
- **Fuentes:** especificación §4.4 y §4.6–§4.8; arquitectura §5.5 y §6; DTC-005, DTC-013, DTC-016, DTC-018, DTC-024, DTC-031–032, DTC-035, DTC-037–045 y DTC-079; hoja de ruta F3 y trazabilidad vigente.

F3 convierte el estado haplotípico faseado que F2 mantiene durante la ejecución
en FASTQ mínimo para su análisis externo, incluido el mapeo que la persona
usuaria quiera realizar con minimap2 u otra herramienta. `truth.jsonl` es la
evidencia independiente publicada de esa verdad, pero no se relee ni se
materializa completo como fuente operativa de cada lectura. El consumidor real
es un único evento técnico de lecturas cortas Illumina de extremo único y sin
contaminación ni degradaciones avanzadas. Reutiliza referencias, identidades,
aleatoriedad, verdad y publicación existentes: no crea otra vía biológica ni un
mapeador.

## 2. Alcance y exclusiones

**Incluye** lecturas WGS mínimas desde haplotipos faseados; exactamente uno
entre cobertura o número de lecturas; reparto por dosis y ploidía; N muestreable
con Phred 0; un evento técnico explícito y su subsemilla; FASTQ, comprobación
interna por flujo y ensamblado mínimo. Cada lectura conserva durante su
procesamiento ID técnico, individuo, copia/haplotipo, secuencia, intervalo y
orientación; se descarta tras escribirla y comprobarla. Al cierre solo quedan
contadores agregados y los artefactos con sus hashes. Las pruebas usan verdad
sintética independiente para proteger origen, dosis/ploidía, cobertura, N e
identidades.

**Excluye:** minimap2, otros mapeadores, SAM/BAM/CRAM e índices; paired-end, otros perfiles, PacBio, contaminación, errores, duplicados, adaptadores, sesgo GC y quimeras; múltiples referencias/eventos y sus interacciones; modo directo; procedencia durable por lectura; paquete, OCI, hardening y plataforma. El mapeo y alineamientos se realizan fuera del simulador. F3 no crea escenarios, ejecutor de evidencia, protocolos de capacidades ni sesión persistente.

El orden es consumidor de verdad, FASTQ por flujo, validación y ensamblado. No
elige crates, formatos, perfiles cuantitativos validados ni herramientas
externas nuevas.

## 3. Dependencias y bloqueos

- F2 debe estar cerrado y proporcionar referencia válida, estado haplotípico
  directo faseado durante la ejecución, `truth.jsonl`, VCF de verdad e
  identidades, subsemillas y publicación reutilizables.
- F3 adopta únicamente el contrato sintético `illumina-se-minimal-v1`: lectura
  single-end de 150 bases, sin errores ni artefactos; Phred 40 para A/C/G/T y
  Phred 0 para N. No es un perfil cuantitativo validado y F4 conserva esa
  responsabilidad.
- No existe cobertura predeterminada: cada evento declara exactamente una
  cobertura WGS o un número de lecturas. La cobertura se expresa respecto a la
  longitud haploide elegible; su conversión a número de lecturas y el reparto
  entre copias son reproducibles y estocásticos. Las secuencias sin un intervalo
  elegible de 150 bases se excluyen del sorteo y se diagnostican si no queda
  ninguna.

## 4. Comprobaciones comunes

Cada hito ejecutado conserva una prueba focal y las puertas aprobadas:

    cargo fmt --all -- --check
    cargo clippy --workspace --all-targets --all-features --locked -- -D warnings
    cargo test --workspace --all-targets --all-features --locked
    cargo test --workspace --doc --all-features --locked
    cargo doc --workspace --no-deps --all-features --locked
    git diff --check

F3 no invoca validadores de BAM/CRAM, índices, minimap2 ni un SNP caller. La verdad independiente se crea antes de FASTQ.

## 5. F3.1 — Lecturas desde haplotipos

**Resultado:** un evento mínimo genera FASTQ final single-end desde copias faseadas, con identidad y origen suficientes para comprobarlo y exportarlo.

- ✅ F3.1.1. Admitir el único evento técnico mínimo y resolver valores efectivos sin escenarios ni capacidades anticipadas.
  - Depende de: F2 cerrada y contrato `illumina-se-minimal-v1` fijado en §3.
  - Procede de: REQ-CONF-64/69/101–103, REQ-DAT-15–18 y DTC-005/016/034/079.
  - Cambia: admisión TOML, configuración efectiva e identidad para un evento
    real; exige exactamente uno entre cobertura y número de lecturas, y registra
    longitud y calidad efectivas del contrato mínimo. Rechaza contaminación,
    pares, referencias múltiples y campos F4.
  - Se comprueba con: cargo test --test f3_event_admission --locked; fixture
    válido produce valores efectivos estables y fixtures rechazan un ID que
    colisione con un evento directo, doble o nulo objetivo WGS y campos F4.
    F3 usa una tabla única `[technical_event]`, por lo que no existe una lista
    técnica donde pueda aparecer un ID duplicado; tecnología y perfil se fijan
    por contrato y no son campos de autoría. puertas comunes.

- ✅ F3.1.2. Derivar ámbito aleatorio e identidades técnicas mínimas desde causas estables, no desde orden de iteración.
  - Depende de: F3.1.1 y catálogo de subsemillas F2 vigente.
  - Procede de: REQ-CONF-26–27, REQ-DAT-18/68, REQ-RG-03, REQ-MAN-72 y DTC-016/037–039.
  - Cambia: catálogo versionado y derivación solo para evento, biblioteca, tanda, grupo de lectura y plantillas que F3 genera. La frontera del evento recibe `ReferenceId` tipado. Batch, biblioteca y grupo de lectura conservan BLAKE3 y registro de colisiones; por excepción explícita de DTC-039, `TemplateId` se compone como `TPL_<library_id>_<template_ordinal_decimal>`, único por construcción sin retener un registro por plantilla.
  - Se comprueba con: cargo test --test f3_technical_identity --locked; semilla fija reproduce ID/subsemilla y el vector histórico `technical_event`; cambiar cada causa altera su ámbito, texto/ruta inválido no cruza la frontera tipada, añadir o reordenar individuo no cambia plantilla existente y derivar muchos ordinales no incrementa estado residente; puertas comunes.

- ☐ F3.1.3. Resolver el número estocástico de lecturas y muestrear sus orígenes
  desde la vista haplotípica faseada, sin reconstruir `truth.jsonl`.
  - Depende de: F3.1.1–F3.1.2 y estado haplotípico F2.
  - Procede de: REQ-DAT-15/24/68, REQ-WGS-01/04/07/09, DTC-018/031–032 y arquitectura §7.3.
  - Cambia: simulación WGS mínima de solo lectura del estado haplotípico; usa
    longitud haploide elegible como denominador y reparte por copias en
    proporción a esa longitud. Por lectura produce transitoriamente plantilla,
    individuo, copia/haplotipo, intervalo, orientación y secuencia verdadera;
    no serializa procedencia F4.
  - Se comprueba con: cargo test --test f3_wgs_sampling --locked; verdad
    sintética independiente verifica dosis/ploidía, denominador haploide,
    exclusión de extremos y secuencias cortas, semilla repetible, cobertura cero,
    cobertura estocástica e invariancia ante referencias de mapeo; puertas
    comunes.

- ☐ F3.1.4. Materializar y comprobar FASTQ final por flujo, conservando N y
  Phred 0, sin degradaciones ni acumulación de lecturas.
  - Depende de: F3.1.3.
  - Procede de: REQ-DAT-17–18/72–73, REQ-RG-06, REQ-SAL-08–09, DTC-007/031–032 y arquitectura §6.3.
  - Cambia: escritor temporal y comprobación interna de un grupo por
    individuo/evento; emite, comprueba y descarta una lectura cada vez. Conserva
    solo contadores agregados y hashes; nunca reconstruye FASTQ desde
    alineamientos.
  - Se comprueba con: cargo test --test f3_fastq_output --locked; fixture con N
    interior/terminal confirma base literal, Phred 0, número/ID, uso constante
    de memoria y correspondencia con verdad sintética; puertas comunes.

> **PAUSA-F3.1.** Confirmar que FASTQ procede de haplotipos faseados, cobertura no incorpora ploidía y no hay mapeador, perfiles, degradaciones, escenarios o persistencia general.

## 6. F3.5 — Conjunto exportable mínimo

**Resultado:** un directorio de evento real reúne referencia, verdad, FASTQ, configuración efectiva, comprobaciones internas y manifiesto.

- ☐ F3.5.1. Validar FASTQ por flujo desde vistas de solo lectura de verdad y
  reproducción determinista de los orígenes, sin procedencia durable por lectura.
  - Depende de: PAUSA-F3.1.
  - Procede de: REQ-VAL-01–05/10–12, REQ-DAT-72–73, REQ-SAL-09/13, DTC-004/024/031–032 y arquitectura §6.3.
  - Cambia: relee el FASTQ y reproduce el muestreo desde la subsemilla y el
    estado haplotípico para comprobar ID, origen, secuencia, calidad y N una
    lectura cada vez. Informa cobertura sin tolerancias ni SNP caller, mapeador
    o alineamiento.
  - Se comprueba con: cargo test --test f3_fastq_validation --locked; verdad
    sintética distinta del productor detecta origen alterado y conserva como
    informe, no fallo, desviación aleatoria de cobertura; puertas comunes.

- ☐ F3.5.2. Definir, validar y escribir el primer manifiesto JSON Schema de un
  evento técnico real.
  - Depende de: F3.5.1.
  - Procede de: REQ-MAN-01/72–73/83, REQ-SAL-04–05/17–18/21–25, DTC-005–006/013/035/045 y arquitectura §5.5/§7.4.
  - Cambia: esquema y escritor de `events/<event_id>/manifest.json` para un
    evento: ID/tipo, perfil efectivo, individuos, subsemilla, referencia,
    artefactos, hashes y estado. El manifiesto general lo enlaza; no introduce
    conjuntos genéricos, sesiones ni archivos de procedencia F4.
  - Se comprueba con: cargo test --test f3_event_manifest --locked; casos
    completo e inválido pasan o fallan frente al esquema propio y al validador
    JSON Schema existente; puertas comunes.

- ☐ F3.5.3. Ensamblar y publicar el primer evento real desde sus artefactos
  validados.
  - Depende de: F3.5.2.
  - Procede de: REQ-MAN-01/72–73/83, REQ-SAL-04–05/17–18/21–25,
    DTC-005–006/013/035/045 y arquitectura §5.5/§7.4.
  - Cambia: directorio de evento con rutas relativas, estados calculados desde
    artefactos validados y enlace desde el manifiesto general; no conjuntos
    vacíos, archivo exportable ni sesiones.
  - Se comprueba con: cargo test --test f3_event_assembly --locked; caso
    completo y fallo controlado prueban que solo artefactos validados son
    completos e independientes conservan estado; puertas comunes.

- ☐ F3.5.4. Cerrar recorrido extremo a extremo con verdad conocida pequeña y FASTQ verificable.
  - Depende de: F3.5.3.
  - Procede de: REQ-VAL-01–03/11, REQ-MAN-01/37, REQ-SAL-05/08–09/13/18 y ACE-17/34/46/64/69.
  - Cambia: prueba integrada y nota breve de cierre; no amplía modelo técnico ni publica paquete/imagen.
  - Se comprueba con: cargo test --test f3_end_to_end --locked y validadores JSON Schema existentes. Demuestra configuración → referencia → individuo → verdad → FASTQ y que verdad no procede de evidencia; puertas comunes.

> **PAUSA-F3.5 / Revisión de F3.** Usar el primer conjunto FASTQ extremo a extremo antes de ampliar biología o técnica. Confirmar utilidad como entrada externa de mapeo/detección, hechos reales en estados/manifiestos y ausencia de escenarios, protocolos, persistencia, empaquetado u hardening F9.

## 10. Matriz de evidencia y riesgos

| Hito | Resultado construido | Evidencia | Riesgo protegido |
|---|---|---|---|
| F3.1 | FASTQ real por flujo desde haplotipos | pruebas focales de admisión, identidad, muestreo y FASTQ; puertas Cargo | dosis/ploidía, N, cobertura, procedencia o memoria incorrectas sin fallo visible |
| F3.5 | Validación FASTQ y evento exportable mínimo | verdad sintética independiente, recorrido extremo a extremo y JSON Schema de evento | aceptar incoherencia, tratar azar como error o declarar completa una salida parcial |

## 11. Riesgos pendientes

- F3 mínimo no equivale a perfil universal validado: F4 es propietario de perfiles cuantitativos, pares y degradaciones.
- Antes del cierre, una verdad manual debe demostrar que cobertura WGS básica y final FASTQ se distinguen; la cobertura alineada es externa al simulador.
- El manifiesto de evento extiende el patrón F1 con hechos reales; no sustituye DTPD-108 ni F9.
