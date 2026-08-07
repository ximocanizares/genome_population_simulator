# Registro de decisiones técnicas del simulador

- **Versión:** 0.1
- **Estado:** borrador para revisión
- **Fecha:** 7 de agosto de 2026
- **Especificación de origen:**
  [especificaciones-0.1.md](./especificaciones-0.1.md)
- **Arquitectura relacionada:**
  [arquitectura-0.1.md](./arquitectura-0.1.md)
- **Hoja de ruta:** [hoja-de-ruta-0.1.md](./hoja-de-ruta-0.1.md)

## 1. Finalidad

Este registro separa las decisiones ya confirmadas de las propuestas de
arquitectura, las elecciones todavía pendientes y el trabajo futuro. Su
finalidad es impedir que una decisión aparezca de forma implícita durante la
implementación.

No es un catálogo de dependencias ni un plan de trabajo. Nombrar una opción no
autoriza su adopción. Una decisión pendiente solo cambia de estado después de
documentar la evaluación y aprobar sus consecuencias.

## 2. Clasificación

- **Confirmada:** procede directamente de una fuente de diseño vigente y no se
  reabre en este documento.
- **Propuesta:** forma parte del borrador de arquitectura y necesita aprobación
  antes de considerarse vinculante.
- **Pendiente:** requiere comparar opciones o aportar evidencia antes de
  planificar o implementar la fase indicada.
- **Futura:** está fuera del alcance inmediato o depende de evidencia de fases
  posteriores. No bloquea la siguiente fase.

Cada entrada identifica su fuente, motivo, consecuencias y fase afectada. En las
pendientes se incluyen opciones, recomendación razonada y fase que bloquea.

## 3. Discrepancia entre fuentes

**DOC-001 — Lenguaje del núcleo.** La hoja de ruta 0.1 declara «Lenguaje del
núcleo: Rust» y F0 depende de la decisión de usar Rust. En cambio, las secciones
4.1, 4.7 y 9 de la especificación 0.1 todavía afirman que el lenguaje no se ha
elegido. Para esta arquitectura prevalece la decisión posterior y explícita de
la hoja de ruta, sin modificar la especificación ni deducir de ella crates,
formatos o herramientas. Conviene corregir la redacción contradictoria en la
próxima revisión de la especificación.

## 4. Decisiones confirmadas

| ID | Decisión y fuente | Motivo confirmado | Consecuencias |
|---|---|---|---|
| DTC-001 | El núcleo se diseña en Rust. Hoja de ruta, cabecera y F0 | Es una decisión explícita de la hoja de ruta | Los contratos aprovecharán tipos, errores explícitos y propiedad de datos; no selecciona crates ni versión de Rust |
| DTC-002 | La configuración y los manifiestos son estructurados y legibles por programa. REQ-CONF-01, REQ-CONF-10 y REQ-MAN-02 | Permitir validación, reproducción y trazabilidad | El formato y el esquema siguen pendientes |
| DTC-003 | Existen tres vías biológicas que convergen en la misma verdad faseada. Secciones 4.3.2–4.3.4 | Evitar modelos incompatibles por vía | Reproducción, haplotipos, variantes e identidades deben reutilizar contratos comunes |
| DTC-004 | Verdad biológica y evidencia técnica están separadas. REQ-VER-32 a REQ-VER-38 | Los artefactos técnicos no pueden alterar el genotipo verdadero | VCF y verdad se construyen sin inferir desde FASTQ/BAM/CRAM |
| DTC-005 | Cada evento técnico es independiente y usa una tecnología, un perfil, un mapeador y una subsemilla. REQ-CONF-101 a REQ-CONF-103 y REQ-TEC-01 | Representar experimentos distintos sin mezclar datos | Los FASTQ no se reutilizan entre eventos y cada evento tiene manifiesto propio |
| DTC-006 | Cada conjunto pertenece a una referencia de mapeo y sus variantes son locales. REQ-DAT-06, REQ-VER-05 y REQ-VER-10 | Evitar correspondencias falsas entre referencias | No habrá identificador global ni traducción registro a registro entre VCF |
| DTC-007 | El modo realista es predeterminado y el directo es opcional. REQ-CONF-67 y REQ-DAT-17 a REQ-DAT-20 | El uso principal necesita FASTQ finales y el directo cubre casos controlados | Ambos modos comparten reglas funcionales, pero no la misma ruta técnica |
| DTC-008 | CRAM es el alineamiento predeterminado, BAM la alternativa; pueden desactivarse. REQ-CONF-64 y REQ-FMT-01 | Decisión funcional explícita | La versión, biblioteca y herramientas de escritura/indexación siguen pendientes |
| DTC-009 | minimap2 es el mapeador predeterminado y sustituible. REQ-CONF-69 y REQ-DAT-21 | Decisión funcional explícita | La versión, invocación, catálogo de alternativas y validación siguen pendientes |
| DTC-010 | La configuración usa coordenadas base 1 e intervalos cerrados; cada formato externo conserva su convenio. REQ-CONF-28 y REQ-IND-38/39 | Evitar posiciones ambiguas | La representación interna y las conversiones concretas siguen pendientes |
| DTC-011 | La reproducibilidad biológica exige misma semilla, configuración, entradas y versión. REQ-CONF-26/27 y REQ-REP-22/23 | Poder repetir decisiones sin prometer compatibilidad indefinida | No confirma identidad bit a bit ni garantías entre versiones o plataformas |
| DTC-012 | La primera versión prioriza casos pequeños y no impone límites funcionales generales. REQ-NOF-01 a REQ-NOF-04 | Permitir validar primero y escalar con evidencia | Los límites prácticos y objetivos cuantitativos deben medirse y acordarse |
| DTC-013 | Los procesos y formatos externos se validan antes de declarar éxito. REQ-FMT-04/07, REQ-VAL-11 y REQ-SAL-05 | Evitar resultados plausibles pero incoherentes | Un artefacto fallido o incompleto nunca se publica como completo |
| DTC-014 | No se ejecuta código personalizado desde la configuración. REQ-CONF-108 a REQ-CONF-110 | Mantener reglas validables y reproducibles | El catálogo declarativo solo admite operaciones conocidas por el programa |

## 5. Decisiones arquitectónicas propuestas

Estas entradas describen el borrador
[arquitectura-0.1.md](./arquitectura-0.1.md). Su estado no implica que exista
implementación.

| ID | Propuesta | Requisitos que satisface | Consecuencias y revisión necesaria |
|---|---|---|---|
| DTP-001 | Arquitectura orientada al dominio con adaptadores para configuración, formatos y procesos externos | REQ-CONF, REQ-REF, REQ-DAT, REQ-VER | Aísla biología de crates y herramientas; aprobar fronteras antes de F1 |
| DTP-002 | Configuración declarada y configuración efectiva son objetos distintos e inmutables en su ámbito | Precedencias REQ-CONF-18, 29–31, 34 y 106 | Permite auditar predeterminados y normalizaciones; debe probarse que ninguna fase muta retroactivamente su entrada |
| DTP-003 | Identidades tipadas por categoría y registro central de unicidad | REQ-IND-25/26, REQ-EVO-53/54, REQ-RG-04 y REQ-MAN-16 | Impide confundir individuo, evento, conjunto o variante; la sintaxis queda pendiente |
| DTP-004 | Coordenadas etiquetadas con referencia, secuencia y convenio, convertidas por un único servicio | REQ-IND-38 a REQ-IND-48 y REQ-REC-14 | Hace explícita toda conversión; falta decidir la representación canónica interna |
| DTP-005 | Verdad genómica basada en segmentos trazables y eventos ordenados, materializable como secuencia | REQ-IND-28 a REQ-IND-52 y REQ-REC-11 a REQ-REC-29 | Permite fase y linaje sin exigir copias completas; debe medirse antes de fijar la estructura física |
| DTP-006 | Flujos aleatorios con ámbitos nombrados para biología, muestreo y eventos técnicos | REQ-EVO-12/62 y REQ-CONF-27 | Evita que una observación cambie la evolución; el algoritmo y la derivación de subsemillas siguen pendientes |
| DTP-007 | Resultados provisionales que solo pasan a completos después de sus validaciones | REQ-SAL-05, REQ-VAL-11 y REQ-MAN-73 | Evita publicar parciales; quedan pendientes la atomicidad física, reanudación y limpieza |
| DTP-008 | Reproducción, recombinación y selección declarativa compartidas por evolución y familias | REQ-CRU-06/17/31/32 y REQ-EVO-27/29 | Evita rutas biológicas divergentes; exige contratos comunes antes de F7/F8 |
| DTP-009 | Recursos comunes de referencia y ámbitos separados por evento y conjunto | REQ-MAN-72 a REQ-MAN-75 y REQ-SAL-17 a REQ-SAL-19 | Evita duplicación y mezcla; nombres y disposición física siguen pendientes |
| DTP-010 | Validación como consumidor de solo lectura, sin corregir verdad ni evidencia | REQ-VAL-01 a REQ-VAL-13 | Las desviaciones se informan, no se ocultan mediante reparaciones posteriores |

## 6. Decisiones pendientes que bloquean F1

### DTPD-001. Formato y esquema de configuración

- **Opciones:** un formato legible orientado a configuración; un formato de
  intercambio general; o separar formato de autoría y representación canónica.
  Entre las familias que deberán evaluarse están TOML, YAML y JSON, sin que su
  mención constituya una selección.
- **Recomendación:** comparar expresividad, mensajes de error, validación de
  esquemas, soporte de comentarios, estabilidad, seguridad y soporte Rust con
  configuraciones mínimas de F1 y estructuras futuras de F7/F8. Considerar una
  representación canónica efectiva distinta del fichero escrito por la persona.
- **Consecuencias:** determina inclusión de ficheros, sintaxis de referencias,
  versionado, diagnósticos y facilidad de migración. Elegir solo para F1 sin
  comprobar listas, ámbitos y grafos futuros podría obligar a romper el formato.
- **Bloquea:** F1.2 y F1.4; también condiciona F2–F8.

### DTPD-002. Esquema y formato de manifiestos

- **Opciones:** reutilizar el formato de configuración; usar un formato de
  intercambio distinto; o separar manifiesto general, eventos y tablas grandes
  en artefactos enlazados.
- **Recomendación:** evaluar primero volumen, lectura parcial, versionado y
  capacidad de enlazar artefactos. No exigir que configuración y manifiesto
  compartan formato si sus necesidades difieren.
- **Consecuencias:** afecta trazabilidad, recuperación, comparación entre
  ejecuciones y compatibilidad futura. Un manifiesto monolítico puede crecer en
  exceso; uno fragmentado exige integridad referencial explícita.
- **Bloquea:** F1.4.

### DTPD-003. Versionado y compatibilidad de esquemas

- **Opciones:** compatibilidad solo con la versión exacta; lectura de varias
  versiones con migración; o esquema canónico con conversores explícitos.
- **Recomendación:** exigir desde F1 un identificador de versión y rechazo claro
  de versiones desconocidas; decidir migraciones solo cuando exista una segunda
  versión real.
- **Consecuencias:** evita interpretar silenciosamente campos con otro
  significado y reduce abstracciones prematuras.
- **Bloquea:** F1.2 y F1.4.

### DTPD-004. Generador aleatorio y estrategia de subsemillas

- **Opciones:** flujo secuencial único; generador divisible jerárquicamente; o
  generador basado en contador y claves de ámbito. También debe elegirse el
  algoritmo y su crate tras evaluar estabilidad y distribución.
- **Recomendación:** descartar el flujo único porque incumple la independencia
  de muestreos y eventos. Probar ámbitos nombrados y derivación determinista sin
  fijar todavía si se implementan mediante división o contador.
- **Consecuencias:** determina reproducibilidad, paralelismo futuro, orden de
  iteración y compatibilidad entre versiones. Cambiarlo altera decisiones
  aleatorias aunque la semilla sea igual.
- **Bloquea:** F1.3. Su extensión por procesos bloquea F6–F9.

### DTPD-005. Garantía exacta de reproducibilidad

- **Opciones:** solo misma versión y plataforma; misma versión entre plataformas
  admitidas; identidad lógica pero no binaria; o estabilidad declarada de
  secuencias aleatorias entre versiones concretas.
- **Recomendación:** para F1 prometer únicamente la garantía funcional ya
  confirmada de misma versión, entradas y configuración. Diseñar pruebas que
  separen decisiones biológicas de metadatos y posponer garantías
  multiplataforma hasta contar con una matriz real.
- **Consecuencias:** afecta tipos numéricos, orden de colecciones, concurrencia,
  herramientas externas y expectativas de usuarios. No debe prometerse
  identidad bit a bit sin evidencia.
- **Bloquea:** la redacción del criterio de cierre F1.3; garantías mayores no
  bloquean F1.

### DTPD-006. Identificadores

- **Opciones:** identificadores jerárquicos legibles; identificadores aleatorios
  globales; contadores locales con ámbito; o combinaciones de identificador
  estable y etiqueta humana.
- **Recomendación:** comparar estabilidad, detección de colisiones, longitud en
  cabeceras genómicas y capacidad de reproducir identificadores lógicos. Separar
  siempre identidad estable de nombre presentado.
- **Consecuencias:** afecta todas las relaciones, `SM`, `RG`, rutas, manifiestos
  y estabilidad de VCF entre eventos. La sintaxis no debe codificar estado
  mutable.
- **Bloquea:** F1.4.

### DTPD-007. Identidad de entradas y sumas de comprobación

- **Opciones:** algoritmos criptográficos habituales como SHA-256 o alternativas
  aprobadas; suma de contenido completo o estrategia adicional para ficheros muy
  grandes.
- **Recomendación:** elegir un algoritmo ampliamente interoperable y registrar
  siempre nombre y valor. Medir el coste de lectura, sin sustituir la suma de
  exportación por metadatos débiles como tamaño y fecha.
- **Consecuencias:** permite comprobar FASTA exacto y paquetes; el cálculo puede
  añadir E/S significativa.
- **Bloquea:** la identidad verificable de F1.5. El cierre completo de
  exportación corresponde a F9.6.

### DTPD-008. Lectura, escritura y validación de FASTA

- **Opciones:** crate Rust; biblioteca nativa enlazada; herramienta externa; o
  combinación con validación independiente.
- **Recomendación:** comparar aceptación estricta de IUPAC, conservación de
  identificadores, indexación futura, consumo de memoria, licencias y
  distribución. El contrato de referencia debe quedar separado del adaptador.
- **Consecuencias:** afecta copia, validación, materialización e índices. Una
  herramienta externa añade gestión de procesos y versiones.
- **Bloquea:** F1.5.

### DTPD-009. Plataforma y distribución iniciales

- **Opciones:** una plataforma de referencia; varias plataformas nativas; o
  distribución encapsulada además del binario.
- **Recomendación:** declarar primero una plataforma verificable coherente con
  minimap2 y las herramientas genómicas previstas, y ampliar solo con una matriz
  automatizada. Evaluar arquitectura de CPU, sistema de ficheros y disponibilidad
  de procesos externos.
- **Consecuencias:** condiciona compilación, rutas, ejecución de herramientas,
  reproducibilidad y soporte. Este registro no selecciona sistema operativo ni
  arquitectura.
- **Bloquea:** F1.1 y la definición exacta de reproducibilidad de F1.3.

### DTPD-010. Política de crates y dependencias nativas

- **Opciones:** crates Rust puros; bibliotecas nativas con enlaces; procesos
  externos; o combinaciones por formato.
- **Recomendación:** decidir por requisito y no por uniformidad. Evaluar
  corrección del estándar, mantenimiento, seguridad, licencia, rendimiento,
  estabilidad de API y coste de distribución para cada dependencia.
- **Consecuencias:** una biblioteca nativa puede mejorar interoperabilidad pero
  complica compilación; un proceso externo aísla fallos pero exige instalación y
  control de versiones.
- **Bloquea:** F1.1 y las decisiones de adaptadores que use F1.

### DTPD-011. Modelo de errores y diagnósticos

- **Opciones:** categorías jerárquicas tipadas; códigos planos; o errores
  internos enriquecidos y una capa estable de diagnóstico de usuario.
- **Recomendación:** separar causa técnica de categoría y alcance funcional,
  conservar cadenas causales y definir desde F1 códigos estables para entradas,
  configuración, recursos, procesos y defectos internos.
- **Consecuencias:** permite pruebas y mensajes comprensibles sin exponer
  detalles accidentales de crates.
- **Bloquea:** F1.2 y el criterio de rechazo limpio de F1.

### DTPD-012. Publicación atómica, reanudación y residuos

- **Opciones:** ejecución no reanudable con área provisional; reanudación por
  artefacto validado; o registro transaccional más completo.
- **Recomendación:** en F1 asegurar que lo incompleto no figure como completo y
  conservar diagnóstico. Posponer la reanudación automática hasta medir el
  coste de repetir fases reales.
- **Consecuencias:** afecta nombres temporales, manifiestos, recuperación tras
  interrupción y uso de almacenamiento.
- **Bloquea:** F1.4 para la semántica de completitud; la reanudación no bloquea
  F1.

### DTPD-013. Escenarios y métricas iniciales de rendimiento

- **Opciones:** casos sintéticos por tamaño de referencia; matriz de referencia,
  individuos y salidas; o recorridos específicos por componente.
- **Recomendación:** definir al menos un caso pequeño de aceptación y uno de
  observación de crecimiento para F1, midiendo tiempo máximo observado, memoria
  máxima y almacenamiento. No fijar umbrales de producto sin medidas.
- **Consecuencias:** aporta la evidencia de F0.4 y una línea base, pero no
  convierte el resultado local en una garantía universal.
- **Bloquea:** cierre de F0.4 y la aceptación del plan F1. No hay objetivo
  cuantitativo confirmado todavía.

### DTPD-014. Herramientas de prueba y validación independiente

- **Opciones:** validadores externos oficiales o de referencia; bibliotecas
  alternativas; comprobadores propios limitados; o combinación.
- **Recomendación:** usar una fuente independiente para cada invariante crítica
  cuando exista y evitar validar un formato solo con el mismo escritor. Registrar
  herramienta y versión en las pruebas.
- **Consecuencias:** añade dependencias de desarrollo y posiblemente de sistema,
  pero protege contra errores plausibles del propio adaptador.
- **Bloquea:** F1.1 y evidencia de F1.5.

### DTPD-015. Versión y política del conjunto de herramientas de Rust

- **Opciones:** fijar una versión exacta del compilador; declarar una versión
  mínima admitida y probar también una estable actual; o seguir únicamente el
  canal estable sin garantía de versión mínima.
- **Recomendación:** elegir después de comprobar compatibilidad de las
  dependencias aprobadas y registrar de forma reproducible compilador, edición y
  componentes de formato y análisis. No declarar una versión mínima que no se
  pruebe continuamente.
- **Consecuencias:** afecta construcción, características disponibles,
  distribución, reproducibilidad y política de actualización. Este borrador no
  elige versión, edición ni canal.
- **Bloquea:** F1.1.

## 7. Decisiones pendientes para fases posteriores

| ID | Decisión | Opciones | Recomendación razonada | Consecuencias | Bloquea |
|---|---|---|---|---|---|
| DTPD-101 | Representación interna de coordenadas | Base 1 cerrada; base 0 semiabierta; valor etiquetado con varias representaciones | Evaluar operaciones de segmentos y evitar conversiones implícitas; mantener tipos de frontera aunque se elija una canónica | Afecta errores de borde, normalización y compatibilidad con formatos | F2.1/F2.2 |
| DTPD-102 | Estado haplotípico y fichero obligatorio de verdad | Secuencias completas; lista de variantes; segmentos persistentes; modelo híbrido | Probar reconstrucción, memoria y composición con SNP, indels y SV antes de elegir; separar modelo interno y formato de intercambio | Es el contrato central de fase, dosis y linaje | F2.1/F2.4 |
| DTPD-103 | VCF, normalización y complemento | Versiones VCF candidatas; normalizador propio o externo; complemento tabular, documental o binario | Seleccionar por capacidad de poliploidía, BND, repeticiones, BGZF/índices e interoperabilidad; demostrar reconstrucción sin pérdida | Afecta campos, cabecera, identificadores locales y validadores | F2.5 y F5.5 |
| DTPD-104 | Adaptadores SAM/BAM/CRAM e índices | Crate; biblioteca nativa; herramientas externas; combinación | Comparar cumplimiento de estándar, CRAM con referencia exacta, ordenación, índices, cabeceras y distribución | Determina E/S, procesos, errores y rendimiento | F3.3 |
| DTPD-105 | Ejecución de minimap2 y otros mapeadores | Proceso hijo local; servicio encapsulado; capa de herramientas configurable | Empezar con un contrato de proceso local verificable sin ocultar comando, versión o código de salida; no aprobar otros mapeadores sin compatibilidad | Afecta instalación, seguridad, logs y reproducibilidad | F3.2 |
| DTPD-106 | Perfiles WGS cuantitativos | Perfiles incluidos versionados; artefactos externos; perfiles calibrados | Validar cada perfil con datos adecuados y distinguir claramente predeterminado, sustituciones y empírico | Fija distribuciones de longitud, error, calidad y artefactos | F4.1–F4.3 |
| DTPD-107 | Calibración empírica | Implementación Rust; herramienta auxiliar aprobada; flujo híbrido | Elegir por capacidad para separar error, variante y mapeo, y registrar método/versión; Python solo si se aprueba como auxiliar | Afecta distribución, fiabilidad y esquema del perfil reutilizable | F4.7 |
| DTPD-108 | Verdad y procedencia de lecturas | Registro por lectura; registros por bloques; combinación con agregados | Mantener detalle solo donde lo exige la especificación y agregados para descartes; medir volumen antes de fijar formato | Puede dominar el almacenamiento y la validación | F4.4/F4.5 |
| DTPD-109 | Fórmulas de F y Ne reproductivo | Fórmulas publicadas compatibles con ploidía; variantes por diseño reproductivo | Evaluar contra el significado funcional confirmado y casos simulados; documentar supuestos sin introducir modelos polisómicos/disómicos no aprobados | Cambia estadísticas, avisos y pruebas, no el genotipo | F6.3 y F7.3 |
| DTPD-110 | Distribuciones y tolerancias poblacionales | Catálogos de distribuciones; precisión flotante; racionales o enteros escalados en ámbitos concretos | Elegir por semántica y estabilidad numérica, con casos dentro/fuera de tolerancia y sin sustituir `0,5/C` | Afecta normalización, frecuencias y reproducibilidad | F6.1/F6.2 |
| DTPD-111 | Mapas y posiciones de recombinación | Formatos tabulares candidatos; generadores de Poisson; interpolación de mapas | Evaluar cobertura completa, monotonicidad, conversiones y prueba de intervalos colineales antes de seleccionar | Afecta bivalentes, linaje y casos estructurales | F7.1 |
| DTPD-112 | Representación de grafos y reglas declarativas | Estructuras anidadas; tablas de nodos/aristas; lenguaje declarativo limitado sin código | Elegir la sintaxis más validable que represente Y/O/NO, ámbitos y dependencias sin ejecutar código | Afecta mensajes de ciclo, migración, fundación y familias | F7.5 y F8.1 |
| DTPD-113 | Paralelismo reproducible | Por referencia, individuo, evento, región o tareas; ejecución secuencial de referencia | Posponer la elección hasta perfilar y conservar ámbitos aleatorios y orden lógico independientes del planificador | Puede cambiar uso de memoria y orden de salida; no debe cambiar verdad | F9.3 |
| DTPD-114 | Formato del paquete exportado | Directorio autocontenido; archivo contenedor; formatos archivados comprimidos | Elegir por interoperabilidad, preservación de permisos/rutas e integridad antes/después, con eventos grandes | Afecta portabilidad y almacenamiento duplicado | F9.6 |
| DTPD-115 | Objetivos y presupuestos cuantitativos | Umbrales por escenario; presupuesto configurable; perfiles de recursos | Fijarlos solo después de medir escenarios representativos y separar advertencia de rechazo por presupuesto | Define expectativas de escala y criterios de publicación | F9.1/F9.4 |

## 8. Decisiones futuras que no bloquean v1

| ID | Tema futuro | Motivo para aplazarlo | Contrato actual que debe preservarse |
|---|---|---|---|
| DTF-001 | Oxford Nanopore y PacBio CLR | Perfiles y errores más complejos, fuera de alcance | Tecnología y perfil identificados y versionados por evento |
| DTF-002 | Paneles, amplicones, GBS y otros protocolos | v1 se limita a WGS de DNA genómico | Objetivo y preparación pertenecen a la biblioteca |
| DTF-003 | Ploidía individual, cromosomas sexuales y aneuploidía | v1 exige la misma ploidía por secuencia en todos los individuos | Ploidía explícita por secuencia, sin codificarla en algoritmos diploides |
| DTF-004 | Ploidías impares, multivalentes y doble reducción | Reproducción sexual v1 usa ploidía par y bivalentes | Contrato de reproducción separado del estado genotípico general 1–10 |
| DTF-005 | Conversión génica e interferencia | v1 solo modela entrecruzamientos sin interferencia | Eventos de transmisión y segmentos con identidad y procedencia |
| DTF-006 | Selección poligénica, ambiental, estabilizadora o disruptiva | v1 usa un locus, valor genético igual al carácter y selección direccional | Servicio de selección separado del estado genotípico |
| DTF-007 | Transposición espontánea y estados activos de TE | v1 solo inserta eventos solicitados | Catálogo y eventos TE identificados sin motor autónomo implícito |
| DTF-008 | Relación de variantes entre referencias | v1 prohíbe correspondencia registro a registro | Linaje de referencias separado de identificadores locales de variante |
| DTF-009 | Compatibilidad reproducible entre versiones | La especificación no la garantiza | Manifiesto registra versión, entradas, perfiles y dependencias |
| DTF-010 | Reanudación distribuida y ejecución remota | No hay requisito ni plataforma aprobada | Estados de resultado distinguibles y contratos de E/S explícitos |

## 9. Reglas para resolver una decisión

Una decisión pendiente pasa a confirmada o se sustituye por otra entrada solo
cuando consten:

1. requisito y fase que la necesitan;
2. opciones realmente evaluadas y criterios aplicados;
3. evidencia reproducible, incluidos prototipos o medidas cuando correspondan;
4. compatibilidad con Rust, formatos externos y plataformas admitidas;
5. costes de distribución, mantenimiento, seguridad y licencias para cada
   dependencia;
6. consecuencias sobre reproducibilidad, datos existentes y trabajo futuro;
7. pruebas que demostrarán la decisión;
8. aprobación explícita y fecha.

No se añadirá un crate porque sea habitual, ni se fijará un formato porque sea
legible, ni se prometerá rendimiento o reproducibilidad porque un caso aislado
lo haya conseguido.

## 10. Resumen de bloqueos

F1 no debe planificarse como ejecutable hasta resolver DTPD-001 a DTPD-015 en la
medida indicada por cada entrada. Algunas decisiones contienen una parte mínima
bloqueante y otra ampliación aplazable; por ejemplo, F1 necesita semántica de
completitud, pero no necesita reanudación automática.

Las decisiones DTPD-101 a DTPD-115 deben resolverse antes de planificar el hito
que señalan. No bloquean la redacción del plan de F1 salvo que su evaluación
revele una incompatibilidad con un contrato base de configuración, identidad,
referencia o reproducibilidad.
