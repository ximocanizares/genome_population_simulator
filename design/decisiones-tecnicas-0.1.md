# Registro de decisiones técnicas del simulador

- **Versión:** 0.2
- **Estado:** F1 cerrado; DTC-089 cierra las decisiones técnicas F5;
  DTC-094 a DTC-097 cierran el deslizamiento biológico de los microsatélites;
  DTC-098 conecta la inserción de tractos en una referencia suministrada;
  DTC-099 hace configurable el multiplicador de la PCR en los homopolímeros
- **Fecha:** 1 de septiembre de 2026
- **Especificación de origen:**
  [especificaciones-0.1.md](./especificaciones-0.1.md)
- **Arquitectura relacionada:**
  [arquitectura-0.1.md](./arquitectura-0.1.md)
- **Hoja de ruta:** [hoja-de-ruta-0.1.md](./hoja-de-ruta-0.1.md)

Este registro es derivado: ninguna paráfrasis técnica amplía los requisitos
normativos. Ante una divergencia prevalece el REQ y debe corregirse este
registro.
El detalle recuperable y el destino de cada DTF se mantienen en
[registro-ampliaciones-futuras-0.1.md](./registro-ampliaciones-futuras-0.1.md).

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

**DOC-001 — Lenguaje del núcleo (resuelto).** La especificación, la arquitectura
y la hoja de ruta confirman Rust para el núcleo. La versión del lenguaje, los
*crates*, formatos y herramientas continúan separados como decisiones de
planificación; ya no existe una discrepancia entre fuentes.

## 4. Decisiones confirmadas

| ID | Decisión y fuente | Motivo confirmado | Consecuencias |
|---|---|---|---|
| DTC-001 | El núcleo se diseña en Rust. Hoja de ruta, cabecera y F0 | Es una decisión explícita de la hoja de ruta | Los contratos aprovecharán tipos, errores explícitos y propiedad de datos; no selecciona crates ni versión de Rust |
| DTC-002 | La configuración y los manifiestos son estructurados y legibles por programa. REQ-CONF-01, REQ-CONF-10 y REQ-MAN-02 | Permitir validación, reproducción y trazabilidad | DTC-034 concreta la configuración de autoría y DTC-035 los manifiestos y la configuración efectiva serializada |
| DTC-003 | Existen tres vías biológicas que convergen en la misma verdad faseada. Secciones 4.3.2–4.3.4 | Evitar modelos incompatibles por vía | Reproducción, haplotipos, variantes e identidades deben reutilizar contratos comunes |
| DTC-004 | Verdad biológica y evidencia técnica están separadas. REQ-VER-32 a REQ-VER-38 | Los artefactos técnicos no pueden alterar el genotipo verdadero | VCF y verdad se construyen sin inferir desde FASTQ/BAM/CRAM |
| DTC-005 | Cada evento técnico es independiente y usa una tecnología, un perfil, un mapeador y una subsemilla. REQ-CONF-101 a REQ-CONF-103 y REQ-TEC-01 | Representar experimentos distintos sin mezclar datos | Los FASTQ no se reutilizan entre eventos y cada evento tiene manifiesto propio |
| DTC-006 | Cada conjunto pertenece a una referencia de mapeo y sus variantes son locales. REQ-DAT-06, REQ-VER-05 y REQ-VER-10 | Evitar correspondencias falsas entre referencias | No habrá identificador global ni traducción registro a registro entre VCF |
| DTC-007 | El modo realista es predeterminado y el directo es opcional. REQ-CONF-67 y REQ-DAT-17 a REQ-DAT-20 | El uso principal necesita FASTQ finales y el directo cubre casos controlados | Ambos modos comparten reglas funcionales, pero no la misma ruta técnica |
| DTC-008 | CRAM es el alineamiento predeterminado, BAM la alternativa; pueden desactivarse. REQ-CONF-64 y REQ-FMT-01 | Decisión funcional explícita | La versión, biblioteca y herramientas de escritura/indexación siguen pendientes |
| DTC-009 | minimap2 es el mapeador predeterminado y sustituible. REQ-CONF-69 y REQ-DAT-21 | Decisión funcional explícita | La versión, invocación, catálogo de alternativas y validación siguen pendientes |
| DTC-010 | La configuración usa coordenadas base 1 e intervalos cerrados; cada formato externo conserva su convenio. REQ-CONF-28 y REQ-IND-38/39 | Evitar posiciones ambiguas | DTC-056 resuelve la representación interna base 0 semiabierta, los tipos distintos de base/punto/intervalo y las conversiones de frontera |
| DTC-011 | La reproducibilidad biológica exige misma semilla, configuración, entradas y versión. REQ-CONF-26/27 y REQ-REP-22/23 | Poder repetir decisiones sin prometer compatibilidad indefinida | DTC-038 concreta los cuatro niveles y no confirma identidad bit a bit ni garantías entre versiones o plataformas no verificadas |
| DTC-012 | La primera versión prioriza casos pequeños y no impone límites funcionales generales. REQ-NOF-01 a REQ-NOF-04 | Permitir validar primero y escalar con evidencia | Los límites prácticos y objetivos cuantitativos deben medirse y acordarse |
| DTC-013 | Los procesos y formatos externos se validan antes de declarar éxito. REQ-FMT-04/07, REQ-VAL-11 y REQ-SAL-05 | Evitar resultados plausibles pero incoherentes | Un artefacto fallido o incompleto nunca se publica como completo |
| DTC-014 | No se ejecuta código personalizado desde la configuración. REQ-CONF-108 a REQ-CONF-110 | Mantener reglas validables y reproducibles | El catálogo declarativo solo admite operaciones conocidas por el programa |
| DTC-015 | En v1, una unidad reproductiva solo contiene individuos basados en la misma referencia estructural. Pregunta 4 del informe adversarial (RAD-012); REQ-IND-53, REQ-MIG-09 y REQ-CRU-44 | No se ha definido la identidad de loci, los mapas ni la referencia de la descendencia entre estructuras distintas | Subpoblaciones y ramas aisladas pueden usar referencias distintas; el mapeo contra varias referencias sigue siendo independiente de la reproducción |
| DTC-016 | En v1, cada evento que genera FASTQ, BAM o CRAM es su único experimento de secuenciación, contiene una tanda y exactamente una biblioteca y un grupo de lectura por individuo. Pregunta 7 del informe adversarial (RAD-027); REQ-CONF-102, REQ-TEC-01 a REQ-TEC-04 y REQ-RG-01 a REQ-RG-07 | Mantener una unidad técnica simple, explícita y verificable por evento | Otra biblioteca del individuo, tanda, perfil, tecnología o control exige otro evento; los eventos solo de verdad o VCF carecen de entidades técnicas |
| DTC-017 | En v1, una sustitución por identificador individual solo puede dirigirse a un individuo directo, fundador o parental expreso existente durante la validación. Pregunta 8 del informe adversarial (RAD-026); REQ-CONF-17 y REQ-CONF-18 | Evitar configurar entidades que todavía no existen y aplicaciones retroactivas | Los descendientes usan ámbitos poblacional, temporal y espacial; un identificador inexistente es error |
| DTC-018 | En F4, una única cobertura objetivo usa la longitud haploide de la referencia estructural de cada individuo y se distingue de las coberturas WGS básica y final del FASTQ. Pregunta 10 del informe adversarial (RAD-004 y RAD-017); REQ-WGS-04, REQ-WGS-11 y REQ-DAT-24 | Un mismo FASTQ no puede tener objetivos distintos por referencia de mapeo | Cada secuencia se suma una vez sin ploidía; haplotipos no cambian el denominador. El evento no admite objetivos por número ni sustituciones individuales, por secuencia o región; la variación realizada se informa sin forzar una profundidad exacta |
| DTC-019 | V1 no modela elementos transponibles como clase funcional. Pregunta 12 del informe adversarial (RAD-006); sección 3 y REQ-TE-01 a REQ-TE-14 retirados | Reducir el alcance y eliminar la contradicción entre referencias derivadas y eventos TE | No hay catálogos, familias, mecanismos ni metadatos TE; las copias e inserciones se clasifican mediante operaciones genéricas |
| DTC-020 | La generación por eventos aplica primero los explícitos y después los aleatorios en el orden duplicaciones, deleciones, inversiones, translocaciones, repeticiones en tándem, *indels* y SNP. Pregunta 13 del informe adversarial (RAD-008); REQ-IND-54 a REQ-IND-57 | Hacer único y comprobable el resultado ante colisiones y agotamiento del espacio válido | Cada evento se elige del espacio válido actualizado; no hay reintentos sobre inválidos ni límites de intentos, y la insuficiencia hace fallar la generación conservando los recuentos sorteado, aceptado y activo. El muestreo inicial desde una matriz SNV de DTC-052 no es otro evento de esta secuencia |
| DTC-021 | Cada conjunto solicitado produce exactamente un VCF multimuesta estándar; un caso de verdad usa uno o varios registros y, solo si pierde detalle, una única entrada complementaria. Excepción desde el 3 de octubre de 2026 (REQ-VER-25, DTC-065): un caso de efecto neto nulo no aporta ALT propio y puede no usar ningún registro. Pregunta 14 del informe adversarial (RAD-011 y RAD-023); REQ-VER-04 y REQ-VER-19 a REQ-VER-26 | Mantener interoperabilidad VCF sin perder la verdad de casos complejos | El VCF sigue siendo válido por sí mismo; la reconstrucción completa usa el par enlazado, BND conserva su cardinalidad estándar y referencias distintas permanecen independientes |
| DTC-022 | V1 usa exclusivamente Poisson para recuentos de variantes aleatorias obtenidos desde una media y el modelo escalonado simple de una unidad para mutaciones de repeticiones en tándem. Pregunta 17 del informe adversarial (RAD-020); REQ-CONF-25, REQ-CONF-61/62 y REQ-IND-18/35 | Cerrar el catálogo mínimo y hacer observables las configuraciones válidas | Eventos explícitos, número exacto de loci de una matriz y modo controlado no son distribuciones alternativas; el espectro alélico y el muestreo de genotipos tampoco son recuentos Poisson. La tasa y el sesgo expansión/contracción siguen configurables y cualquier otra distribución, modelo o salto se rechaza antes de simular |
| DTC-023 | V1 aplaza completamente la calibración empírica automática desde datos reales. Pregunta 18 del informe adversarial (RAD-024); REQ-CONF-85 a REQ-CONF-87 y REQ-VAL-49 | Evitar una capacidad amplia sin contrato verificable y mantener cerrado el alcance de perfiles | FASTQ y referencia de calibración, solicitudes de estimación y perfiles empíricos reutilizables se rechazan antes de simular; se conservan los perfiles incorporados de Illumina y sus sustituciones explícitas. PacBio HiFi queda en DTF-049 |
| DTC-024 | V1 separa validaciones deterministas de informes estadísticos y no usa tolerancias de aceptación para magnitudes aleatorias. Pregunta 19 del informe adversarial (RAD-025); REQ-VAL-01 y REQ-VAL-10 a REQ-VAL-13 | Evitar que una realización aleatoria válida cambie de éxito a fallo según un umbral no definido | Las incoherencias exactas detienen el ámbito; cobertura, soporte y tasas realizadas se informan sin opción de elevarlas a error; una medida no calculable conserva estado y causa |
| DTC-025 | N objetivo es el número de intentos reproductivos y v1 no admite Ne objetivo ni límites de descendencia parental. Pregunta 20 del informe adversarial (RAD-010 y REQ-REP-07); REQ-CONF-40 y REQ-REP-08 a REQ-REP-16 | Dar efecto observable a N sin introducir un optimizador de Ne no definido | Cada éxito crea un descendiente, los fallos aleatorios conservan causa y aviso, los parentales son reutilizables; Ne solo se calcula como resultado reproductivo evaluable |
| DTC-026 | Las tasas evolutivas se aplican por oportunidad válida y generación en todas las copias haplotípicas actuales. Pregunta 21 del informe adversarial (RAD-014); REQ-MUT-01 y REQ-MUT-11 a REQ-MUT-15 | Hacer reproducible el denominador después de ploidía, recombinación y cambios estructurales | Cada clase cuenta su espacio inicial, obtiene media tasa×oportunidades y sortea una vez; destinos condicionados no multiplican el denominador y `N` se rige por DTC-032 |
| DTC-027 | La selección de viabilidad de v1 usa un único SNP bialélico protegido y valores completos por dosis. Pregunta 22 del informe adversarial (RAD-015); REQ-SEL-01 a REQ-SEL-17 | Evitar identidades ambiguas y normalización indefinida ante máximo cero | Eventos incompatibles no pueden afectar al locus; identidad/dosis/aptitud inválida o ausencia de aptitud positiva detiene antes de supervivencia y no equivale a extinción |
| DTC-028 | Hardy-Weinberg y el modelo automático F son opciones distintas: HWE usa la distribución binomial o multinomial en ploidías 1–10 y F solo se admite en diploides con las fórmulas confirmadas. Pregunta 23 del informe adversarial (RAD-021); REQ-EVO-33/34 y REQ-EVO-46 a REQ-EVO-52 | Evitar inferir una distribución de dosis poliploide que los extremos de F no determinan | F en P=1 o P>2 se rechaza antes de simular; una distribución explícita prevalece localmente; F realizado solo se calcula en diploides y figura como «no aplicable» en las demás ploidías |
| DTC-029 | La reproducción de v1 exige ploidía par e igual y el mismo estado estructural canónico en todas las copias y parentales de cada secuencia. Pregunta 24 del informe adversarial (RAD-022); REQ-CRU-23 a REQ-CRU-31 y REQ-REC-03 a REQ-REC-08 | Sustituir la compatibilidad parcial no definida por un invariante observable | SNP, *indels* pequeños y repeticiones pueden segregar; toda diferencia estructural invalida al individuo o pareja; una SV fijada idéntica se admite aunque conserve otra historia y una SV de novo hace incompatible al portador para usos posteriores |
| DTC-030 | V1 aplica un único orden generacional: demografía/N, migración, viabilidad, reproducción, segregación/recombinación, mutación de novo, registro de N/Ne/frecuencias y muestreo/salidas. Pregunta 25; REQ-CONF-44/45 y REQ-REP-24/25 | Evitar trayectorias distintas para la misma configuración | Una etapa opcional se omite en su posición; otro orden se rechaza antes de simular; varios eventos demográficos conservan su orden interno explícito y el manifiesto registra lo aplicado y omitido |
| DTC-031 | V1 aplica un único contrato obligatorio de trazabilidad, sin niveles configurables ni reducción automática por escala. Pregunta 26 del informe adversarial (RAD-030); REQ-CRU-15, REQ-REC-11 a REQ-REC-13, REQ-DAT-40, REQ-SAL-02 y REQ-SAL-16 | Distinguir el detalle reconstruible exigido de los resúmenes suficientes sin permitir resultados degradados silenciosamente | Evolución materializa estados completos inicial, final y solicitados y resume el resto; familias materializa todos los individuos; cada lectura final conserva detalle y los descartes solo agregados; un presupuesto insuficiente rechaza antes de simular |
| DTC-032 | `N` representa una base desconocida que se conserva exactamente y nunca se resuelve. Pregunta 27 del informe adversarial (RAD-007); REQ-REF-14/23, REQ-IND-58/59, REQ-DAT-72/73 y REQ-VER-47 | Evitar que operaciones distintas inventen verdades incompatibles desde una base no determinada | SNP, selección y motivos repetidos excluyen `N`; SV puede transportarla; WGS la lee con Phred 0, GC la excluye y cobertura la cuenta; VCF no emite SNP con `REF=N` y usa complemento cuando haga falta |
| DTC-033 | Las sustituciones se resuelven por coincidencia de todos sus selectores y dominancia parcial entre ejes, sin prioridad implícita entre ejes distintos. Pregunta 3 del informe adversarial (RAD-009); REQ-CONF-116 a REQ-CONF-123 | Obtener un único valor efectivo auditable sin inventar una prioridad global que cambie según la familia | Una única regla máxima gana; máximas incomparables, duplicados exactos, regiones solapadas y selectores no autorizados se rechazan antes de simular; los eventos técnicos nunca sustituyen parámetros biológicos |
| <a id="dtc-034"></a>DTC-034 | TOML es el único formato de configuración escrito por la persona usuaria en v1, tanto para el fichero general como para los escenarios. Resolución de DTPD-001; REQ-CONF-01 a REQ-CONF-14, REQ-CONF-51/52, REQ-CONF-101/105 y REQ-CONF-116 a REQ-CONF-123 | Mantener un único formato de autoría legible y hacer explícitas identidad, orden y procedencia sin depender de la jerarquía física del fichero | Todo fichero exige versión de esquema; se conserva una copia exacta del TOML; `DeclaredConfig` representa lo escrito con su ubicación y `EffectiveConfig`, la configuración validada, heredada y resuelta e inmutable en su ámbito. Se rechazan campos desconocidos, claves o campos duplicados y combinaciones inválidas con diagnósticos localizables. Fases, eventos y cruces son listas de tablas con identificadores explícitos; su identidad y precedencia no se deducen de la jerarquía física. DTC-036 fija la compatibilidad; DTC-035, la serialización efectiva. No selecciona crate ni versión de TOML o Rust |
| <a id="dtc-035"></a>DTC-035 | Los manifiestos de v1 son JSON conforme a RFC 8259: existe un `manifest.json` general en la raíz y otro por evento, validados con esquemas JSON Schema Draft 2020-12 distintos que reutilizan definiciones comunes. La configuración efectiva se escribe aparte en `effective-config.json`, versionada y enlazada por el manifiesto general. Resolución de DTPD-002; REQ-MAN-01/02, REQ-MAN-72 a REQ-MAN-88 y REQ-SAL-03/04/17/18 | Separar el documento de autoría TOML de salidas de intercambio verificables y mantener completo el contrato de trazabilidad sin convertir artefactos independientes en contenido duplicado del manifiesto | Cada manifiesto exige `schema_version` y `manifest_kind`, se escribe de forma legible y contiene directamente todo lo que REQ-MAN obliga a registrar o enumerar. Los artefactos independientes exigidos se enlazan mediante identificador, ruta relativa, tipo, tamaño y suma de comprobación. Las listas con semántica de conjunto se ordenan establemente por identificador, sin prometer identidad binaria general. V1 no fragmenta manifiestos en JSONL ni tablas auxiliares; el preflight incluye su tamaño y nunca reduce contenido obligatorio. La configuración original permanece como copia TOML exacta. DTC-036 fija la compatibilidad; DTC-040, el algoritmo de suma; DTC-045, la publicación; DTC-043, la política de dependencias; y DTC-048, toolchain y lock |
| <a id="dtc-036"></a>DTC-036 | Todo documento estructurado propio declara su clase y `schema_version`: `config_kind` distingue configuración general y escenario, `manifest_kind` distingue los manifiestos, `effective-config.json` usa el discriminador equivalente de su esquema y la salida JSON de diagnósticos declara `diagnostic_kind`. La primera versión de esquema es la cadena exacta `"1"`; las posteriores se añaden por extensión (F2 la `"2"`, F3 la configuración efectiva `"3"`, F4 la `"4"` de la configuración efectiva y del manifiesto general, y las repeticiones en tándem la `"6"` de ambos). Resolución de DTPD-003; REQ-CONF-01/10/27/111 a REQ-CONF-115, REQ-MAN-02/69/73 y REQ-SAL-03/04 | Impedir que un documento se interprete con la clase o semántica de otra versión y separar la evolución del esquema de la versión del programa o del dialecto de validación | El programa solo lee versiones expresamente conocidas, y hoy son: configuración efectiva `"4"` y `"6"` (`effective_schema_version`); manifiesto general `"4"` y `"6"` (`is_supported_general_manifest_version`); manifiesto de evento `"2"` y `"3"`; TOML general `"2"` y `"4"` (el `"4"` es la extensión que puede declarar repeticiones; el despacho admite además el `"3"` de F5, y el escenario F5 el `"1"`); diagnósticos `"1"`. Los lectores de las versiones `"1"` y `"2"` de la configuración efectiva y del manifiesto general, y el de la `"3"` de la efectiva, se retiraron al cerrar F4.1 (DTC-082), antes de esta revisión; el código muerto que quedaba de la `"1"` se borró el 11/10/2026. Una versión ausente, posterior o desconocida se rechaza. No hay migración automática ni modificación silenciosa; el programa escribe siempre la versión actual de cada clase. Versión del programa, `schema_version` y JSON Schema Draft 2020-12 son conceptos distintos. Los esquemas de las versiones retiradas se conservan en `schemas/` como ficheros congelados, sin lector activo: el campo técnico F3 no era una extensión inocua del contrato cerrado de la efectiva `"2"`. DTC-044 concreta la clase y el esquema del documento de diagnósticos |
| <a id="dtc-037"></a>DTC-037 | La aleatoriedad de v1 usa generadores ChaCha12 nombrados y locales, nunca `StdRng` ni un flujo global. La semilla raíz y cada subsemilla tienen 256 bits y se representan en TOML, configuración efectiva y manifiestos como exactamente 64 hexadecimales minúsculos. Si falta la raíz, se obtienen 32 bytes del CSPRNG del sistema operativo antes de simular y se registran valor y origen. Resolución de DTPD-004; REQ-CONF-26/27, REQ-CONF-101/111 a REQ-CONF-115 y REQ-EVO-12/62 | Aislar biología, muestreos y eventos para que el consumo de un proceso o la planificación no alteren decisiones ajenas y permitir paralelismo reproducible sin compartir estado mutable | BLAKE3 deriva cada subsemilla con el contexto ASCII exacto `genome_population_simulator/v1/subseed`. El material binario canónico contiene, en este orden, la raíz cruda de 32 bytes, la versión de derivación como `u32` big-endian, el número de componentes como `u32` big-endian y, para cada componente en el orden normativo del catálogo de ámbitos, nombre de campo ASCII y valor UTF-8, ambos precedidos por su longitud `u32` big-endian. La identidad de flujo nunca procede de posición accidental, iteración, hilo o planificación. Se separan al menos población inicial, generaciones y procesos biológicos, muestreos y cada evento técnico, afinando por unidades independientes. Cada consumidor solo avanza su generador. Se registran raíz, origen, ChaCha12, versión de derivación y subsemillas de evento; los ámbitos internos se reconstruyen desde IDs y reglas versionadas. Ninguna fila publicada de un catálogo, sus campos/orden ni sus vectores se reescribe: añadir ámbitos exige una nueva versión de catálogo y puede conservar la versión de derivación si codificación/contexto no cambian; alterar una fila existente exige además nueva versión de derivación. Cambiar ChaCha12, BLAKE3, esta codificación o un algoritmo de muestreo cambia la reproducibilidad. DTC-038 concreta la garantía; DTC-039 separa el dominio de IDs; DTC-043 rige dependencias; DTC-048 fija el lock; DTC-040 reserva SHA-256 para integridad publicada |
| <a id="dtc-038"></a>DTC-038 | V1 declara por separado cuatro niveles de reproducibilidad: biología; decisiones y lecturas técnicas; equivalencia semántica de salidas; e identidad binaria. Resolución de DTPD-005; REQ-CONF-26/27 y REQ-CONF-111 a REQ-CONF-115, REQ-REP-22/23 y REQ-MAN-69 | Permitir repetir y comparar cada clase de resultado sin deducir una garantía más fuerte desde otra más débil | La biología se garantiza con la misma versión del simulador, plataforma inicialmente admitida, semilla raíz, configuración efectiva, contenidos de entrada, perfiles y esquema y algoritmos aleatorios versionados. El nivel técnico exige además la misma subsemilla de evento y versiones y condiciones de las dependencias participantes. La equivalencia semántica compara contenido lógico y, ante herramientas externas, exige sus mismas versiones, parámetros y condiciones verificadas. V1 no garantiza identidad binaria por defecto para ninguna clase; solo podrá declararse por formato y condiciones tras una prueba expresa. No se promete reproducibilidad entre versiones del simulador ni entre plataformas no verificadas. El manifiesto registra por nivel `garantizado`, `no garantizado` o `no aplicable`, condiciones, versiones, perfiles, herramientas y dependencias. Toda semilla generada se registra y puede reutilizarse. Se mantienen comparadores separados para verdad lógica, contenido técnico, semántica y bytes. DTC-042 fija la plataforma oficial inicial; crates, toolchain y formatos adicionales permanecen separados |
| <a id="dtc-039"></a>DTC-039 | Cada entidad separa un identificador estable de una etiqueta humana Unicode opcional y no única. Los identificadores aportados son ASCII, sensibles a mayúsculas y cumplen conceptualmente `[A-Za-z0-9][A-Za-z0-9._-]{0,63}`; no admiten espacios, rutas ni otros caracteres. En el dominio se usan tipos distintos como `IndividualId`, `EventId` o `ReferenceId`, nunca `String` genérico en una frontera. Resolución de DTPD-006; REQ-CONF-01/38/39/51/52/101, REQ-IND-25/26, REQ-RG-03/04/07, REQ-DAT-68 a REQ-DAT-71 y REQ-VER-10/12 | Evitar confundir categorías, etiquetas, rutas o ámbitos y mantener identidades reproducibles compatibles con formatos genómicos | Los IDs automáticos se derivan con BLAKE3 y el contexto ASCII exacto `genome_population_simulator/v1/id`, usando el material canónico de DTC-037 y el catálogo normativo de campos. Toman los primeros 128 bits del resultado, los codifican en Base32 RFC 4648 mayúscula sin relleno `=` y anteponen el prefijo ASCII de clase del catálogo. El ID de ejecución incluye como causas mínimas raíz y digest de configuración; los demás incluyen tipo e IDs causales estables según catálogo, nunca hilo, ruta, orden físico o reloj. Los catálogos de campos y prefijos y sus vectores son artefactos normativos de F1 y no cambian sin nueva versión de derivación. Toda colisión se detecta y hace fallar, también la de IDs aportados. **Excepción acotada F3.1.2:** `TemplateId` no usa BLAKE3: se compone exactamente como `TPL_<library_id>_<template_ordinal_decimal>` con decimal canónico sin ceros iniciales. `library_id` es único en el evento y el ordinal es único dentro de la biblioteca; por tanto el ID es único por construcción, conserva sintaxis y longitud compatibles con SAM, VCF y FASTQ. Por esa inyectividad demostrable, `TemplateId` no requiere registro ni detección de colisiones mediante estado retenido. Se respetan los ámbitos normativos: variantes locales al conjunto/VCF, `RG` único en la ejecución, lectura o plantilla única en el evento y mates con el mismo ID de plantilla y roles distintos; el ID de individuo se conserva en `SM`, VCF y manifiesto. Rutas y nombres físicos se derivan o escapan sin adquirir autoridad. El manifiesto relaciona ID opaco, clase, ámbito, causa y etiqueta. Los IDs automáticos permanecen compactos y válidos para SAM, VCF y FASTQ; no se fija aquí el layout de directorios |
| <a id="dtc-040"></a>DTC-040 | SHA-256 es la única suma de comprobación publicada por v1. Se identifica como `sha256` y su digest se representa mediante 64 caracteres hexadecimales minúsculos calculados sobre todos los bytes exactos del fichero, sin normalización. Resolución de DTPD-007; REQ-REF-17/18/21, REQ-MAN-39/76 y REQ-SAL-19 | Proporcionar una identidad de contenido interoperable y verificable para cada entrada, salida y fichero enlazado | Cada manifiesto registra algoritmo, digest, tamaño en bytes y ruta relativa de cada fichero. La suma se calcula o verifica al copiar entradas, al cerrar salidas, antes de exportar y después de exportar; copia y original deben coincidir. Puede calcularse incrementalmente durante la E/S, pero las comprobaciones exigidas no se sustituyen por tamaño, fecha, metadatos ni muestreo parcial. Cada fichero conserva su digest; una suma adicional del paquete no reemplaza las individuales. MD5 solo se emite cuando un estándar genómico lo exige, separado y nunca como integridad principal. BLAKE3 queda reservado a derivación de subsemillas e IDs con dominios distintos. DTC-045 concreta la publicación; la futura dependencia SHA-256 debe superar DTC-043 y quedar fijada por el lock de DTC-048 |
| <a id="dtc-041"></a>DTC-041 | V1 admite FASTA de texto sin comprimir y encapsula `noodles-fasta` detrás de un contrato propio como adaptador Rust puro para lectura y escritura por flujo, FAI y acceso indexado. El dominio valida al menos una secuencia no vacía, primer token de cabecera no vacío y único y alfabeto A/C/G/T/N en cualquier caja, sin huecos ni otros símbolos IUPAC. Resolución de DTPD-008; REQ-REF-01 a REQ-REF-18, REQ-REF-21/23, REQ-VAL-45, REQ-MAN-86 y REQ-SAL-01 | Mantener el contrato de referencia independiente de una biblioteca concreta y admitir ficheros grandes sin perder el original ni relajar las reglas biológicas | Un FASTA suministrado se copia byte a byte sin reescribirlo, se valida y calcula su SHA-256 durante el flujo, deja intacto el original y solo se publica tras éxito. Las bases se normalizan a mayúsculas únicamente en la representación interna; se conservan copia exacta y cabecera cruda, con primer token como ID y el resto como descripción. Una descripción no UTF-8 se conserva sin pérdida mediante bytes codificados reversiblemente y una marca de codificación en salida y manifiesto. Los FASTA simulados o derivados se escriben canónicamente sin compresión, con A/C/G/T/N mayúsculas, LF, ancho fijo y orden estable; el ancho se versiona si afecta a bytes. Los índices solo se crean cuando se necesitan. `samtools faidx` es validador independiente de pruebas. La compresión FASTA se diagnostica como no admitida. Evaluación final de F1 a 2026-08-10: `noodles-fasta` 0.66.0, MIT, MSRV 1.89.0 compatible con Rust 1.97.1 y edición 2024 queda fijada con características mínimas síncronas y sin `async`. Su grafo no compila ni enlaza código nativo; los cuatro `build.rs` transitivos de Crossbeam y Rayon fueron auditados y solo emiten configuración para Cargo, sin ejecutar procesos, usar la red, escribir fuentes ni generar o compilar código. La dependencia transitiva de BGZF incorpora compresión y paralelismo aunque F1 use FASTA sin comprimir; el coste se acepta y se vigila. La API upstream evoluciona con frecuencia; adaptador propio, versión exacta, lock y pruebas mitigan el riesgo conforme a DTC-043 y DTC-048 |
| <a id="dtc-042"></a>DTC-042 | La plataforma nativa oficialmente validada para v1 es Ubuntu 24.04 LTS x86-64 GNU, con triple `x86_64-unknown-linux-gnu`, configuración de CPU portable y sin `target-cpu=native`. Otras distribuciones x86-64 GNU pueden funcionar, pero no reciben garantía hasta incorporarse a una matriz. Ubuntu 24.04 sobre WSL2, siempre en el sistema de ficheros Linux/ext4 y no en `/mnt/c`, es el entorno de desarrollo admitido y actual, no una plataforma de publicación independiente. Resolución de DTPD-009; REQ-CONF-111 a REQ-CONF-115, REQ-REP-22/23, REQ-MAN-05 y REQ-SAL-19; hoja de ruta F1.1 y F9.7 | Disponer de un entorno rápido de desarrollo y una referencia de publicación verificable sin convertir el triple GNU en una promesa genérica ni heredar particularidades de WSL2 | Las pruebas ordinarias se ejecutan en WSL2 sin dependencias exclusivas. Antes de publicar se valida en CI Ubuntu 24.04 nativo y en un entorno OCI cuya base y receta se fijan por digest; se prueban atomicidad, publicación no-clobber en su sistema de ficheros, señales, permisos, recursos y procesos externos. V1 distribuye un `.tar.gz` versionado con ejecutable, esquemas, licencias y ejemplos, cuyo SHA-256 se calcula sobre los bytes exactos del archivo, y una imagen OCI de referencia identificada por el digest `sha256:` de su descriptor de índice o manifiesto, conservando y verificando los descriptores de configuración y capas. Un archivo que exporte la imagen OCI tiene además su propio SHA-256 independiente. Se registran plataforma, arquitectura, versión del simulador, referencia y digest de imagen y herramientas externas. macOS, Windows nativo y ARM64 no tienen soporte oficial; ampliarlos exige dependencias y matriz de pruebas. No se promete identidad binaria entre compilaciones. DTC-043 conserva dependencias, DTC-045 la publicación y DTC-048 el toolchain |
| <a id="dtc-043"></a>DTC-043 | V1 prefiere crates Rust puros para configuración, serialización, hashes, RNG y adaptadores de formatos cuando satisfagan correctamente el estándar. No enlaza bibliotecas nativas C/C++ salvo necesidad imprescindible evaluada por separado y aprobada expresamente. Las herramientas genómicas completas y consolidadas, como minimap2, se ejecutan como procesos externos detrás de adaptadores propios. Resolución de DTPD-010; REQ-CONF-69, REQ-DAT-21/22, REQ-MAN-34 a REQ-MAN-36 y REQ-VAL-11; hoja de ruta F1.1 | Reducir la superficie de distribución y acoplamiento sin reimplementar herramientas especializadas ni filtrar decisiones accidentales de una dependencia al dominio | Toda dependencia queda detrás de un contrato propio con tipos y errores del sistema. El adaptador de proceso controla versión, parámetros, comando, entradas y salidas, código de salida, stdout, stderr y fallos. Antes de admitir una dependencia se documentan propósito y requisito, corrección del estándar, mantenimiento, seguridad, licencia, rendimiento, estabilidad de API, características activadas, código nativo, scripts de construcción y coste de distribución. La mera presencia de `build.rs` no implica rechazo: cada script se audita por su función, entradas, salidas y efectos. Se rechaza si ejecuta procesos externos o accede a la red, modifica fuentes, genera o compila código nativo, o produce efectos no necesarios o no acotados, salvo decisión expresa que justifique y controle la excepción. Se minimizan dependencias y *features* y no se añaden por conveniencia anticipada. Las entregas fijan versiones y el manifiesto registra dependencias y herramientas efectivas; DTC-048 fija el toolchain y el lock de los crates aceptados. `noodles-fasta` permanece aprobado por DTC-041, sin prejuzgar otros crates. DTC-047 conserva la validación independiente y DTPD-105, la concreción de su adaptador de proceso externo |
| <a id="dtc-044"></a>DTC-044 | V1 usa errores internos tipados y jerárquicos, separados de una capa pública estable de diagnósticos. Sus categorías mínimas son configuración o esquema, entrada, formato, validación de dominio o biología, recursos, proceso o herramienta externa, E/S o publicación y defecto interno. Resolución de DTPD-011; REQ-CONF-09/14, REQ-FMT-04/07 y REQ-VAL-01 a REQ-VAL-51; hoja de ruta F1.2 | Ofrecer rechazos automatizables y mensajes accionables sin acoplar el contrato público a textos o tipos accidentales de dependencias | Todo diagnóstico público contiene código estable del catálogo versionado `GPS-{categoría}-NNN`, donde la categoría es `CONF`, `INPUT`, `FMT`, `DOMAIN`, `RES`, `TOOL`, `IO` o `INTERNAL`; además contiene categoría, mensaje accionable y contexto disponible: ruta y ubicación, campo e IDs de entidad, evento o artefacto. El catálogo debe enumerar los códigos de F1 antes de implementarla y ningún código se reutiliza ni cambia de significado. La cadena causal técnica queda para detalle o log, no como mensaje principal. Los fallos esperables retornan sin `panic`; los defectos internos se distinguen. Por defecto se renderiza texto humano a stderr. El modo JSON emite exactamente un documento RFC 8259 por invocación, nunca JSONL, validado por esquema y con `diagnostic_kind: "gps_diagnostics"`, `schema_version: "1"` y la lista `diagnostics`. Salidas: 0 éxito o solo avisos; 1 fallo de ejecución o dominio con configuración válida; 2 uso, configuración o entrada inválida; 3 recursos; 4 herramienta externa; 5 E/S o publicación; 70 defecto interno. Los contratos públicos no filtran tipos de crates; DTC-036 rige clase y versión y DTC-048 fija dependencias cuando se elijan |
| <a id="dtc-045"></a>DTC-045 | Toda ejecución iniciada de v1 trabaja en un directorio provisional hermano del destino, en el mismo sistema de ficheros y marcado `INCOMPLETE`. Cada artefacto se cierra, valida e integra primero desde su propio temporal. Un directorio final de éxito o fallo controlado solo aparece cuando sus manifiestos y catálogo de estados son coherentes y se publica mediante una operación Linux atómica no-clobber equivalente a `renameat2(..., RENAME_NOREPLACE)`, detrás de un adaptador propio. Resolución de DTPD-012; REQ-VAL-11, REQ-MAN-73/76 y REQ-SAL-05/19/21 a REQ-SAL-26 | Impedir que una ejecución parcial parezca completa o sobrescriba otra concurrente sin prometer durabilidad física | Se distinguen tres resultados. Un rechazo de admisión o pre-simulación no publica destino: solo deja provisional o diagnóstico externo conforme a REQ-SAL-26. Tras comenzar, un fallo controlado retira o ignora temporales incompletos, conserva y declara exportables los artefactos o conjuntos independientes ya validados, cierra y valida el manifiesto general con estado `fallido`, diagnóstico y catálogo de unidades, y publica atómicamente el directorio final como fallido, nunca completo. Una interrupción no controlable, como SIGKILL o caída, deja el provisional `INCOMPLETE`; conserva diagnóstico solo cuando fue posible escribirlo. Si el destino existe, la operación devuelve `EEXIST` y se rechaza sin sobrescribir, mezclar ni reparar; si el sistema de ficheros no soporta la garantía, se rechaza de forma identificable antes o durante publicación, sin fallback check-then-rename. V1 no reanuda ni reutiliza parciales. Al exportar a otro sistema de ficheros se copia a un provisional allí, se verifican tamaño y digest y se publica con la misma garantía no-clobber; no hay atomicidad entre sistemas. Se prueban concurrencia e interrupciones en WSL2/ext4, Ubuntu nativo y el sistema de ficheros OCI admitido según DTC-042. Publicación lógica y `fsync` o durabilidad física son garantías distintas |
| <a id="dtc-046"></a>DTC-046 | La línea base inicial usa tres escenarios sintéticos del mismo recorrido: A, aceptación de 100 kbp; B, crecimiento con referencias de 1, 10 y 100 Mb; y C, observación de 1 Gb manual o programada. **DTC-078 limita F1 a ejecutar A/B/C en el entorno local disponible; la comparación nativa de publicación y OCI pertenece a F9.** Resolución de DTPD-013; REQ-NOF-01 a REQ-NOF-08; hoja de ruta F0.4, F1 y F9 | Obtener evidencia comparable de corrección y crecimiento sin convertir medidas iniciales en promesas de producto | Se miden tiempo total y por etapa, CPU cuando sea fiable, RSS máximo, bytes de entrada, temporales y salida, MB/s y número de ficheros o artefactos. Se usa una construcción optimizada y portable, sin `target-cpu=native`, configuración y semilla fijas, salida limpia y repeticiones suficientes. F1 registra el entorno local efectivo y conserva la línea base; F9 repite la aceptación de publicación en Ubuntu 24.04 nativo y OCI. No hay inicialmente objetivos de producto ni fallos por pequeñas variaciones temporales: corrección e integridad siguen siendo obligatorias. Umbrales y presupuestos permanecen en DTPD-115. Estos tamaños son cargas sintéticas, no una promesa de escala biológica futura |
| <a id="dtc-047"></a>DTC-047 | Los validadores independientes de v1 son dependencias de pruebas, CI y aceptación de publicaciones; no se ejecutan sobre cada fichero de una ejecución normal ni requieren instalación por la persona usuaria. El runtime valida internamente los artefactos e invariantes que corresponden a su fase. **DTC-078 asigna a F1 únicamente los oráculos de TOML, FASTA, SHA-256 y JSON necesarios para sus salidas; listado/extracción de `.tar.gz` e inspección/arranque OCI se ejecutan en F9.** Resolución de DTPD-014; REQ-FMT-04/07, REQ-VAL-01 a REQ-VAL-51 y REQ-SAL-05; hoja de ruta F1 y F9.5 | Obtener un oráculo distinto del productor sin convertir herramientas de prueba en dependencias operativas ni duplicar todo el coste de validación en cada ejecución | FASTA se contrasta con `samtools faidx`, SHA-256 con `sha256sum` y documentos JSON con una implementación JSON Schema distinta. En F9, `.tar.gz` se contrasta mediante listado y extracción estándar y la imagen OCI mediante inspección y arranque de prueba. En sus fases, VCF se contrasta con bcftools y BAM/CRAM con samtools u otra referencia independiente. Los casos incluyen válidos, truncados, corruptos, campos desconocidos y errores de borde; CI fija y registra versiones y plataforma. El mismo escritor o lector nunca es el único oráculo. Validez formal y corrección biológica se prueban por separado. Una auditoría externa opcional en runtime queda fuera de v1. DTC-048 rige la fijación de versiones y las DTPD de formatos posteriores, sus validadores concretos |
| <a id="dtc-048"></a>DTC-048 | El toolchain oficial inicial es exactamente Rust 1.97.1 stable, verificado como estable actual el 16 de julio de 2026. Al crear F1, `rust-toolchain.toml` fija `channel = "1.97.1"`, `profile = "minimal"`, los componentes `rustfmt` y `clippy` y el target `x86_64-unknown-linux-gnu`; la edición `2024` se declara en cada `Cargo.toml` o en el workspace, no en el fichero de toolchain. `Cargo.lock` se versiona para la aplicación y CI y entregas usan `--locked`. Resolución de DTPD-015; REQ-CONF-111 a REQ-CONF-115, REQ-MAN-05 y REQ-REP-22/23; hoja de ruta F1.1 | Hacer reproducible y auditable la construcción sin declarar una compatibilidad mínima que no se prueba ni trasladar Rust a las personas usuarias del binario o imagen | No existe otra MSRV inicial: el compilador admitido es el fijado. Manifiesto y metadatos registran rustc y Cargo, edición, lock o resolución efectiva y herramientas, sin prometer identidad binaria. Las versiones concretas de crates se eligen al implementar F1 después de comprobar compatibilidad y DTC-043 y quedan fijadas por el lock; una candidata evaluada no cuenta aún como fijada. Toda actualización de Rust o dependencias es deliberada y aislada y exige pruebas, validadores independientes, vectores de RNG y reproducibilidad y escenarios de rendimiento; después se vuelven a declarar garantías sin asumir equivalencia. Rust 1.97.1 sustituye la propuesta accidental basada en Rust 1.94.0 local. Las personas usuarias de binario o imagen no necesitan Rust |
| <a id="dtc-049"></a>DTC-049 | Los valores numéricos configurables de probabilidad o frecuencia y sus tolerancias generales se representan como `f64`. Todo valor y tolerancia debe ser finito; una probabilidad o frecuencia pertenece al intervalo cerrado `[0,1]`. La tolerancia absoluta general es configurable en el intervalo cerrado `[0, 1e-6]` y vale `1e-9` por defecto. Resolución general de DTPD-110, completada para F5 por DTC-089; REQ-CONF-106/107, REQ-VAL-40, REQ-MAN-81 y ACE-105 | Obtener una validación y normalización únicas, auditables y reproducibles sin convertir el orden físico o el paralelismo en parte accidental del resultado | Los valores se recorren en orden canónico estable por identidad y categoría y se suman con Neumaier, sin reducción paralela. Se calcula `delta = abs(suma - 1)` y el borde inclusivo acepta `delta <= tolerancia`. Si se acepta y la suma es positiva, finita y normalizable, cada valor se divide por la suma y se recomputan suma y delta con el mismo algoritmo y orden; si excede la tolerancia, la suma no es positiva o finita, la normalización no produce valores válidos o una entrada es inválida, se rechaza antes de simular. Configuración efectiva y manifiesto registran tolerancia efectiva, orden canónico, algoritmo y versión, suma y delta, y valores antes y después. Toda serialización decimal debe permitir recuperar exactamente el mismo `f64`; las garantías de reproducibilidad exigen conservar estas condiciones. Es una validación determinista de configuración: no altera DTC-024, no sustituye `0,5/C` ni ninguna tolerancia específica y no selecciona crate |

| <a id="dtc-050"></a>DTC-050 | V1 incorpora en F1.6 una referencia reducida derivada de un FASTA real ya admitido por F1.5. Usa exactamente un objetivo por bases o porcentaje, selección reproducible sin reemplazo ni solapamiento, orden relativo/orientación conservados, como máximo un residual, `N` literal y publicación conjunta de FASTA canónico, mapa completo, SHA-256 y manifiesto. REQ-CONF-124 a REQ-CONF-127, REQ-REF-24 a REQ-REF-36, REQ-VAL-52/53, REQ-MAN-89/90 y REQ-SAL-27 | Decisión funcional explícita del 10 de agosto de 2026, completada el 14 de agosto de 2026 | `reduction_percent` expresa el porcentaje eliminado, pertenece a `(0,100)` y retiene `floor(S × (100 − p) / 100)`. El parser conserva el literal decimal TOML y efectúa el cálculo racional con aritmética comprobada; convertir antes a `f64` está prohibido. Todo objetivo debe ser menor que `S`. Cromosomas y *scaffolds* son elegibles, conservan ID y clase y nunca se mezclan. DTC-056/DTC-057 resuelven coordenadas y selección. La reducción no es captura dirigida ni variante estructural biológica |
| <a id="dtc-051"></a>DTC-051 | Toda capacidad configurable de v1 debe estar disponible en Marimo además de TOML/CLI. La aplicación invoca `validate`/`run` sin shell; Rust, esquemas, `EffectiveConfig` y manifiestos conservan autoridad. F2 usa un TOML general y cada fase añade los controles de sus capacidades reales; los escenarios aparecen con su primer consumidor. REQ-UI-01 a REQ-UI-15 y REQ-VAL-54 a REQ-VAL-56 | Decisión funcional explícita del 10 de agosto de 2026, simplificada por DTC-079 | No autoriza Python dentro del núcleo, ejecución automática, reglas duplicadas, paquetes de resultados ni carga completa de artefactos grandes. DTC-079 concreta la interfaz local proporcional y F9 conserva distribución/endurecimiento. |
| <a id="dtc-052"></a>DTC-052 | Los SNV heredables iniciales de F2 proceden de una matriz de frecuencias común a la población, suministrada o generada desde número/densidad, espacio, sustituciones y un espectro alélico explícito. Los loci son marginalmente independientes; Hardy-Weinberg es el modelo predeterminado y `F` solo se admite en diploides. REQ-CONF-128 a REQ-CONF-130, REQ-IND-60 a REQ-IND-65, REQ-VAL-55, REQ-MAN-91 y REQ-UI-14 | Decisión funcional explícita del 10 de agosto de 2026, ampliada tras revisar espectros empíricos | Separa creación de la matriz y muestreo de genotipos; no usa densidades SNP independientes por individuo, no inventa LD y no convierte 1000 Genomes o ExAC en un perfil universal. En F2, HWE cubre P=1..10 y F solo P=2; una distribución genotípica explícita no tiene esquema/controles y cierra REQ-CONF-130 en F6.2. DTC-054 resuelve los formatos, algoritmos y perfiles de DTPD-118 |
| <a id="dtc-053"></a>DTC-053 | **Superada por DTC-089 para F5.** Conserva solo historia documental del panel externo; no es contrato ni entrada vigente. | Evitar retorno del modelo retirado | F5 sintético; F6 SV; F7.1 meiosis. |
| <a id="dtc-054"></a>DTC-054 | DTPD-118 se resuelve con una matriz canónica `population_snvs.tsv` + `population_snvs.metadata.json`, ligada por SHA-256 a una referencia y con `alt_af` como única frecuencia. Una orden prepara entradas tabulares y otra genera el mismo recurso reutilizable. La generación produce candidatos bialélicos sin reemplazo; usa espacio uniforme o regiones no solapadas ponderadas, 12 sustituciones dirigidas y un espectro no solapado fijo/lineal/logarítmico en modo estocástico o exacto por mayores restos. REQ-CONF-128/129, REQ-IND-60 a REQ-IND-65, REQ-VAL-55, REQ-MAN-91, REQ-UI-14 y ACE-116 | Decisiones explícitas aprobadas el 12 de agosto de 2026 | Poblaciones y subsemillas son independientes; 0 y 1 son válidos, se eliminan sin reemplazo los `ALT` globalmente ausentes y se conservan los fijados; una suma ALT menor o igual que 1 conserva `REF=1-suma`, un exceso dentro de tolerancia se normaliza a 1 y otro mayor se rechaza. La multialelia aleatoria y conversión VCF se aplazan a DTF-039/038. `uniform_snv12_v1` es el valor de sustitución predeterminado; los perfiles humanos aproximados de DTC-068 son ejemplos explícitos, nunca implícitos ni tasas universales. Genotipos incompatibles y matriz vacía solo pueden continuar mediante autorizaciones específicas trazadas, no `--force` general. DTC-049 rige `f64`, Neumaier y tolerancia |
| <a id="dtc-055"></a>DTC-055 | **Superada por DTC-089 para F5.** Paquete y recombinación fundadora retirados. | Evitar dos contratos F5 | DTC-089 rige; F7.1 define meiosis. |
| <a id="dtc-056"></a>DTC-056 | DTPD-101 se resuelve con una única representación interna base 0 semiabierta. Una posición de base, un punto interbase y un intervalo no son intercambiables; cada valor viaja etiquetado por `reference_id` y `sequence_id`. Configuración, mapa de reducción y demás contratos propios externos conservan base 1 cerrada, y cada formato ajeno mantiene su convenio oficial. REQ-REF-29/35/36, REQ-IND-38 a REQ-IND-40, REQ-VAL-53 y ACE-113 | Evitar desplazamientos silenciosos y disponer de un contrato común para reducción, variantes, segmentos, anclas y composición | Las magnitudes usan `u64` y aritmética comprobada. Un intervalo interno `[start,end)` exige `start < end` y longitud `end-start`; un punto interbase admite `0..=sequence_length`. La conversión externa `[s,e]`↔interna `[s−1,e)` es exclusiva de las fronteras, exacta y reversible; desbordamiento, referencia/secuencia incompatible o valor fuera de rango se rechaza sin recorte. Los DTO externos no atraviesan el dominio y F2 extiende el mismo servicio sin crear otra convención |
| <a id="dtc-057"></a>DTC-057 | DTPD-116 usa un selector directo proporcional. Reparte las regiones completas entre secuencias según longitud, con desempate FASTA y capacidad comprobada, y sitúa bloques no solapados por composición directa de huecos; el residual se reserva en una secuencia factible y queda último en el orden de selección. REQ-CONF-126/127, REQ-REF-24 a REQ-REF-36, REQ-VAL-52/53, REQ-MAN-89/90 y ACE-112/113 | Obtener tamaño exacto, coste proporcional al número de bloques y azar reproducible sin reintentos | Preflight rechaza una disposición inviable antes de consumir RNG. El mapa conserva intervalos, clase y orden; materialización y mapa se ordenan por FASTA/coordenada. Un vector congelado y reconstrucción independiente verifican el contrato. |
| <a id="dtc-058"></a>DTC-058 | **Superada para F2 por DTC-079.** Conserva la historia de la concreción inicial de Marimo, pero sus protocolos, almacén de sesión, endurecimiento y puertas experimentales no son obligaciones de F2. | Evitar que una decisión técnica histórica mantenga infraestructura sin consumidor | Las versiones efectivas quedan fijadas por `pyproject.toml` y `uv.lock`; F2 prueba instalación desde el repositorio, arranque local y recorrido con el binario. Distribución, plataforma oficial y endurecimiento se verifican en F9. |
| <a id="dtc-059"></a>DTC-059 | En F2 la longitud de cada inserción aleatoria y la de cada deleción aleatoria son parámetros exactos, independientes y configurables, con 1 pb como valor predeterminado para ambas. Cada valor debe cumplir `1 <= longitud < umbral_indel_estructural_efectivo`; no se sortea una longitud ni se introduce una distribución de longitudes. REQ-CONF-55, REQ-IND-40/41 y ACE-29 | Formalizar la simplificación funcional aprobada sin convertir el valor predeterminado en una longitud fija obligatoria | Configuración efectiva y manifiesto conservan cada longitud y su procedencia. Cero, un valor igual o superior al umbral o una distribución solicitada se rechazan antes de simular. Una ampliación probabilística queda fuera de F2 |
| <a id="dtc-060"></a>DTC-060 | `prepare-snv-matrix` obtiene la identidad y los bytes autoritativos de referencia desde un manifiesto general F1 validado: recibe `--reference-manifest <manifest.json>` y `--reference-id <id>`, resuelve exactamente la entrada `references[]` indicada y verifica su ruta, tamaño y SHA-256 antes de preparar la matriz. RAF2-019; REQ-CONF-128, REQ-IND-60, REQ-VAL-55 y REQ-MAN-91 | Evitar deducir una identidad estable desde una ruta FASTA o aceptar ID y bytes desacoplados | No admite un `--reference-id` aislado ni un FASTA alternativo. Rechaza ID ausente/distinto, mismos bytes bajo otro ID, mismo ID con otros bytes, manifiesto no válido y referencia derivada incoherente; una referencia derivada válida conserva su propio ID y procedencia F1 |
| <a id="dtc-061"></a>DTC-061 | DTPD-102 se resuelve con un estado interno híbrido y un fichero externo obligatorio `truth.jsonl` independiente de los tipos Rust. Cada copia interna mantiene una cadena ordenada, lógicamente inmutable y compartible de segmentos: los ordinarios referencian intervalos de referencia/fundador sin copiar bases; SNP e *indels* pequeños son modificaciones tipadas, y solo una inserción sin fuente almacena sus propios bytes. Segmentos/modificaciones conservan procedencia, copia/fase y coordenadas DTC-056; la secuencia completa se materializa solo bajo demanda. `truth.jsonl` —v4, o v5 cuando publica algún alelo de repetición; las v1 y v3 son históricas— usa UTF-8 sin compresión y JSON Lines versionado. La primera línea separa el SHA-256 de bytes del artefacto fuente en `source_digests` de `reference_content_sha256`, huella canónica de la referencia materializada. Cada copia conserva segmentos heredados que forman el estado inicial, los fragmentos vigentes tipados (`source` o `inserted`) y el historial cronológico de eventos activos en coordenadas vigentes; un indel estructural activo no duplica obligatoriamente una modificación heredada. Esta ampliación representa inserciones anidadas, deleciones que atraviesan bases insertadas y SNP eliminados sin forzar una coordenada de fuente inexistente. La matriz canónica es obligatoria y conserva también loci todo-REF; los genotipos SNP e indel se emiten por caso/locus con alelos tipados, dosis y copias; **desde el 2 de octubre de 2026 el locus de un caso indel está en coordenadas de la referencia** aunque su evento se declarase en coordenadas vigentes, con la regla de anclaje que fija `schemas/truth-v4.md`, y el validador semántico lo recalcula al replicar el historial. Una deleción es un alelo estructural `deletion`, nunca el texto centinela `<deleted>` dentro de un SNP. **Desde el 30 de septiembre de 2026 el alelo estructural enlaza su evento fuente y sus fragmentos vigentes y no copia la secuencia final de la copia portadora**, que se reconstruye desde `current_pieces` y la referencia: las versiones `"1"` y `"3"` la copiaban en `final_haplotype_bases`, lo que hacía crecer el fichero como *indels* × longitud del genoma, y pasan a ser históricas, sin lector; las vigentes son `"4"` (sin alelos de repetición) y `"5"` (con alguno publicado), descritas en `schemas/truth-v4.md`, con la misma regla de promoción por contenido. Lo que aquel campo contrastaba de rebote se conserva en memoria: al publicar, el productor compara para cada copia la huella SHA-256 de su materialización con la de la reconstrucción independiente, sin escribir nada en la verdad; **desde el 3 de octubre de 2026 el juez semántico de los casos genotípicos corre además siempre al publicar, con VCF o sin él**, en esa misma pasada (DT-019). **Desde ese mismo día exige además el conjunto completo de casos, y ninguno más**: por individuo, uno por fila de la matriz de su población aunque sea todo-REF, uno por cambio de base vigente de sus copias, uno por deleción del historial y por inserción que conserve bases, y uno por tracto de repetición realizado; dos casos con la misma fuente son un duplicado. El juez exige además los alelos exactos y en orden de cada caso SNP y de repetición, la población de cada copia y el identificador de cada caso, y que el estado heredado de las copias conserve la referencia entera. **Desde el 9 de octubre de 2026 la conserva una sola vez y en orden** (DT-048): en las versiones `"4"` y `"5"` los segmentos de cada copia, en su orden, recorren su secuencia de la referencia de principio a fin sin huecos, solapes ni saltos atrás; el productor nunca escribió otra cosa y la versión no cambia. Una versión futura con variantes estructurales (duplicaciones, inversiones, translocaciones; REQ-IND-46 a 49 y REQ-REC-20) tendrá que sustituir esa regla. `truth.jsonl` no lleva las filas de la matriz —van en su propio artefacto, ligado por SHA-256—, así que el juez las recibe al proyectar, de la matriz en memoria, y en el juez de los recursos publicados, que juzga los bytes del paquete con la matriz publicada y la población de la configuración efectiva, y exige ahí exactamente las copias de la simulación (individuos × secuencias × ploidía); un lector de solo `truth.jsonl` comprueba el resto sin ellas. Un documento recién proyectado que incumple el contrato estructural es defecto interno, no error de formato (DT-020, DT-021; el conjunto, en `schemas/truth-v4.md`). La verdad F5 `"2"` no cambia: nunca publica alelos estructurales (DT-014). REQ-IND-28/29, REQ-VER-27/31, REQ-SAL-02 y ACE-10/11/28/34/62 | Evitar duplicar genomas, mantener una verdad reconstruible e interoperable y escalar por flujo sin acoplar el contrato externo al contenedor interno | Los esquemas por tipo fijan claves y dominios; el orden canónico es referencia, individuo, secuencia y copia, seguido por las claves internas acordadas. La cabecera declara recuentos, digests de fuentes/grupos y la huella de contenido de referencia; el manifiesto conserva longitud y SHA-256 del fichero completo. Escritura/lectura son por flujo. Faltas, duplicados, desorden, recuentos/digests discordantes o relaciones rotas fallan. Un validador semántico/reconstructor independiente materializa secuencias y exige que cada genotipo coincida exactamente con sus haplotipos. VCF deriva de `truth.jsonl`, nunca al contrario. Recombinación/SV futuras componen el mismo modelo, con un contrato de la verdad que sustituya la colinealidad de los segmentos; compresión queda fuera de F2 y no cambia su lógica |
| <a id="dtc-062"></a>DTC-062 | La única versión VCF de salida de v1 es VCF 4.5 y su primera línea de cabecera es exactamente `##fileformat=VCFv4.5`. Se usan campos y representaciones estándar de VCF 4.5; cuando varias formas sean válidas se evita una novedad innecesaria si existe una representación ampliamente interoperable. No se producen alternativas VCF 4.2, 4.3 o 4.4. Una herramienta externa basada en htslib forma parte de la validación. Resuelve la versión dentro de DTPD-103; REQ-FMT-05/06, REQ-VER-06 y ACE-19 | Ofrecer un solo contrato de salida moderno sin multiplicar variantes de compatibilidad y conservar validación independiente | La preferencia interoperable no autoriza apartarse de VCF 4.5 ni inventar campos no documentados. DTC-066 integra esta decisión en el cierre de DTPD-103 |
| <a id="dtc-063"></a>DTC-063 | La salida VCF predeterminada es `.vcf.gz` comprimida con BGZF y lleva exactamente un índice: TBI si todas las secuencias y posiciones caben en su límite normativo, aproximadamente 2^29 bases por secuencia, o CSI si alguna lo supera. Nunca se publican ambos índices. Como alternativa explícita se admite `.vcf` de texto sin compresión y sin índice. Resuelve compresión e índice dentro de DTPD-103; REQ-FMT-05/06, REQ-SAL-12, REQ-MAN-58 y ACE-61 | Entregar una salida consultable e interoperable por defecto, sin imponer el límite de TBI a referencias grandes ni crear combinaciones ambiguas de artefactos | `EffectiveConfig` y el manifiesto registran modo, tipo de índice y rutas efectivas. VCF BGZF e índice se validan antes de hacerse visibles y se publican como una unidad atómica; el fallo de cualquiera impide publicar el conjunto. En modo texto se registra la ausencia deliberada de índice. DTC-066 integra esta decisión en el cierre de DTPD-103 |

**Concreción F2.5.5 (25 de agosto de 2026).** `vcf_output` es un bloque
opcional del TOML general y queda desactivado si no se declara, conforme a
REQ-SAL-06. Cuando se activa, `mode = "bgzf"` o `mode = "text"` es explícito;
la configuración efectiva conserva esa elección y el manifiesto enlaza el VCF
único, la referencia, muestras, modo y el único índice aplicable. Tamaño y
SHA-256 de VCF e índice se conservan exclusivamente en el catálogo físico.
| <a id="dtc-064"></a>DTC-064 | La proyección VCF aplica una normalización interna determinista y versionada a partir de los haplotipos finales de `truth.jsonl`: recorte de prefijo/sufijo común sin perder el anclaje requerido, alineamiento a la izquierda contra la referencia y orden por catálogo de referencia y coordenada. **Excepción, decidida el 23 de septiembre de 2026: los registros de repetición en tándem no se recortan ni se alinean.** La decisión S-09 exige publicarlos «con la secuencia completa en `REF` y `ALT`», y las dos reglas no pueden cumplirse a la vez sobre el mismo registro: recortarlos dejaba un `ALT` de una base anotado con los campos del estándar como si contuviera noventa, es decir, una anotación falsa. La excepción alcanza solo a esos registros; todos los demás siguen recortándose y alineándose igual. Consecuencia sobre el oráculo externo: `bcftools norm` **sí** recorta los registros de repetición, de modo que la comprobación «el normalizador externo no altera el fichero» los excluye, porque la representación buscada no es la canónica, y en su lugar una prueba propia comprueba tres cosas: que el juez externo **acepta** el registro, que lo **recorta de verdad**, y que recortarlo **no cambia ninguna secuencia de alelo**, reconstruyéndolas contra la referencia y comparándolas como multiconjunto. Esa prueba no afirma —ni puede afirmar, porque su única entrada es el fichero ya publicado— que los alelos publicados sean los que declara la verdad: **esa igualdad la garantiza el oráculo semántico**, que reconstruye los alelos desde `truth.jsonl` y los compara uno a uno. `bcftools norm` es solo un oráculo externo de prueba, no una dependencia de ejecución. En F2 existe como máximo un registro por `CHROM+POS`; todos los alelos finales compatibles se integran como `ALT` de un único registro multialélico y no se descomponen automáticamente. **Fusión de variantes solapadas entre copias, decidida el 30 de septiembre de 2026 y precisada el 2 de octubre.** Dos variantes **se solapan** si sus intervalos de `REF` en el VCF, `[POS, POS+len(REF))` con la base de ancla de los *indels* incluida, comparten al menos una base de referencia; la relación se cierra por transitividad —si A pisa a B y B pisa a C, sale un solo registro aunque A y C no se toquen— y lo que solo es adyacente no se funde. Las variantes solapadas de copias distintas forman **un único registro complejo** que cubre la unión de sus intervalos: cada `ALT` es el alelo final completo de una copia sobre ese tramo, obtenido aplicando a la referencia los alelos ya normalizados de esa copia, y cada copia recibe en `GT` el índice de su alelo, de modo que ninguna copia queda anotada con un alelo que no tiene y no hace falta el alelo `*`. Consecuencias de la definición: una inserción dentro del tramo borrado de una deleción, o en su borde derecho —su ancla es la última base borrada—, se funde con ella; una en el borde izquierdo comparte ancla y `POS` y ya era el mismo registro multialélico; una situada antes del ancla de la deleción solo es adyacente. Un registro fundido se recorta y alinea a la izquierda **entero**, como un registro multialélico de `bcftools norm`: mientras todos sus alelos acaben en la misma base se quita esa base, si alguno se vacía se antepone a todos la base anterior de la referencia, y al final se quita el prefijo común mientras todos conserven dos bases; si al alinearse pisa otro registro, ambos se funden y se reconstruyen hasta que ninguno pisa a otro. Un grupo sin tracto cuyos alelos comparten `POS` no cambia. **Antes de fundir, la normalización de cada copia** recompone en un solo alelo los componentes de esa copia que chocan: si el `REF` normalizado de un alelo pisa el de otro alelo de la copia, si pisa un componente en el que **esa misma copia** tiene alelo —un tramo que la copia conserva como la referencia no la obliga a nada—, o si al alinearlo atraviesa entero otro componente. Y **cada inserción pertenece a una sola ventana**: al anclar un *indel* la ventana gana una base por la izquierda, pero no la inserción que la copia tenga en ese punto, que es de la ventana de su propio componente; hasta el 2 de octubre de 2026 esa inserción salía dos veces. **Con un microsatélite**, decidido el 2 de octubre de 2026: las variantes que pisan un registro de repetición se funden en él y el registro fundido **hereda la excepción** de arriba —tracto entero, sin recortar ni alinear—; un grupo con tracto se construye así aunque sus alelos compartan `POS`, y entonces sale el mismo registro que el agrupador construía. **Un alelo cuya alineación libre acabaría dentro de un tracto que alguna copia realiza** se alinea de nuevo con el límite en la última base del tracto, y así se funde con su registro; **uno que cruza el tracto entero lo cruza**, y conserva la representación de siempre. Pero si ese alelo, escrito delante de un tracto del registro, tenía su tramo original detrás de él, **el registro de repetición llega hasta el final de ese tramo** (añadido el 2 de octubre de 2026 tras la exploración final): sin eso, el registro acababa con el tracto, el alelo quedaba escrito dentro de él y el `ALT` de la copia que conserva el tracto llevaba menos bases de repetición de las que su anotación, por posición, cuenta —un `(CA)7` anotado con 14 bases en un `ALT` que solo llevaba 10—. Si el ancla de un *indel* es una base de un tracto que la propia copia realiza, la ventana arranca en el inicio del tracto, porque esa base no es la de la referencia. **Anotación de repetición de cada `ALT`**, corregida el 2 de octubre de 2026 tras la segunda revisión adversarial: se anota **por posición, tracto a tracto**. Por cada tracto del registro, la copia que lo realiza lleva los bloques que la verdad declara, y la que conserva el de la referencia —que le llega intacto, por la regla del microsatélite invadido— lleva **la estructura en bloques del tracto de la referencia**, con todas sus bases. Se anota lo que la copia lleva en el sitio del tracto, nunca lo que sobra del `ALT`: un cambio de base o una inserción en el flanco no son repetición. La versión anterior derivaba el bloque de las bases sobrantes del `ALT` y contaba en él el flanco y la variante: un `(CA)9` con un cambio de base al lado se anotaba con 19 bases en vez de 18, y un `(CA)15` con otro, con 32 en vez de 30. **Solo se deriva de las bases** cuando dos copias tienen **las mismas bases y anotaciones distintas**; entonces las dos llevan la derivada, porque una anotación describe lo que está escrito. Derivar es quitar del `ALT` lo que coincide con el `REF` por delante —como mucho las bases del registro anteriores al tracto— y por detrás —como mucho las posteriores—, y en lo que queda contar, con el modelo de las repeticiones, solo las bases que de verdad repiten el motivo: un bloque **de la primera a la última copia completa del motivo**, con sus interrupciones interiores, o ninguno si no hay ni una copia completa. Con un locus **compuesto** —cuya estructura en bloques no se conoce de un solo motivo—, o con varios tractos en un registro donde haya que derivar, la ejecución se rechaza con un diagnóstico de configuración que nombra el locus. La regla la aplican el productor y el oráculo semántico por separado; el oráculo rechaza además dos registros cuyos intervalos se pisan, igual que dos con el mismo `CHROM+POS`, y **reconstruye cada copia** aplicando el `GT` de los registros a la referencia para compararla con la que describen las piezas de `truth.jsonl`. **Para descomponer un registro fundido** en registros simples hay que usar `bcftools norm -a`, que escribe el alelo `*` donde una copia no tiene la base; `-m-` separa los `ALT` en registros bialélicos y vuelve a dar `0` a la copia que no lleva ese `ALT` aunque no tenga la base, que es justo la anotación falsa que la fusión evita. Resuelve normalización/agrupación dentro de DTPD-103; REQ-FMT-05/06, REQ-VER-15/16/17/30/44 y ACE-25/29/30 | Derivar una representación canónica repetible de la verdad haplotípica sin introducir una segunda autoridad ni perder la relación entre alelos del mismo sitio | Ante cambios solapados se reconstruye primero el alelo final de cada haplotipo y se usa un intervalo `REF` común suficiente para representar alelos de longitudes distintas antes de recortar y alinear. La configuración efectiva v6 (`vcf_output.normalization_version`) y el manifiesto general v6 (`vcf_outputs[].normalization_version`) registran la versión del algoritmo, `vcf-small-allele-normalization-v2` desde la fusión del 2 de octubre de 2026, y sus lectores rechazan cualquier otra; las versiones 4 de ambos documentos, ya cerradas, no la llevan. Un escenario sin microsatélites publica la configuración efectiva v6, si su documento general es v4, pero el manifiesto **sigue en v4**: la versión del algoritmo va entonces solo en la configuración efectiva. Hasta esa fecha la versión existía como constante del código y **no se publicaba en ningún documento**, contra lo que esta celda decía. La restricción de un registro por posición se limita a SNP e *indels* pequeños de F2: SV/translocaciones futuras podrán requerir varios BND relacionados, incluso en una misma posición, conforme al estándar. DTC-066 integra esta decisión en el cierre de DTPD-103 |
| <a id="dtc-065"></a>DTC-065 | Todo registro VCF declara y usa `INFO/GPS_CASE` con `Number=.` y `Type=String`. Contiene los `case_id` originales únicos, en orden canónico, de todos los casos de `truth.jsonl` que aportan algún ALT; no se crea un caso compuesto y un caso aparece una sola vez aunque aporte varios ALT. `truth.jsonl` conserva para cada alelo final la relación ALT↔lista ordenada de casos fuente; si un ALT procede de varios casos, conserva todos. Un caso que abarca varios registros repite su ID en cada uno. **Desde el 3 de octubre de 2026, un caso que no aporta ningún ALT porque en cada copia portadora su efecto neto es nulo —la copia es la referencia alrededor de su locus— no tiene registro propio** (excepción de REQ-VER-25): la verdad, que registra eventos, lo conserva, y el VCF, que registra secuencia, no tiene nada que publicar por él. Normalmente no aparece en ningún registro; **sí aparece** en el `GPS_CASE` de un registro en que su copia lleva otro ALT cuyo tramo compuesto o fundido contiene el locus, porque la relación ALT↔casos se establece por tramo. El agrupador y el oráculo independiente solo admiten un caso sin ALT tras comprobar que su copia es la referencia ahí: una ventana compuesta en que esa copia queda igual a la referencia contiene el locus entero, o ningún cambio de la copia, leído de sus piezas finales, toca el locus, extremos incluidos. Un caso cuya copia difiere y no tiene ALT sigue siendo un defecto interno. Un registro fundido por solapamiento entre copias (DTC-064) lleva los `case_id` de todos los casos de todas las copias cuyos alelos contiene. Solo si el VCF pierde detalle añade `INFO/GPS_LOSS`, código cerrado/versionado y nunca texto libre, e `INFO/GPS_COMPLEMENT`, ID de la entrada enlazada. Resuelve RAF2-030 y `INFO`/complemento dentro de DTPD-103; REQ-VER-19/23 a REQ-VER-26, REQ-MAN-27/31, REQ-SAL-12 y ACE-33/108 | Mantener un VCF válido y una trazabilidad muchos-a-muchos explícita sin perder identidades originales, duplicarlas ni fabricar un caso compuesto | La cabecera declara los tres campos con número, tipo y descripción conformes a VCF 4.5. F2 prueba dos ALT/dos casos, un ALT/dos casos, un caso/varios registros, reordenación de entrada, duplicados e ID inexistente; como todos sus casos son representables, `GPS_LOSS`, `GPS_COMPLEMENT`, entradas y fichero permanecen ausentes. Cada conjunto genera como máximo un `vcf-complement.jsonl` versionado y solo ante pérdida real. La primera prueba positiva corresponde a F6.5 con un caso complejo/SV/BND real; allí se fijará la cardinalidad positiva multicaso del complemento sin alterar la lista `GPS_CASE`. Una pérdida literal descubierta en F2 puede probarse sin convertirse en puerta. DTC-066 integra esta decisión en el cierre de DTPD-103 |
| <a id="dtc-066"></a>DTC-066 | DTPD-103 queda resuelta integrando DTC-062 a DTC-065 con una puerta externa de aceptación/CI fijada a `bcftools` 1.24 y `tabix`/HTSlib 1.24. Estas herramientas nunca son dependencia de ejecución del simulador. VCF de texto y BGZF deben superar lectura completa, comprobación de normalización contra el FASTA autoritativo y, para BGZF, apertura/consulta regional del único TBI o CSI exigido por DTC-063. EBI `vcf-validator` no es autoridad porque su soporte publicado no cubre VCF 4.5. REQ-FMT-05 a REQ-FMT-07, REQ-VER-06/15 a 17/23 a 26/44, REQ-VAL-28, REQ-SAL-12 y ACE-19/25/30/33/61 | Cerrar un contrato VCF único con aceptación reproducible y dos validaciones independientes, sin convertir herramientas externas en productoras o reparadoras de la salida | Un oráculo semántico independiente, sin llamar ni compartir la proyección productora, lee FASTA, `truth.jsonl`, VCF y el `vcf-complement.jsonl` opcional; comprueba REF/ALT/GT, dosis, fase, muestras, orden, unicidad `CHROM+POS` en F2, ausencia de registros cuyos intervalos `REF` se pisen (DTC-064), que cada copia reconstruida desde el `GT` de los registros sea exactamente la de `truth.jsonl` (DTC-064), lista `GPS_CASE` canónica/sin duplicados y relaciones ALT↔casos existentes, ausencia F2 de `GPS_LOSS/GPS_COMPLEMENT`, cardinalidad caso↔registros y campos prohibidos. El corpus incluye dos ALT/dos casos, un ALT/dos casos, un caso/varios registros, orden de entrada, duplicado e ID inexistente, además de positivos texto, BGZF+TBI/BGZF+CSI y mutaciones individuales de cada invariante. La cardinalidad positiva del complemento se valida en F6.5, su primer caso complejo/SV/BND real. Las invocaciones exactas y la comparación sin reparación quedan fijadas en la resolución de DTPD-103 |
| <a id="dtc-067"></a>DTC-067 | DTPD-120 se resuelve mediante el algoritmo normativo `poisson-count-v1`. Admite `lambda` `f64` finita en `0.0..=2^52` (`4_503_599_627_370_496`), derivada conforme a DTC-022/026. `lambda=0` devuelve 0 sin consumir RNG; `0<lambda<10` usa inversión/Knuth y `lambda>=10` PTRS de Hörmann, con frontera inclusiva exacta. Cada uniforme abierto usa exactamente un `RngCore::next_u64()` del ChaCha12 local DTC-037 y `u = (((x >> 12) as f64) + 0.5) × 2^-52`. `exp`, `log`, `lgamma`, `floor` y `sqrt` proceden exclusivamente de `libm = =0.2.16`, nunca de matemática estándar. REQ-CONF-23 a REQ-CONF-25, REQ-IND-16 a REQ-IND-18, REQ-CONF-129 y ACE-10/109 | Obtener recuentos Poisson reproducibles entre Linux y WSL2, con consumo RNG y límites explícitos, sin hacer de la media un límite biológico | `lambda>2^52` se rechaza antes de consumir RNG con valor/límite/ámbito. Todo candidato fuera de `0..=2^53-1` o no representable exactamente como entero falla numéricamente; no hay cast saturante, recorte ni resorteo. Un recuento mayor que las oportunidades falla por la política existente, sin cambiar `lambda`. Se registran `algorithm_version`, `lambda`, ámbito/subsemilla y `sampled_count`. Vectores congelados cubren cero, subnormales, ambos lados de 10, medias representativas, *next-down*/2^52/*next-up* y máximo; un *hook* de prueba fuerza la rama rara fuera de rango. Oráculo ajeno y pruebas PMF/media/varianza —con `rand_distr` solo opcional/test-only— no reutilizan el productor. `libm` debe superar DTC-043, incluida licencia, auditoría y `build.rs` |
| <a id="dtc-068"></a>DTC-068 | DTPD-121 se resuelve con dos perfiles humanos compactos, aproximados, explícitos y nunca predeterminados. `human_1000g_global_sfs_approx_v1` resume los recuentos publicados finales de 1000 Genomes con pesos 64:12:8 para `[1/5008,0.005)`, `[0.005,0.05)` y `[0.05,1)`; los dos primeros usan `log_uniform`, el tercero `linear_uniform` y DTC-049 normaliza. El límite superior 1 pertenece a la cobertura, pero el muestreo continuo nunca lo devuelve; `alt_af=1` exige `fixed=1`. `1/5008` procede de una copia entre las 5.008 copias autosómicas de 2.504 muestras y pertenece a esta receta, no es un límite universal. `human_harris2017_snv12_approx_v1` deriva de los 96 contextos de Harris y Pritchard (2017): normaliza por separado AFR/AMR/EUR/EAS/SAS, promedia las cinco filas con el mismo peso, agrega contextos a seis clases colapsadas por hebra, refleja cada clase a su sustitución complemento inversa con mitad del peso y normaliza los 12 cambios dirigidos. REQ-CONF-129, REQ-IND-61, REQ-VAL-55, REQ-MAN-91 y ACE-116 | Ofrecer dos ejemplos pequeños, reproducibles y descriptivos sin descargar ni procesar el VCF poblacional, inventar perfiles exactos por superpoblación o presentar una aproximación como tasa humana universal | Fuentes fijadas: artículo final 1000 Genomes `https://www.nature.com/articles/nature15393`; artículo eLife `https://elifesciences.org/articles/24284`; datos fuente DOI `https://doi.org/10.7554/eLife.24284.004`; fichero `https://cdn.elifesciences.org/articles/24284/elife-24284-fig1-data1-v2.txt`, 3.644 bytes, SHA-256 `970cc71578fc008ef4e0d75c7c460ab874073da7c5b5226c526b387627c97a7f`. Artículo y datos eLife se atribuyen bajo CC BY. La receta pequeña versiona literales, parser, orden, normalizaciones, complementos inversos y SHA-256 de los recursos incorporados; F2.2.6 la regenera y contrasta con un oráculo independiente. F2 no depende de `bcftools` para estos perfiles ni descarga/procesa el VCF Phase 3. Perfiles suministrados por población siguen admitidos. Materialización exacta del VCF y perfiles empíricos exactos pasan a DTF-016/post-v1 |
| <a id="dtc-069"></a>DTC-069 | **Aplazada por DTC-079.** El conjunto multidocumento autocontenido se implementará cuando exista el primer consumidor real de escenarios. En F2, Marimo importa y exporta un único TOML general de la vía directa. | No anticipar un grafo de documentos, copias privadas y defensas TOCTOU que F2 no utiliza | El TOML debe seguir siendo válido para el binario, conservar IDs, orden y rutas relativas existentes y no sobrescribirse accidentalmente. La fase consumidora de escenarios recuperará los requisitos multidocumento que sigan siendo necesarios. |
| <a id="dtc-070"></a>DTC-070 | RAF2-029 se resuelve para el modo `exact` calculando primero las cuotas enteras por fila mediante mayores restos y desempate por el orden estable del espectro. Después, cada población permuta reproduciblemente el conjunto completo de `variant_id` canónicos con una subsemilla exclusiva de asignación y aplica las cuotas siguiendo esa permutación. Las poblaciones permutan independientemente y se prohíbe asociar filas por coordenada u orden canónico de variante. REQ-CONF-129, REQ-IND-61, REQ-IND-63, REQ-VAL-55 y ACE-116 | Conservar las cantidades exactas sin crear correlación espacial artificial ni acoplar poblaciones o el valor sorteado a la asignación de intervalos | El valor dentro del intervalo se obtiene después con un flujo RNG separado por población×`variant_id`. Una fila `fixed` participa en la asignación pero no consume RNG para el valor; una fila con cuota cero no consume RNG. Los oráculos demuestran cuotas exactas, estabilidad ante reordenación de entrada, ausencia de asociación con coordenadas, aislamiento entre poblaciones, consumo cero en ambos casos y separación asignación/valor. DTC-071 fija los nombres, componentes, orden y versión concretos de ambos ámbitos sin alterar esta semántica |
| <a id="dtc-071"></a>DTC-071 | RAF2-028 se resuelve con `subseed-scopes-v4.tsv`, catálogo normativo v4 y superconjunto estricto de v3. Conserva BLAKE3 y las seis filas/vectores v3 literalmente con derivación 1; añade once ámbitos F2 con derivación 2 y `scope_name` antes del número de componentes. `scope_derivation_versions` registra por fila, de forma legible por máquina, cuál de ambos materiales aplica. El sorteo de ALT usa `scenario_id`, `reference_id`, `sequence_id`, `position_1based` y `ref`, todos conocidos antes de ALT; `variant_id` solo se deriva después. En la vía directa sin documento de escenario, `scenario_id` es el literal interno reservado `"direct"`: no representa biología, no lo escribe la persona usuaria y ningún escenario real puede usarlo; queda registrado en `EffectiveConfig`, manifiesto y contexto RNG. REQ-CONF-111 a REQ-CONF-115, REQ-IND-63, REQ-VAL-55 y ACE-23/116 | Eliminar la dependencia circular de ALT y hacer que todo consumidor F2 tenga identidad estable, componentes ordenados y aislamiento verificable sin cambiar resultados F1 | Ruta, reloj, hilo, orden de iteración y posición accidental siguen prohibidos. Valores enteros de componentes usan decimal ASCII canónico sin signo ni ceros iniciales salvo `0`; IDs usan su forma validada. El catálogo y el manifiesto registran catálogo 4 y el conjunto de derivaciones `[1,2]`; cada ámbito enlaza algoritmo consumidor. Vectores BLAKE3 con raíz común y mutación aislada de cada componente demuestran que v1–v3 permanecen byte a byte inmutables |
| <a id="dtc-072"></a>DTC-072 | RAF2-033 se resuelve con la familia normativa de consumidores RNG v1 de §3.2: `open-unit-f64-v1`, `uniform-u64-below-v1`, `categorical-f64-v1`, `fisher-yates-v1`, `weighted-without-replacement-v1`, `dynamic-indel-selection-v1`, `linear-uniform-v1`, `log-uniform-v1`, `hwe-composition-v1`, `f-diploid-genotype-v1` y `phase-permutation-v1`. DTC-067 continúa siendo el único Poisson. REQ-CONF-111 a REQ-CONF-115, REQ-CONF-129, REQ-IND-16 a REQ-IND-18, REQ-IND-61 a REQ-IND-63 y ACE-23/107/116 | Fijar transformaciones, fronteras y consumo RNG para que cambiar una utilidad o la plataforma no altere silenciosamente matriz, indels, genotipo o fase | Cada algoritmo declara preflight, orden canónico y consumo cero; registra versión y ámbito. Todo intervalo es `[min,max)`: puede declarar `max=1`, pero nunca produce 1, y la corrección de redondeo lineal impide producir 0; solo `fixed` produce exactamente 0 o 1. La selección SNV es estática; la de indels reenumera el estado vigente después de cada evento aceptado. `log`/`exp` usan exclusivamente `libm = =0.2.16` ya sometido a DTC-043. Vectores literales, RNG instrumentado, enumeración exhaustiva pequeña y oráculos independientes preceden a todo consumidor |

**Concreción F2.3.4 de DTC-072 (23 de agosto de 2026).** La compatibilidad de
un SNP explícito no se compara con la realización RNG: se admite sin
autorización si el genotipo completo tiene probabilidad positiva según las
frecuencias de la matriz y HWE o `F` diploide; con `F=1` solo permanecen los
homocigotos de alelos con frecuencia positiva. `--allow-incompatible-genotypes`
solo conserva un caso de probabilidad cero y nunca alelos, referencia,
coordenadas, ploidía o la huella matriz–fase incompatibles. La fase enlaza la
identidad de referencia y SHA-256 de TSV canónico de la matriz antes de aplicar
eventos; los cambios `direct_events` se preflightan antes de abrir publicación.
| <a id="dtc-073"></a>DTC-073 | **Simplificada por DTC-079.** Marimo muestra y transmite únicamente opciones y autorizaciones implementadas por el binario F2 y descritas por su esquema/interfaz vigente. | Evitar anunciar opciones futuras sin crear un protocolo paralelo de capacidades | `--allow-empty-snv-matrix` y `--allow-incompatible-genotypes` aparecen cuando sus recorridos Rust existen; una edición o cambio de autorización invalida la validación. No se mantienen listas separadas `catalogued`/`available`. |
| <a id="dtc-074"></a>DTC-074 | RAF2-035 se resuelve con evidencia oficial de plataforma detectada, no autoetiquetada, en F9. El verificador registra `/etc/os-release`, kernel, arquitectura, triple Rust, marcadores WSL2, punto/tipo de montaje del árbol y proveedor/job/commit de CI; compara esos hechos con la clase solicitada y rechaza cualquier discordancia. Ubuntu 24.04 nativo exige `ID=ubuntu`, `VERSION_ID=24.04`, x86-64 y ausencia de marcadores WSL; WSL2 exige esos marcadores y un sistema de ficheros Linux fuera de `/mnt`. REQ-UI-13, REQ-CONF-111 a REQ-CONF-115 y ACE-115 | Evitar que dos ejecuciones en el mismo WSL2 se presenten como evidencia de dos plataformas distintas | F2 conserva solo actas locales de instalación/arranque y no afirma compatibilidad ni distribución oficial. F9 enlaza el job nativo identificable y el WSL2 detectado; ninguno sustituye al otro. La evidencia automática conserva versiones de Python/Marimo/uv/navegador y las actas manuales se separan. Una etiqueta de CLI es solo una expectativa que debe concordar con la detección |
| <a id="dtc-075"></a>DTC-075 | F2 conserva intervalos de frecuencia semiabiertos, separación entre selección estática de loci y dinámica de indels, y consumo explícito del recuento candidato. La división en un productor, expediente y paso independiente para cada primitiva queda superada por DTC-079. | Mantener la reproducibilidad científica sin convertir utilidades internas en productos separados | Los algoritmos pueden agruparse por módulo y probarse con vectores y casos de borde. Solo los contratos observables y los cambios capaces de alterar resultados requieren versión y trazabilidad propias. |
| <a id="dtc-076"></a>DTC-076 | **Superada por DTC-078.** Propuso que F1 definiera un `ResolutionDomain` sintético para contextos de fases posteriores. REQ-CONF-117/118/123 y REQ-MAN-88/90 | Conserva la historia de la revisión sin mantener una obligación prematura | No es operativa: cada fase propietaria introducirá sus relaciones reales junto con el primer consumidor |
| <a id="dtc-077"></a>DTC-077 | **Superada por DTC-078.** Propuso publicar en F1 una evidencia independiente de resolución, fijar globalmente IDs de fase y anticipar el orden de parámetros futuros. REQ-CONF-117/118/123 y REQ-MAN-88/90 | Conserva la historia de la revisión sin convertir el mecanismo provisional en contrato | `resolution-evidence` deja de pertenecer a F1; las fases propietarias fijarán evidencia, identidad y orden cuando existan datos reales |
| <a id="dtc-078"></a>DTC-078 | F1 adopta un modelo de confianza local mínimo. Valida como no confiables los TOML y FASTA aportados; después confía en los ficheros creados dentro de su provisional privado frente a manipulación activa y los protege frente a fallos propios, interrupciones y corrupción accidental con estados, tamaño, SHA-256, `INCOMPLETE` y publicación atómica no-clobber. La raíz TOML activa de F1 contiene solo clase, versión, semilla opcional, referencia y referencias reducidas. No admite todavía vía, ploidía, reglas, entidades poblacionales, eventos ni escenarios. La composición simulada A/C/G/T usa un validador específico con tolerancia efectiva fija `1e-9`, sin conservar el sistema genérico de normalización poblacional. El cierre estricto de F1 comprende 54 REQ completos (4/1/6/22/21 en F1.2–F1.6); cualquier comprobación F1 de un REQ reasignado es evidencia parcial, no cierre. | Separar reproducibilidad científica e integridad accidental de un modelo de atacante local activo, y evitar construir biología futura sin consumidor | Población↔subpoblación, fase↔generación, locus↔posición, vía, ploidía, reglas, eventos técnicos, normalización poblacional y sus manifiestos se implementan en sus fases propietarias. `ResolutionDomain` sintético, `resolution-evidence`, validación exhaustiva del árbol contra symlink/hardlink/inodos, ventanas TOCTOU y revalidación combinatoria del paquete salen de F1; el endurecimiento frente a atacante local activo se evalúa en F9/CI según el despliegue. Se conservan rutas relativas canónicas, permisos privados razonables, ausencia de sobrescritura y validadores externos pequeños de TOML, FASTA, SHA y repetibilidad. Los IDs exactos y sus destinos pertenecen a la matriz de trazabilidad vigente. |
| <a id="dtc-079"></a>DTC-079 | F2 adopta una implementación proporcional al uso local: repositorio, persona usuaria y provisional privado se consideran confiables frente a manipulación activa. Marimo configura el TOML general vigente, invoca directamente `validate`/`run` y muestra diagnósticos y resúmenes. F2 no crea escenarios sin consumidor, un ejecutor general de evidencias, protocolos paralelos de capacidades, almacenamiento persistente de sesiones ni defensas contra symlinks, inyección deliberada del entorno o carreras hostiles. El aprovisionamiento inicial explícito con `uv sync --frozen` puede usar Internet para descargar las dependencias resueltas por `uv.lock`; el arranque posterior y todo uso normal son locales, sin red ni servicios remotos. | Aplicar a toda F2 la simplificación que permitió cerrar F1 y evitar que la revisión técnica expanda el modelo de amenaza | Las entradas externas siguen validándose; Rust conserva autoridad semántica, reproducibilidad, publicación no-clobber y pruebas de invariantes biológicas. Cada hito usa pruebas Cargo y validadores de formato pertinentes, con una evidencia de cierre sencilla. La instalación distribuible sin descarga previa, `.tar.gz`, OCI, plataforma oficial, auditoría completa de accesibilidad y endurecimiento del despliegue quedan en F9. DTC-058/069/073 se interpretan según esta decisión; DTC-075 mantiene semántica pero no su partición documental. |
| <a id="dtc-080"></a>DTC-080 | Si un cambio directo elimina la coordenada de un locus SNV de la matriz en una copia fundadora, F2.3.5 la registra como `missing_by_deletion`. Sus recuentos y frecuencias REF/ALT usan solo las copias que aún cubren el locus (`observed_copies`); la deleción no se rechaza ni se corrige. | Hacer observable la variación estructural ya admitida por F2.3.4 sin inventar un alelo ni alterar el sorteo inicial. | El resumen conserva `total_copies`, `observed_copies`, `missing_by_deletion`, frecuencias objetivo, recuentos/frecuencias realizados y genotipos completos observables. Los genotipos con una copia ausente no entran en su histograma; la procedencia declara esta regla versión 1. |
| <a id="dtc-081"></a>DTC-081 | F4 introduce el catálogo RNG `f4-technical-v1`, derivado con DTC-037 y ChaCha12 local por decisión. Sus componentes se codifican exactamente en el orden `event_id`, `library_id`, `stage`, `template_ordinal` o `pre_hop_instance_id`, y después los específicos de la decisión; `stage` es ASCII y distingue `fragment`, `chimera`, `quality`, `error`, `contamination_percent`, `contamination_donors`, `contamination_remove`, `index_hop` e `final_id`. La cobertura y el fragmento pertenecen a una plantilla estable; calidad y error añaden `mate` y `cycle`; contaminación añade `recipient_id` y la clave candidata; *index hopping* usa exclusivamente biblioteca anterior al salto e identidad pre-salto. Todo sorteo se reinicia desde su clave: no existe consumo acumulado entre decisiones. `bernoulli-v1` no sortea para 0 o 1 y, en otro caso, compara un uniforme de 53 bits en `[0,1)` con `p`; `uniform-index-v1` usa rechazo sobre `u64` hasta obtener un intervalo cuyo cardinal sea múltiplo de `n`; listas sin reemplazo se ordenan canónicamente por su clave y aplican Fisher–Yates parcial con ese uniforme. Los casos sin candidatos o con cardinalidad imposible se rechazan antes de sortear. Cada cambio de fila, orden de componentes, algoritmo o vector crea una nueva versión de catálogo y, si altera la codificación DTC-037, una nueva versión de derivación. | Cerrar reproducibilidad, regeneración bajo demanda y paralelismo de F4 sin mantener estados RNG ni colecciones de lecturas | F4 publica en configuración efectiva/manifiesto el nombre y versión del catálogo, y prueba vectores literales, permutación de bibliotecas, ejecución paralela y ramas 0/1. La distribución cuantitativa del perfil y su normal truncada siguen siendo responsabilidad de DTPD-106 y no se inventan aquí. |
| <a id="dtc-082"></a>DTC-082 | F4 sustituye de forma atómica los cuatro contratos que consume al cerrar F4.1.3: TOML general `"2"`, `effective-config.json` `"4"`, manifiesto general `"4"` y manifiesto de evento `"2"`. F4.1.1 y F4.1.2 solo preparan el recurso, perfil e identidades internos; mientras F3 siga siendo el único productor publicado conserva exclusivamente sus contratos v1/v3. Al activar PE real, el runtime admite solo las versiones F4, los documentos y *fixtures* v1–v3 se actualizan o eliminan, sin migración automática ni lector heredado. | Evitar publicar un perfil PE sobre FASTQ SE y, al mismo tiempo, evitar matrices de compatibilidad que no aportan valor científico ni operativo | La transición ocurre junto al primer productor PE verificable. Entonces los cuatro escritores emiten sus versiones activas y los lectores rechazan versión ausente, anterior, posterior o desconocida. Un TOML anterior no puede activar F4. Diagnósticos y formatos ajenos a esos cuatro contratos no cambian por esta decisión. |
| <a id="dtc-083"></a>DTC-083 | F4 acepta exclusivamente UDI explícitos por biblioteca: `i7` e `i5` son cada uno ocho bases `ACGT`, en orientación de secuenciación, y cada componente es único dentro de la tanda. La cabecera FASTQ escribe literalmente `i7+i5`, sin complementar ni invertir `i5`. Ausencia de UDI, marcador `0`, índices simples o combinatorios, longitudes/alfabetos inválidos y duplicados se rechazan antes de generar. | Fijar un contrato pequeño, inequívoco y verificable para la única familia de perfiles F4, sin construir soporte anticipado para modos no usados | Un *index hopping* transfiere solo al UDI completo ya declarado de la biblioteca final; no crea combinaciones nuevas. Una combinación inesperada en la frontera se descarta defensivamente y solo incrementa un recuento agregado. `0` queda reservado para un perfil futuro no indexado. Añadir otro formato o longitud requiere nueva decisión, esquema y pruebas. |
| <a id="dtc-084"></a>DTC-084 | Una quimera F4 conserva el prefijo de longitud `b` de su plantilla de inserto `L`, con `b` uniforme en `1..L-1`, y concatena un segundo intervalo de longitud `L-b`. Sus candidatas son los intervalos válidos del mismo individuo/biblioteca, ordenados por copia, secuencia e inicio; en la copia fuente se excluyen los que solapan o tocan el prefijo. Se elige uniformemente una candidata, ambos segmentos se leen en orientación directa y la disposición PE se aplica al inserto compuesto. | Convertir la noción de «dos segmentos no contiguos» en un algoritmo completo y reproducible sin inventar un modelo de quimeras entre muestras | La selección usa las claves de quimera vigentes de DTC-086; una candidata ausente para una pareja seleccionada falla antes de publicar, sin reintento ni reducción silenciosa de la tasa. Solo se retienen coordenadas/orientaciones/punto durante la construcción y validación; el manifiesto publica tasa y recuento agregados. |
| <a id="dtc-085"></a>DTC-085 | `f4-technical-v2` sustituye al catálogo F4 v1 aún no implementado y añade siempre el componente ASCII `decision` después de `stage`. Fragmento: `[event_id,library_id,fragment,template_ordinal,decision]`, con `insert_length`, `copy` o `start`. Quimera: `[event_id,library_id,chimera,template_ordinal,decision]`, con `bernoulli`, `breakpoint` o `second_source`; esta última añade el `breakpoint` ya sorteado. Calidad/error: `[event_id,library_id,quality|error,template_ordinal,mate,cycle,decision]`, con `quality` o `substitution`. Porcentaje contaminante: `[event_id,recipient_library_id,contamination_percent,recipient_id,decision=percent]`; donantes/retirada: `[event_id,recipient_library_id,contamination_donors|contamination_remove,recipient_id,decision]`; *index hopping*: `[event_id,pre_hop_library_id,index_hop,pre_hop_instance_id,decision]`, con `bernoulli` o `destination`; ID final: `[event_id,final_library_id,final_id,pre_hop_instance_id,decision=technical_pair_id]`. | Completar DTC-081 con claves implementables para cada sorteo y versionar su cambio de codificación | Cada clave usa ChaCha12 nuevo y el consumidor DTC-081 correspondiente; 0/1 no consumen, reintentos de la normal permanecen dentro de `insert_length`, y un candidato ausente falla antes de publicar. `f4-technical-v1` no se publica ni implementa; cualquier cambio posterior crea v3 y nuevos vectores. |
| <a id="dtc-086"></a>DTC-086 | `f4-technical-v3` sustituye a v2, ambas aún no publicadas, y es el único catálogo F4 normativo. Fragmento usa `[event_id,library_id,fragment,template_ordinal,decision]` con `insert_length`, `copy` o `start`; quimera usa `[event_id,library_id,chimera,template_ordinal,decision]` con `bernoulli`, `breakpoint` o `second_source`, que añade el `breakpoint` ya sorteado; *index hopping* usa `[event_id,pre_hop_library_id,index_hop,pre_hop_instance_id,decision]` con `bernoulli` o `destination`. Calidad es determinista: la curva y sus ajustes producen el Phred efectivo y no consumen RNG. Error usa `[event_id,library_id,error,template_ordinal,mate,cycle,decision]`, con `decision=bernoulli` para decidir el error y `decision=substitution` solo cuando ocurre; el primero usa `bernoulli-v1` y el segundo `uniform-index-v1` sobre las otras tres bases. El porcentaje contaminante usa `[event_id,recipient_library_id,contamination_percent,recipient_id,decision=percent]` y `uniform-closed-53-bit-v1`: para límites finitos `a <= b`, si `a=b` devuelve `a` sin sortear; en otro caso toma un entero uniforme `k∈[0,2^53-1]` mediante `uniform-index-v1` y devuelve `a+(b-a)×k/(2^53-1)`, incluyendo ambos extremos. Donantes/retirada usan `[event_id,recipient_library_id,contamination_donors|contamination_remove,recipient_id,decision]`, respectivamente `decision=select_donors` y `decision=select_recipients`. `TechnicalPairId` no es un sorteo ni una excepción: es un ID automático BLAKE3 de DTC-039 con clase `TechnicalPairId`, prefijo `TP_`, versión de derivación `2` y componentes canónicos `class=TechnicalPairId,event_id,final_library_id,pre_hop_instance_id`, en ese orden. `id-catalog-v2` añade solo esa clase F4 y no modifica ni rederiva IDs v1; toma los primeros 128 bits, Base32 RFC 4648 mayúscula sin relleno y detecta colisiones como DTC-039. | Cerrar todos los consumidores RNG F4 y separar identidad determinista de aleatoriedad | Las claves v1/v2 no se publican ni implementan. F4 publica v3 y sus vectores cubren extremos del porcentaje, error sí/no, sustitución y todas las decisiones literales. Cambiar cualquier clave, literal, algoritmo o catálogo de `TechnicalPairId` crea una nueva versión y vectores; no hay consumo compartido ni reintento fuera de `insert_length`. |
| <a id="dtc-087"></a>DTC-087 | `id-catalog-v3` es un superconjunto estricto del fichero publicado `catalogs/id-classes-v2.tsv`: conserva literalmente sus filas, vectores, contexto y derivación v1 para ejecución, referencia, individuo, evento, tanda, biblioteca, grupo de lectura y la excepción de plantilla. Añade solo `PreHopInstanceId` y `TechnicalPairId`, ambas con derivación v2, BLAKE3 y codificación canónica de DTC-037/039, 128 bits y Base32 RFC 4648 mayúscula. `PreHopInstanceId` tiene prefijo `PH_`: una instancia propia usa `class=PreHopInstanceId,event_id,library_id,origin=own,template_id`; una contaminante usa `class=PreHopInstanceId,event_id,library_id=recipient_library_id,origin=contamination,source_pre_hop_instance_id`. `TechnicalPairId` tiene prefijo `TP_` y causas `class=TechnicalPairId,event_id,final_library_id,pre_hop_instance_id`. La propuesta homónima `id-catalog-v2` de DTC-086 es antecedente nulo y no se implementa. Colisiones de las nuevas clases fallan. | Hacer derivables, distintos y acotados los IDs de instancia antes/después de contaminación y demultiplexado sin alterar IDs publicados | El `PreHopInstanceId` propio se crea en F4.1; el contaminante se crea al copiar en F4.4; el `TechnicalPairId`, solo tras el salto en F4.5. Configuración efectiva y manifiestos registran `id-catalog-v3` y el mapa `{v1: clases heredadas, v2: PH_/TP_}`. Cada clase tiene vectores con caracteres permitidos, longitudes máximas, copias entre receptoras, reordenación y colisión. |
| <a id="dtc-088"></a>DTC-088 | Las primitivas de F4 son normativas junto a DTC-086. uniform-index-v1(n), con n>0, toma x:u64 de ChaCha12, calcula en aritmética exacta L=2^64-(2^64 mod n), repite mientras x≥L y devuelve x mod n. uniform-open-53-v1 toma k=uniform-index-v1(2^53) y devuelve (k+0.5)/2^53; bernoulli-v1(p) no consume para 0/1 y en otro caso compara uniform-index-v1(2^53)/2^53 < p. Fisher–Yates parcial recorre el prefijo solicitado de la lista canónica y, en el paso i, intercambia con i+uniform-index-v1(n-i). normal_truncated_box_muller_v1 toma dos uniformes abiertos u1,u2, calcula z=sqrt(-2·ln(u1))·cos(2·π·u2) con IEEE-754 f64 y π=0x1.921fb54442d18p+1, forma μ+σz, redondea al entero más cercano con empate alejado de cero y repite el par completo solo hasta caer en los límites inclusivos; no conserva variate. | Evitar que la implementación elija conversiones, rechazo o variante normal al crear vectores F4 | DTC-086 conserva catálogo, ámbitos y consumidores; DTC-088 define sus primitivas y reemplaza para F4 toda dependencia operativa de DTC-081. Los vectores de F4 cubren rechazo, extremos, 0/1, Fisher–Yates y un inserto aceptado/rechazado. |
| <a id="dtc-090"></a>DTC-090 | **SUSTITUIDA POR DTC-093 (27/09/2026)**, que estrena `stutter-v2` y añade la dependencia con la longitud del tracto. Lo que sigue vigente de esta decisión: el predictor, la estructura del paso, la ausencia de reintentos, el catálogo de ámbitos, el encaje en la fase y que la verdad del individuo no cambia. Lo que cambia es **solo la forma de la probabilidad**. Texto original: `stutter-v1` es el modelo normativo del tartamudeo de base de la polimerasa sobre un microsatélite. La probabilidad de resbalón **por ciclo** es `p = clamp(a + b·L, 0, 1)`, donde `L` es la **tirada ininterrumpida más larga** del motivo contenida en la molécula —no el número total de unidades del tracto—, `a` es la fracción calibrada por longitud de motivo y `b` es el coeficiente de la tirada, **que arranca en cero y se publica como no calibrado**. La forma lineal es la de la probabilidad por ciclo; la curva con meseta que describe la literatura es la acumulación de esa probabilidad sobre muchos ciclos, de modo que las dos son ciertas a la vez. La estructura del paso **no se configura**: es exactamente `{−3,−2,−1,+1}` en unidades enteras, y una configuración que declare `+2` o un paso en bases se rechaza; los **pesos relativos** dentro de esa estructura son obligatorios cuando el tartamudeo está activo y no tienen valor por omisión. **Sin reintentos**: los pasos inviables —los que vaciarían la tirada o dejarían la molécula más corta que una lectura— se excluyen **antes** de sortear, y con un solo candidato el sorteo no consume aleatoriedad. El catálogo `f4-technical-v4` extiende `f4-technical-v3` con la etapa `stutter` y los componentes `[event_id, library_id, stutter, template_ordinal, locus_id, decision]`, con `decision` en `bernoulli` y `step`, indexada por **identificador de locus** y nunca por la posición del tracto dentro de la molécula. Las trece decisiones heredadas conservan el ámbito de derivación `f4-technical-v3`, de modo que añadir el tartamudeo **no desplaza ninguna decisión anterior**; una ejecución publica el catálogo del que salen sus decisiones. **DTC-067 no aplica aquí**: obliga a que los **recuentos** aleatorios salgan de `poisson-count-v1`, y aquí no hay recuento de sucesos sobre oportunidades sino un ensayo binario y una elección entre alternativas, cuyos consumidores normativos ya existen —`bernoulli-v1` (DTC-088) y `categorical-f64-v1` (DTC-072)—. El tartamudeo se aplica **después de la quimera y antes de calidad y errores**, una decisión por molécula y tracto sobre la molécula completa, y **no cambia la verdad del individuo**. La molécula cambia de longitud en `paso × longitud del motivo` y no se compensa ni se descarta. Este artefacto **no es una preparación amplificada**: no introduce duplicados y el perfil se sigue registrando como sin PCR. | Fijar el modelo, su catálogo de ámbitos y su encaje en una fase cerrada sin reescribir ninguna subsemilla publicada | Vectores literales de subsemilla de la etapa nueva; los de las trece heredadas no cambian. Un fixture impuro con el total de unidades igualado distingue el predictor de la tirada más larga del número total. Un generador instrumentado confirma que con un solo paso viable el sorteo categórico no consume. Una ejecución sin loci declarados publica exactamente los mismos bytes que antes de que el tartamudeo existiera. |
| <a id="dtc-091"></a>DTC-091 | `pcr-amplification-v1` es el modelo normativo de la preparación de biblioteca **amplificada por PCR**. **No se simula la población de moléculas**, que crece como `2^ciclos`: se muestrea **hacia atrás** desde las parejas que se van a publicar, abriendo **familias de amplificación** consecutivas por ordinal de plantilla hasta cubrir las `N` parejas que la cobertura ya fija, y cada familia cuesta un número de operaciones acotado por constantes de contrato. **El usuario declara la fracción de duplicados esperada, no el parámetro del modelo.** El parámetro interno `m` —parejas por molécula— se **deriva por inversión numérica** de `f(m) = 1 − (1 − e^(−m))/m`, continua y estrictamente creciente, por **bisección sobre `[0, M_MAX]` con exactamente 64 iteraciones fijas**, sin criterio de parada por tolerancia. La exponencial procede **exclusivamente de `libm = =0.2.16`** conforme a DTC-067 y DTC-072: **las dos mitades hacen falta**, porque un número fijo de iteraciones garantiza el mismo **coste** y no el mismo **resultado** —la bisección compara `f(mid)` con el objetivo, y si el último bit de la exponencial varía entre plataformas la comparación cae del otro lado—. `f = 0` devuelve `m = 0` **por retorno anticipado, sin iterar**: sin él la bisección deja un residuo del orden de `10⁻¹⁸` que rompe la identidad byte a byte del caso sin amplificar. Una `f` por encima de la que `M_MAX` puede producir se **rechaza**, nombrando el valor y el máximo alcanzable; no se pinza en silencio. La reproducibilidad se ancla en la `m` **publicada**, no en volver a derivarla. La distribución de copias tras los ciclos explícitos se calcula por **programación dinámica exacta** con el triángulo de Pascal ponderado, sin factoriales ni `libm`, y se sortea **sesgada por tamaño** (`p*(z) ∝ z·p(z)`), porque una familia observada es más probable cuantas más copias tenga su molécula. Se construye **una vez por clase de eficiencia** al resolver el plan, no una vez por familia. El tamaño de familia es una **Poisson de media `λ` condicionada a `K ≥ 1` y pinzada** en el máximo declarado, con el último peso recogiendo **toda** la masa de la cola, normalizada con la suma compensada de DTC-049 y sorteada con **un solo** `categorical-f64-v1`. **No se usa `poisson-count-v1`** porque su rechazo no acota el número de sorteos y su cola no está pinzada; que conste para que nadie lea DTC-067 como una prohibición. La **composición de bases actúa al ELEGIR la molécula**, con `G = ((1 + e·g(GC)) / (1 + e))^c` sobre la fracción de GC del **fragmento completo** contada sobre A/C/G/T, con rechazo y **tope de intentos**; agotado el tope se **acepta el último intento y se cuenta**. El factor entra **además** en la media del tamaño de familia. Con curva neutra `G` vale **exactamente** `1.0` y el ensayo no consume aleatoriedad. El catálogo `f4-technical-v5` extiende los anteriores con la etapa `amplification` y los componentes `[event_id, library_id, amplification, family_ordinal, decision]`, con `decision` en `efficiency_class`, `early_copies`, `gc_accept` —que añade `attempt`— y `family_size`; y con las tres claves de `fragment` que añaden `attempt` para los intentos `a ≥ 1`. El **intento cero conserva literalmente las claves heredadas**, y eso es lo que hace posible la identidad byte a byte. Las trece decisiones heredadas siguen en `f4-technical-v3` y las dos del tartamudeo en `f4-technical-v4`: **ninguna decisión publicada se mueve**. Una ejecución publica el catálogo más alto cuya etapa abre. Las **constantes de contrato** —30 ciclos, 6 ciclos explícitos, 256 de tamaño máximo de familia, 16 intentos, 4 clases de eficiencia y 64 parejas por molécula— viven en **una sola constante de Rust cada una** y un juez las compara con sus copias en el recurso del perfil, en los tres esquemas y en la interfaz. La **cota de recursos se publica, no se impone**: la fase de secuenciación nunca ha tenido control de recursos. | Fijar el modelo, su catálogo de ámbitos, su método de inversión y sus constantes sin reescribir ninguna subsemilla publicada, y dejar por escrito los dos nombres normativos que el programa publica | Vectores literales de subsemilla de la etapa nueva; los de las heredadas no cambian. Una búsqueda negativa sobre el módulo confirma que la exponencial no sale de la matemática estándar, que es lo único que distingue las dos implementaciones en una sola plataforma. Los pesos de la Poisson pinzada se comparan con el cálculo a mano y con el valor que daría la variante sin pinzar. La distribución de los primeros ciclos se compara con una enumeración exhaustiva de caminos escrita aparte. Una huella **congelada** del FASTQ de un escenario sin amplificación ancla la reproducibilidad: comparar dos ejecuciones del mismo binario es una simetría, no una identidad. |
| <a id="dtc-092"></a>DTC-092 | **SUSTITUIDA POR DTC-093 (27/09/2026)**, que estrena `clonal-stutter-v2`. **La derivación de esta decisión no cambia**: lo único que cambia es que sus dos columnas de entrada llegan ya moduladas por la longitud del tracto, de modo que la probabilidad por ciclo deja de ser función de la longitud de motivo sola. Todo lo demás —la genealogía, las guardas de viabilidad, el orden de aplicación, el catálogo y la marca de tosquedad— sigue vigente tal cual. Texto original: `clonal-stutter-v1` es el modelo normativo del **tartamudeo clonal heredable por linaje**. Con preparación amplificada el deslizamiento puede ocurrir en **cualquier ciclo** de la PCR, y una copia que resbala transmite su número de unidades a **toda su descendencia**: un resbalón en el primer ciclo alcanza a media familia y uno en el último, a una sola copia. De ahí sale el fenómeno que el modelo existe para reproducir: **una familia de duplicados puede mostrar dos o tres alelos en proporciones que imitan a un heterocigoto**, procediendo de una sola molécula. **La genealogía es mínima y se declara como APROXIMACIÓN**, no como el coalescente exacto del proceso de ramificación, del que no se afirma tener forma cerrada: la copia `0` de cada familia tiene linaje propio completo, de un paso por ciclo; la copia `i ≥ 1` elige un **progenitor uniforme** entre las anteriores con `uniform-index-v1(i)` —con `i = 1` hay un solo candidato y no consume—, y un **ciclo de divergencia** `d ∈ 0..=c` con `P(d = j) = (1−q)·q^j` para `j < c` y **toda** la masa restante en `d = c`, con `q = 1/(1+e)`; **hereda los resbalones del progenitor ANTERIORES a la divergencia y sortea los posteriores**. **La probabilidad de resbalón por ciclo se DESPEJA, no se suma.** Las dos columnas de la calibración son **totales, no componentes**: la de sin PCR ya contiene el mosaicismo somático, el deslizamiento de la amplificación en puente y el de la síntesis, y la de con PCR contiene todo eso **más** la PCR de biblioteca, de modo que sumarlas contaría dos veces el mismo fenómeno. Se despeja `r_pcr = 1 − (1−r_con)/(1−r_sin)` y se reparte entre los ciclos con `p = 1 − (1 − r_pcr)^(1/c)`, con la potencia **exclusivamente de `libm`** conforme a DTC-067. Si la fracción con PCR no supera a la de sin PCR, `p = 0` y no se consume aleatoriedad. La composición devuelve `1 − (1−r_sin)(1−r_pcr) = r_con`, exactamente la fracción medida, y de ahí sale el orden forzoso: **una preparación amplificada entregada sin el tartamudeo de base subestima lo medido entre un 14 % y un 28 % en motivos de uno a cinco pares de base**, por lo que la configuración **exige el tartamudeo de base activo** cuando hay loci declarados. **Las tasas por ciclo derivadas son TOSCAS y se publican marcadas como tales.** Las dos columnas proceden de **dos individuos distintos**, de proyectos distintos y a coberturas distintas, de modo que la resta puede llevar dentro efectos de individuo y de lote además del de la amplificación. El manifiesto y la configuración efectiva publican **las dos fracciones de entrada, el número de ciclos, la fórmula y la marca de tosquedad con su motivo**. El autor del escenario puede **sustituir** la probabilidad por ciclo por longitud de motivo, y su sustitución prevalece y se publica **marcada como suya**. **Las guardas de viabilidad se aplican al PASO NETO ACUMULADO**, no a cada paso por separado: un resbalón entra en el sorteo solo si deja al menos una unidad completa en la tirada tal como está tras los resbalones anteriores de esa misma copia, y si deja la molécula al menos tan larga como una lectura. Se excluye **antes** de sortear y **no se reintenta**, con la misma regla —una sola función— que el tartamudeo de base de DTC-090. **Orden de aplicación:** sobre la molécula ya quimérica si lo es, **primero** el paso neto heredado del linaje y **después** la decisión de base por pareja. Los tractos se procesan en orden canónico de inicio sobre la molécula, cada paso se aplica en la **tirada ininterrumpida más larga**, y esa tirada **se vuelve a leer entre un paso y el siguiente**, de modo que el segundo puede caer en otra. El paso neto se sortea sobre la molécula **original de la familia** y se **reverifica** contra la molécula **que se va a emitir**, que con una quimera puede no ser la misma; si no cabe, no se escribe. El catálogo **`f4-technical-v6`** extiende los anteriores con las cuatro decisiones del linaje, con los componentes literales `[event_id, library_id, amplification, family_ordinal, decision, copy_index]` para `lineage_parent` y `lineage_divergence`, y `[event_id, library_id, amplification, family_ordinal, decision, copy_index, cycle, locus_id]` para `lineage_slip` y `lineage_step`. **Se indexa por identificador de locus, nunca por la posición del tracto dentro de la molécula.** Nace catálogo nuevo por la misma razón que el v5 no cupo en el v4: el nombre **es el ámbito de derivación**, y una biblioteca amplificada **sin** loci declarados no abre ninguna decisión de linaje y sigue publicando el v5 **byte a byte**. Las trece decisiones heredadas siguen en el v3, las dos del tartamudeo de base en el v4 y las cinco de la amplificación en el v5: **ninguna decisión publicada se mueve**. **La verdad del individuo no cambia**: el tartamudeo clonal vive solo en las moléculas que se secuencian. | Fijar el modelo, su genealogía declarada como aproximación, su derivación y su catálogo sin reescribir ninguna subsemilla publicada, y dejar por escrito los dos nombres normativos que el programa publica | Vectores literales de subsemilla de las cuatro claves nuevas; los de las heredadas no cambian. Un fixture **saturado** con la divergencia forzada al último ciclo exige que todas las copias de una familia muestren el mismo número de unidades, y con la divergencia forzada a cero que sean independientes; las dos variantes juntas cazan además una herencia que copiara los resbalones posteriores a la divergencia. Las seis probabilidades por ciclo se comparan con el cálculo hecho **a mano** desde las dos columnas, y un control aparte comprueba que usar la columna con PCR como si fuera el incremento sobrepasa lo medido. Una prueba de sistema con **6000 moléculas** y tamaño muestral calculado —α = 0,05 a dos colas, potencia del 90 %— comprueba que la fracción que un oráculo cuenta **sobre los FASTQ** reproduce la medida con PCR, 0,17 y 0,038, y la de sin PCR, 0,027 y 0,011; **es tautológica respecto del modelo y discriminante respecto de la implementación**, y se presenta así. La prueba que separa el orden usa un tracto **corto** con los pesos cargados a la contracción, que es el régimen en el que la guarda muerde y las dos partes dejan de conmutar. Un validador independiente rederiva **las familias enteras** desde la semilla —copias tras los primeros ciclos, tamaño de familia, progenitor y divergencia— sin llamar al productor, y ata cada predicción a **su** lectura por el identificador de pareja técnica, que también rederiva. Una huella congelada de un escenario sin amplificación sigue anclando la reproducibilidad. |
| <a id="dtc-093"></a>DTC-093 | `stutter-v2` sustituye a `stutter-v1` y añade la **dependencia con la longitud del tracto**, que es el efecto más fuerte que hay sobre el tartamudeo y el que el modelo tenía apagado. La probabilidad de resbalón **por ciclo** pasa a ser `p = clamp(a·m(T) + b·L, 0, 1)`, donde `a` es la fracción calibrada de la longitud de motivo, `L` es la **tirada ininterrumpida más larga** del motivo contenida en la molécula y `b` su coeficiente —que **sigue arrancando en cero y publicándose como no calibrado**, sin cambio alguno respecto de DTC-090—, y `m` es el **factor de modulación** de la longitud del tracto, que **toma esa MISMA `L`, la tirada ininterrumpida más larga, y NO el número total de unidades del tracto**. No es un detalle: para un tracto puro las dos cifras coinciden, pero para uno **impuro** no, y sin decirlo el documento admitiría dos lecturas distintas del modelo. Se toma la tirada más larga por dos razones. **La primera es de la fuente**: Fungtammasan y col. declaran en sus métodos que «solo se consideraron STR ininterrumpidos (perfectos)», de modo que el eje de su figura 1A es la longitud de una tirada **sin interrupciones**, y aplicar esa curva al total de unidades de un tracto impuro le atribuiría una medida que no hizo. **La segunda es de coherencia**: es el predictor que DTC-090 ya había fijado para la parte lineal, de modo que el modelo no tiene dos predictores distintos para el mismo fenómeno. **La forma anterior no servía y conviene decir por qué**: `a + b·L` supone crecimiento en línea recta, y de 5 a 15 pares de base la fracción medida sube unas **19 000 veces** mientras que una recta daría **3**. Por eso `m` es una **tabla interpolada en escala logarítmica** y no un coeficiente. `m` procede de Fungtammasan y col. (2015), figura 1A, CC BY-NC, digitalizada en §9 de la calibración con una incertidumbre de unas ±0,15 unidades logarítmicas. **`m` MODULA, NO SUMA, y está NORMALIZADA**: las tasas por longitud de motivo de Gymrek ya son promedios sobre los tractos de aquella muestra y **ya llevan dentro** el efecto de la longitud, de modo que multiplicarlas por la curva sin normalizar contaría dos veces el mismo fenómeno —el error contra el que avisa DTC-092—. La **longitud de tracto de referencia**, donde `m` vale exactamente `1.0` y la tasa medida sale intacta, es el **tope del tramo medido**: **15 pares de base en mononucleótidos y 34 en dinucleótidos**. No sale de ninguna fuente porque **ninguna la publica**: Gymrek no da ninguna longitud de tracto, la distribución de longitudes de Fungtammasan es de otra población de loci cuyo agregado propio es 706 veces menor, e invertir la curva sobre las tasas de Gymrek da cuatro longitudes incompatibles y **ninguna** para los homopolímeros con PCR. El responsable la fijó el 27/09/2026 por tres razones: es la **única que no satura** —cualquier ancla por debajo de ~11,7 pb llevaría la probabilidad de los homopolímeros con PCR a 1 dentro del tramo medido—, es el punto que §9 ya usó para su contraste, y deja la modulación actuando **solo hacia abajo**. **Su defecto se publica igual de claro: un tracto de 30 unidades no tartamudea más que uno de 15.** **FUERA DEL TRAMO MEDIDO NO SE EXTRAPOLA**: por debajo y por encima se mantiene el valor del extremo. **SOLO HAY CURVA MEDIDA PARA MONONUCLEÓTIDOS Y DINUCLEÓTIDOS**: los motivos de tres a seis pares de base **no tienen dependencia ninguna**, y eso se declara con todas las letras en la documentación y en el documento publicado en vez de extrapolar desde las otras dos curvas. **La MISMA forma se aplica a las DOS columnas**, con PCR y sin ella, tomada como media logarítmica de ambas, de modo que el cociente entre columnas **no depende de la longitud**: es de ese cociente de donde DTC-092 despeja su probabilidad por ciclo, y hacerlo dependiente de la longitud cambiaría un modelo cerrado apoyándose además en la parte de la digitalización que el texto del propio artículo **no** corrobora. Con ello la composición sigue cerrando **exactamente** a cualquier longitud: `1 − (1−r_sin·m)(1−r_pcr) = r_con·m`. En consecuencia `clonal-stutter-v2` sustituye a `clonal-stutter-v1`: **su derivación no cambia**, solo pasa a recibir las dos columnas ya moduladas, y su probabilidad por ciclo deja de ser función de la longitud de motivo sola. **La SUSTITUCIÓN del autor del escenario PREVALECE y NO se modula**, tanto la tasa por locus como la probabilidad por ciclo: quien declara la tasa de un locus la ha elegido para ese locus, con su longitud dentro, y modularla encima contaría dos veces lo mismo. Se publica marcada como suya. **La dependencia viaja marcada como TOSCA con su motivo entero por los tres caminos de siempre**: constante de Rust, guarda que rechaza el documento si falta, y valor obligatorio en los esquemas publicados. El documento publica además **la curva entera**, un factor por número de unidades, para que la modulación pueda rehacerse sin este programa. **Ningún catálogo de aleatoriedad cambia**: no se abre ni se cierra ninguna decisión, no se mueve ninguna subsemilla, y una ejecución sin loci de repetición publica **exactamente los mismos bytes**. **Enmienda del 06/10/2026:** el cociente entre columnas de los homopolímeros es configurable desde DTC-099; la forma, la normalización y la composición no cambian. | Encender el efecto más fuerte del modelo sin contar dos veces lo que las tasas vigentes ya contienen, y sin reescribir ninguna subsemilla publicada | La tabla se rederiva desde los valores publicados en §9 —media logarítmica, interpolación logarítmica y resta del ancla— escritos a mano en la prueba, no importados del módulo. A la longitud de referencia la tasa medida se reproduce **exactamente**, bit a bit. Una prueba de sistema con los dos tractos puestos en sus longitudes de referencia comprueba sobre los FASTQ que las cuatro fracciones medidas —0,17 y 0,038 con PCR, 0,027 y 0,011 sin ella— salen intactas: es la demostración de que no se cuenta dos veces. Un barrido sobre todas las longitudes comprueba que la composición cierra con error relativo por debajo de `1e-9`, que es el redondeo del viaje de ida y vuelta por la potencia de `libm`. Dos tractos del mismo motivo, uno en su ancla y otro en el extremo corto, se separan por más de un orden de magnitud sobre las lecturas publicadas. Las huellas congeladas de un escenario sin loci de repetición **no se mueven**. |

### 3.1. Catálogo normativo de subsemillas F2 v4 y extensión v5

El fichero que publicará F2.2.8 conserva primero, sin cambio y en el mismo
orden, las seis filas de `subseed-scopes-v3.tsv`. Añade después estas filas; el
orden de componentes mostrado es normativo:

| Ámbito | Componentes ordenados | Consumidor |
|---|---|---|
| `matrix_candidate_count` | `scenario_id`, `reference_id` | `poisson-count-v1`; no se usa para `exact_loci` |
| `matrix_locus_selection` | `scenario_id`, `reference_id` | `weighted-without-replacement-v1`, exclusivamente sobre el conjunto SNV estático |
| `matrix_alternate_allele` | `scenario_id`, `reference_id`, `sequence_id`, `position_1based`, `ref` | `categorical-f64-v1`; identidad completa anterior a ALT |
| `matrix_exact_assignment` | `scenario_id`, `reference_id`, `population_id` | `fisher-yates-v1` sobre todos los `variant_id` canónicos |
| `matrix_frequency_component` | `scenario_id`, `reference_id`, `population_id`, `variant_id` | `categorical-f64-v1`, solo modo estocástico |
| `matrix_frequency_value` | `scenario_id`, `reference_id`, `population_id`, `variant_id` | `linear-uniform-v1` o `log-uniform-v1`; `fixed` no lo consume |
| `indel_count` | `scenario_id`, `reference_id`, `individual_id`, `sequence_id`, `variant_class`, `spatial_scope_id` | `poisson-count-v1` |
| `indel_position` | `scenario_id`, `reference_id`, `individual_id`, `sequence_id`, `variant_class`, `spatial_scope_id` | `dynamic-indel-selection-v1`, que reenumera el estado vigente |
| `insertion_content` | `scenario_id`, `reference_id`, `individual_id`, `sequence_id`, `spatial_scope_id`, `event_index` | `categorical-f64-v1` sobre A/C/G/T por base |
| `initial_genotype` | `scenario_id`, `reference_id`, `population_id`, `individual_id`, `sequence_id`, `position_1based`, `ref`, `alternate_alleles` | `hwe-composition-v1` o `f-diploid-genotype-v1`; `alternate_alleles` es la lista ALT ASCII ordenada del locus, por lo que un sorteo multialélico no depende de un ALT representativo |
| `initial_phase` | `scenario_id`, `reference_id`, `population_id`, `individual_id`, `variant_id` | Fila histórica v4, no consumida por F2.3.3 |
| `initial_phase_locus` | `scenario_id`, `reference_id`, `population_id`, `individual_id`, `sequence_id`, `position_1based`, `ref`, `alternate_alleles` | `phase-permutation-v1` sobre REF+ALT repetidos por dosis |

`variant_class` solo admite `insertion` o `deletion`. `spatial_scope_id` es el
ID estable DTC-039 de la regla espacial efectiva, con `global` como literal
reservado cuando no existe región. `event_index` es el ordinal base 0 del evento
aceptado dentro del grupo canónico; la primera inserción usa `0`. La posición
base 1 y el ordinal se codifican como decimal ASCII canónico. La selección de
ALT precede a la derivación DTC-039 de `variant_id`, que sí incorpora ALT.

### 3.2. Consumidores RNG normativos F2

**Corrección aprobada de DTC-037/DTC-071 (21 de agosto de 2026).** Los
catálogos v1--v3 y `derivation_version=1` permanecen inmutables. El catálogo
v4 es mixto y autocontenido: sus seis filas heredadas conservan
`derivation_version=1` y sus once filas F2 usan `derivation_version=2`. Cada
fila queda ligada de forma legible por máquina en
`scope_derivation_versions`; no existe una versión escalar que pueda ocultar
esa mezcla. El material v2 incorpora `scope_name` ASCII, precedido por su
longitud `u32` big-endian, después de la versión y antes del número de
componentes. Esta separación es necesaria porque cuatro pares v4 tienen los
mismos componentes; no altera ninguna subsemilla F1.
El literal reservado `scenario_id="direct"` sigue siendo interno a identidad y
RNG, no un escenario ni un campo TOML/UI.

**Corrección aprobada de F2.3.3 (23 de agosto de 2026).** La fila
`initial_phase` de v4 queda inmutable y sin consumidor: una clave por ALT no
puede describir una asignación multialélica disjunta. V5 es superconjunto
estricto de v4: conserva literalmente sus filas y añade
`initial_phase_locus` con la misma codificación v2. El productor construye en
orden canónico la lista REF seguida de ALT ASCII, repite cada índice por su
dosis, aplica una sola `phase-permutation-v1` y asigna el resultado a las copias
en el orden de fase 0..P-1. Así cada copia recibe exactamente un alelo y las
dosis se conservan sin depender del orden de filas ni de un ALT representativo.
Desde esta corrección, la metadata de reproducibilidad registra catálogo 5;
v4 sigue disponible y verificable para toda fila histórica.

**Evaluación DTC-043 de `libm 0.2.16`.** Se admite con pin exacto para
`poisson-count-v1` y `log-uniform-v1`: licencia MIT, MSRV 1.63 compatible con
Rust 1.97.1, API pública estable y sin dependencias transitivas normales ni
enlaces nativos. Su mantenimiento activo, pruebas aguas arriba y la cobertura
de las funciones C matemáticas necesarias se revisan como suficientes para esta
dependencia acotada. Se activa solo con `default-features = false`; el coste de
distribución es un crate Rust puro, sin paquetes, herramientas ni bibliotecas
nativas adicionales. El rendimiento aceptado es el de implementaciones
deterministas de `exp`/`log`/`lgamma`/`floor`/`sqrt`, preferido a las rutinas de
la plataforma para conservar reproducibilidad. Su `build.rs` solo emite
configuración Cargo; no ejecuta procesos ni red, ni escribe fuentes, genera o
compila código nativo. La auditoría bloqueada y los vectores de frontera
protegen seguridad, corrección y estabilidad de la integración.

**Evaluación DTC-043 de `serde 1.0.228` (F2.4.3).** Se admite con pin exacto,
`default-features = false` y solo `derive,std`, como dependencia directa del
lector independiente de `truth.jsonl`. La usa exclusivamente para deserializar
un modelo *wire* cerrado (`deny_unknown_fields` y discriminadores explícitos);
los tipos de dominio y la proyección no derivan serialización ni delegan su
semántica al crate. `serde_json` activa además `arbitrary_precision`: los
enteros JSON se conservan como léxico bruto y se comprueban exactamente contra
el dominio `u64`, incluidas las grafías integrales `2.0` y `2e0`, sin pasar por
`f64`. Es Rust puro, con licencia dual MIT o Apache-2.0, sin
enlace nativo ni proceso externo; `serde_derive` queda como macro de
compilación fijada por `Cargo.lock`. `cargo tree --locked -i serde` confirma
que en esta resolución solo lo solicita directamente el simulador. La revisión
de seguridad ejecutó `cargo audit 0.22.2`, pero el sandbox no pudo adquirir el
lock de la base local de RustSec, por lo que no se afirma un resultado de
vulnerabilidades; CI debe repetirla con una base accesible. La alternativa de
un parser manual duplicaba el formato y debilitaba el rechazo de duplicados.

- `open-unit-f64-v1` consume un `next_u64`, descarta los 12 bits bajos y devuelve
  `(((x >> 12) as f64) + 0.5) × 2^-52`, siempre en `(0,1)`.
- `uniform-u64-below-v1(n)` rechaza `n=0` antes de consumir; para `n=1` devuelve
  0 sin consumir. Para `n>1` calcula
  `threshold = n.wrapping_neg() % n` en aritmética `u64`, consume hasta
  obtener `x >= threshold` y devuelve `x mod n`.
- `categorical-f64-v1` recibe pesos efectivos finitos no negativos en orden
  canónico, sumados/normalizados por DTC-049. Cero pesos se omiten; cero
  componentes positivos es error previo y uno solo se elige sin consumir. Con
  varios consume un `open-unit-f64-v1`, usa el primer acumulado estrictamente
  mayor que `u × suma` y, solo por redondeo final, el último positivo.
- `fisher-yates-v1` recorre índices desde `n-1` hasta 1 y permuta cada posición
  con `uniform-u64-below-v1(i+1)`. Con menos de dos elementos no consume.
- `weighted-without-replacement-v1` preflighta que haya suficientes candidatos
  de peso positivo. Por cada selección recalcula en orden canónico los pesos de
  candidatos restantes con DTC-049, aplica `categorical-f64-v1` y retira el
  elegido. Solo se aplica a un conjunto estático de loci SNV: cero selecciones
  no consume; no reintenta ni repone. Una estructura optimizada solo es conforme
  si reproduce exactamente esta secuencia.
- `dynamic-indel-selection-v1` recibe el recuento ya sorteado y un estado lógico
  transaccional. Prevalida sin consumir únicamente configuración, longitud,
  clase y ámbito. Para cada `event_index` desde 0 reenumera desde cero las
  oportunidades del estado vigente en orden canónico de coordenadas DTC-056,
  recalcula pesos efectivos con DTC-049, usa `categorical-f64-v1` sobre el mismo
  ChaCha12 local `indel_position` y entrega la candidata. El orquestador obtiene
  el contenido de una inserción mediante el ámbito separado `insertion_content`,
  aplica el evento al estado transaccional y solo entonces solicita el ordinal
  siguiente. El selector dinámico nunca consume el flujo de contenido. No
  conserva una colección inicial ni asume que las candidatas solo disminuyen.
  Recuento cero no enumera ni consume. Si antes de completar el recuento no hay
  peso positivo, falla en ese ordinal, registra `sampled_count`,
  `accepted_before_failure`, `active_before_rollback` y espacio vigente,
  descarta el estado transaccional y registra `published_count=0`; no publica
  prefijo ni vuelve a sortear Poisson o el evento.
- `linear-uniform-v1` exige `0<=min<max<=1`, siempre como `[min,max)`, y que el
  intervalo contenga al menos un `f64` positivo representable; lo contrario
  falla antes de consumir. Consume un `open-unit-f64-v1` y evalúa, sin FMA,
  `min + u × (max-min)`. Si el redondeo alcanza `max`, devuelve el `f64`
  inmediatamente anterior. Si el resultado es exactamente 0, devuelve
  `f64::from_bits(1)`, el menor `f64` positivo, que el preflight garantiza menor
  que `max`; así un intervalo continuo nunca produce los extremos exactos 0 o 1.
- `log-uniform-v1` exige además `min > 0`, calcula con `libm` y sin FMA
  `exp(log(min) + u × (log(max)-log(min)))` y aplica la misma corrección del
  extremo superior; si el redondeo devolviera un valor menor que `min`, devuelve
  `min`. También usa `[min,max)` y nunca devuelve 0 ni 1 cuando 1 es `max`.
- `hwe-composition-v1` realiza P sorteos `categorical-f64-v1` independientes en
  orden de copia y ordena después los índices de alelo para obtener la
  composición sin fase; una distribución fijada no consume.
- `f-diploid-genotype-v1` calcula las probabilidades ya aprobadas por
  REQ-EVO-47/48, enumera genotipos no ordenados por índices de alelo y aplica un
  único `categorical-f64-v1`; una distribución degenerada no consume.
- `phase-permutation-v1` parte de la composición ordenada. Una ploidía menor que
  2 o todos los alelos iguales no consume; en otro caso aplica
  `fisher-yates-v1` a las P copias, conservando exactamente el vector obtenido.

Cada vector congela entrada, bytes `next_u64` consumidos, salida, número de
consumos y errores previos. Los oráculos no llaman al productor: enumeran casos
pequeños, fuerzan rechazos de `uniform-u64-below-v1`, comprueban PMF, cero,
subnormales, `next_down(1)`, `max=1` y redondeos a ambos extremos, y comparan
HWE/F con probabilidades calculadas aparte. Para indels, un modelo de cadenas
independiente cubre un evento que crea oportunidades, otro que las elimina,
agotamiento antes de `sampled_count`, orden canónico, consumo continuo del
mismo flujo y descarte íntegro del estado al fallar.

**Evaluación DTC-043 de JSON Schema en F1.6 (14 de agosto de 2026).**
`jsonschema` 0.49.9 queda fijado para validar en runtime el esquema general
Draft 2020-12 publicado. Es MIT, declara MSRV 1.85.0 —compatible con Rust
1.97.1— y mantiene actividad y pruebas de conformidad. Se activa con
`default-features = false`: no incorpora resolución HTTP/fichero, TLS, macros ni
asincronía. El adaptador privado `manifest_schema` registra exclusivamente los
dos esquemas embebidos y compila una vez el validador; sus tipos y errores no
atraviesan el contrato del simulador. El crate y sus subcrates propios declaran
`build = false`; el grafo no compila ni enlaza código nativo. Los `build.rs`
transitivos observados son las sondas acotadas de `ahash`/`version_check` y
`num-traits`/`autocfg`: no ejecutan herramientas externas ni red, no modifican
fuentes y no generan código nativo. `cargo audit --no-fetch` no encuentra
avisos en las 149 dependencias locked. El coste aceptado es 84 paquetes locked
adicionales; se evita el grafo HTTP/TLS y se cachea la compilación del esquema.
No existía una línea base de rendimiento específica de validación de manifiesto
y no se reconstruye retrospectivamente. Un corpus diferencial completo/fallido
contrasta controles, saltos y anclajes con Python/jsonschema; el esquema
publicado hace explícito el fin exacto para que una nueva línea no se acepte por
la semántica particular de `$` de un motor. `boon` 0.6.1 se evaluó y descartó:
no cerraba por sí solo esa diferencia y no forma parte de `Cargo.lock`.

**DTC-055: registro histórico cerrado.** La propuesta de panel externo y
recombinación fundadora fue retirada íntegramente el 31 de agosto de 2026. Se
conserva únicamente esta constancia de sustitución; no se conservan aquí
formatos, algoritmos ni ejemplos ejecutables. DTC-089 es la única semántica F5
activa. Una posible importación de panel requeriría una especificación futura
nueva y la meiosis se diseña exclusivamente en F7.1.

**Concreción aprobada de F1.5 para DTC-037.** El catálogo normativo de ámbitos
se versiona separadamente de la codificación. El catálogo 2 añade el ámbito
`reference_generation`, con el único componente `reference_id`, para generar la
referencia inicial. La codificación y la versión de derivación permanecen en 1
y los vectores existentes del catálogo 1 se conservan literalmente; el catálogo
2 añade su vector independiente y su propia huella SHA-256.

**Concreción aprobada de F1.5 para DTC-041/DTC-047.** El FASTA generado
canónico v1 usa 80 bases por línea, mayúsculas, LF y orden de catálogo. Una
ejecución normal de F1 produce y publica cero índices. La consulta indexada de
aceptación es el único consumidor F1 autorizado para crear un FAI bajo demanda;
lo hace fuera del resultado publicado y lo contrasta con `samtools faidx`. Ni el
FAI ni `samtools` pasan a ser salidas normales o dependencias del runtime.

## 5. Decisiones arquitectónicas propuestas

Estas entradas describen el borrador
[arquitectura-0.1.md](./arquitectura-0.1.md). Su estado no implica que exista
implementación.

| ID | Propuesta | Requisitos que satisface | Consecuencias y revisión necesaria |
|---|---|---|---|
| DTP-001 | Arquitectura orientada al dominio con adaptadores para configuración, formatos y procesos externos | REQ-CONF, REQ-REF, REQ-DAT, REQ-VER | Aísla biología de crates y herramientas; aprobar fronteras antes de F1 |
| DTP-002 | Configuración declarada y configuración efectiva son objetos distintos e inmutables en su ámbito | REQ-CONF-18, REQ-CONF-29 a REQ-CONF-31, REQ-CONF-106 y REQ-CONF-116 a REQ-CONF-123 | Permite auditar reglas base, dominancia y normalizaciones; DTC-034 confirma la separación y DTC-035 concreta `effective-config.json` sin mutar retroactivamente la entrada |
| DTP-003 | Identidades tipadas por categoría y registro central de unicidad | REQ-IND-25/26, REQ-EVO-53/54, REQ-RG-04 y REQ-MAN-16 | Impide confundir individuo, evento, conjunto o variante; DTC-039 concreta sintaxis, generación automática y etiquetas |
| DTP-004 | Coordenadas etiquetadas con referencia, secuencia y convenio, convertidas por un único servicio | REQ-IND-38 a REQ-IND-48 y REQ-REC-14 | Hace explícita toda conversión; falta decidir la representación canónica interna |
| DTP-005 | Verdad genómica basada en segmentos trazables y eventos ordenados, materializable como secuencia | REQ-IND-28 a REQ-IND-52 y REQ-REC-11 a REQ-REC-29 | Permite fase y linaje sin exigir copias completas; debe medirse antes de fijar la estructura física |
| DTP-006 | Flujos aleatorios con ámbitos nombrados para biología, muestreo y eventos técnicos | REQ-EVO-12/62 y REQ-CONF-27 | Evita que una observación cambie la evolución; DTC-037 concreta ChaCha12 y la derivación BLAKE3 de subsemillas |
| DTP-007 | Resultados provisionales que solo pasan a completos después de sus validaciones | REQ-SAL-05, REQ-VAL-11 y REQ-MAN-73 | DTC-045 concreta publicación lógica atómica, residuos y ausencia de reanudación v1; la durabilidad física y una reanudación futura exigen evidencia propia |
| DTP-008 | Reproducción, recombinación y selección declarativa compartidas por evolución y familias | REQ-CRU-06/17/31/32 y REQ-EVO-27/29 | Evita rutas biológicas divergentes; exige contratos comunes antes de F7/F8 |
| DTP-009 | Recursos comunes de referencia y ámbitos separados por evento y conjunto | REQ-MAN-72 a REQ-MAN-75 y REQ-SAL-17 a REQ-SAL-19 | Evita duplicación y mezcla; nombres y disposición física siguen pendientes |
| DTP-010 | Validación como consumidor de solo lectura, sin corregir verdad ni evidencia | REQ-VAL-01 a REQ-VAL-13 | Los invariantes deterministas pueden fallar; las magnitudes aleatorias se informan sin tolerancia de aceptación ni reparación posterior |

## 6. Decisiones que bloqueaban F1, ya resueltas

### DTPD-001. Formato y esquema de configuración (resuelta)

Resuelta el 8 de agosto de 2026 mediante [DTC-034](#dtc-034). Se evaluó separar
el formato de autoría de las representaciones de salida frente a reutilizar un
único formato de intercambio. V1 adopta TOML únicamente para la configuración
escrita por la persona usuaria, conserva el original exacto y separa
`DeclaredConfig` de `EffectiveConfig`. La serialización de la configuración
efectiva queda concretada por [DTC-035](#dtc-035); la compatibilidad de v1, por
[DTC-036](#dtc-036).

### DTPD-002. Esquema y formato de manifiestos (resuelta)

Resuelta el 8 de agosto de 2026 mediante [DTC-035](#dtc-035). Se evaluó
reutilizar TOML, usar un formato de intercambio distinto y fragmentar los
manifiestos. V1 adopta JSON conforme a RFC 8259 con JSON Schema Draft 2020-12,
un manifiesto general y otro por evento. Conserva directamente el contenido
obligatorio y enlaza los artefactos independientes, sin JSONL ni tablas
auxiliares. La compatibilidad de v1 queda resuelta por
[DTC-036](#dtc-036); las sumas quedan resueltas por
[DTC-040](#dtc-040) y la publicación atómica por [DTC-045](#dtc-045).

### DTPD-003. Versionado y compatibilidad de esquemas (resuelta)

Resuelta el 8 de agosto de 2026 mediante [DTC-036](#dtc-036). Se evaluó aceptar
solo la versión exacta, leer varias versiones con migración o mantener
conversores explícitos. V1 usa la cadena `"1"`, solo admite versiones conocidas,
rechaza las demás sin interpretación parcial y no migra ni modifica documentos
de entrada. Un conversor explícito solo podrá evaluarse cuando exista una
segunda versión real.

### DTPD-004. Generador aleatorio y estrategia de subsemillas (resuelta)

Resuelta el 8 de agosto de 2026 mediante [DTC-037](#dtc-037). Se descartó el
flujo secuencial único y se eligieron generadores ChaCha12 locales derivados por
ámbitos nombrados, con subsemillas BLAKE3 de 256 bits y separación de dominios.
La raíz y las subsemillas usan 64 hexadecimales minúsculos; la codificación
binaria, el contexto de dominio y el catálogo normativo de ámbitos quedan
fijados en DTC-037, junto con la obtención y registro de una raíz ausente.
La garantía por niveles queda resuelta en [DTC-038](#dtc-038), la política de
dependencias en [DTC-043](#dtc-043) y el toolchain y lock en
[DTC-048](#dtc-048).

### DTPD-005. Garantía exacta de reproducibilidad (resuelta)

Resuelta el 8 de agosto de 2026 mediante [DTC-038](#dtc-038). Se eligió declarar
y comprobar por separado biología, contenido técnico, equivalencia semántica e
identidad binaria. V1 garantiza los dos primeros solo bajo las condiciones
confirmadas de versión, plataforma admitida, entradas, configuración, perfiles,
semillas y dependencias; no promete compatibilidad entre versiones o plataformas
no verificadas ni identidad binaria por defecto.

### DTPD-006. Identificadores (resuelta)

Resuelta el 8 de agosto de 2026 mediante [DTC-039](#dtc-039). Se eligió separar
un identificador estable y tipado de una etiqueta humana opcional. Los IDs
aportados usan el subconjunto ASCII confirmado y los automáticos se derivan de
forma determinista con el contexto BLAKE3, material binario, prefijo de clase y
Base32 RFC 4648 sin relleno concretados en DTC-039. Los ámbitos siguen las reglas
específicas de cada entidad, los catálogos y vectores son normativos para F1 y
toda colisión hace fallar.

### DTPD-007. Identidad de entradas y sumas de comprobación (resuelta)

Resuelta el 8 de agosto de 2026 mediante [DTC-040](#dtc-040). Se eligió SHA-256
sobre los bytes exactos de cada fichero, con representación y nombre de
algoritmo únicos para v1. La verificación cubre copia, cierre y exportación; no
se sustituye por metadatos ni muestreo parcial. MD5 queda limitado a exigencias
de estándares genómicos y BLAKE3 conserva exclusivamente sus dominios de
derivación. La futura dependencia SHA-256 debe seguir
[DTC-043](#dtc-043) y quedar fijada por el lock de [DTC-048](#dtc-048). La
publicación queda resuelta por [DTC-045](#dtc-045).

### DTPD-008. Lectura, escritura y validación de FASTA (resuelta)

Resuelta el 8 de agosto de 2026 mediante [DTC-041](#dtc-041). Se eligió FASTA
de texto sin comprimir y `noodles-fasta` encapsulado tras el contrato propio de
referencia. El dominio conserva la validación de cabeceras, contenido y
alfabeto; la copia suministrada permanece idéntica, mientras que los FASTA
generados usan una escritura canónica. Los índices se crean bajo demanda y
`samtools faidx` solo interviene como validador independiente de pruebas. La
evaluación final acepta `noodles-fasta` 0.66.0 con características mínimas
síncronas y la fija en `Cargo.lock` después de auditar su grafo y sus cuatro
scripts transitivos. La política queda en [DTC-043](#dtc-043) y el
toolchain y lock, en [DTC-048](#dtc-048).

### DTPD-009. Plataforma y distribución iniciales (resuelta)

Resuelta el 8 de agosto de 2026 mediante [DTC-042](#dtc-042). Ubuntu 24.04 LTS
x86-64 GNU es la plataforma nativa oficialmente validada y Ubuntu 24.04 sobre
WSL2, usando el sistema de ficheros Linux, es el entorno de desarrollo admitido.
Las entregas usan CPU portable y se validan en Ubuntu nativo y un entorno OCI
fijado por digest antes de publicar. V1 distribuye un `.tar.gz` con SHA-256 de
sus bytes y una imagen OCI identificada por los digests de sus descriptores; un
archivo OCI exportado tiene además su propio SHA-256 independiente.
WSL2 no define otro comportamiento oficial; macOS, Windows nativo y ARM64 no
tienen soporte oficial en v1. Dependencias, atomicidad y toolchain conservan sus
decisiones separadas [DTC-043](#dtc-043), [DTC-045](#dtc-045) y
[DTC-048](#dtc-048).

### DTPD-010. Política de crates y dependencias nativas (resuelta)

Resuelta el 8 de agosto de 2026 mediante [DTC-043](#dtc-043). V1 prefiere Rust
puro, evita enlaces nativos C/C++ sin necesidad y aprobación expresas y ejecuta
las herramientas genómicas completas como procesos externos encapsulados. Toda
dependencia supera una evaluación explícita, queda detrás de contratos propios
y se limita a las características necesarias. `noodles-fasta` conserva su
aprobación previa sin convertir esta política en un catálogo anticipado. La
validación independiente, las versiones/toolchain y la ejecución concreta de
minimap2 quedan separadas en [DTC-047](#dtc-047), [DTC-048](#dtc-048) y
DTPD-105.

### DTPD-011. Modelo de errores y diagnósticos (resuelta)

Resuelta el 8 de agosto de 2026 mediante [DTC-044](#dtc-044). Los componentes
mantienen errores tipados y causales, mientras una capa pública estable produce
diagnósticos con código `GPS-*`, categoría, mensaje accionable y contexto. El
mismo diagnóstico se presenta como texto humano o JSON. Los fallos esperables
retornan sin `panic` y los defectos internos se distinguen expresamente. DTC-044
fija un documento JSON único por invocación, el catálogo versionado y las
salidas 0, 1, 2, 3, 4, 5 y 70. Los crates de errores siguen DTC-043 y quedan
fijados por el lock de DTC-048 cuando se seleccionen en F1.

### DTPD-012. Publicación atómica, reanudación y residuos (resuelta)

Resuelta el 8 de agosto de 2026 mediante [DTC-045](#dtc-045). V1 distingue el
rechazo previo sin destino, el fallo controlado que publica atómicamente una
ejecución final `fallido` con manifiesto coherente y unidades independientes
válidas, y la interrupción no controlable que deja un provisional `INCOMPLETE`
con diagnóstico solo cuando fue posible. Cada artefacto se finaliza antes en su
propio temporal. La publicación usa un adaptador Linux no-clobber equivalente a
`RENAME_NOREPLACE`; no hay check-then-rename ni fallback con carrera. V1 no
reanuda o reutiliza parciales. Una exportación entre sistemas de ficheros crea y
verifica otro provisional en el destino. Durabilidad física y reanudación por
artefacto no forman parte de esta garantía lógica.

### DTPD-013. Escenarios y métricas iniciales de rendimiento (resuelta)

Resuelta el 8 de agosto de 2026 mediante [DTC-046](#dtc-046). F1 mantiene un
caso de aceptación de 100 kbp por cambio, una serie de crecimiento de 1, 10 y
100 Mb y una observación grande de 1 Gb manual o programada. Todos comparten
recorrido, configuración y semilla fijas y registran contexto, repeticiones,
tiempo, CPU fiable, RSS, bytes, caudal y artefactos. WSL2 mide tendencias de
desarrollo; Linux nativo u OCI es la referencia de publicación. Las cargas son
sintéticas y no fijan objetivos de producto o escala biológica; DTPD-115
conserva futuros umbrales y presupuestos.

### DTPD-014. Herramientas de prueba y validación independiente (resuelta)

Resuelta el 8 de agosto de 2026 mediante [DTC-047](#dtc-047). Los validadores
independientes pertenecen a pruebas, CI y aceptación de publicaciones, no al
recorrido normal ni a la instalación de usuario. F1 contrasta FASTA, SHA-256,
documentos JSON, `.tar.gz` e imagen OCI con herramientas distintas del producto;
las fases de VCF y BAM/CRAM añaden sus referencias correspondientes. CI conserva
versión y plataforma y cubre casos válidos, corruptos y de borde. Corrección
biológica y validez formal mantienen oráculos separados. Versiones exactas y
su registro siguen DTC-048; los validadores futuros permanecen en las DTPD de
cada formato.

### DTPD-015. Versión y política del conjunto de herramientas de Rust (resuelta)

Resuelta el 8 de agosto de 2026 mediante [DTC-048](#dtc-048). F1 fijará en
`rust-toolchain.toml` el canal Rust 1.97.1 stable, perfil `minimal`, componentes
`rustfmt` y `clippy` y target `x86_64-unknown-linux-gnu`; la edición 2024 se
declara únicamente en cada `Cargo.toml` o en el workspace. Se versionará
`Cargo.lock` y se usará `--locked` en CI y entregas. No se declara otra MSRV.
Las versiones de crates se seleccionarán tras DTC-043 y quedarán fijadas por el
lock. Toda actualización será aislada y volverá a ejecutar pruebas,
validadores, vectores reproducibles y escenarios de rendimiento antes de
reformular garantías. Esto sustituye la propuesta accidental basada en Rust
1.94.0 local y no exige Rust a quien usa el binario o imagen.

### DTPD-101. Representación interna común de coordenadas (resuelta)

Resuelta el 14 de agosto de 2026 mediante [DTC-056](#dtc-056). El dominio usa
intervalos base 0 semiabiertos y separa posiciones de base y puntos interbase,
siempre ligados a referencia y secuencia. Configuración y mapas propios siguen
base 1 cerrada y los formatos ajenos conservan su convenio oficial. Las
conversiones ocurren solo en fronteras, usan aritmética comprobada y nunca
recortan un valor inválido. F1.6 crea el servicio común y F2 lo amplía sin
reabrir la representación.

### DTPD-102. Estado haplotípico y fichero obligatorio de verdad (resuelta)

Resuelta el 15 de agosto de 2026 mediante [DTC-061](#dtc-061). El estado interno
es híbrido, segmentado, compartible y materializable bajo demanda. Su proyección
externa obligatoria es `truth.jsonl` (v1 al resolverse; v4 y v5 desde el 30 de
septiembre de 2026, según DTC-061), UTF-8 sin comprimir, independiente de
los tipos Rust y ordenada canónicamente. Combina esquemas por tipo de registro
con recuentos/digests de cabecera y un validador semántico/reconstructor
independiente; cualquier relación ausente, duplicada, desordenada o
inconsistente falla. Los genotipos deben coincidir exactamente con las copias
haplotípicas y todo VCF se deriva de `truth.jsonl`.

### DTPD-103. Contrato VCF y validación independiente (resuelta)

Resuelta el 15 de agosto de 2026 mediante [DTC-062](#dtc-062) a
[DTC-066](#dtc-066). El contrato integrado fija VCF 4.5; `.vcf.gz` BGZF con un
único TBI/CSI o texto explícito sin índice; normalización y agrupación
multialélica internas; `GPS_CASE`, `GPS_LOSS`, `GPS_COMPLEMENT` y el complemento
condicional; y aceptación externa con `bcftools`/HTSlib 1.24 más un oráculo
semántico independiente. No se adopta EBI `vcf-validator` como autoridad VCF
4.5.

La puerta registra la versión y ejecuta estas plantillas de argumentos, siempre
sobre copias temporales y sin modificar la salida publicada:

```text
bcftools --version
tabix --version
bcftools view --no-version --output-type v --output <parsed.vcf> <input.vcf[.gz]>
bcftools norm --no-version --check-ref e --fasta-ref <reference.fa> --output-type v --output <normalized.vcf> <input.vcf[.gz]>
bcftools index --stats <input.vcf.gz>
bcftools view --no-version --regions <region> --output-type v --output <region.vcf> <input.vcf.gz>
tabix -l <input.vcf.gz>
tabix <input.vcf.gz> <region>
bcftools index --tbi --output <validation-copy.vcf.gz.tbi> <validation-copy.vcf.gz>
bcftools index --csi --output <validation-copy.vcf.gz.csi> <validation-copy.vcf.gz>
```

Solo se ejecuta la última plantilla que corresponde al tipo decidido para el
fixture; se exige ausencia del otro índice. El flujo normalizado debe ser
equivalente registro a registro al original según el comparador independiente:
`bcftools norm` no reescribe ni repara el artefacto. Las consultas de
`bcftools` y `tabix` se comparan con el intervalo esperado del corpus. El corpus
incluye positivos de texto, BGZF+TBI y BGZF+CSI, además de una mutación aislada
por cada invariante sintáctico y semántico.

### DTPD-116. Selección reproducible de intervalos para referencia reducida (resuelta)

Resuelta el 14 de agosto de 2026 mediante [DTC-057](#dtc-057). El selector
directo proporcional reparte bloques completos por secuencia y
una terminación completa en el espacio conjunto, selecciona regiones completas
y después un único residual y materializa solo tras ordenar por
secuencia/coordenada. El preflight rechaza configuraciones inviables sin
reintentos o resultados parciales. El vector canónico vigente no promete
compatibilidad con el selector transitorio anterior; preflight, vector congelado
y reconstrucción origen+mapa son los oráculos obligatorios.

### DTPD-118. Formatos y generación de espectros/matrices SNV (resuelta)

Resuelta el 12 de agosto de 2026 mediante [DTC-054](#dtc-054). El contrato
runtime es un par TSV+JSON versionado, íntegro y ligado a una referencia exacta;
solo admite `alt_af`. Una preparación explícita convierte la tabla mínima a ese
par sin modificar el original y la generación integrada o independiente produce
los mismos bytes y permite reutilización. Los espectros de frecuencia usan
componentes no solapados `fixed`, `linear_uniform` o `log_uniform`; el modo
exacto aplica mayores restos y el estocástico es predeterminado. Posiciones,
alelos, frecuencias, genotipos y fase poseen ámbitos RNG separados. DTC-070
completa el modo exacto con una permutación independiente de todos los
`variant_id` por población y un flujo distinto para el valor. DTC-071 resuelve
RAF2-028 al fijar nombres, componentes, orden y versión del catálogo RNG F2 sin
circularidad respecto a `ALT`.

La generación espacial es uniforme o usa regiones no solapadas con peso y
`outside_weight`; genera solo loci bialélicos. Las sustituciones usan 12 cambios
dirigidos, con perfil uniforme predeterminado, tabla propia, atajo Ti/Tv o
el ejemplo aproximado Harris 2017 seleccionado expresamente. El ejemplo
alélico global aproximado 1000 Genomes y el ejemplo SNV12 conservan fuente,
receta, licencia, transformación y SHA-256 y nunca se cargan implícitamente. No
existen perfiles incorporados exactos por superpoblación; cada población puede
recibir un perfil suministrado. Las poblaciones automáticas se generan independientemente; una
relación entre ellas debe venir en una matriz suministrada. Los valores 0 y 1
son válidos: un `ALT` globalmente ausente se elimina y cuenta como candidato sin
reemplazo, mientras uno fijado se conserva. VCF→TSV y generación multialélica
quedan en DTF-038/039.

### DTPD-119. Haplotipos sintéticos F5 (resuelta)

Resuelta por [DTC-089](#dtc-089). F5 construye internamente una distribución
finita de haplotipos SNP-only por bloque×subpoblación y asigna fundadores
intactos. No existe panel externo, paquete de preparación, construcción
recombinante, mosaico ni fuente donante. F6 añade SV/variantes complejas y F7.1
introduce la primera meiosis. Los contratos físicos, la búsqueda reproducible y
los diagnósticos pertenecen exclusivamente a DTC-089.

### DTPD-120. Implementación reproducible de Poisson (resuelta)

Resuelta el 15 de agosto de 2026 mediante [DTC-067](#dtc-067). El algoritmo
normativo `poisson-count-v1` usa inversión/Knuth para `0<lambda<10`, PTRS de
Hörmann para `10<=lambda<=2^52` y cero sin consumo RNG. La conversión uniforme,
`libm` 0.2.16, los límites numéricos, el fallo por espacio insuficiente, los
metadatos, vectores y oráculos quedan fijados; ningún consumidor puede cambiar
la media, recortar el resultado o reutilizar Poisson para frecuencias, alelos,
genotipos o fase.

### DTPD-121. Receta de perfiles humanos compactos (resuelta)

Resuelta el 15 de agosto de 2026 y sustituida ese mismo día mediante
[DTC-068](#dtc-068). F2 incorpora exclusivamente dos aproximaciones compactas:
`human_1000g_global_sfs_approx_v1` y
`human_harris2017_snv12_approx_v1`. No descarga ni procesa el VCF Phase 3 y no
usa `bcftools` para construir estos recursos.

La receta SFS toma los recuentos publicados 64:12:8, usa el mínimo positivo
`1/5008` derivado de las 5.008 copias autosómicas observables del panel y crea:

```text
[1/5008, 0.005)  weight=64  sampling=log_uniform
[0.005, 0.05)    weight=12  sampling=log_uniform
[0.05, 1)        weight=8   sampling=linear_uniform
```

El límite 1 completa la cobertura, pero el intervalo es semiabierto y no lo
produce; `alt_af=1` solo puede proceder de `fixed=1`. No hay masa fija en cero ni frecuencias positivas inferiores a `1/5008`; un
perfil propio puede representarlas. El mínimo es parte de esta aproximación y
no una restricción general del simulador.

La receta SNV12 descarga o usa una copia verificada de
`elife-24284-fig1-data1-v2.txt` (3.644 bytes; SHA-256
`970cc71578fc008ef4e0d75c7c460ab874073da7c5b5226c526b387627c97a7f`),
normaliza cada fila AFR/AMR/EUR/EAS/SAS, calcula su media aritmética con peso
`1/5`, agrega los 96 contextos a seis clases colapsadas por hebra y divide por
igual cada peso entre la sustitución canónica y su complemento inverso para
obtener los 12 cambios dirigidos; DTC-049 normaliza al final. El recurso se
etiqueta como aproximación descriptiva, no como tasa universal.

F2.2.6 fija un generador pequeño y un oráculo independiente que verifican
fuentes, licencia CC BY, tamaño/SHA-256 de entrada, literales 64:12:8,
transformaciones, orden y SHA-256 de ambas salidas incorporadas. La
materialización exacta desde VCF y los perfiles empíricos exactos quedan en
DTF-016/post-v1 y no bloquean F2.2.

## 7. Decisiones de planificación posteriores

> **Corrección de numeración F5/F6 (31 de agosto de 2026).** La población
> inicial, incluidos F diploide, perfiles LD y fundadores intactos, pertenece a
> F5 (F5.0–F5.6) bajo DTC-089. Las referencias
> estructurales, SV, complemento VCF y multirreferencia pertenecen a F6
> (F6.1–F6.6). Las menciones históricas previas a `F6.4` como consumidor del
> panel se consideran retiradas; F5.3 materializa subpoblaciones antes de la
> publicación y F5.5 solo audita la inmutabilidad del paquete.

### 7.1. Decisiones de planificación

| ID | Decisión | Opciones | Recomendación razonada | Consecuencias | Bloquea |
|---|---|---|---|---|---|
| DTPD-104 | Adaptadores SAM/BAM/CRAM e índices | Crate; biblioteca nativa; herramientas externas; combinación | Comparar cumplimiento de estándar, CRAM con referencia exacta, ordenación, índices, cabeceras y distribución | Determina E/S, procesos, errores y rendimiento | F6.6 |
| DTPD-105 | Concreción del adaptador de proceso externo para minimap2 y otros mapeadores | Descubrimiento explícito o ruta configurada; política de invocación y sandbox; compatibilidad de versiones; clasificación de códigos, señales y fallos | Partir del proceso externo ya decidido por DTC-043 y concretar descubrimiento, comando reproducible, aislamiento, compatibilidad y diagnóstico sin ocultar versión, parámetros, stdout, stderr ni código de salida | Afecta instalación, seguridad, logs y reproducibilidad; no reabre servicio frente a proceso | F6.6 |
| DTPD-106 (resuelta 2026-08-28) | Perfiles WGS cuantitativos incorporados | Perfil inicial propio `gps-hiseqx-pcrfree-pe150-v1`, recurso `resources/f4/gps-hiseqx-pcrfree-pe150-v1.json`, SHA-256 `3d8878ab920237dae03705276a26dbb0823941830d531374d35aa1c540f6deea`: aproximación sintética explícita, no calibración empírica. R1 usa Q36 (ciclos 1–10), interpolación lineal redondeada a Q32 (11–100) y a Q28 (101–150); R2 usa Q35, Q30 y Q25 en los mismos tramos. Admite longitud común 1–150 por prefijo y rechaza una mayor. El inserto usa `normal_truncated_box_muller_v1`, media 350 pb, desviación 20 pb, límites inclusivos 250–450 pb, uniformes abiertos de 53 bits, redondeo al entero más próximo con empate alejado de cero y sin variate normal en caché. ART `HSXn` es solo referencia comparativa externa, no se incorpora ni redistribuye. Tras cada error, la base alternativa se elige uniformemente entre las otras tres. PCR no aplicable, ópticos/quimeras/GC/*index hopping* a cero y contaminación `0–0 %`; sin adaptadores por no existir inserto menor que lectura. El orden es fragmento, quimera opcional, calidad/error, contaminación e *index hopping*; quimeras e *index hopping* se sortean Bernoulli e independientemente por pareja, el primero antes y el segundo después de contaminación. Las quimeras no cambian el número de parejas y *index hopping* las transfiere uniformemente a otra biblioteca. F4 no genera UDI inesperados, aunque los descarta defensivamente si llegan a la frontera; admite por evento desplazamiento de calidad y escala de error independientes para R1/R2, con valores neutros `0 Phred` y `1`, aplicados antes de recalcular coherentemente la calidad FASTQ final, codificada Sanger `Phred+33`, redondeada al entero más próximo y validada en `0..=93`; no hay selectores por individuo o región y siempre se publican los valores efectivos | Aprobación explícita del usuario en la ronda de diseño F4 del 2026-08-28; perfil literal y algoritmo aprobados como decisión técnica directa el 2026-08-28; evidencia: recurso versionado, especificación, arquitectura, plan y trazabilidad | El ID de perfil se conserva porque su semántica cuantitativa no cambió; solo el enlace al catálogo técnico no publicado pasa a v3. Sustituirlo por datos reales exige un nuevo ID/recurso/digest. DTF-053–056 conservan lo excluido | — |
| DTPD-107 | Calibración empírica (retirada de v1) | Sin elección técnica en v1 | La decisión funcional DTC-023 aplaza la capacidad completa; cualquier diseño del estimador y del artefacto corresponde a DTF-016 | No bloquea ninguna fase de v1 | — |
| DTPD-108 (resuelta 2026-08-28) | Procedencia agregada de lecturas | Registro por lectura; registros por bloques; resumen agregado | Conservar únicamente en manifiesto/resumen la semilla, perfil efectivo, cobertura objetivo/real, contaminación objetivo/real y aportes por donante, recuentos de quimeras, errores por mate con oportunidades/tasa realizada, duplicados y descartes, y hashes/tamaños FASTQ | Aprobación explícita del usuario en la ronda de diseño F4 del 2026-08-28; evidencia: contrato agregado incorporado en especificación, plan, trazabilidad y `git diff --check` superado | No existe fichero ni procedencia individual por lectura; reduce volumen y no cambia el FASTQ ni la verdad biológica | — |
| DTPD-109 (resuelta 2026-09-01) | Fórmulas de F diploide y Ne reproductivo realizado | DTC-089 fija las probabilidades genotípicas multialélicas, `F_realized=1-H_observed/H_expected` y el Ne por contribuciones `Ne=(N_p*k_bar-1)/(k_bar-1+V_k/k_bar)` | Aplicar las fórmulas solo donde sus denominadores y datos están definidos; usar `not_applicable`/`not_evaluable` con causa y no adaptar F a otras ploidías | F5 informa F por bloque y Ne no aplicable; F7.3 calcula Ne sin convertirlo en objetivo ni modificar descendientes | — |
| DTPD-110 (resuelta 2026-09-01) | Distribuciones y tolerancias poblacionales F5 | DTC-089 fija los cuatro catálogos distributivos, sus dominios, reducciones y las tolerancias inclusivas `max(q)`, `0.05` y la controlada efectiva por bloque×subpoblación en `[0,0.5/C]`, predeterminada `0.5/C`; DTC-049 conserva normalización general | Las fases posteriores solo pueden definir distribuciones nuevas mediante su propia decisión, sin reabrir ni reinterpretar DTPD-110 | Configuración inviable se rechaza antes del RNG propio de la etapa y búsqueda agotada se distingue de imposibilidad demostrada | — |
| DTPD-111 | Mapas, posiciones y comparación estructural para recombinación meiótica | Formato de mapa, generación de puntos, interpolación, bivalentes y representación canónica de segmentos por definir | Diseñar en F7.1 un motor meiótico propio que consuma copias fundadoras/evolucionadas y emita mosaicos reconstruibles; validar cobertura, monotonicidad, orientación, puntos de ruptura y equivalencia estructural | Afecta solo F7.1 y posteriores. No reutiliza DTC-055, no presupone panel, no bloquea F5 y no introduce recombinación fundadora | F7.1 |
| DTPD-112 | Representación de grafos y reglas declarativas | Estructuras anidadas; tablas de nodos/aristas; lenguaje declarativo limitado sin código | Elegir la sintaxis más validable que represente Y/O/NO, ámbitos y dependencias sin ejecutar código | Afecta mensajes de ciclo, migración, fundación y familias | F7.5 y F8.1 |
| DTPD-113 | Paralelismo reproducible | Por referencia, individuo, evento, región o tareas; ejecución secuencial de referencia | Posponer la elección hasta perfilar y conservar ámbitos aleatorios y orden lógico independientes del planificador | Puede cambiar uso de memoria y orden de salida; no debe cambiar verdad | F9.3 |
| DTPD-114 | Formato del paquete exportado | Directorio autocontenido; archivo contenedor; formatos archivados comprimidos | Elegir por interoperabilidad, preservación de permisos/rutas e integridad antes/después, con eventos grandes | Afecta portabilidad y almacenamiento duplicado | F9.6 |
| DTPD-115 | Objetivos y presupuestos cuantitativos | Umbrales por escenario; presupuesto configurable; perfiles de recursos | Fijarlos después de medir escenarios representativos; todo presupuesto debe incluir el mínimo, la trazabilidad única y las salidas activadas, y su exceso rechaza sin degradación automática | Define expectativas de escala y criterios de publicación | F9.1/F9.4 |

### 7.2. DTPD-117. Concreción técnica de Marimo (simplificada)

DTPD-117 queda resuelta para F2 mediante DTC-079. Marimo es una interfaz local
auxiliar y no una nueva capa de dominio.

1. `pyproject.toml` y `uv.lock` fijan el entorno reproducible; el entorno
   creado localmente no se versiona. Su creación explícita mediante
   `uv sync --frozen` puede descargar de Internet las dependencias resueltas por
   el lock; después, el arranque y el uso ordinario no acceden a la red.
2. El servidor escucha en `127.0.0.1`; el lanzador muestra la URL, admite un
   puerto no privilegiado y no ejecuta al abrir.
3. F2 edita un único TOML general. Escenarios y conjuntos multidocumento se
   implementarán con su primer consumidor real.
4. `validate` y `run` se invocan directamente con argumentos separados. Rust
   conserva toda autoridad semántica.
5. Una edición invalida la última validación y la ejecución requiere una acción
   explícita y un destino visible.
6. La UI muestra diagnósticos, progreso, manifiestos y resúmenes sin cargar
   artefactos genómicos grandes.
7. F2 prueba instalación desde repositorio, arranque local, controles básicos y
   paridad UI→TOML→Rust. F9 verifica distribución, plataforma oficial,
   accesibilidad de entrega y endurecimiento.

No forman parte de F2 un protocolo de capacidades, hashes de contratos entre
componentes del mismo árbol, copia privada de grafos TOML, almacén persistente de
sesiones, supervisión adversaria de procesos, navegador fijado como artefacto ni
defensas frente a manipulación local activa.
## 8. Decisiones futuras que no bloquean v1

Esta tabla es un resumen derivado. El registro canónico con fuentes, estado,
destino y dependencias es
[registro-ampliaciones-futuras-0.1.md](./registro-ampliaciones-futuras-0.1.md).

Las entradas DTF-001 a DTF-048 se definen una sola vez en ese registro. Este
documento conserva únicamente las decisiones técnicas DTC, DTP y DTPD.

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

DTPD-001 a DTPD-015 están resueltas mediante DTC-034 a DTC-048. Quedan así
documentadas las decisiones base que bloqueaban la planificación ejecutable de
F1. Esto no significa que F1 esté planificada, implementada ni cerrada. Algunas
decisiones conservan ampliaciones aplazadas; por ejemplo, F1 tiene semántica de
completitud, pero no reanudación automática.

DTC-049 resuelve la representación, suma, tolerancia general y normalización de
REQ-CONF-106/107. DTPD-110 está cerrada por DTC-089; las distribuciones nuevas
de fases posteriores exigen decisiones propias y no reabren ese contrato.

Las decisiones de la tabla DTPD que no llevan la marca «resuelta» deben
resolverse antes de ejecutar el hito que señalan. DTPD-106 y DTPD-108 quedan
formalmente resueltas el 2026-08-28. DTPD-117, DTPD-118 y DTPD-119 ya están resueltas
mediante DTC-079, DTC-054 y DTC-089, respectivamente. No bloquean
la redacción del plan de F1 salvo que su evaluación
revele una incompatibilidad con un contrato base de configuración, identidad,
referencia o reproducibilidad.

DTPD-101 y DTPD-116 están resueltas mediante DTC-056/DTC-057. F1.6 implementa
el servicio común de coordenadas y la selección versionada; F2.1/F2.2 consumen
el mismo servicio para segmentos, posiciones y anclas sin decisiones locales.

DTC-050 a DTC-079 confirman la capacidad de reducción, la obligación Marimo,
las dos fuentes poblacionales, las longitudes exactas de *indels* F2 y el
anclaje de matriz a referencia, la verdad `truth.jsonl`, VCF 4.5 como única
versión de salida, el catálogo/algoritmos RNG F2, las capacidades efectivas y la
evidencia de plataforma. Dentro de ese grupo, DTC-054 cierra los formatos
y algoritmos de la matriz SNV. F1 está cerrada en `e14705a` con cierre
documental `efca265`; F2.0 queda definido por la frontera local proporcional de
DTC-079, sin puerta experimental separada; F2.2/F2.3 ya no están bloqueadas por
DTPD-118; F5 ya no está bloqueada por DTPD-119, F6 no hereda su contrato y F7.1
continúa bloqueada por DTPD-111. Las
recomendaciones pendientes no se consideran aprobadas hasta que consten
revisión y fecha.

### Ampliación DTC-039 para F2.2.4

`VariantId` usa el prefijo `VAR_`, es común entre poblaciones y se deriva de
`reference_id`, `sequence_id`/contig, `position_1based`, `ref` y `alt`. No
incorpora fila, población, ruta, orden físico ni semilla. El vector y la
codificación BLAKE3 quedan fijados en `catalogs/variant-id-v1.tsv`.
<a id="dtc-089"></a>

## DTC-089 — F5 sintético reproducible

Esta decisión es la autoridad técnica única de F5 y sustituye por completo
DTC-053 y DTC-055. F5 genera internamente haplotipos de SNP bialélicos A/C/G/T
por bloque y subpoblación y asigna fundadores intactos. Rechaza paneles,
haplotipos suministrados, `alt_af` intrabloque, *indels*, repeticiones, SV,
mosaicos, recombinación, mapas, genealogía y `HaplotypeProcess`. F6 añade
variación compleja y F7.1 define por primera vez la meiosis.

Revisión de admisión aprobada el 4 de septiembre de 2026: el algoritmo v7
reemplazó los intermedios v5/v6 y pasó a ser la única autoridad de recursos y
binding para F5.1.2. La revisión F5.3, del 12 de septiembre de 2026, lo
sustituye por v8, que contabiliza además el repositorio de identidades
fundadoras; la autoridad es siempre la versión vigente, que es v8.

### Escenario TOML canónico

F5 es un conjunto de dos documentos de autoría cerrados, no un TOML aislado.
El TOML general F5 (`config_kind="general"`, `schema_version="3"`) es dueño de
`seed`, `ploidy`, `reference_id`, `reference_path`, `reference_sha256` y de la única tabla
`[f5_scenario]`, que contiene exactamente `scenario_id` y `scenario_path`. Una
ejecución F5 exige esa tabla exactamente una vez: no admite una lista, un enlace
ausente ni varios escenarios. Esta cardinalidad corresponde al único recorrido
F5 y a su único conjunto de salidas directas; una ampliación multiescenario
necesitará otra versión y rutas e identidades de salida propias. La tabla del
general es exactamente:

```toml
[f5_scenario]
scenario_id = "scenario"
scenario_path = "scenarios/scenario.toml"
```

El escenario siguiente solo define la política F5. Se rechazan
campos desconocidos, claves duplicadas y la ausencia de un campo obligatorio:

```toml
config_kind = "f5_scenario"
schema_version = "1"
scenario_id = "scenario"
founder_count = 1

[f5]
profile_metadata_path = "profiles/profile.metadata.json"
profile_metadata_sha256 = "<64 hex minúsculas>"
inbreeding_f = 0.0
max_layout_attempts = 1000
max_haplotype_attempts = 100
max_swap_proposals = 100000
max_haplotype_cells = 16777216
max_operation_units = 4294967296
max_working_memory_bytes = 8589934592
max_output_bytes = 2147483648
coverage_fraction = 0.10

[[subpopulations]]
population_id = "P1"
individual_count = 1
sampling = "probabilistic"
founder_keys = ["f001"]

# Sustitución opcional de una unidad bloque×subpoblación.
#[[sampling_overrides]]
#population_id = "P1"
#block_index = 0
#mode = "controlled"
#controlled_max_abs_error = 0.0001
```

La ruta de `--config` se resuelve una vez contra el directorio de trabajo de la
CLI y se normaliza; el directorio de ese TOML general es su origen.
`reference_path` es obligatoria, relativa a ese origen y no admite `..`; sus
bytes FASTA se validan contra `reference_id` y `reference_sha256` antes del
preflight de recursos.
`f5_scenario.scenario_path` se resuelve contra ese origen, y los recursos del escenario
contra el directorio del escenario ya resuelto; ninguna ruta admite `..`.
`scenario_id` del enlace y del escenario deben ser idénticos. No hay búsqueda
por workspace, variable de entorno ni directorio actual posterior. El general
es la única autoridad TOML para semilla, ploidía y referencia; por tanto el
escenario no puede redeclararlos. La semilla resuelta aparece exactamente en
`effective-config.json#/seed`, `effective-config.json#/f5/seed` y
`manifest.json#/f5/seed`; los tres campos son cadenas con patrón
`^[0-9a-f]{64}$` y sus valores son idénticos. Si el TOML general declara
`seed`, coinciden además con ese literal; si la omite, la raíz CSPRNG se resuelve
antes de F5 y se copia sin transformación en los tres campos. Las demás copias
resueltas en `effective-config.json` y `manifest.json` también son obligatorias
e idénticas: los hashes son las mismas cadenas hexadecimales, `ploidy` el mismo entero, y referencia el mismo par
ID+SHA-256. El hash de perfil es propiedad del escenario y sus copias efectiva
y manifiesto también deben coincidir. F5 crea contratos externos nuevos:
`config-general-v3.schema.json`, `effective-config-v5.schema.json` y
`manifest-general-v5.schema.json`; el escenario usa
`config-f5-scenario-v1.schema.json`. F5.0 publica estos cuatro esquemas completos
junto con `profile-metadata-f5-v1.schema.json`, `truth-v2.schema.json`, los
contratos completos de los siete TSV y `f5-failure-v1.schema.json`; ninguna
etapa posterior añade campos ni extiende esos esquemas. Por tanto, las versiones no se confunden:
general `"3"`, escenario F5 `"1"`, envoltura efectiva `"5"` con objeto F5
`"1"`, y manifiesto general `"5"` con objeto F5 `"1"`. Los contratos general
v2, efectivo v4 y manifiesto general v4 son históricos e inmutables: F5 no los
extiende, reescribe ni acepta como documentos F5. `truth-v1.schema.json` y su
contrato textual también son históricos e inmutables; F5 usa exclusivamente la
nueva v2 descrita en «IDs, RNG y formatos de salida».

**Excepción expresa del 12/10/2026 (decisión del responsable, DT-041).** El
esquema del manifiesto general v4 se amplió **una vez** para admitir, en cada
evento técnico, el perfil con PCR (`gps-hiseqx-pcr-pe150-v1`, con su huella)
con el catálogo `f4-technical-v5`, atados entre sí: sin PCR sigue exigiendo
exactamente el perfil sin PCR con `f4-technical-v3`. El motivo: desde que
existe la preparación amplificada, **toda** ejecución amplificada sin loci de
repetición publica esta versión del manifiesto con ese perfil y ese catálogo, y
el contrato v4 la rechazaba. Ampliarlo **no mueve ningún byte publicado**; la
alternativa —que esas ejecuciones publicaran otra versión— sí los movía. La
excepción se limita a eso: no reabre el manifiesto general v4 para nada más, ni
toca el general v2, el efectivo v4 ni la verdad v1, que siguen inmutables.

`scenario_id`, `population_id` y los elementos de `founder_ids`/`founder_keys`
cumplen DTC-039. La ploidía efectiva del general pertenece a `1..=10`,
`founder_count>0` y es exactamente la suma comprobada de `individual_count>0`.
Cada subpoblación usa exactamente uno de `founder_ids` o `founder_keys`, cuya
longitud es `individual_count`; los IDs son globalmente únicos y las claves,
únicas dentro de la subpoblación. `inbreeding_f` de `[f5]` es el F de población;
cada subpoblación puede sustituirlo. Ambos son `f64` finitos en `[0,1]`, solo
pueden aparecer con ploidía 2 y el global vale `0.0` si se omite. Con otra
ploidía ambos deben estar ausentes y el valor efectivo es `null`/no aplicable;
las copias se muestrean independientemente. `sampling` se declara por
subpoblación y admite exactamente `probabilistic` —predeterminado— o
`controlled`; es el modo solicitado predeterminado para sus bloques. Una fila
`[[sampling_overrides]]` puede sustituirlo en exactamente un
`population_id`×`block_index`. No puede haber dos filas para la misma pareja y
cada índice debe pertenecer a `0..B-1`, donde B es el número de bloques resuelto
antes del layout; se refiere al orden canónico final FASTA/inicio/fin. Cada fila
contiene `mode` y, solo si es `controlled`, puede contener
`controlled_max_abs_error`. El modo probabilístico prohíbe esa tolerancia. En
controlado, la tolerancia efectiva es la declarada o, si se omite, `0.5/C` para
esa subpoblación; debe ser finita y cumplir `0<=t<=0.5/C`. Esta tolerancia se
aplica por haplotipo y no sustituye DTC-049. En diploides, cada unidad controlada
exige F efectivo 0. Las rutas son relativas al escenario, no contienen `..`, y cada ruta tiene
su SHA-256 compañero. `excluded_regions_path` y
`excluded_regions_sha256` son opcionales como pareja inseparable.

La opcionalidad exacta del TOML es esta; ninguna tabla ni campo adicional está
permitido:

| Tabla | Campo | Tipo/dominio | Presencia/valor efectivo |
|---|---|---|---|
| raíz | `config_kind` | cadena literal `f5_scenario` | obligatorio |
| raíz | `schema_version` | cadena literal `1` | obligatorio |
| raíz | `scenario_id` | ID DTC-039 | obligatorio; igual al enlace del general |
| raíz | `founder_count` | entero `u64>0` | obligatorio; igual a la suma de individuos |
| `[f5]` | `profile_metadata_path`, `profile_metadata_sha256` | ruta relativa normalizada y 64 hex minúsculas | obligatorios |
| `[f5]` | `inbreeding_f` | `f64` finito `[0,1]` | F de población; solo P=2, opcional, predeterminado `0.0` |
| `[f5]` | `max_layout_attempts`, `max_haplotype_attempts`, `max_swap_proposals`, `max_haplotype_cells` | enteros positivos con los máximos del párrafo siguiente | opcionales, con los cuatro predeterminados del párrafo siguiente |
| `[f5]` | `max_operation_units`, `max_working_memory_bytes`, `max_output_bytes` | enteros positivos no mayores que los topes F5 v1 | opcionales; presupuestos efectivos del escenario |
| `[f5]` | `excluded_regions_path`, `excluded_regions_sha256` | ruta y 64 hex minúsculas | ambos ausentes o ambos presentes |
| `[f5]` | `block_count`, `covered_bases`, `coverage_fraction` | dominios del párrafo siguiente | exactamente uno presente; no hay predeterminado |
| `[[subpopulations]]` | `population_id` | ID DTC-039 | obligatorio; único |
| `[[subpopulations]]` | `individual_count` | entero `u64>0` | obligatorio |
| `[[subpopulations]]` | `sampling` | `probabilistic` o `controlled` | opcional; predeterminado `probabilistic`; efectivo por bloque×subpoblación |
| `[[subpopulations]]` | `inbreeding_f` | `f64` finito `[0,1]` | solo P=2; opcional; sustituye el F de población; en otra ploidía debe estar ausente |
| `[[subpopulations]]` | `founder_ids`, `founder_keys` | lista de exactamente `individual_count` IDs DTC-039 | exactamente una presente; no hay predeterminado |
| `[[sampling_overrides]]` | `population_id`, `block_index` | ID DTC-039 y `u32` en `0..B-1` | pareja obligatoria y única; la población debe existir |
| `[[sampling_overrides]]` | `mode` | `probabilistic` o `controlled` | obligatorio; sustituye solo esa pareja |
| `[[sampling_overrides]]` | `controlled_max_abs_error` | `f64` finito en `[0,0.5/C]` | opcional solo con `controlled`; predeterminado `0.5/C` |

Debe existir al menos una `[[subpopulations]]`. El general puede expresar su
ploidía global predeterminada y futuras sustituciones por secuencia de
REQ-CONF-06, pero F5.0–F5.6 solo consume una ploidía efectiva homogénea: si el
general resuelve valores distintos en las secuencias elegibles, F5 rechaza el
escenario sin fingir soporte. El resolvedor completo y la prueba de
REQ-CONF-06 se difieren a F6.1; F5 no puede declararlo cerrado. Un `founder_id` aportado se usa
como `IndividualId`; una `founder_key` deriva el `IndividualId` mediante el
contrato ya vigente población×clave de REQ-IND-25, sin ordinal físico. Ni la
omisión de `seed` en el general ni un predeterminado efectivo autorizan omitirlo de
`effective-config.json` o del manifiesto.

Los cuatro límites algorítmicos son enteros positivos. Sus predeterminados y máximos son,
respectivamente, `1000/100000`, `100/10000`, `100000/10000000` y
`16777216/67108864`; excederlos es configuración inválida, no una petición de
recursos. Los tres presupuestos F5 v1 tienen como valor predeterminado y máximo,
respectivamente, `4294967296` unidades de operación (2^32), `8589934592` bytes
de memoria de trabajo (2^33) y `2147483648` bytes publicados (2^31). La persona
usuaria puede reducirlos, nunca ampliarlos dentro de esta versión.

Los dos primeros techos subieron desde 2^28 y 2^30 el 2026-09-13 por decisión
expresa del usuario, y el predeterminado subió con el máximo. El motivo está
medido: la cota de operaciones no acota la biología sino el **tamaño del
genoma**, porque el modelo de coste carga dos pasadas completas sobre el FASTA
más las bases recorridas. Con 2^28, el tomate reducido al 10 % —79 230 297 B y
78 252 003 bases— consumía 236 712 597 unidades, el 88 % del techo, antes de
abrir un solo bloque; el tomate completo (795 562 190 B) costaba 2 374 M, 8,8
veces el techo. Un `run` real de 100 individuos sobre la referencia reducida se
rechazaba en la fase `reference_and_mask_scan`, con la reserva prospectiva del
FASTA y sin llegar a abrirlo. 2^32 admite el genoma entero con unos 1,9 G de
margen para el trabajo real, y 2^33 B —8 GiB— da holgura de memoria para
genomas mayores: el pico medido en ese `run` fue de 7 897 938 B, el 0,09 % del
techo nuevo.

Predeterminado y máximo siguen siendo la misma cifra a propósito. El techo no
reserva nada: es un dispositivo de rechazo que aborta antes de abrir la
referencia, de modo que un predeterminado bajo no abarata una ejecución cara,
solo adelanta su negativa —y la negativa medida se disparaba al **leer el
fichero**, no al hacer biología—. Separar predeterminado y máximo obligaría
además a declarar el presupuesto en todo escenario que trabaje sobre un genoma
real, que es exactamente la regla que esta corrección elimina. Quien quiera el
rechazo temprano conserva la capacidad de pedirlo: basta declarar
`max_operation_units = 268435456`.

Lo que **no** cambia: el valor efectivo lo declara el escenario, no la
máquina. No se consulta RAM ni CPU detectadas, ni perfiles externos, ni
DTPD-115, para admitir. La razón es de validabilidad, no de gusto: el
presupuesto se publica en `resource_preflight` del efectivo y del manifiesto y
ACE-117 exige que un oráculo independiente recalcule las cotas; si el techo
dependiera del equipo, el validador dejaría de poder validar y la misma
configuración admitiría en una máquina y no en otra. El techo de bytes
publicados sigue en 2^31.

Se exige exactamente uno de `block_count: u32>0`,
`covered_bases: u64>0` o `coverage_fraction`. Esta última se escribe como literal
decimal TOML sin signo ni exponente, con 1 a 18 cifras decimales, en `(0,1]`.
El lector conserva el literal y forma el racional exacto `n/10^k`; está prohibido
convertirlo primero a `f64`.

`block_index`, `site_index`, `haplotype_rank`, `attempt_index` de layout y de
haplotipos empiezan en 0 y son consecutivos en su ámbito. Los límites de
intentos son cardinalidades: se ejecutan índices `0..limit-1`; `proposals` y
`accepted` son recuentos, no el índice de la última propuesta.

### Paquete de perfiles canónico

`profile.metadata.json` es RFC 8259, usa JSON Schema 2020-12 con
`additionalProperties:false` y contiene exactamente
`profile_kind="f5_generation_profile"`, `schema_version="1"`, `profile_id`,
`provenance_kind` (`built_in`, `user_supplied`), `source_citation`,
`source_license` —cadena o `null`—, `source_sha256` —cadena o `null`—, y
`files`. `files` contiene exactamente cinco objetos,
ordenados como se enumeran abajo, con `artifact_kind`,
`schema_version="1"`, `path`, `size_bytes` y `sha256`:

`profile_id` cumple DTC-039; `source_citation` es una cadena UTF-8 no vacía;
`source_sha256` usa 64 hex minúsculas cuando existe y `size_bytes` es `u64`.
Cada `path` es el nombre de fichero literal de la primera columna, relativo al
metadata, sin `..`; `sha256` son 64 hex minúsculas y los cinco paths son únicos.

| Fichero / `artifact_kind` | Cabecera TSV exacta | Dominio |
|---|---|---|
| `block_lengths.tsv` / `f5_block_length_profile` | `length_min_bp\tlength_max_bp\tweight` | enteros `2<=min<=max`; intervalos disjuntos cuya unión es continua; peso finito `>0` |
| `snp_gaps.tsv` / `f5_snp_gap_profile` | `gap_min_bp\tgap_max_bp\tweight` | enteros `1<=min<=max`; intervalos disjuntos; peso finito `>0` |
| `snp_spectrum.tsv` / `f5_snp_spectrum` | `alt_af_min\talt_af_max\tweight` | `f64` finitos `0<=min<=max<=1`; intervalos disjuntos; peso finito `>0` |
| `haplotype_weights.tsv` / `f5_haplotype_weights` | `haplotype_rank\tweight` | rangos decimales consecutivos `0..h_max-1`, `2<=h_max<=4096`; peso finito `>0` |
| `ld_profile.tsv` / `f5_ld_profile` | `population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error` | IDs DTC-039; bins cerrados contiguos desde 1, objetivos `[0,1]`; `max_abs_error` es exactamente `0.05` en v1 |

Todos los TSV usan UTF-8 sin BOM, tabulador, LF final, cabecera única, campos no
vacíos, sin comillas, comentarios, líneas vacías ni columnas adicionales. Los
enteros son decimales canónicos sin signo ni ceros iniciales; los `f64` usan el
decimal mínimo de ida y vuelta exacta. Las filas se ordenan por mínimo de
intervalo en los tres primeros TSV, por `haplotype_rank` y por
`population_profile_id`+`distance_start_bp`, respectivamente; duplicados o
desorden son error. Los `weight` son probabilidades, no pesos relativos: cada
conjunto completo de filas de longitudes, gaps, espectro o rangos haplotípicos
debe sumar 1 bajo DTC-049 en su orden canónico. Una desviación dentro de la
tolerancia general se normaliza, conserva suma/valores antes y después en la
configuración efectiva y manifiesto; otra desviación se rechaza. Una
selección ponderada obtiene U con 53 bits en `[0,1)`, recorre q en
orden de fila y elige la primera suma acumulada estrictamente mayor que U; la
última fila absorbe solo el residuo numérico posible. U se compara con el
prefijo tal cual, sin reescalar por la suma de pesos, porque los `q` son los
efectivos ya normalizados por DTC-049; un conjunto cuyo total se aparte de uno
más allá de la tolerancia se rechaza antes de sortear. La única excepción es el selector de
enteros del layout, que reescala U por la masa acumulada siempre, porque
comparte implementación con la ruta de longitudes condicionada al residual: ahí
el soporte se restringe a `[remaining*Lmin, remaining*Lmax]` y la masa factible
es estrictamente menor que uno, de modo que no reescalar sesgaría la elección
hacia el último valor. Las otras dos rutas de ese selector —los gaps y las
longitudes en modo `block_count`— reciben el soporte completo, así que allí el
reescalado es un factor uno y no cambia la semántica. Cada subpoblación debe
tener exactamente un `population_profile_id` igual a su `population_id`; no
existe remapeo implícito. Las rutas de `files` son relativas al JSON de metadata.

Cada población debe cubrir mediante sus bins todas las distancias enteras desde
1 hasta `Lmax-1`, aunque un layout concreto no las realice.

El perfil incorporado reservado `human_ld_ragsdale_approx_v1` conserva la
elección funcional previa: cada media `mu` digitalizada de `r²` de Ragsdale y
Gravel se serializa directamente como `target_mean_r2=mu`. La antigua beta con
`kappa=8.0` y sus cinco cuantiles queda como historia de la fuente, no participa
en F5 ni cambia la media. El recurso conserva cita, bytes de digitalización,
receta, valores y hashes; no es predeterminado y no se presenta como estimación
exacta de 1000 Genomas. Queda **reservado y no incorporado** desde el 13 de
septiembre de 2026, porque la fuente publica la caída de `r²_÷` en centimorgans
—desde 0,03 cM— y `ld_profile.tsv` la necesita en pares de base desde 1, de modo
que completar el rango exigiría una tasa de recombinación supuesta y una
extrapolación por debajo del primer punto medido; el registro está en F5.6.2 de
[plan-implementacion-f5-0.1.md](./plan-implementacion-f5-0.1.md).

La máscara opcional `excluded_regions.tsv` usa cabecera exacta
`region_id\tsequence_id\tstart_1based\tend_1based`, coordenadas cerradas válidas
y filas ordenadas por orden FASTA, inicio, fin e ID. Sus intervalos no se
solapan. Las bases elegibles son exclusivamente A/C/G/T fuera de máscara. Si E
es su recuento, la fracción produce `T=floor(E*n/10^k)` con enteros comprobados.
Antes de consumir RNG, un único recorrido secuencial de referencia y máscara
procesa exactamente `reference_bases_processed=sum_s sequence_length_s` bases y
produce la lista ordenada de tramos elegibles maximales; `tract_count` es su
cardinalidad exacta. El recorrido mantiene solo un búfer acotado, el intervalo
de máscara vigente y el tramo en curso: no carga la máscara completa, no
construye un bitset ni otra colección por base. Como los intervalos de máscara
son no vacíos y no se solapan, su número no puede superar las bases procesadas.
La única materialización dependiente de la fragmentación es la lista de tramos.
T=0 se rechaza. `covered_bases` fija T directamente. `block_count` fija B y deja
aleatoria la suma de longitudes; los otros modos calculan el menor B positivo
que satisface `B*Lmin<=T<=B*Lmax`. Para `block_count`, el preflight exige
`B*Lmin<=E` y que `sum_t floor(tract_length_t/Lmin)>=B`, con tramos elegibles
maximales ordenados. Para un objetivo T exige `T<=E` y
`sum_t tract_length_t>=T`; la fragmentación restante pertenece a la búsqueda
acotada, porque longitudes variables pueden cambiar el empaquetado. Si falla
una de esas desigualdades, no existe B, se excede `u64` o el perfil no contiene
un soporte continuo, el preflight devuelve `GPS-DOMAIN-101` sin consumir RNG.
El resultado queda vinculado a `scenario_id`, identidad y SHA-256 de referencia,
perfil y objetivo declarado: `block_count`, `covered_bases` o el racional
reducido exacto de `coverage_fraction`. Un consumidor no puede reutilizarlo con
otro objetivo aunque produzca el mismo B o T.

### Layout y marginales acotados

Para un objetivo de bases, cada longitud se sortea condicionada a que el resto
permanezca en `remaining_blocks*Lmin..remaining_blocks*Lmax`; en el último
bloque la longitud es el residual exacto. El muestreo conserva los pesos
renormalizados sobre valores factibles: el peso base de un entero v es
`row_weight/(row_max-row_min+1)` y después se normalizan con Neumaier en v
creciente; la elección categórica sigue la regla U53 anterior. Gaps usan la
misma regla sin condicionamiento residual. Para `block_count`, todas las
longitudes se sortean sin esa condición. El selector de longitudes recorre los
enteros factibles en orden creciente dos veces: primero para el total Neumaier
y después para el prefijo que cruza U53; no colapsa filas ni sustituye esa
semántica por una elección por banda. Su cardinalidad máxima es
`length_support_cardinality=Lmax-Lmin+1` y se presupuesta antes de RNG. En cada
intento se ordenan de mayor a menor, con `block_index` como
desempate, y cada bloque elige uniformemente un inicio entre los tramos elegibles
que admiten su longitud, actualizados tras cada colocación. «Uniformemente»
significa indexar todos los inicios válidos de todos esos tramos en orden
FASTA/coordenada, obtener un único entero acotado por su recuento total y
resolverlo por prefijos de recuentos; no se elige primero un tramo. Esa
resolución usa un recorrido directo completo de la lista ordenada por cada
bloque y por cada intento: no construye un índice auxiliar de tramos. Por ello
v3 fue sustituida por v4: las consultas de tramo se presuponen como
`2*L*(B*tract_count+B*(B-1)/2)`. El primer sumando cuenta los tramos de
preflight y el segundo los residuos que puede crear una división por bloque;
una pasada calcula el total de inicios y la otra resuelve el prefijo del entero
uniforme global, aunque un intento o la búsqueda terminen antes. V4 además
presupuesta de forma separada la comprobación de capacidad, la validación de la
lista, la copia por intento y los desplazamientos de `remove`/`insert`; por
ello no atribuye esos costes indebidamente a las dos consultas de selección. La lista se ordena
por orden FASTA e inicio. Adyacencia es válida; solape, N o máscara no lo son.
Un fallo descarta el intento completo. Tras éxito, los bloques se identifican en
orden FASTA/inicio/fin; los SNP empiezan en la primera posición y sortean gaps
hasta que el siguiente no cabe. Menos de dos SNP descarta el layout. La suma de
longitudes debe ser exactamente T cuando T existe. `ref` es la base de la
referencia y `site_alt` elige uniformemente una de las otras tres bases en orden
ASCII; no existe transición/transversión implícita. Un bloque que produciría
S>4096 sitios se rechaza con `GPS-RES-101`; junto con H<=4096 y el límite H×S,
esto acota también las estadísticas S².

Agotar `max_layout_attempts` devuelve `GPS-DOMAIN-106`, conserva número de
bloques colocados y bases pendientes del mejor intento. El mejor minimiza
`(unplaced_blocks, remaining_target_bases, placed_layout_bytes, attempt_index)`,
donde `placed_layout_bytes` serializa longitud, orden FASTA e inicio como enteros
big-endian. El error significa «no se halló layout dentro del límite», no
imposibilidad matemática. Se distingue de
`GPS-DOMAIN-101`, que solo se emite por una demostración exacta del preflight.

El espectro elige una fila por peso y, dentro de un intervalo no degenerado, un
`f64` uniforme `min+(max-min)U` con U de 53 bits en `[0,1)`; un intervalo
degenerado no consume el uniforme de valor. Produce `p_i` por
subpoblación×bloque×sitio. Para cada intento haplotípico se crea una matriz
binaria H×S sin enumerar `2^H`: para cada columna no fijada se aplica una
permutación Fisher–Yates independiente de las H filas, se recorren sus pesos
acumulados y se elige entre el prefijo anterior y el que cruza `p_i` el que
minimiza `|p_hat_i-p_i|`; un empate usa el prefijo anterior. `p_i=0/1` crea una
columna constante sin RNG. Por construcción
`marginal_max_abs_error=max_h(q_h)` y se exige
`|p_hat_i-p_i|<=marginal_max_abs_error`; con pesos uniformes equivale a `1/H`.
La única cota marginal activa es `max_h(q_h)`, y es **por bloque**: cuando H
varía entre bloques los pesos efectivos de cada uno son distintos y su máximo
también, de modo que cada fila publica su propia cota y la verifica contra sus
propias columnas.

Antes de sortear marginales se exige `h_max<=2^S` —para S>=12
queda implícito por `h_max<=4096`—; después, con V columnas cuyo objetivo no es
0/1, se exige además `h_max<=2^V`. Incumplir cualquiera devuelve
`GPS-DOMAIN-105`; la segunda comprobación ocurre tras esos sorteos pero antes de
abrir un intento haplotípico y **antes de sortear ningún recuento
haplotípico**. Las dos cotas usan `h_max` y no el H de ningún bloque: con la
cota superior el veredicto de admisión no puede depender de un sorteo, y por
tanto un mismo escenario no puede resultar válido con una semilla e inválido
con otra. El constructor repite la comprobación con el H realmente construido,
pero solo como red: la cota superior ya la implica.

#### Recuento haplotípico por bloque

`H` es por omisión el número de filas de `haplotype_weights.tsv`, igual en
todas las unidades. Un perfil con metadata `schema_version="2"` declara además
un objeto cerrado `haplotype_count_distribution` con
`algorithm="poisson-count-v1"`, `lambda` finita en `(0, 2^52]` —el dominio que
DTC-067 fija para ese muestreador, y cuyo extremo superior se rechaza **antes
de consumir RNG**, en el lector del perfil—, `h_min` y `h_max`,
y entonces **cada** bloque×subpoblación sortea su H abriendo el ámbito
`block_haplotype_count`, muestreando `Poisson(lambda)` y **pinzando** el
resultado en `[h_min, h_max]`.

La versión del metadata decide el conjunto cerrado de claves: una v1 con la
clave se rechaza y una v2 sin ella también, de modo que ningún perfil heredado
puede adquirir una distribución por descuido. Los límites son los mismos que ya
impone la tabla de pesos: `2<=h_min<=h_max<=4096`. La tabla de pesos **es** el
vector de tamaño `h_max`: se exige que traiga exactamente `h_max` filas.

Los pesos efectivos de un bloque con H haplotipos son los `H` primeros pesos
declarados **por rango**, renormalizados por su suma compensada. Cuando `H`
coincide con el tamaño declarado el vector se usa literal, sin dividir: el
perfil ya lo normalizó bajo DTC-049 y dividirlo por una suma compensada que no
vale exactamente uno movería el último bit de cada `q`. No hay un segundo
objeto aleatorio que sortear ni publicar.

El truncamiento es por pinzamiento y no por rechazo con reintentos porque el
rechazo no acota el número de sorteos y el presupuesto de operaciones sí. La
consecuencia estadística está declarada: el pinzamiento acumula en los extremos
toda la masa que la Poisson deja fuera del intervalo —`P(H=h_min)` vale
`P(N<=h_min)` y `P(H=h_max)` vale `P(N>=h_max)`—, de modo que un perfil que
quiera una cola superior rara debe elegir `lambda` holgadamente por debajo de
`h_max`.

El oráculo independiente rederiva ese sorteo con su propia reimplementación
de `poisson-count-v1`. Queda dicho que esa reimplementación es una
**transcripción literal** del muestreador productor —mismas constantes, mismas
ramas, mismo consumo—, de modo que no aporta independencia algebraica: una
constante mal transcrita se habría copiado igual. La independencia del
muestreador viene de donde ya venía, DTC-067: los vectores congelados de
`src/randomness.rs` y el oráculo ajeno de `tests/f2_rng.rs`, que no comparten
código con el productor. Lo que el oráculo F5 sí aporta por su cuenta es el
material de subsemilla, el ChaCha12 y la comparación contra los bytes
publicados.

El H de cada bloque se publica explícito en `manifest.json` como
`/f5/haplotype_counts`, junto a `/f5/haplotype_count_distribution`; las dos
claves existen o faltan juntas y solo aparecen cuando el perfil declara la
distribución. Aunque el H sea recontable desde `f5-haplotypes.tsv`, publicarlo
permite que el oráculo compare su propio sorteo contra una cifra declarada en
vez de reconstruirla. Dos filas con el mismo vector alélico son un duplicado inválido:
participan en la puntuación,
deben desaparecer mediante swaps y nunca se fusionan ni suman silenciosamente.

### LD, búsqueda y terminación

Toda aritmética binary64 de DTC-089 sigue el redondeo y la prohibición de FMA
fijados abajo para los consumidores RNG. En orden canónico `h`, `i`, `j`, y con
suma de Neumaier sin reducción paralela:

`p_hat_i=Σ_h q_h x_hi`, `p_hat_ij=Σ_h q_h x_hi x_hj`,
`D_ij=p_hat_ij-p_hat_i*p_hat_j` y
`r²_ij=D_ij²/[p_hat_i(1-p_hat_i)p_hat_j(1-p_hat_j)]`.

Denominador cero da `not_evaluable`, nunca cero. Cada bin agrega todas las
parejas intrabloque de un bloque×subpoblación cuya distancia cae en él; la media
se reduce por `(site_index_i,site_index_j)`. Cero parejas evaluables da
`not_applicable`; si todos los bins lo son, el intento es inválido. El primer
SNP tiene polaridad `positive`; cada posterior usa Bernoulli(0,5), cuyo éxito
es `positive` y cuyo fracaso es `negative`. Esa correspondencia fija la
secuencia reproducible de polaridades de una semilla, de modo que pertenece al
contrato y no a la implementación. En una pareja
evaluable se exige D>0 para polaridades iguales y D<0 para distintas; D=0
incumple. No hay cuota de signos.

Cada propuesta consume, en este orden, cuatro muestreos enteros sin sesgo para
`a in 0..H`, `b in 0..H-1`, `i in 0..S`, `j in 0..S-1`; `b>=a` y `j>=i` se
incrementan, y después se ordena cada pareja. Solo `10/01↔01/10` es candidata.
La propuesta se acepta únicamente si mejora estrictamente la tupla
`(filas_duplicadas, bins_fuera_de_tolerancia, signos_incumplidos, max_error_r2,
suma_error_r2, max_error_marginal, matrix_bytes)`. Los errores `f64` se calculan
en orden de bin y población, con Neumaier para la suma; `matrix_bytes` empaqueta
bits por filas, sitio creciente y bit más significativo primero. Un mapa de
conteos de hashes de fila se actualiza solo para a/b y contrasta bytes ante
colisión; así `filas_duplicadas=Σ_v max(0,count(v)-1)` cuesta O(S) por propuesta.
No se acepta
empate ni empeoramiento. Un intento termina al alcanzar cero en los tres campos
de error discreto —duplicados, bins y signos— y todas las tolerancias, o al
consumir `max_swap_proposals`. Se ejecutan
exactamente los intentos necesarios hasta el primer éxito o
`max_haplotype_attempts`.

El mejor intento es el mínimo lexicográfico de esa tupla, con menor
`attempt_index` como último desempate. El agotamiento devuelve
`GPS-DOMAIN-104`, conserva ese intento en diagnóstico y no afirma
imposibilidad. Un fallo exacto de dimensiones, `H*S>max_haplotype_cells`, cuota
controlada o ausencia demostrada de toda pareja LD evaluable devuelve respectivamente
`GPS-RES-101`, `GPS-DOMAIN-103` o `GPS-DOMAIN-102` antes
de fundadores. Ningún fallo reduce H/S, cambia perfil/tolerancia, fusiona bins o
acepta menos bloques.

El coste por bloque×subpoblación×intento es
`O(H*S² + max_swap_proposals*S)` tiempo: evaluar todos los pares ponderados
requiere recorrer H filas por cada pareja de sitios; las actualizaciones 2×2
aceptadas mantienen las estadísticas en `O(S)`. El contrato de admisión no deja
las cotas de recursos en notación asintótica.

#### Admisión físicamente acotada y acumulativa

`f5-resource-accounting-v9` sustituye por completo a v8 y, con él, a v7 y a los
parches intermedios v5 y v6. Frente a v8 corrige **dos defectos de la
contabilidad de salida**, ambos del lado de la cota y ninguno del techo: la
**multiplicidad** de la cota de verdad, que contaba cada SNP una vez por copia
**y secuencia**, y los **sobres por fila/elemento**, que eran constantes y por
tanto no acotaban las dos clases de registro cuyo tamaño depende de una
dimensión declarada. La admisión es una sola máquina de estados; no existe una carga previa que
pueda abrir o materializar perfil, referencia o máscara al margen de ella. Sus
estados, en este orden, son `bootstrap`, `profile_metadata`,
`profile_table_block_lengths`, `profile_table_snp_gaps`,
`profile_table_snp_spectrum`, `profile_table_haplotype_weights`,
`profile_table_ld`, `reference_and_mask_scan`, `final_formula` y `sealed`. Una
transición solo ocurre si operaciones, pico de memoria y salida siguen dentro
de sus límites inclusivos. Un exceso u overflow deja el estado en `rejected` y
prohíbe RNG, creación de destino y apertura de entradas posteriores.

El dispatcher CLI no decide F5 mediante una sonda. Abre el TOML general una vez
y lo lee completo con límite físico 65536 bytes; el byte 65537 solo sirve para
confirmar el exceso. Clasifica y valida el documento después de EOF sobre ese
mismo búfer. Esta cota se aplica antes de conocer la versión y, por tanto, a la
ruta transparente de configuración general. El escenario tiene la misma cota.
El bootstrap máximo son esos dos documentos pequeños: de él se obtiene
`max_working_memory_bytes` sin cargar perfil ni recorrer referencia. El coste
real del bootstrap se contrasta inmediatamente con los presupuestos leídos del
escenario.

Toda lectura usa un `BoundedInput` interno con contrato mínimo: rol, ruta
normalizada desde el documento propietario, límite, handle abierto, contador
`u64` y SHA-256 incremental. Lectura, medida y hash observan los mismos bytes.
Antes de ampliar un búfer se carga en el ledger su nueva reserva; al llegar al
límite se intenta leer como máximo un byte adicional. `metadata` no forma parte
de la autoridad y no permite omitir esa comprobación. Al alcanzar EOF se emite
un `PhysicalInputIdentity` con rol, ruta normalizada, tamaño observado y
SHA-256; para enlaces declarados añade tamaño y hash esperados y exige igualdad.
El orden y los roles exactos son `general_toml`, `scenario_toml`,
`profile_metadata`, `profile_block_lengths`, `profile_snp_gaps`,
`profile_snp_spectrum`, `profile_haplotype_weights`, `profile_ld`,
`reference_fasta` y, solo si existe, `mask`. Los atributos del sistema de
ficheros no forman parte de la identidad portable; si se consultan para mejorar
un diagnóstico, se obtienen del handle abierto y nunca sustituyen tamaño ni
digest observados.

`profile.metadata.json` se abre solo tras admitir el bootstrap. Con
`operation_remaining` y `memory_headroom` calculados mediante resta comprobada,
su límite es
`min(floor(operation_remaining/2),floor(memory_headroom/2))`, porque lectura y
parser y sus dos representaciones transitorias cargan dos unidades por byte. Un
límite cero rechaza sin abrir. El búfer se libera después de construir el objeto
tipado. Antes de abrir un TSV se comprueba el
coste total conocido de los cinco tamaños declarados. Cada tabla se abre
después, por orden contractual, con límite `size_bytes+1`; lectura, hash y
parseo usan el mismo búfer y este se libera antes de la tabla siguiente. Solo
persisten filas de dominio y auditorías de normalización. Un tamaño o hash
distinto es `GPS-INPUT-002`, no una oportunidad para recalcular el presupuesto.

Antes de cualquier apertura el ledger reserva todos los costes futuros
obligatorios ya conocidos. Un coste que depende de EOF, filas o geometría se
carga incrementalmente antes de consumir el siguiente byte o retener el
siguiente elemento. En particular, referencia y máscara no se abren si el
perfil y las dimensiones ya resueltas dejan sin margen alguna dimensión.

Antes de leer cada byte del primer recorrido FASTA el ledger reserva también la
visita del segundo recorrido obligatorio. La referencia y la máscara se
procesan en streaming; antes de retener un ID, fila o tramo se comprueba su
reserva. El primer recorrido produce geometría elegible e identidad física de
referencia y máscara. Después del layout, la recuperación REF abre otro handle,
lee con límite `reference.size_bytes+1` y exige el mismo tamaño y SHA-256 antes
de aceptar las bases. Este vínculo cubre cambios accidentales entre pasadas,
pero no promete seguridad frente a sustitución hostil concurrente. No se
materializa la referencia completa y la máscara no se relee.

El resultado `SealedF5Admission` posee el escenario efectivo, perfil, tramos,
dimensiones, ledger cerrado y un `F5AdmissionReceipt`. El receipt contiene
`algorithm`, identidad semántica, lista de `PhysicalInputIdentity` por rol,
límites, aportaciones `other_*`, fases y totales. Geometría solo acepta este
objeto sellado; no acepta un receipt suelto ni escalares reconstruidos por el
llamador. Así dos TOML semánticamente equivalentes pero con comentarios,
padding o tamaños distintos no comparten preflight.

Antes del primer consumo RNG v8 calcula con enteros comprobados y nombres
dimensionalmente distintos:

- `B`, número exacto de bloques resuelto desde el objetivo; `N_pop`, número de
  subpoblaciones; `H`, filas de `haplotype_weights.tsv`, que es a la vez `h_max`
  porque el lector exige esa igualdad; `A`, `Q` y `L`, los
  tres límites de intentos/propuestas; `C_total=founder_count*ploidy`;
- `L_max`, mayor longitud admitida; `gap_min`, menor gap admitido, y
  `S_bound=1+floor((L_max-1)/gap_min)`. Esta es la cota de sitios de cualquier
  bloque y se valida contra 4096 y `max_haplotype_cells` antes de sortear una
  longitud o un gap; nunca se sustituye por el S realizado después del layout;
- `gap_support_cardinality=sum_r(gap_max_r-gap_min_r+1)`, cardinalidad exacta
  de los intervalos disjuntos soportados, sin rellenar sus huecos. Cada bloque
  puede solicitar como máximo `S_bound` gaps por intento: hasta `S_bound-1`
  que caben y el gap terminal que no cabe. La cota conservadora del selector
  literal es por tanto
  `gap_selector_operation_units=2*L*B*S_bound*gap_support_cardinality`; cuenta
  completas tanto la pasada Neumaier como la pasada de prefijo, aunque esta
  última pueda terminar antes o una alternativa única no la necesite. Se
  calcula con aritmética comprobada antes de abrir cualquier ámbito RNG y no
  modifica la distribución ni su secuencia de consumo;
- `reference_bases_processed`, suma exacta de las longitudes de todas las
  secuencias resueltas que recorre el escáner; `tract_count`, número de tramos
  elegibles maximales que materializa; `sequence_count`; y
  `retained_sequence_ids`, IDs FASTA retenidos en el índice de unicidad que
  detecta duplicados. Esta última magnitud sirve para dos unidades distintas y
  esa doble lectura es deliberada: como pico simultáneo multiplica los 256
  bytes de memoria, y como total acumulado aporta una unidad de operación por
  identificador admitido. El conjunto solo crece, así que pico y total
  coinciden;
- `fasta_bytes_processed`, tamaño exacto de los bytes crudos del FASTA, que se
  reserva dos veces; y `mask_bytes_processed`, bytes crudos que atraviesa una
  vez el lector/hash de máscara (cero sin máscara);
  `K_bound=B*S_bound`; `pair_bound=S_bound*(S_bound-1)/2`;
- `bin_count(p)`, número de filas LD del `population_id=p`, y
  `bin_count_total=sum_p bin_count(p)`. Por tanto, la cota de filas de bins es
  `B*bin_count_total`, no `B*N_pop*pair_bound`;
- `profile_rows`, suma de filas de los cinco TSV de perfil, y
  `normalization_rows`, suma de filas de los cuatro conjuntos normalizados;
  ambos incluyen cada fila física aunque dos valores sean iguales;
- `general_toml_bytes`, `scenario_toml_bytes`, `profile_metadata_bytes` y
  `profile_tsv_bytes`, suma de los cinco tamaños: bytes observados por sus
  lectores, no por una consulta previa al sistema de ficheros; cada tamaño
  individual permanece en `input_identities`;
- `X=C_total*sequence_count`, número de copias por secuencia, y
  `G_bound=founder_count*K_bound`, número máximo de genotipos SNP;
- `haplotype_rows_bound=B*N_pop*h_max`, filas de `f5-haplotypes.tsv`, y
  `site_summary_rows_bound=B*N_pop*S_bound`, resúmenes por sitio. Son las dos
  únicas clases de registro cuyo ancho depende de una dimensión declarada y no
  de una longitud fija del contrato, de modo que un sobre constante no puede
  acotarlas. No son escalares publicados nuevos: se recomputan exactamente como
  `tsv_rows_bound` y `json_rows_bound`, de los que ya son sumandos. El
  presupuesto se deriva **sin ejecutar nada**, de modo que el preflight no
  puede conocer el H de ningún bloque: contabiliza `h_max`, que domina
  cualquier sorteo. Con distribución declarada la cota deja holgura —un bloque
  pinzado a `h_min` publica menos filas de las contabilizadas— y sigue siendo
  cota superior. La fórmula no cambia: `H` ya vale `h_max` por la igualdad que
  el lector del perfil exige entre el número de pesos y `h_max`.

La cota de filas de los siete TSV es
`tsv_rows_bound=B+K_bound+B*N_pop*H+B*N_pop*S_bound+
B*bin_count_total+B*C_total+B*N_pop*A`. La cota de elementos repetidos de los
objetos JSON F5 es
`json_rows_bound=founder_count+N_pop+3*B*N_pop+2*B*N_pop*H+
B*N_pop*S_bound+2*normalization_rows+9`: fundadores; subpoblaciones; muestreo,
empates y F por bloque×subpoblación; rangos candidatos/seleccionados de cada
empate en su máximo; resúmenes por sitio; valores antes/después de las cuatro
normalizaciones; siete enlaces TSV, `truth` y cabecera F5.

El preflight común entrega tres magnitudes distintas: `other_operation_units`,
`other_working_memory_bytes` y `other_output_bytes`, que cubren configuración
copiada y toda salida obligatoria o activada ajena a F5. Ninguna se reutiliza
como si tuviera otra unidad. Las cotas conservadoras F5 son:

- `operation_units = other_operation_units +
  2*(general_toml_bytes+scenario_toml_bytes+profile_metadata_bytes+
  profile_tsv_bytes) +
  2*fasta_bytes_processed +
  reference_bases_processed + mask_bytes_processed +
  retained_sequence_ids + tract_count + 2*tract_count +
  2*L*(B*tract_count+B*(B-1)/2) +
  L*(tract_count+3*B*tract_count+3*B*(B-1)/2) + profile_rows +
  normalization_rows + 2*L*B*length_support_cardinality +
  2*L*B*S_bound*gap_support_cardinality +
  L*B*(B+S_bound) +
  B*N_pop*A*(H*S_bound*S_bound+Q*S_bound) + B*C_total +
  X*K_bound + G_bound + B*N_pop*S_bound + tsv_rows_bound +
  json_rows_bound`; el factor dos de los documentos cuenta lectura/hash y
  recorrido del parser sobre el búfer, y el factor dos del FASTA cuenta sus dos
  handles; no se añade una visita de `metadata` inexistente;
- `working_memory_bytes=max(bootstrap_peak, profile_metadata_peak,
  profile_table_peaks, reference_and_mask_scan_peak,
  final_formula_peak)`. Cada registro de fase contiene
  `retained_before_bytes`, `transient_bytes`, `retained_after_bytes` y
  `phase_peak_bytes`, y cumple
  `phase_peak_bytes=other_working_memory_bytes+
  max(retained_before_bytes,retained_after_bytes)+transient_bytes`. Los
  documentos cargan dos veces sus
  bytes observados mientras coexisten búfer y parser; los cinco TSV son fases
  separadas. `final_formula_peak` conserva la cota productiva
  `other_working_memory_bytes +
  256*(4*B+2*tract_count+retained_sequence_ids+profile_rows+
  normalization_rows) +
  B*N_pop*(ceil(H*S_bound/8)+(160+S_bound)*H+256*S_bound+192*pair_bound) +
  256*B*bin_count_total + 192*X + 64*B*C_total +
  128*identity_registry_entries + 65536`, con
  `identity_registry_entries=X+B*N_pop*H+B*C_total`. Solo se suman estructuras
  que coexisten; el ledger carga prospectivamente cada crecimiento y registra
  los picos para que un lector independiente repita el máximo. Los cuatro
  términos fundadores se leen así, y v8 los corrige respecto de v7 porque
  ninguno acotaba lo que la fase productiva retiene de verdad:

  - `192*X` sustituye a `128*C_total`. El repositorio de identidades conserva
    una copia por individuo **y secuencia**, que es exactamente `X`, no
    `C_total`: contar copias del individuo dejaba fuera el segundo eje de la
    identidad y con una sola secuencia ya consumía tres cuartos del término.
    Cada entrada mide 72 bytes de copia, 24 de su índice de pertenencia y 48
    del texto `FCP_`, que es una asignación por identidad;
  - `64*B*C_total` presupuesta los asientos de asignación por
    bloque×subpoblación, que miden 32 bytes cada uno y no aparecían en ninguna
    fórmula. El factor dos cubre la holgura del vector de cada unidad y el
    texto de los `IndividualId`, que es una asignación por fundador y no por
    asiento;
  - `128*identity_registry_entries` presupuesta el registro de colisiones.
    DTC-039 exige que «toda colisión se detecta», de modo que el registro
    conserva una entrada por identidad derivada —`FCP_`, `HAP_` y `FHA_`— y no
    puede soltarla cuando la fila que la nombraba ya se ha escrito. Es el
    término dominante de la fase fundadora, porque `B*C_total` crece con el
    número de bloques aunque las tablas retenidas no lo hagan. Cada entrada es
    la pareja `(clase, texto)` de 32 bytes más la asignación de sus 30 o 31
    caracteres; el resto cubre el relleno parcial de los nodos y la cabecera
    del asignador;
  - `(160+S_bound)*H` sustituye a `64*H`: una fila haplotípica retenida cuesta
    16 bytes de identificador, 48 del texto `HAP_` y la parte proporcional de
    la unidad y su índice, y `S_bound` cubre además su fila de la matriz
    ganadora, que se materializa por filas completas de `ceil(S/8)` bytes y no
    en el empaquetado exacto que `ceil(H*S_bound/8)` describe;
- `truth_bound=65536+4096*X+2048*C_total*K_bound+4096*G_bound` y
  `output_bytes = other_output_bytes + 1048576 +
  1024*tsv_rows_bound + S_bound*haplotype_rows_bound +
  1024*json_rows_bound + 512*(ploidy+1)*site_summary_rows_bound +
  truth_bound`. Los cinco sobres son máximos del contrato F5 v1 de
  contabilidad; necesitarlos mayores es un defecto interno, no permiso para
  rebasar el presupuesto. Cada uno se deriva de los campos que el esquema
  obliga a publicar y de sus longitudes máximas —identificador de 64
  caracteres, `u64` de 20 dígitos, literal decimal de un `f64` de 24 bytes
  **por construcción**, enumerados cerrados— y **no** de ninguna ejecución:

  - `2048` por SNP dentro de una copia. **Un SNP se publica tres veces en el
    mismo registro** —modificación del segmento, SNP vigente del fragmento e
    hito del historial—, y las tres llevan su procedencia completa. El peor
    caso contractual de las tres vistas juntas es **1659 bytes**; el sobre deja
    **1,23×**. No es un sobre generoso: es ajustado, y quien lo mida contra una
    sola de las tres vistas obtendrá un tercio del coste real;
  - `4096` por registro de copia, sin sus SNP: peor caso **1566 bytes**,
    margen **2,62×**;
  - `4096` por genotipo: su lista `carriers` tiene exactamente `ploidy`
    entradas, de modo que el peor caso crece 289 bytes por copia —1091 con
    ploidía 1, **3692 con la ploidía máxima 10**— y el margen es **1,11×**.
    Tampoco es generoso, y por eso el sobre queda **atado a la ploidía máxima
    del contrato**: subirla a 12 llevaría el peor caso a 4270 bytes, y la
    implementación lo convierte en un error de compilación en vez de en una
    cota rota en silencio;
  - `1024` por fila TSV y por elemento JSON repetido. La fila más ancha es la
    de `f5-founders.tsv`, con siete identificadores y un índice: **476 bytes**,
    publicada una sola vez, margen **2,15×**. El elemento JSON más ancho de
    ancho fijo es el de muestreo por bloque×subpoblación, **419 bytes**, pero
    `sampling`, `normalizations` y `subpopulations` se publican en **los dos
    documentos** —efectivo y manifiesto— mientras `json_rows_bound` los cuenta
    una sola vez: el coste real por elemento contado es **838 bytes** y el
    margen, **1,22×**. Los resúmenes por sitio, `f_realized` y los nueve
    enlaces solo aparecen en el manifiesto, de modo que su multiplicidad es
    uno;
  - `S_bound` adicional por fila de `f5-haplotypes.tsv`. Su columna `alleles`
    tiene **exactamente `S` caracteres**, con `S_bound` hasta 4096, de modo que
    ningún sobre constante la acota: su parte constante son 242 bytes —tres
    identificadores, un `u64`, un `f64`, cinco TAB y un LF— que los 1024 cubren
    con margen 4,23×, y el resto crece byte a byte con el sitio;
  - `512*(ploidy+1)` adicional por resumen por sitio. Publica **siete** matrices
    de dosis de `ploidy+1` entradas de hasta 47 bytes, es decir `329*(ploidy+1)`
    bytes que crecen con la ploidía, más 735 bytes constantes que los 1024
    cubren con margen 1,39×. Su peor caso contractual es 1722 bytes con ploidía
    2 y **4354 con ploidía 10**, que desbordaba el sobre constante de 4096.

  **La anchura decimal de un `f64` no es una propiedad del tipo: es una
  construcción.** `f64::to_string` es **posicional y nunca exponencial**, de
  modo que su anchura no la fijan los diecisiete dígitos significativos que
  basta conservar para la ida y vuelta, sino el exponente:
  `f64::MIN_POSITIVE` ocupa **326 bytes** y `f64::MAX`, 309. Los 24 bytes de las
  derivaciones de arriba solo se garantizan por construcción, y de dos maneras
  distintas según el documento:

  - en **JSON** los escribe el serializador en notación exponencial en cuanto la
    posicional sería larga, y su anchura máxima es la del propio formato: signo,
    un dígito, el punto, los dieciséis decimales de la ida y vuelta, la `e`, su
    signo y tres dígitos de exponente;
  - en **TSV** los escribe el *literal decimal canónico F5*, que es el mínimo de
    ida y vuelta **de anchura acotada**: la forma posicional mientras quepa en 24
    bytes y la exponencial en cuanto no quepa. Las dos recuperan exactamente el
    mismo `f64` y la elección es función determinista del valor. Es la misma
    forma que la admisión exige a los perfiles y que el oráculo independiente
    reexige sobre los bytes publicados, de modo que acota además la anchura
    decimal **declarable** en una entrada. Sin ella una sola fila de
    `f5-attempts.tsv`, con sus tres columnas `f64`, podía ocupar 1243 bytes y la
    cota de 1024 dejaba de serlo.

  **Por qué la verdad se cuenta por copia y no por copia y secuencia.** Existe
  una `haplotype_copy` por copia y secuencia poblada —de ahí el término
  `4096*X`—, pero esa copia solo lleva los sitios de los bloques colocados en
  **su** secuencia. Los bloques se reparten entre secuencias sin repetirse, así
  que sumando todas las secuencias de una misma copia-de-individuo se recorren
  como máximo los `K_bound` sitios del escenario **una sola vez**. El número de
  copias-de-individuo es `C_total`, no `X`, y por tanto las incidencias SNP
  publicadas están acotadas por `C_total*K_bound`. Contar `X*K_bound`
  multiplicaba cada SNP por `sequence_count`, que es exactamente el factor que
  v8 sobrestimaba: con trece secuencias, trece veces. Quien vea el término y
  crea que falta un `X` está leyendo el contrato al revés: añadir una secuencia
  añade **registros de copia**, nunca incidencias nuevas.

  Los dos términos `fasta_bytes_processed` cubren, respectivamente, todos los bytes que el escáner clasifica y todos los que la
  pasada secuencial posterior recorre para recuperar REF y comprobar SHA-256;
  incluyen cabeceras, descripciones y separadores, no solo bases. El término
  `reference_bases_processed` cubre el avance del cursor de máscara por base y
  `mask_bytes_processed` la lectura/hash única de sus filas e IDs. No se
  materializa la referencia completa. `tract_count` añade una unidad por tramo emitido. Los
  `2*L*B*length_support_cardinality` cubren las dos enumeraciones literales de
  longitudes por intento y bloque; no se admiten recorridos no presupuestados.
  Del mismo modo,
  `2*L*B*S_bound*gap_support_cardinality` cubre conservadoramente las dos
  enumeraciones literales del soporte discontinuo de gaps por cada sorteo
  potencial, incluido el gap terminal que no cabe; no altera el selector ni
  la secuencia RNG.
  V4 añade `2*tract_count` para el recorrido de capacidad y la validación del
  request antes de RNG. Los factores
  `2*L*(B*tract_count+B*(B-1)/2)` son las dos consultas obligatorias del
  recorrido directo: una calcula el total de inicios y otra resuelve su
  prefijo, incluyendo como máximo un nuevo residuo por colocación previa. El
  término `L*(tract_count+3*B*tract_count+3*B*(B-1)/2)` cubre la copia de
  tramos por intento y los desplazamientos de `remove` y hasta dos `insert`.
  No se usa un índice auxiliar. Por tanto, con todos los demás escalares
  iguales, pasar de un tramo a 4096 aumenta operaciones exactamente en
  `4095*(3+L+5*L*B)`, memoria de trabajo en `512*4095` bytes y no cambia
  `output_bytes`, porque la lista no se publica y sus escalares caben en el
  sobre fijo. Los 256 bytes por tramo cubren la lista transitoria y los 256
  bytes por ID retenido cubren el índice de detección de duplicados; el búfer
  secuencial y el intervalo de máscara vigente están incluidos en los 65536
  bytes fijos. La lista no se publica: efectivo y manifiesto serializan sus
  tres escalares de auditoría del escaneo dentro del sobre fijo de 1048576
  bytes, por lo que no añade una fila repetida a `json_rows_bound`.

Todas las sumas, productos, `ceil`, reservas y el cálculo de `S_bound` se
ejecutan en el orden escrito con aritmética entera comprobada. El efectivo y el
manifiesto conservan `reference_bases_processed`, `fasta_bytes_processed`,
`mask_bytes_processed`, `general_toml_bytes`, `scenario_toml_bytes`,
`profile_metadata_bytes`, `profile_tsv_bytes`, las identidades físicas,
los registros de fase, `tract_count`, `retained_sequence_ids`,
`site_bound_per_block=S_bound`, `bin_count_total`, `profile_rows`,
`normalization_rows`, `tsv_rows_bound`, `json_rows_bound`,
`gap_support_cardinality`, `gap_selector_operation_units`, las tres
aportaciones `other_*` y los tres resultados; así un lector puede repetir la
cuenta sin RNG. Los receipts v5/v6/v7/v8 no son entradas admitidas ni se migran: deben
repetir la admisión v9 desde los ficheros fuente. General v3 y escenario v1 no
cambian de significado; el objeto prerelease `resource_preflight` de efectivo y
manifiesto v5 cambia de forma y exige actualización atómica de ambos esquemas y
sus lectores en F5.1.2.5, junto con el sellado del receipt que produce
`input_identities` y `phases`. La atomicidad se exige entre esquema, lectores y
productor, que es lo que evita dos formas en circulación; hasta ese hito el
esquema publica solo la proyección escalar del algoritmo v8, y su
`additionalProperties: false` es deliberado: rechaza cualquier forma intermedia.

La prueba de frontera de `max_working_memory_bytes` no se construye modificando
un TOML después de calcularlo, porque el propio literal forma parte de la
identidad y del pico de bootstrap. Usa una línea de anchura fija
`max_working_memory_bytes = <entero> # pad:<espacios>`: el entero es decimal
canónico sin ceros iniciales y el comentario compensa su anchura. El generador
itera serialización y cálculo v8 hasta un punto fijo, comprueba que dos
iteraciones consecutivas producen los mismos bytes y conserva igual longitud
física al construir la variante límite−1. El fixture no puede quedar junto a
una potencia de diez que impida esa compensación. La prueba recalcula desde los
bytes finales: el límite exacto admite, límite−1 rechaza, y cualquier edición de
padding cambia SHA-256 y obliga a repetir la admisión.

Cada resultado se compara inclusivamente con
`max_operation_units`, `max_working_memory_bytes` y `max_output_bytes` y con los
topes duros F5 v1 anteriores. `GPS-RES-101` registra algoritmo, fase, rol,
dimensión, acumulado/reserva, límite y primer término; un exceso físico indica
además bytes observados antes del rechazo. `GPS-INPUT-001` conserva fallos de
apertura/lectura y `GPS-INPUT-002`, desigualdad de tamaño o SHA-256. Los errores
TOML/JSON/TSV solo se emiten después de la lectura acotada. Esta admisión
determinista no ejecuta perfiles externos,
no consulta RAM/CPU disponible y no promete tiempo de pared o RSS. Es distinta
de DTPD-115: F9 podrá fijar objetivos de rendimiento y publicación, pero no es
prerrequisito de F5 ni puede relajar estos topes de dominio. Está prohibido
materializar o recorrer el conjunto potencia de las H filas.

Los bytes usados por desempates, hashes o vectores no dependen del serializador
TSV. `placed_layout_bytes` es la concatenación, por bloque en orden
FASTA/inicio/fin/ID, de cuatro `u64` big-endian sin signo:
`sequence_ordinal`, `start_1based`, `end_1based`, `length_bp`. `matrix_bytes`
empieza con `H` y `S` como `u64` big-endian y continúa con H filas en
`haplotype_rank` creciente; cada fila empaqueta `site_index` creciente, primer
bit en el bit 7 de cada octeto, y rellena con ceros los bits bajos del último
octeto. `alleles` es exactamente la vista ASCII de esos S bits (`0` o `1`), sin
separador ni final de línea. Toda matriz, layout o vector de ID usa estas
codificaciones exactas antes de SHA-256 o BLAKE3; cambiar anchura, endianidad,
orden o relleno exige nueva versión de algoritmo y de artefacto.

### Fundadores, F y factibilidad controlada

El modo probabilístico muestrea una fila por copia con q. En diploides, primero
muestrea la copia 0; con probabilidad F la copia 1 reutiliza exactamente esa fila
y, si no, se muestrea independientemente. Esto produce, para categorías i/j,
`P(ii)=q_i²+F*q_i*(1-q_i)` y `P(ij)=2*(1-F)*q_i*q_j`, i<j; en el caso bialélico
son `p²+Fpq`, `2pq(1-F)` y `q²+Fpq`. F=0 coincide con HWE y F=1 produce
homocigosis completa. P=1 o P>2 con F es error. `founder_f` se consume antes
que `founder_sample` de la copia 1; si decide reutilizar, ese segundo flujo no
se abre. F=0/1 no consume `founder_f`, conforme a la regla Bernoulli del
catálogo.

Controlado calcula `C=individual_count*ploidy`, cuotas ideales `C*q_h`, pisos y
reparte el residual por mayor resto. Los restos binary64 se ordenan de mayor a
menor; igualdad significa el mismo valor finito tras el cálculo normativo. Si
un grupo de restos iguales cruza la frontera del residual, sus
rangos se ordenan primero y se aplica una única Fisher–Yates al grupo con
`controlled_quota_tie`, seleccionando el prefijo necesario. Grupos incluidos o
excluidos por completo no consumen RNG. Antes de ese desempate se comprueba la
cuota fija de cada rango y, para el grupo fronterizo, los dos estados posibles
—piso o piso+1—; como todos comparten resto, esa comprobación O(H) demuestra la
cota para todas las combinaciones sin enumerarlas. Para la pareja
bloque×subpoblación se exige suma C y
`|n_h/C-q_h|<=controlled_max_abs_error_effective` en ambos estados aplicables;
su valor predeterminado es `0.5/C` y uno declarado solo puede ser más estricto.
La cuota ideal normativa es `fl(C*q_h)` —el producto evaluado en binary64 sobre
el `q_h` ya normalizado, no el real `C*q_h`— y la desigualdad se evalúa **en
copias**: `|n_h - fl(C*q_h)| <= t*C`, donde el predeterminado aporta el valor
`0.5` **exacto** y no `(0.5/C)*C`, que no vuelve a valer `0.5` para todo `C`.
El motivo es que `0.5/C` es exactamente el umbral de redondeo al entero más
próximo, y ese umbral solo es satisfacible en aritmética exacta: en cuanto `q`
es una normalización binary64, un ideal semientero cae un ulp fuera y perfiles
corrientes se quedan sin ninguna cuota legal. Está verificado con
`q={0.7,0.3}`, cinco fundadores haploides y la tolerancia predeterminada: en
decimales exactos la cuota `[3,2]` cumple la cota con igualdad —`|3/5-7/10|` y
`|2/5-3/10|` valen ambos `1/10`—, mientras que evaluando `|n_h/C-q_h|` sobre los
pesos binary64 **ninguna** de las seis cuotas enteras la cumple. La consecuencia
se declara sin suavizarla: un validador independiente que compruebe la
desigualdad literal `|n_h/C-q_h|<=t` sobre los pesos binary64 puede encontrar
cuotas publicadas que la exceden en un ulp; no es un incumplimiento, sino que la
desigualdad normativa es la de copias y quien valide debe usar esa.
Si falla, devuelve
`GPS-DOMAIN-103` antes de RNG. Una permutación Fisher–Yates de las C identidades de
copia, inicialmente ordenadas por `individual_id` y `copy_index`, asigna el
multiconjunto ordenado por `haplotype_rank`. No se prueba, corrige ni relaja otra
cuota y F>0 es incompatible.

Por población×bloque diploide se informa
`F_realized=1-H_observed/H_expected`, donde
`H_observed=heterozygous_individuals/individual_count`,
`q_realized_h=n_h/C` y `H_expected=1-Σ_h q_realized_h²`; se usa la frecuencia
realizada, no el objetivo q. Si `H_expected=0`, es `not_applicable`. F5 no
reproduce y el Ne reproductivo **no le aplica**: no lo calcula ni lo registra,
y ningún documento F5 lo publica. La versión anterior de esta frase decía que
F5 «registra Ne reproductivo como `not_applicable`», lo que era irrealizable y
queda corregido aquí: la enumeración del objeto `f5` del manifiesto general,
fijada por esta misma decisión, no contiene ningún campo donde escribirlo;
`manifest-general-v5.schema.json` es `additionalProperties:false` y tampoco lo
tiene; y REQ-MAN-92, que enumera lo que F5 debe registrar, no lo pide. La
fórmula de Ne por contribuciones pertenece a F7.3, que es la primera fase con
parentales y contribuciones gaméticas. En F7.3, con `N_p>=2`
parentales elegibles, contribuciones gaméticas realizadas `k_j`,
`k_bar=Σk_j/N_p` y varianza poblacional
`V_k=Σ(k_j-k_bar)²/N_p`, se usa
`Ne=(N_p*k_bar-1)/(k_bar-1+V_k/k_bar)`. `k_bar<=0`, denominador no positivo o
datos incompletos dan `not_evaluable` con causa; Ne nunca es objetivo ni cambia
descendientes. Con ello DTPD-109 queda resuelta.

### IDs, RNG y formatos de salida

Todos los IDs automáticos usan DTC-039, no concatenación textual. Cada fila del
catálogo tiene exactamente estos componentes con nombre y valor UTF-8, en este
orden:

| Prefijo / `class` | Componentes exactos tras `class` |
|---|---|
| `BLK_` / `f5_block` | `execution_id`, `scenario_id`, `reference_id`, `block_index_u64` |
| `SITE_` / `f5_site` | `execution_id`, `scenario_id`, `reference_id`, `block_index_u64`, `site_index_u64` |
| `HAP_` / `f5_haplotype` | `execution_id`, `scenario_id`, `reference_id`, `population_id`, `block_index_u64`, `haplotype_rank_u64`, `alleles` |
| `FCP_` / `f5_founder_copy` | `execution_id`, `scenario_id`, `reference_id`, `population_id`, `individual_id`, `sequence_id`, `copy_index_u32` |
| `FHA_` / `f5_block_haplotype_assignment` | `execution_id`, `scenario_id`, `reference_id`, `population_id`, `block_index_u64`, `founder_copy_id`, `haplotype_id` |

`class` es el literal de la segunda columna, no el prefijo. Los enteros son
decimales canónicos UTF-8;
cada nombre y valor conserva su prefijo de longitud DTC-037, por lo que campos
adyacentes nunca resultan ambiguos. Los prefijos son `BLK_`, `SITE_`, `HAP_`,
`FCP_` y `FHA_`; se toman 128 bits BLAKE3 y Base32 sin relleno. `IndividualId`
identifica al fundador una sola vez en toda la ejecución; `FCP_` identifica su
copia biológica por secuencia y fase a través de todos los bloques de esa
secuencia. `FHA_` no es una copia: es la asignación por bloque a un `HAP_` y
nunca puede sustituir la identidad persistente de copia. Orden físico, ruta, hilo,
intento y RNG no forman parte del ID. `alleles` solo se incluye para `HAP_`; el
rango canónico distingue filas y el vector protege el contenido, mientras la
validación impide publicar duplicados. Colisión detectada reutiliza
`GPS-DOMAIN-002`, cuyo significado general ya es colisión de identidad.
`id-catalog-v4` es el superconjunto estricto de v3 que añade únicamente estas
cinco clases F5, sus componentes y vectores; conserva literalmente todas las
filas anteriores y usa versión de derivación 1 mientras no cambien contexto ni
codificación DTC-039. Antes de F5.1 publica por cada clase al menos un vector
completo con componentes, material binario codificado, digest BLAKE3 de 32
bytes, primeros 128 bits, Base32 y ID final prefijado.

`f5-rng-catalog-v1` usa DTC-037/ChaCha12 y exactamente estos ámbitos;
`f5-rng-catalog-v2` conserva los doce literales y añade el decimotercero,
`block_haplotype_count`:

| decisión | componentes adicionales en orden | sal |
|---|---|---|
| `layout_length` | `scenario_id`, `stage`, `attempt_index_u64`, `block_index_u64` | vacía |
| `layout_position` | `scenario_id`, `stage`, `attempt_index_u64`, `block_index_u64` | vacía |
| `layout_gap` | `scenario_id`, `stage`, `attempt_index_u64`, `block_index_u64` | vacía |
| `site_alt` | `scenario_id`, `block_id`, `site_id`, `stage` | vacía |
| `marginal_spectrum` | `scenario_id`, `population_id`, `block_id`, `site_id`, `stage` | vacía |
| `initial_column` | `scenario_id`, `population_id`, `block_id`, `stage`, `attempt_index_u64` | `site_index_u64` |
| `polarity` | `scenario_id`, `population_id`, `block_id`, `stage`, `attempt_index_u64` | `site_index_u64` |
| `swap_proposal` | `scenario_id`, `population_id`, `block_id`, `stage`, `attempt_index_u64` | vacía; un flujo secuencial por intento |
| `founder_sample` | `scenario_id`, `population_id`, `block_id`, `individual_id`, `copy_index_u32`, `stage` | vacía |
| `founder_f` | `scenario_id`, `population_id`, `block_id`, `individual_id`, `stage` | vacía |
| `controlled_quota_tie` | `scenario_id`, `population_id`, `block_id`, `stage` | vacía; un flujo secuencial por bloque |
| `founder_controlled_permutation` | `scenario_id`, `population_id`, `block_id`, `stage` | vacía |
| `block_haplotype_count` | `scenario_id`, `population_id`, `block_id`, `stage` | vacía; solo en v2 |

`block_haplotype_count` **no** lleva `attempt_index_u64`: H es una dimensión de
la unidad, no del intento, y un H que cambiara entre intentos daría matrices de
alturas distintas para el mismo bloque×subpoblación, con lo que ningún rango
haplotípico sería estable para los fundadores ni para los artefactos.

Una ejecución publica el catálogo del que salen **sus** decisiones, no el más
alto que el binario implemente: un perfil que no declara la distribución del
recuento haplotípico no puede abrir `block_haplotype_count` y publica
`f5-rng-catalog-v1`, con su ruta, su tamaño y su huella intactos. Es lo que
mantiene byte a byte el `manifest.json` y el `effective-config.json` de todos
los paquetes anteriores a esa ampliación.

`stage` es literalmente el nombre de decisión de la primera columna. Los
enteros de componente y sal son decimales canónicos UTF-8 dentro de DTC-037.
Una sal no vacía se codifica como el último componente DTC-037, con nombre ASCII
exacto `salt.<nombre_de_la_columna>` —aquí `salt.site_index_u64`— y su valor;
una sal vacía no añade componente.
Cada combinación distinta de componentes y sal deriva una subsemilla; una fila
sin sal abre un único flujo secuencial en el orden de consumo descrito, y una
fila con sal abre un flujo independiente por valor de sal. No se añade una sal
textual redundante. Los
enteros
acotados usan rechazo desde `next_u64`: se descartan valores `>=floor(2^64/n)*n`
y se devuelve `value%n`. Un U53 toma `bits53=next_u64()>>11` y
`U=bits53/2^53` en binary64. Fisher–Yates recorre exactamente `i=len-1..1` y
sortea `j` en `0..=i` con ese entero acotado. Bernoulli compara `bits53` con
`floor(p*2^53)`; p=0/1 no consume RNG. Toda operación binary64 se evalúa en el
orden escrito con redondeo IEEE 754 *round-to-nearest ties-to-even* entre
operaciones y sin contracción FMA.

El consumo por decisión también es cerrado: `layout_length`/`layout_gap`
consumen un U53 por elección entre dos o más valores factibles;
`layout_position` y `site_alt`, un entero acotado; `marginal_spectrum`, un U53
para fila y otro para valor solo si el intervalo no es degenerado;
`initial_column`, H−1 enteros por permutación no constante; `polarity`, un U53
por sitio posterior al primero; `swap_proposal`, los cuatro enteros ya
enumerados por propuesta; `founder_sample`, un U53 por muestreo independiente;
`founder_f`, un U53 solo para 0<F<1; y la permutación controlada, C−1 enteros.
`controlled_quota_tie` consume g−1 enteros solo para el único grupo de tamaño g
que cruza la frontera. `block_haplotype_count` consume un único muestreo
`poisson-count-v1` por bloque×subpoblación, y solo si el perfil declara la
distribución. Una elección con una sola alternativa no consume RNG.
El catálogo debe publicar al menos un vector completo por fila —material
codificado, subsemilla y primeras cuatro palabras `u64`— antes de F5.1.

Los siete TSV exactos de un éxito son:

| Fichero / `artifact_id` / `artifact_kind` | Cabecera exacta | Clave/orden |
|---|---|---|
| `f5-block-layout.tsv` / `f5-01-block-layout-v1` / `f5_block_layout` | `block_id\tsequence_id\tstart_1based\tend_1based\tlength_bp` | orden FASTA, inicio, fin, ID |
| `f5-sites.tsv` / `f5-02-sites-v1` / `f5_sites` | `block_id\tsite_index\tsite_id\tsequence_id\tposition_1based\tref\talt` | bloque, índice |
| `f5-haplotypes.tsv` / `f5-03-haplotypes-v1` / `f5_haplotypes` | `population_id\tblock_id\thaplotype_id\thaplotype_rank\tq\talleles` | población, bloque, rango |
| `f5-frequencies.tsv` / `f5-04-frequencies-v1` / `f5_frequencies` | `population_id\tblock_id\tsite_id\tpolarity\ttarget_alt_af\trealized_alt_af\tabs_error` | población, bloque, sitio |
| `f5-ld-evaluation.tsv` / `f5-05-ld-evaluation-v1` / `f5_ld_evaluation` | `population_id\tblock_id\tdistance_start_bp\tdistance_end_bp\tstatus\tpair_count\tmean_r2\ttarget_mean_r2\tmax_abs_error\tabs_error` | población, bloque, inicio de bin |
| `f5-founders.tsv` / `f5-06-founders-v1` / `f5_founders` | `population_id\tindividual_id\tsequence_id\tcopy_index\tblock_id\tfounder_copy_id\tblock_assignment_id\thaplotype_id` | población, individuo, secuencia, copia, bloque |
| `f5-attempts.tsv` / `f5-07-attempts-v1` / `f5_attempts` | `population_id\tblock_id\tattempt_index\tproposals\taccepted\tduplicate_rows\tbins_outside\tsigns_failed\tmax_r2_error\tsum_r2_error\tmax_marginal_error\tselected_best` | población, bloque, intento |

Rigen el mismo contrato byte a byte de los perfiles. `alleles` contiene
exactamente S caracteres `0`/`1`; la `polarity` propia de cada población es
`positive` o `negative` y por ello reside en `f5-frequencies.tsv`, no en la
geometría común de `f5-sites.tsv`;
`status` es `evaluated` o `not_applicable`. Solo en este último caso
`mean_r2` y `abs_error` contienen el literal `NA`; ningún otro campo queda vacío.
`selected_best` es `true/false` y existe una sola fila `true` por
bloque×población. Frecuencias, errores y medias usan decimal `f64` de ida y vuelta.
Todos los IDs referenciados existen en su tabla propietaria y tienen la clase
indicada arriba; índices, recuentos, coordenadas y longitudes son enteros
decimales canónicos no negativos, con coordenadas/longitudes positivas,
`length_bp=end-start+1`, `site_index` y `haplotype_rank` consecutivos desde 0,
`copy_index` en `0..ploidy-1`, `q` finito en `(0,1]` y suma 1 por
población×bloque, frecuencias/errores/medias en `[0,1]`, y cada contador de
intento cabe en los límites efectivos. `pair_count=0` si y solo si el estado es
`not_applicable`; si es `evaluated`, es positivo y su error respeta el perfil.

El `effective-config.json` vigente añade exactamente este objeto `f5`; todos los
objetos rechazan propiedades adicionales:

| Clave | Tipo y dominio |
|---|---|
| `schema_version` | cadena exacta `"1"` |
| `scenario_id`, `reference_id` | ID DTC-039 |
| `reference_sha256`, `seed` | 64 hex minúsculas |
| `ploidy`, `founder_count` | enteros validados |
| `subpopulations` | lista ordenada por `population_id` de objetos exactos `population_id`, `individual_count`, `inbreeding_f`, `inbreeding_f_source`, `sampling`, `founders`; `inbreeding_f` es `f64` en P=2 y `null` en otra ploidía; la fuente es `population` o `subpopulation`; `founders` contiene objetos `individual_id` y `founder_key` cadena o `null`, ordenados por ID |
| `profile` | objeto exacto `profile_id`, `path`, `size_bytes`, `sha256` |
| `mask` | `null` u objeto exacto `path`, `size_bytes`, `sha256` |
| `target` | objeto exacto `kind`, `declared_value`, `numerator`, `denominator`, `target_bases`; `kind` es `block_count`, `covered_bases` o `coverage_fraction`; campos no aplicables son `null` |
| `sampling` | lista ordenada por `population_id`+`block_id`, con una fila exacta por bloque×subpoblación: `population_id`, `block_id`, `block_index`, `requested_mode`, `effective_mode`, `source` (`subpopulation_default` o `block_override`), `controlled_max_abs_error_declared` y `controlled_max_abs_error_effective`; las tolerancias son `null` en probabilístico y la efectiva es finita en `[0,0.5/C]` en controlado, con `C=individual_count*ploidy` de esa subpoblación |
| `normalizations` | cuatro objetos exactos en orden `block_lengths`, `snp_gaps`, `snp_spectrum`, `haplotype_weights`, cada uno con `set_kind`, `canonical_order`, `algorithm="neumaier-f64-v1"`, `tolerance`, `input_values`, `input_sum`, `input_delta`, `effective_values`, `output_sum`, `output_delta`; los números recuperan el mismo `f64` y materializan físicamente la auditoría DTC-049 |
| `limits` | objeto exacto con los cuatro límites algorítmicos y los tres presupuestos efectivos TOML |
| `resource_preflight` | objeto exacto `algorithm="f5-resource-accounting-v9"`, `input_identities`, `phases`, `reference_bases_processed`, `fasta_bytes_processed`, `mask_bytes_processed`, `general_toml_bytes`, `scenario_toml_bytes`, `profile_metadata_bytes`, `profile_tsv_bytes`, `tract_count`, `retained_sequence_ids`, `gap_support_cardinality`, `gap_selector_operation_units`, `site_bound_per_block`, `bin_count_total`, `profile_rows`, `normalization_rows`, `tsv_rows_bound`, `json_rows_bound`, `other_operation_units`, `other_working_memory_bytes`, `other_output_bytes`, `operation_units`, `working_memory_bytes`, `output_bytes`, `operation_limit`, `working_memory_limit_bytes`, `output_limit_bytes`; `input_identities` se ordena por rol contractual y cada objeto exacto contiene `role`, `normalized_path`, `observed_size_bytes`, `observed_sha256`, `declared_size_bytes` y `declared_sha256`, estos dos últimos `null` cuando no existen en el documento propietario; `phases` sigue el orden cerrado de la máquina de admisión y cada objeto contiene `phase`, `operation_units_cumulative`, `retained_before_bytes`, `transient_bytes`, `retained_after_bytes`, `phase_peak_bytes` y `output_bytes_reserved`; todos los recuentos son `u64`, los valores calculados no exceden sus límites y efectivo/manifiesto son iguales escalar y posicionalmente |
| `algorithms` | objeto exacto `layout="f5-layout-v1"`, `marginals="f5-spectrum-v1"`, `haplotypes="weighted-prefix-swap-v1"`, `ld="phased-r2-neumaier-v1"`, `controlled="largest-remainder-v1"`, `inbreeding="diploid-f-v1"`, `truth_projection="f5-founder-truth-projection-v2"` |
| `rng_catalog`, `id_catalog` | objetos exactos `catalog_id`, `path`, `size_bytes`, `sha256` |

En `normalizations`, `canonical_order` vale exactamente
`length_min_bp`, `gap_min_bp`, `alt_af_min` o `haplotype_rank` según
`set_kind`. `input_values` y `effective_values` son listas de objetos exactos
`key`+`value` y conservan todas las claves canónicas y todos los pesos en ese
orden; no se sustituyen por un digest ni por un resumen. El manifiesto copia
estos cuatro objetos con igualdad semántica exacta respecto al efectivo.

El objeto `f5` del `manifest.json` general contiene exactamente
`schema_version`, `scenario_id`, `seed`, `reference_id`, `reference_sha256`, `profile`,
`mask`, `eligible_bases:u64`, `target` —igual al efectivo—, `ploidy:u8`,
`founder_count:u64`, `subpopulations`, `algorithms`, `rng_catalog`, `id_catalog`, `limits`,
`resource_preflight`, `normalizations`,
`attempt_summary`, `sampling`, `quota_ties`, `f_realized`,
`initialization_summaries`, `truth` y `artifacts`. `attempt_summary`
contiene exactamente `layout_attempts_used:u64`,
`haplotype_attempts_total:u64`, `max_attempt_index:u64`; `f_realized` se ordena
por población+bloque y contiene objetos exactos `population_id`, `block_id`,
`target_f` (`f64` o `null`), `status` (`evaluated`/`not_applicable`) y `value`
(`f64` o `null`); P distinto de 2 usa ambos valores `null` y `not_applicable`.
`seed` cumple `^[0-9a-f]{64}$` y es idéntica a
`effective-config.json#/seed` y `effective-config.json#/f5/seed`. `profile`,
`mask`, `target`, `algorithms`, `rng_catalog`, `id_catalog`,
`limits`, `resource_preflight` y `normalizations` reutilizan exactamente los
objetos del efectivo. En particular, cada escalar de `resource_preflight` —no
solo los tres totales— debe ser idéntico entre efectivo y manifiesto; ambos se
rechazan si una copia discrepa aunque la fórmula pudiera producir el mismo
total. `subpopulations` se
ordena por ID y cada elemento contiene exactamente `population_id`,
`individual_count`, `inbreeding_f` e `inbreeding_f_source` (`f64` o `null` con
la misma regla de ploidía); no repite claves fundadoras ni el vector de
individuos. `sampling` reutiliza la lista efectiva bloque×subpoblación y
`schema_version` es
la cadena `"1"`. `quota_ties` es una lista ordenada por población+bloque;
queda vacía en probabilístico y, en controlado, contiene un objeto por frontera
empatada con exactamente `population_id`, `block_id`, `remainder` binary64,
`candidate_ranks` en orden ascendente, `select_count:u64` y `selected_ranks` en
el orden sorteado. `candidate_ranks`+`select_count` representa de forma
compacta todas las combinaciones óptimas; está prohibido enumerar sus
subconjuntos.
`initialization_summaries` se ordena por población, bloque y sitio y contiene
exactamente `population_id`, `block_id`, `site_id`, `ploidy`, `model`,
`requested_mode`, `effective_mode`, `input_alt_af`, `realized_alt_af`,
`expected_genotype_probabilities`, `expected_genotype_frequencies`,
`expected_genotype_counts`, `realized_genotype_frequencies`,
`realized_genotype_counts`, `recounted_genotype_frequencies`,
`recounted_genotype_counts`, `recount_source` y `recount_matches`.
`model` es `hwe`, `diploid_f` o `controlled_hwe_conditioned_counts`. Cada vector
es una lista de objetos `{dosage,value}` ordenada por dosis `0..ploidy`; los
valores de probabilidad/frecuencia son `f64`, los recuentos esperados son
`f64=N*probability` y los realizados/recontados son `u64` que suman
`individual_count`. En controlado, la probabilidad esperada por dosis es la
hipergeométrica inducida por permutar sin reemplazo las C copias con el recuento
ALT fijado; no se vuelve a presentar como binomial no condicionada.
`recount_source="outputs/truth.jsonl#f5-founder-truth-projection-v2"`; el oráculo que
produce `recounted_*` solo lee la proyección publicada v2. Agrupa por el par
estructurado `site_id`×`population_id`, exige exactamente los individuos y la
ploidía declarados en la cabecera v2 y calcula cada dosis desde los portadores;
no abre los siete TSV ni usa manifiesto, contadores, asignaciones o acumuladores
del productor. `recount_matches` debe ser `true`.
Así REQ-MAN-70 conserva probabilidades, frecuencias y recuentos esperados,
realizados y recontados de forma independiente por locus×subpoblación.

`truth` es el enlace obligatorio a la proyección fundadora canónica y contiene
exactamente `artifact_id="truth-jsonl"`, `artifact_kind="truth_jsonl"`,
`schema_version="2"`, `projection_version="f5-founder-truth-projection-v2"`,
`path="outputs/truth.jsonl"`, `size_bytes` y `sha256`.

`truth-v2.schema.json` es un esquema nuevo y autocontenido; no modifica ni usa
como extensión abierta `truth-v1.schema.json`, cuyo `additionalProperties:false`
y semántica permanecen byte a byte históricos. V2 conserva los registros y
significados reconstruibles de DTC-061 y exige además los siguientes objetos
cerrados para la proyección F5:

- la cabecera `schema_version="2"` incorpora `projection`, con exactamente
  `kind="f5_founder"`, `projection_version="f5-founder-truth-projection-v2"` y
  `scenario_id`, y `populations`, ordenada por `population_id`, con objetos
  exactos `population_id`, `individual_count` y `ploidy`;
- cada `haplotype_copy` incorpora `founder_membership`, con exactamente
  `population_id` y `founder_copy_id`. Su `identity` conserva la tupla
  `(reference_id,individual_id,sequence_id,copy_index)`, por lo que el lector
  puede unir cada portador a una única subpoblación sin inferirla del ID;
- cada `genotype` SNP incorpora `f5_locus`, con exactamente `population_id`,
  `block_id`, `site_id`, `ref` y `alt`. La población coincide con la membresía
  de todas sus copias portadoras, el sitio pertenece al bloque y `ref`/`alt`
  coinciden con locus y alelos tipados.

La proyección emite una `haplotype_copy` por tupla FCP y un genotipo SNP por
fundador y sitio a partir de FHA→HAP. Las tuplas FCP son exactamente el producto
de las subpoblaciones, sus fundadores, las secuencias que el layout **puebla** y
`copy_index` en `0..P-1`: una secuencia que la referencia retiene pero que no
recibe ningún bloque no aporta ninguna copia, igual que no aporta ninguna fila a
`f5-founders.tsv`. La consecuencia es deliberada y no se suaviza: en una
ejecución cuya referencia retiene 25 secuencias y cuyo layout puebla 24, ningún
fundador tiene copia de la secuencia 25 en `truth.jsonl`, aunque la referencia la
retenga y las fases posteriores la hereden; un lector que espere una copia por
secuencia retenida está leyendo mal el contrato. La alternativa —materializar
también las copias vacías— obligaría a que la admisión sellara los
identificadores de secuencia y no solo su recuento, y rechazar la divergencia
tumbaría configuraciones legítimas, porque con muchas secuencias y pocos bloques
lo normal es que la mayoría no reciba ninguno. Se ordenan primero las copias por
referencia, población, individuo, secuencia y fase, y después los genotipos por
referencia, población, bloque, sitio e individuo. El manifiesto no puede
declarar éxito si falta el fichero, su hash no coincide o un reconstructor que
solo lee `truth.jsonl` v2 no recupera exactamente identidad, pertenencia, fase,
alelos, dosis y recuentos por locus×subpoblación. Publicar v1 con una proyección
F5 es error de versión, no compatibilidad implícita.
`artifacts` contiene exactamente siete enlaces en el orden anterior —que es
también el orden lexicográfico de `artifact_id`—, cada uno con exactamente
`artifact_id`, `artifact_kind`, `schema_version="1"`, `path`, `size_bytes` y
`sha256`; sus siete valores literales son los de la tabla y no duplica filas
TSV. Los tamaños son `u64`; toda ruta es relativa y normalizada.

El catálogo público F5 añade sin reutilizar significado:
`GPS-CONF-101` escenario/one-of/F/P/límite inválido; `GPS-FMT-101` perfil,
máscara, TSV, JSON, versión, orden o hash inválido; `GPS-DOMAIN-101` layout
demostrablemente inviable; `GPS-DOMAIN-102` ninguna pareja LD evaluable;
`GPS-DOMAIN-103` cuota controlada inviable; `GPS-DOMAIN-104` búsqueda agotada;
`GPS-DOMAIN-105` H vectores únicos imposibles para S;
`GPS-DOMAIN-106` búsqueda de layout agotada;
`GPS-RES-101` presupuesto/dimensión excedido, incluido el overflow comprobado
de una fórmula de preflight provocado por dimensiones declaradas; y
`GPS-INTERNAL-101` overflow o resultado numérico no finito en una operación
posterior que el preflight ya había demostrado acotada. Cada diagnóstico incluye
escenario y, cuando existen, población, bloque, sitio, campo, ruta, intento,
límite y valor efectivo. La corrupción de un artefacto publicado usa el código
de formato vigente y las colisiones siguen usando `GPS-DOMAIN-002`.

Un rechazo de preflight devuelve diagnóstico JSON DTC-044 y no abre ni publica
destino. Un fallo controlado posterior publica, conforme a DTC-045, un manifiesto
`fallido` con un objeto `f5_failure` que contiene exactamente `code`, `stage`,
`scenario_id`, `population_id`, `block_id`, `limits`, `counters`, `best_layout`,
`best_attempt` y `artifacts`; los tres IDs, límites, contadores y mejores
intentos usan los mismos tipos del fichero diagnóstico descrito a continuación.
`artifacts` contiene exactamente un enlace no biológico a
`f5-failure.json`; no contiene `truth` ni lista ninguno de los siete TSV, aunque existieran
temporales. `f5-failure.json` tiene `diagnostic_kind="f5_failure"`,
`schema_version="1"` y exactamente `code`, `stage`, `scenario_id`,
`population_id`/`block_id` —cadena o `null`—, `limits`, `counters`,
`best_layout`, `best_attempt` y `artifacts=[]`. `limits` repite exactamente los
cuatro límites algorítmicos y los tres presupuestos efectivos; `counters` contiene `layout_attempts_used`,
`haplotype_attempts_used`, `swap_proposals` y `swap_accepted`, todos `u64`.
`best_layout` es `null` o contiene `attempt_index`, `placed_blocks`,
`unplaced_blocks`, `placed_bases`, `remaining_target_bases` y
`layout_bytes_sha256`. `best_attempt` es `null` o contiene `attempt_index`,
`duplicate_rows`, `bins_outside`, `signs_failed`, `max_r2_error`,
`sum_r2_error`, `max_marginal_error` y `matrix_bytes_sha256`. Su
enlace usa exactamente `artifact_id="f5-00-failure-v1"`,
`artifact_kind="f5_failure"`, `schema_version="1"`, ruta, tamaño y SHA-256.
Un éxito publica conjuntamente los siete TSV y la proyección fundadora
`outputs/truth.jsonl`, o falla atómicamente. El oráculo independiente relee esos
bytes, reconstruye H×S, q, marginales, D, r², bins, signos, puntuaciones,
cuotas, F, IDs, hashes, mejor intento y estado fundador sin importar estado del
productor. Marimo solo importa/exporta este TOML, invoca Rust sin *shell* y
muestra su respuesta.

Esta decisión cierra DTPD-110: los únicos catálogos distributivos F5 son
intervalos enteros ponderados uniformes para longitudes/gaps, intervalos `f64`
ponderados uniformes para marginales, pesos categóricos para H y medias `r²`
por bin; las tolerancias son DTC-049 para normalización,
`max_h(q_h)` para marginales, `0.05` inclusivo para LD y, para controlado, la
tolerancia inclusiva efectiva de cada bloque×subpoblación en `[0,0.5/C]`, con
`0.5/C` predeterminada. Una combinación controlada no factible se rechaza antes
del RNG propio de cuotas. Las
distribuciones de reproducción posteriores se deciden en sus fases y no
bloquean F5.

<a id="dtc-094"></a>

## DTC-094 — Modelo escalonado de deslizamiento como servicio reutilizable

El simulador tiene **un solo** modelo de deslizamiento de la polimerasa sobre
un microsatélite, y se declara aquí una vez para los dos consumidores que lo
usan: el deslizamiento **biológico**, que decide la verdad, y el **tartamudeo**,
que ensucia las lecturas. Que el mecanismo físico sea el mismo no los convierte
en el mismo componente, y DTC-096 fija por qué no pueden compartir azar.

**Cada evento añade o elimina exactamente una unidad completa** del motivo,
nunca una fracción ni un salto de varias unidades. El modelo se identifica por
el literal `stepwise` y el tamaño de paso por el literal `1`; cualquier otro
modelo y cualquier otro tamaño de paso **se rechazan antes de simular**,
nombrando el parámetro y el valor rechazado. El salto de varias unidades queda
fuera por DTF-015 y no se reabre aquí.

La unidad se gana o se pierde **en la tirada ininterrumpida más larga** del
tracto, con desempate por la izquierda, que es la misma regla que fija
`TandemRepeatTract::apply_unit_step` y la misma que aplica el tartamudeo. No se
escribe dos veces: es una sola operación del dominio con dos consumidores.

**Regla del borde.** Un paso que dejaría **menos de una unidad completa en la
tirada más larga** no es viable, y tampoco lo es uno que rebasaría el **tope de
unidades** del locus. Cuando falta alguno de los dos pasos admitidos, el linaje
ha llegado a un **borde** y **se detiene ahí**: el evento **no se aplica**, **no
consume aleatoriedad** y **se cuenta aparte**, y los eventos que le quedaban
tampoco se aplican. No es un error: un microsatélite que ha llegado a su mínimo
es un resultado biológico legítimo.

**El borde no empuja, y esto se decidió con las dos opciones delante.** La
primera implementación aplicaba el paso que quedaba, de modo que un linaje que
tocaba el tope se iba contrayendo paso a paso, todos sus eventos se contaban
como aplicados, ninguno como no aplicado, y el documento publicaba un sesgo
intermedio cuando no había habido ni una expansión: **el borde decidía la
biología en lugar del sesgo**. El responsable lo resolvió el 27 de septiembre de
2026 y descartó expresamente dejar que el borde empujara. Rige **en el tope y en
el suelo**, con el signo cambiado.

*Cuidado con el fixture.* Con un sesgo de exactamente 0 o exactamente 1 las dos
conductas son **indistinguibles**, porque el paso que el borde deja fuera es
justo el que el sesgo ya anulaba. Toda comprobación de esta regla necesita un
**sesgo intermedio**.

El recuento de eventos no aplicados **se publica**. Sin él, un servicio que
ignorase en silencio los eventos sobrantes y otro que los contase producirían el
mismo tracto final y ninguna comprobación podría distinguirlos.

El ensayo que elige entre expandir y contraer es `bernoulli-v1` sobre el sesgo
de expansión declarado, y **las probabilidades exactas 0 y 1 no consumen
aleatoriedad**: un sesgo extremo no puede desplazar la secuencia de azar de
nada.

<a id="dtc-095"></a>

## DTC-095 — Generador del abanico de alelos por deslizamiento biológico

Es **condición de cierre** de la simulación de microsatélites: mientras no
exista, los alelos de un locus y sus frecuencias los escribe el usuario a mano.
El algoritmo lo adoptó el responsable el 20 de septiembre de 2026 y es el
siguiente, sin variantes.

Se simulan **tantos linajes independientes como copias haplotípicas** tenga la
población: una copia, un linaje. Cada linaje parte del **alelo que el genoma de
referencia lleva en el locus** y recibe un **número de pasos sorteado de una
Poisson de media declarada**, con `poisson-count-v1`, que es el único recuento
normativo del proyecto (DTC-067). Cada paso lo da el servicio escalonado de
DTC-094 con el sesgo declarado. **La frecuencia de cada alelo es la proporción
de linajes que terminan en él**, y las frecuencias se auditan con la suma
compensada de DTC-049 en vez de darse por uno.

El alelo de la referencia entra en el abanico **si algún linaje termina en él**,
y **su presencia no se fuerza**. Que desaparezca depende de los parámetros que
se den, no del modelo: con una media de pasos baja la mayoría de los linajes no
se mueven y el alelo está; con una media alta puede perderse, que es lo que
ocurre en una población real.

Es **una sola tirada de inicialización de la diversidad**, no una simulación de
generaciones. El modelo escalonado aplicado generación tras generación pertenece
a la fase de evolución y este trabajo no la invade.

**Alcance.** El generador solo admite un locus de **un motivo**: sin fronteras de
bloque declaradas no hay forma de partir el tracto de la referencia en bloques,
y un microsatélite compuesto se rechaza nombrando el locus. El tramo de la
referencia tiene que contener **al menos una copia completa** del motivo, o no
hay alelo ancestral del que partir. Y el escenario tiene que declarar **al menos
una copia haplotípica**, porque sin linajes no hay abanico.

**El tramo de la referencia tiene que caber bajo el tope del locus**, y se
rechaza si no cabe, nombrando el parámetro y el valor. Mirando el tope solo al
expandir se convertía en un «no crezcas más» y dejaba pasar un tramo que ya lo
rebasaba: el mismo alelo escrito a mano se rechazaba desde la entrega 1 y
generado se admitía. Peor aún, desde encima del tope la expansión nunca es
viable, de modo que el generador quedaba **unidireccional** y el sesgo declarado
dejaba de significar nada.

**La media de pasos tiene techo declarado**, mil por linaje, en el código y en el
esquema publicado. `poisson-count-v1` admite hasta `2^52` como primitiva, pero
aquí el parámetro lo escribe el usuario: sin techo se va de las manos el coste
—cada paso clona el tracto— y el sentido, porque un tracto tiene a lo sumo
`max_units` unidades y un linaje se detiene al tocar un borde, de modo que por
encima de unos pocos cientos de pasos la media deja de cambiar el resultado y
solo cambia lo que se tarda. Mil está holgadamente por encima de cualquier
tracto que este simulador admita: **solo rechaza, y no cambia el resultado de
ningún escenario aceptado**.

**Excluyente con los alelos declarados.** Un locus o declara sus alelos con sus
frecuencias, o declara el generador. Las dos cosas a la vez dejarían sin
respuesta la pregunta que el manifiesto contesta por cada alelo: de dónde salió.

**Tres limitaciones declaradas, que son decisiones biológicas y no detalles.**
Viajan en la documentación del generador **y en cada documento publicado**,
palabra por palabra:

1. **No hay ancestría compartida.** Cada copia deriva por su cuenta desde el
   ancestro, de modo que dos individuos emparentados no comparten historia en
   este locus. El espectro de frecuencias resultante es **más plano** que el de
   una población real, donde la deriva concentra unos pocos alelos. Importa si
   el conjunto se usa para probar herramientas de estructura poblacional.
2. **El número de pasos no es un proceso en el tiempo.** Cada copia acumula un
   número aleatorio de mutaciones sin decir en cuántas generaciones. Es coherente
   con que esto sea una sola tirada de inicialización.
3. **El genoma de referencia puede llevar un alelo que ningún individuo tiene.**
   Es real y ocurre, pero conviene saberlo antes de encadenar un análisis que lo
   dé por imposible.

**Cuándo se ejecuta.** El alelo ancestral son bases de la referencia, de modo que
el generador corre **después de preparar la referencia y antes de publicar
nada**: un rechazo suyo no deja un solo byte escrito. Corre igual en `validate`
y en `run`, para que los dos vean lo mismo.

**El juez rehace el abanico y lo contrasta**, y esto es una obligación y no una
descripción. El juez que relee los documentos publicados **vuelve a recorrer los
linajes por su cuenta** —con su propio filtro de viabilidad, su propia regla del
borde, su propio recuento por alelo y su propio orden canónico— desde la
referencia, los parámetros y la semilla **publicados**, y compara alelo a alelo.
Reconstruir la configuración llamando al productor **no cuenta**: eso compara el
abanico consigo mismo. Se midió por qué importa: emparejando cada alelo con la
frecuencia del simétrico de la lista, el espectro se publicaba del revés y toda
la batería seguía en verde.

La versión del algoritmo es `biological-slippage-v1` y se publica con cada
abanico.

<a id="dtc-096"></a>

## DTC-096 — Ámbito de azar propio del deslizamiento biológico

El deslizamiento biológico y el tartamudeo comparten mecanismo físico y **no
pueden compartir ámbito de azar**, porque uno decide la **verdad** y el otro las
**lecturas**. Si los alelos verdaderos salieran del modelo de tartamudeo,
cambiar el número de ciclos de PCR cambiaría el genotipo real de los individuos:
un parámetro de laboratorio decidiendo la biología.

El ámbito del generador es `initial_repeat_slippage`, nace en la versión **10**
del catálogo de ámbitos —`catalogs/subseed-scopes-v10.tsv`, que conserva
literalmente las filas de la versión 9— y se identifica por el **evento**, que es
el nombre del propio ámbito, y por el **identificador del locus**, junto al
escenario y la referencia. **Ni la posición del tracto ni el orden en que los
loci se declararon entran.** La consecuencia observable es la que importa:
insertar un locus que ordene antes que los existentes **no cambia el abanico de
los previos**.

Es **un ámbito por locus**, no por copia: el abanico y sus frecuencias son una
propiedad de la población entera, y el reparto entre copias lo decide después
`initial_repeat_allele`, que sí es por copia.

**Lo que no demuestra la separación, y conviene dejarlo escrito porque parece
que sí.** La comprobación «cambiar los ciclos de PCR no cambia ni un alelo de la
verdad» **no caza un ámbito compartido**: se cumple sola. Compartir ámbito en
este proyecto significa derivar la misma subsemilla, y DTC-081 exige que cada
decisión arranque su propio generador desde su clave y que no haya consumo
acumulado, de modo que un consumidor no desplaza al otro; además la verdad se
genera **antes** que las lecturas. Lo que sí caza un ámbito compartido es
**comparar las dos subsemillas** para la misma terna de identidades. La
comprobación antigua se conserva como prueba de sistema —es buena contra un
generador que leyera la configuración del tartamudeo por error— y su
documentación dice que no es la prueba del ámbito.

<a id="dtc-097"></a>

## DTC-097 — Publicación del origen de cada alelo de repetición

La configuración efectiva y el manifiesto general publican, **por alelo y no por
locus**, su **origen**: `declared` cuando lo escribió el usuario y `generated`
cuando lo produjo el generador de DTC-095. Va por alelo porque un documento con
loci de los dos tipos tiene que poder distinguirlos sin volver a la
configuración, y porque quien encadene un análisis necesita saber si una
frecuencia la escribió una persona o la produjo un modelo. El origen del alelo y
el del locus tienen que concordar: un locus con generador no publica alelos
declarados, ni al revés.

Un locus generado publica además su bloque `slippage` con, exactamente: el
modelo y el tamaño de paso —literales fijos—, la media de pasos y el sesgo de
expansión declarados, la **versión del algoritmo**, las de los dos algoritmos de
azar que consume, el **nombre del ámbito de azar** y la **versión del catálogo de
ámbitos**, las **tres limitaciones** del modelo palabra por palabra, y la
contabilidad de la generación: linajes, pasos sorteados, pasos aplicados y
**pasos no aplicados al agotar el tracto**. Los pasos aplicados más los no
aplicados tienen que sumar exactamente los sorteados, y el lector de los
documentos publicados lo comprueba, porque sumar dos campos y compararlos con un
tercero no es expresable en JSON Schema.

Las frecuencias **declarada** y **realizada** conservan su significado de
siempre: la declarada es la que el modelo o la persona fijó, y la realizada la
proporción de copias que el sorteo dio a ese alelo. Para un alelo generado la
declarada es la proporción de linajes que terminaron en él, y el origen que
viaja a su lado deja claro que nadie la escribió.

La **cota de bases de tracto** se recomputa después de generar: sale de los
alelos, y los alelos acaban de nacer. Dejarla como estaba publicaría una cota
falsa y el juez independiente encontraría lo realizado por encima de lo acotado.
El tope declarado se vuelve a aplicar en ese momento, que sigue siendo antes de
publicar nada.

**La parte VCF de la cota se mide, no se estima**, decidido el 2 de octubre
de 2026. Era el tracto de la referencia más la longitud de cada alelo
distinto, y suponía un registro de repetición con un `ALT` por alelo. Desde que
las variantes que pisan ese registro se funden en él (DTC-064), el registro se
ensancha y gana `ALT` de maneras que ninguna fórmula sobre lo declarado acota
sin ser absurda: un primer término por indel se quedaba corto, otro más amplio
también —al alinear se recorre el flanco que continúa el motivo, y se
publicaron 246 bases frente a una cota de 185—, y el más amplio rechazaba con
una deleción de 2300 bases una configuración que antes se admitía.

Ahora la publicación del VCF, **después de construir las copias y antes de
escribir nada**, cuenta las bases de tracto exactas de los registros que va a
escribir —`REF` más cada `ALT` de los registros de repetición, por locus— y esa
cifra es la que la configuración efectiva y el manifiesto publican donde iba la
cota del VCF; el total es la cota de la verdad más lo medido. El juez de la
contabilidad exige que esa cifra sea la que el VCF publicado lleva, en total y
**locus a locus**: cada registro de repetición se atribuye al primer locus, por
posición, cuyo intervalo contiene, y cada locus tiene que publicar lo que
llevan sus registros —sin esto, atribuirlo todo al primer locus de la
secuencia conservaría la suma—. La cota de la verdad sigue siendo la de antes,
estimada sobre lo declarado.

**El tope declarado** se aplica en dos momentos, y no con la misma cifra:

1. **Al admitir la configuración, por precaución**, sobre la cota de la verdad
   sola. Es una cota superior: una ejecución que la supere **podría caber**,
   porque lo realizado suele quedar por debajo; se rechaza igual, porque es lo
   único que se sabe antes de materializar a nadie.
2. **Al publicar el VCF**, sobre lo que la ejecución va a publicar de verdad:
   la verdad **ya realizada**, que está construida en ese momento, más lo
   medido en el VCF. El diagnóstico es de configuración, nombra el tope, la
   cifra con sus dos sumandos y el locus que más aporta, no escribe el VCF y
   deja un paquete fallido coherente. Usar aquí la cota de la verdad en lugar
   de lo realizado rechazaría ejecuciones que caben y daría en el mensaje una
   cifra que nadie iba a publicar.

De ahí una asimetría que conviene saber: **una configuración que `validate`
admite puede ser rechazada por `run`** al publicar el VCF, porque lo medido no
se conoce hasta entonces.

Consecuencia, dicha en voz alta: **la cifra publicada puede cambiar** en un
escenario con microsatélites, haya o no indels directas. Sin indels, la
estimación contaba también los alelos declarados que nadie sacó y el de la
referencia, que no produce `ALT`; lo medido no los cuenta. Por ejemplo, un
`(CA)50` con alelos de 50 y 45 unidades publicaba 290 y ahora publica lo que el
VCF lleva, 190. **No cambia en todos**: en la muestra de la segunda revisión
adversarial, 109 de 474 escenarios con microsatélites publicaron la misma
cifra que antes. Lo admite la regla S-26, porque la versión 6 aún no se ha
publicado fuera.

## DTC-098 — La referencia derivada por inserción es la de la simulación

Decisión del responsable del 01/10/2026, que concreta REQ-REF-37 y la decisión
S-03 del plan de microsatélites.

Un locus sobre referencia **suministrada** puede declarar `insert_after = P`
(base uno, en coordenadas del genoma de partida, `1 ≤ P ≤ longitud − 1`) con
`reference_units` obligatorio y sin `start`/`end`; la referencia `[reference]`
declara entonces `inserted_reference_id`. El tracto de cada locus que inserta
es su motivo, o sus motivos en orden, repetido tantas veces como digan las
unidades enteras de cada bloque.

**La referencia final es la de toda la ejecución.** La reducción es un
subproducto de preparación y nunca se simula sobre ella; la inserción no: la
matriz de SNV generada, los individuos, la verdad, el VCF, el evento técnico y
el generador de alelos por deslizamiento reciben la referencia final y su
FASTA, y los sorteos que se siembran con el identificador de referencia usan el
de la final. La traducción de coordenadas ocurre **una sola vez**, al resolver
la configuración: una posición `x` de partida pasa a `x` más la suma de las
longitudes de los tractos insertados con punto menor que `x` en su secuencia,
y un intervalo se traduce por sus extremos y se rechaza si un punto de
inserción cae dentro. Así se traducen los intervalos de los loci con tracto
existente y las regiones del generador de la matriz.

Se rechazan, antes de simular y nombrando las piezas: dos inserciones en el
mismo punto, una inserción dentro del tracto de otro locus, `insert_after` con
referencia simulada, y junto a cualquier inserción los cambios directos (sus
inserciones y deleciones se declaran en coordenadas vigentes, que no tienen
traducción sin un segundo modelo de coordenadas), la matriz de SNV
suministrada o reutilizada (están atadas por huella al genoma de partida) y la
reducción de referencia (se haría sobre la suministrada, sin el tracto).

**Publicación.** Se publican la copia suministrada, como siempre; el FASTA
canónico final en `resources/references/<id>.fa`; y el mapa
`reference_coordinate_map` v2 en `resources/references/<id>.map.json`. La
configuración efectiva v6 añade la sección opcional `inserted_reference`, y en
cada locus, cuando hay inserción, el intervalo **final** en `start`/`end` y lo
declarado en `declared`; el manifiesto general v6 añade el enlace de la final a
`references` y la sección opcional `reference_insertion` con la referencia de
origen, los artefactos y los cambios aplicados (REQ-MAN-07). Sin inserciones no
aparece nada nuevo y los bytes publicados no cambian. No hay versión de esquema
nueva: las tres versiones están en curso (S-26) y el mapa v2 ya bastaba.

**Juez.** El validador de los recursos publicados rehace el genoma final desde
la copia suministrada publicada y lo declarado en la configuración efectiva,
con un empalme propio y sin llamar a `insert_tracts`; lo compara byte a byte
con el FASTA final; valida el mapa por su contenido, porque su esquema admite
huecos, intervalos invertidos e identificadores iguales; comprueba cada
intervalo publicado de los loci, en la efectiva y en el manifiesto; y traduce
por su cuenta las regiones del generador de la matriz desde las del TOML y
exige las mismas que recibió el generador. Después, los jueces que ya
trabajaban sobre las bases de la referencia —el abanico de deslizamiento, la
matriz de SNV y el evento técnico— reciben la final ya verificada. **La verdad y
el VCF no se releen**: se comprueba, con una pasada barata, que nombran la
referencia de la simulación —la cabecera de la verdad con su identificador, su
huella del contenido y la huella del FASTA; el identificador y la longitud de
cada fragmento que los lleva; `##reference` y `##contig` del VCF— y que la
sección de loci del manifiesto también. Lo mismo vale sin inserción, contra la
suministrada. Corregido el 01/10/2026 tras la revisión adversarial: la primera
redacción decía que la verdad y el VCF se juzgaban contra la final, y antes de
esa comprobación de nombres no los leía nadie.


## DTC-099 — Multiplicador configurable de la PCR en los homopolímeros

Decisión del responsable del 06/10/2026, que cierra DT-009 y concreta
REQ-CONF-141, REQ-MAN-97 y REQ-VAL-61.

Las dos fuentes de la calibración discrepan en cuánto multiplica la PCR el
tartamudeo de los homopolímeros: Gymrek (2016), tabla 1, mide 0,17 frente a
0,027, 6,30 veces; Fungtammasan y col. (2015) escriben en su texto que la
diferencia se estabiliza en unas 2,5 veces a partir de unos 11 pb. En
dinucleótidos coinciden, 3,45 frente a unas 3,6. **Manda Gymrek por omisión, y
el autor del escenario puede cambiarlo.**

**Qué se configura.** Un solo número, `homopolymer_pcr_fold` en
`[technical_event.pcr.stutter]`: el cociente entre la fracción con PCR y la sin
PCR de los homopolímeros. Es exactamente el cociente entre columnas que DTC-093
mantiene independiente de la longitud del tracto, de modo que la modulación
sigue multiplicando las dos columnas por el mismo factor y la composición sigue
cerrando exactamente: `1 − (1 − r_sin·m)(1 − r_pcr(m)) = k·r_sin·m`. Se aplica
sobre la **constante** sin PCR de la calibración, 0,027, nunca sobre la tasa de
base que el escenario haya declarado, por la misma razón por la que hoy se usan
las dos constantes: `r_pcr` es lo que añade la amplificación. La columna sin
PCR no cambia, de modo que **ninguna ejecución sin PCR cambia** con ningún
valor. Consecuencia que la descripción del campo dice: con una tasa de base de
0,05 y `k = 2,5`, la fracción con PCR de extremo a extremo no es 0,125 sino
`1 − (1 − 0,05)(1 − 0,0416) ≈ 0,0895`.

**Qué no se configura.** Las longitudes de motivo de dos a seis: solo el
homopolímero está en disputa (DTF-063).

**Dominio.** `1 ≤ k ≤ 37`, finito. Con `k = 1` la probabilidad por ciclo del
homopolímero vale cero y no se consume aleatoriedad; con `k < 1` se rechaza,
porque la derivación solo despeja un incremento no negativo y lo pinzaría a
cero sin avisar; 37 es el mayor entero con `37 × 0,027 < 1`. Si la constante sin
PCR cambiara, el tope se recalcula; una prueba lo comprueba. Declararlo junto a
una sustitución directa de la probabilidad por ciclo de los homopolímeros se
rechaza: serían dos valores para lo mismo. El dominio y la combinación se
comprueban **aunque no haya loci declarados**, donde el campo no tendría
efecto. Diagnóstico `GPS-CONF-008`, nombrando
`technical_event.pcr.stutter.homopolymer_pcr_fold` y el valor.

**Ausente no es 6,30.** Sin el campo se usa la constante con PCR, 0,17, tal
cual; el programa no multiplica nada y publica los mismos bytes que antes. El
campo publicado aparece solo cuando se declara.

**Publicación y juez.** La configuración efectiva v6 y el manifiesto de evento
v3 añaden, en el bloque del tartamudeo clonal, la propiedad opcional
`homopolymer_pcr_fold`, y la fila de los homopolímeros publica la fracción con
PCR resultante, con su probabilidad marcada `derived`. El algoritmo sigue
siendo `clonal-stutter-v2`, porque la derivación no cambia, y ningún catálogo
de ámbitos se mueve. No hay versión de esquema nueva: las versiones están en
curso (S-26). El juez de los recursos publicados ata las fracciones publicadas
a la calibración y al multiplicador (REQ-VAL-61) por dos comparaciones: el
bloque de la configuración efectiva tiene que ser el que resuelve el TOML
copiado —esa comparación ya existía y es la que lo ata—, y el del manifiesto de
evento tiene que ser idéntico al de la efectiva. La segunda es nueva: el
08/10/2026 se reprodujo que, sin ella, un manifiesto de evento con la fracción
con PCR del homopolímero y su probabilidad cambiadas a la vez, y con sus
huellas rehechas en el manifiesto general, pasaba el juez; el contenido del
bloque solo lo comprobaba el esquema. La sospecha original —que la efectiva aceptara esa
misma manipulación— **no** se reprodujo: la rechaza la comparación con el TOML.
**Estas dos comparaciones sustituyen a las guardas campo a campo** que proponía
el estudio previo —fracciones sin PCR y con PCR exactamente iguales a las
constantes, multiplicador en dominio y nunca junto a una fila declarada—, que
habrían repetido lo que la comparación con el TOML ya exige; el responsable lo
aprobó el 08/10/2026. La igualdad cubre **el bloque entero** del tartamudeo
clonal, no el resto de la sección de amplificación del manifiesto de evento:
sus ciclos fuera del bloque y los recuentos de linaje siguen sin atarse a nada
(DT-043).

**Consecuencia sobre una cifra de DTC-092.** El «entre un 14 % y un 28 %» que
perdería una preparación amplificada sin tartamudeo de base vale **con la
calibración por omisión**. Con `k = 2,5` el homopolímero perdería un 38 %, y con
`k = 1`, todo. La guarda que exige el tartamudeo de base no cambia.
