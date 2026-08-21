# Registro de decisiones técnicas del simulador

- **Versión:** 0.2
- **Estado:** F1 cerrado; decisiones F2 simplificadas y pendientes de revisión
- **Fecha:** 21 de agosto de 2026
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
| DTC-018 | La cobertura objetivo usa la longitud haploide de la referencia estructural del individuo y se distingue de las coberturas WGS básica, final del FASTQ y alineada por referencia. Pregunta 10 del informe adversarial (RAD-004 y RAD-017); REQ-WGS-04, REQ-WGS-11 y REQ-DAT-24 | Un mismo FASTQ no puede tener objetivos distintos por referencia de mapeo | Cada secuencia se suma una vez sin ploidía; haplotipos y referencias de mapeo no cambian el denominador; la cobertura alineada es solo una observación |
| DTC-019 | V1 no modela elementos transponibles como clase funcional. Pregunta 12 del informe adversarial (RAD-006); sección 3 y REQ-TE-01 a REQ-TE-14 retirados | Reducir el alcance y eliminar la contradicción entre referencias derivadas y eventos TE | No hay catálogos, familias, mecanismos ni metadatos TE; las copias e inserciones se clasifican mediante operaciones genéricas |
| DTC-020 | La generación por eventos aplica primero los explícitos y después los aleatorios en el orden duplicaciones, deleciones, inversiones, translocaciones, repeticiones en tándem, *indels* y SNP. Pregunta 13 del informe adversarial (RAD-008); REQ-IND-54 a REQ-IND-57 | Hacer único y comprobable el resultado ante colisiones y agotamiento del espacio válido | Cada evento se elige del espacio válido actualizado; no hay reintentos sobre inválidos ni límites de intentos, y la insuficiencia hace fallar la generación conservando los recuentos sorteado, aceptado y activo. El muestreo inicial desde una matriz SNV de DTC-052 no es otro evento de esta secuencia |
| DTC-021 | Cada conjunto solicitado produce exactamente un VCF multimuesta estándar; un caso de verdad usa uno o varios registros y, solo si pierde detalle, una única entrada complementaria. Pregunta 14 del informe adversarial (RAD-011 y RAD-023); REQ-VER-04 y REQ-VER-19 a REQ-VER-26 | Mantener interoperabilidad VCF sin perder la verdad de casos complejos | El VCF sigue siendo válido por sí mismo; la reconstrucción completa usa el par enlazado, BND conserva su cardinalidad estándar y referencias distintas permanecen independientes |
| DTC-022 | V1 usa exclusivamente Poisson para recuentos de variantes aleatorias obtenidos desde una media y el modelo escalonado simple de una unidad para mutaciones de repeticiones en tándem. Pregunta 17 del informe adversarial (RAD-020); REQ-CONF-25, REQ-CONF-61/62 y REQ-IND-18/35 | Cerrar el catálogo mínimo y hacer observables las configuraciones válidas | Eventos explícitos, número exacto de loci de una matriz y modo controlado no son distribuciones alternativas; el espectro alélico y el muestreo de genotipos tampoco son recuentos Poisson. La tasa y el sesgo expansión/contracción siguen configurables y cualquier otra distribución, modelo o salto se rechaza antes de simular |
| DTC-023 | V1 aplaza completamente la calibración empírica automática desde datos reales. Pregunta 18 del informe adversarial (RAD-024); REQ-CONF-85 a REQ-CONF-87 y REQ-VAL-49 | Evitar una capacidad amplia sin contrato verificable y mantener cerrado el alcance de perfiles | FASTQ y referencia de calibración, solicitudes de estimación y perfiles empíricos reutilizables se rechazan antes de simular; se conservan los perfiles incorporados de Illumina y PacBio HiFi y sus sustituciones explícitas |
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
| <a id="dtc-036"></a>DTC-036 | Todo documento estructurado propio declara su clase y `schema_version`: `config_kind` distingue configuración general y escenario, `manifest_kind` distingue los manifiestos, `effective-config.json` usa el discriminador equivalente de su esquema y la salida JSON de diagnósticos declara `diagnostic_kind`. La primera versión de esquema es la cadena exacta `"1"`. Resolución de DTPD-003; REQ-CONF-01/10/27/111 a REQ-CONF-115, REQ-MAN-02/69/73 y REQ-SAL-03/04 | Impedir que un documento se interprete con la clase o semántica de otra versión y separar la evolución del esquema de la versión del programa o del dialecto de validación | V1 solo lee versiones expresamente conocidas: una versión ausente, anterior, posterior o desconocida se rechaza antes de interpretar parcialmente el contenido, con clase, versión encontrada, versiones admitidas y ruta. No hay migración automática ni modificación silenciosa; el programa escribe siempre la versión actual. Versión del programa, `schema_version` y JSON Schema Draft 2020-12 son conceptos distintos. Un conversor solo podrá aprobarse cuando exista una segunda versión real: producirá una copia, conservará el original y registrará origen, destino y resultado; no se diseña ni implementa ahora. DTC-044 concreta la clase y el esquema del documento de diagnósticos |
| <a id="dtc-037"></a>DTC-037 | La aleatoriedad de v1 usa generadores ChaCha12 nombrados y locales, nunca `StdRng` ni un flujo global. La semilla raíz y cada subsemilla tienen 256 bits y se representan en TOML, configuración efectiva y manifiestos como exactamente 64 hexadecimales minúsculos. Si falta la raíz, se obtienen 32 bytes del CSPRNG del sistema operativo antes de simular y se registran valor y origen. Resolución de DTPD-004; REQ-CONF-26/27, REQ-CONF-101/111 a REQ-CONF-115 y REQ-EVO-12/62 | Aislar biología, muestreos y eventos para que el consumo de un proceso o la planificación no alteren decisiones ajenas y permitir paralelismo reproducible sin compartir estado mutable | BLAKE3 deriva cada subsemilla con el contexto ASCII exacto `genome_population_simulator/v1/subseed`. El material binario canónico contiene, en este orden, la raíz cruda de 32 bytes, la versión de derivación como `u32` big-endian, el número de componentes como `u32` big-endian y, para cada componente en el orden normativo del catálogo de ámbitos, nombre de campo ASCII y valor UTF-8, ambos precedidos por su longitud `u32` big-endian. La identidad de flujo nunca procede de posición accidental, iteración, hilo o planificación. Se separan al menos población inicial, generaciones y procesos biológicos, muestreos y cada evento técnico, afinando por unidades independientes. Cada consumidor solo avanza su generador. Se registran raíz, origen, ChaCha12, versión de derivación y subsemillas de evento; los ámbitos internos se reconstruyen desde IDs y reglas versionadas. Ninguna fila publicada de un catálogo, sus campos/orden ni sus vectores se reescribe: añadir ámbitos exige una nueva versión de catálogo y puede conservar la versión de derivación si codificación/contexto no cambian; alterar una fila existente exige además nueva versión de derivación. Cambiar ChaCha12, BLAKE3, esta codificación o un algoritmo de muestreo cambia la reproducibilidad. DTC-038 concreta la garantía; DTC-039 separa el dominio de IDs; DTC-043 rige dependencias; DTC-048 fija el lock; DTC-040 reserva SHA-256 para integridad publicada |
| <a id="dtc-038"></a>DTC-038 | V1 declara por separado cuatro niveles de reproducibilidad: biología; decisiones y lecturas técnicas; equivalencia semántica de salidas; e identidad binaria. Resolución de DTPD-005; REQ-CONF-26/27 y REQ-CONF-111 a REQ-CONF-115, REQ-REP-22/23 y REQ-MAN-69 | Permitir repetir y comparar cada clase de resultado sin deducir una garantía más fuerte desde otra más débil | La biología se garantiza con la misma versión del simulador, plataforma inicialmente admitida, semilla raíz, configuración efectiva, contenidos de entrada, perfiles y esquema y algoritmos aleatorios versionados. El nivel técnico exige además la misma subsemilla de evento y versiones y condiciones de las dependencias participantes. La equivalencia semántica compara contenido lógico y, ante herramientas externas, exige sus mismas versiones, parámetros y condiciones verificadas. V1 no garantiza identidad binaria por defecto para ninguna clase; solo podrá declararse por formato y condiciones tras una prueba expresa. No se promete reproducibilidad entre versiones del simulador ni entre plataformas no verificadas. El manifiesto registra por nivel `garantizado`, `no garantizado` o `no aplicable`, condiciones, versiones, perfiles, herramientas y dependencias. Toda semilla generada se registra y puede reutilizarse. Se mantienen comparadores separados para verdad lógica, contenido técnico, semántica y bytes. DTC-042 fija la plataforma oficial inicial; crates, toolchain y formatos adicionales permanecen separados |
| <a id="dtc-039"></a>DTC-039 | Cada entidad separa un identificador estable de una etiqueta humana Unicode opcional y no única. Los identificadores aportados son ASCII, sensibles a mayúsculas y cumplen conceptualmente `[A-Za-z0-9][A-Za-z0-9._-]{0,63}`; no admiten espacios, rutas ni otros caracteres. En el dominio se usan tipos distintos como `IndividualId`, `EventId` o `ReferenceId`, nunca `String` genérico en una frontera. Resolución de DTPD-006; REQ-CONF-01/38/39/51/52/101, REQ-IND-25/26, REQ-RG-03/04/07, REQ-DAT-68 a REQ-DAT-71 y REQ-VER-10/12 | Evitar confundir categorías, etiquetas, rutas o ámbitos y mantener identidades reproducibles compatibles con formatos genómicos | Los IDs automáticos se derivan con BLAKE3 y el contexto ASCII exacto `genome_population_simulator/v1/id`, usando el material canónico de DTC-037 y el catálogo normativo de campos. Toman los primeros 128 bits del resultado, los codifican en Base32 RFC 4648 mayúscula sin relleno `=` y anteponen el prefijo ASCII de clase del catálogo. El ID de ejecución incluye como causas mínimas raíz y digest de configuración; los demás incluyen tipo e IDs causales estables según catálogo, nunca hilo, ruta, orden físico o reloj. Los catálogos de campos y prefijos y sus vectores son artefactos normativos de F1 y no cambian sin nueva versión de derivación. Toda colisión se detecta y hace fallar. Se respetan los ámbitos normativos: variantes locales al conjunto/VCF, `RG` único en la ejecución, lectura o plantilla única en el evento y mates con el mismo ID de plantilla y roles distintos; el ID de individuo se conserva en `SM`, VCF y manifiesto. Rutas y nombres físicos se derivan o escapan sin adquirir autoridad. El manifiesto relaciona ID opaco, clase, ámbito, causa y etiqueta. Los IDs automáticos permanecen compactos y válidos para SAM, VCF y FASTQ; no se fija aquí el layout de directorios |
| <a id="dtc-040"></a>DTC-040 | SHA-256 es la única suma de comprobación publicada por v1. Se identifica como `sha256` y su digest se representa mediante 64 caracteres hexadecimales minúsculos calculados sobre todos los bytes exactos del fichero, sin normalización. Resolución de DTPD-007; REQ-REF-17/18/21, REQ-MAN-39/76 y REQ-SAL-19 | Proporcionar una identidad de contenido interoperable y verificable para cada entrada, salida y fichero enlazado | Cada manifiesto registra algoritmo, digest, tamaño en bytes y ruta relativa de cada fichero. La suma se calcula o verifica al copiar entradas, al cerrar salidas, antes de exportar y después de exportar; copia y original deben coincidir. Puede calcularse incrementalmente durante la E/S, pero las comprobaciones exigidas no se sustituyen por tamaño, fecha, metadatos ni muestreo parcial. Cada fichero conserva su digest; una suma adicional del paquete no reemplaza las individuales. MD5 solo se emite cuando un estándar genómico lo exige, separado y nunca como integridad principal. BLAKE3 queda reservado a derivación de subsemillas e IDs con dominios distintos. DTC-045 concreta la publicación; la futura dependencia SHA-256 debe superar DTC-043 y quedar fijada por el lock de DTC-048 |
| <a id="dtc-041"></a>DTC-041 | V1 admite FASTA de texto sin comprimir y encapsula `noodles-fasta` detrás de un contrato propio como adaptador Rust puro para lectura y escritura por flujo, FAI y acceso indexado. El dominio valida al menos una secuencia no vacía, primer token de cabecera no vacío y único y alfabeto A/C/G/T/N en cualquier caja, sin huecos ni otros símbolos IUPAC. Resolución de DTPD-008; REQ-REF-01 a REQ-REF-18, REQ-REF-21/23, REQ-VAL-45, REQ-MAN-86 y REQ-SAL-01 | Mantener el contrato de referencia independiente de una biblioteca concreta y admitir ficheros grandes sin perder el original ni relajar las reglas biológicas | Un FASTA suministrado se copia byte a byte sin reescribirlo, se valida y calcula su SHA-256 durante el flujo, deja intacto el original y solo se publica tras éxito. Las bases se normalizan a mayúsculas únicamente en la representación interna; se conservan copia exacta y cabecera cruda, con primer token como ID y el resto como descripción. Una descripción no UTF-8 se conserva sin pérdida mediante bytes codificados reversiblemente y una marca de codificación en salida y manifiesto. Los FASTA simulados o derivados se escriben canónicamente sin compresión, con A/C/G/T/N mayúsculas, LF, ancho fijo y orden estable; el ancho se versiona si afecta a bytes. Los índices solo se crean cuando se necesitan. `samtools faidx` es validador independiente de pruebas. La compresión FASTA se diagnostica como no admitida. Evaluación final de F1 a 2026-08-10: `noodles-fasta` 0.66.0, MIT, MSRV 1.89.0 compatible con Rust 1.97.1 y edición 2024 queda fijada con características mínimas síncronas y sin `async`. Su grafo no compila ni enlaza código nativo; los cuatro `build.rs` transitivos de Crossbeam y Rayon fueron auditados y solo emiten configuración para Cargo, sin ejecutar procesos, usar la red, escribir fuentes ni generar o compilar código. La dependencia transitiva de BGZF incorpora compresión y paralelismo aunque F1 use FASTA sin comprimir; el coste se acepta y se vigila. La API upstream evoluciona con frecuencia; adaptador propio, versión exacta, lock y pruebas mitigan el riesgo conforme a DTC-043 y DTC-048 |
| <a id="dtc-042"></a>DTC-042 | La plataforma nativa oficialmente validada para v1 es Ubuntu 24.04 LTS x86-64 GNU, con triple `x86_64-unknown-linux-gnu`, configuración de CPU portable y sin `target-cpu=native`. Otras distribuciones x86-64 GNU pueden funcionar, pero no reciben garantía hasta incorporarse a una matriz. Ubuntu 24.04 sobre WSL2, siempre en el sistema de ficheros Linux/ext4 y no en `/mnt/c`, es el entorno de desarrollo admitido y actual, no una plataforma de publicación independiente. Resolución de DTPD-009; REQ-CONF-111 a REQ-CONF-115, REQ-REP-22/23, REQ-MAN-05 y REQ-SAL-19; hoja de ruta F1.1 y F9.7 | Disponer de un entorno rápido de desarrollo y una referencia de publicación verificable sin convertir el triple GNU en una promesa genérica ni heredar particularidades de WSL2 | Las pruebas ordinarias se ejecutan en WSL2 sin dependencias exclusivas. Antes de publicar se valida en CI Ubuntu 24.04 nativo y en un entorno OCI cuya base y receta se fijan por digest; se prueban atomicidad, publicación no-clobber en su sistema de ficheros, señales, permisos, recursos y procesos externos. V1 distribuye un `.tar.gz` versionado con ejecutable, esquemas, licencias y ejemplos, cuyo SHA-256 se calcula sobre los bytes exactos del archivo, y una imagen OCI de referencia identificada por el digest `sha256:` de su descriptor de índice o manifiesto, conservando y verificando los descriptores de configuración y capas. Un archivo que exporte la imagen OCI tiene además su propio SHA-256 independiente. Se registran plataforma, arquitectura, versión del simulador, referencia y digest de imagen y herramientas externas. macOS, Windows nativo y ARM64 no tienen soporte oficial; ampliarlos exige dependencias y matriz de pruebas. No se promete identidad binaria entre compilaciones. DTC-043 conserva dependencias, DTC-045 la publicación y DTC-048 el toolchain |
| <a id="dtc-043"></a>DTC-043 | V1 prefiere crates Rust puros para configuración, serialización, hashes, RNG y adaptadores de formatos cuando satisfagan correctamente el estándar. No enlaza bibliotecas nativas C/C++ salvo necesidad imprescindible evaluada por separado y aprobada expresamente. Las herramientas genómicas completas y consolidadas, como minimap2, se ejecutan como procesos externos detrás de adaptadores propios. Resolución de DTPD-010; REQ-CONF-69, REQ-DAT-21/22, REQ-MAN-34 a REQ-MAN-36 y REQ-VAL-11; hoja de ruta F1.1 | Reducir la superficie de distribución y acoplamiento sin reimplementar herramientas especializadas ni filtrar decisiones accidentales de una dependencia al dominio | Toda dependencia queda detrás de un contrato propio con tipos y errores del sistema. El adaptador de proceso controla versión, parámetros, comando, entradas y salidas, código de salida, stdout, stderr y fallos. Antes de admitir una dependencia se documentan propósito y requisito, corrección del estándar, mantenimiento, seguridad, licencia, rendimiento, estabilidad de API, características activadas, código nativo, scripts de construcción y coste de distribución. La mera presencia de `build.rs` no implica rechazo: cada script se audita por su función, entradas, salidas y efectos. Se rechaza si ejecuta procesos externos o accede a la red, modifica fuentes, genera o compila código nativo, o produce efectos no necesarios o no acotados, salvo decisión expresa que justifique y controle la excepción. Se minimizan dependencias y *features* y no se añaden por conveniencia anticipada. Las entregas fijan versiones y el manifiesto registra dependencias y herramientas efectivas; DTC-048 fija el toolchain y el lock de los crates aceptados. `noodles-fasta` permanece aprobado por DTC-041, sin prejuzgar otros crates. DTC-047 conserva la validación independiente y DTPD-105, la concreción de su adaptador de proceso externo |
| <a id="dtc-044"></a>DTC-044 | V1 usa errores internos tipados y jerárquicos, separados de una capa pública estable de diagnósticos. Sus categorías mínimas son configuración o esquema, entrada, formato, validación de dominio o biología, recursos, proceso o herramienta externa, E/S o publicación y defecto interno. Resolución de DTPD-011; REQ-CONF-09/14, REQ-FMT-04/07 y REQ-VAL-01 a REQ-VAL-51; hoja de ruta F1.2 | Ofrecer rechazos automatizables y mensajes accionables sin acoplar el contrato público a textos o tipos accidentales de dependencias | Todo diagnóstico público contiene código estable del catálogo versionado `GPS-{categoría}-NNN`, donde la categoría es `CONF`, `INPUT`, `FMT`, `DOMAIN`, `RES`, `TOOL`, `IO` o `INTERNAL`; además contiene categoría, mensaje accionable y contexto disponible: ruta y ubicación, campo e IDs de entidad, evento o artefacto. El catálogo debe enumerar los códigos de F1 antes de implementarla y ningún código se reutiliza ni cambia de significado. La cadena causal técnica queda para detalle o log, no como mensaje principal. Los fallos esperables retornan sin `panic`; los defectos internos se distinguen. Por defecto se renderiza texto humano a stderr. El modo JSON emite exactamente un documento RFC 8259 por invocación, nunca JSONL, validado por esquema y con `diagnostic_kind: "gps_diagnostics"`, `schema_version: "1"` y la lista `diagnostics`. Salidas: 0 éxito o solo avisos; 1 fallo de ejecución o dominio con configuración válida; 2 uso, configuración o entrada inválida; 3 recursos; 4 herramienta externa; 5 E/S o publicación; 70 defecto interno. Los contratos públicos no filtran tipos de crates; DTC-036 rige clase y versión y DTC-048 fija dependencias cuando se elijan |
| <a id="dtc-045"></a>DTC-045 | Toda ejecución iniciada de v1 trabaja en un directorio provisional hermano del destino, en el mismo sistema de ficheros y marcado `INCOMPLETE`. Cada artefacto se cierra, valida e integra primero desde su propio temporal. Un directorio final de éxito o fallo controlado solo aparece cuando sus manifiestos y catálogo de estados son coherentes y se publica mediante una operación Linux atómica no-clobber equivalente a `renameat2(..., RENAME_NOREPLACE)`, detrás de un adaptador propio. Resolución de DTPD-012; REQ-VAL-11, REQ-MAN-73/76 y REQ-SAL-05/19/21 a REQ-SAL-26 | Impedir que una ejecución parcial parezca completa o sobrescriba otra concurrente sin prometer durabilidad física | Se distinguen tres resultados. Un rechazo de admisión o pre-simulación no publica destino: solo deja provisional o diagnóstico externo conforme a REQ-SAL-26. Tras comenzar, un fallo controlado retira o ignora temporales incompletos, conserva y declara exportables los artefactos o conjuntos independientes ya validados, cierra y valida el manifiesto general con estado `fallido`, diagnóstico y catálogo de unidades, y publica atómicamente el directorio final como fallido, nunca completo. Una interrupción no controlable, como SIGKILL o caída, deja el provisional `INCOMPLETE`; conserva diagnóstico solo cuando fue posible escribirlo. Si el destino existe, la operación devuelve `EEXIST` y se rechaza sin sobrescribir, mezclar ni reparar; si el sistema de ficheros no soporta la garantía, se rechaza de forma identificable antes o durante publicación, sin fallback check-then-rename. V1 no reanuda ni reutiliza parciales. Al exportar a otro sistema de ficheros se copia a un provisional allí, se verifican tamaño y digest y se publica con la misma garantía no-clobber; no hay atomicidad entre sistemas. Se prueban concurrencia e interrupciones en WSL2/ext4, Ubuntu nativo y el sistema de ficheros OCI admitido según DTC-042. Publicación lógica y `fsync` o durabilidad física son garantías distintas |
| <a id="dtc-046"></a>DTC-046 | La línea base inicial usa tres escenarios sintéticos del mismo recorrido: A, aceptación de 100 kbp; B, crecimiento con referencias de 1, 10 y 100 Mb; y C, observación de 1 Gb manual o programada. **DTC-078 limita F1 a ejecutar A/B/C en el entorno local disponible; la comparación nativa de publicación y OCI pertenece a F9.** Resolución de DTPD-013; REQ-NOF-01 a REQ-NOF-08; hoja de ruta F0.4, F1 y F9 | Obtener evidencia comparable de corrección y crecimiento sin convertir medidas iniciales en promesas de producto | Se miden tiempo total y por etapa, CPU cuando sea fiable, RSS máximo, bytes de entrada, temporales y salida, MB/s y número de ficheros o artefactos. Se usa una construcción optimizada y portable, sin `target-cpu=native`, configuración y semilla fijas, salida limpia y repeticiones suficientes. F1 registra el entorno local efectivo y conserva la línea base; F9 repite la aceptación de publicación en Ubuntu 24.04 nativo y OCI. No hay inicialmente objetivos de producto ni fallos por pequeñas variaciones temporales: corrección e integridad siguen siendo obligatorias. Umbrales y presupuestos permanecen en DTPD-115. Estos tamaños son cargas sintéticas, no una promesa de escala biológica futura |
| <a id="dtc-047"></a>DTC-047 | Los validadores independientes de v1 son dependencias de pruebas, CI y aceptación de publicaciones; no se ejecutan sobre cada fichero de una ejecución normal ni requieren instalación por la persona usuaria. El runtime valida internamente los artefactos e invariantes que corresponden a su fase. **DTC-078 asigna a F1 únicamente los oráculos de TOML, FASTA, SHA-256 y JSON necesarios para sus salidas; listado/extracción de `.tar.gz` e inspección/arranque OCI se ejecutan en F9.** Resolución de DTPD-014; REQ-FMT-04/07, REQ-VAL-01 a REQ-VAL-51 y REQ-SAL-05; hoja de ruta F1 y F9.5 | Obtener un oráculo distinto del productor sin convertir herramientas de prueba en dependencias operativas ni duplicar todo el coste de validación en cada ejecución | FASTA se contrasta con `samtools faidx`, SHA-256 con `sha256sum` y documentos JSON con una implementación JSON Schema distinta. En F9, `.tar.gz` se contrasta mediante listado y extracción estándar y la imagen OCI mediante inspección y arranque de prueba. En sus fases, VCF se contrasta con bcftools y BAM/CRAM con samtools u otra referencia independiente. Los casos incluyen válidos, truncados, corruptos, campos desconocidos y errores de borde; CI fija y registra versiones y plataforma. El mismo escritor o lector nunca es el único oráculo. Validez formal y corrección biológica se prueban por separado. Una auditoría externa opcional en runtime queda fuera de v1. DTC-048 rige la fijación de versiones y las DTPD de formatos posteriores, sus validadores concretos |
| <a id="dtc-048"></a>DTC-048 | El toolchain oficial inicial es exactamente Rust 1.97.1 stable, verificado como estable actual el 16 de julio de 2026. Al crear F1, `rust-toolchain.toml` fija `channel = "1.97.1"`, `profile = "minimal"`, los componentes `rustfmt` y `clippy` y el target `x86_64-unknown-linux-gnu`; la edición `2024` se declara en cada `Cargo.toml` o en el workspace, no en el fichero de toolchain. `Cargo.lock` se versiona para la aplicación y CI y entregas usan `--locked`. Resolución de DTPD-015; REQ-CONF-111 a REQ-CONF-115, REQ-MAN-05 y REQ-REP-22/23; hoja de ruta F1.1 | Hacer reproducible y auditable la construcción sin declarar una compatibilidad mínima que no se prueba ni trasladar Rust a las personas usuarias del binario o imagen | No existe otra MSRV inicial: el compilador admitido es el fijado. Manifiesto y metadatos registran rustc y Cargo, edición, lock o resolución efectiva y herramientas, sin prometer identidad binaria. Las versiones concretas de crates se eligen al implementar F1 después de comprobar compatibilidad y DTC-043 y quedan fijadas por el lock; una candidata evaluada no cuenta aún como fijada. Toda actualización de Rust o dependencias es deliberada y aislada y exige pruebas, validadores independientes, vectores de RNG y reproducibilidad y escenarios de rendimiento; después se vuelven a declarar garantías sin asumir equivalencia. Rust 1.97.1 sustituye la propuesta accidental basada en Rust 1.94.0 local. Las personas usuarias de binario o imagen no necesitan Rust |
| <a id="dtc-049"></a>DTC-049 | Los valores numéricos configurables de probabilidad o frecuencia y sus tolerancias generales se representan como `f64`. Todo valor y tolerancia debe ser finito; una probabilidad o frecuencia pertenece al intervalo cerrado `[0,1]`. La tolerancia absoluta general es configurable en el intervalo cerrado `[0, 1e-6]` y vale `1e-9` por defecto. Resolución parcial de DTPD-110; REQ-CONF-106/107, REQ-VAL-40, REQ-MAN-81 y ACE-105 | Obtener una validación y normalización únicas, auditables y reproducibles sin convertir el orden físico o el paralelismo en parte accidental del resultado | Los valores se recorren en orden canónico estable por identidad y categoría y se suman con Neumaier, sin reducción paralela. Se calcula `delta = abs(suma - 1)` y el borde inclusivo acepta `delta <= tolerancia`. Si se acepta y la suma es positiva, finita y normalizable, cada valor se divide por la suma y se recomputan suma y delta con el mismo algoritmo y orden; si excede la tolerancia, la suma no es positiva o finita, la normalización no produce valores válidos o una entrada es inválida, se rechaza antes de simular. Configuración efectiva y manifiesto registran tolerancia efectiva, orden canónico, algoritmo y versión, suma y delta, y valores antes y después. Toda serialización decimal debe permitir recuperar exactamente el mismo `f64`; las garantías de reproducibilidad exigen conservar estas condiciones. Es una validación determinista de configuración: no altera DTC-024, no sustituye `0,5/C` ni ninguna tolerancia específica y no selecciona crate |

| <a id="dtc-050"></a>DTC-050 | V1 incorpora en F1.6 una referencia reducida derivada de un FASTA real ya admitido por F1.5. Usa exactamente un objetivo por bases o porcentaje, selección reproducible sin reemplazo ni solapamiento, orden relativo/orientación conservados, como máximo un residual, `N` literal y publicación conjunta de FASTA canónico, mapa completo, SHA-256 y manifiesto. REQ-CONF-124 a REQ-CONF-127, REQ-REF-24 a REQ-REF-36, REQ-VAL-52/53, REQ-MAN-89/90 y REQ-SAL-27 | Decisión funcional explícita del 10 de agosto de 2026, completada el 14 de agosto de 2026 | `reduction_percent` expresa el porcentaje eliminado, pertenece a `(0,100)` y retiene `floor(S × (100 − p) / 100)`. El parser conserva el literal decimal TOML y efectúa el cálculo racional con aritmética comprobada; convertir antes a `f64` está prohibido. Todo objetivo debe ser menor que `S`. Cromosomas y *scaffolds* son elegibles, conservan ID y clase y nunca se mezclan. DTC-056/DTC-057 resuelven coordenadas y selección. La reducción no es captura dirigida ni variante estructural biológica |
| <a id="dtc-051"></a>DTC-051 | Toda capacidad configurable de v1 debe estar disponible en Marimo además de TOML/CLI. La aplicación invoca `validate`/`run` sin shell; Rust, esquemas, `EffectiveConfig` y manifiestos conservan autoridad. F2 usa un TOML general y cada fase añade los controles de sus capacidades reales; los escenarios aparecen con su primer consumidor. REQ-UI-01 a REQ-UI-15 y REQ-VAL-54 a REQ-VAL-56 | Decisión funcional explícita del 10 de agosto de 2026, simplificada por DTC-079 | No autoriza Python dentro del núcleo, ejecución automática, reglas duplicadas, paquetes de resultados ni carga completa de artefactos grandes. DTC-079 concreta la interfaz local proporcional y F9 conserva distribución/endurecimiento. |
| <a id="dtc-052"></a>DTC-052 | Los SNV heredables iniciales de F2 proceden de una matriz de frecuencias común a la población, suministrada o generada desde número/densidad, espacio, sustituciones y un espectro alélico explícito. Los loci son marginalmente independientes; Hardy-Weinberg es el modelo predeterminado y `F` solo se admite en diploides. REQ-CONF-128 a REQ-CONF-130, REQ-IND-60 a REQ-IND-65, REQ-VAL-55, REQ-MAN-91 y REQ-UI-14 | Decisión funcional explícita del 10 de agosto de 2026, ampliada tras revisar espectros empíricos | Separa creación de la matriz y muestreo de genotipos; no usa densidades SNP independientes por individuo, no inventa LD y no convierte 1000 Genomes o ExAC en un perfil universal. En F2, HWE cubre P=1..10 y F solo P=2; una distribución genotípica explícita no tiene esquema/controles y cierra REQ-CONF-130 en F6.2. DTC-054 resuelve los formatos, algoritmos y perfiles de DTPD-118 |
| <a id="dtc-053"></a>DTC-053 | Una fase posterior admite un panel de bloques de LD locales para crear individuos. Su paquete canónico contiene exactamente `panel_blocks.tsv`, `panel_sites.tsv`, `panel_haplotype_alleles.tsv`, `panel_frequencies.tsv` y `panel.metadata.json`. El JSON declara clase/versión/`panel_id`, referencia y SHA-256, coordenadas base 1 cerradas y clase/ruta/tamaño/SHA-256 de cada TSV, más procedencia/herramienta cuando exista, sin duplicar tablas. Cada bloque declara `block_id`/`reference_id`/`sequence_id`/`start`/`end`; cada sitio ordenado, `block_id`/`site_id`/posición/`REF`/`ALT`. Los intervalos no se solapan y cada haplotipo se identifica por `block_id` + `haplotype_id`. V1 solo admite SNP suministrados bialélicos o multialélicos con `REF`/`ALT` de una base A/C/G/T. `panel_haplotype_alleles.tsv` es densa: una asignación por sitio y haplotipo, 0=`REF`, 1..N=`ALT`, sin ausencias, duplicados ni alelos fuera de catálogo. `panel_frequencies.tsv` es densa por población, bloque y haplotipo, incluidos 0/1, sin frecuencia global ni herencia entre poblaciones; sus sumas siguen DTC-049. Dentro del bloque se hereda la combinación local; fuera de sus intervalos, los SNP proceden de la matriz `alt_af`. Ningún locus de la matriz puede caer dentro de un intervalo de bloque: el conflicto detiene antes de crear individuos, identifica bloque/intervalo y locus/fila y no admite prioridad ni sobrescritura. El modo intacto conserva cada combinación intrabloque; la alternativa recombinante forma copias fundadoras directamente por individuo y bloque, reutiliza el motor común y conserva procedencia por segmentos. REQ-EVO-74 a REQ-EVO-77, REQ-VAL-56, REQ-MAN-92 y REQ-UI-15 | Decisión funcional explícita del 10 de agosto de 2026, precisada el 12–13 de agosto de 2026 | F6.4 incorpora paquete y ambas alternativas de inicialización; F7.1 reutiliza el motor común para meiosis después de resolver DTPD-111. La alternativa fundadora no es una meiosis histórica. No se infieren haplotipos desde frecuencias marginales, no se imputa `REF` o frecuencia 0 por ausencia ni se promete conservar LD después de recombinar. Otras clases de variante quedan en DTF-040 |
| <a id="dtc-054"></a>DTC-054 | DTPD-118 se resuelve con una matriz canónica `population_snvs.tsv` + `population_snvs.metadata.json`, ligada por SHA-256 a una referencia y con `alt_af` como única frecuencia. Una orden prepara entradas tabulares y otra genera el mismo recurso reutilizable. La generación produce candidatos bialélicos sin reemplazo; usa espacio uniforme o regiones no solapadas ponderadas, 12 sustituciones dirigidas y un espectro no solapado fijo/lineal/logarítmico en modo estocástico o exacto por mayores restos. REQ-CONF-128/129, REQ-IND-60 a REQ-IND-65, REQ-VAL-55, REQ-MAN-91, REQ-UI-14 y ACE-116 | Decisiones explícitas aprobadas el 12 de agosto de 2026 | Poblaciones y subsemillas son independientes; 0 y 1 son válidos, se eliminan sin reemplazo los `ALT` globalmente ausentes y se conservan los fijados; una suma ALT menor o igual que 1 conserva `REF=1-suma`, un exceso dentro de tolerancia se normaliza a 1 y otro mayor se rechaza. La multialelia aleatoria y conversión VCF se aplazan a DTF-039/038. `uniform_snv12_v1` es el valor de sustitución predeterminado; los perfiles humanos aproximados de DTC-068 son ejemplos explícitos, nunca implícitos ni tasas universales. Genotipos incompatibles y matriz vacía solo pueden continuar mediante autorizaciones específicas trazadas, no `--force` general. DTC-049 rige `f64`, Neumaier y tolerancia |
| <a id="dtc-055"></a>DTC-055 | DTPD-119 se resuelve con un paquete canónico de bloques haplotípicos locales, preparado por una única orden Rust, y dos construcciones por población×bloque: intacta o recombinante fundadora. REQ-EVO-74 a REQ-EVO-77, REQ-VAL-56, REQ-MAN-92, REQ-UI-15 y ACE-117 | Permitir LD local y recombinación fundadora reproducible sin confundir bloques con cromosomas, F con reproducción ni la fuente poblacional con sus resultados | Los cuatro TSV tienen esquema y orden exactos y una única referencia; v1 solo admite SNP A/C/G/T y fuentes del panel original. Intacto admite muestreo probabilístico/controlado y HWE, F diploide o distribución explícita. Recombinante admite solo muestreo probabilístico y HWE o F diploide, usa productores directos, el constructor común uniforme cerrado aquí, puntos x\|x+1, IDs DTC-039 y RNG DTC-037, y rechaza distribución explícita; DTPD-111 solo conecta mapas/meiosis en F7.1. Los bloques se muestrean independientemente sin LD interbloque objetivo, mientras DTF-044 reserva esa dependencia posterior sin retirar la recombinación fundadora intrabloque. Tres TSV de salida con esquemas exactos y el manifiesto conservan asignaciones, segmentos, frecuencias y origen/proceso reconstruibles. DTF-040 a DTF-044 delimitan ampliaciones futuras |
| <a id="dtc-056"></a>DTC-056 | DTPD-101 se resuelve con una única representación interna base 0 semiabierta. Una posición de base, un punto interbase y un intervalo no son intercambiables; cada valor viaja etiquetado por `reference_id` y `sequence_id`. Configuración, mapa de reducción y demás contratos propios externos conservan base 1 cerrada, y cada formato ajeno mantiene su convenio oficial. REQ-REF-29/35/36, REQ-IND-38 a REQ-IND-40, REQ-VAL-53 y ACE-113 | Evitar desplazamientos silenciosos y disponer de un contrato común para reducción, variantes, segmentos, anclas y composición | Las magnitudes usan `u64` y aritmética comprobada. Un intervalo interno `[start,end)` exige `start < end` y longitud `end-start`; un punto interbase admite `0..=sequence_length`. La conversión externa `[s,e]`↔interna `[s−1,e)` es exclusiva de las fronteras, exacta y reversible; desbordamiento, referencia/secuencia incompatible o valor fuera de rango se rechaza sin recorte. Los DTO externos no atraviesan el dominio y F2 extiende el mismo servicio sin crear otra convención |
| <a id="dtc-057"></a>DTC-057 | DTPD-116 se resuelve mediante `reduced_reference_selection_v1`. Para objetivo `T` y región `R`, calcula `q=T/R` regiones completas y `r=T%R`; preflight demuestra que existe una disposición completa. Selecciona primero las `q` regiones y después, si `r>0`, un único residual. En cada paso son candidatos exclusivamente los inicios no solapantes que preservan al menos una terminación completa, y todos esos inicios del conjunto de cromosomas y *scaffolds* tienen la misma probabilidad. REQ-CONF-126/127, REQ-REF-24 a REQ-REF-36, REQ-VAL-52/53, REQ-MAN-89/90 y ACE-112/113 | Obtener tamaño exacto y azar por oportunidad genómica sin ponderación artificial por secuencia, reintentos arbitrarios ni resultados parciales | Los inicios seguros se agrupan en rangos disjuntos ordenados por el orden FASTA y coordenada; se cuentan con aritmética comprobada y un índice entero uniforme sin sesgo selecciona una oportunidad mediante el RNG local. El orden del sorteo se registra, pero materialización y mapa se ordenan después por secuencia/coordenada. Cada secuencia derivada conserva ID y clase; el mapa JSON `reference_reduction_map` versión `1` publica intervalos origen/derivados base 1 cerrados, orientación directa, orden y clase completa/residual. Un caso inviable falla antes de consumir selección o publicar. Cambiar oportunidades, consumo RNG o muestreador exige nueva versión. Fixtures exhaustivos, vectores congelados, prueba estadística predeclarada y reconstrucción independiente actúan como oráculos; no se añade una dependencia por esta decisión |
| <a id="dtc-058"></a>DTC-058 | **Superada para F2 por DTC-079.** Conserva la historia de la concreción inicial de Marimo, pero sus protocolos, almacén de sesión, endurecimiento y puertas experimentales no son obligaciones de F2. | Evitar que una decisión técnica histórica mantenga infraestructura sin consumidor | Las versiones efectivas quedan fijadas por `pyproject.toml` y `uv.lock`; F2 prueba instalación desde el repositorio, arranque local y recorrido con el binario. Distribución, plataforma oficial y endurecimiento se verifican en F9. |
| <a id="dtc-059"></a>DTC-059 | En F2 la longitud de cada inserción aleatoria y la de cada deleción aleatoria son parámetros exactos, independientes y configurables, con 1 pb como valor predeterminado para ambas. Cada valor debe cumplir `1 <= longitud < umbral_indel_estructural_efectivo`; no se sortea una longitud ni se introduce una distribución de longitudes. REQ-CONF-55, REQ-IND-40/41 y ACE-29 | Formalizar la simplificación funcional aprobada sin convertir el valor predeterminado en una longitud fija obligatoria | Configuración efectiva y manifiesto conservan cada longitud y su procedencia. Cero, un valor igual o superior al umbral o una distribución solicitada se rechazan antes de simular. Una ampliación probabilística queda fuera de F2 |
| <a id="dtc-060"></a>DTC-060 | `prepare-snv-matrix` obtiene la identidad y los bytes autoritativos de referencia desde un manifiesto general F1 validado: recibe `--reference-manifest <manifest.json>` y `--reference-id <id>`, resuelve exactamente la entrada `references[]` indicada y verifica su ruta, tamaño y SHA-256 antes de preparar la matriz. RAF2-019; REQ-CONF-128, REQ-IND-60, REQ-VAL-55 y REQ-MAN-91 | Evitar deducir una identidad estable desde una ruta FASTA o aceptar ID y bytes desacoplados | No admite un `--reference-id` aislado ni un FASTA alternativo. Rechaza ID ausente/distinto, mismos bytes bajo otro ID, mismo ID con otros bytes, manifiesto no válido y referencia derivada incoherente; una referencia derivada válida conserva su propio ID y procedencia F1 |
| <a id="dtc-061"></a>DTC-061 | DTPD-102 se resuelve con un estado interno híbrido y un fichero externo obligatorio `truth.jsonl` independiente de los tipos Rust. Cada copia interna mantiene una cadena ordenada, lógicamente inmutable y compartible de segmentos: los ordinarios referencian intervalos de referencia/fundador sin copiar bases; SNP e *indels* pequeños son modificaciones tipadas, y solo una inserción sin fuente almacena sus propios bytes. Segmentos/modificaciones conservan procedencia, copia/fase y coordenadas DTC-056; la secuencia completa se materializa solo bajo demanda. `truth.jsonl` v1 usa UTF-8 sin compresión y JSON Lines versionado. La primera línea es una cabecera con clase/versión, ejecución, referencia, convenio, recuentos y SHA-256; siguen registros tipados de copias haplotípicas con cadena ordenada de segmentos/modificaciones/fase/procedencia y registros de genotipo por caso/locus con alelos, dosis y copias. REQ-IND-28/29, REQ-VER-27/31, REQ-SAL-02 y ACE-10/11/28/34/62 | Evitar duplicar genomas, mantener una verdad reconstruible e interoperable y escalar por flujo sin acoplar el contrato externo al contenedor interno | Los esquemas por tipo fijan claves y dominios; el orden canónico es referencia, individuo, secuencia y copia, seguido por las claves internas acordadas. La cabecera declara recuentos y digests de fuentes/grupos; el manifiesto conserva longitud y SHA-256 del fichero completo. Escritura/lectura son por flujo. Faltas, duplicados, desorden, recuentos/digests discordantes o relaciones rotas fallan. Un validador semántico/reconstructor independiente materializa secuencias y exige que cada genotipo coincida exactamente con sus haplotipos. VCF deriva de `truth.jsonl`, nunca al contrario. Recombinación/SV futuras componen el mismo modelo; compresión queda fuera de F2 y no cambia su lógica |
| <a id="dtc-062"></a>DTC-062 | La única versión VCF de salida de v1 es VCF 4.5 y su primera línea de cabecera es exactamente `##fileformat=VCFv4.5`. Se usan campos y representaciones estándar de VCF 4.5; cuando varias formas sean válidas se evita una novedad innecesaria si existe una representación ampliamente interoperable. No se producen alternativas VCF 4.2, 4.3 o 4.4. Una herramienta externa basada en htslib forma parte de la validación. Resuelve la versión dentro de DTPD-103; REQ-FMT-05/06, REQ-VER-06 y ACE-19 | Ofrecer un solo contrato de salida moderno sin multiplicar variantes de compatibilidad y conservar validación independiente | La preferencia interoperable no autoriza apartarse de VCF 4.5 ni inventar campos no documentados. DTC-066 integra esta decisión en el cierre de DTPD-103 |
| <a id="dtc-063"></a>DTC-063 | La salida VCF predeterminada es `.vcf.gz` comprimida con BGZF y lleva exactamente un índice: TBI si todas las secuencias y posiciones caben en su límite normativo, aproximadamente 2^29 bases por secuencia, o CSI si alguna lo supera. Nunca se publican ambos índices. Como alternativa explícita se admite `.vcf` de texto sin compresión y sin índice. Resuelve compresión e índice dentro de DTPD-103; REQ-FMT-05/06, REQ-SAL-12, REQ-MAN-58 y ACE-61 | Entregar una salida consultable e interoperable por defecto, sin imponer el límite de TBI a referencias grandes ni crear combinaciones ambiguas de artefactos | `EffectiveConfig` y el manifiesto registran modo, tipo de índice y rutas efectivas. VCF BGZF e índice se validan antes de hacerse visibles y se publican como una unidad atómica; el fallo de cualquiera impide publicar el conjunto. En modo texto se registra la ausencia deliberada de índice. DTC-066 integra esta decisión en el cierre de DTPD-103 |
| <a id="dtc-064"></a>DTC-064 | La proyección VCF aplica una normalización interna determinista y versionada a partir de los haplotipos finales de `truth.jsonl`: recorte de prefijo/sufijo común sin perder el anclaje requerido, alineamiento a la izquierda contra la referencia y orden por catálogo de referencia y coordenada. `bcftools norm` es solo un oráculo externo de prueba, no una dependencia de ejecución. En F2 existe como máximo un registro por `CHROM+POS`; todos los alelos finales compatibles se integran como `ALT` de un único registro multialélico y no se descomponen automáticamente. Resuelve normalización/agrupación dentro de DTPD-103; REQ-FMT-05/06, REQ-VER-15/16/17/30/44 y ACE-25/29/30 | Derivar una representación canónica repetible de la verdad haplotípica sin introducir una segunda autoridad ni perder la relación entre alelos del mismo sitio | Ante cambios solapados se reconstruye primero el alelo final de cada haplotipo y se usa un intervalo `REF` común suficiente para representar alelos de longitudes distintas antes de recortar y alinear. Configuración efectiva/manifiesto registran la versión del algoritmo. La restricción de un registro por posición se limita a SNP e *indels* pequeños de F2: SV/translocaciones futuras podrán requerir varios BND relacionados, incluso en una misma posición, conforme al estándar. DTC-066 integra esta decisión en el cierre de DTPD-103 |
| <a id="dtc-065"></a>DTC-065 | Todo registro VCF declara y usa `INFO/GPS_CASE` con `Number=.` y `Type=String`. Contiene los `case_id` originales únicos, en orden canónico, de todos los casos de `truth.jsonl` que aportan algún ALT; no se crea un caso compuesto y un caso aparece una sola vez aunque aporte varios ALT. `truth.jsonl` conserva para cada alelo final la relación ALT↔lista ordenada de casos fuente; si un ALT procede de varios casos, conserva todos. Un caso que abarca varios registros repite su ID en cada uno. Solo si el VCF pierde detalle añade `INFO/GPS_LOSS`, código cerrado/versionado y nunca texto libre, e `INFO/GPS_COMPLEMENT`, ID de la entrada enlazada. Resuelve RAF2-030 y `INFO`/complemento dentro de DTPD-103; REQ-VER-19/23 a REQ-VER-26, REQ-MAN-27/31, REQ-SAL-12 y ACE-33/108 | Mantener un VCF válido y una trazabilidad muchos-a-muchos explícita sin perder identidades originales, duplicarlas ni fabricar un caso compuesto | La cabecera declara los tres campos con número, tipo y descripción conformes a VCF 4.5. F2 prueba dos ALT/dos casos, un ALT/dos casos, un caso/varios registros, reordenación de entrada, duplicados e ID inexistente; como todos sus casos son representables, `GPS_LOSS`, `GPS_COMPLEMENT`, entradas y fichero permanecen ausentes. Cada conjunto genera como máximo un `vcf-complement.jsonl` versionado y solo ante pérdida real. La primera prueba positiva corresponde a F5.5 con un caso complejo/SV/BND real; allí se fijará la cardinalidad positiva multicaso del complemento sin alterar la lista `GPS_CASE`. Una pérdida literal descubierta en F2 puede probarse sin convertirse en puerta. DTC-066 integra esta decisión en el cierre de DTPD-103 |
| <a id="dtc-066"></a>DTC-066 | DTPD-103 queda resuelta integrando DTC-062 a DTC-065 con una puerta externa de aceptación/CI fijada a `bcftools` 1.24 y `tabix`/HTSlib 1.24. Estas herramientas nunca son dependencia de ejecución del simulador. VCF de texto y BGZF deben superar lectura completa, comprobación de normalización contra el FASTA autoritativo y, para BGZF, apertura/consulta regional del único TBI o CSI exigido por DTC-063. EBI `vcf-validator` no es autoridad porque su soporte publicado no cubre VCF 4.5. REQ-FMT-05 a REQ-FMT-07, REQ-VER-06/15 a 17/23 a 26/44, REQ-VAL-28, REQ-SAL-12 y ACE-19/25/30/33/61 | Cerrar un contrato VCF único con aceptación reproducible y dos validaciones independientes, sin convertir herramientas externas en productoras o reparadoras de la salida | Un oráculo semántico independiente, sin llamar ni compartir la proyección productora, lee FASTA, `truth.jsonl`, VCF y el `vcf-complement.jsonl` opcional; comprueba REF/ALT/GT, dosis, fase, muestras, orden, unicidad `CHROM+POS` en F2, lista `GPS_CASE` canónica/sin duplicados y relaciones ALT↔casos existentes, ausencia F2 de `GPS_LOSS/GPS_COMPLEMENT`, cardinalidad caso↔registros y campos prohibidos. El corpus incluye dos ALT/dos casos, un ALT/dos casos, un caso/varios registros, orden de entrada, duplicado e ID inexistente, además de positivos texto, BGZF+TBI/BGZF+CSI y mutaciones individuales de cada invariante. La cardinalidad positiva del complemento se valida desde F5.5. Las invocaciones exactas y la comparación sin reparación quedan fijadas en la resolución de DTPD-103 |
| <a id="dtc-067"></a>DTC-067 | DTPD-120 se resuelve mediante el algoritmo normativo `poisson-count-v1`. Admite `lambda` `f64` finita en `0.0..=2^52` (`4_503_599_627_370_496`), derivada conforme a DTC-022/026. `lambda=0` devuelve 0 sin consumir RNG; `0<lambda<10` usa inversión/Knuth y `lambda>=10` PTRS de Hörmann, con frontera inclusiva exacta. Cada uniforme abierto usa exactamente un `RngCore::next_u64()` del ChaCha12 local DTC-037 y `u = (((x >> 12) as f64) + 0.5) × 2^-52`. `exp`, `log`, `lgamma`, `floor` y `sqrt` proceden exclusivamente de `libm = =0.2.16`, nunca de matemática estándar. REQ-CONF-23 a REQ-CONF-25, REQ-IND-16 a REQ-IND-18, REQ-CONF-129 y ACE-10/109 | Obtener recuentos Poisson reproducibles entre Linux y WSL2, con consumo RNG y límites explícitos, sin hacer de la media un límite biológico | `lambda>2^52` se rechaza antes de consumir RNG con valor/límite/ámbito. Todo candidato fuera de `0..=2^53-1` o no representable exactamente como entero falla numéricamente; no hay cast saturante, recorte ni resorteo. Un recuento mayor que las oportunidades falla por la política existente, sin cambiar `lambda`. Se registran `algorithm_version`, `lambda`, ámbito/subsemilla y `sampled_count`. Vectores congelados cubren cero, subnormales, ambos lados de 10, medias representativas, *next-down*/2^52/*next-up* y máximo; un *hook* de prueba fuerza la rama rara fuera de rango. Oráculo ajeno y pruebas PMF/media/varianza —con `rand_distr` solo opcional/test-only— no reutilizan el productor. `libm` debe superar DTC-043, incluida licencia, auditoría y `build.rs` |
| <a id="dtc-068"></a>DTC-068 | DTPD-121 se resuelve con dos perfiles humanos compactos, aproximados, explícitos y nunca predeterminados. `human_1000g_global_sfs_approx_v1` resume los recuentos publicados finales de 1000 Genomes con pesos 64:12:8 para `[1/5008,0.005)`, `[0.005,0.05)` y `[0.05,1)`; los dos primeros usan `log_uniform`, el tercero `linear_uniform` y DTC-049 normaliza. El límite superior 1 pertenece a la cobertura, pero el muestreo continuo nunca lo devuelve; `alt_af=1` exige `fixed=1`. `1/5008` procede de una copia entre las 5.008 copias autosómicas de 2.504 muestras y pertenece a esta receta, no es un límite universal. `human_harris2017_snv12_approx_v1` deriva de los 96 contextos de Harris y Pritchard (2017): normaliza por separado AFR/AMR/EUR/EAS/SAS, promedia las cinco filas con el mismo peso, agrega contextos a seis clases colapsadas por hebra, refleja cada clase a su sustitución complemento inversa con mitad del peso y normaliza los 12 cambios dirigidos. REQ-CONF-129, REQ-IND-61, REQ-VAL-55, REQ-MAN-91 y ACE-116 | Ofrecer dos ejemplos pequeños, reproducibles y descriptivos sin descargar ni procesar el VCF poblacional, inventar perfiles exactos por superpoblación o presentar una aproximación como tasa humana universal | Fuentes fijadas: artículo final 1000 Genomes `https://www.nature.com/articles/nature15393`; artículo eLife `https://elifesciences.org/articles/24284`; datos fuente DOI `https://doi.org/10.7554/eLife.24284.004`; fichero `https://cdn.elifesciences.org/articles/24284/elife-24284-fig1-data1-v2.txt`, 3.644 bytes, SHA-256 `970cc71578fc008ef4e0d75c7c460ab874073da7c5b5226c526b387627c97a7f`. Artículo y datos eLife se atribuyen bajo CC BY. La receta pequeña versiona literales, parser, orden, normalizaciones, complementos inversos y SHA-256 de los recursos incorporados; F2.2.6 la regenera y contrasta con un oráculo independiente. F2 no depende de `bcftools` para estos perfiles ni descarga/procesa el VCF Phase 3. Perfiles suministrados por población siguen admitidos. Materialización exacta del VCF y perfiles empíricos exactos pasan a DTF-016/post-v1 |
| <a id="dtc-069"></a>DTC-069 | **Aplazada por DTC-079.** El conjunto multidocumento autocontenido se implementará cuando exista el primer consumidor real de escenarios. En F2, Marimo importa y exporta un único TOML general de la vía directa. | No anticipar un grafo de documentos, copias privadas y defensas TOCTOU que F2 no utiliza | El TOML debe seguir siendo válido para el binario, conservar IDs, orden y rutas relativas existentes y no sobrescribirse accidentalmente. La fase consumidora de escenarios recuperará los requisitos multidocumento que sigan siendo necesarios. |
| <a id="dtc-070"></a>DTC-070 | RAF2-029 se resuelve para el modo `exact` calculando primero las cuotas enteras por fila mediante mayores restos y desempate por el orden estable del espectro. Después, cada población permuta reproduciblemente el conjunto completo de `variant_id` canónicos con una subsemilla exclusiva de asignación y aplica las cuotas siguiendo esa permutación. Las poblaciones permutan independientemente y se prohíbe asociar filas por coordenada u orden canónico de variante. REQ-CONF-129, REQ-IND-61, REQ-IND-63, REQ-VAL-55 y ACE-116 | Conservar las cantidades exactas sin crear correlación espacial artificial ni acoplar poblaciones o el valor sorteado a la asignación de intervalos | El valor dentro del intervalo se obtiene después con un flujo RNG separado por población×`variant_id`. Una fila `fixed` participa en la asignación pero no consume RNG para el valor; una fila con cuota cero no consume RNG. Los oráculos demuestran cuotas exactas, estabilidad ante reordenación de entrada, ausencia de asociación con coordenadas, aislamiento entre poblaciones, consumo cero en ambos casos y separación asignación/valor. DTC-071 fija los nombres, componentes, orden y versión concretos de ambos ámbitos sin alterar esta semántica |
| <a id="dtc-071"></a>DTC-071 | RAF2-028 se resuelve con `subseed-scopes-v4.tsv`, catálogo normativo v4 y superconjunto estricto de v3. Conserva BLAKE3, contexto, codificación y `derivation_version=1`, copia literalmente todas las filas/vectores v3 y añade once ámbitos F2 definidos en §3.1. El sorteo de ALT usa `scenario_id`, `reference_id`, `sequence_id`, `position_1based` y `ref`, todos conocidos antes de ALT; `variant_id` solo se deriva después. En la vía directa sin documento de escenario, `scenario_id` es el literal interno reservado `"direct"`: no representa biología, no lo escribe la persona usuaria y ningún escenario real puede usarlo; queda registrado en `EffectiveConfig`, manifiesto y contexto RNG. REQ-CONF-111 a REQ-CONF-115, REQ-IND-63, REQ-VAL-55 y ACE-23/116 | Eliminar la dependencia circular de ALT y hacer que todo consumidor F2 tenga identidad estable, componentes ordenados y aislamiento verificable sin cambiar resultados F1 | Ruta, reloj, hilo, orden de iteración y posición accidental siguen prohibidos. Valores enteros de componentes usan decimal ASCII canónico sin signo ni ceros iniciales salvo `0`; IDs usan su forma validada. El catálogo y el manifiesto registran versiones 4/1 y cada ámbito enlaza algoritmo consumidor. Vectores BLAKE3 con raíz común y mutación aislada de cada componente demuestran que v1–v3 permanecen byte a byte inmutables |
| <a id="dtc-072"></a>DTC-072 | RAF2-033 se resuelve con la familia normativa de consumidores RNG v1 de §3.2: `open-unit-f64-v1`, `uniform-u64-below-v1`, `categorical-f64-v1`, `fisher-yates-v1`, `weighted-without-replacement-v1`, `dynamic-indel-selection-v1`, `linear-uniform-v1`, `log-uniform-v1`, `hwe-composition-v1`, `f-diploid-genotype-v1` y `phase-permutation-v1`. DTC-067 continúa siendo el único Poisson. REQ-CONF-111 a REQ-CONF-115, REQ-CONF-129, REQ-IND-16 a REQ-IND-18, REQ-IND-61 a REQ-IND-63 y ACE-23/107/116 | Fijar transformaciones, fronteras y consumo RNG para que cambiar una utilidad o la plataforma no altere silenciosamente matriz, indels, genotipo o fase | Cada algoritmo declara preflight, orden canónico y consumo cero; registra versión y ámbito. Todo intervalo es `[min,max)`: puede declarar `max=1`, pero nunca produce 1, y la corrección de redondeo lineal impide producir 0; solo `fixed` produce exactamente 0 o 1. La selección SNV es estática; la de indels reenumera el estado vigente después de cada evento aceptado. `log`/`exp` usan exclusivamente `libm = =0.2.16` ya sometido a DTC-043. Vectores literales, RNG instrumentado, enumeración exhaustiva pequeña y oráculos independientes preceden a todo consumidor |
| <a id="dtc-073"></a>DTC-073 | **Simplificada por DTC-079.** Marimo muestra y transmite únicamente opciones y autorizaciones implementadas por el binario F2 y descritas por su esquema/interfaz vigente. | Evitar anunciar opciones futuras sin crear un protocolo paralelo de capacidades | `--allow-empty-snv-matrix` y `--allow-incompatible-genotypes` aparecen cuando sus recorridos Rust existen; una edición o cambio de autorización invalida la validación. No se mantienen listas separadas `catalogued`/`available`. |
| <a id="dtc-074"></a>DTC-074 | RAF2-035 se resuelve con evidencia oficial de plataforma detectada, no autoetiquetada, en F9. El verificador registra `/etc/os-release`, kernel, arquitectura, triple Rust, marcadores WSL2, punto/tipo de montaje del árbol y proveedor/job/commit de CI; compara esos hechos con la clase solicitada y rechaza cualquier discordancia. Ubuntu 24.04 nativo exige `ID=ubuntu`, `VERSION_ID=24.04`, x86-64 y ausencia de marcadores WSL; WSL2 exige esos marcadores y un sistema de ficheros Linux fuera de `/mnt`. REQ-UI-13, REQ-CONF-111 a REQ-CONF-115 y ACE-115 | Evitar que dos ejecuciones en el mismo WSL2 se presenten como evidencia de dos plataformas distintas | F2 conserva solo actas locales de instalación/arranque y no afirma compatibilidad ni distribución oficial. F9 enlaza el job nativo identificable y el WSL2 detectado; ninguno sustituye al otro. La evidencia automática conserva versiones de Python/Marimo/uv/navegador y las actas manuales se separan. Una etiqueta de CLI es solo una expectativa que debe concordar con la detección |
| <a id="dtc-075"></a>DTC-075 | F2 conserva intervalos de frecuencia semiabiertos, separación entre selección estática de loci y dinámica de indels, y consumo explícito del recuento candidato. La división en un productor, expediente y paso independiente para cada primitiva queda superada por DTC-079. | Mantener la reproducibilidad científica sin convertir utilidades internas en productos separados | Los algoritmos pueden agruparse por módulo y probarse con vectores y casos de borde. Solo los contratos observables y los cambios capaces de alterar resultados requieren versión y trazabilidad propias. |
| <a id="dtc-076"></a>DTC-076 | **Superada por DTC-078.** Propuso que F1 definiera un `ResolutionDomain` sintético para contextos de fases posteriores. REQ-CONF-117/118/123 y REQ-MAN-88/90 | Conserva la historia de la revisión sin mantener una obligación prematura | No es operativa: cada fase propietaria introducirá sus relaciones reales junto con el primer consumidor |
| <a id="dtc-077"></a>DTC-077 | **Superada por DTC-078.** Propuso publicar en F1 una evidencia independiente de resolución, fijar globalmente IDs de fase y anticipar el orden de parámetros futuros. REQ-CONF-117/118/123 y REQ-MAN-88/90 | Conserva la historia de la revisión sin convertir el mecanismo provisional en contrato | `resolution-evidence` deja de pertenecer a F1; las fases propietarias fijarán evidencia, identidad y orden cuando existan datos reales |
| <a id="dtc-078"></a>DTC-078 | F1 adopta un modelo de confianza local mínimo. Valida como no confiables los TOML y FASTA aportados; después confía en los ficheros creados dentro de su provisional privado frente a manipulación activa y los protege frente a fallos propios, interrupciones y corrupción accidental con estados, tamaño, SHA-256, `INCOMPLETE` y publicación atómica no-clobber. La raíz TOML activa de F1 contiene solo clase, versión, semilla opcional, referencia y referencias reducidas. No admite todavía vía, ploidía, reglas, entidades poblacionales, eventos ni escenarios. La composición simulada A/C/G/T usa un validador específico con tolerancia efectiva fija `1e-9`, sin conservar el sistema genérico de normalización poblacional. El cierre estricto de F1 comprende 54 REQ completos (4/1/6/22/21 en F1.2–F1.6); cualquier comprobación F1 de un REQ reasignado es evidencia parcial, no cierre. | Separar reproducibilidad científica e integridad accidental de un modelo de atacante local activo, y evitar construir biología futura sin consumidor | Población↔subpoblación, fase↔generación, locus↔posición, vía, ploidía, reglas, eventos técnicos, normalización poblacional y sus manifiestos se implementan en sus fases propietarias. `ResolutionDomain` sintético, `resolution-evidence`, validación exhaustiva del árbol contra symlink/hardlink/inodos, ventanas TOCTOU y revalidación combinatoria del paquete salen de F1; el endurecimiento frente a atacante local activo se evalúa en F9/CI según el despliegue. Se conservan rutas relativas canónicas, permisos privados razonables, ausencia de sobrescritura y validadores externos pequeños de TOML, FASTA, SHA y repetibilidad. Los IDs exactos y sus destinos pertenecen a la matriz de trazabilidad vigente. |
| <a id="dtc-079"></a>DTC-079 | F2 adopta una implementación proporcional al uso local: repositorio, persona usuaria y provisional privado se consideran confiables frente a manipulación activa. Marimo configura el TOML general vigente, invoca directamente `validate`/`run` y muestra diagnósticos y resúmenes. F2 no crea escenarios sin consumidor, un ejecutor general de evidencias, protocolos paralelos de capacidades, almacenamiento persistente de sesiones ni defensas contra symlinks, inyección deliberada del entorno o carreras hostiles. El aprovisionamiento inicial explícito con `uv sync --frozen` puede usar Internet para descargar las dependencias resueltas por `uv.lock`; el arranque posterior y todo uso normal son locales, sin red ni servicios remotos. | Aplicar a toda F2 la simplificación que permitió cerrar F1 y evitar que la revisión técnica expanda el modelo de amenaza | Las entradas externas siguen validándose; Rust conserva autoridad semántica, reproducibilidad, publicación no-clobber y pruebas de invariantes biológicas. Cada hito usa pruebas Cargo y validadores de formato pertinentes, con una evidencia de cierre sencilla. La instalación distribuible sin descarga previa, `.tar.gz`, OCI, plataforma oficial, auditoría completa de accesibilidad y endurecimiento del despliegue quedan en F9. DTC-058/069/073 se interpretan según esta decisión; DTC-075 mantiene semántica pero no su partición documental. |

### 3.1. Catálogo normativo de subsemillas F2 v4

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
| `initial_genotype` | `scenario_id`, `reference_id`, `population_id`, `individual_id`, `variant_id` | `hwe-composition-v1` o `f-diploid-genotype-v1` |
| `initial_phase` | `scenario_id`, `reference_id`, `population_id`, `individual_id`, `variant_id` | `phase-permutation-v1` |

`variant_class` solo admite `insertion` o `deletion`. `spatial_scope_id` es el
ID estable DTC-039 de la regla espacial efectiva, con `global` como literal
reservado cuando no existe región. `event_index` es el ordinal base 0 del evento
aceptado dentro del grupo canónico; la primera inserción usa `0`. La posición
base 1 y el ordinal se codifican como decimal ASCII canónico. La selección de
ALT precede a la derivación DTC-039 de `variant_id`, que sí incorpora ALT.

### 3.2. Consumidores RNG normativos F2

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

**Concreción aprobada de DTC-055 para el panel y los modos de construcción.** Una única
orden Rust, conceptualmente `prepare-haplotype-panel`, recibe referencia y los
cuatro TSV, valida sus relaciones, los copia en orden canónico, calcula SHA-256,
genera `panel.metadata.json` y publica atómicamente sin modificar originales ni
sobrescribir destinos. El runtime solo acepta el paquete completo y vuelve a
verificar hashes, referencia, clases, versiones y rutas; Marimo invoca la misma
orden. La conversión desde VCF queda en DTF-041.

El paquete enlaza exactamente una referencia. Todo `reference_id` de
`panel_blocks.tsv` coincide con `metadata.reference.reference_id`; el SHA-256
recalculado del FASTA autoritativo coincide con
`metadata.reference.sha256`; y cada `sequence_id` existe en ese FASTA. Mezclar
referencias, declarar otro ID o digest o referir una secuencia inexistente se
rechaza antes de interpretar individuos.

Los cuatro TSV usan UTF-8 sin BOM, LF final y tabuladores, sin comillas, vacíos,
líneas vacías ni campos extra. Sus cabeceras exactas son, respectivamente,
`block_id/reference_id/sequence_id/start/end`,
`block_id/site_id/position/ref/alt`,
`block_id/haplotype_id/site_id/allele_index` y
`population_id/block_id/haplotype_id/frequency`. Las barras solo separan columnas
en esta descripción. `alt` contiene ALT A/C/G/T separados por comas sin espacios;
su orden define 1..N y `ref`, 0. Enteros y `f64` usan las formas decimales
canónicas vigentes. El orden estable y las claves únicas siguen bloque, orden
FASTA/posición, HAP e identificadores ASCII conforme al esquema detallado en
REQ-EVO-74.

Cada vector completo de `allele_index` es único dentro de su bloque, incluso si
su frecuencia es 0; los duplicados se rechazan informando todos los IDs, sin
fusionarlos ni sumar frecuencias. Cada bloque se trata como un locus multialélico
independiente. HWE/multinomial es predeterminado para P=1..10; F reutiliza las
fórmulas multialélicas solo en diploides y, en construcción intacta, una
distribución explícita puede prevalecer localmente. En modo probabilístico cada copia se sortea con reemplazo.
En modo controlado, un vector entero por HAP que suma C forma un multiconjunto
que se asigna sin alterar totales; la incompatibilidad con HWE, F o la
distribución explícita se rechaza sin ajustar cantidades ni relajar el modelo.

Cada bloque exige al menos dos sitios SNP y dos vectores HAP distintos. Las
frecuencias no reducen el catálogo: una o todas las poblaciones pueden fijar un
único HAP y conservar los demás con frecuencia 0. Un sitio fijado en el mismo ALT
en todos los HAP se conserva; un sitio `REF` en todos se rechaza como redundante.
La configuración separa `haplotype_construction` (`intact`/`recombinant`) de
`frequency_sampling` (`probabilistic`/`controlled`). La construcción se resuelve
por valor general y sustituciones población×bloque conforme a DTC-033; una unidad
no mezcla modos. `intact` prohíbe el parámetro de recombinación y `recombinant`
exige `recombination_probability_per_interbase` finita en (0,1]. K=0 puede
ocurrir por azar sin convertir la unidad en intacta. La construcción intacta
admite ambos muestreos; la recombinante solo admite `probabilistic` y rechaza
`controlled` antes de individuos, identificando campos y ámbito. El control
exacto de fuentes o mosaicos recombinantes queda en DTF-042.

El modo intacto selecciona exactamente un HAP por
población×individuo×copia×bloque, copia conjuntamente su vector y registra
población, bloque, copia y HAP fundador. Un HAP de frecuencia 0 no se selecciona.
El estado previo a mutaciones es idéntico al panel; las mutaciones de novo se
trazan aparte y no cambian la identidad fundadora.

El modo recombinante inicializa directamente cada individuo y bloque y no
representa una meiosis histórica. Un productor selecciona una fuente primaria
según el panel y sortea las L−1 uniones con
`recombination_probability_per_interbase=r` finita en (0,1]. K=0 devuelve la
primaria; K≥1 sortea una secundaria y el constructor común cerrado por DTC-055
crea el mosaico con puntos únicos, ordenados e internos. DTPD-111 no bloquea
F6.4: decidirá para F7.1 cómo conectar mapas y semántica meiótica a ese mismo
constructor. El punto x es la unión x|x+1; la salida
empieza con primaria en `start`, alterna en cada punto y usa segmentos *forward*
cerrados, sin desplazar el cruce a un SNP ni introducir inversiones. Se conservan
en la traza cruces alélicamente invisibles y pares con efecto neto cancelado. No existe depósito prototipo ni segundo
muestreo, y no se corrigen o descartan mosaicos.

Primaria y secundaria proceden siempre del panel original y de las frecuencias
de la población. Los recombinantes ya creados no se incorporan como fuentes ni
actualizan probabilidades; el procesamiento es independiente del orden. La
recombinación acumulada y las fuentes descendientes quedan para F7/F8.

En diploides, Bernoulli(F) decide si se ejecuta un productor y se duplica su
vector en ambas copias o si se ejecutan dos productores independientes. Así se
aplican exactamente las fórmulas multialélicas: F=0 es HWE y F=1 homocigosis
completa, sin autofecundación ni pedigrí. Para P distinta de 2 se ejecutan P
productores independientes y se rechaza F. Las C=N×P copias finales se agrupan y
su frecuencia realizada es `n_h/C`; las frecuencias del panel son el objetivo de
las fuentes.

Un vector coincidente reutiliza el HAP local del panel. Uno nuevo recibe mediante
DTC-039 un ID tipado determinista derivado de `panel_id` + `block_id` + vector
canónico; no incluye individuo, orden, subsemilla ni evento. Se agregan vectores
idénticos, pero cada productor/copia conserva fuentes, puntos, individuo, copia y
subsemilla. Se distinguen productores con cruce y copias finales derivadas —una
salida duplicada por F cuenta uno y dos, respectivamente—.

Toda decisión aleatoria reutiliza DTC-037 con ámbitos por población, bloque,
individuo y, cuando aplica, copia, separando F/autozigosis, primaria, puntos y
secundaria. La misma raíz+panel+configuración+versión reproduce; reordenar o
añadir otras identidades y cambiar `r` en otro bloque no perturban el resultado.
Una salida duplicada por F mantiene un productor/procedencia y dos copias
enlazadas. Se registran versión de derivación, ámbitos y subsemillas.

Las salidas obligatorias son `founder_block_assignments.tsv` por copia y bloque,
`recombination_segments.tsv` por proceso y segmento y
`realized_haplotype_frequencies.tsv` por población, bloque y HAP final. El
manifiesto las enlaza por clase/versión/ruta/tamaño/SHA-256 y resume configuración
efectiva, F, r, RNG, productores con cruce, cruces totales, vectores nuevos y
coincidencias con panel, sin duplicar tablas.

Todo productor posee `process_id`, incluso con K=0. `origin_kind` vale únicamente
`intact` si K=0 o `recombined` si K≥1; una duplicación por F produce dos
asignaciones al mismo proceso y origen. Los procesos intactos tienen un segmento
primario que cubre el bloque completo. La tabla de frecuencias no agrega un
origen único: contiene recuento total, intacto y recombinado, exige que el total
sea su suma y calcula `frequency=copy_count/C`.

Las tres salidas usan exactamente las cabeceras tabuladas
`population_id/individual_id/copy_id/block_id/final_haplotype_id/origin_kind/process_id`,
`process_id/segment_index/block_id/start/end/source_haplotype_id` y
`population_id/block_id/haplotype_id/copy_count/intact_copy_count/recombined_copy_count/frequency`;
las barras solo separan columnas en esta descripción. Comparten UTF-8 sin BOM,
LF final, ausencia de comillas, vacíos y columnas extra, tipos canónicos y orden
y claves población+individuo+copia+bloque, proceso+índice y
población+bloque+HAP. Clase/versión se declaran en el enlace del manifiesto y
esquema, no como columnas.

La construcción intacta admite HWE, F diploide o una distribución genotípica
explícita. La recombinante admite solo HWE o F diploide; una distribución
explícita se rechaza antes de individuos con población, bloque y campos. Su
condicionamiento sobre resultados recombinantes queda en DTF-043.

Marimo expone una sección poblacional que elige directorio, invoca las acciones
Rust de preparar/validar, resume el panel, mapea poblaciones y edita la tabla
población×bloque. Muestra los controles compatibles por construcción, r tipado y
F solo diploide; previsualiza métricas devueltas por Rust, exige validación
vigente antes de `run` y después enlaza las tres tablas sin cargarlas completas.

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
externa obligatoria es `truth.jsonl` v1, UTF-8 sin comprimir, independiente de
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

Resuelta el 14 de agosto de 2026 mediante [DTC-057](#dtc-057). El algoritmo
`reduced_reference_selection_v1` pondera uniformemente cada inicio que preserve
una terminación completa en el espacio conjunto, selecciona regiones completas
y después un único residual y materializa solo tras ordenar por
secuencia/coordenada. El preflight rechaza configuraciones inviables sin
reintentos o resultados parciales. Vectores, enumeración exhaustiva pequeña,
prueba estadística y reconstrucción origen+mapa son oráculos obligatorios.

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

### DTPD-119. Formato y muestreo del panel de bloques haplotípicos (resuelta)

Resuelta el 13 de agosto de 2026 mediante [DTC-055](#dtc-055). El panel usa
bloques locales, cuatro TSV canónicos y `panel.metadata.json`, con preparación,
integridad, catálogos densos, identidad por contenido y exclusión completa frente
a loci `alt_af` dentro de los intervalos. La inicialización intacta y la
recombinante se configuran por población×bloque con ejes separados. La
recombinante usa exclusivamente muestreo probabilístico, HWE o F diploide,
productores directos con fuentes del panel original, motor común, RNG DTC-037 y
salidas reconstruibles. Distribuciones explícitas solo se admiten en intacto;
DTF-040 a DTF-044 conservan las ampliaciones aplazadas.

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

### 7.1. Decisiones pendientes

| ID | Decisión | Opciones | Recomendación razonada | Consecuencias | Bloquea |
|---|---|---|---|---|---|
| DTPD-104 | Adaptadores SAM/BAM/CRAM e índices | Crate; biblioteca nativa; herramientas externas; combinación | Comparar cumplimiento de estándar, CRAM con referencia exacta, ordenación, índices, cabeceras y distribución | Determina E/S, procesos, errores y rendimiento | F3.3 |
| DTPD-105 | Concreción del adaptador de proceso externo para minimap2 y otros mapeadores | Descubrimiento explícito o ruta configurada; política de invocación y sandbox; compatibilidad de versiones; clasificación de códigos, señales y fallos | Partir del proceso externo ya decidido por DTC-043 y concretar descubrimiento, comando reproducible, aislamiento, compatibilidad y diagnóstico sin ocultar versión, parámetros, stdout, stderr ni código de salida | Afecta instalación, seguridad, logs y reproducibilidad; no reabre servicio frente a proceso | F3.2 |
| DTPD-106 | Perfiles WGS cuantitativos incorporados | Perfiles versionados de Illumina y PacBio HiFi; sustituciones explícitas permitidas | Validar cada perfil con datos adecuados y distinguir claramente valores incorporados y sustituciones | Fija distribuciones de longitud, error, calidad y artefactos; no incluye estimación ni perfiles empíricos | F4.1–F4.3 |
| DTPD-107 | Calibración empírica (retirada de v1) | Sin elección técnica en v1 | La decisión funcional DTC-023 aplaza la capacidad completa; cualquier diseño del estimador y del artefacto corresponde a DTF-016 | No bloquea ninguna fase de v1 | — |
| DTPD-108 | Verdad y procedencia de lecturas | Registro por lectura; registros por bloques; combinación con agregados | Mantener para cada lectura o pareja final el origen y la diferencia verdadera/observada, y agregados solo para descartes o no asignadas; medir volumen antes de fijar formato, sin cambiar el contrato funcional | Puede dominar el almacenamiento y la validación | F4.4/F4.5 |
| DTPD-109 | Fórmulas de F diploide y Ne reproductivo realizado | Fórmulas publicadas aplicables a la heterocigosidad diploide y a las contribuciones observadas; estado no evaluable cuando corresponda | Evaluar contra el significado funcional confirmado y casos simulados; no adaptar F a P=1 o P>2 ni introducir un objetivo u optimización de Ne | Cambia estadísticas y pruebas, no los descendientes ni el genotipo; F fuera de diploides conserva «no aplicable» | F6.3 y F7.3 |
| DTPD-110 | Distribuciones y tolerancias poblacionales posteriores | Catálogos de distribuciones aún no elegidos; tolerancias específicas de modelos poblacionales o estadísticos posteriores | Aplicar [DTC-049](#dtc-049) a toda probabilidad o frecuencia configurable y decidir aquí solo las distribuciones y tolerancias específicas cuando su fase defina la semántica; ninguna puede sustituir `0,5/C` | No reabre `f64`, suma, normalización ni tolerancia general ya resueltos; afecta modelos poblacionales posteriores y sus pruebas | F6.1/F6.2 |
| DTPD-111 | Mapas, posiciones y comparación estructural para recombinación meiótica | Formatos tabulares candidatos; generación meiótica de puntos; interpolación de mapas; representación canónica comparable | Evaluar cobertura completa, monotonicidad, conversiones y equivalencia de segmentos, orden, orientación, copias y puntos de ruptura sin exigir el mismo ID histórico; conectar el resultado al constructor común DTC-055 | Afecta bivalentes, linaje y casos estructurales de F7.1; no reabre ni bloquea la variante uniforme fundadora de F6.4 y no puede crear otro constructor o introducir compatibilidad parcial | F7.1 |
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
REQ-CONF-106/107. DTPD-110 queda restringida a distribuciones y tolerancias
poblacionales específicas posteriores y no reabre ese contrato general.

Las decisiones todavía abiertas de la tabla DTPD deben resolverse antes de
ejecutar el hito que señalan. DTPD-117, DTPD-118 y DTPD-119 ya están resueltas
mediante DTC-079, DTC-054 y DTC-055, respectivamente. No bloquean
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
DTPD-118; F6.4 ya no
está bloqueada por DTPD-119; y F7.1 continúa bloqueada por DTPD-111. Las
recomendaciones pendientes no se consideran aprobadas hasta que consten
revisión y fecha.
