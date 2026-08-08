# Registro de decisiones técnicas del simulador

- **Versión:** 0.1
- **Estado:** borrador para revisión
- **Fecha:** 8 de agosto de 2026
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
| DTC-010 | La configuración usa coordenadas base 1 e intervalos cerrados; cada formato externo conserva su convenio. REQ-CONF-28 y REQ-IND-38/39 | Evitar posiciones ambiguas | La representación interna y las conversiones concretas siguen pendientes |
| DTC-011 | La reproducibilidad biológica exige misma semilla, configuración, entradas y versión. REQ-CONF-26/27 y REQ-REP-22/23 | Poder repetir decisiones sin prometer compatibilidad indefinida | DTC-038 concreta los cuatro niveles y no confirma identidad bit a bit ni garantías entre versiones o plataformas no verificadas |
| DTC-012 | La primera versión prioriza casos pequeños y no impone límites funcionales generales. REQ-NOF-01 a REQ-NOF-04 | Permitir validar primero y escalar con evidencia | Los límites prácticos y objetivos cuantitativos deben medirse y acordarse |
| DTC-013 | Los procesos y formatos externos se validan antes de declarar éxito. REQ-FMT-04/07, REQ-VAL-11 y REQ-SAL-05 | Evitar resultados plausibles pero incoherentes | Un artefacto fallido o incompleto nunca se publica como completo |
| DTC-014 | No se ejecuta código personalizado desde la configuración. REQ-CONF-108 a REQ-CONF-110 | Mantener reglas validables y reproducibles | El catálogo declarativo solo admite operaciones conocidas por el programa |
| DTC-015 | En v1, una unidad reproductiva solo contiene individuos basados en la misma referencia estructural. Pregunta 4 del informe adversarial (RAD-012); REQ-IND-53, REQ-MIG-09 y REQ-CRU-44 | No se ha definido la identidad de loci, los mapas ni la referencia de la descendencia entre estructuras distintas | Subpoblaciones y ramas aisladas pueden usar referencias distintas; el mapeo contra varias referencias sigue siendo independiente de la reproducción |
| DTC-016 | En v1, cada evento que genera FASTQ, BAM o CRAM es su único experimento de secuenciación, contiene una tanda y exactamente una biblioteca y un grupo de lectura por individuo. Pregunta 7 del informe adversarial (RAD-027); REQ-CONF-102, REQ-TEC-01 a REQ-TEC-04 y REQ-RG-01 a REQ-RG-07 | Mantener una unidad técnica simple, explícita y verificable por evento | Otra biblioteca del individuo, tanda, perfil, tecnología o control exige otro evento; los eventos solo de verdad o VCF carecen de entidades técnicas |
| DTC-017 | En v1, una sustitución por identificador individual solo puede dirigirse a un individuo directo, fundador o parental expreso existente durante la validación. Pregunta 8 del informe adversarial (RAD-026); REQ-CONF-17 y REQ-CONF-18 | Evitar configurar entidades que todavía no existen y aplicaciones retroactivas | Los descendientes usan ámbitos poblacional, temporal y espacial; un identificador inexistente es error |
| DTC-018 | La cobertura objetivo usa la longitud haploide de la referencia estructural del individuo y se distingue de las coberturas WGS básica, final del FASTQ y alineada por referencia. Pregunta 10 del informe adversarial (RAD-004 y RAD-017); REQ-WGS-04, REQ-WGS-11 y REQ-DAT-24 | Un mismo FASTQ no puede tener objetivos distintos por referencia de mapeo | Cada secuencia se suma una vez sin ploidía; haplotipos y referencias de mapeo no cambian el denominador; la cobertura alineada es solo una observación |
| DTC-019 | V1 no modela elementos transponibles como clase funcional. Pregunta 12 del informe adversarial (RAD-006); sección 3 y REQ-TE-01 a REQ-TE-14 retirados | Reducir el alcance y eliminar la contradicción entre referencias derivadas y eventos TE | No hay catálogos, familias, mecanismos ni metadatos TE; las copias e inserciones se clasifican mediante operaciones genéricas |
| DTC-020 | La generación aplica primero los eventos explícitos y después los aleatorios en el orden duplicaciones, deleciones, inversiones, translocaciones, repeticiones en tándem, *indels* y SNP. Pregunta 13 del informe adversarial (RAD-008); REQ-IND-54 a REQ-IND-57 | Hacer único y comprobable el resultado ante colisiones y agotamiento del espacio válido | Cada evento se elige del espacio válido actualizado; no hay reintentos sobre inválidos ni límites de intentos, y la insuficiencia hace fallar la generación conservando los recuentos sorteado, aceptado y activo |
| DTC-021 | Cada conjunto solicitado produce exactamente un VCF multimuesta estándar; un caso de verdad usa uno o varios registros y, solo si pierde detalle, una única entrada complementaria. Pregunta 14 del informe adversarial (RAD-011 y RAD-023); REQ-VER-04 y REQ-VER-19 a REQ-VER-26 | Mantener interoperabilidad VCF sin perder la verdad de casos complejos | El VCF sigue siendo válido por sí mismo; la reconstrucción completa usa el par enlazado, BND conserva su cardinalidad estándar y referencias distintas permanecen independientes |
| DTC-022 | V1 usa exclusivamente Poisson para recuentos de variantes aleatorias y el modelo escalonado simple de una unidad para mutaciones de repeticiones en tándem. Pregunta 17 del informe adversarial (RAD-020); REQ-CONF-25, REQ-CONF-61/62 y REQ-IND-18/35 | Cerrar el catálogo mínimo y hacer observables las configuraciones válidas | Eventos explícitos y modo controlado no son distribuciones alternativas; la tasa y el sesgo expansión/contracción siguen configurables y cualquier otra distribución, modelo o salto se rechaza antes de simular |
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
| <a id="dtc-037"></a>DTC-037 | La aleatoriedad de v1 usa generadores ChaCha12 nombrados y locales, nunca `StdRng` ni un flujo global. La semilla raíz y cada subsemilla tienen 256 bits y se representan en TOML, configuración efectiva y manifiestos como exactamente 64 hexadecimales minúsculos. Si falta la raíz, se obtienen 32 bytes del CSPRNG del sistema operativo antes de simular y se registran valor y origen. Resolución de DTPD-004; REQ-CONF-26/27, REQ-CONF-101/111 a REQ-CONF-115 y REQ-EVO-12/62 | Aislar biología, muestreos y eventos para que el consumo de un proceso o la planificación no alteren decisiones ajenas y permitir paralelismo reproducible sin compartir estado mutable | BLAKE3 deriva cada subsemilla con el contexto ASCII exacto `genome_population_simulator/v1/subseed`. El material binario canónico contiene, en este orden, la raíz cruda de 32 bytes, la versión de derivación como `u32` big-endian, el número de componentes como `u32` big-endian y, para cada componente en el orden normativo del catálogo de ámbitos, nombre de campo ASCII y valor UTF-8, ambos precedidos por su longitud `u32` big-endian. La identidad de flujo nunca procede de posición accidental, iteración, hilo o planificación. Se separan al menos población inicial, generaciones y procesos biológicos, muestreos y cada evento técnico, afinando por unidades independientes. Cada consumidor solo avanza su generador. Se registran raíz, origen, ChaCha12, versión de derivación y subsemillas de evento; los ámbitos internos se reconstruyen desde IDs y reglas versionadas. El catálogo de ámbitos, campos y orden y sus vectores son artefactos normativos de F1 y no cambian sin nueva versión de derivación. Cambiar ChaCha12, BLAKE3, esta codificación o un algoritmo de muestreo cambia la reproducibilidad. DTC-038 concreta la garantía; DTC-039 separa el dominio de IDs; DTC-043 rige dependencias; DTC-048 fija el lock; DTC-040 reserva SHA-256 para integridad publicada |
| <a id="dtc-038"></a>DTC-038 | V1 declara por separado cuatro niveles de reproducibilidad: biología; decisiones y lecturas técnicas; equivalencia semántica de salidas; e identidad binaria. Resolución de DTPD-005; REQ-CONF-26/27 y REQ-CONF-111 a REQ-CONF-115, REQ-REP-22/23 y REQ-MAN-69 | Permitir repetir y comparar cada clase de resultado sin deducir una garantía más fuerte desde otra más débil | La biología se garantiza con la misma versión del simulador, plataforma inicialmente admitida, semilla raíz, configuración efectiva, contenidos de entrada, perfiles y esquema y algoritmos aleatorios versionados. El nivel técnico exige además la misma subsemilla de evento y versiones y condiciones de las dependencias participantes. La equivalencia semántica compara contenido lógico y, ante herramientas externas, exige sus mismas versiones, parámetros y condiciones verificadas. V1 no garantiza identidad binaria por defecto para ninguna clase; solo podrá declararse por formato y condiciones tras una prueba expresa. No se promete reproducibilidad entre versiones del simulador ni entre plataformas no verificadas. El manifiesto registra por nivel `garantizado`, `no garantizado` o `no aplicable`, condiciones, versiones, perfiles, herramientas y dependencias. Toda semilla generada se registra y puede reutilizarse. Se mantienen comparadores separados para verdad lógica, contenido técnico, semántica y bytes. DTC-042 fija la plataforma oficial inicial; crates, toolchain y formatos adicionales permanecen separados |
| <a id="dtc-039"></a>DTC-039 | Cada entidad separa un identificador estable de una etiqueta humana Unicode opcional y no única. Los identificadores aportados son ASCII, sensibles a mayúsculas y cumplen conceptualmente `[A-Za-z0-9][A-Za-z0-9._-]{0,63}`; no admiten espacios, rutas ni otros caracteres. En el dominio se usan tipos distintos como `IndividualId`, `EventId` o `ReferenceId`, nunca `String` genérico en una frontera. Resolución de DTPD-006; REQ-CONF-01/38/39/51/52/101, REQ-IND-25/26, REQ-RG-03/04/07, REQ-DAT-68 a REQ-DAT-71 y REQ-VER-10/12 | Evitar confundir categorías, etiquetas, rutas o ámbitos y mantener identidades reproducibles compatibles con formatos genómicos | Los IDs automáticos se derivan con BLAKE3 y el contexto ASCII exacto `genome_population_simulator/v1/id`, usando el material canónico de DTC-037 y el catálogo normativo de campos. Toman los primeros 128 bits del resultado, los codifican en Base32 RFC 4648 mayúscula sin relleno `=` y anteponen el prefijo ASCII de clase del catálogo. El ID de ejecución incluye como causas mínimas raíz y digest de configuración; los demás incluyen tipo e IDs causales estables según catálogo, nunca hilo, ruta, orden físico o reloj. Los catálogos de campos y prefijos y sus vectores son artefactos normativos de F1 y no cambian sin nueva versión de derivación. Toda colisión se detecta y hace fallar. Se respetan los ámbitos normativos: variantes locales al conjunto/VCF, `RG` único en la ejecución, lectura o plantilla única en el evento y mates con el mismo ID de plantilla y roles distintos; el ID de individuo se conserva en `SM`, VCF y manifiesto. Rutas y nombres físicos se derivan o escapan sin adquirir autoridad. El manifiesto relaciona ID opaco, clase, ámbito, causa y etiqueta. Los IDs automáticos permanecen compactos y válidos para SAM, VCF y FASTQ; no se fija aquí el layout de directorios |
| <a id="dtc-040"></a>DTC-040 | SHA-256 es la única suma de comprobación publicada por v1. Se identifica como `sha256` y su digest se representa mediante 64 caracteres hexadecimales minúsculos calculados sobre todos los bytes exactos del fichero, sin normalización. Resolución de DTPD-007; REQ-REF-17/18/21, REQ-MAN-39/76 y REQ-SAL-19 | Proporcionar una identidad de contenido interoperable y verificable para cada entrada, salida y fichero enlazado | Cada manifiesto registra algoritmo, digest, tamaño en bytes y ruta relativa de cada fichero. La suma se calcula o verifica al copiar entradas, al cerrar salidas, antes de exportar y después de exportar; copia y original deben coincidir. Puede calcularse incrementalmente durante la E/S, pero las comprobaciones exigidas no se sustituyen por tamaño, fecha, metadatos ni muestreo parcial. Cada fichero conserva su digest; una suma adicional del paquete no reemplaza las individuales. MD5 solo se emite cuando un estándar genómico lo exige, separado y nunca como integridad principal. BLAKE3 queda reservado a derivación de subsemillas e IDs con dominios distintos. DTC-045 concreta la publicación; la futura dependencia SHA-256 debe superar DTC-043 y quedar fijada por el lock de DTC-048 |
| <a id="dtc-041"></a>DTC-041 | V1 admite FASTA de texto sin comprimir y encapsula `noodles-fasta` detrás de un contrato propio como adaptador Rust puro para lectura y escritura por flujo, FAI y acceso indexado. El dominio valida al menos una secuencia no vacía, primer token de cabecera no vacío y único y alfabeto A/C/G/T/N en cualquier caja, sin huecos ni otros símbolos IUPAC. Resolución de DTPD-008; REQ-REF-01 a REQ-REF-18, REQ-REF-21/23, REQ-VAL-45, REQ-MAN-86 y REQ-SAL-01 | Mantener el contrato de referencia independiente de una biblioteca concreta y admitir ficheros grandes sin perder el original ni relajar las reglas biológicas | Un FASTA suministrado se copia byte a byte sin reescribirlo, se valida y calcula su SHA-256 durante el flujo, deja intacto el original y solo se publica tras éxito. Las bases se normalizan a mayúsculas únicamente en la representación interna; se conservan copia exacta y cabecera cruda, con primer token como ID y el resto como descripción. Una descripción no UTF-8 se conserva sin pérdida mediante bytes codificados reversiblemente y una marca de codificación en salida y manifiesto. Los FASTA simulados o derivados se escriben canónicamente sin compresión, con A/C/G/T/N mayúsculas, LF, ancho fijo y orden estable; el ancho se versiona si afecta a bytes. Los índices solo se crean cuando se necesitan. `samtools faidx` es validador independiente de pruebas. La compresión FASTA se diagnostica como no admitida. Evaluación revisable a 2026-08-08: candidata `noodles-fasta` 0.66.0, MIT, MSRV 1.89.0 compatible con Rust 1.97.1, edición 2024 y `build = false`; sus dependencias directas son `bstr`, `memchr`, `noodles-bgzf`, `noodles-core` y `tokio` solo para la característica asíncrona. F1 usaría características mínimas y síncronas, sin activar async. La API upstream aún se considera experimental y evoluciona con frecuencia; adaptador propio, lock y pruebas mitigan el riesgo. La versión es candidata evaluada, no fijada aún en `Cargo.lock`: antes de aceptarla en F1 se revalidan seguridad y licencias del grafo final y ausencia de código nativo y scripts de construcción conforme a DTC-043 y DTC-048 |
| <a id="dtc-042"></a>DTC-042 | La plataforma nativa oficialmente validada para v1 es Ubuntu 24.04 LTS x86-64 GNU, con triple `x86_64-unknown-linux-gnu`, configuración de CPU portable y sin `target-cpu=native`. Otras distribuciones x86-64 GNU pueden funcionar, pero no reciben garantía hasta incorporarse a una matriz. Ubuntu 24.04 sobre WSL2, siempre en el sistema de ficheros Linux/ext4 y no en `/mnt/c`, es el entorno de desarrollo admitido y actual, no una plataforma de publicación independiente. Resolución de DTPD-009; REQ-CONF-111 a REQ-CONF-115, REQ-REP-22/23, REQ-MAN-05 y REQ-SAL-19; hoja de ruta F1.1 y F9.7 | Disponer de un entorno rápido de desarrollo y una referencia de publicación verificable sin convertir el triple GNU en una promesa genérica ni heredar particularidades de WSL2 | Las pruebas ordinarias se ejecutan en WSL2 sin dependencias exclusivas. Antes de publicar se valida en CI Ubuntu 24.04 nativo y en un entorno OCI cuya base y receta se fijan por digest; se prueban atomicidad, publicación no-clobber en su sistema de ficheros, señales, permisos, recursos y procesos externos. V1 distribuye un `.tar.gz` versionado con ejecutable, esquemas, licencias y ejemplos, cuyo SHA-256 se calcula sobre los bytes exactos del archivo, y una imagen OCI de referencia identificada por el digest `sha256:` de su descriptor de índice o manifiesto, conservando y verificando los descriptores de configuración y capas. Un archivo que exporte la imagen OCI tiene además su propio SHA-256 independiente. Se registran plataforma, arquitectura, versión del simulador, referencia y digest de imagen y herramientas externas. macOS, Windows nativo y ARM64 no tienen soporte oficial; ampliarlos exige dependencias y matriz de pruebas. No se promete identidad binaria entre compilaciones. DTC-043 conserva dependencias, DTC-045 la publicación y DTC-048 el toolchain |
| <a id="dtc-043"></a>DTC-043 | V1 prefiere crates Rust puros para configuración, serialización, hashes, RNG y adaptadores de formatos cuando satisfagan correctamente el estándar. No enlaza bibliotecas nativas C/C++ salvo necesidad imprescindible evaluada por separado y aprobada expresamente. Las herramientas genómicas completas y consolidadas, como minimap2, se ejecutan como procesos externos detrás de adaptadores propios. Resolución de DTPD-010; REQ-CONF-69, REQ-DAT-21/22, REQ-MAN-34 a REQ-MAN-36 y REQ-VAL-11; hoja de ruta F1.1 | Reducir la superficie de distribución y acoplamiento sin reimplementar herramientas especializadas ni filtrar decisiones accidentales de una dependencia al dominio | Toda dependencia queda detrás de un contrato propio con tipos y errores del sistema. El adaptador de proceso controla versión, parámetros, comando, entradas y salidas, código de salida, stdout, stderr y fallos. Antes de admitir una dependencia se documentan propósito y requisito, corrección del estándar, mantenimiento, seguridad, licencia, rendimiento, estabilidad de API, características activadas, código nativo o scripts de construcción y coste de distribución. Se minimizan dependencias y *features* y no se añaden por conveniencia anticipada. Las entregas fijan versiones y el manifiesto registra dependencias y herramientas efectivas; DTC-048 fija el toolchain y el lock de los crates aceptados. `noodles-fasta` permanece aprobado por DTC-041, sin prejuzgar otros crates. DTC-047 conserva la validación independiente y DTPD-105, la concreción de su adaptador de proceso externo |
| <a id="dtc-044"></a>DTC-044 | V1 usa errores internos tipados y jerárquicos, separados de una capa pública estable de diagnósticos. Sus categorías mínimas son configuración o esquema, entrada, formato, validación de dominio o biología, recursos, proceso o herramienta externa, E/S o publicación y defecto interno. Resolución de DTPD-011; REQ-CONF-09/14, REQ-FMT-04/07 y REQ-VAL-01 a REQ-VAL-51; hoja de ruta F1.2 | Ofrecer rechazos automatizables y mensajes accionables sin acoplar el contrato público a textos o tipos accidentales de dependencias | Todo diagnóstico público contiene código estable del catálogo versionado `GPS-{categoría}-NNN`, donde la categoría es `CONF`, `INPUT`, `FMT`, `DOMAIN`, `RES`, `TOOL`, `IO` o `INTERNAL`; además contiene categoría, mensaje accionable y contexto disponible: ruta y ubicación, campo e IDs de entidad, evento o artefacto. El catálogo debe enumerar los códigos de F1 antes de implementarla y ningún código se reutiliza ni cambia de significado. La cadena causal técnica queda para detalle o log, no como mensaje principal. Los fallos esperables retornan sin `panic`; los defectos internos se distinguen. Por defecto se renderiza texto humano a stderr. El modo JSON emite exactamente un documento RFC 8259 por invocación, nunca JSONL, validado por esquema y con `diagnostic_kind: "gps_diagnostics"`, `schema_version: "1"` y la lista `diagnostics`. Salidas: 0 éxito o solo avisos; 1 fallo de ejecución o dominio con configuración válida; 2 uso, configuración o entrada inválida; 3 recursos; 4 herramienta externa; 5 E/S o publicación; 70 defecto interno. Los contratos públicos no filtran tipos de crates; DTC-036 rige clase y versión y DTC-048 fija dependencias cuando se elijan |
| <a id="dtc-045"></a>DTC-045 | Toda ejecución iniciada de v1 trabaja en un directorio provisional hermano del destino, en el mismo sistema de ficheros y marcado `INCOMPLETE`. Cada artefacto se cierra, valida e integra primero desde su propio temporal. Un directorio final de éxito o fallo controlado solo aparece cuando sus manifiestos y catálogo de estados son coherentes y se publica mediante una operación Linux atómica no-clobber equivalente a `renameat2(..., RENAME_NOREPLACE)`, detrás de un adaptador propio. Resolución de DTPD-012; REQ-VAL-11, REQ-MAN-73/76 y REQ-SAL-05/19/21 a REQ-SAL-26 | Impedir que una ejecución parcial parezca completa o sobrescriba otra concurrente sin prometer durabilidad física | Se distinguen tres resultados. Un rechazo de admisión o pre-simulación no publica destino: solo deja provisional o diagnóstico externo conforme a REQ-SAL-26. Tras comenzar, un fallo controlado retira o ignora temporales incompletos, conserva y declara exportables los artefactos o conjuntos independientes ya validados, cierra y valida el manifiesto general con estado `fallido`, diagnóstico y catálogo de unidades, y publica atómicamente el directorio final como fallido, nunca completo. Una interrupción no controlable, como SIGKILL o caída, deja el provisional `INCOMPLETE`; conserva diagnóstico solo cuando fue posible escribirlo. Si el destino existe, la operación devuelve `EEXIST` y se rechaza sin sobrescribir, mezclar ni reparar; si el sistema de ficheros no soporta la garantía, se rechaza de forma identificable antes o durante publicación, sin fallback check-then-rename. V1 no reanuda ni reutiliza parciales. Al exportar a otro sistema de ficheros se copia a un provisional allí, se verifican tamaño y digest y se publica con la misma garantía no-clobber; no hay atomicidad entre sistemas. Se prueban concurrencia e interrupciones en WSL2/ext4, Ubuntu nativo y el sistema de ficheros OCI admitido según DTC-042. Publicación lógica y `fsync` o durabilidad física son garantías distintas |
| <a id="dtc-046"></a>DTC-046 | La línea base inicial de F1 usa tres escenarios sintéticos del mismo recorrido: A, aceptación de 100 kbp con configuración mínima representativa, semilla, manifiestos e integridad en cada cambio; B, crecimiento con referencias deterministas de 1, 10 y 100 Mb; y C, observación de 1 Gb manual o programada, no ejecutada en cada cambio. Resolución de DTPD-013; REQ-NOF-01 a REQ-NOF-08; hoja de ruta F0.4 y F1 | Obtener evidencia comparable de corrección y crecimiento sin convertir medidas iniciales en promesas de producto | Se miden tiempo total y por etapa, CPU cuando sea fiable, RSS máximo, bytes de entrada, temporales y salida, MB/s y número de ficheros o artefactos. Se usa una construcción optimizada y portable, sin `target-cpu=native`, configuración y semilla fijas, salida limpia y repeticiones suficientes. Se registran versión, toolchain, plataforma, CPU y memoria, sistema de ficheros, WSL2, Ubuntu 24.04 LTS nativo u OCI y número de repeticiones. WSL2 sirve para tendencias de desarrollo; Ubuntu nativo u OCI es la referencia de publicación. No hay inicialmente objetivos de producto ni fallos por pequeñas variaciones temporales: se conservan línea base y crecimiento observado, mientras corrección e integridad siguen siendo obligatorias. Umbrales y presupuestos permanecen en DTPD-115. Estos tamaños son cargas sintéticas de F1, no una promesa de escala biológica futura |
| <a id="dtc-047"></a>DTC-047 | Los validadores independientes de v1 son dependencias de pruebas, CI y aceptación de publicaciones; no se ejecutan sobre cada fichero de una ejecución normal ni requieren instalación por la persona usuaria. El runtime valida internamente todos los artefactos e invariantes, tamaño e integridad y la publicación. Resolución de DTPD-014; REQ-FMT-04/07, REQ-VAL-01 a REQ-VAL-51 y REQ-SAL-05; hoja de ruta F1 y F9.5 | Obtener un oráculo distinto del productor sin convertir herramientas de prueba en dependencias operativas ni duplicar todo el coste de validación en cada ejecución | En F1 se contrastan FASTA con `samtools faidx`, SHA-256 con `sha256sum`, manifiestos y `effective-config.json` con una implementación JSON Schema distinta, `.tar.gz` mediante listado y extracción estándar e imagen OCI mediante inspección y arranque de prueba. En sus fases, VCF se contrasta con bcftools y BAM/CRAM con samtools u otra referencia independiente. Los casos incluyen válidos, truncados, corruptos, campos desconocidos y errores de borde; CI fija y registra versiones y plataforma. El mismo escritor o lector nunca es el único oráculo. Validez formal y corrección biológica se prueban por separado, esta última con casos conocidos. Una auditoría externa opcional en runtime queda fuera de v1. DTC-048 rige la fijación de versiones y las DTPD de formatos posteriores, sus validadores concretos |
| <a id="dtc-048"></a>DTC-048 | El toolchain oficial inicial es exactamente Rust 1.97.1 stable, verificado como estable actual el 16 de julio de 2026. Al crear F1, `rust-toolchain.toml` fija `channel = "1.97.1"`, `profile = "minimal"`, los componentes `rustfmt` y `clippy` y el target `x86_64-unknown-linux-gnu`; la edición `2024` se declara en cada `Cargo.toml` o en el workspace, no en el fichero de toolchain. `Cargo.lock` se versiona para la aplicación y CI y entregas usan `--locked`. Resolución de DTPD-015; REQ-CONF-111 a REQ-CONF-115, REQ-MAN-05 y REQ-REP-22/23; hoja de ruta F1.1 | Hacer reproducible y auditable la construcción sin declarar una compatibilidad mínima que no se prueba ni trasladar Rust a las personas usuarias del binario o imagen | No existe otra MSRV inicial: el compilador admitido es el fijado. Manifiesto y metadatos registran rustc y Cargo, edición, lock o resolución efectiva y herramientas, sin prometer identidad binaria. Las versiones concretas de crates se eligen al implementar F1 después de comprobar compatibilidad y DTC-043 y quedan fijadas por el lock; una candidata evaluada no cuenta aún como fijada. Toda actualización de Rust o dependencias es deliberada y aislada y exige pruebas, validadores independientes, vectores de RNG y reproducibilidad y escenarios de rendimiento; después se vuelven a declarar garantías sin asumir equivalencia. Rust 1.97.1 sustituye la propuesta accidental basada en Rust 1.94.0 local. Las personas usuarias de binario o imagen no necesitan Rust |
| <a id="dtc-049"></a>DTC-049 | Los valores numéricos configurables de probabilidad o frecuencia y sus tolerancias generales se representan como `f64`. Todo valor y tolerancia debe ser finito; una probabilidad o frecuencia pertenece al intervalo cerrado `[0,1]`. La tolerancia absoluta general es configurable en el intervalo cerrado `[0, 1e-6]` y vale `1e-9` por defecto. Resolución parcial de DTPD-110; REQ-CONF-106/107, REQ-VAL-40, REQ-MAN-81 y ACE-105 | Obtener una validación y normalización únicas, auditables y reproducibles sin convertir el orden físico o el paralelismo en parte accidental del resultado | Los valores se recorren en orden canónico estable por identidad y categoría y se suman con Neumaier, sin reducción paralela. Se calcula `delta = abs(suma - 1)` y el borde inclusivo acepta `delta <= tolerancia`. Si se acepta y la suma es positiva, finita y normalizable, cada valor se divide por la suma y se recomputan suma y delta con el mismo algoritmo y orden; si excede la tolerancia, la suma no es positiva o finita, la normalización no produce valores válidos o una entrada es inválida, se rechaza antes de simular. Configuración efectiva y manifiesto registran tolerancia efectiva, orden canónico, algoritmo y versión, suma y delta, y valores antes y después. Toda serialización decimal debe permitir recuperar exactamente el mismo `f64`; las garantías de reproducibilidad exigen conservar estas condiciones. Es una validación determinista de configuración: no altera DTC-024, no sustituye `0,5/C` ni ninguna tolerancia específica y no selecciona crate |

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
evaluación revisable considera `noodles-fasta` 0.66.0 candidata compatible, con
características mínimas síncronas, pero no la fija todavía en `Cargo.lock`; el
grafo final deberá revalidarse. La política queda en [DTC-043](#dtc-043) y el
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

## 7. Decisiones pendientes para fases posteriores

| ID | Decisión | Opciones | Recomendación razonada | Consecuencias | Bloquea |
|---|---|---|---|---|---|
| DTPD-101 | Representación interna de coordenadas | Base 1 cerrada; base 0 semiabierta; valor etiquetado con varias representaciones | Evaluar operaciones de segmentos y evitar conversiones implícitas; mantener tipos de frontera aunque se elija una canónica | Afecta errores de borde, normalización y compatibilidad con formatos | F2.1/F2.2 |
| DTPD-102 | Estado haplotípico y fichero obligatorio de verdad | Secuencias completas; lista de variantes; segmentos persistentes; modelo híbrido | Probar reconstrucción, memoria y composición con SNP, indels y SV antes de elegir; separar modelo interno y formato de intercambio | Es el contrato central de fase, dosis y linaje | F2.1/F2.4 |
| DTPD-103 | Representación técnica de VCF y complemento | Versión oficial VCF; normalizador propio o externo; nombres de campos `INFO`; complemento tabular, documental o binario | Seleccionar por capacidad de poliploidía, BND, repeticiones, `N`, BGZF/índices e interoperabilidad; demostrar validez independiente del VCF y reconstrucción sin pérdida, sin emitir SNP con `REF=N` ni inventar bases | Afecta nombres de campos, cabecera, esquema complementario y validadores, no el reparto funcional de verdad | F2.5 y F5.5 |
| DTPD-104 | Adaptadores SAM/BAM/CRAM e índices | Crate; biblioteca nativa; herramientas externas; combinación | Comparar cumplimiento de estándar, CRAM con referencia exacta, ordenación, índices, cabeceras y distribución | Determina E/S, procesos, errores y rendimiento | F3.3 |
| DTPD-105 | Concreción del adaptador de proceso externo para minimap2 y otros mapeadores | Descubrimiento explícito o ruta configurada; política de invocación y sandbox; compatibilidad de versiones; clasificación de códigos, señales y fallos | Partir del proceso externo ya decidido por DTC-043 y concretar descubrimiento, comando reproducible, aislamiento, compatibilidad y diagnóstico sin ocultar versión, parámetros, stdout, stderr ni código de salida | Afecta instalación, seguridad, logs y reproducibilidad; no reabre servicio frente a proceso | F3.2 |
| DTPD-106 | Perfiles WGS cuantitativos incorporados | Perfiles versionados de Illumina y PacBio HiFi; sustituciones explícitas permitidas | Validar cada perfil con datos adecuados y distinguir claramente valores incorporados y sustituciones | Fija distribuciones de longitud, error, calidad y artefactos; no incluye estimación ni perfiles empíricos | F4.1–F4.3 |
| DTPD-107 | Calibración empírica (retirada de v1) | Sin elección técnica en v1 | La decisión funcional DTC-023 aplaza la capacidad completa; cualquier diseño del estimador y del artefacto corresponde a DTF-016 | No bloquea ninguna fase de v1 | — |
| DTPD-108 | Verdad y procedencia de lecturas | Registro por lectura; registros por bloques; combinación con agregados | Mantener para cada lectura o pareja final el origen y la diferencia verdadera/observada, y agregados solo para descartes o no asignadas; medir volumen antes de fijar formato, sin cambiar el contrato funcional | Puede dominar el almacenamiento y la validación | F4.4/F4.5 |
| DTPD-109 | Fórmulas de F diploide y Ne reproductivo realizado | Fórmulas publicadas aplicables a la heterocigosidad diploide y a las contribuciones observadas; estado no evaluable cuando corresponda | Evaluar contra el significado funcional confirmado y casos simulados; no adaptar F a P=1 o P>2 ni introducir un objetivo u optimización de Ne | Cambia estadísticas y pruebas, no los descendientes ni el genotipo; F fuera de diploides conserva «no aplicable» | F6.3 y F7.3 |
| DTPD-110 | Distribuciones y tolerancias poblacionales posteriores | Catálogos de distribuciones aún no elegidos; tolerancias específicas de modelos poblacionales o estadísticos posteriores | Aplicar [DTC-049](#dtc-049) a toda probabilidad o frecuencia configurable y decidir aquí solo las distribuciones y tolerancias específicas cuando su fase defina la semántica; ninguna puede sustituir `0,5/C` | No reabre `f64`, suma, normalización ni tolerancia general ya resueltos; afecta modelos poblacionales posteriores y sus pruebas | F6.1/F6.2 |
| DTPD-111 | Mapas, posiciones y comparación estructural para recombinación | Formatos tabulares candidatos; generadores de Poisson; interpolación de mapas; representación canónica comparable | Evaluar cobertura completa, monotonicidad, conversiones y equivalencia de segmentos, orden, orientación, copias y puntos de ruptura sin exigir el mismo ID histórico | Afecta bivalentes, linaje y casos estructurales; no puede introducir compatibilidad parcial | F7.1 |
| DTPD-112 | Representación de grafos y reglas declarativas | Estructuras anidadas; tablas de nodos/aristas; lenguaje declarativo limitado sin código | Elegir la sintaxis más validable que represente Y/O/NO, ámbitos y dependencias sin ejecutar código | Afecta mensajes de ciclo, migración, fundación y familias | F7.5 y F8.1 |
| DTPD-113 | Paralelismo reproducible | Por referencia, individuo, evento, región o tareas; ejecución secuencial de referencia | Posponer la elección hasta perfilar y conservar ámbitos aleatorios y orden lógico independientes del planificador | Puede cambiar uso de memoria y orden de salida; no debe cambiar verdad | F9.3 |
| DTPD-114 | Formato del paquete exportado | Directorio autocontenido; archivo contenedor; formatos archivados comprimidos | Elegir por interoperabilidad, preservación de permisos/rutas e integridad antes/después, con eventos grandes | Afecta portabilidad y almacenamiento duplicado | F9.6 |
| DTPD-115 | Objetivos y presupuestos cuantitativos | Umbrales por escenario; presupuesto configurable; perfiles de recursos | Fijarlos después de medir escenarios representativos; todo presupuesto debe incluir el mínimo, la trazabilidad única y las salidas activadas, y su exceso rechaza sin degradación automática | Define expectativas de escala y criterios de publicación | F9.1/F9.4 |

## 8. Decisiones futuras que no bloquean v1

Esta tabla es un resumen derivado. El registro canónico con fuentes, estado,
destino y dependencias es
[registro-ampliaciones-futuras-0.1.md](./registro-ampliaciones-futuras-0.1.md).

Las entradas DTF-001 a DTF-036 se definen una sola vez en ese registro. Este
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

Las decisiones DTPD-101 a DTPD-115 deben resolverse antes de planificar el hito
que señalan. No bloquean la redacción del plan de F1 salvo que su evaluación
revele una incompatibilidad con un contrato base de configuración, identidad,
referencia o reproducibilidad.
