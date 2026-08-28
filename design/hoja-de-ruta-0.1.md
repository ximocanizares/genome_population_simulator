# Hoja de ruta del simulador de genomas y poblaciones

- **Versión:** 0.2
- **Estado:** vigente; F0 aprobada, F1 cerrada en `e14705a`/`efca265`; F2
  simplificada por DTC-079 y pendiente de revisión
- **Fecha:** 21 de agosto de 2026
- **Especificación de origen:**
  [especificaciones-0.1.md](./especificaciones-0.1.md)
- **Lenguaje del núcleo:** Rust

Esta hoja es un instrumento derivado de planificación y no crea obligaciones
funcionales distintas de la especificación. La asignación exhaustiva de REQ,
ACE, fase y prueba prevista está en
[trazabilidad-requisitos-0.1.md](./trazabilidad-requisitos-0.1.md).
Las capacidades excluidas de v1 y su destino post-v1 se mantienen en
[registro-ampliaciones-futuras-0.1.md](./registro-ampliaciones-futuras-0.1.md).

## 1. Finalidad

Esta hoja de ruta divide el desarrollo en fases que produzcan resultados
verificables y mantengan controlado el crecimiento de la aplicación. No es
todavía un plan de implementación: no decide módulos, crates, formatos de
configuración, algoritmos concretos ni archivos de código.

Cada fase deberá transformarse en un plan técnico detallado únicamente cuando
la arquitectura y las decisiones de las que dependa estén acordadas. El plan de
una fase dividirá sus hitos en cambios pequeños, con pruebas y pausas de revisión.

## 2. Principios para controlar el crecimiento

1. **Mantener un recorrido completo.** Desde las primeras fases debe existir un
   camino pequeño pero funcional desde la configuración hasta resultados que
   puedan validarse.
2. **Ampliar una dimensión cada vez.** Añadir variantes, modelos poblacionales o
   realismo técnico sin introducir simultáneamente varias capacidades no
   relacionadas.
3. **Cerrar cada fase antes de apoyar otra sobre ella.** Una fase termina con
   pruebas, evidencias y una revisión explícita.
4. **No anticipar trabajo futuro.** Una interfaz podrá dejar una extensión
   prevista, pero no se implementará una abstracción sin un caso confirmado.
5. **Conservar la verdad conocida.** Toda ampliación debe mantener referencia,
   coordenadas, ploidía, fase, dosis, linaje y procedencia sin pérdidas.
6. **Separar verdad biológica y evidencia técnica.** FASTQ, BAM, CRAM, errores y
   contaminación nunca modificarán el genotipo verdadero.
7. **Medir pronto, optimizar con evidencia.** Cada fase observará tiempo y
   memoria en casos representativos. Las optimizaciones complejas requerirán
   pruebas de rendimiento reproducibles.
8. **Planificar solo la siguiente fase.** No se redactará de una vez el plan
   detallado de toda la aplicación, porque las decisiones posteriores dependen
   de lo aprendido en fases anteriores.
9. **Mantener paridad Marimo.** Desde F2.0, ninguna capacidad configurable de
   F2–F9 cierra hasta que la aplicación permita configurarla, validarla y
   lanzarla mediante el mismo conjunto TOML y binario y supere equivalencia y
   accesibilidad.

## 3. Niveles de planificación

- **Fase:** capacidad amplia del producto y resultado que permite avanzar.
- **Hito:** resultado verificable dentro de una fase.
- **Paso:** cambio pequeño de implementación, revisión y prueba independiente.

Esta hoja de ruta define fases e hitos. Los pasos se escribirán en los planes de
implementación posteriores.

## 4. Dependencias generales

```text
F0 Arquitectura
 └─ F1 Núcleo reproducible
     └─ F2 Verdad genómica mínima
         └─ F3 Recorrido técnico mínimo
             └─ F4 Realismo técnico de secuenciación
                 └─ F5 Variación y referencias complejas
                     └─ F6 Población inicial
                         └─ F7 Evolución poblacional
                             └─ F8 Familias y segregantes
                                 └─ F9 Escalabilidad y consolidación
```

El orden expresa la secuencia recomendada. Una tarea aislada de preparación
podrá adelantarse si no introduce contratos ni decisiones prematuras, pero una
fase no se considerará terminada hasta satisfacer todas sus dependencias.

## 5. Fases e hitos

### F0. Arquitectura y decisiones técnicas

**Estado:** cerrada y aprobada el 8 de agosto de 2026. Queda autorizada la
ejecución del plan detallado de F1.

**Objetivo:** convertir la especificación funcional en una arquitectura Rust
trazable y resolver las decisiones necesarias para planificar el primer código.

**Depende de:** especificación 0.1 y decisión de usar Rust.

**Hitos:**

- **F0.1. Arquitectura del sistema.** Definir componentes, responsabilidades,
  flujos, contratos, errores y dirección de dependencias.
- **F0.2. Modelo de dominio.** Definir la representación conceptual de
  referencias, coordenadas, segmentos, haplotipos, individuos, poblaciones,
  eventos, lecturas y resultados.
- **F0.3. Decisiones técnicas iniciales.** Elegir mediante evaluación las
  opciones necesarias para configuración, manifiestos, aleatoriedad,
  serialización, formatos genómicos, procesos externos y pruebas.
- **F0.4. Estrategia de rendimiento.** Definir escenarios pequeños y grandes,
  métricas, límites observables y pruebas iniciales de memoria y tiempo.
- **F0.5. Plan de la fase F1.** Dividir el núcleo reproducible en pasos pequeños
  con dependencias y pruebas concretas.

**Resultado de salida:** arquitectura 0.1, decisiones técnicas justificadas y
plan ejecutable de F1.

**Evidencia de cierre:** todos los requisitos necesarios para F1 tienen un
componente responsable y no quedan decisiones de diseño ocultas en su plan.

> Revisión de F0: aprobar la arquitectura y el plan de F1 antes de crear el
> núcleo de producción.

### F1. Núcleo reproducible

**Objetivo:** obtener una aplicación Rust mínima que valide una configuración,
cree una ejecución identificada y produzca resultados trazables y repetibles.

**Estado:** replanificación aprobada tras la auditoría sistemática del commit
`66ad8b2`. La implementación anterior queda como checkpoint recuperable, no
como base que deba conservarse íntegramente.

**Depende de:** F0.

**Hitos:**

- **F1.1. Base del proyecto.** Disponer del proyecto Cargo, comprobaciones de
  formato, análisis estático, pruebas y política inicial de dependencias.
- **F1.2. Configuración validada.** Leer la configuración mínima y resolver solo
  los valores y selectores consumidos por la referencia inicial o reducida.
  Las matrices poblacionales, temporales y espaciales se incorporan con su
  primer consumidor real en F2/F6/F7 y se auditan conjuntamente en F9.7.
- **F1.3. Infraestructura de aleatoriedad reproducible.** Registrar la semilla,
  validar los vectores normativos, aislar flujos locales y preparar comparadores
  sin presentar valores sintéticos como una ejecución o un artefacto real.
- **F1.4. Identidad y resultados.** Crear los identificadores realmente
  producidos, el provisional privado, la configuración efectiva y un
  manifiesto general mínimo con tamaños/SHA-256 y publicación atómica
  no-clobber. Los manifiestos de evento y conjuntos aparecen con esas unidades
  en F3/F4; el hardening contra atacante local activo se evalúa en F9.
- **F1.5. Referencia y repetición mínimas.** Aceptar una referencia pequeña
  simulada o suministrada, conservar su FASTA e identidad verificable y repetir
  el recorrido ejecutable con semilla, configuración, entrada, versión y
  manifiesto mínimos idénticos.
- **F1.6. Referencia reducida reproducible.** Derivar de un FASTA real validado
  una referencia de tamaño exacto mediante regiones aleatorias no solapadas,
  conservar orden/orientación y `N`, y publicar FASTA, mapa, integridad,
  subsemilla y manifiesto reconstruibles. Las cuestiones funcionales y
  DTPD-101/DTPD-116 están cerradas mediante DTC-050/DTC-056/DTC-057. El código
  útil del checkpoint puede reutilizarse, pero el hito se vuelve a validar tras
  retirar `ResolutionDomain`, `resolution-evidence`, eventos sintéticos y el
  hardening físico que DTC-078 saca de F1.

**Resultado de salida:** ejecución mínima repetible con configuración, semilla,
referencia completa o reducida y manifiesto; una reducción incluye su mapa.

**Evidencia de cierre:** dos ejecuciones equivalentes producen el mismo
resultado lógico y una configuración inválida no deja una ejecución marcada
como completa. Una reducción se reconstruye exactamente desde origen+mapa y
conserva tamaño, orden, clase, orientación y `N`.

> Revisión de F1: comprobar reproducibilidad, diagnósticos y trazabilidad antes
> de representar individuos.

### F2. Verdad genómica mínima

**Objetivo:** crear individuos directamente configurados y conservar una verdad
genotípica faseada respecto a una única referencia.

**Depende de:** F1.

**Alcance inicial deliberado:** referencia pequeña, una ploidía uniforme,
individuos directos, SNP e inserciones o deleciones pequeñas no ambiguas.

**Hitos:**

- **F2.0. Aplicación Marimo base.** Exponer a Marimo únicamente la configuración
  F1 ya implementada; ploidía, población, variantes y `scenario_id = "direct"`
  aparecen después con su primer consumidor F2. F2 no crea documentos de
  escenario. Incorporar un entorno `uv` bloqueado, un
  lanzador en `127.0.0.1` y controles para editar, importar y exportar un único
  TOML general. Invocar `validate`/`run` mediante argumentos separados y mostrar
  diagnósticos, progreso y resúmenes. Probar instalación desde el repositorio,
  arranque local sin red y paridad UI→TOML→Rust. El aprovisionamiento inicial
  explícito con `uv sync --frozen` puede descargar lo resuelto por `uv.lock`;
  el uso normal no accede a la red. DTC-079 excluye protocolos paralelos,
  almacenes de sesión, defensas contra atacante local y expedientes por paso;
  la instalación distribuible sin descarga previa, `.tar.gz`, OCI, plataforma
  oficial y endurecimiento se cierran en F9.
- **F2.1. Estado haplotípico.** Representar copias cromosómicas, fase, dosis y
  procedencia mediante DTC-056 y el modelo híbrido DTC-061, sin duplicar
  innecesariamente secuencias completas; conservar `N` y preparar la proyección
  obligatoria `truth.jsonl` sin acoplarla a tipos Rust.
- **F2.2. Variantes pequeñas y matriz SNV.** Aplicar SNP explícitos e *indels*
  exactos o aleatorios; preparar, admitir, generar y reutilizar el par canónico
  TSV+JSON de DTC-054 ligado por SHA-256 a la referencia. Generar candidatos
  bialélicos sin reemplazo desde número exacto/media/densidad, espacio uniforme
  o ponderado, 12 sustituciones y espectro `alt_af` no solapado en modo
  estocástico o exacto por mayores restos. `poisson-count-v1` de DTC-067 solo
  rige recuentos desde una media, no alelos. Regenerar y validar mediante
  DTC-068 los ejemplos compactos `human_1000g_global_sfs_approx_v1` y
  `human_harris2017_snv12_approx_v1`, incluida fuente eLife fijada por
  tamaño/SHA-256, receta, CC BY y digests de salida. No procesar el VCF Phase 3
  ni depender de `bcftools` para estos builtins; nunca seleccionarlos en runtime
  de forma implícita. Rechazar SNV sobre `N`, VCF o
  semánticas de frecuencia alternativas; conservar `N` en intervalos que no
  necesitan interpretarla. Todos los intervalos de frecuencia son `[min,max)`;
  0/1 exactos solo proceden de `fixed`. DTC-071 fija los once ámbitos RNG F2 sin circularidad,
  DTC-072 sus consumidores no Poisson y DTC-075 separa selección SNV estática
  de selección indel dinámica con reenumeración tras cada evento. REQ-CONF-08/21 aportan aquí evidencia
  parcial y cierran en F5.2/F7.2 según su alcance. DTPD-118 está resuelta por
  DTC-054.
- **F2.3. Individuos directos desde población de origen.** Construir individuos
  asignados a una única población y muestrear la matriz mediante Hardy-Weinberg
  para P=1..10 o `F` solo diploide, con fase y subsemillas independientes. La
  distribución genotípica explícita no tiene fichero, esquema ni control Marimo
  en F2; REQ-CONF-130 aporta evidencia parcial y cierra en F6.2. Integrar restricciones
  explícitas, incluidas las dos autorizaciones limitadas y trazadas para
  genotipos incompatibles o matriz activa vacía; integrar los demás cambios y
  reconstruir cada haplotipo. REQ-IND-02/03 aportan aquí evidencia parcial sobre
  una referencia estructural por individuo y cierran en F5.2. Activar preparación/generación/reutilización,
  perfiles, población y modelo en la pestaña `Individuos`.
- **F2.4. Verdad estructurada.** Publicar `truth.jsonl` v1 obligatorio, UTF-8
  sin comprimir, con cabecera, copias/genotipos tipados, orden e integridad;
  reconstruirlo y contrastar genotipos/haplotipos independientemente del VCF.
- **F2.5. VCF de verdad mínimo.** Producir exactamente un VCF 4.5 multimuesta
  estándar por conjunto, ordenado, derivado de `truth.jsonl` y sin campos
  inferidos de lecturas. Publicar por defecto `.vcf.gz` BGZF con exactamente un
  índice TBI o CSI según el límite de TBI, o `.vcf` de texto sin índice cuando
  se solicite explícitamente; validar y publicar VCF+índice de forma atómica.
  Normalizar de forma determinista/versionada desde los haplotipos finales y
  agrupar en F2 todos los alelos compatibles en un único registro multialélico
  por `CHROM+POS`, sin descomposición automática. Todo registro usa `GPS_CASE`;
  solo una pérdida añade `GPS_LOSS`, `GPS_COMPLEMENT` y una entrada por caso en
  el `vcf-complement.jsonl` condicional del conjunto. No emitir SNP con `REF=N`
  ni inventar una base. Aceptar en CI con bcftools/HTSlib-tabix 1.24 y un
  oráculo semántico independiente sobre texto, BGZF+TBI y BGZF+CSI.

**Resultado de salida:** referencia, individuos faseados, `truth.jsonl` y VCF
4.5 verdadero coherentes.

**Evidencia de cierre:** casos sintéticos independientes reconstruyen las
secuencias esperadas y demuestran coordenadas, alelos, fase y dosis. Incluyen un
SNP rechazado sobre `N`, una inserción aleatoria sin `N`, una explícita que la
conserva y la ausencia de `vcf-complement.jsonl` para todos los casos F2
representables, sin resolver `N` como A/C/G/T. F2 no fabrica una pérdida. Deja evidencia
parcial de REQ-SAL-02 en F2.4.1/F2.4.3/F2.4.4 y de REQ-SAL-12 en
F2.5.4/F2.5.5/F2.5.6; sus cierres completos siguen en F8.5 y F5.5,
respectivamente.

> Revisión de F2: validar la verdad biológica antes de generar una sola lectura.

### F3. Recorrido técnico mínimo

**Objetivo:** producir el primer conjunto pequeño útil para probar externamente
un detector de SNP.

**Depende de:** F2.

**Alcance inicial deliberado:** lecturas cortas tipo Illumina, un único perfil,
un evento, sin contaminación ni degradaciones avanzadas y minimap2 como
mapeador.

**Hitos:**

- **F3.1. Lecturas desde haplotipos.** Generar FASTQ con cobertura o número de
  lecturas y contribución proporcional a dosis y ploidía. Las regiones con `N`
  son muestreables y la lectura la conserva con Phred 0.
- **F3.5. Conjunto exportable mínimo.** Reunir referencia, verdad, FASTQ,
  configuración y manifiesto, y comprobar el origen de las lecturas sin inferir
  genotipos ni ejecutar programas detectores de SNP.

El mapeo con minimap2 u otra herramienta y la creación de BAM o CRAM se realizan
fuera del simulador a partir de los FASTQ exportados; F3 no los invoca ni los
valida.

**Resultado de salida:** recorrido completo:

```text
configuración → referencia → individuo → verdad → FASTQ
```

**Evidencia de cierre:** el FASTQ y sus metadatos se corresponden con una verdad
biológica conocida; pueden entregarse a un mapeador externo sin que este forme
parte del simulador.

> Revisión de F3: usar el primer conjunto extremo a extremo antes de ampliar el
> modelo biológico o técnico.

### F4. Realismo técnico de secuenciación

**Objetivo:** incorporar pronto los perfiles, degradaciones y eventos técnicos
de WGS confirmados, sin alterar la verdad genotípica.

**Depende de:** F3.

**Alcance inicial deliberado:** usar individuos directos y variantes pequeñas
de F2. Las fases biológicas posteriores deberán demostrar que sus individuos y
poblaciones pueden alimentar este mismo recorrido sin crear otra implementación.

**Hitos:**

- **F4.1. Perfiles tecnológicos.** Completar Illumina *paired-end* para una
  población y un evento técnico únicos. PacBio HiFi se aplaza a DTF-049 y los
  recorridos con varias poblaciones o eventos a DTF-050.
- **F4.2. Cobertura, longitud, error y calidad.** Aplicar perfiles validados y
  sustituciones configurables por los ámbitos permitidos.
- **F4.3. Artefactos WGS.** Simular quimeras y rechazar configuraciones de
  duplicados, adaptadores, sesgo GC no neutro y anomalías *paired-end*, que se
  aplazan en DTF-054–056. El contenido GC
  excluye `N`; una ventana sin bases conocidas no recibe corrección y se registra.
- **F4.4. Contaminación biológica.** Sustituir lecturas según fuentes, unidades y
  matrices dirigidas compatibles.
- **F4.5. Demultiplexado técnico.** Representar reasignaciones y descartes por
  índices y conservar solo sus recuentos agregados, sin duplicar secuencias.
- **F4.6. Alineamientos y BAM/CRAM.** Aplazados en DTF-055. F4 termina en FASTQ
  final; el modo directo permanece además en DTF-051.
- **F4.7 (retirado de v1). Calibración desde datos reales.** No se implementa en
  esta versión. La admisión rechaza antes de simular FASTQ o referencia de
  calibración, solicitudes de estimación y perfiles empíricos reutilizables.

**Resultado de salida:** un evento WGS Illumina configurable y trazable para
una población, disponible desde una etapa temprana del proyecto.

**Evidencia de cierre:** validaciones internas cuantifican lo solicitado y lo
realizado, y demuestran que las degradaciones solo afectan a la evidencia
técnica.

> Revisión de F4: aprobar los perfiles y artefactos antes de ampliar las clases
> de variación y los modelos poblacionales.

La calibración empírica automática queda en la hoja de trabajo futura, sin
bloquear F4. Antes de recuperarla deberán acordarse el contrato mínimo de
entrada, la separación entre variantes y errores, el esquema reutilizable, las
métricas y umbrales de fiabilidad y las reglas de reutilización. Los perfiles
incorporados de Illumina y sus sustituciones explícitas sí forman parte de F4.
DTF-049 reserva PacBio HiFi. DTF-016 reserva también la materialización exacta de VCF
poblacionales y perfiles empíricos exactos; no sustituye ni bloquea los dos
ejemplos compactos aproximados de F2.2.

### F5. Referencias y variación complejas

**Objetivo:** ampliar la verdad genómica manteniendo reconstrucción, fase y
separación estricta entre referencias.

**Depende de:** F4.

**Hitos:**

- **F5.1. Referencias estructurales.** Crear referencias derivadas mediante
  duplicaciones, deleciones, inversiones y translocaciones permitidas.
- **F5.2. Variantes estructurales individuales.** Representar cambios
  estructurales en haplotipos y conservar sus efectos finales, incluidas las `N`
  internas de intervalos copiados, eliminados, invertidos o desplazados.
- **F5.3. Ploidía general.** Admitir la ploidía confirmada por cromosoma o
  *scaffold* sin diferencias entre individuos.
- **F5.4. Repeticiones en tándem.** Añadir repeticiones en tándem dentro del
  alcance y su modelo escalonado simple de expansión o contracción de una unidad
  completa por mutación. Los elementos transponibles no son una clase funcional
  de v1.
- **F5.5. Variantes solapadas y complemento.** Resolver cambios secuenciales,
  variantes eliminadas y alelos finales complejos. Cuando un caso pierda detalle
  en VCF, enlazar sus registros con una única entrada complementaria y demostrar
  la reconstrucción completa del par. Esta fase aporta la primera prueba
  positiva obligatoria VCF+complemento mediante un caso complejo/SV/BND real.
  **Subpaso diferido F5.5.a — variación aleatoria directa por individuo:**
  incorporar en la creación inicial directa las cantidades medias o densidades
  globales y sus sustituciones por individuo, distribuciones y regiones, y la
  exclusión explícita respecto de la matriz SNV canónica; demostrar que cada
  sustitución afecta solo a su individuo. Este subpaso conserva separados los
  cambios *de novo* posteriores a la inicialización (F7.2.a) y cubre
  REQ-CONF-20, REQ-CONF-22, REQ-IND-15, REQ-DIR-02 y la parte inicial de
  REQ-IND-64.
- **F5.6. Varios conjuntos de referencia.** Mapear las mismas muestras contra
  referencias diferentes, producir BAM/CRAM y VCF independientes y demostrar
  que el manifiesto del evento y el manifiesto raíz enlazan todos los resultados
  reales con su evento, conjunto, referencia y artefacto conforme a REQ-EVO-78.
  **Subpaso diferido F5.6.a — matrices SNV por referencia:** permitir y validar
  enlaces de matriz canónica separados para cada referencia estructural, su
  selección trazable y su publicación sin mezclar matrices entre referencias;
  cierra la parte multirreferencia de REQ-CONF-128.

**Resultado de salida:** individuos complejos reconstruibles y conjuntos
independientes por referencia.

**Evidencia de cierre:** pruebas con verdad conocida demuestran transformaciones,
límites, solapamientos y ausencia de mezcla entre referencias, incluida una SV
que contiene `N` y la conserva donde corresponda.

> Revisión de F5: confirmar que la complejidad estructural no rompe la verdad
> mínima ni el realismo técnico establecido en F4.

### F6. Población inicial

**Objetivo:** crear poblaciones iniciales reproducibles sin ejecutar todavía su
evolución a lo largo de generaciones.

**Depende de:** F5.

**Hitos:**

- **F6.1. Creación probabilística y controlada.** Generar alelos, genotipos,
  dosis y recuentos según el modo seleccionado.
- **F6.2. Frecuencias y Hardy-Weinberg.** Aplicar y comprobar la distribución
  binomial bialélica y la multinomial multialélica en ploidías 1–10, tanto en
  modo probabilístico como controlado, reutilizando el muestreo matricial
  probabilístico de F2.3 y extendiéndolo a subpoblaciones, recuentos controlados
  y distribución genotípica explícita. Aquí cierran REQ-CONF-130 y
  REQ-EVO-35/36/49.
- **F6.3. Endogamia inicial.** Aplicar las fórmulas confirmadas de F automático
  solo en diploides, probar F=0, valores intermedios y F=1, rechazar F en P=1 o
  P>2 y registrar F realizado como «no aplicable» fuera de diploides, sin crear
  una segunda implementación distinta de la usada para matrices en F2.3.
- **F6.4. Paneles haplotípicos y desequilibrio de ligamiento.** Admitir un panel
  canónico de cinco ficheros, preparado por una única orden Rust, con bloques de
  LD no solapados, SNP A/C/G/T, haplotipos locales únicos y frecuencias densas por
  población y bloque. Exigir una única referencia: cada ID de bloque coincide
  con el ID del metadato, el SHA-256 recalculado del FASTA coincide con su digest
  declarado y cada secuencia existe en ese FASTA. En construcción intacta, tratar cada bloque como locus
  multialélico y reutilizar HWE, F diploide o distribución explícita y los modos
  probabilístico o controlado comunes. Crear combinaciones intrabloque intactas sin recombinación,
  ligarlas a copias biológicas y conservar el HAP fundador aunque mutaciones
  posteriores formen otra capa. Separar `haplotype_construction`
  intacto/recombinante de `frequency_sampling`
  probabilístico/controlado y resolver el primero por población×bloque conforme
  a DTC-033. Admitir ambos muestreos y HWE/F/distribución explícita en intacto;
  en recombinante, solo muestreo probabilístico y HWE o F diploide. En esta alternativa, ejecutar productores directamente por
  individuo×bloque: en diploides Bernoulli(F) elige uno duplicado o dos
  independientes; en otras ploidías, P productores independientes. Cada
  productor usa solo fuentes del panel original, sortea las uniones y devuelve
  primaria intacta o mosaico; interpretar x como unión x|x+1, comenzar en la
  primaria y alternar segmentos *forward* sin ajustar a SNP, conservando cruces
  invisibles. Implementar en F6.4 la variante uniforme y el constructor común de
  segmentos cerrados por DTC-055; F7.1 conectará mapas y meiosis al mismo
  constructor tras DTPD-111, que no bloquea F6.4. Agrupar las C copias finales
  por ID de contenido y calcular `n_h/C`, distinguiendo productores con cruce de
  copias derivadas.
  Incorporar en Marimo la sección poblacional para preparar/validar el panel,
  mapear poblaciones, editar la tabla por bloque, previsualizar métricas y enlazar
  las tres salidas obligatorias sin cargar artefactos grandes. Publicar estas
  salidas con las cabeceras, claves, orden y serialización TSV exactos de
  DTC-055, dejando clase/versión en sus enlaces de manifiesto.
  Permitir HAP distintos entre bloques y obtener los SNP exteriores desde
  `alt_af`. Rechazar cualquier locus de matriz dentro de un bloque, sin prioridad
  ni sobrescritura. Consumir DTC-055, que cierra DTPD-119.

  **Encaje de planificación pendiente para generación automática.** La
  generación por el programa de bloques, SNP, haplotipos y frecuencias que
  produzcan ese mismo paquete se ubica como preparación previa dentro de F6.4,
  no en F2. Todavía no amplía el alcance ejecutable de este hito: la
  especificación vigente solo obliga a admitir un panel suministrado y faltan
  sus REQ/ACE, arquitectura y decisión técnica. Las decisiones cerradas, las
  puertas y los bloqueos se mantienen en el
  [preplan de generación automática de bloques de LD](./preplan-generacion-automatica-bloques-ld-0.1.md).
  La preparación automática producirá el mismo paquete canónico de cinco
  ficheros sin ampliarlo y enlazará desde su manifiesto cuatro informes TSV
  externos de geometría, frecuencias, LD e intentos. La publicación opcional de
  todas las parejas de `r²` queda reservada a DTF-047.
  Cada perfil reutilizable de entrada será, tanto si está incorporado como si lo
  suministra el usuario, un paquete validado de exactamente cuatro ficheros:
  `block_lengths.tsv`, `snp_gaps.tsv`, `ld_profile.tsv` y
  `profile.metadata.json`. La máscara de exclusión ligada a la referencia se
  mantendrá fuera y se declarará y verificará por separado.
  `block_lengths.tsv` tendrá la cabecera exacta
  `length_min_bp<TAB>length_max_bp<TAB>weight`, serialización TSV canónica y
  filas ordenadas, únicas y no solapadas; validará la suma según DTC-049 sin
  normalizarla automáticamente.
  `snp_gaps.tsv` aplicará el mismo contrato con cabecera exacta
  `gap_min_bp<TAB>gap_max_bp<TAB>weight`; su tabla 1–5000 documentada es solo un
  ejemplo ilustrativo y nunca un valor predeterminado.
  `ld_profile.tsv` tendrá una fila canónica por población e intervalo, con
  `population_profile_id`, límites, mínimo evaluable, tolerancia absoluta y los
  cuantiles p10/p25/p50/p75/p90; todas las poblaciones compartirán exactamente
  las mismas fronteras contiguas desde 1. Los valores CEU/YRI documentados serán
  ilustrativos.
  `profile.metadata.json` será un objeto cerrado que identificará el recurso,
  su origen builtin o supplied, su derivación aproximada o empírica, el modelo
  de cuantiles y la procedencia obligatoria, y enlazará exactamente los tres TSV
  por nombre canónico, tamaño y SHA-256 recalculado, sin duplicar sus tablas.
  La máscara opcional `excluded_regions.tsv` permanecerá fuera del paquete,
  declarará intervalos cerrados base 1 no solapados sobre secuencias del FASTA y
  quedará enlazada por SHA-256 junto con la identidad de la referencia. Su
  ausencia no desactivará la exclusión automática de `N`, loci reservados o
  bloques aceptados.
  La configuración usará `[initial_population.ld_blocks]` con
  `source = "generated"`, una única fuente `profile|profile_path`, una única
  medida `block_count|covered_bases|coverage_fraction` y máscara opcional. Cada
  población simulada se mapeará en su subsección a un
  `population_profile_id` existente. Las sustituciones DTC-033 solo cambiarán
  ese mapeo y nunca parchearán tablas o parámetros del builtin; la geometría
  común se resolverá una sola vez.
  La geometría será común e inmutable. Frecuencias, cuantización, inicialización,
  cálculo de cuantiles Type 7, puntuación y recocido se resolverán
  independientemente por población, con `population_attempt_index` propio. El score poblacional minimizará, en
  orden canónico y con suma Neumaier, los excesos cuadrados de los cinco
  cuantiles sobre sus tolerancias para sus intervalos activos. Solo aceptará
  `score = 0.0` junto con la validación individual de cuantiles y mínimos;
  temperatura, mejor score y mejores cuantiles serán también poblacionales.
  El recocido usará `T0 = max(score_inicial, 1e-12)`, aceptará siempre mejoras y
  empates y los empeoramientos con `U < exp(-delta_score/T)`. Enfriará
  `T *= 0.95` tras cada 1000 propuestas, independientemente de cuántas se
  acepten, y se detendrá inmediatamente al alcanzar un estado válido.
  Para cada población se calculará
  `population_haplotype_cells = sum(C * snp_count)` sobre sus bloques y su
  intento quedará limitado a
  `min(10_000_000, max(100_000, 100 * population_haplotype_cells))` movimientos
  propuestos. `max_haplotype_attempts = 1000` se aplicará por población. Agotar
  un intento reintentará solo esa población; agotar su límite hará fallar el
  panel sin resortear las demás y sin convertir el agotamiento en imposibilidad.
  Podrán conservarse resultados válidos internamente, pero la publicación será
  atómica después de que todas las poblaciones superen el validador común. Solo
  entonces se formará la unión densa final de HAP por bloque.
  La derivación RNG 2 separará exactamente `ld_geometry`, `ld_site_alleles`,
  `ld_population_frequencies`, `ld_haplotype_initialization`,
  `ld_annealing_proposals` y `ld_annealing_acceptance`; este último consumirá
  uniformes solo para empeoramientos y no desplazará propuestas posteriores.
  Las derivaciones usarán exclusivamente IDs estables y
  `population_attempt_index` donde corresponda, nunca ruta, reloj, hilo,
  *scheduler* u orden accidental. Los ámbitos poblacionales incorporarán
  `population_id`; cambiar, añadir o reordenar una población no alterará la
  geometría ni frecuencias, haplotipos o recocido de las demás. Los catálogos v1
  y v2 seguirán intactos y legibles. El futuro
  `catalogs/subseed-scopes-v3.tsv` declarará `catalog_version = 3` y
  `derivation_version = 2`, contendrá los ámbitos anteriores necesarios y los
  seis ámbitos LD con vectores normativos propios. La generación automática
  elegirá explícitamente catálogo 3 y derivación 2, sin mezclas silenciosas;
  configuración efectiva y manifiesto registrarán catálogo, derivación,
  contexto, ámbitos, componentes y subsemillas conforme a DTC-037/DTC-038. El
  catálogo real solo se publicará tras cerrar las puertas normativas.
  Los cuantiles empíricos usarán siempre Hyndman–Fan Type 7: sobre valores
  ordenados, `t=(n-1)q`, interpolación entre `floor(t)` y `ceil(t)`, y `x[0]`
  cuando `n=1`, para q 0.10/0.25/0.50/0.75/0.90. La misma regla regirá
  generador, preflight, oráculo e informes.
  El preflight enumerará, para cada pareja no fijada con recuentos `a`, `b` y
  `C`, todos los recuentos conjuntos enteros entre `max(0,a+b-C)` y `min(a,b)`.
  Derivará cotas relajadas de cada cuantil aplicando Type 7 a los mínimos y
  máximos por pareja y solo declarará imposibilidad si la banda solicitada no
  las interseca o faltan parejas evaluables. Ignorará compatibilidades
  multilocus, no admitirá `--force` y se protegerá con enumeración exhaustiva en
  fixtures pequeños para demostrar ausencia de falsos imposibles.
  La representabilidad de `covered_bases` o `coverage_fraction` se comprobará
  aparte de la colocación mediante un grafo exacto de residuos módulo la menor
  longitud admisible que quepa en algún segmento elegible. Caminos mínimos
  obtendrán la menor suma por residuo y predecesores: un objetivo será
  representable si y solo si la menor suma de su residuo no lo supera, y se
  reconstruirá una combinación testigo exacta con aritmética comprobada. Ese
  testigo será evidencia del preflight y no sustituirá el muestreo ponderado.
  Para `block_count`, solo exceder
  `sum(floor(segment_length/min_allowed_length))` demostrará imposibilidad. Pasar
  estas pruebas no demostrará colocabilidad conjunta; el layout agotado y la
  imposibilidad de ejecutar el cálculo por recursos tendrán diagnósticos
  separados y nunca se reclasificarán como suma imposible. Fixtures pequeños
  usarán enumeración o programación dinámica independiente, incluidos intervalos
  degenerados, huecos y recomposición del testigo.
  Los oráculos de aislamiento cambiarán el perfil, reordenarán o añadirán una
  población y exigirán que los artefactos poblacionales no relacionados, o sus
  filas/proyecciones canónicas cuando compartan fichero, permanezcan idénticos.
  También comprobarán que el fallo de una población no publica un panel parcial
  ni resortea las demás.
  Podrá ofrecer, siempre mediante selección explícita, el perfil
  `human_ld_ragsdale_approx_v1`: una aproximación sintética y trazable a la
  tendencia humana `r_div^2` de Ragsdale y Gravel (2020), nunca presentada como
  cuantiles empíricos ni usada por defecto. Sus cinco cuantiles procederán de
  una distribución Beta por población e intervalo, con media dada por la curva
  y `kappa = 8.0` global, no configurable y versionado con el builtin, nunca
  inferido ni ajustado durante la ejecución. Un perfil empírico exacto podrá
  añadirse posteriormente sin bloquear esta primera capacidad. El `profile_id`
  queda ratificado como identificador físico definitivo en todos los contratos.
- **F6.5. Grupos y subpoblaciones.** Asignar identidades, pertenencias y
  parámetros iniciales sin ambigüedades.
- **F6.6. Fundadores.** Identificar individuos y copias haplotípicas fundadoras
  de forma estable.

**Resultado de salida:** población inicial completa, faseada y trazable, apta
para evolución o para producir eventos técnicos.

**Evidencia de cierre:** pruebas estadísticas y controladas validan recuentos,
frecuencias, fase, F, bloques y reproducibilidad. Al menos una población inicial
produce un evento técnico mediante el recorrido de F4 sin una ruta alternativa.

> Revisión de F6: aprobar el estado inicial antes de introducir transmisión y
> generaciones.

### F7. Evolución poblacional

**Objetivo:** simular generaciones discretas manteniendo verdad, linaje y orden
configurado de los procesos.

**Depende de:** F6.

**Hitos:**

- **F7.1. Reproducción y recombinación.** Validar ploidía par e igual y estado
  estructural canónico idéntico, formar bivalentes al azar y transmitir la mitad
  de copias de cada parental con fase y procedencia reconstruibles. Incluir
  *scaffolds* sin recombinación y portadores de SV de novo incompatibles después.
  Reutilizar el mismo motor, mapas y procedencia por segmentos que la
  inicialización fundadora con panel de F6.4, sin convertir aquella creación en
  una meiosis histórica ni introducir un algoritmo alternativo.
- **F7.2. Mutación heredable.** Introducir SNP, *indels*, variantes estructurales
  y repeticiones según sus modelos y momentos confirmados. Contar por clase las
  oportunidades válidas en todas las copias descendientes actuales, registrar
  tasa, espacio, media Poisson y recuentos sorteado, aceptado y activo. Excluir
  `N` de SNP, motivos repetidos y anclas que exigen base concreta sin excluirla
  del interior de intervalos puramente estructurales. **Subpaso diferido F7.2.a
  — mutaciones *de novo*:** aplicar y trazar exclusivamente los cambios
  posteriores a la inicialización, sin reinterpretar la matriz o la variación
  aleatoria inicial directa; cierra esa parte de REQ-IND-64.
- **F7.3. Tamaño poblacional.** Procesar N intentos reproductivos, registrar
  éxitos, fallos y N realizado, y calcular Ne reproductivo realizado cuando sea
  evaluable. Ne objetivo y el control de contribuciones parentales quedan fuera
  de v1.
- **F7.4. Selección y deriva.** Aplicar selección direccional de viabilidad sobre
  un único SNP bialélico protegido, validar todas las dosis y aptitudes antes de
  supervivencia, rechazar un locus sobre `N` y demostrar los fallos por identidad
  o máximo cero. La deriva
  surge del muestreo aleatorio en el orden configurado.
- **F7.5. Migración y subpoblaciones.** Ejecutar migraciones simultáneas,
  divisiones, fusiones, fundaciones y extinciones.
- **F7.6. Fases evolutivas.** Variar parámetros y encadenar episodios como
  panmixia, cuellos de botella o migración sin alterar el orden generacional
  fijo. Registrar cada etapa aplicada u omitida y el orden interno de varios
  eventos demográficos. Completar en `Evolución` la autoría de uno o varios TOML
  de escenario y su equivalencia con el conjunto manual.
- **F7.7. Observaciones.** Muestrear generaciones sin alterar la trayectoria y
  materializar estado completo y linaje detallado para población inicial,
  final y generaciones o muestreos solicitados; conservar solo resúmenes de las
  generaciones restantes, manteniendo internamente lo necesario para reconstruir
  los individuos materializados.

**Resultado de salida:** poblaciones evolucionadas reproducibles con linaje
hasta los fundadores.

**Evidencia de cierre:** casos pequeños permiten reconstruir cada transmisión y
casos estadísticos comprueban que los procesos realizados concuerdan con sus
parámetros mediante informes de esperado, realizado, tamaño observado y
desviación, sin usar una realización aleatoria para aceptar o rechazar la
ejecución. Una muestra evolutiva puede producir los resultados técnicos de F4
sin modificar su verdad genotípica.
También debe demostrar las ocho etapas completas, la omisión de una etapa
opcional en su posición, el rechazo previo de un orden alternativo y la
trazabilidad de varios eventos demográficos ordenados dentro de su etapa.
Debe incluir una generación materializada y otra no materializada y demostrar
que esta última solo se enumera mediante su resumen sin impedir reconstruir el
linaje solicitado.

> Revisión de F7: comprobar exactitud y consumo de recursos antes de añadir
> grafos de cruces complejos.

### F8. Familias y poblaciones segregantes

**Objetivo:** crear descendencia mediante grafos declarativos de cruces y
conservar todos los individuos y haplotipos intermedios.

**Depende de:** F7 para reproducción, recombinación, mutación y linaje.

**Hitos:**

- **F8.1. Grafo de cruces.** Validar parentales, ciclos, dependencias y reglas
  declarativas antes de generar descendencia, y completar en `Familias` la
  autoría de sus escenarios enlazados sin alterar IDs, orden ni rutas relativas.
- **F8.2. Cruces básicos.** Implementar autofecundación, cruzamiento sexual,
  F2 y retrocruce como composiciones verificables que rechazan parentales con
  ploidía impar o desigual o con diferencias estructurales.
- **F8.3. Selección de descendientes.** Elegir por identificador, azar,
  pertenencia, genotipo, dosis o haplotipo mediante el catálogo común.
- **F8.4. Diseños complejos.** Expresar pedigrís, líneas de introgresión y MAGIC
  mediante el mismo modelo, sin código ejecutable en la configuración.
- **F8.5. Poblaciones intermedias.** Materializar el estado faseado completo y
  los metadatos de todos sus individuos en todos los eventos; generar FASTQ,
  BAM/CRAM o VCF únicamente cuando se solicite.

**Resultado de salida:** familias y poblaciones segregantes reconstruibles desde
sus fundadores y cruces.

**Evidencia de cierre:** casos diploides y tetraploides con estructura idéntica,
variantes pequeñas permitidas, cada clase de SV heterocigota, una SV fijada,
ploidías inválidas, un *scaffold* sin recombinación y una SV de novo seguida de
otro cruce demuestran parentesco, compatibilidad, segregación, fase,
recombinación y selección de descendientes.
La evidencia debe comprobar que una familia conserva y enumera todos sus
individuos, incluidos los intermedios que no alimentan una salida técnica.

> Revisión de F8: validar los grafos, los cruces y su integración con los
> eventos técnicos antes de la consolidación.

### F9. Escalabilidad y consolidación

**Objetivo:** convertir las capacidades anteriores en una primera versión
robusta, distribuible y capaz de manejar poblaciones grandes dentro de los
recursos disponibles.

**Depende de:** todas las fases anteriores incluidas en la versión que vaya a
publicarse.

**Hitos:**

- **F9.1. Escenarios de rendimiento.** Fijar casos representativos y medir CPU,
  memoria, almacenamiento, individuos, generaciones y lecturas por segundo con
  el contrato único de trazabilidad y todas las salidas activadas.
- **F9.2. Representación compacta.** Optimizar, con evidencia, haplotipos,
  segmentos, genotipos, historia y materialización de secuencias.
- **F9.3. Paralelismo reproducible.** Escalar por unidades independientes sin
  cambiar resultados con la misma semilla.
- **F9.4. Presupuesto de recursos.** Estimar memoria, almacenamiento y tiempo,
  emitir avisos y rechazar antes de simular ejecuciones que no puedan cumplir el
  mínimo, la trazabilidad y las salidas activadas. Nunca reducir detalle de
  forma automática; el escenario deberá pedir menos materializaciones,
  muestreos, individuos, lecturas u otras salidas.
- **F9.5. Interoperabilidad.** Validar versiones, cabeceras, compresión, índices
  y compatibilidad de FASTA, FASTQ, BAM, CRAM y VCF.
- **F9.6. Exportación autocontenida.** Empaquetar eventos con referencias,
  índices, configuraciones, manifiestos y sumas de comprobación.
- **F9.7. Preparación de publicación.** Completar documentación, ejemplos,
  pruebas de regresión y proceso reproducible de construcción; construir y
  validar oficialmente `.tar.gz` y OCI, y ejecutar la evidencia de plataforma
  y distribución detectada de DTC-074.

**Resultado de salida:** primera versión candidata a publicación.

**Evidencia de cierre:** la batería funcional, genómica, estadística, de
interoperabilidad y rendimiento supera los criterios acordados sin hallazgos
bloqueantes abiertos.

> Revisión de F9: decidir si la versión está lista para publicarse o necesita un
> ciclo adicional de corrección.

DTF-037 reserva para después de F9.7 una decisión funcional separada sobre
lotes reproducibles de ejecuciones completas. No es un hito de v1: solo podrá
planificarse cuando existan las vías F3–F8 y la consolidación F9, y entonces
deberá incluir la misma paridad Marimo exigida a cualquier capacidad
configurable.

## 6. Matriz resumida de comprobación

| Fase | Resultado construido | Evidencia principal | Riesgo protegido |
|---|---|---|---|
| F0 | Arquitectura y plan de F1 | Trazabilidad completa y decisiones cerradas | Diseño improvisado durante la implementación |
| F1 | Ejecución mínima reproducible y referencia reducida | Repetición con semilla, reconstrucción origen+mapa y rechazo limpio | Resultados sin identidad, configuración efectiva o procedencia de referencia |
| F2 | Aplicación base, matriz SNV e individuos faseados | Equivalencia Marimo→TOML, matriz suministrada/generada, HWE/F y reconstrucción independiente de haplotipos/VCF | Perfil implícito, loci distintos por individuo, genotipos/fase incorrectos o semántica duplicada en UI |
| F3 | Primer recorrido extremo a extremo | FASTQ, compresión, hashes y comprobación contra la verdad biológica | Lecturas o metadatos técnicos incoherentes |
| F4 | Realismo técnico WGS | Comparación de parámetros solicitados y realizados | Alterar la verdad con artefactos técnicos |
| F5 | Variación y referencias complejas | Casos estructurales, VCF más complemento y aislamiento por referencia | Mezcla de referencias o pérdida de linaje o detalle complejo |
| F6 | Población inicial y panel haplotípico intacto | Recuentos, frecuencias, reconstrucción del panel y pruebas estadísticas reproducibles | Frecuencias, fase, LD o F mal representados |
| F7 | Evolución y recombinación meiótica | Reconstrucción de transmisiones, mosaicos y medidas realizadas | Orden evolutivo, puntos de cruce, procedencia o genealogía incorrectos |
| F8 | Familias y segregantes | Casos mendelianos y grafos validados | Cruces o selección de descendientes erróneos |
| F9 | Versión consolidada | Pruebas completas y escenarios de rendimiento | Fallos de escala, distribución o interoperabilidad |

## 7. Reglas de avance

Una fase puede comenzar a planificarse cuando sus dependencias estén cerradas.
Puede comenzar a implementarse cuando disponga de arquitectura suficiente, un
plan aprobado y criterios de salida verificables.

Una fase no se cerrará si:

- existe un defecto bloqueante o grave sin resolución o destino aprobado;
- falta una prueba independiente para una invariante crítica;
- los archivos producidos no superan sus validaciones de formato;
- la documentación afirma un resultado que no se ha ejecutado;
- una decisión pendiente cambiaría los contratos de la fase siguiente;
- una opción configurable implementada desde F2 carece de control Marimo,
  equivalencia UI→TOML→`EffectiveConfig` o comprobación accesible.

Todo cierre termina con una pausa para revisar evidencias y decidir si se
continúa, se corrige o se modifica la hoja de ruta.

## 8. Trabajo posterior a esta hoja de ruta

F0 queda cerrada con la arquitectura Rust 0.1 y las decisiones base. F1.6 y
F2.0 son ampliaciones posteriores: no se consideran implementadas ni aprobadas
por aparecer en esta hoja y deben respetar los planes, puertas y pausas.

Las fases F3–F9 no deben dividirse todavía en archivos, módulos, crates ni
commits. Su desglose se hará progresivamente con la evidencia obtenida en las
fases anteriores; F2 ya dispone de su plan 0.1, todavía pendiente de aprobación.
