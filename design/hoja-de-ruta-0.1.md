# Hoja de ruta del simulador de genomas y poblaciones

- **Versión:** 0.2
- **Estado:** vigente; F0 aprobada, F1 cerrada en `e14705a`/`efca265`; F2
  simplificada por DTC-079 y F5 reconstruida documentalmente
- **Fecha:** 1 de septiembre de 2026
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
                 └─ F5 Población inicial y haplotipos sintéticos
                     └─ F6 Variación y referencias complejas
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

> **Reapertura acotada, 23/09/2026.** Esta fase estaba cerrada y se le añadió
> capacidad: el generador de referencia simulada admite **regiones declaradas**
> además del relleno al azar, y una referencia real admite la **inserción de
> tractos** con su tabla de correspondencia de coordenadas. Lo que **no** se
> reabrió: ninguna subsemilla, ningún vector congelado y ningún esquema
> publicado. Con la lista de regiones vacía, los bytes son los de antes. El
> registro completo —qué no se reabre, qué cambia y por decisión de quién,
> dónde vive cada literal, las pruebas con su variante defectuosa y la
> evidencia repetida— está en la **sección 14** de
> `design/plan-implementacion-f1-0.1.md`.

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
  parcial y cierran en F6.2/F7.2 según su alcance. DTPD-118 está resuelta por
  DTC-054.
- **F2.3. Individuos directos desde población de origen.** Construir individuos
  asignados a una única población y muestrear la matriz mediante Hardy-Weinberg
  para P=1..10 o `F` solo diploide, con fase y subsemillas independientes. La
  distribución genotípica explícita no tiene fichero, esquema ni control Marimo
  en F2; REQ-CONF-130 aporta evidencia parcial y cierra en F6.2. Integrar restricciones
  explícitas, incluidas las dos autorizaciones limitadas y trazadas para
  genotipos incompatibles o matriz activa vacía; integrar los demás cambios y
  reconstruir cada haplotipo. REQ-IND-02/03 aportan aquí evidencia parcial sobre
  una referencia estructural por individuo y cierran en F6.2. Activar preparación/generación/reutilización,
  perfiles, población y modelo en la pestaña `Individuos`.
- **F2.4. Verdad estructurada.** Publicar `truth.jsonl` obligatorio (v1 en F2;
  v4, o v5 con alelos de repetición, desde el 30/09/2026), UTF-8
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
F2.5.4/F2.5.5/F2.5.6; sus cierres completos siguen en F8.5 y F6.5,
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

> **Reapertura acotada, 24/09/2026.** Esta fase estaba cerrada con su evidencia
> aprobada y se le añadió un paso: el **tartamudeo de base** de la polimerasa
> sobre los microsatélites declarados, intercalado entre la quimera y calidad y
> errores. **No es una preparación amplificada**, no introduce duplicados y el
> perfil se sigue registrando como sin PCR: el deslizamiento ocurre en la
> amplificación en puente sobre la celda y en la síntesis, que existen en el
> perfil vigente. **La verdad del individuo no cambia.** Lo que **no** se
> reabrió: ninguna clave aleatoria heredada —la etapa nueva estrena catálogo y
> solo ella lo usa—, ningún esquema publicado y ningún artefacto de un
> escenario sin loci declarados, que sigue dando los mismos bytes. El registro
> completo está en la **sección 16** de
> `design/plan-implementacion-f4-0.1.md`.

> **Segunda reapertura acotada, 25/09/2026.** La misma fase recibió un segundo
> paso: la **preparación de biblioteca amplificada por PCR**, por elección
> expresa y con la preparación sin PCR como predeterminada. Con ella, el
> recorrido deja de producir una pareja por fragmento y pasa a producir
> **familias**: de ahí salen los duplicados, la sobrerrepresentación de unas
> moléculas frente a otras y el sesgo por composición. El usuario declara **la
> fracción de duplicados** que ve en su informe de calidad, no el parámetro del
> modelo, que se deriva por inversión numérica. **La verdad del individuo sigue
> sin cambiar.** Lo que **no** se reabrió: el recurso del perfil sin PCR y su
> huella —la amplificación estrena un perfil **hermano**—, ninguna clave
> aleatoria heredada —la etapa `amplification` estrena el catálogo
> `f4-technical-v5` y solo ella lo usa—, ningún esquema publicado, y ningún
> artefacto de un escenario que no la declare, anclado desde ahora por una
> **huella congelada** y no solo por comparar dos ejecuciones entre sí. Lo que
> **quedaba pendiente** era el **tartamudeo heredable por linaje**, de la
> entrega 4, y mientras no existió la preparación amplificada **rechazaba** los
> escenarios que declaraban loci de repetición. **Ese rechazo se levantó el 26
> de septiembre de 2026**; ver la nota de esa fecha, más abajo. El modelo queda
> fijado por **DTC-091**.

> **Nota del 27 de septiembre de 2026 — generador de alelos por deslizamiento
> biológico, y cierre de la simulación de microsatélites.** Hasta ahora los
> alelos de un microsatélite y sus frecuencias **los escribía el usuario a
> mano**. El simulador puede generarlos: parte del alelo que el genoma de
> referencia lleva en el locus, simula **un linaje independiente por copia
> haplotípica** de la población, da a cada uno un número de pasos sorteado de una
> Poisson de media declarada, y **la frecuencia de cada alelo es la proporción de
> linajes que terminan en él**. Va en la **verdad**, no en las lecturas.
>
> **Es una sola tirada para inicializar la diversidad**, no una simulación de
> generaciones: el modelo escalonado aplicado generación tras generación
> pertenece a la fase de evolución y no se ha invadido.
>
> **El alelo de la referencia no tiene la presencia garantizada.** Entra en el
> abanico si algún linaje termina en él; que desaparezca depende de los
> parámetros que se den, no del modelo, y con una media alta perderse es lo que
> ocurre en una población real.
>
> **Tres limitaciones viajan declaradas con cada ejecución**, porque son
> decisiones biológicas y no detalles: no hay ancestría compartida, de modo que
> el espectro sale más plano que el de una población real y eso importa si el
> conjunto se usa para probar herramientas de estructura poblacional; el número
> de pasos no es un proceso en el tiempo; y el genoma de referencia puede llevar
> un alelo que ningún individuo tiene.
>
> **El deslizamiento biológico y el tartamudeo no comparten azar**, aunque
> compartan mecanismo físico: uno decide la verdad y el otro las lecturas. Si los
> alelos verdaderos salieran del modelo de tartamudeo, cambiar el número de
> ciclos de PCR cambiaría el genotipo real de los individuos.
>
> Las decisiones son **DTC-094** a **DTC-097**. Con esto se cierra la regla de
> cierre de la interfaz para este trabajo: toda opción configurable que ha
> añadido tiene control en Marimo, y una comprobación recorre la lista entera.
>
> **Lo que el cambio mueve en documentos ya publicados** está medido y
> registrado en la sección 14 del plan de microsatélites y amplificación: un
> escenario **sin** microsatélites no mueve un solo byte, y uno con
> microsatélites **declarados a mano** gana únicamente el marcador de origen de
> cada alelo.

> **Nota del 26 de septiembre de 2026 — tartamudeo clonal heredable por
> linaje.** Con preparación amplificada, el deslizamiento puede ocurrir en
> **cualquier ciclo** de la PCR, y una copia que resbala transmite su número de
> unidades a **toda su descendencia**. Un resbalón en el primer ciclo alcanza a
> media familia; uno en el último, a una sola copia. De ahí sale el fenómeno que
> importa: **una familia de duplicados puede mostrar dos o tres alelos en
> proporciones que imitan a un heterocigoto**, cuando procede de una sola
> molécula. En genética forense ese artefacto se confunde con un genotipo real.
>
> La genealogía es **mínima y declarada como aproximación**, no como el
> coalescente exacto del proceso de ramificación: progenitor uniforme entre las
> copias anteriores, ciclo de divergencia geométrico pinzado en el último, y
> herencia de los resbalones **anteriores** a la divergencia. Las tasas por
> ciclo se **despejan** de las dos columnas de la calibración —restando, nunca
> sumando, porque cada columna es un total— y se publican **marcadas como
> toscas**: las dos columnas proceden de dos individuos distintos, de proyectos
> distintos y a coberturas distintas.
>
> **El rechazo temporal de la entrega 3 queda levantado.** En su lugar hay una
> condición permanente: con loci declarados, la preparación amplificada **exige
> el tartamudeo de base activo**, porque el clonal solo despeja el incremento de
> la amplificación y sin la parte de base la biblioteca subestimaría lo medido
> entre un 14 % y un 28 %.
>
> Lo que **no** se reabrió: ninguna clave aleatoria heredada —la etapa de
> linaje estrena el catálogo `f4-technical-v6` y **solo ella** lo usa, de modo
> que una biblioteca amplificada sin microsatélites sigue publicando el v5 byte
> a byte—, ningún esquema de una fase cerrada, y ningún artefacto de un
> escenario que no declare loci. **La verdad del individuo sigue sin cambiar**:
> el tartamudeo vive solo en las moléculas que se secuencian. El modelo queda
> fijado por **DTC-092**.

**Hitos:**

- **F4.0. Catálogo RNG técnico.** Fijar y probar las claves, consumidores y
  vectores reproducibles del recorrido F4 antes de producir pares o artefactos.
- **F4.1. Perfiles tecnológicos.** Completar Illumina *paired-end* para una
  población y un evento técnico únicos. PacBio HiFi se aplaza a DTF-049 y los
  recorridos con varias poblaciones o eventos a DTF-050.
- **F4.2. Cobertura, longitud, error y calidad.** Aplicar perfiles validados y
  sustituciones configurables por los ámbitos permitidos.
- **F4.3. Artefactos WGS.** Simular quimeras y rechazar configuraciones de
  duplicados, adaptadores, sesgo GC no neutro y anomalías *paired-end*, que se
  aplazan en DTF-054–056. Las ventanas y correcciones GC no neutras pertenecen
  a DTF-056.
- **F4.4. Contaminación biológica.** Sustituir lecturas completas mediante un
  porcentaje sorteado por receptora y donantes aleatorios de su misma población;
  matrices dirigidas y fuentes externas quedan en DTF-052.
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

### F5. Población inicial y haplotipos sintéticos

**Objetivo:** generar una población fundadora faseada con bloques de LD sintéticos,
sin panel externo, SV ni meiosis. **Depende de:** F4 cerrada. F6 añade SV; F7.1
inicia meiosis.

- **F5.0. Contratos y preflight.** Cerrar los esquemas completos general v3,
  escenario F5 v1, efectivo v5, manifiesto v5, perfiles/TSV/fallo y verdad v2
  sin modificar v2/v4 ni verdad v1; admitir exactamente un escenario y
  resolver referencia, subpoblaciones, perfiles, IDs, RNG, modo/tolerancia por
  bloque×subpoblación, `f5.seed`, auditorías DTC-049, `S_bound`, bins, filas,
  presupuestos dimensionalmente separados y diagnóstico.
- **F5.1. Geometría y marginales.** Generar layout compartido, elegibilidad,
  cobertura exacta y marginales por bloque×subpoblación.
- **F5.2. Haplotipos y LD.** Aplicar DTC-089, r2(d), signos y mejor intento.
- **F5.3. Fundadores intactos.** Muestrear probabilístico/controlado, F
  diploide y verdad fundadora sin recombinación.
- **F5.4. Procedencia y oráculo.** Publicar y validar los siete TSV y la
  proyección fundadora `truth.jsonl` v2, incluidos MAN-70 recontado solo desde
  verdad por locus×subpoblación, vectores, lector
  independiente, cotas de recursos y manifiestos atómicos.
- **F5.5. Auditoría de identidad y publicación.** Comprobar que identidad y
  pertenencia, materializadas antes de los artefactos en F5.3, no se modifican
  después del paquete atómico.
- **F5.6. Marimo y cierre.** Comprobar paridad y gates de cierre.

No hay panel suministrado, construcción recombinante fundadora ni dependencia
F5→F6; la trazabilidad antigua se migra por REQ-EVO-74–77 y ACE-117.

### F6. Referencias y variación complejas

**Objetivo:** ampliar la verdad F5 sin alterar sus fundadores intactos. **Depende
de:** F5 cerrada. **F6.1. Referencias estructurales.** **F6.2. Variantes
estructurales.** **F6.3. Ploidía generalizada.** **F6.4. Repeticiones.**
**F6.5. Variantes complejas y complemento VCF.** **F6.6. Conjuntos y matrices
por referencia.** F6 no hereda
panel ni recombinación de F5. Su cierre demuestra reconstrucción, límites,
solapamientos y aislamiento entre referencias.

**Hito nuevo, sin desarrollar (24/09/2026):**

- **F6.4.a. Localizar los microsatélites de un genoma suministrado.** Barrer una
  referencia real y **descubrir** los tractos que ya contiene, en vez de exigir
  que el usuario los declare uno a uno. Hoy el simulador **no lo hace**: solo
  sabe comprobar que un locus declarado existe de verdad donde se dice (paso
  F6.4.2, hecho) y, cuando se desarrolle F6.4.6, generar loci nuevos al azar por
  regiones. Descubrir los que ya están es la capacidad hermana de esa última y
  falta por completo.

  **Por qué importa más allá de la comodidad.** Sin ella, todo microsatélite es
  declarado por el usuario, y por eso el tartamudeo tiene un interruptor
  explícito: quien declara el locus decide también si tartamudea y aporta sus
  parámetros. Con descubrimiento automático aparecería una clase nueva de locus
  —el que el simulador encontró y el usuario no pidió— y habría que decidir su
  régimen por omisión. El responsable lo anticipó el 24 de septiembre de 2026.

  **Lo que hay que decidir cuando se desarrolle**, y que este hito no resuelve:
  el criterio de descubrimiento —qué pureza y qué tirada mínima hacen que un
  tramo cuente como microsatélite—, qué pasa si el usuario declara uno que el
  barrido no encontró o al revés, el coste de barrer un genoma completo, y el
  régimen por omisión del tartamudeo para los descubiertos.

**Subpaso diferido:**

- **F6.5.a. Variación aleatoria directa.** Configurar cantidad media o densidad,
  unidad, distribuciones y regiones globales y por individuo, y comprobar que la
  creación inicial directa consume las sustituciones individuales sin alterar la
  matriz poblacional.

> Revisión de F6: confirmar que la complejidad estructural no rompe la verdad
> fundadora de F5 ni el realismo técnico establecido en F4.

### F7. Evolución poblacional

**Objetivo:** simular generaciones discretas manteniendo verdad, linaje y orden
configurado de los procesos.

**Depende de:** F6.

**Hitos:**

- **F7.1. Reproducción y recombinación.** Validar ploidía par e igual y estado
  estructural canónico idéntico, formar bivalentes al azar y transmitir la mitad
  de copias de cada parental con fase y procedencia reconstruibles. Incluir
  *scaffolds* sin recombinación y portadores de SV de novo incompatibles después.
  Definir el motor, mapas y procedencia por segmentos en F7.1; no reutiliza un
  panel ni una inicialización recombinante F5, que no existen.
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
| F5 | Población inicial y haplotipos sintéticos | Siete TSV, reconstrucción H×S/q/D/`r²`, fundadores y pruebas reproducibles | Frecuencias, fase, LD o F mal representados |
| F6 | Variación y referencias complejas posteriores | Casos estructurales, VCF más complemento y aislamiento por referencia | Mezcla de referencias o pérdida de linaje o detalle complejo |
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
commits salvo en su plan de fase aprobado. Su desglose se hace progresivamente
con la evidencia obtenida en las fases anteriores. F5 dispone del plan 0.1
jerárquico F5.0–F5.6 y no crea identificadores `REQ-F5-*` paralelos.
