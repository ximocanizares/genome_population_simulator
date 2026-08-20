# Plan de implementación de F2 0.1

## 1. Cabecera y estado

- **Título:** plan detallado de F2, verdad genómica mínima.
- **Versión del plan:** 0.1.
- **Fecha:** 20 de agosto de 2026.
- **Estado:** corregido después del cierre simplificado de F1; pendiente de una
  nueva revisión adversarial y de aprobación explícita. La revisión F2 0.4 es
  histórica y no aprueba esta versión, porque precede a DTC-078 y a esta
  corrección. F2 no se ha iniciado.
- **Ámbito temporal:** únicamente F2.0–F2.5.
- **Fuentes normativas y de diseño:**
  [especificaciones-0.1.md](./especificaciones-0.1.md),
  [arquitectura-0.1.md](./arquitectura-0.1.md),
  [decisiones-tecnicas-0.1.md](./decisiones-tecnicas-0.1.md),
  [hoja-de-ruta-0.1.md](./hoja-de-ruta-0.1.md),
  [trazabilidad-requisitos-0.1.md](./trazabilidad-requisitos-0.1.md),
  [registro-ampliaciones-futuras-0.1.md](./registro-ampliaciones-futuras-0.1.md)
  y [plan-reimplementacion-f1-0.2.md](./plan-reimplementacion-f1-0.2.md).
- **Decisiones confirmadas aplicables:** DTC-001, DTC-003, DTC-004, DTC-006,
  DTC-010, DTC-011, DTC-013, DTC-020 a DTC-022, DTC-024, DTC-031 a
  DTC-049 y DTC-051 a DTC-075, DTC-078.
- **Decisiones resueltas y puertas internas:** DTPD-102, DTPD-103, DTPD-120 y
  DTPD-121 están resueltas por DTC-061 y DTC-062 a DTC-068. REQ-CONF-130 queda
  limitado en F2 a HWE y F diploide; su distribución explícita cierra en F6.2.
  DTC-069 resuelve la unidad visible de importación/exportación Marimo.
  DTC-058 resuelve DTPD-117 sin alternativas técnicas abiertas, pero F2.0 queda
  bloqueada hasta que F2.0.2 verifique y apruebe las versiones candidatas, el
  entorno local/offline y el contrato versionado con el binario. La construcción
  y validación oficial de `.tar.gz`, OCI y plataforma/distribución se trasladan
  a F9. DTC-054
  resuelve DTPD-118 y fija para
  F2.2/F2.3 los formatos de matriz/espectro, el algoritmo reproducible, los
  perfiles versionados y las excepciones limitadas. DTC-055 resuelve DTPD-119 y
  deja preparada la inicialización fundadora desde paneles para F6.4; no se
  consume dentro de F2. DTC-055 ya cierra para F6.4 el constructor común y su
  variante uniforme; F7.1 lo reutilizará después para la recombinación meiótica
  y DTPD-111 solo resolverá mapas y semántica meiótica.
- **Base cerrada y corrección pendiente:** F1 está cerrada en `e14705a` y su
  cierre documental en `efca265`. RAF2-027 queda satisfecho por esa base;
  RAF2-028 a RAF2-044 conservan sus decisiones. Esta corrección sustituye las
  premisas F1 retiradas; no cierra una nueva revisión adversarial de F2.

Todas las casillas permanecen abiertas. El documento ordena trabajo futuro: no
  inicia F2, no resuelve decisiones pendientes, no selecciona
crates, formatos ni herramientas externas no confirmados y no autoriza commits.
Las versiones candidatas de DTC-058 siguen sujetas a F2.0.2. F2 continúa sin
aprobarse ni autorizada para iniciar.

## 2. Resultado y alcance

F2 termina con un recorrido directo pequeño pero completo que, partiendo de una
referencia de F1, construye individuos de ploidía uniforme, conserva copias
haplotípicas separadas, aplica SNP e inserciones/deleciones pequeñas inequívocas
y publica `truth.jsonl` canónico junto con exactamente un VCF de verdad
multimuestra por conjunto. La fase, la dosis, la procedencia, las coordenadas y
los alelos deben reconstruirse sin consultar lecturas, alineamientos ni el propio
VCF como fuente de verdad.

### 2.1. Incluido

- Una referencia pequeña y una ploidía uniforme configurable entre 1 y 10 para
  todas sus secuencias en el recorrido de F2; la variación de ploidía por
  secuencia sigue reservada a F5.3.
- Estado haplotípico híbrido DTC-061 por individuo, secuencia y copia, con fase,
  dosis y procedencia: cadena ordenada/persistente de segmentos que referencia
  fuentes, modificaciones pequeñas tipadas y materialización exacta solo bajo
  demanda, sin duplicar secuencias completas.
- SNP e inserciones/deleciones pequeñas inequívocas, explícitas y aleatorias,
  sobre toda una secuencia o regiones declaradas. Los SNV heredables iniciales
  proceden de una matriz común a la población, suministrada o generada desde
  número/densidad, espacio, espectro de sustituciones y espectro alélico; no se
  generan por densidades independientes por individuo. Solo un recuento pedido
  como media o densidad usa Poisson: ni frecuencias, alelos ni genotipos se
  muestrean mediante Poisson. La longitud de los *indels* aleatorios de F2 es un
  valor exacto configurable por clase, predeterminado a 1 pb y limitado al
  intervalo desde 1 hasta `umbral estructural efectivo - 1`.
- Coordenadas externas basadas en 1 y cerradas, inserción anclada en la base de
  referencia precedente y solo interna, y deleción interna que deja al menos
  una base a cada lado.
- Conservación literal de `N`: un SNP nunca tiene origen ni destino `N`, una
  inserción aleatoria usa solo `A/C/G/T`, y una inserción explícita o copiada
  puede conservar `N`.
- Configuración declarada, validada y efectiva del recorrido directo, ampliando
  la frontera mínima F1 solo con los tipos, relaciones e IDs que tengan su
  primer consumidor F2; reutiliza sintaxis TOML, derivación RNG, integridad,
  manifiesto general y publicación F1, sin suponer reglas ni resolvedores F1.
- Verdad estructurada obligatoria e independiente del VCF, y un VCF de verdad
  mínimo con GT faseado, dosis correcta, orden de referencia, alelos derivados
  de la verdad y validación externa independiente.
- Aplicación Marimo base con nueve pestañas, cobertura completa de F1/F1.6,
  importación/exportación del conjunto TOML e invocación segura del binario; durante F2 se
  amplía para individuos, `truth.jsonl` y VCF.

### 2.2. Excluido

- Lecturas, FASTQ, BAM/CRAM, mapeo y cualquier error técnico; comienzan en F3.
- Variantes estructurales, repeticiones, referencias derivadas, ploidía distinta
  por secuencia y la cobertura VCF/complemento completa de casos complejos; F5.
- Dinámica de poblaciones, evolución, recombinación, selección, familias y
  cruces; F6–F8. F2 solo usa la población de origen necesaria para interpretar
  la matriz de frecuencias y crear individuos directos, sin simularla.
- Paneles de haplotipos, desequilibrio de ligamiento y creación fundadora
  intacta o recombinante desde ellos; DTC-055 destina ambas alternativas a
  F6.4. F7.1 reutiliza el motor común en la reproducción, pero no inicializa el
  conjunto fundador.
- Distribuciones no Poisson para recuentos, distribuciones probabilísticas de
  longitud de *indels*, gVCF con posiciones invariantes y cualquier resolución
  silenciosa de `N`. Las distribuciones configurables de longitud corresponden
  a REQ-MUT-08 en F7.2; F2 solo admite la longitud exacta acordada.
- DTF-003, DTF-008, DTF-013, DTF-014, DTF-023, DTF-033, DTF-034, DTF-036 y
  DTF-037,
  además de las demás ampliaciones post-v1 del registro vigente.
- Afirmar el cierre total de requisitos asignados a fases posteriores. F2 puede
  dejar evidencia parcial reutilizable para, por ejemplo, REQ-IND-58,
  REQ-IND-59, REQ-VER-01, REQ-VER-24 a REQ-VER-26, REQ-FMT-07 y REQ-SAL-02,
  pero no los marcará satisfechos en su matriz de cierre.
- Reimplementar en Marimo validación, dominancia, biología, generación aleatoria
  o publicación; la aplicación solo configura e invoca al binario.

## 3. Orden, reutilización y fronteras

1. Partir exclusivamente de F1 cerrado en `e14705a` y de su cierre documental
   `efca265`; la autoridad de planificación F1 es `plan-reimplementacion-f1-0.2`.
2. Validar DTC-058 en F2.0.2 y cerrar F2.0 antes de F2.1. Consumir sin reabrir
   DTC-061 para implementar el estado haplotípico y `truth.jsonl`; consumir
   DTC-062 a DTC-066 como resolución completa de DTPD-103 antes de implementar
   el VCF.
3. Reutilizar y extender el único servicio de coordenadas decidido e iniciado
   en F1.6. Las fronteras traducen la base 1 cerrada configurada; ningún
   adaptador ni algoritmo realiza conversiones implícitas ni abre otra
   representación interna.
4. Extender `DeclaredConfig` y `EffectiveConfig` desde su superficie mínima F1.
   F2 crea el primer catálogo/formas de reglas y el normalizador solo cuando
   haya consumidores reales; no presupone `ParameterKey`, `ResolutionContext`,
   `ResolutionDomain` ni `resolution-evidence` existentes. Reutiliza solo la
   sintaxis TOML, derivación RNG, detección de colisiones de IDs, hashes y
   publicación realmente cerrados en F1.
5. Crear primero verdad de dominio y materialización; derivar después el fichero
   estructurado y, por último, el VCF. Ninguna salida técnica se convierte en
   estado canónico.
6. Aplicar eventos explícitos en orden declarado y después clases aleatorias en
   el orden normativo aprobado. Los ámbitos aleatorios se registran antes de
   que exista un consumidor y no dependen del orden físico de ejecución.
7. Introducir cada dependencia mediante DTC-043 y cada cambio de algoritmo,
   catálogo o versión como cambio aislado que repite vectores, reproducibilidad
   y validadores.
8. Terminar cada hito con una pausa obligatoria. Un hito consumidor no comienza
   mientras la revisión o las evidencias del productor estén pendientes.

### 3.1. Línea base cerrada y frontera de extensión

F1 está cerrada en `e14705a` y documentada en `efca265`. Su base reutilizable
es deliberadamente mínima: TOML general para semilla, referencia y reducciones;
RNG v1–v3 para `reference_generation` y `reference_reduction`; IDs, SHA-256,
manifiesto general y publicación local. No existen aún vía, `scenario`,
ploidía, reglas, `ParameterKey`, `ResolutionContext`, dominio de resolución,
evidencia de resolución, eventos ni catálogos de locus/posición. F2 no supone
dependencias para Poisson, JSONL o VCF antes de evaluarlas.

## 4. Prerrequisitos, decisiones y comprobaciones comunes

### 4.1. Puerta de entrada a F2

Antes de F2.0.1 deben constar: F1.4, F1.5 y F1.6 terminadas; revisión de cierre de F1
aprobada; pruebas y evidencias F1 verdes; hashes/esquemas/catálogos F1 congelados;
árbol de trabajo atribuido; y aprobación explícita de este plan. Si F1 cambia
un contrato que F2 pretendía reutilizar, se actualiza primero este plan y su
trazabilidad.

### 4.2. Decisiones que detienen solo a sus consumidores

- DTPD-101 está resuelta mediante DTC-056. F1.6 implementa representación común,
  tipos de frontera y vectores de intervalos origen↔derivada; F2.1/F2.2 deben
  auditarlos y extender el mismo servicio para segmentos, anclas y composición
  sin reabrir la decisión.
- DTC-061 resuelve DTPD-102: modelo interno híbrido de segmentos/modificaciones
  y `truth.jsonl` v1 obligatorio, UTF-8 sin comprimir, ordenado y validado por
  esquema por registro más reconstructor semántico independiente. F2 consume
  ese contrato sin reabrir formato, compresión ni dirección de derivación.
- DTC-062 fija VCF 4.5 como única versión y la cabecera exacta
  `##fileformat=VCFv4.5`, con validación externa basada en htslib. DTC-062 a
  DTC-066 ya completan normalización, representación, campos `INFO`, complemento
  y versiones/comandos de validadores; F2.5 los consume y verifica sin reabrir
  DTPD-103 ni producir VCF 4.2/4.3/4.4.
- DTC-063 fija `.vcf.gz` BGZF con exactamente un índice TBI o CSI como salida
  predeterminada, seleccionado por el límite normativo de TBI, y `.vcf` de
  texto sin índice como alternativa explícita. La configuración efectiva y el
  manifiesto registran modo, tipo y rutas; VCF e índice se validan y publican
  atómicamente.
- DTC-064 fija la normalización interna determinista/versionada desde los
  haplotipos finales, agrupación de todos los alelos compatibles en un único
  registro multialélico por `CHROM+POS` y `bcftools norm` solo como oráculo.
  Los cambios solapados se reconstruyen antes sobre un intervalo `REF` común;
  la excepción de varios BND se reserva a SV/translocaciones futuras.
- DTC-065 fija los `INFO` declarados `GPS_CASE`, `GPS_LOSS` y
  `GPS_COMPLEMENT`, y un único `vcf-complement.jsonl` versionado por conjunto,
  presente solo si existe pérdida. `GPS_CASE` es una lista canónica de casos
  originales y `truth.jsonl` conserva por ALT sus casos fuente; F2 no crea casos
  compuestos ni complementos. La primera cardinalidad positiva del complemento
  se cierra en F5.5 sobre una pérdida real.
- DTC-066 cierra DTPD-103 con `bcftools` 1.24 y HTSlib/tabix 1.24 solo en
  aceptación/CI, comandos fijados y un oráculo semántico independiente. El
  corpus positivo cubre texto, BGZF+TBI y BGZF+CSI; mutaciones aisladas prueban
  cada invariante. EBI `vcf-validator` no es autoridad VCF 4.5.
- DTC-058 resuelve DTPD-117 y concreta versión compatible de Marimo, entorno,
  dependencias, licencia, seguridad, empaquetado Ubuntu/WSL2, contratos y
  adaptador de proceso sin *shell*. La pareja candidata
  exacta es CPython 3.14.7 + Marimo 0.23.16; todavía no está aprobada ni puede
  incorporarse antes de superar la puerta de F2.0.2. `uv` será el único gestor:
  dependencias en `pyproject.toml`, resolución exacta en `uv.lock`, entorno
  local `.venv` excluido de Git y reproducción mediante `uv sync --frozen`,
  sin paquetes globales ni instalación automática de paquetes desde Marimo.
  La app, `pyproject.toml` y `uv.lock` se incluyen en el repositorio junto al lanzador único
  `genome-population-simulator-ui`. Su suborden `install` crea `.venv` y `run`
  no descarga ni actualiza; el entorno no se distribuye y Marimo no se incrusta
  en Rust. La construcción/validación oficial de `.tar.gz` e imagen OCI se
  traslada a F9. El
  enlace es exclusivamente `127.0.0.1`, sin `0.0.0.0`, interfaces LAN, acceso
  remoto ni configuración para ampliar la superficie; WSL2 accede mediante
  `localhost` desde el navegador de Windows. El puerto predeterminado es 2718;
  el lanzador admite `--port <n>` en `1024..=65535`, muestra el puerto efectivo
  y falla ante una colisión con un ejemplo para elegir otro, sin selección
  automática. El
  lanzamiento interactivo imprime siempre la URL y después intenta abrir el
  navegador predeterminado —el de Windows bajo WSL2—; un fallo de apertura solo
  avisa y `--no-open` omite el intento. Abrir la aplicación nunca ejecuta una
  simulación. Estos parámetros quedan fuera del TOML biológico y de
  `EffectiveConfig`. La puerta de versiones debe demostrar resolución y *lock*
  exactos, inventario de dependencias, *builds* y licencias, auditoría de
  seguridad, arranque y controles en Ubuntu 24.04 nativo y WSL2, proceso falso,
  binario Rust real y pruebas del lanzador. Si falla, se reabre DTPD-117 sin
  degradar ni cambiar versiones silenciosamente.
  La concreción incluye además descubrimiento del binario sin *fallback*,
  contratos locales emparejados por SHA-256, protocolo de proceso v1,
  cancelación escalonada, seguridad local, estado de formulario en memoria,
  instantánea validada exacta, registros acotados, accesibilidad WCAG 2.2 AA y
  pruebas automatizadas y manuales. F2.0.2–F2.0.18 materializan y verifican cada parte.
- DTC-054 resuelve DTPD-118 y debe consumirse sin reinterpretación: solo admite
  `alt_af`; separa el espectro, la matriz larga canónica y sus metadatos; fija
  generación bialélica, poblaciones independientes, selección sin reemplazo,
  perfiles explícitos, ámbitos RNG, preparación/reutilización y excepciones
  limitadas. F2.2.1 debe congelar su versión documental antes de implementar sus
  consumidores. La capacidad de panel resuelta por DTC-055 permanece fuera de
  este alcance.
- DTC-052 fija la separación funcional posterior a esa puerta: F2.2 produce o
  admite la matriz poblacional y F2.3 muestrea genotipos HWE o `F` diploide.
  Poisson solo puede obtener un recuento de loci/eventos desde una media o
  densidad; no gobierna frecuencias, alelos, genotipos ni fase.
- DTC-067 resuelve DTPD-120 con `poisson-count-v1`: cero sin consumo, Knuth para
  `0<lambda<10`, PTRS para `10<=lambda<=2^52`, uniformes/consumo exactos y
  `libm` 0.2.16 sujeto a DTC-043. Candidatos no exactos o fuera de `0..=2^53-1`
  fallan; vectores y oráculo independiente cubren límites y Linux/WSL2.
- DTC-059 formaliza las longitudes de *indels* F2: inserción y deleción usan
  parámetros exactos, independientes, predeterminados a 1 pb y dentro de
  `1 <= longitud < umbral estructural efectivo`. No consumen RNG de longitud ni
  admiten una distribución probabilística.
- DTC-068 resuelve DTPD-121 con dos perfiles humanos compactos aproximados. El
  SFS usa los literales 64:12:8 y mínimo `1/5008`; el SNV12 transforma los 96
  contextos Harris 2017 desde el fichero eLife fijado por tamaño/SHA-256. F2 no
  procesa el VCF Phase 3 ni usa `bcftools` para estos builtins. F2.2.14 regenera
  y revisa los recursos pequeños antes de incorporarlos; es una puerta de
  evidencia interna, no una decisión abierta ni un bloqueo para iniciar F2.2.

Que una puerta quede bloqueada no impide preparar matrices, casos literales o
validadores que no dependan de su resultado. Sí impide escribir el contrato de
producción que la consuma.

### 4.3. Puertas de calidad futuras

Estos comandos son criterios de ejecución futura y no se presentan como
ejecutados por este documento:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --workspace --all-targets --locked
cargo test --workspace --doc --locked
cargo build --release --locked --target x86_64-unknown-linux-gnu
cargo metadata --locked --format-version 1
cargo tree --locked -e features
uv sync --frozen
uv run --frozen pytest
```

Cada validador externo elegido debe registrarse con versión y comando exactos.
Los oráculos de aceptación no pueden llamar al productor para calcular el
resultado esperado ni leer el VCF para reconstruir la verdad con la que se lo
compara.

### 4.4. Evidencia ejecutable por paso

F2.0.1 crea `scripts/verify-f2-step` y el registro
`tests/f2/step-evidence.toml`. Cada casilla del plan incluye después una línea
`Evidencia ejecutable` con el comando exacto. El ejecutor corre los comandos y
oráculos registrados para ese identificador, conserva stdout/stderr/versiones y
escribe `evidence/f2/<paso>/result.json` conforme a esquema. El resultado
esperado mínimo es código 0, `status = "passed"` y ninguna comprobación omitida.
Un paso detenido por una decisión pendiente nombrada usa `status = "blocked"`;
una orden ausente, no ejecutada o con esperado no declarado hace fallar el
paso. Un expediente de cobertura parcial conserva el paso como `passed`, pero
usa `status = "partial"` en el registro del REQ, identifica su destino final y
nunca permite superar una pausa como si fuera cierre completo de ese REQ.

## 5. F2.0 — Aplicación Marimo base

**Resultado del hito:** una aplicación local versionada permite configurar,
validar y lanzar todas las capacidades cerradas en F1/F1.6 mediante el mismo
TOML y binario, sin duplicar semántica y con las nueve pestañas preparadas para
las ampliaciones F2–F9.

La preparación de evidencia F2.0.1 se ejecuta antes de F2.0.0 aunque su
numeración conserve el identificador histórico: ningún paso exige
`verify-f2-step` hasta que F2.0.1 lo haya creado y validado.

- ☐ F2.0.0. Ampliar en Rust la frontera mínima F1 para el conjunto de configuración F2 antes de invocar `validate`.
  - Depende de: F2.0.1, F1 cerrado en `e14705a`, cierre documental `efca265` y aprobación futura de inicio F2.
  - Procede de: DTC-034, DTC-036, DTC-037, DTC-071 y DTC-078; REQ-CONF-01, REQ-CONF-07, REQ-CONF-09/10 y REQ-UI-06/07.
  - Cambia: crea `ConfigKind::Scenario`, enlaces general→0..N escenarios, identidad, ruta relativa y digest de cada documento, y preservación opaca; para la vía directa sin documento de escenario fija internamente `scenario_id = "direct"`, reservado, no biológico, no aportable por la persona usuaria ni utilizable por un escenario real, y lo registra en `EffectiveConfig`, manifiesto general y contexto RNG. No activa fases, cruces, eventos ni reglas sin consumidor.
  - Se comprueba con: pruebas Rust para 0/1/N escenarios, enlaces, ruta/digest, campos opacos, rechazo de `scenario_id = "direct"` externo o escenario real con ese ID, y paridad de `EffectiveConfig`/manifiesto/contexto RNG; ninguna ruta alcanza semántica de fases, cruces, eventos o reglas.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.0.0`; debe terminar con código 0 y escribir `evidence/f2/F2.0.0/result.json` con `status = "passed"`.

- ☐ F2.0.1. Crear el ejecutor y registro uniforme de evidencia de F2.
  - Depende de: cierre y commit base de F1, árbol atribuido y aprobación de este plan; no inicia producto F2.
  - Procede de: RAF2-023 y reglas de evidencia de §4.4.
  - Cambia: `scripts/verify-f2-step`, el registro versionado `tests/f2/step-evidence.toml` y el esquema de `evidence/f2/<paso>/result.json`; cada entrada enumera comandos, fixtures, oráculos, versiones y resultados esperados y falla si alguno falta o no se ejecuta.
  - Se comprueba con: `cargo test --test f2_plan_evidence_contract --locked` y `./scripts/verify-f2-step F2.0.1`; ambos terminan 0, el segundo crea `evidence/f2/F2.0.1/result.json` válido y una entrada sin comando, fixture o esperado falla.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.0.1`; debe terminar con código 0 y escribir `evidence/f2/F2.0.1/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.0.2. Evaluar y aprobar experimentalmente las candidatas de DTC-058.
  - Depende de: F2.0.1.
  - Procede de: DTC-043, DTC-048 y DTC-058 (resuelve DTPD-117).
  - Cambia: expediente sin dependencia de producción que evalúa exactamente CPython 3.14.7, Marimo 0.23.16 y `uv` 0.12.0, resolución, *builds*, licencias, mantenimiento, CVE-2026-39987 y compatibilidad básica; termina en aprobación fechada o reapertura de DTPD-117, nunca en sustitución silenciosa.
  - Se comprueba con: `./scripts/evaluate-marimo-toolchain --python 3.14.7 --marimo 0.23.16 --uv 0.12.0 --output evidence/f2/F2.0.2` y `./scripts/verify-f2-step F2.0.2`; inventario, auditoría y decisión fechada son válidos o el paso queda bloqueado sin modificar dependencias de producción.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.0.2`; debe terminar con código 0 y escribir `evidence/f2/F2.0.2/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.0.3. Incorporar el proyecto Python y su entorno bloqueado.
  - Depende de: F2.0.2 aprobada.
  - Procede de: DTC-043, DTC-048 y DTC-058.
  - Cambia: `pyproject.toml`, `uv.lock`, `.python-version`, exclusión de `.venv`, instalación `genome-population-simulator-ui install` con `uv sync --frozen` y prohibición de paquetes globales o instalación automática de Marimo.
  - Se comprueba con: `uv sync --frozen --python 3.14.7`, `uv run --frozen python -VV`, `git check-ignore .venv` y `./scripts/verify-f2-step F2.0.3`; lock divergente falla sin actualizar y las versiones coinciden exactamente.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.0.3`; debe terminar con código 0 y escribir `evidence/f2/F2.0.3/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.0.4. Crear el esqueleto local de la aplicación en el repositorio.
  - Depende de: F2.0.3.
  - Procede de: DTC-042, DTC-045 y DTC-058; REQ-UI-01 y REQ-UI-12.
  - Cambia: aplicación vacía, directorio de contratos y lanzador `genome-population-simulator-ui` en el repositorio; excluye `.venv`, CPython y `uv`, sin incorporar aún controles ni protocolo de F2.0.8/F2.0.14. No construye una entrega oficial.
  - Se comprueba con: instalación desde un clon local mediante `uv sync --frozen`, inspección del árbol versionado y `./scripts/verify-f2-step F2.0.4`; confirma contenidos y exclusiones sin producir `.tar.gz` u OCI.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.0.4`; debe terminar con código 0 y escribir `evidence/f2/F2.0.4/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.0.5. Implementar el ciclo de arranque local y navegador.
  - Depende de: F2.0.3–F2.0.4.
  - Procede de: DTC-058; REQ-UI-08, REQ-UI-09 y REQ-UI-13.
  - Cambia: `genome-population-simulator-ui run` usa `marimo run`, token obligatorio, solo `127.0.0.1`, puerto 2718 o `--port 1024..=65535`, colisión sin selección automática, URL visible, navegador predeterminado/Windows WSL2 y `--no-open`; abrir nunca simula.
  - Se comprueba con: `./scripts/test-ui-launcher --bind 127.0.0.1 --port 2718 --evidence evidence/f2/F2.0.5` y `./scripts/verify-f2-step F2.0.5`; cubre bordes de puerto, colisión, `0.0.0.0` prohibido, apertura/fallo/`--no-open` y ausencia de `run` Rust.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.0.5`; debe terminar con código 0 y escribir `evidence/f2/F2.0.5/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.0.6. Demostrar arranque local sin red después de instalar desde el repositorio.
  - Depende de: F2.0.4–F2.0.5.
  - Procede de: DTC-043 y DTC-058.
  - Cambia: prueba reproducible posterior a instalación desde el repositorio que bloquea red durante el arranque; no construye ni inspecciona OCI.
  - Se comprueba con: `./scripts/test-ui-offline --launcher ./genome-population-simulator-ui --evidence evidence/f2/F2.0.6` y `./scripts/verify-f2-step F2.0.6`; toda descarga o actualización durante el arranque falla.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.0.6`; debe terminar con código 0 y escribir `evidence/f2/F2.0.6/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.0.7. Comprobar de forma local el arranque en Ubuntu/WSL2.
  - Depende de: F2.0.2–F2.0.6.
  - Procede de: DTC-058; REQ-UI-13.
  - Cambia: actas locales separadas de instalación desde repositorio, arranque, navegador, permisos y binario falso/real en Ubuntu Linux y WSL2/ext4. No constituye evidencia oficial de plataforma, distribución ni publicación; DTC-074 y ese gate quedan en F9.
  - Se comprueba con: los recorridos locales identifican su entorno real y registran resultado/limitaciones; `./scripts/verify-f2-step F2.0.7` exige acta automática y manual sin etiquetarla como validación oficial.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.0.7`; debe terminar con código 0 y escribir `evidence/f2/F2.0.7/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.0.8. Crear el primer contrato versionado aplicación↔esquema↔catálogos↔binario.
  - Depende de: F2.0.0–F2.0.7.
  - Procede de: REQ-UI-01, REQ-UI-09 y REQ-UI-12; DTC-036, DTC-051 y DTC-073 (resuelve RAF2-034).
  - Cambia: primer directorio versionado con manifiesto, esquemas, catálogo y `capabilities --format json`, versiones y SHA-256; no extiende un catálogo F1 inexistente. `--simulator-bin` es una selección estricta. El protocolo `"1"` separa `catalogued_authorizations` de `available_authorizations`; en F2.0 ambas opciones están catalogadas y ninguna disponible. Marimo solo muestra la intersección y registra ruta, versión y SHA-256 en sesión —nunca TOML/`EffectiveConfig`—.
  - Se comprueba con: rutas/fallos sin *fallback*; `available` no subconjunto, duplicado, opción desconocida o disponibilidad sin parser/semántica fallan; una opción solo catalogada no aparece ni se transmite. También cubre protocolo frente a *semver*, contratos compatible/anterior/posterior/ausente/corrupto/intercambiado y modo diagnóstico sin búsqueda externa ni mezcla.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.0.8`; debe terminar con código 0 y escribir `evidence/f2/F2.0.8/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.0.9. Crear el armazón accesible de las nueve pestañas y su visibilidad por vía/capacidad.
  - Depende de: F2.0.1–F2.0.8.
  - Procede de: REQ-UI-02 y REQ-UI-13; DTC-051.
  - Cambia: aplicación servida solo con `marimo run`, token obligatorio y sin editor, terminal, IA, telemetría, compartición o instalador; contiene `General/Ejecución`, `Genoma`, `Individuos`, `Poblaciones`, `Evolución`, `Familias`, `Eventos técnicos`, `Salidas` y `Revisar y ejecutar`, con indisponibilidad explicada.
  - Se comprueba con: presencia/orden estable, vías direct/evolution/families, token y rutas negativas `edit`/terminal/IA/red, navegación por teclado y comprobación de que una pestaña no disponible no pierde su destino trazable.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.0.9`; debe terminar con código 0 y escribir `evidence/f2/F2.0.9/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.0.10. Mantener un único estado de formulario tipado e invalidar toda validación tras una edición.
  - Depende de: F2.0.8–F2.0.9.
  - Procede de: REQ-UI-03, REQ-UI-07 y REQ-UI-09; arquitectura §§6.1 y 7.1.
  - Cambia: modelo de autoría multidocumento en memoria, sin autoguardado, recuperación ni estado biológico del navegador; `Validar` crea una instantánea TOML visible en temporal privado 0700 y ficheros 0600, la asocia a su digest y `Ejecutar` reutiliza esos bytes; solo exportar persiste configuración.
  - Se comprueba con: editar/añadir/quitar/reordenar/renombrar invalida instantánea y ejecución; reiniciar pierde cambios no exportados; permisos, mismo UID y limpieza segura de temporales; la app nunca calcula regla ganadora, valor efectivo o verdad.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.0.10`; debe terminar con código 0 y escribir `evidence/f2/F2.0.10/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.0.11. Cubrir con controles tipados solo referencia, secuencias y reducciones cerradas en F1.
  - Depende de: F2.0.8–F2.0.10.
  - Procede de: REQ-UI-01, REQ-UI-03 y REQ-CONF-124 a REQ-CONF-127; DTC-050/DTC-051.
  - Cambia: controles de clase/versión, seed, referencia simulada/suministrada/reducida, IDs, clases, tamaños, `region_size_bases`, unión de objetivo y opciones de salida F1, con unidad normativa en bases y ayuda en kb. Vía, reglas y eventos permanecen indisponibles.
  - Se comprueba con: matriz opción F1/F1.6→pestaña/control/campo TOML, límites y exclusividad de reducción, conversión visual kb↔bases sin pérdida y ausencia de opción implementada solo en TOML.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.0.11`; debe terminar con código 0 y escribir `evidence/f2/F2.0.11/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.0.12. Implementar tablas editables solo para secuencias y reducciones F1.
  - Depende de: F2.0.9–F2.0.11.
  - Procede de: REQ-UI-04, REQ-UI-05 y DTC-051.
  - Cambia: tablas con añadir/eliminar/reordenar, IDs y orden visibles de secuencias y reducciones. Vía, reglas, selectores, eventos y sus listas no existen aún y se muestran como indisponibles hasta su primer consumidor; no concede precedencia por orden accidental de celdas.
  - Se comprueba con: alta/baja/reordenación, duplicados y comparación del TOML F1 resultante con fixtures manuales equivalentes; las rutas de regla, selector y evento permanecen ausentes.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.0.12`; debe terminar con código 0 y escribir `evidence/f2/F2.0.12/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.0.13. Cerrar la importación/exportación atómica del directorio autocontenido DTC-069.
  - Depende de: F2.0.13.1–F2.0.13.4.
  - Procede de: REQ-UI-06 y REQ-UI-12; DTC-034, DTC-036, DTC-051 y DTC-069; RAF2-017/036 resueltos.
  - Cambia: integra inventario, copia aislada, preservación opaca y exportación no-clobber como una única puerta de aceptación, sin volver a mezclar sus responsabilidades.
  - Se comprueba con: ejecución conjunta de los cuatro expedientes hijos y prueba de que ningún fallo adopta una sesión ni publica un directorio parcial.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.0.13`; debe terminar con código 0 y escribir `evidence/f2/F2.0.13/result.json` con `status = "passed"` y los cuatro expedientes hijos enlazados.

- ☐ F2.0.13.1. Inventariar y validar el grafo autocontenido y sus límites.
  - Depende de: F2.0.8 y F2.0.10–F2.0.12.
  - Procede de: DTC-069; contrato de importación de REQ-UI-06/12.
  - Cambia: selector de un directorio, inventario canónico de exactamente un TOML general y 0..N escenarios alcanzables, validación de enlaces relativos y límites antes de abrir contenido.
  - Se comprueba con: general ausente/múltiple, enlace ausente/duplicado, ruta absoluta/`..`/no canónica, ciclo, escape, extra y cada borde de cantidad/tamaño/profundidad DTC-069.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.0.13.1`; debe terminar con código 0 y escribir `evidence/f2/F2.0.13.1/result.json` con inventario, límites ejercitados y `status = "passed"`.

- ☐ F2.0.13.2. Copiar el conjunto validado al espacio privado con controles TOCTOU proporcionales de DTC-069 y sin escapes.
  - Depende de: F2.0.13.1.
  - Procede de: DTC-069; aislamiento local de DTC-051.
  - Cambia: apertura segura, rechazo de symlink/hardlink que rompa aislamiento, copia privada atómica y comparación de identidad/hash entre inventario, copia y adopción.
  - Se comprueba con: sustitución y carrera entre inventario/apertura/copia/hash, enlace inseguro, cambio de inode/contenido y permisos; todo fallo descarta el provisional y conserva intacta la sesión anterior.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.0.13.2`; debe terminar con código 0 y escribir `evidence/f2/F2.0.13.2/result.json` con carreras reproducibles y ausencia de adopción parcial.

- ☐ F2.0.13.3. Preservar opacamente documentos y campos válidos sin controles.
  - Depende de: F2.0.13.2.
  - Procede de: DTC-069; autoría multidocumento F2.0.10.
  - Cambia: modelo de round-trip que conserva bytes/documentos válidos aún no editables y solo invalida el digest cuando el usuario realiza un cambio explícito.
  - Se comprueba con: 0/1/N escenarios, espacios, campos desconocidos válidos y documentos alcanzables sin control; importación-exportación sin edición es exacta y una edición invalida el digest sin perder el resto opaco.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.0.13.3`; debe terminar con código 0 y escribir `evidence/f2/F2.0.13.3/result.json` con hashes antes/después y `status = "passed"`.

- ☐ F2.0.13.4. Exportar otro directorio completo de forma atómica y no-clobber.
  - Depende de: F2.0.13.1–F2.0.13.3.
  - Procede de: DTC-069; publicación segura DTC-045.
  - Cambia: exportación mediante directorio provisional, rutas relativas canónicas, revalidación, hashes y publicación no-clobber; no admite ZIP/TAR ni selección manual múltiple.
  - Se comprueba con: destino inexistente/existente, carrera de creación, fallo de escritura/revalidación/hash y grafo alterado; no queda salida parcial y el original importado nunca se modifica.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.0.13.4`; debe terminar con código 0 y escribir `evidence/f2/F2.0.13.4/result.json` con publicación atómica demostrada.

- ☐ F2.0.14. Invocar `validate` sin shell y asociar diagnósticos a campos sin reinterpretarlos.
  - Depende de: F2.0.1–F2.0.8 y F2.0.13.
  - Procede de: REQ-UI-07, REQ-UI-09 y REQ-VAL-54; DTC-044 y DTC-051.
  - Cambia: protocolo `"1"` sin *shell* para `validate <general> --diagnostics json [autorizaciones]`, cwd del general, un proceso activo, `stdin` cerrado y sin reintentos; el digest incluye documentos, contratos y autorizaciones ordenadas. Confirmar o retirar una autorización invalida y repite `validate`. `stdout`/`stderr` se drenan concurrentemente; el primero queda vacío o contiene un JSON y el segundo solo log.
  - Se comprueba con: proceso falso registra argv/cwd y rechaza opción no catalogada/general; espacios/metacaracteres, códigos 0/1/2/3/4/5/70, JSON ausente solo en éxito limpio, JSON inválido/doble/incompatible, ambos flujos llenos/intercalados sin bloqueo y prueba de que stderr/texto nunca se ejecuta ni adquiere semántica.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.0.14`; debe terminar con código 0 y escribir `evidence/f2/F2.0.14/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.0.15. Seleccionar el destino de forma segura y lanzar `run` solo por acción explícita.
  - Depende de: F2.0.10, F2.0.14 y publicación F1.4.
  - Procede de: REQ-UI-08 a REQ-UI-10; DTC-045 y DTC-051.
  - Cambia: destino/no-clobber y botón habilitado solo para el digest vigente; lanza una vez `run <general> --diagnostics json --output <destino> [autorizaciones]` con exactamente documentos, contrato y autorizaciones validados. Crea un grupo dedicado, no reutiliza hijo y cancelar envía solo a ese grupo SIGINT/10 s, SIGTERM/5 s y SIGKILL; registra etapas y nunca completa salida parcial.
  - Se comprueba con: destino vacío/existente/válido, doble clic, edición/recarga y autorización añadida/retirada; proceso falso/real, descendientes que ignoran señales, drenaje hasta EOF, espera sin huérfanos y proceso ajeno intacto; manifiesto final ausente tras código 0 falla y Rust conserva no-clobber.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.0.15`; debe terminar con código 0 y escribir `evidence/f2/F2.0.15/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.0.16. Presentar estado, logs, diagnósticos y resultados sin cargar artefactos grandes.
  - Depende de: F2.0.14–F2.0.15 y manifiestos F1.4/F1.6.
  - Procede de: REQ-UI-11 y DTC-051.
  - Cambia: ficheros 0600 bajo estado XDG, cola UI de 1 MiB y paginación; límites 16 MiB stdout/64 MiB stderr por invocación y 256 MiB total, exceso termina sin truncar, token redactado y limpieza explícita/evicción de sesiones terminadas; resumen solo desde manifiestos/metadatos.
  - Se comprueba con: bordes y exceso de los tres límites, proceso falso que llena ambos flujos intercalados sin bloquear, sesión activa nunca eliminada, token ausente, proceso largo/fallido, manifiesto válido/roto y artefacto señuelo grande que demuestra lectura parcial.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.0.16`; debe terminar con código 0 y escribir `evidence/f2/F2.0.16/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.0.17. Cerrar accesibilidad y funcionamiento local Ubuntu/WSL2.
  - Depende de: F2.0.17.1–F2.0.17.4.
  - Procede de: REQ-UI-13; DTC-042, DTC-051, DTC-058 y DTC-074; RAF2-035/036 resueltos.
  - Cambia: reúne cuatro expedientes independientes —navegador, accesibilidad automática, teclado/foco y lector de pantalla/plataformas— bajo una puerta WCAG 2.2 AA de controles propios.
  - Se comprueba con: los cuatro expedientes hijos aprobados, limitaciones de Marimo registradas y guía coherente con `127.0.0.1`, navegador Windows desde WSL2 y alternativa manual.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.0.17`; debe terminar con código 0 y escribir `evidence/f2/F2.0.17/result.json` con `status = "passed"` y enlaces a F2.0.17.1–F2.0.17.4.

- ☐ F2.0.17.1. Fijar y verificar el navegador automatizado reproducible.
  - Depende de: F2.0.2 y F2.0.9–F2.0.16.
  - Procede de: DTC-058 y DTC-074; puerta de navegador de REQ-UI-13.
  - Cambia: Playwright y axe fijados en `uv.lock`; un Chromium único con versión, URL, artefacto y SHA-256 inventariados, más receta Ubuntu exacta y ejecución posterior sin red.
  - Se comprueba con: instalación reproducible de Chromium y dependencias desde receta aprobada, segunda ejecución offline, rechazo de hash/versión distintos y registro de la plataforma detectada.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.0.17.1`; debe terminar con código 0 y escribir `evidence/f2/F2.0.17.1/result.json` con inventario y repetición offline.

- ☐ F2.0.17.2. Ejecutar la auditoría automática de accesibilidad.
  - Depende de: F2.0.17.1.
  - Procede de: REQ-UI-13; objetivo WCAG 2.2 AA de controles propios.
  - Cambia: corpus automatizado de etiquetas/ayuda, errores asociados, contraste, zoom 200 %, estados no basados solo en color y rutas principales de la aplicación.
  - Se comprueba con: Playwright+axe sobre estados vacío, válido, inválido, validando, ejecutando y fallido; toda infracción o ruta no visitada bloquea.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.0.17.2`; debe terminar con código 0 y escribir `evidence/f2/F2.0.17.2/result.json` con reglas, rutas, versión y resultados completos.

- ☐ F2.0.17.3. Revisar manualmente teclado, orden de foco y recuperación tras errores.
  - Depende de: F2.0.17.2.
  - Procede de: REQ-UI-13; controles y estados implementados en F2.0.9–F2.0.16.
  - Cambia: acta manual separada para navegación sin ratón, foco visible, orden lógico, modales, cancelación y retorno del foco tras diagnóstico.
  - Se comprueba con: guion versionado ejecutado en Ubuntu 24.04 nativo y WSL2/navegador Windows, con fecha, revisor, versiones, pasos y resultado; una omisión o limitación no aceptada bloquea.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.0.17.3`; debe terminar con código 0 y escribir `evidence/f2/F2.0.17.3/result.json` que enlace ambas actas y su plataforma detectada.

- ☐ F2.0.17.4. Revisar manualmente lector de pantalla y funcionamiento por plataforma.
  - Depende de: F2.0.17.1–F2.0.17.3 y F2.0.7.
  - Procede de: DTC-042 y DTC-074; REQ-UI-13.
  - Cambia: actas independientes de lector de pantalla/navegador para Ubuntu nativo y WSL2 con navegador Windows, sin presentar Windows nativo ni acceso remoto como soportados.
  - Se comprueba con: anuncio de etiquetas, ayuda, errores, cambios de estado y resultados; cada acta registra lector/navegador/versiones, evidencia de plataforma detectada, fecha, revisor y limitaciones.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.0.17.4`; debe terminar con código 0 y escribir `evidence/f2/F2.0.17.4/result.json` con las dos actas enlazadas y `status = "passed"`.

- ☐ F2.0.18. Completar equivalencia y trazabilidad de F2.0.
  - Depende de: F2.0.1–F2.0.17.
  - Procede de: REQ-UI-01 a REQ-UI-13, REQ-VAL-54 y ACE-114/115.
  - Cambia: 13 casos `TEST-REQ-*` de cierre exclusivo, un expediente parcial para REQ-UI-13 con destino F9.7, matriz opción→control, pruebas UI→TOML→binario y expediente DTPD-117; inventario del repositorio confirma la UI y no marca controles futuros como implementados. La entrega oficial queda en F9.
  - Se comprueba con: par UI/manual para referencia, secuencias y reducciones F1, y general+0..N escenarios con igual `EffectiveConfig`/diagnóstico; puertas Python/Rust/navegador, arranque local/offline, recuento de 13 REQ exclusivos sin huérfanos y expediente REQ-UI-13 `status = "partial"` con destino F9.7. `.tar.gz`/OCI y evidencia oficial quedan fuera de F2.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.0.18`; debe terminar con código 0 y escribir `evidence/f2/F2.0.18/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

> **PAUSA-F2.0 — cierre obligatorio de F2.0.** No iniciar F2.1 hasta superar y aprobar la puerta
> experimental F2.0.2 de DTC-058, la equivalencia de referencia, secuencias y reducciones F1, el
> emparejamiento de contratos, la invocación sin *shell* y la cancelación, la
> autoría/exportación autoconsistente del conjunto TOML, la invalidación tras
> cambiar cualquier documento/enlace, el destino no-clobber, los límites de
> registros, la seguridad local y la accesibilidad. Las entregas diferenciadas
> `.tar.gz`/OCI y la evidencia oficial de DTC-074 pertenecen a F9. RAF2-001 y RAF2-027 están resueltos; este último mediante el
> cierre de F1 en `e14705a`. RAF2-004 y RAF2-007 se resuelven
> mediante las reasignaciones y pruebas negativas explícitas de este plan.

## 6. F2.1 — Estado haplotípico mínimo

**Resultado del hito:** coordenadas y estado propios representan cada copia
haplotípica, fase, dosis, procedencia y `N`, y permiten materializar secuencias
exactas desde la referencia cerrada en F1.

- ☐ F2.1.1. Congelar y auditar la línea base que F2 consumirá.
  - Depende de: PAUSA-F2.0 aprobada, cierre explícito de F1.4, F1.5 y F1 y commit base identificado.
  - Procede de: dependencia F1→F2 de la hoja de ruta; DTC-043, DTC-045 y DTC-048.
  - Cambia: evidencia de preflight de F2 con commit base, toolchain, lock, esquemas, catálogos, comandos verdes y contratos F1 reutilizados; no cambia producción.
  - Se comprueba con: revisión del acta de cierre F1, `git status --short`, `rustc --version`, `cargo metadata --locked --format-version 1` y las puertas de §4.3 sobre el commit base.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.1.1`; debe terminar con código 0 y escribir `evidence/f2/F2.1.1/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.1.2. Auditar y congelar para F2 la DTPD-101 resuelta en F1.6.
  - Depende de: F2.1.1.
  - Procede de: DTPD-101; DTC-010; REQ-CONF-28, REQ-IND-38 y REQ-IND-39.
  - Cambia: expediente de consumo F2 que enlaza la decisión y los tipos/servicio comunes de F1.6, completa vectores de inserciones ancladas, deleciones internas y composición y prohíbe otra representación; no reabre opciones ni cambia DTPD-101.
  - Se comprueba con: huella/versión coincidentes con F1.6, vectores literales de primera/última base y rangos vacíos, ida/vuelta en fronteras y revisión independiente que detectaría cualquier conversión o tipo local nuevo.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.1.2`; debe terminar con código 0 y escribir `evidence/f2/F2.1.2/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.1.3. Congelar y auditar DTC-061 como resolución completa de DTPD-102.
  - Depende de: F2.1.1 y F2.1.2 para interpretar coordenadas del prototipo.
  - Procede de: DTPD-102 y DTC-061; arquitectura §§5.3, 6.2, 9 y 10; REQ-IND-28, REQ-IND-29 y REQ-VER-31.
  - Cambia: expediente de consumo que congela modelo híbrido y `truth.jsonl` v1 sin comparar alternativas: UTF-8 sin comprimir, cabecera, tipos de registro, orden, recuentos/digests, enlaces a fuentes, integridad, flujo y dirección única verdad→VCF.
  - Se comprueba con: recorrido DTC-061→arquitectura→pasos, prototipo JSONL manual, esquemas por tipo y validador semántico ajeno que detectan copia/genotipo discordantes, relaciones ausentes/duplicadas/desordenadas y recuentos/digests falsos; ninguna prueba reelige formato.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.1.3`; debe terminar con código 0 y escribir `evidence/f2/F2.1.3/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.1.4. Evaluar cada dependencia necesaria para las decisiones F2.1.2–F2.1.3 antes de incorporarla.
  - Depende de: F2.1.2–F2.1.3; se omite si los contratos aprobados no necesitan dependencia nueva.
  - Procede de: DTC-043, DTC-047 y DTC-048.
  - Cambia: evaluación individual con propósito, REQ, versión exacta, licencia, mantenimiento, seguridad, características, árbol, código nativo, scripts y coste de distribución; después, solo las dependencias aprobadas y `Cargo.lock`.
  - Se comprueba con: revisión de la evaluación, `cargo metadata --locked --format-version 1`, `cargo tree --locked -e features` y repetición de §4.3.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.1.4`; debe terminar con código 0 y escribir `evidence/f2/F2.1.4/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.1.5. Extender el servicio común de coordenadas de F1.6 a los casos de F2.
  - Depende de: F2.1.2 y F2.1.4 cuando corresponda.
  - Procede de: DTC-010, DTC-013 y DTC-039; REQ-CONF-28, REQ-IND-38 y REQ-IND-39.
  - Cambia: los tipos y el servicio ya usados por el mapa de reducción ganan constructores validados para posición, intervalo, ancla de inserción y deleción interna; no expone tipos de crates ni crea un segundo módulo/convenio.
  - Se comprueba con: pruebas de tabla compartidas con F1.6 desde la semántica externa base 1 cerrada, propiedades de ida/vuelta y fallos explícitos en cero, fuera de rango, desbordamiento y límites terminales.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.1.5`; debe terminar con código 0 y escribir `evidence/f2/F2.1.5/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.1.6. Definir bases y vistas de referencia que preserven `N` literalmente.
  - Depende de: F2.1.1, F2.1.5 y contrato FASTA cerrado en F1.5.
  - Procede de: DTC-032; REQ-REF-14 y REQ-IND-28.
  - Cambia: tipos de dominio para base y acceso inmutable por secuencia/coordenada, reutilizando la referencia F1 sin reinterpretar ni sustituir `N`.
  - Se comprueba con: referencia sintética que contiene `N` en principio, interior y final, comparación byte a byte con la referencia F1 y casos de error de secuencia o coordenada inexistente.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.1.6`; debe terminar con código 0 y escribir `evidence/f2/F2.1.6/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.1.7. Modelar individuo, secuencia y copia haplotípica con ploidía uniforme explícita.
  - Depende de: DTC-061, F2.1.5–F2.1.6 y catálogo de IDs cerrado en F1.4.
  - Procede de: DTC-003, DTC-039 y DTC-061; REQ-IND-28, REQ-IND-29 y REQ-VER-27; alcance deliberado de F2.
  - Cambia: tipos de estado con exactamente una copia por índice de fase entre 1 y la ploidía 1–10; cada copia posee una cadena ordenada, lógicamente inmutable y compartible de segmentos, sin habilitar ploidía distinta por secuencia ni fijar un crate/contenedor concreto.
  - Se comprueba con: constructores y pruebas para ploidías 1, 2 y 10, rechazo de 0 y 11, unicidad de IDs/índices, cardinalidad exacta y compartición de segmentos intactos sin clonar sus bases.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.1.7`; debe terminar con código 0 y escribir `evidence/f2/F2.1.7/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.1.8. Representar ediciones y procedencia sin duplicación innecesaria de la referencia.
  - Depende de: DTC-061 y F2.1.5–F2.1.7.
  - Procede de: DTC-003, DTC-031, DTC-039 y DTC-061; REQ-IND-28 y REQ-IND-29.
  - Cambia: segmentos ordinarios referencian intervalos de referencia/fundador; SNP e *indels* pequeños son modificaciones tipadas con alelo verdadero, origen explícito/aleatorio, copia/fase y coordenadas DTC-056. Solo una inserción sin fuente conserva sus propios bytes; no añade clases fuera de F2.
  - Se comprueba con: copia intacta, SNP, inserción con/sin fuente y deleción; identidad/procedencia/coordenadas por modificación, ausencia de bytes duplicados y prototipo de composición futura que conserva modificaciones al cortar/concatenar segmentos.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.1.8`; debe terminar con código 0 y escribir `evidence/f2/F2.1.8/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.1.9. Materializar una copia y calcular fase y dosis desde el estado canónico.
  - Depende de: F2.1.7–F2.1.8.
  - Procede de: DTC-003 y DTC-031; REQ-IND-28, REQ-IND-29 y REQ-VER-27.
  - Cambia: operaciones de materialización determinista y consulta genotípica por locus; fase y dosis se derivan de copias ordenadas, no de una salida VCF.
  - Se comprueba con: oráculo manual que aplica ediciones a cadenas cortas, casos haploide/diploide/poliploide y propiedad de que dosis equivale al número exacto de copias portadoras.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.1.9`; debe terminar con código 0 y escribir `evidence/f2/F2.1.9/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.1.10. Completar la evidencia y trazabilidad de F2.1.
  - Depende de: F2.1.5–F2.1.9 y puertas de calidad verdes.
  - Procede de: F2.1; REQ-REF-14, REQ-IND-28, REQ-IND-29 y REQ-VER-27; ACE-03, ACE-10, ACE-11, ACE-28 y ACE-66.
  - Cambia: pruebas de integración, matriz de evidencia F2.1, casos de aceptación y documentación de límites frente a F5.3.
  - Se comprueba con: §4.3, recuento automático de tres REQ exclusivos sin huérfanos, verificación de `REQ-REF-14` como invariante heredada de F1.5 y revisión cruzada de materialización, fase, dosis, procedencia y conservación de `N`.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.1.10`; debe terminar con código 0 y escribir `evidence/f2/F2.1.10/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

> **PAUSA-F2.1 — cierre obligatorio de F2.1.** Revisar el consumo sin reapertura de DTPD-101 y
> DTPD-102, auditar sin reabrir el modelo interno ni `truth.jsonl` de DTC-061,
> la extensión de tipos de coordenada, representación
> haplotípica, consumo de memoria, oráculos y los
> tres REQ exclusivos y la invariante heredada `REQ-REF-14`. No iniciar ningún paso de F2.2, F2.3 o F2.4 si la representación no reconstruye
> exactamente cada copia o si una conversión de coordenadas queda implícita.

## 7. F2.2 — Variantes pequeñas y matriz de frecuencias SNV

**Resultado del hito:** el motor admite o genera una matriz poblacional de
frecuencias SNV reproducible y validada, aplica SNP explícitos e *indels*
pequeños inequívocos y nunca usa Poisson para muestrear frecuencias, alelos o
genotipos.

- ☐ F2.2.1. Congelar como entrada de implementación DTC-054, que resuelve DTPD-118.
  - Depende de: PAUSA-F2.1 aprobada, F2.1.2 y F2.1.3.
  - Procede de: DTC-052, DTC-054 y DTPD-118 resuelta; REQ-CONF-128, REQ-CONF-129, REQ-IND-60 y REQ-IND-61.
  - Cambia: expediente de consumo que enlaza la versión aprobada y enumera formato, algoritmos, perfiles, límites, ámbitos aleatorios y excepciones; declara fuera de F2 VCF→TSV, generación multialélica, estructura/evolución poblacional y la capacidad de panel de DTC-055.
  - Se comprueba con: revisión requisito→cláusula DTC-054→paso que no deja comportamiento aprobado sin consumidor ni reabre DTPD-118 al implementar.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.1`; debe terminar con código 0 y escribir `evidence/f2/F2.2.1/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.2. Congelar DTC-067 y DTC-059 antes de sus consumidores.
  - Depende de: PAUSA-F2.1 aprobada y F2.1.2–F2.1.3; es independiente de la elección de formato F2.2.1.
  - Procede de: DTPD-120 resuelta; DTC-022, DTC-037, DTC-043, DTC-059, DTC-067 y DTC-072; REQ-CONF-23 a REQ-CONF-25 y REQ-IND-16 a REQ-IND-18.
  - Cambia: expediente de consumo y evaluación DTC-043 de `libm = =0.2.16`; implementa únicamente `poisson-count-v1`, con ramas/fronteras, límites, errores y metadatos exactos, y deja disponible la dependencia matemática aprobada para F2.2.8.4. Congela DTC-059 sin RNG de longitud.
  - Se comprueba con: licencia/auditoría/`build.rs`/lock de libm; vectores Poisson, hook fuera de rango, PMF/media/varianza y bytes iguales Linux/WSL2; longitudes 1/2/`umbral-1` y rechazos 0/`>=umbral`. Los vectores de `log-uniform-v1` pertenecen exclusivamente a F2.2.8.4.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.2`; debe terminar con código 0 y escribir `evidence/f2/F2.2.2/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.3. Publicar los esquemas versionados de espectro y matriz SNV.
  - Depende de: F2.2.1.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-CONF-128, REQ-CONF-129 y REQ-IND-60.
  - Cambia: esquema del espectro TSV por intervalos y filas `fixed`; entrada suministrada mínima `contig/position/ref/alt/population_id/alt_af` con `source_id` opcional; pareja canónica `population_snvs.tsv` + `population_snvs.metadata.json`, sin entrada VCF.
  - Se comprueba con: fixtures manuales válidos e inválidos validados por una implementación de esquema separada, rechazo de columnas de frecuencia distintas de `alt_af` y prueba negativa de que `run` no acepta VCF ni TSV sin metadatos canónicos.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.3`; debe terminar con código 0 y escribir `evidence/f2/F2.2.3/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.4. Fijar el contrato de metadatos, integridad y referencia de la pareja canónica.
  - Depende de: F2.2.3 y contratos SHA-256 de F1.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-CONF-128, REQ-IND-60, REQ-VAL-55 y REQ-MAN-91.
  - Cambia: JSON con versión, referencia esperada, SHA-256 del FASTA, coordenadas `1-based-closed`, poblaciones, procedencia sintética/empírica y SHA-256 del TSV; muestra/fuente/filtros desconocidos se declaran `unknown`, y configuración/manifiesto enlazan el SHA-256 del JSON sin digest circular.
  - Se comprueba con: pareja escrita a mano, `sha256sum`, intercambio independiente de TSV/JSON/FASTA, contig literal `1` frente a `chr1` y rechazo previo a individuos ante cualquier discordancia.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.4`; debe terminar con código 0 y escribir `evidence/f2/F2.2.4/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.5. Especificar la autoría TOML excluyente de matriz preparada, reutilizada o generada y de los demás eventos F2.
  - Depende de: F2.2.2–F2.2.4.
  - Procede de: REQ-CONF-08, REQ-CONF-21, REQ-CONF-24, REQ-CONF-25, REQ-CONF-28, REQ-CONF-128 y REQ-CONF-129; DTC-034 y DTC-052.
  - Cambia: esquema documental para fuente canónica o generación, referencia/poblaciones, `exact_loci`/`mean_loci`/`density_per_mb`, espacio, sustituciones y espectro explícito, más SNP e *indels*; elimina densidad SNP inicial por individuo.
  - Se comprueba con: ejemplos de las tres rutas y rechazo de fuentes simultáneas, perfil alélico ausente, más de un modo de recuento, densidad SNP individual, VCF y campos futuros.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.5`; debe terminar con código 0 y escribir `evidence/f2/F2.2.5/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.6. Introducir la autoría F2 de tolerancia general y las primeras reglas/formas necesarias.
  - Depende de: F2.2.3–F2.2.5.
  - Procede de: DTC-034, DTC-036, DTC-044 y DTC-052; REQ-CONF-128, REQ-CONF-129 y REQ-DIR-03.
  - Cambia: `DeclaredConfig`, adaptador TOML y diagnósticos para tolerancia general, matriz, metadatos, espectros, regiones, indels y población de origen; conserva literales y ubicaciones sin construir verdad, calcular digests ni normalizar silenciosamente, y no presupone un resolvedor o `ParameterKey` F1.
  - Se comprueba con: parsing de cada forma, claves desconocidas/duplicadas, clase/versión/ruta errónea y snapshots de documento, campo, literal y ubicación.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.6`; debe terminar con código 0 y escribir `evidence/f2/F2.2.6/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.6.1. Crear el primer normalizador real de DTC-049 para colecciones F2.
  - Depende de: F2.2.6.
  - Procede de: DTC-049, DTC-078, REQ-CONF-106/107, REQ-VAL-40 y REQ-MAN-81.
  - Cambia: crea el único normalizador F2 de probabilidades/frecuencias con orden canónico, Neumaier, tolerancia efectiva, antes/después y diagnósticos; es el primer consumidor real de DTC-049 y registra en `EffectiveConfig` y manifiesto la tolerancia, suma, delta, decisión y valores antes/después. F6.2 reutilizará este contrato para colecciones poblacionales posteriores, sin atribuirlo retrospectivamente a F1.
  - Se comprueba con: sumas dentro/fuera de tolerancia, cero/no finita, borde inclusivo, registro completo de `EffectiveConfig`/manifiesto y oráculo numérico independiente; expediente parcial de los cuatro REQ con destino F6.2.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.6.1`; debe terminar con código 0 y escribir `evidence/f2/F2.2.6.1/result.json` con `status = "passed"`.

- ☐ F2.2.7. Crear el primer catálogo y dominio real de locus, posición y espacio.
  - Depende de: F2.2.5–F2.2.6.1.
  - Procede de: DTC-006, DTC-033, DTC-049 y DTC-052; REQ-IND-10, REQ-IND-13, REQ-IND-14 y REQ-CONF-129.
  - Cambia: crea el catálogo/dominio/resolvedor real de locus, posición y espacio para recuentos excluyentes, sustituciones, perfil, modo `exact/stochastic` y longitudes exactas de inserción/deleción, con parámetros efectivos y procedencia; no introduce jerarquías futuras ni `resolution-evidence`.
  - Se comprueba con: catálogo completo, dominios/unidades, cero explícito, precedencia espacial, valor efectivo y rechazo de ejes no admitidos o fallback empírico implícito.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.7`; debe terminar con código 0 y escribir `evidence/f2/F2.2.7/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.8. Cerrar el catálogo RNG v4 y la familia de consumidores F2.
  - Depende de: F2.2.8.1–F2.2.8.6.
  - Procede de: DTC-011, DTC-037, DTC-038, DTC-052 y DTC-070 a DTC-075; REQ-IND-61, REQ-IND-63 y REQ-VAL-55; RAF2-028/033/041–043 resueltos.
  - Cambia: integra catálogo, primitivas, selectores, distribuciones, genotipo/fase y oráculos como una puerta única, conservando un solo productor por algoritmo.
  - Se comprueba con: aprobación de los seis expedientes hijos y prueba de que ningún consumidor se compila o habilita contra una versión parcial de catálogo/algoritmos.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.8`; debe terminar con código 0 y escribir `evidence/f2/F2.2.8/result.json` con `status = "passed"` y enlaces a F2.2.8.1–F2.2.8.6.

- ☐ F2.2.8.1. Publicar el catálogo `subseed-scopes-v4.tsv` y sus vectores de derivación.
  - Depende de: F2.2.1, F2.2.7 y catálogo RNG cerrado en F1.3.
  - Procede de: DTC-011, DTC-037, DTC-038, DTC-070, DTC-071 y DTC-075.
  - Cambia: conserva literalmente las seis filas/vectores v3 y añade los once ámbitos F2 con componentes exactos, `catalog_version=4`, `derivation_version=1` y algoritmo consumidor; ALT usa escenario/referencia/secuencia/posición/REF y `variant_id` nace después.
  - Se comprueba con: parser ajeno, igualdad byte a byte v3⊂v4, vectores BLAKE3 por ámbito/componente, mutación aislada y rechazo de componente/orden/entero no canónico; ningún ámbito depende de ALT antes de sortearlo.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.8.1`; debe terminar con código 0 y escribir `evidence/f2/F2.2.8.1/result.json` con catálogo, vectores y `status = "passed"`.

- ☐ F2.2.8.2. Implementar el uniforme abierto y el entero acotado normativos.
  - Depende de: F2.2.8.1.
  - Procede de: DTC-037, DTC-071, DTC-072 y DTC-075.
  - Cambia: productores únicos `open-unit-f64-v1` y `uniform-u64-below-v1`, con conversión de bits, rechazo, fronteras y consumo cero exactos.
  - Se comprueba con: vectores literales de `next_u64`, 0/1/máximo, `n=0/1`, rechazo forzado, consumos exactos y oráculo independiente de aritmética.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.8.2`; debe terminar con código 0 y escribir `evidence/f2/F2.2.8.2/result.json` con vectores y `status = "passed"`.

- ☐ F2.2.8.3. Implementar categórico, permutación y selectores estático/dinámico.
  - Depende de: F2.2.6.1 y F2.2.8.1–F2.2.8.2.
  - Procede de: DTC-020, DTC-049, DTC-071, DTC-072 y DTC-075; ACE-107.
  - Cambia: productores únicos `categorical-f64-v1`, `fisher-yates-v1`, `weighted-without-replacement-v1` solo para SNV estáticos y `dynamic-indel-selection-v1` para reenumeración del estado vigente; este último consume solo `indel_position`, nunca `insertion_content`.
  - Se comprueba con: enumeración exhaustiva pequeña, pesos degenerados, permutaciones y selección SNV; para indels, cadenas manuales que ganan/pierden oportunidades, agotamiento, orden canónico, mismo flujo continuo, consumo cero, aislamiento respecto a `insertion_content` y descarte transaccional sin prefijo publicado.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.8.3`; debe terminar con código 0 y escribir `evidence/f2/F2.2.8.3/result.json` con vectores de ambos selectores y `status = "passed"`.

- ☐ F2.2.8.4. Implementar las distribuciones continuas semiabiertas.
  - Depende de: F2.2.2, F2.2.6.1 y F2.2.8.2.
  - Procede de: DTC-043, DTC-049, DTC-068, DTC-072 y DTC-075; RAF2-041/043 resueltos.
  - Cambia: productores únicos `linear-uniform-v1` y `log-uniform-v1`; ambos consumen `[min,max)`, permiten `max=1` sin devolverlo, impiden 0 continuo y la variante logarítmica usa exclusivamente `libm = =0.2.16` aprobado en F2.2.2.
  - Se comprueba con: vectores de cero, subnormales, intervalo sin `f64` positivo interior, `next_down(1)`, `max=1`, redondeo inferior/superior y extremos inválidos; oráculo independiente confirma que solo `fixed` produce 0/1 y que F2.2.2 no implementa ni prueba otra copia logarítmica.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.8.4`; debe terminar con código 0 y escribir `evidence/f2/F2.2.8.4/result.json` con fronteras, versión de `libm` y `status = "passed"`.

- ☐ F2.2.8.5. Implementar los consumidores de composición genotípica y fase.
  - Depende de: F2.2.8.2–F2.2.8.3.
  - Procede de: DTC-052, DTC-071, DTC-072 y DTC-075; REQ-IND-62/63.
  - Cambia: productores únicos `hwe-composition-v1`, `f-diploid-genotype-v1` y `phase-permutation-v1`, sin adelantar su integración biológica F2.3.
  - Se comprueba con: enumeración HWE/F ajena, ploidías 1/2/10, distribuciones degeneradas, orden de composición/fase y consumo exacto independiente entre genotipo y fase.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.8.5`; debe terminar con código 0 y escribir `evidence/f2/F2.2.8.5/result.json` con vectores y `status = "passed"`.

- ☐ F2.2.8.6. Verificar conjuntamente versiones, consumos y aislamiento de la familia RNG.
  - Depende de: F2.2.8.1–F2.2.8.5.
  - Procede de: DTC-038, DTC-071, DTC-072 y DTC-075; REQ-IND-63 y REQ-VAL-55.
  - Cambia: oráculo de integración ajeno a los productores y expediente de propiedad única que enlaza cada ámbito con exactamente un algoritmo, módulo, vector y consumidor.
  - Se comprueba con: mutación aislada de componentes/algoritmos, reordenación, adición de población/individuo/clase y búsqueda de productores duplicados; Linux/WSL2 producen bytes y consumos idénticos.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.8.6`; debe terminar con código 0 y escribir `evidence/f2/F2.2.8.6/result.json` con matriz ámbito→productor→consumidor y `status = "passed"`.

- ☐ F2.2.9. Enumerar bases elegibles antes de determinar o seleccionar loci.
  - Depende de: F2.1.5–F2.1.6 y F2.2.7.
  - Procede de: REQ-CONF-21, REQ-CONF-28, REQ-IND-13, REQ-IND-14, REQ-IND-19 y REQ-DIR-03.
  - Cambia: servicio que cuenta posiciones `A/C/G/T` elegibles tras inclusiones, exclusiones, posiciones protegidas/ocupadas y pesos cero; excluye `N` y no sortea todavía.
  - Se comprueba con: enumeración ajena al productor en secuencias cortas, contigs/regiones literales, espacio vacío y denominador de bases elegibles por cromosoma.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.9`; debe terminar con código 0 y escribir `evidence/f2/F2.2.9/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.10. Implementar `uniform` y `weighted_regions` sin solapamientos.
  - Depende de: F2.2.3, F2.2.7 y F2.2.9.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-CONF-129, REQ-IND-19 y REQ-VAL-55.
  - Cambia: perfil incorporado único `uniform` con peso 1; tabla propia `1-based-closed` con pesos finitos no negativos, `outside_weight` explícito, regiones no solapadas y probabilidad por base proporcional al peso.
  - Se comprueba con: tabla manual por longitud×peso, casos `outside_weight` 0/1, peso 0/0.5/1/2, frontera contigua, solapamiento rechazado y requisito de al menos una base con peso positivo.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.10`; debe terminar con código 0 y escribir `evidence/f2/F2.2.10/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.11. Determinar el número de loci candidatos sin confundirlo con loci activos.
  - Depende de: F2.2.2, F2.2.7, F2.2.8.1 y F2.2.9–F2.2.10.
  - Procede de: REQ-CONF-129, REQ-IND-17, REQ-IND-18 y REQ-IND-61; DTC-022, DTC-052, DTC-071 y DTC-075; RAF2-044 resuelto.
  - Cambia: `exact_loci` fija candidatos sin consumir `matrix_candidate_count`; `mean_loci` y `density_per_mb × bases_elegibles/1_000_000` fijan una media y llaman una vez a `poisson-count-v1` con ese ámbito. Los tres modos son excluyentes; se registran algoritmo, media, ámbito/subsemilla, consumo y recuento aunque después se eliminen ALT globalmente ausentes.
  - Se comprueba con: exactos con consumo cero/ausencia de subsemilla sorteada, vectores Poisson para media/densidad, denominador manual, tasa cero, espacio vacío/insuficiente, aislamiento frente a `matrix_locus_selection`/ALT/frecuencias y manifiesto que reconstruye el recuento sin ambigüedad.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.11`; debe terminar con código 0 y escribir `evidence/f2/F2.2.11/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.12. Validar espectros alélicos por intervalos no solapados.
  - Depende de: F2.2.3, F2.2.6.1–F2.2.7 y F2.2.8.4.
  - Procede de: DTC-054, DTC-068, DTC-072 y DTC-075; REQ-CONF-129, REQ-IND-61 y REQ-VAL-55; RAF2-041 resuelto.
  - Cambia: lector `f64` de decimal/científica para intervalos siempre `[min_af,max_af)`, incluido `max_af=1`, y filas `fixed`; admite huecos, rechaza solapamiento, porcentajes, NaN/infinito y pesos inválidos; solo `fixed` genera exactamente 0/1.
  - Se comprueba con: corpus de fronteras 0/1, subnormales, `next_down(1)`, intervalos contiguos/solapados/con huecos, pesos no normalizados, todo cero, `1e-3`, porcentaje, NaN/infinito y oráculo que prueba que el máximo 1 cubre pero nunca se sortea.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.12`; debe terminar con código 0 y escribir `evidence/f2/F2.2.12/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.13. Implementar los repartos `stochastic` y `exact` por población.
  - Depende de: F2.2.8.1–F2.2.8.3 y F2.2.11–F2.2.12.
  - Procede de: DTC-054, DTC-070, DTC-071 y DTC-072 (resuelven DTPD-118 y RAF2-028/029/033); REQ-CONF-129, REQ-IND-61, REQ-IND-63 y REQ-VAL-55.
  - Cambia: sorteo `categorical-f64-v1` de fila en modo estocástico. En `exact`, calcula cuotas por mayores restos y cada población aplica `fisher-yates-v1` con `matrix_exact_assignment` a todos los `variant_id`; no usa coordenadas ni orden canónico para asociar filas.
  - Se comprueba con: 11 loci y pesos 60/25/10/5→7/3/1/0, suma exacta por población, empate literal, estabilidad al reordenar entrada, ausencia de asociación con coordenada, aislamiento entre poblaciones, cuota cero sin consumo RNG, separación respecto al flujo de valor y prueba estadística prefijada del modo estocástico.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.13`; debe terminar con código 0 y escribir `evidence/f2/F2.2.13/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.14. Aprobar e incorporar los dos perfiles humanos compactos aproximados.
  - Depende de: F2.2.14.1–F2.2.14.4.
  - Procede de: DTC-054 y DTC-068 (resuelven DTPD-118/121); REQ-CONF-129, REQ-IND-61 y ACE-116; RAF2-018/036 resueltos.
  - Cambia: publica conjuntamente fuente fijada, perfil SFS, perfil SNV12 y reproducción independiente, manteniendo cada responsabilidad en un expediente revisable.
  - Se comprueba con: los cuatro expedientes hijos aprobados, digests canónicos de ambos recursos y prueba de selección nominal; ningún builtin se activa por defecto.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.14`; debe terminar con código 0 y escribir `evidence/f2/F2.2.14/result.json` con `status = "passed"` y enlaces a F2.2.14.1–F2.2.14.4.

- ☐ F2.2.14.1. Fijar e inventariar la fuente pequeña Harris 2017.
  - Depende de: F2.2.1 y DTC-068.
  - Procede de: DTC-068; atribución y licencia de ACE-116.
  - Cambia: copia de prueba versionada de `elife-24284-fig1-data1-v2.txt`, URL/DOI, Nature/eLife, licencia CC BY, tamaño 3.644 bytes y SHA-256 `970cc71578fc008ef4e0d75c7c460ab874073da7c5b5226c526b387627c97a7f`.
  - Se comprueba con: tamaño/hash antes de parsear, licencia/atribución presentes, cinco filas continentales y 96 contextos exactos; una diferencia bloquea sin descargar ni sustituir la fuente.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.14.1`; debe terminar con código 0 y escribir `evidence/f2/F2.2.14.1/result.json` con inventario y `status = "passed"`.

- ☐ F2.2.14.2. Generar el perfil SFS global aproximado versionado.
  - Depende de: F2.2.3, F2.2.8.4 y F2.2.12–F2.2.13.
  - Procede de: DTC-054 y DTC-068; `human_1000g_global_sfs_approx_v1`.
  - Cambia: generador offline pequeño para pesos 64:12:8, intervalos `[1/5008,0.005)`, `[0.005,0.05)` y `[0.05,1)`, muestreos log/log/lineal, literales, orden y digest del recurso.
  - Se comprueba con: límites y pesos exactos, normalización DTC-049, máximo 1 nunca sorteado, ausencia de masa fija en cero y de frecuencias positivas inferiores a `1/5008`, metadatos que lo identifican como aproximación global no universal.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.14.2`; debe terminar con código 0 y escribir `evidence/f2/F2.2.14.2/result.json` con recurso candidato y SHA-256.

- ☐ F2.2.14.3. Generar el perfil SNV12 Harris 2017 aproximado versionado.
  - Depende de: F2.2.6.1, F2.2.14.1 y F2.2.3.
  - Procede de: DTC-049 y DTC-068; `human_harris2017_snv12_approx_v1`.
  - Cambia: parser y receta que consumen el único normalizador F2 de DTC-049 para normalizar AFR/AMR/EUR/EAS/SAS por fila, promedian con igual peso, agregan 96 contextos a seis clases por hebra, reflejan complementos inversos con igual peso y normalizan 12 cambios dirigidos.
  - Se comprueba con: `cargo +1.97.1 run --locked --bin gps-build-human-snv-profiles -- --harris-source tests/fixtures/profiles/elife-24284-fig1-data1-v2.txt --output-dir target/f2-profile-candidate`; recuentos, sumas intermedias, complementos, orden, atribución y digest son exactos.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.14.3`; debe terminar con código 0 y escribir `evidence/f2/F2.2.14.3/result.json` con todas las sumas intermedias y SHA-256 final.

- ☐ F2.2.14.4. Reproducir independientemente y autorizar la incorporación de ambos perfiles.
  - Depende de: F2.2.14.1–F2.2.14.3.
  - Procede de: DTC-068; puerta reproducible de RAF2-018.
  - Cambia: oráculo ajeno al productor que recalcula recursos/digests y puerta que los incorpora solo tras coincidencia exacta; prohíbe VCF Phase 3, `bcftools`, perfiles inventados por superpoblación y carga predeterminada.
  - Se comprueba con: segundo cálculo de implementación independiente, mutaciones de fuente/literal/orden/licencia y selección nominal; perfiles suministrados por población siguen siendo válidos.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.14.4`; debe terminar con código 0 y escribir `evidence/f2/F2.2.14.4/result.json` con comparación byte a byte, hashes y aprobación explícita.

- ☐ F2.2.15. Resolver el espectro dirigido de las 12 sustituciones.
  - Depende de: F2.2.3 y F2.2.7; solo el builtin Harris 2017 depende además de DTC-068/F2.2.14.
  - Procede de: DTC-054 y DTC-068; REQ-CONF-129 y REQ-IND-61.
  - Cambia: `uniform_snv12_v1` predeterminado, tabla propia, atajo `ti_tv_ratio=R` convertido a peso `2R` para la transición y 1 para cada transversión, y perfil explícito `human_harris2017_snv12_approx_v1`; la configuración efectiva conserva siempre los 12 pesos.
  - Se comprueba con: las tres ALT por cada REF, R=0/1/2, rechazo de R negativo/infinito y razón esperada independiente; para Harris se normaliza cada fila AFR/AMR/EUR/EAS/SAS, se promedia a igual peso, se colapsan los 96 contextos a seis clases y cada una se refleja con mitad de peso a su complemento inverso. El oráculo no comparte productor y el builtin nunca se carga por defecto ni se presenta como tasa universal.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.15`; debe terminar con código 0 y escribir `evidence/f2/F2.2.15/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.16. Crear IDs canónicos y orden independiente del orden de entrada.
  - Depende de: F2.2.3–F2.2.4 y servicios de identidad F1.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-IND-60, REQ-IND-61 y REQ-VAL-55.
  - Cambia: `site_id` por referencia/contig/posición/REF y `variant_id` por referencia/contig/posición/REF/ALT; conserva `source_id` como alias y ordena por catálogo FASTA, posición, REF, ALT y población.
  - Se comprueba con: permutación de filas/poblaciones, mismo SNV en referencias distintas, sitio multialélico y comparación literal de IDs y bytes ordenados.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.16`; debe terminar con código 0 y escribir `evidence/f2/F2.2.16/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.17. Preparar una matriz suministrada sin modificar el original.
  - Depende de: F2.2.3–F2.2.4, F2.2.9 y F2.2.16.
  - Procede de: DTC-054 y DTC-060; REQ-CONF-128, REQ-IND-60, REQ-VAL-55 y REQ-MAN-91.
  - Cambia: orden `prepare-snv-matrix --reference-manifest <manifest.json> --reference-id <id> --input <TSV> --output <destino>` que valida el manifiesto F1, resuelve exactamente su entrada de referencia, verifica ruta/tamaño/SHA-256 y después calcula IDs/orden, comprueba contigs/coordenadas/REF y publica la pareja nueva.
  - Se comprueba con: ID correcto, ID distinto con mismos bytes, mismo ID con bytes distintos, referencia derivada válida/incoherente y manifiesto corrupto; fixture mínimo original inalterado, fallo sin pareja parcial y equivalencia Marimo «Vincular y verificar referencia».
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.17`; debe terminar con código 0 y escribir `evidence/f2/F2.2.17/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.18. Validar filas completas y poblaciones de una matriz suministrada.
  - Depende de: F2.2.3–F2.2.4 y F2.2.16–F2.2.17.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-CONF-128, REQ-IND-60 y REQ-VAL-55.
  - Cambia: exige una fila por `ALT × población`, `alt_af` explícito en `[0,1]`, población declarada y claves únicas; rechaza vacío/NA, MAF/DAF, AC/AN, inferencias y alias de población, pero admite y valida poblaciones no seleccionadas.
  - Se comprueba con: corpus que omite una combinación, usa cero explícito, añade población válida no usada, solicita una inexistente e intenta cada semántica de frecuencia prohibida.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.18`; debe terminar con código 0 y escribir `evidence/f2/F2.2.18/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.19. Validar loci multialélicos suministrados con suma estable y tolerancia aprobada.
  - Depende de: F2.2.6.1 y F2.2.18.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-IND-60 y REQ-VAL-55.
  - Cambia: agrupa por sitio/población, calcula suma ALT con Neumaier en orden canónico y `ref_af=1-suma`; acepta suma ≤1, normaliza solo un exceso `delta <= tolerancia_efectiva` de DTC-049 y registra tolerancia, suma/delta y valores antes/después; rechaza exceso mayor sin recortar valores individuales.
  - Se comprueba con: tolerancia efectiva 0, `1e-9` y un valor sustituido válido; para cada una prueba `delta == tolerancia` y el siguiente `f64` mayor, sumas 0.5/1/1.01, orden permutado y lector independiente de frecuencias/normalización/diagnóstico.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.19`; debe terminar con código 0 y escribir `evidence/f2/F2.2.19/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.20. Normalizar alelos globalmente ausentes sin reponer candidatos.
  - Depende de: F2.2.18–F2.2.19.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-IND-60, REQ-IND-61, REQ-VAL-55 y REQ-MAN-91.
  - Cambia: elimina cada ALT con `alt_af=0` en todas las poblaciones y el sitio si pierde todos sus ALT; conserva `alt_af=1`, cuenta candidatos/ALT descartados/loci activos por separado y nunca selecciona sustitutos.
  - Se comprueba con: bialélico global 0/1, multialélico con un ALT global 0 y otro activo, matriz suministrada/generada y recuentos antes/después calculados fuera del normalizador.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.20`; debe terminar con código 0 y escribir `evidence/f2/F2.2.20/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.21. Implementar autorización explícita para una matriz activa vacía.
  - Depende de: F2.2.20 y diagnósticos F1.
  - Procede de: DTC-054 y DTC-073 (resuelven DTPD-118 y RAF2-034); REQ-VAL-55 y REQ-MAN-91.
  - Cambia: incorpora conjuntamente parser, validación y semántica de `--allow-empty-snv-matrix` y solo entonces la añade a `available_authorizations`. El primer `validate` avisa; confirmarla crea otro digest y obliga a validar de nuevo. Solo `run` con los mismos documentos, contrato y autorizaciones continúa; queda en configuración efectiva/manifiesto y no relaja errores.
  - Se comprueba con: antes del paso está catalogada pero no disponible/visible; después, capabilities/parser/validate/run coinciden. Cubre matriz vacía inicial/tras normalización, confirmación→revalidación, edición/flag que bloquea, SNP explícitos y negativos de digest, REF, frecuencia, opción desconocida o disponibilidad sin semántica.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.21`; debe terminar con código 0 y escribir `evidence/f2/F2.2.21/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.22. Seleccionar las posiciones candidatas sin reemplazo.
  - Depende de: F2.2.8.1–F2.2.8.3 y F2.2.9–F2.2.11.
  - Procede de: DTC-054 y DTC-071/072 (resuelven DTPD-118 y RAF2-028/033); REQ-CONF-129, REQ-IND-03, REQ-IND-19 y REQ-IND-61.
  - Cambia: `weighted-without-replacement-v1` sobre bases elegibles en orden canónico obtiene exactamente el recuento candidato, sin duplicados, y falla en preflight si faltan posiciones; cero candidatos no consume y no promete prefijo al cambiar N.
  - Se comprueba con: referencia exhaustiva `A/C/G/T/N`, enumeración literal de pesos/acumulados/elecciones, N=0/1/máximo/exceso, forzado de rechazo entero, semilla fija, consumo exacto y comparación con un enumerador independiente.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.22`; debe terminar con código 0 y escribir `evidence/f2/F2.2.22/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.23. Asignar un único ALT generado por locus candidato.
  - Depende de: F2.2.8.1–F2.2.8.3, F2.2.15–F2.2.16 y F2.2.22.
  - Procede de: DTC-054 y DTC-071/072 (resuelven DTPD-118 y RAF2-028/033); REQ-CONF-129 y REQ-IND-61.
  - Cambia: generador bialélico que lee REF del FASTA y aplica `categorical-f64-v1` a sus tres ALT mediante `matrix_alternate_allele`, identificado antes de ALT; deriva `variant_id` solo después. No genera multialelia ni contexto trinucleotídico.
  - Se comprueba con: cada REF A/C/G/T, pesos degenerados con consumo cero, perfiles uniforme/propio/builtin, vectores de scope sin ALT y prueba de que jamás produce REF, `N`, segundo ALT o cambio fuera del perfil.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.23`; debe terminar con código 0 y escribir `evidence/f2/F2.2.23/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.24. Asignar frecuencias independientes a todas las poblaciones.
  - Depende de: F2.2.8.1, F2.2.8.3–F2.2.8.4, F2.2.12–F2.2.13, F2.2.16 y F2.2.23.
  - Procede de: DTC-054, DTC-070, DTC-071 y DTC-072 (resuelven DTPD-118 y RAF2-028/029/033); REQ-CONF-129, REQ-IND-61, REQ-IND-63 y REQ-VAL-55.
  - Cambia: en modo estocástico aplica `categorical-f64-v1`; en `exact` consume la fila asociada por F2.2.13. Después usa `linear-uniform-v1` o `log-uniform-v1` en el ámbito separado por población×`variant_id`. Admite 0/1 solo mediante `fixed`, conserva versiones/subsemillas y no consume RNG para `fixed` ni cuota cero.
  - Se comprueba con: dos poblaciones/perfiles, frecuencias `fixed` 0/1 e intervalos `[min,max)` que nunca producen 0/1 aunque `max=1`, reordenación, consumo cero, vectores lineal/log de F2.2.8.4, separación asignación/valor e independencia sin exigir correlación.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.24`; debe terminar con código 0 y escribir `evidence/f2/F2.2.24/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.25. Publicar la matriz generada mediante el mismo contrato canónico.
  - Depende de: F2.2.4, F2.2.16 y F2.2.20–F2.2.24.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-CONF-129, REQ-IND-61 y REQ-MAN-91.
  - Cambia: escritor determinista de TSV+JSON con origen generado, parámetros/perfiles efectivos, candidatos/descartes/activos, algoritmo de recuento, `matrix_candidate_count`/subsemilla/consumo cuando corresponda, fuente/versión y digests, usando publicación segura. En `exact_loci` registra modo y recuento sin fingir consumo ni subsemilla sorteada.
  - Se comprueba con: dos generaciones de igual referencia/configuración/semilla con mismos bytes, rutas exacta/Poisson, hashes externos, fallo inyectado sin pareja parcial y lector independiente que reconstruye todos los metadatos y distingue ausencia normativa de ámbito en `exact_loci`.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.25`; debe terminar con código 0 y escribir `evidence/f2/F2.2.25/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.26. Ofrecer generación independiente y reutilización sin reescribir la matriz.
  - Depende de: F2.2.17 y F2.2.25.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-CONF-128, REQ-CONF-129, REQ-IND-60 y REQ-IND-61.
  - Cambia: orden `generate-snv-matrix config.toml --output ...` y consumo posterior de la pareja canónica; el recurso conserva bytes/digest/origen y el manifiesto consumidor registra su reutilización.
  - Se comprueba con: generación autónoma e integrada equivalentes byte a byte, segunda ejecución con individuos distintos que no modifica `mtime`/hash del recurso y rechazo de pareja parcial o mezclada.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.26`; debe terminar con código 0 y escribir `evidence/f2/F2.2.26/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.27. Validar y aplicar SNP explícitos como restricciones trazables, no como otra fuente aleatoria inicial.
  - Depende de: F2.1.8–F2.1.9 y F2.2.6–F2.2.7.
  - Procede de: REQ-IND-05, REQ-IND-11 y REQ-REC-24; DTC-020, DTC-032 y DTC-052.
  - Cambia: evento SNP explícito con referencia esperada, ALT `A/C/G/T`, copia y procedencia; prepara su contraste posterior con la matriz y rechaza origen/destino `N`.
  - Se comprueba con: transición/transversión, desacuerdo, copia ausente, borde y `N`, más caso negativo que no acepta una frecuencia poblacional como valor literal individual.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.27`; debe terminar con código 0 y escribir `evidence/f2/F2.2.27/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.28. Publicar ámbitos RNG propios de inserciones y deleciones aleatorias.
  - Depende de: F2.2.2, F2.2.7, F2.2.8.1 y F2.2.8.3.
  - Procede de: DTC-037, DTC-038, DTC-059, DTC-071, DTC-072 y DTC-075; RAF2-016/028/033/042 resueltos.
  - Cambia: consume `indel_count`, `indel_position` e `insertion_content` v4 con componentes exactos; enlaza `indel_position` exclusivamente a `dynamic-indel-selection-v1`. Deleción no tiene contenido; clases, matriz, recuento, posición y bases no comparten estado.
  - Se comprueba con: vectores BLAKE3/componentes/orden, `event_index=0/1`, `global`, región tipada, tasa cero sin avance y cambios aislados de matriz/clase/región/individuo/secuencia.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.28`; debe terminar con código 0 y escribir `evidence/f2/F2.2.28/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.29. Resolver medias y regiones de indels sin reutilizar el selector de matriz.
  - Depende de: F2.2.5–F2.2.7 y F2.2.28.
  - Procede de: REQ-CONF-20 a REQ-CONF-22, REQ-IND-10, REQ-IND-13 a REQ-IND-15; DTC-033 y DTC-059.
  - Cambia: parámetros independientes de inserción/deleción para media o densidad, distribución uniforme y regiones efectivas por individuo; crea ámbitos espaciales disjuntos según DTC-033 y calcula `media = densidad × oportunidades_validas`, sin `weighted_regions` de matriz.
  - Se comprueba con: valores global/región/individuo, densidad cero, regiones solapadas/incomparables, denominador manual y prueba de que modificar `outside_weight` o un perfil de matriz no cambia ninguna media de indel.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.29`; debe terminar con código 0 y escribir `evidence/f2/F2.2.29/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.30. Enumerar oportunidades de indels sobre el estado vigente.
  - Depende de: F2.1.5–F2.1.9 y F2.2.29.
  - Procede de: REQ-IND-13/14/19/21/40/41 y DTC-020/DTC-032/DTC-075; RAF2-042 resuelto.
  - Cambia: enumeradores separados de anclas de inserción y deleciones flanqueadas para la longitud exacta efectiva; excluyen destinos inválidos, respetan `N` según clase y recalculan después de cada evento aceptado sin reintentar candidatos descartados.
  - Se comprueba con: enumeración exhaustiva ajena sobre cadenas cortas, bordes, `N`, posiciones ocupadas y estado tras inserción/deleción previa; espacio vacío no consume selección.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.30`; debe terminar con código 0 y escribir `evidence/f2/F2.2.30/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.31. Sortear una vez el recuento Poisson local y seleccionar secuencialmente.
  - Depende de: DTC-067, F2.2.2, F2.2.8.3 y F2.2.28–F2.2.30.
  - Procede de: DTC-020, DTC-022, DTC-037, DTC-067, DTC-071, DTC-072 y DTC-075; REQ-CONF-23 a REQ-CONF-25 y REQ-IND-16 a REQ-IND-18; RAF2-042 resuelto.
  - Cambia: cada clase/individuo/secuencia/ámbito llama una vez a `poisson-count-v1` con `indel_count` y entrega `sampled_count` a `dynamic-indel-selection-v1` en `indel_position`. Por ordinal reenumera/pondera/aplica sobre el estado vigente; registra versión, lambda, subsemillas, oportunidades iniciales/vigentes, sorteado, aceptado y activo antes del rollback. Agotamiento descarta la transacción, registra `published_count=0` y no hace nuevo Poisson, recorte, reintento ni publicación parcial.
  - Se comprueba con: vectores DTC-067/072, lambda cero sin consumo, límite/next-up, cadenas independientes donde un evento crea o elimina oportunidades, un único candidato degenerado, espacio vacío inicial/tardío, orden canónico, flujo continuo y fallo inyectado que deja cero eventos publicados y conserva diagnóstico completo.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.31`; debe terminar con código 0 y escribir `evidence/f2/F2.2.31/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.32. Generar contenido de inserción aislado y verificar independencia completa.
  - Depende de: F2.2.8.3 y F2.2.28–F2.2.31.
  - Procede de: DTC-032, DTC-037, DTC-059 y DTC-071/072; REQ-IND-40 y RAF2-016/028/033 resueltos.
  - Cambia: cada inserción aplica un categórico uniforme A/C/G/T por base en `insertion_content` con su `event_index`; consume exactamente la longitud efectiva. Deleción y longitud no consumen ese ámbito.
  - Se comprueba con: RNG instrumentado y vectores de bases/consumos para longitudes 1/2/`umbral - 1`, ordinales 0/1, ausencia de `N` y aislamiento frente a matriz/individuo/región/clase.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.32`; debe terminar con código 0 y escribir `evidence/f2/F2.2.32/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.33. Validar y aplicar inserciones pequeñas explícitas y aleatorias.
  - Depende de: F2.1.5, F2.1.8–F2.1.9, F2.2.2, F2.2.6–F2.2.7 y F2.2.28–F2.2.32.
  - Procede de: REQ-IND-05, REQ-IND-40 y REQ-REC-24; DTC-032 y DTC-059.
  - Cambia: inserción interna anclada en la base precedente; la explícita conserva `N` literal y la aleatoria genera solo `A/C/G/T` con longitud exacta efectiva.
  - Se comprueba con: oráculo manual, terminal imposible, `N` explícita, longitudes 1/2/`umbral - 1` y barrido aleatorio sin `N`.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.33`; debe terminar con código 0 y escribir `evidence/f2/F2.2.33/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.34. Validar y aplicar deleciones pequeñas solo con bases flanqueantes.
  - Depende de: F2.1.5, F2.1.8–F2.1.9, F2.2.2, F2.2.6–F2.2.7 y F2.2.28–F2.2.31.
  - Procede de: REQ-IND-05, REQ-IND-41 y REQ-REC-24; DTC-059.
  - Cambia: deleción con intervalo/REF explícitos o longitud exacta aleatoria; rechaza cero, desacuerdo, solapamiento inválido y cualquier deleción terminal.
  - Se comprueba con: cadenas cortas, longitudes 1/2/`umbral - 1`, intentos sobre primera/última/toda la secuencia y estado intacto tras error.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.34`; debe terminar con código 0 y escribir `evidence/f2/F2.2.34/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.35. Componer eventos pequeños sobre el estado vigente con conflictos explícitos.
  - Depende de: F2.2.27–F2.2.34.
  - Procede de: DTC-020 y DTC-031; REQ-REC-24 y REQ-IND-21.
  - Cambia: aplicación secuencial de explícitos y demás clases aleatorias en orden normativo, con coordenada de referencia, alelo vigente y diagnóstico de objetivo inexistente; la matriz todavía no se muestrea.
  - Se comprueba con: SNP→deleción, deleción→SNP, inserciones vecinas y repetidas con resultado manual y prueba de ausencia de edición parcial.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.35`; debe terminar con código 0 y escribir `evidence/f2/F2.2.35/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.36. Cerrar diagnósticos y preflight de las entradas F2.2.
  - Depende de: F2.2.6–F2.2.35.
  - Procede de: DTC-044; REQ-CONF-25, REQ-CONF-28, REQ-IND-40, REQ-IND-41 y REQ-VAL-55.
  - Cambia: catálogo estable y agregado para formato/digest, referencia, población, coordenada, alelo, espectro/intervalo, región/peso, recuentos, `N`, *indel*, matriz vacía y conflicto.
  - Se comprueba con: tabla diagnóstico→caso, texto/JSON equivalentes, orden determinista, flags limitados y configuración inválida que no crea artefactos finales.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.36`; debe terminar con código 0 y escribir `evidence/f2/F2.2.36/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.37. Construir un validador independiente de los formatos, digests y normalización.
  - Depende de: F2.2.3–F2.2.4; no depende de F2.2.17–F2.2.26 para calcular el esperado.
  - Procede de: REQ-VAL-55 y ACE-116; DTC-047 y DTC-054 (resuelve DTPD-118).
  - Cambia: corpus escrito a mano y lector separado para TSV+JSON, referencia, filas completas, IDs, orden, multialelia, sumas Neumaier, tolerancia y eliminación de ALT cero.
  - Se comprueba con: mutación deliberada de cada campo/digest/invariante y reconstrucción de resultados literales sin llamar al parser, normalizador ni escritor de producción.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.37`; debe terminar con código 0 y escribir `evidence/f2/F2.2.37/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.38. Construir oráculos independientes para selección, espectros y perfiles.
  - Depende de: F2.2.8–F2.2.15; no depende de F2.2.22–F2.2.25 para producir el esperado.
  - Procede de: REQ-VAL-55 y ACE-107/116; DTC-047, DTC-054, DTC-070 y DTC-075; DTPD-118 y RAF2-029/041/042 resueltos.
  - Cambia: casos sintéticos y cálculos separados de elegibilidad/selección SNV estática, selección indel dinámica, Poisson, cuotas, permutación `exact`, intervalos semiabiertos, ti/tv, SFS 64:12:8 y Harris 96→6→12.
  - Se comprueba con: vectores y semillas fijas que cubren cantidades, reordenación, aislamiento, consumo cero, extremos 0/1, creación/pérdida dinámica de oportunidades, agotamiento transaccional y separación de flujos; pruebas estadísticas y mutaciones del productor se detectan sin reutilizar sus funciones.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.38`; debe terminar con código 0 y escribir `evidence/f2/F2.2.38/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.39. Integrar fuente, matriz efectiva, digest y procedencia sin crear individuos.
  - Depende de: F2.2.17–F2.2.26 y F2.2.36–F2.2.38.
  - Procede de: REQ-CONF-128, REQ-CONF-129, REQ-IND-06, REQ-IND-60 y REQ-IND-61; DTC-040 y DTC-052.
  - Cambia: resultado F2.2 inmutable con referencia, poblaciones disponibles/seleccionadas, matriz preparada/generada/reutilizada, perfiles, candidatos/descartes/activos, algoritmos, ámbitos/subsemillas y consumos aplicables —incluido `matrix_candidate_count` o su ausencia normativa—, versiones y SHA-256 para F2.3.
  - Se comprueba con: las tres rutas, referencias cruzadas, `sha256sum`, igualdad con semilla fija, poblaciones adicionales válidas y ausencia demostrada de genotipos o SNP individuales sorteados.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.39`; debe terminar con código 0 y escribir `evidence/f2/F2.2.39/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.2.40. Completar evidencia y trazabilidad de F2.2 con oráculos ajenos al productor.
  - Depende de: F2.2.8–F2.2.39 y puertas de calidad verdes.
  - Procede de: los 22 REQ de cierre y seis de evidencia parcial asignados a F2.2 en §15/§15.1; ACE-01, ACE-09, ACE-10, ACE-21, ACE-22, ACE-28, ACE-85, ACE-95, ACE-103, ACE-104, ACE-109, ACE-110, ACE-111 y ACE-116.
  - Cambia: pruebas unitarias, de propiedades, estadísticas e integración y matriz de evidencia para los 22 cierres exclusivos, REQ-CONF-08/21 parciales sobre matriz e indels y REQ-CONF-106/107, REQ-VAL-40 y REQ-MAN-81 parciales sobre el normalizador único.
  - Se comprueba con: §4.3, recuento automático de 22 cierres y seis expedientes parciales —REQ-CONF-08/21 hacia F5.2/F7.2 y los cuatro de DTC-049 hacia F6.2—, revisión de F2.2.37/F2.2.38 y repetición exacta; ningún parcial aparece como cierre total ni atribuye creación de individuos a F2.2.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.2.40`; debe terminar con código 0 y escribir `evidence/f2/F2.2.40/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

> **PAUSA-F2.2 — cierre obligatorio de F2.2.** Confirmar el consumo sin reinterpretación de
> DTC-054 y el cierre de DTPD-118, la pareja TSV+JSON y sus digests, perfiles solo explícitos,
> candidatos/descartes/activos, poblaciones independientes, validadores,
> Poisson limitado a recuentos, espacios elegibles y casos con `N`. No iniciar
> F2.3 si un oráculo se deriva del productor, una densidad SNP individual sigue
> activa, un flag relaja otra validación o un rechazo deja cambios.

## 8. F2.3 — Individuos directos desde una población de origen

**Resultado del hito:** la vía directa muestrea una matriz común mediante HWE o
`F` diploide, conserva fase y valores esperados/realizados, integra los demás
cambios y reconstruye cada haplotipo sin simular dinámica poblacional.

- ☐ F2.3.1. Extender la autoría HWE/F y trasladar la distribución explícita a F6.2.
  - Depende de: PAUSA-F2.2 aprobada y F2.2.1–F2.2.6.
  - Procede de: REQ-CONF-20, REQ-CONF-22, REQ-CONF-130, REQ-DIR-01 y REQ-DIR-02; DTC-034 y DTC-052.
  - Cambia: F2 admite HWE para P=1..10 y F únicamente en diploides; P!=2 rechaza F y usa HWE. No añade esquema, fichero ni controles Marimo para distribución genotípica explícita, cuyo cierre completo queda en F6.2 con REQ-EVO-35/36/49.
  - Se comprueba con: HWE en ploidías 1/2/10, F=0/intermedio/1 solo P=2, rechazo de F en P=1/3/10 y búsqueda negativa de distribución explícita en esquema/CLI/Marimo F2; expediente `REQ-CONF-130.status=partial` con destino F6.2.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.3.1`; debe terminar con código 0 y escribir `evidence/f2/F2.3.1/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.3.2. Ampliar el catálogo F2 con individuo y población de origen.
  - Depende de: F2.3.1 y contratos F1 de identidad/configuración.
  - Procede de: DTC-033, DTC-039 y DTC-044; REQ-CONF-20, REQ-CONF-22, REQ-CONF-38, REQ-CONF-130 y evidencia compartida de REQ-IND-02/03/25 y REQ-MAN-10.
  - Cambia: amplía el catálogo/dominio creado en F2.2.7 con individuo y población de origen, y configura/valida/resuelve sus IDs tipados y procedencia completa. Cuando falta un ID de individuo, exige `founder_key` único y estable dentro de la población; deriva el ID con DTC-039 desde raíz, población y `founder_key`, comprueba colisión y lo reserva en el registro de ejecución. Un ID no se reutiliza dentro del recorrido disponible y el orden físico nunca aporta una causa. No introduce subpoblaciones, generaciones ni jerarquías futuras.
  - Se comprueba con: ID aportado, omitido con `founder_key` y derivado, dos anónimos equivalentes sin claves distintas, `founder_key` duplicado, duplicado/colisión de ID, estabilidad ante reordenación y repetición, ausencia de reutilización, IDs desconocidos, orden, sustitución individual, HWE omitido/predeterminado y aislamiento de valores efectivos.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.3.2`; debe terminar con código 0 y escribir `evidence/f2/F2.3.2/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.3.3. Versionar subsemillas de genotipo y fase separadas de la generación de matriz.
  - Depende de: F2.2.8.1, F2.2.8.5–F2.2.8.6 y F2.3.1–F2.3.2.
  - Procede de: REQ-IND-63 y REQ-VAL-55; DTC-011, DTC-037, DTC-038, DTC-052 y DTC-071/072.
  - Cambia: consume `initial_genotype` e `initial_phase` del catálogo v4 con escenario, referencia, población, individuo y variante en orden normativo; enlaza `hwe-composition-v1`/`f-diploid-genotype-v1` y `phase-permutation-v1` sin RNG mutable compartido.
  - Se comprueba con: vectores literales de ambos ámbitos/algoritmos, componentes omitidos/reordenados rechazados, mismo resultado con distinto orden físico y estabilidad al añadir otro ID.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.3.3`; debe terminar con código 0 y escribir `evidence/f2/F2.3.3/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.3.4. Construir el estado inicial faseado de cada individuo desde la referencia.
  - Depende de: F2.1.6–F2.1.9 y F2.3.2.
  - Procede de: DTC-003 y DTC-006; REQ-DIR-01, REQ-DIR-02, REQ-IND-28 y REQ-IND-29; evidencia compartida de REQ-IND-02/03.
  - Cambia: copias intactas ordenadas por secuencia/ploidía y enlaces a referencia, población, individuo y procedencia antes de asignar alelos.
  - Se comprueba con: copias 1–10 iguales byte a byte a la referencia, cardinalidad/IDs exactos y ausencia de escenarios o poblaciones simuladas.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.3.4`; debe terminar con código 0 y escribir `evidence/f2/F2.3.4/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.3.5. Muestrear alelos HWE por locus y copia conservando fase.
  - Depende de: F2.2.8.5, F2.2.39 y F2.3.3–F2.3.4.
  - Procede de: REQ-CONF-130 y REQ-IND-62; DTC-052 y DTC-072.
  - Cambia: `hwe-composition-v1` realiza P categóricos independientes, ordena la composición y deja la fase a F2.3.7; `alt_af=0/1` es determinista sin consumo y ningún recuento se corrige para acercar frecuencias realizadas.
  - Se comprueba con: corpus bialélico/multialélico y P=1/2/10, 0/intermedias/1, vectores de `next_u64`, distribución degenerada sin consumo, oráculo de probabilidades independiente y ausencia de Poisson/autocorrección.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.3.5`; debe terminar con código 0 y escribir `evidence/f2/F2.3.5/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.3.6. Muestrear genotipos diploides con `F` usando las fórmulas normativas.
  - Depende de: F2.2.8.5, F2.2.39 y F2.3.3–F2.3.4.
  - Procede de: REQ-CONF-130, REQ-IND-62 y REQ-EVO-47/48; DTC-052 y DTC-072.
  - Cambia: `f-diploid-genotype-v1` enumera genotipos no ordenados por índices de alelo, aplica las fórmulas normativas y un categórico; entrega composición a F2.3.7 y no generaliza F.
  - Se comprueba con: probabilidades ajenas para F=0/intermedio/1, orden literal de genotipos, suma/normalización, distribución degenerada sin consumo, multialelia y rechazo no diploide.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.3.6`; debe terminar con código 0 y escribir `evidence/f2/F2.3.6/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.3.7. Asignar fase independiente después de sortear la composición genotípica.
  - Depende de: F2.2.8.3, F2.2.8.5, F2.3.3 y F2.3.5–F2.3.6.
  - Procede de: REQ-IND-62, REQ-IND-63 y REQ-VAL-55; DTC-052, DTC-054 y DTC-072 (resuelve DTPD-118/RAF2-033).
  - Cambia: `phase-permutation-v1` aplica Fisher-Yates a la composición ordenada para heterocigotos diploides, multialélicos y poliploides; P<2 o todos iguales no consume. La fase conocida no representa LD.
  - Se comprueba con: genotipos 0/1, 1/1, 1/2 y poliploides, vectores incluido rechazo de entero acotado, consumos cero y prueba de que cambiar solo fase no vuelve a sortear genotipo ni correlaciona loci.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.3.7`; debe terminar con código 0 y escribir `evidence/f2/F2.3.7/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.3.8. Contrastar SNP explícitos con matriz y modelo antes de modificar copias.
  - Depende de: F2.2.27, F2.3.2 y F2.3.5–F2.3.7.
  - Procede de: REQ-IND-64 y REQ-DIR-02; DTC-020 y DTC-052.
  - Cambia: preflight transaccional que acepta restricciones compatibles y diagnostica probabilidad cero, ALT ajeno al sitio, dosis/genotipo contradictorio, población ausente y densidad SNP individual; no modifica la matriz ni sus frecuencias.
  - Se comprueba con: tabla matriz/modelo/restricción para `alt_af=0/1`, F=1, ALT nuevo, dosis imposible y contradicción tardía sin individuo parcial.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.3.8`; debe terminar con código 0 y escribir `evidence/f2/F2.3.8/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.3.9. Aplicar solo la excepción `--allow-incompatible-genotypes` aprobada.
  - Depende de: F2.3.8 y catálogo de diagnósticos F1.
  - Procede de: DTC-054 y DTC-073 (resuelven DTPD-118 y RAF2-034); REQ-IND-64, REQ-VAL-55 y REQ-MAN-91.
  - Cambia: incorpora conjuntamente parser, validación y semántica de `--allow-incompatible-genotypes` y solo entonces la añade a `available_authorizations`. Confirmarla invalida el digest y obliga a validar; `run` idéntico conserva y registra el genotipo sin relajar otras invariantes.
  - Se comprueba con: antes catalogada/no disponible y después capabilities/parser/validate/run coherentes; confirmación→revalidación, incompatibilidades con/sin flag, edición/conjunto distinto, opción desconocida/disponibilidad falsa y errores no anulables; manifiesto incluye el forzado.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.3.9`; debe terminar con código 0 y escribir `evidence/f2/F2.3.9/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.3.10. Integrar *indels* y los demás cambios directos después de fijar genotipos SNV.
  - Depende de: F2.2.33–F2.2.36, F2.3.4 y F2.3.8–F2.3.9.
  - Procede de: REQ-IND-15, REQ-DIR-01, REQ-DIR-02 y REQ-REC-24; DTC-020 y DTC-031.
  - Cambia: orquestación ordenada por individuo/copia sobre el estado faseado, con transacción lógica y RNG locales para clases no SNV iniciales.
  - Se comprueba con: caso manual matriz+SNP restringido+inserción+deleción, permutación controlada y fallo tardío sin publicación parcial.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.3.10`; debe terminar con código 0 y escribir `evidence/f2/F2.3.10/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.3.11. Materializar y validar cada secuencia haplotípica final contra su historial.
  - Depende de: F2.3.5–F2.3.10 y F2.1.9.
  - Procede de: DTC-003 y DTC-031; REQ-IND-15 y REQ-DIR-01.
  - Cambia: resultado con secuencia final, eventos, fase, dosis y procedencia sin convertir bytes materializados en una segunda verdad.
  - Se comprueba con: aplicación de eventos y ensamblado independiente que producen los mismos bytes/longitudes y conservan `N` por copia.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.3.11`; debe terminar con código 0 y escribir `evidence/f2/F2.3.11/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.3.12. Integrar múltiples individuos sin acoplar matriz, genotipos ni cambios locales.
  - Depende de: F2.3.2–F2.3.11.
  - Procede de: REQ-IND-63, REQ-CONF-22 y REQ-DIR-01; DTC-006 y DTC-038.
  - Cambia: recorrido ordenado de varios individuos con matriz común inmutable y valores/verdad separados por identidad.
  - Se comprueba con: al menos tres individuos, adición/eliminación/reordenación, ejecución física distinta y comparación exacta de matriz e IDs previos.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.3.12`; debe terminar con código 0 y escribir `evidence/f2/F2.3.12/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.3.13. Calcular objetivos y realizados sin autocorregir el muestreo finito.
  - Depende de: F2.3.5–F2.3.12.
  - Procede de: REQ-IND-65 y REQ-VAL-55; DTC-049 y DTC-052.
  - Cambia: agregados canónicos por población/locus de `alt_af` objetivo, copias totales/ALT, `alt_af` realizada y genotipos esperados/realizados, incluyendo individuos forzados.
  - Se comprueba con: cálculo manual sobre corpus pequeño, orden/suma canónicos, desviación finita informada, extremos deterministas y prueba de que ningún genotipo se modifica para acercarla al objetivo.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.3.13`; debe terminar con código 0 y escribir `evidence/f2/F2.3.13/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.3.14. Integrar procedencia de matriz y muestreo en configuración efectiva y manifiesto.
  - Depende de: F2.2.39 y F2.3.2–F2.3.13.
  - Procede de: REQ-MAN-91, REQ-IND-65 y evidencia compartida de REQ-MAN-10; DTC-040, DTC-052 y DTC-075.
  - Cambia: manifiesto con modo/ruta/SHA-256 de matriz y JSON, referencia, poblaciones, parámetros/perfiles/fuente/versión, catálogo/algoritmos RNG, `matrix_candidate_count`/subsemilla/consumo solo cuando hubo Poisson, HWE/F, autorizaciones y valores objetivo/realizados.
  - Se comprueba con: inspección campo a campo, `sha256sum`, enlaces resolubles y ejecuciones preparada/generada exacta/Poisson/reutilizada/forzada; `exact_loci` conserva modo/recuento pero no inventa consumo o subsemilla.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.3.14`; debe terminar con código 0 y escribir `evidence/f2/F2.3.14/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.3.15. Validar estadísticamente HWE/`F` con un corpus y cálculos independientes.
  - Depende de: F2.3.5–F2.3.7 y F2.3.13; no usa F2.3.11 para calcular el esperado.
  - Procede de: REQ-VAL-55 y ACE-116; DTC-047 y DTC-052.
  - Cambia: pruebas de aceptación con probabilidades y umbrales predeclarados para HWE/`F`, fase y aislamiento de subsemillas, separadas de la validez formal F2.2.37.
  - Se comprueba con: muestras reproducibles de tamaños pequeño/grande, comparación de esperado/realizado sin exigir igualdad exacta, límites F=0/1 y fallo deliberado detectado por el oráculo.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.3.15`; debe terminar con código 0 y escribir `evidence/f2/F2.3.15/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.3.16. Completar evidencia y trazabilidad de F2.3.
  - Depende de: F2.3.2–F2.3.15 y puertas de calidad verdes.
  - Procede de: los 13 REQ exclusivos y cinco compartidos asignados a F2.3 en §15/§15.1; ACE-01, ACE-10, ACE-12, ACE-22, ACE-23, ACE-28, ACE-104, ACE-111, ACE-114 y ACE-116.
  - Cambia: matriz de evidencia y escenario sintético reproducible con matriz, genotipos, fase, haplotipos, métricas y manifiesto calculados independientemente.
  - Se comprueba con: §4.3, recuento automático de 13 exclusivos y cinco compartidos sin huérfanos; REQ-IND-02/03 prueban una referencia estructural por individuo, `REQ-CONF-130` queda parcial hacia F6.2, el ID estable omitido se cubre y no se reclama dinámica poblacional, lecturas ni formatos técnicos.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.3.16`; debe terminar con código 0 y escribir `evidence/f2/F2.3.16/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.3.17. Activar matriz, población y creación directa completa en Marimo.
  - Depende de: PAUSA-F2.0 aprobada y F2.3.1–F2.3.16.
  - Procede de: REQ-UI-01 a REQ-UI-07, REQ-UI-13, REQ-UI-14 y DTC-051/DTC-052.
  - Cambia: pestaña `Individuos` para preparar/vincular/verificar referencia, generar/previsualizar/reutilizar matriz, elegir población/perfiles/modos/regiones/HWE/F y confirmar matriz vacía o genotipos incompatibles; cada confirmación forma el conjunto cerrado de autorizaciones, invalida el digest y relanza `validate` antes de habilitar `run` con el mismo conjunto.
  - Se comprueba con: matriz campo→control, SHA/REF y resúmenes por cromosoma/sustitución/población/intervalo, accesibilidad e igualdad UI→TOML→`EffectiveConfig`→artefactos; proceso falso verifica flags idénticos en revalidación/ejecución y edición, retirada o flag general bloquea.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.3.17`; debe terminar con código 0 y escribir `evidence/f2/F2.3.17/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

> **PAUSA-F2.3 — cierre obligatorio de F2.3.** Reconstruir manualmente todos los haplotipos,
> aprobar HWE/`F`, fase, métricas, manifiesto y aislamiento de subsemillas. No
> iniciar F2.4 si Poisson aparece en muestreo de alelos, añadir un individuo
> altera a otro, una excepción limitada oculta otra validación o F2.3.17 deja una
> opción F2.3 solo en TOML. La pausa no puede aprobarse mientras el alcance de
> REQ-CONF-130 se marque como cierre F2 en lugar de evidencia parcial hacia F6.2.

## 9. F2.4 — Verdad estructurada independiente del VCF

**Resultado del hito:** una ejecución directa válida publica `truth.jsonl` v1,
UTF-8 sin comprimir, que basta para reconstruir referencia, individuos, copias,
fase, dosis, eventos y procedencia aunque no se produzca ningún VCF.

- ☐ F2.4.1. Publicar los esquemas por tipo de registro de `truth.jsonl` v1 aprobados por DTC-061.
  - Depende de: PAUSA-F2.3 aprobada y F2.1.3.
  - Procede de: DTPD-102; DTC-031, DTC-035, DTC-036, DTC-061 y DTC-065; REQ-VER-25, REQ-VER-31 y REQ-SAL-06; RAF2-030 resuelto.
  - Cambia: esquemas independientes del modelo Rust para cabecera, copia haplotípica y genotipo; fijan clase/versión, ejecución, referencia, convenio, recuentos/digests, segmentos/modificaciones, fase/procedencia, caso/locus, alelos/dosis/copias y, para cada alelo final proyectable, ALT y lista canónica de casos fuente; fijan también la política ante claves/versiones desconocidas.
  - Se comprueba con: una línea manual válida/inválida de cada tipo, recorrido campo→REQ y validación por implementación de esquema independiente; se rechazan cabecera no inicial, tipo desconocido, clave ausente/duplicada y JSON multilínea.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.4.1`; debe terminar con código 0 y escribir `evidence/f2/F2.4.1/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.4.2. Implementar una proyección pura desde la verdad canónica al modelo externo.
  - Depende de: F2.3.7–F2.3.11 y F2.4.1.
  - Procede de: DTC-003, DTC-004 y DTC-031; REQ-VER-31 y evidencia compartida de REQ-IND-26.
  - Cambia: transformación a registros lógicos `truth.jsonl` que incluye referencia/huella, IDs, copias, cadenas ordenadas de segmentos/modificaciones, ploidía, fase, genotipos, dosis, alelos, coordenadas, `N`, procedencia y relación ALT↔lista ordenada de casos fuente, sin consultar VCF, crear casos compuestos ni copiar el contenedor Rust.
  - Se comprueba con: prueba exhaustiva campo a campo sobre el escenario sintético y propiedad de que desactivar la salida VCF deja idéntica la proyección de verdad.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.4.2`; debe terminar con código 0 y escribir `evidence/f2/F2.4.2/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.4.3. Escribir `truth.jsonl` por flujo, en orden canónico y mediante la publicación F1.
  - Depende de: F2.4.1–F2.4.2 y cierre de publicación/integridad F1.4.
  - Procede de: DTC-035, DTC-040 y DTC-045; REQ-SAL-06.
  - Cambia: escritor UTF-8 sin compresión con cabecera inicial y registros ordenados por referencia/individuo/secuencia/copia y claves internas; calcula por flujo longitudes, recuentos/digests de grupos y SHA-256 del fichero, y publica desde temporal de forma atómica/no-clobber.
  - Se comprueba con: dos escrituras lógicas iguales con bytes/hash iguales, corpus deliberadamente reordenado que se canonicaliza, cabecera/recuentos/digests contrastados, fallo inyectado sin destino completo y `sha256sum` comparado con manifiesto.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.4.3`; debe terminar con código 0 y escribir `evidence/f2/F2.4.3/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.4.4. Construir el validador semántico/reconstructor independiente de `truth.jsonl`.
  - Depende de: F2.4.1; no depende de F2.4.2–F2.4.3 para calcular resultados esperados.
  - Procede de: DTC-047; REQ-VER-31 y REQ-SAL-06; ACE-34 y ACE-62.
  - Cambia: herramienta o implementación de prueba separada que lee por flujo, valida cada esquema, cabecera, orden, unicidad, relaciones, recuentos/digests, reconstruye cada haplotipo y vuelve a calcular genotipos/dosis por caso/locus; registra versión/comando si es externa.
  - Se comprueba con: corpus manual y mutaciones de relación ausente/duplicada/desordenada, segmento, fase, alelo, dosis, copia, recuento y digest; compara secuencias/genotipos precomputados y detecta toda discordancia sin llamar al productor.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.4.4`; debe terminar con código 0 y escribir `evidence/f2/F2.4.4/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.4.5. Demostrar que la verdad permanece disponible sin ninguna salida VCF.
  - Depende de: F2.4.2–F2.4.4.
  - Procede de: DTC-004 y DTC-031; REQ-VER-31 y REQ-SAL-06.
  - Cambia: opción efectiva de salida y recorrido que publica verdad/manifiesto aun con VCF desactivado, sin cambiar el genotipo ni eliminar detalles necesarios.
  - Se comprueba con: ejecuciones gemelas VCF activado/desactivado que tienen la misma verdad lógica y reconstrucción, y ausencia verificada de dependencia desde F2.4.2 hacia cualquier módulo VCF.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.4.5`; debe terminar con código 0 y escribir `evidence/f2/F2.4.5/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.4.6. Integrar `truth.jsonl` en el manifiesto general y el primer conjunto real F2.
  - Depende de: F2.4.3 y contratos cerrados de F1.4/F1.5.
  - Procede de: DTC-039, DTC-040 y DTC-045; REQ-VER-31, REQ-SAL-06 y evidencia compartida de REQ-IND-26.
  - Cambia: extiende el manifiesto general F1 y su catálogo de artefactos para enlazar `truth.jsonl`, sus esquemas, longitud/SHA-256, referencia exacta, configuración efectiva y la primera unidad real de conjunto. No crea manifiestos de evento: los eventos técnicos y sus manifiestos pertenecen a F3/F4.
  - Se comprueba con: recorrido bidireccional de cada ID/ruta/hash, `sha256sum`, rechazo de referencia o artefacto intercambiado y revisión de que no existen rutas huérfanas.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.4.6`; debe terminar con código 0 y escribir `evidence/f2/F2.4.6/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.4.7. Cubrir corrupción, truncado, versión y pérdida de campos sin degradación silenciosa.
  - Depende de: F2.4.3–F2.4.6.
  - Procede de: DTC-031, DTC-036, DTC-044 y DTC-045; REQ-VER-31 y REQ-SAL-06.
  - Cambia: validación y diagnósticos de lectura/publicación para bytes corruptos, JSON incompleto o equivalente del formato elegido, versión desconocida, duplicados y relaciones rotas.
  - Se comprueba con: mutación sistemática del corpus, códigos/ubicaciones estables y prueba de que ningún fallo se acepta como verdad parcial ni deja una ejecución marcada completa.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.4.7`; debe terminar con código 0 y escribir `evidence/f2/F2.4.7/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.4.8. Completar evidencia y trazabilidad de F2.4.
  - Depende de: F2.4.1–F2.4.7 y puertas de calidad verdes.
  - Procede de: REQ-VER-31 y REQ-SAL-06; evidencia compartida de REQ-IND-26 y evidencia parcial directa de REQ-SAL-02; ACE-14, ACE-34 y ACE-62.
  - Cambia: corpus dorado `truth.jsonl`, pruebas de integración, informe del reconstructor independiente y matriz de evidencia F2.4 que separa cierre, evidencia compartida y expediente parcial `REQ-SAL-02 → F8.5`.
  - Se comprueba con: §4.3, recuento automático de dos REQ de cierre, un REQ compartido y REQ-SAL-02 parcial sin huérfanos; validador independiente, coherencia del mismo ID y reconstrucción exacta de todas las copias sin VCF. El expediente parcial no puede cerrar REQ-SAL-02 ni la pausa final asignada a F8.5.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.4.8`; debe terminar con código 0 y escribir `evidence/f2/F2.4.8/result.json` con el paso `status = "passed"`, `requirements.REQ-SAL-02.status = "partial"` y destino `F8.5`, sin comprobaciones omitidas.

- ☐ F2.4.9. Añadir controles y resumen enlazado de verdad estructurada a Marimo.
  - Depende de: PAUSA-F2.0 aprobada y F2.4.1–F2.4.8.
  - Procede de: REQ-UI-01, REQ-UI-07 a REQ-UI-13 y DTC-051; F2.4 de la hoja de ruta.
  - Cambia: `Salidas` muestra `truth.jsonl` obligatorio y no configurable en formato/compresión; `Revisar y ejecutar` resume clase/versión/recuentos/longitud/SHA desde manifiesto sin cargar el fichero completo.
  - Se comprueba con: configuración UI/manual equivalente con VCF activado/desactivado, validación invalidada tras cambio, manifiesto grande simulado y prueba de que la UI solo lee metadatos/enlaces.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.4.9`; debe terminar con código 0 y escribir `evidence/f2/F2.4.9/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

> **PAUSA-F2.4 — cierre obligatorio de F2.4.** Entregar el fichero y la referencia a una revisión
> que no use el productor y exigir reconstrucción exacta. No iniciar F2.5 si el
> VCF es necesario para interpretar la verdad, si falta fase/dosis/procedencia o
> si una corrupción puede parecer una ejecución completa. F2.4.9 debe demostrar que
> la app resume y enlaza la verdad sin cargarla ni interpretarla.

## 10. F2.5 — VCF de verdad mínimo

**Resultado del hito:** cada conjunto publica exactamente un VCF 4.5 de verdad
multimuestra, ordenado por referencia, derivado de `truth.jsonl` y
validado de forma independiente, con GT faseado y sin campos observacionales
inventados.

- ☐ F2.5.1. Congelar y auditar DTC-062 a DTC-066 como resolución completa de DTPD-103.
  - Depende de: PAUSA-F2.4 aprobada.
  - Procede de: DTPD-103; DTC-021 y DTC-061 a DTC-066; REQ-FMT-05 a REQ-FMT-07, REQ-VER-06 y REQ-VER-34.
  - Cambia: expediente de consumo que congela VCF 4.5, compresión/índice, normalización/agrupación, `GPS_*`/complemento y validadores de DTC-066 sin comparar alternativas ni reabrir decisiones.
  - Se comprueba con: recorrido DTC-062..066→arquitectura→pasos, prototipos de todos los modos y enlaces, plantillas exactas de la resolución DTPD-103 y matriz de corpus/mutaciones; ninguna herramienta externa produce o repara la salida.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.5.1`; debe terminar con código 0 y escribir `evidence/f2/F2.5.1/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.5.2. Evaluar e incorporar solo dependencias y validadores aprobados para VCF.
  - Depende de: F2.5.1.
  - Procede de: DTC-043, DTC-047, DTC-048 y DTC-066.
  - Cambia: evaluación de cualquier crate productor bajo DTC-043/048 y provisión CI de `bcftools`/HTSlib-tabix 1.24 como herramientas externas, no dependencias runtime; registra versión y cada argv de DTPD-103. EBI `vcf-validator` queda excluido como autoridad 4.5.
  - Se comprueba con: revisión DTC-043, `cargo metadata --locked --format-version 1`, `cargo tree --locked -e features`, `bcftools --version`, `tabix --version`, §4.3 y ejecución íntegra de las plantillas aprobadas con versión 1.24.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.5.2`; debe terminar con código 0 y escribir `evidence/f2/F2.5.2/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.5.3. Definir la proyección VCF desde verdad final y referencia exacta.
  - Depende de: F2.5.1–F2.5.2 y F2.4.2.
  - Procede de: DTC-003, DTC-004, DTC-021, DTC-061 y DTC-062; REQ-VER-06, REQ-VER-15, REQ-VER-16, REQ-VER-17 y REQ-VER-30.
  - Cambia: modelo intermedio VCF 4.5 con contigs, muestras, posiciones, REF, ALT, genotipos y relación ALT↔lista canónica de casos fuente derivados únicamente de `truth.jsonl` validado y la referencia; no existe proyección inversa VCF→verdad, caso compuesto ni observación inferida.
  - Se comprueba con: fixtures construidos a mano, prueba de dependencia de módulos y comparación de cada registro con referencia/verdad mediante un oráculo distinto.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.5.3`; debe terminar con código 0 y escribir `evidence/f2/F2.5.3/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.5.4. Cerrar la normalización y orden de SNP e indels mediante DTC-064.
  - Depende de: F2.5.4.1–F2.5.4.3.
  - Procede de: REQ-FMT-05, REQ-FMT-06, REQ-VER-15 y REQ-VER-16; DTC-013, DTC-021 y DTC-064; RAF2-036 resuelto.
  - Cambia: integra reconstrucción, recorte anclado y alineación/orden en una puerta del algoritmo interno versionado, sin usar `bcftools norm` para producir la salida.
  - Se comprueba con: aprobación de los tres expedientes hijos y ejecución conjunta de los vectores congelados contra el oráculo externo.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.5.4`; debe terminar con código 0 y escribir `evidence/f2/F2.5.4/result.json` con `status = "passed"` y enlaces a F2.5.4.1–F2.5.4.3.

- ☐ F2.5.4.1. Reconstruir alelos finales y su intervalo REF común.
  - Depende de: F2.5.1 y F2.5.3.
  - Procede de: DTC-061 y DTC-064; verdad haplotípica final.
  - Cambia: reconstrucción por haplotipo y cálculo del intervalo de referencia común suficiente para cambios solapados y alelos de longitudes distintas.
  - Se comprueba con: SNP, inserción, deleción, sustituciones solapadas y combinaciones con longitudes distintas; un reconstructor independiente obtiene las mismas cadenas e intervalo.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.5.4.1`; debe terminar con código 0 y escribir `evidence/f2/F2.5.4.1/result.json` con alelos/intervalos esperados y observados.

- ☐ F2.5.4.2. Recortar contexto común conservando el anclaje VCF.
  - Depende de: F2.5.4.1.
  - Procede de: DTC-064; representación VCF 4.5 DTC-062.
  - Cambia: recorte común determinista que conserva al menos la base de anclaje necesaria y nunca crea REF/ALT vacío o coordenada fuera de referencia.
  - Se comprueba con: primera/interior/última zona admitida, prefijos/sufijos comunes, alelos de distinta longitud y repeticiones; vectores de bordes demuestran anclaje y REF exactos.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.5.4.2`; debe terminar con código 0 y escribir `evidence/f2/F2.5.4.2/result.json` con vectores congelados de recorte.

- ☐ F2.5.4.3. Alinear a la izquierda, ordenar y congelar la versión del algoritmo.
  - Depende de: F2.5.4.2 y catálogo de secuencias F1.
  - Procede de: DTC-013, DTC-021 y DTC-064.
  - Cambia: alineación izquierda respecto al FASTA autoritativo y orden canónico por catálogo/posición, con versión de normalización registrada.
  - Se comprueba con: repeticiones locales y permutaciones de entrada; comparación independiente con `bcftools norm` contra el mismo FASTA sin llamarlo desde el productor ni runtime.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.5.4.3`; debe terminar con código 0 y escribir `evidence/f2/F2.5.4.3/result.json` con versión, vectores y comparación externa.

- ☐ F2.5.5. Agrupar alelos por sitio y calcular GT faseado/dosis para todas las muestras.
  - Depende de: F2.5.3–F2.5.4.
  - Procede de: DTC-003, DTC-021, DTC-064 y DTC-065; REQ-VER-02, REQ-VER-12, REQ-VER-17, REQ-VER-23, REQ-VER-25, REQ-VER-30, REQ-VER-32, REQ-VER-33 y REQ-VER-34; RAF2-030 resuelto.
  - Cambia: agrupación de todos los alelos finales compatibles en exactamente un registro multialélico por `CHROM+POS`, sin descomposición automática ni caso compuesto, con ALT fijos incluidos, índices estables, relación por ALT con la lista ordenada de casos fuente, una columna por individuo y GT con una entrada por copia en orden de fase.
  - Se comprueba con: dos ALT/dos casos, un ALT/dos casos, un caso/varios registros, entrada reordenada, caso duplicado e ID inexistente; además loci bialélicos/multialélicos, cambios solapados compatibles, rechazo de segundo registro, dosis 0..P y ploidías 1/2/10. El oráculo de copia y la permutación de muestras no cambian alelos, coordenadas ni relaciones. Los fixtures F2 no usan BND.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.5.5`; debe terminar con código 0 y escribir `evidence/f2/F2.5.5/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.5.6. Omitir invariantes y alelos eliminados sin perder la verdad canónica.
  - Depende de: F2.5.3–F2.5.5.
  - Procede de: REQ-REC-26, REQ-VER-30, REQ-VER-36, REQ-VER-37, REQ-VER-39, REQ-VER-40 y REQ-VER-41; DTC-021.
  - Cambia: filtros de proyección que excluyen posiciones invariantes, variantes ausentes del conjunto y SNP borrados por una deleción posterior, manteniendo el historial completo en el fichero estructurado.
  - Se comprueba con: verdad manual que contiene alelo fijo, evento revertido o eliminado y sitio invariante; el VCF contiene solo estados finales representables y los registros de modificaciones/procedencia de `truth.jsonl` conservan todo el historial exigido.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.5.6`; debe terminar con código 0 y escribir `evidence/f2/F2.5.6/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.5.7. Aplicar la regla de `N` y emitir el complemento mínimo aprobado cuando el VCF pierda detalle.
  - Depende de: F2.5.1, F2.5.3–F2.5.6.
  - Procede de: DTC-021, DTC-032 y DTC-065; REQ-VER-19, REQ-VER-24 a REQ-VER-26, REQ-VER-34, REQ-VER-47 y REQ-MAN-27; cierre F2 de la hoja de ruta.
  - Cambia: implementa el mecanismo condicional de `vcf-complement.jsonl`, pero F2 no fabrica un caso con pérdida. Rechaza SNP con `REF=N`, conserva `N` representable y exige ausencia de `GPS_LOSS/GPS_COMPLEMENT`, entradas y fichero para todos los casos F2 plenamente representables. Una pérdida literal auténtica descubierta puede añadirse sin ser requisito previo.
  - Se comprueba con: corpus F2 positivo representable, incluida inserción explícita con `N` cuando VCF 4.5 la conserva, inserción aleatoria sin `N`, ausencia exacta de complemento/links y prueba negativa de que `N` nunca se sustituye por A/C/G/T; unitarias del esquema/cardinalidad rechazan relaciones rotas o sobrantes sin inventar una simulación con pérdida. La primera reconstrucción positiva VCF+complemento queda en F5.5.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.5.7`; debe terminar con código 0 y escribir `evidence/f2/F2.5.7/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.5.8. Escribir exactamente un VCF multimuestra por conjunto con cabecera y campos mínimos correctos.
  - Depende de: F2.5.2–F2.5.7.
  - Procede de: DTC-021, DTC-062 y DTC-065; REQ-FMT-05, REQ-FMT-06, REQ-VER-06, REQ-VER-23 a REQ-VER-25, REQ-VER-43, REQ-VER-44 y REQ-MAN-57.
  - Cambia: adaptador VCF 4.5 con `GPS_CASE` `Number=.`/`Type=String` en todo registro, IDs originales únicos en orden canónico y sin caso compuesto; `GPS_LOSS` de catálogo cerrado y `GPS_COMPLEMENT` solo ante pérdida futura; `QUAL=.`, `FILTER=PASS`, GT faseado y ausencia de campos observacionales no justificados. La cabecera declara número/tipo/descripción de cada `GPS_*`.
  - Se comprueba con: parser independiente de cabecera/registros; dos ALT/dos casos, un ALT/dos casos, un caso/varios registros, orden de entrada, duplicados e ID inexistente; ausencia F2 exacta de `GPS_LOSS/GPS_COMPLEMENT` y complemento; campos prohibidos, cardinalidad por conjunto y orden de muestras/contigs.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.5.8`; debe terminar con código 0 y escribir `evidence/f2/F2.5.8/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.5.9. Publicar BGZF e índice por defecto y texto solo cuando se solicite explícitamente.
  - Depende de: F2.5.1–F2.5.2, F2.5.8 y publicación F1.4 cerrada.
  - Procede de: REQ-SAL-12 y REQ-MAN-58 como fuentes directas de BGZF, índice y modo texto; REQ-SAL-21 a REQ-SAL-23 y evidencia compartida de REQ-MAN-83 para los estados y la publicación; DTC-040, DTC-045 y DTC-063.
  - Cambia: artefactos temporales `.vcf.gz` y exactamente un `.tbi` si todo cabe en el límite TBI o un `.csi` si alguna secuencia/posición lo supera; alternativa explícita `.vcf` sin compresión ni índice; hashes, `EffectiveConfig` y manifiestos con modo, tipo y rutas efectivos. VCF+índice se validan y promueven como una unidad coherente o no se publica ninguno.
  - Se comprueba con: validador externo que abre BGZF y el único índice, referencias a ambos lados del límite TBI, consultas por región, ausencia del índice alternativo, modo texto sin índice, `sha256sum` y fallos inyectados en VCF/índice sin conjunto parcialmente publicado.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.5.9`; debe terminar con código 0 y escribir `evidence/f2/F2.5.9/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.5.10. Ejecutar validación sintáctica independiente y validación semántica contra la verdad.
  - Depende de: F2.5.7–F2.5.9.
  - Procede de: DTC-047 y DTC-066; REQ-VAL-28, REQ-VER-02, REQ-VER-12, REQ-VER-23 a REQ-VER-26, REQ-VER-43, REQ-VER-44 y REQ-VER-47.
  - Cambia: aceptación/CI ejecuta lectura, `bcftools norm` contra FASTA y consultas regionales con bcftools/tabix 1.24; el comprobador semántico separado lee FASTA, `truth.jsonl`, VCF y complemento opcional sin llamar ni compartir código con la proyección productora.
  - Se comprueba con: corpus positivo texto/BGZF+TBI/BGZF+CSI; dos ALT/dos casos, un ALT/dos casos, un caso/varios registros y orden de entrada; consultas esperadas; ausencia del índice alternativo; mutaciones individuales de REF/ALT/GT, dosis/fase/muestras, orden, unicidad `CHROM+POS`, `GPS_CASE` duplicado/no canónico/ID inexistente/relación ALT↔casos alterada, presencia indebida F2 de `GPS_LOSS/GPS_COMPLEMENT` o complemento y campos prohibidos, cada una detectada por el oráculo independiente.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.5.10`; debe terminar con código 0 y escribir `evidence/f2/F2.5.10/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.5.11. Verificar el recorrido completo configuración→verdad→VCF y sus invariantes cruzadas.
  - Depende de: F2.5.3–F2.5.10.
  - Procede de: F2.5 y cierre F2; DTC-003, DTC-004, DTC-021 y DTC-031.
  - Cambia: escenario de aceptación con varias muestras, ploidías del alcance, SNP, inserción/deleción, sitio multialélico, alelo fijo y `N`, más inventario de artefactos y métricas pequeñas de tiempo/memoria.
  - Se comprueba con: reconstrucción independiente de todos los haplotipos, validación VCF externa, contraste fase/dosis/alelos/coordenadas y repetición con igual semilla que conserva verdad lógica e integridad esperada.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.5.11`; debe terminar con código 0 y escribir `evidence/f2/F2.5.11/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

- ☐ F2.5.12. Completar evidencia y trazabilidad de F2.5 y preparar la revisión de cierre F2.
  - Depende de: F2.5.1–F2.5.11 y puertas de calidad verdes.
  - Procede de: los 25 REQ de cierre asignados a F2.5 en §15, evidencia compartida de REQ-MAN-83 y evidencia parcial directa de REQ-SAL-12, REQ-FMT-07, REQ-VER-23 a REQ-VER-26, REQ-MAN-31 y REQ-MAN-58; ACE-19, ACE-25, ACE-28, ACE-29, ACE-30, ACE-60, ACE-61, ACE-62, ACE-63, ACE-64, ACE-65, ACE-66, ACE-67, ACE-93, ACE-94, ACE-100, ACE-103 y ACE-108.
  - Cambia: matriz de evidencia F2.5, informe de validación independiente, inventario de limitaciones y dossier de cierre; separa el REQ compartido y los ocho expedientes parciales directos de los cierres exclusivos.
  - Se comprueba con: §4.3, recuento automático de 25 REQ de cierre, un REQ compartido y ocho parciales directos sin huérfanos; comandos de validadores y revisión cruzada de §15. Ningún expediente de formato/VCF/complemento/manifiesto cierra antes de F5.5 o F9.5.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.5.12`; debe terminar con código 0 y escribir `evidence/f2/F2.5.12/result.json` con el paso `status = "passed"`, los ocho `requirements.<REQ>.status = "partial"` y sus destinos F5.5/F9.5, sin comprobaciones omitidas.

- ☐ F2.5.13. Incorporar todas las opciones VCF F2.5 a Marimo y cerrar la paridad F2.
  - Depende de: PAUSA-F2.0 aprobada y F2.5.1–F2.5.12.
  - Procede de: REQ-UI-01, REQ-UI-03, REQ-UI-07 a REQ-UI-13 y DTC-051; F2.5 de la hoja de ruta.
  - Cambia: controles `Salidas` para VCF BGZF/indexado predeterminado o texto, revisión de artefactos esperados, validación/lanzamiento explícitos y resumen de manifiesto/VCF sin cargar registros completos.
  - Se comprueba con: pares UI/manual para ambos modos con igual `EffectiveConfig`, destino existente, edición tras validar, ejecución real mínima y matriz de todas las opciones configurables F2 sin huecos Marimo.
  - Evidencia ejecutable: `./scripts/verify-f2-step F2.5.13`; debe terminar con código 0 y escribir `evidence/f2/F2.5.13/result.json` con `status = "passed"`, o `status = "blocked"` únicamente mientras subsista una decisión pendiente nombrada, sin comprobaciones omitidas.

> **PAUSA-F2.5 — cierre obligatorio de F2.5 y de F2.** No autorizar F3 hasta que una
> revisión independiente reconstruya todos los haplotipos, confirme
> coordenadas/alelos/fase/dosis, valide VCF/BGZF/índice y demuestre los cuatro
> casos de `N`: SNP rechazado, inserción aleatoria sin `N`, inserción explícita
> con `N` conservada, ausencia del complemento no necesario y ninguna resolución
> inventada. La primera prueba positiva VCF+complemento se exige en F5.5 sobre
> un caso complejo/SV/BND real. Aprobar por
> separado el cierre de F2; completar F2.5.12/F2.5.13 no lo concede automáticamente y
> ninguna opción configurable F2 puede quedar solo en TOML.

## 11. Dependencias posteriores confirmadas, fuera de la ejecución de F2

Esta sección no añade hitos ni casillas a F2. Conserva el destino explícito de
DTC-053 sin anticipar la planificación completa de F6 o F7:

| Destino | Contrato que deberá consumir | Puertas previas | Evidencia mínima futura |
|---|---|---|---|
| F5.5, primer complemento positivo | DTC-065/066; REQ-VER-19 y REQ-VER-23 a REQ-VER-26 | Mecanismo condicional y negativas F2.5 cerrados; clase compleja/SV/BND real implementada | Producir un caso auténtico que pierda detalle, enlazar todos sus registros mediante `GPS_CASE/GPS_LOSS/GPS_COMPLEMENT` a exactamente una entrada, reconstruir fase/haplotipo/genotipo/coordenadas/cambio y probar ausencia de entradas sobrantes; no fabricar pérdida ni reinterpretar un caso F2 representable |
| F6.4, panel fundador intacto o recombinante | REQ-EVO-74 a REQ-EVO-77; REQ-VAL-56, REQ-MAN-92 y REQ-UI-15 | Consumir íntegramente DTC-055: paquete/TSV, referencia única, mínimos, ejes/modelos, productor directo, constructor común uniforme, puntos x\|x+1, F, IDs DTC-039, RNG DTC-037, origen/proceso y salidas canónicas | Probar catálogo/ejes y rechazos de combinaciones, incluida distribución explícita solo en intacto; ejecutar Bernoulli(F), productores y fuentes solo del panel; probar referencia/secuencias, K, puntos no SNP/invisibles; verificar C, ID, `n_h/C`, recuentos por origen y los tres esquemas/órdenes/claves TSV; reconstruir con oráculo independiente y corrupción deliberada; publicar/enlazar las tres tablas; probar Marimo, validación vigente y paridad TOML |
| F7.1, recombinación meiótica | REQ-REC y requisitos reproductivos aplicables | F6.4 cerrado y DTPD-101/DTPD-111 resueltas | Conectar los mapas y la semántica meiótica decididos por DTPD-111 al constructor común ya implementado en F6.4 y demostrar que no reinterpreta la inicialización fundadora como evento histórico ni crea una segunda implementación |

DTPD-119 está resuelta mediante DTC-055. Quedan decididos el paquete y su preparación, el
alcance SNP, bloques e identidades locales, tablas y vectores densos/únicos, los
modelos y modos del muestreo intacto, la asignación a copias, el uso de la matriz
para loci exteriores y el rechazo de cualquier locus `alt_af` dentro de un
bloque. La alternativa recombinante se ha trasladado a la formación del conjunto
fundador F6.4 y ya tiene productor directo, puntos, F, IDs, RNG, tres salidas y
Marimo definidos, incluida la sintaxis TSV, la semántica de origen/proceso y el
rechazo de distribuciones genotípicas explícitas en recombinante.
Ningún paso F2.2/F2.3
puede aceptar un panel,
inferir haplotipos desde la matriz marginal ni implementar recombinación local.
Marimo solo expondrá REQ-UI-15 cuando cierren las capacidades correspondientes.

## 12. Riesgos genómicos silenciosos y oráculos obligatorios

| Riesgo silencioso | Prevención en producción | Oráculo independiente obligatorio | Evidencia de fallo que debe detectarse |
|---|---|---|---|
| Desfase base 0/base 1 o extremo abierto/cerrado | Tipos/servicio DTPD-101 iniciados en F1.6 y extendidos en F2.1.5 | Vectores compartidos de mapa y cadenas etiquetadas por posición | Mapa, SNP o indel desplazado, ancla incorrecta o acceso a posición 0 |
| Longitud de *indel* mal resuelta o convertida en distribución implícita | Parámetros exactos F2.2.2/F2.2.7 y constructores F2.2.33–F2.2.34 | Casos 1/2/`umbral - 1` y RNG instrumentado | Longitud distinta, consumo RNG para longitud o clasificación estructural incorrecta |
| Inserción/deleción terminal inválida | Constructores validados F2.2.33–F2.2.34 | Enumeración manual de todos los indels de una cadena corta | Pérdida del ancla o secuencia sin base a uno de los lados |
| Sustitución silenciosa de `N` | Dominio F2.1.6 y reglas F2.2.9/F2.2.17/F2.2.22–F2.2.34/F2.5.7 | Corpus con `N` en borde/interior y barrido de salidas | SNP `REF=N`, ALT `N` o inserción aleatoria con `N` |
| Fase o dosis alterada | Copias ordenadas F2.1.7/F2.1.9 y GT F2.5.5 | Conteo manual por copia, no parsing circular del VCF | GT sin P alelos, separador no faseado o dosis distinta |
| Evento perdido por solapamiento | Composición F2.2.35 e historial F2.4.2 | Aplicador de cadenas corto escrito desde fixtures | SNP superviviente tras ser borrado o conflicto aceptado |
| Tasa cero que consume RNG | Rutas explícitas F2.2.11/F2.2.20/F2.2.24 | RNG instrumentado y comparación de ámbitos vecinos | Cambios en otras clases o individuos al activar cero |
| Poisson incorrecto o aplicado a alelos/genotipos | Decisión F2.2.2 y servicio F2.2.11/F2.3.5–F2.3.6 | Vectores o implementación ajena y test estadístico prefijado | Media/varianza fuera del umbral, Poisson fuera del recuento o resultado no reproducible |
| Intervalo continuo produce 0/1 o interpreta 1 como cerrado | DTC-072/075 y F2.2.8.4/F2.2.12/F2.2.24 | Vectores independientes de cero, subnormales, `next_down(1)` y redondeos | `alt_af` extremo sin componente `fixed` o perfil SFS incompatible |
| Selector estático aplicado a indels sobre estado mutable | `dynamic-indel-selection-v1` en F2.2.8.3/F2.2.30–F2.2.31 | Modelo de cadenas que gana/pierde oportunidades y agota el espacio | Preflight falso, candidato obsoleto, prefijo publicado o nuevo Poisson |
| Matriz o espectro incompatibles con referencia/población | Esquemas/validación DTC-054 y F2.2.3–F2.2.4/F2.2.12/F2.2.17–F2.2.21 | Corpus/lector F2.2.37 calculado fuera del productor | REF/ALT, suma, duplicado, población, perfil o digest incorrectos aceptados |
| Recuento candidato confundido con loci activos | Separación F2.2.11/F2.2.20/F2.2.25 | Tabla manual candidatos→ALT cero→activos | Reposición de un locus globalmente ausente o `exact_loci` interpretado como activos |
| Distribución o perfil implícito/incorrecto | Perfiles F2.2.10/F2.2.12–F2.2.15 y selección F2.2.22–F2.2.24 | Oráculos F2.2.38 y regeneración compacta DTC-068 | Solapamiento aceptado, mínimo `1/5008` alterado, promedio/complemento incorrecto o builtin aproximado cargado por defecto |
| Identidad o azar cambia por orden/adición | IDs F2.2.16 y ámbitos F2.2.8/F2.3.3 | Permutaciones y vectores literales | Reordenar filas/poblaciones o añadir individuo altera resultados existentes |
| Excepción limitada relaja otra validación | Flags F2.2.21/F2.3.9 y preflight F2.2.36 | Matriz independiente error×flag | Digest, REF, frecuencia, ploidía/dosis o no-clobber inválidos aceptados |
| HWE/`F`, fase o frecuencia realizada incorrectos | Muestreo F2.3.5–F2.3.9 y agregados F2.3.13 | Probabilidades/corpus independiente F2.3.15 | Genotipo fuera del modelo, fase cambiada o muestra finita autocorregida |
| Verdad y secuencia materializada divergen | Una verdad canónica F2.1.8–F2.1.9/F2.3.11 | Segundo reconstructor de F2.4.4 | Bytes o longitud distintos para la misma copia |
| VCF convertido en fuente de verdad | Dirección de dependencias F2.4.2/F2.5.3 | Prueba sin VCF y auditoría de módulos | No se puede reconstruir o cambia la verdad al desactivar VCF |
| REF/ALT/normalización incorrectos | Proyección F2.5.3–F2.5.4 | Normalizador/validador independiente aprobado | Registro válido sintácticamente pero no representa la verdad |
| Pérdida de detalle al agrupar multialélicos | Agrupación F2.5.5 y complemento F2.5.7 | Tabla manual sitio×copia y reconstructor combinado | ALT o copia portadora desaparece, o se inventa una base |
| Artefactos mezclados entre ejecuciones | IDs/hashes/publicación F1 reutilizados | `sha256sum`, referencias cruzadas y fallos inyectados | Índice, VCF, verdad o referencia de ejecuciones distintas |
| UI y conjunto TOML producen semánticas distintas, rompen enlaces o ejecutan una revisión obsoleta | Modelo/digest multidocumento F2.0.10/F2.0.13, binario único F2.0.14/F2.0.15 y pasos F2.3.17/F2.4.9/F2.5.13 | Conjunto manual gemelo, escenarios faltante/duplicado y proceso falso que registra argumentos | `EffectiveConfig` distinta, ruta/orden perdido, ejecución automática o validación antigua reutilizada |

Una prueba que serializa y vuelve a leer con el mismo código solo demuestra
consistencia interna. No satisface por sí sola ninguna fila de esta tabla.

## 13. Cambios pequeños, revisión y política de bloqueo

- Cada paso se implementa y revisa de forma aislada; una revisión puede dividirlo
  aún más, pero no fusionarlo con un paso dependiente para ocultar su evidencia.
- Si F2.0.2, F2.0.13, F2.1.2, F2.1.3, F2.2.2, F2.2.14, F2.3.1 o F2.5.1 no obtienen aprobación, o F2.2.1 no confirma sin
  contradicción DTC-054 y el cierre de DTPD-118, se detienen sus consumidores y se
  actualiza este plan con el bloqueo. Se puede continuar únicamente con corpus,
  matrices o validadores que no presupongan la opción elegida.
- F2.0.2, F2.2.1 y F2.2.14 son pasos internos ejecutables después de autorizar
  el inicio de F2, no prerrequisitos de esa autorización: F2.0.2 bloquea
  F2.0.3 y sus demás consumidores, F2.2.1 bloquea únicamente los consumidores
  declarados y F2.2.14 no incorpora los perfiles hasta materializar y revisar
  los literales SFS, tamaño/SHA-256 de la fuente eLife, cinco filas, 96
  contextos, transformaciones y SHA-256 de salida.
- Si una decisión aprobada contradice especificación, arquitectura o DTC
  confirmada, no se adapta el requisito al código: se pausa y se eleva la
  contradicción para decisión explícita.
- Las mediciones pequeñas observan regresiones y consumo, pero no prometen los
  objetivos de escalabilidad de F9 ni justifican optimizaciones especulativas.
- Ninguna pausa se supera solo porque compila. Requiere revisión de contratos,
  pruebas negativas, oráculo independiente y matriz REQ/ACE actualizada.

## 14. Validación mecánica del documento

Antes de aprobar el plan se ejecutarán, como mínimo:

```bash
git diff --check -- design/plan-implementacion-f2-0.1.md
rg -n '^- ☐ F2\.[0-5]\.[0-9]+(\.[0-9]+)*\.' design/plan-implementacion-f2-0.1.md
rg -n '^> \*\*PAUSA-F2\.[0-5]' design/plan-implementacion-f2-0.1.md
rg -n 'DTC-052|DTC-053|DTC-054|DTC-055|DTC-058|DTC-059|DTC-060|DTC-061|DTC-062|DTC-063|DTC-064|DTC-065|DTC-066|DTC-067|DTC-068|DTC-069|DTC-070|DTC-071|DTC-072|DTC-073|DTC-074|DTC-075|DTPD-101|DTPD-102|DTPD-103|DTPD-117|DTPD-118|DTPD-119|DTPD-120|DTPD-121|F1\.5|F1\.6' design/plan-implementacion-f2-0.1.md
rg -n 'RAF2-001|RAF2-004|RAF2-007|RAF2-012|RAF2-017' design/plan-implementacion-f2-0.1.md
```

Además, un comprobador mecánico debe verificar que cada casilla contiene una y
solo una línea `Depende de`, `Procede de`, `Cambia`, `Se comprueba con` y
`Evidencia ejecutable`; que existen exactamente 130 casillas; que la
matriz siguiente contiene una sola vez los 78 REQ cuyo cierre está asignado
exclusivamente a F2; que la lista separada de evidencia compartida contiene
los otros catorce REQ cuyo campo de hito incluye F2 junto con fases posteriores;
y que la lista directa contiene nueve obligaciones parciales cuyo hito final
no incluye F2. Los dos conjuntos parciales deben permanecer fuera de los 78
cierres y declarar `status = "partial"` y destino final.

## 15. Matriz de trazabilidad REQ/ACE/DTC

La lista de REQ de cierre exclusivo es exhaustiva para la asignación F2 vigente:
13 + 3 + 22 + 13 + 2 + 25 = **78**. Se añaden por separado **catorce REQ compartidos**
y **nueve REQ parciales directos** que exigen evidencia en F2 y cierre posterior. Los ACE son la unión de
los criterios enlazados por esas filas de trazabilidad. Una mención no acredita
ejecución ni cierre. REQ-CONF-08, REQ-CONF-21, REQ-IND-02, REQ-IND-03 y
REQ-CONF-130 figuran como evidencia parcial con destino posterior; RAF2-004 y
RAF2-012 quedan resueltos sin ampliar F2. La auditoría de la matriz maestra asigna los doce REQ
nuevos como cuatro a F2.2 y ocho a F2.3; no se adopta el reparto 5/7 porque
contradiría esas filas sin una revisión explícita de trazabilidad.

| Hito | REQ asignados | ACE cubiertos | Decisiones/puertas principales | Pasos de cierre |
|---|---|---|---|---|
| F2.0 | REQ-UI-01; REQ-UI-02; REQ-UI-03; REQ-UI-04; REQ-UI-05; REQ-UI-06; REQ-UI-07; REQ-UI-08; REQ-UI-09; REQ-UI-10; REQ-UI-11; REQ-UI-12; REQ-VAL-54; REQ-UI-13 parcial local en §15.1 | ACE-114, ACE-115 | DTC-034, DTC-036, DTC-043, DTC-044, DTC-045, DTC-048, DTC-051, DTC-058, DTC-069 y DTC-073; DTC-074/RAF2-035 cierran en F9 | F2.0.0–F2.0.18 |
| F2.1 | REQ-IND-28; REQ-IND-29; REQ-VER-27 | ACE-03, ACE-10, ACE-11, ACE-28, ACE-66 | DTC-003, DTC-010, DTC-013, DTC-031, DTC-032, DTC-039, DTC-056 y DTC-061; DTPD-102 resuelta; REQ-REF-14 ya cerrado en F1.5 se conserva como invariante heredada | F2.1.2–F2.1.10 |
| F2.2 | REQ-CONF-24; REQ-CONF-25; REQ-CONF-28; REQ-CONF-128; REQ-CONF-129; REQ-IND-05; REQ-IND-06; REQ-IND-10; REQ-IND-11; REQ-IND-13; REQ-IND-14; REQ-IND-17; REQ-IND-18; REQ-IND-19; REQ-IND-21; REQ-IND-38; REQ-IND-40; REQ-IND-41; REQ-IND-60; REQ-IND-61; REQ-DIR-03; REQ-REC-24 | ACE-01, ACE-09, ACE-10, ACE-21, ACE-22, ACE-28, ACE-85, ACE-95, ACE-103, ACE-104, ACE-109, ACE-110, ACE-111, ACE-116 | DTC-006, DTC-011, DTC-020, DTC-022, DTC-032, DTC-033, DTC-037, DTC-038, DTC-044, DTC-049, DTC-052, DTC-054, DTC-056, DTC-059/DTC-060, DTC-067/DTC-068 y DTC-070 a DTC-075; DTPD-118/120/121 y RAF2-028/029/033/041–044 resueltas; puerta compacta F2.2.14 sin VCF | F2.2.1–F2.2.40 |
| F2.3 | REQ-CONF-20; REQ-CONF-22; REQ-CONF-38; REQ-IND-15; REQ-IND-62; REQ-IND-63; REQ-IND-64; REQ-IND-65; REQ-DIR-01; REQ-DIR-02; REQ-UI-14; REQ-VAL-55; REQ-MAN-91 | ACE-01, ACE-10, ACE-12, ACE-22, ACE-23, ACE-28, ACE-104, ACE-111, ACE-114, ACE-116 | DTC-003, DTC-006, DTC-011, DTC-020, DTC-037, DTC-038, DTC-039, DTC-049, DTC-051, DTC-052, DTC-054 y DTC-071 a DTC-073; DTPD-118 y RAF2-028/033/034 resueltas; REQ-CONF-130 parcial hacia F6.2 | F2.3.1–F2.3.17 |
| F2.4 | REQ-VER-31; REQ-SAL-06 | ACE-34, ACE-62 | DTC-003, DTC-004, DTC-031, DTC-035, DTC-036, DTC-040, DTC-045, DTC-047, DTC-051 y DTC-061; DTPD-102 resuelta | F2.4.1–F2.4.9 |
| F2.5 | REQ-FMT-05; REQ-FMT-06; REQ-REC-26; REQ-VER-02; REQ-VER-06; REQ-VER-12; REQ-VER-15; REQ-VER-16; REQ-VER-17; REQ-VER-30; REQ-VER-32; REQ-VER-33; REQ-VER-34; REQ-VER-36; REQ-VER-37; REQ-VER-39; REQ-VER-40; REQ-VER-41; REQ-VER-43; REQ-VER-44; REQ-VER-47; REQ-VAL-28; REQ-MAN-27; REQ-MAN-57; REQ-SAL-21 | ACE-19, ACE-25, ACE-28, ACE-29, ACE-30, ACE-60, ACE-61, ACE-62, ACE-63, ACE-64, ACE-66, ACE-67, ACE-93, ACE-94, ACE-100, ACE-103 | DTC-003, DTC-004, DTC-013, DTC-021, DTC-032, DTC-040, DTC-043, DTC-045, DTC-047, DTC-051, DTC-061 a DTC-066; DTPD-103 y RAF2-030 resueltas | F2.5.1–F2.5.13 |

### 15.1. Evidencia compartida obligatoria en F2

Estas catorce filas tienen un destino compuesto en la matriz maestra. F2 debe
producir la evidencia concreta indicada, pero no declarará el REQ completamente
satisfecho hasta que cierren también sus fases posteriores:

| REQ compartido | Evidencia parcial exigida en F2 | Pasos |
|---|---|---|
| REQ-CONF-08 (F2.2/F5.2/F7.2) | Configurar cambios concretos/aleatorios solo para matriz SNV e *indels* F2; SV y mutación heredable quedan para los destinos indicados | F2.2.5, F2.2.25, F2.2.33–F2.2.35, F2.2.40 |
| REQ-CONF-21 (F2.2/F5.2/F7.2) | Probar distribución espacial de matriz e *indels* sin reclamar SV ni mutación heredable | F2.2.10, F2.2.29, F2.2.33–F2.2.34, F2.2.40 |
| REQ-CONF-106 (F2.2/F6.2) | Normalizar las primeras colecciones reales F2 y registrar tolerancia, sumas, delta, decisión y valores antes/después; F6.2 amplía el mismo contrato a frecuencias poblacionales | F2.2.6.1, F2.2.19, F2.2.40 |
| REQ-CONF-107 (F2.2/F6.2) | Distinguir normalización general de tolerancias específicas y conservar la decisión efectiva; F6.2 la aplica a poblaciones | F2.2.6.1, F2.2.19, F2.2.40 |
| REQ-VAL-40 (F2.2/F6.2) | Oráculo independiente verifica suma, tolerancia, normalización y rechazo F2; F6.2 añade el modo poblacional controlado | F2.2.6.1, F2.2.19, F2.2.40 |
| REQ-MAN-81 (F2.2/F6.2) | `EffectiveConfig` y manifiesto F2 conservan suma, valores originales, tolerancia, decisión y valores normalizados; F6.2 completa población | F2.2.6.1, F2.2.19, F2.2.40 |
| REQ-IND-02 (F2.3/F5.2) | Crear individuos sobre la única referencia estructural F2; múltiples referencias estructurales cierran en F5.2 | F2.3.2, F2.3.4, F2.3.10–F2.3.11, F2.3.16 |
| REQ-IND-03 (F2.3/F5.2) | Conservar ID/huella de la referencia estructural F2 por muestra sin reclamar todas las referencias disponibles | F2.3.2, F2.3.4, F2.3.11, F2.3.14, F2.3.16 |
| REQ-CONF-130 (F2.3/F6.2) | HWE P=1..10 y F solo diploide; rechazo fuera de P=2 y ausencia de distribución explícita en esquema/CLI/Marimo | F2.3.1, F2.3.4, F2.3.16 |
| REQ-IND-25 (F2.3/F6.1) | Un creador directo real asigna un ID estable cuando falta, comprueba unicidad y no reutilización dentro del recorrido disponible | F2.3.2, F2.3.16 |
| REQ-IND-26 (F2.4/F3.5/F7.7/F8.5) | El mismo ID aparece sin traducción en estado genotípico, verdad estructurada y metadatos F2; observaciones y linaje quedan para sus fases | F2.4.2, F2.4.6, F2.4.8 |
| REQ-MAN-10 (F2.3/F5.5/F9.7) | El manifiesto distingue valores globales, específicos y efectivos consumidos por individuos/regiones/clases de F2 | F2.3.2, F2.3.14, F2.3.16 |
| REQ-MAN-83 (F2.5/F3.5) | El primer conjunto real y sus artefactos VCF/índice manifiestan estados y diagnósticos coherentes con la ejecución, el conjunto y los artefactos; el estado de evento se añade y comprueba en F3.5 | F2.5.9, F2.5.12 |
| REQ-UI-13 (F2.0/F9.7) | Probar instalación desde repositorio, arranque local offline, accesibilidad y paridad UI; F9.7 acredita distribución y plataforma oficial | F2.0.2–F2.0.6, F2.0.18 |

### 15.2. Evidencia parcial directa con cierre posterior

F2 debe registrar además nueve expedientes `status = "partial"` directamente
reutilizables. No cambian el hito de cierre de los REQ ni permiten aprobarlos:

| REQ | Evidencia parcial ejecutable en F2 | Pasos | Cierre final |
|---|---|---|---|
| REQ-SAL-02 | Publicar para creación directa el fichero de verdad F2, reconstruir cada copia sin VCF y enlazar identidad/integridad. `./scripts/verify-f2-step F2.4.8` debe producir el expediente parcial. | F2.4.2, F2.4.5, F2.4.8 | F8.5, tras cubrir todos los modos y estados exigidos |
| REQ-SAL-12 | Publicar el VCF F2 solicitado como BGZF con índice, abrirlo y consultarlo con el validador fijado; `./scripts/verify-f2-step F2.5.12` debe producir el expediente parcial. | F2.5.9, F2.5.10, F2.5.12 | F5.5, tras cubrir el contrato VCF completo |
| REQ-FMT-07 | Leer completamente VCF texto/BGZF, normalizar contra FASTA y consultar regionalmente con las herramientas fijadas. | F2.5.2, F2.5.10, F2.5.12 | F9.5, tras cubrir todas las vías |
| REQ-VER-23 | Emitir y verificar `GPS_CASE` multivaluado, único, ordenado y enlazado con `truth.jsonl`. | F2.5.5, F2.5.8, F2.5.10, F2.5.12 | F5.5, al incorporar campos de repeticiones/SV |
| REQ-VER-24 | Demostrar que todos los casos representables F2 omiten complemento, sin fabricar pérdida. | F2.5.7, F2.5.10, F2.5.12 | F5.5, con el primer caso positivo real complejo/SV/BND |
| REQ-VER-25 | Verificar ALT↔casos, ausencia condicional de pérdida y cardinalidades negativas. | F2.5.5, F2.5.7, F2.5.8, F2.5.10, F2.5.12 | F5.5, con cardinalidad positiva de complemento real |
| REQ-VER-26 | Validar el mecanismo de reconstrucción y su ausencia en casos sin pérdida, sin reclamar caso positivo F2. | F2.5.7, F2.5.10, F2.5.12 | F5.5, al reconstruir un caso complejo real perdido |
| REQ-MAN-31 | Validar que ruta/tamaño/SHA de `vcf-complement.jsonl` solo existen cuando el fichero existe. | F2.5.7, F2.5.9, F2.5.12 | F5.5, con todas las clases complejas |
| REQ-MAN-58 | Registrar versión VCF, normalización, modo, índice, rutas, validación y publicación atómica. | F2.5.8–F2.5.10, F2.5.12 | F9.5 |

Los casos de aceptación de F2 conservan también evidencia adyacente sobre `N`,
inserciones aleatorias y salida directa, pero no incorporan como cerrados
REQ-IND-58, REQ-IND-59 ni REQ-VER-01.

## 16. Criterio de salida del plan

F1 está cerrado y la validación mecánica de esta corrección debe detectar
errores —incluidos los 78 REQ de cierre únicos, los catorce compartidos, los nueve
parciales directos, las seis pausas y las cinco líneas de cada uno de los 130
pasos—. El inicio de F2 no puede aprobarse todavía: requiere nueva revisión
adversarial sobre este texto corregido y aprobación explícita del usuario.

F2.0.2, F2.2.1 y la regeneración/revisión compacta DTC-068 dentro de F2.2.14 no
son prerrequisitos de inicio: son puertas de ejecución internas con los
consumidores delimitados en §13. DTPD-103/120/121 y
RAF2-001/004/007/012/017 quedan resueltos por las decisiones y reasignaciones
aprobadas. F2 solo podrá cerrarse cuando F2.5.12/F2.5.13, la pausa final y la
aprobación independiente de cierre consten por separado.
