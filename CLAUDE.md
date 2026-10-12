# Guía para Claude Code

Autoridad única sobre política de agentes y forma de trabajo en este
repositorio. **Sustituye a `AGENTS.md`**, que queda como histórico de la etapa
de Codex y no debe consultarse para decidir nada.

Las secciones 1 a 5 son genéricas y **exportables a otros proyectos** tal cual.
De la 6 en adelante es específico de este repositorio.

---

## 1. Enrutamiento de modelos

El agente raíz elige **explícitamente** modelo y esfuerzo de cada subagente. No
se hereda el modelo por inercia.

Antes de crear un subagente:

1. comprobar si basta con una búsqueda, un comando o una edición pequeña; en ese
   caso **no** crear el subagente;
2. acotar una responsabilidad independiente y verificable;
3. elegir el nivel **mínimo** que preserve la calidad necesaria;
4. enviar un contexto autocontenido y reducido;
5. decir al usuario qué nivel se eligió cuando la elección sea material.

| Tipo de trabajo | Modelo | Esfuerzo |
|---|---|---|
| Estado, inventario, búsquedas acotadas, ejecución de comandos, formato, documentación mecánica | Haiku 4.5 | bajo o medio |
| Exploración concreta del código | Haiku 4.5, o Sonnet 5 si hay que razonar sobre lo hallado | medio |
| Implementación ordinaria, correcciones ya diagnosticadas, pruebas, planificación, documentación técnica | Sonnet 5 | medio |
| Implementación compleja entre módulos, corrección de un defecto grave | Opus 5 | alto |
| Arquitectura material, decisiones biológicas o estadísticas difíciles, **revisión adversarial**, revisión independiente de cierre | Opus 5 | alto |
| Solo cuando Opus 5 en alto se ha quedado corto **de forma demostrada** | Fable 5.1 | alto |

**Sobre Fable 5.1.** Es el modelo más capaz ampliamente disponible, para el
razonamiento más exigente y el trabajo agéntico de horizonte largo. **Cuesta el
doble que Opus 5** —10 $/50 $ por millón de tokens frente a 5 $/25 $—, así que
no es una mejora gratuita: se justifica por un fracaso concreto de Opus, no por
precaución. Particularidades: el pensamiento está siempre activo y la
profundidad se controla con el esfuerzo; el uso forzado de herramientas
(`tool_choice` `any`/`tool`) devuelve error.

`xhigh`, `max` o equivalentes solo si el usuario lo pide o si una evaluación
reproducible demuestra que hacen falta. No se eleva por precaución genérica.

**Escalado.** Subir de Sonnet a Opus solo si se cumple alguna condición:
contradicción entre especificación, arquitectura y código; invariante científica,
corrupción o pérdida de datos en juego; una corrección grave ha fallado o
reaparecido tras una revisión; el plan exige revisión independiente de cierre; o
el usuario prioriza expresamente la calidad sobre la cuota.

Resuelto el punto crítico, **volver a bajar** para correcciones mecánicas,
documentación, puertas y seguimiento.

## 2. Contexto y concurrencia

- No heredar el chat entero. Enviar una tarea autocontenida con ficheros,
  contratos, límites y resultado esperado.
- Historial completo solo si la tarea depende de decisiones conversacionales que
  no estén registradas en el repositorio.
- **No más de dos subagentes Opus activos a la vez.** Confirmado por el
  responsable el 22 de septiembre de 2026; hasta entonces el documento decía
  uno, heredado literalmente de la política de Codex.
- **Y nunca dos sobre el mismo estado mutable.** No es lo mismo que el límite
  anterior y protege de otra cosa: dos agentes que trabajen a la vez sobre el
  árbol de trabajo se pisan. El 21 de septiembre un revisor y un implementador
  compartieron el directorio temporal y **por poco cruzan sus arneses de
  mutación**, lo que habría hecho perseguir un fantasma. La forma de evitarlo
  que funcionó: el que solo lee trabaja sobre una copia extraída de un commit
  con `git archive`, con su propio directorio de compilación, y el encargo se
  lo dice explícitamente.
- No lanzar dos agentes para la misma pregunta. Las revisiones independientes se
  ejecutan en secuencia y sobre el estado ya corregido.
- No crear un agente solo para esperar, consultar estado o ejecutar un comando.
- **No reutilizar un implementador para una tarea ajena**: su contexto acumulado
  se cobra entero en cada reanudación. Para un asunto nuevo, agente nuevo.

## 3. Registro mínimo de cada encargo

Todo encargo a un subagente debe llevar:

- objetivo y propiedad exactos;
- modelo y esfuerzo elegidos, y por qué;
- contexto heredado mínimo;
- ficheros o módulos en alcance;
- prohibición de revertir cambios ajenos;
- puertas, pausas y acciones externas permitidas;
- condición de escalado o de **detención** —qué debe hacer que se pare e informe
  en vez de improvisar.

## 4. Ciclo de trabajo por hito

1. **Implementador** (Opus 5, esfuerzo alto).
2. **Revisor adversarial** independiente, que no modifica el repositorio y
   **reproduce cada hallazgo ejecutando**. Si no puede reproducir algo, lo
   declara como sospecha; si no puede ejecutar una puerta, lo dice en vez de
   darla por buena.
3. Los hallazgos vuelven **al mismo implementador**, que los aplica él.
4. El agente raíz verifica **por su cuenta** las puertas y la discriminación de
   las regresiones, mutando y restaurando.
5. Commit y push **solo con autorización explícita**.

El implementador puede rechazar un arreglo propuesto por el revisor si demuestra
que no basta; ha ocurrido y tenía razón. Lo que no puede es rechazarlo sin
demostrarlo.

### Regla de corte: qué reabre una ronda

Decidida por el responsable el 12 de octubre de 2026, tras varias tandas de
cobertura que se alargaron a dos y tres rondas por hallazgos menores.

- **Solo un hallazgo medio o grave reabre ronda**: un defecto real del
  programa, o una prueba que consagra algo contrario al contrato o que no puede
  fallar por lo que dice cubrir.
- **Los menores no reabren ronda**: comentarios o textos del registro
  inexactos, límites que otra comprobación ya caza, equivalentes sin anotar.
  El agente raíz los aplica al hacer el commit si son triviales, o se anotan
  como deuda.
- **Una sola revisión adversarial por tanda o hito.** Si la ronda siguiente solo
  corrige hallazgos menores, basta la verificación del agente raíz.
- **Un objetivo sin final natural se cierra con un criterio escrito**, no cuando
  ya no queda nada que encontrar. Para DT-001: se da por cerrada al terminar la
  tanda 9 del juez de los recursos publicados, con lo que quede anotado como
  inalcanzable.

### El defecto recurrente número uno

**Pruebas que no pueden fallar por el defecto que dicen cubrir.** Apareció diez
veces durante F5, y después dos guardas del oráculo resultaron estar
descubiertas por la batería entera. Sus formas: aserción implicada por una
condición anterior; fixture degenerado; oráculo circular; identidad escrita a
mano junto al bucle en vez de contabilizada dentro; mutación que no muta.

De ahí la disciplina, exigida en cada encargo y verificada de nuevo por el
agente raíz:

1. Toda regresión se verifica **fallando** contra la variante defectuosa.
2. Se comprueba **que la mutación mutó** (SHA antes y después).
3. La copia de seguridad va **fuera de git**. Nunca `git checkout --`: ya
   descartó una corrección sin commitear.
4. Si una aserción no puede ejercitarse, se dice **en su docstring**, sin
   adornos, en vez de fingir que discrimina.

Un fixture ordinario a menudo **no** ejercita la guarda, porque otra
comprobación caza antes la discrepancia. Cuando pase, hace falta un fixture
**saturado**, que haga indistinguible todo lo demás.

### Prohibiciones permanentes

- **Bucles de sondeo** (`until ...; do :; done`). Consumen un núcleo, falsean
  cualquier medida de tiempo o CPU y rebotan notificaciones obsoletas. Ha
  ocurrido dos veces; una falseó un rendimiento en un 66 %. **También con
  `sleep` dentro** (`until <condición>; do sleep 30; done`): no consume CPU,
  pero si la condición nunca se cumple espera para siempre. El 12 de octubre de
  2026 uno esperó 13 horas a un fichero que se quedó en 38 de 40 líneas y dejó
  parado al agente que lo esperaba. Para esperar, se ejecuta el trabajo y se
  aguarda su notificación.
- **Commits o push sin autorización explícita.** No hay permiso permanente.
- **Etiquetar una puerta con un comando que no se ejecutó.**

---

## 5. Cómo escribir en el chat

El usuario ha pedido expresamente el lenguaje **más claro posible**. A veces le
resulta difícil seguir lo que se le cuenta, y casi siempre la causa es la misma:
nombrar las cosas por su código.

- **No nombrar nada solo por su código, su número o su identificador.** Ni
  «DTC-089», ni «F5.3», ni «GPS-RES-101», ni un hash de commit a secas, ni el
  nombre de una función o una prueba sin más.
- Si el código hace falta —y a veces hace falta, porque es lo que permite
  buscarlo—, **va siempre acompañado de qué es**: «DTC-089, el contrato que fija
  cuántos recursos puede consumir una ejecución».
- Lo mismo con ficheros y funciones: decir **qué hacen**, no solo cómo se
  llaman.
- Al recuperar un asunto de antes, recordar en una frase de qué iba, sin dar por
  supuesto que el código lo evoca.
- Si el usuario pide que algo se explique más claro, reexplicarlo entero en
  llano, sin hacerle sentir que está repitiendo la pregunta.

Esto vale para el chat. En los documentos y en el código, los identificadores
normativos se usan como siempre: allí son la referencia exacta.

---

## 6. Este proyecto

Simulador de poblaciones genómicas en Rust 2024, con `unsafe_code = "forbid"`
declarado en los dos raíces de crate (`src/lib.rs` y `src/main.rs`).

- **Castellano** en documentos, planes, informes, comentarios y conversación.
  **Inglés** en el código y en los nombres de las pruebas.
- **Contra la sobreingeniería.** «¿No sería más fácil...?» es una corrección, no
  una duda.

### Dónde está la autoridad

| Qué | Dónde |
|---|---|
| Fases, hitos y reglas de avance | `design/hoja-de-ruta-0.1.md` |
| Decisiones técnicas normativas (DTC-xxx) | `design/decisiones-tecnicas-0.1.md` |
| Requisitos funcionales | `design/especificaciones-0.1.md` |
| Arquitectura | `design/arquitectura-0.1.md` |
| Plan de la fase en curso | `design/plan-implementacion-fN-0.1.md` |
| Ampliaciones aplazadas (DTF-xxx) | `design/registro-ampliaciones-futuras-0.1.md` |
| **Deudas técnicas del código actual (DT-xxx)** | `design/registro-deudas-tecnicas-0.1.md` |
| Skills del proyecto (diez, en castellano) | `.claude/skills/*/SKILL.md` |

Ante un conflicto entre el código y un documento normativo, **manda el
documento**, y la contradicción se resuelve antes de seguir. Una cifra o una
propiedad afirmada en un documento normativo **se comprueba**: ya vivió meses
una premisa falsa («el literal `f64` mínimo de ida y vuelta mide 24 bytes»,
cuando `f64::to_string` puede dar 326).

### Skills disponibles

`implementar-funcionalidad-rust`, `revisar-codigo-rust`,
`revisar-rendimiento-rust`, `aplicar-correcciones-rust`,
`auditar-clonaciones-rust`, `ejecutar-plan-rust`, `planificar-implementaciones`,
`redactar-arquitectura`, `redactar-especificaciones`, `redactar-texto-tecnico`.

Viven en `.claude/skills/`, que es donde Claude Code las descubre: se ofrecen
solas y se invocan por su nombre. Hasta el 19 de septiembre de 2026 estuvieron
en `.agents/skills/`, la ruta de Codex, donde **no se cargaban**; había que
abrirlas a mano como ficheros y citar su ruta en cada encargo. Si alguna vez
dejan de aparecer en la lista de skills, ese es el primer sitio donde mirar.

## 7. Puertas

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --workspace --all-targets --locked
cargo test --workspace --doc --locked
cargo build --release --locked --target x86_64-unknown-linux-gnu
./scripts/check-independent-validators.sh <directorio-nuevo>   # 5 PASS; falla si ya existe
./scripts/check-guard-coverage.sh <directorio-nuevo>           # cobertura de los jueces; falla si ya existe
MARIMO_SKIP_UPDATE_CHECK=1 uv run --frozen marimo check --strict ui/gps_app.py
uv run --frozen python -m unittest discover -s tests_ui
```

**Dependencias externas de las puertas.** Este proyecto tiene un problema
recurrente con ellas: sus pruebas de Rust **lanzan procesos externos**
—`python3` para el validador de esquemas, `bcftools` y `tabix` como jueces del
formato VCF— y el resultado de la batería depende de qué encuentre el camino de
búsqueda. Ha mordido tres veces. Antes de dar por buena una cifra de pruebas,
comprobar que se midió con el entorno correcto.

- **Hay que instalar también las dependencias de desarrollo**, no solo las de
  uso: `uv sync --all-groups` o equivalente. Sin ellas faltan 41 pruebas,
  porque el validador de esquemas se lanza como proceso externo y no encuentra
  su módulo. La versión está **fijada a propósito** en la 4.10.3: las modernas
  retiran la pieza que el validador usa.
- **La puerta de los validadores independientes se ejecuta con el `python3`
  que encuentre el camino de búsqueda**, que es el del entorno virtual del
  proyecto (3.12.7) si está activado y el del sistema (3.12.3) si no. Desde el
  22 de septiembre de 2026 **las dos están aprobadas**, revalidadas sobre el
  caso corrupto con resultado idéntico; antes solo lo estaba la del sistema y
  la puerta fallaba según desde dónde se lanzara. Y el script **crea el
  directorio de evidencia antes de comprobar las versiones**, así que un fallo
  deja el directorio creado y obliga a usar un nombre nuevo.
- **`bcftools` y `tabix` deben ser exactamente 1.24**, que es lo que DTC-066
  —la decisión que fija la puerta externa del formato VCF— exige. Se instalaron
  compiladas en `/usr/local` el 20 de septiembre de 2026, que va por delante de
  `/usr/bin` sin tocar el `.bashrc` ni desinstalar las del sistema, de las que
  depende `vcftools`. Con ellas, la prueba que comparaba la normalización
  contra ese juez externo dejó de estar ignorada.

**La puerta de cobertura de los jueces.** `check-guard-coverage.sh` mide, en
los módulos que verifican que los datos publicados son correctos, cuántas de
sus comprobaciones internas ejercita de verdad alguna prueba. Existe desde el
23 de septiembre de 2026, cuando una medición encontró que **dos tercios no
las ejercitaba nadie**: comprobaciones escritas, ejecutadas y jamás probadas,
de modo que si una estuviera mal escrita nadie lo notaría y el simulador daría
por bueno un fichero incorrecto.

Vigila **catorce módulos**. Los cinco primeros: el oráculo de aceptación de F5
(`src/f5_oracle.rs`), el oráculo semántico del VCF
(`src/vcf_semantic_validation.rs`), la auditoría de identidad y publicación
(`src/f5_closure.rs`), el juez de los recursos publicados, de la configuración
efectiva y del manifiesto general (`src/results.rs`) y el del manifiesto de
evento (`src/event_manifest.rs`); los dos últimos entraron el 29 de septiembre
de 2026. Los nueve restantes entraron el 5 de octubre de 2026 y son los que
juzgan lo publicado o lo leído: el juez semántico de la verdad
(`src/truth_validation.rs`) y su contrato (`src/truth_schema.rs`), el del FASTQ
final (`src/fastq_validation.rs`), la validación de la matriz de SNV
(`src/population_snvs.rs`), la verificación del paquete F5
(`src/f5_publication.rs`), la del FASTA suministrado (`src/reference_fasta.rs`),
la del mapa de reducción (`src/reference_reduction_io.rs`), la identidad de
artefactos (`src/integrity.rs`) y la del perfil técnico
(`src/technical_profile.rs`). **La autoridad es el array `judges` de
`scripts/check-guard-coverage.sh`**, que fija cada módulo y los nombres de sus
constructoras de error; la fila de `scripts/guard-coverage-baseline.tsv` es
obligatoria pero no basta: una fila sin entrada en el array se ignora sin
avisar. Para añadir un módulo hay que tocar los dos ficheros, este párrafo y la
lista de `tests/guard_constructors.rs`, que impide que un rechazo se salte la
constructora. En los nueve últimos, cada rechazo pasa por `rejected(..)`
(`src/guard_trace.rs`) o por una constructora propia con `#[track_caller]`; en
`ok_or` hay que escribir `ok_or_else(|| rejected(..))`, porque
`ok_or(rejected(..))` anotaría siempre. Las comprobaciones de estos módulos
subieron de golpe al entrar (565 nuevas, 446 sin ejercitar): se hizo visible lo
que no se veía, no es una regresión.

Compara contra `scripts/guard-coverage-baseline.tsv`, que declara por juez
cuántas comprobaciones tiene y cuántas pueden quedar sin ejercitar. **Se
comprueban las dos columnas**: sin la primera, borrar una comprobación bajaría
las dos cifras a la vez y la puerta lo tomaría por una mejora.

**Lo que esta puerta mide y lo que no**, porque confundirlo sería peor que no
tenerla: mide si una comprobación **llega a dispararse**, no si una prueba
**fallaría** en caso de que estuviera mal escrita. Es un suelo, no un techo. Lo
que demuestra que una comprobación está probada de verdad sigue siendo
**mutarla y ver morir la prueba**. Ver DT-001 en el registro de deudas.

**Las pruebas de Python viven en `tests_ui/`**, que es lo único que recorre la
puerta canónica (`unittest`). `pytest` no está instalado ni declarado: no se
usa. Un fichero de pruebas Python fuera de `tests_ui/` **no lo ejecuta nadie**;
así pasó con `tests/test_process_adapter.py` del 25/08 al 06/10/2026 (DT-004,
cerrada) y sigue pasando con `tests/test_population_snvs_editor.py` (DT-037).

Fallos preexistentes y ajenos que **no** hay que perseguir: dos en `tests_ui`
(F2.5 VCF UI y F4.1 evento técnico) y nueve pruebas ignoradas históricas de
F3/F4 —eran doce hasta que la de normalización con `bcftools` se encendió el
20 de septiembre de 2026, y once hasta que el 06/10/2026 se borraron dos que
estaban vaciadas (DT-002)—.

## 8. Entorno

Disco WSL. **El 21 de septiembre de 2026 el volumen se amplió**: hoy son 393 GB
con 143 libres, al 62 %. Antes eran 256 GiB al 96 %, y esa cifra circuló en
informes hasta que una revisión la encontró caducada — conviene volver a
medirla antes de citarla. Sigue valiendo lo demás: el volumen **no crece solo**
más allá de su techo, y liberar dentro no devuelve espacio a Windows sin
compactar el `.vhdx`. `target/debug` ocupa 16 GB; `cargo clean` solo adelanta
trabajo, porque se rellena en la siguiente compilación.
