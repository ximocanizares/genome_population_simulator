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
- **No más de un subagente Opus activo a la vez.**
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
  ocurrido dos veces; una falseó un rendimiento en un 66 %.
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
| Skills del proyecto (diez, en castellano) | `.agents/skills/*/SKILL.md` |

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

Cada encargo debe citar las que apliquen.

## 7. Puertas

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --workspace --all-targets --locked
cargo test --workspace --doc --locked
cargo build --release --locked --target x86_64-unknown-linux-gnu
./scripts/check-independent-validators.sh <directorio-nuevo>   # 5 PASS; falla si ya existe
MARIMO_SKIP_UPDATE_CHECK=1 uv run --frozen marimo check --strict ui/gps_app.py
uv run --frozen python -m unittest discover -s tests_ui
```

**Trampa de `pytest`:** no está instalado ni declarado en `pyproject.toml`. La
puerta canónica es `unittest`, y bajo ella `tests/test_process_adapter.py`
aborta con `ModuleNotFoundError` desde el commit `6139fe9` (25/08/2026): lleva
sin ejecutarse desde entonces. Solo corre con `uv run --with pytest`. Deuda
registrada en el plan de F5, sin destino asignado.

Fallos preexistentes y ajenos que **no** hay que perseguir: dos en `tests_ui`
(F2.5 VCF UI y F4.1 evento técnico) y doce pruebas ignoradas históricas de
F3/F4.

## 8. Entorno

Disco WSL con techo duro de **256 GiB**, al 96 %. No crece solo, y liberar
dentro no devuelve espacio a Windows sin compactar el `.vhdx`. `target/debug`
ocupa 16 GB; `cargo clean` solo adelanta trabajo, porque se rellena en la
siguiente compilación.
