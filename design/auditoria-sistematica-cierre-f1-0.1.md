# Auditoría sistemática de cierre de F1 0.1

## 1. Identificación, alcance y veredicto

- **Fecha de revisión:** 15 de agosto de 2026.
- **Revisión congelada:** `66ad8b2f96b7b32b7ccc78a4845d26de6d0f1967`.
- **Rama remota contrastada:** `origin/main` apuntaba al mismo commit.
- **Árbol inicial:** limpio.
- **Cambios concurrentes excluidos:** durante la redacción aparecieron cambios
  no confirmados ajenos a esta auditoría en `design/arquitectura-0.1.md`,
  `design/decisiones-tecnicas-0.1.md` y `design/especificaciones-0.1.md`. No se
  revisaron como parte del snapshot ni se modificaron desde esta tarea.
- **Base del cambio revisado:** `8dd676a49e89379d9c5ea408f61be4a9c8d9e143`.
- **Tamaño del cambio:** 55 ficheros, 22.208 inserciones y 2.030 eliminaciones.
- **Ámbito funcional:** F1.1–F1.6, con atención especial a configuración,
  reglas, `EffectiveConfig`, referencias completas y reducidas, RNG,
  publicación, manifiestos, esquemas, catálogos, pruebas y evidencia.
- **Fuera de alcance:** implementar F2 o corregir código. F2 se consideró solo
  para comprobar que F1 no asumiera autoridad o biología que pertenece a una
  fase posterior.

**Veredicto sobre el alcance actual: BLOQUEAR el cierre de F1.** Se han
reproducido tres defectos graves: un evento declarado en un escenario no puede
usarse donde el propio contrato lo admite; la procedencia de un escenario no
queda ligada a su enlace autoritativo; y un paquete fallido temprano acepta
evidencia de resolución semánticamente falsificada si se recalculan sus hashes.
Los gates verdes no detectan estas tres contradicciones.

Este veredicto no obliga a reparar individualmente los tres defectos. La
evaluación de alcance concluye que dos de los mecanismos afectados y gran parte
del tercero anticipan fases futuras. Si se adopta la frontera F1 mínima de §11,
deben eliminarse o trasladarse primero los caminos que dejan de pertenecer a
F1; un hallazgo desaparece legítimamente cuando desaparece toda la capacidad y
su promesa documental, no cuando solo se suprime su prueba.

No se modificó código, pruebas, planes ni evidencia durante la auditoría. Las
tres reproducciones se crearon y ejecutaron exclusivamente bajo `/tmp`.

## 2. Intención funcional reconstruida

F1 pretende producir un núcleo Rust que admita TOML v1, construya una
configuración efectiva inmutable, derive decisiones aleatorias reproducibles,
prepare una referencia simulada o suministrada y publique un paquete coherente
sin sobrescribir un destino. F1.6 añade una reducción reproducible de un FASTA.
No pretende crear individuos, poblaciones, variantes, VCF, lecturas ni cruces.

La implementación revisada va bastante más lejos que ese núcleo: materializa
el retículo de ocho familias de parámetros sobre objetos aún inexistentes,
publica evidencia de resolución separada, crea manifiestos para eventos
técnicos sin consumidor, revalida un grafo documental y físico exhaustivo y
defiende parte de la publicación frente a un atacante local activo. Esta
extensión explica tanto su cobertura amplia como los tres huecos de integración
encontrados.

## 3. Modelo de amenaza y criterio A/B/C/D

Las categorías usadas en este informe son:

- **A — esencial para F1 mínima:** sin el mecanismo no se puede ejecutar F1 de
  forma científicamente reproducible o detectar corrupción accidental normal.
- **B — necesario, pero propiedad de una fase posterior:** el contrato tiene
  sentido cuando exista su consumidor real y debe probarse allí.
- **C — hardening opcional:** defensa de CI/publicación condicionada al modelo de
  amenaza, especialmente frente a un atacante local activo.
- **D — redundante o simplificable:** duplicación de autoridad, derivado usado
  para probarse a sí mismo o producto sin consumidor que puede eliminarse sin
  perder los requisitos funcionales actuales.

Se distinguen tres objetivos que no deben mezclarse:

1. **Reproducibilidad científica:** misma versión, semilla, configuración,
   entradas, algoritmos y orden lógico producen la misma verdad lógica.
2. **Integridad accidental:** una escritura truncada, copia distinta, fallo de
   proceso, destino existente o corrupción de bytes se detecta mediante estado,
   tamaño, hash y publicación atómica.
3. **Atacante local activo:** otro proceso con acceso al árbol intenta sustituir
   rutas o bytes mediante symlinks, hardlinks, alias de inodo o carreras TOCTOU.
   Defenderse de él requiere un contrato de seguridad explícito, descriptores y
   aislamiento; no se deduce de la reproducibilidad científica.

## 4. Estado de gates y herramientas

| Comprobación | Resultado en `66ad8b2` | Observación |
|---|---:|---|
| `cargo fmt --all -- --check` | pasa | Sin diferencias de formato. |
| `cargo clippy --workspace --all-targets --all-features --locked -- -D warnings` | pasa | Cero warnings. |
| `cargo test --workspace --all-targets --all-features --locked` | pasa | 196 pruebas; cero fallos. |
| `cargo test --workspace --doc --all-features --locked` | pasa | Dos `compile_fail`; cero fallos. |
| `cargo doc --workspace --no-deps --all-features --locked` | pasa | Documentación generada. |
| `cargo build --release --workspace --all-features --locked --target x86_64-unknown-linux-gnu` | pasa | Target fijado por F1. |
| `cargo metadata --locked --format-version 1 --no-deps` | pasa | Rust 1.97.1, edición 2024, cero features propias. |
| `cargo tree --locked -e features` | pasa | Árbol inspeccionado; incluye `jsonschema` y sus transitivas. |
| `cargo audit --no-fetch` | pasa | 1.216 avisos locales cargados; cero vulnerabilidades informadas. La actualización de índices no fue posible desde el entorno de solo lectura. |
| `git diff --check 8dd676a..HEAD` | pasa | Sin errores de espacios del cambio. |
| validadores negativos FASTA, SHA-256, `.tar.gz`, JSON Schema | pasan | Las cuatro corrupciones deliberadas fueron rechazadas. |
| validador OCI | **no ejecutado** | Cliente Docker 29.6.1 presente, pero acceso a `/var/run/docker.sock` denegado. No se declara superado. |

Entorno observado: Ubuntu 24.04.4 sobre WSL2, kernel
`6.18.33.2-microsoft-standard-WSL2`, `x86_64`, Rust/Cargo 1.97.1. La evidencia
histórica OCI del repositorio no sustituye la ejecución de esta revisión.

## 5. Inventario de invariantes y cobertura de autoridad

| Dato o invariante | Fuente autoritativa | Derivados principales | Validador | Oráculo/prueba independiente | Clase |
|---|---|---|---|---|---|
| Clase y versión TOML | bytes del TOML admitido | `DeclaredConfig` | frontera de cabecera antes del parseo semántico | fixtures de clase/versión y spans del parser | A |
| Configuración general y escenarios | TOML general y escenarios enlazados | modelo declarado y copias del paquete | parser + validación cruzada | reordenación semántica y copia byte a byte; falta el vínculo de ruta del hallazgo H2 | A para general; B para semántica futura de escenarios |
| Valores efectivos F1 | TOML + reglas admitidas + semilla | `EffectiveConfig`, JSON efectivo | resolvedor y reconstrucción desde copias | expectativas literales de reglas y permutaciones | A para valores consumidos; B para ejes futuros |
| Dominio/solicitudes de resolución | `resolution-evidence.json` creada por el mismo productor | valores efectivos y digest de ejecución | esquema + parser; se cruza solo si existe `effective-config` | corpus correlacionado; queda autofirmada en fallo temprano (H3) | D en F1 mínima; B/A cuando una fase aporte dominio externo real |
| Regla ganadora | conjunto de reglas declarado + contexto probado | valor efectivo, selectores, procedencia | máxima única del orden parcial | retículo manual, reglas permutadas y conflictos | A solo para reglas activas; B para matriz futura |
| FASTA simulado | secuencias/composición TOML + subsemilla de referencia | bases y FASTA canónico | dominio DNA y escritor | vectores, conteo y `samtools faidx` | A |
| FASTA suministrado | bytes del fichero de entrada | copia autoritativa y modelo de referencia | parser interno, catálogo y relectura de copia | `cmp`, SHA-256 y `samtools faidx` | A |
| Reducción de referencia | FASTA copiado + parámetros exactos + subsemilla | plan, FASTA reducido y mapa | preflight, materialización y relector del mapa | enumeración exhaustiva pequeña y reconstructor ajeno | B, trasladable a F2 |
| Porcentaje exacto | literal TOML decimal | número de bases retenidas | parser decimal exacto y `floor` único | casos de frontera y precisión arbitraria | B con reducción |
| Coordenadas | referencia, secuencia y convenio tipado | intervalos internos y mapa | tipos base-1/base-0 e intervalos | vectores de primera/última base y límites | B en F2; ampliación estructural en F5 |
| Semilla raíz | literal TOML o CSPRNG inyectada antes de resolver | `RootSeed` efectivo | parser de 32 bytes | vectores literales y fuente de entropía de prueba | A |
| Subsemillas | raíz + ámbito/versionado + componentes canónicos | ChaCha12 local | catálogo y codificador | vectores BLAKE3 literales y aislamiento | A para referencia; B para consumidores futuros |
| IDs automáticos | clase + causa + campos canónicos | IDs tipados | registro de identidad | vectores de 128 bits/Base32 y colisiones | A para clases F1; B para individuo/población/eventos futuros |
| Copia de configuración | bytes admitidos | `inputs/config/*` | tamaño/hash y relectura | comparación byte a byte | A |
| `effective-config.json` | modelo efectivo | documento JSON | esquema + semántica | reconstrucción desde entradas y evidencia | A, simplificando su grafo |
| Manifiesto general | estado de ejecución + catálogo cerrado | `manifest.json` | esquema, grafo y relectura física | oráculo semántico separado y mutaciones | A en forma mínima |
| Manifiesto de evento | evento explícito sin consumidor F1 | `events/*/manifest.json` | esquema/grafo | matriz 0/1/2 eventos | B para F3/F4; D mientras esté vacío |
| Catálogo de artefactos | escrituras cerradas del workspace | enlaces dedicados y unidades | identidad/clase/ruta/estado/hash | mutaciones por campo | A en forma simple; D si triplica el mismo enlace |
| Estado complete/failed | máquina de estados de ejecución/unidad | visibilidad y exportabilidad | transiciones y cierre fallido | fallos inyectados por etapa | A |
| Publicación | provisional hermano + `INCOMPLETE` + rename no-clobber | destino visible | validación final y `renameat2` | carrera de dos procesos/SIGKILL | A para atomicidad accidental |
| Symlink/hardlink/inodo/TOCTOU | árbol físico justo antes de publicar | autorización de publicación | `NOFOLLOW`, `nlink`, inventario e inodos | mutaciones de enlaces y hooks | C |
| Normalizaciones numéricas | valores TOML + orden canónico | frecuencias/probabilidades efectivas | finitud, dominio, Neumaier, tolerancia | bits y sumas literales | B hasta el primer consumidor biológico; A solo para un campo F1 real |
| Diagnóstico | error estructurado único | texto humano o JSON | catálogo código/categoría/contexto | equivalencia de renderizados | A |
| Plataforma/dependencias | toolchain, lock y detección de entorno | campos de manifiesto/evidencia | metadata y scripts | CI nativo/WSL2/OCI | C/F9 salvo versión necesaria para reproducción |

## 6. Grafo de autoridad

El grafo pretendido para una ejecución completa es:

```text
TOML general ──enlaza──> TOML escenario(s)
      │                         │
      └────────┬────────────────┘
               v
     DeclaredConfig + ResolutionDomain/requests
               │
               ├──> resolution-evidence.json
               └──> EffectiveConfig ──> effective-config.json
                            │
FASTA origen + semilla ─────┼──> referencia(s) / mapa(s)
                            │
                            v
        catálogo de artefactos + unidades + eventos
                            │
                            v
                      manifest.json
```

La revalidación relee hashes y luego intenta invertir parte del grafo:

```text
manifest ─> copias TOML + resolution-evidence + effective-config
                                      │
                                      v
                           EffectiveConfig reconstruida
```

Resultados de la auditoría del grafo:

- La separación TOML declarado → efectivo es correcta y evita reinterpretación
  mutable.
- FASTA suministrado se vuelve autoritativo desde la copia y no desde la ruta
  externa después de publicarse.
- Los hashes verifican bytes, no su significado, y se contrastan antes de la
  publicación.
- `resolution-evidence` no tiene una fuente externa en el CLI F1: la genera el
  mismo objeto que genera `effective-config`. Por tanto no rompe por sí sola la
  autofirma; H3 demuestra que, sin el derivado, queda solo autenticada por su
  propio hash actualizado.
- La ruta de un escenario aparece en el TOML general, en la procedencia de la
  entrada y en la ruta interna de copia. Solo las dos últimas se usan para
  reparsear; falta cruzar la primera con la segunda (H2).
- Eventos se integran para identidad y manifiestos desde general+escenarios,
  pero el catálogo del resolvedor y las solicitudes usa solo el general (H1).
- El mismo artefacto aparece como entrada especializada, fila de `artifacts` y
  en ocasiones unidad/enlace adicional. La igualdad se comprueba ampliamente,
  pero la multiplicidad aumenta combinaciones sin crear otra autoridad.

## 7. Matriz exhaustiva de reglas y objetos

`✓` significa permitido por el catálogo actual; `—`, prohibido. `P/S` reúne
población y subpoblación; `F/G`, fase y generación; `S/R/L`, secuencia, región y
locus. La regla global/base existe sin selector y no cuenta como eje.

| Parámetro | Dominio | Individuo | P/S | F/G | Secuencia | Región | Locus | Grupo | Bloque | Evento técnico |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `ploidy` | entero 1..10 | — | — | — | ✓ | — | — | — | — | — |
| `variant_rate` | real finito ≥0 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | — |
| `selfing_rate` | probabilidad [0,1] | — | ✓ | ✓ | — | — | — | — | — | — |
| `recombination_rate` | real finito ≥0 | — | — | — | ✓ | — | — | ✓ | — | — |
| `coverage` | real finito ≥0 | ✓ | — | — | ✓ | ✓ | — | — | — | ✓ |
| `error_rate` | probabilidad [0,1] | ✓ | — | — | ✓ | ✓ | — | — | — | ✓ |
| `contamination_rate` | probabilidad [0,1] | ✓ | — | — | — | — | — | — | — | ✓ |
| `output_enabled` | booleano | ✓ | — | — | — | — | — | — | — | ✓ |

Matriz de jerarquías y objetos:

| Familia | Prueba de pertenencia exigida | Contexto completo | Incompleto/incoherente | Estado |
|---|---|---|---|---|
| global | `object_id=global`, salvo ploidía | sin selectores | cualquier objeto extranjero | cubierto |
| individuo | ID directo declarado | individuo = objeto | ID inexistente o categoría distinta | cubierto para IDs generales; consumidor futuro B |
| población/subpoblación | `ResolutionDomain.population_memberships` | subpoblación pertenece exactamente a población | lista plana, padre ausente o padre distinto | cubierto con proveedor sintético; B |
| fase/generación | intervalo del dominio | generación dentro de la fase | fase ausente, fuera de intervalo o fase duplicada global | cubierto sintéticamente; B |
| referencia/secuencia | catálogo de referencia | secuencia pertenece a referencia | referencia/secuencia cruzada | cubierto y A |
| región | referencia+secuencia+intervalo cerrado válido | posición dentro del intervalo | vacío, fuera de límites o solape prohibido | cubierto; A solo para reducción, B general |
| locus | dominio con referencia+secuencia+posición | los cuatro coinciden | locus huérfano o coordenada distinta | cubierto sintéticamente; B |
| grupo/bloque | ID declarado | objeto coincide | ID ausente/categoría incorrecta | cubierto sintéticamente; B |
| evento técnico | ID explícito declarado | evento coincide | ausente/categoría incorrecta | **incompleto para eventos de escenario (H1)**; B |

La dominancia implementada es un orden parcial por refinamiento dentro de cada
eje. Una regla combinada domina sus bases; reglas específicas en ejes distintos
son incomparables; dos máximas incomparables fallan incluso si tienen el mismo
valor. No existe desempate por orden físico. Se cubren base obligatoria,
combinadas, permutaciones, regiones adyacentes, duplicados exactos y máximas
incomparables. El hallazgo H1 no está en el retículo: ocurre antes, al construir
un catálogo de identidades incompleto.

## 8. Matriz de paquetes y mutaciones

### 8.1. Paquetes recorridos por la suite

| Eje | Casos observados | Resultado de cobertura |
|---|---|---|
| estado | `complete`, `failed`, provisional interrumpido | cubierto |
| escenarios | 0, 1 y varios documentos | positivos cubiertos; vínculo de ruta no cubierto (H2) |
| reducciones | 0, 1 y 2 | complete/failed, fallo en FASTA/mapa/metadatos y unidades posteriores |
| eventos | 0, 1, 2 y varios | enlaces/orden/IDs cubiertos; catálogo de evento de escenario no cubierto (H1) |
| entradas | general único, escenarios 0..N | ausencia, duplicado, extra, reorden y clase cruzada cubiertos; procedencia correlacionada sobrevive (H2) |
| artefactos | configuración, evidencia, efectivo, FASTA, mapa y manifiesto de evento | ruta, tamaño, hash, clase, estado, validación y propietario mutados |
| estados/unidades | ejecución, reducción y evento | prematuro, huérfano, cruzado, complete/failed |

### 8.2. Operadores de mutación ejecutados

| Familia | Mutaciones ejecutadas | Resultado |
|---|---|---|
| hash/tamaño/ruta | un byte, truncado, digest falso, tamaño falso, ruta ajena | rechazadas |
| mutación correlacionada | contenido + hash recalculado + enlace actualizado | rechazadas en paquetes completos; **sobrevive en evidencia de fallo temprano (H3)** |
| identidad/clase | ID existente de clase incorrecta, ID cruzado, propietario reutilizado, huérfano | rechazadas |
| cardinalidad | añadir, eliminar y duplicar entradas, artefactos, unidades, referencias y eventos | rechazadas en el corpus cubierto |
| orden | TOML, reglas, solicitudes, eventos, secuencias y artefactos reordenados | se conserva semántica o se rechaza cuando el orden es normativo |
| sistema de ficheros | fichero extra, symlink intermedio/final, hardlink externo, alias de inodo | rechazadas antes del rename |
| TOCTOU | cambio de evidencia/FASTA entre preparación y validación final | rechazado por hooks de prueba; no demuestra inmunidad ante toda carrera activa |
| gramática | TOML duplicado/desconocido, JSON truncado/campos extra/clase/versión, decimal con exponentes/signos/separadores | rechazadas según contrato |
| números/coordenadas | cero, máximos, overflow, no finitos, tolerancias, límites base-1/base-0, intervalos vacíos/solapados | rechazados o aceptados en el borde literal previsto |
| UTF-8/bytes | etiquetas Unicode, escapes de diagnóstico, cabeceras FASTA no UTF-8 y CRLF | cubierto sin reinterpretar bases |
| fuentes externas | cambio del FASTA entre validar/copiar, `samtools`, `sha256sum`, JSON Schema Python | cubierto; OCI no ejecutable en esta sesión |
| integración nueva de auditoría | evento solo en escenario; `original_path` correlacionado; evidencia forjada en fallo temprano | tres supervivientes, H1–H3 |

## 9. Independencia de los oráculos y pasadas temáticas

### 9.1. Oráculos

- Los vectores SHA-256, BLAKE3, Base32, RNG y porcentajes tienen expectativas
  literales y no llaman al productor para calcular el esperado.
- La selección de intervalos se compara con enumeraciones exhaustivas pequeñas;
  la reconstrucción de FASTA desde mapa recorre otra ruta.
- JSON Schema se contrasta con una implementación Python Draft 2020-12 además
  del validador Rust.
- El oráculo semántico de manifiestos usa parser/tablas/orden y aritmética
  propios en pruebas, y ejercita válidos aceptados e inválidos rechazados.
- No obstante, muchos documentos positivos siguen siendo producidos por el
  mismo código auditado. La independencia solo es real para la expectativa que
  luego se muta, no para relaciones que nunca se enumeraron. H1 y H2 son dos
  ejemplos.
- En H3, esquema, parser y hashes son independientes en forma, pero ninguno
  dispone de una fuente semántica con la que contrastar el dominio/solicitudes.
  Varias validaciones de la misma afirmación no crean autoridad independiente.

### 9.2. Pasadas completas

| Pasada | Resultado resumido |
|---|---|
| autoridad/procedencia | H2 y H3; copias y FASTA autoritativo correctos |
| reglas/ejes | H1; retículo y matriz aislada correctos |
| referencias/coordenadas | sin defecto funcional nuevo; reducción bien cubierta, pero clasificada B |
| esquemas/grafo documental | estructura robusta; exceso de enlaces duplicados D y hueco semántico H3 |
| publicación/FS/TOCTOU | integridad accidental robusta; defensa activa parcial y clase C |
| diagnósticos/estados/evidencia | renderizados coherentes; fallo temprano no liga evidencia (H3) |
| reproducibilidad/semillas | vectores y aislamiento correctos; ámbitos futuros B |
| dependencias/pruebas | gates verdes; `jsonschema` amplía mucho el árbol; audit local verde; OCI no ejecutado |

### 9.3. Funciones no triviales modificadas

Se agruparon por causa para evitar cientos de observaciones repetidas:

| Grupo funcional | Clase de entrada adicional desafiada | Resultado |
|---|---|---|
| admisión/parsers/builders de configuración | documento reordenado, cruzado, duplicado, escenario externo y evento solo en escenario | H1; resto rechazado/canónico |
| catálogos/identidades/solicitudes | misma cadena en otra clase, objeto extranjero, jerarquía no probada | rechazo correcto salvo la omisión asimétrica H1 |
| matching/dominancia/solapes | combinaciones correlacionadas, máximas iguales incomparables, fronteras máximas | correcto |
| decimal/coordenadas/reducción | precisión mayor que ancho fijo, residual, huecos pequeños, límites y overflow | correcto en corpus exhaustivo pequeño |
| FASTA/generación/copia | byte no UTF-8 en cabecera, `N`, CRLF, cambio concurrente de origen | correcto |
| serializadores/relectores de resultados | campo semántico alterado con hash recalculado y procedencia cruzada | H2 y H3 |
| catálogo/grafo de manifiesto | 0/1/N, clase/propietario cruzados, añadir/quitar/reordenar | correcto donde el grafo tiene fuente |
| workspace/publicador | temporal truncado, destino existente, carrera, extra, symlink y hardlink | correcto bajo los hooks; carreras activas generales quedan C |
| diagnóstico/RNG/integridad | control chars, ámbito desconocido, consumidor reordenado, chunks de hash | correcto |

## 10. Hallazgos priorizados

### H1 — Grave — confianza alta

- **Ubicación:** `src/application.rs:683`, `src/application.rs:1258` y
  `src/application.rs:1475`.
- **Contrato afectado:** el TOML general y sus escenarios forman la unidad de
  configuración; F1 admite cero/uno/varios eventos explícitos, y un selector o
  solicitud `technical_event` debe resolver una identidad declarada.
- **Evidencia:** `build_validated_summary` une correctamente eventos del general
  y de todos los escenarios para `explicit_event_ids`. Sin embargo,
  `validate_declared_rules` añade al `SelectorIdentityCatalog` solo
  `general.events()`, y `validate_resolution_request` vuelve a consultar solo
  esa lista.
- **Reproducción independiente:** en `/tmp` se creó un general de vía evolutiva,
  un escenario con `event-scenario`, una regla base de cobertura y otra regla de
  cobertura seleccionada por ese evento. `load_and_validate_config` devolvió
  `GPS-CONF-013`, `expected_category="technical_event"`,
  `found_category="missing"`.
- **Mecanismo:** tres consumidores construyen el mismo catálogo con conjuntos
  distintos; manifiesto/ID ve el evento, reglas/solicitudes no.
- **Impacto:** una configuración válida según el formato y la propia identidad
  efectiva se rechaza antes de ejecutar. No hay corrupción silenciosa, pero la
  infraestructura general declarada como completa no funciona para uno de sus
  orígenes admitidos.
- **Corrección si se conserva el alcance:** construir una única colección
  validada de eventos general+escenarios y pasarla a identidad, reglas y
  solicitudes; probar ambos orígenes, duplicados cruzados y categoría errónea.
- **Alternativa preferida al simplificar:** retirar de F1 los eventos técnicos y
  sus reglas/manifiestos hasta F3/F4. En ese caso no debe repararse código que se
  elimina.
- **Prueba de regresión:** evento de escenario aceptado por regla y solicitud;
  evento ausente o duplicado general↔escenario rechazado.

### H2 — Grave — confianza alta

- **Ubicación:** `src/application.rs:1512` y `src/results.rs:1480`.
- **Contrato afectado:** la copia de escenario y su procedencia deben
  corresponder al enlace declarado por el TOML general; el paquete debe permitir
  reconstruir qué unidad de configuración se ejecutó.
- **Evidencia:** `validate_route_and_scenarios` cruza cardinalidad, ID y vía, pero
  no la ruta del enlace. Al revalidar un paquete,
  `validate_effective_against_configuration_inputs` usa como origen la ruta
  tomada del propio manifiesto para parsear los bytes copiados, sin compararla
  con `general.scenarios[*].path` resuelta desde la ubicación del general.
- **Reproducción independiente:** se produjo un paquete válido y se cambió solo
  `configuration_inputs[kind=scenario].original_path` de `scenario.toml` a
  `different.toml`. Los bytes, tamaño y SHA-256 de la copia no cambiaron. El
  validador devolvió `Ok(())`.
- **Mecanismo:** dos afirmaciones de procedencia se validan localmente, pero no
  existe la arista que las hace referirse al mismo fichero.
- **Impacto:** el contenido simulado no cambia, pero el paquete afirma una
  procedencia que no fue la declarada. Se pierde trazabilidad reproducible y la
  importación/exportación multidocumento futura podría confiar en una ruta falsa.
- **Corrección si se conserva el alcance:** resolver cada ruta relativa desde
  la fuente general con una política léxica/canónica definida, y exigir
  correspondencia ID↔ruta↔copia tanto al cargar como al revalidar; cubrir alias,
  espacios y complete/failed.
- **Alternativa preferida al simplificar:** en F1, tratar escenarios solo como
  bytes enlazados y posponer semántica multidocumento a F2.0; o admitir solo
  ruta directa. Eliminar entonces las promesas de revalidación completa de esa
  procedencia en F1.
- **Prueba de regresión:** mutar `original_path` sin cambiar bytes debe fallar;
  ruta relativa equivalente solo se acepta si la política aprobada la declara
  equivalente.

### H3 — Grave — confianza alta

- **Ubicación:** `src/results.rs:1332`, especialmente el retorno temprano de
  `src/results.rs:1352`.
- **Contrato afectado:** DTC-077 y REQ-MAN-88 presentan
  `resolution-evidence.json` como fuente autoritativa anterior e independiente,
  también en corpus fallidos; una salida derivada no debe autenticarse a sí
  misma.
- **Evidencia:** el validador comprueba ruta, hash, esquema y sintaxis de la
  evidencia. Si el manifiesto fallido no contiene `effective_config`, retorna
  inmediatamente y nunca contrasta dominio/solicitudes con las copias TOML, un
  catálogo de objetos o el ID de ejecución.
- **Reproducción independiente:** se inyectó el fallo
  `AfterConfigurationCopy`, se añadió a la evidencia una solicitud válida en
  esquema para `ghost-object`, se recalcularon tamaño y SHA-256 tanto en el
  enlace como en `artifacts`, y se revalidó el manifiesto. El resultado fue
  `Ok(())`.
- **Mecanismo:** el hash demuestra identidad con la afirmación modificada, no su
  veracidad. En ausencia del derivado, ninguna fuente externa queda disponible
  y la evidencia se autofirma.
- **Impacto:** un paquete `failed` publicado puede pasar validación mientras
  atribuye al productor un dominio o solicitudes que nunca admitió. Afecta la
  cadena de auditoría, aunque no convierte el paquete en `complete` ni altera
  una simulación terminada.
- **Corrección si se conserva el alcance:** validar evidencia contra entradas
  independientes incluso sin configuración efectiva, o registrar una prueba
  autoritativa externa que permita hacerlo; rederivar también la identidad de
  ejecución cuando sea aplicable.
- **Alternativa preferida al simplificar:** eliminar la evidencia separada de
  F1 mientras dominio/solicitudes no sean entradas reales del CLI. Registrar en
  `effective-config` los valores F1 reconstruibles desde TOML y reintroducir la
  evidencia en la primera fase que reciba un `ResolutionDomain` externo real.
  En ese caso H3 desaparece junto al mecanismo.
- **Prueba de regresión:** paquete fallido antes y después del efectivo, con
  dominio, solicitud, object ID y hash mutados correlacionadamente.

No se encontraron otros defectos bloqueantes o graves tras agrupar las
supervivencias por causa. Las ventanas generales frente a un atacante local no
se elevan a hallazgo porque el proyecto no ha declarado ese modelo de amenaza;
se clasifican C.

## 11. Evaluación de alcance y frontera mínima propuesta

### 11.1. Clasificación por mecanismo

| Mecanismo actual | Clase | Conservar/mover/simplificar |
|---|---:|---|
| Rust 1.97.1, edición 2024, lock, fmt/clippy/test/release | A | conservar en F1 |
| TOML v1, cabecera, unknown/duplicate, `DeclaredConfig`/`EffectiveConfig` | A | conservar; reducir campos a consumidores existentes |
| ocho parámetros × todos los ejes y jerarquías sintéticas | B | mover por familias a F2/F5/F6/F7 y F3/F4 |
| regla base, coincidencia y máxima única para valores usados en F1 | A | conservar núcleo pequeño |
| normalización general de frecuencias/probabilidades | B | mover al primer consumidor F2/F6; conservar solo validación numérica local necesaria |
| raíz RNG, derivación BLAKE3 y ChaCha12 local | A | conservar |
| ámbitos RNG de individuos, población, mutaciones y eventos futuros | B | mover con su consumidor; mantener vectores históricos al incorporarlos |
| IDs tipados de ejecución/config/referencia/reducción/artefacto | A | conservar los realmente emitidos |
| IDs de individuo, población, fase, bloque y evento sin consumidor | B | mover con la fase propietaria |
| FASTA simulado/suministrado, copia y validación | A | conservar |
| reducción, decimal exacto, servicio de coordenadas y mapa | B | mover a F2 como preparación de genomas pequeños; F5 solo extiende a referencias estructurales biológicas |
| `effective-config.json` y manifiesto general mínimo | A | conservar |
| `resolution-evidence.json` sin proveedor externo real | D/B | eliminar de F1; reintroducir cuando F2/F6/F7 entregue dominio/solicitudes reales |
| manifiestos de eventos sin artefactos técnicos | D/B | eliminar de F1; implementar en F3/F4 |
| enlace dedicado + fila `artifacts` + unidad para el mismo objeto | D | usar una representación canónica y referencias por ID |
| estado failed/complete, provisional, marcador y rename no-clobber | A | conservar |
| relectura hash/tamaño de todos los artefactos | A | conservar para corrupción accidental |
| rechazo exhaustivo de symlink/hardlink/inodo y hooks TOCTOU | C | mover a suite de seguridad/publicación F9 o mantener sin hacerlo puerta F1 |
| JSON Schema publicado | A | conservar contrato estructural mínimo |
| validación JSON Schema externa dentro del recorrido normal | C/D | ejecutar en CI/aceptación; producción valida entradas y contratos propios |
| corpus combinatorio completo de grafo físico/documental | C/D | reducir en F1 a pares de frontera; ampliar por fase propietaria y en F9 |
| `samtools`, SHA-256 y reconstructor ajenos | A/C | SHA/reconstructor A; `samtools` como CI C, no dependencia runtime |
| tar/OCI, plataforma nativa, rendimiento y auditoría de distribución | C/B | F9 |

### 11.2. F1 mínima concreta

F1 debería cerrar con este recorrido único:

1. `validate` y `run` admiten un TOML general v1 y, como composición documental
   básica, cero o varios escenarios copiados sin interpretar biología futura.
2. Se valida únicamente la referencia inicial, la semilla y los valores que el
   recorrido F1 consume.
3. Se construye una `EffectiveConfig` inmutable y pequeña.
4. Se deriva una identidad de ejecución desde contenido/configuración y se
   genera o copia un FASTA.
5. Se publican copia(s) TOML, FASTA, configuración efectiva y un manifiesto
   general con estado, tamaño y SHA-256.
6. La publicación usa provisional hermano, `INCOMPLETE` y rename no-clobber.
7. Dos ejecuciones equivalentes se comparan lógicamente; entradas inválidas no
   crean destino completo; fallos no se presentan como éxito.

Código/pruebas/documentos que deben conservarse en F1: admisión y parser TOML,
tipos mínimos de configuración, diagnóstico, raíz RNG y derivación usada,
identidad activa, referencia/FASTA, integridad, workspace/publicación atómica,
manifiesto general mínimo y sus pruebas directas. Deben simplificarse los
cross-links y eliminarse expectativas que repitan el productor.

### 11.3. Traslados sin perder requisitos

- **A F2:** reducción de referencia y Marimo que la configura; coordenadas
  comunes; decimal exacto; reglas de variante/ploidía consumidas; loci,
  individuos y población de origen; normalización de frecuencias; catálogos RNG
  e IDs de esos consumidores. F2.0 posee además la importación/exportación
  multidocumento y debe cerrar el vínculo de rutas de escenarios si lo usa.
- **A F5:** referencias derivadas por variantes estructurales, múltiples
  referencias, mapas estructurales y extensiones de coordenadas. La mera
  reducción para obtener un genoma pequeño no debe confundirse con una variante
  estructural individual.
- **A F3/F4:** evento técnico, reglas de cobertura/error/contaminación y
  manifiesto de evento, cuando existan FASTQ/alineamiento y perfiles reales.
- **A F6/F7:** población/subpoblación, autofecundación, fase/generación y su
  `ResolutionDomain` real.
- **A F9:** exportación autocontenida, OCI/tar, validación de plataforma,
  presupuestos/rendimiento, hardening contra atacante local y matrices
  combinatorias de publicación.

La documentación de F1 debe dejar de acreditar como «implementada» cualquier
fila B trasladada. Las decisiones técnicas pueden conservarse como contrato
futuro, pero su estado y dueño deben cambiar; las pruebas existentes pueden
moverse con el módulo o mantenerse como tests de biblioteca no bloqueantes,
sin usarse para afirmar que una fase biológica inexistente está cerrada.

## 12. Aspectos correctos y huecos no ejecutables

Aspectos especialmente sólidos:

- tipos de coordenadas y decimal exacto con fronteras explícitas;
- selección de intervalos sin solape comparada con un oráculo exhaustivo;
- copia autoritativa de FASTA y conservación de `N` y cabeceras crudas;
- aislamiento RNG y vectores normativos literales;
- orden parcial sin desempates accidentales;
- estados y publicación no-clobber bajo fallos y carreras ordinarias;
- rechazo físico de extras, symlinks y hardlinks antes del rename;
- esquemas con discriminadores, campos cerrados y doble implementación;
- diagnósticos humano/JSON derivados del mismo valor estructurado.

Huecos que esta revisión no declara superados:

- OCI actual, por falta de acceso al daemon;
- Ubuntu 24.04 nativo, porque el entorno actual es WSL2;
- inmunidad completa ante atacante local activo: existe una ventana entre
  validación final y rename y otra tras publicación; cerrarla exigiría diseñar
  la operación completa sobre descriptores/aislamiento;
- escalabilidad, paralelismo y grandes genomas, que corresponden a F9;
- toda verdad biológica F2–F8, que no debe inferirse de pruebas sintéticas del
  resolvedor.

## 13. Conclusión

El árbol `66ad8b2` supera todos los gates Rust ejecutables y contiene trabajo de
alta calidad en fronteras numéricas, referencia, RNG y publicación. Sin
embargo, el alcance actual de F1 no puede aprobarse: sus tres contradicciones de
integración y autoridad están reproducidas y son graves.

La vía recomendada no es añadir otra ronda de parches al grafo completo. Primero
debe aprobarse y aplicar la frontera F1 mínima. Después:

- reparar solo los hallazgos que sobrevivan en mecanismos A;
- trasladar los mecanismos B con sus consumidores y oráculos reales;
- decidir explícitamente si C forma parte del modelo de amenaza de publicación;
- eliminar duplicaciones D y sus afirmaciones documentales.

Tras esa simplificación debe ejecutarse una nueva auditoría limpia sobre el
nuevo diff. Esta revisión es la primera revisión sistemática y no constituye la
segunda revisión limpia exigida para aprobar F1.
