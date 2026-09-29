# Plan de implementación F6 0.1 — referencias y variación complejas

## 1. Estado, autoridad y alcance

- **Fecha:** 14 de septiembre de 2026.
- **Estado:** borrador planificable **con decisiones bloqueantes declaradas**.
  F5 está cerrada con sus evidencias aprobadas, de modo que F6 puede
  planificarse conforme a la sección 7 de la hoja de ruta. Pero **no es
  ejecutable de principio a fin** hasta que se cierren las decisiones de la
  sección 3: a diferencia de F5, que llegó a su plan con DTC-089 como
  autoridad técnica única, F6 no dispone de ninguna decisión que fije sus
  contratos externos. Cada paso que dependa de una de ellas lo declara con la
  marca «bloqueado por». Este plan no cierra ninguna de esas decisiones ni las
  sustituye por una elección propia.
- **Fuentes normativas:** REQ-REF-08–11/19–22, REQ-IND-01–03/06–08/12–20/
  28–29/31–59, REQ-DIR-01–05, REQ-CONF-06/08/15/20–25/28–31/55–65/96/116–123,
  REQ-REC-14–20/24–29, REQ-MUT-11, REQ-VER-01–47, REQ-VAL-27/47/48, REQ-MAN-
  07/09–12/27–31/79, REQ-SAL-11/12/17/23, REQ-EVO-42 y los REQ de datos,
  formato, grupos de lectura y alineamiento asignados a F6.6 en la matriz de
  trazabilidad; ACE-06/07/08/15/16/18/20/30/31/32/58/59/60/61/87/95/96/99/
  100/103/107/109/110/111; DTC-036/037/039/043/045/056/061/062–067/089;
  arquitectura §5.2, §5.3, §6, §7.2, §8, §9 y §11; hoja de ruta §F6 y §7;
  registro de ampliaciones DTF-003/008/011/015/021/033/034/040/048/051/055.
- **Modelo de forma:** el plan de F5 (`plan-implementacion-f5-0.1.md`) y el de
  F1. De ellos toma la estructura, no el contenido.

F6 amplía la verdad de F5 **sin alterar sus fundadores intactos**: las copias
fundadoras, los siete TSV y la proyección `truth.jsonl` v2 que F5 publica no
cambian ni un byte por el hecho de que F6 exista, y una ejecución F6 sin
eventos estructurales publica exactamente lo que F5 publica. F6 añade
referencias estructurales derivadas, variantes estructurales y repeticiones en
tándem sobre individuos, ploidía por secuencia, la primera representación VCF
compleja con complemento, y la separación de conjuntos y matrices por
referencia de mapeo.

**Dos fronteras que no se cruzan.** F6 no hereda panel ni recombinación de
F5, porque no existen: la recombinación fundadora se retiró íntegramente el 31
de agosto de 2026 y no se conservaron formatos ni algoritmos. Meiosis, mosaicos
y recombinación son propiedad exclusiva de F7.1 (DTPD-111). Y ningún hito de
F6 reescribe fundadores: si un paso pareciera exigirlo, es una señal de alarma
que detiene el paso y va al informe, no una licencia.

### 1.1. Cómo leer este plan

- Los pasos heredan la coordenada de la hoja de ruta: F6.1.1, F6.1.2, …; el
  subpaso diferido usa F6.5.a.1, F6.5.a.2, …
- «DP-F6-nn» es una **etiqueta local de este plan** para una decisión
  pendiente. No es un identificador normativo: cerrar una DP exige una fila
  DTC o DTPD en `decisiones-tecnicas-0.1.md`, o una corrección de la
  especificación o de la matriz, según su naturaleza.
- Toda prueba que este plan pide dice **contra qué variante defectuosa debe
  fallar**. Una aserción que no pueda ejercitarse se declara así en su
  docstring, sin adornos (CLAUDE.md §4). Las mutaciones se aplican sobre una
  copia fuera de git, se comprueba por SHA que la mutación mutó y que la
  restauración devolvió el fichero byte a byte.

## 2. Criterios de orden y puertas comunes

1. **Contratos primero, dentro de cada hito.** No se reabre ningún esquema
   publicado por F1–F5: F6 crea versiones nuevas (general, efectivo, manifiesto
   general, verdad) y las completa **dentro de la misma versión** mientras no
   estén publicadas externamente, con la misma justificación que el plan de F5
   §1 dio para `resource_preflight`. Los esquemas heredados se comparan byte a
   byte en cada hito.
2. **Núcleo puro antes que ficheros.** El motor estructural (F6.1.2) se prueba
   sobre secuencias en memoria, sin FASTA ni TOML, antes de derivar una
   referencia o tocar una cadena de segmentos. Lo mismo para el modelo de
   repetición (F6.4.1) y para el modelo escalonado (F6.4.4).
3. **Referencia antes que individuo.** Las cuatro operaciones estructurales se
   construyen una vez y se aplican primero a una referencia (sin fase ni dosis)
   y después a copias haplotípicas. REQ-IND-48 exige que translocaciones sirvan
   en ambos usos con el mismo motor.
4. **Reutilizar antes que duplicar.** Servicio de coordenadas
   (`src/coordinates.rs`, DTC-056), cadena de segmentos y materialización
   (`src/segments.rs`, `src/materialization.rs`, DTC-061), mapa de reducción y
   escritor FASTA canónico (`src/reference_reduction_io.rs`,
   `src/reference_fasta.rs`), `poisson-count-v1` y
   `dynamic-indel-selection-v1` (`src/indel_application.rs`), proyección y
   validación de verdad (`src/truth_projection.rs`, `src/truth_validation.rs`),
   agrupación y oráculo VCF (`src/vcf_grouping.rs`,
   `src/vcf_semantic_validation.rs`), publicación atómica (`src/publication.rs`).
5. **Riesgo silencioso aislado.** Coordenadas tras una inversión o una
   duplicación, orientación, fase de un alelo estructural, dosis con ploidía
   por secuencia, correspondencia origen↔derivada y aislamiento entre
   referencias se contrastan siempre con un reconstructor u oráculo que **no
   comparte código** con el productor (arquitectura §11).
6. **Cerrar un hito antes de apoyar otro.** Cada hito termina con pausa y
   evidencia; F6.2 no empieza sin F6.1 cerrada, F6.5 no empieza sin F6.2 y
   F6.4 cerradas.

Cada paso ejecuta sus pruebas específicas y, antes de su pausa, las puertas de
`CONTRIBUTING.md` y CLAUDE.md §7, restringiendo las de prueba al subconjunto
pertinente; el cierre de cada hito ejecuta la batería completa:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --workspace --all-targets --locked
cargo test --workspace --doc --locked
cargo build --release --locked --target x86_64-unknown-linux-gnu
./scripts/check-independent-validators.sh <directorio-nuevo>
MARIMO_SKIP_UPDATE_CHECK=1 uv run --frozen marimo check --strict ui/gps_app.py
uv run --frozen python -m unittest discover -s tests_ui
```

`pytest` no está instalado ni declarado: la puerta canónica de Python es
`unittest`. Los dos fallos preexistentes de `tests_ui` (F2.5 VCF UI y F4.1
evento técnico) y las doce pruebas ignoradas históricas de F3/F4 no se
persiguen desde F6.

**Regla de cierre Marimo.** Toda opción configurable que F6 implemente
—referencias estructurales, eventos SV, umbral indel/SV, ploidía por secuencia,
loci y parámetros de repetición, cantidades y regiones aleatorias, referencias
de mapeo, formato y mapeador— exige control en la interfaz y equivalencia
UI→TOML→`EffectiveConfig` antes de cerrar su hito (hoja de ruta §7).

**Techo implícito heredado.** DTC-089 fija 65536 bytes inclusivos para el
general y el escenario. Una lista larga de eventos explícitos o de loci de
repetición en el TOML choca con ese techo mucho antes que en F5. F6.1.1 lo
comprueba con un documento legítimo próximo al límite; ampliarlo reabre
DTC-089 y no pertenece a F6.

## 3. Requisitos previos y decisiones bloqueantes

Ninguna de estas decisiones la toma este plan. Cada una dice qué hay que
decidir, qué opciones hay a la vista y qué pasos bloquea.

| DP | Qué hay que decidir | Opciones a la vista | Bloquea |
|---|---|---|---|
| DP-F6-01 | **Contrato TOML y versiones de F6.** Dónde se declaran referencias estructurales, eventos SV explícitos, loci de repetición y parámetros aleatorios: en el general, en un escenario F6 enlazado, o en una extensión de `f5_scenario`. Qué versiones nuevas nacen (general v4, efectivo v6, manifiesto general v6, verdad v3, manifiesto de evento v3) y cómo coexisten la vía directa (general v2) y la fundadora (general v3, que hoy **exige** `f5_scenario`) | (a) general v4 único que admite ambas vías y bloques F6 opcionales; (b) escenario F6 propio enlazado desde el general, al estilo F5; (c) extender `direct_events` del v2 y el escenario F5 por separado | F6.1.1 y todo lo que publica configuración o esquema |
| DP-F6-02 | **Contrato propio de la extensión sobre fundadores F5** (REQ-EVO-42, DTF-040): en qué orden se aplican los eventos F6 respecto al estado F2/F5 ya existente, qué ocurre con un SNP fundador dentro de una deleción, duplicación o inversión F6, y si un SNP o indel previo puede quedar bajo un punto de ruptura posterior | (a) F6 se aplica después de todo el estado previo y un SNP cubierto pasa a `missing_by_deletion`/copiado/invertido con procedencia; (b) F6 rechaza toda colisión con el estado previo; (c) F6 solo se aplica sobre regiones sin variación previa | F6.2.5, F6.2.8, F6.5.a |
| DP-F6-03 | **Esquema `truth.jsonl` v3** para segmentos estructurales: el contenido lo fijan REQ-IND-44/47/48, REQ-REC-14–17/25 y DTC-061; la forma literal (claves, tipos de registro, cómo se expresan orientación, copia, correspondencia y estado activo/modificado/eliminado) no la fija nadie | (a) nuevos tipos de registro junto a los de v2; (b) ampliar los fragmentos vigentes `source`/`inserted` con orientación y familia de copia | F6.2.7, F6.4.3, F6.5 |
| DP-F6-04 | **Catálogos RNG e ID de F6.** Ámbitos literales para recuento, posición, origen, destino y orientación por clase estructural y por repetición (patrón DTC-071/086/089); clases de identificador para referencia derivada, evento estructural, copia duplicada, inserción nueva y caso VCF | (a) `subseed-scopes-v9` + `id-catalog-v5` con componentes al estilo F5; (b) catálogo F6 propio por decisión, como `f4-technical-v3` | F6.1.4, F6.2.1, F6.4.6, F6.5.a.2 |
| DP-F6-05 | **Representación VCF por clase de SV y catálogo de pérdidas.** Símbolico `<DEL>/<INV>/<DUP>` con `SVTYPE/END/SVLEN` frente a alelos explícitos por tamaño; cómo se expresa el destino de una duplicación insertada; BND con `MATEID`/`EVENT` para translocaciones; el catálogo cerrado y versionado de códigos `GPS_LOSS`; el esquema `vcf-complement-v1` (hoy solo una constante reservada en `src/vcf_representation.rs`); la cardinalidad positiva multicaso | (a) símbolico por defecto y explícito bajo un umbral; (b) explícito siempre que quepa; (c) BND para todo lo que no sea DEL/INV puras | F6.5.1–F6.5.3 |
| DP-F6-06 | **Campos VCF de repeticiones y criterio de «tamaño manejable»** (REQ-VER-20/22/23, ACE-31) | (a) los campos de repetición en tándem que VCF 4.4 introdujo y 4.5 conserva (`RN`, `RUS`, `RUL`, `RUC`, `RB`, …); (b) campos `INFO/GPS_*` propios declarados; el umbral de tamaño en pares de base | F6.5.1 en su parte de repeticiones |
| DP-F6-07 | **Realización de pureza e interrupciones.** La especificación define pureza como proporción del tracto que sigue el motivo, pero no dice cómo se colocan las interrupciones ni qué bases toman | (a) posiciones uniformes sobre las unidades y sustitución uniforme A/C/G/T ≠ base del motivo; (b) interrupción por unidad completa; cualquiera exige decisión | F6.4.1, F6.4.2 |
| DP-F6-08 | **Regiones «altamente homólogas»** (REQ-REF-20, REQ-CONF-96): las duplicaciones producen copias idénticas; una copia divergente exige un modelo de divergencia que nadie ha fijado | (a) F6.1 cubre solo repetitivas/duplicadas y la divergencia va a un DTF nuevo; (b) divergencia como porcentaje de sustituciones uniformes | F6.1.5 en su parte de homología |
| DP-F6-09 | **Contenido de F6.3 y forma de la ploidía por secuencia.** La hoja de ruta titula F6.3 «Ploidía generalizada»; la matriz asigna REQ-CONF-06 a F6.1 y deja en F6.3 solo REQ-EVO-49. Forma TOML de las sustituciones por cromosoma o *scaffold*. Y si un escenario F5 con ploidía heterogénea pasa a admitirse (DTC-089 hoy lo rechaza) | (a) mover REQ-CONF-06 a F6.3 en la matriz; (b) mantener F6.1 y dar a F6.3 otro contenido; F5 heterogénea exige reabrir DTC-089 | F6.3.1 |
| DP-F6-10 | **Distribución genotípica explícita** (REQ-EVO-35/36/37/49, REQ-CONF-130, REQ-CONF-32): a qué fase y a qué fuente poblacional pertenece. F5 la rechaza por contrato cerrado; la matriz F2 es HWE/`F`; ninguna fuente actual la consume | (a) F6.3 con un consumidor que hoy no existe; (b) F7 con la vía evolutiva; (c) decisión funcional que la aplace | F6.3.3 (no se planifica) |
| DP-F6-11 | **Filas varadas de la renumeración.** REQ-EVO-04/05/07/12 y REQ-IND-25 («ciclo de vida demográfico») en F6.1; REQ-CONF-39 en F6.5; REQ-CONF-106/107, REQ-VAL-40 («0,5/C controlado») y REQ-MAN-81 («frecuencias poblacionales») en F6.2. Son requisitos poblacionales o evolutivos sin relación con referencias ni SV | Reasignar a F7 o a F5 cerrada según corresponda, en una corrección de la matriz que no toque las 46 filas que `tests/f5_6_closure.rs` relee | Sección 4 (no se planifica) |
| DP-F6-12 | **DTPD-104 y DTPD-105**, abiertas y declaradas bloqueantes de F6.6: adaptadores SAM/BAM/CRAM e índices; proceso externo minimap2 (descubrimiento, versión, comando, aislamiento, diagnóstico) | Las que enumeran sus filas en `decisiones-tecnicas-0.1.md` §7.1 | F6.6.4–F6.6.6 |
| DP-F6-13 | **Proyección de la verdad a una referencia de mapeo distinta de la estructural** (REQ-DAT-03, REQ-IND-20, REQ-VER-09). La arquitectura nombra al servicio de coordenadas como única frontera, pero no dice cómo se componen las diferencias de un individuo sobre A con la correspondencia A↔B para construir el VCF sobre B sin «transformar registros» | (a) composición de cadenas de segmentos por la correspondencia estructural, con pérdida declarada donde no hay homología; (b) restringir v1 a referencia de mapeo igual a la estructural y registrar un DTF | F6.6.1 en su parte cruzada |
| DP-F6-14 | **Herramientas de la puerta externa DTC-066.** El entorno tiene `bcftools` 1.21, `tabix`/HTSlib 1.19, `samtools` 1.21 y `minimap2` 2.26; la puerta exige `bcftools` y `tabix` 1.24 y lleva pendiente desde F2 | (a) instalar 1.24; (b) reabrir DTC-066 para fijar 1.21 tras evaluarla; la puerta no se etiqueta como superada con otra versión | F6.5.5, F6.6.5 |
| DP-F6-15 | **Semántica de «subpaso diferido» F6.5.a**: si su cierre es condición de cierre de F6 o puede quedar abierto al revisar la fase | (a) obligatorio para cerrar F6, ordenado el último; (b) diferible con destino explícito | Pausa de F6.5 y matriz de cierre |
| DP-F6-16 | **DTF-055** registra BAM/CRAM y mapeo real como ampliación con destino F6.6, y el registro exige «nueva especificación funcional» para cerrar un DTF; pero los REQ de BAM/CRAM ya existen. Solo las anomalías PE (REQ-DAT-65/66) carecen de especificación | Aclarar que DTF-055 cubre únicamente las anomalías PE y que BAM/CRAM se rigen por los REQ vigentes | F6.6.5 |

**Consecuencia sobre la ejecutabilidad.** F6.1, F6.2 y F6.4 son planificables
en su núcleo y quedan bloqueadas solo en sus contratos externos (DP-01 a
DP-04, DP-07, DP-08). F6.3 depende de DP-09 y no incluye DP-10. F6.5 depende
de DP-05/06 y de la herramienta DP-14. F6.6 es planificable en conjuntos VCF
por referencia y matrices por referencia, y **no es planificable al nivel de
paso** en BAM/CRAM y mapeo hasta cerrar DP-12, DP-13 y DP-16.

## 4. Requisitos asignados a F6 por la matriz que este plan no planifica

Se enumeran para que ninguno quede sin destino. Ninguno pertenece a
«referencias y variación complejas»; su asignación actual procede de la
numeración anterior al 31 de agosto de 2026 y requiere DP-F6-10 o DP-F6-11.

| Requisito | Hito en la matriz | Por qué no entra aquí |
|---|---|---|
| REQ-EVO-04, REQ-EVO-05, REQ-EVO-07, REQ-EVO-12 | F6.1 | Propiedades de escenario evolutivo, continuidad entre generaciones y aislamiento del muestreo: vía 2, sin consumidor antes de F7 |
| REQ-IND-25 (parte «ciclo de vida demográfico») | F6.1 | La no reutilización de un ID tras dejar la población exige generaciones (F7) |
| REQ-CONF-39 | F6.5 | Unicidad de subpoblaciones evolutivas: vía 2 |
| REQ-CONF-130, REQ-EVO-35, REQ-EVO-36, REQ-EVO-37, REQ-CONF-32 | F6.2 | Distribución genotípica explícita y restricciones individuales sobre un objetivo poblacional: sin fuente que las consuma (DP-F6-10) |
| REQ-EVO-49 | F6.3 | Prevalencia local de esa misma distribución (DP-F6-10) |
| REQ-CONF-106, REQ-CONF-107, REQ-VAL-40, REQ-MAN-81 (parte F6.2) | F6.2 | «Frecuencias poblacionales» y «0,5/C controlado» son F5 cerrada o F7; F6 no crea frecuencias |
| REQ-EVO-30, REQ-CRU-18 | F6.4 | Aplicación del modelo escalonado **durante la evolución o la transmisión familiar**: F6.4 entrega el servicio y lo prueba en creación directa; el consumo generacional es F7.2/F8 |
| REQ-IND-51 | — | La estimación de recursos por duplicaciones no está asignada a F6; el presupuesto es DTPD-115 (F9). F6.2.3 solo aplica el máximo local y la advertencia |
| REQ-IND-53, REQ-IND-19 en su parte poblacional | — | Mezcla de referencias en una unidad reproductiva: F7/F8 |

## 5. F6.1 — Referencias estructurales

Resultado: una o más referencias derivadas de una referencia disponible por
duplicaciones, deleciones, inversiones y translocaciones solicitadas, cada una
con FASTA canónico propio, identidad, linaje y mapa reconstruible, y un motor
estructural puro reutilizable por F6.2.

- ☐ **F6.1.1.** Crear los contratos externos de F6 en versiones nuevas y
  dejar intactos los heredados. **Bloqueado por DP-F6-01.**
  - Depende de: cierre de F5.
  - Procede de: REQ-CONF-13/14, REQ-REF-08–11, REQ-MAN-07, DTC-036 y
    arquitectura §9.3.
  - Cambia: esquemas nuevos en `schemas/`, lectores en
    `src/configuration.rs`, `DeclaredConfig`/`EffectiveConfig`, escritor de
    manifiesto general. Ningún esquema v1–v5 ni `truth-v1/v2` cambia.
  - Se comprueba con: pareja mínima válida; rechazo de versiones anteriores
    como F6; campo desconocido rechazado (`additionalProperties: false`);
    comparación byte a byte de todos los esquemas heredados; documento
    legítimo próximo a 65536 bytes admitido y su límite más uno rechazado;
    validación externa con `python3 -m jsonschema` (Draft 2020-12, ajena al
    productor). **Debe fallar contra:** un lector que acepte `schema_version`
    anterior como F6; un esquema heredado con una propiedad añadida.
- ☐ **F6.1.2.** Implementar el motor estructural puro sobre una colección de
  secuencias en memoria: duplicación segmentaria (intra e intersecuencia,
  orientación directa o invertida), deleción interna, inversión interna y
  translocación recíproca y no recíproca como evento compuesto con subeventos.
  - Depende de: ninguno (sin E/S ni TOML).
  - Procede de: REQ-IND-40–48, REQ-IND-52, REQ-IND-58, REQ-REC-18/19,
    ACE-87/95 y arquitectura §5.2 (DTC-056/061).
  - Cambia: módulo nuevo de dominio sobre `src/coordinates.rs`
    (`InterbasePosition`, `GenomicInterval`, `InsertionAnchor`,
    `InternalDeletion`) sin crear otra convención de coordenadas.
  - Se comprueba con: una secuencia conocida de 40 bases por operación y su
    resultado escrito a mano; reglas terminales —inserción antes de la primera
    o después de la última base, deleción o inversión que toca un extremo o
    abarca la secuencia completa, punto de ruptura terminal, secuencia vacía—
    rechazadas con el evento identificado; fuente de duplicación que toca un
    extremo o abarca la secuencia completa **aceptada**; destino de duplicación
    parcial dentro de su fuente rechazado, adyacente aceptado, y la única
    excepción de fuente completa aceptada; `N` interior copiada o movida sin
    interpretarse; dos eventos no solapados y dos adyacentes aceptados;
    solapamiento o anidamiento en la misma secuencia rechazado. **Debe fallar
    contra:** una inversión que no complemente (solo invierta orden); una
    duplicación invertida que complemente sin invertir; un destino interpretado
    sobre el estado posterior en vez del anterior (REQ-IND-52: el caso usa una
    deleción previa que desplaza el destino); una regla terminal que admita
    `start == 1`; una excepción de fuente completa que deje pasar una fuente
    parcial de longitud `L-1`.
- ☐ **F6.1.3.** Producir la correspondencia origen↔derivada por segmentos con
  orientación, solo donde la transformación la hace significativa, y la
  reconstrucción inversa allí donde existe.
  - Depende de: F6.1.2.
  - Procede de: REQ-REF-10, REQ-IND-52, REQ-REC-14, DTC-056 y arquitectura
    §5.2 («nunca se inventan para regiones sin homología»).
  - Cambia: servicio de coordenadas y transformaciones (extensión, no
    duplicado).
  - Se comprueba con: ida y vuelta base a base sobre las regiones con
    homología; ausencia de correspondencia para bases eliminadas y para la
    copia insertada nueva; orientación invertida en la copia invertida.
    **Debe fallar contra:** una correspondencia que asigne a la base eliminada
    la base vecina; una que ignore la orientación en la inversión (la ida y
    vuelta sobre una secuencia palindrómica **no** discrimina: el fixture usa
    una secuencia sin palíndromos, y la docstring lo dice).
- ☐ **F6.1.4.** Derivar una referencia estructural: ID propio, FASTA canónico
  v1 (80 bases por línea), SHA-256, clase de secuencia conservada, linaje
  origen+cambios en el catálogo (grafo dirigido sin ciclos), mapa estructural
  como documento JSON cerrado con clase y versión, y comprobador independiente
  que reconstruye el FASTA desde origen+mapa antes de publicar.
  **Bloqueado por DP-F6-04** en la derivación del ID.
  - Depende de: F6.1.1, F6.1.3.
  - Procede de: REQ-REF-03/08–11/22, REQ-MAN-07, DTC-039/045 y arquitectura
    §7.2 (patrón de la reducción: plan → escritor → comprobador → publicación
    conjunta).
  - Cambia: `src/reference.rs` (`ReferenceOrigin` gana una variante de
    derivación estructural junto a `Reduced`), escritor FASTA canónico
    reutilizado, documento de mapa nuevo con el patrón de
    `src/reference_reduction_io.rs`, publicación atómica.
  - Se comprueba con: dos referencias derivadas de la misma inicial y una
    tercera derivada de una derivada (REQ-REF-11); `samtools faidx` sobre cada
    FASTA; `sha256sum` externo; reconstrucción exacta desde origen+mapa;
    identidades separadas de origen y derivada aunque un cambio vacío produzca
    bytes iguales (REQ-REF-34 por analogía; un cambio vacío se rechaza);
    ciclo de linaje rechazado; fallo de reconstrucción impide publicar FASTA,
    mapa y manifiesto como unidad. **Debe fallar contra:** un comprobador que
    lea el FASTA derivado en vez de reconstruirlo; un catálogo que acepte una
    referencia cuyo origen es ella misma.
- ☐ **F6.1.5.** Admitir regiones repetitivas o duplicadas destinadas a baja
  mapeabilidad mediante duplicaciones explícitas sobre la referencia.
  **Bloqueado por DP-F6-08** en su parte «altamente homólogas».
  - Depende de: F6.1.4.
  - Procede de: REQ-REF-20, REQ-CONF-96, ACE-52.
  - Cambia: nada más que casos y configuración; no crea otra vía.
  - Se comprueba con: una duplicación intersecuencia declarada como región de
    baja mapeabilidad produce dos copias idénticas verificadas por `grep -c`
    sobre el FASTA canónico. **Debe fallar contra:** una duplicación que
    altere una base de la copia.
- ☐ **F6.1.6.** Registrar en configuración efectiva y manifiesto general cada
  referencia disponible, su origen, sus cambios en orden y el enlace al mapa.
  - Depende de: F6.1.4.
  - Procede de: REQ-MAN-07, REQ-REF-07/10, REQ-MAN-75.
  - Se comprueba con: validación externa de esquema; lista de referencias
    igual entre efectivo y manifiesto; mapa enlazado con ruta, tamaño y
    SHA-256. **Debe fallar contra:** un manifiesto que registre los cambios
    en otro orden (el caso usa dos cambios cuya permutación produce otro
    FASTA).
- ☐ **F6.1.7.** Marimo, pestaña `Genoma`: tabla editable de referencias
  estructurales y de sus cambios ordenados; importar/exportar sin alterar IDs,
  orden ni rutas; paridad UI→TOML→`EffectiveConfig`.
  - Depende de: F6.1.6.
  - Procede de: REQ-UI-01/04/05/06/13 y hoja de ruta §7.
  - Cambia: `ui/config_editor.py` o editor propio al estilo `ui/f5_editor.py`,
    `ui/gps_app.py`, `tests_ui/`.
  - Se comprueba con: `unittest` que renderiza el TOML desde el formulario y
    vuelve a analizarlo byte a byte; reordenar filas cambia el TOML; `marimo
    check --strict`. **Debe fallar contra:** un renderizador que emita los
    cambios en el orden de un diccionario en vez del orden de la tabla.

> Revisión del hito F6.1: comprobar los cambios y sus evidencias, incluidas
> las dos referencias derivadas reconstruidas y el motor puro. No continuar
> hasta confirmar el resultado.

## 6. F6.2 — Variantes estructurales

Resultado: individuos cuyas copias haplotípicas contienen duplicaciones,
deleciones estructurales, inversiones y translocaciones explícitas, con
identificadores estables por copia e inserción, fase, dosis y trazabilidad
sin pérdida, sobre la vía directa y sobre los fundadores de F5 intactos.

- ☐ **F6.2.1.** Extender la cadena de segmentos DTC-061 con clases de segmento
  estructural: copia duplicada con ID estable propio y vínculo al intervalo
  fundador, segmento invertido, segmento trasladado, intervalo eliminado
  conservado y secuencia insertada sin fuente con ID propio; materialización
  bajo demanda con complemento inverso y `N` conservada. **Bloqueado por
  DP-F6-04** en las clases de identificador.
  - Depende de: F6.1.2, F6.1.3.
  - Procede de: REQ-IND-28/29/44, REQ-REC-14–17, DTC-061, ACE-87.
  - Cambia: `src/segments.rs` (`Segment`, `SegmentChain`),
    `src/materialization.rs` (`materialize_copy`, `materialize_copy_interval`),
    `src/identity.rs` (`IdClass`).
  - Se comprueba con: materialización de una copia con cada clase frente a la
    secuencia escrita a mano; dos copias del mismo intervalo con IDs distintos
    y el mismo vínculo fundador; la deleción conserva el intervalo ausente y
    distingue «eliminado» de «nunca perteneció». **Debe fallar contra:** una
    materialización que ignore la orientación del segmento; un ID de copia
    derivado del ordinal de aplicación en vez de sus causas estables (el caso
    aplica dos eventos en orden inverso y exige los mismos IDs).
- ☐ **F6.2.2.** Aplicar en la copia haplotípica las reglas de composición:
  eventos estructurales no solapados y adyacentes admitidos; solapamiento,
  anidamiento o modificación de una estructura previa rechazados con los
  eventos identificados; haplotipos distintos pueden solaparse; SNP e
  *indels* dentro de segmentos duplicados, invertidos o insertados con
  coordenadas actuales y procedencia; SNP o *indel* que solape o atraviese un
  punto de ruptura rechazado, adyacente admitido; eventos explícitos en orden
  y objetivo eliminado detectado.
  - Depende de: F6.2.1.
  - Procede de: REQ-REC-18/19/20/24/27, REQ-IND-54, ACE-87/103.
  - Cambia: compositor (`src/small_event_composition.rs`) y validación
    estructural del motor de verdad.
  - Se comprueba con: fixture de dos SV adyacentes admitidas y las mismas
    desplazadas una base para solapar, rechazadas; SV solapadas en copias
    distintas admitidas; SNP dentro de una copia duplicada con coordenada
    actual verificada; *indel* que atraviesa un punto de ruptura rechazado y
    el mismo *indel* una base más allá admitido; evento que apunta a una base
    eliminada por otro anterior rechazado nombrando ambos. **Debe fallar
    contra:** una comprobación de solapamiento con intervalo cerrado en vez de
    semiabierto (la adyacencia exacta pasaría a rechazarse); una validación que
    no reenumere el estado tras cada evento explícito (el caso apunta a una
    base que solo existe antes del evento previo).
- ☐ **F6.2.3.** Aplicar el máximo local de copias, predeterminado 10 (original
  más nueve), acumulativo por familia de copias e intervalo fundador y por
  copia haplotípica; valor mayor admitido con advertencia; sin otro máximo
  absoluto.
  - Depende de: F6.2.1.
  - Procede de: REQ-IND-49/50, ACE-96. La estimación de recursos de REQ-IND-51
    no es de F6 (sección 4).
  - Se comprueba con: nueve duplicaciones de la misma fuente admitidas y la
    décima rechazada; la décima admitida con máximo 11 y advertencia
    registrada; dos fuentes distintas no comparten contador; dos copias
    haplotípicas no comparten contador; una duplicación de una copia
    duplicada (misma familia) cuenta contra el mismo límite. **Debe fallar
    contra:** un contador por evento en vez de por familia (el caso encadena
    copia de copia); un límite aplicado por individuo en vez de por copia.
- ☐ **F6.2.4.** Hacer configurable el umbral que separa *indel* de variante
  estructural (predeterminado 50 pb), aplicarlo igual al generar, clasificar y
  contar, y registrarlo en efectivo y manifiesto.
  - Depende de: F6.1.1.
  - Procede de: REQ-CONF-55, REQ-IND-31, REQ-MUT-09, REQ-VER-16, REQ-MAN-28.
  - Cambia: admisión, `src/indel_application.rs`, estadísticas por clase.
  - Se comprueba con: una deleción de 49 pb clasificada *indel* y una de 50
    SV con el predeterminado; con umbral 30, la de 49 pasa a SV; efectivo y
    manifiesto conservan el valor. **Debe fallar contra:** una clasificación
    con `>` en vez de `>=`; una estadística que use el literal 50 en vez del
    umbral efectivo (el caso usa umbral 30 y cuenta clases).
- ☐ **F6.2.5.** Admitir eventos estructurales explícitos en la configuración
  para la vía directa y para los fundadores F5, con referencia estructural,
  copia, coordenadas base 1 cerradas y, por clase, fuente, destino,
  orientación y subeventos. **Bloqueado por DP-F6-01 y DP-F6-02.**
  - Depende de: F6.1.1, F6.2.2.
  - Procede de: REQ-CONF-08, REQ-IND-05/38/39, REQ-DIR-01, REQ-EVO-42.
  - Cambia: admisión (`direct_events` gana clases o su equivalente decidido),
    integración inicial (`src/initial_individual_integration.rs`), y la ruta
    F5 tras `src/f5_founders.rs` sin tocar su muestreo.
  - Se comprueba con: un individuo directo con las cuatro clases; un fundador
    F5 con una inversión; coordenadas convertidas por el servicio común y
    registradas con su convenio; ID inexistente de individuo o copia
    rechazado. **Debe fallar contra:** una conversión base 1→0 aplicada dos
    veces (el caso usa una inversión de una sola base en la posición 1 tras
    una base intacta: el doble desplazamiento la saca de rango o la mueve).
- ☐ **F6.2.6.** Conservar la referencia estructural de cada individuo, permitir
  que individuos distintos usen referencias disponibles distintas y publicar
  esa procedencia por muestra.
  - Depende de: F6.1.4, F6.2.5.
  - Procede de: REQ-IND-02/03, REQ-IND-19, ACE-09/21/22. La regla de mezcla en
    unidades reproductivas (REQ-IND-53) es de F7/F8.
  - Cambia: `Individual` (`src/haplotype.rs`) ya nace de una referencia;
    la verdad y el manifiesto publican su ID y huella.
  - Se comprueba con: dos individuos sobre referencias distintas; la región de
    generación de cada uno se interpreta sobre la suya; el manifiesto
    distingue ambas. **Debe fallar contra:** una integración que resuelva la
    región sobre la referencia inicial para todos (el caso usa una región que
    solo existe en la derivada).
- ☐ **F6.2.7.** Proyectar la verdad estructural a `truth.jsonl` v3 y validarla
  con un reconstructor independiente que materializa cada copia y vuelve a
  obtener genotipos y dosis. **Bloqueado por DP-F6-03.**
  - Depende de: F6.2.1–F6.2.6.
  - Procede de: REQ-IND-48, REQ-REC-14/15/16/17/25, REQ-VER-27, DTC-061,
    REQ-SAL-02.
  - Cambia: `src/truth_schema.rs`, `src/truth_projection.rs`,
    `src/truth_validation.rs`, `schemas/truth-v3.schema.json`; v1 y v2 no
    cambian.
  - Se comprueba con: validación externa de cada línea; reconstrucción de la
    secuencia final de cada copia igual a la materializada por el productor
    **sin usar su código**; una mutación aislada por invariante (orientación,
    ID de copia, intervalo eliminado, estado activo/eliminado, dosis). **Debe
    fallar contra:** un reconstructor que consuma los bytes materializados
    publicados por el productor en vez de recomponerlos; un esquema que admita
    un segmento sin orientación.
- ☐ **F6.2.8.** Demostrar que los fundadores F5 quedan intactos.
  **Bloqueado por DP-F6-02.**
  - Depende de: F6.2.5, F6.2.7.
  - Procede de: REQ-EVO-42, DTC-089, hoja de ruta §F6 y §revisión de F6.
  - Se comprueba con: (1) una ejecución F5 y la misma con la sección F6 vacía
    publican los siete TSV, `truth.jsonl` v2 fundadora y el manifiesto F5
    **byte a byte iguales**; (2) con eventos F6, `f5-haplotypes.tsv`,
    `f5-frequencies.tsv` y `f5-ld.tsv` no cambian y la verdad v3 referencia
    a la fundadora sin reescribirla; (3) la capa F6 nunca muestrea alelos SNP.
    **Debe fallar contra:** una capa F6 que reordene o vuelva a muestrear una
    copia fundadora (mutante que cambia un alelo SNP al aplicar una inversión
    que lo contiene: la matriz H×S publicada debe seguir igual y el alelo
    invertido debe aparecer complementado en la verdad v3, no en los TSV).
- ☐ **F6.2.9.** Marimo, pestaña `Individuos`: tabla de eventos estructurales
  explícitos con selectores de individuo, copia y referencia; umbral
  indel/SV; paridad.
  - Depende de: F6.2.5.
  - Procede de: REQ-UI-01/04/05/13.
  - Se comprueba con: `unittest` de ida y vuelta y de selectores no
    permitidos diagnosticados sin resolverlos en Python. **Debe fallar
    contra:** un control que resuelva la clase *indel*/SV por su cuenta (el
    caso cambia el umbral y exige que el TOML no cambie la clase declarada).

> Revisión del hito F6.2: comprobar reconstrucción independiente, límites,
> solapamientos e intactitud de los fundadores. No continuar hasta confirmar
> el resultado.

## 7. F6.3 — Ploidía generalizada

Resultado: ploidía global predeterminada con sustituciones por cromosoma o
*scaffold*, idéntica para todos los individuos en cada secuencia, con dosis,
fase y `GT` coherentes. La hoja de ruta titula así el hito; la matriz asigna
REQ-CONF-06 a F6.1 (DP-F6-09). Este plan lo sitúa aquí porque es el único
contenido propio del hito; la matriz debe corregirse o el paso moverse, pero
el paso existe una sola vez.

- ☐ **F6.3.1.** Resolver la ploidía efectiva por secuencia: valor global
  predeterminado y sustituciones por cromosoma o *scaffold*, cada uno entero en
  `1..=10`, misma ploidía para todos los individuos en cada secuencia, secuencia
  sin valor efectivo o fuera de dominio rechazada identificando la secuencia.
  **Bloqueado por DP-F6-09** en la forma TOML y en la admisión F5 heterogénea.
  - Depende de: F6.1.1.
  - Procede de: REQ-CONF-06/65, REQ-IND-06/07/08/37, REQ-MAN-05, ACE-11.
  - Cambia: admisión; `Individual::from_reference` (`src/haplotype.rs`) pasa
    a recibir ploidía por secuencia; `Ploidy` se conserva como tipo.
  - Se comprueba con: global 2 con una secuencia a 4 produce dos copias en
    una y cuatro en otra en **todos** los individuos; sustitución sobre una
    secuencia inexistente rechazada; 0 y 11 rechazados; efectivo y manifiesto
    conservan global y sustituciones. **Debe fallar contra:** una resolución
    que aplique la sustitución solo al primer individuo; una que lea la
    sustitución por posición de la lista en vez de por `sequence_id`.
- ☐ **F6.3.2.** Expresar estado alélico, dosis y `GT` con la ploidía de cada
  secuencia en verdad y VCF.
  - Depende de: F6.3.1, F6.2.7.
  - Procede de: REQ-IND-12, REQ-VER-27/36/44.
  - Cambia: `src/truth_projection.rs`, `src/vcf_projection.rs`.
  - Se comprueba con: un individuo con ploidía 2 en `chr1` y 4 en `chr2`
    publica `GT` con dos y cuatro alelos faseados respectivamente; el oráculo
    semántico independiente (`src/vcf_semantic_validation.rs`) recalcula la
    dosis desde los haplotipos; `bcftools view` lee el fichero. **Debe fallar
    contra:** una proyección que tome la ploidía global para todas las
    secuencias.
- ☐ **F6.3.3.** Prevalencia local de una distribución genotípica explícita
  (REQ-EVO-49). **No se planifica**: depende de DP-F6-10.
- ☐ **F6.3.4.** Marimo: control de ploidía global y tabla de sustituciones por
  secuencia en `Genoma`; `F` de F5 solo visible con ploidía 2 en todas las
  secuencias elegibles, como exige REQ-UI-15.
  - Depende de: F6.3.1.
  - Procede de: REQ-UI-01/03/15.
  - Se comprueba con: `unittest` de ida y vuelta y de ocultación de `F`.
    **Debe fallar contra:** un formulario que muestre `F` con una sustitución
    a 4.

> Revisión del hito F6.3: comprobar ploidía por secuencia en verdad, VCF e
> interfaz. No continuar hasta confirmar el resultado.

## 8. F6.4 — Repeticiones en tándem

Resultado: loci de repetición en tándem con motivo de longitud libre, copias o
longitud total, pureza e interrupciones, definidos o generados, en referencias
simuladas y en haplotipos; el modelo escalonado simple como servicio
reutilizable; mutaciones internas desactivadas por defecto; rechazo de todo lo
relativo a elementos transponibles.

> **Nota del 27 de septiembre de 2026 — cinco pasos de este hito se han
> adelantado y están hechos.** El plan
> `design/plan-implementacion-repeticiones-pcr-0.1.md`, que cruza tres fases,
> ha implementado **F6.4.1, F6.4.2, F6.4.3, F6.4.4 y F6.4.8**, y ha abierto
> coordenadas nuevas **F6.4.11 a F6.4.13** para el generador de alelos por
> deslizamiento biológico, que aquel plan llama su entrega 5. Este documento
> **no los replanifica**: lo que aquí queda de F6.4 es lo que sigue sin hacer.
>
> Dos avisos para quien lea las casillas de abajo sin haber leído aquel plan:
>
> 1. **F6.4.4 está hecho con menos alcance del que esta casilla describe.** El
>    servicio escalonado existe, con su sesgo, su regla del borde y su recuento
>    de eventos no aplicados, pero **la tasa y el sesgo se declaran por locus y
>    no por región**, y su consumo desde F7.2 y F8 sigue sin ocurrir porque esas
>    fases no existen. Lo que sí está fijado es el contrato, en **DTC-094**.
> 2. **El generador de alelos es una sola tirada de inicialización.** Aplicar el
>    modelo escalonado **generación tras generación** sigue perteneciendo a la
>    fase de evolución y **no** se ha invadido. Lo fija **DTC-095**.
>
> Las decisiones normativas del trabajo adelantado son **DTC-094** —el servicio
> escalonado—, **DTC-095** —el generador—, **DTC-096** —su ámbito de azar
> propio, separado del tartamudeo— y **DTC-097** —la publicación del origen de
> cada alelo—.

- ☐ **F6.4.1.** Modelo de dominio del locus de repetición: motivo A/C/G/T no
  vacío, longitud de unidad coincidente con el motivo, número de copias o
  longitud total, pureza y patrón de interrupciones, y su secuencia
  realizada; validación previa con locus y error identificados; motivo con
  `N` rechazado. **Bloqueado por DP-F6-07** en la realización de
  interrupciones.
  - Depende de: ninguno.
  - Procede de: REQ-CONF-57/60, REQ-IND-32/33, REQ-REF-23, ACE-31/110.
  - Se comprueba con: `(ACG)×5` produce 15 bases; longitud declarada 4 para
    `ACG` rechazada; motivo vacío y motivo `ACN` rechazados; pureza 1 no
    produce interrupciones; secuencia realizada solo A/C/G/T. **Debe fallar
    contra:** un validador que compare la longitud declarada con la del
    motivo en mayúsculas pero acepte `acn`; una realización de pureza 1 que
    introduzca una interrupción.
- ☐ **F6.4.2.** Incluir repeticiones en una referencia simulada, en loci
  concretos o generados por configuración aleatoria y regional, reflejadas en
  el FASTA canónico.
  - Depende de: F6.4.1, F6.1.4.
  - Procede de: REQ-REF-19, REQ-CONF-56.
  - Cambia: `src/reference_generation.rs`.
  - Se comprueba con: el FASTA contiene la secuencia realizada del locus en
    la coordenada declarada (`grep -ob` externo); misma semilla, mismo FASTA;
    otra semilla mueve los loci aleatorios y no los concretos. **Debe fallar
    contra:** una escritura que coloque el locus una base más allá.
- ☐ **F6.4.3.** Representar alelos de repetición por haplotipo como
  modificación tipada de la cadena de segmentos, con fase, copias o longitud
  y secuencia realizada, y permitir que difieran entre haplotipos, individuos
  y poblaciones. **Bloqueado por DP-F6-03** en la forma de la verdad.
  - Depende de: F6.4.1, F6.2.1.
  - Procede de: REQ-IND-32/33, REQ-CONF-58/59, ACE-103.
  - Se comprueba con: dos copias de un diploide con 5 y 7 unidades
    materializan longitudes distintas; la verdad v3 asigna cada alelo a su
    copia; reconstrucción independiente. **Debe fallar contra:** una
    materialización que use el número de copias de la referencia para ambas.
- ☐ **F6.4.4.** Implementar el modelo escalonado simple como servicio: cada
  evento añade o elimina exactamente una unidad completa según el sesgo
  efectivo, con tasa y sesgo configurables globalmente y por locus o región;
  rechazo previo de otro modelo o de un salto de varias unidades con parámetro
  y valor; aplicable en creación directa; F7.2 y F8 lo consumirán sin
  reimplementarlo.
  - Depende de: F6.4.3.
  - Procede de: REQ-CONF-61/62, REQ-IND-35, REQ-DIR-04, REQ-MUT-11,
    REQ-VAL-48, ACE-32/109; DTF-015 mantiene fuera el salto múltiple.
  - Se comprueba con: 1000 eventos con sesgo 1 son todos expansiones y con
    sesgo 0 todas contracciones (0 y 1 no consumen RNG, DTC-088
    `bernoulli-v1`); cada evento cambia la longitud exactamente en la unidad;
    `model = "geometric"` o `step_units = 2` rechazados antes de simular;
    misma semilla reproduce la secuencia de decisiones. **Debe fallar
    contra:** un paso que cambie dos unidades; un sesgo interpretado al revés
    (el caso con sesgo 1 exige solo expansiones).
- ☐ **F6.4.5.** Mutaciones internas de interrupciones, pureza o motivo
  desactivadas por defecto y aplicables solo con activación expresa, sobre un
  haplotipo concreto y trazadas. **Bloqueado por DP-F6-07.**
  - Depende de: F6.4.4.
  - Procede de: REQ-CONF-63, REQ-IND-36, REQ-DIR-05, REQ-MAN-30.
  - Se comprueba con: configuración sin la clave no produce ninguna mutación
    interna aunque la tasa global sea alta; activada, cada evento se registra
    aparte del paso unitario. **Debe fallar contra:** un predeterminado
    activado.
- ☐ **F6.4.6.** Generar loci de repetición al azar y por regiones con las
  mismas reglas globales, regionales, individuales y reproducibles que las
  demás clases: recuento por `poisson-count-v1`, selección dinámica sobre el
  estado vigente con `dynamic-indel-selection-v1` o su equivalente para esta
  clase, exclusión de candidatos cuyo motivo contendría `N`, fallo si el
  espacio se agota. **Bloqueado por DP-F6-04** en los ámbitos RNG.
  - Depende de: F6.4.3.
  - Procede de: REQ-CONF-56, REQ-IND-34, REQ-IND-56/57, DTC-067, ACE-107.
  - Cambia: `src/indel_application.rs` o módulo hermano; `src/randomness.rs`.
  - Se comprueba con: vectores congelados del recuento; media 0 no consume
    RNG; candidato sobre `N` excluido (fixture con una única ventana válida y
    otra con `N`: 1000 sorteos nunca eligen la segunda); espacio insuficiente
    falla y registra sorteado/aceptado/activo. **Debe fallar contra:** un
    selector que reintente candidatos inválidos (el caso cuenta consumos RNG:
    reintentar consume más); un fallo por agotamiento convertido en
    advertencia.
- ☐ **F6.4.7.** Rechazar toda clave, catálogo, salida o metadato de elementos
  transponibles como capacidad retirada, con diagnóstico identificable.
  - Depende de: F6.1.1.
  - Procede de: REQ-TE-01–14 (retirados), REQ-MAN-67, REQ-SAL-14, DTF-007.
  - Se comprueba con: un TOML con una tabla `transposable_elements` rechazado
    por campo desconocido; ausencia de cualquier campo TE en efectivo,
    manifiesto y verdad (`grep -ci transposable` sobre los tres devuelve 0).
    **Debe fallar contra:** un esquema con `additionalProperties: true` en
    la tabla afectada.
- ☐ **F6.4.8.** Registrar en efectivo y manifiesto loci y parámetros
  configurados, procedencia concreta o aleatoria, motivo, unidad, copias o
  longitud, objetivos y valores realizados de pureza e interrupciones, y por
  mutación el modelo, tasa y sesgo efectivos, locus, haplotipo y dirección.
  - Depende de: F6.4.4–F6.4.6.
  - Procede de: REQ-MAN-29/30, REQ-VAL-48.
  - Se comprueba con: validación externa; recuento de mutaciones registradas
    igual al número de eventos aplicados; cada registro «añade» o «elimina»
    exactamente una unidad. **Debe fallar contra:** un manifiesto que registre
    la tasa declarada como efectiva cuando una región la sustituye.
- ☐ **F6.4.9.** Marimo, `Individuos` y `Genoma`: tablas de loci concretos,
  parámetros aleatorios y regionales, tasa y sesgo, activación explícita de
  mutaciones internas; sin control de salto múltiple ni de otro modelo.
  - Depende de: F6.4.8.
  - Procede de: REQ-UI-01/03/04/05.
  - Se comprueba con: `unittest` de ida y vuelta y de ausencia de controles
    para lo no admitido, al estilo de
    `tests_ui/test_f5_6_1_scenario_ui.py::test_no_panel_supplied_haplotypes_recombination_or_map_control_exists`.
    **Debe fallar contra:** un control `step_units` presente.

> Revisión del hito F6.4: comprobar loci, alelos por haplotipo, paso unitario,
> generación aleatoria y rechazo de TE. No continuar hasta confirmar el
> resultado.

## 9. F6.5 — Variantes complejas y complemento VCF

Resultado: el VCF de verdad representa SV, translocaciones y repeticiones con
la representación estándar más completa, sigue siendo VCF 4.5 válido, y el
primer caso real con pérdida de detalle produce `GPS_LOSS`, `GPS_COMPLEMENT` y
exactamente una entrada en `vcf-complement.jsonl`, reconstruible junto al VCF.
Después, el subpaso diferido F6.5.a incorpora la variación aleatoria directa.

- ☐ **F6.5.1.** Representar cada clase estructural en VCF: alelos simbólicos
  con `SVTYPE`, `END` y `SVLEN` cuando corresponda, BND relacionados por
  `MATEID`/`EVENT` para cada nueva adyacencia (dos por adyacencia, cuatro para
  una recíproca, sin cardinalidad universal), alelos de repetición explícitos
  o simbólicos según tamaño, y una versión nueva de la normalización DTC-064
  que admite varios BND en una misma posición sin relajar la cardinalidad de
  SNP e *indels*. **Bloqueado por DP-F6-05 y DP-F6-06.**
  - Depende de: F6.2.7, F6.4.3.
  - Procede de: REQ-VER-15–23/42/44/47, DTC-062/064, ACE-30/31.
  - Cambia: `src/vcf_grouping.rs`, `src/vcf_projection.rs`,
    `src/vcf_contract.rs` (cabecera declara cada `INFO` adicional).
  - Se comprueba con: una deleción de 50 pb como SV y una de 49 como *indel*
    explícito; una inversión con `N` interior conserva `N` solo si el estándar
    lo permite sin pérdida y en otro caso usa complemento; una translocación
    recíproca produce cuatro BND con `MATEID` cruzados; `bcftools view` y
    `bcftools norm --check-ref e` (oráculo externo, DP-F6-14) aceptan el
    fichero; nunca un SNP con `REF=N`. **Debe fallar contra:** un `END` igual
    a `POS+SVLEN` en una inserción (semántica oficial distinta); un `MATEID`
    que apunte al propio registro; un agrupador que fusione dos BND de la misma
    posición en un registro multialélico.
- ☐ **F6.5.2.** Publicar el catálogo cerrado y versionado de códigos
  `GPS_LOSS`, el esquema `vcf-complement-v1` y la cardinalidad positiva
  multicaso: uno o varios registros por caso, exactamente una entrada por caso
  con pérdida, `GPS_CASE` intacta. **Bloqueado por DP-F6-05.**
  - Depende de: F6.5.1.
  - Procede de: REQ-VER-24/25/26, DTC-065, REQ-MAN-31.
  - Cambia: catálogo nuevo en `catalogs/`, esquema nuevo, escritor JSONL por
    flujo con publicación atómica junto al VCF e índice.
  - Se comprueba con: texto libre en `GPS_LOSS` rechazado por el oráculo;
    dos registros del mismo caso comparten entrada; un caso representable no
    crea entrada y ningún caso perdido no crea fichero (`vcf-complement.jsonl`
    ausente y no enlazado). **Debe fallar contra:** un escritor que cree el
    fichero vacío; una entrada por registro en vez de por caso.
- ☐ **F6.5.3.** Producir el primer caso positivo real de pérdida y reconstruir
  desde VCF+complemento el cambio, el haplotipo, la fase, las coordenadas y el
  genotipo, comparándolos con `truth.jsonl` v3 mediante un lector que no
  comparte la proyección productora.
  - Depende de: F6.5.2.
  - Procede de: REQ-VER-19/26/34, REQ-MAN-27/31, DTC-066, ACE-33.
  - Se comprueba con: el caso de pérdida que DP-F6-05 determine (candidatos:
    destino de una duplicación insertada; evento compuesto de translocación)
    reconstruido exactamente; `GPS_COMPLEMENT` ausente, duplicado, cruzado o
    sobrante rechazado. **Debe fallar contra:** una entrada complementaria
    que omita la fase (el caso es heterocigoto en la segunda copia y la
    reconstrucción exige la copia exacta).
- ☐ **F6.5.4.** Representar eventos pequeños solapados de un mismo haplotipo
  como un único alelo final complejo normalizado, haplotipos distintos como
  ALT distintos con fase y dosis, y registrar el historial ordenado con marco
  de coordenadas anterior y estado activo/modificado/eliminado, también dentro
  de segmentos duplicados o invertidos.
  - Depende de: F6.2.2, F6.5.1.
  - Procede de: REQ-REC-25/28/29, REQ-VER-45, REQ-MAN-79, ACE-103.
  - Se comprueba con: SNP seguido de deleción que lo cubre → el SNP no se
    emite y consta como eliminado; dos SNP adyacentes en copias distintas →
    dos ALT en un registro con `GT` faseado; SNP dentro de una copia invertida
    emitido en su coordenada actual complementado. **Debe fallar contra:** un
    historial que marque «activo» el SNP cubierto; una proyección que emita el
    SNP invertido sin complementar.
- ☐ **F6.5.5.** Ampliar el oráculo semántico independiente y ejecutar la
  puerta externa DTC-066 sobre texto, BGZF+TBI y BGZF+CSI con SV, BND,
  repeticiones y complemento. **Bloqueado por DP-F6-14.**
  - Depende de: F6.5.1–F6.5.4.
  - Procede de: REQ-VAL-47, REQ-SAL-12, REQ-VER-01/06/23/28, REQ-FMT-05/06/07,
    DTC-063/066.
  - Cambia: `src/vcf_semantic_validation.rs`, `scripts/validate_vcf_external.sh`.
  - Se comprueba con: corpus con un caso sin complemento, uno con registro y
    entrada, una adyacencia con dos BND, una recíproca con cuatro, una
    repetición explícita y una simbólica; mutación aislada de cada invariante
    del oráculo; `bcftools` y `tabix` en la versión que DTC-066 fije, con
    versión observada registrada en la evidencia. **Debe fallar contra:** un
    oráculo que valide la cardinalidad caso↔registro leyendo `GPS_CASE` sin
    contar los registros (el caso publica un caso con dos registros y solo
    uno lleva el ID).
- ☐ **F6.5.6.** Validación estructural de cierre: convenio de coordenadas,
  límites terminales, puntos de ruptura y orientaciones de cada SV, límite
  acumulado de copias, excepción de fuente completa, correspondencias de cada
  estado anterior y posterior, y validación secuencial de eventos explícitos.
  - Depende de: F6.2.2, F6.2.3, F6.5.4.
  - Procede de: REQ-VAL-27, ACE-95/96/107.
  - Se comprueba con: inyección de una verdad v3 con una orientación
    inconsistente, con una copia de más sobre el límite y con un evento cuyo
    objetivo no existía en el estado anterior; cada una produce error y ningún
    conjunto se marca completo. **Debe fallar contra:** un validador que
    consulte el límite de copias en la configuración en vez de contar en la
    verdad.

### F6.5.a — Variación aleatoria directa (subpaso diferido)

La hoja de ruta lo declara diferido; si su cierre condiciona el cierre de F6
es DP-F6-15. Se ordena el último de F6.5 y lleva su propia pausa.

- ☐ **F6.5.a.1.** Configurar, por clase distinta del SNV de matriz, cantidad
  media esperada o densidad por base, kilobase o megabase, distribución
  espacial uniforme o regiones con valor mayor, menor o cero, con valores
  globales y sustituciones por individuo declarado y existente, resueltas por
  la matriz de ejes REQ-CONF-116–121 y registradas con regla ganadora y
  procedencia. **Bloqueado por DP-F6-01.**
  - Depende de: F6.1.1, F6.2.5, F6.4.6.
  - Procede de: REQ-CONF-08/17–23/28–31, REQ-IND-04/09/10/13/14/15/21,
    REQ-DIR-02/03, REQ-MAN-10/12, ACE-111.
  - Se comprueba con: una sustitución individual reemplaza el valor global
    solo para su destinatario; regiones solapadas con iguales selectores
    rechazadas aunque coincidan en valor; individuo inexistente rechazado;
    región frente a individuo sin regla combinada rechazada como máximas
    incomparables. **Debe fallar contra:** una precedencia implícita
    «individuo gana a región» sin regla combinada; una densidad por kilobase
    multiplicada por la longitud en bases sin dividir.
- ☐ **F6.5.a.2.** Orquestar la generación aleatoria: eventos explícitos
  primero y en orden; después, por clase y en el orden fijo duplicaciones,
  deleciones, inversiones, translocaciones, repeticiones, *indels* y SNP, un
  único sorteo Poisson del recuento al inicio de la clase sobre el espacio
  válido del estado vigente, selección solo entre candidatos válidos con
  actualización tras cada evento, sin reintentos ni límite de intentos, fallo
  sin reducción silenciosa y registro de sorteado, aceptado y activo.
  **Bloqueado por DP-F6-02 y DP-F6-04.**
  - Depende de: F6.5.a.1, F6.1.2, F6.4.6.
  - Procede de: REQ-IND-16/17/18/54–57, REQ-CONF-23/24/25, REQ-MUT-13/14,
    DTC-067, ACE-107/109.
  - Se comprueba con: fixture donde una deleción aleatoria reduce el espacio
    de una inversión posterior y el recuento de inversiones no se vuelve a
    sortear (consumo RNG contado); agotamiento en la clase tercera falla y
    conserva aceptados de las dos primeras en el diagnóstico sin publicarlos
    como verdad; orden alterado por configuración rechazado; `poisson` es la
    única distribución y otra se rechaza con parámetro y valor. **Debe fallar
    contra:** un orquestador que procese SNP antes que SV; uno que resortee el
    recuento tras cada evento; uno que publique la verdad parcial tras el
    fallo.
- ☐ **F6.5.a.3.** Hacer excluyentes en una misma población inicial la matriz de
  frecuencias SNV y la generación inicial de SNP por cantidad o densidad, y
  comprobar que la creación directa consume la sustitución individual sin
  alterar la matriz poblacional.
  - Depende de: F6.5.a.2.
  - Procede de: REQ-IND-64, REQ-DIR-02, REQ-CONF-20/22, hoja de ruta §F6.5.a.
  - Se comprueba con: población con matriz y densidad SNP a la vez rechazada;
    matriz preparada byte a byte igual antes y después de generar SV
    aleatorias sobre sus individuos. **Debe fallar contra:** una generación que
    añada filas a la matriz canónica.
- ☐ **F6.5.a.4.** Marimo, `Individuos`: cantidad/densidad con unidad,
  distribución y regiones globales y por individuo, con selectores explícitos.
  - Depende de: F6.5.a.1.
  - Procede de: REQ-UI-01/03/04/05.
  - Se comprueba con: `unittest` de ida y vuelta y de que la unidad no se
    convierte en la interfaz. **Debe fallar contra:** un control que muestre
    kb y guarde bases redondeadas.

> Revisión del hito F6.5: comprobar VCF válido con SV/BND/repeticiones, el
> primer complemento real reconstruido, la puerta externa con versión
> registrada y, si DP-F6-15 lo exige, la variación aleatoria directa. No
> continuar hasta confirmar el resultado.

## 10. F6.6 — Conjuntos y matrices por referencia

Resultado planificable ahora: selección de una o varias referencias de mapeo,
un conjunto independiente por referencia con su propio VCF multimuesta y su
propio espacio de identificadores, matrices SNV canónicas separadas por
referencia, y manifiestos que separan resultados por referencia sin
correspondencias registro a registro. Los pasos de mapeo, BAM/CRAM, índices y
validación de alineamientos **no son planificables al nivel de paso** hasta
cerrar DTPD-104/105 (DP-F6-12), DP-F6-13 y DP-F6-16; se enumeran con el
contenido que los REQ ya fijan para que conserven destino.

- ☐ **F6.6.1.** Seleccionar en el evento una o varias referencias de mapeo
  entre las disponibles y producir un conjunto por referencia con el mismo
  conjunto de muestras, cada uno con exactamente un VCF multimuesta construido
  directamente contra su referencia, identificadores de variante locales y
  sin tabla de correspondencias. **Bloqueado por DP-F6-13** cuando la
  referencia de mapeo difiere de la estructural del individuo.
  - Depende de: F6.1.4, F6.2.6, F6.5.1.
  - Procede de: REQ-CONF-15, REQ-VER-03/04/05/07/08/09/10/28/29, REQ-DAT-06/07
    (parte VCF), REQ-IND-20, ACE-15/18/20/60/99.
  - Cambia: orquestación de eventos, `src/vcf_publication.rs`,
    `src/vcf_projection.rs`, manifiesto de evento.
  - Se comprueba con: dos referencias seleccionadas producen dos VCF con las
    mismas columnas de muestra, `##contig` distintos y registros
    independientes; ningún ID de variante ni campo cruza entre ellos (`grep`
    de los IDs de uno sobre el otro devuelve 0); el oráculo semántico valida
    cada VCF contra su FASTA; con referencia de mapeo igual a la estructural
    el VCF coincide byte a byte con el de F6.5. **Debe fallar contra:** una
    publicación que reutilice el espacio de IDs entre conjuntos; una
    proyección que tome la referencia estructural del individuo en vez de la
    de mapeo para `##contig` (el caso usa una referencia de mapeo derivada
    con una secuencia renombrada).
- ☐ **F6.6.2.** Admitir la coexistencia, selección y separación de matrices
  SNV canónicas para varias referencias: cada matriz ligada por SHA-256 a su
  referencia exacta, cada población inicial usa la matriz de su referencia y
  una matriz de otra referencia se rechaza.
  - Depende de: F6.1.4, F6.2.6.
  - Procede de: REQ-CONF-128 (parte F6.6.a), REQ-IND-60/63, DTC-054, ACE-116.
  - Cambia: `src/population_snvs.rs`, admisión de `[population_snvs]`.
  - Se comprueba con: dos matrices para dos referencias derivadas de la
    misma inicial; cruzar una matriz con la otra referencia se rechaza por
    huella; los individuos de cada referencia muestrean solo su matriz.
    **Debe fallar contra:** una comprobación de huella que compare solo
    `reference_id` y no la huella canónica de secuencias (el caso usa dos
    referencias con el mismo ID de secuencia y bases distintas).
- ☐ **F6.6.3.** Registrar en manifiesto de evento y raíz las referencias de
  mapeo, el conjunto de cada una, la relación exclusiva artefacto↔conjunto,
  los estados de conjunto y la separación de directorios por referencia; un
  conjunto solo se marca completo con todos sus artefactos publicados y
  validados y un fallo lo marca fallido sin publicar parciales.
  - Depende de: F6.6.1.
  - Procede de: REQ-MAN-09/12/57/58/74, REQ-SAL-17/23, REQ-EVO-78, DTC-045.
  - Se comprueba con: validación externa del esquema; un VCF que no supera el
    oráculo deja su conjunto `fallido` y el otro `completo` en la misma
    ejecución; el directorio del evento separa por referencia. **Debe fallar
    contra:** un estado de conjunto derivado de la existencia del fichero en
    vez de la validación.
- ☐ **F6.6.4.** Adaptador de proceso externo del mapeador: exactamente un
  mapeador por evento, minimap2 predeterminado, *preset* automático para
  Illumina *paired-end* con sustituciones validadas y advertidas, versión,
  comando efectivo completo, stdout, stderr y código de salida registrados, los
  mismos FASTQ finales y el mismo mapeador para todas las referencias del
  evento. **Bloqueado por DTPD-105 (DP-F6-12).**
  - Procede de: REQ-CONF-69/70/71, REQ-DAT-21/22/23, REQ-TEC-04, REQ-MAN-36,
    REQ-VAL-25, DTC-043; arquitectura §6.4 (proceso ausente, versión
    incompatible, código no cero, salida inválida, interrupción).
- ☐ **F6.6.5.** BAM o CRAM final por individuo y conjunto: CRAM por omisión,
  BAM por elección, ninguno sin solicitud, solo desde FASTQ final del modo
  realista con FASTQ activado, ordenado por coordenadas con índice compatible,
  cabecera con `@RG` que identifica evento, muestra, biblioteca, tanda y
  plataforma, `SM` con el ID biológico estable, CRAM ligado al FASTA exacto de
  los recursos comunes, solo la versión final con degradaciones, repeticiones
  reflejadas. **Bloqueado por DTPD-104 (DP-F6-12) y DP-F6-16.**
  - Procede de: REQ-CONF-64/67/68, REQ-FMT-01/02/03, REQ-DAT-01/02/03/04/09/
    11/20, REQ-RG-01/04/05/07, REQ-REF-21, REQ-WGS-12, REQ-MAN-53/55/56,
    REQ-SAL-11, ACE-16/56/58/59/100.
- ☐ **F6.6.6.** Validación interna de alineamientos: coherencia con referencia,
  posiciones verdaderas cuando proceda, campos, banderas y grupos de lectura;
  compatibilidad alineamiento↔índice↔referencia antes de marcar completo;
  métricas observadas en resultados de validación y nunca en el VCF; enlace de
  resultados desde el manifiesto. **Bloqueado por DP-F6-12.**
  - Procede de: REQ-VAL-01/08/25, REQ-FMT-04, REQ-VER-38, REQ-MAN-59,
    REQ-SAL-13, ACE-64.
- ☐ **F6.6.7.** Marimo, `Salidas` y `Eventos técnicos`: selección de
  referencias de mapeo (planificable ahora); formato BAM/CRAM y mapeador
  (cuando F6.6.4–F6.6.6 se desbloqueen).
  - Depende de: F6.6.1.
  - Procede de: REQ-UI-01/02/05/11.
  - Se comprueba con: `unittest` de ida y vuelta; el resumen final enlaza
    manifiestos sin cargar VCF ni FASTA. **Debe fallar contra:** un control de
    referencia de mapeo que solo ofrezca la inicial.

> Revisión del hito F6.6 y de la fase F6: comprobar aislamiento entre
> referencias, matrices separadas, estados de conjunto y, si DTPD-104/105
> están cerradas, mapeo y BAM/CRAM. La revisión de F6 confirma además la
> sección 11. No continuar hasta confirmar el resultado.

## 11. Fundadores F5 intactos y realismo F4: pruebas transversales de cierre

La revisión de F6 debe confirmar que la complejidad estructural **no rompe la
verdad fundadora de F5 ni el realismo técnico de F4**. Estas pruebas se
ejecutan al cerrar la fase, además de las de cada hito.

| Invariante | Cómo se comprueba | Debe fallar contra |
|---|---|---|
| F5 sin F6 es F5 | Ejecución F5 de referencia y la misma con los bloques F6 vacíos: siete TSV, `truth.jsonl` v2 y manifiesto F5 byte a byte iguales (`sha256sum`) | Un lector F6 que reescriba un valor predeterminado en el efectivo F5 |
| Fundadores intactos con F6 activo | `f5-haplotypes.tsv`, `f5-frequencies.tsv` y `f5-ld.tsv` no cambian; la verdad v3 referencia la fundadora; el oráculo F5 (`src/f5_oracle.rs`) sigue rederivando H×S, q, D y `r²` | Una capa F6 que aplique una inversión reescribiendo alelos en la tabla H×S |
| F4 sobre individuos F6 | El recorrido F4 completo sobre un individuo con duplicación, inversión y translocación produce FASTQ; el validador por flujo de `src/fastq_validation.rs` regenera cada pareja desde la verdad v3 y coincide literalmente en bases y calidades | Una materialización WGS que lea el intervalo fuente en orientación directa para una copia invertida |
| Denominador de cobertura | La cobertura haploide usa la longitud de la referencia estructural del individuo tras F6 (REQ-WGS-04): una duplicación de 1 kb en un genoma de 100 kb cambia el denominador en 1 % exactamente | Un denominador que use la referencia inicial |
| `N` nunca se resuelve | SNP sobre `N` rechazado; SV con `N` interior conserva la `N`; lectura con `N` a Phred 0; VCF sin `REF=N`; complemento cuando el estándar no la conserva | Una copia estructural que sustituya `N` por `A` |
| Sin meiosis | Ningún módulo F6 forma bivalentes, gametos ni mosaicos; ningún ámbito RNG F6 nombra recombinación (`grep -ri "recombin\|meiosis\|bivalent" src/<módulos F6>` devuelve 0 salvo comentarios que remitan a F7.1) | — (búsqueda mecánica; se declara como comprobación, no como prueba discriminante) |

## 12. Riesgos de error silencioso y sus pruebas independientes

| Riesgo | Dónde aparece | Protección independiente |
|---|---|---|
| Desplazamiento de una base al convertir coordenadas o al aplicar un evento sobre el estado posterior | F6.1.2, F6.2.2, F6.2.5 | Secuencias escritas a mano de 40 bases; casos límite en posición 1 y `L`; fixture de deleción previa que desplaza el destino |
| Inversión sin complemento o complemento sin inversión | F6.1.2, F6.2.1, F6.5.4 | Secuencia no palindrómica; comparación con complemento inverso calculado por una función de prueba independiente |
| Copia duplicada sin identidad propia | F6.2.1, F6.2.3 | Dos copias del mismo intervalo con IDs distintos; contador de familia con copia de copia |
| Alelo estructural en la copia equivocada | F6.2.7, F6.5.3 | Reconstructor de verdad v3 que materializa cada copia y compara con el productor sin compartir código |
| Dosis con ploidía por secuencia | F6.3.2 | Oráculo semántico VCF que recalcula la dosis desde haplotipos; `bcftools view` |
| Mezcla de referencias | F6.6.1–F6.6.3 | `grep` cruzado de IDs; `##contig` contra el FASTA del conjunto; estados de conjunto independientes en una misma ejecución |
| Verdad parcial publicada tras un fallo aleatorio | F6.5.a.2 | Agotamiento inducido en la tercera clase; no existe `truth.jsonl` y el manifiesto es `fallido` |
| Puerta etiquetada sin ejecutarse | F6.5.5, F6.6.5 | La evidencia registra la versión observada de cada herramienta; una versión distinta de la normativa deja la puerta **pendiente**, no superada |

## 13. Matriz de cierre

| Hito | Resultado construido | Evidencia | Riesgo protegido | Decisiones que lo condicionan |
|---|---|---|---|---|
| F6.1 | Motor estructural puro y referencias derivadas con FASTA, mapa y linaje | Reconstrucción origen+mapa, `samtools faidx`, `sha256sum`, casos terminales y de anidamiento | Referencia derivada no reconstruible o coordenadas desplazadas | DP-01, DP-04, DP-08 |
| F6.2 | SV explícitas en copias haplotípicas con IDs, límites, umbral y verdad v3 | Reconstructor independiente, fundadores F5 byte a byte, límite de copias, adyacencia/solapamiento | Pérdida de linaje o de fase; fundadores alterados | DP-01, DP-02, DP-03, DP-04 |
| F6.3 | Ploidía por secuencia en verdad, VCF e interfaz | Dosis recalculada por el oráculo; `GT` con la ploidía de cada secuencia | Dosis o `GT` con ploidía equivocada | DP-09 (DP-10 excluido) |
| F6.4 | Repeticiones definidas y aleatorias, alelos por haplotipo, paso unitario, TE rechazados | Secuencia realizada en FASTA y copias, vectores Poisson, 0/1 sin consumo, rechazo de salto múltiple | Motivo sobre `N`, salto de varias unidades, mutación interna implícita | DP-03, DP-04, DP-07 |
| F6.5 | VCF 4.5 con SV/BND/repeticiones, catálogo `GPS_LOSS`, complemento y primer caso real reconstruido; F6.5.a variación aleatoria directa | Oráculo semántico ampliado, puerta DTC-066 con versión registrada, corpus con mutación por invariante, orden fijo y agotamiento | Pérdida silenciosa de detalle; verdad parcial publicada | DP-05, DP-06, DP-14, DP-15 |
| F6.6 | Conjuntos VCF y matrices por referencia; BAM/CRAM y mapeo cuando se desbloqueen | Dos referencias con VCF independientes, matrices separadas por huella, estados de conjunto | Mezcla de referencias o de identificadores | DP-12, DP-13, DP-16 |
| F6 | Verdad F5 ampliada sin alterar fundadores ni realismo F4 | Sección 11 completa | Complejidad estructural que rompa F5 o F4 | Todas las anteriores cerradas o con destino |

## 14. Fuera de alcance con destino

- Meiosis, mosaicos, recombinación, mapas y bivalentes: F7.1 (DTPD-111).
- Mutación heredable y *de novo* de todas las clases, incluido el consumo
  generacional del modelo escalonado: F7.2 y F7.2.a.
- Mezcla de referencias estructurales en una unidad reproductiva y
  reproducción entre referencias o con diferencias estructurales: F7/F8 y
  DTF-011/021.
- SV anidadas o que modifican una estructura previa; variantes terminales, de
  secuencia completa o copia dentro de su fuente: DTF-033/034.
- Saltos de varias unidades en repeticiones: DTF-015. Elementos transponibles:
  DTF-007.
- Ploidía variable por individuo, cromosomas sexuales y aneuploidía: DTF-003.
- Correspondencia de variantes entre referencias: DTF-008 (F6.6 la prohíbe
  expresamente).
- Varios índices por referencia: DTF-048. Alineamientos sintéticos directos:
  DTF-051. Anomalías *paired-end* observables: DTF-055 (DP-F6-16).
- Estimación de recursos por duplicaciones y presupuestos cuantitativos:
  REQ-IND-51 y DTPD-115, F9.
- Distribución genotípica explícita y filas varadas de la renumeración:
  sección 4 (DP-F6-10/11).
