# Revisión adversarial del plan de implementación de F2 0.1

> **Informe histórico.** Revisa una versión sustituida el 21 de agosto de 2026
> por el plan F2 0.2 y DTC-079; no impone sus ampliaciones de amenaza o
> infraestructura al plan vigente.

## 1. Cabecera, alcance y método

- **Documento revisado:**
  [plan-implementacion-f2-0.1.md](./plan-implementacion-f2-0.1.md).
- **Fecha de revisión:** 15 de agosto de 2026.
- **Estado del informe:** revisión independiente posterior a DTC-058 y a la
  ampliación del plan a F2.0–F2.5, con seguimiento documental de las
  correcciones aprobado el 15 de agosto de 2026.
- **Veredicto:** **no aprobar todavía el plan para ejecución**.
- **Cambios realizados durante la revisión inicial:** únicamente este informe;
  el seguimiento posterior sí corrige plan, especificación, arquitectura,
  decisiones, hoja de ruta y trazabilidad. No se modifica código ni se inicia
  F2.

Durante la revisión inicial se contrastaron completamente los 96 pasos y las
seis pausas de aquella versión del plan con
las secciones aplicables de:

- [especificaciones-0.1.md](./especificaciones-0.1.md), incluidos los contratos
  de configuración, individuos, matriz SNV, HWE/`F`, verdad, VCF, resultados y
  Marimo;
- [arquitectura-0.1.md](./arquitectura-0.1.md), en particular dominio,
  componentes, admisión, proceso UI, fuentes poblacionales, publicación,
  observabilidad y comprobación;
- [decisiones-tecnicas-0.1.md](./decisiones-tecnicas-0.1.md), incluidas
  DTC-049, DTC-052, DTC-054, DTC-056, DTC-058 y DTPD-102/103;
- [hoja-de-ruta-0.1.md](./hoja-de-ruta-0.1.md), la matriz completa de
  [trazabilidad](./trazabilidad-requisitos-0.1.md), el
  [registro de ampliaciones futuras](./registro-ampliaciones-futuras-0.1.md),
  la revisión adversarial anterior y el punto de reanudación;
- el estado real de `Cargo.toml`, `Cargo.lock`, `src/`, `tests/`, esquemas,
  catálogos, CI y `distribution/f1/`.

La inspección mecánica de la versión originalmente auditada confirmó 96 casillas, una línea `Depende de`, `Procede
de`, `Cambia` y `Se comprueba con` en cada una, seis pausas, 84 REQ exclusivos
y cuatro compartidos. Esa regularidad no hace ejecutable el plan: persisten
decisiones abiertas, dependencias equivocadas y comprobaciones sin comando ni
oráculo materializable.

## 2. Veredicto ejecutivo

La versión auditada había mejorado sustancialmente respecto a la revisión del
10 de agosto: incorporaba F2.0, separaba matriz y genotipos, mantenía
verdad→VCF, exigía pausas y proponía numerosos oráculos independientes. No
obstante, todavía obligaba al implementador a tomar decisiones durante la
ejecución.

En la versión originalmente auditada, los bloqueos se agrupaban así:

- **funcionales:** cierres de REQ mayores que el alcance de F2, distribución
  genotípica explícita no resuelta y evidencia positiva de complemento VCF que
  puede ser imposible sin inventar pérdida;
- **técnicos:** F1 no está cerrada, DTPD-102/103 siguen abiertas, Poisson e
  *indels* carecen de decisión canónica, el protocolo UI no puede expresar sus
  propias confirmaciones y faltan contratos ejecutables de importación,
  aleatoriedad, datos 1000 Genomes y navegador;
- **documentales:** numeración contraria a la convención vigente, pasos
  excesivos, casi ninguna comprobación con comando exacto y contadores
  desactualizados.

Tras el seguimiento, las decisiones y correcciones documentales están
integradas en 107 pasos jerárquicos. Para autorizar el inicio de F2 todavía debe
cerrarse y aprobarse F1 y una nueva revisión adversarial del plan corregido debe
confirmar que no quedan bloqueantes. F2.0.2, F2.2.1 y la regeneración/revisión
compacta DTC-068 de F2.2.14 son puertas internas de la ejecución posterior, no
prerrequisitos para iniciarla.

## 3. Hallazgos bloqueantes funcionales

### RAF2-004 — Cuatro REQ se siguen declarando cerrables con un subconjunto de su alcance

**Prioridad:** bloqueante funcional. **Estado:** resuelto el 15 de agosto de
2026.

**Resolución:** REQ-CONF-08/21 pasan a F2.2/F5.2/F7.2 y REQ-IND-02/03 a
F2.2/F5.2. F2 conserva expedientes parciales explícitos; §15 cuenta 79 cierres
exclusivos y nueve compartidos, sin ampliar su alcance.

La matriz y §15 asignan a F2.2:

- REQ-CONF-08, que incluye *indels*, variantes estructurales y mutaciones *de
  novo* posteriores;
- REQ-CONF-21, que exige el modelo espacial por separado para *indels*,
  variantes estructurales y mutaciones *de novo*;
- REQ-IND-02/03, que exigen escoger y conservar cualquiera de varias referencias
  disponibles.

F2 excluye variantes estructurales, mutación evolutiva y varias referencias.
B35 reconoce la discrepancia, pero continúa contando los cuatro entre sus 26
REQ exclusivos y §16 permite aprobar el plan con ellos en esa matriz. Una nota
de cautela no impide un cierre falso.

**Corrección exigida:** reasignar el cierre completo a la primera fase que cubra
todo el REQ y registrar en F2 solo evidencia parcial. No ampliar F2 para hacer
encajar la matriz.

**Condición de cierre:** los cuatro dejan de formar parte de los 84 cierres
exclusivos de F2 y aparecen, si procede, en una tabla explícita de evidencia
parcial con destino final.

### RAF2-007 — La hoja de ruta sigue exigiendo un complemento positivo no demostrado en F2

**Prioridad:** bloqueante funcional. **Estado:** resuelto el 15 de agosto de
2026.

**Resolución:** F2 prueba positivamente la ausencia de complemento para todos
sus casos representables y nunca fabrica pérdida; F5.5 aporta el primer caso
positivo complejo/SV/BND real conforme a DTC-065/066.

E7 aplica correctamente REQ-VER-24 y prohíbe crear complemento para un caso
plenamente representable. Sin embargo, la evidencia de cierre de F2 en la hoja
de ruta continúa exigiendo «un caso VCF/complemento», mientras el alcance solo
incluye SNP e *indels* inequívocos, GT faseado y ploidía 1–10. DTPD-103 aún no
ha demostrado una pérdida real en ese conjunto.

**Corrección exigida:** si DTPD-103 no encuentra un caso literal de F2 con
pérdida, cambiar la evidencia F2 por una prueba negativa de ausencia de
complemento y trasladar el primer caso positivo a F5.5. No crear una entrada
redundante para cumplir la hoja de ruta.

**Condición de cierre:** hoja de ruta, plan, trazabilidad y DTPD-103 piden el
mismo caso observable.

### RAF2-012 — REQ-CONF-130 admite una alternativa que C1 rechaza sin destino explícito

**Prioridad:** bloqueante funcional. **Estado:** resuelto el 15 de agosto de
2026.

**Resolución:** F2 implementa HWE P=1..10 y F solo diploide y registra
REQ-CONF-130 como parcial; la distribución genotípica explícita y el cierre
completo pasan a F6.2 con REQ-EVO-35/36/49.

REQ-CONF-130 dispone que en haploides y poliploides se use HWE **o una
distribución genotípica explícita ya admitida**. C1 limita F2.3 a HWE o `F`
diploide y §15 declara REQ-CONF-130 cerrado exclusivamente en F2.3. El plan no
implementa la distribución explícita, no la identifica como evidencia parcial
y tampoco traslada esa cláusula a otra fase.

**Corrección exigida:** decidir entre incorporarla en F2.3 o mantener el alcance
pequeño y reasignar el cierre completo del REQ, dejando en F2 evidencia parcial
de HWE/`F`.

**Condición de cierre:** la prueba asignada cubre todas las opciones normativas
del REQ o deja de presentarlo como cerrado.

## 4. Hallazgos bloqueantes técnicos

### RAF2-001 — Tres decisiones centrales siguen dentro del plan en lugar de precederlo

**Prioridad:** bloqueante técnico. **Estado:** resuelto documentalmente el 15 de
agosto de 2026.

**Resolución:** DTC-059 y DTC-061 a DTC-068 preceden a sus consumidores; los
pasos F2.1.3, F2.2.2 y F2.5.1 congelan/auditan contratos y ya no los eligen.

En la versión auditada, DTPD-102 continuaba pendiente y A3 pretendía elegir a la vez la representación
interna y el formato obligatorio de verdad. DTPD-103 continuaba pendiente y E1
elegía versión VCF, normalización, campos, complemento y validadores. B2 todavía
debía crear la decisión de implementación Poisson y registrar el contrato ya
acordado de longitudes exactas de *indels*.

Esto contradice la puerta de entrada del proceso de planificación: un plan
ejecutable no puede usar pasos de implementación para decidir sus contratos.
Además, A3 prueba solo copias, SNP, *indels* y `N`, aunque la propia DTPD-102
exige evaluar composición con variantes estructurales para no elegir una
representación incapaz de sostener F5.

**Corrección exigida:** resolver y aprobar antes del plan ejecutable:

1. estructura interna haplotípica y formato externo como decisiones distintas,
   incluyendo un prototipo estructural mínimo;
2. contrato VCF/complemento completo;
3. algoritmo Poisson, consumo RNG, límites y oráculo;
4. parámetros exactos de longitud de inserción/deleción de F2.

A3, B2 y E1 deben convertirse después en comprobaciones de decisiones ya
aprobadas, no en lugares donde se escogen.

### RAF2-013 — F1 sigue abierta mientras su corrección espera revisión independiente

**Prioridad:** bloqueante técnico. **Estado:** nuevo.

El árbol real no coincide con la fotografía de §3.1 del plan. F1.5 está
implementada. Durante esta auditoría H-F16-R13 fue corregido mediante cruce
global exacto y un corpus con cero, una y dos reducciones; el expediente vuelve
a declarar 22/22 REQ de F1.6 y 83/83 de F1. Sin embargo, la revisión
independiente final sigue en curso y ni la pausa F1.6 ni F1 están aprobadas.

F2 reutiliza precisamente integridad, manifiestos y publicación. Empezar M1 o A1
sobre esa base convertiría un defecto conocido en contrato de la interfaz y de
las nuevas salidas.

**Condición de cierre:** revisión independiente de H-F16-R13 verde, pausa F1.6
y F1 aprobadas y commit base congelado.

### RAF2-014 — El protocolo fijo de Marimo no puede validar sus dos autorizaciones limitadas

**Prioridad:** bloqueante técnico/funcional. **Estado:** resuelto
documentalmente el 15 de agosto de 2026.

**Resolución:** las autorizaciones forman parte de la invocación y del digest;
cada confirmación obliga a revalidar y `run` solo admite la instantánea y las
autorizaciones exactamente validadas.

M8 fija exactamente:

```text
validate <general> --diagnostics json
```

y M9 fija:

```text
run <general> --diagnostics json --output <destino>
```

Pero B21 y C9 requieren `--allow-empty-snv-matrix` y
`--allow-incompatible-genotypes`. Sin autorización, una contradicción debe
bloquear antes de crear individuos; con autorización debe volver a validarse,
quedar en `EffectiveConfig`/manifiesto y solo entonces ejecutarse. C17 afirma
que la confirmación produce flags, pero M8/M9 no permiten pasarlos ni asocian
las autorizaciones al digest validado.

**Riesgo:** la UI no puede mostrar el conflicto para confirmarlo, ejecuta una
configuración distinta de la validada o aplica una autorización antigua a un
documento editado.

**Corrección exigida:** versionar una unión cerrada de opciones de invocación
admitidas en `validate` y `run`; incluir flags y valores en el digest de la
instantánea; repetir `validate` tras confirmar; exigir que `run` use exactamente
documentos, contrato y autorizaciones validados. El proceso falso debe rechazar
cualquier flag general o no catalogado.

### RAF2-015 — La tolerancia multialélica de B19 contradice DTC-049

**Prioridad:** bloqueante técnico con resultado biológico silencioso.
**Estado:** resuelto documentalmente el 15 de agosto de 2026.

**Resolución:** F2.2 consume y registra la tolerancia efectiva de DTC-049 y
prueba el valor cero, el predeterminado, un valor sustituido y sus fronteras.

B19 acepta y normaliza solo un exceso `<=1e-9`. DTC-049 establece una tolerancia
general **efectiva y configurable** en `[0,1e-6]`, con `1e-9` únicamente como
predeterminado. REQ-IND-60 remite a esa tolerancia general.

**Riesgo:** dos ejecuciones con tolerancias efectivas válidas distintas serían
admitidas o rechazadas como si ambas usaran el valor predeterminado.

**Corrección exigida:** B19 debe consumir la tolerancia efectiva ya resuelta,
registrar antes/después y probar 0, el predeterminado, un valor sustituido y sus
dos fronteras inclusiva/excesiva.

### RAF2-016 — La generación de *indels* reutiliza dependencias de la matriz sin definir sus propios ámbitos

**Prioridad:** bloqueante técnico/genómico. **Estado:** resuelto
documentalmente el 15 de agosto de 2026.

**Resolución:** F2.2 define ámbitos RNG independientes de recuento, posición y
contenido para inserciones y deleciones, y separa sus oportunidades y medias
del generador de matrices.

B8 solo versiona `matrix/locus_selection`, `matrix/alternate_allele` y
`matrix/frequency` y dice que quedan separados de *indels*. No existe después
un paso que publique ámbitos para recuento, posición y bases insertadas por
clase, individuo y región. B28/B29 dependen de B9–B11, aunque B10/B11 definen el
selector ponderado y el recuento de **loci candidatos de la matriz**, no el
modelo espacial por cantidad/densidad individual de REQ-CONF-20/21/22 y
REQ-IND-13/15/18.

**Riesgo:** cambiar una matriz puede alterar *indels*, una tasa cero consumir
azar, o `weighted_regions` sustituir por accidente el modelo regional de
*indels*. También falta separar la secuencia aleatoria insertada de su posición.

**Corrección exigida:** pasos propios y pequeños para cada clase:

- catálogo de ámbitos de recuento, selección espacial y contenido insertado;
- enumeración de oportunidades válidas sobre el estado vigente;
- media por ámbito y Poisson local;
- selección secuencial sin reintentos y actualización del espacio;
- aislamiento frente a matriz, individuo y demás clase.

### RAF2-017 — M7 no define una importación multidocumento segura ni conservadora

**Prioridad:** bloqueante técnico de F2.0. **Estado:** resuelto el 15 de agosto
de 2026 mediante DTC-069.

**Resolución:** la unidad es un directorio autocontenido; importación y
exportación son atómicas, aisladas, resistentes a TOCTOU y limitadas, y
preservan opacamente lo no editado sin ZIP/TAR ni selección múltiple.

M7 promete importar exactamente un general y cero o varios escenarios, pero no
define si la unidad de entrada es directorio, selección múltiple o archivo; cómo
se preserva la jerarquía relativa; cómo se rechazan `..`, rutas absolutas,
*symlinks* y cambios TOCTOU; ni cómo se copian los bytes a la instantánea
privada. Tampoco exige conservar como contenido opaco los campos de escenarios
Evolución/Familias que F2.0 todavía no sabe editar. «Canonicalizar rutas» puede
convertir una ruta relativa portable en una absoluta de la máquina local.

**Riesgo:** pérdida silenciosa de campos al importar/exportar, enlaces que salen
del conjunto, semántica distinta entre previsualización y binario o acceso a
ficheros no seleccionados.

**Corrección exigida:** fijar la unidad de importación y su modelo de amenaza;
separar literal TOML de ruta resuelta mostrada; copiar de forma segura todos los
documentos enlazados; conservar opacos los no editables; y validar el conjunto
copiado, no rutas mutables del origen.

### RAF2-018 — El perfil 1000 Genomes no tiene todavía una receta ejecutable ni una dependencia correcta

**Prioridad:** bloqueante técnico de B14/B15. **Estado:** resuelto el 15 de
agosto de 2026 mediante DTC-068 y simplificado después a dos perfiles compactos
aproximados, con fuente/receta/hash pequeños y sin VCF ni `bcftools`.

En la versión auditada, DTC-054 fijaba el nombre lógico y filtros generales, pero no el ensamblado exacto,
URL/artefactos de origen, SHA-256 de entradas, campos por superpoblación,
tratamiento de cromosomas/contigs, descomposición multialélica, herramientas ni
comando reproducible. B14 exige regeneración independiente sin que esas entradas
existan. B15 creaba el antiguo perfil SNV12 desde el mismo corpus pero no
dependía de B14 ni de una evaluación de datos común.

**Resolución vigente:** F2.2.14 regenera
`human_1000g_global_sfs_approx_v1` desde 64:12:8 y `1/5008` y
`human_harris2017_snv12_approx_v1` desde el fichero eLife de 3.644 bytes y
SHA-256 fijado. Un oráculo independiente comprueba normalizaciones,
transformación 96→6→12, licencia CC BY y digests. La materialización exacta del
VCF queda en DTF-016/post-v1 y no bloquea F2.

### RAF2-019 — `prepare-snv-matrix` no puede construir los metadatos que B4 exige

**Prioridad:** bloqueante técnico. **Estado:** resuelto documentalmente el 15
de agosto de 2026.

**Resolución:** `prepare-snv-matrix` obtiene `reference_id` y la huella desde el
manifiesto F1 autoritativo asociado al FASTA y contrasta ambos antes de escribir.

B17 propone:

```text
prepare-snv-matrix --reference <FASTA> --input <TSV> --output <destino>
```

B4 exige además el identificador estable de la referencia, no deducible de una
ruta FASTA. Tampoco se indica si se consume un manifiesto F1, un descriptor de
referencia o un `--reference-id` explícito ni cómo se contrasta ese ID con el
FASTA autoritativo.

**Corrección exigida:** fijar una fuente inequívoca de `reference_id` y probar
ID correcto, ID distinto con mismos bytes, mismos ID con bytes distintos y
referencia derivada.

### RAF2-020 — El entorno bloqueado no incluye los binarios que ejecuta Playwright

**Prioridad:** bloqueante técnico de M1/M11. **Estado:** resuelto
documentalmente el 15 de agosto de 2026.

**Resolución:** la puerta F2.0 fija navegador y revisión, inventario y hashes de
artefactos, dependencias de Ubuntu y una ejecución posterior sin red; separa la
evidencia automatizada de la revisión manual WSL2.

`uv.lock` puede fijar los paquetes Python `playwright` y la integración de axe,
pero no fija por sí solo los binarios de Chromium/Firefox/WebKit ni las
bibliotecas del sistema que descarga o instala Playwright. `uv sync --frozen`
no demuestra una prueba de navegador repetible ni un arranque posterior sin red.

**Corrección exigida:** elegir y registrar navegador, revisión, artefacto/hash,
comando de instalación y dependencias Ubuntu; separar pruebas automatizadas de
la revisión manual WSL2; ejecutar ambas sin descarga después de la instalación.

### RAF2-021 — M1 mezcla evaluación, implementación y aceptación en una sola casilla

**Prioridad:** bloqueante estructural del plan. **Estado:** resuelto
documentalmente el 15 de agosto de 2026.

**Resolución:** la antigua M1 se dividió en una puerta experimental previa y
pasos independientes de entorno, empaquetado, lanzador, red, OCI y evidencias de
plataforma.

M1 evalúa tres herramientas, crea entorno y *lock*, incorpora aplicación,
lanzador, contratos y paquete, modifica la entrega `.tar.gz`, comprueba OCI,
prueba red, puertos, navegador, proceso falso, binario real, Ubuntu, WSL2,
licencias y una vulnerabilidad concreta. No es un cambio independiente ni puede
revisarse o revertirse como una unidad. Además, contratos y protocolo se
implementan formalmente en M2/M8, por lo que parte de la aceptación de M1 carece
todavía de productor.

**Corrección exigida:** dividir la puerta experimental de versiones/licencias,
el esqueleto de empaquetado, el lanzador/red, el entorno sin red, OCI y las
evidencias de plataforma. M1 debe terminar en una decisión fechada antes de
incorporar dependencias de producción; la implementación depende de esa puerta.

## 5. Hallazgos mayores no bloqueantes por sí solos

### RAF2-022 — La numeración no cumple la convención jerárquica vigente

**Estado:** resuelto documentalmente el 15 de agosto de 2026.

Los identificadores M1/A1/B1/C1/D1/E1 no heredan F2.0–F2.5. La habilidad vigente
y la directriz incorporada al proyecto exigen `F2.0.1`, `F2.1.1`, etc. El plan
no se ha ejecutado, por lo que aún puede renumerarse sin romper evidencia
histórica. También deben cambiar dependencias, pausas, matrices y comandos de
validación.

### RAF2-023 — Las comprobaciones describen pruebas, pero casi nunca dan comandos ejecutables

**Estado:** resuelto documentalmente el 15 de agosto de 2026.

Las 96 casillas tienen texto `Se comprueba con`, pero solo una fracción mínima
incluye un comando. No se fijan nombres de pruebas, rutas de fixtures, comandos
de validadores, artefactos de evidencia ni resultado esperado material. §4.3
solo aporta puertas globales; no permite demostrar cada paso ni reanudarlo con
trazabilidad.

Debe existir por paso al menos un comando real o un enlace a un script/fixture
planificado, más el fichero de evidencia y el resultado esperado. Los
validadores externos deben llevar versión y comando exactos después de sus
decisiones.

### RAF2-024 — El supervisor de procesos no exige drenaje concurrente ni prueba de bloqueo

**Estado:** resuelto documentalmente el 15 de agosto de 2026.

M8–M10 limitan `stdout`/`stderr`, pero no exigen leer ambos concurrentemente.
Un proceso que llena una tubería mientras la UI espera la otra puede bloquearse
antes de alcanzar el límite. Tampoco se exige cerrar `stdin`, reaprovechar el
hijo ni comprobar que las señales solo alcanzan el grupo creado por la UI.

Añadir pruebas falsas que llenen ambos flujos intercalados, mantengan
descendientes, ignoren cada señal y comprueben cierre, reaprovechamiento y
ausencia de proceso huérfano.

### RAF2-025 — La trazabilidad omite evidencia parcial directa de REQ-SAL-02 y REQ-SAL-12

**Estado:** resuelto documentalmente el 15 de agosto de 2026.

F2.4 implementa para creación directa el fichero obligatorio de REQ-SAL-02 y
F2.5 implementa BGZF/índice de REQ-SAL-12, pero §15.2 solo los menciona como
«evidencia adyacente» sin enlazarlos a D2/D8 y E9/E12 como ocurre con los cuatro
REQ compartidos. No deben cerrarse en F2, pero sí dejar una evidencia parcial
ejecutable y un destino final explícito.

### RAF2-026 — El estado de fuentes y contadores documentales está desactualizado

**Estado:** resuelto documentalmente el 15 de agosto de 2026.

- §3.1 del plan afirma F1.5 pendiente, aunque está implementada y publicada;
  lo pendiente es la revisión y aprobación final de F1.6/F1.
- La cabecera de trazabilidad declara 843 filas, mientras el fichero contiene
  844 filas `REQ`.
- La revisión adversarial anterior aún describe 52 pasos y cinco pausas; este
  informe la sustituye para el estado actual de 96 y seis.
- La arquitectura conserva una frase que dice que la representación interna de
  coordenadas «se decidirá por separado» después de haber fijado DTC-056.

Son correcciones documentales, pero deben resolverse antes de congelar la línea
base para que las puertas mecánicas no acrediten una fotografía falsa.

## 6. Preguntas y decisiones que sí requieren intervención del usuario

Las seis decisiones quedaron resueltas el 15 de agosto de 2026:

1. **DTPD-102:** modelo híbrido y `truth.jsonl` mediante DTC-061.
2. **DTPD-103:** contrato VCF completo mediante DTC-062 a DTC-066.
3. **RAF2-004:** cuatro REQ amplios pasan a evidencia parcial F2 y cierre
   posterior, con 79 cierres exclusivos y nueve compartidos.
4. **RAF2-007:** F2 prueba ausencia de complemento redundante y F5.5 aporta el
   primer positivo real.
5. **REQ-CONF-130/RAF2-012:** HWE/F en F2 y distribución explícita en F6.2.
6. **Importación Marimo/RAF2-017:** directorio autocontenido seguro según
   DTC-069.

## 7. Aspectos correctos y reutilizables

Se deben conservar al corregir el plan:

- la separación matriz común → genotipo/fase por individuo → verdad externa →
  VCF;
- la prohibición de usar Poisson para frecuencias, alelos o genotipos;
- candidatos, `ALT` globalmente descartados y loci activos como cantidades
  distintas, sin reposición;
- `alt_af` en `[0,1]`, alelos fijados conservados, poblaciones independientes y
  perfiles humanos aproximados DTC-068 solo por selección explícita y perfiles
  suministrados distintos por población;
- IDs canónicos, subsemillas por identidad y pruebas de estabilidad al añadir o
  reordenar entidades;
- coordenadas mediante el servicio común DTC-056 y casos de borde de inserción
  y deleción;
- conservación literal de `N` y separación de inserción explícita/aleatoria;
- oráculos que no llaman al productor y prueba de verdad aun sin VCF;
- VCF derivado de la verdad, GT faseado, dosis real y ausencia de campos
  observacionales inventados;
- publicación no-clobber, hashes, fallos sin salida parcial y pausa al final de
  cada hito;
- Marimo sin autoridad biológica, sin *shell*, sin ejecución automática, con
  digest invalidado tras editar y sin cargar artefactos grandes;
- enlace exclusivo a `127.0.0.1`, colisión de puerto explícita y separación
  `.tar.gz` con UI frente a OCI solo CLI.

La inspección del repositorio confirma nombres reutilizables reales:
`genome-population-simulator`, `validate`, `run`, `--diagnostics`, `--output`,
`DeclaredConfig`, `EffectiveConfig`, el resolvedor, los tipos de coordenadas,
SHA-256, publicación y manifiestos. No existe todavía `capabilities`, Marimo,
`pyproject.toml`, `uv.lock`, los dos flags de autorización, matriz SNV, verdad
genotípica ni VCF; el plan debe presentarlos como nuevas capacidades, no como
componentes disponibles.

## 8. Orden de corrección recomendado

1. **Pendiente externo a este informe:** revisar F1.6/F1 y congelar el commit base.
2. **Completado documentalmente:** decisiones RAF2-004/007/012/017 y
   DTPD-102/103/120/121 integradas.
3. **Completado documentalmente:** Poisson e *indels* canónicos.
4. **Completado documentalmente:** autorizaciones UI, importación, receta
   compacta DTC-068 y navegador bloqueado.
5. **Completado documentalmente:** plan jerárquico de 107 pasos.
6. **Completado documentalmente:** 79 cierres, nueve compartidos y dos
   expedientes parciales directos.
7. **Pendiente antes de autorizar el inicio:** cerrar/aprobar F1 y repetir la
   revisión adversarial sobre el plan corregido; debe concluir sin bloqueantes.
   F2.0.2, F2.2.1 y F2.2.14 producirán sus evidencias después, durante F2, y
   bloquearán solo a sus consumidores.

## 9. Comandos de revalidación

Comprobaciones documentales mínimas tras las correcciones:

```bash
git diff --check -- design/plan-implementacion-f2-0.1.md \
  design/especificaciones-0.1.md design/arquitectura-0.1.md \
  design/decisiones-tecnicas-0.1.md design/hoja-de-ruta-0.1.md \
  design/trazabilidad-requisitos-0.1.md

rg -n '^- ☐ F2\.[0-5]\.[0-9]+\.' design/plan-implementacion-f2-0.1.md
rg -n '^> \*\*Pausa obligatoria F2\.[0-5]' \
  design/plan-implementacion-f2-0.1.md
rg -n 'DTPD-102|DTPD-103|RAF2-004|RAF2-007|REQ-CONF-130' \
  design/plan-implementacion-f2-0.1.md \
  design/decisiones-tecnicas-0.1.md design/hoja-de-ruta-0.1.md
rg '\| F2(\.[0-5])?(/| \||$)' design/trazabilidad-requisitos-0.1.md | wc -l
rg '^\| REQ-' design/trazabilidad-requisitos-0.1.md | wc -l
```

Puertas de línea base antes de iniciar F2:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --workspace --all-targets --locked
cargo test --workspace --doc --locked
cargo build --release --locked --target x86_64-unknown-linux-gnu
cargo metadata --locked --format-version 1
cargo tree --locked -e features
git status --short
```

Puertas internas que el plan corregido ejecutará después de iniciar F2 y antes
de habilitar sus respectivos consumidores:

```text
F2.0.2: uv sync --frozen con CPython exacto, arranque Marimo y navegador bloqueado
F2.2.1: auditoría de consumo de DTC-054 sin reabrir la decisión
F2.2.14: regeneración compacta DTC-068 y revisión de fuente, receta y SHA-256
F2.4.1/F2.4.4: validadores formal y semántico de truth.jsonl
F2.5.1/F2.5.2/F2.5.10: bcftools/HTSlib/tabix 1.24 y oráculo VCF independiente
```

No se han ejecutado las puertas Cargo durante esta auditoría porque otro agente
está corrigiendo F1.6 en el mismo árbol compartido. Esta revisión no atribuye ni
certifica esos cambios concurrentes.
