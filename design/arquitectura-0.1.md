# Arquitectura del simulador de genomas y poblaciones

- **Versión:** 0.2
- **Estado:** F1 cerrado en `e14705a`/`efca265`; arquitectura F2 simplificada
  para ejecución local confiable y pendiente de aprobación
- **Fecha:** 21 de agosto de 2026
- **Especificación de origen:**
  [especificaciones-0.1.md](./especificaciones-0.1.md)
- **Hoja de ruta:** [hoja-de-ruta-0.1.md](./hoja-de-ruta-0.1.md)
- **Registro técnico:**
  [decisiones-tecnicas-0.1.md](./decisiones-tecnicas-0.1.md)
- **Lenguaje del núcleo:** Rust, según la especificación y la hoja de ruta 0.1

Este documento es derivado: organiza técnicamente los requisitos normativos sin
ampliarlos ni sustituirlos. La cobertura completa se consulta en
[trazabilidad-requisitos-0.1.md](./trazabilidad-requisitos-0.1.md).
Las capacidades post-v1 se remiten al
[registro de ampliaciones futuras](./registro-ampliaciones-futuras-0.1.md);
las elecciones DTP/DTPD siguen siendo planificación técnica, no funciones.

## 1. Finalidad y estado

**Corrección F2 posterior a F1.** La base cerrada es `e14705a` con cierre
documental `efca265`. F1 no proporciona vía, escenarios, ploidía, reglas ni
dominios poblacionales. F2 añade únicamente los conceptos que usa su recorrido
directo: ploidía, locus, población de origen, individuo, haplotipo y variantes
pequeñas. La vía directa usa el ID interno reservado
`scenario_id = "direct"` solo para identidad y RNG; F2 no crea documentos ni
infraestructura de escenarios. Marimo edita un TOML general, llama a
`validate`/`run` y presenta sus resultados. La instalación y el arranque se
prueban desde el repositorio; `.tar.gz`, OCI y acreditación oficial de
plataforma quedan en F9.

F2 hereda el modelo de confianza local mínimo de F1: valida entradas externas
y protege frente a fallos propios y corrupción accidental, pero confía en la
persona usuaria, el repositorio y el provisional privado frente a manipulación
activa. No construye defensas contra symlinks hostiles, inyección deliberada del
entorno, carreras adversarias ni sustitución concurrente de ficheros. Tampoco
crea un ejecutor de evidencia general: las pruebas normales de cada componente
y las evidencias de cierre del hito son suficientes.

Este documento define cómo se divide el sistema y qué contratos deben mantener
sus partes para satisfacer la especificación 0.1. No define archivos de código,
paquetes, algoritmos concretos ni un plan de implementación salvo las elecciones
técnicas ya confirmadas. TOML está confirmado como formato de configuración
de autoría, y JSON como formato de los manifiestos y de la configuración
efectiva serializada. `noodles-fasta` 0.66.0 está confirmado y fijado para el
adaptador FASTA tras superar la evaluación revisable de F1. Las versiones
concretas de los crates aceptados quedan fijadas por `Cargo.lock`.

La ampliación confirmada incorpora una derivación reducida de FASTA en F1.6 y
una aplicación Marimo local desde F2.0. DTC-056/DTC-057 cierran las coordenadas
comunes y DTC-079 simplifica la frontera Marimo; la distribución oficial se
pospone a F9.

La arquitectura cubre toda la versión 1 para evitar contratos incompatibles
entre fases. Sin embargo, su materialización seguirá la secuencia F1–F9 de la
hoja de ruta. Una capacidad futura no autoriza a implementarla antes de su fase.

La especificación y la hoja de ruta confirman Rust para el núcleo. Esta elección
no decide por sí sola formatos, *crates* ni algoritmos. DTC-048 fija el
toolchain inicial y DTC-042, la plataforma oficial y la distribución.

## 2. Contexto, alcance y límites

### 2.1. Responsabilidad del sistema

El sistema recibe una configuración, referencias y entradas opcionales; crea
una verdad biológica simulada; puede producir evidencia técnica de
secuenciación y alineamiento; valida la coherencia de los resultados; y conserva
la información necesaria para reconstruir qué ocurrió. Satisface principalmente
las secciones 2, 4, 5 y 6 de la especificación.

Las tres vías biológicas comparten el mismo modelo de referencia, individuo,
haplotipo, variante y verdad:

1. creación directa de individuos;
2. población inicial y evolución por generaciones;
3. familias y poblaciones segregantes definidas mediante un grafo de cruces.

Toda salida asociada a una observación se organiza mediante un evento explícito.
Una ejecución puede declarar cero eventos y, en ese caso, no se crea ninguno de
forma implícita. Todo evento tiene identidad, selección de individuos,
referencias aplicables, salidas y subsemilla. Solo los eventos que producen
FASTQ, BAM o CRAM constituyen su único experimento de secuenciación y requieren
tecnología, perfil, exactamente una tanda y una biblioteca por individuo; solo
los que producen alineamientos requieren además un mapeador. Todas las
bibliotecas se procesan juntas dentro de la tanda. Por tanto, un evento puede
materializar únicamente verdad genotípica o VCF sin parámetros ni entidades
técnicas de secuenciación.

### 2.2. Fronteras externas

Son fronteras del sistema:

- la aplicación Marimo local, que crea conjuntos TOML e invoca el binario como proceso
  sin shell, pero no interpreta reglas biológicas;
- los ficheros TOML de configuración general y poblacional;
- `effective-config.json`, `manifest.json` general y los `manifest.json` de
  evento, todos generados por el sistema;
- FASTA de referencia de texto sin comprimir; la compresión FASTA no se admite
  en v1;
- FASTQ suministrados para contaminación;
- los estándares FASTA, FASTQ, SAM/BAM, CRAM y VCF;
- minimap2 como mapeador predeterminado y los mapeadores compatibles que se
  autoricen en el futuro;
- las herramientas externas necesarias para mapear, comprimir, ordenar,
  indexar o validar formatos, cuando la decisión técnica correspondiente las
  seleccione;
- el sistema de ficheros y los recursos de CPU, memoria y almacenamiento.
- las entregas oficiales para Ubuntu 24.04 LTS x86-64 GNU, con triple
  `x86_64-unknown-linux-gnu`: la construcción y validación oficial del paquete
  `.tar.gz` y de la imagen OCI se posponen a F9.

Cada integración externa se realiza mediante un contrato sustituible. El
dominio biológico no conoce comandos, rutas concretas, cabeceras físicas ni
detalles de procesos externos.

El adaptador Rust puro `noodles-fasta` cubre FASTA, FAI y acceso indexado detrás
del contrato propio de referencia. `samtools faidx` solo actúa como validador
independiente en pruebas y no es una dependencia normal de ejecución.

### 2.3. Fuera de alcance arquitectónico

No forman parte del sistema los *SNP callers*, la comparación de sus resultados,
la imputación, el faseado inferido ni la evaluación de herramientas externas.
Tampoco se incluyen las capacidades funcionales excluidas en la sección 3 de la
especificación, entre ellas señales de instrumento, BCL, paneles dirigidos,
Oxford Nanopore, PacBio CLR, aneuploidía individual, reproducción clonal,
conversión génica, calibración empírica automática y todo modelado específico de
elementos transponibles. V1 no define una frontera de entrada para FASTQ y
referencia de calibración ni para perfiles empíricos reutilizables.

Python no forma parte del núcleo. Se usa solo para Marimo y queda fuera del
dominio biológico. Sus dependencias se declaran en `pyproject.toml`, se fijan en
`uv.lock` y se instalan en un entorno local excluido de Git. Marimo no instala
paquetes ni consulta la red durante el uso normal. El aprovisionamiento inicial
explícito mediante `uv sync --frozen` puede acceder a Internet para descargar
exclusivamente las dependencias resueltas por `uv.lock`; después se prueba el
arranque local sin red. Las versiones se aceptan cuando el entorno bloqueado
instala, arranca y ejecuta sus pruebas; F2 no añade una evaluación de
distribución o seguridad separada para cada dependencia. La instalación
distribuible sin descarga previa pertenece a F9.

Marimo es un componente auxiliar separado del binario Rust. La aplicación,
`pyproject.toml` y `uv.lock` se versionan en el repositorio. El entorno local no
se versiona ni se incrusta en Rust. El contenido de las entregas oficiales y de
la imagen OCI se decide y verifica en F9.

La aplicación Marimo escucha exclusivamente en `127.0.0.1`. El lanzador usa un
puerto local predeterminado, permite elegir otro puerto no privilegiado, muestra
la URL y puede omitir la apertura del navegador. Un puerto ocupado produce un
diagnóstico. El puerto no pertenece al TOML biológico. Abrir la aplicación no
ejecuta una simulación.

macOS, Windows nativo y ARM64 no son plataformas oficialmente soportadas en v1.
Su ampliación requiere comprobar compilación, dependencias, herramientas
externas y una matriz de pruebas. WSL2 es un entorno de desarrollo admitido, no
una plataforma de publicación adicional ni una fuente de dependencias
exclusivas. Otras distribuciones x86-64 GNU pueden funcionar, pero no se
consideran oficialmente validadas hasta incorporarlas a esa matriz.

## 3. Principios arquitectónicos

1. **La verdad biológica precede a la evidencia técnica.** Lecturas,
   alineamientos, contaminación, errores y calidades consultan la verdad, pero
   nunca la modifican.
2. **Cada dato lleva identidad y procedencia.** Referencia, sistema de
   coordenadas, individuo, haplotipo, evento y conjunto no se deducen de una
   ruta o de una posición implícita.
3. **Las referencias no se mezclan.** Un conjunto de resultados pertenece a
   una única referencia de mapeo y sus identificadores de variante son locales.
   En v1, una unidad reproductiva contiene únicamente individuos basados en la
   misma referencia estructural, aunque puedan observarse contra varias
   referencias de mapeo.
4. **La configuración se resuelve antes de simular.** El adaptador TOML conserva
   el documento original exacto, comprueba su clase y versión conocida y solo
   entonces construye `DeclaredConfig` con procedencia y ubicación. Cada fase
   incorpora al resolvedor únicamente los parámetros, selectores y relaciones
   que consume realmente y obtiene un `EffectiveConfig` inmutable para ese
   recorrido. F1 no simula dominios poblacionales, temporales o espaciales para
   anticipar la validación de fases posteriores.
5. **Los estados inválidos no cruzan fronteras.** Los identificadores,
   coordenadas, ploidías, proporciones y relaciones se representan mediante
   valores validados o producen errores explícitos.
6. **La aleatoriedad tiene ámbitos independientes.** Cada ámbito usa su propio
   ChaCha12, derivado mediante BLAKE3 y separación de dominios; no existe
   `StdRng`, flujo global ni generador mutable compartido. Población inicial,
   generaciones y procesos biológicos, muestreos y eventos técnicos no
   comparten consumo. La reproducibilidad se declara por separado para biología,
   contenido técnico, equivalencia semántica e identidad binaria; alcanzar un
   nivel no presupone alcanzar el siguiente.
7. **Los estados y la publicación crecen con sus consumidores.** F1 distingue
   ejecución y artefactos como `incompleto`, `completo` o `fallido`, trabaja en
   un provisional privado y publica mediante renombrado atómico no-clobber. Las
   jerarquías de evento, conjunto y exportación se incorporan en las fases que
   creen esas unidades, sin simularlas mediante manifiestos vacíos.
8. **Los recursos se observan antes de optimizar.** La arquitectura expone
   puntos de medida y admite representaciones compactas posteriores, sin fijar
   objetivos cuantitativos no acordados.
9. **Marimo es un adaptador de autoría y ejecución.** Produce el mismo conjunto TOML,
   consume versiones publicadas de esquemas y catálogos e invoca `validate` y
   `run`. Nunca calcula valores efectivos, resuelve reglas ni modifica la
   semántica biológica. Toda capacidad configurable F2–F9 añade su control antes
   de cerrar el hito que la introduce.
10. **El modelo de confianza es proporcional al despliegue.** Los TOML, FASTA y
    demás entradas externas son no confiables y se validan antes de entrar. Los
    ficheros creados por el proceso dentro de su provisional privado se confían
    frente a manipulación local activa y se protegen frente a errores propios,
    interrupciones y corrupción accidental mediante estados, tamaños, SHA-256
    y publicación atómica. La defensa exhaustiva contra symlink, hardlink,
    sustitución concurrente e inodos ajenos pertenece a F9/CI cuando el modelo
    de despliegue la requiera; no condiciona la reproducibilidad científica de
    F1 ni F2.

## 4. Dirección de dependencias

```text
Persona usuaria / aplicación Marimo / entradas y herramientas externas
             │
             ▼
Adaptadores de UI, configuración, formatos y procesos
             │
             ▼
Orquestación de ejecución y de eventos
             │
             ▼
Servicios de dominio biológico y técnico
             │
             ▼
Modelo de dominio e invariantes compartidas

Persistencia y manifiestos implementan contratos definidos hacia el dominio;
el dominio no depende de sus formatos físicos.
```

La dirección permitida es hacia el modelo de dominio. Los componentes de
secuenciación pueden consultar referencias y verdad genotípica. Los componentes
biológicos no pueden depender de FASTQ, BAM, CRAM, VCF, mapeadores ni perfiles
técnicos. La validación puede consultar ambos lados, pero no corregirlos.

Las dependencias externas solo aparecen en la capa de adaptadores. Cuando un
crate Rust puro cumple correctamente el estándar, es la opción preferida para
configuración, serialización, hashes, RNG y formatos. V1 no enlaza bibliotecas
nativas C/C++ salvo necesidad imprescindible, evaluación separada y aprobación
explícita. Las herramientas genómicas completas y consolidadas se ejecutan como
procesos externos; esta política no decide aún la invocación concreta de
minimap2.

Todo adaptador traduce tipos y errores de la dependencia al contrato propio del
sistema. Ninguna API pública o frontera de dominio expone estructuras, errores
o convenciones accidentales de un crate, una biblioteca o un proceso.

## 5. Modelo conceptual del dominio

### 5.1. Identidades

Cada identidad es estable y pertenece a una categoría. Identidades de
categorías distintas no son intercambiables aunque su texto coincida. Las
fronteras de dominio usan tipos distintos —por ejemplo, `IndividualId`,
`EventId` y `ReferenceId`— y nunca una cadena genérica. La etiqueta humana es
Unicode, opcional y no única; no participa en identidad, referencias ni
precedencia.

| Entidad | Identidad y ciclo de vida | Relaciones obligatorias |
|---|---|---|
| Ejecución | Se crea tras validar y resolver la configuración. No se reutiliza | Versión, entradas, configuración efectiva, semilla raíz, estados y manifiesto general |
| Referencia | Estable durante la ejecución | FASTA exacto, secuencias ordenadas, clasificación cromosoma/*scaffold*, origen y cambios si es derivada |
| Derivación reducida | Identidad propia dentro de la ejecución; no sustituye a su origen | Referencia/FASTA de origen, parámetros efectivos, subsemilla, versión de selección, intervalos ordenados, mapa e integridad |
| Secuencia de referencia | Su identificador procede del primer elemento no vacío de la cabecera FASTA y es único dentro de su referencia | Descripción restante de cabecera, longitud, clasificación, ploidía efectiva y convenio de coordenadas |
| Segmento | Identifica una porción orientada y trazable | Estado actual, origen fundador o evento que lo creó, orientación y correspondencias disponibles |
| Variante o evento genómico | Identidad estable en el historial biológico | Clase, orden, marco de coordenadas previo, haplotipo afectado y estado final |
| Espectro de frecuencias alélicas | Identidad y versión estables; puede ser incorporado o suministrado | Definición de frecuencia, valores o intervalos, pesos, población/fuente/filtros y transformación |
| Matriz de frecuencias SNV | Inmutable una vez admitida o generada | Referencia exacta, poblaciones, loci/alelos, frecuencias, perfil y procedencia generadora |
| Panel de bloques haplotípicos | Paquete inmutable de cinco ficheros durante la ejecución | `panel.metadata.json`, referencia exacta, cuatro TSV enlazados, bloques e intervalos no solapados, sitios SNP A/C/G/T, haplotipos de identidad local `block_id` + `haplotype_id`, asignaciones y frecuencias densas por población y bloque |
| Individuo | Único y no reutilizable en la ejecución | Referencia estructural, generación, población, parentales, haplotipos y observaciones |
| Copia haplotípica | Distinta aunque tenga la misma secuencia que otra | Individuo, secuencia, índice de copia, ploidía, segmentos y fundador |
| Población y subpoblación | Identidades no reutilizables | Generación, miembros, estado activa/extinguida y relaciones demográficas |
| Evento biológico | Identifica muestreo, cruce, mutación o cambio demográfico según su categoría | Entradas, orden, resultado y causas de omisión o fallo |
| Evento explícito | Unidad de observación independiente, con directorio y manifiesto propios | Individuos, configuración resuelta, referencias y subsemilla enlazados desde la raíz; las salidas reales se enlazan cuando existen, sin artefactos vacíos; tecnología, perfil, tanda, bibliotecas y mapeador solo cuando procedan |
| Experimento de secuenciación | Condición de un evento que genera FASTQ, BAM o CRAM; es el único experimento del evento y no posee otra identidad | Tecnología, perfil, una tanda, una biblioteca por individuo y, si genera alineamientos, mapeador |
| Tanda | Única dentro de un evento de secuenciación | Evento, tecnología, perfil y todas las bibliotecas procesadas conjuntamente |
| Biblioteca | Única para un individuo dentro del evento | Evento, tanda, individuo, tecnología, perfil, objetivo WGS e índices técnicos |
| Grupo de lectura | Único en toda la ejecución y único para el individuo dentro del evento | Evento, individuo, biblioteca, tanda y archivos asociados |
| Instancia de lectura | Observación concreta con identidad técnica única dentro del evento; una pareja constituye una unidad y comparte el nombre exigido por el formato | Biblioteca actual, miembros cuando sea pareja, origen estable y mecanismo de copia o transferencia |
| Conjunto | Local a un evento y una referencia de mapeo | Mismas muestras que los demás conjuntos del evento y resultados de esa referencia |
| Caso de verdad | Local a un conjunto | Diferencia biológica simple o compuesta y uno o varios registros VCF |
| Registro VCF | Local a un conjunto y enlazado a uno o varios casos de verdad originales | Referencia, representación estándar, lista canónica `GPS_CASE` y, si pierde detalle, el enlace complementario que se defina para ese caso real |
| Entrada complementaria | Única para un caso que pierde detalle en VCF | Identificador local estable, registros enlazados y detalle completo no representable |

Un identificador aportado usa ASCII, distingue mayúsculas de minúsculas y
cumple conceptualmente `[A-Za-z0-9][A-Za-z0-9._-]{0,63}`. No contiene espacios,
rutas ni otros caracteres. Un identificador automático usa BLAKE3 con el
contexto ASCII exacto `genome_population_simulator/v1/id` y el material binario
canónico descrito en 7.1. Toma los primeros 128 bits del resultado, los codifica
en Base32 RFC 4648 mayúscula sin relleno `=` y antepone el prefijo ASCII de
clase definido por el catálogo. El ID de ejecución incluye como causas mínimas
la raíz y el digest de configuración; los demás incluyen su tipo y los IDs
causales estables que fije el catálogo. Nunca participan reloj, hilo, ruta,
orden físico ni posición accidental. Toda colisión se detecta y hace fallar la
operación.

El catálogo de derivación de F1 fija ámbitos, nombres y orden de campos y
prefijos de clase. Junto con sus vectores de subsemillas e IDs es un artefacto
normativo: modificarlo requiere una nueva versión del esquema de derivación.

Los espacios de identidad siguen el contrato de cada entidad. Las variantes son
locales a su conjunto y VCF; `RG` es único en toda la ejecución; una lectura o
plantilla es única dentro del evento; y los mates comparten el identificador de
plantilla y se distinguen por su rol. El individuo conserva su mismo ID biológico
en `SM`, VCF y manifiesto. Los identificadores automáticos deben permanecer
compactos y válidos en SAM, VCF y FASTQ.

### 5.2. Referencias, coordenadas y segmentos

Una referencia posee una secuencia ordenada de cromosomas o *scaffolds* y una
identidad verificable del FASTA exacto. Una referencia derivada enlaza una única
referencia de origen y la secuencia ordenada de cambios estructurales que la
creó. El catálogo de referencias es un grafo dirigido sin ciclos durante una
ejecución.

Una referencia reducida es una derivación de preparación distinta de los
cambios estructurales biológicos. Enlaza una única copia autoritativa de FASTA
suministrado y una lista ordenada de segmentos en orientación directa. Cada
segmento lleva coordenadas de origen y derivadas basadas en 1 y cerradas,
identidad y clase de secuencia, posición en la concatenación y procedencia. La
lista y el FASTA canónico derivado son dos proyecciones del mismo resultado; el
mapa no se reconstruye a posteriori desde nombres físicos.

El espacio elegible contiene cromosomas y *scaffolds*. Cada secuencia derivada
conserva exactamente el `sequence_id` y la clase de su fuente; la identidad de
referencia separa los dos espacios. `reduction_percent` expresa el porcentaje
eliminado y el objetivo retenido es `floor(S × (100 − p) / 100)` calculado desde
el literal decimal TOML exacto. La frontera conserva ese literal y obtiene el
entero mediante aritmética racional comprobada, sin pasar por `f64`. Tanto el
porcentaje como `target_bases` deben producir un objetivo positivo y
estrictamente menor que `S`.

El servicio común de coordenadas posee tipos distintos para posición de base,
punto interbase e intervalo. Internamente usa magnitudes `u64`, intervalos base
0 semiabiertos no vacíos `[start,end)` y valores ligados a `ReferenceId` y
`SequenceId`. Un punto interbase admite `0..=sequence_length`. Los adaptadores
convierten explícitamente `[s,e]` base 1 cerrado en `[s−1,e)` y viceversa, con
aritmética comprobada. Ningún valor inválido se recorta, y un DTO externo nunca
se usa como valor de dominio. F2 amplía este mismo servicio para anclas,
variantes y segmentos.

`reduced_reference_selection_v1` calcula `q=T/R` y `r=T%R`. El preflight
comprueba que existe una disposición completa. Para cada región completa
enumera de forma compacta rangos de inicios no solapantes cuya elección conserva
al menos una terminación posible; todos los inicios seguros del conjunto reciben
el mismo peso. Un índice entero uniforme sin sesgo elige el inicio mediante el
RNG local. Tras las `q` regiones se aplica la misma regla al único residual si
`r>0`. No existen reintentos de disposición ni resultados parciales. El orden
de sorteo se conserva como traza, pero el plan materializado se ordena por orden
FASTA y coordenada.

En un FASTA, el identificador de secuencia es el primer elemento no vacío que
sigue a `>` y termina en el primer espacio o carácter de separación. Debe ser
único dentro de la referencia. El resto de la cabecera se conserva como
descripción completa y no participa en la unicidad.

El FASTA suministrado se trata como bytes de texto sin comprimir. Su copia se
realiza byte a byte sin reescritura, mientras el flujo valida el contenido y
calcula SHA-256; el original queda intacto y la copia solo se publica si termina
válida. Se conservan la cabecera cruda y la copia exacta. Si la descripción no
es UTF-8, sus bytes se representan de forma reversible y la salida estructurada
y el manifiesto declaran esa codificación; esa condición por sí sola no invalida
la referencia. Las bases se normalizan a mayúsculas únicamente en la
representación interna.

Los FASTA simulados o derivados se escriben de forma canónica: texto sin
comprimir, A/C/G/T/N mayúsculas, saltos LF, ancho fijo de línea y orden estable.
La escritura canónica FASTA v1 usa exactamente 80 bases por línea. Ese ancho es
parte de la versión porque afecta a los bytes. Una ejecución F1 normal produce
cero índices. FAI solo se crea bajo demanda en el espacio de trabajo de la
consulta indexada de aceptación que lo consume; ese índice no se publica como
salida normal de la ejecución.

Toda secuencia suministrada o generada usa exclusivamente `A`, `C`, `G`, `T` y
`N`, sin distinguir mayúsculas de minúsculas durante la validación. `N` expresa
una base desconocida o no determinada. Su tratamiento operativo debe decidirse
funcionalmente para cada clase de operación; el sistema no puede resolverla al
azar, sustituirla ni atribuirle otro significado de forma implícita.

Toda coordenada contiene, como mínimo, referencia, secuencia, posición o
intervalo y convenio. La configuración utiliza base 1 e intervalos cerrados;
los formatos externos conservan su convenio oficial. DTC-056 fija la
representación interna base 0 semiabierta y tipos distintos para base, punto
interbase e intervalo. Ningún componente puede convertir coordenadas sin usar
el servicio común de coordenadas y registrar la conversión.

Según DTC-061, un segmento ordinario no conserva una copia de la secuencia:
referencia un intervalo de la referencia o de una fuente fundadora exacta y
conserva cómo materializarlo y su procedencia. Una secuencia insertada se
almacena como dato propio únicamente cuando no existe en ninguna fuente. Una
deleción conserva el intervalo ausente. Segmentos y modificaciones tipadas usan
los tipos de DTC-056 y mantienen referencia/fuente, copia y fase; las
correspondencias entre estados solo existen donde la transformación las hace
significativas y nunca se inventan para regiones sin homología.

### 5.3. Estado genotípico y verdad

El estado genotípico de un individuo agrupa, para cada secuencia, exactamente el
número de copias que indica su ploidía efectiva. Cada copia contiene la cadena
ordenada y persistente de segmentos de DTC-061. SNP e *indels* pequeños son
modificaciones tipadas sobre esos segmentos, con procedencia y coordenadas
explícitas. La estructura es lógicamente inmutable y admite compartir segmentos
entre estados; DTC-061 no obliga a un crate o contenedor físico concreto. La
cadena basta para materializar bajo demanda la secuencia, conservar fase y dosis
y alcanzar la referencia, una fuente fundadora o un evento de inserción. Los
bytes materializados son una vista y nunca otra fuente de verdad.

La futura recombinación compone cadenas de segmentos y las variantes
estructurales añaden transformaciones o clases de segmento trazables. Ninguna de
las dos capacidades puede sustituir el modelo híbrido, aplanarlo de forma
permanente en secuencias completas ni perder modificaciones pequeñas ya
asociadas. DTPD-102 queda completamente resuelta por DTC-061: este mismo modelo
lógico se proyecta al fichero externo obligatorio `truth.jsonl` sin copiar el
contenedor físico Rust.

El motor genómico aplica primero los eventos explícitos en su orden declarado.
Después recibe un recuento sorteado por clase y procesa, en orden fijo,
duplicaciones, deleciones, inversiones, translocaciones, repeticiones en tándem,
*indels* y SNP. Para cada evento consulta el espacio válido del estado actual y
lo actualiza tras aceptarlo. La ausencia de destinos suficientes es un error del
motor, no una advertencia ni una autorización para producir menos eventos.

El recuento sorteado de variantes aleatorias procede exclusivamente de Poisson.
Los eventos explícitos y el modo controlado aplicable entregan resultados
concretos y no implementan otra distribución. Las mutaciones de repeticiones en
tándem usan el modelo escalonado simple: cada evento añade o elimina una unidad
completa según el sesgo efectivo entre ambas direcciones.

El servicio de recuento implementa únicamente `poisson-count-v1` de DTC-067.
Recibe `lambda` ya derivada y el RNG ChaCha12 del ámbito; cero retorna sin
avanzarlo, la rama `0<lambda<10` usa Knuth y `10<=lambda<=2^52` usa PTRS de
Hörmann. Cada uniforme abierto consume un `next_u64` mediante la conversión
normativa de 52 bits. Toda matemática trascendente pasa por `libm` 0.2.16 y el
resultado se valida antes de convertirlo a `u64`. El servicio no conoce ni
recorta el espacio genómico: el consumidor aplica después el fallo por
oportunidades insuficientes y registra versión, media, ámbito/subsemilla y
recuento.

La **verdad biológica** comprende referencias, haplotipos, genotipos,
parentescos, poblaciones, transmisiones, mutaciones y observaciones. Es la única
fuente autorizada para generar el fichero obligatorio de verdad y los VCF.

El fichero obligatorio es `truth.jsonl` v1, UTF-8 sin comprimir. Se escribe y
lee por flujo. La primera línea es siempre la cabecera versionada con clase,
ejecución, referencia, convenio de coordenadas, recuentos y digests; las líneas
siguientes son registros tipados de copia haplotípica o genotipo. Cada copia
conserva su cadena ordenada de segmentos/modificaciones, fase y procedencia; cada
genotipo identifica caso/locus, alelos, dosis y copias portadoras. Los esquemas
por tipo validan cada línea. El orden canónico comienza por referencia,
individuo, secuencia y copia y continúa por las claves internas del registro.
El manifiesto enlaza longitud y SHA-256 del fichero completo.

La validación de conjunto no se reduce a validar líneas: un reconstructor
independiente comprueba por flujo cabecera, orden, unicidad, recuentos/digests y
relaciones, materializa cada secuencia y vuelve a obtener los genotipos desde
los haplotipos. Una ausencia, duplicado, desorden, relación rota o diferencia de
alelos, dosis o copia invalida el artefacto. VCF es exclusivamente una vista
derivada de `truth.jsonl`; ningún VCF puede reconstruir, corregir o sustituir la
verdad. La compresión queda fuera de F2 y una ampliación no podrá cambiar el
modelo lógico.

La vista VCF aplica desde esos haplotipos finales el algoritmo de normalización
versionado por DTC-064: reconstruye primero cada alelo final, amplía los cambios
solapados al intervalo de referencia común necesario, recorta contexto común
sin perder el anclaje, alinea a la izquierda contra la referencia y ordena por
catálogo de secuencias y coordenada. En F2 agrupa todos los alelos compatibles
en un único registro multialélico por `CHROM+POS`, sin descomponerlo. El
materializador no invoca `bcftools`; las pruebas lo usan como oráculo externo.
La excepción para varios BND relacionados, incluso en una misma posición, se
reserva a SV/translocaciones futuras y no relaja la cardinalidad de SNP e
*indels* pequeños de F2.

Cada registro VCF conserva siempre `GPS_CASE`, declarado como `Number=.` y
`Type=String`. Contiene los IDs locales estables, únicos y ordenados
canónicamente de todos los casos originales en `truth.jsonl` que aportan uno o
varios ALT al registro; no se crea un caso compuesto. `truth.jsonl`, como
autoridad canónica, conserva para cada alelo final la relación explícita
ALT↔lista ordenada de casos fuente. Un caso que aparece en varios registros
repite su ID en cada uno. Solo una pérdida real añade `GPS_LOSS`, código
de catálogo cerrado y no texto libre, y `GPS_COMPLEMENT`, ID de la única entrada
correspondiente en `vcf-complement.jsonl`. Este JSON Lines versionado es único
por conjunto y no existe si ningún caso pierde detalle; contiene una entrada por
caso perdido con fase, haplotipo, genotipo, coordenadas, cambio completo y enlace
a la verdad. Varios registros del mismo caso, incluidos BND futuros, comparten
entrada. El manifiesto conserva ruta, tamaño y SHA-256, y el validador rechaza
enlaces ausentes, múltiples, cruzados o entradas sobrantes. Todos los casos de
F2 son representables y, por tanto, omiten `GPS_LOSS`, `GPS_COMPLEMENT` y el
complemento. La cardinalidad positiva multicaso del complemento se definirá en
F5.5 sobre un caso real con pérdida, no en F2.

La **verdad técnica de lecturas** conserva el origen real de fragmentos,
secuencias original y observada, errores, duplicados, contaminaciones,
quimeras, descartes agregados y categoría sintética cuando corresponda. Esta
verdad se deriva de la biológica y del evento técnico, pero no la sustituye.
También distingue cobertura objetivo, WGS básica anterior a artefactos, final
del FASTQ y alineada observada por referencia de mapeo.

El alcance mínimo de los estados genotípicos completos depende de la vía:

- creación directa: todos los individuos creados;
- evolución: poblaciones inicial y final, y todos los individuos de las
  generaciones o muestreos solicitados; las demás generaciones conservan el
  historial resumido obligatorio;
- familias y poblaciones segregantes: todos los individuos de todos los eventos
  de cruce y poblaciones intermedias.

Un historial resumido nunca se presenta como estado genotípico completo. Un
evento explícito puede solicitar una materialización adicional de verdad o VCF
sin activar el recorrido de secuenciación.

### 5.4. Poblaciones y generaciones

Una generación activa contiene individuos que pueden actuar como progenitores.
La descendencia crea otra generación y sustituye a la anterior. Los individuos
históricos mantienen identidad y linaje, pero no reaparecen como miembros de una
generación posterior.

El estado de una subpoblación distingue N objetivo, intentos reproductivos,
éxitos, fallos por causa, N realizado y Ne realizado o «no evaluable» con causa.
Distingue también subpoblación activa, extinción
planificada y extinción por resultado. Una recolonización crea una identidad
nueva.

Las probabilidades y frecuencias configurables cruzan la frontera de admisión
como `f64` finitos dentro de `[0,1]`. Sus colecciones normalizables conservan
identidad y categoría suficientes para establecer el orden canónico; ni su
posición física ni el reparto de trabajo adquieren semántica. La configuración
efectiva consume únicamente los valores ya validados y normalizados junto con
su registro de auditoría.

### 5.5. Eventos y resultados

Un evento de muestreo crea una observación no destructiva del estado biológico.
Todo evento es explícito y puede producir solo verdad, solo VCF o evidencia
técnica. Cuando genera FASTQ, BAM o CRAM constituye además un experimento de
secuenciación, sin otra identidad, con exactamente una tanda y una biblioteca
por individuo. Todas las bibliotecas se procesan juntas en la tanda. Cada evento
produce cero o más conjuntos, uno por referencia
aplicable. Los FASTQ pertenecen al evento y no a un conjunto concreto, porque
dentro del evento las mismas lecturas pueden mapearse contra varias referencias.

Ejecución, evento, conjunto y artefacto tienen como mínimo los estados
`incompleto`, `completo` y `fallido`. La omisión autorizada de un evento
dependiente se registra como resultado y diagnóstico, sin presentarlo como
completo. El estado de una unidad superior se calcula desde sus obligaciones y
unidades subordinadas; no se deduce de que exista algún fichero.

## 6. Componentes y responsabilidades

### 6.1. Catálogo de componentes

| Componente | Responsabilidad única | Entradas y salidas | Datos que posee o consulta |
|---|---|---|---|
| Aplicación Marimo | Editar y revisar configuración, invocar `validate`/`run` y presentar estado sin adquirir semántica | Esquema vigente y acciones humanas → TOML general, argumentos de proceso, diagnósticos, resúmenes y enlaces | Posee el estado del formulario y la última validación; no posee `EffectiveConfig`, reglas biológicas, verdad ni artefactos grandes |
| Admisión de configuración | Leer TOML, conservar su texto exacto, validar y normalizar colecciones numéricas y resolver `DeclaredConfig` en `EffectiveConfig` mediante la matriz de ejes y dominancia | Ficheros TOML y entradas declaradas → configuración efectiva con procedencia, auditoría numérica o diagnósticos localizables | Posee `DeclaredConfig`, ubicaciones, `EffectiveConfig` y registros de suma/normalización; consulta catálogos de capacidades |
| Identidad de ejecución | Validar IDs aportados, derivar IDs automáticos tipados y comprobar sus ámbitos | Configuración validada, raíz, digest y causas estables → identidades, etiquetas y contexto o colisión | Posee registro, catálogo normativo de campos y prefijos y vectores BLAKE3; no posee datos biológicos |
| Aleatoriedad reproducible | Obtener o validar la raíz de 256 bits, derivar subsemillas y proporcionar un ChaCha12 local por ámbito nombrado | Raíz, origen, versión y componentes estables → generador auditable | Posee el esquema binario, catálogo normativo de ámbitos y vectores BLAKE3; no posee un RNG global ni interpreta decisiones biológicas |
| Catálogo de referencias | Validar, copiar, simular y derivar referencias mediante un contrato independiente del adaptador FASTA | FASTA de texto sin comprimir o configuración → referencias disponibles, copia exacta e índices necesarios | Posee identidades, cabeceras crudas, secuencias internas normalizadas, linaje y recursos exactos; usa `noodles-fasta` solo tras el contrato |
| Reducción de referencias | Seleccionar segmentos elegibles con RNG local, proyectar el mapa y materializar una derivación exacta | Copia FASTA validada, objetivo, tamaño de región, subsemilla y versión → segmentos ordenados, mapa y FASTA canónico o rechazo previo | Posee el plan inmutable de segmentos y sus métricas; consulta catálogo, coordenadas, aleatoriedad y publicación, y no modifica el FASTA fuente |
| Coordenadas y transformaciones | Validar intervalos y convertir convenios sin ambigüedad | Coordenadas etiquetadas → coordenadas validadas o correspondencias | Posee reglas de conversión; consulta referencias y segmentos |
| Motor de verdad genómica | Crear y transformar haplotipos y eventos | Referencia, cambios y flujos aleatorios → estado faseado | Posee el estado genotípico, eventos, segmentos y correspondencias |
| Creación directa | Traducir la vía directa a operaciones del motor genómico | Configuración efectiva → individuos | Posee solo el estado de progreso; el motor posee la verdad resultante |
| Fuentes de variación poblacional | Preparar, admitir, generar o reutilizar espectros y matrices SNV canónicas sin crear todavía individuos; posteriormente admitir paneles de haplotipos | Tabla mínima o par TSV+JSON/configuración/semilla → fuente poblacional inmutable, par publicado, auditoría o error | Posee recursos normalizados, IDs, versiones, perfiles, descartes y estadísticas de generación; consulta referencia, coordenadas, SHA-256, aleatoriedad y publicación |
| Inicialización poblacional | Crear fundadores, frecuencias, genotipos, LD, HWE y F diploide inicial | Fuente poblacional y configuración → población inicial o error previo | Posee estadísticas solicitadas/realizadas y el estado «no aplicable» de F fuera de diploides; entrega verdad al repositorio biológico |
| Reproducción y recombinación | Validar ploidía y estructura canónica, formar bivalentes, gametos y descendencia trazable | Progenitores, mapas/tasas y flujo → descendientes o incompatibilidad reproductiva | Posee eventos de transmisión y comparaciones estructurales; consulta verdad parental |
| Dinámica evolutiva | Ejecutar las ocho etapas generacionales en su orden fijo | Población y fases → población final e historial de etapas aplicadas u omitidas | Posee estado generacional y resúmenes; delega cambios genómicos |
| Selección de viabilidad | Evaluar el SNP bialélico protegido y decidir supervivencia | Individuos elegibles, dosis y función → supervivientes o error | Posee evaluaciones y decisiones; consulta identidad única, dosis y valores genéticos |
| Grafo de familias | Validar y ejecutar el grafo declarativo de cruces | Grafo y parentales → poblaciones intermedias y finales | Posee dependencias y resultados del grafo; consulta selección y reproducción |
| Selección declarativa | Evaluar el catálogo común sobre verdad biológica | Reglas y elegibles → seleccionados y auditoría | No posee genotipos; posee evaluación y resultado de selección |
| Muestreo y resúmenes | Crear observaciones no destructivas y ámbitos estadísticos | Generación y regla → muestra, resumen poblacional y muestral | Posee observaciones y resúmenes; consulta verdad e historial |
| Orquestación de eventos | Validar salidas y parámetros condicionales y aislar cada observación | Observación/configuración → plan de evento y conjuntos | Posee estado del evento y relaciones entre artefactos; solo crea un plan técnico si procede |
| Simulación WGS | Generar fragmentos e instancias de lectura desde haplotipos | Verdad, perfil y objetivo → lecturas con identidad y origen | Posee lecturas previas al demultiplexado y verdad técnica asociada |
| Artefactos y demultiplexado | Convertir conjuntamente la contaminación y aplicar errores, duplicados, adaptadores y reasignaciones | Lecturas y configuración → FASTQ finales y procedencia | Posee transformaciones técnicas e identidades de instancia; consulta fuentes compatibles |
| Alineamiento | Ejecutar el mapeador sobre el FASTQ final | FASTQ y referencia → alineamientos observados | Posee invocación, salida, diagnóstico y metadatos del proceso |
| Alineamiento directo | Crear registros sintéticos coherentes sin FASTQ | Verdad y perfil directo → alineamientos sintéticos | Posee categorías y registros sintéticos; no llama a un mapeador |
| Materialización de verdad | Aplicar el único contrato obligatorio de trazabilidad, proyectar/escribir/validar `truth.jsonl` y derivar después VCF/complemento | Estado faseado, historial y referencia de conjunto → `truth.jsonl`, vistas interoperables y resúmenes | No altera el dominio; posee esquemas por registro, escritor por flujo y validación semántica independiente; materializa el detalle exigido por vía sin niveles ni degradación |
| Ensamblado de resultados | Organizar recursos comunes, eventos y conjuntos | Artefactos validados → estructura de resultados | Posee catálogo lógico de artefactos y estados de completitud |
| Validación interna | Comparar configuración, verdad y evidencia | Objetivos y resultados → informe, aviso o error | Posee evidencias de comprobación; no corrige resultados |
| Diagnósticos públicos | Traducir errores internos tipados a un contrato público estable y renderizarlo | Error causal y contexto → texto humano a stderr o documento JSON único y código de salida | Posee catálogo versionado de códigos, esquema y renderizadores; no expone tipos de dependencias ni decide la consecuencia funcional |
| Manifiestos y exportación | Escribir los manifiestos JSON completos, enlazar artefactos independientes y crear paquetes autocontenidos | Catálogo de artefactos y ejecución → `manifest.json` general y por evento, enlaces y paquete | Enumera directamente individuos materializados y de familias; resume generaciones evolutivas no materializadas y posee enlaces e integridad SHA-256 por fichero |
| Construcción y distribución | Producir y validar las entregas oficiales portables para Ubuntu 24.04 LTS x86-64 GNU | Código, esquemas, licencias, ejemplos y versiones aprobadas → paquete `.tar.gz` e imagen OCI identificables | Posee recetas fijadas, toolchain aprobado y metadatos; no altera resultados ni decide nuevas dependencias |
| Estimación y observación de recursos | Estimar y medir tiempo, memoria y almacenamiento | Configuración, trazabilidad y salidas activadas → estimaciones y medidas | Posee mediciones; no decide modelos biológicos ni reduce resultados para ajustarlos al presupuesto |

### 6.2. Contratos comunes de los componentes

Cada componente ofrece una interfaz pequeña expresada mediante conceptos del
dominio. No recibe estructuras físicas de otro componente salvo en las
fronteras de E/S. Toda operación devuelve un resultado válido o un error con:

- categoría y código estable de diagnóstico;
- entidad, requisito y ámbito afectados;
- valor recibido y condición incumplida cuando sea seguro registrarlos;
- causa original y cadena de contexto;
- alcance de la detención: configuración, ejecución, evento, conjunto o rama;
- artefactos provisionales que no deben publicarse como completos.

Los componentes no registran dos veces la misma regla. Admisión decide
precedencias de configuración; coordenadas decide conversiones; identidad
decide unicidad; el motor genómico decide validez estructural; el orquestador
generacional decide el orden; el ensamblador decide completitud.

Admisión posee además el contrato numérico común: representación `f64`, validez
de intervalos, orden canónico, suma compensada de Neumaier, tolerancia absoluta
y normalización. Los consumidores no vuelven a sumar con otro orden o algoritmo
ni normalizan por su cuenta.

Una interfaz de dominio recibe el tipo de identificador de la categoría que
espera. Convertir desde texto solo ocurre en la admisión o en un adaptador; una
etiqueta, ruta o nombre físico nunca se convierte implícitamente en identidad.

### 6.3. Dependencias permitidas y prohibidas

- Marimo puede depender de esquemas y catálogos publicados y del adaptador de
  proceso del binario. No puede depender del motor, del resolvedor interno ni de
  escritores de resultados; su única configuración ejecutable es el conjunto
  TOML que presenta y entrega a `validate`/`run`.
- Reducción de referencias puede consultar catálogo, coordenadas, RNG local e
  integridad. No puede escribir el origen, usar un RNG compartido, reinterpretar
  DTC-050/DTC-056/DTC-057 ni crear una representación local de coordenadas.
- Creación directa, inicialización, evolución y familias pueden depender del
  motor genómico, aleatoriedad, identidad y selección declarativa.
- Todo consumidor de aleatoriedad recibe un generador local y solo avanza ese
  estado. Ningún componente comparte un RNG mutable ni deriva identidad de
  posiciones en listas, orden de iteración, hilos o planificación.
- Reproducción es compartida por evolución y familias; no se duplica por vía.
- Reproducción compara por secuencia la ploidía y el estado estructural canónico
  de todas las copias antes de formar bivalentes. La equivalencia no depende de
  que las historias compartan identificadores, pero no elimina ninguna traza.
- El motor genómico debe consultar la identidad del locus seleccionado para
  rechazar eventos explícitos solapantes y excluirlo de los espacios aleatorios;
  no puede moverlo, copiarlo ni eliminarlo.
- Reproducción y dinámica evolutiva deben rechazar antes de formar descendencia
  cualquier unidad reproductiva que mezcle referencias estructurales. Esta
  validación no restringe las referencias de mapeo de un evento explícito.
- Simulación WGS, artefactos, alineamiento y modo directo pueden consultar
  verdad e identidades, pero no llamar a componentes que cambien genotipos.
- Materialización VCF depende del estado faseado y de la referencia del
  conjunto; no depende de lecturas ni alineamientos.
- Validación depende de vistas de solo lectura de configuración, verdad y
  evidencia. No puede normalizar o reparar archivos después del hecho.
- Manifiestos reciben hechos ya registrados por los responsables. No
  reconstruyen decisiones a partir de nombres de fichero o mensajes de log.
- Ningún componente de dominio depende de crates o herramientas concretas. Los
  adaptadores implementan contratos propios y controlan esa traducción.
- Un proceso externo se invoca únicamente mediante un adaptador que controla
  versión, parámetros, comando, entradas y salidas, código de salida, stdout,
  stderr y fallos; el proceso no adquiere autoridad de dominio.

### 6.4. Comprobación aislada

Cada componente admite dobles de prueba para sus dependencias externas. Las
pruebas del dominio usan referencias y poblaciones pequeñas construidas de forma
controlada. Los adaptadores de formato se prueban con casos válidos e inválidos
del estándar elegido. Los procesos externos se prueban con respuestas de éxito,
salida inválida, código de error, interrupción y versión incompatible.
El runtime realiza siempre sus validaciones internas de artefactos, invariantes,
tamaño, integridad y publicación. Los validadores independientes se reservan a
pruebas, CI y aceptación de publicaciones; una ejecución normal no los invoca ni
exige que la persona usuaria los instale.

| Componente | Errores propios principales | Comprobación aislada mínima |
|---|---|---|
| Admisión de configuración | Sintaxis TOML; clave/campo desconocido o duplicado; versión ausente; precedencia conflictiva; `f64` no finito o fuera de intervalo; tolerancia inválida; suma no positiva/no finita o delta excesivo; referencia inexistente y solicitud fuera de v1 | Resolver casos mínimos y conservar TOML; probar tolerancias 0, `1e-9`, `1e-6` y fuera de rango, `delta == t` y valores inmediatamente dentro/fuera, NaN/infinito, cero total, reorden físico y rechazo; comprobar vectores Neumaier, recomputación, decimal round-trip y registro antes/después |
| Aplicación Marimo | Binario ausente/incompatible; intento de *shell*; estado validado obsoleto; destino existente; salida excesiva | Probar versión declarada, proceso falso/real, argumentos separados, invalidación tras editar, truncado visual de logs, arranque en `127.0.0.1` y controles básicos de teclado/etiquetas/color. La auditoría completa de navegador, plataforma y distribución pertenece a F9 |
| Identidad de ejecución | Sintaxis inválida, categoría o ámbito incorrectos, colisión o reutilización | Probar sintaxis y sensibilidad, etiquetas repetidas, vectores normativos del contexto ID, material binario y Base32 sin relleno; causas mínimas del ID de ejecución; colisiones y cruces de tipos; y rechazo de cambios de catálogo sin nueva versión |
| Aleatoriedad reproducible | Raíz mal codificada, CSPRNG fallido, ámbito desconocido, codificación ambigua, identificador inestable o estado no regenerable | Probar raíces y subsemillas de exactamente 64 hex minúsculos, generación por CSPRNG y registro de origen, vectores normativos del contexto subseed y enteros big-endian; reconstruir ámbitos y demostrar independencia al añadir o reordenar observaciones |
| Catálogo de referencias | FASTA comprimido o inválido, secuencia vacía, identificador de primer elemento vacío o repetido, símbolo distinto de A/C/G/T/N, linaje cíclico o derivación imposible | Copiar sin cambiar bytes mientras se valida y calcula SHA-256; separar identificador, descripción y cabecera cruda; conservar una descripción no UTF-8 reversiblemente; comprobar escritura canónica e índices bajo demanda; conservar `N` con coordenada y longitud; y reconstruir una referencia pequeña desde su origen y cambios |
| Coordenadas y transformaciones | Convenio ausente, intervalo fuera de límites o correspondencia inexistente | Probar extremos, inserciones, inversiones y conversiones de ida y vuelta donde sean definidas |
| Motor de verdad genómica | Evento solapado prohibido, objetivo eliminado, destinos válidos insuficientes, punto de ruptura o ploidía inválidos | Materializar haplotipos conocidos; rechazar SNP o motivo sobre `N`, conservarla dentro de SV e inserciones explícitas y agotar de forma controlada el espacio válido |
| Creación directa | Individuo incompleto o cambio incompatible | Crear individuos pequeños con cambios concretos y aleatorios controlados |
| Fuentes e inicialización poblacional | Paquete/hash/referencia inválido; bloque/sitio/alelo/frecuencia incompatible; locus `alt_af` dentro de bloque; ejes o ámbitos ambiguos; recombinante+controlado o distribución explícita; r fuera de (0,1]; F fuera de diploides; ID, asignación, segmento, recuento o enlace incoherente | Validar corpus y corrupción con oráculo ajeno al productor; reconstruir HAP y mosaicos; probar productor directo, F=0/intermedio/1, puntos conocidos/múltiples/invisibles, fuentes solo del panel, IDs por contenido, suma C y `n_h/C`; repetición, aislamiento por identidad/bloque y paridad Marimo↔TOML |
| Reproducción y recombinación | Referencias distintas, ploidía impar o desigual, estado estructural no idéntico, mapa incompleto o bivalente imposible | Probar diploides y tetraploides idénticos; SNP, *indels* pequeños y repeticiones distintos; cada clase SV heterocigota; SV fijada con historias distintas; *scaffold* intacto; SV de novo seguida de otro intento |
| Dinámica evolutiva | Orden alternativo solicitado, población activa vacía no autorizada, mezcla de referencias estructurales tras migración o evento demográfico imposible | Ejecutar las ocho etapas, omitir una opcional en su posición, rechazar otro orden, ordenar varios eventos demográficos y comprobar el historial efectivo |
| Selección de viabilidad | Locus ausente, múltiple, ambiguo o situado sobre `N`; dosis sin valor; aptitud negativa, no finita o máximo cero | Probar varias ploidías, rechazo de `N`, fallos antes de supervivencia y recombinación que conserva la identidad única |
| Grafo de familias | Ciclo real, fuente futura, cardinalidad o salida imposible | Validar grafos válidos ramificados y ciclos transitivos mínimos |
| Selección declarativa | Regla desconocida, conjunto ambiguo o elegibles insuficientes | Evaluar Y/O/NO y cantidades exactas/proporcionales sobre un catálogo de verdad fijo |
| Muestreo y resúmenes | Duplicado, estratos solapados, tamaño excesivo o ámbito estadístico confundido | Comparar resumen poblacional y muestral y demostrar observación no destructiva |
| Orquestación de eventos | Evento implícito, parámetro técnico indebido, tecnología múltiple, perfil o mapeador incompatible y referencia ajena | Construir eventos solo de verdad y eventos técnicos independientes, y comprobar requisitos condicionales y aislamiento |
| Simulación WGS | Objetivo ausente/doble, perfil incompatible, identidad repetida o fragmento imposible | Comprobar recuentos, cobertura haploide, aporte por copia y muestreo de `N` con Phred 0 antes de sesgos |
| Artefactos y demultiplexado | Conversión conjunta inviable, exceso de contaminación, colisión de identidad o par huérfano | Resolver una asignación entera conocida; reconstruir origen y diferencia verdadera/observada de cada lectura final, y comprobar agregados de descartes |
| Alineamiento | Herramienta ausente, versión incompatible, código no cero o salida inválida | Sustituir el proceso externo y comprobar comando, contexto y propagación de cada fallo |
| Alineamiento directo | Proporciones, MAPQ, banderas, CIGAR o pares incoherentes | Generar un caso por categoría y validarlo contra la referencia y verdad de origen |
| Materialización de verdad | VCF no conforme, SNP con `REF=N`, alelo no reconstruible, cardinalidad ambigua, enlace roto o pérdida de fase | Validar casos sin/con complemento, BND, dos referencias y cambios con `N`; reconstruir cada caso enlazado sin inventar una base y compararlo con la verdad canónica |
| Ensamblado de resultados | Artefacto activado ausente o mezcla de referencias | Intentar cerrar conjuntos completos e incompletos y verificar estados resultantes |
| Validación interna | Evidencia insuficiente, invariante rota o desviación fuera de política | Inyectar resultados conocidos y comprobar conforme, aviso y error sin corrección automática |
| Diagnósticos públicos | Código o categoría desconocidos, contexto perdido, causa de crate filtrada, JSONL o renderizados divergentes | Probar catálogo inmutable y un caso por categoría; humano a stderr; un único documento JSON por invocación con clase, versión y esquema; salidas 0/1/2/3/4/5/70, avisos con 0; conservación causal solo en detalle; y ausencia de `panic` ante fallos esperables |
| Manifiestos y exportación | JSON o esquema inválido, enlace roto, integridad fallida, destino existente, carrera o interrupción | Probar temporales por artefacto; rechazo pre-simulación sin destino; fallo controlado final `fallido` con diagnóstico, catálogo y unidades independientes exportables; SIGKILL con provisional `INCOMPLETE` y diagnóstico opcional; publicación completa; `EEXIST` concurrente; y rechazo sin fallback cuando no exista no-clobber, en WSL2/ext4, Ubuntu nativo y el FS OCI admitido |
| Construcción y distribución | Baseline, triple o CPU no portables; receta OCI mutable; digest incorrecto; contenido incompleto o divergente | Comprobar Ubuntu 24.04 LTS y `x86_64-unknown-linux-gnu`, rechazar `target-cpu=native`, validar SHA-256 y contenido Marimo del `.tar.gz`; validar que OCI sigue siendo CLI y no contiene Python/Marimo/puerto UI, además de sus digests, config, capas y arranque |
| Estimación y observación de recursos | Dato no medible, estimación inválida o presupuesto excedido | Comparar estimación y medida; ejecutar A en cada cambio, la serie B y la observación C según su frecuencia; y rechazar antes de simular un caso insuficiente sin degradar trazabilidad ni salidas |

## 7. Flujos de datos y ejecución

Esta sección describe la **composición técnica** de los componentes. No redefine
el comportamiento funcional: las condiciones, órdenes, cardinalidades y
errores observables pertenecen a los REQ de la especificación, y los ACE
indican cómo comprobarlos.

### 7.1. Admisión y preparación

En F2, Marimo mantiene un formulario para un único TOML general de la vía
directa. Puede abrirlo, editar los campos implementados, previsualizar el texto
resultante y guardarlo en una ruta elegida por la persona usuaria. Los escenarios
y el conjunto multidocumento se incorporarán cuando una fase posterior tenga su
primer consumidor real.

El servidor escucha en `127.0.0.1`. El lanzador muestra la URL, admite un
puerto local no privilegiado y puede omitir la apertura del navegador. Abrir,
importar o editar nunca ejecuta una simulación.

Marimo invoca directamente el binario Rust con argumentos separados:
`validate <general.toml>` y
`run <general.toml> --output <destino>`. El binario conserva toda la
validación y semántica. La UI muestra sus diagnósticos, código de salida, progreso
y enlaces o resúmenes finales. Una edición posterior invalida la última
validación; `run` solo se habilita mediante una acción explícita.

F2 no define un protocolo paralelo de capacidades ni copia contratos a un
almacén privado. La aplicación y el binario del mismo árbol comparten los
esquemas versionados del repositorio. Una versión incompatible o una opción no
implementada produce un diagnóstico y bloquea la acción.

La persona usuaria puede confirmar únicamente las excepciones específicas que
el binario F2 implemente. Estas confirmaciones se transmiten como argumentos
separados y obligan a validar de nuevo. No existe un `--force` general.

Mientras un proceso está activo, la UI muestra un estado de progreso; al
terminar presenta sus logs de forma acotada y los artefactos grandes mediante
sus rutas o resúmenes, sin cargarlos. F2 no necesita un daemon, una base de datos de
sesiones, un protocolo de cancelación de grupos de procesos ni un sistema de
expedientes de seguridad. El cierre y la publicación de resultados siguen
reutilizando los mecanismos locales de F1 frente a fallos accidentales.
### 7.2. Referencias y estado biológico

El catálogo de referencias es responsable de identidad, metadatos y acceso a
secuencias. El servicio de coordenadas y transformaciones es la única frontera
que traduce posiciones entre una referencia estructural, sus estados derivados
y las referencias de mapeo. Los motores biológicos reciben estas interfaces y
producen cambios de estado, nunca rutas de fichero ni registros de un formato
externo.

Para una referencia suministrada, el adaptador `noodles-fasta` lee por flujo
detrás del contrato del catálogo. El mismo recorrido crea la copia exacta,
calcula SHA-256 y entrega al validador de dominio cabecera y secuencia; un fallo
impide publicar la copia. El catálogo conserva los bytes crudos de cabecera y
expone bases internas en mayúsculas. Para referencias generadas, el flujo
`reference_generation` usa solo el ID estable y el escritor aplica el FASTA
canónico v1 de 80 bases por línea. La ejecución normal conserva el catálogo de
índices vacío. La consulta indexada de aceptación crea y consume un FAI bajo
demanda fuera del resultado publicado. Un FASTA comprimido se rechaza como
formato no admitido, no se descomprime implícitamente.

Para una reducción, el preflight termina antes de sortear. El planificador recibe
el tamaño elegible, objetivo entero, tamaño de región, ID y RNG local; aplica
DTC-057 y devuelve todos los segmentos o un error, nunca un prefijo parcial.
Otra operación los ordena por secuencia fuente y coordenada, asigna coordenadas
derivadas acumuladas y crea el mapa inmutable. El escritor FASTA consume ese
plan una vez, copia `N` como cualquier otra base y aplica el formato canónico.
Un comprobador independiente reconstruye el FASTA desde origen+mapa antes de
publicar conjuntamente mapa, FASTA, hashes y manifiesto.

El mapa es un documento JSON cerrado con `document_kind` igual a
`reference_reduction_map` y `schema_version` igual a `1`. Identifica referencias
origen/derivada, algoritmo y cada segmento mediante orden, clase
completa/residual, `sequence_id`, clase de secuencia, intervalos origen y
derivado base 1 cerrados y orientación directa. Sus filas siguen el orden
materializado; el manifiesto conserva además el orden original del sorteo y la
subsemilla. El mapa no usa posiciones base 0 como contrato externo.

La fuente de variación poblacional de DTC-054 materializa un recurso inmutable
formado por `population_snvs.tsv` y `population_snvs.metadata.json`. El JSON
enlaza por SHA-256 el TSV y la referencia exacta; configuración y manifiesto
enlazan el JSON. La preparación recibe una tabla mínima, calcula SHA-256 e IDs,
recibe `reference_id` y el FASTA autoritativo exclusivamente desde la entrada
seleccionada de un manifiesto general F1 validado conforme a DTC-060, comprueba
ruta, tamaño, SHA-256, coordenadas base 1 cerradas y `REF`, filtra `ALT`
globalmente ausentes, ordena y publica el par sin modificar el original. La generación independiente
y la integrada llaman al mismo servicio y producen los mismos bytes, por lo que
el par puede reutilizarse sin volver a sortear posiciones o frecuencias.

El planificador de candidatos recibe exactamente un recuento o una realización
Poisson obtenida de media/densidad y selecciona bases A/C/G/T sin reemplazo. La
semántica es DTC-022 y el algoritmo, límites, consumo RNG y vectores son los de
DTC-067; ningún consumidor implementa una variante local. Un
selector espacial uniforme o por intervalos no solapados asigna pesos por base y
`outside_weight`; otro servicio elige uno de los tres `ALT` mediante los 12
cambios dirigidos. Cada población consume su espectro no solapado de frecuencias
en modo estocástico o exacto. En el segundo se calculan primero las cuotas por
mayores restos; luego una permutación reproducible e independiente por población
de todos los `variant_id` canónicos asocia las filas siguiendo su orden estable.
Esta asignación usa una subsemilla exclusiva y nunca el orden de coordenadas o
el orden canónico de variantes. Un flujo distinto por población×`variant_id`
obtiene después el valor dentro del intervalo; `fixed` no consume RNG para el
valor y una cuota cero no consume RNG. Todo intervalo es `[min,max)`, incluso si
`max=1`; la corrección numérica DTC-072 impide que un intervalo continuo produzca
0 o 1. Los nombres, componentes, orden y versión
literales pertenecen al catálogo RNG v4 con extensión v5 fijado por DTC-071, que cierra RAF2-028
sin modificar DTC-070. Posición, `ALT`, frecuencia, genotipo y fase usan
ámbitos RNG separados y estables. El filtrado posterior de
un `ALT` cero en todas las poblaciones no repone el candidato; un valor fijado se
conserva. V1 solo genera loci bialélicos, aunque el lector canónico admite
multialelia suministrada con suma ALT validada mediante la tolerancia efectiva
de DTC-049, no mediante su valor predeterminado fijo.

Las inserciones y deleciones aleatorias no reutilizan ámbitos ni selectores de
la matriz. Cada clase separa por referencia, individuo y región el sorteo del
recuento, la selección espacial y, solo para inserciones, el contenido de bases;
una media cero no consume RNG. DTC-059 aporta longitudes exactas sin sorteo. El
motor enumera oportunidades sobre el estado vigente, sortea una vez mediante
`poisson-count-v1` y entrega el recuento a `dynamic-indel-selection-v1`. Este
consumidor no preflighta una colección inicial: por ordinal reenumera y pondera
canónicamente el estado actual y consume el mismo flujo `indel_position`. El
orquestador completa una inserción con el flujo aislado `insertion_content`,
aplica el evento y solo entonces pide el ordinal siguiente. El agotamiento antes del recuento falla y descarta el estado
transaccional completo, aunque registra el número aceptado antes del fallo. La
selección estática `weighted-without-replacement-v1` queda reservada a loci SNV.

Los dos perfiles humanos aproximados se incorporan mediante la receta compacta
DTC-068 y nunca requieren el VCF poblacional. El generador de recursos toma los
literales SFS 64:12:8 y el mínimo versionado `1/5008`; para SNV12 consume el
fichero eLife fijado de 3.644 bytes solo después de verificar su SHA-256, normaliza
por fila continental, promedia las cinco filas, agrega los 96 contextos a seis
clases de hebra y refleja cada clase a dos sustituciones dirigidas. Un oráculo
independiente repite los cálculos y compara orden, pesos normalizados y SHA-256
de las dos salidas pequeñas. El catálogo las etiqueta `approximate`, conserva
fuentes, licencia CC BY y transformación y nunca las selecciona por omisión.
No existen builtins exactos por superpoblación; perfiles suministrados por
población siguen usando el contrato general. `bcftools` permanece como validador
VCF de DTC-066, no como dependencia de estos perfiles.

El inicializador consume esa misma matriz para todos los individuos, asigna a
cada uno una única población y muestrea genotipos con Hardy-Weinberg o `F`
diploide; nunca vuelve a elegir loci por individuo. Después permuta los alelos
entre copias con un ámbito de fase independiente. Un genotipo explícito de
probabilidad cero solo atraviesa la admisión con la autorización específica
`--allow-incompatible-genotypes`; una matriz activa vacía exige
`--allow-empty-snv-matrix`. Ambas decisiones se incorporan a configuración
efectiva/manifiesto y no relajan otras invariantes. La matriz conserva solo
frecuencias marginales, por lo que los loci son independientes y cualquier LD
requiere el panel de haplotipos posterior.

F2.3.4 admite `direct_events` únicamente como una vista tipada por individuo y
fase: su adaptador resuelve ID, contig, alelos y coordenadas, y el orquestador
prepara una vez la matriz, muestrea y fasea antes de construir estados lógicos
temporales. La fase conserva la identidad de referencia y SHA-256 canónico de
la matriz; cualquier diferencia de locus o frecuencia rechaza el enlace. Para
un SNP explícito que coincide con la matriz, la compatibilidad se decide por
probabilidad positiva bajo HWE o `F` diploide, no contra el alelo sorteado; la
autorización específica solo conserva un caso de probabilidad cero. El
preflight transaccional finaliza antes de crear el workspace, conserva el orden
de cambios por individuo y no materializa ni publica las copias.
`validate` y `run` llaman al mismo preflight de solo lectura; únicamente el
segundo abre el workspace de publicación.

El panel de bloques haplotípicos cruza la misma frontera como otro recurso
inmutable. El paquete contiene exactamente `panel_blocks.tsv`, `panel_sites.tsv`,
`panel_haplotype_alleles.tsv`, `panel_frequencies.tsv` y `panel.metadata.json`.
El JSON es la raíz del recurso: declara clase, versión, `panel_id`, referencia y
SHA-256, coordenadas base 1 cerradas y, para cada TSV, clase, ruta, tamaño y
SHA-256. Registra procedencia y herramienta cuando existan, pero no duplica las
tablas. El adaptador verifica todos los enlaces antes de interpretar la fuente.

El adaptador tabular comparte un contrato byte a byte para los cuatro TSV:
UTF-8 sin BOM, LF incluido al final, campos separados por tabulador, una cabecera
exacta, sin comillas, vacíos, líneas vacías ni columnas extra. Las cabeceras son
`block_id reference_id sequence_id start end`,
`block_id site_id position ref alt`,
`block_id haplotype_id site_id allele_index` y
`population_id block_id haplotype_id frequency`, donde los espacios de esta
enumeración representan tabuladores en el fichero. Enteros positivos y no
negativos usan decimal canónico sin signo ni ceros iniciales y `frequency` la
representación decimal mínima *round-trip* de DTC-049.

Cada fila de `panel_blocks.tsv` declara `block_id`, `reference_id`,
`sequence_id`, `start` y `end`; sus intervalos no se solapan. Cada fila de
`panel_sites.tsv` declara `block_id`, `site_id`, posición, `REF` y los `ALT`
ordenados, y todos los sitios quedan dentro de su bloque.
El adaptador solo admite una referencia por paquete: contrasta todos los
`block.reference_id` con `metadata.reference.reference_id`, recalcula el
SHA-256 del FASTA autoritativo y lo contrasta con
`metadata.reference.sha256`, y comprueba que cada `sequence_id` existe en ese
FASTA. Cualquier mezcla, ID o digest divergente o secuencia inexistente falla en
la frontera antes de entregar el panel al inicializador.
Cada haplotipo pertenece a un único bloque y su identidad normativa combina
`block_id` y `haplotype_id`; no existe una columna o identidad de donante que
persista entre bloques. V1 admite aquí únicamente SNP suministrados cuyos `REF` y
`ALT` son una base A/C/G/T. La columna `alt` serializa uno o varios ALT distintos
mediante comas sin espacios; el orden textual define 1..N y `ref` define 0. La frontera rechaza las demás clases de variante,
`N`, alelos simbólicos, *breakends* y sitios repetidos o solapados.

`panel_haplotype_alleles.tsv` es una relación densa: cada pareja bloque-haplotipo
contiene exactamente una asignación por sitio, con 0 para `REF` y 1..N para los
`ALT` declarados. El adaptador rechaza huecos, duplicados y valores ajenos al
catálogo y nunca completa un hueco como `REF`. Tras ordenar por los sitios del
bloque, compara el vector completo de `allele_index` de todos sus haplotipos. Si
varios IDs tienen el mismo vector, rechaza el panel e informa el bloque y todos
los IDs equivalentes; no los fusiona ni suma sus frecuencias. La unicidad solo se
exige dentro del bloque y también se aplica a haplotipos con frecuencia 0.
`panel_frequencies.tsv` también es densa: contiene una fila por población,
bloque y haplotipo, incluidos 0 y 1.
No existe frecuencia global ni herencia entre poblaciones; ausencia o duplicado
son errores y cada suma población×bloque usa DTC-049.

La escritura ordena bloques por referencia, orden FASTA de secuencia, límites e
ID; sitios por bloque, posición e ID; alelos por bloque, HAP ASCII y sitio; y
frecuencias por población ASCII, bloque y HAP ASCII. El adaptador exige las
claves únicas de bloque, bloque+sitio, bloque+HAP+sitio y
población+bloque+HAP, respectivamente, antes de calcular hashes.

Cada bloque debe contener como mínimo dos sitios y dos vectores HAP distintos.
El catálogo no se reduce por las frecuencias de las poblaciones: una o todas
pueden concentrar frecuencia 1 en un único HAP y mantener a 0 los restantes. El
validador conserva un sitio fijado en el mismo ALT en todos los HAP, pero rechaza
como redundante un sitio cuyo `allele_index` sea 0 en todos ellos.

F6.4 admite y valida el recurso y trata cada bloque como un locus multialélico
cuyas categorías son sus haplotipos locales. Para cada individuo de ploidía P,
el inicializador asigna P categorías a copias haplotípicas biológicas
identificadas. Reutiliza HWE/multinomial para P=1..10 y exactamente el modelo F
multialélico solo en diploides; F=0 coincide con HWE, F=1 asigna el mismo HAP a
ambas copias y F se rechaza para P=1 o P>2. En construcción intacta, una
distribución explícita de combinaciones puede prevalecer localmente conforme a
REQ-EVO-49.

La configuración mantiene dos ejes separados. `haplotype_construction` elige
`intact` o `recombinant`, mientras `frequency_sampling` elige `probabilistic` o
`controlled`. El primero admite un valor general y sustituciones por población y
bloque mediante el resolvedor DTC-033; la regla conjunta población+bloque
prevalece. Una unidad efectiva no mezcla construcciones, aunque otras unidades
usen otra. `intact` rechaza parámetros de recombinación y `recombinant` exige
`recombination_probability_per_interbase` y `frequency_sampling=probabilistic`.
La combinación recombinante+controlado falla antes de individuos con ambos
campos y el ámbito efectivo. Marimo proyecta esta misma unión discriminada y
oculta el campo de probabilidad cuando no corresponde.

El selector genotípico también forma una unión validada: `intact` admite HWE, F
diploide o distribución explícita; `recombinant` solo HWE o F diploide. La
combinación recombinante+distribución explícita falla antes de individuos e
identifica población, bloque y campos. La arquitectura no reserva una ruta
parcial para condicionar distribuciones explícitas sobre mosaicos.

Cada población×bloque intacto usa el modo probabilístico predeterminado o el
controlado explícito; el recombinante usa únicamente el probabilístico. El
probabilístico sortea copias con reemplazo según las frecuencias.
El controlado construye un multiconjunto desde recuentos enteros por HAP que
suman C y lo asigna sin alterar los totales. Ambos se condicionan al modelo
genotípico. Una incompatibilidad se rechaza antes de crear individuos; nunca se
ajustan recuentos ni se relaja el modelo. Los ámbitos aleatorios separan bloques
y ligan reproduciblemente cada asignación a la identidad de copia, sin expresar
un donante común entre bloques. Se registran esperado, solicitado y realizado.

Los SNP exteriores a los intervalos de bloque se obtienen de
la matriz independiente `alt_af`, de modo que un individuo puede combinar
haplotipos distintos entre bloques sin que las bases intermedias queden forzadas
a la referencia. Antes de inicializar ningún individuo, Fuentes de variación
poblacional comprueba cada posición de la matriz contra los intervalos del
panel. Un locus situado dentro de un bloque produce un error con referencia,
secuencia, posición y fila o identidad de matriz, además de `block_id` e
intervalo. La frontera no expone prioridad, fusión ni sobrescritura y no entrega
fuentes parcialmente admitidas al inicializador.

La inicialización no mantiene estado de transición, donante común ni objetivo de
LD entre bloques. El LD interbloque que aparezca por tamaño muestral se informa
como realizado cuando corresponda, pero no condiciona ni repite asignaciones.
DTF-044 reserva para una versión posterior las dependencias haplotípicas entre
bloques. Esta frontera no altera el productor recombinante local descrito a
continuación: sus fuentes, puntos y segmentos siguen perteneciendo a un único
bloque.

El modo recombinante también pertenece a la inicialización fundadora, no al flujo
de meiosis histórica. No construye un depósito prototipo ni vuelve a muestrear.
Un productor de copia final selecciona un HAP primario con las frecuencias del
panel. Con `recombination_probability_per_interbase=r` finita en (0,1], prueba
independientemente las L−1 uniones interbase: K~Binomial(L−1,r). K=0 devuelve la
fuente primaria intacta. K≥1 selecciona una secundaria con las mismas
frecuencias y delega al constructor común los puntos únicos, ordenados e internos
y la alternancia de segmentos. DTC-055 cierra para F6.4 tanto el sorteo uniforme
de puntos como este constructor: recibe dos fuentes y puntos ordenados y devuelve
segmentos reconstruibles. DTPD-111 no es una dependencia de F6.4; en F7.1
conectará mapas y semántica meiótica al mismo constructor. El punto x representa
la unión x|x+1. El constructor
empieza con la primaria en `start`, alterna en cada punto y produce segmentos
cerrados *forward* sin inversión; 100–300 con puntos 150/240 produce P 100–150,
S 151–240 y P 241–300. No ajusta puntos a sitios SNP. Las fuentes pueden coincidir, se admiten varios
cruces y no se corrigen ni descartan mosaicos. Los mapas posteriores cambian las
probabilidades por posición, no el motor.

Las dos selecciones consultan exclusivamente el panel original y las frecuencias
de su población. El repositorio de verdad recibe las salidas, pero nunca vuelve
a ofrecerlas al productor durante esta inicialización. No se actualizan pesos ni
catálogos, de modo que el orden de procesamiento no cambia la distribución. La
recombinación acumulada y el uso de descendientes como fuentes pertenecen a
F7/F8.

La traza conserva todos los puntos generados, aunque fuentes iguales, regiones
sin SNP o un número par de alternancias produzcan el mismo vector alélico. La
identidad por contenido y la procedencia histórica son capas distintas.

El inicializador invoca ese productor directamente por individuo y bloque. En
diploides sortea primero Bernoulli(F): en autozigosis ejecuta un productor y
asigna el mismo vector a las dos copias; en caso contrario ejecuta dos productores
independientes. Así implementa directamente las fórmulas multialélicas aprobadas;
F=0 es HWE y F=1 produce homocigosis completa sin representar reproducción o
pedigrí. Para P distinta de 2 ejecuta P productores independientes y rechaza F.
Las C=N×P copias finales se agrupan por vector y la frecuencia realizada es
`n_h/C`. Las frecuencias del panel siguen siendo las probabilidades objetivo de
las fuentes; no se crea otra distribución intermedia.

La identidad final depende del contenido. Un vector que coincide con un HAP del
panel reutiliza su identidad local. Para un vector nuevo, el servicio común de
IDs aplica DTC-039 y deriva un ID tipado y determinista desde `panel_id`,
`block_id` y el vector canónico ordenado de `allele_index`. Individuo, orden,
subsemilla y evento no participan en la identidad. Por tanto, vectores iguales
comparten ID y recuento y vectores distintos tienen IDs distintos; las colisiones
se tratan con el contrato común, no con otro algoritmo. La procedencia permanece
separada por ejecución y copia: fuentes, puntos, segmentos, individuo, copia y
subsemilla.

El repositorio distingue el número de productores con K≥1 del número de copias
finales derivadas de ellos; una salida autozigótica duplicada cuenta una ejecución
y dos copias. Registra además recuentos/frecuencias por vector, genotipos y F
realizados y la procedencia intacta o mosaico. El panel no posee el resultado.

El productor solicita todos sus flujos al servicio DTC-037. La clave de ámbito
combina población, bloque, individuo y, salvo en la decisión compartida de
autozigosis, copia; usa decisiones nombradas distintas para F, primaria, uniones
interbase y secundaria. La versión de derivación forma parte del contrato. La
identidad lógica evita que añadir o reordenar individuos/bloques cambie los ya
existentes, y el ámbito de bloque evita que cambiar `r` contamine otros. Una
salida duplicada por F conserva un único registro de productor enlazado a dos
copias biológicas.

La persistencia materializa tres proyecciones obligatorias y versionadas:
`founder_block_assignments.tsv` indexa población/individuo/copia/bloque y enlaza
HAP final, tipo de origen y proceso. Cada productor recibe un ID DTC-039, incluso
con K=0; el origen es exactamente `intact` para K=0 o `recombined` para K≥1. La
duplicación por F enlaza dos asignaciones al mismo proceso y origen, sin una
categoría especial. `recombination_segments.tsv` ordena por
proceso cada segmento reconstruible con bloque, límites y HAP fuente; y
un proceso K=0 conserva un segmento primario para todo el bloque.
`realized_haplotype_frequencies.tsv` agrega población/bloque/HAP final con
recuento total, intacto y recombinado y `n_h/C`; el total debe ser la suma de los
dos orígenes. El manifiesto enlaza cada una por clase,
versión, ruta, tamaño y SHA-256 y solo resume configuración, F, r, RNG,
productores con cruce, cruces totales, vectores nuevos y coincidencias con panel.

Las tres proyecciones usan UTF-8 sin BOM, LF final, tabuladores, una cabecera
exacta y ningún campo vacío, comilla, línea vacía o columna extra. Sus esquemas
son `population_id individual_id copy_id block_id final_haplotype_id origin_kind process_id`,
`process_id segment_index block_id start end source_haplotype_id` y
`population_id block_id haplotype_id copy_count intact_copy_count recombined_copy_count frequency`;
los espacios de esta descripción son tabuladores. Identificadores, ordinales,
coordenadas, recuentos y `frequency` reutilizan las formas canónicas DTC-039 y
DTC-049. El orden y clave son, respectivamente, población+individuo+copia+bloque,
proceso+índice y población+bloque+HAP. Clase y versión residen en el enlace del
manifiesto/esquema, no en columnas añadidas a los TSV.

La pestaña poblacional de Marimo contiene la sección «Bloques haplotípicos en
población». Su adaptador solo selecciona directorio, invoca preparar/validar y
proyecta el resumen devuelto por Rust; no analiza tablas. Mantiene el mapeo entre
poblaciones del panel y simuladas y una tabla de sustituciones población×bloque.
La unión discriminada muestra muestreo probabilístico/controlado en intacto y r
tipado en recombinante, con controlado desactivado y explicado; F solo aparece en
diploides. La previsualización consume frecuencias, `P(K>0)`, homocigosis esperada
y conflictos `alt_af` calculados por el binario. `run` exige validación vigente y
la vista posterior enlaza/resume las tres tablas sin cargarlas completas.

El repositorio de verdad conserva estados faseados, genealogía, procedencia y
capas de cambio. Creación directa, inicialización poblacional, reproducción,
recombinación, dinámica evolutiva y mutación acceden a él mediante operaciones
de dominio explícitas. El orquestador decide qué componente invocar; los
componentes no se llaman entre sí para alterar el orden funcional.

Las fuentes normativas son REQ-REF-01 a REQ-REF-36, REQ-IND-01 a
REQ-IND-65, REQ-DIR-01 a REQ-DIR-05, REQ-EVO-01 a REQ-EVO-77,
REQ-POB-01 a REQ-POB-03, REQ-CRU-01 a REQ-CRU-45, REQ-REP-01 a
REQ-REP-26, REQ-MIG-01 a REQ-MIG-09, REQ-SEL-01 a REQ-SEL-17,
REQ-REC-01 a REQ-REC-29 y REQ-MUT-01 a REQ-MUT-15. Sus recorridos se
aceptan con ACE-02 a ACE-14, ACE-24 a ACE-32, ACE-75 a ACE-91, ACE-95 a
ACE-98, ACE-103 a ACE-107 y ACE-109 a ACE-117 según la asignación individual
de la matriz de trazabilidad.

### 7.3. Observación técnica

La orquestación de eventos recibe una vista de solo lectura de la verdad y crea
un contexto técnico aislado. Simulación WGS, artefactos, demultiplexado y
alineamiento trabajan sobre ese contexto; los adaptadores de procesos externos
solo traducen invocaciones y resultados, y no adquieren autoridad para cambiar
la verdad.

Para una herramienta consolidada como minimap2, el adaptador conserva la
versión observada, los parámetros efectivos, el comando, entradas y salidas,
código de salida, stdout y stderr. DTC-043 ya fija la ejecución como proceso
externo. DTPD-105 concretará descubrimiento e invocación, sandbox, compatibilidad
de versiones y clasificación de fallos, sin reabrir servicio frente a proceso;
DTC-047 conserva los validadores independientes.

Los FASTQ y demás evidencias nacen en el contexto del evento. Las vistas de
verdad interoperable nacen del repositorio canónico y del contexto de una
referencia, no de inferencias sobre las lecturas. Ambas ramas se reúnen
únicamente al ensamblar y validar el conjunto de resultados.

Esta frontera implementa REQ-TEC-01 a REQ-TEC-10, REQ-RG-01 a REQ-RG-07,
REQ-DAT-01 a REQ-DAT-73, REQ-WGS-01 a REQ-WGS-14, REQ-CONT-01 a
REQ-CONT-16 y REQ-VER-01 a REQ-VER-47. La comprobación funcional se
distribuye entre ACE-15 a ACE-21, ACE-30, ACE-31, ACE-37 a ACE-63,
ACE-66 y ACE-68 a ACE-74.

### 7.4. Cierre y publicación

F1 aplica esta frontera a la ejecución y a sus artefactos reales de
configuración y referencia. No crea manifiestos de evento, conjuntos vacíos ni
evidencias sintéticas para unidades que todavía no existen. Los productores de
F2–F8 extenderán el mismo patrón al introducir cada unidad real.

Cada productor escribe primero cada artefacto en su propio temporal. Solo tras
cerrarlo, validarlo y comprobar tamaño y SHA-256 lo finaliza atómicamente dentro del
directorio provisional de ejecución, acompañado de identidad, referencia,
procedencia y métricas técnicas. Un temporal incompleto se retira o se ignora y
nunca se cataloga como artefacto válido. El ensamblador forma conjuntos
coherentes; los validadores consumen vistas de solo lectura y emiten resultados
estructurados. El escritor recibe hechos ya producidos y genera JSON conforme a
RFC 8259: un `manifest.json` general y uno por evento, validados respectivamente
contra esquemas JSON Schema Draft 2020-12 distintos con definiciones comunes.
Cada documento incluye `manifest_kind` y `schema_version` con la cadena exacta
`"1"`, y no reconstruye decisiones desde nombres o mensajes. El escritor siempre
emite la versión actual conocida.

Los manifiestos contienen directamente todo dato que REQ-MAN exige registrar o
enumerar. Verdad, procedencia de lecturas, informes de validación, complemento
VCF y los demás artefactos independientes conservan sus propios ficheros y se
enlazan sin duplicarlos mediante identificador, ruta relativa, tipo, tamaño en
bytes, algoritmo `sha256` y digest hexadecimal minúsculo de 64 caracteres. Las
listas con semántica de conjunto se escriben en orden estable por identificador.
V1 no usa fragmentos JSONL ni tablas auxiliares para
repartir el contenido obligatorio del manifiesto.

SHA-256 se calcula sobre todos los bytes exactos del fichero, sin normalizar
texto, FASTA, JSON ni saltos de línea. Se calcula o verifica al copiar una
entrada, al cerrar una salida, antes de exportarla y después de exportarla; una
copia debe coincidir con su original. El cálculo puede ser incremental durante
la E/S, pero ninguna comprobación exigida se sustituye por tamaño, fecha,
metadatos o muestreo parcial. Cada fichero conserva su propio digest y una suma
adicional del paquete no lo reemplaza. MD5 solo se añade, por separado, cuando
lo exige un estándar genómico; BLAKE3 queda reservado a los dominios separados
de derivación de subsemillas e IDs.

El manifiesto declara para cada uno de los cuatro niveles de reproducibilidad
`garantizado`, `no garantizado` o `no aplicable`, junto con condiciones,
versiones, perfiles, herramientas y dependencias. La semilla raíz se registra
aunque la haya generado el programa y pasa a ser una entrada reutilizable. V1
no declara identidad binaria por defecto para ninguna clase de salida.

El manifiesto relaciona cada identificador opaco con su clase, ámbito, causa y
etiqueta opcional. La etiqueta y el nombre físico no tienen autoridad de
identidad.

La publicación es la última frontera y distingue tres recorridos. Un rechazo
durante admisión o pre-simulación no publica destino; solo puede dejar el
provisional de trabajo o un diagnóstico externo conforme a REQ-SAL-26. Si la
ejecución ya comenzó y falla de forma controlada, el ensamblador conserva los
artefactos y conjuntos independientes ya validados, los declara exportables,
retira o ignora incompletos y cierra un `manifest.json` general válido con estado
`fallido`, diagnóstico y catálogo de estados de unidades. Solo entonces puede
aparecer el destino final como ejecución fallida, nunca completa. Ante SIGKILL,
caída u otra interrupción no controlable queda el provisional marcado
`INCOMPLETE`; habrá diagnóstico únicamente si alcanzó a escribirse.

Los destinos finales completos o `fallido` se publican mediante un adaptador
Linux propio que realiza una operación atómica no-clobber equivalente a
`renameat2(..., RENAME_NOREPLACE)` dentro del mismo sistema de ficheros. Un
destino existente produce `EEXIST` y rechazo. Un sistema de ficheros que no
ofrezca la garantía se rechaza de forma identificable antes o durante la
publicación; nunca se usa check-then-`std::fs::rename` ni otro fallback con
carrera. Esta garantía lógica no implica `fsync` ni durabilidad física. Los
estados y efectos funcionales de cierre están en REQ-VAL-01 a REQ-VAL-56,
REQ-MAN-01 a REQ-MAN-92 y REQ-SAL-01 a REQ-SAL-27, con ACE-22 a ACE-24,
ACE-33, ACE-34, ACE-36, ACE-46 a ACE-67, ACE-78, ACE-86 y ACE-92 a ACE-108.

En el recorrido local ordinario, el provisional pertenece al proceso y no es
una frontera frente a un atacante con acceso de escritura concurrente. F1 exige
rutas relativas canónicas, permisos privados razonables, temporales propios,
hashes y renombrado no-clobber; no enumera exhaustivamente el árbol para
rechazar hardlinks, inodos externos o cada carrera TOCTOU. Esas defensas se
evalúan en F9/CI si el despliegue declara ese modelo de amenaza.

### 7.5. Construcción y publicación del programa (F9)

El recorrido ordinario de desarrollo se ejecuta de forma nativa en Ubuntu
24.04 sobre WSL2 y siempre dentro del sistema de ficheros Linux/ext4, nunca en
`/mnt/c`. La plataforma nativa oficialmente validada es Ubuntu 24.04 LTS
x86-64 GNU y la construcción de entrega usa el triple
`x86_64-unknown-linux-gnu`, CPU portable y prohíbe `target-cpu=native`. Otras
distribuciones x86-64 GNU pueden funcionar, pero no se prometen hasta probarlas
en una matriz. Ningún contrato depende de WSL.

Al comenzar F1, `rust-toolchain.toml` fija `channel = "1.97.1"`,
`profile = "minimal"`, los componentes `rustfmt` y `clippy` y el target
`x86_64-unknown-linux-gnu`. La edición `2024` se declara en cada `Cargo.toml` o
en el workspace, no en `rust-toolchain.toml`. La aplicación versiona
`Cargo.lock`; CI y entregas usan `--locked`. No existe otra MSRV y las entregas
binarias no requieren Rust en el sistema de la persona usuaria.

Antes de publicar, la candidata se comprueba en CI Ubuntu 24.04 LTS nativo y en
el entorno OCI de referencia, cuya base y receta están fijadas por digest. Las
pruebas cubren la operación no-clobber en el sistema de ficheros OCI admitido,
señales, permisos, límites de recursos y procesos externos. El `.tar.gz` publica
SHA-256 de sus bytes exactos. La imagen OCI se identifica por el digest
`sha256:` de su descriptor de índice o manifiesto y conserva y verifica los
descriptores de configuración y capas; si se exporta como archivo OCI, ese
archivo publica además su SHA-256 independiente. Nada de ello implica identidad
binaria entre compilaciones.

## 8. Invariantes de frontera

Los invariantes de esta sección protegen la dirección de dependencias. Su
contenido funcional procede de los REQ citados; aquí solo se asignan
responsables y puntos de control.

| Invariante arquitectónico | Responsable y frontera de control | Fuente funcional |
|---|---|---|
| Toda entidad cruza fronteras con identidad tipada, no inferida de una etiqueta, ruta o nombre físico | Identidad de ejecución, catálogos y ensamblador de resultados | REQ-CONF-01 a REQ-CONF-14, REQ-TEC-01 a REQ-TEC-10 y REQ-MAN-01 a REQ-MAN-92 |
| Toda coordenada viaja asociada a su referencia y convenio | Catálogo de referencias y servicio de coordenadas | REQ-REF-01 a REQ-REF-36 y REQ-VER-01 a REQ-VER-47 |
| Una referencia reducida conserva origen, orientación, clase y mapa completo; FASTA, mapa y manifiesto se publican o fallan como unidad | Reducción de referencias, catálogo, coordenadas y ensamblador | REQ-REF-24 a REQ-REF-36, REQ-VAL-52/53, REQ-MAN-89/90 y REQ-SAL-27 |
| Una matriz SNV canónica e íntegra se fija antes de crear individuos; todos comparten sus loci/frecuencias y solo genotipos/fase se muestrean por identidad. Preparación y generación integrada/independiente comparten productor; los descartes no se reponen. Un panel se conserva intacto y cualquier mosaico exige el motor común de recombinación | Fuentes de variación, inicialización poblacional y recombinación | REQ-CONF-128 a REQ-CONF-130, REQ-IND-60 a REQ-IND-65 y REQ-EVO-74 a REQ-EVO-77 |
| Solo los motores biológicos escriben en el repositorio de verdad | Interfaces de dominio y repositorio de verdad | REQ-IND-01 a REQ-IND-65, REQ-EVO-01 a REQ-EVO-77 y REQ-CRU-01 a REQ-CRU-45 |
| Los componentes técnicos observan instantáneas y no mutan genotipo, genealogía ni trayectoria | Orquestador de eventos y vistas de solo lectura | REQ-TEC-01 a REQ-TEC-10 y REQ-DAT-01 a REQ-DAT-73 |
| `DeclaredConfig` conserva lo escrito y su ubicación; `EffectiveConfig` es inmutable durante cada ámbito de uso y mantiene la procedencia de la resolución | Adaptador TOML y resolvedor de configuración | REQ-CONF-01 a REQ-CONF-130 y REQ-MAN-01 a REQ-MAN-92 |
| Marimo solo entrega un TOML visible y autoconsistente y argumentos de proceso; cualquier cambio invalida validación y ninguna opción implementada queda solo en TOML. Su servidor solo enlaza con `127.0.0.1`; la versión incompatible del binario bloquea las acciones, sin hashes entre componentes del mismo árbol | Aplicación Marimo y adaptador de proceso | REQ-UI-01 a REQ-UI-15 y REQ-VAL-54 a REQ-VAL-56 |
| Toda probabilidad, frecuencia y tolerancia general es un `f64` finito válido; las colecciones se suman y normalizan una sola vez con Neumaier y orden canónico, y conservan auditoría decimal de ida/vuelta exacta | Admisión de configuración y frontera `DeclaredConfig` → `EffectiveConfig` | REQ-CONF-106/107, REQ-VAL-40, REQ-MAN-81 y ACE-105 |
| Todo manifiesto declara versión y clase, cumple su esquema y conserva directamente el contenido REQ-MAN obligatorio | Escritor de manifiestos y validación de publicación | REQ-MAN-01 a REQ-MAN-92 y REQ-SAL-04/17/18/27 |
| Todo documento estructurado propio usa una clase inequívoca y la versión de esquema conocida `"1"`, incluidos los diagnósticos JSON; una versión distinta se rechaza antes de interpretar su contenido | Adaptadores de entrada y escritores de documentos | REQ-CONF-01/10/111 a REQ-CONF-115 y REQ-MAN-02/69/73 |
| Las decisiones aleatorias consumen ChaCha12 locales derivados con los contextos, material binario y catálogo normativos; ningún consumidor afecta el estado de otro | Proveedor de aleatoriedad y motores consumidores | REQ-CONF-03 a REQ-CONF-08, REQ-CONF-126/129, REQ-MUT-01 a REQ-MUT-15, REQ-EVO-01 a REQ-EVO-77, REQ-REP-01 a REQ-REP-26, REQ-REC-01 a REQ-REC-29 y REQ-DAT-01 a REQ-DAT-73 |
| Toda entrega oficial v1 usa la referencia Ubuntu 24.04 LTS x86-64 GNU y CPU portable y supera validación en Ubuntu nativo y OCI antes de publicarse | Construcción y distribución | REQ-CONF-111 a REQ-CONF-115, REQ-REP-22/23, REQ-MAN-05 y REQ-SAL-19 |
| Cada artefacto se finaliza desde su temporal; un destino completo o `fallido` solo aparece tras manifiestos coherentes y mediante no-clobber atómico, mientras una interrupción no controlable conserva `INCOMPLETE` | Productores, validadores, ensamblador y publicador | REQ-VAL-01 a REQ-VAL-56, REQ-MAN-01 a REQ-MAN-92 y REQ-SAL-01 a REQ-SAL-27 |
| Verdad canónica, vistas interoperables y evidencia técnica permanecen como clases distintas | Repositorio de verdad, materializador VCF y productores técnicos | REQ-VER-01 a REQ-VER-47, REQ-DAT-01 a REQ-DAT-73 y REQ-SAL-01 a REQ-SAL-27 |

Ningún adaptador de formato, proceso externo, manifiesto o escritor de salida
puede saltarse estas fronteras. Los casos observables son los ACE asignados a
cada REQ en
[trazabilidad-requisitos-0.1.md](trazabilidad-requisitos-0.1.md).

## 9. Configuración, persistencia y resultados

### 9.1. Clases de información

La arquitectura separa cuatro almacenes lógicos para impedir que una
representación derive autoridad de otra:

| Clase | Propietario | Mutabilidad y uso |
|---|---|---|
| `DeclaredConfig` | Adaptador TOML | Conserva lo escrito, su procedencia y ubicación; la copia TOML exacta se guarda para auditoría y los motores no la consultan |
| `EffectiveConfig` | Resolvedor de configuración | Vista validada, heredada y resuelta e inmutable en su ámbito; la consumen orquestadores y motores y se serializa versionada en `effective-config.json` |
| Estado interno y verdad | Repositorios de dominio | Evoluciona solo mediante operaciones biológicas; usa DTC-061 y conserva lo necesario para reconstrucción |
| Resultados publicados | Ensamblador, validadores y publicador | Incluyen `truth.jsonl` como proyección canónica obligatoria; VCF y demás vistas se hacen visibles solo después de validarlo y cerrar el conjunto |

Los contenidos exigidos a estas clases están en REQ-CONF-01 a REQ-CONF-130,
REQ-VER-01 a REQ-VER-47, REQ-MAN-01 a REQ-MAN-92 y REQ-SAL-01 a
REQ-SAL-27. TOML fija la configuración de autoría; JSON fija
`effective-config.json` y los manifiestos. Los esquemas de manifiesto usan JSON
Schema Draft 2020-12. Todos estos documentos declaran su clase y el campo
`schema_version` con el valor de cadena `"1"`, usando la sintaxis TOML o JSON
que corresponda, sin confundirlo con la versión del programa ni con el dialecto
de JSON Schema. La salida JSON de diagnósticos es otra clase estructurada propia
y declara `diagnostic_kind: "gps_diagnostics"` y la misma versión de esquema.

### 9.2. Espacios lógicos

La raíz de ejecución contiene el `manifest.json` general y
`effective-config.json`; cada espacio de evento contiene su propio
`manifest.json`. Los recursos compartidos, los eventos y las referencias de
mapeo son espacios lógicos distintos. Las identidades tipadas enlazan esos
espacios; salvo esos nombres confirmados, las rutas físicas son responsabilidad
de los adaptadores de persistencia, que las derivan o escapan desde identidades
sin convertirlas en identidad. Ninguna convención de directorios actúa como
clave biológica y esta arquitectura no fija todavía un layout adicional.

### 9.3. Escritura y recuperación

Los productores pueden escribir por flujo o por fragmentos, pero cada artefacto
empieza en un temporal propio. Solo después de cierre, validación, tamaño y
SHA-256 pasa al catálogo validado del provisional de ejecución. Cada artefacto
atraviesa estados explícitos —provisional, validado, fallido o publicado— y el
manifiesto se construye desde esas transiciones. El tamaño completo de los
manifiestos y de `effective-config.json` entra en el preflight; no permite omitir
campos, enumeraciones ni enlaces. El provisional es hermano del destino, reside
en el mismo sistema de ficheros y lleva el marcador `INCOMPLETE`.

Dentro de F1, una salida que el propio proceso acaba de escribir en ese
provisional privado no se vuelve a tratar como entrada hostil. Se comprueban su
cierre, tamaño y SHA-256 y se conservan pruebas independientes pequeñas de los
formatos; la revalidación semántica exhaustiva del paquete completo y las
matrices de mutaciones físicas pertenecen a la consolidación F9.

Hay tres cierres distintos. Un rechazo durante admisión o pre-simulación no
publica destino y deja, como máximo, provisional o diagnóstico externo. Un
fallo controlado tras comenzar elimina o ignora temporales incompletos, conserva
las unidades independientes ya validadas, cierra un manifiesto general válido
con estado `fallido`, diagnóstico, catálogo de estados y declaración expresa de
qué es exportable, retira el marcador y publica el directorio final como fallido.
Un SIGKILL, caída u otra interrupción no controlable deja el provisional
`INCOMPLETE`; el diagnóstico solo existe si pudo escribirse antes. V1 no reanuda
automáticamente ni reutiliza artefactos parciales; repetir crea una ejecución
nueva.

Cada transición de E/S conserva tamaño y SHA-256 por fichero. La copia de una
entrada verifica igualdad con el original; una salida se verifica al cerrar y
el paquete, antes y después de exportar. Estas verificaciones recorren el
contenido exacto aunque el digest también se haya calculado incrementalmente.

Si una exportación cruza a otro sistema de ficheros, se copia a un provisional
del sistema de destino, se verifican allí tamaño y digest y solo entonces se
publica allí mediante una operación Linux no-clobber equivalente a
`renameat2(..., RENAME_NOREPLACE)`. La misma operación publica dentro del sistema
de ficheros original: `EEXIST` rechaza un destino ya existente y la ausencia de
soporte se diagnostica sin fallback check-then-rename. No se declara atomicidad
entre sistemas de ficheros. La publicación lógica por renombrado no implica
por sí misma
durabilidad ante fallo físico: `fsync` y cualquier garantía de persistencia
deben implementarse y verificarse expresamente antes de prometerlas. La
reanudación por artefacto queda como ampliación futura tras medir costes.

`effective-config.json` declara `schema_version` con `"1"` y el discriminador de
clase equivalente definido por su esquema. Todo escritor emite siempre la
versión actual. V1 no migra documentos ni modifica silenciosamente los leídos;
una eventual conversión entre versiones exigiría una decisión posterior y
produciría otra copia conservando el original.

### 9.4. Entregas del programa

El paquete versionado `.tar.gz` contiene el ejecutable para
`x86_64-unknown-linux-gnu`, los esquemas, las licencias y los ejemplos, además
de la aplicación Marimo, `pyproject.toml`, `uv.lock`, el directorio de contratos
y el lanzador único. No incluye `.venv`, CPython ni `uv`; una instalación
explícita mediante `genome-population-simulator-ui install` prepara el entorno y
`genome-population-simulator-ui run` funciona sin red. La imagen
OCI es la referencia encapsulada de ejecución CLI, fija base y receta por digest
y no contiene Marimo, Python, `uv` ni un puerto UI.
El `.tar.gz` publica SHA-256 de todos los bytes exactos del fichero. Una imagen
OCI abstracta no es un único flujo de bytes: se identifica por el digest
`sha256:` del descriptor de índice o manifiesto y conserva y verifica los
descriptores de configuración y capas. Si se exporta como archivo OCI, ese
fichero publica además su propio SHA-256 independiente. Las recetas no activan
optimizaciones específicas de la CPU del constructor. DTC-043 fija dependencias,
DTC-045 publicación y DTC-048 toolchain y lock; las versiones concretas de crates
se aceptan al implementar F1 conforme a DTC-043 y quedan fijadas por
`Cargo.lock`.

## 10. Errores y observabilidad

### 10.1. Estrategia de errores

Los componentes devuelven errores tipados de su propio vocabulario. El
orquestador añade ejecución, etapa y entidad; los adaptadores añaden ruta,
herramienta o posición cuando corresponda. Una capa puede traducir un error al
vocabulario de la capa superior, pero conserva la causa y no repara ni
reclasifica el resultado por su cuenta.

La frontera pública recibe esos errores y crea un único objeto diagnóstico con
código estable del espacio de nombres del proyecto —por ejemplo,
`GPS-CONF-001`—, categoría, mensaje accionable y todo contexto disponible:
ruta y ubicación, campo e identificadores de entidad, evento o artefacto. La
cadena causal técnica se conserva para detalle y registro, pero ni sus textos ni
los tipos de crates se convierten en el mensaje principal o en el contrato
público. Por defecto, el renderizado humano se escribe en stderr. El modo JSON
emite exactamente un documento RFC 8259 por invocación, nunca JSONL, validado
por su esquema y con el sobre `diagnostic_kind: "gps_diagnostics"`,
`schema_version: "1"` y una lista `diagnostics`. Ambos se producen desde los
mismos objetos.

El catálogo versionado usa `GPS-{categoría}-NNN`, con las categorías `CONF`,
`INPUT`, `FMT`, `DOMAIN`, `RES`, `TOOL`, `IO` e `INTERNAL`. Antes de implementar
F1 debe enumerar todos sus códigos. Un código publicado no se reutiliza ni cambia
de significado.

| Categoría pública mínima | Contexto que debe aportar la arquitectura | Fuente de la consecuencia funcional |
|---|---|---|
| Configuración o esquema | Fichero TOML, clase, versión, campo, valor, ruta y ubicación; en duplicados, todas las ubicaciones; en validación numérica, colección, identidad/categoría, tolerancia, orden, suma, delta, algoritmo/versión y condición incumplida | REQ-CONF-01 a REQ-CONF-130 y REQ-VAL-01 a REQ-VAL-56 |
| Entrada o formato | Ruta, formato, posición, identificador, valor encontrado y contrato incumplido | REQ-FMT-01 a REQ-FMT-07, REQ-REF-01 a REQ-REF-36 y REQ-VAL-01 a REQ-VAL-56 |
| Validación de dominio o biología | Entidad, generación o evento, estado previo identificable e invariante de dominio | REQ-IND-01 a REQ-IND-65, REQ-EVO-01 a REQ-EVO-77, REQ-CRU-01 a REQ-CRU-45, REQ-REP-01 a REQ-REP-26, REQ-MIG-01 a REQ-MIG-09, REQ-SEL-01 a REQ-SEL-17, REQ-REC-01 a REQ-REC-29 y REQ-MUT-01 a REQ-MUT-15 |
| Recursos | Estimación, medida o límite, unidad, ámbito y operación impedida | REQ-NOF-01 a REQ-NOF-08 y REQ-VAL-01 a REQ-VAL-56 |
| Proceso o herramienta externa | Evento, herramienta, versión, comando, código de salida, entradas, salidas y ubicación de stdout o stderr | REQ-TEC-01 a REQ-TEC-10, REQ-DAT-01 a REQ-DAT-73, REQ-UI-01 a REQ-UI-15 y REQ-VAL-01 a REQ-VAL-56 |
| E/S o publicación | Conjunto, referencia, artefactos y enlaces; ante integridad fallida, algoritmo, digest esperado y observado, tamaño y etapa | REQ-VER-01 a REQ-VER-47, REQ-MAN-01 a REQ-MAN-92, REQ-SAL-01 a REQ-SAL-27 y REQ-VAL-01 a REQ-VAL-56 |
| Defecto interno | Componente, operación, invariante imposible y contexto reproducible | REQ-VAL-01 a REQ-VAL-56 |

La clasificación observable entre rechazo, fallo, aviso, resultado estadístico y
estado no evaluable es exclusivamente la de REQ-VAL-01 a REQ-VAL-56 y sus
ACE; esta tabla no introduce otra política.

Los fallos esperables recorren retornos tipados y nunca provocan `panic`. Un
`panic` capturado, una aserción imposible o una invariante interna rota se
clasifican como defecto interno y no como dato inválido de la persona usuaria.
La salida del proceso es estable: 0 para éxito y también para una invocación con
solo avisos; 1 para fallo de ejecución o dominio con configuración válida; 2
para uso, configuración o entrada inválida; 3 para recursos; 4 para herramienta
externa; 5 para E/S o publicación; y 70 para defecto interno.

### 10.2. Registro observable

Los componentes emiten eventos estructurados a un colector común. El sobre
técnico contiene identidad de ejecución, fase, componente, entidad, gravedad y
correlación causal. La carga útil procede del componente responsable. Los
mensajes para personas se renderizan desde ese registro y no son la única fuente
de diagnóstico.

La observación de tiempo, memoria, almacenamiento y cantidades de dominio se
sitúa en las fronteras que producen o consumen esos recursos. La selección de
métricas, política de datos sensibles y umbrales de rendimiento se mantiene en
las decisiones de planificación, sin modificar los criterios funcionales.

Marimo presenta `stdout` y `stderr` de forma incremental durante la sesión. No
los interpreta como verdad biológica ni necesita persistirlos en un almacén
propio. Los resultados duraderos son los manifiestos y artefactos producidos
por Rust; la UI los resume a partir de sus metadatos y enlaces, sin cargar
ficheros genómicos grandes completos.

El registro de reproducibilidad conserva la raíz de 256 bits como 64
hexadecimales minúsculos, su origen declarado o CSPRNG del sistema operativo,
el nombre ChaCha12, la versión y catálogo del esquema de derivación y las
subsemillas de evento en el mismo formato. Los ámbitos internos se reconstruyen
desde IDs y reglas versionadas, sin persistir cada estado mutable. Un cambio de
ChaCha12, BLAKE3, contextos, material binario, catálogo o algoritmo consumidor
es un cambio de reproducibilidad y exige nueva versión cuando afecte a la
derivación.

La evidencia de plataforma de DTC-074 procede de detección, no de etiquetas.
Un colector registra `/etc/os-release`, kernel, arquitectura, triple Rust,
marcadores WSL2, montaje del árbol y metadatos de CI. La expectativa solicitada
solo se compara con esos hechos y una discordancia falla. Ubuntu 24.04 nativo
exige un job identificable sin marcadores WSL; WSL2 exige sus marcadores y un
montaje Linux fuera de `/mnt`. Navegador y herramientas quedan inventariados.
Teclado y lector de pantalla producen actas manuales separadas con revisor,
fecha, versiones, comprobaciones y resultado; ninguna sustituye la evidencia
automática.

La biología solo se declara reproducible con la misma versión del simulador,
plataforma inicialmente admitida, semilla raíz, `EffectiveConfig`, contenidos
de entrada, perfiles y esquema y algoritmos aleatorios versionados. El nivel
técnico exige además la misma subsemilla de evento y las mismas versiones y
condiciones de las dependencias participantes. La equivalencia semántica compara
el contenido lógico y, para herramientas externas, requiere mismas versiones,
parámetros y condiciones verificadas. No se garantiza reproducibilidad entre
versiones del simulador ni entre plataformas no verificadas.

El manifiesto y los metadatos de entrega registran plataforma, arquitectura,
versión del simulador, referencia y digest de la imagen OCI cuando intervenga y
dependencias y herramientas externas efectivas con sus versiones. Registran
además rustc, Cargo, edición y `Cargo.lock` o la resolución efectiva relevante.
La presencia de esos datos permite evaluar las condiciones de reproducibilidad,
pero no declara por sí sola identidad binaria entre compilaciones.

## 11. Estrategia de comprobación

La matriz
[trazabilidad-requisitos-0.1.md](trazabilidad-requisitos-0.1.md) es el índice
obligatorio entre cada REQ, sus ACE, el hito y la evidencia prevista. Esta
sección solo organiza las capas técnicas de comprobación.

| Capa | Propósito técnico | Evidencia |
|---|---|---|
| Dominio unitario | Probar tipos, transformaciones e invariantes sin E/S ni procesos externos | Casos TEST de cada REQ biológico y sus ACE |
| Frontera de componente | Probar traducción de identidades, coordenadas, configuración y errores con dobles de dependencias | Vectores normativos BLAKE3 e IDs; vectores Neumaier con orden canónico, bordes inclusivos, normalización y decimal round-trip; tipos no intercambiables; catálogo diagnóstico y JSON único; comparadores de verdad y contenido técnico; TEST de cada familia |
| Adaptador | Comprobar lectura, escritura e invocación contra contratos externos | TOML desconocido o duplicado con ubicación; clase/versión inválida; FASTA con `samtools faidx`; SHA-256 con `sha256sum`; manifiestos, configuración efectiva y diagnósticos con otra implementación JSON Schema; `.tar.gz` por listado/extracción; descriptores OCI por inspección/arranque; VCF/BAM/CRAM en sus fases; casos válidos, truncados, corruptos y de borde |
| Aplicación Marimo | Demostrar paridad local sin duplicar semántica | Pruebas de controles y recorrido UI→TOML→binario, arranque en `127.0.0.1`, acción explícita, invalidación tras editar, teclado/etiquetas y presentación de un resultado pequeño; la acreditación completa de plataforma y distribución queda en F9 |
| Propiedad y estadística | Explorar espacios amplios sin confundir variación aleatoria con un invariante | Vectores de derivación por ámbito, independencia al añadir o reordenar observaciones y casos TEST asociados a REQ-MUT-01 a REQ-MUT-15, REQ-EVO-01 a REQ-EVO-77, REQ-REP-01 a REQ-REP-26, REQ-REC-01 a REQ-REC-29 y REQ-DAT-01 a REQ-DAT-73 |
| Recorrido completo | Demostrar que las fronteras componen un paquete trazable | Rechazo previo sin destino, fallo controlado final y SIGKILL provisional; carreras no-clobber y FS sin soporte; comparadores semánticos y de bytes; pruebas WSL2/ext4 y publicación en Ubuntu 24.04 LTS nativo y OCI fijado por digest, incluido su FS, señales, permisos, recursos y procesos externos |

Las versiones concretas y los restantes validadores independientes se eligen en
la planificación. La prueba del productor nunca sustituye por sí sola a una
comprobación independiente cuando el hito exige interoperabilidad.

Toda incorporación o ampliación de características de una dependencia adjunta
una evaluación revisable: propósito y requisito cubierto, corrección del
estándar, mantenimiento, seguridad, licencia, rendimiento, estabilidad de API,
*features* activadas, código nativo, scripts de construcción y coste de
distribución. La presencia de `build.rs` no basta para rechazarla: se auditan su
función, entradas, salidas y efectos. Sin aprobación expresa, el script no puede
ejecutar procesos externos ni usar la red, modificar fuentes, generar o compilar
código nativo ni realizar acciones no necesarias o no acotadas. Se comprueba
además que el contrato propio impide filtrar tipos y errores externos. No se
incorpora una dependencia para una necesidad futura aún no implementada.

La evaluación final de `noodles-fasta` comprobada el 10 de agosto de 2026 acepta
la versión exacta 0.66.0: licencia MIT, MSRV 1.89.0 compatible con Rust 1.97.1,
edición 2024 y características mínimas síncronas sin `async`. Es Rust puro y
cubre lectura, escritura, FAI e indexación. El grafo transitivo incorpora BGZF,
Rayon y cuatro `build.rs` de Crossbeam/Rayon; estos scripts fueron auditados y
solo emiten configuración para Cargo, sin procesos externos, red, escritura de
fuentes ni generación o compilación de código. El coste transitivo se acepta y se
vigila. La API upstream cambia con frecuencia; adaptador propio, versión exacta,
lock y pruebas mitigan el riesgo.

Una actualización de Rust o del grafo de dependencias se propone como cambio
aislado. Antes de aceptarla ejecuta pruebas, validadores independientes, vectores
de RNG y reproducibilidad y los escenarios de rendimiento; después vuelve a
declarar las garantías aplicables sin suponer equivalencia con la versión
anterior.

La matriz de CI fija y registra versión y plataforma de cada validador
independiente. Ningún escritor o lector del producto es el único oráculo de su
propia salida. Las pruebas separan validez formal de corrección biológica: la
primera usa implementaciones externas; la segunda, verdades pequeñas conocidas
y comparadores de dominio. Una auditoría externa opcional durante el runtime no
forma parte de v1. Los validadores concretos de formatos posteriores se eligen
en sus DTPD; DTC-048 fija el toolchain y el lock de la aplicación.

F1 conserva una batería sintética estable para observar recursos y crecimiento:

| Escenario | Carga y frecuencia | Propósito |
|---|---|---|
| A | Referencia total de 100 kbp; configuración mínima representativa, semilla, manifiestos e integridad; en cada cambio | Aceptación rápida y corrección del recorrido completo |
| B | El mismo recorrido con referencias deterministas de 1, 10 y 100 Mb | Serie de crecimiento comparable |
| C | Referencia de 1 Gb; ejecución manual o programada, no en cada cambio | Observación grande y detección de cambios de orden de magnitud |

Los tres usan construcción optimizada portable sin `target-cpu=native`,
configuración y semilla fijas, salida limpia y repeticiones suficientes. Se
registran versión, toolchain, plataforma, CPU y memoria, sistema de ficheros,
entorno WSL2, Ubuntu 24.04 LTS nativo u OCI y número de repeticiones. WSL2 sirve
para la tendencia de desarrollo; Ubuntu nativo u OCI es la referencia de publicación.
Estos tamaños son cargas sintéticas de F1 y no prometen escala biológica futura.

## 12. Recursos y evolución prevista

Las interfaces de secuencia, haplotipo, lecturas y artefactos admiten acceso por
tramos y escritura incremental. La representación interna, el tamaño de bloque,
la concurrencia y la política de caché se eligen con medidas; no forman parte
del contrato funcional. El estimador y el observador de recursos comparten
unidades, pero permanecen separados para poder comparar estimación y uso real.
La optimización conserva una línea base de CPU portable para las entregas y no
usa características detectadas del equipo de construcción como parte implícita
del contrato.

Cada escenario mide tiempo total y por etapa, CPU cuando sea fiable, RSS máximo,
bytes de entrada, temporales y salida, MB/s y número de ficheros o artefactos.
Inicialmente no existe un objetivo de producto ni una pequeña variación temporal
convierte el caso en fallo; se conservan línea base y crecimiento observado.
Corrección e integridad sí son obligatorias. Los futuros umbrales y presupuestos
permanecen en DTPD-115 y solo se fijarán con evidencia.

Los puntos de extensión se sitúan en adaptadores de formatos y herramientas,
estrategias de dominio expresamente aprobadas y productores de artefactos. Una
extensión no puede entrar por una rama genérica y eludir identidad, validación,
trazabilidad o versionado.

Las capacidades funcionales fuera de v1 se registran únicamente en
[registro-ampliaciones-futuras-0.1.md](registro-ampliaciones-futuras-0.1.md);
las decisiones técnicas aún pendientes se mantienen en la sección 7 de
[decisiones-tecnicas-0.1.md](decisiones-tecnicas-0.1.md). Esta arquitectura no
anticipa esos comportamientos con entidades implícitas ni los presenta como
obligaciones vigentes.

## 13. Trazabilidad de requisitos a responsables

| Familia de requisitos | Responsable principal | Colaboradores y comprobación |
|---|---|---|
| REQ-CONF | Admisión de configuración | Identidad, coordenadas, aleatoriedad y validación prueban precedencias y rechazos |
| REQ-TEC, REQ-RG | Orquestación de eventos técnicos | WGS, artefactos, alineamiento y manifiestos prueban aislamiento e identidad técnica |
| REQ-FMT | Ensamblado de resultados | Adaptadores de formato y validación externa prueban orden, índice y referencia |
| REQ-REF | Catálogo de referencias | Coordenadas, integridad y exportación prueban FASTA, linaje e índices necesarios |
| REQ-IND, REQ-DIR | Motor de verdad genómica y creación directa | Configuración y materialización prueban ploidía, fase, eventos y reconstrucción |
| REQ-EVO, REQ-POB | Dinámica evolutiva e inicialización poblacional | Reproducción, selección, muestreo y validación poblacional |
| REQ-CRU | Grafo de familias | Reproducción, selección declarativa, verdad y validación del grafo |
| REQ-REP | Reproducción y dinámica evolutiva | Historial y validación comprueban N, Ne, orden y generaciones |
| REQ-SEL | Dinámica evolutiva: selección de viabilidad | Motor genómico y validación comprueban locus, dosis, función y supervivencia |
| REQ-MIG | Dinámica evolutiva: migración | Selección declarativa e historial prueban asignación conjunta y pertenencia |
| REQ-MUT | Motor de verdad genómica | Evolución, familias y validación prueban orden, tasas y transmisión |
| REQ-REC | Reproducción y recombinación | Coordenadas, motor genómico y salida detallada prueban bivalentes y segmentos |
| REQ-DAT, REQ-WGS | Simulación WGS y alineamiento | Artefactos, modo directo y validación técnica prueban origen y coherencia |
| REQ-CONT | Artefactos y demultiplexado | Procedencia y validación prueban sustitución, transferencias y totales |
| REQ-VER | Materialización de verdad | Motor genómico, referencias y validación VCF prueban fase y aislamiento |
| REQ-VAL | Validación interna | Todos los productores aportan objetivos, hechos y contexto sin autocorregirse |
| REQ-MAN | Manifiestos y exportación | Identidad y ensamblado aportan el catálogo de hechos y artefactos |
| REQ-SAL | Ensamblado de resultados | Validación y manifiestos deciden completitud y exportabilidad |
| REQ-NOF | Estimación y observación de recursos | Todos los componentes exponen medidas sin imponer límites funcionales fijos |
| REQ-UI | Aplicación Marimo | Admisión, diagnósticos, proceso y manifiestos prueban equivalencia, seguridad y ausencia de semántica duplicada |

Los criterios ACE-01 a ACE-117 se distribuyen siguiendo esta matriz. Antes de
planificar cada fase, sus criterios se convertirán en casos concretos y se
enlazarán con el componente y la frontera que aportan la evidencia.

## 14. Decisiones abiertas y condición de aprobación

Las decisiones confirmadas, propuestas, pendientes y futuras se registran en
[decisiones-tecnicas-0.1.md](./decisiones-tecnicas-0.1.md). TOML para la
configuración de autoría, JSON para manifiestos y configuración efectiva, RFC
8259 y JSON Schema Draft 2020-12 ya están confirmados. V1 usa únicamente la
versión de esquema `"1"`, rechaza otras versiones y no realiza migraciones. La
aleatoriedad usa ChaCha12 y deriva subsemillas con BLAKE3; esto no selecciona
por sí solo crates. DTC-043 fija su política de admisión y DTC-048 fija Rust
1.97.1 stable con perfil minimal, rustfmt, clippy y target GNU en
`rust-toolchain.toml`; la edición 2024 se declara en Cargo, y `Cargo.lock` se
versiona y usa bloqueado. SHA-256 es la única suma publicada por v1 y permanece
separada de BLAKE3. El FASTA de v1 es texto sin comprimir y su adaptador usa
`noodles-fasta` 0.66.0 tras un contrato propio y una evaluación final favorable
de su grafo y scripts de construcción. La plataforma nativa oficial es Ubuntu
24.04 LTS con triple `x86_64-unknown-linux-gnu`; WSL2 es un entorno de desarrollo
y las entregas se verifican en Ubuntu nativo y OCI. El `.tar.gz` publica SHA-256
de sus bytes; la imagen OCI se identifica por sus descriptores y un archivo OCI
exportado tiene su propio SHA-256.
Los diagnósticos, la publicación atómica, los escenarios iniciales y la
validación independiente quedan fijados por DTC-044 a DTC-047. Esta arquitectura
aplica además DTC-049 para `f64`, suma Neumaier, tolerancia absoluta general y
normalización auditada de probabilidades y frecuencias; DTPD-110 solo conserva
distribuciones y tolerancias poblacionales específicas posteriores. No se
seleccionan silenciosamente otros crates, conversores futuros, formatos aún
abiertos, versiones de estándares genómicos, reglas de paralelismo ni objetivos
cuantitativos de rendimiento. Tampoco amplía la garantía a versiones del
simulador o plataformas no verificadas ni promete identidad binaria.

DTC-050 confirma el contrato funcional de reducción, DTC-056 fija las
coordenadas comunes y DTC-057 el algoritmo uniforme/versionado y su mapa.
DTC-051 confirma Marimo como interfaz sin autoridad biológica. DTC-079
sustituye la concreción sobrediseñada de DTC-058 para F2: entorno bloqueado,
servidor local, TOML general, invocación directa y pruebas de recorrido. La
distribución y el endurecimiento permanecen en F9.

DTC-052 asigna a una matriz común las frecuencias SNV iniciales, DTC-054 resuelve
DTPD-118 con el par TSV+JSON, perfiles, preparación y generación reproducible, y
DTC-053 separa el panel posterior de bloques haplotípicos locales. La matriz
aporta los SNP independientes exteriores y el panel, las combinaciones conjuntas
intrabloque; ningún locus de matriz puede caer dentro del intervalo de un bloque
y una colisión se rechaza antes de crear individuos, sin prioridad ni
sobrescritura. Las asignaciones de alelos y frecuencias del panel son densas y no
admiten completar valores ausentes. El paquete, la preparación y el modo intacto
están definidos. La alternativa recombinante fundadora usa productores directos,
F, puntos, IDs, RNG y salidas trazables ya fijados. DTC-055 resuelve DTPD-119 y
rechaza distribuciones genotípicas explícitas con construcción recombinante.
Ninguna de estas fuentes puede poseer o reescribir la verdad resultante.

DTC-061 resuelve DTPD-102 con la representación interna híbrida y su proyección
obligatoria `truth.jsonl` v1, validada por esquemas de registro y reconstrucción
semántica independiente. DTC-062 fija VCF 4.5 como única versión de salida, la
cabecera `##fileformat=VCFv4.5` y validación externa basada en htslib. DTC-063
fija `.vcf.gz` BGZF con exactamente un índice TBI o CSI como salida
predeterminada: TBI solo si toda secuencia y posición cabe en su límite
normativo, aproximadamente 2^29 bases por secuencia, y CSI en caso contrario;
el modo explícito `.vcf` produce texto sin índice. Configuración efectiva y
manifiesto conservan modo, índice y rutas; BGZF e índice se validan y publican
atómicamente. DTC-064 fija la normalización y agrupación multialélica interna,
con `bcftools norm` solo como oráculo externo. DTC-065 fija `GPS_CASE` como
lista canónica de casos originales, su relación ALT↔casos en `truth.jsonl`,
`GPS_LOSS`, `GPS_COMPLEMENT` y el complemento condicional
`vcf-complement.jsonl`. DTC-066 cierra DTPD-103 con `bcftools` 1.24 y
HTSlib/tabix 1.24 solo en aceptación/CI, más un oráculo semántico independiente
que no comparte la proyección productora. El corpus cubre texto, BGZF+TBI,
BGZF+CSI y una mutación aislada de cada invariante; EBI `vcf-validator` no actúa
como autoridad VCF 4.5.

DTC-067 resuelve DTPD-120 con `poisson-count-v1`, límite `lambda<=2^52`,
conversión uniforme/consumo RNG exactos y `libm` 0.2.16 sujeto a DTC-043.
DTC-068 resuelve fuente y receta de los dos ejemplos humanos compactos
`human_1000g_global_sfs_approx_v1` y
`human_harris2017_snv12_approx_v1`; F2.2.6 valida sus literales, fuente eLife
fijada, licencia, transformaciones y digests sin procesar el VCF Phase 3. La
materialización empírica exacta queda en DTF-016/post-v1 y no bloquea F2.2.
DTC-069 queda aplazada por DTC-079 hasta la primera fase que implemente
escenarios reales; F2 importa y exporta un único TOML general.
DTC-070 resuelve RAF2-029: mayores restos fija las cuotas del modo `exact`, una
permutación independiente por población las asocia a todos los `variant_id` y
otro flujo obtiene el valor sin consumo para `fixed` o cuotas cero. DTC-071
publica los literales del catálogo RNG F2 y resuelve RAF2-028. DTC-072 fija
los consumidores aleatorios no Poisson de F2; DTC-073 queda simplificada por
DTC-079 a mostrar únicamente opciones implementadas por el binario; DTC-074
mantiene en F9 la evidencia oficial de plataforma. DTC-075 conserva los bordes
y la separación de selectores, pero no obliga a crear un expediente ni un paso
independiente para cada primitiva.

Las decisiones base DTPD-001 a DTPD-015 están resueltas. Para cerrar F0 debe
aprobarse esta división; esa aprobación habilita planificar F1, pero no significa
que F1 esté implementada ni cerrada. Para cerrar una fase posterior deben
resolverse las decisiones que cambien sus contratos antes de redactar su plan de
implementación.

La arquitectura se considera comprobable cuando:

- cada requisito de la fase tiene un responsable y una evidencia prevista;
- ninguna frontera permite mezclar referencias, identidades o ámbitos;
- verdad y evidencia técnica recorren caminos separados;
- los fallos y resultados parciales tienen un estado no ambiguo;
- un plan de implementación puede dividir el trabajo sin inventar nuevos
  componentes, contratos o decisiones técnicas.
