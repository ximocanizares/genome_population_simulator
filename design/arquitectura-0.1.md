# Arquitectura del simulador de genomas y poblaciones

- **Versión:** 0.1
- **Estado:** borrador para revisión
- **Fecha:** 7 de agosto de 2026
- **Especificación de origen:**
  [especificaciones-0.1.md](./especificaciones-0.1.md)
- **Hoja de ruta:** [hoja-de-ruta-0.1.md](./hoja-de-ruta-0.1.md)
- **Registro técnico:**
  [decisiones-tecnicas-0.1.md](./decisiones-tecnicas-0.1.md)
- **Lenguaje del núcleo:** Rust, según la especificación y la hoja de ruta 0.1

## 1. Finalidad y estado

Este documento define cómo se divide el sistema y qué contratos deben mantener
sus partes para satisfacer la especificación 0.1. No define archivos de código,
paquetes o bibliotecas de Rust, denominados *crates*, formatos de configuración,
algoritmos concretos ni un plan de
implementación.

La arquitectura cubre toda la versión 1 para evitar contratos incompatibles
entre fases. Sin embargo, su materialización seguirá la secuencia F1–F9 de la
hoja de ruta. Una capacidad futura no autoriza a implementarla antes de su fase.

La especificación y la hoja de ruta confirman Rust para el núcleo. Esta elección
no decide formatos, *crates*, versión del lenguaje, algoritmos, arquitectura ni
plataformas.

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
FASTQ, BAM o CRAM requieren tecnología, perfil y bibliotecas; solo los que
producen alineamientos requieren además un mapeador. Por tanto, un evento puede
materializar únicamente verdad genotípica o VCF sin parámetros técnicos de
secuenciación.

### 2.2. Fronteras externas

Son fronteras del sistema:

- los ficheros estructurados de configuración general y poblacional;
- FASTA de referencia y de elementos transponibles;
- FASTQ suministrados para contaminación o calibración;
- perfiles empíricos reutilizables;
- los estándares FASTA, FASTQ, SAM/BAM, CRAM y VCF;
- minimap2 como mapeador predeterminado y los mapeadores compatibles que se
  autoricen en el futuro;
- las herramientas externas necesarias para mapear, comprimir, ordenar,
  indexar o validar formatos, cuando la decisión técnica correspondiente las
  seleccione;
- el sistema de ficheros y los recursos de CPU, memoria y almacenamiento.

Cada integración externa se realiza mediante un contrato sustituible. El
dominio biológico no conoce comandos, rutas concretas, cabeceras físicas ni
detalles de procesos externos.

### 2.3. Fuera de alcance arquitectónico

No forman parte del sistema los *SNP callers*, la comparación de sus resultados,
la imputación, el faseado inferido ni la evaluación de herramientas externas.
Tampoco se incluyen las capacidades funcionales excluidas en la sección 3 de la
especificación, entre ellas señales de instrumento, BCL, paneles dirigidos,
Oxford Nanopore, PacBio CLR, aneuploidía individual, reproducción clonal,
conversión génica y transposición espontánea.

Python no forma parte del núcleo. Solo podrá incorporarse más adelante como
herramienta auxiliar de calibración o análisis si se aprueba una decisión
técnica que defina su frontera, distribución y reproducibilidad.

## 3. Principios arquitectónicos

1. **La verdad biológica precede a la evidencia técnica.** Lecturas,
   alineamientos, contaminación, errores y calidades consultan la verdad, pero
   nunca la modifican.
2. **Cada dato lleva identidad y procedencia.** Referencia, sistema de
   coordenadas, individuo, haplotipo, evento y conjunto no se deducen de una
   ruta o de una posición implícita.
3. **Las referencias no se mezclan.** Un conjunto de resultados pertenece a
   una única referencia de mapeo y sus identificadores de variante son locales.
4. **La configuración se resuelve antes de simular.** Se aplican valores
   predeterminados y precedencias, se validan incompatibilidades y se obtiene
   una configuración efectiva inmutable para la ejecución.
5. **Los estados inválidos no cruzan fronteras.** Los identificadores,
   coordenadas, ploidías, proporciones y relaciones se representan mediante
   valores validados o producen errores explícitos.
6. **La aleatoriedad tiene ámbitos independientes.** La evolución, los
   muestreos y cada evento explícito no comparten consumo accidental de una
   única secuencia aleatoria. La reproducibilidad se declara por niveles y
   alcanzar uno no presupone alcanzar el siguiente.
7. **Los estados y la publicación son jerárquicos.** Ejecución, evento, conjunto
   y artefacto distinguen como mínimo `incompleto`, `completo` y `fallido`. Cada
   artefacto se publica atómicamente solo después de validarse y una unidad
   superior solo está completa cuando lo están todas sus unidades obligatorias.
8. **Los recursos se observan antes de optimizar.** La arquitectura expone
   puntos de medida y admite representaciones compactas posteriores, sin fijar
   objetivos cuantitativos no acordados.

## 4. Dirección de dependencias

```text
Entradas y herramientas externas
             │
             ▼
Adaptadores de configuración, formatos y procesos
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

## 5. Modelo conceptual del dominio

### 5.1. Identidades

Cada identidad es estable y pertenece a una categoría. Identidades de
categorías distintas no son intercambiables aunque su texto coincida.

| Entidad | Identidad y ciclo de vida | Relaciones obligatorias |
|---|---|---|
| Ejecución | Se crea tras validar y resolver la configuración. No se reutiliza | Versión, entradas, configuración efectiva, semilla raíz, estados y manifiesto general |
| Referencia | Estable durante la ejecución | FASTA exacto, secuencias ordenadas, clasificación cromosoma/*scaffold*, origen y cambios si es derivada |
| Secuencia de referencia | Su identificador procede del primer elemento no vacío de la cabecera FASTA y es único dentro de su referencia | Descripción restante de cabecera, longitud, clasificación, ploidía efectiva y convenio de coordenadas |
| Segmento | Identifica una porción orientada y trazable | Estado actual, origen fundador o evento que lo creó, orientación y correspondencias disponibles |
| Variante o evento genómico | Identidad estable en el historial biológico | Clase, orden, marco de coordenadas previo, haplotipo afectado y estado final |
| Individuo | Único y no reutilizable en la ejecución | Referencia estructural, generación, población, parentales, haplotipos y observaciones |
| Copia haplotípica | Distinta aunque tenga la misma secuencia que otra | Individuo, secuencia, índice de copia, ploidía, segmentos y fundador |
| Población y subpoblación | Identidades no reutilizables | Generación, miembros, estado activa/extinguida y relaciones demográficas |
| Evento biológico | Identifica muestreo, cruce, mutación o cambio demográfico según su categoría | Entradas, orden, resultado y causas de omisión o fallo |
| Evento explícito | Unidad de observación independiente | Individuos, referencias, subsemilla y salidas; tecnología, perfil, bibliotecas y mapeador solo cuando procedan |
| Experimento de secuenciación | Condición de un evento que genera FASTQ, BAM o CRAM; no es otra identidad de evento | Tecnología y perfil; bibliotecas y, si genera alineamientos, mapeador |
| Biblioteca, tanda y carril | Únicos dentro del ámbito definido por el esquema futuro | Evento, muestra, tecnología, perfil e índices técnicos |
| Grupo de lectura | Único en toda la ejecución | Evento, individuo, biblioteca, tanda o carril y archivos asociados |
| Instancia de lectura | Observación concreta con identidad técnica única dentro del evento; una pareja constituye una unidad y comparte el nombre exigido por el formato | Biblioteca actual, miembros cuando sea pareja, origen estable y mecanismo de copia o transferencia |
| Conjunto | Local a un evento y una referencia de mapeo | Mismas muestras que los demás conjuntos del evento y resultados de esa referencia |
| Variante VCF | Local a un conjunto | Referencia, registro VCF y complemento de verdad si existe |

La sintaxis y el mecanismo de generación de identificadores quedan pendientes.
El contrato exige unicidad, estabilidad, categoría y ausencia de reutilización.

### 5.2. Referencias, coordenadas y segmentos

Una referencia posee una secuencia ordenada de cromosomas o *scaffolds* y una
identidad verificable del FASTA exacto. Una referencia derivada enlaza una única
referencia de origen y la secuencia ordenada de cambios estructurales que la
creó. El catálogo de referencias es un grafo dirigido sin ciclos durante una
ejecución.

En un FASTA, el identificador de secuencia es el primer elemento no vacío que
sigue a `>` y termina en el primer espacio o carácter de separación. Debe ser
único dentro de la referencia. El resto de la cabecera se conserva como
descripción completa y no participa en la unicidad.

Toda secuencia suministrada o generada usa exclusivamente `A`, `C`, `G`, `T` y
`N`, sin distinguir mayúsculas de minúsculas durante la validación. `N` expresa
una base desconocida o no determinada. Su tratamiento operativo debe decidirse
funcionalmente para cada clase de operación; el sistema no puede resolverla al
azar, sustituirla ni atribuirle otro significado de forma implícita.

Toda coordenada contiene, como mínimo, referencia, secuencia, posición o
intervalo y convenio. La configuración utiliza base 1 e intervalos cerrados;
los formatos externos conservan su convenio oficial. La representación interna
se decidirá por separado. Ningún componente puede convertir coordenadas sin
usar el servicio común de coordenadas y registrar la conversión.

Un segmento conserva la secuencia actual o la forma de materializarla y su
procedencia. Puede apuntar a una copia fundadora, a una copia duplicada o a una
inserción sin origen fundador. Una deleción conserva el intervalo ausente. Las
correspondencias entre estados solo existen donde la transformación las hace
significativas; nunca se inventan para cubrir regiones sin homología.

### 5.3. Estado genotípico y verdad

El estado genotípico de un individuo agrupa, para cada secuencia, exactamente el
número de copias que indica su ploidía efectiva. Cada copia contiene una cadena
ordenada de segmentos y eventos suficiente para materializar su secuencia,
conservar fase y dosis y alcanzar un origen fundador o un evento de inserción.

La **verdad biológica** comprende referencias, haplotipos, genotipos,
parentescos, poblaciones, transmisiones, mutaciones y observaciones. Es la única
fuente autorizada para generar el fichero obligatorio de verdad y los VCF.

La **verdad técnica de lecturas** conserva el origen real de fragmentos,
secuencias original y observada, errores, duplicados, contaminaciones,
quimeras, descartes agregados y categoría sintética cuando corresponda. Esta
verdad se deriva de la biológica y del evento técnico, pero no la sustituye.

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

El estado de una subpoblación distingue N objetivo, N realizado, Ne objetivo o
ausente y Ne realizado. Distingue también subpoblación activa, extinción
planificada y extinción por resultado. Una recolonización crea una identidad
nueva.

### 5.5. Eventos y resultados

Un evento de muestreo crea una observación no destructiva del estado biológico.
Todo evento es explícito y puede producir solo verdad, solo VCF o evidencia
técnica. Cuando genera FASTQ, BAM o CRAM constituye además un experimento de
secuenciación. Cada evento produce cero o más conjuntos, uno por referencia
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
| Admisión de configuración | Leer, combinar, validar y resolver configuración general y poblacional | Ficheros y entradas declaradas → configuración efectiva o diagnósticos | Posee modelos de entrada y configuración efectiva; consulta catálogos de capacidades |
| Identidad de ejecución | Crear y comprobar identidades tipadas y el contexto común | Configuración validada → identidad y contexto | Posee el registro de identidades; no posee datos biológicos |
| Aleatoriedad reproducible | Proporcionar flujos independientes por ámbito | Semilla y claves de ámbito → flujo aleatorio auditable | Posee la jerarquía de ámbitos y el registro de semillas; no interpreta decisiones biológicas |
| Catálogo de referencias | Validar, copiar, simular y derivar referencias | FASTA o configuración → referencias disponibles | Posee identidades, secuencias, linaje y recursos exactos de referencia |
| Coordenadas y transformaciones | Validar intervalos y convertir convenios sin ambigüedad | Coordenadas etiquetadas → coordenadas validadas o correspondencias | Posee reglas de conversión; consulta referencias y segmentos |
| Motor de verdad genómica | Crear y transformar haplotipos y eventos | Referencia, cambios y flujos aleatorios → estado faseado | Posee el estado genotípico, eventos, segmentos y correspondencias |
| Creación directa | Traducir la vía directa a operaciones del motor genómico | Configuración efectiva → individuos | Posee solo el estado de progreso; el motor posee la verdad resultante |
| Inicialización poblacional | Crear fundadores, frecuencias, genotipos, LD y F inicial | Configuración poblacional → población inicial | Posee estadísticas solicitadas/realizadas; entrega verdad al repositorio biológico |
| Reproducción y recombinación | Formar gametos, bivalentes y descendencia trazable | Progenitores, mapas/tasas y flujo → descendientes | Posee eventos de transmisión; consulta verdad parental |
| Dinámica evolutiva | Ejecutar generaciones y procesos en orden válido | Población y fases → población final e historial | Posee estado generacional y resúmenes; delega cambios genómicos |
| Grafo de familias | Validar y ejecutar el grafo declarativo de cruces | Grafo y parentales → poblaciones intermedias y finales | Posee dependencias y resultados del grafo; consulta selección y reproducción |
| Selección declarativa | Evaluar el catálogo común sobre verdad biológica | Reglas y elegibles → seleccionados y auditoría | No posee genotipos; posee evaluación y resultado de selección |
| Muestreo y resúmenes | Crear observaciones no destructivas y ámbitos estadísticos | Generación y regla → muestra, resumen poblacional y muestral | Posee observaciones y resúmenes; consulta verdad e historial |
| Orquestación de eventos | Validar salidas y parámetros condicionales y aislar cada observación | Observación/configuración → plan de evento y conjuntos | Posee estado del evento y relaciones entre artefactos; solo crea un plan técnico si procede |
| Simulación WGS | Generar fragmentos e instancias de lectura desde haplotipos | Verdad, perfil y objetivo → lecturas con identidad y origen | Posee lecturas previas al demultiplexado y verdad técnica asociada |
| Artefactos y demultiplexado | Convertir conjuntamente la contaminación y aplicar errores, duplicados, adaptadores y reasignaciones | Lecturas y configuración → FASTQ finales y procedencia | Posee transformaciones técnicas e identidades de instancia; consulta fuentes compatibles |
| Alineamiento | Ejecutar el mapeador sobre el FASTQ final | FASTQ y referencia → alineamientos observados | Posee invocación, salida, diagnóstico y metadatos del proceso |
| Alineamiento directo | Crear registros sintéticos coherentes sin FASTQ | Verdad y perfil directo → alineamientos sintéticos | Posee categorías y registros sintéticos; no llama a un mapeador |
| Materialización de verdad | Escribir verdad obligatoria, VCF y complemento | Estado faseado y referencia de conjunto → artefactos de verdad | No altera el dominio; posee la representación materializada |
| Ensamblado de resultados | Organizar recursos comunes, eventos y conjuntos | Artefactos validados → estructura de resultados | Posee catálogo lógico de artefactos y estados de completitud |
| Validación interna | Comparar configuración, verdad y evidencia | Objetivos y resultados → informe, aviso o error | Posee evidencias de comprobación; no corrige resultados |
| Manifiestos y exportación | Registrar procedencia y crear paquetes autocontenidos | Catálogo de artefactos y ejecución → manifiestos/paquete | Posee representaciones de auditoría, enlaces e integridad |
| Estimación y observación de recursos | Estimar y medir tiempo, memoria y almacenamiento | Configuración/ejecución → estimaciones y medidas | Posee mediciones; no decide modelos biológicos |

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

### 6.3. Dependencias permitidas y prohibidas

- Creación directa, inicialización, evolución y familias pueden depender del
  motor genómico, aleatoriedad, identidad y selección declarativa.
- Reproducción es compartida por evolución y familias; no se duplica por vía.
- Simulación WGS, artefactos, alineamiento y modo directo pueden consultar
  verdad e identidades, pero no llamar a componentes que cambien genotipos.
- Materialización VCF depende del estado faseado y de la referencia del
  conjunto; no depende de lecturas ni alineamientos.
- Validación depende de vistas de solo lectura de configuración, verdad y
  evidencia. No puede normalizar o reparar archivos después del hecho.
- Manifiestos reciben hechos ya registrados por los responsables. No
  reconstruyen decisiones a partir de nombres de fichero o mensajes de log.
- Ningún componente de dominio depende de crates o herramientas concretas. Los
  adaptadores futuros implementarán los contratos aprobados.

### 6.4. Comprobación aislada

Cada componente admite dobles de prueba para sus dependencias externas. Las
pruebas del dominio usan referencias y poblaciones pequeñas construidas de forma
controlada. Los adaptadores de formato se prueban con casos válidos e inválidos
del estándar elegido. Los procesos externos se prueban con respuestas de éxito,
salida inválida, código de error, interrupción y versión incompatible.

| Componente | Errores propios principales | Comprobación aislada mínima |
|---|---|---|
| Admisión de configuración | Sintaxis, campo desconocido, precedencia conflictiva y referencia inexistente | Resolver una configuración mínima y rechazar conflictos con diagnósticos que señalen ambos ámbitos |
| Identidad de ejecución | Identidad duplicada, categoría incorrecta o reutilización | Crear identidades por categoría y demostrar colisiones y referencias cruzadas inválidas |
| Aleatoriedad reproducible | Ámbito desconocido, clave duplicada o estado no regenerable | Repetir decisiones por ámbito y demostrar que añadir un muestreo no cambia el flujo evolutivo |
| Catálogo de referencias | FASTA inválido, identificador de primer elemento repetido, símbolo distinto de A/C/G/T/N, linaje cíclico o derivación imposible | Separar identificador y descripción, validar el alfabeto y reconstruir una referencia pequeña desde su origen y cambios |
| Coordenadas y transformaciones | Convenio ausente, intervalo fuera de límites o correspondencia inexistente | Probar extremos, inserciones, inversiones y conversiones de ida y vuelta donde sean definidas |
| Motor de verdad genómica | Evento solapado prohibido, objetivo eliminado, punto de ruptura o ploidía inválidos | Materializar haplotipos conocidos y comparar exactamente secuencia, fase, dosis e historial |
| Creación directa | Individuo incompleto o cambio incompatible | Crear individuos pequeños con cambios concretos y aleatorios controlados |
| Inicialización poblacional | Frecuencias, recuentos, F, LD o restricciones incompatibles | Enumerar casos pequeños y comprobar sumas, fase y resultados por modo |
| Reproducción y recombinación | Ploidía impar, parental incompatible, mapa incompleto o bivalente imposible | Usar casos mendelianos con puntos de recombinación conocidos y origen reconstruible |
| Dinámica evolutiva | Orden inválido, población activa vacía no autorizada o evento demográfico imposible | Ejecutar pocas generaciones enumerables y comprobar estado antes/después de cada proceso |
| Grafo de familias | Ciclo real, fuente futura, cardinalidad o salida imposible | Validar grafos válidos ramificados y ciclos transitivos mínimos |
| Selección declarativa | Regla desconocida, conjunto ambiguo o elegibles insuficientes | Evaluar Y/O/NO y cantidades exactas/proporcionales sobre un catálogo de verdad fijo |
| Muestreo y resúmenes | Duplicado, estratos solapados, tamaño excesivo o ámbito estadístico confundido | Comparar resumen poblacional y muestral y demostrar observación no destructiva |
| Orquestación de eventos | Evento implícito, parámetro técnico indebido, tecnología múltiple, perfil o mapeador incompatible y referencia ajena | Construir eventos solo de verdad y eventos técnicos independientes, y comprobar requisitos condicionales y aislamiento |
| Simulación WGS | Objetivo ausente/doble, perfil incompatible, identidad repetida o fragmento imposible | Comprobar recuentos, cobertura haploide, aporte por copia e identidades únicas antes de sesgos |
| Artefactos y demultiplexado | Conversión conjunta inviable, exceso de contaminación, colisión de identidad o par huérfano | Resolver una asignación entera conocida y rastrear copias biológicas, transferencias, duplicados y descartes |
| Alineamiento | Herramienta ausente, versión incompatible, código no cero o salida inválida | Sustituir el proceso externo y comprobar comando, contexto y propagación de cada fallo |
| Alineamiento directo | Proporciones, MAPQ, banderas, CIGAR o pares incoherentes | Generar un caso por categoría y validarlo contra la referencia y verdad de origen |
| Materialización de verdad | Alelo no reconstruible, normalización ambigua o pérdida de fase | Reconstruir cada caso desde VCF y complemento y compararlo con el haplotipo final |
| Ensamblado de resultados | Artefacto activado ausente o mezcla de referencias | Intentar cerrar conjuntos completos e incompletos y verificar estados resultantes |
| Validación interna | Evidencia insuficiente, invariante rota o desviación fuera de política | Inyectar resultados conocidos y comprobar conforme, aviso y error sin corrección automática |
| Manifiestos y exportación | Enlace roto, identidad incoherente o integridad fallida | Recorrer todos los enlaces de un evento mínimo y verificar el paquete antes/después |
| Estimación y observación de recursos | Dato no medible, estimación inválida o presupuesto excedido | Comparar estimación y medida en casos registrados sin presentar la estimación como exacta |

## 7. Flujos de datos y ejecución

### 7.1. Admisión y preparación

1. Se leen los ficheros declarados sin iniciar la simulación.
2. Se validan sintaxis, categorías, identificadores, referencias cruzadas,
   ploidías, ámbitos, reglas y dependencias.
   Para FASTA se separan el primer elemento de cabecera, usado como
   identificador, y su descripción; para toda secuencia se valida el alfabeto
   `A`, `C`, `G`, `T` y `N`. Si una operación aún no tiene definida su conducta
   ante `N`, no debe inventarla durante la validación.
3. Se aplican predeterminados y precedencias para crear una configuración
   efectiva. Se conservan valor original, valor efectivo y procedencia.
4. Se validan y copian las entradas que deben conservarse. La ejecución usa las
   copias aprobadas y no modifica los originales.
5. Se genera o acepta la semilla raíz y se registran los ámbitos aleatorios.
6. Se estiman recursos. Un presupuesto excedido rechaza; una estimación grande
   sin presupuesto excedido advierte.
7. Solo entonces se crea el contexto de ejecución preparado.

Un error de esta fase no produce una ejecución marcada como iniciada o completa.

### 7.2. Referencias y verdad biológica

1. El catálogo crea o incorpora la referencia inicial.
2. Cada referencia derivada se produce desde una referencia ya disponible y
   conserva el orden de transformaciones.
3. La vía elegida crea individuos directamente, inicializa y evoluciona una
   población, o ejecuta el grafo de familias.
4. Todas las vías entregan el mismo contrato de individuos faseados y
   reconstruibles.
5. Se materializa el alcance obligatorio de verdad de cada vía: todos los
   individuos directos; población inicial, final y generaciones o muestreos
   solicitados en evolución; y todos los individuos y poblaciones intermedias de
   familias. Las generaciones evolutivas restantes conservan resúmenes.
6. El estado completo requerido y los historiales resumidos se entregan a
   materialización y persistencia sin pasar por el recorrido técnico.

### 7.3. Generación evolutiva

El orquestador resuelve los eventos demográficos en su orden explícito y aplica
las dependencias funcionales de REQ-REP-24. Migraciones salientes se calculan
desde el mismo estado inicial y se aplican conjuntamente. La selección de
viabilidad precede a la reproducción. La segregación y recombinación preceden a
las mutaciones de novo. Registro y muestreo observan el estado resultante.

Una extinción por resultado autorizada omite de forma atómica los eventos
dependientes y conserva una cadena causal. Si no está autorizada, detiene toda la
ejecución. Las ramas independientes solo continúan en el primer caso.

### 7.4. Admisión de un evento explícito

Todo evento fija identidad, individuos, referencias, salidas y subsemilla. Si
solo solicita verdad o VCF, la orquestación no exige ni inventa tecnología,
perfil, biblioteca o mapeador. Si solicita FASTQ, BAM o CRAM, valida tecnología,
perfil y bibliotecas. Si solicita alineamientos, valida además un único mapeador
y aplica CRAM como formato predeterminado dentro del evento, salvo elección de
BAM. Sin eventos no se ejecuta este flujo.

### 7.5. Evento de secuenciación en modo realista

1. El evento fija individuos, tecnología, perfil, bibliotecas, referencias,
   mapeador y subsemilla.
2. La simulación WGS selecciona fragmentos desde los haplotipos en proporción a
   dosis, ploidía y longitud antes de sesgos.
3. Después de conocer el número limpio realizado de unidades de cada biblioteca,
   todas las contribuciones de contaminación se convierten conjuntamente a
   lecturas o parejas completas. La asignación entera válida minimiza la
   desviación global; su total es el entero más próximo a la suma convertida y
   los empates se resuelven con la subsemilla. Si no existe solución, el evento
   falla antes de modificar los FASTQ finales.
4. Se aplican perfiles, errores, calidades, artefactos, contaminación biológica
   y demultiplexado en el orden funcional acordado. Una copia biológica recibe
   identidad técnica nueva y enlaza su unidad de origen; una transferencia por
   demultiplexado conserva la misma identidad y cambia de biblioteca. Solo una
   duplicación configurada puede crear otra instancia adicional.
5. Se conserva un FASTQ final o pareja por grupo de lectura y un fichero de
   procedencia sin duplicar secuencias ni calidades.
6. El mismo FASTQ final se entrega al único mapeador del evento contra cada
   referencia seleccionada.
7. Cada referencia produce un conjunto independiente. Por individuo y conjunto
   se combina un único alineamiento final con sus grupos de lectura.
8. El resultado se ordena, indexa y valida según las decisiones de formato.

### 7.6. Evento de secuenciación en modo directo

El modo directo omite FASTQ y mapeador. Genera registros sintéticos desde la
verdad conocida, aplica las mismas reglas funcionales de cobertura,
contaminación y reasignación, y conserva categoría, origen y parámetros. Sus
campos, banderas, CIGAR, MAPQ y relaciones de pareja deben ser coherentes con la
categoría. El manifiesto lo identifica inequívocamente como sintético.

Las identidades técnicas de lectura o pareja siguen las mismas reglas de
unicidad, copia y transferencia que en el modo realista, aunque no exista FASTQ.

### 7.7. Verdad por referencia

Para cada conjunto, el materializador compara directamente los haplotipos de
las muestras con la referencia de mapeo de ese conjunto. Crea un VCF multimuesta
local y, solo cuando haga falta, un complemento estructurado. Nunca traduce un
VCF de otra referencia y nunca consulta FASTQ, BAM o CRAM para decidir un
genotipo.

Dos eventos con las mismas muestras, referencia y verdad deben producir el
mismo resultado lógico de verdad aunque cambie su evidencia técnica. La
estrategia para lograr identificadores locales estables queda pendiente.

### 7.8. Validación, cierre y exportación

1. Cada productor construye su artefacto en un área temporal y mantiene su
   estado `incompleto`.
2. Se validan estructura, referencia, identidades técnicas, índices y
   relaciones.
3. Las comprobaciones estadísticas registran objetivo, observado, tolerancia y
   decisión. N y Ne nunca se elevan a error por su desviación aislada.
4. Cada artefacto validado se publica atómicamente y pasa a `completo`; si
   falla, pasa a `fallido` y no se publica una versión parcial.
5. Un conjunto solo queda completo con todos sus artefactos activados. Un evento
   solo queda completo con sus salidas propias y conjuntos obligatorios. La
   ejecución solo queda completa cuando lo están todas las salidas mínimas y
   activadas, eventos y conjuntos obligatorios.
6. Un fallo superior no borra artefactos o conjuntos completos e independientes;
   pueden exportarse declarando el fallo global, sin presentar el paquete como
   ejecución completa.
7. La exportación copia los recursos comunes necesarios y vuelve a validar la
   integridad antes y después del paquete.

## 8. Invariantes de frontera

### 8.1. Identidad y referencia

- Un identificador pertenece a una sola categoría y no se reutiliza.
- Cada dato posicional identifica referencia, secuencia y convenio.
- Cada conjunto contiene una única referencia de mapeo.
- Los conjuntos de un evento contienen las mismas muestras biológicas.
- Los identificadores de variante no enlazan registros entre referencias.
- Todo CRAM se relaciona con el FASTA exacto contra el que se mapeó.

### 8.2. Ploidía, fase y procedencia

- La ploidía efectiva de una secuencia es un entero de 1 a 10 e igual para todos
  los individuos de la ejecución.
- La reproducción sexual genérica solo acepta ploidía par en v1.
- Cada individuo tiene exactamente tantas copias por secuencia como indique su
  ploidía.
- Toda variante se asigna a copias concretas y la dosis deriva de esa
  asignación.
- Cada segmento heredado alcanza una copia fundadora concreta; una inserción
  nueva alcanza su evento creador.
- Una mutación de novo no borra el origen del segmento sobre el que apareció.

### 8.3. Transformaciones genómicas

- Los eventos explícitos se aplican en orden y usan el estado inmediatamente
  anterior como marco.
- Un objetivo eliminado produce error; no se recoloca ni se ignora.
- Las variantes estructurales anidadas o que modifican otra estructura previa en
  el mismo linaje se rechazan en v1.
- Adyacencia no equivale a solapamiento.
- Los límites terminales, puntos de ruptura, orientación y máximo local de
  copias se validan antes de publicar verdad.

### 8.4. Aleatoriedad y reproducibilidad

- Toda decisión aleatoria pertenece a un ámbito identificable.
- Cambiar o retirar un muestreo no cambia la trayectoria evolutiva.
- Dos eventos técnicos tienen flujos independientes aunque compartan individuo.
- La misma semilla, configuración, entradas y versión reproduce las mismas
  decisiones biológicas dentro de la garantía finalmente aprobada.
- No se promete identidad bit a bit, compatibilidad entre versiones,
  reproducibilidad entre plataformas ni independencia del paralelismo hasta
  resolver las decisiones correspondientes.

### 8.5. Verdad y evidencia

- La evidencia técnica nunca modifica verdad, fase o genotipo.
- VCF se deriva exclusivamente de la verdad faseada.
- Cobertura cero, contaminación o mal mapeo no convierten un genotipo verdadero
  en ausente.
- Los campos observacionales no se incorporan como genotipo inferido al VCF de
  verdad.
- Toda desviación técnica intencionada se explica mediante configuración,
  perfil, artefacto o decisión aleatoria registrada.

### 8.6. Integridad de resultados

- Un artefacto activado que falta impide el cierre con éxito.
- Un alineamiento final, su índice y su referencia se validan conjuntamente.
- VCF, compresión, índice, cabecera y referencia se validan conjuntamente.
- Los FASTQ de un evento no se reutilizan en otro evento.
- El fichero de procedencia no duplica secuencias ni calidades.
- Una salida provisional o fallida nunca aparece como completa en el
  manifiesto.

## 9. Configuración, persistencia y resultados

### 9.1. Cuatro clases de información

1. **Configuración declarada:** lo escrito por la persona usuaria, conservado
   sin reinterpretación retrospectiva.
2. **Configuración efectiva:** resultado validado de predeterminados,
   precedencias y normalizaciones; queda fijada antes de su ámbito de uso.
3. **Estado interno:** referencias, verdad, eventos e historial necesarios para
   continuar y reconstruir la simulación.
4. **Resultados publicados:** artefactos validados, manifiestos e informes que
   otras herramientas pueden consumir.

Los formatos físicos de estas cuatro clases pueden ser distintos. La elección
queda pendiente y deberá permitir versionado de esquema y lectura por programa.

### 9.2. Organización lógica

La raíz de la ejecución contiene configuración efectiva, verdad obligatoria,
manifiesto general y recursos comunes. Cada referencia y sus índices se
almacenan una vez como recurso común verificable. Cada evento posee un ámbito
separado y, dentro de él, un ámbito por referencia de mapeo para BAM, CRAM, VCF e
índices. Los FASTQ y su procedencia permanecen en el ámbito del evento y grupo
de lectura.

Esta organización es lógica. Los nombres de directorio, extensiones concretas,
esquemas y formato de paquete no quedan decididos por este documento.

### 9.3. Escritura y recuperación

Las escrituras de artefactos grandes pueden ser incrementales, pero su
publicación es atómica desde la perspectiva del manifiesto: o el artefacto
validado está disponible o figura como fallido/incompleto. La política concreta
de archivos temporales, reanudación y limpieza queda pendiente.

## 10. Errores y observabilidad

### 10.1. Categorías

| Categoría | Ejemplos | Efecto |
|---|---|---|
| Configuración | Identidad duplicada, ploidía inválida, ciclo, orden ambiguo, proporciones incompatibles | Rechazo antes de simular |
| Entrada | FASTA/FASTQ ilegible, símbolos inválidos, pares incoherentes | Rechazo del ámbito antes de generar resultados dependientes |
| Dominio biológico | Objetivo eliminado, estructura incompatible, parentales insuficientes | Detención del evento biológico o ejecución según el requisito |
| Resultado aleatorio previsto | N o Ne realizado distinto del objetivo | Conservación del resultado y advertencia; N=0 aplica su política específica |
| Extinción autorizada | Dependencia desaparecida por N=0 permitido | Omisión atómica y causal de la rama dependiente |
| Recursos | Presupuesto configurado excedido o agotamiento durante ejecución | Rechazo previo o fallo explícito; nunca resultado parcial completo |
| Proceso externo | Ejecutable ausente, versión incompatible, código no cero, salida inválida | Fallo del evento/conjunto y de la ejecución si la salida era obligatoria |
| Formato e integridad | Cabecera, índice, referencia o suma incompatibles | El artefacto o conjunto no se marca completo |
| Estadístico | Desviación fuera de tolerancia | Advertencia predeterminada o error configurable, salvo excepciones N/Ne |
| Defecto interno | Invariante imposible rota | Detención segura con contexto; no se intenta continuar silenciosamente |

### 10.2. Registro observable

El sistema emite eventos estructurados de progreso y diagnóstico, además de
mensajes comprensibles. Cada registro contiene ejecución, fase, componente,
entidad, gravedad y contexto. Las invocaciones externas conservan herramienta,
versión, parámetros efectivos, código de salida y diagnósticos, con una política
pendiente para datos sensibles o excesivamente voluminosos.

Se observan, como mínimo, tiempo, memoria, almacenamiento, recuentos de
individuos, haplotipos, variantes, lecturas y artefactos por las fronteras que
los producen. No se fijan todavía umbrales cuantitativos ni objetivos de
rendimiento.

## 11. Estrategia de comprobación

### 11.1. Pruebas de dominio

- ejemplos pequeños que reconstruyen cada haplotipo base por base;
- propiedades de ploidía, fase, dosis, orden y unicidad;
- transformación y composición permitida de variantes;
- casos mendelianos y poliploides controlados;
- poblaciones pequeñas con resultados enumerables;
- invariancia de trayectoria al añadir o retirar muestreos;
- selección declarativa sobre verdad, nunca sobre evidencia;
- aislamiento entre referencias y entre eventos técnicos.

### 11.2. Pruebas estadísticas

Las pruebas estadísticas comparan distribuciones y agregados en múltiples
réplicas y registran semilla, tamaño del caso, estadístico y tolerancia. No
exigen exactitud a una realización estocástica. Las tolerancias y tamaños de
muestra se decidirán antes de implementar cada modelo.

### 11.3. Pruebas de contratos externos

Cada formato se comprueba con un validador independiente cuando exista una
opción aprobada. Se prueban cabeceras, coordenadas, orden, compresión, índices,
referencias, grupos de lectura, pares y casos no mapeados. La selección de
validadores y sus versiones queda en el registro técnico.

### 11.4. Recorridos completos

La hoja de ruta define el crecimiento del recorrido: F1 prueba configuración,
semilla, identidad y referencia; F2 añade verdad faseada; F3 añade FASTQ,
mapeo y alineamiento; F4 amplía realismo técnico; F5–F8 amplían biología; F9
consolida escala e interoperabilidad.

Cada cierre demuestra una verdad independiente del productor. Un archivo que
solo puede validarse volviendo a leer el mismo estado o mediante el mismo
algoritmo no basta para detectar errores biológicamente plausibles.

## 12. Recursos y evolución prevista

Las secuencias y haplotipos se exponen mediante contratos que permiten evitar
copias completas y materializar por tramos, pero la representación concreta se
elige con medidas. La E/S admite escritura secuencial y procesamiento por
unidades independientes sin exigir que todo resida en memoria. Estas son
capacidades de frontera, no una promesa de rendimiento.

La arquitectura deja puntos de extensión para nuevos perfiles, mapeadores,
formatos, modelos mutacionales y tecnologías. Cada extensión requiere una
capacidad funcional aprobada, un adaptador o estrategia con identidad versionada
y pruebas de compatibilidad. No se construyen hoy contratos para tecnologías o
modelos excluidos de v1.

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
| REQ-TE | Motor genómico y catálogo de referencias | Materialización de verdad y manifiestos prueban secuencia, evento y linaje |
| REQ-REC | Reproducción y recombinación | Coordenadas, motor genómico y salida detallada prueban bivalentes y segmentos |
| REQ-DAT, REQ-WGS | Simulación WGS y alineamiento | Artefactos, modo directo y validación técnica prueban origen y coherencia |
| REQ-CONT | Artefactos y demultiplexado | Procedencia y validación prueban sustitución, transferencias y totales |
| REQ-VER | Materialización de verdad | Motor genómico, referencias y validación VCF prueban fase y aislamiento |
| REQ-VAL | Validación interna | Todos los productores aportan objetivos, hechos y contexto sin autocorregirse |
| REQ-MAN | Manifiestos y exportación | Identidad y ensamblado aportan el catálogo de hechos y artefactos |
| REQ-SAL | Ensamblado de resultados | Validación y manifiestos deciden completitud y exportabilidad |
| REQ-NOF | Estimación y observación de recursos | Todos los componentes exponen medidas sin imponer límites funcionales fijos |

Los criterios ACE-01 a ACE-106 se distribuyen siguiendo esta matriz. Antes de
planificar cada fase, sus criterios se convertirán en casos concretos y se
enlazarán con el componente y la frontera que aportan la evidencia.

## 14. Decisiones abiertas y condición de aprobación

Las decisiones confirmadas, propuestas, pendientes y futuras se registran en
[decisiones-tecnicas-0.1.md](./decisiones-tecnicas-0.1.md). En particular, esta
arquitectura no decide silenciosamente crates, formatos o esquemas, versiones de
estándares, plataformas, algoritmo aleatorio, garantías entre plataformas,
reglas de paralelismo ni objetivos cuantitativos de rendimiento.

Para cerrar F0 deben aprobarse esta división y, como mínimo, resolverse las
decisiones marcadas como bloqueantes de F1. Para cerrar una fase posterior deben
resolverse las decisiones que cambien sus contratos antes de redactar su plan de
implementación.

La arquitectura se considera comprobable cuando:

- cada requisito de la fase tiene un responsable y una evidencia prevista;
- ninguna frontera permite mezclar referencias, identidades o ámbitos;
- verdad y evidencia técnica recorren caminos separados;
- los fallos y resultados parciales tienen un estado no ambiguo;
- un plan de implementación puede dividir el trabajo sin inventar nuevos
  componentes, contratos o decisiones técnicas.
