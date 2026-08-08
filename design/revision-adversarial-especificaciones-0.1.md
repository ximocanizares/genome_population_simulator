# Revisión adversarial de las especificaciones 0.1

- **Documento revisado:** `design/especificaciones-0.1.md`
- **Revisión de la especificación:** 0.1
- **Fecha de la revisión:** 7 de agosto de 2026
- **Commit de referencia:** `c9892a8`
- **SHA-256 del documento revisado:** `8d72785a27cc6903a6979aac46f7a0e1d966d0ab7e769549188eac9301d4d18b`
- **Estado:** informe para revisión; no modifica ni sustituye la especificación

Este informe es diagnóstico y derivado. La cobertura REQ→ACE→fase→prueba se
mantiene en
[trazabilidad-requisitos-0.1.md](./trazabilidad-requisitos-0.1.md).
Las exclusiones funcionales recuperables se centralizan en
[registro-ampliaciones-futuras-0.1.md](./registro-ampliaciones-futuras-0.1.md).
Los ataques, alternativas y decisiones históricas que siguen explican cómo se
cerraron ambigüedades, pero no constituyen requisitos ni compiten con la sección
4 de la especificación. Si una formulación histórica difiere del contrato
vigente, se corrige este informe y prevalece el REQ correspondiente.

## 1. Objetivo y método

Esta revisión intenta romper el contrato funcional antes de implementarlo. No
evalúa si la idea del simulador es adecuada, sino si dos equipos independientes
podrían implementar la especificación y obtener comportamientos compatibles.

Se han buscado:

- contradicciones entre requisitos;
- términos que admiten resultados observables diferentes;
- casos límite sin resultado definido;
- obligaciones que no pueden verificarse con los criterios actuales;
- dependencias funcionales presentadas como decisiones meramente técnicas;
- diferencias entre la especificación y los documentos de diseño vigentes.

La severidad usada es:

- **Bloqueante:** debe resolverse antes de cerrar F0 o antes de congelar el
  contrato afectado. Permite implementaciones funcionalmente incompatibles o
  impide determinar qué salida es correcta.
- **Alta:** puede aplazarse hasta antes de la fase indicada, pero no debe llegar
  a su implementación sin una decisión explícita.
- **Media:** no suele impedir el primer código, pero debilita la trazabilidad,
  la validación o la evolución del documento.

## 2. Resultado general

La especificación contiene una base funcional extensa y útil. Tras incorporar
las 27 decisiones y auditar de nuevo los hallazgos, no permanece abierta ninguna
cuestión funcional para v1.

La comprobación mecánica ofrece resultados positivos:

- existen 788 definiciones `REQ` y no hay identificadores duplicados;
- toda mención con forma `REQ-…` corresponde a un requisito definido;
- los criterios `ACE-01` a `ACE-111` forman una secuencia completa, sin saltos
  ni duplicados;
- la separación entre verdad genotípica y evidencia técnica está expresada de
  forma consistente en la mayor parte del documento;
- la independencia de los conjuntos referidos a genomas distintos está
  reiterada y resulta inequívoca como intención general.

Los riesgos funcionales ya tienen contratos observables. RAD-028 queda cerrado
tras auditar semánticamente las 788 filas de trazabilidad, RAD-031 por el
registro auditado de ampliaciones y RAD-029 tras subordinar los resúmenes y
documentos derivados a los REQ sin mantener contratos paralelos.

## 3. Hallazgos bloqueantes

### RAD-001. La elección de Rust no está incorporada a la fuente funcional

**Evidencia:** las secciones 4.1, 4.7 y 9 afirman que el lenguaje todavía no se
ha elegido. La hoja de ruta 0.1 y la arquitectura 0.1 declaran Rust como lenguaje
del núcleo.

**Ataque:** una planificación que tome la especificación como fuente principal
puede volver a abrir la elección de lenguaje o justificar decisiones de formato
como si Rust aún no estuviera confirmado.

**Decisión necesaria:** actualizar la especificación para registrar Rust como
decisión confirmada y conservar como pendientes únicamente las decisiones que
realmente sigan abiertas. Distinguir la versión 0.1 del documento de la versión
1 del producto.

**Bloquea:** F0.

### RAD-002. La sección 8 cierra decisiones funcionales que siguen abiertas

**Evidencia:** la sección 8 declara que no queda ninguna decisión funcional
abierta. Sin embargo, REQ-REP-07 exige límites o distribuciones de descendencia sin fijar su
semántica mínima; REQ-CONF-44 permite órdenes alternativos sin enumerar todas
las dependencias; y diversos requisitos exigen «otras distribuciones» u «otros
modelos admitidos» sin definir un catálogo mínimo.

**Ataque:** dos implementaciones pueden cumplir literalmente el texto y, ante el
mismo escenario, una continuar con una advertencia, otra fallar y otra rechazar
una opción por no pertenecer a su catálogo.

**Decisión necesaria:** reabrir la sección 8 y trasladar a ella toda decisión
que cambie entradas aceptadas, salidas o errores observables. Mantener en la
sección 9 solo elecciones que no alteren el contrato funcional.

**Estado: resuelto.** La decisión 17 cerró los catálogos probabilísticos y de
repeticiones citados; la decisión 20 retiró los límites o distribuciones de
descendencia parental; y la decisión 25 eliminó los órdenes generacionales
alternativos de v1. La sección 8 conserva explícitas las cuestiones que siguen
abiertas y la sección 9 se limita a elecciones técnicas.

**Bloquea:** cierre de la especificación y F0.

### RAD-003. Cero eventos técnicos es incompatible con los alineamientos predeterminados

**Evidencia:** REQ-CONF-101 admite cero eventos de muestreo con salidas técnicas.
REQ-CONF-64 y el objetivo general activan CRAM y el modo realista de forma
predeterminada. A la vez, tecnología, perfil, biblioteca, individuos, mapeador y
subsemilla pertenecen a cada evento.

**Ataque:** una ejecución directa con individuos, sin eventos declarados y sin
desactivar alineamientos no tiene dónde obtener la selección de muestras ni los
parámetros técnicos necesarios para producir el CRAM predeterminado.

**Decisión necesaria:** elegir una regla única entre estas posibilidades:

1. toda salida técnica requiere al menos un evento explícito;
2. el sistema crea un evento implícito con valores definidos;
3. cero eventos desactiva las salidas técnicas aunque el valor general sea
   predeterminado.

También debe definirse si un conjunto puede existir solo con VCF cuando los
alineamientos estén desactivados.

**Bloquea:** F0 y F3.

### RAD-004. La cobertura depende de una «referencia seleccionada» incompatible con FASTQ compartidos

**Evidencia:** REQ-WGS-04 define cobertura como bases de lecturas divididas por
la longitud haploide de la referencia seleccionada. REQ-DAT-22 exige usar los
mismos FASTQ contra todas las referencias de mapeo seleccionadas. Las referencias
pueden tener longitudes distintas por deleciones, duplicaciones o
translocaciones, y la referencia estructural del individuo también puede ser
distinta.

**Ataque:** un único FASTQ no puede alcanzar simultáneamente 30× respecto a dos
referencias de longitudes diferentes si 30× determina el número de bases que se
generan. Tampoco está claro si el denominador es la referencia estructural del
individuo, la referencia inicial o cada referencia de mapeo.

**Decisión necesaria:** fijar un denominador biológico único anterior al mapeo
y definir cómo se calcula cuando los haplotipos o cromosomas tienen longitudes y
ploidías diferentes. La cobertura observada contra cada referencia de mapeo debe
tratarse como una medida realizada, no como varios objetivos simultáneos, salvo
que se generen eventos y FASTQ distintos.

**Bloquea:** F0, F3 y F4.

**Estado: resuelto.** El objetivo usa como denominador la longitud haploide de
la referencia estructural del individuo, sumando cada secuencia una sola vez y
sin ploidía. Haplotipos y referencias de mapeo no lo cambian. El mismo FASTQ se
mapea contra todas las referencias y su cobertura alineada es una observación,
no otro objetivo.

### RAD-005. El alcance del fichero obligatorio de verdad es contradictorio

**Evidencia:** REQ-SAL-02 y la tabla de salidas hablan del estado completo de
«los individuos y poblaciones que correspondan» o incluso de «todos los
individuos y poblaciones». REQ-REP-19 indica que, en evolución, solo se conservan
por defecto los genotipos completos de la población inicial, final y
generaciones solicitadas. REQ-CRU-15, en cambio, obliga a conservar todos los
individuos intermedios de los cruces.

**Ataque:** una implementación puede guardar todos los genotipos históricos y
otra únicamente los hitos indicados por REQ-REP-19; ambas pueden alegar
conformidad. El tamaño y la capacidad de reconstrucción del fichero de verdad
serían radicalmente distintos.

**Decisión necesaria:** enumerar por vía y por clase de entidad qué estados
completos son obligatorios, cuáles son resúmenes y qué significa exactamente
«reconstruible» cuando el estado completo no se conserva. Diferenciar verdad
biológica, historial de eventos y materialización completa del genoma.

**Bloquea:** F0, F1, F7 y F8.

### RAD-006. Las referencias derivadas no pueden a la vez cambiar solo por cuatro SV e incorporar TE

**Evidencia:** REQ-REF-09 limita la diferencia de una referencia estructural a
duplicaciones, deleciones, inversiones y translocaciones. REQ-TE-03 exige poder
insertar elementos transponibles al construir referencias. Una inserción de TE
no pertenece necesariamente a ninguna de las cuatro clases enumeradas.

**Ataque:** una referencia derivada con un TE insertado es obligatoria según
REQ-TE-03 e inválida según REQ-REF-09. Clasificar toda inserción TE como
duplicación tampoco funciona cuando la secuencia procede de un catálogo externo
o es nueva.

**Decisión necesaria:** aclarar si los TE solo pueden formar parte de una
referencia inicial simulada, si pueden crear referencias derivadas o si
«referencia estructural» debe admitir también inserciones configuradas. La lista
de operaciones debe ser cerrada y coherente en todas las secciones.

**Bloquea:** F2 y F5.

**Estado: resuelto.** V1 elimina los elementos transponibles como clase
funcional. No admite catálogos, familias, mecanismos, destinos, identidad ni
metadatos TE. Una región preexistente copiada es una duplicación segmentaria
genérica; una secuencia nueva es una inserción genérica; un FASTA real se trata
como bases sin anotación TE. Todo modelado específico se aplaza.

### RAD-007. [RESUELTO] No se define la biología de los símbolos IUPAC ambiguos

**Evidencia:** REQ-REF-14 admite R, Y, S, W, K, M, B, D, H, V y N en un FASTA
real. Después, el simulador debe crear haplotipos, mutaciones, lecturas, alelos
REF/ALT y validaciones desde esa secuencia.

**Ataque:** ante una `R`, una implementación puede conservar `R` en las lecturas,
resolverla una vez como A/G, excluir la posición o resolverla de nuevo para cada
lectura. Esas opciones generan verdades, errores y VCF diferentes. Una matriz de
mutación tampoco puede elegir las tres bases alternativas sin una base de origen
concreta.

**Decisión necesaria:** definir para cada operación si las bases ambiguas se
preservan, se enmascaran, se excluyen del espacio válido o se resuelven mediante
una decisión aleatoria trazable. Debe quedar definido qué REF se permite emitir
y cómo afectan a cobertura y variantes.

**Estado: resuelto.** V1 solo admite A/C/G/T/N y la pregunta 27 fija `N` como
base desconocida que siempre se conserva. SNP, selección y motivos repetidos la
excluyen; las operaciones estructurales pueden transportarla; WGS la representa
con Phred 0, GC la excluye y cobertura la cuenta; VCF nunca emite un SNP con
`REF=N` ni inventa una base y usa complemento cuando haga falta.

**Bloquea:** F2 y F3.

### RAD-008. Falta una política para colisiones durante la generación aleatoria

**Evidencia:** los eventos estructurales no pueden solaparse o anidarse en el
mismo linaje; SNP e *indels* no pueden tocar puntos de ruptura; algunas variantes
pequeñas sí pueden solaparse y los eventos explícitos se aplican en orden. No se
define qué ocurre cuando dos eventos aleatorios, o uno aleatorio y otro
explícito, compiten por una posición incompatible.

**Ataque:** para una muestra de Poisson igual a 20 en una región con pocos sitios
válidos, una implementación puede remuestrear hasta obtener 20, reducir el
número realizado, fallar o cambiar el orden de las clases. Cada política altera
la distribución, la semilla consumida y la verdad final.

**Decisión necesaria:** fijar el orden entre clases y entre eventos explícitos y
aleatorios; el algoritmo funcional de rechazo o remuestreo; el límite de
intentos; y la consecuencia cuando no existen suficientes sitios válidos. Debe
indicarse qué número se registra como muestreado, aceptado y finalmente activo.

**Bloquea:** F2 y F5.

**Estado: resuelto.** En v1 se aplican primero los eventos explícitos en su
orden declarado. Después se sortean los recuentos por clase y se procesan, en
orden fijo, duplicaciones, deleciones, inversiones, translocaciones,
repeticiones en tándem, *indels* y SNP. Cada evento se elige solo entre destinos
válidos en el estado actual y actualiza ese espacio; no se reintentan candidatos
inválidos ni existe un límite de intentos. Si faltan destinos para completar el
recuento, la generación falla sin reducción silenciosa ni advertencia y registra
los números sorteado, aceptado antes del fallo y finalmente activo.

### RAD-009. Precedencia entre ámbitos configurables

**Evidencia:** se definen precedencias global, por secuencia, región e individuo,
y en otros lugares por población, subpoblación, fase, generación y locus. Se
afirma que se aplicarán las precedencias «correspondientes», pero no se define
un orden total ni cómo resolver ejes cruzados. Tampoco se fija la precedencia
entre el fichero general y el poblacional.

**Ataque:** un individuo de una subpoblación tiene una sustitución individual,
su locus está dentro de una región específica y la fase redefine el mismo
parámetro. No puede saberse cuál prevalece. Dos esquemas de configuración pueden
resolverlo de forma distinta sin contradecir una regla local.

**Decisión aceptada:** REQ-CONF-116 a REQ-CONF-123 definen regla base,
selectores por familia, jerarquías internas y dominancia parcial. No hay
prioridad implícita entre ejes: varias máximas incomparables se rechazan aunque
coincidan en valor. También se rechazan duplicados exactos, regiones solapadas y
selectores no permitidos; REQ-VAL-51 y ACE-111 hacen observable el contrato.

**Estado:** Resuelto.

### RAD-010. N objetivo y Ne objetivo pueden quedar funcionalmente sin efecto

**Evidencia:** REQ-REP-08 exige intentar crear N descendientes, pero REQ-REP-26
acepta cualquier N realizado inferior causado por un resultado aleatorio y solo
advierte. REQ-REP-14 exige aproximar Ne objetivo, mientras REQ-REP-15 establece
que cualquier desviación de Ne solo advierte y nunca invalida.

**Ataque:** una implementación defectuosa podría producir sistemáticamente una
fracción pequeña de N o ignorar Ne objetivo y limitarse a registrar advertencias.
Los criterios actuales no fijan cuánto esfuerzo constituye «intentar» ni cómo
distinguir una desviación legítima de no aplicar el objetivo.

**Decisión necesaria:** definir condiciones verificables de intento y de
imposibilidad, junto con invariantes del generador reproductivo. Para Ne debe
definirse una política de aproximación comprobable, aunque la desviación final
siga sin ser un error. Los avisos no deben servir para ocultar que el algoritmo
no aplicó el objetivo.

**Bloquea:** arquitectura del motor poblacional y F7.

**Estado: resuelto.** N objetivo fija exactamente los intentos reproductivos;
cada éxito crea un descendiente y cada fallo conserva una causa. Los parentales
pueden reutilizarse sin límite. Una imposibilidad previa se rechaza y una
aleatoria conserva los descendientes, la desviación y el aviso; N=0 aplica su
política específica. Ne objetivo se rechaza en v1 y Ne realizado se registra
cuando es evaluable.

### RAD-011. El VCF debe contener toda la verdad y a la vez puede delegarla al complemento

**Evidencia:** REQ-VER-22 exige que el propio VCF exponga sin pérdida motivo,
copias por alelo o haplotipo, longitud, pureza e interrupciones. REQ-VER-24 a
REQ-VER-26 permiten que un complemento conserve la información que no pueda
representarse completamente en VCF.

**Ataque:** una repetición compleja puede considerarse inválida porque no cabe
por sí sola en VCF o válida porque el complemento contiene el detalle. Ambas
lecturas están respaldadas por requisitos obligatorios.

**Decisión necesaria:** definir qué subconjunto mínimo debe estar siempre en el
VCF, qué puede residir exclusivamente en el complemento y qué información debe
duplicarse. Sustituir «el propio VCF» por una regla inequívoca sobre el par
VCF–complemento si esa es la intención.

**Bloquea:** F2, F3 y F5.

**Estado: resuelto.** El VCF conserva una representación estándar válida e
interoperable y expone todo lo que admita sin pérdida la versión oficial
seleccionada. Cuando falta detalle, todos los registros aplicables enlazan
mediante campos `INFO` declarados la única entrada complementaria del caso. El
complemento, no el manifiesto, conserva lo perdido; la verdad estructurada
obligatoria sigue siendo la fuente canónica completa.

### RAD-012. No está definido si una población puede mezclar referencias estructurales

**Evidencia:** REQ-IND-02 permite que cada individuo corresponda a cualquiera de
las referencias disponibles. Selección, LD, migración, mapas de recombinación y
cruces referencian loci o coordenadas de una referencia estructural. Los cruces
solo exigen referencias «compatibles», término no definido funcionalmente.

**Ataque:** una subpoblación contiene individuos basados en dos referencias con
una inversión o translocación diferencial. No está definido qué referencia
describe a la descendencia, cómo se identifica el mismo locus seleccionado,
qué mapa se usa ni qué genotipo satisface una regla declarativa.

**Decisión necesaria:** elegir entre prohibir en v1 la mezcla de referencias
estructurales dentro de una unidad reproductiva o definir completamente su
compatibilidad, la referencia de la descendencia, la identidad de loci y los
mapas aplicables. La posibilidad de mapear las mismas muestras contra varias
referencias no obliga a mezclar referencias biológicas en una población.

**Bloquea:** F0, F6, F7 y F8.

### RAD-013. No se define la atomicidad de una ejecución fallida

**Evidencia:** muchos requisitos ordenan detener la ejecución e impedir que un
conjunto se marque como completo. Otros permiten varios eventos y conjuntos y
la exportación independiente. La validación de un FASTA debe fallar «sin
generar resultados», aunque el flujo copia configuraciones y otras entradas en
pasos anteriores.

**Ataque:** falla el último conjunto después de que otros ya estén completos. No
se sabe si la ejecución entera es inválida, si los conjuntos previos son
exportables, qué archivos temporales se conservan ni qué manifiesto explica el
fallo. En otro caso, una copia temprana de configuración contradice la promesa
de no generar resultados ante un FASTA inválido.

**Decisión necesaria:** definir estados de ejecución, evento y conjunto; qué
unidad es atómica; qué artefactos parciales pueden conservarse; cómo se marcan;
y qué puede exportarse tras un fallo. Distinguir un área temporal de las salidas
publicadas.

**Bloquea:** F0 y F1.

## 4. Hallazgos de severidad alta

### RAD-014. El espacio de oportunidades de las tasas de mutación no está definido

REQ-MUT-01 expresa tasas por base y generación, pero no indica si se aplican a
la longitud haploide, a todas las copias, a las longitudes actuales después de
SV, o solo a posiciones válidas para cada clase. Tampoco define cómo cambia el
denominador tras deleciones y duplicaciones. Debe fijarse antes de F7 para que la
media esperada sea reproducible y comparable.

**Estado: resuelto.** Cada clase cuenta sus
oportunidades válidas en todas las copias haplotípicas descendientes actuales y
usa tasa×oportunidades como media Poisson. Ploidía, duplicaciones y deleciones
modifican el espacio; destinos condicionados no multiplican el denominador. La
pregunta 27 define además cuándo `N` queda excluida sin reabrir el denominador.

### RAD-015. La selección tiene casos matemáticos y de identidad de locus sin resultado

REQ-SEL-13 normaliza por la aptitud máxima, pero no resuelve el caso en que todas
las aptitudes sean cero. Además, un locus ligado a una referencia estructural
puede estar eliminado, duplicado o ser ambiguo en individuos basados en otra
estructura. Deben definirse la respuesta a máximo cero y la semántica de
genotipo o dosis del locus en esos casos antes de F7.

**Estado: resuelto.** V1 usa un único SNP bialélico identificado y protegido,
presente una vez por copia aplicable, y exige valores para todas las dosis. Los
eventos incompatibles se rechazan o excluyen del espacio aleatorio. Una
identidad, dosis o aptitud inválida, o un máximo cero, detiene antes de
supervivencia y no se interpreta como extinción.

### RAD-016. Las unidades de contaminación no tienen conversión ni redondeo común

REQ-CONF-75 admite proporción, cobertura equivalente o número de lecturas o
parejas por fuente. REQ-CONT-02 y REQ-CONT-03 validan una fracción total de
0–100 %. Falta el momento y denominador de conversión, el redondeo a unidades
enteras, la resolución conjunta de varias unidades y la relación exacta entre
diagonal y valores fuera de la diagonal de la matriz. Debe resolverse antes de
F4.

### RAD-017. No se distingue siempre el objetivo anterior de la cobertura final

Los descartes técnicos reducen la cobertura; la contaminación biológica es
sustitutiva; los duplicados y el sesgo GC pueden alterar bases y distribución.
El documento usa «cobertura objetivo», «realizada» y «final» sin una tabla que
indique en qué etapa se mide cada una. Deben definirse al menos cobertura WGS
básica, posterior a sesgos, posterior a contaminación/demultiplexado y observada
tras mapeo, sin convertirlas en un único valor.

**Estado: resuelto.** Se registran por separado la cobertura objetivo, la WGS
básica anterior a sesgos y artefactos, la final del FASTQ después de todos los
efectos técnicos con bases propias y contaminantes separadas y la alineada
observada para cada referencia de mapeo.

### RAD-018. La identidad de lecturas puede colisionar o reutilizarse sin regla

Los FASTQ se separan por grupo, después se combinan en un BAM por individuo. La
contaminación biológica copia lecturas sin reducir al donante y la técnica las
transfiere. No se exige unicidad de identificadores de lectura en el evento ni
se define si una lectura donante puede copiarse a varias receptoras. Deben
fijarse el ámbito de unicidad de `QNAME`, la política de renombrado, el vínculo
con el identificador original y los límites de reutilización antes de F4.

### RAD-019. La garantía de reproducibilidad técnica es ambigua

REQ-CONF-27 garantiza la misma simulación biológica, no archivos idénticos.
Otros requisitos exigen que errores, lecturas, contaminación y alineamientos
sintéticos usen semillas reproducibles. El mapeador externo puede además variar
por versión, plataforma, hilos o metadatos. Debe separarse expresamente:

- identidad de decisiones biológicas;
- identidad de lecturas y artefactos técnicos;
- equivalencia semántica de BAM, CRAM y VCF;
- identidad binaria de archivos.

También debe decirse qué cambios no relacionados deben mantener estables las
subsemillas. Este contrato debe cerrarse antes de F1 y F3.

### RAD-020. «Otros modelos admitidos» no define capacidad mínima

REQ-CONF-25, REQ-DIR-04, REQ-CONF-61 y otros puntos obligan a poder seleccionar
alternativas, pero no especifican al menos una alternativa que v1 deba ofrecer.
Una implementación solo con Poisson o con el modelo escalonado puede alegar que
no admite ninguna otra; otra puede considerar obligatorio un catálogo amplio.
Cada capacidad debe indicar si v1 exige una lista mínima, un mecanismo de
extensión futura o ambas cosas.

**Estado: resuelto.** V1 admite únicamente Poisson para los recuentos de
variantes aleatorias y el modelo escalonado simple de una unidad para mutaciones
de repeticiones en tándem. La tasa y el sesgo entre expansión y contracción
siguen siendo configurables. Otra distribución, otro modelo o un salto de varias
unidades se rechazan antes de simular; las alternativas se aplazan expresamente.

### RAD-021. Hardy-Weinberg y F en poliploides no determinan una distribución única

La especificación exige un modelo para ploidías mayores que dos, declara que no
distingue herencia polisómica o disómica y fija los extremos F=0 y F=1. Esas
condiciones no determinan por sí solas las probabilidades de todas las dosis
intermedias ni el comportamiento de F entre cero y uno. Elegir esa fórmula
cambia genotipos observables y no es solo una decisión interna. Debe aprobarse
el modelo funcional y ejemplos numéricos antes de F6.

**Estado: resuelto.** V1 separa Hardy-Weinberg del modelo automático F. HWE usa
la distribución binomial bialélica o multinomial multialélica para ploidías
1–10. F automático se limita a diploides, con fórmulas bialélicas y
multialélicas explícitas; se rechaza en ploidía 1 o mayor que 2. Una distribución
explícita prevalece localmente. F realizado solo se calcula en diploides y se
registra como «no aplicable» en las demás ploidías. Los modelos F poliploides se
aplazan.

### RAD-022. «Compatibilidad estructural» y transmisión de secuencias no reproductivas son vagas

La primera versión prohíbe reproducir heterocigotos para translocaciones, pero
permite tratar inversiones, duplicaciones y deleciones heterocigotas mediante
intervalos colineales. No se define el algoritmo funcional de compatibilidad ni
qué ocurre si no pueden formarse todos los bivalentes. Tampoco queda claro si un
*scaffold* con ploidía impar es reproductivo, cómo se transmite o si provoca el
rechazo general de REQ-CRU-24. Debe existir una tabla por tipo estructural,
ploidía y clase de secuencia antes de F5–F8.

**Estado: resuelto.** V1 solo reproduce secuencias de ploidía par 2–10, igual en
los parentales, y exige que todas sus copias tengan los mismos segmentos
estructurales, orden, orientación, número de copias y puntos de ruptura. SNP,
*indels* pequeños y repeticiones pueden diferir. Cualquier SV heterocigota o
diferencia estructural invalida al individuo o pareja completos; una SV fijada
idéntica se admite aunque proceda de historias distintas. Los bivalentes se
forman al azar y los *scaffolds* sin recombinación transmiten copias intactas.
Una SV de novo se conserva, pero hace incompatible a su portador para un uso
parental posterior mientras no cumpla el invariante. La compatibilidad parcial
queda futura.

### RAD-023. La relación entre un evento complejo, sus registros VCF y el complemento no tiene cardinalidad

Una translocación puede requerir varios BND y una variante compleja puede
producir varios registros. REQ-VER-25 habla de un identificador compartido con
«su registro VCF» y REQ-VER-44 intenta reunir alelos en un único registro cuando
sea posible. Falta definir relaciones uno-a-uno, uno-a-varios y varios-a-uno,
además de una normalización canónica suficiente para cumplir REQ-VER-46 entre
eventos. Debe cerrarse antes de F3 y F5.

**Estado: resuelto.** Un caso de verdad produce uno o varios registros VCF y,
solo si pierde detalle, exactamente una entrada complementaria. Todos los
registros que necesiten ese detalle enlazan la misma entrada y cada registro
enlazado pertenece a un solo caso. Una nueva adyacencia usa normalmente dos BND;
una translocación recíproca puede contener dos adyacencias y necesitar cuatro,
sin una cardinalidad universal por translocación.

### RAD-024. La calibración empírica es obligatoria como capacidad, pero no tiene resultado mínimo verificable

REQ-DAT-41 exige estimar longitudes, calidades, errores y «otros parámetros
útiles» e «intentar» separar variantes y errores de mapeo. No define qué campos
deben estimarse siempre, qué entradas mínimas se aceptan, cuándo la fiabilidad
insuficiente es advertencia o error ni qué significa que un perfil sea válido
para reutilizarlo. Es una capacidad funcional amplia que necesita un contrato
mínimo antes de F4.

**Estado: resuelto.** La calibración empírica automática se aplaza por completo
fuera de v1. No se aceptan FASTQ y referencia de calibración ni perfiles
empíricos reutilizables, no se infieren tecnología o parámetros y no se genera
el artefacto. Esos campos se rechazan antes de simular. Los perfiles
incorporados de Illumina y PacBio HiFi y sus sustituciones explícitas se
mantienen.

### RAD-025. Las validaciones estadísticas aún no pueden aceptar ni rechazar una ejecución

Expresiones como «aproximadamente Poisson», «compatible dentro de tolerancias»,
«tamaño manejable», «cuando pueda medirse» y «cuando sea posible» aparecen en
requisitos obligatorios. La sección 9 aplaza métricas y tolerancias, pero estas
determinan advertencias, errores y éxito observable. Debe aprobarse por cada
comprobación el ámbito, estadístico, tamaño mínimo, tolerancia y consecuencia
antes de implementar la fase correspondiente.

**Estado: resuelto.** V1 separa invariantes deterministas de magnitudes
aleatorias. Las incoherencias exactas detienen el ámbito afectado. Cobertura,
soporte alélico y tasas realizadas se informan sin tolerancias de aceptación ni
opción de elevarlas a error; una medida no calculable figura como «no evaluable»
con causa. Tras la decisión 20, solo N conserva su aviso específico no
invalidante; Ne no tiene objetivo en v1.

### RAD-026. Las sustituciones por individuo abarcan entidades que todavía no existen

REQ-CONF-17 afirma que todos los parámetros aleatorios pueden definirse por
individuo. En evolución y cruces, gran parte de los individuos se crea durante
la ejecución y sus identificadores pueden ser generados por el sistema. No se
aclara qué parámetros admiten realmente sustitución individual, cómo se nombra
un descendiente futuro ni si la regla pretende limitarse a fundadores,
individuos directos o selecciones declarativas. Debe acotarse en F1, F7 y F8.

**Estado: resuelto.** En v1, una sustitución por identificador individual solo
puede dirigirse a un individuo directo, fundador o parental expreso ya existente
durante la validación. Un identificador inexistente es un error. Los
descendientes futuros usan sus ámbitos poblacional, temporal y espacial, y
seleccionarlos después no modifica retroactivamente su creación. Las reglas para
descendientes futuros y la preasignación de sus identificadores se aplazan.

## 5. Hallazgos de severidad media

### RAD-027. «Evento», «experimento», «tanda» y «carril» necesitan un modelo cardinal explícito

El texto suele tratar un evento técnico como un experimento, pero también habla
de experimentos que contienen varias bibliotecas, tandas o carriles. Una tabla
de entidades y cardinalidades evitaría diferencias en directorios, subsemillas,
perfiles, grupos de lectura y manifiestos.

**Estado: resuelto.** En v1, el evento que genera FASTQ, BAM o CRAM es también
su único experimento de secuenciación y no recibe otra identidad. Contiene una
única tanda y exactamente una biblioteca y un grupo de lectura por individuo.
Todas las bibliotecas se procesan conjuntamente. Otra biblioteca, tanda, perfil,
tecnología o control requiere otro evento. Los carriles y las cardinalidades
mayores se aplazan.

### RAD-028. Los criterios de aceptación no trazan cada obligación

Los 770 requisitos se condensan en 109 criterios cualitativos. Muchos `ACE`
agrupan decenas de obligaciones y no indican los `REQ` que cubren. Conviene crear
una matriz `REQ → ACE → fase → prueba`, permitir varios casos por requisito y
detectar requisitos sin evidencia. Los `ACE` actuales pueden conservarse como
escenarios de alto nivel.

**Estado: resuelto.**
[trazabilidad-requisitos-0.1.md](./trazabilidad-requisitos-0.1.md) cubre los 788
REQ mediante 788 filas explícitas, sin huecos mecánicos. Dos auditorías
complementarias revisaron 412 y 376 filas y corrigieron asignaciones ACE, hitos
y evidencias falsas o genéricas. Cada fila conserva un TEST planificado; ninguna
se presenta como prueba ejecutada.

### RAD-029. El mismo contrato está repetido en demasiadas secciones normativas

Objetivo, requisitos, flujo funcional, tablas de entradas y salidas y criterios
de aceptación repiten reglas extensas. Ya hay pequeñas diferencias de alcance y
terminología. Debe declararse qué sección es normativa cuando dos redacciones
difieran y convertir el flujo y las tablas en resúmenes referenciados, no en una
segunda fuente de requisitos.

**Estado: resuelto.** La sección 4 de la especificación es la única fuente
normativa. Sus secciones 5–6 son ahora mapas conceptuales derivados que remiten
a rangos REQ y ACE sin repetir valores, cardinalidades, órdenes, errores ni
excepciones. Las secciones 7–12 de arquitectura conservan componentes,
fronteras, estados, dirección de dependencias y estrategia técnica, y remiten al
contrato funcional exacto cuando describen un efecto observable. Este informe
se declara histórico y no normativo. Los ACE conservan solo evidencia de
aceptación y los DTC la justificación técnica, ambos subordinados al REQ.

### RAD-030. [RESUELTO] La trazabilidad exhaustiva puede contradecir la prioridad de escala

El manifiesto relaciona cada individuo histórico; la verdad técnica conserva
información por lectura; los muestreos trazan segmentos hasta fundadores; y los
cruces guardan todos los individuos intermedios. Todo es coherente con ejemplos
pequeños, pero puede dominar memoria y almacenamiento en poblaciones grandes.
La pregunta 26 fija un único contrato obligatorio: evolución materializa los
estados completos inicial, final y solicitados y resume las generaciones
restantes; familias materializa todos sus individuos; cada lectura final
conserva detalle y los descartes solo agregados. La representación interna puede
optimizarse, pero la escala no reduce el contrato: un presupuesto insuficiente
se rechaza antes de simular.

### RAD-031. Los requisitos futuros no siempre tienen un destino recuperable

La sección 3 enumera ampliaciones futuras, pero no asigna identificadores ni una
fase o documento de destino a todas ellas. Conviene crear un registro de trabajo
futuro que preserve el motivo del aplazamiento y los requisitos afectados, en
especial ploidías impares reproductivas, cromosomas sexuales, multivalentes,
selección adicional, tecnologías aplazadas y protocolos no WGS.

**Estado: resuelto.**
[registro-ampliaciones-futuras-0.1.md](./registro-ampliaciones-futuras-0.1.md)
conserva DTF-001..036 con fuente, motivo/estado, destino post-v1 y dependencias;
distingue estas ampliaciones funcionales de DTP/DTPD. La auditoría clasificó
291 apariciones: 142 DTF, 16 DTP, 36 DTPD y 97 históricas o diagnósticas. Los
trece grupos antes huérfanos se incorporaron como DTF-024..036; quedan cero
exclusiones funcionales futuras huérfanas.

### RAD-032. La semántica de identificador FASTA e integridad de entradas queda demasiado tarde

Se exigen identificadores de secuencia únicos, pero no se define si el
identificador es el primer token de la cabecera o toda la descripción. Tampoco
se exige aún una huella de las entradas en el contrato de reproducibilidad,
aunque el manifiesto prevé datos de integridad. Esta decisión puede ser técnica,
pero debe cerrarse antes de aceptar FASTA, FASTQ, mapas o catálogos reales en
F1–F4.

### 5.1. Auditoría vigente de los 32 hallazgos

Esta clasificación contrasta cada ataque con los requisitos y criterios
vigentes. «Resuelto» exige una regla observable; una pregunta marcada como
hecha no basta por sí sola.

| RAD | Estado auditado | Evidencia vigente breve |
|---|---|---|
| RAD-001 | Resuelto | La cabecera de la especificación y DTC-001 confirman Rust para el núcleo y distinguen documento 0.1, producto v1 y decisiones técnicas pendientes |
| RAD-002 | Resuelto | La sección 8 declara cerradas las 27 decisiones funcionales y limita las elecciones restantes a planificación sin cambiar el contrato |
| RAD-003 | Resuelto | REQ-CONF-101 exige evento explícito para toda observación, interpreta cero eventos como cero salidas asociadas y REQ-CONF-102/103 permite verdad o VCF sin entidades técnicas |
| RAD-004 | Resuelto | REQ-WGS-04/11 fija la longitud haploide estructural como único denominador y separa la cobertura alineada por referencia; REQ-DAT-22 reutiliza el mismo FASTQ |
| RAD-005 | Resuelto | REQ-SAL-02 enumera el estado completo por vía; REQ-REC-11/12 limita el linaje evolutivo exportado y REQ-CRU-15 materializa todas las familias |
| RAD-006 | Resuelto | REQ-TE-01..14 están retirados y DTC-019 reclasifica copias e inserciones mediante operaciones genéricas |
| RAD-007 | Resuelto | REQ-REF-14/23, REQ-IND-58/59, REQ-DAT-72/73, REQ-WGS-14 y REQ-VER-47 fijan alfabeto y semántica completa de `N` |
| RAD-008 | Resuelto | REQ-IND-54..57 fija explícitos primero, orden de clases, espacio válido actualizado, ausencia de reintentos y fallo por insuficiencia |
| RAD-009 | Resuelto | REQ-CONF-116..123 fija matriz por familias, jerarquías internas, dominancia, conflictos, solapamientos, separación biológica/técnica y procedencia; REQ-VAL-51 y ACE-111 lo verifican |
| RAD-010 | Resuelto | REQ-REP-08..16/26 convierte N objetivo en intentos trazables, rechaza imposibilidad previa y retira Ne objetivo; Ne realizado conserva resultado o causa |
| RAD-011 | Resuelto | REQ-VER-19..26 separa VCF estándar interoperable, complemento único cuando hay pérdida y verdad estructurada canónica |
| RAD-012 | Resuelto | REQ-IND-53, REQ-MIG-09 y REQ-CRU-44 exigen una referencia estructural común por unidad reproductiva y rechazan mezclas antes de reproducir |
| RAD-013 | Resuelto | REQ-SAL-21..26 define estados, publicación atómica, conservación de conjuntos independientes y ausencia de resultados publicados tras fallo previo |
| RAD-014 | Resuelto | REQ-MUT-01 y REQ-MUT-13..15 definen oportunidades actuales por clase, copia y ámbito, incluidos los casos de `N` |
| RAD-015 | Resuelto | REQ-SEL-01..17 fija identidad, dosis, protección, aptitud válida y error previo ante máximo cero o locus ambiguo |
| RAD-016 | Resuelto | REQ-CONT-15/16 fija momento, denominador, unidades completas, redondeo conjunto, empates y fallo si no existe asignación válida |
| RAD-017 | Resuelto | REQ-WGS-11 separa objetivo, WGS básica, FASTQ final y cobertura alineada por referencia, con bases propias y contaminantes |
| RAD-018 | Resuelto | REQ-DAT-68..71 fija unicidad por evento, copia biológica, transferencia técnica y rechazo de colisiones sin renombrado incoherente |
| RAD-019 | Resuelto | REQ-CONF-111..115 separa reproducción biológica, técnica, equivalencia semántica e identidad binaria y REQ-MAN-69 registra sus condiciones |
| RAD-020 | Resuelto | REQ-CONF-25/61 y REQ-VAL-48 limitan v1 a Poisson y salto unitario, y rechazan alternativas antes de simular |
| RAD-021 | Resuelto | REQ-EVO-33/34 y REQ-EVO-46..52 fijan HWE en ploidía 1–10, F automático solo diploide y F realizado no aplicable fuera de diploides |
| RAD-022 | Resuelto | REQ-CRU-23..31 y REQ-REC-03..08 exigen ploidía par igual y estado estructural canónico idéntico, sin compatibilidad parcial |
| RAD-023 | Resuelto | REQ-VER-24..26 y REQ-VER-42 fijan caso→uno o varios registros→cero o una entrada, incluida la cardinalidad BND por adyacencia |
| RAD-024 | Resuelto | REQ-CONF-85..87, REQ-DAT-41..45 y REQ-SAL-10 retirados rechazan calibración y perfiles empíricos en v1 |
| RAD-025 | Resuelto | REQ-VAL-10..13 separa error determinista de informe estadístico sin tolerancias de aceptación y conserva «no evaluable» con causa |
| RAD-026 | Resuelto | REQ-CONF-17/18 limita sustituciones individuales a entidades preexistentes y prohíbe aplicación retroactiva a descendientes |
| RAD-027 | Resuelto | REQ-CONF-102, REQ-TEC-01..04 y REQ-RG-01..07 fijan evento/experimento, tanda, biblioteca, grupo de lectura y exclusión de carriles |
| RAD-028 | Resuelto | Las 788 filas REQ explícitas fueron auditadas en lotes de 412 y 376; no hay huecos, ACE/hitos inexistentes ni pruebas presentadas como superadas |
| RAD-029 | Resuelto | Especificación §§5–6 son mapas derivados; arquitectura §§7–12 conserva solo composición técnica con referencias REQ/ACE; el informe histórico se declara no normativo |
| RAD-030 | Resuelto | REQ-SAL-02/16, REQ-DAT-40, REQ-MAN-16/26/84 y REQ-NOF-05/06 fijan un contrato único y rechazo previo sin degradación |
| RAD-031 | Resuelto | DTF-001..036 y la auditoría de 291 apariciones dejan 0 exclusiones funcionales futuras huérfanas y separan DTP/DTPD técnicos |
| RAD-032 | Resuelto | REQ-REF-13 y REQ-VAL-45 usan el primer elemento de cabecera; REQ-MAN-76 exige sumas de los ficheros enlazados, incluidas copias de entradas |

No queda ninguna cuestión funcional abierta ni ningún hallazgo RAD pendiente.
RAD-028, RAD-029 y RAD-031 quedan resueltos con evidencia documental auditable.

## 6. Casos adversariales mínimos que deberían añadirse

Estos casos no eligen la solución. Fuerzan a que la decisión quede observable.

1. Un evento de 30× con un individuo de referencia A, mapeado contra A y una
   referencia B cuya longitud difiere un 20 %.
2. Una ejecución directa con individuos y cero eventos técnicos, sin indicar si
   se desactivan los alineamientos.
3. Un FASTA con `R`, `N` y una región solicitada para SNP, lectura y VCF.
4. Una región corta cuya muestra de Poisson solicita más SV no solapadas que
   posiciones válidas.
5. Un SNP aleatorio que cae dentro de una deleción explícita y otro que cae en
   su punto de ruptura.
6. Un parámetro definido a la vez globalmente, por región, individuo,
   subpoblación y fase.
7. Una población con individuos basados en referencias que difieren por una
   inversión y un locus seleccionado dentro de ella.
8. Una función de selección que asigna aptitud cero a todos los genotipos.
9. Una generación que completa sus N intentos con N descendientes; otra con
   éxitos parciales y fallos aleatorios trazados; una imposibilidad detectable
   antes de simular; un N realizado igual a 0 bajo cada política; y una
   configuración con Ne objetivo rechazada, mientras Ne realizado se registra
   o queda «no evaluable» con causa.
10. Un caso completamente representable sin complemento y una repetición
    compleja representada por un registro VCF más una entrada complementaria;
    todos sus campos `INFO` adicionales están declarados en la cabecera.
11. Una nueva adyacencia representada por dos BND más una entrada complementaria
    y una translocación recíproca representada por los cuatro BND de sus dos
    adyacencias más una sola entrada. El mismo evento contra dos referencias
    produce VCF y espacios de identificadores separados.
12. Contaminación combinada como porcentaje, cobertura y número de parejas, con
    redondeos que compiten por las últimas lecturas disponibles.
13. La misma lectura biológica copiada a dos receptoras y posteriormente
    afectada por reasignación técnica de índices.
14. Un fallo de validación en el último conjunto de una ejecución con otro
    conjunto ya generado correctamente.
15. Una generación evolutiva no guardada por completo que sea ancestro de un
    individuo muestreado varias generaciones después.
16. Una configuración de recuentos aleatorios con Poisson reproducible y otra
    que solicita binomial negativa y se rechaza identificando parámetro y valor.
17. Dos mutaciones de repetición en tándem, una expansión y una contracción de
    una unidad, y una tercera que solicita un salto de dos unidades y se rechaza.
18. Tres configuraciones que aportan, respectivamente, FASTQ y referencia de
    calibración, una solicitud de estimación y un perfil empírico reutilizable:
    todas se rechazan antes de simular identificando el campo. Una ejecución con
    un perfil incorporado y sustituciones permitidas completa sin generar ni
    enlazar ningún artefacto de perfil empírico.
19. Una cobertura aleatoria muy desviada se informa y permite completar; una
    lectura con procedencia imposible hace fallar su ámbito; una métrica sin
    unidades suficientes figura como «no evaluable» con causa; y N desviado
    produce su aviso no invalidante.
20. Para una misma tasa, comparar una, dos y cuatro copias haplotípicas; una
    región duplicada y otra eliminada; dos regiones configuradas solapadas; y un
    evento con destino o extremos condicionados. Comprobar además que dos copias
    actuales del mismo locus repetido aportan dos oportunidades.
21. Un SNP bialélico seleccionado válido en varias ploidías; un evento explícito
    solapante rechazado y uno aleatorio que nunca lo elige; máximo cero y
    aptitudes negativa o no finita; locus ausente o duplicado; y una
    recombinación que transmite el locus conservando su identidad única.
22. Inicializaciones HWE haploide, diploide y poliploide, bialélicas y
    multialélicas; F diploide igual a 0, intermedio e igual a 1; F rechazado en
    P=1 y P>2; distribución explícita poliploide que prevalece; y F realizado
    «no aplicable» fuera de diploides.
23. Parentales diploides y tetraploides estructuralmente idénticos; SNP,
    *indels* pequeños y repeticiones distintos que se admiten; una deleción,
    duplicación, inversión, inserción estructural y translocación heterocigotas
    que se rechazan; una SV fijada idéntica que se admite; ploidía impar y
    desigual; un *scaffold* sin recombinación; y una SV de novo cuyo portador se
    intenta usar después como parental.
24. Una generación con las ocho etapas aplicadas; otra con una etapa opcional
    omitida en su posición; una configuración con orden alternativo rechazada;
    varios eventos demográficos con orden interno explícito; y un historial que
    distingue etapas aplicadas y omitidas y conserva ambos órdenes.
25. Una evolución con una generación solicitada y otra no materializada; una
    familia cuyos individuos intermedios se conservan íntegramente; una lectura
    final reconstruible y una descartada representada solo en agregados; y un
    preflight cuyo presupuesto insuficiente rechaza sin reducir trazabilidad ni
    salidas. El manifiesto enumera los individuos materializados y de familias,
    pero no listas completas de generaciones evolutivas resumidas.
26. Un SNP explícito sobre `N` y un motivo repetido con `N` rechazados; una SV
    que contiene `N`; inserciones aleatoria y explícita; herencia y recombinación;
    una lectura con `N` y Phred 0; una ventana compuesta solo por `N`; cobertura
    que la cuenta; un caso VCF/complemento; y la comprobación de que ninguna
    operación resolvió ni creó `N` silenciosamente.

## 7. Estado del orden de resolución

Las acciones funcionales y documentales recomendadas por este informe están
incorporadas. Las fases de implementación deben convertir la evidencia
planificada de la matriz de trazabilidad en pruebas sin alterar los REQ.

## 8. Conclusión

No se recomienda reducir el alcance de manera automática a partir de este
informe. Las contradicciones funcionales identificadas ya tienen una decisión
observable para v1.

La especificación puede seguir siendo la base del proyecto. Para considerarla
cerrada deberá permitir que, ante cada caso adversarial anterior, exista una
única clase de resultados válidos o un error definido, observable y trazable.

Tras contrastar este informe con las decisiones posteriores, se ha corregido
RAD-001 en la especificación. También se han corregido la declaración de cierre
de RAD-002 y la ausencia de prioridad normativa de RAD-029. Las decisiones 1 a
27 de la sección siguiente se han aceptado e
incorporado a la especificación; RAD-002, RAD-003, RAD-004, RAD-005, RAD-012, RAD-013,
RAD-006, RAD-007, RAD-008, RAD-011, RAD-016, RAD-017, RAD-018,
RAD-009, RAD-010, RAD-014, RAD-015, RAD-019, RAD-020, RAD-021, RAD-022, RAD-023, RAD-024, RAD-025, RAD-026, RAD-027, RAD-030 y RAD-032 quedan resueltos.
RAD-028 y RAD-031 quedan resueltos mediante la matriz auditada y el registro de
expresiones diferidas. RAD-029 queda resuelto al sustituir los focos extensos
por mapas derivados, limitar la arquitectura a decisiones técnicas y declarar
este historial no normativo. No queda ninguna cuestión funcional abierta. Las
decisiones aceptadas se conservan a continuación para mantener su trazabilidad.

## Registro de decisiones valoradas

Las preguntas se ordenan por la primera fase que bloqueaban. Todas están
marcadas como **HECHA** y conservan la decisión aceptada incorporada a la
especificación. La numeración se mantiene estable para conservar trazabilidad.

### Antes de F0 y F1

1. **[HECHA] ¿Qué ocurre cuando se solicitan salidas técnicas sin declarar ningún
   evento?** (RAD-003; REQ-CONF-64, REQ-CONF-67, REQ-CONF-101,
   REQ-CONF-103, REQ-SAL-05; ACE-93).
   **Contexto:** los alineamientos están activados de forma predeterminada, pero
   sus parámetros pertenecen a eventos que pueden no existir. **Opciones:**
   exigir al menos un evento explícito; crear un evento implícito con valores
   definidos; o interpretar cero eventos como desactivación de toda salida
   técnica. También debe decidirse si puede existir un conjunto solo con VCF.
   **Consecuencias:** cambia qué configuraciones se aceptan y qué salidas se
   generan por omisión. **Propuesta:** exigir eventos explícitos para cualquier
   salida técnica y permitir verdad o VCF sin alineamientos cuando se soliciten.
   **Decisión aceptada:** toda salida de observación exige un evento explícito;
   cero eventos no crea salidas de evento. Un evento puede generar solo verdad o
   VCF sin tecnología, biblioteca, FASTQ, BAM, CRAM ni mapeador. Si solicita
   alineamientos, CRAM es el formato predeterminado dentro de ese evento.

2. **[HECHA] ¿Qué estados genotípicos completos debe materializar el fichero
   obligatorio de verdad?** (RAD-005; REQ-SAL-02, REQ-REP-19, REQ-CRU-15;
   ACE-14, ACE-24 y ACE-78).
   **Contexto:** evolución conserva por defecto ciertos hitos, mientras familias
   conserva todos los intermedios. **Opciones:** materializar todos los
   individuos históricos; materializar solo estados definidos por vía y guardar
   historial suficiente para reconstruir el resto; o hacer configurable la
   materialización sobre un mínimo obligatorio. **Consecuencias:** afecta tamaño,
   reconstrucción y consulta. **Propuesta:** fijar un mínimo distinto por vía y
   separar estado completo, historial de eventos y resúmenes.
   **Decisión aceptada:** en la vía directa se materializan todos los individuos
   creados; en evolución, la población inicial, la final y las generaciones o
   muestreos solicitados, con historial resumido del resto; en familias, todos
   los individuos de todos los eventos y poblaciones intermedias.

3. **[HECHA] ¿Cuál es la precedencia total entre ámbitos configurables?** (RAD-009;
   REQ-CONF-18, REQ-CONF-29 a REQ-CONF-34, REQ-CONF-62, REQ-CONF-73 y
   REQ-CONF-84; ACE-10 y ACE-17).
   **Contexto:** pueden coincidir valores globales, poblacionales, por
   subpoblación, fase, generación, individuo, secuencia, región y locus.
   **Opciones:** un orden total común; una matriz de precedencias por familia de
   parámetros; o prohibir combinaciones de ejes no enumeradas. **Consecuencias:**
   cambia el valor efectivo o el rechazo de una misma configuración.
   **Decisión aceptada:** cada familia enumera sus ejes permitidos. Dentro de un
   eje rigen las jerarquías de REQ-CONF-117, sin prioridad entre ejes distintos.
   Las reglas aplicables se resuelven por dominancia: una única máxima gana y
   varias máximas incomparables se rechazan, aunque coincidan en valor. También
   se rechazan duplicados, regiones solapadas y selectores no permitidos. Una
   regla combinada puede dominar reglas parciales. La configuración efectiva y
   el manifiesto conservan regla ganadora, selectores y procedencia; un evento
   técnico nunca sustituye parámetros biológicos.

4. **[HECHA] ¿Puede una unidad reproductiva mezclar individuos basados en referencias
   estructurales distintas?** (RAD-012; REQ-IND-02 y requisitos de selección,
   recombinación, migración y cruces; ACE-82, ACE-83, ACE-85, ACE-88 y
   ACE-97).
   **Contexto:** mapear contra varias referencias no define cómo reproducir
   estructuras diferentes. **Opciones:** prohibir la mezcla en v1; admitir solo
   clases de diferencias declaradas compatibles; o definir referencia de la
   descendencia, identidad de loci y mapas para cualquier mezcla. **Consecuencias:**
   afecta poblaciones válidas, genotipos y reproducción. **Propuesta:** prohibir
   en v1 la mezcla dentro de una unidad reproductiva y mantener independiente
   el mapeo contra varias referencias.
   **Decisión aceptada:** en v1, todos los individuos de una unidad reproductiva
   deben usar la misma referencia estructural. La regla se aplica a los
   conjuntos reproductores de las subpoblaciones, los cruces y los eventos
   demográficos o migratorios que los formen. Las unidades reproductivamente
   aisladas pueden usar referencias distintas y un individuo puede mapearse
   contra varias referencias. La reproducción entre referencias estructurales
   distintas se aplaza a futuras implementaciones.

5. **[HECHA] ¿Cuál es la unidad atómica ante un fallo?** (RAD-013; REQ-FMT-04,
   REQ-FMT-07, REQ-SAL-05 y REQ-MAN-73; ACE-94).
   **Contexto:** pueden fallar ejecución, evento, conjunto o artefacto después de
   que otros resultados estén validados. **Opciones:** invalidar siempre toda la
   ejecución; permitir eventos o conjuntos completos exportables; o definir la
   atomicidad por clase de fallo. **Consecuencias:** determina qué se publica,
   conserva y exporta. **Propuesta:** estados jerárquicos explícitos, publicación
   atómica por artefacto y conjunto, y éxito global solo si todas las salidas
   obligatorias están completas.
   **Decisión aceptada:** ejecución, evento, conjunto y artefacto tienen estados
   jerárquicos explícitos; los artefactos se publican atómicamente después de
   validarse y el éxito global requiere todas las salidas obligatorias completas.
   Los conjuntos completos e independientes se conservan y pueden exportarse
   declarando el fallo global.

6. **[HECHA] ¿Qué garantía exacta de reproducibilidad se ofrece?** (RAD-019;
   REQ-CONF-27 y requisitos de semillas técnicas; ACE-23, ACE-93 y ACE-101).
   **Contexto:** reproducir decisiones biológicas no implica FASTQ o archivos
   binariamente idénticos. **Opciones:** garantizar solo verdad biológica;
   garantizar además lecturas y decisiones técnicas dentro de una versión; o
   prometer identidad lógica o binaria de todas las salidas bajo condiciones
   concretas. **Consecuencias:** condiciona subsemillas, paralelismo, herramientas
   externas y pruebas. **Propuesta:** separar cuatro niveles —biológico,
   lecturas, equivalencia semántica e identidad binaria— y prometer cada uno
   solo bajo condiciones verificadas.
   **Decisión aceptada:** se distinguen reproducibilidad biológica, de lecturas y
   decisiones técnicas aleatorias, equivalencia semántica e identidad binaria.
   Los dos primeros niveles se garantizan bajo sus entradas, semillas, versiones
   y dependencias aplicables; la equivalencia semántica solo se declara bajo
   condiciones verificadas y no existe una garantía binaria general.

7. **[HECHA] ¿Qué cardinalidades relacionan evento, experimento, biblioteca, tanda y
   carril?** (RAD-027; REQ-CONF-101 a REQ-CONF-103, REQ-TEC-01 a REQ-TEC-04 y
   requisitos REQ-RG; ACE-93 y ACE-100).
   **Contexto:** el texto equipara a veces evento y experimento, pero admite
   varias bibliotecas, tandas o carriles. **Opciones:** hacer evento y
   experimento sinónimos con jerarquía interna fija; separarlos como entidades;
   o reducir v1 a una cardinalidad limitada. **Consecuencias:** cambia
   configuración, identificadores, grupos de lectura y manifiestos.
   **Propuesta:** usar una sola entidad técnica superior por evento y declarar
   cardinalidades exactas de sus componentes.
   **Decisión aceptada:** en v1, un evento que genera FASTQ, BAM o CRAM es
   también su único experimento de secuenciación, sin identidad separada. Tiene
   exactamente una tanda y una biblioteca por individuo; todas las bibliotecas
   se tratan conjuntamente y comparten tecnología, perfil y, cuando proceda,
   mapeador. Cada individuo tiene un grupo de lectura y, en modo realista, un
   FASTQ o pareja. Otra biblioteca, tanda, perfil, tecnología o control exige
   otro evento. Un evento solo de verdad o VCF no es un experimento y no tiene
   entidades técnicas. Los carriles y las cardinalidades mayores se aplazan.

8. **[HECHA] ¿A qué individuos pueden aplicarse sustituciones individuales?**
   (RAD-026; REQ-CONF-17, REQ-CONF-18 y requisitos dependientes; ACE-10).
   **Contexto:** muchos descendientes todavía no existen ni tienen identificador
   al leer la configuración. **Opciones:** limitar sustituciones a individuos
   declarados o fundadores; permitir reglas declarativas sobre descendientes;
   o exigir identificadores futuros preasignados. **Consecuencias:** cambia la
   expresividad y la resolución previa de la configuración. **Propuesta:**
   distinguir sustitución por identificador existente de reglas declarativas
   aplicables a descendientes futuros.
   **Decisión aceptada:** en v1, una sustitución por identificador individual
   solo puede aplicarse a individuos directos, fundadores o parentales expresos
   declarados y existentes antes de iniciar cualquier simulación. Un
   identificador inexistente durante la validación es un error. Los descendientes
   futuros usan los valores de sus ámbitos poblacional, temporal y espacial;
   seleccionarlos después no cambia retroactivamente su creación. Las reglas
   específicas sobre descendientes futuros y la preasignación de identificadores
   se aplazan.

9. **[HECHA] ¿Qué parte de una cabecera FASTA constituye el identificador de
   secuencia?** (RAD-032; REQ-REF-13; ACE-03).
   **Contexto:** la unicidad depende de si se usa el primer token o toda la
   descripción. **Opciones:** primer token separado por espacio; texto completo
   de la cabecera; o identificador elegido mediante una regla configurable.
   **Consecuencias:** cambia qué FASTA se acepta y cómo se enlazan coordenadas y
   formatos externos. **Propuesta:** usar el primer token como identificador y
   conservar por separado la descripción completa.
   **Decisión aceptada:** el primer elemento no vacío después de `>` y antes del
   primer separador es el identificador; el resto se conserva como descripción
   completa y no participa en la comprobación de unicidad.

### Antes de F2 a F5

10. **[HECHA] ¿Qué denominador fija la cobertura objetivo y qué coberturas realizadas
    se informan?** (RAD-004 y RAD-017; REQ-CONF-72 a REQ-CONF-74,
    REQ-WGS-04 y REQ-DAT-22; ACE-17 y ACE-20).
    **Contexto:** un FASTQ compartido no puede satisfacer el mismo objetivo
    frente a referencias de longitudes distintas, y los sesgos cambian el
    resultado por etapas. **Opciones:** longitud de la referencia estructural;
    suma de longitudes de haplotipos ponderadas por ploidía; referencia inicial;
    o FASTQ distintos por referencia de mapeo. **Consecuencias:** cambia el
    número de bases y las medidas posteriores. **Propuesta:** un denominador
    biológico único anterior al mapeo y medidas separadas antes de sesgos,
    después de artefactos y después del mapeo.
    **Decisión aceptada:** cada biblioteca usa la longitud haploide de la
    referencia estructural del individuo: cada cromosoma o *scaffold* se suma una
    vez, sin ploidía. Haplotipos y referencias de mapeo no cambian el denominador;
    las regiones usan la referencia estructural. Se registran por separado el
    objetivo, la cobertura WGS básica antes de sesgos y artefactos, la final del
    FASTQ con bases propias y contaminantes separadas y la alineada observada por
    referencia. El mismo FASTQ se mapea contra todas las referencias y la
    cobertura alineada nunca es otro objetivo.

11. **[HECHA] ¿Qué alfabeto de DNA admite el simulador?**
    (RAD-007; REQ-REF-14 y requisitos de generación, lectura y VCF; ACE-03).
    **Contexto original:** preservar, excluir o resolver símbolos ambiguos
    produce verdades distintas. **Decisión aceptada:** todas las secuencias
    suministradas o generadas solo pueden contener A, C, G, T y N, sin distinguir
    mayúsculas de minúsculas. Se rechazan los demás códigos IUPAC, los huecos y
    cualquier otro símbolo. Esta decisión cierra el alfabeto; la semántica de N
    en cada operación se resuelve en la decisión 27.

12. **[HECHA] ¿Pueden las inserciones de elementos transponibles crear referencias
    derivadas?** (RAD-006; REQ-REF-09 y REQ-TE-03; ACE-07 y ACE-80).
    **Contexto:** una referencia estructural solo admite cuatro clases de cambio,
    mientras que los TE pueden insertarse “al construir referencias”.
    **Opciones:** limitar los TE a la referencia inicial simulada; admitirlos
    como quinta clase de cambio de una referencia derivada; o crear una categoría
    de referencia distinta. **Consecuencias:** cambia qué referencias son
    válidas y cómo se registra su linaje. **Propuesta:** admitir inserciones TE
    solicitadas como cambio trazable de una referencia derivada y actualizar de
    forma conjunta la definición, REQ-REF-09, el flujo y ACE-07.
    **Decisión aceptada:** v1 elimina los elementos transponibles como clase
    funcional propia. No hay catálogo FASTA de TE, familias, mecanismos de copiar
    o cortar y pegar, reglas de destino, identidad ni metadatos específicos. Una
    región preexistente que se copie se trata como duplicación segmentaria
    genérica; una secuencia nueva puede ser una inserción genérica; un FASTA real
    se interpreta como bases sin anotación TE. Todo modelado específico se aplaza.

13. **[HECHA] ¿Cómo se resuelven colisiones y ausencia de sitios válidos durante la
    generación aleatoria?** (RAD-008; REQ-IND-54 a REQ-IND-57; ACE-10, ACE-87,
    ACE-103 y ACE-107).
    **Contexto:** eventos explícitos y aleatorios pueden competir o agotar una
    región. **Opciones:** remuestrear hasta un límite y fallar; aceptar menos
    eventos con advertencia; o aplicar una política configurable por clase.
    Debe fijarse también el orden entre clases. **Consecuencias:** cambia verdad,
    consumo aleatorio y número realizado.
    **Decisión aceptada:** primero se aplican los eventos explícitos en su orden
    declarado. Después se sortea el recuento aleatorio por clase y se procesan,
    en orden fijo, duplicaciones, deleciones, inversiones, translocaciones,
    repeticiones en tándem, *indels* y SNP. Cada evento se elige exclusivamente
    entre destinos válidos en el estado inmediatamente anterior y actualiza el
    espacio tras aceptarse; no hay reintentos sobre candidatos inválidos ni
    límite arbitrario. Si faltan destinos para completar el recuento sorteado,
    la generación falla sin reducirlo silenciosamente. Se registran los números
    sorteado, aceptado antes del fallo y finalmente activo.

14. **[HECHA] ¿Cómo se reparte la verdad entre VCF y complemento y qué cardinalidad
    une sus registros?** (RAD-011 y RAD-023; REQ-VER-04, REQ-VER-19 a
    REQ-VER-26, REQ-VER-42 y REQ-VER-44 a REQ-VER-46; ACE-33, ACE-99,
    ACE-102, ACE-103 y ACE-108).
    **Contexto:** algunos eventos no caben sin pérdida en un registro VCF y una
    translocación puede requerir varios BND. **Opciones:** exigir toda la
    información en VCF; definir un mínimo VCF y detalle exclusivo en complemento;
    o tratar VCF y complemento como una unidad inseparable. **Consecuencias:**
    afecta validez, reconstrucción e identificadores uno-a-varios o
    varios-a-uno. **Propuesta:** definir un núcleo VCF interoperable y un
    complemento enlazado con cardinalidades y normalización canónicas.
    **Decisión aceptada:** cada conjunto solicitado genera exactamente un VCF
    multimuesta referido a una sola referencia, sin separarlo por individuo,
    clase ni complejidad. Cumple estrictamente la versión oficial que se elija
    en planificación y usa la representación estándar más completa. Un caso
    produce uno o varios registros y, solo si pierde detalle,
    una única entrada complementaria; todos los registros que necesitan ese
    detalle incluyen campos `INFO` declarados con una descripción breve y el
    mismo identificador local. El VCF sigue siendo válido e interoperable por sí
    mismo, mientras la reconstrucción completa del caso usa el par enlazado. Una
    nueva adyacencia usa normalmente dos BND relacionados y una translocación
    recíproca puede necesitar cuatro para sus dos adyacencias. El fichero
    obligatorio de verdad sigue siendo la fuente canónica completa.

15. **[HECHA] ¿Cómo se convierten y redondean las unidades de contaminación?**
    (RAD-016; REQ-CONF-75, REQ-CONT-02 y REQ-CONT-03; ACE-42 y ACE-43).
    **Contexto:** proporción, cobertura equivalente y número de lecturas o pares
    deben formar una mezcla sustitutiva común. **Opciones:** convertir todo a
    recuentos antes del muestreo; resolver cada fuente en su unidad y ajustar al
    final; o prohibir unidades mezcladas en una muestra. **Consecuencias:** cambia
    el número exacto de lecturas propias y contaminantes. **Propuesta:** convertir
    conjuntamente a unidades enteras del tipo de biblioteca mediante una regla
    única de redondeo y desempate con semilla.
    **Decisión aceptada:** todas las fuentes se convierten conjuntamente, a
    partir del tamaño limpio realizado, a unidades indivisibles de la biblioteca:
    lectura individual para *single-end* o larga y pareja completa para
    *paired-end*. La asignación entera minimiza la desviación global, conserva el
    total entero más próximo, resuelve empates con la subsemilla y registra
    objetivos, conversión, resultado y desviaciones.

16. **[HECHA] ¿Qué identidad y reutilización se permiten para las lecturas?**
    (RAD-018; REQ-DAT-26, REQ-DAT-34 a REQ-DAT-36, REQ-CONT-01,
    REQ-CONT-10 y REQ-RG-04; ACE-100).
    **Contexto:** una lectura puede copiarse por contaminación biológica y
    transferirse por demultiplexado. **Opciones:** QNAME único en todo el evento;
    único solo por FASTQ o grupo; o reutilizable con un identificador técnico
    adicional. También debe limitarse la copia a varias receptoras. **Consecuencias:**
    afecta colisiones, BAM combinado y trazabilidad. **Propuesta:** identidad
    técnica única por instancia de lectura, con vínculo estable a la lectura de
    origen y política explícita de reutilización.
    **Decisión aceptada:** cada instancia de lectura o pareja tiene identidad
    técnica única dentro del evento. La contaminación biológica crea una nueva
    instancia vinculada a un origen estable sin retirar la donante; la
    reasignación técnica conserva la instancia y la retira de la donante. No se
    copia dos veces el mismo origen en una biblioteca receptora salvo mediante
    un mecanismo de duplicación configurado, y toda colisión es un error.

17. **[HECHA] ¿Qué significa que v1 admita “otras distribuciones” u “otros modelos”?**
    (RAD-020; REQ-CONF-25, REQ-CONF-61/62, REQ-DIR-04 y REQ-VAL-48;
    ACE-10, ACE-32, ACE-79 y ACE-109).
    **Contexto:** no se sabe si una alternativa concreta es obligatoria.
    **Opciones:** fijar un catálogo mínimo de v1; exigir solo un mecanismo de
    extensión futura; o eliminar la obligación de alternativas hasta una fase
    posterior. **Consecuencias:** cambia las configuraciones que toda
    implementación conforme debe aceptar. **Propuesta:** enumerar por capacidad
    el catálogo mínimo de v1 y separar explícitamente extensiones futuras.
    **Decisión aceptada:** el recuento realizado de variantes aleatorias usa
    exclusivamente Poisson. Los eventos explícitos y el modo controlado
    aplicable expresan resultados concretos, no distribuciones alternativas. Las
    mutaciones de repeticiones en tándem usan exclusivamente el modelo escalonado
    simple: cada evento añade o elimina una unidad completa y conserva
    configurable el sesgo entre expansión y contracción. Otra distribución,
    otro modelo o un salto de varias unidades se rechazan antes de simular e
    identifican parámetro y valor. Distribuciones adicionales y saltos mayores
    quedan para versiones futuras.

18. **[HECHA] ¿Cuál es el contrato mínimo de calibración empírica?** (RAD-024;
    REQ-CONF-85 a REQ-CONF-87, REQ-DAT-41 a REQ-DAT-45 y
    REQ-VAL-49; ACE-48).
    **Contexto:** “otros parámetros útiles” e “intentar separar” no permiten
    validar un perfil. **Opciones:** exigir un conjunto mínimo fijo de
    estimaciones; perfiles parciales con campos y fiabilidad declarados; o
    aplazar la calibración de v1. **Consecuencias:** cambia entradas válidas,
    resultados obligatorios y reutilización del perfil. **Propuesta:** admitir
    perfiles parciales solo con campos mínimos, métricas de fiabilidad y causas
    de parámetros no estimables.
    **Decisión aceptada:** la calibración empírica automática se aplaza por
    completo fuera de v1. Se rechazan antes de simular los FASTQ y referencias
    de calibración, las solicitudes de estimación o inferencia y los perfiles
    empíricos reutilizables, identificando el campo y la exclusión. V1 no infiere
    tecnología ni parámetros, no separa variantes y errores desde esos datos y
    no genera artefactos empíricos. Se mantienen los perfiles incorporados de
    Illumina y PacBio HiFi y sus sustituciones explícitas. Contrato de entrada,
    esquema, fiabilidad, separación y reutilización quedan para el futuro.

19. **[HECHA] ¿Qué reglas cuantitativas convierten una validación estadística en éxito,
    advertencia o error?** (RAD-025; REQ-VAL-01, REQ-VAL-04 a REQ-VAL-06,
    REQ-VAL-10, REQ-VAL-12 y REQ-VAL-13; ACE-64 y ACE-65).
    **Contexto:** expresiones como “aproximadamente” o “cuando sea posible” no
    deciden el resultado. **Opciones:** políticas fijas por comprobación;
    umbrales configurables con predeterminados; o informes sin capacidad de
    rechazo. **Consecuencias:** afecta el estado final de ejecuciones idénticas.
    **Propuesta:** definir antes de cada fase ámbito, estadístico, tamaño mínimo,
    tolerancia y consecuencia predeterminada, permitiendo elevar advertencias a
    error cuando ya esté autorizado.
    **Decisión aceptada:** v1 separa validaciones deterministas y resultados
    estadísticos. Toda incoherencia exacta detiene e impide completar su ámbito.
    Cobertura, soporte alélico, errores, contaminación, duplicados y demás
    magnitudes aleatorias se registran con ámbito, esperado, realizado, tamaño
    observado, desviaciones absoluta y relativa y explicación, sin tolerancia de
    aceptación ni opción de elevarlas a error. Una medida no calculable figura
    como «no evaluable» con causa. N conserva su aviso específico no invalidante
    y N=0 aplica su política propia; Ne no tiene objetivo en v1. Umbrales, intervalos,
    categorías y escalamiento quedan para versiones futuras.

### Antes de F6 a F8

20. **[HECHA] ¿Qué condiciones demuestran que N y Ne objetivo se intentaron aplicar?**
    (RAD-010 y la indefinición de REQ-REP-07 citada por RAD-002; REQ-REP-08,
    REQ-REP-14, REQ-REP-15 y REQ-REP-26; ACE-76 y ACE-77).
    **Contexto:** una advertencia no distingue imposibilidad aleatoria de un
    algoritmo que ignora el objetivo. **Opciones:** intentos y límites mínimos;
    invariantes del reparto reproductivo; o tolerancias de aproximación con
    diagnóstico de imposibilidad. **Consecuencias:** cambia cuándo una desviación
    es válida o defecto. **Propuesta:** fijar invariantes y evidencia de intento,
    mantener la desviación aleatoria como aviso y tratar como error no aplicar
    el objetivo o sus restricciones.
    **Decisión aceptada:** N objetivo es el número exacto de intentos
    reproductivos por generación y subpoblación; cada éxito crea un descendiente
    y los parentales pueden reutilizarse sin límite. Se procesan los N intentos o
    se registra inequívocamente la imposibilidad común de los restantes, con
    éxitos, fallos y causas. Una imposibilidad previa se rechaza; una aleatoria
    conserva descendientes, N realizado, desviación y aviso, y N=0 aplica su
    política. Ne objetivo absoluto o Ne/N se rechaza antes de simular. Ne
    realizado se calcula solo como Ne reproductivo cuando es evaluable, sin
    objetivo, desviación ni aviso; en otro caso conserva la causa. El control de
    Ne y los límites o distribuciones parentales quedan para el futuro.

21. **[HECHA] ¿Cuál es el espacio de oportunidades de cada tasa de mutación?**
    (RAD-014; REQ-MUT-01 y requisitos relacionados; ACE-79).
    **Contexto:** no se sabe si la tasa usa longitud haploide, todas las copias,
    longitud actual tras SV o solo posiciones válidas. **Opciones:** bases de
    referencia; bases actuales por haplotipo; o sitios válidos específicos de
    cada clase. **Consecuencias:** cambia la media esperada y su evolución tras
    deleciones o duplicaciones. **Propuesta:** definir por clase el conjunto de
    sitios válidos en el estado haplotípico inmediatamente anterior y registrar
    su tamaño.
    **Decisión aceptada:** una tasa evolutiva se aplica por oportunidad válida y
    generación en todas las copias haplotípicas descendientes actuales, justo
    antes de cada clase. Ploidía y duplicaciones aportan copias; las deleciones
    no. Los ámbitos efectivos no solapan oportunidades. La media es
    tasa×oportunidades y el recuento se sortea una vez; después se actualizan los
    destinos válidos sin volver a sortear. Cada clase usa candidatos propios y,
    si necesita varios extremos o fuente/destino, solo cuenta el ancla inicial.
    Cada locus repetido actual es una oportunidad. Se registran tasa, espacio,
    media y recuentos sorteado, aceptado y activo. Conforme a Q27, `N` se excluye
    cuando la clase necesita una base concreta, pero no por estar dentro de un
    intervalo puramente estructural.

22. **[HECHA] ¿Qué ocurre cuando la aptitud máxima es cero o el locus de selección no
    tiene una identidad inequívoca?** (RAD-015; REQ-SEL-13 y requisitos de locus;
    ACE-82).
    **Contexto:** no se puede normalizar por cero y una estructura puede
    eliminar o duplicar el locus. **Opciones:** error; supervivencia uniforme;
    extinción; o política configurable, junto con reglas de identidad de locus.
    **Consecuencias:** cambia supervivencia, N y trayectoria. **Propuesta:**
    rechazar antes de reproducir cuando la función o la identidad del locus no
    produzcan aptitudes definidas para todos los elegibles.
    **Decisión aceptada:** v1 usa un único SNP bialélico identificado por
    referencia estructural, secuencia, coordenada e ID, presente exactamente una
    vez por copia aplicable y protegido frente a todo cambio solapante. Para
    ploidías 1–10 se declaran valores para todas las dosis alternativas. Los
    eventos explícitos incompatibles se rechazan y los aleatorios excluyen el
    locus. Antes de supervivencia, cada elegible debe tener identidad, dosis,
    carácter y aptitud inequívocos; las aptitudes son finitas, no negativas y al
    menos una es positiva. Máximo cero o dato ausente, múltiple, ambiguo,
    negativo o no finito detiene con diagnóstico y no equivale a supervivencia
    uniforme ni extinción. La recombinación puede transmitir el locus sin perder
    su identidad. Otras clases de locus quedan para versiones futuras.

23. **[HECHA] ¿Qué distribución de dosis representa Hardy-Weinberg con F en
    poliploides?** (RAD-021; REQ-EVO-33, REQ-EVO-34 y REQ-EVO-46 a
    REQ-EVO-52; ACE-88, ACE-90 y ACE-104).
    **Contexto:** los extremos F=0 y F=1 no determinan los valores intermedios en
    ploidías mayores que dos. **Opciones:** aprobar una fórmula explícita;
    restringir F a diploides en v1; o exigir distribuciones de dosis explícitas
    para poliploides. **Consecuencias:** cambia los genotipos iniciales.
    **Propuesta:** no inferir una fórmula; restringir temporalmente F a los casos
    con modelo aprobado y usar distribución explícita en los demás.
    **Decisión aceptada:** Hardy-Weinberg y F automático son modelos distintos.
    HWE se admite en ploidías 1–10 y usa la distribución binomial por dosis en
    loci bialélicos o la multinomial por composición en multialélicos, sin
    sustituir los modos probabilístico y controlado ya definidos. F automático
    solo se admite en diploides: para dos alelos usa
    `P(AA)=p²+Fpq`, `P(AB)=2pq(1-F)` y `P(BB)=q²+Fpq`; para varios alelos usa
    `P(ii)=p_i²+F p_i(1-p_i)` y `P(ij)=2(1-F)p_i p_j` con `i<j`. F debe estar
    entre 0 y 1; F=0 equivale a HWE y F=1 produce homocigosis completa. F en
    P=1 o P>2 se rechaza antes de simular con ámbito y ploidía. Los poliploides
    usan HWE sin F o una distribución explícita, que prevalece localmente. F
    realizado solo se calcula en diploides y figura como «no aplicable» en las
    demás ploidías. Los modelos automáticos F poliploides, la herencia
    polisómica o disómica y sus estimadores se aplazan.

24. **[HECHA] ¿Qué combinaciones estructurales y ploidías pueden reproducirse?**
    (RAD-022; requisitos REQ-REC y REQ-CRU, incluido REQ-CRU-24; ACE-85 y
    ACE-98).
    **Contexto:** “compatibilidad estructural”, intervalos colineales y
    *scaffolds* no reproductivos no determinan todos los bivalentes.
    **Opciones:** una tabla cerrada por tipo estructural y ploidía; restringir v1
    a estructuras idénticas entre parentales; o definir un algoritmo completo
    de compatibilidad. **Consecuencias:** cambia qué cruces se aceptan y cómo se
    transmiten secuencias. **Propuesta:** aprobar una tabla cerrada de casos v1 y
    rechazar cualquier combinación no enumerada.
    **Decisión aceptada:** cruces, autofecundación y evolución reproductiva solo
    admiten ploidía par 2, 4, 6, 8 o 10 por cromosoma o *scaffold*, igual entre
    parentales. Todas las copias reproductivas deben compartir segmentos
    estructurales, orden, orientación, número de copias y puntos de ruptura; no
    se exige el mismo ID histórico si el estado canónico es inequívocamente
    idéntico. SNP, *indels* pequeños y repeticiones pueden diferir. Toda
    deleción, duplicación, inversión, inserción estructural o translocación debe
    estar fijada idénticamente en todas las copias y parentales; cualquier
    heterocigosis estructural invalida la opción reproductiva. Los bivalentes se
    forman al azar, cada parental aporta la mitad y un *scaffold* sin
    recombinación transmite copias intactas. Una incompatibilidad explícita se
    rechaza antes de simular; una realizada en evolución excluye la opción y
    causa los fallos de intento previstos por Q20. Una SV de novo se conserva,
    pero su portador no puede reproducirse después mientras sea incompatible.
    La compatibilidad parcial y sus consecuencias quedan futuras.

25. **[HECHA] ¿Qué órdenes alternativos de procesos generacionales son válidos?**
    (RAD-002 y RAD-009; REQ-CONF-44, REQ-CONF-45 y REQ-REP-24; ACE-84).
    **Contexto:** permitir reordenación sin enumerar dependencias admite
    trayectorias incompatibles. **Opciones:** solo el orden predeterminado en v1;
    catálogo cerrado de órdenes; o grafo de dependencias que valide cualquier
    orden topológico. **Consecuencias:** cambia la población resultante y las
    configuraciones aceptadas. **Propuesta:** empezar con el orden
    predeterminado y admitir únicamente alternativas enumeradas y justificadas.
    **Decisión aceptada:** v1 solo admite, en este orden, demografía y N objetivo;
    migración; selección de viabilidad; elección parental y reproducción;
    segregación y recombinación; mutación de novo; cálculo y registro de N
    realizado, Ne realizado y frecuencias; y muestreos y salidas. La
    configuración no puede reordenarlas y cualquier solicitud alternativa se
    rechaza antes de simular identificando el valor y la regla. Una etapa
    opcional desactivada se omite sin mover las demás. Varios eventos
    demográficos conservan un orden interno total explícito dentro de la primera
    etapa. Historial y manifiesto registran etapas aplicadas u omitidas y el
    orden efectivo. Otros órdenes y grafos de dependencias quedan futuros.

### Antes de F9

26. **[HECHA] ¿Qué trazabilidad detallada es obligatoria en escenarios grandes?**
    (RAD-030; REQ-CRU-15, REQ-REC-11 a REQ-REC-13, REQ-DAT-40,
    REQ-MAN-16 y REQ-MAN-26; ACE-14, ACE-86 y ACE-94).
    **Contexto:** conservar detalle por lectura, individuo y segmento puede
    dominar memoria y almacenamiento. **Opciones:** mantener siempre el detalle
    completo; permitir niveles configurables con un mínimo obligatorio; o usar
    agregados para entidades no observadas conservando detalle solo para
    muestreos. **Consecuencias:** afecta reconstrucción, auditoría y escala.
    **Propuesta:** definir un núcleo de verdad reconstruible obligatorio y
    niveles de materialización que nunca eliminen la trazabilidad exigida para
    resultados publicados.
    **Decisión aceptada:** v1 ofrece un único contrato obligatorio y no niveles
    configurables. En evolución materializa estado completo de la población
    inicial, final y generaciones o muestreos solicitados; las demás
    generaciones exportan resúmenes. Solo los individuos de esas
    materializaciones exportan recombinación y origen por segmento hasta copia
    fundadora, con mutaciones de novo en capas, aunque la ejecución conserva
    internamente lo necesario para reconstruirlo. En familias materializa el
    estado faseado y metadatos de todos los individuos de todos los eventos y
    poblaciones intermedias. Cada lectura o pareja del FASTQ final conserva
    origen y diferencia entre secuencia verdadera y observada; las descartadas
    o no asignadas solo conservan recuentos agregados por muestra y causa u
    origen conocido. El manifiesto enumera individuos materializados y de
    familias, y representa las demás generaciones evolutivas mediante
    resúmenes. Si el preflight determina que el presupuesto no permite el
    mínimo, esta trazabilidad y las salidas activadas, rechaza antes de simular;
    nunca degrada el detalle. La persona usuaria debe reducir solicitudes.
    Formatos y esquemas permanecen pendientes de planificación.

### Cuestión residual tras fijar el alfabeto

27. **[HECHA] ¿Cómo debe tratar cada operación las posiciones que contienen N?**
    (cuestión residual de RAD-007; REQ-REF-14, REQ-REF-23 y requisitos de
    generación, lectura y VCF; ACE-03).
    **Contexto:** N es el único símbolo admitido para una base desconocida, pero
    aceptar la secuencia no determina si una posición con N puede ser origen o
    destino de una variante, formar parte de un alelo REF/ALT, generar una
    lectura o contar en el espacio de oportunidades y la cobertura.
    **Opciones:** excluir N de todas las operaciones que requieran una base
    concreta; conservarla en lecturas y regiones sin variantes; o aprobar una
    política distinta y explícita para cada operación. **Consecuencias:** afecta
    mutaciones, variantes, lecturas, cobertura y VCF. **Propuesta:**
    definir una tabla normativa por operación y prohibir toda resolución
    aleatoria implícita de N. Esta cuestión se revisará junto con las demás
    decisiones pendientes.
    **Decisión aceptada:** `N` es una base desconocida que se conserva, ocupa
    coordenada y longitud y nunca se resuelve. No puede ser origen o destino SNP,
    locus seleccionado ni motivo repetido; sí puede permanecer dentro de una
    operación puramente estructural y de una inserción explícita o copiada. Las
    inserciones aleatorias contienen solo A/C/G/T. Se transmite sin cambios por
    herencia y recombinación. WGS puede muestrearla: aparece como `N` con Phred 0
    y se propaga al alineamiento; el modelo de error no la resuelve ni crea otras
    `N`. GC usa `(G+C)/(A+C+G+T)` y una ventana sin bases conocidas no recibe
    corrección. `N` cuenta en longitud y cobertura y se resume por separado. VCF
    no emite SNP con `REF=N`; otros casos usan `N` solo si el estándar elegido lo
    permite sin pérdida y, si no, el complemento de Q14. Verdad, validación y
    manifiesto conservan `N`, sus exclusiones y la ausencia de resolución.
