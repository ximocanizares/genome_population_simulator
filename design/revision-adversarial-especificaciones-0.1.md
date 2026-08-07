# Revisión adversarial de las especificaciones 0.1

- **Documento revisado:** `design/especificaciones-0.1.md`
- **Revisión de la especificación:** 0.1
- **Fecha de la revisión:** 7 de agosto de 2026
- **Commit de referencia:** `c9892a8`
- **SHA-256 del documento revisado:** `8d72785a27cc6903a6979aac46f7a0e1d966d0ab7e769549188eac9301d4d18b`
- **Estado:** informe para revisión; no modifica ni sustituye la especificación

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

La especificación contiene una base funcional extensa y útil, pero todavía no
puede considerarse un contrato funcional cerrado. La afirmación de la sección 8
de que no quedan decisiones funcionales abiertas es demasiado fuerte.

La comprobación mecánica ofrece resultados positivos:

- existen 730 definiciones `REQ` y no hay identificadores duplicados;
- toda mención con forma `REQ-…` corresponde a un requisito definido;
- los criterios `ACE-01` a `ACE-106` forman una secuencia completa, sin saltos
  ni duplicados;
- la separación entre verdad genotípica y evidencia técnica está expresada de
  forma consistente en la mayor parte del documento;
- la independencia de los conjuntos referidos a genomas distintos está
  reiterada y resulta inequívoca como intención general.

Los riesgos se concentran en la interacción entre requisitos. Los casos más
urgentes son la cobertura cuando un evento se mapea contra varias referencias,
el alcance del fichero obligatorio de verdad, la generación aleatoria de
variantes incompatibles, la precedencia entre configuraciones, el significado
operativo de N y Ne, y los límites entre VCF y complemento de verdad.

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
abierta. Sin embargo, requisitos como REQ-TE-12 dejan elegir entre advertencia y
error; REQ-REP-07 exige límites o distribuciones de descendencia sin fijar su
semántica mínima; REQ-CONF-44 permite órdenes alternativos sin enumerar todas
las dependencias; y diversos requisitos exigen «otras distribuciones» u «otros
modelos admitidos» sin definir un catálogo mínimo.

**Ataque:** dos implementaciones pueden cumplir literalmente el texto y, ante el
mismo escenario, una continuar con una advertencia, otra fallar y otra rechazar
una opción por no pertenecer a su catálogo.

**Decisión necesaria:** reabrir la sección 8 y trasladar a ella toda decisión
que cambie entradas aceptadas, salidas o errores observables. Mantener en la
sección 9 solo elecciones que no alteren el contrato funcional.

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

### RAD-007. No se define la biología de los símbolos IUPAC ambiguos

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
y cómo afectan a cobertura, variantes y calibración.

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

### RAD-009. No existe un orden total de precedencia entre todos los ámbitos configurables

**Evidencia:** se definen precedencias global, por secuencia, región e individuo,
y en otros lugares por población, subpoblación, fase, generación y locus. Se
afirma que se aplicarán las precedencias «correspondientes», pero no se define
un orden total ni cómo resolver ejes cruzados. Tampoco se fija la precedencia
entre el fichero general y el poblacional.

**Ataque:** un individuo de una subpoblación tiene una sustitución individual,
su locus está dentro de una región específica y la fase redefine el mismo
parámetro. No puede saberse cuál prevalece. Dos esquemas de configuración pueden
resolverlo de forma distinta sin contradecir una regla local.

**Decisión necesaria:** crear una matriz normativa por familia de parámetros con
origen permitido, precedencia, herencia temporal, regla de empate y error por
conflicto. Definir también qué significa que dos valores del mismo nivel sean
«compatibles».

**Bloquea:** F0 y F1.

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

### RAD-015. La selección tiene casos matemáticos y de identidad de locus sin resultado

REQ-SEL-13 normaliza por la aptitud máxima, pero no resuelve el caso en que todas
las aptitudes sean cero. Además, un locus ligado a una referencia estructural
puede estar eliminado, duplicado o ser ambiguo en individuos basados en otra
estructura. Deben definirse la respuesta a máximo cero y la semántica de
genotipo o dosis del locus en esos casos antes de F7.

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

### RAD-021. Hardy-Weinberg y F en poliploides no determinan una distribución única

La especificación exige un modelo para ploidías mayores que dos, declara que no
distingue herencia polisómica o disómica y fija los extremos F=0 y F=1. Esas
condiciones no determinan por sí solas las probabilidades de todas las dosis
intermedias ni el comportamiento de F entre cero y uno. Elegir esa fórmula
cambia genotipos observables y no es solo una decisión interna. Debe aprobarse
el modelo funcional y ejemplos numéricos antes de F6.

### RAD-022. «Compatibilidad estructural» y transmisión de secuencias no reproductivas son vagas

La primera versión prohíbe reproducir heterocigotos para translocaciones, pero
permite tratar inversiones, duplicaciones y deleciones heterocigotas mediante
intervalos colineales. No se define el algoritmo funcional de compatibilidad ni
qué ocurre si no pueden formarse todos los bivalentes. Tampoco queda claro si un
*scaffold* con ploidía impar es reproductivo, cómo se transmite o si provoca el
rechazo general de REQ-CRU-24. Debe existir una tabla por tipo estructural,
ploidía y clase de secuencia antes de F5–F8.

### RAD-023. La relación entre un evento complejo, sus registros VCF y el complemento no tiene cardinalidad

Una translocación puede requerir varios BND y una variante compleja puede
producir varios registros. REQ-VER-25 habla de un identificador compartido con
«su registro VCF» y REQ-VER-44 intenta reunir alelos en un único registro cuando
sea posible. Falta definir relaciones uno-a-uno, uno-a-varios y varios-a-uno,
además de una normalización canónica suficiente para cumplir REQ-VER-46 entre
eventos. Debe cerrarse antes de F3 y F5.

### RAD-024. La calibración empírica es obligatoria como capacidad, pero no tiene resultado mínimo verificable

REQ-DAT-41 exige estimar longitudes, calidades, errores y «otros parámetros
útiles» e «intentar» separar variantes y errores de mapeo. No define qué campos
deben estimarse siempre, qué entradas mínimas se aceptan, cuándo la fiabilidad
insuficiente es advertencia o error ni qué significa que un perfil sea válido
para reutilizarlo. Es una capacidad funcional amplia que necesita un contrato
mínimo antes de F4.

### RAD-025. Las validaciones estadísticas aún no pueden aceptar ni rechazar una ejecución

Expresiones como «aproximadamente Poisson», «compatible dentro de tolerancias»,
«tamaño manejable», «cuando pueda medirse» y «cuando sea posible» aparecen en
requisitos obligatorios. La sección 9 aplaza métricas y tolerancias, pero estas
determinan advertencias, errores y éxito observable. Debe aprobarse por cada
comprobación el ámbito, estadístico, tamaño mínimo, tolerancia y consecuencia
antes de implementar la fase correspondiente.

### RAD-026. Las sustituciones por individuo abarcan entidades que todavía no existen

REQ-CONF-17 afirma que todos los parámetros aleatorios pueden definirse por
individuo. En evolución y cruces, gran parte de los individuos se crea durante
la ejecución y sus identificadores pueden ser generados por el sistema. No se
aclara qué parámetros admiten realmente sustitución individual, cómo se nombra
un descendiente futuro ni si la regla pretende limitarse a fundadores,
individuos directos o selecciones declarativas. Debe acotarse en F1, F7 y F8.

## 5. Hallazgos de severidad media

### RAD-027. «Evento», «experimento», «tanda» y «carril» necesitan un modelo cardinal explícito

El texto suele tratar un evento técnico como un experimento, pero también habla
de experimentos que contienen varias bibliotecas, tandas o carriles. Una tabla
de entidades y cardinalidades evitaría diferencias en directorios, subsemillas,
perfiles, grupos de lectura y manifiestos.

### RAD-028. Los criterios de aceptación no trazan cada obligación

Los 730 requisitos se condensan en 106 criterios cualitativos. Muchos `ACE`
agrupan decenas de obligaciones y no indican los `REQ` que cubren. Conviene crear
una matriz `REQ → ACE → fase → prueba`, permitir varios casos por requisito y
detectar requisitos sin evidencia. Los `ACE` actuales pueden conservarse como
escenarios de alto nivel.

### RAD-029. El mismo contrato está repetido en demasiadas secciones normativas

Objetivo, requisitos, flujo funcional, tablas de entradas y salidas y criterios
de aceptación repiten reglas extensas. Ya hay pequeñas diferencias de alcance y
terminología. Debe declararse qué sección es normativa cuando dos redacciones
difieran y convertir el flujo y las tablas en resúmenes referenciados, no en una
segunda fuente de requisitos.

### RAD-030. La trazabilidad exhaustiva puede contradecir la prioridad de escala

El manifiesto relaciona cada individuo histórico; la verdad técnica conserva
información por lectura; los muestreos trazan segmentos hasta fundadores; y los
cruces guardan todos los individuos intermedios. Todo es coherente con ejemplos
pequeños, pero puede dominar memoria y almacenamiento en poblaciones grandes.
Antes de F9 debe distinguirse qué información debe estar indexada y consultable,
qué puede transmitirse como flujo, y qué agregados son suficientes sin perder
la verdad obligatoria.

### RAD-031. Los requisitos futuros no siempre tienen un destino recuperable

La sección 3 enumera ampliaciones futuras, pero no asigna identificadores ni una
fase o documento de destino a todas ellas. Conviene crear un registro de trabajo
futuro que preserve el motivo del aplazamiento y los requisitos afectados, en
especial ploidías impares reproductivas, cromosomas sexuales, multivalentes,
selección adicional, tecnologías aplazadas y protocolos no WGS.

### RAD-032. La semántica de identificador FASTA e integridad de entradas queda demasiado tarde

Se exigen identificadores de secuencia únicos, pero no se define si el
identificador es el primer token de la cabecera o toda la descripción. Tampoco
se exige aún una huella de las entradas en el contrato de reproducibilidad,
aunque el manifiesto prevé datos de integridad. Esta decisión puede ser técnica,
pero debe cerrarse antes de aceptar FASTA, FASTQ, mapas o catálogos reales en
F1–F4.

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
9. Un Ne objetivo incompatible con los límites de descendencia, y otro posible
   pero difícil de aproximar.
10. Una repetición compleja cuya pureza e interrupciones no caben completamente
    en un único registro VCF.
11. Una translocación que necesita varios BND y un único caso complementario.
12. Contaminación combinada como porcentaje, cobertura y número de parejas, con
    redondeos que compiten por las últimas lecturas disponibles.
13. La misma lectura biológica copiada a dos receptoras y posteriormente
    afectada por reasignación técnica de índices.
14. Un fallo de validación en el último conjunto de una ejecución con otro
    conjunto ya generado correctamente.
15. Una generación evolutiva no guardada por completo que sea ancestro de un
    individuo muestreado varias generaciones después.

## 7. Orden recomendado de resolución

Antes de cerrar F0:

1. RAD-001, RAD-002 y RAD-029: fijar las fuentes normativas y el estado real de
   las decisiones.
2. RAD-003, RAD-004, RAD-005, RAD-009, RAD-012 y RAD-013: cerrar el modelo de
   ejecución, referencias, verdad, configuración y fallos.
3. RAD-019 y RAD-027: fijar el contrato de reproducibilidad y las entidades
   técnicas.

Antes de F2–F5:

1. RAD-006, RAD-007 y RAD-008: cerrar referencias, ambigüedad y colisiones.
2. RAD-011, RAD-018 y RAD-023: cerrar la relación entre verdad, lecturas, VCF y
   complemento.
3. RAD-016, RAD-017, RAD-024 y RAD-025: definir magnitudes y aceptación del
   realismo técnico.

Antes de F6–F8:

1. RAD-010, RAD-014, RAD-015, RAD-021 y RAD-022: cerrar el significado
   cuantitativo de población, selección, recombinación y poliploidía.
2. RAD-020 y RAD-026: delimitar catálogos y ámbitos configurables de v1.

## 8. Conclusión

No se recomienda reducir el alcance de manera automática a partir de este
informe. Sí se recomienda corregir primero las contradicciones que impiden saber
qué resultado es válido y convertir las decisiones funcionales restantes en
preguntas explícitas.

La especificación puede seguir siendo la base del proyecto. Para considerarla
cerrada deberá permitir que, ante cada caso adversarial anterior, exista una
única clase de resultados válidos o un error definido, observable y trazable.

Tras contrastar este informe con las decisiones posteriores, se ha corregido
RAD-001 en la especificación. También se han corregido la declaración de cierre
de RAD-002 y la ausencia de prioridad normativa de RAD-029. Las decisiones 1,
2, 5, 6, 9, 11, 15 y 16 de la sección siguiente se han aceptado e incorporado a
la especificación; RAD-003, RAD-005, RAD-013, RAD-016, RAD-018, RAD-019 y la
semántica funcional del identificador FASTA de RAD-032 quedan resueltos.
RAD-007 queda resuelto respecto al alfabeto permitido, pero mantiene pendiente
el tratamiento operativo de N. RAD-028 y RAD-031 continúan como mejoras
documentales: crear una matriz de trazabilidad y un registro recuperable de
ampliaciones futuras no exige elegir ahora un comportamiento del simulador. Los
demás hallazgos que cambian resultados, entradas aceptadas o errores se plantean
a continuación sin presuponer una respuesta.

## Decisiones pendientes de valoración

Las preguntas se ordenan por la primera fase que bloquean. Las marcadas como
**HECHA** contienen una decisión aceptada e incorporada a la especificación; las
demás conservan una **propuesta** que todavía no es una decisión vigente. La
numeración se mantiene estable para conservar la trazabilidad de la revisión.

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

3. **¿Cuál es la precedencia total entre ámbitos configurables?** (RAD-009;
   REQ-CONF-18, REQ-CONF-29 a REQ-CONF-34, REQ-CONF-62, REQ-CONF-73 y
   REQ-CONF-84; ACE-10 y ACE-17).
   **Contexto:** pueden coincidir valores globales, poblacionales, por
   subpoblación, fase, generación, individuo, secuencia, región y locus.
   **Opciones:** un orden total común; una matriz de precedencias por familia de
   parámetros; o prohibir combinaciones de ejes no enumeradas. **Consecuencias:**
   cambia el valor efectivo o el rechazo de una misma configuración.
   **Propuesta:** aprobar una matriz por familia, con procedencia permitida,
   herencia temporal, empate y conflicto explícitos.

4. **¿Puede una unidad reproductiva mezclar individuos basados en referencias
   estructurales distintas?** (RAD-012; REQ-IND-02 y requisitos de selección,
   recombinación, migración y cruces; ACE-82, ACE-85, ACE-88 y ACE-97).
   **Contexto:** mapear contra varias referencias no define cómo reproducir
   estructuras diferentes. **Opciones:** prohibir la mezcla en v1; admitir solo
   clases de diferencias declaradas compatibles; o definir referencia de la
   descendencia, identidad de loci y mapas para cualquier mezcla. **Consecuencias:**
   afecta poblaciones válidas, genotipos y reproducción. **Propuesta:** prohibir
   en v1 la mezcla dentro de una unidad reproductiva y mantener independiente
   el mapeo contra varias referencias.

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

7. **¿Qué cardinalidades relacionan evento, experimento, biblioteca, tanda y
   carril?** (RAD-027; REQ-CONF-101 a REQ-CONF-103, REQ-TEC-01 a REQ-TEC-04 y
   requisitos REQ-RG; ACE-93 y ACE-100).
   **Contexto:** el texto equipara a veces evento y experimento, pero admite
   varias bibliotecas, tandas o carriles. **Opciones:** hacer evento y
   experimento sinónimos con jerarquía interna fija; separarlos como entidades;
   o reducir v1 a una cardinalidad limitada. **Consecuencias:** cambia
   configuración, identificadores, grupos de lectura y manifiestos.
   **Propuesta:** usar una sola entidad técnica superior por evento y declarar
   cardinalidades exactas de sus componentes.

8. **¿A qué individuos pueden aplicarse sustituciones individuales?**
   (RAD-026; REQ-CONF-17, REQ-CONF-18 y requisitos dependientes; ACE-10).
   **Contexto:** muchos descendientes todavía no existen ni tienen identificador
   al leer la configuración. **Opciones:** limitar sustituciones a individuos
   declarados o fundadores; permitir reglas declarativas sobre descendientes;
   o exigir identificadores futuros preasignados. **Consecuencias:** cambia la
   expresividad y la resolución previa de la configuración. **Propuesta:**
   distinguir sustitución por identificador existente de reglas declarativas
   aplicables a descendientes futuros.

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

10. **¿Qué denominador fija la cobertura objetivo y qué coberturas realizadas
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

11. **[HECHA] ¿Qué alfabeto de DNA admite el simulador?**
    (RAD-007; REQ-REF-14 y requisitos de generación, lectura y VCF; ACE-03).
    **Contexto original:** preservar, excluir o resolver símbolos ambiguos
    produce verdades distintas. **Decisión aceptada:** todas las secuencias
    suministradas o generadas solo pueden contener A, C, G, T y N, sin distinguir
    mayúsculas de minúsculas. Se rechazan los demás códigos IUPAC, los huecos y
    cualquier otro símbolo. Esta decisión cierra el alfabeto; la semántica de N
    en cada operación se separa como decisión 27.

12. **¿Pueden las inserciones de elementos transponibles crear referencias
    derivadas?** (RAD-006; REQ-REF-09 y REQ-TE-03; ACE-07 y ACE-80).
    **Contexto:** una referencia estructural solo admite cuatro clases de cambio,
    mientras que los TE pueden insertarse “al construir referencias”.
    **Opciones:** limitar los TE a la referencia inicial simulada; admitirlos
    como quinta clase de cambio de una referencia derivada; o crear una categoría
    de referencia distinta. **Consecuencias:** cambia qué referencias son
    válidas y cómo se registra su linaje. **Propuesta:** admitir inserciones TE
    solicitadas como cambio trazable de una referencia derivada y actualizar de
    forma conjunta la definición, REQ-REF-09, el flujo y ACE-07.

13. **¿Cómo se resuelven colisiones y ausencia de sitios válidos durante la
    generación aleatoria?** (RAD-008 y la política de REQ-TE-12 citada por
    RAD-002; requisitos estructurales; ACE-10, ACE-81, ACE-87 y ACE-103).
    **Contexto:** eventos explícitos y aleatorios pueden competir o agotar una
    región. **Opciones:** remuestrear hasta un límite y fallar; aceptar menos
    eventos con advertencia; o aplicar una política configurable por clase.
    Debe fijarse también el orden entre clases. **Consecuencias:** cambia verdad,
    consumo aleatorio y número realizado. **Propuesta:** aplicar primero eventos
    explícitos, definir un orden de clases y un límite reproducible de intentos,
    con fallo predeterminado cuando no pueda cumplirse el contrato solicitado.

14. **¿Cómo se reparte la verdad entre VCF y complemento y qué cardinalidad
    une sus registros?** (RAD-011 y RAD-023; REQ-VER-22, REQ-VER-24 a
    REQ-VER-26 y REQ-VER-44 a REQ-VER-46; ACE-80, ACE-99, ACE-102 y ACE-103).
    **Contexto:** algunos eventos no caben sin pérdida en un registro VCF y una
    translocación puede requerir varios BND. **Opciones:** exigir toda la
    información en VCF; definir un mínimo VCF y detalle exclusivo en complemento;
    o tratar VCF y complemento como una unidad inseparable. **Consecuencias:**
    afecta validez, reconstrucción e identificadores uno-a-varios o
    varios-a-uno. **Propuesta:** definir un núcleo VCF interoperable y un
    complemento enlazado con cardinalidades y normalización canónicas.

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

17. **¿Qué significa que v1 admita “otras distribuciones” u “otros modelos”?**
    (RAD-020; REQ-CONF-25, REQ-DIR-04, REQ-CONF-61 y requisitos análogos;
    ACE-10 y ACE-79).
    **Contexto:** no se sabe si una alternativa concreta es obligatoria.
    **Opciones:** fijar un catálogo mínimo de v1; exigir solo un mecanismo de
    extensión futura; o eliminar la obligación de alternativas hasta una fase
    posterior. **Consecuencias:** cambia las configuraciones que toda
    implementación conforme debe aceptar. **Propuesta:** enumerar por capacidad
    el catálogo mínimo de v1 y separar explícitamente extensiones futuras.

18. **¿Cuál es el contrato mínimo de calibración empírica?** (RAD-024;
    REQ-DAT-41 a REQ-DAT-45; ACE-48).
    **Contexto:** “otros parámetros útiles” e “intentar separar” no permiten
    validar un perfil. **Opciones:** exigir un conjunto mínimo fijo de
    estimaciones; perfiles parciales con campos y fiabilidad declarados; o
    aplazar la calibración de v1. **Consecuencias:** cambia entradas válidas,
    resultados obligatorios y reutilización del perfil. **Propuesta:** admitir
    perfiles parciales solo con campos mínimos, métricas de fiabilidad y causas
    de parámetros no estimables.

19. **¿Qué reglas cuantitativas convierten una validación estadística en éxito,
    advertencia o error?** (RAD-025; REQ-VAL-01, REQ-VAL-04 a REQ-VAL-06,
    REQ-VAL-10, REQ-VAL-12 y REQ-VAL-13; ACE-64 y ACE-65).
    **Contexto:** expresiones como “aproximadamente” o “cuando sea posible” no
    deciden el resultado. **Opciones:** políticas fijas por comprobación;
    umbrales configurables con predeterminados; o informes sin capacidad de
    rechazo. **Consecuencias:** afecta el estado final de ejecuciones idénticas.
    **Propuesta:** definir antes de cada fase ámbito, estadístico, tamaño mínimo,
    tolerancia y consecuencia predeterminada, permitiendo elevar advertencias a
    error cuando ya esté autorizado.

### Antes de F6 a F8

20. **¿Qué condiciones demuestran que N y Ne objetivo se intentaron aplicar?**
    (RAD-010 y la indefinición de REQ-REP-07 citada por RAD-002; REQ-REP-08,
    REQ-REP-14, REQ-REP-15 y REQ-REP-26; ACE-76 y ACE-77).
    **Contexto:** una advertencia no distingue imposibilidad aleatoria de un
    algoritmo que ignora el objetivo. **Opciones:** intentos y límites mínimos;
    invariantes del reparto reproductivo; o tolerancias de aproximación con
    diagnóstico de imposibilidad. **Consecuencias:** cambia cuándo una desviación
    es válida o defecto. **Propuesta:** fijar invariantes y evidencia de intento,
    mantener la desviación aleatoria como aviso y tratar como error no aplicar
    el objetivo o sus restricciones.

21. **¿Cuál es el espacio de oportunidades de cada tasa de mutación?**
    (RAD-014; REQ-MUT-01 y requisitos relacionados; ACE-79).
    **Contexto:** no se sabe si la tasa usa longitud haploide, todas las copias,
    longitud actual tras SV o solo posiciones válidas. **Opciones:** bases de
    referencia; bases actuales por haplotipo; o sitios válidos específicos de
    cada clase. **Consecuencias:** cambia la media esperada y su evolución tras
    deleciones o duplicaciones. **Propuesta:** definir por clase el conjunto de
    sitios válidos en el estado haplotípico inmediatamente anterior y registrar
    su tamaño.

22. **¿Qué ocurre cuando la aptitud máxima es cero o el locus de selección no
    tiene una identidad inequívoca?** (RAD-015; REQ-SEL-13 y requisitos de locus;
    ACE-82).
    **Contexto:** no se puede normalizar por cero y una estructura puede
    eliminar o duplicar el locus. **Opciones:** error; supervivencia uniforme;
    extinción; o política configurable, junto con reglas de identidad de locus.
    **Consecuencias:** cambia supervivencia, N y trayectoria. **Propuesta:**
    rechazar antes de reproducir cuando la función o la identidad del locus no
    produzcan aptitudes definidas para todos los elegibles.

23. **¿Qué distribución de dosis representa Hardy-Weinberg con F en
    poliploides?** (RAD-021; REQ-EVO-33, REQ-EVO-34 y REQ-EVO-46 a
    REQ-EVO-48; ACE-88 y ACE-90).
    **Contexto:** los extremos F=0 y F=1 no determinan los valores intermedios en
    ploidías mayores que dos. **Opciones:** aprobar una fórmula explícita;
    restringir F a diploides en v1; o exigir distribuciones de dosis explícitas
    para poliploides. **Consecuencias:** cambia los genotipos iniciales.
    **Propuesta:** no inferir una fórmula; restringir temporalmente F a los casos
    con modelo aprobado y usar distribución explícita en los demás.

24. **¿Qué combinaciones estructurales y ploidías pueden reproducirse?**
    (RAD-022; requisitos REQ-REC y REQ-CRU, incluido REQ-CRU-24; ACE-85 y
    ACE-98).
    **Contexto:** “compatibilidad estructural”, intervalos colineales y
    *scaffolds* no reproductivos no determinan todos los bivalentes.
    **Opciones:** una tabla cerrada por tipo estructural y ploidía; restringir v1
    a estructuras idénticas entre parentales; o definir un algoritmo completo
    de compatibilidad. **Consecuencias:** cambia qué cruces se aceptan y cómo se
    transmiten secuencias. **Propuesta:** aprobar una tabla cerrada de casos v1 y
    rechazar cualquier combinación no enumerada.

25. **¿Qué órdenes alternativos de procesos generacionales son válidos?**
    (RAD-002 y RAD-009; REQ-CONF-44, REQ-CONF-45 y REQ-REP-24; ACE-84).
    **Contexto:** permitir reordenación sin enumerar dependencias admite
    trayectorias incompatibles. **Opciones:** solo el orden predeterminado en v1;
    catálogo cerrado de órdenes; o grafo de dependencias que valide cualquier
    orden topológico. **Consecuencias:** cambia la población resultante y las
    configuraciones aceptadas. **Propuesta:** empezar con el orden
    predeterminado y admitir únicamente alternativas enumeradas y justificadas.

### Antes de F9

26. **¿Qué trazabilidad detallada es obligatoria en escenarios grandes?**
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

### Cuestión residual tras fijar el alfabeto

27. **¿Cómo debe tratar cada operación las posiciones que contienen N?**
    (cuestión residual de RAD-007; REQ-REF-14, REQ-REF-23 y requisitos de
    generación, lectura y VCF; ACE-03).
    **Contexto:** N es el único símbolo admitido para una base desconocida, pero
    aceptar la secuencia no determina si una posición con N puede ser origen o
    destino de una variante, formar parte de un alelo REF/ALT, generar una
    lectura o contar en el espacio de oportunidades y la cobertura.
    **Opciones:** excluir N de todas las operaciones que requieran una base
    concreta; conservarla en lecturas y regiones sin variantes; o aprobar una
    política distinta y explícita para cada operación. **Consecuencias:** afecta
    mutaciones, variantes, lecturas, cobertura, calibración y VCF. **Propuesta:**
    definir una tabla normativa por operación y prohibir toda resolución
    aleatoria implícita de N. Esta cuestión se revisará junto con las demás
    decisiones pendientes.
