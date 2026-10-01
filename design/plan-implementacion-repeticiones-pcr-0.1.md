# Plan de implementación — repeticiones en tándem y tartamudeo, con y sin PCR 0.1

## 1. Estado, autoridad y alcance

- **Fecha:** 19 de septiembre de 2026. **Reescrito el 20 de septiembre de
  2026** contra las decisiones cerradas por el responsable y contra los
  hallazgos de la revisión adversarial de ese mismo día.
- **Estado:** **ejecutable en las entregas 1 a 4; la quinta, parcialmente
  bloqueada** por una decisión que el responsable no llegó a fijar —el
  pendiente **P-4**, cómo el generador de deslizamiento biológico produce una
  distribución de frecuencias y no una sola secuencia de alelos—. Las
  veintiocho decisiones que los dos documentos
  de origen dejaban abiertas **están cerradas**: el responsable las resolvió
  una a una el 20 de septiembre y el registro literal de cada una vive en la
  sección «Decisiones confirmadas por el responsable». Cinco de ellas
  **corrigieron** la recomendación que este plan había adoptado. Los pasos de
  las secciones 5 a 9 están escritos contra lo decidido, no contra lo
  recomendado.
- **Cómo nació y por qué la tabla de supuestos sigue aquí.** La primera versión
  de este plan fue una excepción a la regla que prohíbe a un plan decidir
  diseño: adoptó la recomendación de cada ficha y la marcó como supuesto, para
  que el responsable pudiera confirmarla o cambiarla de un vistazo. Funcionó, y
  por eso la tabla de la sección 2 se conserva: es el **registro de qué se
  propuso**, y la columna de estado dice cuál se aceptó y cuál se corrigió.
  Leerla ya no sirve para decidir —eso está hecho— sino para entender por qué
  un paso dice lo que dice.
- **Por qué existió la excepción.** El plan de la fase de variación compleja
  (`plan-implementacion-f6-0.1.md`) se paró ante cada decisión abierta y dejó
  21 de sus 46 pasos bloqueados. Un plan hueco no sirve para trabajar.
- **Revisión adversarial del 20 de septiembre de 2026.** Declaró el plan no
  integrable: los treinta y ocho pasos seguían escritos contra lo recomendado y
  no contra lo decidido; siete piezas decididas no tenían ningún paso; la
  sección de correcciones a otros documentos proponía consagrar en la
  especificación algo que acababa de rechazarse; cuatro hechos afirmados eran
  falsos; y diez pruebas no discriminaban. Todo ello está aplicado en esta
  versión, y los cuatro hechos falsos se enumeran en la sección 1.1.
- **Fuentes normativas y de diseño:**
  - `design/expediente-decisiones-repeticiones-tandem-0.1.md` — corte mínimo
    para tener microsatélites que varíen entre individuos, dieciocho fichas de
    decisión y el borrador del requisito del tartamudeo.
  - `design/diseno-bibliotecas-pcr-0.1.md` — muestreo por familias de
    amplificación, duplicados, sesgo por composición, tartamudeo clonal, cota
    de recursos, cambio de contrato y once decisiones.
  - `design/calibracion-tartamudeo-0.1.md` — valores medidos con y sin PCR,
    **declarados humanos**, con sus fuentes y su situación de licencia.
  - `design/especificaciones-0.1.md`, `design/decisiones-tecnicas-0.1.md`,
    `design/hoja-de-ruta-0.1.md` (§7, reglas de avance),
    `design/plan-implementacion-f6-0.1.md` §8 (hito de repeticiones),
    `design/plan-implementacion-f4-0.1.md` (fase de realismo técnico, cerrada).
  - Requisitos existentes que el trabajo toca: REQ-CONF-55/57/59/60/63,
    REQ-CONF-88/89/94/95, REQ-IND-31/32/33/34, REQ-REF-19/23, REQ-TEC-07/08,
    REQ-DAT-17/46/47/52/53/56, REQ-MAN-29/45/48/60, REQ-VER-20/22,
    ACE-31/49/51/68/103/110, DTC-049, DTC-056, DTC-061, DTC-064, DTC-067,
    DTC-072, DTC-081, DTC-086, DTC-087, DTC-088, DTC-089, DTF-055, DTF-056.
- **Modelo de forma:** el plan de la fase de población fundadora
  (`plan-implementacion-f5-0.1.md`), y en particular su sección 11, que es el
  precedente de cómo este proyecto documenta un **cambio de contrato posterior
  al cierre de una fase**. La sección 10 de este plan copia esa forma.
- **Qué entrega.** Dos cosas que el usuario pidió, en este orden:
  1. **Microsatélites que varían entre individuos**: loci de repetición
     declarados sobre una referencia, con unos pocos alelos por locus —cada uno
     con **su secuencia**, no con un contador— y sus **frecuencias**, que el
     simulador sortea entre los individuos y entre las dos copias de cada uno,
     publicados en la verdad y en el VCF.
  2. **Tartamudeo visible en las lecturas**: el cambio del número de unidades
     observado por deslizamiento de la polimerasa, primero en la biblioteca sin
     amplificar y después en la amplificada por PCR, con sus duplicados y su
     herencia por linaje.

  Y una tercera que el responsable añadió como **condición de cierre**, no como
  ampliación aplazada: el **generador de alelos por deslizamiento biológico**,
  que produce el abanico de alelos y sus frecuencias en vez de exigir que el
  usuario los escriba. Mientras no exista, la simulación de microsatélites no
  puede darse por cerrada.
- **Esto es un hito propio, no un adelanto.** La primera versión lo presentaba
  como «adelantar pasos que ya existen». Ya no basta: tres piezas decididas el
  20 de septiembre **no existían en ningún plan** y dos de ellas tocan fases
  cerradas —el generador de referencia simulada con regiones definidas y la
  publicación de la correspondencia de coordenadas—. El plan conserva las
  coordenadas de los pasos que sí anticipa, para no perder la trazabilidad, y
  abre coordenadas nuevas donde no hay nada que anticipar.
- **Qué no entrega, y dónde queda:** sección 13.
- **Prohibiciones que este plan hereda y no relaja.** No autoriza *commit* ni
  publicación externa. No revierte trabajo ajeno. No se etiqueta una puerta con
  un comando que no se haya ejecutado. Las copias de seguridad de cualquier
  mutación van **fuera de git**; nunca `git checkout --`.

### 1.1. Cuatro hechos que este plan afirmaba y son falsos

Se enumeran porque tres de ellos **abarataban o encarecían trabajo**, y uno
cambia lo que puede citarse como fuente. Los dos primeros los he comprobado yo
mismo en el repositorio; los dos últimos proceden de la revisión adversarial,
que los verificó contra las fuentes publicadas, y los recojo atribuidos a ella
porque no he repetido esa comprobación.

1. **«Insertar un tracto sería la primera operación del proyecto que mueve
   coordenadas» — falso, y abarata el trabajo.** La reducción de referencia
   **ya reubica** cada segmento: `src/reference_reduction.rs`, en la función
   que materializa el plan de reducción, lleva un desplazamiento acumulado por
   secuencia y construye un `derived_interval` distinto del `source_interval`;
   y `schemas/reference-reduction-map-v1.schema.json` **ya publica los dos
   intervalos por segmento**, con su orden de materialización y su orientación.
   La tabla de correspondencia que la decisión S-03 exige es por tanto una
   **extensión de un mapa que existe y está publicado**, no maquinaria nueva.
   Comprobado el 20 de septiembre de 2026.
2. **`tabix` es 1.19, no 1.21.** El entorno tiene `bcftools` 1.21, `samtools`
   1.21 y `tabix` (HTSlib) **1.19**. La puerta externa exige 1.24 para las dos.
   Comprobado ejecutando `--version`.
3. **«Los dos conjuntos de la calibración pasaron por las mismas células» —
   falso.** Son **dos individuos distintos**, de proyectos distintos: uno de
   Papuasia Nueva Guinea secuenciado con PCR a 25,9× y otro de Pakistán sin PCR
   a 37,6×, y el primero además en una corrida agrupada. Que cada columna sea un
   **total** sigue en pie, y con ello el orden forzoso de la sección 3.1; lo que
   se debilita es atribuir la diferencia **solo** a la PCR, porque puede llevar
   dentro efectos de individuo y de lote. Las tasas por ciclo derivadas de esa
   resta son, por tanto, más toscas de lo que el texto daba a entender.
   Corregido también en `design/calibracion-tartamudeo-0.1.md`.
4. **La fuente de la dependencia con la longitud del tracto no está tras muro
   de pago.** Fungtammasan y colaboradores (2015) está en acceso abierto en
   Europe PMC con el identificador PMC4417121, y la curva ajustada **se puede
   recuperar** de su figura 1A. El pendiente P-1 deja de ser «no accesible» y
   pasa a ser «accesible, sin digitalizar todavía». Y la licencia, que la
   primera corrección dejaba por comprobar, **se conoce: CC BY-NC**, que
   permite reproducir la curva en este proyecto mientras su uso no sea
   comercial. Corregido también en el documento de calibración.

Y una quinta corrección, menor pero con consecuencias de licencia: **la fuente
de las cifras de duplicados no es CC BY**, declara todos los derechos
reservados; y **sí se sabe de qué plataforma es cada cifra** —el 1,5 % es una
plataforma MGI **con** PCR y el 10 %–13 % es Illumina **sin** PCR—, de modo que
no hay ninguna contradicción interna que explicar: son plataformas distintas.
La conclusión de no fijar un valor por omisión se sostiene, pero **por otra
razón**: el parámetro depende de la masa de ADN de entrada, de los ciclos y del
laboratorio.

---

## 2. Tabla de supuestos

**Cómo se lee, y para qué sirve ya.** Cada fila fue una decisión que los
documentos de origen dejaron abierta y que la primera versión de este plan
resolvió adoptando la recomendación de su documento, marcándola como supuesto.
**El 20 de septiembre de 2026 el responsable las cerró todas**, y la columna
«Opción adoptada» recoge ya el resultado: «Aceptado» cuando la recomendación se
mantuvo, y en mayúsculas cuando cambió. Los cinco cambios son S-02, S-04, S-05,
S-07 más S-14, y S-19; S-03, S-09 y S-29 se reformularon o ampliaron. El
razonamiento completo de cada una está en «Decisiones confirmadas por el
responsable», y **manda ese texto, no esta tabla**, que es un índice.

La columna «Pasos que la usan» sirve para lo mismo que antes: dice qué hay que
rehacer si alguna vez se revisa la decisión.

Las diecisiete primeras proceden del expediente de repeticiones (fichas D-01 a
D-18, menos D-10, que es una confirmación y no una decisión). Las once últimas,
de la sección 11 del diseño de bibliotecas con PCR.

| # | Qué se decide, en llano | Opción adoptada | Ficha | Pasos que la usan |
|---|---|---|---|---|
| **S-01** | Dónde se declaran los loci y los alelos, y qué versión de fichero de configuración nace | Tabla propia en el fichero general que ya se usa, sin documento enlazado nuevo. **Aceptado**, ver «Decisiones confirmadas» | D-01 (a) | F6.1.1.1, F6.4.8 |
| **S-02** | Cómo se guarda el alelo dentro del estado genotípico | ~~Modificación pequeña con el número de unidades~~ → **RESUELTO el 20/09/2026, ver «Decisiones confirmadas»**: variante pequeña nueva que **lleva la secuencia** del tracto, no un contador | D-02, corregida | F6.4.1, F6.4.3, F6.4.3.3, F6.5.1.1 |
| **S-03** | De dónde sale el tracto que hay en la referencia | **Reformulado el 20/09/2026**: las dos vías, con el orden del generador simulado fijado y la inserción en referencia real admitida con tabla de correspondencia, ver «Decisiones confirmadas» | D-03 (c), con (a) primero | F1.7.1, F1.7.2, F6.4.2 |
| **S-04** | Cómo se asignan los alelos a los individuos | **INVERTIDO el 20/09/2026**: frecuencias por locus **ahora**; asignación explícita de individuos concretos, después, ver «Decisiones confirmadas» | D-04 (a) ahora, (c) después | F6.1.1.1, F6.4.3.2, F6.4.11 |
| **S-05** | Qué pasa con un alelo que se separa 50 pb o más de la referencia | **CORREGIDO el 20/09/2026**: las repeticiones **no** se someten al umbral de 50 pb, con corrección de la especificación y un tope propio configurable, ver «Decisiones confirmadas» | D-05 (a) | F6.4.3.1, F6.4.8 |
| **S-06** | Pureza e interrupciones | **RESUELTO el 20/09/2026, ver «Decisiones confirmadas»**: la pureza siempre se calcula, nunca se declara; las interrupciones vienen de la referencia si es real, o se escriben si es simulada | D-06, reformulada | F6.4.1, F6.4.2, F6.4.8 |
| **S-07** | Por qué extremo se ganan o se pierden unidades | **CORREGIDO el 20/09/2026**: el resbalón ocurre en la **tirada ininterrumpida más larga**, no en el extremo derecho, ver «Decisiones confirmadas» | D-07 (a) | F6.4.1.1, F4.8.2, F4.8.4, F4.10.4 |
| **S-08** | Forma de la verdad publicada | Caso de inserción o deleción con metadatos de repetición, en versión nueva del esquema de verdad. **Aceptado**, ver «Decisiones confirmadas» | D-08 (a) | F6.4.3.3 |
| **S-09** | Campos del VCF y umbral de «tamaño manejable» | **RESUELTO el 20/09/2026, ver «Decisiones confirmadas»**: campos de repetición del estándar, alelos explícitos siempre, **sin umbral de tamaño** | D-09, corregida | F6.5.1.1 |
| **S-11** | Ámbitos de aleatoriedad e identificadores | Identificador del locus declarado por el usuario y validado por unicidad; catálogo técnico nuevo para el tartamudeo. **Aceptado**, ver «Decisiones confirmadas» | D-11 | F4.8.1, F6.4.3.2, F4.10.1, F6.4.11 |
| **S-12** | Si cabe un artefacto de PCR en una fase que declara no usar PCR | Nombrarlo por su mecanismo en el caso sin amplificar. **Aceptado**, ver «Decisiones confirmadas» | D-12 (a), ampliada | F4.8.1, sección 10 |
| **S-13** | En qué punto exacto del recorrido se aplica el tartamudeo | Después de la quimera y antes de calidad y errores. **Aceptado y verificado en el código**, ver «Decisiones confirmadas» | D-13 | F4.8.4 |
| **S-14** | El modelo de probabilidad | **CORREGIDO el 20/09/2026**: forma lineal confirmada por la literatura, pero el predictor pasa a ser la tirada más larga; estructura del paso, publicada, ver «Decisiones confirmadas» | D-14 | F4.8.2, F4.8.3 |
| **S-15** | Qué pasa con la longitud del fragmento | La molécula cambia de longitud y el inserto con ella, sin compensar. **Aceptado**, ver «Decisiones confirmadas» | D-15 (a) | F4.8.5 |
| **S-16** | Compatibilidad con lo ya publicado | Congelación acotada a escenarios sin loci declarados. **Aceptado**, ver «Decisiones confirmadas» | D-16, corregida | F4.8.6, F4.9.10 |
| **S-17** | Si la interfaz Marimo es condición de cierre | Interfaz como condición de cierre, para tartamudeo **y** loci. **Aceptado, con alcance mayor del previsto**, ver «Decisiones confirmadas» | D-17, corregida | F6.4.9.1, F4.8.9, F4.9.12 |
| **S-18** | Si esto es un hito propio o adelanta el hito de repeticiones | Adelanta pasos existentes y añade hitos nuevos. **Aceptado, con la advertencia de que hoy se parece más a un hito propio**, ver «Decisiones confirmadas» | D-18 (b) | Sección 4, sección 10 |
| **S-19** | Cómo se expresa la cantidad de material de partida de la PCR | **CAMBIADO el 20/09/2026**: el usuario declara la **fracción de duplicados**, obligatoria y sin valor por defecto; el parámetro interno se deriva, ver «Decisiones confirmadas» | PCR 1 (a) | F4.9.3, F4.9.4, F4.9.13 |
| **S-20** | Cómo actúa el sesgo por composición | Rechazo con tope de intentos, aceptación del último **contada**, aplicado al elegir la molécula. **Aceptado**, ver «Decisiones confirmadas» | PCR 2 (a) | F4.9.6 |
| **S-21** | Genealogía dentro de la familia | Progenitor uniforme y divergencia geométrica, declarada como aproximación. **Aceptado**, ver «Decisiones confirmadas» | PCR 3 (a) | F4.10.2 |
| **S-22** | Heterogeneidad de eficiencia entre moléculas | Clases de eficiencia desde el principio, opcionales y sin valores por defecto. **Aceptado, y respaldado por Best y col.**, ver «Decisiones confirmadas» | PCR 4 (a) | F4.9.3, F4.9.7 |
| **S-23** | Quimeras con PCR | Quimeras por pareja, la heredada para más adelante. **Aceptado, y respaldado por Kebschull y Zador**, ver «Decisiones confirmadas» | PCR 5 (a) | F4.9.7 |
| **S-24** | Cómo se dan las tasas de resbalón por ciclo | Derivadas de las dos columnas de la calibración y del número de ciclos. **Aceptado**, ver «Decisiones confirmadas» | PCR 6 (a) | F4.10.3 |
| **S-25** | Perfil tecnológico | Recurso hermano, dejando intacto el de sin PCR y su huella. **Aceptado**, ver «Decisiones confirmadas» | PCR 7 (a) | F4.9.1 |
| **S-26** | Versiones de esquema | Las versiones nuevas las abre este trabajo y la fase de variación compleja las completa mientras no se publiquen fuera. **Aceptado**, ver «Decisiones confirmadas» | PCR 8 (a) | F6.1.1.1, F4.9.3 |
| **S-27** | La cota de recursos | Se publica en el manifiesto, sin rechazo. **Aceptado y verificado**: la fase de secuenciación nunca ha tenido control de recursos, ver «Decisiones confirmadas» | PCR 9 (a) | F4.9.9 |
| **S-28** | Cobertura de moléculas distintas en el manifiesto | Se añade la cobertura de moléculas distintas. **Aceptado**, ver «Decisiones confirmadas» | PCR 10 (a) | F4.9.8 |
| **S-29** | Reparto en entregas | **AMPLIADO el 20/09/2026**: cuatro entregas más una quinta, el generador de alelos por deslizamiento biológico, ver «Decisiones confirmadas» | PCR 11 (a), refinada | Sección 3 |

**Por qué falta S-10.** La numeración de los supuestos sigue a la de las fichas
del expediente: S-01 corresponde a D-01, S-02 a D-02, y así hasta S-18. **No
hay S-10** porque la ficha D-10 no es una decisión, sino la confirmación que
sigue. Los supuestos S-19 a S-29 corresponden, en orden, a las once decisiones
numeradas 1 a 11 de la sección 11 del diseño de bibliotecas con PCR.

**Confirmación registrada, que no es un supuesto.** La ficha D-10 no decide
nada: DTC-064, la decisión que fija cómo se normaliza la proyección al VCF, ya
obliga a alinear los *indels* a la izquierda, de modo que un alelo de
repetición quedará anclado en el extremo izquierdo del tracto, y ya obliga a un
único registro por cromosoma y posición con varios `ALT`. Un diploide con 5 y 7
unidades produce **un** registro con dos alternativas y genotipo faseado. No
hay que tocar nada. Se anota para que no sorprenda en F6.5.1.1.

**Aviso de recuento.** El encargo cuenta veintiséis decisiones abiertas: quince
del expediente más once del diseño con PCR. El expediente anuncia en su parte B
que «tres» de sus dieciocho fichas son confirmaciones, pero **solo marca una**,
D-10. Por prudencia, este plan trata como supuesto todo lo que no está marcado
como confirmación, y de ahí salen veintiocho filas en vez de veintiséis. Si el
usuario identifica las otras dos confirmaciones, esas filas se bajan de la
tabla sin más consecuencia.

### 2.1. Decisiones **sin recomendación**, que este plan deja pendientes

No se inventan. El plan trabaja alrededor de ellas y dice cómo.

| # | Qué falta | Por qué no hay recomendación | Cómo se planifica alrededor |
|---|---|---|---|
| **P-1** | **Cómo crece la probabilidad de resbalón con la longitud de la tirada ininterrumpida.** Es el coeficiente de la pendiente del modelo lineal | **Corregido el 20/09/2026:** la primera versión decía que la fuente estaba tras muro de pago. **Es falso**: Fungtammasan y colaboradores (2015) está en acceso abierto en Europe PMC, identificador PMC4417121, y la curva se puede recuperar de su figura 1A. Lo que falta ya no es acceso ni licencia —es **CC BY-NC**, que no impide reproducirla en un uso no comercial— sino **digitalizar la curva**. **Actualizado el 25/09/2026:** se intentó y el pendiente queda **bloqueado por dos cosas, no por una**. El historial de Galaxy donde el artículo publicaba los valores exactos **ha desaparecido** —comprobado por tres vías, con el detalle en §9 de la calibración—, de modo que la tabla leída a ojo es hoy lo mejor que hay de esa fuente. Y sobre todo: **el artículo no dice cuántos ciclos de PCR usó**, así que aunque se digitalizara la figura con precisión perfecta, sus cifras de extremo a extremo **no se pueden convertir** en la pendiente por ciclo que el modelo necesita. Registrado como **DTF-061** para cuando aparezca una fuente utilizable. Y la forma sí está publicada y adoptada: Raz y col. (2019) fijan la dependencia **lineal por ciclo**. **CERRADO EN SU PARTE PRINCIPAL el 27/09/2026, por otra vía.** El responsable adoptó la curva como **tabla normalizada en el tope del tramo medido** —15 pb en mononucleótidos, 34 en dinucleótidos—, que usa solo su **forma** y no su nivel. Con eso el número de ciclos que el artículo no declara **deja de hacer falta**: se midió probando de 6 a 30 ciclos y la pendiente con la longitud sale idéntica, de modo que los ciclos afectan al nivel y no a la forma. Y se confirmó que la palabra «ciclo» aparece **cero veces** en el texto completo del artículo, con lo que esa pregunta queda cerrada del todo. Queda en **DTC-093** y en §11 de la calibración. **Lo que sigue abierto es otra cosa, y conviene no confundirla**: el **coeficiente lineal de la tirada ininterrumpida** sigue en cero y sigue publicándose como no calibrado. Es un parámetro distinto del factor de longitud de tracto | El coeficiente existe y es configurable, **arrancando en cero**, con lo que la probabilidad se reduce a la fracción medida por longitud de motivo, que sí está calibrada. El manifiesto lo publica como no calibrado. F4.8.2 prueba con un fixture que la función **sí** depende de la tirada más larga cuando el coeficiente es positivo, y F4.8.3 que con el coeficiente en cero reproduce la tabla |
| **P-2** | **Qué versión de `bcftools` y `tabix` vale para la puerta externa del VCF.** La decisión vigente exige **1.24 para las dos** | **RESUELTO el 20/09/2026 instalando lo que el contrato pedía**, en vez de decidir entre las opciones. Se compilaron HTSlib 1.24 y `bcftools` 1.24 en `/usr/local`, que precede a `/usr/bin` tanto en el camino por defecto como en el de `sudo`; no se desinstaló nada, porque `vcftools` depende del paquete del sistema | Ninguno. La puerta externa **sí** puede ejecutarse y etiquetarse como superada. Además se encendió la prueba de normalización contra `bcftools` que llevaba ignorada desde que se escribió DTC-066: la batería pasa de 1126 a 1127 y de doce ignoradas a once |
| **P-3** | **Los valores de las probabilidades de contracción y expansión** por longitud de motivo | **RESUELTO el 20/09/2026 en su forma, no en sus valores.** Se buscó fuente para los pesos y **no hay ninguna adoptable**: las cifras publicadas están repartidas entre trabajos distintos, sobre marcadores humanos concretos, y juntarlas en un vector sería **síntesis propia disfrazada de fuente**. La estructura sigue fijada por Raz y col. (2019) | Los pesos siguen **obligatorios y sin valor por omisión**, pero la documentación lleva un **ejemplo trabajado** marcado como ilustrativo, con sus fuentes y su advertencia de procedencia. Ver «El ejemplo trabajado de los pesos» más abajo |

| **P-4** | **Cómo produce el generador de deslizamiento biológico un abanico de alelos *con sus frecuencias*.** Aplicar un servicio de pasos N veces produce **una** secuencia de alelos, no una distribución | **RESUELTO el 20/09/2026.** El responsable adopta la forma candidata tal cual. Sobre la pregunta 4, su criterio: que el alelo de la referencia se pierda o no **depende de los parámetros** que se den —con media de pasos baja la mayoría de linajes no se mueven—, de modo que no se fuerza su presencia | Ninguno. F6.4.11 a F6.4.13 **desbloqueados**; la entrega 5 es ejecutable |
| **P-5** | **El valor por omisión del tope propio de tamaño del tracto** | **RESUELTO el 20/09/2026: 200 unidades.** El responsable confirma el valor propuesto. Razón registrada: un microsatélite de trabajo rara vez pasa de 30 o 40 unidades, de modo que 200 no limita ningún uso previsto y actúa solo como red contra un error de escritura; y el tracto se materializa en memoria y se publica en la verdad y en el VCF, así que un tope enorme permitiría generar ficheros desmesurados por descuido. Quien quiera simular expansiones patológicas —cientos o miles de unidades— debe subirlo explícitamente | Ninguno |

**Un desajuste menor que se registra en vez de taparse.** Los predictores de Raz
y col. incluyen el **tipo** de motivo —distinguen `A` de `C`, y `AC` de `AG` y
de `AT`—, mientras que la calibración de este proyecto es **por longitud de
motivo**. Es un desajuste conocido, no un error: se anota aquí y en el
manifiesto, y no se corrige inventando una tabla por tipo de motivo que nadie
ha medido para este uso.

### 2.2. Contradicciones encontradas al adoptar las recomendaciones

Se declaran aquí porque son hallazgos, no obstáculos. Las tres primeras cambian
lo que hay que construir.

**C-1. «Tartamudeo desactivado por omisión» contra «un escenario con loci debe
tartamudear sin PCR».** La ficha D-16 y los requisitos propuestos REQ-DAT-80 y
REQ-CONF-132 piden el tartamudeo desactivado por omisión y la identidad byte a
byte con la ejecución anterior. El diseño con PCR corrige esa congelación el
mismo día en que se escribe: vale para los escenarios que **no declaran loci de
repetición**, que son todos los que existen hoy; uno que sí los declare debe
tartamudear también sin amplificar, porque el fenómeno existe sin PCR y está
medido en la columna sin PCR de la calibración. Y hay una razón aritmética
además de una empírica: si el caso sin PCR se deja en cero, el caso con PCR
queda solo con el incremento despejado y **subestima** el tartamudeo entre un
14 % y un 28 % según el motivo.

*Resolución adoptada:* manda la corrección. La congelación se acota a los
escenarios sin loci declarados.

> **Nota del 24 de septiembre de 2026: la segunda mitad de esta resolución
> quedó superada.** Decía que con loci declarados el tartamudeo de base estaría
> **activo** por omisión y que lo desactivado sería el heredable de PCR. **El
> responsable lo descartó**, y en su lugar el tartamudeo de base se activa con
> un **interruptor explícito del evento técnico, apagado por omisión**: quien
> declara el microsatélite decide también si tartamudea y aporta sus
> parámetros.
>
> Las dos razones, registradas en la nota de REQ-CONF-132 de la
> especificación. Con activación implícita no existiría «este locus, **sin**
> tartamudeo», y con esa configuración desaparece la única prueba que el paso
> F4.8.1.1 exige para distinguir una guarda que mira las dos configuraciones de
> una que mira solo el locus. Y los pesos del paso son **obligatorios y sin
> valor por omisión**, de modo que todo escenario con un microsatélite quedaría
> obligado a declarar cuatro cifras que nadie ha calibrado, incluidos los de la
> entrega 1, que hoy no las declaran.
>
> **La primera mitad sigue vigente**: la congelación acotada a los escenarios
> sin loci declarados es la que se comprueba y la que la evidencia archivada
> sostiene. Y **REQ-DAT-80 y REQ-CONF-132 ya están escritos y aprobados**; la
> fila que los proponía en la sección 11 queda tachada como superada.

**C-2. Dos formas distintas del mismo modelo.** La calibración §5 propone la
forma de tres parámetros de HipSTR —probabilidad de añadir, probabilidad de
quitar y un parámetro geométrico del tamaño—. La ficha D-14 propone
`stutter-v1`: un ensayo de Bernoulli y un sorteo categórico sobre pasos
declarados. Además, D-14 afirma que «no existe en el repositorio ninguna
calibración de tartamudeo» y por eso recomienda valores predeterminados
**cero**; esa premisa **ya no es cierta**: `design/calibracion-tartamudeo-0.1.md`
existe, con fuentes verificadas y licencia CC BY 4.0.

*Resolución adoptada:* se conserva el mecanismo de D-14, porque un categórico
sobre pasos declarados **puede** expresar la forma geométrica de HipSTR sin
atarse a ella, y porque encaja con las primitivas normativas que el proyecto ya
tiene. Se abandona la parte de D-14 que recomendaba valores cero: las tasas por
longitud de motivo salen de la calibración, con su etiqueta «humano, Illumina»
viajando con ellas. Los pesos del tamaño del paso siguen sin valor (P-3).

**C-3. D-12 recomienda no reabrir la ampliación aplazada; el usuario la ha
reabierto.** La ficha D-12 recomienda nombrar el tartamudeo por su mecanismo
para evitar tocar DTF-056, la ampliación que aplaza la preparación amplificada.
El diseño con PCR parte de que el usuario ya decidió simular bibliotecas
amplificadas, y por tanto **divide** DTF-056. No es un choque real, sino un
solapamiento de fechas: la recomendación (a) de D-12 sigue siendo la forma
correcta de nombrar el tartamudeo **sin** PCR, que existe en los dos perfiles, y
la reapertura cubre la preparación amplificada. El plan usa las dos cosas y lo
dice en la sección 10.

**C-4. La ficha D-18 nombra menos coordenadas de las que el corte toca.** Dice
«adelantar F6.4.1 y F6.4.3». El corte necesita además el paso de la referencia
con el tracto (F6.4.2), el del registro en configuración efectiva y manifiesto
(F6.4.8), el de la interfaz (F6.4.9), un subconjunto del paso de contratos
(F6.1.1) y un subconjunto del paso de proyección VCF (F6.5.1). Se resuelve
nombrándolos todos en la sección 4.

**C-5. Un límite conocido que la recomendación original alcanzaba pronto, y que
la decisión del responsable elimina.** DTC-089 fija un techo duro de 65536
bytes al fichero general. Declarar cada alelo a mano cuesta de 70 a 90 bytes;
con 100 individuos diploides y 5 loci son entre 70 y 90 kB, es decir, **se
rompe el techo antes de llegar a un estudio de tamaño normal**. Esta
contradicción es una de las razones por las que el responsable **invirtió**
S-04: se declaran unos pocos alelos por locus con sus **frecuencias**, y el
simulador sortea los individuos, con lo que el coste del fichero deja de crecer
con el número de individuos. **Queda resuelta**, y el paso F6.1.1.1 la
comprueba con un documento próximo al techo para que no vuelva por otra puerta.

**C-6. La asimetría de D-17 deja la interfaz inservible para lo que pide el
usuario.** D-17 recomienda control Marimo para el tartamudeo y lo difiere para
los loci, con el argumento de que el catálogo de loci se escribe una vez y los
parámetros del tartamudeo se tocan muchas. El argumento es bueno pero el
resultado no funciona: **el tartamudeo solo actúa sobre loci declarados**, de
modo que una interfaz con controles de tartamudeo y sin forma de declarar un
locus no permite producir ni una sola lectura tartamudeada. El plan corrige la
ficha: la interfaz incluye una tabla mínima de loci y de alelos (F6.4.9.1), y
la tabla completa —parámetros aleatorios y regionales— sigue difiriéndose al
cierre del hito de repeticiones de la fase de variación compleja.

---


---

### 2.3. Decisiones confirmadas por el responsable

Lo que sigue **ya no es supuesto**: está decidido y sustituye a la fila
correspondiente de la tabla de arriba.

#### S-02 — el alelo guarda su secuencia, no un contador (20/09/2026)

**Decidido.** El alelo de un microsatélite en una copia haplotípica se guarda
como una **clase nueva del mundo de las variantes pequeñas** —junto al cambio
de base, la inserción y la deleción, y **no** como segmento estructural— que
lleva el **intervalo de referencia** que sustituye, la **secuencia realizada**
del tracto, y su anotación: motivo, número de unidades y las interrupciones.
El número de unidades y la pureza se **derivan** de la secuencia; no se
almacenan por separado, para que no puedan contradecirla.

**Por qué se corrigió la recomendación.** Guardar solo el número de unidades no
puede expresar un tracto **impuro** —uno con bases que rompen el motivo—, y la
mayoría de los microsatélites reales lo son. Además contradecía la decisión
D-06, ya adoptada, que hace al usuario declarar las interrupciones una a una:
con un contador no hay dónde ponerlas. La propia ficha D-06 planteaba un tercer
hueco, qué le pasa al patrón de interrupciones cuando el alelo cambia de
longitud, **que con un contador no tiene respuesta posible** y con la secuencia
no necesita ninguna regla: el patrón es el que la secuencia diga.

**Lo que lo hace barato.** `InsertionModification` (`src/segments.rs`) **ya
guarda secuencia** (`bases: Arc<[u8]>`), de modo que la estructura no estrena
capacidad. Y el trío intervalo de referencia + secuencia alternativa +
anotación **es el modelo del propio VCF** para una repetición en tándem, así
que la proyección es directa y no una traducción.

**Beneficio científico, no solo de representación.** Los tractos interrumpidos
**tartamudean menos** que los puros, y la posición de la interrupción importa.
Con la secuencia eso es modelable; con un contador no hay nada que mirar. Como
el tartamudeo es el objetivo del encargo, pesa.

**Declaración compacta, almacenamiento exacto.** Declarar cada alelo como una
cadena de bases en el TOML reventaría el límite de tamaño del fichero general
(contradicción C-5). Por eso se separan las dos cosas: por locus se declara un
puñado de alelos —motivo, número de unidades y su lista de interrupciones si la
tiene— y a los individuos se les **asigna** un alelo por índice o por
frecuencia. El programa materializa la secuencia desde esa definición.

**Frontera de alcance, confirmada.** Este trabajo **varía tractos que ya existen
y están declarados**. La aparición de microsatélites **nuevos** por inserción o
duplicación es trabajo estructural futuro y queda fuera.

**Consecuencia pendiente de ajustar:** el presupuesto de tamaño de salida debe
contar las bases del tracto en la verdad y en el VCF. Es un término pequeño,
pero el contrato de recursos es estricto y no puede descubrirse a mitad.


#### S-09 — campos del estándar, alelos explícitos, y ningún umbral (20/09/2026)

**Decidido.** El VCF usa los **campos de repetición en tándem del propio
estándar** —secuencia del motivo, longitud del motivo, número de copias, bases
totales y número de unidades distintas— declarados en la cabecera, y publica
los alelos **siempre de forma explícita**, con la secuencia completa en `REF` y
`ALT`. **No se fija ningún umbral de tamaño**, ni por longitud del tracto ni
por tamaño de la diferencia. No hay entrada complementaria, porque no hay
pérdida.

**Por qué se retira el umbral de 50 pb que la ficha proponía.** Se leyó la
especificación de VCF 4.5 (25 de febrero de 2026, 50 páginas) y **no contiene
ningún umbral de tamaño** para esta elección. Su criterio, literal en la
sección de repeticiones en tándem, es otro:

> «When the exact sequence is known, the variant can be represented as a
> non-symbolic ALT allele. […] When the exact sequence is not known, or when
> reporting tandem repeat "summary" information, the variant can be represented
> as `<CNV:TR>` copy number variants […]»

Es decir, el estándar distingue por **si se conoce la secuencia exacta**, no
por su longitud. La única cifra que aparece en ese contexto es un ejemplo con
un tracto de diez mil bases, representado de forma simbólica **por ser un
resumen, no por ser largo**.

**Un simulador siempre conoce la secuencia exacta: la ha generado él.** Nunca
hay incertidumbre que resumir, de modo que la forma explícita aplica siempre y
la simbólica nunca. El umbral de 50 pb era una convención de variantes
estructurales importada de fuera del estándar, que además no aplica a este
caso.

**Lo que el estándar no tiene, y cómo se cubre.** No hay campo de pureza ni de
interrupciones. El estándar expresa un tracto impuro **descomponiéndolo en
unidades consecutivas**, y de sus longitudes se calcula la pureza. Aun así, la
afirmación de la ficha de que «usando el estándar no hay pérdida de detalle»
solo es cierta **gracias a [S-02](#)**: como el alelo lleva su secuencia,
`REF` y `ALT` contienen el tracto base por base y no puede perderse nada; los
campos de repetición son anotación encima. Con la versión descartada —guardar
solo un contador— sí habría habido pérdida real.

**Consecuencia:** el hito de variantes complejas y complemento del VCF sigue
fuera del alcance, que era lo que esta ficha decidía.

**Fuente:** [especificación VCF 4.5](https://samtools.github.io/hts-specs/VCFv4.5.pdf),
sección de repeticiones en tándem. Sustituye a la nota del plan que declaraba
sin verificar los nombres y tipos de estos campos: quedan verificados.


#### S-06 — la pureza siempre se calcula; las interrupciones vienen de donde venga la referencia (20/09/2026)

**Decidido, y reformulado respecto a la ficha.** La ficha planteaba elegir entre
declarar las interrupciones, declarar una pureza con colocación determinista, o
declarar una pureza con colocación sorteada. El responsable señaló que la
pregunta está mal puesta, porque **hay dos situaciones distintas según de dónde
venga la referencia**, y el código ya las distingue: `ReferenceOrigin`
(`src/reference.rs`) tiene las variantes `Generated`, `Supplied` y `Reduced`.

**Con referencia real** —suministrada o reducida a partir de un genoma real—
los microsatélites **ya están ahí** con las interrupciones que tengan. No se
declaran: se declara **la posición del locus y su motivo**, y la secuencia se
**lee de la referencia**. La pureza es entonces un valor observado.

> **Validación que esto obliga a añadir.** Al declarar un locus con su motivo
> hay que **comprobar que en esa posición hay realmente un tracto de ese
> motivo**, y rechazar la configuración indicando el locus si no lo hay. Sin
> esa guarda se anotaría como microsatélite un trozo cualquiera del genoma.

**Con referencia simulada** el tracto no existe hasta que se crea, así que se
declara qué insertar, **pudiendo dar la secuencia exacta** con sus impurezas.

> **Corrección del 01/10/2026: la guarda de admisión corre también con
> referencia simulada.** Hasta ese día se saltaba con el argumento de arriba,
> «el tracto no existe hasta que se crea». Desde que la preparación de la
> referencia coloca el tracto (S-03, corrección del mismo día), el tracto
> **sí existe** antes de simular, y saltarse la guarda dejaba publicar tractos
> que el propio programa rechaza en cualquier otro sitio: `ACG` sobre ocho
> bases coloca `ACGACGAC`, de pureza 0,75, que se publicaba aunque un alelo
> declarado idéntico se rechazara con el mínimo de 0,8 y aunque el mismo FASTA
> usado como referencia suministrada no pasara la guarda; y con el generador,
> `AGAT` sobre siete bases publicaba como alelo ancestral `AGATAGA`, de pureza
> 0,57, sin aviso. Ahora la misma guarda —pureza mínima y tirada más larga
> mínima, con los mínimos del locus o los generales— juzga el tracto colocado,
> antes de simular y nombrando el locus. Un tracto de un motivo sin copia
> parcial, y uno compuesto con sus bloques enteros, tienen pureza 1 y la
> tirada que declaran, de modo que solo se rechaza lo que de verdad no cumple.

**Lo común, que es la decisión de fondo: la pureza es siempre una salida, nunca
una entrada.** En el caso real se calcula de lo que hay; en el simulado, de lo
que se escribió. En ninguno hace falta inventar una regla de colocación, que
era justo el problema que la ficha intentaba resolver y que así **desaparece**.
Una configuración que declare una pureza objetivo se rechaza antes de simular,
indicando el locus y el parámetro.

**Por qué no la colocación sorteada.** Tiene un defecto silencioso: si las
interrupciones se sortean, cambiar la semilla **mueve la secuencia de la
referencia**, y con ella todo lo demás.

**Argumento biológico a favor de declarar la secuencia y no un número.** Dos
alelos con el **mismo número de unidades** pero interrupciones en distinto
sitio son alelos **distintos**, y en genética de poblaciones y en forense esa
diferencia se usa precisamente para separar alelos que por longitud serían
indistinguibles. Un valor de pureza no puede expresarlo; una secuencia sí.

**Frontera con el simulador.** La referencia fija **el tracto**; el **número de
cambios de ese microsatélite entre los individuos lo genera el simulador**. Son
dos cosas separadas: una es el punto de partida, la otra la variación que se
simula sobre él.

**Cuestión que se traslada al paso del tartamudeo.** Cuando un alelo gana o
pierde una unidad, ¿dónde se añade o se quita respecto a las interrupciones?
Con la secuencia guardada la pregunta tiene respuesta posible; queda anotada
para ese paso y no se resuelve aquí.

**Nota sobre el umbral de tamaño, cerrada por el responsable.** Cómo un programa
de análisis posterior interprete un cambio —variación del número de
repeticiones o rotura de secuencia— **no es asunto de este simulador**. Se
producen las lecturas con las secuencias correctas y el VCF con la verdad;
acertar al interpretarlas es precisamente lo que un simulador sirve para poner
a prueba.


#### Bloque «cómo se declara todo» — S-01, S-03, S-04, S-11 y S-26 (20/09/2026)

**Aceptados sin cambio:** dónde se declaran los loci (tabla propia en el fichero
general que ya se usa), el identificador del locus (lo pone el usuario, se
valida que no se repita, sin clase de identificador automático nueva) y las
versiones de esquema (las abre este trabajo, la fase de variación compleja las
completa mientras no estén publicadas fuera).

##### S-03 — de dónde sale el tracto, reformulado

**Las dos vías, y con el orden fijado.** Con **referencia simulada** el orden es:
número de cromosomas y tamaño → **regiones definidas con su secuencia, que es
donde entran los microsatélites** → relleno al azar con los parámetros de
composición. **Los cambios estructurales quedan fuera de este trabajo** y van
con el hito estructural, **después** de generar la secuencia.

> **Aviso de coste, comprobado.** El generador **de partida**
> (`src/reference_generation.rs`) **solo sabía rellenar al azar**: no tenía
> regiones definidas ni orden. Lo anterior es un generador **nuevo**, en una
> fase cerrada (F1), no una reordenación de lo existente.
>
> *Cumplido el 23/09/2026 con la entrega 1.* El generador ya coloca regiones
> declaradas antes del relleno, y con la lista vacía produce los mismos bytes
> que antes. La reapertura está registrada en la sección 14 del plan del núcleo
> reproducible.
>
> **Corrección del 01/10/2026: el generador existía, pero el producto no lo
> usaba.** Hasta ese día la preparación de la referencia
> (`prepare_reference`, en `src/execution.rs`) llamaba al generador **sin
> regiones**: el generador con regiones y su tipo de región declarada solo
> tenían llamadores en sus propias pruebas. Con referencia simulada, el FASTA
> publicado llevaba **bases al azar** en el tramo de cada microsatélite
> declarado, todas las copias recibían un alelo alternativo y el alelo de la
> referencia no lo podía tener nadie. Lo encontró la revisión adversarial del
> 30/09/2026. Desde el 01/10/2026 la preparación coloca el tracto de cada
> locus **de un solo motivo** antes del relleno: el tramo declarado, de
> `start` a `end` inclusive, se rellena con el motivo repetido desde su primera
> base, y si su longitud no es múltiplo del motivo la última copia queda
> parcial. Es decisión del responsable de ese día; el tramo ya fija cuánto mide
> el tracto, así que no hace falta campo nuevo. **Los loci compuestos
> (`motifs`) declaran `reference_units`**, también por decisión del
> responsable del mismo día: la lista de unidades de cada bloque en la
> referencia, en el orden de los motivos, **obligatoria** en un compuesto
> sobre referencia simulada. Cada valor es un entero de al menos una unidad, y
> la suma de unidades por longitud de cada motivo tiene que ser **exactamente**
> la longitud del tramo; si no, la configuración se rechaza antes de simular
> nombrando el locus, la suma y el tramo. Cada bloque se coloca con su motivo
> repetido, en orden. El campo **se rechaza** en un locus de un solo motivo y
> con referencia suministrada, donde no haría nada. Las unidades son enteras
> en todos los bloques por la misma razón por la que lo son en los alelos
> compuestos: el número de unidades por bloque de un alelo se declara entero
> (`DeclaredRepeatAllele::units` en `src/configuration.rs`, leído como entero
> no negativo en `src/adapters.rs`) y el tracto de bloques se materializa con
> `TandemRepeatTract::pure_blocks`, que recibe copias enteras y rechaza el
> bloque vacío. Una unidad fraccionaria solo aparece **derivada** de unas bases
> declaradas, y en un compuesto las bases no se pueden declarar. El campo se
> publica en la configuración efectiva v6, de donde lo lee el juez que
> regenera la referencia simulada, y lo ofrece el editor de loci de la
> interfaz. Cabe en la general v4 y la efectiva v6 sin versión nueva por la
> decisión S-26.
>
> **La guarda de admisión juzga también el tracto colocado.** Un tramo que no
> es múltiplo del motivo deja una copia parcial que baja la pureza; si la baja
> por debajo del mínimo, o si el tramo no llega a la tirada mínima, la
> configuración se rechaza antes de simular, como con una referencia real. Ver
> la corrección de S-06 del mismo día.
>
> **Lo que cambia en el consumo de azar**, dicho porque es observable: el
> relleno consume una palabra por base **no declarada**, y el flujo de azar es
> uno solo para todas las secuencias, en su orden. Un tracto de `L` bases hace
> que todo lo que viene **después** de él —el resto de su secuencia y las
> secuencias siguientes— salga desplazado `L` palabras respecto de la
> referencia de antes; lo que viene **antes** no cambia. Una referencia
> simulada **sin** microsatélites no coloca nada y sale idéntica byte a byte a
> la de antes de la corrección.

**Con referencia real** se admite **insertar** el tracto en la posición
indicada, además del caso en que ya esté. Dos consecuencias que hay que
sostener:

1. **Las posiciones que declara el usuario se refieren al genoma de partida**,
   y el sistema las traduce a las del genoma final y **publica esa
   correspondencia**. Sin esa regla, declarar dos inserciones seguidas es
   ambiguo. Es válido porque el resultado es **un solo genoma final** y todo lo
   posterior habla en sus coordenadas.
2. **~~Sería la primera operación del proyecto que mueve coordenadas.~~
   Corregido el 20/09/2026: es FALSO, y abarata el trabajo.** La reducción de
   referencia **ya reubica** cada segmento: en `src/reference_reduction.rs`, la
   función que materializa el plan de reducción lleva un desplazamiento
   acumulado por secuencia y construye un `derived_interval` distinto del
   `source_interval`; y `schemas/reference-reduction-map-v1.schema.json` **ya
   publica los dos intervalos por segmento**, con su orden de materialización,
   su clase y su orientación. Es verdad que `src/coordinates.rs` no tiene un
   tipo de correspondencia entre dos versiones de una referencia, pero el
   **documento** que hace ese papel existe, está publicado y tiene esquema. La
   tabla que esta decisión exige es por tanto una **extensión de ese mapa**, no
   maquinaria nueva. El paso F1.7.2 la construye así.
3. El genoma final **no es el fichero suministrado**, así que se publica como
   **referencia derivada con identidad propia**, igual que las reducciones. De
   otro modo los resultados quedarían atribuidos a un genoma que no existe.

##### S-04 — frecuencias ahora, individuos concretos después

Se **invierte** la recomendación. Por locus se declaran unos pocos alelos con
sus **frecuencias**, y el simulador **sortea** los individuos. Declarar a mano
cien individuos por dos copias y por locus no es generar, es transcribir, y
revienta el límite de tamaño del fichero general (contradicción C-5).

La asignación **explícita de individuos concretos con combinaciones de alelos
específicas** queda para más adelante.

##### El deslizamiento biológico, y por qué no es el tartamudeo

El responsable aclaró que la variación entre individuos nace de **deslizamiento
biológico** —el de la replicación del ADN a lo largo de las generaciones—, no
del tartamudeo. **Son el mismo mecanismo físico en sitios distintos y con
consecuencias opuestas:** uno crea variación real heredable y va en la
**verdad**; el otro solo ensucia la medida y va en las **lecturas**.

La prueba de que no pueden confundirse: si los alelos verdaderos salieran del
modelo de tartamudeo, **cambiar el número de ciclos de PCR cambiaría el
genotipo real de los individuos**. Un parámetro de laboratorio decidiendo la
biología.

> ### Condición de cierre, no ampliación aplazada
>
> Se acuerda un **generador de alelos por deslizamiento biológico**: parte del
> alelo de la referencia, aplica pasos de una unidad y produce el abanico de
> alelos con sus frecuencias, que después se sortean entre los individuos. Es
> biología, va en la verdad y tiene su **propio ámbito de azar**, separado del
> tartamudeo. Es **una sola tirada para inicializar la diversidad**, no una
> simulación de generaciones: el modelo escalonado aplicado generación tras
> generación pertenece a la fase de evolución y no se invade.
>
> **Se hace en una segunda vuelta, pero la simulación de microsatélites NO
> puede darse por cerrada hasta tenerlo.** No es una ampliación aplazada —que
> podría no hacerse nunca— sino una **condición de cierre**. Hasta entonces los
> alelos y sus frecuencias los declara el usuario.


#### S-30 — un locus admite **varios motivos**, en bloques (20/09/2026)

**Decidido después de la primera parada de la entrega 1**, al hacerse efectivos
los mínimos de admisión y descubrirse que rechazaban un caso legítimo.

Un locus pasa a admitir una **lista ordenada de bloques**, cada uno con su motivo
y su número de unidades, en vez de un motivo único. Un locus de un solo motivo es
el caso de un bloque y **se comporta exactamente igual que antes**.

**El caso que lo motiva es real.** Un microsatélite compuesto del tipo
`(TCTA)n(TCTG)n(TCTA)n` —la forma de muchos marcadores forenses— declarado con un
solo motivo sobre el intervalo entero mide pureza 40/60 y lo rechaza el mínimo por
omisión de 0,8, cuando en realidad es perfectamente regular: lo que pasa es que su
regularidad tiene tres tramos.

**Y acerca el modelo al estándar en vez de alejarlo.** El VCF 4.5 ya lo contempla:
su ejemplo literal es `(CAG)10(TG)7(CAGG)3`, con un campo que cuenta cuántas
secuencias repetidas tiene el alelo y **listas paralelas** de secuencia del
motivo, longitud del motivo, número de copias y bases totales. Estábamos usando un
subconjunto de esos campos sin aprovechar para qué fueron diseñados. La
configuración efectiva los publica ya con esa forma, de modo que la proyección al
VCF sea un volcado y no una traducción.

**Las tres reglas:**

1. **La tirada ininterrumpida más larga se mide sobre todos los bloques y gana la
   más larga, sea del motivo que sea.** Es lo biológicamente correcto: la
   polimerasa resbala donde hay más repetición seguida, no donde el usuario puso
   el primer motivo.
2. **El paso del tartamudeo gana o pierde una unidad del bloque donde ocurre el
   resbalón**, es decir, del motivo de la tirada más larga.
3. **La pureza se sigue midiendo sobre bases** y sobre el tracto completo,
   sumando las copias completas de todos los bloques.

**Lo que no cambia:** la descomposición voraz, con su dirección fijada de
izquierda a derecha, sigue valiendo **dentro de cada bloque**.

**Una frontera que la decisión impone y conviene ver.** La **estructura en
bloques se declara; no puede derivarse**. Dadas unas bases sueltas no hay forma de
saber dónde acaba un motivo y empieza el siguiente. Por eso la forma «secuencia
exacta» del alelo solo se admite con **un** motivo: con varios habría que adivinar
los cortes. Todo lo numérico —copias, pureza, interrupciones, tiradas— se sigue
derivando de las bases dentro del tramo de cada bloque.

#### S-31 — los mínimos de admisión, también por locus (20/09/2026)

**Decidido junto con S-30.** La pureza mínima y la tirada mínima pueden
declararse **dentro de la tabla del locus**, y los generales quedan como **valor
de partida** para los que no los declaren. Los valores por omisión siguen siendo
**0,8 y 3 unidades**.

**El motivo: el mismo número no significa lo mismo según el motivo.** Tres
unidades de un mononucleótido son tres bases, que aparecen por azar
constantemente; tres unidades de un tetranucleótido son doce, que no. Con un
umbral único, o deja pasar basura en los motivos cortos o rechaza loci legítimos
en los largos.

**Ampliación posible, registrada y no hecha:** en el futuro podría ajustarse
**por longitud de motivo**, de modo que el suelo saliera de una tabla en vez de
declararse locus a locus. Hoy no se hace.

#### S-32 — una modificación pequeña dentro de un microsatélite rechaza el escenario (22/09/2026)

**Decidido.** Si un cambio de base, una inserción o una deleción declarados caen
**dentro del intervalo de un locus de repetición**, el escenario se **rechaza
antes de simular**, nombrando el locus y:

- **la modificación exacta**, cuando se declara en coordenadas de la
  referencia: los cambios de base y las filas de la matriz de variantes;
- **los candidatos**, cuando se declara en coordenadas vigentes: las
  inserciones y las deleciones. El diagnóstico dice que **una** inserción o
  deleción anterior rompió el tramo y **enumera las declaradas sobre esa
  copia**, sin afirmar cuál de ellas fue.

La precisión es del 22 de septiembre de 2026 y está justificada más abajo, en
«por qué con un indel se nombran los candidatos».

**Por qué.** El alelo de un microsatélite describe **el tracto tal como está en
la referencia**: su motivo, su número de unidades, su pureza y su tirada más
larga se derivan de esas bases. Si otra modificación ya tocó ese tramo, el alelo
seguiría describiendo un tracto que **ya no existe**, y la anotación publicada
en la verdad y en el VCF sería **falsa**. No sería un error visible: sería un
dato incorrecto que alguien se cree, que es la clase de fallo que este proyecto
persigue.

**Consecuencia práctica, dicha en voz alta.** Un usuario que quiera introducir
una variante puntual dentro de un microsatélite **no puede hacerlo declarando un
cambio de base**: tiene que expresarlo como parte del alelo del locus, con su
secuencia o con su lista de interrupciones. Es más restrictivo, y es
deliberado — la alternativa es publicar anotaciones que no describen lo que hay.

**Cómo se aplica.** En **dos sitios que cubren mitades distintas**, y ninguno
basta solo:

1. **En la frontera de configuración**, sobre lo que el usuario escribió y
   antes de materializar a nadie. Ahí **no hay sorteo**, de modo que el rechazo
   no puede depender de la semilla. Juzga los cambios de base y las filas de la
   matriz, suministrada o generada.
2. **En la cadena de segmentos**, consultada **para toda copia, se aplique o no
   el alelo**, exigiendo que el intervalo del tracto llegue intacto. Juzga las
   inserciones y las deleciones.

La segunda condición es la que quita la dependencia de la semilla de raíz:
comprobar solo al aplicar el alelo hacía que el escenario se aceptara o se
rechazara **según qué alelo saliera sorteado**, y cuando se aceptaba la verdad
publicaba que la copia «conserva el tracto de la referencia» sin conservarlo.

**Por qué con un indel se nombran los candidatos.** Porque **no todas las
modificaciones se cuentan con la misma regla de medir**. Un cambio de base
declara su posición en **coordenadas de la referencia**; una inserción o una
deleción declaran su punto en las **vigentes**, las de la secuencia después de
los cambios anteriores. Una inserción cambia la longitud, así que a partir de
la segunda las dos reglas se separan.

Comparar el punto de un indel contra el intervalo del locus es entonces
comparar dos reglas distintas, y **se equivoca en las dos direcciones**:
reproducido el 22 de septiembre de 2026, con una inserción previa de 60 bases,
un indel declarado en el punto 210 cae de verdad en la 150 —dentro del tracto—
y la comparación lo dejaba pasar; y otro declarado en el punto 150 cae de
verdad en la 90 —fuera— y la comparación **rechazaba un escenario válido**
diciéndole al usuario que su cambio estaba dentro del microsatélite.

De ahí que los indeles no se juzguen en la frontera, sino en la cadena de
segmentos, que mira el resultado ya montado y **no tiene deriva posible**. El
precio es que ahí ya no se distingue cuál de las modificaciones rompió el
tramo, y por eso se enumeran los candidatos.

**La alternativa que se descartó, y por qué.** Traducir el punto declarado a
coordenadas de la referencia descontando los indeles anteriores. Es posible,
pero significa **una segunda aritmética de coordenadas** conviviendo con la que
la cadena de segmentos ya sostiene, y dos implementaciones del mismo concepto
acaban separándose. La ronda de revisión anterior encontró justo eso en este
mismo trabajo: la pregunta «¿lleva este documento alelos de repetición?»,
escrita dos veces y respondiendo cosas distintas.

**Ampliación posible, registrada y no hecha:** nombrar el indel exacto
recomponiendo la copia modificación a modificación y comprobando la integridad
después de cada una, para quedarse con la primera que la rompe. Solo costaría
en el camino de error. Hoy no se hace: el usuario que recibe el diagnóstico
tiene delante su propio fichero con los indeles que declaró.

**Origen.** La regla la propuso el implementador al construir el alelo dentro
del estado haplotípico, la declaró para confirmar por no estar en el plan, y el
responsable la confirmó el 22 de septiembre de 2026.

#### Bloque «cómo se publica la verdad» — S-08 y S-05 (20/09/2026)

**S-08, aceptado sin cambio.** El alelo se publica como un caso de inserción o
deleción con metadatos de repetición —locus, motivo y número de unidades—, en
una versión nueva del esquema de verdad. Las versiones anteriores no cambian ni
un byte. Encaja con [S-02](#): como el alelo lleva su secuencia, el caso de
inserción o deleción ya la contiene y los metadatos solo la anotan.

##### S-05 — las repeticiones no se someten al umbral de 50 pb

**Precisión previa, porque se confundieron dos cosas.** Esta ficha **no** trata
la representación en el VCF, que quedó resuelta en [S-09](#) —siempre
explícita, sin umbral, porque el criterio del estándar es si se conoce la
secuencia—. Trata una **regla interna del proyecto**: REQ-CONF-55 y REQ-IND-31
fijan un umbral configurable, 50 pb por omisión, que separa inserciones y
deleciones pequeñas de variantes estructurales, y obligan a aplicarlo igual al
generar, al clasificar y al contar.

**Decidido: la opción (c), con corrección de la especificación**, en vez de la
(a) que la ficha recomendaba.

**Primera razón, biológica.** Un alelo `(CA)40` frente a uno `(CA)10` difiere en
60 bases, pero **no es una variante estructural en ningún sentido**: es un
polimorfismo de longitud de un microsatélite. El umbral existe para separar
inserciones de secuencia **nueva** de reordenamientos, y ahí tiene sentido.
Aplicarlo a un tracto repetido es un **error de categoría**: mide bases cuando
lo que define al fenómeno es que la misma secuencia se repite más veces.

**Segunda razón: la opción (a) aprieta más de lo que su ficha admitía.** Para
que el umbral muerda hace falta una diferencia respecto a la referencia de
**25 unidades en un dinucleótido**, 17 en un trinucleótido y 13 en un
tetranucleótido. La ficha afirmaba que «sobra»; para un tetranucleótido sí,
pero **25 unidades en un dinucleótido con rango alélico amplio no es una cifra
imposible** en una población diversa.

**Redacción a corregir en la especificación:** el umbral se aplica a
inserciones y deleciones de **secuencia nueva**; un alelo de repetición en
tándem **declarado sobre un locus declarado** no se somete a él. Queda así
cubierta la exigencia de coherencia de REQ-IND-31, que pedía aplicarlo igual en
los tres momentos: se aplica igual, y a las repeticiones no se aplica en
ninguno.

**Salvaguarda que sustituye al umbral.** Las repeticiones llevan su **tope
propio, configurable**, para que un error de dedo no genere un tracto
desmesurado. Rechaza antes de simular, identificando locus, individuo y copia.

> **Nota del plan, que no cambia la decisión.** La frase «identificando locus,
> individuo y copia» se escribió cuando los alelos se declaraban individuo a
> individuo. Con **frecuencias** —decisión S-04, tomada el mismo día— el alelo
> que rebasa el tope se conoce **al leer la configuración**, antes de que
> exista ningún individuo al que asignárselo, de modo que el diagnóstico solo
> puede nombrar **el locus y el alelo**. El paso F6.4.3.1 lo recoge así. Si el
> responsable prefiere conservar además individuo y copia, habría que mover la
> comprobación a después del sorteo, que es más tarde y menos útil.


#### Bloque «mecánica del tartamudeo» — S-07, S-12, S-13, S-14, S-15 y S-24 (20/09/2026)

**Aceptados sin cambio:** el punto del recorrido (después de la quimera y antes
de calidad y errores, **verificado en el código**: la quimera rematerializa el
fragmento y lo compara byte a byte, de modo que tocarlo antes la rompería); que
la molécula **cambie de longitud** y con ella el inserto observado, sin
compensar, porque así es como el fenómeno se detecta en datos reales; cómo se
derivan las tasas por ciclo de las dos columnas de la calibración; y nombrar el
artefacto por su mecanismo en el caso sin amplificar.

##### Dos correcciones, ambas con respaldo publicado

Se buscó si existía un modelo publicado en vez de inventar la forma, y existe.

**Fuente principal: Raz O., Biezuner T., Spiro A. y col. (2019),** *Short tandem
repeat stutter model inferred from direct measurement of in vitro stutter
noise*, Nucleic Acids Research 47(5):2436–2445, doi:10.1093/nar/gky1318,
**CC BY 4.0**. Infieren el modelo midiendo el ruido *in vitro* directamente. Su
modelo «lineal, un paso arriba y hasta tres abajo» es una **cadena de Markov en
tiempo discreto donde cada ciclo de PCR es un paso**, con probabilidades de
transición **linealmente dependientes de la longitud actual del tracto**.

##### S-14 — la forma lineal se confirma; cambia el predictor

**Se retira una crítica equivocada.** Se objetó que la forma lineal era una
aproximación tosca, porque la literatura describe crecimiento exponencial con
meseta. **Son dos magnitudes distintas y ambas ciertas:** la dependencia lineal
es de la probabilidad **por ciclo**; la curva con meseta es la **fracción
observada de extremo a extremo**, y sale de acumular la primera sobre muchos
ciclos. La forma lineal del plan **es la publicada**.

**Lo que sí cambia: el predictor.** En genética forense está establecido que el
mejor predictor del tartamudeo es la **tirada ininterrumpida más larga** del
motivo, mejor que la longitud total del alelo, y se modela con **regresión
lineal sobre ella**. El plan usaba el número total de unidades. Pasa a usar la
tirada más larga.

**Y la estructura del tamaño del paso deja de ser inventada.** El plan la hacía
obligatoria y sin valor por defecto. La forma publicada es **asimétrica**:
contracciones de hasta tres unidades por ciclo y expansiones de una sola, con
sesgo hacia la contracción. Se adopta esa estructura; los **valores** siguen sin
fuente (ver más abajo).

##### S-07 — el resbalón ocurre en la tirada más larga, no en el extremo derecho

La propuesta era ganar o perder unidades **siempre por el extremo derecho**.
Para un tracto **puro** da igual: el resultado es idéntico. Para un tracto
**impuro** no, y la convención tenía dos efectos falsos: **vacía un bloque
mientras el otro queda intacto**, y mueve la pureza como consecuencia de una
decisión nuestra y no de la biología.

El mismo hallazgo de la tirada más larga lo resuelve: **el resbalón ocurre en la
tirada ininterrumpida más larga del motivo**, que es donde la polimerasa resbala
—las interrupciones **estabilizan** el tracto, y por eso los tractos
interrumpidos tartamudean menos—. En caso de empate, la primera por la
izquierda, para que siga siendo determinista. Se conserva la propiedad buena de
la propuesta original: el resto del tracto queda igual byte a byte y la
diferencia sigue siendo una inserción o deleción limpia.

Esto responde la pregunta que [S-06](#) había dejado aplazada.

##### Tres salvedades que se registran, no se ocultan

1. **Raz y col. no publican los valores de los parámetros**, solo el método de
   ajuste. Tenemos la **forma**, no los números: la dependencia con la longitud
   sigue sin calibrar y su coeficiente arranca en cero.
2. Sus predictores incluyen el **tipo** de motivo —distinguen `A`, `C`, `AC`,
   `AG`, `AT`—, no solo su longitud. La calibración de este proyecto es por
   longitud de motivo. Desajuste menor, anotado.
3. La literatura forense advierte que la tirada más larga **tiene
   contraejemplos** y que la variación en los flancos influye. **No es una ley:
   es el mejor predictor conocido.**
4. **Atribución corregida el 20/09/2026.** La observación de que, con
   secuenciación masiva, un mismo alelo puede producir tartamudeos originados
   en **partes distintas del tracto** —lo que obliga a refinar el concepto de
   «la» tirada más larga— es de **Vilsen y colaboradores (2018)**, no del
   artículo de *Genes* de 2017 al que la primera redacción se la atribuía.

**Fuentes:** Raz y col. (2019), NAR 47(5):2436–2445, CC BY 4.0, para la forma
del modelo; *Flanking Variation Influences Rates of Stutter in Simple Repeats*,
Genes 2017, PMC5704242, para la tirada más larga como predictor y la influencia
de los flancos; y Vilsen y col. (2018) para el matiz de los orígenes múltiples
en secuenciación masiva.


#### Bloque «la amplificación» — S-19 a S-28 (20/09/2026)

**Siete aceptados**, tres de ellos ahora con respaldo de las referencias
verificadas el 19 de septiembre:

- **Clases de eficiencia** entre moléculas, opcionales y sin valores por
  defecto. Best y col. (2015) encuentran que esa heterogeneidad es
  **independiente de la secuencia** y de las condiciones globales, de modo que
  modelarla como clase aleatoria por molécula no es arbitrario.
- **Quimeras por pareja**, la heredada por la familia para más adelante.
  Kebschull y Zador (2015) midieron que los cambios de plantilla son **raros y
  limitados a números de copia bajos**: el aplazamiento está respaldado por
  medición, no por comodidad.
- **Sesgo por composición** al elegir la molécula, con rechazo y tope de
  intentos. Benjamini y Speed (2012) establecen que lo que manda es la
  composición **del fragmento completo**. Se acepta el atajo de aceptar el
  último intento **porque queda contado**: el usuario puede ver cuántas veces
  se agotó el tope y juzgar si el sesgo le afecta.
- **Genealogía** con progenitor uniforme y divergencia geométrica, declarada
  como aproximación del coalescente y no como forma exacta.
- **Perfil tecnológico** como recurso hermano, dejando intacto el de sin PCR y
  su huella.
- **Cota de recursos publicada, no impuesta.** Verificado: la fase de
  secuenciación **nunca ha tenido** control de recursos —los límites viven solo
  en la vía de haplotipos—, así que imponerlo aquí estrenaría algo que no le
  corresponde.
- **Cobertura de moléculas distintas** en el manifiesto: es lo que un analista
  compara tras eliminar duplicados.

##### S-19 — el usuario declara la fracción de duplicados, no «parejas por molécula»

**Cambiado.** La propuesta pedía declarar `mean_pairs_per_molecule`. Esa es la
variable del modelo, no la del usuario: un biólogo piensa en el porcentaje de
duplicados que ve en su informe de calidad. Se invierte: **se declara la
fracción de duplicados** y el sistema deriva el parámetro interno y **publica la
fracción realizada**. Es además el patrón que el proyecto ya usa en todas
partes: declarar un objetivo y publicar lo realizado con su tolerancia.

##### Por qué no se fija un valor por defecto, pese a haberlo buscado

Se buscó uno defendible y **no lo hay**, por una razón que conviene dejar
escrita para que nadie lo intente de nuevo sin saberlo.

Shen H., Liu P., Li Z. y col. (2019), *Advanced Whole Genome Sequencing Using an
Entirely PCR-free Massively Parallel Sequencing Workflow*, bioRxiv
doi:10.1101/2019.12.20.885517, miden a 30× con tres réplicas: bibliotecas **con
PCR, ~1,5 %**; bibliotecas **sin PCR en plataformas Illumina, 10,43 %–13,62 %**.
Duplicados identificados con Picard.

> **Dos correcciones del 20 de septiembre, y una de ellas invalida el
> razonamiento original.**
>
> **Primera, de licencia:** este preprint **no es CC BY**; declara todos los
> derechos reservados. Se cita, no se reproduce. La primera redacción le
> atribuía CC BY 4.0 y es falso.
>
> **Segunda, de hecho:** la primera redacción decía que «no se puede determinar
> de qué plataforma es el 1,5 %» y construía sobre eso una contradicción
> aparente —«las bibliotecas sin PCR salen con más duplicados que las que sí la
> llevan»— que resolvía conjeturando cómo marca Picard. **Sí se puede
> determinar, y no hay contradicción**: el 1,5 % es una plataforma **MGI con**
> PCR y el 10 %–13 % es **Illumina sin** PCR. Son plataformas distintas, así
> que las cifras no son comparables entre sí y no hay nada que explicar. El
> razonamiento original era erróneo.

**La conclusión se mantiene, por la razón buena.** No se fija valor por omisión
**porque el parámetro depende de la masa de ADN de entrada, del número de
ciclos y del laboratorio**, no porque las cifras publicadas midan otra cosa. Y
sigue siendo cierto y útil distinguir en la documentación entre la fracción
**atribuible a la amplificación**, que es la nuestra, y la que marca una
herramienta de análisis contando parejas que caen en el mismo sitio, que
incluye además las copias de la química de agrupación: quien compare ambas sin
saberlo no las verá cuadrar nunca.

**Lo que sí queda establecido:** la tasa **no depende apreciablemente de la
profundidad** en el rango normal —midieron de 10× a 46× y se mantiene—, así que
un valor por omisión no necesitaría ir atado a una profundidad. Ese obstáculo
queda descartado.

**El ejemplo trabajado, con una salvedad que apareció al revisarlo.** La fuente
da **dos** valores para bibliotecas amplificadas, **1,5 % y 0,9 %**, no uno. Por
tanto la diferencia «medio punto atribuible a la amplificación» que se anotó en
una versión anterior **no puede darse por buena sin comprobar a qué condición
corresponde cada uno de los dos valores y cuál es el término sin amplificar con
el que se compara**. El paso que escriba la documentación tiene esa
comprobación como requisito previo: **o cita la condición exacta de cada cifra,
o no usa el ejemplo**. Y en cualquier caso viajará con dos salvedades: que es
de una plataforma **MGI y no Illumina**, y que es un ejemplo, no un valor por
omisión. Consagrar otra cifra sin comprobar sería repetir el error que este
plan lleva dos rondas corrigiendo.

**Decisión:** parámetro **obligatorio y sin valor por omisión**, con ese ejemplo
trabajado en la documentación.

**La mejor fuente para este número son los informes de calidad del propio
laboratorio**, no una media de genomas humanos ajenos.


#### Bloque «proceso y alcance» — S-16, S-17, S-18 y S-29 (20/09/2026)

**Aceptados en el fondo**, pero las decisiones del 20 de septiembre **han
cambiado su tamaño**, y eso se registra aquí para que el plan no se lea como
algo menor de lo que es.

**S-16.** Congelación byte a byte acotada a los escenarios que **no** declaran
loci de repetición. Sin cambio.

**S-17 — la interfaz mínima es el doble de lo previsto.** Cuando se escribió la
ficha, «mínima» significaba una tabla de loci y alelos. Con lo decidido hoy hace
falta: declarar loci con **motivo y posición**, declarar alelos con sus
**interrupciones**, y declarar **frecuencias** por locus. Sigue siendo
razonable, pero hay que presupuestarlo por lo que es.

**S-18 — esto ya no es solo «adelantar pasos que existen».** La ficha decía que
se aprovechaban pasos del plan de la fase de variación compleja, y en parte
sigue siendo cierto. Pero hoy se han añadido **tres piezas que no existían en
ningún plan**, y las tres tocan fases cerradas:

1. el **generador de referencia simulada con regiones definidas** —el de
   partida solo sabía rellenar al azar (`src/reference_generation.rs`);
   **cumplido el 23/09/2026**—;
2. la **tabla de correspondencia de coordenadas** entre el genoma suministrado
   y el genoma final con el tracto insertado. **Corregido el 20/09/2026:** la
   primera redacción decía que la reducción «nunca ha movido coordenadas», y es
   falso —ya publica intervalo de origen y derivado por segmento—, de modo que
   esta pieza es una **extensión de un mapa que existe**, y es la más barata de
   las tres, no la más cara;
3. la **corrección de la especificación** sobre el umbral de 50 pb.

El trabajo se parece más a un **hito propio** que a un adelanto, y el plan debe
usar esa palabra para que nadie lo subestime al leerlo.

**S-29 — cinco entregas, no cuatro.** Las cuatro del plan se mantienen:
microsatélites que varían entre individuos; tartamudeo de base sin amplificar;
preparación amplificada; y tartamudeo clonal heredable. Se añade una **quinta**:
el **generador de alelos por deslizamiento biológico**, que es condición de
cierre y no tenía sitio en ninguna de las cuatro. Se pone como entrega propia
para que se vea que existe en vez de quedar escondida dentro de otra.


### El ejemplo trabajado de los pesos del paso (20/09/2026)

**No es una calibración y no debe usarse como tal.** Es una forma de partida
para que quien escriba su primer escenario no invente de cero. Va en la
documentación del parámetro, **nunca** como valor por omisión: los pesos siguen
siendo obligatorios.

**Forma ilustrativa**, tomando como unidad el paso de una unidad hacia abajo:

| Paso | Peso relativo |
|---|---|
| **−1** unidad | 100 |
| **+1** unidad | 20 |
| **−2** unidades | 10 |
| **−3** unidades | 1 |

**De dónde sale cada orden de magnitud**, y por qué no es una calibración:

- El paso de **una unidad hacia arriba** se sitúa en torno al **1–2 %** de la
  altura del alelo y el de **dos unidades hacia abajo** en el **0,5–1 %**, según
  Aponte, Gettings, Duewer, Coble y Vallone, *Sequence-based Analysis of Stutter
  at STR Loci: Characterization and Utility*, Instituto Nacional de Estándares y
  Tecnología de Estados Unidos, sobre dos marcadores comparando secuenciación
  masiva y electroforesis capilar.
- El de **tres unidades hacia abajo** va de **0,04 % a 0,35 %** según el
  marcador, en la literatura forense de filtrado de tartamudeo.
- El de **una unidad hacia abajo**, dominante, suele situarse en el **5–10 %**
  en tetranucleótidos.

**Las tres advertencias, que pesan más que la tabla:**

1. **Cada cifra viene de un trabajo distinto**, sobre marcadores distintos.
   Componerlas en un solo vector es una **operación propia**, no un resultado
   publicado. Por eso la tabla da órdenes de magnitud relativos y no
   probabilidades.
2. **Son tetranucleótidos humanos forenses**, y este parámetro es **por longitud
   de motivo**. El desajuste es el mismo que impidió fijar un valor por omisión
   a la fracción de duplicados.
3. La propia fuente advierte que **la variación depende del motivo concreto** y
   que la tirada ininterrumpida **no lo explica todo**: dentro de una misma
   longitud de motivo hay dispersión entre marcadores. Un valor por omisión
   ocultaría esa dispersión.

**Lo que la tabla sí sostiene, y es lo útil:** que el paso de una unidad domina
en las dos direcciones; que **bajar una es unas cinco veces más frecuente que
subir una**, lo que realiza el sesgo de contracción que la calibración exige; y
que los pasos múltiples **caen muy deprisa**, un orden de magnitud por unidad
adicional. Esa forma sí está respaldada por las tres fuentes a la vez, aunque
los números concretos no.


## 3. El orden es forzoso, y por qué. El reparto en entregas

### 3.1. La aritmética manda el orden

Las dos columnas de la calibración —«con PCR» y «sin PCR»— son **totales, no
componentes**. Cada conjunto pasó por una celda de flujo con su amplificación
en puente y por la síntesis de la secuenciación, de modo que la columna sin PCR
**ya contiene** el mosaicismo somático de la muestra, el deslizamiento del
puente y el de la síntesis, y la columna con PCR contiene todo eso **más** la
PCR de biblioteca. Restar, y no sumar, es por tanto lo correcto.

> **Precisión que corrige una afirmación falsa de la primera versión.** El plan
> decía que «los dos conjuntos pasaron por las mismas células». **No es cierto:
> son dos individuos distintos**, de proyectos distintos —uno de Papuasia Nueva
> Guinea con PCR a 25,9× y otro de Pakistán sin PCR a 37,6×, el primero además
> en una corrida agrupada—. Que cada columna sea un total no depende de eso y
> sigue en pie, igual que el orden forzoso que se deduce abajo. Lo que sí se
> debilita es atribuir la diferencia **solo** a la PCR: la resta puede llevar
> dentro efectos de individuo y de lote. En consecuencia, **las tasas por ciclo
> derivadas son toscas**, el manifiesto debe publicarlas marcadas como tales, y
> el paso F4.10.3 lo dice en su propio texto.

El diseño con PCR no suma nada a la columna con PCR: **despeja** el incremento
propio de la amplificación,

```
r_pcr = 1 − (1 − r_con) / (1 − r_sin)
```

y lo compone después con la parte de base, de modo que
`1 − (1−r_sin)(1−r_pcr) = r_con`, exactamente la fracción medida.

**Consecuencia que fija el orden de este plan:** una preparación amplificada
entregada **sin** el tartamudeo de base se queda solo en el incremento
despejado y subestima el tartamudeo en un 14 % en homopolímeros, 28 % en
dinucleótidos, 26 % en trinucleótidos, 23 % en tetranucleótidos y 20 % en
pentanucleótidos. No es una preferencia de orden: es que la aritmética no
cierra. La preparación amplificada va **después** del tartamudeo de base.

> **El rango «14 % a 28 %» vale para motivos de 1 a 5 pares de base.** Los
> hexanucleótidos dan **7 %**, no porque el fenómeno sea distinto, sino porque
> su fila es la menos fiable de la tabla: el propio autor advierte que la
> diferencia entre sus dos columnas procede probablemente del ruido de haber
> analizado pocos loci. Siempre que este plan cite el rango, lo acota así.

### 3.2. Qué no depende de que existan los microsatélites, y qué sí

| Pieza | ¿Necesita loci de repetición? |
|---|---|
| Perfil amplificado, variante `Pcr`, configuración, esquemas y rechazos | **No** |
| Muestreo por familias, copias de los primeros ciclos, tamaño de familia | **No** |
| Duplicados y su oráculo por bytes | **No** |
| Aceptación por composición de bases | **No** |
| Recuentos del manifiesto, fracción esperada y realizada, cobertura de moléculas distintas, cota de recursos | **No** |
| Tartamudeo de base (agrupación y síntesis) | **Sí** |
| Genealogía de la familia y tartamudeo heredable | **Sí**, y además necesita el tartamudeo de base |
| Generador de alelos por deslizamiento biológico | **Sí**, y además necesita que el alelo lleve su secuencia |

### 3.3. El reparto: cinco entregas, no dos

El diseño con PCR recomienda **dos** entregas: primero la amplificación con sus
duplicados y su composición, que no depende de nada de la fase de variación
compleja; después el tartamudeo clonal. Este plan **respeta esa separación pero
reordena y subdivide**, y el responsable añadió una quinta el 20 de septiembre.
Quedan cinco:

| Entrega | Contenido | Coordenadas | Depende de |
|---|---|---|---|
| **1** | Microsatélites que varían entre individuos: referencia con el tracto, loci, alelos con su secuencia, sorteo por frecuencias, verdad, VCF, manifiesto, interfaz | F6.4.0, F1.7.1–F1.7.2, F6.4.1–F6.4.1.1, F6.1.1.1, F6.4.2, F6.4.3–F6.4.3.3, F6.5.1.1, F6.4.8–F6.4.8.1, F6.4.9.1 | Nada |
| **2** | Tartamudeo de base, sin amplificar | F4.8.1–F4.8.9 | Entrega 1 |
| **3** | Preparación amplificada: familias, duplicados, composición, cota | F4.9.1–F4.9.13 | Nada técnicamente; **ordenada aquí por la aritmética** |
| **4** | Tartamudeo clonal heredable por linaje | F4.10.1–F4.10.6 | Entregas 1, 2 y 3 |
| **5** | Generador de alelos por deslizamiento biológico, **condición de cierre**. **F6.4.11 a F6.4.13 desbloqueados el 20/09/2026 (antes bloqueados por P-4**; F6.4.4 no | F6.4.4, F6.4.11–F6.4.13 | Entrega 1; independiente de las entregas 2 a 4 |

**Por qué cinco y no dos.**

1. **Porque lo que el usuario pidió cabe en las dos primeras.** Microsatélites
   que varían entre individuos y tartamudeo visible en las lecturas se entregan
   al terminar la entrega 2, sin haber tocado el contrato de duplicados ni de
   sesgo por composición. Una entrega única que empezara por la amplificación
   retrasaría el resultado pedido detrás de un cambio de contrato grande.
2. **Porque la aritmética prohíbe el orden contrario.** Poner la amplificación
   antes del tartamudeo de base dejaría existir, aunque fuera unos días, una
   biblioteca amplificada que subestima el tartamudeo entre un 14 % y un 28 %.
3. **Porque la primera versión del esquema la abre quien la necesita.** El
   supuesto S-26 dice que las versiones nuevas las abre el trabajo que llegue
   primero. Con este orden, las abre la entrega 1 —que es la que tiene el
   contrato de entrada más grande— y las entregas 3 y 4 las completan dentro de
   la misma versión mientras no estén publicadas fuera.
4. **Porque la entrega 3 sigue siendo verificable sola**, que es la virtud que
   el diseño reclamaba para ella: sus trece pasos no citan ni un locus.
5. **Porque la quinta no cabía dentro de ninguna de las otras y no es
   opcional.** El generador de alelos por deslizamiento biológico es lo que
   convierte «el usuario escribe los alelos y sus frecuencias» en «el simulador
   los genera». El responsable lo declaró **condición de cierre**, no ampliación
   aplazada: la simulación de microsatélites no puede darse por terminada sin
   él. Se le da entrega propia precisamente para que se vea, en vez de quedar
   escondido como un paso más dentro de la primera.

**Por qué la quinta va la última y no la primera.** No es una preferencia: para
generar un abanico de alelos por deslizamiento hace falta que exista el alelo,
que lleve su secuencia y que haya dónde publicarlo, y eso es justo lo que
entrega la primera. Y va después de las entregas 2 a 4 solo por prioridad del
usuario —ver el tartamudeo antes—, no por dependencia técnica: si se
prefiriera, puede ejecutarse inmediatamente después de la entrega 1 sin tocar
nada de lo demás.

**Qué se pierde y cómo se compensa.** Si el usuario prefiriera el orden del
diseño —amplificación primero—, el plan sigue siendo válido intercambiando las
entregas 3 y 2, **con una condición inexcusable**: mientras falte el tartamudeo
de base, la configuración debe **rechazar** la combinación de preparación
amplificada con loci de repetición declarados, con un diagnóstico que nombre el
paso pendiente. Ese rechazo está planificado de todos modos en F4.9.10, porque
también hace falta en el orden recomendado mientras la entrega 4 no exista.

---

## 4. Cómo se numera, y puertas comunes

### 4.1. Numeración

Este trabajo cruza **tres** fases: la del núcleo reproducible y la del realismo
técnico, ambas cerradas, y la de variación compleja, sin empezar. La decisión
S-18 pide que se le llame **hito propio** y no adelanto, porque tres de sus
piezas no existen en ningún plan. A la vez, perder la trazabilidad con los
pasos que sí existen sería un retroceso. Se hacen las dos cosas:

- Donde este trabajo **anticipa** un paso que ya está escrito en el plan de la
  fase de variación compleja, conserva su coordenada: F6.4.1, F6.4.2, F6.4.3,
  F6.4.4, F6.4.8 y F6.4.9. Cuando solo necesita **una parte** de un paso mayor
  de aquel plan, usa una subdivisión: **F6.1.1.1** es el subconjunto de
  repeticiones del paso de contratos externos y **F6.5.1.1** el de la
  proyección al VCF.
- Donde **no hay nada que anticipar**, abre coordenada nueva: **F6.4.0** es la
  línea base de evidencias; **F6.4.11** a **F6.4.13** son el generador de
  alelos por deslizamiento biológico; **F1.7.1** y **F1.7.2** son el trabajo de
  referencia que cae en la fase del núcleo reproducible, cuyo último hito hoy
  es F1.6; y **F4.8**, **F4.9** y **F4.10** son hitos nuevos de la fase de
  realismo técnico, cuyo último hito hoy es F4.7.

**Tres coordenadas caen en fases cerradas** —F1.7.1, F1.7.2 y todo F4.8 a
F4.10— y por tanto exigen el cambio de contrato de la sección 10. No hay forma
de evitarlo: el tracto vive en la referencia, que es de la primera fase, y el
tartamudeo vive en el motor de lecturas, que es de la cuarta.

**Correcciones de dependencia que este plan declara** y que hay que aplicar al
plan de la fase de variación compleja (sección 11):

- **F6.4.3 ya no depende del primer paso del hito de variantes estructurales.**
  La cadena de segmentos (`src/segments.rs`) ya sostiene modificaciones tipadas
  por copia haplotípica, y una de ellas, la inserción, **ya guarda secuencia**
  en un `Arc<[u8]>`. El alelo de repetición es una clase hermana de esas, no un
  segmento estructural.
- **F6.4.2 ya no depende de derivar una referencia estructural.**
- **F6.5.1.1 no arrastra el fichero complementario del VCF**, porque los campos
  de repetición del estándar más el alelo explícito no pierden nada.
- **F6.4.4 se adelanta** y se usa para lo que la decisión S-04 exige: generar
  el abanico de alelos **una sola vez**, para inicializar la diversidad. Su
  aplicación generación tras generación sigue siendo de la fase de evolución y
  no se invade.

### 4.2. Puertas

Cada paso ejecuta sus pruebas específicas. Antes de cada pausa de hito se
ejecuta la batería completa, con los comandos literales:

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

- **Base de partida: se mide, no se cita.** Circulan dos recuentos distintos de
  la batería de Rust —1124 y 1126— en documentos distintos del proyecto, y este
  plan **no fija ninguno**: el paso F6.4.0 ejecuta la batería, registra el
  recuento exacto en la línea base y es contra ese número contra el que se
  comparan los hitos. Citar de memoria un recuento que nadie ha vuelto a medir
  es exactamente el error que este proyecto ya ha pagado.
- `./scripts/check-independent-validators.sh` **falla si el directorio ya
  existe** y debe dar **5 PASS**. Cada hito usa un directorio propio:
  `evidence/repeats-f6-4`, `evidence/stutter-f4-8`, `evidence/pcr-f4-9`,
  `evidence/clonal-f4-10` y `evidence/slippage-f6-4-11`.
- `uv run --frozen python -m unittest discover -s tests_ui` da hoy **96
  ejecutadas con 2 fallos preexistentes y ajenos** (interfaz del VCF de la fase
  2 y del evento técnico de la fase 4). No se persiguen. Tampoco las doce
  pruebas ignoradas históricas.
- **`pytest` no es puerta.** No está instalado ni declarado en
  `pyproject.toml`. Ningún paso de este plan lo cita.
- **La puerta externa del VCF no se puede etiquetar como superada.** Exige
  `bcftools` y `tabix` en la versión 1.24; el entorno tiene `bcftools` 1.21,
  `samtools` 1.21 y `tabix` (HTSlib) **1.19**. Es el pendiente P-2.
- **Regla de cierre de la interfaz** (hoja de ruta §7): toda opción
  configurable que este plan implemente exige control en Marimo y equivalencia
  interfaz → TOML → configuración efectiva antes de cerrar su hito.

### 4.3. Disciplina de pruebas

Cada paso dice **contra qué variante defectuosa concreta** debe fallar la
prueba que pide. Es la respuesta al defecto recurrente número uno del proyecto:
pruebas que no pueden fallar por el defecto que dicen cubrir, que apareció diez
veces durante la fase de población fundadora, que dejó al descubierto dos
guardas del verificador frente a la batería entera, y que la revisión
adversarial del 20 de septiembre volvió a encontrar **diez veces** en la
primera versión de este plan.

Reglas, sin excepción:

1. Toda regresión se verifica **fallando** contra la variante defectuosa.
2. Se comprueba que **la mutación mutó**: SHA del fichero antes y después.
3. La copia de seguridad va **fuera de git**. Nunca `git checkout --`.
4. Si una aserción **no puede ejercitarse**, se dice en su docstring, sin
   adornos, en vez de fingir que discrimina.
5. Cuando un fixture ordinario no ejercita la guarda —porque otra comprobación
   caza antes la discrepancia—, hace falta un fixture **saturado**. Los pasos
   que lo necesitan lo dicen y explican qué tiene que saturar.
6. **Una prueba que solo comprueba que la implementación hace lo que el diseño
   dice no valida el diseño.** Cuando una comprobación sea de ese tipo
   —tautológica—, se declara así en su docstring en vez de presentarla como
   validación del modelo. Afecta a F4.10.5.

---

## 5. Entrega 1 — Microsatélites que varían entre individuos

Resultado: una referencia con tractos de microsatélite en posiciones conocidas,
unos pocos alelos por locus **con su secuencia completa** y sus frecuencias,
sorteados entre los individuos y entre las dos copias de cada uno, y publicados
en la verdad, en el VCF y en el manifiesto, con control en la interfaz. Sin
tartamudeo todavía.

**Lo que cambia respecto de la primera versión de este plan, en una frase:** el
alelo ya **no es un número de unidades** sino una secuencia con su anotación;
los individuos ya **no se declaran uno a uno** sino que se sortean de unas
frecuencias; y **no hay umbral de 50 pb**, sino un tope propio configurable.

- ☐ **F6.4.0.** Línea base de evidencias, medida **antes** de tocar una sola
  línea de código y guardada **fuera de git**.
  - Depende de: ninguno.
  - Procede de: CLAUDE.md §4; prueba 1 de la sección 10 del diseño con PCR, que
    recuerda que medir los hashes después de tocar el código es circular.
  - Cambia: nada del repositorio. Crea `/tmp/gps-baseline-<fecha>/` con:
    `sha256sum` de los FASTQ finales, la verdad, el VCF de verdad, la
    configuración efectiva y los dos manifiestos de al menos dos ejecuciones de
    referencia; la salida completa de
    `cargo test --workspace --all-targets --locked` **con su recuento exacto de
    pruebas que pasan, fallan y se ignoran**; y la identidad del árbol de
    trabajo.
  - Se comprueba con: el fichero de hashes existe, tiene una línea por
    artefacto y ninguna ruta bajo el árbol de git; el recuento de pruebas
    registrado es el que devolvió el comando, no una cifra copiada de un
    documento.
  - **Cómo se identifica la línea base, corregido dos veces.** La primera
    versión anotaba el `HEAD`, que **no discrimina**: un cambio sin *commit*
    deja exactamente el mismo `HEAD`. La segunda propuso el **hash del árbol**
    del *commit*, `git rev-parse HEAD^{tree}`, y **tampoco discrimina**: ese
    hash es función del *commit* y de nada más, no mira el directorio de
    trabajo. Comprobado el 20 de septiembre de 2026 con dos ficheros
    modificados sin preparar: `git rev-parse HEAD^{tree}` y `git write-tree`
    devolvieron **el mismo hash** que antes de tocarlos.
  - **Lo que sí discrimina**, y es lo que el paso registra: **(1)** la salida de
    `git status --porcelain`, que **debe estar vacía**; **(2)** el `HEAD`; y
    **(3)** un **hash del contenido real de los ficheros en disco**, calculado
    sobre el árbol de trabajo y no sobre ningún objeto de git —por ejemplo
    `find src schemas resources ui scripts -type f | sort | xargs sha256sum |
    sha256sum`—. El tercero es el único que cambia cuando cambia un fichero sin
    prepararlo ni confirmarlo.
  - **Debe fallar contra:** una línea base capturada con el árbol ya modificado
    —la (1) y la (3) la cazan—; y, si alguien retirase la (1) creyéndola
    redundante, una capturada tras un cambio sin *commit*, que solo la (3)
    caza. *Se comprueba que la (3) discrimina* modificando un byte de un
    fichero cualquiera bajo `src/`, recalculando y exigiendo que el hash
    cambie, y restaurándolo después desde la copia de seguridad fuera de git.

- ☐ **F1.7.1.** Generador de referencia simulada **con regiones definidas**, en
  este orden: número de secuencias y tamaño → **regiones declaradas con su
  secuencia exacta, que es donde entran los microsatélites** → relleno al azar
  con los parámetros de composición que ya existen.
  - Depende de: ninguno.
  - Procede de: REQ-REF-19, REQ-CONF-56; decisión **S-03**.
  - Cambia: `src/reference_generation.rs`; escribe por el escritor FASTA
    canónico `src/reference_fasta.rs`, que no cambia.
  - **Por qué es trabajo nuevo y no una reordenación.** El generador **de
    partida** solo sabía rellenar al azar: no tenía regiones declaradas ni orden
    entre ellas. Cae en la fase del núcleo reproducible, que está cerrada, y por
    eso lleva coordenada propia y entra en el cambio de contrato de la sección
    10. *Cumplido el 23/09/2026; la reapertura está registrada en la sección 14
    del plan del núcleo reproducible.* **Corregido el 01/10/2026: cumplido en
    el generador, no en el producto.** Hasta ese día ninguna ejecución pasaba
    regiones al generador, y el tramo de un microsatélite en una referencia
    simulada salía al azar; ver la corrección del mismo día en S-03. Desde
    entonces la preparación de la referencia coloca el tracto de cada locus de
    un solo motivo, y el juez que regenera la referencia simulada desde la
    configuración efectiva coloca los mismos tractos por su cuenta, de modo
    que un producto que dejara de pasarlos se rechaza al publicar. Los loci
    compuestos se colocan bloque a bloque con las unidades que declara
    `reference_units`, también desde el 01/10/2026.
  - **Qué queda fuera, dicho aquí para que no se cuele:** los cambios
    estructurales —duplicaciones, inversiones, translocaciones— **no** entran.
    Van con el hito estructural y **después** de generar la secuencia.
  - Se comprueba con: una región declarada aparece **exactamente** en su
    coordenada, localizada con `grep -ob` externo al productor; el relleno
    aleatorio respeta las fracciones de base configuradas fuera de las regiones
    declaradas; misma semilla, mismo FASTA byte a byte; dos regiones declaradas
    que se solapan se rechazan antes de generar, nombrando las dos; una región
    que rebasa el final de la secuencia se rechaza; y —**esta es la que
    discrimina el orden**— un generador instrumentado confirma que el relleno
    consume exactamente **la longitud de la secuencia menos las bases de las
    regiones declaradas**, no la longitud entera.
  - **Debe fallar contra:** un generador que **rellene primero y sobrescriba
    después**. *Corregido:* las cinco comprobaciones de secuencia y rechazo que
    la versión anterior de este paso pedía **las cumple ese mutante sin
    despeinarse**, porque el resultado final es idéntico: la región queda en su
    sitio y la composición observable fuera de ella es la misma. El **único**
    observable que separa los dos órdenes es el **consumo de aleatoriedad**, y
    por eso la prueba del contador es obligatoria y no un adorno. Va en la
    docstring. También debe fallar contra una escritura que coloque la región
    una base más allá, que el `grep -ob` externo caza.

- ☐ **F1.7.2.** **Inserción de un tracto en una referencia real** y publicación
  del resultado como **referencia derivada con identidad propia**, con la
  **tabla de correspondencia** entre las coordenadas del genoma suministrado y
  las del genoma final.
  - Depende de: F1.7.1 (comparte el contrato de región declarada).
  - Procede de: REQ-REF-08–11 y REQ-MAN-07 por analogía con la referencia
    reducida; decisión **S-03**, puntos 1 y 3.
  - Cambia: `src/reference.rs`, `src/reference_reduction_io.rs` o un módulo
    hermano, y una **versión nueva** del esquema del mapa de coordenadas.
  - **Cuánto cuesta realmente, corregido.** La primera versión de este plan
    decía que sería «la primera operación del proyecto que mueve coordenadas».
    **Es falso.** La reducción de referencia ya reubica cada segmento con un
    desplazamiento acumulado por secuencia
    (`src/reference_reduction.rs`), y `schemas/reference-reduction-map-v1.schema.json`
    **ya publica** `source_interval` y `derived_interval` por segmento, con su
    orden de materialización y su orientación. Esta pieza es una **extensión de
    ese mapa a un segundo tipo de operación**, no maquinaria nueva. Es la más
    barata de las tres piezas que la decisión S-18 llama nuevas.
  - **Las dos reglas que hay que sostener.** Primera: las posiciones que
    declara el usuario se refieren **al genoma de partida**, y el sistema las
    traduce a las del genoma final y publica la correspondencia; sin esa regla,
    declarar dos inserciones seguidas es ambiguo. Segunda: el genoma final **no
    es el fichero suministrado**, así que se publica con identidad propia,
    igual que las reducciones, para que ningún resultado quede atribuido a un
    genoma que no existe.
  - Se comprueba con: tras insertar un tracto de 30 bases en la posición 1000
    de una secuencia, una posición declarada en 2000 sobre el genoma de partida
    aparece en 2030 en el mapa publicado, y la base que hay allí es la misma
    que había en 2000 del original, comparada byte a byte; **dos** inserciones
    seguidas en la misma secuencia dan un mapa en el que la segunda posición
    declarada sigue refiriéndose al original y no al intermedio; el mapa valida
    contra su esquema con `python3 -m jsonschema`; la referencia derivada tiene
    identificador propio y el manifiesto la nombra como derivada de la
    suministrada; el mapa de la **reducción** existente no cambia ni un byte.
    **Debe fallar contra:** una implementación que interprete la segunda
    posición declarada en coordenadas del genoma ya modificado —la prueba de
    las dos inserciones es la única que lo caza, y una sola inserción **no
    discrimina**, lo que va en la docstring—; una que publique el resultado con
    el identificador del FASTA suministrado.

- ☐ **F6.4.1.** Modelo de dominio del locus de repetición, puro: intervalo en
  la referencia, motivo A/C/G/T no vacío y sin `N`, **secuencia realizada del
  tracto**, y la anotación **derivada de esa secuencia** —número de unidades,
  pureza, lista de interrupciones y **tirada ininterrumpida más larga** del
  motivo—.
  - Depende de: ninguno. Es código puro sobre secuencias en memoria.
  - Procede de: REQ-CONF-57/59/60, REQ-IND-32/33/34, REQ-REF-23, ACE-31/110;
    decisiones **S-02** y **S-06**.
  - Cambia: módulo nuevo `src/tandem_repeat.rs`.
  - **La regla que gobierna el módulo entero:** nada de lo anotado se almacena
    por separado de la secuencia; **todo se deriva de ella**, para que no pueda
    contradecirla. La pureza es siempre una **salida**, nunca una entrada.
  - **El intervalo NO tiene por qué ser múltiplo de la longitud del motivo, y
    el número de unidades puede ser fraccionario.** Esto corrige una
    contradicción interna de la versión anterior de este plan, que aquí
    rechazaba «un motivo cuya longitud no divide el intervalo» mientras el paso
    del VCF exigía que «un tracto con unidad parcial produzca un valor no
    entero». No podían ser ciertas las dos.
  - **Y la que manda es esta, por biología.** Las interrupciones reales suelen
    ser de **una sola base** —una `T` dentro de un `(CA)n`—, con lo que el
    intervalo deja de ser múltiplo del motivo. Con la regla anterior, la guarda
    de admisión de F6.4.2 habría admitido el tracto real impuro y este modelo
    lo habría rechazado acto seguido. El propio estándar del VCF trae como
    ejemplo `(CAG)5(CA)1(CAG)4`, que un modelo de múltiplos exactos no puede
    representar.
  - **Definiciones exactas, para que no haya dos implementaciones razonables.**
    El **número de unidades** es la longitud del intervalo dividida por la del
    motivo, y puede tener parte fraccionaria. La **pureza** es la proporción de
    bases del intervalo que pertenecen a copias **completas y exactas** del
    motivo. La **tirada ininterrumpida más larga** se cuenta en copias completas
    y consecutivas. La descomposición en unidades que los campos del VCF piden
    la hace F6.5.1.1 a partir de esta misma secuencia.
  - Se comprueba con: `(ACG)×5` produce quince bases y anota cinco unidades,
    pureza 1, ninguna interrupción y tirada más larga 5; un tracto
    `AC AC AC GG AC AC` anota doce bases, **seis** unidades, pureza **10/12**
    sobre bases, una interrupción y tirada más larga **3**, no 6; **el ejemplo
    del estándar `(CAG)5(CA)1(CAG)4` se admite** y anota veintinueve bases,
    **29/3 = 9,667 unidades**, pureza **27/29** y tirada más larga **5**; un
    tracto `CACACATCACA` con una interrupción de **una sola base** se admite y
    anota **5,5 unidades**; motivo vacío y motivo con `N` se rechazan
    identificando el locus; una configuración que declare una **pureza
    objetivo** se rechaza nombrando el locus y el parámetro; dos tractos con el
    **mismo número de unidades** y la interrupción en sitio distinto producen
    anotaciones distintas y secuencias distintas.
  - **Debe fallar contra:** **un modelo que exija que el intervalo sea múltiplo
    del motivo** —los dos fixtures de unidad parcial lo cazan, y son los que
    esta corrección añade—; un número de unidades redondeado a entero, que el
    de 5,5 caza; una pureza contada sobre unidades en vez de sobre bases, que
    el fixture de la interrupción de una sola base distingue —contada sobre
    unidades daría 5/6, sobre bases 10/11—; un modelo que guarde la pureza
    recibida en vez de calcularla; un cálculo de la tirada más larga que
    devuelva el número **total** de unidades, que el fixture impuro distingue;
    una anotación que colapse dos alelos de igual longitud y distinta
    interrupción en el mismo objeto.

- ☐ **F6.4.1.1.** La operación «aplicar un paso de ±k unidades a un tracto»:
  las unidades se ganan o se pierden **en la tirada ininterrumpida más larga**
  del motivo, con desempate por la **primera por la izquierda**, y el resto del
  tracto queda **byte a byte igual**, de modo que la diferencia sigue siendo
  una inserción o una deleción contigua y limpia.
  - Depende de: F6.4.1.
  - Procede de: decisión **S-07**, que corrige la convención del extremo
    derecho; Raz y col. (2019) y la literatura forense para el hecho de que la
    polimerasa resbala donde el motivo corre sin interrupción.
  - Cambia: `src/tandem_repeat.rs`.
  - **Por qué se retiró el extremo derecho, y por qué la prueba de la primera
    versión era falsa.** Con la convención del extremo derecho, un tracto de
    cinco unidades con una interrupción en la tercera se compone de dos tiradas,
    [1–2] y [4–5]; quitar dos unidades **vacía una de ellas y deja la otra
    intacta**, mueve la pureza por una decisión nuestra y no por biología. La
    primera versión de este plan pedía además comprobar que «las bases
    anteriores al extremo derecho son idénticas», aserción que bajo la
    convención nueva es **falsa por construcción**: con las tiradas empatadas
    gana la izquierda, y las bases que cambian son justamente esas. Se retira.
  - Se comprueba con: en un tracto puro de diez unidades, quitar una deja las
    nueve restantes idénticas y la diferencia es una deleción contigua de una
    unidad; en un tracto con tiradas de 3 y 5 unidades, el paso se aplica en la
    de 5 y la de 3 queda **byte a byte intacta**, comprobado sobre las bases;
    con tiradas **empatadas** de 4 y 4, el paso cae en la de la izquierda, y el
    resultado es reproducible; aplicar `+1` y luego `−1` devuelve la secuencia
    original byte a byte. **Debe fallar contra:** la convención del extremo
    derecho —el fixture de tiradas 3 y 5, con la más larga a la izquierda, la
    caza—; un desempate no determinista; una implementación que regenere el
    tracto desde el motivo, que recolocaría las interrupciones. *Cuidado, en la
    docstring:* con un tracto **puro** todas las convenciones coinciden y
    **ninguna prueba discrimina**; los fixtures de este paso son impuros a
    propósito, y el de tiradas empatadas es el único que ejercita el desempate.

- ☐ **F6.1.1.1.** Contrato de entrada: esquema general **v4** y configuración
  efectiva **v6**. Declara, por locus: identificador, secuencia, **inicio y
  final del intervalo**, y motivo; una **lista corta de alelos**, cada uno con
  su número de unidades y sus interrupciones si las tiene, o con su secuencia
  exacta cuando la referencia es simulada; las **frecuencias** de esos alelos;
  y el **tope propio configurable** de tamaño del tracto.
  - **Ampliado el 01/10/2026, dentro de la misma versión (S-26):** un locus
    compuesto sobre referencia simulada declara `reference_units`, las
    unidades enteras de cada bloque en la referencia, y la configuración
    efectiva v6 las publica en el locus. Ver la corrección de S-03 de ese día.
  - **El usuario declara el intervalo entero, inicio y final, y no solo el
    inicio.** Hacía falta decirlo: la pureza y la tirada más larga se miden
    **dentro de un intervalo**, y sin final no hay dónde medirlas. La
    alternativa —que el programa extienda el tracto desde el inicio mientras el
    motivo siga apareciendo— es un algoritmo que nadie ha fijado, tendría que
    decidir cuándo parar ante una interrupción y **queda fuera de este
    trabajo**.
  - Depende de: F6.4.1.
  - Procede de: REQ-CONF-60, REQ-MAN-29, DTC-089; decisiones **S-01**,
    **S-02**, **S-04**, **S-05** y **S-26**; contradicción **C-5**, que esta
    forma resuelve.
  - Cambia: `schemas/config-general-v4.schema.json` y
    `schemas/effective-config-v6.schema.json` nuevos; `src/configuration.rs`,
    `src/adapters.rs`, `src/results.rs`. **Ningún esquema v1 a v5 cambia, ni
    `truth-v1`, ni `truth-v2`.**
  - **Declaración compacta, almacenamiento exacto.** Es la razón de que se
    separen las dos cosas: escribir la cadena de bases de cada alelo de cada
    individuo reventaría el techo de 65536 bytes del fichero general, mientras
    que un puñado de definiciones de alelo por locus más sus frecuencias cabe
    de sobra y **no crece con el número de individuos**. El programa
    materializa la secuencia desde la definición.
  - Se comprueba con: pareja mínima válida admitida; clave desconocida
    rechazada (`additionalProperties: false`); identificador de locus repetido
    rechazado; frecuencias que no suman 1 dentro de la tolerancia de DTC-049
    rechazadas nombrando el locus; una frecuencia negativa rechazada; un alelo
    sin frecuencia y una frecuencia sin alelo rechazados; **un escenario con
    500 individuos y 5 loci cabe holgadamente en el techo de 65536 bytes**, que
    es la comprobación de que C-5 está resuelta; un documento legítimo de 65536
    bytes admitido y el de 65537 rechazado; validación externa con
    `python3 -m jsonschema`; comparación byte a byte de los esquemas heredados;
    y **un fichero general v2 sin loci sigue publicando `effective-config-v4`,
    no v6**. **Debe fallar contra:** un lector que promueva todo documento a
    v6; un esquema heredado con una propiedad añadida; una validación de
    frecuencias que acepte cualquier reparto.

- ☐ **F6.4.2.** **Guarda de admisión con referencia real:** al declarar un
  locus sobre una referencia suministrada o reducida, comprobar que en esa
  posición **hay realmente un tracto de ese motivo**, y rechazar la
  configuración identificando el locus si no lo hay.
  - Depende de: F6.1.1.1, F1.7.2.
  - Procede de: decisión **S-06**, que introduce expresamente esta guarda:
    «sin ella se anotaría como microsatélite un trozo cualquiera del genoma».
  - Cambia: `src/reference.rs`, `src/tandem_repeat.rs`.
  - **El criterio, que la decisión no fijaba y sin el cual no hay guarda.** Una
    guarda que exija el motivo **perfecto** rechazaría tractos reales
    legítimos, que casi siempre son impuros —era el defecto de la primera
    versión de este plan, que rechazaba «una base cambiada dentro del tracto»—.
    Una que no exija nada no es guarda. El criterio es por tanto un **mínimo de
    pureza configurable**, con la pureza medida como en F6.4.1, más **un mínimo
    de unidades completas en la tirada más larga**. Ambos tienen valor por
    omisión, se publican en la configuración efectiva y **son lo que hace que
    exista una variante defectuosa posible**.
  - Se comprueba con: un tracto puro se admite; un tracto real impuro por
    encima del mínimo de pureza se admite —este es el caso que la primera
    versión rechazaba por error—; una región **sin ninguna repetición del
    motivo** se rechaza nombrando el locus, la pureza medida y el mínimo; un
    tracto justo por debajo del mínimo se rechaza y justo por encima se admite,
    lo que demuestra que el umbral se lee y no está codificado; bajar el mínimo
    en la configuración convierte el rechazo en admisión. **Debe fallar
    contra:** una guarda que exija pureza 1, que rechazaría el caso real
    legítimo; una guarda ausente, que admitiría la región sin repeticiones; un
    mínimo codificado en el programa en vez de leído.

- ☐ **F6.4.3.** El alelo por copia haplotípica: **clase nueva del mundo de las
  variantes pequeñas**, junto al cambio de base, la inserción y la deleción,
  que lleva el **intervalo de referencia** que sustituye, la **secuencia
  realizada** del tracto y su anotación —motivo, unidades, interrupciones—
  derivada de esa secuencia.
  - Depende de: F6.4.2.
  - Procede de: REQ-IND-32/33, REQ-CONF-58/59, ACE-103, DTC-061; decisión
    **S-02**. Corrección de dependencia: no depende del hito de variantes
    estructurales.
  - Cambia: `src/segments.rs`, `src/materialization.rs`,
    `src/small_event_composition.rs`.
  - **Lo que lo hace barato.** `InsertionModification` ya guarda secuencia en
    un `Arc<[u8]>`, de modo que la estructura no estrena capacidad. Y el trío
    intervalo de referencia + secuencia alternativa + anotación **es el modelo
    del propio VCF** para una repetición en tándem, así que la proyección de
    F6.5.1.1 es directa y no una traducción.
  - **Frontera de alcance.** Este trabajo **varía tractos que ya existen y
    están declarados**. La aparición de microsatélites **nuevos** por inserción
    o duplicación es trabajo estructural futuro y queda fuera.
  - Se comprueba con: las dos copias de un diploide con alelos de 5 y 7
    unidades materializan secuencias distintas; **dos alelos con el mismo
    número de unidades y la interrupción en sitio distinto materializan
    secuencias distintas y se distinguen en el estado**, que es lo que un
    contador no podía expresar; dos modificaciones sobre el mismo tracto se
    rechazan por solapamiento, como ya hace
    `SegmentError::OverlappingModifications`; un lector independiente
    reconstruye la secuencia y, de ella, el número de unidades y la pureza.
    **Debe fallar contra:** una materialización que use la secuencia de la
    referencia para las dos copias; **una representación que guarde solo el
    número de unidades**, que colapsaría los dos alelos de igual longitud y
    distinta interrupción en el mismo objeto —es la variante defectuosa que
    sustituye a la de «perder el número de unidades», que con la secuencia
    guardada ya no tiene sentido porque no hay número que perder—; una
    anotación que se almacene aparte y pueda contradecir la secuencia, que se
    caza mutando la anotación y exigiendo que el modelo la recalcule.

- ☐ **F6.4.3.1.** **Tope propio configurable** del tamaño del tracto, que
  sustituye al umbral de 50 pb: rechazo antes de simular de todo alelo cuya
  longitud, o cuya diferencia con la referencia, rebase el tope declarado,
  identificando **el locus y el alelo**.
  - Depende de: F6.4.3.
  - Procede de: decisión **S-05**, que retira expresamente el umbral de 50 pb
    para las repeticiones y lo sustituye por esta salvaguarda «para que un
    error de dedo no genere un tracto desmesurado».
  - Cambia: `src/configuration.rs`, `src/tandem_repeat.rs`,
    `schemas/config-general-v4.schema.json`.
  - **Por qué no el umbral de 50 pb.** Un alelo `(CA)40` frente a uno `(CA)10`
    difiere en 60 bases y **no es una variante estructural en ningún sentido**:
    es un polimorfismo de longitud de un microsatélite. El umbral existe para
    separar inserciones de secuencia **nueva** de reordenamientos; aplicarlo a
    un tracto repetido es un error de categoría. Y aprieta más de lo que su
    ficha admitía: muerde a partir de 25 unidades de diferencia en un
    dinucleótido, cifra que una población diversa alcanza.
  - **Por qué se identifica el locus y el alelo, y no el individuo.** Con
    alelos declarados a mano el rechazo podía nombrar individuo y copia. Con
    **frecuencias** (S-04) ya no: el alelo que rebasa el tope se conoce **al
    leer la configuración**, antes de que exista ningún individuo al que
    asignárselo. Nombrar individuo y copia aquí era un resto del modelo
    anterior.
  - **Su contrato, que hay que fijar aquí porque no existe en ningún sitio:**
    una clave propia en la tabla del locus, llamada explícitamente en
    **unidades** —`max_units`—, de dominio entero y positivo, publicada en la
    configuración efectiva. **Valor por omisión propuesto: 200 unidades**,
    holgado para cualquier microsatélite real y estrecho para un error de dedo;
    es el pendiente **P-5** y está a confirmar. El comportamiento no depende de
    ese número: los fixtures configuran el tope explícitamente.
  - Se comprueba con: con el tope en 60 unidades, un alelo de 60 se admite y
    uno de 61 se rechaza **nombrando el locus y el alelo**; **con el tope en 20
    unidades, el alelo de 60 pasa a rechazarse**, lo que demuestra que se lee
    el valor configurado y no una constante; una clave `max_bases` se rechaza
    **como campo desconocido**, porque el esquema es cerrado; un alelo de 30
    unidades de un dinucleótido —60 pares de base, que el umbral retirado
    habría rechazado— **se admite**. **Debe fallar contra:** un tope codificado
    en el programa —la prueba del tope a 20 lo caza—; la supervivencia del
    umbral de 50 pares de base en cualquiera de los tres momentos en que
    REQ-IND-31 lo exigía, que la prueba del dinucleótido de 60 pares de base
    caza; un esquema abierto, que la clave `max_bases` caza.
  - *Corregido:* la versión anterior pedía que «un tope expresado en pares de
    base se rechace nombrando el parámetro». **No es realizable**: un entero no
    lleva unidad consigo, de modo que el programa no puede distinguir 60
    unidades de 60 pares de base. Lo único comprobable es el **nombre de la
    clave**, y por eso la prueba es la del campo desconocido.

- ☐ **F6.4.3.2.** **Sorteo de los genotipos desde las frecuencias del locus**,
  con su ámbito de azar propio: cada copia haplotípica de cada individuo recibe
  un alelo del locus sorteado de las frecuencias declaradas.
  - Depende de: F6.4.3.1.
  - Procede de: decisión **S-04**, que invierte la recomendación original;
    REQ-DAT-15, DTC-037, DTC-072 y DTC-081 para la disciplina de aleatoriedad.
  - Cambia: `src/randomness.rs` (ámbito nuevo), `src/initial_genotypes.rs` o
    módulo hermano.
  - **Esto consume azar, y la primera versión de este plan no lo cubría.** Con
    alelos declarados a mano no se consumía ni un número aleatorio; con
    frecuencias sí. Hace falta por tanto: un **ámbito de semilla propio**,
    identificado por el evento y las identidades estables —individuo, copia y
    **identificador de locus**, nunca la posición ni el orden físico—; y un
    **algoritmo normativo**, que es `categorical-f64-v1` (DTC-072), ya
    existente, sin primitiva nueva. Los pesos se normalizan con la suma
    compensada de DTC-049.
  - Se comprueba con: vectores congelados del sorteo para una terna
    individuo/copia/locus; misma semilla, mismos genotipos; **insertar un
    individuo al principio de la lista no cambia el genotipo de ninguno de los
    demás**, comparado por identificador y no por posición, y **reordenar la
    lista de individuos tampoco**; **añadir un locus en un cromosoma que
    ordena antes que los existentes no cambia los genotipos de los loci
    previos**; con una frecuencia de 1 para un alelo, todos los individuos lo
    reciben y **no se consume aleatoriedad**, contado con un generador
    instrumentado; con dos alelos al **0,9 y 0,1**, sobre mil copias y semilla
    congelada, la proporción observada cae dentro del intervalo binomial
    declarado alrededor de 0,9.
  - **Debe fallar contra:** un ámbito indexado por el **orden** de los
    individuos o por la **posición** del locus —las pruebas de insertar al
    principio, reordenar y añadir un locus que ordena antes son las que lo
    cazan—; un sorteo que consuma con frecuencia 1; **un reparto uniforme que
    ignore los pesos**, que la proporción de 0,9 caza.
  - *Corregido, dos veces en el mismo paso.* La versión anterior usaba **dos
    alelos al 50 %**, y ese fixture es **degenerado**: una elección uniforme
    que ignore los pesos produce exactamente la misma distribución que la
    ponderada, de modo que el mutante sobrevivía. Y pedía **«añadir un
    individuo al final»**, que **no caza** un ámbito indexado por orden, porque
    añadir al final deja intactos todos los ordinales anteriores; hay que
    **insertar al principio o reordenar**. Lo mismo con el locus: añadirlo en
    un cromosoma que ordene **después** no desplaza nada. Va en la docstring.

- ☐ **F6.4.3.3.** Verdad publicada, versión **v3**: el alelo se publica como un
  caso de inserción o deleción con metadatos de repetición —locus, motivo y
  número de unidades—, y como el alelo lleva su secuencia, el caso ya la
  contiene y los metadatos solo la anotan. `truth-v1` y `truth-v2` no cambian
  ni un byte.
  - Depende de: F6.4.3.2.
  - Procede de: REQ-VER-22, REQ-IND-33; decisión **S-08**. Cierra la decisión
    pendiente DP-F6-03 en su parte de repeticiones.
  - Cambia: `schemas/truth-v3.schema.json` nuevo, `src/truth_schema.rs`,
    `src/truth_projection.rs`, `src/truth_validation.rs`.
  - **Cuándo se promueve la versión, que hacía falta fijar.** La verdad pasa a
    la versión 3 **cuando se publica al menos un alelo de repetición**, no por
    el mero hecho de que la configuración declare loci. Un escenario que
    declare un locus cuyo único alelo es el de la referencia no produce ninguna
    variante que publicar, y su verdad sigue siendo la versión 2, byte a byte
    igual a la de antes. Es lo que hace posible la comprobación (c) de F4.8.6,
    y es coherente con la regla de que un fichero general v2 siga publicando
    configuración efectiva v4.
  - Se comprueba con: el registro de un diploide con dos alelos distintos lleva
    los dos con su fase, su locus, su motivo y su número de unidades, y la
    secuencia de cada uno; **dos alelos de igual longitud y distinta
    interrupción producen registros distinguibles**; el validador de verdad
    existente sigue aceptando un documento v2 sin cambios; `sha256sum` de la
    verdad de una ejecución **sin** loci, idéntico al de la línea base; los dos
    esquemas heredados comparados byte a byte; **un escenario con un locus
    declarado cuyo único alelo es el de la referencia publica verdad en la
    versión 2 y byte a byte igual a la línea base**. **Debe fallar contra:**
    una proyección que publique solo el número de unidades y no la secuencia,
    que la prueba de los dos alelos de igual longitud caza; un cambio que
    promueva a v3 las ejecuciones sin repeticiones; **una promoción disparada
    por la declaración del locus en vez de por la publicación de un alelo**,
    que la prueba del alelo de referencia único caza y que dejaría sin base la
    comprobación (c) de F4.8.6.

- ☐ **F6.5.1.1.** Proyección al VCF: campos de repetición en tándem **del
  estándar `VCFv4.5`** declarados en la cabecera, y alelos **siempre
  explícitos**, con la secuencia completa en `REF` y `ALT`. **Sin ningún umbral
  de tamaño** y sin entrada complementaria.
  - Depende de: F6.4.3.3.
  - Procede de: REQ-VER-20/22, ACE-31, DTC-064, DTC-066; decisión **S-09**;
    confirmación D-10. Cierra la decisión pendiente DP-F6-06 y mantiene fuera
    del alcance el hito de variantes complejas y complemento del VCF.
  - Cambia: `src/vcf_contract.rs`, `src/vcf_projection.rs`,
    `src/vcf_representation.rs`, `src/vcf_semantic_validation.rs`.
  - **Por qué no hay umbral.** El criterio del estándar no es la longitud sino
    **si se conoce la secuencia exacta**: la forma explícita cuando se conoce,
    la simbólica cuando se está resumiendo. **Un simulador siempre la conoce,
    porque la ha generado él**, de modo que la forma explícita aplica siempre.
    El umbral de 50 pb era una convención de variantes estructurales importada
    de fuera del estándar.
  - **Tipos, ya verificados contra la especificación**, con dos precisiones.
    Primera: el campo que expresa el número de copias es de tipo **decimal
    (`Float`)**, no entero, porque un tracto puede terminar en una unidad
    parcial —lo que F6.4.1 ahora admite—. Declararlo entero produciría un VCF
    que una herramienta conforme rechaza. Segunda, **corregida el 23/09/2026**:
    la regla de omisión **no es del número de copias sino del número de
    secuencias repetidas**. El estándar dice, de ese campo y no del otro, que
    si se omite se asume uno para los alelos simbólicos de repetición y **cero
    para los demás**, de modo que omitirlo no es «no decir nada» sino
    **afirmar cero**. El proyecto lo emite **explícito para cada alelo
    alternativo**; el número de copias no tiene regla de omisión propia.

> **Retirado el 22/09/2026: la exigencia sobre el alelo de referencia.** Este
> paso pedía que el número de copias se emitiera también para el **alelo de
> referencia**. **No es expresable** con un locus de varios motivos: el número
> de copias por bloque del tracto que hay en la referencia exigiría una
> estructura en bloques, y la decisión S-30 declara que esa estructura **se
> declara y no puede derivarse**. El tracto de la referencia no la declara en
> ninguna parte. Queda por tanto **solo para los alelos alternativos**, que son
> los que el usuario declara con su estructura. Falta además confirmar contra
> la especificación si el campo que cuenta las secuencias repetidas es «uno por
> alelo alternativo», en cuyo caso el alelo de referencia sencillamente no se
> describe con estos campos y la exigencia sobraba por partida doble.
  - Se comprueba con: un diploide con dos alelos produce **un** registro con
    dos alternativas y genotipo faseado, alineado a la izquierda del tracto
    como DTC-064 obliga; la cabecera declara cada campo usado, con su número y
    su tipo, y ninguno más; el campo de número de copias se emite como decimal
    y un tracto con unidad parcial produce un valor no entero; **cada alelo
    alternativo lleva su número de copias explícito** —el de referencia **no**,
    por lo dicho en el recuadro de arriba—, y el oráculo comprueba que cada
    valor **se corresponde con la longitud del alelo al que acompaña**, en su
    mismo orden; el oráculo
    semántico propio reconstruye la secuencia del alelo desde `REF`/`ALT` y la
    compara con la verdad; **ningún registro usa representación simbólica**,
    comprobado con una búsqueda negativa con control positivo, ni siquiera con
    un tracto de diez mil bases. **Debe fallar contra:** una cabecera que
    declare el número de copias como entero; una proyección que recaiga en la
    forma simbólica por encima de algún tamaño, que la prueba del tracto de
    diez mil bases caza; dos registros para el mismo sitio en vez de uno
    multialélico; **un campo de número de copias omitido para algún alelo**,
    que el estándar interpretaría como cero y que la comprobación de
    correspondencia por orden caza; una lista de valores desalineada respecto
    de la de alelos, que es el error silencioso más probable de este campo.
  - **Pendiente P-2, dicho sin adornos:** la puerta externa exige `bcftools` y
    `tabix` 1.24, y el entorno tiene 1.21 y **1.19**. Este paso ejecuta la
    comprobación con lo instalado, **registra las dos versiones en la
    evidencia** y **no etiqueta la puerta como superada**.

- ☐ **F6.4.8.** Registro en configuración efectiva y manifiesto: loci y sus
  parámetros, procedencia de la referencia, motivo, alelos con su secuencia y
  sus frecuencias declaradas, **frecuencias realizadas** tras el sorteo, tope
  configurado, mínimos de la guarda de admisión, y los **valores realizados**
  de pureza e interrupciones por alelo.
  - Depende de: F6.5.1.1.
  - Procede de: REQ-MAN-29, REQ-VAL-48; decisiones **S-01**, **S-04**,
    **S-05**, **S-06**.
  - Cambia: `src/results.rs`, `schemas/manifest-general-v6.schema.json`.
    **Corregido el 23/09/2026:** la primera redacción nombraba también
    `src/event_manifest.rs`, que describe un **evento técnico de
    secuenciación** —bibliotecas, cobertura, quimeras—. Los loci de repetición
    pertenecen a la referencia y a la población inicial, y una ejecución que
    los declare puede no publicar ningún evento técnico: atar su registro a un
    documento que a menudo no existe lo dejaría sin publicar. No se toca.
  - Se comprueba con: validación externa de los dos documentos; la pureza
    publicada por alelo coincide con la calculada a mano sobre su secuencia;
    **la frecuencia realizada difiere de la declarada** en un corpus pequeño y
    converge a ella en uno grande, que es lo que distingue publicar lo realizado
    de copiar lo declarado; el tope y los mínimos publicados son los
    configurados, no los de omisión, cuando la configuración los sustituye.
    **Debe fallar contra:** un manifiesto que publique la frecuencia declarada
    como realizada —la prueba del corpus pequeño lo caza, y con un corpus
    grande **no discrimina**, lo que va en la docstring—; un manifiesto que
    publique los valores de omisión cuando la configuración los sustituye.

- ☐ **F6.4.8.1.** **Ajustar el presupuesto de tamaño de salida** para que
  cuente las bases del tracto que la verdad y el VCF ahora publican.
  - Depende de: F6.4.8.
  - Procede de: DTC-089 y la consecuencia que la decisión **S-02** dejó
    anotada: «el presupuesto de tamaño de salida debe contar las bases del
    tracto en la verdad y en el VCF; es un término pequeño, pero el contrato de
    recursos es estricto y no puede descubrirse a mitad».
  - Cambia: el modelo de coste de salida donde viva, junto a las constantes de
    DTC-089.
  - Se comprueba con: un escenario con loci declarados produce una cota de
    bytes de salida **estrictamente mayor** que el mismo sin ellos, y la
    diferencia es al menos las bases de los alelos publicados; la cota sigue
    siendo **superior o igual** a los bytes realmente escritos, comprobado
    sobre un escenario con alelos largos; un escenario sin loci da exactamente
    la cota de antes. **Debe fallar contra:** una cota que ignore el término
    nuevo —el escenario con alelos largos la caza, porque los bytes escritos la
    rebasarían—; una que lo cuente también cuando no hay loci, que la
    comparación con la línea base caza.

- ☐ **F6.4.9.1.** Interfaz Marimo: declarar loci **con motivo y posición**,
  declarar alelos **con sus interrupciones**, y declarar **frecuencias** por
  locus, con equivalencia interfaz → TOML → configuración efectiva.
  - Depende de: F6.4.8.1.
  - Procede de: REQ-UI-01/03/04/05, hoja de ruta §7; decisión **S-17**, que
    advierte de que «mínima» hoy significa **el doble** de lo que significaba
    cuando se escribió la ficha; contradicción **C-6**.
  - Cambia: `ui/config_editor.py`, `ui/ui_panels.py`, `ui/gps_app.py`; prueba
    nueva en `tests_ui/`.
  - Se comprueba con: ida y vuelta de la interfaz al TOML y del TOML a la
    configuración efectiva; las interrupciones declaradas sobreviven la ida y
    la vuelta; las frecuencias se normalizan y se rechaza un reparto inválido
    **en la interfaz**, no solo en el lector; ausencia de controles para lo que
    el corte no admite —pureza objetivo, generación aleatoria de loci—;
    `marimo check --strict` sin avisos. **Debe fallar contra:** un control de
    pureza objetivo presente; una ida y vuelta que pierda las interrupciones,
    que es justo lo que se perdería si la interfaz siguiera pensando en número
    de unidades.

> **Revisión del hito de la entrega 1.** Comprobar la referencia con el tracto
> y su mapa de coordenadas, el modelo del locus con la secuencia, la operación
> sobre la tirada más larga, la guarda de admisión con su criterio, el sorteo
> por frecuencias y su ámbito de azar, la verdad v3, el VCF sin umbral, el
> presupuesto de salida y la interfaz. Ejecutar la batería completa con
> `./scripts/check-independent-validators.sh evidence/repeats-f6-4`. Confirmar
> que los artefactos de una ejecución **sin** loci coinciden byte a byte con la
> línea base de F6.4.0. Aplicar la parte de la sección 10 que corresponde a la
> fase del núcleo reproducible. No continuar hasta confirmar el resultado.

---

## 6. Entrega 2 — Tartamudeo de base, sin amplificar (F4.8)

Resultado: una molécula que contiene un tracto declarado puede mostrar un
número de unidades distinto del que el individuo tiene, por deslizamiento de la
polimerasa durante la amplificación en puente sobre la celda de flujo y durante
la síntesis de la secuenciación. **La verdad del individuo no cambia.**

**Cómo se nombra, y por qué no contradice el contrato de la fase.** Una
biblioteca «sin PCR» evita la amplificación en la preparación, pero no evita
toda polimerasa: la formación de agrupaciones en la celda es amplificación en
puente, y la secuenciación por síntesis es copia con polimerasa. El requisito
nuevo debe decir que **este artefacto no es una preparación amplificada**, que
no introduce duplicados y que el perfil sigue registrándose como sin PCR
(decisión **S-12**).

- ☐ **F4.8.1.** Decisión técnica nueva y catálogo de ámbitos aleatorios
  **`f4-technical-v4`**, que extiende el vigente con la etapa `stutter` y los
  componentes `[event_id, library_id, stutter, template_ordinal, locus_id,
  decision]`, indexada por **identificador de locus** y nunca por la posición
  del tracto dentro de la molécula.
  - Depende de: cierre de la entrega 1.
  - Procede de: REQ-DAT-15, DTC-037, DTC-081, DTC-086, DTC-088; decisiones
    **S-11** y **S-12**. Fila nueva propuesta: **DTC-090** (provisional).
  - Cambia: `src/randomness.rs`.
  - **Aclaración que la decisión debe llevar por escrito:** DTC-067 obliga a
    que los **recuentos** aleatorios salgan de `poisson-count-v1`. Aquí no hay
    recuento de sucesos sobre oportunidades: hay un ensayo binario y una
    elección entre alternativas, y los consumidores que corresponden ya existen
    —`bernoulli-v1` y `uniform-index-v1` (DTC-088), `categorical-f64-v1`
    (DTC-072)—. Decirlo evita que alguien lea la regla al revés.
  - Se comprueba con: vectores congelados de derivación de subsemilla para una
    terna evento/biblioteca/locus; el catálogo vigente sigue produciendo
    exactamente las mismas subsemillas que antes, comparadas con la línea base.
    **Debe fallar contra:** una extensión que renumere o reordene las etapas
    existentes.

- ☐ **F4.8.1.1.** **Guarda de admisión: un locus que vaya a tartamudear exige
  una tirada mínima efectiva de al menos 2 unidades.**
  - Depende de: F4.8.1.
  - Procede de: **decisión del responsable del 21/09/2026**, tomada sobre un
    hallazgo de la segunda revisión adversarial de la entrega 1.
  - Cambia: la admisión que cruza el evento técnico con los loci declarados.
    **No** `src/tandem_repeat.rs`: el modelo no se toca.

  **Qué exige.** Cuando un locus vaya a sufrir tartamudeo, su mínimo de tirada
  **efectivo** —el suyo si lo declara, y si no el general— debe ser **2 o más**.
  Si es 1, la configuración se rechaza nombrando el locus, el valor declarado y
  el mínimo exigido.

  **Por qué.** El paso de unidades no permite que una contracción vacíe su
  tirada (ver la nota de `apply_unit_step` en `src/tandem_repeat.rs`). La
  consecuencia es que **un tracto cuya tirada más larga sea de una sola copia no
  puede contraerse nunca, mientras que expandirse sí puede**: una clase entera de
  tractos queda congelada frente a la contracción con la expansión libre, lo que
  **sesga el modelo hacia el alargamiento**. Hasta la decisión S-31 esa clase no
  llegaba a entrar, porque la tirada mínima general es 3; desde que cada locus
  puede declarar la suya, uno con tirada mínima 1 la admite. **La decisión S-31
  abrió esa puerta sin querer y esta guarda la cierra.**

  El responsable eligió esta salida en vez de tocar el modelo, y conviene
  registrar por qué:

  - **Lo que cierra:** desaparece la clase de microsatélites que no pueden
    contraerse nunca. Todo alelo admitido podrá perder al menos una unidad.
  - **Lo que no cierra, y está bien:** con una tirada de dos copias siguen
    prohibidos los pasos de dos y tres unidades, porque la vaciarían. Esa
    asimetría residual **es físicamente correcta** —en la realidad tampoco se
    pueden perder tres unidades de una tirada que solo tiene dos— y además los
    tramos cortos apenas tartamudean, porque la tasa crece con la longitud. La
    opción convierte un **sesgo artificial del modelo** en una **restricción
    física real**.
  - **A quién afecta:** a casi nadie. El valor por omisión es 3. Solo toca a
    quien baje deliberadamente el mínimo de su locus a 1, cosa que el esquema
    hoy permite.

  **Dónde hay que comprobarlo, y por qué no se pudo hacer antes.** El
  microsatélite se declara **junto a la referencia** y el tartamudeo se
  configura **en el evento técnico**: son dos sitios distintos del fichero. La
  condición «si va a tartamudear, exige 2» **no se puede comprobar al leer el
  locus**, porque en ese momento todavía no se sabe si habrá tartamudeo. Hay que
  comprobarla **cuando las dos configuraciones están resueltas**. Por eso la
  guarda no se escribió durante la entrega 1 pese a haberse decidido allí: una
  guarda escrita entonces habría quedado sin efecto.

  - Se comprueba con: un locus con tirada mínima 1 **y** tartamudeo activo se
    rechaza nombrando el locus y los dos valores; el mismo locus **sin**
    tartamudeo se admite, que es lo que demuestra que la guarda mira las dos
    configuraciones y no solo el locus; con tirada mínima 2 y tartamudeo activo
    se admite; y un locus que no declara mínimo propio hereda el general de 3 y
    entra sin tocar nada.
  - **Debe fallar contra:** una guarda que mire solo el locus y rechace la
    tirada mínima 1 aunque no haya tartamudeo —la segunda comprobación la caza—;
    una que no mire el mínimo **efectivo** sino solo el declarado, que dejaría
    pasar un locus sin mínimo propio bajo un general de 1; y la ausencia de
    guarda, que la primera comprobación caza.

- ☐ **F4.8.2.** Modelo `stutter-v1`, puro, con la **forma publicada**: la
  probabilidad de resbalón por ciclo depende **linealmente de la longitud de la
  tirada ininterrumpida más larga** del motivo contenida en la molécula, y el
  tamaño del paso sigue una estructura **asimétrica** —contracciones de hasta
  tres unidades, expansiones de una sola, con sesgo hacia la contracción—.
  **Sin reintentos**: los pasos inviables se excluyen **antes** de sortear.
  - Depende de: F4.8.1, F6.4.1.
  - Procede de: decisión **S-14**; Raz y col. (2019), NAR 47(5):2436–2445,
    CC BY 4.0, que infieren el modelo midiendo el ruido *in vitro* y fijan la
    dependencia lineal **por ciclo**; DTC-072, DTC-088; pendientes **P-1** y
    **P-3**.
  - Cambia: módulo nuevo `src/stutter.rs`.
  - **El predictor es la tirada más larga, no el número total de unidades.** Es
    lo que cambió respecto de la primera versión de este plan: en genética
    forense está establecido que la tirada ininterrumpida más larga predice el
    tartamudeo mejor que la longitud total del alelo, y es coherente con que
    las interrupciones **estabilicen** el tracto. El valor lo calcula F6.4.1 y
    este paso lo consume; no se recalcula aquí.
  - **La forma lineal no es una aproximación tosca, y conviene dejarlo escrito
    porque parece serlo.** La literatura describe la fracción observada de
    extremo a extremo con una curva que satura. Son **dos magnitudes
    distintas**: la dependencia lineal es de la probabilidad **por ciclo**, y
    la curva con meseta sale de acumular esa probabilidad sobre muchos ciclos.
    Ambas son ciertas y no se contradicen.
  - Se comprueba con:
    **(a) el predictor.** Dos tractos con el **mismo número total de unidades**
    y **distinta tirada más larga** —por ejemplo diez unidades en una sola
    tirada frente a diez unidades partidas en tiradas de seis y cuatro— dan
    probabilidades **distintas**, y mayor la del primero, cuando el coeficiente
    es positivo.
    **(b) la longitud del motivo.** Un mononucleótido y un tetranucleótido con
    la **misma longitud en pares de base** dan mayor probabilidad al
    mononucleótido.
    **(c) la asimetría.** El conjunto de pasos admitidos es exactamente
    `{−3, −2, −1, +1}`; una configuración que declare `+2` se rechaza.
    **(d) el filtrado previo.** Un generador instrumentado confirma que,
    **cuando solo queda un paso viable**, el sorteo categórico **no consume**
    aleatoriedad. El fixture lo consigue con un tracto de una sola unidad
    completa, donde las tres contracciones dejarían el tracto sin ninguna y
    solo la expansión de una unidad sobrevive al filtrado.
    **Debe fallar contra:** un modelo cuyo predictor sea el **número total de
    unidades** —prueba (a), que es la única que lo caza—; un modelo que dependa
    solo de los pares de base y no de la longitud del motivo —prueba (b)—; una
    estructura simétrica de pasos —prueba (c)—; una implementación que sortee
    entre los cuatro pasos y **descarte** los inviables en vez de filtrarlos
    antes —prueba (d), que cuenta el consumo, no el resultado.
  - *Corregido:* la versión anterior pedía comprobar que no se consume
    aleatoriedad «cuando ningún paso deja el tracto viable». **Ese caso no
    existe y no tiene fixture posible**: una expansión de una unidad es
    **siempre** viable, porque ni vacía el tracto ni acorta la molécula. La
    aserción no podía fallar nunca. La comprobación que sí discrimina es la del
    candidato único.
  - *Cuidado, y va en la docstring de la prueba (a):* **con tractos puros el
    número total de unidades y la tirada más larga coinciden**, de modo que un
    fixture puro **no discrimina** entre los dos predictores. El fixture tiene
    que ser impuro y con el total igualado a propósito. Era el defecto exacto
    de la primera versión de este plan.

- ☐ **F4.8.3.** Tasas por longitud de motivo tomadas de la **columna sin PCR**
  de la calibración, que con el coeficiente de la tirada en cero son
  directamente la probabilidad por longitud de motivo. La procedencia —humano,
  Illumina, con su fuente y su licencia— viaja con los valores en la
  configuración efectiva y en el manifiesto.
  - Depende de: F4.8.2.
  - Procede de: REQ-CONF-131/132 propuestos; contradicciones **C-1** y **C-2**;
    pendientes **P-1** y **P-3**.
  - Cambia: `src/stutter.rs`, `src/technical_profile.rs`, los dos esquemas
    nuevos.
  - **Valores de partida** (fracción de lecturas afectadas, columna sin PCR):
    motivo de 1 pb, 0,027; 2 pb, 0,011; 3 pb, 0,0029; 4 pb, 0,0016; 5 pb,
    0,0015; 6 pb, 0,0010, **este último marcado como el menos fiable**, porque
    el propio autor advierte que procede del ruido de analizar pocos loci.
  - **Qué se configura y qué no.** La **estructura** del paso no se configura:
    la fija el algoritmo (F4.8.2). Los **pesos relativos** dentro de esa
    estructura son **obligatorios** cuando el tartamudeo está activo, sin valor
    por omisión, porque la fuente publica el método de ajuste pero no los
    valores. El **coeficiente de la tirada** arranca en cero y se publica como
    no calibrado (P-1).
  - Se comprueba con: la configuración efectiva publica los seis valores, la
    cita de la fuente, la etiqueta de especie y tecnología y la marca de no
    calibrado del coeficiente; con el coeficiente en cero, la probabilidad
    calculada **reproduce exactamente** la tabla por longitud de motivo; una
    configuración con tartamudeo activo y **sin** pesos declarados se rechaza
    nombrando el parámetro; un peso negativo, una tasa fuera de `[0,1]` y un
    paso expresado en bases se rechazan antes de simular. **Debe fallar
    contra:** unos valores embebidos en el programa sin la etiqueta de
    procedencia; una configuración sin pesos aceptada con un reparto implícito;
    un coeficiente distinto de cero por omisión, que afirmaría una calibración
    que no existe.

- ☐ **F4.8.4.** `stutter::transform_paired_end`, con la misma forma que los dos
  transformadores que ya existen, aplicado **después de la quimera y antes de
  calidad y errores**. Una única decisión por molécula y tracto, aplicada a la
  molécula completa —reutilizando la operación de F6.4.1.1, es decir, **en la
  tirada más larga**— antes de volver a emitir R1 y R2.
  - Depende de: F4.8.3.
  - Procede de: REQ-DAT-74/75/76 propuestos; decisiones **S-07** y **S-13**.
  - Cambia: `src/stutter.rs`, `src/execution.rs`, `src/wgs_sampling.rs` si el
    punto de inserción lo exige.
  - **Por qué ahí y no antes.** `src/chimera.rs` vuelve a materializar el
    prefijo del fragmento desde la copia haplotípica y lo compara **byte a
    byte**; si no coincide devuelve `ChimeraError::PairSourceMismatch`. Un
    fragmento ya modificado rompería esa comprobación. Después de la quimera la
    pareja todavía lleva su fragmento completo, de modo que la molécula sigue
    disponible.
  - Se comprueba con:
    **(a) coherencia del par**, sobre un fixture **saturado**: copia
    haplotípica lo bastante corta para forzar el inserto mínimo del perfil, 250
    pb, de modo que R1 y R2 —de 150 pb— compartan las 50 bases centrales, y un
    tracto de 30 pb **íntegramente dentro de ese solapamiento**; se cuentan las
    unidades en R1 y en el complemento inverso de R2 y se exige que coincidan.
    **(b) quimera**: un corpus con tasa de quimera y de tartamudeo altas a la
    vez no produce ni un `ChimeraError::PairSourceMismatch`.
    **(c) dónde cae el paso**: sobre un tracto impuro con tiradas de 3 y 5
    unidades **y la tirada de 3 a la derecha**, la molécula tartamudeada
    conserva esa tirada de 3 byte a byte. *La posición importa y va en la
    docstring:* si la tirada larga estuviera a la derecha, la convención de la
    tirada más larga y la del extremo derecho coincidirían y la prueba no
    discriminaría nada.
    **Debe fallar contra:** un tartamudeo aplicado **por lectura** (a); una
    implementación que modifique el fragmento **antes** de la quimera (b); una
    que aplique el paso por el extremo derecho en vez de en la tirada más larga
    (c). *Cuidado, en la docstring:* un fixture en el que solo R1 cubra el
    tracto **no discrimina** —la aserción se cumple sola— y con la media de
    inserto de 350 pb la mayoría de las parejas están en ese caso; y sin tasa
    de quimera alta, la prueba (b) no ejercita el camino que dice cubrir.

- ☐ **F4.8.5.** La molécula cambia de longitud en el número de unidades por la
  longitud del motivo, y con ella el inserto observado. No se compensa ni se
  descarta la molécula. Las dos lecturas conservan la longitud del perfil.
  - Depende de: F4.8.4.
  - Procede de: REQ-DAT-78 propuesto; decisión **S-15**.
  - Cambia: `src/stutter.rs`, `src/paired_end.rs`.
  - Se comprueba con: una molécula que pierde una unidad de un motivo de 4 pb
    produce un fragmento exactamente 4 pb más corto, observado **sobre la
    molécula**; las dos lecturas siguen midiendo 150 pb; una contracción que
    dejaría la molécula más corta que la lectura queda excluida por el filtrado
    de F4.8.2 y no produce ningún error duro. **Debe fallar contra:** una
    implementación que recorte o rellene para conservar la longitud sorteada.
    *Cuidado:* observar las lecturas **no discrimina**, porque siguen midiendo
    lo mismo pase lo que pase; la aserción tiene que mirar la molécula.

- ☐ **F4.8.6.** Compatibilidad con lo ya publicado, **acotada**: una ejecución
  de un escenario que **no declara loci de repetición** publica exactamente los
  mismos bytes que antes de que el tartamudeo existiera.
  - Depende de: F4.8.5.
  - Procede de: REQ-DAT-80 propuesto con la acotación de **C-1**, DTC-081;
    decisión **S-16**.
  - **Recordatorio:** declarar un locus activa el tartamudeo de base, y por eso
    la guarda de F4.8.1.1 —tirada mínima efectiva de al menos 2 cuando el locus
    vaya a tartamudear— tiene que estar en pie antes de que este paso pueda dar
    por buena la acotación.
  - Cambia: nada; es una prueba de sistema y su evidencia.
  - Se comprueba con:
    **(a)** `sha256sum` de los dos FASTQ, la verdad, el VCF, la configuración
    efectiva y los dos manifiestos de un escenario **sin** loci, comparados con
    la línea base de F6.4.0.
    **(b)** un generador instrumentado confirma **cero consumo** en la etapa
    `stutter` en ese escenario.
    **(c)** el mismo escenario **con un locus declarado cuyo único alelo es el
    de la referencia, a frecuencia 1**, produce **FASTQ distintos** y **la
    verdad y el VCF de verdad idénticos** a los de (a), byte a byte. *Los
    artefactos se enumeran uno a uno a propósito:* la **configuración
    efectiva** y el **manifiesto general** **tienen que diferir**, porque el
    fichero general pasa a la versión 4 al declarar el locus, la configuración
    efectiva sube en consecuencia a la versión 6 y REQ-MAN-29 obliga a
    registrar el locus. Exigir su identidad, como hacía la versión anterior de
    este paso, era pedir lo imposible. La verdad y el VCF sí pueden ser
    idénticos, y lo son gracias a la regla de promoción que fija F6.4.3.3:
    sin alelo de repetición publicado, la verdad no sube de versión.
    **(d)** dos ejecuciones idénticas salvo en la tasa de tartamudeo conservan
    **la verdad y el VCF de verdad**, salvo por el identificador de ejecución,
    que cambia por contrato al cambiar los bytes del documento. *Corregido el
    24/09/2026:* la redacción anterior decía «difieren **solo** en los FASTQ y
    en los recuentos agregados», y eso es falso — también difieren la
    configuración efectiva, que publica la tasa, y el manifiesto general, que
    resume las huellas. Lo que hay que enumerar es lo que **no** cambia.
    **Debe fallar contra:** una implementación que consuma aleatoriedad con la
    tasa a cero o que desplace el ordinal de plantilla (a y b); una que aplique
    el tartamudeo a la **copia haplotípica** en vez de a la molécula, es decir,
    que cambie el genotipo (c); una lectura de la congelación como absoluta,
    que dejaría el caso con loci en cero tartamudeo (c, que exige que los FASTQ
    **difieran**).
  - *Por qué el alelo de referencia a frecuencia 1, y no un locus polimórfico.*
    La primera versión de este plan pedía «el mismo escenario con un locus
    declarado» y exigía verdad idéntica. **Es imposible de cumplir** si el
    locus es polimórfico: entonces hay alelos que publicar y la verdad cambia
    legítimamente. Con el alelo de referencia a frecuencia 1 el tracto existe
    —y por tanto tartamudea— pero no hay ninguna variante que publicar, que es
    justo lo que hace la comparación posible. Va en la docstring.

- ☐ **F4.8.7.** Configuración efectiva y manifiesto del tartamudeo: activación,
  parámetros efectivos por longitud de motivo y por locus, versión del
  algoritmo, versión del catálogo de ámbitos y, como **recuentos agregados**,
  moléculas evaluadas, moléculas que tartamudearon y distribución de los pasos
  aplicados por locus. Nada por lectura ni por pareja.
  - Depende de: F4.8.6.
  - Procede de: REQ-MAN-93 propuesto, REQ-DAT-56; decisión **S-14**;
    pendientes **P-1** y **P-3**.
  - Cambia: `src/event_manifest.rs`, `schemas/manifest-event-v3.schema.json`.
  - Se comprueba con: validación externa del manifiesto; la suma de la
    distribución de pasos por locus es igual al número de moléculas que
    tartamudearon; el manifiesto marca los parámetros con su calibración de
    origen, el coeficiente de la tirada como no calibrado, y **el desajuste
    conocido** de que la fuente del modelo distingue el *tipo* de motivo
    mientras la calibración del proyecto va por longitud; ningún campo nombra
    una lectura o una pareja concreta, comprobado con una búsqueda negativa con
    control positivo. **Debe fallar contra:** un manifiesto que publique el
    número **esperado** de moléculas tartamudeadas en vez del realizado; una
    procedencia por pareja filtrada; la ausencia de la marca de no calibrado,
    que haría pasar por medido un coeficiente que nadie ha medido.

- ☐ **F4.8.8.** Validador independiente del tartamudeo: reproduce las
  decisiones desde la semilla **sin reutilizar el código productor**, y
  comprueba sobre los FASTQ finales que las dos lecturas de cada pareja que
  cubren un mismo tracto muestran el mismo número de unidades.
  - Depende de: F4.8.7.
  - Procede de: REQ-VAL-57 propuesto, arquitectura §11.
  - Cambia: validador nuevo bajo `tests/`.
  - Se comprueba con: el validador coincide con los recuentos agregados del
    manifiesto sobre el fixture saturado de F4.8.4;
    `./scripts/check-independent-validators.sh evidence/stutter-f4-8` da **5
    PASS**. **Debe fallar contra:** un validador que lea los bytes publicados
    por el productor en vez de recomponerlos —se detecta mutando un recuento
    del manifiesto y exigiendo que el validador lo rechace.

- ☐ **F4.8.9.** Interfaz Marimo del tartamudeo: activación, parámetros por
  longitud de motivo, pesos de los pasos y sustituciones por locus, con
  equivalencia interfaz → TOML → configuración efectiva.
  - Depende de: F4.8.8, F6.4.9.1.
  - Procede de: REQ-UI-01/03/04/05, hoja de ruta §7; decisión **S-17**;
    contradicción **C-6**.
  - Cambia: `ui/config_editor.py`, `ui/ui_panels.py`; prueba en `tests_ui/`.
  - Se comprueba con: ida y vuelta completa; no hay control para un paso `+2`
    ni para un modelo distinto del admitido; `marimo check --strict` sin
    avisos; y **la comprobación que C-6 motiva**: partiendo solo de la
    interfaz, se puede construir un escenario completo —referencia, locus,
    alelos, frecuencias y tasa— que produzca lecturas tartamudeadas. **Debe
    fallar contra:** una interfaz con parámetros de tartamudeo y sin forma de
    declarar un locus con su motivo, sus alelos y sus frecuencias.

> **Revisión del hito de la entrega 2.** Comprobar el predictor sobre el
> fixture de total igualado, el punto de inserción, la coherencia del par sobre
> el fixture saturado, la congelación acotada, la verdad intacta, el validador
> independiente y la interfaz. Ejecutar la batería completa con
> `./scripts/check-independent-validators.sh evidence/stutter-f4-8`. Aplicar el
> cambio de contrato de la sección 10 en su parte de tartamudeo y **repetir la
> evidencia de cierre de la fase de realismo técnico**. No continuar hasta
> confirmar el resultado.

---

## 7. Entrega 3 — Preparación amplificada por PCR (F4.9)

Resultado: un modo de preparación `pcr`, alternativo al único que hoy existe,
en el que el recorrido técnico deja de producir una pareja por fragmento y pasa
a producir **familias de amplificación**: grupos de parejas que descienden de
la misma molécula original. De ahí salen los duplicados, la
sobrerrepresentación de unas moléculas frente a otras y el sesgo por
composición. **Ningún paso de este hito cita un locus de repetición.**

**Cómo se acota lo inacotable.** No se simula la población de moléculas, que
crece como `2^ciclos`. Se muestrea **hacia atrás** desde las parejas que se van
a publicar: como el número de parejas `N` de una biblioteca ya es determinista
desde la cobertura, se generan familias hasta cubrirlas, y cada familia cuesta
un número de operaciones acotado por constantes de contrato.

**Lo que cambia respecto de la primera versión de este plan:** el usuario ya no
declara «parejas por molécula», que es la variable del modelo, sino la
**fracción de duplicados**, que es la que un biólogo ve en su informe de
calidad. El parámetro interno se **deriva por inversión numérica**, y eso
exige un paso propio que antes no existía.

- ☐ **F4.9.1.** Perfil hermano `resources/f4/gps-hiseqx-pcr-pe150-v1.json`,
  idéntico al vigente en lectura, inserto, calidad y error, con preparación
  `pcr` y un bloque de amplificación con sus dominios; y variante `Pcr` en
  `LibraryPreparation`.
  - Depende de: cierre de la entrega 2 (por el orden aritmético de §3.1).
  - Procede de: REQ-CONF-82, REQ-CONF-89 reescrito; decisión **S-25**.
  - Cambia: recurso nuevo; `src/technical_profile.rs`.
  - Se comprueba con: el recurso vigente y su SHA-256 no cambian, comparados
    byte a byte con la línea base; el perfil nuevo se carga y declara su
    preparación. **Debe fallar contra:** una edición del recurso vigente para
    meter las dos preparaciones en uno, que rompería los valores fijos de los
    tres esquemas publicados y la huella congelada.

- ☐ **F4.9.2.** Constantes de contrato en una sola constante de Rust cada una
  —30 ciclos, 6 ciclos explícitos, 256 de tamaño máximo de familia, 16 intentos
  de aceptación, 4 clases de eficiencia— y un **juez** que las compare con los
  máximos de los esquemas y de la interfaz.
  - Depende de: F4.9.1.
  - Procede de: precedente del juez de techos de la sección 11.4 del plan de la
    fase de población fundadora.
  - Cambia: `src/technical_profile.rs`, los esquemas nuevos, `ui/`.
  - Se comprueba con: una prueba que compare cada constante con su máximo en
    cada esquema y en la interfaz, **con control positivo de que las
    comparaciones ocurrieron**. **Debe fallar contra:** un máximo cambiado en
    un solo esquema. *Motivo:* una copia decimal silenciosa es lo que sobrevivió
    a un cambio de contrato la vez anterior.

- ☐ **F4.9.3.** Configuración del modo amplificado y sus rechazos. El usuario
  declara **la fracción de duplicados esperada**, los **ciclos** y la
  **eficiencia por ciclo**, los tres **obligatorios y sin valor por omisión**;
  y opcionalmente ciclos explícitos, tamaño máximo de familia, intentos de
  composición, curva de composición y clases de eficiencia.
  - Depende de: F4.9.2.
  - Procede de: REQ-CONF-134/135/136 propuestos; decisiones **S-19** y
    **S-22**; decisión **S-26** para las versiones.
  - Cambia: los dos esquemas nuevos, `src/configuration.rs`, `src/adapters.rs`,
    `src/application.rs`.
  - **Qué significa exactamente el número que el usuario declara, y hay que
    escribirlo en el contrato.** Es la fracción **esperada para una biblioteca
    homogénea**, no la realizada, y la realizada es **mayor** por **tres**
    motivos acumulativos, no dos: las **clases de eficiencia**, la **curva de
    composición** y —el que faltaba— la **desigualdad de los primeros ciclos**,
    que concentra las parejas en menos familias **aunque no haya ni clases ni
    curva**. Basta con que los ciclos explícitos sean más de cero, que es el
    caso por omisión. Quien espere que las dos fracciones coincidan no las verá
    cuadrar nunca. El manifiesto publica la realizada aparte (F4.9.9).
  - **Y qué no es.** No es la fracción que marca una herramienta de análisis
    contando parejas que caen en el mismo sitio: esa incluye además las copias
    de la química de agrupación, que no son PCR. La documentación lo dice con
    el ejemplo trabajado de la decisión S-19 —medio punto porcentual atribuible
    a la amplificación dentro de una misma plataforma—, marcado como ejemplo y
    con la salvedad de que no es Illumina.
  - Se comprueba con: la batería de rechazos, cada uno citando parámetro y
    valor —preparación sin PCR con tabla de amplificación presente; `pcr` sin
    fracción de duplicados; sin ciclos; 31 ciclos; eficiencia 0; fracción de
    duplicados negativa o igual o mayor que 1; curva con eficiencia 0 o con
    nudos decrecientes; clases cuyos pesos no suman 1 dentro de la tolerancia
    de DTC-049—. **Debe fallar contra:** un esquema con
    `additionalProperties: true` en la tabla nueva; un dominio abierto en
    cualquiera de los parámetros; un valor por omisión materializado para la
    fracción de duplicados, que afirmaría una calibración inexistente.

- ☐ **F4.9.4.** **Inversión numérica** de la fracción de duplicados: dada la
  fracción esperada `f` que el usuario declara, derivar el parámetro interno
  `m` —parejas por molécula— que la produce en una biblioteca homogénea, donde
  `f(m) = 1 − (1 − e^(−m))/m`.
  - Depende de: F4.9.3.
  - Procede de: decisión **S-19**. **Este paso no existía**: la primera versión
    del plan declaraba `m` directamente y por eso no hacía falta. El diseño con
    PCR había descartado esta opción **precisamente porque exigía una inversión
    numérica cuya reproducibilidad habría que fijar**; ahora es obligatoria, y
    fijar esa reproducibilidad es el contenido del paso.
  - Cambia: `src/pcr_amplification.rs`.
  - **El método, fijado aquí, con sus dos mitades.** `f` es continua y
    **estrictamente creciente** en `m`, con `f(0) = 0` y límite 1, de modo que
    la inversa existe y es única. Se calcula por **bisección sobre el intervalo
    `[0, M_MAX]` con `M_MAX = 64` y exactamente 64 iteraciones**, sin criterio
    de parada por tolerancia. Agotar la resolución de un `f64` en ese intervalo
    cuesta entre 53 y 68 pasos según el tope, así que 64 bastan y sobran. No se
    usa Newton, cuyo número de pasos depende del punto de partida.
  - **La exponencial procede exclusivamente de `libm = =0.2.16`, nunca de la
    matemática estándar.** Esta mitad no es una elección de este plan: es
    política vigente del proyecto, literal en DTC-067 —«`exp`, `log`, `lgamma`,
    `floor` y `sqrt` proceden exclusivamente de `libm = =0.2.16`, nunca de
    matemática estándar»— y repetida en DTC-072 para `log` y `exp`. La versión
    fijada ya está en `Cargo.toml` y ya superó la evaluación de dependencias.
  - **Por qué hacen falta las dos mitades, y por qué la afirmación anterior era
    falsa.** La versión anterior de este paso decía que «un número fijo de
    iteraciones da el mismo resultado bit a bit en toda plataforma». **No es
    cierto**: la bisección compara `f(mid)` con el objetivo, y `f` contiene una
    exponencial; si esa exponencial es la de la biblioteca del sistema, su
    último bit varía entre plataformas, la comparación cae del otro lado y la
    bisección se va por otra rama **aunque el número de iteraciones sea el
    mismo**. El número fijo garantiza el mismo **coste**, no el mismo
    **resultado**. La revisión lo midió: **194 875 de dos millones** de
    exponenciales difieren en bits entre implementaciones, y **2 591 de 20 000**
    fracciones objetivo acaban dando una `m` distinta con 64 iteraciones fijas.
    La redacción correcta es: **mismo resultado en toda plataforma porque la
    exponencial es la de la librería fijada y el número de pasos es fijo.**
  - **Los dos casos límite, y el primero no es un adorno.** `f = 0` devuelve
    `m = 0` **por retorno anticipado, sin iterar**. Sin ese retorno la
    bisección devuelve un residuo del orden de `10⁻¹⁸` en vez de cero
    —comprobado al reproducirlo—, y ese residuo **rompe la identidad byte a
    byte** de F4.9.11, porque un `m` que no es exactamente cero hace que el
    tamaño de familia deje de ser trivial y se consuma aleatoriedad. Es la
    clase de defecto que no da error y cambia todos los FASTQ. El segundo: una
    `f` por encima de la que `M_MAX` puede producir se **rechaza antes de
    simular**, nombrando el parámetro, el valor y el máximo alcanzable; no se
    pinza en silencio.
  - **Lo que se publica.** La configuración efectiva publica la `f` declarada,
    la `m` derivada con precisión de ida y vuelta, el número de iteraciones y
    el nombre del método. La reproducibilidad de una ejecución se ancla en la
    `m` publicada, no en volver a derivarla.
  - Se comprueba con: la inversión de las fracciones que el diseño ya tiene
    calculadas en el sentido directo —`m = 0,02` da 1 %, `m = 0,1` da 4,8 % y
    `m = 0,5` da 21 %— devuelve esos `m` dentro de la precisión declarada,
    comparada contra el cálculo directo hecho **a mano**, no contra la propia
    inversión; `f = 0` devuelve `m = 0` exacto y el contador de iteraciones
    vale cero; una `f` inalcanzable se rechaza; la función directa es
    estrictamente creciente sobre una malla de cien puntos, que es la
    precondición que hace única la inversa; el número de iteraciones es
    **constante** para toda entrada distinta de cero; y —**esta es la que
    protege el determinismo**— una búsqueda negativa sobre el módulo,
    con control positivo, confirma que **no se llama a la exponencial de la
    matemática estándar** en ninguna ruta de la inversión, solo a la de
    `libm = =0.2.16`.
  - **Debe fallar contra:** una bisección con criterio de parada por tolerancia
    en vez de iteraciones fijas —la prueba del contador constante la caza—;
    **una implementación que use `f64::exp` de la matemática estándar**, que
    solo la búsqueda negativa caza, porque en una sola plataforma **los
    resultados coinciden y ninguna prueba numérica la distingue**; **la
    supresión del retorno anticipado de fracción cero**, que devolvería un
    residuo del orden de `10⁻¹⁸` y que la prueba de `m = 0` exacto y la de
    identidad byte a byte de F4.9.11 cazan; un pinzamiento silencioso de la
    fracción inalcanzable; una inversión validada contra sí misma, que se evita
    comparando con el cálculo directo a mano.
  - *Cuidado, y va en la docstring:* una prueba de reproducibilidad ejecutada
    **dos veces en la misma máquina** no discrimina nada sobre el determinismo
    entre plataformas, porque la exponencial del sistema es determinista dentro
    de un sistema. La versión anterior de este paso se apoyaba en esa prueba, y
    por eso no protegía lo que decía proteger.

- ☐ **F4.9.5.** Distribución del número de copias tras los primeros ciclos, por
  programación dinámica exacta, y su versión **sesgada por tamaño**, que es la
  que se sortea porque una familia observada es más probable cuantas más copias
  tenga su molécula.
  - Depende de: F4.9.4.
  - Procede de: diseño con PCR §5.5; Sun (1995) y Weiss y von Haeseler (1995).
  - Cambia: `src/pcr_amplification.rs`.
  - Se comprueba con: con 3 ciclos explícitos y eficiencia 0,5, la distribución
    en `1..8` coincide con una **enumeración a mano** de los ocho caminos, y la
    sesgada por tamaño tiene la media que le corresponde. **Debe fallar
    contra:** una programación dinámica que olvide el término «no se copia»;
    una que no sesgue por tamaño. *Cuidado:* la enumeración se escribe a mano;
    reusar la propia programación dinámica como oráculo sería circular.

- ☐ **F4.9.6.** Tamaño de familia: Poisson de media `λ` **condicionada a al
  menos una pareja** y **pinzada** en el máximo declarado, normalizada con la
  suma compensada de DTC-049 y sorteada con **un solo** `categorical-f64-v1`.
  - Depende de: F4.9.5.
  - Procede de: DTC-049, DTC-072; REQ-CONF-137 propuesto.
  - **Por qué no `poisson-count-v1`:** su muestreo por rechazo no acota el
    número de sorteos y su cola no está pinzada, y el proyecto ya fijó que un
    recuento se **pinza** en vez de rechazarse, «porque el rechazo no acota el
    número de sorteos y el presupuesto de operaciones sí». La decisión técnica
    nueva debe decirlo para que nadie lea DTC-067 como una prohibición.
  - Cambia: `src/pcr_amplification.rs`.
  - Se comprueba con: para `λ = 0,3` y máximo 4, los pesos coinciden con los
    calculados a mano y el último recoge **toda** la masa de la cola; con
    `λ = 0` no se consume aleatoriedad, contado con un generador instrumentado.
    **Debe fallar contra:** unos pesos normalizados sin pinzar la cola —la
    comparación con los pesos a mano lo caza—; un `λ = 0` que consuma.

- ☐ **F4.9.7.** Aceptación por composición de bases **al elegir la molécula**,
  con rechazo y **tope de intentos**; agotado el tope, **se acepta el último
  intento y se cuenta**.
  - Depende de: F4.9.6.
  - Procede de: REQ-CONF-94/95 y REQ-DAT-52 reescritos; decisión **S-20**;
    Benjamini y Speed (2012), que establecen que lo que manda es la composición
    **del fragmento completo**.
  - **Por qué al elegir y no al repartir copias.** Con pocas parejas por
    molécula —el caso normal— casi todas las familias tienen tamaño 1 tanto en
    regiones ricas como pobres en GC. La infrarrepresentación real viene de las
    moléculas que **no llegan a observarse**, y eso solo puede expresarse al
    elegir cuál da lugar a la familia.
  - Cambia: `src/pcr_amplification.rs`, `src/wgs_sampling.rs`.
  - Se comprueba con:
    **(a)** referencia con una mitad al 30 % de GC y otra al 70 %, curva que da
    factor 1 a la primera y 0,2 a la segunda, 10 ciclos y fracción de
    duplicados pequeña: la proporción de familias en la segunda mitad cae muy
    por debajo de 0,5, al orden que la fórmula predice.
    **(b)** con una curva que rechaza casi siempre y un tope de **2** intentos,
    el recuento de **familias aceptadas al agotar el tope es estrictamente
    mayor que cero**, y con 16 intentos es **menor** que con 2.
    **Debe fallar contra:** una implementación que aplique la composición solo
    al tamaño de familia, que con fracción de duplicados pequeña daría
    aproximadamente 0,5 (a); un rechazo **sin tope**, que nunca contaría
    ninguna familia aceptada al agotarlo y dejaría el recuento en cero (b).
  - *Cuidado, corregido:* la primera versión declaraba que la variante
    defectuosa de (b) era «un mutante que no termina», y **no es cierto**: con
    una curva que rechaza «casi siempre» el bucle **sí termina**, solo que
    tarde, de modo que la prueba quedaba dependiendo del reloj y era
    intermitente. El discriminador real es el **recuento de aceptadas al agotar
    el tope**, que un rechazo sin tope deja necesariamente en cero. Va en la
    docstring. Con fracción de duplicados grande la variante defectuosa de (a)
    también produce sesgo y no se distingue, así que el fixture la exige
    pequeña.

- ☐ **F4.9.8.** Muestreo por familias integrado en `WgsPairedEndSampler`,
  conservando su contrato externo: sigue produciendo exactamente `N` parejas
  con ordinales `0..N-1`. Cada pareja de una familia tiene **identidad técnica
  propia**. La quimera se sigue decidiendo **por pareja** sobre la molécula sin
  modificar, de modo que su comparación byte a byte se conserva intacta.
  - Depende de: F4.9.7.
  - Procede de: REQ-DAT-18, REQ-DAT-47 reescrito, REQ-DAT-81 propuesto,
    DTC-087; decisiones **S-22** y **S-23**.
  - Cambia: `src/wgs_sampling.rs`, `src/pcr_amplification.rs`,
    `src/randomness.rs` (etapa `amplification`, sin las decisiones de linaje).
  - Se comprueba con:
    **(a)** con errores de secuenciación desactivados, quimera a cero y sin
    loci, el número de parejas cuyas secuencias R1 y R2 coinciden con otra,
    **contado por un oráculo que solo lee el FASTQ**, es igual al recuento de
    duplicados del manifiesto; y con una fracción de duplicados alta es
    **mayor que cero**.
    **(b)** las parejas de una familia **no** aparecen adyacentes en el fichero
    publicado, porque la escritura ordena por identificador de instancia.
    **(c)** añadir una biblioteca no cambia las familias de las demás.
    **Debe fallar contra:** un manifiesto que publique el recuento esperado en
    vez del realizado; una implementación que sortee el tamaño de familia pero
    emita una sola pareja (a, cuya exigencia de «mayor que cero» evita que pase
    por el motivo equivocado); un ámbito indexado por orden de biblioteca (c).
    *Cuidado:* con los errores activados los duplicados dejan de ser idénticos
    y (a) no discrimina; el fixture los desactiva y la docstring lo dice.

- ☐ **F4.9.9.** Recuentos agregados del manifiesto: familias, parejas, parejas
  duplicadas, **fracción realizada** de duplicados, histograma de tamaños de
  familia, familias truncadas al final, intentos de composición totales y
  agotados, y la **cobertura de moléculas distintas**. Nada por pareja.
  - Depende de: F4.9.8.
  - Procede de: REQ-MAN-45 reescrito, REQ-WGS-11, REQ-DAT-56; decisión
    **S-28**.
  - Cambia: `src/event_manifest.rs`, `schemas/manifest-event-v3.schema.json`.
  - Se comprueba con: la suma del histograma ponderada por su índice es igual
    al número de parejas; la fracción realizada es uno menos familias entre
    parejas; **el par de ejecuciones que aísla el mecanismo**: con **cero
    ciclos explícitos**, curva neutra y sin clases, la fracción realizada
    coincide con la declarada dentro de la tolerancia muestral; con **tres
    ciclos explícitos** y todo lo demás igual, es **estrictamente mayor**; la
    cobertura de moléculas distintas se deriva de la realizada; validación
    externa del manifiesto.
  - **Debe fallar contra:** un manifiesto que copie la fracción declarada en el
    campo de la realizada —el par de ejecuciones lo caza, porque la segunda
    tendría que subir y no subiría—; una implementación que no propague la
    desigualdad de los primeros ciclos al tamaño de familia; un histograma que
    cuente familias en vez de parejas en alguna celda. *Cuidado, en la
    docstring:* una **sola** ejecución no discrimina, porque declarada y
    realizada quedan cerca; y la versión anterior de este paso usaba las clases
    de eficiencia como único mecanismo, cuando el de los primeros ciclos actúa
    ya por omisión y es el que hay que aislar.

- ☐ **F4.9.10.** Cota de recursos **publicada, no impuesta**: operaciones y
  memoria derivadas de la configuración y de la geometría conocida antes de
  abrir el primer ámbito aleatorio, publicadas como valores **calculados sin
  ejecutar** y no como medidas.
  - Depende de: F4.9.9.
  - Procede de: REQ-MAN-94 propuesto; decisión **S-27**, verificada: la fase de
    secuenciación **nunca ha tenido** control de recursos, de modo que imponer
    aquí un rechazo estrenaría algo que no le corresponde.
  - Cambia: `src/pcr_amplification.rs`, `src/event_manifest.rs`.
  - **Lo que no se acota, dicho expresamente.** Los sorteos por rechazo que el
    recorrido ya usa —la normal truncada del inserto y el índice uniforme— no
    tienen cota superior de palabras de generador en el peor caso, y el
    proyecto ya lo reconoce. Este trabajo no añade ningún sorteo por rechazo
    nuevo; sí multiplica por el tope de intentos los sorteos de inserto de una
    familia. La cota declara además que el recorrido regenera cada biblioteca
    **dos veces**.
  - Se comprueba con: un contador de unidades instrumentado, solo en pruebas,
    sobre un fixture con intentos, tamaño máximo de familia, ciclos explícitos
    y número de tractos **en sus máximos de contrato** y una curva que rechace
    casi siempre: las unidades reales son menores o iguales que la cota, con
    control positivo de que el contador contó. **Debe fallar contra:** un
    mutante que suprima el término de los intentos de la fórmula. *Cuidado:*
    con curva neutra el tope de intentos **no se ejercita** y el mutante
    sobrevive; sin el fixture saturado la prueba no discrimina.

- ☐ **F4.9.11.** Reducción al caso sin amplificar, y **rechazo temporal** de la
  combinación de preparación amplificada con loci de repetición declarados.
  - Depende de: F4.9.10.
  - Procede de: REQ-DAT-84 propuesto; decisión **S-16**; restricción aritmética
    de §3.1.
  - Cambia: `src/configuration.rs`, `src/diagnostics.rs`.
  - **El rechazo, y por qué.** Mientras la entrega 4 no exista, una biblioteca
    amplificada sobre un escenario con loci declarados tartamudearía **solo con
    la parte de base** y quedaría por debajo de lo medido —entre un 14 % y un
    28 % en motivos de 1 a 5 pares de base—. Antes que publicar una cifra que
    se sabe baja, la configuración **rechaza** esa combinación con un
    diagnóstico que nombra el paso pendiente. Se levanta en F4.10.6.
  - Se comprueba con: **(a)** con preparación amplificada, **fracción de
    duplicados cero**, curva neutra y tasas a cero, los FASTQ son byte a byte
    iguales a los de la preparación sin PCR, y el manifiesto difiere solo en
    los campos declarativos; **(b)** un generador instrumentado confirma cero
    consumo en la etapa de amplificación con esos parámetros; **(c)**
    preparación amplificada más un locus declarado se rechaza citando el paso
    pendiente. **Debe fallar contra:** una implementación que consuma
    aleatoriedad con parámetros neutros o que desplace los ordinales de
    plantilla (a y b); una inversión numérica que devuelva un `m` distinto de
    cero para fracción cero, que rompería (a); una que acepte en silencio la
    combinación (c).

- ☐ **F4.9.12.** Validador independiente de la amplificación: rederiva desde la
  semilla y la configuración efectiva las familias, las copias de los primeros
  ciclos y los tamaños de familia, **sin llamar al productor**.
  - Depende de: F4.9.11.
  - Procede de: REQ-VAL-58 propuesto, arquitectura §11.
  - Cambia: validador nuevo bajo `tests/`.
  - Se comprueba con: coincidencia con el manifiesto sobre dos fixtures;
    `./scripts/check-independent-validators.sh evidence/pcr-f4-9` da **5 PASS**.
    **Debe fallar contra:** un oráculo que lea el histograma publicado en vez
    de rederivarlo —se detecta mutando el histograma y exigiendo que el
    validador lo rechace. *Nota:* la independencia algebraica de la Poisson
    pinzada la da la prueba de pesos a mano de F4.9.6, no este oráculo; y la de
    la inversión numérica, la comparación con el cálculo directo de F4.9.4.

- ☐ **F4.9.13.** Interfaz Marimo de la preparación amplificada: selector de
  preparación, **fracción de duplicados**, ciclos, eficiencia, ciclos
  explícitos, tamaño máximo de familia, intentos, curva de composición y clases
  de eficiencia; y la **`m` derivada** mostrada al usuario junto a la fracción
  que ha declarado.
  - Depende de: F4.9.12.
  - Procede de: REQ-UI-01/03/04/05, hoja de ruta §7; decisión **S-19**.
  - Cambia: `ui/config_editor.py`, `ui/ui_panels.py`; prueba en `tests_ui/`.
  - Se comprueba con: ida y vuelta completa; la tabla de amplificación no
    aparece con preparación sin PCR; la `m` que la interfaz muestra coincide
    con la que la configuración efectiva publica para la misma fracción; los
    límites de los controles coinciden con las constantes de contrato,
    comprobado por el juez de F4.9.2; `marimo check --strict` sin avisos.
    **Debe fallar contra:** una interfaz que calcule la `m` con un método
    distinto del de F4.9.4, que la comparación caza; un control cuyo máximo no
    coincida con la constante.

> **Revisión del hito de la entrega 3.** Comprobar la inversión numérica y su
> determinismo, las familias, los duplicados contados por bytes, la composición
> al elegir la molécula, el tope de intentos por su recuento, la cota
> publicada, la reducción al caso sin PCR y el rechazo temporal. Ejecutar la
> batería completa con
> `./scripts/check-independent-validators.sh evidence/pcr-f4-9`. Aplicar el
> cambio de contrato de la sección 10 en su parte de preparación amplificada y
> **repetir la evidencia de cierre de la fase de realismo técnico con
> preparación sin PCR**, que debe dar los mismos bytes. No continuar hasta
> confirmar el resultado.

---

## 8. Entrega 4 — Tartamudeo clonal heredable por linaje (F4.10)

Resultado: con preparación amplificada, el deslizamiento puede ocurrir en
**cualquier ciclo**, y una copia que resbala transmite su número de unidades a
**toda su descendencia**. Un resbalón en el primer ciclo alcanza a media
familia; uno en el último, a una sola copia. Una familia de duplicados puede
así mostrar dos o tres alelos en proporciones que imitan a un heterocigoto,
cuando procede de una sola molécula.

- ☐ **F4.10.1.** Extensión del catálogo de ámbitos aleatorios con las
  decisiones de linaje: progenitor, ciclo de divergencia, resbalón por ciclo y
  tracto, y paso; indexadas por **identificador de locus**, nunca por posición.
  - Depende de: cierre de las entregas 2 y 3.
  - Procede de: REQ-DAT-85 propuesto, DTC-081, DTC-086; decisión **S-11**.
  - Cambia: `src/randomness.rs`.
  - Se comprueba con: vectores congelados; las etapas existentes producen las
    mismas subsemillas que antes; **añadir un locus no cambia los resbalones de
    los loci previos**, ni en otro cromosoma ni —y esta es la que importa— **a
    la izquierda del anterior dentro de la misma molécula**. **Debe fallar
    contra:** un ámbito indexado por la posición del tracto dentro de la
    molécula.
  - **Corrección del 26/09/2026.** Este paso decía que «la prueba del locus
    añadido en otro cromosoma es la única que lo caza». **Es falso, y se midió
    mutando el ámbito para que usara la posición del tracto dentro de la
    molécula: esa prueba pasó sin inmutarse.** Dos tractos de cromosomas
    distintos nunca comparten molécula, de modo que una implementación por
    posición le sigue dando el índice cero a los dos. La que lo caza es la del
    locus añadido **a la izquierda**, que es donde el tracto anterior deja de
    ser el primero. Las dos pruebas existen; la afirmación de arriba queda
    corregida para que nadie borre la buena creyendo que sobra.

- ☐ **F4.10.2.** Genealogía mínima de la familia: la primera pareja tiene
  linaje propio completo de un paso por ciclo; la pareja siguiente elige un
  **progenitor uniforme** entre las anteriores y un **ciclo de divergencia**
  geométrico, hereda los resbalones del progenitor **anteriores** a la
  divergencia y sortea los posteriores.
  - Depende de: F4.10.1.
  - Procede de: REQ-DAT-82 propuesto; decisión **S-21**; diseño con PCR §5.7.
  - **Declarado como aproximación**, no como el coalescente exacto del proceso
    de ramificación, del que no se afirma tener forma cerrada.
  - Cambia: `src/pcr_amplification.rs` o módulo hermano.
  - Se comprueba con: **fixture saturado** —un locus dentro de la molécula,
    fracción de duplicados alta para forzar familias grandes, divergencia
    forzada al último ciclo para todas las parejas y probabilidad de resbalón 1
    en el primer ciclo y 0 en los demás, por sustitución solo en prueba—:
    **todas** las parejas de cada familia muestran el **mismo** número de
    unidades, distinto del verdadero. Segunda variante con divergencia forzada
    a cero: las parejas son **independientes**. **Debe fallar contra:** un
    tartamudeo sorteado por pareja o por lectura; una herencia que copie los
    resbalones **posteriores** a la divergencia en vez de los anteriores.
    *Cuidado, en la docstring:* sin forzar la divergencia y sin probabilidad 1,
    algunas discordancias dentro de la familia son legítimas y la aserción no
    discrimina.

- ☐ **F4.10.3.** Derivación de la probabilidad de resbalón por ciclo desde las
  dos columnas de la calibración y el número de ciclos, con sustitución directa
  admitida. Si la fracción con PCR no supera a la de sin PCR, la probabilidad
  vale cero.
  - Depende de: F4.10.2.
  - Procede de: REQ-CONF-138 propuesto; decisión **S-24**.
  - Cambia: `src/pcr_amplification.rs`, `src/technical_profile.rs`.
  - **Honestidad sobre lo que esta resta puede y no puede afirmar.** Las dos
    columnas proceden de **dos individuos distintos**, de proyectos distintos y
    con coberturas distintas —no de las mismas células, como este plan afirmó
    por error—. La resta sigue siendo el modo correcto de componer, porque cada
    columna es un total; pero puede llevar dentro efectos de individuo y de
    lote, de modo que **las tasas por ciclo derivadas son toscas**. El
    manifiesto las publica marcadas como tales, y la documentación lo dice.
  - Se comprueba con: con diez ciclos, las probabilidades por ciclo coinciden a
    tres cifras con el cálculo hecho **a mano** desde las dos columnas; la
    configuración efectiva publica las dos fracciones de entrada, el número de
    ciclos, la fórmula y la marca de tosquedad; una sustitución directa del
    usuario prevalece y se marca como suya; con fracción con PCR igual o menor
    que la de sin PCR, la probabilidad vale cero y no se consume aleatoriedad.
    **Debe fallar contra:** una implementación que **sume** la columna sin PCR
    a la de con PCR, que contaría el mismo fenómeno dos veces; una que use la
    fracción de extremo a extremo directamente como probabilidad por ciclo; una
    publicación sin la marca de tosquedad, que presentaría como medida una
    cifra que arrastra efectos de individuo y de lote.

- ☐ **F4.10.4.** Aplicación en el recorrido: sobre la molécula ya quimérica si
  lo es, **primero** el paso neto heredado del linaje de PCR y **después** la
  decisión de base por pareja de la entrega 2. Los tractos se procesan en orden
  canónico de inicio sobre la molécula, y cada paso se aplica **en la tirada
  ininterrumpida más larga** mediante la operación de F6.4.1.1. Las guardas de
  viabilidad se aplican al **paso neto acumulado**, excluyendo antes de sortear
  y sin reintentar.
  - Depende de: F4.10.3.
  - Procede de: REQ-DAT-76 y REQ-DAT-82 propuestos; decisiones **S-07**,
    **S-13**, **S-14**, **S-15**.
  - Cambia: `src/stutter.rs`, `src/pcr_amplification.rs`, `src/execution.rs`.
  - Se comprueba con: con el paso de PCR fijado y la tasa de base a 1 con paso
    fijado, el resultado es la **suma** de los dos; un paso neto que dejaría el
    tracto sin ninguna unidad completa queda excluido y no consume; un corpus
    con tasa de quimera alta no produce `ChimeraError::PairSourceMismatch`; tras
    dos pasos sucesivos, la **tirada más larga se recalcula** entre uno y otro,
    de modo que el segundo puede caer en otra tirada. **Debe fallar contra:**
    un orden invertido; una guarda aplicada a cada paso por separado en vez de
    al neto; una implementación que calcule la tirada más larga **una sola vez**
    y aplique los dos pasos en la misma, que la prueba de los dos pasos caza.

- ☐ **F4.10.5.** **Comprobación aritmética de la composición**: con las dos
  partes activas, la fracción de lecturas afectadas por tartamudeo que un
  oráculo cuenta sobre los FASTQ reproduce la fracción **medida** con PCR.
  - Depende de: F4.10.4.
  - Procede de: la restricción aritmética de §3.1, que es la razón del orden de
    todo este plan.
  - Cambia: prueba de sistema nueva bajo `tests/`.
  - **Qué valida exactamente, dicho sin adornos y obligatorio en su docstring.**
    La recomposición `1 − (1−r_sin)(1−r_pcr) = r_con` es **exacta por
    construcción**: la derivación de la probabilidad por ciclo se hizo
    despejando justamente eso. Por tanto esta prueba comprueba que **la
    implementación compone como el diseño dice**, no que el modelo sea
    correcto. Es **tautológica respecto del modelo** y discriminante respecto
    de la implementación, y se presenta así. Lo que sí caza, y es su razón de
    ser, son los dos errores de implementación que el orden de este plan
    existe para evitar.
  - **Tamaño muestral, calculado y no estimado a ojo.** Para distinguir la
    fracción correcta de la que produciría una entrega sin la parte de base
    hace falta un número de moléculas que cubran el tracto igual a
    `n ≈ (z_{α/2}·√(p₀(1−p₀)) + z_β·√(p₁(1−p₁)))² / (p₀−p₁)²`. Con α = 0,05 a
    dos colas y potencia del 90 %: para el homopolímero, `p₀ = 0,17` frente a
    `p₁ = 0,147`, **n ≈ 2700**; para el dinucleótido, `p₀ = 0,038` frente a
    `p₁ = 0,0273`, **n ≈ 3000**. La docstring escribe α, la potencia, el `n`
    empleado y la tolerancia que de ahí sale; **no** se escribe una tolerancia
    a ojo, que es lo que hacía la primera versión y que habría dejado pasar las
    dos variantes defectuosas.
  - Se comprueba con: un corpus con un homopolímero y un dinucleótido
    declarados, con al menos ese número de moléculas cubriendo cada tracto y
    semilla congelada: con preparación amplificada la fracción observada cae
    dentro del intervalo alrededor de 0,17 y 0,038; con preparación sin PCR,
    alrededor de 0,027 y 0,011. **Debe fallar contra:** una entrega de la
    preparación amplificada **sin** la parte de base, que daría 0,147 y 0,0273;
    una implementación que sume las dos columnas, que daría por encima de lo
    medido.

- ☐ **F4.10.6.** Cierre: recuentos de resbalones y de resbalones heredados por
  locus en el manifiesto; validador independiente que rederiva los pasos netos
  de linaje desde la semilla; control en la interfaz; y **levantar el rechazo
  temporal** de F4.9.11.
  - Depende de: F4.10.5.
  - Procede de: REQ-MAN-45 y REQ-VAL-58 propuestos, hoja de ruta §7.
  - Cambia: `src/event_manifest.rs`, `ui/`, validador bajo `tests/`,
    `src/configuration.rs`.
  - Se comprueba con: los resbalones heredados son mayores que cero en el
    fixture de familias grandes y **exactamente cero** con familias de tamaño
    1; el validador coincide sin llamar al productor;
    `./scripts/check-independent-validators.sh evidence/clonal-f4-10` da **5
    PASS**; la combinación antes rechazada se acepta y produce la fracción de
    F4.10.5; ida y vuelta de la interfaz. **Debe fallar contra:** un recuento
    de heredados que cuente **todos** los resbalones —la prueba de las familias
    de tamaño 1 lo caza, porque ahí no puede heredarse nada—; un levantamiento
    del rechazo sin la prueba de F4.10.5 verde.

> **Revisión del hito de la entrega 4.** Comprobar la genealogía sobre el
> fixture saturado, la derivación de probabilidades y su marca de tosquedad, el
> orden de aplicación con recálculo de la tirada, la composición aritmética con
> su tamaño muestral y el levantamiento del rechazo. Ejecutar la batería
> completa con `./scripts/check-independent-validators.sh
> evidence/clonal-f4-10`. Cerrar el cambio de contrato de la sección 10. No
> continuar hasta confirmar el resultado.

---

## 9. Entrega 5 — Generador de alelos por deslizamiento biológico

**Condición de cierre, no ampliación aplazada.** Mientras esta entrega no
exista, los alelos y sus frecuencias los escribe el usuario, y **la simulación
de microsatélites no puede darse por cerrada**. Se le da entrega propia para
que se vea, en vez de quedar escondida como un paso más.

> **Esta entrega ya no está bloqueada.** *(Corregido el 27 de septiembre de
> 2026: este párrafo decía que lo estaba, y contradecía a la tabla de
> pendientes del propio plan, que la desbloqueó nueve días antes.)* El
> pendiente **P-4** —cómo un servicio de pasos produce una distribución de
> frecuencias y no una sola secuencia de alelos— **quedó resuelto el 20 de
> septiembre de 2026**: el responsable adoptó la forma candidata tal cual y
> fijó que el alelo de la referencia se pierda o no **según los parámetros que
> se den**, sin forzar su presencia. Los pasos F6.4.11 a F6.4.13 están
> desbloqueados, y F6.4.4 nunca lo estuvo.
>
> Las entregas 1 a 4 no dependían de P-4 en absoluto: mientras no existió, los
> alelos y sus frecuencias los declaraba el usuario, que es lo que la entrega 1
> ya entrega.

**Qué es y qué no es.** Es **biología, no artefacto**: parte del alelo de la
referencia, aplica pasos de una unidad y produce el abanico de alelos con sus
frecuencias, que después se sortean entre los individuos. Va en la **verdad**.
Es **una sola tirada para inicializar la diversidad**, no una simulación de
generaciones: el modelo escalonado aplicado generación tras generación
pertenece a la fase de evolución y no se invade.

**Por qué no puede confundirse con el tartamudeo, aunque el mecanismo físico
sea el mismo.** Uno crea variación real heredable y va en la verdad; el otro
solo ensucia la medida y va en las lecturas. La prueba de que no son
intercambiables: si los alelos verdaderos salieran del modelo de tartamudeo,
**cambiar el número de ciclos de PCR cambiaría el genotipo real de los
individuos** —un parámetro de laboratorio decidiendo la biología—.

- ☑ **F6.4.4.** Modelo escalonado simple como **servicio reutilizable**: cada
  evento añade o elimina exactamente **una unidad completa** según el sesgo
  efectivo, con tasa y sesgo configurables; rechazo previo de otro modelo o de
  un salto de varias unidades, con parámetro y valor.
  - Depende de: cierre de la entrega 1.
  - Procede de: REQ-CONF-61/62, REQ-IND-35, REQ-DIR-04, REQ-MUT-11,
    REQ-VAL-48, ACE-32/109; DTF-015 mantiene fuera el salto múltiple. Es el
    paso F6.4.4 del plan de la fase de variación compleja, **adelantado**.
  - Cambia: módulo nuevo o hermano de `src/tandem_repeat.rs`; reutiliza la
    operación de paso de F6.4.1.1, **incluida la regla de la tirada más
    larga**.
  - **Qué hace el servicio cuando el tracto se queda sin unidades**, que
    hacía falta decir porque mil contracciones sobre un tracto finito lo
    agotan. Se aplica **la misma regla que el tartamudeo** y por la misma
    razón: un paso que dejaría menos de una unidad completa en la tirada más
    larga **no es viable**, se excluye **antes** de sortear y no se reintenta;
    si no queda ningún paso viable, el evento **no se aplica**, no consume
    aleatoriedad y se **cuenta aparte** como evento no aplicado. No es un
    error: un microsatélite que ha llegado a su mínimo es un resultado
    biológico legítimo.
  - Se comprueba con: mil eventos con sesgo 1 son todos expansiones y con sesgo
    0 todas contracciones, sin consumo de aleatoriedad en los extremos, como
    `bernoulli-v1` exige; cada evento **aplicado** cambia la longitud
    exactamente en una unidad; **mil eventos con sesgo 0 sobre un tracto de
    cinco unidades terminan con el tracto en una unidad, sin error, y el
    recuento de eventos no aplicados vale exactamente 996**; un modelo distinto
    del admitido y un salto de dos unidades se rechazan antes de simular; misma
    semilla reproduce la secuencia de decisiones. **Debe fallar contra:** un
    paso que cambie dos unidades; un sesgo interpretado al revés —el caso con
    sesgo 1 exige solo expansiones—; un extremo que consuma aleatoriedad; **un
    servicio que falle, que sature en silencio o que deje el tracto en cero
    unidades al agotarlo**, que la prueba del recuento exacto de no aplicados
    caza. *Cuidado:* sin fijar ese recuento, un servicio que simplemente
    ignorase los eventos sobrantes y otro que los contase darían el mismo
    tracto final y la prueba no discriminaría.

- ☑ **F6.4.11.** **Generador del abanico de alelos**, con su **ámbito de azar
  propio** separado del tartamudeo. **Algoritmo adoptado el 20/09/2026**; ya no
  hay bloqueo.
  - Depende de: F6.4.4, F6.4.3.2.
  - Procede de: decisión **S-04** y su recuadro «Condición de cierre, no
    ampliación aplazada»; REQ-DAT-15, DTC-037, DTC-072, DTC-081.
  - Cambia: `src/randomness.rs` (ámbito nuevo), el módulo del generador,
    `schemas/config-general-v4.schema.json`.

  - **Algoritmo, adoptado el 20 de septiembre de 2026.** Se simulan **tantos
    linajes independientes como copias haplotípicas** tenga la población. Cada
    linaje parte del **alelo de la referencia** y recibe un **número de pasos
    sorteado de una Poisson de media declarada**, con `poisson-count-v1`, que
    es el único recuento normativo del proyecto; cada paso lo da el servicio
    escalonado de F6.4.4 con su sesgo. La **frecuencia de cada alelo es la
    proporción de linajes que terminan en él**. El alelo de la referencia entra
    en el abanico **si algún linaje termina en él**, y no se fuerza su
    presencia. Un linaje que rebasaría el tope de F6.4.3.1 **se detiene ahí**,
    igual que el servicio hace al agotar el tracto.

  - **Criterio del responsable sobre la pérdida del alelo original.** Que el
    alelo de la referencia desaparezca del abanico **depende de los parámetros
    que se den**, no del modelo: con una media de pasos baja la mayoría de los
    linajes no se mueven y el alelo está; con una media alta puede perderse,
    que es lo que ocurre en una población real. Por eso no se fuerza.

  - **Tres limitaciones declaradas, porque son decisiones biológicas y no
    detalles.** Deben ir en la documentación del generador y en el manifiesto:
    1. **No hay ancestría compartida.** Cada copia deriva por su cuenta desde
       el ancestro, de modo que dos individuos emparentados no comparten
       historia en este locus. El espectro de frecuencias resultante es **más
       plano** que el de una población real, donde la deriva concentra unos
       pocos alelos. Importa si el conjunto se usa para probar herramientas de
       estructura poblacional.
    2. **El número de pasos no es un proceso en el tiempo.** Es un atajo: cada
       copia acumula un número aleatorio de mutaciones sin decir en cuántas
       generaciones. Es coherente con que esto sea **una sola tirada de
       inicialización** y no una simulación de generaciones, que pertenece a la
       fase de evolución.
    3. **El genoma de referencia puede llevar un alelo que ningún individuo
       tiene.** Es real y ocurre, pero conviene saberlo antes de encadenar un
       análisis que lo dé por imposible.

  - **Dos ámbitos separados, y esta parte no está bloqueada.** El deslizamiento
    biológico y el tartamudeo comparten mecanismo físico y **no pueden
    compartir ámbito de azar**, porque uno decide la verdad y el otro las
    lecturas. El ámbito de este paso se identifica por el evento y el
    **identificador de locus**, nunca por posición ni por orden.
  - Se comprueba con: **vectores congelados de derivación de subsemilla** para
    el ámbito nuevo, al estilo de los de F4.8.1; **una aserción directa de que
    la subsemilla del deslizamiento biológico y la del tartamudeo son
    distintas** para la misma terna de identidades; misma semilla, mismo
    abanico; el abanico no contiene alelos con fracciones de unidad;
    **insertar un locus que ordena antes que los existentes no cambia el
    abanico de los previos**; con cero eventos el abanico es solo el alelo de
    la referencia y **no se consume aleatoriedad**.
  - **Debe fallar contra:** un generador que **comparta ámbito** con el
    tartamudeo, que solo la aserción de subsemillas distintas y los vectores
    congelados cazan; un abanico con alelos fraccionarios; un ámbito indexado
    por posición.
  - *Corregido, y es un error instructivo.* La versión anterior de este paso
    ofrecía como discriminador «cambiar los ciclos de PCR no cambia ni un alelo
    de la verdad». **No caza un ámbito compartido**, y conviene entender por
    qué: en este proyecto compartir ámbito significa **derivar la misma
    subsemilla**, y como DTC-081 exige que cada decisión arranque su propio
    generador desde su clave y que no haya consumo acumulado, un consumidor no
    desplaza al otro; además la verdad se genera **antes** que las lecturas, de
    modo que ningún parámetro de estas puede alcanzarla. La aserción se cumple
    sola. Lo que sí caza es comparar las dos subsemillas. La comprobación
    antigua se conserva en la pausa del hito como prueba de sistema —es buena
    contra un generador que leyera la configuración del tartamudeo por error—,
    pero **no** como prueba del ámbito, y su docstring lo dice.

- ☑ **F6.4.12.** Publicación: la configuración efectiva y el manifiesto
  registran el **origen de cada alelo** —declarado por el usuario o generado
  por deslizamiento—, los parámetros del generador, su versión de algoritmo y
  su versión de catálogo de ámbitos, y las frecuencias **declaradas y
  realizadas**.
  - Depende de: F6.4.11. Su bloqueo por P-4 **quedó levantado el 20/09/2026**.
  - Procede de: REQ-MAN-29/30, REQ-VAL-48; decisión **S-04**.
  - Cambia: `src/event_manifest.rs`, `src/results.rs`, los esquemas nuevos.
  - Se comprueba con: validación externa; un escenario con alelos declarados a
    mano los marca como declarados y uno generado los marca como generados, en
    el mismo documento si hay loci de los dos tipos; el número de eventos de
    deslizamiento registrado coincide con el configurado. **Debe fallar
    contra:** un manifiesto que marque todos los alelos con el mismo origen
    —el escenario mixto lo caza, y uno de un solo tipo **no discrimina**, lo
    que va en la docstring.

- ☑ **F6.4.13.** Interfaz Marimo del generador y **cierre de la simulación de
  microsatélites**: control para elegir entre declarar los alelos o generarlos,
  con los parámetros del generador; y comprobación de que ya no queda ninguna
  opción configurable de este trabajo sin control.
  - Depende de: F6.4.12. Su bloqueo por P-4 **quedó levantado el 20/09/2026** en la parte del
    generador; la comprobación de recorrido de opciones no está bloqueada.
  - Procede de: REQ-UI-01/03/04/05, hoja de ruta §7; decisión **S-17**.
  - Cambia: `ui/config_editor.py`, `ui/ui_panels.py`; prueba en `tests_ui/`.
  - Se comprueba con: ida y vuelta completa; elegir «generar» oculta la tabla
    de alelos y muestra los parámetros, y al revés; **una comprobación que
    recorra la lista de opciones configurables publicadas en la configuración
    efectiva por este trabajo y exija que cada una tenga control en la
    interfaz**, con control positivo de que la lista no estaba vacía;
    `marimo check --strict` sin avisos. **Debe fallar contra:** una opción
    añadida en cualquiera de las cinco entregas sin control en la interfaz —la
    comprobación de recorrido es la única que lo caza de forma sistemática, y
    es la que permite declarar cerrada la regla de cierre de la hoja de ruta.

> **Revisión del hito de la entrega 5, y cierre de la simulación de
> microsatélites.** Comprobar el servicio escalonado, el abanico de alelos, la
> separación demostrada de los dos ámbitos de azar, la publicación del origen
> de cada alelo y el recorrido completo de opciones en la interfaz. Ejecutar la
> batería completa con
> `./scripts/check-independent-validators.sh evidence/slippage-f6-4-11`. Solo
> al cerrar esta pausa puede declararse cerrada la simulación de
> microsatélites. No continuar hasta confirmar el resultado.

---

## 10. Cambio de contrato posterior al cierre de dos fases

Con la forma de la sección 11 del plan de la fase de población fundadora. **No
es un añadido silencioso y no debe presentarse como tal.**

La primera versión de este plan hablaba de reabrir el contrato de **una** fase,
la de realismo técnico. Las decisiones del 20 de septiembre añaden **una
segunda**: el generador de referencia simulada con regiones definidas y la
publicación de la correspondencia de coordenadas caen en la fase del **núcleo
reproducible**, también cerrada.

### 10.1. Qué **no** se reabre

Ningún hito de las dos fases se reabre, ni se revoca ninguna de sus decisiones,
ni se rehace ninguna de sus evidencias salvo las que la sección 10.4 obliga a
repetir. Lo que cambia son contratos posteriores y acotados, por decisión
expresa del responsable del 19 y del 20 de septiembre de 2026.

### 10.2. Qué cambia, y por qué no es una contradicción

- **En la fase del núcleo reproducible.** El generador de referencia simulada
  gana **regiones declaradas con su secuencia**, y la referencia con un tracto
  insertado se publica como **derivada con identidad propia**, con su mapa de
  correspondencia. Esto último es una **extensión** del mapa que la referencia
  reducida ya publica —con intervalo de origen y derivado por segmento—, no
  maquinaria nueva; el contrato que se amplía es el del mapa, no el de las
  coordenadas.
- **El tartamudeo sin PCR no contradice la preparación sin PCR.** Los cinco
  requisitos que declaran la preparación exclusivamente sin PCR hablan de la
  **preparación de biblioteca**. El deslizamiento ocurre además en la
  amplificación en puente sobre la celda y en la síntesis, que existen en el
  perfil vigente. El requisito nuevo debe decirlo para que nadie lo lea como
  una grieta.
- **La preparación amplificada sí reabre una ampliación aplazada.** DTF-056 se
  **divide**: salen de ella la preparación amplificada, los duplicados de PCR,
  el sesgo de composición **configurable**, la sobrerrepresentación estocástica
  y el tartamudeo heredable; **permanecen** el perfil de sesgo GC calibrado con
  datos reales, los duplicados ópticos y de agrupación, y la combinación con
  excepciones regionales. DTF-055 no cambia.
- **La especificación se corrige en un punto de fondo**, no solo se amplía: el
  umbral de 50 pb deja de aplicarse a los alelos de repetición declarados sobre
  loci declarados. Está en la sección 11.

### 10.3. Dónde vive cada literal que hay que tocar

Una copia silenciosa de un literal es lo que sobrevive a un cambio de contrato.
Los sitios conocidos, comprobados en el árbol de partida:

| Sitio | Papel | Qué dice hoy |
|---|---|---|
| `src/technical_profile.rs`, líneas 28–30 | canónico | `LibraryPreparation` con una sola variante, `PcrFree` |
| `src/application.rs`, líneas 776–786 | constantes | duplicados de PCR no aplicables y ópticos a cero |
| `src/event_manifest.rs`, líneas 1516–1531 | contrato | exige preparación sin PCR, PCR no aplicable y composición neutra |
| `schemas/effective-config-v4.schema.json`, líneas 312–313 | contrato | valores fijos del identificador del perfil, su ruta, su SHA-256, la preparación, la aplicabilidad de la PCR y el modo de composición |
| `schemas/manifest-event-v2.schema.json`, líneas 115–118 y 157 | contrato | los mismos valores fijos |
| `schemas/config-general-v2.schema.json`, `$defs/technical_event` | contrato | objeto cerrado de once claves; ninguna nombra la PCR, de modo que hoy una clave de PCR se rechaza **como campo desconocido**, no como capacidad aplazada |
| `src/reference_generation.rs` | canónico | generador que sabe **rellenar al azar y colocar regiones declaradas**, en ese orden. *Corregido el 24/09/2026: decía «solo sabe rellenar al azar», y dejó de ser cierto al cerrar la entrega 1; ver la sección 14 del plan del núcleo reproducible* |
| `src/reference_reduction.rs` y `schemas/reference-reduction-map-v1.schema.json` | contrato | el mapa de coordenadas que hay que extender, con `source_interval` y `derived_interval` por segmento |
| `ui/config_editor.py`, `ui/gps_app.py`, `ui/ui_panels.py` | contrato | controles del evento técnico y de la referencia. *Corregido el 24/09/2026: la fila nombraba dos ficheros y son tres; el cuaderno `gps_app.py` es el que reconstruye el formulario, y omitirlo de esta lista es exactamente cómo el tartamudeo llegó al lector del TOML y no a la interfaz* |
| `resources/f4/gps-hiseqx-pcrfree-pe150-v1.json` | recurso | perfil único, con su huella congelada |

**Regla:** los esquemas publicados **no se editan**. Las versiones nuevas
—general v4, efectivo v6, manifiesto de evento v3, manifiesto general v6,
verdad v3 y mapa de coordenadas v2— llevan la capacidad, y los valores fijos de
las versiones heredadas se comparan byte a byte en cada hito. Un fichero
general v2 sigue publicando efectivo v4 y manifiesto de evento v2.

### 10.4. Evidencia que hay que repetir

- Al cerrar la **entrega 1**, la evidencia de cierre de la fase del núcleo
  reproducible que toque la generación de referencias, con un escenario **sin
  regiones declaradas**, que debe dar los mismos bytes que la línea base.
- Al cerrar la **entrega 2** y otra vez al cerrar la **entrega 3**, la
  evidencia de cierre de la fase de realismo técnico **con preparación sin PCR
  y sin loci declarados**, que debe dar los mismos bytes.

Esas repeticiones son las que sostienen la afirmación de que los contratos se
han ampliado sin romperse.

### 10.5. Dónde se registra

- Sección **16 nueva** del plan de la fase de realismo técnico, titulada
  «Cambio de contrato posterior al cierre — tartamudeo y preparación
  amplificada (fecha)», con: qué no se reabre, qué cambia y por decisión de
  quién, dónde vive cada literal, las pruebas nuevas con sus variantes
  defectuosas, los fixtures que se reabren y la repetición de la evidencia.
- Sección **14 nueva** del plan de la fase del núcleo reproducible, con la
  misma forma, para el generador con regiones y el mapa de coordenadas.
- Una nota fechada en las dos secciones correspondientes de la hoja de ruta.

---

## 11. Correcciones que este plan **propone** a otros documentos

**Estado al 26 de septiembre de 2026.** Esta sección se escribió cuando ninguna
corrección estaba aplicada, y esa frase ha caducado. Hoy:

- las del **documento de calibración** están aplicadas desde el principio, y lo
  estaban por un motivo distinto del de las demás: sus hechos falsos eran del
  propio proyecto, y dejar publicada una afirmación que se sabe falsa no es una
  propuesta pendiente sino un defecto;
- las de las **entregas 2 y 3** están aplicadas al cerrar cada hito, con sus
  decisiones normativas escritas —el modelo del tartamudeo de base y el de la
  preparación amplificada—;
- las de la **entrega 4** se aplican al cerrar este hito, con su decisión
  normativa propia y las notas fechadas en la hoja de ruta y en el plan de la
  fase de realismo técnico.

Lo que sigue pendiente es lo que cada fila diga, no la sección entera.

| Documento | Dónde | Qué | Origen |
|---|---|---|---|
| `especificaciones-0.1.md` | §4.3, glosario de la repetición en tándem | La pureza es **siempre** una salida calculada, nunca una entrada; el alelo se define por su **secuencia** y no por un número de unidades; dos alelos de igual longitud y distinta interrupción son alelos distintos | S-02, S-06 |
| `especificaciones-0.1.md` | §4.3, mismo párrafo | El resbalón y el paso de deslizamiento ocurren en la **tirada ininterrumpida más larga** del motivo, con desempate por la izquierda. **Sustituye a la redacción «por el extremo derecho» que la primera versión de este plan proponía y que la decisión S-07 rechazó** | S-07 |
| **`especificaciones-0.1.md`** | **REQ-CONF-55 y REQ-IND-31** | **Corregir**: el umbral que separa inserciones y deleciones pequeñas de variantes estructurales se aplica a secuencia **nueva**; **un alelo de repetición en tándem declarado sobre un locus declarado no se somete a él**, ni al generar, ni al clasificar, ni al contar. La coherencia que REQ-IND-31 exige se conserva: se aplica igual en los tres momentos, y a las repeticiones no se aplica en ninguno | **S-05** |
| **`especificaciones-0.1.md`** | **REQ-VER-20** | **Corregir**: hoy exige alelos explícitos «cuando su tamaño sea manejable» y **representación simbólica para tramos largos**. Para las repeticiones en tándem eso contradice lo decidido: el criterio del estándar es **si se conoce la secuencia exacta**, y un simulador siempre la conoce. Debe decir que los alelos de repetición se publican **siempre** de forma explícita, sin umbral de tamaño | **S-09** |
| `especificaciones-0.1.md` | §4.4, frase de orden normativo, y REQ-DAT-17 | Insertar el tartamudeo entre la quimera y calidad y errores, y las familias de amplificación antes de la quimera | S-13, diseño con PCR §8 |
| `especificaciones-0.1.md` | §4.4, tras REQ-DAT-73 | Añadir REQ-DAT-74 a 80 (tartamudeo) y REQ-DAT-81 a 85 (amplificación) | Expediente §10.2, diseño §9.2 |
| `especificaciones-0.1.md` | **REQ-DAT-80** | **Reescribir**: la identidad byte a byte vale para los escenarios **que no declaran loci de repetición** | C-1 |
| `especificaciones-0.1.md` | §4.1, tras REQ-CONF-130 | Añadir REQ-CONF-131 a 133 (tartamudeo) y 134 a 138 (amplificación) | Expediente §10.2, diseño §9.2 |
| ~~`especificaciones-0.1.md`~~ | ~~**REQ-CONF-132**~~ | ~~**Reescribir**: lo desactivado por omisión es el tartamudeo **heredable de PCR**; el de base se activa al declarar loci~~ — **SUPERADA el 24/09/2026 por decisión del responsable**: «si el usuario mete el microsatélite, tiene que encender el tartamudeo y aportar los parámetros». El tartamudeo de base lleva **interruptor explícito** y declarar un locus **no** lo activa. El registro de la decisión, con su justificación, está junto a REQ-CONF-132 en la especificación | C-1 |
| `especificaciones-0.1.md` | REQ-CONF-88/89/94/95, REQ-TEC-07/08, REQ-DAT-46/47/52/53, REQ-MAN-45/48/60 | Las redacciones de la sección 9.1 del diseño con PCR, **con el cambio de que el parámetro declarado es la fracción de duplicados y no las parejas por molécula** | S-19, diseño §9.1 |
| `especificaciones-0.1.md` | REQ-REF-19 y los de referencia derivada | Admitir regiones declaradas en la referencia simulada y la referencia derivada por inserción de tracto, con su mapa de correspondencia | S-03 |
| `especificaciones-0.1.md` | §4.7, §4.6, §7 | Añadir REQ-MAN-93/94, REQ-VAL-57/58, ACE-118 y ACE-119; reescribir ACE-49/51/68 | Ambos documentos |
| `decisiones-tecnicas-0.1.md` | Tabla | **DTC-090** (provisional): `stutter-v1` con el predictor de la **tirada más larga** y la estructura asimétrica del paso, catálogo `f4-technical-v4`, y la aclaración de que DTC-067 no prohíbe usar `bernoulli-v1` ni `categorical-f64-v1` aquí | S-07, S-11, S-14 |
| `decisiones-tecnicas-0.1.md` | Tabla | **DTC-091** (provisional): `pcr-amplification-v1`, la programación dinámica de los primeros ciclos, la Poisson pinzada, el rechazo con tope, la genealogía, las constantes y su juez, la derivación de la probabilidad por ciclo; y **el método de inversión numérica de la fracción de duplicados: bisección sobre `[0, 64]` con 64 iteraciones fijas, exponencial exclusivamente de `libm = =0.2.16` conforme a DTC-067 y DTC-072, retorno anticipado exacto para fracción cero y rechazo explícito de la fracción inalcanzable** | S-19 a S-28 |
| `decisiones-tecnicas-0.1.md` | Tabla | Decisión sobre la realización del tracto: pureza como salida, alelo por secuencia, guarda de admisión con mínimo de pureza y de tirada | S-02, S-06 |
| `decisiones-tecnicas-0.1.md` | DTC-089 | Anotar que el presupuesto de bytes de salida cuenta ahora las bases del tracto publicadas en la verdad y en el VCF | S-02 |
| `plan-implementacion-f6-0.1.md` | §8, F6.4.3 | Corregir: **ya no depende** del primer paso del hito de variantes estructurales | S-02 |
| `plan-implementacion-f6-0.1.md` | §8, F6.4.2 | Corregir: **no depende** de derivar una referencia estructural | S-03 |
| `plan-implementacion-f6-0.1.md` | §8, F6.4.1 y F6.4.2 | Corregir el criterio de prueba: la pureza no se declara, y una realización con pureza 1 no es el caso que discrimina | S-06 |
| `plan-implementacion-f6-0.1.md` | §3 | DP-F6-03, DP-F6-06 y DP-F6-07 quedan cerradas **en su parte de repeticiones**; DP-F6-04 anota que el catálogo técnico nuevo lleva las etapas `stutter` y `amplification`, y que hay ámbitos nuevos para el sorteo por frecuencias y para el deslizamiento biológico | S-06/08/09/11, S-04 |
| `plan-implementacion-f6-0.1.md` | §8 | Añadir un paso F6.4.10, «Tartamudeo heredable con PCR», que remite a la entrega 4 de este plan; y anotar que F6.4.4 se adelanta en la entrega 5 | S-18, S-29 |
| `plan-implementacion-f4-0.1.md` | Sección 16 nueva | El cambio de contrato de la sección 10 de este plan | S-18 |
| `plan-implementacion-f1-0.1.md` | Sección 14 nueva | El cambio de contrato de la fase del núcleo reproducible: generador con regiones y mapa de coordenadas | S-03, S-18 |
| `hoja-de-ruta-0.1.md` | Fases del núcleo reproducible, de realismo técnico y de variación compleja | Registrar el hito propio, con fecha, y la reapertura acotada de los dos contratos | S-18 |
| `registro-ampliaciones-futuras-0.1.md` | DTF-056 | Dividir según §10.2; añadir **DTF-059**: errores de polimerasa heredables, quimeras de PCR heredadas, sesgo de longitud y genealogía exacta | Diseño §9.4 |
| `calibracion-tartamudeo-0.1.md` | §2, §3, §3.1, §4, §5, §6, §7 y §8 | **Aplicado ya**: los dos conjuntos comparados **no** son las mismas células sino dos individuos distintos, y se retira el residuo que decía que la PCR era lo único que los distinguía; la fuente de la dependencia con la longitud del tracto **no** está tras muro de pago y su licencia es **CC BY-NC**; la forma del modelo pasa a citar a Raz y col. (2019) junto a la de tres parámetros; y la reducción media se corrige de 4,63× a **4,64×** | Hechos falsos 3 y 4 de §1.1, y menores |
| `diseno-bibliotecas-pcr-0.1.md` | §5.2, §5.6, §5.9, §11 decisión 1 | **Propuesto**: el parámetro declarado es la **fracción de duplicados**, no las parejas por molécula; añadir la inversión numérica y su contrato de determinismo; corregir la licencia y la interpretación de las cifras de duplicados | S-19, hecho falso 5 |
| `diseno-bibliotecas-pcr-0.1.md` | §1, §5.7 | **Propuesto**: corregir «los dos conjuntos pasaron por las mismas células», que se copió del documento de calibración | Hecho falso 3 |
| `expediente-decisiones-repeticiones-tandem-0.1.md` | D-02, D-04, D-05, D-07, D-09, D-14, D-16, D-17 | **Propuesto**: anotar en cada ficha la decisión que la corrigió, sin reescribirlas, remitiendo a la sección «Decisiones confirmadas» de este plan | Todas las correcciones |
| `expediente-decisiones-repeticiones-tandem-0.1.md` | Parte B, encabezado | **Propuesto**: identificar **cuáles** son las tres confirmaciones; hoy solo D-10 lleva la marca | Aviso de recuento, §2 |
| `trazabilidad-requisitos-0.1.md` | Matriz | Filas nuevas para los requisitos nuevos; las que hoy apuntan a DTF-056 pasan a la sección 16 de la fase de realismo técnico | Ambos documentos |

---

## 12. Matriz de cierre

| Entrega | Resultado construido | Evidencia | Riesgo protegido |
|---|---|---|---|
| **1** | Referencia con tractos y su mapa de coordenadas; loci con alelos que llevan **su secuencia**; sorteo por frecuencias; verdad v3; VCF explícito sin umbral; presupuesto de salida ajustado; interfaz | Dos inserciones seguidas con el mapa correcto; tractos impuros con tiradas distintas; guarda de admisión con su mínimo de pureza; dos alelos de igual longitud y distinta interrupción distinguibles; frecuencia realizada distinta de la declarada en corpus pequeño; hashes sin loci iguales a la línea base | Que el alelo se reduzca a un contador y pierda las interrupciones; que la guarda rechace tractos reales o no rechace nada; que el umbral de 50 pb sobreviva; que el sorteo dependa del orden de los individuos |
| **2** | Tartamudeo de base después de la quimera, por molécula, en la tirada más larga, con tasas calibradas | Fixture de total igualado y tirada distinta; fixture saturado con el tracto en el solapamiento; corpus con quimera alta sin error de comparación; locus con alelo de referencia a frecuencia 1 | Que el predictor sea el número total de unidades; que el par salga incoherente; que el tartamudeo cambie la verdad; que la congelación se lea como absoluta |
| **3** | Familias de amplificación, duplicados, composición al elegir la molécula, **fracción de duplicados invertida numéricamente**, cota publicada | Inversión comparada con el cálculo directo a mano y con contador de iteraciones constante; duplicados contados por un oráculo que solo lee el FASTQ; recuento de aceptadas al agotar el tope mayor que cero; identidad byte a byte con fracción cero | Que la inversión no sea reproducible entre plataformas; que la composición se aplique al tamaño de familia; que se publique lo esperado en vez de lo realizado; que el rechazo no tenga tope |
| **4** | Tartamudeo heredable por linaje, y la fracción compuesta que reproduce lo medido | Fixture con divergencia forzada; probabilidades por ciclo contra el cálculo a mano; fracción observada con tamaño muestral calculado (≈2700 y ≈3000 moléculas) | Que el tartamudeo clonal se sortee por pareja; que se sumen las dos columnas; que la biblioteca amplificada subestime entre un 14 % y un 28 % en motivos de 1 a 5 pb |
| **5** | Servicio escalonado con su regla de agotamiento; generador de alelos por deslizamiento biológico con ámbito propio; **cierre de la simulación de microsatélites**. Bloqueada en su algoritmo por **P-4** | Recuento exacto de eventos no aplicados al agotar el tracto; vectores congelados del ámbito nuevo y **aserción de que su subsemilla difiere de la del tartamudeo**; origen de cada alelo publicado; recorrido completo de opciones con control en la interfaz | Que la biología y el artefacto compartan ámbito de azar, de modo que un parámetro de laboratorio decidiera el genotipo real; que el servicio sature en silencio al agotar el tracto; que quede una opción configurable sin control |

---

## 13. Lo que queda fuera a propósito, con destino

| Queda fuera | Destino |
|---|---|
| Duplicaciones, deleciones estructurales, inversiones y translocaciones | Hitos de referencias estructurales y de variantes estructurales de la fase de variación compleja. **No entran tampoco en el generador de referencia con regiones**: van después de generar la secuencia |
| **Microsatélites nuevos** que aparecen por inserción o duplicación | Trabajo estructural futuro. Este plan **varía tractos que ya existen y están declarados** |
| Ploidía generalizada por secuencia | Hito de ploidía generalizada |
| El modelo escalonado aplicado **generación tras generación** | Fase de evolución. La entrega 5 hace **una sola tirada** para inicializar la diversidad y no invade eso |
| Mutaciones internas de pureza, interrupciones o motivo | Paso F6.4.5. REQ-CONF-63 las exige desactivadas por omisión |
| Generación aleatoria de loci por regiones | Paso F6.4.6 |
| **Localizar los microsatélites que un genoma suministrado ya contiene**, en vez de exigir que se declaren | **Hito F6.4.a**, anotado en la hoja de ruta el 24/09/2026 a instancia del responsable. Este plan **comprueba** que un locus declarado existe de verdad donde se dice (F6.4.2), pero no descubre ninguno. Con descubrimiento automático aparecería una clase de locus que el usuario no pidió, y habría que decidir su régimen de tartamudeo por omisión |
| **Mínimos de admisión ajustados por longitud de motivo**, en vez de declarados locus a locus | Registrado por el responsable el 20/09/2026 junto a la decisión S-31. Hoy el suelo se declara por locus, con el general como valor de partida; derivarlo de la longitud del motivo es una ampliación posible y **no se hace** |
| Asignación **explícita de individuos concretos** con combinaciones de alelos específicas | Después de este trabajo; la decisión S-04 la invirtió respecto de las frecuencias, que son lo primero |
| Fichero complementario del VCF y catálogo de códigos de pérdida | Hito de variantes complejas. No hace falta: el alelo explícito con su secuencia no pierde nada |
| Conjuntos y matrices por referencia de mapeo, BAM, CRAM y mapeo real | Hito de conjuntos y matrices por referencia |
| Errores de polimerasa heredables por familia y quimeras de PCR heredadas | DTF-059 nuevo. Kebschull y Zador (2015) midieron que los primeros apenas afectan a la distribución y que las segundas son raras |
| Sesgo de longitud de la amplificación y genealogía exacta del proceso de ramificación | DTF-059 |
| Duplicados ópticos y de agrupación, y marcado de duplicados en BAM/CRAM | DTF-056 y DTF-055, sin cambios |
| Perfil de sesgo GC **calibrado con datos reales** | DTF-056, que conserva esa parte |
| **Mosaicismo somático** modelado aparte | No se modela: se queda dentro de la tasa sin PCR, porque la medición no permite separarlo. Si algún día se modelara, pertenecería a la **verdad** y no a las lecturas, y por tanto no cabe en una fase cuyo contrato prohíbe que los eventos técnicos alteren la verdad |
| Una tabla de tasas **por tipo de motivo** —distinguiendo `A` de `C`, `AC` de `AG`— como la que usa la fuente del modelo | Queda como desajuste conocido y anotado. La calibración del proyecto es por longitud de motivo y no se inventa una tabla por tipo que nadie ha medido para este uso |
| **Digitalizar la curva** de la dependencia con la longitud de la tirada | Pendiente P-1. La fuente **sí** es accesible; falta extraer la curva de su figura y comprobar la licencia antes de reproducir nada |
| Imposición de la cota de recursos como rechazo previo | Fase de presupuesto de recursos. Aquí se **publica**, no se impone |
| Optimizar el coste del tope de intentos y el factor dos de las dos pasadas | Fase de rendimiento. Se declaran y no se optimizan |
| La deuda de `tests/test_process_adapter.py`, que aborta bajo la puerta canónica | Sigue registrada en el plan de la fase de población fundadora, sin destino asignado. Este plan no la toca ni la cita como puerta |

---

## 14. Cambio de contrato posterior al cierre — el origen de cada alelo de repetición (27/09/2026)

Esta sección existe porque la **entrega 1 ya estaba cerrada y publicada** cuando
la entrega 5 añadió el generador de alelos, y esa adición **cambia bytes** en
documentos que escenarios ya existentes venían publicando. No es un añadido
silencioso y no debe presentarse como tal. Tiene la forma de la sección 14 del
plan del núcleo reproducible y de la sección 16 del plan de realismo técnico.

**Quién lo decide.** El responsable, el 27 de septiembre de 2026, con la
medición delante. El implementador de la entrega 5 había planteado la tensión
entre dos instrucciones: una cláusula de reproducibilidad pedía que un escenario
que **no** use el generador publicara exactamente los mismos bytes, y el paso
F6.4.12 de este plan exige marcar el **origen de cada alelo** y ofrece el
escenario mixto como discriminador. Las dos no caben a la vez. **Manda el plan**,
que es la autoridad, y la cláusula queda acotada por esta sección.

### 14.1. Qué **no** se reabre

- **Ninguna subsemilla, ningún vector congelado y ningún catálogo de azar
  publicado.** El ámbito nuevo nace en la versión 10 del catálogo de ámbitos,
  que conserva literalmente las filas de la 9, y la versión del catálogo **no se
  publica en ningún documento**: solo la comparan dos aserciones de prueba.
- **Las huellas congeladas de `tests/pcr_frozen_identity.rs` no se mueven.** Son
  el ancla absoluta de la reproducibilidad del recorrido técnico y siguen
  valiendo lo que valían.
- **Ningún esquema publicado se edita.** La capacidad va en las versiones que ya
  la llevaban abiertas —general v4, efectivo v6 y manifiesto general v6—, que
  son las que este mismo trabajo estrenó.
- **Ninguna decisión de las entregas 1 a 4 se revoca**, ni se rehace ninguna de
  sus evidencias.
- **El significado de las frecuencias no cambia.** La declarada sigue siendo la
  que se fijó y la realizada la proporción que el sorteo dio.

### 14.2. Qué cambia, y por decisión de quién

| Qué | Por qué no contradice lo cerrado | Decisión |
|---|---|---|
| Cada alelo de repetición publica su **origen**, `declared` o `generated`, en la configuración efectiva y en el manifiesto | Va **por alelo y no por locus** porque un documento con loci de los dos tipos tiene que poder distinguirlos sin volver a la configuración. Un escenario de un solo tipo **no discrimina** un productor que marcara todos igual, y ese es justamente el discriminador que F6.4.12 exige | **S-04** y F6.4.12; responsable, 27/09/2026 |
| Un locus generado publica además su bloque `slippage` con parámetros, versiones, ámbito, las tres limitaciones y la contabilidad | El bloque **solo existe** cuando el abanico lo generó el simulador; un locus con alelos declarados no lo lleva y su documento no cambia por esto | **DTC-097** |
| La cota de bases de tracto se **recomputa** después de generar | Sale de los alelos, y los alelos acaban de nacer. Dejarla como estaba publicaría una cota falsa y el juez independiente encontraría lo realizado por encima de lo acotado | **DTC-097** |

### 14.3. Qué bytes se mueven exactamente, medido

No se razona: se midió el 27/09/2026 compilando el binario del commit anterior
en un árbol aparte extraído con `git archive`, ejecutando los mismos escenarios
con los dos binarios y comparando los ficheros publicados.

| Escenario | Diferencias |
|---|---|
| **Sin loci de repetición** | **Ninguna. Byte a byte idéntico**, todos los ficheros. Es lo que sostiene que el contrato se ha ampliado sin romperse |
| **Con loci de repetición declarados a mano** | **Dos campos `origin`, uno por alelo**, en la configuración efectiva y los mismos dos en el manifiesto; y, derivados de ellos, el `size_bytes` y el `sha256` del artefacto de la configuración efectiva dentro del manifiesto. **Nada más** |
| **Con un locus generado** | No aplica: ese escenario no existía antes de la entrega 5 |

**La medición se repitió entera el 27 de septiembre de 2026**, después de las
correcciones de la revisión de cierre, y **da lo mismo hasta el byte**. Conviene
decir por qué dos de esas correcciones no mueven nada publicado:

- **La regla del borde** —que un linaje detenido no se deje empujar por el tope
  ni por el suelo— solo cambia el resultado de un **locus generado**, y esos no
  existían antes de esta entrega.
- **El recuento de Poisson que ahora repite en vez de abortar** (DT-010) solo
  toma la rama nueva **donde la anterior devolvía error**, de modo que ninguna
  ejecución que antes terminara consume palabras distintas ni obtiene un
  resultado distinto. Lo que sí cambió es una **prueba** que había congelado la
  tasa del defecto —veintisiete fallos numéricos sobre un corpus de 4096
  semillas— y que ahora exige cero; el cambio y su motivo están en DT-010.

### 14.4. Dónde vive cada literal

| Sitio | Papel | Qué dice hoy |
|---|---|---|
| `src/application.rs`, `RepeatAlleleOrigin` | canónico | las dos únicas procedencias, `declared` y `generated`, y su texto publicado |
| `src/repeat_slippage.rs` | canónico | el servicio escalonado, el generador, el ámbito y las tres limitaciones publicadas |
| `src/results.rs`, `repeat_allele_json` | contrato | el alelo publicado, compartido por la configuración efectiva y el manifiesto para que los dos no puedan decir cosas distintas |
| `src/results.rs`, `validate_repeat_slippage` | juez | lo que el esquema no puede comprobar: que los pasos aplicados más los no aplicados sumen los sorteados, y que el origen del alelo concuerde con el del locus |
| `schemas/effective-config-v6.schema.json` y `schemas/manifest-general-v6.schema.json` | contrato | `origin` obligatorio en cada alelo y el bloque `slippage` opcional del locus |
| `catalogs/subseed-scopes-v10.tsv` | contrato | la fila del ámbito nuevo, con las de la versión 9 intactas |

### 14.5. Pruebas nuevas y su variante defectuosa

Las de la entrega 5 están en `tests/repeat_slippage.rs`, en las pruebas del
módulo y en `tests/subseed_derivation.rs`. Las que sostienen **esta sección** en
concreto:

| Prueba | Contra qué muere |
|---|---|
| Un escenario **mixto** marca cada alelo con su origen, en los dos documentos | un productor que marque todos los alelos con el mismo origen; un escenario de un solo tipo **no lo caza**, y así se dice en su docstring |
| El lector del manifiesto rechaza un origen inventado y uno que contradiga al del locus | borrar la comprobación de coherencia, o borrar el bloque del generador entero |
| El lector rechaza una contabilidad en la que los aplicados y los no aplicados no suman los sorteados | un servicio que pierda eventos por el camino o que los cuente dos veces |
| Las huellas congeladas del recorrido técnico siguen en su sitio | cualquier cambio que desplace la secuencia de azar de una ejecución sin microsatélites |

### 14.6. Fixtures que se reabren

Ninguno. Los fixtures de la entrega 1 siguen valiendo: sus aserciones miran
campos concretos y ninguna exige el conjunto cerrado de claves del alelo, de
modo que el marcador nuevo entra sin tocarlos. El conjunto cerrado sí lo exige
el **juez** del manifiesto, y ahí el campo está declarado.

### 14.7. Evidencia que hay que repetir

Al cerrar la entrega 5: la comparación de bytes de la sección 14.3, que es la
que sostiene la afirmación de que el contrato se ha ampliado sin romperse. Se
ejecutó el 27/09/2026 con el resultado que esa tabla recoge, y **debe repetirse
si alguna vez se toca el alelo publicado**, porque es la única comprobación que
distingue una ampliación acotada de una reescritura silenciosa.
