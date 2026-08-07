# Hoja de ruta del simulador de genomas y poblaciones

- **Versión:** 0.1
- **Estado:** borrador inicial para revisión
- **Fecha:** 7 de agosto de 2026
- **Especificación de origen:**
  [especificaciones-0.1.md](./especificaciones-0.1.md)
- **Lenguaje del núcleo:** Rust

## 1. Finalidad

Esta hoja de ruta divide el desarrollo en fases que produzcan resultados
verificables y mantengan controlado el crecimiento de la aplicación. No es
todavía un plan de implementación: no decide módulos, crates, formatos de
configuración, algoritmos concretos ni archivos de código.

Cada fase deberá transformarse en un plan técnico detallado únicamente cuando
la arquitectura y las decisiones de las que dependa estén acordadas. El plan de
una fase dividirá sus hitos en cambios pequeños, con pruebas y pausas de revisión.

## 2. Principios para controlar el crecimiento

1. **Mantener un recorrido completo.** Desde las primeras fases debe existir un
   camino pequeño pero funcional desde la configuración hasta resultados que
   puedan validarse.
2. **Ampliar una dimensión cada vez.** Añadir variantes, modelos poblacionales o
   realismo técnico sin introducir simultáneamente varias capacidades no
   relacionadas.
3. **Cerrar cada fase antes de apoyar otra sobre ella.** Una fase termina con
   pruebas, evidencias y una revisión explícita.
4. **No anticipar trabajo futuro.** Una interfaz podrá dejar una extensión
   prevista, pero no se implementará una abstracción sin un caso confirmado.
5. **Conservar la verdad conocida.** Toda ampliación debe mantener referencia,
   coordenadas, ploidía, fase, dosis, linaje y procedencia sin pérdidas.
6. **Separar verdad biológica y evidencia técnica.** FASTQ, BAM, CRAM, errores y
   contaminación nunca modificarán el genotipo verdadero.
7. **Medir pronto, optimizar con evidencia.** Cada fase observará tiempo y
   memoria en casos representativos. Las optimizaciones complejas requerirán
   pruebas de rendimiento reproducibles.
8. **Planificar solo la siguiente fase.** No se redactará de una vez el plan
   detallado de toda la aplicación, porque las decisiones posteriores dependen
   de lo aprendido en fases anteriores.

## 3. Niveles de planificación

- **Fase:** capacidad amplia del producto y resultado que permite avanzar.
- **Hito:** resultado verificable dentro de una fase.
- **Paso:** cambio pequeño de implementación, revisión y prueba independiente.

Esta hoja de ruta define fases e hitos. Los pasos se escribirán en los planes de
implementación posteriores.

## 4. Dependencias generales

```text
F0 Arquitectura
 └─ F1 Núcleo reproducible
     └─ F2 Verdad genómica mínima
         └─ F3 Recorrido técnico mínimo
             └─ F4 Realismo técnico de secuenciación
                 └─ F5 Variación y referencias complejas
                     └─ F6 Población inicial
                         └─ F7 Evolución poblacional
                             └─ F8 Familias y segregantes
                                 └─ F9 Escalabilidad y consolidación
```

El orden expresa la secuencia recomendada. Una tarea aislada de preparación
podrá adelantarse si no introduce contratos ni decisiones prematuras, pero una
fase no se considerará terminada hasta satisfacer todas sus dependencias.

## 5. Fases e hitos

### F0. Arquitectura y decisiones técnicas

**Objetivo:** convertir la especificación funcional en una arquitectura Rust
trazable y resolver las decisiones necesarias para planificar el primer código.

**Depende de:** especificación 0.1 y decisión de usar Rust.

**Hitos:**

- **F0.1. Arquitectura del sistema.** Definir componentes, responsabilidades,
  flujos, contratos, errores y dirección de dependencias.
- **F0.2. Modelo de dominio.** Definir la representación conceptual de
  referencias, coordenadas, segmentos, haplotipos, individuos, poblaciones,
  eventos, lecturas y resultados.
- **F0.3. Decisiones técnicas iniciales.** Elegir mediante evaluación las
  opciones necesarias para configuración, manifiestos, aleatoriedad,
  serialización, formatos genómicos, procesos externos y pruebas.
- **F0.4. Estrategia de rendimiento.** Definir escenarios pequeños y grandes,
  métricas, límites observables y pruebas iniciales de memoria y tiempo.
- **F0.5. Plan de la fase F1.** Dividir el núcleo reproducible en pasos pequeños
  con dependencias y pruebas concretas.

**Resultado de salida:** arquitectura 0.1, decisiones técnicas justificadas y
plan ejecutable de F1.

**Evidencia de cierre:** todos los requisitos necesarios para F1 tienen un
componente responsable y no quedan decisiones de diseño ocultas en su plan.

> Revisión de F0: aprobar la arquitectura y el plan de F1 antes de crear el
> núcleo de producción.

### F1. Núcleo reproducible

**Objetivo:** obtener una aplicación Rust mínima que valide una configuración,
cree una ejecución identificada y produzca resultados trazables y repetibles.

**Depende de:** F0.

**Hitos:**

- **F1.1. Base del proyecto.** Disponer del proyecto Cargo, comprobaciones de
  formato, análisis estático, pruebas y política inicial de dependencias.
- **F1.2. Configuración validada.** Leer la configuración mínima, aplicar valores
  predeterminados y rechazar estados incompatibles antes de simular.
- **F1.3. Aleatoriedad reproducible.** Registrar la semilla y demostrar la
  repetición de un caso mínimo con entradas y versión idénticas.
- **F1.4. Identidad y resultados.** Crear identificadores, directorios,
  configuración efectiva y manifiesto mínimo sin confundir referencias,
  muestras o eventos.
- **F1.5. Referencia mínima.** Aceptar una referencia pequeña simulada o
  suministrada y conservar su FASTA e identidad verificable.

**Resultado de salida:** ejecución mínima repetible con configuración, semilla,
referencia y manifiesto.

**Evidencia de cierre:** dos ejecuciones equivalentes producen el mismo
resultado lógico y una configuración inválida no deja una ejecución marcada
como completa.

> Revisión de F1: comprobar reproducibilidad, diagnósticos y trazabilidad antes
> de representar individuos.

### F2. Verdad genómica mínima

**Objetivo:** crear individuos directamente configurados y conservar una verdad
genotípica faseada respecto a una única referencia.

**Depende de:** F1.

**Alcance inicial deliberado:** referencia pequeña, una ploidía uniforme,
individuos directos, SNP e inserciones o deleciones pequeñas no ambiguas.

**Hitos:**

- **F2.1. Estado haplotípico.** Representar copias cromosómicas, fase, dosis y
  procedencia sin duplicar innecesariamente secuencias completas.
- **F2.2. Variantes pequeñas.** Aplicar SNP e *indels* exactos y aleatorios con
  coordenadas y orden definidos.
- **F2.3. Individuos directos.** Construir individuos desde la configuración y
  reconstruir la secuencia de cada haplotipo.
- **F2.4. Verdad estructurada.** Guardar el genotipo real y la información
  necesaria para reconstruirlo independientemente del VCF.
- **F2.5. VCF de verdad mínimo.** Producir un VCF multimuesta ordenado, derivado
  del estado real y sin campos inferidos de lecturas.

**Resultado de salida:** referencia, individuos faseados, verdad estructurada y
VCF verdadero coherentes.

**Evidencia de cierre:** casos sintéticos independientes reconstruyen las
secuencias esperadas y demuestran coordenadas, alelos, fase y dosis.

> Revisión de F2: validar la verdad biológica antes de generar una sola lectura.

### F3. Recorrido técnico mínimo

**Objetivo:** producir el primer conjunto pequeño útil para probar externamente
un detector de SNP.

**Depende de:** F2.

**Alcance inicial deliberado:** lecturas cortas tipo Illumina, un único perfil,
un evento, sin contaminación ni degradaciones avanzadas y minimap2 como
mapeador.

**Hitos:**

- **F3.1. Lecturas desde haplotipos.** Generar FASTQ con cobertura o número de
  lecturas y contribución proporcional a dosis y ploidía.
- **F3.2. Mapeo.** Ejecutar minimap2 contra la referencia elegida, registrar
  versión y parámetros y detectar cualquier fallo del proceso.
- **F3.3. Alineamiento final.** Producir BAM o CRAM ordenado e indexado con
  cabeceras, muestra y grupos de lectura coherentes.
- **F3.4. Validación interna.** Comparar lecturas y alineamientos con su origen
  sin inferir genotipos ni ejecutar programas detectores de SNP.
- **F3.5. Conjunto exportable mínimo.** Reunir referencia, verdad, FASTQ,
  alineamientos, índices, configuración y manifiesto.

**Resultado de salida:** recorrido completo:

```text
configuración → referencia → individuo → verdad → FASTQ → minimap2 → BAM/CRAM
```

**Evidencia de cierre:** el conjunto supera validaciones de formatos y su BAM o
CRAM refleja de forma coherente la verdad biológica conocida.

> Revisión de F3: usar el primer conjunto extremo a extremo antes de ampliar el
> modelo biológico o técnico.

### F4. Realismo técnico de secuenciación

**Objetivo:** incorporar pronto los perfiles, degradaciones y eventos técnicos
de WGS confirmados, sin alterar la verdad genotípica.

**Depende de:** F3.

**Alcance inicial deliberado:** usar individuos directos y variantes pequeñas
de F2. Las fases biológicas posteriores deberán demostrar que sus individuos y
poblaciones pueden alimentar este mismo recorrido sin crear otra implementación.

**Hitos:**

- **F4.1. Perfiles tecnológicos.** Completar Illumina *single-end* y
  *paired-end* y añadir PacBio HiFi en eventos separados.
- **F4.2. Cobertura, longitud, error y calidad.** Aplicar perfiles validados y
  sustituciones configurables por los ámbitos permitidos.
- **F4.3. Artefactos WGS.** Simular duplicados, adaptadores, sesgo GC, quimeras y
  anomalías *paired-end* incluidas en la primera versión.
- **F4.4. Contaminación biológica.** Sustituir lecturas según fuentes, unidades y
  matrices dirigidas compatibles.
- **F4.5. Demultiplexado técnico.** Representar reasignaciones y descartes por
  índices y conservar la procedencia sin duplicar secuencias.
- **F4.6. Modo directo.** Producir alineamientos sintéticos coherentes sin FASTQ
  y demostrar concordancia funcional con el modo realista donde corresponda.
- **F4.7. Calibración desde datos reales.** Estimar perfiles reutilizables desde
  FASTQ y referencia suministrados sin copiar las lecturas concretas.

**Resultado de salida:** eventos WGS configurables, trazables y separados por
tecnología, referencia y muestra, disponibles desde una etapa temprana del
proyecto.

**Evidencia de cierre:** validaciones internas cuantifican lo solicitado y lo
realizado, y demuestran que las degradaciones solo afectan a la evidencia
técnica.

> Revisión de F4: aprobar los perfiles y artefactos antes de ampliar las clases
> de variación y los modelos poblacionales.

### F5. Referencias y variación complejas

**Objetivo:** ampliar la verdad genómica manteniendo reconstrucción, fase y
separación estricta entre referencias.

**Depende de:** F4.

**Hitos:**

- **F5.1. Referencias estructurales.** Crear referencias derivadas mediante
  duplicaciones, deleciones, inversiones y translocaciones permitidas.
- **F5.2. Variantes estructurales individuales.** Representar cambios
  estructurales en haplotipos y conservar sus efectos finales.
- **F5.3. Ploidía general.** Admitir la ploidía confirmada por cromosoma o
  *scaffold* sin diferencias entre individuos.
- **F5.4. Repeticiones y elementos transponibles.** Añadir repeticiones en
  tándem y eventos configurados de elementos transponibles dentro del alcance.
- **F5.5. Variantes solapadas.** Resolver cambios secuenciales, variantes
  eliminadas y alelos finales complejos.
- **F5.6. Varios conjuntos de referencia.** Mapear las mismas muestras contra
  referencias diferentes y producir BAM/CRAM y VCF independientes.

**Resultado de salida:** individuos complejos reconstruibles y conjuntos
independientes por referencia.

**Evidencia de cierre:** pruebas con verdad conocida demuestran transformaciones,
límites, solapamientos y ausencia de mezcla entre referencias.

> Revisión de F5: confirmar que la complejidad estructural no rompe la verdad
> mínima ni el realismo técnico establecido en F4.

### F6. Población inicial

**Objetivo:** crear poblaciones iniciales reproducibles sin ejecutar todavía su
evolución a lo largo de generaciones.

**Depende de:** F5.

**Hitos:**

- **F6.1. Creación probabilística y controlada.** Generar alelos, genotipos,
  dosis y recuentos según el modo seleccionado.
- **F6.2. Frecuencias y Hardy-Weinberg.** Aplicar y comprobar los modelos
  confirmados para las ploidías admitidas.
- **F6.3. Endogamia inicial.** Aplicar F en los ámbitos permitidos y registrar
  valores esperados y realizados.
- **F6.4. Desequilibrio de ligamiento.** Crear bloques no solapados mediante
  haplotipos y frecuencias conjuntas, conservando la fase.
- **F6.5. Grupos y subpoblaciones.** Asignar identidades, pertenencias y
  parámetros iniciales sin ambigüedades.
- **F6.6. Fundadores.** Identificar individuos y copias haplotípicas fundadoras
  de forma estable.

**Resultado de salida:** población inicial completa, faseada y trazable, apta
para evolución o para producir eventos técnicos.

**Evidencia de cierre:** pruebas estadísticas y controladas validan recuentos,
frecuencias, fase, F, bloques y reproducibilidad. Al menos una población inicial
produce un evento técnico mediante el recorrido de F4 sin una ruta alternativa.

> Revisión de F6: aprobar el estado inicial antes de introducir transmisión y
> generaciones.

### F7. Evolución poblacional

**Objetivo:** simular generaciones discretas manteniendo verdad, linaje y orden
configurado de los procesos.

**Depende de:** F6.

**Hitos:**

- **F7.1. Reproducción y recombinación.** Formar descendencia, bivalentes y
  entrecruzamientos con fase y procedencia reconstruibles.
- **F7.2. Mutación heredable.** Introducir SNP, *indels*, variantes estructurales
  y repeticiones según sus modelos y momentos confirmados.
- **F7.3. Tamaño poblacional.** Gestionar N y Ne objetivos y realizados,
  contribuciones parentales y avisos.
- **F7.4. Selección y deriva.** Aplicar selección direccional de viabilidad y
  muestreo aleatorio en el orden configurado.
- **F7.5. Migración y subpoblaciones.** Ejecutar migraciones simultáneas,
  divisiones, fusiones, fundaciones y extinciones.
- **F7.6. Fases evolutivas.** Variar parámetros y encadenar episodios como
  panmixia, cuellos de botella o migración.
- **F7.7. Observaciones.** Muestrear generaciones sin alterar la trayectoria y
  guardar historia, genotipos y salidas opcionales.

**Resultado de salida:** poblaciones evolucionadas reproducibles con linaje
hasta los fundadores.

**Evidencia de cierre:** casos pequeños permiten reconstruir cada transmisión y
casos estadísticos comprueban que los procesos realizados concuerdan con sus
parámetros dentro de las tolerancias acordadas. Una muestra evolutiva puede
producir los resultados técnicos de F4 sin modificar su verdad genotípica.

> Revisión de F7: comprobar exactitud y consumo de recursos antes de añadir
> grafos de cruces complejos.

### F8. Familias y poblaciones segregantes

**Objetivo:** crear descendencia mediante grafos declarativos de cruces y
conservar todos los individuos y haplotipos intermedios.

**Depende de:** F7 para reproducción, recombinación, mutación y linaje.

**Hitos:**

- **F8.1. Grafo de cruces.** Validar parentales, ciclos, dependencias y reglas
  declarativas antes de generar descendencia.
- **F8.2. Cruces básicos.** Implementar autofecundación, cruzamiento sexual,
  F2 y retrocruce como composiciones verificables.
- **F8.3. Selección de descendientes.** Elegir por identificador, azar,
  pertenencia, genotipo, dosis o haplotipo mediante el catálogo común.
- **F8.4. Diseños complejos.** Expresar pedigrís, líneas de introgresión y MAGIC
  mediante el mismo modelo, sin código ejecutable en la configuración.
- **F8.5. Poblaciones intermedias.** Conservar sus genotipos y generar FASTQ,
  BAM/CRAM o VCF únicamente cuando se solicite.

**Resultado de salida:** familias y poblaciones segregantes reconstruibles desde
sus fundadores y cruces.

**Evidencia de cierre:** casos mendelianos y poliploides controlados demuestran
parentesco, segregación, fase, recombinación y selección de descendientes.

> Revisión de F8: validar los grafos, los cruces y su integración con los
> eventos técnicos antes de la consolidación.

### F9. Escalabilidad y consolidación

**Objetivo:** convertir las capacidades anteriores en una primera versión
robusta, distribuible y capaz de manejar poblaciones grandes dentro de los
recursos disponibles.

**Depende de:** todas las fases anteriores incluidas en la versión que vaya a
publicarse.

**Hitos:**

- **F9.1. Escenarios de rendimiento.** Fijar casos representativos y medir CPU,
  memoria, almacenamiento, individuos, generaciones y lecturas por segundo.
- **F9.2. Representación compacta.** Optimizar, con evidencia, haplotipos,
  segmentos, genotipos, historia y materialización de secuencias.
- **F9.3. Paralelismo reproducible.** Escalar por unidades independientes sin
  cambiar resultados con la misma semilla.
- **F9.4. Presupuesto de recursos.** Estimar memoria, almacenamiento y tiempo,
  emitir avisos y rechazar ejecuciones que superen límites configurados.
- **F9.5. Interoperabilidad.** Validar versiones, cabeceras, compresión, índices
  y compatibilidad de FASTA, FASTQ, BAM, CRAM y VCF.
- **F9.6. Exportación autocontenida.** Empaquetar eventos con referencias,
  índices, configuraciones, manifiestos y sumas de comprobación.
- **F9.7. Preparación de publicación.** Completar documentación, ejemplos,
  pruebas de regresión y proceso reproducible de construcción.

**Resultado de salida:** primera versión candidata a publicación.

**Evidencia de cierre:** la batería funcional, genómica, estadística, de
interoperabilidad y rendimiento supera los criterios acordados sin hallazgos
bloqueantes abiertos.

> Revisión de F9: decidir si la versión está lista para publicarse o necesita un
> ciclo adicional de corrección.

## 6. Matriz resumida de comprobación

| Fase | Resultado construido | Evidencia principal | Riesgo protegido |
|---|---|---|---|
| F0 | Arquitectura y plan de F1 | Trazabilidad completa y decisiones cerradas | Diseño improvisado durante la implementación |
| F1 | Ejecución mínima reproducible | Repetición con semilla y rechazo limpio | Resultados sin identidad o configuración efectiva |
| F2 | Verdad genómica mínima | Reconstrucción independiente de haplotipos y VCF | Genotipos, fase o coordenadas incorrectos |
| F3 | Primer recorrido extremo a extremo | Validación FASTQ, BAM/CRAM y procedencia | Alineamientos incoherentes con la verdad |
| F4 | Realismo técnico WGS | Comparación de parámetros solicitados y realizados | Alterar la verdad con artefactos técnicos |
| F5 | Variación y referencias complejas | Casos estructurales y aislamiento por referencia | Mezcla de referencias o pérdida de linaje |
| F6 | Población inicial | Recuentos y pruebas estadísticas reproducibles | Frecuencias, LD o F mal representados |
| F7 | Evolución poblacional | Reconstrucción de transmisiones y medidas realizadas | Orden evolutivo o genealogía incorrectos |
| F8 | Familias y segregantes | Casos mendelianos y grafos validados | Cruces o selección de descendientes erróneos |
| F9 | Versión consolidada | Pruebas completas y escenarios de rendimiento | Fallos de escala, distribución o interoperabilidad |

## 7. Reglas de avance

Una fase puede comenzar a planificarse cuando sus dependencias estén cerradas.
Puede comenzar a implementarse cuando disponga de arquitectura suficiente, un
plan aprobado y criterios de salida verificables.

Una fase no se cerrará si:

- existe un defecto bloqueante o grave sin resolución o destino aprobado;
- falta una prueba independiente para una invariante crítica;
- los archivos producidos no superan sus validaciones de formato;
- la documentación afirma un resultado que no se ha ejecutado;
- una decisión pendiente cambiaría los contratos de la fase siguiente.

Todo cierre termina con una pausa para revisar evidencias y decidir si se
continúa, se corrige o se modifica la hoja de ruta.

## 8. Trabajo posterior a esta hoja de ruta

El siguiente documento debe ser la arquitectura Rust 0.1 correspondiente a F0.
Después se resolverán las decisiones técnicas estrictamente necesarias para F1
y se escribirá su plan de implementación detallado.

Las fases F2–F9 no deben dividirse todavía en archivos, módulos, crates ni
commits. Su desglose se hará progresivamente con la evidencia obtenida en las
fases anteriores.
