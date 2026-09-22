# Registro de deudas técnicas 0.1

## 1. Finalidad, y en qué se distingue de los otros registros

Registro único de **deficiencias conocidas del código y de las pruebas que ya
existen**. Cada entrada describe algo que **está mal o sin cubrir hoy**, no algo
que falte por construir.

Dónde va cada cosa, para que nadie vuelva a dudarlo:

| Si es… | Va a |
|---|---|
| Una capacidad funcional excluida de la versión 1 | `registro-ampliaciones-futuras-0.1.md` |
| Una decisión de representación o implementación | `decisiones-tecnicas-0.1.md` |
| **Una deficiencia del código o de las pruebas actuales** | **este documento** |
| Un paso pendiente de una fase concreta | el plan de esa fase |

Hasta ahora las deudas vivían **dispersas dentro de los planes de fase**, lo que
hacía que una deuda transversal no tuviera sitio y que nadie supiera cuántas hay.
Este registro las reúne.

Una entrada se cierra cuando la deficiencia deja de existir, **no** cuando se
documenta mejor.

---

## 2. Registro

| ID | Deuda | Gravedad | Detectada | Estado |
|---|---|---|---|---|
| DT-001 | Dos tercios de las comprobaciones de los verificadores independientes no las ejercita ninguna prueba | **Grave** | 21/09/2026 | Abierta |
| DT-002 | Dos pruebas ignoradas están vaciadas, no solo apagadas | Media | 21/09/2026 | Abierta |
| DT-003 | Ocho pruebas ignoradas corresponden a contratos que ya no existen | Baja | 20/09/2026 | Abierta |
| DT-004 | La puerta de Python no corre desde el 25/08/2026 | Media | 19/09/2026 | Abierta |

---

## DT-001. Los verificadores independientes están casi sin probar

**Gravedad: grave.** Es la deuda más importante del proyecto hoy.

### Qué pasa

Los tres módulos que actúan como **jueces independientes** —los que permiten
afirmar que la verdad publicada es correcta— tienen la mayoría de sus
comprobaciones **sin ejercitar por ninguna prueba**:

| Verificador | Comprobaciones | Nunca se disparan |
|---|---|---|
| Oráculo de aceptación del paquete de haplotipos (`src/f5_oracle.rs`) | 366 | **253 (69 %)** |
| Oráculo semántico del formato VCF (`src/vcf_semantic_validation.rs`) | 58 | **36 (62 %)** |
| Auditoría de identidad y publicación (`src/f5_closure.rs`) | 42 | 11 (26 %) |

De las 253 del primero, 21 son autocomprobaciones contra catálogos compilados
—inalcanzables a propósito— y 4 están declaradas en el código como defensa
inalcanzable. **Quedan 228 comprobaciones de contrato ordinarias que nadie ha
probado nunca.**

### Por qué importa más que un fallo corriente

Estos módulos **no producen datos: los verifican**. Si una de esas
comprobaciones estuviera mal escrita —invertida, mirando la columna equivocada,
comparando algo consigo mismo—, la batería entera no lo notaría. El resultado no
sería un error visible, sino **un dato incorrecto que alguien se cree**.

Las peor cubiertas son precisamente las que vigilan la verdad publicada:

| Función | Puntos de rechazo | Muertos |
|---|---|---|
| `recount_f5_man70`, recuento independiente desde la verdad | 30 | **30** |
| `read_truth_copy`, validación de cada copia haplotípica | 26 | **26** |
| `recount_genotype`, validación de cada genotipo | 25 | **24** |
| `read_truth_header`, cabecera de la verdad | 13 | **13** |

Y en el verificador del formato VCF está muerta **la comprobación que ata la
columna de referencia al genoma real** (`validate_reference_allele`, el rechazo
«REF no coincide con el FASTA autoritativo»). Verificado de forma independiente:
esa frase **solo aparece en su propio módulo**; ninguna prueba la provoca.
También están muertas las dos direcciones de «el VCF dice exactamente lo que
dice la verdad, ni más ni menos», y **todo el contrato de cabecera del VCF**.

Además, el mecanismo que impone que cada unidad bloque×subpoblación aparezca
**exactamente una vez y en orden ascendente** en cuatro de los siete ficheros
publicados tiene sus cuatro rechazos muertos. La batería de corrupciones sí
intercambia filas, pero siempre **dentro de la misma unidad**, que es justo lo
que ese mecanismo no vigila.

### Cómo se midió, y cómo rehacerlo

No se leyeron siete mil líneas buscando sospechas: **se midió**. Los tres
verificadores canalizan todos sus rechazos por dos o tres funciones
constructoras. Se les puso la anotación `#[track_caller]` de Rust —que hace que
la función sepa desde qué línea la han llamado— y una línea que registra esa
llamada. Con eso, **una sola ejecución de la batería dice qué guardas llegan a
dispararse y cuáles no**.

Después se confirmó borrándolas: desactivar **trece** comprobaciones a la vez
—entre ellas «el sitio cae fuera de las coordenadas de su bloque», «los pesos
haplotípicos no suman 1» y «el recuento independiente halla otro número de
copias»— deja la batería en **1198 pasan, 0 fallan**, cifra por cifra idéntica a
la línea base. Lo mismo con cuatro del verificador del VCF.

### Qué hacer

**Lo primero, y es barato: convertir el barrido en una puerta.** Una anotación,
una función de cinco líneas y una ejecución de la batería bastan para medir en
cada cambio qué fracción de las comprobaciones del verificador se ejercita de
verdad. Sin esa puerta, la cobertura volverá a caer sin que nadie lo note.

**Después, cubrir.** Conviene empezar por el **verificador del formato VCF**, que
es el más pequeño —58 comprobaciones— y contiene la que ata la referencia al
genoma real.

Un aviso del propio hallazgo: la auditoría de identidad y publicación es **el
módulo mejor probado de los tres**, y no por casualidad. Tiene tres pruebas que
leen su propio texto fuente **con control positivo**, y una de ellas —la que
prohíbe que el módulo contenga cualquier vía de escritura en disco— **cazó la
instrumentación del revisor** en cuanto la introdujo. Ese es el patrón a imitar.

### Lo que quedó sin medir

Las fases 3 y 4 al completo —salida y validación de lecturas, errores y
calidades, quimeras, contaminación, multiplexado, pares y muestreo—; la admisión
del productor; y el contrato de la primera versión de la verdad, que no tiene un
punto único de construcción del error y exigiría otro enfoque.

---

## DT-002. Dos pruebas ignoradas están vaciadas, no solo apagadas

**Gravedad: media.** En `src/execution.rs`, dos pruebas marcadas como ignoradas
llevan además `#[allow(unreachable_code)]` y un `return;` **a media función**.
Todo lo que viene detrás —las aserciones sobre el manifiesto fallido, el estado
del evento, las referencias descartadas— es código muerto **incluso si alguien
las reactivara**.

No son pruebas apagadas: son **pruebas vaciadas que conservan el nombre de lo
que decían comprobar**. Quien las encienda creerá haber recuperado cobertura que
no existe.

O se restauran enteras, o se borran. Conservarlas así es peor que no tenerlas.

---

## DT-003. Ocho pruebas ignoradas son de contratos que ya no existen

**Gravedad: baja.** De las once pruebas ignoradas del proyecto:

- **ocho** corresponden a contratos que la fase 4 dejó obsoletos al pasar de
  lecturas sueltas a lecturas en pareja: una en `src/publication.rs`, dos en
  `src/execution.rs` —las de DT-002—, una en `tests/f3_event_assembly.rs` y
  cuatro en `tests/f3_event_admission.rs`;
- **dos** son experimentos de rendimiento que se ejecutan a mano, en
  `tests/f3_wgs_sampling.rs` y `tests/f4_demultiplexing.rs`;
- **una** es un ayudante de subproceso, en `tests/publication_process.rs`.

Las tres últimas están bien así. Las ocho primeras **deberían borrarse**: una
prueba apagada porque su contrato desapareció no es inofensiva, porque nadie
sabe si lo es sin volver a averiguarlo — y en dos de los ocho casos la respuesta
resulta ser que **ni siquiera dirían nada si se encendieran**.

---

## DT-004. La puerta de Python no corre desde el 25 de agosto de 2026

**Gravedad: media.** `tests/test_process_adapter.py` importa `pytest`, que no
está declarado en `pyproject.toml`. Bajo la puerta canónica —`unittest`— el
módulo aborta al importar y se cuenta como error; solo corre invocándolo con la
dependencia añadida a mano.

Viene del commit `6139fe9`. Registrada también en el plan de la fase 5; se
recoge aquí porque es transversal. Quien la cierre debe decidir entre declarar
`pytest` como dependencia de desarrollo o reescribir el módulo con `unittest`.

---

## 3. Lección transversal

Las cuatro deudas comparten una forma: **algo que parece cubierto y no lo
está**. Una guarda que existe pero no se ejercita, una prueba que conserva su
nombre pero no su cuerpo, una prueba apagada que nadie sabe si importa, una
puerta que se cita en los informes y no se ejecuta.

De ahí la regla que el proyecto ya aplica a las regresiones y que conviene
extender a la verificación entera: **no basta con que exista la comprobación;
hay que demostrar que puede fallar**.
