# Revisión de F1.3

- Fecha: 9 de agosto de 2026.
- Alcance: C1–C8, sin iniciar F1.4 ni integrar aleatoriedad en la generación de
  referencia.
- Método: revisión completa del diff conforme a `revisar-codigo-rust`, seguida
  de correcciones conservadoras con `aplicar-correcciones-rust` y repetición de
  todas las puertas.

## Resultado por frente

- **Raíz y entropía:** `RootSeed` encapsula exactamente 32 bytes y solo analiza
  o emite 64 hexadecimales minúsculos. La entropía es inyectable, el CSPRNG se
  consulta una vez solo cuando falta la raíz y su fallo conserva categoría
  `RES`, código y salida 3. La configuración efectiva registra valor y origen.
- **Derivación y aislamiento:** el material usa raíz cruda, versión y número de
  componentes `u32` big-endian, nombres ASCII normativos y valores UTF-8 con
  longitud `u32` big-endian. Los cuatro vectores A5 pasan contra un oráculo
  independiente. Cada subsemilla crea un `ChaCha12Rng` local no clonable; añadir,
  retirar o reordenar consumidores ajenos no cambia el flujo observado.
- **Metadatos y comparación:** raíz, origen, ChaCha12, versión de derivación,
  catálogo y subsemillas de evento quedan en un contrato propio inmutable. Los
  ámbitos internos se reconstruyen desde IDs. Los cuatro comparadores devuelven
  estado y condiciones por separado; igualdad semántica no implica bytes y la
  identidad binaria nunca se infiere por defecto.
- **Dependencias:** versiones exactas y rasgos mínimos constan en el lock y en
  la evaluación C1. La revisión de licencias, código nativo, scripts, árbol y
  RustSec no deja hallazgos pendientes.

## Hallazgos y correcciones

| ID | Prioridad | Decisión inicial | Corrección y regresión | Estado |
| --- | --- | --- | --- | --- |
| R-F13-01 | media | Aplicar | Un artefacto presente en una observación y ausente en la otra se clasificaba como `no aplicable`; se corrigió a `no garantizado`, reservando `no aplicable` para ausencia en ambas, y `TEST-REQ-VAL-44` cubre ambos casos. | Corregido |
| R-F13-02 | baja | Aplicar | Clippy detectó un `must_use` redundante sobre un iterador ya marcado; se retiró y la puerta `-D warnings` volvió a pasar. | Corregido |
| R-F13-03 | baja | Aplicar | La evaluación C1 contenía un error tipográfico que oscurecía el inventario de transitivas; se corrigió y se añadió el resultado final de árbol, auditoría y vectores. | Corregido |
| R-F13-04 | media | Aplicar | Una firma de componentes imposible para un ámbito conocido reutilizaba el diagnóstico de ámbito desconocido; ahora se trata como invariante interna y se conserva `GPS-CONF-018` exclusivamente para un ámbito realmente ajeno al catálogo. | Corregido |
| R-F13-05 | media | Aplicar | Los comparadores trataban el orden físico de dependencias, herramientas y subsemillas como significativo; ahora comparan multiconjuntos sin orden y las regresiones reordenan dependencias y herramientas sin retirar la garantía. | Corregido |

Balance: **5 confirmados, 5 corregidos, 0 descartados y 0 pendientes**.

## Revisión independiente H-F13

Los tres hallazgos se verificaron literalmente antes de corregir. Para
H-F13-01 se añadió primero una regresión que falló: cambiar
`rand_chacha=0.10.0` por `rand_chacha=changed` todavía producía
`Guaranteed` en el nivel biológico.

| ID | Gravedad | Decisión inicial | Estado final | Corrección y regresión |
| --- | --- | --- | --- | --- |
| H-F13-01 | grave/alta | Aplicar | Aplicado | La base común compara sin orden las dependencias RNG y exige algoritmo, versión de derivación y catálogo. Una regresión cambia individualmente rand_chacha, BLAKE3, derivación y catálogo; todos retiran la garantía, mientras el mero reordenamiento no. Las dependencias exclusivamente técnicas quedan en otro campo y nivel. |
| H-F13-02 | grave/alta | Aplicar con adaptación | Aplicado con adaptación | Se reauditaron literalmente los diez REQ: solo REQ-REP-23 está completamente demostrado; 3 pasan a F1.5 y 6 a F4.6/F9.7. C6 se estrecha a metadatos RNG en `EffectiveConfig` y C8 a la auditoría de infraestructura. Las nueve pruebas sintéticas se renombran y conservan como evidencia parcial, sin afirmar manifiesto, biología, lecturas ni artefactos. |
| H-F13-03 | grave, incluido funcionalmente en H-F13-02 | Aplicar con adaptación | Aplicado con adaptación | El booleano se sustituyó por `BinaryIdentityEvidence`, que identifica clase de artefacto, política e ID de evidencia no vacíos. El comparador sigue siendo preparatorio y REQ-CONF-115 queda reasignado a F9.7 hasta existir verificación sobre artefactos reales. |
| H-F13-02-R1 | grave residual | Aplicar | Aplicado | C3 todavía afirmaba que F1.3 introducía la raíz en un manifiesto y sugería una prueba CLI del fallo CSPRNG inexistente. El paso limita el resultado a `EffectiveConfig`, difiere serialización/publicación a D7/E8/F1.5 y distingue el doble fallido de aplicación del rechazo CLI de semilla inválida. |

Balance de la revisión independiente: **4 confirmados, 4 aplicados, 0
disputados y 0 pendientes**. Balance acumulado de F1.3: **9 confirmados, 9
corregidos o aplicados con adaptación, 0 pendientes**.

La inspección adversarial final no encontró dependencia de ruta, reloj, hilo,
mapa, posición accidental o planificador; tampoco `StdRng`, flujo global o
compartido, truncado de raíces, prioridad de consumo ni inferencia entre
niveles. La pausa de F1.3 queda alcanzada y pendiente de una nueva aprobación
independiente.
