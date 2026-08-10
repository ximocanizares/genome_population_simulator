# Revisión Rust y aplicación de correcciones de F1.4

Fecha: 10 de agosto de 2026.

## Alcance y fuentes de autoridad

Se revisaron el diff y los archivos nuevos de D1–D13: dependencias y lock,
tipos de identidad, SHA-256, esquemas y escritores JSON, publicación atómica,
integración CLI, pruebas, trazabilidad y evidencias. Se contrastaron
`especificaciones-0.1.md`, `arquitectura-0.1.md`,
`decisiones-tecnicas-0.1.md`, `trazabilidad-requisitos-0.1.md`, la hoja de ruta,
el registro de ampliaciones y el plan F1 aprobado. No se revisó ni modificó el
plan F2 ajeno al alcance.

## Veredicto

Sin hallazgos pendientes en F1.4. Los diez hallazgos de la implementación y las
revisiones fueron corregidos, recibieron prueba de regresión y superaron de
nuevo las puertas completas. La repetición posterior en GitHub Actions sobre
`ubuntu-24.04` cerró D12, D13 y la pausa obligatoria; F1.5 queda autorizado.

## Preguntas de revisión

1. ¿Puede un ID estable depender de ruta, reloj, hilo u orden accidental? No:
   el ID de ejecución se prueba con dos rutas distintas y materializa un digest
   encuadrado de los bytes admitidos, ordenados por clase e ID.
2. ¿Puede aparecer un final completo con un artefacto abierto, truncado o sin
   manifiesto? No: estados monotónicos, catálogo terminal, rechazo de
   `*.gps-tmp` y requisito de manifiesto general validado bloquean la
   publicación.
3. ¿Se conserva DTC-049 sin pérdida numérica? Sí: el manifiesto se vuelve a
   parsear y se compara mediante `to_bits()` para todos los `f64` exigidos.

## Hallazgos y correcciones aplicadas

| ID | Severidad | Hallazgo verificable | Corrección | Regresión/validación | Estado |
| --- | --- | --- | --- | --- | --- |
| H-F14-01 | Alta | `serde_json` solo con `std` recuperaba un vector DTC-049 a dos ULP del original. | Se aprobó y activó `float_roundtrip`; se actualizó la evaluación D1. | `manifest_preserves_every_dtc049_value_and_exact_f64_roundtrip` compara `to_bits()`. | Corregido |
| H-F14-02 | Alta | El primer ID de ejecución usaba el hash de `effective-config.json`, que contiene rutas; además la copia releía archivos mutables tras la admisión. | `EffectiveConfig` conserva los bytes admitidos; el digest usa encuadre canónico de clase, ID y bytes, nunca ruta; las copias consumen la instantánea. | `admitted_configuration_bytes_are_frozen_and_execution_identity_ignores_source_paths`. | Corregido |
| H-F14-03 | Alta | `ExecutionWorkspace::publish` podía ser llamado sin manifiesto general catalogado y no detectaba un temporal truncado residual. | La puerta exige un `manifest_general` validado y rechaza recursivamente cualquier `*.gps-tmp`. | `truncated_temporary_is_never_cataloged_or_published` y casos completo/fallido. | Corregido |
| H-F14-04 | Media | Las APIs de estado permitían invalidar un artefacto ya catalogado o modificar el catálogo de una unidad completa. | Transiciones terminales y monotónicas; completo/fallido no admiten cambios posteriores. | `artifact_and_unit_state_machines_reject_every_premature_completion`. | Corregido |
| H-F14-05 | Media | Las entradas de configuración del manifiesto se reordenaban por ruta física. | Orden estable por clase e ID de la copia; IDs de escenario y eventos se validan antes del provisional. | escenario real, rutas distintas y JSON semántico estable. | Corregido |
| H-F14-06 | Baja | El CLI contenía un `expect` alcanzable solo si se rompía la invariante del parser, contrario al criterio de no usar panic en fronteras. | Se sustituyó por `GPS-INTERNAL-001` tipado con salida 70. | Clippy, suite CLI y búsqueda de `panic`/`expect` en producción. | Corregido |
| H-F14-R1 | Alta | Decisión inicial: aplicar. `ExecutionWorkspace::publish` retiraba `INCOMPLETE` antes de entrar en `publish_no_replace`; un SIGKILL en esa ventana dejaba un provisional sin marcador. | El marcador permanece dentro del directorio durante `renameat2`; solo se retira desde el destino publicado tras el éxito. Ante fallo se preserva o restaura, y cualquier fallo de retirada/restauración se propaga con el invariante incompleto visible. No se añadió comprobación previa del destino ni fallback. | La regresión falló inicialmente al observar marcador ausente. `incomplete_marker_is_present_during_atomic_publication_and_survives_publisher_failure` y `marker_restoration_and_removal_failures_preserve_the_incomplete_invariant` pasan tras la corrección. | Corregido |
| H-F14-R2 | Alta | Decisión inicial: aplicar. Los objetos y elementos anidados de los tres esquemas aceptaban propiedades inesperadas y apenas restringían tipos o campos requeridos. | Se definieron contratos Draft 2020-12 cerrados para programa, plataforma, toolchain, dependencias, comprobaciones, entradas, eventos, salidas, unidades, diagnósticos, frecuencias controladas, valores efectivos, selectores y procedencia. Las colecciones aún futuras quedan vacías. | La regresión falló inicialmente porque `program.unexpected_nested_field` era válido. `schemas_reject_unknown_fields_in_every_nested_f1_4_contract` usa `python-jsonschema` y rechaza propiedades nuevas, requeridos ausentes y tipos/rangos corruptos en los tres documentos. | Corregido |
| H-F14-R3 | Media/documental | Decisión inicial: aplicar. E13 y la pausa F1.5 conservaban los recuentos anteriores 24/66/26, en contradicción con la matriz vigente. | Se corrigieron únicamente esos residuos a 27 REQ de F1.5, 61 filas vigentes de F1 y 31 reasignaciones a F2–F9. | Búsqueda literal sin coincidencias de `24 REQ asignados a F1.5`, `66 filas asignadas`, `26 casos reasignados`, `16/1/25/24` o `16 + 1 + 25 + 24`. | Corregido |
| H-F14-R2-R1 | Alta | Decisión inicial: aplicar. El cierre sintáctico de H-F14-R2 aún aceptaba contradicciones: un manifiesto raíz `failed` sin diagnóstico, dos niveles `technical` sin `biological`, combinaciones raíz/unidad incompatibles y un catálogo de dependencias que omitía una entrada duplicando otra con distinta versión. | Los contratos Draft 2020-12 modelan únicamente los estados emitidos en F1.4: un nivel exacto de cada clase y con sus condiciones reales; catálogo exacto de las siete dependencias bloqueadas; ramas condicionales `complete`/`failed` con configuración efectiva, diagnósticos y unidades coherentes; y manifiesto de evento con su única unidad completa y exportable. | La regresión nueva falló primero porque el esquema anterior aceptó `complete`→`failed` manteniendo `diagnostics=[]`. `schemas_reject_semantically_contradictory_f1_4_documents` pasa después y hace que `python-jsonschema` rechace de forma independiente estados, diagnósticos, exportabilidad, códigos, nivel ausente/duplicado y dependencia ausente/duplicada o con versión distinta del lock. | Corregido |

Contabilidad: 10 hallazgos recibidos, 10 corregidos, 0 descartados, 0 aplazados y
0 pendientes. No se aplicaron ciegamente soluciones sugeridas: cada corrección
se contrastó con el contrato y se verificó sobre el código real.

## Riesgos y límites restantes

- La validación nativa Ubuntu 24.04 quedó registrada en la ejecución
  `31372725778`, job `93405082155`, para el commit exacto `5d12e58`.
- Ocho requisitos permanecen correctamente reasignados porque sus consumidores
  reales pertenecen a F1.5 o F2–F9.
- No se implementó referencia, FASTA, verdad, VCF, FASTQ, BAM/CRAM ni biología
  de F1.5 o posteriores.

## Cierre de la pausa

La ejecución
[`31372725778`](https://github.com/ximocanizares/genome_population_simulator/actions/runs/31372725778)
y su job
[`93405082155`](https://github.com/ximocanizares/genome_population_simulator/actions/runs/31372725778/job/93405082155)
repitieron satisfactoriamente la matriz de publicación en `ubuntu-24.04`.
D12–D13 y F1.4 quedan aprobados; se conservan las ocho reasignaciones con sus
destinos y se autoriza E1/F1.5.
