# Punto de reanudación

## Estado actual — 24 de agosto de 2026 (F2.4 cerrada localmente, pendiente de commit)

- F2.3.1–F2.3.7 está cerrada: configuración de fundadores, genotipos HWE/`F`,
  fase multialélica, eventos explícitos, materialización temporal, estadísticas
  y controles Marimo locales. El cierre está publicado en `2348f37`.
- `founder_keys` se exige cuando faltan IDs explícitos; es única dentro de cada
  población y deriva IDs estables desde población×clave.
- La matriz multialélica conserva una representación preparada única de REF y
  ALT. Cubre normalización con REF=0 y ALT globalmente ausentes, tanto en
  matriz suministrada como reutilizada.
- `validate` y `run` comparten un preflight que materializa temporalmente;
  ambos descartan copias y `run` publica solamente el resumen estadístico.
- F2.3 superó revisión independiente y revisión adversarial final. Marimo
  permanece exclusivamente local; el gate reproducible es
  `MARIMO_SKIP_UPDATE_CHECK=1 uv run --frozen marimo check --strict ui/gps_app.py`.
- F2.4 fija el esquema externo v1 e independiente de `truth.jsonl` en
  `src/truth_schema.rs` y `schemas/truth-v1.schema.json`, incluyendo
  fragmentos vigentes, historial y alelos estructurales tipados REF/inserción/
  deleción. `src/truth_projection.rs` proyecta el estado F2.3 a registros
  ordenados y publica JSONL UTF-8 por flujo mediante el provisional F1.
- `src/truth_validation.rs` valida esquema, orden, unicidad, recuentos,
  digests y relaciones, y reconstruye copias y genotipos usando únicamente
  `truth.jsonl` y la referencia. La integración registra el artefacto con
  ruta, tamaño y SHA-256 en el manifiesto; Marimo presenta ese resumen sin
  cargar el JSONL completo. La evidencia de cierre está en `evidence/f2.4/cierre.md`.
- La matriz canónica es una entrada obligatoria de F2.4: se valida contra la
  referencia y exige digests tanto de matriz suministrada/reutilizada como de
  grupo, incluso cuando un locus queda todo-REF o la matriz vacía fue
  autorizada. Su identidad se deriva de la huella canónica preparada, separada
  de los `source_id` heterogéneos de sus filas. Una copia cuyo SNP está
  cubierto por deleción no recibe `<deleted>`: el caso indel tipado conserva
  sus portadores y dosis 0..P, el evento fuente, enlaces a fragmentos vigentes
  y los bytes materializados finales de la copia. La revisión independiente de
  F2.4.4 y las puertas de cierre F2.4 se han completado localmente; no se ha
  creado commit.
- Adaptación aprobada de F2.4.1: `truth.jsonl` v1 conserva además fragmentos
  vigentes de fuente o inserción, el historial cronológico con procedencia
  explícita o aleatoria y casos indel activos.
  Por ello representa inserciones anidadas, deleciones sobre bases insertadas
  y SNP eliminados sin restringir el estado válido de F2.3.
- Existe un cambio local ajeno y sin publicar en
  `tests_ui/test_process_adapter.py`; no incluirlo, revertirlo ni atribuirlo a
  F2.3 sin revisarlo explícitamente.

## Estado histórico — 21 de agosto de 2026 (cierre F2.2.3)

- F1 sigue cerrada en `e14705a`.
- F2.0, F2.1 y F2.2.1–F2.2.3 están implementadas, validadas y revisadas de forma
  independiente; este commit las publica por primera vez.
- F2.1 aporta coordenadas de variantes pequeñas, ploidía 1..10, estado
  haplotípico por segmentos compartibles, materialización bajo demanda y dosis
  o fase derivadas del estado canónico.
- F2.2.1–F2.2.2 fijan la matriz SNV TSV+JSON, modos Marimo/TOML, normalización
  Neumaier y tolerancia efectiva. F2.2.3 incorpora RNG v4: v1 se conserva para
  F1 y v2 añade el nombre de ámbito para F2, con `scenario_id = "direct"`.
- `libm = =0.2.16` fue evaluada conforme DTC-043 y queda fijada sin *features*
  predeterminadas. Los vectores, consumos, PTRS/Poisson y sus oráculos quedan
  cubiertos por pruebas de regresión.
- La re-revisión final de F2.2.3 no encontró hallazgos. El siguiente paso es
  F2.2.4: validar y preparar matrices suministradas; no iniciar F2.2.5 ni
  generación hasta cerrar ese paso.

## Contexto histórico de simplificación

- F1 está cerrada en `e14705a` y su cierre documental en `efca265`.
- El plan F2 0.1 aprobado en `b64767b` resultó sobrediseñado al comenzar
  F2.0.1: creó un ejecutor general de evidencia y las revisiones ampliaron el
  modelo de amenaza a manipulación local activa.
- La persona usuaria ordenó simplificar toda F2 con el mismo criterio aplicado a
  F1.
- El agente de corrección quedó interrumpido antes de continuar el bucle.
- La implementación no aprobada de `verify-f2-step`, su esquema, fixtures,
  pruebas y `evidence/f2/F2.0.1` se retiró mediante `apply_patch`; nunca se
  comprometió.
- Se conservan los ajustes de trazabilidad F1→F2.2/F6.2 en
  `tests/f1_2_traceability.rs` y
  `evidence/f1.2/test-req-matrix.tsv`.

## Nueva autoridad F2

La simplificación está documentada en:

- `design/especificaciones-0.1.md`;
- `design/arquitectura-0.1.md`, versión 0.2;
- DTC-079 en `design/decisiones-tecnicas-0.1.md`;
- `design/hoja-de-ruta-0.1.md`;
- `design/plan-implementacion-f2-0.1.md`, ahora versión 0.2;
- `design/trazabilidad-requisitos-0.1.md`.

El plan pasa de 130 a 41 pasos y conserva F2.0–F2.5 y sus seis pausas.

## Criterio acordado

- Uso exclusivamente local en `127.0.0.1`, sin servidor remoto.
- Repositorio, persona usuaria y provisional privado confiables frente a
  manipulación activa.
- Entradas externas validadas; protección frente a fallos propios, interrupción
  y corrupción accidental.
- Marimo solo configura, invoca `validate`/`run` y muestra diagnósticos,
  estadísticas y enlaces.
- F2 usa un único TOML general y el ID interno `scenario_id = "direct"`; los
  escenarios se implementan con su primer consumidor real.
- Sin ejecutor de evidencia general, protocolo paralelo de capacidades,
  almacenamiento persistente de sesiones ni defensas frente a symlinks,
  inyección deliberada del entorno o carreras hostiles.
- Distribución `.tar.gz`/OCI, plataforma oficial, accesibilidad de entrega y
  endurecimiento quedan en F9.
- El aprovisionamiento inicial explícito con `uv sync --frozen` puede descargar
  las dependencias resueltas por `uv.lock`; el arranque posterior y el uso
  normal son locales y sin red. La instalación distribuible offline queda en
  F9.
- Se conservan las decisiones biológicas: matriz SNV común, perfiles
  aproximados, HWE/`F`, estado haplotípico, `truth.jsonl` y VCF 4.5.

## Siguiente acción

1. Revisar la pausa F2.4 y confirmar el cierre local antes de iniciar F2.5.
   No comenzar VCF hasta esa aprobación.

La revisión directa de simplificación comprobó 41 pasos bien formados, seis
pausas, 78 REQ exclusivos, 14 compartidos y 9 parciales. Corrigió restos en
F2.0, el esquema documental de escenarios, el comentario Rust y el resumen de
DTPD-117. Las puertas Cargo completas y `git diff --check` terminaron con código
0. No quedan preguntas biológicas abiertas en esta revisión.

La revisión documental de simplificación queda cerrada. No reintroducir
`verify-f2-step` ni la infraestructura retirada durante la implementación.
