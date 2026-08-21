# Punto de reanudación

## Estado actual — 21 de agosto de 2026

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

1. Tras aprobación y commit/push, comenzar el nuevo F2.0.1: exponer a Marimo la
   configuración F1 existente sin introducir conceptos biológicos F2.

La revisión directa de simplificación comprobó 41 pasos bien formados, seis
pausas, 78 REQ exclusivos, 14 compartidos y 9 parciales. Corrigió restos en
F2.0, el esquema documental de escenarios, el comentario Rust y el resumen de
DTPD-117. Las puertas Cargo completas y `git diff --check` terminaron con código
0. No quedan preguntas biológicas abiertas en esta revisión.

La revisión documental de simplificación queda cerrada. No reintroducir
`verify-f2-step` ni la infraestructura retirada durante la implementación.
