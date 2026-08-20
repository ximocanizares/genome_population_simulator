# Punto de reactivación: PAUSA-F1.2 aprobada

Este documento es el punto autoritativo para la próxima sesión. El historial
extenso de `punto-reanudacion.md` conserva decisiones anteriores, pero no cambia
la siguiente acción.

## Estado confirmado

- Base recuperable: `66ad8b2`, publicada en `origin/main`.
- Plan vigente: `plan-reimplementacion-f1-0.2.md`.
- F1.1 completada y revisada.
- F1.2 reimplementada, corregida y validada con 134 pruebas y dos doctests.
- F1.2.1–F1.2.4 están completados y la tercera revisión independiente aprobó
  PAUSA-F1.2 sin hallazgos nuevos.
- F1.3 y F2 no han comenzado.
- F2 está documentalmente aprobable según la revisión adversarial 0.4, pero solo
  después de cerrar y fijar F1.

## Cierre alcanzado

La primera revisión confirmó cinco defectos:

1. reducción admitida erróneamente sobre referencia simulada;
2. pruebas CLI de no publicación que usan el subcomando o la ubicación TOML
   equivocados y, por tanto, no acreditan el contrato;
3. testigos CONF-04/05 de una sola secuencia que no cubren el requisito completo;
4. esquema documental de escenarios que atribuye capacidades futuras a F1.2;
5. regresión introducida en el contrato mínimo del manifiesto F1.4: esquema y
   validador aceptan estructuras inválidas tras retirar sus pruebas negativas.

Todos quedaron corregidos. La segunda revisión encontró dos residuos en rutas,
esquemas y evidencia histórica; también quedaron corregidos. La tercera revisión
no encontró hallazgos nuevos y aprobó la pausa.

## Primera acción al reanudar

Solicitar o comprobar la aprobación explícita del usuario para comenzar F1.3.
Sin esa aprobación, no ejecutar F1.3.1 ni modificar RNG.

## Restricciones de continuidad

- No recuperar `ResolutionDomain`, `resolution-evidence`, reglas, poblaciones,
  eventos o escenarios en el runtime F1.
- No ampliar el modelo de amenaza a manipulación local activa.
- No descartar el árbol sucio ni restaurar ficheros en bloque.
- No hacer commit/push ni iniciar F2 sin orden del usuario.
- Preservar la sintaxis compacta `[reference.composition]` salvo que una revisión
  demuestre una contradicción normativa.

La copia de seguridad del hilo está en
`docs/chat-reimplementacion-f1-2026-08-15.md`.
