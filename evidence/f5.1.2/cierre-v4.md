# Cierre F5.1.2/v4

## Alcance y base

Esta evidencia se refiere únicamente al cierre de F5.1.2 sobre la base F5
`c73a14e` y el árbol compartido que contiene sus correcciones posteriores. No
atribuye a F5 el commit `fe5e3ab`: ese commit cambia el plan VE y es ajeno a
este cierre. Tampoco modifica ni certifica `TO_DO.md` ni
`design/plan-validacion-externa-fastq-vcf-0.1.md`; sus cambios presentes en el
árbol quedan fuera del diff y de los resultados siguientes.

| Hallazgo | Gravedad registrada | Decisión | Estado | Evidencia verificable |
| --- | --- | --- | --- | --- |
| H1 | Antecedente histórico | No reconstruir un informe no persistido | Antecedente resumido | El informe detallado no existe en el repositorio. El cierre actual conserva el código y gates F5, pero no reclasifica ni inventa su resultado. |
| H2 | Antecedente histórico | No reconstruir un informe no persistido | Antecedente resumido | Igual que H1. |
| H3 | Antecedente histórico | No reconstruir un informe no persistido | Antecedente resumido | Igual que H1. |
| H4 | Grave | Aplicar | Aplicado | `F5ResourcePreflight` liga identidad semántica de admisión (referencia, perfil, máscara, objetivo, ploidía, fundadores, poblaciones, overrides y límites) y `F5LayoutRequest` exige las mismas aportaciones `other_*`. `f5_resource_preflight` cubre máscaras de igual capacidad con tramos distintos, límites y recursos externos. |
| H5 | Grave | Aplicar con adaptación | Aplicado con adaptación | El layout recupera únicamente bases REF seleccionadas por streaming FASTA y verifica SHA-256, sin construir `Reference`; la prueba de FASTA de 2 MiB conserva la cota de memoria admitida. |
| H6 | Antecedente histórico | No reconstruir un informe no persistido | Antecedente resumido | El informe detallado no existe en el repositorio; no se afirma una causa, severidad ni prueba que no estén conservadas. |
| H7 | Antecedente histórico | No reconstruir un informe no persistido | Antecedente resumido | Igual que H6. |
| H8 | Antecedente histórico | No reconstruir un informe no persistido | Antecedente resumido | Igual que H6. |
| H9 | Antecedente histórico | No reconstruir un informe no persistido | Antecedente resumido | Igual que H6. |
| H10 | Grave | Aplicar | Aplicado | `f5_layout::tests::production_v4_instrumentation_covers_fragmentation_interior_last_residues_and_retry` usa carga/configuración, preflight real, binding, streaming FASTA y layout de producción; no usa `from_test_parts`. |
| H11 | Grave documental | Aplicar con adaptación | Aplicado con adaptación | Plan F5, arquitectura, DTC vigente y esta evidencia emplean `f5-resource-accounting-v4`; se actualizó la fórmula/delta y se eliminó la referencia residual v3 del recorrido F5. |

## H10 — recorrido instrumental v4

El fixture fija `T=3` tramos, `B=4` bloques, `L=2` intentos y semilla `01`×32.
El primer intento se descarta y el segundo completa el layout. La traza exige,
en el intento fallido, un ganador interior con residuo izquierdo y derecho y
una selección del último tramo inicial. El mismo fixture primero admite
`max_operation_units=operation_units` y después rechaza
`max_operation_units=operation_units-1` con `GPS-RES-101`, antes de layout.

| Contador/pico instrumentado | Comparación que ejecuta la prueba |
| --- | --- |
| Escaneo de capacidad | Igualdad con `T`: una visita por tramo en el preflight real. |
| Escaneo del request | Igualdad con `T`: una visita por tramo antes de RNG. |
| Copias iniciales | Igualdad con `L*T`. |
| Reservas | Dos por intento y `L*(T+2*B)` slots solicitados. |
| Escaneos total y prefijo | Cada uno `<= L*(B*T+B*(B-1)/2)`. |
| Crecimiento y división | Cada uno `<= L*B`; una colocación solo puede añadir un tramo neto y dividir una vez. |
| Desplazamientos `remove` | `<= L*(B*T+B*(B-1)/2)`. |
| Desplazamientos `insert` | `<= 2*L*(B*T+B*(B-1)/2)`. |
| Copias+`remove`+`insert` | `<= L*(T+3*B*T+3*B*(B-1)/2)`. |
| Reintentos | Igualdad con `L-1`; la transacción fallida no se publica. |
| Pico de slots/memoria | `slots <= 2*T+4*B` y `bytes=256*slots <= 256*(2*T+4*B)`. |

Estas son las cotas que alimentan v4: las dos visitas iniciales de tramos, las
dos pasadas de selección, la copia por intento y hasta tres desplazamientos por
colocación. Para `tract_count` que pasa de 1 a 4096, el delta confirmado por
la fórmula vigente es `4095*(3+L+5*L*B)` operaciones y `512*4095` bytes de
trabajo, sin cambio de bytes de salida.

## Gates ejecutados

- `cargo fmt --check` — correcto.
- `cargo clippy --all-targets --all-features -- -D warnings` — correcto.
- `cargo test --test f5_resource_preflight` — 25 correctas.
- `cargo test --test f5_0_1_external_schemas` — 6 correctas.
- `cargo test --lib production_v4_instrumentation_covers_fragmentation_interior_last_residues_and_retry` — 1 correcta.
- `cargo test --all-targets --all-features` — correcto (código de salida 0; las
  pruebas marcadas `ignored` permanecen ignoradas).
- `cargo doc --no-deps` — correcto.
- `git diff --check` — correcto.
