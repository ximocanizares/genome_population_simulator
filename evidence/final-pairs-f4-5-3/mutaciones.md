# Mutaciones del productor frente a la batería de REQ-VAL-03 (27/10/2026)

Evidencia del paso F4.5.3 del plan de F4 y de DTC-101: la comparación pareja a
pareja del FASTQ final con la verdad la garantiza la batería obligatoria
`tests/final_pairs_against_truth.rs`, juzgada por el oráculo independiente de
`tests/final_pairs_oracle/`. Esta tabla demuestra que la batería **falla**
cuando el productor se equivoca en cada etapa.

## Cómo se hizo

Cada mutación se aplicó al código de producción en el árbol real, con copia de
seguridad fuera de git, y se ejecutó
`cargo test --locked --test final_pairs_against_truth` (doce escenarios).
Después se restauró la copia y se comprobó que la huella SHA-256 del fichero
volvía a la de antes. Las huellas van abreviadas a 12 caracteres: antes → mutado.
Las treinta se restauraron con la huella original.

«Aborta el propio programa» quiere decir que una guarda del simulador detiene
la ejecución antes que el oráculo. Para esas cinco hay una variante con la
guarda desactivada o esquivada (sufijo `s`) que demuestra que el oráculo también las caza
por su cuenta.

## Tabla

| Código | Mutación | Escenarios que fallan | Quién la caza | SHA-256 |
|---|---|---|---|---|
| M01 | inicio del fragmento desplazado una base | PCR sola, clonal con quimeras, contaminación con salto, diploide, dos secuencias con quimeras, haploide, quimeras y microsatélites en dos secuencias, reproducción de DT-060, saturado, tartamudeo clonal, tartamudeo de base, tetraploide | oráculo | `src/wgs_sampling.rs` 554e125aef72 → 1cc25c7a31b9 |
| M02 | copia cromosómica equivocada (la siguiente candidata) | clonal con quimeras, contaminación con salto, diploide, dos secuencias con quimeras, quimeras y microsatélites en dos secuencias, saturado, tartamudeo clonal, tartamudeo de base, tetraploide | oráculo | `src/wgs_sampling.rs` 554e125aef72 → d70445f47966 |
| M03 | copias recorridas en orden de secuencia inverso | dos secuencias con quimeras, quimeras y microsatélites en dos secuencias | oráculo | `src/wgs_sampling.rs` 554e125aef72 → 2c073f07fd87 |
| M04 | media del inserto desplazada una base | PCR sola, clonal con quimeras, contaminación con salto, diploide, dos secuencias con quimeras, haploide, quimeras y microsatélites en dos secuencias, reproducción de DT-060, saturado, tartamudeo clonal, tartamudeo de base, tetraploide | oráculo | `src/wgs_sampling.rs` 554e125aef72 → 5b1b70c93634 |
| M05 | R2 invertido sin complementar | PCR sola, clonal con quimeras, contaminación con salto, diploide, dos secuencias con quimeras, haploide, quimeras y microsatélites en dos secuencias, reproducción de DT-060, saturado, tartamudeo clonal, tartamudeo de base, tetraploide | oráculo | `src/paired_end.rs` 718e1f9a890f → fea3f3676238 |
| M06 | punto de corte de la quimera desplazado | clonal con quimeras, dos secuencias con quimeras, quimeras y microsatélites en dos secuencias, reproducción de DT-060, saturado | oráculo | `src/chimera.rs` fa88db1cfb72 → b22fa7355503 |
| M07 | candidatas del segundo origen en orden de referencia | dos secuencias con quimeras, quimeras y microsatélites en dos secuencias | oráculo | `src/chimera.rs` fa88db1cfb72 → 96fd981cfb2e |
| M08 | error sorteado pero no aplicado | contaminación con salto, diploide, saturado | oráculo | `src/error_quality.rs` cae0bbf7cfc6 → 87cd40220c0f |
| M09 | error decidido con la clave del ciclo siguiente | contaminación con salto, diploide, saturado | oráculo | `src/error_quality.rs` cae0bbf7cfc6 → 1a1558795314 |
| M10 | calidad desplazada un punto Phred | PCR sola, clonal con quimeras, contaminación con salto, diploide, dos secuencias con quimeras, haploide, quimeras y microsatélites en dos secuencias, reproducción de DT-060, saturado, tartamudeo clonal, tartamudeo de base, tetraploide | oráculo | `src/error_quality.rs` cae0bbf7cfc6 → aad21b9ec3a3 |
| M11 | calidad de N distinta de cero | diploide, haploide | aborta el propio programa | `src/error_quality.rs` cae0bbf7cfc6 → d54c7f27ddd7 |
| M12 | contaminación que no sustituye (repite la propia retirada) | contaminación con salto, saturado | oráculo | `src/biological_contamination.rs` b7674399e2ea → 695f5a9dba3d |
| M13 | contaminación con otra donante (orden de candidatas invertido) | contaminación con salto, saturado | oráculo | `src/biological_contamination.rs` b7674399e2ea → 1493553219ea |
| M14 | donantes de cualquier población | contaminación con salto | aborta el propio programa | `src/biological_contamination.rs` b7674399e2ea → e6f134855c4f |
| M15 | salto de índice que no mueve la pareja | contaminación con salto, saturado | aborta el propio programa | `src/demultiplexing.rs` e39000ac44ce → a20ba026da6b |
| M16 | salto de índice que duplica la pareja | contaminación con salto, saturado | aborta el propio programa | `src/demultiplexing.rs` e39000ac44ce → f660cddee344; `src/execution.rs` 1d2f807f448e → 1a5f04a79002 |
| M17 | PCR sin duplicados (familias de tamaño 1) | PCR sola, clonal con quimeras, saturado, tartamudeo clonal | oráculo | `src/wgs_sampling.rs` 554e125aef72 → 5f3a319cbfe1 |
| M18 | PCR sin aceptación por composición | PCR sola, saturado | oráculo | `src/wgs_sampling.rs` 554e125aef72 → 3f6dadd9ca61 |
| M19 | linaje: herencia un ciclo más corta | saturado, tartamudeo clonal | oráculo | `src/clonal_stutter.rs` 7bf0c0092e49 → 0f2888fe7a80 |
| M20 | tartamudeo: paso elegido simétrico al sorteado | clonal con quimeras, quimeras y microsatélites en dos secuencias, reproducción de DT-060, saturado, tartamudeo clonal, tartamudeo de base | oráculo | `src/stutter.rs` 3e610853ee36 → fd1b344d33f5 |
| M21 | DT-060 reintroducido (tractos con la copia del prefijo) | clonal con quimeras, quimeras y microsatélites en dos secuencias, reproducción de DT-060, saturado | oráculo | `src/stutter.rs` 3e610853ee36 → 5a0bd680698f |
| M22 | linaje aplicado también al segundo tramo | clonal con quimeras | oráculo | `src/stutter.rs` 3e610853ee36 → 42da1f20f10b |
| M11s | calidad de N distinta de cero, sin la guarda del programa | diploide, haploide | oráculo | `src/error_quality.rs` cae0bbf7cfc6 → d54c7f27ddd7; `src/fastq_validation.rs` 7f5a054c67fe → 11e80d217186 |
| M14s | donantes de cualquier población, sin la guarda del programa | contaminación con salto | oráculo | `src/biological_contamination.rs` b7674399e2ea → e6f134855c4f; `src/results.rs` 96df83dedc7e → f399fe26402e |
| M15s | salto de índice sorteado que no se aplica | contaminación con salto, saturado | oráculo | `src/demultiplexing.rs` e39000ac44ce → 35d05f64aaf6 |
| M16s | salto que duplica la pareja, sin la guarda del programa | contaminación con salto, saturado | oráculo (pareja repetida) | `src/demultiplexing.rs` e39000ac44ce → f660cddee344; `src/execution.rs` 1d2f807f448e → 1a5f04a79002; `src/fastq_validation.rs` 7f5a054c67fe → e1dab222b19d |
| M23 | tractos del sufijo buscados con la secuencia del prefijo (revisión) | quimeras y microsatélites en dos secuencias | oráculo | `src/stutter.rs` 3e610853ee36 → 2170c4b64c42 |
| M24 | calidad de N cambiada solo con errores activos | diploide | aborta el propio programa | `src/error_quality.rs` cae0bbf7cfc6 → a968156e222f |
| M24s | ídem, sin la guarda del programa | diploide | oráculo | `src/error_quality.rs` cae0bbf7cfc6 → a968156e222f; `src/fastq_validation.rs` 7f5a054c67fe → 11e80d217186 |
| M25 | desempate de la tirada más larga por la derecha | tartamudeo de base | oráculo | `src/tandem_repeat.rs` dbb2bb975e13 → a55b3729830a |

M23 es la mutación de la revisión adversarial del 27/10/2026: buscar los
tractos del sufijo con la secuencia del prefijo. Con la batería de entonces,
once escenarios, sobrevivía; la caza el escenario de quimeras y microsatélites
en dos secuencias, añadido por esa revisión. M24 y M25 también salen de esa
revisión: sobrevivían antes de añadir una tanda de `N` con errores activos y un
tracto impuro con tiradas empatadas.

## Las dos mitades de la reverificación del linaje

Medidas sobre `tests/clonal_lineage.rs::a_net_step_that_does_not_fit_the_chimeric_molecule_is_not_written`,
con las cifras congeladas tras DT-060 (243 contracciones y 248 expansiones de
300):

| Mutación en `src/stutter.rs` | Resultado |
|---|---|
| Sin la comprobación de que el paso cabe | la ejecución aborta; la prueba falla |
| Sin la comprobación del tracto ausente | 249 expansiones frente a 248; la prueba falla |

## Bytes que no cambia el arreglo de DT-060

Con el binario anterior al arreglo y el corregido, sobre los mismos ficheros y
la misma ruta: idénticas las siete entradas congeladas (una de
`evidence/stutter-f4-8-congelacion/` y seis de
`evidence/dt-009-multiplicador-congelacion/`), una ejecución con tartamudeo y
sin quimeras y una con quimeras y sin microsatélites. Cambia solo la que
combina quimeras, microsatélites y tartamudeo.
