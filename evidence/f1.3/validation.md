# Validación de F1.3

Entorno observado: Rust/Cargo 1.97.1, objetivo
`x86_64-unknown-linux-gnu`, Linux WSL2 x86_64; `cargo-audit` 0.22.2 y
`jsonschema` 4.10.3 sobre CPython 3.12.3.

Puertas ejecutadas después de la revisión y las correcciones:

```text
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --workspace --all-targets --locked
cargo test --workspace --doc --locked
cargo doc --workspace --no-deps --locked
cargo build --release --locked --target x86_64-unknown-linux-gnu
cargo metadata --locked --no-deps --format-version 1
cargo tree --locked -e features
cargo audit
python3 -m jsonschema -i evidence/f1.3/diagnostic-invalid-seed.json schemas/diagnostics-v1.schema.json
git diff --check
```

Todas terminaron con código 0. Pasaron 83 pruebas de integración: 20 nuevas
para F1.3 (4 de raíz/CSPRNG, 4 de materiales/subsemillas, 10 casos auditados
`TEST-REQ-*`, 1 regresión específica de condiciones RNG y 1 comprobación
automática de la matriz),
sin fallos ni casos ignorados; los objetivos unitarios y documentales no
contienen todavía casos propios. La documentación y la construcción release
portable terminaron sin avisos.

El árbol bloqueado fija `getrandom 0.4.3`, `blake3 1.8.6` con `pure` y
`rand_chacha 0.10.0`, sin características predeterminadas en las tres
dependencias directas. Se comprobó el `build.rs` selector de BLAKE3 y getrandom,
la ausencia de compilación C/ensamblador con `pure`, las licencias y el grafo de
rasgos completo. No aparece `StdRng`, un estado global, reloj, ruta, hilo,
planificador ni orden accidental en la identidad de flujo.

`cargo audit` cargó 1198 avisos RustSec y examinó 34 dependencias de
`Cargo.lock`; terminó con código 0 y no informó vulnerabilidades. Los cuatro
materiales y subsemillas A5 coinciden byte a byte con el catálogo. El oráculo
BLAKE3 independiente de `tests/normative_catalogs.rs` y la implementación de
producción obtienen los mismos digests.

El CLI aceptó tanto la raíz declarada como la generación desde el CSPRNG. La
raíz mayúscula fue rechazada en modo humano y JSON con `GPS-CONF-017`, salida 2
y sin crear el destino solicitado. El documento conservado en este directorio
fue aceptado de forma independiente por JSON Schema Draft 2020-12.

El fallo de entropía se ejercita mediante el proveedor inyectable en la frontera
de aplicación y produce `GPS-RES-001`, categoría `RES` y salida 3; no existe ni
se afirma una inyección de fallo CSPRNG a través del binario CLI. F1.3 registra
la raíz y su origen en `EffectiveConfig`; su serialización/publicación en
manifiesto permanece diferida a D7/E8/F1.5.

La matriz `test-req-matrix.tsv` audita exactamente los diez requisitos que se
habían asignado a F1.3: solo REQ-REP-23 queda `pass`; 3 se trasladan a F1.5 y 6
a F4.6/F9.7. Las nueve pruebas restantes se conservan y enlazan expresamente
como evidencia parcial, sin sustituir simulaciones, manifiestos, lecturas,
salidas o artefactos todavía inexistentes.

La regresión de H-F13-01 demuestra que reordenar dependencias RNG no cambia el
estado, mientras que cambiar por separado `rand_chacha`, BLAKE3, la versión de
derivación o el catálogo retira la garantía biológica preparatoria. Las
dependencias técnicas permanecen separadas y solo endurecen el nivel técnico y
el binario. La comparación binaria preparatoria requiere clase de artefacto,
política e identificador de evidencia concretos; no se presenta como validación
de REQ-CONF-115.
