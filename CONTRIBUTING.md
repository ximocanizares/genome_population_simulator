# Desarrollo

El entorno fijado es Rust 1.97.1 con edición 2024. El desarrollo ordinario se
realiza en Ubuntu 24.04 nativo o en Ubuntu 24.04 sobre WSL2, siempre sobre el
sistema de ficheros Linux (ext4) y nunca bajo `/mnt/c`. La entrega usa el target
`x86_64-unknown-linux-gnu` y CPU portable; `target-cpu=native` está prohibido.

## Puertas Cargo rápidas

Ejecutar antes de entregar un cambio:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --workspace --all-targets --locked
cargo test --workspace --doc --locked
cargo build --release --locked --target x86_64-unknown-linux-gnu
```

Las pruebas unitarias y de integración de Rust forman las puertas rápidas. Los
validadores independientes de formatos, integridad, paquetes e imágenes OCI se
ejecutan únicamente en pruebas, CI o aceptación de una publicación. Las
observaciones manuales de recursos se registran aparte y no convierten una
variación temporal pequeña en un fallo funcional.

La disponibilidad y el rechazo deliberado de entradas corruptas se comprueban
por separado con:

```bash
./scripts/check-independent-validators.sh
```

Sus versiones, comandos literales y formato de evidencia están fijados en
`design/evaluacion-validadores-independientes-f1.md`.

Toda dependencia nueva requiere primero una evaluación versionada conforme a
DTC-043. El archivo `Cargo.lock` de la aplicación se versiona y todos los usos
de CI y entrega emplean `--locked`.
