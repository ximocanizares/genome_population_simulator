# Evidencia A1–A3: toolchain, lock y puertas Cargo

- Fecha: 8 de agosto de 2026.
- Plataforma: Ubuntu 24.04 sobre WSL2, x86_64, Linux
  `6.18.33.2-microsoft-standard-WSL2`.
- Sistema de ficheros: ext4 bajo `/home/ximo/devel`, fuera de `/mnt/c`.
- `rustc`: `rustc 1.97.1 (8bab26f4f 2026-07-14)`.
- `cargo`: `cargo 1.97.1 (c980f4866 2026-06-30)`.
- Target: `x86_64-unknown-linux-gnu`.
- SHA-256 de `Cargo.lock`:
  `f872ca0dfa2527a527b1de24a049a89da4518c7aee07711c3401da32e36f67bb`.
- Árbol de características: solo
  `genome-population-simulator v0.1.0`; cero dependencias externas.

## Salida de las cinco puertas

### Formato

```text
$ cargo fmt --all -- --check
[sin salida; código 0]
```

### Clippy

```text
$ cargo clippy --workspace --all-targets --locked -- -D warnings
Finished `dev` profile [unoptimized + debuginfo] target(s)
```

### Pruebas de todos los targets

```text
$ cargo test --workspace --all-targets --locked
unit tests src/lib.rs: 0 passed; 0 failed
unit tests src/main.rs: 0 passed; 0 failed
tests/normative_catalogs.rs: 6 passed; 0 failed
```

Las seis pruebas comprueban el catálogo diagnóstico, las representaciones
humana/JSON, la asignación de rechazos B–E, los cuatro vectores de subsemilla y
los cuatro vectores tipados de ID.

### Pruebas de documentación

```text
$ cargo test --workspace --doc --locked
Doc-tests genome_population_simulator
test result: ok. 0 passed; 0 failed
```

### Construcción portable

```text
$ cargo build --release --locked --target x86_64-unknown-linux-gnu
Finished `release` profile [optimized] target(s)
```

## Resolución bloqueada

```text
$ cargo metadata --locked --format-version 1
workspace_members: genome-population-simulator 0.1.0
dependencies: []
features: {}

$ cargo tree --locked -e features
genome-population-simulator v0.1.0
```

Todos los comandos terminaron con código 0. No se configuró
`target-cpu=native`.
