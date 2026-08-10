# Contraste de distribución F1 (E12)

- Fecha: 10 de agosto de 2026.
- Naturaleza: prueba de base de distribución; no es publicación F9 ni declara
  una entrega oficial v1.
- Construcción: `cargo build --release --locked --target
  x86_64-unknown-linux-gnu`, sin `target-cpu=native`.

## `.tar.gz`

Orden: `scripts/package-f1-distribution.sh
/tmp/gps-f1-distribution-e12`. El paquete incluyó ejecutable, cinco esquemas,
licencia, ejemplo y receta OCI. GNU tar 1.35 lo listó y extrajo; el ejecutable
extraído respondió `genome-population-simulator 0.1.0` a `--version`.

SHA-256 independiente del archivo, regenerado tras H-F15-R01–R05:
`61010ac25197a9f3868f027a0e8cf5e169941da338f54a5a5dc7406e6fe5968a`.

## OCI

La receta `distribution/f1/Dockerfile` fija el índice Ubuntu 24.04
`sha256:561618e2c15bf2397621dd04f96926663a3b5616c189cf7e38db7e82f5c538ea`;
el manifiesto linux/amd64 observado fue
`sha256:019e8eb29a85e74d64925745884f2ec79aa27e3feab36353d24656f4d6b89467`.

La construcción Docker 29.6.1 produjo:

- índice: `sha256:c6680cf33973f07121cf27f08c559e8a884f2618e6db38d5e99cf30eb8367c23`;
- manifiesto: `sha256:a1ab31ebdd28c2084d9f83f42b66d73eede33b3ff5ce012b62ff013482f2c7dc`;
- configuración: `sha256:552d57525ca7925ca657a83e6272743712cb725aada8e28ea8de503e028421f2`;
- cinco capas inspeccionadas mediante `docker image inspect`.

El comando aprobado `docker run --rm genome-population-simulator:f1-e12
--version` terminó con código 0 y mostró la versión esperada. La prueba se
ejecutó en el daemon OCI disponible desde Ubuntu 24.04/WSL2; no sustituye una
publicación F9 ni altera decisiones pendientes.
