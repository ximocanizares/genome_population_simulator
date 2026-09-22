# Evaluación de validadores independientes de F1

- **Versión del registro:** 2.
- **Fecha de evaluación:** 8 de agosto de 2026.
- **Revalidación:** 22 de septiembre de 2026 (CPython; ver más abajo).
- **Ámbito:** pruebas, CI y aceptación; nunca el recorrido normal del producto.
- **Plataforma observada:** Ubuntu 24.04 sobre WSL2, x86_64, repositorio en el
  sistema de ficheros Linux.

## Herramientas aprobadas

| Clase | Herramienta exacta | Licencia y mantenimiento | Comando literal aprobado | Uso y límites |
| --- | --- | --- | --- | --- |
| FASTA/FAI | samtools 1.21 con htslib 1.21 | MIT/BSD; proyecto HTSlib mantenido | `samtools faidx "$fasta"` | Valida estructura indexable y consultas FAI en pruebas. No sustituye la validación interna del alfabeto ni se instala para ejecutar el producto. |
| Integridad | GNU `sha256sum` 9.4 | GPL-3.0-or-later; GNU coreutils mantenido por la distribución Ubuntu | `sha256sum --check "$checksums"` | Oráculo externo de SHA-256 sobre bytes completos. |
| `.tar.gz` | GNU `tar` 1.35 | GPL-3.0-or-later; GNU tar mantenido por la distribución Ubuntu | `tar --list --gzip --file "$archive"` y `tar --extract --gzip --file "$archive" --directory "$destination"` | Comprueba listado y extracción; no define todavía el contenido de publicación de E12. |
| JSON Schema | CPython 3.12.3 **o** 3.12.7, con `jsonschema` 4.10.3 | MIT; paquete mantenido por Ubuntu | `python3 -m jsonschema -i "$instance" "$schema"` | Implementación Draft 2020-12 ajena al futuro productor Rust. Es dependencia de prueba, no de Cargo ni del runtime. |
| OCI | Docker CLI 29.6.1, build 8900f1d | Apache-2.0; CLI mantenida por Docker | `docker image inspect "$image"` y `docker run --rm "$image" --version` | Inspección y arranque de aceptación. Requiere un daemon disponible en CI/publicación; no entra en el binario. |

Las cinco opciones son herramientas externas y no filtran tipos ni errores al
dominio. No añaden código nativo, scripts de construcción, características ni
licencias al grafo Cargo, y su coste se limita al entorno de pruebas. Las
versiones se registran en cada evidencia; un cambio obliga a revalidar los casos
válidos y corruptos.

La evaluación de seguridad limita las herramientas a entradas sintéticas o
artefactos producidos por el proyecto, ejecuta los procesos sin privilegios y
no concede red al validador. Docker solo consulta el daemon de aceptación y las
imágenes fijadas; no ejecuta una imagen no identificada por digest. Las cinco
interfaces CLI están estabilizadas para los comandos seleccionados. Su
rendimiento no afecta al runtime ni a las personas usuarias y se medirá solo
como coste de CI si llegara a dominarla. Las instalaciones usadas proceden de
Ubuntu o de la distribución oficial ya disponible en el entorno de aceptación;
no se enlazan al producto. No se detecta una necesidad de código nativo dentro
del binario, una característica opcional o un script de construcción Cargo.

La receta OCI no existe en F1.1 y no se anticipa: E12 la creará junto con la
imagen. En ese mismo cambio la primera línea fijará la base Ubuntu 24.04 mediante
`FROM ...@sha256:<digest>` y la evidencia conservará el digest de la receta, del
descriptor de índice o manifiesto, de la configuración y de cada capa. No se
aceptará una etiqueta mutable ni un digest pendiente al cerrar E12.

## Revalidación de CPython, 22 de septiembre de 2026

El registro aprobaba una sola versión de CPython, la 3.12.3 que trae Ubuntu
24.04. La puerta empezó a fallar con «unexpected CPython version: 3.12.7», y la
causa **no era una actualización del sistema**: el intérprete del sistema sigue
siendo el 3.12.3. Lo que ocurre es que `python3` resuelve al del **entorno
virtual del proyecto**, que es 3.12.7, cuando ese entorno está activado, y al
del sistema cuando no lo está. La puerta pasaba o fallaba según desde dónde se
lanzara.

Conforme a la regla de este mismo registro —«un cambio obliga a revalidar los
casos válidos y corruptos»—, el caso corrupto de esquema se ejecutó con las dos
versiones sobre la misma entrada:

| | CPython 3.12.3 | CPython 3.12.7 |
|---|---|---|
| Código de salida | 1 | 1 |
| Salida estándar | vacía | vacía |
| Diagnóstico | `Additional properties are not allowed ('unknown_field' was unexpected)` | idéntico |

Y la puerta completa se ejecutó por los dos caminos: **5 PASS** en ambos, con
los `stderr` de las cuatro clases no-OCI idénticos byte a byte.

Se aprueban por tanto **las dos versiones**, como lista explícita de dos y no
como patrón de familia: una 3.12 futura no queda aprobada de antemano. Los
metadatos de cada evidencia siguen registrando cuál corrió, de modo que la
evidencia no pierde precisión.

**Aviso para quien la ejecute.** El script crea el directorio de evidencia
*antes* de comprobar las versiones, así que una ejecución que falle por versión
deja el directorio creado y la siguiente se niega a continuar porque ya existe.
Hay que usar un nombre nuevo.

## Formato obligatorio de evidencia

Cada ejecución de un validador conserva, como mínimo:

| Campo | Contenido |
| --- | --- |
| `evidence_version` | `1` |
| `validator_class` | `fasta`, `sha256`, `tar_gz`, `json_schema` u `oci` |
| `tool` y `tool_version` | Nombre y salida de versión exactos |
| `platform` y `filesystem` | Sistema, arquitectura, WSL2/nativo/OCI y FS |
| `command` | Orden literal, sin ocultar parámetros efectivos |
| `fixture` | Ruta e identidad del válido o corrupto usado |
| `expected` y `observed` | Aceptación/rechazo esperado, código de salida y resumen |
| `stdout_sha256` y `stderr_sha256` | Identidad de las salidas conservadas |
| `artifact_sha256` | SHA-256 del artefacto sometido a prueba |

El script `scripts/check-independent-validators.sh` constituye la prueba
negativa reproducible de A7 y exige como único argumento un directorio de
evidencia nuevo. Por cada clase conserva allí copias estables de las entradas,
los ficheros completos `<clase>.stdout` y `<clase>.stderr` y un
`<clase>.metadata` con todos los campos de la tabla anterior. El script no
sobrescribe un directorio previo; por ejemplo:

```bash
./scripts/check-independent-validators.sh evidence/f1.1/validators-corrupt-20260809
```

Las herramientas se ejecutan sobre nombres relativos dentro de la copia de
evidencia, de modo que stdout y stderr no contienen rutas temporales. Un fallo
de infraestructura —por ejemplo, un daemon OCI inaccesible— no cuenta como
rechazo correcto del artefacto, no genera metadatos de éxito para esa clase e
invalida la ejecución completa. Para cada clase también se comprueban el código
de salida y el diagnóstico semántico esperado; un código distinto o un error
ajeno a la corrupción hace fallar la prueba.

## Resultado de la evaluación

Las herramientas se aprueban para el uso limitado anterior. La prueba negativa
de A7 exige que cada clase rechace su caso corrupto: FASTA sin cabecera,
checksum discordante, falso `.tar.gz`, documento con un campo no permitido e
identidad OCI inexistente con digest explícito. Las pruebas positivas completas
se incorporarán con los productores de B–E; A7 no anticipa esos productores.
