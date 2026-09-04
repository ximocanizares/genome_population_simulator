# Nota de rediseño de la admisión F5.1.2 v7

Fecha: 2026-09-04
Estado: diseño aprobado; implementación y cierre pendientes

## Autoridad y alcance

Esta nota traza la decisión de sustituir las correcciones incrementales
`f5-resource-accounting-v5` y v6 por una sola admisión
`f5-resource-accounting-v7`. La autoridad normativa queda en DTC-089 y en la
arquitectura §15; esta nota no las sustituye ni acredita código ya corregido.

El rediseño conserva la biología, RNG, IDs, geometría, dos pasadas FASTA, una
pasada de máscara y contratos de perfil. Cambia la frontera física de entrada,
la contabilidad de recursos, el binding del preflight y el objeto
`resource_preflight` prerelease de efectivo/manifiesto v5.

## Hallazgos cubiertos

| Problema de revisión | Causa de raíz | Regla v7 | Evidencia de cierre requerida |
|---|---|---|---|
| H13: una sonda CLI de 4096 bytes podía derivar a `fs::read` | clasificación basada en prefijo y lectura completa previa a la cota | TOML general completo bajo 65536+byte de exceso; clasificación solo tras EOF | discriminadores posteriores a 4096, 65536/65537 y búsqueda negativa de la sonda y `fs::read` F5 |
| H13: TOML/perfil podían materializarse antes del presupuesto | la carga precedía al escenario o usaba tamaños externos | bootstrap acotado; después, ledger acumulativo antes de cada apertura y reserva | espías de aperturas, exacto/exacto−1 y primer rechazo sin recursos posteriores |
| H12: la segunda pasada FASTA podía recorrer un fichero cambiado sin cota física | el SHA detectaba el cambio solo tras EOF | límite `tamaño_recibido+1`, tamaño y SHA idénticos antes de aceptar REF | truncado, crecimiento y sustitución entre pasadas |
| Grave: preflight reutilizable entre documentos físicamente distintos | binding solo semántico | identidad por rol+ruta normalizada+tamaño observado+SHA de todas las entradas | cruce de TOML con distinto padding e intercambio de TSV |
| Grave: memoria calculada como suma/pico global impreciso | no existía propiedad ni vida por fase | máximo de picos de fase con retenido/transitorio y crecimiento prospectivo | reconstrucción independiente de cada fase y del máximo |
| Grave: tamaños obtenidos con `metadata` se trataban como límites | observación y lectura ocurrían en handles distintos | contador, hash y EOF en el mismo `BoundedInput`; atributos solo diagnósticos | mutación entre consulta de atributos y lectura |
| Grave: trabajo físico de entrada no era acumulativo | costes se incorporaban tarde o una sola vez | documentos cuentan lectura/hash y parser; FASTA reserva dos visitas; máscara una | oráculo de bytes/operaciones y frontera exacta |
| Grave: el fixture de memoria se autocambiaba al editar el límite | el literal presupuestario altera tamaño, digest y pico | línea de anchura fija e iteración hasta punto fijo; variante −1 conserva longitud | dos iteraciones idénticas y recálculo desde bytes finales |

## Orden obligatorio

1. Leer, medir, hashear y clasificar el general bajo cota.
2. Leer y parsear el escenario bajo cota; crear el ledger y admitir el
   bootstrap.
3. Leer metadata de perfil bajo el restante; reservar los cinco tamaños
   declarados y procesar cada TSV secuencialmente.
4. Resolver dimensiones y salidas conocidas antes del escaneo.
5. Escanear FASTA y máscara con cargo incremental, reservando la segunda visita
   FASTA antes de consumir la primera.
6. Cerrar fórmulas y emitir `SealedF5Admission` antes de RNG o destino.
7. Ejecutar layout y recuperar REF con el límite e identidad del recibo.

## Compatibilidad

- Los receipts y valores `resource_preflight` v5/v6 no se aceptan ni migran; se
  repite la admisión desde las entradas originales.
- General v3 y escenario v1 conservan campos y significado. La ruta CLI
  transparente queda físicamente acotada a 65536 bytes antes de clasificar; un
  documento mayor se rechaza, no se deriva al cargador histórico.
- Efectivo y manifiesto v5 son contratos prerelease y deben actualizarse de
  forma atómica para publicar algoritmo v7, identidades y fases. No puede
  mezclarse un esquema v7 con un productor/lector v5 o v6.
- No se añaden crates. El modelo de amenaza sigue siendo local: se detectan
  cambios accidentales, no se garantiza resistencia a sustituciones hostiles
  concurrentes.

## Gates que faltan

La nota no cierra H12/H13 ni los graves. El cierre requiere las tres pausas
F5.1.2-A/B/C del plan, pruebas independientes, schemas coherentes,
`cargo fmt --check`, Clippy estricto, suite dirigida y completa, validación de
JSON Schema y búsquedas negativas de rutas v5/v6, sonda 4096 y lectura completa
no acotada.
