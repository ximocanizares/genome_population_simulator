# Manifiesto de evento v1

`events/<event_id>/manifest.json` usa `manifest_kind = "gps_event_manifest"`
y `schema_version = "1"`. Describe exclusivamente el primer evento técnico
F3: su tecnología explícita `illumina-short-read-v1`, perfil Illumina single-end
mínimo, tanda, bibliotecas, referencia común y FASTQ. Las rutas de los ficheros
enlazados son relativas a la raíz de la ejecución; el ensamblado y el enlace
desde `manifest.json` general pertenecen a F3.5.3.

La v1 no modela conjuntos de mapeo, mapeadores, sesiones, paquetes exportables
ni procedencia durable por lectura. Un artefacto `complete` conserva su ruta,
tamaño y SHA-256; uno `incomplete` o `failed` no inventa una identidad física y
enlaza un diagnóstico ya declarado por el manifiesto. Existe exactamente un
resultado FASTQ por biblioteca seleccionada; `artifacts` se ordena por
`individual_id`, `library_id`, `id` y `diagnostics` por `code`. Un evento
`complete` solo contiene FASTQ completos y ningún diagnóstico; uno `incomplete`
o `failed` contiene al menos un FASTQ no completo y un diagnóstico.
