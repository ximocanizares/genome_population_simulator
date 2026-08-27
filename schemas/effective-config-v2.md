# Configuración efectiva v2 (histórica)

F2 publicó `effective-config.json` con `effective_config_kind =
"gps_effective_config"` y `schema_version = "2"`. Su forma ejecutable está en
`effective-config-v2.schema.json`.

La v2 añadió campos obligatorios de F2 —ploidía, normalización, autorizaciones,
poblaciones iniciales y, cuando corresponda, eventos directos, matriz SNV y
VCF— que no pueden presentarse como una extensión inocua de la v1. Por ello la
v1 permanece congelada en `effective-config-v1.schema.json`. F3 conserva este
esquema como registro, pero no mantiene un lector v2 activo: la configuración
efectiva vigente es v3 para impedir que `technical_event` se interprete como
una extensión silenciosa de v2. El lector actual conserva la ruta archivada v1
y no migra ni reescribe ese documento.
