# Configuración efectiva v2

F2 publica `effective-config.json` con `effective_config_kind =
"gps_effective_config"` y `schema_version = "2"`. Su forma ejecutable está en
`effective-config-v2.schema.json`.

La v2 añade campos obligatorios de F2 —ploidía, normalización, autorizaciones,
poblaciones iniciales y, cuando corresponda, eventos directos, matriz SNV y
VCF— que no pueden presentarse como una extensión inocua de la v1. Por ello la
v1 permanece congelada en `effective-config-v1.schema.json`; los lectores del
programa aceptan explícitamente ambas versiones y nunca migran ni reescriben
un documento v1.
