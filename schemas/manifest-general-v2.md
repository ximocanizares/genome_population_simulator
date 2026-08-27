# Manifiesto general v2

F2 publica `manifest.json` con `manifest_kind = "gps_execution_manifest"` y
`schema_version = "2"`. Su forma ejecutable está en
`manifest-general-v2.schema.json`. F3.5.3 conserva esta versión para lectores
existentes y añade `manifest-general-v3.schema.json` para los enlaces de evento.

La v2 incorpora la verdad obligatoria, los enlaces de matriz SNV, las
autorizaciones y las salidas VCF. El manifiesto F1 permanece congelado en
`manifest-general-v1.schema.json` y se valida por su lector v1 dedicado. No
hay migración automática: la versión declarada selecciona el contrato que se
lee.
