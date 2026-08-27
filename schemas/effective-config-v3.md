# Configuración efectiva v3

F3 publica `effective-config.json` con `effective_config_kind =
"gps_effective_config"` y `schema_version = "3"`. Su forma ejecutable está en
`effective-config-v3.schema.json`.

La v3 añade el campo opcional `technical_event` de F3 a los campos efectivos de
F2 —ploidía, normalización, autorizaciones, poblaciones iniciales y, cuando
corresponda, eventos directos, matriz SNV y VCF—. No redefine la semántica de
v2: ese esquema permanece congelado como registro histórico; el lector vigente
acepta v1 y v3, nunca migra ni reescribe un documento v1.

Cuando el TOML declara `[technical_event]`, la configuración efectiva conserva
su ID, individuos y objetivo bajo `technical_event`, y fija explícitamente el
perfil `illumina-se-minimal-v1`, `read_length = 150` y las calidades efectivas
`acgt_phred = 40` y `n_phred = 0`. El objetivo se expresa como exactamente uno
de `target.coverage_wgs` o `target.read_count`; no materializa perfiles F4,
pares, contaminación, referencias de mapeo ni alineamientos.
