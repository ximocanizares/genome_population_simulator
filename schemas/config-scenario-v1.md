# Esquema de autoría TOML: envoltura de escenario v1

Todo escenario declara primero:

```toml
config_kind = "scenario"
schema_version = "1"
```

Las claves raíz admitidas en F1.2 son `config_kind`, `schema_version`, `id`,
`route`, `phases`, `events` y `crosses`. `id` es obligatorio y `route` es
exactamente `evolution` o `families`. `phases`, `events` y `crosses` son listas
de tablas con `id` explícito; su orden solo tiene efecto donde el requisito
funcional lo declara normativo. No se admiten campos desconocidos, claves
duplicadas ni IDs repetidos.

La envoltura no ejecuta todavía evolución ni cruces: F1.2 solo conserva su
clase, versión, identidad, vía, ubicación y listas declaradas para validarlas
antes de las fases biológicas consumidoras.
