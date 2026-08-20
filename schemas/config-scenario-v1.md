# Borrador documental para F2: envoltura de escenario v1

Este documento conserva el contrato previsto para F2.
No forma parte del contrato ejecutable de F1. El binario de esa fase solo admite
`config_kind = "general"`, rechaza cualquier enlace de escenario en la raíz y
no abre ficheros de escenario.

La envoltura prevista declara primero:

```toml
config_kind = "scenario"
schema_version = "1"
```

Las claves raíz previstas para F2 son `config_kind`, `schema_version`, `id`,
`route`, `phases`, `events` y `crosses`. `id` es obligatorio y `route` es
exactamente `evolution` o `families`. `phases`, `events` y `crosses` son listas
de tablas con `id` explícito; su orden solo tiene efecto donde el requisito
funcional lo declara normativo. No se admiten campos desconocidos, claves
duplicadas ni IDs repetidos.

Este texto no habilita todavía evolución ni cruces. Su admisión, validación y
consumo se implementarán en F2 con las fases biológicas correspondientes.
