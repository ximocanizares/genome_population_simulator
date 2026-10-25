# ssr-sin-pcr: el mismo escenario sin PCR

Idéntico a `ssr-pcr` (misma semilla, misma población, mismo microsatélite) pero con
`library_preparation = "pcr_free"` y sin la tabla `[technical_event.pcr]`. La
verdad y el VCF son los mismos; cambian las lecturas.

```bash
gps run general.toml --output ../../../../mi-salida/ssr-sin-pcr
```

Valores comprobados: en `I001` cruzan el tracto 17 de 400 lecturas (11 con 21
unidades y 6 con 24), frente a 28 con PCR. No hay duplicados; el tartamudeo
(`artifacts.stutter` del manifiesto del evento) afectó a 10 de los tractos
evaluados.
