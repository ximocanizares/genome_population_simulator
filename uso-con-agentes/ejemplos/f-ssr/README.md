# f-ssr: endogamia (F = 0,5) con SNV y un microsatélite

200 individuos diploides con `genotype_model = "f_diploid"` e
`inbreeding_coefficient = 0.5`, cuatro SNV (`chr1:200`, `500`, `1500`, `1800`, con
frecuencia alternativa 0,5) y un microsatélite (CA) con cuatro alelos explícitos
de 22 a 25 unidades, equiprobables. Solo verdad y VCF, sin lecturas.

```bash
gps run general.toml --output ../../../../mi-salida/f-ssr
```

Resultado con la semilla incluida (F pedida 0,5, con ruido de muestreo de 200
individuos): los SNV dan F_IS de 0,64, 0,56, 0,50 y 0,47, y el microsatélite
0,49. Hasta el 25/10/2026 el microsatélite ignoraba F y daba -0,02; si al
ejecutarlo sale cerca de 0, el binario es anterior al arreglo y hay que
recompilarlo.

F se aplica locus a locus e independientemente (ver `GUIA.md`).
