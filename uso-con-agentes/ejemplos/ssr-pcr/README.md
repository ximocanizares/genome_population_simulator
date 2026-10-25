# ssr-pcr: microsatélite con PCR

20 individuos diploides (equilibrio de Hardy-Weinberg), una referencia de 2000 pb
(`chr1`), un microsatélite (CA)24 en `chr1:1001-1048` cuyos alelos se generan por
deslizamiento, dos SNV (`chr1:300` y `chr1:1600`), lecturas a 30X con PCR (10
ciclos, eficiencia 0,9, 70 % de duplicados esperados) y tartamudeo activo.

```bash
gps validate general.toml
gps run general.toml --output ../../../../mi-salida/ssr-pcr   # el destino no debe existir
```

Comprobado con la semilla incluida (dura 1 s y 17 MB con el binario optimizado).
Valores que deben salir:

- 7 alelos en el microsatélite; Ho = 0,80 y He = 0,77.
- `events/reads-ssr/manifest.json`: 2817 parejas duplicadas, 1183 familias,
  fracción realizada de duplicados 0,704.
- El individuo `I001` es 21/24 unidades. Entre sus 400 lecturas, 28 cruzan el tracto
  con flancos: 16 con 21 unidades, 11 con 24 y 1 con 22 (un paso de tartamudeo).

Sirve para comparar con `ssr-sin-pcr`, que tiene la misma semilla y la misma
verdad. Para ver tartamudeo heredado por familias haría falta un homopolímero,
menos ciclos o más cobertura.
