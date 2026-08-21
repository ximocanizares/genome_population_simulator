# Revisión adversarial final del plan de implementación F2 0.4

> **Informe histórico.** Su aprobación se refería al plan anterior. El plan F2
> 0.2 simplificado necesita una nueva revisión limitada por DTC-079.

## 1. Alcance y método

- **Documento revisado:**
  [plan-implementacion-f2-0.1.md](./plan-implementacion-f2-0.1.md), en el estado
  del árbol compartido del 15 de agosto de 2026, después de DTC-075 y de las
  correcciones RAF2-041 a RAF2-044.
- **Tipo de revisión:** adversarial final, repetida desde cero y sin modificar
  el plan ni el código.
- **Fuentes contrastadas:** especificaciones, arquitectura, DTC-001 a DTC-075,
  hoja de ruta, matriz maestra de trazabilidad, registro de ampliaciones
  futuras, plan F1, punto de reanudación y las revisiones adversariales F2
  0.1, 0.2 y 0.3. El código y la CLI actuales se inspeccionaron únicamente
  para distinguir prerrequisitos F1 de trabajo futuro F2.
- **Escala:** un hallazgo `bloqueante` impide ejecutar el plan con un contrato
  único; uno `mayor` exige corregir el plan antes de aprobarlo; uno `menor` no
  invalida el recorrido, pero deja una ambigüedad verificable.

## 2. Veredicto ejecutivo

**Aprobar el plan F2 condicionado únicamente al cierre y aprobación de F1.**

No quedan hallazgos bloqueantes, mayores ni menores propios del plan F2. La
presente revisión satisface la revisión adversarial final que el propio plan
exige. Las puertas F2.0.2, F2.2.1 y F2.2.14 son pasos internos de ejecución y
no condiciones documentales adicionales para aprobar el plan.

RAF2-027 permanece abierto como **bloqueante externo**: F1 todavía no dispone
de un veredicto independiente de cierre favorable y de un commit base aprobado
e identificable. No debe ejecutarse F2.0.1 mientras falte esa condición. Esta
conclusión no convierte defectos de F1 en hallazgos propios de F2 ni obliga a
rediseñar F2 antes de conocer la línea base F1 definitiva.

## 3. Hallazgos clasificados

| Clase | Hallazgos propios abiertos | Resultado |
|---|---:|---|
| Bloqueantes | 0 | El contrato F2 es unívoco y ejecutable cuando exista una base F1 aprobada. |
| Mayores | 0 | Los pasos, propiedades y dependencias tienen responsables verificables. |
| Menores | 0 | No se detectaron ambigüedades documentales nuevas que requieran corrección. |
| Bloqueante externo | 1 | RAF2-027: cierre, aprobación y commit base F1 aún pendientes. |

### 3.1. Decisiones que requieren al usuario

**Ninguna.** No se ha encontrado una alternativa biológica abierta ni una
decisión funcional que deba resolverse antes de aprobar el plan.

### 3.2. Correcciones técnicas directas pendientes

**Ninguna en F2.** Las correcciones de RAF2-041 a RAF2-044 ya están incorporadas
y verificadas. La reparación y nueva revisión de F1 pertenecen a su propio
proceso de cierre.

## 4. Estructura, dependencias y evidencia

| Comprobación | Resultado | Valoración |
|---|---:|---|
| Pasos jerárquicos | 128 | Correcto |
| Pasos con exactamente los cinco campos obligatorios | 128/128 | Correcto |
| Campos obligatorios | `Depende de`, `Procede de`, `Cambia`, `Se comprueba con`, `Evidencia ejecutable` | Completos |
| Identificadores de paso duplicados | 0 | Correcto |
| Referencias o extremos de rango F2 inexistentes en `Depende de` | 0 | Correcto |
| Pausas `PAUSA-F2.0` a `PAUSA-F2.5` | 6 | Correcto |
| REQ de cierre exclusivo | 79 = 14+4+22+12+2+25 | Correcto |
| REQ compartidos | 9 | Correcto |
| REQ parciales directos | 9 | Correcto |

La distribución de pasos es F2.0=26, F2.1=10, F2.2=50, F2.3=17, F2.4=9 y
F2.5=16. El primer trabajo consumidor de cada hito depende de la pausa anterior
y no existe una ruta que eluda las seis revisiones obligatorias.

Cada paso declara un comando ejecutable, un expediente JSON esperado y una
regla explícita de éxito, fallo o bloqueo. Los expedientes parciales conservan
`status = "partial"`, indican el hito final y no se confunden con los 79 cierres
exclusivos. Los oráculos críticos reconstruyen o mutan artefactos sin usar el
productor como autoridad.

## 5. Trazabilidad 79+9+9

La lista de 79 REQ exclusivos coincide exactamente con las filas cuyo hito es
solo F2 en la matriz maestra; no contiene duplicados ni omisiones. Los nueve REQ
compartidos coinciden también por identidad con la matriz maestra:

`REQ-CONF-08`, `REQ-CONF-21`, `REQ-CONF-130`, `REQ-IND-02`, `REQ-IND-03`,
`REQ-IND-25`, `REQ-IND-26`, `REQ-MAN-10` y `REQ-MAN-83`.

Los nueve expedientes parciales directos son `REQ-SAL-02`, `REQ-SAL-12`,
`REQ-FMT-07`, `REQ-VER-23`, `REQ-VER-24`, `REQ-VER-25`, `REQ-VER-26`,
`REQ-MAN-31` y `REQ-MAN-58`. Todos indican productor F2 y cierre posterior.
Los tres conjuntos permanecen separados, de modo que ninguna evidencia parcial
puede aprobar un requisito con alcance posterior.

## 6. DTC-001 a DTC-075 y contratos aleatorios

El registro contiene exactamente una entrada para cada DTC-001 a DTC-075. Se
revisaron las 75 decisiones, incluidas las que no son consumidoras directas de
F2, y no se encontró una contradicción vigente con el plan. Las decisiones
aplicables están citadas por los pasos y por la matriz de trazabilidad; una DTC
de otra fase no se incorpora artificialmente al alcance F2.

### 6.1. Frecuencias y extremos

- Todo intervalo continuo es semiabierto `[min,max)`.
- Puede declararse `max=1`, pero un componente continuo nunca produce 1.
- La corrección de redondeo de `linear-uniform-v1` impide producir 0; su
  corrección superior impide alcanzar `max`.
- `log-uniform-v1` exige `min>0`, usa exclusivamente `libm = =0.2.16` y aplica
  las correcciones inferior y superior fijadas.
- Solo una fila `fixed` produce exactamente 0 o 1. Esto se refiere a la
  generación desde espectro; una matriz suministrada conserva los `alt_af`
  literales válidos en `[0,1]`.
- El perfil compacto SFS usa `[1/5008,0.005)`, `[0.005,0.05)` y `[0.05,1)`.

### 6.2. Selección SNV e indels

`weighted-without-replacement-v1` pertenece únicamente a la selección SNV
estática. `dynamic-indel-selection-v1` recibe el recuento ya sorteado, enumera
de nuevo las oportunidades después de cada evento aceptado y consume solo el
ámbito `indel_position`. Si se agota el espacio, revierte la transacción,
registra el recuento sorteado y el estado observado, publica cero eventos y no
realiza otro Poisson, recorte, reintento ni publicación parcial.

### 6.3. Propiedad única de algoritmos

F2.2.2 posee únicamente `poisson-count-v1` y evalúa la dependencia matemática
`libm`. F2.2.8 se divide en seis hijos con propiedad exclusiva: catálogo,
uniformes/enteros, selectores, distribuciones continuas, genotipo/fase y
oráculos integrados. F2.2.8.4 es el único productor de `log-uniform-v1` y
depende explícitamente de F2.2.2. Los consumidores dependen del hijo concreto
que necesitan y F2.2.8.6 busca productores duplicados.

### 6.4. `matrix_candidate_count`

`exact_loci` fija el recuento y no consume este ámbito. `mean_loci` y la
densidad convertida a media llaman una sola vez a `poisson-count-v1` con
`matrix_candidate_count`. Resultado F2.2, metadatos y manifiesto registran
algoritmo, media, ámbito, subsemilla y consumo cuando corresponden, o su
ausencia normativa en `exact_loci`.

## 7. RNG v4, modo `exact` y genotipos

El catálogo `subseed-scopes-v4.tsv` está planificado como superconjunto estricto
de v3: conserva literalmente sus seis filas y añade once ámbitos F2 con nombres,
componentes ordenados, consumidor y versión. `matrix_alternate_allele` usa solo
datos conocidos antes del ALT; `variant_id` se deriva después, por lo que no
existe circularidad. Los vectores exigen estabilidad v1-v3 y aislamiento ante
la mutación de cada componente.

La familia DTC-072 fija conversión uniforme abierta, entero acotado por rechazo,
categoría, Fisher-Yates, selectores estático y dinámico, uniformes lineal y
logarítmico, HWE, genotipo diploide con F y fase. Orden canónico, preflight,
consumo cero y versión son parte de cada contrato.

En modo `exact`, mayores restos calcula primero cuotas enteras; después, cada
población permuta independientemente el conjunto completo de `variant_id`. Un
flujo separado por población y variante obtiene el valor dentro del intervalo.
Una cuota cero y un valor `fixed` no consumen RNG. Las pruebas impiden asociar
filas por coordenada y desacoplan asignación, frecuencia, genotipo y fase.

## 8. Verdad, VCF y `GPS_CASE`

`truth.jsonl` es la verdad estructurada independiente y conserva la relación
ALT↔casos fuente. El VCF 4.5 es una representación derivada, única por conjunto,
normalizada y con una sola línea por `CHROM+POS` para todos los alelos
compatibles.

`INFO/GPS_CASE` tiene `Number=.` y `Type=String`; contiene los `case_id`
originales, únicos y en orden canónico. No se crea un caso compuesto. El plan
cubre dos ALT/dos casos, un ALT/dos casos y un caso repartido en varios
registros. F2 no fabrica pérdida ni complemento: `GPS_LOSS`, `GPS_COMPLEMENT` y
`vcf-complement.jsonl` deben estar ausentes en sus casos representables. La
primera prueba positiva real queda asignada a F5.5.

## 9. Capacidades, plataforma, Marimo y perfiles compactos

- **Capacidades:** `available_authorizations` es subconjunto de
  `catalogued_authorizations`. F2.0 no anuncia soporte efectivo. Cada
  autorización se vuelve disponible solo junto con parser, validación y
  semántica en F2.2.21 o F2.3.9; Marimo solo muestra la intersección.
- **Plataforma:** la evidencia detecta distribución, versión, kernel,
  arquitectura, triple, WSL2, montaje y CI. Ubuntu 24.04 nativo, WSL2 sobre un
  sistema de ficheros Linux y las actas manuales son expedientes distintos; una
  etiqueta discordante falla.
- **Marimo:** se mantiene como adaptador local, no autoridad biológica. Usa
  `127.0.0.1`, token obligatorio, protocolo local versionado, procesos sin
  *shell*, drenaje concurrente, cancelación escalonada y un grafo
  multidocumento importado/exportado de forma atómica y autocontenida.
- **Perfiles compactos:** los ejemplos SFS y SNV12 son recursos pequeños,
  versionados, con fuente, licencia, receta, tamaño y SHA-256. Nunca se cargan
  por omisión y no materializan ni procesan el VCF Phase 3 de 1000 Genomes.

## 10. Estado de todos los hallazgos RAF2 registrados

| Grupo | Estado final | Comprobación principal |
|---|---|---|
| RAF2-001, 004, 007 y 012 | Cerrados | Decisiones previas, evidencia parcial y destinos posteriores explícitos. |
| RAF2-013 a RAF2-026 | Cerrados en sus causas propias | Protocolo, multialelia, RNG de indels, importación, perfiles, identidad de referencia, entorno, pasos, comandos, supervisor, trazabilidad y contadores corregidos. |
| RAF2-027 | Abierto, externo | F1 debe cerrarse, aprobarse y fijar su commit base. |
| RAF2-028 a RAF2-040 | Cerrados | RNG v4, `exact`, `GPS_CASE`, parciales, pausas, algoritmos, capacidades, plataforma, subdivisión, productores, estados y nombres de artefacto verificados. |
| RAF2-041 | Cerrado | Intervalos `[min,max)`, 0/1 solo por `fixed` en generación y redondeos definidos. |
| RAF2-042 | Cerrado | Selector SNV estático separado del selector dinámico de indels. |
| RAF2-043 | Cerrado | F2.2.8.1–F2.2.8.6, dependencia explícita y productor único. |
| RAF2-044 | Cerrado | `matrix_candidate_count` consumido, aislado y registrado explícitamente. |

No existen hallazgos RAF2-002/003/005/006/008–011: la numeración histórica
conserva solo los identificadores que llegaron a registrarse y no debe
rellenarse con entradas ficticias.

## 11. Condición única antes de ejecutar F2

La revisión adversarial final del plan queda superada. Para ejecutar F2.0.1
solo falta que el proceso independiente de F1:

1. corrija y vuelva a revisar sus hallazgos;
2. apruebe F1.6 y el cierre global de F1;
3. publique un commit base identificable con pruebas, evidencias, esquemas y
   catálogos coherentes.

Si ese cierre cambia un contrato reutilizado por F2, debe revisarse el impacto
antes de ejecutar el primer paso. Esto es una salvaguarda de compatibilidad, no
una decisión F2 abierta.

## 12. Validaciones realizadas

- Lectura completa del plan F2 y de las fuentes de diseño aplicables.
- Auditoría de 128 identificadores, cinco campos por paso, seis pausas y todas
  las referencias simples o extremos de rango de `Depende de`.
- Comparación por conjuntos de los 79 REQ exclusivos y nueve compartidos con la
  matriz maestra, más revisión de los nueve parciales directos.
- Recuento y revisión de DTC-001 a DTC-075, una entrada por identificador.
- Contraste manual de DTC-070 a DTC-075 con especificación, arquitectura,
  pasos productores, consumidores, manifiestos y oráculos.
- Revisión de todos los RAF2 registrados en los informes 0.1 a 0.3.
- Comprobación final de Markdown y espacios con `git diff --check`.
