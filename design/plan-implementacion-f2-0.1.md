# Plan de implementación de F2 0.2

## 1. Cabecera y estado

- **Título:** implementación simplificada de F2, verdad genómica mínima.
- **Versión:** 0.2.
- **Fecha:** 21 de agosto de 2026.
- **Estado:** simplificado tras detectar que la versión 0.1 repetía el
  sobrediseño retirado de F1; pendiente de revisión y aprobación.
- **Base:** F1 cerrada en `e14705a` y documentada en `efca265`.
- **Alcance:** F2.0–F2.5.
- **Autoridad:** especificación 0.2, arquitectura 0.2, DTC-079, hoja de ruta y
  trazabilidad vigentes.

Este plan sustituye los 130 pasos de la versión 0.1. Conserva las decisiones
biológicas, estadísticas y de formato ya aprobadas, pero elimina infraestructura
anticipada, expedientes por paso y defensas frente a manipulación local activa.
Los cambios F2.0.1 aún no comprometidos del plan anterior no se consideran
implementación aprobada y deben retirarse o adaptarse al ejecutar el primer paso.

## 2. Resultado de F2

F2 termina con un recorrido directo local que:

- parte de una referencia cerrada en F1;
- crea individuos con ploidía uniforme entre 1 y 10;
- mantiene copias haplotípicas separadas y reconstruibles;
- admite o genera una matriz poblacional de frecuencias SNV;
- aplica SNP e indels pequeños explícitos o aleatorios;
- publica `truth.jsonl` y exactamente un VCF de verdad multimuestra por
  conjunto;
- permite configurar, validar, ejecutar y consultar resultados desde Marimo en
  `127.0.0.1`.

La verdad, la fase, la dosis, la procedencia y las coordenadas se derivan del
estado genómico, nunca del VCF ni de lecturas.

## 3. Límites de simplificación

### 3.1. Incluido

- Un único TOML general para la vía directa. El identificador interno reservado
  `scenario_id = "direct"` se usa para identidad y RNG, no como entidad
  biológica ni documento configurable.
- Marimo como interfaz local de configuración y resultados: controles de las
  capacidades implementadas, importación/exportación del TOML, `validate`,
  `run`, diagnósticos, progreso y enlaces o resúmenes.
- Estado haplotípico compacto, matriz SNV, perfiles aproximados aprobados,
  individuos HWE o `F` diploide, verdad estructurada y VCF 4.5.
- Validación de entradas externas, errores explícitos, reproducibilidad,
  publicación no-clobber y pruebas de invariantes biológicas y de formatos.

### 3.2. Excluido

- Escenarios, configuración multidocumento y preservación opaca de capacidades
  futuras; se incorporan con su primer consumidor real.
- Ejecutores o registros generales de evidencia, un `result.json` por paso,
  protocolos paralelos de capacidades y almacenes persistentes de sesión UI.
- Defensa frente a una persona con escritura simultánea en el repositorio o el
  provisional, symlinks hostiles, funciones de entorno inyectadas y carreras
  adversarias.
- Auditoría completa de navegador, plataforma oficial, `.tar.gz`, OCI y
  endurecimiento de despliegue; pertenecen a F9.
- Lecturas, alineamientos, variantes estructurales, dinámica poblacional,
  familias, cruces, recombinación y paneles haplotípicos; pertenecen a fases
  posteriores.
- Distribuciones probabilísticas de longitud de indel; F2 usa una longitud
  exacta por clase.

## 4. Forma de trabajo y evidencia

Cada paso produce código o documentación que utiliza el siguiente paso. No se
crea un sistema adicional para demostrarlo. La evidencia es:

1. prueba focal del comportamiento;
2. puertas Cargo/Python aplicables;
3. validadores externos solo para interoperabilidad real;
4. revisión del diff y una nota breve de cierre por hito.

Puertas comunes:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --all-features --locked -- -D warnings
cargo test --workspace --all-targets --all-features --locked
cargo test --workspace --doc --all-features --locked
cargo doc --workspace --no-deps --all-features --locked
git diff --check
```

Marimo añade las pruebas declaradas en `pyproject.toml` mediante `uv run`.
Los validadores VCF se ejecutan en F2.5. No se exige que cada paso cree un
artefacto de evidencia independiente.

## 5. F2.0 — Marimo local y configuración directa

**Resultado:** la interfaz local arranca, edita el TOML general vigente, invoca
el binario y presenta sus diagnósticos y resultados.

- ☐ F2.0.1. Exponer a Marimo la configuración F1 ya implementada.
  - Depende de: F1 cerrada.
  - Procede de: REQ-UI-01 a REQ-UI-03, DTC-051 y DTC-079.
  - Cambia: una frontera de lectura del esquema y los diagnósticos F1 existentes
    para que Marimo pueda representarlos; no añade ploidía, población, vía ni
    escenarios.
  - Se comprueba con: los TOML válidos e inválidos F1 producen mediante la UI
    el mismo documento o diagnóstico que por CLI, y los campos F2 siguen
    rechazados hasta su hito propietario.

- ☐ F2.0.2. Incorporar el entorno Marimo bloqueado y un lanzador local pequeño.
  - Depende de: F2.0.1.
  - Procede de: REQ-UI-01, REQ-UI-13, DTC-051 y DTC-079.
  - Cambia: `pyproject.toml`, `uv.lock`, aplicación y lanzador que enlaza
    `127.0.0.1`, muestra la URL y admite puerto local y `--no-open`.
  - Se comprueba con: aprovisionamiento inicial mediante `uv sync --frozen`,
    que puede descargar únicamente las dependencias resueltas por `uv.lock`;
    arranque posterior sin red, puerto ocupado y confirmación de que abrir no
    ejecuta. La instalación distribuible sin descarga previa se valida en F9.

- ☐ F2.0.3. Implementar edición e importación/exportación del TOML general.
  - Depende de: F2.0.1–F2.0.2.
  - Procede de: REQ-UI-02 a REQ-UI-06, DTC-051 y DTC-079.
  - Cambia: pestañas previstas, controles solo para capacidades implementadas,
    tablas necesarias y previsualización TOML; no crea escenarios.
  - Se comprueba con: ciclo TOML→UI→TOML→Rust equivalente, orden estable,
    nombres con espacios y errores del esquema visibles.

- ☐ F2.0.4. Conectar `validate` y `run` con argumentos separados.
  - Depende de: F2.0.3.
  - Procede de: REQ-UI-07 a REQ-UI-12 y DTC-079.
  - Cambia: comprobación simple de versión mediante `--version`, acción
    explícita, invalidación tras editar, destino visible, diagnósticos, progreso
    y enlaces/resúmenes finales; no crea un protocolo de capacidades ni hashes
    entre componentes del mismo árbol.
  - Se comprueba con: proceso falso y binario real; edición posterior a validar,
    versión incompatible, destino existente, fallo y artefacto grande no
    cargado.

- ☐ F2.0.5. Cerrar la paridad y usabilidad local de F2.0.
  - Depende de: F2.0.2–F2.0.4.
  - Procede de: REQ-UI-01 a REQ-UI-13, REQ-VAL-54 y ACE-114/115.
  - Cambia: pruebas de recorrido, etiquetas, teclado, foco y estados no basados
    solo en color; nota de cierre sin afirmar distribución oficial.
  - Se comprueba con: configuración F1/F2 representativa produce el mismo
    `EffectiveConfig` por UI y TOML manual; arranque local y puertas comunes.

> **PAUSA-F2.0.** Revisar que Marimo siga siendo una interfaz local del binario,
> que no existan escenarios anticipados ni infraestructura de seguridad o
> evidencia general. No iniciar F2.1 sin aprobar esta pausa.

## 6. F2.1 — Estado haplotípico mínimo

**Resultado:** cada copia haplotípica conserva fase, procedencia y ediciones y
puede materializarse exactamente desde la referencia.

- ☐ F2.1.1. Extender el servicio común de coordenadas a variantes pequeñas.
  - Depende de: PAUSA-F2.0 y coordenadas F1.6.
  - Procede de: DTC-010, DTC-013, DTC-056; REQ-IND-38/39.
  - Cambia: posición, intervalo, ancla de inserción y deleción interna, todos
    etiquetados por referencia y secuencia.
  - Se comprueba con: vectores base 1 cerrada↔base 0 semiabierta, extremos,
    cero, desbordamiento y secuencia incorrecta.

- ☐ F2.1.2. Modelar referencia, individuo y copias con ploidía uniforme.
  - Depende de: F2.1.1.
  - Procede de: DTC-003, DTC-032, DTC-039, DTC-061;
    REQ-IND-28/29 y REQ-VER-27.
  - Cambia: configuración y tipos de base, individuo, copia y fase para ploidía
    1..=10, preservando `N`; es el primer consumidor real de ploidía F2.
  - Se comprueba con: ploidías 1/2/10, cardinalidad, IDs únicos, rechazo 0/11 y
    referencia con `N` en bordes e interior.

- ☐ F2.1.3. Representar segmentos, modificaciones y procedencia sin copiar la referencia.
  - Depende de: F2.1.1–F2.1.2.
  - Procede de: DTC-031 y DTC-061; REQ-IND-28/29.
  - Cambia: segmentos que referencian fuentes y modificaciones tipadas SNP,
    inserción y deleción; solo contenido nuevo de inserción posee bytes.
  - Se comprueba con: copia intacta y cada tipo de edición, origen, coordenadas,
    compartición de segmentos y composición corta conocida.

- ☐ F2.1.4. Materializar copias y derivar fase y dosis.
  - Depende de: F2.1.3.
  - Procede de: DTC-003, DTC-031 y DTC-061; REQ-VER-27.
  - Cambia: materializador determinista y consultas genotípicas desde el estado
    canónico.
  - Se comprueba con: cadenas manuales haploides, diploides y poliploides; dosis
    igual al número de copias portadoras.

- ☐ F2.1.5. Integrar y cerrar el estado haplotípico.
  - Depende de: F2.1.1–F2.1.4.
  - Procede de: REQ-REF-14, REQ-IND-28/29, REQ-VER-27 y ACE-03/10/11/28/66.
  - Cambia: pruebas de integración y documentación del límite frente a F5.3.
  - Se comprueba con: reconstrucción exacta, conservación de `N`, consumo de
    memoria razonable y puertas comunes.

> **PAUSA-F2.1.** No continuar si una copia no puede reconstruirse exactamente,
> si fase/dosis proceden de una salida o si existen dos convenios internos de
> coordenadas.

## 7. F2.2 — Matriz SNV y variantes pequeñas

**Resultado:** F2 admite o genera una matriz común reproducible y aplica SNP e
indels pequeños sin confundir frecuencias, alelos, genotipos o recuentos.

- ☐ F2.2.1. Definir la configuración y los formatos mínimos de la matriz.
  - Depende de: PAUSA-F2.1.
  - Procede de: DTC-052, DTC-054, DTC-060; REQ-CONF-128/129 y REQ-IND-60.
  - Cambia: `population_snvs.tsv`, metadatos JSON y TOML excluyente para
    matriz suministrada, reutilizada o generada.
  - Se comprueba con: fixtures válidos/inválidos, SHA-256, referencia,
    poblaciones, `alt_af` única y rechazo de VCF como entrada.

- ☐ F2.2.2. Incorporar normalización y parámetros efectivos usados por F2.
  - Depende de: F2.2.1.
  - Procede de: DTC-049 y DTC-079; REQ-CONF-106/107, REQ-VAL-40 y REQ-MAN-81.
  - Cambia: normalización canónica de sumas alélicas, tolerancia, procedencia y
    registro en configuración efectiva/manifiesto.
  - Se comprueba con: sumas dentro/fuera del borde, no finitos, cero y oráculo
    numérico independiente.

- ☐ F2.2.3. Añadir los ámbitos y utilidades RNG realmente consumidos por F2.
  - Depende de: F2.2.1–F2.2.2 y RNG F1.
  - Procede de: DTC-037/038, DTC-067, DTC-070 a DTC-072, DTC-075/079.
  - Cambia: catálogo v4 y módulos agrupados para Poisson, selección de loci,
    espectros e indels; introduce aquí el `scenario_id = "direct"` interno. Las
    utilidades de genotipo y fase esperan a F2.3.
  - Se comprueba con: vectores de subsemilla, semillas fijas, bordes, consumo
    cero, aislamiento entre población/individuo/clase y propiedades estadísticas.

- ☐ F2.2.4. Validar y preparar matrices suministradas.
  - Depende de: F2.2.1–F2.2.3.
  - Procede de: DTC-054 y DTC-060; REQ-IND-60/61 y REQ-VAL-55.
  - Cambia: orden, IDs, filas multialélicas, REF/ALT, poblaciones y eliminación
    de ALT globalmente ausentes sin modificar el original.
  - Se comprueba con: duplicados, reordenación, referencia errónea, sumas,
    matriz activa vacía y autorización específica.

- ☐ F2.2.5. Generar loci y frecuencias de una matriz nueva.
  - Depende de: F2.2.1–F2.2.3.
  - Procede de: DTC-052, DTC-054, DTC-067/070/072;
    REQ-CONF-24/25/129 y REQ-IND-17 a REQ-IND-19/61.
  - Cambia: bases elegibles, regiones ponderadas, recuento exacto o Poisson,
    selección sin reemplazo, ALT y frecuencias por población.
  - Se comprueba con: `N` excluida, espacio vacío/insuficiente, tasa cero,
    mayores restos, independencia poblacional y semilla repetible.

- ☐ F2.2.6. Incorporar espectros y perfiles aproximados aprobados.
  - Depende de: F2.2.2–F2.2.5.
  - Procede de: DTC-068 y DTC-070/072; REQ-CONF-129 y ACE-116.
  - Cambia: espectros suministrados, `uniform_snv12_v1`,
    `human_1000g_global_sfs_approx_v1` y
    `human_harris2017_snv12_approx_v1`.
  - Se comprueba con: intervalos semiabiertos/fijos, 12 sustituciones, recursos
    pequeños versionados y reproducción independiente de sus literales.

- ☐ F2.2.7. Publicar y reutilizar la matriz canónica.
  - Depende de: F2.2.4–F2.2.6.
  - Procede de: DTC-054; REQ-MAN-91 y REQ-VAL-55.
  - Cambia: generación y preparación producen el mismo par TSV+JSON y la
    ejecución puede reutilizarlo sin reescribirlo.
  - Se comprueba con: round-trip semántico, SHA-256, orden canónico y ejecución
    repetida con iguales bytes/resultados.

- ☐ F2.2.8. Aplicar SNP explícitos y derivados de la matriz.
  - Depende de: F2.1.3–F2.1.4 y F2.2.4–F2.2.7.
  - Procede de: REQ-CONF-08/21/28, REQ-IND-05/06/10/11/13/14.
  - Cambia: constructores y aplicación trazable de SNP sobre posiciones
    elegibles, sin `REF=N` ni otra fuente aleatoria por individuo.
  - Se comprueba con: REF/ALT, regiones, conflictos, procedencia y matriz común.

- ☐ F2.2.9. Aplicar inserciones y deleciones pequeñas.
  - Depende de: F2.1.3–F2.1.4 y F2.2.3.
  - Procede de: DTC-059 y DTC-075; REQ-IND-38/40/41 y REQ-REC-24.
  - Cambia: longitudes exactas, oportunidades sobre estado vigente, contenido
    A/C/G/T y aplicación con anclas/bases flanqueantes.
  - Se comprueba con: longitudes 1, 2 y umbral−1; terminales rechazados,
    reenumeración dinámica, `N`, agotamiento y tasa cero.

- ☐ F2.2.10. Componer eventos pequeños y cerrar diagnósticos.
  - Depende de: F2.2.8–F2.2.9.
  - Procede de: DTC-031/032/044; REQ-DIR-03 y ACE-01/09/21/22/28.
  - Cambia: orden declarado, conflictos explícitos, preflight y ausencia de
    publicación parcial.
  - Se comprueba con: SNP eliminado por deleción, inserciones vecinas,
    solapamientos, entrada incompatible y cadenas cortas con verdad manual.

- ☐ F2.2.11. Exponer matriz y variantes pequeñas en Marimo.
  - Depende de: F2.2.7–F2.2.10.
  - Procede de: REQ-UI-14, DTC-051/079 y ACE-114/116.
  - Cambia: controles de fuente, generación, perfiles, regiones, SNP/indels,
    previsualización y autorización de matriz vacía.
  - Se comprueba con: UI→TOML→Rust equivalente y previsualización coherente.

- ☐ F2.2.12. Cerrar integración, oráculos y trazabilidad F2.2.
  - Depende de: F2.2.1–F2.2.11.
  - Procede de: requisitos y ACE asignados a F2.2.
  - Cambia: pruebas integradas y nota de cierre; no añade otro ejecutor.
  - Se comprueba con: matriz suministrada/generada, cambios explícitos/aleatorios,
    semillas, oráculos independientes pertinentes y puertas comunes.

> **PAUSA-F2.2.** Revisar especialmente coordenadas, `N`, sumas, aislamiento
> RNG, matriz común y composición de indels. No iniciar individuos si la matriz
> no es reproducible y validable de forma independiente.

## 8. F2.3 — Individuos directos

**Resultado:** la población de origen produce individuos faseados con HWE o
`F` diploide, manteniendo objetivos y realizaciones separados.

- ☐ F2.3.1. Extender configuración e identidades de población e individuo.
  - Depende de: PAUSA-F2.2.
  - Procede de: DTC-003/039/052; REQ-CONF-20/22/38 y REQ-IND-15/25.
  - Cambia: población de origen, número de individuos, ploidía, HWE/`F` e IDs
    estables; no crea estructura poblacional evolutiva.
  - Se comprueba con: IDs aportados/derivados, unicidad, P=1..10 y `F` solo P=2.

- ☐ F2.3.2. Muestrear composición genotípica HWE o `F`.
  - Depende de: F2.2.3–F2.2.7 y F2.3.1.
  - Procede de: DTC-052/071/072; REQ-IND-62/63.
  - Cambia: utilidad RNG y composición por locus e individuo desde `alt_af`,
    sin Poisson.
  - Se comprueba con: casos degenerados, fórmulas independientes, ploidías
    1/2/10 y estadísticas prefijadas.

- ☐ F2.3.3. Asignar fase con un ámbito independiente.
  - Depende de: F2.3.2.
  - Procede de: DTC-037/071/072; REQ-IND-63.
  - Cambia: utilidad RNG y permutación de copias posterior a la composición.
  - Se comprueba con: dosis conservada, aislamiento genotipo/fase y semilla fija.

- ☐ F2.3.4. Integrar SNP explícitos e indels en las copias individuales.
  - Depende de: F2.2.8–F2.2.10 y F2.3.2–F2.3.3.
  - Procede de: REQ-DIR-01/02 y REQ-IND-64/65.
  - Cambia: validación de incompatibilidades y autorización específica, seguida
    por cambios directos locales.
  - Se comprueba con: incompatibilidad rechazada/confirmada, ninguna corrección
    silenciosa y ausencia de `--force`.

- ☐ F2.3.5. Materializar individuos y calcular estadísticas solicitadas/realizadas.
  - Depende de: F2.3.4.
  - Procede de: REQ-IND-62 a REQ-IND-65 y REQ-MAN-91.
  - Cambia: secuencias finales, procedencia, frecuencias realizadas y resumen.
  - Se comprueba con: reconstructor manual, varios individuos, adición/reordenación
    aislada y frecuencias finitas sin autocorrección.

- ☐ F2.3.6. Añadir los controles de individuos a Marimo.
  - Depende de: F2.3.1–F2.3.5.
  - Procede de: REQ-UI-14 y DTC-051/079.
  - Cambia: población, N, ploidía, HWE/`F`, autorización y resúmenes.
  - Se comprueba con: paridad UI→TOML→Rust y diagnósticos localizados.

- ☐ F2.3.7. Cerrar integración y trazabilidad F2.3.
  - Depende de: F2.3.1–F2.3.6.
  - Procede de: requisitos y ACE asignados a F2.3.
  - Cambia: pruebas y nota de cierre.
  - Se comprueba con: HWE/`F`, fase, dosis, múltiples individuos, estadísticas
    y puertas comunes.

> **PAUSA-F2.3.** No continuar si la fase o dosis no se reconstruyen desde las
> copias, si `F` se aplica fuera de diploides o si Marimo obtiene resultados
> distintos del TOML manual.

## 9. F2.4 — Verdad estructurada

**Resultado:** `truth.jsonl` permite reconstruir cada copia sin consultar VCF.

- ☐ F2.4.1. Definir el esquema externo mínimo de `truth.jsonl`.
  - Depende de: PAUSA-F2.3.
  - Procede de: DTC-061; REQ-VER-31 y REQ-SAL-06.
  - Cambia: cabecera y registros de copias/genotipos con versión, orden,
    identidad, fase, dosis y procedencia.
  - Se comprueba con: fixtures manuales válidos, campos ausentes/duplicados,
    orden, versiones y relaciones rotas.

- ☐ F2.4.2. Proyectar y escribir la verdad por flujo.
  - Depende de: F2.4.1 y estado F2.3.
  - Procede de: DTC-004/031/061.
  - Cambia: proyección pura y escritor JSONL mediante publicación F1.
  - Se comprueba con: fichero determinista, recuentos, SHA-256, truncado y fallo.

- ☐ F2.4.3. Reconstruir y validar la verdad independientemente.
  - Depende de: F2.4.2.
  - Procede de: DTC-047/061 y ACE-34/62.
  - Cambia: validador de esquema y reconstructor semántico separado de la
    proyección.
  - Se comprueba con: secuencias y genotipos conocidos, mutaciones aisladas y
    funcionamiento con VCF desactivado.

- ☐ F2.4.4. Integrar verdad, manifiesto y presentación Marimo.
  - Depende de: F2.4.2–F2.4.3.
  - Procede de: REQ-SAL-02, REQ-IND-26 y DTC-051/079.
  - Cambia: enlace de artefacto, resumen y acceso desde la UI sin carga completa.
  - Se comprueba con: identidad/ruta/tamaño/SHA, enlace roto y resumen pequeño.

- ☐ F2.4.5. Cerrar integración y trazabilidad F2.4.
  - Depende de: F2.4.1–F2.4.4.
  - Procede de: requisitos y ACE asignados a F2.4.
  - Cambia: pruebas y nota de cierre.
  - Se comprueba con: reconstrucción sin VCF, corrupción y puertas comunes.

> **PAUSA-F2.4.** No iniciar VCF si `truth.jsonl` no es independiente,
> reconstruible y derivado únicamente del estado canónico.

## 10. F2.5 — VCF de verdad mínimo

**Resultado:** exactamente un VCF 4.5 multimuestra representa la verdad F2 con
GT faseado, dosis y procedencia correctos.

- ☐ F2.5.1. Incorporar el contrato VCF y sus validadores.
  - Depende de: PAUSA-F2.4.
  - Procede de: DTC-062/063/066; REQ-FMT-05/06/07.
  - Cambia: dependencia/escritor aprobados y validación externa
    `bcftools`/HTSlib solo en pruebas.
  - Se comprueba con: evaluación de dependencia, texto, BGZF+TBI y BGZF+CSI.

- ☐ F2.5.2. Proyectar y normalizar alelos desde la verdad.
  - Depende de: F2.4 y F2.5.1.
  - Procede de: DTC-064; REQ-VER-15/16/17/30/44.
  - Cambia: intervalo REF común, recorte, alineamiento izquierdo y orden
    determinista.
  - Se comprueba con: SNP, inserción, deleción, solapamiento y comparación con
    normalizador externo sin repararlo.

- ☐ F2.5.3. Agrupar alelos y calcular GT faseado.
  - Depende de: F2.5.2.
  - Procede de: DTC-064/065; REQ-VER-02/06/12/32 a REQ-VER-41/47.
  - Cambia: un registro multialélico por posición F2, ALT y GT de todas las
    muestras, dosis y `GPS_CASE`.
  - Se comprueba con: varios ALT/casos/copias, ploidías 1/2/10, orden, duplicados
    y dosis manual.

- ☐ F2.5.4. Aplicar `N` y el complemento condicional.
  - Depende de: F2.5.2–F2.5.3.
  - Procede de: DTC-032/065; REQ-VER-23 a REQ-VER-26.
  - Cambia: omisión o representación válida y ausencia del complemento cuando
    F2 no pierde detalle.
  - Se comprueba con: `N` en REF/ALT, casos representables y campos/fichero
    complementario ausentes.

- ☐ F2.5.5. Escribir, indexar y publicar el VCF.
  - Depende de: F2.5.3–F2.5.4.
  - Procede de: DTC-063/066; REQ-MAN-27/57/58 y REQ-SAL-12/21.
  - Cambia: exactamente un VCF, BGZF+TBI/CSI por defecto o texto explícito,
    manifiesto y unidad de publicación.
  - Se comprueba con: lectura completa, consulta regional, índice único,
    destino existente, fallo de índice y SHA-256.

- ☐ F2.5.6. Validar semánticamente VCF contra verdad y referencia.
  - Depende de: F2.5.5.
  - Procede de: DTC-066; REQ-VAL-28.
  - Cambia: comparador independiente de FASTA+`truth.jsonl`+VCF.
  - Se comprueba con: REF/ALT/GT/fase/dosis/muestras/`GPS_CASE`, mutaciones
    aisladas y prueba de que VCF nunca alimenta la verdad.

- ☐ F2.5.7. Integrar VCF en Marimo y cerrar F2.
  - Depende de: F2.5.1–F2.5.6.
  - Procede de: DTC-051/079, requisitos y ACE asignados a F2.5.
  - Cambia: controles de salida, resumen/enlace, recorrido completo y nota de
    cierre.
  - Se comprueba con: UI/TOML equivalentes, configuración→verdad→VCF,
    validadores externos y puertas comunes.

> **PAUSA-F2.5 / CIERRE F2.** Revisar de forma independiente las invariantes
> genómicas, la trazabilidad y el recorrido local completo. No convertir la
> revisión en una auditoría frente a atacante local ni reintroducir
> infraestructura retirada por DTC-079.

## 11. Trazabilidad de cierre

La asignación por hito permanece sin cambiar el significado biológico:

| Hito | REQ de cierre exclusivo | Evidencia principal |
|---|---:|---|
| F2.0 | 13 | Arranque local, TOML directo, paridad UI→Rust y ejecución explícita |
| F2.1 | 3 | Copias, fase, dosis, procedencia y reconstrucción |
| F2.2 | 22 | Matriz SNV, perfiles, SNP/indels y reproducibilidad |
| F2.3 | 13 | Individuos HWE/`F`, fase y estadísticas |
| F2.4 | 2 | `truth.jsonl` independiente y reconstruible |
| F2.5 | 25 | VCF 4.5 multimuestra validado contra la verdad |
| **Total** | **78** | |

Los catorce REQ compartidos conservan sus destinos: REQ-CONF-08,
REQ-CONF-21, REQ-CONF-106, REQ-CONF-107, REQ-VAL-40, REQ-MAN-81, REQ-IND-02,
REQ-IND-03, REQ-CONF-130, REQ-IND-25, REQ-IND-26, REQ-MAN-10, REQ-MAN-83 y
REQ-UI-13. F2 aporta solo la parte observable de su alcance.

Los nueve REQ con cierre posterior conservan evidencia parcial: REQ-SAL-02,
REQ-SAL-12, REQ-FMT-07, REQ-VER-23, REQ-VER-24, REQ-VER-25, REQ-VER-26,
REQ-MAN-31 y REQ-MAN-58. Su prueba en F2 no autoriza marcarlos como cerrados.

La matriz maestra mantiene la autoridad sobre los identificadores concretos y
los destinos finales. La simplificación reduce tareas y expedientes, no rebaja
ningún requisito biológico ni de interoperabilidad.

## 12. Criterio de salida del plan

El plan puede aprobarse cuando:

- especificación, arquitectura, decisiones y hoja de ruta expresan el mismo
  modelo local confiable;
- no queda ninguna referencia operativa al ejecutor `verify-f2-step`;
- los pasos no anticipan escenarios ni componentes de fases posteriores;
- cada paso tiene dependencia, origen, cambio y prueba;
- existen seis pausas, una por hito;
- la trazabilidad conserva 78 cierres exclusivos, 14 compartidos y 9 parciales;
- una revisión comprueba proporcionalidad y contratos científicos, sin ampliar
  otra vez el modelo de amenaza.
