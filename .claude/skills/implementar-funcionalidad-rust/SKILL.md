---
name: implementar-funcionalidad-rust
description: Implementar funcionalidades nuevas del simulador genómico en Rust con tipos de dominio, errores explícitos, pruebas y atención al rendimiento. Usar al añadir módulos, comportamientos, interfaces, algoritmos poblacionales, transformaciones genómicas, generación de lecturas o escritura de formatos, siempre a partir de una especificación y una arquitectura ya acordadas.
---

# Implementar funcionalidad Rust

## Objetivo

Construir código Rust correcto, legible y verificable para el núcleo del
simulador. Priorizar la fidelidad biológica y la conservación de invariantes;
optimizar después los recorridos que una medida identifique como costosos.

## Preparación

1. Leer la especificación, la arquitectura y el paso del plan aplicables.
2. Inspeccionar `Cargo.toml`, la organización de módulos, las convenciones y las
   pruebas existentes.
3. Confirmar entradas, salidas, orden, cardinalidad, errores, casos límite,
   memoria y paralelismo esperados.
4. Localizar tipos, funciones y bibliotecas que deban reutilizarse.
5. Exponer cualquier supuesto no resuelto. Detenerse solo si la elección cambia
   el comportamiento o la arquitectura; resolver localmente y registrar las
   decisiones menores.

Presentar un plan breve antes de editar cuando la tarea sea amplia. No pedir
confirmación si el usuario ya aprobó un plan suficiente y no existe una decisión
bloqueante.

## Modelado del dominio

- Representar con tipos distintos identidades, unidades y coordenadas que no
  puedan intercambiarse: referencia, individuo, evento, cromosoma, posición,
  longitud, ploidía, cobertura y semilla, entre otras.
- Usar enumeraciones para estados mutuamente excluyentes y constructores que
  validen rangos e invariantes.
- Evitar combinaciones de booleanos capaces de representar estados imposibles.
- Mantener explícitas referencia, ploidía, fase, dosis, haplotipo y procedencia
  de los segmentos.
- Impedir por tipo o validación temprana que se mezclen datos de referencias
  distintas.
- Distinguir coordenadas de configuración y de formatos externos. Centralizar
  las conversiones y probar sus extremos.
- Mantener privados los campos cuya construcción directa pueda romper una
  invariante.

## Organización del código

- Separar lógica biológica, orquestación, configuración y entrada/salida.
- Definir cada concepto compartido junto a su contrato, no dentro de una fase
  que solo lo produce.
- Mantener la visibilidad mínima: privado primero, `pub(crate)` cuando baste y
  `pub` únicamente para la interfaz prevista.
- Usar módulos hoja planos; crear directorios de módulos cuando realmente haya
  submódulos.
- Mantener una dirección de dependencias coherente con la arquitectura.
- Evitar abstracciones anticipadas. Crear un rasgo cuando exista una frontera
  estable o más de una implementación real, no solo para ocultar un tipo.

## Propiedad, memoria y rendimiento

- Preferir préstamos a clonaciones cuando el propietario y la vida útil sean
  claros.
- No forzar vidas útiles complejas para evitar una copia irrelevante; medir los
  cambios costosos de diseño.
- Representar poblaciones e historiales con estructuras compactas. Evitar una
  secuencia completa por individuo cuando puedan conservarse segmentos y
  cambios respecto a referencias o fundadores.
- Procesar por bloques o flujos cuando materializar todos los individuos,
  lecturas o registros aumente innecesariamente la memoria.
- Preasignar colecciones cuando exista un tamaño conocido o una cota útil.
- Diseñar el paralelismo con subsemillas reproducibles y sin estado mutable
  compartido oculto. El resultado biológico no debe depender del orden de los
  hilos.
- Añadir complejidad de bajo nivel solo después de medir un recorrido
  representativo en modo optimizado.

## Errores

- Usar errores tipados en las fronteras internas cuando el llamador deba
  distinguir causas.
- Añadir contexto en las fronteras de aplicación y procesos externos.
- No usar `unwrap`, `expect` o `panic!` en rutas alcanzables por entradas del
  usuario, salvo que una invariante demostrada haga el estado imposible y quede
  documentada.
- No ignorar resultados, errores parciales ni fallos de procesos externos.
- Conservar en cada diagnóstico el parámetro, archivo, referencia, muestra o
  evento necesario para corregir el problema.

## Pruebas obligatorias

Incluir, según corresponda:

- recorrido normal;
- límites, entradas vacías y valores mínimo y máximo;
- cada error estable que forme parte del contrato;
- regresiones para defectos encontrados;
- propiedades de ida y vuelta para formatos y transformaciones;
- casos sintéticos con verdad conocida para coordenadas, fase y genotipos;
- reproducibilidad con semilla y con distinto número de hilos cuando aplique;
- pruebas de integración para fronteras de archivos y procesos;
- prueba de rendimiento para un requisito cuantitativo o un recorrido crítico.

Una prueba debe usar una verdad independiente. No demostrar la corrección de un
algoritmo comparándolo únicamente con otra salida del mismo algoritmo.

## Flujo de implementación

1. Añadir o ajustar los tipos y contratos necesarios.
2. Escribir primero una prueba que exponga el comportamiento cuando sea
   práctico.
3. Implementar el cambio mínimo que cumpla el contrato.
4. Ejecutar las pruebas específicas y corregirlas sin debilitar sus
   comprobaciones.
5. Ejecutar la validación amplia disponible.
6. Revisar el diff para eliminar cambios ajenos, clonaciones innecesarias y
   supuestos no documentados.
7. Informar de cambios, pruebas, resultados reales y riesgos residuales.

## Validación

Usar los comandos y características definidos por el proyecto. Como base,
cuando estén disponibles:

```text
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-targets --all-features
cargo doc --no-deps
```

Ejecutar además las pruebas específicas del cambio. Usar `cargo audit` si está
instalado y forma parte de la política del proyecto. No instalar herramientas
ni afirmar resultados sin haberlos ejecutado.

## Finalización

La funcionalidad está terminada cuando:

- satisface los requisitos citados;
- conserva los invariantes genómicos y poblacionales;
- incluye pruebas precisas que pasan;
- no introduce cambios públicos ni dependencias sin justificar;
- supera las comprobaciones disponibles;
- deja explícitos supuestos, decisiones menores y trabajo posterior;
- el directorio de trabajo no contiene cambios accidentales.

No crear commits ni publicar cambios salvo autorización expresa del usuario o
del flujo de trabajo aprobado.
