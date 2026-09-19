# Invariantes del simulador genómico

## Referencias e identidad

- Asociar cada secuencia, coordenada, alineamiento, variante e índice con una
  referencia identificada y versionada.
- Impedir mezclar BAM, CRAM o VCF correspondientes a referencias diferentes.
- Mantener identificadores biológicos estables y separar identificadores
  técnicos de evento, biblioteca, tanda, carril y grupo de lectura.
- Verificar cabeceras, diccionarios de secuencia, orden e índices antes de marcar
  una salida como completa.

## Coordenadas y estructura

- Centralizar la conversión entre intervalos de configuración y convenciones de
  FASTA, SAM/BAM/CRAM y VCF.
- Probar primera base, última base, intervalo de una base, inserciones, límites
  terminales permitidos y solapamientos.
- Interpretar cada cambio respecto al estado que indique la especificación y
  conservar su correspondencia con la referencia o segmento fundador.
- No recolocar silenciosamente una variante cuyo objetivo ya no existe.

## Genotipos, fase y ploidía

- Conservar por cromosoma la ploidía común de la ejecución y validar el límite
  configurado.
- Mantener dosis, fase y copia haplotípica durante reproducción, recombinación y
  mutación.
- Reconstruir el linaje hasta fundadores sin depender del VCF.
- Construir el VCF de verdad desde el estado genotípico final, nunca desde
  lecturas o alineamientos.
- No incluir en el VCF verdadero medidas inferidas como profundidad o calidad de
  genotipo.

## Evolución y aleatoriedad

- Aplicar los procesos en el orden total configurado.
- Mantener separadas N objetivo, N realizado, Ne objetivo y Ne realizado.
- Ejecutar migración, selección, reproducción, mutación y eventos demográficos
  con reglas y ámbitos explícitos.
- Derivar subsemillas de forma estable por fase, generación, población, evento y
  muestra según la estrategia acordada.
- Confirmar que observaciones opcionales y el número de hilos no cambian la
  trayectoria biológica.

## Lecturas y alineamientos

- Generar lecturas desde haplotipos en proporción a dosis y ploidía antes de
  aplicar sesgos.
- Mantener juntas las parejas `paired-end` en sustitución, contaminación y
  demultiplexado.
- Separar contaminación biológica de asignación errónea técnica.
- No reutilizar FASTQ entre eventos que deban ser independientes.
- Verificar que el BAM o CRAM refleja razonablemente la verdad sin convertir esa
  evidencia técnica en genotipo verdadero.

## Trazabilidad y resultados

- Registrar configuración solicitada y efectiva, semilla, versiones,
  referencias, herramientas, parámetros y archivos.
- Conservar estados parciales o avisos sin presentar una ejecución incompleta
  como válida.
- Hacer reproducibles los sorteos y registrar cantidades esperadas y realizadas.
- Validar sumas de comprobación y recursos de exportaciones autocontenidas
  cuando esa capacidad se implemente.
