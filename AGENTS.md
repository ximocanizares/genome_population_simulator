> **HISTÓRICO. No es autoridad.**
>
> Este fichero corresponde a la etapa en que el proyecto se trabajaba con
> Codex, y sus niveles (`gpt-5.6-luna`, `terra`, `sol`) ya no se usan. La
> política de agentes vigente está en **`CLAUDE.md`**, que lo sustituye por
> completo. Se conserva solo como registro de cómo se decidía antes.

# Política de agentes del proyecto

Estas instrucciones se aplican a todos los agentes que trabajen en este
repositorio. Complementan las skills del proyecto; no sustituyen sus puertas,
autoridades ni pausas obligatorias.

## Enrutamiento de modelos y cuota

El agente raíz debe elegir de forma explícita el modelo y el esfuerzo de cada
subagente. No se usa siempre el modelo heredado.

Antes de crear un subagente:

1. comprobar si la tarea puede resolverse directamente con una búsqueda, un
   comando o una edición pequeña; en ese caso no crear el subagente;
2. acotar una responsabilidad independiente y verificable;
3. elegir el nivel mínimo que preserve la calidad necesaria;
4. enviar un contexto autocontenido y reducido;
5. indicar al usuario qué nivel se eligió cuando la elección sea material.

| Tipo de trabajo | Modelo predeterminado | Esfuerzo predeterminado |
|---|---|---|
| Estado, inventario, búsquedas acotadas, ejecución de comandos, formato y documentación mecánica | `gpt-5.6-luna` | `low` o `medium` |
| Exploración concreta del código | `gpt-5.6-luna` | `medium` |
| Implementación ordinaria, correcciones verificadas, pruebas, planificación y documentación técnica | `gpt-5.6-terra` | `medium` |
| Implementación compleja entre módulos o corrección de un defecto grave | `gpt-5.6-terra` | `high` |
| Arquitectura material, decisiones biológicas/estadísticas difíciles, revisión adversarial o revisión independiente de cierre | `gpt-5.6-sol` | `high` |

`xhigh`, `max`, `ultra` o modos equivalentes solo se usan si el usuario lo pide
o si una evaluación reproducible demuestra una mejora necesaria. No se elevan
por precaución genérica.

### Escalado

Escalar de Terra a Sol únicamente cuando ocurra al menos una condición:

- existe una contradicción entre especificación, arquitectura y código;
- está en juego una invariante científica, corrupción o pérdida de datos;
- una corrección grave ha fallado o reaparecido tras una revisión;
- el plan exige una revisión independiente de cierre;
- el usuario prioriza expresamente máxima calidad sobre cuota y latencia.

Tras resolver el punto crítico, volver a Terra o Luna para correcciones,
documentación, gates y seguimiento.

## Contexto y concurrencia

- No heredar todo el chat por defecto. Usar `fork_turns="none"` o el menor
  número de turnos que permita trabajar, acompañado de una tarea
  autocontenida con archivos, contratos, límites y resultado esperado.
- Usar historial completo solo cuando la tarea dependa realmente de decisiones
  conversacionales que no estén registradas en el repositorio.
- No mantener más de un subagente Sol activo a la vez.
- No lanzar dos agentes para la misma pregunta. Las revisiones independientes
  consecutivas se ejecutan secuencialmente y sobre el estado corregido.
- No crear un agente únicamente para esperar, consultar estado o ejecutar un
  único comando.

## Modo de conservación de cuota

Codex no expone al agente el porcentaje semanal de cuota. El agente raíz usa el
último porcentaje comunicado por el usuario:

- con más del 50 % restante, aplicar la tabla normal;
- entre 25 % y 50 %, reservar Sol para revisiones obligatorias, decisiones
  materiales y bloqueos; usar Terra para el resto;
- con menos del 25 %, no iniciar agentes opcionales y pedir confirmación antes
  de abrir un nuevo agente Sol que no sea obligatorio para una pausa ya
  aprobada.

La conservación de cuota no permite omitir pruebas, revisiones obligatorias ni
reglas de seguridad. Puede reducir paralelismo, contexto y nivel de modelo, o
aplazar trabajo opcional.

## Registro mínimo

Al crear un subagente, el agente raíz debe dejar en su encargo:

- objetivo y propiedad exactos;
- modelo y esfuerzo elegidos;
- contexto heredado mínimo;
- archivos o módulos en alcance;
- prohibición de revertir cambios ajenos;
- gates, pausa y acciones externas permitidas;
- condición de escalado o detención.

Si el modelo elegido no está disponible, no sustituir silenciosamente una tarea
rutinaria por Sol alto. Usar el siguiente nivel eficiente disponible o informar
al usuario cuando la sustitución cambie materialmente cuota, calidad o plazo.
