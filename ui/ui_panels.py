"""Composición reutilizable de los paneles disponibles de la interfaz local."""


def individuals_panel(mo, ploidy_control):
    """Controles disponibles para la pestaña de individuos en el alcance actual."""
    return mo.vstack(
        [
            mo.md("## Individuos"),
            ploidy_control,
        ]
    )
