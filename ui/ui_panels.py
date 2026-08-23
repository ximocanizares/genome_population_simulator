"""Composición reutilizable de los paneles disponibles de la interfaz local."""


_INITIAL_PLACEHOLDER_DEFAULTS = {
    "individual_count": 1,
    "ploidy": 2,
    "genotype_model": "hwe",
    "inbreeding_coefficient": "",
    "individual_ids": "",
    "founder_keys": "",
    "individual_ids_specified": False,
}
_DIRECT_EVENT_PLACEHOLDER_DEFAULTS = {
    "individual_id": "",
    "phase": 0,
    "kind": "snp",
    "contig": "",
    "position_1based": 1,
    "ref": "A",
    "alt": "G",
    "point": "",
    "bases": "",
    "start": "",
    "end": "",
}


def is_pristine_placeholder(row, identifier_field):
    """Only omit a row explicitly marked as a still-unedited UI placeholder."""
    if not bool(row.get("omit_if_blank", False)) or str(row.get(identifier_field, "")).strip():
        return False
    defaults = (
        _INITIAL_PLACEHOLDER_DEFAULTS
        if identifier_field == "population_id"
        else _DIRECT_EVENT_PLACEHOLDER_DEFAULTS
        if identifier_field == "event_id"
        else {}
    )
    return all(row.get(field, expected) == expected for field, expected in defaults.items())


def is_explicit_snp_row(row):
    """Return whether an authored row requests an explicit SNP."""
    return not is_pristine_placeholder(row, "event_id") and str(row.get("kind", "")) == "snp"


def authored_editor_rows(rows, identifier_field):
    """Drop only rows explicitly marked as blank UI placeholders."""
    return [
        row
        for row in rows
        if not is_pristine_placeholder(row, identifier_field)
    ]


def effective_incompatible_genotype_authorization(rows, requested):
    """Keep the narrow CLI authorization off unless an explicit SNP can consume it."""
    return bool(requested) and any(is_explicit_snp_row(row) for row in rows)


def explicit_snp_authorization_summary(rows):
    """Describe the bounded set that the narrow authorization could affect."""
    snps = [row for row in rows if is_explicit_snp_row(row)]
    if not snps:
        return "No hay SNP explícitos: la autorización de genotipos incompatibles no se aplicará."
    entries = []
    for row in snps[:20]:
        entries.append(
            "{individual}@{contig}:{position}".format(
                individual=str(row.get("individual_id", "")).strip() or "(individual_id vacío)",
                contig=str(row.get("contig", "")).strip() or "(contig vacío)",
                position=str(row.get("position_1based", "")).strip() or "(posición vacía)",
            )
        )
    suffix = f"; … y {len(snps) - 20} más" if len(snps) > 20 else ""
    return "SNP explícitos afectados por esta autorización: " + ", ".join(entries) + suffix


def individuals_panel(
    mo,
    ploidy_control,
    initial_populations_editor=None,
    direct_events_editor=None,
    allow_incompatible_genotypes_control=None,
    incompatible_genotypes_summary=None,
):
    """Render F2.3 controls only; Rust remains the biological authority."""
    children = [mo.md("## Individuos"), ploidy_control]
    if initial_populations_editor is not None:
        children.extend(
            [
                mo.md("### Cohortes iniciales"),
                initial_populations_editor,
                mo.callout(
                    "Use `individual_ids` o, para IDs derivados, una lista de "
                    "`founder_keys` separada por comas; ambas columnas son excluyentes.",
                    kind="info",
                ),
            ]
        )
    if direct_events_editor is not None:
        children.extend(
            [
                mo.md("### Cambios explícitos por individuo"),
                direct_events_editor,
            ]
        )
    if incompatible_genotypes_summary is not None:
        children.append(mo.callout(incompatible_genotypes_summary, kind="info"))
    if allow_incompatible_genotypes_control is not None:
        children.append(allow_incompatible_genotypes_control)
    return mo.vstack(children)
