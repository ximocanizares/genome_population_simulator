import marimo

__generated_with = "0.23.16"
app = marimo.App(width="medium")

@app.cell
def _():
    import marimo as mo

    from config_editor import (
        BASES,
        DirectEventRow,
        EditorError,
        GeneralForm,
        InitialPopulationRow,
        PopulationSnvsForm,
        ReductionRow,
        ReferenceForm,
        SequenceRow,
        default_form,
        direct_event_from_editor_row,
        editor_integer,
        initial_population_from_editor_row,
        ordered_editor_rows,
        parse_general_toml,
        population_snvs_population_from_editor_row,
        render_general_toml,
    )
    from process_adapter import (
        DEFAULT_BINARY,
        PROJECT_ROOT,
        ProcessAdapterError,
        diagnostic_summary,
        initial_population_statistics_summary,
        truth_artifact_summary,
        vcf_artifact_summary,
        run_toml,
        validate_toml,
    )
    from ui_panels import (
        authored_editor_rows,
        effective_incompatible_genotype_authorization,
        explicit_snp_authorization_summary,
        individuals_panel,
    )

    return (
        BASES,
        authored_editor_rows,
        DirectEventRow,
        DEFAULT_BINARY,
        EditorError,
        GeneralForm,
        InitialPopulationRow,
        PopulationSnvsForm,
        PROJECT_ROOT,
        ProcessAdapterError,
        ReductionRow,
        ReferenceForm,
        SequenceRow,
        default_form,
        direct_event_from_editor_row,
        diagnostic_summary,
        initial_population_statistics_summary,
        truth_artifact_summary,
        vcf_artifact_summary,
        editor_integer,
        effective_incompatible_genotype_authorization,
        explicit_snp_authorization_summary,
        initial_population_from_editor_row,
        individuals_panel,
        mo,
        ordered_editor_rows,
        parse_general_toml,
        population_snvs_population_from_editor_row,
        render_general_toml,
        run_toml,
        validate_toml,
    )


@app.cell
def _(mo):
    config_upload = mo.ui.file(
        filetypes=[".toml"],
        max_size=2_000_000,
        label="Importar configuración TOML general",
    )
    return (config_upload,)


@app.cell
def _(EditorError, config_upload, default_form, parse_general_toml):
    import_error = None
    imported_name = "config.toml"
    imported_form = default_form()
    if config_upload.value:
        imported_name = config_upload.value[0].name
        try:
            imported_text = config_upload.value[0].contents.decode("utf-8")
            imported_form = parse_general_toml(imported_text)
        except (UnicodeDecodeError, EditorError) as error:
            import_error = str(error)
    return import_error, imported_form, imported_name


@app.cell
def _(BASES, imported_form, mo):
    seed_control = mo.ui.text(
        value=imported_form.seed or "",
        label="Semilla raíz opcional (64 hexadecimales minúsculos)",
        full_width=True,
    )
    ploidy_control = mo.ui.number(
        value=imported_form.ploidy,
        start=1,
        stop=10,
        step=1,
        label="Ploidía global uniforme (1–10)",
    )
    normalization_tolerance_control = mo.ui.text(
        value=imported_form.normalization_tolerance,
        label="Tolerancia general de normalización (literal decimal, 0–1e-6)",
        full_width=True,
    )
    selected_population_snvs = imported_form.population_snvs
    vcf_output_mode_control = mo.ui.dropdown(
        options=["none", "text", "bgzf"],
        value=imported_form.vcf_output_mode or "none",
        label="Vista VCF de verdad opcional",
    )
    population_snvs_mode_control = mo.ui.dropdown(
        options=["none", "supplied", "reuse", "generate"],
        value=(selected_population_snvs.mode if selected_population_snvs else "none"),
        label="Fuente de matriz SNV",
    )
    population_snvs_input_path_control = mo.ui.text(
        value=(selected_population_snvs.input_path if selected_population_snvs else "")
        or "",
        label="Ruta de tabla SNV suministrada",
        full_width=True,
    )
    population_snvs_metadata_path_control = mo.ui.text(
        value=(selected_population_snvs.metadata_path if selected_population_snvs else "")
        or "",
        label="Ruta de metadatos JSON de matriz SNV",
        full_width=True,
    )
    population_snvs_locus_count_control = mo.ui.number(
        value=((selected_population_snvs.generate or {}).get("locus_count", 10) if selected_population_snvs else 10), start=0, step=1, label="Loci candidatos exactos"
    )
    generate = (selected_population_snvs.generate or {}) if selected_population_snvs else {}
    population_snvs_count_mode_control = mo.ui.dropdown(options=["exact", "poisson"], value=generate.get("count_mode", "exact"), label="Modo de recuento")
    population_snvs_poisson_parameter_control = mo.ui.dropdown(options=["mean_loci", "mean_density_per_base"], value=("mean_density_per_base" if "mean_density_per_base" in generate else "mean_loci"), label="Parámetro Poisson")
    population_snvs_mean_loci_control = mo.ui.text(value=str(generate.get("mean_loci", "10")), label="Media Poisson de loci")
    population_snvs_density_control = mo.ui.text(value=str(generate.get("mean_density_per_base", "0.01")), label="Densidad media por base")
    population_snvs_spatial_mode_control = mo.ui.dropdown(options=["uniform", "weighted_regions"], value=generate.get("spatial_mode", "uniform"), label="Espacio de loci")
    population_snvs_substitution_profile_control = mo.ui.dropdown(
        options=["uniform_snv12_v1", "human_harris2017_snv12_approx_v1"],
        value=generate.get("substitution_profile", "uniform_snv12_v1"),
        label="Perfil de sustituciones SNV12",
    )
    population_snvs_outside_weight_control = mo.ui.text(value=str(generate.get("outside_weight", "1")), label="Peso fuera de regiones")
    population_snvs_regions_editor = mo.ui.data_editor([{"order": index, **row} for index, row in enumerate(generate.get("regions", [{"contig": "", "start": 1, "end": 1, "weight": 1.0}]), start=1)], label="Regiones ponderadas (fila vacía se omite)", pagination=False)
    population_snvs_populations_editor = mo.ui.data_editor([{"order": index, **row} for index, row in enumerate(generate.get("populations", [{"population_id": "population-1", "frequency_model": "uniform", "min_alt_af": 0.01, "max_alt_af": 0.99}]), start=1)], label="Poblaciones y perfiles de frecuencia (modo solo para spectrum/1000G)", pagination=False)
    population_snvs_spectra_editor = mo.ui.data_editor(
        [{"order": index, "population_id": population.get("population_id", ""), **component} for population in generate.get("populations", []) for index, component in enumerate(population.get("spectrum", []), start=1)] or [{"order": 1, "population_id": "population-1", "kind": "interval", "min_alt_af": 0.01, "max_alt_af": 0.99, "weight": 1.0, "sampling": "linear_uniform"}],
        label="Componentes de espectro por población (solo spectrum)",
        pagination=False,
    )
    reference_kind_control = mo.ui.dropdown(
        options=["simulated", "supplied"],
        value=imported_form.reference.kind,
        label="Origen de la referencia",
    )
    reference_id_control = mo.ui.text(
        value=imported_form.reference.id,
        label="ID de referencia",
        full_width=True,
    )
    fasta_control = mo.ui.text(
        value=imported_form.reference.fasta or "",
        label="Ruta FASTA relativa o absoluta",
        full_width=True,
    )
    composition_controls = {
        base: mo.ui.text(
            value=(imported_form.reference.composition or {}).get(base, "0.25"),
            label=f"Fracción {base} (literal decimal)",
        )
        for base in BASES
    }
    sequence_editor = mo.ui.data_editor(
        [
            {
                "order": order,
                "id": row.id,
                "class": row.sequence_class,
                "length": row.length,
            }
            for order, row in enumerate(imported_form.reference.sequences, start=1)
        ],
        label="Secuencias (edite order para reordenar)",
        pagination=False,
    )
    initial_reduction_rows = [
        {
            "order": order,
            "id": row.id,
            "source_reference_id": row.source_reference_id,
            "region_size_bases": row.region_size_bases,
            "target_kind": row.target_kind,
            "target_literal": row.target_literal,
        }
        for order, row in enumerate(imported_form.reduced_references, start=1)
    ] or [
        {
            "order": 1,
            "id": "",
            "source_reference_id": imported_form.reference.id,
            "region_size_bases": 1,
            "target_kind": "target_bases",
            "target_literal": "",
        }
    ]
    reduction_editor = mo.ui.data_editor(
        initial_reduction_rows,
        label="Referencias reducidas (la fila vacía no se exporta)",
        pagination=False,
    )
    initial_population_rows = [
        {
            "order": order,
            "population_id": row.population_id,
            "individual_count": row.individual_count,
            "ploidy": row.ploidy,
            "genotype_model": row.genotype_model,
            "inbreeding_coefficient": row.inbreeding_coefficient or "",
            "individual_ids": ", ".join(row.individual_ids or []),
            "founder_keys": ", ".join(row.founder_keys or []),
            "individual_ids_specified": row.individual_ids is not None,
            "omit_if_blank": False,
        }
        for order, row in enumerate(imported_form.initial_populations, start=1)
    ] or [{
        "order": 1,
        "population_id": "",
        "individual_count": 1,
        "ploidy": 2,
        "genotype_model": "hwe",
        "inbreeding_coefficient": "",
        "individual_ids": "",
        "founder_keys": "",
        "individual_ids_specified": False,
        "omit_if_blank": True,
    }]
    initial_populations_editor = mo.ui.data_editor(
        initial_population_rows,
        label="Cohortes iniciales (`omit_if_blank` permite omitir solo el borrador vacío)",
        pagination=False,
    )
    direct_event_rows = [
        {
            "order": order,
            "event_id": row.event_id,
            "individual_id": row.individual_id,
            "phase": row.phase,
            "kind": row.kind,
            "contig": row.contig,
            "position_1based": row.position_1based or "",
            "ref": row.reference or "",
            "alt": row.alternate or "",
            "point": row.point if row.point is not None else "",
            "bases": row.bases or "",
            "start": row.start if row.start is not None else "",
            "end": row.end if row.end is not None else "",
            "omit_if_blank": False,
        }
        for order, row in enumerate(imported_form.direct_events, start=1)
    ] or [{
        "order": 1,
        "event_id": "",
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
        "omit_if_blank": True,
    }]
    direct_events_editor = mo.ui.data_editor(
        direct_event_rows,
        label="Cambios directos (`omit_if_blank` permite omitir solo el borrador vacío; complete las columnas de kind)",
        pagination=False,
    )
    return (
        composition_controls,
        fasta_control,
        reduction_editor,
        reference_id_control,
        reference_kind_control,
        seed_control,
        ploidy_control,
        normalization_tolerance_control,
        population_snvs_input_path_control,
        population_snvs_metadata_path_control,
        population_snvs_mode_control,
        population_snvs_locus_count_control,
        population_snvs_count_mode_control,
        population_snvs_poisson_parameter_control,
        population_snvs_mean_loci_control,
        population_snvs_density_control,
        population_snvs_spatial_mode_control,
        population_snvs_substitution_profile_control,
        population_snvs_outside_weight_control,
        population_snvs_regions_editor,
        population_snvs_populations_editor,
        population_snvs_spectra_editor,
        initial_populations_editor,
        direct_events_editor,
        vcf_output_mode_control,
        sequence_editor,
    )


@app.cell
def _(
    BASES,
    authored_editor_rows,
    EditorError,
    GeneralForm,
    PopulationSnvsForm,
    ReductionRow,
    ReferenceForm,
    SequenceRow,
    composition_controls,
    direct_event_from_editor_row,
    direct_events_editor,
    editor_integer,
    fasta_control,
    reduction_editor,
    reference_id_control,
    reference_kind_control,
    render_general_toml,
    ploidy_control,
    normalization_tolerance_control,
    population_snvs_input_path_control,
    population_snvs_metadata_path_control,
    population_snvs_mode_control,
    population_snvs_locus_count_control,
    population_snvs_count_mode_control,
    population_snvs_poisson_parameter_control,
    population_snvs_mean_loci_control,
    population_snvs_density_control,
    population_snvs_spatial_mode_control,
    population_snvs_substitution_profile_control,
    population_snvs_outside_weight_control,
    population_snvs_regions_editor,
    population_snvs_populations_editor,
    population_snvs_spectra_editor,
    initial_population_from_editor_row,
    initial_populations_editor,
    seed_control,
    sequence_editor,
    vcf_output_mode_control,
    ordered_editor_rows,
    population_snvs_population_from_editor_row,
):
    preview_error = None
    preview_toml = ""
    def build_population_snvs_generate():
        count_mode = str(population_snvs_count_mode_control.value)
        generated = {
            "count_mode": count_mode,
            "spatial_mode": str(population_snvs_spatial_mode_control.value),
            "substitution_profile": str(population_snvs_substitution_profile_control.value),
            "populations": [],
        }
        components_by_population = {}
        for row in ordered_editor_rows(population_snvs_spectra_editor.value, "population_snvs.generate.populations.spectrum"):
            population_id = str(row.get("population_id", "")).strip()
            if population_id:
                component = {key: row[key] for key in ("kind", "fixed_alt_af", "min_alt_af", "max_alt_af", "weight", "sampling") if key in row and row[key] is not None}
                components_by_population.setdefault(population_id, []).append(component)
        for row in ordered_editor_rows(population_snvs_populations_editor.value, "population_snvs.generate.populations"):
            population_id = str(row["population_id"]).strip()
            if not population_id:
                continue
            generated["populations"].append(
                population_snvs_population_from_editor_row(row, components_by_population)
            )
        if count_mode == "exact":
            generated["locus_count"] = editor_integer(population_snvs_locus_count_control.value, "population_snvs.generate.locus_count")
        else:
            if str(population_snvs_poisson_parameter_control.value) == "mean_loci":
                generated["mean_loci"] = float(population_snvs_mean_loci_control.value)
            else:
                generated["mean_density_per_base"] = float(population_snvs_density_control.value)
        if generated["spatial_mode"] == "weighted_regions":
            generated["outside_weight"] = float(population_snvs_outside_weight_control.value)
            generated["regions"] = [
                {"contig": str(row["contig"]), "start": editor_integer(row["start"], "population_snvs.generate.regions.start"), "end": editor_integer(row["end"], "population_snvs.generate.regions.end"), "weight": float(row["weight"])}
                for row in ordered_editor_rows(population_snvs_regions_editor.value, "population_snvs.generate.regions")
                if str(row["contig"]).strip()
            ]
        return generated
    try:
        ordered_sequences = ordered_editor_rows(
            sequence_editor.value, "reference.sequences"
        )
        sequence_rows = [
            SequenceRow(
                id=str(row["id"]),
                sequence_class=str(row["class"]),
                length=editor_integer(row["length"], "reference.sequences.length"),
            )
            for row in ordered_sequences
        ]
        active_reductions = [row for row in reduction_editor.value if str(row["id"]).strip()]
        ordered_reductions = ordered_editor_rows(
            active_reductions, "reduced_references"
        )
        authored_reductions = [
            ReductionRow(
                id=str(row["id"]),
                source_reference_id=str(row["source_reference_id"]),
                region_size_bases=editor_integer(
                    row["region_size_bases"], "reduced_references.region_size_bases"
                ),
                target_kind=str(row["target_kind"]),
                target_literal=str(row["target_literal"]),
            )
            for row in ordered_reductions
        ]
        active_initial_populations = authored_editor_rows(
            initial_populations_editor.value, "population_id"
        )
        authored_initial_populations = [
            initial_population_from_editor_row(row, f"initial_populations[{index}]")
            for index, row in enumerate(
                ordered_editor_rows(active_initial_populations, "initial_populations"), start=0
            )
        ]
        active_direct_events = authored_editor_rows(
            direct_events_editor.value, "event_id"
        )
        authored_direct_events = [
            direct_event_from_editor_row(row, f"direct_events[{index}]")
            for index, row in enumerate(
                ordered_editor_rows(active_direct_events, "direct_events"), start=0
            )
        ]
        selected_kind = str(reference_kind_control.value)
        selected_population_mode = str(population_snvs_mode_control.value)
        population_snvs = (
            None
            if selected_population_mode == "none"
            else PopulationSnvsForm(
                mode=selected_population_mode,
                input_path=(
                    population_snvs_input_path_control.value
                    if selected_population_mode == "supplied"
                    else None
                ),
                metadata_path=(
                    population_snvs_metadata_path_control.value
                    if selected_population_mode == "reuse"
                    else None
                ),
                generate=(build_population_snvs_generate() if selected_population_mode == "generate" else None),
            )
        )
        edited_form = GeneralForm(
            seed=seed_control.value or None,
            ploidy=editor_integer(ploidy_control.value, "ploidy"),
            normalization_tolerance=normalization_tolerance_control.value,
            reference=ReferenceForm(
                kind=selected_kind,
                id=reference_id_control.value,
                fasta=fasta_control.value if selected_kind == "supplied" else None,
                composition=(
                    {base: composition_controls[base].value for base in BASES}
                    if selected_kind == "simulated"
                    else None
                ),
                sequences=sequence_rows,
            ),
            reduced_references=authored_reductions,
            population_snvs=population_snvs,
            initial_populations=authored_initial_populations,
            direct_events=authored_direct_events,
            vcf_output_mode=(
                None
                if str(vcf_output_mode_control.value) == "none"
                else str(vcf_output_mode_control.value)
            ),
        )
        preview_toml = render_general_toml(edited_form)
    except (EditorError, KeyError, TypeError, ValueError) as error:
        preview_error = str(error)
    return preview_error, preview_toml


@app.cell
def _(import_error, preview_error, preview_toml):
    authoring_error = import_error or preview_error
    authoring_toml = "" if import_error else preview_toml
    return authoring_error, authoring_toml


@app.cell
def _(DEFAULT_BINARY, PROJECT_ROOT, authoring_error, mo):
    allow_empty_snv_matrix_control = mo.ui.checkbox(
        value=False,
        label="Autorizar una matriz SNV activa vacía",
    )
    allow_incompatible_genotypes_control = mo.ui.checkbox(
        value=False,
        label="Autorizar genotipos explícitos incompatibles con el modelo",
    )
    binary_path_control = mo.ui.text(
        value=str(DEFAULT_BINARY),
        label="Binario Rust",
        full_width=True,
    )
    base_directory_control = mo.ui.text(
        value=str(PROJECT_ROOT),
        label="Directorio base de las rutas relativas del TOML",
        full_width=True,
    )
    output_directory_control = mo.ui.text(
        value=str(PROJECT_ROOT / "results" / "marimo-run"),
        label="Directorio de resultados (debe no existir)",
        full_width=True,
    )
    validate_button = mo.ui.run_button(
        label="Validar con Rust",
        kind="success",
        disabled=bool(authoring_error),
        keyboard_shortcut="Ctrl-Shift-V",
    )
    return (
        allow_empty_snv_matrix_control,
        allow_incompatible_genotypes_control,
        base_directory_control,
        binary_path_control,
        output_directory_control,
        validate_button,
    )


@app.cell
def _(
    ProcessAdapterError,
    allow_empty_snv_matrix_control,
    allow_incompatible_genotypes_control,
    authoring_toml,
    base_directory_control,
    binary_path_control,
    direct_events_editor,
    effective_incompatible_genotype_authorization,
    mo,
    validate_button,
    validate_toml,
):
    from pathlib import Path as _Path

    validation_error = None
    validation_result = None
    validated_session = None
    if validate_button.value:
        try:
            _allow_incompatible_genotypes = effective_incompatible_genotype_authorization(
                direct_events_editor.value,
                allow_incompatible_genotypes_control.value,
            )
            with mo.status.spinner(title="Validando con el binario Rust…"):
                validation_result, validated_session = validate_toml(
                    _Path(binary_path_control.value),
                    authoring_toml,
                    _Path(base_directory_control.value),
                    allow_empty_snv_matrix=allow_empty_snv_matrix_control.value,
                    allow_incompatible_genotypes=_allow_incompatible_genotypes,
                )
        except ProcessAdapterError as error:
            validation_error = str(error)
    return validation_error, validation_result, validated_session


@app.cell
def _(authoring_error, mo, validated_session):
    run_button = mo.ui.run_button(
        label="Ejecutar simulación",
        kind="warn",
        disabled=bool(authoring_error) or validated_session is None,
        keyboard_shortcut="Ctrl-Shift-R",
    )
    return (run_button,)


@app.cell
def _(
    ProcessAdapterError,
    allow_empty_snv_matrix_control,
    allow_incompatible_genotypes_control,
    authoring_toml,
    base_directory_control,
    binary_path_control,
    direct_events_editor,
    effective_incompatible_genotype_authorization,
    mo,
    output_directory_control,
    run_button,
    run_toml,
    validated_session,
):
    from pathlib import Path as _Path

    run_error = None
    run_result = None
    if run_button.value and validated_session is not None:
        try:
            _allow_incompatible_genotypes = effective_incompatible_genotype_authorization(
                direct_events_editor.value,
                allow_incompatible_genotypes_control.value,
            )
            with mo.status.spinner(title="Ejecutando el simulador local…"):
                run_result = run_toml(
                    validated_session,
                    authoring_toml,
                    _Path(binary_path_control.value),
                    _Path(base_directory_control.value),
                    _Path(output_directory_control.value),
                    allow_empty_snv_matrix=allow_empty_snv_matrix_control.value,
                    allow_incompatible_genotypes=_allow_incompatible_genotypes,
                )
        except ProcessAdapterError as error:
            run_error = str(error)
    return run_error, run_result


@app.cell
def _(
    allow_empty_snv_matrix_control,
    allow_incompatible_genotypes_control,
    authoring_error,
    authoring_toml,
    base_directory_control,
    binary_path_control,
    composition_controls,
    config_upload,
    effective_incompatible_genotype_authorization,
    fasta_control,
    import_error,
    imported_name,
    explicit_snp_authorization_summary,
    initial_population_statistics_summary,
    initial_populations_editor,
    individuals_panel,
    mo,
    output_directory_control,
    reduction_editor,
    reference_id_control,
    reference_kind_control,
    run_button,
    run_error,
    run_result,
    seed_control,
    ploidy_control,
    normalization_tolerance_control,
    population_snvs_input_path_control,
    population_snvs_metadata_path_control,
    population_snvs_mode_control,
    population_snvs_locus_count_control,
    population_snvs_count_mode_control,
    population_snvs_poisson_parameter_control,
    population_snvs_mean_loci_control,
    population_snvs_density_control,
    population_snvs_spatial_mode_control,
    population_snvs_substitution_profile_control,
    population_snvs_outside_weight_control,
    population_snvs_regions_editor,
    population_snvs_populations_editor,
    population_snvs_spectra_editor,
    sequence_editor,
    direct_events_editor,
    validate_button,
    validation_error,
    validation_result,
    diagnostic_summary,
    vcf_artifact_summary,
    vcf_output_mode_control,
):
    visible_error = authoring_error
    visible_preview = authoring_toml
    general_panel = mo.vstack(
        [
            mo.md("## General/Ejecución"),
            config_upload,
            seed_control,
            (
                mo.callout(import_error, kind="danger", title="Error de importación")
                if import_error
                else mo.callout(
                    "Importar o editar no valida ni ejecuta la simulación.", kind="info"
                )
            ),
        ]
    )
    selected_source_controls = (
        fasta_control
        if reference_kind_control.value == "supplied"
        else mo.hstack([composition_controls[base] for base in ("A", "C", "G", "T")])
    )
    genome_panel = mo.vstack(
        [
            mo.md("## Genoma"),
            mo.hstack([reference_kind_control, reference_id_control]),
            selected_source_controls,
            sequence_editor,
            reduction_editor,
        ]
    )
    population_panel_mode = str(population_snvs_mode_control.value)
    population_path_control = (
        population_snvs_input_path_control
        if population_panel_mode == "supplied"
        else (
            population_snvs_metadata_path_control
            if population_panel_mode == "reuse"
            else mo.vstack([
                population_snvs_count_mode_control,
                population_snvs_poisson_parameter_control,
                population_snvs_locus_count_control,
                population_snvs_mean_loci_control,
                population_snvs_density_control,
                population_snvs_spatial_mode_control,
                population_snvs_substitution_profile_control,
                population_snvs_outside_weight_control,
                population_snvs_regions_editor,
                population_snvs_populations_editor,
                population_snvs_spectra_editor,
            ])
            if population_panel_mode == "generate"
            else mo.callout("No se ha seleccionado una matriz SNV.", kind="neutral")
        )
    )
    populations_panel = mo.vstack(
        [
            mo.md("## Poblaciones"),
            population_snvs_mode_control,
            allow_empty_snv_matrix_control,
            normalization_tolerance_control,
            population_path_control,
            mo.callout(
                "La interfaz solo conserva la forma estructural; Rust valida la semántica.",
                kind="info",
            ),
        ]
    )
    vcf_mode = str(vcf_output_mode_control.value)
    outputs_panel = mo.vstack(
        [
            mo.md("## Salidas"),
            vcf_output_mode_control,
            mo.callout(
                "No se solicitará VCF; truth.jsonl sigue siendo la verdad obligatoria."
                if vcf_mode == "none"
                else (
                    "Se solicitará un único truth.vcf sin índice."
                    if vcf_mode == "text"
                    else "Se solicitará un único truth.vcf.gz con exactamente un índice TBI o CSI."
                ),
                kind="info",
            ),
        ]
    )
    unavailable = lambda phase: mo.callout(
        f"Sin controles todavía: esta capacidad se incorpora en {phase}.", kind="neutral"
    )
    if validation_error:
        validation_status = mo.callout(
            validation_error, kind="danger", title="No se pudo validar"
        )
    elif validation_result is None:
        validation_status = mo.callout(
            "Pendiente de validación con Rust.", kind="neutral"
        )
    elif validation_result.succeeded:
        validation_status = mo.callout(
            "Configuración validada por el binario Rust.", kind="success"
        )
    else:
        validation_status = mo.callout(
            validation_result.stdout or validation_result.stderr,
            kind="danger",
            title=f"Validación rechazada (código {validation_result.returncode})",
        )
    if run_error:
        run_status = mo.callout(run_error, kind="danger", title="No se pudo ejecutar")
    elif run_result is None:
        run_status = mo.callout("Ejecución no iniciada.", kind="neutral")
    elif run_result.succeeded:
        run_status = mo.callout(
            "Ejecución completada. "
            f"Resultados: {run_result.output_directory}. "
            f"Manifiesto: {run_result.output_directory / 'manifest.json'}.",
            kind="success",
        )
    else:
        run_status = mo.callout(
            run_result.stdout or run_result.stderr,
            kind="danger",
            title=f"Ejecución fallida (código {run_result.returncode})",
        )
    process_logs = ""
    for process_result in (validation_result, run_result):
        if process_result is not None:
            process_logs += process_result.stdout + process_result.stderr
    if len(process_logs) > 20_000:
        process_logs = "[logs truncados a los últimos 20000 caracteres]\n" + process_logs[-20_000:]
    localized_diagnostics = "\n".join(
        summary
        for summary in (
            diagnostic_summary(
                validation_result.diagnostics if validation_result is not None else None
            ),
            diagnostic_summary(run_result.diagnostics if run_result is not None else None),
        )
        if summary
    )
    founders_summary = (
        initial_population_statistics_summary(run_result.output_directory)
        if run_result is not None and run_result.succeeded and run_result.output_directory is not None
        else "El resumen F2.3.5 aparecerá tras una ejecución correcta."
    )
    truth_summary = (
        truth_artifact_summary(run_result.output_directory)
        if run_result is not None and run_result.succeeded and run_result.output_directory is not None
        else "El resumen de truth.jsonl aparecerá tras una ejecución correcta."
    )
    vcf_summary = (
        vcf_artifact_summary(run_result.output_directory)
        if run_result is not None and run_result.succeeded and run_result.output_directory is not None
        else "El resumen del VCF aparecerá tras una ejecución correcta."
    )
    incompatible_authorization = (
        allow_incompatible_genotypes_control
        if effective_incompatible_genotype_authorization(
            direct_events_editor.value,
            True,
        )
        else None
    )
    incompatible_authorization_summary = explicit_snp_authorization_summary(
        direct_events_editor.value
    )
    review_panel = mo.vstack(
        [
            mo.md("## Revisar y ejecutar"),
            (
                mo.callout(visible_error, kind="danger", title="Error de autoría")
                if visible_error
                else mo.callout(
                    "Previsualización estructural correcta. Use «Validar con Rust» "
                    "antes de ejecutar.",
                    kind="success",
                )
            ),
            mo.ui.text_area(
                value=visible_preview,
                label="TOML general previsualizado",
                rows=20,
                disabled=True,
                full_width=True,
            ),
            mo.download(
                visible_preview,
                filename=imported_name or "config.toml",
                mimetype="application/toml",
                disabled=bool(visible_error),
                label="Exportar TOML general",
            ),
            mo.md("### Validar y ejecutar"),
            binary_path_control,
            base_directory_control,
            output_directory_control,
            mo.hstack([validate_button, run_button]),
            validation_status,
            run_status,
            mo.ui.text_area(
                value=founders_summary,
                label="Resumen de fundadores F2.3.5 (sin cargar secuencias)",
                rows=5,
                disabled=True,
                full_width=True,
            ),
            mo.ui.text_area(
                value=truth_summary,
                label="Resumen de truth.jsonl (sin cargar el JSONL)",
                rows=3,
                disabled=True,
                full_width=True,
            ),
            mo.ui.text_area(
                value=vcf_summary,
                label="Resumen del VCF del manifiesto (sin cargar el VCF)",
                rows=3,
                disabled=True,
                full_width=True,
            ),
            mo.ui.text_area(
                value=localized_diagnostics,
                label="Diagnósticos localizados por pestaña y campo",
                rows=5,
                disabled=True,
                full_width=True,
            ),
            mo.ui.text_area(
                value=process_logs,
                label="Diagnósticos y logs del proceso",
                rows=10,
                disabled=True,
                full_width=True,
            ),
        ]
    )
    application_tabs = mo.ui.tabs(
        {
            "General/Ejecución": general_panel,
            "Genoma": genome_panel,
            "Individuos": individuals_panel(
                mo,
                ploidy_control,
                initial_populations_editor,
                direct_events_editor,
                incompatible_authorization,
                incompatible_authorization_summary,
            ),
            "Poblaciones": populations_panel,
            "Evolución": unavailable("F7"),
            "Familias": unavailable("F8"),
            "Eventos técnicos": unavailable("F3–F5"),
            "Salidas": outputs_panel,
            "Revisar y ejecutar": review_panel,
        }
    )
    mo.vstack(
        [
            mo.md("# Genome Population Simulator"),
            mo.callout(
                "Interfaz local en desarrollo. Rust conserva toda autoridad semántica.",
                kind="info",
            ),
            application_tabs,
        ]
    )
    return


if __name__ == "__main__":
    app.run()
