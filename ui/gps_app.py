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
        REDUCTION_OBJECTIVE_KEEP_BASES,
        ReductionRow,
        ReferenceForm,
        RepeatAdmissionForm,
        SequenceRow,
        default_form,
        direct_event_from_editor_row,
        editor_integer,
        initial_population_from_editor_row,
        ordered_editor_rows,
        parse_general_toml,
        population_snvs_population_from_editor_row,
        read_fasta_sequence_rows,
        reduction_editor_row,
        reduction_from_editor_row,
        render_general_toml,
        repeat_allele_editor_rows,
        repeat_interruption_editor_rows,
        repeat_loci_from_editor_rows,
        repeat_locus_editor_rows,
    )
    from f5_editor import (
        F5GeneralForm,
        F5Pair,
        F5ScenarioForm,
        F5_SAMPLING_MODES,
        F5_TARGET_KINDS,
        available_sampling_modes,
        default_f5_pair,
        f5_override_editor_rows,
        f5_override_from_editor_row,
        f5_pair_documents,
        f5_subpopulation_editor_rows,
        f5_subpopulation_from_editor_row,
        import_f5_pair,
        inbreeding_is_available,
        pair_without_inbreeding,
    )
    from process_adapter import (
        DEFAULT_BINARY,
        PROJECT_ROOT,
        ProcessAdapterError,
        diagnostic_summary,
        f5_package_summary,
        initial_population_statistics_summary,
        truth_artifact_summary,
        vcf_artifact_summary,
        run_f5_pair,
        run_toml,
        validate_f5_pair,
        validate_toml,
    )
    from ui_panels import (
        authored_editor_rows,
        effective_incompatible_genotype_authorization,
        explicit_snp_authorization_summary,
        individuals_panel,
        repeats_panel,
    )

    return (
        BASES,
        authored_editor_rows,
        DirectEventRow,
        DEFAULT_BINARY,
        EditorError,
        F5GeneralForm,
        F5Pair,
        F5ScenarioForm,
        F5_SAMPLING_MODES,
        F5_TARGET_KINDS,
        available_sampling_modes,
        default_f5_pair,
        f5_override_editor_rows,
        f5_override_from_editor_row,
        f5_package_summary,
        f5_pair_documents,
        f5_subpopulation_editor_rows,
        f5_subpopulation_from_editor_row,
        import_f5_pair,
        inbreeding_is_available,
        pair_without_inbreeding,
        run_f5_pair,
        validate_f5_pair,
        GeneralForm,
        InitialPopulationRow,
        PopulationSnvsForm,
        REDUCTION_OBJECTIVE_KEEP_BASES,
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
        RepeatAdmissionForm,
        repeats_panel,
        repeat_allele_editor_rows,
        repeat_interruption_editor_rows,
        repeat_loci_from_editor_rows,
        repeat_locus_editor_rows,
        parse_general_toml,
        population_snvs_population_from_editor_row,
        read_fasta_sequence_rows,
        reduction_editor_row,
        reduction_from_editor_row,
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
    load_default_config_button = mo.ui.run_button(
        label="Cargar configuración predeterminada del proyecto",
    )
    return config_upload, load_default_config_button


@app.cell
def _(PROJECT_ROOT, mo):
    base_directory_control = mo.ui.text(
        value=str(PROJECT_ROOT),
        label="Directorio base de las rutas relativas del TOML",
        full_width=True,
    )
    base_directory_browser = mo.ui.file_browser(
        initial_path=PROJECT_ROOT,
        selection_mode="directory",
        multiple=False,
        label="O seleccione el directorio base local",
    )
    return base_directory_browser, base_directory_control


@app.cell
def _(
    EditorError,
    PROJECT_ROOT,
    config_upload,
    default_form,
    load_default_config_button,
    parse_general_toml,
):
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
    elif load_default_config_button.value:
        default_config_path = PROJECT_ROOT / "tests/fixtures/config/general-minimal.toml"
        imported_name = str(default_config_path.relative_to(PROJECT_ROOT))
        try:
            imported_form = parse_general_toml(default_config_path.read_text())
        except (OSError, EditorError) as error:
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
    selected_technical_event = imported_form.technical_event
    technical_event_mode_control = mo.ui.dropdown(
        options=["none", "configured"],
        value="configured" if selected_technical_event else "none",
        label="Evento técnico paired-end F4",
    )
    technical_event_id_control = mo.ui.text(
        value=selected_technical_event.event_id if selected_technical_event else "",
        label="ID del evento técnico",
        full_width=True,
    )
    technical_event_individual_ids_control = mo.ui.text(
        value=", ".join(selected_technical_event.individual_ids)
        if selected_technical_event
        else "",
        label="IDs de individuos (separados por comas)",
        full_width=True,
    )
    technical_event_coverage_control = mo.ui.text(
        value=selected_technical_event.coverage_wgs if selected_technical_event else "0",
        label="Cobertura WGS objetivo (× haploide)",
        full_width=True,
    )
    technical_event_read_length_control = mo.ui.text(
        value=str(selected_technical_event.read_length_bases)
        if selected_technical_event and selected_technical_event.read_length_bases is not None
        else "",
        label="Longitud de lectura opcional (1–150; vacío = 150)",
        full_width=True,
    )
    technical_library_rows = [
        {
            "order": order,
            "individual_id": library.individual_id,
            "i7": library.i7,
            "i5": library.i5,
            "omit_if_blank": False,
        }
        for order, library in enumerate(
            selected_technical_event.libraries if selected_technical_event else [], start=1
        )
    ] or [{
        "order": 1,
        "individual_id": "",
        "i7": "",
        "i5": "",
        "omit_if_blank": True,
    }]
    technical_libraries_editor = mo.ui.data_editor(
        technical_library_rows,
        label="Bibliotecas UDI (i7/i5: exactamente 8 bases A/C/G/T; library_id se deriva en Rust)",
        pagination=False,
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
    fasta_browser = mo.ui.file_browser(
        initial_path=PROJECT_ROOT,
        filetypes=[".fa", ".fasta", ".fna", ".fas"],
        selection_mode="file",
        multiple=False,
        label="O seleccione un FASTA local (puede navegar fuera del proyecto)",
    )
    composition_controls = {
        base: mo.ui.text(
            value=(imported_form.reference.composition or {}).get(
                base, "0" if base == "N" else "0.25"
            ),
            label=(
                "Fracción N (posiciones ambiguas, independiente)"
                if base == "N"
                else f"Fracción {base}"
            ),
        )
        for base in (*BASES, "N")
    }
    read_fasta_headers_button = mo.ui.run_button(
        label="Leer cabeceras y longitudes del FASTA",
    )
    initial_reduction_rows = [
        reduction_editor_row(row, order)
        for order, row in enumerate(imported_form.reduced_references, start=1)
    ] or [
        {
            "order": 1,
            "id": "",
            "source_reference_id": imported_form.reference.id,
            "region_size_bases": 1,
            "Objetivo": REDUCTION_OBJECTIVE_KEEP_BASES,
            "Valor": "",
        }
    ]
    reduction_editor = mo.ui.data_editor(
        initial_reduction_rows,
        label="Reducción de referencia (la fila sin ID no se exporta)",
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
        fasta_browser,
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
        technical_event_mode_control,
        technical_event_id_control,
        technical_event_individual_ids_control,
        technical_event_coverage_control,
        technical_event_read_length_control,
        technical_libraries_editor,
        vcf_output_mode_control,
        read_fasta_headers_button,
    )


@app.cell
def _(
    PROJECT_ROOT,
    EditorError,
    base_directory_control,
    fasta_browser,
    fasta_control,
    imported_form,
    mo,
    read_fasta_headers_button,
    read_fasta_sequence_rows,
):
    fasta_headers_error = None
    _fasta_sequence_rows = imported_form.reference.sequences
    if read_fasta_headers_button.value:
        declared_path = (
            fasta_browser.path(0)
            if fasta_browser.value
            else fasta_control.value.strip()
        )
        from pathlib import Path

        fasta_path = Path(declared_path).expanduser()
        if not fasta_path.is_absolute():
            base_directory = Path(base_directory_control.value.strip() or PROJECT_ROOT)
            fasta_path = base_directory / fasta_path
        try:
            _fasta_sequence_rows = read_fasta_sequence_rows(fasta_path)
        except EditorError as error:
            fasta_headers_error = str(error)
    sequence_editor = mo.ui.data_editor(
        [
            {
                "order": order,
                "id": row.id,
                "class": row.sequence_class,
                "length": row.length,
            }
            for order, row in enumerate(_fasta_sequence_rows, start=1)
        ],
        label="Secuencias (edite order para reordenar)",
        pagination=False,
    )
    return fasta_headers_error, sequence_editor


@app.cell
def _(
    imported_form,
    mo,
    repeat_allele_editor_rows,
    repeat_interruption_editor_rows,
    repeat_locus_editor_rows,
):
    _imported_loci = imported_form.reference.repeat_loci
    _imported_admission = imported_form.reference.repeat_admission
    repeat_loci_editor = mo.ui.data_editor(
        repeat_locus_editor_rows(_imported_loci)
        or [
            {
                "order": 1,
                "id": "",
                "sequence": "",
                "start": 1,
                "end": 1,
                "motifs": "",
                "max_units": "",
                "min_purity": "",
                "min_longest_run_units": "",
            }
        ],
        label=(
            "Loci de repetición (motivo y posición; una fila sin `id` se omite."
            " `motifs` admite varios separados por comas)"
        ),
        pagination=False,
    )
    repeat_alleles_editor = mo.ui.data_editor(
        repeat_allele_editor_rows(_imported_loci)
        or [
            {
                "order": 1,
                "locus_id": "",
                "units": "",
                "sequence": "",
                "frequency": "1.0",
            }
        ],
        label=(
            "Alelos por locus (`units` con un valor por bloque, o `sequence`;"
            " son excluyentes)"
        ),
        pagination=False,
    )
    repeat_interruptions_editor = mo.ui.data_editor(
        repeat_interruption_editor_rows(_imported_loci)
        or [
            {
                "order": 1,
                "locus_id": "",
                "allele_order": 1,
                "offset": 0,
                "bases": "",
            }
        ],
        label=(
            "Interrupciones declaradas (`allele_order` es el `order` del alelo"
            " dentro de su locus)"
        ),
        pagination=False,
    )
    repeat_min_purity_control = mo.ui.text(
        value=(
            ""
            if _imported_admission is None or _imported_admission.min_purity is None
            else _imported_admission.min_purity
        ),
        label="Pureza mínima general (vacío usa el valor de omisión)",
    )
    repeat_min_longest_run_control = mo.ui.text(
        value=(
            ""
            if _imported_admission is None
            or _imported_admission.min_longest_run_units is None
            else str(_imported_admission.min_longest_run_units)
        ),
        label="Tirada mínima general en unidades (vacío usa el valor de omisión)",
    )
    repeat_max_tract_output_control = mo.ui.text(
        value=(
            ""
            if imported_form.reference.max_tract_output_bases is None
            else str(imported_form.reference.max_tract_output_bases)
        ),
        label=(
            "Tope declarado de bases de tracto en la salida"
            " (vacío no declara tope)"
        ),
    )
    return (
        repeat_alleles_editor,
        repeat_interruptions_editor,
        repeat_loci_editor,
        repeat_max_tract_output_control,
        repeat_min_longest_run_control,
        repeat_min_purity_control,
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
    TechnicalEventForm,
    TechnicalLibraryRow,
    composition_controls,
    direct_event_from_editor_row,
    direct_events_editor,
    technical_event_mode_control,
    technical_event_id_control,
    technical_event_individual_ids_control,
    technical_event_coverage_control,
    technical_event_read_length_control,
    technical_libraries_editor,
    editor_integer,
    fasta_control,
    fasta_browser,
    reduction_editor,
    reference_id_control,
    reference_kind_control,
    render_general_toml,
    RepeatAdmissionForm,
    repeat_alleles_editor,
    repeat_interruptions_editor,
    repeat_loci_editor,
    repeat_loci_from_editor_rows,
    repeat_max_tract_output_control,
    repeat_min_longest_run_control,
    repeat_min_purity_control,
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
    reduction_from_editor_row,
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
        authored_sequence_rows = [
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
            reduction_from_editor_row(row, f"reduced_references[{index}]")
            for index, row in enumerate(ordered_reductions)
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
        _selected_technical_event = str(technical_event_mode_control.value)
        active_technical_libraries = authored_editor_rows(
            technical_libraries_editor.value, "individual_id"
        )
        technical_event = (
            None
            if _selected_technical_event == "none"
            else TechnicalEventForm(
                event_id=str(technical_event_id_control.value).strip(),
                individual_ids=[
                    individual_id.strip()
                    for individual_id in str(technical_event_individual_ids_control.value).split(",")
                    if individual_id.strip()
                ],
                coverage_wgs=str(technical_event_coverage_control.value).strip(),
                libraries=[
                    TechnicalLibraryRow(
                        individual_id=str(row["individual_id"]).strip(),
                        i7=str(row["i7"]).strip(),
                        i5=str(row["i5"]).strip(),
                    )
                    for row in ordered_editor_rows(
                        active_technical_libraries, "technical_event.libraries"
                    )
                ],
                read_length_bases=(
                    None
                    if not str(technical_event_read_length_control.value).strip()
                    else editor_integer(
                        technical_event_read_length_control.value,
                        "technical_event.read_length_bases",
                    )
                ),
            )
        )
        authored_repeat_loci = repeat_loci_from_editor_rows(
            repeat_loci_editor.value,
            repeat_alleles_editor.value,
            repeat_interruptions_editor.value,
        )
        _declared_min_purity = str(repeat_min_purity_control.value).strip()
        _declared_min_run = str(repeat_min_longest_run_control.value).strip()
        authored_repeat_admission = (
            None
            if not (_declared_min_purity or _declared_min_run)
            else RepeatAdmissionForm(
                min_purity=_declared_min_purity or None,
                min_longest_run_units=(
                    editor_integer(
                        _declared_min_run,
                        "reference.repeat_admission.min_longest_run_units",
                    )
                    if _declared_min_run
                    else None
                ),
            )
        )
        _declared_tract_cap = str(repeat_max_tract_output_control.value).strip()
        authored_max_tract_output_bases = (
            editor_integer(
                _declared_tract_cap, "reference.max_tract_output_bases"
            )
            if _declared_tract_cap
            else None
        )
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
                fasta=(
                    str(fasta_browser.path(0))
                    if fasta_browser.value
                    else fasta_control.value.strip()
                )
                if selected_kind == "supplied"
                else None,
                composition=(
                    {base: composition_controls[base].value for base in (*BASES, "N")}
                    if selected_kind == "simulated"
                    else None
                ),
                sequences=authored_sequence_rows,
                repeat_loci=authored_repeat_loci,
                repeat_admission=authored_repeat_admission,
                max_tract_output_bases=authored_max_tract_output_bases,
            ),
            reduced_references=authored_reductions,
            population_snvs=population_snvs,
            initial_populations=authored_initial_populations,
            direct_events=authored_direct_events,
            technical_event=technical_event,
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
    output_directory_control = mo.ui.text(
        value=str(PROJECT_ROOT / "results" / "marimo-run"),
        label="Directorio de resultados (debe no existir)",
        full_width=True,
    )
    output_parent_directory_browser = mo.ui.file_browser(
        initial_path=PROJECT_ROOT,
        selection_mode="directory",
        multiple=False,
        label="O seleccione la carpeta padre local de resultados",
    )
    output_name_control = mo.ui.text(
        value="marimo-run",
        label="Nombre nuevo del resultado al usar la carpeta padre",
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
        binary_path_control,
        output_name_control,
        output_directory_control,
        output_parent_directory_browser,
        validate_button,
    )


@app.cell
def _(
    ProcessAdapterError,
    allow_empty_snv_matrix_control,
    allow_incompatible_genotypes_control,
    authoring_toml,
    base_directory_control,
    base_directory_browser,
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
                    (
                        _Path(base_directory_browser.path(0))
                        if base_directory_browser.value
                        else _Path(base_directory_control.value.strip())
                    ),
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
    base_directory_browser,
    binary_path_control,
    direct_events_editor,
    effective_incompatible_genotype_authorization,
    mo,
    output_directory_control,
    output_name_control,
    output_parent_directory_browser,
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
                    (
                        _Path(base_directory_browser.path(0))
                        if base_directory_browser.value
                        else _Path(base_directory_control.value.strip())
                    ),
                    (
                        _Path(output_parent_directory_browser.path(0))
                        / output_name_control.value.strip()
                        if output_parent_directory_browser.value
                        else _Path(output_directory_control.value.strip())
                    ),
                    allow_empty_snv_matrix=allow_empty_snv_matrix_control.value,
                    allow_incompatible_genotypes=_allow_incompatible_genotypes,
                )
        except ProcessAdapterError as error:
            run_error = str(error)
    return run_error, run_result


@app.cell
def _(mo):
    # F5.6.1. La importación F5 es de **pareja**: exactamente un general v3 y
    # exactamente un escenario v1. Cero o varios son errores de autoría.
    f5_pair_upload = mo.ui.file(
        filetypes=[".toml"],
        multiple=True,
        max_size=2_000_000,
        label="Importar la pareja F5: general v3 y escenario v1",
    )
    return (f5_pair_upload,)


@app.cell
def _(EditorError, default_f5_pair, f5_pair_upload, import_f5_pair):
    f5_import_error = None
    imported_f5_pair = default_f5_pair()
    if f5_pair_upload.value:
        try:
            imported_f5_pair = import_f5_pair(
                [
                    (uploaded.name, uploaded.contents.decode("utf-8"))
                    for uploaded in f5_pair_upload.value
                ]
            )
        except (UnicodeDecodeError, EditorError) as error:
            f5_import_error = str(error)
    return f5_import_error, imported_f5_pair


@app.cell
def _(
    F5_SAMPLING_MODES,
    F5_TARGET_KINDS,
    f5_override_editor_rows,
    f5_subpopulation_editor_rows,
    imported_f5_pair,
    mo,
):
    _general = imported_f5_pair.general
    _scenario = imported_f5_pair.scenario
    f5_seed_control = mo.ui.text(
        value=_general.seed or "",
        label="Semilla raíz F5 (64 hexadecimales; vacío la resuelve Rust)",
        full_width=True,
    )
    f5_ploidy_control = mo.ui.number(
        start=1, stop=10, step=1, value=_general.ploidy, label="Ploidía F5"
    )
    f5_reference_id_control = mo.ui.text(
        value=_general.reference_id, label="ID de la referencia"
    )
    f5_reference_path_control = mo.ui.text(
        value=_general.reference_path,
        label="Ruta relativa del FASTA de referencia",
        full_width=True,
    )
    f5_reference_sha256_control = mo.ui.text(
        value=_general.reference_sha256,
        label="SHA-256 de la referencia",
        full_width=True,
    )
    f5_general_scenario_id_control = mo.ui.text(
        value=_general.scenario_id, label="ID del escenario enlazado por el general"
    )
    f5_general_scenario_path_control = mo.ui.text(
        value=_general.scenario_path,
        label="Ruta relativa del único escenario F5 enlazado",
        full_width=True,
    )
    f5_scenario_link_controls = mo.vstack(
        [f5_general_scenario_id_control, f5_general_scenario_path_control]
    )
    f5_scenario_id_control = mo.ui.text(
        value=_scenario.scenario_id, label="ID declarado por el escenario"
    )
    f5_founder_count_control = mo.ui.text(
        value=str(_scenario.founder_count), label="Fundadores del escenario"
    )
    f5_profile_metadata_path_control = mo.ui.text(
        value=_scenario.profile_metadata_path,
        label="Ruta relativa del metadata del perfil versionado",
        full_width=True,
    )
    f5_profile_metadata_sha256_control = mo.ui.text(
        value=_scenario.profile_metadata_sha256,
        label="SHA-256 del metadata del perfil",
        full_width=True,
    )
    f5_excluded_regions_path_control = mo.ui.text(
        value=_scenario.excluded_regions_path or "",
        label="Ruta relativa de la máscara de regiones excluidas (opcional)",
        full_width=True,
    )
    f5_excluded_regions_sha256_control = mo.ui.text(
        value=_scenario.excluded_regions_sha256 or "",
        label="SHA-256 de la máscara (obligatorio si hay máscara)",
        full_width=True,
    )
    f5_inbreeding_control = mo.ui.text(
        value=_scenario.inbreeding_f or "",
        label="F del escenario (solo ploidía 2; vacío lo omite)",
    )
    f5_max_layout_attempts_control = mo.ui.text(
        value="" if _scenario.max_layout_attempts is None else str(_scenario.max_layout_attempts),
        label="max_layout_attempts",
    )
    f5_max_haplotype_attempts_control = mo.ui.text(
        value="" if _scenario.max_haplotype_attempts is None else str(_scenario.max_haplotype_attempts),
        label="max_haplotype_attempts",
    )
    f5_max_swap_proposals_control = mo.ui.text(
        value="" if _scenario.max_swap_proposals is None else str(_scenario.max_swap_proposals),
        label="max_swap_proposals",
    )
    f5_max_haplotype_cells_control = mo.ui.text(
        value="" if _scenario.max_haplotype_cells is None else str(_scenario.max_haplotype_cells),
        label="max_haplotype_cells",
    )
    f5_max_operation_units_control = mo.ui.text(
        value="" if _scenario.max_operation_units is None else str(_scenario.max_operation_units),
        label="max_operation_units",
    )
    f5_max_working_memory_bytes_control = mo.ui.text(
        value="" if _scenario.max_working_memory_bytes is None else str(_scenario.max_working_memory_bytes),
        label="max_working_memory_bytes",
    )
    f5_max_output_bytes_control = mo.ui.text(
        value="" if _scenario.max_output_bytes is None else str(_scenario.max_output_bytes),
        label="max_output_bytes",
    )
    f5_target_kind_control = mo.ui.dropdown(
        options=list(F5_TARGET_KINDS),
        value=_scenario.target_kind,
        label="Objetivo exclusivo del escenario",
    )
    f5_target_value_control = mo.ui.text(
        value=_scenario.target_value, label="Valor del objetivo"
    )
    f5_subpopulations_editor = mo.ui.data_editor(
        f5_subpopulation_editor_rows(_scenario),
        label=(
            "Subpoblaciones: population_id, individual_count, founder_kind "
            "(founder_ids|founder_keys), founders separados por espacios, "
            f"sampling ({'|'.join(F5_SAMPLING_MODES)}) y F"
        ),
        pagination=False,
    )
    f5_sampling_overrides_editor = mo.ui.data_editor(
        f5_override_editor_rows(_scenario)
        or [
            {
                "order": 1,
                "population_id": "",
                "block_index": 0,
                "mode": "probabilistic",
                "controlled_max_abs_error": "",
            }
        ],
        label=(
            "Sustituciones por bloque×subpoblación: modo y tolerancia "
            "(fila con población vacía se omite)"
        ),
        pagination=False,
    )
    return (
        f5_excluded_regions_path_control,
        f5_excluded_regions_sha256_control,
        f5_founder_count_control,
        f5_general_scenario_id_control,
        f5_general_scenario_path_control,
        f5_inbreeding_control,
        f5_max_haplotype_attempts_control,
        f5_max_haplotype_cells_control,
        f5_max_layout_attempts_control,
        f5_max_operation_units_control,
        f5_max_output_bytes_control,
        f5_max_swap_proposals_control,
        f5_max_working_memory_bytes_control,
        f5_ploidy_control,
        f5_profile_metadata_path_control,
        f5_profile_metadata_sha256_control,
        f5_reference_id_control,
        f5_reference_path_control,
        f5_reference_sha256_control,
        f5_sampling_overrides_editor,
        f5_scenario_id_control,
        f5_scenario_link_controls,
        f5_seed_control,
        f5_subpopulations_editor,
        f5_target_kind_control,
        f5_target_value_control,
    )


@app.cell
def _(
    EditorError,
    F5GeneralForm,
    F5Pair,
    F5ScenarioForm,
    editor_integer,
    f5_excluded_regions_path_control,
    f5_excluded_regions_sha256_control,
    f5_founder_count_control,
    f5_general_scenario_id_control,
    f5_general_scenario_path_control,
    f5_inbreeding_control,
    f5_max_haplotype_attempts_control,
    f5_max_haplotype_cells_control,
    f5_max_layout_attempts_control,
    f5_max_operation_units_control,
    f5_max_output_bytes_control,
    f5_max_swap_proposals_control,
    f5_max_working_memory_bytes_control,
    f5_override_from_editor_row,
    f5_pair_documents,
    f5_ploidy_control,
    f5_profile_metadata_path_control,
    f5_profile_metadata_sha256_control,
    f5_reference_id_control,
    f5_reference_path_control,
    f5_reference_sha256_control,
    f5_sampling_overrides_editor,
    f5_scenario_id_control,
    f5_seed_control,
    f5_subpopulation_from_editor_row,
    f5_subpopulations_editor,
    f5_target_kind_control,
    f5_target_value_control,
    imported_f5_pair,
    ordered_editor_rows,
    pair_without_inbreeding,
):
    f5_preview_error = None
    f5_general_toml = ""
    f5_scenario_toml = ""

    def _f5_optional(control, name):
        text = str(control.value).strip()
        return editor_integer(text, f"f5.{name}") if text else None

    def _f5_authored(rows, path):
        # Una fila sin `population_id` es el hueco vacío del editor, no una
        # subpoblación ni una sustitución: se omite antes de ordenar.
        return ordered_editor_rows(
            [row for row in rows if str(row.get("population_id", "")).strip()], path
        )

    try:
        # Los bytes importados viajan con el par: `f5_pair_documents` los
        # reemite mientras sigan describiendo el formulario vigente y deja de
        # hacerlo en cuanto un control cambia (REQ-UI-15).
        _f5_pair = F5Pair(
            general=F5GeneralForm(
                seed=str(f5_seed_control.value).strip() or None,
                ploidy=editor_integer(f5_ploidy_control.value, "ploidy"),
                reference_id=str(f5_reference_id_control.value).strip(),
                reference_path=str(f5_reference_path_control.value).strip(),
                reference_sha256=str(f5_reference_sha256_control.value).strip(),
                scenario_id=str(f5_general_scenario_id_control.value).strip(),
                scenario_path=str(f5_general_scenario_path_control.value).strip(),
            ),
            scenario=F5ScenarioForm(
                scenario_id=str(f5_scenario_id_control.value).strip(),
                founder_count=editor_integer(
                    str(f5_founder_count_control.value).strip(), "founder_count"
                ),
                profile_metadata_path=str(
                    f5_profile_metadata_path_control.value
                ).strip(),
                profile_metadata_sha256=str(
                    f5_profile_metadata_sha256_control.value
                ).strip(),
                target_kind=str(f5_target_kind_control.value),
                target_value=str(f5_target_value_control.value).strip(),
                subpopulations=[
                    f5_subpopulation_from_editor_row(row, f"subpopulations[{index}]")
                    for index, row in enumerate(
                        _f5_authored(f5_subpopulations_editor.value, "subpopulations")
                    )
                ],
                sampling_overrides=[
                    f5_override_from_editor_row(row, f"sampling_overrides[{index}]")
                    for index, row in enumerate(
                        _f5_authored(
                            f5_sampling_overrides_editor.value, "sampling_overrides"
                        )
                    )
                ],
                excluded_regions_path=str(
                    f5_excluded_regions_path_control.value
                ).strip()
                or None,
                excluded_regions_sha256=str(
                    f5_excluded_regions_sha256_control.value
                ).strip()
                or None,
                inbreeding_f=str(f5_inbreeding_control.value).strip() or None,
                max_layout_attempts=_f5_optional(
                    f5_max_layout_attempts_control, "max_layout_attempts"
                ),
                max_haplotype_attempts=_f5_optional(
                    f5_max_haplotype_attempts_control, "max_haplotype_attempts"
                ),
                max_swap_proposals=_f5_optional(
                    f5_max_swap_proposals_control, "max_swap_proposals"
                ),
                max_haplotype_cells=_f5_optional(
                    f5_max_haplotype_cells_control, "max_haplotype_cells"
                ),
                max_operation_units=_f5_optional(
                    f5_max_operation_units_control, "max_operation_units"
                ),
                max_working_memory_bytes=_f5_optional(
                    f5_max_working_memory_bytes_control, "max_working_memory_bytes"
                ),
                max_output_bytes=_f5_optional(
                    f5_max_output_bytes_control, "max_output_bytes"
                ),
            ),
            general_source=imported_f5_pair.general_source,
            scenario_source=imported_f5_pair.scenario_source,
        )
        f5_general_toml, f5_scenario_toml = f5_pair_documents(
            pair_without_inbreeding(_f5_pair)
        )
    except EditorError as error:
        f5_preview_error = str(error)
    return f5_general_toml, f5_preview_error, f5_scenario_toml


@app.cell
def _(f5_import_error, f5_preview_error):
    f5_authoring_error = f5_import_error or f5_preview_error
    return (f5_authoring_error,)


@app.cell
def _(PROJECT_ROOT, f5_authoring_error, mo):
    f5_authoring_directory_control = mo.ui.text(
        value=str(PROJECT_ROOT),
        label="Directorio de autoría donde se escriben los dos documentos F5",
        full_width=True,
    )
    f5_authoring_directory_browser = mo.ui.file_browser(
        initial_path=PROJECT_ROOT,
        selection_mode="directory",
        multiple=False,
        label="O seleccione el directorio de autoría local",
    )
    f5_general_name_control = mo.ui.text(
        value="general.toml",
        label="Nombre del general: se publica en effective-config.json",
    )
    f5_overwrite_control = mo.ui.checkbox(
        value=False,
        label="Sobrescribir los documentos existentes del directorio de autoría",
    )
    f5_output_directory_control = mo.ui.text(
        value=str(PROJECT_ROOT / "results" / "marimo-f5-run"),
        label="Directorio del paquete F5 (debe no existir)",
        full_width=True,
    )
    f5_validate_button = mo.ui.run_button(
        label="Validar la pareja F5 con Rust",
        kind="success",
        disabled=bool(f5_authoring_error),
    )
    return (
        f5_authoring_directory_browser,
        f5_authoring_directory_control,
        f5_general_name_control,
        f5_output_directory_control,
        f5_overwrite_control,
        f5_validate_button,
    )


@app.cell
def _(
    ProcessAdapterError,
    binary_path_control,
    f5_authoring_directory_browser,
    f5_authoring_directory_control,
    f5_general_name_control,
    f5_general_toml,
    f5_overwrite_control,
    f5_scenario_toml,
    f5_validate_button,
    mo,
    validate_f5_pair,
):
    from pathlib import Path as _F5Path

    f5_validation_error = None
    f5_validation_result = None
    f5_validated_session = None
    if f5_validate_button.value:
        try:
            with mo.status.spinner(title="Validando la pareja F5 con Rust…"):
                f5_validation_result, f5_validated_session = validate_f5_pair(
                    _F5Path(binary_path_control.value),
                    f5_general_toml,
                    f5_scenario_toml,
                    (
                        _F5Path(f5_authoring_directory_browser.path(0))
                        if f5_authoring_directory_browser.value
                        else _F5Path(f5_authoring_directory_control.value.strip())
                    ),
                    general_name=str(f5_general_name_control.value).strip(),
                    overwrite=f5_overwrite_control.value,
                )
        except ProcessAdapterError as error:
            f5_validation_error = str(error)
    return f5_validation_error, f5_validation_result, f5_validated_session


@app.cell
def _(f5_authoring_error, f5_validated_session, mo):
    # REQ-UI-15: una validación vigente es obligatoria antes de ejecutar.
    f5_run_button = mo.ui.run_button(
        label="Ejecutar el escenario F5",
        kind="warn",
        disabled=bool(f5_authoring_error) or f5_validated_session is None,
    )
    return (f5_run_button,)


@app.cell
def _(
    ProcessAdapterError,
    binary_path_control,
    f5_authoring_directory_browser,
    f5_authoring_directory_control,
    f5_general_name_control,
    f5_general_toml,
    f5_output_directory_control,
    f5_run_button,
    f5_scenario_toml,
    f5_validated_session,
    mo,
    run_f5_pair,
):
    from pathlib import Path as _F5RunPath

    f5_run_error = None
    f5_run_result = None
    if f5_run_button.value and f5_validated_session is not None:
        try:
            with mo.status.spinner(title="Ejecutando el escenario F5…"):
                f5_run_result = run_f5_pair(
                    f5_validated_session,
                    f5_general_toml,
                    f5_scenario_toml,
                    _F5RunPath(binary_path_control.value),
                    (
                        _F5RunPath(f5_authoring_directory_browser.path(0))
                        if f5_authoring_directory_browser.value
                        else _F5RunPath(f5_authoring_directory_control.value.strip())
                    ),
                    _F5RunPath(f5_output_directory_control.value.strip()),
                    general_name=str(f5_general_name_control.value).strip(),
                )
        except ProcessAdapterError as error:
            f5_run_error = str(error)
    return f5_run_error, f5_run_result


@app.cell
def _(
    allow_empty_snv_matrix_control,
    allow_incompatible_genotypes_control,
    authoring_error,
    authoring_toml,
    f5_authoring_directory_browser,
    f5_authoring_directory_control,
    f5_authoring_error,
    f5_excluded_regions_path_control,
    f5_excluded_regions_sha256_control,
    f5_founder_count_control,
    f5_general_name_control,
    f5_general_toml,
    f5_inbreeding_control,
    f5_max_haplotype_attempts_control,
    f5_max_haplotype_cells_control,
    f5_max_layout_attempts_control,
    f5_max_operation_units_control,
    f5_max_output_bytes_control,
    f5_max_swap_proposals_control,
    f5_max_working_memory_bytes_control,
    f5_output_directory_control,
    f5_overwrite_control,
    f5_package_summary,
    f5_pair_upload,
    f5_ploidy_control,
    f5_profile_metadata_path_control,
    f5_profile_metadata_sha256_control,
    f5_reference_id_control,
    f5_reference_path_control,
    f5_reference_sha256_control,
    f5_run_button,
    f5_run_error,
    f5_run_result,
    f5_sampling_overrides_editor,
    f5_scenario_id_control,
    f5_scenario_link_controls,
    f5_scenario_toml,
    f5_seed_control,
    f5_subpopulations_editor,
    f5_target_kind_control,
    f5_target_value_control,
    f5_validate_button,
    f5_validation_error,
    f5_validation_result,
    f5_validated_session,
    available_sampling_modes,
    inbreeding_is_available,
    base_directory_control,
    base_directory_browser,
    binary_path_control,
    composition_controls,
    config_upload,
    load_default_config_button,
    effective_incompatible_genotype_authorization,
    fasta_control,
    fasta_browser,
    fasta_headers_error,
    import_error,
    imported_name,
    explicit_snp_authorization_summary,
    initial_population_statistics_summary,
    initial_populations_editor,
    individuals_panel,
    mo,
    repeats_panel,
    repeat_alleles_editor,
    repeat_interruptions_editor,
    repeat_loci_editor,
    repeat_max_tract_output_control,
    repeat_min_longest_run_control,
    repeat_min_purity_control,
    output_directory_control,
    output_name_control,
    output_parent_directory_browser,
    reduction_editor,
    read_fasta_headers_button,
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
    technical_event_mode_control,
    technical_event_id_control,
    technical_event_individual_ids_control,
    technical_event_coverage_control,
    technical_event_read_length_control,
    technical_libraries_editor,
):
    visible_error = authoring_error
    visible_preview = authoring_toml
    general_panel = mo.vstack(
        [
            mo.md("## General/Ejecución"),
            config_upload,
            load_default_config_button,
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
        mo.vstack(
            [
                fasta_control,
                fasta_browser,
                read_fasta_headers_button,
                mo.callout(fasta_headers_error, kind="danger", title="No se pudo leer el FASTA")
                if fasta_headers_error
                else mo.md(
                    "Al leerlo se rellenan ID y longitud en orden de FASTA; "
                    "todas las filas empiezan como `scaffold`."
                ),
            ]
        )
        if reference_kind_control.value == "supplied"
        else mo.vstack([
            mo.md("A/C/G/T deben sumar 1. `N` es una fracción adicional de posiciones ambiguas."),
            mo.hstack([composition_controls[base] for base in ("A", "C", "G", "T", "N")]),
        ])
    )
    reduction_controls = (
        mo.vstack(
            [
                mo.callout(
                    "La selección es global entre cromosomas y scaffolds elegibles y forma bloques no solapados de `region_size_bases`. "
                    "«Eliminar porcentaje» usa `reduction_percent`: valor 10 elimina el 10 % y retiene el 90 %. "
                    "Para conservar el 10 %, indique «Eliminar porcentaje» y valor 90.",
                    kind="info",
                    title="Cómo crear una referencia reducida",
                ),
                reduction_editor,
            ]
        )
        if reference_kind_control.value == "supplied"
        else mo.callout(
            "La reducción de referencia solo está disponible para una referencia FASTA suministrada.",
            kind="neutral",
        )
    )
    genome_panel = mo.vstack(
        [
            mo.md("## Genoma"),
            mo.hstack([reference_kind_control, reference_id_control]),
            selected_source_controls,
            sequence_editor,
            reduction_controls,
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
    if f5_validation_error:
        f5_validation_status = mo.callout(
            f5_validation_error, kind="danger", title="No se pudo validar la pareja F5"
        )
    elif f5_validation_result is None:
        f5_validation_status = mo.callout(
            "Pendiente de validar la pareja F5 con Rust.", kind="neutral"
        )
    elif f5_validation_result.succeeded:
        f5_validation_status = mo.callout(
            "Pareja F5 admitida: el preflight de recursos, el perfil versionado "
            "y la referencia ya se resolvieron.",
            kind="success",
        )
    else:
        f5_validation_status = mo.callout(
            f5_validation_result.stdout or f5_validation_result.stderr,
            kind="danger",
            title=f"Validación F5 rechazada (código {f5_validation_result.returncode})",
        )
    if f5_run_error:
        f5_run_status = mo.callout(
            f5_run_error, kind="danger", title="No se pudo ejecutar el escenario F5"
        )
    elif f5_run_result is None:
        f5_run_status = mo.callout("Ejecución F5 no iniciada.", kind="neutral")
    elif f5_run_result.succeeded:
        f5_run_status = mo.callout(
            f"Paquete F5 publicado en {f5_run_result.output_directory}.", kind="success"
        )
    else:
        f5_run_status = mo.callout(
            f5_run_result.stdout or f5_run_result.stderr,
            kind="danger",
            title=f"Ejecución F5 fallida (código {f5_run_result.returncode})",
        )
    f5_package_summary_text = (
        f5_package_summary(f5_run_result.output_directory)
        if f5_run_result is not None
        and f5_run_result.succeeded
        and f5_run_result.output_directory is not None
        else "El resumen del paquete F5 aparecerá tras una ejecución correcta."
    )
    # F5.6.1. El escenario poblacional sintético: un general v3 y un único
    # escenario F5 v1. Sin panel externo, haplotipos suministrados,
    # construcción recombinante ni mapas (REQ-UI-15).
    f5_status = (
        mo.callout(f5_authoring_error, kind="danger", title="Error de autoría F5")
        if f5_authoring_error
        else mo.callout(
            "Pareja F5 estructuralmente completa. Valide con Rust antes de ejecutar.",
            kind="success",
        )
    )
    f5_inbreeding_controls = (
        mo.vstack([f5_inbreeding_control])
        if inbreeding_is_available(int(f5_ploidy_control.value))
        else mo.callout(
            "F solo se edita con ploidía 2; con otra ploidía no se publica.",
            kind="neutral",
        )
    )
    f5_offered_sampling_modes = available_sampling_modes(
        str(f5_inbreeding_control.value).strip() or None
    )
    f5_budget_controls = mo.vstack(
        [
            f5_max_layout_attempts_control,
            f5_max_haplotype_attempts_control,
            f5_max_swap_proposals_control,
            f5_max_haplotype_cells_control,
            f5_max_operation_units_control,
            f5_max_working_memory_bytes_control,
            f5_max_output_bytes_control,
        ]
    )
    f5_panel = mo.vstack(
        [
            mo.md("### Escenario poblacional sintético F5"),
            f5_pair_upload,
            f5_status,
            mo.md("#### General v3"),
            f5_seed_control,
            f5_ploidy_control,
            f5_reference_id_control,
            f5_reference_path_control,
            f5_reference_sha256_control,
            f5_scenario_link_controls,
            mo.md("#### Escenario F5 v1"),
            f5_scenario_id_control,
            f5_founder_count_control,
            f5_profile_metadata_path_control,
            f5_profile_metadata_sha256_control,
            f5_excluded_regions_path_control,
            f5_excluded_regions_sha256_control,
            f5_inbreeding_controls,
            mo.md("#### Objetivo exclusivo, límites y presupuestos"),
            f5_target_kind_control,
            f5_target_value_control,
            f5_budget_controls,
            mo.md("#### Subpoblaciones y sustituciones por bloque×subpoblación"),
            f5_subpopulations_editor,
            f5_sampling_overrides_editor,
            mo.callout(
                "Modos de muestreo ofrecidos con el F declarado: "
                + ", ".join(f5_offered_sampling_modes)
                + (
                    "."
                    if len(f5_offered_sampling_modes) > 1
                    else "; `controlled` se desactiva porque F es mayor que cero."
                ),
                kind="info",
            ),
            mo.ui.text_area(
                value=f5_general_toml,
                label="General v3 exportable",
                rows=12,
                disabled=True,
                full_width=True,
            ),
            mo.ui.text_area(
                value=f5_scenario_toml,
                label="Escenario F5 v1 exportable",
                rows=16,
                disabled=True,
                full_width=True,
            ),
            mo.hstack(
                [
                    mo.download(
                        f5_general_toml,
                        filename=str(f5_general_name_control.value).strip()
                        or "general.toml",
                        mimetype="application/toml",
                        disabled=bool(f5_authoring_error),
                        label="Exportar el general v3",
                    ),
                    mo.download(
                        f5_scenario_toml,
                        filename="scenario.toml",
                        mimetype="application/toml",
                        disabled=bool(f5_authoring_error),
                        label="Exportar el escenario F5 v1",
                    ),
                ]
            ),
            mo.md("#### Preflight y ejecución F5"),
            f5_authoring_directory_control,
            f5_authoring_directory_browser,
            f5_general_name_control,
            f5_overwrite_control,
            f5_output_directory_control,
            mo.hstack([f5_validate_button, f5_run_button]),
            f5_validation_status,
            f5_run_status,
            mo.ui.text_area(
                value=f5_package_summary_text,
                label="Resumen del paquete F5 publicado (solo desde el manifiesto)",
                rows=10,
                disabled=True,
                full_width=True,
            ),
        ]
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
            f5_panel,
        ]
    )
    vcf_mode = str(vcf_output_mode_control.value)
    outputs_panel = mo.vstack(
        [
            mo.md("## Salidas"),
            technical_event_mode_control,
            (
                mo.vstack([
                    technical_event_id_control,
                    technical_event_individual_ids_control,
                    technical_event_coverage_control,
                    technical_event_read_length_control,
                    technical_libraries_editor,
                ])
                if str(technical_event_mode_control.value) == "configured"
                else mo.callout("No se solicitará evento técnico.", kind="neutral")
            ),
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
            base_directory_browser,
            output_directory_control,
            output_parent_directory_browser,
            output_name_control,
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
            "Repeticiones": repeats_panel(
                mo,
                repeat_loci_editor,
                repeat_alleles_editor,
                repeat_interruptions_editor,
                mo.vstack(
                    [
                        repeat_min_purity_control,
                        repeat_min_longest_run_control,
                        repeat_max_tract_output_control,
                    ]
                ),
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
