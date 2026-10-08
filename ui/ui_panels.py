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


# Los dos valores del selector de procedencia de los alelos. Viven aquí, junto
# al panel que los usa, para que la interfaz y sus pruebas no puedan discrepar.
REPEAT_ALLELE_SOURCE_DECLARED = "Declarar los alelos"
REPEAT_ALLELE_SOURCE_GENERATED = "Generarlos por deslizamiento"
REPEAT_ALLELE_SOURCES = (
    REPEAT_ALLELE_SOURCE_DECLARED,
    REPEAT_ALLELE_SOURCE_GENERATED,
)


def repeats_panel(
    mo,
    repeat_loci_editor,
    repeat_allele_source_control,
    repeat_alleles_editor,
    repeat_interruptions_editor,
    repeat_slippage_editor,
    repeat_admission_controls,
):
    """Render the tandem repeat authoring controls; Rust stays the authority.

    **El selector de procedencia decide qué se ve.** Con «declarar» se muestran
    la tabla de alelos y la de interrupciones y **no** la del generador; con
    «generar» se muestra la del generador y **no** las otras dos. No es un
    adorno: declarar los alelos y pedir además que se generen es excluyente, y
    una interfaz que enseñara las dos cosas a la vez invitaría a escribir un
    documento que el lector rechaza.

    Lo que **no** hay aquí es tan normativo como lo que hay: no existe control
    de pureza objetivo, porque la pureza no se declara nunca —se calcula sobre
    la secuencia—, ni control para generar loci al azar, que pertenece a un
    paso posterior. Una casilla que los ofreciera prometería algo que el
    simulador rechaza.
    """
    generating = getattr(repeat_allele_source_control, "value", None) == (
        REPEAT_ALLELE_SOURCE_GENERATED
    )
    children = [
        mo.md("## Repeticiones en tándem"),
        mo.callout(
            "Declare cada locus con su motivo y su posición. Después elija si "
            "sus alelos los escribe usted o los genera el simulador partiendo "
            "del alelo que la referencia lleva en ese locus. La pureza no se "
            "declara nunca: se calcula sobre la secuencia de cada alelo.",
            kind="info",
        ),
        mo.md("### Loci"),
        repeat_loci_editor,
        mo.md("### Procedencia de los alelos"),
        repeat_allele_source_control,
    ]
    if generating:
        children.extend(
            [
                mo.md("### Generador por deslizamiento"),
                repeat_slippage_editor,
                mo.callout(
                    "Cada copia haplotípica deriva por su cuenta desde el alelo "
                    "de la referencia, con un número de pasos sorteado de una "
                    "Poisson de media `mean_events`. La frecuencia de cada alelo "
                    "es la proporción de copias que terminan en él. El alelo de "
                    "la referencia entra si alguna copia termina en él, y su "
                    "presencia no se fuerza.",
                    kind="info",
                ),
            ]
        )
    else:
        children.extend(
            [
                mo.md("### Alelos por locus"),
                repeat_alleles_editor,
                mo.callout(
                    "Cada alelo declara `units` —un valor por bloque, separados "
                    "por comas— o `sequence` con sus bases exactas; son "
                    "excluyentes. Las frecuencias de un locus deben sumar uno "
                    "dentro de la tolerancia de normalización.",
                    kind="info",
                ),
                mo.md("### Interrupciones declaradas"),
                repeat_interruptions_editor,
                mo.callout(
                    "Una interrupción sustituye las bases del alelo desde "
                    "`offset`. Solo se admite en alelos declarados por `units`: "
                    "uno declarado por su secuencia ya las lleva dentro.",
                    kind="info",
                ),
            ]
        )
    children.extend([mo.md("### Mínimos generales de admisión"), repeat_admission_controls])
    return mo.vstack(children)


def pcr_panel(
    mo,
    library_preparation_control,
    pcr_duplicate_fraction_control,
    pcr_derived_mean_pairs,
    pcr_cycles_control,
    pcr_efficiency_control,
    pcr_early_cycles_control,
    pcr_max_family_size_control,
    pcr_gc_max_attempts_control,
    pcr_gc_curve_editor,
    pcr_efficiency_classes_editor,
    pcr_per_cycle_slip_editor,
    pcr_homopolymer_fold_control,
):
    """Render the amplified-preparation controls; Rust stays the authority.

    Lo importante de este panel es lo que **advierte**: la fracción que el
    usuario declara es la esperada de una biblioteca homogénea, y la realizada
    será mayor. Sin ese aviso, quien compare las dos no las verá cuadrar nunca
    y creerá que el simulador se equivoca.

    La `m` derivada se muestra junto a la fracción declarada porque es el
    parámetro que de verdad entra en el modelo: enseñarla evita que el usuario
    pida una fracción sin ver qué está pidiendo. **Es orientativa**: la calcula
    la interfaz con el mismo algoritmo pero con la exponencial del sistema, y
    la que gobierna la ejecución es la que publica la configuración efectiva.
    """
    return mo.vstack(
        [
            mo.md("## Preparación de biblioteca"),
            library_preparation_control,
            mo.callout(
                "Sin PCR, cada fragmento da una pareja. Con PCR, cada molécula "
                "da una FAMILIA de parejas: de ahí salen los duplicados, la "
                "sobrerrepresentación de unas moléculas frente a otras y el "
                "sesgo por composición de bases.",
                kind="info",
            ),
            mo.md("### Fracción de duplicados"),
            pcr_duplicate_fraction_control,
            pcr_derived_mean_pairs,
            mo.callout(
                "Es la fracción ESPERADA de una biblioteca homogénea, no la "
                "realizada. La realizada es MAYOR por tres motivos que se suman: "
                "las clases de eficiencia, la curva de composición y la "
                "desigualdad de los primeros ciclos, que actúa ya por omisión. "
                "El manifiesto publica la realizada aparte. Y tampoco es la "
                "fracción que marca una herramienta de análisis contando parejas "
                "que caen en el mismo sitio: esa incluye además las copias de la "
                "química de agrupación, que no son PCR.",
                kind="warn",
            ),
            mo.callout(
                "No hay valor por omisión, y no es un descuido: depende de la "
                "masa de ADN de entrada, del número de ciclos y del laboratorio. "
                "La mejor fuente son los informes de calidad del propio "
                "laboratorio, no una media de genomas ajenos.",
                kind="info",
            ),
            mo.md("### Ciclos y eficiencia"),
            pcr_cycles_control,
            pcr_efficiency_control,
            pcr_early_cycles_control,
            mo.callout(
                "Los ciclos explícitos son los que se simulan uno a uno para "
                "capturar la desigualdad de los primeros; los posteriores entran "
                "como un factor determinista. Con pocos ciclos explícitos ya se "
                "alcanza casi toda la desigualdad.",
                kind="info",
            ),
            mo.md("### Topes de la familia y de la composición"),
            pcr_max_family_size_control,
            pcr_gc_max_attempts_control,
            mo.callout(
                "Agotado el tope de intentos se acepta el último y se CUENTA: el "
                "manifiesto publica cuántas veces ocurrió, para que el usuario "
                "juzgue si el sesgo le afecta.",
                kind="info",
            ),
            mo.md("### Curva de composición"),
            pcr_gc_curve_editor,
            mo.callout(
                "Eficiencia relativa por fracción de GC, en tramos lineales. La "
                "curva neutra —eficiencia 1 en todo el intervalo— es el valor de "
                "omisión y no consume aleatoriedad. Ninguna curva calibrada se "
                "incorpora.",
                kind="info",
            ),
            mo.md("### Clases de eficiencia"),
            pcr_efficiency_classes_editor,
            mo.callout(
                "Opcionales y sin valores por omisión. Si se declaran, sus pesos "
                "tienen que sumar 1 dentro de la tolerancia de normalización.",
                kind="info",
            ),
            mo.callout(
                "OJO: si declara clases, la eficiencia por ciclo de arriba DEJA DE "
                "INTERVENIR. Cada familia usa la de su clase. La escalar sigue "
                "siendo obligatoria y se publica, pero ya no cambia nada.",
                kind="warn",
            ),
            mo.md("### Tartamudeo clonal heredable"),
            pcr_homopolymer_fold_control,
            mo.callout(
                "Multiplicador de la PCR en los homopolímeros: cuántas veces "
                "multiplica la PCR la fracción de lecturas con tartamudeo en un "
                "motivo de 1 pb. En blanco rige Gymrek (2016), 0,17 frente a "
                "0,027; Fungtammasan y col. (2015) dan 2.5. Escribir 6.30 NO "
                "equivale a dejarlo en blanco. Se aplica sobre la constante sin "
                "PCR de la calibración, 0,027, no sobre la tasa de base de la "
                "pestaña de tartamudeo. Entre 1 y 37, y no junto a una "
                "sustitución de la longitud de motivo 1.",
                kind="info",
            ),
            pcr_per_cycle_slip_editor,
            mo.callout(
                "Con PCR, el deslizamiento puede ocurrir en CUALQUIER ciclo, y "
                "una copia que resbala transmite su número de unidades a toda su "
                "descendencia. De ahí sale el artefacto que importa: una familia "
                "de duplicados puede mostrar dos o tres alelos en proporciones "
                "que imitan a un heterocigoto, procediendo de una sola molécula.",
                kind="info",
            ),
            mo.callout(
                "Lo que se declara aquí es la probabilidad de UN CICLO, no la "
                "fracción de extremo a extremo: son dos magnitudes distintas. Con "
                "diez ciclos, la del homopolímero vale 0,0158 y su fracción de "
                "extremo a extremo 0,17. En blanco se despeja de las dos columnas "
                "de la calibración, que es lo normal.",
                kind="warn",
            ),
            mo.callout(
                "Las tasas por ciclo derivadas son TOSCAS, y el manifiesto las "
                "publica marcadas como tales: las dos columnas de la calibración "
                "proceden de dos individuos distintos, de proyectos distintos y a "
                "coberturas distintas, de modo que la resta puede llevar dentro "
                "efectos de individuo y de lote además del de la amplificación.",
                kind="warn",
            ),
            mo.callout(
                "Solo tiene efecto con loci de repetición declarados. Y con ellos, "
                "la preparación amplificada EXIGE el tartamudeo de base: sin él, "
                "con la calibración por omisión, la biblioteca subestimaría el "
                "tartamudeo entre un 14 % y un 28 % según la longitud del motivo.",
                kind="info",
            ),
        ]
    )


def stutter_panel(
    mo,
    stutter_mode_control,
    stutter_run_length_coefficient_control,
    stutter_motif_rates_editor,
    stutter_steps_editor,
    stutter_locus_rates_editor,
):
    """Render the base stutter controls; Rust stays the authority.

    Lo que **no** hay aquí es tan normativo como lo que hay: no existe control
    para un paso de dos unidades ni para ningún modelo distinto del admitido,
    porque la estructura del paso la fija el algoritmo y declarar otra cosa se
    rechaza. Tampoco hay reparto sugerido de pesos: son obligatorios y no
    tienen valor por omisión, y rellenarlos aquí afirmaría una calibración que
    nadie ha hecho.
    """
    return mo.vstack(
        [
            mo.md("## Tartamudeo de base"),
            mo.callout(
                "Una molécula que contiene un microsatélite puede mostrar un "
                "número de repeticiones distinto del que el individuo tiene, "
                "porque la polimerasa resbala al copiarla. Es un artefacto de la "
                "lectura: la verdad del individuo no cambia. No es una "
                "preparación amplificada, no introduce duplicados y el perfil se "
                "sigue registrando como sin PCR.",
                kind="info",
            ),
            stutter_mode_control,
            mo.md("### Tasa por longitud de motivo"),
            stutter_motif_rates_editor,
            mo.callout(
                "Las cifras que vienen puestas proceden de la columna sin PCR de "
                "la calibración publicada, para humano e Illumina. Solo se cubren "
                "longitudes de 1 a 6 pares de base: un locus con un motivo más "
                "largo se rechaza en vez de extrapolarse.",
                kind="info",
            ),
            stutter_run_length_coefficient_control,
            mo.callout(
                "El coeficiente multiplica la tirada ininterrumpida más larga, "
                "que es el predictor del resbalón. Arranca en cero y se publica "
                "como no calibrado: la curva que lo fijaría está publicada pero "
                "sin digitalizar.",
                kind="warn",
            ),
            mo.md("### Pesos de los cuatro pasos admitidos"),
            stutter_steps_editor,
            mo.callout(
                "Contracciones de hasta tres unidades y expansión de una sola. "
                "La estructura no se configura: solo sus pesos relativos, que son "
                "obligatorios con el tartamudeo activo.",
                kind="info",
            ),
            mo.md("### Sustitución de la tasa por locus"),
            stutter_locus_rates_editor,
            mo.callout(
                "Una sustitución tiene que nombrar un locus declarado en la "
                "referencia; un identificador que no exista se rechaza.",
                kind="info",
            ),
        ]
    )
