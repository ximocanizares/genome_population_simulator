"""F5.6.1. Controles F5 en Marimo y paridad con el TOML manual.

El término independiente de la paridad es un TOML **escrito a mano** en este
fichero y ejecutado por `subprocess` sin pasar por el adaptador. Comparar la
salida de la interfaz con la de la interfaz no demostraría nada.
"""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
import sys
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "ui"))

from config_editor import EditorError  # noqa: E402
from f5_editor import (  # noqa: E402
    F5_BUDGET_FIELDS,
    F5_GENERAL_FIELD_CONTROLS,
    F5_INTEGER_BOUNDS,
    F5_ID_PATTERN,
    F5_RELATIVE_PATH_PATTERN,
    F5_SCENARIO_FIELD_CONTROLS,
    F5_SHA_PATTERN,
    F5SamplingOverrideRow,
    F5SubpopulationRow,
    available_sampling_modes,
    controlled_is_available,
    default_f5_pair,
    f5_control_names,
    f5_override_editor_rows,
    f5_override_from_editor_row,
    f5_pair_documents,
    f5_subpopulation_editor_rows,
    f5_subpopulation_from_editor_row,
    import_f5_pair,
    inbreeding_is_available,
    pair_without_inbreeding,
    parse_f5_general_toml,
    parse_f5_scenario_toml,
    render_f5_general_toml,
    render_f5_scenario_toml,
)
from process_adapter import (  # noqa: E402
    SUPPORTED_BINARY_NAME,
    ProcessAdapterError,
    f5_package_summary,
    run_f5_pair,
    validate_f5_pair,
)

SEED = "a1a2a3a4a5a6a7a8b1b2b3b4b5b6b7b8c1c2c3c4c5c6c7c8d1d2d3d4d5d6d7d8"

PROFILE_TSV = (
    ("block_lengths.tsv", "f5_block_length_profile", "length_min_bp\tlength_max_bp\tweight\n3\t3\t1\n"),
    ("snp_gaps.tsv", "f5_snp_gap_profile", "gap_min_bp\tgap_max_bp\tweight\n1\t1\t1\n"),
    ("snp_spectrum.tsv", "f5_snp_spectrum", "alt_af_min\talt_af_max\tweight\n0.25\t0.75\t1\n"),
    (
        "haplotype_weights.tsv",
        "f5_haplotype_weights",
        "haplotype_rank\tweight\n0\t0.5\n1\t0.5\n",
    ),
    (
        "ld_profile.tsv",
        "f5_ld_profile",
        "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2"
        "\tmax_abs_error\nP1\t1\t1\t1\t0.05\nP1\t2\t2\t1\t0.05\n",
    ),
)

REFERENCE_FASTA = ">chr1\nACGTAC\n"


def sha256(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def write_supporting_files(root: Path) -> tuple[str, str]:
    """Escribe referencia y perfil versionado; devuelve sus dos SHA-256."""
    profiles = root / "scenarios" / "profiles"
    profiles.mkdir(parents=True, exist_ok=True)
    (root / "references").mkdir(parents=True, exist_ok=True)
    (root / "references" / "ref.fa").write_text(REFERENCE_FASTA, encoding="utf-8")
    links = []
    for name, kind, contents in PROFILE_TSV:
        (profiles / name).write_text(contents, encoding="utf-8")
        links.append(
            '{"artifact_kind":"%s","schema_version":"1","path":"%s",'
            '"size_bytes":%d,"sha256":"%s"}'
            % (kind, name, len(contents.encode("utf-8")), sha256(contents))
        )
    metadata = (
        '{"profile_kind":"f5_generation_profile","schema_version":"1",'
        '"profile_id":"profile","provenance_kind":"user_supplied",'
        '"source_citation":"test","source_license":null,"source_sha256":null,'
        '"files":[%s]}' % ",".join(links)
    )
    (profiles / "profile.metadata.json").write_text(metadata, encoding="utf-8")
    return sha256(REFERENCE_FASTA), sha256(metadata)


def manual_general(reference_sha: str) -> str:
    """General v3 escrito a mano. Es el término independiente de la paridad."""
    return (
        'config_kind = "general"\n'
        'schema_version = "3"\n'
        f'seed = "{SEED}"\n'
        "ploidy = 2\n"
        'reference_id = "ref"\n'
        'reference_path = "references/ref.fa"\n'
        f'reference_sha256 = "{reference_sha}"\n'
        "\n"
        "[f5_scenario]\n"
        'scenario_id = "scenario"\n'
        'scenario_path = "scenarios/scenario.toml"\n'
    )


def manual_scenario(profile_sha: str) -> str:
    """Escenario F5 v1 escrito a mano, con presupuestos y una sustitución."""
    return (
        'config_kind = "f5_scenario"\n'
        'schema_version = "1"\n'
        'scenario_id = "scenario"\n'
        "founder_count = 1\n"
        "\n"
        "[f5]\n"
        'profile_metadata_path = "profiles/profile.metadata.json"\n'
        f'profile_metadata_sha256 = "{profile_sha}"\n'
        "max_layout_attempts = 8\n"
        "max_haplotype_attempts = 32\n"
        "max_swap_proposals = 256\n"
        "max_haplotype_cells = 65536\n"
        "block_count = 2\n"
        "\n"
        "[[subpopulations]]\n"
        'population_id = "P1"\n'
        "individual_count = 1\n"
        'founder_ids = ["I1"]\n'
        "\n"
        "[[sampling_overrides]]\n"
        'population_id = "P1"\n'
        "block_index = 1\n"
        'mode = "controlled"\n'
        "controlled_max_abs_error = 0.05\n"
    )


def tree_digests(package: Path) -> dict[str, str]:
    """SHA-256 de cada fichero publicado, calculado aquí y no leído del manifiesto."""
    return {
        str(path.relative_to(package)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(package.rglob("*"))
        if path.is_file()
    }


class F5FieldControlMatrixTests(unittest.TestCase):
    """ACE-114: toda opción implementada tiene control, y ninguna futura lo tiene."""

    def test_the_matrix_covers_both_schemas_and_every_control_exists(self) -> None:
        source = (ROOT / "ui" / "gps_app.py").read_text(encoding="utf-8")
        for schema_name, matrix in (
            ("config-general-v3.schema.json", F5_GENERAL_FIELD_CONTROLS),
            ("config-f5-scenario-v1.schema.json", F5_SCENARIO_FIELD_CONTROLS),
        ):
            schema = json.loads((ROOT / "schemas" / schema_name).read_text())
            declared = set(schema["properties"]) - {"config_kind", "schema_version"}
            if schema_name.endswith("scenario-v1.schema.json"):
                # `f5` es un objeto: sus propiedades son las que se editan.
                declared -= {"f5"}
                declared |= set(schema["$defs"]["f5"]["properties"])
            with self.subTest(schema=schema_name):
                self.assertEqual(declared, set(matrix))
        for control in f5_control_names():
            with self.subTest(control=control):
                self.assertIn(control, source)

    def test_no_panel_supplied_haplotypes_recombination_or_map_control_exists(self) -> None:
        source = (ROOT / "ui" / "gps_app.py").read_text(encoding="utf-8")
        editor = (ROOT / "ui" / "f5_editor.py").read_text(encoding="utf-8")
        for forbidden in (
            "panel_control",
            "panel_path",
            "supplied_haplotypes",
            "haplotype_panel",
            "recombination",
            "genetic_map",
            "recombination_probability_per_interbase",
        ):
            with self.subTest(forbidden=forbidden):
                self.assertNotIn(forbidden, source)
                self.assertNotIn(forbidden, editor)
        # Control positivo: la búsqueda sí encuentra lo que F5.6.1 sí añade.
        self.assertIn("f5_subpopulations_editor", source)
        self.assertIn("f5_sampling_overrides_editor", source)

    def test_the_integer_bounds_are_the_ones_the_two_schemas_declare(self) -> None:
        general = json.loads((ROOT / "schemas" / "config-general-v3.schema.json").read_text())
        scenario = json.loads(
            (ROOT / "schemas" / "config-f5-scenario-v1.schema.json").read_text()
        )
        nodes = {
            "ploidy": general["properties"]["ploidy"],
            "founder_count": scenario["properties"]["founder_count"],
            "individual_count": scenario["$defs"]["subpopulation"]["properties"][
                "individual_count"
            ],
            "block_index": scenario["$defs"]["override"]["properties"]["block_index"],
        }
        for name in (*F5_BUDGET_FIELDS, "block_count", "covered_bases"):
            nodes[name] = scenario["$defs"]["f5"]["properties"][name]
        self.assertEqual(set(nodes), set(F5_INTEGER_BOUNDS))
        for name, node in nodes.items():
            with self.subTest(field=name):
                self.assertEqual(
                    (node["minimum"], node["maximum"]), F5_INTEGER_BOUNDS[name]
                )


    def test_the_three_patterns_are_the_ones_the_two_schemas_declare(self) -> None:
        general = json.loads((ROOT / "schemas" / "config-general-v3.schema.json").read_text())
        scenario = json.loads(
            (ROOT / "schemas" / "config-f5-scenario-v1.schema.json").read_text()
        )
        declared = {
            "id": general["$defs"]["id"]["pattern"],
            "sha": general["$defs"]["sha"]["pattern"],
            "relative_path": general["$defs"]["relative_path"]["pattern"],
        }
        # Los dos esquemas declaran las mismas tres, y el editor una sola copia.
        self.assertEqual(declared["id"], scenario["$defs"]["id"]["pattern"])
        self.assertEqual(declared["sha"], scenario["$defs"]["sha"]["pattern"])
        self.assertEqual(declared["relative_path"], scenario["$defs"]["path"]["pattern"])
        for name, compiled in (
            ("id", F5_ID_PATTERN),
            ("sha", F5_SHA_PATTERN),
            ("relative_path", F5_RELATIVE_PATH_PATTERN),
        ):
            with self.subTest(pattern=name):
                # El editor ancla con `fullmatch`, de modo que su expresión es
                # la del esquema sin los anclajes `^`/`$`.
                self.assertEqual(declared[name], f"^{compiled.pattern}$")

    def test_the_three_patterns_accept_and_reject_the_same_strings(self) -> None:
        cases = (
            (F5_ID_PATTERN, "P1", True),
            (F5_ID_PATTERN, "-P1", False),
            (F5_ID_PATTERN, "A" * 65, False),
            (F5_SHA_PATTERN, "a" * 64, True),
            (F5_SHA_PATTERN, "A" * 64, False),
            (F5_RELATIVE_PATH_PATTERN, "scenarios/scenario.toml", True),
            (F5_RELATIVE_PATH_PATTERN, "/scenarios/scenario.toml", False),
            (F5_RELATIVE_PATH_PATTERN, "../scenario.toml", False),
            (F5_RELATIVE_PATH_PATTERN, "a\\b", False),
        )
        for pattern, value, accepted in cases:
            with self.subTest(value=value):
                self.assertEqual(pattern.fullmatch(value) is not None, accepted)


class F5ImportExportTests(unittest.TestCase):
    """REQ-UI-15: exactamente un general v3 y exactamente un escenario v1."""

    def setUp(self) -> None:
        self.general = manual_general("a" * 64)
        self.scenario = manual_scenario("b" * 64)

    def test_the_pair_round_trips_byte_for_byte(self) -> None:
        pair = import_f5_pair(
            [("general.toml", self.general), ("scenario.toml", self.scenario)]
        )
        exported_general, exported_scenario = f5_pair_documents(pair)
        self.assertEqual(exported_general, self.general)
        self.assertEqual(exported_scenario, self.scenario)

    def test_the_order_of_the_two_documents_does_not_matter(self) -> None:
        forward = import_f5_pair(
            [("general.toml", self.general), ("scenario.toml", self.scenario)]
        )
        backward = import_f5_pair(
            [("scenario.toml", self.scenario), ("general.toml", self.general)]
        )
        self.assertEqual(forward, backward)

    def test_two_legal_key_orders_survive_the_round_trip_byte_for_byte(self) -> None:
        """REQ-UI-15 pide conservar **el mismo** TOML, no uno equivalente.

        El renderizador impone un orden de claves único; un escenario legal
        escrito en otro orden tiene que volver **igual**, porque el
        `execution_id` sale de los SHA-256 observados y un documento
        reordenado es otra ejecución.
        """
        reversed_general = (
            "ploidy = 2\n"
            'schema_version = "3"\n'
            'config_kind = "general"\n'
            'reference_sha256 = "' + "a" * 64 + '"\n'
            'reference_path = "references/ref.fa"\n'
            'reference_id = "ref"\n'
            f'seed = "{SEED}"\n'
            "\n"
            "[f5_scenario]\n"
            'scenario_path = "scenarios/scenario.toml"\n'
            'scenario_id = "scenario"\n'
        )
        # `sampling` después de `founder_ids`, el objetivo antes de los
        # presupuestos y `[f5]` después de las subpoblaciones: las tres cosas
        # que el renderizador emite al revés.
        reordered_scenario = (
            'config_kind = "f5_scenario"\n'
            'schema_version = "1"\n'
            "founder_count = 1\n"
            'scenario_id = "scenario"\n'
            "\n"
            "[[subpopulations]]\n"
            'population_id = "P1"\n'
            'founder_ids = ["I1"]\n'
            'sampling = "probabilistic"\n'
            "individual_count = 1\n"
            "\n"
            "[f5]\n"
            "block_count = 2\n"
            "max_layout_attempts = 8\n"
            'profile_metadata_sha256 = "' + "b" * 64 + '"\n'
            'profile_metadata_path = "profiles/profile.metadata.json"\n'
        )
        pair = import_f5_pair(
            [("g.toml", reversed_general), ("s.toml", reordered_scenario)]
        )
        exported_general, exported_scenario = f5_pair_documents(pair)
        self.assertEqual(exported_general, reversed_general)
        self.assertEqual(exported_scenario, reordered_scenario)
        # Y el orden importado **no** es el que el renderizador produce: sin
        # conservación de bytes esta prueba no diría nada.
        self.assertNotEqual(render_f5_general_toml(pair.general), reversed_general)
        self.assertNotEqual(
            render_f5_scenario_toml(pair.scenario), reordered_scenario
        )

    def test_editing_the_form_drops_the_imported_bytes(self) -> None:
        pair = import_f5_pair(
            [("general.toml", self.general), ("scenario.toml", self.scenario)]
        )
        edited = replace(pair, scenario=replace(pair.scenario, founder_count=2))
        exported_general, exported_scenario = f5_pair_documents(edited)
        self.assertEqual(exported_general, self.general)
        self.assertNotEqual(exported_scenario, self.scenario)
        self.assertIn("founder_count = 2\n", exported_scenario)
        self.assertEqual(parse_f5_scenario_toml(exported_scenario), edited.scenario)

    def test_a_non_diploid_pair_drops_the_preserved_scenario_bytes(self) -> None:
        scenario = self.scenario.replace(
            "\n[f5]\n", "\n[f5]\ninbreeding_f = 0.25\n"
        )
        pair = import_f5_pair(
            [("general.toml", self.general), ("scenario.toml", scenario)]
        )
        self.assertEqual(f5_pair_documents(pair)[1], scenario)
        tetraploid = replace(pair, general=replace(pair.general, ploidy=4))
        exported = f5_pair_documents(pair_without_inbreeding(tetraploid))[1]
        self.assertNotIn("inbreeding_f", exported)
        self.assertEqual(parse_f5_general_toml(f5_pair_documents(tetraploid)[0]).ploidy, 4)

    def test_zero_several_and_foreign_versions_are_rejected(self) -> None:
        v2 = self.general.replace('schema_version = "3"', 'schema_version = "2"')
        v4 = self.general.replace('schema_version = "3"', 'schema_version = "4"')
        scenario_v2 = self.scenario.replace('schema_version = "1"', 'schema_version = "2"')
        cases = {
            "sin documentos": [],
            "solo el general": [("general.toml", self.general)],
            "solo el escenario": [("scenario.toml", self.scenario)],
            "dos generales": [
                ("a.toml", self.general),
                ("b.toml", self.general),
                ("scenario.toml", self.scenario),
            ],
            "dos escenarios": [
                ("general.toml", self.general),
                ("a.toml", self.scenario),
                ("b.toml", self.scenario),
            ],
            "general v2": [("general.toml", v2), ("scenario.toml", self.scenario)],
            "general v4": [("general.toml", v4), ("scenario.toml", self.scenario)],
            "escenario v2": [("general.toml", self.general), ("s.toml", scenario_v2)],
        }
        for label, documents in cases.items():
            with self.subTest(case=label):
                with self.assertRaises(EditorError):
                    import_f5_pair(documents)
        # Control positivo: la pareja correcta sí se admite.
        self.assertIsNotNone(
            import_f5_pair(
                [("general.toml", self.general), ("scenario.toml", self.scenario)]
            )
        )

    def test_v2_and_v4_are_rejected_by_their_own_version_and_not_by_a_later_field(
        self,
    ) -> None:
        for version in ("2", "4"):
            with self.subTest(version=version):
                document = self.general.replace(
                    'schema_version = "3"', f'schema_version = "{version}"'
                )
                with self.assertRaises(EditorError) as raised:
                    import_f5_pair(
                        [("general.toml", document), ("scenario.toml", self.scenario)]
                    )
                self.assertIn(f'"{version}"', str(raised.exception))


class F5SubstitutionAndBudgetTests(unittest.TestCase):
    """Modo y tolerancia por bloque×subpoblación, y los siete presupuestos."""

    def test_every_budget_and_both_substitution_modes_round_trip(self) -> None:
        budgets = {
            "max_layout_attempts": 7,
            "max_haplotype_attempts": 11,
            "max_swap_proposals": 13,
            "max_haplotype_cells": 17,
            "max_operation_units": 19,
            "max_working_memory_bytes": 23,
            "max_output_bytes": 29,
        }
        # Ningún valor coincide con otro: una prueba cuyos presupuestos fueran
        # todos iguales no distinguiría un campo escrito en el sitio del otro.
        self.assertEqual(len(set(budgets.values())), len(budgets))
        base = default_f5_pair().scenario
        form = replace(
            base,
            sampling_overrides=[
                F5SamplingOverrideRow("population-1", 0, "probabilistic"),
                F5SamplingOverrideRow("population-1", 1, "controlled", "0.125"),
            ],
            **budgets,
        )
        document = render_f5_scenario_toml(form)
        for name, value in budgets.items():
            with self.subTest(budget=name):
                self.assertIn(f"{name} = {value}\n", document)
        self.assertIn('mode = "probabilistic"\n', document)
        self.assertIn("controlled_max_abs_error = 0.125\n", document)
        self.assertEqual(parse_f5_scenario_toml(document), form)

    def test_probabilistic_cannot_carry_a_tolerance(self) -> None:
        form = replace(
            default_f5_pair().scenario,
            sampling_overrides=[
                F5SamplingOverrideRow("population-1", 0, "probabilistic", "0.05")
            ],
        )
        with self.assertRaises(EditorError):
            render_f5_scenario_toml(form)

    def test_editor_rows_round_trip_through_the_substitution_table(self) -> None:
        form = replace(
            default_f5_pair().scenario,
            sampling_overrides=[
                F5SamplingOverrideRow("PA", 3, "controlled", "0.05"),
                F5SamplingOverrideRow("PB", 4, "probabilistic"),
            ],
        )
        rows = f5_override_editor_rows(form)
        rebuilt = [
            f5_override_from_editor_row(row, f"sampling_overrides[{index}]")
            for index, row in enumerate(rows)
        ]
        self.assertEqual(rebuilt, form.sampling_overrides)

    def test_editor_rows_round_trip_through_the_subpopulation_table(self) -> None:
        form = default_f5_pair().scenario
        rows = f5_subpopulation_editor_rows(form)
        rebuilt = [
            f5_subpopulation_from_editor_row(row, f"subpopulations[{index}]")
            for index, row in enumerate(rows)
        ]
        self.assertEqual(rebuilt, form.subpopulations)


class F5PresentationRuleTests(unittest.TestCase):
    """REQ-UI-15: F solo con ploidía 2; `controlled` fuera con F mayor que cero."""

    def test_inbreeding_is_offered_only_at_ploidy_two(self) -> None:
        self.assertTrue(inbreeding_is_available(2))
        for ploidy in (1, 3, 10):
            with self.subTest(ploidy=ploidy):
                self.assertFalse(inbreeding_is_available(ploidy))

    def test_a_non_diploid_pair_does_not_publish_inbreeding(self) -> None:
        pair = default_f5_pair()
        pair.scenario.inbreeding_f = "0.25"
        pair.scenario.subpopulations[0].inbreeding_f = "0.5"
        diploid = pair_without_inbreeding(pair)
        self.assertEqual(diploid.scenario.inbreeding_f, "0.25")
        pair.general.ploidy = 4
        haploidish = pair_without_inbreeding(pair)
        self.assertIsNone(haploidish.scenario.inbreeding_f)
        self.assertIsNone(haploidish.scenario.subpopulations[0].inbreeding_f)
        self.assertNotIn("inbreeding_f", render_f5_scenario_toml(haploidish.scenario))

    def test_a_document_with_controlled_and_a_positive_f_is_not_emitted(self) -> None:
        """«Desactivarse» no es anunciarse: el documento no llega a existir."""
        base = default_f5_pair().scenario
        cases = {
            "sustitución con F de escenario": replace(
                base,
                inbreeding_f="0.25",
                sampling_overrides=[
                    F5SamplingOverrideRow("population-1", 0, "controlled")
                ],
            ),
            "sustitución con F de subpoblación": replace(
                base,
                subpopulations=[
                    F5SubpopulationRow(
                        "population-1", 1, "founder_ids", ["individual-1"],
                        inbreeding_f="1",
                    )
                ],
                sampling_overrides=[
                    F5SamplingOverrideRow("population-1", 0, "controlled", "0.05")
                ],
            ),
            "subpoblación controlled con F": replace(
                base,
                inbreeding_f="0.5",
                subpopulations=[
                    F5SubpopulationRow(
                        "population-1", 1, "founder_ids", ["individual-1"],
                        sampling="controlled",
                    )
                ],
            ),
        }
        for label, form in cases.items():
            with self.subTest(case=label):
                with self.assertRaises(EditorError):
                    render_f5_scenario_toml(form)
        # Control positivo: con F igual a cero, `controlled` sí se emite, y
        # con F positivo `probabilistic` tampoco estorba.
        allowed = replace(
            base,
            inbreeding_f="0",
            sampling_overrides=[
                F5SamplingOverrideRow("population-1", 0, "controlled", "0.05")
            ],
        )
        self.assertIn('mode = "controlled"', render_f5_scenario_toml(allowed))
        tolerated = replace(
            base,
            inbreeding_f="0.25",
            sampling_overrides=[
                F5SamplingOverrideRow("population-1", 0, "probabilistic")
            ],
        )
        self.assertIn('mode = "probabilistic"', render_f5_scenario_toml(tolerated))

    def test_a_scenario_without_subpopulations_is_not_emitted(self) -> None:
        empty = replace(default_f5_pair().scenario, subpopulations=[])
        with self.assertRaises(EditorError):
            render_f5_scenario_toml(empty)

    def test_controlled_is_disabled_by_a_positive_inbreeding(self) -> None:
        self.assertTrue(controlled_is_available(None))
        self.assertTrue(controlled_is_available("0"))
        self.assertFalse(controlled_is_available("0.0001"))
        self.assertEqual(
            available_sampling_modes("0"), ("probabilistic", "controlled")
        )
        self.assertEqual(available_sampling_modes("0.5"), ("probabilistic",))


class F5ProcessTests(unittest.TestCase):
    """El adaptador ejecuta el mismo `validate`/`run` de Rust."""

    BINARY = ROOT / "target" / "debug" / SUPPORTED_BINARY_NAME

    def test_a_later_edit_invalidates_the_validated_session(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference_sha, profile_sha = write_supporting_files(root)
            general = manual_general(reference_sha)
            scenario = manual_scenario(profile_sha)
            validation, session = validate_f5_pair(self.BINARY, general, scenario, root)
            self.assertTrue(validation.succeeded, validation.stdout + validation.stderr)
            self.assertIsNotNone(session)
            edits = {
                "escenario editado": (
                    general,
                    scenario.replace("block_count = 2", "block_count = 3"),
                    "general.toml",
                ),
                "general editado": (
                    general.replace("ploidy = 2", "ploidy = 3"),
                    scenario,
                    "general.toml",
                ),
                "nombre del general editado": (general, scenario, "otro.toml"),
            }
            for label, (edited_general, edited_scenario, name) in edits.items():
                with self.subTest(case=label):
                    with self.assertRaises(ProcessAdapterError):
                        run_f5_pair(
                            session,
                            edited_general,
                            edited_scenario,
                            self.BINARY,
                            root,
                            root / f"out-{name}",
                            general_name=name,
                        )
            # Control positivo: sin edición, la misma sesión sí ejecuta.
            execution = run_f5_pair(
                session, general, scenario, self.BINARY, root, root / "out"
            )
            self.assertTrue(execution.succeeded, execution.stdout + execution.stderr)

    def test_an_existing_destination_is_refused_before_invoking_rust(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference_sha, profile_sha = write_supporting_files(root)
            general = manual_general(reference_sha)
            scenario = manual_scenario(profile_sha)
            _, session = validate_f5_pair(self.BINARY, general, scenario, root)
            (root / "taken").mkdir()
            with self.assertRaises(ProcessAdapterError):
                run_f5_pair(session, general, scenario, self.BINARY, root, root / "taken")

    def test_the_adapter_only_writes_inside_the_authoring_directory(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "authoring"
            root.mkdir()
            reference_sha, profile_sha = write_supporting_files(root)
            general = manual_general(reference_sha)
            scenario = manual_scenario(profile_sha)
            for name in ("..", "../escapado.toml", "../../escapado.toml"):
                with self.subTest(name=name):
                    with self.assertRaises(ProcessAdapterError):
                        validate_f5_pair(
                            self.BINARY, general, scenario, root, general_name=name
                        )
            self.assertFalse((Path(directory) / "escapado.toml").exists())
            for path in ("../s.toml", "/tmp/s.toml", "sub/../../s.toml"):
                with self.subTest(scenario_path=path):
                    with self.assertRaises(ProcessAdapterError):
                        validate_f5_pair(
                            self.BINARY,
                            general,
                            scenario,
                            root,
                            scenario_relative_path=path,
                        )
            # Control positivo: el nombre por omisión sí se escribe.
            validation, _ = validate_f5_pair(self.BINARY, general, scenario, root)
            self.assertTrue(validation.succeeded, validation.stdout + validation.stderr)

    def test_an_unpublishable_general_name_is_rejected_by_the_rust_admission(self) -> None:
        """La regla vive en Rust, una sola vez, y `validate` la aplica."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference_sha, profile_sha = write_supporting_files(root)
            general = manual_general(reference_sha)
            scenario = manual_scenario(profile_sha)
            validation, session = validate_f5_pair(
                self.BINARY, general, scenario, root, general_name="general con espacios.toml"
            )
            self.assertFalse(validation.succeeded)
            self.assertIsNone(session)
            self.assertIn("GPS-CONF-101", validation.stdout)
            self.assertIn("general_path", validation.stdout)
            # El adaptador no repite la regla: escribió el documento y dejó
            # que Rust lo juzgara.
            self.assertTrue((root / "general con espacios.toml").is_file())

    def test_a_foreign_document_is_not_overwritten_without_authorization(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference_sha, profile_sha = write_supporting_files(root)
            general = manual_general(reference_sha)
            scenario = manual_scenario(profile_sha)
            # El caso real: `<repo>/general.toml` es el documento **v2** que
            # las otras ocho pestañas editan y el destino que la interfaz
            # ofrece por omisión. Comparte `config_kind` con el general F5.
            legacy = (
                'config_kind = "general"\n'
                'schema_version = "2"\n'
                'seed = "' + "0" * 64 + '"\n'
            )
            (root / "general.toml").write_text(legacy, encoding="utf-8")
            with self.assertRaises(ProcessAdapterError):
                validate_f5_pair(self.BINARY, general, scenario, root)
            self.assertEqual(
                (root / "general.toml").read_text(encoding="utf-8"), legacy
            )
            (root / "general.toml").unlink()

            foreign = 'config_kind = "general"\nschema_version = "2"\n'
            (root / "scenarios" / "scenario.toml").write_text(foreign, encoding="utf-8")
            with self.assertRaises(ProcessAdapterError):
                validate_f5_pair(self.BINARY, general, scenario, root)
            self.assertEqual(
                (root / "scenarios" / "scenario.toml").read_text(encoding="utf-8"),
                foreign,
            )
            # Una autoría F5 anterior sí se reemplaza sin autorización: es el
            # documento que la propia interfaz escribió.
            (root / "scenarios" / "scenario.toml").write_text(
                'config_kind = "f5_scenario"\nschema_version = "1"\n', encoding="utf-8"
            )
            validation, session = validate_f5_pair(self.BINARY, general, scenario, root)
            self.assertTrue(validation.succeeded, validation.stdout + validation.stderr)
            self.assertIsNotNone(session)
            # Y con autorización explícita, también el documento ajeno.
            (root / "scenarios" / "scenario.toml").write_text(foreign, encoding="utf-8")
            authorized, _ = validate_f5_pair(
                self.BINARY, general, scenario, root, overwrite=True
            )
            self.assertTrue(authorized.succeeded, authorized.stdout + authorized.stderr)

    def test_the_package_summary_reads_the_manifest_and_notices_a_broken_link(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference_sha, profile_sha = write_supporting_files(root)
            general = manual_general(reference_sha)
            scenario = manual_scenario(profile_sha)
            _, session = validate_f5_pair(self.BINARY, general, scenario, root)
            execution = run_f5_pair(
                session, general, scenario, self.BINARY, root, root / "out"
            )
            self.assertTrue(execution.succeeded, execution.stdout + execution.stderr)
            summary = f5_package_summary(root / "out")
            self.assertIn("estado=complete", summary)
            for name in (
                "f5-block-layout.tsv",
                "f5-sites.tsv",
                "f5-haplotypes.tsv",
                "f5-frequencies.tsv",
                "f5-ld-evaluation.tsv",
                "f5-founders.tsv",
                "f5-attempts.tsv",
                "truth.jsonl",
            ):
                with self.subTest(artifact=name):
                    self.assertIn(name, summary)
            # Control positivo del propio resumen: si el enlace deja de
            # corresponder al fichero, el resumen lo dice en vez de repetirlo.
            founders = root / "out" / "outputs" / "f5-founders.tsv"
            founders.write_text("truncado\n", encoding="utf-8")
            self.assertIn("no corresponde", f5_package_summary(root / "out"))


class F5ParityTests(unittest.TestCase):
    """La interfaz no es otra vía de entrada: es un editor del mismo TOML."""

    BINARY = ROOT / "target" / "debug" / SUPPORTED_BINARY_NAME

    def _manual_run(self, root: Path, general: str, scenario: str, destination: Path):
        """Ejecuta el par manual sin pasar por el adaptador de la interfaz."""
        (root / "general.toml").write_text(general, encoding="utf-8")
        (root / "scenarios" / "scenario.toml").write_text(scenario, encoding="utf-8")
        return subprocess.run(
            [
                str(self.BINARY),
                "run",
                str(root / "general.toml"),
                "--diagnostics",
                "json",
                "--output",
                str(destination),
            ],
            cwd=root,
            capture_output=True,
            text=True,
            check=False,
        )

    def test_the_interface_publishes_the_same_effective_manifest_and_digests(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            manual_root = Path(directory) / "manual"
            ui_root = Path(directory) / "ui"
            manual_root.mkdir(parents=True)
            ui_root.mkdir(parents=True)
            reference_sha, profile_sha = write_supporting_files(manual_root)
            self.assertEqual((reference_sha, profile_sha), write_supporting_files(ui_root))
            general = manual_general(reference_sha)
            scenario = manual_scenario(profile_sha)

            manual = self._manual_run(
                manual_root, general, scenario, manual_root / "package"
            )
            self.assertEqual(manual.returncode, 0, manual.stdout + manual.stderr)

            pair = import_f5_pair(
                [("general.toml", general), ("scenario.toml", scenario)]
            )
            ui_general, ui_scenario = f5_pair_documents(pair)
            # Importar y exportar conserva exactamente el mismo TOML.
            self.assertEqual(ui_general, general)
            self.assertEqual(ui_scenario, scenario)
            validation, session = validate_f5_pair(
                self.BINARY, ui_general, ui_scenario, ui_root
            )
            self.assertTrue(validation.succeeded, validation.stdout + validation.stderr)
            execution = run_f5_pair(
                session, ui_general, ui_scenario, self.BINARY, ui_root, ui_root / "package"
            )
            self.assertTrue(execution.succeeded, execution.stdout + execution.stderr)

            manual_package = manual_root / "package"
            ui_package = ui_root / "package"
            # El recibo de admisión publica la ruta absoluta de cada entrada:
            # lo único que puede diferir entre dos autorías idénticas es el
            # directorio que eligió quien las escribió. Se normaliza ese
            # prefijo y nada más, y se comprueba que la normalización tocó
            # algo: si no, estaría comparando dos textos ya iguales por otra
            # razón y no diría lo que dice.
            for document in ("effective-config.json", "manifest.json"):
                with self.subTest(document=document):
                    raw_manual = (manual_package / document).read_text(encoding="utf-8")
                    raw_ui = (ui_package / document).read_text(encoding="utf-8")
                    self.assertIn(str(manual_root), raw_manual)
                    self.assertIn(str(ui_root), raw_ui)
                    self.assertEqual(
                        raw_manual.replace(str(manual_root), "<root>"),
                        raw_ui.replace(str(ui_root), "<root>"),
                    )
            manual_digests = tree_digests(manual_package)
            ui_digests = tree_digests(ui_package)
            # Diez ficheros: siete TSV, la verdad y los dos documentos raíz.
            self.assertEqual(len(manual_digests), 10)
            outputs = {
                name: digest
                for name, digest in manual_digests.items()
                if name.startswith("outputs/")
            }
            self.assertEqual(len(outputs), 8)
            self.assertEqual(
                outputs,
                {
                    name: digest
                    for name, digest in ui_digests.items()
                    if name.startswith("outputs/")
                },
            )

    def test_a_rejected_pair_produces_the_same_diagnostic_through_both_ways(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            manual_root = Path(directory) / "manual"
            ui_root = Path(directory) / "ui"
            manual_root.mkdir(parents=True)
            ui_root.mkdir(parents=True)
            reference_sha, _ = write_supporting_files(manual_root)
            write_supporting_files(ui_root)
            general = manual_general(reference_sha)
            # Un hash de perfil que no corresponde: rechazo antes de la
            # referencia, con un diagnóstico estable.
            scenario = manual_scenario("c" * 64)

            (manual_root / "general.toml").write_text(general, encoding="utf-8")
            (manual_root / "scenarios" / "scenario.toml").write_text(
                scenario, encoding="utf-8"
            )
            manual = subprocess.run(
                [
                    str(self.BINARY),
                    "validate",
                    str(manual_root / "general.toml"),
                    "--diagnostics",
                    "json",
                ],
                cwd=manual_root,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertNotEqual(manual.returncode, 0)

            pair = import_f5_pair(
                [("general.toml", general), ("scenario.toml", scenario)]
            )
            ui_general, ui_scenario = f5_pair_documents(pair)
            validation, session = validate_f5_pair(
                self.BINARY, ui_general, ui_scenario, ui_root
            )
            self.assertFalse(validation.succeeded)
            self.assertIsNone(session)
            self.assertEqual(manual.returncode, validation.returncode)

            def normalize(text: str, root: Path) -> str:
                return re.sub(re.escape(str(root)), "<root>", text)

            self.assertEqual(
                normalize(manual.stdout, manual_root),
                normalize(validation.stdout, ui_root),
            )
            self.assertTrue(manual.stdout.strip(), "el diagnóstico no puede estar vacío")
            self.assertIn("GPS-", manual.stdout)


if __name__ == "__main__":
    unittest.main()
