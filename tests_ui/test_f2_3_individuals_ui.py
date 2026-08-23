from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "ui"))

from config_editor import (  # noqa: E402
    DirectEventRow,
    EditorError,
    InitialPopulationRow,
    PopulationSnvsForm,
    default_form,
    direct_event_from_editor_row,
    initial_population_from_editor_row,
    parse_general_toml,
    render_general_toml,
)
from process_adapter import (  # noqa: E402
    SUPPORTED_BINARY_NAME,
    ValidatedSession,
    diagnostic_summary,
    initial_population_statistics_summary,
    run_toml,
    validate_toml,
)
from ui_panels import (  # noqa: E402
    authored_editor_rows,
    effective_incompatible_genotype_authorization,
    explicit_snp_authorization_summary,
)


class F23IndividualsUiTests(unittest.TestCase):
    BINARY = ROOT / "target" / "debug" / SUPPORTED_BINARY_NAME

    def test_authorization_participates_in_validation_session_identity(self) -> None:
        session = ValidatedSession(
            "config", Path("/binary"), "0.1.0", Path("/base"), False, True
        )

        self.assertTrue(session.matches("config", Path("/binary"), Path("/base"), False, True))
        self.assertFalse(session.matches("config", Path("/binary"), Path("/base"), False, False))

    def test_authorization_is_effective_only_for_non_placeholder_explicit_snps(self) -> None:
        insertion = {"omit_if_blank": False, "kind": "insertion", "event_id": "i1"}
        snp = {
            "omit_if_blank": True,
            "kind": "snp",
            "event_id": "s1",
            "individual_id": "sample-A",
            "contig": "chr1",
            "position_1based": 42,
        }
        placeholder_snp = {"omit_if_blank": True, "kind": "snp", "event_id": ""}

        self.assertFalse(effective_incompatible_genotype_authorization([insertion], True))
        self.assertFalse(effective_incompatible_genotype_authorization([], True))
        self.assertFalse(effective_incompatible_genotype_authorization([placeholder_snp], True))
        self.assertTrue(effective_incompatible_genotype_authorization([snp], True))
        self.assertIn("sample-A@chr1:42", explicit_snp_authorization_summary([snp]))
        self.assertIn("no se aplicará", explicit_snp_authorization_summary([insertion]))

        session = ValidatedSession("config", Path("/binary"), "0.1.0", Path("/base"))
        self.assertTrue(session.matches("config", Path("/binary"), Path("/base"), False, False))
        self.assertFalse(session.matches("config", Path("/binary"), Path("/base"), False, True))

    def test_only_explicitly_marked_blank_placeholders_are_filtered(self) -> None:
        rows = [
            {"population_id": "", "omit_if_blank": True},
            {"population_id": "", "omit_if_blank": False, "individual_count": 1},
            {"population_id": "", "omit_if_blank": True, "individual_count": 2},
        ]

        self.assertEqual(
            authored_editor_rows(rows, "population_id"),
            [
                {"population_id": "", "omit_if_blank": False, "individual_count": 1},
                {"population_id": "", "omit_if_blank": True, "individual_count": 2},
            ],
        )

    def test_empty_imported_ids_are_indexed_errors_or_reach_rust_without_filtering(self) -> None:
        source = """config_kind = "general"
schema_version = "1"
[reference]
kind = "simulated"
id = "reference-main"
[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 1000
[[initial_populations]]
population_id = ""
individual_count = 1
ploidy = 2
genotype_model = "hwe"
founder_keys = ["founder-empty"]
[[direct_events]]
event_id = ""
individual_id = "sample-A"
phase = 0
kind = "insertion"
contig = "chr1"
point = 1
bases = "A"
"""
        imported = parse_general_toml(source)
        population_row = {
            "population_id": imported.initial_populations[0].population_id,
            "individual_count": 1,
            "ploidy": 2,
            "genotype_model": "hwe",
            "inbreeding_coefficient": "",
            "individual_ids": "",
            "founder_keys": "founder-empty",
            "individual_ids_specified": False,
            "omit_if_blank": False,
        }
        event_row = {
            "event_id": imported.direct_events[0].event_id,
            "individual_id": "sample-A",
            "phase": 0,
            "kind": "insertion",
            "contig": "chr1",
            "point": 1,
            "bases": "A",
            "omit_if_blank": False,
        }
        form = default_form()
        form.initial_populations = [
            initial_population_from_editor_row(
                authored_editor_rows([population_row], "population_id")[0],
                "initial_populations[0]",
            )
        ]
        form.direct_events = [
            direct_event_from_editor_row(
                authored_editor_rows([event_row], "event_id")[0], "direct_events[0]"
            )
        ]
        rendered = render_general_toml(form)

        self.assertIn('population_id = ""', rendered)
        self.assertIn('event_id = ""', rendered)
        validation, session = validate_toml(self.BINARY, rendered, ROOT)
        self.assertFalse(validation.succeeded)
        self.assertIsNone(session)
        with self.assertRaises(EditorError) as raised:
            parse_general_toml(source.replace('genotype_model = "hwe"', 'genotype_model = "hwe"\nindividual_ids = []'))
        self.assertEqual(raised.exception.path, "initial_populations[0].individual_ids")

    def test_editor_rows_cover_hwe_f_and_every_direct_event_kind(self) -> None:
        hwe = initial_population_from_editor_row(
            {
                "population_id": "P1",
                "individual_count": 1,
                "ploidy": 2,
                "genotype_model": "hwe",
                "inbreeding_coefficient": "",
                "individual_ids": "sample-A",
                "founder_keys": "",
            },
            "initial_populations[0]",
        )
        f_diploid = initial_population_from_editor_row(
            {
                "population_id": "P2",
                "individual_count": 1,
                "ploidy": 2,
                "genotype_model": "f_diploid",
                "inbreeding_coefficient": "0.25",
                "individual_ids": "sample-B",
                "founder_keys": "",
            },
            "initial_populations[1]",
        )
        events = [
            direct_event_from_editor_row({"event_id": "s", "individual_id": "sample-A", "phase": 0, "kind": "snp", "contig": "chr1", "position_1based": 1, "ref": "A", "alt": "G"}, "direct_events[0]"),
            direct_event_from_editor_row({"event_id": "i", "individual_id": "sample-A", "phase": 0, "kind": "insertion", "contig": "chr1", "point": 1, "bases": "A"}, "direct_events[1]"),
            direct_event_from_editor_row({"event_id": "d", "individual_id": "sample-A", "phase": 0, "kind": "deletion", "contig": "chr1", "start": 1, "end": 2}, "direct_events[2]"),
        ]
        form = default_form()
        form.initial_populations = [hwe, f_diploid]
        form.direct_events = events

        rendered = render_general_toml(form)

        self.assertEqual(parse_general_toml(rendered), form)
        for kind in ("snp", "insertion", "deletion"):
            self.assertIn(f'kind = "{kind}"', rendered)

    def test_blank_individual_ids_requires_founder_keys_for_derived_ids(self) -> None:
        imported = parse_general_toml(
            """config_kind = "general"
schema_version = "1"
[reference]
kind = "simulated"
id = "reference-main"
[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 10
[[initial_populations]]
population_id = "P1"
individual_count = 1
ploidy = 2
genotype_model = "hwe"
individual_ids = ["sample-A"]
"""
        )
        imported_population = imported.initial_populations[0]
        edited_import_row = {
            "population_id": imported_population.population_id,
            "individual_count": imported_population.individual_count,
            "ploidy": imported_population.ploidy,
            "genotype_model": imported_population.genotype_model,
            "inbreeding_coefficient": "",
            "individual_ids": "",  # La persona usuaria borra la lista importada.
            "founder_keys": "founder-P1",
            "individual_ids_specified": True,
        }
        form = default_form()
        form.initial_populations = [
            initial_population_from_editor_row(
                edited_import_row,
                "initial_populations[0]",
            )
        ]

        rendered = render_general_toml(form)

        self.assertNotIn("individual_ids", rendered)
        self.assertIn('founder_keys = ["founder-P1"]', rendered)
        self.assertEqual(
            parse_general_toml(rendered).initial_populations[0].founder_keys,
            ["founder-P1"],
        )

    def test_individual_diagnostics_are_localized_to_the_individuals_tab(self) -> None:
        rendered = diagnostic_summary(
            {
                "diagnostics": [
                    {
                        "code": "GPS-TEST",
                        "message": "error",
                        "context": {"field": "direct_events.phase"},
                    }
                ]
            }
        )

        self.assertIn("Individuos · direct_events.phase", rendered)

    def test_founder_summary_reads_only_small_effective_configuration(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            (output / "effective-config.json").write_text(
                json.dumps(
                    {
                        "population_snvs": {
                            "initial_population_statistics": {
                                "populations": [
                                    {
                                        "population_id": "P1",
                                        "individual_count": 2,
                                        "ploidy": 2,
                                        "genotype_model": "hwe",
                                        "loci": [{"variant_id": "v1"}],
                                    }
                                ]
                            }
                        }
                    }
                )
            )

            summary = initial_population_statistics_summary(output)

        self.assertEqual(summary, "P1: N=2, ploidía=2, modelo=hwe, loci=1")

    def test_founder_summary_rejects_an_oversized_effective_configuration(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            (output / "effective-config.json").write_bytes(b" " * 2_000_001)

            summary = initial_population_statistics_summary(output)

        self.assertIn("supera 2 MB", summary)

    @unittest.skipUnless(BINARY.exists(), "compile the Rust binary before this integration test")
    def test_ui_toml_with_hwe_f_and_all_direct_events_runs_and_exposes_bounded_summary(self) -> None:
        form = default_form()
        form.reference.composition = {"A": "1", "C": "0", "G": "0", "T": "0"}
        form.population_snvs = PopulationSnvsForm(
            "generate",
            None,
            None,
            {
                "count_mode": "exact",
                "locus_count": 0,
                "spatial_mode": "uniform",
                "substitution_profile": "uniform_snv12_v1",
                "populations": [
                    {
                        "population_id": "P1",
                        "frequency_model": "uniform",
                        "min_alt_af": 0.1,
                        "max_alt_af": 0.9,
                    },
                    {
                        "population_id": "P2",
                        "frequency_model": "uniform",
                        "min_alt_af": 0.1,
                        "max_alt_af": 0.9,
                    },
                ],
            },
        )
        initial_rows = [
            {
                "population_id": "P1",
                "individual_count": 1,
                "ploidy": 2,
                "genotype_model": "hwe",
                "inbreeding_coefficient": "",
                "individual_ids": "sample-A",
                "founder_keys": "",
                "individual_ids_specified": True,
                "omit_if_blank": True,
            },
            {
                "population_id": "P2",
                "individual_count": 1,
                "ploidy": 2,
                "genotype_model": "f_diploid",
                "inbreeding_coefficient": "0.25",
                "individual_ids": "sample-B",
                "founder_keys": "",
                "individual_ids_specified": True,
                "omit_if_blank": True,
            },
        ]
        event_rows = [
            {
                "event_id": "event-snp",
                "individual_id": "sample-A",
                "phase": 0,
                "kind": "snp",
                "contig": "chr1",
                "position_1based": 1,
                "ref": "A",
                "alt": "G",
                "point": "",
                "bases": "",
                "start": "",
                "end": "",
                "omit_if_blank": True,
            },
            {
                "event_id": "event-insertion",
                "individual_id": "sample-A",
                "phase": 0,
                "kind": "insertion",
                "contig": "chr1",
                "position_1based": "",
                "ref": "",
                "alt": "",
                "point": 1,
                "bases": "A",
                "start": "",
                "end": "",
                "omit_if_blank": True,
            },
            {
                "event_id": "event-deletion",
                "individual_id": "sample-A",
                "phase": 0,
                "kind": "deletion",
                "contig": "chr1",
                "position_1based": "",
                "ref": "",
                "alt": "",
                "point": "",
                "bases": "",
                "start": 1,
                "end": 2,
                "omit_if_blank": True,
            },
        ]
        form.initial_populations = [
            initial_population_from_editor_row(row, f"initial_populations[{index}]")
            for index, row in enumerate(authored_editor_rows(initial_rows, "population_id"))
        ]
        form.direct_events = [
            direct_event_from_editor_row(row, f"direct_events[{index}]")
            for index, row in enumerate(authored_editor_rows(event_rows, "event_id"))
        ]
        toml = render_general_toml(form)

        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            validation, session = validate_toml(
                self.BINARY,
                toml,
                base,
                allow_empty_snv_matrix=True,
                allow_incompatible_genotypes=True,
            )
            self.assertTrue(validation.succeeded, validation.stdout or validation.stderr)
            self.assertIsNotNone(session)
            execution = run_toml(
                session,
                toml,
                self.BINARY,
                base,
                base / "result",
                allow_empty_snv_matrix=True,
                allow_incompatible_genotypes=True,
            )
            self.assertTrue(execution.succeeded, execution.stdout or execution.stderr)
            self.assertIn(
                "P1: N=1, ploidía=2, modelo=hwe, loci=0",
                initial_population_statistics_summary(base / "result"),
            )


if __name__ == "__main__":
    unittest.main()
