#!/usr/bin/env python3
"""Oráculo ajeno al productor Rust para el manifiesto mínimo F1."""

import json
import pathlib
import sys

import jsonschema


def require_kind(catalog: dict[str, dict], artifact_id: str, kind: str) -> None:
    assert artifact_id in catalog
    assert catalog[artifact_id]["kind"] == kind


schema_path = pathlib.Path(sys.argv[1])
instance_path = pathlib.Path(sys.argv[2])
diagnostic_catalog_path = pathlib.Path(sys.argv[3])
schema = json.loads(schema_path.read_text(encoding="utf-8"))
instance = json.loads(instance_path.read_text(encoding="utf-8"))
jsonschema.Draft202012Validator.check_schema(schema)
jsonschema.Draft202012Validator(schema).validate(instance)

diagnostic_catalog: dict[str, str] = {}
for line in diagnostic_catalog_path.read_text(encoding="utf-8").splitlines():
    if not line or line.startswith("#"):
        continue
    fields = line.split("|")
    assert len(fields) == 6
    code, _category, meaning, _context, _exit_code, _steps = fields
    assert code not in diagnostic_catalog
    diagnostic_catalog[code] = meaning

artifacts = instance["artifacts"]
catalog = {artifact["id"]: artifact for artifact in artifacts}
assert len(catalog) == len(artifacts)
assert len({artifact["path"] for artifact in artifacts}) == len(artifacts)

configuration = instance["configuration_input"]
if configuration is not None:
    require_kind(catalog, configuration["artifact_id"], "config_general_toml")
effective_id = instance["effective_config_artifact_id"]
if effective_id is not None:
    require_kind(catalog, effective_id, "effective_config")
references = {reference["id"]: reference["fasta_artifact_id"] for reference in instance["references"]}
assert len(references) == len(instance["references"])
for reference in instance["references"]:
    require_kind(catalog, reference["fasta_artifact_id"], "reference_fasta")
population_snvs = instance.get("population_snvs")
population_minimum_ids: set[str] = set()
if population_snvs is not None:
    assert population_snvs["mode"] in {"supplied", "generated", "reuse"}
    artifact_id = population_snvs.get("metadata_artifact_id")
    if population_snvs["mode"] == "reuse":
        assert artifact_id is None
    else:
        assert isinstance(artifact_id, str) and artifact_id
        require_kind(catalog, artifact_id, "population_snvs_metadata")
        assert population_snvs["metadata_sha256"] == catalog[artifact_id]["checksum"]["digest"]
        require_kind(catalog, "population-snvs-tsv", "population_snvs_tsv")
        population_minimum_ids.update({artifact_id, "population-snvs-tsv"})
reduction_ids: set[str] = set()
requested_ids: set[str] = set()
for reduction in instance["reference_reductions"]:
    assert reduction["id"] not in reduction_ids
    reduction_ids.add(reduction["id"])
    links = (
        ("fasta_artifact_id", "reference_fasta"),
        ("map_artifact_id", "reference_reduction_map"),
        ("source_fasta_artifact_id", "reference_fasta"),
    )
    for field, kind in links:
        artifact_id = reduction[field]
        if artifact_id is not None:
            require_kind(catalog, artifact_id, kind)
            if field != "source_fasta_artifact_id":
                requested_ids.add(artifact_id)
    if reduction["state"] == "complete":
        assert all(reduction[field] is not None for field, _ in links)
    if (
        reduction["fasta_artifact_id"] is not None
        and reduction["source_fasta_artifact_id"] is not None
    ):
        assert reduction["fasta_artifact_id"] != reduction["source_fasta_artifact_id"]
    if reduction["fasta_artifact_id"] is not None:
        assert references[reduction["id"]] == reduction["fasta_artifact_id"]
    if reduction["source_fasta_artifact_id"] is not None:
        assert references[reduction["source_reference_id"]] == reduction["source_fasta_artifact_id"]

classified: list[str] = []
for category in ("minimum", "defaulted", "requested", "disabled"):
    classified.extend(instance["output_categories"][category])
assert len(classified) == len(set(classified))
assert set(classified) == set(catalog)
minimum_ids: set[str] = set()
if configuration is not None:
    minimum_ids.add(configuration["artifact_id"])
if effective_id is not None:
    minimum_ids.add(effective_id)
for reference_id, artifact_id in references.items():
    if reference_id not in reduction_ids:
        minimum_ids.add(artifact_id)
minimum_ids |= population_minimum_ids
assert set(instance["output_categories"]["minimum"]) == minimum_ids
assert set(instance["output_categories"]["requested"]) == requested_ids
assert not instance["output_categories"]["defaulted"]
assert not instance["output_categories"]["disabled"]
assert minimum_ids | requested_ids == set(catalog)

if instance["state"] == "complete":
    assert configuration is not None
    assert effective_id is not None
    assert instance["references"]
    assert not instance["diagnostics"]
    assert all(reduction["state"] == "complete" for reduction in instance["reference_reductions"])
else:
    assert instance["diagnostics"]

for diagnostic in instance["diagnostics"]:
    assert diagnostic["code"] in diagnostic_catalog
    assert diagnostic["message"] == diagnostic_catalog[diagnostic["code"]]
