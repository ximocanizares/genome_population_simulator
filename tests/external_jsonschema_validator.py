#!/usr/bin/env python3
"""Validador JSON Schema Draft 2020-12 independiente del productor Rust."""

import json
import pathlib
import sys

import jsonschema


def reject_duplicate_keys(pairs):
    document = {}
    for key, value in pairs:
        if key in document:
            raise ValueError(f"duplicate JSON key: {key}")
        document[key] = value
    return document


schema_path = pathlib.Path(sys.argv[1])
schema = json.loads(
    schema_path.read_text(encoding="utf-8"),
    object_pairs_hook=reject_duplicate_keys,
)
instance = json.loads(
    pathlib.Path(sys.argv[2]).read_text(encoding="utf-8"),
    object_pairs_hook=reject_duplicate_keys,
)
jsonschema.Draft202012Validator.check_schema(schema)
# Los contratos versionados pueden reutilizar definiciones cerradas de otro
# artefacto publicado; fijar la base física evita que el cwd del validador
# cambie el significado de una referencia relativa.
resolver = jsonschema.RefResolver(schema_path.resolve().as_uri(), schema)
validator = jsonschema.Draft202012Validator(schema, resolver=resolver)

if len(sys.argv) == 4 and sys.argv[3] == "--cases":
    for index, case in enumerate(instance):
        errors = list(validator.iter_errors(case["instance"]))
        observed = not errors
        if observed != case["expected"]:
            detail = "valid" if observed else errors[0].message
            raise AssertionError(
                f"case {index}: expected valid={case['expected']}, observed {detail}"
            )
else:
    validator.validate(instance)
