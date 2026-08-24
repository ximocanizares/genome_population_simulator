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


schema = json.loads(
    pathlib.Path(sys.argv[1]).read_text(encoding="utf-8"),
    object_pairs_hook=reject_duplicate_keys,
)
instance = json.loads(
    pathlib.Path(sys.argv[2]).read_text(encoding="utf-8"),
    object_pairs_hook=reject_duplicate_keys,
)
jsonschema.Draft202012Validator.check_schema(schema)
jsonschema.Draft202012Validator(schema).validate(instance)
