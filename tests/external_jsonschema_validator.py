#!/usr/bin/env python3
"""Validador JSON Schema Draft 2020-12 independiente del productor Rust."""

import json
import pathlib
import sys

import jsonschema


schema = json.loads(pathlib.Path(sys.argv[1]).read_text(encoding="utf-8"))
instance = json.loads(pathlib.Path(sys.argv[2]).read_text(encoding="utf-8"))
jsonschema.Draft202012Validator.check_schema(schema)
jsonschema.Draft202012Validator(schema).validate(instance)
