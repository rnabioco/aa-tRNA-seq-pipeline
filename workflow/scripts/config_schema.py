"""
Validate the merged pipeline config against workflow/schemas/config.schema.yaml.

The schema closes every section, so a misspelt or retired key is an error rather
than a silent fall-back to a default in `config.get(key, default)`. Every problem
is reported at once, and an unknown key names its nearest valid neighbour:

    charging.ml_treshold: unknown key (did you mean 'ml_threshold'?)
"""

from __future__ import annotations

import difflib
from pathlib import Path

import yaml
from jsonschema import Draft7Validator

SCHEMA_PATH = Path(__file__).resolve().parent.parent / "schemas" / "config.schema.yaml"

# Keys that were renamed or removed, with what to do instead. Checked before the
# generic unknown-key message so the user is told the replacement, not guessed at.
RETIRED = {
    "opts.escpod_align": "renamed opts.eschalign when alignment moved out of escpod",
}


def load_schema(path: Path = SCHEMA_PATH) -> dict:
    return yaml.safe_load(path.read_text())


def _dotted(parts) -> str:
    return ".".join(str(p) for p in parts) or "<top level>"


def _unknown_keys(error) -> list[str]:
    """Messages for the properties `additionalProperties: false` rejected."""
    props = error.schema.get("properties", {})
    extra = sorted(set(error.instance) - set(props), key=str)
    out = []
    for key in extra:
        path = _dotted([*error.absolute_path, key])
        if path in RETIRED:
            out.append(f"{path}: {RETIRED[path]}")
            continue
        near = difflib.get_close_matches(str(key), [str(p) for p in props], n=1)
        hint = f" (did you mean {near[0]!r}?)" if near else ""
        out.append(f"{path}: unknown key{hint}")
    return out


def config_problems(config: dict, schema: dict | None = None) -> list[str]:
    """Every way `config` departs from the schema, as readable one-liners."""
    validator = Draft7Validator(schema or load_schema())
    problems = []
    for error in sorted(validator.iter_errors(config), key=lambda e: list(e.path)):
        if error.validator == "additionalProperties":
            problems.extend(_unknown_keys(error))
        elif error.validator == "required":
            problems.append(f"{_dotted(error.absolute_path)}: {error.message}")
        elif error.validator in ("oneOf", "anyOf"):
            problems.append(
                f"{_dotted(error.absolute_path)}: {error.instance!r} is not valid "
                f"({error.schema.get('description', 'see config.schema.yaml')})"
            )
        else:
            problems.append(f"{_dotted(error.absolute_path)}: {error.message}")
    return problems


def validate_config(config: dict) -> None:
    """Exit with every problem listed, or return when the config is sound."""
    problems = config_problems(config)
    if problems:
        import sys

        sys.exit(
            "Invalid configuration (schema: workflow/schemas/config.schema.yaml):\n  "
            + "\n  ".join(problems)
        )
