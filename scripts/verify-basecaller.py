#!/usr/bin/env python3
"""Prove the vendored basecalling model IS the one the charging bundle names.

The DAG-time check in `common.smk` compares NAMES, which is cheap and catches
the configuration mistake. This proves byte identity against the `model_sha256`
the bundle declares, which catches the other thing: a model directory that is
named right and is not the same bytes — a partial download, a re-fetch of a
retagged upstream model, an edited config.toml.

It hashes the whole model (~300 MB, 159 files for rna004_130bps_sup@v5.3.0),
which is why it is a task and not part of DAG construction.

The scheme is upstream's, reproduced exactly — sha256 over every file in the
directory, name-sorted, each contributing its POSIX relative path then its bytes
(`escapepod-models/scripts/charging/verify_basecaller.py`). Any deviation
produces a hash that can never match, which would read as tampering.

    pixi run verify-basecaller

Exit 0 when the declared hash matches or the bundle declares nothing; 1 on a
mismatch or a missing model.
"""

import os
import sys
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "workflow" / "scripts"))

from basecaller_compat import bundle_basecaller, model_digest  # noqa: E402


def resolve(path):
    """Config paths are repo-relative so they stay portable."""
    path = str(path)
    return Path(path) if os.path.isabs(path) else REPO_ROOT / path


def main():
    config = yaml.safe_load((REPO_ROOT / "config" / "config-base.yml").read_text())

    bundle = resolve(config["charging"]["model"])
    model = resolve(config["base_calling_model"])
    print(f"bundle: {bundle.name}")
    print(f"model : {model.name}")

    declared = bundle_basecaller(bundle)
    if declared is None:
        print(
            "\nThis bundle declares no `basecaller` block, so there is nothing to\n"
            "verify against. Bundles published before escapepod-models#106 do not\n"
            "carry one — that is silence, not a claim of compatibility."
        )
        return 0

    if not model.is_dir():
        print(f"\nFAILED: basecalling model directory not found: {model}")
        print("Run `pixi run setup` to download it.")
        return 1

    if model.name != declared["model"]:
        print(
            f"\nFAILED: name mismatch — the bundle names {declared['model']}, "
            f"config basecalls with {model.name}."
        )
        print("The DAG-time check reports this too; see charging.basecaller_check.")
        return 1

    got, count = model_digest(model)
    want = declared["model_sha256"]
    print(f"files : {count}")
    print(f"sha256: {got}")

    if got != want:
        print(f"\nFAILED: declared {want}")
        print(
            "\nSame name, different bytes. That is a partial download, a re-fetch "
            "of a\nretagged upstream model, or a locally edited file — not a "
            "configuration\nmistake, and the name check cannot see it. Re-download "
            "the model with\n`pixi run setup` and re-run this."
        )
        return 1

    print(f"\nOK: byte-identical to the model {bundle.name} was trained on.")
    print(f"    (bundle built with dorado {declared['dorado_version']}; "
          f"config pins {config.get('dorado_version')})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
