#!/usr/bin/env python3
"""Report vendored escapepod model bundles that upstream has since superseded.

`verify-demux-model.py` answers "is this copy byte-for-byte what upstream
released?". This answers the other half: "is that still the release we want?".
A bundle can verify perfectly and be two retrains out of date.

It is a REPORT, never an upgrade. Newest is not automatically best here, and the
shipped ldx16 pin is the standing example: upstream's own v0.1.1 rebuild scores
-0.0038 exact match against the v0.1.0 we vendor, upstream is still bounding
that gap as training-seed variance, and v0.1.1 deliberately did not become the
shipped bundle. Moving a barcode model is also not a neutral act — nbc16 to
ldx16 changed ~8% of calls — so a human decides, and records the decision in
`resources/models/pins.yml`.

Network, and therefore NOT part of the workflow. Upstream is a private repo and
compute nodes have no route to GitHub, which is why the bundles are vendored in
the first place; run this from a login or orchestration node:

    pixi run check-models

Exit status: 0 when every bundle is current or deliberately pinned, 1 when one
has drifted with no recorded reason (so this can gate CI). Unreachable upstream
is reported and exits 0 — the check is advisory by default — unless --strict.
"""

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

import yaml

REPO = "rnabioco/escapepod-models"
REPO_ROOT = Path(__file__).resolve().parents[1]
MODEL_DIRS = [REPO_ROOT / "resources" / "models" / d for d in ("demux", "charging")]
PINS = REPO_ROOT / "resources" / "models" / "pins.yml"

# A vendored bundle directory, and an upstream release tag, are the same shape:
# `<family>@v<semver>`. That is what makes this comparison a lookup rather than
# a mapping table someone has to keep in step.
BUNDLE = re.compile(r"^(?P<family>.+)@v(?P<version>\d+\.\d+\.\d+)$")


def version_key(version):
    return tuple(int(part) for part in version.split("."))


def vendored_bundles():
    """Every vendored bundle, as {family: (version, path)}."""
    found = {}
    for models in MODEL_DIRS:
        if not models.is_dir():
            continue
        for path in sorted(models.iterdir()):
            match = BUNDLE.match(path.name)
            if path.is_dir() and match:
                found[match["family"]] = (match["version"], path)
    return found


def upstream_latest():
    """Newest release per family, as {family: version}.

    Raises RuntimeError with the underlying message rather than exiting, so the
    caller decides whether being offline is fatal.
    """
    try:
        out = subprocess.run(
            ["gh", "release", "list", "--repo", REPO, "--limit", "300",
             "--json", "tagName"],
            capture_output=True, text=True, check=True,
        ).stdout
    except FileNotFoundError as exc:
        raise RuntimeError("the `gh` CLI is not installed") from exc
    except subprocess.CalledProcessError as exc:
        detail = (exc.stderr or "").strip().splitlines()
        raise RuntimeError(detail[-1] if detail else f"gh exited {exc.returncode}")

    latest = {}
    for release in json.loads(out):
        match = BUNDLE.match(release["tagName"])
        if not match:
            continue  # repo-level tags like `v0.2.0`, which name no bundle
        family, version = match["family"], match["version"]
        if family not in latest or version_key(version) > version_key(latest[family]):
            latest[family] = version
    return latest


def load_pins():
    """Deliberate holds, as {family: {"version": ..., "reason": ...}}."""
    if not PINS.is_file():
        return {}
    return yaml.safe_load(PINS.read_text()) or {}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--strict", action="store_true",
        help="treat an unreachable upstream as a failure rather than a notice",
    )
    args = parser.parse_args()

    vendored = vendored_bundles()
    if not vendored:
        print("No vendored bundles found under resources/models/.")
        return 0

    try:
        latest = upstream_latest()
    except RuntimeError as exc:
        print(f"Could not read {REPO} releases: {exc}")
        print("\nThis check needs network and a GitHub login, so it belongs on a")
        print("login or orchestration node — not on a compute node, which has no")
        print("route to GitHub. Vendored bundles are unaffected either way.")
        return 1 if args.strict else 0

    pins = load_pins()
    drifted = []
    for family, (version, path) in sorted(vendored.items()):
        newest = latest.get(family)
        pin = pins.get(family) or {}

        if newest is None:
            status = "no upstream release found — is the family renamed?"
            drifted.append(family)
        elif version_key(version) >= version_key(newest):
            status = "current"
        elif pin.get("version") == version:
            status = f"pinned at v{version} (upstream v{newest}) — deliberate"
        elif pin:
            status = (
                f"BEHIND: v{version} vendored, v{newest} upstream; "
                f"pins.yml pins v{pin.get('version')}, which is NOT what is here"
            )
            drifted.append(family)
        else:
            status = f"BEHIND: v{version} vendored, v{newest} upstream"
            drifted.append(family)

        print(f"{family}")
        print(f"  {status}")
        if pin.get("reason") and pin.get("version") == version:
            reason = " ".join(pin["reason"].split())
            print(f"  reason: {reason}")

    if drifted:
        print(
            f"\n{len(drifted)} bundle(s) need a decision: "
            + ", ".join(drifted)
        )
        print(
            "Either vendor the newer release (and re-run "
            "`pixi run verify-demux-model`), or record why not in "
            f"{PINS.relative_to(REPO_ROOT)}."
        )
        print(
            "Do not upgrade a barcode model casually: nbc16 -> ldx16 changed "
            "~8% of calls, so a switch mid-project makes runs incomparable."
        )
        return 1

    print("\nEvery vendored bundle is current or deliberately pinned.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
