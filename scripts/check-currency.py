#!/usr/bin/env python3
"""Report pins that upstream has superseded: the vendored model bundles under
resources/models/, and the tool versions config-base.yml pins (escpod, dorado).

`verify-demux-model.py` answers "is this copy byte-for-byte what upstream
released?". This answers the other half: "is that still the release we want?".
A bundle can verify perfectly and be two retrains out of date; an escpod pin can
sit on a release for weeks after the one that fixed the bug we are working
around shipped -- which is how the pipeline was still on 0.19.0 the day after
0.20.0 came out.

It is a REPORT, never an upgrade. Newest is not automatically best here, and the
shipped ldx16 pin is the standing example: upstream's own v0.1.1 rebuild scores
-0.0038 exact match against the v0.1.0 we vendor, upstream is still bounding
that gap as training-seed variance, and v0.1.1 deliberately did not become the
shipped bundle. Moving a barcode model is also not a neutral act -- nbc16 to
ldx16 changed ~8% of calls -- and a dorado bump is a basecaller change the
charging bundle was not trained on. So a human decides, and records the decision
in `resources/models/pins.yml` (model families at the top level, tools under
`tools:`). A recorded hold must name the version actually pinned, so a later
bump that leaves the file untouched is reported as stale rather than silently
inheriting the old justification.

Two more things are checked for escpod, because a bump is two edits that have
to agree: the SHA256 constants in scripts/setup-tools.sh are compared against
the pinned release's own SHA256SUMS.txt, so a pin moved without its checksums
(or checksums moved without the pin) fails here rather than at the next
`pixi run setup`.

Network, and therefore NOT part of the workflow. escapepod-models is a private
repo and compute nodes have no route to GitHub, which is why the bundles are
vendored in the first place; run this from a login or orchestration node, or
let .github/workflows/currency.yml run it weekly:

    pixi run check-currency               # everything
    pixi run check-models                 # bundles only
    pixi run check-escpod                 # tool pins only (public repos)

Exit status: 0 when every pin is current or deliberately held, 1 when one has
drifted with no recorded reason or a checksum disagrees with the release (so
this can gate CI). Unreachable upstream is reported and exits 0 -- the check is
advisory by default -- unless --strict.
"""

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
MODELS_REPO = "rnabioco/escapepod-models"
MODEL_DIRS = [REPO_ROOT / "resources" / "models" / d for d in ("demux", "charging")]
PINS = REPO_ROOT / "resources" / "models" / "pins.yml"
CONFIG = REPO_ROOT / "config" / "config-base.yml"
SETUP = REPO_ROOT / "scripts" / "setup-tools.sh"

# A vendored bundle directory, and an upstream release tag, are the same shape:
# `<family>@v<semver>`. That is what makes this comparison a lookup rather than
# a mapping table someone has to keep in step.
BUNDLE = re.compile(r"^(?P<family>.+)@v(?P<version>\d+\.\d+\.\d+)$")
# A tool release tag: `v<semver>`. The newest is the highest version, not the
# most recently published -- dorado keeps two lines alive (1.3.x, 2.1.x) and
# publishes patch releases on the old one after the new one.
SEMVER_TAG = re.compile(r"^v(?P<version>\d+\.\d+\.\d+)$")

# The tools config-base.yml pins, and where upstream publishes them. escpod's
# checksums map setup-tools.sh constants to the artifact each one guards.
TOOLS = {
    "escpod": {
        "repo": "rnabioco/escapepod-rs",
        "config_key": "escpod_version",
        "checksums": {
            "ESCPOD_SHA256_x86_64_linux": "escpod-v{v}-x86_64-unknown-linux-musl.tar.gz",
            "ESCPOD_SHA256_aarch64_linux": "escpod-v{v}-aarch64-unknown-linux-musl.tar.gz",
            "ESCPOD_SHA256_x86_64_darwin": "escpod-v{v}-x86_64-apple-darwin.tar.gz",
            "ESCPOD_SHA256_aarch64_darwin": "escpod-v{v}-aarch64-apple-darwin.tar.gz",
            "ESCPOD_SHA256_x86_64_linux_gpu": "escpod-v{v}-x86_64-unknown-linux-gnu-gpu.tar.gz",
        },
    },
    "dorado": {
        "repo": "nanoporetech/dorado",
        "config_key": "dorado_version",
    },
}


def version_key(version):
    return tuple(int(part) for part in str(version).split("."))


def gh(*args, text=True):
    """Run a `gh` command; RuntimeError with the underlying message on failure.

    Raised rather than exiting so the caller decides whether being offline is
    fatal (--strict) or a notice.
    """
    try:
        return subprocess.run(
            ["gh", *args], capture_output=True, text=text, check=True
        ).stdout
    except FileNotFoundError as exc:
        raise RuntimeError("the `gh` CLI is not installed") from exc
    except subprocess.CalledProcessError as exc:
        detail = (exc.stderr or "").strip().splitlines()
        raise RuntimeError(
            detail[-1] if detail else f"gh exited {exc.returncode}"
        ) from exc


# ----------------------------------------------------------------- models ---


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


def upstream_latest_bundles():
    """Newest release per family in escapepod-models, as {family: version}."""
    out = gh(
        "release", "list", "--repo", MODELS_REPO, "--limit", "300", "--json", "tagName"
    )
    latest = {}
    for release in json.loads(out):
        match = BUNDLE.match(release["tagName"])
        if not match:
            continue  # repo-level tags like `v0.2.0`, which name no bundle
        family, version = match["family"], match["version"]
        if family not in latest or version_key(version) > version_key(latest[family]):
            latest[family] = version
    return latest


# ------------------------------------------------------------------ tools ---


def upstream_latest_tag(repo):
    """The highest `v<semver>` release of a tool repo, or None."""
    out = gh(
        "release",
        "list",
        "--repo",
        repo,
        "--limit",
        "100",
        "--json",
        "tagName,isPrerelease,isDraft",
    )
    versions = [
        SEMVER_TAG.match(r["tagName"])["version"]
        for r in json.loads(out)
        if SEMVER_TAG.match(r["tagName"]) and not r["isPrerelease"] and not r["isDraft"]
    ]
    return max(versions, key=version_key) if versions else None


def release_checksums(repo, version):
    """{artifact: sha256} from the release's SHA256SUMS.txt."""
    out = gh(
        "release",
        "download",
        f"v{version}",
        "--repo",
        repo,
        "--pattern",
        "SHA256SUMS.txt",
        "--output",
        "-",
    )
    sums = {}
    for line in out.splitlines():
        parts = line.split()
        if len(parts) == 2:
            sums[parts[1]] = parts[0]
    return sums


def setup_checksums():
    """{constant name: sha256} as written in scripts/setup-tools.sh."""
    text = SETUP.read_text()
    return dict(
        re.findall(r'^(ESCPOD_SHA256_\w+)="([0-9a-f]{64})"', text, re.MULTILINE)
    )


def config_versions():
    cfg = yaml.safe_load(CONFIG.read_text()) or {}
    return {
        tool["config_key"]: str(cfg.get(tool["config_key"], ""))
        for tool in TOOLS.values()
    }


# ------------------------------------------------------------------ shared ---


def load_pins():
    """Recorded holds: ({family: {version, reason}}, {tool: {version, reason}})."""
    if not PINS.is_file():
        return {}, {}
    pins = yaml.safe_load(PINS.read_text()) or {}
    tools = pins.pop("tools", None) or {}
    return pins, tools


def judge(have, newest, pin):
    """(status line, drifted?) for one pin against upstream and its recorded hold."""
    if newest is None:
        return "no upstream release found — is it renamed?", True
    if version_key(have) >= version_key(newest):
        return "current", False
    if pin.get("version") is not None and str(pin["version"]) == str(have):
        return f"pinned at v{have} (upstream v{newest}) — deliberate", False
    if pin:
        return (
            f"BEHIND: v{have} here, v{newest} upstream; pins.yml holds "
            f"v{pin.get('version')}, which is NOT what is here"
        ), True
    return f"BEHIND: v{have} here, v{newest} upstream", True


def report(name, have, status, pin):
    print(f"{name}")
    print(f"  {status}")
    if (
        pin.get("reason")
        and str(pin.get("version")) == str(have)
        and "deliberate" in status
    ):
        print(f"  reason: {' '.join(pin['reason'].split())}")


def offline_notice(what, exc):
    print(f"Could not read {what}: {exc}")
    print("  This check needs network and a GitHub login, so it belongs on a")
    print("  login or orchestration node -- not on a compute node, which has no")
    print("  route to GitHub. Pinned copies are unaffected either way.")


def check_models(pins):
    """Returns (drifted families, unreachable?)."""
    vendored = vendored_bundles()
    if not vendored:
        print("No vendored bundles found under resources/models/.")
        return [], False
    try:
        latest = upstream_latest_bundles()
    except RuntimeError as exc:
        offline_notice(f"{MODELS_REPO} releases", exc)
        return [], True

    drifted = []
    for family, (version, _path) in sorted(vendored.items()):
        pin = pins.get(family) or {}
        status, bad = judge(version, latest.get(family), pin)
        report(family, version, status, pin)
        if bad:
            drifted.append(family)
    return drifted, False


def check_tools(tool_pins):
    """Returns (drifted tools, checksum problems, unreachable?)."""
    versions = config_versions()
    drifted, problems, unreachable = [], [], False
    for name, tool in TOOLS.items():
        have = versions[tool["config_key"]]
        pin = tool_pins.get(name) or {}
        try:
            newest = upstream_latest_tag(tool["repo"])
        except RuntimeError as exc:
            offline_notice(f"{tool['repo']} releases", exc)
            unreachable = True
            continue
        status, bad = judge(have, newest, pin)
        report(f"{name} ({tool['config_key']}: {have})", have, status, pin)
        if bad:
            drifted.append(name)

        if "checksums" in tool:
            try:
                published = release_checksums(tool["repo"], have)
            except RuntimeError as exc:
                offline_notice(f"{tool['repo']} v{have} SHA256SUMS.txt", exc)
                unreachable = True
                continue
            ours = setup_checksums()
            for const, artifact in tool["checksums"].items():
                want = published.get(artifact.format(v=have))
                got = ours.get(const)
                if want is None:
                    problems.append(
                        f"{const}: release v{have} ships no {artifact.format(v=have)}"
                    )
                elif got != want:
                    problems.append(
                        f"{const} in {SETUP.name} is {got or 'missing'}, but "
                        f"release v{have} publishes {want} for {artifact.format(v=have)}"
                    )
            if not [p for p in problems if p.startswith("ESCPOD")]:
                print(f"  setup-tools.sh checksums match the v{have} release")
    return drifted, problems, unreachable


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--strict",
        action="store_true",
        help="treat an unreachable upstream as a failure rather than a notice",
    )
    parser.add_argument(
        "--only",
        choices=("models", "tools"),
        help="check only the vendored bundles, or only the tool pins",
    )
    args = parser.parse_args()

    model_pins, tool_pins = load_pins()
    drifted, problems, unreachable = [], [], False

    if args.only != "tools":
        print("== vendored model bundles (resources/models)")
        d, u = check_models(model_pins)
        drifted += d
        unreachable |= u
        print()

    if args.only != "models":
        print("== tool pins (config/config-base.yml)")
        d, p, u = check_tools(tool_pins)
        drifted += d
        problems += p
        unreachable |= u
        print()

    if problems:
        print("CHECKSUMS DISAGREE WITH THE RELEASE:")
        for p in problems:
            print(f"  {p}")
        print(
            "A version pin and its checksums move together. Copy the values from "
            "the release's SHA256SUMS.txt into scripts/setup-tools.sh."
        )
    if drifted:
        print(f"{len(drifted)} pin(s) need a decision: " + ", ".join(drifted))
        print(
            "Either move the pin (bundles: vendor the release and re-run "
            "`pixi run verify-demux-model`; escpod: bump escpod_version AND the "
            "checksums, then `pixi run setup`; dorado: only with a charging "
            f"bundle trained on it), or record why not in {PINS.relative_to(REPO_ROOT)}."
        )
    if problems or drifted:
        return 1
    if unreachable:
        return 1 if args.strict else 0
    print("Every pin is current or deliberately held.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
