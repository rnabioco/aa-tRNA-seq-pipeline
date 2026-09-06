#!/usr/bin/env python3
"""Verify a vendored escapepod demux bundle against upstream's own attestation.

Deliberately does NOT trust any checksum file we could have written, for the
files upstream can attest. The two hashes it checks there are the ones the
RELEASE declares inside its own metadata:

    metadata.json   boundary.sha256   -> the boundary detector ONNX
    provenance.json sha256            -> the CRF ONNX

A bundle that has been edited locally fails here. That is the point: patching a
vendored bundle (as #107/#108 did to nbc16, adding boundary.margin and
clamp_max_shift) left it failing its own SHA256SUMS for eleven days, and then
looked like evidence that a later upstream release had removed those keys.

A bundle's own SHA256SUMS.txt is checked too when the release shipped one.

THE SIDECARS ARE THE GAP, and they are covered separately. `metadata.json` is
where the ONNX hashes live, so nothing upstream can attest it without being
self-referential; `provenance.json` is the same. An edit to either — adding
window rules, say — passes every check above. Upstream cannot close this and
compute nodes have no route to GitHub, so the only offline control is a digest
recorded when the bundle was vendored, in sidecars.sha256. That buys
tamper-evidence, not authenticity: it proves a bundle has not changed since we
vendored it, not that it is upstream's. Vendor from a verified source, then
record.

    scripts/verify-demux-model.py --record            # bundles this repo vendors
    scripts/verify-demux-model.py --record --local    # bundles only YOU vendor

The two-file split matters. `sidecars.sha256` is tracked and describes the
bundles in git. A deployment that vendors extra bundles by hand -- ldx32, say,
which this repo does not ship -- records them in `sidecars.local.sha256`, which
is gitignored. Without that split, pulling an updated tracked manifest would
drop the local bundles' digests and fail them closed at the worst moment.

Until escapepod-models#127/#128 (2026-09-06) no upstream CRF bundle declared
boundary.margin or clamp_max_shift, and refusing any bundle that did was the
cheap way to catch a local edit. Bundles now declare them legitimately —
barcode_crf_ldx32_rna004@v0.2.2 is the first — so that heuristic is retired and
the sidecar digest is what replaces it.
"""

import hashlib
import json
import sys
from pathlib import Path

DEMUX_DIR = Path("resources/models/demux")
RECORD_NAME = "sidecars.sha256"
LOCAL_RECORD_NAME = "sidecars.local.sha256"
SIDECARS = ("metadata.json", "provenance.json")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_record(path):
    if not path.is_file():
        return {}
    recorded = {}
    for line in path.read_text().split("\n"):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        want, name = line.split(None, 1)
        recorded[name.strip()] = want
    return recorded


def load_record(demux_dir):
    """Tracked digests, extended by the gitignored local overlay."""
    recorded = parse_record(demux_dir / RECORD_NAME)
    recorded.update(parse_record(demux_dir / LOCAL_RECORD_NAME))
    return recorded


def check(label, path, expected, failures):
    if expected is None:
        print(f"  {label}: no hash declared upstream, skipped")
        return
    if not path.is_file():
        print(f"  {label}: MISSING {path.name}")
        failures.append(label)
        return
    actual = sha256(path)
    ok = actual == expected
    print(f"  {label}: {'OK' if ok else 'FAILED'}  {path.name}")
    if not ok:
        print(f"      declared {expected}\n      actual   {actual}")
        failures.append(label)


def verify(bundle):
    bundle = Path(bundle)
    print(f"{bundle}")
    if not bundle.is_dir():
        print("  not a directory")
        return False

    meta = json.loads((bundle / "metadata.json").read_text())
    prov = json.loads((bundle / "provenance.json").read_text())
    failures = []

    boundary = meta.get("boundary", {})
    check("boundary onnx (metadata.boundary.sha256)",
          bundle / boundary.get("onnx", "adapter_rna004.onnx"),
          boundary.get("sha256"), failures)
    check("crf onnx (provenance.sha256)",
          bundle / meta.get("onnx", ""), prov.get("sha256"), failures)

    sums = bundle / "SHA256SUMS.txt"
    if sums.is_file():
        for line in sums.read_text().split("\n"):
            if not line.strip():
                continue
            want, name = line.split()
            check(f"SHA256SUMS: {name}", bundle / name, want, failures)
    else:
        print("  SHA256SUMS.txt: not shipped by this release, skipped")

    # The sidecars, which no upstream hash covers. See the module docstring.
    recorded = load_record(bundle.parent)
    for name in SIDECARS:
        key = f"{bundle.name}/{name}"
        if key not in recorded:
            print(f"  sidecar {name}: NO DIGEST RECORDED in {RECORD_NAME} or "
                  f"{LOCAL_RECORD_NAME} — a local edit to it would pass "
                  f"unnoticed. Vendor from a verified source, then "
                  f"`scripts/verify-demux-model.py --record` (a bundle this "
                  f"repo ships) or `--record --local` (one you vendored).")
            failures.append(key)
            continue
        check(f"sidecar {name}", bundle / name, recorded[key], failures)

    # Window rules are legitimate upstream content as of escapepod-models
    # #127/#128; report them so a run's flags can be compared against them.
    declared = {k: boundary[k] for k in ("margin", "clamp_max_shift")
                if k in boundary}
    if declared:
        print(f"  window rules declared by the bundle: {declared}")

    return not failures


def record(targets, local=False):
    """Rewrite the sidecar digest manifest from what is on disk.

    With `local`, only bundles absent from the tracked manifest are written,
    into the gitignored overlay — so a deployment's own bundles are recorded
    without touching the file git owns.
    """
    by_dir = {}
    for t in targets:
        bundle = Path(t)
        if bundle.is_dir():
            by_dir.setdefault(bundle.parent, []).append(bundle)

    for demux_dir, bundles in sorted(by_dir.items()):
        if local:
            tracked = parse_record(demux_dir / RECORD_NAME)
            bundles = [b for b in bundles
                       if f"{b.name}/{SIDECARS[0]}" not in tracked]
        path = demux_dir / (LOCAL_RECORD_NAME if local else RECORD_NAME)
        lines = [
            "# Digests of the bundle sidecars, which no upstream hash covers.",
            "# Tamper-evidence only: recorded at vendoring time, so it proves a",
            "# bundle has not changed since, not that it is upstream's.",
            "# Regenerate with `scripts/verify-demux-model.py --record"
            + (" --local`." if local else "`."),
        ]
        for bundle in sorted(bundles, key=lambda b: b.name):
            for name in SIDECARS:
                f = bundle / name
                if f.is_file():
                    lines.append(f"{sha256(f)}  {bundle.name}/{name}")
        path.write_text("\n".join(lines) + "\n")
        print(f"recorded {len(lines) - 4} digests in {path}")


if __name__ == "__main__":
    flags = {"--record", "--local"}
    args = [a for a in sys.argv[1:] if a not in flags]
    targets = args or sorted(
        str(p) for p in DEMUX_DIR.iterdir() if p.is_dir()
    )
    if "--record" in sys.argv[1:]:
        record(targets, local="--local" in sys.argv[1:])
        sys.exit(0)
    # The list is deliberate: it verifies EVERY bundle before deciding,
    # so one bad bundle does not hide the others. `all(generator)`
    # short-circuits, which is what a ruff C419 autofix would give you.
    ok = all([verify(t) for t in targets])  # noqa: C419
    print("\nall bundles verify" if ok else "\nVERIFICATION FAILED")
    sys.exit(0 if ok else 1)
