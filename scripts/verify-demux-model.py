#!/usr/bin/env python3
"""Verify a vendored escapepod demux bundle against upstream's own attestation.

Deliberately does NOT trust any checksum file we could have written. The two
hashes it checks are the ones the RELEASE declares inside its own metadata:

    metadata.json   boundary.sha256   -> the boundary detector ONNX
    provenance.json sha256            -> the CRF ONNX

A bundle that has been edited locally fails here. That is the point: patching a
vendored bundle (as #107/#108 did to nbc16, adding boundary.margin and
clamp_max_shift) left it failing its own SHA256SUMS for eleven days, and then
looked like evidence that a later upstream release had removed those keys.
Overrides belong in config-base.yml, not in the model.

A bundle's own SHA256SUMS.txt is checked too when the release shipped one.
"""

import hashlib
import json
import sys
from pathlib import Path


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


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

    # A locally added margin/clamp is the specific failure mode this guards.
    for key in ("margin", "clamp_max_shift"):
        if key in boundary:
            print(f"  boundary.{key} is declared in the bundle — no upstream "
                  f"CRF bundle declares it, so this copy has been edited. "
                  f"Set ldx.{key if key != 'margin' else 'boundary_margin'} "
                  f"in config-base.yml instead.")
            failures.append(f"boundary.{key}")

    return not failures


if __name__ == "__main__":
    targets = sys.argv[1:] or sorted(
        str(p) for p in Path("resources/models/demux").iterdir()
        if p.is_dir()
    )
    ok = all([verify(t) for t in targets])
    print("\nall bundles verify" if ok else "\nVERIFICATION FAILED")
    sys.exit(0 if ok else 1)
