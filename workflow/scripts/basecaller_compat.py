"""Check this run's basecaller against the one a charging bundle was trained on.

Since escapepod-models#106 a charging bundle declares a `basecaller` block —
`{model, model_sha256, dorado_version}` — and escpod states it at load without
enforcing it. This is the enforcing half.

It matters because of what the model actually reads. The charging feature set is
`mean + z-scored k-mer residual`, and the EXPECTED level is predicted from the
read's own basecall, so a charging model substantially detects *how the
basecaller fails* at the aminoacyl adduct. Change the basecaller and the
dominant feature means something else. Upstream measured it: the same reads
called two ways lose ~0.0097 AUROC and 3.0-3.2 pp of TPR, gain ~0.4 pp of FPR
(so no threshold recovers it), and flip **3.9% of per-read calls** — while the
aggregate charged fraction moves 0.04 pp. The one statistic anyone would check
when changing basecaller reads "no change" while one read in 26 answers
differently. That is the failure this exists to make visible.

Two rules, deliberately of different strength:

  MODEL IDENTITY is the rule. The bundle names the basecalling model it was
  trained against and the declaration's own note says inference data should be
  basecalled with *that model*. A different one is a domain shift, so a
  mismatch is an error by default.

  DORADO VERSION is provenance, and secondary — the same model weights run by a
  later dorado are still the same weights. But a MAJOR version difference is a
  different implementation of the basecaller, which can move basecalls with the
  weights unchanged, so that is worth saying out loud. It warns; it never
  blocks, because blocking on it would strand every run whose bundle was built
  by an older dorado than the one currently pinned.

A bundle that declares nothing (anything published before #106) produces no
findings. "Cannot tell" is not "invalid" — the same rule
`bundle_barcode_names()` follows next door.

`model_digest()` reproduces upstream's hash EXACTLY — sha256 over every file in
the directory, name-sorted, each contributing its relative path then its bytes
(`escapepod-models/scripts/charging/verify_basecaller.py`). It reads the whole
model (~300 MB), so it is not run while the DAG is built; `pixi run
verify-basecaller` is where it lives.
"""

import hashlib
import re
from pathlib import Path

# `{model, model_sha256, dorado_version}` — upstream refuses to build a charging
# bundle without all three, so a block carrying only some of them is malformed
# rather than partial.
REQUIRED_KEYS = ("model", "model_sha256", "dorado_version")


def bundle_basecaller(bundle_dir):
    """The `basecaller` block a charging bundle declares, or None.

    None means the bundle predates escapepod-models#106 and simply does not say
    — not that it is basecaller-agnostic. Reads `metadata.json` only, so this is
    cheap enough for DAG construction.
    """
    meta = Path(bundle_dir) / "metadata.json"
    if not meta.is_file():
        return None
    import json

    with open(meta) as fh:
        block = json.load(fh).get("basecaller")
    if not isinstance(block, dict) or not all(block.get(k) for k in REQUIRED_KEYS):
        return None
    return block


def parse_version(version):
    """Leading numeric components of a version string, as a tuple of ints.

    `1.4.0+ba44a013` -> `(1, 4, 0)`. The build metadata after `+` identifies the
    commit, not the release, and comparing it would make every rebuild look like
    a different version. Returns () when nothing numeric leads, so an unparseable
    version compares equal to nothing and is reported rather than guessed at.
    """
    match = re.match(r"\s*v?(\d+(?:\.\d+)*)", str(version or ""))
    return tuple(int(p) for p in match.group(1).split(".")) if match else ()


def model_digest(model_dir):
    """(sha256, file_count) over a dorado model directory.

    Upstream's scheme, reproduced exactly so the two sides can be compared at
    all: name-sorted, each file contributing its POSIX relative path and then
    its bytes. Any deviation — skipping dotfiles, sorting differently — silently
    produces a hash that can never match.
    """
    model_dir = Path(model_dir)
    digest = hashlib.sha256()
    count = 0
    for path in sorted(p for p in model_dir.rglob("*") if p.is_file()):
        digest.update(path.relative_to(model_dir).as_posix().encode())
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                digest.update(chunk)
        count += 1
    return digest.hexdigest(), count


def check_basecaller(bundle_dir, basecall_model, dorado_version):
    """Compare a run's basecaller against a charging bundle's declaration.

    `basecall_model` is the model actually handed to dorado (a path or a name;
    only its final component is compared), and `dorado_version` the pinned
    binary version.

    Returns a list of `(level, message)` with level "error" or "warn", most
    severe first. Empty means compatible, or that the bundle does not say.
    """
    declared = bundle_basecaller(bundle_dir)
    if declared is None:
        return []

    findings = []
    want_model = declared["model"]
    got_model = Path(str(basecall_model)).name

    if got_model != want_model:
        findings.append(
            (
                "error",
                f"charging bundle {Path(bundle_dir).name} was trained on basecalls "
                f"from {want_model}, but this run basecalls with {got_model}.\n"
                "The charging feature set is mean + z-scored k-mer residual, and the "
                "expected level is predicted from the read's own basecall — so the "
                "model substantially detects HOW THE BASECALLER FAILS at the "
                "aminoacyl adduct, and a different one is a domain shift rather than "
                "a detail. Upstream measured ~0.0097 AUROC and 3.9% of per-read "
                "calls flipped, while the aggregate charged fraction moved 0.04 pp: "
                "the number you would check reads 'no change'.\n"
                "Point `base_calling_model` at the declared model, or take a "
                "charging bundle trained on this one. To run anyway — arm-to-arm "
                "contrasts survive a shared shift, absolute charged fractions do "
                "not — set `charging.basecaller_check: warn`.",
            )
        )

    want_dorado = parse_version(declared["dorado_version"])
    got_dorado = parse_version(dorado_version)
    if not want_dorado or not got_dorado:
        findings.append(
            (
                "warn",
                f"could not compare dorado versions (bundle declares "
                f"{declared['dorado_version']!r}, config pins {dorado_version!r}).",
            )
        )
    elif want_dorado[0] != got_dorado[0]:
        findings.append(
            (
                "warn",
                f"dorado major version differs: the bundle was built with "
                f"{declared['dorado_version']}, this run pins {dorado_version}.\n"
                "Model identity is the rule that governs charging calls, so this "
                "never blocks a run on its own. But a major bump is a different "
                "implementation of the basecaller and can move basecalls with the "
                "weights unchanged — worth revalidating before trusting absolute "
                "charged fractions across it.",
            )
        )
    return findings
