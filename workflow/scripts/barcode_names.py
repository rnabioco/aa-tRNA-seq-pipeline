"""Translate between the barcode names a demux model emits and the names this
project uses.

Every CRF bundle emits whatever `metadata.json` calls its references, and that
vocabulary is upstream's, not ours. Two panels are in play and they run in
OPPOSITE directions, which is why this cannot be one conditional rename:

    panel      samples YAML   model emits   project-facing
    LDX        nbc01          nbc01         ldx01      <- upstream name is configured
    WDX-CRF    barcode03      bc03          barcode03  <- OUR name is configured

So a sample has two names, and which of them is "upstream" depends on the panel.
Callers should say which one they want:

    emitted_to_label   what to call a barcode the model just emitted
    label_to_emitted   what to look for in the model's own output

Both directions are pure renames of a numbered code (documented in
config/README.md and resources/models/demux/README.md), so nothing is lost.

The renames are applied with `fullmatch` against a fixed prefix table, so no
caller has to ask which backend is enabled.

Scope note: `label_to_emitted` is meaningful only on the escpod path. WarpDemuX
names (`barcode07`) are indistinguishable from this project's WDX-CRF labels --
that is the point, the CRF panel is the same physical kit -- so on a WarpDemuX
run, where reads are routed by `pod5 filter` rather than by escpod's own output
filenames, do not call it. `emitted_to_label` is safe everywhere: it is the
identity on every WarpDemuX name.
"""

import json
import re
from pathlib import Path

# emitted prefix -> project-facing prefix, for a name of the form <prefix><digits>.
# Order matters only in that the emitted prefixes must be checked longest-first,
# so `nbc01` is not read as `n` + `bc01`; `_EMITTED_TO_LABEL` is built sorted.
_PREFIXES = {
    "nbc": "ldx",
    "bc": "barcode",
}
_EMITTED_TO_LABEL = sorted(_PREFIXES.items(), key=lambda kv: -len(kv[0]))
_LABEL_TO_EMITTED = sorted(
    ((label, emitted) for emitted, label in _PREFIXES.items()),
    key=lambda kv: -len(kv[0]),
)


def _rename(name, table):
    """Rewrite <prefix><digits> using the first matching prefix in `table`."""
    if not name:
        return name
    for src, dst in table:
        match = re.fullmatch(rf"{re.escape(src)}(\d+)", name)
        if match:
            return f"{dst}{match.group(1)}"
    return name


def emitted_to_label(name):
    """The project-facing name for a barcode as the model emits it.

    `nbc01` -> `ldx01`, `bc03` -> `barcode03`. Anything else -- including
    `unclassified`, an already-project-facing name, and None -- is returned
    unchanged.
    """
    return _rename(name, _EMITTED_TO_LABEL)


def label_to_emitted(name):
    """The name the model emits for a project-facing barcode.

    The inverse of `emitted_to_label`, and the one that matters at run time:
    `escpod demux` writes `barcode_<emitted>.pod5`, so a sample configured as
    `barcode03` has to be looked up as `bc03`.

    A name already in the emitted vocabulary is returned unchanged, so this is
    safe on a samples file written either way. See the module docstring for why
    it must not be used on a WarpDemuX run.
    """
    if not name:
        return name
    for emitted, _ in _EMITTED_TO_LABEL:
        if re.fullmatch(rf"{re.escape(emitted)}(\d+)", name):
            return name  # already emitted-side
    return _rename(name, _LABEL_TO_EMITTED)


def bundle_barcode_names(bundle_dir):
    """The barcode names a CRF bundle declares, as a set.

    Reads `metadata.json` only -- no ONNX is loaded -- so this is cheap enough to
    call during DAG construction. Returns an empty set when the path is not a
    CRF bundle directory (a GBM/DTW model is a single JSON file and carries its
    barcode set differently), letting callers treat "cannot tell" as "do not
    validate" rather than as "no barcodes are valid".
    """
    meta = Path(bundle_dir) / "metadata.json"
    if not meta.is_file():
        return set()
    with open(meta) as f:
        data = json.load(f)
    return {bc["name"] for bc in data.get("barcodes", []) if "name" in bc}
