"""Translate between the barcode names a demux model emits and the names this
project uses.

Every CRF bundle emits whatever `metadata.json` calls its references, and only
ONE panel still needs translating:

    panel      samples YAML   model emits   project-facing
    LDX        ldx01          ldx01         ldx01      <- no rename at all
    WDX-CRF    barcode03      bc03          barcode03  <- OUR name is configured

The LDX row used to read `nbc01` in the first two columns, because the shipped
bundle was `barcode_crf_nbc16_rna004` and upstream called its references `nbcNN`.
That family is retired here -- the pipeline runs `barcode_crf_ldx16_rna004`,
whose class names ARE `ldx01`..`ldx16` -- so translating the LDX panel is not
merely unnecessary, it is wrong: rewriting a configured `ldx01` to `nbc01`
produces a name that appears nowhere in an ldx16 `classifications.csv`, and
every sample on the run then fails with "No reads were assigned". `nbc` is
therefore absent from the prefix table on purpose. Do not add it back to support
an nbc bundle; point `ldx.model` at an ldx one.

Callers should say which direction they want:

    emitted_to_label   what to call a barcode the model just emitted
    label_to_emitted   what to look for in the model's own output
    resolve_to_bundle  the same question asked of a SPECIFIC bundle

The first two are pure renames of a numbered code against a fixed prefix table
(documented in config/README.md and resources/models/demux/README.md), applied
with `fullmatch` so no caller has to ask which backend is enabled. They are
guesses about a vocabulary, though, and a wrong guess is silent until the run
dies hours later -- which is exactly how the nbc rename survived the switch to
ldx16. `resolve_to_bundle` asks the bundle instead, and is what the workflow
uses; prefer it wherever a bundle path is in hand.

Scope note: `label_to_emitted` and `resolve_to_bundle` are meaningful only on
the escpod path. WarpDemuX names (`barcode07`) are indistinguishable from this
project's WDX-CRF labels -- that is the point, the CRF panel is the same
physical kit -- so on a WarpDemuX run, where reads are routed by `pod5 filter`
rather than by escpod's own output filenames, do not call them.
`emitted_to_label` is safe everywhere: it is the identity on every WarpDemuX
name.
"""

import json
import re
from pathlib import Path

# emitted prefix -> project-facing prefix, for a name of the form <prefix><digits>.
# One entry, and the LDX panel is deliberately not in it: `ldx16` emits `ldxNN`
# already. Tables are still built sorted longest-first, because a second entry
# whose prefix ends in an existing one would otherwise be mis-split.
_PREFIXES = {
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

    `bc03` -> `barcode03`. Anything else -- including `ldx01`, `unclassified`,
    an already-project-facing name, and None -- is returned unchanged.
    """
    return _rename(name, _EMITTED_TO_LABEL)


def label_to_emitted(name):
    """The name the model emits for a project-facing barcode.

    The inverse of `emitted_to_label`: `escpod demux` writes `bc03` into
    `classifications.csv`, so a sample configured as `barcode03` has to be
    looked up as `bc03`. An LDX name is returned unchanged, because the shipped
    ldx16 bundle emits it verbatim.

    A name already in the emitted vocabulary is returned unchanged, so this is
    safe on a samples file written either way. See the module docstring for why
    it must not be used on a WarpDemuX run, and prefer `resolve_to_bundle` when
    the bundle is known.
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


class BarcodeVocabularyError(ValueError):
    """A configured barcode is not in the bundle's declared reference set."""


def resolve_to_bundle(name, bundle_dir):
    """The name to match in `bundle_dir`'s `classifications.csv` for `name`.

    Asks the bundle what it calls its references rather than guessing from a
    prefix, and is the form the workflow uses. The prefix tables above still do
    the work -- this only decides which of a name's spellings the bundle
    actually recognises, and says so out loud when none of them does.

    Resolution order is "believe the config first": a configured name the bundle
    declares is returned as-is, and only then are the two renames tried. That
    matters because the panels are not disjoint in shape -- `bc03` is a real
    ldx-era spelling and a real WDX4 one -- so a name that is already correct
    must never be rewritten into a different panel's vocabulary.

    Raises BarcodeVocabularyError when nothing resolves, naming the bundle's own
    vocabulary. This is the whole point of reading the bundle: the alternative is
    a name that appears nowhere in the model's output, which surfaces as "No
    reads were assigned" only AFTER the demux pass and the whole-run basecall
    have been paid for. Raised during DAG construction instead.

    A path that is not a CRF bundle directory declares nothing, and then this
    falls back to `label_to_emitted` without validating -- "cannot tell" is not
    "nothing is valid". See the module docstring for why neither may be used on
    a WarpDemuX run.
    """
    declared = bundle_barcode_names(bundle_dir)
    if not declared:
        return label_to_emitted(name)
    for candidate in (name, label_to_emitted(name), emitted_to_label(name)):
        if candidate in declared:
            return candidate
    raise BarcodeVocabularyError(
        f"barcode {name!r} is not declared by the demux bundle "
        f"{Path(bundle_dir).name}.\n"
        f"That bundle emits: {', '.join(sorted(declared))}.\n"
        "Fix the samples file to use one of those names, or point `ldx.model` "
        "at the bundle this run was demultiplexed with. (Names of the form "
        "`nbcNN` are the retired nbc16 panel: its successor is ldx16, whose "
        "references are already called `ldxNN`.)"
    )
