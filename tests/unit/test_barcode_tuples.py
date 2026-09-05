"""
Tests for `_check_barcode_tuples` in workflow/rules/common.smk.

The uniqueness check refuses two samples on one run that would receive
identical reads. It has to agree with the sibling check in
`workflow/scripts/run_config.py` (`pixi run check-run-config`), or the
validator passes configs the pipeline then refuses at parse time -- which is
exactly what shipped in v0.7.0: the tuple was built from (ldx, fdx) and left
`edx` out, so any design fanning one LDX code across several EDX adapters died
with "have identical barcode assignments" (#161).

common.smk cannot be imported: it mixes Snakemake directives with Python, and
its module level touches `config` and the filesystem. `_check_barcode_tuples`
is pure, though, so the loader below lifts that one function out with `ast` and
executes it alone. That keeps the test on the shipped source -- editing
common.smk changes what runs here -- without needing a Snakemake workflow.
"""

import ast
import re
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).parent.parent.parent
COMMON_SMK = REPO_ROOT / "workflow" / "rules" / "common.smk"

# Snakemake directives are not valid Python, so only the prefix before the
# first one is parsed. Every function definition lives above it.
_DIRECTIVE = re.compile(
    r"^(rule|checkpoint|wildcard_constraints|onstart|onsuccess|onerror"
    r"|ruleorder|localrules|include|configfile|container)\b"
)


def _load(func_name):
    """Return one top-level function from common.smk, executed in isolation."""
    lines = COMMON_SMK.read_text().splitlines(keepends=True)
    cut = next(
        (i for i, line in enumerate(lines) if _DIRECTIVE.match(line)), len(lines)
    )
    tree = ast.parse("".join(lines[:cut]))
    node = next(
        n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == func_name
    )
    namespace = {"sys": sys}
    exec(
        compile(ast.Module(body=[node], type_ignores=[]), str(COMMON_SMK), "exec"),
        namespace,
    )
    return namespace[func_name]


check_barcode_tuples = _load("_check_barcode_tuples")


def sample(barcode, fdx=None, edx=None, run_id="run1"):
    """A samples-dict entry shaped like parse_samples_yaml builds it."""
    return {"barcode": barcode, "fdx": fdx, "edx": edx, "run_id": run_id}


class TestAcceptsDistinctTuples:
    """Designs that give each sample a different set of reads."""

    def test_one_ldx_across_two_edx(self):
        # The v0.7.0 regression: legal, and was refused.
        samples = {
            "a_edx01": sample("ldx01", edx="edx01"),
            "a_edx02": sample("ldx01", edx="edx02"),
        }
        check_barcode_tuples(samples, "samples.yml")

    def test_one_ldx_fanned_across_seven_edx(self):
        # The shape of the GlnRS timecourse that surfaced this.
        samples = {f"s{i}": sample("ldx01", edx=f"edx{i:02d}") for i in range(1, 8)}
        check_barcode_tuples(samples, "samples.yml")

    def test_one_ldx_across_two_fdx(self):
        samples = {
            "a_fdx01": sample("ldx01", fdx="fdx01"),
            "a_fdx02": sample("ldx01", fdx="fdx02"),
        }
        check_barcode_tuples(samples, "samples.yml")

    def test_same_tuple_on_different_runs(self):
        samples = {
            "a": sample("ldx01", edx="edx01", run_id="run1"),
            "b": sample("ldx01", edx="edx01", run_id="run2"),
        }
        check_barcode_tuples(samples, "samples.yml")

    def test_distinct_ldx_codes(self):
        samples = {"a": sample("ldx01"), "b": sample("ldx02")}
        check_barcode_tuples(samples, "samples.yml")

    def test_unbarcoded_samples_are_skipped(self):
        samples = {"a": sample(None), "b": sample(None)}
        check_barcode_tuples(samples, "samples.yml")


class TestRefusesOverlap:
    """Designs where two samples would receive the same reads."""

    def test_identical_ldx_alone(self):
        samples = {"a": sample("ldx01"), "b": sample("ldx01")}
        with pytest.raises(SystemExit, match="identical"):
            check_barcode_tuples(samples, "samples.yml")

    def test_identical_ldx_and_edx(self):
        samples = {
            "a": sample("ldx01", edx="edx01"),
            "b": sample("ldx01", edx="edx01"),
        }
        with pytest.raises(SystemExit, match="identical"):
            check_barcode_tuples(samples, "samples.yml")

    def test_identical_ldx_and_fdx(self):
        samples = {
            "a": sample("ldx01", fdx="fdx01"),
            "b": sample("ldx01", fdx="fdx01"),
        }
        with pytest.raises(SystemExit, match="identical"):
            check_barcode_tuples(samples, "samples.yml")

    def test_bare_ldx_would_swallow_a_dual_index_sample(self):
        # Caught by the earlier all-or-none `fdx:` check, not the tuple check.
        samples = {
            "bare": sample("ldx01"),
            "dual": sample("ldx01", fdx="fdx01"),
        }
        with pytest.raises(SystemExit, match="fdx"):
            check_barcode_tuples(samples, "samples.yml")


class TestAgreesWithRunConfig:
    """
    The two checks must refuse the same designs. run_config.py keys its tuple
    on (fdx, edx) within an LDX group; common.smk keys on (ldx, fdx, edx)
    across the run. Same axes, so the verdicts have to match -- a divergence
    here is what let #161 reach a release.
    """

    @pytest.mark.parametrize(
        "codes",
        [
            [{"ldx": "ldx01", "edx": "edx01"}, {"ldx": "ldx01", "edx": "edx02"}],
            [{"ldx": "ldx01", "fdx": "fdx01"}, {"ldx": "ldx01", "fdx": "fdx02"}],
            [{"ldx": "ldx01", "edx": "edx01"}, {"ldx": "ldx01", "edx": "edx01"}],
            [{"ldx": "ldx01"}, {"ldx": "ldx01"}],
        ],
    )
    def test_same_verdict(self, codes):
        names = [f"s{i}" for i in range(len(codes))]

        smk_samples = {
            name: sample(c["ldx"], fdx=c.get("fdx"), edx=c.get("edx"))
            for name, c in zip(names, codes, strict=True)
        }
        try:
            check_barcode_tuples(smk_samples, "samples.yml")
            smk_refused = False
        except SystemExit:
            smk_refused = True

        # The same grouping run_config.py applies, on the same inputs.
        by_ldx = {}
        for name, c in zip(names, codes, strict=True):
            by_ldx.setdefault(c["ldx"], []).append((name, c))
        cfg_refused = False
        for _, group in by_ldx.items():
            tuples = [(c.get("fdx"), c.get("edx")) for _, c in group]
            if len(set(tuples)) != len(tuples):
                cfg_refused = True

        assert smk_refused == cfg_refused
