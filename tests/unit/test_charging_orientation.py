"""
Tests for `get_charging_orientation_arg` in workflow/rules/common.smk.

`escpod classify` decides the move-table frame from the data by default, which
needs >= 50 informative reads and is a hard ERROR below that -- so on a heavily
multiplexed run a sparse sample (~20 anchored reads at 224 samples per flow
cell) fails, and one such sample fails `rule all` for the whole corpus.

The knob exists to force the frame for those. Two properties matter:

  * `auto` must emit NOTHING, so a run that never sets it produces exactly the
    command line it produced before the knob existed.
  * a bad value must stop the DAG, not 448 individual jobs -- a forced frame
    that is wrong does not error, it silently mis-anchors every feature.

common.smk cannot be imported (Snakemake directives, and its module level
touches `config` and the filesystem), so the loader lifts out the one function
with `ast`. Same approach as tests/unit/test_barcode_tuples.py.
"""

import ast
import re
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).parent.parent.parent
COMMON_SMK = REPO_ROOT / "workflow" / "rules" / "common.smk"

_DIRECTIVE = re.compile(
    r"^(rule|checkpoint|wildcard_constraints|onstart|onsuccess|onerror"
    r"|ruleorder|localrules|include|configfile|container)\b"
)


def _load(func_name, config):
    """Lift one function out of common.smk and run it against `config`."""
    lines = COMMON_SMK.read_text().splitlines(keepends=True)
    cut = next(
        (i for i, line in enumerate(lines) if _DIRECTIVE.match(line)), len(lines)
    )
    tree = ast.parse("".join(lines[:cut]))
    wanted = {func_name, "_ORIENTATIONS"}
    body = [
        n
        for n in tree.body
        if (isinstance(n, ast.FunctionDef) and n.name in wanted)
        or (
            isinstance(n, ast.Assign)
            and any(
                isinstance(t, ast.Name) and t.id in wanted for t in n.targets
            )
        )
    ]
    namespace = {"sys": sys, "config": config}
    exec(compile(ast.Module(body=body, type_ignores=[]), str(COMMON_SMK), "exec"), namespace)
    return namespace[func_name]


def orientation_arg(value=...):
    charging = {} if value is ... else {"orientation": value}
    return _load("get_charging_orientation_arg", {"charging": charging})()


def fallback(orientation=..., fallback_value=...):
    charging = {}
    if orientation is not ...:
        charging["orientation"] = orientation
    if fallback_value is not ...:
        charging["orientation_fallback"] = fallback_value
    return _load("get_charging_orientation_fallback", {"charging": charging})()


class TestDefaultIsSilent:
    def test_unset_emits_nothing(self):
        assert orientation_arg() == ""

    def test_explicit_auto_emits_nothing(self):
        """`auto` IS escpod's default; restating it on the command line would
        say nothing and would change the shell string for every existing run."""
        assert orientation_arg("auto") == ""

    def test_missing_charging_block_is_fine(self):
        assert _load("get_charging_orientation_arg", {})() == ""


class TestForcedFrames:
    @pytest.mark.parametrize("value", ["time", "reversed"])
    def test_forced_frame_is_passed_through(self, value):
        assert orientation_arg(value) == f"--orientation {value}"


class TestRejectsBadValues:
    @pytest.mark.parametrize(
        "value", ["Auto", "AUTO", "forward", "reverse", "time ", "", None, True, 1]
    )
    def test_bad_value_stops_the_dag(self, value):
        with pytest.raises(SystemExit, match="charging.orientation"):
            orientation_arg(value)

    def test_message_names_the_valid_values(self):
        with pytest.raises(SystemExit) as exc:
            orientation_arg("sideways")
        message = str(exc.value)
        for value in ("auto", "time", "reversed"):
            assert value in message

    def test_reverse_is_not_silently_accepted_as_reversed(self):
        """The near-miss a person actually types."""
        with pytest.raises(SystemExit):
            orientation_arg("reverse")


class TestFallbackResolution:
    """
    `orientation_fallback` supplies a frame for a sample too thin for `auto`,
    instead of failing it. It is a fallback rather than a forced default so that
    detection keeps running on every deep sample -- the frame's two inputs are
    already pinned (`base_calling_model`) and enforced
    (`charging.basecaller_check: error`), so the consensus escpod computes is a
    free tripwire for the case where those pins stop being true.
    """

    def test_default_is_reversed(self):
        """The measured frame: 8,057,646 votes across 446 samples, no dissent."""
        assert fallback() == "reversed"

    def test_explicit_value_is_used(self):
        assert fallback(fallback_value="time") == "time"

    @pytest.mark.parametrize("off", ["none", None, False])
    def test_can_be_switched_off(self, off):
        """Restores v0.7.2: an underpowered sample fails the run."""
        assert fallback(fallback_value=off) == ""

    @pytest.mark.parametrize("forced", ["time", "reversed"])
    def test_unused_when_a_frame_is_already_forced(self, forced):
        """Nothing to fall back from, so the retry branch stays off."""
        assert fallback(orientation=forced, fallback_value="reversed") == ""

    def test_applies_when_orientation_is_explicitly_auto(self):
        assert fallback(orientation="auto", fallback_value="time") == "time"

    @pytest.mark.parametrize("bad", ["auto", "Reversed", "reverse", "", 1])
    def test_bad_value_stops_the_dag(self, bad):
        with pytest.raises(SystemExit, match="orientation_fallback"):
            fallback(fallback_value=bad)

    def test_auto_is_not_a_valid_fallback(self):
        """`auto` is what failed; naming it as its own fallback is a loop."""
        with pytest.raises(SystemExit):
            fallback(fallback_value="auto")
