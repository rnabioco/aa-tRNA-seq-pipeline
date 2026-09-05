"""
Unit tests for run_config.py: the scaffolder and checker for a run's config pair.

Validation runs against the real repository (its vendored bundles and its
committed fixtures), so these also pin that the shipped test configs pass the
checker they document.
"""

from pathlib import Path

import pytest
import yaml

from run_config import (
    PlanError,
    load_base_config,
    load_pair,
    parse_spec,
    render_config,
    render_samples,
    validate_plan,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
LDX_RUN = str(REPO_ROOT / ".tests" / "fixtures" / "ldx-demux" / "run")
RAW_FA = str(REPO_ROOT / ".tests" / "fixtures" / "ldx-demux" / "collapsed.fa")
PLAIN = "GGCTTCTTCTTGCTCTTAGGAAAAAAAAAA"


@pytest.fixture(scope="module")
def base():
    return load_base_config(REPO_ROOT)


def plan(samples, **kw):
    p = {
        "name": "t",
        "output_dir": "results/t",
        "runs": [{"path": LDX_RUN, "samples": samples}],
        "adapters": None,
        "reference": {"mode": "build", "raw_fasta": RAW_FA},
        "ldx_gpu": False,
    }
    p.update(kw)
    return p


def errors(problems):
    return [m for lv, m in problems if lv == "error"]


class TestParseSpec:
    @pytest.mark.parametrize(
        "spec,expected",
        [
            ("ldx01", {"ldx": "ldx01"}),
            ("ldx01+fdx01", {"ldx": "ldx01", "fdx": "fdx01"}),
            ("fdx02+ldx04", {"ldx": "ldx04", "fdx": "fdx02"}),
            ("ldx01/edx01", {"ldx": "ldx01", "edx": "edx01"}),
            ("ldx01+fdx01/edx07", {"ldx": "ldx01", "fdx": "fdx01", "edx": "edx07"}),
            ("wdx:barcode03", {"wdx": "barcode03"}),
            ("-", {}),
            ("", {}),
        ],
    )
    def test_shapes(self, spec, expected):
        assert parse_spec(spec) == expected

    @pytest.mark.parametrize("spec", ["nbc01", "ldx01+ldx02", "ldx01+", "wdx:", "barcode03"])
    def test_rejects(self, spec):
        with pytest.raises(PlanError):
            parse_spec(spec)


class TestValidatePlan:
    def test_clean_dual_index_plan(self, base):
        p = plan(
            {"a": {"ldx": "ldx01", "fdx": "fdx01"}, "b": {"ldx": "ldx04", "fdx": "fdx02"}},
            adapters={"five_prime": "CCTAAGAGCAAGAAGAAGCCTGG", "three_prime": PLAIN},
        )
        assert errors(validate_plan(p, base, REPO_ROOT)) == []

    def test_unknown_codes_are_named_with_the_bundle(self, base):
        p = plan({"a": {"ldx": "ldx17"}, "b": {"ldx": "ldx01", "fdx": "fdx09"}})
        errs = errors(validate_plan(p, base, REPO_ROOT))
        assert any("ldx17" in e and "barcode_crf_ldx16" in e for e in errs)
        assert any("fdx09" in e and "barcode_crf_fdx4" in e for e in errs)

    def test_fdx_needs_ldx(self, base):
        p = plan({"a": {"fdx": "fdx01"}})
        assert any("fdx" in e and "ldx" in e for e in errors(validate_plan(p, base, REPO_ROOT)))

    def test_mixed_axis_sets_on_one_ldx_code(self, base):
        p = plan({"a": {"ldx": "ldx01"}, "b": {"ldx": "ldx01", "fdx": "fdx01"}})
        assert any("all or none" in e for e in errors(validate_plan(p, base, REPO_ROOT)))

    def test_identical_tuples(self, base):
        p = plan({"a": {"ldx": "ldx01", "fdx": "fdx01"}, "b": {"ldx": "ldx01", "fdx": "fdx01"}})
        assert any("identical" in e for e in errors(validate_plan(p, base, REPO_ROOT)))

    def test_edx_must_be_a_declared_adapter(self, base):
        p = plan({"a": {"ldx": "ldx01", "edx": "edx09"}})
        assert any("edx09" in e for e in errors(validate_plan(p, base, REPO_ROOT)))
        p = plan(
            {"a": {"ldx": "ldx01", "edx": "edx07"}},
            adapters={"three_prime": [{"name": "edx07", "seq": "GGCTTCTTCTTGCTCTTAGGAAGGCTGACAGTCTCAAAAAAAAAA"}]},
        )
        assert errors(validate_plan(p, base, REPO_ROOT)) == []

    def test_three_prime_adapter_rules(self, base):
        p = plan(
            {"a": {"ldx": "ldx01"}},
            adapters={"three_prime": [{"name": "x", "seq": "GGCAAAA"}, {"name": "y", "seq": "GGCAAAAAAAA"}]},
        )
        assert any("length" in e for e in errors(validate_plan(p, base, REPO_ROOT)))
        p = plan({"a": {"ldx": "ldx01"}}, adapters={"three_prime": "TTCTTCTTGCTCTT"})
        assert any("GGC" in e for e in errors(validate_plan(p, base, REPO_ROOT)))

    def test_missing_run_and_reference(self, base):
        p = plan({"a": {"ldx": "ldx01"}}, reference={"mode": "build", "raw_fasta": "/nope/x.fa"})
        p["runs"][0]["path"] = "/nope/run"
        errs = errors(validate_plan(p, base, REPO_ROOT))
        assert any("run directory" in e for e in errs)
        assert any("raw reference" in e for e in errs)

    def test_unbarcoded_sample_cannot_share_a_run(self, base):
        p = plan({"a": {}, "b": {"ldx": "ldx01"}})
        assert any("unbarcoded" in e for e in errors(validate_plan(p, base, REPO_ROOT)))

    def test_wdx_and_escpod_do_not_mix(self, base):
        p = plan({"a": {"wdx": "barcode03"}, "b": {"ldx": "ldx01"}})
        assert any("mutually exclusive" in e for e in errors(validate_plan(p, base, REPO_ROOT)))


class TestRender:
    def test_samples_yaml_round_trips(self):
        p = plan({"a": {"ldx": "ldx01", "fdx": "fdx01", "edx": "edx07"}, "u": {"ldx": "ldx02"}})
        text, ext = render_samples(p)
        assert ext == "yml"
        data = yaml.safe_load(text)
        (run,) = data["runs"]
        assert run["path"] == LDX_RUN
        assert run["samples"]["a"] == {"ldx": "ldx01", "fdx": "fdx01", "edx": "edx07"}
        assert run["samples"]["u"] == {"ldx": "ldx02"}

    def test_unbarcoded_plan_becomes_tsv(self):
        p = plan({"only": {}})
        text, ext = render_samples(p)
        assert ext == "tsv"
        assert text == f"only\t{LDX_RUN}\n"

    def test_config_enables_only_the_axes_used(self):
        p = plan({"a": {"ldx": "ldx01", "fdx": "fdx01"}}, adapters={"three_prime": PLAIN})
        cfg = yaml.safe_load(render_config(p, "config/samples-t.yml"))
        assert cfg["ldx"]["enabled"] is True
        assert cfg["fdx"]["enabled"] is True
        assert "edx" not in cfg and "warpdemux" not in cfg
        assert cfg["adapters"]["three_prime"] == PLAIN
        assert cfg["reference"] == {"mode": "build", "raw_fasta": RAW_FA}
        assert cfg["samples"] == "config/samples-t.yml"
        assert cfg["output_directory"] == "results/t"

    def test_config_for_unbarcoded_tsv_plan_has_no_demux_block(self):
        p = plan({"only": {}}, reference=None)
        cfg = yaml.safe_load(render_config(p, "config/samples-t.tsv"))
        assert set(cfg) == {"samples", "output_directory"}


class TestLoadPair:
    @pytest.mark.parametrize("name", ["config-ldx-test.yml", "config-fdx-test.yml", "config-test.yml"])
    def test_shipped_test_configs_pass(self, name):
        cfg = REPO_ROOT / "config" / name
        p, merged, flags = load_pair(cfg, REPO_ROOT)
        problems = flags + validate_plan(p, merged, REPO_ROOT)
        assert errors(problems) == [], problems

    def test_fdx_samples_without_fdx_enabled_is_flagged(self, temp_dir):
        samples = temp_dir / "s.yml"
        samples.write_text(
            "runs:\n  - path: %s\n    samples:\n      a: {ldx: ldx01, fdx: fdx01}\n" % LDX_RUN
        )
        cfg = temp_dir / "config-x.yml"
        cfg.write_text(
            f"samples: {samples}\noutput_directory: out\nldx:\n  enabled: true\n"
            f"reference:\n  mode: build\n  raw_fasta: {RAW_FA}\n"
        )
        p, merged, flags = load_pair(cfg, REPO_ROOT)
        assert any("fdx.enabled" in m for _, m in flags)
