"""
Unit tests for config_schema.py and workflow/schemas/config.schema.yaml

The schema closes every config section, so a misspelt or retired key is an error
instead of a silent fall-back to a default.
"""

import copy
import glob
from pathlib import Path

import pytest
import yaml
from config_schema import config_problems, validate_config

REPO = Path(__file__).resolve().parents[2]
BASE = yaml.safe_load((REPO / "config" / "config-base.yml").read_text())


def merged(override):
    """config-base.yml with a run's override applied, one level deep like Snakemake."""
    out = copy.deepcopy(BASE)
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k].update(v)
        else:
            out[k] = v
    return out


def valid(**extra):
    return {
        **copy.deepcopy(BASE),
        "samples": "s.tsv",
        "output_directory": "out",
        **extra,
    }


@pytest.mark.parametrize(
    "path", sorted(glob.glob(str(REPO / "config" / "config-*.yml"))), ids=Path
)
def test_every_committed_config_is_valid(path):
    if path.endswith("config-base.yml"):
        cfg = {**BASE, "samples": "s.tsv", "output_directory": "out"}
    else:
        cfg = merged(yaml.safe_load(Path(path).read_text()))
    assert config_problems(cfg) == []


def test_base_alone_lacks_only_the_run_keys():
    assert {p.split(": ", 1)[1] for p in config_problems(BASE)} == {
        "'samples' is a required property",
        "'output_directory' is a required property",
    }


def test_typo_names_the_nearest_key():
    cfg = valid()
    cfg["charging"]["ml_treshold"] = 220
    assert config_problems(cfg) == [
        "charging.ml_treshold: unknown key (did you mean 'ml_threshold'?)"
    ]


def test_unknown_top_level_key():
    assert config_problems(valid(dorado_verison="2.1.1")) == [
        "dorado_verison: unknown key (did you mean 'dorado_version'?)"
    ]


def test_retired_key_is_named():
    cfg = valid()
    cfg["opts"]["escpod_align"] = "--x"
    (problem,) = config_problems(cfg)
    assert problem.startswith("opts.escpod_align: renamed opts.eschalign")


def test_every_problem_is_reported_at_once():
    cfg = valid()
    cfg["ldx"]["min_crf_margn"] = 1.0
    cfg["charging"]["ml_threshold"] = 300
    assert len(config_problems(cfg)) == 2


def test_wrong_type_and_range():
    cfg = valid()
    cfg["charging"]["ml_threshold"] = 256
    cfg["ldx"]["gpu"] = "yes"
    problems = config_problems(cfg)
    assert any(p.startswith("charging.ml_threshold:") for p in problems)
    assert any(p.startswith("ldx.gpu:") for p in problems)


def test_odds_ratio_threshold_may_be_null_or_a_cut():
    assert config_problems(valid(odds_ratios={"ml_threshold": None})) == []
    assert config_problems(valid(odds_ratios={"ml_threshold": 180})) == []


def test_three_prime_accepts_string_or_named_list():
    cfg = valid()
    cfg["adapters"]["three_prime"] = "GGCTTCTTCTTGC"
    assert config_problems(cfg) == []
    cfg["adapters"]["three_prime"] = [{"name": "edx01"}]
    assert config_problems(cfg)


def test_cleanup_intermediates_bool_or_list():
    assert config_problems(valid(cleanup_intermediates=True)) == []
    assert config_problems(valid(cleanup_intermediates=["cascade"])) == []
    assert config_problems(valid(cleanup_intermediates="cascade"))


def test_validate_config_exits_listing_problems():
    cfg = valid(dorado_verison="x")
    with pytest.raises(SystemExit) as exc:
        validate_config(cfg)
    assert "dorado_verison" in str(exc.value)
