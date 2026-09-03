"""
Unit tests for basecaller_compat.py

The charging bundle declares the basecaller it was trained against; these cover
the two rules that declaration implies, which are deliberately of different
strength — model identity is an error, dorado version is a warning.
"""

import json
from pathlib import Path

import pytest

from basecaller_compat import (
    bundle_basecaller,
    check_basecaller,
    model_digest,
    parse_version,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
CHARGING_BUNDLE = (
    REPO_ROOT
    / "resources"
    / "models"
    / "charging"
    / "charging_feature_nn_sup6_rna004@v0.1.0"
)

needs_bundle = pytest.mark.skipif(
    not CHARGING_BUNDLE.is_dir(), reason="charging bundle is not vendored in this repo"
)


def write_bundle(tmp_path, basecaller):
    """A minimal bundle directory declaring (or omitting) a basecaller block."""
    bundle = tmp_path / "bundle@v1.0.0"
    bundle.mkdir()
    meta = {"format": "charging", "classes": ["uncharged", "charged"]}
    if basecaller is not None:
        meta["basecaller"] = basecaller
    (bundle / "metadata.json").write_text(json.dumps(meta))
    return bundle


GOOD = {
    "model": "rna004_130bps_sup@v5.3.0",
    "model_sha256": "0" * 64,
    "dorado_version": "1.4.0+ba44a013",
}


class TestParseVersion:
    @pytest.mark.parametrize(
        "raw,expected",
        [
            ("1.4.0+ba44a013", (1, 4, 0)),
            ("2.1.1", (2, 1, 1)),
            ("v0.19.0", (0, 19, 0)),
            ("2.1.1+d66c17c", (2, 1, 1)),
            ("3", (3,)),
        ],
    )
    def test_numeric_lead_is_taken(self, raw, expected):
        """Build metadata after `+` names the commit, not the release."""
        assert parse_version(raw) == expected

    @pytest.mark.parametrize("raw", [None, "", "unknown", "latest"])
    def test_unparseable_is_empty(self, raw):
        """Empty compares equal to nothing, so it is reported, not guessed at."""
        assert parse_version(raw) == ()


class TestBundleBasecaller:
    def test_declared_block_is_returned(self, tmp_path):
        assert bundle_basecaller(write_bundle(tmp_path, GOOD)) == GOOD

    def test_absent_block_is_none(self, tmp_path):
        """A bundle published before escapepod-models#106 simply does not say."""
        assert bundle_basecaller(write_bundle(tmp_path, None)) is None

    @pytest.mark.parametrize("missing", ["model", "model_sha256", "dorado_version"])
    def test_partial_block_is_none(self, tmp_path, missing):
        """Upstream refuses to build one of these, so a partial block is malformed."""
        block = {k: v for k, v in GOOD.items() if k != missing}
        assert bundle_basecaller(write_bundle(tmp_path, block)) is None

    def test_missing_metadata_is_none(self, tmp_path):
        assert bundle_basecaller(tmp_path / "nope") is None

    @needs_bundle
    def test_vendored_bundle_declares_one(self):
        block = bundle_basecaller(CHARGING_BUNDLE)
        assert block is not None
        assert block["model"] == "rna004_sup@v6.0.0"


class TestCheckBasecaller:
    def test_matching_model_and_version_is_silent(self, tmp_path):
        bundle = write_bundle(tmp_path, GOOD)
        assert check_basecaller(bundle, "rna004_130bps_sup@v5.3.0", "1.4.0") == []

    def test_model_is_compared_by_name_not_path(self, tmp_path):
        """`base_calling_model` is a path; only its final component names a model."""
        bundle = write_bundle(tmp_path, GOOD)
        found = check_basecaller(
            bundle, "resources/models/rna004_130bps_sup@v5.3.0", "1.4.0"
        )
        assert found == []

    def test_wrong_model_is_an_error(self, tmp_path):
        bundle = write_bundle(tmp_path, GOOD)
        found = check_basecaller(bundle, "rna004_sup@v6.0.0", "1.4.0")
        assert [level for level, _ in found] == ["error"]
        assert "rna004_sup@v6.0.0" in found[0][1]
        assert "rna004_130bps_sup@v5.3.0" in found[0][1]

    def test_same_major_dorado_is_silent(self, tmp_path):
        """Same weights run by a later patch/minor are still the same weights."""
        bundle = write_bundle(tmp_path, GOOD)
        assert check_basecaller(bundle, "rna004_130bps_sup@v5.3.0", "1.9.3") == []

    def test_major_dorado_difference_only_warns(self, tmp_path):
        """The live case: model matches, dorado 1.x vs 2.x. Must not block."""
        bundle = write_bundle(tmp_path, GOOD)
        found = check_basecaller(bundle, "rna004_130bps_sup@v5.3.0", "2.1.1")
        assert [level for level, _ in found] == ["warn"]

    def test_unparseable_version_warns_rather_than_passing(self, tmp_path):
        bundle = write_bundle(tmp_path, GOOD)
        found = check_basecaller(bundle, "rna004_130bps_sup@v5.3.0", "")
        assert [level for level, _ in found] == ["warn"]

    def test_both_findings_report_error_first(self, tmp_path):
        bundle = write_bundle(tmp_path, GOOD)
        found = check_basecaller(bundle, "rna004_sup@v6.0.0", "2.1.1")
        assert [level for level, _ in found] == ["error", "warn"]

    def test_version_warning_does_not_claim_the_model_matches(self, tmp_path):
        """It fires alongside a model mismatch, so it must not assert the opposite."""
        bundle = write_bundle(tmp_path, GOOD)
        _, warning = check_basecaller(bundle, "rna004_sup@v6.0.0", "2.1.1")[1]
        assert "model matches" not in warning

    def test_undeclared_bundle_produces_nothing(self, tmp_path):
        """'Cannot tell' is not 'invalid' — the rule bundle_barcode_names follows."""
        bundle = write_bundle(tmp_path, None)
        assert check_basecaller(bundle, "anything@v9", "99.0.0") == []


class TestModelDigest:
    def test_path_and_contents_both_contribute(self, tmp_path):
        """Renaming a file must change the digest, or a swap would go unseen."""
        a = tmp_path / "a"
        a.mkdir()
        (a / "one.tensor").write_bytes(b"xyz")
        b = tmp_path / "b"
        b.mkdir()
        (b / "two.tensor").write_bytes(b"xyz")
        assert model_digest(a)[0] != model_digest(b)[0]

    def test_counts_files_recursively(self, tmp_path):
        model = tmp_path / "m"
        (model / "sub").mkdir(parents=True)
        (model / "config.toml").write_bytes(b"a")
        (model / "sub" / "w.tensor").write_bytes(b"b")
        assert model_digest(model)[1] == 2

    def test_stable_across_calls(self, tmp_path):
        model = tmp_path / "m"
        model.mkdir()
        (model / "config.toml").write_bytes(b"a")
        assert model_digest(model)[0] == model_digest(model)[0]

    def test_matches_upstreams_declared_hash(self):
        """The scheme is upstream's; if this drifts, nothing can ever verify.

        Runs only where the basecalling model is present — it is gitignored and
        downloaded by `pixi run setup`, so CI has neither it nor the 300 MB read.
        """
        block = bundle_basecaller(CHARGING_BUNDLE) if CHARGING_BUNDLE.is_dir() else None
        if block is None:
            pytest.skip("charging bundle is not vendored in this repo")
        model = REPO_ROOT / "resources" / "models" / block["model"]
        if not model.is_dir():
            pytest.skip("basecalling model is not downloaded (pixi run setup)")
        assert model_digest(model)[0] == block["model_sha256"]
