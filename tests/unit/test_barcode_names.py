"""
Unit tests for barcode_names.py

Tests translation between the barcode names a demux model emits and the names
this project uses, in both directions and for both panels.
"""

from pathlib import Path

import pytest

from barcode_names import bundle_barcode_names, emitted_to_label, label_to_emitted

REPO_ROOT = Path(__file__).resolve().parents[2]
DEMUX_MODELS = REPO_ROOT / "resources" / "models" / "demux"
WDX_BUNDLE = DEMUX_MODELS / "barcode_crf_wdx4_rna004@v0.2.0"
NBC_BUNDLE = DEMUX_MODELS / "barcode_crf_nbc16_rna004@v0.2.0"

# Guard, not an expectation of absence: the WDX4 bundle is vendored, so these
# run. It stays because `bundle_barcode_names` returns an empty set for a
# missing bundle, which would make every `not in` assertion below true for the
# wrong reason -- including the bc11 coverage check, whose whole job is to fail
# when a barcode is outside the panel.
needs_wdx_bundle = pytest.mark.skipif(
    not WDX_BUNDLE.is_dir(), reason="WDX4 CRF bundle is not vendored in this repo"
)


class TestEmittedToLabel:
    """Model output name -> project-facing name."""

    @pytest.mark.parametrize(
        "emitted,label",
        [("nbc01", "ldx01"), ("nbc08", "ldx08"), ("nbc16", "ldx16")],
    )
    def test_ldx_panel(self, emitted, label):
        """nbcNN is renamed to ldxNN."""
        assert emitted_to_label(emitted) == label

    @pytest.mark.parametrize(
        "emitted,label",
        [
            ("bc03", "barcode03"),
            ("bc04", "barcode04"),
            ("bc05", "barcode05"),
            ("bc07", "barcode07"),
        ],
    )
    def test_wdx_panel(self, emitted, label):
        """bcNN is renamed to barcodeNN."""
        assert emitted_to_label(emitted) == label

    def test_nbc_is_not_read_as_n_plus_bc(self):
        """The longer prefix wins, or nbc01 would become barcode01."""
        assert emitted_to_label("nbc01") == "ldx01"

    @pytest.mark.parametrize(
        "name", ["barcode03", "barcode07", "barcode11", "ldx01", "unclassified"]
    )
    def test_already_project_facing_is_unchanged(self, name):
        """Idempotent: a name already in our vocabulary passes through.

        This is also what keeps the legacy WarpDemuX path untouched -- every
        WarpDemuX barcode name is already project-facing.
        """
        assert emitted_to_label(name) == name

    @pytest.mark.parametrize("name", [None, "", "bcXX", "nbc", "bc"])
    def test_non_matching_input_passes_through(self, name):
        """Anything that is not <known prefix><digits> is returned unchanged."""
        assert emitted_to_label(name) == name


class TestLabelToEmitted:
    """Project-facing name -> the name to look for in the model's output."""

    @pytest.mark.parametrize(
        "label,emitted", [("ldx01", "nbc01"), ("ldx16", "nbc16")]
    )
    def test_ldx_panel(self, label, emitted):
        assert label_to_emitted(label) == emitted

    @pytest.mark.parametrize(
        "label,emitted",
        [
            ("barcode03", "bc03"),
            ("barcode04", "bc04"),
            ("barcode05", "bc05"),
            ("barcode07", "bc07"),
        ],
    )
    def test_wdx_panel(self, label, emitted):
        """This is the rename that makes escpod's barcode_bc03.pod5 findable."""
        assert label_to_emitted(label) == emitted

    @pytest.mark.parametrize("name", ["nbc01", "bc03", "bc07"])
    def test_already_emitted_is_unchanged(self, name):
        """A samples file written in the emitted vocabulary still works."""
        assert label_to_emitted(name) == name

    @pytest.mark.parametrize("name", [None, "", "unclassified", "bcXX"])
    def test_non_matching_input_passes_through(self, name):
        assert label_to_emitted(name) == name


class TestRoundTrip:
    """The two directions must compose back to identity on emitted names."""

    @pytest.mark.parametrize(
        "emitted", ["nbc01", "nbc16", "bc03", "bc04", "bc05", "bc07"]
    )
    def test_emitted_round_trip(self, emitted):
        assert label_to_emitted(emitted_to_label(emitted)) == emitted

    @pytest.mark.parametrize(
        "label", ["ldx01", "ldx16", "barcode03", "barcode07"]
    )
    def test_label_round_trip(self, label):
        assert emitted_to_label(label_to_emitted(label)) == label


class TestBundleBarcodeNames:
    """Reading a bundle's declared references out of metadata.json."""

    @needs_wdx_bundle
    def test_wdx_bundle_declares_four_codes(self):
        """The WDX CRF panel covers 4 of WarpDemuX's 12 codes."""
        assert bundle_barcode_names(WDX_BUNDLE) == {"bc03", "bc04", "bc05", "bc07"}

    def test_nbc_bundle_declares_sixteen_codes(self):
        names = bundle_barcode_names(NBC_BUNDLE)
        assert len(names) == 16
        assert {"nbc01", "nbc16"} <= names

    def test_non_bundle_path_returns_empty_set(self):
        """'Cannot tell' must read as 'do not validate', not 'nothing is valid'."""
        assert bundle_barcode_names(DEMUX_MODELS) == set()

    def test_missing_path_returns_empty_set(self):
        assert bundle_barcode_names(DEMUX_MODELS / "does-not-exist") == set()


class TestPanelCoverage:
    """The check that turns an uncovered barcode into a parse-time error."""

    @needs_wdx_bundle
    def test_covered_wdx_sample_validates(self):
        """A barcode03 sample resolves to a code the WDX bundle declares."""
        assert label_to_emitted("barcode03") in bundle_barcode_names(WDX_BUNDLE)

    @needs_wdx_bundle
    def test_uncovered_wdx_code_is_detectable(self):
        """barcode11 is a real WarpDemuX code the CRF panel does not cover.

        This is the case the pipeline must reject at DAG construction rather
        than discover after a multi-hour demux.
        """
        assert label_to_emitted("barcode11") not in bundle_barcode_names(WDX_BUNDLE)

    def test_wrong_bundle_is_detectable(self):
        """WDX samples against the nbc16 default bundle must not validate."""
        assert label_to_emitted("barcode03") not in bundle_barcode_names(NBC_BUNDLE)
