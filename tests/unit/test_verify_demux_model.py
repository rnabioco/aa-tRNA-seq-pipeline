"""Unit tests for scripts/verify-demux-model.py.

The script is hyphenated and lives outside workflow/scripts, so it is loaded by
path rather than imported by name.
"""

import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).parent.parent.parent
SCRIPT = REPO_ROOT / "scripts" / "verify-demux-model.py"


def _load():
    spec = importlib.util.spec_from_file_location("verify_demux_model", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


vdm = _load()


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


@pytest.fixture
def bundle(temp_dir):
    """A minimal bundle whose ONNX hashes agree with its own metadata."""
    demux = temp_dir / "demux"
    b = demux / "barcode_crf_test_rna004@v1.0.0"
    b.mkdir(parents=True)

    crf = b"crf-weights"
    boundary = b"boundary-weights"
    (b / "barcode_crf_test_rna004.onnx").write_bytes(crf)
    (b / "adapter_rna004.onnx").write_bytes(boundary)

    (b / "metadata.json").write_text(json.dumps({
        "onnx": "barcode_crf_test_rna004.onnx",
        "boundary": {
            "onnx": "adapter_rna004.onnx",
            "sha256": _sha256(boundary),
        },
    }))
    (b / "provenance.json").write_text(json.dumps({"sha256": _sha256(crf)}))
    return b


def _set_metadata(bundle, **keys):
    meta = json.loads((bundle / "metadata.json").read_text())
    meta["boundary"].update(keys)
    (bundle / "metadata.json").write_text(json.dumps(meta))


class TestSidecarDigests:
    def test_unrecorded_bundle_fails_closed(self, bundle):
        """No recorded digest must FAIL, never silently skip.

        This is the whole value of the check: it replaced a blanket refusal of
        any bundle declaring window rules, so a bundle nobody recorded has to
        be refused rather than waved through.
        """
        assert vdm.verify(bundle) is False

    def test_recorded_bundle_verifies(self, bundle):
        vdm.record([str(bundle)])
        assert vdm.verify(bundle) is True

    def test_edited_metadata_fails(self, bundle):
        """The nbc16 failure mode: window rules patched in locally.

        Every ONNX hash still agrees — metadata.json is where those hashes
        live, so nothing upstream covers it. Only the recorded digest catches
        this.
        """
        vdm.record([str(bundle)])
        _set_metadata(bundle, margin=0, clamp_max_shift=300)
        assert vdm.verify(bundle) is False

    def test_edited_provenance_fails(self, bundle):
        vdm.record([str(bundle)])
        prov = bundle / "provenance.json"
        prov.write_text(prov.read_text().replace("sha256", "sha256 "))
        assert vdm.verify(bundle) is False


class TestUpstreamWindowRules:
    def test_declared_window_rules_are_accepted(self, bundle):
        """escapepod-models#127/#128 ship these legitimately as of 2026-09-06.

        The old guard refused any bundle declaring them, which made the current
        `latest` of two families unverifiable.
        """
        _set_metadata(bundle, margin=200, clamp_max_shift=0)
        vdm.record([str(bundle)])
        assert vdm.verify(bundle) is True


class TestLocalOverlay:
    def test_local_record_does_not_touch_tracked_manifest(self, bundle):
        demux = bundle.parent
        tracked = demux / vdm.RECORD_NAME

        vdm.record([str(bundle)])
        before = tracked.read_text()

        extra = demux / "barcode_crf_other_rna004@v1.0.0"
        extra.mkdir()
        for name in ("metadata.json", "provenance.json"):
            (extra / name).write_text((bundle / name).read_text())
        for name in ("barcode_crf_test_rna004.onnx", "adapter_rna004.onnx"):
            (extra / name).write_bytes((bundle / name).read_bytes())

        vdm.record([str(extra)], local=True)

        assert tracked.read_text() == before
        assert (demux / vdm.LOCAL_RECORD_NAME).is_file()
        assert vdm.verify(extra) is True

    def test_recording_one_bundle_keeps_the_others(self, bundle):
        """Recording by path MERGES; it must not drop the other bundles.

        Replacing would fail every other local bundle closed the moment a
        deployment recorded a new one by path — the same trap the
        tracked/local split exists to prevent, through a different door.
        """
        demux = bundle.parent
        other = demux / "barcode_crf_other_rna004@v1.0.0"
        other.mkdir()
        for name in ("metadata.json", "provenance.json"):
            (other / name).write_text((bundle / name).read_text())
        for name in ("barcode_crf_test_rna004.onnx", "adapter_rna004.onnx"):
            (other / name).write_bytes((bundle / name).read_bytes())

        vdm.record([str(bundle)], local=True)
        vdm.record([str(other)], local=True)

        assert vdm.verify(bundle) is True
        assert vdm.verify(other) is True

    def test_record_drops_entries_whose_bundle_is_gone(self, bundle, temp_dir):
        """A manifest should describe what is on disk, not accumulate."""
        import shutil

        demux = bundle.parent
        stale = demux / "barcode_crf_stale_rna004@v1.0.0"
        stale.mkdir()
        for name in ("metadata.json", "provenance.json"):
            (stale / name).write_text((bundle / name).read_text())

        vdm.record([str(bundle), str(stale)])
        assert stale.name in (demux / vdm.RECORD_NAME).read_text()

        shutil.rmtree(stale)
        vdm.record([str(bundle)])
        assert stale.name not in (demux / vdm.RECORD_NAME).read_text()
        assert bundle.name in (demux / vdm.RECORD_NAME).read_text()

    def test_local_record_skips_bundles_already_tracked(self, bundle):
        demux = bundle.parent
        vdm.record([str(bundle)])
        vdm.record([str(bundle)], local=True)
        overlay = (demux / vdm.LOCAL_RECORD_NAME).read_text()
        assert bundle.name not in overlay
