"""
Integration tests for the dual-index (LDX + FDX) demultiplexing path.

Same two tiers as test_ldx_demux.py:

  * Fixture and config tests run ANYWHERE, including CI, and hold the committed
    fixture and the configs that reference it to each other.

  * Output tests need a completed run (`.tests/outputs-fdx`) and SKIP without
    one, since basecalling needs a GPU. Produce them with:

      pixi run test-fdx

What is special about this fixture: the fdx label of every read comes from the
LDX call through the donor pool's design (each library = one 5' code x three 3'
codes), an independent channel from the 5' signal the fdx model reads. And the
donor flowcell (Run2) is not in the shipped fdx model's training corpus, so the
per-sample recovery measured here is a held-out number.
"""

import gzip
from collections import Counter
from pathlib import Path

import pysam
import pytest
import yaml

REPO_ROOT = Path(__file__).parent.parent.parent
FIXTURE = REPO_ROOT / ".tests" / "fixtures" / "fdx-demux"
MANIFEST = FIXTURE / "fixture_manifest.tsv"
OUTPUTS = REPO_ROOT / ".tests" / "outputs-fdx"
CONFIG = REPO_ROOT / "config" / "config-fdx-test.yml"
SAMPLES = REPO_ROOT / "config" / "samples-fdx-test.yml"
RUN_ID = "run"

# sample -> (ldx, fdx)
EXPECTED_SAMPLES = {
    "fdx01_ldx01": ("ldx01", "fdx01"),
    "fdx02_ldx04": ("ldx04", "fdx02"),
    "fdx03_ldx07": ("ldx07", "fdx03"),
}
# The plain 3' adapter, recovered from the donor run's soft-clipped tails
# (escapepod-models dev-notes/fdx-three-prime-chemistry.md).
PLAIN_3P = "GGCTTCTTCTTGCTCTTAGGAAAAAAAAAA"

needs_run = pytest.mark.skipif(
    not (OUTPUTS / "bam" / "final").exists(),
    reason="no completed FDX run; see module docstring",
)


def manifest_populations():
    counts = Counter()
    with open(MANIFEST) as fh:
        next(fh)
        for line in fh:
            _read_id, selected_as = line.rstrip("\n").split("\t")
            counts[selected_as] += 1
    return counts


def codes_of(selected_as):
    """`ldx01:fdx01:aligned` -> (`ldx01`, `fdx01`); `unclassified` -> None."""
    parts = selected_as.split(":")
    return (parts[0], parts[1]) if len(parts) >= 2 else None


# --------------------------------------------------------------------------
# Fixture integrity — runs in CI
# --------------------------------------------------------------------------


class TestFixture:
    def test_fixture_present(self):
        pod5s = list((FIXTURE / "run" / "pod5").glob("*.pod5"))
        assert pod5s, "no POD5 in the fixture run directory"
        assert sum(p.stat().st_size for p in pod5s) > 1_000_000

    def test_manifest_matches_documented_composition(self):
        assert manifest_populations() == {
            "ldx01:fdx01:aligned": 100,
            "ldx04:fdx02:aligned": 100,
            "ldx07:fdx03:aligned": 100,
            "ldx10:fdx04:unclaimed": 40,
            "unclassified": 25,
        }

    def test_read_ids_unique(self):
        ids = [ln.split("\t")[0] for ln in MANIFEST.read_text().splitlines()[1:]]
        assert len(ids) == len(set(ids)) == 365

    def test_reference_is_raw_trna(self):
        """reference.mode is `build`, so collapsed.fa must be UNadapted."""
        fa = (FIXTURE / "collapsed.fa").read_text()
        assert fa.startswith(">")
        seqs = fa.split(">")[1:]
        assert len(seqs) == 47
        assert "CCTAAGAGCAAGAAGAAGCCTGG" not in fa
        bodies = ["".join(s.split("\n")[1:]) for s in seqs]
        assert all(b.endswith("CCA") for b in bodies)


# --------------------------------------------------------------------------
# Config agreement — runs in CI
# --------------------------------------------------------------------------


@pytest.fixture(scope="module")
def config():
    return yaml.safe_load(CONFIG.read_text())


@pytest.fixture(scope="module")
def samples():
    return yaml.safe_load(SAMPLES.read_text())


@pytest.fixture(scope="module")
def base_config():
    return yaml.safe_load((REPO_ROOT / "config" / "config-base.yml").read_text())


class TestConfig:
    def test_samples_point_at_the_fixture(self, samples):
        (run,) = samples["runs"]
        assert (REPO_ROOT / run["path"]).is_dir()
        assert set(run["samples"]) == set(EXPECTED_SAMPLES)

    def test_assignments_match_fixture(self, samples):
        present = {codes_of(s) for s in manifest_populations()} - {None}
        (run,) = samples["runs"]
        for name, val in run["samples"].items():
            ldx, fdx = EXPECTED_SAMPLES[name]
            assert (val["ldx"], val["fdx"]) == (ldx, fdx)
            assert (ldx, fdx) in present, f"{name} claims {ldx}+{fdx}, absent from fixture"

    def test_both_axes_enabled_and_wdx_not(self, config):
        assert config["ldx"]["enabled"] is True
        assert config["fdx"]["enabled"] is True
        assert "warpdemux" not in config, "backends are mutually exclusive"

    def test_fdx_model_and_gate_come_from_base(self, config, base_config):
        """The test must exercise what ships, so a bump is covered not bypassed."""
        assert "model" not in config["fdx"]
        assert "min_crf_margin" not in config["fdx"]
        bundle = REPO_ROOT / base_config["fdx"]["model"]
        assert (bundle / "metadata.json").is_file()
        assert base_config["fdx"]["min_crf_margin"] == 3.5

    def test_fdx_codes_are_declared_by_the_bundle(self, base_config):
        import json

        meta = json.loads((REPO_ROOT / base_config["fdx"]["model"] / "metadata.json").read_text())
        declared = {b["name"] for b in meta["barcodes"]}
        for _, fdx in EXPECTED_SAMPLES.values():
            assert fdx in declared
        # and the bundle is the read-end-anchored kind the two-pass rule assumes
        assert meta["signal"]["anchor"] == "read_end"

    def test_plain_three_prime_adapter(self, config):
        """
        Every library in the fixture used the plain 3' adapter, which is the
        EDX scaffold truncated before the barcode. A single string, since the
        reference build requires one 3' adapter length.
        """
        assert config["adapters"]["three_prime"] == PLAIN_3P
        assert PLAIN_3P.startswith("GGC"), "charging anchors on the CCA|GGC junction"

    def test_reference_build_path_resolves(self, config):
        assert config["reference"]["mode"] == "build"
        assert (REPO_ROOT / config["reference"]["raw_fasta"]).is_file()


# --------------------------------------------------------------------------
# Pipeline outputs — needs a completed run
# --------------------------------------------------------------------------


def read_summary(path):
    with gzip.open(path, "rt") as fh:
        next(fh)
        return {ln.split("\t")[0]: int(ln.split("\t")[1]) for ln in fh if ln.strip()}


@needs_run
class TestDemuxRouting:
    def test_ldx_axis_broadly_reproduces_the_fixture(self):
        """The 3' axis, same bounded-concordance shape as the LDX fixture test."""
        observed = read_summary(OUTPUTS / "demux" / "read_ids" / RUN_ID / "demux_summary.tsv.gz")
        expected = Counter()
        for label, n in manifest_populations().items():
            codes = codes_of(label)
            expected[codes[0] if codes else "unclassified"] += n
        assert sum(observed.values()) == sum(expected.values()), "reads went missing"
        for code, n in expected.items():
            if code == "unclassified":
                continue
            assert observed.get(code, 0) >= 0.85 * n, f"{code}: {observed.get(code, 0)} of {n}"

    def test_fdx_axis_agrees_with_the_library_design(self):
        """
        The 5' axis against the LDX-derived truth. Held out: the donor flowcell
        is not in this model's corpus. The bundle publishes balanced recall
        0.938 at its gate on its own held-out split; assert the shape, not the
        decimal.
        """
        observed = read_summary(
            OUTPUTS / "demux" / "read_ids" / RUN_ID / "fdx" / "demux_summary.tsv.gz"
        )
        expected = Counter()
        for label, n in manifest_populations().items():
            codes = codes_of(label)
            if codes:
                expected[codes[1]] += n
        assert sum(observed.values()) == 365, "reads went missing"
        for code, n in expected.items():
            assert observed.get(code, 0) >= 0.80 * n, f"{code}: {observed.get(code, 0)} of {n}"

    def test_join_keeps_most_of_each_claimed_pair(self):
        summary = OUTPUTS / "demux" / "read_ids" / RUN_ID / "assigned_summary.tsv"
        rows = [ln.split("\t") for ln in summary.read_text().splitlines()[1:]]
        assigned = {r[0]: int(r[1]) for r in rows}
        assert set(assigned) == set(EXPECTED_SAMPLES)
        for name, n in assigned.items():
            # both axes must agree, so this is the product of two recoveries
            assert n >= 75, f"{name}: only {n} of 100 selected reads survived the join"

    def test_unclaimed_pair_produces_no_sample(self):
        assert not list((OUTPUTS / "bam" / "final").glob("*ldx10*"))
        assert not list((OUTPUTS / "bam" / "final").glob("*fdx04*"))


@needs_run
class TestFinalBamTags:
    @pytest.mark.parametrize("sample,codes", list(EXPECTED_SAMPLES.items()))
    def test_dual_index_barcode_tag_on_every_read(self, sample, codes):
        """BC carries both codes, 3' first, joined with `-` (the SAM dual-index form)."""
        ldx, fdx = codes
        bam = OUTPUTS / "bam" / "final" / sample / f"{sample}.bam"
        with pysam.AlignmentFile(bam, "rb") as fh:
            tags = {r.get_tag("BC") if r.has_tag("BC") else None for r in fh}
            rgs = fh.header.to_dict().get("RG", [])
        assert tags == {f"{ldx}-{fdx}"}
        assert [rg.get("SM") for rg in rgs] == [sample]
        assert {rg.get("BC") for rg in rgs} == {f"{ldx}-{fdx}"}

    @pytest.mark.parametrize("sample", list(EXPECTED_SAMPLES))
    def test_reads_are_charge_called(self, sample):
        bam = OUTPUTS / "bam" / "final" / sample / f"{sample}.bam"
        with pysam.AlignmentFile(bam, "rb") as fh:
            n = sum(1 for _ in fh)
            fh.reset()
            called = sum(1 for r in fh if r.has_tag("cl"))
        assert n > 0
        assert called > 0.5 * n, f"{sample}: {called} of {n} reads carry cl"
