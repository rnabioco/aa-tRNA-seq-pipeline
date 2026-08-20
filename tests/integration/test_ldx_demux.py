"""
Integration tests for the LDX demultiplexing path.

Two tiers, deliberately:

  * Fixture and config tests run ANYWHERE, including CI. They assert that the
    committed fixture and the configs that reference it agree with each other.
    These are what stop issue #120 recurring — that bug was a config pointing at
    data which could not satisfy it, and it was invisible because `dry-run-demux`
    only builds the DAG.

  * Output tests need a completed run (`.tests/outputs-ldx`) and SKIP without
    one, since basecalling needs a GPU that CI does not have. Produce them with:

      pixi run snakemake --configfile=config/config-ldx-test.yml --cores 8

Note on determinism: barcode routing is decided from raw signal before
basecalling and reproduces exactly, so it is asserted exactly. 3' adapter
detection runs on basecalls, and GPU basecalling is not bit-reproducible, so
adapter counts are asserted with tolerance rather than pinned.
"""

import gzip
import subprocess
from collections import Counter
from pathlib import Path

import pysam
import pytest
import yaml

REPO_ROOT = Path(__file__).parent.parent.parent
FIXTURE = REPO_ROOT / ".tests" / "fixtures" / "ldx-demux"
MANIFEST = FIXTURE / "fixture_manifest.tsv"
OUTPUTS = REPO_ROOT / ".tests" / "outputs-ldx"
CONFIG = REPO_ROOT / "config" / "config-ldx-test.yml"
SAMPLES = REPO_ROOT / "config" / "samples-ldx-test.yml"

# sample -> (barcode, edx assignment or None)
EXPECTED_SAMPLES = {
    "ldx01_edx01": ("nbc01", "edx01"),
    "ldx02_edx02": ("nbc02", "edx02"),
    "ldx08_pool": ("nbc08", None),
}

needs_run = pytest.mark.skipif(
    not (OUTPUTS / "bam" / "final").exists(),
    reason="no completed LDX run; see module docstring",
)


def manifest_populations():
    """selected_as label -> read count, from the committed fixture manifest."""
    counts = Counter()
    with open(MANIFEST) as fh:
        next(fh)
        for line in fh:
            _read_id, selected_as = line.rstrip("\n").split("\t")
            counts[selected_as] += 1
    return counts


def barcode_of(selected_as):
    """`nbc01:edx01:on_target` -> `nbc01`; `unclassified` -> `unclassified`."""
    return selected_as.split(":")[0]


# --------------------------------------------------------------------------
# Fixture integrity — runs in CI
# --------------------------------------------------------------------------


class TestFixture:
    def test_fixture_present(self):
        pod5s = list((FIXTURE / "run" / "pod5").glob("*.pod5"))
        assert pod5s, "no POD5 in the fixture run directory"
        # find_raw_inputs globs pod5/, so an empty or missing dir would surface
        # as an empty DAG rather than an error
        assert sum(p.stat().st_size for p in pod5s) > 1_000_000

    def test_manifest_matches_documented_composition(self):
        counts = manifest_populations()
        assert counts == {
            "nbc01:edx01:on_target": 100,
            "nbc01:other_adapter": 25,
            "nbc02:edx02:on_target": 100,
            "nbc02:other_adapter": 25,
            "nbc08:pool_aligned": 100,
            "nbc05:unclaimed": 40,
            "unclassified": 25,
        }

    def test_read_ids_unique(self):
        ids = [ln.split("\t")[0] for ln in MANIFEST.read_text().splitlines()[1:]]
        assert len(ids) == len(set(ids)) == 415

    def test_reference_is_raw_trna(self):
        """reference.mode is `build`, so collapsed.fa must be UNadapted."""
        fa = (FIXTURE / "collapsed.fa").read_text()
        assert fa.startswith(">")
        seqs = [s for s in fa.split(">")[1:]]
        assert len(seqs) == 47
        # the 5' adapter must not already be present, or build would double it
        assert "CCTAAGAGCAAGAAGAAGCCTGG" not in fa

    def test_off_target_reads_exist(self):
        """
        The filtered samples must carry wrong-adapter reads, or
        filter_{fastq,pod5}_by_edx is a no-op that cannot fail.
        """
        counts = manifest_populations()
        assert counts["nbc01:other_adapter"] > 0
        assert counts["nbc02:other_adapter"] > 0


# --------------------------------------------------------------------------
# Config agreement — runs in CI, and is the direct regression test for #120
# --------------------------------------------------------------------------


@pytest.fixture(scope="module")
def config():
    return yaml.safe_load(CONFIG.read_text())


@pytest.fixture(scope="module")
def samples():
    return yaml.safe_load(SAMPLES.read_text())


class TestConfig:
    def test_samples_point_at_the_fixture(self, samples):
        (run,) = samples["runs"]
        assert (REPO_ROOT / run["path"]).is_dir()
        assert set(run["samples"]) == set(EXPECTED_SAMPLES)

    def test_barcode_assignments_match_fixture(self, samples):
        """Every `ldx:` assignment must name a barcode the fixture contains."""
        present = {barcode_of(s) for s in manifest_populations()}
        (run,) = samples["runs"]
        for name, val in run["samples"].items():
            barcode, edx = EXPECTED_SAMPLES[name]
            assert val["ldx"] == barcode
            assert val.get("edx") == edx
            assert barcode in present, f"{name} claims {barcode}, absent from fixture"

    def test_ldx_enabled_and_wdx_not(self, config):
        assert config["ldx"]["enabled"] is True
        assert "warpdemux" not in config, "backends are mutually exclusive"

    def test_edx_adapters_are_overridden(self, config):
        """
        The exact failure in #120: config-base.yml defines edx01/edx02 with the
        sacCer3 dual-adapter sequences. Same names, different molecules. Without
        this override every fixture read detects as `none`.
        """
        base = yaml.safe_load((REPO_ROOT / "config" / "config-base.yml").read_text())
        base_seqs = {a["name"]: a["seq"] for a in base["adapters"]["three_prime"]}
        test_seqs = {a["name"]: a["seq"] for a in config["adapters"]["three_prime"]}
        assert test_seqs["edx01"] != base_seqs["edx01"]
        # all seven of the donor run's adapters must be listed, so that
        # off-target reads resolve to an adapter instead of falling to `none`
        assert set(test_seqs) == {f"edx0{i}" for i in range(1, 8)}

    def test_every_assigned_edx_is_a_declared_adapter(self, config, samples):
        declared = {a["name"] for a in config["adapters"]["three_prime"]}
        (run,) = samples["runs"]
        for name, val in run["samples"].items():
            if val.get("edx"):
                assert val["edx"] in declared

    def test_reference_build_path_resolves(self, config):
        assert config["reference"]["mode"] == "build"
        assert (REPO_ROOT / config["reference"]["raw_fasta"]).is_file()


# --------------------------------------------------------------------------
# Pipeline outputs — needs a completed run
# --------------------------------------------------------------------------


@needs_run
class TestDemuxRouting:
    def test_routing_reproduces_fixture_exactly(self):
        """
        Barcode assignment is made from raw signal, so a rerun must place every
        read in the barcode it was selected for.
        """
        (summary,) = (OUTPUTS / "demux" / "read_ids").glob("*/demux_summary.tsv.gz")
        with gzip.open(summary, "rt") as fh:
            next(fh)
            observed = {
                ln.split("\t")[0]: int(ln.split("\t")[1]) for ln in fh if ln.strip()
            }
        expected = Counter()
        for label, n in manifest_populations().items():
            expected[barcode_of(label)] += n
        assert observed == dict(expected)

    def test_unclaimed_barcode_produces_no_sample(self):
        """nbc05 is in the fixture but claimed by no sample."""
        assert not (OUTPUTS / "bam" / "final" / "ldx05").exists()
        assert not list((OUTPUTS / "bam" / "final").glob("*nbc05*"))


@needs_run
class TestEdxFiltering:
    @pytest.mark.parametrize("sample,target", [("ldx01_edx01", "edx01"), ("ldx02_edx02", "edx02")])
    def test_filter_keeps_only_the_target_adapter(self, sample, target):
        kept = set(
            (OUTPUTS / "demux" / "edx" / sample / f"{sample}.edx_read_ids.txt")
            .read_text()
            .split()
        )
        detected = {}
        with gzip.open(
            OUTPUTS / "demux" / "edx" / sample / f"{sample}.edx_adapters.tsv.gz", "rt"
        ) as fh:
            next(fh)
            for ln in fh:
                f = ln.split("\t")
                detected[f[0]] = f[1]
        assert kept, "filter kept nothing"
        assert {detected[r] for r in kept} == {target}

    @pytest.mark.parametrize("sample", ["ldx01_edx01", "ldx02_edx02"])
    def test_filter_actually_removed_reads(self, sample):
        """A filter that drops nothing would pass every other assertion here."""
        detected_n = 0
        with gzip.open(
            OUTPUTS / "demux" / "edx" / sample / f"{sample}.edx_adapters.tsv.gz", "rt"
        ) as fh:
            next(fh)
            detected_n = sum(1 for ln in fh if ln.strip())
        kept_n = len(
            (OUTPUTS / "demux" / "edx" / sample / f"{sample}.edx_read_ids.txt")
            .read_text()
            .split()
        )
        assert kept_n < detected_n, "EDX filter was a no-op"
        # ~25 off-target reads were planted; allow drift from basecall variation
        assert 10 <= detected_n - kept_n <= 45

    def test_unfiltered_sample_is_not_filtered(self):
        """ldx08_pool has no `edx:` key, so no filtering artefacts may exist."""
        assert not (OUTPUTS / "demux" / "edx" / "ldx08_pool").exists()

    def test_concordance_covers_filtered_samples(self):
        with gzip.open(OUTPUTS / "summary" / "edx" / "edx_concordance.tsv.gz", "rt") as fh:
            header = next(fh).rstrip("\n").split("\t")
            rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
        assert header[:3] == ["sample", "edx_adapter", "n_reads"]
        assert {r[0] for r in rows} == {"ldx01_edx01", "ldx02_edx02"}


@needs_run
class TestFinalBamTags:
    @pytest.mark.parametrize("sample,barcode", [(s, b) for s, (b, _) in EXPECTED_SAMPLES.items()])
    def test_barcode_tag_on_every_read(self, sample, barcode):
        """
        BC carries OUR name for the barcode (ldx01), not upstream's (nbc01),
        and it must be on every read — a partially tagged BAM is worse than an
        untagged one.
        """
        expected = barcode.replace("nbc", "ldx")
        bam = OUTPUTS / "bam" / "final" / sample / f"{sample}.bam"
        with pysam.AlignmentFile(bam, "rb") as fh:
            tags = [r.get_tag("BC") if r.has_tag("BC") else None for r in fh]
        assert tags, "empty BAM"
        assert set(tags) == {expected}

    @pytest.mark.parametrize("sample", list(EXPECTED_SAMPLES))
    def test_read_groups_resolve(self, sample):
        """Regression test for the dangling @RG that PR #121 fixed."""
        bam = OUTPUTS / "bam" / "final" / sample / f"{sample}.bam"
        with pysam.AlignmentFile(bam, "rb") as fh:
            declared = {rg["ID"] for rg in fh.header.to_dict().get("RG", [])}
            used = {r.get_tag("RG") for r in fh if r.has_tag("RG")}
        assert declared, "no @RG in header"
        assert used, "no read carries RG"
        assert used <= declared, f"dangling RG: {used - declared}"

    @pytest.mark.parametrize("sample,barcode", [(s, b) for s, (b, _) in EXPECTED_SAMPLES.items()])
    def test_upstream_barcode_recorded_in_comment(self, sample, barcode):
        bam = OUTPUTS / "bam" / "final" / sample / f"{sample}.bam"
        with pysam.AlignmentFile(bam, "rb") as fh:
            comments = fh.header.to_dict().get("CO", [])
        assert f"aa-tRNA-seq:upstream_barcode={barcode}" in comments

    @pytest.mark.parametrize("sample", list(EXPECTED_SAMPLES))
    def test_sample_name_in_read_group(self, sample):
        bam = OUTPUTS / "bam" / "final" / sample / f"{sample}.bam"
        with pysam.AlignmentFile(bam, "rb") as fh:
            rgs = fh.header.to_dict().get("RG", [])
        assert [rg.get("SM") for rg in rgs] == [sample]

    def test_bam_is_valid_sam(self):
        """quickcheck passed even with the dangling @RG; this is stricter."""
        bam = OUTPUTS / "bam" / "final" / "ldx01_edx01" / "ldx01_edx01.bam"
        p = subprocess.run(
            ["samtools", "view", "-h", "--no-PG", str(bam)],
            capture_output=True, text=True,
        )
        assert p.returncode == 0, p.stderr


@needs_run
class TestAttrition:
    def test_attrition_starts_from_the_whole_fixture(self):
        with gzip.open(OUTPUTS / "summary" / "read_attrition.tsv.gz", "rt") as fh:
            header = next(fh).rstrip("\n").split("\t")
            rows = [dict(zip(header, ln.rstrip("\n").split("\t"))) for ln in fh if ln.strip()]
        assert rows
        assert int(rows[0]["entered"]) == 415
        # each stage must hand its retained count to the next
        for prev, nxt in zip(rows, rows[1:]):
            assert int(prev["retained"]) == int(nxt["entered"])
