"""
Integration tests that validate pre-computed pipeline outputs.

These tests verify that expected output files exist and have the correct structure.
They use pre-computed GPU outputs as fixtures.
"""

import gzip
import pandas as pd
import pysam
import pytest
from pathlib import Path

# Test data paths
REPO_ROOT = Path(__file__).parent.parent.parent
TEST_OUTPUTS_DIR = REPO_ROOT / ".tests" / "outputs"


def test_outputs_exist():
    """Verify pre-computed test outputs directory exists."""
    if not TEST_OUTPUTS_DIR.exists():
        pytest.skip("Test outputs not available. Run 'pixi run dl-test-data'")
    assert TEST_OUTPUTS_DIR.is_dir()


@pytest.mark.skipif(
    not (TEST_OUTPUTS_DIR / "bam" / "final").exists(),
    reason="Pre-computed outputs not available",
)
class TestFinalBamOutput:
    """Tests for final BAM output files."""

    @pytest.fixture
    def final_bam_path(self):
        """Get path to sample1 final BAM."""
        path = TEST_OUTPUTS_DIR / "bam" / "final" / "sample1.bam"
        if not path.exists():
            pytest.skip("Final BAM not found")
        return path

    def test_final_bam_is_valid(self, final_bam_path):
        """Final BAM should be a valid BAM file."""
        with pysam.AlignmentFile(str(final_bam_path), "rb") as bam:
            # Should be able to iterate without errors
            reads = list(bam.fetch())
            assert len(reads) > 0

    def test_final_bam_has_cl_tag(self, final_bam_path):
        """Final BAM reads should have CL (charging likelihood) tag."""
        with pysam.AlignmentFile(str(final_bam_path), "rb") as bam:
            reads_with_cl = 0
            total_reads = 0
            for read in bam.fetch():
                total_reads += 1
                if read.has_tag("CL"):
                    reads_with_cl += 1
                    # CL should be in valid range
                    cl_val = read.get_tag("CL")
                    if isinstance(cl_val, (list, tuple)):
                        cl_val = cl_val[0]
                    assert 0 <= cl_val <= 255

            # Most reads should have CL tag
            assert reads_with_cl > 0
            assert reads_with_cl / total_reads > 0.5

    def test_final_bam_has_pt_tag(self, final_bam_path):
        """Final BAM reads should have PT (adapter positions) tag."""
        with pysam.AlignmentFile(str(final_bam_path), "rb") as bam:
            reads_with_pt = 0
            for read in bam.fetch():
                if read.has_tag("PT"):
                    reads_with_pt += 1
                    pt_val = read.get_tag("PT")
                    # PT should have expected format
                    assert isinstance(pt_val, str)
                    # Should contain adapter annotations
                    if "|" in pt_val:
                        # Both adapters present
                        assert "5p_adapter" in pt_val
                        assert "3p_adapter" in pt_val
                    else:
                        # Single adapter
                        assert "adapter" in pt_val

            assert reads_with_pt > 0

    def test_final_bam_index_exists(self, final_bam_path):
        """Final BAM should have accompanying index."""
        index_path = Path(str(final_bam_path) + ".bai")
        assert index_path.exists()


@pytest.mark.skipif(
    not (TEST_OUTPUTS_DIR / "summary" / "tables").exists(),
    reason="Pre-computed outputs not available",
)
class TestChargingTableOutput:
    """Tests for charging probability and CPM tables."""

    @pytest.fixture
    def charging_prob_path(self):
        """Get path to charging probability table."""
        path = (
            TEST_OUTPUTS_DIR
            / "summary"
            / "tables"
            / "sample1"
            / "sample1.charging_prob.tsv.gz"
        )
        if not path.exists():
            pytest.skip("Charging probability table not found")
        return path

    @pytest.fixture
    def charging_cpm_path(self):
        """Get path to charging CPM table."""
        path = (
            TEST_OUTPUTS_DIR
            / "summary"
            / "tables"
            / "sample1"
            / "sample1.charging.cpm.tsv.gz"
        )
        if not path.exists():
            pytest.skip("Charging CPM table not found")
        return path

    def test_charging_prob_format(self, charging_prob_path):
        """Charging probability table should have correct format."""
        df = pd.read_csv(charging_prob_path, sep="\t")

        # Required columns
        assert "read_id" in df.columns
        assert "tRNA" in df.columns
        assert "charging_likelihood" in df.columns

        # Should have data
        assert len(df) > 0

        # Likelihood should be in valid range
        assert df["charging_likelihood"].min() >= 0
        assert df["charging_likelihood"].max() <= 255

    def test_charging_cpm_format(self, charging_cpm_path):
        """Charging CPM table should have correct format."""
        df = pd.read_csv(charging_cpm_path, sep="\t", index_col=0)

        # Should have count columns
        has_charged = "counts_charged" in df.columns
        has_uncharged = "counts_uncharged" in df.columns
        assert has_charged or has_uncharged

        # Should have CPM columns
        has_cpm_charged = "cpm_charged" in df.columns
        has_cpm_uncharged = "cpm_uncharged" in df.columns
        assert has_cpm_charged or has_cpm_uncharged

        # Should have tRNA entries
        assert len(df) > 0

    def test_charging_cpm_values(self, charging_cpm_path):
        """CPM values should be non-negative and reasonable."""
        df = pd.read_csv(charging_cpm_path, sep="\t", index_col=0)

        if "cpm_charged" in df.columns:
            assert df["cpm_charged"].min() >= 0

        if "cpm_uncharged" in df.columns:
            assert df["cpm_uncharged"].min() >= 0

        if "counts_charged" in df.columns:
            assert df["counts_charged"].min() >= 0

        if "counts_uncharged" in df.columns:
            assert df["counts_uncharged"].min() >= 0


@pytest.mark.skipif(
    not (TEST_OUTPUTS_DIR / "summary" / "tables").exists(),
    reason="Pre-computed outputs not available",
)
class TestAlignmentStatsOutput:
    """Tests for alignment statistics output."""

    @pytest.fixture
    def align_stats_path(self):
        """Get path to alignment stats file."""
        path = (
            TEST_OUTPUTS_DIR
            / "summary"
            / "tables"
            / "sample1"
            / "sample1.align_stats.txt"
        )
        if not path.exists():
            pytest.skip("Alignment stats not found")
        return path

    def test_align_stats_exists(self, align_stats_path):
        """Alignment stats file should exist and have content."""
        assert align_stats_path.exists()
        content = align_stats_path.read_text()
        assert len(content) > 0

    def test_align_stats_has_key_metrics(self, align_stats_path):
        """Alignment stats should contain expected metrics."""
        content = align_stats_path.read_text()

        # Check for expected metric names
        expected_metrics = [
            "total_alignments",
            "passed_alignments",
        ]

        for metric in expected_metrics:
            assert metric in content, f"Missing metric: {metric}"


@pytest.mark.skipif(
    not (TEST_OUTPUTS_DIR / "reference").exists(),
    reason="Pre-computed outputs not available",
)
class TestReferenceOutput:
    """Tests for reference validation/build output."""

    @pytest.fixture
    def validated_ref_path(self):
        """Get path to validated reference."""
        # Look for any .fa file in reference output
        ref_dir = TEST_OUTPUTS_DIR / "reference"
        if not ref_dir.exists():
            pytest.skip("Reference output not found")

        fa_files = list(ref_dir.glob("**/*.fa"))
        if not fa_files:
            pytest.skip("No FASTA files in reference output")
        return fa_files[0]

    def test_validated_reference_is_fasta(self, validated_ref_path):
        """Validated reference should be valid FASTA."""
        with open(validated_ref_path) as f:
            content = f.read()

        # Should have FASTA header
        assert content.startswith(">")

        # Should have sequences
        lines = content.strip().split("\n")
        seq_lines = [l for l in lines if not l.startswith(">")]
        assert len(seq_lines) > 0
