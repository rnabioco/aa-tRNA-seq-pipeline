"""
Integration tests that validate pre-computed pipeline outputs.

These tests verify that expected output files exist and have the correct structure.
They use pre-computed GPU outputs as fixtures.
"""

import gzip
import json
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
        path = TEST_OUTPUTS_DIR / "bam" / "final" / "sample1" / "sample1.bam"
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
        """Final BAM reads should have the cl (charging likelihood) tag.

        Lowercase: two-letter lowercase tags are the SAM spec's local-use
        space, and tags are case-sensitive. A read with NO cl tag is a
        deliberate no-call from the charging model, not a failure, so this
        asserts a majority rather than all.
        """
        with pysam.AlignmentFile(str(final_bam_path), "rb") as bam:
            reads_with_cl = 0
            total_reads = 0
            for read in bam.fetch():
                total_reads += 1
                if read.has_tag("cl"):
                    reads_with_cl += 1
                    # cl should be in valid range
                    cl_val = read.get_tag("cl")
                    if not isinstance(cl_val, (int, float)):
                        cl_val = cl_val[0]
                    assert 0 <= cl_val <= 255

            assert reads_with_cl > 0
            assert reads_with_cl / total_reads > 0.5

    def test_final_bam_has_pt_tag(self, final_bam_path):
        """Final BAM reads should have the pt (adapter positions) tag."""
        with pysam.AlignmentFile(str(final_bam_path), "rb") as bam:
            reads_with_pt = 0
            for read in bam.fetch():
                if read.has_tag("pt"):
                    reads_with_pt += 1
                    pt_val = read.get_tag("pt")
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

    def test_final_bam_preserves_modbase_tags(self, final_bam_path):
        """Dorado's MM/ML modbase tags must survive charging classification.

        `escpod classify` adds `cl` alongside them rather than writing
        into ML/MM the way Remora did, and modkit reads MM/ML downstream. A
        regression here would silently empty the modification pileups.
        """
        with pysam.AlignmentFile(str(final_bam_path), "rb") as bam:
            reads_with_mods = sum(
                1 for read in bam.fetch() if read.has_tag("MM") and read.has_tag("ML")
            )
        assert reads_with_mods > 0

    def test_final_bam_has_no_cm_tag(self, final_bam_path):
        """The cm tag is retired along with the Remora ML/MM round-trip."""
        with pysam.AlignmentFile(str(final_bam_path), "rb") as bam:
            assert not any(read.has_tag("cm") for read in bam.fetch())

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
            / "sample1.align_stats.tsv.gz"
        )
        if not path.exists():
            pytest.skip("Alignment stats not found")
        return path

    def test_align_stats_exists(self, align_stats_path):
        """Alignment stats file should exist and have content."""
        assert align_stats_path.exists()
        df = pd.read_csv(align_stats_path, sep="\t")
        assert len(df) > 0

    def test_align_stats_has_key_metrics(self, align_stats_path):
        """Alignment stats should contain expected columns."""
        df = pd.read_csv(align_stats_path, sep="\t")

        expected_columns = [
            "n_reads",
            "pct_mapped",
            "mapped_reads",
        ]

        for col in expected_columns:
            assert col in df.columns, f"Missing column: {col}"


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


@pytest.mark.skipif(
    not (TEST_OUTPUTS_DIR / "summary" / "tables").exists(),
    reason="Pre-computed outputs not available",
)
class TestBaseCallingErrorOutput:
    @pytest.fixture
    def bcerror_path(self):
        path = TEST_OUTPUTS_DIR / "summary" / "tables" / "sample1" / "sample1.bcerror.tsv.gz"
        if not path.exists():
            pytest.skip("bcerror output not found")
        return path

    def test_bcerror_columns(self, bcerror_path):
        df = pd.read_csv(bcerror_path, sep="\t")
        for col in ["Reference", "Position", "Spanning_Reads", "MismatchFreq",
                     "InsertionFreq", "DeletionFreq", "BCErrorFreq"]:
            assert col in df.columns, f"Missing column: {col}"

    def test_freq_values_in_range(self, bcerror_path):
        df = pd.read_csv(bcerror_path, sep="\t")
        for col in ["MismatchFreq", "InsertionFreq", "DeletionFreq", "BCErrorFreq"]:
            assert df[col].min() >= 0, f"{col} has negative values"
            assert df[col].max() <= 1, f"{col} has values > 1"

    def test_positions_1_indexed(self, bcerror_path):
        df = pd.read_csv(bcerror_path, sep="\t")
        assert df["Position"].min() >= 1


@pytest.mark.skipif(
    not (TEST_OUTPUTS_DIR / "summary" / "modkit").exists(),
    reason="Pre-computed outputs not available",
)
class TestModkitOutputs:
    @pytest.fixture
    def pileup_path(self):
        path = TEST_OUTPUTS_DIR / "summary" / "modkit" / "sample1" / "sample1.pileup.bed.gz"
        if not path.exists():
            pytest.skip("Pileup bed not found")
        return path

    @pytest.fixture
    def mod_calls_path(self):
        path = TEST_OUTPUTS_DIR / "summary" / "modkit" / "sample1" / "sample1.mod_calls.tsv.gz"
        if not path.exists():
            pytest.skip("mod_calls not found")
        return path

    @pytest.fixture
    def mod_full_path(self):
        path = TEST_OUTPUTS_DIR / "summary" / "modkit" / "sample1" / "sample1.mod_full.tsv.gz"
        if not path.exists():
            pytest.skip("mod_full not found")
        return path

    def test_pileup_bedmethyl_format(self, pileup_path):
        df = pd.read_csv(pileup_path, sep="\t", header=None)
        assert df.shape[1] >= 10, "bedMethyl should have at least 10 columns"
        assert len(df) > 0

    def test_mod_calls_has_required_columns(self, mod_calls_path):
        df = pd.read_csv(mod_calls_path, sep="\t", nrows=5)
        for col in ["read_id", "ref_position", "chrom", "call_code"]:
            assert col in df.columns, f"Missing column: {col}"

    def test_mod_full_has_header(self, mod_full_path):
        df = pd.read_csv(mod_full_path, sep="\t", nrows=5)
        assert "read_id" in df.columns
        assert len(df) > 0


@pytest.mark.skipif(
    not (TEST_OUTPUTS_DIR / "summary" / "tables").exists(),
    reason="Pre-computed outputs not available",
)
class TestCoverageBedgraphOutput:
    @pytest.fixture
    def counts_path(self):
        path = TEST_OUTPUTS_DIR / "summary" / "tables" / "sample1" / "sample1.counts.bg.gz"
        if not path.exists():
            pytest.skip("counts bedgraph not found")
        return path

    @pytest.fixture
    def cpm_path(self):
        path = TEST_OUTPUTS_DIR / "summary" / "tables" / "sample1" / "sample1.cpm.bg.gz"
        if not path.exists():
            pytest.skip("cpm bedgraph not found")
        return path

    def test_bedgraph_4_columns(self, counts_path):
        df = pd.read_csv(counts_path, sep="\t", header=None)
        assert df.shape[1] == 4, "bedGraph should have exactly 4 columns"

    def test_values_non_negative(self, counts_path, cpm_path):
        for path in [counts_path, cpm_path]:
            df = pd.read_csv(path, sep="\t", header=None)
            assert df.iloc[:, 3].min() >= 0

    def test_positions_0_based(self, counts_path):
        df = pd.read_csv(counts_path, sep="\t", header=None)
        assert df.iloc[:, 1].min() >= 0

    def test_same_chroms_in_both(self, counts_path, cpm_path):
        counts_df = pd.read_csv(counts_path, sep="\t", header=None)
        cpm_df = pd.read_csv(cpm_path, sep="\t", header=None)
        assert set(counts_df.iloc[:, 0]) == set(cpm_df.iloc[:, 0])


@pytest.mark.skipif(
    not (TEST_OUTPUTS_DIR / "manifest.json").exists(),
    reason="Pre-computed outputs not available",
)
class TestManifestOutput:
    @pytest.fixture
    def manifest(self):
        path = TEST_OUTPUTS_DIR / "manifest.json"
        with open(path) as f:
            return json.load(f)

    def test_valid_json(self, manifest):
        assert isinstance(manifest, dict)

    def test_required_top_level_keys(self, manifest):
        for key in ["manifest_version", "pipeline", "execution", "config", "samples", "tools"]:
            assert key in manifest, f"Missing key: {key}"

    def test_status_success(self, manifest):
        assert manifest["execution"]["status"] == "success"

    def test_sample_count(self, manifest):
        assert manifest["samples"]["count"] == 2


@pytest.mark.skipif(
    not (TEST_OUTPUTS_DIR / "squiggy-session.json").exists(),
    reason="Pre-computed outputs not available",
)
class TestSquiggySessionOutput:
    @pytest.fixture
    def session(self):
        path = TEST_OUTPUTS_DIR / "squiggy-session.json"
        with open(path) as f:
            return json.load(f)

    def test_valid_json(self, session):
        assert isinstance(session, dict)

    def test_has_samples(self, session):
        samples = session.get("samples", {})
        assert "sample1" in samples
        assert "sample2" in samples

    def test_sample_has_paths(self, session):
        for sample_name in ["sample1", "sample2"]:
            sample = session["samples"][sample_name]
            assert "bamPath" in sample
            assert "pod5Paths" in sample
            assert "fastaPath" in sample


@pytest.mark.skipif(
    not (TEST_OUTPUTS_DIR / "summary" / "qc" / "reference_similarity.tsv").exists(),
    reason="Pre-computed outputs not available",
)
class TestReferenceSimilarityOutput:
    @pytest.fixture
    def sim_matrix(self):
        path = TEST_OUTPUTS_DIR / "summary" / "qc" / "reference_similarity.tsv"
        return pd.read_csv(path, sep="\t", index_col=0)

    def test_square_matrix(self, sim_matrix):
        assert sim_matrix.shape[0] == sim_matrix.shape[1]

    def test_diagonal_100(self, sim_matrix):
        for i in range(sim_matrix.shape[0]):
            assert sim_matrix.iloc[i, i] == pytest.approx(100.0)

    def test_symmetric(self, sim_matrix):
        for i in range(sim_matrix.shape[0]):
            for j in range(i + 1, sim_matrix.shape[1]):
                assert sim_matrix.iloc[i, j] == pytest.approx(sim_matrix.iloc[j, i], abs=0.01)

    def test_values_in_range(self, sim_matrix):
        assert sim_matrix.min().min() >= 0
        assert sim_matrix.max().max() <= 100


@pytest.mark.skipif(
    not (TEST_OUTPUTS_DIR / "reference" / "trna_only.fa").exists(),
    reason="Pre-computed outputs not available",
)
class TestTrnaOnlyReference:
    @pytest.fixture
    def trna_only_path(self):
        return TEST_OUTPUTS_DIR / "reference" / "trna_only.fa"

    def test_valid_fasta(self, trna_only_path):
        with open(trna_only_path) as f:
            content = f.read()
        assert content.startswith(">")
        seq_lines = [l for l in content.strip().split("\n") if not l.startswith(">")]
        assert len(seq_lines) > 0

    def test_no_adapter_substrings(self, trna_only_path):
        adapter_5p = "CCTAAGAGCAAGAAGAAGCCTGG"
        adapter_3p_prefix = "GGCTTCTTCTTGCTCTT"
        with open(trna_only_path) as f:
            for line in f:
                if line.startswith(">"):
                    continue
                seq = line.strip().upper()
                assert adapter_5p not in seq, "Found 5' adapter in tRNA-only reference"
                assert adapter_3p_prefix not in seq, "Found 3' adapter in tRNA-only reference"

    def test_sequences_end_with_cca(self, trna_only_path):
        sequences = {}
        current = None
        with open(trna_only_path) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    current = line[1:].split()[0]
                    sequences[current] = ""
                else:
                    sequences[current] += line.upper()
        for name, seq in sequences.items():
            assert seq.endswith("CCA"), f"{name} does not end with CCA: ...{seq[-5:]}"


@pytest.mark.skipif(
    not (TEST_OUTPUTS_DIR / "bam" / "final").exists()
    or not (TEST_OUTPUTS_DIR / "summary" / "tables").exists(),
    reason="Pre-computed outputs not available",
)
class TestMultiSampleConsistency:
    SAMPLES = ["sample1", "sample2"]

    def test_both_samples_have_final_bam(self):
        for sample in self.SAMPLES:
            path = TEST_OUTPUTS_DIR / "bam" / "final" / sample / f"{sample}.bam"
            assert path.exists(), f"Missing final BAM for {sample}"

    def test_both_samples_have_charging_tables(self):
        for sample in self.SAMPLES:
            path = (
                TEST_OUTPUTS_DIR / "summary" / "tables" / sample
                / f"{sample}.charging_prob.tsv.gz"
            )
            assert path.exists(), f"Missing charging table for {sample}"

    def test_both_samples_have_bcerror(self):
        for sample in self.SAMPLES:
            path = (
                TEST_OUTPUTS_DIR / "summary" / "tables" / sample
                / f"{sample}.bcerror.tsv.gz"
            )
            assert path.exists(), f"Missing bcerror for {sample}"

    def test_both_samples_have_modkit(self):
        for sample in self.SAMPLES:
            path = (
                TEST_OUTPUTS_DIR / "summary" / "modkit" / sample
                / f"{sample}.mod_calls.tsv.gz"
            )
            assert path.exists(), f"Missing modkit output for {sample}"

    def test_column_names_match_between_samples(self):
        for table in ["bcerror", "charging_prob", "charging.cpm", "align_stats"]:
            dfs = {}
            for sample in self.SAMPLES:
                path = (
                    TEST_OUTPUTS_DIR / "summary" / "tables" / sample
                    / f"{sample}.{table}.tsv.gz"
                )
                if not path.exists():
                    pytest.skip(f"{table} not found for {sample}")
                dfs[sample] = pd.read_csv(path, sep="\t", nrows=0)

            cols1 = set(dfs["sample1"].columns)
            cols2 = set(dfs["sample2"].columns)
            assert cols1 == cols2, f"Column mismatch in {table}: {cols1.symmetric_difference(cols2)}"
