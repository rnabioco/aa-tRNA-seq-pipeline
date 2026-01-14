"""
Shared pytest fixtures for aa-tRNA-seq pipeline tests.
"""

import os
import sys
import tempfile
from array import array
from pathlib import Path

import pysam
import pytest

# Add workflow/scripts to path for imports
REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT / "workflow" / "scripts"))

# Test data paths
TEST_DATA_DIR = REPO_ROOT / ".tests"
TEST_OUTPUTS_DIR = TEST_DATA_DIR / "outputs"
TEST_INPUTS_DIR = TEST_DATA_DIR / "inputs"

# Default adapter sequences (matching add_adapter_tags.py)
ADAPTER_5P = "CCTAAGAGCAAGAAGAAGCCTGG"
ADAPTER_3P = "GGCTTCTTCTTGCTCTTCCAACCTTGCCTTAAAAAAAAAA"


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test outputs."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_bam(temp_dir):
    """Create a minimal BAM file with a few reads for testing."""
    bam_path = temp_dir / "test.bam"

    # Create a minimal header
    header = {
        "HD": {"VN": "1.0"},
        "SQ": [
            {"SN": "tRNA-Ala-AGC-1-1", "LN": 140},
            {"SN": "tRNA-Gly-GCC-1-1", "LN": 135},
        ],
    }

    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as outf:
        # Create a mapped read
        read1 = pysam.AlignedSegment()
        read1.query_name = "read1"
        read1.query_sequence = ADAPTER_5P + "GCGGCTATAGCTCAGTTGGTA" + "CCA" + ADAPTER_3P
        read1.flag = 0
        read1.reference_id = 0
        read1.reference_start = 0
        read1.mapping_quality = 60
        read1.cigartuples = [(0, len(read1.query_sequence))]
        read1.query_qualities = pysam.qualitystring_to_array(
            "I" * len(read1.query_sequence)
        )
        outf.write(read1)

        # Create a read with truncated 5' adapter
        read2 = pysam.AlignedSegment()
        read2.query_name = "read2"
        read2.query_sequence = "AAGCCTGG" + "GCGGCTATAGCTCAGTTGGTA" + "CCA" + ADAPTER_3P
        read2.flag = 0
        read2.reference_id = 0
        read2.reference_start = 15  # Starts after some adapter
        read2.mapping_quality = 60
        read2.cigartuples = [(0, len(read2.query_sequence))]
        read2.query_qualities = pysam.qualitystring_to_array(
            "I" * len(read2.query_sequence)
        )
        outf.write(read2)

        # Create an unmapped read
        read3 = pysam.AlignedSegment()
        read3.query_name = "read3"
        read3.query_sequence = "ACGTACGTACGT"
        read3.flag = 4  # Unmapped
        read3.reference_id = -1
        read3.reference_start = -1
        read3.mapping_quality = 0
        read3.query_qualities = pysam.qualitystring_to_array(
            "I" * len(read3.query_sequence)
        )
        outf.write(read3)

    # Index the BAM
    pysam.index(str(bam_path))

    return bam_path


@pytest.fixture
def sample_bam_with_ml_tag(temp_dir):
    """Create a BAM file with ML tags for charging classification testing."""
    bam_path = temp_dir / "test_ml.bam"

    header = {
        "HD": {"VN": "1.0"},
        "SQ": [
            {"SN": "tRNA-Ala-AGC-1-1", "LN": 140},
            {"SN": "tRNA-Gly-GCC-1-1", "LN": 135},
        ],
    }

    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as outf:
        # Charged read (ML >= 200)
        read1 = pysam.AlignedSegment()
        read1.query_name = "charged_read"
        read1.query_sequence = "A" * 100
        read1.flag = 0
        read1.reference_id = 0
        read1.reference_start = 0
        read1.mapping_quality = 60
        read1.cigartuples = [(0, 100)]
        read1.query_qualities = pysam.qualitystring_to_array("I" * 100)
        read1.set_tag("ML", array("B", [220]))
        outf.write(read1)

        # Uncharged read (ML < 200)
        read2 = pysam.AlignedSegment()
        read2.query_name = "uncharged_read"
        read2.query_sequence = "A" * 100
        read2.flag = 0
        read2.reference_id = 0
        read2.reference_start = 0
        read2.mapping_quality = 60
        read2.cigartuples = [(0, 100)]
        read2.query_qualities = pysam.qualitystring_to_array("I" * 100)
        read2.set_tag("ML", array("B", [150]))
        outf.write(read2)

        # Read on different tRNA
        read3 = pysam.AlignedSegment()
        read3.query_name = "gly_charged"
        read3.query_sequence = "A" * 100
        read3.flag = 0
        read3.reference_id = 1
        read3.reference_start = 0
        read3.mapping_quality = 60
        read3.cigartuples = [(0, 100)]
        read3.query_qualities = pysam.qualitystring_to_array("I" * 100)
        read3.set_tag("ML", array("B", [255]))
        outf.write(read3)

    pysam.index(str(bam_path))
    return bam_path


@pytest.fixture
def sample_fasta(temp_dir):
    """Create a sample FASTA file with adapted tRNA sequences."""
    fasta_path = temp_dir / "test.fa"

    # Create properly adapted tRNA sequence
    trna_seq = "GCGGCTATAGCTCAGTTGGTAGAGCGCTTGCTTAGCATGCAAGAGGTCAGCGGTTCGATCCCGCTATAGCCGCCA"
    adapted_seq = ADAPTER_5P + "G" + trna_seq + ADAPTER_3P

    with open(fasta_path, "w") as f:
        f.write(">tRNA-Ala-AGC-1-1\n")
        f.write(adapted_seq + "\n")
        f.write(">tRNA-Gly-GCC-1-1\n")
        f.write(ADAPTER_5P + "G" + "GCATGCATGCATGCATGCATCCA" + ADAPTER_3P + "\n")

    return fasta_path


@pytest.fixture
def raw_trna_fasta(temp_dir):
    """Create a sample FASTA file with raw (unadapted) tRNA sequences."""
    fasta_path = temp_dir / "raw_trnas.fa"

    with open(fasta_path, "w") as f:
        # tRNA with CCA
        f.write(">tRNA-Ala-AGC-1\n")
        f.write("GCGGCTATAGCTCAGTTGGTACCA\n")
        # tRNA without CCA (should be added)
        f.write(">tRNA-Gly-GCC-1\n")
        f.write("GCATGCATGCATGCATGCAT\n")

    return fasta_path


@pytest.fixture
def test_outputs_available():
    """Check if pre-computed test outputs are available."""
    return TEST_OUTPUTS_DIR.exists() and any(TEST_OUTPUTS_DIR.iterdir())


@pytest.fixture
def test_final_bam():
    """Path to pre-computed final BAM with all tags."""
    bam_path = TEST_OUTPUTS_DIR / "bam" / "final" / "sample1.bam"
    if bam_path.exists():
        return bam_path
    pytest.skip("Pre-computed test outputs not available. Run 'pixi run dl-test-data'")


@pytest.fixture
def test_charging_table():
    """Path to pre-computed charging probability table."""
    table_path = (
        TEST_OUTPUTS_DIR
        / "summary"
        / "tables"
        / "sample1"
        / "sample1.charging_prob.tsv.gz"
    )
    if table_path.exists():
        return table_path
    pytest.skip("Pre-computed test outputs not available. Run 'pixi run dl-test-data'")


# Utility functions for tests


def create_bam_with_reads(path, reads, header=None):
    """
    Helper to create a BAM file with specified reads.

    Args:
        path: Output BAM path
        reads: List of dicts with read attributes
        header: Optional header dict, or default will be used
    """
    if header is None:
        header = {
            "HD": {"VN": "1.0"},
            "SQ": [{"SN": "ref1", "LN": 1000}],
        }

    with pysam.AlignmentFile(str(path), "wb", header=header) as outf:
        for read_dict in reads:
            read = pysam.AlignedSegment()
            read.query_name = read_dict.get("name", "read")
            read.query_sequence = read_dict.get("seq", "ACGT")
            read.flag = read_dict.get("flag", 0)
            read.reference_id = read_dict.get("ref_id", 0)
            read.reference_start = read_dict.get("ref_start", 0)
            read.mapping_quality = read_dict.get("mapq", 60)
            read.cigartuples = read_dict.get(
                "cigar", [(0, len(read.query_sequence))]
            )
            read.query_qualities = pysam.qualitystring_to_array(
                "I" * len(read.query_sequence)
            )

            # Add any tags
            for tag, value in read_dict.get("tags", {}).items():
                if isinstance(value, list):
                    read.set_tag(tag, value, "B")
                else:
                    read.set_tag(tag, value)

            outf.write(read)

    pysam.index(str(path))
    return path
