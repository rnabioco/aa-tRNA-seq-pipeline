"""
Unit tests for anchor_coverage.py

The anchor is located per reference from the adapter prefix shared by every
3' adapter, so contigs carrying a different EDX adapter are all found.
"""

import os

import pysam
import pytest

from anchor_coverage import locate_anchors

EDX01 = "GGCTTCTTCTTGCTCTTCCAACCTTGCCTTAAAAAAAAAA"
EDX02 = "GGCTTCTTCTTGCTCTTATGGAAGGTAGGCAAAAAAAAAA"
COMMON = os.path.commonprefix([EDX01, EDX02])
BODY = "GCCGAGGTAGCTCAGTTGGTAGAGCATGGGACTGAAAATCCCATGGGTCGGGGGTTCGAATCCCTCCCTCGGCACCA"


@pytest.fixture
def dual_adapter_fasta(tmp_path):
    path = tmp_path / "ref.fa"
    path.write_text(f">tRNA-A-edx01\n{BODY}{EDX01}\n>tRNA-A-edx02\n{BODY}{EDX02}\n")
    pysam.faidx(str(path))
    return path


def test_common_prefix_is_shared_stem():
    assert COMMON == "GGCTTCTTCTTGCTCTT"


def test_first_adapter_prefix_misses_other_adapter(dual_adapter_fasta):
    spans = locate_anchors(dual_adapter_fasta, EDX01[:25])
    assert set(spans) == {"tRNA-A-edx01"}


def test_common_prefix_finds_every_contig(dual_adapter_fasta):
    spans = locate_anchors(dual_adapter_fasta, COMMON)
    assert spans == {
        "tRNA-A-edx01": (len(BODY) - 3, len(BODY)),
        "tRNA-A-edx02": (len(BODY) - 3, len(BODY)),
    }
    assert BODY[-3:] == "CCA"
