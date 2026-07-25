"""Unit tests for get_mismatch_calls.py."""

import gzip

import pandas as pd
import pysam
import pytest

from conftest import create_bam_with_reads
from get_mismatch_calls import load_sites, write_mismatch_calls


@pytest.fixture
def simple_fasta(temp_dir):
    """Create a simple indexed FASTA."""
    fa_path = temp_dir / "ref.fa"
    fa_path.write_text(">ref1\nACGTACGT\n")
    pysam.faidx(str(fa_path))
    return fa_path


@pytest.fixture
def adapted_fasta(temp_dir):
    """Create FASTA with adapter regions for trimming tests."""
    fa_path = temp_dir / "ref.fa"
    # 4bp 5' adapter + 8bp tRNA + 4bp 3' adapter = 16bp total
    fa_path.write_text(">ref1\nAAAAACGTACGTTTTT\n")
    pysam.faidx(str(fa_path))
    return fa_path


def read_calls(path):
    with gzip.open(path, "rt") as fh:
        return pd.read_csv(fh, sep="\t")


def run_calls(temp_dir, fasta, reads, ref_len=8, **kwargs):
    """Run the per-read output, keeping matches so position-level tests read clearly."""
    kwargs.setdefault("include_matches", True)
    bam_path = temp_dir / "test.bam"
    out_path = temp_dir / "calls.tsv.gz"
    header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "ref1", "LN": ref_len}]}
    create_bam_with_reads(bam_path, reads, header=header)
    write_mismatch_calls(str(bam_path), str(fasta), str(out_path), **kwargs)
    return read_calls(out_path)


def run_counts(temp_dir, fasta, reads, ref_len=8, **kwargs):
    """Run the per-site error-by-charging counts output."""
    bam_path = temp_dir / "test.bam"
    out_path = temp_dir / "calls.tsv.gz"
    counts_path = temp_dir / "counts.tsv.gz"
    header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "ref1", "LN": ref_len}]}
    create_bam_with_reads(bam_path, reads, header=header)
    write_mismatch_calls(
        str(bam_path), str(fasta), str(out_path), counts_tsv=str(counts_path), **kwargs
    )
    return read_calls(counts_path)


class TestWriteMismatchCalls:
    def test_perfect_alignment_is_all_match(self, temp_dir, simple_fasta):
        """A read matching the reference gets a '-' at every position."""
        reads = [
            {"name": "r1", "seq": "ACGTACGT", "ref_start": 0, "cigar": [(0, 8)]},
        ]
        df = run_calls(temp_dir, simple_fasta, reads)

        assert len(df) == 8
        assert set(df["call_code"]) == {"-"}
        assert list(df["ref_position"]) == list(range(1, 9))
        assert set(df["read_id"]) == {"r1"}
        assert set(df["chrom"]) == {"ref1"}

    def test_single_mismatch(self, temp_dir, simple_fasta):
        """One substituted base is called X at that position only."""
        reads = [
            {"name": "r1", "seq": "AGGTACGT", "ref_start": 0, "cigar": [(0, 8)]},
        ]
        df = run_calls(temp_dir, simple_fasta, reads)

        mismatched = df[df["call_code"] == "X"]
        assert list(mismatched["ref_position"]) == [2]
        assert len(df[df["call_code"] == "-"]) == 7

    def test_deletion(self, temp_dir, simple_fasta):
        """Deleted reference positions are called D."""
        reads = [
            {
                "name": "r1",
                "seq": "ACGTCGT",
                "ref_start": 0,
                "cigar": [(0, 4), (2, 1), (0, 3)],
            },
        ]
        df = run_calls(temp_dir, simple_fasta, reads)

        deleted = df[df["call_code"] == "D"]
        assert list(deleted["ref_position"]) == [5]

    def test_insertion(self, temp_dir, simple_fasta):
        """An insertion is called I at the following reference position."""
        reads = [
            {
                "name": "r1",
                "seq": "ACGTNNACGT",
                "ref_start": 0,
                "cigar": [(0, 4), (1, 2), (0, 4)],
            },
        ]
        df = run_calls(temp_dir, simple_fasta, reads)

        inserted = df[df["call_code"] == "I"]
        assert list(inserted["ref_position"]) == [5]

    def test_reads_are_reported_separately(self, temp_dir, simple_fasta):
        """Each read gets its own row per position, which is the point."""
        reads = [
            {"name": "r1", "seq": "ACGTACGT", "ref_start": 0, "cigar": [(0, 8)]},
            {"name": "r2", "seq": "AGGTACGT", "ref_start": 0, "cigar": [(0, 8)]},
        ]
        df = run_calls(temp_dir, simple_fasta, reads)

        assert len(df) == 16
        at_two = df[df["ref_position"] == 2].set_index("read_id")["call_code"]
        assert at_two["r1"] == "-"
        assert at_two["r2"] == "X"

    def test_partial_alignment_omits_uncovered_positions(
        self, temp_dir, simple_fasta
    ):
        """Positions the read does not span are absent, not reported as matches."""
        reads = [
            {"name": "r1", "seq": "ACGT", "ref_start": 2, "cigar": [(0, 4)]},
        ]
        df = run_calls(temp_dir, simple_fasta, reads)

        assert list(df["ref_position"]) == [3, 4, 5, 6]

    def test_unmapped_and_reverse_reads_are_skipped(self, temp_dir, simple_fasta):
        """Only forward, mapped alignments contribute, as in get_bcerror_freqs."""
        reads = [
            {"name": "fwd", "seq": "ACGTACGT", "ref_start": 0, "cigar": [(0, 8)]},
            {
                "name": "rev",
                "seq": "ACGTACGT",
                "flag": 16,
                "ref_start": 0,
                "cigar": [(0, 8)],
            },
        ]
        df = run_calls(temp_dir, simple_fasta, reads)

        assert set(df["read_id"]) == {"fwd"}

    def test_sites_restrict_output(self, temp_dir, simple_fasta):
        """A site list keeps only the requested positions."""
        sites_path = temp_dir / "sites.tsv"
        sites_path.write_text("ref\tpos\nref1\t2\nref1\t5\n")

        reads = [
            {"name": "r1", "seq": "AGGTACGT", "ref_start": 0, "cigar": [(0, 8)]},
        ]
        df = run_calls(
            temp_dir, simple_fasta, reads, sites=load_sites(str(sites_path))
        )

        assert list(df["ref_position"]) == [2, 5]
        assert list(df["call_code"]) == ["X", "-"]

    def test_sites_skip_absent_references(self, temp_dir, simple_fasta):
        """References with no requested site produce no rows."""
        sites_path = temp_dir / "sites.tsv"
        sites_path.write_text("ref\tpos\nother-ref\t2\n")

        reads = [
            {"name": "r1", "seq": "ACGTACGT", "ref_start": 0, "cigar": [(0, 8)]},
        ]
        df = run_calls(
            temp_dir, simple_fasta, reads, sites=load_sites(str(sites_path))
        )

        assert len(df) == 0

    def test_offsets_convert_to_trna_coordinates(self, temp_dir, adapted_fasta):
        """Adapter positions are dropped and the rest renumbered from 1."""
        reads = [
            {
                "name": "r1",
                "seq": "AAAAACGTACGTTTTT",
                "ref_start": 0,
                "cigar": [(0, 16)],
            },
        ]
        df = run_calls(
            temp_dir,
            adapted_fasta,
            reads,
            ref_len=16,
            trim_5p=4,
            trim_3p=4,
        )

        assert list(df["ref_position"]) == list(range(1, 9))
        assert set(df["call_code"]) == {"-"}

    def test_offsets_apply_to_sites(self, temp_dir, adapted_fasta):
        """Site positions are interpreted in tRNA-only coordinates."""
        sites_path = temp_dir / "sites.tsv"
        sites_path.write_text("ref\tpos\nref1\t1\n")

        reads = [
            {
                "name": "r1",
                "seq": "AAAAGCGTACGTTTTT",
                "ref_start": 0,
                "cigar": [(0, 16)],
            },
        ]
        df = run_calls(
            temp_dir,
            adapted_fasta,
            reads,
            ref_len=16,
            sites=load_sites(str(sites_path)),
            trim_5p=4,
            trim_3p=4,
        )

        # tRNA position 1 is reference position 5, which is substituted.
        assert list(df["ref_position"]) == [1]
        assert list(df["call_code"]) == ["X"]

    def test_output_schema_matches_modkit(self, temp_dir, simple_fasta):
        """Downstream code reads modkit and mismatch calls interchangeably."""
        reads = [
            {"name": "r1", "seq": "ACGTACGT", "ref_start": 0, "cigar": [(0, 8)]},
        ]
        df = run_calls(temp_dir, simple_fasta, reads)

        assert list(df.columns) == [
            "read_id",
            "ref_position",
            "chrom",
            "within_alignment",
            "call_code",
        ]
        assert set(df["within_alignment"]) == {True}


class TestMatchRowsAreOmitted:
    def test_matches_are_dropped_by_default(self, temp_dir, simple_fasta):
        """Match rows are ~70% of positions and carry nothing individually."""
        reads = [
            {"name": "r1", "seq": "AGGTACGT", "ref_start": 0, "cigar": [(0, 8)]},
        ]
        df = run_calls(temp_dir, simple_fasta, reads, include_matches=False)

        assert list(df["call_code"]) == ["X"]
        assert list(df["ref_position"]) == [2]

    def test_include_matches_restores_them(self, temp_dir, simple_fasta):
        reads = [
            {"name": "r1", "seq": "AGGTACGT", "ref_start": 0, "cigar": [(0, 8)]},
        ]
        df = run_calls(temp_dir, simple_fasta, reads, include_matches=True)

        assert len(df) == 8
        assert sorted(df["call_code"].unique()) == ["-", "X"]

    def test_a_read_with_no_errors_emits_nothing(self, temp_dir, simple_fasta):
        reads = [
            {"name": "r1", "seq": "ACGTACGT", "ref_start": 0, "cigar": [(0, 8)]},
        ]
        df = run_calls(temp_dir, simple_fasta, reads, include_matches=False)

        assert len(df) == 0


class TestChargingErrorCounts:
    def test_counts_fill_the_2x2(self, temp_dir, simple_fasta):
        """Every cell of each site's contingency table is recorded."""
        reads = [
            # charged, error at position 2
            {"name": "r1", "seq": "AGGTACGT", "ref_start": 0, "cigar": [(0, 8)],
             "tags": {"cl": 240}},
            # charged, no error
            {"name": "r2", "seq": "ACGTACGT", "ref_start": 0, "cigar": [(0, 8)],
             "tags": {"cl": 240}},
            # uncharged, error at position 2
            {"name": "r3", "seq": "AGGTACGT", "ref_start": 0, "cigar": [(0, 8)],
             "tags": {"cl": 40}},
            # uncharged, no error
            {"name": "r4", "seq": "ACGTACGT", "ref_start": 0, "cigar": [(0, 8)],
             "tags": {"cl": 40}},
        ]
        df = run_counts(temp_dir, simple_fasta, reads)
        at_two = df[df["pos"] == 2].iloc[0]

        assert at_two["err_charged"] == 1
        assert at_two["err_uncharged"] == 1
        assert at_two["match_charged"] == 1
        assert at_two["match_uncharged"] == 1

        # A position where nobody errs still records its coverage.
        at_one = df[df["pos"] == 1].iloc[0]
        assert at_one["err_charged"] == 0
        assert at_one["err_uncharged"] == 0
        assert at_one["match_charged"] == 2
        assert at_one["match_uncharged"] == 2

    def test_ml_threshold_splits_charging(self, temp_dir, simple_fasta):
        reads = [
            {"name": "r1", "seq": "AGGTACGT", "ref_start": 0, "cigar": [(0, 8)],
             "tags": {"cl": 200}},
        ]

        at_default = run_counts(temp_dir, simple_fasta, reads)
        assert at_default[at_default["pos"] == 2].iloc[0]["err_charged"] == 1

        at_strict = run_counts(
            temp_dir, simple_fasta, reads, ml_threshold=201
        )
        assert at_strict[at_strict["pos"] == 2].iloc[0]["err_uncharged"] == 1

    def test_zero_tag_counts_as_uncharged(self, temp_dir, simple_fasta):
        """A cl of 0 is a confident uncharged call, not a missing one."""
        reads = [
            {"name": "r1", "seq": "AGGTACGT", "ref_start": 0, "cigar": [(0, 8)],
             "tags": {"cl": 0}},
        ]
        df = run_counts(temp_dir, simple_fasta, reads)

        assert df[df["pos"] == 2].iloc[0]["err_uncharged"] == 1

    def test_reads_without_a_tag_are_excluded(self, temp_dir, simple_fasta):
        reads = [
            {"name": "r1", "seq": "AGGTACGT", "ref_start": 0, "cigar": [(0, 8)],
             "tags": {"cl": 240}},
            {"name": "r2", "seq": "AGGTACGT", "ref_start": 0, "cigar": [(0, 8)]},
        ]
        df = run_counts(temp_dir, simple_fasta, reads)
        at_two = df[df["pos"] == 2].iloc[0]

        assert at_two["err_charged"] == 1
        assert at_two["err_uncharged"] == 0

    def test_counts_respect_sites(self, temp_dir, simple_fasta):
        sites_path = temp_dir / "sites.tsv"
        sites_path.write_text("ref\tpos\nref1\t2\n")

        reads = [
            {"name": "r1", "seq": "AGGTACGT", "ref_start": 0, "cigar": [(0, 8)],
             "tags": {"cl": 240}},
        ]
        df = run_counts(
            temp_dir, simple_fasta, reads, sites=load_sites(str(sites_path))
        )

        assert list(df["pos"]) == [2]

    def test_counts_totals_match_read_depth(self, temp_dir, simple_fasta):
        """The four cells must sum to the number of tagged reads covering a site."""
        reads = [
            {"name": f"r{i}", "seq": "ACGTACGT", "ref_start": 0, "cigar": [(0, 8)],
             "tags": {"cl": 240 if i % 2 else 40}}
            for i in range(6)
        ]
        df = run_counts(temp_dir, simple_fasta, reads)
        totals = (
            df["err_charged"] + df["err_uncharged"]
            + df["match_charged"] + df["match_uncharged"]
        )

        assert set(totals) == {6}


class TestLoadSites:
    def test_returns_none_without_a_path(self):
        assert load_sites(None) is None

    def test_accepts_bcerror_column_names(self, temp_dir):
        """Site lists may come straight from a bcerror table."""
        path = temp_dir / "sites.tsv"
        path.write_text("Reference\tPosition\nref1\t3\nref1\t7\n")

        assert load_sites(str(path)) == {"ref1": {3, 7}}

    def test_accepts_gzipped_input(self, temp_dir):
        path = temp_dir / "sites.tsv.gz"
        with gzip.open(path, "wt") as fh:
            fh.write("ref\tpos\nref1\t3\n")

        assert load_sites(str(path)) == {"ref1": {3}}

    def test_errors_without_position_columns(self, temp_dir):
        path = temp_dir / "sites.tsv"
        path.write_text("foo\tbar\n1\t2\n")

        with pytest.raises(ValueError, match="reference and position"):
            load_sites(str(path))
