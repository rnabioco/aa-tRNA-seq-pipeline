"""
Unit tests for get_charging_table.py

Tests ML tag extraction from BAM files, and the tie set (`XA`) each read's
count is split across.
"""

import gzip
from array import array

import pandas as pd
import pysam
import pytest

from get_charging_table import extract_tag, reference_weights, tie_references
from get_trna_charging_cpm import per_read_charging


class TestExtractTag:
    """Tests for extract_tag function."""

    def test_extracts_ml_tag(self, temp_dir):
        """Should extract ML tag values to TSV."""
        input_bam = temp_dir / "input.bam"
        output_tsv = temp_dir / "output.tsv"

        header = {
            "HD": {"VN": "1.0"},
            "SQ": [{"SN": "tRNA-Ala-AGC-1-1", "LN": 100}],
        }

        with pysam.AlignmentFile(str(input_bam), "wb", header=header) as outf:
            read = pysam.AlignedSegment()
            read.query_name = "read1"
            read.query_sequence = "A" * 100
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 0
            read.cigartuples = [(0, 100)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            read.set_tag("ML", array("B", [220]))
            outf.write(read)

        pysam.index(str(input_bam))

        extract_tag(str(input_bam), str(output_tsv), "ML")

        # Check output
        with open(output_tsv) as f:
            lines = f.readlines()

        assert len(lines) == 2  # Header + 1 read
        assert "read_id\ttRNA\tcharging_likelihood\ttie_refs\n" == lines[0]
        # a read with no XA has an empty tie set
        assert "read1\ttRNA-Ala-AGC-1-1\t220\t\n" == lines[1]

    def test_handles_gzip_output(self, temp_dir):
        """Should write gzipped output when filename ends in .gz."""
        input_bam = temp_dir / "input.bam"
        output_tsv = temp_dir / "output.tsv.gz"

        header = {
            "HD": {"VN": "1.0"},
            "SQ": [{"SN": "tRNA-Ala-AGC-1-1", "LN": 100}],
        }

        with pysam.AlignmentFile(str(input_bam), "wb", header=header) as outf:
            read = pysam.AlignedSegment()
            read.query_name = "read1"
            read.query_sequence = "A" * 100
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 0
            read.cigartuples = [(0, 100)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            read.set_tag("ML", array("B", [200]))
            outf.write(read)

        pysam.index(str(input_bam))

        extract_tag(str(input_bam), str(output_tsv), "ML")

        # Check gzipped output
        with gzip.open(str(output_tsv), "rt") as f:
            lines = f.readlines()

        assert len(lines) == 2
        assert "read1" in lines[1]

    def test_multiple_reads(self, temp_dir):
        """Should handle multiple reads."""
        input_bam = temp_dir / "input.bam"
        output_tsv = temp_dir / "output.tsv"

        header = {
            "HD": {"VN": "1.0"},
            "SQ": [
                {"SN": "tRNA-Ala-AGC-1-1", "LN": 100},
                {"SN": "tRNA-Gly-GCC-1-1", "LN": 100},
            ],
        }

        with pysam.AlignmentFile(str(input_bam), "wb", header=header) as outf:
            for i, (ref_id, ml_val) in enumerate([(0, 220), (0, 150), (1, 255)]):
                read = pysam.AlignedSegment()
                read.query_name = f"read{i}"
                read.query_sequence = "A" * 100
                read.flag = 0
                read.reference_id = ref_id
                read.reference_start = 0
                read.cigartuples = [(0, 100)]
                read.query_qualities = pysam.qualitystring_to_array("I" * 100)
                read.set_tag("ML", array("B", [ml_val]))
                outf.write(read)

        pysam.index(str(input_bam))

        extract_tag(str(input_bam), str(output_tsv), "ML")

        with open(output_tsv) as f:
            lines = f.readlines()

        assert len(lines) == 4  # Header + 3 reads

    def test_skips_unmapped_reads(self, temp_dir):
        """Unmapped reads should be excluded."""
        input_bam = temp_dir / "input.bam"
        output_tsv = temp_dir / "output.tsv"

        header = {
            "HD": {"VN": "1.0"},
            "SQ": [{"SN": "tRNA-Ala-AGC-1-1", "LN": 100}],
        }

        with pysam.AlignmentFile(str(input_bam), "wb", header=header) as outf:
            # Mapped read
            read1 = pysam.AlignedSegment()
            read1.query_name = "mapped"
            read1.query_sequence = "A" * 100
            read1.flag = 0
            read1.reference_id = 0
            read1.reference_start = 0
            read1.cigartuples = [(0, 100)]
            read1.query_qualities = pysam.qualitystring_to_array("I" * 100)
            read1.set_tag("ML", array("B", [200]))
            outf.write(read1)

            # Unmapped read with ML tag
            read2 = pysam.AlignedSegment()
            read2.query_name = "unmapped"
            read2.query_sequence = "A" * 50
            read2.flag = 4
            read2.reference_id = -1
            read2.reference_start = -1
            read2.query_qualities = pysam.qualitystring_to_array("I" * 50)
            read2.set_tag("ML", array("B", [150]))
            outf.write(read2)

        pysam.index(str(input_bam))

        extract_tag(str(input_bam), str(output_tsv), "ML")

        with open(output_tsv) as f:
            lines = f.readlines()

        # Only mapped read should appear (plus header)
        assert len(lines) == 2
        assert "mapped" in lines[1]

    def test_retains_ml_zero_reads(self, temp_dir):
        """ML==0 is a valid, maximally-confident *uncharged* call and must be
        written to the table. Regression test for a bug where the write gate
        used `if tag_value` (falsy at 0), silently dropping ML==0 reads and
        biasing charging fraction upward."""
        input_bam = temp_dir / "input.bam"
        output_tsv = temp_dir / "output.tsv"

        header = {
            "HD": {"VN": "1.0"},
            "SQ": [{"SN": "tRNA-Ala-AGC-1-1", "LN": 100}],
        }

        with pysam.AlignmentFile(str(input_bam), "wb", header=header) as outf:
            # One confidently-uncharged read (ML=0) and one charged (ML=255)
            for name, ml_val in [("uncharged_zero", 0), ("charged_max", 255)]:
                read = pysam.AlignedSegment()
                read.query_name = name
                read.query_sequence = "A" * 100
                read.flag = 0
                read.reference_id = 0
                read.reference_start = 0
                read.cigartuples = [(0, 100)]
                read.query_qualities = pysam.qualitystring_to_array("I" * 100)
                read.set_tag("ML", array("B", [ml_val]))
                outf.write(read)

        pysam.index(str(input_bam))

        extract_tag(str(input_bam), str(output_tsv), "ML")

        with open(output_tsv) as f:
            lines = f.readlines()

        # Header + both reads (the ML=0 read must NOT be dropped)
        assert len(lines) == 3
        assert any(
            line.startswith("uncharged_zero\t") and line.rstrip().endswith("\t0")
            for line in lines[1:]
        )

    def test_skips_and_warns_on_multielement_tag(self, temp_dir, capsys):
        """A multi-element tag array (e.g. the dorado mod-base ML tag) is not a
        single charging score and is skipped, but the skip must be reported on
        stderr rather than dropped silently (silent drops bias the charging
        fraction like the ML==0 bug did)."""
        input_bam = temp_dir / "input.bam"
        output_tsv = temp_dir / "output.tsv"

        header = {
            "HD": {"VN": "1.0"},
            "SQ": [{"SN": "tRNA-Ala-AGC-1-1", "LN": 100}],
        }

        with pysam.AlignmentFile(str(input_bam), "wb", header=header) as outf:
            # One scalar-usable read and one multi-element read (2 mod probs)
            scalar_read = pysam.AlignedSegment()
            scalar_read.query_name = "scalar"
            scalar_read.query_sequence = "A" * 100
            scalar_read.flag = 0
            scalar_read.reference_id = 0
            scalar_read.reference_start = 0
            scalar_read.cigartuples = [(0, 100)]
            scalar_read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            scalar_read.set_tag("ML", array("B", [200]))
            outf.write(scalar_read)

            multi_read = pysam.AlignedSegment()
            multi_read.query_name = "multi"
            multi_read.query_sequence = "A" * 100
            multi_read.flag = 0
            multi_read.reference_id = 0
            multi_read.reference_start = 0
            multi_read.cigartuples = [(0, 100)]
            multi_read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            multi_read.set_tag("ML", array("B", [10, 250]))
            outf.write(multi_read)

        pysam.index(str(input_bam))

        extract_tag(str(input_bam), str(output_tsv), "ML")

        with open(output_tsv) as f:
            lines = f.readlines()

        # Header + only the scalar read; the multi-element read is skipped
        assert len(lines) == 2
        assert lines[1].startswith("scalar\t")
        assert not any(line.startswith("multi\t") for line in lines[1:])

        # ...but the skip is reported, not silent
        err = capsys.readouterr().err
        assert "skipped 1" in err

    def test_different_tag(self, temp_dir):
        """Should work with different tag names."""
        input_bam = temp_dir / "input.bam"
        output_tsv = temp_dir / "output.tsv"

        header = {
            "HD": {"VN": "1.0"},
            "SQ": [{"SN": "tRNA-Ala-AGC-1-1", "LN": 100}],
        }

        with pysam.AlignmentFile(str(input_bam), "wb", header=header) as outf:
            read = pysam.AlignedSegment()
            read.query_name = "read1"
            read.query_sequence = "A" * 100
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 0
            read.cigartuples = [(0, 100)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            read.set_tag("CL", array("B", [180]))
            outf.write(read)

        pysam.index(str(input_bam))

        extract_tag(str(input_bam), str(output_tsv), "CL")

        with open(output_tsv) as f:
            lines = f.readlines()

        assert len(lines) == 2
        assert "180" in lines[1]


REFS = ["tRNA-Ala-AGC-1-1", "tRNA-Ala-AGC-2-1", "tRNA-Ala-AGC-3-1", "tRNA-Gly-GCC-1-1"]


def _write_tied_bam(path, reads):
    """reads: (name, ref_index, cl, xa or None). Coordinate-sorted and indexed."""
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": ref, "LN": 100} for ref in REFS],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as outf:
        for name, ref_id, cl, xa in sorted(reads, key=lambda r: r[1]):
            read = pysam.AlignedSegment(outf.header)
            read.query_name = name
            read.query_sequence = "A" * 50
            read.flag = 0
            read.reference_id = ref_id
            read.reference_start = 0
            read.cigartuples = [(0, 50)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 50)
            read.mapping_quality = 0 if xa else 60
            read.set_tag("cl", cl, value_type="C")
            if xa:
                read.set_tag("XA", xa)
            outf.write(read)
    pysam.index(str(path))


def _xa(*refs):
    return "".join(f"{ref},+3,50M,2;" for ref in refs)


class TestTieReferences:
    def _read(self, temp_dir, xa):
        bam = temp_dir / "one.bam"
        _write_tied_bam(bam, [("r", 0, 200, xa)])
        with pysam.AlignmentFile(str(bam)) as fh:
            return next(fh.fetch())

    def test_parses_escpod_xa(self, temp_dir):
        read = self._read(temp_dir, _xa(REFS[1], REFS[2]))
        assert tie_references(read) == [REFS[1], REFS[2]]

    def test_no_xa_is_no_ties(self, temp_dir):
        assert tie_references(self._read(temp_dir, None)) == []

    def test_primary_and_repeats_are_not_extra_references(self, temp_dir):
        """bwa-style XA can list another position on the same reference."""
        read = self._read(temp_dir, _xa(REFS[0], REFS[1], REFS[1]))
        assert tie_references(read) == [REFS[1]]


class TestReferenceWeights:
    @pytest.mark.parametrize("n", [0, 1, 2, 3])
    def test_n_ties_give_one_over_n_plus_one_each(self, n):
        weights = reference_weights(REFS[0], REFS[1 : 1 + n])
        assert set(weights) == set(REFS[: 1 + n])
        for w in weights.values():
            assert w == pytest.approx(1 / (n + 1))
        assert sum(weights.values()) == pytest.approx(1.0)

    def test_untied_read_weighs_one_on_its_reference(self):
        assert reference_weights(REFS[2]) == {REFS[2]: 1.0}


class TestTieSplitCounts:
    """charging_prob (extract_tag) -> charging.cpm (per_read_charging)."""

    def _counts(self, temp_dir, reads):
        bam = temp_dir / "tied.bam"
        prob = temp_dir / "prob.tsv.gz"
        cpm = temp_dir / "cpm.tsv"
        _write_tied_bam(bam, reads)
        extract_tag(str(bam), str(prob), "cl")
        per_read_charging(str(prob), str(cpm), 200)
        return pd.read_csv(prob, sep="\t"), pd.read_csv(cpm, sep="\t", index_col=0)

    def test_tied_read_split_evenly_untied_read_whole(self, temp_dir):
        reads = [
            ("tied3", 0, 250, _xa(REFS[1], REFS[2])),  # charged, 1/3 each
            ("tied2", 1, 10, _xa(REFS[3])),  # uncharged, 1/2 each
            ("unique", 3, 220, None),  # charged, 1 on Gly
        ]
        prob, cpm = self._counts(temp_dir, reads)
        # one row per scored read, primary reference in tRNA
        assert len(prob) == 3
        assert set(prob["read_id"]) == {"tied3", "tied2", "unique"}

        charged = cpm["counts_charged"]
        uncharged = cpm["counts_uncharged"]
        assert charged[REFS[0]] == pytest.approx(1 / 3)
        assert charged[REFS[1]] == pytest.approx(1 / 3)
        assert charged[REFS[2]] == pytest.approx(1 / 3)
        assert charged[REFS[3]] == pytest.approx(1.0)
        assert uncharged[REFS[1]] == pytest.approx(1 / 2)
        assert uncharged[REFS[3]] == pytest.approx(1 / 2)

    def test_totals_sum_to_scored_read_count(self, temp_dir):
        reads = [
            ("a", 0, 250, _xa(REFS[1], REFS[2], REFS[3])),
            ("b", 0, 100, _xa(REFS[1])),
            ("c", 2, 200, None),
            ("d", 3, 0, _xa(REFS[0], REFS[2])),
            ("e", 1, 199, None),
        ]
        prob, cpm = self._counts(temp_dir, reads)
        total = cpm["counts_charged"].sum() + cpm["counts_uncharged"].sum()
        assert total == pytest.approx(len(reads))
        assert len(prob) == len(reads)
        cpm_total = cpm["cpm_charged"].sum() + cpm["cpm_uncharged"].sum()
        assert cpm_total == pytest.approx(1e6)

    def test_table_without_tie_column_counts_each_read_once(self, temp_dir):
        """charging_prob tables written before the tie_refs column."""
        prob = temp_dir / "old.tsv"
        prob.write_text(
            "read_id\ttRNA\tcharging_likelihood\n"
            f"r1\t{REFS[0]}\t250\nr2\t{REFS[0]}\t10\n"
        )
        cpm = temp_dir / "cpm.tsv"
        per_read_charging(str(prob), str(cpm), 200)
        df = pd.read_csv(cpm, sep="\t", index_col=0)
        assert df.loc[REFS[0], "counts_charged"] == 1
        assert df.loc[REFS[0], "counts_uncharged"] == 1
