"""
Unit tests for transfer_tags.py

Tests tag transfer between BAM files and tag renaming.
"""

from array import array

import pysam
import pytest

from transfer_tags import transfer_tags, parse_tag_items, parse_set_tags


class TestParseTagItems:
    """Tests for parse_tag_items function."""

    def test_single_rename(self):
        """Single tag rename should parse correctly."""
        result = parse_tag_items(["ML=CL"])
        assert result == {"ML": "CL"}

    def test_multiple_renames(self):
        """Multiple renames should all be captured."""
        result = parse_tag_items(["ML=CL", "MM=CM"])
        assert result == {"ML": "CL", "MM": "CM"}

    def test_empty_list(self):
        """Empty list should return empty dict."""
        result = parse_tag_items([])
        assert result == {}

    def test_whitespace_handling(self):
        """Whitespace around = should be stripped."""
        result = parse_tag_items(["ML = CL"])
        assert result == {"ML": "CL"}


class TestTransferTags:
    """Tests for transfer_tags function."""

    def test_basic_transfer(self, temp_dir):
        """Tags should be transferred from source to target BAM."""
        source_bam = temp_dir / "source.bam"
        target_bam = temp_dir / "target.bam"
        output_bam = temp_dir / "output.bam"

        header = {
            "HD": {"VN": "1.0"},
            "SQ": [{"SN": "ref", "LN": 100}],
        }

        # Create source BAM with ML tag
        with pysam.AlignmentFile(str(source_bam), "wb", header=header) as outf:
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

        # Create target BAM without ML tag
        with pysam.AlignmentFile(str(target_bam), "wb", header=header) as outf:
            read = pysam.AlignedSegment()
            read.query_name = "read1"
            read.query_sequence = "A" * 100
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 0
            read.cigartuples = [(0, 100)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            read.set_tag("XY", "target_tag")  # Different tag
            outf.write(read)

        transfer_tags(
            tags=["ML"],
            rename=[],
            source_bam=str(source_bam),
            target_bam=str(target_bam),
            output_bam=str(output_bam),
        )

        # Verify output
        with pysam.AlignmentFile(str(output_bam), "rb") as bam:
            read = next(bam)
            assert read.has_tag("ML")
            ml_val = read.get_tag("ML")
            # Single-element arrays are unwrapped to scalar by transfer_tags
            assert ml_val == 220
            # Original target tag should also be present
            assert read.has_tag("XY")

    def test_tag_renaming(self, temp_dir):
        """Tags should be renamed during transfer."""
        source_bam = temp_dir / "source.bam"
        target_bam = temp_dir / "target.bam"
        output_bam = temp_dir / "output.bam"

        header = {
            "HD": {"VN": "1.0"},
            "SQ": [{"SN": "ref", "LN": 100}],
        }

        # Create source with ML tag
        with pysam.AlignmentFile(str(source_bam), "wb", header=header) as outf:
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

        # Create target
        with pysam.AlignmentFile(str(target_bam), "wb", header=header) as outf:
            read = pysam.AlignedSegment()
            read.query_name = "read1"
            read.query_sequence = "A" * 100
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 0
            read.cigartuples = [(0, 100)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            outf.write(read)

        # Transfer with rename ML -> CL
        transfer_tags(
            tags=["ML"],
            rename=["ML=CL"],
            source_bam=str(source_bam),
            target_bam=str(target_bam),
            output_bam=str(output_bam),
        )

        with pysam.AlignmentFile(str(output_bam), "rb") as bam:
            read = next(bam)
            # Should have CL, not ML
            assert read.has_tag("CL")
            assert not read.has_tag("ML")

    def test_skips_secondary_alignments(self, temp_dir):
        """Secondary alignments should be excluded from output."""
        source_bam = temp_dir / "source.bam"
        target_bam = temp_dir / "target.bam"
        output_bam = temp_dir / "output.bam"

        header = {
            "HD": {"VN": "1.0"},
            "SQ": [{"SN": "ref", "LN": 100}],
        }

        # Create source with tags
        with pysam.AlignmentFile(str(source_bam), "wb", header=header) as outf:
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

        # Create target with primary and secondary alignments
        with pysam.AlignmentFile(str(target_bam), "wb", header=header) as outf:
            # Primary alignment
            read1 = pysam.AlignedSegment()
            read1.query_name = "read1"
            read1.query_sequence = "A" * 100
            read1.flag = 0  # Primary
            read1.reference_id = 0
            read1.reference_start = 0
            read1.cigartuples = [(0, 100)]
            read1.query_qualities = pysam.qualitystring_to_array("I" * 100)
            outf.write(read1)

            # Secondary alignment
            read2 = pysam.AlignedSegment()
            read2.query_name = "read1"
            read2.query_sequence = "A" * 100
            read2.flag = 256  # Secondary flag
            read2.reference_id = 0
            read2.reference_start = 10
            read2.cigartuples = [(0, 100)]
            read2.query_qualities = pysam.qualitystring_to_array("I" * 100)
            outf.write(read2)

        transfer_tags(
            tags=["ML"],
            rename=[],
            source_bam=str(source_bam),
            target_bam=str(target_bam),
            output_bam=str(output_bam),
        )

        with pysam.AlignmentFile(str(output_bam), "rb") as bam:
            reads = list(bam)

        # Only primary alignment should be in output
        assert len(reads) == 1
        assert not reads[0].is_secondary

    def test_only_outputs_reads_with_tags(self, temp_dir):
        """Reads without matching tags in source should not be in output."""
        source_bam = temp_dir / "source.bam"
        target_bam = temp_dir / "target.bam"
        output_bam = temp_dir / "output.bam"

        header = {
            "HD": {"VN": "1.0"},
            "SQ": [{"SN": "ref", "LN": 100}],
        }

        # Source only has read1
        with pysam.AlignmentFile(str(source_bam), "wb", header=header) as outf:
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

        # Target has read1 and read2
        with pysam.AlignmentFile(str(target_bam), "wb", header=header) as outf:
            for name in ["read1", "read2"]:
                read = pysam.AlignedSegment()
                read.query_name = name
                read.query_sequence = "A" * 100
                read.flag = 0
                read.reference_id = 0
                read.reference_start = 0
                read.cigartuples = [(0, 100)]
                read.query_qualities = pysam.qualitystring_to_array("I" * 100)
                outf.write(read)

        transfer_tags(
            tags=["ML"],
            rename=[],
            source_bam=str(source_bam),
            target_bam=str(target_bam),
            output_bam=str(output_bam),
        )

        with pysam.AlignmentFile(str(output_bam), "rb") as bam:
            reads = list(bam)

        # Only read1 (which has tags transferred) should be in output
        assert len(reads) == 1
        assert reads[0].query_name == "read1"

    def test_multiple_tags(self, temp_dir):
        """Multiple tags should all be transferred."""
        source_bam = temp_dir / "source.bam"
        target_bam = temp_dir / "target.bam"
        output_bam = temp_dir / "output.bam"

        header = {
            "HD": {"VN": "1.0"},
            "SQ": [{"SN": "ref", "LN": 100}],
        }

        # Source with multiple tags
        with pysam.AlignmentFile(str(source_bam), "wb", header=header) as outf:
            read = pysam.AlignedSegment()
            read.query_name = "read1"
            read.query_sequence = "A" * 100
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 0
            read.cigartuples = [(0, 100)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            read.set_tag("ML", array("B", [220]))
            read.set_tag("MM", "A+a.,0;")
            outf.write(read)

        # Target
        with pysam.AlignmentFile(str(target_bam), "wb", header=header) as outf:
            read = pysam.AlignedSegment()
            read.query_name = "read1"
            read.query_sequence = "A" * 100
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 0
            read.cigartuples = [(0, 100)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            outf.write(read)

        transfer_tags(
            tags=["ML", "MM"],
            rename=[],
            source_bam=str(source_bam),
            target_bam=str(target_bam),
            output_bam=str(output_bam),
        )

        with pysam.AlignmentFile(str(output_bam), "rb") as bam:
            read = next(bam)
            assert read.has_tag("ML")
            assert read.has_tag("MM")

    def test_multithreaded_output_matches(self, temp_dir):
        """Using threads > 1 should produce identical results."""
        source_bam = temp_dir / "source.bam"
        target_bam = temp_dir / "target.bam"
        output_st = temp_dir / "output_st.bam"
        output_mt = temp_dir / "output_mt.bam"

        header = {
            "HD": {"VN": "1.0"},
            "SQ": [{"SN": "ref", "LN": 100}],
        }

        # Create source BAM with tags
        with pysam.AlignmentFile(str(source_bam), "wb", header=header) as outf:
            for i in range(10):
                read = pysam.AlignedSegment()
                read.query_name = f"read{i}"
                read.query_sequence = "A" * 100
                read.flag = 0
                read.reference_id = 0
                read.reference_start = i * 10
                read.cigartuples = [(0, 100)]
                read.query_qualities = pysam.qualitystring_to_array("I" * 100)
                read.set_tag("ML", array("B", [200 + i]))
                outf.write(read)

        # Create target BAM
        with pysam.AlignmentFile(str(target_bam), "wb", header=header) as outf:
            for i in range(10):
                read = pysam.AlignedSegment()
                read.query_name = f"read{i}"
                read.query_sequence = "A" * 100
                read.flag = 0
                read.reference_id = 0
                read.reference_start = i * 10
                read.cigartuples = [(0, 100)]
                read.query_qualities = pysam.qualitystring_to_array("I" * 100)
                outf.write(read)

        # Single-threaded
        transfer_tags(
            tags=["ML"],
            rename=[],
            source_bam=str(source_bam),
            target_bam=str(target_bam),
            output_bam=str(output_st),
            threads=1,
        )

        # Multi-threaded
        transfer_tags(
            tags=["ML"],
            rename=[],
            source_bam=str(source_bam),
            target_bam=str(target_bam),
            output_bam=str(output_mt),
            threads=4,
        )

        # Compare outputs
        with pysam.AlignmentFile(str(output_st), "rb") as st_bam:
            st_reads = list(st_bam)
        with pysam.AlignmentFile(str(output_mt), "rb") as mt_bam:
            mt_reads = list(mt_bam)

        assert len(st_reads) == len(mt_reads) == 10
        for st_read, mt_read in zip(st_reads, mt_reads):
            assert st_read.query_name == mt_read.query_name
            # Single-element arrays are unwrapped to scalar by transfer_tags
            assert st_read.get_tag("ML") == mt_read.get_tag("ML")


class TestParseSetTags:
    """Tests for parse_set_tags."""

    def test_string_tag(self):
        assert parse_set_tags(["BC:Z:ldx04"]) == [("BC", "ldx04", "Z")]

    def test_integer_tag_is_typed(self):
        """`i` values must become ints, or pysam writes them as strings."""
        assert parse_set_tags(["nn:i:7"]) == [("nn", 7, "i")]

    def test_value_may_contain_colons(self):
        """Split on the first two colons only."""
        assert parse_set_tags(["co:Z:a:b:c"]) == [("co", "a:b:c", "Z")]

    def test_multiple(self):
        assert parse_set_tags(["BC:Z:ldx04", "xx:Z:y"]) == [
            ("BC", "ldx04", "Z"),
            ("xx", "y", "Z"),
        ]

    def test_malformed_raises(self):
        with pytest.raises(ValueError):
            parse_set_tags(["BC:ldx04"])


def _write_bam(path, header, tags_per_read, names=("read1",)):
    """Write a minimal single-ref BAM with one read per name."""
    with pysam.AlignmentFile(str(path), "wb", header=header) as outf:
        for name in names:
            read = pysam.AlignedSegment()
            read.query_name = name
            read.query_sequence = "A" * 100
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 0
            read.cigartuples = [(0, 100)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            for tag, val in tags_per_read.items():
                read.set_tag(tag, val)
            outf.write(read)


class TestSetTag:
    """Constant tags stamped onto every output read."""

    HEADER = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "ref", "LN": 100}]}

    def test_set_tag_applied(self, temp_dir):
        src, tgt, out = (temp_dir / f"{n}.bam" for n in ("s", "t", "o"))
        _write_bam(src, self.HEADER, {"mv": "x"})
        _write_bam(tgt, self.HEADER, {})

        transfer_tags(
            tags=[],
            rename=[],
            source_bam=str(src),
            target_bam=str(tgt),
            output_bam=str(out),
            all_tags=True,
            set_tags=["BC:Z:ldx04"],
        )

        with pysam.AlignmentFile(str(out), "rb") as bam:
            read = next(bam)
            assert read.get_tag("BC") == "ldx04"

    def test_applied_to_reads_with_no_source_match(self, temp_dir):
        """--all-tags emits unmatched target reads; they must be tagged too.

        Otherwise a sample's BAM carries the barcode on only some of its reads,
        which is worse than not carrying it at all.
        """
        src, tgt, out = (temp_dir / f"{n}.bam" for n in ("s", "t", "o"))
        _write_bam(src, self.HEADER, {"mv": "x"}, names=("read1",))
        _write_bam(tgt, self.HEADER, {}, names=("read1", "orphan"))

        transfer_tags(
            tags=[],
            rename=[],
            source_bam=str(src),
            target_bam=str(tgt),
            output_bam=str(out),
            all_tags=True,
            set_tags=["BC:Z:ldx04"],
        )

        with pysam.AlignmentFile(str(out), "rb") as bam:
            tagged = {r.query_name: r.get_tag("BC") for r in bam}
        assert tagged == {"read1": "ldx04", "orphan": "ldx04"}

    def test_absent_when_not_requested(self, temp_dir):
        """Unbarcoded samples must get no BC at all, not an empty one."""
        src, tgt, out = (temp_dir / f"{n}.bam" for n in ("s", "t", "o"))
        _write_bam(src, self.HEADER, {"mv": "x"})
        _write_bam(tgt, self.HEADER, {})

        transfer_tags(
            tags=[],
            rename=[],
            source_bam=str(src),
            target_bam=str(tgt),
            output_bam=str(out),
            all_tags=True,
        )

        with pysam.AlignmentFile(str(out), "rb") as bam:
            assert not next(bam).has_tag("BC")


class TestReadGroupRepair:
    """The output header must declare the read groups its reads reference."""

    SRC_HEADER = {
        "HD": {"VN": "1.0"},
        "SQ": [{"SN": "ref", "LN": 100}],
        "RG": [
            {
                "ID": "runid_model",
                "PL": "ONT",
                "PU": "FLOW1",
                "SM": "ont_run_name",
                "LB": "ont_lib",
            }
        ],
    }
    TGT_HEADER = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "ref", "LN": 100}]}

    def _run(self, temp_dir, **kwargs):
        src, tgt, out = (temp_dir / f"{n}.bam" for n in ("s", "t", "o"))
        _write_bam(src, self.SRC_HEADER, {"RG": "runid_model"})
        _write_bam(tgt, self.TGT_HEADER, {})
        transfer_tags(
            tags=[],
            rename=[],
            source_bam=str(src),
            target_bam=str(tgt),
            output_bam=str(out),
            all_tags=True,
            **kwargs,
        )
        return out

    def test_rg_header_is_restored(self, temp_dir):
        """Regression: bwa drops @RG, --all-tags copies RG:Z back -> invalid SAM."""
        out = self._run(temp_dir)
        with pysam.AlignmentFile(str(out), "rb") as bam:
            declared = {rg["ID"] for rg in bam.header.to_dict().get("RG", [])}
            used = {r.get_tag("RG") for r in bam}
        assert used, "test setup: reads should carry RG"
        assert used <= declared, f"reads reference undeclared @RG: {used - declared}"

    def test_id_preserved_but_identity_overwritten(self, temp_dir):
        """ID must survive (RG:Z points at it); SM/LB/BC become ours."""
        out = self._run(
            temp_dir, rg_sample="ldx04", rg_library="run123", rg_barcode="ldx04"
        )
        with pysam.AlignmentFile(str(out), "rb") as bam:
            rg = bam.header.to_dict()["RG"][0]
        assert rg["ID"] == "runid_model"
        assert rg["PU"] == "FLOW1" and rg["PL"] == "ONT"  # provenance kept
        assert rg["SM"] == "ldx04"
        assert rg["LB"] == "run123"
        assert rg["BC"] == "ldx04"

    def test_comment_recorded(self, temp_dir):
        out = self._run(temp_dir, comments=["aa-tRNA-seq:upstream_barcode=nbc04"])
        with pysam.AlignmentFile(str(out), "rb") as bam:
            assert "aa-tRNA-seq:upstream_barcode=nbc04" in bam.header.to_dict()["CO"]

    def test_no_source_rg_is_not_an_error(self, temp_dir):
        """Non-demux/older BAMs may have no @RG at all."""
        src, tgt, out = (temp_dir / f"{n}.bam" for n in ("s", "t", "o"))
        _write_bam(src, self.TGT_HEADER, {"mv": "x"})
        _write_bam(tgt, self.TGT_HEADER, {})
        transfer_tags(
            tags=[],
            rename=[],
            source_bam=str(src),
            target_bam=str(tgt),
            output_bam=str(out),
            all_tags=True,
            rg_sample="s1",
        )
        with pysam.AlignmentFile(str(out), "rb") as bam:
            assert "RG" not in bam.header.to_dict()
