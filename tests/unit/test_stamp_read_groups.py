"""
Unit tests for stamp_read_groups.py

The @RG/@CO header lines bwa_align hands to `bwa mem -H`, so that the per-read
RG:Z: that `bwa mem -C` copies through from dorado's uBAM resolves to a
declared read group in the aligned BAM.
"""

import pysam

from stamp_read_groups import (
    format_header_lines,
    read_groups_from_bam,
    stamp_read_groups,
)

DORADO_RG = {
    "ID": "8b08a11a_rna004_sup@v6.0.0",
    "PU": "PAU09453",
    "PM": "ont-p2-jh",
    "DT": "2024-07-24T20:41:37.673000+00:00",
    "PL": "ONT",
    "DS": "runid=8b08a11a basecall_model=rna004_sup@v6.0.0",
    "LB": "JMW_510_28C",
}


class TestStampReadGroups:
    def test_overwrites_only_identity_fields(self):
        (rg,) = stamp_read_groups(
            [DORADO_RG], sample="s1", library="run1", barcode="ldx04"
        )
        assert rg["SM"] == "s1"
        assert rg["LB"] == "run1"
        assert rg["BC"] == "ldx04"
        # dorado's ID is what every read's RG:Z: points at; it must not move
        assert rg["ID"] == DORADO_RG["ID"]
        for key in ("PU", "PM", "DT", "PL", "DS"):
            assert rg[key] == DORADO_RG[key]

    def test_absent_barcode_writes_no_bc(self):
        """Absence means 'no demultiplexing', not 'unknown barcode'."""
        (rg,) = stamp_read_groups([DORADO_RG], sample="s1")
        assert "BC" not in rg
        # and dorado's own LB survives when no run id is given
        assert rg["LB"] == DORADO_RG["LB"]

    def test_input_is_not_mutated(self):
        original = dict(DORADO_RG)
        stamp_read_groups([DORADO_RG], sample="s1", barcode="bc")
        assert DORADO_RG == original

    def test_every_read_group_is_stamped(self):
        rgs = [dict(DORADO_RG), {**DORADO_RG, "ID": "other_run"}]
        stamped = stamp_read_groups(rgs, sample="s1")
        assert [rg["SM"] for rg in stamped] == ["s1", "s1"]
        assert [rg["ID"] for rg in stamped] == [DORADO_RG["ID"], "other_run"]


class TestFormatHeaderLines:
    def test_rg_line_is_tab_separated_sam(self):
        (line,) = format_header_lines([{"ID": "x", "SM": "s1"}])
        assert line == "@RG\tID:x\tSM:s1"

    def test_comments_follow_read_groups(self):
        lines = format_header_lines(
            [{"ID": "x"}], comments=["aa-tRNA-seq:upstream_barcode=bc03"]
        )
        assert lines == ["@RG\tID:x", "@CO\taa-tRNA-seq:upstream_barcode=bc03"]

    def test_no_read_groups_gives_no_lines(self):
        assert format_header_lines([]) == []

    def test_lines_parse_back_as_a_valid_header(self):
        """What bwa -H inserts must be something htslib accepts as a header."""
        lines = format_header_lines(
            stamp_read_groups([DORADO_RG], sample="s1", barcode="ldx04"),
            comments=["note"],
        )
        header = pysam.AlignmentHeader.from_text(
            "@HD\tVN:1.6\tSO:unsorted\n" + "\n".join(lines) + "\n"
        ).to_dict()
        assert header["RG"][0]["SM"] == "s1"
        assert header["RG"][0]["BC"] == "ldx04"
        assert header["CO"] == ["note"]


class TestReadGroupsFromBam:
    def test_reads_rg_from_an_unaligned_bam(self, temp_dir):
        ubam = temp_dir / "reads.bam"
        header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "RG": [DORADO_RG]}
        with pysam.AlignmentFile(str(ubam), "wb", header=header) as fh:
            read = pysam.AlignedSegment(fh.header)
            read.query_name = "r1"
            read.query_sequence = "ACGT"
            read.query_qualities = pysam.qualitystring_to_array("IIII")
            read.flag = 4
            read.set_tag("RG", DORADO_RG["ID"])
            fh.write(read)
        (rg,) = read_groups_from_bam(str(ubam))
        assert rg["ID"] == DORADO_RG["ID"]
        assert rg["DS"] == DORADO_RG["DS"]

    def test_bam_without_rg_gives_empty_list(self, temp_dir):
        ubam = temp_dir / "reads.bam"
        with pysam.AlignmentFile(
            str(ubam), "wb", header={"HD": {"VN": "1.6", "SO": "unsorted"}}
        ):
            pass
        assert read_groups_from_bam(str(ubam)) == []
