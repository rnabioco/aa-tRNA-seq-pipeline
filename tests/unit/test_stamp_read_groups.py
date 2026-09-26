"""
Unit tests for stamp_read_groups.py

The uBAM rewrite escpod_align runs before `escpod align`: SM/LB/BC stamped on
dorado's @RG (so the per-read RG:Z: escpod copies through still resolves to a
declared read group) and, on a barcoded sample, a constant BC:Z: on every
record -- both in one pass.
"""

from array import array

import pysam
import pytest

from stamp_read_groups import (
    format_header_lines,
    read_groups_from_bam,
    stamp_bam,
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
        """The stamped header lines must be something htslib accepts."""
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


def _write_ubam(path, read_groups=(DORADO_RG,), n_reads=3, extra_tags=()):
    """A small dorado-shaped uBAM: unmapped records with RG, a move table and
    MM/ML, the tags escpod align must see unchanged."""
    header = {"HD": {"VN": "1.6", "SO": "unknown"}, "PG": [{"ID": "basecaller"}]}
    if read_groups:
        header["RG"] = [dict(rg) for rg in read_groups]
    with pysam.AlignmentFile(str(path), "wb", header=header) as fh:
        for i in range(n_reads):
            read = pysam.AlignedSegment(fh.header)
            read.query_name = f"r{i}"
            read.query_sequence = "ACGTACGT"
            read.query_qualities = pysam.qualitystring_to_array("IIIIIIII")
            read.flag = 4
            if read_groups:
                read.set_tag("RG", read_groups[0]["ID"])
            read.set_tag("mv", array("b", [6, 1, 0, 1, 0, 0, 1]))
            read.set_tag("ns", 42 + i)
            read.set_tag("MM", "A+a?,0;")
            read.set_tag("ML", array("B", [200]))
            for tag, value in extra_tags:
                read.set_tag(tag, value)
            fh.write(read)


def _read_all(path):
    with pysam.AlignmentFile(str(path), "rb", check_sq=False) as fh:
        return fh.header.to_dict(), [r for r in fh.fetch(until_eof=True)]


class TestStampBam:
    """The combined reheader + tag-stamp pass escpod_align runs."""

    @pytest.mark.parametrize("compress", [False, True])
    def test_barcoded_sample_gets_rg_fields_and_bc_on_every_record(
        self, temp_dir, compress
    ):
        ubam, out = temp_dir / "in.bam", temp_dir / "out.bam"
        _write_ubam(ubam)
        n = stamp_bam(
            str(ubam),
            str(out),
            sample="s1",
            library="run1",
            barcode="ldx04-fdx01",
            compress=compress,
        )
        header, reads = _read_all(out)
        assert n == 3 and len(reads) == 3
        (rg,) = header["RG"]
        assert rg["SM"] == "s1"
        assert rg["LB"] == "run1"
        assert rg["BC"] == "ldx04-fdx01"
        assert [r.get_tag("BC") for r in reads] == ["ldx04-fdx01"] * 3

    def test_rg_identity_fields_untouched(self, temp_dir):
        ubam, out = temp_dir / "in.bam", temp_dir / "out.bam"
        _write_ubam(ubam)
        stamp_bam(str(ubam), str(out), sample="s1", library="run1", barcode="bc")
        header, reads = _read_all(out)
        (rg,) = header["RG"]
        for key in ("ID", "PU", "PM", "DT", "PL", "DS"):
            assert rg[key] == DORADO_RG[key]
        # every read still points at the (unchanged) declared read group
        assert {r.get_tag("RG") for r in reads} == {DORADO_RG["ID"]}

    def test_unbarcoded_sample_gets_no_bc_anywhere(self, temp_dir):
        """Absence means 'no demultiplexing', on the @RG and on every record."""
        ubam, out = temp_dir / "in.bam", temp_dir / "out.bam"
        _write_ubam(ubam)
        stamp_bam(str(ubam), str(out), sample="s1")
        header, reads = _read_all(out)
        (rg,) = header["RG"]
        assert "BC" not in rg
        assert rg["SM"] == "s1"
        assert rg["LB"] == DORADO_RG["LB"]
        assert not any(r.has_tag("BC") for r in reads)

    def test_dorado_tags_and_records_pass_through(self, temp_dir):
        ubam, out = temp_dir / "in.bam", temp_dir / "out.bam"
        _write_ubam(ubam)
        stamp_bam(str(ubam), str(out), sample="s1", barcode="bc")
        _, before = _read_all(ubam)
        _, after = _read_all(out)
        assert [r.query_name for r in after] == [r.query_name for r in before]
        for a, b in zip(after, before):
            assert a.query_sequence == b.query_sequence
            assert a.flag == b.flag
            for tag in ("mv", "ns", "MM", "ML", "RG"):
                assert a.get_tag(tag) == b.get_tag(tag)

    def test_existing_bc_tag_is_replaced_not_duplicated(self, temp_dir):
        ubam, out = temp_dir / "in.bam", temp_dir / "out.bam"
        _write_ubam(ubam, extra_tags=[("BC", "upstream")])
        stamp_bam(str(ubam), str(out), sample="s1", barcode="ldx04")
        _, reads = _read_all(out)
        for read in reads:
            assert read.get_tag("BC") == "ldx04"
            assert [t for t, _ in read.get_tags()].count("BC") == 1

    def test_comments_and_other_header_lines_kept(self, temp_dir):
        ubam, out = temp_dir / "in.bam", temp_dir / "out.bam"
        _write_ubam(ubam)
        stamp_bam(
            str(ubam),
            str(out),
            sample="s1",
            barcode="ldx04",
            comments=["aa-tRNA-seq:upstream_barcode=bc04"],
        )
        header, _ = _read_all(out)
        assert header["CO"] == ["aa-tRNA-seq:upstream_barcode=bc04"]
        assert header["PG"][0]["ID"] == "basecaller"
        assert header["HD"]["VN"] == "1.6"

    def test_every_read_group_stamped(self, temp_dir):
        ubam, out = temp_dir / "in.bam", temp_dir / "out.bam"
        _write_ubam(ubam, read_groups=(DORADO_RG, {**DORADO_RG, "ID": "other"}))
        stamp_bam(str(ubam), str(out), sample="s1", barcode="bc")
        header, _ = _read_all(out)
        assert [rg["SM"] for rg in header["RG"]] == ["s1", "s1"]
        assert [rg["ID"] for rg in header["RG"]] == [DORADO_RG["ID"], "other"]

    def test_ubam_without_rg_warns_and_still_writes(self, temp_dir, capsys):
        ubam, out = temp_dir / "in.bam", temp_dir / "out.bam"
        _write_ubam(ubam, read_groups=())
        assert stamp_bam(str(ubam), str(out), sample="s1") == 3
        header, _ = _read_all(out)
        assert "RG" not in header
        assert "declares no @RG" in capsys.readouterr().err
