"""
Unit tests for split_ubam_by_sample.py

One pass over a run-level uBAM must give each sample exactly the reads
`samtools view -N` would, in the same order.
"""

import pysam
import pytest
from split_ubam_by_sample import load_routes, split


def _write_ubam(path, ids):
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}})
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for rid in ids:
            rec = pysam.AlignedSegment(header)
            rec.query_name = rid
            rec.query_sequence = "ACGT"
            rec.flag = 4
            out.write(rec)


def _ids(path):
    with pysam.AlignmentFile(str(path), "rb", check_sq=False) as bam:
        return [r.query_name for r in bam.fetch(until_eof=True)]


@pytest.fixture
def run(tmp_path):
    bam = tmp_path / "run.bam"
    _write_ubam(bam, ["r1", "r2", "r3", "r4", "r5"])
    (tmp_path / "a.txt").write_text("r1\nr4\nr_absent\n")
    (tmp_path / "b.txt").write_text("r3\nr2\n\n")
    return tmp_path, bam


def test_each_sample_gets_its_reads_in_run_order(run):
    tmp, bam = run
    routes, outs = load_routes(
        [f"a={tmp}/a.txt={tmp}/o/a.bam", f"b={tmp}/b.txt={tmp}/o/b.bam"]
    )
    counts = split(bam, routes, outs, threads=1)
    assert counts == {"a": 2, "b": 2}
    assert _ids(tmp / "o/a.bam") == ["r1", "r4"]
    assert _ids(tmp / "o/b.bam") == ["r2", "r3"]


def test_read_listed_for_two_samples_goes_to_both(run):
    tmp, bam = run
    (tmp / "c.txt").write_text("r1\n")
    routes, outs = load_routes(
        [f"a={tmp}/a.txt={tmp}/a.bam", f"c={tmp}/c.txt={tmp}/c.bam"]
    )
    split(bam, routes, outs, threads=1)
    assert "r1" in _ids(tmp / "a.bam") and _ids(tmp / "c.bam") == ["r1"]


def test_duplicate_sample_rejected(run):
    tmp, _ = run
    with pytest.raises(SystemExit):
        load_routes([f"a={tmp}/a.txt=x.bam", f"a={tmp}/b.txt=y.bam"])
