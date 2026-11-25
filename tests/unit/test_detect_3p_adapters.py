"""Unit tests for detect_3p_adapters.py and the unassigned columns of
edx_concordance.py (#223).

The panel mimics the one in the issue: the plain `edx00` (shared 20-nt stem +
A13), three 45-nt `edx` adapters and twelve 60-nt `ndx` adapters. Barcodes are
seeded random; only the geometry (stem, barcode, poly-A) matters here.
"""

import gzip
import random
import subprocess
import sys
from pathlib import Path

import pytest

from add_adapter_tags import create_scoring_matrix
from detect_3p_adapters import UNASSIGNED, classify_read, self_score

STEM = "GGCTTCTTCTTGCTCTTAGG"
EDX00 = STEM + "A" * 13
GATE = dict(min_score_frac=0.6, min_margin_frac=0.05, residual=("edx00",))
SCRIPTS = Path(__file__).parents[2] / "workflow" / "scripts"


def _barcode(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def build_panel(seed=7):
    rng = random.Random(seed)
    panel = [("edx00", EDX00)]
    panel += [(f"edx0{i}", STEM + _barcode(rng, 12) + "A" * 13) for i in (2, 3, 4)]
    panel += [
        (f"ndx{i:02d}", STEM + _barcode(rng, 25) + "A" * 15) for i in range(1, 13)
    ]
    return panel


PANEL = build_panel()
SEQS = dict(PANEL)
assert [len(SEQS[n]) for n in ("edx00", "edx02", "ndx01")] == [33, 45, 60]
MATRIX = create_scoring_matrix()


def noisy(seq, rng, rate=0.10):
    """Substitutions, insertions and deletions at `rate` per base in total."""
    out = []
    for base in seq:
        r = rng.random()
        if r < rate * 0.6:
            out.append(rng.choice([b for b in "ACGT" if b != base]))
        elif r < rate * 0.8:
            out.append(base)
            out.append(rng.choice("ACGT"))
        elif r < rate:
            continue
        else:
            out.append(base)
    return "".join(out)


def make_read(adapter_seq, rng, truncate=0, rate=0.10):
    """A read body followed by the adapter, with `truncate` 3' nt of it lost."""
    body = _barcode(rng, 60)
    kept = adapter_seq[: len(adapter_seq) - truncate]
    return body + noisy(kept, rng, rate)


def call(seq, **kw):
    return classify_read(seq, PANEL, MATRIX, **kw)


class TestFullLength:
    def test_detect_3p_full_length_adapters_called_correctly(self):
        rng = random.Random(1)
        wrong = []
        for name, seq in PANEL:
            for _ in range(25):
                got = call(make_read(seq, rng), **GATE)["adapter_3p"]
                if got != name:
                    wrong.append((name, got))
        # 10 % error on a 33-nt adapter is rough; allow the odd stray, not a
        # systematic miscall. The gate must not be what removes true calls.
        assert len(wrong) <= 0.03 * 25 * len(PANEL), wrong


class TestTruncation:
    def test_detect_3p_truncated_barcode_is_unassigned_not_forced(self):
        rng = random.Random(2)
        for name, seq in PANEL:
            if not name.startswith("ndx"):
                continue
            for _ in range(10):
                c = call(make_read(seq, rng, truncate=25), **GATE)
                assert c["adapter_3p"] == UNASSIGNED, (name, c)
                assert c["reason"] in {"low_score", "low_margin", "stem_only"}

    def test_detect_3p_truncated_edx_does_not_fall_into_edx00(self):
        rng = random.Random(3)
        for name in ("edx02", "edx03", "edx04"):
            for _ in range(15):
                c = call(make_read(SEQS[name], rng, truncate=15), **GATE)
                assert c["adapter_3p"] != "edx00", (name, c)
                assert c["adapter_3p"] in (UNASSIGNED, name), (name, c)

    def test_stem_only_is_edx00_when_polya_tail_is_present(self):
        rng = random.Random(4)
        clean = [
            call(make_read(EDX00, rng, rate=0.0), **GATE)["adapter_3p"]
            for _ in range(30)
        ]
        assert set(clean) == {"edx00"}
        noisy5 = [
            call(make_read(EDX00, rng, rate=0.05), **GATE)["adapter_3p"]
            for _ in range(30)
        ]
        # At 5 % error the tail alignment sometimes drops under the 0.8 bar; a
        # read lost to `unassigned` is the cost of refusing to guess, a read
        # given to another adapter would be the bug.
        assert noisy5.count("edx00") >= 24
        assert set(noisy5) <= {"edx00", UNASSIGNED}


class TestColumns:
    def test_hit_position_and_aligned_bases_are_reported(self):
        # One adapter, so the call is about the columns, not about which of the
        # panel a truncated read resembles.
        rng = random.Random(5)
        seq = SEQS["ndx05"]
        only = [("ndx05", seq)]
        read = make_read(seq, rng, truncate=25, rate=0.0)
        c = classify_read(read, only, MATRIX)
        # Error-free read cut 25 nt short: the hit reaches the end of the read
        # and only the 35 adapter bases that are there aligned.
        assert c["read_len"] == len(read)
        assert len(read) - c["hit_end"] <= 2
        assert 33 <= c["n_aligned"] <= 35
        assert float(c["score_frac"]) < 0.6
        full = classify_read(make_read(seq, rng, rate=0.0), only, MATRIX)
        assert full["n_aligned"] == 60 and full["adapter_3p"] == "ndx05"
        assert float(full["score_frac"]) == pytest.approx(1.0)

    def test_truncated_ndx_resembles_a_shorter_adapter(self):
        """Why the gate is on the score and not on the nearest adapter's name.

        An ndx read cut 25 nt short is stem + a bit of barcode, which a 45-nt
        adapter fits as well as its own. The nearest adapter is then arbitrary;
        the gate is what stops it being believed.
        """
        rng = random.Random(5)
        read = make_read(SEQS["ndx05"], rng, truncate=25, rate=0.0)
        c = call(read, **GATE)
        assert c["adapter_3p"] == UNASSIGNED and c["reason"] == "stem_only"

    def test_gate_off_keeps_nearest_adapter_call(self):
        rng = random.Random(6)
        read = make_read(SEQS["ndx05"], rng, truncate=25)
        c = call(read)  # no gate, no residual
        assert c["adapter_3p"] == c["adapter_nearest"] != UNASSIGNED

    def test_no_adapter_is_none(self):
        assert call("ACGT" * 5, **GATE)["adapter_3p"] == "none"
        assert call("", **GATE)["adapter_3p"] == "none"

    def test_self_score_scales_with_length(self):
        assert self_score(EDX00) == 66
        assert self_score(SEQS["ndx01"]) == 120


def _write_adapter_tsv(path, rows):
    header = "read_id\tadapter_3p\tscore_best\tscore_second\tmargin\n"
    with gzip.open(path, "wt") as f:
        f.write(header)
        for i, (adapter, margin) in enumerate(rows):
            f.write(f"r{i}\t{adapter}\t50\t40\t{margin}\n")


def test_edx_concordance_reports_unassigned(tmp_path):
    tsv = tmp_path / "pool.tsv.gz"
    _write_adapter_tsv(
        tsv,
        [("edx02", 20)] * 6
        + [("edx03", 2)] * 2
        + [(UNASSIGNED, 1)] * 3
        + [("none", "")],
    )
    out = tmp_path / "conc.tsv.gz"
    subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "edx_concordance.py"),
            "--tsvs",
            str(tsv),
            "--samples",
            "pool",
            "--min-margin",
            "5",
            "--output",
            str(out),
        ],
        check=True,
        capture_output=True,
    )
    with gzip.open(out, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = [
            dict(zip(header, line.rstrip("\n").split("\t"), strict=True)) for line in f
        ]
    assert {"n_unassigned", "pct_unassigned"} <= set(header)
    by = {r["edx_adapter"]: r for r in rows}
    # 12 reads in the pool; 3 unassigned, repeated on every row.
    assert all(r["n_unassigned"] == "3" and r["pct_unassigned"] == "25.0" for r in rows)
    assert by[UNASSIGNED]["n_reads"] == "3"
    # Unassigned reads are never "clean", and do not dilute the clean fractions.
    assert by[UNASSIGNED]["n_clean"] == "0"
    assert by["edx02"]["n_clean"] == "6" and by["edx03"]["n_clean"] == "0"
    assert by["edx02"]["pct_clean"] == "100.0"
