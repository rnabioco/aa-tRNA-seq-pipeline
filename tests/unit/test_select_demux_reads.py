"""
Unit tests for select_demux_reads.py

The per-sample read selection that joins escpod classifications -- one CSV per
axis from separate passes, or one fused CSV with per-axis columns -- into the
read list a sample owns.
"""

import pytest

from select_demux_reads import (
    SelectionError,
    add_split_children,
    main,
    parse_axes,
    parse_gates,
    parse_samples,
    select_reads,
)


def write_csv(path, header, rows):
    with open(path, "w") as fh:
        fh.write(",".join(header) + "\n")
        for row in rows:
            fh.write(",".join(str(v) for v in row) + "\n")


@pytest.fixture
def two_axes(temp_dir):
    """Two single-model CSVs, as two escpod passes write them."""
    ldx = temp_dir / "ldx.csv"
    fdx = temp_dir / "fdx.csv"
    write_csv(
        ldx,
        ["read_id", "barcode", "confidence", "crf_logp", "crf_margin"],
        [
            ("r1", "ldx01", 14, -0.1, 12.0),
            ("r2", "ldx01", 14, -0.1, 12.0),
            ("r3", "ldx01", 14, -0.1, 12.0),
            ("r4", "ldx04", 14, -0.1, 12.0),
            ("r5", "ldx04", 14, -0.1, 12.0),
            ("r6", "ldx10", 14, -0.1, 12.0),
            ("r7", "unclassified", 0, "", ""),
        ],
    )
    write_csv(
        fdx,
        ["read_id", "barcode", "confidence", "crf_logp", "crf_margin"],
        [
            ("r1", "fdx01", 16, -0.2, 18.0),
            ("r2", "fdx01", 16, -0.2, 4.0),
            ("r3", "fdx02", 16, -0.2, 18.0),  # ldx01 read whose 5' index disagrees
            ("r4", "fdx02", 16, -0.2, 2.0),  # below the fdx gate of 3.5
            ("r5", "unclassified", 0, "", ""),
            ("r6", "fdx04", 16, -0.2, 18.0),
            ("r7", "fdx01", 16, -0.2, 18.0),
        ],
    )
    return [("ldx", str(ldx)), ("fdx", str(fdx))]


@pytest.fixture
def fused_csv(temp_dir):
    """One fused CSV with per-axis columns, the same reads as `two_axes`."""
    path = temp_dir / "fused.csv"
    header = [
        "read_id", "ldx", "ldx_confidence", "ldx_crf_logp", "ldx_crf_margin",
        "fdx", "fdx_confidence", "fdx_crf_logp", "fdx_crf_margin",
    ]
    write_csv(
        path,
        header,
        [
            ("r1", "ldx01", 14, -0.1, 12.0, "fdx01", 16, -0.2, 18.0),
            ("r2", "ldx01", 14, -0.1, 12.0, "fdx01", 16, -0.2, 4.0),
            ("r3", "ldx01", 14, -0.1, 12.0, "fdx02", 16, -0.2, 18.0),
            ("r4", "ldx04", 14, -0.1, 12.0, "fdx02", 16, -0.2, 2.0),
            ("r5", "ldx04", 14, -0.1, 12.0, "unclassified", 0, "", ""),
            ("r6", "ldx10", 14, -0.1, 12.0, "fdx04", 16, -0.2, 18.0),
            ("r7", "unclassified", 0, "", "", "fdx01", 16, -0.2, 18.0),
        ],
    )
    return str(path)


DUAL = {"a": {"ldx": "ldx01", "fdx": "fdx01"}, "b": {"ldx": "ldx04", "fdx": "fdx02"}}


class TestParsing:
    def test_axes_keep_order_and_primary(self):
        axes = parse_axes(["ldx=a.csv", "fdx=b.csv"])
        assert axes == [("ldx", "a.csv"), ("fdx", "b.csv")]

    def test_axis_needs_name_and_path(self):
        with pytest.raises(SelectionError):
            parse_axes(["a.csv"])

    def test_gate_names_a_declared_axis(self):
        axes = parse_axes(["ldx=a.csv", "fdx=b.csv"])
        assert parse_gates(["fdx=3.5"], axes) == {"fdx": 3.5}
        with pytest.raises(SelectionError):
            parse_gates(["wdx=1"], axes)
        with pytest.raises(SelectionError):
            parse_gates(["fdx=high"], axes)

    def test_sample_must_name_primary_axis(self):
        axes = parse_axes(["ldx=a.csv", "fdx=b.csv"])
        with pytest.raises(SelectionError, match="primary axis"):
            parse_samples(["s:fdx=fdx01"], axes)

    def test_sample_cannot_name_unknown_axis(self):
        axes = parse_axes(["ldx=a.csv"])
        with pytest.raises(SelectionError, match="no --axis"):
            parse_samples(["s:ldx=ldx01,wdx=barcode03"], axes)

    def test_prefix_tuples_are_refused(self):
        """ldx01 alone would swallow every read of ldx01+fdx01."""
        axes = parse_axes(["ldx=a.csv", "fdx=b.csv"])
        with pytest.raises(SelectionError, match="different axis sets"):
            parse_samples(["a:ldx=ldx01", "b:ldx=ldx01,fdx=fdx01"], axes)

    def test_identical_tuples_are_refused(self):
        axes = parse_axes(["ldx=a.csv", "fdx=b.csv"])
        with pytest.raises(SelectionError, match="identical"):
            parse_samples(["a:ldx=ldx01,fdx=fdx01", "b:ldx=ldx01,fdx=fdx01"], axes)

    def test_same_primary_different_secondary_is_fine(self):
        axes = parse_axes(["ldx=a.csv", "fdx=b.csv"])
        samples = parse_samples(["a:ldx=ldx01,fdx=fdx01", "b:ldx=ldx01,fdx=fdx02"], axes)
        assert samples["a"] == {"ldx": "ldx01", "fdx": "fdx01"}
        assert samples["b"] == {"ldx": "ldx01", "fdx": "fdx02"}


class TestSelectReads:
    def test_single_axis_matches_the_old_behaviour(self, two_axes):
        (ldx,) = two_axes[:1]
        assigned = select_reads([ldx], {"a": {"ldx": "ldx01"}, "d": {"ldx": "ldx04"}})
        assert assigned == {"a": {"r1", "r2", "r3"}, "d": {"r4", "r5"}}

    def test_join_requires_every_axis_to_agree(self, two_axes):
        assigned = select_reads(two_axes, DUAL)
        # r3 is ldx01 but fdx02; r5 is ldx04 but unclassified on fdx
        assert assigned == {"a": {"r1", "r2"}, "b": {"r4"}}

    def test_gate_drops_low_margin_calls(self, two_axes):
        assigned = select_reads(two_axes, DUAL, gates={"fdx": 3.5})
        # r4's fdx margin is 2.0: gated out; r2's 4.0 survives
        assert assigned == {"a": {"r1", "r2"}, "b": set()}

    def test_unclaimed_codes_are_routed_nowhere(self, two_axes):
        assigned = select_reads(two_axes, {"a": {"ldx": "ldx01", "fdx": "fdx01"}})
        assert "r6" not in assigned["a"] and "r7" not in assigned["a"]

    def test_secondary_axis_can_split_one_primary_code(self, two_axes):
        assigned = select_reads(
            two_axes,
            {"a1": {"ldx": "ldx01", "fdx": "fdx01"}, "a2": {"ldx": "ldx01", "fdx": "fdx02"}},
        )
        assert assigned == {"a1": {"r1", "r2"}, "a2": {"r3"}}

    def test_fused_csv_gives_the_same_answer(self, two_axes, fused_csv):
        """The same CSV handed for both axes, read by per-axis column."""
        fused = [("ldx", fused_csv), ("fdx", fused_csv)]
        for gates in (None, {"fdx": 3.5}, {"ldx": 1.0, "fdx": 3.5}):
            assert select_reads(fused, DUAL, gates) == select_reads(two_axes, DUAL, gates)

    def test_csv_without_a_usable_column_is_an_error(self, temp_dir):
        bad = temp_dir / "bad.csv"
        bad.write_text("read_id,code\nr1,ldx01\n")
        with pytest.raises(SelectionError, match="barcode"):
            select_reads([("ldx", str(bad))], {"a": {"ldx": "ldx01"}})

    def test_gate_without_margins_is_an_error(self, temp_dir):
        csv_path = temp_dir / "nomargin.csv"
        csv_path.write_text("read_id,barcode,confidence\nr1,ldx01,14\n")
        with pytest.raises(SelectionError, match="ref-scores"):
            select_reads([("ldx", str(csv_path))], {"a": {"ldx": "ldx01"}}, {"ldx": 1.0})


class TestSplitChildren:
    def test_children_inherit_their_parents_sample(self, temp_dir):
        parents = temp_dir / "split_parents.tsv"
        parents.write_text("c1\tr1\nc2\tr4\nc3\tr9\n")
        assigned = {"a": {"r1"}, "b": {"r4"}}
        added = add_split_children(assigned, str(parents))
        assert assigned == {"a": {"r1", "c1"}, "b": {"r4", "c2"}}
        assert dict(added) == {"a": 1, "b": 1}


class TestMain:
    def test_writes_sorted_union_and_summary(self, two_axes, temp_dir):
        out = temp_dir / "reads.txt"
        summary = temp_dir / "summary.tsv"
        main(
            [
                *(f"--axis={a}={p}" for a, p in two_axes),
                "--sample=a:ldx=ldx01,fdx=fdx01",
                "--sample=b:ldx=ldx04,fdx=fdx02",
                f"--output={out}",
                f"--summary={summary}",
            ]
        )
        assert out.read_text().split() == ["r1", "r2", "r4"]
        assert summary.read_text().splitlines() == [
            "sample\tassigned\tsplit_children",
            "a\t2\t0",
            "b\t1\t0",
        ]

    def test_empty_sample_fails(self, two_axes, temp_dir):
        out = temp_dir / "reads.txt"
        with pytest.raises(SystemExit) as exc:
            main(
                [
                    *(f"--axis={a}={p}" for a, p in two_axes),
                    "--sample=z:ldx=ldx16,fdx=fdx01",
                    f"--output={out}",
                ]
            )
        assert "no reads were assigned" in str(exc.value)
