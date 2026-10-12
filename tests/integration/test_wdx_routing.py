"""
The WarpDemuX path joins the shared raw-run tail (#207).

parse_warpdemux writes `classifications.csv` with `read_id,barcode`, and
select_demux_reads.py selects a sample's reads from it on a `wdx` axis, with no
bundle to resolve codes against and no lattice margin to gate on. These tests
run anywhere; the DAG itself is covered by `pixi run dry-run-demux`.
"""

import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).parent.parent.parent
SCRIPT = REPO_ROOT / "workflow" / "scripts" / "select_demux_reads.py"

CALLS = [
    ("r1", "barcode03"),
    ("r2", "barcode04"),
    ("r3", "barcode04"),
    ("r4", "unclassified"),
    ("r5", "barcode07"),  # a code no sample claims
]


def _select(tmp_path, *samples, parents=None):
    csv_path = tmp_path / "classifications.csv"
    csv_path.write_text("read_id,barcode\n" + "".join(f"{r},{b}\n" for r, b in CALLS))
    out = tmp_path / "ids.txt"
    cmd = [
        sys.executable,
        str(SCRIPT),
        "--axis",
        f"wdx={csv_path}",
        "--output",
        str(out),
    ]
    for s in samples:
        cmd += ["--sample", s]
    if parents:
        cmd += ["--parents", str(parents)]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return set(out.read_text().split())


def test_selects_reads_of_the_configured_barcode(tmp_path):
    assert _select(tmp_path, "s4:wdx=barcode04") == {"r2", "r3"}


def test_all_samples_union_excludes_unclaimed_and_unclassified(tmp_path):
    ids = _select(tmp_path, "s3:wdx=barcode03", "s4:wdx=barcode04")
    assert ids == {"r1", "r2", "r3"}


def test_split_children_inherit_the_parent_barcode(tmp_path):
    parents = tmp_path / "split_parents.tsv"
    parents.write_text("r2a\tr2\nr2b\tr2\nr1a\tr1\n")
    ids = _select(tmp_path, "s4:wdx=barcode04", parents=parents)
    assert ids == {"r2", "r3", "r2a", "r2b"}


# --- sidecar migration: the table may only be deleted once this passes ---

VERIFY = REPO_ROOT / "workflow" / "scripts" / "verify_wdx_sidecar.py"


def _run_verify(tmp_path, sidecar_labels):
    pa = pytest.importorskip("pyarrow")
    pod5_dir = tmp_path / "run" / "pod5"
    pod5_dir.mkdir(parents=True)
    table = pa.table(
        {
            "read_id": [f"s{i}" for i in range(len(sidecar_labels))],
            "barcode": sidecar_labels,
        }
    )
    with (
        pa.OSFile(str(tmp_path / "run" / "pod5.p5s"), "wb") as sink,
        pa.ipc.new_file(sink, table.schema) as writer,
    ):
        writer.write_table(table)
    csv_path = tmp_path / "classifications.csv"
    csv_path.write_text("read_id,barcode\n" + "".join(f"{r},{b}\n" for r, b in CALLS))
    return subprocess.run(
        [sys.executable, str(VERIFY), "--table", str(csv_path), str(pod5_dir)],
        capture_output=True,
        text=True,
    )


def test_verify_accepts_a_sidecar_that_carries_the_table(tmp_path):
    # the sidecar also labels reads the table never named: extra `unclassified`
    labels = [b for _, b in CALLS] + ["unclassified"] * 3
    assert _run_verify(tmp_path, labels).returncode == 0


def test_verify_refuses_a_sidecar_missing_a_read(tmp_path):
    labels = [b for _, b in CALLS][:-1]  # barcode07 never landed
    result = _run_verify(tmp_path, labels)
    assert result.returncode == 1
    assert "barcode07" in result.stderr


def test_verify_refuses_a_missing_sidecar(tmp_path):
    pod5_dir = tmp_path / "run" / "pod5"
    pod5_dir.mkdir(parents=True)
    csv_path = tmp_path / "classifications.csv"
    csv_path.write_text("read_id,barcode\nr1,barcode03\n")
    result = subprocess.run(
        [sys.executable, str(VERIFY), "--table", str(csv_path), str(pod5_dir)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 1
    assert "no sidecar" in result.stderr
