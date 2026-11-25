"""Unit tests for generate_squiggy_session.py."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

from generate_squiggy_session import (
    generate_squiggy_session,
    read_pod5_map,
    read_samples_file,
)

SCRIPT = (
    Path(__file__).parents[2] / "workflow" / "scripts" / "generate_squiggy_session.py"
)


def run_script(tmp_path, *extra, name="session.json"):
    out = tmp_path / name
    cmd = [
        sys.executable,
        str(SCRIPT),
        "--output-dir",
        str(tmp_path / "out"),
        "--fasta",
        str(tmp_path / "ref.fa"),
        "--output",
        str(out),
        "--no-checksums",
        *extra,
    ]
    return subprocess.run(cmd, capture_output=True, text=True), out


def session_of(out):
    data = json.loads(out.read_text())
    data.pop("timestamp")
    return data


class TestFileInputs:
    def test_files_match_flags(self, tmp_path):
        """Same samples and POD5s give the same JSON by flag or by file."""
        samples = ["a", "b", "c"]
        pod5 = {"a": ["/r/1.pod5", "/r/2.pod5"], "b": ["/r/1.pod5"]}

        flags = []
        for s, paths in pod5.items():
            for p in paths:
                flags += ["--pod5", f"{s}={p}"]
        r1, out1 = run_script(
            tmp_path, "--samples", *samples, *flags, name="flags.json"
        )
        assert r1.returncode == 0, r1.stderr

        (tmp_path / "samples.txt").write_text("".join(f"{s}\n" for s in samples))
        (tmp_path / "map.tsv").write_text(
            "".join(f"{s}\t{p}\n" for s, ps in pod5.items() for p in ps)
        )
        r2, out2 = run_script(
            tmp_path,
            "--samples-file",
            str(tmp_path / "samples.txt"),
            "--pod5-map",
            str(tmp_path / "map.tsv"),
            name="files.json",
        )
        assert r2.returncode == 0, r2.stderr

        a, b = session_of(out1), session_of(out2)
        assert a == b
        assert list(b["samples"]) == samples
        assert b["samples"]["a"]["pod5Paths"] == ["/r/1.pod5", "/r/2.pod5"]

    def test_no_samples_is_an_error(self, tmp_path):
        r, _ = run_script(tmp_path)
        assert r.returncode != 0
        assert "no samples" in r.stderr

    def test_bad_map_row_is_an_error(self, tmp_path):
        (tmp_path / "samples.txt").write_text("a\n")
        (tmp_path / "map.tsv").write_text("a /no/tab/here.pod5\n")
        r, _ = run_script(
            tmp_path,
            "--samples-file",
            str(tmp_path / "samples.txt"),
            "--pod5-map",
            str(tmp_path / "map.tsv"),
        )
        assert r.returncode != 0
        assert "expected SAMPLE<TAB>PATH" in r.stderr


class TestReaders:
    def test_samples_file_skips_blanks_and_comments(self, tmp_path):
        f = tmp_path / "s.txt"
        f.write_text("# header\na\n\n  b  \n")
        assert read_samples_file(f) == ["a", "b"]

    def test_pod5_map_groups_in_file_order(self, tmp_path):
        f = tmp_path / "m.tsv"
        f.write_text("a\t/2\nb\t/1\na\t/1\n")
        assert read_pod5_map(f) == {"a": ["/2", "/1"], "b": ["/1"]}


class TestLargeRun:
    """The failure was `Argument list too long` at 884 libraries x 108 POD5s."""

    def test_884_libraries_108_pod5(self, tmp_path):
        samples = [f"lib{i:04d}" for i in range(884)]
        pod5s = [f"/data/run/pod5_pass/file_{j:03d}.pod5" for j in range(108)]
        (tmp_path / "samples.txt").write_text("".join(f"{s}\n" for s in samples))
        (tmp_path / "map.tsv").write_text(
            "".join(f"{s}\t{p}\n" for s in samples for p in pod5s)
        )
        r, out = run_script(
            tmp_path,
            "--samples-file",
            str(tmp_path / "samples.txt"),
            "--pod5-map",
            str(tmp_path / "map.tsv"),
        )
        assert r.returncode == 0, r.stderr
        data = json.loads(out.read_text())
        assert len(data["samples"]) == 884
        assert len(data["samples"]["lib0883"]["pod5Paths"]) == 108

    def test_flag_form_would_overflow(self):
        """Documents why the files exist.

        The rule ran `python ... --pod5 ...` through `bash -c`, so the whole
        command is ONE argument, capped at MAX_ARG_STRLEN (32 pages = 128 KiB on
        Linux) -- far below ARG_MAX, which is why a 477-library run passed.
        """
        samples = [f"lib{i:04d}" for i in range(884)]
        pod5s = [f"/data/run/pod5_pass/file_{j:03d}.pod5" for j in range(108)]
        old = " ".join(f"--pod5 {s}={p}" for s in samples for p in pod5s)
        assert len(old) > 32 * 4096
