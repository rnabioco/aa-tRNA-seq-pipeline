"""
Unit tests for workflow/scripts/dorado_basecall_resume.sh

The wrapper's whole job is to survive a job being killed: keep the partial
basecall somewhere Snakemake will not reap, hand it back to dorado as
``--resume-from``, and refuse it when continuing would be wrong. None of that
needs a GPU, so it is tested against a stub `dorado` that mimics the parts of
the real one the wrapper depends on -- an ``@PG ID:basecaller`` header carrying
the command line, and copy-through of every record in the resume file.

The stub's behaviour is not a guess -- see the header of the script under test
for the dorado sources it was read off, and PR #150 for the run against the real
binary that confirmed them.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
WRAPPER = REPO_ROOT / "workflow" / "scripts" / "dorado_basecall_resume.sh"

# The wrapper reads BAM headers to decide whether a checkpoint is resumable, so
# samtools is a hard requirement rather than a convenience. `pixi run -e test`
# provides it; a bare `pytest` off an unactivated shell may not.
pytestmark = pytest.mark.skipif(
    shutil.which("samtools") is None, reason="samtools not on PATH"
)

# The stub dorado. Writes a uBAM to stdout the way the real one does: a @PG
# record with ID:basecaller and the full command line in CL, then the reads
# named in the resume file (copied through verbatim), then the reads it was
# asked to basecall that the resume file did not already cover.
#
# FAKE_DORADO_ARGV      - file to record argv into, so tests can assert on it
# FAKE_DORADO_NREADS    - how many reads the "run" contains (default 4)
# FAKE_DORADO_STOP_AFTER- write this many records, then exit 1 (a killed job)
STUB = '''#!{python}
import os, random, sys
import pysam

argv = sys.argv[1:]
if argv and argv[0] == "basecaller":
    argv = argv[1:]

with open(os.environ["FAKE_DORADO_ARGV"], "w") as fh:
    fh.write("\\n".join(argv))

resume = None
for i, a in enumerate(argv):
    if a == "--resume-from":
        resume = argv[i + 1]

carried = []
if resume:
    with pysam.AlignmentFile(resume, "rb", check_sq=False) as fin:
        carried = [r.query_name for r in fin]

nreads = int(os.environ.get("FAKE_DORADO_NREADS", "4"))
todo = [f"read{{i}}" for i in range(nreads) if f"read{{i}}" not in carried]
stop = os.environ.get("FAKE_DORADO_STOP_AFTER")
names = carried + todo
if stop is not None:
    names = names[: int(stop)]

header = pysam.AlignmentHeader.from_dict(
    {{
        "HD": {{"VN": "1.6", "SO": "unknown"}},
        "PG": [
            {{
                "ID": "basecaller",
                "PN": "dorado",
                "VN": "0.0-stub",
                "CL": "dorado " + " ".join(sys.argv[1:]),
            }}
        ],
    }}
)
rng = random.Random(0)
with pysam.AlignmentFile("-", "wb", header=header) as out:
    for name in names:
        # Long and incompressible on purpose: the wrapper compares checkpoints
        # by byte size, so a record has to outweigh the header text the way a
        # real basecall does.
        seq = "".join(rng.choice("ACGT") for _ in range(2000))
        rec = pysam.AlignedSegment(header)
        rec.query_name = name
        rec.query_sequence = seq
        rec.flag = 4
        rec.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
        out.write(rec)

sys.exit(1 if stop is not None else 0)
'''


@pytest.fixture
def rig(tmp_path, monkeypatch):
    """A stub dorado on PATH plus the paths the wrapper derives from an output."""
    bindir = tmp_path / "bin"
    bindir.mkdir()
    stub = bindir / "dorado"
    stub.write_text(STUB.format(python=sys.executable))
    stub.chmod(0o755)

    argv_log = tmp_path / "argv.txt"
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_DORADO_ARGV", str(argv_log))

    outdir = tmp_path / "bam"
    outdir.mkdir()
    final = outdir / "run.rbc.bam"

    class Rig:
        pass

    rig = Rig()
    rig.final = final
    rig.argv_log = argv_log
    rig.part = outdir / f".{final.name}.partial"
    rig.prev = outdir / f".{final.name}.resume"
    rig.args = outdir / f".{final.name}.args"

    def run(*dorado_args, env=None):
        return subprocess.run(
            ["bash", str(WRAPPER), str(final), *dorado_args],
            capture_output=True,
            text=True,
            env={**os.environ, **(env or {})},
        )

    rig.run = run
    rig.argv = lambda: argv_log.read_text().splitlines()
    return rig


def names_in(bam):
    with pysam.AlignmentFile(str(bam), "rb", check_sq=False) as fh:
        return [r.query_name for r in fh]


DORADO_ARGS = ("--models-directory", "/models", "modelA", "/pod5")


class TestFreshRun:
    def test_produces_output_and_leaves_no_checkpoint(self, rig):
        res = rig.run(*DORADO_ARGS)

        assert res.returncode == 0, res.stderr
        assert names_in(rig.final) == [f"read{i}" for i in range(4)]
        assert not rig.part.exists()
        assert not rig.prev.exists()
        assert not rig.args.exists()

    def test_does_not_pass_resume_from(self, rig):
        rig.run(*DORADO_ARGS)

        assert "--resume-from" not in rig.argv()


class TestKilledRun:
    def test_partial_survives_for_the_next_attempt(self, rig):
        res = rig.run(*DORADO_ARGS, env={"FAKE_DORADO_STOP_AFTER": "2"})

        # The rule fails, so Snakemake reaps {output} -- but the checkpoint is
        # not an output, which is the entire point.
        assert res.returncode != 0
        assert not rig.final.exists()
        assert names_in(rig.part) == ["read0", "read1"]
        assert rig.args.exists()

    def test_second_attempt_resumes_and_finishes(self, rig):
        rig.run(*DORADO_ARGS, env={"FAKE_DORADO_STOP_AFTER": "2"})
        res = rig.run(*DORADO_ARGS)

        assert res.returncode == 0, res.stderr
        argv = rig.argv()
        assert "--resume-from" in argv
        # Rotated: dorado must never read the file it is writing.
        assert argv[argv.index("--resume-from") + 1] == str(rig.prev)
        # Self-contained, and nothing counted twice.
        assert names_in(rig.final) == [f"read{i}" for i in range(4)]
        assert not rig.prev.exists()
        assert not rig.args.exists()

    def test_resume_survives_repeated_kills(self, rig):
        rig.run(*DORADO_ARGS, env={"FAKE_DORADO_STOP_AFTER": "1"})
        rig.run(*DORADO_ARGS, env={"FAKE_DORADO_STOP_AFTER": "3"})
        res = rig.run(*DORADO_ARGS)

        assert res.returncode == 0, res.stderr
        assert names_in(rig.final) == [f"read{i}" for i in range(4)]


class TestRefusesTheWrongCheckpoint:
    def test_changed_arguments_discard_it(self, rig):
        rig.run(*DORADO_ARGS, env={"FAKE_DORADO_STOP_AFTER": "2"})
        res = rig.run("--models-directory", "/models", "modelB", "/pod5")

        # A checkpoint from another basecalling model would be rejected by
        # dorado itself, wedging every retry. It is dropped here instead.
        assert res.returncode == 0, res.stderr
        assert "--resume-from" not in rig.argv()
        assert "discarding checkpoint" in res.stderr
        assert not rig.prev.exists()

    def test_headerless_partial_is_dropped(self, rig):
        rig.run(*DORADO_ARGS, env={"FAKE_DORADO_STOP_AFTER": "2"})
        rig.part.write_bytes(b"\x1f\x8b" + os.urandom(64))  # killed before the header

        res = rig.run(*DORADO_ARGS)

        assert res.returncode == 0, res.stderr
        assert "--resume-from" not in rig.argv()
        assert names_in(rig.final) == [f"read{i}" for i in range(4)]

    def test_never_resumes_from_the_shorter_of_two(self, rig):
        # A job killed inside dorado's copy-through leaves a partial shorter than
        # the file it resumed from. Taking the newest blindly would lose reads.
        rig.run(*DORADO_ARGS, env={"FAKE_DORADO_STOP_AFTER": "3"})
        rig.part.rename(rig.prev)
        rig.run(*DORADO_ARGS, env={"FAKE_DORADO_STOP_AFTER": "1"})
        assert rig.part.stat().st_size < rig.prev.stat().st_size

        res = rig.run(*DORADO_ARGS)

        assert res.returncode == 0, res.stderr
        assert "--resume-from" in rig.argv()
        assert names_in(rig.final) == [f"read{i}" for i in range(4)]


class TestUsage:
    def test_requires_an_output_and_a_command(self, rig):
        res = subprocess.run(
            ["bash", str(WRAPPER), str(rig.final)],
            capture_output=True,
            text=True,
        )

        assert res.returncode == 2
        assert "usage:" in res.stderr
