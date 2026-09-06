"""
Tests for workflow/scripts/escpod_classify_fallback.sh.

The script retries `escpod classify` with a supplied move-table frame when
detection came back underpowered, and MUST NOT retry on anything else -- a
second attempt that hides an OOM or a missing POD5 would turn a loud failure
into a quiet one.

Driven against a fake `escpod` on PATH, so the real binary, a GPU and a POD5
are all unnecessary.
"""

import os
import shutil
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).parent.parent.parent
SCRIPT = REPO_ROOT / "workflow" / "scripts" / "escpod_classify_fallback.sh"

UNDERPOWERED = (
    "Error: orientation check underpowered: 19 informative reads (need 50); "
    "pass --orientation to override for small batches"
)


def fake_escpod(tmp_path, behaviour):
    """Put a fake `escpod` on PATH that fails in a chosen way.

    It appends every invocation to calls.txt, so a test can tell one attempt
    from two without parsing the log.
    """
    bindir = tmp_path / "bin"
    bindir.mkdir(exist_ok=True)
    calls = tmp_path / "calls.txt"
    (bindir / "escpod").write_text(
        "#!/usr/bin/env bash\n"
        f'echo "$@" >>"{calls}"\n'
        f"{behaviour}\n"
    )
    (bindir / "escpod").chmod(0o755)
    env = dict(os.environ, PATH=f"{bindir}:{os.environ['PATH']}")
    return env, calls


# Fails only while no frame is forced -- exactly escpod's real behaviour.
RECOVERS = (
    'if [[ "$*" == *--orientation* ]]; then echo "classified"; exit 0; fi\n'
    f'echo "{UNDERPOWERED}"; exit 1'
)
ALWAYS_OOM = 'echo "Error: out of memory allocating arena"; exit 1'
ALWAYS_OK = 'echo "classified"; exit 0'


def run(tmp_path, behaviour, fallback):
    env, calls = fake_escpod(tmp_path, behaviour)
    log = tmp_path / "logs" / "sample1"
    proc = subprocess.run(
        ["bash", str(SCRIPT), fallback, str(log), "--bam", "s.bam"],
        env=env,
        capture_output=True,
        text=True,
    )
    attempts = calls.read_text().splitlines() if calls.exists() else []
    return proc.returncode, attempts, (log.read_text() if log.exists() else "")


@pytest.mark.skipif(shutil.which("bash") is None, reason="needs bash")
class TestRetriesOnlyWhenUnderpowered:
    def test_underpowered_recovers_with_the_fallback(self, tmp_path):
        rc, attempts, log = run(tmp_path, RECOVERS, "reversed")
        assert rc == 0
        assert len(attempts) == 2
        assert "--orientation reversed" in attempts[1]

    def test_underpowered_without_a_fallback_still_fails(self, tmp_path):
        """`orientation_fallback: none` restores the v0.7.2 behaviour."""
        rc, attempts, _ = run(tmp_path, RECOVERS, "")
        assert rc == 1
        assert len(attempts) == 1

    def test_a_real_failure_is_never_retried(self, tmp_path):
        """The property that matters: an OOM must not be masked by a retry."""
        rc, attempts, _ = run(tmp_path, ALWAYS_OOM, "reversed")
        assert rc == 1
        assert len(attempts) == 1

    def test_success_runs_once(self, tmp_path):
        rc, attempts, _ = run(tmp_path, ALWAYS_OK, "reversed")
        assert rc == 0
        assert len(attempts) == 1


@pytest.mark.skipif(shutil.which("bash") is None, reason="needs bash")
class TestTheRunRecordsWhatHappened:
    def test_a_supplied_frame_is_announced_in_the_log(self, tmp_path):
        """Otherwise a forced call is indistinguishable from a detected one."""
        _, _, log = run(tmp_path, RECOVERS, "reversed")
        assert "SUPPLIED rather than detected" in log
        assert "--orientation reversed" in log
        assert "depth guard" in log

    def test_first_attempt_output_is_kept(self, tmp_path):
        _, _, log = run(tmp_path, RECOVERS, "reversed")
        assert "orientation check underpowered" in log

    def test_nothing_is_announced_when_no_retry_happened(self, tmp_path):
        _, _, log = run(tmp_path, ALWAYS_OK, "reversed")
        assert "SUPPLIED" not in log


@pytest.mark.skipif(shutil.which("bash") is None, reason="needs bash")
class TestUsage:
    def test_too_few_arguments_is_a_usage_error(self, tmp_path):
        proc = subprocess.run(
            ["bash", str(SCRIPT), "reversed"], capture_output=True, text=True
        )
        assert proc.returncode == 2
        assert "usage" in proc.stderr.lower()

    def test_log_directory_is_created(self, tmp_path):
        _, _, log = run(tmp_path, ALWAYS_OK, "reversed")
        assert log == "classified\n"
