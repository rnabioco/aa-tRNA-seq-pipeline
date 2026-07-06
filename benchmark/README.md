# Benchmark harness — migration result-parity

Verifies that the **dorado 2.0.1 / escpod / leech** migration (commit `ca528e0`)
did not change pipeline results. It answers two different questions with two
different standards of proof:

| Layer | Change | Expectation | Tool |
|-------|--------|-------------|------|
| **Isolation** | `pod5` CLI → `escpod` (merge/filter) | **bit-lossless** (same reads + signal) | `equivalence/pod5_equivalence.py` |
| **Isolation** | `remora` → `leech` (same `cca_classifier.pt`) | **near-exact** per-read charging calls | `equivalence/classifier_equivalence.py` |
| **End-to-end** | dorado 1.4.0→2.0.1, model v5.0.0→v5.3.0 | basecalls change; **aggregate** results hold | `fingerprint.py` + `compare.py` |

Why split it: the dorado bump legitimately changes basecalls, so an end-to-end
diff can never be bit-identical — only aggregate biological conclusions (mapping
rate, per-tRNA CPM, charged fraction) should survive. The escpod and leech swaps,
by contrast, *should* be lossless, so they get their own exact checks that are
not muddied by basecaller drift.

Everything runs in the `benchmark` pixi environment:

```bash
pixi install            # picks up the new benchmark env (adds pod5, pyyaml)
```

There are two ways to drive the harness: **Snakemake targets** (integrated with
the pipeline DAG, below) or the **standalone scripts / pixi tasks** (sections 1–3
further down). Use whichever fits.

## Snakemake targets (run separately from `all`)

The rules live in `workflow/rules/benchmark.smk` and are **not** part of `rule
all`, so a normal pipeline run never triggers them. Name the target explicitly.
Put the target **before** `--configfile` (Snakemake's `--configfile` greedily
consumes following args):

```bash
# fingerprint this run (builds the pipeline, then snapshots it)
pixi run snakemake benchmark_fingerprint --configfile=config/config-test.yml --cores 4

# isolation checks
pixi run snakemake benchmark_pod5_equivalence_all --configfile=config/config-test.yml --cores 4
pixi run snakemake benchmark_classifier_equivalence_all --profile cluster/slurm   # GPU (leech)
pixi run snakemake benchmark_isolation_all --profile cluster/slurm                 # both isolation checks

# compare two fingerprints (paths from config; see benchmark_compare below)
pixi run snakemake benchmark_compare --configfile=config/config-test.yml --cores 1 \
    --config benchmark_baseline=benchmark/snapshots/old/fingerprint \
             benchmark_candidate=benchmark/snapshots/new/fingerprint \
             benchmark_profile=aggregate
```

Reports land under `{output_directory}/benchmark/`. The isolation and compare
rules always write their report (with a `RESULT: PASS/FAIL` line) rather than
failing the Snakemake job on divergence — read the report to see the verdict. For
CI-style gating on exit code, call `benchmark/compare.py` directly (section 1).

The two ends of an old→new comparison are still two separate pipeline runs, so
generate each fingerprint with its own `benchmark_fingerprint` invocation (one per
ref/config, e.g. via `run_ref.sh`), then diff them with `benchmark_compare`.

## 1. End-to-end: fingerprint + compare

A **fingerprint** is a compact, diffable snapshot of a run's result-bearing
outputs (charging tables, CPM, alignment stats, base-calling errors, modkit
pileups) plus the tool-version manifest. `compare.py` diffs two fingerprints
under a named tolerance profile in `tolerances.yml`:

- **`strict`** — for the isolation swaps (run old vs new with dorado held fixed).
- **`aggregate`** — for the full old→new comparison including the dorado bump.

### Generate the two sides

`run_ref.sh` checks out a git ref into a throwaway worktree, runs the pipeline
into an isolated output dir, and fingerprints it. GPU rules (`rebasecall`,
`classify_charging`) need a GPU — pass a cluster `--profile`.

```bash
# baseline = pre-migration commit
benchmark/run_ref.sh --ref 403e755 --config config/config-test.yml \
    --label old --profile cluster/slurm --setup

# candidate = migrated branch
benchmark/run_ref.sh --ref migrate/dorado-2.0.1-escpod --config config/config-test.yml \
    --label new --profile cluster/slurm --setup
```

`--setup` runs `pixi run setup` inside each worktree so each ref fetches its own
dorado/escpod/leech versions. Snapshots land in `benchmark/snapshots/<label>/`
(git-ignored).

### Compare

```bash
pixi run -e benchmark python benchmark/compare.py \
    benchmark/snapshots/old/fingerprint \
    benchmark/snapshots/new/fingerprint \
    --profile aggregate --json benchmark/snapshots/report.json
```

Exit status is non-zero if any checked metric exceeds tolerance (use in CI).
Metrics missing from either side are **skipped**, not failed. `--report-only`
prints the report without gating.

If you already have two output directories (however produced), skip `run_ref.sh`:

```bash
pixi run bench-fingerprint -- /path/to/old/outputs benchmark/snapshots/old/fingerprint --label old
pixi run bench-fingerprint -- /path/to/new/outputs benchmark/snapshots/new/fingerprint --label new
pixi run bench-compare    -- benchmark/snapshots/old/fingerprint benchmark/snapshots/new/fingerprint --profile aggregate
```

## 2. Isolation: escpod vs pod5

Run both tools on the same input and compare read-id set + raw signal per read.
The `pod5` CLI and the `pod5` python oracle both ship in the benchmark env; only
`escpod` must be added to PATH (via `scripts/setup-env.sh` or `~/.cargo/bin`):

```bash
pixi run -e benchmark bash -c '
  source scripts/setup-env.sh            # or: export PATH="$HOME/.cargo/bin:$PATH"

  # merge losslessness
  python benchmark/equivalence/pod5_equivalence.py merge --workdir benchmark/snapshots/pod5eq \
      .tests/sample1/pod5_pass/1.pod5 .tests/sample1/pod5_pass/2.pod5

  # filter losslessness
  python benchmark/equivalence/pod5_equivalence.py filter --ids read_ids.txt \
      --workdir benchmark/snapshots/pod5eq .tests/sample1/pod5_pass/1.pod5
'
```

Verified on the test data: `escpod merge` of the two sample1 pod5s is
signal-identical to `pod5 merge` for all 200 reads.

## 3. Isolation: leech vs remora

Both engines load the same `cca_classifier.pt`; charging scores should agree.

```bash
# run both engines on a tagged BAM + pod5, then diff (needs GPU for leech)
pixi run bench-classifier-eq -- run \
    --pod5  <output>/pod5/sample1/sample1.pod5 \
    --bam   <output>/bam/tagged/sample1/sample1.tagged.bam \
    --model resources/models/cca_classifier.pt \
    --workdir benchmark/snapshots/clfeq

# or diff two already-classified charging BAMs
pixi run bench-classifier-eq -- compare --leech leech.charging.bam --remora remora.charging.bam
```

Reports charged/uncharged call agreement (threshold 200) and the `|Δcl|`
distribution; fails if agreement < 0.995, max `|Δcl|` > 5, or read overlap < 0.99.

## Files

```
benchmark/
├── README.md
├── tolerances.yml                     # strict / aggregate profiles
├── lib.py                             # output-schema readers
├── fingerprint.py                     # OUTPUT_DIR -> fingerprint/
├── compare.py                         # two fingerprints -> report + exit code
├── run_ref.sh                         # run a git ref -> snapshot
└── equivalence/
    ├── pod5_equivalence.py            # escpod vs pod5 (exact)
    └── classifier_equivalence.py      # leech vs remora (near-exact)
```

## Tuning tolerances

Edit `tolerances.yml`. Set any field to `null` to disable that check. After a
first real old→new run, inspect the `aggregate` report and tighten limits to
just above the observed drift so future regressions stand out.
```
