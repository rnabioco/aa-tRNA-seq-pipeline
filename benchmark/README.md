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

Everything runs in the **default pixi env** — `pixi run setup` installs `pod5`
(the escpod oracle) and the analysis deps (pandas/pysam/pyarrow/pyyaml) are
already there. No separate environment.

There are two ways to drive the harness: **Snakemake targets** (integrated with
the pipeline DAG, below) or the **standalone scripts / pixi tasks** (sections 1–4
further down). Use whichever fits.

For comparing **dorado versions specifically** (the change that actually moves
results), the harness installs multiple dorado versions side-by-side and compares
their outputs directly — no git worktrees needed. See section 4.

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

# dorado version comparison: basecall each sample's pod5 with every configured
# dorado version and diff the basecalls (GPU). See section 4.
pixi run snakemake benchmark_basecall_compare_all --profile cluster/slurm

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
pileups), the tool-version manifest, and **per-rule wall time + peak RSS**
(from Snakemake `benchmark:` directives on merge_pods, rebasecall, and the
classifiers). `compare.py` diffs two fingerprints under a named tolerance
profile in `tolerances.yml`, and prints an informational runtime section
(baseline → candidate seconds, % faster/slower) so you can see whether the
migration sped up or regressed each step. The leech-vs-remora check also reports
each engine's wall time (GPU vs CPU) directly.

- **`strict`** — for the isolation swaps (run old vs new with dorado held fixed).
- **`aggregate`** — for the full old→new comparison including the dorado bump.

### Generate the two sides

**By dorado version (recommended for tool-version parity):** run the pipeline
once per version into separate output dirs, then fingerprint each. Same code,
same env — only the basecaller changes.

Use `dorado_opts_override="--emit-moves"` (canonical basecalls, no
`--modified-bases`) for BOTH runs. Modification models differ across model
versions (e.g. v5.0.0 has no 2′-O-methyl mods), so mods are not comparable
across versions — and they don't affect charging or alignment anyway (same
simplex sequence/moves). The comparison then covers the mod-independent
readouts: charging, CPM, mapping rate, base-call errors.

```bash
# old basecaller (canonical)
pixi run snakemake benchmark_fingerprint \
    --configfile=config/config-test.yml --profile cluster/slurm \
    --config output_directory=.tests/cmp-1.4.0 \
             dorado_version=1.4.0 dorado_model=rna004_130bps_sup@v5.0.0 \
             base_calling_model=resources/models/rna004_130bps_sup@v5.0.0 \
             dorado_opts_override=--emit-moves benchmark_label=v1.4.0
# new basecaller (canonical; version/model default to config)
pixi run snakemake benchmark_fingerprint \
    --configfile=config/config-test.yml --profile cluster/slurm \
    --config output_directory=.tests/cmp-2.0.1 \
             dorado_opts_override=--emit-moves benchmark_label=v2.0.1
```

**By git ref (for arbitrary code changes, not just tool versions):**
`run_ref.sh` checks out a ref into a throwaway worktree, runs the pipeline, and
fingerprints it. Heavier (worktree + optional per-ref `--setup`); use it when the
difference is in pipeline *code*, not just the dorado version.

```bash
benchmark/run_ref.sh --ref 403e755 --config config/config-test.yml \
    --label old --profile cluster/slurm --setup
```

Snapshots land in `benchmark/snapshots/<label>/` (git-ignored).

### Compare

```bash
pixi run python benchmark/compare.py \
    .tests/cmp-1.4.0/benchmark/fingerprint \
    .tests/cmp-2.0.1/benchmark/fingerprint \
    --profile aggregate --json report.json
```

Modkit is skipped (no shared modified sites under canonical basecalls); the
report covers charging, CPM, mapping, and base-call errors.

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
The `pod5` CLI and the `pod5` python oracle are in the default env; only `escpod`
must be on PATH (via `scripts/setup-env.sh` or `~/.cargo/bin`):

```bash
source scripts/setup-env.sh            # or: export PATH="$HOME/.cargo/bin:$PATH"

# merge losslessness
pixi run python benchmark/equivalence/pod5_equivalence.py merge \
    --workdir benchmark/snapshots/pod5eq \
    .tests/sample1/pod5_pass/1.pod5 .tests/sample1/pod5_pass/2.pod5

# filter losslessness
pixi run python benchmark/equivalence/pod5_equivalence.py filter --ids read_ids.txt \
    --workdir benchmark/snapshots/pod5eq .tests/sample1/pod5_pass/1.pod5
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

## 4. Basecaller: compare dorado versions

Config `dorado_variants` lists the `(version, model)` pairs to install and
compare; `pixi run setup` installs each dorado binary + base model side-by-side
under `resources/tools/dorado/<version>/`:

```yaml
dorado_variants:
  - version: "1.4.0"
    model: "rna004_130bps_sup@v5.0.0"
  - version: "2.0.1"
    model: "rna004_130bps_sup@v5.3.0"
```

`benchmark_basecall_compare_all` basecalls each sample's merged pod5 with every
configured version (canonical basecalls, no modified bases) and diffs them per
sample — read-id overlap, exact sequence match, mean similarity, and length/quality
deltas — writing `{output_directory}/benchmark/basecall/<sample>.compare.txt`:

```bash
pixi run snakemake benchmark_basecall_compare_all --profile cluster/slurm
```

Or diff already-basecalled uBAMs directly:

```bash
pixi run bench-basecall-compare -- --bam 1.4.0=a.ubam --bam 2.0.1=b.ubam
```

This isolates the basecaller's sequence/quality change. Whether those changes move
the biological conclusions (charging, CPM, mods) is the end-to-end fingerprint
comparison in section 1 — run per version via config override (see "Generate the
two sides").

## Files

```
benchmark/
├── README.md
├── tolerances.yml                     # strict / aggregate profiles
├── lib.py                             # output-schema readers
├── fingerprint.py                     # OUTPUT_DIR -> fingerprint/
├── compare.py                         # two fingerprints -> report + exit code
├── basecall_compare.py                # dorado version basecall diff
├── run_ref.sh                         # run a git ref -> snapshot (arbitrary code changes)
└── equivalence/
    ├── pod5_equivalence.py            # escpod vs pod5 (exact)
    └── classifier_equivalence.py      # leech vs remora (near-exact)
```

## Tuning tolerances

Edit `tolerances.yml`. Set any field to `null` to disable that check. After a
first real old→new run, inspect the `aggregate` report and tighten limits to
just above the observed drift so future regressions stand out.
```
