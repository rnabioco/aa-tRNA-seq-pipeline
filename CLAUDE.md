# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Snakemake pipeline for processing Oxford Nanopore Technologies (ONT) aa-tRNA-seq data. The pipeline distinguishes between charged (aminoacylated) and uncharged tRNA molecules using a machine learning model — run by `escpod classify` — trained on nanopore signal data over the CCA 3' end of tRNA molecules.

## Setup and Environment

### Initial Setup

```bash
# Install all dependencies
pixi install

# One-time setup: downloads dorado, basecalling models, escpod, and WarpDemuX
# IMPORTANT: Run this once before using the pipeline, from a single node only
pixi run setup

# Download test data (optional, for testing only)
pixi run dl-test-data
```

**Note:** The `pixi run setup` command installs tools that are not available via conda (dorado, escpod, WarpDemuX). Run this once from a single node before submitting cluster jobs to avoid race conditions on shared filesystems.

### Running the Pipeline

```bash
# Dry run with test config
pixi run dry-run

# Run locally with test data (4 cores)
pixi run test

# Run on LSF cluster
pixi run test-lsf

# Run on SLURM cluster
pixi run test-slurm

# Run preprint pipeline on cluster
pixi run run-preprint
```

### Direct Snakemake Commands

```bash
pixi run snakemake --configfile=config/config-test.yml --cores 8
```

### Cluster Execution

The pipeline supports both LSF and SLURM schedulers. Key files:

**LSF:**
- `run-test.sh`: Test data execution on LSF
- `run-preprint.sh`: Full preprint data execution on LSF
- `cluster/lsf/config.yaml`: LSF-specific resource configurations

**SLURM:**
- `cluster/slurm/config.yaml`: SLURM-specific resource configurations (customize partition/account for your cluster)

GPU-intensive rules (rebasecall) automatically request GPU resources via queue/partition configuration. `classify_charging` is CPU-only — `escpod classify` has no GPU path.

## Architecture

### Snakemake Workflow Structure

The workflow is modular with rules split across multiple files:

```
workflow/
├── Snakefile                          # Main entry point, includes all rule modules
├── rules/
│   ├── common.smk                     # Sample parsing, helper functions, outputs definition
│   ├── aatrnaseq-process.smk          # Core processing: pod5 merge → basecalling → alignment
│   ├── aatrnaseq-charging.smk         # Charging classification outputs (2 rules)
│   ├── aatrnaseq-qc.smk               # QC metrics: base calling errors, alignment stats (3 rules)
│   ├── aatrnaseq-modifications.smk    # Modification calling: coverage, modkit outputs (4 rules)
│   └── demux.smk                      # Barcode demux (WDX or LDX) + EDX concordance (conditionally loaded)
├── scripts/                           # Python scripts called by rules
│   ├── detect_3p_adapters.py           # Detect 3' adapter identity per read (for EDX splitting)
│   └── generate_squiggy_session.py    # Generate Squiggy/Positron session JSON
└── envs/
    └── aatrnaseqpipe-env.yml          # Conda environment (legacy)
```

**Key Architectural Details:**

- **Sample Management**: `workflow/rules/common.smk` contains `parse_samples()` which reads `config/samples.tsv` and `find_raw_inputs()` which recursively searches for pod5 files in specified directories
- **Tool Management**: Modkit is managed by pixi. Dorado and escpod are downloaded by `pixi run setup` into `resources/tools/`, and put on PATH by the Snakefile's `onstart` shell prefix
- **Output Aggregation**: `pipeline_outputs()` in `common.smk` defines all final output files for the `rule all` target

### Pipeline Flow

```
POD5 files → stage_pod5 (symlinks) → rebasecall (Dorado) → bwa_align (dorado tags carried through) →
classify_charging (escpod) → add_adapter_tags → finalize_bam → Summary tables
```

On an LDX run the first two steps are replaced (there is no per-sample POD5 to
stage or basecall), and the flow rejoins at `bwa_align`:
```
raw POD5 → escapepod_demux (--annotate, writes .p5s) → rebasecall_ldx_run (whole run, one dorado pass)
         → ldx_split_parent_map + extract_ldx_sample_reads → split_ldx_ubam → (as above)
```

For EDX samples (dual barcoding), 3' adapter detection happens on the uBAM before alignment, and the resulting read-id list is the whole filter: `bwa_align` aligns only those reads (`samtools view -N`), and `classify_charging` only ever touches reads the BAM names, so no filtered FASTQ or POD5 is written:
```
rebasecall → detect_edx_adapters → extract_edx_read_ids → bwa_align → classify_charging → ...
```

### Core Processing Pipeline (aatrnaseq-process.smk)

1. **stage_pod5**: Lay a directory of symlinks over the sample's raw POD5 files (`pod5/{sample}/<run>/<pod5_pass|pod5_fail|pod5>/`). Nothing copies the signal: dorado and `escpod classify` both take a directory. Replaced `merge_pods`, which wrote a full second copy of every run per sample
2. **rebasecall**: Use dorado (`--recursive` over the staged directory) to rebasecall with move tables (required by the charging model). Runs through `workflow/scripts/dorado_basecall_resume.sh`, as does `rebasecall_ldx_run`: the partial BAM of a killed attempt is kept beside the output (where Snakemake will not reap it) and fed back as dorado's `--resume-from`, so an overrun costs the tail of a basecall rather than the whole flowcell. `escapepod_demux` has no equivalent and does not resume
3. **bwa_align**: Stream the uBAM through `samtools fastq -T '*'` into `bwa mem -C`, so dorado's tags (the `mv`/`ns`/`ts` move table, MM/ML modbase calls, RG) ride the FASTQ comment onto the aligned records. `-H` inserts the uBAM's own @RG lines with SM/LB/BC stamped (`stamp_read_groups.py`), and on demux runs a constant `BC` tag goes onto every read. No FASTQ is written and there is no tag-injection step (`inject_ubam_tags` is retired). `-K 100000000` pins bwa's batch so the tag payload (~13x the read) costs a few GB, not the OOM PR #86 saw with the default per-thread batch. Output is primary forward alignments (`-F 2324`)
4. **classify_charging**: Run `escpod classify` to classify charged vs uncharged reads. Writes a `cl` tag onto the records it scored and passes every other record through unchanged, so dorado's MM/ML modbase tags survive and no tag round-trip is needed. Also emits a per-read calls TSV with a `reason` for every read it did not score. It is handed the sample's whole signal store (staged directory, WarpDemuX split POD5, or the raw LDX run) and looks reads up by the BAM's ids, so no EDX-filtered POD5 exists
5. **add_adapter_tags**: Detect adapter positions and add pt tags with 5'/3' boundaries
6. **finalize_bam**: Hardlink adapter-tagged BAM as final output

### Summary Generation

After classification, generates (split across three rule files):

**aatrnaseq-charging.smk:**
- Charging probability tables (`cl` tag values per read)
- CPM (counts per million) for charged/uncharged tRNA
- Per-read charging calls with a no-call `reason` (`{sample}.charging_calls.tsv.gz`)

**aatrnaseq-qc.smk:**
- Base calling error frequencies
- Alignment statistics (the `classified` row counts reads carrying a `cl` tag)
- CCA anchor coverage
- Read attrition (`summary/read_attrition.tsv.gz`), which folds in the classifier's no-call reasons

**aatrnaseq-modifications.smk:**
- Coverage bedGraph files (counts and CPM)
- Modkit modification pileups
- Modkit per-read modification calls

**common.smk:**
- Squiggy session JSON (`squiggy-session.json`) for loading outputs in Positron IDE

## Configuration

### Writing a new run's config

Do not hand-write the pair. `pixi run new-run-config -- --name <name> --run <dir>
--sample <name>=<spec> ... --reference-raw <fa> --three-prime <name>=<seq>
--dry-run` writes `config/config-<name>.yml` and its samples file with every
override explained, and refuses what the pipeline would refuse (codes not
emitted by the vendored bundles, `fdx` without `ldx`, overlapping barcode
tuples, `edx` names not in `adapters.three_prime`, 3' adapters of unequal
length or not starting with GGC, run dirs without POD5, missing references).
`pixi run check-run-config -- config/config-<name>.yml --dry-run` re-applies the
same checks to an existing pair. The sample spec grammar and the questions to
ask first are in `.claude/skills/new-run-config/SKILL.md`; the logic is
`workflow/scripts/run_config.py`.

### Main Config Files

- `config/config-base.yml`: Base configuration included by Snakefile
  - Base calling model path
  - Reference fasta
  - Charging model bundle and operating point (`charging`)
  - Dorado and escpod versions for download
  - Command-line options for tools (dorado, bwa, filters)

- `config/samples.tsv`: Two-column TSV (no header)
  - Column 1: Unique sample ID
  - Column 2: Path to sequencing run folder containing pod5_pass/pod5_fail/pod5 subdirectories

- `config/config-test.yml`: Overrides base config for test data

### Important Config Parameters

- **opts.bam_filter**: Controls full-length read filtering (`-5 24 -3 23 -s` requires 24bp 5' adapter, 23bp 3' adapter, positive strand)
- **opts.dorado**: Includes `--modified-bases m5C_2OmeC inosine_m6A_2OmeA pseU_2OmeU 2OmeG --emit-moves` for modification calling and move tables
- **opts.bwa**: RNA-optimized alignment parameters (`-W 13 -k 6 -T 20 -x ont2d`)
- **charging.ml_threshold**: 200 (200-255 = charged, <200 = uncharged). Config, not hardcoded. It is the bundle's declared operating point, measured against ligation chemistry at FPR 1.74% / TPR 0.939; its precision depends on the sample's own charged fraction, so low-charging samples need a higher value (see `docs/troubleshooting/faq.md`)
- **cleanup_intermediates**: Opt-in auto-deletion of large regenerable intermediates during a run, via `temp()`. Accepts a bool or a list of tier names (`cascade`, `basecall`, `demux_scratch`, `split_pod5`) resolved by `maybe_temp()` / `_enabled_cleanup_tiers()` in `common.smk`; the retired `fastq` and `merged_pod5` tiers are accepted and ignored (nothing writes a FASTQ or a merged POD5 any more), and an unknown name is an error. `bam/final` is always kept. `split_pod5` deletes the WarpDemuX split POD5 once classification and the session file are done; it is that path's classification input, so re-running classification afterwards means re-splitting from the raw run. Unbarcoded and LDX samples have no `demux/pod5` (their store is the raw run, via symlinks for unbarcoded samples), so the tier is inert there. The on-demand `clean` rule (`rules/clean.smk`) remains the catch-all superset for reclaiming space on already-completed runs. See `config/README.md` for tier→directory mapping.

## Demultiplexing (Optional)

The pipeline supports optional signal-level barcode demultiplexing for
pooled/multiplexed runs, via one of two mutually exclusive backends:

| Backend | Config key | Barcodes | Tool | Shape |
|---|---|---|---|---|
| WarpDemuX | `warpdemux.enabled` | WDX (`barcode04`) | `warpdemux` | Classify to a table → parse to a read→barcode mapping → `escpod filter` per sample → basecall each split POD5 |
| escapepod | `ldx.enabled` | LDX (`ldx01`) | `escpod demux --annotate` | One fused pass records each read's barcode in a `.p5s` sidecar beside the raw POD5 → basecall the run once → split the uBAM per sample |
| escapepod, second axis | `fdx.enabled` (needs `ldx.enabled`) | FDX (`fdx01`), the 5' index | `escpod demux --model fdx=…` | A dual-index sample is `{ldx: ldx01, fdx: fdx01}` and owns the reads on which BOTH calls agree; `select_demux_reads.py` joins the axes per read |

Both converge on the same per-sample uBAM (`bam/rebasecall/{sample}/`), and
everything downstream is identical. Enabling both is rejected at parse time.

LDX writes no POD5 of its own: the raw run plus a few-MB sidecar per POD5
directory is the whole signal store, where the WarpDemuX path leaves a second
full copy of every read. The classifiers are handed the raw run directory —
`escpod classify` walks a directory recursively — so a run split across
pod5_pass/pod5_fail needs no merged copy either.
That is also why the LDX path basecalls per run rather than per sample — there
is no per-sample POD5 to hand dorado — and why the signal classifiers are given
the raw run and rely on the per-sample BAM to bound what they read.

LDX is the successor path; see `config/README.md` for the LDX sample-file
format, the self-describing model bundle, the sidecar's semantics, and its CPU
cost. The vendored models live in `resources/models/demux/` (see the README
there for why they are committed rather than fetched).

**FDX** is a 5' index on the same molecule (a 27-nt code in a 77-nt DNA adapter,
read off the READ END, where LDX is read off the adapter boundary). It is a
second axis of the escpod path, not a backend: `fdx.enabled` requires
`ldx.enabled`, a sample names both codes, and `select_demux_reads.py` assigns a
read to a sample only when every axis agrees. escpod can call both models in
one sweep (`fdx.fused: true`); escpod 0.19.0 and 0.20.0 refused the LDX
boundary flags whenever a read-end model was in the run. 0.21.0 scopes them per
head (escapepod-rs#323), so fusing is legal on the pinned version, but the
default stays a second pass (`escapepod_demux_fdx`) until it is measured here —
see the `fdx` block in `config-base.yml`. The FDX
gate (3.5 nats, the bundle's declared operating point) is load-bearing: the CRF
snaps reads without a 5' index onto a code, and only the gate refuses them.

### Testing the demux path

Only the LDX backend has an end-to-end test, because it is the only one with a
committed barcoded fixture (`.tests/fixtures/ldx-demux`, 415 reads of a real
pooled run — see the README there):

```bash
pixi run dry-run-ldx      # DAG only, no GPU, runs in CI
pixi run test-ldx         # full run against the fixture (needs a GPU for dorado)
pixi run test-ldx-resume  # basecall resume, #149 (needs a GPU; rebuilds .tests/outputs-ldx)
pixi run dry-run-fdx      # dual-index (LDX + FDX) DAG, runs in CI
pixi run test-fdx         # dual-index run against .tests/fixtures/fdx-demux (needs a GPU)
```

The FDX fixture (`.tests/fixtures/fdx-demux`, 365 reads) is Run2 of the FDX
pilot, which the shipped fdx4 model was NOT trained on, so its per-sample
recovery is a held-out number; `tests/integration/test_fdx_demux.py` asserts
its shape.

`test-ldx-resume` (`.tests/run_ldx_resume_e2e.sh`) is the only thing that
exercises `dorado_basecall_resume.sh` against real dorado — run it after any
change to that wrapper. Read its header before editing it: the fixture
basecalls faster than dorado starts, so the kill and the resume have to be
proved separately, and it carries a **control** run whose only job is to
measure the pipeline's own nondeterminism (bwa picks arbitrarily among
equal-scoring near-identical tRNA isodecoders, which moves per-reference
summary rows) so that instability is never mistaken for a resume bug. The
header also records two traps that make a naive version of this test pass while
testing nothing: SIGTERM to snakemake alone makes it *wait* for the running job
rather than kill it, and deleting a rule's output does not make it re-run when
its consumers are up to date.

`config/config-demux-test.yml` (WarpDemuX) **cannot complete a run** and is a
dry-run target only: it points at unbarcoded sacCer3 data relabelled as
barcoded, so adapter detection finds nothing and bwa maps none of the reads
WarpDemuX routes to those barcodes. See issue #120 and that file's header.

The rest of this section describes the WarpDemuX backend.

### Enabling Demultiplexing

1. **Install WarpDemuX**: `pixi run setup` (installs WarpDemuX along with other tools)
2. **Create YAML sample file** with barcode assignments (see `config/samples-demux-example.yml`)
3. **Enable in config**: Set `warpdemux.enabled: true`

### Sample File Formats

**TSV format (existing, no demux):**
```
sample1    /path/to/run1
sample2    /path/to/run2
```

**YAML format (with demux):**
```yaml
runs:
  - path: /path/to/pooled/run
    barcode_kit: "WDX4_rna004_v1_0"
    samples:
      charged_sample: "barcode04"
      uncharged_sample: "barcode05"
```

**YAML format (dual barcoding — WDX + EDX):**
```yaml
runs:
  - path: /path/to/pooled/run
    barcode_kit: "WDX4_tRNA_rna004_v1_0"
    samples:
      sample_bc03:
        wdx: "barcode03"
        edx: "edx01"       # must match adapter name from adapters.three_prime config
      sample_bc04:
        wdx: "barcode04"
        edx: "edx02"       # must match adapter name from adapters.three_prime config
```

EDX values must directly match adapter names from `adapters.three_prime` in the config (e.g., `edx01`, `edx02`). EDX splitting happens before alignment — the pipeline detects 3' adapter identity on the uBAM, then aligns only the matching reads; the classifier reads only what the BAM names, so no filtered FASTQ or POD5 is written.

When `edx.enabled: true` in config, the `edx_concordance` rule produces `summary/edx/edx_concordance.tsv.gz` — a concordance table of WDX assignment vs EDX (3' adapter) identity per sample.

### Running with Demux

```bash
# Dry run with demux config
pixi run snakemake -n --configfile=config/config-demux-test.yml

# Execute with demux
pixi run snakemake --configfile=config/config-demux-test.yml --cores 8
```

## Charged vs Uncharged Classification

The pipeline runs `escpod classify` against a vendored ONNX model bundle:

- **Model Location**: `charging.model` config parameter — a bundle **directory**,
  `resources/models/charging/charging_feature_nn_sup6_rna004@v0.1.0`. It is
  self-describing (anchor, feature recipe, k-mer table pinned by sha256, abstain
  rule, operating point), so nothing about the recipe is passed as a flag.
  See the README beside it.
- **Anchor**: the CCA|adapter junction in reference coordinates — the adapter G at
  `index(CCAGGC) + 3` — mapped through the move table. Features run over offsets
  -8..+24 around it
- **`cl` Tag**: `round(P(charged) * 255)`, 0-255, written directly onto the aligned
  records. There is no `cm` tag
- **Threshold**: `cl` >= 200 = charged (`charging.ml_threshold`, matching the bundle's
  declared `operating_point.cl`)
- **Abstention**: reads where the aligner placed no common-arm base get **no `cl`
  tag**, not a default class. This is charging-correlated, so a charging fraction
  over called reads alone is an UNDERESTIMATE — always report the no-call rate
  beside it (`{sample}.charging_calls.tsv.gz`, `read_attrition.tsv.gz`)
- **Runtime pinning**: the bundle needs escpod >= 0.19.0; older refuses it with
  ``unknown field `basecaller` `` (the schema is `deny_unknown_fields` on purpose).
  `escpod_version` is pinned alongside the model, and the two move together
- **Basecaller pairing**: the bundle names the basecalls it was trained on —
  `rna004_sup@v6.0.0`, dorado 2.1.1+d66c17c, matching `base_calling_model` and the
  installed binary byte-for-byte. escpod states this at load and does NOT check it,
  so the pipeline does: `charging.basecaller_check` errors on a model mismatch at
  DAG construction, and `pixi run verify-basecaller` proves byte identity. A
  different basecalling model is a domain shift on the k-mer residual — arm-to-arm
  contrasts survive it, absolute charged fractions do not
- **Two bundles are vendored, one per basecalling model** and they are NOT
  interchangeable: `charging_feature_nn_sup6_rna004@v0.1.0` for `rna004_sup@v6.0.0`
  (the default), `charging_feature_nn_rna004@v0.1.1` for data already basecalled with
  `rna004_130bps_sup@v5.3.0`. Switching means setting `base_calling_model`,
  `dorado_model` and `charging.model` together

## Development

### Adding New Rules

When adding new Snakemake rules:
- Place processing rules in `aatrnaseq-process.smk`
- Place charging analysis rules in `aatrnaseq-charging.smk`
- Place QC/statistics rules in `aatrnaseq-qc.smk`
- Place modification/coverage rules in `aatrnaseq-modifications.smk`
- Add helper functions to `common.smk`
- Reference Python scripts should go in `workflow/scripts/`
- Update `pipeline_outputs()` if rule produces final outputs

### Testing Changes

```bash
# Always test with dry run first
pixi run dry-run

# Run specific rule
pixi run snakemake <rule_name> --configfile=config/config-test.yml

# Force rerun of specific rule
pixi run snakemake <rule_name> --forcerun <rule_name> --configfile=config/config-test.yml

# Demultiplexing: dry-run the DAG, then run it against the committed fixture
pixi run dry-run-ldx
pixi run test-ldx

# Basecall resume, after any change to dorado_basecall_resume.sh
pixi run test-ldx-resume
```

A dry-run only builds the DAG — it never looks at the data, which is why
`config-demux-test.yml` passed `dry-run-demux` for as long as it did while being
unable to execute (#120). `tests/integration/test_ldx_demux.py` closes that gap:
its fixture and config tests run without a GPU and assert that the configs and
the fixture actually agree, while the output tests skip unless `.tests/outputs-ldx`
exists from a real run.

### Cluster Resource Configuration

**LSF** - Modify `cluster/lsf/config.yaml` to adjust:
- Memory requirements per rule (mem_mb)
- GPU queue assignments (lsf_queue)
- LSF project tags (lsf_project)
- Maximum concurrent jobs

Rules requiring GPU (rebasecall) must set:
- lsf_queue: "gpu"
- lsf_extra: "-gpu num=1:j_exclusive=yes"
- ngpu: 1

**SLURM** - Modify `cluster/slurm/config.yaml` to adjust:
- Memory requirements per rule (mem_mb)
- GPU partition (slurm_partition)
- Account/allocation (slurm_account)
- Runtime limits (runtime, in minutes)
- Maximum concurrent jobs

Rules requiring GPU (rebasecall) must set:
- slurm_partition: "gpu" (or your cluster's GPU partition)
- gres: "gpu:1"

## Important Notes

- The pipeline requires Snakemake 8.0+
- Modkit is managed by pixi; dorado and escpod are downloaded by `pixi run setup`
- Pins are checked against upstream weekly (`.github/workflows/currency.yml`, by hand `pixi run check-currency`): vendored model bundles against escapepod-models releases, `escpod_version`/`dorado_version` against escapepod-rs/dorado releases, and the escpod checksums in `scripts/setup-tools.sh` against the pinned release's SHA256SUMS. A red run is a decision, not an order to bump: hold a pin by recording why in `resources/models/pins.yml` (families at the top level, tools under `tools:`)
- The charging and demux model bundles are vendored under `resources/models/`
  rather than fetched: upstream `rnabioco/escapepod-models` is private and
  compute nodes have no route to GitHub
- The pipeline tracks git commit ID for reproducibility (see `get_pipeline_commit()`)
- CUDA_VISIBLE_DEVICES is passed through to dorado if set
- Pod5 files are searched recursively in pod5_pass/pod5_fail/pod5 subdirectories
- The charging threshold lives in `charging.ml_threshold` and must match the model bundle's declared operating point
- `shell.prefix()` in the Snakefile's `onstart` REPLACES snakemake's default prefix, which is bash strict mode — the `set -euo pipefail` in it is load-bearing, do not drop it

## Key Outputs

Outputs go to directory specified by `output_dir` in config. Test outputs: `.tests/outputs/`

Key outputs per sample:
- `summary/tables/{sample}/{sample}.charging.cpm.tsv.gz` - CPM-normalized charging counts
- `summary/tables/{sample}/{sample}.charging_prob.tsv.gz` - Per-read charging probabilities
- `summary/tables/{sample}/{sample}.charging_calls.tsv.gz` - Per-read calls, with a `reason` for every read the model did not score
- `bam/final/{sample}/{sample}.bam` - Final BAM with cl (charging), pt (adapter positions) and BC (barcode, demux runs only) tags, plus dorado's MM/ML modbase tags and a valid @RG

Pipeline-level outputs:
- `squiggy-session.json` - Squiggy session file for loading samples in Positron
