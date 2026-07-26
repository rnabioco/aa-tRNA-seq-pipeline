# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Snakemake pipeline for processing Oxford Nanopore Technologies (ONT) aa-tRNA-seq data. The pipeline distinguishes between charged (aminoacylated) and uncharged tRNA molecules using Remora machine learning models trained on nanopore signal data over the CCA 3' end of tRNA molecules.

## Setup and Environment

### Initial Setup

```bash
# Install all dependencies
pixi install

# One-time setup: dorado + basecalling models, remora, escpod (built from source
# with demux), the escapepod python package, leech, demux models, and WarpDemuX
# IMPORTANT: Run this once before using the pipeline, from a single node only
pixi run setup

# Download test data (optional, for testing only)
pixi run dl-test-data
```

**Note:** The `pixi run setup` command installs tools that are not available via conda. Run this once from a single node before submitting cluster jobs to avoid race conditions on shared filesystems.

Setup has two host prerequisites beyond pixi:

- **Rust toolchain >= 1.95** — `escpod` is built from source (`scripts/install-escpod.sh`, or `pixi run install-escpod`) into `resources/tools/escapepod/<version>/bin`. This is mandatory, not a convenience: the published escapepod-rs release binaries are default-features only and do **not** contain a working `escpod demux`. The build uses `--features cnn-detect`, which implies `demux` and adds `--method cnn`.
- **Authenticated `gh` CLI** — leech is installed from the GitHub release wheels of the private `rnabioco/leech` repo (`scripts/install-leech.sh`, or `pixi run install-leech`), pinned by `leech_version` / `leech_core_version` in `config/config-base.yml`. Set `LEECH_WHEEL_DIR` to a local directory of wheels to bypass `gh`. leech requires Python >= 3.12.

There is no longer a `resources/leech` git submodule (and no `.gitmodules`) — do not add `git submodule update` steps.

Optional per-component tasks: `pixi run install-escpod`, `pixi run install-leech`, `pixi run install-demux-models`, `pixi run install-warpdemux`.

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

GPU-intensive rules (rebasecall, classify_charging) automatically request GPU resources via queue/partition configuration.

## Architecture

### Snakemake Workflow Structure

The workflow is modular with rules split across multiple files:

```
workflow/
├── Snakefile                          # Main entry point, includes all rule modules
├── rules/
│   ├── common.smk                     # Sample parsing, helper functions, outputs definition
│   ├── aatrnaseq-process.smk          # Core processing: escpod merge → basecalling → alignment
│   ├── aatrnaseq-charging.smk         # Charging classification outputs (2 rules)
│   ├── aatrnaseq-qc.smk               # QC metrics: base calling errors, alignment stats (3 rules)
│   ├── aatrnaseq-modifications.smk    # Modification calling: coverage, modkit outputs (4 rules)
│   ├── demux.smk                      # Demux entry point (conditionally loaded): backend
│   │                                  #   dispatch + EDX splitting/concordance rules
│   ├── demux-escpod.smk               # escpod demux backend (default), included by demux.smk
│   └── demux-warpdemux.smk            # WarpDemuX backend, included by demux.smk
├── scripts/                           # Python scripts called by rules
│   ├── detect_3p_adapters.py           # Detect 3' adapter identity per read (for EDX splitting)
│   └── generate_squiggy_session.py    # Generate Squiggy/Positron session JSON
└── envs/
    └── aatrnaseqpipe-env.yml          # Conda environment (legacy)
```

**Key Architectural Details:**

- **Sample Management**: `workflow/rules/common.smk` contains `parse_samples()` which reads `config/samples.tsv` and `find_raw_inputs()` which recursively searches for pod5 files in specified directories
- **Tool Management**: Modkit is managed by pixi. Everything else is installed by `scripts/setup-tools.sh` (`pixi run setup`): dorado + models, remora, the `escpod` CLI (source build, `scripts/install-escpod.sh`), the `escapepod` python package, leech (`scripts/install-leech.sh`), the escpod demux models (`scripts/install-demux-models.sh`), and WarpDemuX. `scripts/setup-env.sh` runs on `pixi shell` activation and puts `resources/tools/dorado/<version>/bin` and `resources/tools/escapepod/<version>/bin` on `PATH`
- **POD5 I/O**: All POD5 manipulation goes through `escpod` (escapepod-rs), not the ONT `pod5` CLI — `escpod merge` in `merge_pods`, `escpod filter` in `split_pod5` and `filter_pod5_by_edx`. It is 3-9x faster on these operations and writes crash-safely (output is staged to a temp file and renamed, so an interrupted run never leaves a corrupt POD5). On the python side, `escapepod` (PyPI, pinned by `escapepod_version`) is a mostly-drop-in replacement for the ONT `pod5` python package and is what leech uses internally. The ONT `pod5` python package is still installed for exactly one optional path: `workflow/scripts/extract_signal_metrics.py` hands pod5 objects to remora's `io` API and so cannot use escapepod's reader; that script only runs when `remora_kmer_table` is set (default `null`)
- **`escpod filter` behavior differences** vs ONT `pod5 filter`, worth remembering when editing these rules: it takes a **single positional input** (a file, or a directory it walks recursively) — this is why `split_pod5` passes the run directory rather than a list of `pod5_pass`/`pod5_fail`/`pod5` dirs; there is no `--missing-ok` flag, since read IDs absent from the input are a warning rather than an error (the flag was dropped); and an **empty read-ID file is a hard error** ("No read IDs found") where ONT `pod5 --missing-ok` produced an empty POD5, so a barcode or EDX adapter matching zero reads now fails that job instead of silently emitting an empty file
- **Output Aggregation**: `pipeline_outputs()` in `common.smk` defines all final output files for the `rule all` target

### Pipeline Flow

```
POD5 files → merge_pods → rebasecall (Dorado) → ubam_to_fastq → bwa_align →
classify_charging (Remora) → transfer_bam_tags → add_adapter_tags → finalize_bam → Summary tables
```

For EDX samples (dual barcoding), 3' adapter detection and FASTQ/POD5 splitting happens before alignment:
```
rebasecall → detect_edx_adapters → extract_edx_read_ids
                                     ├── filter_fastq_by_edx → bwa_align → ...
                                     └── filter_pod5_by_edx → classify_charging
```

### Core Processing Pipeline (aatrnaseq-process.smk)

1. **merge_pods**: Merge all pod5 files per sample into single pod5 (`escpod merge`)
2. **rebasecall**: Use dorado to rebasecall with move tables (required for Remora)
3. **ubam_to_fastq**: Extract reads from unmapped BAM to FASTQ
4. **bwa_align**: Align reads to tRNA + adapter reference with BWA MEM
5. **classify_charging**: Use Remora model to classify charged vs uncharged reads (adds ML tag to BAM)
6. **transfer_bam_tags**: Transfer alignment tags back to classified BAM (ML→cl, MM→cm)
7. **add_adapter_tags**: Detect adapter positions and add pt tags with 5'/3' boundaries
8. **finalize_bam**: Symlink adapter-tagged BAM as final output (EDX filtering now happens before alignment)

### Summary Generation

After classification, generates (split across three rule files):

**aatrnaseq-charging.smk:**
- Charging probability tables (ML tag values per read)
- CPM (counts per million) for charged/uncharged tRNA

**aatrnaseq-qc.smk:**
- Base calling error frequencies
- Alignment statistics
- Remora signal metrics (if kmer table provided)

**aatrnaseq-modifications.smk:**
- Coverage bedGraph files (counts and CPM)
- Modkit modification pileups
- Modkit per-read modification calls

**common.smk:**
- Squiggy session JSON (`squiggy-session.json`) for loading outputs in Positron IDE

## Configuration

### Main Config Files

- `config/config-base.yml`: Base configuration included by Snakefile
  - Base calling model path
  - Reference fasta
  - Remora models and kmer tables
  - Dorado version for download
  - Pinned tool versions: `escapepod_version` (escpod CLI source build + `escapepod` PyPI package), `leech_version`, `leech_core_version`
  - Demux backend settings: `warpdemux.backend`, `warpdemux.method`, `warpdemux.barcode_model`, `warpdemux.adapter_model`
  - Command-line options for tools (dorado, bwa, filters)

- `config/samples.tsv`: Two-column TSV (no header)
  - Column 1: Unique sample ID
  - Column 2: Path to sequencing run folder containing pod5_pass/pod5_fail/pod5 subdirectories

- `config/config-test.yml`: Overrides base config for test data

### Important Config Parameters

- **opts.bam_filter**: Controls full-length read filtering (`-5 24 -3 23 -s` requires 24bp 5' adapter, 23bp 3' adapter, positive strand)
- **opts.dorado**: Includes `--modified-bases m5C_2OmeC inosine_m6A_2OmeA pseU_2OmeU 2OmeG --emit-moves` for modification calling and move tables
- **opts.bwa**: RNA-optimized alignment parameters (`-W 13 -k 6 -T 20 -x ont2d`)
- **ml-threshold**: Currently hardcoded in `get_cca_trna_cpm` rule (200-255 = charged, <200 = uncharged)
- **cleanup_intermediates**: Opt-in auto-deletion of large regenerable intermediates during a run, via `temp()`. Accepts a bool or a list of tier names (`cascade`, `basecall`, `fastq`, `merged_pod5`, `demux_scratch`, `split_pod5`) resolved by `maybe_temp()` / `_enabled_cleanup_tiers()` in `common.smk`. `demux_scratch` also covers the escpod `classifications.csv`; `split_pod5` also covers `demux/escpod_output`. `bam/final` and `demux/edx/pod5` (the classification-input POD5) are always kept, as is `demux/read_ids/{run_id}/demux_summary.tsv.gz`. Only enable `split_pod5` for all-EDX runs (where `demux/edx/pod5` is the leaf classification input); for non-EDX/mixed runs `demux/pod5` must be kept. The on-demand `clean` rule (`rules/clean.smk`) remains the catch-all superset for reclaiming space on already-completed runs. See `config/README.md` for tier→directory mapping.

## Demultiplexing (Optional)

The pipeline supports optional signal-level barcode demultiplexing for pooled/multiplexed sequencing runs. The config section is still named `warpdemux` for backward compatibility, but there are now two interchangeable backends selected by `warpdemux.backend`.

### Enabling Demultiplexing

1. **Install tools/models**: `pixi run setup` (builds `escpod` with demux, installs the demux models, and installs WarpDemuX)
2. **Create YAML sample file** with barcode assignments (see `config/samples-demux-example.yml`)
3. **Enable in config**: Set `warpdemux.enabled: true`

### Backends (`warpdemux.backend`)

Both backends produce the same downstream contract — `demux/read_ids/{run_id}/barcode_mapping.tsv.gz` and `demux/pod5/{sample}/{sample}.pod5` — so everything downstream is backend-agnostic.

- **`escpod` (default, `demux-escpod.smk`)**: one fused `escpod demux` pass (detect → fingerprint → classify → split) that emits per-read classifications *and* one POD5 per barcode. It therefore does not need the separate read-ID extraction + POD5 filter pass, and runs strictly fewer jobs. Rules: `escpod_demux`, `parse_escpod_demux`, `collect_escpod_pod5` (symlinks the per-barcode POD5 into the per-sample path).
- **`warpdemux` (`demux-warpdemux.smk`)**: the original python implementation, still fully supported. It only classifies, so POD5 splitting is a second `escpod filter` pass. Rules: `warpdemux`, `parse_warpdemux`, `extract_sample_reads`, `split_pod5`.

**The two backends do NOT produce identical barcode calls.** `escpod demux` cannot load a WarpDemuX kit: the tRNA kits (`WDX4_tRNA_rna004_v1_0`, `WDX4b_tRNA_rna004_v1_0`) are `warpdemux.models.fpt_boost.Fpt_Boost` objects wrapping a CatBoost classifier, and every converter shipped with escapepod-rs requires either a scikit-learn SVC or a scikit-learn `HistGradientBoostingClassifier`. Conversion is not possible. Instead escpod uses `barcode_wdx4_rna004`, a GBM *distilled* from the `WDX4_tRNA_rna004_v1_0` teacher's high-confidence calls (confidence >= 0.9), trained and served on the same Rust `--warpdemux-compat` fingerprint. It scores 0.971 balanced recall against that teacher (0.988 at 95% recovery) and its label mapping (`{0:3, 1:4, 2:5, 3:7}`) is identical to the tRNA kit's, so barcode *numbering* is unchanged — escpod labels them `BC03` where the pipeline uses `barcode03`, and `wdx_to_escpod_barcode()` / `parse_escpod_demux` translate.

Two semantic deltas to keep in mind:

1. The shipped GBM carries **no per-class confidence thresholds**, so escpod assigns every read with a usable adapter boundary to some barcode. The `unclassified` fraction drops sharply versus WarpDemuX, which rejected low-confidence reads. Filter on the `confidence` column of the classifications table to restore WarpDemuX-like rejection.
2. `bc05` <-> `bc07` is the residual confusable pair.

The 0.971 figure is agreement with the WarpDemuX teacher on held-out *confident* calls, not independent ground truth; read-level concordance on real aa-tRNA-seq data has not been measured. To validate the swap, run both backends on one run and compare their `demux_summary.tsv.gz` / `barcode_mapping.tsv.gz` outputs (`escapepod-rs` ships `scripts/compare_demux_results.py` for a per-read concordance table).

`warpdemux.method` selects boundary detection for the escpod backend: `cnn` (default, matches how the barcode model was trained) or `llr` (needs no ONNX model, but agrees with the CNN only ~82% within +/-200 samples, which shifts calls).

### Demux Models

`bash scripts/install-demux-models.sh` (or `pixi run install-demux-models`) installs the models into `resources/models/demux/` (gitignored), verifying sha256 against the upstream `MANIFEST.json`. It sources them from a local `rnabioco/escapepod-models` checkout (`ESCAPEPOD_MODELS_DIR`, or a sibling directory) or, failing that, a GitHub release. **Caveat:** the barcode GBM models are not yet published as releases on `rnabioco/escapepod-models` (only `adapter_rna004@v1.0.1` is), so today the local-checkout path is required; the script prints the `scripts/release_model.sh` command needed to publish them. Config keys: `warpdemux.barcode_model`, `warpdemux.adapter_model`.

### Demux Outputs

`demux/read_ids/{run_id}/demux_summary.tsv.gz` (per-barcode read counts) is a first-class pipeline output whenever demux is enabled — it used to be cleanup-tiered scratch. With the escpod backend it also carries a `mean_confidence` column.

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

EDX values must directly match adapter names from `adapters.three_prime` in the config (e.g., `edx01`, `edx02`). EDX splitting happens before alignment — the pipeline detects 3' adapter identity on the uBAM, then filters FASTQ and POD5 so downstream rules only process matching reads.

When `edx.enabled: true` in config, the `edx_concordance` rule produces `summary/edx/edx_concordance.tsv.gz` — a concordance table of WDX assignment vs EDX (3' adapter) identity per sample.

### Running with Demux

```bash
# Dry run with demux config
pixi run snakemake -n --configfile=config/config-demux-test.yml

# Execute with demux
pixi run snakemake --configfile=config/config-demux-test.yml --cores 8
```

## Charged vs Uncharged Classification

The pipeline uses Remora machine learning to classify charging state:

- **Model Location**: `remora_cca_classifier` config parameter (resources/models/cca_classifier.pt)
- **Signal Region**: 6-nucleotide kmer spanning CCA 3' end + first 3 adapter bases (CCAGGC)
- **ML Tag**: Classification score stored in BAM ML tag (0-255 scale)
- **Threshold**: ML ≥ 200 = charged, ML < 200 = uncharged (adjustable in get_cca_trna_cpm rule)
- **Filtering**: Only full-length tRNA reads with proper 5'/3' adapters are classified

### leech (alternative classifier)

Setting `classifier: leech` routes charging classification through the `classify_charging_leech` rule (GPU-accelerated), and `classify_aa.enabled` / `aa_identity.enabled` use leech for amino-acid identity. Points to remember:

- leech writes the `ML` (score) and `MP` (query position) tags on the binary charging path — note **`MP` where remora writes `MM`**, so downstream tag transfer differs by classifier.
- leech is at v0.4.1 (was ~v0.2.0), pinned by `leech_version`; `leech-core` (`leech_core_version`) is a separate optional Rust extension wheel. Without `leech-core`, leech falls back to a slower pure-python extraction backend and `--backend rust` is unavailable.
- v0.4.x CLI changes that affected the pipeline: `--reference-anchored` is deprecated and hidden in favor of `--anchor {basecall,reference}`, and `--anchor` now **defaults to `reference`**; `classify_charging_leech` passes `--anchor reference` explicitly. `--raw` also changed meaning: `pn` tags are always written now, and `ac`/`am`/`pp` default to compact uint8 encoding unless `--raw` requests floats.
- leech requires Python >= 3.12.

## Development

### Adding New Rules

When adding new Snakemake rules:
- Place processing rules in `aatrnaseq-process.smk`
- Place charging analysis rules in `aatrnaseq-charging.smk`
- Place QC/statistics rules in `aatrnaseq-qc.smk`
- Place modification/coverage rules in `aatrnaseq-modifications.smk`
- Place backend-agnostic demux/EDX rules in `demux.smk`; backend-specific rules go in `demux-escpod.smk` or `demux-warpdemux.smk` (a new backend rule must still produce `demux/read_ids/{run_id}/barcode_mapping.tsv.gz` and `demux/pod5/{sample}/{sample}.pod5`)
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
```

### Cluster Resource Configuration

**LSF** - Modify `cluster/lsf/config.yaml` to adjust:
- Memory requirements per rule (mem_mb)
- GPU queue assignments (lsf_queue)
- LSF project tags (lsf_project)
- Maximum concurrent jobs

Rules requiring GPU (rebasecall, classify_charging) must set:
- lsf_queue: "gpu"
- lsf_extra: "-gpu num=1:j_exclusive=yes"
- ngpu: 1

**SLURM** - Modify `cluster/slurm/config.yaml` to adjust:
- Memory requirements per rule (mem_mb)
- GPU partition (slurm_partition)
- Account/allocation (slurm_account)
- Runtime limits (runtime, in minutes)
- Maximum concurrent jobs

Rules requiring GPU (rebasecall, classify_charging) must set:
- slurm_partition: "gpu" (or your cluster's GPU partition)
- gres: "gpu:1"

## Important Notes

- The pipeline requires Snakemake 8.0+
- Modkit is managed by pixi; remora, escpod/escapepod, leech, dorado and WarpDemuX come from `pixi run setup`
- Building `escpod` requires a Rust toolchain (>= 1.95); installing leech requires an authenticated `gh` (or `LEECH_WHEEL_DIR`) and Python >= 3.12
- There is no `resources/leech` submodule any more — leech comes from release wheels
- `Snakefile`'s `onstart` hook validates tool locations and, when the escpod demux backend is selected, warns if `escpod demux` is unavailable (i.e. a default-features build)
- The pipeline tracks git commit ID for reproducibility (see `get_pipeline_commit()`)
- CUDA_VISIBLE_DEVICES is passed through to dorado if set
- Pod5 files are searched recursively in pod5_pass/pod5_fail/pod5 subdirectories
- The ML threshold for charging classification is currently hardcoded in the `get_cca_trna_cpm` rule

## Key Outputs

Outputs go to directory specified by `output_dir` in config. Test outputs: `.tests/outputs/`

Key outputs per sample:
- `summary/tables/{sample}/{sample}.charging.cpm.tsv.gz` - CPM-normalized charging counts
- `summary/tables/{sample}/{sample}.charging_prob.tsv.gz` - Per-read charging probabilities
- `bam/final/{sample}/{sample}.bam` - Final BAM with cl/cm (charging) and pt (adapter positions) tags

Pipeline-level outputs:
- `squiggy-session.json` - Squiggy session file for loading samples in Positron
- `demux/read_ids/{run_id}/demux_summary.tsv.gz` - per-barcode read counts (when demux is enabled; includes `mean_confidence` with the escpod backend)
