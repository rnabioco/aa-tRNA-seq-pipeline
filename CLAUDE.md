# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Snakemake pipeline for processing Oxford Nanopore Technologies (ONT) aa-tRNA-seq data. The pipeline distinguishes between charged (aminoacylated) and uncharged tRNA molecules using Remora machine learning models trained on nanopore signal data over the CCA 3' end of tRNA molecules.

## Setup and Environment

### Initial Setup

```bash
# Install all dependencies (modkit, remora, and other tools)
pixi install

# Enter the environment (downloads dorado on first activation)
pixi shell

# Download test data (first time only)
pixi run dl-test-data
```

### Running the Pipeline

```bash
# Dry run with test config
pixi run dry-run

# Run locally with test data (4 cores)
pixi run test

# Run on LSF cluster
pixi run test-lsf

# Run preprint pipeline on cluster
pixi run run-preprint
```

### Direct Snakemake Commands

```bash
pixi run snakemake --configfile=config/config-test.yml --cores 8
```

### Cluster Execution

The pipeline is optimized for LSF scheduler. Key files:
- `run-test.sh`: Test data execution on LSF
- `run-preprint.sh`: Full preprint data execution on LSF
- `cluster/lsf/config.yaml`: LSF-specific resource configurations

GPU-intensive rules (rebasecall, classify_charging) automatically request GPU resources via LSF queue configuration.

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
│   └── warpdemux.smk                  # WarpDemuX demultiplexing (conditionally loaded)
├── scripts/                           # Python scripts called by rules
└── envs/
    └── aatrnaseqpipe-env.yml          # Conda environment (legacy)
```

**Key Architectural Details:**

- **Sample Management**: `workflow/rules/common.smk` contains `parse_samples()` which reads `config/samples.tsv` and `find_raw_inputs()` which recursively searches for pod5 files in specified directories
- **Tool Management**: Modkit and Remora are managed by pixi. Dorado is downloaded on first `pixi shell` activation via `scripts/setup-dorado.sh`
- **Output Aggregation**: `pipeline_outputs()` in `common.smk` defines all final output files for the `rule all` target

### Pipeline Flow

```
POD5 files → merge_pods → rebasecall (Dorado) → ubam_to_fastq → bwa_align →
classify_charging (Remora) → transfer_bam_tags → Summary tables
```

### Core Processing Pipeline (aatrnaseq-process.smk)

1. **merge_pods**: Merge all pod5 files per sample into single pod5
2. **rebasecall**: Use dorado to rebasecall with move tables (required for Remora)
3. **ubam_to_fastq**: Extract reads from unmapped BAM to FASTQ
4. **bwa_align**: Align reads to tRNA + adapter reference with BWA MEM
5. **classify_charging**: Use Remora model to classify charged vs uncharged reads (adds ML tag to BAM)
6. **transfer_bam_tags**: Transfer alignment tags back to classified BAM

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

## Configuration

### Main Config Files

- `config/config-base.yml`: Base configuration included by Snakefile
  - Base calling model path
  - Reference fasta
  - Remora models and kmer tables
  - Dorado version for download
  - Command-line options for tools (dorado, bwa, filters)

- `config/samples.tsv`: Two-column TSV (no header)
  - Column 1: Unique sample ID
  - Column 2: Path to sequencing run folder containing pod5_pass/pod5_fail/pod5 subdirectories

- `config/config-test.yml`: Overrides base config for test data

### Important Config Parameters

- **opts.bam_filter**: Controls full-length read filtering (`-5 24 -3 23 -s` requires 24bp 5' adapter, 23bp 3' adapter, positive strand)
- **opts.dorado**: Includes `--modified-bases pseU m5C inosine_m6A --emit-moves` for modification calling and move tables
- **opts.bwa**: RNA-optimized alignment parameters (`-W 13 -k 6 -T 20 -x ont2d`)
- **ml-threshold**: Currently hardcoded in `get_cca_trna_cpm` rule (200-255 = charged, <200 = uncharged)

## WarpDemuX Demultiplexing (Optional)

The pipeline supports optional barcode demultiplexing using WarpDemuX for pooled/multiplexed sequencing runs.

### Enabling Demultiplexing

1. **Install demux environment**: `pixi install -e demux && pixi run -e demux install-warpdemux`
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

### Running with Demux

```bash
# Dry run with demux config
pixi run -e demux snakemake -n --configfile=config/config-demux-test.yml

# Execute with demux
pixi run -e demux snakemake --configfile=config/config-demux-test.yml --cores 8
```

## Charged vs Uncharged Classification

The pipeline uses Remora machine learning to classify charging state:

- **Model Location**: `remora_cca_classifier` config parameter (resources/models/cca_classifier.pt)
- **Signal Region**: 6-nucleotide kmer spanning CCA 3' end + first 3 adapter bases (CCAGGC)
- **ML Tag**: Classification score stored in BAM ML tag (0-255 scale)
- **Threshold**: ML ≥ 200 = charged, ML < 200 = uncharged (adjustable in get_cca_trna_cpm rule)
- **Filtering**: Only full-length tRNA reads with proper 5'/3' adapters are classified

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
```

### Cluster Resource Configuration

Modify `cluster/lsf/config.yaml` to adjust:
- Memory requirements per rule (mem_mb)
- GPU queue assignments
- LSF project tags
- Maximum concurrent jobs

Rules requiring GPU (rebasecall, classify_charging) must set:
- lsf_queue: "gpu"
- lsf_extra: "-gpu num=1:j_exclusive=yes"
- ngpu: 1

## Important Notes

- The pipeline requires Snakemake 8.0+
- Modkit and Remora are managed by pixi (bioconda and pypi-dependencies)
- Dorado is downloaded automatically on first `pixi shell` activation
- The pipeline tracks git commit ID for reproducibility (see `get_pipeline_commit()`)
- CUDA_VISIBLE_DEVICES is passed through to dorado if set
- Pod5 files are searched recursively in pod5_pass/pod5_fail/pod5 subdirectories
- The ML threshold for charging classification is currently hardcoded in the `get_cca_trna_cpm` rule

## Key Outputs

Outputs go to directory specified by `output_dir` in config. Test outputs: `.tests/outputs/`

Key outputs per sample:
- `summary/tables/{sample}/{sample}.charging.cpm.tsv.gz` - CPM-normalized charging counts
- `summary/tables/{sample}/{sample}.charging_prob.tsv.gz` - Per-read charging probabilities
- `bam/final/{sample}.bam` - Final BAM with CL/CM charging tags
