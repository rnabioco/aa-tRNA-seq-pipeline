# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Snakemake pipeline for processing Oxford Nanopore Technologies (ONT) aa-tRNA-seq data. The pipeline distinguishes between charged (aminoacylated) and uncharged tRNA molecules using Remora machine learning models trained on nanopore signal data over the CCA 3' end of tRNA molecules.

## Setup and Environment

### Initial Setup

```bash
# Create conda environment
mamba env create -f workflow/envs/aatrnaseqpipe-env.yml
mamba activate aatrnaseqpipe

# Download test data (first time only)
bash .tests/dl_test_data.sh

# Install dorado and modkit (first time only)
snakemake setup_dorado dorado_model setup_modkit
```

### Running the Pipeline

```bash
# Dry run with test config
snakemake -n --configfile=config/config-test.yml

# Execute locally (specify cores)
snakemake --cores 12 --configfile=config/config-test.yml

# Run on LSF cluster (uses cluster/lsf profile)
bsub < run-test.sh
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
│   ├── tool_setup.smk                 # Dorado and modkit installation
│   ├── aatrnaseq-process.smk          # Core processing: pod5 merge → basecalling → alignment
│   └── aatrnaseq-summaries.smk        # Summary statistics and output tables
├── scripts/                           # Python scripts called by rules
└── envs/
    └── aatrnaseqpipe-env.yml          # Conda environment
```

**Key Architectural Details:**

- **Sample Management**: `workflow/rules/common.smk` contains `parse_samples()` which reads `config/samples.tsv` and `find_raw_inputs()` which recursively searches for pod5 files in specified directories
- **Dynamic PATH Setup**: The Snakefile `onstart` handler dynamically adds dorado and modkit binaries to PATH based on configured versions
- **Output Aggregation**: `pipeline_outputs()` in `common.smk` defines all final output files for the `rule all` target

### Core Processing Pipeline (aatrnaseq-process.smk)

The main data flow through the pipeline:

1. **merge_pods**: Merge all pod5 files per sample into single pod5
2. **rebasecall**: Use dorado to rebasecall with move tables (required for Remora)
3. **ubam_to_fastq**: Extract reads from unmapped BAM to FASTQ
4. **bwa_align**: Align reads to tRNA + adapter reference with BWA MEM
5. **filter_reads**: Filter for full-length tRNA reads with proper adapter boundaries
6. **classify_charging**: Use Remora model to classify charged vs uncharged reads (adds ML tag to BAM)
7. **transfer_bam_tags**: Transfer alignment tags back to classified BAM

### Summary Generation (aatrnaseq-summaries.smk)

After classification, generates:
- Charging probability tables (ML tag values per read)
- CPM (counts per million) for charged/uncharged tRNA
- Base calling error frequencies
- Alignment statistics
- Remora signal metrics (if kmer table provided)
- Modkit modification calls and pileups

## Configuration

### Main Config Files

- `config/config-base.yml`: Base configuration included by Snakefile
  - Base calling model path
  - Reference fasta
  - Remora models and kmer tables
  - Tool versions (dorado, modkit)
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
- Place summary/analysis rules in `aatrnaseq-summaries.smk`
- Add helper functions to `common.smk`
- Reference Python scripts should go in `workflow/scripts/`
- Update `pipeline_outputs()` if rule produces final outputs

### Testing Changes

```bash
# Always test with dry run first
snakemake -n --configfile=config/config-test.yml

# Run specific rule
snakemake <rule_name> --configfile=config/config-test.yml

# Force rerun of specific rule
snakemake <rule_name> --forcerun <rule_name> --configfile=config/config-test.yml
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
- Dorado and modkit are installed by the pipeline (not via conda) to specific versions
- The pipeline tracks git commit ID for reproducibility (see `get_pipeline_commit()`)
- CUDA_VISIBLE_DEVICES is passed through to dorado if set
- Pod5 files are searched recursively in pod5_pass/pod5_fail/pod5 subdirectories
- The ML threshold for charging classification is currently hardcoded in the `get_cca_trna_cpm` rule
