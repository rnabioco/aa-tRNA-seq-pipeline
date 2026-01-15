# aa-tRNA-seq-pipeline

[![CI](https://github.com/rnabioco/aa-tRNA-seq-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/rnabioco/aa-tRNA-seq-pipeline/actions/workflows/ci.yml)
[![Lint](https://github.com/rnabioco/aa-tRNA-seq-pipeline/actions/workflows/lint.yml/badge.svg)](https://github.com/rnabioco/aa-tRNA-seq-pipeline/actions/workflows/lint.yml)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://rnabioco.github.io/aa-tRNA-seq-pipeline)

A Snakemake pipeline to process ONT aa-tRNA-seq data.

## Prerequisites

This pipeline uses [Pixi](https://pixi.sh) for environment management. Install it first:

```bash
curl -fsSL https://pixi.sh/install.sh | sh
```

See the [Pixi installation guide](https://pixi.prefix.dev/latest/installation/) for alternative methods (Homebrew, Windows, etc.).

## Usage

The pipeline can be configured by editing the `config/config.yml` file. The config file specifications will
run a small example dataset through the pipeline.

```bash
git clone https://github.com/rnabioco/aa-tRNA-seq-pipeline.git
cd aa-tRNA-seq-pipeline

# Install environment
pixi install

# Download test data and setup tools (first time only)
pixi run dl-test-data
pixi run setup-tools

# Dry run
pixi run dry-run

# Run pipeline locally with test data
pixi run test
```

## Configuration

To use on your own samples, create a config file and sample file in `config/`.

### Standard Usage (Non-Multiplexed Runs)

Create a TSV sample file with sample IDs and run paths:

```
sample1    /path/to/run1
sample2    /path/to/run2
```

Then create a config file pointing to it:

```yaml
samples: config/samples.tsv
output_directory: "results"
```

### Multiplexed Runs (WarpDemuX Demultiplexing)

For pooled sequencing runs with WarpDemuX barcodes, use a YAML sample file:

```yaml
runs:
  - path: /path/to/pooled/run
    barcode_kit: "WDX4_tRNA_rna004_v1_0"
    samples:
      charged_sample: "barcode03"
      uncharged_sample: "barcode04"
```

Enable demultiplexing in your config:

```yaml
samples: config/samples-demux.yml
output_directory: "results"

warpdemux:
    enabled: true
    barcode_kit: "WDX4_tRNA_rna004_v1_0"
```

Run with the demux environment:

```bash
pixi run -e demux snakemake --configfile=config/config-demux.yml --cores 8
```

See [README.md in the config directory](https://github.com/rnabioco/aa-tRNA-seq-pipeline/tree/main/config) for additional details on all configuration options.

## Workflow

```mermaid
flowchart TD
    subgraph Input
        POD5[POD5 files]
    end

    subgraph Demux [Optional Demultiplexing]
        W[warpdemux<br/>barcode classification]
    end

    subgraph Processing
        A[merge_pods] --> B[rebasecall<br/>Dorado + move tables]
        B --> C[ubam_to_fastq]
        C --> D[bwa_align<br/>tRNA + adapter reference]
        D --> E[filter_reads<br/>full-length tRNAs only]
    end

    subgraph Classification
        E --> F[classify_charging<br/>Remora ML model]
        B -.-> F
        A -.-> F
        F --> G[transfer_bam_tags]
    end

    subgraph Outputs
        G --> H[charging_prob<br/>per-read ML scores]
        G --> I[get_cca_trna_cpm<br/>CPM counts]
        G --> J[bcerror<br/>basecalling errors]
        G --> K[align_stats]
        G --> L[modkit pileups]
    end

    POD5 -.-> W
    W -.-> A
    POD5 --> A
```

Given a directory of POD5 files, this pipeline:

1. **(Optional) Demultiplexes** pooled runs using WarpDemuX barcode classification
2. **Merges** all POD5 files per sample into a single file
3. **Rebasecalls** with Dorado to generate unmapped BAM with move tables (required for Remora)
4. **Converts** BAM to FASTQ and **aligns** to tRNA + adapter reference with BWA MEM
5. **Filters** for full-length tRNA reads with proper adapter boundaries
6. **Classifies** charged vs. uncharged reads using a Remora model trained on nanopore signal over the CCA 3' end

The classification generates ML tag values (0-255) indicating the likelihood of aminoacylation. By default, ML values of 200-255 are treated as charged, and values <200 as uncharged. This threshold can be adjusted via the `ml-threshold` parameter in the `get_cca_trna_cpm` rule.

The final steps of the pipeline calculate a number of outputs that may be useful for analysis and visualization, including normalized counts for charged and uncharged tRNA (`get_cca_trna_cpm`), basecalling error values (`bcerror`), alignment statistics (`align_stats`) and information on raw nanopore signal from Remora (`remora_signal_stats`).

### Remora classification

A few notes about Remora classification for charged vs. uncharged tRNA reads

1. this step retains only full length tRNA reads (with an allowance for signal loss at the 5´ end of nanopore direct RNA sequencing)
2. Additionally, due to the iterative nature of sequencing method development, the present approach does not rely on differences in adapter sequences attached to charged vs. uncharged tRNA molecules (though these sequences are retained as separate entries in the alignment reference and downstream files). While we anticipate being able to leverage this information in the future, the current pipeline relies exclusively on signal data over a 6-nt modification kmer spanning the universal CCA 3′ end of tRNA and the first three nucleotides of the 3′ adapter (CCAGGC) to distinguish charged and uncharged reads.

## Cluster execution

The pipeline includes cluster profiles for LSF and SLURM schedulers.

```bash
# Run test data on LSF cluster
pixi run test-lsf

# Run full preprint analysis on cluster
pixi run run-preprint
```

For more details on configuring HPC jobs, see `cluster/lsf/config.yaml` or `cluster/generic/config.yaml`.

## Citation

To cite this work, see:

> White LK, Radakovic A, Sajek MP, Dobson K, Riemondy KA, Del Pozo S, Szostak JW, Hesselberth JR. Nanopore sequencing of intact aminoacylated tRNAs. Nat Commun. 2025 Aug 20;16(1):7781. doi: 10.1038/s41467-025-62545-9. PMID: 40835813; PMCID: PMC12368100.

