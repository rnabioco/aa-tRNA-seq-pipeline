# Workflow Overview

This page provides a high-level view of the aa-tRNA-seq pipeline architecture and data flow.

## Pipeline Architecture

The pipeline is organized into modular Snakemake rule files:

```mermaid
flowchart TB
    subgraph rules[Rule Files]
        A[aatrnaseq-process.smk<br/>Core processing]
        B[aatrnaseq-charging.smk<br/>Charging analysis]
        C[aatrnaseq-qc.smk<br/>Quality control]
        D[aatrnaseq-modifications.smk<br/>Modification calling]
        OR[aatrnaseq-odds-ratios.smk<br/>Odds ratio analysis]
        R[aatrnaseq-report.smk<br/>QC report]
        E[warpdemux.smk<br/>Demultiplexing<br/><i>conditional</i>]
    end

    subgraph common[common.smk]
        F[Sample parsing]
        G[Helper functions]
        H[Output definitions]
    end

    I[Snakefile<br/>Main entry] --> rules
    I --> common
```

## Complete Pipeline Flow

### Standard Pipeline (No Demultiplexing)

```mermaid
flowchart TB
    subgraph Input
        A[POD5 files<br/>per sample]
    end

    subgraph Processing[aatrnaseq-process.smk]
        B[merge_pods<br/>Merge POD5s]
        C[rebasecall<br/>Dorado basecalling]
        D[ubam_to_fastq<br/>Extract FASTQ]
        E[bwa_align<br/>Align to reference]
        F[classify_charging<br/>Remora ML]
        G[transfer_bam_tags<br/>Rename tags]
        G2[add_adapter_tags<br/>PT tags]
        G3[finalize_bam<br/>Symlink final BAM]
    end

    subgraph Charging[aatrnaseq-charging.smk]
        H[get_cca_trna<br/>Extract probabilities]
        I[get_cca_trna_cpm<br/>Calculate CPM]
    end

    subgraph QC[aatrnaseq-qc.smk]
        J[base_calling_error<br/>Error metrics]
        K[align_stats<br/>Read statistics]
        L[remora_signal_stats<br/>Signal metrics]
    end

    subgraph Mods[aatrnaseq-modifications.smk]
        M[bam_to_coverage<br/>Coverage tracks]
        N[modkit_pileup<br/>Modification consensus]
        O[modkit_extract_calls<br/>Per-read mods]
        P[modkit_extract_full<br/>Full export]
    end

    subgraph OddsRatios[aatrnaseq-odds-ratios.smk]
        Q[compute_odds_ratios<br/>Pairwise OR]
    end

    subgraph Report[aatrnaseq-report.smk]
        R[render_combined_qc_report<br/>QC report]
    end

    A --> B --> C --> D --> E --> F --> G --> G2 --> G3

    G3 --> H --> I
    G3 --> J
    G3 --> K
    G3 --> L
    G3 --> M
    G3 --> N
    G3 --> O
    G3 --> P
    O --> Q
    H --> Q
    J --> R
    H --> R
    I --> R
    K --> R
```

### With Demultiplexing (WarpDemuX + EDX)

```mermaid
flowchart TB
    subgraph Input
        A[Pooled POD5 files<br/>per run]
    end

    subgraph WDX[WDX Demultiplexing]
        B[warpdemux<br/>Barcode prediction]
        C[parse_warpdemux<br/>Create mapping]
        D[extract_sample_reads<br/>Per-sample IDs]
        E[split_pod5<br/>Split by WDX barcode]
    end

    subgraph Standard[Standard Pipeline]
        F[rebasecall]
    end

    subgraph EDX[EDX Early Splitting]
        G[detect_edx_adapters<br/>3' adapter ID per read]
        H[extract_edx_read_ids]
        I[filter_fastq_by_edx]
        J[filter_pod5_by_edx]
    end

    subgraph Downstream[Downstream Processing]
        K[bwa_align → classify_charging → ...]
    end

    A --> B --> C --> D --> E --> F --> G --> H
    H --> I --> K
    H --> J --> K
```

## Rule Categories

### Processing Rules

Core data processing from raw signal to classified reads:

| Rule | Purpose | GPU |
|------|---------|-----|
| `merge_pods` | Combine POD5 files per sample | No |
| `rebasecall` | Basecall with Dorado | Yes |
| `ubam_to_fastq` | Extract reads for alignment | No |
| `bwa_idx` | Build BWA index | No |
| `bwa_align` | Align reads to reference | No |
| `classify_charging` | ML charging classification | No |
| `transfer_bam_tags` | Rename ML→CL tags | No |
| `add_adapter_tags` | Add PT tags for adapter positions | No |
| `finalize_bam` | Symlink final BAM | No |

### Charging Analysis Rules

Extract and summarize charging classification:

| Rule | Purpose |
|------|---------|
| `get_cca_trna` | Extract per-read charging scores |
| `get_cca_trna_cpm` | Calculate CPM-normalized counts |

### Quality Control Rules

Generate QC metrics and statistics:

| Rule | Purpose |
|------|---------|
| `compute_reference_similarity` | Pairwise reference sequence similarity matrix |
| `base_calling_error` | Per-position error frequencies |
| `align_stats` | Read counts through pipeline |
| `remora_signal_stats` | Raw signal metrics |

### Modification Rules

RNA modification calling with Modkit:

| Rule | Purpose |
|------|---------|
| `bam_to_coverage` | Generate coverage tracks |
| `modkit_pileup` | Per-site modification consensus |
| `modkit_extract_calls` | Per-read modification calls |
| `modkit_extract_full` | Comprehensive modification export |

### Odds Ratio Rules

Per-tRNA pairwise modification odds ratios:

| Rule | Purpose |
|------|---------|
| `compute_odds_ratios` | Pairwise modification odds ratios per tRNA |

### Report Rules

QC report generation:

| Rule | Purpose |
|------|---------|
| `render_combined_qc_report` | Combined Quarto QC report with per-sample tabs |

### Demultiplexing Rules

Optional WarpDemuX barcode demultiplexing:

| Rule | Purpose |
|------|---------|
| `warpdemux` | Run WDX barcode prediction |
| `parse_warpdemux` | Parse predictions to mapping |
| `extract_sample_reads` | Filter reads by WDX barcode |
| `split_pod5` | Create per-sample WDX POD5s |
| `detect_edx_adapters` | Detect 3' adapter identity per read |
| `extract_edx_read_ids` | Extract matching read IDs for EDX |
| `filter_fastq_by_edx` | Create EDX-filtered FASTQ |
| `filter_pod5_by_edx` | Create EDX-filtered POD5 |
| `edx_concordance` | WDX vs EDX concordance table |

## Key Processing Steps

### 1. POD5 Merging

Individual POD5 files from a sequencing run are merged into a single file per sample:

```
run1/pod5_pass/*.pod5  ─┐
run1/pod5_fail/*.pod5  ─┼──► sample.pod5
run2/pod5/*.pod5       ─┘
```

### 2. Basecalling

Dorado re-basecalls with:

- Move tables (`--emit-moves`) required for Remora
- Modification calling (`--modified-bases pseU m5C inosine_m6A`)
- High-accuracy model (rna004_130bps_sup@v5.3.0)

### 3. Alignment

BWA MEM with RNA-optimized parameters:

- `-x ont2d` preset for ONT reads
- `-F 20`: Unmapped and reverse-strand read removal

### 4. Charging Classification

Remora ML model analyzes signal at CCA 3' end:

- Input: POD5 (signal) + BAM (alignment)
- Output: ML tag (0-255 score)
- Threshold: ≥200 = charged

!!! info "Remora Classification Details"

    A few notes about Remora classification for charged vs. uncharged tRNA reads:

    1. This step retains only full-length tRNA reads (with an allowance for signal loss at the 5' end of nanopore direct RNA sequencing)

    2. The current approach does not rely on differences in adapter sequences attached to charged vs. uncharged tRNA molecules (though these sequences are retained as separate entries in the alignment reference). The pipeline relies exclusively on signal data over a **6-nucleotide modification kmer** spanning the universal CCA 3' end of tRNA and the first three nucleotides of the 3' adapter (**CCAGGC**) to distinguish charged and uncharged reads.

### 5. Tag Transfer

Original Remora tags are renamed to avoid conflicts:

- `ML` → `CL` (charging likelihood)
- `MM` → `CM` (charging metadata)

### 6. Adapter Position Tagging

The `add_adapter_tags` rule adds PT tags with adapter boundaries:

- Uses parasail Smith-Waterman alignment
- Detects 5' and 3' adapter positions
- Can infer 5' adapter from alignment position when truncated

## Resource Requirements

### GPU Rules

These rules require GPU access:

| Rule | Typical Runtime | Memory |
|------|-----------------|--------|
| `rebasecall` | 30-60 min/sample | 24 GB |

### CPU-Intensive Rules

| Rule | Threads | Memory |
|------|---------|--------|
| `merge_pods` | 12 | 16 GB |
| `bwa_align` | 12 | 24 GB |
| `classify_charging` | 8 | 24 GB |
| `modkit_extract_full` | 12 | 48 GB |

### Memory-Intensive Rules

| Rule | Memory |
|------|--------|
| `modkit_extract_calls` | 96 GB |
| `warpdemux` | 32 GB |
| `remora_signal_stats` | 24 GB |

## Configuration Points

Key parameters that affect pipeline behavior:

| Parameter | Affects |
|-----------|---------|
| `opts.dorado` | Basecalling modifications |
| `opts.bwa` | Alignment sensitivity |
| `opts.bam_filter` | Full-length read filtering |
| `modkit.mod_thresholds` | Modification calling stringency |
| `warpdemux.barcode_kit` | Demultiplexing model |

## Next Steps

- [Rules Reference](rules-reference.md) - Detailed rule documentation
- [Scripts Reference](scripts-reference.md) - Python scripts documentation
- [Demultiplexing](demultiplexing.md) - WarpDemuX setup guide
