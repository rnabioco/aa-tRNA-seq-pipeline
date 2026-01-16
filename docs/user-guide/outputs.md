# Output Files

This guide documents all output files produced by the pipeline.

## Output Directory Structure

```
{output_directory}/
├── pod5/                    # Merged POD5 files
├── bam/                     # BAM files at each stage
├── fq/                      # Extracted FASTQ files
├── summary/                 # Analysis outputs
│   ├── tables/             # Tabular summaries
│   └── modkit/             # Modification calling
├── demux/                   # Demultiplexing outputs (if enabled)
└── logs/                    # Rule execution logs
```

## Data Flow and Outputs

```mermaid
flowchart TB
    subgraph Input
        A[POD5 files]
    end

    subgraph Processing
        B[pod5/{sample}/{sample}.pod5<br/>Merged POD5]
        C[bam/rebasecall/{sample}/{sample}.rbc.bam<br/>Basecalled]
        D[fq/{sample}/{sample}.fq.gz<br/>FASTQ]
        E[bam/aln/{sample}/{sample}.aln.bam<br/>Aligned]
        F[bam/charging/{sample}/{sample}.charging.bam<br/>Classified]
        G[bam/final/{sample}/{sample}.bam<br/>Final BAM]
    end

    subgraph Outputs
        H[summary/tables/<br/>Charging & Stats]
        I[summary/modkit/<br/>Modifications]
    end

    A --> B --> C --> D --> E --> F --> G
    G --> H
    G --> I
```

## Core Outputs

### Final BAM

`bam/final/{sample}/{sample}.bam`

The final BAM file with charging classification and adapter position tags.

**Tags:**

| Tag | Type | Description |
|-----|------|-------------|
| `CL` | `B:C` | Charging likelihood (0-255 scale) |
| `CM` | `Z` | Charging model metadata |
| `PT` | `Z` | Adapter positions (5' and 3' boundaries) |

**View tags:**

```bash
samtools view results/bam/final/sample1/sample1.bam | head -1 | tr '\t' '\n' | grep -E "^(CL|CM|PT):"
```

!!! note "Tag Renaming"
    The original Remora tags ML/MM are renamed to CL/CM to avoid conflicts with standard SAM modification tags.

### Charging Probability Table

`summary/tables/{sample}/{sample}.charging_prob.tsv.gz`

Per-read charging likelihood scores.

| Column | Description |
|--------|-------------|
| `read_id` | Nanopore read identifier |
| `tRNA` | Aligned tRNA reference |
| `charging_likelihood` | ML score (0-255) |

**Interpretation:**

- Score ≥ 200: Charged (aminoacylated)
- Score < 200: Uncharged

**Example:**

```bash
zcat results/summary/tables/sample1/sample1.charging_prob.tsv.gz | head
```

```
read_id                                 tRNA                    charging_likelihood
00a1b2c3-4567-89ab-cdef-0123456789ab   tRNA-Ala-AGC-1-1        245
00a1b2c3-4567-89ab-cdef-0123456789ac   tRNA-Gly-GCC-2-1        87
```

### Charging CPM Table

`summary/tables/{sample}/{sample}.charging.cpm.tsv.gz`

Per-tRNA aggregated charging counts, normalized to CPM (counts per million).

| Column | Description |
|--------|-------------|
| `tRNA` | tRNA reference name |
| `counts_charged` | Number of charged reads |
| `counts_uncharged` | Number of uncharged reads |
| `cpm_charged` | Charged CPM |
| `cpm_uncharged` | Uncharged CPM |

**Example:**

```bash
zcat results/summary/tables/sample1/sample1.charging.cpm.tsv.gz | column -t | head
```

```
tRNA              counts_charged  counts_uncharged  cpm_charged  cpm_uncharged
tRNA-Ala-AGC-1-1  1523            234               15230.5      2340.2
tRNA-Gly-GCC-2-1  892             1456              8920.3       14560.8
```

## Quality Control Outputs

### Alignment Statistics

`summary/tables/{sample}/{sample}.align_stats.tsv.gz`

Read counts through pipeline stages.

| Column | Description |
|--------|-------------|
| `bam_file` | BAM file path |
| `id` | Sample identifier |
| `info` | Pipeline stage |
| `n_reads` | Total reads |
| `pct_mapped` | Percent mapped |
| `mapped_reads` | Number mapped |
| `pos_reads` | Positive strand reads |
| `mapq0_reads` | MAPQ 0 reads |
| `mean_length` | Mean read length |
| `mean_bq` | Mean base quality |
| `mean_mapq` | Mean mapping quality |

### Base Calling Errors

`summary/tables/{sample}/{sample}.bcerror.tsv.gz`

Per-position base calling error metrics.

| Column | Description |
|--------|-------------|
| `Position` | Reference position |
| `Coverage` | Read coverage |
| `A_Freq`, `T_Freq`, `G_Freq`, `C_Freq` | Base frequencies |
| `MismatchFreq` | Mismatch frequency |
| `InsertionFreq` | Insertion frequency |
| `DeletionFreq` | Deletion frequency |
| `BCErrorFreq` | Combined error frequency |
| `MeanQual` | Mean base quality |

### Coverage Tracks

`summary/tables/{sample}/{sample}.{cpm,counts}.bg.gz`

BedGraph coverage tracks for visualization.

- `.cpm.bg.gz` - CPM-normalized coverage
- `.counts.bg.gz` - Raw count coverage

**Load in IGV:**

```bash
gunzip -c results/summary/tables/sample1/sample1.cpm.bg.gz > sample1.cpm.bg
# Load sample1.cpm.bg in IGV
```

## Modification Outputs

### Modification Pileup

`summary/modkit/{sample}/{sample}.pileup.bed.gz`

Per-site modification consensus in BED format.

| Column | Description |
|--------|-------------|
| `chrom` | Reference name |
| `start` | Start position |
| `end` | End position |
| `mod` | Modification type |
| `score` | Modification score |
| `strand` | Strand |
| Additional | Modkit-specific columns |

### Per-Read Modification Calls

`summary/modkit/{sample}/{sample}.mod_calls.tsv.gz`

Individual modification calls per read.

### Full Modification Export

`summary/modkit/{sample}/{sample}.mod_full.tsv.gz`

Comprehensive modification information including all modkit fields.

## Intermediate Files

These files are produced but typically not used directly:

### Merged POD5

`pod5/{sample}/{sample}.pod5`

Merged POD5 file containing all raw signal data for the sample.

### Rebasecalled BAM

`bam/rebasecall/{sample}/{sample}.rbc.bam`

Dorado output with basecalls and move tables. Protected output (not deleted).

### Aligned BAM

`bam/aln/{sample}/{sample}.aln.bam`

BWA MEM alignment output.

### Charging BAM

`bam/charging/{sample}/{sample}.charging.bam`

Remora classification output with ML/MM tags (before renaming).

### FASTQ

`fq/{sample}/{sample}.fq.gz`

Extracted reads for alignment.

## Demultiplexing Outputs

When demultiplexing is enabled:

### Barcode Mapping

`demux/read_ids/{run_id}/barcode_mapping.tsv.gz`

Read ID to barcode assignments.

### Per-Sample Read Lists

`demux/read_ids/{sample}.txt`

Read IDs belonging to each sample.

### Split POD5

`demux/pod5/{sample}.pod5`

Per-sample POD5 files after demultiplexing.

## Signal Metrics (Optional)

If `remora_kmer_table` is configured:

`summary/tables/{sample}/{sample}.remora.tsv.gz`

Remora signal metrics per read per position.

## Log Files

`logs/{rule}/{sample}.log`

Standard output and error for each rule execution.

## File Sizes

Approximate file sizes for a typical sample:

| File | Size |
|------|------|
| Merged POD5 | 5-50 GB |
| Final BAM | 100-500 MB |
| Charging CPM | 10-50 KB |
| Charging Prob | 1-10 MB |
| Modkit pileup | 1-5 MB |

## Cleanup

Remove intermediate files to save space:

```bash
# Remove intermediate BAMs (keep final)
rm -rf results/bam/rebasecall results/bam/aln results/bam/charging

# Remove FASTQ (can be regenerated)
rm -rf results/fq
```

!!! warning "Keep Final Outputs"
    Do not delete `bam/final/`, `summary/`, or `pod5/` directories - these are primary outputs.

## Next Steps

- [Workflow Overview](../workflow/overview.md) - Understand the pipeline stages
- [Rules Reference](../workflow/rules-reference.md) - Detailed rule documentation
