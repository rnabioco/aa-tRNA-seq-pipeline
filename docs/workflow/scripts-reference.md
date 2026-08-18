# Scripts Reference

Documentation for Python scripts in `workflow/scripts/`.

## Overview

| Script | Purpose |
|--------|---------|
| `transfer_tags.py` | Transfer BAM tags between files |
| `get_charging_table.py` | Extract charging likelihood per read |
| `get_trna_charging_cpm.py` | Calculate CPM-normalized counts |
| `get_charging_summary.py` | Generate charging statistics |
| `get_bcerror_freqs.py` | Calculate basecalling error metrics |
| `get_align_stats.py` | Summarize alignment statistics |
| `read_attrition.py` | Report where the run's reads were lost |
| `anchor_coverage.py` | Count aligned reads spanning the CCA anchor |
| `filter_reads.py` | Filter BAM by quality criteria |
| `generate_squiggy_session.py` | Generate Squiggy session JSON for Positron |
| `compute_odds_ratios.py` | Compute per-tRNA pairwise modification odds ratios |
| `compute_seq_similarity.py` | Compute pairwise reference sequence similarity |
| `collapse_gtrndb_fasta.py` | Collapse redundant GtRNAdb FASTA sequences |

---

## transfer_tags.py

Transfer specified BAM tags from one file to another with optional renaming.

### Usage

```bash
python transfer_tags.py \
    --tags ML MM \
    --rename ML=CL MM=CM \
    --source source.bam \
    --target target.bam \
    --output output.bam
```

### Arguments

| Argument | Description |
|----------|-------------|
| `--tags` | Tags to transfer (space-separated) |
| `--rename` | Tag renaming (OLD=NEW format) |
| `--source` | Source BAM with tags |
| `--target` | Target BAM to add tags to |
| `--output` | Output BAM path |

### Behavior

1. Reads source BAM, stores specified tags by read ID
2. Iterates target BAM, transfers matching tags
3. Outputs only primary alignments (filters secondary/supplementary)

---

## get_charging_table.py

Extract charging likelihood tag values from BAM to TSV format.

### Usage

```bash
python get_charging_table.py --tag CL input.bam output.tsv.gz
```

### Arguments

| Argument | Description |
|----------|-------------|
| `--tag` | Tag to extract (default: CL) |
| `input` | Input BAM file |
| `output` | Output TSV file (gzipped) |

### Output Format

```
read_id                                  tRNA              charging_likelihood
00a1b2c3-4567-89ab-cdef-0123456789ab    tRNA-Ala-AGC-1-1  245
```

| Column | Description |
|--------|-------------|
| read_id | Nanopore read identifier |
| tRNA | Reference name (aligned position) |
| charging_likelihood | Tag value (0-255) |

---

## get_trna_charging_cpm.py

Aggregate per-read charging data to tRNA-level CPM counts.

### Usage

```bash
python get_trna_charging_cpm.py \
    --input charging_prob.tsv.gz \
    --output charging.cpm.tsv.gz \
    --ml-threshold 200
```

### Arguments

| Argument | Description |
|----------|-------------|
| `--input` | Input charging probability TSV |
| `--output` | Output CPM TSV |
| `--ml-threshold` | Threshold for charged classification (default: 200) |

### Classification Logic

```python
if charging_likelihood >= threshold:
    classification = "charged"
else:
    classification = "uncharged"
```

### Output Format

| Column | Description |
|--------|-------------|
| tRNA | Reference name |
| counts_charged | Number of charged reads |
| counts_uncharged | Number of uncharged reads |
| cpm_charged | (counts_charged / total) × 1,000,000 |
| cpm_uncharged | (counts_uncharged / total) × 1,000,000 |

---

## get_charging_summary.py

Generate tRNA-gene level charging statistics.

### Usage

```bash
python get_charging_summary.py \
    -b sample.bam \
    -t trna_table.txt \
    -o output.tsv
```

### Arguments

| Argument | Description |
|----------|-------------|
| `-b` | Input BAM file |
| `-t` | tRNA lookup table |
| `-o` | Output TSV file |

### Lookup Table Format

4 whitespace-delimited columns (no header):

```
uncharged_seq    charged_seq    isodecoder    gene_id
tRNA-Ala-AGC-u   tRNA-Ala-AGC-c Ala-AGC       tRNA-Ala-1
```

### Output Format

| Column | Description |
|--------|-------------|
| tRNA-gene | Gene identifier |
| charged_read_counts | Charged reads |
| uncharged_read_counts | Uncharged reads |
| percent_charged | (charged / total) × 100 |

---

## get_bcerror_freqs.py

Calculate per-position basecalling error frequencies.

### Usage

```bash
python get_bcerror_freqs.py sample.bam reference.fa output.tsv.gz
```

### Arguments

| Position | Description |
|----------|-------------|
| 1 | Input BAM file |
| 2 | Reference FASTA |
| 3 | Output TSV file |

### Metrics Calculated

| Metric | Formula |
|--------|---------|
| Coverage | Reads spanning position |
| Bases_Mapped | Aligned bases at position |
| MismatchFreq | mismatches / bases_mapped |
| InsertionFreq | insertions / bases_mapped |
| DeletionFreq | deletions / coverage |
| BCErrorFreq | (mismatch + ins + del) / coverage |
| MeanQual | Mean base quality score |

### Output Format

| Column | Description |
|--------|-------------|
| Reference | Reference name |
| Position | 1-based position |
| RefBase | Reference base |
| Coverage | Read coverage |
| Bases_Mapped | Bases aligned |
| A_Freq, T_Freq, G_Freq, C_Freq, N_Freq | Base frequencies |
| MismatchFreq | Mismatch frequency |
| InsertionFreq | Insertion frequency |
| DeletionFreq | Deletion frequency |
| BCErrorFreq | Combined error frequency |
| MeanQual | Mean base quality |

---

## get_align_stats.py

Summarize alignment statistics across multiple BAM files.

### Usage

```bash
python get_align_stats.py \
    -o stats.tsv.gz \
    -a unmapped aligned classified \
    -i sample_id \
    -b ubam.bam aligned.bam classified.bam
```

### Arguments

| Argument | Description |
|----------|-------------|
| `-o` | Output TSV file |
| `-a` | Labels for each BAM (space-separated) |
| `-i` | Sample identifier |
| `-b` | BAM files (space-separated) |

### Classes

**CountArray**: Memory-efficient integer count tracking

- `push(index)`: Increment count at index
- `nth(n)`: Get nth value
- `mean()`: Calculate mean
- `quantile(p)`: Calculate percentile

**ReadStats**: Aggregate read metrics

- Tracks read length, base quality, MAPQ distributions
- Returns summary with mean/median for each metric

### Output Format

| Column | Description |
|--------|-------------|
| bam_file | Source BAM path |
| id | Sample identifier |
| info | Pipeline stage label |
| n_reads | Total reads |
| pct_mapped | Percent mapped |
| mapped_reads | Mapped read count |
| pos_reads | Positive strand reads |
| mapq0_reads | MAPQ 0 reads |
| mapq_pass_reads | MAPQ > 0 reads |
| mean_length | Mean read length |
| median_length | Median read length |
| mean_bq | Mean base quality |
| median_bq | Median base quality |
| mean_mapq | Mean mapping quality |
| median_mapq | Median mapping quality |

---

## filter_reads.py

Filter BAM alignments by quality criteria.

!!! note "Not used in main pipeline"
    This script is available for custom filtering but not called by default rules.

### Usage

```bash
python filter_reads.py \
    -i input.bam \
    -o output.bam \
    -5 24 \
    -3 23 \
    -s
```

### Arguments

| Argument | Description |
|----------|-------------|
| `-i` | Input BAM |
| `-o` | Output BAM |
| `-5` | 5' truncation threshold |
| `-3` | 3' truncation threshold |
| `-s` | Require positive strand |
| `-q` | Minimum MAPQ |
| `--multi` | Filter multi-mapping reads |
| `--adapter-dist` | Maximum adapter edit distance |

### Filter Flags

The script adds a `zf` tag with bitwise filter flags:

| Bit | Filter |
|-----|--------|
| 1 | 5' truncation |
| 2 | 3' truncation |
| 4 | Low MAPQ |
| 8 | Negative strand |
| 16 | Secondary alignment |
| 32 | Supplementary alignment |
| 64 | Multi-mapped |
| 128 | Adapter distance |

### FilterStats Class

Tracks filtering reasons with methods:

- `add(reason)`: Record filter reason
- `summary()`: Get filter statistics
- `write_stats(path)`: Write statistics to file

---

## generate_squiggy_session.py

Generate a Squiggy session JSON file for loading pipeline outputs in Positron IDE.

### Usage

```bash
python generate_squiggy_session.py \
    --samples sample1 sample2 \
    --output-dir /path/to/output \
    --fasta /path/to/reference.fa \
    --output squiggy-session.json \
    [--session-name "My session"] \
    [--no-checksums]
```

### Arguments

| Argument | Description |
|----------|-------------|
| `--samples` | Sample names to include (space-separated) |
| `--output-dir` | Pipeline output directory (paths are relative to this) |
| `--fasta` | Path to reference FASTA file |
| `--output` | Output JSON file path |
| `--session-name` | Optional session name (defaults to `aa-tRNA-seq: {directory}`) |
| `--no-checksums` | Skip computing MD5 checksums (faster but no integrity verification) |

### Output Format

JSON with the following structure:

| Key | Description |
|-----|-------------|
| `version` | Schema version (`1.0.0`) |
| `timestamp` | ISO 8601 generation time |
| `sessionName` | Display name for the session |
| `samples` | Per-sample `pod5Paths`, `bamPath`, `fastaPath` (absolute paths) |
| `plotOptions` | Default plot settings (eventalign mode, z-normalization) |
| `fileChecksums` | MD5, size, and last-modified per file (unless `--no-checksums`) |

---

## compute_odds_ratios.py

Compute per-tRNA pairwise modification odds ratios from modkit calls and charging data.

### Usage

```bash
python compute_odds_ratios.py \
    --modkit mod_calls.tsv.gz \
    --charging charging_prob.tsv.gz \
    --output odds_ratios.tsv.gz \
    --ml-threshold 200 \
    --min-coverage 10
```

### Arguments

| Argument | Description |
|----------|-------------|
| `--modkit` | Path to modkit extract calls TSV (.gz) |
| `--charging` | Path to charging probability TSV (.gz) |
| `--output` | Output path (.tsv.gz) |
| `--ml-threshold` | CL tag threshold for charged classification (default: 200) |
| `--min-coverage` | Minimum reads per tRNA or position pair (default: 10) |

### Behavior

1. Loads modkit per-read modification calls and binarizes (canonical vs modified)
2. Loads charging probabilities and binarizes at threshold
3. For each tRNA, pivots to per-read matrix (rows=reads, columns=positions)
4. Adds charging as position 999
5. Tests all pairwise combinations via 2x2 contingency tables
6. Applies Haldane correction for zero cells, Fisher's exact test
7. Applies BH correction across all results

### Output Format

| Column | Description |
|--------|-------------|
| tRNA | Reference tRNA name |
| pos1 | First position |
| pos2 | Second position (999 = charging) |
| n00, n01, n10, n11 | Contingency table counts |
| total_obs | Total observations |
| odds_ratio | Odds ratio |
| log_odds_ratio | Log odds ratio |
| se_log_or | Standard error of log OR |
| ci_lower, ci_upper | 95% confidence interval |
| fisher_or | Fisher's exact test OR |
| p_value | Fisher's exact test p-value |
| p_adjusted | BH-adjusted p-value |

---

## compute_seq_similarity.py

Compute pairwise sequence similarity matrix for reference FASTA sequences.

### Usage

```bash
python compute_seq_similarity.py reference.fa output.tsv
```

### Arguments

| Position | Description |
|----------|-------------|
| 1 | Input reference FASTA file |
| 2 | Output TSV file |

### Behavior

1. Reads all sequences from FASTA using pysam
2. Computes all-vs-all pairwise Needleman-Wunsch global alignment using parasail
3. Calculates percent identity = matches / max(len_seq1, len_seq2) * 100
4. Writes square TSV similarity matrix

### Output Format

Square TSV matrix with sequence names as row and column headers. Values are percent identity (0-100). Diagonal entries are always 100.0.

---

## collapse_gtrndb_fasta.py

Collapse redundant sequences in a GtRNAdb FASTA reference.

### Usage

```bash
python collapse_gtrndb_fasta.py \
    -i ecoliK12-mature-tRNAs.fa \
    -o ecoliK12-collapsed.fa \
    -m ecoliK12-mapping.tsv \
    [--keep-unparsed]
```

### Arguments

| Argument | Description |
|----------|-------------|
| `-i`, `--input` | Input GtRNAdb FASTA file |
| `-o`, `--output` | Output collapsed FASTA file |
| `-m`, `--mapping` | Output TSV mapping report |
| `--keep-unparsed` | Pass through sequences with non-GtRNAdb headers (default: error) |

### Behavior

1. Parses GtRNAdb FASTA headers using end-anchored regex to extract tRNA identity
2. Groups sequences by isodecoder family (amino acid + anticodon)
3. Strips trailing CCA for sequence comparison
4. Identifies redundant gene copies with identical sequences
5. Writes collapsed FASTA (first occurrence as representative) and TSV mapping

### Mapping Output Format

| Column | Description |
|--------|-------------|
| collapsed_name | Representative sequence name |
| original_name | Original GtRNAdb header |
| isodecoder | Isodecoder family key |
| amino_acid | Amino acid identity |
| anticodon | Anticodon sequence |
| family_num | Family number |
| copy_num | Gene copy number |
| is_representative | Whether this is the representative sequence |
| sequence_length | Sequence length |

---

## Common Patterns

### Reading BAM Files

```python
import pysam

with pysam.AlignmentFile(bam_path, "rb") as bam:
    for read in bam.fetch():
        # Access tags
        tag_value = read.get_tag("CL")
        # Access alignment info
        ref_name = read.reference_name
        read_id = read.query_name
```

### Writing Compressed Output

```python
import gzip

with gzip.open(output_path, "wt") as f:
    f.write("column1\tcolumn2\n")
    f.write(f"{value1}\t{value2}\n")
```

### Using pandas

```python
import pandas as pd

df = pd.read_csv(input_path, sep="\t", compression="gzip")
df.to_csv(output_path, sep="\t", index=False, compression="gzip")
```
