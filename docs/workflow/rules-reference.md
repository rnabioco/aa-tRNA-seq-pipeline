# Rules Reference

Complete documentation for all Snakemake rules in the pipeline.

## Processing Rules

These rules form the core data processing pipeline.

### merge_pods

Merge all POD5 files for a sample into a single file.

**File:** `workflow/rules/aatrnaseq-process.smk`

| Property | Value |
|----------|-------|
| Input | All POD5 files for sample |
| Output | `pod5/{sample}/{sample}.pod5` |
| Threads | 12 |
| GPU | No |

**Command:**
```bash
pod5 merge -t {threads} -f -o {output} {input}
```

---

### rebasecall

Re-basecall POD5 files with Dorado, emitting move tables for Remora.

**File:** `workflow/rules/aatrnaseq-process.smk`

| Property | Value |
|----------|-------|
| Input | Merged POD5, mod model sentinel |
| Output | `bam/rebasecall/{sample}/{sample}.rbc.bam` |
| GPU | Yes |
| Parameters | `base_calling_model`, `opts.dorado`, `models_dir` |

**Command:**
```bash
dorado basecaller --models-directory {models_dir} {opts.dorado} {model} {input} > {output}
```

**Notes:**

- Depends on `download_mod_models` rule to pre-download modification models
- Respects `CUDA_VISIBLE_DEVICES` environment variable
- Default options include `--modified-bases pseU m5C inosine_m6A --emit-moves`

---

### ubam_to_fastq

Extract reads from unmapped BAM to FASTQ format.

**File:** `workflow/rules/aatrnaseq-process.smk`

| Property | Value |
|----------|-------|
| Input | Rebasecalled BAM |
| Output | `fq/{sample}/{sample}.fq.gz` |

**Command:**
```bash
samtools fastq -T "*" {input} | gzip > {output}
```

---

### bwa_idx

Build BWA index for reference FASTA.

**File:** `workflow/rules/aatrnaseq-process.smk`

| Property | Value |
|----------|-------|
| Input | Reference FASTA |
| Output | `.amb`, `.ann`, `.bwt`, `.pac`, `.sa` files |

**Command:**
```bash
bwa index {input}
```

**Notes:**

- Only runs once per reference
- Index files are stored alongside the FASTA

---

### bwa_align

Align reads to tRNA reference with BWA MEM.

**File:** `workflow/rules/aatrnaseq-process.smk`

| Property | Value |
|----------|-------|
| Input | FASTQ (EDX-filtered for EDX samples), BWA index |
| Output | `bam/aln/{sample}/{sample}.aln.bam`, `.bai` |
| Threads | 12 |
| Parameters | `fasta`, `opts.bwa` |

**Command:**
```bash
bwa mem -C -t {threads} {opts.bwa} {index} {reads} \
    | samtools view -F 20 -Sb - \
    | samtools sort -o {output}
samtools index {output}
```

**Filtering:**

- `-F 20`: Remove unmapped reads (`0x4`) and reverse-strand reads (`0x10`)

**Default BWA options:**

- `-W 13 -k 6 -T 20 -x ont2d` (RNA-optimized)

---

### classify_charging

Run Remora ML model to classify charged vs uncharged reads.

**File:** `workflow/rules/aatrnaseq-process.smk`

| Property | Value |
|----------|-------|
| Input | POD5, aligned BAM |
| Output | `bam/charging/{sample}/{sample}.charging.bam`, `.bai` |
| Threads | 8 |
| GPU | No (CPU) |
| Parameters | `remora_cca_classifier` |

**Command:**
```bash
remora infer from_pod5_and_bam {pod5} {bam} \
    --model {model} \
    --out-bam {output} \
    --reference-anchored \
    --num-extract-alignment-workers 2 \
    --num-prepare-read-workers 2 \
    --num-prepare-nn-input-workers 2 \
    --num-post-process-workers 2
samtools sort -@ {threads} {output} > {temp}
samtools index {output}
```

**Output tags:**

- `ML`: Modification likelihood (0-255)
- `MM`: Modification metadata

---

### transfer_bam_tags

Transfer and rename charging tags from Remora output to classified BAM.

**File:** `workflow/rules/aatrnaseq-process.smk`

| Property | Value |
|----------|-------|
| Input | Charging BAM, aligned BAM |
| Output | `bam/classified/{sample}/{sample}.bam`, `.bai` |

**Command:**
```bash
python transfer_tags.py \
    --tags ML MM \
    --rename ML=CL MM=CM \
    --source {charging_bam} \
    --target {aligned_bam} \
    --output {output}
samtools index {output}
```

**Tag renaming:**

- `ML` → `CL`: Charging likelihood
- `MM` → `CM`: Charging metadata

This prevents interference with standard SAM modification tags.

---

### add_adapter_tags

Detect adapter positions using parasail alignment and add PT tags to create final BAM.

**File:** `workflow/rules/aatrnaseq-process.smk`

| Property | Value |
|----------|-------|
| Input | Classified BAM |
| Output | `bam/adapter_tagged/{sample}/{sample}.bam`, `.bai` |
| Parameters | `adapters.*` config options |

**Command:**
```bash
python add_adapter_tags.py \
    -i {classified_bam} \
    -o {output} \
    --adapter-5p "{adapter_5p}" \
    --adapter-3p "{adapter_3p}" \
    --min-score-5p {min_score_5p} \
    --min-score-3p {min_score_3p}
samtools index {output}
```

**PT tag format:**
```
PT:Z:start;end;strand;type|start;end;strand;type
Example: PT:Z:0;24;+;5p_adapter|118;135;+;3p_adapter
```

**Notes:**

- Uses parasail Smith-Waterman alignment to find adapter positions
- Can infer 5' adapter presence from alignment position when adapter is truncated
- Output goes to `bam/adapter_tagged/`; the downstream `finalize_bam` rule produces the final BAM at `bam/final/`

---

### finalize_bam

Produce the final BAM for downstream analysis. Symlinks the adapter-tagged BAM as the final output. EDX filtering now happens early in the pipeline (before alignment) via the `detect_edx_adapters` / `filter_fastq_by_edx` / `filter_pod5_by_edx` rules.

**File:** `workflow/rules/aatrnaseq-process.smk`

| Property | Value |
|----------|-------|
| Input | `bam/adapter_tagged/{sample}/{sample}.bam`, `.bai` |
| Output | `bam/final/{sample}/{sample}.bam`, `.bai` |

**Notes:**

- Creates symlinks to the adapter-tagged BAM (zero-copy passthrough)
- This is the final BAM with all tags: CL/CM (charging) and PT (adapters)

---

## Charging Analysis Rules

### get_cca_trna

Extract charging probability (CL tag) per read to TSV.

**File:** `workflow/rules/aatrnaseq-charging.smk`

| Property | Value |
|----------|-------|
| Input | Final BAM |
| Output | `summary/tables/{sample}/{sample}.charging_prob.tsv.gz` |

**Command:**
```bash
python get_charging_table.py --tag CL {input} {output}
```

**Output columns:**

| Column | Description |
|--------|-------------|
| read_id | Nanopore read identifier |
| tRNA | Aligned reference name |
| charging_likelihood | CL tag value (0-255) |

---

### get_cca_trna_cpm

Calculate CPM-normalized charging counts per tRNA.

**File:** `workflow/rules/aatrnaseq-charging.smk`

| Property | Value |
|----------|-------|
| Input | Charging probability TSV |
| Output | `summary/tables/{sample}/{sample}.charging.cpm.tsv.gz` |
| Parameters | `ml_thresh=200` (hardcoded) |

**Command:**
```bash
python get_trna_charging_cpm.py \
    --input {charging_prob} \
    --output {output} \
    --ml-threshold 200
```

**Classification:**

- Score ≥ 200: Charged
- Score < 200: Uncharged

**Output columns:**

| Column | Description |
|--------|-------------|
| tRNA | Reference name |
| counts_charged | Charged read count |
| counts_uncharged | Uncharged read count |
| cpm_charged | Charged CPM |
| cpm_uncharged | Uncharged CPM |

---

## Quality Control Rules

### compute_reference_similarity

Compute pairwise sequence similarity matrix for the reference FASTA.

**File:** `workflow/rules/aatrnaseq-qc.smk`

| Property | Value |
|----------|-------|
| Input | Reference FASTA |
| Output | `summary/qc/reference_similarity.tsv`, `summary/qc/reference_similarity.clusters.tsv` |
| Threads | 4 |
| GPU | No |

**Command:**
```bash
python compute_seq_similarity.py {fasta} {output} \
    --threads {threads} --max-mismatch {max_mismatch} --clusters {clusters}
```

Identical sequences are collapsed before aligning (lossless). Gated by
`qc.reference_similarity` and skipped above `qc.reference_similarity_max_seqs`;
`qc.reference_similarity_max_mismatch` additionally collapses near-identical
sequences by Hamming distance.

**Notes:**

- Uses Needleman-Wunsch global alignment to compute all-vs-all pairwise similarity
- Percent identity = matches / max(len_seq1, len_seq2) * 100
- Identifies potential cross-mapping issues from homologous tRNA sequences

---

### base_calling_error

Extract per-position basecalling error metrics.

**File:** `workflow/rules/aatrnaseq-qc.smk`

| Property | Value |
|----------|-------|
| Input | Final BAM, reference FASTA |
| Output | `summary/tables/{sample}/{sample}.bcerror.tsv.gz` |

**Command:**
```bash
python get_bcerror_freqs.py {bam} {fasta} {output}
```

**Output columns:**

| Column | Description |
|--------|-------------|
| Position | Reference position |
| Coverage | Read coverage |
| A_Freq, T_Freq, G_Freq, C_Freq | Base frequencies |
| MismatchFreq | Mismatch rate |
| InsertionFreq | Insertion rate |
| DeletionFreq | Deletion rate |
| BCErrorFreq | Combined error rate |
| MeanQual | Mean base quality |

---

### align_stats

Summarize alignment statistics across pipeline stages.

**File:** `workflow/rules/aatrnaseq-qc.smk`

| Property | Value |
|----------|-------|
| Input | Unmapped BAM, aligned BAM, classified BAM |
| Output | `summary/tables/{sample}/{sample}.align_stats.tsv.gz` |

**Command:**
```bash
python get_align_stats.py \
    -o {output} \
    -a unmapped aligned classified \
    -i {sample} \
    -b {unmapped} {aligned} {classified}
```

**Output columns:**

| Column | Description |
|--------|-------------|
| bam_file | Source BAM path |
| id | Sample identifier |
| info | Pipeline stage |
| n_reads | Total reads |
| pct_mapped | Percent mapped |
| mean_length | Mean read length |
| mean_bq | Mean base quality |
| mean_mapq | Mean mapping quality |

---

### remora_signal_stats

Extract raw signal metrics using Remora API.

**File:** `workflow/rules/aatrnaseq-qc.smk`

| Property | Value |
|----------|-------|
| Input | Final BAM, POD5 |
| Output | `summary/tables/{sample}/{sample}.remora.tsv.gz` |
| Parameters | `remora_kmer_table`, `opts.remora` |

**Command:**
```bash
python extract_signal_metrics.py \
    --pod5_dir {pod5} \
    --bam {bam} \
    --kmer {kmer_table} \
    --sample_name {sample} \
    | gzip > {output}
```

**Notes:**

- Only runs if `remora_kmer_table` is configured
- Uses custom Remora fork for metrics extraction

---

## Modification Rules

### bam_to_coverage

Generate BedGraph coverage tracks.

**File:** `workflow/rules/aatrnaseq-modifications.smk`

| Property | Value |
|----------|-------|
| Input | Final BAM |
| Output | `summary/tables/{sample}/{sample}.{cpm,counts}.bg.gz` (protected) |
| Threads | 4 |
| Parameters | `opts.coverage` |

**Command:**
```bash
bamCoverage -b {bam} -o {cpm} --normalizeUsing CPM -bs 1 {opts}
bamCoverage -b {bam} -o {counts} -bs 1 {opts}
gzip {outputs}
```

---

### modkit_pileup

Generate per-site modification consensus.

**File:** `workflow/rules/aatrnaseq-modifications.smk`

| Property | Value |
|----------|-------|
| Input | Final BAM |
| Output | `summary/modkit/{sample}/{sample}.pileup.bed.gz` |
| Parameters | `fasta`, modkit thresholds |

**Command:**
```bash
modkit pileup --ref {fasta} {threshold_opts} {bam} - | gzip > {output}
```

---

### modkit_extract_calls

Extract per-read modification calls.

**File:** `workflow/rules/aatrnaseq-modifications.smk`

| Property | Value |
|----------|-------|
| Input | Final BAM |
| Output | `summary/modkit/{sample}/{sample}.mod_calls.tsv.gz` |
| Memory | 96 GB |
| Parameters | `fasta`, modkit thresholds |

**Command:**
```bash
modkit extract calls \
    --bgzf \
    --reference {fasta} \
    --edge-filter 10 \
    --mapped --pass \
    {threshold_opts} \
    {bam} {output}
```

---

### modkit_extract_full

Export comprehensive modification information.

**File:** `workflow/rules/aatrnaseq-modifications.smk`

| Property | Value |
|----------|-------|
| Input | Final BAM |
| Output | `summary/modkit/{sample}/{sample}.mod_full.tsv.gz` |
| Threads | 12 |
| Memory | 48 GB |

**Command:**
```bash
modkit extract full \
    --bgzf \
    --threads 12 \
    --reference {fasta} \
    --edge-filter 10 \
    --mapped \
    {threshold_opts} \
    {bam} {output}
```

---

## Utility Rules

### generate_squiggy_session

Generate a Squiggy session JSON file for loading pipeline outputs in Positron IDE.

**File:** `workflow/rules/common.smk`

| Property | Value |
|----------|-------|
| Input | All final BAMs, all merged POD5s, reference FASTA |
| Output | `squiggy-session.json` (at output root) |
| GPU | No |

**Command:**
```bash
python generate_squiggy_session.py \
    --samples {sample_names} \
    --output-dir {output_dir} \
    --fasta {fasta} \
    --output {output.session}
```

**Notes:**

- Generates absolute paths to POD5, BAM, and FASTA files for each sample
- Computes MD5 checksums for file integrity verification
- Includes default plot options for the Squiggy viewer (eventalign mode, z-normalization)

---

## Odds Ratio Rules

### compute_odds_ratios

Compute per-tRNA pairwise modification odds ratios.

**File:** `workflow/rules/aatrnaseq-odds-ratios.smk`

| Property | Value |
|----------|-------|
| Input | Modkit extract calls TSV, charging probability TSV |
| Output | `summary/tables/{sample}/{sample}.odds_ratios.tsv.gz` |
| Parameters | `odds_ratios.ml_threshold` (default: 200), `odds_ratios.min_coverage` (default: 10) |

**Command:**
```bash
python compute_odds_ratios.py \
    --modkit {modkit_calls} \
    --charging {charging_prob} \
    --output {output} \
    --ml-threshold {ml_thresh} \
    --min-coverage {min_cov}
```

**Notes:**

- For each tRNA, tests whether modification at position X is correlated with modification at position Y (and with charging status) via 2x2 contingency tables
- Uses Haldane correction for zero cells and Fisher's exact test
- Applies BH correction across all results
- Charging status is represented as position 999

**Output columns:**

| Column | Description |
|--------|-------------|
| `tRNA` | Reference tRNA name |
| `pos1` | First position |
| `pos2` | Second position (999 = charging) |
| `n00`, `n01`, `n10`, `n11` | Contingency table counts |
| `total_obs` | Total observations |
| `odds_ratio` | Odds ratio |
| `log_odds_ratio` | Log odds ratio |
| `se_log_or` | Standard error of log OR |
| `ci_lower`, `ci_upper` | 95% confidence interval |
| `fisher_or` | Fisher's exact test OR |
| `p_value` | Fisher's exact test p-value |
| `p_adjusted` | BH-adjusted p-value |

---

## Report Rules

### render_combined_qc_report

Render a combined Quarto QC report with per-sample tabs.

**File:** `workflow/rules/aatrnaseq-report.smk`

| Property | Value |
|----------|-------|
| Input | Alignment stats, charging probabilities, charging CPM, basecalling errors (all samples) |
| Output | `reports/qc_report.html` |
| Parameters | `ml-threshold`, `report.custom_include` |

**Command:**
```bash
quarto render qc-report.qmd \
    -P config_file:{config} \
    -P ml_threshold:{threshold} \
    --output-dir {output_dir}
```

**Notes:**

- Requires the `report` pixi environment: `pixi run -e report snakemake render_combined_qc_report`
- Generates faceted QC plots with per-sample patchwork tabs
- Supports optional custom Quarto include via `report.custom_include` config

---

## Demultiplexing Rules

See [Demultiplexing](demultiplexing.md) for detailed documentation.

### warpdemux

Run WarpDemuX barcode prediction.

| Output | `demux/warpdemux_output/{run_id}/` |

### parse_warpdemux

Parse predictions to barcode mapping file.

| Output | `demux/read_ids/{run_id}/barcode_mapping.tsv.gz` |

### extract_sample_reads

Filter read IDs for specific sample's barcode.

| Output | `demux/read_ids/{sample}.txt` |

### split_pod5

Split merged POD5 by sample using read ID list.

| Output | `demux/pod5/{sample}.pod5` |

### detect_edx_adapters

Detect 3' adapter identity per read on unaligned BAM (before alignment). Only runs for EDX samples.

| Output | `demux/edx/{sample}/{sample}.edx_adapters.tsv.gz` |

### extract_edx_read_ids

Extract read IDs matching the sample's EDX adapter assignment.

| Output | `demux/edx/{sample}/{sample}.edx_read_ids.txt` |

### filter_fastq_by_edx

Extract FASTQ for reads matching the sample's EDX adapter.

| Output | `demux/edx/fq/{sample}/{sample}.fq.gz` |

### filter_pod5_by_edx

Filter POD5 to keep only reads matching the sample's EDX adapter.

| Output | `demux/edx/pod5/{sample}/{sample}.pod5` |

### edx_concordance

Build concordance table of WDX vs EDX adapter identity from adapter detection TSVs.

| Output | `summary/edx/edx_concordance.tsv.gz` |

---

## Rule Dependencies

```mermaid
flowchart LR
    merge_pods --> rebasecall
    rebasecall --> ubam_to_fastq
    rebasecall --> detect_edx_adapters
    detect_edx_adapters --> extract_edx_read_ids
    extract_edx_read_ids --> filter_fastq_by_edx
    extract_edx_read_ids --> filter_pod5_by_edx
    ubam_to_fastq -.-> bwa_align
    filter_fastq_by_edx -.-> bwa_align
    bwa_align --> inject_ubam_tags
    inject_ubam_tags --> classify_charging
    filter_pod5_by_edx -.-> classify_charging
    classify_charging --> transfer_bam_tags
    transfer_bam_tags --> add_adapter_tags
    add_adapter_tags --> finalize_bam
    finalize_bam --> get_cca_trna
    finalize_bam --> base_calling_error
    finalize_bam --> align_stats
    finalize_bam --> bam_to_coverage
    finalize_bam --> modkit_pileup
    finalize_bam --> modkit_extract_calls
    get_cca_trna --> get_cca_trna_cpm
    get_cca_trna --> compute_odds_ratios
    modkit_extract_calls --> compute_odds_ratios
    finalize_bam --> generate_squiggy_session
    merge_pods --> generate_squiggy_session
    finalize_bam --> compute_reference_similarity
    align_stats --> render_combined_qc_report
    get_cca_trna --> render_combined_qc_report
    get_cca_trna_cpm --> render_combined_qc_report
    base_calling_error --> render_combined_qc_report
```

**Note:** Dashed lines (-.->`) indicate conditional paths. For EDX samples, `filter_fastq_by_edx` feeds `bwa_align` and `filter_pod5_by_edx` feeds `classify_charging`. For non-EDX samples, `ubam_to_fastq` feeds `bwa_align` directly.
