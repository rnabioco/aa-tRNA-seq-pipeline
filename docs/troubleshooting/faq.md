# FAQ

Frequently asked questions about the aa-tRNA-seq pipeline.

## General Questions

### What species can I analyze?

The pipeline is designed for *Saccharomyces cerevisiae* tRNAs by default. To use with other species:

1. Create a reference FASTA with your tRNA sequences + adapters
2. Train or adapt the Remora charging model for your species
3. Update the config to point to your reference

### What is the charging threshold?

The default threshold is **ML ≥ 200 = charged**. This is currently hardcoded in the `get_cca_trna_cpm` rule. Values range from 0-255.

### Can I change the charging threshold?

Yes, but requires editing the rule. In `workflow/rules/aatrnaseq-charging.smk`:

```python
params:
    ml_thresh=200,  # Change this value
```

### What modifications are detected?

The pipeline calls these modifications via Dorado:

- **pseU** - Pseudouridine
- **m5C** - 5-methylcytosine
- **m6A** - N6-methyladenosine
- **inosine** - Inosine

### How accurate is charging classification?

Accuracy depends on:

- Signal quality
- Full-length read retention
- Model training data match

The 6-nucleotide kmer (CCAGGC) provides the classification signal.

---

## Data Requirements

### What input format is required?

**POD5 format** from Oxford Nanopore sequencing. The pipeline does not accept:

- FAST5 (convert with `pod5 convert`)
- FASTQ (signal data required)
- BAM (unless from Dorado with move tables)

### How should I organize my data?

```
run_directory/
├── pod5_pass/    # Passed reads (required)
│   └── *.pod5
├── pod5_fail/    # Failed reads (optional)
│   └── *.pod5
└── pod5/         # Alternative name
    └── *.pod5
```

### Can I merge multiple runs per sample?

Yes! List the same sample ID multiple times in your samples file:

```
sample1    /data/run1
sample1    /data/run2
```

POD5 files from both runs will be merged.

### What's the minimum read count?

There's no hard minimum, but:

- **100+ reads** per tRNA recommended for reliable CPM estimates
- **1000+ reads** recommended for statistical analysis

---

## Execution Questions

### How long does the pipeline take?

Typical times per sample (GPU):

| Step | Time |
|------|------|
| Merge POD5 | 5-15 min |
| Rebasecall | 30-60 min |
| Alignment | 5-10 min |
| Classification | 10-30 min |
| Summaries | 5-10 min |
| **Total** | **1-2 hours** |

CPU-only is 10-100x slower.

### Can I run without a GPU?

Yes, but not recommended. Dorado and Remora will fall back to CPU:

```bash
export CUDA_VISIBLE_DEVICES=""
pixi run snakemake --cores 12 --configfile=config/config.yml
```

Expect significant slowdown.

### How much disk space is needed?

Per sample (approximate):

| Data | Size |
|------|------|
| Merged POD5 | 5-50 GB |
| Rebasecalled BAM | 1-5 GB |
| Final BAM | 100-500 MB |
| Summary tables | 10-100 MB |

Plan for ~50-100 GB per sample during processing. Clean intermediate files after.

### Can I resume a failed run?

Yes! Snakemake automatically resumes:

```bash
pixi run snakemake --profile cluster/lsf --configfile=config/config.yml
```

Only incomplete rules will re-run.

---

## Output Questions

### What does the charging_prob table contain?

Per-read charging likelihood scores:

| Column | Description |
|--------|-------------|
| read_id | Nanopore read ID |
| tRNA | Aligned tRNA reference |
| charging_likelihood | Score 0-255 (≥200 = charged) |

### What does the CPM table contain?

Per-tRNA aggregated counts:

| Column | Description |
|--------|-------------|
| tRNA | Reference name |
| counts_charged | Charged read count |
| counts_uncharged | Uncharged read count |
| cpm_charged | Charged CPM |
| cpm_uncharged | Uncharged CPM |

### Why are CL/CM tags used instead of ML/MM?

The original Remora tags (ML/MM) are renamed to CL/CM to avoid conflicts with standard SAM modification tags used by Dorado.

### How do I visualize results?

Load output files in:

- **IGV**: BAM files, BedGraph coverage
- **R/Python**: TSV tables for statistical analysis
- **Downstream repo**: [github.com/rnabioco/aa-tRNA-seq](https://github.com/rnabioco/aa-tRNA-seq)

---

## Configuration Questions

### What's the difference between config files?

| File | Purpose |
|------|---------|
| `config-base.yml` | Default parameters (always loaded) |
| `config-test.yml` | Test data paths |
| `config-preprint.yml` | Full preprint analysis |
| `config-demux-test.yml` | Demultiplexing test |

Your config only needs to override what's different.

### How do I add a new sample?

1. Add to samples.tsv:
   ```
   new_sample    /path/to/data
   ```

2. Run pipeline:
   ```bash
   pixi run snakemake --configfile=config/config.yml
   ```

Only new samples will be processed.

### Can I use a custom reference?

Yes, set in your config:

```yaml
fasta: "path/to/your/reference.fa"
```

The reference should include adapter sequences matching your library prep.

---

## Demultiplexing Questions

### When should I use demultiplexing?

When samples were **pooled with WarpDemuX barcodes** during library preparation.

### Which barcode kit should I use?

| Kit | Barcodes | Recommendation |
|-----|----------|----------------|
| `WDX4_tRNA_rna004_v1_0` | 03, 04, 05, 07 | **Recommended** |
| `WDX4b_tRNA_rna004_v1_0` | 04, 05, 07, 11 | Alternative |

### Why are my demux reads unbalanced?

Common causes:

1. Library concentration differences
2. Barcode ligation efficiency variation
3. Loading bias

Check `demux_summary.tsv.gz` for distribution.

### Can I mix demux and direct samples?

Yes, in YAML format use `~` for direct samples:

```yaml
runs:
  - path: /data/pooled
    samples:
      pooled_sample: "barcode03"
  - path: /data/direct
    samples:
      direct_sample: ~
```

---

## Troubleshooting Questions

### Why did my job fail?

1. Check the log file:
   ```bash
   cat results/logs/<rule>/<sample>
   ```

2. Check cluster job output:
   ```bash
   # LSF
   bjobs -l <job_id>
   ```

### How do I force re-run a rule?

```bash
pixi run snakemake --forcerun <rule_name> --configfile=config/config.yml
```

### How do I clean up and start fresh?

```bash
# Remove all outputs
rm -rf results/

# Or remove specific intermediate files
rm -rf results/bam/rebasecall results/bam/aln results/fq
```

### Where can I get help?

- **GitHub Issues**: [github.com/rnabioco/aa-tRNA-seq-pipeline/issues](https://github.com/rnabioco/aa-tRNA-seq-pipeline/issues)
- **Documentation**: This site
- **Downstream analysis**: [github.com/rnabioco/aa-tRNA-seq](https://github.com/rnabioco/aa-tRNA-seq)
