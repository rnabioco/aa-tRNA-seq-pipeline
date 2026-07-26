# Quick Start

Run the pipeline with test data to verify your installation.

## Prerequisites

Complete the [Installation](installation.md) guide first:

```bash
pixi install
pixi run setup
pixi run dl-test-data
```

`pixi run setup` builds the `escpod` CLI from source (needs a Rust toolchain >= 1.95)
and installs `leech` from private release wheels (needs an authenticated `gh`, or
`LEECH_WHEEL_DIR`). Neither the git clone nor setup involves submodules.

## Dry Run

Preview what the pipeline will do without executing:

```bash
pixi run dry-run
```

Expected output:

```
Building DAG of jobs...
Job stats:
job                      count
---------------------  -------
all                          1
bwa_align                    2
classify_charging            2
get_cca_trna                 2
get_cca_trna_cpm             2
merge_pods                   2
rebasecall                   2
transfer_bam_tags            2
...
total                       XX
```

## Run Test Pipeline

Execute the pipeline locally with 4 cores:

```bash
pixi run test
```

This processes two test samples and takes approximately 15-30 minutes depending on your hardware (GPU required).

!!! note "GPU Required"
    The `rebasecall` and `classify_charging` rules require GPU access. Ensure CUDA is available.

## Expected Output

After completion, outputs are in `.tests/outputs/`:

```
.tests/outputs/
├── pod5/
│   └── sample1/
│       └── sample1.pod5              # Merged POD5
├── bam/
│   ├── rebasecall/sample1/
│   │   └── sample1.rbc.bam           # Basecalled BAM
│   ├── aln/sample1/
│   │   └── sample1.aln.bam           # Aligned BAM
│   ├── charging/sample1/
│   │   └── sample1.charging.bam      # Remora classification
│   └── final/sample1/
│       └── sample1.bam               # Final BAM with CL/CM/PT tags
├── fq/
│   └── sample1/
│       └── sample1.fq.gz             # Extracted FASTQ
└── summary/
    ├── tables/sample1/
    │   ├── sample1.charging_prob.tsv.gz  # Per-read charging
    │   ├── sample1.charging.cpm.tsv.gz   # CPM counts
    │   ├── sample1.bcerror.tsv.gz        # Base errors
    │   └── sample1.align_stats.tsv.gz    # Alignment stats
    └── modkit/sample1/
        ├── sample1.pileup.bed.gz         # Modification pileup
        └── sample1.mod_calls.tsv.gz      # Per-read mods
```

## Inspect Key Outputs

### Charging CPM Table

View per-tRNA charging counts:

```bash
zcat .tests/outputs/summary/tables/sample1/sample1.charging.cpm.tsv.gz | head
```

Output columns:

| Column | Description |
|--------|-------------|
| tRNA | tRNA reference name |
| counts_charged | Number of charged reads |
| counts_uncharged | Number of uncharged reads |
| cpm_charged | Charged counts per million |
| cpm_uncharged | Uncharged counts per million |

### Per-Read Charging Probabilities

View individual read classifications:

```bash
zcat .tests/outputs/summary/tables/sample1/sample1.charging_prob.tsv.gz | head
```

Output columns:

| Column | Description |
|--------|-------------|
| read_id | Nanopore read ID |
| tRNA | Aligned tRNA reference |
| charging_likelihood | ML score (0-255; ≥200 = charged) |

### Final BAM Tags

The final BAM contains charging classification in tags:

```bash
samtools view .tests/outputs/bam/final/sample1/sample1.bam | head -1 | tr '\t' '\n' | grep -E "^(CL|CM|PT):"
```

- `CL:B:C` - Charging likelihood (ML tag renamed to avoid conflict)
- `CM:Z` - Charging model metadata (MM tag renamed)
- `PT:Z` - Adapter positions (5' and 3' adapter boundaries)

## Run on Cluster

For cluster execution with LSF:

```bash
pixi run test-lsf
```

Or manually with a profile:

```bash
pixi run snakemake --profile cluster/lsf --configfile=config/config-test.yml
```

See [LSF Setup](../cluster/lsf-setup.md) for cluster configuration.

## Cleanup

Remove test outputs:

```bash
rm -rf .tests/outputs
```

## Next Steps

- [First Analysis](first-analysis.md) - Use your own data
- [Configuration](../user-guide/configuration.md) - Customize parameters
- [Output Files](../user-guide/outputs.md) - Detailed output documentation
