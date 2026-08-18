# First Analysis

This guide walks through analyzing your own aa-tRNA-seq data.

## Prerequisites

- Completed [Installation](installation.md)
- POD5 files from your Oxford Nanopore sequencing run
- Reference tRNA sequences with adapters (or use the provided S. cerevisiae reference)

## Step 1: Organize Your Data

The pipeline expects POD5 files organized by sequencing run:

```
/path/to/your/data/
├── run1/
│   ├── pod5_pass/     # Passed reads (required)
│   │   ├── file1.pod5
│   │   └── file2.pod5
│   └── pod5_fail/     # Failed reads (optional)
│       └── file3.pod5
└── run2/
    └── pod5/          # Alternative directory name
        └── file4.pod5
```

The pipeline searches for POD5 files in `pod5_pass/`, `pod5_fail/`, or `pod5/` subdirectories.

## Step 2: Create Sample File

Create a tab-separated sample file listing your samples:

=== "samples.tsv"

    ```
    sample1	/path/to/your/data/run1
    sample2	/path/to/your/data/run2
    sample1	/path/to/your/data/run1_replicate
    ```

!!! note "Multiple Runs per Sample"
    You can list the same sample ID multiple times with different run paths. POD5 files will be merged before processing.

For multiplexed samples with barcodes, use YAML format instead. See [Sample Files](../user-guide/samples-format.md).

## Step 3: Create Configuration

Create a configuration file for your analysis:

=== "config/config-myproject.yml"

    ```yaml
    # Sample file path
    samples: config/samples-myproject.tsv

    # Output directory
    output_directory: "results/myproject"

    # Optional: Override base configuration
    # fasta: "path/to/custom/reference.fa"
    # charging:
    #   model: "path/to/custom/bundle"
    ```

The configuration inherits defaults from `config/config-base.yml`. Override any parameter as needed.

## Step 4: Dry Run

Preview the pipeline execution:

```bash
pixi run snakemake -n --configfile=config/config-myproject.yml
```

Verify:

- All samples are detected
- Expected number of jobs
- No errors in DAG construction

## Step 5: Run Pipeline

### Local Execution

For small datasets or testing:

```bash
pixi run snakemake --cores 12 --configfile=config/config-myproject.yml
```

Adjust `--cores` based on your system.

### Cluster Execution (Recommended)

For production runs, use cluster execution:

=== "LSF"

    ```bash
    pixi run snakemake --profile cluster/lsf --configfile=config/config-myproject.yml
    ```

=== "SLURM"

    ```bash
    pixi run snakemake --profile cluster/slurm --configfile=config/config-myproject.yml
    ```

See [Cluster Setup](../cluster/lsf-setup.md) for detailed cluster configuration.

### Background Execution

For long-running jobs, run in background:

```bash
nohup pixi run snakemake --profile cluster/lsf \
    --configfile=config/config-myproject.yml \
    > pipeline.log 2>&1 &
```

Monitor progress:

```bash
tail -f pipeline.log
```

## Step 6: Monitor Progress

### Check Job Status

```bash
# For LSF
bjobs -u $USER

# For SLURM
squeue -u $USER
```

### View Snakemake Progress

```bash
# Show running jobs
pixi run snakemake --profile cluster/lsf \
    --configfile=config/config-myproject.yml \
    --summary
```

### Check Logs

Rule-specific logs are in the output directory:

```bash
ls results/myproject/logs/
```

## Step 7: Verify Outputs

After completion, check your outputs:

```bash
ls -la results/myproject/summary/tables/
```

### Key Files to Check

1. **Charging CPM** - Per-tRNA charging quantification:
   ```bash
   zcat results/myproject/summary/tables/sample1/sample1.charging.cpm.tsv.gz | column -t | head
   ```

2. **Alignment Stats** - Read counts through pipeline:
   ```bash
   zcat results/myproject/summary/tables/sample1/sample1.align_stats.tsv.gz | column -t
   ```

3. **Final BAM** - Verify charging tags:
   ```bash
   samtools view results/myproject/bam/final/sample1/sample1.bam | head -1
   ```

## Troubleshooting

### No POD5 Files Found

```
Error: No POD5 files found for sample 'sample1'
```

**Solution**: Verify POD5 files exist in `pod5_pass/`, `pod5_fail/`, or `pod5/` subdirectories.

### GPU Memory Error

```
CUDA out of memory
```

**Solution**: Reduce batch size or run one GPU job at a time. Adjust cluster profile `ngpu` resource.

### Rule Failed

```
Error in rule rebasecall
```

**Solution**: Check rule-specific log:

```bash
cat results/myproject/logs/rebasecall/sample1.log
```

See [Troubleshooting](../troubleshooting/common-errors.md) for more solutions.

## Next Steps

- [Output Files](../user-guide/outputs.md) - Understand all output files
- [Configuration](../user-guide/configuration.md) - Customize parameters
- [Workflow Overview](../workflow/overview.md) - Understand the pipeline stages
