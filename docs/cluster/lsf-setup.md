# LSF Setup

Configure the pipeline for LSF (IBM Spectrum LSF) cluster execution.

## Overview

The pipeline includes a pre-configured LSF profile at `cluster/lsf/config.yaml` that handles:

- Job submission via `bsub`
- Memory and queue allocation
- GPU job routing
- Resource limits

## Quick Start

```bash
# Run with LSF profile
pixi run snakemake --profile cluster/lsf --configfile=config/config.yml

# Or use the test shortcut
pixi run test-lsf
```

## Profile Configuration

### Location

```
cluster/lsf/config.yaml
```

### Full Configuration

```yaml
executor: lsf
jobs: 300

default-resources:
  - 'mem_mb=8'
  - 'lsf_queue=rna'
  - 'lsf_project=aatrnaseq'
  - 'lsf_extra=""'

resources:
  - ngpu=12

set-resources:
  - rebasecall:lsf_queue="gpu"
  - rebasecall:lsf_extra="-gpu num=1:j_exclusive=yes"
  - rebasecall:ngpu=1
  - rebasecall:mem_mb=24
  - classify_charging:lsf_queue="gpu"
  - classify_charging:lsf_extra="-gpu num=1:j_exclusive=yes"
  - classify_charging:ngpu=1
  - classify_charging:mem_mb=24
  - remora_signal_stats:mem_mb=24
  - bwa_align:mem_mb=24
  - modkit_extract_calls:mem_mb=96
  - warpdemux:mem_mb=32
  - parse_warpdemux:mem_mb=8

printshellcmds: True
show-failed-logs: True
latency-wait: 15
```

## Configuration Options

### Global Settings

| Option | Value | Description |
|--------|-------|-------------|
| `executor` | `lsf` | Use LSF executor |
| `jobs` | `300` | Maximum concurrent jobs |
| `latency-wait` | `15` | Seconds to wait for file sync |

### Default Resources

Applied to all rules unless overridden:

| Resource | Value | Description |
|----------|-------|-------------|
| `mem_mb` | `8` | Memory in GB |
| `lsf_queue` | `rna` | Default LSF queue |
| `lsf_project` | `aatrnaseq` | Project for accounting |
| `lsf_extra` | `""` | Additional bsub options |

### GPU Resource Limit

```yaml
resources:
  - ngpu=12
```

Limits total concurrent GPU jobs to 12. Adjust based on your cluster's GPU availability.

## Per-Rule Resources

### GPU Rules

These rules are automatically submitted to the GPU queue:

| Rule | Queue | Memory | GPU |
|------|-------|--------|-----|
| `rebasecall` | gpu | 24 GB | 1 (exclusive) |
| `classify_charging` | gpu | 24 GB | 1 (exclusive) |

### Memory-Intensive Rules

| Rule | Memory |
|------|--------|
| `modkit_extract_calls` | 96 GB |
| `warpdemux` | 32 GB |
| `remora_signal_stats` | 24 GB |
| `bwa_align` | 24 GB |

!!! note "escpod demux backend"
    The profile ships resources for the `warpdemux` / `parse_warpdemux` rules. If
    `warpdemux.backend` is `escpod` (the default), add equivalent entries for
    `escpod_demux` and `parse_escpod_demux` — the fused `escpod demux` pass replaces
    both `warpdemux` and `split_pod5`.

## Customization

### Change Default Queue

Edit `lsf_queue` in default-resources:

```yaml
default-resources:
  - 'lsf_queue=your_queue'
```

### Change Project Tag

For job accounting:

```yaml
default-resources:
  - 'lsf_project=your_project'
```

### Adjust Max Jobs

Reduce if you're filling up the queue:

```yaml
jobs: 100
```

### Increase Memory for a Rule

Add or modify in `set-resources`:

```yaml
set-resources:
  - your_rule:mem_mb=64
```

### Change GPU Queue Name

If your GPU queue has a different name:

```yaml
set-resources:
  - rebasecall:lsf_queue="your_gpu_queue"
  - classify_charging:lsf_queue="your_gpu_queue"
```

## Monitoring Jobs

### View All Jobs

```bash
bjobs -u $USER
```

### View Job Details

```bash
bjobs -l <job_id>
```

### View Job History

```bash
bhist -l <job_id>
```

### Kill All Jobs

```bash
bkill 0  # Kill all your jobs
```

### View Queue Status

```bash
bqueues
```

## Submit Scripts

For long-running pipelines, use a submit script:

=== "run-analysis.sh"

    ```bash
    #!/bin/bash
    #BSUB -J aa-tRNA-seq
    #BSUB -o logs/pipeline.%J.out
    #BSUB -e logs/pipeline.%J.err
    #BSUB -q rna
    #BSUB -n 1
    #BSUB -R "rusage[mem=4000]"

    mkdir -p logs

    pixi run snakemake --profile cluster/lsf \
        --configfile=config/config.yml
    ```

Submit:

```bash
bsub < run-analysis.sh
```

## Troubleshooting

### Jobs Pending Too Long

Check queue limits:

```bash
bqueues -l rna
```

Reduce concurrent jobs:

```yaml
jobs: 50
```

### Memory Errors

Increase memory for the failing rule:

```yaml
set-resources:
  - failing_rule:mem_mb=128
```

### GPU Jobs Not Starting

Check GPU queue availability:

```bash
bqueues -l gpu
```

Verify GPU resource syntax for your cluster.

### File Sync Errors

Increase latency wait:

```yaml
latency-wait: 60
```

## Next Steps

- [SLURM Setup](slurm-setup.md) - For SLURM clusters
- [GPU Configuration](gpu-configuration.md) - GPU-specific settings
