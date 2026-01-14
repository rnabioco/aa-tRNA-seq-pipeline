# SLURM Setup

Configure the pipeline for SLURM or other generic cluster schedulers.

## Overview

The pipeline includes a generic cluster profile at `cluster/generic/config.yaml` that uses Snakemake's `cluster-generic` executor. While configured for LSF's `bsub` command by default, it can be adapted for SLURM.

## Generic Profile (LSF-based)

### Location

```
cluster/generic/config.yaml
```

### Configuration

```yaml
executor: cluster-generic
cluster-generic-submit-cmd:
  bsub
    -o "{log}.out"
    -e "{log}.err"
    -J "{rule}-{wildcards}"
    -R "select[mem>{resources.mem_mb}] rusage[mem={resources.mem_mb}] span[hosts=1]"
    -n {threads}
    -q "{resources.queue}"
    "{resources.gpu_opts}"

default-resources:
  - mem_mb=8
  - queue="rna"
  - gpu_opts=""

jobs: 50

resources:
  - ngpu=8

set-resources:
  - rebasecall:queue="gpu"
  - rebasecall:gpu_opts="-gpu num=1:j_exclusive=yes"
  - rebasecall:ngpu=1
  - rebasecall:mem_mb=24
  - classify_charging:queue="gpu"
  - classify_charging:gpu_opts="-gpu num=1:j_exclusive=yes"
  - classify_charging:ngpu=1
  - classify_charging:mem_mb=24
  - remora_signal_stats:mem_mb=24
  - bwa_align:mem_mb=24
  - modkit_extract_calls:mem_mb=96
  - modkit_extract_full:mem_mb=48

printshellcmds: True
show-failed-logs: True
latency-wait: 60
cluster-generic-cancel-cmd: "bkill"
```

## SLURM Adaptation

To use with SLURM, create a modified profile:

### Create SLURM Profile

```bash
mkdir -p cluster/slurm
```

=== "cluster/slurm/config.yaml"

    ```yaml
    executor: cluster-generic
    cluster-generic-submit-cmd:
      sbatch
        --output="{log}.out"
        --error="{log}.err"
        --job-name="{rule}-{wildcards}"
        --mem={resources.mem_mb}M
        --cpus-per-task={threads}
        --partition="{resources.partition}"
        {resources.gpu_opts}
        --wrap

    default-resources:
      - mem_mb=8000
      - partition="compute"
      - gpu_opts=""

    jobs: 50

    resources:
      - ngpu=8

    set-resources:
      - rebasecall:partition="gpu"
      - rebasecall:gpu_opts="--gres=gpu:1"
      - rebasecall:ngpu=1
      - rebasecall:mem_mb=24000
      - classify_charging:partition="gpu"
      - classify_charging:gpu_opts="--gres=gpu:1"
      - classify_charging:ngpu=1
      - classify_charging:mem_mb=24000
      - remora_signal_stats:mem_mb=24000
      - bwa_align:mem_mb=24000
      - modkit_extract_calls:mem_mb=96000
      - modkit_extract_full:mem_mb=48000

    printshellcmds: True
    show-failed-logs: True
    latency-wait: 60
    cluster-generic-cancel-cmd: "scancel"
    ```

### Key Differences from LSF

| Feature | LSF | SLURM |
|---------|-----|-------|
| Submit command | `bsub` | `sbatch` |
| Cancel command | `bkill` | `scancel` |
| Queue/Partition | `-q queue` | `--partition=partition` |
| Memory | `-R "rusage[mem=X]"` | `--mem=XM` |
| Threads | `-n X` | `--cpus-per-task=X` |
| GPU | `-gpu num=1` | `--gres=gpu:1` |

## Usage

### With Generic Profile

```bash
pixi run snakemake --profile cluster/generic --configfile=config/config.yml
```

### With Custom SLURM Profile

```bash
pixi run snakemake --profile cluster/slurm --configfile=config/config.yml
```

## Configuration Options

### Max Concurrent Jobs

```yaml
jobs: 50
```

Adjust based on your cluster's fair share policy.

### Latency Wait

```yaml
latency-wait: 60
```

Increase for network file systems with slow sync.

### GPU Limits

```yaml
resources:
  - ngpu=8
```

Limits concurrent GPU jobs. Set to your available GPUs.

## Per-Rule Resources

### Memory Requirements

| Rule | Memory (MB) |
|------|-------------|
| `rebasecall` | 24000 |
| `classify_charging` | 24000 |
| `modkit_extract_calls` | 96000 |
| `modkit_extract_full` | 48000 |
| `remora_signal_stats` | 24000 |
| `bwa_align` | 24000 |

### GPU Rules

| Rule | Partition | GPU |
|------|-----------|-----|
| `rebasecall` | gpu | 1 |
| `classify_charging` | gpu | 1 |

## Monitoring Jobs

### SLURM Commands

```bash
# View your jobs
squeue -u $USER

# View job details
scontrol show job <job_id>

# Cancel job
scancel <job_id>

# Cancel all your jobs
scancel -u $USER

# View partition status
sinfo
```

### LSF Commands (Generic Profile)

```bash
# View your jobs
bjobs -u $USER

# View job details
bjobs -l <job_id>

# Cancel job
bkill <job_id>

# View queue status
bqueues
```

## Submit Scripts

### SLURM Submit Script

=== "run-slurm.sh"

    ```bash
    #!/bin/bash
    #SBATCH --job-name=aa-tRNA-seq
    #SBATCH --output=logs/pipeline.%j.out
    #SBATCH --error=logs/pipeline.%j.err
    #SBATCH --partition=compute
    #SBATCH --mem=4G
    #SBATCH --cpus-per-task=1

    mkdir -p logs

    pixi run snakemake --profile cluster/slurm \
        --configfile=config/config.yml
    ```

Submit:

```bash
sbatch run-slurm.sh
```

## Troubleshooting

### Jobs Not Starting

Check partition limits:

```bash
sinfo -p compute
```

### Memory Errors

SLURM uses different memory units. Ensure values are in MB:

```yaml
set-resources:
  - rule:mem_mb=24000  # 24 GB
```

### GPU Not Detected

Verify GPU resource specification for your cluster:

```yaml
set-resources:
  - rebasecall:gpu_opts="--gres=gpu:1"
```

Or for specific GPU types:

```yaml
set-resources:
  - rebasecall:gpu_opts="--gres=gpu:v100:1"
```

## Next Steps

- [LSF Setup](lsf-setup.md) - For LSF clusters
- [GPU Configuration](gpu-configuration.md) - GPU-specific settings
