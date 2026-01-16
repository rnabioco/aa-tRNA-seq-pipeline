# Installation

This guide covers installing the aa-tRNA-seq pipeline and its dependencies.

## Prerequisites

- **Operating System**: Linux (tested on CentOS/RHEL 9, Ubuntu 20.04+)
- **Python**: 3.10+
- **GPU**: NVIDIA GPU with CUDA support (required for basecalling and charging classification)
- **Storage**: ~50GB for tools, models, and test data

## Install Pixi

The pipeline uses [Pixi](https://pixi.sh) for environment management.

=== "Linux/macOS"

    ```bash
    curl -fsSL https://pixi.sh/install.sh | sh
    ```

=== "macOS (Homebrew)"

    ```bash
    brew install pixi
    ```

=== "Windows"

    ```powershell
    powershell -ExecutionPolicy Bypass -c "irm -useb https://pixi.sh/install.ps1 | iex"
    ```

After installation, restart your shell or run:

```bash
source ~/.bashrc  # or ~/.zshrc for Zsh
```

For additional installation options, see the [official Pixi installation guide](https://pixi.prefix.dev/latest/installation/).

## Clone the Repository

```bash
git clone https://github.com/rnabioco/aa-tRNA-seq-pipeline.git
cd aa-tRNA-seq-pipeline
```

## Install Dependencies

Install all Python dependencies via Pixi:

```bash
pixi install
```

This creates a `.pixi` directory with all required packages including:

- Snakemake 8.0+
- pysam
- pandas
- pod5
- remora
- samtools
- bwa
- deeptools

## Install External Tools

The pipeline requires several external tools. Install them with a single command:

```bash
pixi run setup
```

This downloads and installs:

- **Dorado** v1.3.1 - Oxford Nanopore basecaller
- **Dorado model** - `rna004_130bps_sup@v5.1.0` basecalling model
- **Remora** - ONT signal analysis for charging classification
- **WarpDemuX** - Barcode demultiplexing (optional, for multiplexed samples)

Dorado and models are installed to `resources/tools/` and `resources/models/`. Modkit is managed by pixi (installed via conda).

## Download Test Data (Optional)

To run the test pipeline, download the test dataset:

```bash
pixi run dl-test-data
```

This downloads ~1GB of test POD5 files to `.tests/`.

## Verify Installation

Verify everything is installed correctly:

```bash
# Check Snakemake version
pixi run snakemake --version

# Check Dorado installation
resources/tools/dorado/1.3.1/bin/dorado --version

# Check Modkit installation (managed by pixi)
pixi run modkit --version

# Dry run with test config
pixi run dry-run
```

## Directory Structure After Installation

```
aa-tRNA-seq-pipeline/
├── .pixi/                    # Pixi environment (includes modkit, remora)
├── resources/
│   ├── tools/
│   │   ├── dorado/1.3.1/    # Dorado binaries
│   │   └── WarpDemuX/       # WarpDemuX (if demux enabled)
│   ├── models/
│   │   ├── rna004_130bps_sup@v5.1.0/  # Basecalling model
│   │   └── cca_classifier.pt          # Remora charging model
│   ├── ref/                  # Reference sequences
│   └── kmers/               # Kmer level tables
├── .tests/                  # Test data (if downloaded)
├── workflow/                # Snakemake workflow
├── config/                  # Configuration files
└── cluster/                 # Cluster profiles
```

## Updating

To update the pipeline:

```bash
git pull
pixi install  # Update dependencies if pixi.lock changed
```

To update external tools, modify the version in `config/config-base.yml` and rerun:

```bash
pixi run setup
```

## Troubleshooting

### Pixi Installation Issues

If Pixi fails to install, ensure you have:

- curl installed
- Write permissions to `~/.pixi`
- Internet access to download packages

### GPU Not Detected

If Dorado fails to detect GPU:

1. Check CUDA is installed: `nvidia-smi`
2. Verify CUDA_VISIBLE_DEVICES is set correctly
3. Ensure GPU drivers are up to date

### Remora Installation Issues

If Remora fails to install with CUDA/PyTorch errors:

```bash
# Manually specify CUDA version
CUDA_VERSION=cu121 pixi run setup
```

## Next Steps

- [Quick Start](quickstart.md) - Run the test pipeline
- [First Analysis](first-analysis.md) - Analyze your own data
