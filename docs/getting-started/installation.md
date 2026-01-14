# Installation

This guide covers installing the aa-tRNA-seq pipeline and its dependencies.

## Prerequisites

- **Operating System**: Linux (tested on CentOS/RHEL 9, Ubuntu 20.04+)
- **Python**: 3.10+
- **GPU**: NVIDIA GPU with CUDA support (required for basecalling and charging classification)
- **Storage**: ~50GB for tools, models, and test data

## Install Pixi

The pipeline uses [Pixi](https://pixi.sh) for environment management. Install it with:

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

After installation, restart your shell or run:

```bash
source ~/.bashrc  # or ~/.zshrc
```

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

The pipeline requires Dorado (ONT basecaller) and Modkit (modification toolkit). Install them with:

```bash
pixi run setup-tools
```

This downloads and installs:

- **Dorado** v0.9.1 - Oxford Nanopore basecaller
- **Modkit** v0.4.3 - Modification calling toolkit

Tools are installed to `resources/tools/` and automatically added to PATH when running the pipeline.

## Download Basecalling Model

Download the Dorado basecalling model:

```bash
pixi run snakemake dorado_model --cores 1
```

This downloads `rna004_130bps_sup@v5.1.0` to `resources/models/`.

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
resources/tools/dorado/0.9.1/bin/dorado --version

# Check Modkit installation
resources/tools/modkit/0.4.3/bin/modkit --version

# Dry run with test config
pixi run dry-run
```

## Directory Structure After Installation

```
aa-tRNA-seq-pipeline/
├── .pixi/                    # Pixi environment
├── resources/
│   ├── tools/
│   │   ├── dorado/0.9.1/    # Dorado binaries
│   │   └── modkit/0.4.3/    # Modkit binaries
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
pixi run setup-tools
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

### Modkit Build Fails

Modkit is built from source and requires Rust. If installation fails:

```bash
# Install Rust manually
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source ~/.cargo/env

# Retry setup
pixi run setup-tools
```

## Next Steps

- [Quick Start](quickstart.md) - Run the test pipeline
- [First Analysis](first-analysis.md) - Analyze your own data
