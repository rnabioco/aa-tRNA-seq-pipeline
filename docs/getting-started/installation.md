# Installation

This guide covers installing the aa-tRNA-seq pipeline and its dependencies.

## Prerequisites

- **Operating System**: Linux (tested on CentOS/RHEL 9, Ubuntu 20.04+)
- **Python**: 3.10+ (3.12+ required if you use the `leech` classifier)
- **GPU**: NVIDIA GPU with CUDA support (required for basecalling and charging classification)
- **Storage**: ~50GB for tools, models, and test data
- **Rust toolchain**: >= 1.95, required by `pixi run setup` to build the `escpod` CLI
  from source (see [escpod](#escpod-pod5-io-and-demultiplexing) below)
- **GitHub CLI (`gh`)**, authenticated, to install `leech` from its private release
  wheels — or a local wheel directory via `LEECH_WHEEL_DIR`

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

### Install Rust (for `escpod`)

If `cargo` is not already available:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

## Clone the Repository

```bash
git clone https://github.com/rnabioco/aa-tRNA-seq-pipeline.git
cd aa-tRNA-seq-pipeline
```

!!! note "No submodules"
    Earlier versions of the pipeline vendored `leech` as a `resources/leech` git
    submodule. That submodule (and `.gitmodules`) has been removed — leech is now
    installed from release wheels, so a plain `git clone` is all you need. Ignore any
    older instructions to run `git submodule update --init --recursive`.

## Install Dependencies

Install all Python dependencies via Pixi:

```bash
pixi install
```

This creates a `.pixi` directory with all required packages including:

- Snakemake 8.0+
- pysam
- pandas
- modkit
- samtools
- bwa
- deeptools

## Install External Tools

The pipeline requires several external tools. Install them with a single command:

```bash
pixi run setup
```

This downloads and installs:

- **Dorado** - Oxford Nanopore basecaller (version set in `config/config-base.yml`)
- **Dorado model** - `rna004_130bps_sup@v5.3.0` basecalling model, plus modification models
- **Remora** - ONT signal analysis for charging classification
- **escpod** - the escapepod-rs CLI, **built from source** (POD5 merge/filter and demultiplexing)
- **escapepod** - the matching python package for POD5 I/O
- **leech** - GPU charging / amino-acid classifiers, from private GitHub release wheels
- **escpod demux models** - barcode + adapter-boundary models (best effort; see below)
- **WarpDemuX** - the original python demultiplexer (alternative demux backend)

Dorado, escpod and leech wheels are installed under `resources/tools/`, models under `resources/models/`. Modkit is managed by pixi (installed via conda).

Individual components can be (re)installed on their own:

```bash
pixi run install-escpod          # build the escpod CLI with the demux feature
pixi run install-leech           # install leech + leech-core from release wheels
pixi run install-demux-models    # install the escpod barcode/adapter models
pixi run install-warpdemux       # clone and install the python WarpDemuX
```

### escpod (POD5 I/O and demultiplexing)

All POD5 manipulation uses `escpod` rather than the ONT `pod5` CLI — it is 3-9x faster
on merge/filter and writes crash-safely (output is staged to a temp file and renamed,
so an interrupted run never leaves a corrupt POD5).

`escpod` **must** be built from source: the published escapepod-rs release binaries are
compiled with default features only and do not contain a working `escpod demux`. The
build requires Rust >= 1.95 and uses `--features cnn-detect` (which implies `demux` and
adds `--method cnn`). The binary lands in
`resources/tools/escapepod/<version>/bin`, which `pixi shell` puts on `PATH`.

The version is pinned by `escapepod_version` in `config/config-base.yml`, which also
pins the `escapepod` python package installed from PyPI. The ONT `pod5` python package
is still installed, but only for the optional remora signal-metrics QC script
(`workflow/scripts/extract_signal_metrics.py`), which needs remora's `io` API and so
cannot use escapepod's reader. That path only runs when `remora_kmer_table` is set
(default `null`).

### leech (optional classifier)

leech lives in the private `rnabioco/leech` repo and is not on PyPI, so
`pixi run setup` fetches its release wheels with `gh` — you need an authenticated
`gh` CLI with read access, or `LEECH_WHEEL_DIR` pointing at a local directory of
wheels. Versions are pinned by `leech_version` and `leech_core_version` in
`config/config-base.yml`.

`leech-core` is a separate, optional Rust extension wheel built per interpreter. If no
matching wheel is available, leech still installs but falls back to a slower
pure-python extraction backend and `--backend rust` is unavailable. leech requires
Python >= 3.12.

### escpod demux models

Only needed when demultiplexing with the default `escpod` backend:

```bash
pixi run install-demux-models
```

This installs `barcode_wdx4_rna004.gbm.json` and `adapter_rna004.onnx` into
`resources/models/demux/` (gitignored), verifying sha256 against the upstream
`MANIFEST.json`. Models are sourced from a local `rnabioco/escapepod-models` checkout
(`ESCAPEPOD_MODELS_DIR`, or a sibling directory of this repo), or failing that, a
GitHub release.

!!! warning "Barcode models are not released yet"
    The barcode GBM models are not yet published as releases on
    `rnabioco/escapepod-models` (only `adapter_rna004@v1.0.1` is), so the
    local-checkout path is currently required. The script prints the
    `scripts/release_model.sh` command needed to publish them. `pixi run setup`
    treats a failure here as a warning, since it only affects the escpod demux
    backend.

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
resources/tools/dorado/*/bin/dorado --version

# Check Modkit installation (managed by pixi)
pixi run modkit --version

# Check escpod, and that it was built with demux (needed for the escpod demux backend)
pixi run escpod --version
pixi run escpod demux --help

# Check leech (only if using classifier: leech or classify_aa)
pixi run leech --version

# Dry run with test config
pixi run dry-run
```

## Directory Structure After Installation

```
aa-tRNA-seq-pipeline/
├── .pixi/                    # Pixi environment (includes modkit, remora)
├── resources/
│   ├── tools/
│   │   ├── dorado/<version>/     # Dorado binaries
│   │   ├── escapepod/<version>/  # escpod CLI (source build)
│   │   ├── leech/<version>/      # downloaded leech wheels
│   │   └── WarpDemuX/            # WarpDemuX (warpdemux demux backend)
│   ├── models/
│   │   ├── rna004_130bps_sup@v5.3.0/  # Basecalling model
│   │   ├── cca_classifier.pt          # Remora charging model
│   │   └── demux/                     # escpod barcode + adapter models (gitignored)
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

### `cargo: not found` During Setup

The `escpod` build needs a Rust toolchain (>= 1.95):

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source "$HOME/.cargo/env"
pixi run install-escpod
```

### leech Wheels Cannot Be Downloaded

`gh release download` fails when `gh` is missing or unauthenticated, since
`rnabioco/leech` is private:

```bash
gh auth login
pixi run install-leech

# or install from a local wheel directory instead
LEECH_WHEEL_DIR=/path/to/wheels pixi run install-leech
```

leech is only needed for `classifier: leech` and the amino-acid classification rules;
the default Remora path works without it.

## Next Steps

- [Quick Start](quickstart.md) - Run the test pipeline
- [First Analysis](first-analysis.md) - Analyze your own data
