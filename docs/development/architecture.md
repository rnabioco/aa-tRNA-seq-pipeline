# Architecture

Technical architecture of the aa-tRNA-seq pipeline.

## Directory Structure

```
aa-tRNA-seq-pipeline/
├── workflow/
│   ├── Snakefile              # Main entry point
│   ├── rules/                 # Modular rule files
│   │   ├── common.smk         # Utilities and helpers
│   │   ├── aatrnaseq-process.smk    # Core processing
│   │   ├── aatrnaseq-charging.smk   # Charging analysis
│   │   ├── aatrnaseq-qc.smk         # Quality control
│   │   ├── aatrnaseq-modifications.smk  # Modkit rules
│   │   └── warpdemux.smk      # Demultiplexing (conditional)
│   └── scripts/               # Python processing scripts
├── config/
│   ├── config-base.yml        # Default configuration
│   ├── config-test.yml        # Test configuration
│   └── samples-*.tsv/yml      # Sample definitions
├── cluster/
│   ├── lsf/                   # LSF cluster profile
│   └── generic/               # Generic cluster profile
├── resources/
│   ├── ref/                   # Reference sequences
│   ├── models/                # ML models
│   ├── kmers/                 # Kmer tables
│   └── tools/                 # External tools (dorado, modkit)
├── .tests/                    # Test data and scripts
├── pixi.toml                  # Pixi environment definition
└── pixi.lock                  # Locked dependencies
```

## Snakemake Organization

### Entry Point

`workflow/Snakefile` is the main entry point:

```python
# Load base configuration
configfile: "config/config-base.yml"

# Dynamic PATH setup
onstart:
    # Add dorado and modkit to PATH
    dorado_path = f"resources/tools/dorado/{config['dorado_version']}/bin"
    os.environ["PATH"] = f"{dorado_path}:{os.environ['PATH']}"

# Include rule modules
include: "rules/common.smk"
include: "rules/aatrnaseq-process.smk"
include: "rules/aatrnaseq-charging.smk"
include: "rules/aatrnaseq-qc.smk"
include: "rules/aatrnaseq-modifications.smk"

# Conditionally include demux rules
if is_demux_enabled():
    include: "rules/warpdemux.smk"

# Target rule
rule all:
    input:
        pipeline_outputs()
```

### Rule Modules

```mermaid
flowchart TB
    A[Snakefile] --> B[common.smk<br/>Utilities]
    A --> C[aatrnaseq-process.smk<br/>7 rules]
    A --> D[aatrnaseq-charging.smk<br/>2 rules]
    A --> E[aatrnaseq-qc.smk<br/>3 rules]
    A --> F[aatrnaseq-modifications.smk<br/>4 rules]
    A -.->|conditional| G[warpdemux.smk<br/>5 rules]
```

### Rule Categories

| Module | Purpose | Rules |
|--------|---------|-------|
| `common.smk` | Utilities, sample parsing, output definition | - |
| `aatrnaseq-process.smk` | Core pipeline | 7 |
| `aatrnaseq-charging.smk` | Charging analysis | 2 |
| `aatrnaseq-qc.smk` | Quality control | 3 |
| `aatrnaseq-modifications.smk` | Modification calling | 4 |
| `warpdemux.smk` | Demultiplexing | 5 |

## Key Functions in common.smk

### Sample Management

```python
def parse_samples(fl):
    """Auto-detect and parse sample file (TSV or YAML)."""

def parse_samples_tsv(fl):
    """Parse TSV format: sample_id<TAB>path"""

def parse_samples_yaml(fl):
    """Parse YAML format with barcode assignments."""

def find_raw_inputs():
    """Recursively find POD5 files in run directories."""
```

### Configuration Helpers

```python
def get_modkit_threshold_opts():
    """Build modkit threshold options from config."""

def get_pipeline_commit():
    """Get git commit ID for reproducibility."""

def is_demux_enabled():
    """Check if WarpDemuX is enabled in config."""
```

### Output Definition

```python
def pipeline_outputs():
    """Define all final outputs for rule all."""
    outputs = []
    for sample in samples:
        outputs.extend([
            f"{outdir}/summary/tables/{sample}/{sample}.charging.cpm.tsv.gz",
            f"{outdir}/summary/tables/{sample}/{sample}.charging_prob.tsv.gz",
            f"{outdir}/summary/tables/{sample}/{sample}.bcerror.tsv.gz",
            # ... more outputs
        ])
    return outputs
```

## Data Flow

### Input Discovery

```mermaid
flowchart LR
    A[samples.tsv/yml] --> B[parse_samples]
    B --> C[find_raw_inputs]
    C --> D[samples dict<br/>path, barcode, run_id]
    D --> E[get_raw_inputs<br/>per sample]
```

### Processing Pipeline

```mermaid
flowchart TD
    subgraph Per Sample
        A[POD5 files] --> B[stage_pod5]
        B --> C[rebasecall]
        C --> E[bwa_align<br/>tags carried through]
        E --> F[classify_charging]
        F --> G[add_adapter_tags<br/>finalize_bam]
    end

    subgraph Summaries
        G --> H[get_cca_trna]
        G --> I[base_calling_error]
        G --> J[modkit_pileup]
    end
```

## Global Variables

Defined in `common.smk`:

```python
# Script directory path
SCRIPT_DIR = os.path.join(SNAKEFILE_DIR, "scripts")

# Output directory from config
outdir = config.get("output_directory", "results")

# Parsed samples dictionary
samples = parse_samples(config["samples"])
find_raw_inputs()  # Populates samples with POD5 paths
```

## Configuration System

### Hierarchy

```mermaid
flowchart TB
    A[config-base.yml<br/>Defaults] --> B[Your config.yml<br/>Overrides]
    B --> C[Command line<br/>--config key=value]
    C --> D[Final config dict]
```

### Key Config Sections

```yaml
# Tool paths and versions
base_calling_model: "resources/models/..."
dorado_version: 2.1.1

# Reference files
fasta: "resources/ref/..."
charging:
    model: "resources/models/charging/..."

# Command options
opts:
    dorado: "..."
    bwa: "..."
    bam_filter: "..."

# Modkit thresholds
modkit:
    filter_threshold: 0.5
    mod_thresholds: {...}

# Optional demux
warpdemux:
    enabled: false
```

## External Tools

### Tool Management

Tools are installed to `resources/tools/`:

```
resources/tools/
├── dorado/
│   └── <version>/
│       └── bin/dorado
```

### PATH Setup

Tools are added to PATH dynamically in `Snakefile`:

```python
onstart:
    dorado_path = f"resources/tools/dorado/{config['dorado_version']}/bin"
    os.environ["PATH"] = f"{dorado_path}:{os.environ['PATH']}"
```

## Python Scripts

Located in `workflow/scripts/`:

| Script | Called By | Purpose |
|--------|-----------|---------|
| `stamp_read_groups.py` | `bwa_align` | @RG/@CO header lines for `bwa mem -H` |
| `get_charging_table.py` | `get_cca_trna` | Extract CL tag |
| `get_trna_charging_cpm.py` | `get_cca_trna_cpm` | Calculate CPM |
| `get_bcerror_freqs.py` | `base_calling_error` | Error metrics |
| `get_align_stats.py` | `align_stats` | Read statistics |

Scripts are called via `SCRIPT_DIR`:

```python
params:
    src=SCRIPT_DIR,
shell:
    "python {params.src}/script.py ..."
```

## Reproducibility

### Git Commit Tracking

```python
def get_pipeline_commit():
    """Return current git commit ID."""
    try:
        repo = Repo(SNAKEFILE_DIR, search_parent_directories=True)
        return repo.head.commit.hexsha[:8]
    except:
        return "unknown"
```

### Locked Dependencies

`pixi.lock` ensures reproducible environment:

```bash
# Install exact versions
pixi install

# Update lock file
pixi update
```

## Next Steps

- [Adding Rules](adding-rules.md) - Add new functionality
- [Contributing](contributing.md) - Contribution guidelines
