"""
Rules for barcode demultiplexing and EDX (3' adapter barcode) concordance.
Only loaded when warpdemux.enabled is true in config.

Two backends produce the same downstream contract — a per-run
``barcode_mapping.tsv.gz`` and a per-sample ``demux/pod5/{sample}/{sample}.pod5``:

  escpod    ``escpod demux`` (Rust). One pass over the raw POD5s emits both the
            per-read classifications and the per-barcode POD5 files, so no
            separate filter pass is needed.
  warpdemux The original python implementation, followed by an ``escpod filter``
            pass to split the POD5s by read ID.

See ``warpdemux.backend`` in config/config-base.yml for why the two do not
produce identical calls.
"""

import pandas as pd
from pathlib import Path
import gzip


def get_demux_backend():
    """Demultiplexing backend: 'escpod' (default) or 'warpdemux'."""
    backend = config.get("warpdemux", {}).get("backend", "escpod")
    if backend not in ("escpod", "warpdemux"):
        sys.exit(f"Unknown warpdemux.backend '{backend}'. Use 'escpod' or 'warpdemux'.")
    return backend


def wdx_to_escpod_barcode(barcode):
    """Map a samples-file barcode name to the label escpod emits.

    escpod formats the integer barcode id from the model's label_mapper as
    ``BC{id:02}``, where WarpDemuX and this pipeline use ``barcode{id:02}``. The
    digits are the same, so this is a prefix rename.
    """
    if barcode is None:
        return None
    if barcode.startswith("barcode"):
        return "BC" + barcode[len("barcode") :]
    if barcode.startswith("WDX_bc"):
        return "BC" + barcode[len("WDX_bc") :]
    if barcode.startswith("BC"):
        return barcode
    sys.exit(
        f"Cannot map barcode '{barcode}' to an escpod label. Expected a name "
        "like 'barcode03', 'WDX_bc03', or 'BC03'."
    )


# WarpDemuX helper functions


def get_run_ids():
    """Get unique run_ids for samples that require demultiplexing."""
    run_ids = set()
    for sample, info in samples.items():
        if info.get("barcode") and info.get("run_id"):
            run_ids.add(info["run_id"])
    return list(run_ids)


def get_run_path(run_id):
    """Get the path for a run_id."""
    for sample, info in samples.items():
        if info.get("run_id") == run_id:
            # Return first path from the set
            return list(info["path"])[0]
    return None


def get_samples_for_run(run_id):
    """Get list of sample names assigned to a run."""
    return [
        sample
        for sample, info in samples.items()
        if info.get("run_id") == run_id and info.get("barcode")
    ]


def get_barcodes_for_run(run_id):
    """Get barcode→sample mapping for a run."""
    return {
        info["barcode"]: sample
        for sample, info in samples.items()
        if info.get("run_id") == run_id and info.get("barcode")
    }


def get_barcode_kit_for_run(run_id):
    """Get the barcode kit for a run."""
    for sample, info in samples.items():
        if info.get("run_id") == run_id:
            return info.get("barcode_kit")
    return config.get("warpdemux", {}).get("barcode_kit")


def get_run_raw_inputs(wildcards):
    """Get all POD5 files for a run directory (for Snakemake DAG tracking)."""
    POD5_DIRS = ["pod5_pass", "pod5_fail", "pod5"]
    ext = ".pod5"
    run_path = get_run_path(wildcards.run_id)

    raw_fls = []
    for subdir in POD5_DIRS:
        data_path = os.path.join(run_path, subdir, "*" + ext)
        fls = glob.glob(data_path)
        raw_fls += fls

    if len(raw_fls) == 0:
        sys.exit(
            f"No input files found for run: {wildcards.run_id} at {run_path}. "
            "Please check the path in the samples file"
        )
    return raw_fls


def get_run_pod5_dirs(run_id):
    """Get existing POD5 subdirectories for a run (for shell commands)."""
    POD5_DIRS = ["pod5_pass", "pod5_fail", "pod5"]
    run_path = get_run_path(run_id)
    return [
        os.path.join(run_path, d)
        for d in POD5_DIRS
        if os.path.isdir(os.path.join(run_path, d))
    ]


def get_sample_run_raw_inputs(wildcards):
    """Get all raw POD5 files for a sample's run (for Snakemake DAG tracking)."""
    run_id = samples[wildcards.sample]["run_id"]
    POD5_DIRS = ["pod5_pass", "pod5_fail", "pod5"]
    ext = ".pod5"
    run_path = get_run_path(run_id)

    raw_fls = []
    for subdir in POD5_DIRS:
        data_path = os.path.join(run_path, subdir, "*" + ext)
        fls = glob.glob(data_path)
        raw_fls += fls

    if len(raw_fls) == 0:
        sys.exit(
            f"No input files found for run: {run_id} at {run_path}. "
            "Please check the path in the samples file"
        )
    return raw_fls


def get_sample_barcode_mapping(wildcards):
    """Get the barcode mapping file for a sample's run."""
    run_id = samples[wildcards.sample]["run_id"]
    return os.path.join(outdir, "demux", "read_ids", run_id, "barcode_mapping.tsv.gz")


# --- Backend-specific demultiplexing rules ---
#
# Both backends produce demux/read_ids/{run_id}/barcode_mapping.tsv.gz and
# demux/pod5/{sample}/{sample}.pod5, so everything below and downstream is
# backend-agnostic.


if get_demux_backend() == "escpod":

    include: "demux-escpod.smk"

else:

    include: "demux-warpdemux.smk"


# --- EDX (3' adapter barcode) early demultiplexing ---


def get_edx_samples():
    """Get sample names that have EDX adapter assignments."""
    return [s for s, info in samples.items() if info.get("edx")]


rule detect_edx_adapters:
    """
    Detect 3' adapter identity per read on the unaligned BAM (before alignment).
    Produces a TSV mapping each read_id to its best-matching 3' adapter name.
    """
    input:
        bam=lambda wildcards: os.path.join(
            outdir,
            "bam",
            "rebasecall",
            wildcards.sample,
            f"{wildcards.sample}.rbc.bam",
        ),
    output:
        tsv=os.path.join(
            outdir, "demux", "edx", "{sample}", "{sample}.edx_adapters.tsv.gz"
        ),
    log:
        os.path.join(outdir, "logs", "detect_edx_adapters", "{sample}"),
    params:
        src=SCRIPT_DIR,
        adapter_3p_args=lambda wc: " ".join(
            f'--adapter-3p "{name}:{seq}"' for name, seq in get_adapter_3p_list()
        ),
        min_score_3p=config["adapters"]["min_score_3p"],
    shell:
        """
        python {params.src}/detect_3p_adapters.py \
            --bam {input.bam} \
            {params.adapter_3p_args} \
            --min-score-3p {params.min_score_3p} \
            --output {output.tsv} \
            2>&1 | tee {log}
        """


rule extract_edx_read_ids:
    """
    Extract read IDs matching this sample's EDX adapter assignment.
    """
    input:
        tsv=rules.detect_edx_adapters.output.tsv,
    output:
        read_ids=maybe_temp(
            os.path.join(
                outdir, "demux", "edx", "{sample}", "{sample}.edx_read_ids.txt"
            ),
            tier="demux_scratch",
        ),
    params:
        edx_adapter_name=get_sample_edx,
    run:
        import gzip

        with (
            gzip.open(input.tsv, "rt") as f_in,
            open(output.read_ids, "w") as f_out,
        ):
            header = f_in.readline()  # skip header
            for line in f_in:
                read_id, adapter = line.rstrip("\n").split("\t", 1)
                if adapter == params.edx_adapter_name:
                    f_out.write(f"{read_id}\n")


rule filter_fastq_by_edx:
    """
    Extract FASTQ for reads matching this sample's EDX adapter.
    """
    input:
        bam=lambda wildcards: os.path.join(
            outdir,
            "bam",
            "rebasecall",
            wildcards.sample,
            f"{wildcards.sample}.rbc.bam",
        ),
        read_ids=rules.extract_edx_read_ids.output.read_ids,
    output:
        fq=maybe_temp(
            os.path.join(outdir, "demux", "edx", "fq", "{sample}", "{sample}.fq.gz"),
            tier="fastq",
        ),
    log:
        os.path.join(outdir, "logs", "filter_fastq_by_edx", "{sample}"),
    shell:
        """
        samtools view -N {input.read_ids} {input.bam} \
            | samtools fastq - \
            | gzip >{output.fq} \
                2>&1 | tee {log}
        """


rule filter_pod5_by_edx:
    """
    Filter POD5 to keep only reads matching this sample's EDX adapter.
    """
    input:
        pod5=get_sample_pod5,
        read_ids=rules.extract_edx_read_ids.output.read_ids,
    output:
        pod5=os.path.join(outdir, "demux", "edx", "pod5", "{sample}", "{sample}.pod5"),
    log:
        os.path.join(outdir, "logs", "filter_pod5_by_edx", "{sample}"),
    threads: 8
    shell:
        """
        escpod filter {input.pod5} \
            --ids {input.read_ids} \
            --threads {threads} \
            --force \
            --output {output.pod5} 2>&1 | tee {log}
        """


# --- EDX (3' adapter barcode) concordance analysis ---


rule edx_concordance:
    """
    Build concordance table of WDX sample assignment vs EDX adapter identity.
    Uses pre-alignment adapter detection TSVs (which contain ALL reads with their
    detected adapter) rather than final BAMs (which only contain matching reads).
    """
    input:
        tsvs=lambda wildcards: expand(
            os.path.join(
                outdir, "demux", "edx", "{sample}", "{sample}.edx_adapters.tsv.gz"
            ),
            sample=get_edx_samples(),
        ),
    output:
        concordance=os.path.join(outdir, "summary", "edx", "edx_concordance.tsv.gz"),
    log:
        os.path.join(outdir, "logs", "edx", "edx_concordance.log"),
    params:
        src=SCRIPT_DIR,
        sample_names=lambda wildcards: " ".join(get_edx_samples()),
    shell:
        """
        python {params.src}/edx_concordance.py \
            --tsvs {input.tsvs} \
            --samples {params.sample_names} \
            --output {output.concordance} \
            2>&1 | tee {log}
        """
