"""
Rules for WarpDemuX barcode demultiplexing and EDX (3' adapter barcode) concordance.
Only loaded when warpdemux.enabled is true in config.
"""

import pandas as pd
from pathlib import Path
import gzip

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


rule warpdemux:
    """
Run WarpDemuX barcode demultiplexing directly on raw POD5 files.
"""
    input:
        get_run_raw_inputs,
    output:
        outdir=maybe_temp(
            directory(os.path.join(outdir, "demux", "warpdemux_output", "{run_id}")),
            tier="demux_scratch",
        ),
        done=maybe_temp(
            os.path.join(outdir, "demux", "warpdemux_output", "{run_id}", ".done"),
            tier="demux_scratch",
        ),
    log:
        os.path.join(outdir, "logs", "warpdemux", "{run_id}"),
    threads: config.get("warpdemux", {}).get("threads", 16)
    params:
        model=lambda wildcards: get_barcode_kit_for_run(wildcards.run_id),
        save_boundaries=lambda wildcards: (
            "true"
            if config.get("warpdemux", {}).get("save_boundaries", True)
            else "false"
        ),
        pod5_dirs=lambda wildcards: " ".join(get_run_pod5_dirs(wildcards.run_id)),
    shell:
        """
        warpdemux demux \
            -i {params.pod5_dirs} \
            -o {output.outdir} \
            -m {params.model} \
            -j {threads} \
            --save_boundaries {params.save_boundaries} 2>&1 | tee {log}
        touch {output.done}
        """


rule parse_warpdemux:
    """
Parse WarpDemuX predictions and create a barcode mapping file per run.
"""
    input:
        demux_done=os.path.join(
            outdir, "demux", "warpdemux_output", "{run_id}", ".done"
        ),
        demux_dir=os.path.join(outdir, "demux", "warpdemux_output", "{run_id}"),
    output:
        mapping=maybe_temp(
            os.path.join(
                outdir, "demux", "read_ids", "{run_id}", "barcode_mapping.tsv.gz"
            ),
            tier="demux_scratch",
        ),
        summary=maybe_temp(
            os.path.join(
                outdir, "demux", "read_ids", "{run_id}", "demux_summary.tsv.gz"
            ),
            tier="demux_scratch",
        ),
    log:
        os.path.join(outdir, "logs", "parse_warpdemux", "{run_id}"),
    run:
        # WarpDemuX creates a timestamped subdirectory, find it
        demux_base = Path(input.demux_dir)
        subdirs = [
            d
            for d in demux_base.iterdir()
            if d.is_dir() and d.name.startswith("warpdemux_")
        ]
        if subdirs:
            pred_dir = subdirs[0] / "predictions"
        else:
            pred_dir = demux_base / "predictions"
        all_predictions = []
        for pred_file in pred_dir.glob("*.csv.gz"):
            with gzip.open(pred_file, "rt") as f:
                df = pd.read_csv(f, comment=None)
                # Handle #read_id column name
                df.columns = [c.lstrip("#") for c in df.columns]
                all_predictions.append(df)
        predictions = (
            pd.concat(all_predictions, ignore_index=True)
            if all_predictions
            else pd.DataFrame(columns=["read_id", "predicted_barcode"])
        )
        # Convert numeric barcode (e.g., 7) to barcode format (e.g., "barcode07")
        predictions["predicted_barcode"] = predictions["predicted_barcode"].apply(
            lambda x: f"barcode{int(x):02d}" if x != -1 else "unclassified"
        )
        # Write the full mapping file
        predictions[["read_id", "predicted_barcode"]].to_csv(
            output.mapping, sep="\t", index=False, compression="gzip"
        )
        # Write summary statistics
        summary_data = (
            predictions.groupby("predicted_barcode")
            .size()
            .reset_index(name="n_reads")
        )
        summary_data.to_csv(
            output.summary, sep="\t", index=False, compression="gzip"
        )


def get_sample_barcode_mapping(wildcards):
    """Get the barcode mapping file for a sample's run."""
    run_id = samples[wildcards.sample]["run_id"]
    return os.path.join(outdir, "demux", "read_ids", run_id, "barcode_mapping.tsv.gz")


rule extract_sample_reads:
    """
Extract read IDs for a specific sample based on its barcode assignment.
"""
    input:
        mapping=get_sample_barcode_mapping,
    output:
        read_ids=maybe_temp(
            os.path.join(outdir, "demux", "read_ids", "{sample}", "{sample}.txt"),
            tier="demux_scratch",
        ),
    log:
        os.path.join(outdir, "logs", "extract_sample_reads", "{sample}"),
    params:
        barcode=lambda wildcards: samples[wildcards.sample]["barcode"],
    run:
        # Read the mapping file and filter for this sample's barcode
        mapping = pd.read_csv(input.mapping, sep="\t", compression="gzip")
        sample_reads = mapping[mapping["predicted_barcode"] == params.barcode][
            "read_id"
        ]
        with open(output.read_ids, "w") as f:
            for read_id in sample_reads:
                f.write(f"{read_id}\n")


def get_sample_read_ids(wildcards):
    """Get the read ID file for a sample."""
    return os.path.join(
        outdir, "demux", "read_ids", wildcards.sample, f"{wildcards.sample}.txt"
    )


def get_sample_pod5_dirs(wildcards):
    """Get POD5 directories for a sample's run (for shell commands)."""
    run_id = samples[wildcards.sample]["run_id"]
    return " ".join(get_run_pod5_dirs(run_id))


rule split_pod5:
    """
Filter raw POD5 files by sample using read IDs from demultiplexing.
"""
    input:
        pod5=get_sample_run_raw_inputs,
        read_ids=get_sample_read_ids,
    output:
        maybe_temp(
            os.path.join(outdir, "demux", "pod5", "{sample}", "{sample}.pod5"),
            tier="split_pod5",
        ),
    log:
        os.path.join(outdir, "logs", "split_pod5", "{sample}"),
    params:
        pod5_dirs=get_sample_pod5_dirs,
    shell:
        """
        pod5 filter {params.pod5_dirs} --ids {input.read_ids} --missing-ok --output {output} 2>&1 | tee {log}
        """


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
        pod5=os.path.join(
            outdir, "demux", "edx", "pod5", "{sample}", "{sample}.pod5"
        ),
    log:
        os.path.join(outdir, "logs", "filter_pod5_by_edx", "{sample}"),
    shell:
        """
        pod5 filter {input.pod5} --ids {input.read_ids} --missing-ok --output {output.pod5} 2>&1 | tee {log}
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
