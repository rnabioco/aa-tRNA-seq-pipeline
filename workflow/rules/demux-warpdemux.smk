"""
Barcode demultiplexing with WarpDemuX (the original python implementation).

Included from demux.smk when ``warpdemux.backend`` is ``warpdemux``. Produces the
same downstream contract as the escpod backend:

  demux/read_ids/{run_id}/barcode_mapping.tsv.gz   per-read barcode assignment
  demux/pod5/{sample}/{sample}.pod5                per-sample POD5

Unlike the escpod backend, WarpDemuX only classifies — splitting the POD5s is a
second pass, done here with ``escpod filter`` over the run's read IDs.
"""


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
        # Not tiered for cleanup: this is the durable record of what the
        # demultiplexer called, and it is small.
        summary=os.path.join(
            outdir, "demux", "read_ids", "{run_id}", "demux_summary.tsv.gz"
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


def get_sample_run_dir(wildcards):
    """Get the run directory for a sample (for shell commands).

    ``escpod filter`` takes a single positional input, expanding a directory
    recursively, so the run root is passed rather than the individual
    pod5_pass/pod5_fail/pod5 subdirectories. Reads outside the ID list are
    dropped regardless, so a wider scan only costs walk time.
    """
    run_id = samples[wildcards.sample]["run_id"]
    return get_run_path(run_id)


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
    threads: 8
    params:
        run_dir=get_sample_run_dir,
    shell:
        """
        escpod filter {params.run_dir} \
            --ids {input.read_ids} \
            --threads {threads} \
            --force \
            --output {output} 2>&1 | tee {log}
        """
