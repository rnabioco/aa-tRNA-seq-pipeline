"""
Barcode demultiplexing with ``escpod demux`` (Rust).

Included from demux.smk when ``warpdemux.backend`` is ``escpod``. Provides the
same downstream contract as the warpdemux backend:

  demux/read_ids/{run_id}/barcode_mapping.tsv.gz   per-read barcode assignment
  demux/pod5/{sample}/{sample}.pod5                per-sample POD5

The fused ``escpod demux`` pipeline (detect -> fingerprint -> classify -> split)
makes a single pass over the raw POD5s and emits both the classifications table
and one POD5 per barcode, so the separate read-ID extraction and POD5 filter
pass that the warpdemux backend needs are not required here.
"""


def get_escpod_demux_dir(run_id):
    return os.path.join(outdir, "demux", "escpod_output", run_id)


def get_escpod_barcode_model():
    model = config.get("warpdemux", {}).get("barcode_model")
    if not model:
        sys.exit(
            "warpdemux.barcode_model is not set. Run "
            "'bash scripts/install-demux-models.sh' and check config."
        )
    return model


def get_escpod_min_confidence():
    """Minimum classifier confidence for a read to keep its barcode call.

    The shipped barcode GBM carries no per-class thresholds, so escpod assigns
    every read with a usable adapter boundary to some barcode -- on a test run
    of 210 non-barcoded reads it produced zero `unclassified` calls at a median
    confidence of 0.57. WarpDemuX instead rejected low-confidence reads. This
    threshold restores that behavior; 0.0 keeps escpod's native output.
    """
    value = config.get("warpdemux", {}).get("min_confidence", 0.0)
    if not 0.0 <= float(value) <= 1.0:
        sys.exit(f"warpdemux.min_confidence must be between 0 and 1, got {value}.")
    return float(value)


def get_escpod_detect_opts():
    """Boundary-detection flags for the fused demux pipeline.

    The barcode GBM was trained behind the CNN adapter-boundary detector, so
    'cnn' is the parity setting. 'llr' needs no ONNX model but agrees with the
    CNN only ~82% within +/-200 samples, which shifts calls.
    """
    wdx = config.get("warpdemux", {})
    method = wdx.get("method", "cnn")
    if method == "llr":
        return "--method llr"
    if method != "cnn":
        sys.exit(f"Unknown warpdemux.method '{method}'. Use 'cnn' or 'llr'.")
    adapter_model = wdx.get("adapter_model")
    if not adapter_model:
        sys.exit(
            "warpdemux.method is 'cnn' but warpdemux.adapter_model is not set. "
            "Run 'bash scripts/install-demux-models.sh', or set "
            "warpdemux.method to 'llr'."
        )
    return f"--method cnn --cnn-model {adapter_model}"


rule escpod_demux:
    """
    Classify barcodes and split POD5s in one pass over the raw signal.
    """
    input:
        pod5=get_run_raw_inputs,
        barcode_model=get_escpod_barcode_model(),
    output:
        outdir=maybe_temp(
            directory(os.path.join(outdir, "demux", "escpod_output", "{run_id}")),
            tier="split_pod5",
        ),
        # Deliberately outside the directory() output above: a file output nested
        # inside a directory output is removed along with the directory when
        # Snakemake cleans up before a rerun.
        classifications=maybe_temp(
            os.path.join(
                outdir, "demux", "read_ids", "{run_id}", "classifications.csv"
            ),
            tier="demux_scratch",
        ),
    log:
        os.path.join(outdir, "logs", "escpod_demux", "{run_id}"),
    threads: config.get("warpdemux", {}).get("threads", 16)
    params:
        detect_opts=get_escpod_detect_opts(),
    shell:
        """
        mkdir -p {output.outdir} $(dirname {output.classifications})

        escpod demux {input.pod5} \
            --model {input.barcode_model} \
            {params.detect_opts} \
            --output-dir {output.outdir} \
            --classifications {output.classifications} \
            --prefix barcode \
            --threads {threads} 2>&1 | tee {log}
        """


rule parse_escpod_demux:
    """
    Convert the escpod classifications CSV to the pipeline's barcode mapping table.

    escpod labels barcodes ``BC{id:02}``; the pipeline and its samples files use
    ``barcode{id:02}``. Reads with no usable adapter boundary are ``unclassified``.
    """
    input:
        classifications=rules.escpod_demux.output.classifications,
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
        os.path.join(outdir, "logs", "parse_escpod_demux", "{run_id}"),
    params:
        min_confidence=get_escpod_min_confidence(),
    run:
        predictions = pd.read_csv(input.classifications)
        predictions.columns = [c.lstrip("#") for c in predictions.columns]
        # `classify --model` names the column predicted_barcode; the fused
        # pipeline names it barcode.
        barcode_col = (
            "barcode" if "barcode" in predictions.columns else "predicted_barcode"
        )
        predictions["predicted_barcode"] = predictions[barcode_col].apply(
            lambda x: (
                "unclassified"
                if pd.isna(x) or str(x) == "unclassified"
                else (
                    "barcode" + str(x)[len("BC") :]
                    if str(x).startswith("BC")
                    else str(x)
                )
            )
        )
        # Reject low-confidence calls to `unclassified`, standing in for the
        # per-class thresholds the shipped GBM does not carry.
        if params.min_confidence > 0:
            if "confidence" not in predictions.columns:
                sys.exit(
                    "warpdemux.min_confidence is set but the escpod "
                    "classifications table has no `confidence` column."
                )
            below = predictions["confidence"] < params.min_confidence
            logger.info(
                f"{wildcards.run_id}: rejecting {int(below.sum())} of "
                f"{len(predictions)} reads below confidence "
                f"{params.min_confidence}"
            )
            predictions.loc[below, "predicted_barcode"] = "unclassified"
        predictions[["read_id", "predicted_barcode"]].to_csv(
            output.mapping, sep="\t", index=False, compression="gzip"
        )
        if "confidence" in predictions.columns:
            summary_data = (
                predictions.groupby("predicted_barcode")
                .agg(
                    n_reads=("read_id", "size"),
                    mean_confidence=("confidence", "mean"),
                )
                .reset_index()
            )
        else:
            summary_data = (
                predictions.groupby("predicted_barcode")
                .size()
                .reset_index(name="n_reads")
            )
        summary_data.to_csv(
            output.summary, sep="\t", index=False, compression="gzip"
        )


def get_escpod_sample_pod5(wildcards):
    """The per-barcode POD5 that escpod demux wrote for this sample."""
    run_id = samples[wildcards.sample]["run_id"]
    label = wdx_to_escpod_barcode(samples[wildcards.sample]["barcode"])
    return os.path.join(get_escpod_demux_dir(run_id), f"barcode_{label}.pod5")


def get_escpod_run_dir(wildcards):
    """The escpod demux output directory for this sample's run."""
    return get_escpod_demux_dir(samples[wildcards.sample]["run_id"])


# How the per-sample POD5 is produced depends on whether a confidence threshold
# is in force. Without one, escpod's per-barcode POD5 is already exactly the
# sample's read set and can be linked. With one, the low-confidence reads have to
# be removed from it, which needs a filter pass over the (small) per-barcode file.
if get_escpod_min_confidence() > 0:

    rule extract_escpod_sample_reads:
        """
        Read IDs for this sample after confidence thresholding.
        """
        input:
            mapping=get_sample_barcode_mapping,
        output:
            read_ids=maybe_temp(
                os.path.join(outdir, "demux", "read_ids", "{sample}", "{sample}.txt"),
                tier="demux_scratch",
            ),
        params:
            barcode=lambda wildcards: samples[wildcards.sample]["barcode"],
        run:
            mapping = pd.read_csv(input.mapping, sep="\t", compression="gzip")
            sample_reads = mapping[mapping["predicted_barcode"] == params.barcode][
                "read_id"
            ]
            with open(output.read_ids, "w") as f:
                for read_id in sample_reads:
                    f.write(f"{read_id}\n")

    rule collect_escpod_pod5:
        """
        Drop below-threshold reads from the escpod per-barcode POD5.
        """
        input:
            demux_dir=get_escpod_run_dir,
            read_ids=rules.extract_escpod_sample_reads.output.read_ids,
        output:
            maybe_temp(
                os.path.join(outdir, "demux", "pod5", "{sample}", "{sample}.pod5"),
                tier="split_pod5",
            ),
        log:
            os.path.join(outdir, "logs", "collect_escpod_pod5", "{sample}"),
        threads: 4
        params:
            src=get_escpod_sample_pod5,
        shell:
            """
            escpod filter {params.src} \
                --ids {input.read_ids} \
                --threads {threads} \
                --force \
                --output {output} 2>&1 | tee {log}
            """

else:

    rule collect_escpod_pod5:
        """
        Expose the escpod per-barcode POD5 under the per-sample path.

        escpod demux already wrote one POD5 per barcode, so this is a symlink
        rather than a second filter pass. The link and its target share the
        split_pod5 cleanup tier so they are removed together.
        """
        input:
            demux_dir=get_escpod_run_dir,
        output:
            maybe_temp(
                os.path.join(outdir, "demux", "pod5", "{sample}", "{sample}.pod5"),
                tier="split_pod5",
            ),
        log:
            os.path.join(outdir, "logs", "collect_escpod_pod5", "{sample}"),
        params:
            src=get_escpod_sample_pod5,
        shell:
            """
            if [ ! -f "{params.src}" ]; then
                echo "Error: {params.src} not found." >&2
                echo "escpod demux emits one POD5 per class in its barcode" >&2
                echo "model; this sample's barcode is not among them. Check" >&2
                echo "that the samples file barcodes match the model in" >&2
                echo "warpdemux.barcode_model." >&2
                exit 1
            fi
            ln -sf "$(realpath {params.src})" {output} 2>&1 | tee {log}
            """
