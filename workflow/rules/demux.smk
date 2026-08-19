"""
Rules for signal-level barcode demultiplexing and EDX (3' adapter barcode)
concordance. Loaded when either demux backend is enabled in config.

Two backends converge on the same per-sample split POD5, by different routes:

  warpdemux : WarpDemuX DTW/fingerprint classifier, WDX barcodes ("barcode04").
              Classifies into a table, which is then parsed into a read->barcode
              mapping and applied with `pod5 filter` to cut the POD5 per sample.
  ldx       : escapepod CTC-CRF barcode basecaller, LDX barcodes ("nbc01").
              A single fused pass detects, basecalls, matches and routes each
              read into its barcode's POD5, so the split already exists when the
              command returns and only needs renaming barcode -> sample.

Only one may be enabled (see is_demux_enabled); a ruleorder picks which rule
supplies the split POD5.
"""

import pandas as pd
from pathlib import Path
from snakemake.exceptions import WorkflowError
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


# --- escapepod CTC-CRF (LDX / nbc) demultiplexing ---


def get_ort_dylib():
    """Absolute path to the CUDA-enabled libonnxruntime for `ldx.gpu`.

    Pinned to 1.27.x rather than whatever is newest: `ort 2.0.0-rc.13` enables
    `api-27`, so it refuses to load an older onnxruntime with a BadVersion
    error. conda-forge's CUDA builds stop at 1.26 today, which is why this is a
    vendored tarball from the onnxruntime GitHub release rather than a package.
    """
    root = os.path.join(PIPELINE_DIR, "resources", "tools", "onnxruntime")
    hits = sorted(glob.glob(os.path.join(root, "*", "lib", "libonnxruntime.so.1.27.*")))
    if not hits:
        sys.exit(
            "ldx.gpu is true but no CUDA libonnxruntime 1.27.x was found under "
            f"{root}.\nInstall it with `pixi run install-ort-gpu` (or set "
            "ldx.gpu: false to run demux on CPU)."
        )
    # cuDNN is a separate requirement: the onnxruntime tarball ships the CUDA
    # provider but not libcudnn.so.9, and without it the provider fails to
    # register and onnxruntime falls back to CPU with only a warning — the run
    # succeeds, just without a GPU, which is the worst failure mode here.
    cudnn = glob.glob(
        os.path.join(PIPELINE_DIR, ".pixi", "envs", "gpu", "lib", "libcudnn.so.9*")
    )
    if not cudnn:
        sys.exit(
            "ldx.gpu is true but libcudnn.so.9 was not found in the `gpu` pixi "
            "environment.\nInstall it with `pixi install -e gpu` — without it "
            "the CUDA provider silently falls back to CPU."
        )
    return hits[0]


rule escapepod_demux:
    """
    Demultiplex LDX (nbc) barcodes with escapepod's fused CTC-CRF pipeline.

    One pass over the raw POD5 does detect -> prep -> basecall -> match -> route,
    writing the per-barcode POD5 files directly. There is deliberately no follow-up
    rule to derive the split or a read->barcode mapping: the WarpDemuX path needs
    those because its classifier only emits a table, whereas here the routed POD5
    *is* the product and re-deriving it would be a second full pass for an
    identical result.

    No --barcodes or --method is passed: the bundle carries its own references and
    pins the boundary detector it was calibrated against, and overriding either
    silently degrades the calls.

    The classifications CSV is kept because it is the only per-read record of the
    call and its confidence margin — the POD5 routing preserves which barcode won,
    but not by how much, and that margin is what separates a confident call from a
    near-tie when auditing demux quality.
    """
    input:
        get_run_raw_inputs,
    output:
        outdir=directory(os.path.join(outdir, "demux", "escapepod", "{run_id}")),
        classifications=os.path.join(
            outdir, "demux", "read_ids", "{run_id}", "classifications.csv"
        ),
        summary=os.path.join(
            outdir, "demux", "read_ids", "{run_id}", "demux_summary.tsv.gz"
        ),
    log:
        os.path.join(outdir, "logs", "escapepod_demux", "{run_id}"),
    threads: demux_config().get("threads", 16)
    params:
        model=get_demux_model(),
        min_margin=demux_config().get("min_margin", 0),
        # Overrules the model bundle's declared `boundary.margin` — the samples
        # of adapter_end a read needs beyond the model's chunk before the CRF
        # will decode it. Unset (the default) leaves the bundle in charge, which
        # is where this belongs; set it only to evaluate a change the bundle has
        # not adopted yet. Needs escpod with the flag (escapepod-rs#193).
        boundary_margin=lambda wildcards: (
            f"--boundary-margin {demux_config()['boundary_margin']}"
            if demux_config().get("boundary_margin") is not None
            else ""
        ),
        # The sibling of --boundary-margin, for reads whose adapter ends BEFORE
        # chunk, where relaxing the margin cannot help because the window would
        # start before sample 0. Same null-means-the-bundle-decides contract.
        clamp_max_shift=lambda wildcards: (
            f"--clamp-max-shift {demux_config()['clamp_max_shift']}"
            if demux_config().get("clamp_max_shift") is not None
            else ""
        ),
        # --gpu runs the CRF encoder and the boundary CNN through onnxruntime's
        # CUDA provider; the lattice decode stays on the CPU either way. It is
        # opt-in because the *released* escpod has neither GPU feature compiled
        # in and would reject the flag — see config-base.yml `ldx.gpu`.
        # `--method cnn` is passed alongside --gpu on purpose. The bundle already
        # pins cnn, but a pinned detector runs on the CPU and makes --gpu a
        # no-op; only an explicit --method engages the CUDA detection path.
        # Naming the same detector the bundle pins is not an override, so this
        # cannot silently downgrade to LLR.
        gpu=lambda wildcards: (
            "--gpu --method cnn" if demux_config().get("gpu", False) else ""
        ),
        # ort dlopens onnxruntime at run time from ORT_DYLIB_PATH, and the CUDA
        # execution provider sits beside it, so its directory must also be on
        # LD_LIBRARY_PATH or the core library loads and the provider then fails
        # to. Exported inside the rule rather than relying on the submitting
        # shell's environment surviving into the batch job.
        ort_env=lambda wildcards: (
            f"export ORT_DYLIB_PATH={get_ort_dylib()}; "
            f"export LD_LIBRARY_PATH={os.path.dirname(get_ort_dylib())}:"
            f"{os.path.join(PIPELINE_DIR, '.pixi', 'envs', 'gpu', 'lib')}:"
            f"$LD_LIBRARY_PATH; "
            if demux_config().get("gpu", False)
            else ""
        ),
        pod5_dirs=lambda wildcards: " ".join(
            os.path.join(d, "*.pod5") for d in get_run_pod5_dirs(wildcards.run_id)
        ),
        summarize_awk=os.path.join(SCRIPT_DIR, "summarize_demux.awk"),
    shell:
        """
        # escpod opens --classifications with a plain create(), so a missing
        # parent directory is a bare ENOENT ("No such file or directory") with
        # nothing naming the path. Snakemake does not reliably pre-create it
        # here: when a previous attempt fails it removes this rule's outputs,
        # taking demux/read_ids/<run_id>/ with them, and the retry then dies on
        # the very first write. Creating them up front makes the rule
        # re-runnable after any failure.
        mkdir -p $(dirname {output.classifications}) \
            $(dirname {output.summary}) \
            {output.outdir}

        {params.ort_env}escpod demux {params.pod5_dirs} \
            --model {params.model} \
            --output-dir {output.outdir} \
            --classifications {output.classifications} \
            --min-margin {params.min_margin} \
            {params.boundary_margin} \
            {params.clamp_max_shift} \
            {params.gpu} \
            --threads {threads} 2>&1 | tee {log}

        # escpod prints its per-barcode tally to the log only, so tabulate the
        # same counts into a file the QC report can read. The awk program lives
        # in workflow/scripts/ rather than inline; see the note at the top of it.
        awk -F, -f {params.summarize_awk} {output.classifications} \
            | gzip >{output.summary}
        """


def get_sample_escapepod_dir(wildcards):
    """Locate the escapepod demux output directory for a sample's run."""
    run_id = samples[wildcards.sample]["run_id"]
    return os.path.join(outdir, "demux", "escapepod", run_id)


rule link_demux_pod5:
    """
    Adopt escapepod's per-barcode POD5 as the sample's split POD5.

    escapepod already routed every read into its barcode's file with a block-level
    copy during the demux pass, so re-deriving the same split with `pod5 filter`
    would be a second full pass for an identical result. This rule only renames
    barcode -> sample.

    The file to adopt is named for the barcode the MODEL emits, which is not
    necessarily the name the samples file used: the WDX panel is configured as
    `barcode03` but emitted as `bc03`. label_to_emitted() is what bridges that —
    see workflow/scripts/barcode_names.py.
    """
    input:
        demux_dir=get_sample_escapepod_dir,
    output:
        pod5=maybe_temp(
            os.path.join(outdir, "demux", "pod5", "{sample}", "{sample}.pod5"),
            tier="split_pod5",
        ),
    log:
        os.path.join(outdir, "logs", "link_demux_pod5", "{sample}"),
    params:
        barcode=lambda wildcards: samples[wildcards.sample]["barcode"],
        emitted=lambda wildcards: label_to_emitted(samples[wildcards.sample]["barcode"]),
    run:
        # escapepod names outputs <prefix>_<barcode>.pod5 and skips barcodes
        # with no reads, so a missing file means this barcode got nothing.
        src = Path(input.demux_dir) / f"barcode_{params.emitted}.pod5"
        if not src.exists():
            named = (
                f"'{params.barcode}'"
                if params.emitted == params.barcode
                else f"'{params.barcode}' (emitted as '{params.emitted}')"
            )
            raise WorkflowError(
                f"No reads were assigned to barcode {named} "
                f"(sample '{wildcards.sample}'): {src} was not written.\n"
                f"Check the per-barcode counts in "
                f"{Path(input.demux_dir).parent.parent}/read_ids/"
                f"{samples[wildcards.sample]['run_id']}/demux_summary.tsv.gz"
            )
        dest = Path(output.pod5)
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists() or dest.is_symlink():
            dest.unlink()
        # Hardlink rather than symlink: both are free, but a hardlink keeps the
        # POD5 reachable if demux/escapepod is cleaned out from under it, which
        # a symlink would turn into a dangling path that only fails much later
        # at classify time. Falls back to a symlink across filesystems.
        try:
            os.link(src, dest)
        except OSError:
            os.symlink(os.path.relpath(src.resolve(), dest.parent), dest)
        with open(log[0], "w") as f:
            f.write(f"{src} -> {dest}\n")


# Route the two backends. They overlap on exactly two outputs — the per-sample
# split POD5 and the per-run demux summary — so each needs an explicit winner.
# is_demux_enabled() has already rejected a config that turns on both backends.
if is_crf_demux_enabled():

    ruleorder: link_demux_pod5 > split_pod5
    ruleorder: escapepod_demux > parse_warpdemux

else:

    ruleorder: split_pod5 > link_demux_pod5
    ruleorder: parse_warpdemux > escapepod_demux


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
        # Several run directories in one call: multi-input filter landed in
        # escapepod-rs 0.8.1 (#196), which config-base.yml now pins. Before that
        # escpod took a single <INPUT> and would have filtered against only the
        # first directory — quietly, which is why this call site waited for a
        # release rather than an unreleased commit.
        escpod filter {params.pod5_dirs} -i {input.read_ids} -o {output} 2>&1 | tee {log}
        """


# --- EDX (3' adapter barcode) early demultiplexing ---


def get_edx_samples():
    """Sample names to run 3' adapter detection on.

    Normally that is the samples carrying an explicit `edx:` assignment, since
    detection exists to serve the EDX filter. With `edx.detect_all`, it is every
    sample instead: an adapter-ligation QC run pools several 3' adapters behind
    one signal barcode, so no sample has a single `edx:` to filter to, yet the
    adapter composition of each barcode is the whole measurement.
    """
    if config.get("edx", {}).get("detect_all", False):
        return list(samples.keys())
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
                # read_id, adapter_3p, score_best, score_second, margin
                fields = line.rstrip("\n").split("\t")
                read_id, adapter = fields[0], fields[1]
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
    shell:
        """
        escpod filter {input.pod5} -i {input.read_ids} -o {output.pod5} 2>&1 | tee {log}
        """


# --- EDX (3' adapter barcode) concordance analysis ---


rule edx_concordance:
    """
    Build the signal-barcode x EDX (3' adapter) contingency table.

    Uses pre-alignment adapter detection TSVs, which contain ALL reads with their
    detected adapter, rather than final BAMs, which only contain matching reads —
    the reads that ended up under the wrong adapter are exactly what this measures,
    so a filtered input would hide the signal.
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
        min_margin=config.get("edx", {}).get("min_margin", 5),
    shell:
        """
        python {params.src}/edx_concordance.py \
            --tsvs {input.tsvs} \
            --samples {params.sample_names} \
            --min-margin {params.min_margin} \
            --output {output.concordance} \
            2>&1 | tee {log}
        """
