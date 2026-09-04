"""
Rules for signal-level barcode demultiplexing and EDX (3' adapter barcode)
concordance. Loaded when either demux backend is enabled in config.

The two backends split the run at different points:

  warpdemux : WarpDemuX DTW/fingerprint classifier, WDX barcodes ("barcode04").
              Classifies into a table, which is then parsed into a read->barcode
              mapping and applied with `escpod filter` to cut the POD5 per
              sample. Each split POD5 is then basecalled on its own.
  ldx       : escapepod CTC-CRF barcode basecaller, LDX barcodes ("nbc01").
              One fused pass detects, basecalls and matches each read, recording
              the assignment in a `.p5s` sidecar beside the raw POD5. No POD5 is
              written at all: the run is basecalled once, whole, and the split
              happens on the resulting uBAM.

The LDX route exists to keep exactly one copy of the signal on disk. Routing
reads into per-barcode POD5s duplicates the entire run for the life of the
analysis; a sidecar is a few MB and leaves the raw POD5 untouched, which also
lets a re-demux (new model, new margin) cost nothing but the pass itself.

Both backends converge on the per-sample unaligned BAM
(`bam/rebasecall/{sample}/{sample}.rbc.bam`), and everything downstream of that
is identical. Only one may be enabled (see is_demux_enabled); a ruleorder picks
which rule supplies it.
"""

import pandas as pd
import re
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
    """Get barcode→sample mapping for a run, keyed as the MODEL emits it.

    The keys are matched against the `barcode` column of classifications.csv,
    which carries the bundle's own vocabulary. A WDX4 sample configured as
    `barcode03` is recorded there as `bc03`, so keying on the configured string
    would match nothing and every sample would come back empty; an ldx16 sample
    configured `ldx01` is already the bundle's own name and passes through.

    Resolved against the CONFIGURED bundle rather than by prefix, so a samples
    file and a model that disagree are caught here, while the DAG is being
    built, instead of by the "No reads were assigned" guard in
    rebasecall_ldx_run — which fires only after the demux pass has run.
    """
    bundle = get_ldx_model()
    return {
        resolve_to_bundle(info["barcode"], bundle): sample
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


# --- escapepod CTC-CRF (LDX) demultiplexing ---


def get_ldx_model():
    """Path to the CRF bundle directory, resolved against the pipeline dir.

    Config carries a repo-relative path so it stays portable, but Snakemake runs
    with the working directory set by the caller, so relative paths cannot be
    handed to the shell as-is.
    """
    model = config.get("ldx", {}).get("model")
    if not model:
        sys.exit(
            "ldx.enabled is true but ldx.model is unset. Point it at a CRF "
            "bundle directory, e.g. "
            "resources/models/demux/barcode_crf_ldx16_rna004@v0.1.0"
        )
    if not os.path.isabs(model):
        model = os.path.join(PIPELINE_DIR, model)
    if not os.path.isdir(model):
        sys.exit(
            f"ldx.model is not a directory: {model}\n"
            "escapepod CRF models are self-describing BUNDLES (metadata.json "
            "plus the ONNX graphs it names), not a single file."
        )
    # The nbc16 family is retired here, and refused rather than merely
    # undocumented. It is not a naming preference: nbc16 and ldx16 are separate
    # RETRAINS whose calls differ on ~8% of reads, so a run demultiplexed with
    # one is not comparable to a run demultiplexed with the other, and nothing
    # downstream records which was used. The vocabulary is also no longer
    # translated (see workflow/scripts/barcode_names.py), so an nbc bundle would
    # now fail in resolve_to_bundle anyway — this says why, at the point the
    # model is chosen.
    if os.path.basename(model).startswith("barcode_crf_nbc"):
        sys.exit(
            f"ldx.model points at a retired nbc bundle: {os.path.basename(model)}\n"
            "This pipeline runs the ldx16 panel, whose references are named "
            "`ldx01`..`ldx16`. nbc16 is a different retrain — ~8% of calls "
            "differ — so switching is not a rename and results are not "
            "comparable across it.\n"
            "Use resources/models/demux/barcode_crf_ldx16_rna004@v0.1.0."
        )
    return model


def get_escpod_bin():
    """The escpod binary `escapepod_demux` should run.

    GPU placement needs a build carrying the crf-gpu/cnn-gpu features, which
    the portable musl artifact does not have. Since escapepod-rs 0.17.1 that
    build is PUBLISHED (`...-x86_64-unknown-linux-gnu-gpu.tar.gz`) and
    `pixi run setup` downloads it into `<version>-gpu/`; before 0.17.1 the same
    path had to be produced by a source build of the private repo. The path
    convention is unchanged, so only the remediation in the error below moved.

    Resolved HERE rather than by pointing `escpod_version` at `-gpu`, because
    demux is the only rule with a GPU path and the GPU artifact is x86_64 Linux
    only. Moving the global pin would hand every other rule a dynamically
    linked, single-platform binary just to run `escpod merge` and
    `escpod classify`, and would break `pixi run setup`, which derives
    its download URL from that same string and would ask GitHub for a
    `v<version>-gpu` release that does not exist.
    """
    if not config.get("ldx", {}).get("gpu", False):
        return "escpod"
    version = config.get("escpod_version", ESCPOD_VERSION)
    binary = os.path.join(
        PIPELINE_DIR, "resources", "tools", "escpod", f"{version}-gpu", "bin", "escpod"
    )
    if not os.path.isfile(binary):
        sys.exit(
            f"ldx.gpu is true but no GPU-enabled escpod was found at {binary}.\n"
            "Install it with `pixi run setup`, which downloads the published "
            f"GPU artifact for escpod {version} (x86_64 Linux only).\n"
            "It also needs a CUDA libonnxruntime (`pixi run install-ort-gpu`) "
            "and cuDNN (`pixi install -e gpu`).\n"
            "Set ldx.gpu: false to run demux on the CPU instead."
        )
    return binary


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

    One pass over the raw POD5 does detect -> prep -> basecall -> match, and
    `--annotate` records each read's barcode in a `.p5s` sidecar. With no
    -d/--output-dir this is the only routing output: no per-barcode POD5 is
    written, so the run is never duplicated. The split is taken later, on the
    basecalled uBAM (see split_ldx_ubam).

    Since escpod 0.19.0 a DIRECTORY argument produces ONE collection sidecar for
    the whole directory — `<run>/pod5/` gets `<run>/pod5.p5s` — rather than a
    `*.pod5.p5s` beside every member. A run that wrote fifty POD5s produced one
    set of barcode calls, not fifty, and this is the shape that says so. It is
    the argument that selects it, not the flag: naming files individually still
    writes one sidecar each, which is why params.pod5_dirs passes the
    directories.

    The sidecar lives next to the raw data rather than under this run's results
    because escpod binds it to the POD5s' footer UUIDs and sizes, and reads it
    from the adjacent path. A member gets rows only when both match its entry, so
    a file dropped into the directory afterwards is told it has no sidecar rather
    than inheriting a neighbour's labels. Three consequences worth knowing:

      - Re-demuxing the same directory replaces the `barcode` column in place, so
        changing the model or the margin is safe and needs no cleanup.
      - A POD5 replaced under the same name leaves a collection describing reads
        that are no longer there. escpod detects this from the footer and refuses
        to read it; recover by deleting the `.p5s` and re-running. (Deleting a
        `.p5s` by itself does not re-trigger this rule, since its outputs are
        still present — force it with `--forcerun escapepod_demux`.)
      - An escpod older than 0.19.0 refuses a collection sidecar BY NAME rather
        than reporting a missing column, so a tree written by this pipeline and
        read by an older binary fails loudly instead of silently.
      - A per-file `*.pod5.p5s` left behind by a pre-0.19.0 demux SHADOWS the
        collection: lookup merges the two per column with the file's own winning,
        so a stale `barcode` column would outrank the fresh one. Nothing here
        reads a sidecar — the split is driven by classifications.csv — so this
        cannot affect the pipeline's own output, but delete the old per-file
        sidecars before pointing `escpod demux split --sidecar` or the Python
        Reader at a run demuxed under both versions.

    No --barcodes or --method is passed: the bundle carries its own references and
    pins the boundary detector it was calibrated against, and overriding either
    silently degrades the calls.

    The classifications CSV is kept for two reasons. It is the read->barcode source
    the basecall and split rules downstream read, so it must outlive the demux
    pass -- there is no routed POD5 to recover it from any more. And it is the
    only per-read record of how the call went: under `ldx.ref_scores` (on by
    default) it carries the lattice's own log P(barcode | signal), which is the
    score with enough resolution to audit against and the one `ldx.min_crf_margin`
    gates on. `confidence` alone is an edit-distance margin that a designed panel
    pins to three values.
    """
    input:
        get_run_raw_inputs,
    output:
        classifications=os.path.join(
            outdir, "demux", "read_ids", "{run_id}", "classifications.csv"
        ),
        summary=os.path.join(
            outdir, "demux", "read_ids", "{run_id}", "demux_summary.tsv.gz"
        ),
    log:
        os.path.join(outdir, "logs", "escapepod_demux", "{run_id}"),
    threads: config.get("ldx", {}).get("threads", 16)
    params:
        model=get_ldx_model(),
        escpod=get_escpod_bin(),
        min_margin=config.get("ldx", {}).get("min_margin", 0),
        # The CRF's own confidence. `--ref-scores` restricts the forward
        # recursion to the paths emitting each reference and normalises by the
        # partition function, so classifications.csv gains crf_logp /
        # crf_margin / crf_best / mean_logpost — a continuous per-read score,
        # where `confidence` is an edit-distance margin a designed panel pins
        # to three values. Needs escpod >= 0.12.0 (escapepod-rs#241); an older
        # binary rejects the flag outright rather than ignoring it, which is
        # the failure mode we want given the pin lives in the same config.
        ref_scores=lambda wildcards: (
            "--ref-scores" if config.get("ldx", {}).get("ref_scores", True) else ""
        ),
        # Gates on that score, both unset by default — see config-base.yml for
        # why an unmeasured operating point is not shipped. Either flag implies
        # --ref-scores upstream, so a gate cannot silently do nothing.
        min_crf_margin=lambda wildcards: (
            f"--min-crf-margin {config['ldx']['min_crf_margin']}"
            if config.get("ldx", {}).get("min_crf_margin") is not None
            else ""
        ),
        min_crf_prob=lambda wildcards: (
            f"--min-crf-prob {config['ldx']['min_crf_prob']}"
            if config.get("ldx", {}).get("min_crf_prob") is not None
            else ""
        ),
        # Boundary gating. NO upstream CRF bundle declares `boundary.margin` or
        # `boundary.clamp_max_shift` — build_crf_bundle.py cannot write them —
        # so unset does not mean "the bundle decides", it means escpod's
        # fallback of margin 200 and no clamp. Both are therefore configured
        # explicitly in config-base.yml, which documents the measurements.
        # Needs escpod with the flags (escapepod-rs#193).
        boundary_margin=lambda wildcards: (
            f"--boundary-margin {config['ldx']['boundary_margin']}"
            if config.get("ldx", {}).get("boundary_margin") is not None
            else ""
        ),
        # Reaches the reads --boundary-margin cannot: those whose adapter ends
        # before the chunk, decoded from [0, chunk] instead. 0 disables.
        clamp_max_shift=lambda wildcards: (
            f"--clamp-max-shift {config['ldx']['clamp_max_shift']}"
            if config.get("ldx", {}).get("clamp_max_shift") is not None
            else ""
        ),
        # Device placement. `--device gpu` runs the CRF encoder and the
        # boundary CNN through onnxruntime's CUDA provider; the lattice decode
        # stays on the CPU either way.
        #
        # Passed EXPLICITLY in both directions rather than leaning on the
        # `auto` default, because the two failure modes are not symmetric.
        # `--device gpu` is a requirement: it fails when the feature is
        # missing, no device is visible, or onnxruntime cannot register its
        # CUDA provider. `auto` would quietly run on the CPU instead — a 20x
        # slowdown (2.6 h against 458.9 s on our own flowcell) that still
        # produces correct output, so it looks like success and is caught only
        # by noticing the wall clock. Under `ldx.gpu: false` the explicit
        # `--device cpu` likewise stops a GPU binary from opportunistically
        # using a device the config said not to.
        #
        # Replaces `--gpu`, deprecated in escapepod-rs 0.17.1. That spelling
        # still runs (it warns and continues as `--device gpu`), so this is a
        # migration rather than a break — but the meaning changed underneath
        # it: `--gpu` fell back to the CPU where `--device gpu` fails.
        #
        # `--method cnn` is passed alongside on purpose. The bundle already
        # pins cnn, and naming the same detector the bundle pins is not an
        # override, so this cannot silently downgrade to LLR.
        gpu=lambda wildcards: (
            "--device gpu --method cnn"
            if config.get("ldx", {}).get("gpu", False)
            else "--device cpu"
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
            if config.get("ldx", {}).get("gpu", False)
            else ""
        ),
        # DIRECTORIES, not a `*.pod5` glob, and the distinction is the whole
        # point since escpod 0.19.0: pointed at a directory `--annotate` writes
        # ONE collection sidecar beside it, pointed at files it writes one per
        # file. Naming the directories is what turns fifty sidecars into one.
        pod5_dirs=lambda wildcards: " ".join(get_run_pod5_dirs(wildcards.run_id)),
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
            $(dirname {output.summary})

        {params.ort_env}{params.escpod} demux {params.pod5_dirs} \
            --model {params.model} \
            --annotate \
            --classifications {output.classifications} \
            --min-margin {params.min_margin} \
            {params.ref_scores} \
            {params.min_crf_margin} \
            {params.min_crf_prob} \
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


def ldx_sample_constraint():
    """Regex matching exactly the samples the LDX rules may claim.

    The rules below produce the same per-sample paths as `rebasecall` and
    `extract_sample_reads` and win them by ruleorder, which would otherwise also
    capture the un-barcoded samples of a mixed run — samples with no barcode to
    split on. Constraining the wildcard leaves those to `rebasecall`, and makes
    these rules unmatchable on a WarpDemuX run rather than merely outranked.
    """
    if not is_ldx_enabled():
        return r"$^"
    names = [re.escape(s) for s, info in samples.items() if info.get("barcode")]
    return "|".join(names) if names else r"$^"


def get_sample_run_classifications(wildcards):
    """Locate the demux classifications CSV for a sample's run."""
    run_id = samples[wildcards.sample]["run_id"]
    return os.path.join(outdir, "demux", "read_ids", run_id, "classifications.csv")


def get_sample_run_split_parents(wildcards):
    """Locate the split-read parent map for a sample's run."""
    run_id = samples[wildcards.sample]["run_id"]
    return os.path.join(outdir, "demux", "read_ids", run_id, "split_parents.tsv")


def get_sample_run_ubam(wildcards):
    """Locate the run-level basecall for a sample's run."""
    run_id = samples[wildcards.sample]["run_id"]
    return os.path.join(outdir, "bam", "rebasecall_run", run_id, f"{run_id}.rbc.bam")


rule rebasecall_ldx_run:
    """
    Basecall a whole LDX run in one dorado pass, straight off the raw POD5.

    The WarpDemuX path basecalls each sample's split POD5 separately because that
    split is what its `escpod filter` step produced. Here there is no split POD5 —
    demux left only a sidecar — so the run is basecalled whole and cut up
    afterwards, on the uBAM.

    `-l` restricts the pass to reads that demux actually assigned to one of this
    run's samples. That is the same read universe the per-barcode POD5s used to
    cover: unclassified reads were never in any of them, and basecalling them here
    would be GPU time spent on output nothing reads.
    """
    input:
        pod5=get_run_raw_inputs,
        classifications=os.path.join(
            outdir, "demux", "read_ids", "{run_id}", "classifications.csv"
        ),
        mod_models=rules.download_mod_models.output.sentinel,
    output:
        bam=maybe_temp(
            os.path.join(
                outdir, "bam", "rebasecall_run", "{run_id}", "{run_id}.rbc.bam"
            ),
            tier="basecall",
        ),
        read_ids=maybe_temp(
            os.path.join(
                outdir, "demux", "read_ids", "{run_id}", "assigned_read_ids.txt"
            ),
            tier="demux_scratch",
        ),
    log:
        os.path.join(outdir, "logs", "rebasecall_ldx_run", "{run_id}"),
    params:
        model=config["base_calling_model"],
        dorado_opts=config["opts"]["dorado"],
        models_dir=os.path.join(PIPELINE_DIR, "resources", "models"),
        run_path=lambda wildcards: get_run_path(wildcards.run_id),
        barcodes=lambda wildcards: ",".join(get_barcodes_for_run(wildcards.run_id)),
        resume_sh=os.path.join(SCRIPT_DIR, "dorado_basecall_resume.sh"),
    shell:
        """
        # Exact set membership rather than a regex: barcode names come from the
        # model bundle, and matching them loosely is how nbc1 would swallow
        # nbc16.
        awk -F, -v bcs="{params.barcodes}" \
            'BEGIN {{ n = split(bcs, a, ","); for (i = 1; i <= n; i++) keep[a[i]] = 1 }}
             NR > 1 && ($2 in keep) {{ print $1 }}' \
            {input.classifications} >{output.read_ids}

        if [ ! -s {output.read_ids} ]; then
            echo "ERROR: demux assigned no reads to any barcode of run {wildcards.run_id}." >&2
            exit 1
        fi

        if [[ "${{CUDA_VISIBLE_DEVICES:-}}" ]]; then
            echo "CUDA_VISIBLE_DEVICES $CUDA_VISIBLE_DEVICES"
            export CUDA_VISIBLE_DEVICES
        fi

        # Via the resume wrapper rather than a bare redirect. This is the rule
        # whose walltime is hardest to size -- one job for a whole flowcell, hours
        # of GPU -- and without a resume an overrun destroyed all of it rather
        # than a tail, so the wall had to be padded rather than measured. The
        # partial it leaves behind is now a checkpoint; see the script.
        bash {params.resume_sh} {output.bam} \
            --models-directory {params.models_dir} {params.dorado_opts} \
            {params.model} {params.run_path} --recursive \
            --read-ids {output.read_ids}
        """


rule ldx_split_parent_map:
    """
    Record which basecalled reads are split children, and of whom.

    dorado splits a concatenated signal read into several reads, each with a NEW
    read id and its origin in the `pi` tag. Those ids appear in no classifications
    row — demux only ever saw the parent — so a split that keyed on the CSV alone
    would drop every split child on the floor, silently. (Measured on the 2026-08-06
    run: 0.85% of a sample's reads.) One scan here lets the per-sample extraction
    give each child its parent's barcode.
    """
    input:
        bam=os.path.join(
            outdir, "bam", "rebasecall_run", "{run_id}", "{run_id}.rbc.bam"
        ),
    output:
        parents=maybe_temp(
            os.path.join(outdir, "demux", "read_ids", "{run_id}", "split_parents.tsv"),
            tier="demux_scratch",
        ),
    log:
        os.path.join(outdir, "logs", "ldx_split_parent_map", "{run_id}"),
    threads: 4
    shell:
        """
        samtools view -@ {threads} {input.bam} \
            | awk -v OFS='\\t' \
                '{{ for (i = 12; i <= NF; i++) if ($i ~ /^pi:Z:/) {{ print $1, substr($i, 6); break }} }}' \
                >{output.parents} 2>{log}
        """


rule extract_ldx_sample_reads:
    """
    List the reads belonging to one sample: its barcode's reads, plus any split
    children of those reads.
    """
    input:
        classifications=get_sample_run_classifications,
        parents=get_sample_run_split_parents,
    output:
        read_ids=maybe_temp(
            os.path.join(outdir, "demux", "read_ids", "{sample}", "{sample}.txt"),
            tier="demux_scratch",
        ),
    log:
        os.path.join(outdir, "logs", "extract_ldx_sample_reads", "{sample}"),
    wildcard_constraints:
        sample=ldx_sample_constraint(),
    params:
        # As the MODEL emits it: this is compared against classifications.csv,
        # which speaks the bundle's vocabulary, not the samples file's. Asked of
        # the bundle rather than guessed from the prefix — see resolve_to_bundle.
        barcode=lambda wildcards: resolve_to_bundle(
            samples[wildcards.sample]["barcode"], get_ldx_model()
        ),
    run:
        import csv

        assigned = set()
        with open(input.classifications, newline="") as f:
            for row in csv.DictReader(f):
                if row["barcode"] == params.barcode:
                    assigned.add(row["read_id"])
        if not assigned:
            raise WorkflowError(
                f"No reads were assigned to barcode '{params.barcode}' "
                f"(sample '{wildcards.sample}').\n"
                f"Check the per-barcode counts in "
                f"{os.path.dirname(input.classifications)}/demux_summary.tsv.gz"
            )
        # A split child inherits its parent's barcode. The parent id itself is
        # not in the BAM — dorado replaces the read with its children — so
        # keeping it in the list is harmless and keeps this file a faithful
        # record of what demux assigned.
        children = set()
        with open(input.parents) as f:
            for line in f:
                child, parent = line.rstrip("\n").split("\t")
                if parent in assigned:
                    children.add(child)
        with open(output.read_ids, "w") as f:
            for read_id in sorted(assigned | children):
                f.write(f"{read_id}\n")
        with open(log[0], "w") as f:
            f.write(
                f"{len(assigned)} assigned reads, {len(children)} split children\n"
            )


rule split_ldx_ubam:
    """
    Cut one sample's unaligned BAM out of the run-level basecall.

    Emits the path `rebasecall` would have, so everything downstream — EDX adapter
    detection, FASTQ extraction, alignment, tag transfer — is unchanged and cannot
    tell which backend produced it.
    """
    input:
        bam=get_sample_run_ubam,
        read_ids=rules.extract_ldx_sample_reads.output.read_ids,
    output:
        maybe_temp(
            os.path.join(outdir, "bam", "rebasecall", "{sample}", "{sample}.rbc.bam"),
            tier="basecall",
        ),
    log:
        os.path.join(outdir, "logs", "split_ldx_ubam", "{sample}"),
    wildcard_constraints:
        sample=ldx_sample_constraint(),
    threads: 4
    shell:
        """
        samtools view -b -@ {threads} -N {input.read_ids} {input.bam} \
            >{output} 2>{log}
        """


# Route the two backends. They overlap on three outputs — the per-sample read-id
# list, the per-sample uBAM and the per-run demux summary — so each needs an
# explicit winner. is_demux_enabled() has already rejected a config that turns
# on both backends.
#
# The LDX rules are additionally constrained to barcoded samples of an LDX run
# (ldx_sample_constraint), so on the WarpDemuX path they cannot match at all and
# these orderings are belt and braces.
if is_ldx_enabled():

    ruleorder: extract_ldx_sample_reads > extract_sample_reads
    ruleorder: split_ldx_ubam > rebasecall
    ruleorder: escapepod_demux > parse_warpdemux

else:

    ruleorder: extract_sample_reads > extract_ldx_sample_reads
    ruleorder: rebasecall > split_ldx_ubam
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
