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

The escpod path can carry a SECOND AXIS, `fdx` (`fdx.enabled`): a 5' index on the
same molecule as the 3' LDX code, called by its own pass over the raw POD5
(escapepod_demux_fdx) and joined against the LDX call per read by
select_demux_reads.py. A dual-index sample is `{ldx: ldx01, fdx: fdx01}` and owns
the reads on which both calls agree. Why two passes rather than one fused run is
recorded on that rule and in config-base.yml.
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


# `{run_id}` may only be a run this samples file names. Without this, the fdx
# axis's `demux/read_ids/{run}/fdx/classifications.csv` also matches the LDX
# rule's `demux/read_ids/{run_id}/classifications.csv` with run_id = "run/fdx",
# and Snakemake evaluates that rule's input function (which cannot resolve the
# path) before settling on the right rule -- a logged error on every dry-run,
# and one candidate away from an ambiguous-rule failure.
wildcard_constraints:
    run_id="|".join(re.escape(r) for r in get_run_ids()) or r"$^",


def get_samples_for_run(run_id):
    """Get list of sample names assigned to a run."""
    return [
        sample
        for sample, info in samples.items()
        if info.get("run_id") == run_id and info.get("barcode")
    ]


def get_sample_axis_codes(sample):
    """{axis: code} for an escpod-demultiplexed sample, each code as its MODEL emits it.

    The codes are matched against the `barcode` column of that axis's
    classifications.csv, which carries the bundle's own vocabulary. A WDX4
    sample configured as `barcode03` is recorded there as `bc03`, so keying on
    the configured string would match nothing and every sample would come back
    empty; an ldx16 sample configured `ldx01` is already the bundle's own name
    and passes through, as is an fdx4 one configured `fdx01`.

    Resolved against the CONFIGURED bundle rather than by prefix, so a samples
    file and a model that disagree are caught while the DAG is being built
    (see validate_escpod_barcodes) instead of by the "no reads were assigned"
    guard in select_demux_reads.py, which fires only after the demux pass has
    run. Axis order is the join order: `ldx` first, because that pass is the one
    the run's accounting is keyed on and the one every sample must name.
    """
    info = samples[sample]
    codes = {"ldx": resolve_to_bundle(info["barcode"], get_ldx_model())}
    if info.get("fdx"):
        codes["fdx"] = resolve_to_bundle(info["fdx"], get_fdx_model())
    return codes


def get_sample_downstream_codes(sample):
    """{axis: code} a sample carries that select_demux_reads.py does NOT select on.

    Only `edx`, and only because it is not a signal axis: the 3' adapter is read
    off the uBAM by detect_edx_adapters, long after the demux pass, so there is
    no classifications CSV for it here.

    It is kept OUT of get_sample_axis_codes() deliberately -- that function also
    names which CSVs to read (get_sample_axis_csv_args), and an `edx` entry
    there would look up a path that does not exist. It is passed to the script
    instead as a `--downstream-axis`, where it counts towards sample uniqueness
    and nothing else. Without it, samples that share an LDX code and differ only
    by adapter are indistinguishable to that script's duplicate check and the
    run dies after the demux pass has been paid for (#163).
    """
    if is_edx_enabled() and sample_has_edx(sample):
        return {"edx": samples[sample]["edx"]}
    return {}


def get_downstream_axis_args():
    """`--downstream-axis` flags for the axes resolved later in the pipeline."""
    return "--downstream-axis edx" if is_edx_enabled() else ""


def get_sample_axis_args(sample):
    """The `--sample` argument select_demux_reads.py takes for one sample."""
    codes = {**get_sample_axis_codes(sample), **get_sample_downstream_codes(sample)}
    return f"{sample}:" + ",".join(f"{axis}={code}" for axis, code in codes.items())


def get_run_axis_csvs(run_id):
    """{axis: classifications.csv} for every escpod axis enabled on a run.

    Two passes give two CSVs; a fused pass gives one CSV with per-axis columns,
    and then both axes point at it. select_demux_reads.py reads either shape.
    """
    base = os.path.join(outdir, "demux", "read_ids", run_id)
    axes = {"ldx": os.path.join(base, "classifications.csv")}
    if is_fdx_fused():
        axes["fdx"] = axes["ldx"]
    elif is_fdx_enabled():
        axes["fdx"] = os.path.join(base, "fdx", "classifications.csv")
    return axes


def get_axis_gate_args():
    """`--gate AXIS=NATS` for every axis with a configured lattice gate.

    escpod applies one `--min-crf-margin` to every model of a fused pass, and
    the axes want different operating points (ldx 1.0, fdx 3.5), so the
    per-axis gate is applied in the join as well. Re-applying a gate escpod
    already applied changes nothing, which is what lets the two-pass and fused
    shapes share this code.
    """
    gates = []
    for axis in ("ldx", "fdx"):
        if axis == "fdx" and not is_fdx_enabled():
            continue
        margin = config.get(axis, {}).get("min_crf_margin")
        if margin is not None:
            gates.append(f"--gate {axis}={margin}")
    return " ".join(gates)


def get_escpod_model_args():
    """The `--model` argument(s) of the primary escpod pass.

    Unnamed for a single axis, so the sidecar column stays `barcode` and the
    CSV keeps the single-model shape every existing consumer reads. Named
    (`ldx=`, `fdx=`) for a fused pass, which is what makes the sidecar and CSV
    columns per-axis.
    """
    if is_fdx_fused():
        return f"--model ldx={get_ldx_model()} --model fdx={get_fdx_model()}"
    return f"--model {get_ldx_model()}"


def get_sample_axis_csvs(wildcards):
    """The classifications CSV(s) a sample's selection reads, each path once.

    A fused pass puts both axes in one file, so this is deduplicated rather than
    one-per-axis; the `--axis` arguments are built from the run's axis map
    instead (get_sample_axis_csv_args), which keeps name and path paired.
    """
    axes = get_run_axis_csvs(samples[wildcards.sample]["run_id"])
    return sorted({axes[axis] for axis in get_sample_axis_codes(wildcards.sample)})


def get_sample_axis_csv_args(wildcards):
    """`--axis NAME=CSV` arguments, in the same order as the sample's codes."""
    axes = get_run_axis_csvs(samples[wildcards.sample]["run_id"])
    return " ".join(
        f"--axis {axis}={axes[axis]}"
        for axis in get_sample_axis_codes(wildcards.sample)
    )


def get_run_axis_csv_args(wildcards):
    """`--axis NAME=CSV` for every axis enabled on a run."""
    return " ".join(
        f"--axis {axis}={path}"
        for axis, path in get_run_axis_csvs(wildcards.run_id).items()
    )


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


def get_fdx_model():
    """Path to the FDX (5' index) CRF bundle directory, resolved like get_ldx_model.

    Same shape as the LDX bundle -- `metadata.json` plus the ONNX graph it
    names -- but anchored on the read end rather than on the boundary detector's
    `adapter_end`, which is why it takes none of the boundary flags and why it
    is demultiplexed in its own pass. See the `fdx` block in config-base.yml.
    """
    model = config.get("fdx", {}).get("model")
    if not model:
        sys.exit(
            "fdx.enabled is true but fdx.model is unset. Point it at a CRF "
            "bundle directory, e.g. "
            "resources/models/demux/barcode_crf_fdx4_rna004@v0.2.0"
        )
    if not os.path.isabs(model):
        model = os.path.join(PIPELINE_DIR, model)
    if not os.path.isdir(model):
        sys.exit(
            f"fdx.model is not a directory: {model}\n"
            "escapepod CRF models are self-describing BUNDLES (metadata.json "
            "plus the ONNX graphs it names), not a single file."
        )
    return model


def validate_escpod_barcodes():
    """Fail at DAG time if a sample names a code its axis's bundle does not emit.

    The alternative is a name that appears nowhere in the model's output, which
    surfaces as "no reads were assigned" only AFTER the demux pass (and, for
    the LDX axis, the whole-run basecall) have been paid for.
    """
    for sample, info in samples.items():
        if info.get("barcode") and sample_is_ldx(sample):
            get_sample_axis_codes(sample)


def escpod_device_args(boundary_anchored=True):
    """Device placement for an escpod demux pass.

    `--device gpu` runs the CRF encoder (and, for a boundary-anchored model,
    the boundary CNN) through onnxruntime's CUDA provider; the lattice decode
    stays on the CPU either way.

    Passed EXPLICITLY in both directions rather than leaning on the `auto`
    default, because the two failure modes are not symmetric. `--device gpu` is
    a requirement: it fails when the feature is missing, no device is visible,
    or onnxruntime cannot register its CUDA provider. `auto` would quietly run
    on the CPU instead — a 20x slowdown (2.6 h against 458.9 s on our own
    flowcell) that still produces correct output, so it looks like success and
    is caught only by noticing the wall clock. Under `ldx.gpu: false` the
    explicit `--device cpu` likewise stops a GPU binary from opportunistically
    using a device the config said not to.

    Replaces `--gpu`, deprecated in escapepod-rs 0.17.1. That spelling still
    runs (it warns and continues as `--device gpu`), so this is a migration
    rather than a break — but the meaning changed underneath it: `--gpu` fell
    back to the CPU where `--device gpu` fails.

    `--method cnn` is passed alongside on purpose for a boundary-anchored
    bundle. The bundle already pins cnn, and naming the same detector the
    bundle pins is not an override, so this cannot silently downgrade to LLR.
    A read-end-anchored bundle (fdx) consumes no detector and REFUSES the flag.
    """
    if not config.get("ldx", {}).get("gpu", False):
        return "--device cpu"
    return "--device gpu --method cnn" if boundary_anchored else "--device gpu"


def escpod_ort_env():
    """Shell prefix that lets a GPU escpod find its CUDA onnxruntime.

    ort dlopens onnxruntime at run time from ORT_DYLIB_PATH, and the CUDA
    execution provider sits beside it, so its directory must also be on
    LD_LIBRARY_PATH or the core library loads and the provider then fails to.
    Exported inside the rule rather than relying on the submitting shell's
    environment surviving into the batch job.
    """
    if not config.get("ldx", {}).get("gpu", False):
        return ""
    return (
        f"export ORT_DYLIB_PATH={get_ort_dylib()}; "
        f"export LD_LIBRARY_PATH={os.path.dirname(get_ort_dylib())}:"
        f"{os.path.join(PIPELINE_DIR, '.pixi', 'envs', 'gpu', 'lib')}:"
        f"$LD_LIBRARY_PATH; "
    )


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
        # `--model <ldx>` alone, or `--model ldx=<ldx> --model fdx=<fdx>` for a
        # fused dual-index pass (fdx.fused); see get_escpod_model_args.
        models=get_escpod_model_args(),
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
        # Device placement and the CUDA runtime it needs; see the two helpers
        # for why both directions are spelled out and why `--method cnn` rides
        # along on a boundary-anchored pass but is refused once the read-end
        # fdx model is in the same run.
        gpu=lambda wildcards: escpod_device_args(boundary_anchored=not is_fdx_fused()),
        ort_env=lambda wildcards: escpod_ort_env(),
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
            {params.models} \
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


if is_fdx_enabled() and is_fdx_fused():

    rule summarize_fdx_axis:
        """
        The per-code tally of the fdx axis of a FUSED pass.

        With `fdx.fused`, escapepod_demux calls both axes in one sweep and its
        CSV carries them as per-axis columns; this pulls the fdx column into the
        same summary shape the two-pass rule writes, at the same path, so the
        report and the tests read one layout whichever way the run was demuxed.
        """
        input:
            classifications=os.path.join(
                outdir, "demux", "read_ids", "{run_id}", "classifications.csv"
            ),
        output:
            summary=os.path.join(
                outdir, "demux", "read_ids", "{run_id}", "fdx", "demux_summary.tsv.gz"
            ),
        log:
            os.path.join(outdir, "logs", "summarize_fdx_axis", "{run_id}"),
        params:
            summarize_awk=os.path.join(SCRIPT_DIR, "summarize_demux.awk"),
        shell:
            """
            mkdir -p $(dirname {output.summary})
            awk -F, -v col=fdx -f {params.summarize_awk} {input.classifications} \
                2>{log} | gzip >{output.summary}
            """


if is_fdx_enabled() and not is_fdx_fused():

    rule escapepod_demux_fdx:
        """
        Call the 5' FDX index on every read of a run: the second escpod axis,
        as its own pass over the raw POD5.

        escpod can call both axes in ONE sweep (`--model ldx=... --model fdx=...`,
        `fdx.fused: true`), and that is the shape to want on a 600 GB flowcell
        where the sweep is hours of IO. This rule exists because escpod (0.19.0 and
        0.20.0, both measured) cannot yet run that pass the way the LDX axis needs it: in a fused run
        every model shares one `--boundary-margin` / `--clamp-max-shift`, and
        escpod refuses those flags outright when any axis anchors on the read
        end -- which this bundle does ("--boundary-margin is not applicable:
        this model anchors its window on the read end"). boundary_margin 0 is
        worth ~14% of a run's LDX reads, so until upstream scopes the flags to
        the axis that has a boundary, each axis runs on its own terms and
        select_demux_reads.py joins them per read. The join is the same in both
        shapes, so flipping `fdx.fused` changes nothing downstream.

        `--model fdx=` names the axis, which is what the sidecar column is called
        (the LDX pass, unnamed, writes `barcode`), so both calls live in the one
        collection sidecar beside the raw POD5. The lattice columns (`crf_*`) are
        shared and belong to whichever pass ran last; the per-axis classifications
        CSVs are the record the pipeline reads. The LDX pass's CSV is an input
        here for ORDERING only: both passes rewrite that sidecar, and two
        concurrent rewrites would lose one axis. No `--method`: this bundle
        consumes no boundary detector and rejects the flag. Its own lattice gate
        is passed here AND re-applied in the join, so a fused run gates the same
        way.
        """
        input:
            pod5=get_run_raw_inputs,
            ldx_done=os.path.join(
                outdir, "demux", "read_ids", "{run_id}", "classifications.csv"
            ),
        output:
            classifications=os.path.join(
                outdir, "demux", "read_ids", "{run_id}", "fdx", "classifications.csv"
            ),
            summary=os.path.join(
                outdir, "demux", "read_ids", "{run_id}", "fdx", "demux_summary.tsv.gz"
            ),
        log:
            os.path.join(outdir, "logs", "escapepod_demux_fdx", "{run_id}"),
        threads: config.get("ldx", {}).get("threads", 16)
        params:
            model=get_fdx_model(),
            escpod=get_escpod_bin(),
            # This axis's own gate; the bundle's declared operating point by default.
            # Not optional in practice: a CRF snaps a read carrying none of its codes
            # onto the nearest one, so an ungated FDX call on a read with no 5'
            # index is a confident wrong answer. See config-base.yml `fdx`.
            min_crf_margin=lambda wildcards: (
                f"--min-crf-margin {config['fdx']['min_crf_margin']}"
                if config.get("fdx", {}).get("min_crf_margin") is not None
                else ""
            ),
            min_crf_prob=lambda wildcards: (
                f"--min-crf-prob {config['fdx']['min_crf_prob']}"
                if config.get("fdx", {}).get("min_crf_prob") is not None
                else ""
            ),
            gpu=lambda wildcards: escpod_device_args(boundary_anchored=False),
            ort_env=lambda wildcards: escpod_ort_env(),
            pod5_dirs=lambda wildcards: " ".join(get_run_pod5_dirs(wildcards.run_id)),
            summarize_awk=os.path.join(SCRIPT_DIR, "summarize_demux.awk"),
        shell:
            """
            mkdir -p $(dirname {output.classifications})

            {params.ort_env}{params.escpod} demux {params.pod5_dirs} \
                --model fdx={params.model} \
                --annotate \
                --classifications {output.classifications} \
                --ref-scores \
                {params.min_crf_margin} \
                {params.min_crf_prob} \
                {params.gpu} \
                --threads {threads} 2>&1 | tee {log}

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


def get_sample_run_split_parents(wildcards):
    """Locate the split-read parent map for a sample's run."""
    run_id = samples[wildcards.sample]["run_id"]
    return os.path.join(outdir, "demux", "read_ids", run_id, "split_parents.tsv")


def get_sample_run_ubam(wildcards):
    """Locate the run-level basecall for a sample's run."""
    run_id = samples[wildcards.sample]["run_id"]
    return os.path.join(outdir, "bam", "rebasecall_run", run_id, f"{run_id}.rbc.bam")


rule ldx_run_read_ids:
    """
    The reads demux assigned to ANY of this run's samples, across every axis.

    This is the `-l` list rebasecall_ldx_run hands dorado, so it bounds the GPU
    work to reads some sample will actually own: unclassified reads, reads of a
    code no sample claims, and on a dual-index run reads whose two calls do not
    make a configured pair are never basecalled. It used to be an awk over the
    `barcode` column inside the basecall rule; with a second axis to join it is
    the same script the per-sample selection uses, run over all samples at once,
    so the two cannot disagree about what "assigned" means. Every sample's count
    is recorded beside it, which is the per-sample view of where the join lost
    reads (the LDX pass's demux_summary is the per-code one).
    """
    input:
        csvs=lambda wildcards: sorted(set(get_run_axis_csvs(wildcards.run_id).values())),
    output:
        read_ids=maybe_temp(
            os.path.join(
                outdir, "demux", "read_ids", "{run_id}", "assigned_read_ids.txt"
            ),
            tier="demux_scratch",
        ),
        summary=os.path.join(
            outdir, "demux", "read_ids", "{run_id}", "assigned_summary.tsv"
        ),
    log:
        os.path.join(outdir, "logs", "ldx_run_read_ids", "{run_id}"),
    params:
        src=SCRIPT_DIR,
        axes=get_run_axis_csv_args,
        gates=get_axis_gate_args(),
        samples=lambda wildcards: " ".join(
            f"--sample {get_sample_axis_args(s)}"
            for s in get_samples_for_run(wildcards.run_id)
        ),
        downstream=get_downstream_axis_args(),
    shell:
        """
        python {params.src}/select_demux_reads.py \
            {params.axes} \
            {params.gates} \
            {params.downstream} \
            {params.samples} \
            --output {output.read_ids} \
            --summary {output.summary} \
            2>&1 | tee {log}
        """


rule rebasecall_ldx_run:
    """
    Basecall a whole LDX run in one dorado pass, straight off the raw POD5.

    The WarpDemuX path basecalls each sample's split POD5 separately because that
    split is what its `escpod filter` step produced. Here there is no split POD5 —
    demux left only a sidecar — so the run is basecalled whole and cut up
    afterwards, on the uBAM.

    `-l` restricts the pass to reads that demux actually assigned to one of this
    run's samples (ldx_run_read_ids). That is the same read universe the
    per-barcode POD5s used to cover: unclassified reads were never in any of
    them, and basecalling them here would be GPU time spent on output nothing
    reads.
    """
    input:
        pod5=get_run_raw_inputs,
        read_ids=rules.ldx_run_read_ids.output.read_ids,
        mod_models=rules.download_mod_models.output.sentinel,
    output:
        bam=maybe_temp(
            os.path.join(
                outdir, "bam", "rebasecall_run", "{run_id}", "{run_id}.rbc.bam"
            ),
            tier="basecall",
        ),
    log:
        os.path.join(outdir, "logs", "rebasecall_ldx_run", "{run_id}"),
    params:
        model=config["base_calling_model"],
        dorado_opts=config["opts"]["dorado"],
        models_dir=os.path.join(PIPELINE_DIR, "resources", "models"),
        run_path=lambda wildcards: get_run_path(wildcards.run_id),
        resume_sh=os.path.join(SCRIPT_DIR, "dorado_basecall_resume.sh"),
    shell:
        """
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
            --read-ids {input.read_ids}
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
    List the reads belonging to one sample: the reads every axis the sample
    names agrees on, plus any split children of those reads.

    One axis (`ldx`) on a single-index library; `ldx` and `fdx` on a dual-index
    one, where the two calls come from two passes and are joined here by read
    id. A split child inherits its parent's assignment: the parent id itself is
    not in the BAM — dorado replaces the read with its children — so keeping it
    in the list is harmless and keeps this file a faithful record of what demux
    assigned. The codes are passed as the MODEL emits them (see
    get_sample_axis_codes), and the script fails loudly when a sample ends up
    with no reads at all.
    """
    input:
        csvs=get_sample_axis_csvs,
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
        src=SCRIPT_DIR,
        axes=get_sample_axis_csv_args,
        gates=get_axis_gate_args(),
        sample=lambda wildcards: get_sample_axis_args(wildcards.sample),
        downstream=get_downstream_axis_args(),
    shell:
        """
        python {params.src}/select_demux_reads.py \
            {params.axes} \
            {params.gates} \
            {params.downstream} \
            --sample {params.sample} \
            --parents {input.parents} \
            --output {output.read_ids} \
            2>&1 | tee {log}
        """


rule split_ldx_ubam:
    """
    Cut one sample's unaligned BAM out of the run-level basecall.

    Emits the path `rebasecall` would have, so everything downstream — EDX adapter
    detection, alignment, classification — is unchanged and cannot tell which
    backend produced it.
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

    # Every configured code against the bundle that has to emit it, now rather
    # than after the demux pass and the whole-run basecall have run.
    validate_escpod_barcodes()

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

    This list is the whole EDX filter: bwa_align aligns only the reads it names
    (`samtools view -N`), and the classifier only ever touches reads the BAM
    names, so neither a filtered FASTQ nor a filtered POD5 is written any more.
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
