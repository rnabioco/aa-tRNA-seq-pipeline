import os
import glob
import re
import sys
import pysam
import yaml
from git import Repo

SCRIPT_DIR = os.path.join(SNAKEFILE_DIR, "scripts")

# The barcode-name crosswalk. Every escpod CRF bundle emits whatever its
# metadata.json calls its references. For the shipped ldx16 panel that is
# already this project's vocabulary (`ldx01`); for WDX4 it is not (`bc03`
# against a configured `barcode03`). Resolved against the bundle rather than
# guessed; see workflow/scripts/barcode_names.py.
sys.path.insert(0, SCRIPT_DIR)
from barcode_names import emitted_to_label, resolve_to_bundle

# Whether this run's basecaller is the one the charging bundle was trained
# against. A charging model reads the k-mer residual against a level predicted
# from the read's own basecall, so the basecaller is part of the feature
# definition rather than tooling; see workflow/scripts/basecaller_compat.py.
from basecaller_compat import check_basecaller

# Cleanup tiers for maybe_temp(). Each large intermediate is assigned a tier so
# they can be deleted or kept independently via the `cleanup_intermediates`
# config key (see maybe_temp / _enabled_cleanup_tiers below).
_CLEANUP_TIERS = {
    "cascade",  # bam/aln, calmd, charging, adapter_tagged (redundant near-copies; bam/final hardlinks the last)
    "basecall",  # bam/rebasecall, bam/rebasecall_run (GPU-hours to regenerate)
    "demux_scratch",  # demux/warpdemux_output, demux/read_ids, edx read_ids
    "split_pod5",  # demux/pod5 (WDX split POD5: that path's signal store and classification input)
}

# Tiers that existed before v0.7 and no longer name anything. `fastq` went when
# alignment started streaming straight from the uBAM (no fq/ is written);
# `merged_pod5` went when stage_pod5 replaced merge_pods (pod5/ now holds
# symlinks, not a copy of the run). Accepted and ignored so an existing config
# keeps working, but named, so the reader learns why the space it expected to
# reclaim is no longer being spent.
_RETIRED_CLEANUP_TIERS = {"fastq", "merged_pod5"}


def _enabled_cleanup_tiers():
    """Resolve `cleanup_intermediates` config into a set of enabled tier names.

    Accepts a bool or a list of tier names:
      - True         -> all tiers enabled
      - False/absent -> no tiers (opt-in default; nothing auto-deleted)
      - list         -> only the named tiers
    """
    cfg = config.get("cleanup_intermediates", False)
    if cfg is True:
        return set(_CLEANUP_TIERS)
    if not cfg:
        return set()
    unknown = set(cfg) - _CLEANUP_TIERS - _RETIRED_CLEANUP_TIERS
    if unknown:
        sys.exit(
            f"cleanup_intermediates names unknown tier(s) {sorted(unknown)}; "
            f"choose from {sorted(_CLEANUP_TIERS)}"
        )
    return set(cfg) & _CLEANUP_TIERS


def maybe_temp(path, tier="cascade"):
    """Mark path as temp() only when its cleanup tier is enabled."""
    return temp(path) if tier in _enabled_cleanup_tiers() else path


def get_charging_model():
    """Absolute path to the charging model BUNDLE directory.

    Relative paths resolve against the pipeline directory, not the invocation
    cwd, so a run launched from elsewhere still finds the vendored bundle.
    Mirrors get_ldx_model() in demux.smk — same failure mode, since both are
    self-describing directories rather than single model files.
    """
    model = config.get("charging", {}).get("model")
    if not model:
        sys.exit(
            "charging.model is unset. Point it at a charging bundle "
            "directory, e.g. "
            "resources/models/charging/charging_feature_nn_sup6_rna004@v0.1.0"
        )
    if not os.path.isabs(model):
        model = os.path.join(PIPELINE_DIR, model)
    if not os.path.isdir(model):
        sys.exit(
            f"charging.model is not a directory: {model}\n"
            "escapepod charging models are self-describing BUNDLES "
            "(metadata.json plus the ONNX graph and k-mer table it names), "
            "not a single file. The Remora .pt classifier this replaced was a "
            "file; a config carried over from before v0.4.0 will still point "
            "at one."
        )
    return model


# escpod's own accepted values for `classify --orientation`.
_ORIENTATIONS = ("auto", "time", "reversed")


def get_charging_orientation_arg():
    """The `--orientation` flag for `escpod classify`, or "" to leave it alone.

    `auto` is escpod's default and passing it explicitly would say nothing, so
    the default emits no flag at all: the command line stays exactly what it was
    before this knob existed, and a run that never sets it is unaffected.

    Forcing a frame is a correctness decision, not a convenience. `auto` errors
    below 50 informative reads, which on a heavily multiplexed run fails an
    otherwise fine sample; a forced frame that is WRONG does not error at all --
    it mis-anchors every feature and the calls become noise. Validated here so a
    typo stops the DAG rather than 448 jobs.

    The frame is a property of the RUN (chemistry plus basecaller), not of a
    sample, so the value belongs to one run and must be re-derived for the next.
    Every sample deep enough for `auto` logs the frame it found and its vote
    counts, which is where the value comes from; see config-base.yml.
    """
    value = config.get("charging", {}).get("orientation", "auto")
    if value not in _ORIENTATIONS:
        sys.exit(
            f"charging.orientation is {value!r}; expected one of "
            f"{', '.join(_ORIENTATIONS)}. `auto` detects the frame per sample "
            "and is right unless a sample is too shallow for it (fewer than 50 "
            "informative reads); force a frame only from a value `auto` agreed "
            "on for the same chemistry and basecaller on a deeper run."
        )
    return "" if value == "auto" else f"--orientation {value}"


def get_charging_orientation_fallback():
    """The frame to supply when `auto` is underpowered, or "" for none.

    Empty whenever a frame is already forced -- there is nothing to fall back
    from -- and whenever the fallback is switched off, which restores v0.7.2's
    behaviour of failing an underpowered sample.

    This is a fallback rather than a forced default because detection is worth
    keeping wherever it works. `base_calling_model` and
    `charging.basecaller_check` already pin and enforce the two things the frame
    depends on, so it is effectively a constant -- but letting escpod reach its
    own consensus on every deep sample keeps a free tripwire for the case where
    that assumption stops holding, and forcing everywhere would silence it.
    """
    charging = config.get("charging", {})
    if charging.get("orientation", "auto") != "auto":
        return ""
    value = charging.get("orientation_fallback", "reversed")
    if value in (None, "none", False):
        return ""
    if value not in ("time", "reversed"):
        sys.exit(
            f"charging.orientation_fallback is {value!r}; expected time, "
            "reversed, or none. It supplies the frame for a sample too thin for "
            "`auto` to detect one (fewer than 50 informative reads) instead of "
            "failing it; `none` restores the failure."
        )
    return value


def get_charging_device_arg():
    """`--device cpu`/`--device gpu` for `escpod classify`, passed explicitly.

    Same reasoning as `escpod_device_args()` in demux.smk: `--device gpu` is a
    requirement, so it fails loudly if the feature or a CUDA device is
    missing, where `auto` would silently keep running on the CPU and look like
    a normal, if slow, success. `--device cpu` under `charging.gpu: false`
    likewise stops an opportunistically-GPU-capable binary from using a device
    the config said not to.

    Only the windowed (TCN) charging bundle has a GPU path; the GBM/feature-
    network bundles ignore this flag with a one-line log note
    (`note_cpu_only`) rather than an error, so it is safe to always pass one
    regardless of which bundle `charging.model` names.
    """
    return "--device gpu" if config["charging"].get("gpu", False) else "--device cpu"


def get_charging_escpod_gpu_prefix():
    """Shell prefix putting a GPU-enabled escpod ahead of the default build on
    PATH for one command, or "" under `charging.gpu: false`.

    `escpod_classify_fallback.sh` calls plain `escpod`, which the Snakefile's
    `onstart` prefix always resolves to the portable (CPU) musl build. GPU
    classify needs the dynamically-linked `-gpu` release artifact instead --
    the same one `ldx.gpu` already downloads via `pixi run setup` into
    `<escpod_version>-gpu/` (see `get_escpod_bin` in demux.smk) -- so this
    shadows it onto PATH for just this invocation rather than moving the
    global pin, for the same reason `get_escpod_bin` gives: the GPU artifact
    is x86_64 Linux only and dynamically linked, and every other rule should
    keep the portable one.

    Unlike `get_escpod_bin`, this does not also need to check for a vendored
    CUDA libonnxruntime or cuDNN: `escpod classify`'s GPU path is
    tract-cuda/cudarc, which dlopens the CUDA driver and compiles its own
    kernels at run time (NVRTC) against the node's own CUDA toolkit, not a
    vendored onnxruntime.
    """
    if not config["charging"].get("gpu", False):
        return ""
    version = config.get("escpod_version", ESCPOD_VERSION)
    bin_dir = os.path.join(
        PIPELINE_DIR, "resources", "tools", "escpod", f"{version}-gpu", "bin"
    )
    if not os.path.isfile(os.path.join(bin_dir, "escpod")):
        sys.exit(
            f"charging.gpu is true but no GPU-enabled escpod was found at "
            f"{bin_dir}.\n"
            "Install it with `pixi run setup`, which downloads the published "
            f"GPU artifact for escpod {version} (x86_64 Linux only).\n"
            "A visible CUDA device is also required at run time; unlike "
            "ldx.gpu, no separate onnxruntime or cuDNN install is needed.\n"
            "Set charging.gpu: false to run classify on the CPU instead."
        )
    return f'export PATH="{bin_dir}:$PATH"; '


def is_warpdemux_enabled():
    """Check if WarpDemuX (WDX) demultiplexing is enabled in config."""
    return config.get("warpdemux", {}).get("enabled", False)


def is_ldx_enabled():
    """Check if escapepod CRF (LDX) demultiplexing is enabled in config."""
    return config.get("ldx", {}).get("enabled", False)


def is_edx_enabled():
    """Whether the 3' adapter (EDX) is in use as a sample-distinguishing axis.

    Unlike LDX and FDX this is NOT a signal axis: the adapter is read off the
    uBAM by detect_edx_adapters, after basecalling. It still separates samples,
    so anything reasoning about whether two samples can collide has to account
    for it -- see get_sample_downstream_codes in demux.smk.
    """
    return config.get("edx", {}).get("enabled", False)


def is_fdx_enabled():
    """Whether the 5' FDX index is demultiplexed as a second escpod axis.

    FDX is a second barcode on the SAME molecule as the LDX one -- a 5' index
    read off the read end, where LDX is a 3' index read off the adapter
    boundary -- so it is an additional axis of the escpod path, not a backend
    of its own. It requires `ldx.enabled`: a library carrying only a 5' index
    is not a case this pipeline has seen yet, and the LDX pass is what the
    read-level accounting (demux_summary, read_attrition) is keyed on.
    """
    return config.get("fdx", {}).get("enabled", False)


def is_fdx_fused():
    """Whether both escpod axes are called in ONE pass over the raw POD5.

    `escpod demux --model ldx=... --model fdx=...` decodes each read's signal
    once and calls both axes, writing one classifications CSV with per-axis
    columns. It is the shape upstream recommends and the one to want on a
    600 GB flowcell, where the POD5 sweep is hours of IO-bound wall. It is off
    by default. escpod 0.19.0 and 0.20.0 (both measured) refused
    `--boundary-margin` / `--clamp-max-shift` whenever a model in the run
    anchored on the read end, which the fdx bundle does, and the LDX axis
    cannot give those flags up. escpod 0.21.0 scopes the flags per head
    (escapepod-rs#323), so the fused path is legal on the pinned version; it
    stays off until someone measures it. See the `fdx` block in
    config-base.yml.
    """
    return is_fdx_enabled() and config.get("fdx", {}).get("fused", False)


def get_sample_barcode_label(sample):
    """The project-facing barcode name for a sample, or None if it has none.

    `samples[s]["barcode"]` holds whatever the samples YAML said. Both live
    panels now configure a name this project already owns — WarpDemuX barcodes
    are ours (`barcode04`), and the ldx16 bundle's own references are
    `ldx01`..`ldx16` — so this is the identity on every current run. It stays
    because the emitted vocabulary is a property of the bundle, not a constant:
    the retired nbc16 panel named the same physical codes `nbcNN`, and a future
    bundle may differ again.

    Each mapping is a pure documented rename (bcNN == barcodeNN; see
    config/README.md and resources/models/demux/README.md), so nothing is lost.
    When one applies, the upstream name is also recorded in an @CO line on the
    BAM, so a reader never has to guess which naming a file is using.

    A dual-index sample gets both codes joined with `-`, the SAM convention for
    a dual index in `BC` (`ldx01-fdx01`), 3' code first because that is the
    axis the run is keyed on.
    """
    barcode = samples.get(sample, {}).get("barcode")
    if not barcode:
        return None
    # Driven by the name itself rather than by which backend is enabled: the
    # WDX4 panel is served by a CRF bundle too, and it emits `bc03` where the
    # samples file says `barcode03`, so the rename is a property of the
    # vocabulary and not of the config. `emitted_to_label` is the identity on
    # every name that is already project-facing, which since the ldx16 switch
    # includes every LDX one, and every FDX one.
    label = emitted_to_label(barcode)
    fdx = samples[sample].get("fdx")
    return f"{label}-{emitted_to_label(fdx)}" if fdx else label


def get_sample_barcode_upstream(sample):
    """The raw barcode name(s) as configured, for provenance. None if unbarcoded."""
    barcode = samples.get(sample, {}).get("barcode")
    if not barcode:
        return None
    fdx = samples[sample].get("fdx")
    return f"{barcode}-{fdx}" if fdx else barcode


def is_demux_enabled():
    """Check if any signal-level barcode demultiplexing backend is enabled.

    The two backends are mutually exclusive: they populate the same per-sample
    `barcode` field and their rules write the same barcode_mapping output, so
    enabling both would make the DAG ambiguous.
    """
    if is_warpdemux_enabled() and is_ldx_enabled():
        sys.exit(
            "Config enables both `warpdemux` and `ldx` demultiplexing. "
            "These are alternative backends for the same step — enable exactly one."
        )
    if is_fdx_enabled() and not is_ldx_enabled():
        sys.exit(
            "Config enables `fdx` without `ldx`. The 5' FDX index is a second axis "
            "of the escpod demux path and is joined against the 3' LDX call per "
            "read; a library carrying only an FDX index is not supported yet. "
            "Enable `ldx` as well, or disable `fdx`."
        )
    return is_warpdemux_enabled() or is_ldx_enabled()


def parse_samples_tsv(fl):
    """Parse TSV-format sample file (existing format)."""
    samples = {}
    with open(fl) as f:
        for l in f:
            line = l.rstrip()
            if not line or line.startswith("#"):
                continue
            try:
                sample, path = line.split()
            except:
                print(
                    "samples file must have 2 columns, sample_id and data_path, separated by whitespace",
                    file=sys.stderr,
                )
                sys.exit(f"found {line}")
            if sample in samples:
                samples[sample]["path"].add(path)
            else:
                samples[sample] = {
                    "path": {path},
                    "barcode": None,
                    "fdx": None,
                    "edx": None,
                    "run_id": None,
                }
    return samples


def parse_samples_yaml(fl):
    """
    Parse YAML-format sample file with barcode assignments.

    YAML format:
    runs:
      - path: /path/to/pooled/run
        barcode_kit: "WDX4_rna004_v1_0"  # optional, uses config default
        samples:
          sample_name: "barcode04"
          another_sample: "barcode05"

      - path: /path/to/dual-index/run
        samples:
          lib_a_rep1: {ldx: "ldx01", fdx: "fdx01"}   # 3' LDX code + 5' FDX code
          lib_a_rep2: {ldx: "ldx02", fdx: "fdx01"}

      - path: /path/to/non-demux/run
        samples:
          direct_sample: ~  # null barcode = no demux

    `fdx:` needs `fdx.enabled` in the config, and within one run every sample
    sharing an `ldx:` code must either all name an `fdx:` or none of them: a
    sample named by `ldx01` alone would otherwise swallow the reads of one named
    by `ldx01` + `fdx01`. Tuples must be unique per run.
    """
    samples = {}
    with open(fl) as f:
        data = yaml.safe_load(f)

    if "runs" not in data:
        sys.exit(f"YAML samples file must contain 'runs' key: {fl}")

    for run_idx, run in enumerate(data["runs"]):
        if "path" not in run:
            sys.exit(f"Run {run_idx} missing 'path' in samples file: {fl}")
        if "samples" not in run:
            sys.exit(f"Run {run_idx} missing 'samples' in samples file: {fl}")

        run_path = run["path"]
        # Create a unique run_id from the path (last directory component)
        run_id = os.path.basename(run_path.rstrip("/"))
        barcode_kit = run.get(
            "barcode_kit", config.get("warpdemux", {}).get("barcode_kit")
        )

        for sample_name, sample_val in run["samples"].items():
            if sample_name in samples:
                sys.exit(f"Duplicate sample name '{sample_name}' in samples file: {fl}")

            # Backward compat: plain string or null = signal barcode only
            if isinstance(sample_val, str) or sample_val is None:
                barcode = sample_val
                fdx = None
                edx = None
            elif isinstance(sample_val, dict):
                # `wdx` (WarpDemuX) and `ldx` (escapepod CRF) both name the
                # signal-level barcode; which one is meaningful depends on the
                # enabled backend, so only one may be given per sample.
                wdx_bc, ldx_bc = sample_val.get("wdx"), sample_val.get("ldx")
                if wdx_bc is not None and ldx_bc is not None:
                    sys.exit(
                        f"Sample '{sample_name}' sets both 'wdx' and 'ldx'. "
                        "These name the same field for different demux backends — "
                        "give exactly one."
                    )
                barcode = wdx_bc if wdx_bc is not None else ldx_bc
                fdx = sample_val.get("fdx")
                edx = sample_val.get("edx")
                if fdx is not None and not is_fdx_enabled():
                    sys.exit(
                        f"Sample '{sample_name}' names an `fdx:` code but `fdx.enabled` "
                        "is off. Enable the FDX axis in the config, or drop the key."
                    )
                if fdx is not None and ldx_bc is None:
                    sys.exit(
                        f"Sample '{sample_name}' names an `fdx:` code without an "
                        "`ldx:` one. The FDX axis is joined against the LDX call per "
                        "read, so a dual-index sample needs both."
                    )
            else:
                sys.exit(f"Invalid sample value for '{sample_name}': {sample_val}")

            samples[sample_name] = {
                "path": {run_path},
                "barcode": barcode,
                "fdx": fdx,
                "edx": edx,
                "run_id": run_id,
                "barcode_kit": barcode_kit,
            }

    _check_barcode_tuples(samples, fl)
    return samples


def _check_barcode_tuples(samples, fl):
    """Refuse sample tuples that would overlap on a run.

    The per-sample read selection (select_demux_reads.py) assigns a read to a
    sample when every axis the sample names agrees. A sample naming `ldx01`
    alone therefore contains every read of a sample naming `ldx01` + `fdx01`,
    and two samples with identical tuples receive identical reads. Either is a
    mistake in the samples file, caught here rather than in a BAM.

    `edx` counts as an axis here even though it is applied later and elsewhere
    -- 3' adapter identity on the uBAM, not a signal-level call -- because it
    still distinguishes the reads a sample ends up with. One LDX code fanned
    across several EDX adapters is a normal design, so leaving `edx` out of the
    tuple refuses a legitimate samples file (#161).
    """
    by_run = {}
    for name, info in samples.items():
        if not info.get("barcode"):
            continue
        by_run.setdefault(info["run_id"], {}).setdefault(info["barcode"], []).append(
            name
        )
    for run_id, by_code in by_run.items():
        for code, names in by_code.items():
            if len(names) < 2:
                continue
            with_fdx = [n for n in names if samples[n].get("fdx")]
            if with_fdx and len(with_fdx) != len(names):
                sys.exit(
                    f"Samples {', '.join(sorted(names))} share barcode {code} on run "
                    f"{run_id}, but only {', '.join(sorted(with_fdx))} name an `fdx:` "
                    "code. A sample named by the LDX code alone would swallow the "
                    "others' reads; give every one of them an `fdx:`, or none. "
                    f"({fl})"
                )
            with_edx = [n for n in names if samples[n].get("edx")]
            if with_edx and len(with_edx) != len(names):
                sys.exit(
                    f"Samples {', '.join(sorted(names))} share barcode {code} on run "
                    f"{run_id}, but only {', '.join(sorted(with_edx))} name an `edx:` "
                    "code. EDX filtering is applied per sample, so one named by the "
                    "LDX code alone is never filtered and swallows the others' reads; "
                    f"give every one of them an `edx:`, or none. ({fl})"
                )
            tuples = [
                (code, samples[n].get("fdx"), samples[n].get("edx")) for n in names
            ]
            if len(set(tuples)) != len(tuples):
                sys.exit(
                    f"Samples {', '.join(sorted(names))} on run {run_id} have identical "
                    f"barcode assignments ({fl})."
                )


def parse_samples(fl):
    """
    Parse sample file, detecting format based on extension.
    .yml/.yaml files use YAML format with barcode support.
    Other files use TSV format (backward compatible).
    """
    if fl.endswith(".yml") or fl.endswith(".yaml"):
        return parse_samples_yaml(fl)
    else:
        return parse_samples_tsv(fl)


def get_pipeline_commit():
    try:
        repo = Repo(PIPELINE_DIR)
        return repo.head.commit
    except Exception:
        return None


def get_pipeline_version():
    """
    Get comprehensive pipeline version information from git.

    Returns dict with commit, tag, branch, and dirty status.
    """
    try:
        repo = Repo(PIPELINE_DIR)
        commit = repo.head.commit

        # Find tags pointing to current commit
        tags = [tag.name for tag in repo.tags if tag.commit == commit]

        # Get branch name (None if detached HEAD)
        try:
            branch = repo.active_branch.name
        except TypeError:
            branch = None

        return {
            "git_commit": str(commit),
            "git_tag": tags[0] if tags else None,
            "git_branch": branch,
            "git_dirty": repo.is_dirty(),
        }
    except Exception:
        return {
            "git_commit": None,
            "git_tag": None,
            "git_branch": None,
            "git_dirty": None,
        }


def format_config_values():
    x = []
    x.append("Config settings:")
    for k, v in config.items():
        if k == "opts":
            x.append(f"\t{k}:")
            for cmd, opts in v.items():
                x.append(f"\t\t{cmd}: {opts}")
        else:
            x.append(f"\t{k}: {v}")
    return "\n".join(x)


def report_metadata():
    from snakemake.logging import logger

    cid = get_pipeline_commit()
    if cid is not None:
        logger.info(f"Pipeline commit: {cid}")
    else:
        logger.warning("Pipeline commit: unable to resolve git commit")
    logger.info(format_config_values())


def find_raw_inputs(sample_dict):
    """
    parse through directories listed in samples.tsv and identify pod5 files to process
    store input files and uuid base file names in dictionary for each sample
    """
    POD5_DIRS = ["pod5_pass", "pod5_fail", "pod5"]
    ext = ".pod5"

    for sample, info in sample_dict.items():
        raw_fls = []
        for path in info["path"]:
            for subdir in POD5_DIRS:
                data_path = os.path.join(path, subdir, "*" + ext)
                fls = glob.glob(data_path)
                raw_fls += fls
        if len(raw_fls) == 0:
            sys.exit(
                f"No input files found for sample: {sample}. Please check the path in the samples.tsv file"
            )
        sample_dict[sample]["raw_files"] = raw_fls

    return sample_dict


# set up global samples dictionary to be used throughout pipeline
outdir = config["output_directory"]

_reuse_from = config.get("reuse_outputs_from")
if _reuse_from:
    _reuse_from = os.path.realpath(_reuse_from)
    # Guard: source must exist
    if not os.path.isdir(_reuse_from):
        sys.exit(f"reuse_outputs_from: directory not found: {_reuse_from}")
    # Guard: must not be same as outdir
    os.makedirs(outdir, exist_ok=True)
    if os.path.realpath(_reuse_from) == os.path.realpath(outdir):
        sys.exit("reuse_outputs_from cannot be the same as output_directory")

    # `pod5/` from a pre-v0.7 run holds a merged copy per sample where this
    # version stages symlinks; both are a directory of POD5 the tools walk, so
    # either shape serves. (No `fq/`: alignment streams from the uBAM now.)
    _REUSE_DIRS = ["pod5", "demux", "bam/rebasecall"]
    for _subdir in _REUSE_DIRS:
        _src = os.path.join(_reuse_from, _subdir)
        _dst = os.path.join(os.path.realpath(outdir), _subdir)
        if not os.path.isdir(_src):
            continue  # skip missing (e.g., demux/ when demux disabled)
        os.makedirs(os.path.dirname(_dst), exist_ok=True)
        if os.path.exists(_dst):
            if os.path.islink(_dst) and os.path.realpath(_dst) == os.path.realpath(
                _src
            ):
                continue  # already linked correctly
            sys.exit(f"reuse_outputs_from: {_dst} already exists. Remove it first.")
        os.symlink(os.path.realpath(_src), _dst)
        print(f"reuse_outputs_from: {_subdir}/ -> {os.path.realpath(_src)}")
    del _reuse_from, _REUSE_DIRS

samples = parse_samples(config["samples"])
samples = find_raw_inputs(samples)


# Define target files for rule all
def get_demux_summaries(wildcards=None):
    """Per-run demux summaries, or nothing when no backend is enabled.

    Defined here rather than in demux.smk because `read_attrition` lives in the
    QC rules and must work either way: demux.smk is only included when a backend
    is on, so referencing its helpers unconditionally breaks every non-demux run.
    """
    if not (
        config.get("warpdemux", {}).get("enabled", False)
        or config.get("ldx", {}).get("enabled", False)
    ):
        return []
    run_ids = {
        info["run_id"]
        for info in samples.values()
        if info.get("barcode") and info.get("run_id")
    }
    return [
        os.path.join(outdir, "demux", "read_ids", rid, "demux_summary.tsv.gz")
        for rid in sorted(run_ids)
    ]


def get_assigned_summaries(wildcards=None):
    """Per-run `assigned_summary.tsv` from the escpod join, or nothing.

    Written by ldx_run_read_ids, so only the LDX path has one. read_attrition
    uses it to split the barcode-assigned -> basecalled gate into "no sample
    claims this read" (unclaimed code, or a dual-index pair the axes did not
    agree on) and "the basecaller could not read it", which that row otherwise
    conflates -- on the LDX fixture the 40 decoy reads looked like basecaller
    loss.
    """
    if not config.get("ldx", {}).get("enabled", False):
        return []
    run_ids = {
        info["run_id"]
        for info in samples.values()
        if info.get("barcode") and info.get("run_id")
    }
    return [
        os.path.join(outdir, "demux", "read_ids", rid, "assigned_summary.tsv")
        for rid in sorted(run_ids)
    ]


def get_fdx_summaries():
    """Per-run tallies of the FDX axis, or nothing when it is off.

    Requested as a pipeline output because nothing downstream consumes them:
    the join reads the classifications, not the tally. In the two-pass shape
    the tally is written beside the fdx classifications anyway; in the fused
    shape it is what makes summarize_fdx_axis run at all.
    """
    if not is_fdx_enabled():
        return []
    run_ids = {
        info["run_id"]
        for info in samples.values()
        if info.get("barcode") and info.get("run_id")
    }
    return [
        os.path.join(outdir, "demux", "read_ids", rid, "fdx", "demux_summary.tsv.gz")
        for rid in sorted(run_ids)
    ]


def pipeline_outputs():
    outs = expand(
        os.path.join(
            outdir, "summary", "tables", "{sample}", "{sample}.charging_prob.tsv.gz"
        ),
        sample=samples.keys(),
    )

    outs += expand(
        os.path.join(
            outdir, "summary", "tables", "{sample}", "{sample}.charging.cpm.tsv.gz"
        ),
        sample=samples.keys(),
    )

    outs += expand(
        os.path.join(
            outdir, "summary", "tables", "{sample}", "{sample}.bcerror.tsv.gz"
        ),
        sample=samples.keys(),
    )

    outs += expand(
        os.path.join(
            outdir, "summary", "tables", "{sample}", "{sample}.align_stats.tsv.gz"
        ),
        sample=samples.keys(),
    )

    outs += expand(
        os.path.join(
            outdir, "summary", "tables", "{sample}", "{sample}.anchor_coverage.tsv.gz"
        ),
        sample=samples.keys(),
    )

    # One table per run saying where the reads went. Always on: every gate's
    # loss is already implicit in some artifact, but only as a difference
    # between rows in different files, which is why a 12% drop at charge-calling
    # went unnoticed for every run before 2026-08-09 (issue #110).
    outs.append(os.path.join(outdir, "summary", "read_attrition.tsv.gz"))

    outs += expand(
        os.path.join(
            outdir, "summary", "tables", "{sample}", "{sample}.{values}.bg.gz"
        ),
        sample=samples.keys(),
        values=["cpm", "counts"],
    )

    # modkit outputs
    outs += expand(
        os.path.join(
            outdir,
            "summary",
            "modkit",
            "{sample}",
            "{sample}.pileup.bed.gz",
        ),
        sample=samples.keys(),
    )

    # outs += expand(
    #     os.path.join(
    #         outdir, "summary", "modkit", "{sample}", "{sample}.mod_calls.tsv.gz"
    #     ),
    #     sample=samples.keys(),
    # )

    # outs += expand(
    #     os.path.join(
    #         outdir, "summary", "modkit", "{sample}", "{sample}.mod_full.tsv.gz"
    #     ),
    #     sample=samples.keys(),
    # )

    # Per-read base-calling error calls, for co-occurrence analysis that does
    # not depend on the modification caller. Off by default: the BAM walk is
    # expensive and only some projects need read-level calls.
    if config.get("mismatch_calls", {}).get("enabled", False):
        outs.append(os.path.join(outdir, "summary", "tables", "bcerror_sites.tsv.gz"))

        outs += expand(
            os.path.join(
                outdir,
                "summary",
                "tables",
                "{sample}",
                "{sample}.mismatch_calls.tsv.gz",
            ),
            sample=samples.keys(),
        )

        outs += expand(
            os.path.join(
                outdir,
                "summary",
                "tables",
                "{sample}",
                "{sample}.charging_error.tsv.gz",
            ),
            sample=samples.keys(),
        )

    # outs += expand(
    #     os.path.join(
    #         outdir, "summary", "tables", "{sample}", "{sample}.odds_ratios.tsv.gz"
    #     ),
    #     sample=samples.keys(),
    # )

    # outs += expand(
    #     os.path.join(
    #         outdir,
    #         "summary",
    #         "tables",
    #         "{sample}",
    #         "{sample}.odds_ratios_filtered.tsv.gz",
    #     ),
    #     sample=samples.keys(),
    # )

    # Per-read charging calls, with a `reason` row for every read the model
    # did NOT score. Not optional: abstention is charging-correlated, so a
    # charging fraction without its no-call rate beside it is biased low.
    outs += expand(
        os.path.join(
            outdir, "summary", "tables", "{sample}", "{sample}.charging_calls.tsv.gz"
        ),
        sample=samples.keys(),
    )

    # tRNA-only reference FASTA (adapters stripped)
    outs.append(get_trna_fasta())

    # Squiggy session file for loading samples in Positron
    outs.append(os.path.join(outdir, "squiggy-session.json"))

    # Reference sequence similarity QC (runs once per pipeline execution)
    if want_reference_similarity():
        outs.append(os.path.join(outdir, "summary", "qc", "reference_similarity.tsv"))

    # EDX (3' adapter barcode) concordance table
    if config.get("edx", {}).get("enabled", False) and get_edx_samples():
        outs.append(os.path.join(outdir, "summary", "edx", "edx_concordance.tsv.gz"))

    # Per-run tally of the FDX (5' index) axis, beside the LDX one
    outs += get_fdx_summaries()

    return outs


wildcard_constraints:
    sample="|".join(samples.keys()),


# various additional helper functions
def get_raw_inputs(wildcards):
    return samples[wildcards.sample]["raw_files"]


def get_basecalling_dir(wildcards):
    return samples[wildcards.sample]["path"]


def get_modified_bases():
    """Parse modification names from the dorado opts string.

    Returns list of modification names, e.g. ["pseU", "m5C", "inosine_m6A"].
    """
    dorado_opts = config.get("opts", {}).get("dorado", "")
    if "--modified-bases" not in dorado_opts:
        return []
    # extract tokens after --modified-bases until the next flag or end of string
    parts = dorado_opts.split("--modified-bases")[1].split()
    mods = []
    for part in parts:
        if part.startswith("--"):
            break
        mods.append(part)
    return mods


def count_fasta_seqs(path):
    """Count records in a FASTA by '>' lines. Returns None if unreadable."""
    try:
        with open(path) as f:
            return sum(1 for line in f if line.startswith(">"))
    except OSError:
        return None


def get_similarity_max_mismatch():
    """Hamming distance for collapsing near-identical reference sequences.

    0 collapses exact duplicates only, which is lossless.
    """
    return config.get("qc", {}).get("reference_similarity_max_mismatch", 0)


def want_reference_similarity():
    """Whether to emit the reference similarity QC matrix.

    Alignment count is quadratic in the number of distinct reference sequences,
    and the downstream heatmap stops being legible well before it gets slow, so
    a large reference is skipped with a warning rather than silently costing
    hours and gigabytes.
    """
    qc = config.get("qc", {})
    if not qc.get("reference_similarity", True):
        return False

    max_seqs = qc.get("reference_similarity_max_seqs", 2000)
    if not max_seqs:
        return True

    # Collapsing is an explicit opt-in to a large reference: the matrix is then
    # reported over cluster representatives, so it does not grow with input size
    if get_similarity_max_mismatch() > 0:
        return True

    n_seqs = count_fasta_seqs(get_raw_reference())
    if n_seqs is None:
        # Reference not readable yet; include the target rather than drop QC
        return True

    if n_seqs > max_seqs:
        logger.warning(
            f"Skipping reference similarity QC: {get_raw_reference()} has "
            f"{n_seqs} sequences (limit {max_seqs}). Raise or null out "
            f"qc.reference_similarity_max_seqs to run it anyway, and consider "
            f"setting qc.reference_similarity_max_mismatch to collapse "
            f"near-identical sequences."
        )
        return False

    return True


def get_modkit_threshold_opts():
    """
    Build modkit threshold options from config.
    Uses optimized thresholds based on ModkitOpt (Sneddon et al. 2025).
    See https://github.com/comprna/modkitopt for threshold optimization.
    """
    opts = []
    modkit_config = config.get("modkit", {})
    if modkit_config:
        # Global filter threshold for canonical base confidence
        filter_thresh = modkit_config.get("filter_threshold")
        if filter_thresh is not None:
            opts.append(f"--filter-threshold {filter_thresh}")
        # Per-modification pass thresholds (mod code or ChEBI ID)
        mod_thresholds = modkit_config.get("mod_thresholds", {})
        if mod_thresholds:
            for mod_code, threshold in mod_thresholds.items():
                if threshold is not None:
                    opts.append(f"--mod-thresholds {mod_code}:{threshold}")
    return " ".join(opts)


# WarpDemuX helper functions (additional functions in demux.smk)


def sample_needs_demux(sample):
    """Check if a sample needs demultiplexing (has barcode assigned)."""
    return is_demux_enabled() and samples[sample].get("barcode") is not None


def sample_has_edx(sample):
    """Check if a sample has an EDX adapter assignment."""
    return samples[sample].get("edx") is not None


def get_sample_edx(wildcards):
    """Return the EDX adapter name for a sample, or None if not set."""
    return samples[wildcards.sample].get("edx")


def sample_is_ldx(sample):
    """Check if a sample is demultiplexed by the escapepod (LDX) backend."""
    return is_ldx_enabled() and samples[sample].get("barcode") is not None


def get_ldx_pod5_source(run_id):
    """Return the one POD5 path to hand a classifier for an LDX run.

    There is no per-sample POD5 on the LDX path — the whole run stays in its raw
    POD5 and each sample is defined by its BAM. `escpod classify` takes one path,
    a POD5 *or* a directory, and it walks a directory RECURSIVELY: a run split
    across pod5_pass/pod5_fail is covered by naming the run directory itself.
    That is why there is no per-run `escpod merge` here any more — a merged copy
    of the run existed only to collapse two paths into one, and it duplicated the
    entire run's signal to do it, which is the cost the LDX path exists to avoid.

    Naming the run is a SUPERSET of its POD5 directories, and that is harmless
    here: classify is driven by the BAM, looking each aligned read up by id, so
    signal it is never asked for is never touched. The narrower path is still
    preferred when a run keeps its reads in one directory.
    """
    dirs = get_run_pod5_dirs(run_id)
    if len(dirs) == 1:
        return dirs[0]
    return get_run_path(run_id)


def get_sample_pod5(wildcards):
    """
    The sample's signal store, as handed to dorado and to `escpod classify`.

    LDX samples have no POD5 of their own: they resolve to the raw run.
    A WarpDemuX sample resolves to the split POD5 `escpod filter` wrote for it.
    Every other sample resolves to the directory of symlinks stage_pod5 laid
    over its raw files -- both tools take a directory, so nothing is copied.

    This is the whole store, never a subset: classification is driven by the
    BAM and looks each aligned read up by id, so an EDX sample's classifier
    reads exactly its own reads out of the shared store. (There used to be an
    EDX-filtered POD5 per sample for this, a copy kept forever; it went with
    merge_pods.)
    """
    if sample_is_ldx(wildcards.sample):
        source = get_ldx_pod5_source(samples[wildcards.sample]["run_id"])
        # A directory would make this input's mtime move whenever demux writes
        # a sidecar beside the POD5s, re-triggering classification for nothing.
        # Track the POD5 files themselves.
        return (
            samples[wildcards.sample]["raw_files"] if os.path.isdir(source) else source
        )
    if sample_needs_demux(wildcards.sample):
        return os.path.join(
            outdir, "demux", "pod5", wildcards.sample, f"{wildcards.sample}.pod5"
        )
    else:
        return os.path.join(outdir, "pod5", wildcards.sample)


def stage_pod5_links(raw_files, dest_dir, log_path):
    """Lay a directory of symlinks over a sample's raw POD5 files (stage_pod5).

    Each link is <dest>/<run>/<pod5_pass|pod5_fail|pod5>/<file>, mirroring where
    the file was found, so two runs pooled into one sample cannot collide on a
    basename. Targets are canonical paths. A link that already points at the
    right file is left alone; one pointing elsewhere is a real collision and is
    refused rather than silently repointed.
    """
    from snakemake.exceptions import WorkflowError

    n_linked = n_kept = 0
    for raw in raw_files:
        target = os.path.realpath(raw)
        subdir = os.path.dirname(raw)
        run = os.path.basename(os.path.normpath(os.path.dirname(subdir)))
        link_dir = os.path.join(dest_dir, run, os.path.basename(subdir))
        link = os.path.join(link_dir, os.path.basename(raw))
        os.makedirs(link_dir, exist_ok=True)
        if os.path.islink(link):
            if os.path.realpath(link) == target:
                n_kept += 1
                continue
            raise WorkflowError(
                f"stage_pod5: {link} already points at {os.path.realpath(link)}, "
                f"not {target}. Two of this sample's runs share a directory name "
                "and a file name; give them distinct run directory names."
            )
        if os.path.exists(link):
            raise WorkflowError(
                f"stage_pod5: {link} exists and is not a symlink; refusing to "
                "replace it (this directory is expected to hold only links)."
            )
        os.symlink(target, link)
        n_linked += 1
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    with open(log_path, "w") as fh:
        fh.write(f"{n_linked} links created, {n_kept} already in place\n")


def get_alignment_read_ids(wildcards):
    """The read-id list that bounds a sample's alignment, or nothing.

    An EDX sample aligns only the reads carrying its own 3' adapter, which
    detect_edx_adapters / extract_edx_read_ids identify on the uBAM before
    alignment. Every other sample aligns its whole uBAM. Returned as a list so
    bwa_align's input is empty rather than absent in the second case.
    """
    if sample_has_edx(wildcards.sample):
        return [
            os.path.join(
                outdir,
                "demux",
                "edx",
                wildcards.sample,
                f"{wildcards.sample}.edx_read_ids.txt",
            )
        ]
    return []


def get_read_group_args(wildcards):
    """Arguments to stamp_read_groups.py: the identity written into the BAM.

    bwa_align is where the sample's identity is written INTO the BAM, because it
    is the first place both demux backends have converged (see the note at the
    top of demux.smk). The @RG gets SM (sample), LB (run id, when the sample has
    one) and BC (barcode). Until this point the barcode lives only in the output
    path, and the per-read record that could recover it (demux/read_ids/) is
    deleted by the `clean` rule and, on the WarpDemuX path, is temp() under the
    demux_scratch tier.

    The @CO records the upstream (escapepod-models) barcode name next to ours
    whenever the two differ. Since the ldx16 switch neither live panel renames --
    a sample is configured as the name its bundle emits -- so this emits nothing
    today. Kept because the emitted vocabulary belongs to the bundle: the
    retired nbc16 panel did differ, and a future one may, and then a BAM tagged
    `ldx04` should still say what it came from.

    Unbarcoded samples get neither BC nor a barcode on the @RG: absence means
    "no demultiplexing", not "unknown barcode".
    """
    args = [f"--sample {wildcards.sample}"]
    library = samples[wildcards.sample].get("run_id")
    if library:
        args.append(f"--library {library}")
    label = get_sample_barcode_label(wildcards.sample)
    if label:
        args.append(f"--barcode {label}")
        upstream = get_sample_barcode_upstream(wildcards.sample)
        if upstream != label:
            args.append(f'--comment "aa-tRNA-seq:upstream_barcode={upstream}"')
    return " ".join(args)


def get_classification_pod5_arg(wildcards, input):
    """Return the POD5 argument to pass a classifier on the command line.

    Identical to `input.pod5` everywhere except the LDX path, where the input is
    the raw run's POD5 files (for dependency tracking) but the tool must be
    given the single directory holding them.
    """
    if sample_is_ldx(wildcards.sample):
        return get_ldx_pod5_source(samples[wildcards.sample]["run_id"])
    return input.pod5


def get_all_final_bams():
    """Return list of all final BAM files for all samples."""
    return expand(
        os.path.join(outdir, "bam", "final", "{sample}", "{sample}.bam"),
        sample=samples.keys(),
    )


def get_sample_pod5_files(sample):
    """The POD5 FILES holding a sample's signal, for the Squiggy session.

    Squiggy opens files, not directories, so this is the file-level view of
    what get_sample_pod5 hands the tools: the split POD5 for a WarpDemuX
    sample, and the raw run's own files for everyone else -- a staged
    directory is only symlinks to those, and an LDX sample has nothing else.
    """
    if sample_needs_demux(sample) and not sample_is_ldx(sample):
        return [os.path.join(outdir, "demux", "pod5", sample, f"{sample}.pod5")]
    return list(samples[sample]["raw_files"])


def get_all_pod5_files():
    """Every POD5 file any sample's session entry names, each once.

    Samples of one LDX run all point at the same raw files, and two unbarcoded
    samples may name the same run, so this is the one place a path can repeat.
    """
    seen = []
    for sample in samples:
        seen.extend(f for f in get_sample_pod5_files(sample) if f not in seen)
    return seen


rule generate_squiggy_session:
    """
    Generate squiggy session JSON file for loading samples in Positron.

    Creates a session file with absolute paths to POD5, BAM, and FASTA files
    along with MD5 checksums for integrity verification.
    """
    input:
        bams=get_all_final_bams(),
        pod5s=get_all_pod5_files(),
        fasta=config["fasta"],
    output:
        session=os.path.join(outdir, "squiggy-session.json"),
    log:
        os.path.join(outdir, "logs", "generate_squiggy_session.log"),
    params:
        src=SCRIPT_DIR,
        samples=" ".join(samples.keys()),
        outdir=outdir,
        pod5_args=" ".join(
            f"--pod5 {sample}={path}"
            for sample in samples
            for path in get_sample_pod5_files(sample)
        ),
    shell:
        """
        python {params.src}/generate_squiggy_session.py \
            --samples {params.samples} \
            {params.pod5_args} \
            --output-dir {params.outdir} \
            --fasta {input.fasta} \
            --output {output.session} \
            --no-checksums \
            2>&1 | tee {log}
        """


def check_basecaller_compatibility():
    """Fail (or warn) when this run's basecaller is not the bundle's.

    Runs while the DAG is built, so a mismatch costs a dry-run rather than a
    basecall plus a classification pass. Reads `metadata.json` only — the
    ~300 MB model digest that would prove byte identity is `pixi run
    verify-basecaller`, not this.

    `charging.basecaller_check` selects the strength: `error` (default),
    `warn`, or `off`. Only the MODEL-IDENTITY finding is governed by it; the
    dorado major-version finding always warns, because the model is what governs
    charging calls and blocking a run over the runtime that produced identical
    weights would strand every bundle built by an older dorado.
    """
    mode = config.get("charging", {}).get("basecaller_check", "error")
    if mode == "off":
        return
    if mode not in ("error", "warn"):
        sys.exit(
            f"charging.basecaller_check must be one of error, warn, off; "
            f"got {mode!r}"
        )

    findings = check_basecaller(
        get_charging_model(),
        config.get("base_calling_model", ""),
        config.get("dorado_version", ""),
    )
    fatal = [msg for level, msg in findings if level == "error" and mode == "error"]
    for level, msg in findings:
        label = "ERROR" if (level == "error" and mode == "error") else "WARNING"
        print(f"{label}: {msg}\n", file=sys.stderr)
    if fatal:
        sys.exit(
            "Refusing to run with a basecaller the charging bundle was not "
            "trained on. Set `charging.basecaller_check: warn` to proceed anyway."
        )


check_basecaller_compatibility()
