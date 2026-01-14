import os
import glob
import sys
import pysam
import yaml
from git import Repo

SCRIPT_DIR = os.path.join(SNAKEFILE_DIR, "scripts")


def is_demux_enabled():
    """Check if WarpDemuX demultiplexing is enabled in config."""
    return config.get("warpdemux", {}).get("enabled", False)


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
                samples[sample] = {"path": {path}, "barcode": None, "run_id": None}
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

      - path: /path/to/non-demux/run
        samples:
          direct_sample: ~  # null barcode = no demux
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
        barcode_kit = run.get("barcode_kit", config.get("warpdemux", {}).get("barcode_kit"))

        for sample_name, barcode in run["samples"].items():
            if sample_name in samples:
                sys.exit(f"Duplicate sample name '{sample_name}' in samples file: {fl}")
            samples[sample_name] = {
                "path": {run_path},
                "barcode": barcode,
                "run_id": run_id,
                "barcode_kit": barcode_kit,
            }

    return samples


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

    repo = Repo(PIPELINE_DIR)
    return repo.head.commit


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
    logger.info(f"Pipeline commit: {cid}")
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

samples = parse_samples(config["samples"])
samples = find_raw_inputs(samples)


# Define target files for rule all
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

    outs += expand(
        os.path.join(
            outdir, "summary", "modkit", "{sample}", "{sample}.mod_calls.tsv.gz"
        ),
        sample=samples.keys(),
    )

    outs += expand(
        os.path.join(
            outdir, "summary", "modkit", "{sample}", "{sample}.mod_full.tsv.gz"
        ),
        sample=samples.keys(),
    )

    outs += [f"{DORADO_DIR}/bin/dorado"]

    outs += [os.path.join("resources/models", config["dorado_model"])]

    outs += [f"{MODKIT_DIR}/bin/modkit"]

    if (
        "remora_kmer_table" in config
        and config["remora_kmer_table"] != ""
        and config["remora_kmer_table"] is not None
    ):
        outs += expand(
            os.path.join(
                outdir, "summary", "tables", "{sample}", "{sample}.remora.tsv.gz"
            ),
            sample=samples.keys(),
        )

    return outs


wildcard_constraints:
    sample="|".join(samples.keys()),


# various additional helper functions
def get_raw_inputs(wildcards):
    return samples[wildcards.sample]["raw_files"]


def get_basecalling_dir(wildcards):
    return samples[wildcards.sample]["path"]


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


def sample_needs_demux(sample):
    """Check if a sample needs demultiplexing (has barcode assigned)."""
    return is_demux_enabled() and samples[sample].get("barcode") is not None


def get_sample_pod5(wildcards):
    """
    Return the correct POD5 path for a sample.
    If demux is enabled and sample has barcode, use split POD5.
    Otherwise, use merged POD5 from merge_pods rule.
    """
    if sample_needs_demux(wildcards.sample):
        return os.path.join(outdir, "demux", "pod5", f"{wildcards.sample}.pod5")
    else:
        return os.path.join(outdir, "pod5", wildcards.sample, f"{wildcards.sample}.pod5")
