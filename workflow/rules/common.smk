import os
import glob
import sys
import pysam
import yaml
from git import Repo

SCRIPT_DIR = os.path.join(SNAKEFILE_DIR, "scripts")


def maybe_temp(path):
    """Mark path as temp() when cleanup_intermediates is enabled."""
    if config.get("cleanup_intermediates", True):
        return temp(path)
    return path


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
                samples[sample] = {
                    "path": {path},
                    "barcode": None,
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
        barcode_kit = run.get(
            "barcode_kit", config.get("warpdemux", {}).get("barcode_kit")
        )

        for sample_name, sample_val in run["samples"].items():
            if sample_name in samples:
                sys.exit(f"Duplicate sample name '{sample_name}' in samples file: {fl}")

            # Backward compat: plain string or null = wdx barcode only
            if isinstance(sample_val, str) or sample_val is None:
                barcode = sample_val
                edx = None
            elif isinstance(sample_val, dict):
                barcode = sample_val.get("wdx")
                edx = sample_val.get("edx")
            else:
                sys.exit(f"Invalid sample value for '{sample_name}': {sample_val}")

            samples[sample_name] = {
                "path": {run_path},
                "barcode": barcode,
                "edx": edx,
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

    _REUSE_DIRS = ["pod5", "demux", "bam/rebasecall", "fq"]
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

    # modkit outputs — only when modified bases are called. Canonical basecalls
    # (e.g. dorado_opts_override without --modified-bases, used for cross-version
    # comparison) carry no MM/ML tags, and `modkit pileup` errors on them.
    if get_modified_bases():
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

    # tRNA-only reference FASTA (adapters stripped)
    outs.append(get_trna_fasta())

    # Squiggy session file for loading samples in Positron
    outs.append(os.path.join(outdir, "squiggy-session.json"))

    # Reference sequence similarity QC (runs once per pipeline execution)
    outs.append(os.path.join(outdir, "summary", "qc", "reference_similarity.tsv"))

    # Amino acid classification outputs (leech multiclass)
    if config.get("classify_aa", {}).get("enabled", False):
        outs += expand(
            os.path.join(
                outdir,
                "summary",
                "tables",
                "{sample}",
                "{sample}.aa_classify.tsv.gz",
            ),
            sample=samples.keys(),
        )

    # EDX (3' adapter barcode) concordance table
    if config.get("edx", {}).get("enabled", False) and get_edx_samples():
        outs.append(os.path.join(outdir, "summary", "edx", "edx_concordance.tsv.gz"))

    # AA identity classification (one-vs-all bundle)
    if config.get("aa_identity", {}).get("enabled", False):
        outs += expand(
            os.path.join(outdir, "bam", "aa_classified", "{sample}", "{sample}.bam"),
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


def get_modified_bases():
    """Parse modification names from the dorado opts string.

    Returns list of modification names, e.g. ["pseU", "m5C", "inosine_m6A"].

    Honors `dorado_opts_override` so a canonical run (no --modified-bases) does
    not try to fetch mod models the active dorado_model lacks.
    """
    dorado_opts = config.get(
        "dorado_opts_override", config.get("opts", {}).get("dorado", "")
    )
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


def get_sample_pod5(wildcards):
    """
    Return the correct POD5 path for a sample.
    If demux is enabled and sample has barcode, use split POD5.
    Otherwise, use merged POD5 from merge_pods rule.
    """
    if sample_needs_demux(wildcards.sample):
        return os.path.join(
            outdir, "demux", "pod5", wildcards.sample, f"{wildcards.sample}.pod5"
        )
    else:
        return os.path.join(
            outdir, "pod5", wildcards.sample, f"{wildcards.sample}.pod5"
        )


def get_alignment_fastq(wildcards):
    """
    Return the correct FASTQ path for alignment.
    EDX samples use the EDX-filtered FASTQ; others use ubam_to_fastq output.
    """
    if sample_has_edx(wildcards.sample):
        return os.path.join(
            outdir,
            "demux",
            "edx",
            "fq",
            wildcards.sample,
            f"{wildcards.sample}.fq.gz",
        )
    return os.path.join(outdir, "fq", wildcards.sample, f"{wildcards.sample}.fq.gz")


def get_classification_pod5(wildcards):
    """
    Return the correct POD5 path for classification/signal analysis.
    EDX samples use the EDX-filtered POD5; others use the WDX-split or merged POD5.

    NOTE: Do NOT use this for rebasecall — rebasecall needs the pre-EDX POD5
    (use get_sample_pod5 instead).
    """
    if sample_has_edx(wildcards.sample):
        return os.path.join(
            outdir,
            "demux",
            "edx",
            "pod5",
            wildcards.sample,
            f"{wildcards.sample}.pod5",
        )
    return get_sample_pod5(wildcards)


def get_all_final_bams():
    """Return list of all final BAM files for all samples."""
    return expand(
        os.path.join(outdir, "bam", "final", "{sample}", "{sample}.bam"),
        sample=samples.keys(),
    )


def get_all_merged_pod5s():
    """Return list of all merged/filtered POD5 files for all samples.

    For EDX samples, returns the EDX-filtered POD5 (subset matching the sample's adapter).
    For WDX-only samples, returns the WDX-split POD5.
    For non-demux samples, returns the merged POD5.
    """
    pod5_paths = []
    for sample in samples.keys():
        if sample_has_edx(sample):
            pod5_paths.append(
                os.path.join(outdir, "demux", "edx", "pod5", sample, f"{sample}.pod5")
            )
        elif sample_needs_demux(sample):
            pod5_paths.append(
                os.path.join(outdir, "demux", "pod5", sample, f"{sample}.pod5")
            )
        else:
            pod5_paths.append(os.path.join(outdir, "pod5", sample, f"{sample}.pod5"))
    return pod5_paths


rule generate_squiggy_session:
    """
Generate squiggy session JSON file for loading samples in Positron.

Creates a session file with absolute paths to POD5, BAM, and FASTA files
along with MD5 checksums for integrity verification.
"""
    input:
        bams=get_all_final_bams(),
        pod5s=get_all_merged_pod5s(),
        fasta=config["fasta"],
    output:
        session=os.path.join(outdir, "squiggy-session.json"),
    log:
        os.path.join(outdir, "logs", "generate_squiggy_session.log"),
    params:
        src=SCRIPT_DIR,
        samples=" ".join(samples.keys()),
        outdir=outdir,
    shell:
        """
        python {params.src}/generate_squiggy_session.py \
            --samples {params.samples} \
            --output-dir {params.outdir} \
            --fasta {input.fasta} \
            --output {output.session} \
            --no-checksums \
            2>&1 | tee {log}
        """
