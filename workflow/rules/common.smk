import os
import glob
import sys
import pysam
from git import Repo

SCRIPT_DIR = os.path.join(SNAKEFILE_DIR, "scripts")


def parse_samples(fl):
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
                samples[sample] = {"path": {path}}
    return samples


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
    parse through directories listed in samples.tsv and identify fast5 or pod5 files to process
    store input files and uuid base file names in dictionary for each sample
    """
    POD5_DIRS = ["pod5_pass", "pod5_fail", "pod5"]
    FAST5_DIRS = ["fast5_pass", "fast5_fail"]
    fmt = config["input_format"]

    data_subdirs = []
    if fmt == "POD5":
        data_subdirs = POD5_DIRS
        ext = ".pod5"
    elif fmt == "FAST5":
        data_subdirs = FAST5_DIRS
        ext = ".fast5"
    else:
        sys.exit("input_format config option must be either FAST5, or POD5")

    for sample, info in sample_dict.items():
        raw_fls = []
        for path in info["path"]:
            for subdir in data_subdirs:
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
rbc_outdir = os.path.join(outdir, "rbc_bams")

samples = parse_samples(config["samples"])
samples = find_raw_inputs(samples)


# Define target files for rule all
def pipeline_outputs():
    outs = expand(
        os.path.join(outdir, "tables", "{sample}", "{sample}.charging_prob.tsv.gz"),
        sample=samples.keys(),
    )

    outs += expand(
        os.path.join(outdir, "tables", "{sample}", "{sample}.charging.cpm.tsv.gz"),
        sample=samples.keys(),
    )

    outs += expand(
        os.path.join(outdir, "tables", "{sample}", "{sample}.bcerror.tsv.gz"),
        sample=samples.keys(),
    )

    outs += expand(
        os.path.join(outdir, "tables", "{sample}", "{sample}.align_stats.tsv.gz"),
        sample=samples.keys(),
    )

    outs += expand(
        os.path.join(outdir, "tables", "{sample}", "{sample}.{values}.bg.gz"),
        sample=samples.keys(),
        values=["cpm", "counts"],
    )

    outs += expand(
        os.path.join(outdir, "squigualizer", "{sample}", "{sample}.tsv"),
        sample=samples.keys(),
    )

    if (
        "remora_kmer_table" in config
        and config["remora_kmer_table"] != ""
        and config["remora_kmer_table"] is not None
    ):
        outs += expand(
            os.path.join(outdir, "tables", "{sample}", "{sample}.remora.tsv.gz"),
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
