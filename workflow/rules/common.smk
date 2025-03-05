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
            fields = line.split()

            if len(fields) == 2:
                # If samples.tsv has the old format, assume aa-tRNA-seq input
                try:
                    sample, path = fields
                except ValueError:
                    print(
                        "samples file must have 2 columns (sample_id and data_path, in which case "
                        "aa-tRNA-seq input will be assumed), or 5 columns (sample_id, data_path, "
                        "sequencing_input, organism, chemistry) separated by whitespace",
                        file=sys.stderr,
                    )
                    sys.exit(f"found {line}")

                sequencing_input = "aa-tRNA"
                organism = "scerevisiae"  # ✅ Fixed typo (was `organisms`)
                chemistry = "RNA004"
                basecall_model = "sup"

            elif len(fields) == 5:
                # New format, use provided values
                try:
                    sample, path, sequencing_input, organism, chemistry = fields
                except ValueError:
                    print(
                        "sample file must have either 2 or 5 columns, separated by whitespace.",
                        file=sys.stderr
                    )
                    sys.exit(f"found {line}")

            else:
                print(
                    "Error: samples file must have either 2 or 5 columns:\n"
                    "2-column format: sample_id, data_path (defaults to scerevisiae RNA004 aa-tRNA)\n"
                    "5-column format: sample_id, data_path, sequencing_input, organism, chemistry",
                    file=sys.stderr,
                )
                sys.exit(f"found {line}")

            if sample in samples:
                print(f"Duplicate sample found: {sample}", file=sys.stderr)
                sys.exit(1)
            else:
                samples[sample] = {
                    "path": path,
                    "sequencing_input": sequencing_input,
                    "organism": organism,
                    "chemistry": chemistry,
                    "basecall_model": basecall_model
                }

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
    Parse directories listed in samples.tsv and identify FAST5 or POD5 files to process.
    Store input files and UUID base file names in a dictionary for each sample.
    """
    POD5_DIRS = ["pod5_pass", "pod5_fail"]
    FAST5_DIRS = ["fast5_pass", "fast5_fail"]
    fmt = config["input_format"]  # Get input format from config file

    print(f"DEBUG: input_format set to {fmt}", file=sys.stderr)  # Print input format for debugging

    # Select correct subdirectories and file extension
    if fmt == "POD5":
        data_subdirs = POD5_DIRS
        ext = ".pod5"
    elif fmt == "FAST5":
        data_subdirs = FAST5_DIRS
        ext = ".fast5"
    else:
        sys.exit("ERROR: input_format in config must be either FAST5 or POD5")

    print(f"DEBUG: Searching for {ext} files in {data_subdirs}", file=sys.stderr)  # Debug file search path

    for sample, info in sample_dict.items():
        raw_fls = []
        print(f"DEBUG: Processing sample {sample}, Looking inside: {info['path']}", file=sys.stderr)

        # Ensure path is stored as a string (not split into characters!)
        if isinstance(info["path"], set):
            path_list = list(info["path"])  # Convert set to list if necessary
        else:
            path_list = [info["path"]]  # Ensure it's a list for iteration

        for path in path_list:
            absolute_path = os.path.abspath(path)  # Ensure correct absolute path
            for subdir in data_subdirs:
                data_path = os.path.join(absolute_path, subdir, "*" + ext)  # Corrected path handling
                fls = glob.glob(data_path)

                print(f"DEBUG: Sample {sample}, Searching {data_path}, Found: {fls}", file=sys.stderr)

                raw_fls += fls

        if len(raw_fls) == 0:
            sys.exit(
                f"ERROR: No input files found for sample: {sample}. Please check the path in the samples.tsv file"
            )

        sample_dict[sample]["raw_files"] = raw_fls

    return sample_dict

# set up global samples dictionary to be used throughout pipeline
outdir = config["output_directory"]
rbc_outdir = os.path.join(outdir, "rbc_bams")

samples = parse_samples(config["samples"])
print("Parsed samples:", samples, file=sys.stderr)
for sample, info in samples.items():
    print(f"Sample: {sample}, Paths: {info['path']}", file=sys.stderr)

samples = find_raw_inputs(samples)
print("Parsed samples dictionary:", samples, file=sys.stderr)

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

    if (
        "remora_kmer_table" in config
        and config["remora_kmer_table"] != ""
        and config["remora_kmer_table"] is not None
    ):
        outs += expand(
            os.path.join(outdir, "tables", "{sample}", "{sample}.remora.tsv.gz"),
            sample=samples.keys(),
        )

    # if "trna_table" in config and config["trna_table"] != "" and config["trna_table"] is not None:
    #    outs += expand(os.path.join(outdir, "tables", "{sample}", "{sample}.charging_status.tsv"),
    #        sample = samples.keys())

    return outs


wildcard_constraints:
    sample="|".join(samples.keys()),


# various additional helper functions
def get_raw_inputs(wildcards):
    return samples[wildcards.sample]["raw_files"]


def get_basecalling_dir(wildcards):
    return samples[wildcards.sample]["path"]
