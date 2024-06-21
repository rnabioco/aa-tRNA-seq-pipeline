import os
import glob
import sys
import pysam

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
                print("samples file must have 2 columns, sample_id and data_path, separated by whitespace", file = sys.stderr)
                sys.exit(f"found {line}")
            if sample in samples:
                samples[sample]["path"].add(path)
            else: 
                samples[sample] = {"path" : {path}}
    return samples

def is_bam_file(fl):
    
    try:
        fo = pysam.AlignmentFile(fl, 'rb', check_sq = False)
        fo.close()
        return True
    except:
        return False

def get_pipeline_commit():
    from git import Repo
    repo = Repo(PIPELINE_DIR)
    return repo.head.commit

def print_config():
    print("Config settings:")
    for k,v in config.items():
        if k == "opts":
            print(f"\t{k}:")
            [print(f"\t\t{cmd}: {opts}")  for cmd,opts in v.items()]
        else:
            print(f"\t{k}: {v}")  

def find_raw_inputs(sample_dict):
    """
    parse through directories listed in samples.tsv and identify fast5 or pod5 files to process
    store input files and uuid base file names in dictionary for each sample
    """
    POD5_DIRS = ["pod5_pass", "pod5_fail"]
    FAST5_DIRS = ["fast5_pass", "fast5_fail"]
    fmt = config["input_format"]

    data_subdirs = [] 
    if fmt == "POD5":
        data_subdirs = POD5_DIRS
        ext = ".pod5"
    elif fmt == "FAST5":
        data_subdirs = FAST5_DIRS
        ext = ".fast5"
    elif fmt == "BAM":
        data_subdirs = []
        ext = ".bam"
    else:
        sys.exit("input_format config option must be either FAST5, POD5, or BAM")
    
    for sample, info in sample_dict.items():
        raw_fls = []
        if fmt == "BAM":
            bam_fl = info["path"]
            if len(bam_fl) > 1:
                sys.exit("BAM input format only accepts one file per sample")
            
            bam_fl = next(iter(bam_fl))

            if is_bam_file(bam_fl):
                raw_fls = [bam_fl]
            else:
                sys.exit(f"{bam_fl} is not a valid BAM file, check format")

        else:
            for path in info["path"]:
                for subdir in data_subdirs:
                    data_path = os.path.join(path, subdir, "*" + ext)
                    fls = glob.glob(data_path)
                    raw_fls += fls

        sample_dict[sample]["raw_files"] = raw_fls

    return sample_dict


# set up global samples dictionary to be used throughout pipeline
outdir = config["output_directory"]
rbc_outdir = config["rebasecalled_bam_directory"]

samples = parse_samples(config["samples"])
samples = find_raw_inputs(samples)

# Define target files for rule all
def pipeline_outputs():  
    outs = expand(os.path.join(outdir, "tables", "{sample}.bwa.bcerror.tsv"),
        sample = samples.keys())

    outs += expand(os.path.join(outdir, "tables", "{sample}.bwa.{values}.bg"),
        sample = samples.keys(),
        values = ["cpm","counts"])

    outs += [os.path.join(outdir, "tables", "align_stats.tsv")]
    return outs

wildcard_constraints:
    sample="|".join(samples.keys()),

# various additional helper functions
def get_raw_inputs(wildcards):
    return samples[wildcards.sample]["raw_files"]

def get_basecalling_dir(wildcards):
    return samples[wildcards.sample]["path"]
