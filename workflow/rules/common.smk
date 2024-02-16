import os
import glob
import sys

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

            samples[sample] = {"path" : path}
    return samples

def find_raw_inputs(sample_dict):
    """
    parse through directories listed in samples.tsv and identify fast5 or pod5 files to process
    store input files and uuid base file names in dictionary for each sample
    """
    POD5_DIR = "pod5_pass"
    FAST5_DIR = "fast5_pass"
    fmt = config["input_format"]

    if fmt == "POD5":
        data_subdir = POD5_DIR
        ext = ".pod5"
    elif fmt == "FAST5":
        data_subdir = FAST5_DIR
        ext = ".fast5"
    else:
        sys.exit("input_format config option must be either FAST5 or POD5")
    
    for sample, info in sample_dict.items():
        data_path = os.path.join(info["path"], data_subdir, "*" + ext)
        sample_dict[sample]["raw_data_path"] = os.path.join(info["path"], data_subdir) 
        fls = glob.glob(data_path)
        sample_dict[sample]["raw_files"] = {}
        for fl in fls:
            fl_id = os.path.basename(fl.replace(ext, ""))
            sample_dict[sample]["raw_files"][fl_id] = fl

    return sample_dict

# set up global samples dictionary to be used throughout pipeline
outdir = config["output_directory"]
rbc_outdir = config["rebasecalled_bam_directory"]
samples = parse_samples(config["samples"])
samples = find_raw_inputs(samples)

# Define target files for rule all:
# As of now it is a linear pipeline so only need to provide final files here 
def pipeline_outputs():  
    outs = expand(os.path.join(outdir, "tables", "{sample}.{aligner}.bcerror.tsv"),
        sample = samples.keys(),
        aligner = config["aligner"])
    sb_outs = [os.path.join(outdir, "tables", "sb_values.tsv")]
    outs.append(sb_outs)
    return outs

# various additional helper functions
def get_basecalling_inputs(wildcards):
    return samples[wildcards.sample]["raw_files"].values()

def get_basecalling_dir(wildcards):
    return samples[wildcards.sample]["raw_data_path"]

def get_rebasecalled_outputs(wildcards):
    partitions = samples[wildcards.sample]["raw_files"].keys()
    outputs = [os.path.join(outdir, "bams", wildcards.sample, p + ".bam") for p in partitions]
    return outputs

def get_aligner_output(wildcards):
    return os.path.join(outdir, "bams", wildcards.sample, wildcards.sample + "." + config["aligner"] + ".bam")
