# AAtRNAseqPipe

A pipeline to process ONT tRNA-seq data built using snakemake.

# Usage

The pipeline can be configued by editing the `config/config.yml` file. The config file specifications will
run a small example dataset through the pipeline. To download these data files:
```
git clone https://github.com/rnabioco/AAtRNAseqPipe.git
cd .test
# download test data 
bash dl_data.sh
```

Then you can test the pipeline by invoking snakemake in the pipeline root directory
```
cd ../
snakemake -c 1 -p
```

# Configuration

To use on your own samples you will need to edit the config.yml and samples.tsv files in the config directory. 
See README.md in the config directory for additional details.

# Cluster execution

To use on `bodhi`, see the `run.sh` script and the `cluster/config.yaml` for configuration details.

# Software requirements

The following tools are required and are expected to be present in your PATH. 

## Dorado

The dorado basecaller can be installed using pre-built binaries available from [github](https://github.com/nanoporetech/dorado?tab=readme-ov-file#installation)

## Other tools

Other CLI tools and python dependencies can be installed using conda/mamba, or manually installed. 

```bash
mamba env create -f environment.yml
conda activate aatrnaseqpipe 
```

If using a macOS arm64 CPU (e.g. M1-3) many bioinformatics dependencies are not yet available via conda for arm64.

To install these you'll need to use the x86 versions. Perform the following to set up the x86 conda enviroment 
and install the dependencies.

```bash
CONDA_SUBDIR=osx-64 mamba env create -f environment.yml 
conda activate aatrnaseqpipe
conda config --env --set subdir osx-64
```
