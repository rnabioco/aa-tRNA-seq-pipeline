# aa-tRNA-seq-pipeline

A pipeline to process ONT aa-tRNA-seq data built using snakemake.

## Usage

The pipeline can be configued by editing the `config/config.yml` file. The config file specifications will
run a small example dataset through the pipeline. To download these data files:

```
git clone https://github.com/rnabioco/AAtRNAseqPipe.git
cd .test
# download test data 
bash dl_data.sh
```

Then you can test the pipeline by invoking snakemake in the pipeline root directory:

```
cd ../
snakemake -c 1 -p
```

## Configuration

To use on your own samples you will need to edit the config.yml and samples.tsv files in the config directory. 
See [README.md in the config directory](https://github.com/rnabioco/aa-tRNA-seq-pipeline/tree/main/config) for additional details.

## Workflow overview

The workflow for the aa-tRNA-seq pipeline is illustrated below. This Directed Acyclic Graph (DAG) provides an overview of the pipeline structure for a single sample.

![Workflow DAG](https://github.com/rnabioco/aa-tRNA-seq-pipeline/blob/main/workflow/workflow_dag.png)

To generate a DAG for your own configuration, use the following command in the pipeline root directory:

`snakemake --dag | dot -Tpng > workflow_dag.png`

You will need [Graphviz](https://graphviz.org) installed to run the above command.

## Cluster execution

The pipeline includes a `run.sh` script optimized for the authors' local compute cluster (`bodhi`) which uses the LSF scheduler. For more details on configuring for HPC jobs, see `cluster/config.yaml`.

## Software requirements

The following tools are required and are expected to be present in your PATH. 

### Dorado

The dorado basecaller can be installed using pre-built binaries available from [github](https://github.com/nanoporetech/dorado?tab=readme-ov-file#installation)

### Other tools

Other CLI tools and python dependencies can be installed using conda/mamba, or manually installed. Make sure to activate
the conda environment before running the pipeline.  

```bash
mamba env create -f environment.yml
conda activate aatrnaseqpipe 
```

## Notes

If using a macOS arm64 CPU (e.g. M1-3) many bioinformatics dependencies are not yet available via conda for arm64.

To install these you'll need to use the x86 versions. Perform the following to set up the x86 conda enviroment 
and install the dependencies.

```bash
CONDA_SUBDIR=osx-64 mamba env create -f environment.yml 
conda activate aatrnaseqpipe
conda config --env --set subdir osx-64
```
