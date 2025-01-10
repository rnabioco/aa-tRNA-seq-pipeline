# aa-tRNA-seq-pipeline

A pipeline to process ONT aa-tRNA-seq data built using snakemake.
Downstream analysis to generate figures for the initial preprint can be found at: (https://github.com/rnabioco/aa-tRNA-seq)

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

## Workflow

Given a directory of pod5 files, this pipeline merges all files from each sample into a single pod5, rebasecalls them to generate an unmapped bam with move table information (for downstream use by Remora), converts the bam into a fastq, and aligns that fastq to a reference containing tRNA + adapter sequences with BWA MEM. The resulting data (pod5s and aligned reads) are then fed to a model trained using Remora to classify charged vs. uncharged reads in the rule `cca_classify`, generating numeric values indicating the likelihood of a read being aminoacylated in the `ML` tag of the BAM file. For classifying charged vs. uncharged reads, we treat ML values of 200-255 as aminoacylated in downstream steps, and values <200 as uncharged. This can be altered by adjusting the `ml-threshold` parameter in the rule `get_cca_trna_cpm`.

The final steps of the pipeline calculate a number of outputs that may be useful for analysis and visualization, including normalized counts for charged and uncharged tRNA (`get_cca_trna_cpm`), basecalling error values (`bcerror`), alignment statistics (`align_stats`) and information on raw nanopore signal from Remora (`remora_signal_stats`).

A few notes about Remora classification for charged vs. uncharged tRNA reads: First, this step retains only full length tRNA reads (with an allowance for signal loss at the 5´ end of nanopore direct RNA sequencing). Additionally, due to the iterative nature of sequencing method development, the present approach does not rely on differences in adapter sequences attached to charged vs. uncharged tRNA molecules (though these sequences are retained as separate entries in the alignment reference and downstream files). While we anticipate being able to leverage this information in the future, the current pipeline relies exclusively on signal data over a 6-nt modification kmer spanning the universal CCA 3′ end of tRNA and the first three nucleotides of the 3′ adapter (CCAGGC) to distinguish charged and uncharged reads.

The Directed Acyclic Graph (DAG) below provides an overview of the pipeline structure for a single sample.

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
