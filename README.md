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

Then you can run the pipeline by invoking snakemake in the pipeline root directory
```
cd ../
snakemake -c 1 -p
```

To use on `bodhi`, see the `run.sh` script and the `cluster/config.yaml` for configuration details.

# Requirements

The following tools are required and are expected to be present in your PATH. The instructions below install the necssary tools for a x86 linux platform. 

## Dorado

The dorado basecaller can be installed using pre-built binaries available from [github](https://github.com/nanoporetech/dorado?tab=readme-ov-file#installation)

```
wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.5.3-linux-x64.tar.gz
tar -zxvf dorado-0.5.3-linux-x64.tar.gz 
```

The pipeline will use the `dorado` executable in the `bin` directory.

Also available via modules on bodhi:
`module load dorado`

## Parasail

The parasail alignment library is compiled from source with instructions provided on [github](https://github.com/jeffdaily/parasail?tab=readme-ov-file#compiling-and-installing)

```
git clone https://github.com/jeffdaily/parasail/tree/master
cd parasail
autoreconf -fi
./configure --prefix=/a/path/to/install/location
# use multiple cores for make
make -j 4
make install
```
The pipeline will use the `parasail_aligner` executable in the `bin` directory.

Also available here on bodhi:
`/beevol/home/riemondy/Projects/AAtRNAseq/ext/parasail/bin`

## Other tools

Other CLI tools and python dependencies can be installed using conda/mamba, or manually installed. 

```
mamba env create -f environment.yml
```

