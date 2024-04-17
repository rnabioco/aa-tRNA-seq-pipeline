
source $CONDA_PREFIX/etc/profile.d/conda.sh

if { mamba env list | grep 'aatrnaseqpipe'; } >/dev/null 2>&1; then
  conda activate aatrnaseqpipe
else
  mamba env create -f environment.yml
fi

cd .test
bash  dl_test_data.sh
cd -

snakemake \
  -c 1 \
  -p \
  --configfile=config/config.yml  \
  --config \
    samples=config/samples_bam.tsv \
    input_format=BAM

