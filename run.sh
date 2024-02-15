#! /usr/bin/env bash
#BSUB -n 1
#BSUB -J aatrnaseq 
#BSUB -e logs/snakemake_%J.err
#BSUB -o logs/snakemake_%J.out

mkdir -p logs

snakemake \
  --configfile=config/config.yml \
  --profile cluster

# run with parasail instead (will overwrite output bams and summary tables)
snakemake \
  --configfile=config/config.yml \
  --config aligner=parasail \
  --profile cluster
