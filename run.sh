#! /usr/bin/env bash
#BSUB -n 1
#BSUB -J aatrnaseq 
#BSUB -e logs/snakemake_%J.err
#BSUB -o logs/snakemake_%J.out

mkdir -p logs

snakemake \
  --configfile=config/config.yml \
  --profile cluster
  