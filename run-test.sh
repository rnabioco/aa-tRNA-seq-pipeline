#! /usr/bin/env bash

#BSUB -n 1
#BSUB -J aatrnaseq-main
#BSUB -e results/logs/aatrnaseq-main_%J.err
#BSUB -o results/logs/aatrnaseq-main_%J.out

mkdir -p results/logs

snakemake \
  --configfile=config/config.yml \
  --profile cluster
