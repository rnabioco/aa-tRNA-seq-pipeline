#! /usr/bin/env bash

#BSUB -n 1
#BSUB -J aatrnaseq-main
#BSUB -e logs/aatrnaseq-main_%J.err
#BSUB -o logs/aatrnaseq-main_%J.out

mkdir -p logs

snakemake \
    --configfile=config/config-preprint.yml \
    --profile cluster \
    --rerun-incomplete
