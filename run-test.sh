#! /usr/bin/env bash

#BSUB -n 1
#BSUB -J aatrnaseq-main-test
#BSUB -e .test/logs/aatrnaseq-main-test_%J.err
#BSUB -o .test/logs/aatrnaseq-main-test_%J.out

mkdir -p .test/logs

snakemake \
  --configfile=config/config-test.yml \
  --profile cluster
