#!/bin/bash

# download and unpack test data for the aa-trna-seq-pipeline Snakemake workflow

set -x -e

wget -q https://aatrnaseq-testdata.s3.amazonaws.com/test_data.tar.gz
tar -zxf test_data.tar.gz
rm test_data.tar.gz
