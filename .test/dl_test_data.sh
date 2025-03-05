#!/bin/bash
set -e

curl -L -o test_data.tar.gz https://aatrnaseq-testdata.s3.amazonaws.com/test_data.tar.gz

tar -zxvf test_data.tar.gz
rm test_data.tar.gz
