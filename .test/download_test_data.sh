#!/bin/bash
set -x -e

wget --progress=bar:force:noscroll https://aatrnaseq-testdata.s3.amazonaws.com/test_data.tar.gz
tar -zxvf test_data.tar.gz
rm test_data.tar.gz
