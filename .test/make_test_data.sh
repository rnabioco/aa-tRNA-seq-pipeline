#! /usr/bin/env bash

# generate test data for snakemake pipeline
# upload the test_data.tar.gz tarball  to aws s3 (s3://aatrnaseq-testdata)
# e.g. aws s3 cp test_data.tar.gz s3://aatrnaseq-testdata/

# get reference fasta
cp /beevol/home/riemondy/Projects/AAtRNAseq/ref/yeast.tRNAnanoref.wsplicing.fa ./


# get dorado model
dorado download --model rna004_130bps_sup@v3.0.1

# extract a subset of reads and make small pod5 and fast5s for testing
rm -rf sample1 sample2
ex1=/beevol/home/whitel/data/AAtRNAseq/biological_tRNA/20240126_S288C_3AT/20240126_1602_P2S-01618-B_PAS98845_337afbe4
ex2=/beevol/home/whitel/data/AAtRNAseq/biological_tRNA/20231211_S288C_ctrl_chemligonly/20231211_1443_MN35252_FAX71838_5dd8ff5a

samtools view -q 30 -F 20 $ex1/20240126_S288C_3AT.bwa.bam | cut -f 1 | head -n 100  > ex1_read_ids_1.txt
samtools view -q 30 -F 20 $ex1/20240126_S288C_3AT.bwa.bam | cut -f 1 | tail -n 100  > ex1_read_ids_2.txt
samtools view -q 30 -F 20 $ex1/20240126_S288C_3AT.bwa.bam | cut -f 1 | head -n 200 | tail -n 10   > ex1_read_ids_3.txt

samtools view -q 30 -F 20 $ex2/20231211_S288C_ctrl_chemligonly.bwa.bam | cut -f 1 | head -n 100  > ex2_read_ids_1.txt
samtools view -q 30 -F 20 $ex2/20231211_S288C_ctrl_chemligonly.bwa.bam | cut -f 1 | tail -n 100  > ex2_read_ids_2.txt
samtools view -q 30 -F 20 $ex2/20231211_S288C_ctrl_chemligonly.bwa.bam | cut -f 1 | head -n 200 | tail -n 10   > ex2_read_ids_3.txt

rm -rf sample1/ sample2/ sample2_1/
od=sample1/pod5_pass
rm -rf $od
mkdir -p $od
pod5 filter $ex1/pod5_pass/*.pod5 --ids ex1_read_ids_1.txt --force-overwrite -o $od/1.pod5
pod5 filter $ex1/pod5_pass/*.pod5 --ids ex1_read_ids_2.txt --force-overwrite -o $od/2.pod5

od=sample1/pod5_fail
rm -rf $od
mkdir -p $od
pod5 filter $ex1/pod5_pass/*.pod5 --ids ex1_read_ids_3.txt --force-overwrite -o $od/1.pod5

od=sample2/pod5_pass
mkdir -p $od
pod5 filter $ex2/pod5_pass/*.pod5 --ids ex2_read_ids_1.txt --force-overwrite -o $od/1.pod5
pod5 filter $ex2/pod5_pass/*.pod5 --ids ex2_read_ids_2.txt --force-overwrite -o $od/2.pod5

od=sample2/pod5_fail
mkdir -p $od
pod5 filter $ex2/pod5_pass/*.pod5 --ids ex2_read_ids_3.txt --force-overwrite -o $od/1.pod5

od=sample1/fast5_pass
mkdir -p $od
pod5 convert to_fast5 -f -o $od sample1/pod5_pass/*.pod5

od=sample1/fast5_fail
mkdir -p $od
pod5 convert to_fast5 -f -o $od sample1/pod5_fail/*.pod5

od=sample2/fast5_pass
mkdir -p $od
pod5 convert to_fast5 -f -o $od sample2/pod5_pass/*.pod5

od=sample2/fast5_fail
mkdir -p $od
pod5 convert to_fast5 -f -o $od sample2/pod5_fail/*.pod5

dorado basecaller --emit-moves -v -r rna004_130bps_sup@v3.0.1 sample1/pod5_pass > sample1/sample1.unmapped.bam
dorado basecaller --emit-moves -v -r rna004_130bps_sup@v3.0.1 sample2/pod5_pass > sample2/sample2.unmapped.bam

# make another "sample2" dataset to test merging multiple runs
# don't duplicate sample2 reads to avoid throwing an error when merging
# pod5s
cp -r sample1 sample2_1

rm ex1*.txt ex2*.txt

tar -czf test_data.tar.gz rna004_130bps_sup@v3.0.1 sample* yeast.tRNAnanoref.wsplicing.fa
