#! /usr/bin/env bash

# generate test data for snakemake pipeline
# upload the test_data.tar.gz tarball  to aws s3 (s3://aatrnaseq-testdata)
# e.g. aws s3 cp test_data.tar.gz s3://aatrnaseq-testdata/

# get reference fasta
cp /beevol/home/riemondy/Projects/AAtRNAseq/ref/sacCer3-mature-tRNAs-dual-adapt-v2.fa ./

# get reference tRNA isodecoder mapping table
cp /beevol/home/riemondy/Projects/AAtRNAseq/ref/sacCer3-mature-tRNAs-dual-adapt-v2.table ./

# get dorado model
dorado download --model rna004_130bps_sup@v5.0.0

# extract a subset of reads and make small pod5 and fast5s for testing
rm -rf sample1 sample2 sample2_1
ex=/beevol/home/riemondy/Projects/AAtRNAseq/data/2024-07-25_Sami_Phizickystrains

samtools view -F 20 $ex/bams/JMW_510_28C/JMW_510_28C.bwa.bam "Ala-AGC" | grep -v "pi:Z:" | cut -f 1 | head -n 100  > ex1_read_ids_1.txt
samtools view -F 20 $ex/bams/JMW_510_28C/JMW_510_28C.bwa.bam "Ala-AGC-uncharged" | grep -v "pi:Z:" | cut -f 1 | tail -n 100  > ex1_read_ids_2.txt
samtools view -s 42.005 -F 20 $ex/bams/JMW_510_28C/JMW_510_28C.bwa.bam | grep -v "pi:Z:" | cut -f 1 | head -n 200 | tail -n 10   > ex1_read_ids_3.txt

samtools view -F 20 $ex/bams/JMW_510_37C/JMW_510_37C.bwa.bam | grep -v "pi:Z:" | cut -f 1 | head -n 100  > ex2_read_ids_1.txt
samtools view -F 20 $ex/bams/JMW_510_37C/JMW_510_37C.bwa.bam | grep -v "pi:Z:" |cut -f 1 | tail -n 100  > ex2_read_ids_2.txt
samtools view -s 42.005 -F 20 $ex/bams/JMW_510_37C/JMW_510_37C.bwa.bam | grep -v "pi:Z:" | cut -f 1 | head -n 200 | tail -n 10   > ex2_read_ids_3.txt

od=sample1/pod5_pass
rm -rf $od
mkdir -p $od
pod5 filter $ex/rbc/JMW_510_28C/JMW_510_28C.pod5 --ids ex1_read_ids_1.txt --force-overwrite -o $od/1.pod5
pod5 filter $ex/rbc/JMW_510_28C/JMW_510_28C.pod5 --ids ex1_read_ids_2.txt --force-overwrite -o $od/2.pod5

od=sample1/pod5_fail
rm -rf $od
mkdir -p $od
pod5 filter $ex/rbc/JMW_510_28C/JMW_510_28C.pod5 --ids ex1_read_ids_3.txt --force-overwrite -o $od/1.pod5

od=sample2/pod5_pass
mkdir -p $od
pod5 filter $ex/rbc/JMW_510_37C/JMW_510_37C.pod5 --ids ex2_read_ids_1.txt --force-overwrite -o $od/1.pod5
pod5 filter $ex/rbc/JMW_510_37C/JMW_510_37C.pod5 --ids ex2_read_ids_2.txt --force-overwrite -o $od/2.pod5

od=sample2/pod5_fail
mkdir -p $od
pod5 filter $ex/rbc/JMW_510_37C/JMW_510_37C.pod5 --ids ex2_read_ids_3.txt --force-overwrite -o $od/1.pod5

# make another "sample2" dataset to test merging multiple runs
# don't duplicate sample2 reads to avoid throwing an error when merging
# pod5s
cp -r sample1 sample2_1

rm ex1*.txt ex2*.txt

# get kmer levels for remora  
wget https://raw.githubusercontent.com/nanoporetech/kmer_models/master/rna004/9mer_levels_v1.txt

# tarball the test data
tar -czf test_data.tar.gz \
  9mer_levels_v1.txt \
  rna004_130bps_sup@v5.0.0 \
  sample* \
  sacCer3-mature-tRNAs-dual-adapt-v2.fa \
  sacCer3-mature-tRNAs-dual-adapt-v2.table

