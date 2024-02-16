import argparse
import pysam
import numpy as np
from pathlib import Path

'''
This script takes the output merged bams from merge_bams.py, which contain mv tags.
It generates a new metric, samples per base, and adds it to the bam file as a new tag, 'sb'.

python add_sb_tags.py  \
    --input /path/to/input.bam \
    --output /path/to/output.bam
'''

def calculate_samples_per_base(mv_array):
    samples_per_base = []
    sample_count = 0
    stride_length = mv_array[0]
    for value in mv_array[1:]:
        if value == 1:
            if sample_count > 0:
                samples_per_base.append(sample_count)
            sample_count = stride_length
        else:
            sample_count += stride_length
    if sample_count > 0:
        samples_per_base.append(sample_count)
    return samples_per_base

def add_samples_per_base_tags(input_bam_path, output_bam_path):
    with pysam.AlignmentFile(input_bam_path, "rb") as infile, \
         pysam.AlignmentFile(output_bam_path, "wb", template=infile) as outfile:
        for read in infile:
            if read.has_tag("mv"):
                mv_array = read.get_tag("mv")
                samples_per_base = calculate_samples_per_base(mv_array)
                read.set_tag("sb", samples_per_base)
            outfile.write(read)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Add samples per base tags to BAM file.')
    parser.add_argument('--input', help='Path to input (merged, sorted) BAM file relative to root')
    parser.add_argument('--output', help='Path to output BAM file for samples per base tags, relative to root')
    
    args = parser.parse_args()
    
    input_bam = Path(args.input)
    output_bam = Path(args.output)
    
    add_samples_per_base_tags(input_bam, output_bam)
    pysam.index(str(output_bam))
