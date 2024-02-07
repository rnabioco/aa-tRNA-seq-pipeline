import argparse
import pysam
import os
from pathlib import Path

'''
This script takes unmapped bams with mv tags and merges them to mapped bams without them.
It uses pathlib to specify a root directory, so that you can
then give it relative paths for the unmapped and mapped bam files, e.g.:

python merge_bams.py --root "/path/to/root/directory/" \
    --ubam "relpath/to/rbc.bam" \
    --mbam "relpath/to/bwa.bam" \
    --merged "relpath/to/merged.bam"

'''

def index_bam(bam_path):
    pysam.index(str(bam_path))

def read_mv_tags(ubam_path):
    mv_tags = {}
    with pysam.AlignmentFile(ubam_path, "rb", check_sq=False) as infile:
        for read in infile:
            if read.has_tag("mv"):
                mv_tags[read.query_name] = read.get_tag("mv")
    return mv_tags

def merge_bams(ubam_path, mbam_path, merged_bam_path, mv_tags):
    with pysam.AlignmentFile(mbam_path, "rb") as infile, \
         pysam.AlignmentFile(merged_bam_path, "wb", template=infile) as outfile:
        for read in infile:
            if read.query_name in mv_tags:
                read.set_tag("mv", mv_tags[read.query_name])
            outfile.write(read)

def sort_bam(unsorted_bam_path, sorted_bam_path):
    pysam.sort("-o", str(sorted_bam_path), str(unsorted_bam_path))

def delete_file(file_path):
    os.remove(file_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge BAM files.')
    parser.add_argument('--root', help='Root directory of the project')
    parser.add_argument('--ubam', help='Path to Unmapped BAM file relative to root')
    parser.add_argument('--mbam', help='Path to Mapped BAM file relative to root')
    parser.add_argument('--merged', help='Path to save Merged BAM file relative to root')

    args = parser.parse_args()

    root = Path(args.root)
    ubam = root / args.ubam
    mbam = root / args.mbam
    merged_bam = root / args.merged
    sorted_merged_bam = merged_bam.with_suffix(".sorted.bam")

    # Index the unmapped BAM
    index_bam(ubam)

    # Read and store "mv" tags from the Unmapped BAM
    mv_tags = read_mv_tags(ubam)

    # Merge BAM files while adding "mv" tags
    merge_bams(ubam, mbam, merged_bam, mv_tags)

    # Sort the merged BAM
    sort_bam(merged_bam, sorted_merged_bam)

    # Index the output above
    index_bam(sorted_merged_bam)

    # Delete unsorted merged BAM
    delete_file(merged_bam)