from pathlib import Path
import pysam
import numpy as np
import pandas as pd
import argparse
from multiprocessing import Pool

'''
Take the sb tag values from a directory of bam files,
filter out only the mapped reads, and processes the
the values of this tag to retain sb values at mapped
nucleotides (match/mismatch positions), and append sb 
values from insertions to the last match/mismatch position.

Uses multiprocessing to parallelize.

Outputs a tsv file of read IDs and the retained sb values.
'''

def process_cigar(cigar_tuples, samples_per_base):
    sb_filtered = []
    sb_index = 0
    read_position = 0  # Initialize the read position

    for op, length in cigar_tuples:
        if op == 0:  # Match or mismatch
            sb_filtered.extend(samples_per_base[sb_index:sb_index + length])
            sb_index += length
            read_position += length  # Update read position
        elif op == 1:  # Insertion
            # Append insertion values to the last match or mismatch position
            sb_filtered.extend([samples_per_base[read_position - 1]] * length)
            sb_index += length
        elif op == 2:  # Deletion
            continue  # Do nothing
    return sb_filtered

# Define a function to process each BAM file
def process_bam_file(bam_file):
    print(f"Processing {bam_file}...")
    results = []
    with pysam.AlignmentFile(bam_file, 'rb') as infile:
        for i, read in enumerate(infile):
            if read.is_unmapped:
                continue
            if read.has_tag("sb"):
                sb_values = np.array(read.get_tag("sb"))
                cigar_tuples = read.cigar
                sb_filtered = process_cigar(cigar_tuples, sb_values)

                # Add to the results list
                results.append({
                    'read_id': read.query_name,
                    'reference_start': read.reference_start,
                    'processed_sb': ','.join(map(str, sb_filtered)),
                    'file': bam_file.name
                })

    return results

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Extract SB values')
    parser.add_argument('--output', help = 'Path to output CSV file', required = True)
    parser.add_argument('--proc', help = 'number of parallel processes', type = int, default = 1)
    parser.add_argument('BAM', help = 'Path to input BAM file(s)', nargs = '+')
    
    args = parser.parse_args()

    # Define the path to your BAM files
    bam_file_paths = [Path(x) for x in args.BAM]

    # Define the number of parallel processes to use
    num_processes = args.proc 

    # Path for saving the output
    output_file = args.output

    # Create a Pool of worker processes
    with Pool(num_processes) as pool:
        # Use the pool to process BAM files in parallel
        all_results = pool.map(process_bam_file, bam_file_paths)

    # Combine results from all processes
    results_df = pd.concat([pd.DataFrame(result) for result in all_results], ignore_index=True)

    # Save DataFrame to TSV file
    results_df.to_csv(output_file, sep='\t', index=False)
