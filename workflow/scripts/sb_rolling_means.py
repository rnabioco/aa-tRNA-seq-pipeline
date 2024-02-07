import pandas as pd
import numpy as np
import sys
import os

'''
given a tsv file of processed sb tags from a bam file
(e.g., produced from process_sb_with_insertions.py)
generate a new tsv file with a column that takes the mean sb
value over a window of 20 nucleotides

call this script as sb_rolling_means.py processed_sb.tsv
to get a tsv file appended with _prepared
'''


# Check if the file name is provided as a command-line argument
if len(sys.argv) < 2:
    print("Usage: python script.py <filename>")
    sys.exit(1)

# First command-line argument is the script name, second is the filename
file_name = sys.argv[1]

# Check if the provided file exists in the current directory
if not os.path.isfile(file_name):
    print(f"Error: File '{file_name}' not found in the current directory.")
    sys.exit(1)

# Read the data
extracted_data = pd.read_csv(file_name, sep='\t', usecols=['read_id', 'reference_start', 'processed_sb', 'file'])

# Process columns as before
extracted_data['processed_sb'] = extracted_data['processed_sb'].apply(lambda x: [int(i) for i in x.split(',')])
extracted_data[['date', 'sample']] = extracted_data['file'].str.extract(r'(\d{8})_(.*).spb.bam')

# Generate 'reference_position_list' and elongate the DataFrame
extracted_data['reference_position_list'] = extracted_data.apply(
    lambda row: np.arange(row['reference_start'], row['reference_start'] + len(row['processed_sb'])), axis=1)

prepared_data = extracted_data.explode(['processed_sb', 'reference_position_list'])

prepared_data['adjpos'] = prepared_data['reference_position_list'].astype(int) - 130 + 30 + 1

# Convert 'processed_sb' to numeric before rolling
prepared_data['processed_sb'] = prepared_data['processed_sb'].astype(int)

# Take rolling mean of 'processed_sb' over a 20 nt window
prepared_data['rolling_mean_sb'] = prepared_data.groupby('read_id')['processed_sb'] \
                                                .transform(lambda x: x.rolling(window=20, min_periods=1).mean())

# Write the prepared_data to a new TSV file, adding '_prepared' suffix before the extension
base_name, extension = os.path.splitext(file_name)
output_file_path = f"{base_name}_prepared{extension}"
prepared_data.to_csv(output_file_path, sep='\t', index=False)

print(f"Processed file written to: {output_file_path}")

