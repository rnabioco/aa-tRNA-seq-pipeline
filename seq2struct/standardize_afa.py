import argparse
import tempfile
import shutil
from Bio import SeqIO

"""
TODO:
1. rename the file output to include 'clean header' for clarity
2. include checks on header to make sure it does not change non-fasta headers
3. run tests to check for bugs


input: 
1. afa: reference file generated from a cmalignment (esl-reformat .sto -> .afa)

x. increase flexibility in input file formats resulting in handling header line differently

"""
class StandardizeAFA:
    
    
    def __init__(self, afa):
        """ Initialize inputs afa """
        self.afa = afa
        

    # this function is now fixed
    def clean_headers(self):
        """ Standardize the FASTA headers, add adapters, and replace afa with standardized afa """
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
            temp_filename = temp_file.name  # Get temp file path

            with open(self.afa) as infile, open(self.output_file, 'w') as outfile:
                for record in SeqIO.parse(infile, "fasta"):

                    # remove most irrelevant info from header
                    parts = record.description.split()
                    type_, amino_acid, anticodon = parts[0], parts[3], parts[4].replace('U', 'T').replace('(', '').replace(')', '')
                    header = f'>{type_}-{amino_acid}-{anticodon}'

                    # add adapters to seq
                    modified_seq = f'CCUAAGAGCAAGAAGAAGCCUGGN{record.seq}GGCUUCUUCUUGCUCUUAGGAAAAAAAAAA'

                    # write to new file
                    temp_file.write(header + '\n' + modified_seq + '\n')
        shutil.move(temp_filename, self.afa)
        

# Argument Parsing
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Standardize AFA files")
    parser.add_argument("input_file", type=str, help="Path to input AFA file")
    parser.add_argument("output_file", type=str, help="Path to output AFA file")
    args = parser.parse_args()

    processor = StandardizeAFA(args.input_file, args.output_file)
    processor.process()
    print(f"Standardized AFA saved to {args.output_file}")