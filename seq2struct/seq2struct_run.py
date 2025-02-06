import subprocess
import sys
"""
Currently, this is two classes that provide structure to making the seq to structure mapping file.
The function of this driver script: 
1. takes an afa file from Modomics with annotated modification sites for various organisms and provides a nucleotide readible in ATGC.
2. Stockholm files of the reference tRNAs are converted to afa files prior to this script. 
these reference files are anchored to a consensus structure. This script will take these anchored reference sequences and 
provide a structural position to each nucleotide.
3. then both the mod.afa and ref.afa will be used to give each nucleotide in the fasta reference a structural position and modification annotation (if known)


Workflow:
1. take a stockholm file reformated to afa and define special characters(?)
2. count non-gap / space positions with corresponding nucleotide
3. save this consensus secondary structure anchored mapping in a new file
4. join the reference with the map to generate a sequence to structure file that can be used for nanopore sequencing conversion of sequence space to structure space.

File inputs:
ref.afa -> reference file anchored to consensus secondary structure converted from .sto to .afa (esl-reformat)
mod.afa -> modomics afa file that contains modification annotations across multiple species (could just start with ecoli? Could get rid of entirely)
ref.fasta -> fasta reference that contains 5' and 3' sequencing adapters

Outputs:
standardized_afa: a afa file that has standardized header format, sequencing adapters added to each sequence, and modification annotations
seq2struct_tsv: map of the secondary structure to reference file with modification annotation at relevant nucleotides.

"""
'---------------------------------------------------------------'  

def main():
    if len(sys.argv) < 2:
        print("Usage: python main.py <script> [args]")
        print("Scripts:")
        print("  standardize_afa  - Standardize AFA files")
        print("  seq2struct_tsv   - Generate Sequence-to-Structure TSV")
        sys.exit(1)

    script = sys.argv[1]
    args = sys.argv[2:]

    if script == "standardize_afa":
        subprocess.run(["python", "standardize_afa.py"] + args)
    elif script == "seq2struct_tsv":
        subprocess.run(["python", "seq2struct_tsv.py"] + args)
    else:
        print(f"Error: Unknown script '{script}'")
        sys.exit(1)

if __name__ == "__main__":
    main()

