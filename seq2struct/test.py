import csv
import gzip
from Bio import SeqIO

## TODO
# make sure the files have adapters (returning correctly)
# make sure the NULL is explained by the sequencing adapters.

def add_adapters(infile, out_file):
    """ 
    Keep headers the same so it is easy to match to reference file.
    add adapters
    Replace old AFA file with new AFA file.
    """
    with open(infile) as infile, open(out_file, 'w') as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            # Keep header the same
            header = record.description

            # Add adapters to afa file sequence
            modified_seq = f'CCUAAGAGCAAGAAGAAGCCUGGN{record.seq}GGCUUCUUCUUGCUCUUAGGAAAAAAAAAA'
            
            # Write to new file
            outfile.write(f">{header}\n{modified_seq}\n")
    
    return outfile

def header_seq_dictonary(file_name):
    """Reads a FASTA file (gzipped or plain) and returns sequences and headers."""
    ref_sequences, ref_headers = {}, []

    with open(file_name, 'rt') as file:
        seq_id = None
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                seq_id = line[1:]
                ref_headers.append(seq_id)
                ref_sequences[seq_id] = ''
            else:
                ref_sequences[seq_id] += line.upper()

    return ref_sequences, ref_headers

def create_mapping(ref_seq, afa_seq):
    """Create a mapping of annotation sequence positions to reference sequence positions."""
    mapping = {}
    ref_index = 0  # Tracks position in ref_seq

    for afa_index, afa_nuc in enumerate(afa_seq):
        if afa_nuc != '-':
            while ref_index < len(ref_seq) and ref_seq[ref_index] == '-':
                ref_index += 1  

            if ref_index < len(ref_seq):
                mapping[afa_index + 1] = ref_index + 1
                ref_index += 1  # Move to the next ref_seq position
            else:
                mapping[afa_index + 1] = None
        else:
            mapping[afa_index + 1] = None

    return mapping

def process(fasta_file, fsa_file, output_file):
    """Process sequences and write mappings to a TSV file."""
    # 1. add adpaters to both the fasta and afa sequences
    fasta_file = add_adapters(fasta_file)
    fsa_file = add_adapters(fsa_file)

    # 2. create a dictionary for reference file and afa file containing 
    # headers (keys) and sequences (dictionary)
    ref_seqs, ref_headers = header_seq_dictonary(fasta_file)
    afa_seqs, afa_headers = header_seq_dictonary(fsa_file)

    # 3. write the tsv file using the create_mappings function
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow(['ref_id', 'afa_ref', 'afa_nt', 'afa_pos', 'struct_pos'])

        for ann_id in afa_headers:
            matching_refs = [ref_id for ref_id in ref_seqs if ref_id == ann_id]

            for ref_id in matching_refs:
                mapping = create_mapping(ref_seqs[ref_id], afa_seqs[ann_id])

                for ann_index, ann_nuc in enumerate(afa_seqs[ann_id]):
                    writer.writerow([ref_id, ann_id, ann_nuc, ann_index + 1, mapping.get(ann_index + 1, '')])

# Files:
fasta_file = "test_data/ecoli-t4-trna.fasta.gz"
fsa_file = "test_data/ecoli-t4-trna.afa"
output_file = "test-seq2struct.tsv"

process(fasta_file, fsa_file, output_file)
