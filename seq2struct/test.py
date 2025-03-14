import csv
import os
import gzip
from Bio import SeqIO

def add_adapters(infile, is_gzipped=False):
    """ 
    Keep headers the same but add adapters to sequences.
    Preserves the original file extension.
    """
    # Extract base name and extension
    base_name, extension = os.path.splitext(os.path.basename(infile))
    if extension == '.gz':
        base_name, extension = os.path.splitext(base_name)
    
    # Generate output filename with original extension
    out_file = f"{base_name}_modified{extension}"
    
    # Handle gzipped or regular files
    open_func = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'
    
    with open_func(infile, mode) as infile, open(out_file, 'w') as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            # Keep header the same
            header = record.description
            
            modified_seq = f'CCTAAGAGCAAGAAGAAGCCTGGN{record.seq}GGCTTCTTCTTGCTCTTAGGAAAAAAAAAA'
            
            # Write to new file
            outfile.write(f">{header}\n{modified_seq}\n")
    
    return out_file

def header_seq_dictionary(file_name, is_gzipped=False):
    """Reads a FASTA file (gzipped or plain) and returns sequences and headers."""
    ref_sequences, ref_headers = {}, []
    
    open_func = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'
    
    with open_func(file_name, mode) as file:
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

def create_mapping(ref_seq, struct_seq):
    """
    Create a bidirectional mapping between reference sequence and structural sequence.
    Also track the actual nucleotides at each position.
    Skip positions where the reference has a gap ('-' or '.').
    """
    # Keep track of both nucleotides and positions
    mapping = {}
    ref_index = 0  # Position in reference sequence (0-indexed)
    ref_pos_counter = 1  # Actual position counter (1-indexed, skips gaps)
    
    for struct_index, struct_nuc in enumerate(struct_seq):
        if struct_nuc not in ['-', '.']:  # If not a gap in structural sequence
            # Skip gaps in reference sequence
            while ref_index < len(ref_seq) and (ref_seq[ref_index] == '-' or ref_seq[ref_index] == '.'):
                ref_index += 1
                
            if ref_index < len(ref_seq):
                # Store all relevant information for this position
                mapping[struct_index + 1] = {
                    'ref_pos': ref_pos_counter,
                    'ref_nt': ref_seq[ref_index] if ref_index < len(ref_seq) else None,
                    'struct_nt': struct_nuc
                }
                ref_index += 1
                ref_pos_counter += 1  # Only increment position counter for non-gap characters
            else:
                # Reached end of reference sequence or part of the adapter
                mapping[struct_index + 1] = {
                    'ref_pos': None,
                    'ref_nt': None,
                    'struct_nt': struct_nuc
                }
        else:
            # Gap in structural sequence
            mapping[struct_index + 1] = {
                'ref_pos': None,
                'ref_nt': None,
                'struct_nt': struct_nuc
            }

    return mapping

def process(fasta_file, fsa_file, output_file):
    """Process sequences and write mappings to a TSV file."""
    # 1. Determine if files are gzipped
    fasta_is_gzipped = fasta_file.endswith('.gz')
    fsa_is_gzipped = fsa_file.endswith('.gz')
    
    # 2. Add adapters to both files
    modified_fasta = add_adapters(fasta_file, fasta_is_gzipped)
    modified_fsa = add_adapters(fsa_file, fsa_is_gzipped)

    # 3. Create dictionaries for reference and structural files
    ref_seqs, ref_headers = header_seq_dictionary(modified_fasta, False)
    struct_seqs, struct_headers = header_seq_dictionary(modified_fsa, False)

    # 4. Write the tsv file with new format
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow(['ref_id', 'struct_id', 'ref_nt', 'struct_nt', 'ref_pos', 'struct_pos'])

        for struct_id in struct_headers:
            matching_refs = [ref_id for ref_id in ref_seqs if ref_id == struct_id]

            for ref_id in matching_refs:
                mapping = create_mapping(ref_seqs[ref_id], struct_seqs[struct_id])

                for struct_pos in sorted(mapping.keys()):
                    map_data = mapping[struct_pos]
                    writer.writerow([
                        ref_id, 
                        struct_id, 
                        map_data['ref_nt'] or '', 
                        map_data['struct_nt'], 
                        map_data['ref_pos'] or '', 
                        struct_pos
                    ])

# Files:
fasta_file = "test_data/ec-t4-trna.fasta"
fsa_file = "test_data/ec-t4-trna.afa"
output_file = "test_result/test-seq2struct.tsv"

process(fasta_file, fsa_file, output_file)