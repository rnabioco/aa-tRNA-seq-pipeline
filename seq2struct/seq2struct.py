import csv
import argparse

class Seq2StructTSV:
    def __init__(self, fasta_file, fsa_file, output_file):
        self.fasta_file = fasta_file
        self.fsa_file = fsa_file
        self.output_file = output_file

    def read_fasta(self, file_name):
        sequences, headers = {}, []
        with open(file_name, 'r') as file:
            seq_id = None
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    seq_id = line[1:]
                    headers.append(seq_id)
                    sequences[seq_id] = ''
                else:
                    sequences[seq_id] += line.upper()
        return sequences, headers

    def create_mapping(self, ref_seq, ann_seq):
        mapping, ref_index = {}, 0
        for ann_index, ann_nuc in enumerate(ann_seq):
            mapping[ann_index + 1] = ref_index + 1 if ann_nuc != '-' else None
            if ann_nuc != '-':
                ref_index += 1
        return mapping

    def process(self):
        ref_seqs, _ = self.read_fasta(self.fasta_file)
        ann_seqs, ann_headers = self.read_fasta(self.fsa_file)
        with open(self.output_file, 'w', newline='') as file:
            writer = csv.writer(file, delimiter='\t')
            writer.writerow(['seq_ref', 'struct_ref', 'struct_nt', 'struct_pos', 'seq_pos'])
            for ann_id in ann_headers:
                parts = ann_id.split('-')
                prefix, amino_acid, anticodon_fsa = parts[0], parts[2], parts[3].replace('T', 'U')
                matching_refs = [ref_id for ref_id in ref_seqs if ref_id.startswith(f'{prefix}-tRNA-{amino_acid}-{anticodon_fsa}')]
                for ref_id in matching_refs:
                    mapping = self.create_mapping(ref_seqs[ref_id], ann_seqs[ann_id])
                    for ann_index, ann_nuc in enumerate(ann_seqs[ann_id]):
                        writer.writerow([ref_id, ann_id, ann_nuc, ann_index + 1, mapping.get(ann_index + 1, '')])

# Argument Parsing
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate Sequence-to-Structure TSV")
    parser.add_argument("fasta_file", type=str, help="Path to reference FASTA file")
    parser.add_argument("fsa_file", type=str, help="Path to annotated FSA file")
    parser.add_argument("output_file", type=str, help="Path to output TSV file")
    args = parser.parse_args()

    processor = Seq2StructTSV(args.fasta_file, args.fsa_file, args.output_file)
    processor.process()
    print(f"Sequence-to-Structure TSV saved to {args.output_file}")