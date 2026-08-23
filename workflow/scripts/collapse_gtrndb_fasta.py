#!/usr/bin/env python
"""
Collapse redundant sequences in a GtRNAdb FASTA reference.

Parses GtRNAdb FASTA headers to extract tRNA identity, groups sequences
by isodecoder family (amino acid + anticodon), identifies redundant gene
copies with identical sequences, and writes a collapsed (non-redundant)
reference FASTA plus a TSV mapping report.

GtRNAdb header format (variable species prefix):
  >Escherichia_coli_str_K-12_substr_MG1655_tRNA-Ala-GGC-1-1 (tRNAscan-SE ID: ...)

The tRNA portion is extracted with an end-anchored regex, making parsing
robust to arbitrary species prefixes containing hyphens.
"""

import argparse
import csv
import re
import sys


def uf_find(parent, x):
    """Union-find root of `x`, path-compressing `parent` in place.

    Takes `parent` as an argument rather than closing over it: the map is
    rebuilt per group, and a closure over a rebound loop variable is the B023
    footgun even when -- as here -- every call happens in the same iteration.
    """
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


TRNA_NAME_RE = re.compile(
    r"((?:pre)?tRNA)"  # group 1: tRNA or pretRNA
    r"-([A-Za-z0-9]+)"  # group 2: amino acid (Ala, Ile2, fMet, SeC, etc.)
    r"-([A-Z]{3}|NNN)"  # group 3: anticodon
    r"-(\d+)"  # group 4: family number
    r"-(\d+)"  # group 5: gene copy number
    r"$"
)


def read_fasta(fasta_path):
    """Read FASTA file and yield (name, sequence) tuples."""
    name = None
    seq_parts = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_parts).upper()
                name = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)

        if name is not None:
            yield name, "".join(seq_parts).upper()


def write_fasta(sequences, output_path):
    """Write sequences to FASTA format (single line per sequence)."""
    with open(output_path, "w") as f:
        for name, seq in sequences:
            f.write(f">{name}\n")
            f.write(seq + "\n")


def parse_trna_name(header):
    """
    Parse a GtRNAdb FASTA header to extract tRNA name components.

    Returns a dict with keys: trna_type, amino_acid, anticodon,
    family_num, copy_num, prefix, short_name.
    Returns None if header does not match GtRNAdb naming convention.
    """
    m = TRNA_NAME_RE.search(header)
    if m is None:
        return None

    trna_type = m.group(1)
    amino_acid = m.group(2)
    anticodon = m.group(3)
    family_num = m.group(4)
    copy_num = m.group(5)

    short_name = f"{trna_type}-{amino_acid}-{anticodon}-{family_num}-{copy_num}"
    prefix = header[: m.start()]
    prefix = prefix.removesuffix("_")

    return {
        "trna_type": trna_type,
        "amino_acid": amino_acid,
        "anticodon": anticodon,
        "family_num": family_num,
        "copy_num": copy_num,
        "prefix": prefix,
        "short_name": short_name,
    }


def strip_trailing_cca(seq):
    """Strip trailing CCA from a sequence for comparison purposes."""
    if seq.upper().endswith("CCA"):
        return seq[:-3]
    return seq


def hamming_distance(s1, s2):
    """Hamming distance between two equal-length strings.

    Returns None if strings have different lengths.
    """
    if len(s1) != len(s2):
        return None
    return sum(a != b for a, b in zip(s1, s2, strict=True))


def collapse_sequences(records, keep_unparsed=False, max_hamming=0):
    """
    Collapse redundant tRNA sequences.

    Args:
        records: list of (header, sequence) tuples from FASTA
        keep_unparsed: if True, pass through sequences with non-GtRNAdb headers;
                       if False, raise ValueError on unparseable headers
        max_hamming: maximum Hamming distance for merging different families
                     within the same isodecoder. 0 means exact match only.

    Returns:
        (output_sequences, mapping_rows) where:
        - output_sequences: list of (short_name, sequence) for collapsed FASTA
        - mapping_rows: list of dicts for the mapping TSV
    """
    parsed_records = []
    unparsed_records = []
    seen_names = set()

    for header, seq in records:
        if header in seen_names:
            raise ValueError(f"Duplicate input name: {header}")
        seen_names.add(header)

        parsed = parse_trna_name(header)
        if parsed is None:
            if not keep_unparsed:
                raise ValueError(
                    f"Could not parse GtRNAdb header: {header}\n"
                    "Use --keep-unparsed to pass through non-GtRNAdb sequences."
                )
            unparsed_records.append((header, seq))
        else:
            parsed_records.append((header, seq, parsed))

    # Group by isodecoder key
    isodecoder_groups = {}
    for header, seq, parsed in parsed_records:
        key = f"{parsed['trna_type']}-{parsed['amino_acid']}-{parsed['anticodon']}"
        if key not in isodecoder_groups:
            isodecoder_groups[key] = []
        isodecoder_groups[key].append((header, seq, parsed))

    output_sequences = []
    mapping_rows = []
    seen_output_names = set()

    # Process each isodecoder group
    for iso_key in sorted(isodecoder_groups):
        members = isodecoder_groups[iso_key]

        # Group by normalized (CCA-stripped) sequence
        seq_groups = {}
        for header, seq, parsed in members:
            norm_seq = strip_trailing_cca(seq)
            if norm_seq not in seq_groups:
                seq_groups[norm_seq] = []
            seq_groups[norm_seq].append((header, seq, parsed))

        # Hamming-based merging of different sequence groups
        if max_hamming > 0 and len(seq_groups) > 1:
            seqs = list(seq_groups.keys())
            # Union-find for transitive merges
            parent = {s: s for s in seqs}

            for i in range(len(seqs)):
                for j in range(i + 1, len(seqs)):
                    dist = hamming_distance(seqs[i], seqs[j])
                    if dist is not None and dist <= max_hamming:
                        ri = uf_find(parent, seqs[i])
                        rj = uf_find(parent, seqs[j])
                        if ri != rj:
                            # Merge later into earlier (by first appearance)
                            parent[rj] = ri

            # Rebuild seq_groups by merging
            merged = {}
            for s in seqs:
                root = uf_find(parent, s)
                if root not in merged:
                    merged[root] = []
                merged[root].extend(seq_groups[s])
            seq_groups = merged

        for group in seq_groups.values():
            # First encountered is the representative
            rep_header, rep_seq, rep_parsed = group[0]
            rep_short = rep_parsed["short_name"]

            if rep_short not in seen_output_names:
                output_sequences.append((rep_short, rep_seq))
                seen_output_names.add(rep_short)

            for header, seq, parsed in group:
                mapping_rows.append(
                    {
                        "collapsed_name": rep_short,
                        "original_name": header,
                        "isodecoder": iso_key,
                        "amino_acid": parsed["amino_acid"],
                        "anticodon": parsed["anticodon"],
                        "family_num": parsed["family_num"],
                        "copy_num": parsed["copy_num"],
                        "is_representative": header == rep_header,
                        "sequence_length": len(seq),
                    }
                )

    # Append unparsed records as-is
    for header, seq in unparsed_records:
        output_sequences.append((header, seq))
        mapping_rows.append(
            {
                "collapsed_name": header,
                "original_name": header,
                "isodecoder": "unparsed",
                "amino_acid": "",
                "anticodon": "",
                "family_num": "",
                "copy_num": "",
                "is_representative": True,
                "sequence_length": len(seq),
            }
        )

    return output_sequences, mapping_rows


def write_mapping(mapping_rows, output_path):
    """Write mapping report as TSV."""
    fieldnames = [
        "collapsed_name",
        "original_name",
        "isodecoder",
        "amino_acid",
        "anticodon",
        "family_num",
        "copy_num",
        "is_representative",
        "sequence_length",
    ]

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(mapping_rows)


def main():
    parser = argparse.ArgumentParser(
        description="Collapse redundant sequences in a GtRNAdb FASTA reference",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  python collapse_gtrndb_fasta.py -i ecoliK12-mature-tRNAs.fa \\
      -o ecoliK12-collapsed.fa -m ecoliK12-mapping.tsv

  python collapse_gtrndb_fasta.py -i reference.fa \\
      -o collapsed.fa -m mapping.tsv --keep-unparsed
""",
    )

    parser.add_argument("-i", "--input", required=True, help="Input GtRNAdb FASTA file")
    parser.add_argument(
        "-o", "--output", required=True, help="Output collapsed FASTA file"
    )
    parser.add_argument(
        "-m", "--mapping", required=True, help="Output TSV mapping report"
    )
    parser.add_argument(
        "--keep-unparsed",
        action="store_true",
        default=False,
        help="Pass through sequences with non-GtRNAdb headers (default: error)",
    )
    parser.add_argument(
        "--max-hamming",
        type=int,
        default=0,
        help="Max Hamming distance for merging same-anticodon families "
        "(0 = exact match only, default: 0)",
    )

    args = parser.parse_args()

    # Read input
    records = list(read_fasta(args.input))
    if not records:
        print("Error: Input FASTA is empty.", file=sys.stderr)
        sys.exit(1)

    # Collapse
    try:
        output_sequences, mapping_rows = collapse_sequences(
            records,
            keep_unparsed=args.keep_unparsed,
            max_hamming=args.max_hamming,
        )
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    # Count stats
    n_input = len(records)
    n_isodecoders = len(
        {r["isodecoder"] for r in mapping_rows if r["isodecoder"] != "unparsed"}
    )
    n_output = len(output_sequences)
    n_collapsed = n_input - n_output

    # Write outputs
    write_fasta(output_sequences, args.output)
    write_mapping(mapping_rows, args.mapping)

    # Summary to stderr
    print(
        f"Input: {n_input} sequences from {n_isodecoders} isodecoder families",
        file=sys.stderr,
    )
    print(f"Output: {n_output} unique sequences", file=sys.stderr)
    print(f"Collapsed: {n_collapsed} redundant gene copies removed", file=sys.stderr)


if __name__ == "__main__":
    main()
