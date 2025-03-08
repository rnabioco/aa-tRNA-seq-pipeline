import argparse
import pysam
import sys


def get_charging_stats(fn):
    fo = pysam.AlignmentFile(fn)
    read_counts = {}
    for read in fo:
        if not read.is_unmapped:
            if read.reference_name not in read_counts:
                read_counts[read.reference_name] = 0
            read_counts[read.reference_name] += 1

    fo.close()
    return read_counts


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
        Generate a table of read counts and percent aminoacylation. 
        """
    )
    parser.add_argument("-b", "--bam", help="Input BAM file", required=True)

    parser.add_argument(
        "-t",
        "--table",
        required=True,
        help="""
        Path to a table with trna isodecoder and charging status information.
        three whitespace deliminated columns, no header: 
           - sequence name of uncharged tRNA: Name of the uncharged tRNA sequence, must match the fasta entry (e.g. tRNA-Ala-AGC-1-1-uncharged)
           - sequence name of charged tRNA: Name of the charged tRNA sequence, must match the fasta entry (e.g. tRNA-Ala-AGC-1-1-charged)
           - isodecoder: The isodecoder family of the tRNA sequence (e.g. Ala-AGC)
        """,
    )

    parser.add_argument("-o", "--out", help="Path to output tsv file", required=False)

    args = parser.parse_args()
    bam_fl = args.bam

    if args.out:
        fout = open(args.out, "w")
    else:
        fout = sys.stdout

    trna_ref_dict = {}

    with open(args.table, "r") as f:
        for line in f:
            try:
                uncharged, charged, _, gene = line.strip().split()
            except ValueError:
                sys.exit(
                    "tRNA table must have 4 columns (no header): uncharged, charged, isodecoder, gene"
                )

            trna_ref_dict[uncharged] = {"gene": gene, "charge_status": "uncharged"}

            trna_ref_dict[charged] = {"gene": gene, "charge_status": "charged"}

    read_counts = get_charging_stats(bam_fl)
    charged_counts = {}

    for ref, count in read_counts.items():
        try:
            chrg = trna_ref_dict[ref]["charge_status"]
            gene = trna_ref_dict[ref]["gene"]
        except KeyError:
            print("Reference not found in tRNA table: ", ref, file=sys.stderr)
            continue

        if gene not in charged_counts:
            charged_counts[gene] = {"charged": 0, "uncharged": 0}

        charged_counts[gene][chrg] += count

    print(
        "tRNA-gene",
        "charged_read_counts",
        "uncharged_read_counts",
        "percent_charged",
        sep="\t",
        file=fout,
    )
    for gene, counts in charged_counts.items():
        total = counts["charged"] + counts["uncharged"]
        pct_charged = 100 * counts["charged"] / total
        print(
            gene,
            counts["charged"],
            counts["uncharged"],
            pct_charged,
            sep="\t",
            file=fout,
        )
    fout.close()
