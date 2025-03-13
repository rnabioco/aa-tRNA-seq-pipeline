import argparse
import pysam
import sys

# Adapter lengths
ADAPTER_5P_LEN = 24  # Total length of the 5' tRNA adapter
ADAPTER_3P_LEN = 40  # Total length of the 3' tRNA adapter
FIVE_P_TRUNCATION = ADAPTER_5P_LEN  # Maximum truncation allowed at the 5' end
THREE_P_TRUNCATION = 23  # Maximum truncation allowed at the 3' end, 23 is the first position of the discriminating 3' adapter region
# Filtering logic to ID "full length" reads that align to adapters but permit 5' end truncation

FILTER_CODES = {
    "5p_trunc": 1 << 0,
    "3p_trunc": 1 << 1,
    "both_trunc": 1 << 2,
    "unmapped": 1 << 3,
    "low_mapq": 1 << 4,
    "negative_strand": 1 << 5,
    "supplemental_secondary": 1 << 6,
    "invalid_multimapping": 1 << 7,
    "too_many_edits_in_adapter": 1 << 8,
}
FILTER_TAG = "zf"


class FilterStats:
    """
    Class to log information about why reads were filtered.

    Attributes:
        n_align (int): Total number of alignments.
        n_filtered (int): Total number of filtered alignments.
        filter_counts (dict): Dictionary to store the count of filtered alignments for each reason.
    """

    n_align = 0
    n_filtered = 0

    def __init__(self):
        self.filter_counts = {k: 0 for k in FILTER_CODES.keys()}

    def log(self, tag):
        """
        Logs the information about the alignment and updates the filter counts.

        Args:
            tag (int): Bit flag representing the reason for filtering the alignment.
        """
        self.n_align += 1
        filtered = False
        for reason, val in FILTER_CODES.items():
            if tag & val:
                filtered = True
                self.filter_counts[reason] += 1
        if filtered:
            self.n_filtered += 1

    def summary(self):
        d = {
            "total_alignments": self.n_align,
            "failed_alignments": self.n_filtered,
            "passed_alignments": self.n_align - self.n_filtered,
        }
        d.update(self.filter_counts)
        out = []
        for k, v in d.items():
            out.append(f"{k} {v}")

        return "\n".join(out)

    def __str__(self):
        return str(self.filter_counts)


def count_adapter_edits(aln, ref_start, ref_end):
    # use get_aligned_pairs to get the aligned positions
    # count number of mismatches, deletions and insertions
    # in the adapter region
    # edits defined as sum of insertions deletions and mismatches

    mismatches = 0
    insertions = 0
    deletions = 0
    aligned_pairs = aln.get_aligned_pairs(with_seq=True, matches_only=False)

    # trim 3' softclipped bases
    for i in reversed(range(len(aligned_pairs))):
        if aligned_pairs[i][1] is not None:
            break
    aligned_pairs = aligned_pairs[: (i + 1)]

    last_ref_pos = aligned_pairs[-1][1]
    # alignment doesn't reach adapter region
    if last_ref_pos < ref_start:
        return None

    in_adapter = False
    for query_pos, ref_pos, seq in aligned_pairs:
        if ref_pos is not None and ref_pos >= ref_end:
            break

        if ref_pos is not None and ref_start <= ref_pos < ref_end:
            in_adapter = True
            if query_pos is None:
                # Deletion in the read
                deletions += 1
            elif seq is not None and seq.islower():
                # Mismatch
                mismatches += 1
        elif ref_pos is None and query_pos is not None and in_adapter:
            # Insertion in the read
            insertions += 1

    # assign soft clipped adapter bases as deletions
    if ref_pos < ref_end:
        n_adapter_sc = ref_end - ref_pos - 1
        deletions += n_adapter_sc

    return insertions + deletions + mismatches


def compatible_secondary_alignments(aln, trna_ref_dict, isodecoder_ref_dict):
    # read doesn't have multiple alignments (at least not in the XA tag)
    # likely due to a larger number of secondary alignments than -h setting in bwa mem
    if not aln.has_tag("XA"):
        return False

    xa_list = aln.get_tag("XA")
    xa_list = xa_list.split(";")
    xa_list = [xa for xa in xa_list if xa != ""]

    primary_aln_charge = trna_ref_dict[aln.reference_name]["charge_status"]
    primary_aln_isodecoder = trna_ref_dict[aln.reference_name]["isodecoder"]
    isodecoder_genes = isodecoder_ref_dict[primary_aln_isodecoder][primary_aln_charge]

    res = True
    for xa in xa_list:
        xa_ref = xa.split(",")[0]

        if trna_ref_dict[xa_ref]["isodecoder"] != primary_aln_isodecoder:
            res = False

        if xa_ref not in isodecoder_genes:
            res = False

    return res


def filter_bam(args):
    five_p_truncation = args.five_p_truncation
    p3_truncation_max = args.three_p_truncation
    min_mapq = args.min_mapq
    only_positive = args.only_positive
    rescue_multi_mappers = args.rescue_multi_mappers

    bamfile = pysam.AlignmentFile(args.input_bam, "rb")
    if args.failed_bam:
        failed_bam = pysam.AlignmentFile(args.failed_bam, "wb", template=bamfile)

    if rescue_multi_mappers:
        # store two dictionaries:
        # isodecoder_ref contains the charged and uncharged tRNAs for each isodecoder
        # trna_ref_dict contains the individual tRNAs and their charging status
        isodecoder_ref = {}
        trna_ref_dict = {}
        if not args.trna_table:
            sys.exit("Rescue multi-mappers (-r) requires a tRNA table (-t)")

        with open(args.trna_table, "r") as f:
            for line in f:
                try:
                    uncharged, charged, isodecoder, _ = line.strip().split()
                except ValueError:
                    sys.exit(
                        "tRNA table must have 4 columns (no header): uncharged, charged, isodecoder, gene"
                    )

                if isodecoder not in isodecoder_ref:
                    isodecoder_ref[isodecoder] = {"charged": set(), "uncharged": set()}
                isodecoder_ref[isodecoder]["charged"].add(charged)
                isodecoder_ref[isodecoder]["uncharged"].add(uncharged)

                trna_ref_dict[uncharged] = {
                    "isodecoder": isodecoder,
                    "charge_status": "uncharged",
                    "key_charge_inverse": charged,
                }

                trna_ref_dict[charged] = {
                    "isodecoder": isodecoder,
                    "charge_status": "charged",
                    "key_charge_inverse": uncharged,
                }

    if args.max_edit_dist:
        max_edit_dist, adapter_region = args.max_edit_dist.split(":")
        max_edit_dist = int(max_edit_dist)
        adapter_region = adapter_region.split("-")
        adapter_region = [int(adapter_region[0]) + 1, int(adapter_region[1])]

    stats = FilterStats()
    with pysam.AlignmentFile(args.output_bam, "wb", template=bamfile) as outfile:
        for read in bamfile:
            tag = 0
            if read.is_unmapped:
                tag |= FILTER_CODES["unmapped"]
                read.set_tag(FILTER_TAG, tag, "i")
                stats.log(tag)
                if args.failed_bam:
                    failed_bam.write(read)
                continue

            if read.is_secondary or read.is_supplementary:
                tag |= FILTER_CODES["supplemental_secondary"]

            if read.mapping_quality < min_mapq:
                tag |= FILTER_CODES["low_mapq"]

            if (min_mapq == 0 and read.mapping_quality == 0) or (
                tag & FILTER_CODES["low_mapq"]
            ):
                if rescue_multi_mappers:
                    if not compatible_secondary_alignments(
                        read, trna_ref_dict, isodecoder_ref
                    ):
                        tag |= FILTER_CODES["invalid_multimapping"]

            if only_positive and read.is_reverse:
                tag |= FILTER_CODES["negative_strand"]

            ref_length = bamfile.get_reference_length(read.reference_name)
            expected_min_end = ref_length - p3_truncation_max

            if five_p_truncation >= 0 and read.reference_start > five_p_truncation:
                tag |= FILTER_CODES["5p_trunc"]

            if p3_truncation_max >= 0 and read.reference_end < expected_min_end:
                tag |= FILTER_CODES["3p_trunc"]
                if tag & FILTER_CODES["5p_trunc"]:
                    tag |= FILTER_CODES["both_trunc"]
                    tag &= ~FILTER_CODES["5p_trunc"]
                    tag &= ~FILTER_CODES["3p_trunc"]

            if args.max_edit_dist:
                edits = count_adapter_edits(
                    read, ref_length - adapter_region[0], ref_length - adapter_region[1]
                )
                if edits is None or edits > max_edit_dist:
                    tag |= FILTER_CODES["too_many_edits_in_adapter"]

            read.set_tag(FILTER_TAG, tag, "i")
            stats.log(tag)
            if not tag:
                outfile.write(read)

            if args.failed_bam and tag:
                failed_bam.write(read)
    print(stats.summary())
    bamfile.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
        Filter reads in a BAM file according to the indicated criteria.

        Note the output BAM file(s) will have a new tag "zf" containing a 
        bitwise flag of the reasons for filtering. The tag values
        are as follows:
        0: Passes filtering
        1: 5' truncation
        2: 3' truncation
        4: both 5' and 3' truncation
        8: unmapped
        16: low MAPQ
        32: negative strand
        64: supplemental or secondary alignment
        128: read has one or more equivalent alignments to charged and uncharged tRNA, or multiple isodecoder families
        256: too many edits in adapter region
        """
    )

    parser.add_argument("-i", "--input_bam", required=True, help="Input BAM file")

    parser.add_argument("-o", "--output_bam", required=True, help="Output BAM file")

    parser.add_argument(
        "-5",
        "--five_p_truncation",
        help="maximum 5' truncation allowed, to -1 to disable",
        type=int,
        default=FIVE_P_TRUNCATION,
    )
    parser.add_argument(
        "-3",
        "--three_p_truncation",
        help="maximum 3' truncation allowed, to -1 to disable",
        type=int,
        default=THREE_P_TRUNCATION,
    )
    parser.add_argument(
        "-q",
        "--min_mapq",
        help="Required MAPQ score for reads to be retained",
        type=int,
        default=0,
    )
    parser.add_argument(
        "-s",
        "--only_positive",
        help="if set, only return reads on positive strand",
        action="store_true",
    )
    parser.add_argument(
        "-r",
        "--rescue_multi_mappers",
        help="""
        if set, parse secondary alignments to rescue reads that map to same tRNA isodecoder family, 
        requires -t --trna-table input
        """,
        action="store_true",
    )
    parser.add_argument(
        "-t",
        "--trna-table",
        help="""
        text file with three whitespace deliminated columns, no header:
          - sequence name of uncharged tRNA: Name of the uncharged tRNA sequence, must match the fasta entry (e.g. tRNA-Ala-AGC-1-1-uncharged)
          - sequence name of charged tRNA: Name of the charged tRNA sequence, must match the fasta entry (e.g. tRNA-Ala-AGC-1-1-charged)
          - isodecoder: The isodecoder family of the tRNA sequence (e.g. Ala-AGC)
        """,
        required=False,
    )

    parser.add_argument(
        "-e",
        "--max-edit-dist",
        help="""
        Exclude reads with higher than max-edit-dist within the adapter region. Edit distance is defined as the # of 
        read mismatches, insertions, and deletions in the adapter region.

        Supplied this as:
        4:30-11 
        [maximum edit distance]:[5' position (1-based, closed) of 3' adapter offset from end of reference]-[5' position (1-based, closed) of 3' adapter offset from end of reference] ]  
        e.g. 
        4:22-10 means exclude reads with more than 4 mismatches in the 3' adapter region from reference length - 22 to reference length - 10.
        internally the 5' position will be converted to 0-based, open to match pysam conventions.
        For v2 charged and uncharged adapter set, use :22-10 to include the 3' adapter region.
        """,
        required=False,
    )

    parser.add_argument(
        "-f",
        "--failed_bam",
        help="if provided write alignements failing filtering to indicated BAM file",
        required=False,
    )
    args = parser.parse_args()

    filter_bam(args)
