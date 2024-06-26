import argparse
import pysam


# Adapter lengths
ADAPTER_5P_LEN = 24  # Total length of the 5' tRNA adapter
ADAPTER_3P_LEN = 30  # Total length of the 3' tRNA adapter
FIVE_P_TRUNCATION = 15  # Maximum truncation allowed at the 5' end (as helicase loses contact with the RNA)
# Filtering logic to ID "full length" reads that align to adapters but permit 5' end truncation

FILTER_CODES = {
    "5p_trunc": 1 << 0,
    "3p_trunc": 1 << 1,
    "both_trunc": 1 << 2,
    "unmapped": 1 << 3,
    "low_mapq": 1 << 4,
    "negative_strand": 1 << 5,
    "supplemental_secondary": 1 << 6,
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
        d = {"total_alignments": self.n_align,
             "failed_alignments": self.n_filtered,
             "passed_alignments": self.n_align - self.n_filtered}
        d.update(self.filter_counts)
        out = []
        for k,v in d.items():
            out.append(f"{k} {v}")   
       
        return  "\n".join(out)

    def __str__(self):
        return str(self.filter_counts)
    
def filter_bam(args):

    five_p_truncation = args.five_p_truncation
    p3_truncation_max = args.three_p_truncation
    min_mapq = args.min_mapq
    only_positive = args.only_positive

    bamfile = pysam.AlignmentFile(args.input_bam, "rb")
    if args.failed_bam:
        failed_bam = pysam.AlignmentFile(args.failed_bam, "wb", template=bamfile)

    stats = FilterStats() 
    with pysam.AlignmentFile(args.output_bam, "wb", template=bamfile) as outfile:
        for read in bamfile:
            tag = 0
            if read.is_unmapped:
                tag |= FILTER_CODES["unmapped"]

            if read.is_secondary or read.is_supplementary:
                tag |= FILTER_CODES["supplemental_secondary"]

            if read.mapping_quality < min_mapq:
                tag |= FILTER_CODES["low_mapq"]
            
            if only_positive and read.is_reverse:
                tag |= FILTER_CODES["negative_strand"]

            ref_length = bamfile.get_reference_length(read.reference_name)
            expected_min_end = ref_length - p3_truncation_max

            if read.reference_start > five_p_truncation:
                tag |= FILTER_CODES["5p_trunc"]
            if read.reference_end < expected_min_end:
                tag |= FILTER_CODES["3p_trunc"]
                if tag & FILTER_CODES["5p_trunc"]:
                    tag |= FILTER_CODES["both_trunc"] 
                    tag &= ~FILTER_CODES["5p_trunc"]
                    tag &= ~FILTER_CODES["3p_trunc"]
            
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
        """
    )

    parser.add_argument(
        "-i",
        "--input_bam",
        help="Input BAM file")
    
    parser.add_argument(
        "-o",
        "--output_bam",
        help="Output BAM file")
    
    parser.add_argument(
        "-5",
        "--five_p_truncation",
        help="maximum 5' truncation allowed",
        type=int,
        default=FIVE_P_TRUNCATION,
    )
    parser.add_argument(
        "-3",
        "--three_p_truncation",
        help="maximum 3' truncation allowed",
        type=int,
        default=ADAPTER_3P_LEN,
    )
    parser.add_argument(
        "-q",
        "--min_mapq",
        help="Required MAPQ score for reads to be retained",
        type=int,
        default=1,
    )
    parser.add_argument(
        "-s",
        "--only_positive",
        help="if set, only return reads on positive strand",
        action="store_true",
    )
    parser.add_argument(
        "-f",
        "--failed_bam",
        help="if provided write alignements failing filtering to indicated BAM file",
        required=False,
    )
    args = parser.parse_args()

    filter_bam(args)
