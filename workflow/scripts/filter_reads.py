import pysam
import sys
import argparse
import pysam



# Adapter lengths
ADAPTER_5P_LEN = 24  # Total length of the 5' tRNA adapter
ADAPTER_3P_LEN = 30  # Total length of the 3' tRNA adapter
FIVE_P_TRUNCATION = 15  # Maximum truncation allowed at the 5' end (as helicase loses contact with the RNA)
# Filtering logic to ID "full length" reads that align to adapters but permit 5' end truncation

def filter_bam(inbam, outbam, p5_truncation_max, p3_truncation_max, min_mapq, only_positive):
    bamfile = pysam.AlignmentFile(inbam, "rb")
    with pysam.AlignmentFile(outbam, "wb", template=bamfile) as outfile:
        for read in bamfile:
            if read.is_unmapped:
                continue
            
            if read.mapping_quality < min_mapq:
                continue
            
            if only_positive and read.is_reverse:
                continue

            ref_length = bamfile.get_reference_length(read.reference_name)
            expected_min_end = ref_length - p3_truncation_max

            if read.reference_start <= p5_truncation_max and read.reference_end >= expected_min_end:
                outfile.write(read)
    bamfile.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_bam", help="Input BAM file")
    parser.add_argument("-o", "--output_bam", help="Output BAM file")
    parser.add_argument("-5",  "--five_p_truncation", help="maximum 5' truncation allowed", type = int, default=FIVE_P_TRUNCATION)
    parser.add_argument("-3",  "--three_p_truncation", help="maximum 3' truncation allowed", type = int, default=ADAPTER_3P_LEN)
    parser.add_argument("-q",  "--min_mapq", help="Required MAPQ score for reads to be retained", type = int, default = 1)
    parser.add_argument("-s",  "--only_positive", help="if set, only return reads on positive strand", action="store_true")
    args = parser.parse_args()

    filter_bam(args.input_bam, 
               args.output_bam, 
               args.five_p_truncation,
               args.three_p_truncation, 
               args.min_mapq, 
               args.only_positive)
