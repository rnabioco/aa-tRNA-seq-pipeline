import pysam
import sys

# CLI arguments
input_bam = sys.argv[1]
output_bam = sys.argv[2]

# Open the BAM file
bamfile = pysam.AlignmentFile(input_bam, "rb")

# Adapter lengths
adapter_5p_len = 24  # Total length of the 5' tRNA adapter
adapter_3p_len = 30  # Total length of the 3' tRNA adapter
truncation_max = 15  # Maximum truncation allowed at the 5' end (as helicase loses contact with the RNA)

# Filtering logic to ID "full length" reads that align to adapters but permit 5' end truncation
with pysam.AlignmentFile(output_bam, "wb", template=bamfile) as outfile:
    for read in bamfile:
        ref_length = bamfile.get_reference_length(read.reference_name)
        expected_min_end = ref_length - adapter_3p_len
        
        if read.reference_start <= truncation_max and read.reference_end >= expected_min_end:
            outfile.write(read)

bamfile.close()