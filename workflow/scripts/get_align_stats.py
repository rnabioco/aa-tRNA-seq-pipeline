
import argparse
import pysam
import os
from pathlib import Path

# tracking unique reads will use alot a memory 
# consider using a bloom filter if this becomes an issue 
def get_read_stats(fn, unmapped_bam = False):
    n_uniq_reads = 0
    n_pos_reads = 0
    seen_qnames = set()
    mean_qlen = 0
    if unmapped_bam:
        fo = pysam.AlignmentFile(fn, "rb", check_sq=False)
    else:
        fo = pysam.AlignmentFile(fn, "rb") 
    
    for read in fo:
        if not unmapped_bam and read.is_unmapped:
            continue

        qname = read.query_name        
        if qname in seen_qnames:
            continue

        n_uniq_reads += 1
        seen_qnames.add(qname)
        
        if not unmapped_bam and not read.is_reverse:
            n_pos_reads += 1

        qlen = read.query_length
        mean_qlen = mean_qlen + (qlen - mean_qlen) / n_uniq_reads


    fo.close()
    return n_uniq_reads, round(mean_qlen, 6), n_pos_reads 


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge BAM files.')
    parser.add_argument('-u', '--ubam', help='Path to Unmapped BAM file', required = True)
    parser.add_argument('-m', '--mbam', help='Path to Mapped BAM file', required = True)
    parser.add_argument('-o', '--out', help='Path to output tsv file', required = True)
    args = parser.parse_args()

    ubam = Path(args.ubam)
    mbam = Path(args.mbam)
    total_reads, um_read_length, _ = get_read_stats(ubam, unmapped_bam=True)
    m_reads, m_read_length, m_pos_strand = get_read_stats(mbam, unmapped_bam=False)

    stats = [ubam, total_reads, um_read_length, mbam, m_reads, m_read_length, 100 * m_reads / total_reads, 100 * m_pos_strand / m_pos_strand]
    cols = ["unmapped_bam", "total_reads", "mean_read_length", "mapped_bam", "mapped_reads", "mean_mapped_read_length", "pct_aligned", "pct_pos_strand"]

    fout = open(Path(args.out), 'w')
    fout.write("\t".join(cols) + "\n")

    out = "\t".join([str(x) for x in stats])
    fout.write(out + "\n")
    fout.close()