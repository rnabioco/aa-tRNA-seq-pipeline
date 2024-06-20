import argparse
import pysam
import sys
from collections import OrderedDict

MAX_ARRAY_LENGTH = 100000

class CountArray:
    """
    A class to represent integer count frequencies. Memory efficent approach
    to store various integer metrics (read lengths, quality scores, etc.)

    Attributes:
        counts (list): A list of counts.
        total (int): The total count.

    Methods:
        push(index): Increment the count at the given index.
        nth(n): Get the nth element from count distribution.
        mean(): Calculate the mean value of count distribution.
        quantile(p): Calculate the quantile value at the given percentile.
    """

    def __init__(self, max = MAX_ARRAY_LENGTH):
        self.counts = [0] * max
        self.total = 0
        self.len = len(self.counts) 

    def push(self, index):
        if index >= self.len:
            return 
        self.counts[index] += 1
        self.total += 1
    
    def nth(self, n):
        """
        Get the nth value from array frequency distribution
        Args:
            n (int): The index of the element to retrieve.

        Returns:
            int or None: The value at the nth index, or None if index is out of range.
        """
        start_idx = 0
        uniq_pos = [i for i,v in enumerate(self.counts) if v > 0]
        for pos in uniq_pos:
            freq = self.counts[pos]
            end_idx = start_idx + freq
            if n < end_idx:
                return pos
            else:
                start_idx = end_idx 
        return None
    
    def mean(self):
        if self.total == 0:
            return 0   
        
        value_total = 0
        for i,v in enumerate(self.counts):
            value_total += i * v
        return value_total / self.total
    
    def quantile(self, p):
        if self.total == 0:
            return 0   
        
        if p < 0 or p > 1:
            raise ValueError("p must be between 0 and 1")
        
        idx = self.total * p
        if idx % 1 == 0:
            i = int(idx)
            val = (self.nth(i - 1) + self.nth(i)) / 2
        else: 
            val = self.nth(int(idx))
        return val

class ReadStats:
    """
    Class to calculate statistics for a set of reads.

    Attributes:
        read_type (str): The type of reads (one of "mapped" or "unmapped").
        n_reads (int): The total number of reads.
        n_pos_reads (int): The number of positive strand reads.
        qlens (CountArray): An instance of CountArray to store read lengths.
        qquals (CountArray): An instance of CountArray to store base qualities.
        mapq (CountArray): An instance of CountArray to store MAPQ scores.

    Methods:
        summary(q=0.5, digits=3): Returns a summary dictionary of read statistics.
    """

    def __init__(self, read_type="mapped", max_qlen=MAX_ARRAY_LENGTH, max_qqual=60, max_mapq=255):
        """
        Initializes a ReadStats object.

        Args:
            read_type (str, optional): The type of reads. Defaults to "mapped".
            max_qlen (int, optional): The maximum read length. Defaults to MAX_ARRAY_LENGTH.
            max_qqual (int, optional): The maximum base quality score. Defaults to 60.
            max_mapq (int, optional): The maximum MAPQ score. Defaults to 255.
        """
        self.read_type = read_type
        self.n_reads = 0
        self.n_pos_reads = 0
        self.qlens = CountArray(max=max_qlen)
        self.qquals = CountArray(max=max_qqual)
        self.mapq = CountArray(max=max_mapq)

    def summary(self, q=0.5, digits=3):
        d = {
            "reads": self.n_reads,
            "pos_reads": self.n_pos_reads,
            "mean_length": round(self.qlens.mean(), digits),
            "median_length": self.qlens.quantile(q),
            "mean_base_quality": round(self.qquals.mean(), digits),
            "median_base_quality": self.qquals.quantile(q),
            "mean_MAPQ": round(self.mapq.mean(), digits),
            "median_MAPQ": self.mapq.quantile(q)
        }
        if self.read_type == "unmapped":
            d.pop("pos_reads")
            d.pop("mean_MAPQ")
            d.pop("median_MAPQ")
        return {f'{self.read_type}_{k}': v for k, v in d.items()}


def get_read_stats(fn, flag = None, id = None):
    n_uniq_reads = 0
    seen_qnames = set()

    m_stats = ReadStats(read_type = "mapped")
    um_stats = ReadStats(read_type = "unmapped")

    fo = pysam.AlignmentFile(fn, check_sq=False)
    
    for read in fo:
        if flag is not None:
            if read.flag & flag != flag:
                continue
        # tracking unique reads will use alot a memory 
        # consider using a bloom filter if this becomes an issue 
        qname = read.query_name        
        if qname in seen_qnames:
            continue

        n_uniq_reads += 1
        seen_qnames.add(qname)
        
        qlen = read.query_length
        mean_qual = round(sum(read.query_alignment_qualities) / len(read.query_alignment_qualities))

        if read.is_unmapped:
            um_stats.n_reads += 1
            um_stats.qlens.push(qlen)
            um_stats.qquals.push(mean_qual)
        else:
            m_stats.n_reads += 1
            m_stats.qlens.push(qlen)
            m_stats.qquals.push(mean_qual)
            m_stats.mapq.push(read.mapping_quality)
            if not read.is_reverse:
                m_stats.n_pos_reads += 1

    um_summary = um_stats.summary()
    m_summary = m_stats.summary()
    read_summary = {**um_summary , **m_summary}
    
    read_summary["id"] = id if id is not None else ""
    read_summary["n_reads"] = n_uniq_reads
    read_summary["bam_file"] = "stdin" if fn == "-" else fn

    first_col_order = ["bam_file", "id", "n_reads"] +  list(m_summary.keys()) + list(um_summary.keys())
    read_summary = {k: read_summary[k] for k in first_col_order}

    fo.close()
    return read_summary

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Collect read and alignment statistics from BAM files."
    )
    parser.add_argument(
        "-f",
        "--flag",
        help="Require that mapped reads have the indicated BAM flag set",
        required=False,
        type=int,
    )
    parser.add_argument("-o", "--out", help="Path to output tsv file", required=False)

    parser.add_argument(
        "-n",
        "--id",
        help="Additional id information to append to output, supplied for each input BAM",
        required=False,
        nargs="+"
    )

    parser.add_argument(
        "-b", 
        "--bam", 
        help="Path to one or more BAM file(s)",
        nargs="+"
    )

    args = parser.parse_args()

    bam_fls = args.bam
    flag = args.flag
    if args.out:
        fout = open(args.out, 'w')
    else: 
        fout = sys.stdout
    
    if args.id and len(args.id) != len(bam_fls):
        sys.exit("Number of IDs must match number of BAM files")

    for i, bam in enumerate(bam_fls):
        id = None
        if args.id:
            id = args.id[i]

        read_summaries = get_read_stats(bam, flag, id)
        if i == 0:
            cols = list(read_summaries.keys())
            fout.write("\t".join(cols) + "\n")

        out = "\t".join([str(x) for x in read_summaries.values()])
        fout.write(out + "\n")

    fout.close()
