import argparse
import pysam
import sys
import gzip

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

    def __init__(self, max=MAX_ARRAY_LENGTH):
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
        uniq_pos = [i for i, v in enumerate(self.counts) if v > 0]
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
        for i, v in enumerate(self.counts):
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

    def __init__(self, max_qlen=MAX_ARRAY_LENGTH, max_qqual=60, max_mapq=255):
        """
        Initializes a ReadStats object.

        Args:
            max_qlen (int, optional): The maximum read length. Defaults to MAX_ARRAY_LENGTH.
            max_qqual (int, optional): The maximum base quality score. Defaults to 60.
            max_mapq (int, optional): The maximum MAPQ score. Defaults to 255.
        """
        self.mapped_reads = 0
        self.n_pos_reads = 0
        self.qlens = CountArray(max=max_qlen)
        self.qquals = CountArray(max=max_qqual)
        self.mapq = CountArray(max=max_mapq)

    def summary(self, q=0.5, digits=3):
        d = {
            "mapped_reads": self.mapped_reads,
            "pos_reads": self.n_pos_reads,
            "mapq0_reads": self.mapq.counts[0],
            "mapq_pass_reads": self.mapq.total - self.mapq.counts[0],
            "mean_length": round(self.qlens.mean(), digits),
            "median_length": self.qlens.quantile(q),
            "mean_base_quality": round(self.qquals.mean(), digits),
            "median_base_quality": self.qquals.quantile(q),
            "mean_MAPQ": round(self.mapq.mean(), digits),
            "median_MAPQ": self.mapq.quantile(q),
        }

        return d


def get_read_stats(fn, flag=None, sample_id=None, sample_info=None):
    n_uniq_reads = 0
    seen_qnames = set()

    stats = ReadStats()

    fo = pysam.AlignmentFile(fn, check_sq=False)

    for read in fo:
        if flag is not None:
            if read.flag & flag != flag:
                continue
            if read.is_secondary or read.is_supplementary:
                continue

        # tracking unique reads will use alot a memory
        # consider using a bloom filter if this becomes an issue
        qname = read.query_name
        if qname in seen_qnames:
            continue

        n_uniq_reads += 1
        seen_qnames.add(qname)

        qlen = read.query_length
        mean_qual = round(
            sum(read.query_alignment_qualities) / len(read.query_alignment_qualities)
        )

        stats.qlens.push(qlen)
        stats.qquals.push(mean_qual)

        if not read.is_unmapped:
            stats.mapped_reads += 1
            stats.mapq.push(read.mapping_quality)
            if not read.is_reverse:
                stats.n_pos_reads += 1

    read_summary = stats.summary()
    read_stat_keys = read_summary.keys()

    read_summary["id"] = sample_id if sample_id is not None else ""
    read_summary["info"] = sample_info if sample_info is not None else ""
    read_summary["n_reads"] = n_uniq_reads
    read_summary["bam_file"] = "stdin" if fn == "-" else fn
    read_summary["pct_mapped"] = 0  # to fill in later

    first_col_order = ["bam_file", "id", "info", "n_reads", "pct_mapped"] + list(
        read_stat_keys
    )
    read_summary = {k: read_summary[k] for k in first_col_order}

    fo.close()
    return read_summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
        Collect read and alignment statistics from unmapped, mapped, and further processed BAM files.
        This tool is designed to summarize read counts across BAM files from a single sample."
        """
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
        "-i",
        "--id",
        help="Sample identifier to append to output",
        required=False,
        default="sample",
    )

    parser.add_argument(
        "-a",
        "--info",
        help="Additional information to append to output, supplied for each input BAM",
        required=False,
        nargs="+",
    )

    parser.add_argument(
        "-b",
        "--bam",
        help="Path to one or more BAM file(s), the first BAM file is assumed to be a uBAM containing only unmapped reads.",
        nargs="+",
    )

    args = parser.parse_args()

    bam_fls = args.bam
    flag = args.flag
    if args.out:
        if args.out.endswith(".gz"):
            fout = gzip.open(args.out, "wt")
        else:
            fout = open(args.out, "w")
    else:
        fout = sys.stdout

    if args.info and len(args.info) != len(bam_fls):
        sys.exit("Number of info fields must match number of BAM files")

    total_reads = 0

    for i, bam in enumerate(bam_fls):
        sample_info = None
        if args.info:
            sample_info = args.info[i]

        read_summaries = get_read_stats(bam, flag, args.id, sample_info)

        if i == 0:
            total_reads = read_summaries["n_reads"]
            cols = list(read_summaries.keys())
            fout.write("\t".join(cols) + "\n")

        if read_summaries["mapped_reads"] > 0:
            read_summaries["pct_mapped"] = round(
                read_summaries["mapped_reads"] / total_reads * 100, 2
            )

        out = "\t".join([str(x) for x in read_summaries.values()])
        fout.write(out + "\n")

    fout.close()
