import argparse
import shutil
import pysam
import os
import sys
import gzip
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.metrics import precision_recall_curve, PrecisionRecallDisplay


'''
This script is an implementation of the pairwise alignment scheme for tRNAs described by:

Yu Sun, Michael Piechotta, Isabel Naarmann-de Vries, Christoph Dieterich, Ann E Ehrenhofer-Murray, 
Detection of queuosine and queuosine precursors in tRNAs by direct RNA sequencing,
Nucleic Acids Research, Volume 51, Issue 20, 10 November 2023, Pages 11197–11212,
https://doi.org/10.1093/nar/gkad826

This python code is a reimplementation of their pipeline provided here:
https://github.com/dieterich-lab/QutRNA
'''

PARASAIL_OPTS =  "-a sw_trace_striped_sse41_128_16 -M 1 -X 1 -e 1 -o 1 -c 20 -x -d -O SAMH".split() 
  
def gzwrite(f, compression = 6):
    if f.endswith(".gz"):
        fh = gzip.open(f, 'wt', compresslevel = compression) 
    else: 
        fh = open(f, 'w')
    return fh  

def filter_bam(bam, outbam, aln_score_threshold = 50):
    """
    Filter bam file based on alignment score cutoff. Note that the alignment
    scores are assumed to be in the MAPQ field
    """
    bam_fh = pysam.AlignmentFile(bam, "rb")
    out_fh = pysam.AlignmentFile(outbam, "wb", template = bam_fh)
    for read in bam_fh.fetch():
        if read.mapping_quality >= aln_score_threshold:
            out_fh.write(read)

    bam_fh.close()
    out_fh.close()

def threshold_alignments(og_bam, randomized_bam, out_prefix, p= 0.95):
    """
    Compare alignment scores between true and random alignments and determine
    score cutoff that results in a precision > p 
    """
    aln_scores = []
    rand_aln_scores = []

    with pysam.AlignmentFile(og_bam, 'r') as og_fh:
        for read in og_fh.fetch():
            aln_scores.append(read.mapping_quality)

    with pysam.AlignmentFile(randomized_bam, 'r') as rnd_fh:
        for read in rnd_fh.fetch():
            rand_aln_scores.append(read.mapping_quality)
    
    n_true = len(aln_scores)
    n_rand = len(rand_aln_scores)

    aln_scores = np.array(aln_scores)
    rand_aln_scores = np.array(rand_aln_scores)
    aln_scores = np.append(aln_scores, rand_aln_scores)

    pos_labels = np.ones(n_true)
    neg_labels = np.zeros(n_rand)
    labels = np.append(pos_labels, neg_labels)

    precision, recall, thresholds = precision_recall_curve(labels, aln_scores)
    
    disp = PrecisionRecallDisplay(precision = precision, recall = recall)
    disp.plot()
    plt.savefig(out_prefix + '.pr_curve.pdf')

    with open(out_prefix + ".pr_values.tsv", 'w') as fo:
        fo.write("precision\trecall\tthreshold\n")
        for i in range(len(precision)):
            # the last precision and recall values do not have a threshold and are 1 and 0 respectively
            try:
                thres = thresholds[i]
            except IndexError:
                thres = "NA"
            fo.write("{}\t{}\t{}\n".format(precision[i], recall[i], thres))
    
    try:
        cutoff = np.min(thresholds[precision[:-1] > p])
    except:
        print("Unable to determine alignment score cutoff",
              file = sys.stderr)
        print("Setting to alignment score cutoff to 0",
              file = sys.stderr)
        cutoff = 0

    return cutoff 

def run_parasail(fastq_fn, bam_out, fasta_ref, threads, mem_budget,
        min_align_score = 10):
    
    cmd = [f"gunzip -c {fastq_fn} | "]
    cmd += ["parasail_aligner"]
    cmd += PARASAIL_OPTS

    input_opts = f"-r {mem_budget} -t {threads} -f {fasta_ref}".split()
    cmd += input_opts

    min_align_score = int(min_align_score)
    filter_cmd = f"\"[AS]>={min_align_score}\""

    cmd += f" | samtools view -e {filter_cmd} -bS > {bam_out}".split() 
    cmd = " ".join(cmd)
    ps = subprocess.run(cmd, shell = True, capture_output = True)
    
    if ps.returncode != 0:
        print(" ".join(cmd))
        print(ps)
        sys.exit("Error occured running parasail_aligner")

def process_alignments(bam_fn, out_fn, min_align_score, max_align_score):
    """
    Post process alignments from parasail. Move alignment scores (AS tag) to MAPQ field, capping
    them at max_align_score. The set of alignments with the highest AS scores are retained for each read.
    """
    bam_fh = pysam.AlignmentFile(str(bam_fn), 'r') 

    out_fh = pysam.AlignmentFile(out_fn, "wb", header=bam_fh.header)
    
    # Iterate through alignments
    read_grp = ""
    best_as = 0
    read_cache = []
    
    for read in bam_fh.fetch(until_eof=True):

        try:
            ascore = read.get_tag("AS")
        except KeyError:
            continue
        
        if ascore > max_align_score:
            ascore = max_align_score
        read.mapping_quality = ascore

        if ascore < min_align_score:
            continue
    
        if read.query_name != read_grp:
            read_grp = read.query_name 
            for r in read_cache:
                out_fh.write(r)
            read_cache.clear()
            read_cache.append(read)
            best_as = read.mapping_quality
        
        else:
            if ascore > best_as:
                read_cache.clear()
                read_cache.append(read)
                best_as = ascore
            elif ascore == best_as:
                read_cache.append(read)
            else:
                continue 

    for r in read_cache:
        out_fh.write(r) 

    bam_fh.close()
    out_fh.close()


def align_reads(fastq_fn, fasta_ref, bam_out, threads = 1, min_align_score = 10, max_align_score = 255, mem_budget = "4GB"):
    """
    Use the parasail pairwise aligner to align fastq_fn to fasta_ref. 
    
    """
    bam_out = Path(bam_out)
    ps_bam_out = bam_out.with_suffix(".ps.bam")

    run_parasail(fastq_fn, ps_bam_out, fasta_ref, threads, mem_budget, min_align_score)
    
    as_bam_out = bam_out.with_suffix(".as.bam")
    process_alignments(ps_bam_out, as_bam_out, min_align_score, max_align_score)

    # add MD tag to bam, sort and index 
    calmd_bam_out = bam_out.with_suffix(".calmd.bam")
    
    cmd = f"samtools calmd -b {as_bam_out} {fasta_ref} > {calmd_bam_out}"
    ps = subprocess.run(cmd, shell = True, capture_output = True)
    
    if ps.returncode != 0:
        print(cmd)
        print(ps)
        sys.exit("Error occured running samtools calmd")
    
    pysam.sort("-o", str(bam_out), str(calmd_bam_out))
    pysam.index(str(bam_out))
    
    os.remove(calmd_bam_out)
    os.remove(as_bam_out)
    os.remove(ps_bam_out)


def generate_test_seqs(fastq_fn, out_fn):
    """
    Reverse the sequence (and qual score) to generate pseudorandom sequence for assessing alignment score
    """
    with pysam.FastxFile(fastq_fn, persist = True) as fh, gzwrite(out_fn) as fout:
        for rec in fh:
            rseq = rec.sequence[::-1]
            rec.sequence = rseq
            rqual = rec.quality[::-1]
            rec.quality = rqual
            fout.write(str(rec) + "\n")

def cleanup(*args):
    for f in args:
        try:
            os.remove(f)
        except FileNotFoundError:
            continue

def check_exe(exe):
    if shutil.which(exe) is None:
        sys.exit(f"Unable to find {exe} executable")

def check_deps():
    deps = ["parasail_aligner", "samtools", "gunzip"]
    for dep in deps:
        check_exe(dep)

def make_output_dir(filepath):
    dn = os.path.dirname(filepath)
    if dn:
        os.makedirs(dn, exist_ok=True)

def get_args():
    parser = argparse.ArgumentParser(description="""
                                     Perform optimized alignment with parasail.
                                     Implements the method described by: 
                                     Yu Sun, Michael Piechotta, Isabel Naarmann-de Vries, Christoph Dieterich, Ann E Ehrenhofer-Murray, 
                                     Detection of queuosine and queuosine precursors in tRNAs by direct RNA sequencing,
                                     Nucleic Acids Research, Volume 51, Issue 20, 10 November 2023, Pages 11197–11212,
                                     https://doi.org/10.1093/nar/gkad826
                                     https://github.com/dieterich-lab/QutRNA
                                     """)
    parser.add_argument('-p', '--precision', help = "precision cutoff for alignment score", default = 0.95, type = float)
    parser.add_argument('-s', '--align-score', help = "cutoff for alignment score, if set will override precision parameter", type = int)
    parser.add_argument('-o', '--output', required = True, help='output prefix (can include a directory)')
    parser.add_argument('-t', '--threads', help = "number of threads for alignment", default = 1, type = int)
    parser.add_argument('-m', '--memory', help = "Memory limit for parasail", default = "8GB")
    parser.add_argument('REF', help='FASTA file with reference sequences')
    parser.add_argument('FASTQ', help='input FASTQ file')

    args = parser.parse_args()

    p = args.precision
    if p < 0 or p > 1:
        sys.exit("-p must be between 0 and 1")
    
    if args.threads < 1:
        sys.exit("-t must be > 0")
    
    return args


def main():
    args = get_args()

    check_deps()
    make_output_dir(args.output)
    
    out_prefix = args.output

    print("running parasail on input reads")
    og_bam = out_prefix + ".og.bam"
    align_reads(args.FASTQ, args.REF, og_bam, mem_budget = args.memory, threads = args.threads)
    
    aln_score_cutoff = 0
    if args.align_score:
        aln_score_cutoff = args.align_score
    else:
        rfq_out = out_prefix + ".rev.fastq.gz"
        print("generating reversed input fastq")
        generate_test_seqs(args.FASTQ, rfq_out)
        
        print("running parasail on rand/rev sequences")
        rev_bam = out_prefix + ".rev.bam"
        align_reads(rfq_out, args.REF, rev_bam, mem_budget = args.memory, threads = args.threads)

        aln_score_cutoff = threshold_alignments(og_bam, rev_bam, out_prefix, p = args.precision)
        print(f"alignment score cutoff determined to be {aln_score_cutoff}") 
        cleanup(rev_bam, rev_bam + "bai", rfq_out)

    print(f"filtering final bam for alignment score of {aln_score_cutoff}")
    final_bam = out_prefix + ".bam"
    filter_bam(og_bam, final_bam, aln_score_cutoff)

    pysam.index(final_bam)

    cleanup(og_bam, og_bam + ".bai")

if __name__ == "__main__": main()
