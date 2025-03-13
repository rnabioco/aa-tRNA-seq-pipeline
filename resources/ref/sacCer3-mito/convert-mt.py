#! /usr/bin/env python

from Bio import SeqIO
import sys
import pdb

adaptor5 = "CCTAAGAGCAAGAAGAAGCCTGGN"
adaptor3 = "CCAACCTTGCCTTAAAAAAAAAA"

for record in SeqIO.parse(sys.argv[1], "fasta"):
    if "pseudo" in record.description: continue

    fs = record.description.split(" ")
    name = f"mt-{fs[3]}-{fs[4].replace('(','').replace(')','')}"
    seq = adaptor5 + str(record.seq) + adaptor3

    print(f">{name}\n{seq}")
