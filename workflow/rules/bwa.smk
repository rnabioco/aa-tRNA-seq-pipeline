
rule bwa_idx:
  input:
    config["fasta"]
  output:
    multiext(config["fasta"], ".amb", ".ann", ".bwt", ".pac", ".sa")
  log:
    os.path.join(outdir, "logs", "bwa_idx", "log") 
  shell:
    """
    bwa index {input}
    """ 
    
rule bwa:
  """
  align reads to tRNA references with bwa mem
  """
  input:
    reads = rules.ubam_to_fq.output,
    idx = rules.bwa_idx.output
  output:
    bam = os.path.join(outdir, "bams", "{sample}", "{sample}.bam"),
    bai = os.path.join(outdir, "bams", "{sample}", "{sample}.bam.bai"),
  params:
    index = config["fasta"],
    src = config["src"]
  log:
    os.path.join(outdir, "logs", "bwa", "{sample}") 
  threads: 4
  shell:
    """
    bwa mem -t {threads} -W 13 -k 6 -T 20 -x ont2d {params.index} {input.reads} \
        | samtools view -F4 -hu - \
        | python {params.src}/filter_reads.py "-" "-" \
        | samtools sort -o {output.bam}

    samtools index {output.bam}
    """


