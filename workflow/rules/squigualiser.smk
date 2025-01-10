
rule blue_crab:
    input:
        rules.merge_pods.output,
    output:
        os.path.join(outdir, "squigualizer", "{sample}", "{sample}.blow5"),
    log:
        os.path.join(outdir, "logs", "blue_crab", "{sample}"),
    shell:
        """
        blue-crab p2s {input} -o {output}
        """


rule resquiggle:
    input:
        fastq=rules.ubam_to_fq.output,
        signal=rules.blue_crab.output,
    output:
        os.path.join(outdir, "squigualizer", "{sample}", "{sample}.tsv"),
    log:
        os.path.join(outdir, "logs", "resquiggle", "{sample}"),
    shell:
        """
        f5c_x86_64_linux_cuda \
            resquiggle \
            # XXX -B depends on memory available
            -B 10 \
            --rna \
            -c {input.fastq} \
            {input.signal} \
            -o {output} 
        """
