rule call_variants:
    input:
        bam = RESULTS_DIR + "/mapping/{sample}.sorted.bam",
        bai = RESULTS_DIR + "/mapping/{sample}.sorted.bam.bai",
        ref = REF
    output:
        bcf = RESULTS_DIR + "/variants/{sample}.bcf",
        vcf = RESULTS_DIR + "/variants/{sample}.vcf.gz"
    log:
        "logs/calling_{sample}.log"
    threads: 4
    shell:
        """
        mkdir -p {RESULTS_DIR}/variants

        bcftools mpileup -f {input.ref} {input.bam} \
            | bcftools call -mv -Ob -o {output.bcf} 2> {log}

        bcftools view {output.bcf} -Oz -o {output.vcf}
        bcftools index {output.vcf}
        """

