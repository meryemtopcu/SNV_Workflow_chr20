
rule index_reference:
    input:
        fasta = REF
    output:
        BWA_INDEX
    log:
        "logs/index_reference.log"
    shell:
        """
        bwa index {input.fasta} 2> {log}
        """


rule bwa_map:
    input:
        fq  = RAW_DIR + "/{sample}.fastq.gz",
        ref = REF,
        idx = BWA_INDEX
    output:
        RESULTS_DIR + "/mapping/{sample}.bam"
    log:
        "logs/bwa_{sample}.log"
    threads: 4
    shell:
        """
        mkdir -p {RESULTS_DIR}/mapping
        bwa mem {input.ref} {input.fq} \
            | samtools view -bS - \
            > {output} 2> {log}
        """


rule sort_bam:
    input:
        RESULTS_DIR + "/mapping/{sample}.bam"
    output:
        RESULTS_DIR + "/mapping/{sample}.sorted.bam"
    log:
        "logs/sort_{sample}.log"
    threads: 4
    shell:
        """
        samtools sort -o {output} {input} 2> {log}
        """


rule index_bam:
    input:
        RESULTS_DIR + "/mapping/{sample}.sorted.bam"
    output:
        RESULTS_DIR + "/mapping/{sample}.sorted.bam.bai"
    log:
        "logs/index_{sample}.log"
    shell:
        """
        samtools index {input} 2> {log}
        """

