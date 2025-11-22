rule fastqc:
    input:
        RAW_DIR + "/{sample}.fastq.gz"
    output:
        html = RESULTS_DIR + "/qc/{sample}_fastqc.html",
        zip  = RESULTS_DIR + "/qc/{sample}_fastqc.zip"
    log:
        "logs/fastqc_{sample}.log"
    shell:
        """
        mkdir -p {RESULTS_DIR}/qc
        fastqc {input} --outdir {RESULTS_DIR}/qc &> {log}
        """

