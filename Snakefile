configfile: "config.yaml"

samples = ["GSE139659"]

rule all:
    input:
        expand("{sample}_trimmed.fq", sample=samples)

rule trim:
    input:
        "{sample}.fastq"
    output:
        "{sample}_trimmed.fq"
    shell:
        """
        trim_galore -q {config[trim][q]} --phred33 --length {config[trim][length]} {input}
        """

rule mapping:
    input:
        "{sample}_trimmed.fq"
    output:
        "{sample}_aligned.bam"
    shell:
        """
        bowtie-build reference.fasta index_reference.fasta
        bowtie -p 4 -S index_reference.fasta {input} | samtools sort -@ 4 > {output}
        """

rule counting:
    input:
        "{sample}_aligned.bam"
    output:
        "stats/{sample}.txt"
    shell:
        """
        featureCounts -a annotation.gtf -o {output} {input}
        """