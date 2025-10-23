configfile: "config.yaml"

samples = ["GSE139659"]

rule all:
    input:
        expand("results/featurecounts/{sample}_counts.txt", sample=samples)

rule download_data:
    output:
        "results/raw-data/{sample}.fastq"
    container:
        "https://zenodo.org/records/17423176/files/sratoolkit-fasterq-dump.sif"
    shell:
        """
        fasterq-dump --threads 4 --progress -O results/raw-data {wildcards.sample}
        """

rule trimming:
    input:
        "results/raw-data/{sample}.fastq"
    output:
        "results/trimming/{sample}_trimmed.fastq"
    container:
        "cutadapt.img"
    shell:
        """
        cutadapt -a {config[trimming][a]} -m {config[trimming][m]} -o {output} {input}
        """


rule reference_genome:
    output:
        "results/Reference_Genome/reference.fasta"
    log:
        "logs/reference_genome.log"
    container:
        ""
    shell:
        """
        wget -q -O {output} "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta" &> {log}
        """

rule genome_annotation:
    output:
        "results/Genome_Annotation/reference.gff"
    log:
        "logs/genome_annotation.log"
    container:
        ""
    shell:
        """
        wget -q -O {output} "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=CP000253.1" &> {log}
        """

rule genome_index:
    input:
        "results/Reference_Genome/reference.fasta"
    output:
        "results/Reference_Genome/index_reference.1.ebwt"
    log:
        "logs/genome_index.log"
    threads: config["threads"]
    container:
        ""
    shell:
        """
        bowtie-build {input} results/Reference_Genome/index_reference &> {log}
        """



rule mapping:
    input:
        "results/trimming/{sample}_trimmed.fastq"
    output:
        "results/mapping/{sample}_aligned.bam"
    shell:
        """
        bowtie-build reference.fasta index_reference.fasta
        bowtie -p 4 -S index_reference.fasta {input} | samtools sort -@ 4 > {output}
        """

rule featurecounts:
    input:
        "results/mapping/{sample}_aligned.bam"
    output:
        "results/featurecounts/{sample}_counts.txt"
    shell:
        """
        featureCounts --extraAttributes Name -t gene -g ID -F GTF -T 4 -a reference.gff -o {output} {input}
        """
