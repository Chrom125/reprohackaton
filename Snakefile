configfile: "config.yaml"

############################################## Script for generating a mapping of sample to SRA IDs ##############################################

sample_to_sra_id = {}
with open(config["sample_table"]) as f:
    next(f)
    for line in f:
        sra_id, sample, label = line.strip().split("\t")
        sample_to_sra_id[sample] = sra_id



########################################################## Snakefile rules ############################################
#Final output: a single counts file with all samples combined
rule all:
    input:
        "results/featurecounts/processed_counts.txt"

rule download_data:
    output:
        "results/raw-data/{sample}.fastq"
    log:
        out = "logs/download_data/{sample}.log"
    container:
        "https://zenodo.org/records/17423176/files/sratoolkit-fasterq-dump.sif"
    params: 
        read_data = lambda wildcards : sample_to_sra_id[wildcards.sample]
    shell:
        """
        fasterq-dump --threads {config[download_data][threads]} --progress -O results/raw-data \
        -o {wildcards.sample}.fastq {params.read_data} 2>{log.out}
        """

rule trimming:
    input:
        "results/raw-data/{sample}.fastq"
    output:
        "results/trimming/{sample}_trimmed.fastq"
    log:
        out = "logs/trimming/{sample}.out",
        err = "logs/trimming/{sample}.err"
    container:
        "https://zenodo.org/records/17424137/files/cutadapt.img?download=1"
    shell:
        """
        cutadapt -a {config[trimming][a]} -m {config[trimming][m]} -q {config[trimming][q]} -o {output} {input}\
        >{log.out} 2>{log.err}
        """


rule reference_genome:
    output:
        "results/Reference_Genome/reference.fasta"
    log:
        out = "logs/Reference_Genome/reference_genome_download.log"
    shell:
        """
        wget -q -O {output} "{config[reference_genome][fasta_url]}" 2>{log.out}
        """

rule genome_annotation:
    output:
        "results/Genome_Annotation/reference.gff"
    log:
        out = "logs/Genome_Annotation/reference_genome_annotation_download.log"
    shell:
        """
        wget -q -O {output} "{config[genome_annotation][annotation_url]}" 2>{log.out}
        """

rule genome_index:
    input:
        "results/Reference_Genome/reference.fasta"
    output:
        expand("results/Reference_Genome/index_reference.{i}.ebwt", i=[1,2,3,4,"rev.1","rev.2"])
    log:
        out = "logs/Reference_Genome/reference_genome_index.out",
        err = "logs/Reference_Genome/reference_genome_index.err"
    container:
        "https://zenodo.org/records/17425965/files/bowtie-samtools.img?download=1"
    shell:
        """
        bowtie-build {input} results/Reference_Genome/index_reference >{log.out} 2>{log.err}
        """

rule mapping:
    input:
        trimmed = "results/trimming/{sample}_trimmed.fastq",
        index_reference = expand("results/Reference_Genome/index_reference.{i}.ebwt", i=[1,2,3,4,"rev.1","rev.2"])
    output:
        "results/mapping/{sample}_aligned.bam"
    log:
        mapping = "logs/mapping/bowtie/{sample}_trimmed_mapping.log",
        sorting = "logs/mapping/samtools-sort/{sample}_aligned_sorting.log"
    container:
        "https://zenodo.org/records/17425965/files/bowtie-samtools.img?download=1"
    shell:
        """
        bowtie -p {config[mapping][threads]} -S results/Reference_Genome/index_reference \
        {input.trimmed} 2>{log.mapping} | samtools sort -@ {config[mapping][threads]} >{output} 2>{log.sorting}
        """

rule featurecounts:
    input:
        mapping = expand("results/mapping/{sample}_aligned.bam", sample= list(sample_to_sra_id.keys())),
        annotation = "results/Genome_Annotation/reference.gff"
    output:
        "results/featurecounts/counts.txt"
    log:
        out = "logs/featurecounts/counts.log"
    container:
        "https://zenodo.org/records/17426103/files/feature-counts.sif?download=1"
    shell:
        """
        featureCounts -t {config[featurecounts][feature_type]} \
        -g {config[featurecounts][attribute_type]} -F {config[featurecounts][annotation_format]} \
        -T  {config[featurecounts][threads]} \
        -a {input.annotation} \
        -s {config[featurecounts][s]} -o {output} {input.mapping} 2>{log.out}
        """

rule processing_counts:
    input : 
        counts= "results/featurecounts/counts.txt"
    output:
        processed_counts="results/featurecounts/processed_counts.txt"
    run:
        with open(input.counts) as infile, open(output.processed_counts, "w") as outfile:
            for line in infile:
                if line.startswith("#"):
                    continue             #ignoring the header lines starting with #
                line = line.strip().split("\t")
                if line[0] == "Geneid":       # handling the table header
                    for i in range(6, len(line)):
                        a = line[i].split("/")
                        a = a[-1].replace("_aligned.bam", "")#Extracting only the sample label from the full path
                        line[i] = a
                if line[0].startswith("gene-"):
                    line[0] = line[0].replace("gene-", "")  #Removing 'gene-' prefix from Gene IDs
                
                #Keeping only Geneid and counts columns
                line = [line[0]] + line[6:]   #Keeping only Geneid and counts columns
                outfile.write("\t".join(line) + "\n")