configfile: "config.yaml"

import pandas as pd

metadata = pd.read_table(config['sample_table'], dtype=str)
sra_id = sample_table['SRA_ID'].unique().tolist()
label_to_sra_id = metadata.set_index("Label")["SRA_ID"].to_dict()

#Final output: a single counts file with all samples combined
rule all:
    input:
        counts_file = "counts.txt"


rule download_data:
    output:
        "{config[output_dir]}/{config[download_data][output_dir]}/{sample}.fastq"
    container:
        "https://zenodo.org/records/17423176/files/sratoolkit-fasterq-dump.sif"
    params: 
        read_data = lambda wildcards.sample : label_to_sra_id[wildcards.sample]
    shell:
        """
        fasterq-dump --threads {config[download_data][threads]} --progress -O {config[output_dir]}/{config[download_data][output_dir]} \
        -o {sample}.fastq {params.read_data}
        """

rule trimming:
    input:
        "{config[output_dir]}/{config[download_data][output_dir]}/{sample}.fastq"
    output:
        "{config[output_dir]}/trimming/{sample}_trimmed.fastq"
    container:
        "https://zenodo.org/records/17424137/files/cutadapt.img?download=1"
    shell:
        """
        cutadapt -a {config[trimming][a]} -m {config[trimming][m]} -q {config[trimming][q]} -o {output} {input}
        """


rule reference_genome:
    output:
        "{config[output_dir]}/Reference_Genome/{config[reference_genome][filename]}"
    shell:
        """
        wget -q -O {output} {config[reference_genome][fasta_url]}
        """

rule genome_annotation:
    output:
        "{config[output_dir]}/Genome_Annotation/{config[reference_genome][annotation_filename]}"
    shell:
        """
        wget -q -O {output} {config[reference_genome][annotation_url]}
        """

rule genome_index:
    input:
        "{config[output_dir]}/Reference_Genome/{config[reference_genome][filename]}"
    output:
        expand("{config[output_dir]}/Reference_Genome/{config[genome_index][index_files_prefix]}.{i}.ebwt", i=[1,2,3,4,"rev.1","rev.2"])
    container:
        "https://zenodo.org/records/17425965/files/bowtie-samtools.img?download=1"
    shell:
        """
        bowtie-build {input} results/Reference_Genome/{config[genome_index][index_files_prefix]}
        """

rule mapping:
    input:
        trimmed = "{config[output_dir]}/trimming/{sample}_trimmed.fastq"
        index_reference = expand("{config[output_dir]}/Reference_Genome/{config[genome_index][index_files_prefix]}.{i}.ebwt", i=[1,2,3,4,"rev.1","rev.2"])
    output:
        "{config[output_dir]}/mapping/{sample}_aligned.bam"
    container:
        "https://zenodo.org/records/17425965/files/bowtie-samtools.img?download=1"
    shell:
        """
        bowtie -p {config[mapping][threads]} -S {input.index_reference} {input.trimmed} | samtools sort -@ {config[mapping][threads]} > {output}
        """

rule featurecounts:
    input:
        expand("{config[output_dir]}/mapping/{sample}_aligned.bam", sample=metadata['Label'].tolist())
    output:
        "{config[output_dir]}/featurecounts/counts.txt"
    container:
        "https://zenodo.org/records/17425965/files/bowtie-samtools.img?download=1"
    shell:
        """
        featureCounts -t {config[featurecounts][feature_type]} \
        -g {config[featurecounts][attribute_type]} -F {config[featurecounts][annotation_format]} \
        -T  {config[featurecounts][threads]} \
        -a {config[featurecounts][annotation_file]} -o {output} {input}
        """
