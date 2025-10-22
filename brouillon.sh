snakemake -s Snakefile --cores 4

conda install -c bioconda snakemake

sudo apt install sra-toolkit # Most recent ok

fasterq-dump --threads 4 --progress GSE139659
gzip *.fastq

conda install -c bioconda fastqc
conda install -c bioconda cutadapt # version 1.11
conda install -c bioconda trim-galore

trim_galore -q 20 --phred33 --length 25 GSE139659.fastq

wget -q -O reference.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta"

wget -O reference.gff "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=CP000253.1"

conda install -c bioconda bowtie # Version 0.12.7

bowtie-build reference.fasta index_reference.fasta

bowtie -p 4 -S index_reference.fasta GSE139659.fastq.gz | samtools sort -@ 4 > GSE139659.bam

conda install -c bioconda subread

# FeatureCounts 1.4.6-p3

featureCounts --extraAttributes Name -t gene -g ID -F GTF -T <#CPUS> -a <GFF> -o counts.txt <BAM FILES>

# DESeq2 version 1.16

samtools index GSE139659.bam

snakemake -s Snakefile --cores 1 --forceall