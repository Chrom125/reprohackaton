# Cloner le projet :
```
git clone https://github.com/Chrom125/reprohackaton
cd reprohackaton/
```

# Pre-requis :

## Apptainer :
```
wget https://github.com/apptainer/apptainer/releases/download/v1.4.3/apptainer_1.4.3_amd64.deb
sudo apt install -y ./apptainer_1.4.3_amd64.deb
rm ./apptainer_1.4.3_amd64.deb
```

## Snakemake :
```
conda init
conda install -c bioconda snakemake
conda activate
```

# Générer le dag graphe :
```
snakemake -S Snakefile --dag | dot -Tpng -o dag.png
```

# Lancer l'analyse :
new
```
snakemake -s Snakefile --use-apptainer --singularity-args --cores <number-of-cores>
```

old
```
snakemake -s Snakefile --use-singularity --singularity-args "--bind /tmp,/home/ubuntu/data" --cores 1
```
