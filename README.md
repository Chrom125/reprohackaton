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
exec bash -l
```

# Générer le dag graphe :
```
sudo apt install graphviz
snakemake -s Snakefile --dag | dot -Tpng -o dag.png
```

# Lancer l'analyse :
new
```
snakemake -s Snakefile --use-apptainer --cores <number-of-cores>
```

new2
'''
snakemake -s Snakefile --use-singularity --singularity-args "--bind $(pwd)" --cores 1
'''