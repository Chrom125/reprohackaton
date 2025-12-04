# Cloner le projet :
```
git clone https://github.com/Chrom125/reprohackaton
cd reprohackaton/
```

# À propos du projet :

Ce projet vise à reproduire l'analyse bioinformatique de l'article : [Intracellular Staphylococcus aureus persisters upon antibiotic exposure](https://www.nature.com/articles/s41467-020-15966-7). 

L'article montre une différence d'expression des gènes de Staphylococcus aureus lorsque les bactéries sont dans un état de persistance, qui leur permet de résister aux antibiotiques. 

L'analyse bioinformatique part des fichiers fasta/fastq obtenus du séquençage RNA‑seq des bactéries dans un état de persistance et dans un état de contrôle. Elle vise à obtenir une liste de gènes différemment exprimés dont certains sont réprimés et d'autres surexprimés, ainsi que plusieurs figures (MA‑plots, PCA, heatmaps, et les "upseq plots" fournis dans le Snakefile2). Le pipeline implémente les étapes classiques : contrôle qualité, trimming, alignement, comptage par gène puis analyse différentielle (ex. DESeq2) — les détails de chaque étape sont définis dans les Snakefiles.

# Pre-requis :

## Apptainer :
```
wget https://github.com/apptainer/apptainer/releases/download/v1.4.3/apptainer_1.4.3_amd64.deb
sudo apt install -y ./apptainer_1.4.3_amd64.deb
rm ./apptainer_1.4.3_amd64.deb
```

## Snakemake :

```
sudo apt install snakemake
```

ou

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

Attention, les fichiers intermédiaires ne sont pas supprimés. Il faut prévoir assez de place pour les stocker.

Le fichier snakefile2 permet de comparer les résultats avec les auteurs. Il permet notamment de tracer les upseq plots.

```
snakemake -s Snakefile --use-singularity --singularity-args "--bind $(pwd)" --cores <number_of_cores>
```

# Description succincte des fichiers et dossiers

- Snakefile  
  Le pipeline principal Snakemake décrivant les étapes (QC, trimming, alignement, comptage, analyse différentielle, génération de figures). C'est le Snakefile à utiliser pour l'exécution standard.

- Snakefile2  
  Variante du pipeline (utilisée pour comparer avec les auteurs ou tracer des figures supplémentaires comme les "upseq plots").

- config.yaml  
  Fichier de configuration contenant les paramètres.

- data/  
  Contient des informations sur les gènes pour l'analyse.

- results/  
  Contient les sorties du pipeline : matrices de comptage, tables DESeq2, figures (MA‑plots, PCA, heatmaps, upseq plots), et listes de gènes différentiellement exprimés.

- scripts/  
  Scripts R et Python appelés par le Snakefile pour l'analyse différentielle et la génération de figures (p.ex. scripts/DESeq2.R, scripts/plot_* .R). Conserver ces scripts pour tracer la provenance des analyses.

- README.md  
  Ce fichier — informations d'usage et contexte du projet.
