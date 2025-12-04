

# Description du projet

Ce projet vise à reproduire l'analyse bioinformatique de l'article : [Intracellular Staphylococcus aureus persisters upon antibiotic exposure](https://www.nature.com/articles/s41467-020-15966-7) dans le cadre de l'unité d'enseignement Hackaton reproductibilité du Master (M2) de biologie computationnelle : AMI2B de l'Université Paris-Saclay.

L'article porte sur l'identification des gènes différentiellement exprimés des bactéries de type ***Staphylococcus aureus*** lorsqu'elles deviennent persistantes à l'issue d'un contact avec des antibiotiques; la persistance leur conférant une tolérance à ceux-ci. Le but est ici de reprendre l'analyse d'expression différentielle qui compare les bactéries "Persister" et les bactéries "Control" en reproduisant deux MA plots, l'un à l'échelle de l'ensemble des gènes étudiés et l'autre restreint aux gènes de la traduction.  

L'analyse bioinformatique réalisée a été conçue dans une optique de reproductibilité computationnelle, en utilisant un système de workflow **Snakemake** pour l'orchestration des différentes tâches et des conteneurs de type **Apptainer** pour l'encapsulation de chaque outil bioinformatique et en assurer la portabilité. 

Le pipeline d'analyse inclut: le trimming des reads issus du séquençage RNA‑seq des bactéries, le téléchargement et l'indexation du génome de référence, le mapping des reads sur le génome de référence, la génération d'une table de comptage, l'analyse d'expression différentielle via le package R DESeq2 avec construction de MA‑plots et volcano plot ainsi que la comparaison des résultats obtenus avec ceux issus de la table de comptage des auteurs (Upset plot et comparaison des distributions de log2FoldChange et de moyenne de comptages normalisés).

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
