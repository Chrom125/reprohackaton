

# Description du projet

Ce projet vise à reproduire l'analyse bioinformatique de l'article : [Intracellular Staphylococcus aureus persisters upon antibiotic exposure](https://www.nature.com/articles/s41467-020-15966-7) dans le cadre de l'unité d'enseignement Hackaton reproductibilité du Master (M2) de biologie computationnelle : AMI2B de l'Université Paris-Saclay.

L'article porte sur l'identification des gènes différentiellement exprimés des bactéries de type ***Staphylococcus aureus*** lorsqu'elles deviennent persistantes à l'issue d'un contact avec des antibiotiques; la persistance leur conférant une tolérance à ceux-ci. Le but est ici de reprendre l'analyse d'expression différentielle qui compare les bactéries "Persister" et les bactéries "Control" en reproduisant deux MA plots, l'un à l'échelle de l'ensemble des gènes étudiés et l'autre restreint aux gènes de la traduction.  

L'analyse bioinformatique réalisée a été conçue dans une optique de reproductibilité computationnelle, en utilisant un système de workflow **Snakemake** pour l'orchestration des différentes tâches et des conteneurs de type **Apptainer** pour l'encapsulation de chaque outil bioinformatique et en assurer la portabilité. 

# Pipeline d'analyse

Le pipeline d'analyse automatisé couvre les étapes suivantes :

* **Trimming des reads** issus du séquençage RNA‑seq bactérien (Nettoyage des données brutes).
* **Préparation du génome de référence** : Téléchargement et indexation du génome de référence.
* **Mapping** : Alignement des reads nettoyés sur le génome de référence.
* **Quantification** : Génération de la table de comptage (*Count table*).
* **Analyse d'expression différentielle** (via le package R **DESeq2**) incluant la génération de graphiques de diagnostic :
    * MA‑plots
    * Volcano plots
* **Comparaison et vérification de la reproductibilité** des résultats par rapport à la table de comptage fournie par les auteurs :
    * Visualisation des intersections (Upset plot).
    * Comparaison des distributions de *log2FoldChange*.
    * Comparaison des moyennes de comptages normalisés (*baseMean*).
 

# Structure du projet

L'organisation des fichiers du projet est la suivante :

```text
.
├── Data/
│   ├── GeneSpecificInformation_NCTC8325.xlsx    # Informations d'annotation des gènes
│   └── metadata.tsv                             # Métadonnées et design expérimental : Correspondance entre les réplicats et leur label "Persister" ou "Control"
│
├── Singularity/Recipes/             # Recettes pour la construction des conteneurs
│   ├── bowtie-samtools.recipe                   # bowtie 0.12.7 et samtools 1.22.1
│   ├── cutadapt.recipe                          # cutadapt 1.11
│   ├── data_viz_tools_suite.def                 # packages R de visualtion graphique
│   ├── deseq.recipe                             # DESeq2 R package 1.16.1
│   ├── feature-counts.def                       # subreads 1.4.6
│   ├── keggrest.def                             # keggrest R package 1.50.0
│   └── sratoolkit-fasterq-dump.def              # sratoolkit current version
│
├── envs/
│   └── bowtie2.yaml                      # Environnement Conda pour bowtie2 2.5.4 et samtools 1.22.1
│
├── scripts/                               # Scripts R exécutés par le pipeline
│   ├── KEGG_functional_annotation.R              # Interrogation de l'API REST de KEGG pour avoir les annotations fonctionnelles
│   ├── comparison_with_authors_results.R         # Comparaison aux résultats obtenus de l'analyse différentielle à partir des comptages des auteurs (Upset plots, boxplots)
│   ├── differential_expression.R                 # Reproduction de l'analyse différentielle
│   ├── differential_expression_authors.R         # Analyse différentielle à partir de la table de comptage des auteurs
│   └── plots_differential_expression_.results.R  # Génération des MA plots, volcano plot
│
├── config.yaml       # Fichier de configuration pour le pipeline de reproduction des résultats de l'article (Trimming, Mapping, Indextation, analyse différentielle, MA plots)
├── config2.yaml      # Fichier de configuration pour le pipeline de comparaison des résultats à ceux obtenus avec la table de comptage des auteurs (Analyses R/Comparaisons)
├── Snakefile         # Workflow Snakemake de reproduction des résultats de l'article
├── Snakefile2        # Workflow Snakemake pour la comparaison des résultats obtenus à ceux obtenus à partir de la table de comptage des auteurs
└── README.md         # Documentation du projet
```

# Protocole d'exécution du workflow
## Cloner le projet
```
git clone https://github.com/Chrom125/reprohackaton
cd reprohackaton/
```

## Préalables:

### Apptainer :
```
wget https://github.com/apptainer/apptainer/releases/download/v1.4.3/apptainer_1.4.3_amd64.deb
sudo apt install -y ./apptainer_1.4.3_amd64.deb
rm ./apptainer_1.4.3_amd64.deb
```

### Snakemake :

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
