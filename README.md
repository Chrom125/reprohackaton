Pre-requis :

Installer apptainer :

wget https://github.com/apptainer/apptainer/releases/download/v1.4.3/apptainer_1.4.3_amd64.deb
sudo apt install -y ./apptainer_1.4.3_amd64.deb
rm ./apptainer_1.4.3_amd64.deb

Installer snakemake :

conda install -c bioconda snakemake

Commande pour lancer le snakefile :

snakemake -s Snakefile --use-singularity --singularity-args "--bind /tmp,/home/ubuntu/data" --cores 1