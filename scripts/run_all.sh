#!/usr/bin/env bash
# Lance tout le projet d'un coup, depuis la racine du projet :
#   conda activate singlecell_pbmc
#   bash scripts/run_all.sh
set -euo pipefail
bash    scripts/00_telecharger_donnees.sh
Rscript scripts/01_controle_qualite.R
Rscript scripts/02_normalisation_clustering.R
Rscript scripts/03_marqueurs_annotation.R
Rscript scripts/04_annotation_automatique_singler.R
echo ">> Projet terminé : toutes les figures et tables sont dans results/"
