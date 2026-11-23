#!/usr/bin/env bash
# Télécharge le jeu de données public PBMC 3k (10x Genomics) :
# ~2 700 cellules mononucléées du sang périphérique (PBMC) d'un donneur sain.
# À lancer depuis la racine du projet :  bash scripts/00_telecharger_donnees.sh
set -euo pipefail

mkdir -p data
cd data

URL="https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"

if [ ! -d "filtered_gene_bc_matrices/hg19" ]; then
  echo ">> Téléchargement de PBMC 3k (~7 Mo)..."
  curl -L --fail -o pbmc3k.tar.gz "$URL"
  echo ">> Décompression..."
  tar -xzf pbmc3k.tar.gz
  rm pbmc3k.tar.gz
else
  echo ">> Les données sont déjà présentes, rien à faire."
fi

echo ">> Contenu attendu : barcodes.tsv, genes.tsv, matrix.mtx"
ls -lh filtered_gene_bc_matrices/hg19/
