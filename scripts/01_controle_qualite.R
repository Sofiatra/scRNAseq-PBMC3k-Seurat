# =============================================================================
# Étape 1 : chargement des données et contrôle qualité (QC) des cellules
# Entrée  : data/filtered_gene_bc_matrices/hg19/
# Sorties : results/01_*.png, results/pbmc_qc.rds
# Lancer depuis la racine du projet :  Rscript scripts/01_controle_qualite.R
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})
dir.create("results", showWarnings = FALSE)
set.seed(42)

# 1. Lecture de la matrice de comptage 10x (gènes x cellules) ----------------
counts <- Read10X(data.dir = "data/filtered_gene_bc_matrices/hg19/")
cat("Dimensions de la matrice brute (gènes x cellules) :", dim(counts), "\n")

# 2. Création de l'objet Seurat ----------------------------------------------
# min.cells = 3    : on garde les gènes détectés dans au moins 3 cellules
# min.features = 200 : on garde les cellules avec au moins 200 gènes détectés
pbmc <- CreateSeuratObject(counts = counts, project = "pbmc3k",
                           min.cells = 3, min.features = 200)
cat("Après création de l'objet :", ncol(pbmc), "cellules,", nrow(pbmc), "gènes\n")

# 3. Pourcentage de lectures mitochondriales -----------------------------------
# Une cellule abîmée ou mourante perd son ARN cytoplasmique mais garde ses
# mitochondries : un % mitochondrial élevé signale une cellule de mauvaise qualité.
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# 4. Visualisation des métriques AVANT filtrage --------------------------------
p_violin <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    ncol = 3, pt.size = 0.1)
ggsave("results/01_qc_violons_avant_filtrage.png", p_violin, width = 12, height = 5, dpi = 150)

p1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("results/01_qc_nuages_de_points.png", p1 + p2, width = 12, height = 5, dpi = 150)

# 5. Filtrage des cellules -----------------------------------------------------
# - moins de 200 gènes : gouttelette vide ou cellule abîmée
# - plus de 2 500 gènes : probable doublet (2 cellules dans une gouttelette)
# - plus de 5 % de mitochondrial : cellule abîmée
seuil_min_genes <- 200
seuil_max_genes <- 2500
seuil_max_mt    <- 5
n_avant <- ncol(pbmc)
pbmc <- subset(pbmc, subset = nFeature_RNA > seuil_min_genes &
                              nFeature_RNA < seuil_max_genes &
                              percent.mt   < seuil_max_mt)
cat("Cellules avant filtrage :", n_avant, "| après filtrage :", ncol(pbmc),
    "| retirées :", n_avant - ncol(pbmc), "\n")

p_violin2 <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                     ncol = 3, pt.size = 0.1)
ggsave("results/01_qc_violons_apres_filtrage.png", p_violin2, width = 12, height = 5, dpi = 150)

# 6. Sauvegarde ------------------------------------------------------------
saveRDS(pbmc, "results/pbmc_qc.rds")
cat("Objet sauvegardé : results/pbmc_qc.rds\n")
