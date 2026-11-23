# =============================================================================
# Étape 2 : normalisation, gènes variables, ACP, clustering et UMAP
# Entrée  : results/pbmc_qc.rds
# Sorties : results/02_*.png, results/pbmc_clusters.rds
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})
set.seed(42)
pbmc <- readRDS("results/pbmc_qc.rds")

# 1. Normalisation : chaque cellule est ramenée à 10 000 lectures puis log1p --
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 1e4)

# 2. Gènes les plus variables (ceux qui portent l'information biologique) ----
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)
cat("Les 10 gènes les plus variables :", paste(top10, collapse = ", "), "\n")
p_hvg <- LabelPoints(plot = VariableFeaturePlot(pbmc), points = top10, repel = TRUE)
ggsave("results/02_genes_variables.png", p_hvg, width = 8, height = 6, dpi = 150)

# 3. Mise à l'échelle (moyenne 0, variance 1 pour chaque gène) ----------------
pbmc <- ScaleData(pbmc, features = rownames(pbmc))

# 4. ACP sur les gènes variables -----------------------------------------------
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), verbose = FALSE)
print(pbmc[["pca"]], dims = 1:3, nfeatures = 5)
ggsave("results/02_acp_pc1_pc2.png", DimPlot(pbmc, reduction = "pca"), width = 7, height = 6, dpi = 150)
ggsave("results/02_acp_heatmap_pc1_6.png",
       DimHeatmap(pbmc, dims = 1:6, cells = 500, balanced = TRUE, fast = FALSE),
       width = 10, height = 10, dpi = 150)

# 5. Combien de composantes garder ? Le "coude" de ce graphique -------------
ggsave("results/02_elbow_plot.png", ElbowPlot(pbmc, ndims = 30), width = 7, height = 5, dpi = 150)
n_pcs <- 10   # valeur classique pour PBMC 3k ; à justifier avec l'elbow plot

# 6. Graphe des plus proches voisins puis clustering (algorithme de Louvain) --
pbmc <- FindNeighbors(pbmc, dims = 1:n_pcs)
pbmc <- FindClusters(pbmc, resolution = 0.5)
cat("Nombre de clusters :", length(levels(Idents(pbmc))), "\n")
print(table(Idents(pbmc)))

# 7. UMAP : projection en 2D pour la visualisation ----------------------------
pbmc <- RunUMAP(pbmc, dims = 1:n_pcs)
p_umap <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("PBMC 3k : clusters (résolution 0.5)")
ggsave("results/02_umap_clusters.png", p_umap, width = 8, height = 6, dpi = 150)

# 8. Effet de la résolution : comparer 0.2, 0.5 et 1.0 -----------------------
plots <- list()
for (res in c(0.2, 0.5, 1.0)) {
  tmp <- FindClusters(pbmc, resolution = res, verbose = FALSE)
  plots[[as.character(res)]] <- DimPlot(tmp, reduction = "umap", label = TRUE) +
    ggtitle(paste("résolution", res)) + NoLegend()
}
ggsave("results/02_umap_effet_resolution.png", wrap_plots(plots, ncol = 3),
       width = 15, height = 5, dpi = 150)

# On revient à la résolution 0.5 pour la suite
pbmc <- FindClusters(pbmc, resolution = 0.5, verbose = FALSE)
saveRDS(pbmc, "results/pbmc_clusters.rds")
cat("Objet sauvegardé : results/pbmc_clusters.rds\n")
