# =============================================================================
# Étape 4 : annotation automatique avec SingleR et comparaison avec l'étape 3
# Entrée  : results/pbmc_annote.rds
# Sorties : results/04_*.png, results/04_*.csv, results/pbmc_final.rds
# Nécessite une connexion internet au premier lancement (téléchargement de la
# référence celldex, ~100 Mo, mise en cache ensuite).
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(SingleCellExperiment)
  library(ggplot2)
  library(pheatmap)
})
set.seed(42)
pbmc <- readRDS("results/pbmc_annote.rds")

# 1. Référence : Human Primary Cell Atlas (profils de types cellulaires purs) -
ref <- celldex::HumanPrimaryCellAtlasData()

# 2. Conversion de l'objet Seurat en SingleCellExperiment ----------------------
sce <- as.SingleCellExperiment(pbmc)

# 3. Prédiction cellule par cellule ------------------------------------------
pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)
pbmc$SingleR <- pred$pruned.labels
cat("Types prédits par SingleR :\n")
print(sort(table(pbmc$SingleR), decreasing = TRUE))

png("results/04_singler_scores_heatmap.png", width = 1400, height = 900, res = 120)
plotScoreHeatmap(pred)
dev.off()

p <- DimPlot(pbmc, reduction = "umap", group.by = "SingleR", label = TRUE, repel = TRUE) +
  ggtitle("Annotation automatique (SingleR, Human Primary Cell Atlas)")
ggsave("results/04_umap_singler.png", p, width = 10, height = 7, dpi = 150)

# 4. Comparaison annotation manuelle vs automatique --------------------------
confusion <- table(manuelle = pbmc$type_cellulaire, SingleR = pbmc$SingleR,
                   useNA = "ifany")   # NA = cellules que SingleR n'a pas su annoter
write.csv(as.data.frame.matrix(confusion), "results/04_comparaison_manuelle_vs_singler.csv")
png("results/04_comparaison_manuelle_vs_singler.png", width = 1400, height = 900, res = 120)
pheatmap(log10(as.matrix(confusion) + 1), display_numbers = as.matrix(confusion),
         cluster_rows = FALSE, cluster_cols = FALSE,
         main = "Annotation manuelle (lignes) vs SingleR (colonnes)")
dev.off()

saveRDS(pbmc, "results/pbmc_final.rds")

# 5. Traçabilité : versions des logiciels utilisés -----------------------------
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
cat("Terminé. Versions des packages : results/sessionInfo.txt\n")
