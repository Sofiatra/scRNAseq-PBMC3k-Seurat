# =============================================================================
# Étape 3 : gènes marqueurs et annotation des types cellulaires
# Entrée  : results/pbmc_clusters.rds
# Sorties : results/03_*.png, results/03_*.csv, results/pbmc_annote.rds
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})
set.seed(42)
pbmc <- readRDS("results/pbmc_clusters.rds")

# 1. Gènes marqueurs de chaque cluster (test de Wilcoxon, gènes surexprimés) -
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "results/03_marqueurs_tous_clusters.csv", row.names = FALSE)

top5 <- markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 5)
write.csv(top5, "results/03_marqueurs_top5_par_cluster.csv", row.names = FALSE)
print(top5, n = 60)

top10 <- markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 10)
ggsave("results/03_heatmap_top10_marqueurs.png",
       DoHeatmap(pbmc, features = top10$gene) + NoLegend(),
       width = 12, height = 12, dpi = 150)

# 2. Gènes marqueurs connus des PBMC (littérature) ----------------------------
signatures <- list(
  "T CD4 naïfs"        = c("IL7R", "CCR7"),
  "T CD4 mémoire"      = c("IL7R", "S100A4"),
  "Monocytes CD14+"    = c("CD14", "LYZ"),
  "Lymphocytes B"      = c("MS4A1", "CD79A"),
  "T CD8"              = c("CD8A", "CD8B"),
  "Monocytes FCGR3A+"  = c("FCGR3A", "MS4A7"),
  "Cellules NK"        = c("GNLY", "NKG7"),
  "Cellules dendritiques" = c("FCER1A", "CST3"),
  "Plaquettes"         = c("PPBP", "PF4")
)
genes_connus <- unique(unlist(signatures))
ggsave("results/03_featureplot_marqueurs_connus.png",
       FeaturePlot(pbmc, features = genes_connus, ncol = 4),
       width = 16, height = 16, dpi = 120)
ggsave("results/03_dotplot_marqueurs_connus.png",
       DotPlot(pbmc, features = genes_connus) + RotatedAxis(),
       width = 12, height = 6, dpi = 150)

# 3. Annotation semi-automatique : score de chaque signature par cluster ----
# Pour chaque cluster, on calcule l'expression moyenne de chaque signature et
# on attribue la signature la plus forte. Cette méthode ne dépend pas de la
# numérotation des clusters (qui peut changer d'une version de Seurat à l'autre).
moy <- AverageExpression(pbmc, features = genes_connus, assays = "RNA",
                         layer = "data")$RNA
# AverageExpression renvoie des moyennes à l'échelle linéaire (exp des données
# log) : sans correction, les gènes très exprimés (LYZ, GNLY, NKG7...) écrasent
# les autres. On repasse en log puis on centre-réduit chaque gène entre clusters
# (z-score), pour que chaque marqueur pèse autant dans le score.
moy <- log1p(as.matrix(moy))
moy_z <- t(scale(t(moy)))
scores <- sapply(signatures, function(g) colMeans(moy_z[g, , drop = FALSE]))
# scores : lignes = clusters, colonnes = types cellulaires
annotation_auto <- colnames(scores)[apply(scores, 1, which.max)]

# Correction d'analyste : T CD4 naïfs vs mémoire --------------------------------
# Les deux signatures partagent IL7R, et S100A4 est surtout exprimé par les
# cellules myéloïdes et NK : le score seul classe les deux clusters CD4 en
# "naïfs". On les départage sur CCR7 (naïfs) moins S100A4 (mémoire) : le cluster
# CD4 le plus "CCR7 haut / S100A4 bas" est naïf, les autres sont mémoire.
# Vérifié avec les marqueurs du cluster : CCR7, LEF1, MAL (naïfs) contre
# AQP3, CD40LG, TRAT1 (mémoire).
annotation <- annotation_auto
cd4 <- which(annotation %in% c("T CD4 naïfs", "T CD4 mémoire"))
if (length(cd4) >= 2) {
  naivete <- moy_z["CCR7", cd4] - moy_z["S100A4", cd4]
  annotation[cd4] <- "T CD4 mémoire"
  annotation[cd4[which.max(naivete)]] <- "T CD4 naïfs"
}

table_annot <- data.frame(cluster = rownames(scores),
                          type_score = annotation_auto,
                          type_final = annotation,
                          round(scores, 2), check.names = FALSE)
print(table_annot)
write.csv(table_annot, "results/03_annotation_par_cluster.csv", row.names = FALSE)

# Les noms renvoyés par AverageExpression peuvent être "g0", "g1"... selon la
# version de Seurat : on retrouve l'identifiant du cluster en retirant le "g".
noms_clusters <- sub("^g", "", rownames(scores))
correspondance <- setNames(annotation, noms_clusters)
pbmc$type_cellulaire <- unname(correspondance[as.character(Idents(pbmc))])

p_annot <- DimPlot(pbmc, reduction = "umap", group.by = "type_cellulaire",
                   label = TRUE, repel = TRUE, pt.size = 0.5) +
  ggtitle("PBMC 3k : annotation par marqueurs connus") + NoLegend()
ggsave("results/03_umap_annote.png", p_annot, width = 9, height = 7, dpi = 150)

# 4. Proportions de chaque type cellulaire -----------------------------------
prop <- as.data.frame(table(pbmc$type_cellulaire))
colnames(prop) <- c("type_cellulaire", "n_cellules")
prop$pourcentage <- round(100 * prop$n_cellules / sum(prop$n_cellules), 1)
print(prop)
write.csv(prop, "results/03_proportions_types_cellulaires.csv", row.names = FALSE)
p_prop <- ggplot(prop, aes(x = reorder(type_cellulaire, -n_cellules), y = pourcentage)) +
  geom_col(fill = "#1F3B57") + labs(x = NULL, y = "% des cellules") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("results/03_proportions_types_cellulaires.png", p_prop, width = 8, height = 5, dpi = 150)

saveRDS(pbmc, "results/pbmc_annote.rds")
cat("Objet sauvegardé : results/pbmc_annote.rds\n")
