#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(anndataR)
  library(SingleCellExperiment)
  library(slingshot)
  library(ggplot2)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
h5ad_file <- args[1]
output_prefix <- args[2]
n_pc_use <- 10L

sce <- read_h5ad(h5ad_file, as = "SingleCellExperiment")

if (inherits(reducedDim(sce, "X_pca"), "LinearEmbeddingMatrix")) {
  reducedDim(sce, "X_pca") <- as.matrix(reducedDim(sce, "X_pca"))
}

n_cells <- ncol(sce)
k <- max(1L, floor(0.2 * n_cells))
set.seed(42)
idx <- sample.int(n_cells, k, replace = FALSE)
sce_fit <- sce[, idx]

ncp <- min(n_pc_use, ncol(as.matrix(reducedDim(sce_fit, "X_pca"))))
rd <- as.matrix(reducedDim(sce_fit, "X_pca"))[, seq_len(ncp), drop = FALSE]

lineages <- getLineages(
  data = rd,
  clusterLabels = colData(sce_fit)$cell_label,
  start.clus = "B IN",
  end.clus = "B Mem"
)

curves <- getCurves(
  lineages,
  approx_points = 100,
  thresh = 0.01,
  stretch = 0.8,
  allow.breaks = FALSE,
  shrink = 0.99
)

## Same PC basis as the fit (required for predict)
rd_full <- as.matrix(reducedDim(sce, "X_pca"))[, seq_len(ncp), drop = FALSE]
cell_ids <- colnames(sce)
rownames(rd_full) <- cell_ids

curves_all <- tryCatch(
  slingshot::predict(curves, rd_full),
  error = function(e) {
    slingshot::predict(as.SlingshotDataSet(curves), rd_full)
  }
)
pt_full <- as.matrix(slingPseudotime(curves_all, na = TRUE))
rownames(pt_full) <- cell_ids

pseudotime_df <- as.data.frame(pt_full)
pseudotime_df$cell_id <- cell_ids
pseudotime_df$cell_label <- colData(sce)$cell_label
write_tsv(pseudotime_df, paste0(output_prefix, "-pseudotime.tsv"))

pca <- as.data.frame(rd_full[, 1:2, drop = FALSE])
colnames(pca) <- c("PC1", "PC2")
pca$cell_label <- colData(sce)$cell_label
pca$pseudotime <- rowMeans(pt_full, na.rm = TRUE)

sc_list <- slingCurves(curves)
if (!length(sc_list)) {
  stop("Slingshot produced no curves (empty slingCurves()).")
}
cv <- sc_list[[1]]
coords <- cv$s[cv$ord, 1:2, drop = FALSE]
curve_df <- data.frame(PC1 = coords[, 1], PC2 = coords[, 2], lineage = 1)

p_pt <- ggplot(pca, aes(PC1, PC2)) +
  geom_point(aes(colour = pseudotime), size = 0.3, alpha = 0.6) +
  scale_colour_viridis_c(na.value = "grey80") +
  geom_path(data = curve_df, aes(group = lineage), linewidth = 1) +
  theme_minimal() +
  labs(
    title = paste("Slingshot trajectory (pseudotime):", output_prefix),
    colour = "Pseudotime"
  )

ggsave(paste0(output_prefix, "-slingshot.pdf"), p_pt, width = 8, height = 6)

p_ct <- ggplot(pca, aes(PC1, PC2)) +
  geom_path(
    data = curve_df,
    aes(PC1, PC2, group = lineage),
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  geom_point(aes(colour = cell_label), size = 0.35, alpha = 0.75) +
  scale_colour_brewer(type = "qual", palette = "Set1") +
  theme_minimal() +
  labs(
    title = paste("Slingshot trajectory (cell type):", output_prefix),
    colour = "Cell type"
  )

ggsave(paste0(output_prefix, "-slingshot-celltype.pdf"), p_ct, width = 8, height = 6)
