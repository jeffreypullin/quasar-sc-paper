#!/usr/bin/env Rscript

# Use the process conda env's Python; otherwise reticulate prefers
# ~/.virtualenvs/r-reticulate, which has no phate.
Sys.setenv(RETICULATE_PYTHON = Sys.which("python"))

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(slingshot)
  library(phateR)
  library(ggplot2)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
input_tsv <- args[1]
output_prefix <- args[2]
n_pca <- 10L

meta_cols <- c("cell_id", "cell_label")
dat <- read_tsv(input_tsv, show_col_types = FALSE)
cell_ids <- dat$cell_id
gene_cols <- setdiff(colnames(dat), meta_cols)

res_mat <- as.matrix(dat[, gene_cols, drop = FALSE])
rownames(res_mat) <- cell_ids

phate_op <- phate(res_mat, npca = n_pca)
rd_mat <- phate_op$embedding
colnames(rd_mat) <- paste0("PHATE_", seq_len(ncol(rd_mat)))
rownames(rd_mat) <- cell_ids

cd <- data.frame(
  cell_label = dat$cell_label,
  row.names = cell_ids,
  stringsAsFactors = FALSE
)

sce <- SingleCellExperiment(
  assays = list(X = matrix(0, nrow = 1L, ncol = length(cell_ids))),
  colData = cd,
  reducedDims = list(PHATE = rd_mat)
)
colnames(sce) <- cell_ids

sds <- slingshot(
  rd_mat,
  clusterLabels = colData(sce)$cell_label,
  start.clus = "B IN",
  end.clus = "B Mem"
)

pt_full <- as.matrix(slingPseudotime(sds, na = TRUE))
rownames(pt_full) <- cell_ids

curve_idx <- 1L
pseudotime <- pt_full[, curve_idx, drop = TRUE]

pseudotime_df <- as.data.frame(pt_full)
colnames(pseudotime_df) <- paste0("Lineage", seq_len(ncol(pseudotime_df)))
pseudotime_df$cell_id <- cell_ids
pseudotime_df$cell_label <- colData(sce)$cell_label
pseudotime_df$pseudotime <- pseudotime
write_tsv(pseudotime_df, paste0(output_prefix, "-pseudotime.tsv"))

plot_df <- as.data.frame(rd_mat[, 1:2, drop = FALSE])
colnames(plot_df) <- c("PHATE_1", "PHATE_2")
plot_df$cell_label <- colData(sce)$cell_label
plot_df$pseudotime <- pseudotime

sc_list <- slingCurves(sds)
cv <- sc_list[[curve_idx]]
coords <- cv$s[cv$ord, 1:2, drop = FALSE]
curve_df <- data.frame(PHATE_1 = coords[, 1], PHATE_2 = coords[, 2], lineage = curve_idx)

p_pt <- ggplot(plot_df, aes(PHATE_2, PHATE_1)) +
  geom_point(aes(colour = pseudotime), size = 0.3, alpha = 0.6) +
  scale_colour_viridis_c(na.value = "grey80") +
  geom_path(data = curve_df, aes(PHATE_2, PHATE_1, group = lineage), linewidth = 1) +
  theme_minimal() +
  labs(
    title = paste("Slingshot trajectory (pseudotime):", output_prefix),
    x = "PHATE 2",
    y = "PHATE 1",
    colour = "Pseudotime"
  )

ggsave(
  paste0(output_prefix, "-slingshot.pdf"),
  p_pt,
  width = 8,
  height = 6
)

p_ct <- ggplot(plot_df, aes(PHATE_2, PHATE_1)) +
  geom_path(
    data = curve_df,
    aes(PHATE_2, PHATE_1, group = lineage),
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  geom_point(aes(colour = cell_label), size = 0.35, alpha = 0.75) +
  scale_colour_brewer(type = "qual", palette = "Set1") +
  theme_minimal() +
  labs(
    title = paste("Slingshot trajectory (cell type):", output_prefix),
    x = "PHATE 2",
    y = "PHATE 1",
    colour = "Cell type"
  )

ggsave(
  paste0(output_prefix, "-slingshot-celltype.pdf"),
  p_ct,
  width = 8,
  height = 6
)
