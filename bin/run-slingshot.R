#!/usr/bin/env Rscript

# Use the process conda env's Python; otherwise reticulate prefers
# ~/.virtualenvs/r-reticulate, which has no phate.
Sys.setenv(RETICULATE_PYTHON = Sys.which("python"))

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(slingshot)
  library(phateR)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
input_tsv <- args[1]
dataset <- args[2]
cell_type <- args[3]
# No fourth argument: the published trajectory (PHATE seed 20261013).
# A fourth argument re-draws PHATE and writes seed-tagged files.
seed <- if (length(args) >= 4L) as.integer(args[4]) else 20261013L
output_prefix <- paste(dataset, cell_type, sep = "-")
if (length(args) >= 4L) {
  output_prefix <- paste(output_prefix, paste0("seed", seed), sep = "-")
}
n_pca <- 10L

meta_cols <- c("cell_id", "cell_label")
dat <- read_tsv(input_tsv, show_col_types = FALSE)
cell_ids <- dat$cell_id
gene_cols <- setdiff(colnames(dat), meta_cols)

res_mat <- as.matrix(dat[, gene_cols, drop = FALSE])
rownames(res_mat) <- cell_ids

phate_op <- phate(res_mat, npca = n_pca, seed = seed)
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

sc_list <- slingCurves(sds)
cv <- sc_list[[curve_idx]]
coords <- cv$s[cv$ord, 1:2, drop = FALSE]

cells <- data.frame(
  kind = "cell",
  PHATE_1 = rd_mat[, 1],
  PHATE_2 = rd_mat[, 2],
  cell_label = as.character(colData(sce)$cell_label),
  pseudotime = pseudotime,
  lineage = NA_integer_,
  stringsAsFactors = FALSE
)
curve <- data.frame(
  kind = "curve",
  PHATE_1 = coords[, 1],
  PHATE_2 = coords[, 2],
  cell_label = NA_character_,
  pseudotime = NA_real_,
  lineage = curve_idx,
  stringsAsFactors = FALSE
)
write_tsv(rbind(cells, curve), paste0(output_prefix, "-slingshot-plot.tsv"))
