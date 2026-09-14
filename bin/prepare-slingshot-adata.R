#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(hdf5r)
  library(Matrix)
  library(Seurat)
  library(readr)
})

N_HVG <- 500L

COMBINED_MEMBERS <- list(
  B_all = c("B IN", "B Mem"),
  T_all = c(
    "CD4 NC", "CD4 ET", "CD4 SOX4",
    "CD8 ET", "CD8 NC", "CD8 S100B"
  )
)

cell_label_keep <- function(labels, cell_type) {
  if (cell_type %in% names(COMBINED_MEMBERS)) {
    return(labels %in% COMBINED_MEMBERS[[cell_type]])
  }
  labels == gsub("_", " ", cell_type, fixed = TRUE)
}

# Old anndata (<0.8) categorical: integer codes + __categories/<col>.
read_obs_labels <- function(h5, col) {
  codes <- h5[[paste0("obs/", col)]][]
  levels <- as.character(h5[[paste0("obs/__categories/", col)]][])
  levels[as.integer(codes) + 1L]
}

args <- commandArgs(trailingOnly = TRUE)
cell_type <- args[[1]]
h5ad_file <- args[[2]]

h5 <- H5File$new(h5ad_file, mode = "r")
on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)

barcodes <- as.character(h5[["obs/barcode"]][])
cell_labels <- read_obs_labels(h5, "cell_label")
pool <- h5[["obs/pool"]][]
percent_mt <- as.numeric(h5[["obs/percent.mt"]][])
gene_ids <- as.character(h5[["var/Geneid"]][])

keep_idx <- which(cell_label_keep(cell_labels, cell_type))
if (!length(keep_idx)) {
  stop(sprintf("No cells matched cell_type=%s", cell_type))
}

indptr <- as.integer(h5[["X/indptr"]][])
x_data <- h5[["X/data"]][]
x_indices <- as.integer(h5[["X/indices"]][])
h5$close_all()
on.exit(NULL)

nnz <- indptr[keep_idx + 1L] - indptr[keep_idx]
sel <- sequence(nnz, from = indptr[keep_idx] + 1L)
counts <- sparseMatrix(
  i = x_indices[sel] + 1L,
  j = rep.int(seq_along(keep_idx), nnz),
  x = as.numeric(x_data[sel]),
  dims = c(length(gene_ids), length(keep_idx)),
  dimnames = list(gene_ids, barcodes[keep_idx])
)
rm(x_data, x_indices, indptr, sel, nnz)
invisible(gc())

meta <- data.frame(
  cell_label = cell_labels[keep_idx],
  pool = factor(pool[keep_idx]),
  percent.mt = percent_mt[keep_idx],
  row.names = barcodes[keep_idx],
  stringsAsFactors = FALSE
)

obj <- CreateSeuratObject(counts = counts, meta.data = meta)
rm(counts, meta)
invisible(gc())

obj <- SCTransform(
  obj,
  vars.to.regress = c("pool", "percent.mt"),
  do.correct.umi = FALSE,
  variable.features.n = N_HVG,
  verbose = TRUE
)

scale_data <- GetAssayData(obj, assay = "SCT", layer = "scale.data")

out <- data.frame(
  cell_id = colnames(obj),
  cell_label = as.character(obj$cell_label),
  t(as.matrix(scale_data)),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

out_path <- paste0(cell_type, "-slingshot-input.tsv.gz")
write_tsv(out, out_path)
