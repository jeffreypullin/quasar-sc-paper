#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(hdf5r)
  library(sctransform)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: compute-sc-sct-counts.R <counts_tsv> <h5ad> <cell_type>")
}

counts_file <- args[[1]]
h5ad_file <- args[[2]]
cell_type <- args[[3]]

# Phenotype TSV: sample_id, cell_id, then gene columns (cells x genes).
counts <- fread(counts_file, sep = "\t", header = TRUE, data.table = TRUE)
if (!all(c("sample_id", "cell_id") %in% names(counts))) {
  stop("Counts TSV must have sample_id and cell_id columns.")
}

sample_ids <- counts$sample_id
cell_ids <- as.character(counts$cell_id)
gene_ids <- setdiff(names(counts), c("sample_id", "cell_id"))
umi_dense <- as.matrix(counts[, ..gene_ids])
storage.mode(umi_dense) <- "double"
rownames(umi_dense) <- cell_ids
# sctransform::vst expects genes as rows, cells as columns.
umi <- Matrix(t(umi_dense), sparse = TRUE)
rm(umi_dense, counts)
invisible(gc(verbose = FALSE))

# Read obs only from the h5ad (do not touch X).
h5 <- H5File$new(h5ad_file, mode = "r")
on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)

read_obs_vec <- function(path) {
  ds <- h5[[path]]
  as.vector(ds[])
}

barcodes <- as.character(read_obs_vec("obs/barcode"))
pool <- read_obs_vec("obs/pool")
percent_mt <- as.numeric(read_obs_vec("obs/percent.mt"))
n_count_rna <- as.numeric(read_obs_vec("obs/nCount_RNA"))
h5$close_all()
on.exit(NULL)

obs_idx <- match(cell_ids, barcodes)
if (anyNA(obs_idx)) {
  n_miss <- sum(is.na(obs_idx))
  stop(sprintf("%d cell_id values from counts TSV not found in h5ad obs/barcode.", n_miss))
}

cell_attr <- data.frame(
  log_umi = log(n_count_rna[obs_idx]),
  percent.mt = percent_mt[obs_idx],
  batch = factor(pool[obs_idx]),
  row.names = cell_ids,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (any(!is.finite(cell_attr$log_umi)) || any(n_count_rna[obs_idx] <= 0)) {
  stop("nCount_RNA must be positive and finite for every cell in the counts TSV.")
}
if (any(!is.finite(cell_attr$percent.mt))) {
  stop("percent.mt must be finite for every cell in the counts TSV.")
}
if (nlevels(cell_attr$batch) < 2L) {
  stop("batch (pool) must have at least 2 levels for sctransform batch_var.")
}

vst_out <- vst(
  umi,
  cell_attr = cell_attr,
  latent_var = c("log_umi", "percent.mt"),
  batch_var = "batch",
  method = "poisson",
  return_cell_attr = TRUE,
  verbosity = 1
)
corrected <- correct_counts(vst_out, umi, verbosity = 1)
# Cells x genes again, matching the SC phenotype layout.
corrected_cells <- t(as.matrix(corrected))
storage.mode(corrected_cells) <- "integer"

out <- data.table(
  sample_id = sample_ids,
  cell_id = cell_ids
)
out <- cbind(out, as.data.table(corrected_cells))

out_path <- paste0(cell_type, "-sct-sc-pheno.tsv")
fwrite(out, out_path, sep = "\t", quote = FALSE)
