#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(hdf5r)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[[1]]
h5ad_file <- args[[2]]
dataset <- args[[3]]
cell_type <- args[[4]]

counts <- fread(counts_file)
sample_ids <- counts$sample_id
cell_ids <- as.character(counts$cell_id)
gene_ids <- setdiff(names(counts), c("sample_id", "cell_id"))
umi <- Matrix(
  t(as.matrix(counts[, ..gene_ids])),
  sparse = TRUE,
  dimnames = list(gene_ids, cell_ids)
)
rm(counts)
gc()

h5 <- H5File$new(h5ad_file, mode = "r")
barcodes <- as.character(h5[["obs/barcode"]][])
pool <- h5[["obs/pool"]][]
percent_mt <- as.numeric(h5[["obs/percent.mt"]][])
h5$close_all()

obs_idx <- match(cell_ids, barcodes)
obj <- CreateSeuratObject(
  counts = umi,
  meta.data = data.frame(
    pool = pool[obs_idx],
    percent.mt = percent_mt[obs_idx],
    row.names = cell_ids
  )
)
rm(umi)

# OneK1K norm.R. These two arguments keep every input gene in both outputs.
obj <- SCTransform(
  object = obj,
  vars.to.regress = c("percent.mt", "pool"),
  return.only.var.genes = FALSE,
  min_cells = 1
)

write_pheno <- function(mat, suffix) {
  cells <- t(as.matrix(mat))
  fwrite(
    cbind(data.table(sample_id = sample_ids, cell_id = cell_ids), as.data.table(cells)),
    paste0(dataset, "-", cell_type, "-", suffix, ".tsv"),
    sep = "\t"
  )
}

write_pheno(GetAssayData(obj, assay = "SCT", slot = "counts"), "sct-sc-pheno")
write_pheno(GetAssayData(obj, assay = "SCT", slot = "scale.data"), "sct-residuals-sc-pheno")
