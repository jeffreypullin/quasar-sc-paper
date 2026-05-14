#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

cell_type <- args[1]
chr <- args[2]
sc_counts_file <- args[3]
n <- as.numeric(args[4])

sc_counts <- read_tsv(sc_counts_file, n_max = 1, show_col_types = FALSE)

genes <- colnames(sc_counts)
genes <- genes[str_sub(genes, 1, 4) == "ENSG"]

genes_split <- split(genes, ceiling(seq_along(genes) / n))
for (i in seq_along(genes_split)) {
  grp <- genes_split[[i]]
  writeLines(grp, paste0(cell_type, "-", chr, "-genes-", i, ".txt"))
}
