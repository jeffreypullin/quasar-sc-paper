#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

ids <- read_csv(args[1], col_names = FALSE, show_col_types = FALSE)[[1]]
pheno <- read_tsv(args[2], show_col_types = FALSE)
anno <- read_tsv(args[3], show_col_types = FALSE)

# NB: This is not the same set of genes as the bulk examples.
out_pheno <- pheno |>
  filter(sample_id %in% ids) |>
  select(sample_id, cell_id, 3:22)

genes <- colnames(out_pheno)[3:22]

out_anno <- anno |>
  filter(phenotype_id %in% genes)

write_tsv(out_pheno, "sc-pheno-n100.tsv")
write_tsv(out_anno, "example-anno.tsv")
