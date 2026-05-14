#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

cell_type <- args[1]
expr_covs_file <- args[2]
geno_pcs_file <- args[3]

expr_covs <- read_tsv(expr_covs_file, show_col_types = FALSE)
geno_pcs <- read_tsv(geno_pcs_file, show_col_types = FALSE)

covs <- left_join(expr_covs, geno_pcs, by = "sample_id") |>
  select(sample_id, sex, age, paste0("PC_", 1:2), paste0("geno_pc", 1:6)) 

write_tsv(covs, file = paste0(cell_type, "-covs.tsv"))
