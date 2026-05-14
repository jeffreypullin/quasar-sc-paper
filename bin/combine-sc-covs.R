#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

cell_type <- args[1]
pb_expr_covs <- read_tsv(args[2], show_col_types = FALSE)
sc_expr_covs <- read_tsv(args[3], show_col_types = FALSE)
geno_pcs  <- read_tsv(args[4], show_col_types = FALSE)

pb_covs <- left_join(pb_expr_covs, geno_pcs, by = "sample_id") |>
  select(sample_id, sex, age, paste0("PC_", 1:5), paste0("geno_pc", 1:6))

covs <- sc_expr_covs |>
  select(sample_id, cell_id,  paste0("scPC_", 1:5), pct_counts_mt, S_score, G2M_score) |>
  left_join(pb_covs, by = join_by(sample_id)) |>
  select(sample_id, cell_id, everything())

write_tsv(covs, file = paste0(cell_type, "-covs.tsv"))
