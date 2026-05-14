#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

cell_type <- args[1]
cov_spec <- args[2]
covs <- read_tsv(args[3], show_col_types = FALSE)

if (cov_spec == "bulk_pca") {
  covs <- covs |>
    select(sample_id, cell_id, sex, age,
           starts_with("geno_pc"), starts_with("PC_"))
} else if (cov_spec == "sc_pca") {
  covs <- covs |>
    select(sample_id, cell_id, sex, age,
           starts_with("geno_pc"), starts_with("scPC_"))
} else if (cov_spec == "bulk_pca+pct_mito") {
  covs <- covs |>
    select(sample_id, cell_id, sex, age,
           starts_with("geno_pc"), starts_with("PC_"), "pct_counts_mt")
} else if (cov_spec == "bulk_pca+cell_cycle") {
  covs <- covs |>
    select(sample_id, cell_id, sex, age,
           starts_with("geno_pc"), starts_with("PC_"), "G2M_score", "S_score")
} else {
  step("Unknown cov_spec.")
}

write_tsv(covs, file = paste0(cell_type, "-", cov_spec, "-covs.tsv"))
