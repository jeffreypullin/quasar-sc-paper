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
  select(sample_id, sex, age, paste0("PC_", 1:2), paste0("geno_pc", 1:6))

int_cov_cols <- c(
  "pseudotime",
  "starcat_Cytotoxic", "starcat_Th22", "starcat_MAIT", "starcat_TEMRA",
  "starcat_CD4_CM", "starcat_CD8_EM", "starcat_CD4_Naive",
  "starcat_Th2_Activated", "starcat_Th2_Resting", "starcat_Th1_Like",
  "starcat_CD8_Trm", "starcat_Th17_Activated", "starcat_Tfh_2",
  "starcat_Tph", "starcat_Exhaustion", "starcat_Tfh_1",
  "starcat_CellCycle_S", "starcat_CellCycle_G2M"
)

covs <- sc_expr_covs |>
  select(sample_id, cell_id, paste0("scPC_", 1:2), pct_counts_mt, S_score, G2M_score,
         any_of(int_cov_cols)) |>
  left_join(pb_covs, by = join_by(sample_id)) |>
  select(sample_id, cell_id, everything())

write_tsv(covs, file = paste0(cell_type, "-covs.tsv"))
