#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

cell_type <- args[[1]]
pb_covs <- read_tsv(args[[2]], show_col_types = FALSE)
sc_expr_covs <- read_tsv(args[[3]], show_col_types = FALSE)
int_cov_data <- read_tsv(args[[4]], show_col_types = FALSE)

if ("Lineage1" %in% names(int_cov_data)) {
  int_cov_data <- int_cov_data |>
    rename(pseudotime = Lineage1) |>
    select(-any_of("cell_label"))
}

pb_int_cov <- sc_expr_covs |>
  select(sample_id, cell_id) |>
  left_join(int_cov_data, by = "cell_id") |>
  group_by(sample_id) |>
  summarise(pseudotime = mean(pseudotime, na.rm = TRUE), .groups = "drop")

out <- pb_covs |>
  left_join(pb_int_cov, by = "sample_id")

write_tsv(out, paste0(cell_type, "-pb-cov-data-with-int-cov.tsv"))
