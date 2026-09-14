#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

cell_type <- args[[1]]
cov_data <- read_tsv(args[[2]], show_col_types = FALSE)
int_cov_data <- read_tsv(args[[3]], show_col_types = FALSE)

if ("pseudotime" %in% names(int_cov_data)) {
  int_cov_data <- int_cov_data |>
    select(cell_id, pseudotime)
} else if ("Lineage1" %in% names(int_cov_data)) {
  int_cov_data <- int_cov_data |>
    rename(pseudotime = Lineage1) |>
    select(-cell_label)
}

out <- cov_data |>
  left_join(int_cov_data, by = "cell_id")

write_tsv(out, paste0(cell_type, "-cov-data-with-int-cov.tsv"))
