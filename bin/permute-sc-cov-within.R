#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
ind <- args[[1]]
cell_type <- args[[2]]
int_cov <- args[[3]]
cov <- args[[4]]

set.seed(as.numeric(ind) * 10 + 64)

cov_data <- read_delim(cov, show_col_types = FALSE)
perm_col_name <- paste0(int_cov, "_perm")

cov_data <- cov_data |>
  group_by(sample_id) |>
  mutate(!!perm_col_name := sample(.data[[int_cov]])) |>
  ungroup()

write_delim(
  cov_data,
  paste0("permute-within-", cell_type, "-", int_cov, "-", ind, ".tsv"),
  delim = "\t",
  col_names = TRUE
)
