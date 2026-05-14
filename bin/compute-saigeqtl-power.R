#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
files_data <- read_tsv(args[1], show_col_types = FALSE)

compute_variant_power <- function(file) {
  if (file == "NA" || is.na(file)) {
    out <- NA
  } else {
    out <- sum(fread(file, select = "p.value")$p.value < 5e-6, na.rm = TRUE)
  }
  out
}

n_sig_data <- files_data |>
  rowwise() |>
  mutate(n_sig_variant = compute_variant_power(variant_file)) |>
  ungroup()

write_tsv(n_sig_data, basename(args[1]))
