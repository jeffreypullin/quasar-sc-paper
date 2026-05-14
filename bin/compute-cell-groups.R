#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

cov_file <- args[[1]]
int_cov <- args[[2]]
k <- as.integer(args[[3]])
cell_type <- args[[4]]

if (is.na(k) || k < 2L) {
  stop("K must be an integer >= 2")
}

covs <- read_tsv(cov_file, show_col_types = FALSE)

x <- covs[[int_cov]]
breaks <- unique(quantile(x, probs = seq(0, 1, length.out = k + 1L), na.rm = TRUE))
bin <- cut(x, breaks = breaks, include.lowest = TRUE, labels = FALSE)

out <- tibble(
  group = paste0(int_cov, "_q", bin),
  cell_id = covs$cell_id
) |>
  filter(!is.na(bin))

write_tsv(out, paste0(cell_type, "-", int_cov, "-K", k, "-cell-groups.tsv"))
