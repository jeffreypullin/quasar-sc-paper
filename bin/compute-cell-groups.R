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

# Quantile bins (equal cell counts per group when continuous)
breaks <- unique(quantile(x, probs = seq(0, 1, length.out = k + 1L), na.rm = TRUE))
bin <- cut(x, breaks = breaks, include.lowest = TRUE, labels = FALSE)

# Equal-width bins on interaction covariate values (e.g. pseudotime)
#x_finite <- x[!is.na(x)]
#if (length(x_finite) == 0) {
#  stop("No finite values in interaction covariate: ", int_cov)
#}
#min_x <- min(x_finite)
#max_x <- max(x_finite)
#if (min_x == max_x) {
#  stop("Cannot form ", k, " groups: ", int_cov, " is constant.")
#}
#breaks <- seq(min_x, max_x, length.out = k + 1L)
#bin <- cut(x, breaks = breaks, include.lowest = TRUE, labels = FALSE)

out <- tibble(
  group = paste0(int_cov, "_q", bin),
  cell_id = covs$cell_id,
  value = x
) |>
  filter(!is.na(bin))

write_tsv(out, paste0(cell_type, "-", int_cov, "-K", k, "-cell-groups.tsv"))
