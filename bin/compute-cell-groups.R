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

# Within-donor residual of the interaction covariate (same construction as
# Quasar's add_bw_covariates() for *_w). Ranking by raw x within donor is
# identical to ranking by x_w; we store x_w as `value` so group_linear scores
# track within-donor position rather than global PT.
#
# Previously: global quantile bins of raw x (equal cell counts across the
# whole dataset). Kept here for reference:
# breaks <- unique(quantile(x, probs = seq(0, 1, length.out = k + 1L), na.rm = TRUE))
# bin <- cut(x, breaks = breaks, include.lowest = TRUE, labels = FALSE)

out <- covs |>
  transmute(
    sample_id,
    cell_id,
    x = .data[[int_cov]]
  ) |>
  filter(is.finite(x)) |>
  group_by(sample_id) |>
  mutate(
    x_w = x - mean(x),
    bin = ntile(x, k)
  ) |>
  ungroup() |>
  filter(!is.na(bin)) |>
  transmute(
    group = paste0(int_cov, "_q", bin),
    cell_id,
    value = x_w
  )

n_groups <- n_distinct(out$group)
if (n_groups < 2L) {
  stop("Within-donor binning produced fewer than 2 groups; check K and cell counts.")
}

write_tsv(out, paste0(cell_type, "-", int_cov, "-K", k, "-cell-groups.tsv"))
