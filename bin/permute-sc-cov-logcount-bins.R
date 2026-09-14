#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

CHUNK_ROWS <- 1000L

read_depth_from_counts <- function(counts_file, chunk_rows = CHUNK_ROWS) {
  acc <- new.env(parent = emptyenv())
  acc$parts <- list()
  acc$i <- 0L

  callback <- SideEffectChunkCallback$new(function(x, pos) {
    acc$i <- acc$i + 1L
    gene_cols <- setdiff(names(x), c("sample_id", "cell_id"))
    totals <- rowSums(x[gene_cols], na.rm = TRUE)
    acc$parts[[acc$i]] <- tibble(
      sample_id = as.character(x$sample_id),
      cell_id = as.character(x$cell_id),
      log_counts = if_else(totals > 0, log(totals), 0)
    )
    NULL
  })

  read_tsv_chunked(
    counts_file,
    callback = callback,
    chunk_size = chunk_rows,
    show_col_types = FALSE,
    progress = FALSE
  )

  bind_rows(acc$parts)
}

args <- commandArgs(trailingOnly = TRUE)
ind <- args[[1]]
cell_type <- args[[2]]
int_cov <- args[[3]]
cov_file <- args[[4]]
counts_file <- args[[5]]
n_bins <- if (length(args) >= 6) as.integer(args[[6]]) else 10L

set.seed(as.numeric(ind) * 10 + 64)

depth <- read_depth_from_counts(counts_file)

cov_data <- read_delim(cov_file, show_col_types = FALSE)
if (!int_cov %in% names(cov_data)) {
  stop("Interaction covariate '", int_cov, "' not found in ", cov_file)
}
cov_data <- cov_data |>
  mutate(
    sample_id = as.character(sample_id),
    cell_id = as.character(cell_id)
  )

merged <- inner_join(cov_data, depth, by = c("sample_id", "cell_id"))
if (nrow(merged) == 0) {
  stop("No overlapping sample_id/cell_id rows between cov and counts")
}

n_finite <- sum(is.finite(merged$log_counts))
if (n_finite == 0) {
  stop("No finite log_counts values available for binning")
}
if (n_finite < n_bins) {
  n_bins <- max(1L, n_finite)
}

perm_col <- paste0(int_cov, "_perm")
merged <- merged |>
  mutate(
    log_count_bin = ntile(rank(log_counts, ties.method = "first"), n_bins),
    !!perm_col := .data[[int_cov]]
  ) |>
  group_by(log_count_bin) |>
  mutate(
    !!perm_col := if (n() >= 2) sample(.data[[int_cov]]) else .data[[int_cov]]
  ) |>
  ungroup() |>
  select(-log_counts, -log_count_bin)

write_delim(
  merged,
  paste0("permute-logcount-bins-", cell_type, "-", int_cov, "-", ind, ".tsv"),
  delim = "\t",
  col_names = TRUE
)
