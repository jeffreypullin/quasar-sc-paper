#!/usr/bin/env Rscript

# Build the three log-scale offset files used by the offset-precision experiment.
#
# percell   : offset_ic = log(n_ic), the per-cell library size quasar computes internally.
# donorflat : offset_ic = log(N_i / C_i), constant within a donor. Total exposure per
#             donor, sum_c exp(offset_ic) = N_i, is preserved exactly, so this is the
#             pseudobulk-equivalent offset spread uniformly over that donor's cells.
# constant  : offset_ic = log(N / C), the same value for every cell, so exposure is
#             proportional to cell count alone.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
})

args <- commandArgs(trailingOnly = TRUE)
dataset <- args[1]
cell_type <- args[2]
sc_pheno_file <- args[3]
prefix <- paste0(dataset, "-", cell_type, "-offset-")

depth <- read_tsv_chunked(
  sc_pheno_file,
  DataFrameCallback$new(function(chunk, pos) {
    tibble(
      sample_id = chunk$sample_id,
      cell_id = chunk$cell_id,
      total = rowSums(select(chunk, -sample_id, -cell_id))
    )
  }),
  chunk_size = 1000,
  col_types = cols(sample_id = col_character(), cell_id = col_character(), .default = col_double())
)

n_zero_cells <- sum(depth$total <= 0)
if (n_zero_cells > 0) {
  message(
    "Warning: ", n_zero_cells, " cells have zero total counts. Assigning them ",
    "offset 0 in the per-cell file, matching quasar's internal behaviour. ",
    "Donor total exposure is therefore not preserved exactly for the affected donors."
  )
}

safe_log <- function(x) if_else(x > 0, log(x), 0)

offsets <- depth |>
  mutate(
    percell = safe_log(total),
    donorflat = safe_log(mean(total)),
    .by = sample_id
  ) |>
  mutate(constant = safe_log(mean(total)))

walk(c("percell", "donorflat", "constant"), \(spec) {
  offsets |>
    select(sample_id, cell_id, offset = all_of(spec)) |>
    write_tsv(paste0(prefix, spec, ".tsv"))
})
