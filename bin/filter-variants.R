#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
chrom <- args[1]
cell_type <- args[2]
variant_data <- read_tsv(args[3], show_col_types = FALSE)
prune_in <- read_lines(args[4])

filt_variant_data <- variant_data |>
  filter(snp_id %in% prune_in)

write_tsv(
  filt_variant_data,
  paste0(chrom, "-", cell_type, "-quasar-cis-variant-filt.tsv")
)
