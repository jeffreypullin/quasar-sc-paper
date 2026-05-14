#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
chr <- args[[1]]
cell_type <- args[[2]]
int_cov <- args[[3]]

clump_files <- list.files(pattern = "*.clumps")
feature_ids <- str_extract(clump_files, "ENSG[0-9]+")

clumps_data <- tibble(file_name = clump_files) |>
  mutate(feature_id = str_extract(file_name, "ENSG[0-9]+")) |>
  rowwise() |>
  mutate(data = list(read_tsv(file_name, show_col_types = FALSE))) |>
  unnest(cols = c(data)) |>
  select(-file_name)

write_tsv(
  clumps_data, 
  paste0(chr, "-", cell_type, "-", int_cov, "-harmonised-clumps.tsv")
)
