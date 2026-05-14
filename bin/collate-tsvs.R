#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
tsvs <- str_remove_all(args, "[,\\[\\]]")
all_data <- bind_rows(lapply(tsvs, function(x) read_tsv(x, show_col_types = FALSE)))

write_tsv(all_data, "saigeqtl-files.tsv")
