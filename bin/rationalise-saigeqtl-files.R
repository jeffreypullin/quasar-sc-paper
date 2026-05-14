#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

input_tsv <- read_tsv(args[[1]], show_col_types = FALSE)

cell_type <- input_tsv$cell_type[[1]]
chrom <- input_tsv$chrom[[1]]
gene_list <- str_replace(
  args[[1]],
  "^(.+?)-\\1-(genes-\\d+\\.txt).*$",
  "\\1-\\2"
)

all_variant_data <- tibble()
all_region_data <- tibble()
error_genes <- c()
for (i in seq_len(nrow(input_tsv))) {

  gene <- input_tsv$gene[[i]]

  if (is.na(input_tsv$variant_file[[i]]) || is.na(input_tsv$region_file[[i]])) {
    error_genes <- c(error_genes, gene)
    next
  }

  variant_data <- read_tsv(
    input_tsv$variant_file[[i]],
    show_col_types = FALSE
  ) |>
    mutate(gene = gene)

  region_data <- read_tsv(input_tsv$region_file[[i]], show_col_types = FALSE) |>
    mutate(gene = gene)

  all_variant_data <- bind_rows(all_variant_data, variant_data)
  all_region_data <- bind_rows(all_region_data, region_data)

}

step1_time <- sum(input_tsv$step1_time)
step2_time <- sum(input_tsv$step2_time)
step3_time <- sum(input_tsv$step3_time)

prefix <- str_remove(args[[1]], "-saigeqtl-files\\.tsv$")

all_variant_file <- paste0("rationalised-", prefix, "-variant.tsv")
all_region_file <- paste0("rationalised-", prefix, "-region.tsv")

write_tsv(all_variant_data, all_variant_file)
write_tsv(all_region_data, all_region_file)

result_tibble <- tibble(
  cell_type = cell_type,
  chrom = chrom,
  gene_list = gene_list,
  step1_time = step1_time,
  step2_time = step2_time,
  step3_time = step3_time,
  variant_file = paste0(getwd(), "/", all_variant_file),
  region_file = paste0(getwd(), "/", all_region_file),
  error_genes = list(error_genes)
)

write_tsv(result_tibble, paste0("rationalised-", args[[1]]))

all_files <- list.files(".")
files_to_delete <- all_files[!str_starts(all_files, "rationalised-")]
file.remove(files_to_delete)
