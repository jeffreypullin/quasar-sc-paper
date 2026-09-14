#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

cell_type <- args[1]
chr <- args[2]
sc_counts_file <- args[3]
covs_file <- args[4]
anno_file <- args[5]
safe_cell_type <- args[6]
int_cov <- args[7]

sc_counts <- read_tsv(sc_counts_file, show_col_types = FALSE)
covs <- read_tsv(covs_file, show_col_types = FALSE)
anno <- read_tsv(anno_file, show_col_types = FALSE)

chr_genes <- anno |>
  filter(`#chr` == as.numeric(str_sub(chr, 4, 5))) |>
  pull(phenotype_id)

gene_cols <- setdiff(colnames(sc_counts), c("sample_id", "cell_id"))
totals <- rowSums(sc_counts[, gene_cols, drop = FALSE])

sc_counts <- sc_counts |>
  mutate(log_cell_read_counts = ifelse(totals > 0, log(totals), 0)) |>
  select(sample_id, cell_id, log_cell_read_counts, any_of(chr_genes))

covs <- covs |>
  select(cell_id, sex, age, paste0("PC_", 1:2), paste0("geno_pc", 1:6), all_of(int_cov))

all_data <- inner_join(sc_counts, covs, by = "cell_id") |>
  filter(!is.na(.data[[int_cov]]), !is.na(sample_id)) |>
  select(
    sample_id, cell_id, sex, age, paste0("PC_", 1:2), paste0("geno_pc", 1:6),
    all_of(int_cov), log_cell_read_counts, any_of(chr_genes)
  )

write_tsv(
  all_data,
  file = paste0(safe_cell_type, "-", chr, "-castie-input.tsv")
)
