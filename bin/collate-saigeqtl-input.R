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
expr_covs_file <- args[4]
geno_pcs_file <- args[5]
anno_file <- args[6]
safe_cell_type <- args[7]

sc_counts <- read_tsv(sc_counts_file, show_col_types = FALSE)
expr_covs <- read_tsv(expr_covs_file, show_col_types = FALSE)
geno_pcs <- read_tsv(geno_pcs_file, show_col_types = FALSE)
anno <- read_tsv(anno_file, show_col_types = FALSE)

chr_genes <- anno |>
  filter(`#chr` == as.numeric(str_sub(chr, 4, 5))) |>
  pull(phenotype_id)

sc_counts <- sc_counts |>
  select(sample_id, any_of(chr_genes))

covs <- left_join(expr_covs, geno_pcs, by = "sample_id") |>
  select(sample_id, sex, age, paste0("PC_", 1:5), paste0("geno_pc", 1:6))

all_data <- left_join(covs, sc_counts, by = "sample_id") |>
  mutate(size_factor = log(rowSums(pick(-(1:14)))))

write_tsv(
  all_data,
  file = paste0(safe_cell_type, "-", chr, "-saigeqtl-input.tsv")
)
