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

gene_cols <- setdiff(colnames(sc_counts), c("sample_id", "cell_id"))
totals <- rowSums(sc_counts[, gene_cols, drop = FALSE])

# Library-size offset from all genes, before restricting to this chromosome.
sc_counts <- sc_counts |>
  mutate(log_cell_read_counts = ifelse(totals > 0, log(totals), 0)) |>
  select(sample_id, log_cell_read_counts, any_of(chr_genes))

covs <- left_join(expr_covs, geno_pcs, by = "sample_id") |>
  select(sample_id, sex, age, paste0("PC_", 1:5), paste0("geno_pc", 1:6))

all_data <- left_join(covs, sc_counts, by = "sample_id")

write_tsv(
  all_data,
  file = paste0(safe_cell_type, "-", chr, "-saigeqtl-input.tsv")
)
