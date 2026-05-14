#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

cell_type <- args[1]
pb_counts_file <- args[2]
annot_bed_file <- args[3]

pb_counts <- read_tsv(pb_counts_file, show_col_types = FALSE)
annot_bed <- read_tsv(annot_bed_file, show_col_types = FALSE)

annot_pheno <- left_join(
  pb_counts,
  annot_bed,
  by = join_by("feature_id" == "phenotype_id")
)

annot_pheno <- annot_pheno |>
  rename(phenotype_id = feature_id) |>
  select(`#chr`, start, end, phenotype_id, everything()) |>
  # Filter out non-autosomal genes which don't have inforamtion.
  filter(!is.na(start))

write_tsv(annot_pheno, file = paste0(cell_type, "-annot-pheno.tsv"))
