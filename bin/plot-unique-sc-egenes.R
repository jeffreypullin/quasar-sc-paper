#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)

pb_data_files <- read_tsv(args[1], show_col_types = FALSE)
sc_data_files <- read_tsv(args[2], show_col_types = FALSE)

cell_types <- c("Plasma", "B_IN", "CD4_NC")

pb_manifest <- pb_data_files |>
  filter(
    model == "nb_glm",
    int_cov == "none",
    cell_frac == 1,
    indiv_frac == 1,
    n_cells_target < 0,
    count_frac == 1,
    cell_type %in% cell_types
  )

sc_manifest <- sc_data_files |>
  filter(
    model == "p_glmm_sc",
    cov_spec == "bulk_pca",
    k == "none",
    cell_frac == 1,
    indiv_frac == 1,
    n_cells_target < 0,
    count_frac == 1,
    cell_type %in% cell_types
  )

pb_genes <- pb_manifest |>
  mutate(method = "nb_glm") |>
  rowwise() |>
  mutate(gene_tbl = list(read_tsv(region_file, show_col_types = FALSE))) |>
  unnest(gene_tbl) |>
  ungroup() |>
  select(cell_type, method, feature_id, region_pvalue = pvalue)

sc_genes <- sc_manifest |>
  mutate(method = "p_glmm_sc") |>
  rowwise() |>
  mutate(gene_tbl = list(read_tsv(region_file, show_col_types = FALSE))) |>
  unnest(gene_tbl) |>
  ungroup() |>
  select(cell_type, method, feature_id, region_pvalue = pvalue)

all_genes <- bind_rows(pb_genes, sc_genes) |>
  mutate(
    region_bh = p.adjust(region_pvalue, method = "BH"),
    .by = c(cell_type, method)
  )

sig_genes <- all_genes |>
  filter(region_bh < 0.05) |>
  select(cell_type, method, feature_id, region_pvalue, region_bh)

unique_sc_egenes <- sig_genes |>
  filter(method == "p_glmm_sc") |>
  anti_join(
    sig_genes |>
      filter(method == "nb_glm") |>
      select(cell_type, feature_id),
    by = c("cell_type", "feature_id")
  ) |>
  select(cell_type, feature_id, region_pvalue, region_bh)

variant_col_types <- cols(
  feature_id = col_character(),
  snp_id = col_character(),
  chrom = col_character(),
  pos = col_double(),
  alt = col_character(),
  ref = col_character(),
  maf = col_double(),
  beta = col_double(),
  se = col_double(),
  pvalue = col_double()
)

sc_variants <- sc_manifest |>
  summarise(
    variant_files = list(variant_file),
    .by = cell_type
  ) |>
  mutate(
    var_data = map(variant_files, function(files) {
      map_dfr(files, function(path) {
        read_tsv(path, show_col_types = FALSE, col_types = variant_col_types) |>
          select(feature_id, snp_id, chrom, pos, maf, beta, se, pvalue)
      })
    })
  ) |>
  select(cell_type, var_data) |>
  unnest(var_data)

lead_variants <- unique_sc_egenes |>
  left_join(sc_variants, by = c("cell_type", "feature_id")) |>
  filter(!is.na(pvalue), pvalue > 0) |>
  group_by(cell_type, feature_id) |>
  slice_min(order_by = pvalue, n = 1, with_ties = FALSE) |>
  ungroup() |>
  arrange(cell_type, region_bh, feature_id) |>
  select(
    cell_type, feature_id, snp_id, chrom, pos, maf, beta, se, pvalue,
    region_pvalue, region_bh
  )

write_tsv(lead_variants, "unique-p-glmm-egenes.tsv")

lead_variants |>
  print(n = Inf)
