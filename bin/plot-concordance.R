#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggupset)
  library(tidyr)
  library(patchwork)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)

pb_data_files <- read_tsv(args[1], show_col_types = FALSE)
sc_data_files <- read_tsv(args[2], show_col_types = FALSE)
saigeqtl_data_files <- read_tsv(args[3], show_col_types = FALSE)

sc_genes <- sc_data_files |>
  filter(cell_frac == 1) |>
  filter(indiv_frac == 1) |>
  mutate(method = "sc") |>
  rowwise() |>
  mutate(gene_tbl = list(read_tsv(region_file, show_col_types = FALSE))) |>
  unnest(gene_tbl) |>
  ungroup()

pb_genes <- pb_data_files |>
  filter(cell_frac == 1) |>
  filter(indiv_frac == 1) |>
  mutate(method = paste0("pb_", model)) |>
  rowwise() |>
  mutate(gene_tbl = list(read_tsv(region_file, show_col_types = FALSE))) |>
  unnest(gene_tbl) |>
  ungroup()

#saigeqtl_genes <- saigeqtl_data_files |>
#  mutate(method = "saigeqtl") |>
#  rowwise() |>
#  mutate(gene_tbl = list(read_tsv(region_file))) |>
#  rename(pvalue = ACAT_p) |>
#  unnest(gene_tbl) |>
#  ungroup()

all_genes <- bind_rows(
  sc_genes,
  pb_genes
) |>
  mutate(
    bh = p.adjust(pvalue, method = "BH"),
    .by = c(cell_type, method)
  )

sig_genes <- all_genes |>
  filter(bh < 0.05) |>
  select(cell_type, method, feature_id)

upset_df <- sig_genes |>
  distinct(method, feature_id) |>
  group_by(feature_id) |>
  summarise(methods = list(sort(unique(method))), .groups = "drop")

p <- upset_df |>
  ggplot(aes(x = methods)) +
  geom_bar() +
  scale_x_upset() +
  labs(x = NULL, y = "Number of significant genes")

ggsave(
  "concordance-plot.pdf",
  p,
  width = 8,
  height = 6
)

# Additional figure: P-GLMM vs LM UpSet, faceted by cell type.
lm_pglmm_genes <- bind_rows(
  sc_data_files |>
    filter(
      model == "p_glmm_sc",
      cov_spec == "bulk_pca",
      k == "none",
      cell_frac == 1,
      indiv_frac == 1
    ) |>
    mutate(method = "p_glmm_sc") |>
    rowwise() |>
    mutate(gene_tbl = list(read_tsv(region_file, show_col_types = FALSE))) |>
    unnest(gene_tbl) |>
    ungroup(),
  pb_data_files |>
    filter(
      model == "lm",
      int_cov == "none",
      cell_frac == 1,
      indiv_frac == 1
    ) |>
    mutate(method = "lm") |>
    rowwise() |>
    mutate(gene_tbl = list(read_tsv(region_file, show_col_types = FALSE))) |>
    unnest(gene_tbl) |>
    ungroup()
) |>
  mutate(
    bh = p.adjust(pvalue, method = "BH"),
    .by = c(cell_type, method)
  ) |>
  filter(bh < 0.05) |>
  mutate(
    method_label = unname(method_lookup[method]),
    cell_type_label = factor(
      coalesce(unname(cell_type_lookup[cell_type]), cell_type),
      levels = unname(cell_type_lookup[c("CD4_NC", "B_IN", "Plasma")])
    )
  )

upset_lm_pglmm <- lm_pglmm_genes |>
  filter(!is.na(cell_type_label)) |>
  distinct(cell_type_label, method_label, feature_id) |>
  group_by(cell_type_label, feature_id) |>
  summarise(methods = list(sort(unique(method_label))), .groups = "drop")

p_lm_pglmm <- upset_lm_pglmm |>
  ggplot(aes(x = methods)) +
  geom_bar() +
  scale_x_upset() +
  facet_wrap(~cell_type_label, scales = "free_y") +
  labs(
    title = "eGene concordance: Pseudobulk LM vs Single-cell P-GLMM",
    x = NULL,
    y = "Number of significant eGenes"
  ) +
  theme_jp() +
  theme(
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 12)
  )

n_facet <- n_distinct(upset_lm_pglmm$cell_type_label)
ggsave(
  "concordance-lm-p_glmm_sc-plot.pdf",
  p_lm_pglmm,
  width = max(8, 4 * n_facet),
  height = 6,
  limitsize = FALSE
)
