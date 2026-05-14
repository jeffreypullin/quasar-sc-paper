#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(data.table)
  library(ggplot2)
  library(forcats)
  library(ggupset)
  library(tidyr)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)

pb_data_files <- read_tsv(args[1], show_col_types = FALSE)
sc_data_files <- read_tsv(args[2], show_col_types = FALSE)
saigeqtl_data_files <- read_tsv(args[3], show_col_types = FALSE)

sc_genes <- sc_data_files |>
  filter(cell_frac == 1) |>
  filter(indiv_frac == 1) |>
  mutate(method = "sc") |>
  rowwise() |>
  mutate(gene_tbl = list(read_tsv(region_file))) |>
  unnest(gene_tbl) |>
  ungroup()

pb_genes <- pb_data_files |>
  filter(cell_frac == 1) |>
  filter(indiv_frac == 1) |>
  mutate(method = paste0("pb_", model)) |>
  rowwise() |>
  mutate(gene_tbl = list(read_tsv(region_file))) |>
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
