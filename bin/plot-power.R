#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(data.table)
  library(ggplot2)
  library(forcats)
  library(patchwork)
  library(qvalue)
  library(purrr)
  library(stringr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

pb_data_files <- read_tsv(args[1], show_col_types = FALSE)
sc_data_files <- read_tsv(args[2], show_col_types = FALSE)
saigeqtl_data_files <- read_tsv(args[3], show_col_types = FALSE)

saigeqtl_n_sig_var_data <- saigeqtl_data_files |>
  filter(variant_file != "NA") |>
  summarise(n_sig_variant = sum(n_sig_variant), .by = cell_type)

saigeqtl_n_sig_gene_data <- saigeqtl_data_files |>
  filter(variant_file != "NA") |>
  rowwise() |>
  mutate(test = list(read_tsv(region_file, show_col_type = FALSE))) |>
  unnest(test) |>
  summarise(
    n_sig_gene = sum(p.adjust(ACAT_p, method = "BH") < 0.05, na.rm = TRUE),
    .by = cell_type
  )

sc_n_sig_gene_data <- sc_data_files |>
  filter(indiv_frac == 1) |>
  rowwise() |>
  mutate(test = list(read_tsv(region_file, show_col_types = FALSE))) |>
  unnest(test) |>
  summarise(
    n_sig_gene = sum(p.adjust(pvalue, method = "BH") < 0.05, na.rm = TRUE),
    .by = c(cell_type, cov_spec, cell_frac)
  )

gene_prop_data <- pb_data_files |>
  rowwise() |>
  mutate(
    gene_data = list(read_tsv(gene_prop_file, show_col_types = FALSE)),
  ) |>
  unnest(gene_data) |>
  distinct(feature_id, cell_type, model, .keep_all = TRUE)

pb_n_sig_gene_data <- pb_data_files |>
  rowwise() |>
  mutate(
    test = list(read_tsv(region_file, show_col_types = FALSE)),
  ) |>
  unnest(test) |>
  select(pvalue, model, cell_type, feature_id, cell_frac) |>
  left_join(
    gene_prop_data |>
      select(model, cell_type, feature_id, pb_mean),
    by = c("cell_type", "feature_id", "model")
  ) |>
  summarise(
    n_sig_gene = sum(p.adjust(pvalue[!is.na(pvalue)], method = "BH") < 0.05, na.rm = TRUE),
    .by = c(cell_type, model, cell_frac)
  )

sc_n_sig_var_data <- sc_data_files |>
  filter(indiv_frac == 1) |>
  rowwise() |>
  mutate(n_sig_variant = list(read_tsv(power_file, show_col_types = FALSE))) |>
  unnest(n_sig_variant) |>
  summarise(n_sig_variant = sum(n_sig_variant), .by = c(cell_type, cell_frac, cov_spec))

pb_n_sig_var_data <- pb_data_files |>
  filter(indiv_frac == 1) |>
  rowwise() |>
  mutate(n_sig_variant = list(read_tsv(power_file, show_col_types = FALSE))) |>
  unnest(n_sig_variant) |>
  summarise(n_sig_variant = sum(n_sig_variant), .by = c(cell_type, model, cell_frac))

p_variant <- bind_rows(
  pb_n_sig_var_data |>
    filter(cell_frac == 1) |>
    summarise(n_sig_variant = sum(n_sig_variant), .by = c(cell_type, model)) |>
    mutate(type = paste0("pb-", model)) |>
    select(-model),
  sc_n_sig_var_data |>
    filter(cell_frac == 1) |>
    filter(cov_spec == "bulk_pca") |>
    summarise(n_sig_variant = sum(n_sig_variant), .by = c(cell_type)) |>
    mutate(type = "sc"),
  saigeqtl_n_sig_var_data |>
    mutate(type = "saigeqtl")
) |>
  filter(cell_type %in% c("Plasma", "B_IN", "CD4_NC")) |>
  mutate(type = case_when(
    type == "sc" ~ "Single-cell quasar",
    type == "pb-nb_glm" ~ "Pseudobulk NB-GLM quasar"
  )) |>
  mutate(cell_type = fct_reorder(factor(cell_type), n_sig_variant)) |>
  ggplot(aes(cell_type, n_sig_variant, fill = type)) +
  geom_col(position = "dodge2") +
  labs(
    y = "Number of signficant variants",
    x = "Cell type"
  ) + 
  scale_fill_manual(values = c(
    "Single-cell quasar" = "#228833", 
    "Pseudobulk NB-GLM quasar" = "#EE6677")
  ) +
  coord_flip() + 
  theme_jp_vgrid()

p_gene <- bind_rows(
  pb_n_sig_gene_data |>
    filter(cell_frac == 1) |>
    mutate(type = paste0("pb-", model)),
  sc_n_sig_gene_data |>
    filter(cell_frac == 1) |>
    filter(cov_spec == "bulk_pca") |>
    mutate(type = "sc"),
  saigeqtl_n_sig_gene_data |>
    mutate(type = "saigeqtl")
) |>
  filter(cell_type %in% c("Plasma", "B_IN", "CD4_NC")) |>
  mutate(cell_type = fct_reorder(factor(cell_type), n_sig_gene)) |>
  ggplot(aes(cell_type, n_sig_gene, fill = type)) +
  geom_col(position = "dodge2") +
  coord_flip() +
  theme_jp_vgrid()

p <- p_variant + p_gene

ggsave(
  "power-method-comparison-plot.pdf",
  plot = p,
  width = 12,
  height = 8
)

#frac_p <- sc_n_sig_var_data |>
#  summarise(
#    n_sig_variant = sum(n_sig_variant),
#    .by = c(cell_type, cov_spec, cell_frac)
#  ) |>
#  filter(cell_type == "CD4_NC") |>
#  filter(cov_spec == "bulk_pca") |>
#  ggplot(aes(cell_frac, n_sig_variant)) +
#  geom_point() +
#  geom_smooth(method = "loess", se = FALSE)

ggsave(
  "power-frac-plot.pdf",
  plot = p,
  width = 12,
  height = 8
)

#non_zero_p <- sc_n_sig_var_data |>
#  filter(cell_frac == 1) |>
#  filter(cov_spec == "bulk_pca") |>
#  summarise(
#    n_sig_variant = sum(n_sig_variant),
#    .by = c(cell_type, nz_bin)
#  ) |>
#  ggplot(aes(nz_bin, n_sig_variant)) +
#  geom_col(position = "dodge2") +
#  facet_wrap(~ cell_type)

ggsave(
  "power-sc-non-zero-frac-plot.pdf",
  plot = p,
  width = 12,
  height = 8
)
