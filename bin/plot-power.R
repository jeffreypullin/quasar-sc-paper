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
  filter(n_cells_target < 0) |>
  filter(count_frac == 1) |>
  rowwise() |>
  mutate(test = list(read_tsv(region_file, show_col_types = FALSE))) |>
  unnest(test) |>
  summarise(
    n_sig_gene = sum(p.adjust(pvalue, method = "BH") < 0.05, na.rm = TRUE),
    .by = c(cell_type, cov_spec, cell_frac, model)
  )

gene_prop_data <- pb_data_files |>
  filter(n_cells_target < 0) |>
  filter(count_frac == 1) |>
  rowwise() |>
  mutate(
    gene_data = list(read_tsv(gene_prop_file, show_col_types = FALSE)),
  ) |>
  unnest(gene_data) |>
  distinct(feature_id, cell_type, model, .keep_all = TRUE)

pb_n_sig_gene_data <- pb_data_files |>
  filter(indiv_frac == 1) |>
  filter(n_cells_target < 0) |>
  filter(count_frac == 1) |>
  filter(int_cov == "none") |>
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
  filter(n_cells_target < 0) |>
  filter(count_frac == 1) |>
  rowwise() |>
  mutate(n_sig_variant = list(read_tsv(power_file, show_col_types = FALSE))) |>
  unnest(n_sig_variant) |>
  summarise(n_sig_variant = sum(n_sig_variant), .by = c(cell_type, cell_frac, cov_spec, model))

pb_n_sig_var_data <- pb_data_files |>
  filter(indiv_frac == 1) |>
  filter(n_cells_target < 0) |>
  filter(count_frac == 1) |>
  rowwise() |>
  mutate(n_sig_variant = list(read_tsv(power_file, show_col_types = FALSE))) |>
  unnest(n_sig_variant) |>
  summarise(n_sig_variant = sum(n_sig_variant), .by = c(cell_type, model, cell_frac))

variant_plot_data <- bind_rows(
  pb_n_sig_var_data |>
    filter(cell_frac == 1) |>
    summarise(n_sig_variant = sum(n_sig_variant), .by = c(cell_type, model)) |>
    mutate(type = paste0("pb-", model)) |>
    select(-model),
  sc_n_sig_var_data |>
    filter(cell_frac == 1) |>
    filter(cov_spec == "bulk_pca") |>
    summarise(n_sig_variant = sum(n_sig_variant), .by = c(cell_type, model)) |>
    mutate(type = paste0("sc-", model)),
  saigeqtl_n_sig_var_data |>
    mutate(type = "saigeqtl")
) |>
  filter(cell_type %in% c("Plasma", "B_IN", "CD4_NC")) |>
  filter(type %in% names(method_lookup)) |>
  mutate(
    cell_type = cell_type_lookup[cell_type],
    type = method_lookup[type]
  ) |>
  mutate(
    cell_type = fct_reorder(factor(cell_type), n_sig_variant, .fun = max),
    type = fct_reorder(factor(type), n_sig_variant, .fun = sum)
  )

method_levels <- levels(variant_plot_data$type)

gene_plot_data <- bind_rows(
  pb_n_sig_gene_data |>
    filter(cell_frac == 1) |>
    mutate(type = paste0("pb-", model)),
  sc_n_sig_gene_data |>
    filter(cell_frac == 1) |>
    filter(cov_spec == "bulk_pca") |>
    mutate(type = paste0("sc-", model)),
  saigeqtl_n_sig_gene_data |>
    mutate(type = "saigeqtl")
) |>
  filter(cell_type %in% c("Plasma", "B_IN", "CD4_NC")) |>
  filter(type %in% names(method_lookup)) |>
  mutate(
    cell_type = cell_type_lookup[cell_type],
    type = factor(method_lookup[type], levels = method_levels)
  ) |>
  mutate(cell_type = fct_reorder(factor(cell_type), n_sig_gene, .fun = max))

p_variant <- variant_plot_data |>
  ggplot(aes(cell_type, n_sig_variant, fill = type)) +
  geom_col(position = "dodge2") +
  scale_fill_manual(
    values = method_col_lookup,
    drop = FALSE,
    guide = guide_legend(reverse = TRUE)
  ) +
  labs(
    y = "Number of significant variants",
    x = "Cell type",
    fill = "Method"
  ) +
  coord_flip() +
  theme_jp_vgrid()

p_gene <- gene_plot_data |>
  ggplot(aes(cell_type, n_sig_gene, fill = type)) +
  geom_col(position = "dodge2") +
  scale_fill_manual(values = method_col_lookup, drop = FALSE) +
  labs(
    y = "Number of significant genes",
    x = NULL,
    fill = "Method"
  ) +
  coord_flip() +
  theme_jp_vgrid() +
  guides(fill = "none")

p <- p_variant + p_gene +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    axis.text = element_text(size = 18)
  )

ggsave(
  "power-method-comparison-plot.pdf",
  plot = p,
  width = 16,
  height = 10
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


