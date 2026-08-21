#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(forcats)
  library(patchwork)
  library(tidyr)
  library(purrr)
})

args <- commandArgs(trailingOnly = TRUE)

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

n_cells_specs <- tribble(
  ~cell_type, ~n_cells_target,
  "CD4_NC", 300L,
  "B_IN", 100L
)

pb_data_files <- read_tsv(args[1], show_col_types = FALSE) |>
  inner_join(n_cells_specs, by = c("cell_type", "n_cells_target"))

sc_data_files <- read_tsv(args[2], show_col_types = FALSE) |>
  inner_join(n_cells_specs, by = c("cell_type", "n_cells_target")) |>
  filter(cov_spec == "bulk_pca")

sc_n_sig_gene_data <- sc_data_files |>
  rowwise() |>
  mutate(test = list(read_tsv(region_file, show_col_types = FALSE))) |>
  unnest(test) |>
  summarise(
    n_sig_gene = sum(p.adjust(pvalue, method = "BH") < 0.05, na.rm = TRUE),
    .by = c(cell_type, model)
  )

pb_n_sig_gene_data <- pb_data_files |>
  filter(int_cov == "none") |>
  rowwise() |>
  mutate(test = list(read_tsv(region_file, show_col_types = FALSE))) |>
  unnest(test) |>
  summarise(
    n_sig_gene = sum(p.adjust(pvalue[!is.na(pvalue)], method = "BH") < 0.05, na.rm = TRUE),
    .by = c(cell_type, model)
  )

sc_n_sig_var_data <- sc_data_files |>
  rowwise() |>
  mutate(n_sig_variant = list(read_tsv(power_file, show_col_types = FALSE))) |>
  unnest(n_sig_variant) |>
  summarise(n_sig_variant = sum(n_sig_variant), .by = c(cell_type, model))

pb_n_sig_var_data <- pb_data_files |>
  rowwise() |>
  mutate(n_sig_variant = list(read_tsv(power_file, show_col_types = FALSE))) |>
  unnest(n_sig_variant) |>
  summarise(n_sig_variant = sum(n_sig_variant), .by = c(cell_type, model))

variant_plot_data <- bind_rows(
  pb_n_sig_var_data |> mutate(type = paste0("pb-", model)),
  sc_n_sig_var_data |> mutate(type = paste0("sc-", model))
) |>
  filter(type %in% names(method_lookup)) |>
  filter(type != "saigeqtl") |>
  mutate(type = method_lookup[type]) |>
  left_join(n_cells_specs, by = "cell_type") |>
  mutate(
    panel_label = paste0(
      cell_type_lookup[cell_type],
      ", donors with >=",
      n_cells_target,
      " cells downsampled to ",
      n_cells_target
    ),
    type = fct_reorder(factor(type), n_sig_variant, .fun = sum)
  )

method_levels <- levels(variant_plot_data$type)

gene_plot_data <- bind_rows(
  pb_n_sig_gene_data |> mutate(type = paste0("pb-", model)),
  sc_n_sig_gene_data |> mutate(type = paste0("sc-", model))
) |>
  filter(type %in% names(method_lookup)) |>
  filter(type != "saigeqtl") |>
  mutate(type = factor(method_lookup[type], levels = method_levels)) |>
  left_join(n_cells_specs, by = "cell_type") |>
  mutate(
    panel_label = paste0(
      cell_type_lookup[cell_type],
      ", donors with >=",
      n_cells_target,
      " cells downsampled to ",
      n_cells_target
    )
  )

panel_levels <- unique(variant_plot_data$panel_label)
variant_plot_data <- variant_plot_data |>
  mutate(panel_label = factor(panel_label, levels = panel_levels))
gene_plot_data <- gene_plot_data |>
  mutate(panel_label = factor(panel_label, levels = panel_levels))

p_variant <- variant_plot_data |>
  ggplot(aes(type, n_sig_variant, fill = type)) +
  geom_col() +
  scale_fill_manual(
    values = method_col_lookup,
    drop = FALSE,
    guide = guide_legend(reverse = TRUE)
  ) +
  facet_wrap(~panel_label, ncol = 1, scales = "free_x") +
  labs(
    y = "Number of significant variants",
    x = NULL,
    fill = "Method"
  ) +
  coord_flip() +
  theme_jp_vgrid()

p_gene <- gene_plot_data |>
  ggplot(aes(type, n_sig_gene, fill = type)) +
  geom_col() +
  scale_fill_manual(values = method_col_lookup, drop = FALSE) +
  facet_wrap(~panel_label, ncol = 1, scales = "free_x") +
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
  "power-n-cells-filter-plot.pdf",
  plot = p,
  width = 16,
  height = 12
)
