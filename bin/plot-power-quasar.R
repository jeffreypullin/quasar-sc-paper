#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(patchwork)
  library(purrr)
  library(readr)
  library(tidyr)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)
pb_data_files <- read_tsv(args[1], show_col_types = FALSE)
sc_data_files <- read_tsv(args[2], show_col_types = FALSE)

is_full_data <- function(data) {
  data |>
    filter(
      as.numeric(cell_frac) == 1,
      as.numeric(indiv_frac) == 1,
      as.numeric(n_cells_target) < 0,
      as.numeric(count_frac) == 1,
      int_cov == "none"
    )
}

count_sig_genes <- function(manifest, extra_by = character()) {
  manifest |>
    mutate(gene_tbl = map(region_file, read_tsv, show_col_types = FALSE)) |>
    unnest(gene_tbl) |>
    summarise(
      n_sig_gene = sum(p.adjust(pvalue, method = "BH") < 0.05, na.rm = TRUE),
      .by = all_of(c("cell_type", extra_by))
    )
}

count_sig_variants <- function(manifest, extra_by = character()) {
  manifest |>
    mutate(power_tbl = map(power_file, read_tsv, show_col_types = FALSE)) |>
    unnest(power_tbl) |>
    summarise(
      n_sig_variant = sum(n_sig_variant),
      .by = all_of(c("cell_type", extra_by))
    )
}

pb_manifest <- pb_data_files |>
  is_full_data() |>
  filter(model == "lm")

sc_manifest <- sc_data_files |>
  is_full_data() |>
  filter(cov_spec == "bulk_pca", as.character(k) == "none")

cell_types <- intersect(unique(pb_manifest$cell_type), unique(sc_manifest$cell_type))
pb_manifest <- pb_manifest |> filter(cell_type %in% cell_types)
sc_manifest <- sc_manifest |> filter(cell_type %in% cell_types)

variant_plot_data <- bind_rows(
  count_sig_variants(pb_manifest) |> mutate(type = "pb-lm"),
  count_sig_variants(sc_manifest, "model") |> mutate(type = paste0("sc-", model))
) |>
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
  count_sig_genes(pb_manifest) |> mutate(type = "pb-lm"),
  count_sig_genes(sc_manifest, "model") |> mutate(type = paste0("sc-", model))
) |>
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

ggsave(
  "power-quasar-comparison-plot.pdf",
  plot = p_variant + p_gene +
    plot_layout(guides = "collect") &
    theme(
      legend.position = "right",
      axis.text = element_text(size = 18)
    ),
  width = 16,
  height = 10
)
