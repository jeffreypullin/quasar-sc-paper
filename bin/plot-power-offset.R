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

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)

offset_lookup <- c(
  "percell"   = "Per-cell depth (native)",
  "donorflat" = "Donor-flattened depth",
  "constant"  = "Constant depth"
)

offset_col_lookup <- c(
  "Per-cell depth (native)" = "#0077BB",
  "Donor-flattened depth"   = "#EE7733",
  "Constant depth"          = "#009988"
)

comparison_lookup <- c(
  "donorflat" = "Donor-flattened vs per-cell",
  "constant"  = "Constant vs per-cell"
)

data_files <- read_tsv(args[1], show_col_types = FALSE)

gene_data <- data_files |>
  rowwise() |>
  mutate(gene_tbl = list(read_tsv(region_file, show_col_types = FALSE))) |>
  ungroup() |>
  unnest(gene_tbl) |>
  select(cell_type, offset_spec, chr, feature_id, pvalue)

variant_plot_data <- data_files |>
  rowwise() |>
  mutate(power_tbl = list(read_tsv(power_file, show_col_types = FALSE))) |>
  ungroup() |>
  unnest(power_tbl) |>
  summarise(n_sig_variant = sum(n_sig_variant), .by = c(cell_type, offset_spec)) |>
  mutate(type = factor(offset_lookup[offset_spec], levels = unname(offset_lookup)))

gene_plot_data <- gene_data |>
  summarise(
    n_sig_gene = sum(p.adjust(pvalue, method = "BH") < 0.05, na.rm = TRUE),
    .by = c(cell_type, offset_spec)
  ) |>
  mutate(type = factor(offset_lookup[offset_spec], levels = unname(offset_lookup)))

p_variant <- variant_plot_data |>
  ggplot(aes(type, n_sig_variant, fill = type)) +
  geom_col() +
  scale_fill_manual(values = offset_col_lookup, drop = FALSE) +
  facet_wrap(~cell_type_lookup[cell_type], ncol = 1, scales = "free_x") +
  labs(y = "Number of significant variants", x = "Offset", fill = "Offset") +
  coord_flip() +
  theme_jp_vgrid() +
  guides(fill = "none")

p_gene <- gene_plot_data |>
  ggplot(aes(type, n_sig_gene, fill = type)) +
  geom_col() +
  scale_fill_manual(values = offset_col_lookup, drop = FALSE) +
  facet_wrap(~cell_type_lookup[cell_type], ncol = 1, scales = "free_x") +
  labs(y = "Number of significant genes", x = "Offset", fill = "Offset") +
  coord_flip() +
  theme_jp_vgrid() +
  guides(fill = "none")

# Paired per-gene comparison against the native per-cell offset.
paired_data <- gene_data |>
  pivot_wider(names_from = offset_spec, values_from = pvalue) |>
  pivot_longer(
    any_of(names(comparison_lookup)),
    names_to = "comparison",
    values_to = "alt_pvalue"
  ) |>
  filter(!is.na(percell), !is.na(alt_pvalue)) |>
  mutate(
    percell_score = -log10(percell),
    alt_score = -log10(alt_pvalue)
  )

p_scatter <- names(comparison_lookup) |>
  map(\(spec) {
    paired_data |>
      filter(comparison == spec) |>
      ggplot(aes(percell_score, alt_score)) +
      geom_abline(slope = 1, intercept = 0, colour = "#cbcbcb", linewidth = 1) +
      geom_point(alpha = 0.3, size = 1.5) +
      facet_wrap(~cell_type_lookup[cell_type], ncol = 1) +
      labs(
        x = bquote(-log[10] ~ "p," ~ .(offset_lookup[["percell"]])),
        y = bquote(-log[10] ~ "p," ~ .(offset_lookup[[spec]]))
      ) +
      theme_jp()
  }) |>
  wrap_plots(nrow = 1)

figure <- wrap_elements(full = p_variant + p_gene) /
  wrap_elements(full = p_scatter) +
  plot_annotation(tag_levels = "a")

ggsave("power-offset-plot.pdf", figure, width = 16, height = 14)
