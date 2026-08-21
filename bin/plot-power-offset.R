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
  "Per-cell depth (native)" = "#882255",
  "Donor-flattened depth"   = "#44AA99",
  "Constant depth"          = "#DDCC77"
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
  labs(y = "Number of significant variants", x = NULL, fill = "Offset") +
  coord_flip() +
  theme_jp_vgrid() +
  guides(fill = "none")

p_gene <- gene_plot_data |>
  ggplot(aes(type, n_sig_gene, fill = type)) +
  geom_col() +
  scale_fill_manual(values = offset_col_lookup, drop = FALSE) +
  facet_wrap(~cell_type_lookup[cell_type], ncol = 1, scales = "free_x") +
  labs(y = "Number of significant genes", x = NULL, fill = "Offset") +
  coord_flip() +
  theme_jp_vgrid() +
  guides(fill = "none")

ggsave(
  "power-offset-plot.pdf",
  plot = p_variant + p_gene & theme(axis.text = element_text(size = 18)),
  width = 16,
  height = 8
)

# Paired per-gene comparison against the native per-cell offset. Under Poisson
# sufficiency the donor-flattened arm should be indistinguishable from it.
paired_data <- gene_data |>
  pivot_wider(names_from = offset_spec, values_from = pvalue) |>
  pivot_longer(
    any_of(names(comparison_lookup)),
    names_to = "comparison",
    values_to = "alt_pvalue"
  ) |>
  filter(!is.na(percell), !is.na(alt_pvalue)) |>
  mutate(
    comparison_label = factor(
      comparison_lookup[comparison],
      levels = unname(comparison_lookup)
    ),
    percell_score = -log10(percell),
    alt_score = -log10(alt_pvalue),
    abs_diff = abs(alt_score - percell_score)
  )

p_scatter <- paired_data |>
  ggplot(aes(percell_score, alt_score)) +
  geom_abline(slope = 1, intercept = 0, colour = "#cbcbcb", linewidth = 1) +
  geom_point(alpha = 0.3, size = 1.5) +
  facet_grid(cell_type_lookup[cell_type] ~ comparison_label) +
  labs(
    x = expression(-log[10] ~ "p, per-cell offset"),
    y = expression(-log[10] ~ "p, alternative offset")
  ) +
  theme_jp()

ggsave("power-offset-scatter-plot.pdf", plot = p_scatter, width = 14, height = 8)

summary_data <- paired_data |>
  summarise(
    n_genes = n(),
    max_abs_diff_neglog10_p = max(abs_diff),
    median_abs_diff_neglog10_p = median(abs_diff),
    n_genes_abs_diff_gt_1e6 = sum(abs_diff > 1e-6),
    .by = c(cell_type, comparison)
  ) |>
  left_join(
    gene_plot_data |> select(cell_type, comparison = offset_spec, n_sig_gene),
    by = c("cell_type", "comparison")
  ) |>
  left_join(
    gene_plot_data |>
      filter(offset_spec == "percell") |>
      select(cell_type, n_sig_gene_percell = n_sig_gene),
    by = "cell_type"
  ) |>
  left_join(
    variant_plot_data |> select(cell_type, comparison = offset_spec, n_sig_variant),
    by = c("cell_type", "comparison")
  ) |>
  left_join(
    variant_plot_data |>
      filter(offset_spec == "percell") |>
      select(cell_type, n_sig_variant_percell = n_sig_variant),
    by = "cell_type"
  )

write_tsv(summary_data, "power-offset-summary.tsv")
