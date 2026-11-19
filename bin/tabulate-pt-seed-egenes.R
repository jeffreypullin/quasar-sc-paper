#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(tidyr)
  library(purrr)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)
manifest <- read_tsv(args[1], show_col_types = FALSE)
gene_map <- read_tsv(args[2], show_col_types = FALSE) |>
  distinct(feature_id, gene_symbol)

genes <- manifest |>
  mutate(analysis = if_else(k == "none", "interaction", "grouped")) |>
  mutate(gene_tbl = map2(region_file, analysis, function(path, analysis) {
    col <- if (analysis == "interaction") {
      "int_pseudotime_acat_pvalue"
    } else {
      "group_linear_acat_pvalue"
    }
    read_tsv(path, show_col_types = FALSE) |>
      select(feature_id, pvalue = all_of(col))
  })) |>
  select(pt_seed, analysis, gene_tbl) |>
  unnest(gene_tbl) |>
  filter(!is.na(feature_id), is.finite(pvalue), pvalue > 0, pvalue <= 1) |>
  slice_min(pvalue, n = 1, with_ties = FALSE, by = c(pt_seed, analysis, feature_id)) |>
  mutate(bh = p.adjust(pvalue, method = "BH"), .by = c(pt_seed, analysis))

sig <- genes |>
  filter(bh < 0.05) |>
  count(analysis, feature_id, name = "n_seeds")

seeds <- sort(unique(genes$pt_seed))
bh_cols <- paste0("bh_", seeds)
p_cols <- paste0("pvalue_", seeds)

table <- genes |>
  inner_join(sig, by = c("analysis", "feature_id")) |>
  pivot_wider(
    id_cols = c(analysis, feature_id, n_seeds),
    names_from = pt_seed,
    values_from = c(bh, pvalue),
    names_glue = "{.value}_{pt_seed}"
  ) |>
  left_join(gene_map, by = "feature_id") |>
  select(analysis, feature_id, gene_symbol, n_seeds, all_of(bh_cols), all_of(p_cols)) |>
  arrange(analysis, desc(n_seeds), feature_id)

write_tsv(table, "pseudotime-seed-egenes.tsv")

plot_df <- genes |>
  inner_join(sig, by = c("analysis", "feature_id")) |>
  left_join(gene_map, by = "feature_id") |>
  mutate(
    seed = factor(pt_seed, levels = seeds),
    neglog = pmin(-log10(bh), 8),
    significant = bh < 0.05
  )

seed_panel <- function(analysis_name, title) {
  data <- plot_df |>
    filter(analysis == analysis_name)
  gene_order <- data |>
    distinct(feature_id, gene_symbol, n_seeds) |>
    arrange(desc(n_seeds), feature_id) |>
    pull(gene_symbol)
  data <- data |>
    mutate(gene_symbol = factor(gene_symbol, levels = gene_order))
  data |>
    ggplot(aes(seed, gene_symbol, fill = neglog)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_point(
      data = filter(data, significant),
      colour = "#CC3311",
      size = 2.4
    ) +
    scale_y_discrete(limits = rev(gene_order)) +
    scale_fill_gradientn(
      colours = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
      limits = c(0, 8),
      breaks = c(0, 2, 4, 6, 8),
      name = expression(-log[10] * " BH, capped at 8")
    ) +
    labs(title = title, x = "PHATE seed", y = NULL) +
    theme_jp() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_text(family = "Helvetica", size = 14, colour = "#222222"),
      legend.position = "bottom"
    )
}

seed_plot <- seed_panel("grouped", "Grouped") /
  seed_panel("interaction", "Interaction") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  "pseudotime-seed-egenes.pdf",
  seed_plot,
  width = 7.5,
  height = 12
)
