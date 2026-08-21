#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(purrr)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)

leads_file <- args[1]
sc_logcounts_file <- args[2]
genotype_file <- args[3]
cell_type <- args[4]

leads <- read_tsv(leads_file, show_col_types = FALSE) |>
  filter(cell_type == !!cell_type) |>
  arrange(region_bh, feature_id) |>
  slice_head(n = 10)

if (nrow(leads) == 0) {
  stop(sprintf("No unique eGenes found for cell type %s", cell_type))
}

sc_logcounts <- read_tsv(sc_logcounts_file, show_col_types = FALSE)
genotypes <- read_tsv(genotype_file, show_col_types = FALSE)

plot_one_pair <- function(gene_id, variant_id, region_bh) {
  plot_df <- sc_logcounts |>
    mutate(log_expr = .data[[gene_id]]) |>
    select(sample_id, cell_id, log_expr) |>
    left_join(
      genotypes |>
        mutate(dosage = .data[[variant_id]]) |>
        select(sample_id, dosage),
      by = "sample_id"
    ) |>
    filter(!is.na(log_expr), !is.na(dosage)) |>
    summarise(
      mean_expr = mean(log_expr, na.rm = TRUE),
      .by = c(sample_id, dosage)
    ) |>
    filter(!is.na(mean_expr), !is.na(dosage))

  ggplot(plot_df, aes(x = dosage, y = mean_expr)) +
    geom_point(
      alpha = 0.35,
      size = 0.8,
      colour = "grey40",
      position = position_jitter(width = 0.06, height = 0)
    ) +
    geom_boxplot(
      aes(group = dosage),
      width = 0.5,
      fill = "white",
      colour = "grey20",
      linewidth = 0.4,
      alpha = 0.9,
      outlier.shape = NA
    ) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, colour = "black") +
    scale_x_continuous(breaks = sort(unique(plot_df$dosage))) +
    labs(
      title = sprintf("%s\n%s", gene_id, variant_id),
      subtitle = sprintf("region BH = %.2e", region_bh),
      x = "Genotype dosage",
      y = "Mean log1p expression"
    ) +
    theme_jp() +
    theme(
      plot.title = element_text(size = 8, lineheight = 1.05),
      plot.subtitle = element_text(size = 7, margin = margin(0, 0, 4, 0)),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 6)
    )
}

plots <- pmap(
  list(leads$feature_id, leads$snp_id, leads$region_bh),
  plot_one_pair
)

cell_label <- coalesce(unname(cell_type_lookup[cell_type]), cell_type)
p <- wrap_plots(plots, ncol = 2) +
  plot_annotation(
    title = sprintf("Top unique P-GLMM eGenes: %s", cell_label),
    theme = theme(plot.title = element_text(size = 11, family = "Helvetica"))
  )

n_plots <- length(plots)
ggsave(
  sprintf("unique-sc-egene-figures-%s.pdf", cell_type),
  p,
  width = 12,
  height = max(8, 3.5 * ceiling(n_plots / 2)),
  limitsize = FALSE
)
