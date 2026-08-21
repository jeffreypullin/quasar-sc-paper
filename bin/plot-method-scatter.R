#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)

pb_data_files <- read_tsv(args[1], show_col_types = FALSE)
sc_data_files <- read_tsv(args[2], show_col_types = FALSE)

cell_types <- c("Plasma", "B_IN", "CD4_NC")

pb_manifest <- pb_data_files |>
  filter(
    model == "lm",
    int_cov == "none",
    cell_frac == 1,
    indiv_frac == 1,
    cell_type %in% cell_types
  )

sc_manifest <- sc_data_files |>
  filter(
    model == "p_glmm_sc",
    cov_spec == "bulk_pca",
    k == "none",
    cell_frac == 1,
    indiv_frac == 1,
    cell_type %in% cell_types
  )

read_region_bh <- function(manifest, method) {
  manifest |>
    mutate(method = method) |>
    rowwise() |>
    mutate(gene_tbl = list(read_tsv(region_file, show_col_types = FALSE))) |>
    unnest(gene_tbl) |>
    ungroup() |>
    mutate(
      region_bh = p.adjust(pvalue, method = "BH"),
      .by = c(cell_type, method)
    ) |>
    select(cell_type, feature_id, region_bh) |>
    distinct()
}

plot_data <- inner_join(
  read_region_bh(pb_manifest, "lm") |> rename(region_bh_lm = region_bh),
  read_region_bh(sc_manifest, "p_glmm_sc") |> rename(region_bh_sc = region_bh),
  by = c("cell_type", "feature_id")
) |>
  mutate(
    neg_log10_bh_lm = -log10(pmax(region_bh_lm, .Machine$double.xmin)),
    neg_log10_bh_sc = -log10(pmax(region_bh_sc, .Machine$double.xmin)),
    delta_neg_log10_bh = neg_log10_bh_lm - neg_log10_bh_sc,
    cell_type_label = factor(
      cell_type_lookup[cell_type],
      levels = unname(cell_type_lookup[cell_types])
    )
  )

p <- ggplot(plot_data, aes(neg_log10_bh_lm, neg_log10_bh_sc)) +
  geom_point(alpha = 0.25, size = 0.6, colour = "grey30") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.4) +
  facet_wrap(~cell_type_label, scales = "free") +
  labs(
    x = expression(-log[10] * " BH (Pseudobulk LM)"),
    y = expression(-log[10] * " BH (Single-cell P-GLMM)")
  ) +
  theme_jp()

ggsave("method-scatter-plot.pdf", p, width = 12, height = 4.5)
