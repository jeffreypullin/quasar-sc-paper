#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(forcats)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

sc_data_files <- read_tsv(args[1], show_col_types = FALSE) |>
  filter(k != "none")

variant_data <- sc_data_files |>
  rowwise() |>
  mutate(var_tbl = list(read_tsv(variant_file, show_col_types = FALSE))) |>
  ungroup() |>
  select(cell_type, int_cov, k, var_tbl) |>
  unnest(var_tbl)

top_pairs <- variant_data |>
  filter(!is.na(group_het_pvalue), !is.na(feature_id), !is.na(snp_id)) |>
  filter(chrom != "6") |>
  group_by(feature_id) |>
  arrange(group_acat_pvalue, .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup() |>
  arrange(group_acat_pvalue) |>
  slice_head(n = 49) |>
  mutate(
    pair_label = paste0(feature_id, "\n", snp_id),
    pair_label = fct_reorder(pair_label, group_acat_pvalue, .desc = TRUE)
  )

plot_data <- top_pairs |>
  select(
    cell_type,
    int_cov,
    k,
    feature_id,
    snp_id,
    pair_label,
    group_het_pvalue,
    group_linear_pvalue,
    group_acat_pvalue,
    matches("_q[0-9]+_(beta|se)$")
  ) |>
  pivot_longer(
    cols = matches("_q[0-9]+_(beta|se)$"),
    names_to = c("quantile", ".value"),
    names_pattern = "(.*)_(beta|se)$"
  ) |>
  filter(!is.na(beta), !is.na(se)) |>
  mutate(
    q_num = as.integer(str_extract(quantile, "(?<=_q)\\d+")),
    quantile = paste0("q", q_num)
  ) |>
  mutate(
    quantile = factor(quantile, levels = paste0("q", sort(unique(q_num))))
  )

p <- plot_data |>
  ggplot(aes(x = quantile, y = beta, colour = cell_type, group = cell_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_errorbar(aes(ymin = beta - se, ymax = beta + se), width = 0.15, alpha = 0.8) +
  geom_line(alpha = 0.8) +
  geom_point(size = 1.8) +
  facet_wrap(~pair_label, ncol = 7, scales = "free_y") +
  labs(
    x = "Quantile",
    y = "Beta (+/- SE)",
    colour = "Cell type",
    title = "Top 9 grouped gene/SNP pairs by group_acat_pvalue"
  ) +
  theme_bw()

ggsave(
  "grouped-sc-output-top9-plot.pdf",
  p,
  width = 11,
  height = 8
)
