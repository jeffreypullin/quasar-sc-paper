#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(stringr)
  library(data.table)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)

neg_log10_p <- function(p) {
  p <- pmin(pmax(p, 1e-300), 1 - 1e-300)
  -log10(p)
}

format_pvalue <- function(p) {
  out <- rep("NA", length(p))
  ok <- !is.na(p)
  out[ok & p[ok] < 1e-4] <- formatC(p[ok & p[ok] < 1e-4], format = "e", digits = 2)
  out[ok & p[ok] >= 1e-4] <- formatC(p[ok & p[ok] >= 1e-4], format = "g", digits = 3)
  out
}

gene_leads <- function(data, pval_col, require_het = FALSE) {
  keep <- data[!is.na(get(pval_col)) & !is.na(feature_id) & !is.na(snp_id)]
  if (require_het) {
    keep <- keep[!is.na(group_het_pvalue)]
  }
  keep[order(get(pval_col)), .SD[1L], by = feature_id]
}

build_plot_data <- function(pairs, pval_col) {
  pairs |>
    mutate(
      pair_label = paste0(feature_id, "\np = ", format_pvalue(.data[[pval_col]])),
      line_label = snp_id
    ) |>
    select(
      cell_type,
      int_cov,
      k,
      feature_id,
      snp_id,
      pair_label,
      line_label,
      matches("_q[0-9]+_(beta|se)$")
    ) |>
    pivot_longer(
      cols = matches("_q[0-9]+_(beta|se)$"),
      names_to = c("quantile", ".value"),
      names_pattern = "(.*)_(beta|se)$"
    ) |>
    mutate(
      q_num = as.integer(str_extract(quantile, "(?<=_q)\\d+")),
      quantile = factor(paste0("q", q_num), levels = paste0("q", sort(unique(q_num)))),
      pair_label = factor(pair_label, levels = unique(pair_label))
    )
}

make_group_value_plot <- function(plot_data, title) {
  plot_data |>
    ggplot(aes(x = quantile, y = beta, colour = line_label, group = line_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_errorbar(aes(ymin = beta - se, ymax = beta + se), width = 0.15, alpha = 0.8) +
    geom_line(alpha = 0.8) +
    geom_point(size = 1.8) +
    facet_wrap(~pair_label, scales = "free_y") +
    labs(
      x = "Quantile",
      y = "Beta (+/- SE)",
      colour = "Lead SNP",
      title = title
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
}

sc_data_files <- read_tsv(args[1], show_col_types = FALSE) |>
  filter(k != "none")

acat_gws_leads <- list()
acat_plot_leads <- list()
linear_plot_leads <- list()
region_genes <- list()

for (i in seq_len(nrow(sc_data_files))) {
  info <- sc_data_files[i, ]
  variant_data <- fread(info$variant_file, showProgress = FALSE)
  variant_data[, `:=`(
    cell_type = info$cell_type,
    int_cov = info$int_cov,
    k = info$k
  )]

  acat_gws_leads[[i]] <- gene_leads(variant_data, "group_acat_pvalue")
  acat_plot_leads[[i]] <- gene_leads(variant_data, "group_acat_pvalue", require_het = TRUE)
  linear_plot_leads[[i]] <- gene_leads(variant_data, "group_linear_pvalue", require_het = TRUE)
  rm(variant_data)

  region_data <- fread(info$region_file, showProgress = FALSE)
  region_data[, `:=`(
    cell_type = info$cell_type,
    int_cov = info$int_cov,
    k = info$k
  )]
  region_genes[[i]] <- region_data[, .(
    cell_type, int_cov, k, feature_id,
    group_het_acat_pvalue, group_linear_acat_pvalue, group_combined_acat_pvalue
  )]
}

lead_variants <- bind_rows(acat_gws_leads) |>
  group_by(cell_type, int_cov, k, feature_id) |>
  arrange(group_acat_pvalue, .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup() |>
  filter(group_acat_pvalue < 5e-8) |>
  arrange(group_acat_pvalue)

write_tsv(lead_variants, "grouped-sc-output-gws-leads.tsv")

egenes <- bind_rows(region_genes) |>
  mutate(
    het_bh = p.adjust(group_het_acat_pvalue, method = "BH"),
    linear_bh = p.adjust(group_linear_acat_pvalue, method = "BH"),
    acat_bh = p.adjust(group_combined_acat_pvalue, method = "BH"),
    .by = c(cell_type, int_cov, k)
  )

linear_plot_leads <- bind_rows(linear_plot_leads) |>
  group_by(cell_type, int_cov, k, feature_id) |>
  arrange(group_linear_pvalue, .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup() |>
  inner_join(
    egenes |> filter(linear_bh < 0.05),
    by = c("cell_type", "int_cov", "k", "feature_id")
  ) |>
  arrange(group_linear_acat_pvalue)

acat_plot_leads <- bind_rows(acat_plot_leads) |>
  group_by(cell_type, int_cov, k, feature_id) |>
  arrange(group_acat_pvalue, .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup() |>
  inner_join(
    egenes |> filter(acat_bh < 0.05),
    by = c("cell_type", "int_cov", "k", "feature_id")
  ) |>
  arrange(group_combined_acat_pvalue)

for (ct in unique(sc_data_files$cell_type)) {
  linear_plot_data <- build_plot_data(
    linear_plot_leads |> filter(cell_type == ct),
    "group_linear_acat_pvalue"
  )
  acat_plot_data <- build_plot_data(
    acat_plot_leads |> filter(cell_type == ct),
    "group_combined_acat_pvalue"
  )

  n_linear <- n_distinct(linear_plot_data$pair_label)
  n_acat <- n_distinct(acat_plot_data$pair_label)

  ggsave(
    paste0("grouped-sc-output-", ct, "-linear-egenes-plot.pdf"),
    make_group_value_plot(
      linear_plot_data,
      paste0("Significant linear eGenes (BH < 0.05) (", ct, ")")
    ),
    width = 12,
    height = 3.5 * ceiling(n_linear / 3),
    limitsize = FALSE
  )
  ggsave(
    paste0("grouped-sc-output-", ct, "-acat-egenes-plot.pdf"),
    make_group_value_plot(
      acat_plot_data,
      paste0("Significant ACAT eGenes (BH < 0.05) (", ct, ")")
    ),
    width = 12,
    height = 3.5 * ceiling(n_acat / 3),
    limitsize = FALSE
  )
}

scatter_p <- egenes |>
  mutate(
    neg_log10_het = neg_log10_p(group_het_acat_pvalue),
    neg_log10_linear = neg_log10_p(group_linear_acat_pvalue)
  ) |>
  ggplot(aes(x = neg_log10_linear, y = neg_log10_het)) +
  geom_point(alpha = 0.35, size = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.3) +
  facet_wrap(vars(cell_type, int_cov)) +
  labs(
    x = expression(-log[10](linear ~ eGene ~ pvalue)),
    y = expression(-log[10](het ~ eGene ~ pvalue))
  ) +
  theme_jp()

ggsave("grouped-sc-output-het-vs-linear-plot.pdf", scatter_p, width = 8, height = 8)
