#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

pc_gwas_data_files <- read_tsv(args[1], show_col_types = FALSE)

pc_col_types <- cols(
  feature_id = col_character(),
  snp_id = col_character(),
  chrom = col_character(),
  pos = col_double(),
  alt = col_character(),
  ref = col_character(),
  maf = col_double(),
  beta = col_double(),
  se = col_double(),
  pvalue = col_double()
)

pc_gwas_data <- pc_gwas_data_files |>
  rowwise() |>
  mutate(var_data = list(read_tsv(
    sig_variant_file,
    show_col_types = FALSE,
    col_types = pc_col_types
  ))) |>
  ungroup() |>
  unnest(var_data) |>
  select(cell_type, feature_id, snp_id, chrom, pos, maf, beta, pvalue)

plot_pc_gwas_manhattan <- function(data, ct, out_file) {
  cell_data <- data |> filter(cell_type == ct)

  chrom_levels <- unique(cell_data$chrom)
  chrom_levels <- chrom_levels[order(
    suppressWarnings(as.integer(gsub("[^0-9]", "", chrom_levels))),
    chrom_levels
  )]

  pc_levels <- unique(cell_data$feature_id)
  pc_levels <- pc_levels[order(
    suppressWarnings(as.integer(gsub("[^0-9]", "", pc_levels))),
    pc_levels
  )]

  manhattan_data <- cell_data |>
    filter(!is.na(pvalue) & pvalue > 0) |>
    mutate(
      chrom = factor(chrom, levels = chrom_levels),
      pc_id = factor(feature_id, levels = pc_levels),
      neg_log10_p = -log10(pvalue)
    )

  manhattan_plot <- manhattan_data |>
    ggplot(aes(x = pos, y = neg_log10_p, colour = chrom)) +
    geom_point(size = 0.8, alpha = 0.7) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed") +
    facet_grid(pc_id ~ chrom, scales = "free", space = "free_x", switch = "x") +
    labs(
      x = "Chromosome",
      y = expression(-log[10](p)),
      title = paste0(ct, " single-cell PC GWAS by principal component")
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing.x = unit(0.05, "lines"),
      strip.background = element_blank()
    )

  n_pcs <- n_distinct(manhattan_data$pc_id)
  ggsave(
    out_file,
    manhattan_plot,
    width = 10,
    height = max(6, 2 * n_pcs),
    limitsize = FALSE
  )
}

plot_pc_gwas_manhattan(pc_gwas_data, "B_IN", "pc-gwas-B_IN-plot.pdf")
plot_pc_gwas_manhattan(pc_gwas_data, "CD4_NC", "pc-gwas-CD4_NC-plot.pdf")

gws_pval_threshold <- 5e-8

lead_variants <- pc_gwas_data |>
  filter(!is.na(pvalue), pvalue > 0) |>
  group_by(cell_type, feature_id, chrom) |>
  slice_min(order_by = pvalue, n = 1, with_ties = FALSE) |>
  ungroup() |>
  filter(pvalue < gws_pval_threshold) |>
  arrange(cell_type, feature_id, chrom, pvalue)

write_tsv(lead_variants, "pc-gwas-gws-leads.tsv")

lead_variants |>
  print(n = Inf)
