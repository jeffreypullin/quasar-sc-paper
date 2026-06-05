#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

csaqtl_data_files <- read_tsv(args[1], show_col_types = FALSE)

csaqtl_col_types <- cols(
  feature_id = col_character(),
  snp_id = col_character(),
  chrom = col_character(),
  pos = col_double(),
  alt = col_character(),
  ref = col_character(),
  maf = col_double(),
  beta = col_double(),
  se = col_double(),
  pvalue = col_double(),
  glm_converged = col_double(),
  phi = col_double(),
  phi_converged = col_double()
)

csaqtl_data <- csaqtl_data_files |>
  rowwise() |>
  mutate(var_data = list(read_tsv(
    sig_variant_file,
    show_col_types = FALSE,
    col_types = csaqtl_col_types
  ))) |>
  ungroup() |>
  unnest(var_data) |>
  select(cell_type, feature_id, geno_chr, snp_id, chrom, pos, maf, beta, pvalue)

csaqtl_data |>
  filter(cell_type == "B_all") |>
  filter(pvalue < 5e-8) |>
  group_by(cell_type, chrom) |>
  slice_min(order_by = pvalue, n = 1, with_ties = FALSE) |>
  arrange(cell_type, chrom) |>
  print(n = Inf) |>
  select(-feature_id) |>
  arrange(pvalue)

2 + "djsl"

cd4_t_all_data <- csaqtl_data |>
  filter(cell_type == "B_all")

chrom_levels <- unique(cd4_t_all_data$chrom)
chrom_levels <- chrom_levels[order(
  suppressWarnings(as.integer(gsub("[^0-9]", "", chrom_levels))),
  chrom_levels
)]

seacell_levels <- unique(cd4_t_all_data$feature_id)
seacell_levels <- seacell_levels[order(
  suppressWarnings(as.integer(gsub("[^0-9]", "", seacell_levels))),
  seacell_levels
)]

manhattan_data <- cd4_t_all_data |>
  filter(!is.na(pvalue) & pvalue > 0) |>
  mutate(
    chrom = factor(chrom, levels = chrom_levels),
    seacell_id = factor(feature_id, levels = seacell_levels),
    neg_log10_p = -log10(pvalue)
  )

manhattan_plot <- manhattan_data |>
  ggplot(aes(x = pos, y = neg_log10_p, colour = chrom)) +
  geom_point(size = 0.8, alpha = 0.7) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed") +
  facet_grid(seacell_id ~ chrom, scales = "free_x", space = "free_x", switch = "x") +
  labs(
    x = "Chromosome",
    y = expression(-log[10](p)),
    title = "CD4_T_all cell-state abundance QTL (csaQTL) by SEACell"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing.x = unit(0.05, "lines"),
    strip.background = element_blank()
  )

n_seacells <- n_distinct(manhattan_data$seacell_id)
ggsave(
  "csaqtl-output-plot.pdf",
  manhattan_plot,
  width = 10,
  height = max(6, 2 * n_seacells),
  limitsize = FALSE
)

csaqtl_data |>
  arrange(pvalue) |>
  print(n = 30)
