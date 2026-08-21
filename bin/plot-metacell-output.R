#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

sc_data_files <- read_tsv(args[1], show_col_types = FALSE) |>
  filter(k == "seacells") |>
  filter(cell_type == "T_all")

seacells_files <- read_tsv(args[2], show_col_types = FALSE)

# Per-metacell UMAP coordinates and sizes for T_all.
metacell_info <- seacells_files |>
  filter(cell_type == "T_all") |>
  slice_head(n = 1) |>
  pull(info_file) |>
  read_tsv(show_col_types = FALSE) |>
  mutate(metacell = as.character(metacell))

metacell_ids <- metacell_info$metacell

# Full per-variant table, including the per-metacell beta/se columns.
variant_data <- sc_data_files |>
  rowwise() |>
  mutate(var_tbl = list(read_tsv(variant_file, show_col_types = FALSE))) |>
  ungroup() |>
  select(var_tbl) |>
  unnest(var_tbl)

plot_feature_id <- "ENSG00000163599"

# Lead variant (min group_het_pvalue) for the gene of interest.
lead_variants <- variant_data |>
  filter(feature_id == plot_feature_id) |>
  filter(!is.na(group_het_pvalue), !is.na(feature_id), !is.na(snp_id)) |>
  slice_min(group_het_pvalue, n = 1, with_ties = FALSE)

lead_variants |>
  select(feature_id, snp_id, group_het_pvalue) |>
  print()

# Reshape per-metacell beta/se into long form and compute a signed z-score
# (beta / se) describing the strength and direction of the effect per metacell.
plot_data <- lead_variants |>
  select(feature_id, snp_id, group_het_pvalue, matches("_(beta|se)$")) |>
  pivot_longer(
    cols = matches("_(beta|se)$"),
    names_to = c("metacell", ".value"),
    names_pattern = "(.*)_(beta|se)$"
  ) |>
  filter(metacell %in% metacell_ids, !is.na(beta), !is.na(se), se > 0) |>
  mutate(z = beta / se) |>
  left_join(metacell_info, by = "metacell")

z_limit <- max(abs(plot_data$z), na.rm = TRUE)

p <- plot_data |>
  ggplot(aes(x = umap_1, y = umap_2, size = n_cells, colour = z)) +
  geom_point(alpha = 0.85) +
  scale_colour_gradient2(
    low = "#2166ac",
    mid = "grey90",
    high = "#b2182b",
    midpoint = 0,
    limits = c(-z_limit, z_limit)
  ) +
  scale_size_continuous(range = c(1, 8)) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
    size = "Cells per\nmetacell",
    colour = "Effect z-score\n(beta / se)",
    title = paste0("T_all metacell effects for ", plot_feature_id, " (", lead_variants$snp_id, ")")
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 7)
  )

ggsave(
  "metacell-output-sizes-plot.pdf",
  p,
  width = 8,
  height = 6
)
