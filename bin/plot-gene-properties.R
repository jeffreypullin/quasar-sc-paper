#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

gene_prop_files <- read_tsv(args[1], show_col_types = FALSE)

gene_prop_data <- gene_prop_files |>
  rowwise() |>
  mutate(gene_props = list(read_tsv(properties_file))) |>
  unnest(gene_props)

p <- gene_prop_data |>
  ggplot(aes(sc_non_zero_frac, pb_non_zero_frac)) +
  geom_point() +
  geom_vline(xintercept = 0.01, linetype = "dashed", color = "red") +
  facet_wrap(~cell_type)

ggsave(
  "sc-non-zero-frac-vs-pb-non-zero-frac-plot.pdf",
  plot = p,
  width = 12,
  height = 8
)

p <- gene_prop_data |>
  ggplot(aes(sc_non_zero_frac, log10(sc_mean + 1))) +
  geom_point() +
  facet_wrap(~cell_type)

ggsave(
  "sc-non-zero-frac-vs-sc-mean-plot.pdf",
  plot = p,
  width = 12,
  height = 8
)

p <- gene_prop_data |>
    mutate(
    sc_nz_pct = sc_non_zero_frac * 100,
    sc_nz_bin = case_when(
      sc_nz_pct == 0 ~ "0%",
      sc_nz_pct > 0 & sc_nz_pct < 0.1 ~ "0-0.1%",
      sc_nz_pct >= 0.1 & sc_nz_pct < 0.5 ~ "0.1-0.5%",
      sc_nz_pct >= 0.5 & sc_nz_pct < 1 ~ "0.5-1%",
      sc_nz_pct >= 1 & sc_nz_pct < 5 ~ "1-5%",
      sc_nz_pct >= 5 & sc_nz_pct < 10 ~ "5-10%",
      sc_nz_pct >= 10 & sc_nz_pct < 50 ~ "10-50%",
      sc_nz_pct >= 50 & sc_nz_pct <= 100 ~ "50-100%"
    ),
    sc_nz_bin = factor(
      sc_nz_bin,
      levels = c("0%", "0-0.1%", "0.1-0.5%", "0.5-1%", "1-5%", "5-10%", "10-50%", "50-100%")
    )
  ) |>
  count(cell_type, sc_nz_bin, name = "n_genes") |>
  complete(cell_type, sc_nz_bin, fill = list(n_genes = 0)) |>
  ggplot(aes(sc_nz_bin, n_genes, fill = cell_type)) +
  geom_col(position = "dodge2") +
  labs(
    x = "sc_non_zero_frac bin", 
    y = "Number of genes",
    fill = "Cell Type"
  )

ggsave(
  "sc-non-zero-frac-bin-gene-counts-plot.pdf",
  plot = p,
  width = 12,
  height = 8
)

