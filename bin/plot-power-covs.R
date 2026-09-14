#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(forcats)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

sc_data_files <- read_tsv(args[1], show_col_types = FALSE) |>
  filter(indiv_frac == 1) |>
  filter(int_cov == "none" | as.character(k) != "none")

sc_n_sig_gene_data <- sc_data_files |>
  rowwise() |>
  mutate(test = list(read_tsv(region_file, show_col_types = FALSE))) |>
  unnest(test) |>
  summarise(
    n_sig_gene = sum(p.adjust(pvalue, method = "BH") < 0.05, na.rm = TRUE),
    .by = c(cell_type, cov_spec, cell_frac)
  )

sc_n_sig_var_data <- sc_data_files |>
  rowwise() |>
  mutate(n_sig_variant = list(read_tsv(power_file, show_col_types = FALSE))) |>
  unnest(n_sig_variant) |>
  summarise(n_sig_variant = sum(n_sig_variant), .by = c(cell_type, cell_frac, cov_spec))

cov_p <- sc_n_sig_var_data |>
  summarise(n_sig_variant = sum(n_sig_variant), .by = c(cell_type, cov_spec, cell_frac)) |>
  filter(cell_type %in% c("Plasma", "B_IN", "CD4_NC")) |>
  filter(cell_frac == 1) |>
  mutate(cell_type = fct_reorder(factor(cell_type), n_sig_variant)) |>
  ggplot(aes(cell_type, n_sig_variant, fill = cov_spec)) +
  geom_col(position = "dodge2") +
  labs(
    x = "Cell type",
    y = "Number of significant variants",
    fill = "Covariate specification"
  ) +
  coord_flip() +
  theme_jp_vgrid()

cov_gene_p <- sc_n_sig_gene_data |>
  filter(cell_type %in% c("Plasma", "B_IN", "CD4_NC")) |>
  filter(cell_frac == 1) |>
  mutate(cell_type = fct_reorder(factor(cell_type), n_sig_gene)) |>
  ggplot(aes(cell_type, n_sig_gene, fill = cov_spec)) +
  geom_col(position = "dodge2") +
  labs(
    x = "Cell type",
    y = "Number of significant genes",
    fill = "Covariate specification"
  ) +
  coord_flip() +
  theme_jp_vgrid()

filter_int_cov_data <- function(df) {
  df |>
    filter(cell_frac == 1) |>
    filter(
      (cell_type == "B_all" & cov_spec %in% c("bulk_pca", "bulk_pca+pseudotime")) |
        (cell_type == "T_all" & (cov_spec == "bulk_pca" | startsWith(cov_spec, "bulk_pca+starcat_")))
    )
}

int_cov_plot_data <- bind_rows(
  filter_int_cov_data(sc_n_sig_var_data) |>
    mutate(level = "Number of eQTLS", n_sig = n_sig_variant),
  filter_int_cov_data(sc_n_sig_gene_data) |>
    mutate(level = "Number of eGenes", n_sig = n_sig_gene)
) |>
  mutate(
    cov_type = if_else(cov_spec == "bulk_pca", "Bulk PCA", "Bulk PCA + additional covariate"),
    cov_type = factor(cov_type, levels = c("Bulk PCA", "Bulk PCA + additional covariate")),
    cell_type = factor(cell_type, levels = c("B_all", "T_all"))
  )

int_cov_p <- int_cov_plot_data |>
  ggplot(aes(cov_type, n_sig, fill = cov_type)) +
  geom_col() +
  facet_grid(level ~ cell_type, scales = "free_y") +
  labs(
    x = NULL,
    y = "Number of significant associations",
    fill = "Covariate specification"
  ) +
  scale_fill_manual(values = c(
    "Bulk PCA" = "#4477AA",
    "Bulk PCA + additional covariate" = "#EE6677"
  )) +
  theme_jp()

ggsave(
  "power-cov-spec-plot.pdf",
  plot = cov_p,
  width = 12,
  height = 8
)

ggsave(
  "power-cov-spec-gene-plot.pdf",
  plot = cov_gene_p,
  width = 12,
  height = 8
)

ggsave(
  "power-int-cov-plot.pdf",
  plot = int_cov_p,
  width = 14,
  height = 8
)
