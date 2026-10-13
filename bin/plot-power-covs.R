#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(forcats)
  library(tidyr)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

sc_data_files <- read_tsv(args[1], show_col_types = FALSE) |>
  filter(
    dataset == "onek1k",
    model == "p_glmm_sc",
    data_type == "counts",
    indiv_frac == 1,
    count_frac == 1,
    as.numeric(n_cells_target) < 0,
    ablate == "none"
  ) |>
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

cov_spec_lookup <- c(
  "bulk_pca" = "Pseudobulk PCs",
  "sc_pca" = "Single-cell PCs",
  "mean_sc_pca" = "Donor-mean single-cell PCs",
  "bulk_pca+pct_mito" = "Pseudobulk PCs + mitochondrial fraction",
  "bulk_pca+cell_cycle" = "Pseudobulk PCs + cell-cycle scores"
)

cov_spec_col_lookup <- c(
  "Pseudobulk PCs" = "#4477AA",
  "Single-cell PCs" = "#EE6677",
  "Donor-mean single-cell PCs" = "#AA3377",
  "Pseudobulk PCs + mitochondrial fraction" = "#228833",
  "Pseudobulk PCs + cell-cycle scores" = "#CCBB44"
)

plot_cov_spec <- function(df, n_sig, y_label, specs) {
  df |>
    filter(cell_type %in% c("Plasma", "B_IN", "CD4_NC"), cell_frac == 1, cov_spec %in% specs) |>
    mutate(
      cell_type = factor(cell_type, levels = c("Plasma", "B_IN", "CD4_NC")),
      cov_spec = factor(cov_spec_lookup[cov_spec], levels = unname(cov_spec_lookup))
    ) |>
    ggplot(aes(cell_type, {{ n_sig }}, fill = cov_spec)) +
    geom_col(position = position_dodge2(reverse = TRUE)) +
    scale_x_discrete(labels = cell_type_lookup) +
    scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = cov_spec_col_lookup) +
    labs(x = "Cell type", y = y_label, fill = "Covariates") +
    coord_flip() +
    theme_jp_vgrid()
}

plot_cov_spec_pair <- function(specs) {
  plot_cov_spec(sc_n_sig_var_data, n_sig_variant, "Number of significant variants", specs) +
    plot_cov_spec(sc_n_sig_gene_data, n_sig_gene, "Number of significant genes", specs) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top")
}

cov_p <- plot_cov_spec_pair(c("bulk_pca", "bulk_pca+pct_mito", "bulk_pca+cell_cycle"))
sc_pca_p <- plot_cov_spec_pair(c("bulk_pca", "sc_pca", "mean_sc_pca"))

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
  width = 16,
  height = 7
)

ggsave(
  "power-sc-pca-plot.pdf",
  plot = sc_pca_p,
  width = 16,
  height = 7
)

if (nrow(int_cov_plot_data) > 0) {
  ggsave(
    "power-int-cov-plot.pdf",
    plot = int_cov_p,
    width = 14,
    height = 8
  )
}
