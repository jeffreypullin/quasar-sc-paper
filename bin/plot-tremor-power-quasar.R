#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(patchwork)
  library(purrr)
  library(readr)
  library(tidyr)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)
pb_data_files <- read_tsv(args[1], show_col_types = FALSE)
sc_data_files <- read_tsv(args[2], show_col_types = FALSE)

is_full_data <- function(data) {
  data |>
    filter(
      dataset == "tremor",
      as.numeric(cell_frac) == 1,
      as.numeric(indiv_frac) == 1,
      as.numeric(min_n_cells) < 0,
      as.numeric(count_frac) == 1,
      int_cov == "none"
    )
}

count_sig_genes <- function(manifest, extra_by = character()) {
  manifest |>
    mutate(gene_tbl = map(region_file, read_tsv, show_col_types = FALSE)) |>
    unnest(gene_tbl) |>
    summarise(
      n_sig_gene = sum(p.adjust(pvalue, method = "BH") < 0.05, na.rm = TRUE),
      .by = all_of(c("cell_type", extra_by))
    )
}

count_sig_variants <- function(manifest, extra_by = character()) {
  manifest |>
    mutate(power_tbl = map(power_file, read_tsv, show_col_types = FALSE)) |>
    unnest(power_tbl) |>
    summarise(
      n_sig_variant = sum(n_sig_variant),
      .by = all_of(c("cell_type", extra_by))
    )
}

pb_manifest <- pb_data_files |>
  is_full_data() |>
  filter(model == "lm", pb_type == "tmm_int")

sc_manifest <- sc_data_files |>
  is_full_data() |>
  filter(
    cov_spec %in% c("bulk_pca", "tqtl"),
    as.character(k) == "none",
    model %in% c("p_glmm_sc", "lmm_sc")
  )

cell_types <- intersect(unique(pb_manifest$cell_type), unique(sc_manifest$cell_type))
pb_manifest <- pb_manifest |> filter(cell_type %in% cell_types)
sc_manifest <- sc_manifest |> filter(cell_type %in% cell_types)

cov_label <- c(bulk_pca = "bulk PCA", tqtl = "TensorQTL")
type_levels <- c(
  "Pseudobulk LM",
  "Single-cell P-GLMM, bulk PCA",
  "Single-cell P-GLMM, TensorQTL",
  "Single-cell LMM, bulk PCA",
  "Single-cell LMM, TensorQTL"
)
type_cols <- c(
  "Pseudobulk LM" = "#AA4499",
  "Single-cell P-GLMM, bulk PCA" = "#882255",
  "Single-cell P-GLMM, TensorQTL" = "#CC6677",
  "Single-cell LMM, bulk PCA" = "#44AA99",
  "Single-cell LMM, TensorQTL" = "#88CCEE"
)

sc_type <- function(model, cov_spec) {
  paste0(method_lookup[paste0("sc-", model)], ", ", cov_label[cov_spec])
}

variant_plot_data <- bind_rows(
  count_sig_variants(pb_manifest) |> mutate(type = "Pseudobulk LM"),
  count_sig_variants(sc_manifest, c("model", "cov_spec")) |>
    mutate(type = sc_type(model, cov_spec))
) |>
  mutate(
    cell_type = cell_type_lookup[cell_type],
    type = factor(type, levels = type_levels)
  ) |>
  mutate(cell_type = fct_reorder(factor(cell_type), n_sig_variant, .fun = max))

cell_type_levels <- levels(variant_plot_data$cell_type)

gene_plot_data <- bind_rows(
  count_sig_genes(pb_manifest) |> mutate(type = "Pseudobulk LM"),
  count_sig_genes(sc_manifest, c("model", "cov_spec")) |>
    mutate(type = sc_type(model, cov_spec))
) |>
  mutate(
    cell_type = factor(cell_type_lookup[cell_type], levels = cell_type_levels),
    type = factor(type, levels = type_levels)
  )

p_variant <- variant_plot_data |>
  ggplot(aes(cell_type, n_sig_variant, fill = type)) +
  geom_col(position = "dodge2") +
  scale_fill_manual(values = type_cols, drop = FALSE) +
  labs(
    y = "Number of significant variants",
    x = "Cell type",
    fill = "Method"
  ) +
  coord_flip() +
  theme_jp_vgrid()

p_gene <- gene_plot_data |>
  ggplot(aes(cell_type, n_sig_gene, fill = type)) +
  geom_col(position = "dodge2") +
  scale_fill_manual(values = type_cols, drop = FALSE) +
  labs(
    y = "Number of significant genes",
    x = NULL,
    fill = "Method"
  ) +
  coord_flip() +
  theme_jp_vgrid() +
  guides(fill = "none")

ggsave(
  "tremor-power-quasar-comparison-plot.pdf",
  plot = p_variant + p_gene +
    plot_layout(guides = "collect") &
    theme(legend.position = "right"),
  width = 18,
  height = 8
)
