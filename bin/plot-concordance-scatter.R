#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)

pb_data_files <- read_tsv(args[1], show_col_types = FALSE)
sc_data_files <- read_tsv(args[2], show_col_types = FALSE)
saigeqtl_data_files <- read_tsv(args[3], show_col_types = FALSE)

cell_types <- c("Plasma", "B_IN", "CD4_NC")
cell_type_levels <- unname(cell_type_lookup[cell_types])

pb_manifest <- pb_data_files |>
  filter(
    model == "lm",
    int_cov == "none",
    cell_frac == 1,
    indiv_frac == 1,
    min_n_cells < 0,
    count_frac == 1,
    cell_type %in% cell_types
  )

pb_id_manifest <- pb_data_files |>
  filter(
    model == "p_glmm_id",
    pb_type == "counts",
    int_cov == "none",
    cell_frac == 1,
    indiv_frac == 1,
    min_n_cells < 0,
    count_frac == 1,
    cell_type %in% cell_types
  )

sc_manifest <- sc_data_files |>
  filter(
    model == "p_glmm_sc",
    data_type == "counts",
    cov_spec == "bulk_pca",
    k == "none",
    int_cov == "none",
    cell_frac == 1,
    indiv_frac == 1,
    min_n_cells < 0,
    count_frac == 1,
    ablate %in% c("none", ""),
    cell_type %in% cell_types
  )

saige_manifest <- saigeqtl_data_files |>
  filter(variant_file != "NA", region_file != "NA", cell_type %in% cell_types)

read_region_bh <- function(manifest, pvalue_col, id_col) {
  if (nrow(manifest) == 0) {
    return(tibble(
      cell_type = character(),
      feature_id = character(),
      region_pvalue = double(),
      region_bh = double()
    ))
  }

  manifest |>
    rowwise() |>
    mutate(gene_tbl = list(
      read_tsv(region_file, show_col_types = FALSE) |>
        transmute(
          feature_id = .data[[id_col]],
          pvalue = .data[[pvalue_col]]
        )
    )) |>
    unnest(gene_tbl) |>
    ungroup() |>
    mutate(
      region_bh = p.adjust(pvalue, method = "BH"),
      .by = cell_type
    ) |>
    select(cell_type, feature_id, region_pvalue = pvalue, region_bh) |>
    distinct()
}

concordance_data <- function(x_data, y_data, value_col = "region_bh") {
  inner_join(
    x_data |> select(cell_type, feature_id, value_x = all_of(value_col)),
    y_data |> select(cell_type, feature_id, value_y = all_of(value_col)),
    by = c("cell_type", "feature_id")
  ) |>
    mutate(
      neg_log10_x = -log10(pmax(value_x, .Machine$double.xmin)),
      neg_log10_y = -log10(pmax(value_y, .Machine$double.xmin)),
      cell_type_label = factor(unname(cell_type_lookup[cell_type]), levels = cell_type_levels)
    )
}

plot_concordance_scatter <- function(data, x_lab, y_lab, equal_axes = FALSE) {
  if (equal_axes) {
    data <- data |>
      mutate(
        axis_limit = max(neg_log10_x, neg_log10_y, na.rm = TRUE),
        .by = cell_type_label
      )
  }

  p <- ggplot(data, aes(neg_log10_x, neg_log10_y)) +
    geom_point(alpha = 0.45, size = 0.6, colour = "grey20") +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.4) +
    facet_wrap(~cell_type_label, scales = "free") +
    labs(x = x_lab, y = y_lab) +
    theme_jp() +
    theme(axis.text = element_text(size = 13))

  if (equal_axes) {
    p <- p +
      geom_blank(aes(x = 0, y = 0)) +
      geom_blank(aes(x = axis_limit, y = axis_limit)) +
      theme(aspect.ratio = 1)
  }

  p
}

save_scatter <- function(data, path, x_lab, y_lab, equal_axes = FALSE) {
  if (nrow(data) == 0) {
    return(invisible(NULL))
  }
  p <- plot_concordance_scatter(data, x_lab, y_lab, equal_axes)
  ggsave(path, p, width = 12, height = 4.5)
  saveRDS(p, sub("\\.pdf$", ".rds", path))
}

pglmm_bh <- read_region_bh(sc_manifest, "pvalue", "feature_id")
lm_bh <- read_region_bh(pb_manifest, "pvalue", "feature_id")
pglmm_id_bh <- read_region_bh(pb_id_manifest, "pvalue", "feature_id")
saige_bh <- read_region_bh(saige_manifest, "ACAT_p", "gene")

save_scatter(
  concordance_data(lm_bh, pglmm_bh),
  "concordance-lm-p_glmm_sc-scatter.pdf",
  expression(-log[10] * " BH (Pseudobulk LM)"),
  expression(-log[10] * " BH (Single-cell P-GLMM)")
)

save_scatter(
  concordance_data(pglmm_id_bh, pglmm_bh),
  "concordance-p_glmm_id-p_glmm_sc-scatter.pdf",
  expression(-log[10] * " BH (Pseudobulk P-GLMM (iid))"),
  expression(-log[10] * " BH (Single-cell P-GLMM)")
)

save_scatter(
  concordance_data(saige_bh, pglmm_bh, value_col = "region_pvalue"),
  "concordance-saigeqtl-p_glmm_sc-scatter.pdf",
  expression(-log[10] * " P-value (SAIGE-QTL)"),
  expression(-log[10] * " P-value (Single-cell P-GLMM)"),
  equal_axes = TRUE
)
