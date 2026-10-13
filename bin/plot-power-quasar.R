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
      as.numeric(cell_frac) == 1,
      as.numeric(indiv_frac) == 1,
      as.numeric(n_cells_target) < 0,
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

label_randolph_cell_type <- function(data) {
  data |>
    mutate(
      cell_type = if_else(
        startsWith(dataset, "randolph_"),
        paste0(
          coalesce(unname(cell_type_lookup[cell_type]), cell_type),
          " (", sub("randolph_", "", dataset), ")"
        ),
        cell_type
      )
    )
}

pb_manifest <- pb_data_files |>
  is_full_data() |>
  label_randolph_cell_type() |>
  filter(model %in% c("lm", "lm_cw", "nb_glm", "p_glmm_id"))

sc_manifest <- sc_data_files |>
  is_full_data() |>
  label_randolph_cell_type() |>
  filter(
    cov_spec == "bulk_pca",
    as.character(k) == "none"
  ) |>
  mutate(model = if_else(ablate == "donor-totals", paste0(model, "_ablate"), model))

exclude_cell_types <- c("B All", "T All")

variant_plot_data <- bind_rows(
  count_sig_variants(pb_manifest, "model") |> mutate(type = paste0("pb-", model)),
  count_sig_variants(sc_manifest, "model") |> mutate(type = paste0("sc-", model))
) |>
  filter(type %in% names(method_lookup)) |>
  mutate(
    cell_type = coalesce(unname(cell_type_lookup[cell_type]), cell_type),
    type = method_lookup[type]
  ) |>
  filter(!cell_type %in% exclude_cell_types) |>
  mutate(
    cell_type = fct_reorder(factor(cell_type), n_sig_variant, .fun = max),
    type = fct_reorder(factor(type), n_sig_variant, .fun = sum)
  )

method_levels <- levels(variant_plot_data$type)

gene_plot_data <- bind_rows(
  count_sig_genes(pb_manifest, "model") |> mutate(type = paste0("pb-", model)),
  count_sig_genes(sc_manifest, "model") |> mutate(type = paste0("sc-", model))
) |>
  filter(type %in% names(method_lookup)) |>
  mutate(
    cell_type = coalesce(unname(cell_type_lookup[cell_type]), cell_type),
    type = factor(method_lookup[type], levels = method_levels)
  ) |>
  filter(!cell_type %in% exclude_cell_types) |>
  mutate(cell_type = fct_reorder(factor(cell_type), n_sig_gene, .fun = max))

make_power_plot <- function(variant_data, gene_data, method_keys) {
  method_levels <- unname(method_lookup[method_keys])
  variant_data <- variant_data |>
    filter(type %in% method_levels) |>
    mutate(type = factor(type, levels = method_levels))
  gene_data <- gene_data |>
    filter(type %in% method_levels) |>
    mutate(type = factor(type, levels = method_levels))

  p_variant <- variant_data |>
    ggplot(aes(cell_type, n_sig_variant, fill = type)) +
    geom_col(position = "dodge2") +
    scale_fill_manual(
      values = method_col_lookup,
      drop = FALSE,
      guide = guide_legend(reverse = TRUE)
    ) +
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    labs(
      y = "Number of significant variants",
      x = "Cell type",
      fill = "Method"
    ) +
    coord_flip() +
    theme_jp_vgrid()

  p_gene <- gene_data |>
    ggplot(aes(cell_type, n_sig_gene, fill = type)) +
    geom_col(position = "dodge2") +
    scale_fill_manual(values = method_col_lookup, drop = FALSE) +
    labs(
      y = "Number of significant genes",
      x = NULL,
      fill = "Method"
    ) +
    coord_flip() +
    theme_jp_vgrid() +
    guides(fill = "none")

  p_variant + p_gene +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")
}

all_methods <- c("lm", "nb_glm", "p_glmm_id", "lmm_sc", "p_glmm_sc")
comparison_methods <- c("lm", "nb_glm", "p_glmm_sc")

ggsave(
  "power-quasar-comparison-plot.pdf",
  plot = make_power_plot(variant_plot_data, gene_plot_data, comparison_methods),
  width = 22,
  height = 16
)

lm_counts <- gene_plot_data |>
  filter(type == method_lookup[["lm"]]) |>
  select(cell_type, n_egene_lm = n_sig_gene)

nb_glm_counts <- gene_plot_data |>
  filter(type == method_lookup[["nb_glm"]]) |>
  select(cell_type, n_egene_nb_glm = n_sig_gene)

egene_tbl <- gene_plot_data |>
  filter(!as.character(cell_type) %in% exclude_cell_types) |>
  filter(type %in% unname(method_lookup[all_methods])) |>
  mutate(method = factor(type, levels = unname(method_lookup[all_methods]))) |>
  left_join(lm_counts, by = "cell_type") |>
  left_join(nb_glm_counts, by = "cell_type") |>
  mutate(
    pct_diff_vs_lm = 100 * (n_sig_gene - n_egene_lm) / n_egene_lm,
    pct_diff_vs_nb_glm = 100 * (n_sig_gene - n_egene_nb_glm) / n_egene_nb_glm
  ) |>
  select(
    cell_type, method, n_egene = n_sig_gene,
    n_egene_lm, pct_diff_vs_lm,
    n_egene_nb_glm, pct_diff_vs_nb_glm
  )

avg_tbl <- egene_tbl |>
  summarise(
    cell_type = "Average",
    n_egene = mean(n_egene),
    n_egene_lm = mean(n_egene_lm),
    pct_diff_vs_lm = mean(pct_diff_vs_lm),
    n_egene_nb_glm = mean(n_egene_nb_glm),
    pct_diff_vs_nb_glm = mean(pct_diff_vs_nb_glm),
    .by = method
  )

total_tbl <- egene_tbl |>
  summarise(
    cell_type = "Total",
    n_egene = sum(n_egene),
    n_egene_lm = sum(n_egene_lm),
    pct_diff_vs_lm = 100 * (sum(n_egene) - sum(n_egene_lm)) / sum(n_egene_lm),
    n_egene_nb_glm = sum(n_egene_nb_glm),
    pct_diff_vs_nb_glm = 100 * (sum(n_egene) - sum(n_egene_nb_glm)) / sum(n_egene_nb_glm),
    .by = method
  )

bind_rows(egene_tbl, avg_tbl, total_tbl) |>
  mutate(
    cell_type = factor(
      as.character(cell_type),
      levels = c(
        setdiff(levels(gene_plot_data$cell_type), exclude_cell_types),
        "Average",
        "Total"
      )
    )
  ) |>
  arrange(cell_type, method) |>
  write_csv("power-quasar-egene-pct-diff.csv")

