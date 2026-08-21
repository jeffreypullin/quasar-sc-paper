#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(forcats)
  library(scales)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)

pb_data_files <- read_tsv(args[1], show_col_types = FALSE)
sc_data_files <- read_tsv(args[2], show_col_types = FALSE)

cell_types <- c("CD4_NC")

pb_manifest <- pb_data_files |>
  filter(
    model == "lm",
    int_cov == "none",
    cell_frac == 1,
    count_frac == 1,
    cell_type %in% cell_types
  )

sc_manifest <- sc_data_files |>
  filter(
    model == "p_glmm_sc",
    cov_spec == "bulk_pca",
    k == "none",
    cell_frac == 1,
    count_frac == 1,
    cell_type %in% cell_types
  )

read_egenes <- function(manifest, method) {
  manifest |>
    mutate(method = method) |>
    rowwise() |>
    mutate(gene_tbl = list(read_tsv(region_file, show_col_types = FALSE))) |>
    unnest(gene_tbl) |>
    ungroup() |>
    mutate(
      region_bh = p.adjust(pvalue, method = "BH"),
      .by = c(cell_type, method, indiv_frac)
    ) |>
    filter(region_bh < 0.05) |>
    select(cell_type, method, indiv_frac, feature_id) |>
    distinct()
}

egenes <- bind_rows(
  read_egenes(pb_manifest, "lm"),
  read_egenes(sc_manifest, "p_glmm_sc")
) |>
  mutate(
    method_label = factor(
      method_lookup[method],
      levels = c(method_lookup[["lm"]], method_lookup[["p_glmm_sc"]])
    ),
    cell_type_label = factor(
      cell_type_lookup[cell_type],
      levels = unname(cell_type_lookup[cell_types])
    )
  )

grid <- bind_rows(
  pb_manifest |> distinct(cell_type, indiv_frac) |> mutate(method = "lm"),
  sc_manifest |> distinct(cell_type, indiv_frac) |> mutate(method = "p_glmm_sc")
) |>
  mutate(
    method_label = factor(
      method_lookup[method],
      levels = c(method_lookup[["lm"]], method_lookup[["p_glmm_sc"]])
    ),
    cell_type_label = factor(
      cell_type_lookup[cell_type],
      levels = unname(cell_type_lookup[cell_types])
    )
  )

count_data <- grid |>
  left_join(
    egenes |> summarise(n_egene = n(), .by = c(cell_type, method, indiv_frac)),
    by = c("cell_type", "method", "indiv_frac")
  ) |>
  mutate(n_egene = replace_na(n_egene, 0L))

method_cols <- method_col_lookup[c(method_lookup[["lm"]], method_lookup[["p_glmm_sc"]])]

p_counts <- count_data |>
  ggplot(aes(indiv_frac, n_egene, colour = method_label, group = method_label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = sort(unique(count_data$indiv_frac))) +
  scale_colour_manual(values = method_cols) +
  facet_wrap(~cell_type_label, scales = "free_y") +
  labs(
    x = "Individual fraction",
    y = "Number of eGenes",
    colour = "Method"
  ) +
  theme_jp() +
  theme(legend.position = "right")

ggsave(
  "power-indiv-frac-counts-plot.pdf",
  p_counts,
  width = 8,
  height = 5
)

truth <- egenes |>
  filter(method == "lm", indiv_frac == 1) |>
  select(cell_type, feature_id) |>
  distinct()

truth_n <- truth |>
  summarise(n_truth = n(), .by = cell_type)

recall_data <- grid |>
  left_join(
    egenes |>
      inner_join(truth, by = c("cell_type", "feature_id")) |>
      summarise(n_recovered = n(), .by = c("cell_type", "method", "indiv_frac")),
    by = c("cell_type", "method", "indiv_frac")
  ) |>
  left_join(truth_n, by = "cell_type") |>
  mutate(
    n_recovered = replace_na(n_recovered, 0L),
    recall = if_else(n_truth > 0, n_recovered / n_truth, NA_real_)
  )

p_recall <- recall_data |>
  ggplot(aes(indiv_frac, recall, colour = method_label, group = method_label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = sort(unique(recall_data$indiv_frac))) +
  scale_y_continuous(limits = c(0, 1), labels = percent_format()) +
  scale_colour_manual(values = method_cols) +
  facet_wrap(~cell_type_label) +
  labs(
    x = "Individual fraction",
    y = "Proportion of full LM eGenes recovered",
    colour = "Method"
  ) +
  theme_jp() +
  theme(legend.position = "right")

ggsave(
  "power-indiv-frac-recall-plot.pdf",
  p_recall,
  width = 8,
  height = 5
)
