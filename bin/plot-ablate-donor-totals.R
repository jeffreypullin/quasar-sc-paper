#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)

ablate_lookup <- c(
  "donor-totals" = "Single-cell data",
  "none" = "Donor totals"
)

ablate_col_lookup <- c(
  "Single-cell data" = "#882255",
  "Donor totals" = "#44AA99"
)

read_time <- function(path) {
  read_lines(path) |>
    str_subset("^real") |>
    str_extract("[0-9.]+") |>
    as.numeric()
}

manifest <- read_tsv(args[1], show_col_types = FALSE) |>
  filter(
    dataset == "onek1k",
    cell_type %in% c("Plasma", "B_IN", "CD4_NC"),
    model == "p_glmm_sc",
    data_type == "counts",
    cov_spec == "bulk_pca",
    as.numeric(cell_frac) == 1,
    as.numeric(indiv_frac) == 1,
    as.numeric(min_n_cells) < 0,
    as.numeric(count_frac) == 1,
    as.character(k) == "none",
    int_cov == "none",
    ablate %in% names(ablate_lookup)
  ) |>
  mutate(cell_type = coalesce(unname(cell_type_lookup[cell_type]), cell_type))

plot_data <- manifest |>
  mutate(gene_tbl = map(region_file, \(f) read_tsv(f, show_col_types = FALSE) |> select(feature_id, pvalue))) |>
  select(cell_type, ablate, gene_tbl) |>
  unnest(gene_tbl) |>
  pivot_wider(names_from = ablate, values_from = pvalue) |>
  filter(!is.na(none), !is.na(`donor-totals`))

p <- plot_data |>
  ggplot(aes(-log10(`donor-totals`), -log10(none))) +
  geom_abline(slope = 1, intercept = 0, colour = "#cbcbcb", linewidth = 1) +
  geom_point(alpha = 0.3, size = 1.5) +
  facet_wrap(~cell_type, ncol = 3) +
  labs(
    x = expression(-log[10] ~ "eGene p, single-cell data"),
    y = expression(-log[10] ~ "eGene p, donor totals")
  ) +
  theme_jp()

ggsave("ablate-donor-totals-plot.pdf", plot = p, width = 12, height = 4.5)

time_data <- manifest |>
  mutate(time = map_dbl(time_file, read_time)) |>
  summarise(hours = sum(time) / 3600, .by = c(cell_type, ablate)) |>
  mutate(
    method = factor(ablate_lookup[ablate], levels = unname(ablate_lookup)),
    cell_type = fct_reorder(cell_type, hours, .fun = max)
  )

p_time <- time_data |>
  ggplot(aes(cell_type, hours, fill = method)) +
  geom_col(position = "dodge2") +
  coord_flip() +
  scale_fill_manual(values = ablate_col_lookup) +
  labs(x = "Cell type", y = "Time (h)", fill = "Method") +
  theme_jp_vgrid()

ggsave("ablate-donor-totals-time-plot.pdf", plot = p_time, width = 10, height = 10)

na_data <- manifest |>
  mutate(
    n_na = map_int(conv_file, \(f) read_tsv(f, show_col_types = FALSE) |>
      filter(is.na(pvalue)) |>
      distinct(feature_id, snp_id) |>
      nrow()),
    n_variant = map_int(variant_file, \(f) nrow(fread(f, select = "feature_id", showProgress = FALSE)))
  ) |>
  summarise(pct_na = 100 * sum(n_na) / sum(n_variant), .by = c(cell_type, ablate)) |>
  mutate(
    method = factor(ablate_lookup[ablate], levels = unname(ablate_lookup)),
    cell_type = fct_reorder(cell_type, pct_na, .fun = max)
  )

p_na <- na_data |>
  ggplot(aes(cell_type, pct_na, fill = method)) +
  geom_col(position = "dodge2") +
  coord_flip() +
  scale_fill_manual(values = ablate_col_lookup) +
  labs(x = "Cell type", y = "Variants with NA p-value (%)", fill = "Method") +
  theme_jp_vgrid()

ggsave("ablate-donor-totals-na-plot.pdf", plot = p_na, width = 10, height = 10)
