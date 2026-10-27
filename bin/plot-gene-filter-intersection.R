#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tidyr)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)
stats <- read_tsv(args[1], show_col_types = FALSE)
donor_summary <- read_tsv(args[2], show_col_types = FALSE)

cell_types <- c(
  "CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_ET", "CD8_NC", "CD8_S100B",
  "NK", "NK_R", "Plasma", "B_Mem", "B_IN", "Mono_C", "Mono_NC", "DC"
)

min_cells_levels <- c(0L, 5L, 10L, 20L)
min_cells_labels <- c(
  "0" = "All donors",
  "5" = ">= 5 cells",
  "10" = ">= 10 cells",
  "20" = ">= 20 cells"
)

# Sparsity cutoffs, including the pipeline 1% rule.
thresholds <- c(0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1)
threshold_labels <- c(
  "0" = "0%",
  "0.001" = "0.1%",
  "0.005" = "0.5%",
  "0.01" = "1%",
  "0.02" = "2%",
  "0.05" = "5%",
  "0.1" = "10%"
)

stats <- stats |>
  filter(cell_type %in% cell_types) |>
  mutate(
    cell_type = factor(
      cell_type_lookup[cell_type],
      levels = unname(cell_type_lookup[cell_types])
    ),
    min_cells_label = factor(
      min_cells_labels[as.character(min_cells)],
      levels = unname(min_cells_labels)
    )
  )

n_retained <- map_dfr(thresholds, function(t) {
  stats |>
    filter(sc_non_zero_frac > t) |>
    count(cell_type, min_cells_label, name = "n_genes") |>
    mutate(threshold = t)
}) |>
  complete(
    cell_type,
    min_cells_label,
    threshold = thresholds,
    fill = list(n_genes = 0L)
  ) |>
  mutate(
    threshold_label = factor(
      threshold_labels[as.character(threshold)],
      levels = unname(threshold_labels)
    ),
    is_pipeline = threshold == 0.01
  )

p_retained <- n_retained |>
  ggplot(aes(threshold_label, min_cells_label, fill = n_genes)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_tile(
    data = ~ filter(.x, is_pipeline),
    fill = NA,
    colour = "#222222",
    linewidth = 1.1
  ) +
  geom_text(
    aes(label = scales::comma(n_genes)),
    size = 3,
    colour = "#222222"
  ) +
  scale_fill_gradient(
    low = "#f7f7f7",
    high = "#882255",
    labels = scales::comma
  ) +
  facet_wrap(~cell_type) +
  labs(
    x = "Sparsity threshold (sc non-zero frac)",
    y = "Min cells per donor",
    fill = "Genes retained"
  ) +
  theme_jp() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    legend.title = element_text(
      family = "Helvetica", size = 22, color = "#222222"
    )
  )

ggsave(
  "gene-filter-n-genes-retained-plot.pdf",
  plot = p_retained,
  width = 18,
  height = 14
)

# Donor counts after each min-cells filter used in the gene heatmap.
donor_summary <- donor_summary |>
  filter(cell_type %in% cell_types) |>
  mutate(
    cell_type = factor(
      cell_type_lookup[cell_type],
      levels = unname(cell_type_lookup[cell_types])
    )
  )

p_donors_filter <- donor_summary |>
  ggplot(aes(min_cells, n_donors)) +
  geom_line(colour = "#44AA99", linewidth = 0.9) +
  geom_point(colour = "#44AA99", size = 2.5) +
  facet_wrap(~cell_type, scales = "free_y") +
  scale_x_continuous(breaks = min_cells_levels) +
  expand_limits(y = 0) +
  labs(
    x = "Min cells per donor",
    y = "Donors retained"
  ) +
  theme_jp()

ggsave(
  "gene-filter-donors-retained-plot.pdf",
  plot = p_donors_filter,
  width = 16,
  height = 10
)
