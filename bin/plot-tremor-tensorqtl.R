#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(scales)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)
sample_files <- args[grepl("sample", basename(args))]
summary_files <- args[grepl("summary", basename(args))]

if (length(sample_files) == 0 || length(summary_files) == 0) {
  stop("Expected *sample*.tsv and *summary*.tsv arguments")
}

samples <- map_dfr(sample_files, read_tsv, show_col_types = FALSE)
summaries <- map_dfr(summary_files, read_tsv, show_col_types = FALSE)

write_tsv(summaries, "tremor-lm-vs-tensorqtl-summary.tsv")

label_cell_type <- function(x) {
  ifelse(x %in% names(cell_type_lookup), unname(cell_type_lookup[x]), x)
}

colour_cell_type <- function(x) {
  ifelse(x %in% names(cell_type_cols), unname(cell_type_cols[x]), "#667788")
}

cell_types <- sort(unique(c(samples$cell_type, summaries$cell_type)))
ct_labels <- setNames(label_cell_type(cell_types), cell_types)
ct_colours <- setNames(colour_cell_type(cell_types), cell_types)

annot <- summaries |>
  mutate(
    cell_type_label = factor(ct_labels[cell_type], levels = unname(ct_labels)),
    label = sprintf(
      "n = %s\nr = %.2f\nsign = %.0f%%",
      comma(n_finite),
      pearson_z,
      100 * sign_concordance
    )
  )

plot_data <- samples |>
  mutate(cell_type_label = factor(ct_labels[cell_type], levels = unname(ct_labels)))

p <- ggplot(plot_data, aes(z_quasar, z_tensorqtl, colour = cell_type)) +
  geom_point(alpha = 0.25, size = 0.45, show.legend = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.4) +
  geom_text(
    data = annot,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.2,
    size = 3.2,
    colour = "#222222"
  ) +
  scale_colour_manual(values = ct_colours) +
  facet_wrap(~cell_type_label, scales = "free") +
  labs(
    title = "Tremor QuASAR LM (TMM+INT) vs TensorQTL",
    x = "QuASAR pseudobulk LM z-score",
    y = "TensorQTL z-score"
  ) +
  theme_jp()

n_ct <- length(cell_types)
ggsave(
  "tremor-lm-vs-tensorqtl.pdf",
  p,
  width = max(8, 4 * min(n_ct, 3)),
  height = max(4, 4 * ceiling(n_ct / 3))
)
