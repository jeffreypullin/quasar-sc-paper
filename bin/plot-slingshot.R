#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
plot_tsv <- args[1]
output_prefix <- args[2]

dat <- read_tsv(plot_tsv, show_col_types = FALSE)
plot_df <- dat |> filter(kind == "cell")
curve_df <- dat |> filter(kind == "curve")

p_pt <- ggplot(plot_df, aes(PHATE_2, PHATE_1)) +
  geom_point(aes(colour = pseudotime), size = 0.3, alpha = 0.6) +
  scale_colour_viridis_c(na.value = "grey80") +
  geom_path(data = curve_df, aes(PHATE_2, PHATE_1, group = lineage), linewidth = 1) +
  theme_minimal() +
  labs(
    title = paste("Slingshot trajectory (pseudotime):", output_prefix),
    x = "PHATE 2",
    y = "PHATE 1",
    colour = "Pseudotime"
  )

ggsave(
  paste0(output_prefix, "-slingshot.pdf"),
  p_pt,
  width = 8,
  height = 6
)

p_ct <- ggplot(plot_df, aes(PHATE_2, PHATE_1)) +
  geom_path(
    data = curve_df,
    aes(PHATE_2, PHATE_1, group = lineage),
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  geom_point(aes(colour = cell_label), size = 0.35, alpha = 0.75) +
  scale_colour_brewer(type = "qual", palette = "Set1") +
  theme_minimal() +
  labs(
    title = paste("Slingshot trajectory (cell type):", output_prefix),
    x = "PHATE 2",
    y = "PHATE 1",
    colour = "Cell type"
  )

ggsave(
  paste0(output_prefix, "-slingshot-celltype.pdf"),
  p_ct,
  width = 8,
  height = 6
)
