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

joint_fracs <- c(0.1, 0.25, 0.5, 0.75, 1.0)
method_levels <- unname(method_lookup[c("lm", "p_glmm_sc")])

pb_manifest <- pb_data_files |>
  filter(
    model == "lm",
    int_cov == "none",
    cell_type == "CD4_NC",
    count_frac == 1,
    cell_frac %in% joint_fracs,
    indiv_frac %in% joint_fracs
  )

sc_manifest <- sc_data_files |>
  filter(
    model == "p_glmm_sc",
    cov_spec == "bulk_pca",
    k == "none",
    int_cov == "none",
    cell_type == "CD4_NC",
    count_frac == 1,
    cell_frac %in% joint_fracs,
    indiv_frac %in% joint_fracs
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
      .by = c(dataset, method, cell_frac, indiv_frac)
    ) |>
    filter(region_bh < 0.05) |>
    select(dataset, method, cell_frac, indiv_frac, feature_id) |>
    distinct()
}

egenes <- bind_rows(
  read_egenes(pb_manifest, "lm"),
  read_egenes(sc_manifest, "p_glmm_sc")
)

grid <- bind_rows(
  pb_manifest |> distinct(dataset, cell_frac, indiv_frac) |> mutate(method = "lm"),
  sc_manifest |> distinct(dataset, cell_frac, indiv_frac) |> mutate(method = "p_glmm_sc")
) |>
  mutate(
    method_label = factor(method_lookup[method], levels = method_levels),
    cell_frac_label = factor(
      percent(cell_frac, accuracy = 1),
      levels = percent(joint_fracs, accuracy = 1)
    ),
    indiv_frac_label = factor(
      percent(indiv_frac, accuracy = 1),
      levels = percent(joint_fracs, accuracy = 1)
    )
  )

count_data <- grid |>
  left_join(
    egenes |>
      summarise(n_egene = n(), .by = c(dataset, method, cell_frac, indiv_frac)),
    by = c("dataset", "method", "cell_frac", "indiv_frac")
  ) |>
  mutate(n_egene = replace_na(n_egene, 0L))

p_heatmap <- count_data |>
  ggplot(aes(indiv_frac_label, cell_frac_label, fill = n_egene)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = n_egene), size = 3.5, colour = "#222222") +
  scale_fill_gradient(low = "#f7f7f7", high = "#882255") +
  facet_grid(dataset ~ method_label) +
  labs(
    x = "Individual fraction",
    y = "Cell fraction",
    fill = "eGenes"
  ) +
  theme_jp() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "right",
    legend.title = element_text(family = "Helvetica", size = 14, color = "#222222")
  )

ggsave(
  "power-both-frac-heatmap.pdf",
  p_heatmap,
  width = 10,
  height = 5
)
