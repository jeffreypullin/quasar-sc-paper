#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(stringr)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

read_time <- function(path) {
  data <- read.delim(path)
  raw_str <- colnames(data)
  str <- str_extract(raw_str, "(?<=real\\.).*")
  as.numeric(str)
}

args <- commandArgs(trailingOnly = TRUE)
pb_data_files <- read_tsv(args[1], show_col_types = FALSE)
sc_data_files <- read_tsv(args[2], show_col_types = FALSE)
cluster_sizes <- read_tsv(args[3], show_col_types = FALSE) |>
  mutate(cell_type = str_replace_all(cell_label, " ", "_"))

method_levels <- unname(method_lookup[c("lm", "lm_cw", "lmm_sc", "p_glmm_sc")])

n_cells_by_type <- cluster_sizes |>
  filter(!cell_type %in% c("B_all", "T_all")) |>
  select(cell_type, n_cells)

pb_manifest <- pb_data_files |>
  filter(
    dataset == "onek1k",
    model %in% c("lm", "lm_cw"),
    pb_type == "logcounts",
    int_cov == "none",
    indiv_frac == 1,
    count_frac == 1,
    n_cells_target < 0
  )

sc_manifest <- sc_data_files |>
  filter(
    dataset == "onek1k",
    (model == "lmm_sc" & data_type == "log_counts") |
      (model == "p_glmm_sc" & data_type == "counts"),
    cov_spec == "bulk_pca",
    as.character(k) == "none",
    int_cov == "none",
    indiv_frac == 1,
    count_frac == 1,
    n_cells_target < 0
  )

time_data <- bind_rows(
  pb_manifest |> mutate(method = model),
  sc_manifest |> mutate(method = model)
) |>
  filter(!cell_type %in% c("B_all", "T_all")) |>
  inner_join(n_cells_by_type, by = "cell_type") |>
  rowwise() |>
  mutate(time = read_time(time_file)) |>
  ungroup() |>
  summarise(
    time = sum(time),
    n_cells = first(n_cells),
    .by = c(cell_type, cell_frac, method)
  ) |>
  mutate(
    method_label = factor(method_lookup[method], levels = method_levels),
    hours = time / 3600,
    n_cells = n_cells * cell_frac
  )

p <- time_data |>
  ggplot(aes(n_cells, hours, colour = method_label)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, alpha = 0.7) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_x_continuous(labels = label_comma()) +
  scale_colour_manual(values = method_col_lookup[method_levels], drop = FALSE) +
  labs(
    x = "Number of cells",
    y = "Time (h)",
    colour = "Method"
  ) +
  theme_jp() +
  theme(legend.position = "right")

ggsave("time-ncells-scatter-plot.pdf", p, width = 8, height = 5)
