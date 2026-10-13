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

indiv_frac_breaks <- c(0.05, 0.1, 0.25, 0.5, 1)
indiv_frac_labels <- percent(indiv_frac_breaks, accuracy = NULL, drop0trailing = TRUE)
method_levels <- unname(method_lookup[c("lm", "lmm_sc", "p_glmm_sc")])

pb_manifest <- pb_data_files |>
  filter(
    model == "lm",
    pb_type == "logcounts",
    int_cov == "none",
    cell_frac == 1,
    count_frac == 1,
    n_cells_target < 0,
    cell_type == "CD4_NC"
  )

sc_manifest <- sc_data_files |>
  filter(
    (model == "lmm_sc" & data_type == "log_counts") |
      (model == "p_glmm_sc" & data_type == "counts"),
    cov_spec == "bulk_pca",
    as.character(k) == "none",
    int_cov == "none",
    cell_frac == 1,
    count_frac == 1,
    n_cells_target < 0,
    cell_type == "CD4_NC"
  )

time_data <- bind_rows(
  pb_manifest |> mutate(method = model),
  sc_manifest |> mutate(method = model)
) |>
  rowwise() |>
  mutate(time = read_time(time_file)) |>
  ungroup() |>
  summarise(time = sum(time), .by = c(indiv_frac, method)) |>
  mutate(
    method_label = factor(method_lookup[method], levels = method_levels),
    hours = time / 3600
  )

p <- time_data |>
  ggplot(aes(indiv_frac, hours, colour = method_label, group = method_label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_x_continuous(
    breaks = indiv_frac_breaks,
    labels = indiv_frac_labels,
    expand = expansion(mult = c(0.08, 0.03))
  ) +
  scale_colour_manual(values = method_col_lookup[method_levels]) +
  labs(
    x = "Donor fraction",
    y = "Time (h)",
    colour = "Method"
  ) +
  theme_jp() +
  theme(legend.position = "right")

ggsave("time-donor-subsample-plot.pdf", p, width = 8, height = 5)
