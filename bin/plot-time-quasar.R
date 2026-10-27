#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
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

method_levels <- unname(method_lookup[c("lm", "p_glmm_id", "p_glmm_sc", "p_glmm_sc_ablate")])

is_full_data <- function(data) {
  data |>
    filter(
      dataset == "onek1k",
      as.numeric(cell_frac) == 1,
      as.numeric(indiv_frac) == 1,
      as.numeric(min_n_cells) < 0,
      as.numeric(count_frac) == 1,
      int_cov == "none"
    )
}

sum_time <- function(manifest) {
  manifest |>
    rowwise() |>
    mutate(time = read_time(time_file)) |>
    ungroup() |>
    summarise(time = sum(time), .by = c(cell_type, model))
}

sc_time_data <- sc_data_files |>
  is_full_data() |>
  filter(
    model == "p_glmm_sc",
    cov_spec == "bulk_pca",
    as.character(k) == "none"
  ) |>
  mutate(model = if_else(ablate == "donor-totals", paste0(model, "_ablate"), model)) |>
  sum_time()

pb_time_data <- pb_data_files |>
  is_full_data() |>
  filter(model %in% c("lm", "p_glmm_id")) |>
  sum_time()

plot_data <- bind_rows(
  sc_time_data |> mutate(type = paste0("sc-", model)),
  pb_time_data |> mutate(type = paste0("pb-", model))
) |>
  filter(type %in% names(method_lookup)) |>
  mutate(
    cell_type = coalesce(unname(cell_type_lookup[cell_type]), cell_type),
    type = factor(method_lookup[type], levels = method_levels),
    hours = time / 3600
  )

p <- plot_data |>
  mutate(cell_type = fct_reorder(cell_type, hours, .fun = max)) |>
  ggplot(aes(cell_type, hours, fill = type)) +
  geom_col(position = "dodge2") +
  coord_flip() +
  scale_y_continuous(breaks = breaks_width(1)) +
  scale_fill_manual(values = method_col_lookup, drop = FALSE) +
  labs(
    x = "Cell type",
    y = "Time (h)",
    fill = "Method"
  ) +
  theme_jp_vgrid() +
  theme(legend.position = "right")

ggsave("time-quasar-plot.pdf", p, width = 14, height = 12)
