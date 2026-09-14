#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

read_time <- function(path) {
  data <- read.delim(path)
  raw_str <- colnames(data)
  str <- str_extract(raw_str, "(?<=real\\.).*")
  as.numeric(str)
}

analysis_levels <- c("Main effect", "Interaction", "Grouped")
cell_type_levels <- c("B_all", "T_all")

args <- commandArgs(trailingOnly = TRUE)

sc_data_files <- read_tsv(args[1], show_col_types = FALSE) |>
  filter(
    cell_type %in% cell_type_levels,
    k != "seacells",
    cell_frac == 1,
    indiv_frac == 1,
    count_frac == 1,
    n_cells_target < 0
  )

time_data <- sc_data_files |>
  rowwise() |>
  mutate(time = read_time(time_file)) |>
  ungroup() |>
  mutate(
    analysis = case_when(
      as.character(k) != "none" ~ "Grouped",
      int_cov == "none" ~ "Main effect",
      TRUE ~ "Interaction"
    )
  ) |>
  summarise(time = sum(time), .by = c(cell_type, analysis)) |>
  complete(
    cell_type = cell_type_levels,
    analysis = analysis_levels
  ) |>
  mutate(
    analysis = factor(analysis, levels = analysis_levels),
    cell_type = factor(
      coalesce(unname(cell_type_lookup[cell_type]), cell_type),
      levels = unname(cell_type_lookup[cell_type_levels])
    ),
    hours = time / (60 * 60)
  )

p <- time_data |>
  ggplot(aes(cell_type, hours, fill = analysis)) +
  geom_col(position = "dodge2") +
  coord_flip() +
  scale_fill_discrete(drop = FALSE) +
  labs(
    x = "Cell type",
    y = "Time (h)",
    fill = NULL
  ) +
  theme_jp_vgrid() +
  theme(legend.position = "right")

ggsave("grouped-int-time-plot.pdf", p, width = 8, height = 4)
