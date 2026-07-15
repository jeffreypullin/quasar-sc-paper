#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(data.table)
  library(ggplot2)
  library(forcats)
  library(patchwork)
  library(qvalue)
  library(purrr)
  library(stringr)
  library(scales)
  library(readxl)
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
saigeqtl_data_files <- read_tsv(args[3], show_col_types = FALSE)

saigeqtl_time_data <- saigeqtl_data_files |>
  filter(variant_file != "NA") |>
  summarise(
    step1_time = sum(step1_time),
    step2_time = sum(step2_time),
    .by = cell_type) |>
  mutate(time = step1_time + step2_time) |>
  mutate(cell_frac = 1, indiv_frac = 1)

sc_time_data <- sc_data_files |>
  filter(indiv_frac == 1) |>
  filter(cov_spec == "bulk_pca") |>
  rowwise() |>
  mutate(time = read_time(time_file)) |>
  ungroup() |>
  summarise(time = sum(time), .by = c(cell_type, cell_frac, indiv_frac))

pb_time_data <- pb_data_files |>
  filter(indiv_frac == 1) |>
  rowwise() |>
  mutate(time = read_time(time_file)) |>
  ungroup() |>
  summarise(time = sum(time), .by = c(cell_type, model, cell_frac, indiv_frac))

plot_data <- bind_rows(
  saigeqtl_time_data |>
    mutate(type = "saigeqtl"),
  sc_time_data |>
    mutate(type = "sc"),
  pb_time_data |>
    mutate(type = paste0("pb-", model)) |>
    select(-model),
)

quasar_p <- plot_data |>
  filter(cell_frac == 1) |>
  filter(indiv_frac == 1) |>
  filter(type != "saigeqtl") |>
  mutate(cell_type = fct_reorder(factor(cell_type), time)) |>
  ggplot(aes(cell_type, time, fill = type)) +
  geom_col(position = "dodge2") +
  coord_flip() +
  scale_y_continuous(labels = label_timespan())

ggsave(
  "time-quasar-plot.pdf",
  quasar_p,
  width = 12,
  height = 10
)

comparison_p <- plot_data |>
  filter(cell_frac == 1) |>
  filter(indiv_frac == 1) |>
  mutate(cell_type = fct_reorder(factor(cell_type), time)) |>
  filter(cell_type %in% c("Plasma", "B_IN", "CD4_NC")) |>
  filter(type %in% c("sc", "saigeqtl")) |>
  mutate(
    type = case_when(
      type == "sc" ~ "quasar",
      type == "saigeqtl" ~ "SAIGE-QTL",
    )
  ) |>
  ggplot(aes(cell_type, time, fill = type)) +
  geom_col(position = "dodge2") +
  coord_flip() +
  scale_y_continuous(
    labels = label_timespan(),
    breaks = breaks_width(24 * 60 * 60)
  ) +
  labs(
    x = "Cell type",
    y = "Time",
    fill = "Method"
  ) +
  scale_fill_manual(
    values = c("quasar" = "#228833", "SAIGE-QTL" = "#66CCEE")
  ) +
  theme_jp_vgrid()

ggsave(
  "time-method-comparison-plot.pdf",
  comparison_p,
  width = 12,
  height = 10
)

ggsave(
  "time-frac-plot.pdf",
  comparison_p,
  width = 12,
  height = 10
)
