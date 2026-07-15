#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(data.table)
  library(ggplot2)
  library(forcats)
  library(patchwork)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

pb_data_files <- read_tsv(args[1], show_col_types = FALSE) |>
  filter(int_cov == "none")
sc_data_files <- read_tsv(args[2], show_col_types = FALSE) |>
  filter(int_cov == "none") |>
  filter(k == "none")
saigeqtl_data_files <- read_tsv(args[3], show_col_types = FALSE)

saigeqtl_data_files

2 + "fdsml"

sc_n_na_data <- sc_data_files |>
  filter(cell_frac == 1) |>
  filter(indiv_frac == 1) |>
  rowwise() |>
  mutate(n_na = read_tsv(conv_file, show_col_types = FALSE) |>
    filter(is.na(pvalue)) |>
    nrow()) |>
  ungroup() |>
  summarise(n_na = sum(n_na), .by = cell_type)

pb_n_na_data <- pb_data_files |>
  filter(cell_frac == 1) |>
  filter(indiv_frac == 1) |>
  rowwise() |>
  mutate(n_na = read_tsv(conv_file, show_col_types = FALSE) |>
    filter(is.na(pvalue)) |>
    nrow()) |>
  ungroup() |>
  summarise(n_na = sum(n_na), .by = c(cell_type, model))

saigeqtl_n_na_data <- saigeqtl_data_files |>
  rowwise() |>
  mutate(n_na = count_na_variant_files(error_genes)) |>
  ungroup() |>
  summarise(n_na = sum(n_na), .by = cell_type)

sc_n_non_converged_data <- sc_data_files |>
  filter(cell_frac == 1) |>
  filter(indiv_frac == 1) |>
  rowwise() |>
  mutate(n_non_converged = read_tsv(conv_file, show_col_types = FALSE) |>
    filter(glmm_converged == 0) |>
    nrow()) |>
  ungroup() |>
  summarise(n_non_converged = sum(n_non_converged), .by = cell_type)

pb_n_non_converged_data <- pb_data_files |>
  filter(cell_frac == 1) |>
  filter(indiv_frac == 1) |>
  filter(model == "nb_glm") |>
  rowwise() |>
  mutate(n_non_converged = read_tsv(conv_file, show_col_types = FALSE) |>
    filter(glm_converged == 0) |>
    nrow()) |>
  ungroup() |>
  summarise(n_non_converged = sum(n_non_converged), .by = c(cell_type, model))

p_na <- bind_rows(
  pb_n_na_data |>
    mutate(type = paste0("pb-", model)),
  sc_n_na_data |>
    mutate(type = "sc"),
  saigeqtl_n_na_data |>
    mutate(type = "saigeqtl")
) |>
  mutate(cell_type = fct_reorder(factor(cell_type), n_na)) |>
  ggplot(aes(cell_type, n_na, fill = type)) +
  geom_col(position = "dodge2") +
  coord_flip()

p_zero <- bind_rows(
  pb_n_zero_data |>
    mutate(type = paste0("pb-", model)),
  sc_n_zero_data |>
    mutate(type = "sc")
) |>
  mutate(cell_type = fct_reorder(factor(cell_type), n_zero)) |>
  ggplot(aes(cell_type, n_zero, fill = type)) +
  geom_col(position = "dodge2") +
  coord_flip()

p_non_converged <- bind_rows(
  pb_n_non_converged_data |>
    mutate(type = paste0("pb-", model)),
  sc_n_non_converged_data |>
    mutate(type = "sc")
) |>
  mutate(cell_type = fct_reorder(factor(cell_type), n_non_converged)) |>
  ggplot(aes(cell_type, n_non_converged, fill = type)) +
  geom_col(position = "dodge2") +
  coord_flip()

p <- p_na + p_zero + p_non_converged +
  plot_layout(guides = "collect")

ggsave(
  "convergence-plot.pdf",
  plot = p,
  width = 8,
  height = 6
)
