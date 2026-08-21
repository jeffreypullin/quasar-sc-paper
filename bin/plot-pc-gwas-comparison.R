#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(forcats)
  library(stringr)
  library(scales)
  library(tidyr)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

read_time <- function(path) {
  data <- read.delim(path)
  raw_str <- colnames(data)
  str <- str_extract(raw_str, "(?<=real\\.).*")
  as.numeric(str)
}

args <- commandArgs(trailingOnly = TRUE)

pc_gwas_files <- read_tsv(args[1], show_col_types = FALSE)
pc_sc_gwas_files <- read_tsv(args[2], show_col_types = FALSE)

pc_data <- bind_rows(
  pc_gwas_files,
  pc_sc_gwas_files
)

time_data <- pc_data |>
  rowwise() |>
  mutate(time = read_time(time_file)) |>
  ungroup() |>
  mutate(
    cell_type = fct_reorder(factor(cell_type), time, .fun = max),
    model = factor(model, levels = c("lm", "lmm_sc"))
  )

p_time <- time_data |>
  ggplot(aes(cell_type, time, fill = model)) +
  geom_col(position = "dodge2") +
  coord_flip() +
  scale_y_continuous(labels = label_timespan()) +
  labs(
    x = "Cell type",
    y = "Time",
    fill = "Model"
  ) +
  theme_jp_vgrid()

ggsave(
  "time-pc-gwas-plot.pdf",
  p_time,
  width = 12,
  height = 10
)

pc_col_types <- cols(
  feature_id = col_character(),
  snp_id = col_character(),
  chrom = col_character(),
  pos = col_double(),
  alt = col_character(),
  ref = col_character(),
  maf = col_double(),
  beta = col_double(),
  se = col_double(),
  pvalue = col_double()
)

gws_pval_threshold <- 5e-8

power_data <- pc_data |>
  rowwise() |>
  mutate(var_data = list(read_tsv(
    sig_variant_file,
    show_col_types = FALSE,
    col_types = pc_col_types
  ))) |>
  ungroup() |>
  unnest(var_data) |>
  summarise(
    n_sig_variant = sum(!is.na(pvalue) & pvalue < gws_pval_threshold),
    .by = c(cell_type, model)
  ) |>
  mutate(
    cell_type = fct_reorder(factor(cell_type), n_sig_variant, .fun = max),
    model = factor(model, levels = c("lm", "lmm_sc"))
  )

p_power <- power_data |>
  ggplot(aes(cell_type, n_sig_variant, fill = model)) +
  geom_col(position = "dodge2") +
  coord_flip() +
  labs(
    x = "Cell type",
    y = "Number of genome-wide significant variants",
    fill = "Model"
  ) +
  theme_jp_vgrid()

ggsave(
  "power-pc-gwas-plot.pdf",
  p_power,
  width = 12,
  height = 10
)
