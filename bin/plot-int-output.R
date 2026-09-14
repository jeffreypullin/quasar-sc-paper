#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(data.table)
  library(purrr)
  library(tidyr)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)

keep_full <- function(data) {
  data <- filter(data, int_cov != "none", cell_frac == 1, indiv_frac == 1)
  if ("count_frac" %in% names(data)) {
    data <- filter(data, count_frac == 1)
  }
  if ("n_cells_target" %in% names(data)) {
    data <- filter(data, n_cells_target < 0)
  }
  data
}

int_acat_col <- function(nms, int_cov) {
  preferred <- paste0("int_", to_snake(int_cov), "_acat_pvalue")
  if (preferred %in% nms) {
    return(preferred)
  }
  if ("int_acat_pvalue" %in% nms) {
    return("int_acat_pvalue")
  }
  NA_character_
}

pb_data_files <- keep_full(read_tsv(args[[1]], show_col_types = FALSE))

region_data <- if (nrow(pb_data_files) == 0) {
  tibble(
    cell_type = character(),
    model = character(),
    int_cov = character(),
    feature_id = character(),
    int_acat_pvalue = double()
  )
} else {
  pb_data_files |>
    mutate(
      gene_tbl = map2(region_file, int_cov, function(path, cov) {
        dt <- fread(path, showProgress = FALSE)
        col <- int_acat_col(names(dt), cov)
        if (is.na(col) || nrow(dt) == 0L) {
          return(tibble(feature_id = character(), int_acat_pvalue = double()))
        }
        tibble(
          feature_id = as.character(dt$feature_id),
          int_acat_pvalue = as.numeric(dt[[col]])
        )
      })
    ) |>
    select(cell_type, model, int_cov, gene_tbl) |>
    unnest(gene_tbl)
}

int_egenes <- region_data |>
  filter(!is.na(feature_id), feature_id != "", !is.na(int_acat_pvalue)) |>
  group_by(cell_type, model, int_cov, feature_id) |>
  slice_min(int_acat_pvalue, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    int_acat_bh = p.adjust(int_acat_pvalue, method = "BH"),
    .by = c(cell_type, model, int_cov)
  ) |>
  filter(int_acat_bh < 0.05) |>
  arrange(cell_type, model, int_cov, int_acat_bh) |>
  select(cell_type, model, int_cov, feature_id, int_acat_pvalue, int_acat_bh)

write_tsv(int_egenes, "int-res-egenes.tsv")
