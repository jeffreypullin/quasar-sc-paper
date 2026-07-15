#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(coloc)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

sc_quasar_file <- args[1]

fcrl5_data <- read_tsv(sc_quasar_file, show_col_types = FALSE) |>
  filter(k == 5) |>
  filter(chr == "chr1") |>
  rowwise() |>
  mutate(var_data = list(read_tsv(variant_file, show_col_types = FALSE))) |>
  unnest(var_data) |>
  filter(feature_id == "ENSG00000143297")

gwas_data <- read_tsv(
  "/home/jp2045/quasar-sc-paper/data/gwas/test.tsv",
  show_col_types = FALSE
)

data1 <- fcrl5_data |>
  mutate(varbeta = se^2) |>
  select(
    snp = snp_id,
    beta = beta,
    varbeta = varbeta,
    position = pos
  ) |>
  as.list()
data1$type <- "quant"
data1$sdY <- 1.1

data2 <- gwas_data |>
  mutate(
    snp = paste0(chromosome, ":", base_pair_location, toupper(effect_allele), "-", toupper(other_allele))
  ) |>
  mutate(varbeta = standard_error^2) |>
  select(
    snp = snp,
    beta = beta,
    varbeta = varbeta,
    position = base_pair_location,
  ) |>
  as.list()
data2$type <- "quant"
data2$sdY <- 1.1

coloc.abf(data1, data2)

write_tsv(sc_data, "coloc-results.tsv")
