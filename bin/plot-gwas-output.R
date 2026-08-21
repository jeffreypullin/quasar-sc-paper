#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(forcats)
  library(patchwork)
  library(qvalue)
  library(purrr)
  library(stringr)
  library(tidyr)
  library(ggh4x)
})

args <- commandArgs(trailingOnly = TRUE)

sc_gwas_data_files <- read_tsv(args[1], show_col_types = FALSE)
pb_gwas_data_files <- read_tsv(args[2], show_col_types = FALSE)

pb_col_types <- cols(
  feature_id = col_character(),
  snp_id = col_character(),
  chrom = col_character(),
  pos = col_double(),
  alt = col_character(),
  ref = col_character(),
  maf = col_double(),
  beta = col_double(),
  se = col_double(),
  pvalue = col_double(),
  glm_converged = col_double(),
  phi = col_double(),
  phi_converged = col_double()
)

sc_col_types <- cols(
  feature_id = col_character(),
  snp_id = col_character(),
  chrom = col_character(),
  pos = col_double(),
  alt = col_character(),
  ref = col_character(),
  maf = col_double(),
  beta = col_double(),
  se = col_double(),
  pvalue = col_double(),
  glmm_converged = col_double(),
  sigma2 = col_double()
)

pb_filtered <- pb_gwas_data_files |>
  filter(model == "nb_glm") |>
  rowwise() |>
  mutate(var_data = list(read_tsv(
    sig_variant_file,
    show_col_types = FALSE,
    col_types = pb_col_types
  ))) |>
  ungroup() |>
  unnest(var_data) |>
  select(feature_id, pheno_chr, geno_chr, snp_id, maf, beta, pvalue, model) |>
  filter(pheno_chr != geno_chr) |>
  filter(geno_chr != "chr6") |>
  # Filter out pseudogene driven DND1.
  filter(feature_id != "ENSG00000256453") |>
  # Filter out COX6A1P2 pseudogene.
  filter(feature_id != "ENSG00000226976") |>
  # Filter out EIF5AL1 with paralog.
  filter(feature_id != "ENSG00000253626")

sc_filtered <- sc_gwas_data_files |>
  rowwise() |>
  mutate(var_data = list(read_tsv(
    sig_variant_file,
    show_col_types = FALSE,
    col_types = sc_col_types
  ))) |>
  ungroup() |>
  unnest(var_data) |>
  select(feature_id, pheno_chr, geno_chr, snp_id, maf, beta, pvalue) |>
  filter(pheno_chr != geno_chr) |>
  filter(geno_chr != "chr6") |>
  # Filter out pseudogene driven DND1.
  filter(feature_id != "ENSG00000256453") |>
  # Filter out COX6A1P2 pseudogene.
  filter(feature_id != "ENSG00000226976") |>
  # Filter out EIF5AL1 with paralog.
  filter(feature_id != "ENSG00000253626")

gwas_data <- bind_rows(
  pb_filtered |>
    mutate(method = paste0("pseudobulk-", model)),
  sc_filtered |>
    mutate(method = "single-cell")
)

plot_data <- gwas_data |>
  expand_grid(threshold = c(5e-8, 5e-11)) |>
  summarise(
    count = sum(pvalue < threshold, na.rm = TRUE), .by = c(threshold, method)
  ) |>
  mutate(threshold = formatC(threshold, format = "e", digits = 0))

p <- plot_data |>
  ggplot(aes(x = threshold, y = count, fill = method)) +
  geom_col(position = "dodge") +
  labs(
    x = "Significance threshold",
    y = "Number of significant variants",
    fill = NULL
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("gwas-sig-counts-plot.pdf", p, width = 7, height = 5)

pb_filtered |>
  arrange(pvalue) |>
  print(n = 20)

sc_filtered |>
  arrange(pvalue) |>
  print(n = 100)
