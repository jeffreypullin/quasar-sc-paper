#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

set.seed(300)

args <- commandArgs(trailingOnly = TRUE)

read_model_variants <- function(files, model_name) {
  files |>
    filter(
      .data$cell_type == "B_IN",
      .data$model == model_name,
      .data$cell_frac == 1,
      .data$indiv_frac == 1,
      .data$n_cells_target < 0,
      .data$count_frac == 1,
      .data$int_cov == "none"
    ) |>
    distinct(
      .data$dataset,
      .data$chr,
      .data$variant_file
    ) |>
    pmap_dfr(function(dataset, chr, variant_file) {
      variant_data <- read_tsv(
        variant_file,
        col_select = all_of(c("feature_id", "snp_id", "chrom", "beta", "se")),
        show_col_types = FALSE
      )
      variant_data$dataset <- dataset
      variant_data$manifest_chr <- chr
      variant_data
    })
}

read_sc_variants <- function(files) {
  if ("data_type" %in% names(files)) {
    files <- filter(files, data_type == "counts")
  } else if ("sc_type" %in% names(files)) {
    files <- filter(files, sc_type == "counts")
  }

  selected_files <- files |>
    filter(
      .data$cell_type == "B_IN",
      .data$model == "p_glmm_sc",
      .data$cov_spec == "bulk_pca",
      .data$cell_frac == 1,
      .data$indiv_frac == 1,
      .data$n_cells_target < 0,
      .data$count_frac == 1,
      .data$int_cov == "none",
      .data$k == "none"
    )

  if (n_distinct(selected_files$dataset) != 1) {
    stop("The SAIGE-QTL manifest lacks dataset IDs; exactly one dataset is required")
  }

  selected_files |>
    distinct(
      .data$dataset,
      .data$chr,
      .data$variant_file
    ) |>
    pmap_dfr(function(dataset, chr, variant_file) {
      variant_data <- read_tsv(
        variant_file,
        col_select = all_of(
          c("feature_id", "snp_id", "chrom", "alt", "ref", "beta", "se")
        ),
        show_col_types = FALSE
      )
      variant_data$dataset <- dataset
      variant_data$manifest_chr <- chr
      variant_data
    })
}

read_saigeqtl_variants <- function(files) {
  files |>
    filter(.data$cell_type == "B_IN") |>
    distinct(
      .data$chrom,
      .data$variant_file
    ) |>
    pmap_dfr(function(chrom, variant_file) {
      variant_data <- read_tsv(
        variant_file,
        col_select = all_of(
          c("gene", "MarkerID", "Allele1", "Allele2", "BETA", "SE")
        ),
        show_col_types = FALSE
      )
      variant_data$manifest_chr <- chrom
      variant_data
    })
}

pb_files <- read_tsv(args[[1]], show_col_types = FALSE)
sc_files <- read_tsv(args[[2]], show_col_types = FALSE)
saigeqtl_files <- read_tsv(args[[3]], show_col_types = FALSE)

nb_variants <- read_model_variants(pb_files, "nb_glm")
names(nb_variants)[names(nb_variants) == "beta"] <- "beta_nb"
names(nb_variants)[names(nb_variants) == "se"] <- "se_nb"

lm_variants <- read_model_variants(pb_files, "lm")
names(lm_variants)[names(lm_variants) == "beta"] <- "beta_lm"
names(lm_variants)[names(lm_variants) == "se"] <- "se_lm"

matched_variants <- inner_join(
  lm_variants,
  nb_variants,
  by = c("dataset", "manifest_chr", "chrom", "feature_id", "snp_id")
) |>
  filter(
    is.finite(.data$beta_lm),
    is.finite(.data$se_lm),
    .data$se_lm > 0,
    is.finite(.data$beta_nb),
    is.finite(.data$se_nb),
    .data$se_nb > 0
  ) |>
  mutate(
    z_lm = .data$beta_lm / .data$se_lm,
    z_nb = .data$beta_nb / .data$se_nb
  ) |>
  filter(is.finite(.data$z_lm), is.finite(.data$z_nb))

if (nrow(matched_variants) == 0) {
  stop("No matching B_IN SNP tests were found for the NB-GLM and LM")
}

n_sample <- min(10000, nrow(matched_variants))
plot_data <- matched_variants |>
  slice_sample(n = n_sample)

p <- ggplot(plot_data, aes(x = .data$z_lm, y = .data$z_nb)) +
  geom_point(alpha = 0.2, size = 0.5, color = cell_type_cols[["B_IN"]]) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.4) +
  coord_equal() +
  labs(
    title = "Pseudobulk model concordance",
    subtitle = paste0("B IN; n = ", scales::comma(n_sample), " SNP tests"),
    x = "Linear model z-score",
    y = "Negative binomial model z-score"
  ) +
  theme_jp()

ggsave(
  "zscore-scatter-B_IN.pdf",
  p,
  width = 8,
  height = 8
)

sc_variants <- read_sc_variants(sc_files) |>
  rename(
    beta_quasar = beta,
    se_quasar = se
  )

saigeqtl_variants <- read_saigeqtl_variants(saigeqtl_files) |>
  rename(
    beta_saigeqtl = BETA,
    se_saigeqtl = SE
  )

matched_saigeqtl <- inner_join(
  sc_variants,
  saigeqtl_variants,
  by = c(
    "manifest_chr",
    "feature_id" = "gene",
    "snp_id" = "MarkerID"
  )
)

allele_mismatches <- matched_saigeqtl |>
  filter(.data$alt != .data$Allele2 | .data$ref != .data$Allele1)

if (nrow(allele_mismatches) > 0) {
  stop(
    "Quasar and SAIGE-QTL alleles disagree for ",
    scales::comma(nrow(allele_mismatches)),
    " matched SNP tests"
  )
}

matched_saigeqtl <- matched_saigeqtl |>
  filter(
    is.finite(.data$beta_quasar),
    is.finite(.data$se_quasar),
    .data$se_quasar > 0,
    is.finite(.data$beta_saigeqtl),
    is.finite(.data$se_saigeqtl),
    .data$se_saigeqtl > 0
  ) |>
  mutate(
    z_quasar = .data$beta_quasar / .data$se_quasar,
    z_saigeqtl = .data$beta_saigeqtl / .data$se_saigeqtl
  ) |>
  filter(is.finite(.data$z_quasar), is.finite(.data$z_saigeqtl))

if (nrow(matched_saigeqtl) == 0) {
  stop("No matching B_IN SNP tests were found for P-GLMM and SAIGE-QTL")
}

n_saigeqtl_sample <- min(10000, nrow(matched_saigeqtl))
saigeqtl_plot_data <- matched_saigeqtl |>
  slice_sample(n = n_saigeqtl_sample)

p_saigeqtl <- ggplot(
  saigeqtl_plot_data,
  aes(x = .data$z_quasar, y = .data$z_saigeqtl)
) +
  geom_point(alpha = 0.2, size = 0.5, color = cell_type_cols[["B_IN"]]) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.4) +
  coord_equal() +
  labs(
    title = "Single-cell model concordance",
    subtitle = paste0(
      "B IN; n = ",
      scales::comma(n_saigeqtl_sample),
      " SNP tests"
    ),
    x = "P-GLMM z-score",
    y = "SAIGE-QTL z-score"
  ) +
  theme_jp()

ggsave(
  "zscore-scatter-saigeqtl-B_IN.pdf",
  p_saigeqtl,
  width = 8,
  height = 8
)
