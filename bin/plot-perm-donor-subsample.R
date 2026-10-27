#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(purrr)
  library(tidyr)
  library(ggh4x)
  library(data.table)
  library(patchwork)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)

sc_data_files <- read_tsv(args[1], show_col_types = FALSE)
pb_data_files <- read_tsv(args[2], show_col_types = FALSE)
sc_obs_files <- read_tsv(args[3], show_col_types = FALSE)
cells_per_indiv <- read_tsv(args[4], show_col_types = FALSE)

cell_type_keep <- "CD4_NC"
models <- c("lm", "p_glmm_sc")
method_levels <- unname(method_lookup[models])

n_total <- cells_per_indiv |>
  filter(cell_label == cell_type_lookup[[cell_type_keep]]) |>
  summarise(n = n_distinct(individual)) |>
  pull(n)

n_donors_at_frac <- function(indiv_frac) {
  k <- floor(indiv_frac * n_total)
  k <- if_else(indiv_frac > 0, pmax(1, k), k)
  pmin(n_total, k)
}

add_n_donors <- function(data) {
  data |>
    mutate(
      indiv_frac = as.numeric(indiv_frac),
      n_donors = n_donors_at_frac(indiv_frac),
      n_label = paste0("n = ", n_donors),
      method_label = factor(method_lookup[method], levels = method_levels)
    )
}

compute_qq_data <- function(files) {
  pvalue <- unlist(map(files, function(x) fread(x, select = "pvalue")$pvalue))
  pvalue <- pvalue[!is.na(pvalue)]
  n <- length(pvalue)
  if (n == 0) {
    return(tibble(
      log_x_bin_mid = numeric(),
      log_y_pvalue = numeric(),
      log_lower_ci = numeric(),
      log_upper_ci = numeric()
    ))
  }

  m <- (1:n) / (n + 1)
  c <- abs(qnorm(0.05 / 2))
  v <- (1:n) * (n - (1:n) + 1) / (n + 1)^2 / (n + 2)
  s <- sqrt(v)
  lower_ci <- m - c * s
  upper_ci <- m + c * s

  log_x_pvalue <- -log10(m)
  y_pvalue <- sort(pvalue)
  x_bin <- cut(
    log_x_pvalue,
    breaks = seq(0, 6, by = 0.1),
    include.lowest = TRUE
  )

  tibble(x_bin, y_pvalue, lower_ci, upper_ci) |>
    summarise(
      log_y_pvalue = -log10(mean(y_pvalue)),
      log_lower_ci = -log10(mean(lower_ci)),
      log_upper_ci = -log10(mean(upper_ci)),
      .by = x_bin
    ) |>
    mutate(log_x_bin_mid = seq(0.05, 5.95, by = 0.1)[as.numeric(x_bin)])
}

build_sc_plot_data <- function(model_name) {
  sc_data_files |>
    filter(
      model == model_name,
      cov_spec == "bulk_pca",
      as.character(k) == "none",
      cell_type == cell_type_keep,
      as.numeric(cell_frac) == 1,
      as.numeric(count_frac) == 1
    ) |>
    mutate(method = model_name) |>
    summarise(file_list = list(variant_file), .by = c(method, indiv_frac)) |>
    rowwise() |>
    mutate(qq_data = list(compute_qq_data(file_list))) |>
    ungroup() |>
    select(-file_list) |>
    unnest(cols = qq_data)
}

build_pb_plot_data <- function(model_name) {
  pb_data_files |>
    filter(
      model == model_name,
      int_cov == "none",
      pb_type == "logcounts",
      cell_type == cell_type_keep,
      as.numeric(cell_frac) == 1,
      as.numeric(count_frac) == 1
    ) |>
    mutate(method = model_name) |>
    summarise(file_list = list(variant_file), .by = c(method, indiv_frac)) |>
    rowwise() |>
    mutate(qq_data = list(compute_qq_data(file_list))) |>
    ungroup() |>
    select(-file_list) |>
    unnest(cols = qq_data)
}

# TEMP: pseudobulk lm dropped; restore with build_pb_plot_data("lm").
perm_plot_data <- build_sc_plot_data("p_glmm_sc") |>
  add_n_donors()

n_label_levels <- perm_plot_data |>
  distinct(n_donors, n_label) |>
  arrange(n_donors) |>
  pull(n_label)

p_perm <- perm_plot_data |>
  mutate(n_label = factor(n_label, levels = n_label_levels)) |>
  ggplot(aes(
    log_x_bin_mid, log_y_pvalue,
    ymin = log_lower_ci, ymax = log_upper_ci
  )) +
  geom_point(alpha = 0.8, size = 0.8) +
  geom_abline(linetype = "dashed") +
  geom_ribbon(linetype = 2, alpha = 0.1) +
  facet_grid2(
    vars(method_label), vars(n_label),
    scales = "free_y",
    independent = "y"
  ) +
  labs(
    x = "Expected -log10(p-value)",
    y = "Observed -log10(p-value)"
  ) +
  theme_jp()

read_egenes <- function(manifest, method) {
  manifest |>
    mutate(method = method) |>
    rowwise() |>
    mutate(gene_tbl = list(read_tsv(region_file, show_col_types = FALSE))) |>
    unnest(gene_tbl) |>
    ungroup() |>
    mutate(
      region_bh = p.adjust(pvalue, method = "BH"),
      .by = c(method, indiv_frac)
    ) |>
    filter(region_bh < 0.05) |>
    select(method, indiv_frac, feature_id) |>
    distinct()
}

sc_manifest <- sc_obs_files |>
  filter(
    model == "p_glmm_sc",
    cov_spec == "bulk_pca",
    as.character(k) == "none",
    cell_type == cell_type_keep,
    as.numeric(cell_frac) == 1,
    as.numeric(count_frac) == 1,
    as.numeric(min_n_cells) < 0,
    int_cov == "none",
    ablate == "none"
  )

egenes <- read_egenes(sc_manifest, "p_glmm_sc")

egene_plot_data <- sc_manifest |>
  distinct(indiv_frac) |>
  mutate(method = "p_glmm_sc") |>
  add_n_donors() |>
  left_join(
    egenes |> summarise(n_egene = n(), .by = c(method, indiv_frac)),
    by = c("method", "indiv_frac")
  ) |>
  mutate(
    n_egene = replace_na(n_egene, 0L),
    pct_egene = 100 * n_egene / n_egene[indiv_frac == 1]
  )

p_egene <- egene_plot_data |>
  ggplot(aes(n_donors, pct_egene)) +
  geom_line(linewidth = 0.8, colour = method_col_lookup[[method_lookup[["p_glmm_sc"]]]]) +
  geom_point(size = 2, colour = method_col_lookup[[method_lookup[["p_glmm_sc"]]]]) +
  scale_x_continuous(breaks = sort(unique(egene_plot_data$n_donors))) +
  labs(
    x = "Number of donors",
    y = "% of full-data eGenes"
  ) +
  theme_jp()

n_frac <- n_distinct(perm_plot_data$n_donors)

ggsave(
  "perm-donor-subsample-plot.pdf",
  p_perm / p_egene + plot_layout(heights = c(2.2, 1)),
  width = max(8, 2.2 * n_frac),
  height = 10,
  limitsize = FALSE
)
