#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(forcats)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)

pb_data_files <- read_tsv(args[1], show_col_types = FALSE)
sc_data_files <- read_tsv(args[2], show_col_types = FALSE)

cell_types <- c("Plasma", "B_IN", "CD4_NC")

pb_manifest <- pb_data_files |>
  filter(
    model == "nb_glm",
    int_cov == "none",
    cell_frac == 1,
    indiv_frac == 1,
    cell_type %in% cell_types
  )

sc_manifest <- sc_data_files |>
  filter(
    model == "p_glmm_sc",
    cov_spec == "bulk_pca",
    k == "none",
    cell_frac == 1,
    indiv_frac == 1,
    cell_type %in% cell_types
  )

read_regions <- function(manifest, method) {
  manifest |>
    mutate(method = method) |>
    rowwise() |>
    mutate(gene_tbl = list(read_tsv(region_file, show_col_types = FALSE))) |>
    unnest(gene_tbl) |>
    ungroup() |>
    select(cell_type, method, feature_id, region_pvalue = pvalue)
}

sig_status <- bind_rows(
  read_regions(pb_manifest, "nb_glm"),
  read_regions(sc_manifest, "p_glmm_sc")
) |>
  mutate(
    region_bh = p.adjust(region_pvalue, method = "BH"),
    .by = c(cell_type, method)
  ) |>
  mutate(is_egene = region_bh < 0.05) |>
  select(cell_type, method, feature_id, is_egene) |>
  pivot_wider(
    names_from = method,
    values_from = is_egene,
    names_prefix = "is_egene_"
  ) |>
  filter(!is.na(is_egene_p_glmm_sc), !is.na(is_egene_nb_glm)) |>
  mutate(
    is_unique_sc = as.integer(is_egene_p_glmm_sc & !is_egene_nb_glm)
  ) |>
  select(cell_type, feature_id, is_unique_sc)

gene_props <- bind_rows(
  pb_manifest |> distinct(cell_type, gene_prop_file),
  sc_manifest |> distinct(cell_type, gene_prop_file)
) |>
  distinct(cell_type, gene_prop_file) |>
  rowwise() |>
  mutate(props = list(read_tsv(gene_prop_file, show_col_types = FALSE))) |>
  unnest(props) |>
  ungroup() |>
  distinct(cell_type, feature_id, .keep_all = TRUE) |>
  select(
    cell_type, feature_id,
    sc_mean, sc_non_zero_frac,
    pb_mean, pb_non_zero_frac
  ) |>
  filter(
    !is.na(sc_mean),
    !is.na(sc_non_zero_frac),
    !is.na(pb_mean),
    !is.na(pb_non_zero_frac)
  ) |>
  mutate(
    log10_sc_mean = as.numeric(scale(log10(sc_mean + 1))),
    log10_pb_mean = as.numeric(scale(log10(pb_mean + 1))),
    sc_non_zero_frac_z = as.numeric(scale(sc_non_zero_frac)),
    pb_non_zero_frac_z = as.numeric(scale(pb_non_zero_frac)),
    .by = cell_type
  ) |>
  select(
    cell_type, feature_id,
    log10_sc_mean, log10_pb_mean,
    sc_non_zero_frac_z, pb_non_zero_frac_z
  )

model_data <- sig_status |>
  inner_join(gene_props, by = c("cell_type", "feature_id")) |>
  mutate(
    cell_type_label = factor(
      cell_type_lookup[cell_type],
      levels = unname(cell_type_lookup[cell_types])
    )
  )

predictor_labels <- c(
  "log10_sc_mean" = "sc mean (log10, scaled)",
  "sc_non_zero_frac_z" = "sc non-zero frac (scaled)",
  "log10_pb_mean" = "pb mean (log10, scaled)",
  "pb_non_zero_frac_z" = "pb non-zero frac (scaled)"
)

fit_logit <- function(df) {
  glm(
    is_unique_sc ~ log10_sc_mean + sc_non_zero_frac_z +
      log10_pb_mean + pb_non_zero_frac_z,
    data = df,
    family = binomial()
  )
}

tidy_logit <- function(fit) {
  est <- coef(fit)
  ci <- confint.default(fit)
  tibble(
    term = names(est),
    estimate = exp(est),
    conf.low = exp(ci[, 1]),
    conf.high = exp(ci[, 2]),
    p.value = summary(fit)$coefficients[, "Pr(>|z|)"]
  )
}

coef_data <- model_data |>
  group_by(cell_type_label) |>
  group_modify(~ {
    tidy_logit(fit_logit(.x)) |>
      filter(term != "(Intercept)") |>
      mutate(term_label = coalesce(predictor_labels[term], term))
  }) |>
  ungroup() |>
  mutate(
    cell_type_label = factor(
      cell_type_label,
      levels = unname(cell_type_lookup[cell_types])
    ),
    term_label = fct_rev(factor(
      term_label,
      levels = unname(predictor_labels)
    ))
  )

write_tsv(coef_data, "egene-sig-model-logit.tsv")

forest_p <- coef_data |>
  ggplot(aes(x = estimate, y = term_label, colour = cell_type_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_pointrange(
    aes(xmin = conf.low, xmax = conf.high),
    linewidth = 0.6,
    size = 0.5
  ) +
  scale_x_log10() +
  scale_colour_manual(
    values = setNames(cell_type_cols, cell_type_lookup[names(cell_type_cols)])
  ) +
  facet_wrap(~cell_type_label, ncol = 1) +
  labs(
    x = "Odds ratio (unique to single-cell P-GLMM eGene)",
    y = NULL,
    colour = "Cell type"
  ) +
  theme_jp_vgrid() +
  theme(legend.position = "none")

ggsave(
  "egene-sig-model-forest-plot.pdf",
  forest_p,
  width = 8,
  height = 9
)
