#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(purrr)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)

neg_log10_p <- function(p) {
  p <- pmin(pmax(p, 1e-300), 1 - 1e-300)
  -log10(p)
}

effect_levels <- c("Main only", "Main in interaction", "Interaction")

egene_effect_specs <- function(int_cov, region_names) {
  if (identical(int_cov, "none")) {
    list(list(col = "pvalue", effect = "Main only"))
  } else {
    specs <- list()
    if ("main_acat_pvalue" %in% region_names) {
      specs <- c(specs, list(list(col = "main_acat_pvalue", effect = "Main in interaction")))
    }
    if ("int_acat_pvalue" %in% region_names) {
      specs <- c(specs, list(list(col = "int_acat_pvalue", effect = "Interaction")))
    }
    specs
  }
}

read_egene_pvalues <- function(manifest) {
  rows <- list()
  row_i <- 1L
  for (i in seq_len(nrow(manifest))) {
    info <- manifest[i, ]
    region <- read_tsv(info$region_file, show_col_types = FALSE)
    specs <- egene_effect_specs(info$int_cov, names(region))
    for (spec in specs) {
      if (!spec$col %in% names(region)) {
        next
      }
      rows[[row_i]] <- region |>
        select(feature_id, pvalue = all_of(spec$col)) |>
        filter(!is.na(feature_id), is.finite(pvalue), pvalue > 0, pvalue <= 1) |>
        mutate(
          cell_type = info$cell_type,
          model = info$model,
          int_cov = info$int_cov,
          effect = spec$effect
        )
      row_i <- row_i + 1L
    }
  }
  bind_rows(rows)
}

sc_data_files <- read_tsv(args[1], show_col_types = FALSE) |>
  filter(
    model %in% c("lmm_sc", "p_glmm_sc"),
    k == "none",
    cell_frac == 1,
    indiv_frac == 1
  )

if ("data_type" %in% names(sc_data_files)) {
  sc_data_files <- filter(
    sc_data_files,
    (model == "p_glmm_sc" & data_type == "counts") |
      (model == "lmm_sc" & data_type == "log_counts")
  )
} else if ("sc_type" %in% names(sc_data_files)) {
  sc_data_files <- filter(
    sc_data_files,
    (model == "p_glmm_sc" & sc_type == "counts") |
      (model == "lmm_sc" & sc_type %in% c("logcounts", "log_counts"))
  )
}
if ("count_frac" %in% names(sc_data_files)) {
  sc_data_files <- filter(sc_data_files, count_frac == 1)
}
if ("n_cells_target" %in% names(sc_data_files)) {
  sc_data_files <- filter(sc_data_files, n_cells_target < 0)
}

lmm_manifest <- filter(sc_data_files, model == "lmm_sc")
pglmm_manifest <- filter(sc_data_files, model == "p_glmm_sc")

if (nrow(lmm_manifest) == 0L || nrow(pglmm_manifest) == 0L) {
  stop("Need both lmm_sc and p_glmm_sc region files.")
}

lmm_genes <- read_egene_pvalues(lmm_manifest)
pglmm_genes <- read_egene_pvalues(pglmm_manifest)

if (nrow(lmm_genes) == 0L) {
  stop("No eGene p-values found for lmm_sc.")
}
if (nrow(pglmm_genes) == 0L) {
  stop("No eGene p-values found for p_glmm_sc.")
}

lmm_genes <- lmm_genes |>
  mutate(
    bh = p.adjust(pvalue, method = "BH"),
    .by = c(cell_type, int_cov, effect)
  )
pglmm_genes <- pglmm_genes |>
  mutate(
    bh = p.adjust(pvalue, method = "BH"),
    .by = c(cell_type, int_cov, effect)
  )

plot_data <- inner_join(
  lmm_genes |>
    select(cell_type, int_cov, effect, feature_id,
           lmm_pvalue = pvalue, lmm_bh = bh),
  pglmm_genes |>
    select(cell_type, int_cov, effect, feature_id,
           pglmm_pvalue = pvalue, pglmm_bh = bh),
  by = c("cell_type", "int_cov", "effect", "feature_id")
)

n_int_covs <- n_distinct(plot_data$int_cov[plot_data$int_cov != "none"])

plot_data <- plot_data |>
  mutate(
    cell_type_label = coalesce(unname(cell_type_lookup[cell_type]), cell_type),
    effect = factor(effect, levels = effect_levels),
    int_label = if_else(
      int_cov == "none" | n_int_covs <= 1,
      as.character(effect),
      paste(as.character(effect), int_cov, sep = " / ")
    ),
    panel = paste(cell_type_label, int_label, sep = " / "),
    neg_log10_lmm = neg_log10_p(lmm_pvalue),
    neg_log10_pglmm = neg_log10_p(pglmm_pvalue)
  ) |>
  arrange(cell_type_label, effect, int_cov)

if (nrow(plot_data) == 0L) {
  stop("No overlapping genes between lmm_sc and p_glmm_sc.")
}

model_summary <- plot_data |>
  summarise(
    n_genes = n(),
    spearman = cor(lmm_pvalue, pglmm_pvalue, method = "spearman"),
    .by = c(panel, cell_type_label, effect, int_cov)
  ) |>
  arrange(cell_type_label, effect, int_cov) |>
  mutate(
    label = paste0(
      panel, "\n",
      format(n_genes, big.mark = ","), " genes; ",
      "Spearman rho = ", sprintf("%.3f", spearman)
    )
  )

print(model_summary |> select(panel, n_genes, spearman))

scatter_data <- plot_data |>
  left_join(select(model_summary, panel, label), by = "panel") |>
  mutate(label = factor(label, levels = model_summary$label))

scatter_p <- scatter_data |>
  ggplot(aes(x = neg_log10_lmm, y = neg_log10_pglmm)) +
  geom_point(alpha = 0.35, size = 0.7, colour = "grey30") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.3) +
  facet_wrap(vars(label), scales = "free") +
  labs(
    x = expression(-log[10](LMM~eGene~pvalue)),
    y = expression(-log[10](P-GLMM~eGene~pvalue))
  ) +
  theme_jp()

n_panels <- n_distinct(scatter_data$label)
n_cols <- min(n_panels, 3)
ggsave(
  "sc-pglmm-vs-lmm-plot.pdf",
  scatter_p,
  width = max(6, 5.5 * n_cols),
  height = 6 * ceiling(n_panels / n_cols),
  limitsize = FALSE
)
