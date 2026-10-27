#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(purrr)
  library(tidyr)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: plot-grouped-vs-int.R <sc_quasar_file> <castie_file>")
}

neg_log10_p <- function(p) {
  p <- pmin(pmax(p, 1e-300), 1 - 1e-300)
  -log10(p)
}

analysis_levels <- c("P-GLMM", "CASTIE", "Linear grouped")

read_region_pvalues <- function(manifest, pval_col) {
  rows <- list()
  for (i in seq_len(nrow(manifest))) {
    info <- manifest[i, ]
    region <- read_tsv(info$region_file, show_col_types = FALSE)
    col <- gsub("{int_cov}", info$int_cov, pval_col, fixed = TRUE)
    if (!col %in% names(region)) {
      next
    }
    rows[[i]] <- region |>
      select(feature_id, pvalue = all_of(col)) |>
      filter(!is.na(feature_id), is.finite(pvalue), pvalue > 0, pvalue <= 1) |>
      mutate(
        cell_type = info$cell_type,
        int_cov = info$int_cov,
        data_type = if ("data_type" %in% names(info)) as.character(info$data_type) else NA_character_
      )
  }
  bind_rows(rows)
}

read_castie_pvalues <- function(manifest) {
  rows <- list()
  for (i in seq_len(nrow(manifest))) {
    info <- manifest[i, ]
    if (is.na(info$region_file) || identical(as.character(info$region_file), "NA")) {
      next
    }
    region <- read_tsv(info$region_file, show_col_types = FALSE)
    if ("int_acat_pvalue" %in% names(region) && "feature_id" %in% names(region)) {
      out <- region |>
        transmute(
          feature_id = as.character(feature_id),
          pvalue = as.numeric(int_acat_pvalue)
        )
    } else if (all(c("Gene", "pval_column", "ACAT_p") %in% names(region))) {
      out <- region |>
        filter(as.character(pval_column) == as.character(info$int_cov)) |>
        transmute(
          feature_id = as.character(Gene),
          pvalue = as.numeric(ACAT_p)
        )
    } else {
      next
    }
    rows[[i]] <- out |>
      filter(!is.na(feature_id), feature_id != "", is.finite(pvalue), pvalue > 0, pvalue <= 1) |>
      mutate(
        cell_type = info$cell_type,
        int_cov = info$int_cov,
        data_type = if ("data_type" %in% names(info)) as.character(info$data_type) else "counts"
      )
  }
  bind_rows(rows)
}

collapse_genes <- function(data) {
  data |>
    slice_min(pvalue, n = 1, with_ties = FALSE, by = c(cell_type, int_cov, data_type, feature_id))
}

add_bh <- function(data) {
  data |>
    mutate(
      bh = p.adjust(pvalue, method = "BH"),
      .by = c(cell_type, int_cov, data_type)
    )
}

panel_label <- function(cell_type, int_cov, data_type = NA_character_) {
  ct <- coalesce(unname(cell_type_lookup[cell_type]), cell_type)
  dt <- ifelse(is.na(data_type) | data_type == "", "counts", data_type)
  paste(ct, int_cov, dt, sep = " / ")
}

sc_data_files <- read_tsv(args[1], show_col_types = FALSE) |>
  filter(
    model == "p_glmm_sc",
    cell_frac == 1,
    indiv_frac == 1,
    int_cov != "none"
  )

if ("data_type" %in% names(sc_data_files)) {
  sc_data_files <- filter(sc_data_files, data_type %in% c("counts", "sct_counts"))
} else {
  sc_data_files$data_type <- "counts"
}
sc_data_files <- mutate(
  sc_data_files,
  data_type = ifelse(is.na(data_type) | data_type == "", "counts", data_type)
)
if ("sc_type" %in% names(sc_data_files)) {
  sc_data_files <- filter(sc_data_files, sc_type == "counts")
}
if ("count_frac" %in% names(sc_data_files)) {
  sc_data_files <- filter(sc_data_files, count_frac == 1)
}
if ("min_n_cells" %in% names(sc_data_files)) {
  sc_data_files <- filter(sc_data_files, min_n_cells < 0)
}

int_manifest <- sc_data_files |>
  filter(k == "none")
grouped_manifest <- sc_data_files |>
  filter(k != "none", k != "seacells")
castie_manifest <- read_tsv(args[2], show_col_types = FALSE)
if (!"data_type" %in% names(castie_manifest)) {
  castie_manifest$data_type <- "counts"
}
castie_manifest <- mutate(
  castie_manifest,
  data_type = ifelse(is.na(data_type) | data_type == "", "counts", data_type)
)

if (nrow(int_manifest) == 0L || nrow(grouped_manifest) == 0L) {
  stop("Need both P-GLMM interaction (k == none) and grouped (k != none) runs.")
}

pglmm_genes <- read_region_pvalues(int_manifest, "int_{int_cov}_acat_pvalue") |>
  mutate(data_type = ifelse(is.na(data_type) | data_type == "", "counts", data_type)) |>
  collapse_genes() |>
  add_bh() |>
  mutate(analysis = "P-GLMM")
linear_genes <- read_region_pvalues(grouped_manifest, "group_linear_acat_pvalue") |>
  mutate(data_type = ifelse(is.na(data_type) | data_type == "", "counts", data_type)) |>
  collapse_genes() |>
  add_bh() |>
  mutate(analysis = "Linear grouped")
castie_genes <- read_castie_pvalues(castie_manifest) |>
  mutate(data_type = ifelse(is.na(data_type) | data_type == "", "counts", data_type)) |>
  collapse_genes() |>
  add_bh() |>
  mutate(analysis = "CASTIE")

if (nrow(pglmm_genes) == 0L) {
  stop("No int_<int_cov>_acat_pvalue values found in P-GLMM interaction region files.")
}
if (nrow(linear_genes) == 0L) {
  stop("No group_linear_acat_pvalue values found in grouped region files.")
}
if (nrow(castie_genes) == 0L) {
  stop("No CASTIE interaction eGene p-values found.")
}

plot_data <- inner_join(
  pglmm_genes |>
    select(cell_type, int_cov, data_type, feature_id, pglmm_pvalue = pvalue),
  linear_genes |>
    select(cell_type, int_cov, data_type, feature_id, grouped_pvalue = pvalue),
  by = c("cell_type", "int_cov", "data_type", "feature_id")
) |>
  mutate(
    panel = panel_label(cell_type, int_cov, data_type),
    neg_log10_pglmm = neg_log10_p(pglmm_pvalue),
    neg_log10_grouped = neg_log10_p(grouped_pvalue)
  )

if (nrow(plot_data) == 0L) {
  stop("No overlapping genes between P-GLMM and grouped analyses.")
}

model_summary <- plot_data |>
  summarise(
    n_genes = n(),
    spearman = cor(pglmm_pvalue, grouped_pvalue, method = "spearman"),
    .by = panel
  ) |>
  arrange(panel) |>
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
  ggplot(aes(x = neg_log10_pglmm, y = neg_log10_grouped)) +
  geom_point(
    colour = method_col_lookup[["Single-cell P-GLMM"]],
    alpha = 0.35, size = 0.7
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.3) +
  facet_wrap(vars(label)) +
  labs(
    x = expression(-log[10](P-GLMM~interaction~eGene~pvalue)),
    y = expression(-log[10](Linear~grouped~eGene~pvalue))
  ) +
  theme_jp()

n_panels <- n_distinct(scatter_data$label)
ggsave(
  "grouped-vs-int-scatter-plot.pdf",
  scatter_p,
  width = max(6, 5.5 * min(n_panels, 2)),
  height = 6 * ceiling(n_panels / min(n_panels, 2)),
  limitsize = FALSE
)

analysis_keys <- c(
  "P-GLMM" = "p_glmm",
  "CASTIE" = "castie",
  "Linear grouped" = "grouped"
)

gene_calls <- bind_rows(pglmm_genes, castie_genes, linear_genes) |>
  mutate(
    analysis = factor(analysis, levels = analysis_levels),
    key = unname(analysis_keys[as.character(analysis)]),
    called = bh < 0.05
  )

called_wide <- gene_calls |>
  select(cell_type, int_cov, data_type, feature_id, key, called) |>
  pivot_wider(
    names_from = key,
    values_from = called,
    names_prefix = "called_",
    values_fill = FALSE
  )

p_wide <- gene_calls |>
  select(cell_type, int_cov, data_type, feature_id, key, pvalue, bh) |>
  pivot_wider(
    names_from = key,
    values_from = c(pvalue, bh)
  )

upset_genes <- called_wide |>
  left_join(p_wide, by = c("cell_type", "int_cov", "data_type", "feature_id")) |>
  mutate(
    called_p_glmm = coalesce(called_p_glmm, FALSE),
    called_castie = coalesce(called_castie, FALSE),
    called_grouped = coalesce(called_grouped, FALSE),
    methods = pmap_chr(
      list(called_p_glmm, called_castie, called_grouped),
      function(p_glmm, castie, grouped) {
        paste(analysis_levels[c(p_glmm, castie, grouped)], collapse = "; ")
      }
    ),
    n_methods = as.integer(called_p_glmm) +
      as.integer(called_castie) +
      as.integer(called_grouped)
  ) |>
  filter(n_methods > 0L) |>
  arrange(cell_type, int_cov, data_type, desc(n_methods), methods, feature_id) |>
  select(
    cell_type, int_cov, data_type, feature_id, methods, n_methods,
    called_p_glmm, called_castie, called_grouped,
    pvalue_p_glmm, bh_p_glmm,
    pvalue_castie, bh_castie,
    pvalue_grouped, bh_grouped
  )

write_tsv(upset_genes, "grouped-vs-int-upset-genes.tsv")
