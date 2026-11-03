#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)

sc_data_files <- read_tsv(args[2], show_col_types = FALSE)
saigeqtl_data_files <- read_tsv(args[3], show_col_types = FALSE)

example_cell_types <- c("Plasma", "B_IN", "CD4_NC")
example_cell_levels <- rev(unname(cell_type_lookup[c("CD4_NC", "B_IN", "Plasma")]))

status_levels <- c("SAIGE-QTL result", "No SAIGE-QTL result")

sc_files <- sc_data_files |>
  filter(cell_type %in% example_cell_types) |>
  filter(cell_frac == 1, indiv_frac == 1, min_n_cells < 0, count_frac == 1) |>
  filter(cov_spec == "bulk_pca", int_cov == "none", as.character(k) == "none") |>
  filter(model == "p_glmm_sc") |>
  filter(if ("data_type" %in% names(sc_data_files)) data_type == "counts" else TRUE) |>
  mutate(ablate = if ("ablate" %in% names(sc_data_files)) coalesce(as.character(ablate), "none") else "none") |>
  filter(ablate %in% c("none", ""))

read_quasar_genes <- function(files) {
  lapply(files, \(f) read_tsv(f, show_col_types = FALSE) |> select(feature_id, pvalue)) |>
    bind_rows() |>
    distinct(feature_id, .keep_all = TRUE)
}

read_saige_genes <- function(files) {
  files <- files[!is.na(files) & files != "NA"]
  lapply(files, \(f) read_tsv(f, show_col_types = FALSE) |> select(gene)) |>
    bind_rows() |>
    distinct(gene) |>
    pull(gene)
}

gene_data <- lapply(example_cell_types, function(ct) {
  quasar_genes <- read_quasar_genes(sc_files |> filter(cell_type == ct) |> pull(region_file))
  saige_genes <- read_saige_genes(saigeqtl_data_files |> filter(cell_type == ct) |> pull(region_file))

  quasar_genes |>
    mutate(
      cell_type = ct,
      status = if_else(feature_id %in% saige_genes, "SAIGE-QTL result", "No SAIGE-QTL result")
    )
}) |>
  bind_rows() |>
  mutate(
    cell_type = factor(unname(cell_type_lookup[cell_type]), levels = example_cell_levels),
    status = factor(status, levels = status_levels)
  )

fail_data <- tibble(cell_type = factor(example_cell_levels, levels = example_cell_levels)) |>
  left_join(
    gene_data |>
      filter(status == "No SAIGE-QTL result") |>
      count(cell_type, name = "n_gene"),
    by = "cell_type"
  ) |>
  mutate(n_gene = coalesce(n_gene, 0L))

p_fail <- fail_data |>
  ggplot(aes(cell_type, n_gene)) +
  geom_col(width = 0.7, fill = "grey35") +
  labs(
    x = "Cell type",
    y = "Non-converging SAIGE-QTL genes"
  ) +
  coord_flip() +
  theme_jp_vgrid()

n_sig_in_missing <- function(files, missing) {
  if (length(missing) == 0) {
    return(0)
  }
  sum(vapply(files, function(f) {
    x <- fread(f, select = c("feature_id", "pvalue"))
    sum(x$feature_id %in% missing & x$pvalue < 5e-6, na.rm = TRUE)
  }, numeric(1)))
}

variant_data <- lapply(example_cell_types, function(ct) {
  missing <- gene_data |>
    filter(cell_type == cell_type_lookup[[ct]], status == "No SAIGE-QTL result") |>
    pull(feature_id)
  tibble(
    cell_type = unname(cell_type_lookup[[ct]]),
    n_sig_variant = n_sig_in_missing(
      sc_files |> filter(cell_type == ct) |> pull(variant_file),
      missing
    )
  )
}) |>
  bind_rows() |>
  mutate(cell_type = factor(cell_type, levels = example_cell_levels))

p_var <- variant_data |>
  ggplot(aes(cell_type, n_sig_variant)) +
  geom_col(width = 0.7, fill = "grey35") +
  scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  labs(
    x = NULL,
    y = "Significant variants in genes\nwith no SAIGE-QTL result"
  ) +
  coord_flip() +
  theme_jp_vgrid()

ggsave(
  "convergence-saigeqtl-vs-quasar-example-celltypes-plot.pdf",
  plot = p_fail + p_var,
  width = 12,
  height = 5
)
