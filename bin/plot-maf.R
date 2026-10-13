#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(purrr)
  library(readr)
  library(tidyr)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)
pb_data_files <- read_tsv(args[1], show_col_types = FALSE)
sc_data_files <- read_tsv(args[2], show_col_types = FALSE)

is_full_data <- function(data) {
  data |>
    filter(
      as.numeric(cell_frac) == 1,
      as.numeric(indiv_frac) == 1,
      as.numeric(n_cells_target) < 0,
      as.numeric(count_frac) == 1,
      int_cov == "none"
    )
}

lookup_method <- function(key) {
  hit <- unname(method_lookup[key])
  ifelse(is.na(hit), key, hit)
}

cov_label <- c(bulk_pca = "bulk PCA", tqtl = "TensorQTL")

pb_manifest <- pb_data_files |>
  is_full_data() |>
  filter(model != "lm_cw") |>
  mutate(method = lookup_method(paste0("pb-", model)))

sc_manifest <- sc_data_files |>
  is_full_data() |>
  filter(as.character(k) == "none") |>
  mutate(multi_cov = n_distinct(cov_spec) > 1, .by = dataset) |>
  mutate(
    base = lookup_method(paste0("sc-", model)),
    cov_name = unname(cov_label[cov_spec]),
    cov_name = if_else(is.na(cov_name), as.character(cov_spec), cov_name),
    method = if_else(multi_cov, paste0(base, ", ", cov_name), base)
  )

onek1k_cell_types <- c("Plasma", "B_IN", "CD4_NC")

manifest <- bind_rows(
  pb_manifest |> transmute(dataset, cell_type, method, region_file, variant_file),
  sc_manifest |> transmute(dataset, cell_type, method, region_file, variant_file)
) |>
  filter(dataset != "onek1k" | cell_type %in% onek1k_cell_types)

sig_genes <- manifest |>
  distinct(dataset, cell_type, method, region_file) |>
  mutate(gene_tbl = map(region_file, function(path) {
    read_tsv(path, show_col_types = FALSE) |>
      select(feature_id, pvalue)
  })) |>
  select(-region_file) |>
  unnest(gene_tbl) |>
  mutate(
    region_bh = p.adjust(pvalue, method = "BH"),
    .by = c(dataset, cell_type, method)
  ) |>
  filter(region_bh < 0.05) |>
  distinct(dataset, cell_type, method, feature_id)

read_leads <- function(paths) {
  leads <- rbindlist(lapply(paths, function(path) {
    dt <- fread(
      path,
      select = c("feature_id", "snp_id", "maf", "pvalue"),
      showProgress = FALSE
    )
    dt <- dt[!is.na(feature_id) & !is.na(pvalue) & pvalue > 0 & !is.na(maf)]
    if (nrow(dt) == 0L) {
      return(data.table(
        feature_id = character(),
        snp_id = character(),
        maf = double(),
        pvalue = double()
      ))
    }
    dt[order(pvalue), .(snp_id = snp_id[1], maf = maf[1], pvalue = pvalue[1]), by = feature_id]
  }))
  if (nrow(leads) == 0L) {
    return(tibble(feature_id = character(), snp_id = character(), maf = double()))
  }
  as_tibble(leads[order(pvalue), .(snp_id = snp_id[1], maf = maf[1]), by = feature_id])
}

plot_data <- manifest |>
  semi_join(sig_genes, by = c("dataset", "cell_type", "method")) |>
  summarise(
    variant_file = list(unique(variant_file)),
    .by = c(dataset, cell_type, method)
  ) |>
  mutate(lead = map(variant_file, read_leads)) |>
  select(-variant_file) |>
  unnest(lead) |>
  inner_join(sig_genes, by = c("dataset", "cell_type", "method", "feature_id"))

label_dataset <- function(x) {
  case_when(
    x == "onek1k" ~ "OneK1K",
    x == "tremor" ~ "Tremor",
    startsWith(x, "randolph") ~ "Randolph",
    TRUE ~ x
  )
}

method_levels <- c(
  "Pseudobulk LM",
  "Pseudobulk NB-GLM",
  "Pseudobulk P-GLMM",
  "Single-cell P-GLMM",
  "Single-cell P-GLMM, bulk PCA",
  "Single-cell P-GLMM, TensorQTL",
  "Single-cell LMM",
  "Single-cell LMM, bulk PCA",
  "Single-cell LMM, TensorQTL"
)

type_cols <- c(
  method_col_lookup,
  "Single-cell P-GLMM, bulk PCA" = "#882255",
  "Single-cell P-GLMM, TensorQTL" = "#CC6677",
  "Single-cell LMM, bulk PCA" = "#44AA99",
  "Single-cell LMM, TensorQTL" = "#88CCEE"
)

label_cell <- function(x) {
  known <- x %in% names(cell_type_lookup)
  out <- x
  out[known] <- unname(cell_type_lookup[x[known]])
  out
}

plot_data <- plot_data |>
  mutate(dataset_label = label_dataset(dataset))

write_tsv(
  plot_data |>
    arrange(dataset, cell_type, method, feature_id) |>
    select(dataset, cell_type, method, feature_id, snp_id, maf),
  "maf-leads.tsv"
)

dataset_order <- c("OneK1K", "Tremor", "Randolph")
dataset_order <- c(
  intersect(dataset_order, unique(plot_data$dataset_label)),
  sort(setdiff(unique(plot_data$dataset_label), dataset_order))
)

make_plot <- function(df, title) {
  cell_labels <- label_cell(df$cell_type)
  cell_levels <- c(
    intersect(unname(cell_type_lookup), unique(cell_labels)),
    sort(setdiff(unique(cell_labels), unname(cell_type_lookup)))
  )
  methods <- intersect(method_levels, unique(df$method))
  methods <- c(methods, sort(setdiff(unique(df$method), methods)))
  cols <- type_cols[methods]
  cols[is.na(cols)] <- "#999999"

  df |>
    mutate(
      cell_type = factor(label_cell(cell_type), levels = cell_levels),
      method = factor(method, levels = methods)
    ) |>
    ggplot(aes(cell_type, maf, fill = method)) +
    geom_boxplot(
      position = position_dodge2(preserve = "single"),
      width = 0.7,
      linewidth = 0.3,
      outlier.size = 0.5
    ) +
    scale_fill_manual(values = cols, drop = FALSE) +
    labs(
      x = "Cell type",
      y = "Lead variant MAF",
      title = title
    ) +
    theme_jp() +
    theme(legend.position = "right")
}

plots <- map(dataset_order, function(ds) {
  make_plot(
    plot_data |> filter(dataset_label == ds),
    title = if (length(dataset_order) > 1) ds else NULL
  )
})

n_cells <- map_int(dataset_order, function(ds) {
  plot_data |>
    filter(dataset_label == ds) |>
    distinct(cell_type) |>
    nrow()
})

ggsave(
  "maf-plot.pdf",
  plot = wrap_plots(plotlist = plots, ncol = 1),
  width = max(8, 2.2 * max(n_cells) + 3),
  height = 4.5 * length(dataset_order)
)
