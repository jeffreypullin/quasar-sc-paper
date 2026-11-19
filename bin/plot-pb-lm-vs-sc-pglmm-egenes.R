#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)

pb_data_files <- read_tsv(args[1], show_col_types = FALSE)
sc_data_files <- read_tsv(args[2], show_col_types = FALSE)
annot <- read_tsv(
  args[3],
  show_col_types = FALSE,
  col_types = cols(
    `#chr` = col_character(),
    start = col_double(),
    end = col_double(),
    phenotype_id = col_character()
  )
) |>
  transmute(
    feature_id = phenotype_id,
    tss = start
  )

cell_types <- c("Plasma", "B_IN", "CD4_NC")
cell_type_levels <- unname(cell_type_lookup[cell_types])

overlap_levels <- c(
  "Unique to Pseudobulk LM",
  "Common to both",
  "Unique to Single-cell P-GLMM"
)
overlap_cols <- c(
  "Unique to Pseudobulk LM" = unname(method_col_lookup[["Pseudobulk LM"]]),
  "Common to both" = "#4477AA",
  "Unique to Single-cell P-GLMM" = unname(method_col_lookup[["Single-cell P-GLMM"]])
)
method_levels <- c("Pseudobulk LM", "Single-cell P-GLMM")

pb_manifest <- pb_data_files |>
  filter(
    dataset == "onek1k",
    model == "lm",
    pb_type == "logcounts",
    int_cov == "none",
    cell_frac == 1,
    indiv_frac == 1,
    min_n_cells < 0,
    count_frac == 1,
    cell_type %in% cell_types
  )

sc_manifest <- sc_data_files |>
  filter(
    dataset == "onek1k",
    model == "p_glmm_sc",
    data_type == "counts",
    cov_spec == "bulk_pca",
    k == "none",
    int_cov == "none",
    cell_frac == 1,
    indiv_frac == 1,
    min_n_cells < 0,
    count_frac == 1,
    ablate %in% c("none", ""),
    cell_type %in% cell_types
  )

if (nrow(pb_manifest) == 0L || nrow(sc_manifest) == 0L) {
  stop("Need both pseudobulk LM and single-cell P-GLMM region files.")
}

read_region_bh <- function(manifest, method) {
  manifest |>
    mutate(method = method) |>
    distinct(cell_type, method, region_file) |>
    mutate(gene_tbl = map(region_file, function(path) {
      read_tsv(path, show_col_types = FALSE) |>
        select(feature_id, pvalue)
    })) |>
    select(-region_file) |>
    unnest(gene_tbl) |>
    filter(!is.na(feature_id), is.finite(pvalue), pvalue > 0, pvalue <= 1) |>
    mutate(
      region_bh = p.adjust(pvalue, method = "BH"),
      .by = c(cell_type, method)
    ) |>
    select(cell_type, method, feature_id, region_pvalue = pvalue, region_bh)
}

read_leads <- function(paths) {
  leads <- rbindlist(lapply(paths, function(path) {
    dt <- fread(
      path,
      select = c("feature_id", "snp_id", "chrom", "pos", "pvalue"),
      showProgress = FALSE
    )
    dt <- dt[!is.na(feature_id) & !is.na(pvalue) & pvalue > 0 & !is.na(pos)]
    if (nrow(dt) == 0L) {
      return(data.table(
        feature_id = character(),
        snp_id = character(),
        chrom = character(),
        pos = double(),
        pvalue = double()
      ))
    }
    dt[order(pvalue), .(
      snp_id = snp_id[1],
      chrom = chrom[1],
      pos = pos[1],
      pvalue = pvalue[1]
    ), by = feature_id]
  }))
  if (nrow(leads) == 0L) {
    return(tibble(
      feature_id = character(),
      snp_id = character(),
      chrom = character(),
      pos = double(),
      lead_pvalue = double()
    ))
  }
  as_tibble(
    leads[order(pvalue), .(
      snp_id = snp_id[1],
      chrom = chrom[1],
      pos = pos[1],
      lead_pvalue = pvalue[1]
    ), by = feature_id]
  )
}

method_leads <- function(manifest, method) {
  manifest |>
    summarise(
      variant_file = list(unique(variant_file)),
      .by = cell_type
    ) |>
    mutate(lead = map(variant_file, read_leads)) |>
    select(-variant_file) |>
    unnest(lead) |>
    mutate(method = method)
}

pb_genes <- read_region_bh(pb_manifest, "Pseudobulk LM")
sc_genes <- read_region_bh(sc_manifest, "Single-cell P-GLMM")

sig_genes <- bind_rows(pb_genes, sc_genes) |>
  filter(region_bh < 0.05)

if (nrow(sig_genes) == 0L) {
  stop("No significant eGenes found.")
}

pb_leads <- method_leads(pb_manifest, "Pseudobulk LM")
sc_leads <- method_leads(sc_manifest, "Single-cell P-GLMM")

sig_with_leads <- sig_genes |>
  left_join(
    bind_rows(pb_leads, sc_leads),
    by = c("cell_type", "method", "feature_id")
  ) |>
  left_join(annot, by = "feature_id") |>
  mutate(tss_distance = abs(pos - tss))

egene_status <- sig_genes |>
  mutate(present = TRUE) |>
  pivot_wider(
    id_cols = c(cell_type, feature_id),
    names_from = method,
    values_from = present,
    values_fill = list(present = FALSE)
  ) |>
  mutate(
    category = case_when(
      `Pseudobulk LM` & `Single-cell P-GLMM` ~ "Common to both",
      `Pseudobulk LM` & !`Single-cell P-GLMM` ~ "Unique to Pseudobulk LM",
      !`Pseudobulk LM` & `Single-cell P-GLMM` ~ "Unique to Single-cell P-GLMM"
    ),
    category = factor(category, levels = overlap_levels)
  ) |>
  select(cell_type, feature_id, category)

gene_props <- bind_rows(
  pb_manifest |> distinct(cell_type, gene_prop_file),
  sc_manifest |> distinct(cell_type, gene_prop_file)
) |>
  distinct(cell_type, gene_prop_file) |>
  mutate(props_tbl = map(gene_prop_file, read_tsv, show_col_types = FALSE)) |>
  select(-gene_prop_file) |>
  unnest(props_tbl) |>
  distinct(cell_type, feature_id, .keep_all = TRUE) |>
  select(cell_type, feature_id, pb_mean)

wide_leads <- sig_with_leads |>
  select(
    cell_type, method, feature_id,
    snp_id, chrom, pos, lead_pvalue, tss_distance,
    region_pvalue, region_bh
  ) |>
  pivot_wider(
    id_cols = c(cell_type, feature_id),
    names_from = method,
    values_from = c(
      snp_id, chrom, pos, lead_pvalue, tss_distance,
      region_pvalue, region_bh
    ),
    names_sep = "__"
  )

rename_method_cols <- function(df) {
  nm <- names(df)
  nm <- gsub("Pseudobulk LM", "lm", nm, fixed = TRUE)
  nm <- gsub("Single-cell P-GLMM", "p_glmm_sc", nm, fixed = TRUE)
  names(df) <- nm
  df
}

comparison <- egene_status |>
  left_join(rename_method_cols(wide_leads), by = c("cell_type", "feature_id")) |>
  left_join(gene_props, by = c("cell_type", "feature_id")) |>
  mutate(
    cell_type_label = factor(
      unname(cell_type_lookup[cell_type]),
      levels = cell_type_levels
    )
  ) |>
  arrange(cell_type, category, feature_id)

write_tsv(comparison, "pb-lm-vs-sc-pglmm-egenes.tsv")

distance_data <- sig_with_leads |>
  inner_join(egene_status, by = c("cell_type", "feature_id")) |>
  filter(is.finite(tss_distance)) |>
  mutate(
    cell_type_label = factor(
      unname(cell_type_lookup[cell_type]),
      levels = cell_type_levels
    ),
    method = factor(method, levels = method_levels),
    category = factor(category, levels = overlap_levels)
  )

distance_p <- distance_data |>
  ggplot(aes(category, tss_distance + 1, fill = method)) +
  geom_boxplot(
    outlier.size = 0.4,
    outlier.alpha = 0.3,
    linewidth = 0.4,
    position = position_dodge2(preserve = "single")
  ) +
  scale_y_log10(
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  ) +
  scale_fill_manual(values = method_col_lookup[method_levels], drop = FALSE) +
  facet_wrap(~cell_type_label, scales = "free_y") +
  labs(
    x = NULL,
    y = "Lead variant distance to TSS + 1 (bp)",
    fill = NULL
  ) +
  theme_jp() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(
  "pb-lm-vs-sc-pglmm-tss-distance-plot.pdf",
  distance_p,
  width = 12,
  height = 5
)

mean_data <- comparison |>
  filter(is.finite(pb_mean))

mean_p <- mean_data |>
  ggplot(aes(category, pb_mean, fill = category)) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.3, linewidth = 0.4) +
  scale_y_log10() +
  scale_fill_manual(values = overlap_cols, guide = "none") +
  facet_wrap(~cell_type_label, scales = "free_y") +
  labs(
    x = NULL,
    y = "Pseudobulk mean count"
  ) +
  theme_jp() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(
  "pb-lm-vs-sc-pglmm-pb-mean-plot.pdf",
  mean_p,
  width = 12,
  height = 5
)
