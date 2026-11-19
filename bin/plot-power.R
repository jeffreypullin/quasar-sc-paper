#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(data.table)
  library(ggplot2)
  library(forcats)
  library(patchwork)
  library(qvalue)
  library(purrr)
  library(stringr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

pb_data_files <- read_tsv(args[1], show_col_types = FALSE)
sc_data_files <- read_tsv(args[2], show_col_types = FALSE)
saigeqtl_data_files <- read_tsv(args[3], show_col_types = FALSE)

saigeqtl_n_sig_var_data <- saigeqtl_data_files |>
  filter(variant_file != "NA") |>
  summarise(n_sig_variant = sum(n_sig_variant), .by = cell_type)

saigeqtl_n_sig_gene_data <- saigeqtl_data_files |>
  filter(variant_file != "NA") |>
  rowwise() |>
  mutate(test = list(read_tsv(region_file, show_col_type = FALSE))) |>
  unnest(test) |>
  summarise(
    n_sig_gene = sum(p.adjust(ACAT_p, method = "BH") < 0.05, na.rm = TRUE),
    .by = cell_type
  )

sc_n_sig_gene_data <- sc_data_files |>
  filter(indiv_frac == 1) |>
  filter(min_n_cells < 0) |>
  filter(count_frac == 1) |>
  filter(
    if ("data_type" %in% names(sc_data_files)) {
      (model == "p_glmm_sc" & data_type == "counts") |
        (model == "lmm_sc" & data_type == "log_counts")
    } else {
      TRUE
    }
  ) |>
  rowwise() |>
  mutate(test = list(read_tsv(region_file, show_col_types = FALSE))) |>
  unnest(test) |>
  summarise(
    n_sig_gene = sum(p.adjust(pvalue, method = "BH") < 0.05, na.rm = TRUE),
    .by = c(cell_type, cov_spec, cell_frac, model)
  )

pb_gene_files <- pb_data_files |>
  filter(indiv_frac == 1, min_n_cells < 0, count_frac == 1, int_cov == "none")

if (nrow(pb_gene_files) == 0) {
  pb_n_sig_gene_data <- tibble(
    cell_type = character(),
    model = character(),
    cell_frac = double(),
    n_sig_gene = integer()
  )
} else {
  gene_prop_data <- pb_gene_files |>
    rowwise() |>
    mutate(gene_data = list(read_tsv(gene_prop_file, show_col_types = FALSE))) |>
    unnest(gene_data) |>
    distinct(feature_id, cell_type, model, .keep_all = TRUE)

  pb_n_sig_gene_data <- pb_gene_files |>
    rowwise() |>
    mutate(test = list(read_tsv(region_file, show_col_types = FALSE))) |>
    unnest(test) |>
    select(pvalue, model, cell_type, feature_id, cell_frac) |>
    left_join(
      gene_prop_data |>
        select(model, cell_type, feature_id, pb_mean),
      by = c("cell_type", "feature_id", "model")
    ) |>
    summarise(
      n_sig_gene = sum(p.adjust(pvalue[!is.na(pvalue)], method = "BH") < 0.05, na.rm = TRUE),
      .by = c(cell_type, model, cell_frac)
    )
}

sc_n_sig_var_data <- sc_data_files |>
  filter(indiv_frac == 1) |>
  filter(min_n_cells < 0) |>
  filter(count_frac == 1) |>
  filter(
    if ("data_type" %in% names(sc_data_files)) {
      (model == "p_glmm_sc" & data_type == "counts") |
        (model == "lmm_sc" & data_type == "log_counts")
    } else {
      TRUE
    }
  ) |>
  rowwise() |>
  mutate(n_sig_variant = list(read_tsv(power_file, show_col_types = FALSE))) |>
  unnest(n_sig_variant) |>
  summarise(n_sig_variant = sum(n_sig_variant), .by = c(cell_type, cell_frac, cov_spec, model))

pb_var_files <- pb_data_files |>
  filter(indiv_frac == 1, min_n_cells < 0, count_frac == 1)

if (nrow(pb_var_files) == 0) {
  pb_n_sig_var_data <- tibble(
    cell_type = character(),
    model = character(),
    cell_frac = double(),
    n_sig_variant = integer()
  )
} else {
  pb_n_sig_var_data <- pb_var_files |>
    rowwise() |>
    mutate(n_sig_variant = list(read_tsv(power_file, show_col_types = FALSE))) |>
    unnest(n_sig_variant) |>
    summarise(n_sig_variant = sum(n_sig_variant), .by = c(cell_type, model, cell_frac))
}

variant_plot_data <- bind_rows(
  pb_n_sig_var_data |>
    filter(cell_frac == 1) |>
    summarise(n_sig_variant = sum(n_sig_variant), .by = c(cell_type, model)) |>
    mutate(type = paste0("pb-", model)) |>
    select(-model),
  sc_n_sig_var_data |>
    filter(cell_frac == 1) |>
    filter(cov_spec == "bulk_pca") |>
    summarise(n_sig_variant = sum(n_sig_variant), .by = c(cell_type, model)) |>
    mutate(type = paste0("sc-", model)),
  saigeqtl_n_sig_var_data |>
    mutate(type = "saigeqtl")
) |>
  filter(cell_type %in% c("Plasma", "B_IN", "CD4_NC")) |>
  filter(type %in% names(method_lookup)) |>
  mutate(
    cell_type = cell_type_lookup[cell_type],
    type = method_lookup[type]
  ) |>
  mutate(
    cell_type = fct_reorder(factor(cell_type), n_sig_variant, .fun = max),
    type = fct_reorder(factor(type), n_sig_variant, .fun = sum)
  )

method_levels <- levels(variant_plot_data$type)

gene_plot_data <- bind_rows(
  pb_n_sig_gene_data |>
    filter(cell_frac == 1) |>
    mutate(type = paste0("pb-", model)),
  sc_n_sig_gene_data |>
    filter(cell_frac == 1) |>
    filter(cov_spec == "bulk_pca") |>
    mutate(type = paste0("sc-", model)),
  saigeqtl_n_sig_gene_data |>
    mutate(type = "saigeqtl")
) |>
  filter(cell_type %in% c("Plasma", "B_IN", "CD4_NC")) |>
  filter(type %in% names(method_lookup)) |>
  mutate(
    cell_type = cell_type_lookup[cell_type],
    type = factor(method_lookup[type], levels = method_levels)
  ) |>
  mutate(cell_type = fct_reorder(factor(cell_type), n_sig_gene, .fun = max))

p_variant <- variant_plot_data |>
  ggplot(aes(cell_type, n_sig_variant, fill = type)) +
  geom_col(position = "dodge2") +
  scale_fill_manual(
    values = method_col_lookup,
    drop = FALSE,
    guide = guide_legend(reverse = TRUE)
  ) +
  labs(
    y = "Number of significant variants",
    x = "Cell type",
    fill = "Method"
  ) +
  coord_flip() +
  theme_jp_vgrid()

p_gene <- gene_plot_data |>
  ggplot(aes(cell_type, n_sig_gene, fill = type)) +
  geom_col(position = "dodge2") +
  scale_fill_manual(values = method_col_lookup, drop = FALSE) +
  labs(
    y = "Number of significant genes",
    x = NULL,
    fill = "Method"
  ) +
  coord_flip() +
  theme_jp_vgrid() +
  guides(fill = "none")

p <- p_variant + p_gene +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(
  "power-method-comparison-plot.pdf",
  plot = p,
  width = 16,
  height = 10
)
#frac_p <- sc_n_sig_var_data |>
#  summarise(
#    n_sig_variant = sum(n_sig_variant),
#    .by = c(cell_type, cov_spec, cell_frac)
#  ) |>
#  filter(cell_type == "CD4_NC") |>
#  filter(cov_spec == "bulk_pca") |>
#  ggplot(aes(cell_frac, n_sig_variant)) +
#  geom_point() +
#  geom_smooth(method = "loess", se = FALSE)

ggsave(
  "power-frac-plot.pdf",
  plot = p,
  width = 12,
  height = 8
)
#non_zero_p <- sc_n_sig_var_data |>
#  filter(cell_frac == 1) |>
#  filter(cov_spec == "bulk_pca") |>
#  summarise(
#    n_sig_variant = sum(n_sig_variant),
#    .by = c(cell_type, nz_bin)
#  ) |>
#  ggplot(aes(nz_bin, n_sig_variant)) +
#  geom_col(position = "dodge2") +
#  facet_wrap(~ cell_type)

ggsave(
  "power-sc-non-zero-frac-plot.pdf",
  plot = p,
  width = 12,
  height = 8
)

example_cell_types <- c("Plasma", "B_IN", "CD4_NC")
example_cell_levels <- rev(unname(cell_type_lookup[c("CD4_NC", "B_IN", "Plasma")]))
onek1k_cell_types <- c(
  "CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_ET", "CD8_NC", "CD8_S100B",
  "NK", "NK_R", "Plasma", "B_Mem", "B_IN", "Mono_C", "Mono_NC", "DC"
)

sc_pglmm_files <- sc_data_files |>
  filter(cell_frac == 1, indiv_frac == 1, min_n_cells < 0, count_frac == 1) |>
  filter(cov_spec == "bulk_pca", int_cov == "none", as.character(k) == "none") |>
  filter(model == "p_glmm_sc") |>
  filter(if ("data_type" %in% names(sc_data_files)) data_type == "counts" else TRUE) |>
  mutate(ablate = if ("ablate" %in% names(sc_data_files)) coalesce(as.character(ablate), "none") else "none") |>
  filter(ablate %in% c("none", ""))

sc_pglmm_variants <- sc_pglmm_files |>
  rowwise() |>
  mutate(n_sig_variant = list(read_tsv(power_file, show_col_types = FALSE))) |>
  unnest(n_sig_variant) |>
  summarise(n_sig_variant = sum(n_sig_variant), .by = cell_type) |>
  mutate(type = "sc-p_glmm_sc")

sc_pglmm_genes <- sc_pglmm_files |>
  rowwise() |>
  mutate(test = list(read_tsv(region_file, show_col_types = FALSE))) |>
  unnest(test) |>
  summarise(
    n_sig_gene = sum(p.adjust(pvalue, method = "BH") < 0.05, na.rm = TRUE),
    .by = cell_type
  ) |>
  mutate(type = "sc-p_glmm_sc")

pb_lm_files <- pb_data_files |>
  filter(
    dataset == "onek1k",
    cell_type %in% example_cell_types,
    model == "lm",
    pb_type == "logcounts",
    int_cov == "none",
    cell_frac == 1,
    indiv_frac == 1,
    min_n_cells < 0,
    count_frac == 1
  )

pb_lm_variants <- pb_lm_files |>
  rowwise() |>
  mutate(n_sig_variant = list(read_tsv(power_file, show_col_types = FALSE))) |>
  unnest(n_sig_variant) |>
  summarise(n_sig_variant = sum(n_sig_variant), .by = cell_type) |>
  mutate(type = "pb-lm")

pb_lm_genes <- pb_lm_files |>
  rowwise() |>
  mutate(test = list(read_tsv(region_file, show_col_types = FALSE))) |>
  unnest(test) |>
  summarise(
    n_sig_gene = sum(p.adjust(pvalue, method = "BH") < 0.05, na.rm = TRUE),
    .by = cell_type
  ) |>
  mutate(type = "pb-lm")

saige_vs_quasar <- function(saige, quasar, pseudobulk, value) {
  bind_rows(
    saige |> mutate(type = "saigeqtl"),
    quasar,
    pseudobulk
  ) |>
    filter(cell_type %in% example_cell_types) |>
    mutate(
      cell_type = factor(unname(cell_type_lookup[cell_type]), levels = example_cell_levels),
      type = factor(
        case_when(
          type == "saigeqtl" ~ "SAIGE-QTL",
          type == "sc-p_glmm_sc" ~ "Single-cell quasar P-GLMM",
          type == "pb-lm" ~ "Pseudobulk LM"
        ),
        levels = c("Pseudobulk LM", "SAIGE-QTL", "Single-cell quasar P-GLMM")
      ),
      value = .data[[value]]
    )
}

saige_vs_quasar_cols <- c(
  "Pseudobulk LM" = method_col_lookup[["Pseudobulk LM"]],
  "SAIGE-QTL" = method_col_lookup[["SAIGE-QTL"]],
  "Single-cell quasar P-GLMM" = method_col_lookup[["Single-cell P-GLMM"]]
)

p_saige_variant <- saige_vs_quasar(
  saigeqtl_n_sig_var_data,
  sc_pglmm_variants,
  pb_lm_variants,
  "n_sig_variant"
) |>
  ggplot(aes(cell_type, value, fill = type)) +
  geom_col(position = "dodge2") +
  scale_fill_manual(
    values = saige_vs_quasar_cols,
    guide = guide_legend(reverse = TRUE)
  ) +
  scale_y_continuous(
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  ) +
  labs(y = "Number of significant variants", x = "Cell type", fill = "Method") +
  coord_flip() +
  theme_jp_vgrid()

p_saige_gene <- saige_vs_quasar(
  saigeqtl_n_sig_gene_data,
  sc_pglmm_genes,
  pb_lm_genes,
  "n_sig_gene"
) |>
  ggplot(aes(cell_type, value, fill = type)) +
  geom_col(position = "dodge2") +
  scale_fill_manual(values = saige_vs_quasar_cols) +
  labs(y = "Number of significant genes", x = NULL, fill = "Method") +
  coord_flip() +
  theme_jp_vgrid() +
  guides(fill = "none")

p_saige_example <- p_saige_variant + p_saige_gene +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

ggsave(
  "power-saigeqtl-vs-quasar-example-celltypes-plot.pdf",
  plot = p_saige_example,
  width = 12,
  height = 5
)

saveRDS(p_saige_example, "power-saigeqtl-vs-quasar-example-celltypes-plot.rds")

sc_lmm_files <- sc_data_files |>
  filter(cell_type %in% onek1k_cell_types) |>
  filter(cell_frac == 1, indiv_frac == 1, min_n_cells < 0, count_frac == 1) |>
  filter(cov_spec == "bulk_pca", int_cov == "none", as.character(k) == "none") |>
  filter(model == "lmm_sc") |>
  filter(if ("data_type" %in% names(sc_data_files)) data_type == "log_counts" else TRUE) |>
  mutate(ablate = if ("ablate" %in% names(sc_data_files)) coalesce(as.character(ablate), "none") else "none") |>
  filter(ablate %in% c("none", ""))

sc_lmm_variants <- sc_lmm_files |>
  rowwise() |>
  mutate(n_sig_variant = list(read_tsv(power_file, show_col_types = FALSE))) |>
  unnest(n_sig_variant) |>
  summarise(n_sig_variant = sum(n_sig_variant), .by = cell_type) |>
  mutate(type = "sc-lmm_sc")

sc_lmm_genes <- sc_lmm_files |>
  rowwise() |>
  mutate(test = list(read_tsv(region_file, show_col_types = FALSE))) |>
  unnest(test) |>
  summarise(
    n_sig_gene = sum(p.adjust(pvalue, method = "BH") < 0.05, na.rm = TRUE),
    .by = cell_type
  ) |>
  mutate(type = "sc-lmm_sc")

discovery_cell_types <- sc_pglmm_genes |>
  filter(cell_type != "B_all", cell_type %in% names(cell_type_lookup)) |>
  arrange(n_sig_gene) |>
  pull(cell_type)

pglmm_lmm_cell_levels <- unname(cell_type_lookup[discovery_cell_types])

pglmm_vs_lmm <- function(pglmm, lmm, value) {
  bind_rows(pglmm, lmm) |>
    filter(cell_type %in% discovery_cell_types) |>
    mutate(
      cell_type = factor(unname(cell_type_lookup[cell_type]), levels = pglmm_lmm_cell_levels),
      type = factor(
        unname(method_lookup[type]),
        levels = c("Single-cell P-GLMM", "Single-cell LMM")
      ),
      value = .data[[value]]
    )
}

pglmm_lmm_cols <- c(
  "Single-cell P-GLMM" = method_col_lookup[["Single-cell P-GLMM"]],
  "Single-cell LMM" = method_col_lookup[["Single-cell LMM"]]
)

p_lmm_variant <- pglmm_vs_lmm(sc_pglmm_variants, sc_lmm_variants, "n_sig_variant") |>
  ggplot(aes(cell_type, value, fill = type)) +
  geom_col(position = "dodge2") +
  scale_fill_manual(values = pglmm_lmm_cols, drop = FALSE) +
  scale_y_continuous(
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  ) +
  labs(y = "Number of significant variants", x = "Cell type", fill = "Method") +
  coord_flip() +
  theme_jp_vgrid()

p_lmm_gene <- pglmm_vs_lmm(sc_pglmm_genes, sc_lmm_genes, "n_sig_gene") |>
  ggplot(aes(cell_type, value, fill = type)) +
  geom_col(position = "dodge2") +
  scale_fill_manual(values = pglmm_lmm_cols, drop = FALSE) +
  labs(y = "Number of significant genes", x = NULL, fill = "Method") +
  coord_flip() +
  theme_jp_vgrid() +
  guides(fill = "none")

p_lmm <- p_lmm_variant + p_lmm_gene +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

ggsave(
  "power-p_glmm_sc-vs-lmm_sc-all-celltypes-plot.pdf",
  plot = p_lmm,
  width = 12,
  height = max(6, 0.45 * n_distinct(p_lmm_variant$data$cell_type) + 1.5),
  limitsize = FALSE
)
saveRDS(p_lmm, "power-p_glmm_sc-vs-lmm_sc-all-celltypes-plot.rds")

b_all_method_levels <- c(
  "P-GLMM (counts)",
  "LMM (log counts)",
  "LMM (SCT residuals)"
)

b_all_files <- sc_data_files |>
  filter(
    dataset == "onek1k",
    cell_type == "B_all",
    cell_frac == 1,
    indiv_frac == 1,
    min_n_cells < 0,
    count_frac == 1,
    cov_spec == "bulk_pca",
    int_cov == "none",
    as.character(k) == "none",
    (model == "p_glmm_sc" & data_type == "counts") |
      (model == "lmm_sc" & data_type == "log_counts") |
      (model == "lmm_sc" & data_type == "sct_residuals")
  ) |>
  mutate(
    ablate = if ("ablate" %in% names(sc_data_files)) coalesce(as.character(ablate), "none") else "none",
    method = case_when(
      model == "p_glmm_sc" ~ "P-GLMM (counts)",
      data_type == "log_counts" ~ "LMM (log counts)",
      data_type == "sct_residuals" ~ "LMM (SCT residuals)"
    )
  ) |>
  filter(ablate %in% c("none", ""))

missing_b_all_methods <- setdiff(b_all_method_levels, unique(b_all_files$method))
if (length(missing_b_all_methods) > 0) {
  stop("Missing B_all power results for: ", paste(missing_b_all_methods, collapse = ", "))
}

b_all_variants <- b_all_files |>
  rowwise() |>
  mutate(power = list(read_tsv(power_file, show_col_types = FALSE))) |>
  unnest(power) |>
  summarise(value = sum(n_sig_variant, na.rm = TRUE), .by = method)

b_all_genes <- b_all_files |>
  rowwise() |>
  mutate(region = list(read_tsv(region_file, show_col_types = FALSE))) |>
  unnest(region) |>
  summarise(
    value = sum(p.adjust(pvalue, method = "BH") < 0.05, na.rm = TRUE),
    .by = method
  )

b_all_cols <- c(
  "P-GLMM (counts)" = method_col_lookup[["Single-cell P-GLMM"]],
  "LMM (log counts)" = method_col_lookup[["Single-cell LMM"]],
  "LMM (SCT residuals)" = "#D55E00"
)

plot_b_all_power <- function(data, y_label) {
  data |>
    mutate(method = factor(method, levels = rev(b_all_method_levels))) |>
    ggplot(aes(method, value, fill = method)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = b_all_cols, drop = FALSE) +
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    labs(x = NULL, y = y_label) +
    coord_flip() +
    theme_jp_vgrid() +
    guides(fill = "none")
}

p_b_all <- plot_b_all_power(b_all_variants, "Number of significant variants") +
  plot_b_all_power(b_all_genes, "Number of significant genes")

ggsave(
  "power-b-all-expression-comparison-plot.pdf",
  plot = p_b_all,
  width = 12,
  height = 4
)

pb_main_files <- pb_data_files |>
  filter(dataset == "onek1k", cell_type %in% onek1k_cell_types) |>
  filter(cell_frac == 1, indiv_frac == 1, min_n_cells < 0, count_frac == 1, int_cov == "none") |>
  filter(
    (model == "lm" & pb_type == "logcounts") |
      (model == "nb_glm" & pb_type == "counts")
  )

count_main_variants <- function(manifest, type_prefix) {
  if (nrow(manifest) == 0) {
    return(tibble(
      cell_type = character(),
      model = character(),
      n_sig_variant = integer(),
      type = character()
    ))
  }
  manifest |>
    rowwise() |>
    mutate(n_sig_variant = list(read_tsv(power_file, show_col_types = FALSE))) |>
    unnest(n_sig_variant) |>
    summarise(n_sig_variant = sum(n_sig_variant), .by = c(cell_type, model)) |>
    mutate(type = paste0(type_prefix, model))
}

count_main_genes <- function(manifest, type_prefix) {
  if (nrow(manifest) == 0) {
    return(tibble(
      cell_type = character(),
      model = character(),
      n_sig_gene = integer(),
      type = character()
    ))
  }
  manifest |>
    rowwise() |>
    mutate(test = list(read_tsv(region_file, show_col_types = FALSE))) |>
    unnest(test) |>
    summarise(
      n_sig_gene = sum(p.adjust(pvalue, method = "BH") < 0.05, na.rm = TRUE),
      .by = c(cell_type, model)
    ) |>
    mutate(type = paste0(type_prefix, model))
}

sc_main_files <- sc_pglmm_files |>
  filter(dataset == "onek1k", cell_type %in% onek1k_cell_types)

main_variants <- bind_rows(
  count_main_variants(pb_main_files, "pb-"),
  count_main_variants(sc_main_files, "sc-")
)
main_genes <- bind_rows(
  count_main_genes(pb_main_files, "pb-"),
  count_main_genes(sc_main_files, "sc-")
)

main_method_keys <- c("pb-lm", "pb-nb_glm", "sc-p_glmm_sc")
main_method_levels <- unname(method_lookup[main_method_keys])
main_cell_order <- main_genes |>
  filter(type == "sc-p_glmm_sc") |>
  arrange(n_sig_gene) |>
  pull(cell_type)
main_cell_order <- c(main_cell_order, setdiff(unique(main_genes$cell_type), main_cell_order))
main_cell_levels <- unname(cell_type_lookup[main_cell_order])

label_main <- function(data, value) {
  data |>
    filter(type %in% main_method_keys) |>
    mutate(
      cell_type = factor(unname(cell_type_lookup[cell_type]), levels = main_cell_levels),
      type = factor(unname(method_lookup[type]), levels = main_method_levels),
      value = .data[[value]]
    )
}

main_cols <- method_col_lookup[main_method_levels]

p_main_variant <- label_main(main_variants, "n_sig_variant") |>
  ggplot(aes(cell_type, value, fill = type)) +
  geom_col(position = "dodge2") +
  scale_fill_manual(values = main_cols, drop = FALSE) +
  scale_y_continuous(
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  ) +
  labs(y = "Number of significant variants", x = "Cell type", fill = "Method") +
  coord_flip() +
  theme_jp_vgrid()

p_main_gene <- label_main(main_genes, "n_sig_gene") |>
  ggplot(aes(cell_type, value, fill = type)) +
  geom_col(position = "dodge2") +
  scale_fill_manual(values = main_cols, drop = FALSE) +
  labs(y = "Number of significant genes", x = NULL, fill = "Method") +
  coord_flip() +
  theme_jp_vgrid() +
  guides(fill = "none")

p_main <- p_main_variant + p_main_gene +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

ggsave(
  "power-lm-nb_glm-p_glmm_sc-lmm_sc-onek1k-plot.pdf",
  plot = p_main,
  width = 14,
  height = max(6, 0.45 * n_distinct(p_main_variant$data$cell_type) + 1.5),
  limitsize = FALSE
)
saveRDS(p_main, "power-lm-nb_glm-p_glmm_sc-lmm_sc-onek1k-plot.rds")

weight_method_keys <- c("pb-lm_cw", "pb-lm", "sc-p_glmm_sc")
weight_method_levels <- unname(method_lookup[weight_method_keys])

pb_weight_files <- pb_data_files |>
  filter(dataset == "onek1k", cell_type %in% onek1k_cell_types) |>
  filter(cell_frac == 1, indiv_frac == 1, min_n_cells < 0, count_frac == 1, int_cov == "none") |>
  filter(model %in% c("lm", "lm_cw"), pb_type == "logcounts")

weight_variants <- bind_rows(
  count_main_variants(pb_weight_files, "pb-"),
  count_main_variants(
    sc_pglmm_files |> filter(dataset == "onek1k", cell_type %in% onek1k_cell_types),
    "sc-"
  )
)
weight_genes <- bind_rows(
  count_main_genes(pb_weight_files, "pb-"),
  count_main_genes(
    sc_pglmm_files |> filter(dataset == "onek1k", cell_type %in% onek1k_cell_types),
    "sc-"
  )
)

missing_weight_methods <- setdiff(weight_method_keys, unique(weight_genes$type))
if (length(missing_weight_methods) > 0) {
  stop("Missing weighted-LM power results for: ", paste(missing_weight_methods, collapse = ", "))
}

label_weight <- function(data, value) {
  labelled <- data |>
    filter(type %in% weight_method_keys) |>
    mutate(
      cell_type = unname(cell_type_lookup[cell_type]),
      type = unname(method_lookup[type]),
      value = .data[[value]]
    )
  cell_levels <- labelled |>
    summarise(max_value = max(value), .by = cell_type) |>
    arrange(max_value) |>
    pull(cell_type)
  labelled |>
    mutate(rank = rank(value, ties.method = "first"), .by = cell_type) |>
    mutate(
      cell_type = factor(cell_type, levels = cell_levels),
      type = factor(type, levels = weight_method_levels)
    )
}

p_weight_variant <- label_weight(weight_variants, "n_sig_variant") |>
  ggplot(aes(cell_type, value, fill = type, group = rank)) +
  geom_col(position = "dodge2") +
  scale_fill_manual(values = method_col_lookup[weight_method_levels], drop = FALSE) +
  scale_y_continuous(
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  ) +
  labs(y = "Number of significant variants", x = "Cell type", fill = "Method") +
  coord_flip() +
  theme_jp_vgrid()

p_weight_gene <- label_weight(weight_genes, "n_sig_gene") |>
  ggplot(aes(cell_type, value, fill = type, group = rank)) +
  geom_col(position = "dodge2") +
  scale_fill_manual(values = method_col_lookup[weight_method_levels], drop = FALSE) +
  labs(y = "Number of significant genes", x = NULL, fill = "Method") +
  coord_flip() +
  theme_jp_vgrid() +
  guides(fill = "none")

p_weight <- p_weight_variant + p_weight_gene +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

ggsave(
  "power-lm-lm_cw-p_glmm_sc-onek1k-plot.pdf",
  plot = p_weight,
  width = 14,
  height = max(6, 0.45 * n_distinct(p_weight_variant$data$cell_type) + 1.5),
  limitsize = FALSE
)

spa_levels <- c("Single-cell P-GLMM", "Single-cell P-GLMM (no SPA)")
spa_cols <- setNames(c(method_col_lookup[["Single-cell P-GLMM"]], "#999999"), spa_levels)

sc_spa_files <- sc_data_files |>
  filter(dataset == "onek1k", cell_type %in% onek1k_cell_types) |>
  filter(cell_frac == 1, indiv_frac == 1, min_n_cells < 0, count_frac == 1) |>
  filter(cov_spec == "bulk_pca", int_cov == "none", as.character(k) == "none") |>
  filter(model == "p_glmm_sc", data_type == "counts", ablate %in% c("none", "spa")) |>
  mutate(model = if_else(ablate == "spa", spa_levels[2], spa_levels[1]))

spa_genes <- count_main_genes(sc_spa_files, "")
spa_cell_levels <- spa_genes |>
  filter(model == spa_levels[1]) |>
  arrange(n_sig_gene) |>
  pull(cell_type)

label_spa <- function(data, value) {
  data |>
    mutate(
      cell_type = factor(unname(cell_type_lookup[cell_type]), levels = unname(cell_type_lookup[spa_cell_levels])),
      type = factor(model, levels = spa_levels),
      value = .data[[value]]
    )
}

p_spa_variant <- label_spa(count_main_variants(sc_spa_files, ""), "n_sig_variant") |>
  ggplot(aes(cell_type, value, fill = type)) +
  geom_col(position = "dodge2") +
  scale_fill_manual(values = spa_cols, drop = FALSE) +
  scale_y_continuous(
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  ) +
  labs(y = "Number of significant variants", x = "Cell type", fill = "Method") +
  coord_flip() +
  theme_jp_vgrid()

p_spa_gene <- label_spa(spa_genes, "n_sig_gene") |>
  ggplot(aes(cell_type, value, fill = type)) +
  geom_col(position = "dodge2") +
  scale_fill_manual(values = spa_cols, drop = FALSE) +
  labs(y = "Number of significant genes", x = NULL, fill = "Method") +
  coord_flip() +
  theme_jp_vgrid() +
  guides(fill = "none")

p_spa <- p_spa_variant + p_spa_gene +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

ggsave(
  "power-p_glmm_sc-spa-vs-no-spa-onek1k-plot.pdf",
  plot = p_spa,
  width = 14,
  height = max(6, 0.45 * n_distinct(p_spa_variant$data$cell_type) + 1.5),
  limitsize = FALSE
)


