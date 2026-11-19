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
  library(scales)
  library(readxl)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

read_time <- function(path) {
  data <- read.delim(path)
  raw_str <- colnames(data)
  str <- str_extract(raw_str, "(?<=real\\.).*")
  as.numeric(str)
}

# multiple is how many times faster the fastest method is than this one.
write_time_tsv <- function(data, path) {
  per_cell <- data |>
    mutate(multiple = time / min(time), .by = cell_type) |>
    arrange(cell_type, time) |>
    select(cell_type, method = type, time_h = hours, multiple)

  overall <- per_cell |>
    summarise(time_h = sum(time_h), .by = method) |>
    mutate(
      cell_type = "Overall",
      multiple = time_h / min(time_h)
    ) |>
    arrange(time_h) |>
    select(cell_type, method, time_h, multiple)

  bind_rows(per_cell, overall) |>
    write_tsv(path)
}

args <- commandArgs(trailingOnly = TRUE)

pb_data_files <- read_tsv(args[1], show_col_types = FALSE)
sc_data_files <- read_tsv(args[2], show_col_types = FALSE)
saigeqtl_data_files <- read_tsv(args[3], show_col_types = FALSE)

saigeqtl_time_data <- saigeqtl_data_files |>
  filter(variant_file != "NA") |>
  summarise(
    step1_time = sum(step1_time),
    step2_time = sum(step2_time),
    step3_time = sum(step3_time),
    .by = cell_type) |>
  mutate(time = step1_time + step2_time + step3_time) |>
  mutate(cell_frac = 1, indiv_frac = 1)

sc_time_data <- sc_data_files |>
  filter(dataset == "onek1k") |>
  filter(indiv_frac == 1) |>
  filter(min_n_cells < 0) |>
  filter(count_frac == 1) |>
  filter(cov_spec == "bulk_pca") |>
  filter(int_cov == "none") |>
  filter(as.character(k) == "none") |>
  mutate(ablate = coalesce(ablate, "none")) |>
  filter(ablate %in% c("none", "donor-totals")) |>
  mutate(model = if_else(ablate == "donor-totals", paste0(model, "_ablate"), model)) |>
  mutate(time = map_dbl(time_file, read_time)) |>
  summarise(time = sum(time), .by = c(cell_type, cell_frac, indiv_frac, model))

pb_time_data <- pb_data_files |>
  filter(indiv_frac == 1) |>
  filter(min_n_cells < 0) |>
  filter(count_frac == 1) |>
  mutate(time = map_dbl(time_file, read_time)) |>
  summarise(time = sum(time), .by = c(cell_type, model, cell_frac, indiv_frac)) |>
  mutate(cell_frac = as.numeric(cell_frac), indiv_frac = as.numeric(indiv_frac))

plot_data <- bind_rows(
  saigeqtl_time_data |>
    mutate(type = "saigeqtl"),
  sc_time_data |>
    mutate(type = paste0("sc-", model)),
  pb_time_data |>
    mutate(type = paste0("pb-", model)) |>
    select(-model),
)

quasar_data <- plot_data |>
  filter(cell_frac == 1) |>
  filter(indiv_frac == 1) |>
  filter(cell_type %in% c("Plasma", "B_IN", "CD4_NC")) |>
  filter(type != "saigeqtl") |>
  filter(type != "pb-nb_glm") |>
  filter(type != "pb-lm_cw") |>
  filter(type %in% names(method_lookup)) |>
  mutate(
    cell_type = factor(
      cell_type_lookup[cell_type],
      levels = cell_type_lookup[c("CD4_NC", "B_IN", "Plasma")]
    ),
    type = method_lookup[type]
  ) |>
  mutate(
    cell_type = fct_reorder(cell_type, time, .fun = max),
    type = fct_reorder(factor(type), time, .fun = sum),
    hours = time / (60 * 60)
  )

write_time_tsv(quasar_data, "time-quasar-plot.tsv")

quasar_p <- quasar_data |>
  ggplot(aes(cell_type, hours, fill = type)) +
  geom_col(position = "dodge2") +
  coord_flip() +
  scale_y_continuous(breaks = breaks_width(1)) +
  labs(
    x = "Cell type",
    y = "Time (h)",
    fill = "Method"
  ) +
  scale_fill_manual(values = method_col_lookup) +
  theme_jp_vgrid() +
  theme(legend.position = "right")

ggsave(
  "time-quasar-plot.pdf",
  quasar_p,
  width = 14,
  height = 6
)

example_method_keys <- c(
  "pb-lm", "pb-p_glmm_id", "sc-p_glmm_sc", "sc-p_glmm_sc_ablate"
)
example_method_levels <- unname(method_lookup[example_method_keys])

example_data <- plot_data |>
  filter(cell_frac == 1) |>
  filter(indiv_frac == 1) |>
  filter(cell_type %in% c("Plasma", "B_IN", "CD4_NC")) |>
  filter(type %in% example_method_keys) |>
  mutate(
    cell_type = factor(
      cell_type_lookup[cell_type],
      levels = cell_type_lookup[c("CD4_NC", "B_IN", "Plasma")]
    ),
    type = factor(method_lookup[type], levels = example_method_levels),
    hours = time / (60 * 60)
  ) |>
  mutate(cell_type = fct_reorder(cell_type, time, .fun = max))

write_time_tsv(
  example_data,
  "time-lm-p_glmm_id-p_glmm_sc-example-celltypes-plot.tsv"
)

example_p <- example_data |>
  ggplot(aes(cell_type, hours, fill = type)) +
  geom_col(position = "dodge2") +
  coord_flip() +
  scale_y_continuous(breaks = breaks_width(1)) +
  labs(
    x = "Cell type",
    y = "Time (h)",
    fill = "Method"
  ) +
  scale_fill_manual(values = method_col_lookup, drop = FALSE) +
  theme_jp_vgrid() +
  theme(legend.position = "right")

ggsave(
  "time-lm-p_glmm_id-p_glmm_sc-example-celltypes-plot.pdf",
  example_p,
  width = 14,
  height = 6
)

comparison_data <- plot_data |>
  filter(cell_frac == 1) |>
  filter(indiv_frac == 1) |>
  filter(cell_type %in% c("Plasma", "B_IN", "CD4_NC")) |>
  filter(type %in% c("saigeqtl", "sc-p_glmm_sc", "sc-p_glmm_sc_ablate")) |>
  mutate(
    cell_type = factor(
      cell_type_lookup[cell_type],
      levels = cell_type_lookup[c("CD4_NC", "B_IN", "Plasma")]
    ),
    type = method_lookup[type]
  ) |>
  mutate(
    cell_type = fct_reorder(cell_type, time, .fun = max),
    type = fct_reorder(factor(type), time, .fun = sum),
    hours = time / (60 * 60)
  )

write_time_tsv(comparison_data, "time-method-comparison-plot.tsv")
write_time_tsv(comparison_data, "time-frac-plot.tsv")

comparison_p <- comparison_data |>
  ggplot(aes(cell_type, hours, fill = type)) +
  geom_col(position = "dodge2") +
  coord_flip() +
  labs(
    x = "Cell type",
    y = "Time (h)",
    fill = "Method"
  ) +
  scale_fill_manual(values = method_col_lookup) +
  theme_jp_vgrid()

ggsave(
  "time-method-comparison-plot.pdf",
  comparison_p,
  width = 16,
  height = 10
)

ggsave(
  "time-frac-plot.pdf",
  comparison_p,
  width = 12,
  height = 10
)

# SAIGE-QTL vs quasar single-cell P-GLMM fitted without donor totals.
saigeqtl_ablate_data <- plot_data |>
  filter(cell_frac == 1) |>
  filter(indiv_frac == 1) |>
  filter(cell_type %in% c("Plasma", "B_IN", "CD4_NC")) |>
  filter(type %in% c("saigeqtl", "sc-p_glmm_sc_ablate")) |>
  mutate(
    cell_type = factor(
      cell_type_lookup[cell_type],
      levels = cell_type_lookup[c("CD4_NC", "B_IN", "Plasma")]
    ),
    type = if_else(
      type == "saigeqtl", "SAIGE-QTL", "Single-cell quasar P-GLMM"
    )
  ) |>
  mutate(
    cell_type = fct_reorder(cell_type, time, .fun = max),
    type = fct_reorder(factor(type), time, .fun = sum),
    hours = time / (60 * 60)
  )

write_time_tsv(saigeqtl_ablate_data, "time-saigeqtl-vs-quasar-example-celltypes-plot.tsv")

saigeqtl_ablate_p <- saigeqtl_ablate_data |>
  ggplot(aes(cell_type, hours, fill = type)) +
  geom_col(position = "dodge2") +
  coord_flip() +
  labs(
    x = "Cell type",
    y = "Time (h)",
    fill = "Method"
  ) +
  scale_fill_manual(
    values = c(
      "SAIGE-QTL" = method_col_lookup[["SAIGE-QTL"]],
      "Single-cell quasar P-GLMM" = method_col_lookup[["Single-cell P-GLMM"]]
    )
  ) +
  theme_jp_vgrid()

ggsave(
  "time-saigeqtl-vs-quasar-example-celltypes-plot.pdf",
  saigeqtl_ablate_p,
  width = 10,
  height = 6
)

saveRDS(saigeqtl_ablate_p, "time-saigeqtl-vs-quasar-example-celltypes-plot.rds")

ablated_pglmm <- plot_data |>
  filter(cell_frac == 1, indiv_frac == 1, type == "sc-p_glmm_sc_ablate") |>
  arrange(time)

pglmm_lmm_cell_levels <- unname(cell_type_lookup[ablated_pglmm$cell_type])

pglmm_lmm_data <- plot_data |>
  filter(cell_frac == 1, indiv_frac == 1) |>
  filter(cell_type %in% ablated_pglmm$cell_type) |>
  filter(type %in% c("sc-p_glmm_sc_ablate", "sc-lmm_sc"))

pglmm_lmm_plot_data <- pglmm_lmm_data |>
  mutate(
    cell_type = factor(unname(cell_type_lookup[cell_type]), levels = pglmm_lmm_cell_levels),
    type = factor(
      unname(method_lookup[type]),
      levels = c("Single-cell P-GLMM (no donor totals)", "Single-cell LMM")
    ),
    hours = time / (60 * 60)
  )

write_time_tsv(pglmm_lmm_plot_data, "time-p_glmm_sc-vs-lmm_sc-all-celltypes-plot.tsv")

pglmm_lmm_p <- pglmm_lmm_plot_data |>
  ggplot(aes(cell_type, hours, fill = type)) +
  geom_col(position = "dodge2") +
  coord_flip() +
  labs(
    x = "Cell type",
    y = "Time (h)",
    fill = "Method"
  ) +
  scale_fill_manual(values = method_col_lookup, drop = FALSE) +
  theme_jp_vgrid()

ggsave(
  "time-p_glmm_sc-vs-lmm_sc-all-celltypes-plot.pdf",
  pglmm_lmm_p,
  width = 10,
  height = max(6, 0.45 * n_distinct(pglmm_lmm_data$cell_type) + 1.5)
)
saveRDS(pglmm_lmm_p, "time-p_glmm_sc-vs-lmm_sc-all-celltypes-plot.rds")
