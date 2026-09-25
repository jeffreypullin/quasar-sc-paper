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
  library(ggh4x)
})

args <- commandArgs(trailingOnly = TRUE)

pb_data_files <- read_tsv(args[1], show_col_types = FALSE) |>
  filter(int_cov == "none")

compute_qq_data_pb <- function(variant_files, prop_file) {

  variant_data <- bind_rows(!!!map(
    variant_files,
    function(x) {
      read_tsv(x, show_col_types = FALSE) |>
        select(feature_id, pvalue, maf)
    })
  )

  prop_data <- read_tsv(prop_file, show_col_types = FALSE)

  pvalue_data <- left_join(
    variant_data,
    prop_data,
    by = "feature_id"
  ) |>
    filter(!is.na(pvalue))

  pvalue <- pvalue_data$pvalue
  n <- length(pvalue)
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
      .by = c(x_bin)
    ) |>
    mutate(log_x_bin_mid = seq(0.05, 5.95, by = 0.1)[as.numeric(x_bin)])
}

pb_overall_plot_data <- pb_data_files |>
  filter(model == "lm") |>
  filter(cell_frac == 1) |>
  filter(indiv_frac == 1) |>
  filter(count_frac == 1) |>
  summarise(
    file_list = list(variant_file),
    prop_file = first(prop_file),
    .by = c(dataset, cell_type)
  ) |>
  rowwise() |>
  mutate(qq_data = list(compute_qq_data_pb(file_list, prop_file))) |>
  ungroup() |>
  select(-file_list) |>
  unnest(cols = qq_data)

pb_overall_p <- pb_overall_plot_data |>
  ggplot(aes(log_x_bin_mid, log_y_pvalue,
             ymin = log_lower_ci, ymax = log_upper_ci,
             colour = cell_type)) +
  geom_point(alpha = 0.8) +
  geom_abline(linetype = "dashed") +
  geom_ribbon(linetype = 2, alpha = 0.1) +
  facet_wrap(~dataset + cell_type) +
  labs(
    x = "Expected -log10(p-value)",
    y = "Observed -log10(p-value)",
    colour = "Cell type"
  )

ggsave(
  "perm-pb-overall-plot.pdf",
  pb_overall_p,
  width = 12,
  height = 10
)

#pb_overall_int_plot_data <- pb_data_files |>
#  filter(cell_frac == 1) |>
#  filter(indiv_frac == 1) |>
#  summarise(file_list = list(variant_file), .by = cell_type) |>
#  rowwise() |>
#  mutate(qq_data = list(compute_qq_data_pb_int(file_list, "age"))) |>
#  ungroup() |>
#  select(-file_list) |>
#  unnest(cols = qq_data)

#pb_overall_int_p <- pb_overall_int_plot_data |>
#  ggplot(aes(log_x_bin_mid, log_y_pvalue,
#             ymin = log_lower_ci, ymax = log_upper_ci,
#             colour = cell_type)) +
#  geom_point(alpha = 0.8) +
#  geom_abline(linetype = "dashed") +
#  geom_ribbon(linetype = 2, alpha = 0.1) +
#  facet_wrap(~cell_type) +
#  labs(
#    x = "Expected -log10(p-value)",
#    y = "Observed -log10(p-value)",
#    colour = "Cell type"
#  )

ggsave(
  "perm-pb-overall-int-plot.pdf",
  pb_overall_p,
  width = 12,
  height = 10
)

#pb_cell_frac_plot_data <- pb_data_files |>
#  filter(indiv_frac == 1) |>
#  mutate(cell_frac = factor(cell_frac)) |>
#  summarise(
#    file_list = list(variant_file),
#    prop_file = first(prop_file),
#    .by = c(cell_frac, cell_type)
#  ) |>
#  rowwise() |>
#  mutate(qq_data = list(compute_qq_data_pb(file_list, prop_file))) |>
#  ungroup() |>
#  select(-file_list) |>
#  unnest(cols = qq_data)

#pb_cell_frac_p <- pb_cell_frac_plot_data |>
#  mutate(cell_type = factor(cell_type,
#                            levels = c("CD4_NC", "B_IN", "Plasma"))) |>
#  ggplot(aes(log_x_bin_mid, log_y_pvalue,
#             ymin = log_lower_ci, ymax = log_upper_ci)) +
#  geom_point(alpha = 0.8) +
#  geom_abline(linetype = "dashed") +
#  geom_ribbon(linetype = 2, alpha = 0.1) +
#  facet_grid2(vars(cell_type), vars(cell_frac),
#              scales = "free_y", independent = "y") +
#  labs(
#    x = "Expected -log10(p-value)",
#    y = "Observed -log10(p-value)"
#  )

ggsave(
  "perm-pb-cell-frac-plot.pdf",
  pb_overall_p,
  width = 14,
  height = 10
)