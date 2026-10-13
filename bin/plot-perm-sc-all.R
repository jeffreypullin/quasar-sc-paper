#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(stringr)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

qq_data <- map(commandArgs(trailingOnly = TRUE), readRDS)

collapse_unique <- function(x) paste(unique(x), collapse = ", ")

cell_type_label <- cell_type_lookup
cell_type_label[c("B_all", "T_all")] <- c("All B", "All T")

plot_qq <- function(data, type) {
  int_cov <- collapse_unique(str_replace_all(data$int_cov, "_", " "))
  title <- if (type == "int") {
    paste0("Genotype × ", int_cov, " interaction")
  } else {
    paste0("Genotype main effect (genotype × ", int_cov, " model)")
  }

  data |>
    mutate(perm = recode(perm, Global = "Genotype")) |>
    ggplot(aes(log_x_bin_mid, log_y_pvalue,
               ymin = log_lower_ci, ymax = log_upper_ci,
               colour = perm, fill = perm)) +
    geom_point(alpha = 0.8) +
    geom_abline(linetype = "dashed") +
    geom_ribbon(linetype = 2, alpha = 0.1, colour = NA) +
    facet_wrap(~perm, labeller = as_labeller(\(x) paste(x, "permutation"))) +
    labs(
      title = title,
      subtitle = paste0(
        collapse_unique(cell_type_label[data$cell_type]), " cells  ·  ",
        collapse_unique(method_lookup[data$model])
      ),
      x = "Expected -log10(p-value)",
      y = "Observed -log10(p-value)"
    ) +
    theme_jp() +
    theme(legend.position = "none")
}

ggsave(
  "perm-sc-all-int-plot.pdf",
  plot_qq(map_dfr(qq_data, "int"), "int"),
  width = 12,
  height = 5
)

ggsave(
  "perm-sc-all-main-plot.pdf",
  plot_qq(map_dfr(qq_data, "main"), "main"),
  width = 12,
  height = 5
)
