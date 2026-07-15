#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)

sc_covs_file <- args[1]
sc_logcounts_file <- args[2]
genotype_file <- args[3]

gene_id <- "ENSG00000143297"
variant_id <- "1:157530620T-C"

sc_covs <- read_tsv(sc_covs_file, show_col_types = FALSE)
sc_logcounts <- read_tsv(sc_logcounts_file, show_col_types = FALSE)
genotypes <- read_tsv(genotype_file, show_col_types = FALSE)

plot_df <- sc_covs |>
  left_join(
    sc_logcounts |>
      mutate(log_expr = .data[[gene_id]]) |>
      select(sample_id, cell_id, log_expr),
    by = c("sample_id", "cell_id")
  ) |>
  left_join(
    genotypes |>
      mutate(dosage = .data[[variant_id]]) |>
      select(sample_id, dosage),
    by = "sample_id"
  ) |>
  group_by(sample_id) |>
  mutate(pseudotime_w = pseudotime - mean(pseudotime, na.rm = TRUE)) |>
  ungroup() |>
  mutate(geno = factor(dosage)) |>
  select(sample_id, cell_id, log_expr, dosage, geno, pseudotime, pseudotime_w) |>
  filter(!is.na(pseudotime), !is.na(log_expr), !is.na(geno))

p1 <- plot_df |>
  ggplot(aes(x = pseudotime, y = log_expr)) +
  geom_bin2d(bins = 50, alpha = 0.85) +
  scale_fill_viridis_c(option = "A", trans = "log10", name = "Cells") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(
    x = "Pseudotime",
    y = "log1p gene expression"
  ) +
  facet_wrap(~geno) +
  theme_bw()

n_bins <- 5
use_quantile <- TRUE

pt_vals <- plot_df$pseudotime
if (length(pt_vals) == 0 || !all(is.finite(range(pt_vals)))) {
  stop("No finite pseudotime values available for binning.")
}
pt_range <- range(pt_vals)
if (diff(pt_range) == 0) {
  stop("Pseudotime is constant; cannot define pseudotime bins.")
}

if (use_quantile) {
  pt_bin_ranges <- plot_df |>
    mutate(pt_bin = ntile(pseudotime, n_bins)) |>
    summarise(
      pt_lo = min(pseudotime),
      pt_hi = max(pseudotime),
      .by = pt_bin
    ) |>
    arrange(pt_bin) |>
    mutate(pt_bin_label = sprintf("%.2f–%.2f", pt_lo, pt_hi))

  plot_df_binned <- plot_df |>
    mutate(pt_bin = ntile(pseudotime, n_bins)) |>
    left_join(
      pt_bin_ranges |> select(pt_bin, pt_bin_label),
      by = "pt_bin"
    ) |>
    mutate(
      pt_bin_label = factor(
        pt_bin_label,
        levels = pt_bin_ranges$pt_bin_label
      )
    )
} else {
  pt_breaks <- seq(pt_range[1], pt_range[2], length.out = n_bins + 1)
  pt_bin_labels <- sprintf(
    "%.2f–%.2f",
    pt_breaks[-length(pt_breaks)],
    pt_breaks[-1]
  )

  plot_df_binned <- plot_df |>
    mutate(
      pt_bin = cut(
        pseudotime,
        breaks = pt_breaks,
        include.lowest = TRUE
      ),
      pt_bin_label = factor(
        sprintf(
          "%.2f–%.2f",
          pt_breaks[as.integer(pt_bin)],
          pt_breaks[as.integer(pt_bin) + 1]
        ),
        levels = pt_bin_labels
      )
    )
}

pt_bin_facet_labels <- plot_df_binned |>
  filter(!is.na(pt_bin_label)) |>
  count(pt_bin_label, name = "n_cells") |>
  mutate(
    facet_label = sprintf(
      "%s\n(n = %s cells)",
      pt_bin_label,
      format(n_cells, big.mark = ",")
    )
  )

bin_plot_df <- plot_df_binned |>
  left_join(
    pt_bin_facet_labels |> select(pt_bin_label, facet_label),
    by = "pt_bin_label"
  ) |>
  mutate(
    pt_bin_label = factor(
      facet_label,
      levels = pt_bin_facet_labels$facet_label
    )
  ) |>
  select(-facet_label) |>
  summarise(
    mean_expr = mean(log_expr, na.rm = TRUE),
    .by = c(sample_id, dosage, pt_bin_label)
  ) |>
  filter(!is.na(pt_bin_label), !is.na(mean_expr), !is.na(dosage))

p2 <- bin_plot_df |>
  ggplot(aes(x = dosage, y = mean_expr)) +
  geom_point(
    alpha = 0.35,
    size = 1.2,
    colour = "grey40",
    position = position_jitter(width = 0.06, height = 0)
  ) +
  geom_boxplot(
    aes(group = dosage),
    width = 0.5,
    fill = "white",
    colour = "grey20",
    linewidth = 0.6,
    alpha = 0.9,
    outlier.shape = NA
  ) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1, colour = "black") +
  scale_x_continuous(breaks = sort(unique(bin_plot_df$dosage))) +
  coord_cartesian(ylim = c(0, 0.75)) +
  facet_wrap(~pt_bin_label, nrow = 1) +
  labs(
    x = "Genotype dosage",
    y = "Mean log1p expression"
  ) +
  theme_bw()

fig_title <- sprintf("%s\n%s", gene_id, variant_id)

ggsave(
  "sc-int-figures.pdf",
  p2,
  width = 14,
  height = 14
)
