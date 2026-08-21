#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(data.table)
  library(ggplot2)
  library(purrr)
  library(stringr)
  library(tidyr)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)

sc_data_files <- read_tsv(args[1], show_col_types = FALSE) |>
  filter(k != "none")

compute_qq_data <- function(pvalue) {
  pvalue <- pvalue[!is.na(pvalue)]
  n <- length(pvalue)
  if (n == 0L) {
    return(tibble(
      log_x_bin_mid = double(),
      log_y_pvalue = double(),
      log_lower_ci = double(),
      log_upper_ci = double()
    ))
  }

  m <- (1:n) / (n + 1)
  z <- abs(qnorm(0.05 / 2))
  v <- (1:n) * (n - (1:n) + 1) / (n + 1)^2 / (n + 2)
  s <- sqrt(v)
  lower_ci <- m - z * s
  upper_ci <- m + z * s

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

stat_id_for_col <- function(col) {
  if (col == "pvalue") {
    return("pvalue")
  }
  if (col == "group_linear_pvalue") {
    return("linear")
  }
  if (col == "group_acat_pvalue") {
    return("acat")
  }
  q <- str_match(col, "_q([0-9]+)_pvalue$")[, 2]
  if (!is.na(q)) {
    return(paste0("q", q))
  }
  col
}

stat_label <- function(stat) {
  if (stat == "pvalue") {
    return("main-effect p-value")
  }
  if (stat == "linear") {
    return("group_linear p-value")
  }
  if (stat == "acat") {
    return("group_acat p-value")
  }
  if (startsWith(stat, "q")) {
    return(paste0(stat, " p-value"))
  }
  stat
}

stat_order <- function(stat) {
  qn <- suppressWarnings(as.integer(sub("^q", "", stat)))
  dplyr::case_when(
    stat == "pvalue" ~ 1L,
    stat == "linear" ~ 2L,
    stat == "acat" ~ 3L,
    !is.na(qn) ~ 10L + qn,
    TRUE ~ 99L
  )
}

collect_group_qq <- function(variant_files, prop_file, int_cov) {
  header <- names(fread(variant_files[[1]], nrows = 0, showProgress = FALSE))
  q_cols <- grep(
    paste0("^", int_cov, "_q[0-9]+_pvalue$"),
    header,
    value = TRUE
  )
  pvalue_cols <- intersect(
    c("pvalue", "group_linear_pvalue", "group_acat_pvalue", q_cols),
    header
  )
  if (length(pvalue_cols) == 0L) {
    return(tibble())
  }

  prop_data <- fread(
    prop_file,
    select = c("feature_id", "pb_non_zero_frac"),
    showProgress = FALSE
  )
  keep_genes <- prop_data[pb_non_zero_frac > 0.1, feature_id]

  variant_data <- rbindlist(lapply(variant_files, function(x) {
    dt <- fread(
      x,
      select = c("feature_id", "maf", pvalue_cols),
      showProgress = FALSE
    )
    dt[maf > 0.1 & feature_id %chin% keep_genes]
  }))

  bind_rows(lapply(pvalue_cols, function(col) {
    compute_qq_data(variant_data[[col]]) |>
      mutate(stat = stat_id_for_col(col))
  }))
}

plot_data <- sc_data_files |>
  filter(cell_frac == 1, indiv_frac == 1) |>
  summarise(
    file_list = list(variant_file),
    prop_file = first(prop_file),
    .by = c(cell_type, int_cov, k)
  ) |>
  rowwise() |>
  mutate(qq_data = list(collect_group_qq(file_list, prop_file, int_cov))) |>
  ungroup() |>
  select(-file_list) |>
  unnest(cols = qq_data) |>
  mutate(
    plot_id = if_else(str_starts(stat, "q"), "quantile", stat),
    quantile = if_else(str_starts(stat, "q"), stat, NA_character_),
    panel_label = paste0(cell_type, " (", int_cov, ")")
  )

make_qq_plot <- function(data, colour_var, xlab, ylab, colour_lab) {
  data |>
    ggplot(aes(
      log_x_bin_mid, log_y_pvalue,
      ymin = log_lower_ci, ymax = log_upper_ci,
      colour = .data[[colour_var]]
    )) +
    geom_point(alpha = 0.8) +
    geom_abline(linetype = "dashed") +
    geom_ribbon(linetype = 2, alpha = 0.1) +
    facet_wrap(vars(panel_label)) +
    labs(
      x = xlab,
      y = ylab,
      colour = colour_lab
    ) +
    theme_jp()
}

single_stats <- c("pvalue", "linear", "acat")
single_stats <- single_stats[single_stats %in% plot_data$plot_id]

for (s in single_stats) {
  p <- make_qq_plot(
    plot_data |> filter(plot_id == s),
    "cell_type",
    paste0("Expected -log10(", stat_label(s), ")"),
    paste0("Observed -log10(", stat_label(s), ")"),
    "Cell type"
  )
  ggsave(
    paste0("perm-sc-grouped-", s, "-plot.pdf"),
    p,
    width = 14,
    height = 10
  )
}

quantile_data <- plot_data |>
  filter(plot_id == "quantile") |>
  mutate(
    quantile = factor(
      quantile,
      levels = unique(quantile)[order(stat_order(unique(quantile)))]
    )
  )

if (nrow(quantile_data) > 0L) {
  quantile_p <- make_qq_plot(
    quantile_data,
    "quantile",
    "Expected -log10(quantile p-value)",
    "Observed -log10(quantile p-value)",
    "Quantile"
  )
  ggsave(
    "perm-sc-grouped-quantile-plot.pdf",
    quantile_p,
    width = 14,
    height = 10
  )
}
