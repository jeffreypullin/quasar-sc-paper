#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(tidyr)
  library(forcats)
})

args <- commandArgs(trailingOnly = TRUE)

pb_data_files <- read_tsv(args[1], show_col_types = FALSE)
sc_data_files <- read_tsv(args[2], show_col_types = FALSE)

load_clumped <- function(data_files, method_col) {
  data_files |>
    filter(cell_frac == 1, indiv_frac == 1) |>
    rowwise() |>
    mutate(
      clump_tbl = list(
        tryCatch(
          read_tsv(clumped_file, show_col_types = FALSE),
          error = function(e) tibble()
        )
      )
    ) |>
    unnest(clump_tbl) |>
    ungroup() |>
    mutate(method = .data[[method_col]]) |>
    select(method, cell_type, feature_id, ID)
}

pb_clumps <- pb_data_files |>
  filter(cell_frac == 1, indiv_frac == 1) |>
  rowwise() |>
  mutate(clump_tbl = list(read_tsv(clumped_file, show_col_types = FALSE))) |>
  unnest(clump_tbl) |>
  ungroup() |>
  select(method = model, cell_type, feature_id, ID)

sc_clumps <- sc_data_files |>
  filter(cell_frac == 1, indiv_frac == 1) |>
  rowwise() |>
  mutate(clump_tbl = list(read_tsv(clumped_file, show_col_types = FALSE))) |>
  unnest(clump_tbl) |>
  ungroup() |>
  select(method = int_cov, cell_type, feature_id, ID)

all_clumps <- bind_rows(
  pb_clumps |> mutate(data_type = "pseudobulk"),
  sc_clumps |> mutate(data_type = "single-cell")
)

all_clumps <- all_clumps |>
  mutate(cell_type = str_replace_all(cell_type, "-", " ")) |>
  filter(cell_type %in% c("Plasma", "B IN", "CD4 NC")) |>
  filter(!(data_type == "pseudobulk" & !method %in% c("nb_glm", "lm"))) |>
  mutate(
    type = if_else(data_type == "pseudobulk", paste0("pb-", method), "single-cell")
  )

signal_counts <- all_clumps |>
  summarise(
    n_signals = n_distinct(ID),
    .by = c(type, cell_type, feature_id)
  )

method_signal_counts <- signal_counts |>
  summarise(
    n_signals = sum(n_signals),
    .by = c(type, cell_type)
  ) |>
  mutate(cell_type = fct_reorder(cell_type, n_signals, .fun = sum))

p_total <- method_signal_counts |>
  ggplot(aes(x = cell_type, y = n_signals, fill = type)) +
  geom_col(position = "dodge2") +
  coord_flip() +
  labs(
    x = NULL,
    y = "Total independent eQTLs",
    fill = "Method",
    title = "Total independent eQTL signals per cell type"
  ) +
  theme_bw()

ggsave(
  "clumped-n-independent-eqtls-plot.pdf",
  plot = p_total,
  width = 10,
  height = 5
)

per_gene_signal_dist <- signal_counts |>
  mutate(
    n_signals_capped = pmin(n_signals, 5L),
    n_signals_label = if_else(
      n_signals_capped == 5L,
      "5+",
      as.character(n_signals_capped)
    )
  ) |>
  count(type, cell_type, n_signals_label) |>
  mutate(
    n_signals_label = factor(n_signals_label, levels = c("1", "2", "3", "4", "5+"))
  )

p_dist <- per_gene_signal_dist |>
  ggplot(aes(x = n_signals_label, y = n, fill = type)) +
  geom_col(position = "dodge2") +
  facet_wrap(~ cell_type, scales = "free_y") +
  labs(
    x = "Number of independent signals",
    y = "Number of genes",
    fill = "Method",
    title = "Distribution of independent eQTL signal counts per gene"
  ) +
  theme_bw() +
  theme(strip.text = element_text(size = 9))

ggsave(
  "clumped-n-signals-plot.pdf",
  plot = p_dist,
  width = 12,
  height = 5
)

multi_signal <- signal_counts |>
  filter(n_signals >= 2) |>
  summarise(
    n_multi_signal_genes = n_distinct(feature_id),
    .by = c(type, cell_type)
  ) |>
  mutate(cell_type = fct_reorder(cell_type, n_multi_signal_genes, .fun = sum))

p_multi <- multi_signal |>
  ggplot(aes(x = cell_type, y = n_multi_signal_genes, fill = type)) +
  geom_col(position = "dodge2") +
  coord_flip() +
  labs(
    x = NULL,
    y = "Genes with \u22652 independent signals",
    fill = "Method",
    title = "Genes with multiple independent eQTL signals"
  ) +
  theme_bw()

ggsave(
  "clumped-multi-signal-plot.pdf",
  plot = p_multi,
  width = 10,
  height = 5
)
