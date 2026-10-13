#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tidyr)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)
pb_data_files <- read_tsv(args[1], show_col_types = FALSE)
sc_data_files <- read_tsv(args[2], show_col_types = FALSE)

# Targets at or above this offset keep the same donors as the equal-cell
# subsample and every cell from those donors.
donor_only_offset <- 1000000

n_cells_specs <- tribble(
  ~cell_type, ~n_cells,
  "Plasma", 5,
  "B_IN", 100,
  "CD4_NC", 300
)

n_cells_targets <- bind_rows(
  n_cells_specs |> transmute(cell_type, n_cells_target = n_cells),
  n_cells_specs |> transmute(
    cell_type,
    n_cells_target = donor_only_offset + n_cells
  )
)

keep_base <- function(data) {
  data |>
    filter(
      cell_type %in% n_cells_specs$cell_type,
      as.numeric(cell_frac) == 1,
      as.numeric(indiv_frac) == 1,
      as.numeric(count_frac) == 1,
      int_cov == "none"
    )
}

select_runs <- function(data) {
  data |>
    mutate(n_cells_target = as.numeric(n_cells_target)) |>
    filter(n_cells_target < 0) |>
    bind_rows(
      data |>
        mutate(n_cells_target = as.numeric(n_cells_target)) |>
        inner_join(n_cells_targets, by = c("cell_type", "n_cells_target"))
    )
}

pb_manifest <- pb_data_files |>
  keep_base() |>
  filter(model == "lm", pb_type == "logcounts") |>
  select_runs()

sc_manifest <- sc_data_files |>
  keep_base() |>
  filter(cov_spec == "bulk_pca", as.character(k) == "none") |>
  select_runs()

count_sig_genes <- function(manifest) {
  manifest |>
    mutate(gene_tbl = map(region_file, read_tsv, show_col_types = FALSE)) |>
    unnest(gene_tbl) |>
    summarise(
      n_sig_gene = sum(p.adjust(pvalue[!is.na(pvalue)], method = "BH") < 0.05, na.rm = TRUE),
      .by = c(cell_type, model, n_cells_target)
    )
}

sample_label <- function(n_cells_target) {
  case_when(
    n_cells_target < 0 ~ "All cells",
    n_cells_target >= donor_only_offset ~ "Same donors",
    TRUE ~ paste0(n_cells_target, " per donor")
  )
}

method_levels <- c(
  "Pseudobulk LM",
  "Single-cell LMM",
  "Single-cell P-GLMM"
)

plot_data <- bind_rows(
  count_sig_genes(pb_manifest) |> mutate(type = paste0("pb-", model)),
  count_sig_genes(sc_manifest) |> mutate(type = paste0("sc-", model))
) |>
  filter(type %in% names(method_lookup)) |>
  mutate(
    type = factor(method_lookup[type], levels = method_levels),
    cell_type = factor(
      cell_type_lookup[cell_type],
      levels = cell_type_lookup[n_cells_specs$cell_type]
    ),
    sample = sample_label(n_cells_target)
  ) |>
  mutate(
    sample = factor(
      sample,
      levels = c(
        "All cells",
        "Same donors",
        paste0(n_cells_specs$n_cells, " per donor")
      )
    )
  )

ggplot(plot_data, aes(sample, n_sig_gene, fill = type)) +
  geom_col(position = "dodge2") +
  scale_fill_manual(values = method_col_lookup, drop = FALSE) +
  facet_wrap(~cell_type, scales = "free") +
  labs(
    y = "Number of eGenes",
    x = NULL,
    fill = "Method"
  ) +
  theme_jp()

ggsave(
  "power-n-cells-filter-plot.pdf",
  width = 16,
  height = 5.5
)
