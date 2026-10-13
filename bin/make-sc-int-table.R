#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
})

args <- commandArgs(trailingOnly = TRUE)

castie_paper <- c("BCL2A1", "HLA-DQA1", "CD96", "UBE2J1", "CYB5A", "FCRL5")

column_lookup <- c(
  "p_glmm_sc counts none" = "P-GLMM (counts)",
  "p_glmm_sc sct_counts none" = "P-GLMM (SCT counts)",
  "p_glmm_sc counts 5" = "P-GLMM K = 5 (counts)",
  "p_glmm_sc sct_counts 5" = "P-GLMM K = 5 (SCT counts)",
  "castie counts none" = "CASTIE (counts)",
  "castie sct_counts none" = "CASTIE (SCT counts)"
)

hits <- read_tsv(args[1], show_col_types = FALSE, col_types = cols(k = "c")) |>
  mutate(column = column_lookup[paste(model, data_type, k)]) |>
  filter(!is.na(column)) |>
  arrange(int_acat_bh)

columns <- c(
  map(set_names(unname(column_lookup)), \(col) filter(hits, column == col)$gene_symbol),
  list("CASTIE (paper)" = castie_paper)
)
columns <- keep(columns, \(genes) length(genes) > 0)
n_rows <- max(lengths(columns))
table <- as_tibble(map(columns, \(genes) c(genes, rep("", n_rows - length(genes)))))

header <- paste0("<th>", names(table), "</th>", collapse = "")
rows <- table |>
  mutate(across(everything(), \(x) paste0("<td>", x, "</td>"))) |>
  unite("row", everything(), sep = "") |>
  pull(row)

html <- c(
  "<!DOCTYPE html>",
  "<html>",
  "<head>",
  "<meta charset=\"utf-8\">",
  "<style>",
  "table { border-collapse: collapse; font-family: Helvetica, Arial, sans-serif; }",
  "th, td { border: 1px solid #d0d0d0; padding: 6px 24px; text-align: left; }",
  "th { background: #f0f0f0; text-align: center; }",
  "td { font-style: italic; }",
  "</style>",
  "</head>",
  "<body>",
  "<table>",
  paste0("<tr>", header, "</tr>"),
  paste0("<tr>", rows, "</tr>"),
  "</table>",
  "</body>",
  "</html>"
)

writeLines(html, args[2])
