#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

to_snake <- function(x) {
  out <- tolower(x)
  out <- gsub("[^a-z0-9]", "_", out)
  out <- gsub("_+", "_", out)
  out <- gsub("_$", "", out)
  out
}

output_cols <- c(
  "cell_type",
  "model",
  "data_type",
  "int_cov",
  "feature_id",
  "snp_id",
  "snp_beta",
  "snp_se",
  "snp_pvalue",
  "int_beta",
  "int_se",
  "int_pvalue",
  "int_acat_pvalue",
  "int_acat_bh"
)

variant_cols <- setdiff(output_cols, c("int_acat_pvalue", "int_acat_bh"))

empty_hits <- tibble::tibble(
  cell_type = character(),
  model = character(),
  data_type = character(),
  int_cov = character(),
  feature_id = character(),
  snp_id = character(),
  snp_beta = double(),
  snp_se = double(),
  snp_pvalue = double(),
  int_beta = double(),
  int_se = double(),
  int_pvalue = double(),
  int_acat_pvalue = double(),
  int_acat_bh = double()
)

read_file_list <- function(path) {
  data <- read_tsv(path, show_col_types = FALSE)
  if (!"k" %in% names(data)) {
    data$k <- "none"
  }
  if (!"data_type" %in% names(data)) {
    data$data_type <- NA_character_
  }
  data |>
    mutate(
      k = ifelse(is.na(k) | k == "", "none", k),
      data_type = ifelse(is.na(data_type) | data_type == "", NA_character_, data_type),
      data_type = ifelse(is.na(data_type) & model == "castie", "counts", data_type)
    ) |>
    filter(k == "none", int_cov != "none")
}

# Quasar region files use int_<cov>_acat_pvalue (or legacy int_acat_pvalue).
# CASTIE step 3 is long format (Gene, pval_column, ACAT_p).
normalise_int_region <- function(dt, int_cov) {
  if (nrow(dt) == 0L) {
    return(NULL)
  }
  if (all(c("Gene", "pval_column", "ACAT_p") %in% names(dt))) {
    keep <- dt[as.character(pval_column) == as.character(int_cov)]
    if (nrow(keep) == 0L) {
      return(NULL)
    }
    out <- keep[, .(
      feature_id = as.character(Gene),
      int_cov = as.character(int_cov),
      int_acat_pvalue = as.numeric(ACAT_p)
    )]
  } else if ("feature_id" %in% names(dt)) {
    acat_cols <- grep("^int_.*_acat_pvalue$", names(dt), value = TRUE)
    if (length(acat_cols) == 0L && "int_acat_pvalue" %in% names(dt)) {
      acat_cols <- "int_acat_pvalue"
    }
    if (length(acat_cols) == 0L) {
      return(NULL)
    }
    out <- rbindlist(lapply(acat_cols, function(col) {
      term <- if (col == "int_acat_pvalue") {
        to_snake(int_cov)
      } else {
        sub("_acat_pvalue$", "", sub("^int_", "", col))
      }
      dt[, .(
        feature_id = as.character(feature_id),
        int_cov = term,
        int_acat_pvalue = as.numeric(get(col))
      )]
    }))
  } else {
    return(NULL)
  }
  out <- out[!is.na(feature_id) & feature_id != "" & !is.na(int_acat_pvalue)]
  if (nrow(out) == 0L) {
    return(NULL)
  }
  out
}

normalise_int_variants <- function(dt) {
  if ("MarkerID" %in% names(dt)) {
    gene_col <- if ("feature_id" %in% names(dt)) {
      "feature_id"
    } else if ("gene" %in% names(dt)) {
      "gene"
    } else if ("Gene" %in% names(dt)) {
      "Gene"
    } else {
      return(dt)
    }
    ge_p <- as.character(dt$pval_ge)
    ge_b <- as.character(dt$Beta_ge)
    ge_s <- as.character(dt$seBeta_ge)
    dt[, `:=`(
      feature_id = as.character(get(gene_col)),
      snp_id = MarkerID,
      snp_beta = BETA,
      snp_se = SE,
      snp_pvalue = `p.value`,
      snp_x_pseudotime_pvalue = as.numeric(sub(",.*", "", ge_p)),
      snp_x_pseudotime_beta = as.numeric(sub(",.*", "", ge_b)),
      snp_x_pseudotime_se = as.numeric(sub(",.*", "", ge_s))
    )]
  }
  dt
}

sc_data_files <- read_file_list(args[[1]])
if (length(args) >= 2) {
  sc_data_files <- bind_rows(sc_data_files, read_file_list(args[[2]]))
}

if (nrow(sc_data_files) == 0) {
  write_tsv(empty_hits, "sc-int-eqtl-hits.tsv")
  quit(save = "no", status = 0)
}

region_list <- list()
for (i in seq_len(nrow(sc_data_files))) {
  info <- sc_data_files[i, ]
  region_data <- fread(info$region_file, showProgress = FALSE)
  region_data <- normalise_int_region(region_data, info$int_cov)
  if (is.null(region_data)) {
    next
  }
  region_data[, `:=`(
    cell_type = info$cell_type,
    model = info$model,
    data_type = info$data_type
  )]
  region_list[[i]] <- region_data[, .(
    cell_type, model, data_type, int_cov, feature_id, int_acat_pvalue
  )]
}

region_list <- region_list[!vapply(region_list, is.null, logical(1))]

if (length(region_list) == 0) {
  write_tsv(empty_hits, "sc-int-eqtl-hits.tsv")
  quit(save = "no", status = 0)
}

all_int_p <- bind_rows(region_list) |>
  group_by(cell_type, model, data_type, int_cov, feature_id) |>
  slice_min(int_acat_pvalue, n = 1, with_ties = FALSE) |>
  ungroup()

# Same BH n for every model in a cell type / context / data type: the largest
# number of genes any method tested. CASTIE chunks of 50 would otherwise get a
# much easier FDR if a run is incomplete or tests a subset.
n_tests <- all_int_p |>
  summarise(n_model = n(), .by = c(cell_type, model, data_type, int_cov)) |>
  summarise(n_tests = max(n_model), .by = c(cell_type, int_cov, data_type))

int_egenes <- all_int_p |>
  left_join(n_tests, by = c("cell_type", "int_cov", "data_type")) |>
  mutate(
    int_acat_bh = p.adjust(int_acat_pvalue, method = "BH", n = n_tests[1]),
    .by = c(cell_type, model, data_type, int_cov)
  ) |>
  filter(int_acat_bh < 0.05) |>
  select(cell_type, model, data_type, int_cov, feature_id, int_acat_pvalue, int_acat_bh)

if (nrow(int_egenes) == 0) {
  write_tsv(empty_hits, "sc-int-eqtl-hits.tsv")
  quit(save = "no", status = 0)
}

gene_int_leads <- function(dt, pval_col) {
  keep <- dt[!is.na(get(pval_col)) & !is.na(feature_id) & !is.na(snp_id)]
  if (nrow(keep) == 0L) {
    return(keep)
  }
  keep[order(get(pval_col)), .SD[1L], by = feature_id]
}

lead_list <- list()

for (i in seq_len(nrow(sc_data_files))) {
  info <- sc_data_files[i, ]
  egene_rows <- int_egenes |>
    filter(
      cell_type == info$cell_type,
      model == info$model,
      (is.na(data_type) & is.na(info$data_type)) |
        (!is.na(data_type) & data_type == info$data_type)
    )
  if (nrow(egene_rows) == 0L) {
    next
  }

  variant_data <- fread(info$variant_file, showProgress = FALSE)
  variant_data <- normalise_int_variants(variant_data)

  for (term in unique(egene_rows$int_cov)) {
    int_p_col <- paste0("snp_x_", term, "_pvalue")
    int_beta_col <- paste0("snp_x_", term, "_beta")
    int_se_col <- paste0("snp_x_", term, "_se")
    if (!all(c(int_p_col, int_beta_col, int_se_col) %in% names(variant_data))) {
      next
    }
    egene_features <- egene_rows |>
      filter(int_cov == term) |>
      pull(feature_id)
    term_data <- variant_data[feature_id %in% egene_features]
    leads <- gene_int_leads(term_data, int_p_col)
    if (nrow(leads) == 0L) {
      next
    }
    leads[, `:=`(
      cell_type = info$cell_type,
      model = info$model,
      data_type = info$data_type,
      int_cov = term,
      int_beta = get(int_beta_col),
      int_se = get(int_se_col),
      int_pvalue = get(int_p_col)
    )]
    lead_list[[length(lead_list) + 1L]] <- leads[, ..variant_cols]
  }
  rm(variant_data)
}

lead_list <- lead_list[!vapply(lead_list, is.null, logical(1))]

if (length(lead_list) == 0) {
  write_tsv(empty_hits, "sc-int-eqtl-hits.tsv")
  quit(save = "no", status = 0)
}

int_hits <- bind_rows(lead_list) |>
  group_by(cell_type, model, data_type, int_cov, feature_id) |>
  arrange(int_pvalue, .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup() |>
  inner_join(
    int_egenes,
    by = c("cell_type", "model", "data_type", "int_cov", "feature_id")
  ) |>
  arrange(cell_type, model, data_type, int_acat_bh, int_pvalue) |>
  select(all_of(output_cols))

write_tsv(int_hits, "sc-int-eqtl-hits.tsv")
