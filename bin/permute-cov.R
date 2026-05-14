#!/usr/bin/env Rscript

library(readr)

args <- commandArgs(trailingOnly = TRUE)
ind <- args[[1]]
cell_type <- args[[2]]
int_cov <- args[[3]]
cov <- args[[4]]

set.seed(as.numeric(ind) * 10 + 64)

cov_data <- read_delim(cov, show_col_types = FALSE)
perm <- order(runif(nrow(cov_data)))

perm_col_name <- paste0(int_cov, "_perm")
cov_data[[perm_col_name]] <- cov_data[[int_cov]][perm]

write_delim(
  cov_data,
  paste0("permute-", cell_type, "-", int_cov, "-", ind, ".tsv"),
  delim = "\t",
  col_names = TRUE
)
