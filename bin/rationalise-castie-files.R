#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

input_tsv <- read_tsv(args[[1]], show_col_types = FALSE)

cell_type <- input_tsv$cell_type[[1]]
chrom <- input_tsv$chrom[[1]]
data_type <- if ("data_type" %in% names(input_tsv)) {
  as.character(input_tsv$data_type[[1]])
} else {
  NA_character_
}
int_cov <- if ("int_cov" %in% names(input_tsv)) {
  as.character(input_tsv$int_cov[[1]])
} else {
  "pseudotime"
}
gene_list <- str_replace(
  args[[1]],
  "^(.+?)-\\1-(genes-\\d+\\.txt).*$",
  "\\1-\\2"
)

read_castie_tsv <- function(path) {
  read_tsv(
    path,
    show_col_types = FALSE,
    col_types = cols(.default = col_character())
  )
}

all_variant_data <- tibble()
error_genes <- c()
for (i in seq_len(nrow(input_tsv))) {
  gene <- input_tsv$gene[[i]]

  if (is.na(input_tsv$variant_file[[i]]) || input_tsv$variant_file[[i]] == "NA") {
    error_genes <- c(error_genes, gene)
    next
  }

  variant_data <- read_castie_tsv(input_tsv$variant_file[[i]])
  if (nrow(variant_data) == 0L) {
    error_genes <- c(error_genes, gene)
    next
  }

  all_variant_data <- bind_rows(
    all_variant_data,
    mutate(variant_data, gene = gene)
  )
}

region_paths <- unique(input_tsv$region_file)
region_paths <- region_paths[!is.na(region_paths) & region_paths != "NA"]
all_region_data <- tibble()
for (path in region_paths) {
  region_data <- read_castie_tsv(path)
  if (nrow(region_data) == 0L) {
    next
  }
  all_region_data <- bind_rows(all_region_data, region_data)
}

if (nrow(all_variant_data) > 0L) {
  all_variant_data <- suppressMessages(type_convert(all_variant_data))
}
if (nrow(all_region_data) > 0L) {
  all_region_data <- suppressMessages(type_convert(all_region_data))
}

step1_time <- sum(input_tsv$step1_time)
step2_time <- sum(input_tsv$step2_time)
step3_time <- sum(input_tsv$step3_time)

prefix <- str_remove(args[[1]], "-castie-files\\.tsv$")

all_variant_file <- paste0("rationalised-", prefix, "-variant.tsv")
all_region_file <- paste0("rationalised-", prefix, "-region.tsv")

write_tsv(all_variant_data, all_variant_file)
write_tsv(all_region_data, all_region_file)

result_tibble <- tibble(
  cell_type = cell_type,
  model = "castie",
  data_type = data_type,
  int_cov = int_cov,
  k = "none",
  chrom = chrom,
  gene_list = gene_list,
  step1_time = step1_time,
  step2_time = step2_time,
  step3_time = step3_time,
  variant_file = paste0(getwd(), "/", all_variant_file),
  region_file = paste0(getwd(), "/", all_region_file),
  error_genes = list(error_genes)
)

write_tsv(result_tibble, paste0("rationalised-", args[[1]]))

all_files <- list.files(".")
files_to_delete <- all_files[!str_starts(all_files, "rationalised-")]
file.remove(files_to_delete)
