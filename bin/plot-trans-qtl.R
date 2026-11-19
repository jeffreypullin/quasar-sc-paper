#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(forcats)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)
pb_gwas_files <- read_tsv(args[1], show_col_types = FALSE)
sc_gwas_files <- read_tsv(args[2], show_col_types = FALSE)
annot <- read_tsv(
  args[3],
  show_col_types = FALSE,
  col_types = cols(
    `#chr` = col_character(),
    start = col_double(),
    end = col_double(),
    phenotype_id = col_character()
  )
)
gene_map <- read_tsv(args[4], show_col_types = FALSE)

read_time <- function(path) {
  data <- read.delim(path)
  raw_str <- colnames(data)
  str <- str_extract(raw_str, "(?<=real\\.).*")
  as.numeric(str)
}

norm_chr <- function(x) {
  str_remove(as.character(x), "^chr")
}

cell_types <- c("Plasma", "B_IN", "CD4_NC")
cell_type_levels <- unname(cell_type_lookup[cell_types])
method_levels <- unname(method_lookup[c("lm", "p_glmm_id")])
artefact_genes <- c("ENSG00000256453", "ENSG00000226976", "ENSG00000253626")
# GRCh38 MHC / HLA region (NCBI RefSeq MHC).
hla_chr <- "6"
hla_start <- 28510120
hla_end <- 33480577
trans_mb <- 5e6
sig_threshold <- 5e-8

in_hla <- function(chr, pos) {
  chr == hla_chr & pos >= hla_start & pos <= hla_end
}

tss <- annot |>
  transmute(
    feature_id = phenotype_id,
    gene_chr = norm_chr(`#chr`),
    tss = start
  )

variant_cols <- cols(
  feature_id = col_character(),
  snp_id = col_character(),
  chrom = col_character(),
  pos = col_double(),
  alt = col_character(),
  ref = col_character(),
  maf = col_double(),
  beta = col_double(),
  se = col_double(),
  pvalue = col_double(),
  .default = col_guess()
)

read_sig <- function(path) {
  read_tsv(path, show_col_types = FALSE, col_types = variant_cols)
}

pb_hits <- pb_gwas_files |>
  filter(
    model %in% c("lm", "p_glmm_id"),
    cell_type %in% cell_types
  ) |>
  mutate(
    geno_chr = norm_chr(geno_chr),
    pheno_chr = norm_chr(pheno_chr)
  ) |>
  rowwise() |>
  mutate(var_data = list(read_sig(sig_variant_file))) |>
  ungroup() |>
  unnest(var_data) |>
  mutate(snp_chr = norm_chr(chrom)) |>
  left_join(tss, by = "feature_id") |>
  mutate(
    cis_distance = if_else(
      snp_chr == gene_chr,
      abs(pos - tss),
      Inf
    ),
    is_trans = cis_distance > trans_mb
  ) |>
  filter(is_trans, pvalue < sig_threshold)

discovery_counts <- pb_hits |>
  summarise(
    n_genes = n_distinct(feature_id),
    n_variants = n(),
    .by = c(cell_type, model)
  ) |>
  mutate(
    cell_type = factor(cell_type_lookup[cell_type], levels = cell_type_levels),
    method = factor(method_lookup[model], levels = method_levels)
  )

write_tsv(discovery_counts, "trans-discovery-counts-plot.tsv")

p_discovery <- discovery_counts |>
  ggplot(aes(cell_type, n_genes, fill = method)) +
  geom_col(position = "dodge2") +
  coord_flip() +
  scale_fill_manual(values = method_col_lookup, drop = FALSE) +
  labs(
    x = "Cell type",
    y = "Trans eGenes (P < 5e-8, > 5 Mb)",
    fill = "Method"
  ) +
  theme_jp_vgrid() +
  theme(legend.position = "right")

ggsave("trans-discovery-counts-plot.pdf", p_discovery, width = 10, height = 5)

# Matched runtime: B IN, phenotype chromosome 21.
pb_time <- pb_gwas_files |>
  filter(
    model == "p_glmm_id",
    cell_type == "B_IN",
    norm_chr(pheno_chr) == "21"
  ) |>
  mutate(time = map_dbl(time_file, read_time)) |>
  summarise(time = sum(time), .by = model) |>
  mutate(method = "p_glmm_id")

sc_time <- sc_gwas_files |>
  filter(
    cell_type == "B_IN",
    norm_chr(pheno_chr) == "21"
  ) |>
  mutate(time = map_dbl(time_file, read_time)) |>
  summarise(time = sum(time)) |>
  mutate(method = "p_glmm_sc_ablate")

time_levels <- unname(method_lookup[c("p_glmm_id", "p_glmm_sc_ablate")])

time_plot_data <- bind_rows(pb_time, sc_time) |>
  mutate(
    method = factor(method_lookup[method], levels = time_levels),
    hours = time / 3600
  )

write_tsv(
  time_plot_data |> select(method, time_h = hours),
  "trans-time-p_glmm_id-vs-sc-plot.tsv"
)

p_time <- time_plot_data |>
  ggplot(aes(method, hours, fill = method)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = method_col_lookup, drop = FALSE) +
  labs(
    x = NULL,
    y = "Time (h)",
    fill = "Method"
  ) +
  theme_jp_vgrid() +
  theme(legend.position = "none")

ggsave("trans-time-p_glmm_id-vs-sc-plot.pdf", p_time, width = 8, height = 4)

# Cross-chromosome p_glmm_id hits.
p_glmm_id_hits <- pb_hits |>
  filter(
    model == "p_glmm_id",
    geno_chr != pheno_chr,
    !in_hla(snp_chr, pos),
    !feature_id %in% artefact_genes
  ) |>
  left_join(gene_map, by = "feature_id") |>
  arrange(pvalue) |>
  transmute(
    cell_type,
    feature_id,
    gene_symbol,
    snp_id,
    snp_chr,
    snp_pos = pos,
    gene_chr,
    tss,
    distance_bp = cis_distance,
    beta,
    pvalue
  )

write_tsv(p_glmm_id_hits, "trans-p_glmm_id-hits.tsv")
