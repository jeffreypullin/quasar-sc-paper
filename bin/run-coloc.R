#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(coloc)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)
pb_quasar_file <- args[[1]]
sc_quasar_file <- args[[2]]
annot_bed <- args[[3]]
finngen_dir <- args[[4]]

gws_p <- 5e-8
window_bp <- 1e6L
min_shared_snps <- 10L
pp4_thresh <- 0.8
n_eqtl <- 981L
eur_bim_path <- "/home/jp2045/quasar-sc-paper/data/reference/EUR.bim"
chain_path <- "/home/jp2045/quasar-sc-paper/data/reference/hg38ToHg19.over.chain.gz"
chain_url <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"

if (!file.exists(chain_path)) {
  download.file(chain_url, chain_path, quiet = TRUE)
}

# EUR.bim (hg19): map chrom:pos:ref:alt -> rsid (both allele orientations)
bim <- fread(eur_bim_path, header = FALSE,
             col.names = c("chr", "rsid", "cm", "pos", "a1", "a2"))
bim <- bim[grepl("^rs", rsid)]
rsid_map <- unique(rbind(
  bim[, .(key = paste(chr, pos, a2, a1, sep = ":"), rsid)],
  bim[, .(key = paste(chr, pos, a1, a2, sep = ":"), rsid)]
), by = "key")

annot <- read_tsv(annot_bed, show_col_types = FALSE) |>
  transmute(feature_id = phenotype_id, gene_chr = as.character(`#chr`), tss = as.integer(start))

fg <- read_tsv(file.path(finngen_dir, "selected_manifest.tsv"), show_col_types = FALSE) |>
  mutate(
    path = file.path(finngen_dir, paste0("finngen_R13_", phenocode, ".gz")),
    trait_s = num_cases / (num_cases + num_controls)
  )

pb <- read_tsv(pb_quasar_file, show_col_types = FALSE) |>
  filter(dataset == "onek1k", model == "lm", pb_type == "logcounts",
         int_cov == "none", cell_frac == 1, indiv_frac == 1,
         min_n_cells < 0, count_frac == 1)

sc <- read_tsv(sc_quasar_file, show_col_types = FALSE) |>
  filter(dataset == "onek1k", model == "p_glmm_sc", data_type == "counts",
         cov_spec == "bulk_pca", k == "none", int_cov == "none",
         cell_frac == 1, indiv_frac == 1, min_n_cells < 0, count_frac == 1,
         ablate %in% c("none", ""))

read_egenes <- function(manifest, method) {
  manifest |>
    distinct(cell_type, region_file) |>
    mutate(genes = map(region_file, \(p) read_tsv(p, show_col_types = FALSE) |>
                         select(feature_id, pvalue))) |>
    select(-region_file) |>
    unnest(genes) |>
    mutate(region_bh = p.adjust(pvalue, method = "BH"), .by = cell_type) |>
    filter(region_bh < 0.05) |>
    transmute(cell_type, method, feature_id)
}

egenes <- bind_rows(read_egenes(pb, "lm"), read_egenes(sc, "p_glmm_sc")) |>
  left_join(annot, by = "feature_id")

var_index <- bind_rows(
  pb |> transmute(method = "lm", cell_type, chr = str_remove(chr, "^chr"), variant_file),
  sc |> transmute(method = "p_glmm_sc", cell_type, chr = str_remove(chr, "^chr"), variant_file)
)

# FinnGen GWS hits -> merged ±1 Mb windows (hg38)
windows <- map_dfr(seq_len(nrow(fg)), function(i) {
  path <- fg$path[[i]]
  pheno <- fg$phenocode[[i]]
  hits <- fread(cmd = sprintf(
    "gzip -cd %s | awk -F'\\t' 'NR>1 && ($7+0)<%s {print $1,$2,$5,$7}'",
    shQuote(path), format(gws_p, scientific = TRUE)
  ), header = FALSE, col.names = c("chr", "pos", "rsid", "pval"))
  hits[, `:=`(
    phenocode = pheno,
    chr = as.character(chr),
    start = pmax(as.integer(pos) - window_bp, 1L),
    end = as.integer(pos) + window_bp
  )]
  setorder(hits, chr, start, end)
  hits[, grp := cumsum(c(TRUE, start[-1] > cummax(end)[-.N])), by = chr]
  hits[, .(
    start = min(start), end = max(end),
    lead_pos = pos[which.min(pval)],
    lead_rsid = rsid[which.min(pval)],
    lead_pval = min(pval)
  ), by = .(phenocode, chr, grp)][, grp := NULL] |>
    as_tibble()
}) |>
  left_join(fg |> select(phenocode, phenotype, path, trait_s), by = "phenocode")

# Lift windows hg38 -> hg19 for TSS overlap
bed_in <- tempfile(fileext = ".bed")
bed_out <- tempfile(fileext = ".bed")
bed_unmap <- tempfile(fileext = ".unmapped")
windows |>
  mutate(
    name = as.character(row_number()),
    chrom = paste0("chr", chr),
    start0 = start - 1L
  ) |>
  select(chrom, start0, end, name) |>
  write_tsv(bed_in, col_names = FALSE)
system2("liftOver", c(bed_in, chain_path, bed_out, bed_unmap))
lifted <- read_tsv(bed_out, col_names = c("chrom", "start0", "end1", "name"),
                   col_types = "ciic", show_col_types = FALSE) |>
  transmute(
    name,
    chr_hg19 = str_remove(chrom, "^chr"),
    start_hg19 = start0 + 1L,
    end_hg19 = end1
  )
windows <- windows |>
  mutate(name = as.character(row_number())) |>
  inner_join(lifted, by = "name") |>
  select(-name)

candidates <- egenes |>
  inner_join(windows, by = c("gene_chr" = "chr_hg19"), relationship = "many-to-many") |>
  filter(tss >= start_hg19, tss <= end_hg19)

# Cache variant files and GWAS windows
eqtl_cache <- new.env(parent = emptyenv())
gwas_cache <- new.env(parent = emptyenv())

get_eqtl <- function(method, cell_type, chrom, feature_id) {
  key <- paste(method, cell_type, chrom, sep = "::")
  if (!exists(key, envir = eqtl_cache, inherits = FALSE)) {
    path <- var_index |>
      filter(method == .env$method, cell_type == .env$cell_type, chr == .env$chrom) |>
      pull(variant_file)
    tbl <- fread(path[[1]], showProgress = FALSE) |>
      as_tibble() |>
      mutate(
        chrom = as.character(chrom),
        key = paste(chrom, pos, ref, alt, sep = ":")
      ) |>
      inner_join(as_tibble(rsid_map), by = "key")
    assign(key, tbl, envir = eqtl_cache)
  }
  get(key, envir = eqtl_cache) |> filter(feature_id == .env$feature_id)
}

get_gwas <- function(path, chrom, start, end) {
  key <- paste(path, chrom, start, end, sep = "::")
  if (!exists(key, envir = gwas_cache, inherits = FALSE)) {
    tbl <- fread(
      cmd = sprintf("tabix %s %s:%d-%d", shQuote(path), chrom, start, end),
      header = FALSE,
      col.names = c("chrom", "pos", "ref", "alt", "rsids", "nearest_genes",
                    "pval", "mlogp", "beta", "sebeta", "af_alt",
                    "af_alt_cases", "af_alt_controls")
    ) |>
      as_tibble() |>
      transmute(
        rsid = str_replace(rsids, ";.*$", ""),
        beta, varbeta = sebeta^2, ref, alt, af_alt
      ) |>
      filter(str_starts(rsid, "rs")) |>
      distinct(rsid, .keep_all = TRUE)
    assign(key, tbl, envir = gwas_cache)
  }
  get(key, envir = gwas_cache)
}

run_coloc <- function(eqtl, gwas, trait_s) {
  shared <- eqtl |>
    transmute(rsid, beta_eqtl = beta, varbeta_eqtl = se^2, maf,
              ref_eqtl = ref, alt_eqtl = alt) |>
    inner_join(
      gwas |> transmute(rsid, beta_gwas = beta, varbeta_gwas = varbeta, af_alt,
                        ref_gwas = ref, alt_gwas = alt),
      by = "rsid"
    ) |>
    mutate(
      flip = ref_eqtl == alt_gwas & alt_eqtl == ref_gwas,
      keep = flip | (ref_eqtl == ref_gwas & alt_eqtl == alt_gwas),
      beta_gwas = if_else(flip, -beta_gwas, beta_gwas)
    ) |>
    filter(keep) |>
    distinct(rsid, .keep_all = TRUE)

  if (nrow(shared) < min_shared_snps) {
    return(NULL)
  }

  fit <- suppressMessages(coloc.abf(
    list(snp = shared$rsid, beta = shared$beta_eqtl, varbeta = shared$varbeta_eqtl,
         MAF = pmin(shared$maf, 1 - shared$maf), N = n_eqtl, type = "quant"),
    list(snp = shared$rsid, beta = shared$beta_gwas, varbeta = shared$varbeta_gwas,
         MAF = pmin(shared$af_alt, 1 - shared$af_alt), type = "cc", s = trait_s)
  ))
  as_tibble(as.list(fit$summary))
}

results <- map(seq_len(nrow(candidates)), function(i) {
  row <- candidates[i, ]
  gwas <- get_gwas(row$path, row$chr, row$start, row$end)
  eqtl <- get_eqtl(row$method, row$cell_type, row$gene_chr, row$feature_id)
  if (nrow(gwas) == 0L || nrow(eqtl) == 0L) {
    return(NULL)
  }
  fit <- run_coloc(eqtl, gwas, row$trait_s)
  if (is.null(fit)) {
    return(NULL)
  }
  bind_cols(
    tibble(
      phenocode = row$phenocode, phenotype = row$phenotype,
      cell_type = row$cell_type, feature_id = row$feature_id, method = row$method,
      window_chr = row$chr, window_start = row$start, window_end = row$end,
      lead_gwas_pos_hg38 = row$lead_pos, lead_gwas_rsid = row$lead_rsid,
      lead_gwas_pval = row$lead_pval
    ),
    fit
  )
}) |>
  bind_rows() |>
  slice_max(PP.H4.abf, n = 1, with_ties = FALSE,
            by = c(phenocode, cell_type, feature_id, method))

write_tsv(results, "coloc-results.tsv")

wide <- results |>
  mutate(
    coloc = PP.H4.abf >= pp4_thresh,
    method = factor(method, levels = c("lm", "p_glmm_sc"))
  ) |>
  select(phenocode, cell_type, feature_id, method, coloc) |>
  pivot_wider(names_from = method, values_from = coloc,
              values_fill = FALSE, names_expand = TRUE) |>
  mutate(category = case_when(
    lm & p_glmm_sc ~ "both",
    p_glmm_sc & !lm ~ "p_glmm_sc only",
    lm & !p_glmm_sc ~ "lm only",
    TRUE ~ "neither"
  ))

cats <- c("both", "p_glmm_sc only", "lm only")
sig <- wide |> filter(category != "neither")

summary_tbl <- bind_rows(
  sig |> count(category, name = "n_loci") |>
    complete(category = cats, fill = list(n_loci = 0L)) |>
    mutate(level = "overall", phenocode = NA_character(), cell_type = NA_character()),
  sig |> count(phenocode, category, name = "n_loci") |>
    complete(phenocode, category = cats, fill = list(n_loci = 0L)) |>
    mutate(level = "phenocode", cell_type = NA_character()),
  sig |> count(cell_type, category, name = "n_loci") |>
    complete(cell_type, category = cats, fill = list(n_loci = 0L)) |>
    mutate(level = "cell_type", phenocode = NA_character())
) |>
  select(level, phenocode, cell_type, category, n_loci)

write_tsv(summary_tbl, "coloc-summary.tsv")

overall <- summary_tbl |> filter(level == "overall") |>
  mutate(category = factor(category, levels = cats))

ggplot(overall, aes(category, n_loci, fill = category)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_loci), vjust = -0.3, size = 5) +
  scale_fill_manual(values = c(
    "both" = "#4477AA",
    "p_glmm_sc only" = unname(method_col_lookup[["Single-cell P-GLMM"]]),
    "lm only" = unname(method_col_lookup[["Pseudobulk LM"]])
  )) +
  theme_jp() +
  theme(legend.position = "none") +
  labs(title = "FinnGen × OneK1K eQTL coloc",
       subtitle = sprintf("PP4 >= %.1f", pp4_thresh),
       x = NULL, y = "Number of loci")
ggsave("coloc-summary.pdf", width = 9, height = 6)
