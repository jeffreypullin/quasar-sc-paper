#!/usr/bin/env Rscript
# Convert one Randolph Seurat cluster → condition-stratified MTX bundles.
#
# Usage:
#   convert-randolph-seurat.R <seurat.rds> <metadata.csv> <cluster> <out_dir>
#
# Writes NI/ and flu/ subdirectories each containing:
#   matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz, obs.tsv.gz

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(readr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop(paste(
    "Usage: convert-randolph-seurat.R",
    "<seurat.rds> <metadata.csv> <cluster> <out_dir>"
  ))
}
seurat_path <- args[1]
meta_path <- args[2]
cluster <- args[3]
out_dir <- args[4]

CELL_LABEL_MAP <- c(
  "monocytes" = "monocytes combined",
  "infected_monocytes" = "monocytes combined",
  "NK" = "NK combined",
  "NK_high_response" = "NK combined",
  "CD4_T" = "CD4 T",
  "CD8_T" = "CD8 T",
  "B" = "B",
  "DC" = "DC",
  "neutrophils" = "neutrophils",
  "NKT" = "NKT"
)

if (!(cluster %in% names(CELL_LABEL_MAP))) {
  stop("Unknown Randolph cluster: ", cluster)
}

meta_names <- names(read_csv(meta_path, n_max = 0, show_col_types = FALSE))
meta_ref <- read_csv(
  meta_path,
  skip = 1,
  col_names = c("sample_id", meta_names),
  show_col_types = FALSE
)
geno_donors <- unique(meta_ref$indiv_ID[!is.na(meta_ref$YRI_GRCh38)])
# Paper: 89 genotyped donors; HMN52545 lacks genotypes.
geno_donors <- setdiff(geno_donors, "HMN52545")

message("Loading Seurat object: ", seurat_path)
obj <- readRDS(seurat_path)
DefaultAssay(obj) <- "RNA"

obs <- as.data.frame(obj@meta.data)
obs$cell_id <- rownames(obs)
obs$cell_label <- unname(CELL_LABEL_MAP[[cluster]])

# Prefer sample_condition from metadata; fall back to common aliases.
sample_col <- intersect(
  c("sample_condition", "infection_ID", "Sample", "sample"),
  colnames(obs)
)[1]
if (is.na(sample_col)) {
  stop("Could not find sample_condition / infection_ID column in Seurat metadata")
}
obs$sample_condition <- as.character(obs[[sample_col]])

# Donor and condition from sample_condition like HMN83551_NI / HMN83551_flu.
obs$individual <- str_replace(obs$sample_condition, "_(NI|flu)$", "")
obs$condition <- str_extract(obs$sample_condition, "(NI|flu)$")
if (any(is.na(obs$condition))) {
  # Alternate: explicit infection_status column.
  status_col <- intersect(c("infection_status", "condition", "Infection"), colnames(obs))[1]
  if (!is.na(status_col)) {
    raw <- as.character(obs[[status_col]])
    obs$condition <- dplyr::case_when(
      raw %in% c("NI", "mock", "Mock", "control") ~ "NI",
      raw %in% c("flu", "IAV", "iav", "infected") ~ "flu",
      TRUE ~ NA_character_
    )
  }
}
if (any(is.na(obs$condition))) {
  stop("Could not derive NI/flu condition for all cells")
}

# Keep only donors with genotypes.
keep <- obs$individual %in% geno_donors
obs <- obs[keep, , drop = FALSE]

# Sex / age: join from published individual metadata when missing on cells.
indiv_meta <- meta_ref |>
  distinct(indiv_ID, .keep_all = TRUE) |>
  transmute(
    individual = indiv_ID,
    sex = dplyr::case_when(
      gender %in% c("Male", "male", "M", "1") ~ 1L,
      gender %in% c("Female", "female", "F", "2") ~ 2L,
      TRUE ~ NA_integer_
    ),
    age = as.numeric(age),
    age_Scale = as.numeric(age_Scale)
  )
obs <- obs |>
  left_join(indiv_meta, by = "individual")

# Male-only cohort in the paper; default sex=1 if still missing.
obs$sex[is.na(obs$sex)] <- 1L
if (all(is.na(obs$age)) && "age" %in% colnames(obj@meta.data)) {
  obs$age <- as.numeric(obj@meta.data[obs$cell_id, "age"])
}

counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
rm(obj)
gc()
counts <- counts[, obs$cell_id, drop = FALSE]

gene_ids <- rownames(counts)

write_condition <- function(cond) {
  idx <- which(obs$condition == cond)
  if (length(idx) == 0) {
    message(cluster, ": no cells for condition ", cond, "; skipping")
    return(invisible(NULL))
  }
  sub_obs <- obs[idx, , drop = FALSE]
  sub_counts <- counts[, idx, drop = FALSE]

  cond_dir <- file.path(out_dir, cond)
  dir.create(cond_dir, recursive = TRUE, showWarnings = FALSE)

  Matrix::writeMM(sub_counts, file.path(cond_dir, "matrix.mtx"))
  system2("gzip", c("-f", file.path(cond_dir, "matrix.mtx")))

  write_tsv(
    tibble(barcode = colnames(sub_counts)),
    file.path(cond_dir, "barcodes.tsv.gz"),
    col_names = FALSE
  )
  write_tsv(
    tibble(gene_id = gene_ids, GeneSymbol = gene_ids),
    file.path(cond_dir, "features.tsv.gz"),
    col_names = FALSE
  )

  out_obs <- sub_obs |>
    transmute(
      cell_id,
      individual,
      cell_label,
      condition,
      sample_condition,
      sex,
      age,
      age_Scale
    )
  write_tsv(out_obs, file.path(cond_dir, "obs.tsv.gz"))
  message(cond, ": ", ncol(sub_counts), " cells, ", nrow(sub_counts), " genes, ",
          n_distinct(sub_obs$individual), " donors")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_condition("NI")
write_condition("flu")
message("Done. MTX bundles written under ", out_dir)
