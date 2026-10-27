
to_snake <- function(x) {
  out <- tolower(x)
  out <- gsub("[^a-z0-9]", "_", out)
  out <- gsub("_+", "_", out)
  out <- gsub("_$", "", out)
  out
}

theme_jp <- function() {

  font <- "Helvetica"

  theme(
    plot.title = element_text(
      family =  "Helvetica", size = 18, color = "#222222"
    ),
    plot.subtitle = element_text(
      family = font, size = 14, margin = ggplot2::margin(0, 0, 10, 0)
    ),
    plot.caption = element_blank(),
    plot.title.position = "plot",
    legend.position = "top",
    legend.box.margin = margin(t = -5),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(family = font, size = 12, color = "#222222", hjust = 0),
    axis.title = element_text(family = font, size = 12, color = "#222222"),
    axis.text = element_text(family = font, size = 11, color = "#222222"),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid.major.y = element_line(color = "#cbcbcb"),
    panel.grid.minor.y = element_line(color = "#cdcdcd"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(family = font, size = 12)
  )
}

theme_jp_vgrid <- function() {
  theme_jp() %+replace%
    theme(
      panel.grid.major.x = element_line(color = "#cbcbcb"),
      panel.grid.minor.x = element_line(color = "#cdcdcd"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
    )
}

method_lookup <- c(
  "saigeqtl" = "SAIGE-QTL",
  "castie" = "CASTIE",
  "lm" = "Pseudobulk LM",
  "pb-lm" = "Pseudobulk LM",
  "nb_glm" = "Pseudobulk NB-GLM",
  "pb-nb_glm" = "Pseudobulk NB-GLM",
  "p_glmm" = "Pseudobulk P-GLMM",
  "pb-p_glmm" = "Pseudobulk P-GLMM",
  "p_glmm_id" = "Pseudobulk P-GLMM (iid)",
  "pb-p_glmm_id" = "Pseudobulk P-GLMM (iid)",
  "p_glmm_sc" = "Single-cell P-GLMM",
  "sc-p_glmm_sc" = "Single-cell P-GLMM",
  "p_glmm_sc_ablate" = "Single-cell P-GLMM (no donor totals)",
  "sc-p_glmm_sc_ablate" = "Single-cell P-GLMM (no donor totals)",
  "lmm_sc" = "Single-cell LMM",
  "sc-lmm_sc" = "Single-cell LMM",
  "lm_cw" = "Pseudobulk LM with cell weights",
  "pb-lm_cw" = "Pseudobulk LM with cell weights",
  "tensorqtl" = "TensorQTL"
)

method_col_lookup <- c(
  "SAIGE-QTL" = "#EE7733",
  "CASTIE" = "#332288",
  "Pseudobulk LM" = "#AA4499",
  "Pseudobulk NB-GLM" = "#DDCC77",
  "Pseudobulk P-GLMM" = "#117733",
  "Pseudobulk P-GLMM (iid)" = "#999933",
  "Single-cell P-GLMM" = "#882255",
  "Single-cell P-GLMM (no donor totals)" = "#CC6677",
  "Single-cell LMM" = "#44AA99",
  "Pseudobulk LM with cell weights" = "#88CCEE",
  "TensorQTL" = "#0077BB"
)

cell_type_cols <- c(
  "Plasma" = "#66CCEE",
  "B_IN" = "#228833",
  "B_Mem" = "#66A61E",
  "B_all" = "#117733",
  "CD4_NC" = "#EE6677",
  "CD4_ET" = "#CC3311",
  "CD4_SOX4" = "#EE7733",
  "CD8_ET" = "#BB5566",
  "CD8_NC" = "#AA4499",
  "CD8_S100B" = "#882255",
  "NK" = "#332288",
  "NK_R" = "#44AA99",
  "Mono_C" = "#0077BB",
  "Mono_NC" = "#33BBEE",
  "DC" = "#999933",
  "T_all" = "#AA4499",
  "Bergmann" = "#CC6677",
  "Microglia" = "#4477AA",
  "Golgi" = "#EE6677",
  "UBC" = "#228833",
  "Endocytes" = "#CCBB44",
  "B" = "#228833",
  "CD4_T" = "#EE6677",
  "CD8_T" = "#CC6677",
  "monocytes_combined" = "#4477AA",
  "NK_combined" = "#CCBB44"
)

cell_type_lookup <- c(
    "Plasma" = "Plasma",
    "B_IN" = "B IN",
    "B_Mem" = "B Mem",
    "CD4_NC" = "CD4 NC",
    "CD4_ET" = "CD4 ET",
    "CD4_SOX4" = "CD4 SOX4",
    "CD8_ET" = "CD8 ET",
    "CD8_NC" = "CD8 NC",
    "CD8_S100B" = "CD8 S100B",
    "NK" = "NK",
    "NK_R" = "NK R",
    "Mono_C" = "Mono C",
    "Mono_NC" = "Mono NC",
    "DC" = "DC",
    "B_all" = "B All",
    "T_all" = "T All",
    "Bergmann" = "Bergmann",
    "Microglia" = "Microglia",
    "Golgi" = "Golgi",
    "UBC" = "UBC",
    "Endocytes" = "Endocytes",
    "B" = "B",
    "CD4_T" = "CD4 T",
    "CD8_T" = "CD8 T",
    "monocytes_combined" = "Monocytes",
    "NK_combined" = "NK"
)