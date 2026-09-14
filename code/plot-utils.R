
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
      family =  "Helvetica", size = 24, color = "#222222"
    ),
    plot.subtitle = element_text(
      family = font, size = 20, margin = ggplot2::margin(0, 0, 10, 0)
    ),
    plot.caption = element_blank(),
    plot.title.position = "plot",
    legend.position = "top",
    legend.box.margin = margin(t = -5),
    legend.text.align = 0,
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(family = font, size = 14, color = "#222222"),
    axis.title = element_text(family = font, size = 14, color = "#222222"),
    axis.text = element_text(family = font, size = 14, color = "#222222"),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid.major.y = element_line(color = "#cbcbcb"),
    panel.grid.minor.y = element_line(color = "#cdcdcd"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(family = font, size = 18)
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
  "p_glmm_sc" = "Single-cell P-GLMM",
  "sc-p_glmm_sc" = "Single-cell P-GLMM",
  "lmm_sc" = "Single-cell LMM",
  "sc-lmm_sc" = "Single-cell LMM"
)

method_col_lookup <- c(
  "SAIGE-QTL" = "#EE7733",
  "CASTIE" = "#332288",
  "Pseudobulk LM" = "#AA4499",
  "Pseudobulk NB-GLM" = "#DDCC77",
  "Pseudobulk P-GLMM" = "#117733",
  "Single-cell P-GLMM" = "#882255",
  "Single-cell LMM" = "#44AA99"
)

cell_type_cols <- c(
  "Plasma" = "#66CCEE",
  "B_IN" = "#228833",
  "CD4_NC" = "#EE6677"
)

cell_type_lookup <- c(
    "Plasma" = "Plasma",
    "B_IN" = "B IN",
    "CD4_NC" = "CD4 NC",
    "B_all" = "B all",
    "T_all" = "T all"
)