#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(forcats)
  library(stringr)
  library(scales)
})

source("/home/jp2045/quasar-sc-paper/code/plot-utils.R")

read_time <- function(path) {
  data <- read.delim(path)
  raw_str <- colnames(data)
  str <- str_extract(raw_str, "(?<=real\\.).*")
  as.numeric(str)
}

args <- commandArgs(trailingOnly = TRUE)

sc_data_files <- read_tsv(args[1], show_col_types = FALSE) |>
  filter(
    k != "seacells",
    cell_frac == 1,
    indiv_frac == 1,
    count_frac == 1,
    min_n_cells < 0,
    is.na(ablate) | ablate == "none",
    int_cov != "none" | cov_spec == "bulk_pca"
  )

if (!"data_type" %in% names(sc_data_files)) {
  sc_data_files$data_type <- "counts"
}
sc_data_files <- mutate(
  sc_data_files,
  data_type = ifelse(is.na(data_type) | data_type == "", "counts", data_type)
)
if ("sc_type" %in% names(sc_data_files)) {
  sc_data_files <- filter(sc_data_files, sc_type == "counts")
}

int_cov_lookup <- c(
  "pseudotime" = "Interaction (pseudotime)",
  "starcat_all" = "Interaction (STARCAT all)",
  "starcat_CD4_Naive" = "Interaction (STARCAT CD4 Naive)"
)

time_data <- sc_data_files |>
  rowwise() |>
  mutate(time = read_time(time_file)) |>
  ungroup() |>
  mutate(
    analysis = case_when(
      as.character(k) != "none" ~ "Grouped",
      int_cov == "none" ~ "Main effect",
      TRUE ~ coalesce(
        unname(int_cov_lookup[as.character(int_cov)]),
        paste0("Interaction (", int_cov, ")")
      )
    )
  ) |>
  summarise(
    time = sum(time),
    .by = c(cell_type, int_cov, analysis, model, data_type)
  )

if (length(args) >= 2) {
  castie_time_data <- read_tsv(args[2], show_col_types = FALSE)
  if (!"data_type" %in% names(castie_time_data)) {
    castie_time_data$data_type <- "counts"
  }
  if (!"int_cov" %in% names(castie_time_data)) {
    castie_time_data$int_cov <- "pseudotime"
  }
  castie_time_data <- castie_time_data |>
    mutate(
      data_type = ifelse(is.na(data_type) | data_type == "", "counts", data_type),
      int_cov = ifelse(is.na(int_cov) | int_cov == "", "pseudotime", int_cov),
      time = coalesce(as.numeric(step1_time), 0) +
        coalesce(as.numeric(step2_time), 0) +
        coalesce(as.numeric(step3_time), 0)
    ) |>
    summarise(time = sum(time), .by = c(cell_type, data_type, int_cov)) |>
    mutate(
      analysis = coalesce(
        unname(int_cov_lookup[as.character(int_cov)]),
        paste0("Interaction (", int_cov, ")")
      ),
      model = "castie"
    )
  time_data <- bind_rows(time_data, castie_time_data)
}

time_data <- time_data |>
  filter(data_type == "counts", cell_type == "B_all") |>
  mutate(
    analysis = fct_reorder(analysis, time, .fun = max),
    model = factor(
      method_lookup[model],
      levels = method_lookup[c("p_glmm_sc", "lmm_sc", "castie")]
    ),
    hours = time / (60 * 60)
  )

p <- time_data |>
  filter(model != method_lookup[["castie"]]) |>
  mutate(model = droplevels(model)) |>
  ggplot(aes(model, hours, fill = analysis)) +
  geom_col(position = "dodge2") +
  coord_flip() +
  labs(
    x = NULL,
    y = "Time (h)",
    fill = NULL
  ) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_jp_vgrid() +
  theme(legend.position = "right")

if (length(args) < 2) {
  ggsave(
    "time-int-plot.pdf",
    p,
    width = 7,
    height = 3
  )
} else {
  panel_levels <- c("Interaction", "Grouped")
  model_levels <- c(
    method_lookup[["p_glmm_sc"]],
    "Grouped",
    method_lookup[["castie"]]
  )
  interaction_pglmm <- time_data |>
    filter(
      analysis == "Interaction (pseudotime)",
      model == method_lookup[["p_glmm_sc"]]
    ) |>
    mutate(panel = "Grouped", model = method_lookup[["p_glmm_sc"]])
  panel_data <- time_data |>
    filter(
      analysis %in% c("Interaction (pseudotime)", "Grouped"),
      model %in% method_lookup[c("p_glmm_sc", "castie")]
    ) |>
    mutate(
      panel = if_else(analysis == "Grouped", "Grouped", "Interaction"),
      model = if_else(analysis == "Grouped", "Grouped", as.character(model))
    ) |>
    bind_rows(interaction_pglmm) |>
    mutate(
      panel = factor(panel, levels = panel_levels),
      model = factor(model, levels = model_levels)
    ) |>
    mutate(
      fold = hours / min(hours),
      label = if_else(fold > 1, paste0(round(fold), "x"), ""),
      .by = panel
    )

  p_panels <- panel_data |>
    ggplot(aes(model, hours, fill = model)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = label), hjust = -0.2) +
    scale_x_discrete(labels = \(x) sub("Single-cell ", "Single-cell\n", x)) +
    scale_fill_manual(
      values = c(method_col_lookup, "Grouped" = "#377EB8"),
      drop = FALSE
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    facet_wrap(~panel, ncol = 1, scales = "free") +
    coord_flip() +
    labs(x = NULL, y = "Time (h)") +
    theme_jp_vgrid() +
    theme(
      strip.text = element_blank(),
      strip.background = element_blank()
    )

  ggsave(
    "time-context-plot.pdf",
    p_panels,
    width = 7,
    height = 4.5
  )
  saveRDS(p_panels, "time-context-plot.rds")
}
