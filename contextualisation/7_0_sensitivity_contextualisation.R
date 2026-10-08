# Run from the repository root in a fresh R session. See README.md.
# EDIT LOCAL PATHS HERE. No study data or results are supplied.
PROJECT_ROOT <- "."
DATA_DIR <- file.path(PROJECT_ROOT, "data")
OUTPUT_ROOT <- file.path(PROJECT_ROOT, "outputs")
MANUSCRIPT_DIR <- file.path(OUTPUT_ROOT, "manuscript")
REFERENCE_DIR <- file.path(DATA_DIR, "reference")
SIGNATURE_SCRIPTS <- file.path(PROJECT_ROOT, "signatures")
CONTEXT_SCRIPTS <- file.path(PROJECT_ROOT, "contextualisation")
dir.create(MANUSCRIPT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Protein-level contextualisation figures
# One row per protein
#
# Cancer grouping:
#   1. Overall and other site-specific cancers:
#      All cancer, breast cancer, prostate cancer, endometrial cancer
#   2. Colorectal cancer
#
# Output:
#   - Preview as low-resolution PNG
#   - Save final as 300 DPI PNG
# ============================================================

library(tidyverse)
library(readxl)
library(openxlsx)
library(patchwork)
library(ggplot2)
library(cowplot)

# ------------------------------------------------------------
# 0. Paths and settings
# ------------------------------------------------------------

script_dir <- CONTEXT_SCRIPTS
out_dir    <- MANUSCRIPT_DIR

context_out_dir <- file.path(out_dir, "protein_contextualisation_by_disease")
dir.create(context_out_dir, showWarnings = FALSE, recursive = TRUE)

sex_script     <- file.path(script_dir, "6_8_1_plot_sexinteraction.R")
smoking_script <- file.path(script_dir, "6_8_2_plotsmokinginteractions.R")
alcohol_script <- file.path(script_dir, "6_8_3_plot_alcoholinteractions.R")
age_script     <- file.path(script_dir, "6_8_5_plot_aging.R")

mr_file <- file.path(out_dir, "MR_results_with sequence_ids_manual_checks.xlsx")

# Use 0.50 for moderate + strong colocalization; use 0.80 for strong only
min_pph4 <- 0.50

# ------------------------------------------------------------
# 1. Text sizes and stronger colour palette
# ------------------------------------------------------------

TXT_SMALL <- 8
TXT_MAIN  <- 10
TXT_TITLE <- 12

pt_to_geom <- function(x) x / ggplot2::.pt

# COL_MALE    <- "#253858"
# COL_FEMALE  <- "#C15A3A"
# 
# COL_AGE50   <- "#C15A3A"
# COL_AGE55   <- "#6B778C"
# COL_AGE60   <- "#253858"
# 
# COL_NEVERSMOKE <- "#253858"
# COL_EVERSMOKE  <- "#C15A3A"
# 
# COL_NEVERALC   <- "#253858"
# COL_OCCALC     <- "#6B778C"
# COL_MODALC     <- "#C15A3A"
# COL_REGALC     <- "#8B0000"

# ------------------------------------------------------------
# Stronger colour palette
# ------------------------------------------------------------

COL_MALE    <- "#003B73"   # deep blue
COL_FEMALE  <- "#D1495B"   # strong rose/red

COL_AGE50   <- "#E66100"   # vivid orange
COL_AGE55   <- "#7A5195"   # purple
COL_AGE60   <- "#003B73"   # deep blue

COL_NEVERSMOKE <- "#003B73"
COL_EVERSMOKE  <- "#D1495B"

COL_NEVERALC   <- "#003B73"
COL_OCCALC     <- "#7A5195"
COL_MODALC     <- "#E66100"
COL_REGALC     <- "#A50026"   # dark red

# ------------------------------------------------------------
# 2. Source existing interaction scripts safely
# ------------------------------------------------------------

source_get <- function(script, object_name) {
  env <- new.env()
  source(script, local = env)
  if (!exists(object_name, envir = env, inherits = FALSE)) {
    stop(sprintf("%s did not create %s. Check the script selected; see README review items.",
                 basename(script), object_name))
  }
  get(object_name, envir = env, inherits = FALSE)
}

plotdat_sex     <- source_get(sex_script, "plotdat_sex")
plotdat_smoking <- source_get(smoking_script, "plotdat_smoking")
plotdat_alcohol <- source_get(alcohol_script, "plotdat_alcohol")
plotdat_age     <- source_get(age_script, "plotdat_age")

# ------------------------------------------------------------
# 3. Read MR + coloc results
# ------------------------------------------------------------

evidence_all <- read_excel(mr_file, sheet = excel_sheets(mr_file)[1]) %>%
  mutate(
    gene = protein,
    
    outcome_clean = case_when(
      disease == "Cancer"      ~ "All cancer",
      disease == "Breast"      ~ "Breast cancer",
      disease == "Prostate"    ~ "Prostate cancer",
      disease == "Colorectal"  ~ "Colorectal cancer",
      disease == "Endometrial" ~ "Endometrial cancer",
      disease == "Lung"        ~ "Lung cancer",
      disease == "Liver"       ~ "Liver cancer",
      disease == "Kidney"      ~ "Kidney cancer",
      disease == "Gastric"     ~ "Gastric cancer",
      disease == "Myeloma"     ~ "Myeloma",
      disease == "T2D"         ~ "T2D",
      disease == "CVD"         ~ "CVD",
      TRUE ~ disease
    ),
    
    disease_group = case_when(
      outcome_clean == "T2D" ~ "T2D",
      outcome_clean == "CVD" ~ "CVD",
      TRUE ~ "Cancer"
    ),
    
    cancer_page = case_when(
      outcome_clean %in% c(
        "All cancer",
        "Breast cancer",
        "Prostate cancer",
        "Endometrial cancer"
      ) ~ "Overall and other site-specific cancers",
      
      outcome_clean == "Colorectal cancer" ~ "Colorectal cancer",
      
      TRUE ~ NA_character_
    ),
    
    OR = or,
    lower_CI = or_lci95,
    upper_CI = or_uci95,
    PP_H4 = PP_H4_final,
    
    coloc_simple = case_when(
      PP_H4 >= 0.80 ~ "Strong",
      PP_H4 >= 0.50 ~ "Moderate",
      TRUE ~ "Limited"
    ),
    
    outcome_clean = factor(
      outcome_clean,
      levels = c(
        "All cancer",
        "Breast cancer",
        "Prostate cancer",
        "Colorectal cancer",
        "Endometrial cancer",
        "Lung cancer",
        "Liver cancer",
        "Kidney cancer",
        "Gastric cancer",
        "Myeloma",
        "T2D",
        "CVD"
      )
    ),
    
    cancer_page = factor(
      cancer_page,
      levels = c(
        "Overall and other site-specific cancers",
        "Colorectal cancer"
      )
    )
  ) %>%
  filter(
    !is.na(pdi_type),
    !is.na(protein_id),
    !is.na(OR),
    !is.na(lower_CI),
    !is.na(upper_CI),
    !is.na(PP_H4)
  )

# ------------------------------------------------------------
# 4. Keep only proteins available in interaction plot data
# ------------------------------------------------------------

available_proteins <- bind_rows(
  plotdat_sex     %>% distinct(pdi_type, gene, protein_id),
  plotdat_smoking %>% distinct(pdi_type, gene, protein_id),
  plotdat_alcohol %>% distinct(pdi_type, gene, protein_id),
  plotdat_age     %>% distinct(pdi_type, gene, protein_id)
) %>%
  distinct()

evidence_plot <- evidence_all %>%
  semi_join(
    available_proteins %>% distinct(pdi_type, protein_id),
    by = c("pdi_type", "protein_id")
  )

# ------------------------------------------------------------
# 5. Contextualisation evidence
# ------------------------------------------------------------

evidence_context <- evidence_plot %>%
  filter(
    PP_H4 >= min_pph4,
    outcome_clean %in% c(
      "All cancer",
      "Breast cancer",
      "Prostate cancer",
      "Colorectal cancer",
      "Endometrial cancer",
      "T2D",
      "CVD"
    )
  ) %>%
  arrange(
    pdi_type,
    disease_group,
    cancer_page,
    outcome_clean,
    gene,
    protein_id,
    desc(PP_H4)
  ) %>%
  distinct(
    pdi_type,
    disease_group,
    cancer_page,
    outcome_clean,
    gene,
    protein_id,
    OR,
    lower_CI,
    upper_CI,
    PP_H4,
    .keep_all = TRUE
  )

print(
  evidence_context %>%
    count(pdi_type, disease_group, cancer_page, outcome_clean)
)

print(
  evidence_context %>%
    distinct(pdi_type, disease_group, cancer_page, gene, protein_id) %>%
    count(pdi_type, disease_group, cancer_page)
)

# ------------------------------------------------------------
# 6. Shared plot theme and helpers
# ------------------------------------------------------------

small_theme <- theme_bw(base_size = TXT_MAIN) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = TXT_MAIN, margin = margin(t = 2)),
    axis.title.y = element_blank(),
    axis.text = element_text(size = TXT_SMALL),
    axis.ticks = element_line(linewidth = 0.25),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.18),
    plot.title = element_text(size = TXT_MAIN, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = TXT_SMALL, hjust = 0.5),
    plot.margin = margin(3, 3, 3, 3)
  )

blank_panel <- function(label = "Not available") {
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = label,
      size = pt_to_geom(TXT_SMALL)
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void()
}

get_one_dat <- function(plotdat, row_meta) {
  plotdat %>%
    filter(
      pdi_type == row_meta$pdi_type,
      protein_id == row_meta$protein_id
    )
}

clean_int_label <- function(x) {
  ifelse(str_detect(x, "included"), "Interaction included", "No interaction")
}

# ------------------------------------------------------------
# 7. Diet-protein panels
# ------------------------------------------------------------

make_sex_panel <- function(row_meta) {
  
  dat <- get_one_dat(plotdat_sex, row_meta)
  if (nrow(dat) == 0) return(blank_panel())
  
  int_label <- clean_int_label(unique(dat$interaction_label)[1])
  
  ggplot(dat, aes(x = diet_value, y = pred, colour = Sex, fill = Sex)) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.10, colour = NA) +
    geom_line(aes(linetype = Sex), linewidth = 0.85) +
    scale_colour_manual(values = c("Female" = COL_FEMALE, "Male" = COL_MALE)) +
    scale_fill_manual(values = c("Female" = COL_FEMALE, "Male" = COL_MALE)) +
    scale_linetype_manual(values = c("Female" = "solid", "Male" = "dotted")) +
    labs(
      title = "Sex",
      x = paste0("Standardized ", row_meta$pdi_type),
      y = NULL,
      subtitle = int_label
    ) +
    small_theme
}

make_age_panel <- function(row_meta) {
  
  dat <- get_one_dat(plotdat_age, row_meta)
  if (nrow(dat) == 0) return(blank_panel())
  
  int_label <- clean_int_label(unique(dat$interaction_label)[1])
  
  ggplot(dat, aes(x = diet_value, y = pred, colour = age_group, fill = age_group)) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.10, colour = NA) +
    geom_line(aes(linetype = age_group), linewidth = 0.85) +
    scale_colour_manual(
      values = c(
        "Age 50" = COL_AGE50,
        "Age 55" = COL_AGE55,
        "Age 60" = COL_AGE60
      )
    ) +
    scale_fill_manual(
      values = c(
        "Age 50" = COL_AGE50,
        "Age 55" = COL_AGE55,
        "Age 60" = COL_AGE60
      )
    ) +
    scale_linetype_manual(
      values = c(
        "Age 50" = "dotted",
        "Age 55" = "dashed",
        "Age 60" = "solid"
      )
    ) +
    labs(
      title = "Age",
      x = paste0("Standardized ", row_meta$pdi_type),
      y = NULL,
      subtitle = int_label
    ) +
    small_theme
}

make_smoking_panel <- function(row_meta) {
  
  dat <- get_one_dat(plotdat_smoking, row_meta)
  if (nrow(dat) == 0) return(blank_panel())
  
  int_label <- clean_int_label(unique(dat$interaction_label)[1])
  
  ggplot(dat, aes(x = diet_value, y = pred, colour = smoke_group, fill = smoke_group)) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.10, colour = NA) +
    geom_line(aes(linetype = smoke_group), linewidth = 0.85) +
    scale_colour_manual(
      values = c(
        "Never smoked" = COL_NEVERSMOKE,
        "Ever Smoked"  = COL_EVERSMOKE
      )
    ) +
    scale_fill_manual(
      values = c(
        "Never smoked" = COL_NEVERSMOKE,
        "Ever Smoked"  = COL_EVERSMOKE
      )
    ) +
    scale_linetype_manual(
      values = c(
        "Never smoked" = "solid",
        "Ever Smoked" = "dotted"
      )
    ) +
    labs(
      title = "Smoking",
      x = paste0("Standardized ", row_meta$pdi_type),
      y = NULL,
      subtitle = int_label
    ) +
    small_theme
}

make_alcohol_panel <- function(row_meta) {
  
  dat <- get_one_dat(plotdat_alcohol, row_meta)
  if (nrow(dat) == 0) return(blank_panel())
  
  int_label <- clean_int_label(unique(dat$interaction_label)[1])
  
  dat <- dat %>%
    mutate(
      alc_group = factor(
        alc_group,
        levels = c("Never", "Occasional", "Moderate", "Regular/Daily")
      )
    )
  
  ggplot(dat, aes(x = diet_value, y = pred, colour = alc_group, fill = alc_group)) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.08, colour = NA) +
    geom_line(aes(linetype = alc_group), linewidth = 0.85) +
    scale_colour_manual(
      breaks = c("Never", "Occasional", "Moderate", "Regular/Daily"),
      values = c(
        "Never"         = COL_NEVERALC,
        "Occasional"    = COL_OCCALC,
        "Moderate"      = COL_MODALC,
        "Regular/Daily" = COL_REGALC
      )
    ) +
    scale_fill_manual(
      breaks = c("Never", "Occasional", "Moderate", "Regular/Daily"),
      values = c(
        "Never"         = COL_NEVERALC,
        "Occasional"    = COL_OCCALC,
        "Moderate"      = COL_MODALC,
        "Regular/Daily" = COL_REGALC
      )
    ) +
    scale_linetype_manual(
      breaks = c("Never", "Occasional", "Moderate", "Regular/Daily"),
      values = c(
        "Never" = "solid",
        "Occasional" = "dashed",
        "Moderate" = "dotted",
        "Regular/Daily" = "dotdash"
      )
    ) +
    labs(
      title = "Alcohol",
      x = paste0("Standardized ", row_meta$pdi_type),
      y = NULL,
      subtitle = int_label
    ) +
    small_theme
}

# ------------------------------------------------------------
# 8. Legend helpers
# ------------------------------------------------------------

legend_theme <- theme_void(base_size = TXT_MAIN) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = TXT_MAIN),
    legend.text = element_text(size = TXT_SMALL),
    legend.key.width = unit(0.7, "cm"),
    legend.key.height = unit(0.4, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  )

extract_legend <- function(p) {
  cowplot::get_legend(p + theme(legend.position = "bottom"))
}

make_sex_legend <- function() {
  
  dat <- tibble(
    x = rep(1:2, 2),
    y = rep(1:2, each = 2),
    Sex = factor(
      c("Male", "Male", "Female", "Female"),
      levels = c("Male", "Female")
    )
  )
  
  p <- ggplot(dat, aes(x, y, colour = Sex, linetype = Sex)) +
    geom_line(linewidth = 0.65) +
    scale_colour_manual(values = c("Male" = COL_MALE, "Female" = COL_FEMALE)) +
    scale_linetype_manual(values = c("Male" = "dotted", "Female" = "solid")) +
    guides(
      colour = guide_legend(title = "Sex", nrow = 1),
      linetype = guide_legend(title = "Sex", nrow = 1)
    ) +
    legend_theme
  
  cowplot::ggdraw(extract_legend(p))
}

make_age_legend <- function() {
  
  dat <- tibble(
    x = rep(1:2, 3),
    y = rep(1:3, each = 2),
    age_group = factor(
      rep(c("Age 50", "Age 55", "Age 60"), each = 2),
      levels = c("Age 50", "Age 55", "Age 60")
    )
  )
  
  p <- ggplot(dat, aes(x, y, colour = age_group, linetype = age_group)) +
    geom_line(linewidth = 0.65) +
    scale_colour_manual(
      values = c(
        "Age 50" = COL_AGE50,
        "Age 55" = COL_AGE55,
        "Age 60" = COL_AGE60
      )
    ) +
    scale_linetype_manual(
      values = c(
        "Age 50" = "dotted",
        "Age 55" = "dashed",
        "Age 60" = "solid"
      )
    ) +
    guides(
      colour = guide_legend(title = "Age", nrow = 1),
      linetype = guide_legend(title = "Age", nrow = 1)
    ) +
    legend_theme
  
  cowplot::ggdraw(extract_legend(p))
}

make_smoking_legend <- function() {
  
  dat <- tibble(
    x = rep(1:2, 2),
    y = rep(1:2, each = 2),
    smoke_group = factor(
      c("Never smoked", "Never smoked", "Ever Smoked", "Ever Smoked"),
      levels = c("Never smoked", "Ever Smoked")
    )
  )
  
  p <- ggplot(dat, aes(x, y, colour = smoke_group, linetype = smoke_group)) +
    geom_line(linewidth = 0.65) +
    scale_colour_manual(
      values = c(
        "Never smoked" = COL_NEVERSMOKE,
        "Ever Smoked"  = COL_EVERSMOKE
      )
    ) +
    scale_linetype_manual(
      values = c(
        "Never smoked" = "solid",
        "Ever Smoked" = "dotted"
      )
    ) +
    guides(
      colour = guide_legend(title = "Smoking", nrow = 1),
      linetype = guide_legend(title = "Smoking", nrow = 1)
    ) +
    legend_theme
  
  cowplot::ggdraw(extract_legend(p))
}

make_alcohol_legend <- function() {
  
  dat <- tibble(
    x = rep(1:2, 4),
    y = rep(1:4, each = 2),
    alc_group = factor(
      rep(c("Never", "Occasional", "Moderate", "Regular/Daily"), each = 2),
      levels = c("Never", "Occasional", "Moderate", "Regular/Daily")
    )
  )
  
  p <- ggplot(dat, aes(x, y, colour = alc_group, linetype = alc_group)) +
    geom_line(linewidth = 0.65) +
    scale_colour_manual(
      values = c(
        "Never"         = COL_NEVERALC,
        "Occasional"    = COL_OCCALC,
        "Moderate"      = COL_MODALC,
        "Regular/Daily" = COL_REGALC
      )
    ) +
    scale_linetype_manual(
      values = c(
        "Never" = "solid",
        "Occasional" = "dashed",
        "Moderate" = "dotted",
        "Regular/Daily" = "dotdash"
      )
    ) +
    guides(
      colour = guide_legend(title = "Alcohol", nrow = 1),
      linetype = guide_legend(title = "Alcohol", nrow = 1)
    ) +
    legend_theme
  
  cowplot::ggdraw(extract_legend(p))
}

# ------------------------------------------------------------
# 9. Text and MR panels
# ------------------------------------------------------------

text_panel <- function(label, size = TXT_MAIN, bold = FALSE) {
  ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0.5,
      label = label,
      hjust = 0,
      vjust = 0.5,
      size = pt_to_geom(size),
      fontface = ifelse(bold, "bold", "plain"),
      lineheight = 0.95
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void() +
    theme(plot.margin = margin(2, 2, 2, 4))
}

make_mr_panel <- function(mr_dat, x_limits = c(0.3, 2.5)) {
  
  if (nrow(mr_dat) == 0) return(blank_panel())
  
  n_rows <- nrow(mr_dat)
  block_spacing <- 1.55
  
  mr_dat <- mr_dat %>%
    arrange(outcome_clean, desc(PP_H4)) %>%
    mutate(
      y_pos = rev(seq(from = 1, by = block_spacing, length.out = n_rows)),
      mr_text = paste0(
        "Outcome: ", as.character(outcome_clean),
        "\nOR (95% CI): ", sprintf("%.2f", OR),
        " (", sprintf("%.2f", lower_CI), "-", sprintf("%.2f", upper_CI), ")",
        "\nPP.H4: ", sprintf("%.2f", PP_H4)
      )
    )
  
  y_max <- max(mr_dat$y_pos) + 0.65
  
  sep_df <- tibble(y = numeric(0))
  if (n_rows > 1) {
    y_sorted <- sort(mr_dat$y_pos)
    sep_df <- tibble(y = head(y_sorted, -1) + block_spacing / 2)
  }
  
  ggplot(mr_dat) +
    annotate(
      "rect",
      xmin = 0, xmax = 1,
      ymin = 0.35, ymax = y_max,
      fill = "white",
      colour = "grey70",
      linewidth = 0.3
    ) +
    geom_text(
      aes(
        x = 0.03,
        y = y_pos,
        label = mr_text
      ),
      hjust = 0,
      vjust = 0.5,
      size = pt_to_geom(TXT_SMALL),
      lineheight = 0.95
    ) +
    geom_hline(
      data = sep_df,
      aes(yintercept = y),
      linewidth = 0.15,
      colour = "grey85"
    ) +
    coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0.35, y_max),
      clip = "off"
    ) +
    labs(title = "Protein to outcome") +
    theme_void(base_size = TXT_MAIN) +
    theme(
      plot.title = element_text(
        size = TXT_MAIN,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 3)
      ),
      plot.margin = margin(3, 3, 3, 3)
    )
}

# ------------------------------------------------------------
# 10. Assemble one protein row
# ------------------------------------------------------------

make_context_row <- function(row_meta, evidence_dat, x_limits = c(0.3, 2.5)) {
  
  page_val <- as.character(row_meta$cancer_page[[1]])
  
  mr_dat <- evidence_dat %>%
    filter(
      pdi_type == row_meta$pdi_type,
      disease_group == row_meta$disease_group,
      protein_id == row_meta$protein_id
    )
  
  if (!is.na(page_val)) {
    mr_dat <- mr_dat %>%
      filter(cancer_page == page_val)
  }
  
  protein_lab <- paste0(
    row_meta$gene,
    "\n",
    row_meta$protein_id
  )
  
  p_label   <- text_panel(protein_lab, size = TXT_MAIN, bold = TRUE)
  p_sex     <- make_sex_panel(row_meta)
  p_age     <- make_age_panel(row_meta)
  p_smoking <- make_smoking_panel(row_meta)
  p_alcohol <- make_alcohol_panel(row_meta)
  p_mr      <- make_mr_panel(mr_dat, x_limits = x_limits)
  
  p_label + p_sex + p_age + p_smoking + p_alcohol + p_mr +
    plot_layout(
      nrow = 1,
      widths = c(0.75, 1.35, 1.35, 1.35, 1.35, 0.95)
    )
}

# ------------------------------------------------------------
# 11. Build contextualisation figure
# ------------------------------------------------------------

make_context_figure <- function(evidence_dat,
                                pdi,
                                disease_group_in,
                                cancer_page_in = NULL,
                                title = NULL,
                                x_limits = c(0.3, 2.5)) {
  
  dat <- evidence_dat %>%
    filter(
      pdi_type == pdi,
      disease_group == disease_group_in
    )
  
  if (!is.null(cancer_page_in)) {
    dat <- dat %>%
      filter(cancer_page == cancer_page_in)
  }
  
  if (nrow(dat) == 0) {
    stop("No rows available for this pdi/disease group/page after filtering.")
  }
  
  row_meta <- dat %>%
    group_by(pdi_type, disease_group, cancer_page, gene, protein_id) %>%
    summarise(
      first_outcome_order = min(as.numeric(outcome_clean), na.rm = TRUE),
      max_pph4 = max(PP_H4, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(first_outcome_order, desc(max_pph4), gene, protein_id)
  
  if (is.null(title)) {
    title <- paste0(
      pdi,
      "-related proteins with colocalized evidence for ",
      disease_group_in
    )
  }
  
  legend_row <- plot_spacer() +
    make_sex_legend() +
    make_age_legend() +
    make_smoking_legend() +
    make_alcohol_legend() +
    plot_spacer() +
    plot_layout(
      nrow = 1,
      widths = c(0.75, 1.3, 1.3, 1.3, 1.90, 1)
    )
  
  rows <- map(
    seq_len(nrow(row_meta)),
    ~ make_context_row(row_meta[.x, ], evidence_dat = dat, x_limits = x_limits)
  )
  
  body <- wrap_plots(plotlist = rows, ncol = 1)
  
  fig <- legend_row / body +
    plot_layout(heights = c(0.45, nrow(row_meta))) +
    plot_annotation(
      title = title,
      subtitle = paste0(
        "Panels show adjusted diet-protein prediction curves by sex, age, smoking, and alcohol intake. ",
        "The right panel summarizes Mendelian randomization odds ratios with 95% CIs and colocalization posterior probabilities. ",
        "Only results with PP.H4 >= ", min_pph4, " are shown."
      ),
      theme = theme(
        plot.title = element_text(size = TXT_TITLE, face = "bold"),
        plot.subtitle = element_text(size = TXT_MAIN, hjust = 0)
      )
    )
  
  attr(fig, "n_rows") <- nrow(row_meta)
  attr(fig, "suggested_width") <- 22
  attr(fig, "suggested_height") <- 1.8 + 2.55 * nrow(row_meta)
  
  fig
}

# ------------------------------------------------------------
# 12. Preview and save helpers
# ------------------------------------------------------------

preview_png <- function(fig, filename = "preview_contextualisation.png", dpi = 120) {
  
  w <- attr(fig, "suggested_width")
  h <- attr(fig, "suggested_height")
  
  if (is.null(w)) w <- 22
  if (is.null(h)) h <- 10
  
  preview_file <- file.path(tempdir(), filename)
  
  ggsave(
    filename = preview_file,
    plot = fig,
    width = w,
    height = h,
    units = "in",
    dpi = dpi,
    bg = "white",
    limitsize = FALSE
  )
  
  shell.exec(preview_file)
}

save_context_figure_png <- function(fig, filename) {
  
  w <- attr(fig, "suggested_width")
  h <- attr(fig, "suggested_height")
  
  if (is.null(w)) w <- 22
  if (is.null(h)) h <- 10
  
  ggsave(
    filename = file.path(context_out_dir, filename),
    plot = fig,
    width = w,
    height = h,
    units = "in",
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
}

# ------------------------------------------------------------
# 13. Create figures
# ------------------------------------------------------------
fig_hpdi_overall_other_cancers <- make_context_figure(
  evidence_dat = evidence_context,
  pdi = "hPDI",
  disease_group_in = "Cancer",
  cancer_page_in = "Overall and other site-specific cancers",
  title = "hPDI-related proteins with colocalized evidence for all cancer and selected site-specific cancers",
  x_limits = c(0.3, 2.5)
)

fig_hpdi_colorectal <- make_context_figure(
  evidence_dat = evidence_context,
  pdi = "hPDI",
  disease_group_in = "Cancer",
  cancer_page_in = "Colorectal cancer",
  title = "hPDI-related proteins with colocalized evidence for colorectal cancer",
  x_limits = c(0.3, 2.5)
)

fig_updi_overall_other_cancers <- make_context_figure(
  evidence_dat = evidence_context,
  pdi = "uPDI",
  disease_group_in = "Cancer",
  cancer_page_in = "Overall and other site-specific cancers",
  title = "uPDI-related proteins with colocalized evidence for all cancer and selected site-specific cancers",
  x_limits = c(0.3, 2.5)
)

fig_updi_colorectal <- make_context_figure(
  evidence_dat = evidence_context,
  pdi = "uPDI",
  disease_group_in = "Cancer",
  cancer_page_in = "Colorectal cancer",
  title = "uPDI-related proteins with colocalized evidence for colorectal cancer",
  x_limits = c(0.3, 2.5)
)

fig_hpdi_t2d <- make_context_figure(
  evidence_dat = evidence_context,
  pdi = "hPDI",
  disease_group_in = "T2D",
  title = "hPDI-related proteins with colocalized evidence for type 2 diabetes",
  x_limits = c(0.3, 2.5)
)

fig_hpdi_cvd <- make_context_figure(
  evidence_dat = evidence_context,
  pdi = "hPDI",
  disease_group_in = "CVD",
  title = "hPDI-related proteins with colocalized evidence for cardiovascular disease",
  x_limits = c(0.7, 1.3)
)

fig_updi_cvd <- make_context_figure(
  evidence_dat = evidence_context,
  pdi = "uPDI",
  disease_group_in = "CVD",
  title = "uPDI-related proteins with colocalized evidence for cardiovascular disease",
  x_limits = c(0.7, 1.3)
)

# ------------------------------------------------------------
# 14. Preview one figure at a time
# ------------------------------------------------------------

preview_png(fig_hpdi_overall_other_cancers, "preview_hpdi_overall_other_cancers.png")

# preview_png(fig_hpdi_colorectal, "preview_hpdi_colorectal.png")
# preview_png(fig_updi_overall_other_cancers, "preview_updi_overall_other_cancers.png")
# preview_png(fig_updi_colorectal, "preview_updi_colorectal.png")
# preview_png(fig_hpdi_t2d, "preview_hpdi_t2d.png")
# preview_png(fig_hpdi_cvd, "preview_hpdi_cvd.png")
# preview_png(fig_updi_cvd, "preview_updi_cvd.png")

# ------------------------------------------------------------
# 15. Save as PNG only after checking previews
# ------------------------------------------------------------

 save_context_figure_png(
   fig_hpdi_overall_other_cancers,
   "Figure_S1a_hPDI_all_breast_prostate_endometrial_contextualisation.png"
 )

 save_context_figure_png(
   fig_hpdi_colorectal,
   "Figure_S1b_hPDI_colorectal_contextualisation.png"
)

save_context_figure_png(
   fig_updi_overall_other_cancers,
   "Figure_S4a_uPDI_all_breast_prostate_endometrial_contextualisation.png"
 )

 save_context_figure_png(
   fig_updi_colorectal,
   "Figure_S4b_uPDI_colorectal_contextualisation.png"
 )

 save_context_figure_png(
 fig_hpdi_t2d,
   "Figure_S2_hPDI_T2D_contextualisation.png"
 )

save_context_figure_png(
  fig_hpdi_cvd,
   "Figure_S3_hPDI_CVD_contextualisation.png"
 )

save_context_figure_png(
   fig_updi_cvd,
   "Figure_S5_uPDI_CVD_contextualisation.png"
 )

# shell.exec(context_out_dir)
