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

#'---- Syntax for ProPlant Project'#
#' Created by Peggy Ler (Pler) - 20th OCt 2025
#' Filename..: 3a_volcano_hPDI.R
#' Updated...: "" by "Pler"
#' Purpose...: Creating outputs from maineffects models with hdpi
#' 
#' Olink data
#' --------------------------------------------------------------
#' Data used...: UKBB on RAP
#' Using R studio 2024.04.02
#' -- Main program ----------------------------------------------
# --- Load data ----------------------------------------------------------------
proj_dir <- MANUSCRIPT_DIR

df <- read.csv(file.path(proj_dir, "limmalist_updi_forsupp.csv"),
                     stringsAsFactors = FALSE)

# --- Packages -----------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggrepel); library(scales)
})


# --- Settings -----------------------------------------------------------------
# --- Settings -----------------------------------------------------------------
alpha_fdr <- 0.05
p_line_fdr <- alpha_fdr
top_n_overall <- 15
# Prep volcano data
# Prep volcano data
volc <- df %>%
  mutate(
    gene_clean = if_else(is.na(gene) | gene == "", protein_id, gene)
  ) %>%
  group_by(gene_clean) %>%
  mutate(
    n_gene_name = n_distinct(protein_id)
  ) %>%
  ungroup() %>%
  mutate(
    Gene_label = if_else(
      n_gene_name > 1,
      paste0(gene_clean, " (", protein_id, ")"),
      gene_clean
    ),
    neglog10_p_adj = -log10(limma_adj_p_val),
    group = case_when(
      limma_adj_p_val < alpha_fdr ~ "Diet significant",
      TRUE ~ "Not significant"
    )
  )
# Print counts
cat("Diet significant proteins:", sum(volc$group == "Diet significant"), "\n")
#619
# Label
#-------
label_df <- volc %>%
  filter(limma_adj_p_val < alpha_fdr) %>%
  arrange(desc(neglog10_p_adj)) %>%
  slice_head(n = top_n_overall)

cat("Labels created:", nrow(label_df), "\n")

# Scales and Aesthetics
group_levels <- c("Not significant", "Diet significant")
volc$group <- factor(volc$group, levels = group_levels)

cols <- c(
  "Diet significant" = "#D27E64",
  "Not significant" = "grey80"
)

# Set the dashed line at FDR = 0.05 as a visual guide
# ---------------------------------------------------
p_line_fdr <- alpha_fdr

# PLOT
#---------
g_volcano <- ggplot(volc, aes(x = limma_log_fc, y = neglog10_p_adj)) +
  geom_hline(
    yintercept = -log10(p_line_fdr),
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_point(
    aes(color = group),
    size = 1.6,
    alpha = 0.7
  ) +
  geom_text_repel(
    data = label_df,
    aes(
      x = limma_log_fc,
      y = neglog10_p_adj,
      label = Gene_label
    ),
    size = 3,
    max.overlaps = Inf,
    min.segment.length = 0,
    box.padding = 0.5,
    point.padding = 0.5,
    seed = 123,
    show.legend = FALSE
  ) +
  scale_color_manual(values = cols, name = NULL) +
  labs(
    x = "Log fold-change",
    y = expression(-log[10](FDR)),
    title = "Volcano plot of uPDI-associated proteins from limma",
    subtitle = "FDR < 0.05 highlighted",
    caption = paste(
      "Dashed horizontal line: FDR =", alpha_fdr,
      "| Significant proteins:",
      sum(volc$group == "Diet significant")
    )
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30"),
    plot.caption = element_text(size = 8, color = "grey50")
  )

g_volcano

ggsave(
  filename = file.path(proj_dir, "Volcano_uPDI_LIMMA.png"),
  plot = g_volcano,
  width = 7,
  height = 6,
  dpi = 600
)
