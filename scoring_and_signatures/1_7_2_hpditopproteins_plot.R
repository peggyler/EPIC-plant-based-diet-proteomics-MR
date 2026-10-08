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
#' Created by Peggy Ler (Pler) - 31 OCt 2025
#' Filename..: Hpdi_protein_plot.R
#' Updated...: "" by "Pler"
#' Purpose...: To visualise the signatures top proteins and the volcano plot
#' 
#' Using R studio 2024.04.02
#' -- Main program ----------------------------------------------

# Top 10 proteins per hPDI signature in one faceted plot
# =====================================================

library(readr)
library(dplyr)
library(forcats)
library(ggplot2)
library(tidytext)
library(scales)

# Paths
# -----
proj_root <- MANUSCRIPT_DIR
out_dir <- proj_root   # change this if your GBM/ENET files are elsewhere

limma_file <- file.path(proj_root, "limmalist_hpdi_forsupp.csv")   
gbm_file   <- file.path(out_dir,   "Hpdi_XGBoost_proteinlist.csv")  
enet_file  <- file.path(out_dir,   "Hpdi_ENET_proteinlist.csv")  

pick_col <- function(dat, candidates, required = TRUE) {
  
  hit <- candidates[candidates %in% names(dat)][1]
  
  if (is.na(hit) && required) {
    stop(
      "Could not find any of these columns: ",
      paste(candidates, collapse = ", "),
      "\nAvailable columns are:\n",
      paste(names(dat), collapse = ", ")
    )
  }
  
  hit
}

parse_num <- function(x) {
  suppressWarnings(readr::parse_number(as.character(x)))
}
# Helper functions
# ----------------
read_signature <- function(file, signature, method = c("limma", "enet", "xgb")) {
  
  method <- match.arg(method)
  dat <- read_csv(file, show_col_types = FALSE)
  
  gene_col <- pick_col(
    dat,
    c("gene", "Gene"),
    required = FALSE
  )
  
  feature_col <- pick_col(
    dat,
    c("protein_id", "seq_id", "seq_id_de", "Feature"),
    required = FALSE
  )
  
  # ----------------------------
  # Signed value + ranking metric
  # ----------------------------
  
  if (method == "limma") {
    
    effect_col <- pick_col(
      dat,
      c("limma_log_fc"),
      required = TRUE
    )
    
    fdr_col <- pick_col(
      dat,
      c("limma_adj_p_val"),
      required = TRUE
    )
    
    dat <- dat %>%
      mutate(
        .limma_log_fc = parse_num(.data[[effect_col]]),
        .limma_fdr    = parse_num(.data[[fdr_col]])
      ) %>%
      filter(
        !is.na(.limma_log_fc),
        !is.na(.limma_fdr),
        .limma_fdr < 0.05
      )
    
    signed_value <- dat$.limma_log_fc
    rank_metric  <- abs(dat$.limma_log_fc)
    
  } else if (method == "enet") {
    
    effect_col <- pick_col(dat, c("elasticnet_coef"), required = TRUE)
    rank_col   <- pick_col(dat, c("elasticnet_abs_coef"), required = TRUE)
    
    signed_value <- parse_num(dat[[effect_col]])
    rank_metric  <- parse_num(dat[[rank_col]])
    
  } else if (method == "xgb") {
    
    rank_col <- pick_col(
      dat,
      c("mean_abs_shap_test", "relative_mean_abs_shap_test"),
      required = TRUE
    )
    
    direction_col <- pick_col(
      dat,
      c("average_shap_direction"),
      required = TRUE
    )
    
    rank_metric <- parse_num(dat[[rank_col]])
    
    direction_sign <- case_when(
      trimws(as.character(dat[[direction_col]])) == "+" ~  1,
      trimws(as.character(dat[[direction_col]])) == "-" ~ -1,
      TRUE ~ NA_real_
    )
    
    signed_value <- rank_metric * direction_sign
  }
  
  # ----------------------------
  # Protein labels
  # ----------------------------
  
  gene <- if (!is.na(gene_col)) {
    as.character(dat[[gene_col]])
  } else {
    rep(NA_character_, nrow(dat))
  }
  
  feature <- if (!is.na(feature_col)) {
    as.character(dat[[feature_col]])
  } else {
    rep(NA_character_, nrow(dat))
  }
  
  gene_label <- coalesce(
    na_if(gene, ""),
    na_if(feature, "")
  )
  
  tibble(
    signature = signature,
    method = method,
    Gene_label = gene_label,
    signed_value = signed_value,
    rank_metric = rank_metric
  )
}

# Read and combine
# ----------------
plot_df <- bind_rows(
  read_signature(limma_file, "hPDI-limma", method = "limma"),
  read_signature(enet_file,  "hPDI-ENET",  method = "enet"),
  read_signature(gbm_file,   "hPDI-XGB",   method = "xgb")
) %>%
  filter(
    !is.na(Gene_label),
    !is.na(signed_value),
    !is.na(rank_metric)
  ) %>%
  distinct(signature, Gene_label, .keep_all = TRUE) %>%
  group_by(signature) %>%
  slice_max(order_by = rank_metric, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    signature_label = case_when(
      signature == "hPDI-limma" ~ "hPDI-limma\nsigned log fold-change",
      signature == "hPDI-ENET"  ~ "hPDI-ENET\nsigned elastic-net coefficient",
      signature == "hPDI-XGB"   ~ "hPDI-XGB\ndirection-coded mean absolute SHAP",
      TRUE ~ signature
    ),
    signature_label = factor(
      signature_label,
      levels = c(
        "hPDI-limma\nsigned log fold-change",
        "hPDI-ENET\nsigned elastic-net coefficient",
        "hPDI-XGB\ndirection-coded mean absolute SHAP"
      )
    )
  ) %>%
  group_by(signature_label) %>%
  mutate(
    rank_rel = ifelse(
      n_distinct(rank_metric) > 1,
      scales::rescale(rank_metric, to = c(1, 10)),
      5
    ),
    Gene_label_plot = reorder_within(Gene_label, signed_value, signature_label)
  ) %>%
  ungroup()

# Plot
# ----
p <- ggplot(
  plot_df,
  aes(
    x = signed_value,
    y = Gene_label_plot,
    size = rank_rel,
    color = signed_value
  )
) +
  geom_point(alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~ signature_label, scales = "free", nrow = 1) +
  scale_y_reordered() +
  scale_color_gradient2(
    low = "#2B3256",
    mid = "#EFEDDF",
    high = "#D27E64",
    midpoint = 0
  ) +
  scale_size_continuous(range = c(3, 7), guide = "none") +
  labs(
    x = "Signed model-specific value",
    y = NULL,
    color = "Direction",
    title = "Top 10 hPDI-related proteins from each signature",
    subtitle = "Proteins were ranked by absolute log fold-change for LIMMA, absolute coefficient for ENET, and mean absolute SHAP in the test set for XGB"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold")
  )

p



# Make direction categorical, so the color legend is clearer
library(patchwork)
library(grid)
max_abs_color <- max(abs(plot_df2$signed_value), na.rm = TRUE)
plot_df2 <- plot_df %>%
  mutate(
    signature = factor(
      signature,
      levels = c("hPDI-limma", "hPDI-ENET", "hPDI-XGB")
    )
  )

make_signature_plot <- function(dat, sig, x_lab) {
  
  dat_sub <- dat %>%
    filter(signature == sig) %>%
    mutate(
      Gene_label_plot = fct_reorder(Gene_label, signed_value)
    )
  
  ggplot(
    dat_sub,
    aes(
      x = signed_value,
      y = Gene_label_plot,
      color = signed_value
    )
  ) +
    geom_point(size = 4, alpha = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
    scale_color_gradient2(
      low = "#2B3256",
      mid = "#EFEDDF",
      high = "#D27E64",
      midpoint = 0,
      limits = c(-max_abs_color, max_abs_color),
      name = NULL,
      breaks = c(-max_abs_color, 0, max_abs_color),
      labels = scales::label_number(accuracy = 0.01),
      guide = guide_colorbar(
        barwidth = unit(7, "cm"),
        barheight = unit(0.35, "cm")
      )
    ) +
    labs(
      title = sig,
      x = x_lab,
      y = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 9)
    )
}

p_limma <- make_signature_plot(
  plot_df2,
  sig = "hPDI-limma",
  x_lab = "Log fold-change"
)

p_enet <- make_signature_plot(
  plot_df2,
  sig = "hPDI-ENET",
  x_lab = "Elastic-net coefficient"
)

p_xgb <- make_signature_plot(
  plot_df2,
  sig = "hPDI-XGB",
  x_lab = "Mean absolute SHAP"
)

p <- (p_limma + p_enet + p_xgb) +
  plot_layout(nrow = 1, guides = "collect") +
  plot_annotation(
    title = "Top 10 hPDI-related proteins from each signature",
    subtitle = "Proteins were ranked by absolute log fold-change for limma, absolute coefficient for ENET, and mean absolute SHAP in the test set for XGB"
  ) &
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

p


# Save
# ----
ggsave(
  filename = file.path(proj_root, "Top10_hPDI_proteins_each_signature_single_plot.png"),
  plot = p,
  width = 12,
  height = 6,
  dpi = 600
)

## create list of top proteins without cis
##-----------------------------------------
library(readr)
library(dplyr)
library(stringr)

# -----------------------------
# Read cis instrument protein list
# -----------------------------

cis_file <- file.path(
  proj_root,
  "Hpdi_proteins_one_cis_IV_per_protein_MR_RAW_v2.csv"
)

cis_df <- read_csv(cis_file, show_col_types = FALSE)

cis_genes <- cis_df %>%
  filter(!is.na(gene), gene != "") %>%
  distinct(gene) %>%
  mutate(has_cis_instrument = TRUE)

top10_cis_check <- plot_df %>%
  left_join(
    cis_genes,
    by = c("Gene_label" = "gene")
  ) %>%
  mutate(
    has_cis_instrument = if_else(
      is.na(has_cis_instrument),
      FALSE,
      has_cis_instrument
    )
  ) %>%
  arrange(signature, desc(rank_metric))

if (interactive()) View(top10_cis_check)
table(top10_cis_check$has_cis_instrument)
# FALSE  TRUE 
# 6    24 

top10_cis_check %>%
  filter(!has_cis_instrument) %>%
  arrange(signature, desc(rank_metric)) %>%
  select(
    signature,
    Gene_label,
    signed_value,
    rank_metric
  )



# A tibble: 6 × 4
# signature  Gene_label signed_value rank_metric
# <chr>      <chr>             <dbl>       <dbl>
# 1 hPDI-ENET  GULP1            0.0558      0.0558
# 2 hPDI-ENET  CAMK1           -0.0500      0.0500
# 3 hPDI-ENET  UBQLN3           0.0437      0.0437
# 4 hPDI-ENET  RND3            -0.0426      0.0426
# 5 hPDI-XGB   TEAD3           -0.0152      0.0152
# 6 hPDI-limma SDC3            -0.0709      0.0709

write_csv(
  top10_cis_check,
  file.path(proj_root, "Top10_hPDI_signature_proteins_cis_instruments_status.csv")
)
