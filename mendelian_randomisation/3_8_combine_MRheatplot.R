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

# =========================================================
# Combined MR heat map for GBM + ENET
# - scans all MR Excel files
# - keeps selected primary GWAS per disease
# - combines GBM and ENET into one reporting table
# - creates ONE heat map
# - adds signature-source annotation:
#     GBM only / ENET only / Both
# - added LIMMA 
# =========================================================
library(readxl)
library(dplyr)
library(purrr)
library(stringr)
library(writexl)

proj_root <- MANUSCRIPT_DIR

mr_files <- tibble(
  file = file.path(proj_root, c("hpdi_MR.xlsx", "Updi_MR.xlsx")),
  pdi_type = c("hPDI", "uPDI")
)

learning_files <- tibble(
  file = file.path(
    proj_root,
    c(
      "Hpdi_proteins_one_cis_IV_per_protein_MR_RAW_v2.csv",
      "Updi_proteins_one_cis_IV_per_protein_MR_MANUAL_v1.csv"
    )
  ),
  pdi_type = c("hPDI", "uPDI")
)

# -----------------------------
# 3. Disease names from MR sheets
# -----------------------------
outcome_lookup <- tibble::tribble(
  ~outcome_label,              ~disease,       ~outcome_name,
  "Colorectal_GCST90255675",   "Colorectal",   "Colorectal cancer (GCST90255675)",
  "Kidney_GCST90320057",       "Kidney",       "Kidney cancer (GCST90320057)",
  "Prostate_GCST90274714",     "Prostate",     "Prostate cancer (GCST90274714)",
  "Liver_GCST90475568",        "Liver",        "Liver cancer (GCST90475568)",
  "Lung_GCST004748",           "Lung",         "Lung cancer (GCST004748)",
  "Myeloma_GCST90624747",      "Myeloma",      "Myeloma (GCST90624747)",
  "Gastric_GCST90477175",      "Gastric",      "Gastric cancer (GCST90477175)",
  "Endom_GCST90477175",        "Endometrial",  "Endometrial cancer (GCST006465)",
  "Breast_ieu-a-1126",         "Breast",       "Breast cancer (ieu-a-1126)",
  "ieu-b-4966",                "Cancer",       "All cancer / Outcome ieu-b-4966",
  "T2D_GCST006867",            "T2D",          "Type 2 diabetes (ebi-a-GCST006867)",
  "CVD_GCST90029019",          "CVD",          "Cardiovascular disease (ebi-a-GCST90029019)"
)

get_disease <- function(sheet_name) {
  disease <- disease_names[str_detect(
    sheet_name,
    regex(paste0("^", disease_names, collapse = "|^"), ignore_case = TRUE)
  )][1]
  
  ifelse(is.na(disease), str_remove(sheet_name, "_GCST.*$"), disease)
}

extract_mr <- function(file, pdi_type) {
  
  sheets <- excel_sheets(file)
  
  map_dfr(sheets, function(sh) {
    
    dat <- read_excel(file, sheet = sh)
    
    names(dat) <- names(dat) %>%
      make.names(unique = TRUE) %>%
      str_to_lower()
    
    if (!all(c("exposure", "p_fdr") %in% names(dat))) {
      return(tibble())
    }
    
    dat %>%
      mutate(
        p_fdr = as.numeric(p_fdr),
        pdi_type = pdi_type,
        outcome_label = sh,
        protein = exposure,
        protein_match = str_to_upper(str_trim(exposure))
      ) %>%
      left_join(outcome_lookup, by = "outcome_label") %>%
      mutate(
        disease = if_else(is.na(disease), outcome_label, disease),
        outcome_name = if_else(is.na(outcome_name), outcome_label, outcome_name)
      ) %>%
      filter(!is.na(p_fdr)) %>%
      relocate(
        pdi_type,
        disease,
        outcome_label,
        outcome_name,
        protein,
        exposure,
        protein_match
      )
  })
}

mr_combined <- pmap_dfr(
  mr_files,
  ~ extract_mr(file = ..1, pdi_type = ..2)
)

# ----------------------
# read learning files
# -----------------------
read_learning_file <- function(file, pdi_type) {
  
  dat <- read.csv(
    file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  names(dat) <- names(dat) %>%
    make.names(unique = TRUE) %>%
    str_to_lower()
  
  # Standardise elasticnet column name to enet
  if ("elasticnet" %in% names(dat) & !"enet" %in% names(dat)) {
    dat <- dat %>% rename(enet = elasticnet)
  }
  
  dat %>%
    mutate(
      pdi_type = pdi_type,
      protein_match = str_to_upper(str_trim(gene)),
      gbm = as.numeric(gbm),
      limma = as.numeric(limma),
      enet = as.numeric(enet)
    ) %>%
    group_by(pdi_type, protein_match) %>%
    summarise(
      gbm = max(coalesce(gbm, 0)),
      limma = max(coalesce(limma, 0)),
      enet = max(coalesce(enet, 0)),
      .groups = "drop"
    )
}

learning_combined <- pmap_dfr(
  learning_files,
  ~ read_learning_file(file = ..1, pdi_type = ..2)
)

mr_with_learning <- mr_combined %>%
  left_join(
    learning_combined,
    by = c("pdi_type", "protein_match")
  ) %>%
  mutate(
    gbm = if_else(is.na(gbm), 0, gbm),
    limma = if_else(is.na(limma), 0, limma),
    enet = if_else(is.na(enet), 0, enet)
  ) %>%
  relocate(
    pdi_type,
    disease,
    protein,
    gbm,
    limma,
    enet,
    exposure,
    outcome
  )

# -----------------------------
# 7. Summaries
# -----------------------------
sig_mr_with_learning<- mr_with_learning %>% filter(p_fdr <0.05)

head(sig_mr_with_learning) #112 significant rows

summary_by_pdi_disease <- sig_mr_with_learning %>%
  count(pdi_type, disease, name = "n_significant") %>%
  arrange(pdi_type, desc(n_significant))

summary_by_pdi_disease
# A tibble: 15 × 3
# pdi_type disease     n_significant
# <chr>    <chr>               <int>
#   1 hPDI     Prostate               23
# 2 hPDI     CVD                    17
# 3 hPDI     Colorectal              8
# 4 hPDI     T2D                     6
# 5 hPDI     Breast                  5
# 6 hPDI     Cancer                  4
# 7 hPDI     Endometrial             2
# 8 hPDI     Liver                   1
# 9 hPDI     Lung                    1
# 10 uPDI     CVD                    18
# 11 uPDI     Colorectal             12
# 12 uPDI     Breast                  7
# 13 uPDI     Prostate                6
# 14 uPDI     Cancer                  1
# 15 uPDI     Endometrial             1

summary_by_learning_method <- sig_mr_with_learning %>%
  summarise(
    n_gbm   = sum(gbm == 1, na.rm = TRUE),
    n_limma = sum(limma == 1, na.rm = TRUE),
    n_enet  = sum(enet == 1, na.rm = TRUE),
    .by = c(pdi_type, disease)
  ) %>%
  arrange(pdi_type, disease)

summary_by_learning_method
# pdi_type disease     n_gbm n_limma n_enet
# <chr>    <chr>       <int>   <int>  <int>
#   1 hPDI     Breast          4       1      0
# 2 hPDI     CVD            14       6      2
# 3 hPDI     Cancer          2       3      1
# 4 hPDI     Colorectal      6       2      0
# 5 hPDI     Endometrial     2       0      0
# 6 hPDI     Liver           0       1      0
# 7 hPDI     Lung            1       0      0
# 8 hPDI     Prostate       16      11      0
# 9 hPDI     T2D             4       2      0
# 10 uPDI     Breast          7       1      1
# 11 uPDI     CVD            18       3      2
# 12 uPDI     Cancer          1       0      0
# 13 uPDI     Colorectal     12       0      0
# 14 uPDI     Endometrial     1       1      0
# 15 uPDI     Prostate        6       0      0

unique_sig_proteins_with_learning <- sig_mr_with_learning %>%
  distinct(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    gbm,
    limma,
    enet
  ) %>%
  arrange(pdi_type, disease, protein)

sig_mr_with_learning <- sig_mr_with_learning %>%
  mutate(
    learning_signature = case_when(
      gbm == 1 & limma == 1 & enet == 1 ~ "gbm + limma + enet",
      gbm == 1 & limma == 1 ~ "gbm + limma",
      gbm == 1 & enet == 1 ~ "gbm + enet",
      limma == 1 & enet == 1 ~ "limma + enet",
      gbm == 1 ~ "gbm",
      limma == 1 ~ "limma",
      enet == 1 ~ "enet",
      TRUE ~ "not found"
    )
  )

table(sig_mr_with_learning$learning_signature, useNA="ifany")
# gbm        gbm + limma gbm + limma + enet              limma 
# 81                  7                  6                 18 

summary_by_learning_signature <- sig_mr_with_learning %>%
  count(pdi_type, disease, learning_signature, name = "n_significant") %>%
  arrange(pdi_type, disease, desc(n_significant))

summary_by_learning_signature

# -----------------------------
# 8. Export
# -----------------------------

write.xlsx(
  list(
    all_mr_with_learning = mr_with_learning,
    significant_MR_with_learning = sig_mr_with_learning,
    summary_by_pdi_disease = summary_by_pdi_disease,
    summary_by_learning_method = summary_by_learning_method,
    summary_by_learning_signature = summary_by_learning_signature
  ),
  file.path(proj_root, "significant_MR_results_hPDI_uPDI_with_learning_indicators.xlsx")
)

## opening file
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(patchwork)

file <- file.path(proj_root, "significant_MR_results_hPDI_uPDI_with_learning_indicators.xlsx")

mr_with_learning_all <- read_excel(file, sheet = "all_mr_with_learning")

names(mr_with_learning_all) <- names(mr_with_learning_all) %>%
  make.names(unique = TRUE) %>%
  str_to_lower()


##PLOTTING ---------------------
# 8) Row-level source annotation within diet_score
# -----------------------------
row_source <- mr_with_learning_all %>%
  filter(!is.na(exposure)) %>%
  distinct(pdi_type, exposure, gbm, limma, enet) %>%
  transmute(
    diet_score = pdi_type,
    exposure,
    XGB   = as.integer(gbm == 1),
    ENET  = as.integer(enet == 1),
    LIMMA = as.integer(limma == 1)
  )


#--------------------------SUPPLEMENT---------------
# extract_mr <- function(file, pdi_type) {
#   
#   sheets <- excel_sheets(file)
#   
#   map_dfr(sheets, function(sh) {
#     
#     dat <- read_excel(file, sheet = sh)
#     
#     names(dat) <- names(dat) %>%
#       make.names(unique = TRUE) %>%
#       str_to_lower()
#     
#     if (!all(c("exposure", "p_fdr") %in% names(dat))) {
#       return(tibble())
#     }
#     
#     dat %>%
#       mutate(
#         p_fdr = as.numeric(p_fdr),
#         pdi_type = pdi_type,
#         outcome_label = sh,
#         protein = exposure,
#         protein_match = str_to_upper(str_trim(exposure))
#       ) %>%
#       left_join(outcome_lookup, by = "outcome_label") %>%
#       mutate(
#         disease = if_else(is.na(disease), outcome_label, disease),
#         outcome_name = if_else(is.na(outcome_name), outcome_label, outcome_name)
#       ) %>%
#       filter(!is.na(p_fdr)) %>%
#       relocate(
#         pdi_type,
#         disease,
#         outcome_label,
#         outcome_name,
#         protein,
#         exposure,
#         protein_match
#       )
#   })
# }
# 
# mr_all <- pmap_dfr(
#   mr_files,
#   ~ extract_mr(file = ..1, pdi_type = ..2)
# )

# ----------------------
# read learning files
# -----------------------
# mr_with_learning_all <- mr_all %>%
#   left_join(
#     learning_combined,
#     by = c("pdi_type", "protein_match")
#   ) %>%
#   mutate(
#     gbm = if_else(is.na(gbm), 0, gbm),
#     limma = if_else(is.na(limma), 0, limma),
#     enet = if_else(is.na(enet), 0, enet)
#   ) %>%
#   relocate(
#     pdi_type,
#     disease,
#     protein,
#     gbm,
#     limma,
#     enet,
#     exposure,
#     outcome
#   )

##keeping all the oucomes for any one protein with less than fdr 0.01
fdr_cutoff <- 0.1   # change to 0.10 if you meant FDR < 0.10

proteins_keep <- mr_with_learning_all %>%
  filter(!is.na(p_fdr), p_fdr < fdr_cutoff) %>%
  distinct(pdi_type, protein)

mr_heat_all_outcomes <- mr_with_learning_all %>%
  semi_join(
    proteins_keep,
    by = c("pdi_type", "protein")
  )

mr_heat_all_outcomes %>%
  count(pdi_type, protein, name = "n_outcomes") %>%
  arrange(pdi_type, protein)

# -----------------------------
# 11) Wrapper to build stacked hPDI/uPDI figure
# -----------------------------
supp_disease_order <- c(
  "All cancer",
  "Breast cancer",
  "Prostate cancer",
  "Colorectal cancer",
  "Kidney cancer",
  "Liver cancer",
  "Lung cancer",
  "Myeloma",
  "Endometrial cancer",
  "Stomach cancer",
  "T2D",
  "CVD"
)


# -----------------------------
#supp
# -----------------------------

heat_dat0 <- mr_heat_all_outcomes %>%
  mutate(
    diet_score = pdi_type,
    exposure = protein,
    beta = b,
    
    disease_label = recode(
      disease,
      "Cancer" = "All cancer",
      "Breast" = "Breast cancer",
      "Prostate" = "Prostate cancer",
      "Colorectal" = "Colorectal cancer",
      "Kidney" = "Kidney cancer",
      "Liver" = "Liver cancer",
      "Lung" = "Lung cancer",
      "Myeloma" = "Myeloma",
      "Endometrial" = "Endometrial cancer",
      "Endom" = "Endometrial cancer",
      "Gastric" = "Stomach cancer",
      "T2D" = "T2D",
      "CVD" = "CVD",
      .default = disease
    ),
    
    XGB   = as.integer(gbm == 1),
    ENET  = as.integer(enet == 1),
    LIMMA = as.integer(limma == 1)
  )

make_clean_heatmap <- function(heat_dat,
                               diet_score_keep,
                               disease_order,
                               title_text = NULL,
                               subtitle_text = "* FDR < 0.05, ^ FDR < 0.1",
                               star_color = "darkred") {
  
  heat_dat <- heat_dat %>%
    filter(
      diet_score == diet_score_keep,
      disease_label %in% disease_order,
      !is.na(exposure),
      !is.na(beta),
      !is.na(pval),
      !is.na(p_fdr)
    )
  
  if (nrow(heat_dat) == 0) {
    return(NULL)
  }
  
  protein_order <- heat_dat %>%
    group_by(exposure) %>%
    summarise(best_p = min(p_fdr, na.rm = TRUE), .groups = "drop") %>%
    arrange(best_p) %>%
    pull(exposure)
  
  source_annot <- heat_dat %>%
    group_by(exposure) %>%
    summarise(
      XGB   = max(XGB, na.rm = TRUE),
      ENET  = max(ENET, na.rm = TRUE),
      LIMMA = max(LIMMA, na.rm = TRUE),
      .groups = "drop"
    )
  
  heat_full <- expand_grid(
    exposure = protein_order,
    disease_label = disease_order
  ) %>%
    left_join(
      heat_dat %>%
        select(exposure, disease_label, beta, pval, p_fdr),
      by = c("exposure", "disease_label")
    ) %>%
    left_join(source_annot, by = "exposure") %>%
    mutate(
      sig_lab = case_when(
        !is.na(p_fdr) & p_fdr < 0.05 ~ "*",
        !is.na(p_fdr) & p_fdr < 0.10 ~ "^",
        TRUE ~ ""
      ),
      exposure = factor(exposure, levels = rev(protein_order)),
      disease_label = factor(disease_label, levels = disease_order),
      XGB   = if_else(is.na(XGB), 0L, XGB),
      ENET  = if_else(is.na(ENET), 0L, ENET),
      LIMMA = if_else(is.na(LIMMA), 0L, LIMMA)
    )
  
  annot_dat <- heat_full %>%
    distinct(exposure, XGB, ENET, LIMMA) %>%
    pivot_longer(
      cols = c(XGB, ENET, LIMMA),
      names_to = "signature",
      values_to = "selected"
    ) %>%
    mutate(
      signature = factor(signature, levels = c("XGB", "ENET", "LIMMA")),
      selected = factor(selected, levels = c(0, 1))
    )
  
  p_annot <- ggplot(
    annot_dat,
    aes(x = signature, y = exposure, fill = selected)
  ) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_manual(
      values = c(
        "0" = "grey75",
        "1" = "grey25"
      ),
      drop = FALSE
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 8) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none",
      plot.margin = margin(16, 0, 16, 16)
    )
  
  p_heat <- ggplot(
    heat_full,
    aes(x = disease_label, y = exposure, fill = beta)
  ) +
    geom_tile(color = "grey85", linewidth = 0.25) +
    geom_text(aes(label = sig_lab), color = star_color, size = 3.2) +
    scale_fill_gradient2(
      low = "#2B3256",
      mid = "grey100",
      high = "#D27E64",
      midpoint = 0,
      na.value = "grey95",
      name = "Beta"
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = title_text,
      subtitle = subtitle_text
    ) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(face = "bold", size = 10),
      plot.margin = margin(16, 16, 16, 0)
    )
  
  final_plot <- p_annot + p_heat + plot_layout(widths = c(0.65, 8))
  
  list(
    plot = final_plot,
    heat_full = heat_full
  )
}

make_stacked_heatmap <- function(heat_dat,
                                 disease_order,
                                 mode = c("main", "supp")) {
  
  mode <- match.arg(mode)
  
  dat_plot <- heat_dat %>%
    filter(disease_label %in% disease_order)
  
  subtitle_text <- "* FDR < 0.05, ^ FDR < 0.01"
  
  hpdi_res <- make_clean_heatmap(
    heat_dat = dat_plot,
    diet_score_keep = "hPDI",
    disease_order = disease_order,
    title_text = "hPDI-linked proteins and associations with cancers, T2D and CVD",
    subtitle_text = subtitle_text
  )
  
  updi_res <- make_clean_heatmap(
    heat_dat = dat_plot,
    diet_score_keep = "uPDI",
    disease_order = disease_order,
    title_text = "uPDI-linked proteins and associations with cancers, T2D and CVD",
    subtitle_text = subtitle_text
  )
  
  stacked_plot <- hpdi_res$plot / updi_res$plot +
    plot_layout(ncol = 1, guides = "collect") &
    theme(legend.position = "right")
  
  list(
    plot = stacked_plot,
    hpdi = hpdi_res,
    updi = updi_res
  )
}

main_res <- make_stacked_heatmap(
  heat_dat = heat_dat0,
  disease_order = supp_disease_order,
  mode = "main"
)


main_res$hpdi$plot
main_res$updi$plot

ggsave(
  filename = file.path(proj_root, "heatmap_HPDI_SUPP.png"),
  plot = main_res$hpdi$plot,
  width = 21,
  height = 27.9,
  units = "cm",
  dpi = 300
)

ggsave(
  filename = file.path(proj_root, "heatmap_UPDI_SUPP.png"),
  plot = main_res$updi$plot,
  width = 21,
  height = 27.9,
  units = "cm",
  dpi = 300
)


# -----------------------------
# 9) Disease order
# -----------------------------
main_disease_order <- c(
  "All cancer",
  "Breast cancer",
  "Prostate cancer",
  "Lung cancer",
  "Colorectal cancer",
  "Liver cancer",
  "Endometrial cancer",
  "T2D",
  "CVD"
)

##keeping all the oucomes for any one protein with less than fdr 0.01
fdr_cutoff <- 0.05   # change to 0.10 if you meant FDR < 0.10

proteins_keep <- mr_with_learning_all %>%
  filter(!is.na(p_fdr), p_fdr < fdr_cutoff) %>%
  distinct(pdi_type, protein)

mr_heat_all_outcomes <- mr_with_learning_all %>%
  semi_join(
    proteins_keep,
    by = c("pdi_type", "protein")
  )

mr_heat_all_outcomes %>%
  count(pdi_type, protein, name = "n_outcomes") %>%
  arrange(pdi_type, protein)



heat_dat1 <- mr_heat_all_outcomes %>%
  mutate(
    diet_score = pdi_type,
    exposure = protein,
    beta = b,
    
    disease_label = recode(
      disease,
      "Cancer" = "All cancer",
      "Breast" = "Breast cancer",
      "Prostate" = "Prostate cancer",
      "Lung" = "Lung cancer",
      "Colorectal" = "Colorectal cancer",
      "Liver" = "Liver cancer",
      "Endometrial" = "Endometrial cancer",
      "Endom" = "Endometrial cancer",
      "T2D" = "T2D",
      "CVD" = "CVD",
      .default = disease
    ),
    
    XGB   = as.integer(gbm == 1),
    ENET  = as.integer(enet == 1),
    LIMMA = as.integer(limma == 1)
  )


main_res1 <- make_stacked_heatmap(
  heat_dat = heat_dat1,
  disease_order = main_disease_order,
  mode = "main"
)


main_res1$hpdi$plot
main_res1$updi$plot

ggsave(
  filename = file.path(proj_root, "heatmap_sig_HPDI.png"),
  plot = main_res1$hpdi$plot,
  width = 21,
  height = 27.9,
  units = "cm",
  dpi = 300
)

ggsave(
  filename = file.path(proj_root, "heatmap_sig_UPDI.png"),
  plot = main_res1$updi$plot,
  width = 21,
  height = 27.9,
  units = "cm",
  dpi = 300
)
