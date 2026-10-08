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

library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(readxl)

## Create MR results tables for supplement
mr_file <- file.path(MANUSCRIPT_DIR, "significant_MR_results_hPDI_uPDI_with_learning_indicators.xlsx")
# First check exact sheet names
excel_sheets(mr_file)
df_mr_sig <- read_excel(
  mr_file,
  sheet = "all_mr_with_learning"
)

colnames(df_mr_sig)
# [1] "pdi_type"      "disease"       "protein"       "gbm"           "limma"         "enet"         
# [7] "exposure"      "outcome"       "outcome_label" "outcome_name"  "protein_match" "id.exposure"  
# [13] "id.outcome"    "method"        "nsnp"          "b"             "se"            "pval"         
# [19] "p_fdr"         "or"            "or_lci95"      "or_uci95"  

hpdi <- read.csv(file.path(MANUSCRIPT_DIR, "Hpdi_proteins_one_cis_IV_per_protein_MR_RAW_v2.csv"))
updi <- read.csv(file.path(MANUSCRIPT_DIR, "Updi_proteins_one_cis_IV_per_protein_MR_MANUAL_v1.csv"))


## ============================================================
## 1. Cancer-only MR results
## ============================================================

mr_cancer <- df_mr_sig %>%
  filter(
    if_any(
      any_of("disease"),
      ~ str_detect(str_to_lower(as.character(.x)), "cancer")
    )
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_match,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  )

## ============================================================
## 3. annotation look up
## =========
hpdi_lookup <- hpdi %>%
dplyr::mutate(
  pdi_type = "hPDI",
  protein = str_trim(as.character(gene))
) %>%
  dplyr::select(
    pdi_type,
    protein,
    protein_id,
    target,
    uniprot,
    rsid
  ) %>%
  dplyr::distinct()

updi_lookup <- updi %>%
  dplyr::mutate(
    pdi_type = "uPDI",
    protein = str_trim(as.character(gene))
  ) %>%
  dplyr::select(
    pdi_type,
    protein,
    protein_id,
    target,
    uniprot,
    rsid
  ) %>%
  dplyr::distinct()

protein_lookup <- dplyr::bind_rows(hpdi_lookup, updi_lookup)

##merge it
mr_cancer_supp <- mr_cancer %>%
  dplyr::left_join(
    protein_lookup,
    by = c("pdi_type", "protein")
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_id,
    target,
    uniprot,
    rsid,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  ) %>%
  dplyr::arrange(pdi_type, protein, method)

##check
mr_cancer_supp %>%
  dplyr::filter(is.na(protein_id) | is.na(target) | is.na(uniprot) | is.na(rsid)) %>%
  dplyr::count(pdi_type, protein, sort = TRUE)

write_csv(
  mr_cancer_supp,
  file.path(MANUSCRIPT_DIR, "Supplement_MR_cancer_hPDI_uPDI.csv")
)

## cvd
## ----
mr_cancer <- df_mr_sig %>%
  filter(
    if_any(
      any_of("disease"),
      ~ str_detect(str_to_lower(as.character(.x)), "cvd")
    )
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_match,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  )

mr_cancer_supp <- mr_cancer %>%
  dplyr::left_join(
    protein_lookup,
    by = c("pdi_type", "protein")
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_id,
    target,
    uniprot,
    rsid,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  ) %>%
  dplyr::arrange(pdi_type, protein, method)

mr_cancer_supp %>%
  dplyr::filter(is.na(protein_id) | is.na(target) | is.na(uniprot) | is.na(rsid)) %>%
  dplyr::count(pdi_type, protein, sort = TRUE)

head(mr_cancer_supp)
write_csv(
  mr_cancer_supp,
  file.path(MANUSCRIPT_DIR, "Supplement_MR_CVD_hPDI_uPDI.csv")
)

## t2d
mr_cancer <- df_mr_sig %>%
  filter(
    if_any(
      any_of("disease"),
      ~ str_detect(str_to_lower(as.character(.x)), "t2d")
    )
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_match,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  )

mr_cancer_supp <- mr_cancer %>%
  dplyr::left_join(
    protein_lookup,
    by = c("pdi_type", "protein")
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_id,
    target,
    uniprot,
    rsid,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  ) %>%
  dplyr::arrange(pdi_type, protein, method)

mr_cancer_supp %>%
  dplyr::filter(is.na(protein_id) | is.na(target) | is.na(uniprot) | is.na(rsid)) %>%
  dplyr::count(pdi_type, protein, sort = TRUE)

head(mr_cancer_supp)
write_csv(
  mr_cancer_supp,
  file.path(MANUSCRIPT_DIR, "Supplement_MR_T2D_hPDI_uPDI.csv")
)

## breast
## ----
mr_cancer <- df_mr_sig %>%
  filter(
    if_any(
      any_of("disease"),
      ~ str_detect(str_to_lower(as.character(.x)), "breast")
    )
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_match,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  )

mr_cancer_supp <- mr_cancer %>%
  dplyr::left_join(
    protein_lookup,
    by = c("pdi_type", "protein")
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_id,
    target,
    uniprot,
    rsid,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  ) %>%
  dplyr::arrange(pdi_type, protein, method)

mr_cancer_supp %>%
  dplyr::filter(is.na(protein_id) | is.na(target) | is.na(uniprot) | is.na(rsid)) %>%
  dplyr::count(pdi_type, protein, sort = TRUE)

head(mr_cancer_supp)
write_csv(
  mr_cancer_supp,
  file.path(MANUSCRIPT_DIR, "Supplement_MR_breast_hPDI_uPDI.csv")
)

## t2d
mr_cancer <- df_mr_sig %>%
  filter(
    if_any(
      any_of("disease"),
      ~ str_detect(str_to_lower(as.character(.x)), "t2d")
    )
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_match,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  )

mr_cancer_supp <- mr_cancer %>%
  dplyr::left_join(
    protein_lookup,
    by = c("pdi_type", "protein")
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_id,
    target,
    uniprot,
    rsid,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  ) %>%
  dplyr::arrange(pdi_type, protein, method)

mr_cancer_supp %>%
  dplyr::filter(is.na(protein_id) | is.na(target) | is.na(uniprot) | is.na(rsid)) %>%
  dplyr::count(pdi_type, protein, sort = TRUE)

head(mr_cancer_supp)
write_csv(
  mr_cancer_supp,
  file.path(MANUSCRIPT_DIR, "Supplement_MR_T2D_hPDI_uPDI.csv")
)

## cvd
## ----
mr_cancer <- df_mr_sig %>%
  filter(
    if_any(
      any_of("disease"),
      ~ str_detect(str_to_lower(as.character(.x)), "cvd")
    )
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_match,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  )

mr_cancer_supp <- mr_cancer %>%
  dplyr::left_join(
    protein_lookup,
    by = c("pdi_type", "protein")
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_id,
    target,
    uniprot,
    rsid,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  ) %>%
  dplyr::arrange(pdi_type, protein, method)

mr_cancer_supp %>%
  dplyr::filter(is.na(protein_id) | is.na(target) | is.na(uniprot) | is.na(rsid)) %>%
  dplyr::count(pdi_type, protein, sort = TRUE)

head(mr_cancer_supp)
write_csv(
  mr_cancer_supp,
  file.path(MANUSCRIPT_DIR, "Supplement_MR_CVD_hPDI_uPDI.csv")
)

## prostate
mr_cancer <- df_mr_sig %>%
  filter(
    if_any(
      any_of("disease"),
      ~ str_detect(str_to_lower(as.character(.x)), "prostate")
    )
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_match,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  )

mr_cancer_supp <- mr_cancer %>%
  dplyr::left_join(
    protein_lookup,
    by = c("pdi_type", "protein")
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_id,
    target,
    uniprot,
    rsid,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  ) %>%
  dplyr::arrange(pdi_type, protein, method)

mr_cancer_supp %>%
  dplyr::filter(is.na(protein_id) | is.na(target) | is.na(uniprot) | is.na(rsid)) %>%
  dplyr::count(pdi_type, protein, sort = TRUE)

head(mr_cancer_supp)

write_csv(
  mr_cancer_supp,
    file.path(MANUSCRIPT_DIR, "Supplement_MR_prostate_hPDI_uPDI.csv")
)

## lung
mr_cancer <- df_mr_sig %>%
  filter(
    if_any(
      any_of("disease"),
      ~ str_detect(str_to_lower(as.character(.x)), "lung")
    )
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_match,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  )

mr_cancer_supp <- mr_cancer %>%
  dplyr::left_join(
    protein_lookup,
    by = c("pdi_type", "protein")
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_id,
    target,
    uniprot,
    rsid,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  ) %>%
  dplyr::arrange(pdi_type, protein, method)

mr_cancer_supp %>%
  dplyr::filter(is.na(protein_id) | is.na(target) | is.na(uniprot) | is.na(rsid)) %>%
  dplyr::count(pdi_type, protein, sort = TRUE)

head(mr_cancer_supp)

write_csv(
  mr_cancer_supp,
  file.path(MANUSCRIPT_DIR, "Supplement_MR_breast_hPDI_uPDI.csv")
)

## lung
mr_cancer <- df_mr_sig %>%
  filter(
    if_any(
      any_of("disease"),
      ~ str_detect(str_to_lower(as.character(.x)), "lung")
    )
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_match,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  )

mr_cancer_supp <- mr_cancer %>%
  dplyr::left_join(
    protein_lookup,
    by = c("pdi_type", "protein")
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_id,
    target,
    uniprot,
    rsid,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  ) %>%
  dplyr::arrange(pdi_type, protein, method)

mr_cancer_supp %>%
  dplyr::filter(is.na(protein_id) | is.na(target) | is.na(uniprot) | is.na(rsid)) %>%
  dplyr::count(pdi_type, protein, sort = TRUE)

head(mr_cancer_supp)

write_csv(
  mr_cancer_supp,
  file.path(MANUSCRIPT_DIR, "Supplement_MR_lung_hPDI_uPDI.csv")
)
## colorectal
mr_cancer <- df_mr_sig %>%
  filter(
    if_any(
      any_of("disease"),
      ~ str_detect(str_to_lower(as.character(.x)), "colorectal")
    )
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_match,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  )

mr_cancer_supp <- mr_cancer %>%
  dplyr::left_join(
    protein_lookup,
    by = c("pdi_type", "protein")
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_id,
    target,
    uniprot,
    rsid,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  ) %>%
  dplyr::arrange(pdi_type, protein, method)

mr_cancer_supp %>%
  dplyr::filter(is.na(protein_id) | is.na(target) | is.na(uniprot) | is.na(rsid)) %>%
  dplyr::count(pdi_type, protein, sort = TRUE)

head(mr_cancer_supp)

write_csv(
  mr_cancer_supp,
  file.path(MANUSCRIPT_DIR, "Supplement_MR_colorectal_hPDI_uPDI.csv")
)

## lung
########
mr_cancer <- df_mr_sig %>%
  filter(
    if_any(
      any_of("disease"),
      ~ str_detect(str_to_lower(as.character(.x)), "kidney")
    )
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_match,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  )

mr_cancer_supp <- mr_cancer %>%
  dplyr::left_join(
    protein_lookup,
    by = c("pdi_type", "protein")
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_id,
    target,
    uniprot,
    rsid,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  ) %>%
  dplyr::arrange(pdi_type, protein, method)

mr_cancer_supp %>%
  dplyr::filter(is.na(protein_id) | is.na(target) | is.na(uniprot) | is.na(rsid)) %>%
  dplyr::count(pdi_type, protein, sort = TRUE)

head(mr_cancer_supp)

write_csv(
  mr_cancer_supp,
  file.path(MANUSCRIPT_DIR, "Supplement_MR_kidney_hPDI_uPDI.csv")
)

## lung
mr_cancer <- df_mr_sig %>%
  filter(
    if_any(
      any_of("disease"),
      ~ str_detect(str_to_lower(as.character(.x)), "liver")
    )
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_match,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  )

mr_cancer_supp <- mr_cancer %>%
  dplyr::left_join(
    protein_lookup,
    by = c("pdi_type", "protein")
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_id,
    target,
    uniprot,
    rsid,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  ) %>%
  dplyr::arrange(pdi_type, protein, method)

mr_cancer_supp %>%
  dplyr::filter(is.na(protein_id) | is.na(target) | is.na(uniprot) | is.na(rsid)) %>%
  dplyr::count(pdi_type, protein, sort = TRUE)

head(mr_cancer_supp)

write_csv(
  mr_cancer_supp,
  file.path(MANUSCRIPT_DIR, "Supplement_MR_liver_hPDI_uPDI.csv")
)

## endometrial
mr_cancer <- df_mr_sig %>%
  filter(
    if_any(
      any_of("disease"),
      ~ str_detect(str_to_lower(as.character(.x)), "endometrial")
    )
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_match,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  )

mr_cancer_supp <- mr_cancer %>%
  dplyr::left_join(
    protein_lookup,
    by = c("pdi_type", "protein")
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_id,
    target,
    uniprot,
    rsid,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  ) %>%
  dplyr::arrange(pdi_type, protein, method)

mr_cancer_supp %>%
  dplyr::filter(is.na(protein_id) | is.na(target) | is.na(uniprot) | is.na(rsid)) %>%
  dplyr::count(pdi_type, protein, sort = TRUE)

head(mr_cancer_supp)

write_csv(
  mr_cancer_supp,
  file.path(MANUSCRIPT_DIR, "Supplement_MR_ENDOM_hPDI_uPDI.csv")
)

## endometrial
mr_cancer <- df_mr_sig %>%
  filter(
    if_any(
      any_of("disease"),
      ~ str_detect(str_to_lower(as.character(.x)), "myeloma")
    )
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_match,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  )

mr_cancer_supp <- mr_cancer %>%
  dplyr::left_join(
    protein_lookup,
    by = c("pdi_type", "protein")
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_id,
    target,
    uniprot,
    rsid,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  ) %>%
  dplyr::arrange(pdi_type, protein, method)

mr_cancer_supp %>%
  dplyr::filter(is.na(protein_id) | is.na(target) | is.na(uniprot) | is.na(rsid)) %>%
  dplyr::count(pdi_type, protein, sort = TRUE)

head(mr_cancer_supp)

write_csv(
  mr_cancer_supp,
  file.path(MANUSCRIPT_DIR, "Supplement_MR_Myeloma_hPDI_uPDI.csv")
)

## Stomach
mr_cancer <- df_mr_sig %>%
  filter(
    if_any(
      any_of("disease"),
      ~ str_detect(str_to_lower(as.character(.x)), "gastric")
    )
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_match,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  )

mr_cancer_supp <- mr_cancer %>%
  dplyr::left_join(
    protein_lookup,
    by = c("pdi_type", "protein")
  ) %>%
  dplyr::select(
    pdi_type,
    disease,
    outcome_label,
    outcome_name,
    protein,
    protein_id,
    target,
    uniprot,
    rsid,
    gbm,
    limma,
    enet,
    exposure,
    method,
    nsnp,
    b,
    se,
    pval,
    p_fdr,
    or,
    or_lci95,
    or_uci95
  ) %>%
  dplyr::arrange(pdi_type, protein, method)

mr_cancer_supp %>%
  dplyr::filter(is.na(protein_id) | is.na(target) | is.na(uniprot) | is.na(rsid)) %>%
  dplyr::count(pdi_type, protein, sort = TRUE)

head(mr_cancer_supp)

write_csv(
  mr_cancer_supp,
  file.path(MANUSCRIPT_DIR, "Supplement_MR_gastric_hPDI_uPDI.csv")
)
