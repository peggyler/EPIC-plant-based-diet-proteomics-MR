# Run from the repository root in a fresh R session. See README.md.
# EDIT LOCAL PATHS HERE. No study data or results are supplied.
PROJECT_ROOT <- "."
DATA_DIR <- file.path(PROJECT_ROOT, "data")
OUTPUT_ROOT <- file.path(PROJECT_ROOT, "outputs")
MANUSCRIPT_DIR <- file.path(OUTPUT_ROOT, "manuscript")
REFERENCE_DIR <- file.path(DATA_DIR, "reference")
SIGNATURE_SCRIPTS <- file.path(PROJECT_ROOT, "signatures")
CONTEXT_SCRIPTS <- file.path(PROJECT_ROOT, "contextualisation")
HPDI_RUN <- "staged_run_20260424_0935"  # choose the manuscript-producing run/date
UPDI_RUN <- "staged_run_20260424_1844"  # choose the manuscript-producing run/date
dir.create(MANUSCRIPT_DIR, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(dplyr))

proj_root    <- MANUSCRIPT_DIR
stable <- read.csv(file.path(OUTPUT_ROOT, "hpdi/stability_bootstrap_ENet_Hpdi_NoVegoils_filtFDR5_lamMIN.csv"))

enet <- read.csv(file.path(proj_root, "Hpdi_ENET_proteinlist.csv"))

stable2 <- stable %>%
  rename(protein_id = feature) 

enet1 <- enet %>%
  select(-any_of(setdiff(names(stable2), "protein_id")),
         -matches("^sel_freq[.][xy]$")) %>%
  left_join(stable2, by = "protein_id")

enet2 <- enet1 %>%
  select(gene, protein_id, target, uniprot, elasticnet_coef, elasticnet_abs_coef, sel_freq,everything()
  )

write.csv(enet2, file.path(proj_root, "Hpdi_ENET_proteinlist.csv"), row.names = FALSE)

##_________
## UPDI
##__________
stable <- read.csv(file.path(OUTPUT_ROOT, "updi/stability_bootstrap_ENet_Updi_NoVegoils_filtFDR5_lamMIN.csv"))

enet <- read.csv(file.path(proj_root, "Updi_ENET_proteinlist.csv"))

stable2 <- stable %>%
  rename(protein_id = feature) 

enet1 <- enet %>%
  select(-any_of(setdiff(names(stable2), "protein_id")),
         -matches("^sel_freq[.][xy]$")) %>%
  left_join(stable2, by = "protein_id")

enet2 <- enet1 %>%
  select(gene, protein_id, target, uniprot, elasticnet_coef, elasticnet_abs_coef, sel_freq,everything()
  )

write.csv(enet2, file.path(proj_root, "Updi_ENET_proteinlist.csv"), row.names = FALSE)

## x gradient boost signatures
## ____________________________
proj_root    <- MANUSCRIPT_DIR
gbm <- read.csv(file.path(OUTPUT_ROOT, file.path("hpdi", "gbm", HPDI_RUN, "signature_SHAP.csv")))
gbm1 <- read.csv(file.path(OUTPUT_ROOT, file.path("hpdi", "gbm", HPDI_RUN, "feature_importance.csv")))
gbm3 <- read.csv(file.path(proj_root, "Hpdi_GBM_proteinlist.csv"))

## rename XGBweights to relative contribution _ in GBM3
supp_xgb <-gbm3 %>%
  rename(
    xgboost_signature                     = gbm,
    relative_mean_abs_shap_train          = gbm_weight_abs,
    signed_relative_shap_train            = gbm_weight_signed,
    mean_signed_shap_train                = gbm_mean_shap_train,
    mean_abs_shap_train                   = gbm_mean_abs_train,
    mean_abs_shap_test                    = gbm_mean_abs_test,
    average_shap_direction                = gbm_direction,
    cumulative_relative_shap_train        = gbm_cumweight
  )

xgb_imp_proteins <- gbm1 %>%
  filter(grepl("^seq\\.", Feature)) %>%
  rename(
    protein_id        = Feature
  )

# ---------------------------------------------------------------------
# Join SHAP summary with XGBoost gain/cover/frequency
# ---------------------------------------------------------------------
supp_xgb <- supp_xgb %>%
  left_join(xgb_imp_proteins, by = "protein_id") %>%
  mutate(
    relative_mean_abs_shap_test =
      mean_abs_shap_test / sum(mean_abs_shap_test, na.rm = TRUE)
  ) %>%
  arrange(desc(mean_abs_shap_test)) %>%
  mutate(
    cumulative_relative_shap_test =
      cumsum(relative_mean_abs_shap_test)
  ) %>%
  select(
    gene,
    protein_id,
    target,
    uniprot,
    seq_id,
    seq_id_decode,
    xgboost_signature,
    mean_abs_shap_test,
    relative_mean_abs_shap_test,
    cumulative_relative_shap_test,
    mean_abs_shap_train,
    relative_mean_abs_shap_train,
    signed_relative_shap_train,
    mean_signed_shap_train,
    average_shap_direction,
    Gain,
    Cover,
    Frequency,
    type,
    dilution,
    organism,
    everything()
  )

if (interactive()) View(supp_xgb)
write.csv(supp_xgb, file.path(proj_root, "Hpdi_XGBoost_proteinlist.csv"), row.names = FALSE)

##______
## UPDI
##______
proj_root    <- MANUSCRIPT_DIR
gbm <- read.csv(file.path(OUTPUT_ROOT, file.path("updi", "gbm", UPDI_RUN, "signature_SHAP.csv")))
gbm1 <- read.csv(file.path(OUTPUT_ROOT, file.path("updi", "gbm", UPDI_RUN, "feature_importance.csv")))
gbm3 <- read.csv(file.path(proj_root, "Updi_GBM_proteinlist.csv"))

## rename XGBweights to relative contribution _ in GBM3
supp_xgb <-gbm3 %>%
  rename(
    xgboost_signature                     = gbm,
    relative_mean_abs_shap_train          = gbm_weight_abs,
    signed_relative_shap_train            = gbm_weight_signed,
    mean_signed_shap_train                = gbm_mean_shap_train,
    mean_abs_shap_train                   = gbm_mean_abs_train,
    mean_abs_shap_test                    = gbm_mean_abs_test,
    average_shap_direction                = gbm_direction,
    cumulative_relative_shap_train        = gbm_cumweight
  )

xgb_imp_proteins <- gbm1 %>%
  filter(grepl("^seq\\.", Feature)) %>%
  rename(
    protein_id        = Feature
  )

# ---------------------------------------------------------------------
# Join SHAP summary with XGBoost gain/cover/frequency
# ---------------------------------------------------------------------
supp_xgb <- supp_xgb %>%
  left_join(xgb_imp_proteins, by = "protein_id") %>%
  mutate(
    relative_mean_abs_shap_test =
      mean_abs_shap_test / sum(mean_abs_shap_test, na.rm = TRUE)
  ) %>%
  arrange(desc(mean_abs_shap_test)) %>%
  mutate(
    cumulative_relative_shap_test =
      cumsum(relative_mean_abs_shap_test)
  ) %>%
  select(
    gene,
    protein_id,
    target,
    uniprot,
    seq_id,
    seq_id_decode,
    xgboost_signature,
    mean_abs_shap_test,
    relative_mean_abs_shap_test,
    cumulative_relative_shap_test,
    mean_abs_shap_train,
    relative_mean_abs_shap_train,
    signed_relative_shap_train,
    mean_signed_shap_train,
    average_shap_direction,
    Gain,
    Cover,
    Frequency,
    type,
    dilution,
    organism,
    everything()
  )

if (interactive()) View(supp_xgb)
write.csv(supp_xgb, file.path(proj_root, "Updi_XGBoost_proteinlist.csv"), row.names = FALSE)
