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

library(dplyr)

hpdi_gbm <- read.csv(file.path(MANUSCRIPT_DIR, "Hpdi_GBM_proteinlist.csv"))

hpdi_gbm_ranked <- hpdi_gbm %>%
  arrange(desc(gbm_mean_abs_test)) %>%
  mutate(
    shap_weight_test = gbm_mean_abs_test / sum(gbm_mean_abs_test, na.rm = TRUE),
    shap_percent_test = shap_weight_test * 100,
    shap_cumweight_test = cumsum(shap_weight_test),
    shap_cumpercent_test = shap_cumweight_test * 100,
    shap_rank_test = row_number(),
    shap_percentile_test = 100 * (1 - (shap_rank_test - 1) / n())
  )

if (interactive()) View(hpdi_gbm_ranked)

#merge with gain frequency cover

moredata <- read.csv(
  file.path(OUTPUT_ROOT, file.path("hpdi", "gbm", HPDI_RUN, "protein_signature.csv"))
)

# Merge Gain / Cover / Frequency
hpdi_gbm_ranked2 <- hpdi_gbm_ranked %>%
  left_join(
    moredata %>%
      select(
        protein_id = Feature,
        Gain,
        Cover,
        Frequency
      ),
    by = "protein_id"
  )

hpdi_gbm_supp <- hpdi_gbm_ranked2 %>%
  select(
    gene,
    protein_id,
    target,
    uniprot,
    
    shap_rank_test,
    gbm_mean_abs_test,
    shap_percent_test,
    shap_cumpercent_test,
    shap_percentile_test,
    
    gbm_direction,
    gbm_mean_abs_train,
    
    Gain,
    Cover,
    Frequency,
    
    everything()
  ) %>%
  rename(
    Gene = gene,
    protein_id = protein_id,
    Protein_target = target,
    UniProt = uniprot,
    
    SHAP_rank_test = shap_rank_test,
    Mean_abs_SHAP_test = gbm_mean_abs_test,
    SHAP_percent_test = shap_percent_test,
    Cumulative_SHAP_percent_test = shap_cumpercent_test,
    SHAP_percentile_test = shap_percentile_test,
    
    SHAP_direction_train = gbm_direction,
    Mean_abs_SHAP_train = gbm_mean_abs_train,
    
    XGBoost_Gain = Gain,
    XGBoost_Cover = Cover,
    XGBoost_Frequency = Frequency
  )


write.csv(
  hpdi_gbm_supp,
  file.path(MANUSCRIPT_DIR, "Hpdi_GBM_proteinlist_SHAP_test_ranked.csv"),
  row.names = FALSE
)

updi_gbm <- read.csv(file.path(MANUSCRIPT_DIR, "Updi_GBM_proteinlist.csv"))

updi_gbm_ranked <- updi_gbm %>%
  arrange(desc(gbm_mean_abs_test)) %>%
  mutate(
    shap_weight_test = gbm_mean_abs_test / sum(gbm_mean_abs_test, na.rm = TRUE),
    shap_percent_test = shap_weight_test * 100,
    shap_cumweight_test = cumsum(shap_weight_test),
    shap_cumpercent_test = shap_cumweight_test * 100,
    shap_rank_test = row_number(),
    shap_percentile_test = 100 * (1 - (shap_rank_test - 1) / n())
  )

if (interactive()) View(updi_gbm_ranked)

#merge with gain frequency cover

moredataupdi <- read.csv(
  file.path(OUTPUT_ROOT, file.path("updi", "gbm", UPDI_RUN, "protein_signature.csv"))
)

# Merge Gain / Cover / Frequency
updi_gbm_ranked2 <- updi_gbm_ranked %>%
  left_join(
    moredataupdi %>%
      select(
        protein_id = Feature,
        Gain,
        Cover,
        Frequency
      ),
    by = "protein_id"
  )

updi_gbm_supp <- updi_gbm_ranked2 %>%
  select(
    gene,
    protein_id,
    target,
    uniprot,
    
    shap_rank_test,
    gbm_mean_abs_test,
    shap_percent_test,
    shap_cumpercent_test,
    shap_percentile_test,
    
    gbm_direction,
    gbm_mean_abs_train,
    
    Gain,
    Cover,
    Frequency,
    
    everything()
  ) %>%
  rename(
    Gene = gene,
    protein_id = protein_id,
    Protein_target = target,
    UniProt = uniprot,
    
    SHAP_rank_test = shap_rank_test,
    Mean_abs_SHAP_test = gbm_mean_abs_test,
    SHAP_percent_test = shap_percent_test,
    Cumulative_SHAP_percent_test = shap_cumpercent_test,
    SHAP_percentile_test = shap_percentile_test,
    
    SHAP_direction_train = gbm_direction,
    Mean_abs_SHAP_train = gbm_mean_abs_train,
    
    XGBoost_Gain = Gain,
    XGBoost_Cover = Cover,
    XGBoost_Frequency = Frequency
  )


write.csv(
  updi_gbm_supp,
  file.path(MANUSCRIPT_DIR, "Updi_GBM_proteinlist_SHAP_test_ranked.csv"),
  row.names = FALSE
)
