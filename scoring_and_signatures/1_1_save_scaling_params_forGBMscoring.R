# Run from the repository root in a fresh R session. See README.md.
# REVIEW: reconstructed factor reference levels must match the fitted GBM; see README.
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

#'---- Syntax for ProPlant Project'#
#' Filename..: save_scaling_params_forGBMscoring.R
#' Purpose...: One-off script to save scaling_params.rds and
#'             factor_levels.rds into existing run folders for
#'             both hPDI and uPDI models.
#'             Run once interactively — does NOT retrain anything.
#' Writes into existing run folders:
#'   {hpdi_run_dir}/scaling_params.rds
#'   {hpdi_run_dir}/factor_levels.rds
#'   {updi_run_dir}/scaling_params.rds
#'   {updi_run_dir}/factor_levels.rds
#' --------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
})

# =================== USER CONFIG =====================================
OUT_BASE_HPDI <- file.path(OUTPUT_ROOT, "hpdi/gbm")
OUT_BASE_UPDI <- file.path(OUTPUT_ROOT, "updi/gbm")
# =====================================================================

# ------------------- Locate existing run folders --------------------
hpdi_run_dir <- file.path(OUTPUT_ROOT, file.path("hpdi", "gbm", HPDI_RUN))
updi_run_dir <- file.path(OUTPUT_ROOT, file.path("updi", "gbm", UPDI_RUN))

# Verify final models exist — fail early if not
stopifnot(file.exists(file.path(hpdi_run_dir, "final_model.rds")))
stopifnot(file.exists(file.path(updi_run_dir, "final_model.rds")))
cat("Both final_model.rds confirmed present.\n\n")

# ------------------- Load raw data once -----------------------------
load(file.path(DATA_DIR, "subcohort_20251106.RData"))
df_raw <- subcohort
cat(sprintf("Loaded subcohort: %d rows\n", nrow(df_raw)))

# ------------------- Create edu2 ------------------------------------
df_raw$edu2 <- dplyr::case_when(
  df_raw$L_School %in% c(0, 1)               ~ "Below secondary",
  df_raw$L_School %in% c(2, 3, 4)            ~ "Secondary/vocational and above",
  df_raw$L_School == 5 | is.na(df_raw$L_School) ~ "Not specified",
  TRUE                                        ~ NA_character_
)
df_raw$edu2 <- factor(
  df_raw$edu2,
  levels = c("Below secondary",
             "Secondary/vocational and above",
             "Not specified")
)

# ------------------- Helper function --------------------------------
save_params <- function(df_raw, outcome_col, run_dir) {
  
  cat(sprintf("\n--- Processing: %s ---\n", outcome_col))
  
  # Complete cases — must match training exactly
  model_covars <- c("Sex", "Age_Recr", outcome_col, "edu2",
                    "alc_freq", "smoke", "PA_total_MET", "Center")
  
  keep_rows <- stats::complete.cases(df_raw[, model_covars])
  df        <- df_raw[keep_rows, , drop = FALSE]
  cat(sprintf("Rows after complete case filter: %d\n", nrow(df)))
  
  # Save scaling params BEFORE any scaling
  scaling_params <- list(
    outcome = list(
      col  = outcome_col,
      mean = mean(df[[outcome_col]], na.rm = TRUE),
      sd   = sd(df[[outcome_col]],   na.rm = TRUE)
    ),
    Age_Recr = list(
      mean = mean(df$Age_Recr,     na.rm = TRUE),
      sd   = sd(df$Age_Recr,       na.rm = TRUE)
    ),
    PA_total_MET = list(
      mean = mean(df$PA_total_MET, na.rm = TRUE),
      sd   = sd(df$PA_total_MET,   na.rm = TRUE)
    )
  )
  
  cat(sprintf("Scaling params:\n"))
  cat(sprintf("  %-20s mean=%.4f, sd=%.4f\n", outcome_col,
              scaling_params$outcome$mean, scaling_params$outcome$sd))
  cat(sprintf("  %-20s mean=%.4f, sd=%.4f\n", "Age_Recr",
              scaling_params$Age_Recr$mean, scaling_params$Age_Recr$sd))
  cat(sprintf("  %-20s mean=%.4f, sd=%.4f\n", "PA_total_MET",
              scaling_params$PA_total_MET$mean, scaling_params$PA_total_MET$sd))
  
  # Save factor levels AFTER droplevels
  df$Center   <- droplevels(as.factor(df$Center))
  df$edu2     <- droplevels(df$edu2)
  df$alc_freq <- droplevels(as.factor(df$alc_freq))
  df$smoke    <- droplevels(as.factor(df$smoke))
  
  factor_levels <- list(
    Center   = levels(df$Center),
    edu2     = levels(df$edu2),
    alc_freq = levels(df$alc_freq),
    smoke    = levels(df$smoke),
    Sex      = levels(as.factor(df$Sex))
  )
  
  cat("Factor levels:\n")
  for (nm in names(factor_levels)) {
    cat(sprintf("  %-10s: %s\n", nm,
                paste(factor_levels[[nm]], collapse = ", ")))
  }
  
  # Save into existing run folder
  saveRDS(scaling_params, file.path(run_dir, "scaling_params.rds"))
  saveRDS(factor_levels,  file.path(run_dir, "factor_levels.rds"))
  
  cat(sprintf("Saved to: %s\n", run_dir))
}

# ------------------- Run for both indices ---------------------------
save_params(df_raw, outcome_col = "Hpdi_NoVegoils", run_dir = hpdi_run_dir)
save_params(df_raw, outcome_col = "Updi_NoVegoils", run_dir = updi_run_dir)

cat("\n========== Complete ==========\n")
cat("Files saved:\n")
cat(sprintf("  %s/scaling_params.rds\n", hpdi_run_dir))
cat(sprintf("  %s/factor_levels.rds\n",  hpdi_run_dir))
cat(sprintf("  %s/scaling_params.rds\n", updi_run_dir))
cat(sprintf("  %s/factor_levels.rds\n",  updi_run_dir))
cat("\nYou can now run 06_score_new_data_hPDI.R and 06_score_new_data_uPDI.R\n")
