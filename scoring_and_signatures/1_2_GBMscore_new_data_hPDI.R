# Run from the repository root in a fresh R session. See README.md.
# REVIEW: use the same run, factor coding and preprocessing as training; see README.
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
#' Filename..: 06_score_new_data_hPDI.R
#' Purpose...: Score an external dataset using the trained hPDI XGBoost model.
#' Corrections from interactive session:
#'   - final_model$feature_names returns NULL in this xgboost version
#'     -> train_cols sourced from 00_data_prep.R via X_train
#'   - cbind of dgCMatrix collapsed to numeric without na.pass
#'     -> na.action = na.pass added to sparse.model.matrix
#'   - All matrices explicitly cast to dgCMatrix before cbind
#' Reads  : {out_dir}/final_model.rds
#'          {out_dir}/scaling_params.rds
#'          {out_dir}/factor_levels.rds
#'          {out_dir}/data_splits.rds
#' Writes : {out_dir}/scores_external.csv
#' --------------------------------------------------------------------

# =================== USER CONFIG =====================================
OUTCOME_COL   <- "Hpdi_NoVegoils"
OUT_BASE      <- file.path(OUTPUT_ROOT, "hpdi/gbm")
SOURCE_DIR    <- file.path(SIGNATURE_SCRIPTS, "gbm/hpdi")
EXTERNAL_DATA <- file.path(DATA_DIR, "analyticaldf_postexcl_20251020.RData")
EXTERNAL_OBJ  <- "df"   # name of dataframe object inside RData
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(xgboost)
  library(data.table)
})

# ------------------- Locate run folder ------------------------------
out_dir <- trimws(readLines(file.path(OUT_BASE, "current_run_dir.txt")))
cat(sprintf("Using run folder: %s\n", out_dir))

# ------------------- Load training artefacts ------------------------
final_model    <- readRDS(file.path(out_dir, "final_model.rds"))
scaling_params <- readRDS(file.path(out_dir, "scaling_params.rds"))
factor_levels  <- readRDS(file.path(out_dir, "factor_levels.rds"))

cat("Training artefacts loaded.\n")
cat(sprintf("Scaling — %s:    mean=%.4f, sd=%.4f\n", OUTCOME_COL,
            scaling_params$outcome$mean,      scaling_params$outcome$sd))
#Scaling — Hpdi_NoVegoils:    mean=52.2040, sd=5.4103
cat(sprintf("Scaling — Age_Recr:     mean=%.4f, sd=%.4f\n",
            scaling_params$Age_Recr$mean,     scaling_params$Age_Recr$sd))
#Scaling — Age_Recr:     mean=51.4349, sd=8.5881
cat(sprintf("Scaling — PA_total_MET: mean=%.4f, sd=%.4f\n",
            scaling_params$PA_total_MET$mean, scaling_params$PA_total_MET$sd))
#Scaling — PA_total_MET: mean=103.4285, sd=61.1102
# ------------------- Get training column names ----------------------
# final_model$feature_names returns NULL in this xgboost version
# Source 00_data_prep.R to rebuild X_train and extract column names
cat("\nRebuilding training column names via 00_data_prep.R...\n")
DATA_FILE <- file.path(DATA_DIR, "subcohort_20251106.RData")
LOAD_SPLITS <- TRUE
source(file.path(SOURCE_DIR, "00_data_prep.R"))
out_dir    <- trimws(readLines(file.path(OUT_BASE, "current_run_dir.txt")))  # reset after source()
train_cols <- colnames(X_train)
cat(sprintf("Training columns: %d\n", length(train_cols)))

# Clean up training objects — keep only what we need
rm(X_full, X_train, X_test, y_train, y_test, dtrain, dtest,
   Xcov_full, Xprot_full, df, splits)

# ------------------- Load external dataset --------------------------
cat(sprintf("\nLoading external dataset: %s\n", EXTERNAL_DATA))
load(EXTERNAL_DATA)
df_new <- get(EXTERNAL_OBJ)
cat(sprintf("External dataset: %d rows, %d cols\n", nrow(df_new), ncol(df_new)))

# ------------------- Education variable (edu2) ----------------------
df_new$edu2 <- dplyr::case_when(
  df_new$L_School %in% c(0, 1)               ~ "Below secondary",
  df_new$L_School %in% c(2, 3, 4)            ~ "Secondary/vocational and above",
  df_new$L_School == 5 | is.na(df_new$L_School) ~ "Not specified",
  TRUE                                        ~ NA_character_
)

cat("\nedu2 distribution in external dataset:\n")
print(table(df_new$edu2, useNA = "ifany"))

# ------------------- Align factor levels to training ----------------
df_new$Center   <- factor(df_new$Center,   levels = factor_levels$Center)
df_new$edu2     <- factor(df_new$edu2,     levels = factor_levels$edu2)
df_new$alc_freq <- factor(df_new$alc_freq, levels = factor_levels$alc_freq)
df_new$smoke    <- factor(df_new$smoke,    levels = factor_levels$smoke)
df_new$Sex      <- factor(df_new$Sex,      levels = factor_levels$Sex)

cat("\nFactor level alignment checks:\n")
cat(sprintf("  Center NAs   : %d\n", sum(is.na(df_new$Center))))
cat(sprintf("  edu2 NAs     : %d\n", sum(is.na(df_new$edu2))))
cat(sprintf("  alc_freq NAs : %d\n", sum(is.na(df_new$alc_freq))))
cat(sprintf("  smoke NAs    : %d\n", sum(is.na(df_new$smoke))))

# ------------------- Scale continuous variables ---------------------
# Always use training mean/SD — not the external dataset's own
df_new$Age_Recr     <- (df_new$Age_Recr     - scaling_params$Age_Recr$mean)     / scaling_params$Age_Recr$sd
df_new$PA_total_MET <- (df_new$PA_total_MET - scaling_params$PA_total_MET$mean) / scaling_params$PA_total_MET$sd

cat("\nScaled covariates (training mean/SD applied):\n")
cat(sprintf("  Age_Recr     mean=%.4f\n", mean(df_new$Age_Recr,     na.rm = TRUE)))
cat(sprintf("  PA_total_MET mean=%.4f\n", mean(df_new$PA_total_MET, na.rm = TRUE)))

# Scale outcome if available — for correlation check only
has_outcome <- OUTCOME_COL %in% names(df_new)
if (has_outcome) {
  df_new[[paste0(OUTCOME_COL, "_raw")]] <- df_new[[OUTCOME_COL]]
  df_new[[OUTCOME_COL]] <- (df_new[[OUTCOME_COL]] - scaling_params$outcome$mean) /
    scaling_params$outcome$sd
  cat(sprintf("  %s mean=%.4f\n", OUTCOME_COL, mean(df_new[[OUTCOME_COL]], na.rm = TRUE)))
}

# ------------------- Protein columns --------------------------------
protein_cols <- grep("^seq", names(df_new), value = TRUE)
is_num       <- vapply(df_new[protein_cols], is.numeric, logical(1))
protein_cols <- protein_cols[is_num]
cat(sprintf("\nProtein columns in external dataset: %d\n", length(protein_cols)))

# ------------------- Build design matrix ----------------------------
# na.action = na.pass prevents sparse.model.matrix dropping rows with
# NA factor levels which causes cbind to silently collapse to numeric
cov_formula <- ~ Age_Recr + Sex + edu2 + alc_freq + smoke + PA_total_MET + Center

Xcov_new <- sparse.model.matrix(cov_formula, data = df_new, na.action = na.pass)
Xcov_new <- Xcov_new[, colnames(Xcov_new) != "(Intercept)", drop = FALSE]
Xcov_new <- methods::as(Xcov_new, "dgCMatrix")

Xprot_new <- methods::as(
  Matrix::Matrix(as.matrix(df_new[, protein_cols, drop = FALSE]), sparse = TRUE),
  "dgCMatrix")
Xprot_new[is.na(Xprot_new)] <- 0

X_new <- cbind(Xcov_new, Xprot_new)
cat(sprintf("X_new built: %s | %d x %d\n", class(X_new), nrow(X_new), ncol(X_new)))

# Verify matrix built correctly before proceeding
stopifnot(is.matrix(X_new) || inherits(X_new, "sparseMatrix"))
stopifnot(nrow(X_new) == nrow(df_new))

# ------------------- Column alignment -------------------------------
missing_cols <- setdiff(train_cols, colnames(X_new))
extra_cols   <- setdiff(colnames(X_new), train_cols)

cat(sprintf("\nColumn alignment:\n"))
cat(sprintf("  Training columns   : %d\n", length(train_cols)))
cat(sprintf("  External columns   : %d\n", ncol(X_new)))
cat(sprintf("  Missing in external: %d\n", length(missing_cols)))
cat(sprintf("  Extra in external  : %d\n", length(extra_cols)))

if (length(missing_cols) > 0) {
  cat("  Adding missing columns as zeros...\n")
  zero_mat <- methods::as(
    Matrix::Matrix(0,
                   nrow = nrow(X_new),
                   ncol = length(missing_cols),
                   dimnames = list(NULL, missing_cols),
                   sparse = TRUE),
    "dgCMatrix")
  X_new <- cbind(X_new, zero_mat)
}

if (length(extra_cols) > 0) {
  X_new <- X_new[, !colnames(X_new) %in% extra_cols, drop = FALSE]
}

X_new <- X_new[, train_cols, drop = FALSE]
stopifnot(identical(colnames(X_new), train_cols))
cat("Column alignment complete.\n")

# ------------------- Score ------------------------------------------
cat("\nScoring external dataset...\n")
dnew <- xgb.DMatrix(X_new)
df_new$GBMfull_Hpdi <- predict(final_model, dnew)

cat("Score distribution:\n")
print(summary(df_new$GBMfull_Hpdi))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.76728 -0.12281  0.02329  0.02915  0.16641  1.25134 

# ------------------- Protein SHAP contribution score ----------------
pred_full <- predict(final_model, dnew)

contrib <- predict(final_model, dnew, predcontrib = TRUE)

# Add column names if missing
if (is.null(colnames(contrib))) {
  colnames(contrib) <- c(colnames(X_new), "BIAS")
}

# Check that SHAP contributions sum to the full prediction
check <- all.equal(
  as.numeric(rowSums(contrib)),
  as.numeric(pred_full),
  tolerance = 1e-6
)

cat("\nSHAP decomposition check:\n")
print(check)

protein_shap_cols <- grep("^seq", colnames(contrib), value = TRUE)

cat(sprintf("Protein SHAP columns: %d\n", length(protein_shap_cols)))

df_new$GBMprotSHAP_Hpdi <- rowSums(
  contrib[, protein_shap_cols, drop = FALSE]
)

cat("Protein SHAP score distribution:\n")
print(summary(df_new$GBMprotSHAP_Hpdi))
# ------------------- Validation -------------------------------------
if (has_outcome) {
  rho <- cor(df_new$GBMfull_Hpdi, df_new[[OUTCOME_COL]],
             method = "spearman", use = "complete.obs")
  r   <- cor(df_new$GBMfull_Hpdi, df_new[[OUTCOME_COL]],
             method = "pearson",  use = "complete.obs")
  cat(sprintf("\nCorrelation with observed %s:\n", OUTCOME_COL))
  cat(sprintf("  Spearman rho : %.4f\n", rho))
  cat(sprintf("  Pearson r    : %.4f\n", r))
}
# Correlation with observed Hpdi_NoVegoils:
# Spearman rho : 0.3390
# Pearson r    : 0.3655

if (has_outcome) {
  rho <- cor(df_new$GBMprotSHAP_Hpdi, df_new[[OUTCOME_COL]],
             method = "spearman", use = "complete.obs")
  r   <- cor(df_new$GBMprotSHAP_Hpdi, df_new[[OUTCOME_COL]],
             method = "pearson",  use = "complete.obs")
  cat(sprintf("\nCorrelation with observed %s:\n", OUTCOME_COL))
  cat(sprintf("  Spearman rho : %.4f\n", rho))
  cat(sprintf("  Pearson r    : %.4f\n", r))
}
# Correlation with observed Hpdi_NoVegoils:
#   Spearman rho : 0.3124
# Pearson r    : 0.3416


# ------------------- Save scores ------------------------------------
save_cols  <- intersect(c("Idepic", "GBMfull_Hpdi", "GBMprotSHAP_Hpdi",
                          paste0(OUTCOME_COL, "_raw")), names(df_new))
scores_out <- df_new[, save_cols, drop = FALSE]

data.table::fwrite(scores_out, file.path(out_dir, "scores_external.csv"))
cat(sprintf("\nScores saved: scores_external.csv (%d individuals)\n", nrow(scores_out)))
cat(sprintf("Run folder: %s\n", out_dir))
cat(sprintf("Completed: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
