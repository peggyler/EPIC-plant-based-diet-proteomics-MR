# Run from the repository root in a fresh R session. See README.md.
# REVIEW: fallback scaling is calculated after filtering, unlike 1.1 training; see README.
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
#' Filename..: score_Enet_Updi.R
#' Purpose...: Score an external dataset using the trained Updi elastic net model.
#'             Uses predict() — handles intercept automatically.
#'             Aligns design matrix columns to training exactly via bundle$colnames_X
#' Key differences from GBM scoring:
#'   - Proteins go in RAW (standardize=FALSE in training)
#'   - PA_total_MET_quint used (factor quintile) not PA_total_MET continuous
#'   - smoke ref = "Never", alc_freq ref = "Never"
#'   - Outcome and Age_Recr scaled using training mean/SD
#' Reads  : elasticnet_model_*.rds   (model bundle)
#'          scaling_params_enet.rds  (training mean/SD for outcome + Age_Recr)
#' Writes : {outdir}/scores_external_enet.csv
#' --------------------------------------------------------------------

# =================== USER CONFIG =====================================
OUTCOME_COL   <- "Updi_NoVegoils"
OUTDIR        <- file.path(OUTPUT_ROOT, "updi")
EXTERNAL_DATA <- file.path(DATA_DIR, "analyticaldf_postexcl_20251020.RData")
EXTERNAL_OBJ  <- "df"
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(glmnet)
  library(data.table)
})

# ------------------- Load model bundle ------------------------------
bundle_path <- list.files(OUTDIR, pattern = "^elasticnet_model_.*\\.rds$", full.names = TRUE)
if (length(bundle_path) != 1) stop("Cannot uniquely find elasticnet_model_*.rds in OUTDIR")
bundle <- readRDS(bundle_path)

cat(sprintf("Model bundle loaded: %s\n", basename(bundle_path)))
cat(sprintf("  alpha      : %.2f\n",  bundle$best_alpha))
cat(sprintf("  lambda     : %.6g\n",  bundle$best_lambda))
cat(sprintf("  Proteins   : %d\n",    length(bundle$protein_cols)))
cat(sprintf("  X columns  : %d\n",    length(bundle$colnames_X)))

# ------------------- Load scaling params ----------------------------
# Saved by 00_save_scaling_params_oneoff.R
# Contains training mean/SD for outcome and Age_Recr
scaling_path <- file.path(OUTDIR, "scaling_params_enet.rds")

if (file.exists(scaling_path)) {
  scaling_params <- readRDS(scaling_path)
  cat(sprintf("Scaling params loaded from: %s\n", scaling_path))
} else {
  # Derive from training data if not saved
  cat("scaling_params_enet.rds not found — deriving from subcohort...\n")
  load(file.path(DATA_DIR, "subcohort_20251106.RData"))
  df_sub <- subcohort
  
  # Complete cases filter matching training
  df_sub$edu2 <- dplyr::case_when(
    df_sub$L_School %in% c(0, 1)               ~ "Below secondary",
    df_sub$L_School %in% c(2, 3, 4)            ~ "Secondary/vocational and above",
    df_sub$L_School == 5 | is.na(df_sub$L_School) ~ "Not specified",
    TRUE                                        ~ NA_character_
  )
  model_covars <- c("Sex", "Age_Recr", OUTCOME_COL, "edu2",
                    "alc_freq", "smoke", "PA_total_MET_quint", "Center")
  df_sub <- df_sub[stats::complete.cases(df_sub[, model_covars]), , drop = FALSE]
  
  scaling_params <- list(
    outcome  = list(col  = OUTCOME_COL,
                    mean = mean(df_sub[[OUTCOME_COL]], na.rm = TRUE),
                    sd   = sd(df_sub[[OUTCOME_COL]],   na.rm = TRUE)),
    Age_Recr = list(mean = mean(df_sub$Age_Recr,       na.rm = TRUE),
                    sd   = sd(df_sub$Age_Recr,          na.rm = TRUE))
  )
  saveRDS(scaling_params, scaling_path)
  cat(sprintf("Scaling params saved to: %s\n", scaling_path))
  rm(df_sub, subcohort)
}

cat(sprintf("Scaling — %s: mean=%.4f, sd=%.4f\n", OUTCOME_COL,
            scaling_params$outcome$mean, scaling_params$outcome$sd))
cat(sprintf("Scaling — Age_Recr: mean=%.4f, sd=%.4f\n",
            scaling_params$Age_Recr$mean, scaling_params$Age_Recr$sd))

# ------------------- Load external dataset --------------------------
cat(sprintf("\nLoading external dataset: %s\n", EXTERNAL_DATA))
load(EXTERNAL_DATA)
df_new <- get(EXTERNAL_OBJ)
cat(sprintf("External dataset: %d rows, %d cols\n", nrow(df_new), ncol(df_new)))

df_new$edu2 <- dplyr::case_when(
  df_new$L_School %in% c(0, 1)               ~ "Below secondary",
  df_new$L_School %in% c(2, 3, 4)            ~ "Secondary/vocational and above",
  df_new$L_School == 5 | is.na(df_new$L_School) ~ "Not specified",
  TRUE                                        ~ NA_character_
)

# ------------------- Check required columns -------------------------
required <- c("Age_Recr", "Sex", "edu2", "alc_freq", "smoke",
              "PA_total_MET_quint", "Center")
missing_vars <- setdiff(required, names(df_new))
if (length(missing_vars) > 0) {
  stop(sprintf("Missing required columns in external dataset: %s",
               paste(missing_vars, collapse = ", ")))
}

# ------------------- Education variable (edu2) ----------------------
df_new$edu2 <- dplyr::case_when(
  df_new$L_School %in% c(0, 1)               ~ "Below secondary",
  df_new$L_School %in% c(2, 3, 4)            ~ "Secondary/vocational and above",
  df_new$L_School == 5 | is.na(df_new$L_School) ~ "Not specified",
  TRUE                                        ~ NA_character_
)
df_new$edu2 <- factor(df_new$edu2,
                      levels = c("Below secondary",
                                 "Secondary/vocational and above",
                                 "Not specified"))

cat("\nedu2 distribution in external dataset:\n")
print(table(df_new$edu2, useNA = "ifany"))
# Below secondary Secondary/vocational and above                  Not specified 
# 7603                           6793                            486 
# ------------------- Factor setup matching training -----------------
# Reference levels must match training exactly
df_new$smoke    <- factor(df_new$smoke)
df_new$alc_freq <- factor(df_new$alc_freq)
df_new$smoke    <- relevel(df_new$smoke,    ref = "Never")
df_new$alc_freq <- relevel(df_new$alc_freq, ref = "Never")

df_new$Center             <- as.factor(df_new$Center)
df_new$Sex                <- as.factor(df_new$Sex)
df_new$PA_total_MET_quint <- as.factor(df_new$PA_total_MET_quint)
df_new$edu2               <- droplevels(df_new$edu2)

cat("\nFactor checks:\n")
cat(sprintf("  smoke levels    : %s\n", paste(levels(df_new$smoke),    collapse = ", ")))
cat(sprintf("  alc_freq levels : %s\n", paste(levels(df_new$alc_freq), collapse = ", ")))
cat(sprintf("  PA_quint levels : %s\n", paste(levels(df_new$PA_total_MET_quint), collapse = ", ")))

# ------------------- Scale continuous variables ---------------------
# Use training mean/SD — NOT external dataset's own
df_new$Age_Recr <- (df_new$Age_Recr - scaling_params$Age_Recr$mean) /
  scaling_params$Age_Recr$sd
cat(sprintf("\nAge_Recr scaled: mean=%.4f\n", mean(df_new$Age_Recr, na.rm = TRUE)))

# Scale outcome if available — for correlation check only
has_outcome <- OUTCOME_COL %in% names(df_new)
if (has_outcome) {
  df_new[[paste0(OUTCOME_COL, "_raw")]] <- df_new[[OUTCOME_COL]]
  df_new[[OUTCOME_COL]] <- (df_new[[OUTCOME_COL]] - scaling_params$outcome$mean) /
    scaling_params$outcome$sd
  cat(sprintf("%s scaled: mean=%.4f\n", OUTCOME_COL,
              mean(df_new[[OUTCOME_COL]], na.rm = TRUE)))
}

# ------------------- Build design matrix ----------------------------
# Must match training formula exactly
# Proteins go in RAW — standardize=FALSE was used in training
cov_formula <- ~ Age_Recr + Sex + edu2 + alc_freq + smoke +
  PA_total_MET_quint + Center

Xcov_new <- model.matrix(cov_formula, data = df_new)[, -1, drop = FALSE]

# Protein matrix — raw, no scaling
prot_cols_new <- intersect(bundle$protein_cols, colnames(df_new))
missing_prots <- setdiff(bundle$protein_cols, colnames(df_new))
cat(sprintf("\nProteins in bundle    : %d\n", length(bundle$protein_cols)))
cat(sprintf("Proteins in external  : %d\n", length(prot_cols_new)))
cat(sprintf("Proteins missing      : %d\n", length(missing_prots)))

Xprot_new <- as.matrix(df_new[, prot_cols_new, drop = FALSE])
storage.mode(Xprot_new) <- "double"

X_new <- cbind(Xcov_new, Xprot_new)
cat(sprintf("X_new built: %d x %d\n", nrow(X_new), ncol(X_new)))

# ------------------- Column alignment -------------------------------
# Align to training column order exactly
train_cols   <- bundle$colnames_X
missing_cols <- setdiff(train_cols, colnames(X_new))
extra_cols   <- setdiff(colnames(X_new), train_cols)

cat(sprintf("\nColumn alignment:\n"))
cat(sprintf("  Training columns   : %d\n", length(train_cols)))
cat(sprintf("  External columns   : %d\n", ncol(X_new)))
cat(sprintf("  Missing in external: %d\n", length(missing_cols)))
cat(sprintf("  Extra in external  : %d\n", length(extra_cols)))

# Add missing columns as zeros
if (length(missing_cols) > 0) {
  cat("  Missing columns (filled with 0):\n")
  cat(paste("   ", head(missing_cols, 10), collapse = "\n"), "\n")
  zero_mat <- matrix(0, nrow = nrow(X_new), ncol = length(missing_cols),
                     dimnames = list(NULL, missing_cols))
  X_new <- cbind(X_new, zero_mat)
}

# Drop extra columns
if (length(extra_cols) > 0) {
  X_new <- X_new[, !colnames(X_new) %in% extra_cols, drop = FALSE]
}

# Reorder to match training exactly
X_new <- X_new[, train_cols, drop = FALSE]
stopifnot(identical(colnames(X_new), train_cols))
cat("Column alignment complete.\n")

# ------------------- Score via predict() ----------------------------
# predict() handles intercept automatically
cat("\nScoring external dataset...\n")
df_new$Enet_Updi <- as.numeric(
  predict(bundle$fit_final, newx = X_new, s = bundle$best_lambda)
)

cat("Score distribution:\n")
print(summary(df_new$Enet_Updi))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.3982 -0.4589 -0.1963 -0.0447  0.3923  1.8090 

# ------------------- ENET decomposition -----------------------------

cat("\nDecomposing ENET prediction into intercept, covariate, and protein parts...\n")

coef_enet <- as.matrix(
  coef(bundle$fit_final, s = bundle$best_lambda)
)

coef_vec <- coef_enet[, 1]
names(coef_vec) <- rownames(coef_enet)

# Intercept
enet_intercept <- unname(coef_vec["(Intercept)"])

# Coefficients excluding intercept
beta <- coef_vec[setdiff(names(coef_vec), "(Intercept)")]

# Align beta to X_new columns
beta_aligned <- beta[colnames(X_new)]

# Safety check
if (any(is.na(beta_aligned))) {
  missing_beta <- colnames(X_new)[is.na(beta_aligned)]
  stop(sprintf(
    "Some X_new columns are missing coefficients: %s",
    paste(head(missing_beta, 20), collapse = ", ")
  ))
}

stopifnot(identical(names(beta_aligned), colnames(X_new)))

# Identify protein and covariate columns
protein_cols_model <- intersect(bundle$protein_cols, colnames(X_new))
cov_cols_model     <- setdiff(colnames(X_new), protein_cols_model)

cat(sprintf("Protein columns in model   : %d\n", length(protein_cols_model))) #16
cat(sprintf("Covariate columns in model : %d\n", length(cov_cols_model))) #29
cat(sprintf("ENET intercept             : %.6f\n", enet_intercept)) #0.9

# Protein-only contribution, no intercept
df_new$EnetProt_Updi <- as.numeric(
  X_new[, protein_cols_model, drop = FALSE] %*%
    beta_aligned[protein_cols_model]
)

# Protein contribution + intercept
df_new$EnetProtInt_Updi <- enet_intercept + df_new$EnetProt_Updi

# Covariate contribution
df_new$EnetCov_Updi <- as.numeric(
  X_new[, cov_cols_model, drop = FALSE] %*%
    beta_aligned[cov_cols_model]
)

# Manual full prediction
df_new$EnetFull_manual_Updi <-
  enet_intercept +
  df_new$EnetCov_Updi +
  df_new$EnetProt_Updi

# Check against predict()
cat("\nCheck manual full prediction vs predict():\n")
print(all.equal(
  df_new$EnetFull_manual_Updi,
  df_new$Enet_Updi,
  tolerance = 1e-8
))

cat("\nScore distributions:\n")
cat("\nFull ENET prediction, from predict():\n")
print(summary(df_new$Enet_Updi))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.3982 -0.4589 -0.1963 -0.0447  0.3923  1.8090 

cat("\nProtein ENET contribution, no intercept:\n")
print(summary(df_new$EnetProt_Updi))
#  Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-1.004185 -0.138600 -0.021642 -0.006474  0.109118  1.073956  

cat("\nProtein ENET score with intercept:\n")
print(summary(df_new$EnetProtInt_Updi))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.09655  0.76904  0.88599  0.90116  1.01675  1.98159 

cat("\nCovariate ENET contribution:\n")

print(summary(df_new$EnetCov_Updi))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -2.0221 -1.3312 -1.0700 -0.9459 -0.4589  0.5640 

# ------------------- Validation -------------------------------------

if (has_outcome) {
  
  score_vars <- c(
    "Enet_Updi",           # full prediction: intercept + covariates + proteins
    "EnetProt_Updi",       # protein contribution only
    "EnetProtInt_Updi",    # protein contribution + intercept
    "EnetCov_Updi"         # covariate contribution only
  )
  
  cat(sprintf("\nCorrelation with observed scaled %s:\n", OUTCOME_COL))
  
  for (score_var in score_vars) {
    
    rho <- cor(df_new[[score_var]], df_new[[OUTCOME_COL]],
               method = "spearman", use = "complete.obs")
    
    r <- cor(df_new[[score_var]], df_new[[OUTCOME_COL]],
             method = "pearson", use = "complete.obs")
    
    cat(sprintf("\n%s:\n", score_var))
    cat(sprintf("  Spearman rho : %.4f\n", rho))
    cat(sprintf("  Pearson r    : %.4f\n", r))
  }
}

# Enet_Updi:
#   Spearman rho : 0.5266
# Pearson r    : 0.5545
# 
# EnetProt_Updi:
#   Spearman rho : 0.1128
# Pearson r    : 0.1092
# 
# EnetProtInt_Updi:
#   Spearman rho : 0.1128
# Pearson r    : 0.1092
# 
# EnetCov_Updi:
#   Spearman rho : 0.4899
# Pearson r    : 0.5273

# ------------------- Save scores ------------------------------------
save_cols  <- intersect(c("Idepic", "Enet_Updi", "EnetProt_Updi", "EnetProtInt_Updi", "EnetCov_Updi",
                          paste0(OUTCOME_COL, "_raw")), names(df_new))
scores_out <- df_new[, save_cols, drop = FALSE]

data.table::fwrite(scores_out, file.path(OUTDIR, "scores_external_enet.csv"))
cat(sprintf("\nScores saved: scores_external_enet.csv (%d individuals)\n",
            nrow(scores_out)))
cat(sprintf("Completed: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
