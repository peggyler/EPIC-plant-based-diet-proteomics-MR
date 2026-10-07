#'---- Syntax for ProPlant Project'#
#' Created by Peggy Ler (Pler) - 31 Oct 2025
#' Filename..: 05_compute_shap_updi.R
#' Updated...: [date] by "Pler"
#' Purpose...: Compute SHAP values for final updi XGBoost model
#'             Uses saved train/test splits from 01_stage_A to ensure
#'             SHAP is computed on identical partitions used in training.
#'             Computes SHAP on both train and test sets.
#' Reads  : current_run_dir.txt -> final_model.rds, data_splits.rds
#' Writes : {out_dir}/shap_train.csv       <- full SHAP matrix (train)
#'          {out_dir}/shap_test.csv        <- full SHAP matrix (test)
#'          {out_dir}/signature_SHAP.csv   <- per-protein summary
#' --------------------------------------------------------------------

# USER PATHS: edit these locally; run from the folder containing README.md.
DATA_FILE <- "data/subcohort_20251106.RData"
OUTPUT_ROOT <- "outputs"

# =================== USER CONFIG =====================================
OUTCOME_COL <- "Updi_NoVegoils"
OUT_BASE <- file.path(OUTPUT_ROOT, "updi", "gbm")
if (!dir.exists(OUT_BASE)) dir.create(OUT_BASE, recursive = TRUE)
OUT_BASE <- normalizePath(OUT_BASE, winslash = "/", mustWork = TRUE)
SOURCE_DIR <- file.path("signatures", "gbm", "updi")
LOAD_SPLITS <- TRUE    # must be TRUE — load same splits used in training
# =====================================================================

# ------------------- Locate run folder ------------------------------
run_pointer <- file.path(OUT_BASE, "current_run_dir.txt")
if (!file.exists(run_pointer)) stop("Run Stage A (01) first.")
out_dir <- trimws(readLines(run_pointer, warn = FALSE))
if (length(out_dir) != 1L || !nzchar(out_dir) || !dir.exists(out_dir)) {
  stop("current_run_dir.txt must name one existing run directory.")
}
cat(sprintf("Using run folder: %s\n", out_dir))

# ------------------- Source shared prep -----------------------------
# Loads data, creates edu2, scales covariates, builds design matrices
# Loads saved train/test splits (LOAD_SPLITS = TRUE)
source(file.path(SOURCE_DIR, "00_data_prep_updi.R"), local = TRUE)

# ------------------- Load final model -------------------------------
model_path <- file.path(out_dir, "final_model.rds")
if (!file.exists(model_path)) stop("final_model.rds not found — run 04 first")
final_model <- readRDS(model_path)
cat("Final model loaded.\n")

# ------------------- Compute SHAP values ----------------------------
# predcontrib = TRUE returns SHAP contributions per feature per observation
# Last column is BIAS (intercept term) — removed below

cat("Computing SHAP values on training set...\n")
shap_train_raw <- predict(final_model, dtrain, predcontrib = TRUE)
shap_train     <- as.data.frame(shap_train_raw)
if ("BIAS" %in% names(shap_train)) shap_train$BIAS <- NULL

cat("Computing SHAP values on test set...\n")
shap_test_raw <- predict(final_model, dtest, predcontrib = TRUE)
shap_test     <- as.data.frame(shap_test_raw)
if ("BIAS" %in% names(shap_test)) shap_test$BIAS <- NULL

# ------------------- Proteins only ----------------------------------
prot_cols <- grep("^seq", colnames(shap_train), value = TRUE)

S_train <- shap_train[, prot_cols, drop = FALSE]
S_test  <- shap_test[,  prot_cols, drop = FALSE]

# Drop proteins with all-zero SHAP across training set
nz <- colSums(abs(S_train), na.rm = TRUE) > 0
S_train <- S_train[, nz, drop = FALSE]
S_test  <- S_test[,  nz, drop = FALSE]   # apply same filter to test

cat(sprintf("Proteins with non-zero SHAP (train): %d\n", ncol(S_train)))
cat(sprintf("Proteins retained in test set:        %d\n", ncol(S_test)))

# ------------------- Per-protein SHAP summary -----------------------
# Computed on training set (larger N, more stable estimates)
# Also report test set mean |SHAP| as cross-check

mean_abs_train <- colMeans(abs(S_train), na.rm = TRUE)
mean_abs_test  <- colMeans(abs(S_test),  na.rm = TRUE)

# Normalised weights (sum to 1) — used for signature scoring downstream
w_abs    <- mean_abs_train / sum(mean_abs_train)

# Sign of the average SHAP contribution; interpretation requires review.
# This sign alone does not establish a protein-dose association direction.
# See README.md.
m_shap   <- colMeans(S_train, na.rm = TRUE)
dir_sign <- sign(m_shap)
w_signed <- as.numeric(dir_sign * w_abs)
names(w_signed) <- names(w_abs)

sig <- data.frame(
  feature          = names(w_abs),
  weight_abs       = as.numeric(w_abs),
  weight_signed    = as.numeric(w_signed),
  mean_shap_train  = as.numeric(m_shap),
  mean_abs_train   = as.numeric(mean_abs_train),
  mean_abs_test    = as.numeric(mean_abs_test[names(mean_abs_train)]),
  direction        = ifelse(m_shap > 0, "+", ifelse(m_shap < 0, "-", "0")),
  stringsAsFactors = FALSE
)

sig <- sig[order(-sig$weight_abs), ]

# Cumulative gain check
sig$cumweight <- cumsum(sig$weight_abs)
cat("\nProteins needed for cumulative weight thresholds:\n")
cat(sprintf("  50%%: %d proteins\n", which(sig$cumweight >= 0.50)[1]))
cat(sprintf("  80%%: %d proteins\n", which(sig$cumweight >= 0.80)[1]))
cat(sprintf("  90%%: %d proteins\n", which(sig$cumweight >= 0.90)[1]))

cat("\nTop 10 proteins by |SHAP|:\n")
print(head(sig[, c("feature", "weight_abs", "mean_shap_train", "direction")], 10))

# ------------------- Save outputs -----------------------------------
# Full SHAP matrices (one row per individual, one col per protein)
S_train_out        <- S_train
S_train_out$split  <- "train"
S_train_out$idx    <- train_idx

S_test_out         <- S_test
S_test_out$split   <- "test"
S_test_out$idx     <- test_idx

data.table::fwrite(S_train_out, file.path(out_dir, "shap_train.csv"))
data.table::fwrite(S_test_out,  file.path(out_dir, "shap_test.csv"))
data.table::fwrite(sig,         file.path(out_dir, "signature_SHAP.csv"))

cat("\n========== SHAP outputs saved ==========\n")
cat(sprintf("  shap_train.csv      (%d obs x %d proteins)\n", nrow(S_train), ncol(S_train)))
cat(sprintf("  shap_test.csv       (%d obs x %d proteins)\n", nrow(S_test),  ncol(S_test)))
cat(sprintf("  signature_SHAP.csv  (%d proteins in signature)\n", nrow(sig)))
cat(sprintf("  Run folder: %s\n", out_dir))
cat(sprintf("Completed: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
