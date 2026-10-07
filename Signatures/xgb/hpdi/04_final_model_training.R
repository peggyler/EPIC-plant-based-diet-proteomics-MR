#'---- Syntax for ProPlant Project'#
#' Filename..: 04_final_model_training.R
#' Purpose...: Train final XGBoost model using tuned parameters.
#'             Evaluate on held-out test set.
#'             Extract protein signature (non-zero gain features).
#'             Save model, signature, predictions, and performance metrics.
#' Reads  : current_run_dir.txt -> {out_dir}/final_params.rds
#'          {out_dir}/data_splits.rds
#' Writes : {out_dir}/final_model.rds
#'          {out_dir}/protein_signature.csv   <- proteins with Gain > 0
#'          {out_dir}/feature_importance.csv  <- full importance table
#'          {out_dir}/test_predictions.csv    <- y_test vs y_pred
#'          {out_dir}/model_performance.csv   <- RMSE, R2, Spearman rho
#' --------------------------------------------------------------------

# USER PATHS: edit these locally; run from the folder containing README.md.
DATA_FILE <- "data/subcohort_20251106.RData"
OUTPUT_ROOT <- "outputs"

# =================== USER CONFIG =====================================
OUTCOME_COL <- "Hpdi_NoVegoils"
OUT_BASE <- file.path(OUTPUT_ROOT, "hpdi", "gbm")
if (!dir.exists(OUT_BASE)) dir.create(OUT_BASE, recursive = TRUE)
OUT_BASE <- normalizePath(OUT_BASE, winslash = "/", mustWork = TRUE)
SOURCE_DIR <- file.path("signatures", "gbm", "hpdi")
LOAD_SPLITS <- TRUE    # always TRUE for scripts 02-04
# =====================================================================

# ------------------- Locate run folder ------------------------------
run_pointer <- file.path(OUT_BASE, "current_run_dir.txt")
if (!file.exists(run_pointer)) stop("Run Stage A (01) first.")
out_dir <- trimws(readLines(run_pointer, warn = FALSE))
if (length(out_dir) != 1L || !nzchar(out_dir) || !dir.exists(out_dir)) {
  stop("current_run_dir.txt must name one existing run directory.")
}
cat(sprintf("Using run folder: %s\n", out_dir))

# ------------------- Load final params ------------------------------
final_params <- readRDS(file.path(out_dir, "final_params.rds"))
cat("\nFinal parameters loaded:\n")
print(unlist(final_params))

# ------------------- Source shared prep -----------------------------
source(file.path(SOURCE_DIR, "00_data_prep.R"), local = TRUE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_finalmodel.txt"))

# ------------------- Build XGB params list --------------------------
xgb_params <- list(
  objective        = "reg:squarederror",
  eval_metric      = "rmse",
  tree_method      = "hist",
  max_bin          = 128,
  nthread          = 16,
  seed             = 258L,
  eta              = final_params$eta,
  subsample        = final_params$subsample,
  max_depth        = final_params$max_depth,
  min_child_weight = final_params$min_child_weight,
  lambda           = final_params$lambda,
  alpha            = final_params$alpha,
  colsample_bytree = final_params$colsample_bytree
)

# ------------------- Train final model ------------------------------
# Train on full training set using nrounds from early stopping
cat(sprintf("\nTraining final model: nrounds=%d\n", final_params$nrounds))

set.seed(258)
final_model <- xgb.train(
  params  = xgb_params,
  data    = dtrain,
  nrounds = final_params$nrounds,
  verbose = 1,
  print_every_n = 100
)

# ------------------- Test set evaluation ----------------------------
y_pred <- predict(final_model, dtest)

rmse_test    <- sqrt(mean((y_test - y_pred)^2))
ss_res       <- sum((y_test - y_pred)^2)
ss_tot       <- sum((y_test - mean(y_test))^2)
r2_test      <- 1 - (ss_res / ss_tot)
spearman_rho <- cor(y_test, y_pred, method = "spearman")

cat("\n========== Test set performance ==========\n")
cat(sprintf("RMSE        : %.4f\n", rmse_test))
cat(sprintf("R2          : %.4f\n", r2_test))
cat(sprintf("Spearman rho: %.4f\n", spearman_rho))

performance <- data.frame(
  outcome      = OUTCOME_COL,
  rmse_test    = rmse_test,
  r2_test      = r2_test,
  spearman_rho = spearman_rho,
  nrounds      = final_params$nrounds,
  n_train      = length(y_train),
  n_test       = length(y_test)
)

# ------------------- Feature importance -----------------------------
imp_full <- xgb.importance(model = final_model)
imp_full <- as.data.frame(imp_full)

# Protein signature: only features with Gain > 0, only seq.* columns
protein_sig <- imp_full[grepl("^seq", imp_full$Feature), ]
protein_sig <- protein_sig[protein_sig$Gain > 0, ]
protein_sig <- protein_sig[order(-protein_sig$Gain), ]

cat(sprintf("\nProtein signature: %d proteins selected (Gain > 0)\n", nrow(protein_sig)))
cat("Top 10 proteins by Gain:\n")
print(head(protein_sig[, c("Feature", "Gain", "Cover", "Frequency")], 10))

# ------------------- Predictions table ------------------------------
preds_df <- data.frame(
  idx    = test_idx,
  y_obs  = y_test,
  y_pred = y_pred
)

# ------------------- Save all outputs -------------------------------
saveRDS(final_model, file.path(out_dir, "final_model.rds"))
data.table::fwrite(protein_sig,  file.path(out_dir, "protein_signature.csv"))
data.table::fwrite(imp_full,     file.path(out_dir, "feature_importance.csv"))
data.table::fwrite(preds_df,     file.path(out_dir, "test_predictions.csv"))
data.table::fwrite(performance,  file.path(out_dir, "model_performance.csv"))

cat("\n========== Files saved ==========\n")
cat(sprintf("  final_model.rds\n"))
cat(sprintf("  protein_signature.csv  (%d proteins)\n", nrow(protein_sig)))
cat(sprintf("  feature_importance.csv (%d features total)\n", nrow(imp_full)))
cat(sprintf("  test_predictions.csv\n"))
cat(sprintf("  model_performance.csv\n"))
cat(sprintf("  Run folder: %s\n", out_dir))
cat(sprintf("\nCompleted: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
