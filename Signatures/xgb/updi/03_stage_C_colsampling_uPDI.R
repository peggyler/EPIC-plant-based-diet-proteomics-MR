#'---- Syntax for ProPlant Project'#
#' Filename..: 03_stage_C_colsampling_uPDI.R
#' Purpose...: Stage C tuning — colsample_bytree
#'             Reads best params from Stages A+B.
#'             Saves final_params.rds for model training.
#' Reads  : current_run_dir.txt -> {out_dir}/best_AB_params.rds
#'          {out_dir}/data_splits.rds
#' Writes : {out_dir}/results_stageC.csv
#'          {out_dir}/final_params.rds
#'          {out_dir}/tuning_summary_all_stages.csv
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
LOAD_SPLITS <- TRUE
# =====================================================================

# ------------------- Locate run folder ------------------------------
run_pointer <- file.path(OUT_BASE, "current_run_dir.txt")
if (!file.exists(run_pointer)) stop("Run Stage A (01) first.")
out_dir <- trimws(readLines(run_pointer, warn = FALSE))
if (length(out_dir) != 1L || !nzchar(out_dir) || !dir.exists(out_dir)) {
  stop("current_run_dir.txt must name one existing run directory.")
}
cat(sprintf("Using run folder: %s\n", out_dir))

# ------------------- Load Stage A+B results -------------------------
best_AB <- readRDS(file.path(out_dir, "best_AB_params.rds"))
cat("Loaded Stage A+B best:\n")
cat(sprintf("  max_depth=%d, min_child_weight=%d\n", best_AB$max_depth, best_AB$min_child_weight))
cat(sprintf("  lambda=%.1f, alpha=%.1f\n", best_AB$lambda, best_AB$alpha))
cat(sprintf("  RMSE after B: %.4f\n", best_AB$rmse_B))

# ------------------- Source shared prep -----------------------------
source(file.path(SOURCE_DIR, "00_data_prep_updi.R"), local = TRUE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_stageC.txt"))

# ------------------- Base XGB parameters ----------------------------
xgb_params <- list(
  objective   = "reg:squarederror",
  eval_metric = "rmse",
  tree_method = "hist",
  max_bin     = 128,
  nthread     = 16,
  seed        = 258L
)

# ------------------- Stage C grid -----------------------------------
# Tuning: colsample_bytree
# uPDI prior best: 0.90 -> explore {0.85, 0.90, 0.95}
grid_C <- expand.grid(
  colsample_bytree = c(0.85, 0.90, 0.95),
  stringsAsFactors = FALSE
)

fixed_C <- list(
  eta              = best_AB$eta,
  subsample        = best_AB$subsample,
  max_depth        = best_AB$max_depth,
  min_child_weight = best_AB$min_child_weight,
  lambda           = best_AB$lambda,
  alpha            = best_AB$alpha
)

cat(sprintf("\nStage C: %d combos | tuning colsample_bytree\n", nrow(grid_C)))
cat("Fixed params: "); print(unlist(fixed_C))

# ------------------- CV function ------------------------------------
cv_one <- function(params_override, combo_idx, stage_tag = "C") {

  params   <- modifyList(xgb_params, params_override)
  log_file <- file.path(out_dir, sprintf("cv_log_stage%s_combo%02d.txt", stage_tag, combo_idx))

  sink(log_file, append = TRUE)
  on.exit(sink(), add = TRUE)

  cat(sprintf("[Stage %s | Combo %d] %s | Start: %s\n",
              stage_tag, combo_idx,
              paste(names(params_override), round(unlist(params_override), 4),
                    sep = "=", collapse = ", "),
              format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

  t0 <- Sys.time()

  cv <- xgb.cv(
    params                = params,
    data                  = dtrain,
    nrounds               = 2000,
    folds                 = folds,
    early_stopping_rounds = 50,
    verbose               = 1,
    print_every_n         = 100,
    stratified            = FALSE
  )

  # Robust best_iter
  best_iter <- if (!is.null(cv$best_iteration) && length(cv$best_iteration) > 0) {
    cv$best_iteration
  } else {
    which.min(cv$evaluation_log$test_rmse_mean)
  }
  if (is.na(best_iter) || best_iter == 0) best_iter <- nrow(cv$evaluation_log)
  cat(sprintf("  best_iter resolved to: %d\n", best_iter))

  bst           <- xgb.train(params, dtrain, nrounds = best_iter, verbose = 0)
  imp           <- xgb.importance(model = bst)
  proteins_used <- nrow(imp)

  cat(sprintf("  Best iter=%d | RMSE=%.4f (SD=%.4f) | Proteins=%d | %.2f min\n\n",
              best_iter,
              cv$evaluation_log$test_rmse_mean[best_iter],
              cv$evaluation_log$test_rmse_std[best_iter],
              proteins_used,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  data.frame(
    stage            = stage_tag,
    combo_idx        = combo_idx,
    as.data.frame(params_override),
    best_iter        = best_iter,
    rmse_mean        = as.numeric(cv$evaluation_log$test_rmse_mean[best_iter]),
    rmse_sd          = as.numeric(cv$evaluation_log$test_rmse_std[best_iter]),
    proteins_used    = proteins_used,
    duration_min     = as.numeric(difftime(Sys.time(), t0, units = "mins")),
    stringsAsFactors = FALSE
  )
}

# ------------------- Run Stage C ------------------------------------
results_list <- lapply(seq_len(nrow(grid_C)), function(i) {
  override <- c(as.list(grid_C[i, , drop = FALSE]), fixed_C)
  cv_one(params_override = override, combo_idx = i)
})

results_C <- data.table::rbindlist(results_list, fill = TRUE)
data.table::setorder(results_C, rmse_mean, proteins_used)

data.table::fwrite(results_C, file.path(out_dir, "results_stageC.csv"))

cat("\n========== Stage C complete ==========\n")
cat("Top 3 results:\n")
print(head(results_C[, .(colsample_bytree, rmse_mean, rmse_sd, proteins_used, best_iter)], 3))

# ------------------- Final params -----------------------------------
final_params <- list(
  eta              = best_AB$eta,
  subsample        = best_AB$subsample,
  max_depth        = best_AB$max_depth,
  min_child_weight = best_AB$min_child_weight,
  lambda           = best_AB$lambda,
  alpha            = best_AB$alpha,
  colsample_bytree = results_C$colsample_bytree[1],
  nrounds          = results_C$best_iter[1]
)

saveRDS(final_params, file.path(out_dir, "final_params.rds"))

# ------------------- Combined tuning summary ------------------------
results_A <- data.table::fread(file.path(out_dir, "results_stageA.csv"))
results_B <- data.table::fread(file.path(out_dir, "results_stageB.csv"))

tuning_summary <- data.table::rbindlist(
  list(results_A, results_B, results_C),
  fill = TRUE
)
data.table::fwrite(tuning_summary, file.path(out_dir, "tuning_summary_all_stages.csv"))

# ------------------- Print final params -----------------------------
cat("\n========== Final optimal parameters ==========\n")
cat(sprintf("Outcome          : %s\n", OUTCOME_COL))
cat(sprintf("eta              : %.3f  [fixed]\n",    final_params$eta))
cat(sprintf("subsample        : %.2f  [fixed]\n",    final_params$subsample))
cat(sprintf("max_depth        : %d    [Stage A]\n",  final_params$max_depth))
cat(sprintf("min_child_weight : %d    [Stage A]\n",  final_params$min_child_weight))
cat(sprintf("lambda           : %.1f  [Stage B]\n",  final_params$lambda))
cat(sprintf("alpha            : %.1f  [Stage B]\n",  final_params$alpha))
cat(sprintf("colsample_bytree : %.2f  [Stage C]\n",  final_params$colsample_bytree))
cat(sprintf("nrounds          : %d    [early stop]\n", final_params$nrounds))
cat(sprintf("\nRMSE progression: A=%.4f -> B=%.4f -> C=%.4f\n",
            best_AB$rmse_A, best_AB$rmse_B, results_C$rmse_mean[1]))
cat(sprintf("\nSaved: final_params.rds, tuning_summary_all_stages.csv\n"))
cat(sprintf("Completed Stage C: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("\nNext step -> run 04_final_model_training_uPDI.R\n"))


# Historical pasted hPDI output below is NOT a uPDI result.
# ========= Stage C complete ==========
#   Top 3 results:
#   colsample_bytree rmse_mean    rmse_sd proteins_used best_iter
# <num>     <num>      <num>         <int>     <int>
# 1:             0.25 0.9655412 0.01677981           676       405
# 2:             0.30 0.9663524 0.01691306           766       470
# 3:             0.20 0.9667941 0.01637805           681       396
#
# ========== Final optimal parameters ==========
#   Outcome          : Hpdi_NoVegoils
# eta              : 0.030  [fixed]
# subsample        : 0.80  [fixed]
# max_depth        : 2    [Stage A]
# min_child_weight : 12    [Stage A]
# lambda           : 8.0  [Stage B]
# alpha            : 10.0  [Stage B]
# colsample_bytree : 0.25  [Stage C]
# nrounds          : 405    [early stop]
#
# RMSE progression: A=0.9656 -> B=0.9655 -> C=0.9655
