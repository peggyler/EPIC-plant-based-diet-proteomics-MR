#'---- Syntax for ProPlant Project'#
#' Filename..: 02_stage_B_regularisation.R
#' Purpose...: Stage B tuning — lambda x alpha
#'             Reads best tree structure from Stage A.
#'             Saves best_AB_params.rds for Stage C to read.
#' Reads  : current_run_dir.txt  -> {out_dir}/best_A_params.rds
#'          {out_dir}/data_splits.rds
#' Writes : {out_dir}/results_stageB.csv
#'          {out_dir}/best_AB_params.rds
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

# ------------------- Load Stage A results ---------------------------
best_A <- readRDS(file.path(out_dir, "best_A_params.rds"))
cat(sprintf("Loaded Stage A best: max_depth=%d, min_child_weight=%d, RMSE=%.4f\n",
            best_A$max_depth, best_A$min_child_weight, best_A$rmse_A))

# ------------------- Source shared prep -----------------------------
source(file.path(SOURCE_DIR, "00_data_prep.R"), local = TRUE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_stageB.txt"))

# ------------------- Base XGB parameters ----------------------------
xgb_params <- list(
  objective   = "reg:squarederror",
  eval_metric = "rmse",
  tree_method = "hist",
  max_bin     = 128,
  nthread     = 16,
  seed        = 258L
)

# ------------------- Stage B grid -----------------------------------
# Tuning: lambda, alpha
# Fixed: best from Stage A + eta, subsample, colsample_bytree
# hPDI prior best: lambda=8, alpha=8 -> explore {6,8,10} x {6,8,10}
grid_B <- expand.grid(
  lambda           = c(6, 8, 10),
  alpha            = c(6, 8, 10),
  stringsAsFactors = FALSE
)

fixed_B <- list(
  eta              = best_A$eta,
  subsample        = best_A$subsample,
  colsample_bytree = best_A$colsample_bytree,
  max_depth        = best_A$max_depth,
  min_child_weight = best_A$min_child_weight
)

cat(sprintf("\nStage B: %d combos | tuning lambda x alpha\n", nrow(grid_B)))
cat("Fixed params: "); print(unlist(fixed_B))

# ------------------- CV function ------------------------------------
cv_one <- function(params_override, combo_idx, stage_tag = "B") {

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

# ------------------- Run Stage B ------------------------------------
results_list <- lapply(seq_len(nrow(grid_B)), function(i) {
  override <- c(as.list(grid_B[i, , drop = FALSE]), fixed_B)
  cv_one(params_override = override, combo_idx = i)
})

results_B <- data.table::rbindlist(results_list, fill = TRUE)
data.table::setorder(results_B, rmse_mean, proteins_used)

data.table::fwrite(results_B, file.path(out_dir, "results_stageB.csv"))

cat("\n========== Stage B complete ==========\n")
cat("Top 3 results:\n")
print(head(results_B[, .(lambda, alpha, rmse_mean, rmse_sd, proteins_used, best_iter)], 3))

# ------------------- Save best AB params ----------------------------
best_AB <- list(
  max_depth        = best_A$max_depth,
  min_child_weight = best_A$min_child_weight,
  eta              = best_A$eta,
  subsample        = best_A$subsample,
  colsample_bytree = best_A$colsample_bytree,
  lambda           = results_B$lambda[1],
  alpha            = results_B$alpha[1],
  best_iter_A      = best_A$best_iter_A,
  best_iter_B      = results_B$best_iter[1],
  rmse_A           = best_A$rmse_A,
  rmse_B           = results_B$rmse_mean[1]
)

saveRDS(best_AB, file.path(out_dir, "best_AB_params.rds"))
cat(sprintf("\nBest Stage B: lambda=%.1f, alpha=%.1f, RMSE=%.4f\n",
            best_AB$lambda, best_AB$alpha, best_AB$rmse_B))
cat(sprintf("Saved: best_AB_params.rds\n"))
cat(sprintf("Completed Stage B: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("\nNext step -> run 03_stage_C_colsampling.R\n"))
