#'---- Syntax for ProPlant Project'#
#' Filename..: 01_stage_A_tree_structure.R
#' Purpose...: Stage A tuning — max_depth x min_child_weight
#'             Creates run folder, train/test splits, and CV folds.
#'             Saves best_A_params.rds for Stage B to read.
#' Reads  : subcohort_20251106.RData
#' Writes : {out_dir}/data_splits.rds
#'          {out_dir}/results_stageA.csv
#'          {out_dir}/best_A_params.rds
#' --------------------------------------------------------------------

# USER PATHS: edit these locally; run from the folder containing README.md.
DATA_FILE <- "data/subcohort_20251106.RData"
OUTPUT_ROOT <- "outputs"
set.seed(258)

# =================== USER CONFIG =====================================
OUTCOME_COL  <- "Hpdi_NoVegoils"
OUT_BASE <- file.path(OUTPUT_ROOT, "hpdi", "gbm")
if (!dir.exists(OUT_BASE)) dir.create(OUT_BASE, recursive = TRUE)
OUT_BASE <- normalizePath(OUT_BASE, winslash = "/", mustWork = TRUE)
SOURCE_DIR <- file.path("signatures", "gbm", "hpdi")
LOAD_SPLITS  <- FALSE    # FALSE = create new splits (always FALSE for script 01)
# =====================================================================

# ------------------- Output folder ----------------------------------
run_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(OUT_BASE, paste0("staged_run_", run_tag))
if (dir.exists(out_dir)) stop("Run folder already exists; retry in a new second.")
dir.create(out_dir, recursive = TRUE)
if (!file.access(out_dir, 2) == 0) stop("No write permission to ", out_dir)
cat(sprintf("Run folder: %s\n", out_dir))

# Save run_tag so subsequent scripts can find this folder
writeLines(out_dir, file.path(OUT_BASE, "current_run_dir.txt"))

# ------------------- Source shared prep -----------------------------
source(file.path(SOURCE_DIR, "00_data_prep.R"), local = TRUE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_stageA.txt"))

# ------------------- Base XGB parameters ----------------------------
xgb_params <- list(
  objective   = "reg:squarederror",
  eval_metric = "rmse",
  tree_method = "hist",
  max_bin     = 128,
  nthread     = 16,
  seed        = 258L
)

# ------------------- Stage A grid -----------------------------------
# Tuning: max_depth, min_child_weight
# Fixed at prior best values (hPDI broad search):
#   eta=0.03, subsample=0.8, colsample_bytree=0.25, lambda=8, alpha=8
grid_A <- expand.grid(
  max_depth        = c(2, 3),
  min_child_weight = c(8, 12, 16),
  stringsAsFactors = FALSE
)

fixed_A <- list(
  eta              = 0.03,
  subsample        = 0.8,
  colsample_bytree = 0.25,
  lambda           = 8,
  alpha            = 8
)

cat(sprintf("Stage A: %d combos | tuning max_depth x min_child_weight\n", nrow(grid_A)))
cat("Fixed params: "); print(unlist(fixed_A))

# ------------------- CV function ------------------------------------
cv_one <- function(params_override, combo_idx, stage_tag = "A") {

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

  # Robust fallback — cv$best_iteration can be NULL in some xgboost versions
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

# ------------------- Run Stage A ------------------------------------
results_list <- lapply(seq_len(nrow(grid_A)), function(i) {
  override <- c(as.list(grid_A[i, , drop = FALSE]), fixed_A)
  cv_one(params_override = override, combo_idx = i)
})

results_A <- data.table::rbindlist(results_list, fill = TRUE)
data.table::setorder(results_A, rmse_mean, proteins_used)

data.table::fwrite(results_A, file.path(out_dir, "results_stageA.csv"))

cat("\n========== Stage A complete ==========\n")
cat("Top 3 results:\n")
print(head(results_A[, .(max_depth, min_child_weight, rmse_mean, rmse_sd, proteins_used, best_iter)], 3))

# ------------------- Save best A params -----------------------------
best_A <- list(
  max_depth        = results_A$max_depth[1],
  min_child_weight = results_A$min_child_weight[1],
  # carry fixed params forward
  eta              = fixed_A$eta,
  subsample        = fixed_A$subsample,
  colsample_bytree = fixed_A$colsample_bytree,
  lambda           = fixed_A$lambda,
  alpha            = fixed_A$alpha,
  best_iter_A      = results_A$best_iter[1],
  rmse_A           = results_A$rmse_mean[1]
)

saveRDS(best_A, file.path(out_dir, "best_A_params.rds"))
cat(sprintf("\nBest Stage A: max_depth=%d, min_child_weight=%d, RMSE=%.4f\n",
            best_A$max_depth, best_A$min_child_weight, best_A$rmse_A))
cat(sprintf("Saved: best_A_params.rds\n"))
cat(sprintf("Completed Stage A: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("\nNext step -> run 02_stage_B_regularisation.R\n"))
