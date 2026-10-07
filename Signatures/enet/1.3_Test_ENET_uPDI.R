#'---- Syntax for ProPlant Project'#
#' Created by Peggy Ler (Pler) - 1st Oct 2025
#' Filename..: 1.3_Test_Lasso_uPDI.R
#' Updated...:  by "Pler"

#' Purpose...: Load saved model, rebuild test matrix, predict and evaluate
#' --------------------------------------------------------------------
#' Data used...: EPIC subcohort (object: subcohort)
#' Using R studio 2024.04.02
#' -- Main program ----------------------------------------------------
#'
# USER PATHS: edit these locally; run from the folder containing README.md.
DATA_FILE <- "data/subcohort_20251106.RData"
OUTPUT_ROOT <- "outputs"
suppressPackageStartupMessages({
  library(glmnet)
  library(splines)
  library(dplyr)
  library(tibble)
  library(readr)
})

# Paths
# -----
outdir    <- file.path(OUTPUT_ROOT, "updi")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# load model bundle
model_rds <- list.files(outdir, pattern = "^elasticnet_model_.*\\.rds$", full.names = TRUE)
if (length(model_rds) != 1) stop("Can't uniquely find model rds in outdir.")
mod <- readRDS(model_rds)

fit_final    <- mod$fit_final
best_alpha   <- mod$best_alpha
best_lambda  <- mod$best_lambda
pf_saved     <- mod$pf
colnames_X   <- mod$colnames_X
covar_terms  <- mod$covar_terms
run_id       <- mod$run_id
train_idx    <- mod$train_idx
test_idx     <- mod$test_idx
protein_cols <- mod$protein_cols
age_df       <- mod$age_df

# Reload data to reconstruct TEST set
# ----------------------------------
data_path <- DATA_FILE
# Load only the expected data object into the analysis workspace.
input_data <- new.env(parent = baseenv())
load(data_path, envir = input_data)
if (!exists("subcohort", envir = input_data, inherits = FALSE) ||
    !is.data.frame(input_data$subcohort)) {
  stop("The RData file must contain a data.frame named subcohort.")
}
subcohort <- input_data$subcohort
rm(input_data)

df <- subcohort

df$edu2 <- dplyr::case_when(
  df$L_School %in% c(0, 1)              ~ "Below secondary",
  df$L_School %in% c(2, 3, 4)           ~ "Secondary/vocational and above",
  df$L_School == 5 | is.na(df$L_School) ~ "Not specified",
  TRUE                                   ~ NA_character_
)
df$edu2 <- factor(
  df$edu2,
  levels = c("Below secondary",
             "Secondary/vocational and above",
             "Not specified")
)

table(df$edu2, useNA = "ifany")
# Below secondary Secondary/vocational and above                  Not specified
# 2242                           2044                             86

cat("\nedu2 distribution (reference = Below secondary):\n")
print(table(df$edu2, useNA = "ifany"))

df$smoke <- factor(df$smoke)
df$alc_freq <- factor(df$alc_freq)
df$smoke <- relevel(df$smoke, ref = "Never")
df$alc_freq <- relevel(df$alc_freq, ref = "Never")

# ---------- covariates / scaling ----------
df$Updi_NoVegoils <- as.numeric(scale(df$Updi_NoVegoils))  # center+scale
df$Age_Recr <- as.numeric(scale(df$Age_Recr))
df$Center <- as.factor(df$Center)
df$PA_total_MET_quint <- as.factor(df$PA_total_MET_quint)
df$Sex <- as.factor(df$Sex)

# Recreate same processing used in 01_prep.R
model_covars <- c("Sex","Age_Recr","Updi_NoVegoils","edu2",
                  "alc_freq","smoke","PA_total_MET_quint","Center")

keep_rows <- stats::complete.cases(df[, model_covars])
df <- df[keep_rows, , drop = FALSE]

# Outcome and age were scaled above, matching the 1.1 preprocessing script.
# Duplicate post-filter scaling removed; see README.md.

df$Center <- droplevels(df$Center)
df$edu2 <- droplevels(df$edu2)
df$alc_freq <- droplevels(df$alc_freq)
df$smoke <- droplevels(df$smoke)
df$PA_total_MET_quint <- droplevels(df$PA_total_MET_quint)

# after complete.cases filtering + droplevels
train_df <- df[train_idx, , drop = FALSE]
test_df  <- df[test_idx,  , drop = FALSE]
# build covariates for TEST
Xcov_test <- model.matrix(
  ~ Age_Recr + Sex + edu2 +
    alc_freq + smoke + PA_total_MET_quint + Center,
  data = test_df
)[, -1, drop = FALSE]

# proteins (already scaled upstream)
test_prot_mat <- as.matrix(test_df[, protein_cols, drop = FALSE])

# combine
X_test_raw <- cbind(Xcov_test, test_prot_mat)

# align columns of TEST to training columns
align_newx <- function(X_new, cols_ref) {
  X_out <- matrix(0, nrow = nrow(X_new), ncol = length(cols_ref))
  colnames(X_out) <- cols_ref
  common_cols <- intersect(cols_ref, colnames(X_new))
  X_out[, common_cols] <- X_new[, common_cols, drop = FALSE]
  X_out
}
X_test <- align_newx(X_test_raw, colnames_X)

# response on test
y_test <- test_df$Updi_NoVegoils

## ----------------- Predict on TEST -----------------
pred_test <- as.numeric(predict(fit_final, newx = X_test))

rmse <- function(y, p) sqrt(mean((y - p)^2))
mae  <- function(y, p) mean(abs(y - p))
r2   <- function(y, p) cor(y, p)^2
spmn <- function(y, p) suppressWarnings(cor(y, p, method = "spearman"))

metrics_test <- tibble(
  split     = "TEST",
  alpha     = best_alpha,
  lambda    = best_lambda,
  RMSE      = rmse(y_test, pred_test),
  MAE       = mae(y_test, pred_test),
  R2        = r2(y_test, pred_test),
  SpearmanR = spmn(y_test, pred_test)
)

print(metrics_test)
# split alpha  lambda  RMSE   MAE    R2 SpearmanR
# <chr> <dbl>   <dbl> <dbl> <dbl> <dbl>     <dbl>
#   1 TEST    0.7 0.00103 0.845 0.665 0.279     0.509

write_csv(
  metrics_test,
  file.path(outdir, paste0("elasticnet_", run_id, "_TEST_metrics.csv"))
)

# Save predictions with IDs if available
id_col <- intersect(c("Idepic"), names(test_df))
id_vec <- if (length(id_col)) test_df[[id_col[1]]] else seq_len(nrow(test_df))

pred_out <- tibble(
  id        = id_vec,
  observed  = y_test,
  predicted = pred_test
)

write_csv(
  pred_out,
  file.path(outdir, paste0("elasticnet_", run_id, "_TEST_predictions.csv"))
)
