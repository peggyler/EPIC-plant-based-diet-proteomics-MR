#'---- Syntax for ProPlant Project'#
#' Created by Peggy Ler (Pler) - 1st Oct 2025
#' Filename..: 1.2_Train_Lasso_hPDI.R
#' Updated...:  by "Pler"

#' Purpose...: Elastic-net protein signatures after training-set LIMMA filtering
#' Protein filter: training-set LIMMA FDR < 0.05
#' 1 - split to training and test set
#' Proteins are assumed preprocessed upstream; no protein scaling here.
#' 3 - create model.matrix
#' 5 - cross validates, pick alpha
#' 6 - refit on training set with best alpha
#' --------------------------------------------------------------------
#' Data used...: controls, no neuro cases
#' Using R studio 2024.04.02
#' -- Main program ----------------------------------------------------
#'
# USER PATHS: edit these locally; run from the folder containing README.md.
DATA_FILE <- "data/subcohort_20251106.RData"
OUTPUT_ROOT <- "outputs"

suppressPackageStartupMessages({
  library(splines)
  library(glmnet)
  library(dplyr)
  library(tibble)
  library(readr)
})

# Paths
# -----
outdir    <- file.path(OUTPUT_ROOT, "hpdi")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
prep_rds  <- list.files(outdir, pattern = "^prep_objects_.*\\.rds$", full.names = TRUE)
if (length(prep_rds) != 1) stop("Can't uniquely find prep_objects_*.rds in outdir.")
prep_obj  <- readRDS(prep_rds)

df           <- prep_obj$df
train_idx    <- prep_obj$train_idx
test_idx     <- prep_obj$test_idx
protein_cols <- prep_obj$prot_keep  # IMPORTANT: filtered proteins only
run_id       <- paste0(prep_obj$run_id, "_lamMIN")  # we'll tag lamMIN

# Split train and test
train_df <- df[train_idx, , drop = FALSE]
test_df  <- df[test_idx,  , drop = FALSE]  # not used for tuning, but we’ll evaluate later

# -----------------------------------------------------------------------------
# Build glmnet matrices
# Build the covariate design using train_df.
# -----------------------------------------------------------------------------
Xcov_train <- model.matrix(
  ~  Age_Recr+ Sex + edu2 +
    alc_freq + smoke + PA_total_MET_quint + Center,
  data = train_df
)[, -1, drop = FALSE]  # drop intercept

covar_terms <- colnames(Xcov_train)

#
train_prot_mat <- as.matrix(train_df[, protein_cols, drop = FALSE])

# Combine covariates (penalty=0) + proteins (penalty=1)
X_train <- cbind(Xcov_train, train_prot_mat)

y_train <- train_df$Hpdi_NoVegoils

# Penalty factors, where 0 for covariates and 1 for proteins
pf <- ifelse(colnames(X_train) %in% covar_terms, 0, 1)


# REVIEW: no CV seed or shared foldid is set in the supplied script.
# See README.md before reproducing a reported fit.
# Cross validate over alphas
# ---------------------------
alphas <- seq(0.1, 0.9, by = 0.2)   # elastic-net mixing grid
cvmods <- vector("list", length(alphas))
cv_err <- numeric(length(alphas))

for (i in seq_along(alphas)) {
  a <- alphas[i]
  fit_cv <- cv.glmnet(
    x = X_train,
    y = y_train,
    family = "gaussian",
    alpha = a,
    penalty.factor = pf,
    nfolds = 10,
    standardize = FALSE,
    type.measure = "mse"
  )
  cvmods[[i]] <- fit_cv
  cv_err[i] <- min(fit_cv$cvm) #error at minimum lambda for that alpha
  cat(sprintf("[alpha=%.2f] lambda.min=%.6g | CV MSE(min)=%.4f\n",
              a, fit_cv$lambda.min, cv_err[i]))
}

# [alpha=0.10] lambda.min=0.053587 | CV MSE(min)=0.8739
# [alpha=0.30] lambda.min=0.031215 | CV MSE(min)=0.8737
# [alpha=0.50] lambda.min=0.0155491 | CV MSE(min)=0.8756
# [alpha=0.70] lambda.min=0.0121894 | CV MSE(min)=0.8733
# [alpha=0.90] lambda.min=0.0125328 | CV MSE(min)=0.8804

best_idx   <- which.min(cv_err)
best_alpha <- alphas[best_idx]
best_cv    <- cvmods[[best_idx]]
best_lambda <- best_cv$lambda.min

cat(sprintf("Chosen alpha=%.2f | lambda.min=%.6g | CV MSE=%.4f\n",
            best_alpha, best_lambda, cv_err[best_idx]))
#Chosen alpha=0.70 | lambda.min=0.0121894 | CV MSE=0.8733

# ----------------------------------------
# Final fit on full training set
# ----------------------------------------
fit_final <- glmnet(
  x = X_train,
  y = y_train,
  family = "gaussian",
  alpha = best_alpha,
  lambda = best_lambda,      # --- CHANGED: lock to lambda.min
  penalty.factor = pf,
  standardize = FALSE
)

#' Extract trained signature
#' --------------------------
coef_train <- as.matrix(coef(fit_final))[,1]
names(coef_train) <- rownames(as.matrix(coef(fit_final)))

sig_train_all <- names(coef_train)[coef_train != 0]
sig_train_prot <- setdiff(sig_train_all, c("(Intercept)", covar_terms))

cat("Number of proteins with nonzero coef in TRAIN model:", length(sig_train_prot), "\n")
#108

coef_tbl <- tibble(
  feature      = setdiff(names(coef_train), "(Intercept)"),
  coef         = as.numeric(coef_train[setdiff(names(coef_train), "(Intercept)")]),
  abs_coef     = abs(coef),
  feature_type = ifelse(feature %in% covar_terms, "covariate", "protein")
) %>% arrange(desc(abs_coef))

write_csv(coef_tbl,
          file.path(outdir, paste0("elasticnet_", run_id, "_coef_full.csv")))
write_csv(
  tibble(selected_protein = sig_train_prot),
  file.path(outdir, paste0("elasticnet_", run_id, "_selected_proteins_train.csv"))
)

##' ---------------------------------------------------------------
##' BOOTSTRAP stabillity of protein signatires in training set only
##' ---------------------------------------------------------------
set.seed(42)
B <- 1000   # adjust if you need faster/slower

feat_names <- colnames(X_train)
sel_counts <- setNames(integer(length(feat_names)), feat_names)

for (b in 1:B) {
  idx_b <- sample.int(nrow(X_train), replace = TRUE)

  fit_b <- glmnet(
    x = X_train[idx_b, , drop = FALSE],
    y = y_train[idx_b],
    family = "gaussian",
    alpha = best_alpha,
    lambda = best_lambda,
    penalty.factor = pf,
    standardize = FALSE
  )

  beta_b <- as.matrix(coef(fit_b))[,1]  # includes intercept
  nz <- names(beta_b)[beta_b != 0]
  nz <- setdiff(nz, "(Intercept)")
  sel_counts[nz] <- sel_counts[nz] + 1L
}

stab_df <- tibble(
  feature      = names(sel_counts),
  sel_freq     = sel_counts / B,
  feature_type = ifelse(feature %in% covar_terms, "covariate", "protein")
) %>%
  arrange(desc(sel_freq))

stable_proteins <- stab_df %>%
  filter(feature_type == "protein", sel_freq >= 0.6) %>%
  arrange(desc(sel_freq))

write_csv(
  stab_df,
  file.path(outdir, paste0("stability_bootstrap_", run_id, ".csv"))
)

write_csv(
  stable_proteins,
  file.path(outdir, paste0("stable_protein_signature_", run_id, ".csv"))
)

cat("Bootstrap complete. Proteins with sel_freq >= 0.6:",
    nrow(stable_proteins), "\n")

#Bootstrap complete. Proteins with sel_freq >= 0.6: 78

## ----------------- Save model bundle for scoring -----------------
model_bundle <- list(
  fit_final     = fit_final,
  best_alpha    = best_alpha,
  best_lambda   = best_lambda,
  pf            = pf,
  colnames_X    = colnames(X_train),  # required column order at prediction time
  covar_terms   = covar_terms,
  run_id        = run_id,
  train_idx     = train_idx,
  test_idx      = test_idx,
  protein_cols  = protein_cols
)

saveRDS(
  model_bundle,
  file = file.path(outdir, paste0("elasticnet_model_", run_id, ".rds"))
)

cat("Training script done. Model bundle saved.\n")
