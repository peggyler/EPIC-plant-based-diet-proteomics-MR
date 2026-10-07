# Shared helper: invoked by the numbered GBM scripts.
# Run from the project root; see README.md.
#'---- Syntax for ProPlant Project'#
#' Filename..: 00_data_prep.R
#' Purpose...: Shared data loading, edu2 creation, preprocessing,
#'             design matrix construction.
#'             Sourced by 01-05 scripts — do not run standalone.
#' Requires in calling script:
#'   DATA_FILE    : character, local path to the prepared RData file
#'   OUTCOME_COL  : character, e.g. "Updi_NoVegoils"
#'   out_dir      : character, path to output folder for this run
#'   LOAD_SPLITS  : logical, TRUE = load saved splits, FALSE = create new ones
#' --------------------------------------------------------------------

# ------------------- Packages ---------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(xgboost)
  library(caret)
  library(data.table)
})

# ------------------- Load data --------------------------------------
# Load only the expected data object into the analysis workspace.
input_data <- new.env(parent = baseenv())
load(DATA_FILE, envir = input_data)
if (!exists("subcohort", envir = input_data, inherits = FALSE) ||
    !is.data.frame(input_data$subcohort)) {
  stop("The RData file must contain a data.frame named subcohort.")
}
subcohort <- input_data$subcohort
rm(input_data)
df <- subcohort
cat(sprintf("Loaded data: %d rows\n", nrow(df)))

# ------------------- Education variable (edu2) ----------------------
# Derived from L_School — avoids buggy L_School_1
# L_School: 0=None, 1=Primary, 2=Technical/professional,
#            3=Secondary, 4=Longer/university, 5=Not specified
# Reference level: "Below secondary"
# table(df$L_School, useNA = "ifany")
# 0    1    2    3    4      5 <NA>
# 659 1583  730  642  672   66   20

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

# ------------------- Complete cases ---------------------------------
model_covars <- c("Sex", "Age_Recr", OUTCOME_COL, "edu2",
                  "alc_freq", "smoke", "PA_total_MET", "Center")
keep_rows <- stats::complete.cases(df[, model_covars])
df        <- df[keep_rows, , drop = FALSE]
cat(sprintf("Rows after complete case filter: %d\n", nrow(df)))

# ------------------- Pre-process ------------------------------------
df$smoke <- factor(df$smoke)
df$alc_freq <- factor(df$alc_freq)
df$smoke <- relevel(df$smoke, ref = "Never")
df$alc_freq <- relevel(df$alc_freq, ref = "Never")

df[[OUTCOME_COL]] <- as.numeric(scale(df[[OUTCOME_COL]]))
df$Age_Recr       <- as.numeric(scale(df$Age_Recr))
df$Center         <- droplevels(as.factor(df$Center))
df$edu2           <- droplevels(df$edu2)
df$alc_freq       <- droplevels(df$alc_freq)
df$smoke          <- droplevels(df$smoke)
df$PA_total_MET   <- as.numeric(scale(df$PA_total_MET))

# ------------------- Protein columns --------------------------------
protein_cols <- grep("^seq", names(df), value = TRUE)
is_num       <- vapply(df[protein_cols], is.numeric, logical(1))
protein_cols <- protein_cols[is_num]
stopifnot(length(protein_cols) > 0)
cat(sprintf("Protein columns: %d\n", length(protein_cols)))

# ------------------- Train/test split & folds -----------------------
splits_path <- file.path(out_dir, "data_splits.rds")

if (!LOAD_SPLITS) {
  # Script 01 creates and saves splits
  set.seed(258)
  n         <- nrow(df)
  train_idx <- sample.int(n, size = floor(0.8 * n))
  test_idx  <- setdiff(seq_len(n), train_idx)

  # Stratified 5-fold CV on training outcome
  y_train_raw <- df[[OUTCOME_COL]][train_idx]
  K           <- 5
  br          <- unique(quantile(y_train_raw, probs = seq(0, 1, length.out = K + 1), na.rm = TRUE))
  if (length(br) < 3) br <- pretty(y_train_raw, n = K)
  y_bins      <- cut(y_train_raw, breaks = br, include.lowest = TRUE, labels = FALSE)
  set.seed(258)
  folds <- caret::createFolds(y_bins, k = K, list = TRUE, returnTrain = FALSE)

  saveRDS(list(train_idx = train_idx, test_idx = test_idx, folds = folds),
          splits_path)
  cat(sprintf("Splits saved: %d train, %d test\n", length(train_idx), length(test_idx)))

} else {
  # Scripts 02-04 load the same splits
  splits    <- readRDS(splits_path)
  train_idx <- splits$train_idx
  test_idx  <- splits$test_idx
  folds     <- splits$folds
  cat(sprintf("Splits loaded: %d train, %d test\n", length(train_idx), length(test_idx)))
}

# ------------------- Design matrices --------------------------------
cov_formula <- ~ Age_Recr + Sex + edu2 + alc_freq + smoke + PA_total_MET + Center
Xcov_full   <- sparse.model.matrix(cov_formula, data = df)
Xcov_full   <- Xcov_full[, colnames(Xcov_full) != "(Intercept)", drop = FALSE]

Xprot_full  <- Matrix::Matrix(as.matrix(df[, protein_cols, drop = FALSE]), sparse = TRUE)
storage.mode(Xprot_full@x) <- "double"
Xprot_full[is.na(Xprot_full)] <- 0

X_full  <- cbind(Xcov_full, Xprot_full)
X_train <- X_full[train_idx, ]
X_test  <- X_full[test_idx, ]
y_train <- df[[OUTCOME_COL]][train_idx]
y_test  <- df[[OUTCOME_COL]][test_idx]

dtrain <- xgb.DMatrix(X_train, label = y_train)
dtest  <- xgb.DMatrix(X_test,  label = y_test)

cat("Data prep complete.\n\n")
