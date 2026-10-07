#'---- Syntax for ProPlant Project'#
#' Created by Peggy Ler (Pler) - 1st Oct 2025
#' Filename..: 1.1_Preprocess_ENET_uPDI_with_filter.R
#' Updated...:  by "Pler"

#' Purpose...: Preparing the data for elastic net
#' Protein filter: training-set LIMMA FDR < 0.05
#' 1 - split to training and test set
#' Proteins are assumed preprocessed upstream; no protein scaling here.
#' 3 - create model.matrix
#' 5 - cross validates, pick alpha
#' 6 - refit on training set with best alpha
#' 7 - evaluate the fitted model on the held-out test set
#' --------------------------------------------------------------------
#' Data used...: Controls only and no neuro cases
#' Using R studio 2024.04.02
#' -- Main program ----------------------------------------------------
#'
# USER PATHS: edit these locally; run from the folder containing README.md.
DATA_FILE <- "data/subcohort_20251106.RData"
OUTPUT_ROOT <- "outputs"

suppressPackageStartupMessages({
  library(limma)
  library(splines)
  library(dplyr)
  library(tibble)
})

# Load data, paths
#-------------
# Load only the expected data object into the analysis workspace.
input_data <- new.env(parent = baseenv())
load(DATA_FILE, envir = input_data)
if (!exists("subcohort", envir = input_data, inherits = FALSE) ||
    !is.data.frame(input_data$subcohort)) {
  stop("The RData file must contain a data.frame named subcohort.")
}
subcohort <- input_data$subcohort
rm(input_data)

outdir <- file.path(OUTPUT_ROOT, "updi")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

df <- subcohort
nrow(df)
#4372

# set seed
# --------
set.seed(258)
train_prop <- 0.8
fdr_thresh <- 0.05  # change if you want stricter/looser filter for proteins

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

# Keep analysis sample
# ---------------------
model_covars <- c("Sex","Age_Recr","Updi_NoVegoils", "edu2",
                  "alc_freq","smoke","PA_total_MET_quint","Center")
keep_rows <- stats::complete.cases(df[, model_covars])
df <- df[keep_rows, , drop = FALSE]

# Drop empty levels
df$Center <- droplevels(df$Center)
df$edu2 <- droplevels(df$edu2)
df$alc_freq <- droplevels(df$alc_freq)
df$smoke <- droplevels(df$smoke)
df$PA_total_MET_quint <- droplevels(df$PA_total_MET_quint)

# Identify protein columns
# -------------------------
protein_cols <- grep("^seq", names(df), value = TRUE)
is_num <- vapply(df[protein_cols], is.numeric, logical(1))
protein_cols <- protein_cols[is_num]
stopifnot(length(protein_cols) > 0)

# Split data in train and test
# ------------------------------
n <- nrow(df)

train_idx <- sample.int(n, size = floor(train_prop * n))
test_idx  <- setdiff(seq_len(n), train_idx)

train_df <- df[train_idx, , drop = FALSE]
test_df  <- df[test_idx,  , drop = FALSE]

# TRAIN-only limma to filter proteins
# --------------------------------------------------------------
X_train_limma <- model.matrix(
  ~ Updi_NoVegoils + Age_Recr + Sex + edu2 +
    alc_freq + smoke + PA_total_MET_quint + Center,
  data = train_df
)

# Proteins matrix for limma: features in rows, samples in columns (TRAIN ONLY)
#---------------------------------------------------------------------------------
prot_mat_train <- as.matrix(train_df[, protein_cols, drop = FALSE])
storage.mode(prot_mat_train) <- "double"
Y_train <- t(prot_mat_train)                          # rows = proteins
rownames(Y_train) <- colnames(prot_mat_train)

# Filter zero-variance in TRAIN
# -----------------------------
keep_var <- apply(Y_train, 1, stats::var, na.rm = TRUE) > 0
Y_train <- Y_train[keep_var, , drop = FALSE]
prot_kept_var <- rownames(Y_train)

# Fit limma on TRAIN
# --------------------
fit_lm <- lmFit(Y_train, X_train_limma)
fit_lm <- eBayes(fit_lm, robust = TRUE)

# Extract TRAIN-only Updi_NoVegoils results
# --------------------------------
coef_name <- "Updi_NoVegoils"
stopifnot(coef_name %in% colnames(X_train_limma))
tt <- topTable(fit_lm, coef = coef_name, number = Inf, sort.by = "none") %>%
  rownames_to_column(var = "protein")

#filter proteins by train only FDR
#----------------------------------------
prot_keep <- tt %>%
  filter(adj.P.Val < fdr_thresh) %>%
  pull(protein) %>%
  unique()

cat("TRAIN-only limma filter | FDR <", fdr_thresh, " | kept proteins: ", length(prot_keep), "\n")
#16 proteins

if (length(prot_keep) == 0) {
  stop("No proteins passed the TRAIN-only limma FDR threshold. Increase fdr_thresh.")
}

# Ensure kept proteins exist in both train and test data columns
# -------------------------------------------------------------------
prot_keep <- intersect(prot_keep, intersect(colnames(train_df), colnames(test_df)))

#' saving script for next step tuning
#' ----------------------------------
run_id <- paste0(
  "ENet_Updi_NoVegoils_filtFDR",
  gsub("\\.", "", as.character(fdr_thresh * 100))
)

saveRDS(
  list(
    df           = df,          # processed df with factors/scaled Updi_NoVegoils & Age_Recr
    train_idx    = train_idx,
    test_idx     = test_idx,
    protein_cols = protein_cols,
    prot_keep    = prot_keep,
    model_covars = model_covars,
    run_id       = run_id,
    fdr_thresh   = fdr_thresh,
    train_prop   = train_prop,
    seed_used    = 258
  ),
  file = file.path(outdir, paste0("prep_objects_", run_id, ".rds"))
)

write.csv(tt,
          file.path(outdir, paste0("trainonly_limma_Updi_NoVegoils_", format(Sys.Date(), "%Y%m%d"), ".csv")),
          row.names = FALSE
)

cat("Prep script done. Saved prep_objects and limma table.\n")
