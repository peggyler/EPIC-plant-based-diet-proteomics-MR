#'---- Syntax for ProPlant Project'#
#' Created by Peggy Ler (Pler) - 20th Oct 2025
#' Filename..: Limma_uPDI.R
#' Updated...: by "Pler"
#' Purpose...: Full-subcohort LIMMA associations, with a linear age term
#' --------------------------------------------------------------------
#' Data used...: EPIC subcohort (object: subcohort)
#' Using R studio 2024.04.02
#' -- Main program ----------------------------------------------------
# USER PATHS: edit these locally; run from the folder containing README.md.
DATA_FILE <- "data/subcohort_20251106.RData"
OUTPUT_ROOT <- "outputs"
outdir <- file.path(OUTPUT_ROOT, "updi")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

suppressPackageStartupMessages({
  library(limma)
  library(data.table)
  library(readr)
  library(dplyr)
  library(tibble)
  library(splines)
})

# load on workbench
# Load only the expected data object into the analysis workspace.
input_data <- new.env(parent = baseenv())
load(DATA_FILE, envir = input_data)
if (!exists("subcohort", envir = input_data, inherits = FALSE) ||
    !is.data.frame(input_data$subcohort)) {
  stop("The RData file must contain a data.frame named subcohort.")
}
subcohort <- input_data$subcohort
rm(input_data)

# identify protein columns
df <- subcohort
nrow(df) #4372

protein_cols <- grep("^seq", names(df), value = TRUE)
#7285

# keep numeric proteins only
is_num <- vapply(df[protein_cols], is.numeric, logical(1))
protein_cols <- protein_cols[is_num]
stopifnot(length(protein_cols) > 0)

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

# ---------- make a single analysis sample (same for all proteins) ----------
model_covars <- c("Sex","Age_Recr","Updi_NoVegoils", "edu2",
                  "alc_freq","smoke","PA_total_MET_quint",
                  "Center")
keep_rows <- stats::complete.cases(df[, model_covars])
df   <- df[keep_rows, , drop = FALSE]
nrow(df)

## create matrix
prot_mat <- as.matrix(df[, protein_cols])
storage.mode(prot_mat) <- "double"
colnames(prot_mat) <- protein_cols

Y <- t(prot_mat) #transpose so that rows are proteins, features in rows instead
rownames(Y)  <- protein_cols
keep <- apply(Y, 1, stats::var, na.rm = TRUE) > 0 #keeping variances that are not missing or 0
Y <- Y[keep, , drop = FALSE]

# Design  model matrix
# -------------------------
X2 <- model.matrix(~ Sex + Age_Recr + Updi_NoVegoils + edu2 +
                     alc_freq + smoke + PA_total_MET_quint + Center,
                   data = df)

# If any columns are non-estimable, drop them up front
ne <- limma::nonEstimable(X2)
if (!is.null(ne)) {
  X2 <- X2[, setdiff(colnames(X2), ne), drop = FALSE]
}

# Sanity check, samples should match
stopifnot(ncol(Y) == nrow(X2))

# ----------------------- limma fit -----------------------
fit2 <- limma::lmFit(Y, X2)
fit2 <- limma::eBayes(fit2, robust = TRUE)

# Extract results for Updi_Novegoils
term_name <- "Updi_NoVegoils"
coef_idx  <- match(term_name, colnames(X2))
stopifnot(!is.na(coef_idx))

tt <- limma::topTable(fit2, coef = coef_idx, number = Inf, sort.by = "none")

tt <- tt %>%
  rownames_to_column(var = "protein") %>%
  mutate(
    q_BH  = p.adjust(P.Value, method = "BH"),
    q_BY  = p.adjust(P.Value, method = "BY"),
    p_bonf = p.adjust(P.Value, method = "bonferroni")
  ) %>%
  relocate(protein)

# SAVE
# ----
datestr <- format(Sys.Date(), "%Y%m%d")
outfile <- paste0("Updinovegoils_limma_", datestr, ".csv")
out_path <- file.path(outdir, outfile)
write.csv(tt, out_path, row.names = FALSE)

sum(tt$adj.P.Val < 0.05) #51 proteins
