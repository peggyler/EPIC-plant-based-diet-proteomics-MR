# Run from the repository root in a fresh R session. See README.md.
# EDIT LOCAL PATHS HERE. No study data or results are supplied.
PROJECT_ROOT <- "."
DATA_DIR <- file.path(PROJECT_ROOT, "data")
OUTPUT_ROOT <- file.path(PROJECT_ROOT, "outputs")
MANUSCRIPT_DIR <- file.path(OUTPUT_ROOT, "manuscript")
REFERENCE_DIR <- file.path(DATA_DIR, "reference")
SIGNATURE_SCRIPTS <- file.path(PROJECT_ROOT, "signatures")
CONTEXT_SCRIPTS <- file.path(PROJECT_ROOT, "contextualisation")
dir.create(MANUSCRIPT_DIR, recursive = TRUE, showWarnings = FALSE)

## sensitivity to test 
## proteamic scores with covariates

library(survival)
library(dplyr)
library(splines)
library(flexsurv)
library(ggplot2)
library(patchwork)
library(tibble)
library(forcats)

load(file.path(DATA_DIR, "df_with_scores_surv_20260425.RData"))
df <-df_sel

# z-scores: hPDI
df$GBMprot_Hpdi_z    <- as.numeric(scale(df$GBMprotSHAP_Hpdi))
df$GBMfull_Hpdi_z    <- as.numeric(scale(df$ GBMfull_Hpdi))
df$Hpdi_NoVegoils_z  <- as.numeric(scale(df$Hpdi_NoVegoils))

# z-scores: uPDI
df$GBMprot_Updi_z    <- as.numeric(scale(df$GBMprotSHAP_Updi))
df$GBMfull_Updi_z    <- as.numeric(scale(df$GBMfull_Updi))
df$Updi_NoVegoils_z  <- as.numeric(scale(df$Updi_NoVegoils))



fit_ca_gbmfull <- flexsurvreg(
  Surv(Age_Recr, Age_Exit_Cancer_1st, Cancer_Status) ~
    Sex + edu2 + alc4 + ever_smoked + PA + Updi_NoVegoils_z + GBMfull_Updi_z + Center_grp,
  data = df, weights = Weights, dist = "weibull"
)

fit_ca_gbmprot <- flexsurvreg(
  Surv(Age_Recr, Age_Exit_Cancer_1st, Cancer_Status) ~
    Sex + edu2 + alc4 + ever_smoked + PA + Updi_NoVegoils_z + GBMprot_Updi_z + Center_grp,
  data = df, weights = Weights, dist = "weibull"
)

fit_ca_gbmfull
# Updi_NoVegoils_z                    -2.28e-16  -8.21e-03  -1.56e-02  -8.42e-04   3.76e-03   9.92e-01   9.85e-01   9.99e-01
# GBMfull_Updi_z                      -1.05e-17  -1.38e-02  -2.45e-02  -3.05e-03   5.48e-03   9.86e-01   9.76e-01   9.97e-01

fit_ca_gbmprot
# Updi_NoVegoils_z                    -2.28e-16  -6.85e-03  -1.41e-02   4.27e-04   3.71e-03   9.93e-01   9.86e-01   1.00e+00
# GBMprot_Updi_z                      -2.04e-17  -1.43e-02  -2.26e-02  -5.96e-03   4.25e-03   9.86e-01   9.78e-01   9.94e-01

##t2d
#____

fit_t2d_gbmfull <- flexsurvreg(
  Surv(Age_Recr, Age_Exit_T2D, T2D_Status) ~
    Sex + edu2 + alc4 + ever_smoked + PA + Updi_NoVegoils_z + GBMfull_Updi_z + Center_grp,
  data = df, weights = Weights, dist = "gengamma"
)

fit_t2d_gbmprot <- flexsurvreg(
  Surv(Age_Recr, Age_Exit_T2D, T2D_Status) ~
    Sex + edu2 + alc4 + ever_smoked + PA + Updi_NoVegoils_z + GBMprot_Updi_z + Center_grp,
  data = df, weights = Weights, dist = "gengamma"
)

fit_t2d_gbmfull
# Updi_NoVegoils_z                    -2.28e-16   3.46e-02   1.83e-02   5.10e-02   8.36e-03   1.04e+00   1.02e+00   1.05e+00
# GBMfull_Updi_z                      -1.05e-17  -2.09e-02  -4.44e-02   2.72e-03   1.20e-02   9.79e-01   9.57e-01   1.00e+00

fit_t2d_gbmprot 
# Updi_NoVegoils_z                    -2.28e-16   3.45e-02   1.81e-02   5.09e-02   8.38e-03   1.04e+00   1.02e+00   1.05e+00
# GBMprot_Updi_z                      -2.04e-17  -1.58e-02  -3.42e-02   2.52e-03   9.36e-03   9.84e-01   9.66e-01   1.00e+00

fit_t2d_gbmfullonly <- flexsurvreg(
  Surv(Age_Recr, Age_Exit_T2D, T2D_Status) ~
    Sex + edu2 + alc4 + ever_smoked + PA + GBMfull_Updi_z + Center_grp,
  data = df, weights = Weights, dist = "gengamma"
)

fit_t2d_gbmprotonly <- flexsurvreg(
  Surv(Age_Recr, Age_Exit_T2D, T2D_Status) ~
    Sex + edu2 + alc4 + ever_smoked + PA + GBMprot_Updi_z + Center_grp,
  data = df, weights = Weights, dist = "gengamma"
)

fit_t2d_gbmfullonly 
#GBMfull_Updi_z                      -1.05e-17   2.03e-03  -1.91e-02   2.32e-02   1.08e-02   1.00e+00   9.81e-01   1.02e+00
fit_t2d_gbmprotonly 
#GBMprot_Updi_z                      -2.04e-17   1.43e-03  -1.50e-02   1.78e-02   8.37e-03   1.00e+00   9.85e-01   1.02e+00
