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
df$Enet_Hpdi_z       <- as.numeric(scale(df$EnetProtInt_Hpdi))
df$GBMfull_Hpdi_z    <- as.numeric(scale(df$GBMprotSHAP_Hpdi))
df$Hpdi_NoVegoils_z  <- as.numeric(scale(df$Hpdi_NoVegoils))

# z-scores: uPDI
df$Updi_NoVegoils_z  <- as.numeric(scale(df$Updi_NoVegoils))
df$Enet_Updi_z       <- as.numeric(scale(df$EnetProtInt_Updi))
df$GBM_Updi_z        <- as.numeric(scale(df$GBMprotSHAP_Updi))

#################
## hPDI models ##
#################

fit_ca_hpdi <- flexsurvreg(
  Surv(Age_Recr, Age_Exit_1_Cvd, Cvd_1_Status) ~
    Sex + edu2 + alc4 + ever_smoked + PA + Hpdi_NoVegoils_z + Center_grp,
  data = df, weights = Weights, dist = "gengamma"
)

fit_ca_enet <- flexsurvreg(
  Surv(Age_Recr, Age_Exit_1_Cvd, Cvd_1_Status) ~
    Sex + edu2 + alc4 + ever_smoked + PA + Hpdi_NoVegoils_z + Enet_Hpdi_z + Center_grp,
  data = df, weights = Weights, dist = "gengamma"
)

fit_ca_gbm_sparse <- flexsurvreg(
  Surv(Age_Recr, Age_Exit_1_Cvd, Cvd_1_Status) ~
    Sex + edu2 + alc4 + ever_smoked + PA + Hpdi_NoVegoils_z + GBMfull_Hpdi_z + Center_grp,
  data = df, weights = Weights, dist = "gengamma"
)

fit_ca_enetonly <- flexsurvreg(
  Surv(Age_Recr, Age_Exit_1_Cvd, Cvd_1_Status) ~
    Sex + edu2 + alc4 + ever_smoked + PA + Enet_Hpdi_z + Center_grp,
  data = df, weights = Weights, dist = "gengamma"
)

fit_ca_gbm_sparseonly <- flexsurvreg(
  Surv(Age_Recr, Age_Exit_1_Cvd, Cvd_1_Status) ~
    Sex + edu2 + alc4 + ever_smoked + PA + GBMfull_Hpdi_z + Center_grp,
  data = df, weights = Weights, dist = "gengamma"
)



fit_ca_limma <- flexsurvreg(
  Surv(Age_Recr, Age_Exit_1_Cvd, Cvd_1_Status) ~
    Sex + edu2 + alc4 + ever_smoked + PA + Hpdi_NoVegoils_z + hPDI_limma_unw_z + Center_grp,
  data = df, weights = Weights, dist = "gengamma")


fit_ca_limmaonly <- flexsurvreg(
  Surv(Age_Recr, Age_Exit_1_Cvd, Cvd_1_Status) ~
    Sex + edu2 + alc4 + ever_smoked + PA + hPDI_limma_unw_z + Center_grp,
  data = df, weights = Weights, dist = "gengamma")


#################
## uPDI models ##
#################

fit_ca_updi <- flexsurvreg(
  Surv(Age_Recr, Age_Exit_1_Cvd, Cvd_1_Status) ~
    Sex + edu2 + alc4 + ever_smoked + PA + Updi_NoVegoils_z + Center_grp,
  data = df, weights = Weights, dist = "gengamma"
)

fit_ca_enetupdi <- flexsurvreg(
  Surv(Age_Recr, Age_Exit_1_Cvd, Cvd_1_Status) ~
    Sex + edu2 + alc4 + ever_smoked + PA + Updi_NoVegoils_z + Enet_Updi_z + Center_grp,
  data = df, weights = Weights, dist = "gengamma"
)

fit_ca_gbm_updi <- flexsurvreg(
  Surv(Age_Recr, Age_Exit_1_Cvd, Cvd_1_Status) ~
    Sex + edu2 + alc4 + ever_smoked + PA + Updi_NoVegoils_z + GBM_Updi_z + Center_grp,
  data = df, weights = Weights, dist = "gengamma"
)

fit_ca_enetupdionly <- flexsurvreg(
  Surv(Age_Recr, Age_Exit_1_Cvd, Cvd_1_Status) ~
    Sex + edu2 + alc4 + ever_smoked + PA + Enet_Updi_z + Center_grp,
  data = df, weights = Weights, dist = "gengamma"
)

fit_ca_gbm_updionly <- flexsurvreg(
  Surv(Age_Recr, Age_Exit_1_Cvd, Cvd_1_Status) ~
    Sex + edu2 + alc4 + ever_smoked + PA + GBM_Updi_z + Center_grp,
  data = df, weights = Weights, dist = "gengamma"
)

fit_ca_limma_updi <- flexsurvreg(
  Surv(Age_Recr, Age_Exit_1_Cvd, Cvd_1_Status) ~
    Sex + edu2 + alc4 + ever_smoked + PA + Updi_NoVegoils_z + uPDI_limma_unw_z + Center_grp,
  data = df, weights = Weights, dist = "gengamma")

fit_ca_limmaupdionly <- flexsurvreg(
  Surv(Age_Recr, Age_Exit_1_Cvd, Cvd_1_Status) ~
    Sex + edu2 + alc4 + ever_smoked + PA + uPDI_limma_unw_z + Center_grp,
  data = df, weights = Weights, dist = "gengamma")

#################
##### PLOTS #####
#################
library(dplyr)
library(tibble)
library(ggplot2)
library(forcats)
library(patchwork)

extract_tr <- function(fit, var_name, divisor = 1) {
  tab <- fit$res
  row <- tab[rownames(tab) == var_name, , drop = FALSE]
  
  if (nrow(row) == 0) {
    stop(paste("Variable", var_name, "not found in model"))
  }
  
  est <- as.numeric(row[, "est"]) / divisor
  se  <- as.numeric(row[, "se"])  / divisor
  
  tibble(
    est     = est,
    se      = se,
    TR      = exp(est),
    TR_low  = exp(est - 1.96 * se),
    TR_high = exp(est + 1.96 * se)
  )
}

# -----------------------------
# LEFT PANEL: diet score effects
# -----------------------------
left_rows <- c(
  "hPDI only",
  "hPDI + limma score",
  "hPDI + ENET score",
  "hPDI + GBM score",
  "pad_mid",
  "uPDI only",
  "uPDI + limma score",
  "uPDI + ENET score",
  "uPDI + GBM score"
)

left_combined_df <- bind_rows(
  extract_tr(fit_ca_hpdi, "Hpdi_NoVegoils_z") %>%
    mutate(row = "hPDI only"),
  extract_tr(fit_ca_limma, "Hpdi_NoVegoils_z") %>%
    mutate(row = "hPDI + limma score"),
  extract_tr(fit_ca_enet, "Hpdi_NoVegoils_z") %>%
    mutate(row = "hPDI + ENET score"),
  extract_tr(fit_ca_gbm_sparse, "Hpdi_NoVegoils_z") %>%
    mutate(row = "hPDI + GBM score"),
  
  tibble(
    est = NA_real_, se = NA_real_, TR = NA_real_,
    TR_low = NA_real_, TR_high = NA_real_,
    row = "pad_mid"
  ),
  
  extract_tr(fit_ca_updi, "Updi_NoVegoils_z") %>%
    mutate(row = "uPDI only"),
  extract_tr(fit_ca_limma_updi, "Updi_NoVegoils_z") %>%
    mutate(row = "uPDI + limma score"),
  extract_tr(fit_ca_enetupdi, "Updi_NoVegoils_z") %>%
    mutate(row = "uPDI + ENET score"),
  extract_tr(fit_ca_gbm_updi, "Updi_NoVegoils_z") %>%
    mutate(row = "uPDI + GBM score")
) %>%
  mutate(row = factor(row, levels = rev(left_rows)))

p_left_combined <- ggplot(
  left_combined_df,
  aes(x = TR, y = row, xmin = TR_low, xmax = TR_high)
) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_pointrange(na.rm = TRUE) +
  scale_x_log10() +
  scale_y_discrete(labels = function(x) ifelse(grepl("^pad", x), "", x)) +
  labs(
    title = "Diet score estimates",
    x = "Time ratio for CVD (per 1 SD higher diet score)",
    y = NULL
  ) +
  theme_bw()

# -----------------------------------
# RIGHT PANEL: protein score effects
# -----------------------------------
right_rows <- c(
  "hPDI dummy",
  "hPDI limma score",
  "hPDI ENET score",
  "hPDI GBM score",
  "pad_mid",
  "uPDI dummy",
  "uPDI limma score",
  "uPDI ENET score",
  "uPDI GBM score"
)

right_combined_df <- bind_rows(
  tibble(
    est = c(NA_real_, NA_real_),
    se = c(NA_real_, NA_real_),
    TR = c(NA_real_, NA_real_),
    TR_low = c(NA_real_, NA_real_),
    TR_high = c(NA_real_, NA_real_),
    row = "hPDI dummy",
    model_type = c("Independent", "Jointly with hPDI")
  ),
  
  extract_tr(fit_ca_limmaonly, "hPDI_limma_unw_z") %>%
    mutate(row = "hPDI limma score", model_type = "Independent"),
  extract_tr(fit_ca_limma, "hPDI_limma_unw_z") %>%
    mutate(row = "hPDI limma score", model_type = "Jointly with hPDI"),
  
  extract_tr(fit_ca_enetonly, "Enet_Hpdi_z") %>%
    mutate(row = "hPDI ENET score", model_type = "Independent"),
  extract_tr(fit_ca_enet, "Enet_Hpdi_z") %>%
    mutate(row = "hPDI ENET score", model_type = "Jointly with hPDI"),
  
  extract_tr(fit_ca_gbm_sparseonly, "GBMfull_Hpdi_z") %>%
    mutate(row = "hPDI GBM score", model_type = "Independent"),
  extract_tr(fit_ca_gbm_sparse, "GBMfull_Hpdi_z") %>%
    mutate(row = "hPDI GBM score", model_type = "Jointly with hPDI"),
  
  tibble(
    est = c(NA_real_, NA_real_),
    se = c(NA_real_, NA_real_),
    TR = c(NA_real_, NA_real_),
    TR_low = c(NA_real_, NA_real_),
    TR_high = c(NA_real_, NA_real_),
    row = "pad_mid",
    model_type = c("Independent", "Jointly with hPDI")
  ),
  
  tibble(
    est = c(NA_real_, NA_real_),
    se = c(NA_real_, NA_real_),
    TR = c(NA_real_, NA_real_),
    TR_low = c(NA_real_, NA_real_),
    TR_high = c(NA_real_, NA_real_),
    row = "uPDI dummy",
    model_type = c("Independent", "Jointly with uPDI")
  ),
  
  extract_tr(fit_ca_limmaupdionly, "uPDI_limma_unw_z") %>%
    mutate(row = "uPDI limma score", model_type = "Independent"),
  extract_tr(fit_ca_limma_updi, "uPDI_limma_unw_z") %>%
    mutate(row = "uPDI limma score", model_type = "Jointly with uPDI"),
  
  extract_tr(fit_ca_enetupdionly, "Enet_Updi_z") %>%
    mutate(row = "uPDI ENET score", model_type = "Independent"),
  extract_tr(fit_ca_enetupdi, "Enet_Updi_z") %>%
    mutate(row = "uPDI ENET score", model_type = "Jointly with uPDI"),
  
  extract_tr(fit_ca_gbm_updionly, "GBM_Updi_z") %>%
    mutate(row = "uPDI GBM score", model_type = "Independent"),
  extract_tr(fit_ca_gbm_updi, "GBM_Updi_z") %>%
    mutate(row = "uPDI GBM score", model_type = "Jointly with uPDI")
) %>%
  mutate(
    row = factor(row, levels = rev(right_rows)),
    model_type = factor(
      model_type,
      levels = c("Independent", "Jointly with hPDI", "Jointly with uPDI")
    )
  )

pd <- position_dodge2(width = 0.5, reverse = TRUE)

p_right_combined <- ggplot(
  right_combined_df,
  aes(x = TR, y = row, xmin = TR_low, xmax = TR_high,
      colour = model_type, shape = model_type)
) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_pointrange(position = pd, na.rm = TRUE) +
  scale_x_log10() +
  scale_y_discrete(
    labels = c(
      "hPDI dummy" = "",
      "hPDI limma score" = "hPDI-limma",
      "hPDI ENET score" = "hPDI-ENET",
      "hPDI GBM score" = "hPDI-GBM",
      "pad_mid" = "",
      "uPDI dummy" = "",
      "uPDI limma score" = "uPDI-limma",
      "uPDI ENET score" = "uPDI-ENET",
      "uPDI GBM score" = "uPDI-GBM"
    )
  ) +
  scale_colour_manual(
    values = c(
      "Independent" = "#2B3256",
      "Jointly with hPDI" = "#D27E64",
      "Jointly with uPDI" = "#D27E64"
    )
  ) +
  labs(
    title = "Protein signature estimates",
    x = "Time ratio for CVD (per 1 unit higher protein signature)",
    y = NULL,
    colour = NULL,
    shape = NULL
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

p_final_2panel <- p_left_combined + p_right_combined +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

p_final_2panel

# save
ggsave(
  file.path(MANUSCRIPT_DIR, "CVD_survival_2panel_withlimma.png"),
  plot = p_final_2panel,
  width = 12, height = 7, units = "in", dpi = 600
)


####################
# Estimates to Excel
####################

library(writexl)
library(dplyr)
library(tibble)

extract_all_tr <- function(fit, model_name, outcome = "CVD") {
  tab <- as.data.frame(fit$res)
  tab$term <- rownames(fit$res)
  
  lcl_col <- names(tab)[grepl("^L95", names(tab))]
  ucl_col <- names(tab)[grepl("^U95", names(tab))]
  
  out <- tab %>%
    filter(
      !term %in% c("mu", "sigma", "Q", "shape", "scale",
                   "(Intercept)", "log(shape)", "log(scale)")
    ) %>%
    mutate(
      outcome = outcome,
      model = model_name,
      beta = est,
      SE = se,
      LCL_beta = .data[[lcl_col]],
      UCL_beta = .data[[ucl_col]],
      TR = exp(beta),
      TR_low = exp(LCL_beta),
      TR_high = exp(UCL_beta),
      p_value = 2 * pnorm(-abs(beta / SE))
    ) %>%
    dplyr::select(
      outcome, model, term, beta, SE, LCL_beta, UCL_beta,
      TR, TR_low, TR_high, p_value
    )
  
  out
}



all_results <- bind_rows(
  extract_all_tr(fit_ca_hpdi, "hPDI only"),
  extract_all_tr(fit_ca_limma, "hPDI + limma"),
  extract_all_tr(fit_ca_enet, "hPDI + Enet"),
  extract_all_tr(fit_ca_gbm_sparse, "hPDI + GBM"),
  extract_all_tr(fit_ca_limmaonly, "limma only (hPDI set)"),
  extract_all_tr(fit_ca_enetonly, "Enet only (hPDI set)"),
  extract_all_tr(fit_ca_gbm_sparseonly, "GBM only (hPDI set)"),
  
  extract_all_tr(fit_ca_updi, "uPDI only"),
  extract_all_tr(fit_ca_limma_updi, "uPDI + limma"),
  extract_all_tr(fit_ca_enetupdi, "uPDI + Enet"),
  extract_all_tr(fit_ca_gbm_updi, "uPDI + GBM"),
  extract_all_tr(fit_ca_limmaupdionly, "limma only (uPDI set)"),
  extract_all_tr(fit_ca_enetupdionly, "Enet only (uPDI set)"),
  extract_all_tr(fit_ca_gbm_updionly, "GBM only (uPDI set)")
)


# view in RStudio
if (interactive()) View(all_results)

# export to Excel
write_xlsx(
  all_results,
  file.path(MANUSCRIPT_DIR, "EPIC_CVD_survmodel_coefficients_limma_gengamma.xlsx")
)
