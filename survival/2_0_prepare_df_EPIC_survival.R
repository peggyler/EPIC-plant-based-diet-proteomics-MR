# Run from the repository root in a fresh R session. See README.md.
# EDIT LOCAL PATHS HERE. No study data or results are supplied.
PROJECT_ROOT <- "."
DATA_DIR <- file.path(PROJECT_ROOT, "data")
OUTPUT_ROOT <- file.path(PROJECT_ROOT, "outputs")
MANUSCRIPT_DIR <- file.path(OUTPUT_ROOT, "manuscript")
REFERENCE_DIR <- file.path(DATA_DIR, "reference")
SIGNATURE_SCRIPTS <- file.path(PROJECT_ROOT, "signatures")
CONTEXT_SCRIPTS <- file.path(PROJECT_ROOT, "contextualisation")
HPDI_RUN <- "staged_run_20260424_0935"  # choose the manuscript-producing run/date
UPDI_RUN <- "staged_run_20260424_1844"  # choose the manuscript-producing run/date
dir.create(MANUSCRIPT_DIR, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(dplyr))



load(file.path(DATA_DIR, "score_tte_analysis_20251126.RData"))
wgts <- readRDS(file.path(DATA_DIR, "Weights_Feb2025.rds"))

df <- df1

wgts <- dplyr::rename(wgts, Idepic = idepic)

df <- df %>%
  left_join(wgts, by = "Idepic") 

summary(df$Weights)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   1.000   1.000   1.000   5.468   9.762  80.680     909 
nrow(df) #14882
df <- df %>%
  filter(!is.na(Weights))

#13973

df$Sex <- factor(df$Sex)
df$Sex <- relevel(df$Sex, ref = "Male")

# quintile cut points
# q <- quantile(df$PA_total_MET, probs = seq(0, 1, 0.2), na.rm = TRUE)
# 
# # create quintile variable
# df$PA <- cut(
#   df$PA_total_MET,
#   breaks = q,
#   include.lowest = TRUE,
#   labels = c("Q1", "Q2", "Q3", "Q4", "Q5")
# )

df$PA <- as.numeric(scale(df$PA_total_MET))

# Alcohol
df$alc4 <- NA
df$alc4[df$alc_freq == "Never"] <- "Never"
df$alc4[df$alc_freq == "Special occasions only"] <- "Occasional"
df$alc4[df$alc_freq %in% c("1–3 times/month", "1–2 times/week")] <- "Moderate"
df$alc4[df$alc_freq %in% c("3–4 times/week", "Daily or almost daily")] <- "Regular/Daily"

df$alc4 <- factor(df$alc4,
                  levels = c("Moderate", "Never", "Occasional", "Regular/Daily"))

df$ever_smoked <- factor(
  ifelse(df$smoke == "Never", "Never smoked", "Ever Smoked"),
  levels = c("Never smoked", "Ever Smoked")
)

# z-scores: hPDI
# df$Enet_Hpdi_z       <- as.numeric(scale(df$Enet_Hpdi))
# df$GBMfull_Hpdi_z    <- as.numeric(scale(df$GBMfull_Hpdi))
# df$GBM_Hpdi_z        <- as.numeric(scale(df$GBM_Hpdi))
# df$Hpdi_NoVegoils_z  <- as.numeric(scale(df$Hpdi_NoVegoils))
# 
# # z-scores: uPDI
# df$Updi_NoVegoils_z  <- as.numeric(scale(df$Updi_NoVegoils))
# df$Enet_Updi_z       <- as.numeric(scale(df$Enet_Updi))
# df$GBM_Updi_z        <- as.numeric(scale(df$GBM_Updi))

# Center groups
df <- df %>%
  mutate(
    Center_num = as.numeric(as.character(Center)),
    Center_grp = case_when(
      Center_num %in% c(21,22,23,24,25) ~ "Italy",
      Center_num %in% c(31,32,33,34,35) ~ "Spain",
      Center_num %in% c(41,42)          ~ "United Kingdom",
      Center_num %in% c(51,52)          ~ "Netherlands",
      TRUE ~ NA_character_
    )
  )

df$Center_grp <- factor(df$Center_grp)
df$Center_grp <- relevel(df$Center_grp, ref = "Spain")


df$edu2 <- dplyr::case_when(
  df$L_School %in% c(0, 1) ~ "Below secondary",
  df$L_School %in% c(2, 3, 4) ~ "Secondary/vocational and above",
  df$L_School == 5 | is.na(df$L_School) ~ "Not specified",
  TRUE ~ NA_character_
)

table(df$edu2, useNA = "ifany")
# Below secondary                  Not specified Secondary/vocational and above 
#7229                            464                           6280 
df$edu2 <- factor(df$edu2)
df$edu2 <- relevel(df$edu2, ref = "Below secondary")

# Cancer variables
df <- df %>%
  mutate(
    moved = is.na(Age_Exit_Cancer_1st) & !is.na(Age_Evt_Cancer_2nd),
    Age_Exit_Cancer_1st = if_else(moved, Age_Exit_Cancer_2nd, Age_Exit_Cancer_1st),
    Age_Evt_Cancer_1st  = if_else(moved, Age_Evt_Cancer_2nd, Age_Evt_Cancer_1st),
    Age_Evt_Cancer_2nd  = if_else(moved, NA_real_, Age_Evt_Cancer_2nd),
    Age_Exit_Cancer_1st = if_else(
      is.na(Age_Exit_Cancer_1st) & !is.na(Age_Exit_Cancer_2nd),
      Age_Exit_Cancer_2nd,
      Age_Exit_Cancer_1st
    )
  )

df$Cancer_Status <- ifelse(!is.na(df$Age_Evt_Cancer_1st), 1, 0)

##-----------
## remove vars
##-----------
vars_to_remove <- c(
  "Enet_Hpdi",
  "Enet60_Hpdi",
  "GBM_sparse_Hpdi",
  "GBM_Hpdi",
  "EnetGBM_Hpdi",
  "EnetGBM_all_Hpdi",
  "GBMfull_sparse_Hpdi",
  "GBMfull_Hpdi",
  "Enet_Updi",
  "uPDI_limma_unw",
  "uPDI_limma_unw_z",
  "uPDI_limma",
  "uPDI_limma_z",
  "hPDI_limma_unw",
  "hPDI_limma_unw_z",
  "hPDI_limma",
  "hPDI_limma_z",
  "GBM_Updi"
)

df_sel <- df %>%
  select(-any_of(vars_to_remove))

DDIR <- OUTPUT_ROOT

limma <- "limma_scores.csv"
gbmu <- file.path("updi", "gbm", UPDI_RUN, "scores_external_updi.csv")
gbmh <- file.path("hpdi", "gbm", HPDI_RUN, "scores_external.csv")
eneth <- "hpdi/scores_external_enet.csv"
enetu <-"updi/scores_external_enet.csv"

## merge
library(dplyr)
library(readr)

# Read files
limma_df <- read_csv(file.path(DDIR, limma), show_col_types = FALSE)
gbmu_df  <- read_csv(file.path(DDIR, gbmu),  show_col_types = FALSE)
gbmh_df  <- read_csv(file.path(DDIR, gbmh),  show_col_types = FALSE)
eneth_df <- read_csv(file.path(DDIR, eneth), show_col_types = FALSE)
enetu_df <- read_csv(file.path(DDIR, enetu), show_col_types = FALSE)

# Check ID column exists
score_dfs <- list(
  limma = limma_df,
  gbmu  = gbmu_df,
  gbmh  = gbmh_df,
  eneth = eneth_df,
  enetu = enetu_df
)

lapply(score_dfs, function(x) stopifnot("Idepic" %in% names(x)))

# Merge with df_sel
df_sel <- df_sel %>%
  left_join(limma_df, by = "Idepic") %>%
  left_join(gbmu_df,  by = "Idepic") %>%
  left_join(gbmh_df,  by = "Idepic") %>%
  left_join(eneth_df, by = "Idepic") %>%
  left_join(enetu_df, by = "Idepic")

df_sel <- df_sel %>%
  select(
    -any_of(c(
      "Hpdi_NoVegoils_raw.x",
      "Hpdi_NoVegoils_raw.y",
      "Updi_NoVegoils_raw.x",
      "Updi_NoVegoils_raw.y"
    ))
  )

save(
  df_sel,
  file = file.path(DATA_DIR, "df_with_scores_surv_20260425.RData")
)
