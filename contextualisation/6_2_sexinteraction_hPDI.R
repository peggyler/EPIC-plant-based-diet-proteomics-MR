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

# ============================================================
# Script: Sex interaction for hPDI-prioritised proteins
# Correct direction:
# protein abundance ~ hPDI_z * Sex + age term + covariates
#
# This script starts from scratch:
# - loads EPIC analytic data
# - loads weights
# - prepares covariates
# - reads selected age df from saved age-linearity file
# - reads selected hPDI df from saved hPDI-linearity file
# - tests sex interaction for each prioritised hPDI protein
# ============================================================

library(dplyr)
library(tibble)
library(purrr)
library(tidyr)
library(splines)
library(openxlsx)

# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

data_file <- file.path(DATA_DIR, "analyticaldf_postexcl_20251020.RData")
weights_file <- file.path(DATA_DIR, "Weights_Feb2025.rds")

out_dir <- MANUSCRIPT_DIR

age_linearity_file <- file.path(
  out_dir,
  "age_linearity_prioritised_proteins_correct_direction_AIC_df2to4.xlsx"
)

hpdi_linearity_file <- file.path(
  out_dir,
  "linearity_hpdi_correct_direction_ageDFplugged.xlsx"
)

output_file <- file.path(
  out_dir,
  "sex_interaction_hpdi_correct_direction_ageHpdiDFplugged.xlsx"
)

# ------------------------------------------------------------
# 1. Load data and weights
# ------------------------------------------------------------

load(data_file)

wgts <- readRDS(weights_file) %>%
  rename(Idepic = idepic)

df <- df %>%
  left_join(wgts, by = "Idepic") %>%
  filter(!is.na(Weights))

# ------------------------------------------------------------
# 2. Prepare covariates
# ------------------------------------------------------------

# Sex
df$Sex <- factor(df$Sex)
df$Sex <- relevel(df$Sex, ref = "Male")

# Physical activity
df$PA_total_MET_quint <- as.factor(df$PA_total_MET_quint)

# Alcohol collapsed into alc4
df$alc4 <- NA_character_

df$alc4[df$alc_freq == "Never"] <- "Never"
df$alc4[df$alc_freq == "Special occasions only"] <- "Occasional"
df$alc4[df$alc_freq %in% c("1–3 times/month", "1–2 times/week")] <- "Moderate"
df$alc4[df$alc_freq %in% c("3–4 times/week", "Daily or almost daily")] <- "Regular/Daily"

df$alc4 <- factor(
  df$alc4,
  levels = c("Moderate", "Never", "Occasional", "Regular/Daily")
)

# Smoking, ever vs never
df$ever_smoked <- factor(
  ifelse(df$smoke == "Never", "Never smoked", "Ever Smoked"),
  levels = c("Never smoked", "Ever Smoked")
)

# Center group
df <- df %>%
  mutate(
    Center_num = as.numeric(as.character(Center)),
    Center_grp = case_when(
      Center_num %in% c(21, 22, 23, 24, 25) ~ "Italy",
      Center_num %in% c(31, 32, 33, 34, 35) ~ "Spain",
      Center_num %in% c(41, 42)             ~ "United Kingdom",
      Center_num %in% c(51, 52)             ~ "Netherlands",
      TRUE ~ NA_character_
    )
  )

df$Center_grp <- factor(df$Center_grp)
df$Center_grp <- relevel(df$Center_grp, ref = "Spain")

# Education
df$edu2 <- case_when(
  df$L_School %in% c(0, 1) ~ "Below secondary",
  df$L_School %in% c(2, 3, 4) ~ "Secondary/vocational and above",
  df$L_School == 5 | is.na(df$L_School) ~ "Not specified",
  TRUE ~ NA_character_
)

df$edu2 <- factor(df$edu2)
df$edu2 <- relevel(df$edu2, ref = "Below secondary")

# Standardised hPDI
df <- df %>%
  mutate(
    Hpdi_z = as.numeric(scale(Hpdi_NoVegoils))
  )

# ------------------------------------------------------------
# 3. hPDI-prioritised proteins
# ------------------------------------------------------------

protein_spec_hpdi_base <- tribble(
  ~gene,    ~protein_id,
  "CTSS",   "seq.3181.50",
  "DNPH1",  "seq.17785.11",
  "PTK7",   "seq.9525.1",
  "GREM1",  "seq.18878.15",
  "CCM2",   "seq.12347.29",
  "POSTN",  "seq.6645.53",
  "CA11",   "seq.8859.51",
  "MMP10",  "seq.10479.18",
  "MMP10",  "seq.8479.4",
  "GSTO1",  "seq.12436.84",
  "LRIG1",  "seq.18831.6",
  "MANBA",  "seq.6382.17",
  "NCAN",   "seq.15573.110",
  "PARVA",  "seq.13434.172",
  "PCSK9",  "seq.5231.79",
  "PCSK9",  "seq.8246.9",
  "PEAR1",  "seq.8275.31",
  "SPP1",   "seq.13113.7"
) %>%
  mutate(pdi_type = "hPDI") %>%
  filter(protein_id %in% names(df))

# ------------------------------------------------------------
# 4. Read selected df from saved files
# ------------------------------------------------------------

age_linearity_summary <- read.xlsx(
  age_linearity_file,
  sheet = "summary"
)

hpdi_linearity_summary <- read.xlsx(
  hpdi_linearity_file,
  sheet = "summary"
)

# Check required columns
required_age_cols <- c(
  "pdi_type",
  "protein_id",
  "selected_age_df",
  "evidence_age_nonlinearity",
  "age_df_selection_reason"
)

required_hpdi_cols <- c(
  "protein_id",
  "selected_hpdi_df",
  "evidence_hpdi_nonlinearity",
  "hpdi_df_selection_reason"
)

missing_age_cols <- setdiff(required_age_cols, names(age_linearity_summary))
missing_hpdi_cols <- setdiff(required_hpdi_cols, names(hpdi_linearity_summary))

if (length(missing_age_cols) > 0) {
  stop("Missing columns in age_linearity_summary: ", paste(missing_age_cols, collapse = ", "))
}

if (length(missing_hpdi_cols) > 0) {
  stop("Missing columns in hpdi_linearity_summary: ", paste(missing_hpdi_cols, collapse = ", "))
}

# ------------------------------------------------------------
# 5. Combine protein list with selected age df and hPDI df
# ------------------------------------------------------------

protein_spec_hpdi_sex <- protein_spec_hpdi_base %>%
  left_join(
    age_linearity_summary %>%
      filter(pdi_type == "hPDI") %>%
      select(
        protein_id,
        selected_age_df,
        evidence_age_nonlinearity,
        age_df_selection_reason
      ),
    by = "protein_id"
  ) %>%
  left_join(
    hpdi_linearity_summary %>%
      select(
        protein_id,
        selected_hpdi_df,
        evidence_hpdi_nonlinearity,
        hpdi_df_selection_reason
      ),
    by = "protein_id"
  ) %>%
  mutate(
    selected_age_df = as.numeric(selected_age_df),
    selected_hpdi_df = as.numeric(selected_hpdi_df),
    selected_age_df = coalesce(selected_age_df, 1),
    selected_hpdi_df = coalesce(selected_hpdi_df, 1),
    use_age_spline = selected_age_df > 1,
    use_hpdi_spline = selected_hpdi_df > 1
  )

if (interactive()) View(protein_spec_hpdi_sex)

# ------------------------------------------------------------
# 6. Function: test sex interaction
# ------------------------------------------------------------

test_hpdi_sex_interaction_one_protein <- function(data,
                                                  protein_col,
                                                  selected_age_df = 1,
                                                  selected_hpdi_df = 1,
                                                  diet_col = "Hpdi_z") {
  
  # Age term selected from age-linearity model
  age_term <- if (selected_age_df > 1) {
    paste0("ns(age_value, df = ", selected_age_df, ")")
  } else {
    "age_value"
  }
  
  # hPDI term selected from hPDI-linearity model
  hpdi_term <- if (selected_hpdi_df > 1) {
    paste0("ns(diet_value, df = ", selected_hpdi_df, ")")
  } else {
    "diet_value"
  }
  
  # Covariates excluding Sex because Sex is part of the interaction test
  covar_rhs_no_sex <- paste(
    age_term,
    "edu2",
    "alc4",
    "ever_smoked",
    "PA_total_MET_quint",
    "Center_grp",
    sep = " + "
  )
  
  vars_needed <- c(
    protein_col,
    diet_col,
    "Weights",
    "Sex",
    "Age_Recr",
    "edu2",
    "alc4",
    "ever_smoked",
    "PA_total_MET_quint",
    "Center_grp"
  )
  
  dat <- data %>%
    select(all_of(vars_needed)) %>%
    mutate(
      protein_value = .data[[protein_col]],
      diet_value = .data[[diet_col]],
      age_value = Age_Recr
    ) %>%
    filter(complete.cases(.))
  
  # Main-effects model
  f_main <- as.formula(
    paste0(
      "protein_value ~ Sex + ",
      hpdi_term, " + ",
      covar_rhs_no_sex
    )
  )
  
  # Sex-interaction model
  f_int <- as.formula(
    paste0(
      "protein_value ~ ",
      hpdi_term, " * Sex + ",
      covar_rhs_no_sex
    )
  )
  
  m_main <- lm(f_main, data = dat, weights = Weights)
  m_int  <- lm(f_int,  data = dat, weights = Weights)
  
  comp <- anova(m_main, m_int)
  
  tibble(
    protein_id = protein_col,
    selected_age_df = selected_age_df,
    selected_hpdi_df = selected_hpdi_df,
    age_form = if_else(
      selected_age_df > 1,
      paste0("Age spline df=", selected_age_df),
      "Linear age"
    ),
    hpdi_form = if_else(
      selected_hpdi_df > 1,
      paste0("hPDI spline df=", selected_hpdi_df),
      "Linear hPDI"
    ),
    n = nobs(m_int),
    AIC_main = AIC(m_main),
    AIC_interaction = AIC(m_int),
    BIC_main = BIC(m_main),
    BIC_interaction = BIC(m_int),
    delta_AIC_interaction = AIC(m_int) - AIC(m_main),
    delta_BIC_interaction = BIC(m_int) - BIC(m_main),
    p_interaction_sex = comp$`Pr(>F)`[2]
  )
}

# ------------------------------------------------------------
# 7. Run sex interaction models
# ------------------------------------------------------------

sex_interaction_hpdi <- protein_spec_hpdi_sex %>%
  mutate(
    result = pmap(
      list(protein_id, selected_age_df, selected_hpdi_df),
      ~ test_hpdi_sex_interaction_one_protein(
        data = df,
        protein_col = ..1,
        selected_age_df = ..2,
        selected_hpdi_df = ..3,
        diet_col = "Hpdi_z"
      ) %>%
        select(-protein_id, -selected_age_df, -selected_hpdi_df)
    )
  ) %>%
  select(
    pdi_type,
    gene,
    protein_id,
    selected_age_df,
    evidence_age_nonlinearity,
    age_df_selection_reason,
    selected_hpdi_df,
    evidence_hpdi_nonlinearity,
    hpdi_df_selection_reason,
    result
  ) %>%
  unnest(result) %>%
  mutate(
    p_interaction_sex_fdr = p.adjust(p_interaction_sex, method = "BH"),
    evidence_sex_interaction = case_when(
      p_interaction_sex_fdr < 0.05 ~ "Strong evidence of sex interaction",
      p_interaction_sex < 0.05 ~ "Nominal evidence of sex interaction",
      TRUE ~ "No clear evidence of sex interaction"
    )
  ) %>%
  arrange(p_interaction_sex)

if (interactive()) View(sex_interaction_hpdi)

# ------------------------------------------------------------
# 8. Save output
# ------------------------------------------------------------

write.xlsx(
  list(
    protein_spec = protein_spec_hpdi_sex,
    sex_interaction = sex_interaction_hpdi
  ),
  file = output_file,
  overwrite = TRUE
)
