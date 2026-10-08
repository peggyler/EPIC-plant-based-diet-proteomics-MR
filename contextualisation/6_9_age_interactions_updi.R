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
# Script: Age interaction for uPDI-prioritised proteins
# Correct direction:
# protein abundance ~ uPDI_z * age + covariates
#
# Uses:
# - selected age df from age-linearity file
# - selected uPDI df from uPDI-linearity file
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

updi_linearity_file <- file.path(
  out_dir,
  "linearity_updi_correct_direction_ageDFplugged.xlsx"
)

output_file <- file.path(
  out_dir,
  "age_interaction_updi_correct_direction_ageUpdiDFplugged.xlsx"
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

df$Sex <- factor(df$Sex)
df$Sex <- relevel(df$Sex, ref = "Male")

df$PA_total_MET_quint <- as.factor(df$PA_total_MET_quint)

df$alc4 <- NA_character_
df$alc4[df$alc_freq == "Never"] <- "Never"
df$alc4[df$alc_freq == "Special occasions only"] <- "Occasional"
df$alc4[df$alc_freq %in% c("1–3 times/month", "1–2 times/week")] <- "Moderate"
df$alc4[df$alc_freq %in% c("3–4 times/week", "Daily or almost daily")] <- "Regular/Daily"

df$alc4 <- factor(
  df$alc4,
  levels = c("Moderate", "Never", "Occasional", "Regular/Daily")
)

df$ever_smoked <- factor(
  ifelse(df$smoke == "Never", "Never smoked", "Ever Smoked"),
  levels = c("Never smoked", "Ever Smoked")
)

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

df$edu2 <- case_when(
  df$L_School %in% c(0, 1) ~ "Below secondary",
  df$L_School %in% c(2, 3, 4) ~ "Secondary/vocational and above",
  df$L_School == 5 | is.na(df$L_School) ~ "Not specified",
  TRUE ~ NA_character_
)

df$edu2 <- factor(df$edu2)
df$edu2 <- relevel(df$edu2, ref = "Below secondary")

df <- df %>%
  mutate(
    Updi_z = as.numeric(scale(Updi_NoVegoils))
  )

# ------------------------------------------------------------
# 3. uPDI-prioritised proteins
# ------------------------------------------------------------

protein_spec_updi_base <- tribble(
  ~gene,       ~protein_id,
  "KRT5",      "seq.11177.16",
  "EGFR",      "seq.2677.1",
  "CPNE1",     "seq.11602.12",
  "PLA2G2A",   "seq.2692.74",
  "SERPINA3",  "seq.2879.9",
  "HHIP",      "seq.10833.64",
  "CCM2",      "seq.12347.29",
  "ACE",       "seq.10714.7",
  "LCT",       "seq.9017.58",
  "PTPN9",     "seq.12633.3",
  "POSTN",     "seq.3457.57",
  "POSTN",     "seq.6645.53",
  "MMP10",     "seq.10479.18",
  "MMP10",     "seq.8479.4",
  "PCSK9",     "seq.8246.9",
  "ACADVL",    "seq.19262.219",
  "NADK",      "seq.13624.17",
  "SPP1",      "seq.13113.7"
) %>%
  mutate(pdi_type = "uPDI") %>%
  filter(protein_id %in% names(df))

# ------------------------------------------------------------
# 4. Read selected df from saved files
# ------------------------------------------------------------

age_linearity_summary <- read.xlsx(
  age_linearity_file,
  sheet = "summary"
)

updi_linearity_summary <- read.xlsx(
  updi_linearity_file,
  sheet = "summary"
)

required_age_cols <- c(
  "pdi_type",
  "protein_id",
  "selected_age_df",
  "evidence_age_nonlinearity",
  "age_df_selection_reason"
)

required_updi_cols <- c(
  "protein_id",
  "selected_updi_df",
  "evidence_updi_nonlinearity",
  "updi_df_selection_reason"
)

missing_age_cols <- setdiff(required_age_cols, names(age_linearity_summary))
missing_updi_cols <- setdiff(required_updi_cols, names(updi_linearity_summary))

if (length(missing_age_cols) > 0) {
  stop("Missing columns in age_linearity_summary: ", paste(missing_age_cols, collapse = ", "))
}

if (length(missing_updi_cols) > 0) {
  stop("Missing columns in updi_linearity_summary: ", paste(missing_updi_cols, collapse = ", "))
}

# ------------------------------------------------------------
# 5. Combine protein list with selected age df and uPDI df
# ------------------------------------------------------------

protein_spec_updi_age <- protein_spec_updi_base %>%
  left_join(
    age_linearity_summary %>%
      filter(pdi_type == "uPDI") %>%
      select(
        protein_id,
        selected_age_df,
        evidence_age_nonlinearity,
        age_df_selection_reason
      ),
    by = "protein_id"
  ) %>%
  left_join(
    updi_linearity_summary %>%
      select(
        protein_id,
        selected_updi_df,
        evidence_updi_nonlinearity,
        updi_df_selection_reason
      ),
    by = "protein_id"
  ) %>%
  mutate(
    selected_age_df = as.numeric(selected_age_df),
    selected_updi_df = as.numeric(selected_updi_df),
    selected_age_df = coalesce(selected_age_df, 1),
    selected_updi_df = coalesce(selected_updi_df, 1),
    use_age_spline = selected_age_df > 1,
    use_updi_spline = selected_updi_df > 1
  )

if (interactive()) View(protein_spec_updi_age)

# ------------------------------------------------------------
# 6. Function: test uPDI x age interaction
# ------------------------------------------------------------

test_updi_age_interaction_one_protein <- function(data,
                                                  protein_col,
                                                  selected_age_df = 1,
                                                  selected_updi_df = 1,
                                                  diet_col = "Updi_z") {
  
  age_term <- if (selected_age_df > 1) {
    paste0("ns(age_value, df = ", selected_age_df, ")")
  } else {
    "age_value"
  }
  
  updi_term <- if (selected_updi_df > 1) {
    paste0("ns(diet_value, df = ", selected_updi_df, ")")
  } else {
    "diet_value"
  }
  
  covar_rhs_no_age_diet <- paste(
    "Sex",
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
    "Age_Recr",
    "Sex",
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
      "protein_value ~ ",
      updi_term, " + ",
      age_term, " + ",
      covar_rhs_no_age_diet
    )
  )
  
  # uPDI x age interaction model
  f_int <- as.formula(
    paste0(
      "protein_value ~ ",
      updi_term, " * ", age_term, " + ",
      covar_rhs_no_age_diet
    )
  )
  
  m_main <- lm(f_main, data = dat, weights = Weights)
  m_int  <- lm(f_int,  data = dat, weights = Weights)
  
  comp <- anova(m_main, m_int)
  
  tibble(
    protein_id = protein_col,
    selected_age_df = selected_age_df,
    selected_updi_df = selected_updi_df,
    age_form = if_else(
      selected_age_df > 1,
      paste0("Age spline df=", selected_age_df),
      "Linear age"
    ),
    updi_form = if_else(
      selected_updi_df > 1,
      paste0("uPDI spline df=", selected_updi_df),
      "Linear uPDI"
    ),
    n = nobs(m_int),
    AIC_main = AIC(m_main),
    AIC_interaction = AIC(m_int),
    BIC_main = BIC(m_main),
    BIC_interaction = BIC(m_int),
    delta_AIC_interaction = AIC(m_int) - AIC(m_main),
    delta_BIC_interaction = BIC(m_int) - BIC(m_main),
    p_interaction_age = comp$`Pr(>F)`[2]
  )
}

# ------------------------------------------------------------
# 7. Run age interaction models
# ------------------------------------------------------------

age_interaction_updi <- protein_spec_updi_age %>%
  mutate(
    result = pmap(
      list(protein_id, selected_age_df, selected_updi_df),
      ~ test_updi_age_interaction_one_protein(
        data = df,
        protein_col = ..1,
        selected_age_df = ..2,
        selected_updi_df = ..3,
        diet_col = "Updi_z"
      ) %>%
        select(-protein_id, -selected_age_df, -selected_updi_df)
    )
  ) %>%
  select(
    pdi_type,
    gene,
    protein_id,
    selected_age_df,
    evidence_age_nonlinearity,
    age_df_selection_reason,
    selected_updi_df,
    evidence_updi_nonlinearity,
    updi_df_selection_reason,
    result
  ) %>%
  unnest(result) %>%
  mutate(
    p_interaction_age_fdr = p.adjust(p_interaction_age, method = "BH"),
    evidence_age_interaction = case_when(
      p_interaction_age_fdr < 0.05 ~ "Strong evidence of age interaction",
      p_interaction_age < 0.05 ~ "Nominal evidence of age interaction",
      TRUE ~ "No clear evidence of age interaction"
    )
  ) %>%
  arrange(p_interaction_age)

if (interactive()) View(age_interaction_updi)

# ------------------------------------------------------------
# 8. Save output
# ------------------------------------------------------------

write.xlsx(
  list(
    protein_spec = protein_spec_updi_age,
    age_interaction = age_interaction_updi
  ),
  file = output_file,
  overwrite = TRUE
)
