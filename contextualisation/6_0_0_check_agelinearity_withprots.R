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
# Script 1: Check age-linearity for prioritized proteins
# Correct direction: protein abundance ~ age + diet + covariates
# Tests age spline df 2–4 and selects df based on AIC
# ============================================================

library(dplyr)
library(tibble)
library(purrr)
library(tidyr)
library(splines)
library(openxlsx)

# ------------------------------------------------------------
# 0. Load data and prepare covariates
# ------------------------------------------------------------

load(file.path(DATA_DIR, "analyticaldf_postexcl_20251020.RData"))
wgts <- readRDS(file.path(DATA_DIR, "Weights_Feb2025.rds"))

wgts <- dplyr::rename(wgts, Idepic = idepic)

df1 <- df %>%
  left_join(wgts, by = "Idepic")

df <- df1 %>%
  filter(!is.na(Weights))

# Sex
df$Sex <- factor(df$Sex)
df$Sex <- relevel(df$Sex, ref = "Male")

# Physical activity
df$PA <- as.numeric(scale(df$PA_total_MET))
df$PA_total_MET_quint <- as.factor(df$PA_total_MET_quint)

# Alcohol, collapsed
df$alc4 <- NA
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
df$edu2 <- dplyr::case_when(
  df$L_School %in% c(0, 1) ~ "Below secondary",
  df$L_School %in% c(2, 3, 4) ~ "Secondary/vocational and above",
  df$L_School == 5 | is.na(df$L_School) ~ "Not specified",
  TRUE ~ NA_character_
)

df$edu2 <- factor(df$edu2)
df$edu2 <- relevel(df$edu2, ref = "Below secondary")

# Diet scores, standardized
df <- df %>%
  mutate(
    Hpdi_z = as.numeric(scale(Hpdi_NoVegoils)),
    Updi_z = as.numeric(scale(Updi_NoVegoils))
  )

# ------------------------------------------------------------
# 1. hPDI prioritized proteins
# ------------------------------------------------------------

protein_spec_hpdi <- tribble(
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
# 2. uPDI prioritized proteins
# ------------------------------------------------------------

protein_spec_updi <- tribble(
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
# 3. Function: test linear age vs spline age
# ------------------------------------------------------------

test_age_linearity_one_protein <- function(data,
                                           protein_col,
                                           diet_col,
                                           pdi_type_label,
                                           spline_dfs = 2:4) {
  
  covar_rhs_no_age <- paste(
    "Sex",
    "diet_value",
    "edu2",
    "alc4",
    "ever_smoked",
    "PA_total_MET_quint",
    "Center_grp",
    sep = " + "
  )
  
  vars_needed <- c(
    protein_col, diet_col, "Weights",
    "Sex", "Age_Recr", "edu2", "alc4", "ever_smoked",
    "PA_total_MET_quint", "Center_grp"
  )
  
  dat <- data %>%
    select(all_of(vars_needed)) %>%
    mutate(
      protein_value = .data[[protein_col]],
      diet_value = .data[[diet_col]],
      age_value = Age_Recr
    ) %>%
    filter(complete.cases(.))
  
  # Linear age model
  f_linear <- as.formula(
    paste0(
      "protein_value ~ age_value + ",
      covar_rhs_no_age
    )
  )
  
  m_linear <- lm(f_linear, data = dat, weights = Weights)
  
  linear_row <- tibble(
    pdi_type = pdi_type_label,
    protein_id = protein_col,
    model = "linear_age",
    age_spline_df = 1,
    n = nobs(m_linear),
    AIC = AIC(m_linear),
    BIC = BIC(m_linear),
    adj_r2 = summary(m_linear)$adj.r.squared,
    p_vs_linear = NA_real_
  )
  
  # Spline age models, df 2–4
  spline_rows <- map_dfr(spline_dfs, function(df_spline) {
    
    f_spline <- as.formula(
      paste0(
        "protein_value ~ ns(age_value, df = ", df_spline, ") + ",
        covar_rhs_no_age
      )
    )
    
    m_spline <- lm(f_spline, data = dat, weights = Weights)
    comp <- anova(m_linear, m_spline)
    
    tibble(
      pdi_type = pdi_type_label,
      protein_id = protein_col,
      model = paste0("age_spline_df", df_spline),
      age_spline_df = df_spline,
      n = nobs(m_spline),
      AIC = AIC(m_spline),
      BIC = BIC(m_spline),
      adj_r2 = summary(m_spline)$adj.r.squared,
      p_vs_linear = comp$`Pr(>F)`[2]
    )
  })
  
  bind_rows(linear_row, spline_rows) %>%
    mutate(
      delta_AIC_vs_linear = AIC - AIC[model == "linear_age"],
      delta_BIC_vs_linear = BIC - BIC[model == "linear_age"]
    )
}

# ------------------------------------------------------------
# 4. Run age-linearity tests for hPDI proteins
# ------------------------------------------------------------

age_linearity_hpdi <- protein_spec_hpdi %>%
  mutate(
    result = map(
      protein_id,
      ~ test_age_linearity_one_protein(
        data = df,
        protein_col = .x,
        diet_col = "Hpdi_z",
        pdi_type_label = "hPDI",
        spline_dfs = 2:4
      ) %>%
        select(-pdi_type, -protein_id)
    )
  ) %>%
  select(pdi_type, gene, protein_id, result) %>%
  unnest(result)

# ------------------------------------------------------------
# 5. Run age-linearity tests for uPDI proteins
# ------------------------------------------------------------

age_linearity_updi <- protein_spec_updi %>%
  mutate(
    result = map(
      protein_id,
      ~ test_age_linearity_one_protein(
        data = df,
        protein_col = .x,
        diet_col = "Updi_z",
        pdi_type_label = "uPDI",
        spline_dfs = 2:4
      ) %>%
        select(-pdi_type, -protein_id)
    )
  ) %>%
  select(pdi_type, gene, protein_id, result) %>%
  unnest(result)

# ------------------------------------------------------------
# 6. Combine all model results
# ------------------------------------------------------------

age_linearity_all <- bind_rows(
  age_linearity_hpdi,
  age_linearity_updi
)

# ------------------------------------------------------------
# 7. Summary table
# ------------------------------------------------------------

age_linearity_summary <- age_linearity_all %>%
  group_by(pdi_type, gene, protein_id) %>%
  summarise(
    n = first(n),
    
    best_df_AIC = age_spline_df[which.min(AIC)],
    best_model_AIC = model[which.min(AIC)],
    best_AIC = min(AIC),
    
    best_df_BIC = age_spline_df[which.min(BIC)],
    best_model_BIC = model[which.min(BIC)],
    best_BIC = min(BIC),
    
    linear_AIC = AIC[model == "linear_age"],
    linear_BIC = BIC[model == "linear_age"],
    
    min_p_vs_linear = min(p_vs_linear, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    p_age_nonlinearity_fdr = p.adjust(min_p_vs_linear, method = "BH"),
    
    evidence_age_nonlinearity = case_when(
      p_age_nonlinearity_fdr < 0.05 ~ "Strong evidence of age non-linearity",
      min_p_vs_linear < 0.05 ~ "Nominal evidence of age non-linearity",
      TRUE ~ "No clear evidence of age non-linearity"
    ),
    
    selected_age_df = case_when(
      evidence_age_nonlinearity == "Strong evidence of age non-linearity" &
        best_model_AIC != "linear_age" ~ best_df_AIC,
      
      TRUE ~ 1
    ),
    
    age_df_selection_reason = case_when(
      evidence_age_nonlinearity == "Strong evidence of age non-linearity" &
        best_model_AIC != "linear_age" ~ "Selected by AIC among age spline df 2–4",
      
      TRUE ~ "Linear age retained"
    )
  ) %>%
  arrange(pdi_type, min_p_vs_linear)

if (interactive()) View(age_linearity_summary)

# ------------------------------------------------------------
# 8. Save
# ------------------------------------------------------------

write.xlsx(
  list(
    hPDI_all_models = age_linearity_hpdi,
    uPDI_all_models = age_linearity_updi,
    all_models = age_linearity_all,
    summary = age_linearity_summary
  ),
  file = file.path(MANUSCRIPT_DIR, "age_linearity_prioritised_proteins_correct_direction_AIC_df2to4.xlsx"),
  overwrite = TRUE
)
