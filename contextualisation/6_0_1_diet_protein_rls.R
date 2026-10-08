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
# Formal diet non-linearity testing for prioritized proteins
#
# Correct direction:
# protein abundance ~ dietary index + age + covariates
#
# Tests:
# M0: protein ~ covariates
# M1: protein ~ diet + covariates
# M2: protein ~ ns(diet, df = 3) + covariates
#
# Main questions:
# 1. Is there a linear diet association?        M0 vs M1
# 2. Is there any flexible diet association?    M0 vs M2
# 3. Is there non-linearity beyond linear?      M1 vs M2
#
# Uses selected age df from previous age-linearity results.
# ============================================================

library(dplyr)
library(tibble)
library(purrr)
library(tidyr)
library(splines)
library(openxlsx)
library(broom)

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

out_file <- file.path(
  out_dir,
  "diet_nonlinearity_formal_test_fixed_df3_ageDFplugged.xlsx"
)

# ------------------------------------------------------------
# 1. Load data and prepare covariates
# ------------------------------------------------------------

load(data_file)

wgts <- readRDS(weights_file) %>%
  rename(Idepic = idepic)

df <- df %>%
  left_join(wgts, by = "Idepic") %>%
  filter(!is.na(Weights))

# Sex
df$Sex <- factor(df$Sex)
df$Sex <- relevel(df$Sex, ref = "Male")

# Physical activity
df$PA <- as.numeric(scale(df$PA_total_MET))
df$PA_total_MET_quint <- as.factor(df$PA_total_MET_quint)

# Alcohol, collapsed
df$alc4 <- NA_character_
df$alc4[df$alc_freq == "Never"] <- "Never"
df$alc4[df$alc_freq == "Special occasions only"] <- "Occasional"
df$alc4[df$alc_freq %in% c("1–3 times/month", "1–2 times/week")] <- "Moderate"
df$alc4[df$alc_freq %in% c("3–4 times/week", "Daily or almost daily")] <- "Regular/Daily"

df$alc4 <- factor(
  df$alc4,
  levels = c("Moderate", "Never", "Occasional", "Regular/Daily")
)

# Smoking
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

# Diet scores, standardized
df <- df %>%
  mutate(
    Hpdi_z = as.numeric(scale(Hpdi_NoVegoils)),
    Updi_z = as.numeric(scale(Updi_NoVegoils))
  )

# ------------------------------------------------------------
# 2. Prioritized proteins
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
  mutate(
    pdi_type = "hPDI",
    diet_col = "Hpdi_z"
  ) %>%
  filter(protein_id %in% names(df))

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
  mutate(
    pdi_type = "uPDI",
    diet_col = "Updi_z"
  ) %>%
  filter(protein_id %in% names(df))

protein_spec_all <- bind_rows(
  protein_spec_hpdi,
  protein_spec_updi
)

# ------------------------------------------------------------
# 3. Add selected age df from previous age-linearity results
# ------------------------------------------------------------

age_linearity_summary <- read.xlsx(
  age_linearity_file,
  sheet = "summary"
)

protein_spec_all <- protein_spec_all %>%
  left_join(
    age_linearity_summary %>%
      select(
        pdi_type,
        protein_id,
        selected_age_df,
        evidence_age_nonlinearity,
        age_df_selection_reason
      ),
    by = c("pdi_type", "protein_id")
  ) %>%
  mutate(
    selected_age_df = as.numeric(selected_age_df),
    selected_age_df = if_else(is.na(selected_age_df), 1, selected_age_df),
    use_age_spline = selected_age_df > 1
  )

# ------------------------------------------------------------
# 4. Function: formal diet shape test
# ------------------------------------------------------------

test_diet_shape_one_protein <- function(data,
                                        protein_col,
                                        gene,
                                        pdi_type_label,
                                        diet_col,
                                        selected_age_df = 1,
                                        diet_spline_df = 3) {
  
  age_term <- if (selected_age_df > 1) {
    paste0("ns(age_value, df = ", selected_age_df, ")")
  } else {
    "age_value"
  }
  
  covar_rhs <- paste(
    "Sex",
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
  
  # M0: no diet model
  f_null <- as.formula(
    paste0(
      "protein_value ~ ",
      covar_rhs
    )
  )
  
  # M1: linear diet model
  f_linear <- as.formula(
    paste0(
      "protein_value ~ diet_value + ",
      covar_rhs
    )
  )
  
  # M2: flexible diet model, fixed df = 3
  f_spline <- as.formula(
    paste0(
      "protein_value ~ ns(diet_value, df = ", diet_spline_df, ") + ",
      covar_rhs
    )
  )
  
  m_null <- lm(f_null, data = dat, weights = Weights)
  m_linear <- lm(f_linear, data = dat, weights = Weights)
  m_spline <- lm(f_spline, data = dat, weights = Weights)
  
  # Model comparisons
  comp_linear <- anova(m_null, m_linear)
  comp_spline <- anova(m_null, m_spline)
  comp_nonlinear <- anova(m_linear, m_spline)
  
  p_linear_overall <- comp_linear$`Pr(>F)`[2]
  p_spline_overall <- comp_spline$`Pr(>F)`[2]
  p_nonlinearity <- comp_nonlinear$`Pr(>F)`[2]
  
  # Linear beta for interpretability
  linear_tidy <- broom::tidy(m_linear)
  
  diet_beta <- linear_tidy %>%
    filter(term == "diet_value") %>%
    summarise(
      beta_linear = estimate[1],
      se_linear = std.error[1],
      t_linear = statistic[1],
      p_linear_term = p.value[1],
      .groups = "drop"
    )
  
  tibble(
    pdi_type = pdi_type_label,
    gene = gene,
    protein_id = protein_col,
    diet_col = diet_col,
    selected_age_df = selected_age_df,
    diet_spline_df = diet_spline_df,
    n = nobs(m_spline),
    
    AIC_null = AIC(m_null),
    AIC_linear = AIC(m_linear),
    AIC_spline = AIC(m_spline),
    
    BIC_null = BIC(m_null),
    BIC_linear = BIC(m_linear),
    BIC_spline = BIC(m_spline),
    
    adj_r2_null = summary(m_null)$adj.r.squared,
    adj_r2_linear = summary(m_linear)$adj.r.squared,
    adj_r2_spline = summary(m_spline)$adj.r.squared,
    
    delta_AIC_linear_vs_null = AIC(m_linear) - AIC(m_null),
    delta_AIC_spline_vs_null = AIC(m_spline) - AIC(m_null),
    delta_AIC_spline_vs_linear = AIC(m_spline) - AIC(m_linear),
    
    delta_BIC_linear_vs_null = BIC(m_linear) - BIC(m_null),
    delta_BIC_spline_vs_null = BIC(m_spline) - BIC(m_null),
    delta_BIC_spline_vs_linear = BIC(m_spline) - BIC(m_linear),
    
    p_linear_overall = p_linear_overall,
    p_spline_overall = p_spline_overall,
    p_nonlinearity = p_nonlinearity
  ) %>%
    bind_cols(diet_beta)
}

# ------------------------------------------------------------
# 5. Run formal diet shape tests
# ------------------------------------------------------------

diet_shape_results <- protein_spec_all %>%
  mutate(
    result = pmap(
      list(
        protein_id,
        gene,
        pdi_type,
        diet_col,
        selected_age_df
      ),
      ~ test_diet_shape_one_protein(
        data = df,
        protein_col = ..1,
        gene = ..2,
        pdi_type_label = ..3,
        diet_col = ..4,
        selected_age_df = ..5,
        diet_spline_df = 3
      )
    )
  ) %>%
  select(
    pdi_type,
    gene,
    protein_id,
    diet_col,
    selected_age_df,
    evidence_age_nonlinearity,
    age_df_selection_reason,
    result
  ) %>%
  select(-pdi_type, -gene, -protein_id, -diet_col, -selected_age_df) %>%
  unnest(result)

# ------------------------------------------------------------
# 6. Multiple-testing correction and classification
# ------------------------------------------------------------

diet_shape_summary <- diet_shape_results %>%
  group_by(pdi_type) %>%
  mutate(
    # Primary correction: within each diet index
    p_linear_overall_fdr_by_pdi = p.adjust(p_linear_overall, method = "BH"),
    p_spline_overall_fdr_by_pdi = p.adjust(p_spline_overall, method = "BH"),
    p_nonlinearity_fdr_by_pdi = p.adjust(p_nonlinearity, method = "BH"),
    p_linear_term_fdr_by_pdi = p.adjust(p_linear_term, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    # Sensitivity correction: across all hPDI + uPDI tests
    p_linear_overall_fdr_global = p.adjust(p_linear_overall, method = "BH"),
    p_spline_overall_fdr_global = p.adjust(p_spline_overall, method = "BH"),
    p_nonlinearity_fdr_global = p.adjust(p_nonlinearity, method = "BH"),
    p_linear_term_fdr_global = p.adjust(p_linear_term, method = "BH")
  ) %>%
  mutate(
    diet_shape_class_by_pdi_fdr = case_when(
      p_spline_overall_fdr_by_pdi < 0.05 &
        p_nonlinearity_fdr_by_pdi < 0.05 &
        p_linear_overall_fdr_by_pdi >= 0.05 ~
        "Non-linear-only diet association",
      
      p_spline_overall_fdr_by_pdi < 0.05 &
        p_nonlinearity_fdr_by_pdi < 0.05 &
        p_linear_overall_fdr_by_pdi < 0.05 ~
        "Linear plus non-linear diet association",
      
      p_linear_overall_fdr_by_pdi < 0.05 &
        p_nonlinearity_fdr_by_pdi >= 0.05 ~
        "Approximately linear diet association",
      
      p_spline_overall_fdr_by_pdi < 0.05 &
        p_nonlinearity_fdr_by_pdi >= 0.05 ~
        "Flexible association, but not clearly non-linear beyond linear",
      
      TRUE ~
        "No clear diet association"
    ),
    
    diet_shape_class_nominal = case_when(
      p_spline_overall < 0.05 &
        p_nonlinearity < 0.05 &
        p_linear_overall >= 0.05 ~
        "Nominal non-linear-only diet association",
      
      p_spline_overall < 0.05 &
        p_nonlinearity < 0.05 &
        p_linear_overall < 0.05 ~
        "Nominal linear plus non-linear diet association",
      
      p_linear_overall < 0.05 &
        p_nonlinearity >= 0.05 ~
        "Nominal approximately linear diet association",
      
      p_spline_overall < 0.05 &
        p_nonlinearity >= 0.05 ~
        "Nominal flexible association, not clearly non-linear",
      
      TRUE ~
        "No nominal diet association"
    ),
    
    selected_diet_df_for_plot = case_when(
      diet_shape_class_by_pdi_fdr %in% c(
        "Non-linear-only diet association",
        "Linear plus non-linear diet association"
      ) ~ diet_spline_df,
      TRUE ~ 1
    ),
    
    diet_df_selection_reason = case_when(
      selected_diet_df_for_plot > 1 ~
        "Formal fixed df=3 test supports diet non-linearity",
      TRUE ~
        "Linear diet retained for formal interpretation"
    )
  ) %>%
  arrange(
    pdi_type,
    p_spline_overall_fdr_by_pdi,
    p_nonlinearity_fdr_by_pdi
  )

if (interactive()) View(diet_shape_summary)

# ------------------------------------------------------------
# 7. Optional quick tables
# ------------------------------------------------------------

# Non-linear-only candidates: most relevant to XGBoost interpretation
nonlinear_only_candidates <- diet_shape_summary %>%
  filter(diet_shape_class_by_pdi_fdr == "Non-linear-only diet association") %>%
  arrange(pdi_type, p_spline_overall_fdr_by_pdi)

# Any formal non-linearity
nonlinear_candidates <- diet_shape_summary %>%
  filter(
    diet_shape_class_by_pdi_fdr %in% c(
      "Non-linear-only diet association",
      "Linear plus non-linear diet association"
    )
  ) %>%
  arrange(pdi_type, p_nonlinearity_fdr_by_pdi)

# Approximately linear candidates
linear_candidates <- diet_shape_summary %>%
  filter(diet_shape_class_by_pdi_fdr == "Approximately linear diet association") %>%
  arrange(pdi_type, p_linear_overall_fdr_by_pdi)

# No clear diet association
no_clear_diet_candidates <- diet_shape_summary %>%
  filter(diet_shape_class_by_pdi_fdr == "No clear diet association") %>%
  arrange(pdi_type, p_spline_overall_fdr_by_pdi)

# ------------------------------------------------------------
# 8. Save output
# ------------------------------------------------------------

write.xlsx(
  list(
    formal_diet_shape_tests = diet_shape_results,
    summary = diet_shape_summary,
    nonlinear_only = nonlinear_only_candidates,
    any_nonlinearity = nonlinear_candidates,
    approximately_linear = linear_candidates,
    no_clear_diet_association = no_clear_diet_candidates
  ),
  file = out_file,
  overwrite = TRUE
)
