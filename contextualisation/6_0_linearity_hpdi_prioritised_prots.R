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

library(dplyr)
library(purrr)
library(tidyr)
library(broom)
library(splines)
library(openxlsx)

##Linear regression with EPIC full cohort
load(file.path(DATA_DIR, "analyticaldf_postexcl_20251020.RData"))
wgts <- readRDS(file.path(DATA_DIR, "Weights_Feb2025.rds"))

wgts <- dplyr::rename(wgts, Idepic = idepic)

df1 <- df %>%
  left_join(wgts, by = "Idepic") 

df <- df1 %>%
  filter(!is.na(Weights))

##linear
df$Sex <- factor(df$Sex)
df$Sex <- relevel(df$Sex, ref = "Male")
df$PA <- as.numeric(scale(df$PA_total_MET))
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
df$edu2 <- factor(df$edu2)
df$edu2 <- relevel(df$edu2, ref = "Below secondary")
df$PA_total_MET_quint <- as.factor(df$PA_total_MET_quint)

# ------------------------------------------------------------
# 1. Scale hPDI
# ------------------------------------------------------------

df <- df %>%
  mutate(
    Hpdi_z = as.numeric(scale(Hpdi_NoVegoils))
  )

# ------------------------------------------------------------
# 2. hPDI prioritized proteins
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
# 3. Add selected age df from age-linearity results
# ------------------------------------------------------------
age_linearity_summary <- openxlsx::read.xlsx(
  file.path(MANUSCRIPT_DIR, "age_linearity_prioritised_proteins_correct_direction_AIC_df2to4.xlsx"),
  sheet = "summary"
)

protein_spec_hpdi <- protein_spec_hpdi %>%
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
  mutate(
    selected_age_df = if_else(is.na(selected_age_df), 1, selected_age_df),
    use_age_spline = selected_age_df > 1
  )

# ------------------------------------------------------------
# 4. Function: test linear hPDI vs spline hPDI
# ------------------------------------------------------------
test_hpdi_linearity_one_protein <- function(data,
                                            protein_col,
                                            selected_age_df = 1,
                                            diet_col = "Hpdi_z",
                                            spline_dfs = 2:4) {
  
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
  
  # Linear hPDI model
  f_linear <- as.formula(
    paste0(
      "protein_value ~ diet_value + ",
      covar_rhs
    )
  )
  
  m_linear <- lm(f_linear, data = dat, weights = Weights)
  
  linear_row <- tibble(
    protein_id = protein_col,
    model = "linear_hPDI",
    hpdi_spline_df = 1,
    selected_age_df = selected_age_df,
    n = nobs(m_linear),
    AIC = AIC(m_linear),
    BIC = BIC(m_linear),
    adj_r2 = summary(m_linear)$adj.r.squared,
    p_vs_linear = NA_real_
  )
  
  # Spline hPDI models, df 2–4
  spline_rows <- map_dfr(spline_dfs, function(df_spline) {
    
    f_spline <- as.formula(
      paste0(
        "protein_value ~ ns(diet_value, df = ", df_spline, ") + ",
        covar_rhs
      )
    )
    
    m_spline <- lm(f_spline, data = dat, weights = Weights)
    comp <- anova(m_linear, m_spline)
    
    tibble(
      protein_id = protein_col,
      model = paste0("hPDI_spline_df", df_spline),
      hpdi_spline_df = df_spline,
      selected_age_df = selected_age_df,
      n = nobs(m_spline),
      AIC = AIC(m_spline),
      BIC = BIC(m_spline),
      adj_r2 = summary(m_spline)$adj.r.squared,
      p_vs_linear = comp$`Pr(>F)`[2]
    )
  })
  
  bind_rows(linear_row, spline_rows) %>%
    mutate(
      delta_AIC_vs_linear = AIC - AIC[model == "linear_hPDI"],
      delta_BIC_vs_linear = BIC - BIC[model == "linear_hPDI"]
    )
}

# ------------------------------------------------------------
# 5. Run hPDI linearity tests
# ------------------------------------------------------------
linearity_hpdi <- protein_spec_hpdi %>%
  mutate(
    result = map2(
      protein_id,
      selected_age_df,
      ~ test_hpdi_linearity_one_protein(
        data = df,
        protein_col = .x,
        selected_age_df = .y,
        diet_col = "Hpdi_z",
        spline_dfs = 2:4
      ) %>%
        select(-protein_id, -selected_age_df)   # remove duplicates
    )
  ) %>%
  select(
    pdi_type,
    gene,
    protein_id,
    selected_age_df,
    evidence_age_nonlinearity,
    age_df_selection_reason,
    result
  ) %>%
  unnest(result)

# ------------------------------------------------------------
# 6. Summarise hPDI linearity results
# ------------------------------------------------------------
linearity_hpdi_summary <- linearity_hpdi %>%
  group_by(
    pdi_type,
    gene,
    protein_id,
    selected_age_df,
    evidence_age_nonlinearity,
    age_df_selection_reason
  ) %>%
  summarise(
    n = first(n),
    
    best_df_AIC = hpdi_spline_df[which.min(AIC)],
    best_model_AIC = model[which.min(AIC)],
    best_AIC = min(AIC),
    
    best_df_BIC = hpdi_spline_df[which.min(BIC)],
    best_model_BIC = model[which.min(BIC)],
    best_BIC = min(BIC),
    
    linear_AIC = AIC[model == "linear_hPDI"],
    linear_BIC = BIC[model == "linear_hPDI"],
    
    min_p_vs_linear = min(p_vs_linear, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    p_hpdi_nonlinearity_fdr = p.adjust(min_p_vs_linear, method = "BH"),
    
    evidence_hpdi_nonlinearity = case_when(
      p_hpdi_nonlinearity_fdr < 0.05 ~ "Strong evidence of hPDI non-linearity",
      min_p_vs_linear < 0.05 ~ "Nominal evidence of hPDI non-linearity",
      TRUE ~ "No clear evidence of hPDI non-linearity"
    ),
    
    selected_hpdi_df = case_when(
      evidence_hpdi_nonlinearity == "Strong evidence of hPDI non-linearity" &
        best_model_AIC != "linear_hPDI" ~ best_df_AIC,
      TRUE ~ 1
    ),
    
    hpdi_df_selection_reason = case_when(
      evidence_hpdi_nonlinearity == "Strong evidence of hPDI non-linearity" &
        best_model_AIC != "linear_hPDI" ~ "Selected by AIC among hPDI spline df 2–4",
      TRUE ~ "Linear hPDI retained"
    )
  ) %>%
  arrange(min_p_vs_linear)

if (interactive()) View(linearity_hpdi_summary)
# ------------------------------------------------------------
# 7. Save
# ------------------------------------------------------------

write.xlsx(
  list(
    all_models = linearity_hpdi,
    summary = linearity_hpdi_summary
  ),
  file = file.path(MANUSCRIPT_DIR, "linearity_hpdi_correct_direction_ageDFplugged.xlsx"),
  overwrite = TRUE
)
