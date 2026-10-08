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
library(ggplot2)
library(purrr)

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

##
df <- df %>%
  mutate(
    Hpdi_z = as.numeric(scale(Hpdi_NoVegoils)),
    Updi_z = as.numeric(scale(Updi_NoVegoils))
  )
##

library(dplyr)
library(tibble)

protein_spec_updi <- tribble(
  ~gene,       ~protein_id,       ~selected_df, ~sex_interaction,
  "KRT5",      "seq.11177.16",    1,  FALSE,
  "EGFR",      "seq.2677.1",      1,  TRUE,
  "CPNE1",     "seq.11602.12",    1,  FALSE,
  "PLA2G2A",   "seq.2692.74",     4,  FALSE,
  "SERPINA3",  "seq.2879.9",      4,  TRUE,
  "HHIP",      "seq.10833.64",    1,  TRUE,
  "CCM2",      "seq.12347.29",    1,  TRUE,
  "ACE",       "seq.10714.7",     1,  FALSE,   # nominal
  "LCT",       "seq.9017.58",     1,  FALSE,
  "PTPN9",     "seq.12633.3",     4,  FALSE,
  "POSTN",     "seq.3457.57",     2,  FALSE,
  "POSTN",     "seq.6645.53",     4,  TRUE,
  "MMP10",     "seq.10479.18",    1,  FALSE,
  "MMP10",     "seq.8479.4",      1,  FALSE,
  "PCSK9",     "seq.8246.9",      5,  FALSE,
  "ACADVL",    "seq.19262.219",   4,  TRUE,
  "NADK",      "seq.13624.17",    2,  FALSE,
  "SPP1",      "seq.13113.7",     1,  FALSE
) %>%
  mutate(
    use_spline = selected_df > 1
  ) %>%
  filter(protein_id %in% names(df))

##plotting function
get_mode <- function(x) {
  x <- x[!is.na(x)]
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

make_ref_value <- function(x) {
  if (is.factor(x)) {
    factor(get_mode(x), levels = levels(x))
  } else if (is.character(x)) {
    get_mode(x)
  } else {
    median(x, na.rm = TRUE)
  }
}

plot_sex_model_one_protein <- function(data,
                                       protein_col,
                                       outcome_col,
                                       selected_df = 1,
                                       sex_interaction = FALSE,
                                       diet_label = "uPDI") {
  
  vars_needed <- c(
    outcome_col, protein_col, "Weights",
    "Sex", "Age_Recr", "edu2", "alc4", "ever_smoked",
    "PA_total_MET_quint", "Center_grp"
  )
  
  dat <- data %>%
    select(all_of(vars_needed)) %>%
    mutate(
      outcome = .data[[outcome_col]],
      prot_value = .data[[protein_col]]
    ) %>%
    filter(complete.cases(.))
  
  protein_term <- if (selected_df > 1) {
    paste0("ns(prot_value, df = ", selected_df, ")")
  } else {
    "prot_value"
  }
  
  protein_sex_term <- if (sex_interaction) {
    paste0(protein_term, " * Sex")
  } else {
    paste0(protein_term, " + Sex")
  }
  
  f_model <- as.formula(
    paste0(
      "outcome ~ ",
      "ns(Age_Recr, df = 3) + ",
      "edu2 + alc4 + ever_smoked + PA_total_MET_quint + Center_grp + ",
      protein_sex_term
    )
  )
  
  fit <- lm(f_model, data = dat, weights = Weights)
  
  x_seq <- seq(
    quantile(dat$prot_value, 0.05, na.rm = TRUE),
    quantile(dat$prot_value, 0.95, na.rm = TRUE),
    length.out = 100
  )
  
  newdat <- expand_grid(
    prot_value = x_seq,
    Sex = factor(levels(dat$Sex), levels = levels(dat$Sex))
  ) %>%
    mutate(
      Age_Recr = median(dat$Age_Recr, na.rm = TRUE),
      edu2 = make_ref_value(dat$edu2),
      alc4 = make_ref_value(dat$alc4),
      ever_smoked = make_ref_value(dat$ever_smoked),
      PA_total_MET_quint = make_ref_value(dat$PA_total_MET_quint),
      Center_grp = make_ref_value(dat$Center_grp)
    )
  
  pred_obj <- stats::predict(
    fit,
    newdata = newdat,
    se.fit = TRUE
  )
  
  pred_fit <- as.numeric(pred_obj$fit)
  pred_se  <- as.numeric(pred_obj$se.fit)
  
  newdat %>%
    mutate(
      pred = pred_fit,
      se   = pred_se,
      lcl  = pred_fit - 1.96 * pred_se,
      ucl  = pred_fit + 1.96 * pred_se,
      protein_form = if_else(selected_df > 1, paste0("spline df=", selected_df), "linear"),
      model_type = if_else(sex_interaction, "Sex interaction included", "Main effects only"),
      diet = diet_label
    )
}
  
## create plot data
plotdat_updi_sex <- protein_spec_updi %>%
  mutate(
    plot_data = pmap(
      list(protein_id, selected_df, sex_interaction),
      ~ plot_sex_model_one_protein(
        data = df,
        protein_col = ..1,
        outcome_col = "Updi_z",
        selected_df = ..2,
        sex_interaction = ..3,
        diet_label = "uPDI"
      ) %>%
        select(-protein_form, -model_type, -diet)
    )
  ) %>%
  unnest(plot_data) %>%
  mutate(
    sex_interaction_label = if_else(
      sex_interaction,
      "Sex interaction included",
      "No sex interaction"
    ),
    facet_label = paste0(
      gene, " | ", protein_id
    )
  )


p_updi_sex <- ggplot(
  plotdat_updi_sex,
  aes(x = prot_value, y = pred, colour = Sex, fill = Sex)
) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.15, colour = NA) +
  geom_line(aes(linetype = Sex), linewidth = 0.8) +
  scale_colour_manual(
    values = c(
      "Female" = "#D27E64",
      "Male"   = "#2B3256"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#D27E64",
      "Male"   = "#2B3256"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Female" = "solid",
      "Male"   = "dashed"
    )
  ) +
  facet_wrap(
    ~ facet_label,
    scales = "free_x",
    ncol = 3
  ) +
  labs(
    x = "Protein abundance",
    y = "Predicted standardized hPDI",
    title = "Protein–uPDI associations by sex"
  ) +
  theme_bw()

p_updi_sex

ggsave(
  filename = file.path(MANUSCRIPT_DIR, "uPDI_sex_interaction_A4.pdf"),
  plot = p_updi_sex,
  width = 8.27,
  height = 11.69,
  units = "in"
)
