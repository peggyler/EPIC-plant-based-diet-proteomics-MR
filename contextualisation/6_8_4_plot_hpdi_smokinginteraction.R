# Run from the repository root in a fresh R session. See README.md.
# Historical hPDI plot with diet as the outcome and manually selected protein terms.
# For the current protein-outcome plots for both diets, use 6_8_2_plotsmokinginteractions.R.
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

protein_spec_hpdi_smoke <- tribble(
  ~gene,    ~protein_id,       ~selected_df, ~smoking_interaction,
  "CTSS",   "seq.3181.50",     1,  TRUE,
  "DNPH1",  "seq.17785.11",    4,  TRUE,
  "PTK7",   "seq.9525.1",      1,  FALSE,  # nominal only
  "GREM1",  "seq.18878.15",    3,  TRUE,
  "CCM2",   "seq.12347.29",    4,  TRUE,
  "POSTN",  "seq.6645.53",     1,  TRUE,
  "CA11",   "seq.8859.51",     4,  FALSE,
  "MMP10",  "seq.10479.18",    4,  TRUE,
  "MMP10",  "seq.8479.4",      3,  FALSE,
  "GSTO1",  "seq.12436.84",    4,  TRUE,
  "LRIG1",  "seq.18831.6",     2,  TRUE,
  "MANBA",  "seq.6382.17",     2,  TRUE,
  "NCAN",   "seq.15573.110",   4,  TRUE,
  "PARVA",  "seq.13434.172",   3,  TRUE,
  "PCSK9",  "seq.5231.79",     2,  TRUE,
  "PCSK9",  "seq.8246.9",      3,  TRUE,
  "PEAR1",  "seq.8275.31",     2,  TRUE,
  "SPP1",   "seq.13113.7",     1,  FALSE
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

plot_smoke_model_one_protein <- function(data,
                                         protein_col,
                                         outcome_col,
                                         selected_df = 1,
                                         smoking_interaction = FALSE,
                                         diet_label = "hPDI",
                                         smoke_var = "ever_smoked") {
  
  vars_needed <- c(
    outcome_col, protein_col, "Weights",
    smoke_var,
    "Sex", "Age_Recr", "edu2", "alc4",
    "PA_total_MET_quint", "Center_grp"
  )
  
  dat <- data %>%
    select(all_of(vars_needed)) %>%
    mutate(
      outcome = .data[[outcome_col]],
      prot_value = .data[[protein_col]],
      smoke_group = .data[[smoke_var]]
    ) %>%
    filter(complete.cases(.))
  
  protein_term <- if (selected_df > 1) {
    paste0("ns(prot_value, df = ", selected_df, ")")
  } else {
    "prot_value"
  }
  
  protein_smoke_term <- if (smoking_interaction) {
    paste0(protein_term, " * smoke_group")
  } else {
    paste0(protein_term, " + smoke_group")
  }
  
  f_model <- as.formula(
    paste0(
      "outcome ~ ",
      "Sex + ns(Age_Recr, df = 3) + ",
      "edu2 + alc4 + PA_total_MET_quint + Center_grp + ",
      protein_smoke_term
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
    smoke_group = factor(levels(dat$smoke_group), levels = levels(dat$smoke_group))
  ) %>%
    mutate(
      Sex = make_ref_value(dat$Sex),
      Age_Recr = median(dat$Age_Recr, na.rm = TRUE),
      edu2 = make_ref_value(dat$edu2),
      alc4 = make_ref_value(dat$alc4),
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
      model_type = if_else(smoking_interaction, "Smoking interaction included", "No smoking interaction"),
      diet = diet_label
    )
}

## create plot data
plotdat_hpdi_smoke <- protein_spec_hpdi_smoke %>%
  mutate(
    plot_data = pmap(
      list(protein_id, selected_df, smoking_interaction),
      ~ plot_smoke_model_one_protein(
        data = df,
        protein_col = ..1,
        outcome_col = "Hpdi_z",
        selected_df = ..2,
        smoking_interaction = ..3,
        diet_label = "hPDI",
        smoke_var = "ever_smoked"
      )
    )
  ) %>%
  unnest(plot_data) %>%
  mutate(
    smoking_interaction_label = if_else(
      smoking_interaction,
      "Smoking interaction included",
      "No smoking interaction"
    ),
    facet_label = paste0(
      gene, " | ", protein_id
    )
  )

p_hpdi_smoke <- ggplot(
  plotdat_hpdi_smoke,
  aes(x = prot_value, y = pred, colour = smoke_group, fill = smoke_group)
) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.15, colour = NA) +
  geom_line(aes(linetype = smoke_group), linewidth = 0.8) +
  scale_colour_manual(
    values = c(
      "Never smoked" = "#2B3256",
      "Ever Smoked"  = "#D27E64"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Never smoked" = "#2B3256",
      "Ever Smoked"  = "#D27E64"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Never smoked" = "solid",
      "Ever Smoked"  = "dotted"
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
    colour = "Smoking status",
    fill = "Smoking status",
    linetype = "Smoking status",
    title = "Protein–hPDI associations by smoking status"
  ) +
  theme_bw()

p_hpdi_smoke

ggsave(
  filename = file.path(MANUSCRIPT_DIR, "hPDI_smoke_interaction_A4.pdf"),
  plot = p_hpdi_smoke,
  width = 8.27,
  height = 11.69,
  units = "in"
)
