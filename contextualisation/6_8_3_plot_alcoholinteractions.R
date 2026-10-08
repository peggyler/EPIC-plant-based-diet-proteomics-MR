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
# Reconstructed from the supplied interaction plotting workflow.
# Model selections come from 6_6_ALC_hPDI.R and 6_7_alcinteraction_uPDI.R.
# Compare regenerated figures with the article; this is not recovered historical code.
# Plot alcohol-specific diet-protein associations
# Correct direction:
# protein abundance ~ dietary index * alcohol + selected age term + covariates
#
# Alcohol interaction is included ONLY when:
# evidence_alcohol_interaction == "Strong evidence of alcohol interaction"
# ============================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(splines)
library(openxlsx)

# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

data_file <- file.path(DATA_DIR, "analyticaldf_postexcl_20251020.RData")
weights_file <- file.path(DATA_DIR, "Weights_Feb2025.rds")

out_dir <- MANUSCRIPT_DIR

hpdi_alcohol_file <- file.path(
  out_dir,
  "alcohol_interaction_hpdi_correct_direction_ageHpdiDFplugged.xlsx"
)

updi_alcohol_file <- file.path(
  out_dir,
  "alcohol_interaction_updi_correct_direction_ageUpdiDFplugged.xlsx"
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
    Hpdi_z = as.numeric(scale(Hpdi_NoVegoils)),
    Updi_z = as.numeric(scale(Updi_NoVegoils))
  )

# ------------------------------------------------------------
# 2. Read alcohol interaction results
# ------------------------------------------------------------

alcohol_hpdi <- read.xlsx(hpdi_alcohol_file, sheet = "alcohol_interaction") %>%
  mutate(
    pdi_type = "hPDI",
    selected_diet_df = selected_hpdi_df,
    evidence_diet_nonlinearity = evidence_hpdi_nonlinearity,
    diet_df_selection_reason = hpdi_df_selection_reason,
    alcohol_interaction = evidence_alcohol_interaction == "Strong evidence of alcohol interaction"
  )

alcohol_updi <- read.xlsx(updi_alcohol_file, sheet = "alcohol_interaction") %>%
  mutate(
    pdi_type = "uPDI",
    selected_diet_df = selected_updi_df,
    evidence_diet_nonlinearity = evidence_updi_nonlinearity,
    diet_df_selection_reason = updi_df_selection_reason,
    alcohol_interaction = evidence_alcohol_interaction == "Strong evidence of alcohol interaction"
  )

alcohol_all <- bind_rows(alcohol_hpdi, alcohol_updi) %>%
  mutate(
    selected_age_df = as.numeric(selected_age_df),
    selected_diet_df = as.numeric(selected_diet_df),
    selected_age_df = coalesce(selected_age_df, 1),
    selected_diet_df = coalesce(selected_diet_df, 1)
  )

# Set TRUE if you only want to plot proteins with strong alcohol interaction
plot_only_strong <- FALSE

alcohol_all_plot <- alcohol_all %>%
  filter(!plot_only_strong | alcohol_interaction)

# ------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# 4. Function to create prediction data
# ------------------------------------------------------------

make_alcohol_plot_data_one_protein <- function(data,
                                               protein_col,
                                               gene,
                                               pdi_type,
                                               selected_age_df = 1,
                                               selected_diet_df = 1,
                                               alcohol_interaction = FALSE) {
  
  diet_col <- ifelse(pdi_type == "hPDI", "Hpdi_z", "Updi_z")
  
  age_term <- if (selected_age_df > 1) {
    paste0("ns(age_value, df = ", selected_age_df, ")")
  } else {
    "age_value"
  }
  
  diet_term <- if (selected_diet_df > 1) {
    paste0("ns(diet_value, df = ", selected_diet_df, ")")
  } else {
    "diet_value"
  }
  
  diet_alc_term <- if (alcohol_interaction) {
    paste0(diet_term, " * alc_group")
  } else {
    paste0(diet_term, " + alc_group")
  }
  
  covar_rhs_no_alcohol <- paste(
    "Sex",
    age_term,
    "edu2",
    "ever_smoked",
    "PA_total_MET_quint",
    "Center_grp",
    sep = " + "
  )
  
  vars_needed <- c(
    protein_col,
    diet_col,
    "Weights",
    "ever_smoked",
    "Sex",
    "Age_Recr",
    "edu2",
    "alc4",
    "PA_total_MET_quint",
    "Center_grp"
  )
  
  dat <- data %>%
    select(all_of(vars_needed)) %>%
    mutate(
      protein_value = .data[[protein_col]],
      diet_value = .data[[diet_col]],
      age_value = Age_Recr,
      alc_group = alc4
    ) %>%
    filter(complete.cases(.))
  
  f_model <- as.formula(
    paste0(
      "protein_value ~ ",
      diet_alc_term, " + ",
      covar_rhs_no_alcohol
    )
  )
  
  fit <- lm(f_model, data = dat, weights = Weights)
  
  x_seq <- seq(
    quantile(dat$diet_value, 0.05, na.rm = TRUE),
    quantile(dat$diet_value, 0.95, na.rm = TRUE),
    length.out = 100
  )
  
  newdat <- expand_grid(
    diet_value = x_seq,
    alc_group = factor(levels(dat$alc_group), levels = levels(dat$alc_group))
  ) %>%
    mutate(
      Sex = make_ref_value(dat$Sex),
      age_value = median(dat$age_value, na.rm = TRUE),
      edu2 = make_ref_value(dat$edu2),
      ever_smoked = make_ref_value(dat$ever_smoked),
      PA_total_MET_quint = make_ref_value(dat$PA_total_MET_quint),
      Center_grp = make_ref_value(dat$Center_grp)
    )
  
  pred_obj <- predict(fit, newdata = newdat, se.fit = TRUE)
  
  pred_fit <- as.numeric(pred_obj$fit)
  pred_se <- as.numeric(pred_obj$se.fit)
  
  newdat %>%
    mutate(
      pred = pred_fit,
      se = pred_se,
      lcl = pred_fit - 1.96 * pred_se,
      ucl = pred_fit + 1.96 * pred_se,
      selected_age_df = selected_age_df,
      selected_diet_df = selected_diet_df,
      age_form = if_else(
        selected_age_df > 1,
        paste0("Age spline df=", selected_age_df),
        "Linear age"
      ),
      diet_form = if_else(
        selected_diet_df > 1,
        paste0(pdi_type, " spline df=", selected_diet_df),
        paste0("Linear ", pdi_type)
      ),
      interaction_label = if_else(
        alcohol_interaction,
        "Alcohol interaction included",
        "No alcohol interaction"
      ),
      facet_label = paste0(.env$gene, " | ", .env$protein_col)
    )
}

# ------------------------------------------------------------
# 5. Create plot data
# ------------------------------------------------------------

plotdat_alcohol <- alcohol_all_plot %>%
  mutate(
    plot_data = pmap(
      list(
        protein_id,
        gene,
        pdi_type,
        selected_age_df,
        selected_diet_df,
        alcohol_interaction
      ),
      ~ make_alcohol_plot_data_one_protein(
        data = df,
        protein_col = ..1,
        gene = ..2,
        pdi_type = ..3,
        selected_age_df = ..4,
        selected_diet_df = ..5,
        alcohol_interaction = ..6
      )
    )
  ) %>%
  select(
    pdi_type,
    gene,
    protein_id,
    alcohol_interaction,
    p_interaction_alcohol,
    p_interaction_alcohol_fdr,
    evidence_alcohol_interaction,
    plot_data
  ) %>%
  unnest(plot_data)

# ------------------------------------------------------------
# 6. Plot hPDI and uPDI
# ------------------------------------------------------------

alcohol_levels <- c("Never", "Occasional", "Moderate", "Regular/Daily")
alcohol_colours <- c("Never" = "#003B73", "Occasional" = "#7A5195",
                     "Moderate" = "#E66100", "Regular/Daily" = "#A50026")
alcohol_linetypes <- c("Never" = "solid", "Occasional" = "dashed",
                       "Moderate" = "dotted", "Regular/Daily" = "dotdash")

make_alcohol_plot <- function(plot_data, diet_label) {
  ggplot(plot_data, aes(x = diet_value, y = pred, colour = alc_group, fill = alc_group)) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.15, colour = NA) +
    geom_line(aes(linetype = alc_group), linewidth = 0.8) +
    scale_colour_manual(values = alcohol_colours, breaks = alcohol_levels) +
    scale_fill_manual(values = alcohol_colours, breaks = alcohol_levels) +
    scale_linetype_manual(values = alcohol_linetypes, breaks = alcohol_levels) +
    facet_wrap(~ facet_label, scales = "free_y", ncol = 3) +
    labs(
      x = paste("Standardized", diet_label),
      y = "Predicted protein abundance",
      colour = "Alcohol frequency", fill = "Alcohol frequency", linetype = "Alcohol frequency",
      title = paste0("Alcohol-specific ", diet_label, "–protein associations")
    ) +
    theme_bw()
}

plotdat_hpdi_alcohol <- plotdat_alcohol %>% filter(pdi_type == "hPDI")
plotdat_updi_alcohol <- plotdat_alcohol %>% filter(pdi_type == "uPDI")
p_hpdi_alcohol <- make_alcohol_plot(plotdat_hpdi_alcohol, "hPDI")
p_updi_alcohol <- make_alcohol_plot(plotdat_updi_alcohol, "uPDI")

ggsave(file.path(out_dir, "hPDI_alcohol_interaction_correct_direction_A4.pdf"),
       plot = p_hpdi_alcohol, width = 8.27, height = 11.69, units = "in")
ggsave(file.path(out_dir, "uPDI_alcohol_interaction_correct_direction_A4.pdf"),
       plot = p_updi_alcohol, width = 8.27, height = 11.69, units = "in")
