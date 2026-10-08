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
# Plot age-specific diet-protein associations
# Correct direction:
# protein abundance ~ dietary index * age + covariates
#
# Age interaction is included ONLY when:
# evidence_age_interaction == "Strong evidence of age interaction"
#
# Predictions are shown at ages 50, 55, and 60
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

hpdi_age_file <- file.path(
  out_dir,
  "age_interaction_hpdi_correct_direction_ageHpdiDFplugged.xlsx"
)

updi_age_file <- file.path(
  out_dir,
  "age_interaction_updi_correct_direction_ageUpdiDFplugged.xlsx"
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
  levels = c("Never", "Occasional", "Moderate", "Regular/Daily")
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
# 2. Read age interaction results
# ------------------------------------------------------------

age_hpdi <- read.xlsx(hpdi_age_file, sheet = "age_interaction") %>%
  mutate(
    pdi_type = "hPDI",
    selected_diet_df = selected_hpdi_df,
    evidence_diet_nonlinearity = evidence_hpdi_nonlinearity,
    diet_df_selection_reason = hpdi_df_selection_reason,
    age_interaction = evidence_age_interaction == "Strong evidence of age interaction"
  )

age_updi <- read.xlsx(updi_age_file, sheet = "age_interaction") %>%
  mutate(
    pdi_type = "uPDI",
    selected_diet_df = selected_updi_df,
    evidence_diet_nonlinearity = evidence_updi_nonlinearity,
    diet_df_selection_reason = updi_df_selection_reason,
    age_interaction = evidence_age_interaction == "Strong evidence of age interaction"
  )

age_all <- bind_rows(age_hpdi, age_updi) %>%
  mutate(
    selected_age_df = as.numeric(selected_age_df),
    selected_diet_df = as.numeric(selected_diet_df),
    selected_age_df = coalesce(selected_age_df, 1),
    selected_diet_df = coalesce(selected_diet_df, 1)
  )

# Set TRUE if you only want to plot proteins with strong age interaction
plot_only_strong <- FALSE

age_all_plot <- age_all %>%
  filter(!plot_only_strong | age_interaction)

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

make_age_plot_data_one_protein <- function(data,
                                           protein_col,
                                           gene,
                                           pdi_type,
                                           selected_age_df = 1,
                                           selected_diet_df = 1,
                                           age_interaction = FALSE,
                                           age_levels = c(50, 55, 60)) {
  
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
  
  diet_age_term <- if (age_interaction) {
    paste0(diet_term, " * ", age_term)
  } else {
    paste0(diet_term, " + ", age_term)
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
  
  f_model <- as.formula(
    paste0(
      "protein_value ~ ",
      diet_age_term, " + ",
      covar_rhs_no_age_diet
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
    age_value = age_levels
  ) %>%
    mutate(
      age_group = factor(
        paste0("Age ", age_value),
        levels = paste0("Age ", age_levels)
      ),
      Sex = make_ref_value(dat$Sex),
      edu2 = make_ref_value(dat$edu2),
      alc4 = make_ref_value(dat$alc4),
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
        age_interaction,
        "Age interaction included",
        "No age interaction"
      ),
      facet_label = paste0(.env$gene, " | ", .env$protein_col)
    )
}

# ------------------------------------------------------------
# 5. Create plot data
# ------------------------------------------------------------

plotdat_age <- age_all_plot %>%
  mutate(
    plot_data = pmap(
      list(
        protein_id,
        gene,
        pdi_type,
        selected_age_df,
        selected_diet_df,
        age_interaction
      ),
      ~ make_age_plot_data_one_protein(
        data = df,
        protein_col = ..1,
        gene = ..2,
        pdi_type = ..3,
        selected_age_df = ..4,
        selected_diet_df = ..5,
        age_interaction = ..6,
        age_levels = c(50, 55, 60)
      )
    )
  ) %>%
  select(
    pdi_type,
    gene,
    protein_id,
    age_interaction,
    p_interaction_age,
    p_interaction_age_fdr,
    evidence_age_interaction,
    plot_data
  ) %>%
  unnest(plot_data)

# ------------------------------------------------------------
# 6. Split hPDI and uPDI
# ------------------------------------------------------------

plotdat_hpdi_age <- plotdat_age %>%
  filter(pdi_type == "hPDI")

plotdat_updi_age <- plotdat_age %>%
  filter(pdi_type == "uPDI")

p_hpdi_age <- ggplot(
  plotdat_hpdi_age,
  aes(x = diet_value, y = pred, colour = age_group, fill = age_group)
) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.10, colour = NA) +
  geom_line(aes(linetype = age_group), linewidth = 0.8) +
  scale_colour_manual(
    values = c(
      "Age 50" = "#D27E64",
      "Age 55" = "#7A869A",
      "Age 60" = "#2B3256"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Age 50" = "#D27E64",
      "Age 55" = "#7A869A",
      "Age 60" = "#2B3256"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Age 50" = "dotted",
      "Age 55" = "dashed",
      "Age 60" = "solid"
    )
  ) +
  facet_wrap(
    ~ facet_label,
    scales = "free_y",
    ncol = 3
  ) +
  labs(
    x = "Standardized hPDI",
    y = "Predicted protein abundance",
    colour = "Age at recruitment",
    fill = "Age at recruitment",
    linetype = "Age at recruitment",
    title = "Age-specific hPDI–protein associations"
  ) +
  theme_bw()

p_hpdi_age

ggsave(
  filename = file.path(out_dir, "hPDI_age_interaction_correct_direction_A4.pdf"),
  plot = p_hpdi_age,
  width = 8.27,
  height = 11.69,
  units = "in"
)

p_updi_age <- ggplot(
  plotdat_updi_age,
  aes(x = diet_value, y = pred, colour = age_group, fill = age_group)
) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.10, colour = NA) +
  geom_line(aes(linetype = age_group), linewidth = 0.8) +
  scale_colour_manual(
    values = c(
      "Age 50" = "#D27E64",
      "Age 55" = "#7A869A",
      "Age 60" = "#2B3256"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Age 50" = "#D27E64",
      "Age 55" = "#7A869A",
      "Age 60" = "#2B3256"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Age 50" = "dotted",
      "Age 55" = "dashed",
      "Age 60" = "solid"
    )
  ) +
  facet_wrap(
    ~ facet_label,
    scales = "free_y",
    ncol = 3
  ) +
  labs(
    x = "Standardized uPDI",
    y = "Predicted protein abundance",
    colour = "Age at recruitment",
    fill = "Age at recruitment",
    linetype = "Age at recruitment",
    title = "Age-specific uPDI–protein associations"
  ) +
  theme_bw()

p_updi_age

ggsave(
  filename = file.path(out_dir, "uPDI_age_interaction_correct_direction_A4.pdf"),
  plot = p_updi_age,
  width = 8.27,
  height = 11.69,
  units = "in"
)
