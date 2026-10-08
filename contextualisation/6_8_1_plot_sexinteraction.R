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
# Plot sex-specific diet-protein associations
# Correct direction:
# protein abundance ~ dietary index * Sex + selected age term + covariates
#
# Sex interaction is included ONLY when:
# evidence_sex_interaction == "Strong evidence of sex interaction"
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

hpdi_sex_file <- file.path(
  out_dir,
  "sex_interaction_hpdi_correct_direction_ageHpdiDFplugged.xlsx"
)

updi_sex_file <- file.path(
  out_dir,
  "sex_interaction_updi_correct_direction_ageUpdiDFplugged.xlsx"
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
# 2. Read sex interaction results
# ------------------------------------------------------------

sex_hpdi <- read.xlsx(hpdi_sex_file, sheet = "sex_interaction") %>%
  mutate(
    pdi_type = "hPDI",
    selected_diet_df = selected_hpdi_df,
    evidence_diet_nonlinearity = evidence_hpdi_nonlinearity,
    diet_df_selection_reason = hpdi_df_selection_reason,
    sex_interaction = evidence_sex_interaction == "Strong evidence of sex interaction"
  )

sex_updi <- read.xlsx(updi_sex_file, sheet = "sex_interaction") %>%
  mutate(
    pdi_type = "uPDI",
    selected_diet_df = selected_updi_df,
    evidence_diet_nonlinearity = evidence_updi_nonlinearity,
    diet_df_selection_reason = updi_df_selection_reason,
    sex_interaction = evidence_sex_interaction == "Strong evidence of sex interaction"
  )

sex_all <- bind_rows(sex_hpdi, sex_updi) %>%
  mutate(
    selected_age_df = as.numeric(selected_age_df),
    selected_diet_df = as.numeric(selected_diet_df),
    selected_age_df = coalesce(selected_age_df, 1),
    selected_diet_df = coalesce(selected_diet_df, 1)
  )

# Optional: set TRUE if you only want to plot proteins with strong sex interaction
plot_only_strong <- FALSE

sex_all_plot <- sex_all %>%
  filter(!plot_only_strong | sex_interaction)

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

make_sex_plot_data_one_protein <- function(data,
                                           protein_col,
                                           gene,
                                           pdi_type,
                                           selected_age_df = 1,
                                           selected_diet_df = 1,
                                           sex_interaction = FALSE) {
  
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
  
  diet_sex_term <- if (sex_interaction) {
    paste0(diet_term, " * Sex")
  } else {
    paste0(diet_term, " + Sex")
  }
  
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
  
  f_model <- as.formula(
    paste0(
      "protein_value ~ ",
      diet_sex_term, " + ",
      covar_rhs_no_sex
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
    Sex = factor(levels(dat$Sex), levels = levels(dat$Sex))
  ) %>%
    mutate(
      age_value = median(dat$age_value, na.rm = TRUE),
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
        sex_interaction,
        "Sex interaction included",
        "No sex interaction"
      ),
      facet_label = paste0(.env$gene, " | ", .env$protein_col)
    )
}

# ------------------------------------------------------------
# 5. Create plot data
# ------------------------------------------------------------

plotdat_sex <- sex_all_plot %>%
  mutate(
    plot_data = purrr::pmap(
      list(
        protein_id,
        gene,
        pdi_type,
        selected_age_df,
        selected_diet_df,
        sex_interaction
      ),
      ~ make_sex_plot_data_one_protein(
        data = df,
        protein_col = ..1,
        gene = ..2,
        pdi_type = ..3,
        selected_age_df = ..4,
        selected_diet_df = ..5,
        sex_interaction = ..6
      )
    )
  ) %>%
  select(
    pdi_type,
    gene,
    protein_id,
    sex_interaction,
    p_interaction_sex,
    p_interaction_sex_fdr,
    evidence_sex_interaction,
    plot_data
  ) %>%
  unnest(plot_data)
# ------------------------------------------------------------
# 6. Split hPDI and uPDI
# ------------------------------------------------------------

plotdat_hpdi_sex <- plotdat_sex %>%
  filter(pdi_type == "hPDI")

plotdat_updi_sex <- plotdat_sex %>%
  filter(pdi_type == "uPDI")

# ------------------------------------------------------------
# 7. Plot hPDI
# ------------------------------------------------------------

p_hpdi_sex <- ggplot(
  plotdat_hpdi_sex,
  aes(x = diet_value, y = pred, colour = Sex, fill = Sex)
) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.15, colour = NA) +
  geom_line(aes(linetype = Sex), linewidth = 0.8) +
  scale_colour_manual(
    values = c(
      "Female" = "#D27E64",
      "Male" = "#2B3256"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#D27E64",
      "Male" = "#2B3256"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Female" = "solid",
      "Male" = "dotted"
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
    colour = "Sex",
    fill = "Sex",
    linetype = "Sex",
    title = "Sex-specific hPDI–protein associations"
  ) +
  theme_bw()

p_hpdi_sex

ggsave(
  filename = file.path(out_dir, "hPDI_sex_interaction_correct_direction_A4.pdf"),
  plot = p_hpdi_sex,
  width = 8.27,
  height = 11.69,
  units = "in"
)

# ------------------------------------------------------------
# 8. Plot uPDI
# ------------------------------------------------------------

p_updi_sex <- ggplot(
  plotdat_updi_sex,
  aes(x = diet_value, y = pred, colour = Sex, fill = Sex)
) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.15, colour = NA) +
  geom_line(aes(linetype = Sex), linewidth = 0.8) +
  scale_colour_manual(
    values = c(
      "Female" = "#D27E64",
      "Male" = "#2B3256"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Female" = "#D27E64",
      "Male" = "#2B3256"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Female" = "solid",
      "Male" = "dotted"
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
    colour = "Sex",
    fill = "Sex",
    linetype = "Sex",
    title = "Sex-specific uPDI–protein associations"
  ) +
  theme_bw()

p_updi_sex

ggsave(
  filename = file.path(out_dir, "uPDI_sex_interaction_correct_direction_A4.pdf"),
  plot = p_updi_sex,
  width = 8.27,
  height = 11.69,
  units = "in"
)
