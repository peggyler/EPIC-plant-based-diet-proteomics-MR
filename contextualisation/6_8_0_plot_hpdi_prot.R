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
# Plot diet-protein associations
# one adjusted for age and sex, the other fully adjusted
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
#    No sex interaction.
#    Compare:
#    Model 1: adjusted for age + sex
#    Model 2: fully adjusted
# ------------------------------------------------------------

make_adjusted_plot_data_one_protein <- function(data,
                                                protein_col,
                                                gene,
                                                pdi_type,
                                                selected_age_df = 1,
                                                selected_diet_df = 1,
                                                model_type = c("age_sex", "full")) {
  
  model_type <- match.arg(model_type)
  
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
  
  base_vars <- c(
    protein_col,
    diet_col,
    "Weights",
    "Sex",
    "Age_Recr"
  )
  
  full_covars <- c(
    "edu2",
    "alc4",
    "ever_smoked",
    "PA_total_MET_quint",
    "Center_grp"
  )
  
  vars_needed <- c(base_vars, full_covars)
  
  dat <- data %>%
    select(all_of(vars_needed)) %>%
    mutate(
      protein_value = .data[[protein_col]],
      diet_value = .data[[diet_col]],
      age_value = Age_Recr
    ) %>%
    filter(complete.cases(.))
  
  covar_rhs <- if (model_type == "age_sex") {
    paste(
      age_term,
      "Sex",
      sep = " + "
    )
  } else {
    paste(
      age_term,
      "Sex",
      "edu2",
      "alc4",
      "ever_smoked",
      "PA_total_MET_quint",
      "Center_grp",
      sep = " + "
    )
  }
  
  f_model <- as.formula(
    paste0(
      "protein_value ~ ",
      diet_term, " + ",
      covar_rhs
    )
  )
  
  fit <- lm(f_model, data = dat, weights = Weights)
  
  x_seq <- seq(
    quantile(dat$diet_value, 0.05, na.rm = TRUE),
    quantile(dat$diet_value, 0.95, na.rm = TRUE),
    length.out = 100
  )
  
  newdat <- tibble(
    diet_value = x_seq
  ) %>%
    mutate(
      age_value = median(dat$age_value, na.rm = TRUE),
      Sex = make_ref_value(dat$Sex)
    )
  
  if (model_type == "full") {
    newdat <- newdat %>%
      mutate(
        edu2 = make_ref_value(dat$edu2),
        alc4 = make_ref_value(dat$alc4),
        ever_smoked = make_ref_value(dat$ever_smoked),
        PA_total_MET_quint = make_ref_value(dat$PA_total_MET_quint),
        Center_grp = make_ref_value(dat$Center_grp)
      )
  }
  
  pred_obj <- predict(fit, newdata = newdat, se.fit = TRUE)
  
  pred_fit <- as.numeric(pred_obj$fit)
  pred_se <- as.numeric(pred_obj$se.fit)
  
  newdat %>%
    mutate(
      pred = pred_fit,
      se = pred_se,
      lcl = pred_fit - 1.96 * pred_se,
      ucl = pred_fit + 1.96 * pred_se,
      n_obs = nrow(dat),
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
      model_type = model_type,
      model_label = case_when(
        model_type == "age_sex" ~ "Age + sex adjusted",
        model_type == "full" ~ "Fully adjusted"
      ),
      facet_label = paste0(.env$gene, " | ", .env$protein_col)
    )
}

# ------------------------------------------------------------
# 5. Create plot data
# ------------------------------------------------------------

sex_all_plot <- sex_all
# No filtering by sex interaction anymore.

plot_model_grid <- bind_rows(
  sex_all_plot %>% mutate(model_type = "age_sex"),
  sex_all_plot %>% mutate(model_type = "full")
)

plotdat_adjusted <- plot_model_grid %>%
  mutate(
    plot_data = pmap(
      list(
        protein_id,
        gene,
        pdi_type,
        selected_age_df,
        selected_diet_df,
        model_type
      ),
      ~ make_adjusted_plot_data_one_protein(
        data = df,
        protein_col = ..1,
        gene = ..2,
        pdi_type = ..3,
        selected_age_df = ..4,
        selected_diet_df = ..5,
        model_type = ..6
      )
    )
  ) %>%
  select(
    pdi_type,
    gene,
    protein_id,
    plot_data
  ) %>%
  unnest(plot_data)

plotdat_adjusted %>%
  distinct(facet_label, model_label, n_obs) %>%
  arrange(facet_label, model_label)

# ------------------------------------------------------------
# 6. Split hPDI and uPDI
# ------------------------------------------------------------

plotdat_hpdi_adjusted <- plotdat_adjusted %>%
  filter(pdi_type == "hPDI")

plotdat_updi_adjusted <- plotdat_adjusted %>%
  filter(pdi_type == "uPDI")

# ----------------
# SET COLORS LINES
#------------------
model_cols <- c(
  "Age + sex adjusted" = "#D27E64",
  "Fully adjusted"    = "#2B3256"
)

model_ltys <- c(
  "Age + sex adjusted" = "solid",
  "Fully adjusted"    = "dotted"
)

# ------------------------------------------------------------
# 7. Plot hPDI
# ------------------------------------------------------------

p_hpdi_adjusted <- ggplot(
  plotdat_hpdi_adjusted,
  aes(x = diet_value, y = pred, colour = model_label, fill = model_label)
) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.15, colour = NA) +
  geom_line(aes(linetype = model_label), linewidth = 0.8) +
  scale_colour_manual(values = model_cols) +
  scale_fill_manual(values = model_cols) +
  scale_linetype_manual(values = model_ltys) +
  facet_wrap(
    ~ facet_label,
    scales = "free_y",
    ncol = 3
  ) +
  labs(
    x = "Standardized hPDI",
    y = "Predicted protein abundance",
    colour = "Model",
    fill = "Model",
    linetype = "Model",
    title = "hPDI–protein associations: age/sex vs fully adjusted models"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

p_hpdi_adjusted

ggsave(
  filename = file.path(out_dir, "hPDI_age_sex_vs_fully_adjusted_A4.pdf"),
  plot = p_hpdi_adjusted,
  width = 8.27,
  height = 11.69,
  units = "in"
)

# ------------------------------------------------------------
# 8. Plot uPDI
# ------------------------------------------------------------

p_updi_adjusted <- ggplot(
  plotdat_updi_adjusted,
  aes(x = diet_value, y = pred, colour = model_label, fill = model_label)
) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.15, colour = NA) +
  geom_line(aes(linetype = model_label), linewidth = 0.8) +
  scale_colour_manual(values = model_cols) +
  scale_fill_manual(values = model_cols) +
  scale_linetype_manual(values = model_ltys) +
  facet_wrap(
    ~ facet_label,
    scales = "free_y",
    ncol = 3
  ) +
  labs(
    x = "Standardized uPDI",
    y = "Predicted protein abundance",
    colour = "Model",
    fill = "Model",
    linetype = "Model",
    title = "uPDI–protein associations: age/sex vs fully adjusted models"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

p_updi_adjusted

ggsave(
  filename = file.path(out_dir, "uPDI_age_sex_vs_fully_adjusted_A4.pdf"),
  plot = p_updi_adjusted,
  width = 8.27,
  height = 11.69,
  units = "in"
)
