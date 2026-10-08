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

##redraw forest plot and stitch with pathchwork
library(readxl)
library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(forcats)
library(grid)
# -----------------------------
# File paths
# -----------------------------
file_cancer <- file.path(MANUSCRIPT_DIR, "EPIC_cancer_survmodel_coefficients_limma_gengamma.xlsx")
file_t2d    <- file.path(MANUSCRIPT_DIR, "EPIC_T2D_survmodel_coefficients_limma_gengamma.xlsx")
file_cvd    <- file.path(MANUSCRIPT_DIR, "EPIC_CVD_survmodel_coefficients_limma_gengamma.xlsx")

res <- read_excel(file_cancer)

# ----------------------------
# Clean / standardize labels
# ----------------------------
res <- res %>%
  mutate(
    # standardize model names
    model = recode(
      model,
      "hPDI + Enet" = "hPDI + ENET",
      "uPDI + Enet" = "uPDI + ENET",
      "Enet only (hPDI set)" = "ENET only (hPDI set)",
      "Enet only (uPDI set)" = "ENET only (uPDI set)",
      "hPDI + GBM" = "hPDI + XGB",
      "uPDI + GBM" = "uPDI + XGB",
      "GBM only (hPDI set)" = "XGB only (hPDI set)",
      "GBM only (uPDI set)" = "XGB only (uPDI set)",
      .default = model
    ),
    
    # standardize term names robustly
    term_clean = case_when(
      term %in% c("Hpdi_NoVegoils_z", "hPDI (per 1 SD)") ~ "hPDI (per 1 SD)",
      term %in% c("Updi_NoVegoils_z", "uPDI (per 1 SD)") ~ "uPDI (per 1 SD)",
      
      term %in% c("hPDI_limma_unw_z", "hPDI_limma_unw", "limma protein signature for hPDI") ~
        "limma protein signature for hPDI",
      term %in% c("uPDI_limma_unw_z", "uPDI_limma_unw", "limma protein signature for uPDI") ~
        "limma protein signature for uPDI",
      
      term %in% c("Enet_Hpdi_z", "Enet_Hpdi", "ENET protein signature for hPDI") ~
        "ENET protein signature for hPDI",
      term %in% c("Enet_Updi_z", "Enet_Updi", "ENET protein signature for uPDI") ~
        "ENET protein signature for uPDI",
      
      term %in% c("GBMfull_Hpdi_z", "GBMfull_Hpdi", "GBM protein signature for hPDI", "XGB protein signature for hPDI") ~
        "XGB protein signature for hPDI",
      term %in% c("GBM_Updi_z", "GBM_Updi", "GBM protein signature for uPDI", "XGB protein signature for uPDI") ~
        "XGB protein signature for uPDI",
      
      TRUE ~ term
    )
  )

# Helper functions
pick_row <- function(data, model_name, term_name, row_label) {
  data %>%
    filter(model == model_name, term_clean == term_name) %>%
    slice(1) %>%
    transmute(
      est = beta,
      se = SE,
      TR = TR,
      TR_low = TR_low,
      TR_high = TR_high,
      row = row_label
    )
}

pick_row_right <- function(data, model_name, term_name, row_label, model_type_label) {
  data %>%
    filter(model == model_name, term_clean == term_name) %>%
    slice(1) %>%
    transmute(
      est = beta,
      se = SE,
      TR = TR,
      TR_low = TR_low,
      TR_high = TR_high,
      row = row_label,
      model_type = model_type_label
    )
}


# -----------------------------
# LEFT PANEL: diet score effects
# -----------------------------
left_rows <- c(
  "hPDI only",
  "hPDI + limma score",
  "hPDI + ENET score",
  "hPDI + XGB score",
  "pad_mid",
  "uPDI only",
  "uPDI + limma score",
  "uPDI + ENET score",
  "uPDI + XGB score"
)

left_combined_df <- bind_rows(
  pick_row(res, "hPDI only",     "hPDI (per 1 SD)", "hPDI only"),
  pick_row(res, "hPDI + limma",  "hPDI (per 1 SD)", "hPDI + limma score"),
  pick_row(res, "hPDI + ENET",   "hPDI (per 1 SD)", "hPDI + ENET score"),
  pick_row(res, "hPDI + XGB",    "hPDI (per 1 SD)", "hPDI + XGB score"),
  
  tibble(
    est = NA_real_, se = NA_real_, TR = NA_real_,
    TR_low = NA_real_, TR_high = NA_real_,
    row = "pad_mid"
  ),
  
  pick_row(res, "uPDI only",     "uPDI (per 1 SD)", "uPDI only"),
  pick_row(res, "uPDI + limma",  "uPDI (per 1 SD)", "uPDI + limma score"),
  pick_row(res, "uPDI + ENET",   "uPDI (per 1 SD)", "uPDI + ENET score"),
  pick_row(res, "uPDI + XGB",    "uPDI (per 1 SD)", "uPDI + XGB score")
) %>%
  mutate(row = factor(row, levels = rev(left_rows)))

# -----------------------------------
# RIGHT PANEL: protein score effects
# -----------------------------------
right_rows <- c(
  "hPDI dummy",
  "hPDI limma score",
  "hPDI ENET score",
  "hPDI XGB score",
  "pad_mid",
  "uPDI dummy",
  "uPDI limma score",
  "uPDI ENET score",
  "uPDI XGB score"
)

right_combined_df <- bind_rows(
  tibble(
    est = c(NA_real_, NA_real_),
    se = c(NA_real_, NA_real_),
    TR = c(NA_real_, NA_real_),
    TR_low = c(NA_real_, NA_real_),
    TR_high = c(NA_real_, NA_real_),
    row = "hPDI dummy",
    model_type = c("Independent", "Jointly with hPDI")
  ),
  
  pick_row_right(res, "limma only (hPDI set)", "limma protein signature for hPDI",
                 "hPDI limma score", "Independent"),
  pick_row_right(res, "hPDI + limma", "limma protein signature for hPDI",
                 "hPDI limma score", "Jointly with hPDI"),
  
  pick_row_right(res, "ENET only (hPDI set)", "ENET protein signature for hPDI",
                 "hPDI ENET score", "Independent"),
  pick_row_right(res, "hPDI + ENET", "ENET protein signature for hPDI",
                 "hPDI ENET score", "Jointly with hPDI"),
  
  pick_row_right(res, "XGB only (hPDI set)", "XGB protein signature for hPDI",
                 "hPDI XGB score", "Independent"),
  pick_row_right(res, "hPDI + XGB", "XGB protein signature for hPDI",
                 "hPDI XGB score", "Jointly with hPDI"),
  
  tibble(
    est = c(NA_real_, NA_real_),
    se = c(NA_real_, NA_real_),
    TR = c(NA_real_, NA_real_),
    TR_low = c(NA_real_, NA_real_),
    TR_high = c(NA_real_, NA_real_),
    row = "pad_mid",
    model_type = c("Independent", "Jointly with hPDI")
  ),
  
  tibble(
    est = c(NA_real_, NA_real_),
    se = c(NA_real_, NA_real_),
    TR = c(NA_real_, NA_real_),
    TR_low = c(NA_real_, NA_real_),
    TR_high = c(NA_real_, NA_real_),
    row = "uPDI dummy",
    model_type = c("Independent", "Jointly with uPDI")
  ),
  
  pick_row_right(res, "limma only (uPDI set)", "limma protein signature for uPDI",
                 "uPDI limma score", "Independent"),
  pick_row_right(res, "uPDI + limma", "limma protein signature for uPDI",
                 "uPDI limma score", "Jointly with uPDI"),
  
  pick_row_right(res, "ENET only (uPDI set)", "ENET protein signature for uPDI",
                 "uPDI ENET score", "Independent"),
  pick_row_right(res, "uPDI + ENET", "ENET protein signature for uPDI",
                 "uPDI ENET score", "Jointly with uPDI"),
  
  pick_row_right(res, "XGB only (uPDI set)", "XGB protein signature for uPDI",
                 "uPDI XGB score", "Independent"),
  pick_row_right(res, "uPDI + XGB", "XGB protein signature for uPDI",
                 "uPDI XGB score", "Jointly with uPDI")
) %>%
  mutate(
    row = factor(row, levels = rev(right_rows)),
    model_type = factor(
      model_type,
      levels = c("Independent", "Jointly with hPDI", "Jointly with uPDI")
    )
  )



## scaling x axis
## --------------
# Common x-axis limits for left and right panels
x_vals <- c(
  left_combined_df$TR_low,
  left_combined_df$TR_high,
  right_combined_df$TR_low,
  right_combined_df$TR_high
)

log_x_limits <- range(log10(x_vals), na.rm = TRUE)

# add a little padding
common_x_limits <- 10^(log_x_limits + c(-0.03, 0.03))

common_x_breaks <- c(0.6, 0.8, 1.0, 1.2, 1.5, 2.0)

##plotting
p_left_combined <- ggplot(
  left_combined_df,
  aes(x = TR, y = row, xmin = TR_low, xmax = TR_high)
) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_pointrange(na.rm = TRUE) +
  scale_x_log10(
    limits = common_x_limits,
    breaks = common_x_breaks
  ) +
  scale_y_discrete(labels = function(x) ifelse(grepl("^pad", x), "", x)) +
  labs(
    title = "Diet score estimates",
    x = "Time ratio for Cancer (per 1 SD higher diet score)",
    y = NULL
  ) +
  theme_bw()


pd <- position_dodge2(width = 0.5, reverse = TRUE)

p_right_combined <- ggplot(
  right_combined_df,
  aes(x = TR, y = row, xmin = TR_low, xmax = TR_high,
      colour = model_type, shape = model_type)
) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_pointrange(position = pd, na.rm = TRUE) +
  scale_x_log10(
    limits = common_x_limits,
    breaks = common_x_breaks
  )+
  scale_y_discrete(
    labels = c(
      "hPDI dummy" = "",
      "hPDI limma score" = "hPDI-limma",
      "hPDI ENET score" = "hPDI-ENET",
      "hPDI XGB score" = "hPDI-XGB",
      "pad_mid" = "",
      "uPDI dummy" = "",
      "uPDI limma score" = "uPDI-limma",
      "uPDI ENET score" = "uPDI-ENET",
      "uPDI XGB score" = "uPDI-XGB"
    )
  ) +
  scale_colour_manual(
    values = c(
      "Independent" = "#2B3256",
      "Jointly with hPDI" = "#D27E64",
      "Jointly with uPDI" = "#D27E64"
    )
  ) +
  labs(
    title = "Protein signature estimates",
    x = "Time ratio for Cancer (per 1 SD higher protein signature)",
    y = NULL,
    colour = NULL,
    shape = NULL
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

# ------------------
# Final combined plot
# ------------------
p_final_2panel <- p_left_combined + p_right_combined +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

p_final_2panel

#_____________________________
#T2d
#-----------------------------
res <- read_excel(file_t2d)

# -----------------------------
# Clean names robustly
# -----------------------------
res <- res %>%
  mutate(
    model_clean = recode(
      model,
      "hPDI + Enet" = "hPDI + ENET",
      "uPDI + Enet" = "uPDI + ENET",
      "Enet only (hPDI set)" = "ENET only (hPDI set)",
      "Enet only (uPDI set)" = "ENET only (uPDI set)",
      "hPDI + GBM" = "hPDI + XGB",
      "uPDI + GBM" = "uPDI + XGB",
      "GBM only (hPDI set)" = "XGB only (hPDI set)",
      "GBM only (uPDI set)" = "XGB only (uPDI set)",
      .default = model
    ),
    
    term_clean = case_when(
      term %in% c("Hpdi_NoVegoils_z", "hPDI (per 1 SD)") ~
        "hPDI (per 1 SD)",
      
      term %in% c("Updi_NoVegoils_z", "uPDI (per 1 SD)") ~
        "uPDI (per 1 SD)",
      
      term %in% c(
        "hPDI_limma_unw_z",
        "hPDI_limma_unw",
        "limma protein signature for hPDI",
        "limma protein signature for hPDI (per 1 SD)"
      ) ~
        "limma protein signature for hPDI",
      
      term %in% c(
        "uPDI_limma_unw_z",
        "uPDI_limma_unw",
        "limma protein signature for uPDI",
        "limma protein signature for uPDI (per 1 SD)"
      ) ~
        "limma protein signature for uPDI",
      
      term %in% c(
        "Enet_Hpdi_z",
        "Enet_Hpdi",
        "ENET protein signature for hPDI",
        "ENET protein signature for hPDI (per 1 SD)"
      ) ~
        "ENET protein signature for hPDI",
      
      term %in% c(
        "Enet_Updi_z",
        "Enet_Updi",
        "ENET protein signature for uPDI",
        "ENET protein signature for uPDI (per 1 SD)"
      ) ~
        "ENET protein signature for uPDI",
      
      term %in% c(
        "GBMfull_Hpdi_z",
        "GBMfull_Hpdi_z_sm",
        "GBMfull_Hpdi",
        "GBM protein signature for hPDI",
        "GBM protein signature for hPDI (per 1 SD)",
        "XGB protein signature for hPDI",
        "XGB protein signature for hPDI (per 1 SD)"
      ) ~
        "XGB protein signature for hPDI",
      
      term %in% c(
        "GBM_Updi_z",
        "GBM_Updi_z_sm",
        "GBM_Updi",
        "GBM protein signature for uPDI",
        "GBM protein signature for uPDI (per 1 SD)",
        "XGB protein signature for uPDI",
        "XGB protein signature for uPDI (per 1 SD)"
      ) ~
        "XGB protein signature for uPDI",
      
      TRUE ~ term
    )
  )

pick_row <- function(data, model_name, term_name, row_label) {
  out <- data %>%
    filter(model_clean == model_name, term_clean == term_name) %>%
    slice(1) %>%
    transmute(
      est = beta,
      se = SE,
      TR = TR,
      TR_low = TR_low,
      TR_high = TR_high,
      row = row_label
    )
  
  if (nrow(out) == 0) {
    warning(paste("Missing:", model_name, "|", term_name))
  }
  
  out
}

pick_row_right <- function(data, model_name, term_name, row_label, model_type_label) {
  out <- data %>%
    filter(model_clean == model_name, term_clean == term_name) %>%
    slice(1) %>%
    transmute(
      est = beta,
      se = SE,
      TR = TR,
      TR_low = TR_low,
      TR_high = TR_high,
      row = row_label,
      model_type = model_type_label
    )
  
  if (nrow(out) == 0) {
    warning(paste("Missing:", model_name, "|", term_name))
  }
  
  out
}

left_rows <- c(
  "hPDI only",
  "hPDI + limma score",
  "hPDI + ENET score",
  "hPDI + XGB score",
  "pad_mid",
  "uPDI only",
  "uPDI + limma score",
  "uPDI + ENET score",
  "uPDI + XGB score"
)

left_combined_df <- bind_rows(
  pick_row(res, "hPDI only",     "hPDI (per 1 SD)", "hPDI only"),
  pick_row(res, "hPDI + limma",  "hPDI (per 1 SD)", "hPDI + limma score"),
  pick_row(res, "hPDI + ENET",   "hPDI (per 1 SD)", "hPDI + ENET score"),
  pick_row(res, "hPDI + XGB",    "hPDI (per 1 SD)", "hPDI + XGB score"),
  
  tibble(
    est = NA_real_, se = NA_real_, TR = NA_real_,
    TR_low = NA_real_, TR_high = NA_real_,
    row = "pad_mid"
  ),
  
  pick_row(res, "uPDI only",     "uPDI (per 1 SD)", "uPDI only"),
  pick_row(res, "uPDI + limma",  "uPDI (per 1 SD)", "uPDI + limma score"),
  pick_row(res, "uPDI + ENET",   "uPDI (per 1 SD)", "uPDI + ENET score"),
  pick_row(res, "uPDI + XGB",    "uPDI (per 1 SD)", "uPDI + XGB score")
) %>%
  mutate(row = factor(row, levels = rev(left_rows)))

p_left_combined <- ggplot(
  left_combined_df,
  aes(x = TR, y = row, xmin = TR_low, xmax = TR_high)
) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_pointrange(na.rm = TRUE) +
  scale_x_log10() +
  scale_y_discrete(labels = function(x) ifelse(grepl("^pad", x), "", x)) +
  labs(
    title = "Diet score estimates",
    x = "Time ratio for T2D (per 1 SD higher diet score)",
    y = NULL
  ) +
  theme_bw()

right_rows <- c(
  "hPDI dummy",
  "hPDI limma score",
  "hPDI ENET score",
  "hPDI XGB score",
  "pad_mid",
  "uPDI dummy",
  "uPDI limma score",
  "uPDI ENET score",
  "uPDI XGB score"
)

right_combined_df <- bind_rows(
  tibble(
    est = c(NA_real_, NA_real_),
    se = c(NA_real_, NA_real_),
    TR = c(NA_real_, NA_real_),
    TR_low = c(NA_real_, NA_real_),
    TR_high = c(NA_real_, NA_real_),
    row = "hPDI dummy",
    model_type = c("Independent", "Jointly with hPDI")
  ),
  
  pick_row_right(res, "limma only (hPDI set)", "limma protein signature for hPDI",
                 "hPDI limma score", "Independent"),
  pick_row_right(res, "hPDI + limma", "limma protein signature for hPDI",
                 "hPDI limma score", "Jointly with hPDI"),
  
  pick_row_right(res, "ENET only (hPDI set)", "ENET protein signature for hPDI",
                 "hPDI ENET score", "Independent"),
  pick_row_right(res, "hPDI + ENET", "ENET protein signature for hPDI",
                 "hPDI ENET score", "Jointly with hPDI"),
  
  pick_row_right(res, "XGB only (hPDI set)", "XGB protein signature for hPDI",
                 "hPDI XGB score", "Independent"),
  pick_row_right(res, "hPDI + XGB", "XGB protein signature for hPDI",
                 "hPDI XGB score", "Jointly with hPDI"),
  
  tibble(
    est = c(NA_real_, NA_real_),
    se = c(NA_real_, NA_real_),
    TR = c(NA_real_, NA_real_),
    TR_low = c(NA_real_, NA_real_),
    TR_high = c(NA_real_, NA_real_),
    row = "pad_mid",
    model_type = c("Independent", "Jointly with hPDI")
  ),
  
  tibble(
    est = c(NA_real_, NA_real_),
    se = c(NA_real_, NA_real_),
    TR = c(NA_real_, NA_real_),
    TR_low = c(NA_real_, NA_real_),
    TR_high = c(NA_real_, NA_real_),
    row = "uPDI dummy",
    model_type = c("Independent", "Jointly with uPDI")
  ),
  
  pick_row_right(res, "limma only (uPDI set)", "limma protein signature for uPDI",
                 "uPDI limma score", "Independent"),
  pick_row_right(res, "uPDI + limma", "limma protein signature for uPDI",
                 "uPDI limma score", "Jointly with uPDI"),
  
  pick_row_right(res, "ENET only (uPDI set)", "ENET protein signature for uPDI",
                 "uPDI ENET score", "Independent"),
  pick_row_right(res, "uPDI + ENET", "ENET protein signature for uPDI",
                 "uPDI ENET score", "Jointly with uPDI"),
  
  pick_row_right(res, "XGB only (uPDI set)", "XGB protein signature for uPDI",
                 "uPDI XGB score", "Independent"),
  pick_row_right(res, "uPDI + XGB", "XGB protein signature for uPDI",
                 "uPDI XGB score", "Jointly with uPDI")
) %>%
  mutate(
    row = factor(row, levels = rev(right_rows)),
    model_type = factor(
      model_type,
      levels = c("Independent", "Jointly with hPDI", "Jointly with uPDI")
    )
  )

pd <- position_dodge2(width = 0.5, reverse = TRUE)

p_right_combined <- ggplot(
  right_combined_df,
  aes(
    x = TR, y = row, xmin = TR_low, xmax = TR_high,
    colour = model_type, shape = model_type
  )
) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_pointrange(position = pd, na.rm = TRUE) +
  scale_x_log10() +
  scale_x_log10(
    limits = common_x_limits,
    breaks = common_x_breaks
  ) +
  scale_y_discrete(
    labels = c(
      "hPDI dummy" = "",
      "hPDI limma score" = "hPDI-limma",
      "hPDI ENET score" = "hPDI-ENET",
      "hPDI XGB score" = "hPDI-XGB",
      "pad_mid" = "",
      "uPDI dummy" = "",
      "uPDI limma score" = "uPDI-limma",
      "uPDI ENET score" = "uPDI-ENET",
      "uPDI XGB score" = "uPDI-XGB"
    )
  ) +
  scale_colour_manual(
    values = c(
      "Independent" = "#2B3256",
      "Jointly with hPDI" = "#D27E64",
      "Jointly with uPDI" = "#D27E64"
    )
  ) +
  labs(
    title = "Protein signature estimates",
    x = "Time ratio for T2D (per 1 SD higher protein signature)",
    y = NULL,
    colour = NULL,
    shape = NULL
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

p_final_t2d_2panel <- p_left_combined + p_right_combined +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

p_final_t2d_2panel

#-----------------
#CVD
#----------------

res <- read_excel(file_cvd)

# -----------------------------
# Clean names robustly
# -----------------------------
res <- res %>%
  mutate(
    term_trim = trimws(term),
    
    model_clean = recode(
      model,
      "hPDI + Enet" = "hPDI + ENET",
      "uPDI + Enet" = "uPDI + ENET",
      "Enet only (hPDI set)" = "ENET only (hPDI set)",
      "Enet only (uPDI set)" = "ENET only (uPDI set)",
      "hPDI + GBM" = "hPDI + XGB",
      "uPDI + GBM" = "uPDI + XGB",
      "GBM only (hPDI set)" = "XGB only (hPDI set)",
      "GBM only (uPDI set)" = "XGB only (uPDI set)",
      .default = model
    ),
    
    term_clean = case_when(
      term_trim %in% c("Hpdi_NoVegoils_z", "hPDI (per 1 SD)") ~
        "hPDI (per 1 SD)",
      
      term_trim %in% c("Updi_NoVegoils_z", "uPDI (per 1 SD)") ~
        "uPDI (per 1 SD)",
      
      term_trim %in% c(
        "hPDI_limma_unw_z",
        "hPDI_limma_unw_z_sm",
        "hPDI_limma_unw",
        "limma protein signature for hPDI",
        "limma protein signature for hPDI (per 1 SD)"
      ) ~
        "limma protein signature for hPDI",
      
      term_trim %in% c(
        "uPDI_limma_unw_z",
        "uPDI_limma_unw_z_sm",
        "uPDI_limma_unw",
        "limma protein signature for uPDI",
        "limma protein signature for uPDI (per 1 SD)"
      ) ~
        "limma protein signature for uPDI",
      
      term_trim %in% c(
        "Enet_Hpdi_z",
        "Enet_Hpdi",
        "ENET protein signature for hPDI",
        "ENET protein signature for hPDI (per 1 SD)"
      ) ~
        "ENET protein signature for hPDI",
      
      term_trim %in% c(
        "Enet_Updi_z",
        "Enet_Updi",
        "ENET protein signature for uPDI",
        "ENET protein signature for uPDI (per 1 SD)"
      ) ~
        "ENET protein signature for uPDI",
      
      term_trim %in% c(
        "GBMfull_Hpdi_z",
        "GBMfull_Hpdi_z_sm",
        "GBMfull_Hpdi",
        "GBM protein signature for hPDI",
        "GBM protein signature for hPDI (per 1 SD)",
        "XGB protein signature for hPDI",
        "XGB protein signature for hPDI (per 1 SD)"
      ) ~
        "XGB protein signature for hPDI",
      
      term_trim %in% c(
        "GBM_Updi_z",
        "GBM_Updi_z_sm",
        "GBM_Updi",
        "GBM protein signature for uPDI",
        "GBM protein signature for uPDI (per 1 SD)",
        "XGB protein signature for uPDI",
        "XGB protein signature for uPDI (per 1 SD)"
      ) ~
        "XGB protein signature for uPDI",
      
      TRUE ~ term_trim
    )
  )

pick_row <- function(data, model_name, term_name, row_label) {
  out <- data %>%
    filter(model_clean == model_name, term_clean == term_name) %>%
    slice(1)
  
  if (nrow(out) == 0) {
    warning(paste("Missing:", model_name, "|", term_name))
    
    return(
      tibble(
        est = NA_real_,
        se = NA_real_,
        TR = NA_real_,
        TR_low = NA_real_,
        TR_high = NA_real_,
        row = row_label
      )
    )
  }
  
  out %>%
    transmute(
      est = beta,
      se = SE,
      TR = TR,
      TR_low = TR_low,
      TR_high = TR_high,
      row = row_label
    )
}

pick_row_right <- function(data, model_name, term_name, row_label, model_type_label) {
  out <- data %>%
    filter(model_clean == model_name, term_clean == term_name) %>%
    slice(1)
  
  if (nrow(out) == 0) {
    warning(paste("Missing:", model_name, "|", term_name))
    
    return(
      tibble(
        est = NA_real_,
        se = NA_real_,
        TR = NA_real_,
        TR_low = NA_real_,
        TR_high = NA_real_,
        row = row_label,
        model_type = model_type_label
      )
    )
  }
  
  out %>%
    transmute(
      est = beta,
      se = SE,
      TR = TR,
      TR_low = TR_low,
      TR_high = TR_high,
      row = row_label,
      model_type = model_type_label
    )
}

left_rows <- c(
  "hPDI only",
  "hPDI + limma score",
  "hPDI + ENET score",
  "hPDI + XGB score",
  "pad_mid",
  "uPDI only",
  "uPDI + limma score",
  "uPDI + ENET score",
  "uPDI + XGB score"
)

left_combined_df <- bind_rows(
  pick_row(res, "hPDI only",     "hPDI (per 1 SD)", "hPDI only"),
  pick_row(res, "hPDI + limma",  "hPDI (per 1 SD)", "hPDI + limma score"),
  pick_row(res, "hPDI + ENET",   "hPDI (per 1 SD)", "hPDI + ENET score"),
  pick_row(res, "hPDI + XGB",    "hPDI (per 1 SD)", "hPDI + XGB score"),
  
  tibble(
    est = NA_real_, se = NA_real_, TR = NA_real_,
    TR_low = NA_real_, TR_high = NA_real_,
    row = "pad_mid"
  ),
  
  pick_row(res, "uPDI only",     "uPDI (per 1 SD)", "uPDI only"),
  pick_row(res, "uPDI + limma",  "uPDI (per 1 SD)", "uPDI + limma score"),
  pick_row(res, "uPDI + ENET",   "uPDI (per 1 SD)", "uPDI + ENET score"),
  pick_row(res, "uPDI + XGB",    "uPDI (per 1 SD)", "uPDI + XGB score")
) %>%
  mutate(row = factor(row, levels = rev(left_rows)))

p_left_combined <- ggplot(
  left_combined_df,
  aes(x = TR, y = row, xmin = TR_low, xmax = TR_high)
) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_pointrange(na.rm = TRUE) +
  scale_x_log10() +
  scale_y_discrete(labels = function(x) ifelse(grepl("^pad", x), "", x)) +
  labs(
    title = "Diet score estimates",
    x = "Time ratio for CVD (per 1 SD higher diet score)",
    y = NULL
  ) +
  theme_bw()

right_rows <- c(
  "hPDI dummy",
  "hPDI limma score",
  "hPDI ENET score",
  "hPDI XGB score",
  "pad_mid",
  "uPDI dummy",
  "uPDI limma score",
  "uPDI ENET score",
  "uPDI XGB score"
)

right_combined_df <- bind_rows(
  tibble(
    est = c(NA_real_, NA_real_),
    se = c(NA_real_, NA_real_),
    TR = c(NA_real_, NA_real_),
    TR_low = c(NA_real_, NA_real_),
    TR_high = c(NA_real_, NA_real_),
    row = "hPDI dummy",
    model_type = c("Independent", "Jointly with hPDI")
  ),
  
  pick_row_right(res, "limma only (hPDI set)", "limma protein signature for hPDI",
                 "hPDI limma score", "Independent"),
  pick_row_right(res, "hPDI + limma", "limma protein signature for hPDI",
                 "hPDI limma score", "Jointly with hPDI"),
  
  pick_row_right(res, "ENET only (hPDI set)", "ENET protein signature for hPDI",
                 "hPDI ENET score", "Independent"),
  pick_row_right(res, "hPDI + ENET", "ENET protein signature for hPDI",
                 "hPDI ENET score", "Jointly with hPDI"),
  
  pick_row_right(res, "XGB only (hPDI set)", "XGB protein signature for hPDI",
                 "hPDI XGB score", "Independent"),
  pick_row_right(res, "hPDI + XGB", "XGB protein signature for hPDI",
                 "hPDI XGB score", "Jointly with hPDI"),
  
  tibble(
    est = c(NA_real_, NA_real_),
    se = c(NA_real_, NA_real_),
    TR = c(NA_real_, NA_real_),
    TR_low = c(NA_real_, NA_real_),
    TR_high = c(NA_real_, NA_real_),
    row = "pad_mid",
    model_type = c("Independent", "Jointly with hPDI")
  ),
  
  tibble(
    est = c(NA_real_, NA_real_),
    se = c(NA_real_, NA_real_),
    TR = c(NA_real_, NA_real_),
    TR_low = c(NA_real_, NA_real_),
    TR_high = c(NA_real_, NA_real_),
    row = "uPDI dummy",
    model_type = c("Independent", "Jointly with uPDI")
  ),
  
  pick_row_right(res, "limma only (uPDI set)", "limma protein signature for uPDI",
                 "uPDI limma score", "Independent"),
  pick_row_right(res, "uPDI + limma", "limma protein signature for uPDI",
                 "uPDI limma score", "Jointly with uPDI"),
  
  pick_row_right(res, "ENET only (uPDI set)", "ENET protein signature for uPDI",
                 "uPDI ENET score", "Independent"),
  pick_row_right(res, "uPDI + ENET", "ENET protein signature for uPDI",
                 "uPDI ENET score", "Jointly with uPDI"),
  
  pick_row_right(res, "XGB only (uPDI set)", "XGB protein signature for uPDI",
                 "uPDI XGB score", "Independent"),
  pick_row_right(res, "uPDI + XGB", "XGB protein signature for uPDI",
                 "uPDI XGB score", "Jointly with uPDI")
) %>%
  mutate(
    row = factor(row, levels = rev(right_rows)),
    model_type = factor(
      model_type,
      levels = c("Independent", "Jointly with hPDI", "Jointly with uPDI")
    )
  )

pd <- position_dodge2(width = 0.5, reverse = TRUE)

p_right_combined <- ggplot(
  right_combined_df,
  aes(
    x = TR, y = row, xmin = TR_low, xmax = TR_high,
    colour = model_type, shape = model_type
  )
) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_pointrange(position = pd, na.rm = TRUE) +
  scale_x_log10() +
  scale_y_discrete(
    labels = c(
      "hPDI dummy" = "",
      "hPDI limma score" = "hPDI-limma",
      "hPDI ENET score" = "hPDI-ENET",
      "hPDI XGB score" = "hPDI-XGB",
      "pad_mid" = "",
      "uPDI dummy" = "",
      "uPDI limma score" = "uPDI-limma",
      "uPDI ENET score" = "uPDI-ENET",
      "uPDI XGB score" = "uPDI-XGB"
    )
  ) +
  scale_colour_manual(
    values = c(
      "Independent" = "#2B3256",
      "Jointly with hPDI" = "#D27E64",
      "Jointly with uPDI" = "#D27E64"
    )
  ) +
  labs(
    title = "Protein signature estimates",
    x = "Time ratio for CVD (per 1 SD higher protein signature)",
    y = NULL,
    colour = NULL,
    shape = NULL
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

p_final_cvd_2panel <- p_left_combined + p_right_combined +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

p_final_cvd_2panel

##____________________________________________________
## some how prefer having them combined with same axes
## PUT all plots together
##____________________________________________________

# -----------------------------
# Function to read + clean each file
# -----------------------------
read_and_clean <- function(file, outcome_label) {
  read_excel(file) %>%
    mutate(
      outcome = outcome_label,
      term_trim = trimws(term),
      
      # harmonize model labels
      model_clean = recode(
        model,
        "hPDI + Enet" = "hPDI + ENET",
        "uPDI + Enet" = "uPDI + ENET",
        "Enet only (hPDI set)" = "ENET only (hPDI set)",
        "Enet only (uPDI set)" = "ENET only (uPDI set)",
        "hPDI + GBM" = "hPDI + XGB",
        "uPDI + GBM" = "uPDI + XGB",
        "GBM only (hPDI set)" = "XGB only (hPDI set)",
        "GBM only (uPDI set)" = "XGB only (uPDI set)",
        .default = model
      ),
      
      # harmonize term labels across Cancer / T2D / CVD
      term_clean = case_when(
        term_trim %in% c("Hpdi_NoVegoils_z", "hPDI (per 1 SD)") ~
          "hPDI (per 1 SD)",
        
        term_trim %in% c("Updi_NoVegoils_z", "uPDI (per 1 SD)") ~
          "uPDI (per 1 SD)",
        
        term_trim %in% c(
          "hPDI_limma_unw_z",
          "hPDI_limma_unw_z",
          "hPDI_limma_unw",
          "limma protein signature for hPDI",
          "limma protein signature for hPDI (per 1 SD)"
        ) ~ "limma protein signature for hPDI",
        
        term_trim %in% c(
          "uPDI_limma_unw_z",
          "uPDI_limma_unw_z",
          "uPDI_limma_unw",
          "limma protein signature for uPDI",
          "limma protein signature for uPDI (per 1 SD)"
        ) ~ "limma protein signature for uPDI",
        
        term_trim %in% c(
          "Enet_Hpdi_z",
          "Enet_Hpdi",
          "ENET protein signature for hPDI",
          "ENET protein signature for hPDI (per 1 SD)"
        ) ~ "ENET protein signature for hPDI",
        
        term_trim %in% c(
          "Enet_Updi_z",
          "Enet_Updi",
          "ENET protein signature for uPDI",
          "ENET protein signature for uPDI (per 1 SD)"
        ) ~ "ENET protein signature for uPDI",
        
        term_trim %in% c(
          "GBMfull_Hpdi_z",
          "GBMfull_Hpdi_z",
          "GBMfull_Hpdi",
          "GBM protein signature for hPDI",
          "GBM protein signature for hPDI (per 1 SD)",
          "XGB protein signature for hPDI",
          "XGB protein signature for hPDI (per 1 SD)"
        ) ~ "XGB protein signature for hPDI",
        
        term_trim %in% c(
          "GBM_Updi_z",
          "GBM_Updi_z",
          "GBM_Updi",
          "GBM protein signature for uPDI",
          "GBM protein signature for uPDI (per 1 SD)",
          "XGB protein signature for uPDI",
          "XGB protein signature for uPDI (per 1 SD)"
        ) ~ "XGB protein signature for uPDI",
        
        TRUE ~ term_trim
      )
    )
}

# -----------------------------
# Read all three
# -----------------------------
res_all <- bind_rows(
  read_and_clean(file_cancer, "Cancer"),
  read_and_clean(file_t2d,    "T2D"),
  read_and_clean(file_cvd,    "CVD")
) %>%
  mutate(
    outcome = factor(outcome, levels = c("Cancer", "T2D", "CVD"))
  )

# Optional check
res_all %>%
  count(outcome, model_clean, term_clean) %>%
  print(n = Inf)

# -----------------------------
# Safe pick helpers
# -----------------------------
pick_row <- function(data, outcome_name, model_name, term_name, row_label) {
  out <- data %>%
    filter(outcome == outcome_name,
           model_clean == model_name,
           term_clean == term_name) %>%
    slice(1)
  
  if (nrow(out) == 0) {
    warning(paste("Missing:", outcome_name, "|", model_name, "|", term_name))
    return(
      tibble(
        outcome = outcome_name,
        est = NA_real_, se = NA_real_,
        TR = NA_real_, TR_low = NA_real_, TR_high = NA_real_,
        row = row_label
      )
    )
  }
  
  out %>%
    transmute(
      outcome = outcome_name,
      est = beta,
      se = SE,
      TR = TR,
      TR_low = TR_low,
      TR_high = TR_high,
      row = row_label
    )
}

pick_row_right <- function(data, outcome_name, model_name, term_name, row_label, model_type_label) {
  out <- data %>%
    filter(outcome == outcome_name,
           model_clean == model_name,
           term_clean == term_name) %>%
    slice(1)
  
  if (nrow(out) == 0) {
    warning(paste("Missing:", outcome_name, "|", model_name, "|", term_name))
    return(
      tibble(
        outcome = outcome_name,
        est = NA_real_, se = NA_real_,
        TR = NA_real_, TR_low = NA_real_, TR_high = NA_real_,
        row = row_label,
        model_type = model_type_label
      )
    )
  }
  
  out %>%
    transmute(
      outcome = outcome_name,
      est = beta,
      se = SE,
      TR = TR,
      TR_low = TR_low,
      TR_high = TR_high,
      row = row_label,
      model_type = model_type_label
    )
}

# -----------------------------
# LEFT PANEL data
# -----------------------------
left_rows <- c(
  "hPDI only",
  "hPDI + limma score",
  "hPDI + ENET score",
  "hPDI + XGB score",
  "pad_mid",
  "uPDI only",
  "uPDI + limma score",
  "uPDI + ENET score",
  "uPDI + XGB score"
)

build_left_df <- function(outcome_name, data) {
  bind_rows(
    pick_row(data, outcome_name, "hPDI only",    "hPDI (per 1 SD)", "hPDI only"),
    pick_row(data, outcome_name, "hPDI + limma", "hPDI (per 1 SD)", "hPDI + limma score"),
    pick_row(data, outcome_name, "hPDI + ENET",  "hPDI (per 1 SD)", "hPDI + ENET score"),
    pick_row(data, outcome_name, "hPDI + XGB",   "hPDI (per 1 SD)", "hPDI + XGB score"),
    
    tibble(
      outcome = outcome_name,
      est = NA_real_, se = NA_real_,
      TR = NA_real_, TR_low = NA_real_, TR_high = NA_real_,
      row = "pad_mid"
    ),
    
    pick_row(data, outcome_name, "uPDI only",    "uPDI (per 1 SD)", "uPDI only"),
    pick_row(data, outcome_name, "uPDI + limma", "uPDI (per 1 SD)", "uPDI + limma score"),
    pick_row(data, outcome_name, "uPDI + ENET",  "uPDI (per 1 SD)", "uPDI + ENET score"),
    pick_row(data, outcome_name, "uPDI + XGB",   "uPDI (per 1 SD)", "uPDI + XGB score")
  )
}

left_combined_df <- bind_rows(
  build_left_df("Cancer", res_all),
  build_left_df("T2D", res_all),
  build_left_df("CVD", res_all)
) %>%
  mutate(
    outcome = factor(outcome, levels = c("Cancer", "T2D", "CVD")),
    row = factor(row, levels = rev(left_rows))
  )

# -----------------------------
# RIGHT PANEL data
# -----------------------------
right_rows <- c(
  "hPDI dummy",
  "hPDI limma score",
  "hPDI ENET score",
  "hPDI XGB score",
  "pad_mid",
  "uPDI dummy",
  "uPDI limma score",
  "uPDI ENET score",
  "uPDI XGB score"
)

build_right_df <- function(outcome_name, data) {
  bind_rows(
    tibble(
      outcome = outcome_name,
      est = c(NA_real_, NA_real_),
      se = c(NA_real_, NA_real_),
      TR = c(NA_real_, NA_real_),
      TR_low = c(NA_real_, NA_real_),
      TR_high = c(NA_real_, NA_real_),
      row = "hPDI dummy",
      model_type = c("Independent", "Jointly with hPDI")
    ),
    
    pick_row_right(data, outcome_name, "limma only (hPDI set)", "limma protein signature for hPDI",
                   "hPDI limma score", "Independent"),
    pick_row_right(data, outcome_name, "hPDI + limma", "limma protein signature for hPDI",
                   "hPDI limma score", "Jointly with hPDI"),
    
    pick_row_right(data, outcome_name, "ENET only (hPDI set)", "ENET protein signature for hPDI",
                   "hPDI ENET score", "Independent"),
    pick_row_right(data, outcome_name, "hPDI + ENET", "ENET protein signature for hPDI",
                   "hPDI ENET score", "Jointly with hPDI"),
    
    pick_row_right(data, outcome_name, "XGB only (hPDI set)", "XGB protein signature for hPDI",
                   "hPDI XGB score", "Independent"),
    pick_row_right(data, outcome_name, "hPDI + XGB", "XGB protein signature for hPDI",
                   "hPDI XGB score", "Jointly with hPDI"),
    
    tibble(
      outcome = outcome_name,
      est = c(NA_real_, NA_real_),
      se = c(NA_real_, NA_real_),
      TR = c(NA_real_, NA_real_),
      TR_low = c(NA_real_, NA_real_),
      TR_high = c(NA_real_, NA_real_),
      row = "pad_mid",
      model_type = c("Independent", "Jointly with hPDI")
    ),
    
    tibble(
      outcome = outcome_name,
      est = c(NA_real_, NA_real_),
      se = c(NA_real_, NA_real_),
      TR = c(NA_real_, NA_real_),
      TR_low = c(NA_real_, NA_real_),
      TR_high = c(NA_real_, NA_real_),
      row = "uPDI dummy",
      model_type = c("Independent", "Jointly with uPDI")
    ),
    
    pick_row_right(data, outcome_name, "limma only (uPDI set)", "limma protein signature for uPDI",
                   "uPDI limma score", "Independent"),
    pick_row_right(data, outcome_name, "uPDI + limma", "limma protein signature for uPDI",
                   "uPDI limma score", "Jointly with uPDI"),
    
    pick_row_right(data, outcome_name, "ENET only (uPDI set)", "ENET protein signature for uPDI",
                   "uPDI ENET score", "Independent"),
    pick_row_right(data, outcome_name, "uPDI + ENET", "ENET protein signature for uPDI",
                   "uPDI ENET score", "Jointly with uPDI"),
    
    pick_row_right(data, outcome_name, "XGB only (uPDI set)", "XGB protein signature for uPDI",
                   "uPDI XGB score", "Independent"),
    pick_row_right(data, outcome_name, "uPDI + XGB", "XGB protein signature for uPDI",
                   "uPDI XGB score", "Jointly with uPDI")
  )
}

right_combined_df <- bind_rows(
  build_right_df("Cancer", res_all),
  build_right_df("T2D", res_all),
  build_right_df("CVD", res_all)
) %>%
  mutate(
    outcome = factor(outcome, levels = c("Cancer", "T2D", "CVD")),
    row = factor(row, levels = rev(right_rows)),
    model_type = factor(
      model_type,
      levels = c("Independent", "Jointly with hPDI", "Jointly with uPDI")
    )
  )

# -----------------------------
# LEFT PLOT
# -----------------------------
common_x_limits <- c(0.85, 1.25)
common_x_breaks <- c(0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25)

x_scale_common <- scale_x_log10(
  limits = common_x_limits,
  breaks = common_x_breaks,
  labels = scales::number_format(accuracy = 0.01),
  expand = expansion(mult = c(0.02, 0.02))
)

common_theme <- theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", linewidth = 0.8),
    panel.spacing.y = unit(0.25, "cm"),
    axis.text.x = element_text(margin = margin(t = 0)),
    axis.ticks.length.x = unit(2, "pt"),
    axis.title.x = element_text(hjust = 0.5, margin = margin(t = 5)),
    plot.margin = margin(8, 8, 8, 8)
  )

disease_labs <- c(
  Cancer = "All cancer",
  T2D    = "Type 2 diabetes",
  CVD    = "Cardiovascular disease"
)

p_left_all <- ggplot(
  left_combined_df,
  aes(x = TR, y = row, xmin = TR_low, xmax = TR_high)
) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkred") +
  geom_pointrange(na.rm = TRUE) +
  x_scale_common +
  scale_y_discrete(labels = function(x) ifelse(grepl("^pad", x), "", x)) +
  facet_grid(
    outcome ~ .,
    switch = "y",
    labeller = labeller(outcome = disease_labs)
  ) +
  labs(
    title = "Diet score estimates",
    x = "Time ratio (per 1 SD higher diet score)",
    y = NULL
  ) +
  common_theme +
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill = "gray75", colour = "gray75", linewidth = 0.8),
    strip.text.y.left = element_text(angle = 90, face = "bold")
  )

# -----------------------------
# RIGHT PLOT
# -----------------------------
pd <- position_dodge2(width = 0.5, reverse = TRUE)

p_right_all <- ggplot(
  right_combined_df,
  aes(
    x = TR, y = row, xmin = TR_low, xmax = TR_high,
    colour = model_type, shape = model_type
  )
) +
  geom_vline(xintercept = 1, linetype = "dashed",  color = "darkred") +
  geom_pointrange(position = pd, na.rm = TRUE) +
  x_scale_common +
  scale_y_discrete(
    labels = c(
      "hPDI dummy" = "",
      "hPDI limma score" = "hPDI-limma",
      "hPDI ENET score" = "hPDI-ENET",
      "hPDI XGB score" = "hPDI-XGB",
      "pad_mid" = "",
      "uPDI dummy" = "",
      "uPDI limma score" = "uPDI-limma",
      "uPDI ENET score" = "uPDI-ENET",
      "uPDI XGB score" = "uPDI-XGB"
    )
  ) +
  facet_grid(outcome ~ .) +
  scale_colour_manual(
    values = c(
      "Independent" = "#2B3256",
      "Jointly with hPDI" = "#D27E64",
      "Jointly with uPDI" = "#D27E64"
    )
  ) +
  labs(
    title = "Protein signature estimates",
    x = "Time ratio (per 1 SD higher protein signature)",
    y = NULL,
    colour = NULL,
    shape = NULL
  ) +
  common_theme +
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill = "gray75", colour = "gray75", linewidth = 0.8),
    strip.text.y.left = element_text(angle = 90, face = "bold")
  )
# -----------------------------
# Final combined figure
# -----------------------------
p_final_combined <- p_left_all + p_right_all +
  plot_layout(ncol = 2, guides = "collect", widths = c(1, 1)) &
  theme(legend.position = "bottom")

p_final_combined

# -----------------------------
# Save
# -----------------------------
ggsave(
  file.path(MANUSCRIPT_DIR, "EPIC_survival_combined_Cancer_T2D_CVD_XGB.png"),
  plot = p_final_combined,
  width = 13,
  height = 16,
  units = "in",
  dpi = 600
)
