# Run from the repository root in a fresh R session. See README.md.
# CHANGED: discovery-cohort correlations now use df_discovery consistently. Rerun and compare.
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

# Load discovery cohort
load(file.path(DATA_DIR, "subcohort_20251106.RData"))
# object: subcohort

# Load full cohort with scores
load(file.path(DATA_DIR, "score_tte_analysis_20251126.RData"))
# object: df

wgts <- readRDS(file.path(DATA_DIR, "Weights_Feb2025.rds"))

df <- df1

wgts <- dplyr::rename(wgts, Idepic = idepic)

df <- df %>%
  left_join(wgts, by = "Idepic") %>%
  filter(!is.na(Weights))

# Keep only rows in df whose IDEPIC is present in subcohort
df_discovery <- df %>%
  semi_join(subcohort %>% distinct(Idepic), by = "Idepic")

# Check numbers
nrow(subcohort)
nrow(df_discovery)
length(intersect(df$Idepic, subcohort$Idepic))

df_discovery$Hpdi_NoVegoils_z  <- as.numeric(scale(df_discovery$Hpdi_NoVegoils))
df_discovery$Updi_NoVegoils_z  <- as.numeric(scale(df_discovery$Updi_NoVegoils))


library(purrr)
library(tibble)

r2_table <- tibble(
  score_var = c("hPDI_limma_unw_z", "uPDI_limma_unw_z"),
  diet_var  = c("Hpdi_NoVegoils_z", "Updi_NoVegoils_z")
) %>%
  mutate(
    n = map2_int(score_var, diet_var, ~sum(complete.cases(df_discovery[, c(.x, .y)]))),
    spearman_r = map2_dbl(score_var, diet_var, ~cor(
      df_discovery[[.x]],
      df_discovery[[.y]],
      method = "spearman",
      use = "complete.obs"
    )),
    spearman_r2 = spearman_r^2
  )

r2_table
# score_var        diet_var             n spearman_r spearman_r2
# <chr>            <chr>            <int>      <dbl>       <dbl>
#   1 hPDI_limma_unw_z Hpdi_NoVegoils_z  4372     0.0706     0.00498
# 2 uPDI_limma_unw_z Updi_NoVegoils_z  4372     0.105      0.0110 

r2_table <- tibble(
  score_var = c("hPDI_limma_z", "uPDI_limma_z"),
  diet_var  = c("Hpdi_NoVegoils_z", "Updi_NoVegoils_z")
) %>%
  mutate(
    n = map2_int(score_var, diet_var, ~sum(complete.cases(df_discovery[, c(.x, .y)]))),
    spearman_r = map2_dbl(score_var, diet_var, ~cor(
      df_discovery[[.x]],
      df_discovery[[.y]],
      method = "spearman",
      use = "complete.obs"
    )),
    spearman_r2 = spearman_r^2
  )

r2_table
# score_var    diet_var             n spearman_r spearman_r2
# <chr>        <chr>            <int>      <dbl>       <dbl>
#   1 hPDI_limma_z Hpdi_NoVegoils_z  4372     0.0739     0.00546
# 2 uPDI_limma_z Updi_NoVegoils_z  4372     0.113      0.0127
