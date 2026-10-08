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
library(tibble)
library(openxlsx)
library(writexl)

out_dir <- MANUSCRIPT_DIR

# ------------------------------------------------------------
# 1. Helper: convert evidence text to short label
# ------------------------------------------------------------

short_interaction_label <- function(x) {
  case_when(
    grepl("^Strong", x) ~ "Strong",
    grepl("^Nominal", x) ~ "Nominal",
    TRUE ~ "No"
  )
}

# ------------------------------------------------------------
# 2. Manual row order, exactly as in your table
# ------------------------------------------------------------

table_order <- tribble(
  ~Diet,  ~Protein,    ~`Protein sequence ID`,
  "hPDI", "CTSS",      "seq.3181.50",
  "hPDI", "DNPH1",     "seq.17785.11",
  "hPDI", "PTK7",      "seq.9525.1",
  "hPDI", "GREM1",     "seq.18878.15",
  "hPDI", "CCM2",      "seq.12347.29",
  "hPDI", "POSTN",     "seq.6645.53",
  "hPDI", "CA11",      "seq.8859.51",
  "hPDI", "MMP10",     "seq.10479.18",
  "hPDI", "MMP10",     "seq.8479.4",
  "hPDI", "GSTO1",     "seq.12436.84",
  "hPDI", "LRIG1",     "seq.18831.6",
  "hPDI", "MANBA",     "seq.6382.17",
  "hPDI", "NCAN",      "seq.15573.110",
  "hPDI", "PARVA",     "seq.13434.172",
  "hPDI", "PCSK9",     "seq.5231.79",
  "hPDI", "PCSK9",     "seq.8246.9",
  "hPDI", "PEAR1",     "seq.8275.31",
  "hPDI", "SPP1",      "seq.13113.7",
  
  "uPDI", "KRT5",      "seq.11177.16",
  "uPDI", "EGFR",      "seq.2677.1",
  "uPDI", "CPNE1",     "seq.11602.12",
  "uPDI", "PLA2G2A",   "seq.2692.74",
  "uPDI", "SERPINA3",  "seq.2879.9",
  "uPDI", "HHIP",      "seq.10833.64",
  "uPDI", "CCM2",      "seq.12347.29",
  "uPDI", "ACE",       "seq.10714.7",
  "uPDI", "LCT",       "seq.9017.58",
  "uPDI", "PTPN9",     "seq.12633.3",
  "uPDI", "POSTN",     "seq.3457.57",
  "uPDI", "POSTN",     "seq.6645.53",
  "uPDI", "MMP10",     "seq.10479.18",
  "uPDI", "MMP10",     "seq.8479.4",
  "uPDI", "PCSK9",     "seq.8246.9",
  "uPDI", "ACADVL",    "seq.19262.219",
  "uPDI", "NADK",      "seq.13624.17",
  "uPDI", "SPP1",      "seq.13113.7"
) %>%
  mutate(row_order = row_number())

# ------------------------------------------------------------
# 3. Read interaction files and make lookup tables
# ------------------------------------------------------------

# Sex
sex_hpdi <- read.xlsx(
  file.path(out_dir, "sex_interaction_hpdi_correct_direction_ageHpdiDFplugged.xlsx"),
  sheet = "sex_interaction"
) %>%
  transmute(
    Diet = "hPDI",
    `Protein sequence ID` = protein_id,
    `Sex interaction` = short_interaction_label(evidence_sex_interaction)
  )

sex_updi <- read.xlsx(
  file.path(out_dir, "sex_interaction_updi_correct_direction_ageUpdiDFplugged.xlsx"),
  sheet = "sex_interaction"
) %>%
  transmute(
    Diet = "uPDI",
    `Protein sequence ID` = protein_id,
    `Sex interaction` = short_interaction_label(evidence_sex_interaction)
  )

sex_lookup <- bind_rows(sex_hpdi, sex_updi)

# Smoking
smoking_hpdi <- read.xlsx(
  file.path(out_dir, "smoking_interaction_hpdi_correct_direction_ageHpdiDFplugged.xlsx"),
  sheet = "smoking_interaction"
) %>%
  transmute(
    Diet = "hPDI",
    `Protein sequence ID` = protein_id,
    `Smoking interaction` = short_interaction_label(evidence_smoking_interaction)
  )

smoking_updi <- read.xlsx(
  file.path(out_dir, "smoking_interaction_updi_correct_direction_ageUpdiDFplugged.xlsx"),
  sheet = "smoking_interaction"
) %>%
  transmute(
    Diet = "uPDI",
    `Protein sequence ID` = protein_id,
    `Smoking interaction` = short_interaction_label(evidence_smoking_interaction)
  )

smoking_lookup <- bind_rows(smoking_hpdi, smoking_updi)

# Alcohol
alcohol_hpdi <- read.xlsx(
  file.path(out_dir, "alcohol_interaction_hpdi_correct_direction_ageHpdiDFplugged.xlsx"),
  sheet = "alcohol_interaction"
) %>%
  transmute(
    Diet = "hPDI",
    `Protein sequence ID` = protein_id,
    `Alcohol consumption interaction` = short_interaction_label(evidence_alcohol_interaction)
  )

alcohol_updi <- read.xlsx(
  file.path(out_dir, "alcohol_interaction_updi_correct_direction_ageUpdiDFplugged.xlsx"),
  sheet = "alcohol_interaction"
) %>%
  transmute(
    Diet = "uPDI",
    `Protein sequence ID` = protein_id,
    `Alcohol consumption interaction` = short_interaction_label(evidence_alcohol_interaction)
  )

alcohol_lookup <- bind_rows(alcohol_hpdi, alcohol_updi)

# Age
age_hpdi <- read.xlsx(
  file.path(out_dir, "age_interaction_hpdi_correct_direction_ageHpdiDFplugged.xlsx"),
  sheet = "age_interaction"
) %>%
  transmute(
    Diet = "hPDI",
    `Protein sequence ID` = protein_id,
    `Age interaction` = short_interaction_label(evidence_age_interaction)
  )

age_updi <- read.xlsx(
  file.path(out_dir, "age_interaction_updi_correct_direction_ageUpdiDFplugged.xlsx"),
  sheet = "age_interaction"
) %>%
  transmute(
    Diet = "uPDI",
    `Protein sequence ID` = protein_id,
    `Age interaction` = short_interaction_label(evidence_age_interaction)
  )

age_lookup <- bind_rows(age_hpdi, age_updi)

# ------------------------------------------------------------
# 4. Join everything in the exact table order
# ------------------------------------------------------------

interaction_cols_ordered <- table_order %>%
  left_join(sex_lookup, by = c("Diet", "Protein sequence ID")) %>%
  left_join(smoking_lookup, by = c("Diet", "Protein sequence ID")) %>%
  left_join(alcohol_lookup, by = c("Diet", "Protein sequence ID")) %>%
  left_join(age_lookup, by = c("Diet", "Protein sequence ID")) %>%
  arrange(row_order) %>%
  mutate(
    `Sex interaction` = coalesce(`Sex interaction`, "No"),
    `Smoking interaction` = coalesce(`Smoking interaction`, "No"),
    `Alcohol consumption interaction` = coalesce(`Alcohol consumption interaction`, "No"),
    `Age interaction` = coalesce(`Age interaction`, "No")
  )

# ------------------------------------------------------------
# 5. Version to paste directly into Excel
# ------------------------------------------------------------

paste_ready <- interaction_cols_ordered %>%
  select(
    `Diet`,
    `Protein sequence ID`,
    `Sex interaction`,
    `Smoking interaction`,
    `Alcohol consumption interaction`,
    `Age interaction`
  )

if (interactive()) View(paste_ready)

# Copy to clipboard for Excel paste
write.table(
  paste_ready,
  file = "clipboard",
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# Also save a checked version with identifiers
write.xlsx(
  list(
    paste_ready = paste_ready,
    with_identifiers = interaction_cols_ordered
  ),
  file = file.path(out_dir, "interaction_strength_columns_for_excel_paste.xlsx"),
  overwrite = TRUE
)
