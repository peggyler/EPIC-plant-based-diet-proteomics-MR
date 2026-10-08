# Run from the repository root in a fresh R session. See README.md.
# EDIT LOCAL PATHS HERE. No study data or results are supplied.
PROJECT_ROOT <- "."
DATA_DIR <- file.path(PROJECT_ROOT, "data")
OUTPUT_ROOT <- file.path(PROJECT_ROOT, "outputs")
MANUSCRIPT_DIR <- file.path(OUTPUT_ROOT, "manuscript")
REFERENCE_DIR <- file.path(DATA_DIR, "reference")
SIGNATURE_SCRIPTS <- file.path(PROJECT_ROOT, "signatures")
CONTEXT_SCRIPTS <- file.path(PROJECT_ROOT, "contextualisation")
HPDI_RUN <- "staged_run_20260424_0935"  # choose the manuscript-producing run/date
LIMMA_DATE <- "20260424"  # choose the manuscript-producing run/date
dir.create(MANUSCRIPT_DIR, recursive = TRUE, showWarnings = FALSE)

#'present protein list limma, gbm, and enet for supplement tables
#'merge with md to get uniprot and genename

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(janitor)
})

proj_root <- PROJECT_ROOT

sig_dir <- file.path(OUTPUT_ROOT, "hpdi"
)

gbm_file <- file.path(
  sig_dir,
  file.path("gbm", HPDI_RUN, "signature_SHAP.csv")
)

elasticnet_file <- file.path(
  sig_dir,
  "elasticnet_ENet_Hpdi_NoVegoils_filtFDR5_lamMIN_coef_full.csv"
)

limma_file <- file.path(
  sig_dir,
  paste0("Hpidnovegoils_limma_", LIMMA_DATE, ".csv")
)

# Metadata file from previous workflow. Adjust this path if your md file is elsewhere.
md_file <- file.path(REFERENCE_DIR, "MetaData_Soma_SMP.rds"
)

out_dir <- file.path(sig_dir, "HPDI_compiled_cis_pqtl")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(gbm_file))
stopifnot(file.exists(elasticnet_file))
stopifnot(file.exists(limma_file))
stopifnot(file.exists(md_file))

message("GBM file:        ", gbm_file)
message("Elastic net file:", elasticnet_file)
message("limma file:      ", limma_file)
message("Metadata file:   ", md_file)

md <- readRDS(md_file)

md_map <- md %>%
  select(
    protein_id = AptName,
    seq_id = SeqId,
    target = Target,
    gene = EntrezGeneSymbol,
    uniprot = UniProt,
    type = Type,
    dilution = Dilution,
    organism = Organism
  ) %>%
  mutate(
    seq_id_decode = str_replace_all(seq_id, "-", "_"),
    uniprot = if_else(is.na(uniprot), NA_character_, sub(";.*$", "", uniprot))
  ) %>%
  distinct(protein_id, .keep_all = TRUE)

# Check metadata map
md_map %>%
  summarise(
    n_rows = n(),
    n_protein_id = n_distinct(protein_id),
    n_seq_id = n_distinct(seq_id),
    n_missing_gene = sum(is.na(gene) | gene == ""),
    n_missing_uniprot = sum(is.na(uniprot) | uniprot == "")
  ) %>%
  print()

limma_raw <- fread(limma_file) %>%
  as_tibble() %>%
  clean_names()

limma_list <- limma_raw %>%
  transmute(
    protein_id = protein,
    limma = 1L,
    limma_log_fc = log_fc,
    limma_ave_expr = ave_expr,
    limma_t = t,
    limma_p_value = p_value,
    limma_adj_p_val = adj_p_val,
    limma_b = b,
    limma_q_bh = q_bh,
    limma_q_by = q_by,
    limma_p_bonf = p_bonf
  ) %>%
  left_join(md_map, by = "protein_id")

#reorder columns
df_hpdi <- limma_list %>%
  select(gene, protein_id, target, uniprot, limma_log_fc, limma_p_value, limma_adj_p_val,everything())

fwrite(
  df_hpdi,
  file.path(MANUSCRIPT_DIR, "limmalist_hpdi_forsupp.csv")
)

## FOR UPDI
## --------
proj_root <- PROJECT_ROOT

sig_dir <- file.path(OUTPUT_ROOT, "updi"
)

limma_file <- file.path(
  sig_dir,
  paste0("Updinovegoils_limma_", LIMMA_DATE, ".csv")
)

limma_raw <- fread(limma_file) %>%
  as_tibble() %>%
  clean_names()

limma_list <- limma_raw %>%
  transmute(
    protein_id = protein,
    limma = 1L,
    limma_log_fc = log_fc,
    limma_ave_expr = ave_expr,
    limma_t = t,
    limma_p_value = p_value,
    limma_adj_p_val = adj_p_val,
    limma_b = b,
    limma_q_bh = q_bh,
    limma_q_by = q_by,
    limma_p_bonf = p_bonf
  ) %>%
  left_join(md_map, by = "protein_id")

#reorder columns
df_updi <- limma_list %>%
  select(gene, protein_id, target, uniprot, limma_log_fc, limma_p_value, limma_adj_p_val,everything())

fwrite(
  df_updi,
  file.path(MANUSCRIPT_DIR, "limmalist_updi_forsupp.csv")
)
