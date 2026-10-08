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

#'---- Compile proteins from GBM, elastic net, and limma signatures ----
#' Project...: ProPlant / HPDI no vegetable oils
#' Created...: Peggy Ler
#' Purpose...: Compile signature proteins from three discovery methods,
#'             annotate them with SomaScan metadata, and create an input file
#'             for cis-pQTL extraction from deCODE ST08.
#'
#' Main outputs:
#'   1. compiled_signature_proteins_by_method.csv
#'   2. compiled_signature_proteins_by_method.rds
#'
#' Notes:
#'   - limma proteins are filtered at FDR < 0.05.
#'   - GBM and elastic net inputs are assumed to already contain selected proteins.
#'   - seq_id is harmonised from SomaScan AptName format, e.g. seq.24954.83 -> 24954_83.
#'   - This script does not perform cis filtering. deCODE ST08 is already a cis-pQTL table.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(janitor)
})

# -------------------------------------------------------------------------
# 0. Paths
# -------------------------------------------------------------------------

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

# -------------------------------------------------------------------------
# Read and standardise SomaScan metadata
# -------------------------------------------------------------------------

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

# n_rows n_protein_id n_seq_id n_missing_gene n_missing_uniprot
# <int>        <int>    <int>          <int>             <int>
#   1   7596         7596     7596            274                34
#missing are mouse targets

# -------------------------------------------------------------------------
# 3. Read GBM signature
# -------------------------------------------------------------------------
gbm_raw <- fread(gbm_file) %>%
  as_tibble() %>%
  clean_names()

# Expected columns include:
# feature, weight_abs, weight_signed, mean_shap_train, mean_abs_train,
# mean_abs_test, direction, cumweight

gbm_sig <- gbm_raw %>%
  transmute(
    protein_id = feature,
    gbm = 1L,
    gbm_weight_abs = weight_abs,
    gbm_weight_signed = weight_signed,
    gbm_mean_shap_train = mean_shap_train,
    gbm_mean_abs_train = mean_abs_train,
    gbm_mean_abs_test = mean_abs_test,
    gbm_direction = direction,
    gbm_cumweight = cumweight
  ) %>%
  left_join(md_map, by = "protein_id")

# Check GBM mapping
cat("\nGBM mapping check:\n")
gbm_sig %>%
  summarise(
    n_rows = n(),
    n_proteins = n_distinct(protein_id),
    n_missing_seq_id = sum(is.na(seq_id)),
    n_missing_gene = sum(is.na(gene) | gene == ""),
    n_missing_uniprot = sum(is.na(uniprot) | uniprot == "")
  ) %>%
  print()

# n_rows n_proteins n_missing_seq_id n_missing_gene n_missing_uniprot
# <int>      <int>            <int>          <int>             <int>
#   1    667        667                0              0                 0

fwrite(
  gbm_sig,
  file.path(MANUSCRIPT_DIR, "Hpdi_GBM_proteinlist.csv")
)

# -------------------------------------------------------------------------
# 4. Read elastic net signature
# -------------------------------------------------------------------------
enet_raw <- fread(elasticnet_file) %>%
  as_tibble() %>%
  clean_names()

# Expected columns include:
# feature, coef, abs_coef, feature_type
# Keep only protein features and non-zero coefficients.

enet_sig <- enet_raw %>%
  filter(feature_type == "protein") %>%
  filter(!is.na(coef), coef != 0) %>%
  transmute(
    protein_id = feature,
    elasticnet = 1L,
    elasticnet_coef = coef,
    elasticnet_abs_coef = abs_coef,
    elasticnet_feature_type = feature_type
  ) %>%
  left_join(md_map, by = "protein_id")

# Check elastic net mapping
cat("\nElastic net mapping check:\n")
enet_sig %>%
  summarise(
    n_rows = n(),
    n_proteins = n_distinct(protein_id),
    n_missing_seq_id = sum(is.na(seq_id)),
    n_missing_gene = sum(is.na(gene) | gene == ""),
    n_missing_uniprot = sum(is.na(uniprot) | uniprot == "")
  ) %>%
  print()
# n_rows n_proteins n_missing_seq_id n_missing_gene n_missing_uniprot
# <int>      <int>            <int>          <int>             <int>
#   1    108        108                0              0                 0

fwrite(
  enet_sig,
  file.path(MANUSCRIPT_DIR, "Hpdi_ENET_proteinlist.csv")
)

# -------------------------------------------------------------------------
# 5. Read limma signature and filter FDR < 0.05
# -------------------------------------------------------------------------

limma_raw <- fread(limma_file) %>%
  as_tibble() %>%
  clean_names()

# Keep FDR < 0.05 using adj_p_val.

limma_sig <- limma_raw %>%
  filter(!is.na(adj_p_val), adj_p_val < 0.05) %>%
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

# Check limma mapping
cat("\nlimma mapping check:\n")
limma_sig %>%
  summarise(
    n_rows = n(),
    n_proteins = n_distinct(protein_id),
    n_missing_seq_id = sum(is.na(seq_id)),
    n_missing_gene = sum(is.na(gene) | gene == ""),
    n_missing_uniprot = sum(is.na(uniprot) | uniprot == "")
  ) %>%
  print()
# n_rows n_proteins n_missing_seq_id n_missing_gene n_missing_uniprot
# <int>      <int>            <int>          <int>             <int>
#   1    619        619                0              0                 0

fwrite(
  limma_sig,
  file.path(MANUSCRIPT_DIR, "Hpdi_LIMMA_proteinlist.csv")
)

# -------------------------------------------------------------------------
# 6. Compile all proteins and mark source method
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Keep protein information + signature marker
# -------------------------------------------------------------------------

gbm_keep <- gbm_sig %>%
  select(
    protein_id,
    seq_id,
    seq_id_decode,
    gene,
    target,
    uniprot
  ) %>%
  distinct() %>%
  mutate(
    gbm = 1L,
    elasticnet = 0L,
    limma = 0L
  )

enet_keep <- enet_sig %>%
  select(
    protein_id,
    seq_id,
    seq_id_decode,
    gene,
    target,
    uniprot
  ) %>%
  distinct() %>%
  mutate(
    gbm = 0L,
    elasticnet = 1L,
    limma = 0L
  )

limma_keep <- limma_sig %>%
  select(
    protein_id,
    seq_id,
    seq_id_decode,
    gene,
    target,
    uniprot
  ) %>%
  distinct() %>%
  mutate(
    gbm = 0L,
    elasticnet = 0L,
    limma = 1L
  )

# -------------------------------------------------------------------------
# Compile into one protein table
# -------------------------------------------------------------------------

compiled_proteins <- bind_rows(
  gbm_keep,
  enet_keep,
  limma_keep
) %>%
  group_by(protein_id) %>%
  summarise(
    seq_id = first(na.omit(seq_id)),
    seq_id_decode = first(na.omit(seq_id_decode)),
    gene = first(na.omit(gene)),
    target = first(na.omit(target)),
    uniprot = first(na.omit(uniprot)),
    
    gbm = max(gbm, na.rm = TRUE),
    elasticnet = max(elasticnet, na.rm = TRUE),
    limma = max(limma, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    n_methods = gbm + elasticnet + limma,
    source_methods = case_when(
      gbm == 1 & elasticnet == 1 & limma == 1 ~ "GBM;elasticnet;limma",
      gbm == 1 & elasticnet == 1 & limma == 0 ~ "GBM;elasticnet",
      gbm == 1 & elasticnet == 0 & limma == 1 ~ "GBM;limma",
      gbm == 0 & elasticnet == 1 & limma == 1 ~ "elasticnet;limma",
      gbm == 1 & elasticnet == 0 & limma == 0 ~ "GBM",
      gbm == 0 & elasticnet == 1 & limma == 0 ~ "elasticnet",
      gbm == 0 & elasticnet == 0 & limma == 1 ~ "limma",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(desc(n_methods), protein_id)

# -------------------------------------------------------------------------
# QC
# -------------------------------------------------------------------------

compiled_proteins %>%
  summarise(
    n_total_proteins = n(),
    n_gbm = sum(gbm == 1),
    n_elasticnet = sum(elasticnet == 1),
    n_limma = sum(limma == 1),
    n_in_1_method = sum(n_methods == 1),
    n_in_2_methods = sum(n_methods == 2),
    n_in_3_methods = sum(n_methods == 3),
    n_missing_seq_id_decode = sum(is.na(seq_id_decode)),
    n_missing_gene = sum(is.na(gene)),
    n_missing_uniprot = sum(is.na(uniprot)),
    n_duplicate_protein_id = sum(duplicated(protein_id))
  )

# A tibble: 1 × 11
# n_total_proteins n_gbm n_elasticnet n_limma n_in_1_method n_in_2_methods n_in_3_methods n_missing_seq_id_decode
# <int> <int>        <int>   <int>         <int>          <int>          <int>                   <int>
#    1206          667      108        619          1055            114             37                       0
# # ℹ 3 more variables: n_missing_gene <int>, n_missing_uniprot <int>, n_duplicate_protein_id <int>

fwrite(
  compiled_proteins,
  file.path(out_dir, "Hpdi_compiled_signature_proteins_for_cis_pqtl.csv")
)
