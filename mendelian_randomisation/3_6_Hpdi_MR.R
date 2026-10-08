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

# =============================================================================
# TWO-SAMPLE MENDELIAN RANDOMISATION — END TO END
# Exposure: cis-pQTLs (GRCh38)
# Outcomes: Local GWAS Catalog files (GRCh38) + IEU OpenGWAS (GRCh37)
# =============================================================================

library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(data.table)
library(openxlsx)

# =============================================================================
# 0) WORKBOOK SETUP
# =============================================================================

wb <- createWorkbook()

add_mr_sheet <- function(wb, sheet_name, res, dat = NULL) {
  res_out <- res %>%
    mutate(
      p_fdr    = p.adjust(pval, method = "BH"),
      OR       = exp(b),
      OR_lci95 = exp(b - 1.96 * se),
      OR_uci95 = exp(b + 1.96 * se)
    ) %>%
    arrange(pval)
  
  sheet_name <- substr(sheet_name, 1, 31)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, x = res_out)
  
  if (!is.null(dat)) {
    sheet2 <- substr(paste0(sheet_name, "_dat"), 1, 31)
    addWorksheet(wb, sheet2)
    writeData(wb, sheet = sheet2, x = dat)
  }
  
  invisible(res_out)
}

# =============================================================================
# 1) EXPOSURE DATA
# =============================================================================

exp_raw <- read.csv(file.path(MANUSCRIPT_DIR, "Hpdi_proteins_one_cis_IV_per_protein_MR_RAW_v2.csv"))

n_raw <- nrow(exp_raw)
cat("=== INSTRUMENT TRACKING ===\n")
cat("Step 0 — Raw instruments loaded:", n_raw, "\n\n")

# ── Step 1: Fix allele whitespace and cap pval=0 ────────────────────────────
exp_raw <- exp_raw %>%
  mutate(
    effect_allele.exposure = trimws(effect_allele.exposure),
    other_allele.exposure  = trimws(other_allele.exposure),
    pval_capped = ifelse(pval.exposure == 0, 1e-300, pval.exposure),
    Z     = qnorm(pval_capped / 2, lower.tail = FALSE),
    se    = abs(beta.exposure / Z),
    Fstat = (beta.exposure / se)^2
  )

# smmary(exp_raw$Fstat)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 36.29   81.11  202.04  434.96  625.53 1379.76 

# ── Step 3: Track each filtering step ───────────────────────────────────────

# 3a. RSID filter — drops chr:pos without an RSID
after_rsid <- exp_raw %>% # Keep if it has a valid rsid OR if it's AMY1A chr:pos
  filter(grepl("^rs\\d+$", trimws(tolower(rsid))) | rsid == "chr1:103682950")
dropped_no_rsid <- exp_raw %>%
  filter(!(grepl("^rs\\d+$", trimws(tolower(rsid))) | rsid == "chr1:103682950")) %>% # ✅
  select(gene, rsid, variant) %>%
  mutate(drop_reason = "No RSID — chr:pos variant absent from dbSNP")

cat("Step 1 — Dropped (no RSID / unresolvable chr:pos):", nrow(dropped_no_rsid),
    sprintf("(%.1f%%)\n", 100 * nrow(dropped_no_rsid) / n_raw))
print(dropped_no_rsid)
cat("\n")

# 3b. Allele filter — drops complex indels with non-ACGT characters
after_allele <- after_rsid %>%
  filter(
    grepl("^[ACGT]+$", effect_allele.exposure),
    grepl("^[ACGT]+$", other_allele.exposure)
  )
dropped_allele <- after_rsid %>%
  filter(!(
    grepl("^[ACGT]+$", effect_allele.exposure) &
      grepl("^[ACGT]+$", other_allele.exposure)
  )) %>%
  select(gene, rsid, effect_allele.exposure, other_allele.exposure) %>%
  mutate(drop_reason = "Complex indel — non-ACGT allele characters")

cat("Step 2 — Dropped (complex indels):", nrow(dropped_allele),
    sprintf("(%.1f%%)\n", 100 * nrow(dropped_allele) / n_raw))
print(dropped_allele)
cat("\n")

# 3c. SE filter — drops instruments where se=0 remains after pval capping
after_se <- after_allele %>% filter(is.finite(se), se > 0)
dropped_se <- after_allele %>%
  filter(!(is.finite(se) & se > 0)) %>%
  select(gene, rsid, beta.exposure, pval.exposure, se) %>%
  mutate(drop_reason = "SE = 0 or non-finite after pval capping")

cat("Step 3 — Dropped (invalid SE):", nrow(dropped_se),
    sprintf("(%.1f%%)\n", 100 * nrow(dropped_se) / n_raw))
print(dropped_se)
cat("\n")

# 3d. Deduplication — keep best instrument per gene
exp_raw2 <- after_se %>%
  mutate(rsid = trimws(tolower(rsid))) %>%
  arrange(pval.exposure) %>%
  distinct(gene, rsid, .keep_all = TRUE)

dropped_dedup <- nrow(after_se) - nrow(exp_raw2)
cat("Step 4 — Dropped (duplicates):", dropped_dedup,
    sprintf("(%.1f%%)\n", 100 * dropped_dedup / n_raw))
#Step 4 — Dropped (duplicates): 44 (9.4%)
cat("Final instruments for MR:", nrow(exp_raw2),
    "across", length(unique(exp_raw2$gene)), "proteins\n\n")
#Final instruments for MR: 414 across 408 proteins

# ── Step 4: F-statistic summary ─────────────────────────────────────────────
cat("F-statistic summary (instruments entering MR):\n")
print(summary(exp_raw2$Fstat))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 36.29  106.45  326.17  551.77 1051.98 1379.76 
cat("Instruments with F < 10:", sum(exp_raw2$Fstat < 10, na.rm = TRUE), "\n\n")
#Instruments with F < 10: 0 
# ── Step 5: Format exposure ─────────────────────────────────────────────────
exp_dat <- format_data(
  exp_raw2,
  type                 = "exposure",
  snp_col              = "rsid",
  beta_col             = "beta.exposure",
  se_col               = "se",
  effect_allele_col    = "effect_allele.exposure",
  other_allele_col     = "other_allele.exposure",
  eaf_col              = "eaf.exposure",
  pval_col             = "pval.exposure",
  phenotype_col        = "gene"
)

# =============================================================================
# 2) OUTCOME PROCESSING FUNCTION
# =============================================================================
# Handles both:
#   - Local GWAS Catalog files (GRCh38, matched by rsid)
#   - IEU OpenGWAS API (GRCh37, queried by rsid)
#   - Outcomes reporting OR + CI instead of beta + SE (e.g. liver)
#
# Parameters:
#   is_or        — set TRUE if the outcome file reports odds ratios
#   beta_col     — when is_or=TRUE, pass the OR column name here
#   ci_upper_col — column name for upper 95% CI bound (used when is_or=TRUE)
#   ci_lower_col — column name for lower 95% CI bound (used when is_or=TRUE)
#   se_col       — only used when is_or=FALSE (direct SE column)

run_mr_for_outcome <- function(
    wb,
    exp_dat,
    outcome_label,
    outcome_name,
    # For local files:
    outcome_file  = NULL,
    snp_col       = "rsid",
    beta_col      = "beta",
    se_col        = "standard_error",
    ea_col        = "effect_allele",
    oa_col        = "other_allele",
    eaf_col       = "effect_allele_frequency",
    pval_col      = "p_value",
    # OR-specific parameters
    is_or         = FALSE,
    ci_upper_col  = "ci_upper",
    ci_lower_col  = "ci_lower",
    # For IEU OpenGWAS:
    ieu_id        = NULL
) {
  
  cat("\n", strrep("=", 60), "\n")
  cat("OUTCOME:", outcome_name, "\n")
  if (is_or) cat("NOTE: Outcome reports OR — converting to log scale\n")
  cat(strrep("=", 60), "\n")
  
  n_exp_snps <- length(unique(exp_dat$SNP))
  
  # ── Load outcome ────────────────────────────────────────────────────────────
  if (!is.null(ieu_id)) {
    
    cat("Source: IEU OpenGWAS —", ieu_id, "\n")
    out_dat <- extract_outcome_data(
      snps     = exp_dat$SNP,
      outcomes = ieu_id
    )
    
    if (is.null(out_dat) || nrow(out_dat) == 0) {
      cat("WARNING: No outcome data returned from IEU OpenGWAS. Skipping.\n")
      return(NULL)
    }
    
    n_found   <- length(unique(out_dat$SNP))
    n_missing <- n_exp_snps - n_found
    cat(sprintf("SNPs requested: %d | Found: %d (%.1f%%) | Missing: %d (%.1f%%)\n",
                n_exp_snps, n_found, 100 * n_found / n_exp_snps,
                n_missing, 100 * n_missing / n_exp_snps))
    
  } else {
    
    cat("Source: Local file —", basename(outcome_file), "\n")
    out_raw <- data.table::fread(outcome_file)
    
    out_std <- out_raw %>%
      mutate(
        SNP = trimws(tolower(.data[[snp_col]])),
        
        # Beta: log-transform OR if needed, otherwise use directly
        beta_outcome = if (is_or) {
          log(as.numeric(.data[[beta_col]]))
        } else {
          as.numeric(.data[[beta_col]])
        },
        
        # SE: derive from CI bounds if OR, otherwise use SE column directly
        se_outcome = if (is_or) {
          (log(as.numeric(.data[[ci_upper_col]])) -
             log(as.numeric(.data[[ci_lower_col]]))) / (2 * 1.96)
        } else {
          as.numeric(.data[[se_col]])
        },
        
        pval_outcome          = as.numeric(.data[[pval_col]]),
        effect_allele.outcome = toupper(.data[[ea_col]]),
        other_allele.outcome  = toupper(.data[[oa_col]]),
        eaf.outcome           = as.numeric(.data[[eaf_col]]),
        outcome_name          = outcome_name
      ) %>%
      filter(grepl("^rs\\d+$", SNP)) %>%
      filter(SNP %in% unique(exp_dat$SNP)) %>%
      filter(is.finite(beta_outcome), is.finite(se_outcome), se_outcome > 0) %>%
      as.data.frame()
    
    n_found   <- length(unique(out_std$SNP))
    n_missing <- n_exp_snps - n_found
    n_absent  <- sum(!unique(exp_dat$SNP) %in% trimws(tolower(out_raw[[snp_col]])))
    
    cat(sprintf("SNPs requested: %d | Found: %d (%.1f%%) | Missing: %d (%.1f%%)\n",
                n_exp_snps, n_found, 100 * n_found / n_exp_snps,
                n_missing, 100 * n_missing / n_exp_snps))
    cat(sprintf("  Of missing: %d absent from outcome GWAS entirely\n", n_absent))
    
    if (nrow(out_std) == 0) {
      cat("WARNING: No SNPs matched. Skipping.\n")
      return(NULL)
    }
    
    out_dat <- format_data(
      out_std,
      type              = "outcome",
      snp_col           = "SNP",
      beta_col          = "beta_outcome",
      se_col            = "se_outcome",
      effect_allele_col = "effect_allele.outcome",
      other_allele_col  = "other_allele.outcome",
      eaf_col           = "eaf.outcome",
      pval_col          = "pval_outcome",
      phenotype_col     = "outcome_name"
    )
  }
  
  # ── Harmonise ───────────────────────────────────────────────────────────────
  dat <- harmonise_data(exp_dat, out_dat)
  
  n_keep     <- sum(dat$mr_keep)
  n_drop_har <- sum(!dat$mr_keep)
  
  cat(sprintf("\nHarmonisation: kept %d | dropped %d\n", n_keep, n_drop_har))
  
  drop_summary <- dat %>%
    filter(!mr_keep) %>%
    mutate(
      drop_reason = case_when(
        ambiguous == TRUE                      ~ "Palindromic — ambiguous strand (intermediate MAF)",
        grepl("incompatible", tolower(remove)) ~ "Incompatible alleles",
        TRUE                                   ~ "Other"
      )
    ) %>%
    count(drop_reason) %>%
    mutate(pct = sprintf("%.1f%%", 100 * n / (n_keep + n_drop_har)))
  
  cat("Drop reasons after harmonisation:\n")
  print(as.data.frame(drop_summary))
  
  # Proteins lost
  proteins_kept <- dat %>% filter(mr_keep) %>% distinct(exposure) %>% pull(exposure)
  proteins_lost <- unique(exp_dat$exposure)[!unique(exp_dat$exposure) %in% proteins_kept]
  
  cat(sprintf("\nProteins with valid instrument: %d / %d (%.1f%%)\n",
              length(proteins_kept),
              length(unique(exp_dat$exposure)),
              100 * length(proteins_kept) / length(unique(exp_dat$exposure))))
  cat(sprintf("Proteins with NO valid instrument: %d (%.1f%%)\n",
              length(proteins_lost),
              100 * length(proteins_lost) / length(unique(exp_dat$exposure))))
  
  if (length(proteins_lost) > 0) {
    lost_detail <- exp_dat %>%
      filter(exposure %in% proteins_lost) %>%
      left_join(
        dat %>% select(SNP, mr_keep, ambiguous, remove),
        by = "SNP"
      ) %>%
      mutate(
        loss_reason = case_when(
          is.na(mr_keep)    ~ "SNP absent from outcome GWAS",
          ambiguous == TRUE ~ "Palindromic — dropped at harmonisation",
          !is.na(remove)    ~ paste("Harmonisation:", remove),
          TRUE              ~ "Other"
        )
      ) %>%
      select(exposure, SNP, eaf.exposure, loss_reason) %>%
      arrange(loss_reason, exposure)
    
    cat("\nLost proteins detail:\n")
    print(as.data.frame(lost_detail))
  } else {
    lost_detail <- NULL
  }
  
  # ── Run MR ─────────────────────────────────────────────────────────────────
  if (n_keep == 0) {
    cat("WARNING: No instruments remain after harmonisation. Skipping MR.\n")
    return(NULL)
  }
  
  res <- mr(dat) %>%
    mutate(
      p_fdr    = p.adjust(pval, method = "BH"),
      OR       = exp(b),
      OR_lci95 = exp(b - 1.96 * se),
      OR_uci95 = exp(b + 1.96 * se)
    ) %>%
    arrange(pval)
  
  cat(sprintf("\nMR complete: %d protein-outcome pairs tested\n", nrow(res)))
  cat(sprintf("Significant after FDR (p_fdr < 0.05): %d\n",
              sum(res$p_fdr < 0.05, na.rm = TRUE)))
  
  # ── Save to workbook ────────────────────────────────────────────────────────
  add_mr_sheet(wb, outcome_label, res, dat)
  
  # Save lost proteins log to workbook
  if (length(proteins_lost) > 0) {
    lost_sheet <- substr(paste0(outcome_label, "_lost"), 1, 31)
    addWorksheet(wb, lost_sheet)
    writeData(wb, sheet = lost_sheet, x = lost_detail)
  }
  
  return(list(res = res, dat = dat, proteins_lost = proteins_lost))
}

# =============================================================================
# 3) RUN MR FOR ALL OUTCOMES
# =============================================================================

# ── Local GWAS Catalog outcomes (GRCh38) ────────────────────────────────────

results_colorectal <- run_mr_for_outcome(
  wb            = wb,
  exp_dat       = exp_dat,
  outcome_label = "Colorectal_GCST90255675",
  outcome_name  = "Colorectal cancer (GCST90255675)",
  outcome_file  = file.path(DATA_DIR, "GWAS/colorectal/GCST90255675.h.tsv.gz")
)

# SNPs requested: 408 | Found: 320 (78.4%) | Missing: 88 (21.6%)
#Of missing: 88 absent from outcome GWAS entirely
# Harmonisation: kept 274 | dropped 49
# Drop reasons after harmonisation:
#   drop_reason  n   pct
# 1                                             Other  3  0.9%
# 2 Palindromic — ambiguous strand (intermediate MAF) 46 14.2%
# 
# Proteins with valid instrument: 272 / 408 (66.7%)
# Proteins with NO valid instrument: 136 (33.3%)
# MR complete: 272 protein-outcome pairs tested
# Significant after FDR (p_fdr < 0.05): 8

results_kidney <- run_mr_for_outcome(
  wb            = wb,
  exp_dat       = exp_dat,
  outcome_label = "Kidney_GCST90320057",
  outcome_name  = "Kidney cancer (GCST90320057)",
  outcome_file  = file.path(DATA_DIR, "GWAS/Kidney_cancer/GCST90320057.tsv.gz"),
)
# MR complete: 334 protein-outcome pairs tested
# Significant after FDR (p_fdr < 0.05): 0

results_prostate <- run_mr_for_outcome(
  wb            = wb,
  exp_dat       = exp_dat,
  outcome_label = "Prostate_GCST90274714",
  outcome_name  = "Prostate cancer (GCST90274714)",
  outcome_file  = file.path(DATA_DIR, "GWAS/prostate/GCST90274714.h.tsv.gz")
)

# MR complete: 342 protein-outcome pairs tested
# Significant after FDR (p_fdr < 0.05): 23

results_liver <- run_mr_for_outcome(
  wb            = wb,
  exp_dat       = exp_dat,
  outcome_label = "Liver_GCST90475568",
  outcome_name  = "Liver cancer (GCST90475568)",
  outcome_file  = file.path(DATA_DIR, "GWAS/liver/GCST90475568.h.tsv.gz"),
  beta_col      = "odds_ratio",  # log-transformed internally
  is_or         = TRUE,
  ci_upper_col  = "ci_upper",
  ci_lower_col  = "ci_lower"
)

# MR complete: 334 protein-outcome pairs tested
# Significant after FDR (p_fdr < 0.05): 1
results_lung <- run_mr_for_outcome(
  wb            = wb,
  exp_dat       = exp_dat,
  outcome_label = "Lung_GCST004748",
  outcome_name  = "Lung cancer (GCST004748)",
  outcome_file  = file.path(DATA_DIR, "GWAS/lung/28604730-GCST004748-EFO_0001071.h.tsv.gz"),
  snp_col       = "hm_rsid",
  ea_col        = "hm_effect_allele",
  oa_col        = "hm_other_allele",
  eaf_col       = "hm_effect_allele_frequency",
  beta = "hm_beta"
)
# MR complete: 272 protein-outcome pairs tested
# Significant after FDR (p_fdr < 0.05): 1

results_myeloma <- run_mr_for_outcome(
  wb            = wb,
  exp_dat       = exp_dat,
  outcome_label = "Myeloma_GCST90624747",
  outcome_name  = "Myeloma (GCST90624747)",
  outcome_file  = file.path(DATA_DIR, "GWAS/myeloma/GCST90624747.h.tsv.gz")
)
# 
# MR complete: 302 protein-outcome pairs tested
# Significant after FDR (p_fdr < 0.05): 0

results_gastric <- run_mr_for_outcome(
  wb            = wb,
  exp_dat       = exp_dat,
  outcome_label = "Gastric_GCST90477175",
  outcome_name  = "Gastric cancer (GCST90477175)",
  outcome_file  = file.path(DATA_DIR, "GWAS/gastric/GCST90477175.h.tsv.gz"),
  beta_col      = "odds_ratio",  # log-transformed internally
  is_or         = TRUE,
  ci_upper_col  = "ci_upper",
  ci_lower_col  = "ci_lower"
)
# MR complete: 334 protein-outcome pairs tested
# Significant after FDR (p_fdr < 0.05): 0
results_endom <- run_mr_for_outcome(
  wb            = wb,
  exp_dat       = exp_dat,
  outcome_label = "Endom_GCST90477175",
  outcome_name  = "Endometrial cancer (GCST90477175)",
  outcome_file  = file.path(DATA_DIR, "GWAS/endom/30093612-GCST006465-EFO_1001512.h.tsv.gz"),
  snp_col       = "hm_rsid",
  ea_col        = "hm_effect_allele",
  oa_col        = "hm_other_allele",
  eaf_col       = "hm_effect_allele_frequency",
  beta = "hm_beta"
)
# MR complete: 297 protein-outcome pairs tested
# Significant after FDR (p_fdr < 0.05): 2

# ── IEU OpenGWAS outcomes (GRCh37 — queried by RSID, build does not matter) ─

results_breast <- run_mr_for_outcome(
  wb            = wb,
  exp_dat       = exp_dat,
  outcome_label = "Breast_ieu-a-1126",
  outcome_name  = "Breast cancer (ieu-a-1126)",
  ieu_id        = "ieu-a-1126"
)
# MR complete: 301 protein-outcome pairs tested
# Significant after FDR (p_fdr < 0.05): 5

results_all_cancer <- run_mr_for_outcome(
  wb            = wb,
  exp_dat       = exp_dat,
  outcome_label = "ieu-b-4966",
  outcome_name  = "Outcome ieu-b-4966",
  ieu_id        = "ieu-b-4966"
)

# MR complete: 321 protein-outcome pairs tested
# Significant after FDR (p_fdr < 0.05): 4

results_t2d <- run_mr_for_outcome(
  wb            = wb,
  exp_dat       = exp_dat,
  outcome_label = "T2D_GCST006867",
  outcome_name  = "Type 2 diabetes (ebi-a-GCST006867)",
  ieu_id        = "ebi-a-GCST006867"
)

# MR complete: 206 protein-outcome pairs tested
# Significant after FDR (p_fdr < 0.05): 6

results_cvd <- run_mr_for_outcome(
  wb            = wb,
  exp_dat       = exp_dat,
  outcome_label = "CVD_GCST90029019",
  outcome_name  = "Cardiovascular disease (ebi-a-GCST90029019)",
  ieu_id        = "ebi-a-GCST90029019"
)
# MR complete: 330 protein-outcome pairs tested
# Significant after FDR (p_fdr < 0.05): 17

# =============================================================================
# 4) GLOBAL INSTRUMENT LOSS SUMMARY ACROSS ALL OUTCOMES
# =============================================================================

all_results <- list(
  Colorectal = results_colorectal,
  Kidney     = results_kidney,
  Prostate   = results_prostate,
  Liver      = results_liver,
  Lung       = results_lung,
  Gastric    = results_gastric,
  Myeloma    = results_myeloma,
  Endometrial = results_endom,
  Breast     = results_breast,
  Cancer = results_all_cancer,
  T2D        = results_t2d,
  CVD        = results_cvd
)

cat("\n", strrep("=", 60), "\n")
cat("GLOBAL SUMMARY — INSTRUMENTS ACROSS ALL OUTCOMES\n")
cat(strrep("=", 60), "\n")

summary_table <- lapply(names(all_results), function(nm) {
  res <- all_results[[nm]]
  if (is.null(res)) return(NULL)
  
  dat <- res$dat
  n_total    <- length(unique(exp_dat$exposure))
  n_kept     <- dat %>% filter(mr_keep) %>% distinct(exposure) %>% nrow()
  n_lost     <- length(res$proteins_lost)
  
  data.frame(
    Outcome              = nm,
    Total_proteins       = n_total,
    Proteins_with_IV     = n_kept,
    Proteins_lost        = n_lost,
    Pct_lost             = sprintf("%.1f%%", 100 * n_lost / n_total),
    Instruments_kept     = sum(dat$mr_keep),
    Palindromic_dropped  = sum(dat$ambiguous == TRUE, na.rm = TRUE),
    stringsAsFactors     = FALSE
  )
}) %>%
  bind_rows()

print(summary_table)

addWorksheet(wb, "Summary")
writeData(wb, sheet = "Summary", x = summary_table)

# =============================================================================
# 5) SAVE WORKBOOK
# =============================================================================

out_path <- file.path(MANUSCRIPT_DIR, "hpdi_MR.xlsx")

saveWorkbook(wb, out_path, overwrite = TRUE)
cat("\nResults saved to:", out_path, "\n")
