# File: 4_1_coloc_all_cancer.R
# Analysis: All cancer.
# Active anchor input: hPDI reviewed instrument table; exact filename is in the input section.
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

# OpenGWAS credentials must be supplied privately via the OPENGWAS_JWT environment variable.
# ieugwasr::get_opengwas_jwt()

# install.packages("BiocManager")
# BiocManager::install("VariantAnnotation")
# remotes::install_github("MRCIEU/gwasvcf")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(readr)
  library(coloc)
  library(gwasvcf)
  library(VariantAnnotation)
  library(SummarizedExperiment)
  library(GenomeInfoDb)
  library(BiocGenerics)
})

set_bcftools("~/miniconda3/envs/bcftools_env/bin/bcftools")
#========================================================
# DIRECTORIES / DEFAULTS
#========================================================
pqtl_dir    <- file.path(DATA_DIR, "SNP_coloc/allcancer")
gwas_dir    <- file.path(DATA_DIR, "GWAS")
results_dir <- file.path(MANUSCRIPT_DIR, "results_coloc_local/allcancer")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

DEFAULT_CIS_WINDOW <- 500000

DEFAULT_MIN_SNPS   <- 50

#========================================================
# ANALYSIS TABLE: ONE ROW = ONE RUN
#========================================================
mr_hits <- tribble(
  ~protein,   ~disease,            ~gwas_id,
  "BPIFA2",     "All cancer",        "ieu-b-4966",
  "CFB",     "All cancer",        "ieu-b-4966",
  "CTSS",     "All cancer",        "ieu-b-4966",
  "VARS1",     "All cancer",        "ieu-b-4966",
  "KRT5",     "All cancer",        "ieu-b-4966"
)

#========================================================
# GWAS MANIFEST
#========================================================
gwas_manifest <- tribble(
  ~gwas_id,                   ~disease,           ~file,                                                   ~rsidx,        ~outcome_type, ~N,      ~ncase,  ~ncontrol,
  "ieu-b-4966",               "All cancer",       file.path(gwas_dir, "ieu-b-4966.vcf.gz"),                NA_character_, "cc",          442239,  70223,   372016
)

#========================================================
# BUILD GENE ANCHOR TABLE
#========================================================
exp_raw1 <- read.csv(file.path(MANUSCRIPT_DIR, "Hpdi_proteins_one_cis_IV_per_protein_MR_MANUAL.csv"))
exp_all <- exp_raw1 %>%
  dplyr::rename(
    chrom        = chrom_variant,
    pos_hg38     = pos_hg38_variant,
    p            = pval.exposure,
    rsid_first   = rsid,
    beta_aligned = beta.exposure
  )

exp_all <- exp_all %>%
  dplyr::mutate(
    chr_use = ifelse(grepl("^chr", chrom), chrom, paste0("chr", chrom)),
    pos_use = as.numeric(pos_hg38)
  )

print(names(exp_all))  # verify before proceeding

exp_all <- exp_all %>%
  dplyr::mutate(
    chr_use = ifelse(grepl("^chr", chrom), chrom, paste0("chr", chrom)),
    pos_use = as.numeric(pos_hg38)
  )

gene_loci <- exp_all %>%
  dplyr::filter(gene %in% unique(mr_hits$protein)) %>%
  dplyr::filter(!is.na(chr_use), !is.na(pos_use), !is.na(p)) %>%
  dplyr::group_by(gene) %>%
  dplyr::slice_min(order_by = p, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::transmute(
    protein = gene,
    chr     = ifelse(grepl("^chr", chr_use), chr_use, paste0("chr", chr_use)),
    start   = pos_use,
    end     = pos_use
  )


print(gene_loci)
print(setdiff(unique(mr_hits$protein), gene_loci$protein))

#========================================================
# INDEX pQTL FILES
#========================================================
pqtl_files <- list.files(
  path     = pqtl_dir,
  pattern  = "\\.txt\\.gz$",
  full.names = TRUE
)

extract_protein_from_filename <- function(x) {
  base  <- basename(x)
  base  <- str_remove(base, "\\.txt\\.gz$")
  parts <- str_split(base, "_", simplify = TRUE)
  if (ncol(parts) >= 3) parts[, 3] else NA_character_
}

file_index <- data.frame(
  file    = pqtl_files,
  protein = extract_protein_from_filename(pqtl_files),
  stringsAsFactors = FALSE
)

#========================================================
# READ ANNOTATION / EXCLUSION
#========================================================
annot    <- fread(file.path(pqtl_dir, "assocvariants.annotated.txt.gz"))
excluded <- fread(file.path(pqtl_dir, "assocvariants.excluded.txt.gz"))

standardize_annot <- function(dat) {
  if ("Chrom"           %in% names(dat)) setnames(dat, "Chrom",           "chr")
  if ("Pos"             %in% names(dat)) setnames(dat, "Pos",             "pos")
  if ("rsids"           %in% names(dat)) setnames(dat, "rsids",           "SNP")
  if ("effectAlleleFreq" %in% names(dat)) setnames(dat, "effectAlleleFreq", "maf")
  # Note: effectAllele and otherAllele must be present in annot for the join to work
  dat
}

standardize_excluded <- function(dat) {
  if ("rsids" %in% names(dat)) setnames(dat, "rsids", "SNP")
  dat
}

annot    <- standardize_annot(annot)
excluded <- standardize_excluded(excluded)

#========================================================
# pQTL HELPERS
#========================================================
standardize_pqtl <- function(dat) {
  if ("Chrom"  %in% names(dat)) setnames(dat, "Chrom",  "chr")
  if ("Pos"    %in% names(dat)) setnames(dat, "Pos",    "pos")
  if ("rsids"  %in% names(dat)) setnames(dat, "rsids",  "SNP")
  if ("Beta"   %in% names(dat)) setnames(dat, "Beta",   "beta")
  if ("Pval"   %in% names(dat)) setnames(dat, "Pval",   "p")
  if ("SE"     %in% names(dat)) setnames(dat, "SE",     "se")
  if ("ImpMAF" %in% names(dat)) setnames(dat, "ImpMAF", "maf")
  dat
}

prepare_pqtl <- function(pqtl_file, annot, excluded) {
  pqtl <- fread(pqtl_file)
  pqtl <- standardize_pqtl(pqtl)
  
  pqtl <- pqtl %>%
    dplyr::filter(
      !is.na(Name), !is.na(SNP), !is.na(beta), !is.na(se), !is.na(p),
      !is.na(chr),  !is.na(pos)
    ) %>%
    # FIX: deduplicate AFTER filtering but BEFORE the annotation join
    # so we don't drop variant entries prematurely
    dplyr::distinct(SNP, .keep_all = TRUE)
  
  if ("Name" %in% names(excluded)) {
    pqtl <- pqtl %>% dplyr::filter(!Name %in% excluded$Name)
  }
  
  # Identify which annotation columns actually exist
  annot_cols_available <- intersect(
    c("Name", "chr", "pos", "SNP", "effectAllele", "otherAllele", "maf"),
    names(annot)
  )
  
  pqtl <- pqtl %>%
    dplyr::left_join(
      annot %>% dplyr::select(dplyr::all_of(annot_cols_available)),
      by   = "Name",
      suffix = c("", ".annot")
    )
  
  # Coalesce only columns that exist after join
  if ("chr.annot"         %in% names(pqtl)) pqtl <- pqtl %>% mutate(chr = coalesce(chr.annot,         chr))
  if ("pos.annot"         %in% names(pqtl)) pqtl <- pqtl %>% mutate(pos = coalesce(pos.annot,         pos))
  if ("SNP.annot"         %in% names(pqtl)) pqtl <- pqtl %>% mutate(SNP = coalesce(SNP.annot,         SNP))
  if ("effectAllele.annot" %in% names(pqtl)) pqtl <- pqtl %>% mutate(effectAllele = coalesce(effectAllele.annot, effectAllele))
  if ("otherAllele.annot"  %in% names(pqtl)) pqtl <- pqtl %>% mutate(otherAllele  = coalesce(otherAllele.annot,  otherAllele))
  if ("maf.annot"         %in% names(pqtl)) pqtl <- pqtl %>% mutate(maf = coalesce(maf.annot,         maf))
  
  pqtl <- pqtl %>%
    dplyr::select(-dplyr::ends_with(".annot")) %>%
    dplyr::filter(!is.na(SNP), !is.na(maf)) %>%
    dplyr::distinct(SNP, .keep_all = TRUE)
  
  pqtl
}

clean_snp_ids <- function(x) {
  x <- as.character(x)
  x <- sub(",.*$", "", x)
  x <- trimws(x)
  x[x == ""] <- NA_character_
  x
}

#========================================================
# GWAS VCF HELPER
# FIX: function is now properly closed and returns a value
#========================================================
prepare_outcome_vcf_by_snp <- function(snps,
                                       gwas_vcf,
                                       rsidx        = NULL,
                                       fallback_ncase = NA,
                                       fallback_N   = NA,
                                       batch_size   = 200) {
  
  snps <- unique(clean_snp_ids(snps))
  snps <- snps[!is.na(snps)]
  snps <- snps[grepl("^rs", snps)]
  
  if (length(snps) == 0) {
    message("No usable rsIDs to query in GWAS VCF")
    return(NULL)
  }
  
  snp_chunks <- split(snps, ceiling(seq_along(snps) / batch_size))
  out_list   <- vector("list", length(snp_chunks))
  
  for (i in seq_along(snp_chunks)) {
    message("GWAS VCF batch ", i, " / ", length(snp_chunks),
            " (", length(snp_chunks[[i]]), " SNPs)")
    
    vcf <- tryCatch(
      {
        if (!is.null(rsidx) && !is.na(rsidx) && file.exists(rsidx)) {
          gwasvcf::query_gwas(gwas_vcf, rsid = snp_chunks[[i]], rsidx = rsidx)
        } else {
          gwasvcf::query_gwas(gwas_vcf, rsid = snp_chunks[[i]])
        }
      },
      error = function(e) {
        message("Batch ", i, " failed: ", e$message)
        return(NULL)
      }
    )
    
    if (is.null(vcf) || nrow(vcf) == 0) next
    
    rr <- SummarizedExperiment::rowRanges(vcf)
    
    out_list[[i]] <- tibble::tibble(
      chr   = as.character(GenomeInfoDb::seqnames(rr)),
      pos   = BiocGenerics::start(rr),
      SNP   = rownames(vcf),
      beta  = as.numeric(VariantAnnotation::geno(vcf)$ES[, 1]),
      se    = as.numeric(VariantAnnotation::geno(vcf)$SE[, 1]),
      lp    = as.numeric(VariantAnnotation::geno(vcf)$LP[, 1]),
      eaf   = as.numeric(VariantAnnotation::geno(vcf)$AF[, 1]),
      N     = dplyr::coalesce(as.numeric(VariantAnnotation::geno(vcf)$SS[, 1]), fallback_N),
      ncase = dplyr::coalesce(as.numeric(VariantAnnotation::geno(vcf)$NC[, 1]), fallback_ncase)
    ) %>%
      dplyr::mutate(
        chr      = ifelse(grepl("^chr", chr), chr, paste0("chr", chr)),
        p        = 10^(-lp),
        maf      = pmin(eaf, 1 - eaf),
        ncontrol = ifelse(!is.na(N) & !is.na(ncase), N - ncase, NA_real_),
        outcome_type = "cc"
      ) %>%
      dplyr::distinct(SNP, .keep_all = TRUE)
  }
  
  # FIX: assemble and return results
  result <- dplyr::bind_rows(Filter(Negate(is.null), out_list))
  if (nrow(result) == 0) return(NULL)
  result
}

#========================================================
# MAIN COLOC FUNCTION
# FIX: renamed argument to gwas_id_input to avoid name collision with column
#========================================================
run_coloc_one_local <- function(protein_name,
                                disease_name,
                                gwas_id_input,       # renamed from gwas_id
                                file_index,
                                annot,
                                excluded,
                                gwas_manifest,
                                gene_loci,
                                pqtl_N     = NULL,
                                min_snps   = 50,
                                cis_window = 500000) {
  
  matched_files <- file_index$file[file_index$protein == protein_name]
  
  if (length(matched_files) == 0) {
    message("No file found for protein: ", protein_name)
    return(NULL)
  }
  
  pqtl_file <- matched_files[1]
  message("Running: ", protein_name, " | ", disease_name, " | ", gwas_id_input)
  
  pqtl <- prepare_pqtl(
    pqtl_file = pqtl_file,
    annot     = annot,
    excluded  = excluded
  )
  
  pqtl <- pqtl %>%
    mutate(SNP = clean_snp_ids(SNP)) %>%
    filter(!is.na(SNP), !is.na(beta), !is.na(se), !is.na(maf),
           !is.na(chr),  !is.na(pos)) %>%
    distinct(SNP, .keep_all = TRUE)
  
  if (nrow(pqtl) < min_snps) {
    message("Too few pQTL SNPs for ", protein_name)
    return(NULL)
  }
  
  if (is.null(pqtl_N)) {
    if ("N" %in% names(pqtl)) {
      pqtl_N <- unique(na.omit(pqtl$N))[1]
    }
    if (is.null(pqtl_N) || is.na(pqtl_N)) {
      message("No pQTL sample size provided or found in file for ", protein_name)
      return(NULL)
    }
  }
  
  loc <- gene_loci %>% filter(protein == protein_name)
  
  if (nrow(loc) == 0) {
    message("No gene_loci entry found for protein: ", protein_name)
    return(NULL)
  }
  
  region_start  <- loc$start[1] - cis_window
  region_end    <- loc$end[1]   + cis_window
  region_chr    <- loc$chr[1]
  region_string <- paste0(region_chr, ":", region_start, "-", region_end)
  
  pqtl_region <- pqtl %>%
    filter(chr == region_chr,
           pos >= region_start,
           pos <= region_end)
  
  message("Region tested: ", region_string)
  message("Number of pQTL SNPs in region: ", nrow(pqtl_region))
  
  if (nrow(pqtl_region) < min_snps) {
    message("Too few SNPs in selected region for ", protein_name)
    return(NULL)
  }
  
  # FIX: use gwas_id_input (no name collision)
  gwas_meta <- gwas_manifest %>% filter(gwas_id == gwas_id_input)
  
  if (nrow(gwas_meta) == 0) {
    message("No GWAS file found in manifest for ", gwas_id_input)
    return(NULL)
  }
  
  gwas_file  <- gwas_meta$file[1]
  gwas_rsidx <- if ("rsidx" %in% names(gwas_meta) && !is.na(gwas_meta$rsidx[1])) {
    gwas_meta$rsidx[1]
  } else {
    NULL
  }
  
  outcome_region <- tryCatch(
    prepare_outcome_vcf_by_snp(
      snps          = pqtl_region$SNP,
      gwas_vcf      = gwas_file,
      rsidx         = gwas_rsidx,
      fallback_ncase = gwas_meta$ncase[1],
      fallback_N    = gwas_meta$N[1],
      batch_size    = 200
    ),
    error = function(e) {
      message("Local GWAS extraction failed for ", protein_name, " / ", disease_name, ": ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(outcome_region) || nrow(outcome_region) == 0) {
    message("No outcome data returned for ", protein_name, " / ", disease_name)
    return(NULL)
  }
  
  merged <- pqtl_region %>%
    inner_join(outcome_region, by = "SNP", suffix = c("", ".outcome")) %>%
    filter(
      !is.na(beta), !is.na(se), !is.na(maf),
      !is.na(beta.outcome), !is.na(se.outcome)
    ) %>%
    distinct(SNP, .keep_all = TRUE)
  
  message("Overlapping SNPs after merge: ", nrow(merged))
  
  if (nrow(merged) < min_snps) {
    message("Too few overlapping SNPs after merge for ", protein_name, " / ", disease_name)
    return(NULL)
  }
  
  dataset1 <- list(
    beta    = merged$beta,
    varbeta = merged$se^2,
    snp     = merged$SNP,
    MAF     = merged$maf,
    N       = pqtl_N,
    type    = "quant"
  )
  
  # FIX: warn explicitly if MAF falls back to pQTL values
  if ("maf.outcome" %in% names(merged)) {
    maf_outcome <- merged$maf.outcome
  } else if ("eaf.outcome" %in% names(merged)) {
    maf_outcome <- pmin(merged$eaf.outcome, 1 - merged$eaf.outcome)
  } else {
    message("Warning: maf.outcome and eaf.outcome not found — falling back to pQTL MAF for dataset2. ",
            "Results for ", protein_name, " / ", disease_name, " should be interpreted with caution.")
    maf_outcome <- merged$maf
  }
  
  outcome_type <- unique(na.omit(merged$outcome_type))
  outcome_type <- if (length(outcome_type) == 0) "cc" else outcome_type[1]
  
  if (outcome_type == "cc") {
    ncase    <- as.numeric(merged$ncase[1])
    ncontrol <- as.numeric(merged$ncontrol[1])
    outcome_N <- as.numeric(merged$N[1])
    
    if (!is.na(ncase) && !is.na(ncontrol)) {
      outcome_N <- ncase + ncontrol
      s_prop    <- ncase / outcome_N
    } else if (!is.na(ncase) && !is.na(outcome_N)) {
      s_prop <- ncase / outcome_N
    } else {
      message("Missing case/control counts for ", protein_name, " / ", disease_name)
      return(NULL)
    }
    
    dataset2 <- list(
      beta    = merged$beta.outcome,
      varbeta = merged$se.outcome^2,
      snp     = merged$SNP,
      MAF     = maf_outcome,
      N       = outcome_N,
      s       = s_prop,
      type    = "cc"
    )
    
  } else {
    outcome_N <- as.numeric(merged$N[1])
    
    if (is.na(outcome_N)) {
      message("Missing quantitative outcome sample size for ", protein_name, " / ", disease_name)
      return(NULL)
    }
    
    dataset2 <- list(
      beta    = merged$beta.outcome,
      varbeta = merged$se.outcome^2,
      snp     = merged$SNP,
      MAF     = maf_outcome,
      N       = outcome_N,
      type    = "quant"
    )
  }
  
  coloc_res <- tryCatch(
    coloc.abf(dataset1, dataset2),
    error = function(e) {
      message("coloc failed for ", protein_name, " / ", disease_name, ": ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(coloc_res)) return(NULL)
  
  summary_vec <- coloc_res$summary
  
  # Save full SNP-level results for plotting
  snp_results <- coloc_res$results
  snp_outfile <- file.path(
    results_dir,
    paste0(sprintf("%03d", hit_index), "_",
           protein_name, "_",
           gsub("[^A-Za-z0-9]+", "_", disease_name),
           "_snp_results.csv")
  )
  
fwrite(snp_results, snp_outfile)
  
  res <- data.frame(
    protein              = protein_name,
    disease              = disease_name,
    gwas_id              = gwas_id_input,
    pqtl_file            = basename(pqtl_file),
    outcome_file         = basename(gwas_file),
    region_type          = "cis_anchor",
    region               = region_string,
    cis_window           = cis_window,
    nsnps_region         = nrow(pqtl_region),
    nsnps_overlap        = nrow(merged),
    PP.H0                = unname(summary_vec["PP.H0.abf"]),
    PP.H1                = unname(summary_vec["PP.H1.abf"]),
    PP.H2                = unname(summary_vec["PP.H2.abf"]),
    PP.H3                = unname(summary_vec["PP.H3.abf"]),
    PP.H4                = unname(summary_vec["PP.H4.abf"]),
    coloc_interpretation = dplyr::case_when(
      unname(summary_vec["PP.H4.abf"]) >= 0.80 ~ "Strong colocalization",
      unname(summary_vec["PP.H4.abf"]) >= 0.50 ~ "Moderate colocalization",
      TRUE                                      ~ "Weak/no colocalization"
    ),
    stringsAsFactors = FALSE
  )
  
  print(coloc_res$summary)
  res
}

#========================================================
# CHOOSE ONE ROW TO RUN
#========================================================
args <- commandArgs(trailingOnly = TRUE)
 
 if (length(args) == 0) {
   stop("Provide one row number from mr_hits, e.g. Rscript coloc_hpdi.R 3")
 }
 
hit_index <- as.integer(args[1])

# TO THIS (to run row 1 interactively):
#hit_index <- 1  # change to 2, 3, 4 as needed

if (is.na(hit_index) || hit_index < 1 || hit_index > nrow(mr_hits)) {
  stop("Invalid row number. Must be between 1 and ", nrow(mr_hits))
}

this_hit <- mr_hits[hit_index, ]

message("Selected row: ", hit_index)
print(this_hit)

res <- run_coloc_one_local(
  protein_name  = this_hit$protein,
  disease_name  = this_hit$disease,
  gwas_id_input = this_hit$gwas_id,   # matches renamed argument
  file_index    = file_index,
  annot         = annot,
  excluded      = excluded,
  gwas_manifest = gwas_manifest,
  gene_loci     = gene_loci,
  pqtl_N        = NULL,
  min_snps      = DEFAULT_MIN_SNPS,
  cis_window    = DEFAULT_CIS_WINDOW
)

outfile <- file.path(
  results_dir,
  paste0(
    sprintf("%03d", hit_index), "_",
    this_hit$protein, "_",
    gsub("[^A-Za-z0-9]+", "_", this_hit$disease),
    "_coloc_local.csv"
  )
)

if (is.null(res)) {
  message("Run failed or returned NULL")
} else {
  fwrite(res, outfile)
  message("Saved result to: ", outfile)
  print(res)
}
