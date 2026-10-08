# Run from the repository root in a fresh R session. See README.md.
# EDIT LOCAL PATHS HERE. No study data or results are supplied.
PROJECT_ROOT <- "."
DATA_DIR <- file.path(PROJECT_ROOT, "data")
OUTPUT_ROOT <- file.path(PROJECT_ROOT, "outputs")
MANUSCRIPT_DIR <- file.path(OUTPUT_ROOT, "manuscript")
REFERENCE_DIR <- file.path(DATA_DIR, "reference")
SIGNATURE_SCRIPTS <- file.path(PROJECT_ROOT, "signatures")
CONTEXT_SCRIPTS <- file.path(PROJECT_ROOT, "contextualisation")
LIMMA_DATE <- "20260424"  # choose the manuscript-producing run/date
dir.create(MANUSCRIPT_DIR, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(dplyr))

#'---- Syntax for ProPlant Project'#
#' Created by Peggy Ler (Pler) - 17th November 2025
#' Filename..: build_weighted_scores.R
#' Updated...:  by "Pler"
#' --------------------------------------------------------------------
#' Using R studio 2024.04.02
#' -- Main program ----------------------------------------------------
#' 
#Load information of selected proteins from all models
#-----------------------------------------------------
# Paths
# -----
#############################################################
##create limma scores ---------------------------------------
#############################################################
# =================== USER CONFIG =====================================
OUTDIR        <- file.path(OUTPUT_ROOT, "updi")
EXTERNAL_DATA <- file.path(DATA_DIR, "analyticaldf_postexcl_20251020.RData")
EXTERNAL_OBJ  <- "df"

proj_root <- OUTPUT_ROOT
weights_file <- file.path(proj_root, paste0("hpdi/Hpidnovegoils_limma_", LIMMA_DATE, ".csv"))   
w <- read.csv(weights_file)

w1 <- filter(w, adj.P.Val<0.05)
nrow(w1) #619

#load full dataset
# Load data
load(file.path(DATA_DIR, "analyticaldf_postexcl_20251020.RData"))
cat(sprintf("Loaded data with %d rows\n", nrow(df)))

#build scores 
#----------------------------------
#keep proteins featured in ENET
common_enet <- intersect(w1$protein, colnames(df))
length(common_enet) #619

colnames(w1)
#[1] "protein"   "logFC"     "AveExpr"   "t"         "P.Value"   "adj.P.Val" "B"         "q_BH"      "q_BY"      "p_bonf" 


w1_use <- w1 %>%
  filter(protein %in% colnames(df)) %>%
  distinct(protein, .keep_all = TRUE)

prot_names <- w1_use$protein

# standardize selected proteins for score construction
X <- scale(df[, prot_names, drop = FALSE])

# weighted score
beta_vec <- w1_use$logFC
names(beta_vec) <- prot_names

df$hPDI_limma <- as.numeric(X %*% beta_vec)
df$hPDI_limma_z <- as.numeric(scale(df$hPDI_limma))

# unweighted aligned score
X_align <- as.data.frame(X)
neg_prots <- prot_names[w1_use$logFC < 0]
X_align[, neg_prots] <- X_align[, neg_prots] * -1

df$hPDI_limma_unw <- rowMeans(X_align)
df$hPDI_limma_unw_z <- as.numeric(scale(df$hPDI_limma_unw))


  
  score_vars <- c(
    "hPDI_limma",           # full prediction: intercept + covariates + proteins
    "hPDI_limma_z",       # protein contribution only
    "hPDI_limma_unw_z",    # protein contribution + intercept
    "hPDI_limma_unw"         # covariate contribution only
  )
  

  
  for (score_var in score_vars) {
    
    rho <- cor(df[[score_var]], df[["Hpdi_NoVegoils"]],
               method = "spearman", use = "complete.obs")
    
    r <- cor(df[[score_var]], df[["Hpdi_NoVegoils"]],
             method = "pearson", use = "complete.obs")
    
    cat(sprintf("\n%s:\n", score_var))
    cat(sprintf("  Spearman rho : %.4f\n", rho))
    cat(sprintf("  Pearson r    : %.4f\n", r))
  }

  # hPDI_limma:
  #   Spearman rho : 0.0695
  # Pearson r    : 0.0712
  # 
  # hPDI_limma_z:
  #   Spearman rho : 0.0695
  # Pearson r    : 0.0712
  # 
  # hPDI_limma_unw_z:
  #   Spearman rho : 0.0654
  # Pearson r    : 0.0666
  # 
  # hPDI_limma_unw:
  #   Spearman rho : 0.0654
  # Pearson r    : 0.0666
  
##----------------------
## NOW UPDI
## ---------------------

  weights_file <- file.path(proj_root, paste0("updi/Updinovegoils_limma_", LIMMA_DATE, ".csv"))   
  w <- read.csv(weights_file)
  
  w1 <- filter(w, adj.P.Val<0.05)
  nrow(w1) #51
  
  #load full dataset
  # Load data
  # Keep the same data frame so the hPDI scores remain available for export.
  cat(sprintf("Loaded data with %d rows\n", nrow(df)))
  
  #build scores 
  #----------------------------------
  #keep proteins featured in ENET
  common_enet <- intersect(w1$protein, colnames(df))
  length(common_enet) #51
  
  colnames(w1)
  #[1] "protein"   "logFC"     "AveExpr"   "t"         "P.Value"   "adj.P.Val" "B"         "q_BH"      "q_BY"      "p_bonf" 
  
  
  w1_use <- w1 %>%
    filter(protein %in% colnames(df)) %>%
    distinct(protein, .keep_all = TRUE)
  
  prot_names <- w1_use$protein
  
  # standardize selected proteins for score construction
  X <- scale(df[, prot_names, drop = FALSE])
  
  # weighted score
  beta_vec <- w1_use$logFC
  names(beta_vec) <- prot_names
  
  df$uPDI_limma <- as.numeric(X %*% beta_vec)
  df$uPDI_limma_z <- as.numeric(scale(df$uPDI_limma))
  
  # unweighted aligned score
  X_align <- as.data.frame(X)
  neg_prots <- prot_names[w1_use$logFC < 0]
  X_align[, neg_prots] <- X_align[, neg_prots] * -1
  
  df$uPDI_limma_unw <- rowMeans(X_align)
  df$uPDI_limma_unw_z <- as.numeric(scale(df$uPDI_limma_unw))
  
  
  
  score_vars <- c(
    "uPDI_limma",           # full prediction: intercept + covariates + proteins
    "uPDI_limma_z",       # protein contribution only
    "uPDI_limma_unw_z",    # protein contribution + intercept
    "uPDI_limma_unw"         # covariate contribution only
  )
  
  
  
  for (score_var in score_vars) {
    
    rho <- cor(df[[score_var]], df[["Updi_NoVegoils"]],
               method = "spearman", use = "complete.obs")
    
    r <- cor(df[[score_var]], df[["Updi_NoVegoils"]],
             method = "pearson", use = "complete.obs")
    
    cat(sprintf("\n%s:\n", score_var))
    cat(sprintf("  Spearman rho : %.4f\n", rho))
    cat(sprintf("  Pearson r    : %.4f\n", r))
  }
  
  # uPDI_limma:
  #   Spearman rho : 0.0978
  # Pearson r    : 0.0967
  # 
  # uPDI_limma_z:
  #   Spearman rho : 0.0978
  # Pearson r    : 0.0967
  # 
  # uPDI_limma_unw_z:
  #   Spearman rho : 0.0842
  # Pearson r    : 0.0841
  # 
  # uPDI_limma_unw:
  #   Spearman rho : 0.0842
  # Pearson r    : 0.0841

##save scores
dfgbm <- df[, c("Idepic", "hPDI_limma_unw", "hPDI_limma_unw_z", "hPDI_limma", "hPDI_limma_z",
                "uPDI_limma_unw", "uPDI_limma_unw_z", "uPDI_limma", "uPDI_limma_z")]

readr::write_csv(
  dfgbm,
  file.path(OUTPUT_ROOT, "limma_scores.csv")
)
