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

#'---- Syntax for ProPlant Project'#
#' Created by Peggy Ler (Pler) - 31 OCt 2025
#' Filename..: Overlap_GB_ENET.R
#' Updated...: "" by "Pler"
#' Purpose...: 
#' 
#' 
#' Using R studio 2024.04.02
#' -- Main program ----------------------------------------------

library(ggplot2)
library(dplyr)
library(forcats)
library(ggrepel)
library(UpSetR)
library(scales)

# Paths
# -----
proj_root <- MANUSCRIPT_DIR

##HPDI
gb_file <- file.path(proj_root, "Hpdi_XGBoost_proteinlist.csv")   
limma_file <- file.path(proj_root, "Hpdi_LIMMA_proteinlist.csv")  
enet_file <- file.path(proj_root, "Hpdi_ENET_proteinlist.csv")  

hPDI_XGB  <- read.csv(gb_file, stringsAsFactors = FALSE)
hPDI_limma  <- read.csv(limma_file, stringsAsFactors = FALSE)
hPDI_ENET  <- read.csv(enet_file, stringsAsFactors = FALSE)


##
    
S1 <- hPDI_limma$gene
S2 <- hPDI_ENET$gene
S3 <- hPDI_XGB$gene 

sets <- list(
  hPDI_limma    = S1,
  hPDI_ENET  = S2,
  hPDI_XGB   = S3
)

# size of each signature
sapply(sets, length)
# hPDI_limma  hPDI_ENET   hPDI_XGB 
# 619        108        667 

#jaccard
# ------
jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

models <- names(sets)
J <- matrix(NA_real_,
            nrow = length(sets),
            ncol = length(sets),
            dimnames = list(models, models))

for (i in seq_along(sets)) {
  for (j in seq_along(sets)) {
    J[i, j] <- jaccard(sets[[i]], sets[[j]])
  }
}

J    
# XGB      limma       ENET
# XGB   1.00000000 0.08792185 0.06506365
# limma 0.08792185 1.00000000 0.16020236
# ENET  0.06506365 0.16020236 1.00000000

# pairwise intersections (counts)
pairwise_intersect_counts <- outer(
  sets, sets,
  Vectorize(function(a, b) length(intersect(a, b)))
)
dimnames(pairwise_intersect_counts) <- list(names(sets), names(sets))
pairwise_intersect_counts
# XGB limma ENET
# XGB   645    99   46
# limma  99   580   95
# ENET   46    95  108

# Visualise
# ---------
J_df <- as.data.frame(as.table(J))  

p<- ggplot(J_df, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = percent(Freq, accuracy = 1)), 
            color = "white", size = 8) +  # 
  scale_fill_gradient(low = "#D27E64", high = "#2B3256") +
  labs(fill = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank()
  )
p

# ggsave(
#   filename = file.path(out_dir, "overlap_percentages_Enet_GBM_full.png"),
#   plot     = p,
#   width    = 9,
#   height   = 6,
#   dpi      = 300
# )

#Upset plots
#-----------
sets <- sets[c("hPDI_limma", "hPDI_ENET", "hPDI_XGB")]
names(sets) <- c("hPDI-limma", "hPDI-ENET", "hPDI-XGB")

upset_obj <- fromList(sets)

upset(
  upset_obj,
  sets = c("hPDI-limma", "hPDI-ENET", "hPDI-XGB"),
  keep.order = TRUE,
  nsets = 3,
  nintersects = 15,
  sets.bar.color = "#D27E64",
  main.bar.color = "#2B3256",
  matrix.color = "grey30"
)

#save
png(
  filename = file.path(proj_root, "intersect_upsets_hPDI_limma_ENET_XGB.png"),
  width = 9,
  height = 6,
  units = "in",
  res = 300
)

upset(
  upset_obj,
  sets = c("hPDI-limma", "hPDI-ENET", "hPDI-XGB"),
  keep.order = TRUE,
  nsets = 3,
  nintersects = 15,
  sets.bar.color = "#D27E64",
  main.bar.color = "#2B3256",
  matrix.color = "grey30"
)

dev.off()

names(sets)

# proteins present in all three sets
overlap_all3 <- Reduce(
  intersect,
  sets[c("hPDI-limma", "hPDI-ENET", "hPDI-XGB")]
)

length(overlap_all3)
overlap_all3_df <- tibble(
  gene= overlap_all3
)

overlap_all3_annotated <- overlap_all3_df %>%
  left_join(
    hPDI_XGB %>% dplyr::select(gene, protein_id, target, uniprot),
    by = "gene"
  ) %>%
  dplyr::select(gene, protein_id, target, uniprot)

overlap_all3_annotated

write.csv(
  overlap_all3_annotated,
  file.path(proj_root, "hPDI_overlap_all3_target_gene.csv"),
  row.names = FALSE
)

## UPDI
## ____
##HPDI
gb_file <- file.path(proj_root, "Updi_XGBoost_proteinlist.csv")   
limma_file <- file.path(proj_root, "Updi_LIMMA_proteinlist.csv")  
enet_file <- file.path(proj_root, "Updi_ENET_proteinlist.csv")  

uPDI_XGB  <- read.csv(gb_file, stringsAsFactors = FALSE)
uPDI_limma  <- read.csv(limma_file, stringsAsFactors = FALSE)
uPDI_ENET  <- read.csv(enet_file, stringsAsFactors = FALSE)


##

S1 <- uPDI_limma$gene
S2 <- uPDI_ENET$gene
S3 <- uPDI_XGB$gene 

sets <- list(
  uPDI_limma    = S1,
  uPDI_ENET  = S2,
  uPDI_XGB   = S3
)

# size of each signature
sapply(sets, length)
# uPDI_limma  uPDI_ENET   uPDI_XGB 
# 51         14        954 

sets <- sets[c("uPDI_limma", "uPDI_ENET", "uPDI_XGB")]
names(sets) <- c("uPDI-limma", "uPDI-ENET", "uPDI-XGB")

upset_obj <- fromList(sets)

upset(
  upset_obj,
  sets = c("uPDI-limma", "uPDI-ENET", "uPDI-XGB"),
  keep.order = TRUE,
  nsets = 3,
  nintersects = 15,
  sets.bar.color = "#D27E64",
  main.bar.color = "#2B3256",
  matrix.color = "grey30"
)

#save
png(
  filename = file.path(proj_root, "intersect_upsets_uPDI_limma_ENET_XGB.png"),
  width = 9,
  height = 6,
  units = "in",
  res = 300
)

upset(
  upset_obj,
  sets = c("uPDI-limma", "uPDI-ENET", "uPDI-XGB"),
  keep.order = TRUE,
  nsets = 3,
  nintersects = 15,
  sets.bar.color = "#D27E64",
  main.bar.color = "#2B3256",
  matrix.color = "grey30"
)

dev.off()

names(sets)

# proteins present in all three sets
overlap_all3 <- Reduce(
  intersect,
  sets[c("uPDI-limma", "uPDI-ENET", "uPDI-XGB")]
)

overlap_all3
length(overlap_all3)
overlap_all3_df <- tibble(
  gene= overlap_all3
)

overlap_all3_annotated <- overlap_all3_df %>%
  left_join(
    uPDI_XGB %>% dplyr::select(gene, protein_id, target, uniprot),
    by = "gene"
  ) %>%
  dplyr::select(gene, protein_id, target, uniprot)

overlap_all3_annotated

write.csv(
  overlap_all3_annotated,
  file.path(proj_root, "uPDI_overlap_all3_target_gene.csv"),
  row.names = FALSE
)

#---------------------------------------------
# comparisons of importance and coefficients
#------------------------------------------
# set direction of effects for elastic net
Enet$direction <- sign(Enet$coef)

## Elastic Net bootstrap 60% model
Enet_60$direction <- sign(Enet_60$coef)

Enet$direction    <- ifelse(Enet$coef    >= 0, "+", "-")
Enet_60$direction <- ifelse(Enet_60$coef >= 0, "+", "-")

# create vector list
eff_gbm_best <- tapply(GBM_best$weight_signed, GBM_best$Gene, mean, na.rm = TRUE)
eff_GBMfull  <- tapply(GBMfull$weight_signed,  GBMfull$Gene,  mean, na.rm = TRUE)
eff_enet     <- tapply(Enet$coef,             Enet$Gene,      mean, na.rm = TRUE)
eff_enet_60  <- tapply(Enet_60$coef,          Enet_60$Gene,   mean, na.rm = TRUE)

effects <- list(
  GBM_best  = eff_gbm_best,
  GBMfull   = eff_GBMfull,
  Enet      = eff_enet,
  Enet_60   = eff_enet_60
)

# pair wise comparisons with sign agreements
# -------------------------------------------
sign_agreement_matrix <- function(effects, drop_zeros = TRUE) {
  models <- names(effects)
  n <- length(effects)
  M <- matrix(NA_real_, nrow = n, ncol = n,
              dimnames = list(models, models))
  
  for (i in seq_along(effects)) {
    for (j in seq_along(effects)) {
      a <- effects[[i]]
      b <- effects[[j]]
      common <- intersect(names(a), names(b))
      if (length(common) == 0) {
        M[i, j] <- NA
      } else {
        sa <- sign(a[common])
        sb <- sign(b[common])
        if (drop_zeros) {
          keep <- (sa != 0 | sb != 0)
          sa <- sa[keep]
          sb <- sb[keep]
        }
        if (length(sa) == 0) {
          M[i, j] <- NA
        } else {
          M[i, j] <- mean(sa == sb)
        }
      }
    }
  }
  M
}

sign_agree <- sign_agreement_matrix(effects)
sign_agree
# GBM_best   GBMfull      Enet   Enet_60
# GBM_best 1.0000000 0.4931973 0.6857143 0.7000000
# GBMfull  0.4931973 1.0000000 0.4166667 0.3571429
# Enet     0.6857143 0.4166667 1.0000000 1.0000000
# Enet_60  0.7000000 0.3571429 1.0000000 1.0000000

df <- as.data.frame(as.table(sign_agree))  

p<- ggplot(df, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = percent(Freq, accuracy = 1)), 
            color = "white", size = 8) +  # 
  scale_fill_gradient(low = "#D27E64", high = "#2B3256") +
  labs(fill = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank()
  )
p

ggsave(
  filename = file.path(out_dir, "sign_agreement_percentages_4models.png"),
  plot     = p,
  width    = 9,
  height   = 6,
  dpi      = 300
)

#pairwise corrlation of effects sizes
#-------------------------------------
effect_correlation_matrix <- function(effects, method = "spearman", abs = FALSE) {
  models <- names(effects)
  n <- length(effects)
  M <- matrix(NA_real_, nrow = n, ncol = n,
              dimnames = list(models, models))
  
  for (i in seq_along(effects)) {
    for (j in seq_along(effects)) {
      a <- effects[[i]]
      b <- effects[[j]]
      common <- intersect(names(a), names(b))
      if (length(common) == 0) {
        M[i, j] <- NA
      } else {
        xa <- a[common]
        xb <- b[common]
        if (abs) {
          xa <- abs(xa)
          xb <- abs(xb)
        }
        M[i, j] <- cor(xa, xb, method = method)
      }
    }
  }
  M
}

# Signed effects (direction included)
effect_cor_signed <- effect_correlation_matrix(effects, method = "spearman", abs = FALSE)
effect_cor_signed
# GBM_best     GBMfull       Enet    Enet_60
# GBM_best  1.00000000 -0.03477404  0.4364146  0.4700779
# GBMfull  -0.03477404  1.00000000 -0.2352638 -0.1450465
# Enet      0.43641457 -0.23526384  1.0000000  1.0000000
# Enet_60   0.47007786 -0.14504652  1.0000000  1.0000000

df <- as.data.frame(as.table(effect_cor_signed))  

p<- ggplot(df, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = percent(Freq, accuracy = 1)), 
            color = "white", size = 8) +  # 
  scale_fill_gradient(low = "#D27E64", high = "#2B3256") +
  labs(fill = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank()
  )
p

ggsave(
  filename = file.path(out_dir, "correlation_matrx_4models.png"),
  plot     = p,
  width    = 9,
  height   = 6,
  dpi      = 300
)

## heavy pruning
## --------------
B75 <- B75 %>%
  left_join(md_map, by = "Feature")

eff_B75  <- tapply(B75$weight_signed,B75$Gene,   mean, na.rm = TRUE)

effects <- list(
  GBM_best  = eff_gbm_best,
  GBMfull   = eff_GBMfull,
  Enet      = eff_enet,
  Enet_60   = eff_enet_60,
  B75 = eff_B75
)

effect_correlation_matrix <- function(effects, method = "spearman", abs = FALSE) {
  models <- names(effects)
  n <- length(effects)
  M <- matrix(NA_real_, nrow = n, ncol = n,
              dimnames = list(models, models))
  
  for (i in seq_along(effects)) {
    for (j in seq_along(effects)) {
      a <- effects[[i]]
      b <- effects[[j]]
      common <- intersect(names(a), names(b))
      if (length(common) == 0) {
        M[i, j] <- NA
      } else {
        xa <- a[common]
        xb <- b[common]
        if (abs) {
          xa <- abs(xa)
          xb <- abs(xb)
        }
        M[i, j] <- cor(xa, xb, method = method)
      }
    }
  }
  M
}

# Signed effects (direction included)
effect_cor_signed <- effect_correlation_matrix(effects, method = "spearman", abs = FALSE)
effect_cor_signed
