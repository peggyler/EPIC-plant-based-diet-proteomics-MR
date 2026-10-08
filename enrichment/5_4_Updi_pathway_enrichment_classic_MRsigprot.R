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
# Pathway enrichment for MR-significant causal proteins
# Peggy
#
# Analyses:
#   1) Hits = MR-significant proteins
#      Universe = all cis proteins from LIMMA ∪ GBM
#
#   2) Hits = MR-significant proteins
#      Universe = all signature proteins from LIMMA ∪ GBM
#
# Output:
#   - duplicate MR protein summary
#   - mapped hit/background tables
#   - enrichment tables
#   - combined dotplots
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
proj_root <- PROJECT_ROOT

# MR-significant proteins (hits)
mr_dir  <- file.path(MANUSCRIPT_DIR)
mr_file <- file.path(mr_dir, "MR_all_source_sigresults_with_coloc.csv")

# cis universe files
cis_file    <- file.path(mr_dir, "Updi_proteins_one_cis_IV_per_protein_MR_MANUAL_v1.csv")


# full signature universe files
# sig_out_dir     <- file.path(proj_root, "outputs/legacy_subcohort/Hpdi_novegoils")
# gbmfull_file    <- file.path(sig_out_dir, "Sig_GB_lowestRMSE_full.csv")
# limmafull_file  <- file.path(sig_out_dir, "Hpidnovegoils_limma_epic_20251107.csv")
# 
# # SomaLogic annotation
# md_file <- file.path(REFERENCE_DIR, "MetaData_Soma_SMP.rds")
# 
# # output directory
# out_dir <- file.path(MANUSCRIPT_DIR, "pathwayenrichment")
# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
sigprot <- read.csv(mr_file, stringsAsFactors = FALSE)
names(sigprot)
head(sigprot)

cis <- read.csv(cis_file, stringsAsFactors = FALSE)
names(cis)
# [1] "seq_id"           "gene"             "protein"          "UniProt"          "rsid"             "chrom"            "pos_hg38"        
# [8] "beta_tab8"        "p"                "pav_r2"           "iv_pav_flag"      "minor_allele"     "major_allele"     "has_bang"        
# [15] "protein_id"       "weight"           "rsid_first"       "name1"            "name2"            "Chrom"            "Pos"             
# [22] "Name"             "effectAllele"     "otherAllele"      "effectAlleleFreq" "beta_aligned"    

# cleaning up MR hit list
mr_hits <- sigprot %>%
  dplyr::filter(pdi_type == "uPDI")%>%
  dplyr::select(Gene = exposure) %>%
  dplyr::filter(!is.na(Gene), Gene != "") %>%
  dplyr::distinct(Gene)

nrow(mr_hits) #41
head(mr_hits)

#cleaning up cis universe
cis_universe <- cis %>%
  dplyr::select(Gene = gene)



nrow(cis_universe) #370
head(cis_universe)


# background mapping - CIS
#-------------------
bg_gene_symbols <- unique(cis_universe$Gene)

entrez_from_symbol <- suppressMessages(
  AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = bg_gene_symbols,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )
)

bg_entrez <- unname(entrez_from_symbol[bg_gene_symbols])
bg_entrez <- unique(bg_entrez[!is.na(bg_entrez)])

length(bg_entrez) #345
head(bg_entrez)

# hits mapping
hit_gene_symbols <- unique(mr_hits$Gene)

hit_entrez <- suppressMessages(
  AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = hit_gene_symbols,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )
)

hit_entrez <- unname(hit_entrez[hit_gene_symbols])
hit_entrez <- unique(hit_entrez[!is.na(hit_entrez)])

length(hit_entrez)#41
head(hit_entrez)

#enrich
#---------
# Function: run enrichment
# ------------------------
run_cis_mr_enrichment <- function(hit_entrez, bg_entrez, signature_name = "MR-cis-hits") {
  
  hit_entrez <- unique(hit_entrez[!is.na(hit_entrez)])
  bg_entrez  <- unique(bg_entrez[!is.na(bg_entrez)])
  
  if (length(hit_entrez) < 5) {
    warning("Too few mapped genes for enrichment")
    return(NULL)
  }
  
  # GO Biological Process
  ego_bp <- enrichGO(
    gene = hit_entrez,
    universe = bg_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )
  
  # GO Cellular Component
  ego_cc <- enrichGO(
    gene = hit_entrez,
    universe = bg_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )
  
  # GO Molecular Function
  ego_mf <- enrichGO(
    gene = hit_entrez,
    universe = bg_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )
  
  # KEGG
  ekegg <- enrichKEGG(
    gene = hit_entrez,
    universe = bg_entrez,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
  
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    ekegg <- setReadable(
      ekegg,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID"
    )
  }
  
  # Reactome
  ereact <- enrichPathway(
    gene = hit_entrez,
    universe = bg_entrez,
    organism = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )
  
  add_source <- function(enrich_obj, source_name) {
    if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
      return(NULL)
    }
    
    as.data.frame(enrich_obj) %>%
      mutate(
        source = source_name,
        signature = signature_name,
        n_hit_entrez = length(hit_entrez),
        n_background_entrez = length(bg_entrez)
      )
  }
  
  bind_rows(
    add_source(ego_bp, "GO_BP"),
    add_source(ego_cc, "GO_CC"),
    add_source(ego_mf, "GO_MF"),
    add_source(ekegg, "KEGG"),
    add_source(ereact, "Reactome")
  )
}

#run enrichment
enrichment_all <- run_cis_mr_enrichment(
  hit_entrez = hit_entrez,
  bg_entrez = bg_entrez,
  signature_name = "MR-cis-hits"
)

enrichment_all <- enrichment_all %>%
  mutate(
    Description = str_replace(
      Description,
      "^([a-z])",
      ~ str_to_upper(.x)
    )
  )

enrichment_sig <- enrichment_all %>%
  filter(p.adjust < 0.05) %>%
  arrange(p.adjust)

table(enrichment_sig$source)

#save
write_csv(enrichment_all,file.path(mr_dir, "Updi_MRsig_cisuniverse_pathway_enrichment_all_results.csv"))

# Helper to convert GeneRatio to numeric
# gene_ratio_to_num <- function(x) {
#   sapply(strsplit(x, "/"), function(z) as.numeric(z[1]) / as.numeric(z[2]))
# }
# 
# cis_top20 <- enrichment_sig %>%
#   arrange(p.adjust) %>%
#   slice_head(n = 20) %>%
#   mutate(
#     GeneRatio_num = gene_ratio_to_num(GeneRatio),
#     logFDR = -log10(p.adjust),
#     Description_short = str_trunc(Description, width = 70),
#     source = factor(
#       source,
#       levels = c("GO_BP", "GO_CC", "GO_MF", "KEGG", "Reactome")
#     )
#   )
# 
# p_cis_top20 <- ggplot(
#   cis_top20,
#   aes(
#     x = GeneRatio_num,
#     y = reorder(Description_short, GeneRatio_num)
#   )
# ) +
#   geom_point(aes(size = logFDR, colour = source)) +
#   scale_size_continuous(name = "-log10(FDR)") +
#   scale_color_manual(
#     name = "Database / Ontology",
#     values = c(
#       "GO_MF" = "#D27E64",
#       "KEGG" = "#D3D0A9",
#       "GO_BP" = "#E3AE9D",
#       "GO_CC" = "#2B3256",
#       "Reactome" = "#8C8C8C"
#     ),
#     labels = c(
#       "GO_BP" = "GO BP",
#       "GO_CC" = "GO CC",
#       "GO_MF" = "GO MF",
#       "KEGG" = "KEGG",
#       "Reactome" = "Reactome"
#     )
#   ) +
#   labs(
#     x = "Gene ratio",
#     y = NULL,
#     title = "MR-significant proteins – Enriched pathways"
#   ) +
#   theme_bw(base_size = 10) +
#   theme(
#     axis.text.y = element_text(size = 10),
#     plot.title = element_text(hjust = 0.5)
#   )
# 
# p_cis_top20
# 
# ggsave(
#   filename = file.path(mr_dir, "MRsig_hPDI_enriched_pathways.png"),
#   plot = p_cis_top20,
#   width = 9,
#   height = 7,
#   dpi = 600
# )
# 
library(openxlsx)


# Make sure ratio columns are character
enrichment_all_xlsx <- enrichment_all %>%
  mutate(
    GeneRatio = as.character(GeneRatio),
    BgRatio   = as.character(BgRatio)
  )

enrichment_sig_xlsx <- enrichment_sig %>%
  mutate(
    GeneRatio = as.character(GeneRatio),
    BgRatio   = as.character(BgRatio)
  )

# Create workbook
wb <- createWorkbook()

addWorksheet(wb, "All_results")
addWorksheet(wb, "FDR_0.05")

writeData(wb, "All_results", enrichment_all_xlsx)
writeData(wb, "FDR_0.05", enrichment_sig_xlsx)

# Force GeneRatio and BgRatio columns to text format
text_style <- createStyle(numFmt = "@")

ratio_cols_all <- which(names(enrichment_all_xlsx) %in% c("GeneRatio", "BgRatio"))
ratio_cols_sig <- which(names(enrichment_sig_xlsx) %in% c("GeneRatio", "BgRatio"))

addStyle(
  wb,
  sheet = "All_results",
  style = text_style,
  cols = ratio_cols_all,
  rows = 1:(nrow(enrichment_all_xlsx) + 1),
  gridExpand = TRUE
)

addStyle(
  wb,
  sheet = "FDR_0.05",
  style = text_style,
  cols = ratio_cols_sig,
  rows = 1:(nrow(enrichment_sig_xlsx) + 1),
  gridExpand = TRUE
)

# Optional: freeze header and adjust width
freezePane(wb, "All_results", firstRow = TRUE)
freezePane(wb, "FDR_0.05", firstRow = TRUE)

setColWidths(wb, "All_results", cols = 1:ncol(enrichment_all_xlsx), widths = "auto")
setColWidths(wb, "FDR_0.05", cols = 1:ncol(enrichment_sig_xlsx), widths = "auto")

# Save
saveWorkbook(
  wb,
  file = file.path(mr_dir, "UPDI_MRsig_CIS_pathway_enrichment_results.xlsx"),
  overwrite = TRUE
)
