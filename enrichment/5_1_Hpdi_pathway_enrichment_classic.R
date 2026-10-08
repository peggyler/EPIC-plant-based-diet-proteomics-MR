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

#'---- Syntax for ProPlant Project'#
#' Created by Peggy Ler (Pler) - 31 Oct 2025
#' Filename..: hPDI_protein_enrichment_classic.R
#' Updated...: "02 Nov 2025" by "Grok-assisted"
#' Purpose...: Protein enrichment with Kegg, reactome, GO, _ over reprepsrentation enrichment instead
#' of gene set enrichment
#' 
#' --------------------------------------------------------------
#' Using R studio 2024.04.02
#' -- Main program ----------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(tidyr)
  library(stringr)
  library(enrichplot)   # for cnetplot, heatplot, emapplot, gseaplot2
  library(pathview)     # for coloured KEGG diagrams (install via BiocManager if missing)
  library(ComplexHeatmap)  # for combined heatmap (install via BiocManager)
})

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

# gene names uniprot, info
md_file   <- file.path(file.path(REFERENCE_DIR, "MetaData_Soma_SMP.rds"))


# Background universe of proteins
# --------------------------------
# The full hPDI LIMMA export is used only as an inventory of assayed protein IDs.
# No hPDI effect estimates, p-values or significance cutoffs define the background.
# The shared background is subsequently mapped to unique non-missing ENTREZ IDs.
allprot <- read.csv(file.path(OUTPUT_ROOT, paste0("hpdi/Hpidnovegoils_limma_", LIMMA_DATE, ".csv")))

head(allprot)
background_proteins <- allprot$protein
# protein        logFC       AveExpr          t   P.Value adj.P.Val         B      q_BH q_BY p_bonf
# 1 seq.10000.28  0.012440795  0.0089120576  0.8670334 0.3859711 0.6386133 -6.319731 0.6386133    1      1
# 2  seq.10001.7 -0.006592828 -0.0118054311 -0.5810348 0.5612469 0.7662267 -6.523687 0.7662267    1      1
# 3 seq.10003.15  0.013748213  0.0281393899  1.0988102 0.2719112 0.5342793 -6.095364 0.5342793    1      1
# 4 seq.10006.25  0.012120783 -0.0006398676  0.9381828 0.3482021 0.6061296 -6.256481 0.6061296    1      1
# 5 seq.10008.43  0.002451172 -0.0015425746  0.2042577 0.8381616 0.9192845 -6.669412 0.9192845    1      1
# 6 seq.10010.10  0.010955506  0.0263789595  0.8609003 0.3893400 0.6409208 -6.324949 0.6409208    1      1
length(background_proteins) #7285

# Output folder
# -------------
out <- file.path(proj_root, "protein_enrichment")
dir.create(out, showWarnings = FALSE, recursive = TRUE)


# Load metadata
# -------------
md <- readRDS(md_file)

md_map <- md %>%
  dplyr::select(AptName, EntrezGeneSymbol, UniProt) %>%
  dplyr::distinct() %>%
  dplyr::rename(
    protein = AptName,
    Gene = EntrezGeneSymbol
  )


# Build background mapping: protein_id -> gene symbol -> ENTREZID
# ---------------------------------------------------------------
bg_genes <- md_map %>%
  dplyr::filter(protein %in% background_proteins)

df_long <- bg_genes %>%
  dplyr::mutate(Gene = dplyr::if_else(is.na(Gene) | Gene == "", protein, Gene)) %>%
  dplyr::mutate(Gene_list = stringr::str_split(Gene, "\\|")) %>%
  tidyr::unnest(Gene_list) %>%
  dplyr::mutate(Gene = stringr::str_trim(Gene_list)) %>%
  dplyr::select(-Gene_list) %>%
  dplyr::distinct(protein, Gene, .keep_all = TRUE)

if ("UniProt" %in% names(df_long)) {
  df_long <- df_long %>%
    dplyr::mutate(
      UniProt = ifelse(
        is.na(UniProt),
        NA_character_,
        sub(";.*$", "", UniProt)
      )
    )
}

# Map gene symbols to ENTREZ IDs
sym_keys <- unique(df_long$Gene)

entrez_from_symbol <- suppressMessages(
  AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = sym_keys,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )
)

df_long$ENTREZID <- unname(entrez_from_symbol[df_long$Gene])

# Fallback: map UniProt to ENTREZ where symbol mapping failed
need_fallback <- is.na(df_long$ENTREZID) & !is.na(df_long$UniProt)

if (any(need_fallback)) {
  up_keys <- unique(df_long$UniProt[need_fallback])
  
  entrez_from_uniprot <- suppressMessages(
    AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = up_keys,
      keytype = "UNIPROT",
      column = "ENTREZID",
      multiVals = "first"
    )
  )
  
  df_long$ENTREZID[need_fallback] <-
    unname(entrez_from_uniprot[df_long$UniProt[need_fallback]])
}

cat(
  "Mapped background proteins to ENTREZ: ",
  round(mean(!is.na(df_long$ENTREZID)) * 100, 2),
  "%\n",
  sep = ""
)

bg_entrez <- df_long %>%
  dplyr::filter(!is.na(ENTREZID)) %>%
  dplyr::pull(ENTREZID) %>%
  unique()

length(bg_entrez)


# Prepare signature list
# ----------------------
# All three files use protein_id, based on your screenshots.
sig_list <- list(
  "hPDI-limma" = hPDI_limma$protein_id,
  "hPDI-ENET"  = hPDI_ENET$protein_id,
  "hPDI-XGB"   = hPDI_XGB$protein_id
)

sig_list <- lapply(sig_list, unique)

sapply(sig_list, length)
# hPDI-limma  hPDI-ENET   hPDI-XGB 
# 619        108        667 

# Function: map signature proteins to ENTREZ
# ------------------------------------------
map_signature_to_entrez <- function(protein_vec, signature_name) {
  
  mapped <- tibble::tibble(protein = unique(protein_vec)) %>%
    dplyr::left_join(
      df_long %>%
        dplyr::select(protein, Gene, UniProt, ENTREZID) %>%
        dplyr::distinct(),
      by = "protein"
    )
  
  cat(
    signature_name, ": ",
    length(unique(protein_vec)), " proteins; ",
    length(unique(mapped$ENTREZID[!is.na(mapped$ENTREZID)])), " unique ENTREZ IDs mapped\n",
    sep = ""
  )
  
  mapped
}


mapped_signatures <- purrr::imap(
  sig_list,
  ~ map_signature_to_entrez(.x, .y)
)


# Function: run enrichment
# ------------------------
run_enrichment <- function(hit_entrez, signature_name) {
  
  hit_entrez <- unique(hit_entrez[!is.na(hit_entrez)])
  
  if (length(hit_entrez) < 5) {
    warning(paste("Too few mapped genes for", signature_name))
    return(NULL)
  }
  
  # GO
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
      dplyr::mutate(
        source = source_name,
        signature = signature_name,
        n_signature_entrez = length(hit_entrez)
      )
  }
  
  dplyr::bind_rows(
    add_source(ego_bp, "GO_BP"),
    add_source(ego_cc, "GO_CC"),
    add_source(ego_mf, "GO_MF"),
    add_source(ekegg, "KEGG"),
    add_source(ereact, "Reactome")
  )
}


# Run enrichment for each hPDI signature
# --------------------------------------
enrichment_all <- purrr::imap_dfr(
  mapped_signatures,
  function(mapped_df, sig_name) {
    hit_entrez <- mapped_df %>%
      dplyr::filter(!is.na(ENTREZID)) %>%
      dplyr::pull(ENTREZID) %>%
      unique()
    
    run_enrichment(hit_entrez, sig_name)
  }
)


# Significant enrichment results
# ------------------------------
enrichment_sig <- enrichment_all %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(signature, p.adjust)


# Save full and significant results
# ---------------------------------
readr::write_csv(
  enrichment_all,
  file.path(out, "hPDI_pathway_enrichment_all_results.csv")
)

readr::write_csv(
  enrichment_sig,
  file.path(out, "hPDI_pathway_enrichment_FDR005.csv")
)

table(enrichment_sig$signature)
# hPDI-ENET hPDI-limma   hPDI-XGB 
# 3         27          4 

saveRDS(
  enrichment_sig,
  file.path(out, "hPDI_pathway_enrichment_FDR005.rds")
)

saveRDS(
  enrichment_all,
  file.path(out, "hPDI_pathway_enrichment_all_results.rds")
)

#'PLOTTING
#'--------------
enrichment_sig <- read.csv(file.path(out, "hPDI_pathway_enrichment_FDR005.csv"))

# Helper to convert GeneRatio to numeric
gene_ratio_to_num <- function(x) {
  sapply(strsplit(x, "/"), function(z) as.numeric(z[1]) / as.numeric(z[2]))
}

hPDI_top20_each <- enrichment_sig %>%
  dplyr::mutate(
    signature = factor(
      signature,
      levels = c("hPDI-limma", "hPDI-ENET", "hPDI-XGB")
    )
  ) %>%
  dplyr::group_by(signature) %>%
  dplyr::arrange(p.adjust, .by_group = TRUE) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    GeneRatio_num = gene_ratio_to_num(GeneRatio),
    logFDR = -log10(p.adjust),
    Description_short = stringr::str_trunc(Description, width = 70)
  )

p_hPDI_top20_each <- ggplot(
  hPDI_top20_each,
  aes(
    x = GeneRatio_num,
    y = reorder(Description_short, GeneRatio_num)
  )
) +
  geom_point(aes(size = logFDR, colour = source)) +
  facet_grid(
    rows = vars(signature),
    scales = "free_y",
    space = "free_y"
  ) +
  scale_size_continuous(name = "-log10(FDR)") +
  scale_color_manual(
    name = "Database / Ontology",
    values = c(
      "GO_MF" = "#D27E64",
      "KEGG" = "#D3D0A9",
      "GO_BP" = "#E3AE9D",
      "GO_CC" = "#2B3256",
      "Reactome" = "#8C8C8C"
    ),
    labels = c(
      "GO_BP" = "GO BP",
      "GO_CC" = "GO CC",
      "GO_MF" = "GO MF",
      "KEGG" = "KEGG",
      "Reactome" = "Reactome"
    )
  ) +
  labs(
    x = "Gene ratio",
    y = NULL,
    title = "hPDI – Top enriched pathways by protein signature"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 10),
    strip.text.y = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5)
  )

p_hPDI_top20_each

ggsave(
  filename = file.path(out, "hPDI_top20_enriched_pathways_by_signature.png"),
  plot = p_hPDI_top20_each,
  width = 9,
  height = 13,
  dpi = 600
)
