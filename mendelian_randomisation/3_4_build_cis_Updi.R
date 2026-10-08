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

## Building cis

library(data.table)
library(dplyr)
library(readxl)
library(janitor)
library(stringr)
library(readr)

# PATHS
# -----
proj_root <- PROJECT_ROOT

compiled_file <- file.path(OUTPUT_ROOT, "updi", "Updi_compiled_cis_pqtl",
                           "Updi_compiled_signature_proteins_for_cis_pqtl.csv")

decode_suppl <- file.path(DATA_DIR,
  "deCODE_Supplement.xlsx"
)

out_dir <- file.path(MANUSCRIPT_DIR
)

compiled_proteins <- fread(compiled_file) %>%
  as_tibble()

head(compiled_proteins)
# A tibble: 6 × 11
# protein_id   seq_id   seq_id_decode gene    target                                          uniprot   gbm elasticnet limma n_methods source_methods      
# <chr>        <chr>    <chr>         <chr>   <chr>                                           <chr>   <int>      <int> <int>     <int> <chr>               
#   1 seq.13496.19 13496-19 13496_19      HMGCS1  Hydroxymethylglutaryl-CoA synthase, cytoplasmic Q01581      1          1     1         3 GBM;elasticnet;limma
# 2 seq.13738.8  13738-8  13738_8       INHBA   Inhibin beta A chain                            P08476      1          1     1         3 GBM;elasticnet;limma
# 3 seq.19622.7  19622-7  19622_7       INHBA   Activin A                                       P08476      1          1     1         3 GBM;elasticnet;limma
# 4 seq.20087.3  20087-3  20087_3       SELENOF 15 kDa selenoprotein                            O60613      1          1     1         3 GBM;elasticnet;limma
# 5 seq.2190.55  2190-55  2190_55       F11     Coagulation Factor XI                           P03951      1          1     1         3 GBM;elasticnet;limma
# 6 seq.22124.94 22124-94 22124_94      DYTN    Dystrotelin                                     A2CJ06      1          1     1         3 GBM;elasticnet;limma
names(compiled_proteins)

compiled_proteins %>%
  summarise(
    n_rows = n(),
    n_proteins = n_distinct(protein_id),
    n_seq_id_decode = n_distinct(seq_id_decode),
    n_gbm = sum(gbm == 1),
    n_elasticnet = sum(elasticnet == 1),
    n_limma = sum(limma == 1)
  )

# n_rows n_proteins n_seq_id_decode n_gbm n_elasticnet n_limma
# <int>      <int>           <int> <int>        <int>   <int>
#   1    969        969             969   954           14      51

compiled_proteins %>%
  dplyr::count(seq_id_decode) %>%
  dplyr::filter(n > 1)

#0

#' read decode
#' -----------
#' read decode
#' -----------
tab8 <- read_excel(
  decode_suppl,
  sheet = "ST08",
  skip = 2
) %>%
  clean_names()

head(tab8)
# p_qtl_id_global region_id_protein chrom_variant pos_hg38_variant variant   ld_class_size gene_encoding_target…¹ seq_id uni_prot full_name_of_target_…² rank_cond_sign
# <chr>           <chr>             <chr>                    <dbl> <chr>             <dbl> <chr>                  <chr>  <chr>    <chr>                           <dbl>
#   1 chr1_1          14151_4_chr1_1    chr1                   1013855 rs2465124            60 ISG15                  14151… P05161   Ubiquitin-like protei…              2
# 2 chr1_1          14151_4_chr1_1    chr1                   1014228 rs1921               31 ISG15                  14151… P05161   Ubiquitin-like protei…              1
# 3 chr1_2          15483_377_chr1_1  chr1                   1019790 rs189343…             1 AGRN                   15483… O00468   Agrin                               5
# 4 chr1_2          15483_377_chr1_1  chr1                   1020990 rs568254…             4 AGRN                   15483… O00468   Agrin                               6
# 5 chr1_2          15483_377_chr1_1  chr1                   1027511 rs4970394            11 AGRN                   15483… O00468   Agrin                               1
# 6 chr1_2          15483_377_chr1_1  chr1                   1042813 rs2710873             1 AGRN                   15483… O00468   Agrin                               2
names(tab8)
# [1] "p_qtl_id_global"                        "region_id_protein"                      "chrom_variant"                         
# [4] "pos_hg38_variant"                       "variant"                                "ld_class_size"                         
# [7] "gene_encoding_target_protein"           "seq_id"                                 "uni_prot"                              
# [10] "full_name_of_target_protein"            "rank_cond_sign"                         "count_cond_sign"                       
# [13] "corre_lated_pav_same_gene"              "corre_lated_pav_diff_gene"              "corre_lated_cis_e_qtl_same_gene"       
# [16] "corre_lated_cis_e_qtl_diff_gene"        "log10_p_adj"                            "log10_p_unadj"                         
# [19] "beta_in_sd_adj"                         "beta_in_sd_unadj"                       "strand"                                
# [22] "distance_tss"                           "maf_percent"                            "info"                                  
# [25] "minor_allele"                           "major_allele"                           "prop_cis_e_qtl_same_direction_as_p_qtl"
# [28] "count_cis_e_qtl_same_gene"              "count_cis_e_qtl_diff_gene"              "cis_e_qtl_pav_genes"                   
# [31] "cis_e_qtl_gene"                         "cis_e_qtl_r2"                           "cis_e_qtl_rs_id"                       
# [34] "cis_e_qtl_pval"                         "cis_e_qtl_effect"                       "cis_e_qtl_tissue"                      
# [37] "cis_e_qtl_effect_type"                  "cis_e_qtl_source"                       "pav_effect"                            
# [40] "pav_class_set"                          "pav_r2"                                 "pav_rs_id"                    

#' clean tab 8
tab8_clean <- tab8 %>%
  mutate(
    seq_id = as.character(seq_id),
    variant = as.character(variant),
    chrom_variant = as.character(chrom_variant),
    
    pos_hg38_variant = parse_number(as.character(pos_hg38_variant)),
    log10_p_adj = parse_number(as.character(log10_p_adj)),
    rank_cond_sign = parse_number(as.character(rank_cond_sign)),
    beta_in_sd_adj = parse_number(as.character(beta_in_sd_adj)),
    
    # deCODE gives MAF as percent
    maf_percent = parse_number(as.character(maf_percent)),
    maf = maf_percent / 100,
    
    minor_allele = toupper(str_trim(as.character(minor_allele))),
    major_allele = toupper(str_trim(as.character(major_allele))),
    
    has_bang = if_else(
      minor_allele == "!" | major_allele == "!",
      TRUE,
      FALSE,
      missing = FALSE
    )
  )

#check
tab8_clean %>%
  summarise(
    n_rows = n(),
    n_seq_id = n_distinct(seq_id),
    n_missing_seq_id = sum(is.na(seq_id)),
    n_missing_log10_p_adj = sum(is.na(log10_p_adj)),
    n_missing_beta = sum(is.na(beta_in_sd_adj)),
    n_missing_maf = sum(is.na(maf)),
    n_has_bang = sum(has_bang, na.rm = TRUE)
  )

# n_rows n_seq_id n_missing_seq_id n_missing_log10_p_adj n_missing_beta n_missing_maf n_has_bang
# <int>    <int>            <int>                 <int>          <int>         <int>      <int>
# 7571     1881                0                     0              0             0       1433

#' lets join
#' ---------
# 1. First match by exact deCODE/SomaScan seq_id
# 2. For proteins not matched by seq_id, recover using UniProt + gene agreement
# 3. Keep downstream column names unchanged: protein_id, seq_id, gene, target, uniprot, etc.

compiled_for_join <- compiled_proteins %>%
  dplyr::rename(
    seq_id_compiled = seq_id,
    seq_id_decode_compiled = seq_id_decode,
    gene_compiled = gene,
    target_compiled = target,
    uniprot_compiled = uniprot
  )

# exact SeqID match: same as your original idea
st08_compiled_exact <- tab8_clean %>%
  dplyr::inner_join(
    compiled_for_join,
    by = c("seq_id" = "seq_id_decode_compiled")
  ) %>%
  dplyr::mutate(
    gene = gene_compiled,
    target = target_compiled,
    uniprot = uniprot_compiled,
    match_type = "exact_seq_id"
  )

# proteins still missing after exact SeqID match
compiled_missing_exact <- compiled_for_join %>%
  dplyr::anti_join(
    st08_compiled_exact %>%
      dplyr::distinct(protein_id),
    by = "protein_id"
  )

# fallback: same UniProt and same gene name
st08_compiled_fallback <- tab8_clean %>%
  dplyr::inner_join(
    compiled_missing_exact,
    by = c("uni_prot" = "uniprot_compiled"),
    relationship = "many-to-many"
  ) %>%
  dplyr::filter(
    toupper(gene_encoding_target_protein) == toupper(gene_compiled)
  ) %>%
  dplyr::mutate(
    gene = gene_compiled,
    target = target_compiled,
    uniprot = uni_prot,
    match_type = "uniprot_gene_fallback"
  )

# final merged file
st08_compiled <- dplyr::bind_rows(
  st08_compiled_exact,
  st08_compiled_fallback
)

st08_compiled %>%
  summarise(
    n_pqtl_rows = n(),
    n_proteins_with_any_cis_pqtl = n_distinct(protein_id),
    n_seq_id = n_distinct(seq_id),
    n_gbm = n_distinct(protein_id[gbm == 1]),
    n_elasticnet = n_distinct(protein_id[elasticnet == 1]),
    n_limma = n_distinct(protein_id[limma == 1]),
    n_has_bang = sum(has_bang, na.rm = TRUE),
    n_exact_seq_id = n_distinct(protein_id[match_type == "exact_seq_id"]),
    n_uniprot_gene_fallback = n_distinct(protein_id[match_type == "uniprot_gene_fallback"])
  )

# n_pqtl_rows n_proteins_with_any_cis_pqtl n_seq_id n_gbm n_elasticnet n_limma n_has_bang n_exact_seq_id n_uniprot_gene_fallback
# <int>                        <int>    <int> <int>        <int>   <int>      <int>          <int>                   <int>
#   1        1567                          370      351   361            7      35        293            305                      65

#proteins that did not appear
compiled_no_st08 <- compiled_proteins %>%
  anti_join(
    st08_compiled %>%
    dplyr::distinct(protein_id),
    by = "protein_id"
  )

compiled_no_st08 %>%
  summarise(
    n_missing_from_st08 = n(),
    n_gbm = sum(gbm == 1),
    n_elasticnet = sum(elasticnet == 1),
    n_limma = sum(limma == 1)
  )

# n_missing_from_st08 n_gbm n_elasticnet n_limma
# <int> <int>        <int>   <int>
#   1                 599   593            7      16

#' select significant cis-pQTLs
#' ----------------------------
# ST08 is already cis, but we still keep the deCODE pQTL significance threshold.

pthresh_log10 <- -log10(1.8e-9)

st08_compiled_sig <- st08_compiled %>%
  filter(
    !is.na(log10_p_adj),
    log10_p_adj >= pthresh_log10
  )

st08_compiled_sig %>%
  summarise(
    n_pqtl_rows = n(),
    n_proteins_with_sig_cis_pqtl = n_distinct(protein_id),
    n_seq_id = n_distinct(seq_id),
    n_gbm = n_distinct(protein_id[gbm == 1]),
    n_elasticnet = n_distinct(protein_id[elasticnet == 1]),
    n_limma = n_distinct(protein_id[limma == 1]),
    n_has_bang = sum(has_bang, na.rm = TRUE)
  )

# n_pqtl_rows n_proteins_with_sig_cis_pqtl n_seq_id n_gbm n_elasticnet n_limma n_has_bang
# <int>                        <int>    <int> <int>        <int>   <int>      <int>
#   1        1565                          370      351   361            7      35        292
#' parse PAV r2
#' ------------

parse_r2_max <- function(x) {
  x <- str_trim(as.character(x))
  x[x %in% c("", "-", NA)] <- NA_character_
  
  out <- sapply(x, function(s) {
    if (is.na(s)) return(NA_real_)
    
    nums <- suppressWarnings(as.numeric(str_split(s, ",", simplify = TRUE)))
    nums <- nums[!is.na(nums)]
    
    if (length(nums) == 0) return(NA_real_)
    
    max(nums)
  })
  
  as.numeric(out)
}

st08_compiled_sig <- st08_compiled_sig %>%
  mutate(
    pav_r2_num = parse_r2_max(pav_r2),
    cis_e_qtl_r2_num = parse_r2_max(cis_e_qtl_r2),
    pav_hi = !is.na(pav_r2_num) & pav_r2_num >= 0.8
  )

#' check has_bang patterns by protein
#' ----------------------------------

bang_by_protein <- st08_compiled_sig %>%
  group_by(protein_id) %>%
  summarise(
    n_candidate_iv = n(),
    n_has_bang = sum(has_bang, na.rm = TRUE),
    n_non_bang = sum(!has_bang, na.rm = TRUE),
    only_has_bang = n_non_bang == 0,
    has_clean_option = n_non_bang > 0,
    gbm = max(gbm, na.rm = TRUE),
    elasticnet = max(elasticnet, na.rm = TRUE),
    limma = max(limma, na.rm = TRUE),
    .groups = "drop"
  )

bang_by_protein %>%
  summarise(
    n_proteins = n(),
    n_with_any_has_bang = sum(n_has_bang > 0),
    n_only_has_bang = sum(only_has_bang),
    n_with_clean_option = sum(has_clean_option)
  )

# n_proteins n_with_any_has_bang n_only_has_bang n_with_clean_option
# <int>               <int>           <int>               <int>
#   1        370                 187              25                 345
#' select one cis-IV per protein
#' -----------------------------
# Rule:
# 1. If a protein has any non-! cis-pQTL, prefer non-! variants.
# 2. If a protein only has ! variants, keep the best ! variant.
# 3. Prefer rank 1 conditional signal.
# 4. Prefer non-PAV.
# 5. Then choose strongest pQTL by log10_p_adj.

st08_iv_one <- st08_compiled_sig %>%
  group_by(protein_id) %>%
  mutate(
    n_candidate_iv = n(),
    n_non_bang = sum(!has_bang, na.rm = TRUE),
    only_has_bang = n_non_bang == 0,
    
    is_rank1 = !is.na(rank_cond_sign) & rank_cond_sign == 1,
    
    iv_priority = case_when(
      # Best: clean, rank 1, non-PAV
      !has_bang & is_rank1 & !pav_hi ~ 1L,
      
      # Clean, rank 1, but PAV
      !has_bang & is_rank1 & pav_hi ~ 2L,
      
      # Clean, non-PAV, but not rank 1
      !has_bang & !pav_hi ~ 3L,
      
      # Any other clean variant
      !has_bang ~ 4L,
      
      # If only bang variants exist, keep them as fallback
      has_bang & only_has_bang & is_rank1 & !pav_hi ~ 5L,
      has_bang & only_has_bang & is_rank1 & pav_hi ~ 6L,
      has_bang & only_has_bang & !pav_hi ~ 7L,
      has_bang & only_has_bang ~ 8L,
      
      TRUE ~ 99L
    )
  ) %>%
  arrange(
    iv_priority,
    desc(log10_p_adj),
    .by_group = TRUE
  ) %>%
  slice_head(n = 1) %>%
  ungroup()

#' selected IV QC
#' --------------
st08_iv_one %>%
  summarise(
    n_iv_rows = n(),
    n_proteins = n_distinct(protein_id),
    n_seq_id = n_distinct(seq_id),
    n_gbm = n_distinct(protein_id[gbm == 1]),
    n_elasticnet = n_distinct(protein_id[elasticnet == 1]),
    n_limma = n_distinct(protein_id[limma == 1]),
    n_duplicate_protein_id = sum(duplicated(protein_id)),
    n_selected_has_bang = sum(has_bang, na.rm = TRUE),
    n_selected_only_has_bang = sum(only_has_bang, na.rm = TRUE),
    n_selected_pav_hi = sum(pav_hi, na.rm = TRUE),
    n_selected_rank1 = sum(rank_cond_sign == 1, na.rm = TRUE)
  )

# n_iv_rows n_proteins n_seq_id n_gbm n_elasticnet n_limma n_duplicate_protein_id n_selected_has_bang n_selected_only_has_bang n_selected_pav_hi n_selected_rank1
# <int>      <int>    <int> <int>        <int>   <int>                  <int>               <int>                    <int>             <int>            <int>
#   1       370        370      350   361            7      35                      0                  25                       25                73              242
# > 
#   > 
#count rows that need filling
iv_needs_fill <- st08_iv_one %>%
  mutate(
    needs_allele_fill = has_bang |
      is.na(minor_allele) | minor_allele == "" |
      is.na(major_allele) | major_allele == "" |
      minor_allele == major_allele
  ) %>%
  filter(needs_allele_fill)

iv_needs_fill %>%
  dplyr::select(
    protein_id,
    seq_id,
    gene,
    target,
    variant,
    chrom_variant,
    pos_hg38_variant,
    minor_allele,
    major_allele,
    has_bang,
    log10_p_adj,
    rank_cond_sign,
    pav_hi,
    gbm,
    elasticnet,
    limma
  ) %>%
  arrange(protein_id)

iv_needs_fill %>%
  summarise(
    n_needs_fill = n(),
    n_proteins_needing_fill = n_distinct(protein_id),
    n_has_bang = sum(has_bang, na.rm = TRUE),
    n_same_allele = sum(minor_allele == major_allele, na.rm = TRUE)
  )

# n_needs_fill n_proteins_needing_fill n_has_bang n_same_allele
# <int>                   <int>      <int>         <int>
#   1           25                      25         25             0

st08_iv_one <- st08_iv_one %>%
  mutate(
    needs_allele_fill = has_bang |
      is.na(minor_allele) | minor_allele == "" |
      is.na(major_allele) | major_allele == "" |
      minor_allele == major_allele
  )

#' prepare MR-ready cis-IV table before allele filling
#' ---------------------------------------------------

st08_iv_mr <- st08_iv_one %>%
  mutate(
    needs_allele_fill = has_bang |
      is.na(minor_allele) | minor_allele == "" |
      is.na(major_allele) | major_allele == "" |
      minor_allele == major_allele,
    
    rsid = str_trim(sub(",.*$", "", variant)),
    
    variant_pos_hg38 = paste0(chrom_variant, ":", pos_hg38_variant),
    
    # useful for matching against files using chr:pos:allele1:allele2
    variant_key_minor_major = if_else(
      !needs_allele_fill,
      paste0(chrom_variant, ":", pos_hg38_variant, ":", minor_allele, ":", major_allele),
      NA_character_
    ),
    
    variant_key_major_minor = if_else(
      !needs_allele_fill,
      paste0(chrom_variant, ":", pos_hg38_variant, ":", major_allele, ":", minor_allele),
      NA_character_
    ),
    
    # deCODE beta is for the minor allele.
    # For rows needing allele fill, leave MR alleles as NA for now.
    effect_allele.exposure = if_else(
      needs_allele_fill,
      NA_character_,
      minor_allele
    ),
    
    other_allele.exposure = if_else(
      needs_allele_fill,
      NA_character_,
      major_allele
    ),
    
    beta.exposure = beta_in_sd_adj,
    eaf.exposure = maf,
    pval.exposure = 10^(-log10_p_adj),
    
    pval.exposure = if_else(
      pval.exposure == 0,
      .Machine$double.xmin,
      pval.exposure
    ),
    
    SNP = rsid,
    exposure = paste0(gene, "_", protein_id),
    id.exposure = protein_id
  ) %>%
  dplyr::select(
    SNP,
    exposure,
    id.exposure,
    
    protein_id,
    seq_id,
    gene,
    target,
    uniprot,
    
    gbm,
    elasticnet,
    limma,
    n_methods,
    source_methods,
    
    chrom_variant,
    pos_hg38_variant,
    variant_pos_hg38,
    variant,
    variant_key_minor_major,
    variant_key_major_minor,
    
    minor_allele,
    major_allele,
    has_bang,
    needs_allele_fill,
    
    beta.exposure,
    effect_allele.exposure,
    other_allele.exposure,
    eaf.exposure,
    pval.exposure,
    
    log10_p_adj,
    beta_in_sd_adj,
    maf,
    maf_percent,
    
    rank_cond_sign,
    count_cond_sign,
    ld_class_size,
    distance_tss,
    
    pav_r2,
    pav_r2_num,
    pav_hi,
    cis_e_qtl_r2,
    cis_e_qtl_r2_num,
    
    everything()
  )
#check
st08_iv_mr %>%
  summarise(
    n_rows = n(),
    n_proteins = n_distinct(protein_id),
    n_snps = n_distinct(SNP),
    n_needs_allele_fill = sum(needs_allele_fill, na.rm = TRUE),
    n_missing_effect_allele = sum(is.na(effect_allele.exposure)),
    n_missing_other_allele = sum(is.na(other_allele.exposure)),
    n_missing_beta = sum(is.na(beta.exposure)),
    n_missing_pval = sum(is.na(pval.exposure)),
    n_missing_eaf = sum(is.na(eaf.exposure))
  )
# n_rows n_proteins n_snps n_needs_allele_fill n_missing_effect_allele n_missing_other_allele n_missing_beta n_missing_pval n_missing_eaf
# <int>      <int>  <int>               <int>                   <int>                  <int>          <int>          <int>         <int>
#   1    370        370    347                  25                      25                     25              0              0             0
 

fwrite(
  st08_iv_mr,
  file.path(out_dir, "Updi_proteins_one_cis_IV_per_protein_MR_RAW_1.csv")
)
