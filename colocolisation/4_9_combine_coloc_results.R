# File: 4_9_combine_coloc_results.R
# Analysis: Combine colocalisation summaries with MR results and create figures.
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
library(purrr)
library(readr)
library(stringr)
library(readxl)

coloc_dir <- file.path(MANUSCRIPT_DIR, "results_coloc_local")

# only CVD coloc-local csv files
files_cvd <- list.files(
  path = coloc_dir,
  pattern = "_CVD_coloc_local.*\\.csv$",
  full.names = TRUE
)

files_cvd
length(files_cvd)

df_coloc_cvd <- map_dfr(files_cvd, function(f) {
  read_csv(f, show_col_types = FALSE) %>%
    mutate(
      source_file = basename(f),
      protein_from_file = str_match(basename(f), "^\\d+_(.+?)_CVD_coloc_local")[,2]
    )
})

df_coloc_cvd

#check
names(df_coloc_cvd)
table(df_coloc_cvd$disease) # CVD 25 rows  
df_coloc_cvd %>% count(protein, sort = TRUE)

## only t2d coloc-local csv files
##--------------------------------
coloc_dir <- file.path(MANUSCRIPT_DIR, "results_coloc_local/t2d")
files_t2d <- list.files(
  path = coloc_dir,
  pattern = "_T2D_coloc_local.*\\.csv$",
  full.names = TRUE
)

files_t2d
length(files_t2d)

df_coloc_t2d <- map_dfr(files_t2d, function(f) {
  read_csv(f, show_col_types = FALSE) %>%
    mutate(
      source_file = basename(f),
      protein_from_file = str_match(basename(f), "^\\d+_(.+?)_T2D_coloc_local")[,2]
    )
})

df_coloc_t2d

#check
names(df_coloc_t2d)
table(df_coloc_t2d$disease) # t2d 12 rows  
df_coloc_t2d %>% count(protein, sort = TRUE)

## only all cancer coloc-local csv files
##--------------------------------
coloc_dir <- file.path(MANUSCRIPT_DIR, "results_coloc_local/allcancer")
files_ca <- list.files(
  path = coloc_dir,
  pattern = "_cancer_coloc_local.*\\.csv$",
  full.names = TRUE
)

files_ca
length(files_ca)

df_coloc_ca <- map_dfr(files_ca, function(f) {
  read_csv(f, show_col_types = FALSE) %>%
    mutate(
      disease = "All cancer",
      source_file = basename(f),
      protein_from_file = str_match(basename(f), "^\\d+_(.+?)_cancer_coloc_local")[,2]
    )
})

df_coloc_ca

#check
names(df_coloc_ca)
table(df_coloc_ca$disease) # t2d 11 rows  
df_coloc_ca %>% count(protein, sort = TRUE)

##LUNG
coloc_dir <- file.path(MANUSCRIPT_DIR, "results_coloc_local/lung")
files_lung <- list.files(
  path = coloc_dir,
  pattern = "_Lung_cancer_coloc_.*\\.csv$",
  full.names = TRUE
)

files_lung
length(files_lung)

df_coloc_lung <- map_dfr(files_lung, function(f) {
  read_csv(f, show_col_types = FALSE) %>%
    mutate(
      source_file = basename(f),
      protein_from_file = str_match(basename(f), "^\\d+_(.+?)_coloc_lung")[,2]
    )
})

df_coloc_lung

#check
names(df_coloc_lung)
table(df_coloc_lung$disease) # t2d 11 rows  
df_coloc_lung %>% count(protein, sort = TRUE)

##breast
##--------
coloc_dir <- file.path(MANUSCRIPT_DIR, "results_coloc_local/breast")
files_breast <- list.files(
  path = coloc_dir,
  pattern = "_Breast_cancer_coloc_local.*\\.csv$",
  full.names = TRUE
)

files_breast
length(files_breast)

df_coloc_breast <- map_dfr(files_breast, function(f) {
  read_csv(f, show_col_types = FALSE) %>%
    mutate(
      source_file = basename(f),
      protein_from_file = str_match(basename(f), "^\\d+_(.+?)_Breast_cancer_coloc_local")[,2]
    )
})

df_coloc_breast

#check
names(df_coloc_breast)
table(df_coloc_breast$disease) # 7 rows  
df_coloc_breast %>% count(protein, sort = TRUE)

##prostate
##--------
coloc_dir <- file.path(MANUSCRIPT_DIR, "results_coloc_local/prostate")
files_prostate <- list.files(
  path = coloc_dir,
  pattern = "_coloc.*\\.csv$",
  full.names = TRUE
)

files_prostate
length(files_prostate)

df_coloc_prostate <- map_dfr(files_prostate, function(f) {
  read_csv(f, show_col_types = FALSE) %>%
    mutate(
      source_file = basename(f),
      protein_from_file = str_match(basename(f), "^\\d+_(.+?)_coloc")[,2]
    )
})

df_coloc_prostate

#check
names(df_coloc_prostate)
table(df_coloc_prostate$disease) # 16 rows  
df_coloc_prostate %>% count(protein, sort = TRUE)

##endometrial cancer
##--------
coloc_dir <- file.path(MANUSCRIPT_DIR, "results_coloc_local/endom")
files_endom <- list.files(
  path = coloc_dir,
  pattern = "_coloc.*\\.csv$",
  full.names = TRUE
)

files_endom 
length(files_endom )

df_coloc_endom <- map_dfr(files_endom, function(f) {
  read_csv(f, show_col_types = FALSE) %>%
    mutate(
      source_file = basename(f),
      protein_from_file = str_match(basename(f), "^\\d+_(.+?)_coloc")[,2]
    )
})

df_coloc_endom

#check
names(df_coloc_endom)
table(df_coloc_endom$disease) # 16 rows  
df_coloc_endom %>% count(protein, sort = TRUE)

##colorectal
##--------
coloc_dir <- file.path(MANUSCRIPT_DIR, "results_coloc_local/colorectal")
files_colo <- list.files(
  path = coloc_dir,
  pattern = "_coloc.*\\.csv$",
  full.names = TRUE
)

files_colo 
length(files_colo)

df_coloc_colorectal <- map_dfr(files_colo, function(f) {
  read_csv(f, show_col_types = FALSE) %>%
    mutate(
      source_file = basename(f),
      protein_from_file = str_match(basename(f), "^\\d+_(.+?)_coloc")[,2]
    )
})

df_coloc_colorectal

##liver cancer
##--------
coloc_dir <- file.path(MANUSCRIPT_DIR, "results_coloc_local/liver")
files_liver <- list.files(
  path = coloc_dir,
  pattern = "_coloc.*\\.csv$",
  full.names = TRUE
)

files_liver 
length(files_liver)

df_coloc_liver <- map_dfr(files_liver, function(f) {
  read_csv(f, show_col_types = FALSE) %>%
    mutate(
      source_file = basename(f),
      protein_from_file = str_match(basename(f), "^\\d+_(.+?)_coloc")[,2]
    )
})

df_coloc_liver

#check
names(df_coloc_liver)
table(df_coloc_liver$disease) # 16 rows  
df_coloc_liver %>% count(protein, sort = TRUE)


##row bind
df_coloc_all <- bind_rows(
  df_coloc_cvd,
  df_coloc_t2d,
  df_coloc_ca,
  df_coloc_lung,
  df_coloc_breast,
  df_coloc_prostate,
  df_coloc_endom,
  df_coloc_liver,
  df_coloc_colorectal
)

head(df_coloc_all)
table(df_coloc_all$coloc_interpretation)
# Moderate colocalization   Strong colocalization  Weak/no colocalization 
# 12                      15                      68 

df_coloc_all %>% filter(coloc_interpretation == "Strong colocalization")
# protein disease gwas_id pqtl_file outcome_file region_type region cis_window nsnps_region nsnps_overlap     PP.H0
# <chr>   <chr>   <chr>   <chr>     <chr>        <chr>       <chr>       <dbl>        <dbl>         <dbl>     <dbl>
#   1 PCSK9   CVD     ebi-a-… 5231_79_… ebi-a-GCST9… cis_anchor  chr1:…     500000         8954          4852 6.00e-113
# 2 LRIG1   T2D     ebi-a-… 18831_6_… ebi-a-GCST0… cis_anchor  chr3:…     500000         9252          1593 0        
# 3 NCAN    T2D     ebi-a-… 15573_11… ebi-a-GCST0… cis_anchor  chr19…     500000         7376           889 8.67e- 93
# 4 MANBA   T2D     ebi-a-… 6382_17_… ebi-a-GCST0… cis_anchor  chr4:…     500000         7218          1438 0        
# 5 KRT5    All ca… ieu-b-… download… ieu-b-4966.… cis_anchor  chr12…     500000         9088          4713 4.25e- 58
# 6 CPNE1   Breast… ieu-a-… download… ieu-a-1126.… cis_anchor  chr20…     500000         6622          2551 0        
# 7 EGFR    Breast… ieu-a-… 2677_1_E… ieu-a-1126.… cis_anchor  chr7:…     500000         9110          4465 1.43e- 30
# 8 MMP10   Endome… GCST00… 8479_4_M… 30093612-GC… cis_anchor  chr11…     500000         8165          3433 1.07e-247
# 9 GSTO1   Endome… GCST00… 12436_84… 30093612-GC… cis_anchor  chr10…     500000         7055          2394 0        
# 10 GREM1   Colore… GCST90… 18878_15… GCST9025567… cis_anchor  chr15…     500000         8029          3050 0        
# 11 CCM2    Colore… GCST90… 12347_29… GCST9025567… cis_anchor  chr7:…     500000         7816          3622 9.90e- 12
# 12 HHIP    Colore… GCST90… 10833_64… GCST9025567… cis_anchor  chr4:…     500000         6324          2360 6.74e- 73
# 13 PTPN9   Colore… GCST90… 12633_3_… GCST9025567… cis_anchor  chr15…     500000         6754          2310 1.13e- 12
# 14 ACE     Colore… GCST90… 10714_7_… GCST9025567… cis_anchor  chr17…     500000         6477          2331 0        
# 15 LCT     Colore… GCST90… 9017_58_… GCST9025567… cis_anchor  chr2:…     500000         5696          2165 0    

df_coloc_all %>% filter(coloc_interpretation == "Moderate colocalization")
# 1 PEAR1   CVD     ebi-a-… 8275_31_… ebi-a-GCST9… cis_anchor  chr1:…     500000         8180          4215 4.95e-104
# 2 SPP1    CVD     ebi-a-… 13113_7_… ebi-a-GCST9… cis_anchor  chr4:…     500000         8086          4209 7.06e-100
# 3 NADK    CVD     ebi-a-… 13624_17… ebi-a-GCST9… cis_anchor  chr1:…     500000         9930          3935 1.05e- 47
# 4 ACADVL  CVD     ebi-a-… 19262_21… ebi-a-GCST9… cis_anchor  chr17…     500000         8112          3964 7.81e-  6
# 5 PARVA   T2D     ebi-a-… 13434_17… ebi-a-GCST0… cis_anchor  chr11…     500000         8516          1947 8.29e- 19
# 6 CTSS    All ca… ieu-b-… 3181_50_… ieu-b-4966.… cis_anchor  chr1:…     500000         6936          2829 0        
# 7 DNPH1   Breast… ieu-a-… 17785_11… ieu-a-1126.… cis_anchor  chr6:…     500000         6958          2722 3.96e- 26
# 8 ASIP    Breast… ieu-a-… 5676_54_… ieu-a-1126.… cis_anchor  chr20…     500000         6020          1889 0        
# 9 PTK7    Prosta… GCST90… 9525_1_P… GCST9027471… cis_anchor  chr6:…     500000         6931          3960 4.36e- 20
# 10 SERPIN… Prosta… GCST90… 2879_9_S… GCST9027471… cis_anchor  chr14…     500000         9468          6143 1.58e-126
# 11 CA11    Colore… GCST90… 8859_51_… GCST9025567… cis_anchor  chr19…     500000         8784          3941 5.39e- 60
# 12 POSTN   Colore… GCST90… 3457_57_… GCST9025567… cis_anchor  chr13…     500000         8795          4427 3.38e-105

##' all significant MR reslts
##' merge with it for forest plot
##' --------------------------------
mr_file <- file.path(MANUSCRIPT_DIR, "significant_MR_results_hPDI_uPDI_with_learning_indicators.xlsx")
# First check exact sheet names
excel_sheets(mr_file)
df_mr_sig <- read_excel(
  mr_file,
  sheet = "significant_MR_with_learning"
)

make_disease_key <- function(x) {
  
  x_clean <- x %>%
    stringr::str_to_lower() %>%
    stringr::str_squish()
  
  dplyr::case_when(
    x_clean %in% c("cancer", "all cancer") ~ "all_cancer",
    
    x_clean %in% c("breast", "breast cancer") ~ "breast",
    x_clean %in% c("colorectal", "colorectal cancer") ~ "colorectal",
    x_clean %in% c("prostate", "prostate cancer") ~ "prostate",
    x_clean %in% c("endometrial", "endometrial cancer") ~ "endometrial",
    x_clean %in% c("liver", "liver cancer") ~ "liver",
    x_clean %in% c("lung", "lung cancer") ~ "lung",
    
    x_clean %in% c("t2d", "type 2 diabetes", "type 2 diabetes mellitus") ~ "t2d",
    x_clean %in% c("cvd", "cardiovascular disease") ~ "cvd",
    
    TRUE ~ x_clean
  )
}

df_mr_sig2 <- df_mr_sig %>%
  dplyr::mutate(
    disease_key = make_disease_key(disease),
    protein_key = stringr::str_to_upper(stringr::str_squish(protein))
  )

df_coloc2 <- df_coloc_all %>%
  dplyr::mutate(
    disease_key = make_disease_key(disease),
    protein_key = stringr::str_to_upper(stringr::str_squish(protein))
  )

df_coloc2 %>%
  count(disease_key, protein_key) %>%
  filter(n > 1)

# disease_key protein_key     n
# <chr>       <chr>       <int>
#   1 breast      RSPO3           2

df_coloc2 <- df_coloc2 %>%
  dplyr::mutate(
    PP_H4_final = dplyr::coalesce(PP.H4, PP.H4.abf),
    PP_H3_final = dplyr::coalesce(PP.H3, PP.H3.abf),
    PP_H2_final = dplyr::coalesce(PP.H2, PP.H2.abf),
    PP_H1_final = dplyr::coalesce(PP.H1, PP.H1.abf),
    PP_H0_final = dplyr::coalesce(PP.H0, PP.H0.abf)
  )

df_coloc2 %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_PP_H4 = sum(!is.na(PP_H4_final)),
    n_strong = sum(PP_H4_final >= 0.80, na.rm = TRUE),
    n_moderate = sum(PP_H4_final >= 0.50 & PP_H4_final < 0.80, na.rm = TRUE),
    n_weak = sum(PP_H4_final < 0.50, na.rm = TRUE),
    n_missing = sum(is.na(PP_H4_final))
  )

# n_rows n_PP_H4 n_strong n_moderate n_weak n_missing
# <int>   <int>    <int>      <int>  <int>     <int>
#   1    110     110       18         14     78         0

df_coloc2_unique <- df_coloc2 %>%
  # Historical selection rule retained: highest PP.H4 per disease/protein.
  # With one result per disease/protein, this leaves the sole result unchanged.
  dplyr::arrange(
    disease_key,
    protein_key,
    dplyr::desc(PP_H4_final),
    dplyr::desc(nsnps_overlap)
  ) %>%
  dplyr::group_by(disease_key, protein_key) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

df_coloc2_unique %>%
  dplyr::count(disease_key, protein_key) %>%
  dplyr::filter(n > 1)

## MERGE
df_mr_coloc <- df_mr_sig2 %>%
  dplyr::left_join(
    df_coloc2_unique %>%
      dplyr::select(
        disease_key,
        protein_key,
        gwas_id_coloc = gwas_id,
        region,
        nsnps_overlap,
        PP_H0_final,
        PP_H1_final,
        PP_H2_final,
        PP_H3_final,
        PP_H4_final,
        dplyr::any_of(c("source_file", "pqtl_file", "outcome_file")),
        coloc_interpretation
      ),
    by = c("disease_key", "protein_key")
  ) %>%
  dplyr::mutate(
    coloc_category = dplyr::case_when(
      is.na(PP_H4_final)      ~ "No coloc result",
      PP_H4_final >= 0.80     ~ "Strong colocalization",
      PP_H4_final >= 0.50     ~ "Moderate colocalization",
      TRUE                    ~ "Weak/no colocalization"
    )
  )

df_mr_coloc %>%
  dplyr::count(coloc_category)
# coloc_category              n
# <chr>                   <int>
#   1 Moderate colocalization    13
# 2 Strong colocalization      20
# 3 Weak/no colocalization     79

if (interactive()) View(df_mr_coloc)

library(readr)

write_csv(
  df_mr_coloc,
  file.path(MANUSCRIPT_DIR, "MR_all_source_sigresults_with_coloc.csv")
)

##read file
mr <- read.csv(file.path(MANUSCRIPT_DIR, "MR_all_source_sigresults_with_coloc.csv"))

# distinct protein-disease pairs in MR results
assoc_all <- mr %>%
  distinct(exposure, disease)
nrow(assoc_all) #99 distinct exposure disease combinations

#distinct disaease - exposure by diet
sum_assoc_diet_disease <- mr %>%
  distinct(pdi_type, exposure, disease) %>%
  count(pdi_type, disease, name = "n_distinct_protein_disease")

sum_assoc_diet_disease

# pdi_type     disease n_distinct_protein_disease
# 1      hPDI      Breast                          5
# 2      hPDI         CVD                         17
# 3      hPDI      Cancer                          4
# 4      hPDI  Colorectal                          8
# 5      hPDI Endometrial                          2
# 6      hPDI       Liver                          1
# 7      hPDI        Lung                          1
# 8      hPDI    Prostate                         23
# 9      hPDI         T2D                          6
# 10     uPDI      Breast                          7
# 11     uPDI         CVD                         18
# 12     uPDI      Cancer                          1
# 13     uPDI  Colorectal                         12
# 14     uPDI Endometrial                          1
# 15     uPDI    Prostate                          6

#number that are colocalised
sum_coloc_diet <- mr %>%
  distinct(pdi_type , exposure, disease, coloc_interpretation) %>%
  group_by(pdi_type) %>%
  summarise(
    n_distinct_protein_disease = n_distinct(paste(exposure, disease)),
    n_colocalised = n_distinct(
      paste(exposure[coloc_interpretation == "Strong colocalization"],
            disease[coloc_interpretation == "Strong colocalization"])
    ),
    .groups = "drop"
  )

sum_coloc_diet
# pdi_type n_distinct_protein_disease n_colocalised
# <chr>                         <int>         <int>
#   1 hPDI                             67             9
# 2 uPDI                             45            12

sum_coloc_diet <- mr %>%
  distinct(pdi_type, exposure, disease, coloc_interpretation) %>%
  group_by(pdi_type) %>%
  summarise(
    n_distinct_protein_disease = n_distinct(paste(exposure, disease)),
    n_colocalised = n_distinct(
      paste(exposure[coloc_interpretation == "Strong colocalization"],
            disease[coloc_interpretation == "Strong colocalization"])
    ),
    .groups = "drop"
  )

sum_coloc_diet
# diet_score n_distinct_protein_disease n_colocalised
# <chr>                           <int>         <int>
#   1 hPDI                               57             7
# 2 uPDI                               33             4

colocalised_proteins <- mr %>%
  filter(coloc_interpretation == "Strong colocalization") %>%
  distinct(pdi_type, disease, exposure) %>%
  arrange(pdi_type, disease, exposure)

colocalised_proteins
# pdi_type     disease exposure
# 1      hPDI         CVD    PCSK9
# 2      hPDI  Colorectal     CCM2
# 3      hPDI  Colorectal    GREM1
# 4      hPDI Endometrial    GSTO1
# 5      hPDI Endometrial    MMP10
# 6      hPDI         T2D    LRIG1
# 7      hPDI         T2D    MANBA
# 8      hPDI         T2D     NCAN
# 9      uPDI      Breast    CPNE1
# 10     uPDI      Breast     EGFR
# 11     uPDI         CVD    PCSK9
# 12     uPDI      Cancer     KRT5
# 13     uPDI  Colorectal      ACE
# 14     uPDI  Colorectal     CCM2
# 15     uPDI  Colorectal    GREM1
# 16     uPDI  Colorectal     HHIP
# 17     uPDI  Colorectal      LCT
# 18     uPDI  Colorectal    PTPN9
# 19     uPDI Endometrial    MMP10

colocalised_proteins <- mr %>%
  filter(coloc_interpretation == "Moderate colocalization") %>%
  distinct(pdi_type, disease, exposure) %>%
  arrange(pdi_type, disease, exposure)

colocalised_proteins
# pdi_type    disease exposure
# 1      hPDI     Breast    DNPH1
# 2      hPDI        CVD    PEAR1
# 3      hPDI        CVD     SPP1
# 4      hPDI     Cancer     CTSS
# 5      hPDI Colorectal     CA11
# 6      hPDI Colorectal    POSTN
# 7      hPDI   Prostate     PTK7
# 8      hPDI        T2D    PARVA
# 9      uPDI        CVD   ACADVL
# 10     uPDI        CVD     NADK
# 11     uPDI        CVD     SPP1
# 12     uPDI Colorectal    POSTN
# 13     uPDI   Prostate SERPINA3

#plot
#___
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(stringr)
library(tidytext)

#-----------------------------
# 1. Prepare plotting data
#-----------------------------
disease_levels <- c(
  "All cancer", "Breast cancer", "Prostate cancer", "Colorectal cancer",
  "Endometrial cancer", "T2D", "CVD"
)

df_plot <- mr %>%
  filter(!is.na(PP_H4_final), PP_H4_final >= 0.50) %>%   # moderate and above
  mutate(
    disease_plot = case_when(
      disease %in% c("Cancer", "All cancer") ~ "All cancer",
      disease %in% c("Breast", "Breast cancer") ~ "Breast cancer",
      disease %in% c("Prostate", "Prostate cancer") ~ "Prostate cancer",
      disease %in% c("Colorectal", "Colorectal cancer") ~ "Colorectal cancer",
      disease %in% c("Endometrial", "Endometrial cancer") ~ "Endometrial cancer",
      disease %in% c("T2D", "Type 2 diabetes") ~ "T2D",
      disease %in% c("CVD") ~ "CVD",
      TRUE ~ disease
    ),
    disease_plot = factor(disease_plot, levels = disease_levels),
    pdi_type = factor(pdi_type, levels = c("hPDI", "uPDI"))
  ) %>%
  arrange(disease_plot, pdi_type, desc(PP_H4_final)) %>%
  mutate(
    # order proteins within each disease x pdi facet
    protein_facet = tidytext::reorder_within(
      protein,
      PP_H4_final,
      interaction(disease_plot, pdi_type)
    )
  )

#-----------------------------
# 2. Source indicator tiles
#-----------------------------
# Put the source tiles on the left of the OR axis
source_positions <- c(XGB = -0.60, ENET = -0.42, limma = -0.24)

df_source <- df_plot %>%
  select(disease_plot, pdi_type, protein, protein_facet, gbm, enet, limma) %>%
  pivot_longer(
    cols = c(gbm, enet, limma),
    names_to = "source",
    values_to = "present"
  ) %>%
  mutate(
    source = recode(
      source,
      "gbm" = "XGB",
      "enet" = "ENET",
      "limma" = "limma"
    ),
    source = factor(source, levels = c("XGB", "ENET", "limma")),
    x = unname(source_positions[as.character(source)])
  )

# PP.H4 label position on the right
x_pph4 <- max(df_plot$or_uci95, na.rm = TRUE) + 0.18

# x-axis breaks
or_breaks <- pretty(c(df_plot$or_lci95, df_plot$or_uci95), n = 4)
or_breaks <- or_breaks[or_breaks > 0]

x_breaks <- c(source_positions["XGB"], source_positions["ENET"], source_positions["limma"], or_breaks)
x_labels <- c("XGB", "ENET", "limma", format(or_breaks, trim = TRUE))

#-----------------------------
# 3. Plot
#-----------------------------
p <- ggplot() +
  # left-hand source tiles
  geom_tile(
    data = df_source,
    aes(x = x, y = protein_facet, fill = factor(present)),
    width = 0.12,
    height = 0.75
  ) +
  scale_fill_manual(
    values = c("0" = "grey75", "1" = "grey25"),
    guide = "none"
  ) +
  
  # divider between source tiles and forest plot
  geom_vline(xintercept = -0.08, colour = "grey60", linewidth = 0.35) +
  
  # null line at OR = 1
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  
  # confidence intervals
  geom_segment(
    data = df_plot,
    aes(
      x = or_lci95,
      xend = or_uci95,
      y = protein_facet,
      yend = protein_facet
    ),
    linewidth = 0.45,
    colour = "grey35"
  ) +
  
  # OR points colored by colocalization strength
  geom_point(
    data = df_plot,
    aes(x = or, y = protein_facet, colour = PP_H4_final),
    size = 2.4
  ) +
  
  # PP.H4 text on the right
  geom_text(
    data = df_plot,
    aes(
      x = x_pph4,
      y = protein_facet,
      label = sprintf("%.2f", PP_H4_final)
    ),
    hjust = 0,
    size = 2.7,
    colour = "grey20"
  ) +
  
  facet_grid(
    disease_plot ~ pdi_type,
    scales = "free_y",
    space = "free_y",
    drop = TRUE,
    switch = "y"
  ) +
  
  tidytext::scale_y_reordered() +
  
  scale_x_continuous(
    breaks = x_breaks,
    labels = x_labels,
    expand = expansion(mult = c(0.02, 0.20))
  ) +
  
  scale_colour_gradient(
    low = "#2B3256",
    high = "#D27E64",
    limits = c(0.50, 1.00),
    name = "PP.H4"
  ) +
  
  labs(
    x = "Odds ratio per 1-SD higher genetically predicted plasma protein abundance (95% CI)",
    y = NULL,
    title = "MR-significant proteins with moderate or strong colocalization",
    subtitle = "Rows = disease, columns = diet type; point color and right-hand text indicate PP.H4"
  ) +
  
  coord_cartesian(clip = "off") +
  
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.margin = margin(5.5, 70, 5.5, 5.5),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold"),
    strip.text.y.right = element_blank()
  )

p

ggsave(
  filename = file.path(MANUSCRIPT_DIR, "forestplot_withcolocalisation.png"),
  plot = p,
  width = 21,
  height = 27.9,
  units = "cm",
  dpi = 300
)

## another version with CI
library(tidytext)
library(ggh4x)

#-----------------------------
# 1. Prepare plotting data
#-----------------------------
disease_levels <- c(
  "All cancer", "Breast cancer", "Prostate cancer", "Colorectal cancer",
  "Endometrial cancer", "T2D", "CVD"
)

pdi_levels <- c("hPDI", "uPDI")

# Start from filtered data
df_base <- mr %>%
  filter(!is.na(PP_H4_final), PP_H4_final >= 0.50) %>%
  mutate(
    disease_plot = case_when(
      disease %in% c("Cancer", "All cancer") ~ "All cancer",
      disease %in% c("Breast", "Breast cancer") ~ "Breast cancer",
      disease %in% c("Prostate", "Prostate cancer") ~ "Prostate cancer",
      disease %in% c("Colorectal", "Colorectal cancer") ~ "Colorectal cancer",
      disease %in% c("Endometrial", "Endometrial cancer") ~ "Endometrial cancer",
      disease %in% c("T2D", "Type 2 diabetes") ~ "T2D",
      disease %in% c("CVD") ~ "CVD",
      TRUE ~ disease
    ),
    disease_plot = factor(disease_plot, levels = disease_levels),
    pdi_type = factor(pdi_type, levels = pdi_levels)
  )

# Create the union of proteins within each disease
# basically create the full list of proteins within each disease
# only has unique combination
protein_union <- df_base %>%
  distinct(disease_plot, protein)

# Expand each disease-level protein list to both hPDI and uPDI
# expands every disease-protein combination into two rows: one for hPDI, one for uPDI
plot_scaffold <- protein_union %>%
  tidyr::crossing(
    pdi_type = factor(pdi_levels, levels = pdi_levels)
  )

# Join actual data back onto the scaffold
df_plot <- plot_scaffold %>%
  left_join(
    df_base,
    by = c("disease_plot", "pdi_type", "protein")
  )

# Determine ordering of proteins within each disease
# (here by strongest PP_H4 across hPDI/uPDI)
protein_order <- df_base %>%
  group_by(disease_plot, protein) %>%
  summarise(
    order_value = max(PP_H4_final, na.rm = TRUE),
    .groups = "drop"
  )

df_plot <- df_plot %>%
  left_join(protein_order, by = c("disease_plot", "protein")) %>%
  mutate(
    protein_facet = tidytext::reorder_within(
      protein,
      order_value,
      disease_plot
    ),
    right_label = ifelse(
      !is.na(or),
      sprintf("%.2f (%.2f, %.2f) | %.2f", or, or_lci95, or_uci95, PP_H4_final),
      ""
    )
  )

#-----------------------------
# 2. Source indicator tiles
#-----------------------------
source_positions <- c(XGB = -0.60, ENET = -0.42, limma = -0.24)

df_source <- df_plot %>%
  select(disease_plot, pdi_type, protein, protein_facet, gbm, enet, limma) %>%
  pivot_longer(
    cols = c(gbm, enet, limma),
    names_to = "source",
    values_to = "present"
  ) %>%
  mutate(
    source = recode(
      source,
      gbm = "XGB",
      enet = "ENET",
      limma = "limma"
    ),
    source = factor(source, levels = c("XGB", "ENET", "limma")),
    x = unname(source_positions[as.character(source)])
  )

#-----------------------------
# 3. Axis positions / labels
#-----------------------------
x_right <- max(df_plot$or_uci95, na.rm = TRUE) +0.15

or_breaks <- pretty(c(df_plot$or_lci95, df_plot$or_uci95), n = 4)
or_breaks <- or_breaks[or_breaks > 0]

x_breaks <- c(source_positions["XGB"], source_positions["ENET"], source_positions["limma"], or_breaks)
x_labels <- c("XGB", "ENET", "limma", format(or_breaks, trim = TRUE))

#-----------------------------
# 3a. Add  point-shape groups
#-----------------------------
highlight_proteins <- c("PCSK9", "MMP10", "EGFR")

df_plot <- df_plot %>%
  mutate(
    n_approaches = rowSums(
      across(c(gbm, enet, limma), ~ as.numeric(.x == 1)),
      na.rm = TRUE
    ),
    point_group = case_when(
      PP_H4_final >= 0.80 & n_approaches >= 2 ~
        "≥2 approaches + strong colocalization",
      PP_H4_final >= 0.80 & n_approaches < 2 ~
        "Other strong colocalization",
      PP_H4_final >= 0.50 & PP_H4_final < 0.80 ~
        "Moderate colocalization",
      TRUE ~ NA_character_
    ),
    point_group = factor(
      point_group,
      levels = c(
        "≥2 approaches + strong colocalization",
        "Other strong colocalization",
        "Moderate colocalization"
      )
    )
  )

#-----------------------------
# 4. Plot
#-----------------------------
p2 <-ggplot() +
  # left-hand source tiles
  geom_tile(
    data = df_source,
    aes(x = x, y = protein_facet, fill = factor(present)),
    width = 0.12,
    height = 0.75
  ) +
  scale_fill_manual(
    values = c("0" = "grey75", "1" = "grey25"),
    na.value = "white",
    guide = "none"
  ) +
  
  # divider between source tiles and forest plot
  geom_vline(xintercept = -0.08, colour = "grey60", linewidth = 0.35) +
  
  # null line at OR = 1
  geom_vline(xintercept = 1, linetype = "dashed", colour = "darkred", linewidth = 0.3) +
  
  # confidence intervals
  geom_segment(
    data = df_plot %>% filter(!is.na(or_lci95), !is.na(or_uci95)),
    aes(
      x = or_lci95,
      xend = or_uci95,
      y = protein_facet,
      yend = protein_facet
    ),
    linewidth = 0.45,
    colour = "grey35"
  ) +
  
  # OR points
  geom_point(
    data = df_plot %>% filter(!is.na(or), !is.na(PP_H4_final)),
    aes(
      x = or,
      y = protein_facet,
      colour = PP_H4_final,
      shape = point_group
    ),
    size = 1.8,
    stroke = 0.7
  ) +
  
  # right-hand label
  geom_text(
    data = df_plot %>% filter(!is.na(or)),
    aes(
      x = x_right,
      y = protein_facet,
      label = right_label
    ),
    hjust = 0,
    size = 2.7,
    colour = "grey20",
    family = "Helvetica"
  ) +
  
  facet_grid(
    disease_plot ~ pdi_type,
    scales = "free_y",
    space = "free_y",
    drop = TRUE,
    switch = "y"
  ) +
  
  tidytext::scale_y_reordered() +
  
  scale_x_continuous(
    breaks = x_breaks,
    labels = x_labels,
    expand = expansion(mult = c(0.02, 0.40))
  ) +
  
  scale_colour_gradient(
    low = "#2B3256",
    high = "#D27E64",
    limits = c(0.50, 1.00),
    name = "PP.H4"
  ) +
  
  scale_shape_manual(
    values = c(
      "≥2 approaches + strong colocalization" = 17,  # triangle
      "Other strong colocalization" = 16,            # circle
      "Moderate colocalization" = 124                 # triangle
    ),
    name = "Evidence tier"
  ) +
  
  labs(
    x = "Odds ratio per 1-SD higher genetically predicted plasma protein abundance",
    y = NULL,
    title = "MR-significant proteins with moderate or strong colocalization",
    subtitle = "Right-hand label shows OR (95% CI) | PP.H4"
  ) +
  
  coord_cartesian(clip = "off") +
  
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.margin = margin(5.5, 120, 5.5, 5.5),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold"),
    strip.text.y.right = element_blank()
  )

p2

ggsave(
  filename = file.path(MANUSCRIPT_DIR, "forestplot_withcolocalisation_CI_evidence.png"),
  plot = p2,
  width = 29.7,
  height = 21,
  units = "cm",
  dpi = 300
)

## revised - plot
##================
#============================================================
# MR + colocalization plot
# Four prioritization sections
#============================================================

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tidytext)
library(patchwork)


#------------------------------------------------------------
# 1. Read data
#------------------------------------------------------------

mr <- read.csv(
  file.path(MANUSCRIPT_DIR, "MR_all_source_sigresults_with_coloc.csv")
)


#------------------------------------------------------------
# 2. Plotting orders
#------------------------------------------------------------

#-----------------------------
# 1. Display definitions
#-----------------------------
disease_levels <- c(
  "All cancer", "Breast cancer", "Prostate cancer", "Colorectal cancer",
  "Endometrial cancer", "T2D", "CVD"
)

pdi_levels <- c("hPDI", "uPDI")

# >>> CHANGE HERE 1: edit this vector to change the protein display order.
priority_proteins <- c("EGFR", "PCSK9", "MMP10", "PTK7", "CCM2", "NCAN")

# >>> CHANGE HERE 2: edit the PP.H4 thresholds if your definitions change.
moderate_pph4_cutoff <- 0.50
strong_pph4_cutoff <- 0.80

# >>> CHANGE HERE 3: edit these labels to rename or reorder the four panels.
section_levels <- c(
  "Identified ≥2 discovery approaches + Strong colocalization",
  "Identified by LIMMA + Strong colocalization",
  "Identified by XGB + Strong colocalization",
  "Moderate colocalization"
)

#-----------------------------
# 2. Prepare plotting data
#-----------------------------
df_base <- mr %>%
  filter(
    !is.na(PP_H4_final),
    PP_H4_final >= moderate_pph4_cutoff
  ) %>%
  mutate(
    disease_plot = case_when(
      disease %in% c("Cancer", "All cancer") ~ "All cancer",
      disease %in% c("Breast", "Breast cancer") ~ "Breast cancer",
      disease %in% c("Prostate", "Prostate cancer") ~ "Prostate cancer",
      disease %in% c("Colorectal", "Colorectal cancer") ~ "Colorectal cancer",
      disease %in% c("Endometrial", "Endometrial cancer") ~ "Endometrial cancer",
      disease %in% c("T2D", "Type 2 diabetes") ~ "T2D",
      disease == "CVD" ~ "CVD",
      TRUE ~ disease
    ),
    disease_plot = factor(disease_plot, levels = disease_levels),
    pdi_type = factor(pdi_type, levels = pdi_levels),
    protein_key = str_to_upper(protein),
    requested_rank = match(protein_key, priority_proteins),
    gbm_flag = as.integer(
      !is.na(gbm) &
        as.character(gbm) %in% c("1", "TRUE", "True", "true")
    ),
    enet_flag = as.integer(
      !is.na(enet) &
        as.character(enet) %in% c("1", "TRUE", "True", "true")
    ),
    limma_flag = as.integer(
      !is.na(limma) &
        as.character(limma) %in% c("1", "TRUE", "True", "true")
    ),
    # >>> CHANGE HERE 4: this counts discovery approaches for each
    # protein–diet-index result. Keep all three terms if all approaches
    # should contribute equally to the count.
    n_approaches = gbm_flag + enet_flag + limma_flag
  )

# Assign each protein–disease combination to exactly one of the four panels.
# The hierarchy is deliberate: multi-approach strong results first, followed
# by LIMMA strong, XGB strong, and finally moderate colocalization.
section_lookup <- df_base %>%
  group_by(protein_key, disease_plot) %>%
  summarise(
    multi_strong = any(
      PP_H4_final >= strong_pph4_cutoff & n_approaches >= 2L,
      na.rm = TRUE
    ),
    limma_strong = any(
      PP_H4_final >= strong_pph4_cutoff & limma_flag == 1L,
      na.rm = TRUE
    ),
    enet_only_strong = any(
      PP_H4_final >= strong_pph4_cutoff &
        enet_flag == 1L & limma_flag == 0L & gbm_flag == 0L,
      na.rm = TRUE
    ),
    xgb_strong = any(
      PP_H4_final >= strong_pph4_cutoff & gbm_flag == 1L,
      na.rm = TRUE
    ),
    max_pph4 = max(PP_H4_final, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # >>> CHANGE HERE 5: these are the four classification rules.
    # If ENET-only strong findings should join the LIMMA panel, change
    # `limma_strong` below to `limma_strong | enet_only_strong`.
    section = case_when(
      multi_strong ~
        "Identified ≥2 discovery approaches + Strong colocalization",
      limma_strong ~
        "Identified by LIMMA + Strong colocalization",
      enet_only_strong ~
        NA_character_,
      xgb_strong ~
        "Identified by XGB + Strong colocalization",
      max_pph4 >= moderate_pph4_cutoff &
        max_pph4 < strong_pph4_cutoff ~
        "Moderate colocalization",
      TRUE ~ NA_character_
    ),
    section = factor(section, levels = section_levels)
  )

unclassified_sections <- section_lookup %>%
  filter(is.na(section))

if (nrow(unclassified_sections) > 0) {
  stop(
    "At least one strong-colocalization protein–disease result could not be ",
    "assigned to the multi-approach, LIMMA, or XGB panels. This may be an ",
    "ENET-only result. See CHANGE HERE 5."
  )
}

df_base <- df_base %>%
  left_join(
    section_lookup %>% select(protein_key, disease_plot, section),
    by = c("protein_key", "disease_plot")
  )

# >>> CHANGE HERE 6: proteins not listed in `priority_proteins` are ordered
# by their strongest PP.H4. Replace `desc(max_pph4)` if another ordering is
# preferred.
other_order <- df_base %>%
  filter(is.na(requested_rank)) %>%
  group_by(protein_key) %>%
  summarise(max_pph4 = max(PP_H4_final, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_pph4), protein_key) %>%
  mutate(other_rank = length(priority_proteins) + row_number()) %>%
  select(protein_key, other_rank)

df_base <- df_base %>%
  left_join(other_order, by = "protein_key") %>%
  mutate(protein_rank = coalesce(requested_rank, other_rank))

missing_priority <- setdiff(priority_proteins, unique(df_base$protein_key))
if (length(missing_priority) > 0) {
  message(
    "These requested proteins are absent after the PP.H4 filter: ",
    paste(missing_priority, collapse = ", ")
  )
}

# Order the nested disease panels according to the first requested protein
# appearing in each panel. This preserves the requested overall hierarchy even
# when the proteins relate to different disease outcomes.
panel_order <- df_base %>%
  distinct(section, disease_plot, protein_rank) %>%
  group_by(section, disease_plot) %>%
  summarise(panel_rank = min(protein_rank, na.rm = TRUE), .groups = "drop") %>%
  arrange(section, panel_rank, disease_plot) %>%
  mutate(panel_key = paste(section, disease_plot, sep = "___"))

panel_levels <- panel_order$panel_key
panel_labels <- setNames(
  as.character(panel_order$disease_plot),
  panel_order$panel_key
)

df_base <- df_base %>%
  mutate(
    panel_key = factor(
      paste(section, disease_plot, sep = "___"),
      levels = panel_levels
    )
  )

# Use the union of proteins within each nested panel so that hPDI and uPDI
# columns have aligned rows.
plot_scaffold <- df_base %>%
  distinct(
    section, panel_key, disease_plot,
    protein, protein_key, protein_rank
  ) %>%
  crossing(
    pdi_type = factor(pdi_levels, levels = pdi_levels)
  )

df_plot <- plot_scaffold %>%
  left_join(
    df_base,
    by = c(
      "section", "panel_key", "disease_plot",
      "protein", "protein_key", "protein_rank", "pdi_type"
    )
  ) %>%
  mutate(
    # Negative rank makes rank 1 appear at the top of a discrete y-axis.
    protein_facet = tidytext::reorder_within(
      protein,
      -protein_rank,
      panel_key
    ),
    right_label = if_else(
      !is.na(or),
      sprintf(
        "%.2f (%.2f, %.2f) | %.2f",
        or, or_lci95, or_uci95, PP_H4_final
      ),
      ""
    )
  )

#-----------------------------
# 3. Discovery-method indicators
#-----------------------------
source_positions <- c(XGB = -0.60, ENET = -0.42, LIMMA = -0.24)

df_source <- df_plot %>%
  select(
    section, panel_key, disease_plot, pdi_type,
    protein, protein_facet, gbm_flag, enet_flag, limma_flag
  ) %>%
  pivot_longer(
    cols = c(gbm_flag, enet_flag, limma_flag),
    names_to = "source",
    values_to = "present"
  ) %>%
  mutate(
    source = recode(
      source,
      gbm_flag = "XGB",
      enet_flag = "ENET",
      limma_flag = "LIMMA"
    ),
    source = factor(source, levels = c("XGB", "ENET", "LIMMA")),
    x = unname(source_positions[as.character(source)])
  )

#-----------------------------
# 4. Axis positions
#-----------------------------
ci_values <- c(df_plot$or_lci95, df_plot$or_uci95)
ci_values <- ci_values[is.finite(ci_values)]

if (length(ci_values) == 0) {
  stop("No finite confidence intervals are available for plotting.")
}

ci_range <- range(ci_values)
ci_span <- diff(ci_range)
if (!is.finite(ci_span) || ci_span == 0) {
  ci_span <- 0.25
}

x_right <- ci_range[2] + 0.10 * ci_span
or_breaks <- pretty(ci_range, n = 4)
or_breaks <- or_breaks[or_breaks > 0]
x_breaks <- c(unname(source_positions), or_breaks)
x_labels <- c(names(source_positions), format(or_breaks, trim = TRUE))

#-----------------------------
# 5. Plot each prioritization section as a separate panel
#-----------------------------
make_section_plot <- function(section_name, show_x_axis = FALSE) {
  section_data <- df_plot %>%
    filter(as.character(section) == section_name)
  
  section_sources <- df_source %>%
    filter(as.character(section) == section_name)
  
  ggplot() +
    geom_tile(
      data = section_sources,
      aes(x = x, y = protein_facet, fill = factor(present)),
      width = 0.12,
      height = 0.75
    ) +
    scale_fill_manual(
      values = c("0" = "grey78", "1" = "grey25"),
      na.value = "white",
      guide = "none"
    ) +
    geom_vline(
      xintercept = -0.08,
      colour = "grey60",
      linewidth = 0.35
    ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      colour = "darkred",
      linewidth = 0.3
    ) +
    geom_segment(
      data = section_data %>%
        filter(!is.na(or_lci95), !is.na(or_uci95)),
      aes(
        x = or_lci95,
        xend = or_uci95,
        y = protein_facet,
        yend = protein_facet
      ),
      linewidth = 0.45,
      colour = "grey35"
    ) +
    geom_point(
      data = section_data %>%
        filter(!is.na(or), !is.na(PP_H4_final)),
      aes(
        x = or,
        y = protein_facet,
        colour = PP_H4_final
      ),
      shape = 16,
      size = 2.1
    ) +
    geom_text(
      data = section_data %>% filter(!is.na(or)),
      aes(
        x = x_right,
        y = protein_facet,
        label = right_label
      ),
      hjust = 0,
      size = 2.7,
      colour = "grey20",
      family = "Helvetica"
    ) +
    facet_grid(
      rows = vars(panel_key),
      cols = vars(pdi_type),
      scales = "free_y",
      space = "free_y",
      drop = TRUE,
      switch = "y",
      labeller = labeller(
        panel_key = as_labeller(panel_labels)
      )
    ) +
    tidytext::scale_y_reordered() +
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      expand = expansion(mult = c(0.02, 0.45))
    ) +
    scale_colour_gradient(
      low = "#2B3256",
      high = "#D27E64",
      limits = c(0.50, 1.00),
      name = "PP.H4"
    ) +
    labs(
      x = if (show_x_axis) {
        "Odds ratio per 1-SD higher genetically predicted plasma protein abundance"
      } else {
        NULL
      },
      y = NULL,
      title = section_name
    ) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 8),
      axis.text.x = if (show_x_axis) {
        element_text(size = 8, angle = 45, hjust = 1)
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x_axis) {
        element_line()
      } else {
        element_blank()
      },
      legend.position = "bottom",
      plot.margin = margin(2, 120, 2, 5.5),
      plot.title = element_text(
        face = "bold",
        size = 10,
        margin = margin(3, 0, 4, 0)
      ),
      plot.title.position = "plot",
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, face = "bold"),
      strip.text.y.right = element_blank()
    )
}

available_sections <- section_levels[
  section_levels %in% as.character(unique(df_plot$section))
]

section_heights <- vapply(
  available_sections,
  function(section_name) {
    df_plot %>%
      filter(as.character(section) == section_name) %>%
      distinct(panel_key, protein) %>%
      nrow()
  },
  numeric(1)
)

# Add a small constant so that the title and disease strips have adequate room
# even in the compact first and second panels.
section_heights <- section_heights + 1.5

section_plots <- lapply(
  seq_along(available_sections),
  function(i) {
    make_section_plot(
      section_name = available_sections[i],
      show_x_axis = i == length(available_sections)
    )
  }
)

p <- patchwork::wrap_plots(
  section_plots,
  ncol = 1,
  heights = section_heights,
  guides = "collect"
) +
  patchwork::plot_annotation(
    title = "MR-significant proteins with moderate or strong colocalization",
    subtitle = paste(
      "Tiles show diet-specific discovery support;",
      "right-hand labels show OR (95% CI) | PP.H4"
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9)
    )
  ) &
  theme(legend.position = "bottom")

p

ggsave(
  filename = paste0(
    PROJECT_ROOT,
    "outputs/manuscript/forestplot_withcolocalisation_revised.png"
  ),
  plot = p,
  width = 21,
  height = 27.9,
  units = "cm",
  dpi = 300
)
