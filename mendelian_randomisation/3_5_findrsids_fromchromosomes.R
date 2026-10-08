# Run from the repository root in a fresh R session. See README.md.
# Optional interactive rsID/liftOver helper with example loci and online queries; review before running.
# EDIT LOCAL PATHS HERE. No study data or results are supplied.
PROJECT_ROOT <- "."
DATA_DIR <- file.path(PROJECT_ROOT, "data")
OUTPUT_ROOT <- file.path(PROJECT_ROOT, "outputs")
MANUSCRIPT_DIR <- file.path(OUTPUT_ROOT, "manuscript")
REFERENCE_DIR <- file.path(DATA_DIR, "reference")
SIGNATURE_SCRIPTS <- file.path(PROJECT_ROOT, "signatures")
CONTEXT_SCRIPTS <- file.path(PROJECT_ROOT, "contextualisation")
dir.create(MANUSCRIPT_DIR, recursive = TRUE, showWarnings = FALSE)

library(httr)
library(data.table)

chrpos_instruments <- data.frame(
  chr = c("12", "12", "6", "8", "11", "4", "3"),
  pos = c(71495073, 69234593, 1719622, 23834609, 48748624, 77607718, 128642712)
)

get_rsid <- function(chr, pos) {
  url <- paste0("https://rest.ensembl.org/overlap/region/human/",
                chr, ":", pos, "-", pos,
                "?feature=variation;content-type=application/json")
  
  response <- GET(url)
  
  if (status_code(response) == 200) {
    hits <- content(response, "parsed")
    if (length(hits) > 0) {
      return(hits[[1]]$id)  # returns RSID
    }
  }
  return(NA)
}

chrpos_instruments$rsid <- mapply(get_rsid, 
                                  chrpos_instruments$chr, 
                                  chrpos_instruments$pos)

print(chrpos_instruments)

# chr      pos         rsid
# 1  12 71495073   rs11178833
# 2  12 69234593   rs72277077
# 3   6  1719622 rs1465912134
# 4   8 23834609   rs10573085
# 5  11 48748624 rs1565143357

#' from hpdi

chrpos_instruments <- data.frame(
  chr = c("1", "4", "1"),
  pos = c(103682950,98262079,158937849)
)

get_rsid <- function(chr, pos) {
  url <- paste0("https://rest.ensembl.org/overlap/region/human/",
                chr, ":", pos, "-", pos,
                "?feature=variation;content-type=application/json")
  
  response <- GET(url)
  
  if (status_code(response) == 200) {
    hits <- content(response, "parsed")
    if (length(hits) > 0) {
      return(hits[[1]]$id)  # returns RSID
    }
  }
  return(NA)
}

chrpos_instruments$rsid <- mapply(get_rsid, 
                                  chrpos_instruments$chr, 
                                  chrpos_instruments$pos)

print(chrpos_instruments)

# chr       pos        rsid
# 1   1 103682950        <NA>
#   2   4  98262079 rs748259104
# 3   1 158937849        <NA>

#'1:158907639 (GRCh37)
#' 1:158937849 is rs2526301806, PYHIN1 on 
#' https://www.ncbi.nlm.nih.gov/snp/?term=1%5BChromosome%5D+AND+158937849%5BBase+Position%5D
#' 
# Install rtracklayer and GenomicRanges separately before running.
library(rtracklayer)
library(GenomicRanges)

chain_url <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"

dir.create(REFERENCE_DIR, recursive = TRUE, showWarnings = FALSE)
chain_gz <- file.path(REFERENCE_DIR, "hg38ToHg19.over.chain.gz")
chain_file <- file.path(REFERENCE_DIR, "hg38ToHg19.over.chain")

# remove old/corrupted files if they exist
if (file.exists(chain_gz)) file.remove(chain_gz)
if (file.exists(chain_file)) file.remove(chain_file)

# download again as binary
download.file(chain_url, destfile = chain_gz, mode = "wb")
R.utils::gunzip(chain_gz, destname = chain_file, overwrite = TRUE, remove = FALSE)
chain <- rtracklayer::import.chain(chain_file)

gr38 <- GRanges(
  seqnames = "chr1",
  ranges = IRanges(start = 103682950, end = 103682950)  # ✅
)

gr37 <- liftOver(gr38, chain)
as.data.frame(unlist(gr37))
# seqnames     start       end width strand
# 1     chr1 104225572 104225572     1      *
