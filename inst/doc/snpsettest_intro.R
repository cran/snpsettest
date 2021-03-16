## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----GWAS summary-------------------------------------------------------------
library(snpsettest)

# Check an example of GWAS summary file (included in this package)
head(exGWAS, 3)

## ----reference data-----------------------------------------------------------
# Path to .bed file
bfile <- system.file("extdata", "example.bed", package = "snpsettest")

# Read a .bed file using bed.matrix-class in gaston package
# Genotypes are retrieved on demand to manage large-scale genotype data
x <- read_reference_bed(bfile, verbose = FALSE)

## ----harmonization------------------------------------------------------------
# Harmonize by SNP IDs
hsumstats1 <- harmonize_sumstats(exGWAS, x)

# Harmonize by genomic position and allele codes
# Reference allele swap will be taken into account (e.g., A/C match C/A)
hsumstats2 <- harmonize_sumstats(exGWAS, x, match_by_id = FALSE)

# Check matching entries by flipping allele codes (e.g., A/C match T/G)
# Ambiguous SNPs will be excluded from harmonization
hsumstats3 <- harmonize_sumstats(exGWAS, x, match_by_id = FALSE, check_strand_flip = TRUE)

## ----snp2gene-----------------------------------------------------------------
# Check gene information from the GENCODE project (included in this package)
head(gene.curated.GRCh37, 3)

# Map SNPs to genes
snp_sets <- map_snp_to_gene(hsumstats1, gene.curated.GRCh37)
str(snp_sets$sets[1:5])

# Allows a certain kb window before/after the gene to be included for SNP mapping
snp_sets_50kb <- map_snp_to_gene(
  hsumstats1, gene.curated.GRCh37, 
  extend_start = 50, extend_end = 50 # default is 20kb
)

## ----tests--------------------------------------------------------------------
# Perform gene-based association tests for the first 5 genes
res <- snpset_test(hsumstats1, x, snp_sets$sets[1:5])

# Show output
res

