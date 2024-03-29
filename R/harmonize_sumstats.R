##' Harmonizing GWAS summary to reference data
##'
##' Finds an intersection of variants between GWAS summary and reference data.
##'
##' Pre-processing of GWAS summary data is required because the sets of variants
##' available in a particular GWAS might be poorly matched to the variants in
##' reference data. SNP matching can be performed either 1) by SNP ID or 2) by
##' chromosome code, base-pair position, and allele codes, while taking into
##' account possible strand flips and reference allele swap. For matched
##' entries, the SNP IDs in GWAS summary data are replaced with the ones in the
##' reference data.
##'
##' @param sumstats A data frame with two columns: "id" and "pvalue".
##' - id = SNP ID (e.g., rs numbers)
##' - pvalue = SNP-level p value
##'
##' If `match_by_id = FALSE`, it requires additional columns: "chr", "pos", "A1"
##' and "A2".
##' - chr =  chromosome
##' - pos =  base-pair position (must be integer)
##' - A1, A2 = allele codes (allele order is not important)
##'
##' @param x A `bed.matrix` object created using the reference data.
##' @param match_by_id If `TRUE`, SNP matching will be performed by SNP IDs
##'   instead of genomic position and allele codes. Default is `TRUE`.
##' @param check_strand_flip Only applies when `match_by_id = FALSE`. If `TRUE`,
##'   the function 1) removes ambiguous A/T and G/C SNPs for which the strand is
##'   not obvious, and 2) attempts to find additional matching entries by
##'   flipping allele codes (i.e., A->T, T->A, C->G, G->A). If the GWAS genotype
##'   data itself is used as the reference data, it would be safe to set
##'   `FALSE`. Default is `FALSE`.
##' @param return_indice Only applied when `match_by_id = FALSE`. If `TRUE`, the
##'   function provides an additional column indicating whether the match is
##'   with swapped alleles. If `check_strand_flip = TRUE`, the function also
##'   provides an additional column indicating whether the match is with flipped
##'   strand. Unnecessary for gene-based tests in this package, but may be
##'   useful for other purposes (e.g., harmonization for meta-analysis that
##'   needs to flip the sign of beta for a match with swapped alleles).
##' @return A data frame with columns: "id", "chr", "pos", "A1", "A2" and
##'   "pvalue". If `return_indice = TRUE`, the data frame includes additional
##'   columns `key_`, `swapped_`, and `flipped_`. `key_` is "chr_pos_A1_A2" in
##'   `sumstat` (the original input before harmonization). `swapped_` contains a
##'   logical vector indicating reference allele swap. `flipped_` contains a
##'   logical vector indicating strand flip.
##'
##' @examples
##' \dontshow{data.table::setDTthreads(1)}
##' ## GWAS summary statistics
##' head(exGWAS)
##'
##' ## Load reference genotype data
##' bfile <- system.file("extdata", "example.bed", package = "snpsettest")
##' x <- read_reference_bed(path = bfile)
##'
##' ## Harmonize by SNP IDs
##' hsumstats1 <- harmonize_sumstats(exGWAS, x)
##'
##' ## Harmonize by genomic position and allele codes
##' ## Reference allele swap will be taken into account
##' hsumstats2 <- harmonize_sumstats(exGWAS, x, match_by_id = FALSE)
##'
##' ## Check matching entries by flipping allele codes
##' ## Ambiguous SNPs will be excluded from harmonization
##' hsumstats3 <- harmonize_sumstats(exGWAS, x, match_by_id = FALSE,
##'                                  check_strand_flip = TRUE)
##' @export
harmonize_sumstats <- function(sumstats, x,
                               match_by_id = TRUE,
                               check_strand_flip = FALSE,
                               return_indice = FALSE
                               ) {

  ## Save the names of data input for error message
  sumstats_name <- deparse(substitute(sumstats))

  ## Assert function arguments
  is_df(sumstats)
  is_bed_matrix(x)
  is_tf(check_strand_flip, "check_strand_flip")
  is_tf(match_by_id, "match_by_id")

  ## Make data.table object
  if (!inherits(sumstats, "data.table")) {
    setDT(sumstats)
  }
  if (!inherits(x@snps, "data.table")) {
    setDT(x@snps)
  }

  message("-----")
  message("Checking the reference data for harmonization...")

  ## Check monomorphic SNPs in the reference data
  ref_monomorphic <- which(x@sigma == 0) # due to possible heterozyguous mono (maf = 0.05)
  message("Found ", length(ref_monomorphic),
          " monomoprhic SNPs in the reference data.")

  if (match_by_id) {

    ## Check duplicate SNP IDs in the reference data
    ## All duplicate-id variants are excluded from harmonization
    ref_dup <- get_duplicate_indice(x@snps$id)
    message("Found ", length(ref_dup),
            " duplicate SNP IDs in the reference data.")

    ## Get the indices of SNPs to exclude from ref-sumstats harmonization
    ref_remove_ind <- unique(c(ref_monomorphic, ref_dup))
    ref_keep_ind <- setdiff(seq_len(ncol(x)), ref_remove_ind)
    message("Excluded ", length(ref_remove_ind),
            " SNPs from the harmonization.")
    message("-----")

    ## Check mandatory columns for `ID` matching
    has_columns(sumstats, c("id", "pvalue"))

    ## Input SNP IDs must be unique
    if (anyDuplicated(sumstats$id) > 0L) {
      stop(
        "Duplicate SNP IDs are found in ",
        "'", sumstats_name, "'. ",
        "If '", sumstats_name, "' ",
        "is being matched by `id`, SNP IDs must be unique.",
        call. = FALSE
      )
    }

    message("Checking the GWAS summary statistics...")
    message(pretty_num(nrow(sumstats)), " variants to be matched.")

    ## Inner join
    sumstats <- x@snps[ref_keep_ind, ][sumstats[, .(id, pvalue)], on = .(id),
                                       nomatch = NULL]

    ## Quick safety check; not strictly necessary
    stopifnot(anyDuplicated(sumstats) == 0L)

    message(pretty_num(nrow(sumstats)),
            " variants have been matched.")

  } else {

    ## Check duplicate SNP IDs in the reference data
    ## Keep one instance of duplicate-ID variants
    ref_dup <- which(duplicated((x@snps[, .(chr, pos, A1, A2)])))
    message("Found ", length(ref_dup),
            " duplicate SNPs in the reference data",
            " by genomic position and alleles codes.")

    ## Get the indices of SNPs to exclude from ref-sumstats harmonization
    ref_remove_ind <- unique(c(ref_monomorphic, ref_dup))
    ref_keep_ind <- setdiff(seq_len(ncol(x)), ref_remove_ind)
    message(
      "Excluded ", length(ref_remove_ind), " SNPs from the harmonization."
    )
    message("-----")

    message("Checking the GWAS summary statistics...")

    ## Check mandatory columns for `chr:pos:A1:A2` matching
    has_columns(sumstats, c("chr", "pos", "A1", "A2", "pvalue"))

    ## Input `chr:pos:A1:A2` combination must be unique
    if (anyDuplicated(sumstats[, .(chr, pos, A1, A2)])) {
      stop(
        "Some variants in ",
        "'", sumstats_name, "' ",
        "share the same genomic position and allele codes. ",
        "Remove duplicate SNPs in ",
        "'", sumstats_name, "'.",
        call. = FALSE
      )
    }

    ## Matching taking into account ref allele swaps (optionally strand flip)
    sumstats <- snp_match(sumstats, x@snps[ref_keep_ind, ], check_strand_flip,
                          return_indice)
  }

  ## Return
  if (!return_indice) {
    sumstats[, .(id, chr, pos, A1, A2, pvalue)]
  } else {
    if (!check_strand_flip) {
      sumstats[, .(id, chr, pos, A1, A2, pvalue, key_, swapped_)]
    } else {
      sumstats[, .(id, chr, pos, A1, A2, pvalue, key_, swapped_, flipped_)]
    }
  }

}

