#' This function can be used to run all steps of rhapsodi
#' 
#' This function runs all steps of rhapsodi by first inputting the data
#' Then calling `phase_donor_haplotypes` to run donor phasing
#' Then calling `impute_gamete_genotypes` to run gamete genotype imputation through the HMM
#' And finally calling `discover_meiotic_recombination` to run meiotic recombination discovery
#' The named list which is returned has 
#' `donor_haps` which is the phased haplotypes as a data frame with column names index, pos (for SNP positions), h1 (haplotype 1), & h2 (haplotype 2)
#' `gamete_haps` which is the filled gametes specifying from which donor haplotype each gamete position originates in data frame form
#' `gamete_genotypes` which is the filled gametes specifying the genotype in (0's and 1's) for each gamete position as a data frame
#' `recomb_breaks`which is a data frame specifying the recombination breakpoints for each gamete
#' 
#' @param input_file a string; the path plus filename for the input sparse gamete genotype data in tabular form. Note the form is different depending on the value of `acgt`. Use NULL if `use_dt` is TRUE
#' @param use_dt a bool; default is FALSE, whether to input a pre-loaded data frame/table rather than using an input file
#' @param input_dt a data frame/table; only necessary if use_dt is TRUE. User-pre-loaded data frame/table. Note the format is different depending on the value of `acgt` 
#' @param acgt a bool; default is FALSE; If TRUE, assumes that the data is not 0/1 encoded and is available in pre-loaded data frame format passed with `input_dt`
#' @param threads an integer; default is 2, number of threads to utilize when we use `mclapply`
#' @param sampleName a string; default is "sampleT", fill in with whatever the sample name is. We assume a single input file is from a single file
#' @param chrom a string; default is "chrT", fill in with whatever the chromosome is. We assume a single input file is from a single chromosome
#' @param seqError_model a numeric; default is 0.005, used in `build_hmm` within `fill_gametes`, the exppected error in sequencing
#' @param avg_recomb_model a numeric; default is 1, used in `build_hmm` within `fill_gametes`, the expected number of average recombination spots per chromosome
#' @param window_length an integer; default is 3000, used in `split_with_overlap` within `impute_donor_haplotypes`, the segment length to use in constructing overlapping windows for phasing
#' @param overlap_denom an integer; default is 2, used in `split_with_overlap` within `impute_donor_haplotypes`, User-input value for denominator in calculation of overlap, or the degree of overlap between segments
#' @param mcstop a bool; used in `stitch_haplotypes` within `impute_donor_haplotypes`, only considered if `stringent_stitch` is TRUE; default is TRUE; this parameter is used to determine whether phasing continues or exits if the mean concordance between two windows is between 0.1 and 0.9. If TRUE, rhapsodi exits. If FALSE, rhapsodi and phasing continues, asking which threshold the concordance is closer to and acting accordingly
#' @param stringent_stitch a bool; used in `stitch_haplotypes` within `impute_donor_haplotypes`, default is TRUE, this parameter is used to determine the threshold values used in determining whether two windows originate from the same donor. If TRUE, the preset thresholds of 0.1 and 0.9 are used.
#' @param stitch_new_min a numeric >0, but <1; default is 0.5; used in `stitch_haplotypes` within `impute_donor_haplotypes`, this parameter is only evaluated if `stringent_stitch` is FALSE and is dually assigned as the `different_max` and `same_min` threshold values when considering the concordance between two windows and therefore which donors they originate from (same or different).
#' @param smooth_imputed_genotypes a bool; default is FALSE; used in `report_gametes` whether to use smoothed data from the HMM or original reads for ending filled gamete genotypes. If TRUE, doesn't replace smoothed data from HMM with original reads when there's a mismatch
#' @param smooth_crossovers a bool; default is TRUE; used in `report_gametes` whether to use smoothed data from the HMM or original reads for recombination finding. If TRUE, doesn't replace smoothed data from HMM with original reads when there's a mismatch
#'
#' @return rhapsodi_out a named list with `donor_haps`, `gamete_haps`, `gamete_genotypes`, `unsmoothed_gamete_haps`,  `unsmoothed_gamete_genotypes`, and `recomb_breaks`
#'
#' @export
#'
rhapsodi_autorun <- function(input_file, use_dt = FALSE, input_dt = NULL, acgt = FALSE, threads=2, sampleName="sampleT", chrom = "chrT", seqError_model = 0.005, avg_recomb_model = 1, 
                             window_length=3000, overlap_denom = 2, mcstop=TRUE, stringent_stitch=TRUE, stitch_new_min=0.5,
                             smooth_imputed_genotypes=FALSE, smooth_crossovers=TRUE){
  input_data <- read_data(input_file, use_dt = use_dt, input_dt = input_dt, acgt = acgt)
  message("Data Input Completed")
  complete_haplotypes <- phase_donor_haplotypes(input_data$dt, input_data$positions, window_length = window_length, overlap_denom = overlap_denom, threads=threads, mcstop=mcstop, stringent_stitch=stringent_stitch, stitch_new_min = stitch_new_min)
  message("Phasing Completed")
  filled_gametes <- impute_gamete_genotypes(input_data$dt, complete_haplotypes, input_data$positions, sequencing_error = seqError_model, avg_recomb = avg_recomb_model, smooth_imputed_genotypes = smooth_imputed_genotypes, threads = threads)
  message("Imputation Completed")
  recomb_breaks <- discover_meiotic_recombination(input_data$dt, complete_haplotypes, filled_gametes, input_data$positions, smooth_crossovers = smooth_crossovers, smooth_imputed_genotypes = smooth_imputed_genotypes, sampleName = sampleName, chrom=chrom, threads = threads)
  message("Discovery Completed")
  if (acgt){
    complete_haplotypes <- record_ref_alt(complete_haplotypes, input_data$ref, input_data$alt)
    filled_gametes$filled_gametes <- record_ref_alt(filled_gametes$filled_gametes, input_data$ref, input_data$alt)
    filed_gametes$unsmoothed_gametes <- record_ref_alt(filled_gametes$unsmoothed_gametes, input_data$ref, input_data$alt)
  }
  rhapsodi_out <- list(donor_haps = complete_haplotypes, gamete_haps = filled_gametes$filled_gametes_haps , gamete_genotypes = filled_gametes$filled_gametes , unsmoothed_gamete_haps = filled_gametes$unsmoothed_gametes_haps, unsmoothed_gamete_genotypes = filled_gametes$unsmoothed_gametes, recomb_breaks = recomb_breaks)
  return (rhapsodi_out)
}