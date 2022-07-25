#' This function can be used to run all steps of rhapsodi
#' 
#' This function runs all steps of rhapsodi by first inputting the data
#' Then calling `phase_donor_haplotypes` to run donor phasing
#' Then calling `impute_gamete_genotypes` to run gamete genotype imputation
#' And finally calling `discover_meiotic_recombination` to run meiotic recombination discovery
#' Input data should be sparse gamete genotype data encoded either as 0/1/NA or as a VCF style input with A/C/G/T/NA. The data is either from a tab-delimited file with a header or a pre-loaded data frame/table.
#' For both input types, the first column should contain SNP positions in integer format. 
#' For ACGT input type, specifically (acgt = TRUE), the second column should be the REF allele and the third column should be the ALT allele. All following columns should be gamete data, with each gamete having its own column. Within these columns, data should be A/C/G/T/NA
#' For 0/1/NA input type, specifically (acgt = FALSE), gamete data starts in the second column and continues for the rest of the columns.
#' After rhapsodi has completed all three tasks, it returns a named list which has 
#' `donor_haps` which is the phased haplotypes as a data frame with column names index, pos (for SNP positions), h1 (haplotype 1), & h2 (haplotype 2) if acgt = FALSE. Otherwise: index, pos, a0, a1, h1, h2
#' `gamete_haps` which is the filled gamete data frame specifying from which donor haplotype each gamete position originates. Column names: index, pos, gamete_names.
#' `gamete_genotypes` which is the filled gamete dataf rame specifying the genotype (in 0's and 1's) for each gamete position. If acgt = FALSE, column names: index, pos, gamete_names. Otherwise: index, pos, a0, a1, gamete_names
#' `unsmoothed_gamete_haps` which is the filled gamete data frame specifying from which donor haplotype each gamete position originates in data frame form, after unsmoothing the data by replacing imputed values with original sequencing reads when there's disagreement between observations and imputation. Column names: index, pos, gamete_names.
#' `unsmoothed_gamete_genotypes` which is the filled gamete data frame specifying the genotype (in 0's and 1's) for each gamete position, after unsmoothing the the data by replacing imputed values with original sequencing reads when there's disagreement between observations and imputation. If acgt = FALSE, column names: index, pos, gamete_names. Otherwise: index, pos, a0, a1, gamete_names
#' `recomb_breaks`which is a data frame specifying the recombination breakpoints for each gamete. Column names: Ident, Genomic_start, Genomic_end
#' 
#' @param input_file a string; the path plus filename for the input sparse gamete genotype data in tabular form. Note the form is different depending on the value of `acgt`. Use NULL if `use_dt` is TRUE
#' @param use_dt a bool; default is FALSE, whether to input a pre-loaded data frame/table rather than using an input file
#' @param input_dt a data frame/table; only necessary if use_dt is TRUE. User-pre-loaded data frame/table. Note the format is different depending on the value of `acgt` 
#' @param acgt a bool; default is FALSE; If TRUE, assumes that the data is not 0/1/NA encoded, rather gamete genotypes are A/C/G/T/NA encoded and the dataframe has ref and alt columns.
#' @param threads an integer; default is 2, number of threads to utilize when we use `mclapply`like functions
#' @param sampleName a string; default is "sampleT", fill in with whatever the sample name is. We assume a single input file is from a single sample/donor
#' @param chrom a string; default is "chrT", fill in with whatever the chromosome is. We assume a single input file is from a single chromosome
#' @param seqError_model a numeric; default is 0.005, used in `build_hmm` within `impute_gamete_genotypes`, the expected error rate in genotyping
#' @param avg_recomb_model a numeric; default is 1, used in `build_hmm` within `impute_gamete_genotypes`, the expected number of average recombination events per chromosome
#' @param window_length an integer; default is 3000, used in `split_with_overlap` within `phase_donor_haplotypes`, the segment length to use in constructing overlapping windows for phasing
#' @param overlap_denom an integer; default is 2, used in `split_with_overlap` within `phase_donor_haplotypes`, User-input value for denominator in calculation of overlap, or the degree of overlap between segments
#' @param calculate_window_size_bool A bool; used in `phase_donor_haplotypes`, whether or not to calculate the window size based on characteristics of the input dataset; default = FALSE
#' @param estimated_coverage a numeric; used in `calculate_window_size` within `phase_donor_haplotypes` only if the user wants rhapsodi to calculate the preferred window size for accurate phasing given characteristics of the data; the estimated sequencing depth of coverage of the input data; default = NULL 
#' @param mcstop a bool; used in `stitch_haplotypes` within `phase_donor_haplotypes`, only considered if `stringent_stitch` is TRUE; default is TRUE; this parameter is used to determine whether phasing continues or exits if the mean concordance between two windows is between 0.1 and 0.9. If TRUE, rhapsodi exits. If FALSE, rhapsodi and phasing continues, asking which threshold the concordance is closer to and acting accordingly
#' @param stringent_stitch a bool; used in `stitch_haplotypes` within `phase_donor_haplotypes`, default is TRUE, this parameter is used to determine the threshold values used in determining whether two windows originate from the same donor. If TRUE, the preset thresholds of 0.1 and 0.9 are used.
#' @param stitch_new_min a numeric >0, but <1; default is 0.5; used in `stitch_haplotypes` within `phase_donor_haplotypes`, this parameter is only evaluated if `stringent_stitch` is FALSE and is dually assigned as the `different_max` and `same_min` threshold values when considering the concordance between two windows and therefore which donors they originate from (same or different).
#' @param smooth_imputed_genotypes a bool; default is FALSE; used in `impute_gamete_genotypes` whether to use smoothed data from the HMM or original reads for the ending filled gamete genotypes, whenever there is disagreement between the two. If TRUE, doesn't replace smoothed data from HMM with original reads when there's a mismatch
#' @param fill_ends a boolean; if TRUE, fills the NAs at the terminal edges of chromosomes with the last known or imputed SNP (for end of chromosome) and the first known or imputed SNP (for beginning of chromosome); if FALSE, leaves these genotypes as NA; default = TRUE 
#' @param smooth_crossovers a bool; default is TRUE; used in `discover_meiotic_recombination` whether to use smoothed data from the HMM or original reads for recombination finding. If TRUE, doesn't replace smoothed data from HMM with original reads when there's a mismatch
#' @param verbose a bool; default is FALSE; if TRUE, prints progress statements after each step is successfully completed
#'
#' @return rhapsodi_out a named list with `donor_haps`, `gamete_haps`, `gamete_genotypes`, `unsmoothed_gamete_haps`,  `unsmoothed_gamete_genotypes`, and `recomb_breaks`
#'
#' @export
#'
rhapsodi_autorun <- function(input_file, use_dt = FALSE, input_dt = NULL, acgt = FALSE, threads=2, sampleName="sampleT", chrom = "chrT", seqError_model = 0.005, avg_recomb_model = 1, 
                             window_length=3000, overlap_denom = 2, calculate_window_size_bool = FALSE, estimated_coverage = NULL,
                             mcstop=TRUE, stringent_stitch=TRUE, stitch_new_min=0.5,
                             smooth_imputed_genotypes=FALSE, fill_ends = TRUE, smooth_crossovers=TRUE, verbose = FALSE){
  input_data <- read_data(input_file, use_dt = use_dt, input_dt = input_dt, acgt = acgt)
  if (verbose){ message("Data Input Completed")}
  complete_haplotypes <- phase_donor_haplotypes(input_data$dt, input_data$positions, window_length = window_length, overlap_denom = overlap_denom, threads=threads, mcstop=mcstop, stringent_stitch=stringent_stitch, stitch_new_min = stitch_new_min, calculate_window_size_bool = calculate_window_size_bool, cov = estimated_coverage, ger = seqError_model, avgr = avg_recomb_model)
  if (verbose){ message("Phasing Completed")}
  filled_gametes <- impute_gamete_genotypes(input_data$dt, complete_haplotypes, input_data$positions, genotyping_error = seqError_model, avg_recomb = avg_recomb_model, smooth_imputed_genotypes = smooth_imputed_genotypes, fill_ends = fill_ends, threads = threads)
  if (verbose){message("Imputation Completed")}
  recomb_breaks <- discover_meiotic_recombination(input_data$dt, complete_haplotypes, filled_gametes, input_data$positions, smooth_crossovers = smooth_crossovers, smooth_imputed_genotypes = smooth_imputed_genotypes, sampleName = sampleName, chrom=chrom, threads = threads)
  if (verbose){message("Discovery Completed")}
  rhapsodi_out <- export_data(input_data, complete_haplotypes, filled_gametes, recomb_breaks, acgt = acgt)
  return (rhapsodi_out)
}