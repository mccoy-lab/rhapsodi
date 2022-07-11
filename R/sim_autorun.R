#' This is the simulation autorun that generates data, calls rhapsodi, and compares rhapsodi's output with the truth/generated data
#' 
#' First this simulation autorun function calls the generative model producing input data for rhapsodi as well as the true fully known data to be used later in assessment of rhapsodi
#' Then it calls the rhapsodi autorun function which passes the generated output through the rhapsodi steps, returning rhapsodi's predictions
#' Then it calls assessment using the function `assess_it` to compare the fully known truth/generated data with the rhapsodi predictions. 
#' which first assesses donor haplotypte phasing, producing a named list with single values for lhs (largest haplotype segment), ser (switch error rate), acc (accuracy), com (completeness)
#' Then second assesses gamete genotype imputation, producing a named list with vectors for lhs (largest haplotype segment), ser (switch error rate), acc (accuracy), com (completeness)
#' Then third assesses recombination discovery producing a named list with single values for precision, recall, accuracy, specificity, fdr (false discovery rate), fpr (false positive rate) f1 (f1 score), 
#' true_n (number of true recombination breakpoints), pred_n (number of predicted recombination breakpoints), tn (true negative), fn (false negative), tp (true positive), fp (false positive)
#' Finally, it returns a list of named lists where `phasing` contains the phasing assessment named list
#' `gam_imputation` contains the gamete genotype imputation assessment named list
#' and `recomb` contains the recombination breakpoint discovery assessment named list
#'
#' @param num_gametes an integer, the number of gametes, used within `sim_run_generative_model`, or the number of columns for the sparse gamete data you want generated
#' @param num_snps an integer, the number of SNPs, used within `sim_run_generative_model`,or the number of rows for the sparse gamete data you want generated. Note: not all of these will be heterozygous due to the coverage and therefore this number won't necessarily equal the number of SNPs following filtering at the end of the generation
#' @param coverage a numeric, input if `input_cov` is TRUE, suggested NULL otherwise, used within `sim_run_generative_model`
#' @param threads an integer, default = 2, the number of cores to use when we use `mclapply` throughout the rhapsodi pipeline
#' @param recomb_lambda a numeric, the average rate of recombination expected for the simulation, used by `sim_find_num_crossovers` within `sim_run_generative_model`
#' @param random_seed an integer, default = 42, the random seed which will be set for the simulation, used within `sim_run_generative_model`
#' @param input_cov a bool, TRUE if coverage (i.e. like 0.01 (x)) will be input rather than missing genotype rate, default = TRUE, used within `sim_run_generative_model`
#' @param input_mgr a bool, TRUE if missing genotype rate (i.e. like 80 (%) or 0.8) will be input rather than coverage, default = FALSE, used within `sim_run_generative_model`
#' @param missing_genotype_rate a numeric, input if `input_mgr` is TRUE and `input_cov` is FALSE, suggested NULL otherwise, default = NULL, used within `sim_run_generative_model`
#' @param add_seq_error a bool, TRUE if you want to simulate sequencing error in the generated data, default = TRUE, used within `sim_run_generative_model`
#' @param seqError_add a numeric, default = 0.005, the sequencing error rate if adding sequencing error to the generated data, used by `sim_add_seq_error` within `sim_run_generative_model`
#' @param add_de_novo_mut a bool, TRUE if you want to add de novo mutations to the generated data, default = FALSE, used within `sim_run_generative_model`
#' @param de_novo_lambda an integer, default = 5, parameterizes a poisson distribution to find the number of de novo mutations (DNM) total, used by `sim_add_de_novo_mut` within `sim_run_generative_model`
#' @param de_novo_alpha a numeric, default=7.5, shape parameter for a gamma distribution to find the number of gametes affected per DNM, used by `sim_add_de_novo_mut` within `sim_run_generative_model`
#' @param de_novo_beta a numeric, default = 10, scale parameter for a gamma distribution to find the number of gametes affected per DNM, used by `sim_add_de_novo_mut` within `sim_run_generative_model`
#' @param cons a bool, default = FALSE, If TRUE, compares recombination breakpoints in a conservative manner such that if two ore more true breakpoints intersect a single prediction, we only consider one intersection to be a tp and the rest to be fn. If FALSE, all are tp. Used by `sim_find_tp` and `sim_find_fn` within `sim_assess_recomb` by `sim_assess_it`
#' @param sampleName a string, default = "sim", only used in reporting recombination breakpoints by rhapsodi. Used by `find_recomb_spots` and `report_gametes` within `rhapsodi_autorun`
#' @param chrom a string, default = "chrS", only used in reporting recombination breakpoints by rhapsodi. Used by `find_recomb_spots` and `report_gametes` within `rhapsodi_autorun`
#' @param seqError_model a numeric, default = 0.005, used in `build_hmm` within `fill_gametes` in `rhapsodi_autorun`, the expected error rate in sequencing, affects emission probabilities
#' @param avg_recomb_model a numeric, default = 1, used in `buil_hmm` within `fill_gametes` in `rhapsodi_autorun`, the expected number of average recombination spots per chromosome, affects transition probabilities
#' @param window_length an integer, default = 3000, used in `split_with_overlap` within `impute_donor_haplotypes` in `rhapsodi_autorun` as the segment length to use in constructing overlapping windows for phasing
#' @param overlap_denom an integer, default = 2, used in `split_with_overlap` within `impute_donor_haplotypes` in `rhapsodi_autorun` as the denominator in calculating the amount of overlap to use in the overlapping windows for phasing
#' @param mcstop a bool, default = FALSE, used in `stitch_haplotypes` within `impute_donor_haplotypes` in `rhapsodi_autorun`, only considered if `stringent_stitch` is TRUE. This parameter is used to determine whether phasing continues or exits if the mean concordance bewteen two windows is between 0.1 and 0.9. We prefer TRUE or continuing despite discordance for the simulation section so rhapsodi continues asking which threshold is closer and acts accordingly
#' @param stringent_stitch a bool, default = TRUE, used in `stitch_haplotypes` within `impute_donor_haplotypes` in `rhapsodi_autorun`, this parameter is used to determine the threshold values used in determining whether two windows originate from the same donor. If TRUE, the preset thresholds of 0.1 and 0.9 are used.
#' @param stitch_new_min a numeric, >0, < 1, default = 0.5, used in `stitch_haplotypes` within `impute_donor_haplotypes` in `rhapsodi_autorun`, only evaluated if `stringent_stitch` is FALSE. This parameter is dually assigned as the `different_max` and `same_min` threshold values when considering the mean concordance between two overlapping windows in phasing
#' @param smooth_imputed_genotypes a bool, default = FALSE, used in `report_gametes` and `unsmooth` within `rhapsodi_autorun`, whether to use smoothed data from the HMM or original reads when there are mismatches for ending/predicted filled gamete genotypes. If TRUE, doesn't replace smoothed data from HMM with original reads
#' @param fill_ends a boolean; if TRUE, fills the NAs at the terminal edges of chromosomes with the last known or imputed SNP (for end of chromosome) and the first known or imputed SNP (for beginning of chromosome); if FALSE, leaves these genotypes as NA; default = TRUE 
#' @param smooth_crossovers a bool, default = TRUE, used in `report_gametes` and `unsmooth` within `rhapsodi_autorun`, whether to use smoothed data from the HMM or original reads when there are mismatches for recombination breakpoint discovery. If TRUE, doesn't replace smoothed data from HMM with original reads
#' @param verbose a bool; default is FALSE; if TRUE, prints progress statements after each step is successfully completed
#' 
#' @return all_metrics a named list of named lists with all the assessment metric values or vectors 
#'  
#' @export 
#'
sim_autorun <- function(num_gametes, num_snps, coverage, 
                        recomb_lambda, random_seed=42, threads=2,
                        input_cov=TRUE, input_mgr=FALSE, missing_genotype_rate=NULL,
                        add_seq_error=TRUE, seqError_add=0.005,
                        add_de_novo_mut=FALSE, de_novo_lambda=5, de_novo_alpha=7.5, de_novo_beta=10,
                        cons=FALSE, sampleName="sim", chrom="chrS", seqError_model=0.005, avg_recomb_model=1,
                        window_length=3000, overlap_denom=2, mcstop = FALSE, stringent_stitch = TRUE, stitch_new_min = 0.5, 
                        smooth_imputed_genotypes=FALSE, fill_ends = TRUE, smooth_crossovers=TRUE, verbose = FALSE){
  #Generate simulated data
  generated_data <- sim_run_generative_model(num_gametes, num_snps, coverage, 
                                             recomb_lambda, random_seed = random_seed, 
                                             input_cov = input_cov, input_mgr =  input_mgr, missing_genotype_rate = missing_genotype_rate,
                                             add_seq_error = add_seq_error, seqError_add=seqError_add, 
                                             add_de_novo_mut=add_de_novo_mut, de_novo_lambda=de_novo_lambda, de_novo_alpha=de_novo_alpha, de_novo_beta=de_novo_beta)
  if (verbose){ message("Input Data Simulated")}
  #Run rhapsodi on sparse simulated data (generated_data$gam_na)
  rhapsodi_out <- rhapsodi_autorun(NULL, use_dt = TRUE, input_dt = generated_data$gam_na, threads=threads, sampleName=sampleName, chrom=chrom, seqError_model=seqError_model, avg_recomb_model=avg_recomb_model,
                                   window_length=window_length, overlap_denom=overlap_denom, mcstop = mcstop, stringent_stitch = stringent_stitch, stitch_new_min = stitch_new_min,
                                   smooth_imputed_genotypes=smooth_imputed_genotypes, fill_ends = fill_ends, smooth_crossovers=smooth_crossovers, verbose = verbose)
  
  if (verbose) {message(" rhapsodi run complete")}
   #Assess how rhapsodi did
  if (smooth_imputed_genotypes){
    if (verbose) { message(" checking metrics")}
    all_metrics <- sim_assess_it(generated_data$donor_haps, rhapsodi_out$donor_haps, generated_data$recomb_spots, rhapsodi_out$recomb_breaks, generated_data$gam_full, rhapsodi_out$gamete_genotypes[,-c(1,2)], cons = cons, verbose = verbose)
  } else{
    if (verbose) { message(" checking metrics")}
    all_metrics <- sim_assess_it(generated_data$donor_haps, rhapsodi_out$donor_haps, generated_data$recomb_spots, rhapsodi_out$recomb_breaks, generated_data$gam_full, rhapsodi_out$unsmoothed_gamete_genotypes[,-c(1,2)], cons=cons, verbose = verbose)
  }
  if (verbose){ message("rhapsodi performance on simulated data assessed")}
  return(all_metrics)
}