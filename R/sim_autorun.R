#'
#'
#'
sim_autorun <- function(num_gametes, num_snps, coverage, 
                        recomb_lambda, random_seed=42, threads=2,
                        input_cov=TRUE, input_mgr=FALSE, missing_genotype_rate=NULL,
                        add_seq_error=TRUE, seqError_add=0.005,
                        add_de_novo_mut=FALSE, de_novo_lambda=5, de_novo_alpha=7.5, de_novo_beta=10,
                        write_out_plot=FALSE, cons=FALSE,
                        sampleName="sim", chrom="chrS", outdir="tmp/", seqError_model=0.005, avg_recomb_model=1,
                        window_length=3000, overlap_denom=2, mcstop = FALSE, stringent_stitch = TRUE, stitch_new_min = 0.5, 
                        smooth_imupted_genotypes=FALSE, smooth_crossovers=TRUE){
  #Generate simulated data
  generated_data <- sim_run_generative_model(num_gametes, num_snps, coverage, 
                                             recomb_lambda, random_seed = random_seed, 
                                             input_cov = input_cov, input_mgr =  input_mgr, missing_genotype_rate = missing_genotype_rate,
                                             add_seq_error = add_seq_error, seqError_add=seqError_add, 
                                             add_de_novo_mut=add_de_novo_mut, de_novo_lambda=de_novo_lambda, de_novo_alpha=de_novo_alpha, de_novo_beta=de_novo_beta)
  #Run rhapsodi on sparse simulated data (generated_data$gam_na)
  rhapsodi_out <- rhapsodi_autorun(threads=threads, sampleName=sampleName, chrom=chrom, outdir=outdir, seqError_model=seqError_model, avg_recomb_model=avg_recomb_model,
                                   window_length=window_length, overlap_denom=overlap_denom, mcstop = mcstop, stringent_stitch = stringent_stitch, stitch_new_min = stitch_new_min,
                                   smooth_imupted_genotypes=smooth_imputed_genotypes, smooth_crossovers=smooth_crossovers)
  
  #Assess how rhapsodi did
  all_metrics <- sim_assess_it(generated_data$donor_haps, rhapsodi_out$donor_haps, generated_data$recomb_spots, rhapsodi_out$recomb_breaks, generated_data$gam_full, rhapsodi_out$gamete_genotypes, write_out_plot = write_out_plot, cons = cons)
  return(all_metrics)
}