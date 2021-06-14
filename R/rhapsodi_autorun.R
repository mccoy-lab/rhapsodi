#'
#'
#'
rhapsodi_autorun <- function(input_file, acgt = FALSE, threads=2, sampleName="sampleT", chrom = "chrT", outdir = "tmp/", seqError_model = 0.005, avg_recomb_model = 1, 
                             window_length=3000, overlap_denom = 2, mcstop=TRUE, stringent_stitch=TRUE, stitch_new_min=0.5,
                             smooth_imputed_genotypes=FALSE, smooth_crossovers=TRUE){
  if (!acgt){
    input_data <- standard_input(input_file)
    dt <- input_data$dt
    positions <- input_data$positions
  } else{
    #will be filled in later
    input_data <- 0
  }
  complete_haplotypes <- impute_donor_haplotypes(dt, positions, window_length = window_length, overlap_denom = overlap_denom, threads=threads, mcstop=mcstop, stringent_stitch=stringent_stitch, stitch_new_min = stitch_new_min)
  filled_gametes <- fill_gametes(dt, complete_haplotypes, sequencing_error = seqError_model, avg_recomb = avg_recomb_model, threads = threads)
  rhapsodi_out <- report_gametes(smooth_crossovers, smooth_imputed_genotypes, complete_haplotypes, dt, filled_gametes, sampleName, chrom)
  return (rhapsodi_out)
}