#'
#'
#'
sim_autorun <- function(){
  #Generate simulated data
  generated_data <- sim_run_generative_model(__lots_of_args__)
  #Run rhapsodi on sparse simulated data (generated_data$gam_na)
  
  #Assess how rhapsodi did
  all_metrics <- sim_assess_it(generated_data$donor_haps, __pred_donor_haps__, generated_data$recomb_spots, __pred_recomb__, generated_data$gam_full, __pred_gam__, write_out_plot = write_out_plot, cons = cons)
  return(all_metrics)
}