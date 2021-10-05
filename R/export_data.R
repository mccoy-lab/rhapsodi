#' A function to export the data from rhapsodi
#' 
#' This function exports the data from rhapsodi by first adding index and pos columns to the gamete output data.
#' Then, continues by checking whether the input data was originally nucleotide encoded (A/C/G/T/NA).
#' If so, the function appends a0 and a1 columns to 0/1 encoded outputs by calling `record_ref_alt`, such that the nucleotide(s) associated with genotype 0 are recorded in column a0
#' and the nucleotide(s) associated with genotype 1 are recorded in column a1
#' Finally, the function returns a named list which has 
#' `donor_haps` which is the phased haplotypes as a data frame with column names index, pos (for SNP positions), h1 (haplotype 1), & h2 (haplotype 2) if acgt = FALSE. Otherwise: index, pos, a0, a1, h1, h2
#' `gamete_haps` which is the filled gamete data frame specifying from which donor haplotype each gamete position originates. Column names: index, pos, gamete_names.
#' `gamete_genotypes` which is the filled gamete dataf rame specifying the genotype (in 0's and 1's) for each gamete position. If acgt = FALSE, column names: index, pos, gamete_names. Otherwise: index, pos, a0, a1, gamete_names
#' `unsmoothed_gamete_haps` which is the filled gamete data frame specifying from which donor haplotype each gamete position originates in data frame form, after unsmoothing the data by replacing imputed values with original sequencing reads when there's disagreement between observations and imputation. Column names: index, pos, gamete_names.
#' `unsmoothed_gamete_genotypes` which is the filled gamete data frame specifying the genotype (in 0's and 1's) for each gamete position, after unsmoothing the the data by replacing imputed values with original sequencing reads when there's disagreement between observations and imputation. If acgt = FALSE, column names: index, pos, gamete_names. Otherwise: index, pos, a0, a1, gamete_names
#' `recomb_breaks`which is a data frame specifying the recombination breakpoints for each gamete. Column names: Ident, Genomic_start, Genomic_end
#'   
#' @param acgt a bool; If TRUE, assumes that the input data was not 0/1/NA encoded, rather gamete genotypes were A/C/G/T/NA encoded and the dataframe had ref and alt columns. Will add these columns back as a0 and a1 for 0/1 encoded output data
#' @param input_data the named list from `read_data`
#' @param complete_haplotypes the dataframe from `phase_donor_haplotypes`
#' @param filled_gametes the named list from `impute_gamete_genotypes`
#' @param recomb_breaks the dataframe from `discover_meitoic_recombination`
#' 
#' @return rhapsodi_out a named list with `donor_haps`, `gamete_haps`, `gamete_genotypes`, `unsmoothed_gamete_haps`,  `unsmoothed_gamete_genotypes`, and `recomb_breaks`
#' 
#' @export
#' 
export_data <- function(acgt, input_data, complete_haplotypes, filled_gametes, recomb_breaks){
  filled_gametes$unsmoothed_gametes_haps <- cbind_pos(filled_gametes$unsmoothed_gametes_haps, input_data$positions)
  filled_gametes$unsmoothed_gametes <- cbind_pos(filled_gametes$unsmoothed_gametes, input_data$positions)
  filled_gametes$filled_gametes_haps <- cbind_pos(filled_gametes$filled_gametes_haps, input_data$positions)
  filled_gametes$filled_gametes <- cbind_pos(filled_gametes$filled_gametes, input_data$positions)
  if (acgt){
    complete_haplotypes <- record_ref_alt(complete_haplotypes, input_data$ref, input_data$alt)
    filled_gametes$filled_gametes <- record_ref_alt(filled_gametes$filled_gametes, input_data$ref, input_data$alt)
    filed_gametes$unsmoothed_gametes <- record_ref_alt(filled_gametes$unsmoothed_gametes, input_data$ref, input_data$alt)
  }
  rhapsodi_out <- list(donor_haps = complete_haplotypes, gamete_haps = filled_gametes$filled_gametes_haps , gamete_genotypes = filled_gametes$filled_gametes , unsmoothed_gamete_haps = filled_gametes$unsmoothed_gametes_haps, unsmoothed_gamete_genotypes = filled_gametes$unsmoothed_gametes, recomb_breaks = recomb_breaks)
  return (rhapsodi_out)
}