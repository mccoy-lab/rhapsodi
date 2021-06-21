#' This function adds de novo mutations (DNM) to the generated data
#' 
#' This function adds de novo mutations (DNM) to the generated data, specifically picking the number of DNMs to add (using a poisson distribution),
#' which donor haplotype the DNM originates from (using a uniform distribution)
#' SNP indices from the diploid donor phased haplotypes after which the new DNM will be added (using a uniform distribution),
#' how many gametes and which gametes could potentially be affected by the DNM because that SNP position originates from the affected donor haplotype,
#' and from that how many gametes actually will be affected by each DNM (using a gamma distribution). Using this info, we construct each new row
#' for the diploid donor haplotypes, giving the originating haplotype the alternate allele, and the other the reference allele
#' for the full gamete data, giving unaffected gametes the reference allele, and the affected gametes the alternate allele,
#' For the sparse gamete data, we replace genotypes with NAs for each new row using the missing genotype rate and a uniform distribution
#' Finally, we track SNP indices and adjust any recombination breakpoints as necessary
#' 
#' @param de_novo_lambda an integer, parameterizes a poisson distribution to find the number of DNMs total
#' @param de_novo_alpha a numeric, shape parameter for a gamma distribution to find the maximum number of gametes affected by each DNM
#' @param de_novo_beta a numeric, scale parameter for a gamma distribution to find the maximum number of gametes affected by each DNM
#' @param num_snps an integer, the number of SNPs or the number of rows, the generated data had before calling this function
#' @param num_gametes an integer, the number of gametes, or the number of columns, the generated data has 
#' @param gam_haps data matrix/frame of the hapltoypes from which each SNP in each gamete originates (encoded as 1's and 2's), necessary to find which gametes potentially could be affected by each DNM 
#' @param gam_mat full data matrix/frame of the genotypes by SNP for each gamete, encoded in 0's and 1's
#' @param gam_mat_with_na sparse data/ matrix/frame of the genotypes by SNP for each gamete, encoded in 0's and 1's and NAs
#' @param donor_haps a data frame with the phased diploid donor haplotypes in two columns `donor1` and `donor2`
#' @param unlist_ci a named vector from unlist with the crossover break points for each gamete 
#' @param missing_genotype_rate a numeric, the missing genotype rate of the simulation
#' 
#' @return out a named list with the adjusted `unlist_ci`, `num_snps`, `gam_haps`, `gam_mat`, `gam_mat_with_na`, `donor_haps`, as well as the new `new_rows` which tracks where the new DNMs are
#'
sim_add_de_novo_mut <- function(de_novo_lambda, de_novo_alpha, de_novo_beta, num_snps, num_gametes, gam_haps, gam_mat, gam_mat_with_na, donor_haps, unlist_ci, missing_genotype_rate){
  out <- list()
  stopifnot(de_novo_beta > 0)
  new_rows <- c()
  num_dnm <- rpois(1, de_novo_lambda) + 1 #ensure this is greater than 0 by adding 1
  num_gametes_affected_per_dnm <- ceiling(rgamma(num_dnm, de_novo_alpha, scale=de_novo_beta))
  donors_with_dnm <- sample(1:2, num_dnm, replace=TRUE)
  for (i in 1:num_dnm){ #for every donor dnm
    donor_with_dnm <- donors_with_dnm[i]
    row_loc <- sample(1:num_snps, 1) #pick random row after which we'll add this new de novo mutation
    gam_can_be_affected <- which(gam_haps[row_loc, ]==donor_with_dnm)
    num_gametes_affected <- min(num_gametes_affected_per_dnm[i], length(gam_can_be_affected))
    message(paste0("Number of gametes affected for de novo mutation ", i, ": ", num_gametes_affected))
    where_locs_greater <- which(new_rows >= (row_loc+1))
    new_rows[where_locs_greater] <- new_rows[where_locs_greater] + 1
    new_rows <- c(new_rows, (row_loc+1))
    message(paste0("new row location: ", row_loc+1))
    message(paste0("new rows vector: ", new_rows, collapse= " ; "))
    donor_haps <- rbind(donor_haps[1:row_loc,], rep(0, 2), donor_haps[(row_loc+1):num_snps, ])
    donor_haps[(row_loc+1), donor_with_dnm] <- 1
    affected_gam <- sort(sample(gam_can_be_affected, num_gametes_affected))
    message(paste0("affected gametes: ", affected_gam, collapse= " ; "))
    new_row_haps <- rep(abs((donor_with_dnm - 1)-1)+1, num_gametes)
    new_row_haps[gam_can_be_affected] <- donor_with_dnm
    new_row_vals <- rep(0, num_gametes)
    new_row_vals[affected_gam] <- 1
    #add sparsity in
    num_nas_to_add <- sim_find_num_nas(num_gametes, 1, missing_genotype_rate)
    if (missing_genotype_rate <= 0.5){
      change_indices <- sample(1:num_gametes, num_nas_to_add)
      new_row_vals[change_indices] <- NA
    } else if (missing_genotype_rate > 0.5){
      sparse_row <- rep(NA, num_gametes)
      keep_indices <- sample(1:num_gametes, num_gametes - num_nas_to_add)
      sparse_row[keep_indices] <- new_row_vals[keep_indices]
      new_row_vals <- sparse_row
    }
    #add in new SNP line to gamete data, automatically changing pseudo positions too for when we later make this a dataframe
    gam_mat <- rbind(gam_mat[1:row_loc, ], new_row_vals, gam_mat[(row_loc+1):num_snps, ])
    gam_mat_with_na <- rbind(gam_mat_with_na[1:row_loc, ], new_row_vals, gam_mat_with_na[(row_loc+1):num_snps, ])
    #add in new haps line
    gam_haps <- rbind(gam_haps[1:row_loc, ], new_row_haps, gam_haps[(row_loc+1):num_snps, ])
    #add in new SNP to number of SNPs
    num_snps <- num_snps + 1 #increase num_snps by one
    #Adjust recombination breakpoints
    bool_adj <- unlist_ci >= row_loc+1
    bool_adj[is.na(bool_adj)] <- FALSE
    unlist_ci[bool_adj] <- unlist_ci[bool_adj] + 1
  }
  out$unlist_ci <- unlist_ci
  out$num_snps <- num_snps
  out$gam_haps <- gam_haps
  out$gam_mat <- gam_mat
  out$gam_mat_with_na <- gam_mat_with_na
  out$new_rows <- new_rows
  out$donor_haps <- donor_haps
  return(out)
}