#'
#'
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
  return(out)
}