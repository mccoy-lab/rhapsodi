#' A function to assign each allele to a haplotype or an NA if not enough information is known
#' 
#' This function reads along each gametes and replaces its read (0 or 1) with the corresponding
#' haplotype (h1 or h2) based on the allele in each donor haplotype; alternatively, the read 
#' may be replaced with an NA if not enough information was known for phasing
#' 
#' @param dt Input matrix of gamete reads
#' @param complete_haplotypes Inferred donor haplotypes
#' 
#' @return dt matrix of gametes coded by haplotype at each position 
#' 
recode_gametes <- function(dt, complete_haplotypes) {
  for (i in 1:ncol(dt)) {
    dt[i][dt[i] == complete_haplotypes$h1] <- "hap1"
    dt[i][dt[i] == complete_haplotypes$h2] <- "hap2"
    dt[c(which(dt[,i] == 0 | dt[,i] == 1)),i] <- NA
  }
  return(dt)
}