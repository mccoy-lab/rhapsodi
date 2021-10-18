#' Simulated Sperm-seq data 0/1/NA encoded
#'
#' Simulated Sperm-seq data from the generative model, with data 0/1/NA encoded,
#' with 50 gametes, 5000 beginning SNPs, 0.1 coverage, 
#' average recombination rate of 1, sequencing error rate of 0.005.
#' Following filtering, there are 4136 hetSNPs. This data originated from the
#' generative model with a random seed of 42. 
#' The first column is the SNP index positions and the following 50 columns
#' are the gamete genotypes which are encoded with 0/1/NA.
#' 
#' @docType data
#' 
#' @usage sim01NA42
#' 
#' @format A dataframe with column names, 4136 rows, and 51 columns:
#' \describe{
#'   \item{positions}{positions, integer, SNP index pre-filtering for hetSNPs only}
#'   \item{gami_}{gami_, 0, 1, or NA, the genotype at that SNP for gamete i}
#' }
#' 
#' @keywords datasets
#' 
"sim01NA42"