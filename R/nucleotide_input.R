#' Function for rhapsodi to convert input of A/C/G/T to 0/1 (ref/alt) 
#' 
#' A function to read in the sparse gamete sequencing data encoded in A/C/G/T tab-delimited format with the first column containing SNP positions in integer format, the second column containing the ref allele, and the third column containing the alt allele
#' 
#' This function is called when data is originally in nucleotide (A/C/G/T) format. The file is expected to be a tab-delimited file and to have the SNP genomic positions in the first column
#' We assume that a single file originates from a single sample and chromosome
#' Finally, we subset to include only heterozygous SNPs or those that have atleast one reference (0) and one alternate allele (1)
#' 
#' @param info_and_gametes dataframe (no header) of genomic position, ref allele, alt allele, and read for each gamete at each position in nucleotide format
#'
#' @return input in a named list format with `gametes_het`, `positions`, `ref`, & `alt` where `gametes_het` is a data frame of the heterozygous sparse gamete data (in 0/1 format) and `positions` is a vector of the SNP genomic positions, `ref` is a vector of the reference allele for each SNP genomic position in A/C/G/T, and `alt` is a vector of the alternate allel for each SNP genomic position in A/C/G/T
#' 
#' @export
#'
nucleotide_input <- function(info_and_gametes) {
  positions <- info_and_gametes[,1] 
  ref <- info_and_gametes[,2] 
  alt <- info_and_gametes[,3] 
  gametes <- info_and_gametes[,4:ncol(info_and_gametes)] 
  gametes <- apply(gametes, 2, function(x) x == alt)
  gametes[gametes==TRUE] <- 1
  # Keep only SNPs that are heterozygous in input data
  gametes_het <- gametes[((rowSums(gametes[,2:ncol(gametes)] == 0, na.rm = TRUE) > 0) & (rowSums(gametes[, 2:ncol(gametes)] == 1, na.rm=TRUE) > 0)),-1]
  return(list(gametes_het = gametes_het, positions = positions, ref=ref, alt=alt))
}