#' Simulated Sperm-seq data A/C/G/T/NA encoded
#'
#' Simulated Sperm-seq data from the generative model, with data A/C/G/T/NA encoded,
#' with 50 gametes, 5000 beginning SNPs, 0.1 coverage, 
#' average recombination rate of 1, sequencing error rate of 0.005.
#' Following filtering, there are 4136 hetSNPs. This data originated from the
#' generative model with a random seed of 42, originally 0/1/NA encoded.
#' Then `dataACGT.R` from the `data-raw` directory was run on the 0/1/NA 
#' encoded data to produce this A/C/G/T/NA encoded version.
#' The first column is the SNP index positions 
#' the second column is the reference allele, the third the alternate allele,
#' and the following 50 columns are the gamete genotypes which are encoded with A/C/G/T/NA.
#' 
#' @docType data
#' 
#' @usage dataACGT
#' 
#' @format A dataframe with column names, 4136 rows, and 53 columns:
#' \describe{
#'   \item{positions}{positions, integer, SNP index pre-filtering for hetSNPs only}
#'   \item{ref}{ref, a character or string, the reference allele at that SNP}
#'   \item{alt}{alt, a character or string, the alternate allele at that SNP}
#'   \item{gami_}{gami_, A, C, G, T, or NA, the genotype at that SNP for gamete i}
#' }
#' 
#' @keywords datasets
#' 
"dataACGT"