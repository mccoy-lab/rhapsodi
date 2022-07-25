#' A function to automatically calculate a preferred window size for phasing 
#' 
#' This function given the number of gametes (`gametes`), coverage (`coverage`), 
#' genotyping error rate (`se`), and average recombination rate (`avgr`),
#' as well as the number of hetSNPs in the data,
#' is used to automatically calculate a preferred window size for phasing.
#' This is accomplished by using the predict function and a pre-trained beta regression model `betaregmodel_20220718`.
#' The user can call this function themselves outside of rhapsodi's main steps, or the user can specify they would rather rhapsodi
#' calculate a window size based on the input data characteristics. 
#' 
#' @param ngam an integer; number of gametes in the dataset; this should be the number of columns in the input dataframe (excluding the positions column)
#' @param cov a numeric; sequencing depth of coverage applicable to the data
#' @param nsnp an integer; number of hetSNPs in the dataset; this should be the number of rows in the input dataframe (after filtering to include hetSNPs only)
#' @param ger a numeric; genotyping error rate applicable to data
#' @param avgr a numeric; average recombination rate
#' 
#' @return window_size an integer; the preferred window size for accurate phasing given the data characteristics
#' 
#' @import betareg
#' @importFrom stats predict
#' 
#' @export

calculate_window_size <- function(ngam, cov, nsnp, ger, avgr){
  load(betaregmodel_20220718)
  window_prop <- predict(betares, data.frame(gametes = ngam, coverage = cov, se = ger, avgr = avgr))
  window_size = as.integer(window_prop * nsnp)
  return(window_size)
}