#' Trained beta regression model for automatic phasing window size calculation
#' 
#' The trained beta regression model to be used for automatic phasing window size calculation.
#' The model was trained with an intercept and predictors of number of gametes, coverage,
#' genotyping error rate, and  average recombination rate. Of those, 
#' the number of gametes, coverage, and average recombination rate were significant predictors. 
#' The response variable was the optimal window proportion (i.e., the optimal window size / total number of SNPs).
#' This trained model can be used to predict new values given all the predictors.
#' 
#' @docType data
#' 
#' @usage betaregmodel_20220718
#' 
#' @format an object of class "betareg", i.e., a list with components `coefficients`, `residuals`, `n`, etc.
#' 
#' @keywords datasets 
#' 