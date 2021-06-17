#' A function to read in the sparse gamete sequencing data encoded in 0/1 tab-delimited format with the first column containing SNP positions in integer format
#' 
#' This function is called when data is already encoded in 0/1 format. The file is expected to be a tab-delimited file
#' And to have the SNP genomic positions in the first column
#' We assume that a single file originates from a single sample and chromosome
#' Finally, we subset to include only heterozygous SNPs or those that have atleast one reference )0) and one alternate allele (1)
#' 
#' @param input_file the path and file name for the tab-delimited 0/1 encoded file
#' 
#' @return input in a named list format with `dt` and `positions` where `dt` is a data frame of the heterozygous sparse gamete data and `positions` is a vector of the SNP genomic positions 
#' 
#' @export
#'
standard_input <- function(input_file){
  stopifnot(file_test("-f", input_file), paste0("The provided file ", input_file, " does not exist."))
  dt <- read_delim(input_file, delim = "\t", col_types = cols(.default = "d")) %>% 
    as.data.frame()
  
  dt <- dt[((rowSums(dt[,2:ncol(dt)] == 0, na.rm = TRUE) > 0) & (rowSums(dt[, 2:ncol(dt)] == 1, na.rm=TRUE) > 0)),]
  
  positions <- dt[,1]
  return(list(dt=dt, positions=positions)) 
}