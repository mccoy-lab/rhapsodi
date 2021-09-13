#' A function to read in the sparse gamete sequencing data encoded in 0/1/NA either from a tab-delimited file with a header or a pre-loaded data frame/table with the first column containing SNP positions in integer format and gamete data in the rest of the columns
#' 
#' This function is called when data is already encoded in 0/1/NA format. If using a file, the file is expected to be a tab-delimited file with a header
#' And to have the SNP genomic positions in the first column and each gamete has its own column for the second and following columns
#' We assume that a single file originates from a single sample and chromosome
#' Finally, we subset to include only heterozygous SNPs or those that have atleast one reference )0) and one alternate allele (1)
#' 
#' @param input_file a string; the path and file name for the tab-delimited 0/1 encoded file
#' @param use_dt a bool; default is FALSE, whether to input a pre-loaded data frame/table rather than using an input file
#' @param input_dt a data frame/table; only necessary if use_dt is TRUE. User-pre-loaded data frame/table with the first column being SNP positions (integers) and all the other columns having 0/1/NA encoded values for gametes. 
#' 
#' @return input in a named list format with `dt` and `positions` where `dt` is a data frame of the heterozygous sparse gamete data and `positions` is a vector of the SNP genomic positions 
#' 
#' @importFrom magrittr %>%
#' @importFrom readr read_delim cols
#' @importFrom utils file_test
#' 
#' @export
#'
standard_input <- function(input_file, use_dt=FALSE, input_dt=NULL){
  if (!use_dt){
    stopifnot(file_test("-f", input_file))
    dt <- read_delim(input_file, delim = "\t", col_names = TRUE, col_types = cols(.default = "d")) %>% 
      as.data.frame()
  } else {
    dt <- input_dt
  }
  
  hetRows <- ((rowSums(dt[,-1] == 0, na.rm = TRUE) > 0) & (rowSums(dt[, -1] == 1, na.rm=TRUE) > 0))
  positions <- dt[hetRows,1]
  dt <- dt[hetRows, -1]

  return(list(dt=dt, positions=positions)) 
}