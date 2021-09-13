#' A function to read in the sparse gamete sequencing data encoded either as 0/1/NA or as a VCF style input with A/C/G/T/NA. The data is either from a tab-delimited file with a header or a pre-loaded data frame/table.
#' For both input types, the first column should contain SNP positions in integer format. For ACGT input type, the second column should be the REF allele and the third column should be the ALT allele.
#' For 0/1/NA gamete data in the rest of the columns
#'
#' This function is called to process input data for rhapsodi. If using a file, the file is expected to be a tab-delimited file with a header
#' For both input types (ACGT and 01), the first column should contain SNP positions in integer format.
#' For ACGT input type, the second column should be the REF allele and the third column should be the ALT allele. All following columns should be gamete data, with each gamete having its own column
#' For 01 input type, each gamete has its own column for the second and following columns
#' We assume that a single file originates from a single sample and chromosome
#' We subset to include only heterozygous SNPs or those that have atleast one reference (0) and one alternate allele (1)
#' Then, if the data is ACGT encoded, we recode it to be 01 encoded
#'
#' @param input_file a string; the path and file name for the tab-delimited file with header
#' @param use_dt a bool; default is FALSE, whether to input a pre-loaded data frame/table rather than using an input file
#' @param input_dt a data frame/table; only necessary if use_dt is TRUE. User-pre-loaded data frame/table with the first column being SNP positions (integers) and all the other columns being 0/1/NA encoded values for gametes if acgt = FALSE.
#' @param acgt a bool; default is FALSE, whether the input data is in ACGT vcf-like format with columns position (int), REF (character), ALT (character), and then gametes (characters), with gamete data stored as A/C/G/T
#'
#' @return input in a named list format with `dt`, `positions`, `ref`, and `alt` where `dt` is a data frame of the heterozygous sparse gamete data encoded as 0/1/NA and `positions` is a vector of the SNP genomic positions, `ref` and `alt` are NULL if `acgt` is FALSE, otherwise, they are character vectors of the REF and ALT alleles for each hetSNP
#'
#' @importFrom magrittr %>%
#' @importFrom readr read_delim cols
#' @importFrom utils file_test
#'
#' @export
#'
read_data <- function(input_file, use_dt=FALSE, input_dt=NULL, acgt = FALSE){
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

  if (acgt){
    ref <- dt[,1]
    alt <- dt[,2]
    dt <- dt[,-c(1,2)]
    dt <- apply(dt, 2, function(x) x == alt)
    dt[dt==TRUE] <- 1
  } else {
    ref <- NULL
    alt <- NULL
  }
  return(list(dt = dt, positions = positions, ref = ref, alt = alt))
}
