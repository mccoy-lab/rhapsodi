#' A function to filter user input data
#' 
#' This function takes in a datatable of user input data, with genomic positions, and either a file of genomic positions that should be included or a file of genomic positions that should be
#' excluded (or both!). It returns list of approved genomic positions. The user can then filter their raw data to only include regions that are within these approved ranges.
#'  
#' 
#' @export
#' 
#' @param input_data_positions Bed-style file (with chr, start, end columns) showing genomic positions in the raw data   
#' @param chrom Number of chromosome (can iterate through e.g., by calling in lapply)
#' @param exclude_genomic_positions Bed-style file (with chr, start, end columns) showing genomic positions that should be excluded (e.g., are in a blacklisted region of the genome) 
#' @param include_genomic_positions Bed-style file (with chr, start, end columns) showing genomic positions that should be included (e.g., are in an accessible region of the genome)
#' 
#' @return approved_genomic_positions Genomic positions that should be included   
#' 
#' @example 
#' R code here showing my function works 
#' to_include_output <- filter_problem_genome(positions_nc26chr21, 21, blacklist, giab_union)
#' new_dt <- input_dt_nc26chr21[input_dt_nc26chr21$positions %in% to_include_output$positions,]
#' 
filter_problem_genome <- function(input_data_positions, chrom, exclude_genomic_positions, include_genomic_positions) {
  col_chr_name <- paste0("chr", chrom)
  # Define notin function for use with exclusion 
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  # rename for later use in case there is no exclude_genomic_positions argument 
  pos_not_excluded <- input_data_positions
  
  # Find regions in the raw data that are in the blacklist 
  if(!missing(exclude_genomic_positions)) {
    exclude_genomic_positions <- exclude_genomic_positions[exclude_genomic_positions$chr == col_chr_name] 
    # Set key for use in foverlaps
    setkey(exclude_genomic_positions, start, end)
    # Get filtered set of points in the raw data that should be excluded  
    blacklisted <- data.table::foverlaps(new_positions, exclude_genomic_positions, type="any", nomatch=NULL)
    
    # Then, remove those rows (i.e., genomic positions) from the list of positions that will be returned
    pos_not_excluded <- new_positions[new_positions$positions %!in% blacklisted$positions,]
  }
  # Find regions in the raw data that are usable 
  if(!missing(include_genomic_positions)) {
    include_genomic_positions <- include_genomic_positions[include_genomic_positions$chr == col_chr_name,] %>% as.data.table()
    
    # Find overlaps between giab and our data table (excluding the blacklisted sites) 
    setkey(include_genomic_positions, start, end)
    pos_included <- data.table::foverlaps(pos_not_excluded, include_genomic_positions, type="any", nomatch=NULL)
    
    # Then, remove those rows (i.e., genomic positions) from the list of positions that will be returned
    pos_not_excluded <- pos_not_excluded[pos_not_excluded$positions %in% pos_included$positions,]
  }
  return (pos_not_excluded)
}