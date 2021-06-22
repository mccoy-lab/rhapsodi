#' A function to separate the gamete SNP positions into overlapping segments 
#' 
#' This function separates each gamete into overlapping windows of size determined by the user. These 
#' overlaps are an input to the distance matrix constructed in `reconstruct_haplotypes`. 
#'
#' Adapated from https://stackoverflow.com/questions/8872376/split-vector-with-overlapping-samples-in-r
#'
#' @param vector_of_positions vector of SNP positions 
#' @param window_size User-input value for length of segment (default = 3000)
#' @param overlap_denom User-input value for denominator in calculation of overlap, or the degree of overlap between segments (default = 2)
#'
#' @return windows list of lists with the number of lists equal to the number of windows and each inner list containing the SNP positions in that overlapping segment
#'
split_with_overlap <- function(vector_of_positions, window_size, overlap_denom) {
  overlap <- window_size %/% overlap_denom #Degree of overlap between segments
  starts = seq(1, length(vector_of_positions), by = window_size - overlap)
  ends   = starts + window_size - 1
  ends[ends > length(vector_of_positions)] = length(vector_of_positions)
  windows <- lapply(1:length(starts), function(i) vector_of_positions[starts[i]:ends[i]])
  
  #merge the last two windows to avoid edge effect
  if (length(windows) > 1){
    combined <- unique(c(windows[[length(windows) - 1]], windows[[length(windows)]]))
    combined <- combined[order(combined)]
    total_combined <- windows[-c((length(windows) - 1), length(windows))]
    total_combined[[length(total_combined) + 1]] <- combined
    windows <- total_combined
  }
  
  return(windows)
}