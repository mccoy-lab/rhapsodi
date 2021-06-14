#'
#'
#'
standard_input <- function(input_file){
  stopifnot(file_test("-f", input_file), paste0("The provided file ", input_file, " does not exist."))
  dt <- read_delim(input_file, delim = "\t", col_types = cols(.default = "d")) %>% 
    as.data.frame()
  
  dt <- dt[((rowSums(dt[,2:ncol(dt)] == 0, na.rm = TRUE) > 0) & (rowSums(dt[, 2:ncol(dt)] == 1, na.rm=TRUE) > 0)),]
  
  positions <- dt[,1]
  return(list(dt=dt, positions=positions)) 
}