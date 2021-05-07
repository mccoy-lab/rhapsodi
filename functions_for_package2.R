# Formatting and commenting functions for `write and upload `rhapsodi` package


# Questions: 
# In `split_with_overlap` 
# How do we communicate that `positions` is something that needs to be globally defined? 



### Part 1: Impute parental haplotypes `impute_parental_haplotypes`
### Internal functions: `get_mode`, `getmode`, invert_bits`, `invertbits`, `split_with_overlap`, `reconstruct_haplotypes`


#' A function to find the more common allele (i.e., 0 or 1) at each SNP 
#'
#' This function gets the mode of a vector after removing the NAs. Explain why we use this 
#'  
#' @param vector A subset of 
#' 
#' @return mode The most frequent value at a position 
#' 
#' @example 
#' R code here showing how my function works  
#' 
get_mode <- function(vector) {
  uniqv <- unique(vector)
  uniqv <- uniqv[!is.na(uniqv)]
  mode <- uniqv[which.max(tabulate(match(vector, uniqv)))]
  return (mode)
}

#' A function to find the more common allele or NA at each SNP
#' 
#' This function gets the mode of a vector for majority voting or returns NA if there is no single mode
#' 
#' @param vector A subset of
#' 
#' @return mode The most frequent value or NA at a position
#' 
#' @example 
#' R code here showing how my function works
#' 
getmode <- function(vector){
  uniqv <- unique(na.omit(vector))
  tabv <- tabulate(match(vector, uniqv))
  if (length(uniqv) != 1 & sum(max(tabv) == tabv) > 1){
    if (is.character(uniqv)) return(NA_character_) else return(NA_real_)
  }
  max_tabv <- tabv == max(tabv)
  return(uniqv[max_tabv])
}


#' A function to invert the values in a data frame
#' 
#' This function replaces 0s with 1s and 1s with 0s in a dataframe. This inverted dataframe is an 
#' input to the reconstruction of parental haplotypes. 
#' 
#' @param input_matrix Matrix of gamete alleles
#' 
#' @return input_matrix Inverted matrix of gamete alleles
#' 
#' @examples 
#' R code here showing my function works  
#' 
invert_bits <- function(input_matrix) {
  input_matrix[input_matrix == 0] <- -1
  input_matrix[input_matrix == 1]  <- 0
  input_matrix[input_matrix == -1] <- 1
  return(input_matrix)
}


#' A function to invert the values in a data frame or matrix, faster than invert_bits
#' 
#' This function replaces 0s with 1s and 1s with 0s in a dataframe or matrix
#' 
#' @param input_data
#' 
#' @return input_data inverted from the actual input
#' 
#' @examples
#' R code here showing how my function works
#' 
invertbits <- function(input_data) {
  return(abs(input_data-1))
}

#' A function to separate the gametes into overlapping segments 
#' 
#' This function separates each gamete into overlapping windows of size determined by the user. These 
#' overlaps are an input to the distance matrix constructed in `reconstruct_haplotypes`. 
#'
#' @param vector_of_positions # print positions and see what it is 
#' @param window_size User-input value for length of segment (default = 2500)
#' @param overlap Degree of overlap between segments (default = window_size/2)
#'
#' @return ? check if need to 
#' 
#' @examples 
#' R code here showing my function works 
#'
# do we need to cite this? 
# overlapping window function from https://stackoverflow.com/questions/8872376/split-vector-with-overlapping-samples-in-r
split_with_overlap <- function(vector_of_positions, window_size, overlap) {
  starts = seq(1, length(vector_of_positions), by = window_size - overlap)
  ends   = starts + window_size - 1
  ends[ends > length(vector_of_positions)] = length(vector_of_positions)
  lapply(1:length(starts), function(i) vector_of_positions[starts[i]:ends[i]])
}


#' A function to reconstruct the parental haplotypes using data from gametes
#' 
#' This function takes gamete data and computes and clusters a distance matrix with a binary method. It 
#' anticipates sparse data and replaces NAs with 0.5. It clusters the tree into two groups, i.e., haplotypes. 
#' It then categorizes the sperm cells falling into each of the two groups. It reconstructs the original 
#' haplotypes by majority vote after inverting the opposite haplotype. 
#' 
#' @param input_dt A dataframe of sparse gamete data coded with reference and alt alleles (i.e., 0, 1, or NA)
#' @param input_positions #print positions to see what this is 
#' @param window_indices Overlapping segments of the gamete data (using segments generated from `split_with_overlap`) 
#' @param threads Number of threads to use for computation (default = 1)
#' 
#' @return inferred_output A tibble with the inferred haplotype sorted by position and window
#' 
#' @example 
#' R code here showing my function works 
#' 
reconstruct_haplotypes <- function(input_dt, input_positions, window_indices, threads = 1) { #why do we need the threads param here?
  window_start <- min(window_indices)
  window_end <- max(window_indices)
  positions_for_window <- input_positions[window_start:window_end]
  # compute a distance matrix
  d <- dist(t(as.matrix(input_dt)[window_start:window_end,]), method = "binary")
  # plug in 0.5 for any NA entries of the distance matrix
  d[is.na(d)] <- 0.5
  # cluster the distance matrix
  tree <- hclust(d, method = "ward.D2")
  # plot(tree, cex = 0.1) # uncomment to plot
  # cut the tree generated by clustering into two groups (haplotypes)
  haplotypes <- cutree(tree, k=2)
  # get the names of the sperm cells falling into the two groups
  h1_sperm <- names(haplotypes[haplotypes == 1])
  h2_sperm <- names(haplotypes[haplotypes == 2])
  # reconstruct the original haplotypes by majority vote after inverting the opposite haplotype
  h1_inferred <- unname(apply(cbind(input_dt[window_start:window_end, h1_sperm],
                                    invertbits(input_dt[window_start:window_end, h2_sperm])),
                              1, function(x) getmode(x)))
  h2_inferred <- unname(apply(cbind(input_dt[window_start:window_end, h2_sperm],
                                    invertbits(input_dt[window_start:window_end, h1_sperm])),
                              1, function(x) getmode(x)))
  inferred_output <- tibble(index = window_indices, pos = positions_for_window, h1 = h1_inferred)
  return(inferred_output)
}


#' A function to impute parental haplotypes 
#' 
#' This function 
#' 
#' @export
#' 
#' @param dt Matrix of gamete alleles 
#' @param window_length Size of window 
#' @param positions Vector of  
#' @param threads 
#' 
#' @return complete_haplotypes
#' 
#' @example 
#' R code here showing my function works 
#' 
impute_parental_haplotypes <- function(dt, window_length, positions, threads) {
  windows <- split_with_overlap(rank(positions), window_length, overlap = window_length / 2)
  
  if (length(windows) > 1){ #merge the last two windows to avoid edge effect
    combined <- unique(c(windows[[length(windows) - 1]], windows[[length(windows)]]))
    combined <- combined[order(combined)]
    total_combined <- windows[-c((length(windows) - 1), length(windows))]
    total_combined[[length(total_combined) + 1]] <- combined
    windows <- total_combined
  }
  
  # infer the haplotypes within the overlapping windows
  inferred_haplotypes <- pbmclapply(1:length(windows), 
                                    function(x) reconstruct_haplotypes(dt, positions, windows[[x]]),
                                    mc.cores = getOption("mc.cores", threads))
  # stitch together the haplotypes
  initial_haplotype <- inferred_haplotypes[[1]]
  for (hap_window in 1:length(windows)) {
    olap_haps <- merge(initial_haplotype, inferred_haplotypes[[hap_window]], by = "index")
    olap_haps_complete <- merge(initial_haplotype, inferred_haplotypes[[hap_window]], by = "index", all = TRUE)
    mean_concordance <- mean(olap_haps$h1.x == olap_haps$h1.y, na.rm=TRUE)
    if (mean_concordance < 0.1) {
      olap_haps_complete$h1.y <- invertbits(olap_haps_complete$h1.y)
    } else if (mean_concordance < 0.9) {
      error(paste0("Haplotypes within overlapping windows are too discordant to merge. Mean: ", mean_concordance))
    }
    initial_haplotype <- tibble(index = olap_haps_complete$index,
                                pos = c(olap_haps_complete[is.na(olap_haps_complete$pos.y),]$pos.x,
                                        olap_haps_complete[!is.na(olap_haps_complete$pos.x) &
                                                             !is.na(olap_haps_complete$pos.y),]$pos.x,
                                        olap_haps_complete[is.na(olap_haps_complete$pos.x),]$pos.y),
                                h1 = c(olap_haps_complete[is.na(olap_haps_complete$pos.y),]$h1.x,
                                       olap_haps_complete[!is.na(olap_haps_complete$pos.x) &
                                                            !is.na(olap_haps_complete$pos.y),]$h1.x,
                                       olap_haps_complete[is.na(olap_haps_complete$pos.x),]$h1.y))
  }
  complete_haplotypes <- initial_haplotype %>%
    mutate(h2 = invertbits(h1))
  return(complete_haplotypes)
}

# Example to call `complete_haplotypes`
# complete_haplotypes <- impute_parental_haplotypes(dt, window_length=2500, positions=input_datatable[, 1], threads=1)

### Part 2: Fill gametes with assignments to one haplotype or the other `fill_gametes`
### Internal functions: `recode_gametes`, `recodegametes`, `build_hmm`, `run_hmm`, `fill_na`

#' A function to assign each allele to a haplotype
#' 
#' This function reads along each gamete and replaces its read (i.e., 0 or 1) with the corresponding
#' haplotype (h1 or h2) based on the allele in each parental haplotype. 
#' 
#' @param dt Input matrix of gametes
#' @param complete_haplotypes Inferred parental haplotypes 
#' 
#' @return dt Matrix of gametes coded by haplotype at each position 
#' 
#' @example 
#' R code here showing my function works 
#' 
recode_gametes <- function(dt, complete_haplotypes) {
  # Going through each gamete, if an allele (0 or 1) in a gamete matches the allele (0 or 1)
  # in h1 at that position, replace the allele with "h1". Do the same for h2.
  for (i in 1:ncol(dt)) {
    dt[i][dt[i] == complete_haplotypes$h1] <- "h1"
    dt[i][dt[i] == complete_haplotypes$h2] <- "h2"
  }
  return(dt)
}


#' A function to assign each allele to a haplotype or an NA if not enough information is known
#' 
#' This function reads along each gametes and replaces its read (0 or 1) with the corresponding
#' haplotype (h1 or h2) based on the allele in each parental haplotype; alternatively, the read 
#' may be replaced with an NA if not enough information was known for phasing
#' 
#' @param dt Input matrix of gamete reads
#' @param complete_haplotypes Inferred parental haplotypes
#' 
#' @return dt matrix of gametes coded by haplotype at each position 
#' 
#' @example
#' R code here showing how my function works
#' 
recodegametes <- function(dt, complete_haplotypes) {
 for (i in 1:ncol(dt)) {
   dt[i][dt[i] == complete_haplotypes$h1] <- "h1"
   dt[i][dt[i] == complete_haplotypes$h2] <- "h2"
   dt[c(which(dt[,i] == 0 | dt[,i] == 1)),i] <- NA
 }
  return(dt)
}


#' Build a model to assign haplotypes to each gamete
#' 
#' This function builds a hidden Markov model that considers sequencing error. Uses the 
#' HMM packages. 
#' 
#' @param dt Input matrix
#' @param sequencing_error User-input for expected error in sequencing (default = 0.005)
#' @param avg_recomb User-input for expected average recombination spots per chromosome
#' 
#' @return hmm The hidden Markov model with transition and emission probabilities set for use. 
#' 
#' @example 
#' R code here showing my function works 
build_hmm <- function(dt, sequencing_error, avg_recomb) {
  # set denominator for transition probability - one recombination event per chromosome
  num_snps <- nrow(dt)
  
  # two states
  states <- c("haplotype1", "haplotype2")
  
  # probability of state at position x+1 given state at position x   
  hap1Prob <- c(1-(avg_recomb/num_snps), avg_recomb/num_snps)
  hap2Prob <- rev(hap1Prob)
  # Probability of transitioning at any given position 
  transProb <- matrix(c(hap1Prob, hap2Prob), 2)
  
  # Two emissions (observations): an allele from h1 or an allele from h2
  emissions <- c("h1","h2")
  
  # Prob of emitting an h1 allele, prob of emitting an h2 allele in state `haplotype1`
  h1ProbEmiss <- c((1-sequencing_error), sequencing_error)
  # Prob of emitting an h1 allele, prob of emitting an h2 allele in state `haplotype2`
  h2ProbEmiss <- rev(h1ProbEmiss)
  emissProb <- matrix(c(h1ProbEmiss, h2ProbEmiss), 2)
  
  #build model with the above inputs
  hmm <- initHMM(States = states,
                 Symbols = emissions,
                 transProbs = transProb,
                 emissionProbs = emissProb)
  return(hmm)
}


#' A function that assigns the position in each gamete to the correct parental haplotype
#' 
#' This function walks along each gamete and replaces the original observation with the inferred state.
#' This replaces any incorrect haplotype assignments (e.g., due to sequencing error).
#'
#' @param gamete_dt Matrix of gametes with 0 and 1
#' @param column_index `Apply` function cycles through gamete_dt to act on each column (i.e., each gamete) 
#' @param complete_haplotypes Inferred parental haplotypes 
#' @param sequencing_error User-input for expected error in sequencing (default = 0.005) 
#' 
#' @return original_obs Replaced observed haplotype with that assigned by the model
#' 
#' @example 
#' R code here showing my function works 
#' 
run_hmm <- function(dt, column_index, complete_haplotypes, sequencing_error, avg_recomb) {
  # build the hmm 
  hmm <- build_hmm(dt, sequencing_error, avg_recomb)
  
  original_obs <- dt[,column_index]
  inferred_state <- viterbi(hmm, na.omit(dt[, column_index]))
  original_obs[!is.na(original_obs)] <- inferred_state
  return(original_obs)
}


# add boolean function to superimpose original 

#' A function to fill in missing data from each gamete 
#' 
#' This function fills in missing data (NAs) on each gamete. For each gamete, it fills the NA values with the nearest haplotype.
#' If the two adjacent haplotypes are not the same (i.e., at a recombination breakpoint), it leaves the values as NA. 
#' It offers the option to avoid oversmoothing by superimposing initial haplotype assignments over each gamete. For example, if an 
#' allele assigned to h1 was changed by the model to h2, this function can fill the NAs to h2, but replace the singular h1
#' at the correct allele. This could be an example of gene conversion or non-crossover. 
#' 
#' @param imputed_gametes Output of `run_hmm` which assigned a parental haplotype to each segment of each gamete
#' @param col_index Each column of `imputed_gametes`, pulled via `apply` function 
#' 
#' @return gamete_sample_imputed Column with each gamete's imputed haplotypes 
#' 
#' @example 
#' R code here showing my function works 
#' 
fill_na <- function(imputed_gametes, col_index) {
  gamete_sample <- imputed_gametes[,col_index] %>%
    rename(gamete = colnames(.)[1]) %>%
    mutate(gamete_up = gamete) %>%
    mutate(gamete_down = gamete) %>%
    fill(gamete_up, .direction = "up") %>%
    fill(gamete_down, .direction = "down") %>%
    mutate(is_match = (gamete_up == gamete_down)) %>%
    replace_na(list(is_match = FALSE))
  gamete_sample$gamete_imputed <- as.character(NA)
  gamete_sample[gamete_sample$is_match == TRUE,]$gamete_imputed <- gamete_sample[gamete_sample$is_match == TRUE,]$gamete_up
  #fill beginning of chromosome NAs
  first <- which(!is.na(gamete_sample$gamete_imputed))[1]
  gamete_sample$gamete_imputed[1:(first-1)] <- gamete_sample$gamete_imputed[first]
  #fill end of chromosome NAs
  gamete_sample$gamete_imputed <- rev(gamete_sample$gamete_imputed)
  first <- which(!is.na(gamete_sample$gamete_imputed))[1]
  gamete_sample$gamete_imputed[1:(first-1)] <- gamete_sample$gamete_imputed[first]
  #reverse chromosome imputation back so it faces the right way
  gamete_sample$gamete_imputed <- rev(gamete_sample$gamete_imputed)
  gamete_sample_imputed <- gamete_sample$gamete_imputed
  return(gamete_sample_imputed)
}


#' A function to assign the haplotypes of each allele on every gamete
#' 
#' This function builds and applies a hidden Markov model to categorize each allele on each gamete. 
#' It then fills the positions missing data with the nearest haplotype assignment. It offers the option to 
#' superimpose the original haplotype. 
#' 
#' @param dt matrix of gametes
#' @param complete_haplotypes Inferred parental haplotypes 
#' @param sequencing_error User-input for expected error in sequencing (default = 0.005) 
#' 
#' @return filled_gametes matrix with same dimensions as dt but with correctly assigned haplotype
#' and NAs filled 
#' 
#' @export
#' 
#' @example 
#' R code here showing my function works 
#' 

# add option for boolean 
fill_gametes <- function(dt, complete_haplotypes, sequencing_error=0.005, threads) { #test to see if pbmclapply complains if a single thread is used
  dt_recoded <- recodegametes(dt, complete_haplotypes)
  
  imputed_gametes <- as_tibble(do.call(cbind, pbmclapply(1:ncol(dt_recoded),
                                                         function(x) run_hmm(dt_recoded, x),
                                                         mc.cores = getOption("mc.cores", threads))))

  filled_gametes <- as_tibble(do.call(cbind, 
                                      pblapply(1:ncol(imputed_gametes),
                                               function(x) fill_na(imputed_gametes, x),
                                               mc.cores = getOption("mc.cores", threads))))
  colnames(filled_gametes) <- colnames(dt)
  return(filled_gametes)
}

# Sample for calling `filled_gametes`
# outcome <- fill_gametes(dt, complete_haplotypes, sequencing_error=0.005, threads=1)


### Part 3: Find recombination spots and write out results (possibly with real reads instead of smoothed reads) `report_gametes`
### Internal functions: `unsmooth`, `find_recomb_spots`, `re_recode_gametes`, `report_gametes`

#' A function to unsmooth or replace original reads for the gametes
#' 
#' This function finds where the resulting haplotype assignments differ from the original reads and replaces the imputed data with
#' the originally observed data, hence unsmoothing the HMM signal.
#' 
#' @param original_gamete_data original gamete data with haplotype by position 
#' @param filled_gamete_data filled gamete data from `fill_gametes`
#'
#' @return filled_gamete_data filled gamete data in tibble form with replaced original reads 
#' 
#' @example
#' R code here explaining what my function does
#' 
unsmooth <- function(original_gamete_df, filled_gamete_data) {
  original_gamete_df[original_gamete_df == "h1"] <- "haplotype1"
  original_gamete_df[original_gamete_df == "h2"] <- "haplotype2"
  original_dt <- as.data.frame(original_gamete_df)
  filled_gamete_data <- as.data.frame(filled_gamete_data)
  filled_gamete_data[!is.na(original_dt)] <- original_dt[!is.na(original_dt)]
  filled_gamete_data <- as_tibble(filled_gamete_data)
  return (filled_gamete_data)
}

#' A function to find recombination breakpoint regions for a single gamete
#' 
#' This function finds recombination breakpoint regions for a single gamete by finding 
#' where adjacent non-NA haplotypes switch. Most likely large regions will be returned 
#' rather than 2bp windows because the sparsity of the original data constrains how 
#' many NAs will remain after filling and therefore we can only say that the true
#' exchange point occurs somewhere between these two non-matching haplotypes
#' 
#' @param input_gamete_data (tibble) SNPs correspond to the rows and gametes correspond to the columns
#' @param x which gamete number to search
#' @param identities gamete identities vector
#' @param genomic_positions Genomic SNP positions 
#' 
#' @return recomb_spots tibble of the recombination spots for gamete x
#' 
#' @example 
#' R code here explaining what my function does
#' 
find_recomb_spots <- function(input_gamete_data, x, identities, genomic_positions) {
  ident <- identities[x]
  input_tibble <- input_gamete_data[, x] %>%
    mutate(., index = row_number()) %>%
    mutate(., positions = genomic_positions)
  complete_cases_tibble <- input_tibble[complete.cases(input_tibble),]
  input_vec <- as.factor(complete_cases_tibble[[1]])
  switch_indices <- which(input_vec[-1] != input_vec[-length(input_vec)])
  switch_indices_input <- complete_cases_tibble[switch_indices,]$index
  crossover_start <- input_tibble[switch_indices_input,]$positions
  rev_input_tibble <- arrange(input_tibble, -index) %>%
    mutate(., index = row_number())
  complete_cases_rev_tibble <- rev_input_tibble[complete.cases(rev_input_tibble),]
  rev_input_vec <- as.factor(complete_cases_rev_tibble[[1]])
  rev_switch_indices <- which(rev_input_vec[-1] != rev_input_vec[-length(rev_input_vec)])
  rev_switch_indices_input <- complete_cases_rev_tibble[rev_switch_indices,]$index
  crossover_end <- rev(rev_input_tibble[rev_switch_indices_input,]$positions)
  recomb_spots <- tibble(Ident = ident, Genomic_start = crossover_start, Genomic_end = crossover_end)
  return (recomb_spots)
}

#' A function to turn a matrix of gamete haplotypes by position back to the 0/1 reads by
#' position, using the phased parental haplotypes as a guide
#' 
#' This function builds a matrix with position by row and gametes by column such that each cell
#' is a 0 or a 1 or an NA based on whether that cell in the input gamete matrix was from haplotype1, haplotype2, 
#' or was an NA. Then the 0 or 1 is found in the complete_haplotypes (phased parentals) input at the corresponding
#' positions 
#'
#' @param dt input gamete haplotype data in tibble form
#' @param complete_haplotypes dataframe of phased parentals with two columns (h1 and h2) and SNP positions as rows 
#'
#' @return to_return gamete data in dataframe form with read data (0's and 1's and NAs) instead of haplotype information
#' 
#' @example
#' R code here explaining what my function does
#'
re_recode_gametes <- function(dt, complete_haplotypes) {
  to_return <- data.frame(matrix(NA_real_, nrow=nrow(dt), ncol=ncol(dt)))
  for (i in 1:ncol(dt)) {
    locs_h1 <- dt[,i] == "haplotype1"
    locs_h1[which(is.na(locs_h1))] <- FALSE
    locs_h2 <- dt[,i] == "haplotype2"
    locs_h2[which(is.na(locs_h2))] <- FALSE
    to_return[locs_h1, i] <- complete_haplotypes$h1[locs_h1]
    to_return[locs_h2, i] <- complete_haplotypes$h2[locs_h2]
  }
  return (to_return)
}

#' A function to drive and report the unsmoothing (if desired) and recombination finding
#' 
#' This function takes as input two booleans controlling whether the reported genotypes are smoothed
#' or unsmoothed (replacing inferred HMM state with the original reads if they disagree) 
#' and whether the data given to recombination finding is smoothed or unsmoothed. Then the function
#' runs recombination finding (and eventually reporting/exporting of the cool data)
#' 
#' @param smooth_crossovers boolean whether to use smoothed data for recombination finding. If TRUE, doesn't replace with original reads
#' @param smooth_imputed_genotypes boolean whether to use smoothed data for ending genotypes. If TRUE, doesn't replace with original reads
#' @param complete_haplotypes dataframe of phased parental genotypes in two columns, one for each parental haplotype
#' @param original_gamete_data 
#' 
#'
#' @export
#' 
#' @example 
#' R code showing how my function works
#'
report_gametes <- function(smooth_crossovers, smooth_imputed_genotypes, complete_haplotypes, original_gamete_data, filled_gamete_data, sampleName, chrom, positions) {
  if (!smooth_crossovers & !smooth_imputed_genotyptes) {
    filled_gamete_forrecomb <- unsmooth(original_gamete_data, filled_gamete_data) #haplotypes
    idents_for_csv <- paste0(paste0(sampleName, "_", chrom, "_"), colnames(filled_gamete_forrecomb))
    recomb_spots_all <- do.call(rbind, pbmclapply(1:ncol(filled_gamete_forrecomb),
                                                  function(x) find_recomb_spots(filled_gamete_forrecomb, x, idents_for_csv, positions),
                                                  mc.cores=getOption("mc.cores", threads))) %>%
      right_join(., tibble(Ident = idents_for_csv), by = "Ident")
    filled_gamete_recode <- re_recode_gametes(filled_gamete_forrecomb, complete_haplotypes) #0's and 1's
  } else if (!smooth_crossovers & smooth_imputed_genotypes) {
    filled_gamete_forrecomb <- unsmooth(original_gamete_data, filled_gamete_data)
    idents_for_csv <- paste0(paste0(sampleName, "_", chrom, "_"), colnames(filled_gamete_forrecomb))
    recomb_spots_all <- do.call(rbind, pbmclapply(1:ncol(filled_gamete_forrecomb),
                                                  function(x) find_recomb_spots(filled_gamete_forrecomb, x, idents_for_csv, positions),
                                                  mc.cores=getOption("mc.cores", threads))) %>%
      right_join(., tibble(Ident = idents_for_csv), by = "Ident")
    filled_gamete_recode <- re_recode_gametes(filled_gamete_data, complete_haplotypes) #0's and 1's
} else if (smooth_crossovers & !smooth_imputed_genotypes){
    idents_for_csv <- paste0(paste0(samplenName, "_", chrom, "_"), colnames(filled_gamete_data))
    recomb_spots_all <- do.call(rbind, pbmclapply(1:ncol(filled_gamete_data),
                                                  function(x) find_recomb_spots(filled_gamete_data, x, idents_for_csv, positions),
                                                  mc.cores=getOption("mc.cores", threads))) %>%
      right_join(., tibble(Ident = idents_for_csv), by = "Ident")
    filled_gamete_data <- unsmooth(original_gamete_data, filled_gamete_data) #haplotypes
    filled_gamete_recode <- re_recode_gamtes(filled_gamete_data, complete_haplotypes) #0's and 1's
    
} else { #both are TRUE
  idents_for_csv <- paste0(paste0(samplenName, "_", chrom, "_"), colnames(filled_gamete_data))
  recomb_spots_all <- do.call(rbind, pbmclapply(1:ncol(filled_gamete_data),
                                                function(x) find_recomb_spots(filled_gamete_data, x, idents_for_csv, positions),
                                                mc.cores=getOption("mc.cores", threads))) %>%
    right_join(., tibble(Ident = idents_for_csv), by = "Ident")
  filled_gamete_recode <- re_recode_gamtes(filled_gamete_data, complete_haplotypes) #0's and 1's
  }
}

### Part 4: Functions for plotting


### Part 5: Functions for simulations 


### Part 6: Functions for testing 
















