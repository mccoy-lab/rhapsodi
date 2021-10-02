#' Build a model to assign haplotypes to each gamete
#' 
#' This function builds a hidden Markov model that considers sequencing error. Uses the 
#' HMM packages. 
#' 
#' @param num_snps number of SNPs, found previously from taking the number of rows in the input data. Used to set the denominator for transition probability
#' @param sequencing_error User-input for expected error in sequencing (default = 0.005)
#' @param avg_recomb User-input for expected average recombination spots per chromosome (default=1). Used to set the numerator for transition probability
#' 
#' @importFrom HMM initHMM
#' 
#' @return hmm The hidden Markov model with transition and emission probabilities set for use. 
#' 
build_hmm <- function(num_snps, sequencing_error, avg_recomb) {
  # two states
  states <- c("h1", "h2")
  
  # probability of state at position x+1 given state at position x   
  hap1Prob <- c(1-(avg_recomb/num_snps), avg_recomb/num_snps)
  hap2Prob <- rev(hap1Prob)
  # Probability of transitioning at any given position 
  transProb <- matrix(c(hap1Prob, hap2Prob), 2)
  
  # Two emissions (observations): an allele from h1 or an allele from h2
  emissions <- c("hap1", "hap2")
  
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