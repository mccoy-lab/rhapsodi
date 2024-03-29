% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fill_na.R
\name{fill_na}
\alias{fill_na}
\title{A function to fill in missing data from each gamete}
\usage{
fill_na(imputed_gametes, col_index, fill_ends = TRUE)
}
\arguments{
\item{imputed_gametes}{Output of \code{run_hmm} which assigned a donor haplotype to each segment of each gamete}

\item{col_index}{Each column of \code{imputed_gametes}, pulled via \code{apply} function}

\item{fill_ends}{a boolean; if TRUE, fills the NAs at the terminal edges of chromosomes with the last known or imputed SNP (for end of chromosome) and the first known or imputed SNP (for beginning of chromosome); if FALSE, leaves these genotypes as NA; default = TRUE}
}
\value{
gamete_sample_imputed Column with each gamete's imputed haplotypes
}
\description{
This function fills in missing data (NAs) on each gamete. For each gamete, it fills the NA values with the nearest haplotype.
If the two adjacent haplotypes are not the same (i.e., at a recombination breakpoint), it leaves the values as NA.
}
