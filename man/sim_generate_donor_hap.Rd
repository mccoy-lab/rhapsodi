% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_generate_donor_hap.R
\name{sim_generate_donor_hap}
\alias{sim_generate_donor_hap}
\title{This function generates the diploid donor haplotypes of length \code{num_snps} with only heterozygous sites}
\usage{
sim_generate_donor_hap(num_snps)
}
\arguments{
\item{num_snps}{an integer, the number of SNP positions for which we want to generate data}
}
\value{
donor_haps a data frame with the truth phased diploid donor haplotypes with two columns \code{donor1} and \code{donor2}
}
\description{
This function generates the diploid donor haplotypes by first generating a single haplotype
using a uniform distribution to pick the reference (a 0) or the alternate (a 1) allele
for each position for \code{num_snps} positions.
Then for the second donor haplotype we assume that each of these SNP positions are heterozygous
and just invert the bits of the first haplotypes to make this so.
}
