% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_mode.R
\name{get_mode}
\alias{get_mode}
\title{A function to find the more common allele or NA at each SNP}
\usage{
get_mode(vector)
}
\arguments{
\item{vector}{A subset of donor haplotypes}
}
\value{
mode The most frequent value or NA at a position
}
\description{
This function gets the mode of a vector for majority voting or returns NA if there is no single mode
}
\details{
Adapted from  #from https://stackoverflow.com/questions/56552709/r-no-mode-and-exclude-na?noredirect=1#comment99692066_56552709
}
