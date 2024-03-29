% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/export_data.R
\name{export_data}
\alias{export_data}
\title{A function to export the data from rhapsodi}
\usage{
export_data(
  input_data,
  complete_haplotypes,
  filled_gametes,
  recomb_breaks,
  acgt = FALSE
)
}
\arguments{
\item{input_data}{the named list from \code{read_data}}

\item{complete_haplotypes}{the dataframe from \code{phase_donor_haplotypes}}

\item{filled_gametes}{the named list from \code{impute_gamete_genotypes}}

\item{recomb_breaks}{the dataframe from \code{discover_meitoic_recombination}}

\item{acgt}{a bool; default is FALSE. If TRUE, assumes that the input data was not 0/1/NA encoded, rather gamete genotypes were A/C/G/T/NA encoded and the dataframe had ref and alt columns. Will add these columns back as a0 and a1 for 0/1 encoded output data}
}
\value{
rhapsodi_out a named list with \code{donor_haps}, \code{gamete_haps}, \code{gamete_genotypes}, \code{unsmoothed_gamete_haps},  \code{unsmoothed_gamete_genotypes}, and \code{recomb_breaks}
}
\description{
This function exports the data from rhapsodi by first adding index and pos columns to the gamete output data.
Then, continues by checking whether the input data was originally nucleotide encoded (A/C/G/T/NA).
If so, the function appends a0 and a1 columns to 0/1 encoded outputs by calling \code{record_ref_alt}, such that the nucleotide(s) associated with genotype 0 are recorded in column a0
and the nucleotide(s) associated with genotype 1 are recorded in column a1
Finally, the function returns a named list which has
\code{donor_haps} which is the phased haplotypes as a data frame with column names index, pos (for SNP positions), h1 (haplotype 1), & h2 (haplotype 2) if acgt = FALSE. Otherwise: index, pos, a0, a1, h1, h2
\code{gamete_haps} which is the filled gamete data frame specifying from which donor haplotype each gamete position originates. Column names: index, pos, gamete_names.
\code{gamete_genotypes} which is the filled gamete dataf rame specifying the genotype (in 0's and 1's) for each gamete position. If acgt = FALSE, column names: index, pos, gamete_names. Otherwise: index, pos, a0, a1, gamete_names
\code{unsmoothed_gamete_haps} which is the filled gamete data frame specifying from which donor haplotype each gamete position originates in data frame form, after unsmoothing the data by replacing imputed values with original sequencing reads when there's disagreement between observations and imputation. Column names: index, pos, gamete_names.
\code{unsmoothed_gamete_genotypes} which is the filled gamete data frame specifying the genotype (in 0's and 1's) for each gamete position, after unsmoothing the the data by replacing imputed values with original sequencing reads when there's disagreement between observations and imputation. If acgt = FALSE, column names: index, pos, gamete_names. Otherwise: index, pos, a0, a1, gamete_names
\code{recomb_breaks}which is a data frame specifying the recombination breakpoints for each gamete. Column names: Ident, Genomic_start, Genomic_end
}
