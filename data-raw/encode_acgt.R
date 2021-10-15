library(tidyverse)
library(here)
library(pbmcapply)

#args <- commandArgs(trailingOnly = TRUE)
#threads <- as.integer(args[1])
threads <- 2L
load("../data/sim01NA42.rda")

set.seed(42)

num_snps <- nrow(sim01NA42)
positions <- sim01NA42$positions
data01 <- sim01NA42[,-1]
possible_alleles <- c("A", "C", "G", "T", "AC", "TG", "GAA", "CAT", "CCCC")
alt_alleles_matched <- c("G", "T", "A", "C", "GT", "CC", "CCT", "TTC", "GAGT")
allele_locs <- sample(1:length(possible_alleles), num_snps, replace = TRUE)
ref <- possible_alleles[allele_locs]
alt <- alt_alleles_matched[allele_locs]

acgt_recode_gametes <- function(dt_col, num_snps, ref, alt){
  locs_ref <- which(dt_col == 0)
  locs_alt <- which(dt_col == 1)
  to_return <- rep(NA, num_snps)
  to_return[locs_ref] <-  ref[locs_ref]
  to_return[locs_alt] <- alt[locs_alt]
  return (to_return)
}

dataACGT <- cbind(positions, ref, alt, do.call(cbind, 
                                       pbmcapply::pbmclapply(1:ncol(data01), function(x) acgt_recode_gametes(data01[,x], num_snps, ref, alt),
                                       mc.cores = getOption("mc.cores", threads)))) %>%  as.data.frame() %>% `colnames<-`(c("positions", "ref", "alt", colnames(data01)))

save(dataACGT, file="../data/dataACGT.rda")

