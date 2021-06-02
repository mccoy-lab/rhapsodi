#'
#'
#'
sim_generate_donor_hap <- function(num_snps){
  hap1 <- data.frame(V1 = sample(c(0,1), size=num_snps, replace=TRUE)) #simulate first donor chromosome
  hap2 <- abs(1-hap1) #switch bits to construct second donor
  donor_haps <- data.frame(cbind(hap1, hap2)) %>% `colnames<-`(c("donor1", "donor2"))
  return (donor_haps)
}