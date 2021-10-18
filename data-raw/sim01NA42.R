library(rhapsodi)
library(here)



generated_data <- rhapsodi::sim_run_generative_model(50, 5000, 0.1, 1)
sim01NA42 <- generated_data$gam_na

save(sim01NA42, file="../data/sim01NA42.rda")
write.csv(sim01NA42, file="simulated_01NA_encoded_rs42.csv")