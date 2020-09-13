k <- commandArgs(traulingONLY = TRUE)
infile <- paste0("/Users/johnshamshoian/Documents/R_projects",
                 "bfcr/sleep/mcmc_output/mcmc_results",
                 k,
                 ".rds")
mcmc_results <- readRDS(infile)
info_criteria <- calculate_waic(mcmc_results)
outfile <- paste0("/Users/johnshamshoian/Documents/R_projects",
                  "bfcr/sleep/info_criteria",
                  k,
                  ".rds")
saveRDS(info_criteria, file = outfile)