library(BayesianConditionalFPCA)
library(tidyverse)
library(mgcv)
library(spam)

load(paste0("/Users/johnshamshoian/Documents/R_projects/", 
            "BayesianConditionalFPCA/sleep/data/",
            "relative_psd.RData"))

num_epochs <- 480
epoch_grid <- 1:num_epochs
basis_df <- 10 / 100 * num_epochs
epoch_basis <- smooth.construct(s(epoch_grid, bs = "ps", k = basis_df),
                                data.frame(epoch_grid), NULL)$X
penalty <- as.matrix(spam::precmat.RW2(basis_df))
mean_penalty <- list(penalty)
var_penalty <- list(penalty)
mean_indices <- c(1)
var_indices <- c(1)
relative_psd_filtered <- relative_psd %>%
  group_by(id) %>%
  filter(n() >= num_epochs, epoch <= num_epochs,
         id %in% 200001:200100) %>%
  ungroup()

num_subjects <- as.numeric(summarise(relative_psd_filtered, n_distinct(id)))
X <- cbind(rep(1, num_subjects))
Z <- cbind(rep(1, num_subjects))
response <- t(matrix(relative_psd_filtered$psd,
                   nrow = num_epochs,
                   ncol = num_subjects))

k <- 15
iter <- 500
burnin <- 100
nchains <- 1
thin <- 1
loglik <- 0
mcmc_results <- mymain(response, epoch_grid, X,
                               Z, epoch_basis, 
                               mean_penalty, var_penalty,
                               mean_indices, var_indices,
                               k, iter,
                               thin)

posterior_intervals <- get_posterior_predictive_bands2(mcmc_results, c(.025, .5, .975))
colnames(posterior_intervals) <- c("ID", "Time", "Y", "Lower_P", "Median_P", "Upper_P", "Lower_M", "Median_M", "Upper_M")
posterior_intervals <- as_tibble(posterior_intervals)
# posterior_intervals %>%
#   filter(ID == sub) %>%
#   ggplot(aes(x = Time, y = Y_no_error)) +
#   geom_point(na.rm = TRUE) +
#   geom_ribbon(aes(ymin = Lower_M, ymax = Upper_M), alpha = 0.3) +
#   theme_bw()
posterior_intervals %>%
  filter(ID == 10) %>%
  ggplot(aes(x = Time, y = Y)) +
  geom_point(na.rm = TRUE) +
  geom_line(aes(y = Median_M)) +
  geom_ribbon(aes(ymin = Lower_M, ymax = Upper_M), alpha = 0.3) +
  theme_bw()

xi <- c(1)
alpha <- .05
means <- get_posterior_means(mcmc_results, xi, alpha)
plot(means[,2], type = "l")
lines(means[,1])
lines(means[,3])
