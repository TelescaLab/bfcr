library(BayesianConditionalFPCA)
library(tidyverse)
library(mgcv)
library(spam)
library(microbenchmark)

load(paste0("/Users/johnshamshoian/Documents/R_projects/", 
            "BayesianConditionalFPCA/sleep/data/",
            "relative_psd.RData"))
sleep_tabulated <- read.csv(paste0("/Users/johnshamshoian/Documents/R_projects/",
            "BayesianConditionalFPCA/sleep/tabulated_data/",
            "shhs1-dataset-0.15.0.csv"))

sleep_tabulated_filtered <- sleep_tabulated %>%
  select(nsrrid, age_s1)
sleep_data <- inner_join(sleep_tabulated_filtered, relative_psd)
sleep_data_filtered <- sleep_data %>%
  group_by(nsrrid) %>%
  filter(n() >= num_epochs, epoch <= num_epochs, 
         nsrrid %in% 200001:200100) %>%
  ungroup()

num_epochs <- 480
epoch_grid <- 1:num_epochs
epoch_df <- ceiling(10 / 100 * num_epochs)
age_grid <- sleep_data_filtered %>%
  group_by(nsrrid) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(age_s1) 
age_df <- ceiling(10 / 100 * (max(age_grid) - min(age_grid)))


epoch_basis <- smooth.construct(s(epoch_grid, bs = "ps", k = epoch_df),
                                data.frame(epoch_grid), NULL)$X
spec_age <- s(age_grid, bs = "ps", k = age_df)
age_dat <- data.frame(age_grid)
spec_age$term <- "age_s1"
age_basis <- smoothCon(object = spec_age,
                       data = data.frame(age_grid))[[1]]$X
penalty <- as.matrix(spam::precmat.RW2(epoch_df))
# penalty <- as.matrix(precmat.IGMRFreglat(age_df, epoch_df), order = 2)
mean_penalty <- list(penalty)
var_penalty <- list(penalty)
mean_indices <- c(1)
var_indices <- c(1)

num_subjects <- as.numeric(summarise(sleep_data_filtered, n_distinct(nsrrid)))
design_mean <- cbind(rep(1, num_subjects))
design_var <- cbind(rep(1, num_subjects))
# design_mean <- age_basis
# design_var <- age_basis
response <- t(matrix(sleep_data_filtered$psd,
                   nrow = num_epochs,
                   ncol = num_subjects))

k <- 15
iter <- 10000
burnin <- 100
nchains <- 1
thin <- 1
loglik <- 0
results <- mymain(response, design_mean,
                  design_var, epoch_basis, 
                  mean_penalty, var_penalty,
                  mean_indices, var_indices,
                  k, iter, thin = 1, var = "unequal")

iter_2 <- 10000
subj <- 70
interval <- 300:480
plot(response[subj,][interval])
for(iter_s in 7500:iter_2) {
  myfit <- epoch_basis %*% results$beta[,,iter] %*% design_mean[subj,]
  # fit_beta <- epoch_basis %*% results$fit_beta
  for(kp in 1:k) {
    myfit <- myfit + results$eta[subj, kp, iter_s] * epoch_basis %*% results$lambda[[iter_s]][,,kp] %*% design_var[subj,]
  }
  
  lines(myfit[interval], col = "red")
}






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

x <- MASS::mvrnorm(50, mu = rep(0, basis_df), Sigma = solve(10*penalty + .1*diag(basis_df)))
matplot(epoch_basis %*% t(x), type = "l")

MASS::mvrnorm(1, mu = rep(0,48), Sigma = 
                solve(results$tau2[1,10000] * var_penalty[[1]] +
                        results$phi_delta[,,1]))

curr_k <- 3
prior_precision <- results$tau2[1,curr_k,5000] * var_penalty[[1]] +
  results$phi_delta[,,curr_k]
precision <- kronecker(t(design_var) %*% diag(results$eta[,curr_k,10000]) %*%
                         diag(results$varphi[,10000]) %*% diag(results$eta[,curr_k,10000]) %*%
                         design_var, t(epoch_basis) %*% epoch_basis)
posterior_precision <- prior_precision + precision
MASS::mvrnorm(1, mu = rep(0, 48), Sigma = solve(posterior_precision))

tau_a <- 1
tau_b <- .0005
update_b <- 0
update_a <- 0
for (curr_k in 1:k) {
  update_a <- update_a + dim(results$lambda[[10000]])[1]
  update_b <- update_b + results$lambda[[10000]][,,curr_k] %*% 
    var_penalty[[1]] %*% 
    results$lambda[[10000]][,,curr_k]
}
post_a <- tau_a + update_a
post_b <- 1 / (tau_b + update_b)


### Check if lambda[,,,1] is still tiny compared to lambda[,,10] with no shrinkage
### Shrinkage isn't doing anything when smoothing is in effect
results$lambda[[10000]][,,10]
