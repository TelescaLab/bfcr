library(BayesianConditionalFPCA)
library(tidyverse)
library(mgcv)
library(spam)
library(plotly)
library(microbenchmark)

load(paste0("/Users/johnshamshoian/Documents/R_projects/", 
            "BayesianConditionalFPCA/sleep/data/",
            "relative_psd.RData"))
sleep_tabulated <- read.csv(paste0("/Users/johnshamshoian/Documents/R_projects/",
            "BayesianConditionalFPCA/sleep/tabulated_data/",
            "shhs1-dataset-0.15.0.csv"), stringsAsFactors = FALSE)

num_epochs <- 60
id_range <- 200001:200100

sleep_tabulated_filtered <- sleep_tabulated %>%
  select(nsrrid, age_s1)
sleep_data <- inner_join(sleep_tabulated_filtered, relative_psd)

sleep_data_filtered <- sleep_data %>%
  group_by(nsrrid) %>%
  filter(n() >= num_epochs, epoch <= num_epochs, 
         nsrrid %in% id_range) %>%
  ungroup()

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
penalty <- as.matrix(precmat.IGMRFreglat(age_df, epoch_df), order = 2)
epoch_marginal_penalty <- as.matrix(precmat.RW2(epoch_df), order = 2)
age_marginal_penalty <- as.matrix(precmat.RW2(age_df), order = 2)
Q1 <- kronecker(diag(age_df), epoch_marginal_penalty)
Q2 <- kronecker(age_marginal_penalty, diag(epoch_df))
age_penalty <- kronecker(penalty_2, diag(epoch_df))
epoch_penalty <- kronecker(diag(age_df), penalty_1)
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
iter <- 15000
burnin <- 5000
nchains <- 1
thin <- 1
loglik <- 0
mcmc_results <- run_mcmc(response, design_mean,
                  design_var, epoch_basis, 
                  epoch_grid,
                  mean_penalty, var_penalty,
                  mean_indices, var_indices,
                  k, iter, burnin, thin = 1,
                  var = "unequal")


subject_bands <- get_posterior_subject_bands(mcmc_results)
mean_bands <- get_posterior_means(mcmc_results, c(1))
eigen_bands <- get_posterior_eigen(mcmc_results, 6, c(1))
subj <- c(1, 10, 20, 30)
subject_bands %>% 
  filter(id %in% subj) %>%
  ggplot() +
  geom_point(aes(x = time, y = response), alpha = .5) +
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper), alpha = 0.5) +
  facet_wrap(. ~ id) +
  theme_bw()

number.labs <- paste0("Eigenfunction ", 1:6, ": ", 100 * round(eigen_bands$prop_var_explained[2,],2), "%",
                      " (", 100 * round(eigen_bands$prop_var_explained[1,], 2), "% - ",
                      100 * round(eigen_bands$prop_var_explained[3,], 2), "%)")
names(number.labs) <- c("1":6)
eigen_bands$eigenfunctions %>%
  ggplot(aes(x=time)) +
  facet_wrap(. ~ number, labeller = labeller(number = number.labs)) +
  geom_line(aes(y=mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .3) +
  labs(x = "Epoch", y = "Value") +
  theme_bw()

plot_ly() %>%
  add_surface(z =~ eigen_bands$surface)

mean_bands %>% 
  ggplot(mapping = aes(time)) +
  geom_line(aes(y=mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .3) +
  labs(x = "Epoch", y = "Relative delta power spectral density") +
  theme_bw()





z <- 1000
tau <- 100000000
prior_precision <- tau * penalty + z * diag(48)
mycov <- solve(prior_precision)
x <- MASS::mvrnorm(100, mu = rep(0, 48), Sigma = mycov)
round(mycov[1:10,1:10],2)
matplot(t(x), type = "l")
