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

num_epochs <- 120
id_range <- 200001:200100

sleep_tabulated_filtered <- sleep_tabulated %>%
  select(nsrrid, age_s1, bmi_s1) %>%
  drop_na()
sleep_data <- inner_join(sleep_tabulated_filtered, relative_psd)

sleep_data_filtered <- sleep_data %>%
  group_by(nsrrid) %>%
  filter(n() >= num_epochs, epoch <= num_epochs,
         nsrrid %in% id_range) %>%
  ungroup()

num_subjects <- as.numeric(summarise(
  sleep_data_filtered, n_distinct(nsrrid)))

epoch_grid <- 1:num_epochs
epoch_df <- ceiling(20 / 100 * num_epochs)
age_grid <- sleep_data_filtered %>%
  group_by(nsrrid) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(age_s1)
age_df <- ceiling(10 / 100 * (max(age_grid) - min(age_grid)))
epoch_spec <- s(epoch_grid, bs = "ps", k = epoch_df)
age_spec <- s(age_grid, bs = "ps", k = age_df)
age_spec$term <- "age_s1"

epoch_basis <- smoothCon(epoch_spec,
                         data = data.frame(epoch_grid))
age_basis <- smoothCon(object = age_spec,
                       data = data.frame(age_grid))
epoch_marginal_penalty <- epoch_basis[[1]]$S[[1]]
age_marginal_penalty <- age_basis[[1]]$S[[1]]
model_penalties <- tensor.prod.penalties(list(age_marginal_penalty,
                                              epoch_marginal_penalty))
mean_penalty <- model_penalties
var_penalty <- model_penalties
mean_indices <- c(1,1)
var_indices <- c(1,1)
design_mean <- age_basis[[1]]$X
design_var <- age_basis[[1]]$X
# design_mean <- cbind(rep(1, num_subjects))
# design_var <- cbind(rep(1, num_subjects))
# mean_penalty <- list(16*epoch_marginal_penalty)
# var_penalty <- list(16*epoch_marginal_penalty)
# reduced_penalty <- matrix(0, nrow = epoch_df, epoch_df)
# reduced_penalty[2:epoch_df, 2:epoch_df] <- epoch_basis[[1]]$S[[1]]
# mean_penalty <- list(reduced_penalty)
# var_penalty <- list(reduced_penalty)
# mean_penalty <- list(matrix(0, 24, 24))
# var_penalty <- list(matrix(0, 24, 24))
# mean_indices <- c(1)
# var_indices <- c(1)

epoch_basis_spline <- epoch_basis[[1]]$X
response <- t(matrix(sleep_data_filtered$psd,
                   nrow = num_epochs,
                   ncol = num_subjects))

k <- 15
iter <- 5000
burnin <- 2500
thin <- 1
loglik <- 0
mcmc_results <- run_mcmc(response, design_mean,
                  design_var, epoch_basis_spline,
                  epoch_grid,
                  mean_penalty, var_penalty,
                  mean_indices, var_indices,
                  k, iter, burnin, thin = 1,
                  var = "unequal")

subject_bands <- get_posterior_subject_bands(mcmc_results)
mean_bands <- get_posterior_means(mcmc_results, design_mean[4,])
eigen_bands <- get_posterior_eigen(mcmc_results, 6, design_mean[5,])
subj <- 40:43
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
