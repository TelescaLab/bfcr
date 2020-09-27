library(BayesianConditionalFPCA)
library(tidyverse)
library(mgcv)
library(spam)
library(plotly)
library(dlnm)
library(fdapace)

load(paste0("/Users/johnshamshoian/Documents/R_projects/",
            "bfcr/sleep/data/",
            "relative_psd.RData"))
sleep_tabulated <- read.csv(
  paste0("/Users/johnshamshoian/Documents/R_projects/",
         "bfcr/sleep/tabulated_data/",
         "shhs1-dataset-0.15.0.csv"), stringsAsFactors = FALSE)

num_epochs <- 120
id_range <- 200001:206000
covariate <- "HTNDerv_s1"

sleep_data <- sleep_tabulated %>%
  filter(EEG1qual == 4) %>%
  select(nsrrid, !!(sym(covariate))) %>%
  drop_na()

sleep_data <- sleep_tabulated %>%
  filter(EEG1qual == 4) %>%
  select(nsrrid, HTNDerv_s1, gender) %>%
  drop_na()
sleep_data <- inner_join(sleep_data, relative_psd)


sleep_data <- sleep_data %>%
  group_by(nsrrid) %>%
  filter(n() >= num_epochs, epoch <= num_epochs,
         nsrrid %in% id_range) %>%
  ungroup()

num_subjects <- as.numeric(summarise(
  sleep_data, n_distinct(nsrrid)))

epoch_grid <- 1:num_epochs
epoch_df <- ceiling(20 / 100 * num_epochs)
covariate_grid <- pull(sleep_data %>%
                         group_by(nsrrid) %>%
                         filter(row_number() == 1) %>%
                         ungroup() %>%
                         select(!!(sym(covariate))))

covariate_grid <- sleep_data %>%
  group_by(nsrrid) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(-nsrrid, -epoch, -psd) %>%
  as.matrix()
covariate_df <- 6
epoch_basis <- ps(epoch_grid, df = epoch_df, intercept = TRUE)
epoch_penalty <- attr(epoch_basis, "S")

# Covariate adjusted
evaluate_spline <- function(smoothCon_list, covariate) {
  c(1, PredictMat(smoothCon_list[[1]], data = data.frame(Bmi = bmi)))
}


# Age adjusted, using mgcv
if (TRUE) {
  # covariate_basis <- smoothCon(s(covariate_grid, bs = "ps",
  #                                k = covariate_df, m = 2),
  #                              data = data.frame(covariate_grid),
  #                              absorb.cons = TRUE)
  # design_mean <- cbind(1, covariate_basis[[1]]$X)
  # design_var <- cbind(1, covariate_basis[[1]]$X)
  # model_penalties <- tensor.prod.penalties(list(covariate_basis[[1]]$S[[1]], epoch_penalty))
  # mean_indices <- c(1, 2, 2)
  # var_indices <- c(1, 2, 2)
  # mean_penalty <- list(epoch_penalty, model_penalties[[1]], model_penalties[[2]])
  # var_penalty <- list(epoch_penalty, model_penalties[[1]], model_penalties[[2]])
  # covariate_quantiles <- quantile(covariate_grid, c(.1, .9))
  # new_covariate <- seq(from = covariate_quantiles[1],
  #                      to = covariate_quantiles[2], length.out = 10)
  # 
  design_mean <- cbind(rep(1, num_subjects), covariate_grid)
  design_var <- cbind(rep(1, num_subjects), covariate_grid)
  mean_indices <- c(1, 2, 3)
  var_indices <- c(1, 2, 3)
  mean_penalty <- list(epoch_penalty, epoch_penalty, epoch_penalty)
  var_penalty <- list(epoch_penalty, epoch_penalty, epoch_penalty)
}

# Not age adjusted
if (FALSE) {
  design_mean <- cbind(rep(1, num_subjects))
  design_var <- cbind(rep(1, num_subjects))
  mean_indices <- c(1)
  var_indices <- c(1)
  mean_penalty <- list(epoch_penalty)
  var_penalty <- list(epoch_penalty)
}
response <- t(matrix(sleep_data$psd,
                     nrow = num_epochs,
                     ncol = num_subjects))

k <- 8
iter <- 5000
burnin <- 2500
thin <- 1
mcmc_results <- run_mcmc(response, design_mean,
                         design_var, epoch_basis,
                         epoch_grid,
                         mean_penalty, var_penalty,
                         mean_indices, var_indices,
                         k, iter, burnin, thin = thin,
                         var = "unequal")
saveRDS(mcmc_results, file = "mcmc_results_systbp.rds")
subject_bands <- get_posterior_subject_bands(mcmc_results)
mean_bands <- get_posterior_means(mcmc_results, mcmc_results$data$design_mean[2,], alpha_level = .05)
evals <- 4
eigen_bands <- get_posterior_eigen(mcmc_results, evals, c(1,1,1))
subj <- 31:34
subject_bands %>%
  filter(id %in% subj) %>%
  ggplot() +
  geom_point(aes(x = time, y = response), alpha = .5) +
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper), alpha = 0.5) +
  facet_wrap(. ~ id) +
  theme_bw()

number.labs <- paste0("Eigenfunction ", 1:evals, ": ", 100 * round(eigen_bands$prop_var_explained[2,],2), "%",
                      " (", 100 * round(eigen_bands$prop_var_explained[1,], 2), "% - ",
                      100 * round(eigen_bands$prop_var_explained[3,], 2), "%)")
names(number.labs) <- c("1":evals)
eigen_bands$eigenfunctions %>%
  ggplot(aes(x=time)) +
  facet_wrap(. ~ number, labeller = labeller(number = number.labs), scales = "free") +
  geom_line(aes(y=mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .3) +
  labs(x = "Epoch", y = "Value") +
  theme_bw()

plot_ly() %>%
  add_surface(z =~ eigen_bands$surface)

L <- list()
L$y <- lapply(1:num_subjects, function(i) response[i,])
L$t <- lapply(1:num_subjects, function(i) 1:num_epochs)
res <- FPCA(L$y, L$t, list(dataType = "Dense", methodMuCovEst = "smooth"))

plot_ly() %>%
  add_surface(z =~ res$fittedCov)
mean_bands %>%
  ggplot(mapping = aes(time)) +
  geom_line(aes(y=mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .3) +
  labs(x = "Epoch", y = "Relative delta power spectral density") +
  theme_bw()
