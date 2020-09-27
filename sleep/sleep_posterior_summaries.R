library(BayesianConditionalFPCA)
library(tidyverse)
library(mgcv)
library(spam)
library(plotly)
library(dlnm)
library(fdapace)

mcmc_results <- readRDS(
  paste0("/Users/johnshamshoian/Documents/R_projects/bfcr",
  "/sleep/mcmc_output/mcmc_results_age.rds")
)

load(
  paste0("/Users/johnshamshoian/Documents/R_projects/",
         "bfcr/sleep/data/",
         "relative_psd.RData"))

sleep_tabulated <- read.csv(
  paste0("/Users/johnshamshoian/Documents/R_projects/",
         "bfcr/sleep/tabulated_data/", 
         "shhs1-dataset-0.15.0.csv"),
  stringsAsFactors = FALSE)

num_epochs <- 120
id_range <- 200001:201000
covariate <- "ahi_a0h3"

sleep_data <- sleep_tabulated %>%
  filter(EEG1qual == 4) %>%
  select(nsrrid, !!(sym(covariate))) %>%
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
covariate_df <- 6
epoch_basis <- ps(epoch_grid, df = epoch_df, intercept = TRUE)
epoch_penalty <- attr(epoch_basis, "S")

# Covariate adjusted
evaluate_spline <- function(smoothCon_list, covariate) {
  c(1, PredictMat(smoothCon_list[[1]],
                  data = data.frame(covariate_grid = covariate)))
}


# Age adjusted, using mgcv
if (TRUE) {
  covariate_basis <- smoothCon(s(covariate_grid, bs = "ps",
                                 k = covariate_df, m = 2),
                               data = data.frame(covariate_grid),
                               absorb.cons = TRUE)
  design_mean <- cbind(1, covariate_basis[[1]]$X)
  design_var <- cbind(1, covariate_basis[[1]]$X)
  model_penalties <- tensor.prod.penalties(list(covariate_basis[[1]]$S[[1]], epoch_penalty))
  mean_indices <- c(1, 2, 2)
  var_indices <- c(1, 2, 2)
  mean_penalty <- list(epoch_penalty, model_penalties[[1]], model_penalties[[2]])
  var_penalty <- list(epoch_penalty, model_penalties[[1]], model_penalties[[2]])
  covariate_quantiles <- quantile(covariate_grid, c(.1, .9))
  new_covariate <- seq(from = covariate_quantiles[1],
                       to = covariate_quantiles[2], length.out = 10)
}


mean_surface <- array(0, dim = c(length(new_covariate), num_epochs, 3))
for (i in seq_along(new_covariate)) {
  new_covariate_value <- new_covariate[i]
  xi <- evaluate_spline(covariate_basis, new_covariate_value)
  mean_bands <- get_posterior_means(
    mcmc_results, 
    xi)
  mean_surface[i,,2] <- mean_bands$mean
  mean_surface[i,,1] <- mean_bands$lower
  mean_surface[i,,3] <- mean_bands$upper
}

plot_ly(showscale = FALSE) %>%
  add_surface(z=~mean_surface[,,2], x=~epoch_grid, y=~new_covariate) %>%
  layout(scene = list(
    xaxis = list(title = "Epoch"),
    yaxis = list(title = "AHI"),
    zaxis = list(title = "Normalized Delta Power")
  )) %>%
  add_surface(z=~mean_surface[,,1], x=~epoch_grid, y=~new_covariate) %>%
  layout(scene = list(
    xaxis = list(title = "Epoch"),
    yaxis = list(title = "AHI"),
    zaxis = list(title = "Normalized Delta Power")
  )) %>%
  add_surface(z=~mean_surface[,,3], x=~epoch_grid, y=~new_covariate) %>%
  layout(scene = list(
    xaxis = list(title = "Epoch"),
    yaxis = list(title = "AHI"),
    zaxis = list(title = "Normalized Delta Power")
  )) 

evals <- 3
eigen_surface <- array(0, dim = c(length(new_covariate), num_epochs, 3, evals))
pve <- array(0, dim = c(length(new_covariate), 3, evals))
magnitude <- matrix(0, 3, length(new_covariate))
cov_surface <- array(0, dim = c(length(new_covariate), num_epochs, num_epochs))
for (i in seq_along(new_covariate)) {
  print(i)
  new_covariate_value <- new_covariate[i]
  zi <- evaluate_spline(covariate_basis, new_covariate_value)
  eigen_bands <- get_posterior_eigen(mcmc_results, evals, zi)
  for (k in 1:evals) {
    eigen_surface[i,,2,k] <- pull(eigen_bands$eigenfunctions %>%
                                  filter(number == k) %>%
                                  select(mean))
    eigen_surface[i,,1,k] <- pull(eigen_bands$eigenfunctions %>%
                                    filter(number == k) %>%
                                    select(lower))
    eigen_surface[i,,3,k] <- pull(eigen_bands$eigenfunctions %>%
                                    filter(number == k) %>%
                                    select(upper))
    cov_surface[i,,] <- eigen_bands$surface
    pve[i,,k] <- eigen_bands$prop_var_explained[,k]
  }
  magnitude[,i] <- eigen_bands$magnitude
}

cov_surface <- matrix(0, num_epochs, num_epochs)
iter < 5000
for (k in 1:12) {
  cov_surface <- cov_surface + 
    tcrossprod(epoch_basis %*% mcmc_results$samples$lambda[[iter]][,,k] %*%
                c(1,rep(0,5)))
}

plot_ly(showscale = FALSE) %>%
  add_surface(z=~eigen_surface[,,2,1], x=~epoch_grid, y=~new_covariate) %>%
  layout(scene = list(
    xaxis = list(title = "Epoch"),
    yaxis = list(title = "Age (years"),
    zaxis = list(title = "Normalized Delta Power")
  ))  
  add_surface(z=~eigen_surface[,,3,1], x=~epoch_grid, y=~new_ages) %>%
  layout(scene = list(
    xaxis = list(title = "Epoch"),
    yaxis = list(title = "Age (years"),
    zaxis = list(title = "Normalized Delta Power")
  ))

magnitude_tibble <- tibble(covariate = new_covariate,
                           lower = magnitude[1,],
                           upper = magnitude[3,],
                           mean = magnitude[2,])

magnitude_tibble %>%
  ggplot(aes(x = new_covariate)) +
  geom_line(aes(y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .3) +
  labs(x = "AHI", y = "Total variance") +
  theme_bw()

