library(BayesianConditionalFPCA)
library(mgcv)
library(tidyverse)
library(fdapace)
library(plotly)
source("/Users/johnshamshoian/Documents/R_projects/bfcr/Peak Alpha Data/Peak_Alpha_Data_Transfer.R")
num_subjects <- pa.dat %>%
  summarise(n_distinct(ID)) %>%
  pull()
freq_grid <- pa.dat %>%
  distinct(func) %>%
  pull()
age_grid <- pa.dat %>% 
  group_by(ID) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(Age) %>%
  pull()
group_grid <- pa.dat %>%
  group_by(ID) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  summarise(group = group - 1) %>%
  select(group) %>%
  pull()

num_times <- pa.dat %>%
  filter(ID == 2, reg == 1) %>%
  count() %>%
  pull()
response <- pa.dat %>%
  filter(reg == 15) %>%
  select(y) %>%
  pull() %>%
  matrix(nrow = num_times, ncol = num_subjects) %>%
  t()

freq_basis_obj <- smoothCon(s(freq_grid, k = 12, bs = "ps", m = 2),
                        data.frame(freq_grid),
                        absorb.cons = FALSE)
# freq_basis <- cbind(1, freq_basis_obj[[1]]$X)
freq_basis <- freq_basis_obj[[1]]$X
# freq_penalty <- matrix(0, 12, 12)
# freq_penalty[-1,-1] <- freq_basis_obj[[1]]$S[[1]] * freq_basis_obj[[1]]$S.scale
freq_penalty <- freq_basis_obj[[1]]$S[[1]] * freq_basis_obj[[1]]$S.scale

{
  age_basis_obj <- smoothCon(s(age_grid, k = 6, bs = "ps", by = as.factor(group_grid)),
                             data = data.frame(age_grid, group_grid),
                             absorb.cons = TRUE)
  age_penalty <- age_basis_obj[[1]]$S[[1]] * age_basis_obj[[1]]$S.scale
  model_penalties <- tensor.prod.penalties(list(age_penalty, freq_penalty))
  design_mean <- cbind(1, group_grid,
                       age_basis_obj[[1]]$X + age_basis_obj[[2]]$X,
                       age_basis_obj[[2]]$X)
  design_var <- cbind(1, group_grid,
                       age_basis_obj[[1]]$X + age_basis_obj[[2]]$X,
                       age_basis_obj[[2]]$X)
  mean_penalty <- list(freq_penalty, freq_penalty,
                       model_penalties[[1]], model_penalties[[2]],
                       model_penalties[[1]], model_penalties[[2]])
  var_penalty <- list(freq_penalty, freq_penalty,
                      model_penalties[[1]], model_penalties[[2]],
                      model_penalties[[1]], model_penalties[[2]])
  var_indices <- c(1, 2, 3, 3, 4, 4)
  mean_indices <- c(1, 2, 3, 3, 4, 4)
  evaluate_basis <- function(obj, age, group) {
    spline_part1 <- 
      PredictMat(obj[[1]], data.frame(age_grid = age ,group_grid = group)) +
      PredictMat(obj[[2]], data.frame(age_grid = age, group_grid = group))
      
    spline_part2 <- 
      PredictMat(obj[[2]], data.frame(age_grid = age, group_grid = group))
    c(1, group, spline_part1, spline_part2)
  }
}

{
  age_basis_obj <- smoothCon(s(age_grid, k = 6, bs = "ps"),
                             data = data.frame(age_grid, group_grid),
                             absorb.cons = TRUE)
  age_penalty <- age_basis_obj[[1]]$S[[1]]
  model_penalties <- tensor.prod.penalties(list(age_penalty, freq_penalty))
  design_mean <- cbind(1, group_grid, age_basis_obj[[1]]$X)
  design_var <- cbind(1, group_grid, age_basis_obj[[1]]$X)
  mean_penalty <- list(freq_penalty, freq_penalty,
                       model_penalties[[1]], model_penalties[[2]])
  var_penalty <- list(freq_penalty, freq_penalty,
                      model_penalties[[1]], model_penalties[[2]])
  var_indices <- c(1, 2, 3, 3)
  mean_indices <- c(1, 2, 3, 3)
  evaluate_basis <- function(age_basis_obj, age, group) {
    spline_part <- PredictMat(age_basis_obj[[1]],
                              data.frame(age_grid = age, group_grid = group))
    c(1, group, spline_part)
  }
}

{
  design_mean <- cbind(rep(1, num_subjects))
  design_var <- cbind(rep(1, num_subjects))
  mean_indices <- c(1)
  var_indices <- c(1)
  mean_penalty <- list(freq_penalty)
  var_penalty <- list(freq_penalty)
}

k <- 6
iter <- 10000
burnin <- 5000
thin <- 5

mcmc_results2 <- run_mcmc(response, design_mean,
                         design_var, freq_basis,
                         freq_grid,
                         mean_penalty, var_penalty,
                         mean_indices, var_indices,
                         k, iter, burnin, thin = thin,
                         var = "pooled")

subject_bands <- get_posterior_subject_bands(mcmc_results)
mean_bands <- get_posterior_means(mcmc_results, mcmc_results$data$design_mean[1,])
evals <- 4
eigen_bands <- get_posterior_eigen(mcmc_results, evals, mcmc_results$data$design_mean[8,])
subj <- 8:10
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
L$t <- lapply(1:num_subjects, function(i) 1:num_times)
res <- FPCA(L$y, L$t, list(dataType = "Dense", methodMuCovEst = "smooth"))

plot_ly() %>%
  add_surface(z =~ res$fittedCov)
plot(res$mu, type = "l")
lines(mean_bands$mean, col = "blue")
mean_bands %>%
  ggplot(mapping = aes(time)) +
  geom_line(aes(y=mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .3) +
  labs(x = "Frequency", y = "Relative delta power spectral density") +
  theme_bw()

N <- 10
num_time <- 33
lower_age <- quantile(age_grid, .1)
upper_age <- quantile(age_grid, .9)
new_age_grid <- seq(from = lower_age, to = upper_age, length.out = N)
mean_surface_td <- matrix(0, N, num_time)
for (i in 1:N) {
  xi <- evaluate_basis(age_basis_obj, new_age_grid[i], 0)
  mean_bands <- get_posterior_means(mcmc_results, xi)
  mean_surface_td[i, ] <- mean_bands$mean
}

plot_ly() %>% add_surface(z =~ mean_surface_td, y =~ new_age_grid, x =~ freq_grid)

N <- 10
num_time <- 33
lower_age <- quantile(age_grid, .1)
upper_age <- quantile(age_grid, .9)
new_age_grid <- seq(from = lower_age, to = upper_age, length.out = N)
mean_surface_asd <- matrix(0, N, num_time)
for (i in 1:N) {
  xi <- evaluate_basis(age_basis_obj, new_age_grid[i], 1)
  mean_bands <- get_posterior_means(mcmc_results, xi)
  mean_surface_asd[i, ] <- mean_bands$mean
}

plot_ly() %>% add_surface(z =~ mean_surface_asd, y =~ new_age_grid, x =~ freq_grid)

eigen_surface_td <- matrix(0, N, num_time)
magnitude_td <- matrix(0, N, 3)
for (i in 1:N) {
  zi <- evaluate_basis(age_basis_obj, new_age_grid[i], 0)
  eigen_bands <- get_posterior_eigen(mcmc_results, 1, zi)
  eigen_surface_td[i,] <- eigen_bands$eigenfunctions$mean
  if (eigen_surface_td[i, 1] < 0) {
    eigen_surface_td[i,] <- -1 * eigen_surface_td[i,]
  }
  magnitude_td[i, ] <- eigen_bands$magnitude
}

plot(new_age_grid, magnitude_td[,1], type = "l", ylim = c(0, .05))
lines(new_age_grid, magnitude_td[,3])
plot_ly() %>% add_surface(z =~ eigen_surface_td, y =~ new_age_grid, x =~ freq_grid)

eigen_surface_asd <- matrix(0, N, num_time)
magnitude_asd <- matrix(0, N, 3)
for (i in 1:N) {
  zi <- evaluate_basis(age_basis_obj, new_age_grid[i], 1)
  eigen_bands <- get_posterior_eigen(mcmc_results, 1, zi)
  eigen_surface_asd[i,] <- eigen_bands$eigenfunctions$mean
  if (eigen_surface_asd[i, 1] < 0) {
    eigen_surface_asd[i,] <- -1 * eigen_surface_asd[i,]
  }
  magnitude_asd[i, ] <- eigen_bands$magnitude
}

plot(new_age_grid, magnitude_asd[,1], type = "l", ylim = c(0, .05))
lines(new_age_grid, magnitude_asd[,3])
plot_ly() %>% add_surface(z =~ eigen_surface_td, y =~ new_age_grid, x =~ freq_grid)


intercept_effect <- numeric(num_time)
group_effect <- numeric(num_time)
age_effect <- numeric(num_time)
interaction_effect <- numeric(num_time)
group_zi <- c(0, 1, rep(0, 10))
age_zi <- c(0, 0, PredictMat(age_basis_obj[[1]], data.frame(age_grid = 60, group_grid = 0)), rep(0, 5))
interaction_zi <- c(rep(0, 7), PredictMat(age_basis_obj[[1]], data.frame(age_grid = 60, group_grid = 0)))
for (iter in 2501:10000) {
  for (kp in 1:k) {
    intercept_effect <- intercept_effect + (freq_basis %*% mcmc_results$samples$lambda[[iter]][,1,kp])^2
    group_effect <- group_effect + (freq_basis %*% mcmc_results$samples$lambda[[iter]][,,kp] %*% group_zi)^2
    age_effect <- age_asd_effect + (freq_basis %*% mcmc_results$samples$lambda[[iter]][,,kp] %*% age_zi)^2
    interaction_effect <- interaction_effect + (freq_basis %*% mcmc_results$samples$lambda[[iter]][,,kp] %*% interaction_zi)^2
  }
}
intercept_effect <- intercept_effect / 7500
age_effect <- age_effect / 7500
group_effect <- group_effect / 7500
interaction_effect <- interaction_effect / 7500

plot(intercept_effect, type = "l")
lines(group_effect, col = "blue")
lines(age_effect, col = "green")
lines(interaction_effect, col = "orange")





fit1 <- numeric(33)
fit1 <- freq_basis %*% mcmc_results$samples$beta[,,5000] %*% evaluate_basis(age_basis_obj,  age_grid[8], group_grid[8])
for (kp in 1:k) {
  fit1 <- fit1 + mcmc_results$samples$eta[1,kp,5000] * freq_basis %*% 
    mcmc_results$samples$lambda[[5000]][,,kp] %*% 
    evaluate_basis(age_basis_obj, age_grid[8] , group_grid[8])
}

plot(response[1,])
lines(fit1)
total_cov <- array(0, dim = c(33, 33, 97))
for (i in 1:97) {
  mycov <- matrix(0, 33, 33)
  for (kp in 1:k) {
    mycov <- mycov + tcrossprod(freq_basis %*% 
                                  mcmc_results$samples$lambda[[5000]][,,kp] %*% 
                                  evaluate_basis(age_basis_obj, age_grid[i], group_grid[i]))
  }
  total_cov[,,i] <- mycov
}


plot_ly() %>%
  add_surface(z =~ mycov2)
l1 <- freq_basis %*% 
  mcmc_results$samples$lambda[[5000]][,,1] %*% 
  evaluate_basis(age_basis_obj, 90 , group_grid[1])
plot(l1, type = "l")
lines(l1)
plot(freq_basis %*% mcmc_results$samples$lambda[[5000]][,1,1], type = "l")
lines(freq_basis %*% mcmc_results$samples$lambda[[5000]][,1,5], col = "red")
plot(freq_basis %*% mcmc_results$samples$lambda[[5000]][,,1] %*% age_td, type = "l")
lines(freq_basis %*% mcmc_results$samples$lambda[[5000]][,,5] %*% age_td, col = "red")
plot(freq_basis %*% mcmc_results$samples$lambda[[5000]][,,1] %*% age_asd, type = "l")
lines(freq_basis %*% mcmc_results$samples$lambda[[5000]][,,5] %*% age_asd, col = "red")

fake_y <- freq_basis %*% mcmc_results$samples$beta[,,6000][,1] %*% rep(1, num_subjects)
for (kp in 1:k) {
  fake_y <- fake_y + freq_basis  %*%mcmc_results$samples$lambda[[6000]][,1,kp] %*% 
    mcmc_results$samples$eta[,kp,6000] 
}
fake_y <- t(fake_y)
matplot(t(fake_y), type = "l")
matplot(t(response), type = "l")


v1 <- get_prediction_error(mcmc_results, .05)
new_design_mean <- design_mean
new_design_var[,2] <- new_design_var[sample(1:97, 97, replace = F), 2]
mcmc_results$data$design_mean <- new_design_var
v2 <- get_prediction_error(mcmc_results, .05)
v2
