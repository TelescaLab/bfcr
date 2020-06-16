set.seed(1)
library(MASS)
library(plot3D)
library(dlnm)
library(BayesianConditionalFPCA)
library(rstan)
library(loo)
library(tidyverse)
library(pracma)
setwd("/Users/johnshamshoian/Rcpp/BayesianConditionalFPCA/Rfuns")
source('simulated_data.r')
#matplot(t, t(Y), type = "l", col = "black", xlab = "Time", ylab = "Response")
par(mfrow = c(1,2))
matplot(t, t(Y[1:(n/2),]), type = "l", xlab = "Time", ylab = "Response", col = "black", main = "Group 1")
matplot(t, t(Y[(n/2 + 1):n,]), type = "l", xlab = "Time", ylab = "Response", col = "black", main = "Group 2")
par(mfrow = c(1,1))
### MCMC ###
K <- 4
Basis <- ps(t, df = 16, intercept = TRUE)
X_red <- cbind(X, X[,2]^2, X[,2]^3)
mcmc_results <- run_mcmc_Morris(Y, t, X, X, Basis, K, iter = 10000, burnin = 30000, nchains = 1, thin = 4, loglik = 0)
results <- numeric(20)
### Visualization ###
# sub <- 4
posterior_intervals <- get_posterior_predictive_bands2(mcmc_results, c(.025, .5, .975))
colnames(posterior_intervals) <- c("ID", "Time", "Y", "Lower_P", "Median_P", "Upper_P", "Lower_M", "Median_M", "Upper_M")
posterior_intervals <- as_tibble(posterior_intervals)
posterior_intervals$Y_no_error <- c(t(Y_no_error))
# posterior_intervals %>%
#   filter(ID == sub) %>%
#   ggplot(aes(x = Time, y = Y_no_error)) +
#   geom_point(na.rm = TRUE) +
#   geom_ribbon(aes(ymin = Lower_P, ymax = Upper_P), alpha = 0.3) +
#   theme_bw()
# posterior_intervals %>%
#   filter(ID == sub) %>%
#   ggplot(aes(x = Time, y = Y)) +
#   geom_point(na.rm = TRUE) +
#   geom_line(aes(y = Median_M)) +
#   theme_bw()
# posterior_intervals %>%
#   filter(ID == sub) %>%
#   ggplot(aes(x = Time, y = Y_no_error)) +
#   geom_point(na.rm = TRUE) +
#   geom_ribbon(aes(ymin = Lower_M, ymax = Upper_M), alpha = 0.3) +
#   theme_bw()
# posterior_intervals %>%
#   group_by(ID) %>%
#   filter(Y_no_error > Lower_M & Y_no_error < Upper_M) %>%
#   summarize(coverage = n() / tmax)
# posterior_intervals %>%
#   filter(Y > Lower_P & Y < Upper_P) %>%
#   summarize(coverage = n() / (tmax * n))

results[1] <- unlist(posterior_intervals %>%
  group_by(ID) %>%
  summarize(mymean = trapz(t,(Median_P - Y_no_error)^2) / 
              trapz(t, (Y_no_error)^2)) %>%
  ungroup() %>%
  summarize(mmean = mean(mymean)))
  
results[2] <- unlist(posterior_intervals %>%
  filter(Y_no_error > Lower_M & Y_no_error < Upper_M) %>%
  summarize(coverage = n() / (n * tmax)))

results[3] <- unlist(posterior_intervals %>%
  filter(Y > Lower_P & Y < Upper_P) %>%
  summarize(coverage = n() / (n * tmax)))

results[4] <- mean(posterior_intervals$Upper_P - posterior_intervals$Lower_P)
results[5] <- mean(posterior_intervals$Upper_M - posterior_intervals$Lower_M)
### Goodness of fit ###
# mcmc_omnibus_fit <- get_omnibus_fit2(mcmc_results)
# plot(mcmc_omnibus_fit$statistic_obs, mcmc_omnibus_fit$statistic_rep, xlab = "Chi2 observed", ylab = "Chi2 repitition")
# abline(a = 0, b = 1)
# sum(mcmc_omnibus_fit$statistic_rep > mcmc_omnibus_fit$statistic_obs) / 5000

### LOO ###
# r_eff <- relative_eff(exp(mcmc_results$log_lik), cores = 4)
# loo_1 <- loo(mcmc_results$log_lik, r_eff = r_eff, cores = 4)
# print(loo_1)

### Posterior mean bands ###
alpha <- .05
# sub <- 13
# xi <- X[sub,]
# xi <- c(0, 1)
# coef_bands <- get_posterior_means(mcmc_results, xi, alpha)
# coef_bands <- cbind(t, coef_bands)
# colnames(coef_bands) <- c("Time", "Lower", "Mean", "Upper")
# coef_bands <- as_tibble(coef_bands)
# coef_bands %>%
#   ggplot(aes(x = Time)) +
#   geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, fill = "blue") +
#   geom_line(aes(y = Mean)) +
#   xlab("Time") +
#   ylab("Response") +
#   ggtitle("Conditional mean bands") +
#   theme_bw()
# 
# plot(t, coef_bands$Upper, type = "l", ylim = c(-.2, 1.2), xlab = "Time", ylab = "Response")
# lines(t, coef_bands$Lower)
# lines(t, Btru%*%Theta1 %*% xi,col="green")

coef_bands <- get_posterior_coefs(mcmc_results, alpha)

Intercept_fn <- Btru %*% Theta1 %*% c(1, 0)
Age_fn <- Btru %*% Theta1 %*% c(0, 1)
coef_bands <- tibble(Frequency = rep(t, times = dim(X)[2]),
                     Covariate = rep(c("Intercept", "Age"), each = dim(Basis)[1]),
                     Lower = c(coef_bands$lower),
                     Mean = c(coef_bands$mean),
                     Upper = c(coef_bands$upper))
Intercept_idx <- which(coef_bands$Covariate == "Intercept")
results[6] <- trapz(t, (coef_bands$Mean[Intercept_idx] - Intercept_fn)^2) / 
  trapz(t, Intercept_fn^2)
results[7] <- trapz(t, (coef_bands$Mean[-Intercept_idx] - Age_fn)^2) /
  trapz(t, Age_fn^2)
results[8] <- mean(as.numeric((Intercept_fn >= coef_bands$Lower[Intercept_idx] &
                 Intercept_fn <= coef_bands$Upper[Intercept_idx]))) * 100
results[9] <- mean(as.numeric((Age_fn >= coef_bands$Lower[-Intercept_idx] &
                                 Age_fn <= coef_bands$Upper[-Intercept_idx]))) *
                                 100
results[10] <- mean(coef_bands$Upper[Intercept_idx] - coef_bands$Lower[Intercept_idx])
results[11] <- mean(coef_bands$Upper[-Intercept_idx] - coef_bands$Lower[-Intercept_idx])
# coef_bands %>%
#   ggplot(aes(x = Frequency)) +
#   geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Covariate), alpha = 0.5) +
#   geom_line(aes(y = Mean)) +
#   geom_hline(yintercept = 0) +
#   facet_wrap(Covariate ~., scales = "free") +
#   ylab("Power") +
#   theme_minimal()

### Some covariance visualization ###
z_seq <- seq(from = -2, to = 2, length.out = 10)
evals <- 2
for(i in 1:length(z_seq)){
  
  zi <- c(1, z_seq[i])
  alpha <- .05
  eigen_bands <- get_posterior_eigen2(mcmc_results, evals, zi, alpha)
  eig_names <- c()
  eig_labs <- c()
  for(k in 1:evals){
    eig_names <- c(eig_names, paste("Eigenfunction", k))
    eig_labs <- c(eig_labs, paste("Eigenfunction", k, " ", round(eigen_bands$eigenval_pve_intervals[1,k], 2), "-", round(eigen_bands$eigenval_pve_intervals[3,k], 2)))
  }
  names(eig_labs) <- eig_names
  eigen_bands_tibble <- tibble(Time = rep(t, evals),
                               number = factor(rep(eig_names, each = length(t))),
                               lower = c(eigen_bands$lower),
                               mean = c(eigen_bands$mean),
                               upper = c(eigen_bands$upper),
                               val_lower = rep(eigen_bands$eigenval_intervals[1,], each = length(t)),
                               val_median = rep(eigen_bands$eigenval_intervals[2,], each = length(t)),
                               val_upper = rep(eigen_bands$eigenval_intervals[3,], each = length(t)))
  
  truecov <- matrix(0, nrow = tmax, ncol = tmax)
  truecov <- Btru %*% Lambda1 %*% outer(zi, zi) %*% t(Lambda1) %*% t(Btru) +
    Btru %*% Lambda2 %*% outer(zi, zi) %*% t(Lambda2) %*% t(Btru)
  eigvec <- eigen(truecov)$vectors
  for(i in 1:2){
    eigvec[,i] <- eigvec[,i] / sqrt(trapz(t, eigvec[,i]^2))
  }
  eigen_bands_tibble <- eigen_bands_tibble %>%
    mutate(eigvec = c(eigvec[,1], eigvec[,2]))
  
  eigen_bands_tibble <- eigen_bands_tibble %>%
    group_by(number) %>%
    mutate(eigvec = case_when(as.numeric(trapz(t, (mean + eigvec)^2) < trapz(t, 
    (mean - eigvec)^2)) == 1 ~ -eigvec, TRUE ~ eigvec))
  
  
  results[11:12] <- 100 * unlist((eigen_bands_tibble %>%
    group_by(number) %>%
    summarize(RISE = trapz(t, (mean - eigvec)^2) / trapz(t, eigvec^2)))[,2]) +
    results[11:12]
  
  results[13:14] <- unlist((eigen_bands_tibble %>%
    group_by(number) %>%
    filter(eigvec >= lower & eigvec <= upper) %>%
    summarize(coverage = 100 * n() / tmax))[,2]) + results[13:14]
  
  results[15:16] <- unlist((eigen_bands_tibble %>%
    group_by(number) %>%
    summarize(mean_l = mean(upper - lower)))[,2]) + results[15:16]
  results[17] <- trapz(t,sapply(1:tmax, function(i) trapz(t, (eigen_bands$surface[,i] -
        truecov[,i])^2))) / trapz(t, sapply(1:tmax, function(i) trapz(t, 
        truecov[,i]^2))) * 100 + results[17]
}
results[11:17] <- results[11:17] / length(z_seq)

plot(eigen_bands_tibble$eigvec[1:50], type = "l")
lines(eigen_bands_tibble$mean[1:50], col = "green")
plot(eigen_bands_tibble$eigvec[51:100], type = "l")
lines(eigen_bands_tibble$mean[51:100], col = "green")
plot(eigen_bands$raw_magnitude, type = "l")

# eigen_bands_tibble %>%
#   ggplot(aes(x = Time)) +
#   geom_ribbon(aes(ymin = lower, ymax = upper, fill = number), alpha = 0.5) +
#   geom_line(aes(y = mean)) +
#   facet_wrap(number ~., labeller = labeller(number = eig_labs)) +
#   theme_bw() +
#   theme(legend.position="none")
# eigen_bands$magnitude
# aX <- list(title = "Time")
# aY <- list(title = "Time")
# aZ <- list(title = "Response")
# plotly::plot_ly(x = t, y = t, z = eigen_bands$surface, type = "surface") %>%
#   plotly::layout(scene = list(xaxis = aX, yaxis = aY, zaxis = aZ, dragmode = "turntable"))
# plot(eigen_bands$raw_magnitude, type = "l")
# truecov <- matrix(0, nrow = tmax, ncol = tmax)
# truecov <- Btru %*% Lambda1 %*% outer(X[sub,], X[sub,]) %*% t(Lambda1) %*% t(Btru) +
#   Btru %*% Lambda2 %*% outer(X[sub,], X[sub,]) %*% t(Lambda2) %*% t(Btru)
# par(mfrow = c(1,2))
# persp3D(t, t, eigen_bands$surface)
# persp3D(t, t, truecov)
# par(mfrow = c(1,1))


# X_seq <- seq(from = -2, to = 2, by = .1)
# Cov_mat <- matrix(0, nrow = length(X_seq), ncol = 3)
# counter <- 1
# for(x in X_seq){
#   zi <- c(1, x)
#  Cov_mat[counter,] <- get_posterior_eigen2(mcmc_results, 4, zi, .05)$magnitude
#   truecov <- matrix(0, nrow = tmax, ncol = tmax)
#   truecov <- Btru %*% Lambda1 %*% outer(zi, zi) %*% t(Lambda1) %*% t(Btru) +
#     Btru %*% Lambda2 %*% outer(zi, zi) %*% t(Lambda2) %*% t(Btru)
#   Cov_mat[counter] <- sum(eigen(truecov)$values)
#   counter <- counter + 1
# }
# plot(Cov_mat,type="l")
# plot(Cov_mat[,1], type = "l")
# lines(Cov_mat[,3])




