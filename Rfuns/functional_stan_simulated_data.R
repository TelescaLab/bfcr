library(MASS)
library(splines)
library(mvtnorm)
library(plot3D)
library(splines)
library(Rcpp)
library(microbenchmark)
library(dlnm)
library(BayesianConditionalFPCA)
library(rstan)
library(loo)
library(tidyverse)
library(pracma)
library(expm)
setwd("E:/Rcpp stuff/BayesianConditionalFPCA/Rfuns")
source('simulated_data.r')



Y[1,1:10] <- NA
Y[2, 25:35] <- NA
Y[4,1:20] <- NA

### MCMC sanity check ###
K <- 2
mcmc_results <- run_mcmc(Y, t, X, Btru, K, iter = 5000, burnin = 2500, nchains = 1, thin = 1, loglik = 0)

### Visualization ###
posterior_intervals <- get_posterior_predictive_bands(mcmc_results, c(.05, .5, .95))
colnames(posterior_intervals) <- c("ID", "Time", "Y", "Lower_P", "Median_P", "Upper_P", "Lower_M", "Median_M", "Upper_M")
posterior_intervals <- as_tibble(posterior_intervals)
posterior_intervals %>%
  filter(ID == 3) %>%
  ggplot(aes(x = Time, y = Y)) +
  geom_point(na.rm = TRUE) +
  geom_ribbon(aes(ymin = Lower_P, ymax = Upper_P), alpha = 0.3) +
  theme_bw()
posterior_intervals %>%
  filter(ID == 1) %>%
  ggplot(aes(x = Time, y = Y)) +
  geom_point(na.rm = TRUE) +
  geom_ribbon(aes(ymin = Lower_M, ymax = Upper_M), alpha = 0.3) +
  theme_bw()
posterior_intervals %>%
  group_by(ID) %>%
  filter(Y > Lower_P & Y < Upper_P) %>%
  summarize(coverage = n() / tmax)
posterior_intervals %>%
  filter(Y > Lower_P & Y < Upper_P) %>%
  summarize(coverage = n() / (tmax * n))

### Goodness of fit ###
mcmc_omnibus_fit <- get_omnibus_fit(mcmc_results)
plot(mcmc_omnibus_fit$statistic_obs, mcmc_omnibus_fit$statistic_rep, xlab = "Chi2 observed", ylab = "Chi2 repitition")
abline(a = 0, b = 1)
sum(mcmc_omnibus_fit$statistic_rep > mcmc_omnibus_fit$statistic_obs) / 5000

### LOO ###
r_eff <- relative_eff(exp(mcmc_results$log_lik), cores = 4)
loo_1 <- loo(mcmc_results$log_lik, r_eff = r_eff, cores = 4)
print(loo_1)

### Posterior mean bands ###
xi <- c(1, 0)
alpha <- .05
coef_bands <- get_posterior_means(mcmc_results, xi, alpha)
coef_bands <- cbind(t, coef_bands)
colnames(coef_bands) <- c("Time", "Lower", "Mean", "Upper")
coef_bands <- as_tibble(coef_bands)
coef_bands %>%
  ggplot(aes(x = Time)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, fill = "blue") +
  geom_line(aes(y = Mean)) +
  xlab("Time") +
  ylab("Response") +
  ggtitle("Conditional mean bands") + 
  theme_bw() 


### Some covariance visualization ###
evals <- 6
xi <- c(1, .2)
alpha <- .05
eigen_bands <- get_posterior_eigen(mcmc_results, evals, xi, alpha)
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

eigen_bands_tibble %>%
  ggplot(aes(x = Time)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = number), alpha = 0.5) +
  facet_wrap(number ~., labeller = labeller(number = eig_labs)) +
  theme_bw() + 
  theme(legend.position="none")
aX <- list(title = "Time")
aY <- list(title = "Time")
aZ <- list(title = "Response")
plotly::plot_ly(x = t, y = t, z = eigen_bands$surface, type = "surface") %>%
  plotly::layout(scene = list(xaxis = aX, yaxis = aY, zaxis = aZ, dragmode = "turntable"))
persp3D(1:50, 1:50, eigen_bands$surface)
trapz1true <- trapz(t, eigen(truecov)$vectors[,1]^2)
trapz2true <- trapz(t, eigen(truecov)$vectors[,2]^2)
boobie <- get_posterior_eigen(mcmc_results, 3, c(1,1), .05)
plot(boobie$lower[,1], type = "l")
lines(boobie$upper[,1])
lines(-1/sqrt(trapz1true)*eigen(truecov)$vectors[,1], col = "red")

plot(boobie$lower[,2], type = "l", ylim = c(-4,4))
lines(boobie$upper[,2])
lines(1/sqrt(trapz2true)*eigen(truecov)$vectors[,2], col = "red")

plot(boobie$eigenbands[,1], type = "l")
lines(boobie$eigenbands[,3])
Btru_bsp <- bs(t, df = 12, intercept = TRUE)
Psi <- matrix(0, nrow = p, ncol = p)
for(i in 1:p){
  for(j in 1:p){
    Psi[i, j] <- trapz(t, Btru[, i] * Btru[, j])
  }
}
Psi_sqrt <- pracma::sqrtm(Psi)$B
Psi_sqrt_inv <- solve(Psi_sqrt)
eigenfns <- extract_eigenfn(mcmc_results$Lambda[[1,5000]], mcmc_results$Delta[[1]][,5000], Psi, Psi_sqrt, Psi_sqrt_inv, Btru, 2, c(1,1))
eigenfns$eigenval
pracma::trapz(t, eigenfns[,11]^2)

temp_vec <- eigenfns[,12] / sqrt(sum(eigenfns[,12]^2))
plot(temp_vec, type = "l")
lines(-eigen(truecov)$vectors[,2])
lines(-eigen(mycov)$vectors[,1], col = "blue")
plot(eigenfns[,12], type = "l")
hist(mcmc_results$Nu[[1]])
truecov <- matrix(0, nrow = tmax, ncol = tmax)
truecov <- Btru %*% Lambda1 %*% outer(X[100,], X[100,]) %*% t(Lambda1) %*% t(Btru) +
  Btru %*% Lambda2 %*% outer(X[100,], X[100,]) %*% t(Lambda2) %*% t(Btru)
persp3D(1:50,1:50, truecov)
mycov <- matrix(0, nrow = 12, ncol = 12)
for(i in 1:5000){
  mycov <- mycov + diag(1/mcmc_results$Delta[[1]][,i])
  for(j in 1:K){
    mycov <- mycov + mcmc_results$Lambda[[1,i]][,,j] %*% outer(X[51,], X[51,]) %*% t(mcmc_results$Lambda[[1,i]][,,j]) 
  }
}
mycov <- Btru %*% mycov %*% t(Btru) / 5000
persp3D(1:50, 1:50, mycov)
persp3D(1:50, 1:50, truecov)
persp3D(1:50,1:50,eigen_bands$surface)
eigen(mycov)$values[2] * trapz2est
eigen(mycov)$values[1:3]
trapz1est <- trapz(t, eigen(mycov)$vectors[,1]^2)
trapz2est <- trapz(t, eigen(mycov)$vectors[,2]^2)
eigen(mycov)$values[1] * trapz1est
eigen(mycov)$values[2] * trapz2est
boobie$eigenval_intervals

cumsum(eigen(mycov)$values)/sum(eigen(mycov)$values)
eigenfns$eigenval[1] / eigenfns$eigenval[2]
eigen(mycov)$values[1]/eigen(mycov)$values[2]
t(eigenfns$eigenfn) %*% Psi %*% mycov %*% Psi %*% eigenfns$eigenfn
plot(eigen(mycov)$vectors[,2], type = "l")
lines(eigen(truecov)$vectors[,2], col = "red")
persp3D(1:50,1:50, mycov)
persp3D(1:50,1:50, truecov)
plot(colMeans(Y), type="l")

## To do ##
# 1. Posterior mean coefficient bands
# 2. Posterior eigenfunction bands
# 5. Sparse data





# 
# ### Stan stuff ###
# options(mc.cores =parallel::detectCores())
# rstan_options(auto_write = TRUE)
# m <- stan_model(file = 'stan_code.stan')
# output_samples <- 1000
# my_data <- list(N = n, P = p, D = 1, Q = 2, X = Intercept_only, tmax = tmax, Y = c(t(Y)), B = Btru)
# fit_vb <- vb(m, data = my_data, output_samples = 1000, iter = 30000, tol_rel_obj = .001)
# list_of_draws <- rstan::extract(fit_vb)
# Beta_init <- list_of_draws$Beta[1000,,]
# Lambda_init <- aperm(list_of_draws$Lambda[1000,,,], c(2, 3, 1))
# Eta_init <- list_of_draws$Eta[1000,,]
# Theta_init <- list_of_draws$Theta[1000,,]
# Tau_init <- matrix(1, nrow = Q + 1, ncol = D)
# Delta_init <- list_of_draws$Phi_t[1000,]
# 
# ### Sanity check ###
# plot(Y[1,])
# lines(Best %*%list_of_draws$Theta[500,,1])
# hist(list_of_draws$tau_lambda[,1])
# sum(list_of_draws$log_lik[1000,])
# x <- 1
# covtruth <- Btru%*%Lambda1%*%outer(x,x)%*%t(Lambda1)%*%t(Btru) + Btru%*%Lambda2%*%outer(x,x)%*%t(Lambda2)%*%t(Btru)
# cov_vb <- Btru %*%list_of_draws$Lambda[1000,1,,1] %*% outer(x, x) %*% t(list_of_draws$Lambda[1000,1,,1]) %*% t(Btru) +
#   Btru %*%list_of_draws$Lambda[1000,2,,1] %*% outer(x, x) %*% t(list_of_draws$Lambda[1000,2,,1]) %*% t(Btru)
# persp3D(1:50,1:50, cov_vb)
# persp3D(1:50, 1:50, covtruth)
# covvb <- list
# 
# ### LOO ###
# log_lik_1 <- extract_log_lik(fit_vb, merge_chains = FALSE)
# r_eff <- relative_eff(exp(log_lik_1), cores = 2) 
# loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 2)
# print(loo_1)


lines(Btru %*% AAA$Beta[[1]][,,3000], col="blue")

mymean <- rep(0, 50)
for(i in 1:5000){
  mymean <- mymean + Btru %*% AAA$Beta[[1]][,,i] / 5000
}
lines(mymean, col= "green")
plot(Y[3,])
lines(Btru %*% AAA$Beta[[1]][,,4000] + Btru %*% AAA$Lambda[[1,4000]][,,1] * AAA$Eta[[1]][3,1,4000] + Btru %*% AAA$Lambda[[1,4000]][,,2] * AAA$Eta[[1]][3,2,4000])
