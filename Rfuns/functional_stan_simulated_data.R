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
setwd("E:/Rcpp stuff/BayesianConditionalFPCA/Rfuns")
options(mc.cores =parallel::detectCores())
rstan_options(auto_write = TRUE)
m <- stan_model(file = 'stan_code.stan')

source('simulated_data.r')
output_samples <- 1000

my_data <- list(N = n, P = p, D = 1, Q = 2, X = Intercept_only, tmax = tmax, Y = c(t(Y)), B = Btru)
fit_vb <- vb(m, data = my_data, output_samples = 1000, iter = 30000, tol_rel_obj = .001)
list_of_draws <- rstan::extract(fit_vb)
Beta_init <- list_of_draws$Beta[1000,,]
Lambda_init <- aperm(list_of_draws$Lambda[1000,,,], c(2, 3, 1))
Eta_init <- list_of_draws$Eta[1000,,]
Theta_init <- list_of_draws$Theta[1000,,]
Tau_init <- matrix(1, nrow = Q + 1, ncol = D)
Delta_init <- list_of_draws$Phi_t[1000,]

### Sanity check ###
plot(Y[1,])
lines(Best %*%list_of_draws$Theta[500,,1])
hist(list_of_draws$tau_lambda[,1])
sum(list_of_draws$log_lik[1000,])
x <- 1
covtruth <- Btru%*%Lambda1%*%outer(x,x)%*%t(Lambda1)%*%t(Btru) + Btru%*%Lambda2%*%outer(x,x)%*%t(Lambda2)%*%t(Btru)
cov_vb <- Btru %*%list_of_draws$Lambda[1000,1,,1] %*% outer(x, x) %*% t(list_of_draws$Lambda[1000,1,,1]) %*% t(Btru) +
  Btru %*%list_of_draws$Lambda[1000,2,,1] %*% outer(x, x) %*% t(list_of_draws$Lambda[1000,2,,1]) %*% t(Btru)
persp3D(1:50,1:50, cov_vb)
persp3D(1:50, 1:50, covtruth)
covvb <- list
### LOO ###
log_lik_1 <- extract_log_lik(fit_vb, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik_1), cores = 2) 
loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 2)
print(loo_1)

### MCMC sanity check ###
K <- 2
AAA <- MCMC(Y, X, Btru, K, 5000, 5000, 1, 1)
P <- posterior_predictive_bands(AAA, c(.05, .5, .95))
colnames(P) <- c("ID", "Time", "Y", "Lower_P", "Median_P", "Upper_P", "Lower_M", "Median_M", "Upper_M")
P <- as_tibble(P)
P %>%
  filter(ID == 7) %>%
  ggplot(aes(x = Time, y = Y)) +
  geom_point() +
  geom_ribbon(aes(ymin = Lower_P, ymax = Upper_P), alpha = 0.3) +
  theme_bw()
P %>%
  group_by(ID) %>%
  filter(Y > Lower_P & Y < Upper_P) %>%
  summarize(coverage = n() / tmax)
P %>%
  filter(Y > Lower_P & Y < Upper_P) %>%
  summarize(coverage = n() / (tmax * n))
A <- omnibus_fit(AAA)
sum(A$statistic_rep > A$statistic_obs) / 5000
truecov <- matrix(0, nrow = tmax, ncol = tmax)
truecov <- Btru %*% Lambda1 %*% outer(X[1,], X[1,]) %*% t(Lambda1) %*% t(Btru) +
  Btru %*% Lambda2 %*% outer(X[1,], X[1,]) %*% t(Lambda2) %*% t(Btru)

mycov <- matrix(0, nrow = 12, ncol = 12)
for(i in 1:5000){
  mycov <- mycov + diag(1/AAA$Delta[[1]][,i]) 
  # mycov <- matrix(0, nrow = tmax, ncol = tmax)
  for(j in 1:K){
    mycov <- mycov + AAA$Lambda[[1,i]][,,j] %*% outer(X[1,], X[1,]) %*% t(AAA$Lambda[[1,i]][,,j])
  }
}
mycov <- mycov / 5000
mycov <- Btru %*% mycov %*% t(Btru)
# plot(eigen(mycov)$vectors[,1], type = "l")
persp3D(1:50,1:50, mycov)
persp3D(1:50,1:50, truecov)
plot(colMeans(Y), type="l")
lines(Btru %*% AAA$Beta[[1]][,,3000], col="blue")

mymean <- rep(0, 50)
for(i in 1:5000){
  mymean <- mymean + Btru %*% AAA$Beta[[1]][,,i] / 5000
}
lines(mymean, col= "green")
plot(Y[3,])
lines(Btru %*% AAA$Beta[[1]][,,4000] + Btru %*% AAA$Lambda[[1,4000]][,,1] * AAA$Eta[[1]][3,1,4000] + Btru %*% AAA$Lambda[[1,4000]][,,2] * AAA$Eta[[1]][3,2,4000])
