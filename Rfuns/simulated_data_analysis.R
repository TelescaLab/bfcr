library(MASS)
library(plot3D)
library(dlnm)
library(BayesianConditionalFPCA)
library(rstan)
library(loo)
library(tidyverse)
setwd("E:/Rcpp stuff/BayesianConditionalFPCA/Rfuns")
source('simulated_data.r')
matplot(t, t(Y), type = "l", col = "black", xlab = "Time", ylab = "Response")

### MCMC sanity check ###
K <- 4
Basis <- ps(t, df = 16, intercept = TRUE)
mcmc_results <- run_mcmc_Morris(Y, t, X, X, Basis, K, iter = 5000, burnin = 10000, nchains = 1, thin = 5, loglik = 1)

### Visualization ###
sub <- 10
posterior_intervals <- get_posterior_predictive_bands2(mcmc_results, c(.025, .5, .975))
colnames(posterior_intervals) <- c("ID", "Time", "Y", "Lower_P", "Median_P", "Upper_P", "Lower_M", "Median_M", "Upper_M")
posterior_intervals <- as_tibble(posterior_intervals)
posterior_intervals %>%
  filter(ID == sub) %>%
  ggplot(aes(x = Time, y = Y)) +
  geom_point(na.rm = TRUE) +
  geom_ribbon(aes(ymin = Lower_P, ymax = Upper_P), alpha = 0.3) +
  theme_bw()
posterior_intervals %>%
  filter(ID == sub) %>%
  ggplot(aes(x = Time, y = Y)) +
  geom_point(na.rm = TRUE) +
  geom_line(aes(y = Median_M)) +
  theme_bw()
posterior_intervals %>%
  filter(ID == sub) %>%
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
mcmc_omnibus_fit <- get_omnibus_fit2(mcmc_results)
plot(mcmc_omnibus_fit$statistic_obs, mcmc_omnibus_fit$statistic_rep, xlab = "Chi2 observed", ylab = "Chi2 repitition")
abline(a = 0, b = 1)
sum(mcmc_omnibus_fit$statistic_rep > mcmc_omnibus_fit$statistic_obs) / 5000

### LOO ###
r_eff <- relative_eff(exp(mcmc_results$log_lik), cores = 4)
loo_1 <- loo(mcmc_results$log_lik, r_eff = r_eff, cores = 4)
print(loo_1)

### Posterior mean bands ###
sub <- 1
xi <- X[sub,]
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

plot(t, coef_bands$Upper, type = "l", ylim = c(-.2, 1.2))
lines(t, coef_bands$Lower)
lines(t, Btru%*%Theta1 %*% X[sub,],col="green")

### Some covariance visualization ###
sub <- 51
evals <- 2
zi <- X[sub, ]
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

eigen_bands_tibble %>%
  ggplot(aes(x = Time)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = number), alpha = 0.5) +
  geom_line(aes(y = mean)) +
  facet_wrap(number ~., labeller = labeller(number = eig_labs)) +
  theme_bw() +
  theme(legend.position="none")
eigen_bands$magnitude
# aX <- list(title = "Time")
# aY <- list(title = "Time")
# aZ <- list(title = "Response")
# plotly::plot_ly(x = t, y = t, z = eigen_bands$surface, type = "surface") %>%
#   plotly::layout(scene = list(xaxis = aX, yaxis = aY, zaxis = aZ, dragmode = "turntable"))
# plot(eigen_bands$raw_magnitude, type = "l")
truecov <- matrix(0, nrow = tmax, ncol = tmax)
truecov <- Btru %*% Lambda1 %*% outer(X[sub,], X[sub,]) %*% t(Lambda1) %*% t(Btru) +
  Btru %*% Lambda2 %*% outer(X[sub,], X[sub,]) %*% t(Lambda2) %*% t(Btru)
par(mfrow = c(1,2))
persp3D(t, t, eigen_bands$surface)
persp3D(t, t, truecov)

