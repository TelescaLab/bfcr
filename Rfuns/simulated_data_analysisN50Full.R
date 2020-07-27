args = commandArgs(trailingOnly = TRUE)
# setwd("/u/home/j/jshamsh1/Documents/BayesianConditionalFPCA_results")
# files_in_dir <- dir()
# numbers_to_run <- which(!(paste0("N100Full", 1:500, ".RData") %in% files_in_dir))
# seed <- numbers_to_run[as.numeric(args)]
seed <- as.numeric(args)
set.seed(seed)
library(splines)
library(MASS)
# library(plot3D)
library(dlnm)
library(BayesianConditionalFPCA)
# library(rstan)
# library(loo)
library(tidyverse)
library(pracma)
library(mgcv)

#####################################################################
################### MEAN AND VARIANCE FUNCTIONS #####################
#####################################################################
n <- 50
tmax <- 100
t <- seq(from = 0, to = 1, length.out = tmax)
p <- 12
B <- bs(t, df = p, intercept = TRUE)
int_function <- 1*(t-.5*t^2)
x_function <- t*sin(2*pi*t)
#x_function <- .2*sqrt(t)
lbasis1_int <- .2*(1-1*(t-.5)^2)
#lbasis1_int <- rep(0, length(t))
#lbasis1_int <- .25*cos(t*pi)
lbasis2_int <- .25*sin(t*pi)
lbasis1_z <- .25*exp(cos(t))
lbasis2_z <- .25*exp(sin(t))
#lbasis1_z <- rep(0, length(t))
#lbasis2_z <- rep(0, length(t))
#####################################################################
####################### DATA GENERATION #############################
#####################################################################
X <- cbind(rep(1,n), rnorm(n, sd = 1))
Age <- X[,2]
AgeBasis <- smooth.construct(s(Age, bs = "ps", k = 10), data.frame(Age), NULL)
noise_sd <- .1
# E <- matrix(rnorm(tmax * n,sd=noise_sd), nrow = n, ncol = tmax)
E <- matrix(0, nrow = n, ncol = tmax)
for(i in 1:n){
  # E[i, ] <- rnorm(tmax, mean = 0, sd = noise_sd + (i %% 5) * noise_sd / 5)
  E[i, ] <- rnorm(tmax, mean = 0, sd = noise_sd)
}
mean_fn <- cbind(int_function, x_function)
mean_grid <- cbind(int_function, x_function) %*% t(X)
basis1_fn <- cbind(lbasis1_int, lbasis1_z)
basis2_fn <- cbind(lbasis2_int, lbasis2_z)
basis1_grid <- cbind(lbasis1_int, lbasis1_z) %*% t(X)
basis2_grid <- cbind(lbasis2_int, lbasis2_z) %*% t(X)

Eta1 <- rnorm(n, sd = sqrt(1))
Eta2 <- rnorm(n, sd = sqrt(1))
Y_mean <- t(mean_grid)
Y_no_error <-  t(mean_grid) + diag(Eta1)%*%t(basis1_grid) + diag(Eta2)%*%t(basis2_grid)
Y <- Y_no_error + E

#####################################################################
########################## MCMC #####################################
#####################################################################
K <- 4
Basis <- ps(t, df = 16, intercept = TRUE)
mcmc_results <- run_mcmc_Morris(Y, t, AgeBasis$X, cbind(rep(1,n)), Basis, K = 2, iter = 20000, burnin = 20000, nchains = 1, thin = 1, loglik = 0)
results <- numeric(21)
#####################################################################
####################### VISUALIZATION ###############################
#####################################################################
# sub <- 4
posterior_intervals <- get_posterior_predictive_bands2(mcmc_results, c(.025, .5, .975))
colnames(posterior_intervals) <- c("ID", "Time", "Y", "Lower_P", "Median_P", "Upper_P", "Lower_M", "Median_M", "Upper_M")
posterior_intervals <- as_tibble(posterior_intervals)
posterior_intervals$Y_no_error <- c(t(Y_no_error))
# posterior_intervals %>%
#   filter(ID == sub) %>%
#   ggplot(aes(x = Time, y = Y_no_error)) +
#   geom_point(na.rm = TRUE) +
#   geom_ribbon(aes(ymin = Lower_M, ymax = Upper_M), alpha = 0.3) +
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

#####################################################################
####################### RECORD RESULTS ##############################
#####################################################################
results[1] <- unlist(posterior_intervals %>%
                       group_by(ID) %>%
                       summarize(mymean = 100 * trapz(t,(Median_P - Y_no_error)^2) / 
                                   trapz(t, (Y_no_error)^2)) %>%
                       ungroup() %>%
                       summarize(mmean = mean(mymean)))

results[2] <- 100 * unlist(posterior_intervals %>%
                             filter(Y_no_error > Lower_M & Y_no_error < Upper_M) %>%
                             summarize(coverage = n() / (n * tmax)))

results[3] <- 100 * unlist(posterior_intervals %>%
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
x_seq <- seq(from = -2, to = 2, length.out = 10)
for(i in 1:length(x_seq)){
  xi <- c(0, x_seq[i])
  coef_bands <- get_posterior_means(mcmc_results, xi, alpha)
  coef_bands <- cbind(t, coef_bands)
  colnames(coef_bands) <- c("Time", "Lower", "Mean", "Upper")
  coef_bands <- as_tibble(coef_bands)
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
  # lines(t, mean_fn %*% xi,col="green")
  results[6] <- trapz(t, (coef_bands$Mean - mean_fn %*% xi)^2) /
    trapz(t, (mean_fn %*% xi)^2) * 100 + results[6]
  results[7] <- mean(as.numeric(mean_fn %*% xi >= coef_bands$Lower & 
                                  mean_fn %*% xi <= coef_bands$Upper)) * 100 +
    results[7]
  results[8] <- mean(coef_bands$Upper - coef_bands$Lower) + results[8]
}
results[6:8] <- results[6:8] / length(x_seq)
coef_bands <- get_posterior_coefs(mcmc_results, alpha)

# Intercept_fn <- Btru %*% Theta1 %*% c(1, 0)
# Age_fn <- Btru %*% Theta1 %*% c(0, 1)
coef_bands <- tibble(Frequency = rep(t, times = dim(X)[2]),
                     Covariate = rep(c("Intercept", "Age"), each = dim(TBasis$X)[1]),
                     Lower = c(coef_bands$lower),
                     Mean = c(coef_bands$mean),
                     Upper = c(coef_bands$upper))
Intercept_idx <- which(coef_bands$Covariate == "Intercept")
results[9] <- 100 * trapz(t, (coef_bands$Mean[Intercept_idx] - int_function)^2) / 
  trapz(t, int_function^2)
results[10] <- 100 * trapz(t, (coef_bands$Mean[-Intercept_idx] - x_function)^2) /
  trapz(t, x_function^2)
results[11] <- mean(as.numeric((int_function >= coef_bands$Lower[Intercept_idx] &
                                  int_function <= coef_bands$Upper[Intercept_idx]))) * 100
results[12] <- mean(as.numeric((x_function >= coef_bands$Lower[-Intercept_idx] &
                                  x_function <= coef_bands$Upper[-Intercept_idx]))) *
  100
results[13] <- mean(coef_bands$Upper[Intercept_idx] - coef_bands$Lower[Intercept_idx])
results[14] <- mean(coef_bands$Upper[-Intercept_idx] - coef_bands$Lower[-Intercept_idx])
# coef_bands %>%
#   ggplot(aes(x = Frequency)) +
#   geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Covariate), alpha = 0.5) +
#   geom_line(aes(y = Mean)) +
#   geom_hline(yintercept = 0) +
#   facet_wrap(Covariate ~., scales = "free") +
#   ylab("Power") +
#   theme_minimal()

### Some covariance visualization ###
mag <- numeric(length(z_seq))
L <- matrix(rnorm(16*5), ncol = 5, nrow = 16)
Bz <- ps(z_seq, df = 5, intercept = TRUE)
for(i in 1:length(z_seq)){
  z <- Bz[i,]
  mag[i] <- Re(sum(eigen(L%*%outer(z, z)%*%t(L))$values))
}
plot(mag, type = "l")
z_seq <- seq(from = -3, to = 3, length.out = 20)
evals <- 2
true_mag <- numeric(length(z_seq))
est_mag <- numeric(length(z_seq))
true_mag <- matrix(0, nrow = 50, ncol = length(z_seq))
est_mag <- matrix(0, nrow = 50, ncol = length(z_seq))
for(i in 1:length(z_seq)){
  zi <- c(1, z_seq[i])
  truecov <- matrix(0, nrow = tmax, ncol = tmax)
  #truecov <- basis1_fn %*% outer(zi, zi) %*% t(basis1_fn) +
  #  basis2_fn %*% outer(zi, zi) %*% t(basis2_fn)
  truecov <- tcrossprod(basis1_fn[,1]) + tcrossprod(basis2_fn[,1])
  e1 <- eigen(truecov)$vectors
  m1 <- sapply(1:2, function(i) trapz(t, e1[,i]^2))
  true_mag[i] <- sum(eigen(truecov)$values * m1)
  #true_mag[,i] <- diag(truecov)
  est_mag[i] <- get_posterior_eigen2(mcmc_results, evals, zi, alpha)$magnitude[2]
  #est_mag[,i] <- diag(get_posterior_eigen2(mcmc_results, evals, zi, alpha)$surface)
}
plot(est_mag[,15], type = "l")
lines(true_mag[,15])
estcov <- get_posterior_eigen2(mcmc_results, evals, c(1,-3), alpha)$surface
persp3D(1:50,1:50, estcov)
persp3D(1:50,1:50, truecov)
plot(true_mag, type = "l")
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
  truecov <- basis1_fn %*% outer(zi, zi) %*% t(basis1_fn) +
    basis2_fn %*% outer(zi, zi) %*% t(basis2_fn)
  # truecov <- Btru %*% Lambda1 %*% outer(zi, zi) %*% t(Lambda1) %*% t(Btru) +
  #  Btru %*% Lambda2 %*% outer(zi, zi) %*% t(Lambda2) %*% t(Btru)
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
  
  
  results[15:16] <- 100 * unlist((eigen_bands_tibble %>%
                                    group_by(number) %>%
                                    summarize(RISE = trapz(t, (mean - eigvec)^2) / trapz(t, eigvec^2)))[,2]) +
    results[15:16]
  
  results[17:18] <- unlist((eigen_bands_tibble %>%
                              group_by(number) %>%
                              filter(eigvec >= lower & eigvec <= upper) %>%
                              summarize(coverage = 100 * n() / tmax))[,2]) + results[17:18]
  
  results[19:20] <- unlist((eigen_bands_tibble %>%
                              group_by(number) %>%
                              summarize(mean_l = mean(upper - lower)))[,2]) + results[19:20]
  results[21] <- trapz(t,sapply(1:tmax, function(i) trapz(t, (eigen_bands$surface[,i] -
                                                                truecov[,i])^2))) / trapz(t, sapply(1:tmax, function(i) trapz(t,
                                                                                                                              truecov[,i]^2))) * 100 + results[21]
  print(trapz(t,sapply(1:tmax, function(i) trapz(t, (eigen_bands$surface[,i] -
                                                       truecov[,i])^2))) / trapz(t, sapply(1:tmax, function(i) trapz(t,
                                                                                                                     truecov[,i]^2))))
  # results[22] <- trapz(t,sapply(1:tmax, function(i) trapz(t, (eigen_bands$surface_cor[,i] -
  #                                                               cov2cor(truecov)[,i])^2))) / trapz(t, sapply(1:tmax, function(i) trapz(t,
  #                                                                                                                             cov2cor(truecov)[,i]^2))) * 100 + results[22]
  
}

results[15:21] <- results[15:21] / length(z_seq)
results
setwd("/u/home/j/jshamsh1/Documents/BayesianConditionalFPCA_results")
print("Changed directory")
save(results, file = paste0("N100Full", seed, ".RData"))
# 
# estcov <- matrix(0, nrow = 50, ncol = 50)
# for(i in 1:20000){
#   for(k in 1:4){
#     estcov <- estcov + mcmc_results$B %*% apple$Lambda[,,k] %*% outer(zi, zi) %*% t(apple$Lambda[,,k]) %*% t(mcmc_results$B)
#   }
# }
# estcov <- estcov / 20000
# persp3D(1:50,1:50, eigen_bands$surface)
# persp3D(1:50,1:50, truecov)
# persp3D(1:50, 1:50, estcov)
# save(results, file = paste0("/Users/johnshamshoian/Rcpp/BayesianConditionalFPCA/simulation/N50Full", seed,".RData"))

# plot(eigen_bands_tibble$eigvec[1:50], type = "l")
# lines(eigen_bands_tibble$mean[1:50], col = "green")
# lines(eigen_bands_tibble$lower[1:50], col = "blue")
# lines(eigen_bands_tibble$upper[1:50], col = "blue")
# plot(eigen_bands_tibble$eigvec[51:100], type = "l")
# lines(eigen_bands_tibble$mean[51:100], col = "green")
# lines(eigen_bands_tibble$lower[51:100], col = "blue")
# lines(eigen_bands_tibble$upper[51:100], col = "blue")
# plot(eigen_bands$raw_magnitude, type = "l")

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
library(spam)
Age <- X[,2]
TBasis <- smooth.construct(s(t, bs = "ps", k = 12), data.frame(t), NULL)
AgeBasis <- smooth.construct(s(Age, bs = "ps", k = 4), data.frame(Age), NULL)
DT <- as.matrix(diff.spam(diag.spam(12), lag = 1, differences = 2))
DAge <- as.matrix(diff.spam(diag.spam(4), lag = 1, differences = 2))
P1 <- t(DT)%*%DT
P2 <- t(DAge) %*% DAge
PP1 <- kronecker(t(DAge)%*%DAge, diag(12))
PP2 <- kronecker(diag(4), t(DT)%*%DT)
Q <- precmat.IGMRFreglat(12, 4, order = 2)

myte_basis<- smooth.construct.tensor.smooth.spec(te(t, Age, bs = "ps", k = c(12,4)), data.frame(t=t, Age=Age), knots = NULL)

matlist <- list()
matlist[[1]] <- PP1
matlist[[2]] <- PP2
matlist[[1]] <- kronecker(myte_basis$XP[[2]] %*% t(DAge)%*%DAge %*%t(myte_basis$XP[[2]]), diag(12))
matlist[[2]] <- kronecker(diag(4), myte_basis$XP[[1]]%*%t(DT)%*%DT%*%t(myte_basis$XP[[1]]))
matlist[[1]] <- as.matrix(Q)
matlistmean <- list()
matlistmean[[1]] <- as.matrix(Q)
matlistvar <- list()
matlistvar[[1]] <- as.matrix(Q)
matlistmean[[1]] <- kron(TBasis$S[[1]], diag(dim(AgeBasis$X)[2]))
matlistmean[[2]] <- kron(diag(dim(TBasis$X)[2]),AgeBasis$S[[1]])
matlistvar[[1]] <- kron(TBasis$S[[1]], diag(dim(AgeBasis$X)[2]))
matlistvar[[2]] <- kron(diag(dim(TBasis$X)[2]),AgeBasis$S[[1]])
matlist[[1]] <- kronecker(TBasis$S[[1]], AgeBasis$S[[1]])
for(i in 1:4){
  matlist[[i]] <- TBasis$S[[1]]
}
mcmc_results <- run_mcmc_Morris_Tensor(Y, t, AgeBasis$X, AgeBasis$X, TBasis$X,
                                       matlist, matlist, c(1), c(1), K = 4, iter = 20000,
                                       burnin = 20000, nchains = 1,
                                       thin = 1, loglik = 0)

matlistmean <- list()
matlistmean[[1]] <- myte_basis$S[[1]]
matlistmean[[2]] <- myte_basis$S[[2]]
matlistvar <- list()
matlistvar[[1]] <- TBasis$S[[1]]
mcmc_results <- run_mcmc_Morris_Tensor(Y, t, AgeBasis$X, cbind(rep(1,n)),
                                       TBasis$X, matlistmean, matlistvar, c(1,1),
                                       c(1), K = 2, iter = 20000, burnin = 20000,
                                       nchains = 1, thin = 1, loglik = 0)
results <- numeric(21)
matlist <- list()
matlist[[1]] <- TBasis$S[[1]]
matlist[[2]] <- TBasis$S[[1]]
mcmc_results <- run_mcmc_Morris_Tensor(Y, t, X, X, TBasis$X,
                                       matlist, matlist, c(1,2), c(1,2),
                                       K = 2, iter = 20000,
                                       burnin = 20000, nchains = 1,
                                       thin = 1, loglik = 0)





alpha <- 0.05
Age_seq <- seq(from = -3, to = 3, length.out = 20)
Age_df <- data.frame(Age = Age_seq)
new_points <- Predict.matrix(AgeBasis, Age_df)
for(i in 1:length(Age_seq)){
  xi <- new_points[i,]
  #xi <- c(1, Age_seq[i])
  xim <- c(1, Age_seq[i])
  coef_bands <- get_posterior_means(mcmc_results, xi, alpha)
  coef_bands <- cbind(t, coef_bands)
  colnames(coef_bands) <- c("Time", "Lower", "Mean", "Upper")
  coef_bands <- as_tibble(coef_bands)
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
  # lines(t, mean_fn %*% xi,col="green")
  results[6] <- trapz(t, (coef_bands$Mean - mean_fn %*% xim)^2) /
    trapz(t, (mean_fn %*% xim)^2) * 100 + results[6]
  results[7] <- mean(as.numeric(mean_fn %*% xim >= coef_bands$Lower & 
                                  mean_fn %*% xim <= coef_bands$Upper)) * 100 +
    results[7]
  results[8] <- mean(coef_bands$Upper - coef_bands$Lower) + results[8]
}
results[6:8] <- results[6:8] / length(Age_seq)

idx <- 1
xi <- new_points[idx,]
#xi <- c(1, Age_seq[idx])
coef_bands <- get_posterior_means(mcmc_results, xi, alpha)
plot(coef_bands[,2], type = "l")
lines(mean_fn%*%c(1, Age_seq[idx]), col = "red")
lines(coef_bands[,1], col = "blue", lty = 2)
lines(coef_bands[,3], col = "blue", lty = 2)

estmeanmat <- matrix(0, nrow = length(Age_seq), ncol = length(t))
truemeanmat <- matrix(0, nrow = length(Age_seq), ncol = length(t))
for(i in 1:length(Age_seq)){
  xi <- new_points[i, ]
  coef_bands <- get_posterior_means(mcmc_results, xi, alpha)
  estmeanmat[i,] <- coef_bands[,2]
  truemeanmat[i,] <- mean_fn%*%c(1, Age_seq[i])
}
par(mfrow = c(1,2))
persp3D(1:length(Age_seq), 1:length(t), estmeanmat)
persp3D(1:length(Age_seq), 1:length(t), truemeanmat)

for(i in 1:length(Age_seq)){
  evals <- 2
  zi <- new_points[i,]
  zim <- c(1, Age_seq[i])
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
  
  #truecov <- matrix(0, nrow = tmax, ncol = tmax)
  #truecov <- tcrossprod(basis1_fn[,1]) + tcrossprod(basis2_fn[,1])
  truecov <- basis1_fn %*% outer(zim, zim) %*% t(basis1_fn) +
    basis2_fn %*% outer(zim, zim) %*% t(basis2_fn)
  #truecov <- Btru %*% Lambda1 %*% outer(zi, zi) %*% t(Lambda1) %*% t(Btru) +
  #  Btru %*% Lambda2 %*% outer(zi, zi) %*% t(Lambda2) %*% t(Btru)
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
  
  
  results[15:16] <- 100 * unlist((eigen_bands_tibble %>%
                                    group_by(number) %>%
                                    summarize(RISE = trapz(t, (mean - eigvec)^2) / trapz(t, eigvec^2)))[,2]) +
    results[15:16]
  
  results[17:18] <- unlist((eigen_bands_tibble %>%
                              group_by(number) %>%
                              filter(eigvec >= lower & eigvec <= upper) %>%
                              summarize(coverage = 100 * n() / tmax))[,2]) + results[17:18]
  
  results[19:20] <- unlist((eigen_bands_tibble %>%
                              group_by(number) %>%
                              summarize(mean_l = mean(upper - lower)))[,2]) + results[19:20]
  results[21] <- trapz(t,sapply(1:tmax, function(i) trapz(t, (eigen_bands$surface[,i] -
                                                                truecov[,i])^2))) / trapz(t, sapply(1:tmax, function(i) trapz(t,
                                                                                                                              truecov[,i]^2))) * 100 + results[21]
  print(100 * trapz(t,sapply(1:tmax, function(i) trapz(t, (eigen_bands$surface[,i] -
                                                       truecov[,i])^2))) / trapz(t, sapply(1:tmax, function(i) trapz(t,
                                                                                                                     truecov[,i]^2)))) 
  # results[22] <- trapz(t,sapply(1:tmax, function(i) trapz(t, (eigen_bands$surface_cor[,i] -
  #                                                               cov2cor(truecov)[,i])^2))) / trapz(t, sapply(1:tmax, function(i) trapz(t,
  #                                                                                                                             cov2cor(truecov)[,i]^2))) * 100 + results[22]
  
}
results[15:21] <- results[15:21] / length(Age_seq)

idx <- 11
evals <- 2
zi <- new_points[idx,]
zim <- c(1, Age_seq[idx])
#zi <- c(1, Age_seq[10])
#truecov <- tcrossprod(lbasis1_int) + tcrossprod(lbasis2_int)
truecov <- basis1_fn %*% outer(zim, zim) %*% t(basis1_fn) +
  basis2_fn %*% outer(zim, zim) %*% t(basis2_fn)
eigen_bands <- get_posterior_eigen2(mcmc_results, evals, zi, alpha)
persp3D(1:50,1:50, truecov)
persp3D(1:50,1:50, eigen_bands$surface)

mag <- numeric(length(Age_seq))
for(i in 1:length(Age_seq)){
  zi <- new_points[i,]
  mag[i] <- sum(eigen(get_posterior_eigen2(mcmc_results, evals, zi, alpha)$surface)$values)
}
truecov <- tcrossprod(basis1_fn[,1]) + tcrossprod(basis2_fn[,1])
persp3D(1:50,1:50, truecov)
persp3D(1:50,1:50, eigen_bands$surface)






L <- diag(5)
Om <- diag(c(10000,10000,1,1,1))
L[2,1] <- -2
L[3,2] <- -2
L[4,3] <- -2
L[5,4] <- -2

L[3,1] <- 1
L[4,2] <- 1
L[5,3] <- 1
L <- diag(4)
L[2,1] <- -1
L[3,2] <- -1
L[4,3] <- -1
Om <- diag(c(10000,1,1,1))
t(L)%*%solve(Om)%*%L
L <- diag(20)
for(i in 1:19){
  L[i + 1, i] <- -1
}
for(i in 1:15){
  L[i + 5, i] <- -1
}
Om <- diag(c(10000, rep(1,19)))
t(L)%*%solve(Om)%*%L
L[2,1] <- -1  
L[6,1] <- -1
L[3,2]
Q <- diag(c(2,3,3,3))


Q <- precmat(5, order = 2)
Q <- precmat.IGMRFreglat(16, 4)
