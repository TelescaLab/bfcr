for(seed in c(51, 371, 385)){
# setwd("/u/home/j/jshamsh1/Documents/BayesianConditionalFPCA_results")
# files_in_dir <- dir()
# numbers_to_run <- which(!(paste0("N100Full", 1:500, ".RData") %in% files_in_dir))
# seed <- numbers_to_run[as.numeric(args)]
# seed <- as.numeric(args)
print(paste0("seed is ", seed))
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


#####################################################################
################### MEAN AND VARIANCE FUNCTIONS #####################
#####################################################################
n <- 100
tmax <- 50
t <- seq(from = 0, to = 1, length.out = tmax)
p <- 12
B <- bs(t, df = p, intercept = TRUE)
int_function <- 1*(t-.5*t^2)
x_function <- t*sin(2*pi*t)
lbasis1_int <- .5*cos(t*pi)
lbasis2_int <- .5*sin(t*pi)
#lbasis1_z <- .25*exp(cos(t))
#lbasis2_z <- .25*exp(sin(t))

#####################################################################
####################### DATA GENERATION #############################
#####################################################################
X <- cbind(rep(1,n), rnorm(n, sd = 1))
noise_sd <- .2
# E <- matrix(rnorm(tmax * n,sd=noise_sd), nrow = n, ncol = tmax)
E <- matrix(0, nrow = n, ncol = tmax)
for(i in 1:n){
  # E[i, ] <- rnorm(tmax, mean = 0, sd = noise_sd + (i %% 5) * noise_sd / 5)
  E[i, ] <- rnorm(tmax, mean = 0, sd = noise_sd)
}
mean_fn <- cbind(int_function, x_function)
mean_grid <- cbind(int_function, x_function) %*% t(X)
# basis1_fn <- cbind(lbasis1_int, lbasis1_z)
# basis2_fn <- cbind(lbasis2_int, lbasis2_z)
basis1_fn <- cbind(lbasis1_int)
basis2_fn <- cbind(lbasis2_int)
# basis1_grid <- cbind(lbasis1_int, lbasis1_z) %*% t(X)
# basis2_grid <- cbind(lbasis2_int, lbasis2_z) %*% t(X)
basis1_grid <- basis1_fn %*% t(rep(1,n))
basis2_grid <- basis2_fn %*% t(rep(1,n))
Eta1 <- rnorm(n, sd = sqrt(1))
Eta2 <- rnorm(n, sd = sqrt(1))
Y_no_error <- t(mean_grid) + diag(Eta1)%*%t(basis1_grid) + diag(Eta2)%*%t(basis2_grid)
Y <- Y_no_error + E


#####################################################################
########################## MCMC #####################################
#####################################################################
K <- 4
Basis <- ps(t, df = 16, intercept = TRUE)
mcmc_results <- run_mcmc_Morris(Y, t, X, X, Basis, K, iter = 20000, burnin = 20000, nchains = 1, thin = 1, loglik = 0)
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
                     Covariate = rep(c("Intercept", "Age"), each = dim(Basis)[1]),
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
  truecov <- basis1_fn %*% outer(1, 1) %*% t(basis1_fn) +
    basis2_fn %*% outer(1, 1) %*% t(basis2_fn)
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
mag <- numeric(30)
z_seq <- seq(from = -3, to = 3, length.out = 30)
estcov <- matrix(0, nrow = 50, ncol = 50)
for(i in 1:length(z_seq)){
  print(i)
  zi <- c(1, z_seq[i])
  mag[i] <- median(get_posterior_eigen2(mcmc_results, evals, zi, alpha)$raw_magnitude)
}
estcov <- estcov / length(z_seq)

trapz(t,sapply(1:tmax, function(i) trapz(t, (estcov[,i] -
                                               truecov[,i])^2))) / trapz(t, sapply(1:tmax, function(i) trapz(t,
                                                                                                             truecov[,i]^2))) * 100
persp3D(1:50,1:50, eigen_bands <- get_posterior_eigen2(mcmc_results, evals, c(1,-2), alpha)$surface)
persp3D(1:50,1:50, truecov)
results[15:21] <- results[15:21] / length(z_seq)

rawmag_2 <- get_posterior_eigen2(mcmc_results, evals, c(1,-2), alpha)$raw_magnitude
results
setwd("/Users/johnshamshoian/Rcpp/results")
print("Changed directory")
save(results, file = paste0("N50NoZFitZ", seed, ".RData"))
}
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




