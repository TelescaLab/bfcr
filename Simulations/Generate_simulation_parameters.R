# This script is used to generate data generating parameters used in the 
# simulation study.
library(plotly)
library(mgcv)
set.seed(8)
n <- 100
k <- 4
df_times <- 10
df_x <- 6
times <- seq(from = 0, to = 1, length.out = 100)
x <- seq(from = 0, to = 1, length.out = n)
tau_smooth_time <- 10
tau_smooth_x <- 100
delta_baseline <- 2
delta_x <- 4
tau_shrink_baseline <- cumprod(rep(delta_baseline, k))
tau_shrink_x <- cumprod(rep(delta_x, k))
time_basis <- smoothCon(s(times, k = df_times, bs = "ps", m = 2),
                        data.frame(times))
x_basis <- smoothCon(s(x, k = df_x, bs = "ps", m = 2), absorb.cons = TRUE,
                     data.frame(x))
time_penalty <- time_basis[[1]]$S.scale * time_basis[[1]]$S[[1]]
x_penalty <- x_basis[[1]]$S.scale * x_basis[[1]]$S[[1]]
S_list <- list(x_penalty, time_penalty)
S <- tensor.prod.penalties(S_list)
precision_beta <- matrix(0, nrow = df_times * df_x, df_times * df_x)
precision_beta[1:df_times, 1:df_times] <- time_penalty
precision_beta[-c(1:df_times), -c(1:df_times)] <- S[[1]] + S[[2]]
precision_beta <- precision_beta + .1 * diag(df_times * df_x)
beta <- MASS::mvrnorm(n=1, mu = rep(0, df_x * df_times),
                      Sigma = solve(precision_beta)) %>%
  matrix(nrow = df_times)
lambda <- array(0, dim = c(df_times, df_x, k))
for (kp in 1:k) {
  precision_lambda <- matrix(0, nrow = df_times * df_x, df_times * df_x)
  precision_lambda[1:df_times, 1:df_times] <- time_penalty + tau_shrink_baseline * diag(df_times)
  precision_lambda[-c(1:df_times), -c(1:df_times)] <- S[[1]] + S[[2]] + tau_shrink_x * diag((df_x - 1) * df_times)
  lambda[,,kp] <- MASS::mvrnorm(n = 1, mu = rep(0, df_x * df_times),
                               Sigma = solve(precision_lambda)) %>%
    matrix(nrow = df_times)
}
eta <- rnorm(n * k) %>% array(dim = c(n, k))
X <- cbind(1, x_basis[[1]]$X)
Y <- time_basis[[1]]$X %*% beta %*% t(X)

# Plot mean curves
matplot(Y, type = "l")
for (kp in 1:k) {
  Y <- Y + time_basis[[1]]$X %*% lambda[,,kp] %*% t(X) %*% diag(eta[,kp])
}

# Plot mean + random deviations from mean
matplot(Y, type = "l")
simparam <- list(beta = beta, lambda = lambda, X = X, k = k, n = n, 
                 df_times = df_times, df_x = df_x, x = x,
                 tau_smooth_time = tau_smooth_time, 
                 delta_base = delta_baseline, delta_x = delta_x,
                 time_basis = time_basis, x_basis = x_basis, times = times)
file_name <- paste0("/Users/johnshamshoian/Documents/R_projects/bfcr/",
                    "Simulations/Simulation_parameters.RData")
save(simparam, file = file_name)
covfs <- array(0, dim = c(length(times), length(times), n))
for (i in 1:n) {
  for (kp in 1:k) {
    covfs[,,i] <- covfs[,,i] + tcrossprod(time_basis[[1]]$X %*% lambda[,,kp] %*% X[i,])
  }
}
covmat <- matrix(0, n, 10)
for (xx in 1:10) {
  covmat[,xx] <- covfs[xx * 10, xx* 10,]
} 

# Plot several elements of within-subject covariance matrix, as a function of
# a covariate
matplot(covmat, type = "l")

# Plot low dimensional g summaries for the intercept and covariate effect
intercept_effect <- numeric(length(times))
x_effect <- numeric(length(times))
for (kp in 1:k) {
  intercept_effect <- intercept_effect + 
    (time_basis[[1]]$X %*% lambda[,1,kp])^2
}
for (kp in 1:k) {
  for (i in 1:10) {
    x_effect <- x_effect + 
      (time_basis[[1]]$X %*% lambda[,-1,kp] %*% x_basis[[1]]$X[i,])^2
  }
}
x_effect <- x_effect / 10
plot(intercept_effect, type = "l", ylim = c(0, max(intercept_effect + .1)))
lines(x_effect, type = "l")



