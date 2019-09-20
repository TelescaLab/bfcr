basis_version <- "
  data{
    int<lower=0> N; // Number of samples
    int<lower=0> P; // Dimension of basis
    int<lower=0> Q; // Dimension of loading matrix
    matrix[N,2] X;
    int<lower=0> tmax;
    vector[tmax * N] Y;
    matrix[tmax, P] B;
  }
  parameters{
    vector[P * N] Theta;
    matrix[P, Q] Lambda; // Loading matrix
    vector<lower=0>[P] Phi_t;
    matrix[2, Q] Beta; // Latent mean matrix
    matrix[2, Q] D; // Latent covariance matrix
    vector[N] Xi; //Latent cov random effect
    real<lower=0> sigma; // error variance
  }
  transformed parameters{
    vector[P * N] Theta_hat;
    vector[N * Q] Eta;
    for(i in 1:N){
      Eta[((i-1)*Q + 1):(i*Q)] = Beta' * X[i]' + D' * X[i]' * Xi[i];
    }
    for(i in 1:N){
      Theta_hat[((i-1)*P + 1):(i*P)] = Lambda * Eta[((i-1)*Q + 1):(i*Q)];
    }
  }
  model{
    for(i in 1:N){
      Y[((i-1)*tmax + 1):(i*tmax)] ~ normal(B * Theta[((i-1)*P + 1):(i*P)], sigma);
    }
    sigma ~ inv_gamma(.001,.001);
    for(i in 1:N){
      for(j in 1:P){
        Theta[(i-1)*P + j] ~ normal(Theta_hat[(i-1)*P + j], Phi_t[j]);
      }
    }
    Phi_t ~ inv_gamma(.001,.001);
    Xi ~ normal(0, 1);
}
"
library(rstan)
library(MASS)
library(dlnm)
library(plot3D)
options(mc.cores =parallel::detectCores())
rstan_options(auto_write = TRUE)
m <- stan_model(model_code = basis_version)
set.seed(1)
n <- 400
p <- 12
q <- 3
x_d <- 2
X <- cbind(rep(1,n), rnorm(n))
phi_t <- matrix(0, nrow = n, ncol = p)
for(i in 1:n){
  phi_t[i,] <- 1*mvrnorm(1, mu = rep(0, p), Sigma = diag(p))
}
Lambda <- matrix(rnorm(p * q), nrow = p, ncol = q)
D <- .5*matrix(rnorm(x_d * q), nrow = x_d, ncol = q)
Beta <- matrix(rnorm(x_d * q), nrow = x_d, ncol = q)
Xi <- rnorm(n)
Eta <- matrix(0, nrow = n, ncol = q)
for(i in 1:n){
  Eta[i,] <- t(Beta)%*%X[i,] + t(D) %*% X[i,] * Xi[i]
}
Theta <- matrix(0, nrow = n, ncol = p)
for(i in 1:n){
  Theta[i,] <- Lambda%*%Eta[i,] + phi_t[i,]
}
Theta_vec <- c(t(Theta))
noise_var <- .1
tmax <- 50
B <- ps(seq(from = 0, to = 1, length.out = tmax), df = p, intercept = TRUE)
Y <- numeric(tmax * n)
for(i in 1:n){
  Y[((i-1)*tmax + 1):(i*tmax)] <- B%*%Theta[i,] + noise_var * rnorm(tmax)
}
output_samples <- 1000
my_data <- list(N = n, P = 12, Q = 3, X = X, tmax = tmax, Y = Y, B = B)
set.seed(4)
fit_vb <- vb(m, data = my_data, par = c("Beta", "D", "Lambda"), output_samples = output_samples, iter = 20000)
posterior <- extract(fit_vb)
#fit_sampling <- sampling(m, data = my_data, iter = 50000)
#posterior_sampling <- extract(fit_sampling)
x <- c(1,.5)
true_mean <- B%*%Theta1%*%x
true_cov <- Btru%*%Lambda1%*%outer(x,x)%*%t(Lambda1)%*%t(Btru)
est_mean <- numeric(tmax)
est_cov <- matrix(0, nrow = tmax, ncol = tmax)
for(i in 1:output_samples){
  est_mean <- est_mean + B%*%posterior$Lambda[i,,]%*%t(posterior$Beta[i,,])%*%x
  est_cov <- est_cov + B%*%(posterior$Lambda[i,,]%*%t(posterior$D[i,,])%*%outer(x,x)%*%posterior$D[i,,]%*%t(posterior$Lambda[i,,]) + diag(p))%*%t(B)
}
est_mean <- est_mean / output_samples
est_cov <- est_cov / output_samples
#est_mean
#true_mean
#est_cov
#true_cov
plot(est_mean, type = "l")
lines(true_mean, col = "red")
for(i in 1:output_samples){
  lines(B%*%posterior$Lambda[i,,]%*%t(posterior$Beta[i,,])%*%x, col = "gray")
}
par(mfrow = c(1,2))
persp3D(1:tmax, 1:tmax, true_cov, theta = 90, phi = 20)
persp3D(1:tmax, 1:tmax, est_cov, theta = 90, phi = 20)
points(true_cov[lower.tri(true_cov)], col = "red")
for(i in 1:output_samples){
  temp <- posterior$Lambda[i,,]%*%t(posterior$D[i,,])%*%outer(x,x)%*%posterior$D[i,,]%*%t(posterior$Lambda[i,,]) + diag(p)
  points(temp[lower.tri(temp)])
}
