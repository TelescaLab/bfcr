basis_version <- "
  data{
    int<lower=0> N; // Number of samples
    int<lower=0> P; // Dimension of basis
    int<lower=0> Q; // Dimension of loading matrix
    vector[p * N] Theta; // Treating basis coefficients as response
    matrix[N,2] X;
  }
  parameters{
    matrix[10, Q] Lambda; // Loading matrix
    vector<lower=0>[p] Phi_t;
    matrix[2, Q] Beta; // Latent mean matrix
    matrix[2, Q] D; // Latent covariance matrix
    vector[N] Xi; //Latent cov random effect
  }
  transformed parameters{
    vector[p * N] Theta_hat;
    vector[N * Q] Eta;
    for(i in 1:N){
      Eta[((i-1)*Q + 1):(i*Q)] = Beta' * X[i]' + D' * X[i]' * Xi[i];
    }
    for(i in 1:N){
      Theta_hat[((i-1)*p + 1):(i*p)] = Lambda * Eta[((i-1)*Q + 1):(i*Q)];
    }
  }
  model{
    for(i in 1:N){
      for(j in 1:p){
        Theta[(i-1)*p + j] ~ normal(Theta_hat[(i-1)*p + j], Phi_t[j]);
      }
    }
    Phi_t ~ inv_gamma(.001,.001);
    Xi ~ normal(0, 1);
    //to_vector(Phi_t) ~ normal(0, 1);
  }
"
library(rstan)
options(mc.cores =parallel::detectCores())
rstan_options(auto_write = TRUE)
m <- stan_model(model_code = basis_version)
set.seed(2)
n <- 400
p <- 10
q <- 3
x_d <- 2
X <- cbind(rep(1,n), rnorm(n))
phi_t <- matrix(0, nrow = n, ncol = p)
for(i in 1:n){
  phi_t[i,] <- .5*mvrnorm(1, mu = rep(0, p), Sigma = diag(p))
}
Lambda <- matrix(rnorm(p * q), nrow = p, ncol = q)
D <- 5*matrix(rnorm(x_d * q), nrow = x_d, ncol = q)
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
output_samples <- 2000
my_data <- list(N = n, p = p, Theta = Theta_vec, X = X)
set.seed(3)
fit_vb <- vb(m, data = my_data, par = c("Beta", "D", "Lambda"), output_samples = output_samples)
posterior <- extract(fit_vb)
fit_sampling <- sampling(m, data = my_data, iter = 50000)
posterior_sampling <- extract(fit_sampling)
x <- c(1,-.5)
true_mean <- Lambda%*%t(Beta)%*%x
true_cov <- Lambda%*%t(D)%*%outer(x,x)%*%D%*%t(Lambda) + diag(p)
est_mean <- numeric(p)
est_cov <- matrix(0, nrow = p, ncol = p)
for(i in 1:output_samples){
  est_mean <- est_mean + posterior$Lambda[i,,]%*%t(posterior$Beta[i,,])%*%x
  est_cov <- est_cov + posterior$Lambda[i,,]%*%t(posterior$D[i,,])%*%outer(x,x)%*%posterior$D[i,,]%*%t(posterior$Lambda[i,,]) + diag(p)
}
est_mean <- est_mean / output_samples
est_cov <- est_cov / output_samples
#est_mean
#true_mean
#est_cov
#true_cov
plot(est_mean)
points(true_mean, col = "red")
plot(est_cov[lower.tri(est_cov)])
points(true_cov[lower.tri(true_cov)], col = "red")
