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
options(mc.cores =parallel::detectCores())
rstan_options(auto_write = TRUE)
m <- stan_model(model_code = basis_version)
{
  library(MASS)
  library(splines)
  library(mvtnorm)
  library(plot3D)
  library(splines)
  library(Rcpp)
  library(microbenchmark)
  library(dlnm)
  library(BayesianConditionalFPCA)
  mu <- function(t, z){
    t + z*sin(t) + (1-z)*cos(t)
  }
  
  tmax <- 50
  zmax <- 30
  t <- seq(from = 0, to = 1, length.out = tmax)
  z <- seq(from = 0, to = 1, length.out = zmax)
  
  meanfunc <- numeric(0)
  for(i in 1:zmax){
    meanfunc <- c(meanfunc, mu(t, z[i]))
  }
  p <- 12
  B <- bs(t, df = p, intercept = TRUE)
  X <- numeric(0)
  X <- cbind(rep(1,zmax), z)
  X <- kronecker(X, B)
  reg1 <- lm(meanfunc ~ X - 1)
  num <- 29
  idx <- (1 + tmax*num):(tmax*(num+1))
  plot(meanfunc[idx],type="l")
  lines(X[idx,]%*%reg1$coefficients,col="blue")
  
  eig1 <- function(t,z){
    -cos(pi*(t+z/2))/sqrt(2) * sqrt(z/9)
  }
  eigfunc1 <- numeric(0)
  for(i in 1:zmax){
    eigfunc1 <- c(eigfunc1, eig1(t,z[i]))
  }
  reg2 <- lm(eigfunc1 ~ X - 1)
  
  eig2 <- function(t,z){
    sin(pi*(t+z/2))/sqrt(2) * sqrt(z/36)
  }
  eigfunc2 <- numeric(0)
  for(i in 1:zmax){
    eigfunc2 <- c(eigfunc2, eig2(t, z[i]))
  }
  reg3 <- lm(eigfunc2 ~ X - 1)
  plot(eigfunc1[idx], type="l")
  lines(X[idx,]%*%reg2$coefficients, col = "blue")
  
  Theta <- matrix(reg1$coefficients, nrow = p, ncol = 2)
  L1 <- matrix(reg2$coefficients, nrow = p, ncol = 2)
  L2 <- matrix(reg3$coefficients, nrow = p, ncol = 2)
}

{
  set.seed(5)
  n <- 500
  tmax <- 50
  p <- 12
  T <- seq(from = 0, to = 1, length.out = tmax)
  library(plot3D)
  library(splines)
  #Btru <- ps(T, df = p)
  Btru <- bs(T, df = p, intercept = TRUE)
  #X <- cbind(rep(1,n))
  X <- cbind(rep(1,n), c(rep(0, n/2), rep(1,n/2)))
  #X <- cbind(rep(1,n),runif(n,min=-1,max=1))
  #X <- cbind(rep(1,n), rnorm(n, sd = 1))
  d <- dim(X)[2]
  Eta1 <- rnorm(n)
  Eta2 <- rnorm(n)
  Y <- matrix(0, nrow = n, ncol = tmax)
  #Lambda1 <- matrix(rnorm(p * dim(X)[2]), nrow = p, ncol = d)
  #Lambda2 <- matrix(rnorm(p * dim(X)[2]), nrow = p, ncol = d)
  #Lambda1 <- L1[,1]
  #Lambda2 <- L2[,1]
  Lambda1 <-  .1*L1
  Lambda2 <-  1*L2
  #Theta1 <- Theta[,1]
  Theta1 <- 1*Theta
  #X <- as.matrix(X[,1])
  #Lambda%*%t(Lambda)
  #Theta <- matrix(rnorm(p * dim(X)[2]), nrow = p, ncol = dim(X)[2])
  noise_sd <- .01
  E <- matrix(rnorm(tmax * n,sd=noise_sd), nrow = n, ncol = tmax)
  Y <- X%*%t(Theta1)%*%t(Btru) + diag(Eta1)%*%X%*%t(Lambda1)%*%t(Btru) + E #+ diag(Eta2)%*%X%*%t(Lambda2)%*%t(Btru)# + E
  inflation <- 5
  Et1 <- matrix(rnorm(tmax * n, sd = inflation), nrow = n, ncol = tmax)
  Yt <- Y + Et1
}
output_samples <- 1000
my_data <- list(N = n, P = p, Q = 10, X = X, tmax = tmax, Y = c(t(Y)), B = Btru)
set.seed(4)
fit_vb <- vb(m, data = my_data, par = c("Beta", "D", "Lambda", "Phi_t", "sigma"), output_samples = output_samples, iter = 20000)
posterior <- extract(fit_vb)
#fit_sampling <- sampling(m, data = my_data, iter = 50000)
#posterior_sampling <- extract(fit_sampling)
x <- c(1,1)
true_mean <- Btru%*%Theta1%*%x
true_cov <- Btru%*%Lambda1%*%outer(x,x)%*%t(Lambda1)%*%t(Btru)
est_mean <- numeric(tmax)
est_cov <- matrix(0, nrow = tmax, ncol = tmax)
for(i in 1:output_samples){
  est_mean <- est_mean + Btru%*%posterior$Lambda[i,,]%*%t(posterior$Beta[i,,])%*%x
  est_cov <- est_cov + Btru%*%(posterior$Lambda[i,,]%*%t(posterior$D[i,,])%*%outer(x,x)%*%posterior$D[i,,]%*%t(posterior$Lambda[i,,]) + diag(posterior$Phi_t[i,]))%*%t(Btru)
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
