library(splines)
#####################################################################
################### MEAN AND VARIANCE FUNCTIONS #####################
#####################################################################
tmax <- 50
t <- seq(from = 0, to = 1, length.out = tmax)
p <- 12
B <- bs(t, df = p, intercept = TRUE)
int_function <- 3*(t-.5*t^2)
x_function <- t*sin(2*pi*t)
lbasis1_int <- .5*cos(t*pi)
lbasis2_int <- .5*sin(t*pi)
lbasis1_z <- .25*exp(cos(t))
lbasis2_z <- .25*exp(sin(t))

#####################################################################
####################### DATA GENERATION #############################
#####################################################################
X <- cbind(rep(1,n), rnorm(n, sd = 1))
noise_sd <- .05
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
Y_no_error <- t(mean_grid) + diag(Eta1)%*%t(basis1) + diag(Eta2)%*%t(basis2)
Y <- Y_no_error + E

# mu <- function(t, z){
#   t + z*sin(t*pi) + (1-z)*cos(t*pi)
# }
# 
# tmax <- 50
# zmax <- 100
# t <- seq(from = 0, to = 1, length.out = tmax)
# z <- seq(from = 0, to = 1, length.out = zmax)
# 
# meanfunc <- numeric(0)
# for(i in 1:zmax){
#   meanfunc <- c(meanfunc, mu(t, z[i]))
# }
# p <- 12
# B <- bs(t, df = p, intercept = TRUE)
# X <- numeric(0)
# X <- cbind(rep(1,zmax), z)
# X <- kronecker(X, B)
# reg1 <- lm(meanfunc ~ X - 1)
# num <- 29
# idx <- (1 + tmax*num):(tmax*(num+1))
# eig1 <- function(t,z){
#   -cos(pi*(t+z/2))/sqrt(2) * sqrt(z / 9)
# }
# eigfunc1 <- numeric(0)
# for(i in 1:zmax){
#   eigfunc1 <- c(eigfunc1, eig1(t,z[i]))
# }
# reg2 <- lm(eigfunc1 ~ X - 1)
# 
# eig2 <- function(t,z){
#   sin(pi*(t+z/2))/sqrt(2) * sqrt(z / 36)
# }
# eigfunc2 <- numeric(0)
# for(i in 1:zmax){
#   eigfunc2 <- c(eigfunc2, eig2(t, z[i]))
# }
# reg3 <- lm(eigfunc2 ~ X - 1)
# 
# Theta <- matrix(reg1$coefficients, nrow = p, ncol = 2)
# L1 <- matrix(reg2$coefficients, nrow = p, ncol = 2)
# L2 <- matrix(reg3$coefficients, nrow = p, ncol = 2)

#####################################################################
####################### DATA GENERATION #############################
#####################################################################
# n <- 40
# tmax <- 50
# p <- 12
# Q <- 2
# D <- 2
# T <- seq(from = 0, to = 1, length.out = tmax)
# Btru <- ps(T, df = p, intercept = TRUE)
# Best <- ps(T, df = ceiling(1.5 * p), intercept = TRUE)
# #X <- cbind(rep(1,n))
# #X <- cbind(rep(1,n), c(rep(0, n/2), rep(1,n/2)))
# #X <- cbind(rep(1,n),runif(n,min=0,max=1))
# X <- cbind(rep(1,n), rnorm(n, sd = 1))
# Intercept_only <- cbind(rep(1, n))
# #X <- Intercept_only
# d <- dim(X)[2]
# Eta1 <- rnorm(n, sd = sqrt(1))
# Eta2 <- rnorm(n, sd = sqrt(1))
# Y <- matrix(0, nrow = n, ncol = tmax)
# Lambda1 <-  L1*3
# Lambda2 <-  L2*3
# Theta1 <- Theta
# #Theta1 <- matrix(0, nrow = p, ncol = 2)
# #Theta1[,2] <- 0
# #Lambda1[,2] <- 0
# #Lambda2[,2] <- 0
# noise_sd <- .05
# # E <- matrix(rnorm(tmax * n,sd=noise_sd), nrow = n, ncol = tmax)
# E <- matrix(0, nrow = n, ncol = tmax)
# for(i in 1:n){
#   # E[i, ] <- rnorm(tmax, mean = 0, sd = noise_sd + (i %% 5) * noise_sd / 5)
#   E[i, ] <- rnorm(tmax, mean = 0, sd = noise_sd)
# }
# Y_no_error <- X%*%t(Theta1)%*%t(Btru) + diag(Eta1)%*%X%*%t(Lambda1)%*%t(Btru) + diag(Eta2)%*%X%*%t(Lambda2)%*%t(Btru)
# Y <- Y_no_error + E