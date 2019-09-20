setwd("/Users/John/Documents/Johnstuff/BayesianConditionalFPCA/Rfuns")
#setwd("E:/Rcpp stuff/BFPCA")
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
  set.seed(6)
  n <- 500
  tmax <- 50
  p <- 12
  T <- seq(from = 0, to = 1, length.out = tmax)
  library(plot3D)
  library(splines)
  #Btru <- ps(T, df = p)
  Btru <- ps(T, df = p, intercept = TRUE)
  #X <- cbind(rep(1,n))
  #X <- cbind(rep(1,n), c(rep(0, n/2), rep(1,n/2)))
  #X <- cbind(rep(1,n),runif(n,min=-1,max=1))
  X <- cbind(rep(1,n), rnorm(n, sd = 1))
  d <- dim(X)[2]
  Eta1 <- rnorm(n)
  Eta2 <- rnorm(n)
  Y <- matrix(0, nrow = n, ncol = tmax)
  #Lambda1 <- matrix(rnorm(p * dim(X)[2]), nrow = p, ncol = d)
  #Lambda2 <- matrix(rnorm(p * dim(X)[2]), nrow = p, ncol = d)
  #Lambda1 <- L1[,1]
  #Lambda2 <- L2[,1]
  Lambda1 <-  1*L1
  Lambda2 <-  1*L2
  #Theta1 <- Theta[,1]
  Theta1 <- 1*Theta
  #X <- as.matrix(X[,1])
  #Lambda%*%t(Lambda)
  #Theta <- matrix(rnorm(p * dim(X)[2]), nrow = p, ncol = dim(X)[2])
  noise_sd <- .1
  E <- matrix(rnorm(tmax * n,sd=noise_sd), nrow = n, ncol = tmax)
  Y <- X%*%t(Theta1)%*%t(Btru) + diag(Eta1)%*%X%*%t(Lambda1)%*%t(Btru) + E + diag(Eta2)%*%X%*%t(Lambda2)%*%t(Btru)# + E
  inflation <- 5
  Et1 <- matrix(rnorm(tmax * n, sd = inflation), nrow = n, ncol = tmax)
  Yt <- Y + Et1
}
dev.off()
plot(Y[1,],type="p")
n_500_high_noise_high_between <- numeric(100)
for(i in 1:100){
  {
    set.seed(6)
    p <- 12
    #B <- bs(T, df = p, intercept = TRUE)
    B <- ps(T, df = p, intercept = TRUE)
    K <- 2
    Xmat <- kronecker(B, X)
    reg <- lm(c(Y) ~ Xmat - 1)
    Theta_init <- t(matrix(reg$coefficients, nrow = 2))
    Theta_init <- Theta1
    #Lambda_init <- cbind(Lambda1, Lambda2)
    Lambda_init <- matrix(rnorm(12 * 2 * K), ncol = 2*K)
    Eta_init <- rbind(Eta1, Eta2)
    param <- cpp_EM(X, Btru, Y, K,Theta_init, 12)
    param_new <- cpp_EM_new(X, Btru, Y, K, Theta_init, Lambda_init, 12)
    print(c("log likelihood to beat is", cpploglik(Theta1, cbind(Lambda1,Lambda2), 1/noise_sd^2, X, Btru, Y, 2, 6)))
#    n_500_high_noise_high_between[i] <- cpploglik(param$Theta, param$Lambda, param$Precision, X, B, Y, K,1)
    
    max_iter <- 5000
    burnin <- 1000
    thin <- 1
    nchain <- 1
    Lambda_init <- array(param$Lambda, dim = c(p, 2, K))
    #Lambda_init <- array(cbind(Lambda1, Lambda2), dim = c(p, 2, K))
    Eta_init <- t(param$EtaM)
    #Eta_init <- cbind(Eta1, Eta2)
    Prec_init <- param$Precision
    set.seed(1)
    bayes_param <- MCMC(Y, X, B, K, max_iter, nchain, thin, Theta_init, Lambda_init, Eta_init, Prec_init)
    #bayes_logliks <- sapply(seq(from = 1, to = max_iter, by = 1), function(i) cpploglik(bayes_param$Theta[[1]][,,i], array(bayes_param$Lambda[[1,i]], dim = c(p,2*K)), bayes_param$Prec[[1]][i], X, B, Y, K, 6))
  }
}

dev.off()
x <- c(1,-1.5)
bayes_mean <- matrix(0, nrow = p, ncol = 2)
for(i in 1:nchain){
  bayes_mean <- bayes_mean + apply(bayes_param$Theta[[i]][,,burnin:max_iter], c(1,2), mean)
}
bayes_mean <- bayes_mean / nchain

plot(T,Btru%*%Theta1%*%x, type = "l", ylab = "Mean, x = -1", xlab = "t")
lines(T,B%*%Theta_init%*%x,col="green")
lines(T,Btru%*%param$Theta%*%x, col = "blue")
lines(T,B%*%bayes_mean%*%x, col = "red")
for(chain in 1:nchain){
  for(i in seq(from = burnin, to = max_iter, by = 10)){
    lines(T,B%*%bayes_param$Theta[[chain]][,,i]%*%x, col = "gray")
  }
}
lines(T,Btru%*%Theta1%*%x,col="red")

sum((B%*%param$Theta%*%x - Btru%*%Theta1%*%x)^2)/sum((Btru%*%Theta1%*%x)^2)*100
sum((B%*%bayes_mean%*%x - Btru%*%Theta1%*%x)^2)/sum((Btru%*%Theta1%*%x)^2)*100

covfreq <- matrix(0, p, p)
for(k in 1:K){
  covfreq <- covfreq + param$Lambda[,(d*(k-1)+1):(d*k)]%*%outer(x,x)%*%t(param$Lambda[,(d*(k-1)+1):(d*k)])
}
covfreq <- B%*%covfreq%*%t(B)

covbayes <- matrix(0, nrow = tmax, ncol = tmax)
for(chain in 1:nchain){
  for(iter in burnin:max_iter){
    for(a in 1:K){
      covbayes <- covbayes + B%*%bayes_param$Lambda[[chain,iter]][,,a]%*%outer(x,x)%*%t(bayes_param$Lambda[[chain,iter]][,,a])%*%t(B) 
    }
  }
}
covbayes <- covbayes / ((max_iter - burnin + 1) * nchain)

covtruth <- Btru%*%Lambda1%*%outer(x,x)%*%t(Lambda1)%*%t(Btru) + Btru%*%Lambda2%*%outer(x,x)%*%t(Lambda2)%*%t(Btru)

sum((covfreq - covtruth)^2)/sum((covtruth)^2)*100
sum((covbayes - covtruth)^2)/sum((covtruth^2))*100
par(mfrow = c(1,3))
persp3D(1:tmax,1:tmax, covtruth, theta=90,phi=10, main = "Truth", colkey = FALSE)
persp3D(1:tmax,1:tmax, covfreq, theta=90,phi=10, zlim = c(min(unlist(covtruth)),max(unlist(covtruth))), main = "EM algorithm", colkey = FALSE)
persp3D(1:tmax,1:tmax, covbayes, theta=90,phi=10, zlim = c(min(unlist(covtruth)),max(unlist(covtruth))), main = "Gibbs sampling", colkey = FALSE)


dev.off()
subj <- 3
iter <- 3000
plot(Y[subj,],type="p")
lines(B%*%param$Theta%*%X[subj,] + B%*%param$Lambda[,1:2]%*%X[subj,] * param$EtaM[1,subj] +
       B%*%param$Lambda[,3:4]%*%X[subj,] * param$EtaM[2,subj],col="blue")
lines(B%*%bayes_param$Theta[[1]][,,iter]%*%X[subj,] + bayes_param$Eta[[1]][subj,1,iter] * B%*%bayes_param$Lambda[[1,iter]][,,1]%*%X[subj,]+
        bayes_param$Eta[[1]][subj,2,iter] * B%*%bayes_param$Lambda[[1,iter]][,,2]%*%X[subj,], col = "red")
sapply(burnin:max_iter, function(iter) lines(B%*%bayes_param$Theta[[1]][,,iter]%*%X[subj,] + bayes_param$Eta[[1]][subj,1,iter] * B%*%bayes_param$Lambda[[1,iter]][,,1]%*%X[subj,]+
                                       bayes_param$Eta[[1]][subj,2,iter] * B%*%bayes_param$Lambda[[1,iter]][,,2]%*%X[subj,], col="gray"))

