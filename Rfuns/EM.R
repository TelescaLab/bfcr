setwd("/Users/John/Documents/Johnstuff/BayesianConditionalFPCA/Rfuns")
{
  library(MASS)
  library(splines)
  library(mvtnorm)
  library(plot3D)
  library(splines)
  library(Rcpp)
  library(microbenchmark)
  library(dlnm)
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
  set.seed(2)
  n <- 200
  tmax <- 50
  p <- 12
  T <- seq(from = 0, to = 1, length.out = tmax)
  library(plot3D)
  library(splines)
  #Btru <- ps(T, df = p)
  Btru <- bs(T, df = p, intercept = TRUE)
  #X <- cbind(rep(1,n))
  #X <- cbind(rep(1,n),runif(n,min=-1,max=1))
  X <- cbind(rep(1,n), rnorm(n, sd = 1))
  d <- dim(X)[2]
  Eta1 <- rnorm(n)
  Eta2 <- rnorm(n)
  Y <- matrix(0, nrow = n, ncol = tmax)
  #Lambda1 <- matrix(rnorm(p * dim(X)[2]), nrow = p, ncol = d)
  #Lambda2 <- matrix(rnorm(p * dim(X)[2]), nrow = p, ncol = d)
  Lambda1 <- L1
  Lambda2 <- L2
  #Lambda%*%t(Lambda)
  #Theta <- matrix(rnorm(p * dim(X)[2]), nrow = p, ncol = dim(X)[2])
  E <- matrix(rnorm(tmax * n,sd=.001), nrow = n, ncol = tmax)
  Y <- X%*%t(Theta)%*%t(Btru) + diag(Eta1)%*%X%*%t(Lambda1)%*%t(Btru) + E + diag(Eta2)%*%X%*%t(Lambda2)%*%t(Btru)# + E
  inflation <- 5
  Et1 <- matrix(rnorm(tmax * n, sd = inflation), nrow = n, ncol = tmax)
  Yt <- Y + Et1
}
plot(Y[1,],type="p")
{
  set.seed(3)
  p <- 25
  #B <- bs(T, df = p, intercept = TRUE)
  B <- ps(T, df = p, diff = 1, intercept = TRUE)
  K <- 2
  Xmat <- kronecker(B, X)
  reg <- lm(c(Y) ~ Xmat - 1)
  Theta_init <- t(matrix(reg$coefficients, nrow = 2))
  param <- cpp_EM(X, B, Y, K,Theta_init, 12)
  cpploglik(Theta, cbind(L1,L2), 1/.001^2, X, Btru, Y, 2, 6)
  cpploglik(param$Theta, param$Lambda, param$Precision, X, B, Y, K,1)

  max_iter <- 15000
  burnin <- 5000
  thin <- 1
  nchain <- 6
  Lambda_init <- array(param$Lambda, dim = c(p, 2, 2))
  Eta_init <- t(param$EtaM)
  bayes_param <- MCMC(Y, X, B, K, max_iter, nchain, thin, param$Theta, Lambda_init, Eta_init)
  bayes_logliks <- sapply(seq(from = 1, to = max_iter, by = 1), function(i) cpploglik(bayes_param$Theta[[1]][,,i], array(bayes_param$Lambda[[1,i]], dim = c(p,2*K)), bayes_param$Prec[[1]][i], X, B, Y, K, 6))
}
plot(bayes_logliks,type="l")
dev.off()
x <- c(1,-1)
bayes_mean <- matrix(0, nrow = p, ncol = 2)
for(i in 1:nchain){
  bayes_mean <- bayes_mean + apply(bayes_param$Theta[[i]][,,burnin:max_iter], c(1,2), mean)
}
bayes_mean <- bayes_mean / nchain

plot(Btru%*%Theta%*%x, type = "l")
lines(B%*%param$Theta%*%x, col = "blue")
lines(B%*%bayes_mean%*%x, col = "red")
for(chain in 1:nchain){
  for(i in seq(from = burnin, to = max_iter, by = 10)){
    lines(B%*%bayes_param$Theta[[chain]][,,i]%*%x, col = "gray")
  }
}
lines(Btru%*%Theta%*%x, type = "l",col="red")

sum((B%*%param$Theta%*%x - Btru%*%Theta%*%x)^2)/sum((Btru%*%Theta%*%x)^2)*100
sum((B%*%bayes_mean%*%x - Btru%*%Theta%*%x)^2)/sum((Btru%*%Theta%*%x)^2)*100

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
persp3D(1:tmax,1:tmax, covtruth, theta=90,phi=10)
persp3D(1:tmax,1:tmax, covfreq, theta=90,phi=10, zlim = c(min(unlist(covtruth)),max(unlist(covtruth))))
persp3D(1:tmax,1:tmax, covbayes, theta=90,phi=10, zlim = c(min(unlist(covtruth)),max(unlist(covtruth))))


dev.off()
subj <- 101
iter <- 3000
plot(Y[subj,],type="p")
lines(B%*%param$Theta%*%X[subj,] + B%*%param$Lambda[,1:2]%*%X[subj,] * param$EtaM[1,subj] +
       B%*%param$Lambda[,3:4]%*%X[subj,] * param$EtaM[2,subj],col="blue")
lines(B%*%bayes_param$Theta[[1]][,,iter]%*%X[subj,] + bayes_param$Eta[[1]][subj,1,iter] * B%*%bayes_param$Lambda[[1,iter]][,,1]%*%X[subj,]+
        bayes_param$Eta[[1]][subj,2,iter] * B%*%bayes_param$Lambda[[1,iter]][,,2]%*%X[subj,], col = "red")
sapply(burnin:max_iter, function(iter) lines(B%*%bayes_param$Theta[[1]][,,iter]%*%X[subj,] + bayes_param$Eta[[1]][subj,1,iter] * B%*%bayes_param$Lambda[[1,iter]][,,1]%*%X[subj,]+
                                       bayes_param$Eta[[1]][subj,2,iter] * B%*%bayes_param$Lambda[[1,iter]][,,2]%*%X[subj,], col="gray"))
