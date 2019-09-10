library(BayesianConditionalFPCA)
library(splines)
library(MASS)
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
A <- numeric(10)

CovMetricsTotal <- array(dim = c(3, 500, 5)) 
MuMetricsTotal <- array(dim = c(3, 500, 5))
Eig1MetricsTotal <- array(dim = c(3, 500, 5))
Eig2MetricsTotal <- array(dim = c(3, 500, 5))
Eig1ValMetricsTotal <- array(dim = c(3, 500, 5))
Eig2ValMetricsTotal <- array(dim = c(3, 500, 5))
setwd("/Users/John/Documents/Johnstuff/BayesianConditionalFPCA/Rfuns")
ni <- 1
i <- 1
for(ni in 1:3){
  N <- c(100, 250, 500)
  for(i in 1:1){
    print(c(ni, i))
    n <- N[ni]
    Eta <- matrix(0, nrow = n, ncol = 2)
    Z <- runif(n)
    Eta[,1] <- sapply(1:n, function(i) rnorm(1, mean = 0, sd = 1))
    Eta[,2] <- sapply(1:n, function(i) rnorm(1, mean = 0, sd = 1))
    error_var <- (.05)^2
    Xd <- cbind(rep(1,n), Z)
    E <- matrix(rnorm(tmax * n, mean = 0, sd = sqrt(error_var)), nrow = n, ncol = tmax)
    Y <- matrix(0, nrow = n, ncol = tmax)
    #Y <- Xd%*%t(Theta)%*%t(B) + diag(Eta[,1])%*%Xd%*%t(L1)%*%t(B) + 
    #  diag(Eta[,2])%*%Xd%*%t(L2)%*%t(B) + E
    for(subj in 1:n){
      Y[subj,] <- mvrnorm(1,mu = B%*%Theta%*%Xd[subj,], Sigma = 
                B%*%L1%*%outer(Xd[subj,],Xd[subj,])%*%t(L1)%*%t(B) + 
                B%*%L2%*%outer(Xd[subj,],Xd[subj,])%*%t(L2)%*%t(B) + 
                diag(rep(.05^2),tmax))
    }
    res <- MCMC(Y, Xd, B, 1, 5000, 1, 1)
    zz <- seq(from = 0, to = 1, length.out = 5)
    m <- numeric(5)
    for(j in 1:5){
      zd <- c(1, zz[j])
      zd <- c(1,0)
      truthcov <- B%*%L1%*%outer(zd, zd)%*%t(L1)%*%t(B) +
        B%*%L2%*%outer(zd,zd)%*%t(L2)%*%t(B)
      trutheig1 <- eigen(truthcov)$vectors[,1]
      trutheig2 <- eigen(truthcov)$vectors[,2]
      cov <- matrix(0, nrow = tmax, ncol = tmax)
      mymean <- numeric(tmax)
      for(iter in 1001:5000){
        mymean <- B%*%res$Theta[[1]][,,iter]%*%zd
        for(a in 1:1){
          cov <- cov + B%*%res$Lambda[[1,iter]][,,a]%*%outer(zd,zd)%*%t(res$Lambda[[1,iter]][,,a])%*%t(B) * 1/res$Sigma[[1]][a,iter]
        }
      }
      cov <- cov / 4000
      CovMetricsTotal[ni, i, j] <- sum((cov-truthcov)^2)/sum(truthcov^2)
      MuMetricsTotal[ni, i, j] <- sum((mymean - B%*%Theta%*%zd)^2)/sum((B%*%Theta%*%zd)^2)
      Eig1MetricsTotal[ni, i, j] <- min(sum((eigen(cov)$vectors[,1] - trutheig1)^2),
                                        sum((eigen(cov)$vectors[,1] + trutheig1)^2))
      Eig2MetricsTotal[ni, i, j] <- min(sum((eigen(cov)$vectors[,2] - trutheig2)^2),
                                        sum((eigen(cov)$vectors[,2] + trutheig2)^2))
      Eig1ValMetricsTotal[ni, i, j] <- (eigen(cov)$values[1] - eigen(truthcov)$values[1])^2
      Eig2ValMetricsTotal[ni, i, j] <- (eigen(cov)$values[2] - eigen(truthcov)$values[2])^2
    }
  }
}
save(CovMetricsTotal, file = "CovMetricsTotal.RData")
save(MuMetricsTotal, file = "MuMetricsTotal.RData")
save(Eig1MetricsTotal, file = "Eig1MetricsTotal.RData")
save(Eig2MetricsTotal, file = "Eig2MetricsTotal.RData")
save(Eig1ValMetricsTotal, file = "Eig1ValMetricsTotal.RData")
save(Eig2ValMetricsTotal, file = "Eig2ValMetricsTotal.RData")
