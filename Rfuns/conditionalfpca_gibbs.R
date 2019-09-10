library(splines)
T <- seq(from = 0, to = 1, length.out = 50)
p <- 12
n <- 300
#X <- matrix(0, nrow = n, ncol = 2)
#X[,1] <- rep(1, n)
#X[1:floor(n/2),2] <- rep(1, floor(n/2))
X <- cbind(rep(1,n), runif(n))
B <- bs(T, df = p, intercept = TRUE)
Lambda1 <- matrix(rnorm(p * 2), nrow = p, ncol = 2)
Lambda2 <- matrix(rnorm(p * 2), nrow = p, ncol = 2)
#L1 <- matrix(rnorm(p * 2), nrow = p, ncol = 2)
#L2 <- matrix(rnorm(p * 2), nrow = p, ncol = 2)
Gamma1 <- diag(rnorm(n))
Gamma2 <- diag(rnorm(n))
mysd <- .05
E <- matrix(rnorm(length(T) * n, mean = 0, sd = mysd), nrow = n, ncol = length(T))
Thetatrue <- matrix(rnorm(p * 2), nrow = p, ncol = 2)
#Thetatrue <- matrix(0, nrow = p, ncol = 2)
Y <- X%*%t(Thetatrue) %*% t(B) + Gamma1 %*% X %*% t(L1)%*%t(B) + E
+ E + Gamma2 %*% X %*% t(L2) %*%t(B)

y <- list()
ot <- list()
BasisS <- list()
Z <- matrix(0, nrow = n, ncol = 2)
for(i in 1:n){
  missing <- sample((1:length(T)), sample(1:50, 1), replace = TRUE)
  #y[[i]] <- Y[i,-i]
  #y[[i]] <- Y[i,]
  y[[i]] <- Y[i, -missing]
  #ot[[i]] <- setdiff(T,T[i])
  ot[[i]] <- setdiff(T, T[missing])
  BasisS[[i]] <- bs(ot[[i]], df = p, intercept = TRUE)
  #BasisS[[i]] <- bs(T, df = p, intercept = TRUE)
  Z[i,] <- c(1,runif(1))
}

ptm <- proc.time()
mod<-MCMC_Impute(y, ot, T, X, B, k, 1000,1,1)
proc.time() - ptm
ptm <- proc.time()
mod2 <- MCMC(Y, X, B, 4, 5000, 1, 1)
proc.time() - ptm
ptm <- proc.time()
mod3 <- MCMC_Sparse(y, X, BasisS, k, 1000, 1, 1)
proc.time() - ptm

x <- c(1,1)
addL <- function(i){
  init <- matrix(0, nrow = tmax, ncol = tmax)
  for(k in 1:4){
    init <- init + B%*%mod2$Lambda[[1,i]][,,k]%*%x%*%t(B%*%mod2$Lambda[[1,i]][,,k]%*%x)
  }
  return(init)
}
LambdaList <- lapply(2500:5000, function(i) addL(i))
truth <- B%*%L1%*%x %*% t(B%*%L1%*%x) + B%*%L2%*%x%*%t(B%*%L2%*%x)
#v <- sapply(1:2500, function(i) LambdaList[[i]][5,5])
#hist(v)
#abline(v=truth[1,3])
persp3D(1:50,1:50, LambdaList[[2400]],theta=90,phi=10)
persp3D(1:50,1:50, truth, theta=90, phi=10)
subj <- 2
plot(ot[[subj]], y[[subj]],ylim=c(-5,5))
myfunc <- function(i){
  mysum <- 0
  for(j in 1:3){
    #mysum <- mysum + B%*%mod2$Lambda[[1,i]][,,j]%*%X[30,]*mod2$Eta[[1]][30,j,i]
    mysum <- mysum + mod3$Eta[[1]][subj,j,i] * BB %*% mod3$Lambda[[1,i]][,,j]%*%X[subj,]
  }
  return(mysum)
}
sapply(800:1000, function(i) lines(seq(from = 0, to = 1, length.out=1000), BB%*%mod3$Theta[[1]][,,i]%*%X[subj,] + myfunc(i), col = "blue"))


subj <- 2
plot(seq(from = 0, to = 1, length.out = 50),Y[subj,],ylim=c(-5,5))
myfunc <- function(i){
  mysum <- 0
  for(j in 1:k){
    mysum <- mysum + BB%*%mod$Lambda[[1,i]][,,j]%*%X[subj,]*mod$Eta[[1]][subj,j,i]
  }
  return(mysum)
}
sapply(500:1000, function(i) lines(seq(from = 0, to = 1, length.out = 1000), BB%*%mod$Theta[[1]][,,i]%*%X[subj,] + myfunc(i), col = "blue"))



for(i in 1:10000){
  m <- sample(1:25, replace = TRUE)
  lines(seq(from = 0, to = 1, length.out = 200), colMeans(Y[m,]),col="blue")
}
Ytilde <- X %*% t(Theta) %*%t(B) + diag(G[,1])%*%X %*%t(L[,,1])%*%t(B)# + diag(G[,2]) %*%X%*%t(L[,,2])%*%t(B) #+ diag(G[,3]) %*%X%*%t(L[,,3])%*%t(B) + diag(G[,4]) %*%X%*%t(L[,,4])%*%t(B) + diag(G[,5]) %*%X%*%t(L[,,5])%*%t(B)
lines(seq(from = 0, to = 1, length.out = 200),Ytilde[2,],col="green")
lines(seq(from = 0, to = 1, length.out = 200),B%*%Theta[,1],type="l")
lines(seq(from = 0, to = 1, length.out = 200), B%*%Thetatrue[,1])
plot(T,Y[25,])
lines(seq(from = 0, to = 1, length.out = 200), colMeans(Y[1:25,]),type="l")
lines(T, B%*%Thetatrue%*%c(1,1))

# Group 1
empth1 <- colMeans(Y[1:25,])
plot(empth1,type="l",ylim=c(-5,5))
lines(B%*%Thetatrue[,1] + B%*%Thetatrue[,2])
# Bootstrapped
sapply(1:1000, function(i){
  m <- sample(1:25, replace = TRUE)
  lines(colMeans(Y[m,]))})
sapply(2500:5000, function(i) lines(B%*%mod$Theta[[1]][,1,i] + B%*%mod$Theta[[1]][,2,i],col="blue"))
sapply(2500:5000, function(i) lines(B%*%mod2$Theta[[1]][,1,i] + B%*%mod2$Theta[[1]][,2,i],col="blue"))

# Group 2
empth2 <- colMeans(Y[26:50,])
plot(empth2, type = "l")
# Bootstrapped
sapply(1:1000, function(i){
  m <- sample(26:50, replace = TRUE)
  lines(colMeans(Y[m,]))})
sapply(500:1000, function(i) lines(B%*%mod$Theta[[1]][,1,i] ,col="blue"))
sapply(500:1000, function(i) lines(B%*%mod2$Theta[[1]][,1,i],col="red"))


subject <- 2
plot(Y[subject,],ylim=c(-5,5))
myfunc <- function(i){
  mysum <- 0
  for(j in 1:k){
    mysum <- mysum + B%*%mod2$Lambda[[1,i]][,,j]%*%X[subject,]*mod2$Eta[[1]][subject,j,i]
  }
  return(mysum)
}
sapply(3000:5000, function(i) lines(B%*%mod2$Theta[[1]][,,i]%*%X[subject,] + myfunc(i), col = "blue"))


ptm <- proc.time()
for(i in 1:1000){
  updateLambda2(Y, L, R, G, X, B, 20, Theta)
  
}
proc.time() - ptm

ptm <- proc.time()
for(i in 1:1000){
  updateTheta(Y, L, c(1,1), G, X, B, 20, Theta)
}
proc.time() - ptm
ptm <- proc.time()
for(i in 1:1000){
  updateEta(Y, L, rep(1,5), G, X, B, 20, Theta)
}
proc.time() - ptm



ptm <- proc.time()
for(i in 1:1000){
  updateLambda(Y, L, rep(1,5), G, X, B, 1/.5^2)
}
proc.time() - ptm

ptm <- proc.time()
for(i in 1:1000){
  prec <- updatePrec(Y, L, G, X, B, Theta)
}
proc.time() - ptm