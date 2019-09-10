library(splines)
library(MASS)
library(fdapace)
library(doParallel)
library(doSNOW)
library(coda)
library(MCMCglmm)
library(LFBayes)
#setwd("/Users/John/Downloads/LongFunc Code/ChenCode")
setwd("/Users/John/Documents/Johnstuff/LFBayes/Rfuns")

source("MarginalFPCA.R")
source("ProductFPCA.R")

errorvar <- .025
SS <- 20
TT <- 20
t <- seq(from = 0, to = 1, length.out = TT)
s <- seq(from = 0, to = 1, length.out = SS)
n <- 150
tt <- list()
tt[[1]] <- 1:(TT*SS)
tt <- rep(tt, n)
p1 <- 10
p2 <- 10
q1 <- 2
q2 <- 2
Bt <- bs(t, df = p1, intercept = TRUE)
Bs <- bs(s, df = p2, intercept = TRUE)
Bt1 <- bs(t, df = 12, intercept = TRUE)
Bs1 <- bs(s, df = 12, intercept = TRUE)
productcov <- function(eig, scores){
  n <- dim(scores)[1]
  neig <- dim(eig)[2]
  v <- numeric(dim(eig)[1])
  mycov <- matrix(0, nrow = length(v), ncol = length(v))
  for(l in 1:neig){
    v = eig[,l]
    mycov <- mycov + outer(v, v) * 1/n * as.numeric(scores[,l]%*%scores[,l])
  }
  mycov
}
marginalcov <- function(eig, scores){
  n <- dim(scores)[1]
  neig <- dim(eig)[2]
  v <- numeric(dim(eig)[1])
  mycov <- matrix(0, nrow = length(v), ncol = length(v))
  for(j in 1:neig){
    v <- eig[,j]
    mycov <- mycov + outer(v,v) * 1/n * as.numeric(scores[,j]%*%scores[,j])
  }


  mycov
}
Loading.Brown.Bridge <- function(t, p, k){
  B <- bs(t, df = p, intercept = TRUE)
  Loading <- matrix(nrow = p, ncol = k)
  for(i in 1:k){
    eigval <- 1/(i^2*pi^2)
    #if(i %% 2 == 0){
    #  psi <- sqrt(2) * -sqrt(2)*cos((i+1)*pi*t)
    #}else{
      psi <- sqrt(2)*sin(i*pi*t)
    #}

    Loading[,i] <- sqrt(eigval)*solve(t(B)%*%B)%*%t(B)%*%psi
  }
  return(Loading)
}
Loading.Brown.Motion <- function(t, p, k){
  B <- bs(t, df = p, intercept = TRUE)
  Loading <- matrix(nrow = p, ncol = k)
  for(i in 1:k){
    eigval <- 1/((k - 1/2)^2*pi^2)
    psi <- sqrt(2)*sin((k-1/2)*pi*t)
    Loading[,i] <- sqrt(eigval)*solve(t(B)%*%B)%*%t(B)%*%psi
  }
  return(Loading)
}
Matern.Cov <- function(s){
  rho = 0.5
  Matern.Cov <- matrix(nrow = length(s), ncol = length(s))
  for(i in 1:length(s)){
    for(j in 1:length(s)){
      d <- abs(s[i] - s[j])
      Matern.Cov[i, j] <- (1 + sqrt(5) * d / rho + 5 * d^2 / (3*rho^2)) * exp(-sqrt(5) * d / rho)
    }
  }
  Matern.Cov
}

Loading.Matern <- function(s, p, k, B){
  rho = 0.5
  Matern.Cov <- matrix(nrow = length(s), ncol = length(s))
  for(i in 1:length(s)){
    for(j in 1:length(s)){
      d <- abs(s[i] - s[j])
      Matern.Cov[i, j] <- (1 + sqrt(5) * d / rho + 5 * d^2 / (3*rho^2)) * exp(-sqrt(5) * d / rho)
    }
  }
  Loading <- matrix(nrow = p, ncol = k)
  evec <- eigen(Matern.Cov)$vectors
  eval <- eigen(Matern.Cov)$values
  for(i in 1:k){
    Loading[,i] <- sqrt(eval[i]) * solve(t(B)%*%B)%*%t(B)%*%evec[,i]
  }
  Loading
}


Loading.CosCov <- function(s, p, q, B){
  alpha <- 2
  CosCov <- matrix(0,nrow=length(s),ncol=length(s))
  for(k in 1:50){
    CosCov <- k^(-2*alpha)*outer(cos(k*pi*s),cos(k*pi*s)) + CosCov
  }
  Loading <- matrix(nrow = p, ncol = q)
  evec <- eigen(CosCov)$vectors
  eval <- eigen(CosCov)$values
  for(i in 1:q){
    print(i)
    Loading[,i] <- sqrt(eval[i]) * solve(t(B)%*%B)%*%t(B)%*%evec[,i]
  }
  Loading
}
GenerateNonsep <- function(s,t){
  Cov <- matrix(0, nrow = length(t)*length(s), ncol = length(t)*length(s))
  S <- length(s)
  T <- length(t)
  for(i in 1:S){
    for(ii in 1:S){
      for(j in 1:T){
        for(jj in 1:T){
          Cov[(i-1)*T + j, (ii - 1)*T + jj] <- 1/((t[j]-t[jj])^2+1) * exp(-((s[i]-s[ii])^2)/((t[j]-t[jj])^2+1))
        }
      }
    }
  }
  Cov
}
GenerateH <- function(q1,q2){
  H <- matrix(0,q1, q2)
  for(i in 1:q1){
    for(j in 1:q2){
      H[i,j] <- exp(-(sqrt(.01*i) + sqrt(.1*j)))
    }
  }
  H <- diag(c(H))
  H
}
GenerateMu1 <- function(s,t){
  mu <- matrix(0, nrow= length(t),ncol=length(s))
  for(i in 1:length(t)){
    for(j in 1:length(s)){
      mu[i,j] <- sqrt(1/(5*sqrt(s[j]+1)))*sin(5*t[i])
    }
  }
  c(mu)
}
GenerateMu2 <- function(s,t){
  mu <- matrix(0, nrow = length(t),ncol=length(s))
  for(i in 1:length(t)){
    for(j in 1:length(s)){
      mu[i,j] <- 5*sqrt(1-(s[j]-.5)^2-(t[i]-.5)^2)
    }
  }
  c(mu)
}
GenerateMu3 <- function(s,t){
  mu <- matrix(0, nrow = length(t),ncol=length(s))
  for(i in 1:length(t)){
    for(j in 1:length(s)){
      mu[i,j] <- sqrt(1+sin(pi*s[j]) + cos(pi*t[i]))
    }
  }
  c(mu)
}
H <- GenerateH(q1, q2)

mu1 <- GenerateMu1(s,t)
Lambda <- Loading.Matern(t, p1, q1, Bt)
Gamma <- Loading.Brown.Bridge(s, p2, q2)
Cov <- kronecker(Bs%*%Gamma, Bt%*%Lambda)%*%H%*%t(kronecker(Bs%*%Gamma, Bt%*%Lambda)) + errorvar * diag(SS * TT)

#image(Cov[,200:1],col=heat.colors(100))
#mu2 <- GenerateMu2(s,t)
#Lambda <- Loading.Brown.Bridge(t, p1, q1)
#Gamma <- Loading.CosCov(s,p2,q2,Bs)
#Cov <- kronecker(Bs%*%Gamma, Bt%*%Lambda)%*%H%*%t(kronecker(Bs%*%Gamma, Bt%*%Lambda)) + errorvar * diag(SS * TT)

#mu3 <- GenerateMu3(s,t)
#Cov <- GenerateNonsep(s,t)
pc.j = NULL
pc.k = NULL
fpca.op1 = list(dataType = "Sparse", maxK = pc.j, FVEthreshold = .9999, nRegGrid = TT)
fpca.op2 = list(dataType = "Sparse", maxK = pc.k, FVEthreshold = .9999, nRegGrid = SS)

iter <- 5000 # Number of iterations
burnin <- 1000 # Burnin iterations
thin <- 1 # Thinning for each chain
nchain <- 1 # Number of chains
neig <- 3 # Number of eigenfunctions for inference
info_thin <- 10 # Thinning for information criteria (shortens computation time)

splinenum <- 10 # Number of splines used in estimation
q1s <- 4 # Number of latent factors for functional dimension
q2s <- 4 # Number of latent factors for longitudinal dimension
Bt1 <- bs(t, df = splinenum, intercept = TRUE)
Bs1 <- bs(s, df = splinenum, intercept = TRUE)
set.seed(1)
setwd("/Users/John/Documents/Johnstuff/splines")
iterations <- 1
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cl <- makeCluster(5)
registerDoSNOW(cl)
system.time(simulation <-foreach(index=1:iterations,.combine=rbind, .packages = c("MASS", "LFBayes", "fdapace"), .options.snow = opts)%dopar%{
  print(index)
  #x <- mvrnorm(n, mu  = rep(0, TT*SS), Sigma = Cov.Weak)
  x <- mvrnorm(n, mu  = as.vector(mu1), Sigma = Cov)
  sx <- sd(x)
  mx <- mean(x)
  x <- (x-mx)/sx
  Smooth_scaled_cov <- (Cov - errorvar * diag(SS * TT)) / sx^2
  mu <- (mu1 - mx)/sx
  Marg.Long <- getMarginalLong(Smooth_scaled_cov,SS,TT)
  Marg.Func <- getMarginalFunc(Smooth_scaled_cov,SS,TT)
  #scaled_cov <- Cov/sx^2
  y <- lapply(1:n, function(i) x[i,])
  missing <- list()
  for(i in 1:n){
    missing[[i]] <- numeric(0)
  }
  X <- rep(1,n)
  dim(X) <- c(n,1)

  MCMC <- mcmcWeakChains(y, missing, X, Bs1, Bt1, q1s, q2s, iter, thin, burnin, nchain)
  MCMC_eigen <- eigenLFChains(Bs1, Bt1, MCMC, neig, iter, burnin, nchain)
  infoD <- calculate_DIC(t(x), X, MCMC, Bs1, Bt1, burnin, info_thin)
  infoB <- calculate_BIC(t(x), X, MCMC, Bs1, Bt1, burnin, info_thin)
  
  resMarginal <- MarginalFPCA(x, n, SS, TT, fpca.op1, fpca.op2, pc.j, rep(pc.k, pc.j))
  resMarginalCov <- marginalcov(resMarginal$eig, resMarginal$scores)
  resProduct <- ProductFPCA(x, n, SS, TT, fpca.op1, fpca.op2, pc.j, rep(pc.k, pc.j))
  resProductCov <- productcov(resProduct$eig, resProduct$scores)
  resPACE <- FPCA(y, tt, list(useBinnedData = "OFF", FVEthreshold = .9999, dataType = "Dense"))
  yt <- t(x)
  EmpCov <- 1 / (n-1)  * (yt - rowMeans(yt))%*%t(yt - rowMeans(yt))
  EmpMean <- colMeans(x)
  results <- numeric(24)
  results[1] <- infoD$DIC
  results[2] <- infoB$BIC1
  results[3] <- infoB$BIC2
  results[4] <- sum((Smooth_scaled_cov - resProductCov)^2) / sum(Smooth_scaled_cov^2)
  results[5] <- sum((Smooth_scaled_cov - resMarginalCov)^2) / sum(Smooth_scaled_cov^2)
  results[6] <- sum((Smooth_scaled_cov - resPACE$smoothedCov)^2) / sum(Smooth_scaled_cov^2)
  results[7] <- sum((Smooth_scaled_cov - (EmpCov - resPACE$sigma2 * diag(SS * TT)))^2) / sum(Smooth_scaled_cov^2)
  results[8] <- sum((Smooth_scaled_cov - MCMC_eigen$postcov)^2) / sum(Smooth_scaled_cov^2)
  
  m1 <- eigen(Marg.Long)$vectors[,1:3] # Marginal longitudinal eigenvectors
  
  results[9] <- min(sum((resProduct$phi[,1] - m1[,1])^2), sum((resProduct$phi[,1] + m1[,1])^2))
  results[10] <-min(sum((resProduct$phi[,2] - m1[,2])^2), sum((resProduct$phi[,2] + m1[,2])^2))
  results[11] <-min(sum((resProduct$phi[,3] - m1[,3])^2), sum((resProduct$phi[,3] + m1[,3])^2))
  results[12] <-min(sum((MCMC_eigen$eigvecLongmean[,3] - m1[,1])^2), sum((MCMC_eigen$eigvecLongmean[,3] + m1[,1])^2))
  results[13] <-min(sum((MCMC_eigen$eigvecLongmean[,2] - m1[,2])^2), sum((MCMC_eigen$eigvecLongmean[,2] + m1[,2])^2))
  results[14] <-min(sum((MCMC_eigen$eigvecLongmean[,1] - m1[,3])^2), sum((MCMC_eigen$eigvecLongmean[,1] + m1[,3])^2))

  m2 <- eigen(Marg.Func)$vectors[,1:3] # Marginal functional eigenvectors

  results[15] <-min(Re(sum((resProduct$psi[,1] - m2[,1])^2)), Re(sum((resProduct$psi[,1] + m2[,1])^2)))
  results[16] <-min(Re(sum((resProduct$psi[,2] - m2[,2])^2)), Re(sum((resProduct$psi[,2] + m2[,2])^2)))
  results[17] <-min(Re(sum((resProduct$psi[,3] - m2[,3])^2)), Re(sum((resProduct$psi[,3] + m2[,3])^2)))
  results[18] <-min(Re(sum((MCMC_eigen$eigvecFuncmean[,3] - m2[,1])^2)), Re(sum((MCMC_eigen$eigvecFuncmean[,3] + m2[,1])^2)))
  results[19] <-min(Re(sum((MCMC_eigen$eigvecFuncmean[,2] - m2[,2])^2)), Re(sum((MCMC_eigen$eigvecFuncmean[,2] + m2[,2])^2)))
  results[20] <-min(Re(sum((MCMC_eigen$eigvecFuncmean[,1] - m2[,3])^2)), Re(sum((MCMC_eigen$eigvecFuncmean[,1] + m2[,3])^2)))
  
  # Compare marginal covariances
  
  results[21] <- sum((getMarginalLong(MCMC_eigen$postcov,SS, TT) - Marg.Long)^2) / sum(Marg.Long^2)
  results[22] <- sum((getMarginalLong(resProductCov, SS, TT) - Marg.Long)^2) / sum(Marg.Long^2)
  results[23] <- sum((getMarginalFunc(MCMC_eigen$postcov,SS,TT) - Marg.Func)^2) / sum(Marg.Func^2)
  results[24] <- sum((getMarginalFunc(resProductCov,SS,TT) - Marg.Func)^2) / sum(Marg.Func^2)
  
  results[25] <- sum((mu - EmpMean)^2)/ sum(mu^2)
  results[26] <- sum((mu - as.numeric(MCMC_eigen$postmean))^2) / sum(mu^2)
  save(results, file= paste0("splinenum",splinenum,"_20long20splines",index,".RData"))
  gc()
  results
})[3]
stopCluster(cl)
save(simulation, file = "simulation.RData")

splinenum <- 10
Bt1 <- bs(t, df = splinenum, intercept = TRUE)
Bs1 <- bs(s, df = splinenum, intercept = TRUE)
set.seed(1)
setwd("/Users/John/Documents/Johnstuff/splines")
iterations <- 1000
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cl <- makeCluster(5)
registerDoSNOW(cl)
system.time(MCMC10<-foreach(index=1:iterations,.combine=rbind, .packages = c("MASS", "LFBayes", "fdapace"), .options.snow = opts)%dopar%{
  print(index)
  #x <- mvrnorm(n, mu  = rep(0, TT*SS), Sigma = Cov.Weak)
  x <- mvrnorm(n, mu  = as.vector(mu2), Sigma = Cov)
  sx <- sd(x)
  mx <- mean(x)
  x <- (x-mx)/sx
  Smooth_scaled_cov <- (Cov - errorvar * diag(SS * TT)) / sx^2
  mu <- (mu2 - mx)/sx
  Marg.Long <- getMarginalLong(Smooth_scaled_cov,SS,TT)
  Marg.Func <- getMarginalFunc(Smooth_scaled_cov,SS,TT)
  #scaled_cov <- Cov/sx^2
  y <- lapply(1:n, function(i) x[i,])
  missing <- list()
  for(i in 1:n){
    missing[[i]] <- numeric(0)
  }
  X <- rep(1,n)
  dim(X) <- c(n,1)
  
  q1s <- 4
  q2s <- 4
  mcmc <- mcmcWeakChains(y, missing, X, Bs1, Bt1, q1s, q2s, 5000, 1, 0,1)
  infoD <- calculate_DIC(t(x), X, mcmc, Bs1, Bt1, 1000, 10)
  infoB <- calculate_BIC(t(x), X, mcmc, Bs1, Bt1, 1000, 10)
  mcmc$postcov <- infoD$postcov
  mcmc$postmean <- infoD$meanvec
  #resMarginal <- MarginalFPCA(x, n, SS, TT, fpca.op1, fpca.op2, pc.j, rep(pc.k, pc.j))
  #resMarginalCov <- marginalcov(resMarginal$eig, resMarginal$scores)
  #resProduct <- ProductFPCA(x, n, SS, TT, fpca.op1, fpca.op2, pc.j, rep(pc.k, pc.j))
  #resProductCov <- productcov(resProduct$eig, resProduct$scores)
  #resPACE <- FPCA(y, tt, list(useBinnedData = "OFF", FVEthreshold = .9999, dataType = "Dense"))
  #yt <- t(x)
  #EmpCov <- 1 / (n-1)  * (yt - rowMeans(yt))%*%t(yt - rowMeans(yt))
  #EmpMean <- colMeans(x)
  results <- numeric(24)
  results[1] <- infoD$DIC
  results[2] <- infoB$BIC1
  results[3] <- infoB$BIC2
  #results[1] <- sum((Smooth_scaled_cov - resProductCov)^2) / sum(Smooth_scaled_cov^2)
  #results[2] <- sum((Smooth_scaled_cov - resMarginalCov)^2) / sum(Smooth_scaled_cov^2)
  #results[3] <- sum((Smooth_scaled_cov - resPACE$smoothedCov)^2) / sum(Smooth_scaled_cov^2)
  #results[4] <- sum((Smooth_scaled_cov - (EmpCov - resPACE$sigma2 * diag(SS * TT)))^2) / sum(Smooth_scaled_cov^2)
  results[5] <- sum((Smooth_scaled_cov - mcmc$postcov)^2) / sum(Smooth_scaled_cov^2)
  #m1 <- eigen(Marg.Long)$vectors[,1:3]
  #results[7] <- min(sum((resProduct$phi[,1] - m1[,1])^2), sum((resProduct$phi[,1] + m1[,1])^2))
  #results[8] <-min(sum((resProduct$phi[,2] - m1[,2])^2), sum((resProduct$phi[,2] + m1[,2])^2))
  #results[9] <-min(sum((resProduct$phi[,3] - m1[,3])^2), sum((resProduct$phi[,3] + m1[,3])^2))
  #results[10] <-min(sum((mcmc$eigvecLongmean[,3] - m1[,1])^2), sum((mcmc$eigvecLongmean[,3] + m1[,1])^2))
  #results[11] <-min(sum((mcmc$eigvecLongmean[,2] - m1[,2])^2), sum((mcmc$eigvecLongmean[,2] + m1[,2])^2))
  #results[12] <-min(sum((mcmc$eigvecLongmean[,1] - m1[,3])^2), sum((mcmc$eigvecLongmean[,1] + m1[,3])^2))
  
  #m2 <- eigen(Marg.Func)$vectors[,1:3]
  
  #results[13] <-min(Re(sum((resProduct$psi[,1] - m2[,1])^2)), Re(sum((resProduct$psi[,1] + m2[,1])^2)))
  #results[14] <-min(Re(sum((resProduct$psi[,2] - m2[,2])^2)), Re(sum((resProduct$psi[,2] + m2[,2])^2)))
  #results[15] <-min(Re(sum((resProduct$psi[,3] - m2[,3])^2)), Re(sum((resProduct$psi[,3] + m2[,3])^2)))
  #results[16] <-min(Re(sum((mcmc$eigvecFuncmean[,3] - m2[,1])^2)), Re(sum((mcmc$eigvecFuncmean[,3] + m2[,1])^2)))
  #results[17] <-min(Re(sum((mcmc$eigvecFuncmean[,2] - m2[,2])^2)), Re(sum((mcmc$eigvecFuncmean[,2] + m2[,2])^2)))
  #results[18] <-min(Re(sum((mcmc$eigvecFuncmean[,1] - m2[,3])^2)), Re(sum((mcmc$eigvecFuncmean[,1] + m2[,3])^2)))
  
  results[19] <- sum((getMarginalLong(mcmc$postcov,SS, TT) - Marg.Long)^2) / sum(Marg.Long^2)
  #results[20] <- sum((getMarginalLong(resProductCov, SS, TT) - Marg.Long)^2) / sum(Marg.Long^2)
  results[21] <- sum((getMarginalFunc(mcmc$postcov,SS,TT) - Marg.Func)^2) / sum(Marg.Func^2)
  #results[22] <- sum((getMarginalFunc(resProductCov,SS,TT) - Marg.Func)^2) / sum(Marg.Func^2)
  
  #results[23] <- sum((mu - resPACE$mu)^2)/ sum(mu^2)
  results[24] <- sum((mu - as.numeric(mcmc$postmean))^2) / sum(mu^2)
  save(results, file= paste0("splinenum",splinenum,"_20long20splines",index,".RData"))
  gc()
  results
})[3]
stopCluster(cl)
save(MCMC10, file = "MCMC10.RData")

splinenum <- 15
Bt1 <- bs(t, df = splinenum, intercept = TRUE)
Bs1 <- bs(s, df = splinenum, intercept = TRUE)
set.seed(1)
setwd("/Users/John/Documents/Johnstuff/splines")
iterations <- 1000
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cl <- makeCluster(5)
registerDoSNOW(cl)
system.time(MCMC15<-foreach(index=1:iterations,.combine=rbind, .packages = c("MASS", "LFBayes", "fdapace"), .options.snow = opts)%dopar%{
  print(index)
  #x <- mvrnorm(n, mu  = rep(0, TT*SS), Sigma = Cov.Weak)
  x <- mvrnorm(n, mu  = as.vector(mu2), Sigma = Cov)
  sx <- sd(x)
  mx <- mean(x)
  x <- (x-mx)/sx
  Smooth_scaled_cov <- (Cov - errorvar * diag(SS * TT)) / sx^2
  mu <- (mu2 - mx)/sx
  Marg.Long <- getMarginalLong(Smooth_scaled_cov,SS,TT)
  Marg.Func <- getMarginalFunc(Smooth_scaled_cov,SS,TT)
  #scaled_cov <- Cov/sx^2
  y <- lapply(1:n, function(i) x[i,])
  missing <- list()
  for(i in 1:n){
    missing[[i]] <- numeric(0)
  }
  X <- rep(1,n)
  dim(X) <- c(n,1)
  
  q1s <- 4
  q2s <- 4
  mcmc <- mcmcWeakChains(y, missing, X, Bs1, Bt1, q1s, q2s, 5000, 1, 0,1)
  infoD <- calculate_DIC(t(x), X, mcmc, Bs1, Bt1, 1000, 10)
  infoB <- calculate_BIC(t(x), X, mcmc, Bs1, Bt1, 1000, 10)
  mcmc$postcov <- infoD$postcov
  mcmc$postmean <- infoD$meanvec
  #resMarginal <- MarginalFPCA(x, n, SS, TT, fpca.op1, fpca.op2, pc.j, rep(pc.k, pc.j))
  #resMarginalCov <- marginalcov(resMarginal$eig, resMarginal$scores)
  #resProduct <- ProductFPCA(x, n, SS, TT, fpca.op1, fpca.op2, pc.j, rep(pc.k, pc.j))
  #resProductCov <- productcov(resProduct$eig, resProduct$scores)
  #resPACE <- FPCA(y, tt, list(useBinnedData = "OFF", FVEthreshold = .9999, dataType = "Dense"))
  #yt <- t(x)
  #EmpCov <- 1 / (n-1)  * (yt - rowMeans(yt))%*%t(yt - rowMeans(yt))
  #EmpMean <- colMeans(x)
  results <- numeric(24)
  results[1] <- infoD$DIC
  results[2] <- infoB$BIC1
  results[3] <- infoB$BIC2
  #results[1] <- sum((Smooth_scaled_cov - resProductCov)^2) / sum(Smooth_scaled_cov^2)
  #results[2] <- sum((Smooth_scaled_cov - resMarginalCov)^2) / sum(Smooth_scaled_cov^2)
  #results[3] <- sum((Smooth_scaled_cov - resPACE$smoothedCov)^2) / sum(Smooth_scaled_cov^2)
  #results[4] <- sum((Smooth_scaled_cov - (EmpCov - resPACE$sigma2 * diag(SS * TT)))^2) / sum(Smooth_scaled_cov^2)
  results[5] <- sum((Smooth_scaled_cov - mcmc$postcov)^2) / sum(Smooth_scaled_cov^2)
  #m1 <- eigen(Marg.Long)$vectors[,1:3]
  #results[7] <- min(sum((resProduct$phi[,1] - m1[,1])^2), sum((resProduct$phi[,1] + m1[,1])^2))
  #results[8] <-min(sum((resProduct$phi[,2] - m1[,2])^2), sum((resProduct$phi[,2] + m1[,2])^2))
  #results[9] <-min(sum((resProduct$phi[,3] - m1[,3])^2), sum((resProduct$phi[,3] + m1[,3])^2))
  #results[10] <-min(sum((mcmc$eigvecLongmean[,3] - m1[,1])^2), sum((mcmc$eigvecLongmean[,3] + m1[,1])^2))
  #results[11] <-min(sum((mcmc$eigvecLongmean[,2] - m1[,2])^2), sum((mcmc$eigvecLongmean[,2] + m1[,2])^2))
  #results[12] <-min(sum((mcmc$eigvecLongmean[,1] - m1[,3])^2), sum((mcmc$eigvecLongmean[,1] + m1[,3])^2))
  
  #m2 <- eigen(Marg.Func)$vectors[,1:3]
  
  #results[13] <-min(Re(sum((resProduct$psi[,1] - m2[,1])^2)), Re(sum((resProduct$psi[,1] + m2[,1])^2)))
  #results[14] <-min(Re(sum((resProduct$psi[,2] - m2[,2])^2)), Re(sum((resProduct$psi[,2] + m2[,2])^2)))
  #results[15] <-min(Re(sum((resProduct$psi[,3] - m2[,3])^2)), Re(sum((resProduct$psi[,3] + m2[,3])^2)))
  #results[16] <-min(Re(sum((mcmc$eigvecFuncmean[,3] - m2[,1])^2)), Re(sum((mcmc$eigvecFuncmean[,3] + m2[,1])^2)))
  #results[17] <-min(Re(sum((mcmc$eigvecFuncmean[,2] - m2[,2])^2)), Re(sum((mcmc$eigvecFuncmean[,2] + m2[,2])^2)))
  #results[18] <-min(Re(sum((mcmc$eigvecFuncmean[,1] - m2[,3])^2)), Re(sum((mcmc$eigvecFuncmean[,1] + m2[,3])^2)))
  
  results[19] <- sum((getMarginalLong(mcmc$postcov,SS, TT) - Marg.Long)^2) / sum(Marg.Long^2)
  #results[20] <- sum((getMarginalLong(resProductCov, SS, TT) - Marg.Long)^2) / sum(Marg.Long^2)
  results[21] <- sum((getMarginalFunc(mcmc$postcov,SS,TT) - Marg.Func)^2) / sum(Marg.Func^2)
  #results[22] <- sum((getMarginalFunc(resProductCov,SS,TT) - Marg.Func)^2) / sum(Marg.Func^2)
  
  #results[23] <- sum((mu - resPACE$mu)^2)/ sum(mu^2)
  results[24] <- sum((mu - as.numeric(mcmc$postmean))^2) / sum(mu^2)
  save(results, file= paste0("splinenum",splinenum,"_20long20splines",index,".RData"))
  gc()
  results
})[3]
stopCluster(cl)
save(MCMC15, file = "MCMC15.RData")


round(mean(n30_50[,20]),3)
short_results <- matrix(0,nrow=10,ncol=20)
for(i in c(1:10)){
  load(paste0("n30_ccc",i,".RData"))
  short_results[i,] <- results
}
short_results<-short_results[c(1:10),]
n30_results <- matrix(0,nrow=250,ncol = 20)
for(i in 1:250){
  load(paste0("n30_",i,".RData"))
  n30_results[i,] <- results
}
n30_results <- matrix(as.numeric(n30_results),nrow=250,ncol=20)

round(mean(n30_results[,5]),3)


n60_results <- matrix(0,nrow=250,ncol = 20)
for(i in c(1:19,21,23:30,32:42,44:209,211:250)){
  load(paste0("n60_",i,".RData"))
  n60_results[i,] <- results
}
n60_results <- matrix(as.numeric(n60_results),nrow=250,ncol=20)
rows <- apply(n60_results,1,function(i) all(i==0))
n60_results<-n60_results[!rows,]
round(mean(n60_results[,5]),3)

n90_results <- matrix(0,nrow=250,ncol = 20)
for(i in c(1:48,50:67,69:87,89:91,93:99,101:118,120:136,138,140:165,
           167:171,173,175:180,182:185,187:199,201,203:207,209:213,215:220,222:250)){
  load(paste0("n90_",i,".RData"))
  n90_results[i,] <- results
}
n90_results <- matrix(as.numeric(n90_results),nrow=250,ncol=20)
rows <- apply(n90_results,1,function(i) all(i==0))
n90_results<-n90_results[!rows,]

round(mean(n90_results[,5]),3)


mlongb <- getMarginalLong(mcmc$postcov,20,20)
mfuncb <- getMarginalFunc(mcmc$postcov,20,20)
mlongt <- getMarginalLong(Smooth_scaled_cov,20,20)
mfunct <- getMarginalFunc(Smooth_scaled_cov,20,20)
mlongf <- getMarginalLong(resProductCov,20,20)
mfuncf <- getMarginalFunc(resProductCov,20,20)
mlonge <- getMarginalLong(EmpCov,20,20)
mfunce <- getMarginalFunc(EmpCov,20,20)
image(mlongt[,20:1],col = heat.colors(100))
image(mlongb[,20:1],col = heat.colors(100))
image(mlongf[,20:1],col = heat.colors(100))
image(mlonge[,20:1],col = heat.colors(100))
plot(eigen(mlongf)$vectors[,1],type="l")
lines(resProduct$phi[,1],col="red")
plot(eigen(mlongf)$vectors[,2],type="l")
lines(-resProduct$phi[,2],col="red")
plot(eigen(mlongf)$vectors[,3],type="l")
lines(resProduct$phi[,3],col="red")
plot(eigen(mfuncf)$vectors[,1],type="l")
lines(-resProduct$psi[,1],col="red")
plot(eigen(mfuncf)$vectors[,2],type="l")
lines(-resProduct$psi[,2],col="red")
plot(eigen(mfuncf)$vectors[,3],type="l")
lines(resProduct$psi[,3],col="red")

j <-11
plot(mlongt[j,],type="l")
lines(mlongb[j,],col="blue")
lines(mlongf[j,],col="red")
lines(mlonge[j,],col="green")

j <- 197
plot(Smooth_scaled_cov[j,],type="l")
lines(mcmc$postcov[j,],col="blue")
lines(resProductCov[j,],col="red")
lines(EmpCov[j,],col="green")

j <- 2
plot(mfunct[j,],type="l")
lines(mfuncb[j,],col="blue")
lines(mfuncf[j,],col="red")
lines(mfunce[j,],col="green")

mylist <- list()
mylist[[1]] <- mfunct
mylist[[2]] <- mfuncb
mylist[[3]] <- mfuncf
mylist[[4]] <- mfunce
image(mfunct, col= heat.colors(100),zlim=c(min(unlist(mylist)),max(unlist(mylist))))
image(mfuncb, col= heat.colors(100),zlim=c(min(unlist(mylist)),max(unlist(mylist))))
image(mfuncf, col= heat.colors(100),zlim=c(min(unlist(mylist)),max(unlist(mylist))))
image(mfunce, col= heat.colors(100),zlim=c(min(unlist(mylist)),max(unlist(mylist))))

image(mlongt, col= heat.colors(100))
image(mlongb, col= heat.colors(100))
image(mlongf, col= heat.colors(100))
image(mlonge, col= heat.colors(100))




errorvar <- .01
Cov.Weak <- kronecker(Bs%*%Gamma, Bt%*%Lambda)%*%H%*%t(kronecker(Bs%*%Gamma, Bt%*%Lambda)) + errorvar*diag(length(s) * length(t))
iterations <- 500
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cl <- makeCluster(6)
registerDoSNOW(cl)
system.time(n90_s.01_Sep<-foreach(index=1:iterations,.combine=rbind, .packages = c("MASS", "LFBayes", "fdapace"), .options.snow = opts)%dopar%{
  print(index)
  #x <- mvrnorm(n, mu  = rep(0, length(t)*length(s)), Sigma = Cov.Weak)
  x <- mvrnorm(30, mu  = as.vector(mu1), Sigma = Cov.Weak)
  sx <- sd(x)
  mx <- mean(x)
  x <- (x-mx)/sx
  Smooth_scaled_cov <- (Cov.Weak - errorvar * diag(length(s) * length(t))) / sx^2
  mu <- (mu1 - mu)/sx
  y <- lapply(1:n, function(i) x[i,])
  missing <- list()
  for(i in 1:n){
    missing[[i]] <- numeric(0)
  }
  X <- rep(1,n)
  dim(X) <- c(n,1)
  
  q <- 8
  mcmc <- mcmcWeak(y, missing, X, Bs1, Bt1, q, q, 25500, 1, 5000)
  
  resMarginal <- MarginalFPCA(x, n, length(s), length(t), fpca.op1, fpca.op2, pc.j, rep(pc.k, pc.j))
  resMarginalCov <- marginalcov(resMarginal$eig, resMarginal$scores)
  resProduct <- ProductFPCA(x, n, length(s), length(t), fpca.op1, fpca.op2, pc.j, rep(pc.k, pc.j))
  resProductCov <- productcov(resProduct$eig, resProduct$scores)
  resPACE <- FPCA(y, tt, list(useBinnedData = "OFF", FVEthreshold = .9999, dataType = "Dense"))
  yt <- t(x)
  EmpCov <- 1 / (n-1)  * (yt - rowMeans(yt))%*%t(yt - rowMeans(yt))
  EmpMean <- colMeans(x)
  results <- numeric(19)
  
  results[1] <- max(abs(eigen(Smooth_scaled_cov - resProductCov)$values)) / (length(s) * length(t))
  results[2] <- max(abs(eigen(Smooth_scaled_cov - resMarginalCov)$values)) / (length(s) * length(t))
  results[3] <- max(abs(eigen(Smooth_scaled_cov - resPACE$smoothedCov)$values)) / (length(s) * length(t))
  results[4] <- max(abs(eigen(Smooth_scaled_cov - (EmpCov - resPACE$sigma2 * diag(length(s) * length(t))))$values)) / (length(s) * length(t))
  results[5] <- max(abs(eigen(Smooth_scaled_cov - mcmc$postcov)$values)) / (length(s) * length(t))
  
  m1 <- eigen(Brown.Motion.Cov)$vectors[,1:3]
  results[6] <- min(sum((resProduct$phi[,1] - m1[,1])^2), sum((resProduct$phi[,1] + m1[,1])^2))
  results[7] <-min(sum((resProduct$phi[,2] - m1[,2])^2), sum((resProduct$phi[,2] + m1[,2])^2))
  results[8] <-min(sum((resProduct$phi[,3] - m1[,3])^2), sum((resProduct$phi[,3] + m1[,3])^2))
  results[9] <-min(sum((mcmc$eigvecLongmean[,3] - m1[,1])^2), sum((mcmc$eigvecLongmean[,3] + m1[,1])^2))
  results[10] <-min(sum((mcmc$eigvecLongmean[,2] - m1[,2])^2), sum((mcmc$eigvecLongmean[,2] + m1[,2])^2))
  results[11] <-min(sum((mcmc$eigvecLongmean[,1] - m1[,3])^2), sum((mcmc$eigvecLongmean[,1] + m1[,3])^2))
  
  m2 <- eigen(Matern.Cov)$vectors[,1:3]
  
  results[12] <-min(sum((resProduct$psi[,1] - m2[,1])^2), sum((resProduct$psi[,1] + m2[,1])^2))
  results[13] <-min(sum((resProduct$psi[,2] - m2[,2])^2), sum((resProduct$psi[,2] + m2[,2])^2))
  results[14] <-min(sum((resProduct$psi[,3] - m2[,3])^2), sum((resProduct$psi[,3] + m2[,3])^2))
  results[15] <-min(sum((mcmc$eigvecFuncmean[,3] - m2[,1])^2), sum((mcmc$eigvecFuncmean[,3] + m2[,1])^2))
  results[16] <-min(sum((mcmc$eigvecFuncmean[,2] - m2[,2])^2), sum((mcmc$eigvecFuncmean[,2] + m2[,2])^2))
  results[17] <-min(sum((mcmc$eigvecFuncmean[,1] - m2[,3])^2), sum((mcmc$eigvecFuncmean[,1] + m2[,3])^2))
  
  results[18] <- sum((mu - resPACE$mu)^2)/(length(s) * length(t))
  results[19] <- sum((mu - as.numeric(mcmc$postmean))^2)/(length(s) * length(t))
  gc()
  results
})[3]
stopCluster(cl)
save(n90_s.01_Sep, file = "n90_s.01_Sep.RData")


errorvar <- .001
Cov.Weak <- kronecker(Bs%*%Gamma, Bt%*%Lambda)%*%H%*%t(kronecker(Bs%*%Gamma, Bt%*%Lambda)) + errorvar*diag(length(s) * length(t))
iterations <- 500
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cl <- makeCluster(6)
registerDoSNOW(cl)
system.time(n90_s.001_Sep<-foreach(index=1:iterations,.combine=rbind, .packages = c("MASS", "LFBayes", "fdapace"), .options.snow = opts)%dopar%{
  print(index)
  #x <- mvrnorm(n, mu  = rep(0, length(t)*length(s)), Sigma = Cov.Weak)
  x <- mvrnorm(n, mu  = as.vector(mu1), Sigma = Cov.Weak)
  sx <- sd(x)
  mx <- mean(x)
  x <- (x-mx)/sx
  Smooth_scaled_cov <- (Cov.Weak - errorvar * diag(length(s) * length(t))) / sx^2
  mu <- (mu1 - mx)/sx
  y <- lapply(1:n, function(i) x[i,])
  missing <- list()
  for(i in 1:n){
    missing[[i]] <- numeric(0)
  }
  X <- rep(1,n)
  dim(X) <- c(n,1)
  
  q <- 8
  mcmc <- mcmcWeak(y, missing, X, Bs1, Bt1, q, q, 25000, 1, 5000)
  
  resMarginal <- MarginalFPCA(x, n, length(s), length(t), fpca.op1, fpca.op2, pc.j, rep(pc.k, pc.j))
  resMarginalCov <- marginalcov(resMarginal$eig, resMarginal$scores)
  resProduct <- ProductFPCA(x, n, length(s), length(t), fpca.op1, fpca.op2, pc.j, rep(pc.k, pc.j))
  resProductCov <- productcov(resProduct$eig, resProduct$scores)
  resPACE <- FPCA(y, tt, list(useBinnedData = "OFF", FVEthreshold = .9999, dataType = "Dense"))
  yt <- t(x)
  EmpCov <- 1 / (n-1)  * (yt - rowMeans(yt))%*%t(yt - rowMeans(yt))
  EmpMean <- colMeans(x)
  results <- numeric(19)
  
  results[1] <- max(abs(eigen(Smooth_scaled_cov - resProductCov)$values)) / (length(s) * length(t))
  results[2] <- max(abs(eigen(Smooth_scaled_cov - resMarginalCov)$values)) / (length(s) * length(t))
  results[3] <- max(abs(eigen(Smooth_scaled_cov - resPACE$smoothedCov)$values)) / (length(s) * length(t))
  results[4] <- max(abs(eigen(Smooth_scaled_cov - (EmpCov - resPACE$sigma2 * diag(length(s) * length(t))))$values)) / (length(s) * length(t))
  results[5] <- max(abs(eigen(Smooth_scaled_cov - mcmc$postcov)$values)) / (length(s) * length(t))
  
  m1 <- eigen(Brown.Motion.Cov)$vectors[,1:3]
  results[6] <- min(sum((resProduct$phi[,1] - m1[,1])^2), sum((resProduct$phi[,1] + m1[,1])^2))
  results[7] <-min(sum((resProduct$phi[,2] - m1[,2])^2), sum((resProduct$phi[,2] + m1[,2])^2))
  results[8] <-min(sum((resProduct$phi[,3] - m1[,3])^2), sum((resProduct$phi[,3] + m1[,3])^2))
  results[9] <-min(sum((mcmc$eigvecLongmean[,3] - m1[,1])^2), sum((mcmc$eigvecLongmean[,3] + m1[,1])^2))
  results[10] <-min(sum((mcmc$eigvecLongmean[,2] - m1[,2])^2), sum((mcmc$eigvecLongmean[,2] + m1[,2])^2))
  results[11] <-min(sum((mcmc$eigvecLongmean[,1] - m1[,3])^2), sum((mcmc$eigvecLongmean[,1] + m1[,3])^2))
  
  m2 <- eigen(Matern.Cov)$vectors[,1:3]
  
  results[12] <-min(sum((resProduct$psi[,1] - m2[,1])^2), sum((resProduct$psi[,1] + m2[,1])^2))
  results[13] <-min(sum((resProduct$psi[,2] - m2[,2])^2), sum((resProduct$psi[,2] + m2[,2])^2))
  results[14] <-min(sum((resProduct$psi[,3] - m2[,3])^2), sum((resProduct$psi[,3] + m2[,3])^2))
  results[15] <-min(sum((mcmc$eigvecFuncmean[,3] - m2[,1])^2), sum((mcmc$eigvecFuncmean[,3] + m2[,1])^2))
  results[16] <-min(sum((mcmc$eigvecFuncmean[,2] - m2[,2])^2), sum((mcmc$eigvecFuncmean[,2] + m2[,2])^2))
  results[17] <-min(sum((mcmc$eigvecFuncmean[,1] - m2[,3])^2), sum((mcmc$eigvecFuncmean[,1] + m2[,3])^2))
  
  results[18] <- sum((mu - resPACE$mu)^2)/(length(s) * length(t))
  results[19] <- sum((mu - as.numeric(mcmc$postmean))^2)/(length(s) * length(t))
  gc()
  results
})[3]
stopCluster(cl)
save(n90_s.001_Sep, file = "n90_s.001_Sep.RData")

halft <- function(sig){
  nu <- 1
  A <- 10
  (1+1/nu*(sig/A)^2)^(-(nu+1)/2)
}



A <- diag(15)
Xsp <- matrix(nrow = 900, ncol = 20)
for(j in 1:60){
  for(i in 1:15){
    print((j-1)*15+i)
    Xsp[(j-1) * 15 + i,] <- Bt1 %*% mcmc$theta[[50]][,,j] %*% t(Bs1) %*% A[,i]
  }
}
plot(x[2,],type="l")
lines(21:40, Xsp[17,],col="red")

evecsp <- eigen(cov(Xsp))$vectors
plot(evecemp[,1],type="l")
lines(resProduct$psi[,1],col="red")
lines(evecsp[,1],col="blue")
lines(evecth[,1],col="green")
plot(evecemp[,2],type="l")
lines(resProduct$psi[,2],col="red")
lines(evecsp[,2],col="blue")
lines(evecth[,2],col="green")
plot(evecemp[,3],type="l")
lines(resProduct$psi[,3],col="red")
lines(evecsp[,3],col="blue")
lines(-evecth[,3],col="green")
