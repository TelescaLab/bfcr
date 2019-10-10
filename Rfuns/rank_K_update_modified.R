splines <- "
	data{
    int<lower=0> N; // Number of samples
    int<lower=0> num_basis; // Dimension of basis matrix
    int<lower=0> t; // number of time points
    vector[t * N] Y;
    matrix[N,2] X;
    matrix[t, num_basis] B;
    int<lower=0> K;
  }
  parameters{
    matrix[num_basis,2] a;
    real<lower=0> precision;
    vector[N] Z[K];
    matrix[num_basis, 2] Lam[K];
    matrix<lower=0>[2, 2] smooth_inv;
    vector<lower=0>[num_basis] phi;
    vector[N * num_basis] projectionp;
  }
  transformed parameters{
    vector[t * N] Y_hat;
    real<lower=0> sigma;
    vector<lower=0>[num_basis] phi_inv_sqrt;
    matrix<lower=0>[2, 2] smooth;
    vector[num_basis * N] projection;
    projection = to_vector(a*X');
    for(k in 1:K){
      projection = projection + to_vector(Lam[k]*X'*diag_matrix(Z[k]));
    }
    sigma = inv(sqrt(precision));
    for(j in 1:num_basis){
      phi_inv_sqrt[j] = 1 / sqrt(phi[j]);
    }
    //phi_inv_sqrt = 1 ./ sqrt(phi);
    smooth = 1 ./ smooth_inv;
  }
  model{
    precision ~ gamma(.0001, .0001);
    to_vector(phi) ~ gamma(1, .5);
    for(k in 1:K){
      Z[k] ~ normal(0,1);
    }
    to_vector(smooth_inv) ~ gamma(1, .00005);
    for(k in 1:K){
      for(i in 1:num_basis){
        for(j in 1:2){
          if(i == 1){
            Lam[k, i, j] ~ normal(0, smooth[2,j]);
            a[i, j] ~ normal(0, smooth[1,j]);
          }
          else if(i == 2){
            Lam[k, i, j] ~ normal(2*Lam[k, i-1,j], smooth[2,j]);
            a[i, j] ~ normal(2*a[i-1,j], smooth[1,j]);
          }
          else {
            Lam[k, i, j] ~ normal(2*Lam[k, i-1,j] - Lam[k,i-2,j], smooth[2,j]);
            a[i, j] ~ normal(2*a[i-1,j] - a[i-2,j], smooth[1,j]);
    
          }
        }
      }
    }
    for(i in 1:N){
      for(j in 1:num_basis){
        projectionp[((i-1) * num_basis + j)] ~ normal(projection[((i-1) * num_basis + j)], phi_inv_sqrt[j]);
      }
    }
    for(i in 1:N){
      Y[(i-1) * t + 1:(t * i)] ~ normal(B * projectionp[((i-1) * num_basis + 1):(num_basis * i)], sigma);
    }
}
"

library(dlnm)
library(rstan)
options(mc.cores =parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=core2')
set.seed(3)
K <- 2
tmax <- 50
N = 200
num_basis_true <- 12
sigma = .0001
tp <- seq(from=0, to=1, length.out = tmax) # generating inputs
B_true <- bs(tp, df = num_basis_true, intercept = TRUE)
#alpha <- matrix(rnorm(num_basis_true*2), nrow = num_basis_true, ncol = 2)
alpha <- Theta1
X <- cbind(rep(1,N), rnorm(N))
Z <- matrix(rnorm(N*2), ncol = 2)
#L <- matrix(rnorm(num_basis_true*2), num_basis_true, 2)
#Y_true <- B_true%*%alpha%*%t(X) + B_true%*%L%*%t(X)%*%diag(Z)
Y_true <- B_true%*%alpha%*%t(X) + B_true%*%Lambda1%*%t(X)%*%diag(Z[,1]) + B_true%*%Lambda2%*%t(X)%*%diag(Z[,2])
Y <- Y_true + matrix(rnorm(N * length(t)), nrow = length(tp), ncol = N) * sigma
Ycol <- c(Y)
plot(Y[,1],type="p")
plot(Y[,2], type = "p")

#data = list(N = N, num_basis = num_basis, Y = Y, X = X, B = B)
num_basis_est <- 15
B_est <- ps(tp, df = num_basis_est, intercept = TRUE)
data1 = list(N = N, num_basis = num_basis_est, t = tmax, Y = Ycol, X = X, B = B_est, K = 2)

reg <- lm(Ycol ~ kronecker(X, B_est)-1)
init <- list(list(a = matrix(reg$coefficients, ncol = 2)))
fitvb <- vb(m, data = data1, algorithm = "meanfield",  init = init[[1]], iter = 200000, par = c("a", "Lam", "Z", "sigma", "smooth_inv", "phi_inv_sqrt", "projectionp"))
posterior <- extract(fitvb)


#fit <- stan(file = "stan_file.stan", data = data1, cores = 1)

m <- stan_model(model_code = splines)
fit <- sampling(m, data = data1, chains=1, iter=1000, thin = 1, par = c("a", "Lam", "Z", "sigma"), init = init2, warmup = 500)

posterior <- extract(fit)
dim(posterior$Lam)

subj <- 5
plot(tp,Ycol[((subj-1) * tmax + 1):(subj*tmax)])
for(i in 1:1000){
  y_hat <- B_est%*%posterior$a[i,,]%*%X[subj,]
  for(k in 1:K){
    y_hat <- y_hat + B_est%*%posterior$Lam[i,k,,]%*%X[subj,]*posterior$Z[i,k,subj]
  }
  lines(tp, y_hat, col="gray")
}
lines(tp,c(Y_true)[((subj-1) * tmax + 1):(subj*tmax)], col = "red")

subj <- 1
plot(tp,Ycol[((subj-1) * tmax + 1):(subj*tmax)])
for(i in 1:1000){
  lines(tp, B_est%*%posterior$projectionp[i, ((subj-1) * num_basis_est + 1):(subj*num_basis_est)], col = "gray")
}
lines(tp,c(Y_true)[((subj-1) * tmax + 1):(subj*tmax)], col = "red")


x <- c(1,.5)
cov_true <- B_true%*%Lambda1%*%outer(x,x)%*%t(Lambda1)%*%t(B_true) + B_true%*%Lambda2%*%outer(x,x)%*%t(Lambda2)%*%t(B_true)
cov_est <- matrix(0,nrow = length(tp), ncol = length(tp))
for(i in 1:500){
  for(k in 1:K){
    cov_est <- B_est%*%posterior$Lam[i,k,,]%*%outer(x,x)%*%t(posterior$Lam[i,k,,])%*%t(B_est) + cov_est
  }
}
cov_est <- cov_est / 500
 cov1 <- matrix(0, nrow = tmax, ncol = tmax)
 j <- 500
 for(k in 1:K){
   cov1 <- B_est%*%posterior$Lam[j,k,,]%*%outer(x,x)%*%t(posterior$Lam[j,k,,])%*%t(B_est) + cov1
 }
library(plot3D)
par(mfrow = c(1,3))
persp3D(tp,tp, cov_est, phi = 10, theta = 90, zlim = c(min(cov_true), max(cov_true)))
persp3D(tp,tp, cov_true, phi = 10, theta = 90, zlim = c(min(cov_true), max(cov_true)))
persp3D(tp,tp, cov1, phi = 10, theta = 90, zlim = c(min(cov_true), max(cov_true)))
dev.off()
plot(B_true%*%alpha%*%x,type="l",col="red")
for(i in 1:500){
  lines(B_est%*%posterior$a[i,,]%*%x,col="gray")
  
}
lines(B_true%*%alpha%*%x,type="l",col="red")

data1 = list(N = N, num_basis = num_basis_est, t = tmax, Y = Ycol, X = X, B = B_est, K = 2)
fitvb <- vb(m, data = data1, algorithm = "meanfield",  init = init[[1]], iter = 200000, par = c("a", "Lam", "Z", "sigma", "smooth_inv"))
posterior <- extract(fitvb)
init2 <- list(list(a = apply(posterior$a, c(2,3), mean), Lam = apply(posterior$Lam, c(2,3,4), mean), smooth_inv = apply(posterior$smooth_inv, c(2,3), mean), sigma = mean(posterior$sigma), Z = apply(posterior$Z, c(2,3), mean)))

dim(summary(fitvb)$c_summary)
avb <- t(matrix((summary(fitvb)$c_summary[1:24,1,1]), nrow = 2))
Lamvb <- t(matrix((summary(fitvb)$c_summary[25:48,1,1]), nrow = 2))
Zvb <- summary(fitvb)$c_summary[49:148,1,1]
sigmavb <- summary(fitvb)$c_summary[149,1,1] 

#init <- list(a = avb, Lam = Lamvb, Z = Zvb, sigma = sigmavb))
fitopt <- optimizing(m, data = data1, init = init[[1]], as_vector = FALSE)
fitopt$value
#plot(B_true%*%alpha%*%x,type="l",col="red")
lines(B_est%*%fitopt$par$a%*%x)

x <- c(1,0)
cov_true <- B_true%*%L%*%outer(x,x)%*%t(L)%*%t(B_true)
cov_est <- matrix(0, nrow = tmax, ncol = tmax)
for(k in 1:K){
  cov_est <- cov_est + B_est%*%fitopt$par$Lam[k,,]%*%outer(x,x)%*%t(fitopt$par$Lam[k,,])%*%t(B_est)
}


library(plot3D)
par(mfrow = c(1,2))
persp3D(tp,tp, cov_est, phi = 10, zlim = c(min(cov_true), max(cov_true)))
persp3D(tp,tp, cov_true, phi = 10, zlim = c(min(cov_true), max(cov_true)))
persp3D(tp, tp,B_est%*%posterior$Lam[i,,]%*%outer(x,x)%*%t(posterior$Lam[1,,])%*%t(B_est), zlim = c(min(cov_true), max(cov_true)), phi = 10)
dev.off()
plot(B_true%*%alpha%*%x,type="l",col="red")
lines(B_est%*%fitopt$par$a%*%x)
plot(tp, Ycol[1:100])
Y_hat <- B_est%*%fitopt$par$a%*%X[1,]
for(k in 1:K){
  Y_hat <- Y_hat + B_est%*%fitopt$par$Lam[k,,]%*%X[1,]%*%fitopt$par$Z[k,1]
}
lines(tp, Y_hat)
lines(tp,Y_true[,1],col="red")
