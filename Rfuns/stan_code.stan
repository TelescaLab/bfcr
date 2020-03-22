data {
  int<lower=0> N; // Number of samples
  int<lower=0> P; // Dimension of basis
  int<lower=0> Q; // Dimension of loading matrix
  int<lower=0> D; // Number of covariates
  int<lower=0> tmax;
  matrix[N, D] X;
  vector[tmax * N] Y;
  matrix[tmax, P] B;
}
parameters {
  matrix[P, N] Theta;
  matrix[P, D] Lambda[Q]; // Coefficient loading matrix
  vector<lower=0>[P] Phi_t;
  matrix[P, D] Beta; // Coefficient mean matrix
  matrix[N, Q] Eta; //Latent cov random effect
  //real<lower=0> sigma; // error variance
  vector<lower=0>[D] tau_beta;
  vector<lower=0>[D] tau_lambda;
  //vector[D] SigBeta;
  //vector<lower=0>[N] sigma;
  vector<lower=0>[N] Gamma;
  real<lower=0> alpha;
  //real<lower=0> Phiii;
}
transformed parameters {
  matrix[P, N] Theta_hat;
  matrix[tmax, N] y_hat;
  vector[N] sigma_hat;
  //vector[N] alpha;
  //vector[N] beta;
  for (i in 1:N){
    Theta_hat[1:P, i] = Beta * X[i]';
	for (q in 1:Q){
	  Theta_hat[1:P, i] = col(Theta_hat, i) + Lambda[q] * X[i]' * Eta[i, q];
	}
  }
  for (i in 1:N) {
    y_hat[1:tmax, i] = B * col(Theta, i);
    //sigma_hat[i] = exp(X[i] * SigBeta);
    sigma_hat[i] = alpha^2 * Gamma[i];
  }
  //alpha = sigma_hat .* sigma_hat ./ Phiii;
  //beta = sigma_hat ./ Phiii; 
}
model {
  for (i in 1:N) {
    Y[((i-1)*tmax + 1):(i*tmax)] ~ normal(y_hat[1:tmax, i], sigma_hat[i]);
  }
  //sigma ~ inv_gamma(.001,.001);
  for (i in 1:N){
    Theta[1:P, i] ~ normal(col(Theta_hat, i), Phi_t);
  }
  Phi_t ~ inv_gamma(.001,.001);
  Gamma ~ inv_gamma(1.0/2.0, 1.0 * .0001 / 2.0);
  //Phiii ~ inv_gamma(.01, .01);
  for (q in 1:Q){
    for(i in 1:N){
	Eta[i, q] ~ normal(0, 1);
     }
  }
  for (d in 1:D) {
    for (p in 2:P) {
	Beta[p, d] ~ normal(Beta[p - 1, d], tau_beta[d]);
    }
  }
  for (q in 1:Q) {
    for (d in 1:D) {
      for (p in 2:P) {
	Lambda[q][p, d] ~ normal(Lambda[q][p - 1, d], tau_lambda[d]);
      }
    }
  }
  tau_beta ~ inv_gamma(1, .005);
  tau_lambda ~ inv_gamma(1, .005);
}
generated quantities {
  vector[tmax * N] log_lik;
  int idx;
  for (i in 1:N) {
    for (j in 1:tmax) {
	  idx = (i - 1) * tmax + j;
	  log_lik[idx] = normal_lpdf(Y[idx] | y_hat[j, i], sigma_hat[i]);
	}
  }
}
