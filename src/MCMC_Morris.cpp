#include <RcppArmadillo.h>
// #ifdef _OPENMP
// #include <omp.h>
// #endif
#include "updateParam.h"
#include "Utility.h"
// [[Rcpp::plugins(openmp)]]
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

arma::vec Psi2;
arma::cube Lambda2;
arma::mat Beta2;
arma::mat Eta2;
arma::vec Prec2;
arma::mat BtB2;
arma::mat BtY2;
arma::mat fit2;
arma::vec Phi2;
arma::vec Tau2;
double Tausq2;
double alpha2;
double Nu2;
double b2;
arma::mat P2;
arma::cube log_lik2;

void updatePsi2(){
  double deg_freedom_a = 1;
  double deg_freedom_b = .0005;
  double n = Eta2.n_rows;
  arma::uword K = Eta2.n_cols;
  for(arma::uword k = 0; k < K; k++){
    Psi2(k) = R::rgamma(deg_freedom_a / 2 + n, 1.0 / (deg_freedom_b + 1.0 / 2.0 *
      arma::accu(arma::square(Eta2.col(k)))));
  }
}
void updatePrec2(arma::mat& Y){
  double a = 0.001;
  double b = 0.001;
  //double my_sum = arma::sum(arma::square(Y.t() - B * Proj.t()));
  double my_sum = arma::accu(arma::square(Y - fit2));
  double p = R::rgamma(a + Y.n_elem/2, 1.0 / (b + 1.0/2.0 * my_sum));
  Prec2.fill(p);
}

void updateAlphaBeta2(arma::mat& Y){
  // Phi2 has precision parameterization
  // Output is on the precision scale
  arma::vec norm_vec(Y.n_rows);

  for(arma::uword i = 0; i < Y.n_rows; i++){
    norm_vec(i) = arma::sum(arma::square(Y.row(i) - fit2.row(i))) * Phi2(i);
  }
  alpha2 = R::rgamma(b2 * (Y.n_elem)/2 + 1 - b2, 1.0 / (b2 * (1.0/2.0 * arma::sum(norm_vec))));
}
void updatePhiBeta2(arma::mat& Y){
  // alpha is on precision scale
  // output on precision scale
  for(arma::uword i = 0; i < Phi2.n_elem; i++){
    Phi2(i) = R::rgamma(b2 * (Nu2 + Y.n_cols) / 2 + 1 - b2,
        1.0 / (b2 * (Nu2 * Tausq2 + arma::dot(Y.row(i) -
          fit2.row(i), Y.row(i) - fit2.row(i)) * alpha2) / 2));
  }
}
void updateEtaBeta2(arma::mat Y, arma::mat& X, arma::mat& Z, arma::mat B){
  arma::uword n = X.n_rows;
  arma::mat Xtilde(B.n_cols, Lambda2.n_slices);
  for(arma::uword i = 0; i < n; i++){
    for(arma::uword k = 0; k < Lambda2.n_slices; k++){
      Xtilde.col(k) = Lambda2.slice(k) * Z.row(i).t();
    }
    arma::vec g = Prec2(i) * Xtilde.t() * (B.t() * Y.row(i).t() - BtB2 * Beta2 * X.row(i).t());
    
    //arma::mat C = (Prec2(i) * Xtilde.t() * BtB2 * Xtilde + arma::eye(Lambda2.n_slices, Lambda2.n_slices));
    arma::mat C = (Prec2(i) * Xtilde.t() * BtB2 * Xtilde + arma::diagmat(Psi2));
    arma::mat mychol = arma::chol(C, "lower");
    arma::vec w = solve(arma::trimatl(mychol), g);
    arma::vec mu = solve(arma::trimatu(mychol.t()), w);
    arma::vec z = 1/sqrt(b2) * arma::randn<arma::vec>(Lambda2.n_slices);
    arma::vec v = arma::solve(arma::trimatu(mychol.t()), z);
    Eta2.row(i) = (mu + v).t();
  }
}
void updateTauBeta2(){
  double t_alpha = 1;
  double t_beta =  .0005;
  //double t_beta = 1;
  arma::uword K = Lambda2.n_slices;
  arma::uword D1 = Beta2.n_cols;
  arma::uword D2 = Lambda2.n_cols;
  
  for(arma::uword d = 0; d < D1; d++){
    Tau2(d) = R::rgamma(b2 * (t_alpha + (Lambda2.n_rows - 1.0) / 2.0) + 1 - b2, 1.0 / (b2 * (t_beta + 1.0 / 2.0 * arma::as_scalar(Beta2.col(d).t() * P2 * Beta2.col(d)))));
  }
  double temp_beta;
  for(arma::uword d = 0; d < D2; d++){
    temp_beta = 0;
    for(arma::uword k = 0; k < K; k++){
      temp_beta = temp_beta + 1.0 / 2.0 * arma::as_scalar(Lambda2.slice(k).col(d).t() * P2 * Lambda2.slice(k).col(d));
    }
    
    //Tau2(D1 + d) = R::rgamma(b2 * (t_alpha + K * (Lambda2.n_rows - 1.0) / 2.0) + 1 - b2, 1.0 / (b2 * (t_beta + temp_beta)));
    for(arma::uword k = 0; k < K; k++){
      //Tau2(D1 + k * D2 + d) = Tau2(D1 + d);
    }
  }
}
void updateThetaLambdaBeta2(arma::mat& Y, arma::mat& X, arma::mat& Z, arma::mat& B){
  arma::uword D1 = X.n_cols;
  arma::uword D2 = Z.n_cols;
  arma::uword p = Beta2.n_rows;
  arma::uword K = Lambda2.n_slices;
  arma::uword n = X.n_rows;
  arma::mat X_eta(n, D1 + K * D2);
  arma::mat C, mychol;
  arma::vec g, w, mu, z, v, result;
  X_eta.cols(0, D1 - 1) = X;
  for(arma::uword k = 0; k < K; k++){
    //X_eta.cols(d*(k+1), d*(k+2)-1) = arma::diagmat(Eta2.col(k)) * Z;
    X_eta.cols(D1 + k * D2, D1 + (k + 1) * D2 - 1) = arma::diagmat(Eta2.col(k)) * Z;
  }
  g = arma::vectorise(B.t() * Y.t() * arma::diagmat(Prec2) * X_eta);
  //C = (arma::kron( X_eta.t() * X_eta, arma::diagmat(Delta)) + arma::kron(arma::diagmat(arma::vectorise(Tau2.t())), P2));
  C = (arma::kron(X_eta.t() * arma::diagmat(Prec2) * X_eta, BtB2) + arma::kron(arma::diagmat(Tau2), P2));
  mychol = arma::chol(C, "lower");
  w = solve(arma::trimatl(mychol), g);
  mu = solve(arma::trimatu(mychol.t()), w);
  z = 1 / sqrt(b2) * arma::randn<arma::vec>(p * (D1 + K * D2));
  v = arma::solve(arma::trimatu(mychol.t()), z);
  result = mu + v;
  
  
  Beta2 = arma::reshape(result.rows(0, p*D1-1), p, D1);
  for(arma::uword k = 0; k < K; k++){
    Lambda2.slice(k) = arma::reshape(result.rows(p*D2*(k+1), p*D2*(k+2) - 1), p, D2);
  }
  
  //return(C);
}
void updateBeta2(arma::mat& Y, arma::mat& X, arma::mat& Z, arma::mat& B){
  /*
  arma::uword D1 = X.n_cols;
  arma::uword p = Beta2.n_rows;
  arma::vec g, w, mu, z, v, result;
  arma::mat PriorPrecision = arma::kron(arma::diagmat(Tau2.head(D1)), P2);
  arma::mat Precision = arma::kron(BtB2, X.t() * arma::diagmat(Prec2) * X) + PriorPrecision;
  arma::mat mychol = arma::chol(Precision, "lower");
  arma::mat Cond = arma::zeros<arma::mat>(Lambda2.n_rows, Y.n_rows);
  for(arma::uword k = 0; k < Lambda2.n_slices; k++){
    Cond = Cond + Lambda2.slice(k) * Z.t() * arma::diagmat(Eta2.col(k));
  }
  g = arma::vectorise((B.t() * Y.t() - BtB2 * Cond) * arma::diagmat(Prec2) * X);
  w = solve(arma::trimatl(mychol), g);
  mu = solve(arma::trimatu(mychol.t()), w);
  z = 1 / sqrt(b2) * arma::randn<arma::vec>(p * D1);
  v = arma::solve(arma::trimatu(mychol.t()), z);
  result = mu + v;
  Beta2 = arma::reshape(result, p, D1);
   */
  arma::uword D1 = X.n_cols;
  arma::uword D2 = Z.n_cols;
  arma::uword p = Beta2.n_rows;
  arma::uword K = Lambda2.n_slices;
  arma::uword n = X.n_rows;
  arma::mat X_eta(n, D1);
  arma::mat C, mychol;
  arma::vec g, w, mu, z, v, result;
  X_eta = X;
  arma::mat Cond = arma::zeros<arma::mat>(Lambda2.n_rows, Y.n_rows);
  for(arma::uword k = 0; k < Lambda2.n_slices; k++){
    Cond = Cond + Lambda2.slice(k) * Z.t() * arma::diagmat(Eta2.col(k));
  }
  g = arma::vectorise((B.t() * Y.t() - BtB2 * Cond) * arma::diagmat(Prec2) * X_eta);
  //C = (arma::kron( X_eta.t() * X_eta, arma::diagmat(Delta)) + arma::kron(arma::diagmat(arma::vectorise(Tau2.t())), P2));
  C = (arma::kron(X_eta.t() * arma::diagmat(Prec2) * X_eta, BtB2) + arma::kron(arma::diagmat(Tau2.head(D1)), P2));
  mychol = arma::chol(C, "lower");
  w = solve(arma::trimatl(mychol), g);
  mu = solve(arma::trimatu(mychol.t()), w);
  z = 1 / sqrt(b2) * arma::randn<arma::vec>(p * (D1));
  v = arma::solve(arma::trimatu(mychol.t()), z);
  result = mu + v;
  
  
  Beta2 = arma::reshape(result.rows(0, p*D1-1), p, D1);
}
void updateLambda2(arma::mat& Y, arma::mat& X, arma::mat& Z, arma::mat B){
  arma::uword D1 = X.n_cols;
  arma::uword D2 = Z.n_cols;
  arma::uword p = Beta2.n_rows;
  arma::uword K = Lambda2.n_slices;
  arma::uword n = X.n_rows;
  arma::mat X_eta(n, K * D2);
  arma::mat C, mychol;
  arma::vec g, w, mu, z, v, result;
  for(arma::uword k = 0; k < K; k++){
    //X_eta.cols(d*(k+1), d*(k+2)-1) = arma::diagmat(Eta2.col(k)) * Z;
    X_eta.cols(k * D2, (k + 1) * D2 - 1) = arma::diagmat(Eta2.col(k)) * Z;
  }
  g = arma::vectorise((B.t() * Y.t() - BtB2 * Beta2 * X.t()) * arma::diagmat(Prec2) * X_eta);
  //C = (arma::kron( X_eta.t() * X_eta, arma::diagmat(Delta)) + arma::kron(arma::diagmat(arma::vectorise(Tau2.t())), P2));
  C = (arma::kron(X_eta.t() * arma::diagmat(Prec2) * X_eta, BtB2) + arma::kron(arma::diagmat(Tau2.tail(K * D2)), P2));
  mychol = arma::chol(C, "lower");
  w = solve(arma::trimatl(mychol), g);
  mu = solve(arma::trimatu(mychol.t()), w);
  z = 1 / sqrt(b2) * arma::randn<arma::vec>(p * (K * D2));
  v = arma::solve(arma::trimatu(mychol.t()), z);
  result = mu + v;
  
  for(arma::uword k = 0; k < K; k++){
    Lambda2.slice(k) = arma::reshape(result.rows(p*D2*(k), p*D2*(k+1) - 1), p, D2);
  }
  
}
/*
void updateDeltaBeta(arma::mat& X, arma::mat& Z){
  arma::uword p = Theta.n_cols;
  arma::uword n = Theta.n_rows;
  arma::uword K = Lambda2.n_slices;
  double a_delta = .1;
  double b_delta = .1;
  arma::vec my_sum = arma::zeros<arma::vec>(p);
  arma::vec ptm;
  for(arma::uword i = 0; i < n; i++){
    ptm = Beta2 * X.row(i).t();
    for(arma::uword k = 0; k < K; k++){
      ptm = ptm + Lambda2.slice(k) * Z.row(i).t() * Eta2(i, k); 
    }
    my_sum = my_sum + arma::square(Theta.row(i).t() - ptm);
  }
  for(arma::uword lp = 0; lp < p; lp++){
    Delta(lp) = R::rgamma(b2 * (a_delta + n/2) + 1 - b2, 1.0 / (b2 * (b_delta + 1.0/2.0 * my_sum(lp))));
  }
}
*/
double updateTausq2(){
  return(R::rgamma(b2 * Phi2.n_elem * Nu2 / 2 + 1 - b2, 1.0 / (b2 * (Nu2 / 2.0 * arma::sum(Phi2)))));
}

void updateNu2(){
  arma::uword n = Phi2.n_elem;
  double nu_alpha = 5;
  double nu_beta = 1;
  double Nu_prop = std::exp(log(Nu2) + R::rnorm(0, 1));

  double log_curr = 0;
  double log_prop = 0;
  for(arma::uword i = 0; i < n; i++){
    log_curr = log_curr + b2 * R::dgamma(Phi2(i), Nu2 / 2, 1.0 / (Nu2 * (Tausq2) / 2), true);
    log_prop = log_prop + b2 * R::dgamma(Phi2(i), Nu_prop / 2, 1.0 / (Nu_prop * (Tausq2) / 2), true);
  }
  log_curr = log_curr + b2 * R::dgamma(Nu2, nu_alpha, 1.0 / nu_beta, true) + std::log(Nu2);
  log_prop = log_curr + b2 * R::dgamma(Nu_prop, nu_alpha, 1.0 / nu_beta, true) + std::log(Nu_prop);
  if(std::log(R::runif(0, 1)) < log_prop - log_curr){
    Nu2 = Nu_prop;
  }
}

// [[Rcpp::export]]
arma::uvec armadillo_modulus3(arma::uvec indicies, arma::uword n){
  return(indicies - n * arma::floor(indicies / n));
}

// [[Rcpp::export]]
void completeY2(arma::mat& Y, arma::uvec missing_sub, arma::uvec missing_time){
  for(arma::uword i = 0; i < missing_sub.n_elem; i++){
    Y(missing_sub(i), missing_time(i)) = fit2(missing_sub(i), missing_time(i)) + R::rnorm(0, std::pow(Prec2(missing_sub(i)), -.5));
  }
}

arma::rowvec update_log_lik2(arma::mat& Y, arma::uword loglik){
  arma::vec SDHolder(Y.n_cols);
  arma::vec current_log_lik;
  
  if(loglik == 1){
    current_log_lik = arma::vec(Y.n_rows * Y.n_cols);
    for(arma::uword i = 0; i < Y.n_rows; i++){
      SDHolder.fill(std::pow(Prec2(i), -.5));
      current_log_lik.subvec(Y.n_cols * i, Y.n_cols * (i + 1) - 1) = arma::log_normpdf(Y.row(i).t(), fit2.row(i).t(), SDHolder);
    }
  }
  if(loglik == 2){
    current_log_lik = arma::vec(Y.n_rows);
    for(arma::uword i = 0; i < Y.n_rows; i++){
      SDHolder.fill(std::pow(Prec2(i), -.5));
      current_log_lik(i) = arma::as_scalar(arma::sum(arma::log_normpdf(Y.row(i).t(), fit2.row(i).t(), SDHolder)));
      
    }
  }
  return(current_log_lik.t());
}
// [[Rcpp::export]]
List run_mcmc_Morris(arma::mat Y, arma::vec Time, arma::mat X, arma::mat Z, arma::mat B, arma::uword K, arma::uword iter, arma::uword burnin, arma::uword nchains, arma::uword thin, arma::uword loglik){//, arma::mat Beta_init, arma::cube Lambda_init, arma::mat Theta_init, arma::mat Eta_init, double Prec_init){
  int p = B.n_cols;
  int N = Y.n_rows;
  int D1 = X.n_cols;
  int D2 = Z.n_cols;
  P2 = getPenalty2(B.n_cols, 2);
  arma::uvec missing = arma::find_nonfinite(Y);
  arma::uvec missing_sub = armadillo_modulus3(missing, Y.n_rows);
  arma::uvec missing_time = arma::floor(missing / Y.n_rows);
  arma::mat Ypred = Y;
  Ypred.elem(missing).fill(0);
  BtB2 = B.t() * B;
  BtY2 = B.t() * Y.t();
  if(loglik == 1){
    log_lik2 = arma::cube(iter, nchains, Y.n_rows * Y.n_cols);
  }
  if(loglik == 2){
    log_lik2 = arma::cube(iter, nchains, Y.n_rows);
  }
  arma::field<arma::cube> LambdaF(nchains, iter);
  arma::field<arma::cube> EtaF(nchains);
  arma::field<arma::mat> TauF(nchains);
  arma::field<arma::cube> BetaF(nchains);
  arma::field<arma::mat> SigmaF(nchains);
  arma::field<arma::mat> PrecF(nchains);
  arma::field<arma::mat> PhiF(nchains);
  arma::field<arma::vec> TausqF(nchains);
  arma::field<arma::vec> NuF(nchains);
  arma::field<arma::vec> AlphaF(nchains);
  for(arma::uword u = 0; u < nchains; u++){
    for(arma::uword i = 0; i < iter; i++){
      LambdaF(u, i) = arma::cube(p, D2, K);
    }
    EtaF(u) = arma::cube(N, K, iter);
    //TauF(u) = arma::cube(K + 1, D, iter);
    TauF(u) = arma::mat(D1 + K * D2, iter);
    BetaF(u) = arma::cube(p, D1, iter);
    SigmaF(u) = arma::mat(K, iter);
    PrecF(u) = arma::mat(N, iter);
    PhiF(u) = arma::mat(N, iter);
    TausqF(u) = arma::vec(iter);
    NuF(u) = arma::vec(iter);
    AlphaF(u) = arma::vec(iter);
  }
  arma::vec Sigma(K);
// #ifdef _OPENMP
//   omp_set_num_threads(12);
//   #pragma omp parallel for shared(LambdaF, BetaF, AlphaF, NuF, TausqF, PhiF, EtaF, PrecF, TauF) schedule(auto)
// #endif  
  for(arma::uword u = 0; u < nchains; u++){
    // Set initial values
    
    
    //arma::cube Lambda2(p, D, K);
    Psi2 = arma::ones<arma::vec>(K);
    Lambda2 = arma::cube(p, D2, K);
    Eta2 = arma::mat(N, K);
    Tau2 = arma::vec(D1 + K * D2);
    Beta2 = arma::mat(p, D1);
    Prec2 = arma::vec(N);
    Phi2 = arma::vec(N);
    arma::cube Lambda_infer(p, D2, K);
    arma::mat Eta_infer(N, K);
    alpha2 = 1;
    Tausq2 = 1;
    Nu2 = 5;
    Lambda2.zeros();
    Eta2.randn();
    Beta2.randn();
    Tau2.ones();
    Prec2.ones();
    Phi2.ones();
    //Lambda2 = Lambda_init;
    //Eta2 = Eta_init;
    arma::mat Y_burnin = Ypred * B * arma::inv(BtB2);
    /*
    Lambda2 = Lambda_init;
    Eta2 = Eta_init;
    Beta2 = Beta_init;
    Theta = Theta_init;
    Prec2 = Prec_init;
    */
    Rcpp::Rcout << "Using newest version" << std::endl;
    Rcpp::Rcout << "Starting burn in..." << std::endl;
    
    for(arma::uword i = 0; i < burnin; i++){
      b2 = .5;
      if(double(i) > double(burnin) / 3.0){
        b2 = .75;
      }
      if(double(i) > double(burnin) * 2.0 / 3.0){
        b2 = 1.0;
      }
      //BtY = Ypred * B;
      fit2 = X * Beta2.t() * B.t();
      for(arma::uword k = 0; k < K; k++){
        fit2 = fit2 + arma::diagmat(Eta2.col(k)) * Z * Lambda2.slice(k).t() * B.t();
      }
      completeY2(Ypred, missing_sub, missing_time);
      //updateAlphaBeta2(Ypred);
      //alpha2 = 1.0;
      //updatePhiBeta2(Ypred);
      //Prec2 = Phi2 * alpha2;
      updatePsi2();
      updatePrec2(Ypred);
      //updateNu2();
      Tausq2 = updateTausq2();
      updateTauBeta2();
      updateLambda2(Ypred, X, Z, B);
      updateBeta2(Ypred, X, Z, B);
      updateEtaBeta2(Ypred, X, Z, B);
      updatePsi2();
    }
    Rcpp::Rcout << "Starting MCMC..." << std::endl;
    double percent = 10;
    for(arma::uword i = 0; i < iter; i++){
      b2 = 1.0;
      
      if(100 * double(i) / iter > percent){
        Rcpp::Rcout << percent << "%" << std::endl;
        percent = 10 + percent;
      }
      for(arma::uword j = 0; j < thin; j++){
        fit2 = X * Beta2.t() * B.t();
        for(arma::uword k = 0; k < K; k++){
          fit2 = fit2 + arma::diagmat(Eta2.col(k)) * Z * Lambda2.slice(k).t() * B.t();
        }
        completeY2(Ypred, missing_sub, missing_time);
        //updateAlphaBeta2(Ypred);
        //updatePhiBeta2(Ypred);
        //Prec2 = Phi2 * alpha2;
        updatePsi2();
        updatePrec2(Ypred);
        //updateNu2();
        Tausq2 = updateTausq2();
        updateTauBeta2();
        updateBeta2(Ypred, X, Z, B);
        updateLambda2(Ypred, X, Z, B);
        updateEtaBeta2(Ypred, X, Z, B);
        if((loglik == 1) || (loglik == 2)){
          log_lik2.tube(i, u) = update_log_lik2(Ypred, loglik);
        }
      }
      
      for(arma::uword k = 0; k < K; k++){
        Lambda_infer.slice(k) = std::pow(Psi2(k), -.5) * Lambda2.slice(k);
        Eta_infer.col(k) = std::pow(Psi2(k), .5) * Eta2.col(k);
      }
      
      LambdaF(u, i) = Lambda_infer;
      EtaF(u).slice(i) = Eta_infer;
      PrecF(u).col(i) = Prec2;
      TauF(u).col(i) = Tau2;
      BetaF(u).slice(i) = Beta2;
      PhiF(u).col(i) = Phi2;
      TausqF(u)(i) = Tausq2;
      NuF(u)(i) = Nu2;
      AlphaF(u)(i) = alpha2;
    }
  }
  Rcpp::Rcout << "All done!";
  List mod = List::create(Named("Lambda", LambdaF), Named("Beta", BetaF),
                          Named("Eta", EtaF),
                          Named("Prec", PrecF), Named("Tau", TauF),
                          Named("Tausq", TausqF), Named("Y", Y),
                          Named("B", B), Named("Phi", PhiF),
                          Named("Nu", NuF), Named("Alpha", AlphaF),
                          Named("log_lik", log_lik2),
                          Named("Time", Time),
                          Named("fit", fit2),
                          Named("X", X),
                          Named("Z", Z));
  return(mod);
}
