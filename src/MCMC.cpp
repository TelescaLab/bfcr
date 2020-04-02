#include <RcppArmadillo.h>
#include <omp.h>
#include "updateParam.h"
#include "Utility.h"
// [[Rcpp::plugins(openmp)]]
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

arma::cube Lambda;
arma::mat Beta;
arma::mat Theta;
arma::mat Eta;
arma::vec Delta;
arma::vec Prec;
arma::mat BtB;
arma::mat BtY;
arma::mat fit;
arma::vec Phi;
arma::vec Tau;
double Tausq;
double alpha;
double Nu;
double b;
arma::mat P;
arma::cube log_lik;

void updateProjBeta2(arma::mat& Y, arma::mat& X, arma::mat& Z, arma::mat& B){
  arma::uword p = B.n_cols;
  for(arma::uword i = 0; i < Y.n_rows; i++){
    arma::vec g = Prec(i) * BtY.row(i).t() + arma::diagmat(Delta) * Beta * X.row(i).t();
    for(arma::uword k = 0; k < Lambda.n_slices; k++){
      g = g + arma::diagmat(Delta) * Lambda.slice(k) * Z.row(i).t() * Eta(i, k);
    }
    arma::mat C = (Prec(i) * BtB + arma::diagmat(Delta));
    arma::mat mychol = arma::chol(C, "lower");
    arma::vec w = solve(arma::trimatl(mychol), g);
    arma::vec mu = solve(arma::trimatu(mychol.t()), w);
    arma::vec z = 1 / sqrt(b) * arma::randn<arma::vec>(p);
    arma::vec v = arma::solve(arma::trimatu(mychol.t()), z);
    Theta.row(i) = (mu + v).t();
  }
}
void updateAlphaBeta(arma::mat& Y){
  // Phi has precision parameterization
  // Output is on the precision scale
  arma::vec norm_vec(Y.n_rows);
  for(arma::uword i = 0; i < Y.n_rows; i++){
    norm_vec(i) = arma::sum(arma::square(Y.row(i) - fit.row(i))) * Phi(i);
  }
  alpha = R::rgamma(b * (Y.n_elem)/2 + 1 - b, 1.0 / (b * (1.0/2.0 * arma::sum(norm_vec))));
}
void updatePhiBeta(arma::mat& Y){
  // alpha is on precision scale
  // output on precision scale
  for(arma::uword i = 0; i < Phi.n_elem; i++){
    Phi(i) = R::rgamma(b * (Nu + Y.n_cols) / 2 + 1 - b,
        1.0 / (b * (Nu * Tausq + arma::dot(Y.row(i) -
          fit.row(i), Y.row(i) - fit.row(i)) * alpha) / 2));
  }
}
void updateEtaPBeta(arma::mat& X, arma::mat& Z){
  arma::uword n = Theta.n_rows;
  arma::mat Xtilde(Theta.n_cols, Lambda.n_slices);
  for(arma::uword i = 0; i < n; i++){
    for(arma::uword k = 0; k < Lambda.n_slices; k++){
      Xtilde.col(k) = Lambda.slice(k) * Z.row(i).t();
    }
    arma::vec g = Xtilde.t() * arma::diagmat(Delta) * (Theta.row(i).t() - Beta * X.row(i).t());
    arma::mat C = (Xtilde.t() * arma::diagmat(Delta) * Xtilde + arma::eye(Lambda.n_slices, Lambda.n_slices));
    arma::mat mychol = arma::chol(C, "lower");
    arma::vec w = solve(arma::trimatl(mychol), g);
    arma::vec mu = solve(arma::trimatu(mychol.t()), w);
    arma::vec z = 1/sqrt(b) * arma::randn<arma::vec>(Lambda.n_slices);
    arma::vec v = arma::solve(arma::trimatu(mychol.t()), z);
    Eta.row(i) = (mu + v).t();
  }
}
void updateTauBeta(){
  double t_alpha = 1;
  double t_beta =  .0005;
  //double t_beta = 1;
  arma::uword K = Lambda.n_slices;
  arma::uword D1 = Beta.n_cols;
  arma::uword D2 = Lambda.n_cols;
  
  /*
  for(arma::uword r = 0; r < R; r++){
  for(arma::uword d = 0; d < D; d++){
  Tau(r + 1, d) = R::rgamma(b * (t_alpha + (Lambda.n_rows - 1.0) / 2.0) + 1 - b, 1.0 / (b * (t_beta + 1.0 / 2.0 * arma::as_scalar(Lambda.slice(r).col(d).t() * P * Lambda.slice(r).col(d)))));
  }
  }
  */
  /*
  for(arma::uword d = 0; d < D; d++){
    Tau(0, d) = R::rgamma(b * (t_alpha + (Lambda.n_rows - 1.0) / 2.0) + 1 - b, 1.0 / (b * (t_beta + 1.0 / 2.0 * arma::as_scalar(Beta.col(d).t() * P * Beta.col(d)))));
  }
  double temp_beta;
  for(arma::uword d = 0; d < D; d++){
    temp_beta = 0;
    for(arma::uword r = 1; r < R; r++){
      temp_beta = temp_beta + 1.0 / 2.0 * arma::as_scalar(Lambda.slice(r).col(d).t() * P * Lambda.slice(r).col(d));
    }
    Tau(1, d) = R::rgamma(b * (t_alpha + R * (Lambda.n_rows - 1.0) / 2.0) + 1 - b, 1.0 / (b * (t_beta + temp_beta)));
    for(arma::uword r = 0; r < R; r++){
      Tau(r + 1, d) = Tau(1, d);
    }
  }*/
  for(arma::uword d = 0; d < D1; d++){
    Tau(d) = R::rgamma(b * (t_alpha + (Lambda.n_rows - 1.0) / 2.0) + 1 - b, 1.0 / (b * (t_beta + 1.0 / 2.0 * arma::as_scalar(Beta.col(d).t() * P * Beta.col(d)))));
  }
  double temp_beta;
  for(arma::uword d = 0; d < D2; d++){
    temp_beta = 0;
    for(arma::uword k = 0; k < K; k++){
      temp_beta = temp_beta + 1.0 / 2.0 * arma::as_scalar(Lambda.slice(k).col(d).t() * P * Lambda.slice(k).col(d));
    }
    
    Tau(D1 + d) = R::rgamma(b * (t_alpha + K * (Lambda.n_rows - 1.0) / 2.0) + 1 - b, 1.0 / (b * (t_beta + temp_beta)));
    for(arma::uword k = 0; k < K; k++){
      Tau(D1 + k * D2 + d) = Tau(D1 + d);
    }
  }
}
void updateThetaLambdaBeta(arma::mat& X, arma::mat& Z){
  arma::uword D1 = X.n_cols;
  arma::uword D2 = Z.n_cols;
  arma::uword p = Beta.n_rows;
  arma::uword K = Lambda.n_slices;
  arma::uword n = Theta.n_rows;
  arma::mat X_eta(n, D1 + K * D2);
  arma::mat C, mychol;
  arma::vec g, w, mu, z, v, result;
  X_eta.cols(0, D1 - 1) = X;
  for(arma::uword k = 0; k < K; k++){
    //X_eta.cols(d*(k+1), d*(k+2)-1) = arma::diagmat(Eta.col(k)) * Z;
    X_eta.cols(D1 + k * D2, D1 + (k + 1) * D2 - 1) = arma::diagmat(Eta.col(k)) * Z;
  }
  g = arma::vectorise(arma::diagmat(Delta) * Theta.t() * X_eta);
  //C = (arma::kron( X_eta.t() * X_eta, arma::diagmat(Delta)) + arma::kron(arma::diagmat(arma::vectorise(Tau.t())), P));
  C = (arma::kron( X_eta.t() * X_eta, arma::diagmat(Delta)) + arma::kron(arma::diagmat(Tau), P));
  mychol = arma::chol(C, "lower");
  w = solve(arma::trimatl(mychol), g);
  mu = solve(arma::trimatu(mychol.t()), w);
  z = 1 / sqrt(b) * arma::randn<arma::vec>(p * (D1 + K * D2));
  v = arma::solve(arma::trimatu(mychol.t()), z);
  result = mu + v;
  
  
  Beta = arma::reshape(result.rows(0, p*D1-1), p, D1);
  for(arma::uword k = 0; k < K; k++){
    Lambda.slice(k) = arma::reshape(result.rows(p*D2*(k+1), p*D2*(k+2) - 1), p, D2);
  }
  
  //return(C);
}

void updateBeta(arma::mat& X, arma::mat& Z){
  arma::uword D1 = X.n_cols;
  arma::uword D2 = Z.n_cols;
  arma::uword p = Beta.n_rows;
  arma::uword K = Lambda.n_slices;
  arma::uword n = Theta.n_rows;
  arma::mat X_eta(n, D1);
  arma::mat C, mychol;
  arma::vec g, w, mu, z, v, result;
  X_eta = X;
  arma::mat Cond = arma::zeros<arma::mat>(p, n);
  for(arma::uword k = 0; k < K; k++){
    Cond = Cond + Lambda.slice(k) * Z.t() * arma::diagmat(Eta.col(k));
  }
  g = arma::vectorise(arma::diagmat(Delta) * (Theta.t() - Cond) * X_eta);
  //C = (arma::kron( X_eta.t() * X_eta, arma::diagmat(Delta)) + arma::kron(arma::diagmat(arma::vectorise(Tau.t())), P));
  C = (arma::kron( X_eta.t() * X_eta, arma::diagmat(Delta)) + arma::kron(arma::diagmat(Tau.head(D1)), P));
  mychol = arma::chol(C, "lower");
  w = solve(arma::trimatl(mychol), g);
  mu = solve(arma::trimatu(mychol.t()), w);
  z = 1 / sqrt(b) * arma::randn<arma::vec>(p * (D1));
  v = arma::solve(arma::trimatu(mychol.t()), z);
  result = mu + v;
  
  Beta = arma::reshape(result.rows(0, p*D1-1), p, D1);
  
}

void updateLambda(arma::mat& X, arma::mat& Z){
  arma::uword D1 = X.n_cols;
  arma::uword D2 = Z.n_cols;
  arma::uword p = Beta.n_rows;
  arma::uword K = Lambda.n_slices;
  arma::uword n = Theta.n_rows;
  arma::mat X_eta(n, K * D2);
  arma::mat C, mychol;
  arma::vec g, w, mu, z, v, result;
  for(arma::uword k = 0; k < K; k++){
    //X_eta.cols(d*(k+1), d*(k+2)-1) = arma::diagmat(Eta.col(k)) * Z;
    X_eta.cols(k * D2, (k + 1) * D2 - 1) = arma::diagmat(Eta.col(k)) * Z;
  }
  g = arma::vectorise(arma::diagmat(Delta) * (Theta.t() - Beta * X.t()) * X_eta);
  //C = (arma::kron( X_eta.t() * X_eta, arma::diagmat(Delta)) + arma::kron(arma::diagmat(arma::vectorise(Tau.t())), P));
  C = (arma::kron(X_eta.t() * X_eta, arma::diagmat(Delta)) + arma::kron(arma::diagmat(Tau.tail(D2 * K)), P));
  mychol = arma::chol(C, "lower");
  w = solve(arma::trimatl(mychol), g);
  mu = solve(arma::trimatu(mychol.t()), w);
  z = 1 / sqrt(b) * arma::randn<arma::vec>(p * (K * D2));
  v = arma::solve(arma::trimatu(mychol.t()), z);
  result = mu + v;
  
  
  for(arma::uword k = 0; k < K; k++){
    Lambda.slice(k) = arma::reshape(result.rows(p*D2*(k), p*D2*(k+1) - 1), p, D2);
  }
}
void updateDeltaBeta(arma::mat& X, arma::mat& Z){
  arma::uword p = Theta.n_cols;
  arma::uword n = Theta.n_rows;
  arma::uword K = Lambda.n_slices;
  double a_delta = .1;
  double b_delta = .1;
  arma::vec my_sum = arma::zeros<arma::vec>(p);
  arma::vec ptm;
  for(arma::uword i = 0; i < n; i++){
    ptm = Beta * X.row(i).t();
    for(arma::uword k = 0; k < K; k++){
      ptm = ptm + Lambda.slice(k) * Z.row(i).t() * Eta(i, k); 
    }
    my_sum = my_sum + arma::square(Theta.row(i).t() - ptm);
  }
  for(arma::uword lp = 0; lp < p; lp++){
    Delta(lp) = R::rgamma(b * (a_delta + n/2) + 1 - b, 1.0 / (b * (b_delta + 1.0/2.0 * my_sum(lp))));
  }
}
double updateTausq(){
  return(R::rgamma(b * Phi.n_elem * Nu / 2 + 1 - b, 1.0 / (b * (Nu / 2.0 * arma::sum(Phi)))));
}

void updateNu(){
  arma::uword n = Phi.n_elem;
  double nu_alpha = 10;
  double nu_beta = 1;
  double Nu_prop = std::exp(log(Nu) + R::rnorm(0, 1));
  //Rcout << "Nu current: " << Nu << std::endl << "Nu proposal: " << Nu_prop << std::endl;
  
  double log_curr = 0;
  double log_prop = 0;
  for(arma::uword i = 0; i < n; i++){
    log_curr = log_curr + b * R::dgamma(Phi(i), Nu / 2, 1.0 / (Nu * (Tausq) / 2), true);
    log_prop = log_prop + b * R::dgamma(Phi(i), Nu_prop / 2, 1.0 / (Nu_prop * (Tausq) / 2), true);
  }
  log_curr = log_curr + b * R::dgamma(Nu, nu_alpha, 1.0 / nu_beta, true) + std::log(Nu);
  log_prop = log_curr + b * R::dgamma(Nu_prop, nu_alpha, 1.0 / nu_beta, true) + std::log(Nu_prop);
  //Rcout << log_curr << std::endl << log_prop << std::endl;
  if(std::log(R::runif(0, 1)) < log_prop - log_curr){
    Nu = Nu_prop;
  }
}

// [[Rcpp::export]]
arma::uvec armadillo_modulus(arma::uvec indicies, arma::uword n){
  return(indicies - n * arma::floor(indicies / n));
}

// [[Rcpp::export]]
void completeY(arma::mat& Y, arma::uvec missing_sub, arma::uvec missing_time){
  for(arma::uword i = 0; i < missing_sub.n_elem; i++){
    Y(missing_sub(i), missing_time(i)) = fit(missing_sub(i), missing_time(i)) + R::rnorm(0, std::pow(Prec(missing_sub(i)), -.5));
  }
}

arma::rowvec update_log_lik(arma::mat& Y, arma::uword loglik){
  arma::vec SDHolder(Y.n_cols);
  arma::vec current_log_lik;

  if(loglik == 1){
    current_log_lik = arma::vec(Y.n_rows * Y.n_cols);
    for(arma::uword i = 0; i < Y.n_rows; i++){
      SDHolder.fill(std::pow(Prec(i), -.5));
      current_log_lik.subvec(Y.n_cols * i, Y.n_cols * (i + 1) - 1) = arma::log_normpdf(Y.row(i).t(), fit.row(i).t(), SDHolder);
    }
  }
  if(loglik == 2){
    current_log_lik = arma::vec(Y.n_rows);
    for(arma::uword i = 0; i < Y.n_rows; i++){
      SDHolder.fill(std::pow(Prec(i), -.5));
      current_log_lik(i) = arma::as_scalar(arma::sum(arma::log_normpdf(Y.row(i).t(), fit.row(i).t(), SDHolder)));
      
    }
  }
  return(current_log_lik.t());
}
// [[Rcpp::export]]
List run_mcmc(arma::mat Y, arma::vec Time, arma::mat X, arma::mat Z, arma::mat B, int K, arma::uword iter, arma::uword burnin, arma::uword nchains, arma::uword thin, arma::uword loglik){//, arma::mat Beta_init, arma::cube Lambda_init, arma::mat Theta_init, arma::mat Eta_init, double Prec_init){
  int p = B.n_cols;
  int N = Y.n_rows;
  int D1 = X.n_cols;
  int D2 = Z.n_cols;
  P = getPenalty2(B.n_cols, 2);
  arma::uvec missing = arma::find_nonfinite(Y);
  arma::uvec missing_sub = armadillo_modulus(missing, Y.n_rows);
  arma::uvec missing_time = arma::floor(missing / Y.n_rows);
  arma::mat Ypred = Y;
  Ypred.elem(missing).fill(0);
  BtB = B.t() * B;
  if(loglik == 1){
    log_lik = arma::cube(iter, nchains, Y.n_rows * Y.n_cols);
  }
  if(loglik == 2){
    log_lik = arma::cube(iter, nchains, Y.n_rows);
  }
  arma::field<arma::cube> LambdaF(nchains, iter);
  arma::field<arma::cube> EtaF(nchains);
  arma::field<arma::mat> TauF(nchains);
  arma::field<arma::cube> BetaF(nchains);
  arma::field<arma::mat> SigmaF(nchains);
  arma::field<arma::mat> PrecF(nchains);
  arma::field<arma::cube> ThetaF(nchains);
  arma::field<arma::mat> DeltaF(nchains);
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
    ThetaF(u) = arma::cube(N, p, iter);
    DeltaF(u) = arma::mat(p, iter);
    PhiF(u) = arma::mat(N, iter);
    TausqF(u) = arma::vec(iter);
    NuF(u) = arma::vec(iter);
    AlphaF(u) = arma::vec(iter);
  }
  arma::vec Sigma(K);
  //omp_set_num_threads(12);
  //#pragma omp parallel for shared(LambdaF, ThetaF, EtaF, PrecF, TauF) schedule(auto)
  for(arma::uword u = 0; u < nchains; u++){
    // Set initial values
    
    
    //arma::cube Lambda(p, D, K);
    Lambda = arma::cube(p, D2, K);
    Eta = arma::mat(N, K);
    Tau = arma::vec(D1 + K * D2);
    Beta = arma::mat(p, D1);
    Delta = arma::vec(p);
    Theta = arma::mat(N, p);
    Prec = arma::vec(N);
    Phi = arma::vec(N);
    
    alpha = 1;
    Tausq = 1;
    Nu = 5;
    Theta.randn();
    Lambda.randn();
    Eta.randn();
    Beta.randn();
    Tau.ones();
    Delta.ones();
    Prec.ones();
    Phi.ones();
    /*
    Lambda = Lambda_init;
    Eta = Eta_init;
    Beta = Beta_init;
    Theta = Theta_init;
    Prec = Prec_init;
    */
    Rcpp::Rcout << "Starting burn in..." << std::endl;
    for(arma::uword i = 0; i < burnin; i++){
      b = .5;
      if(double(i) > double(burnin) / 3.0){
        b = .75;
      }
      if(double(i) > double(burnin) * 2.0 / 3.0){
        b = 1.0;
      }
      BtY = Ypred * B;
      updateProjBeta2(Ypred, X, Z, B);
      fit = Theta * B.t();
      completeY(Ypred, missing_sub, missing_time);
      updateAlphaBeta(Ypred);
      updatePhiBeta(Ypred);
      updateNu();
      Prec = Phi * alpha;
      Tausq = updateTausq();
      updateEtaPBeta(X, Z);
      updateTauBeta();
      //updateThetaLambdaBeta(X, Z);
      updateBeta(X, Z);
      updateLambda(X, Z);
      updateDeltaBeta(X, Z);

    }
    Rcpp::Rcout << "Starting MCMC..." << std::endl;
    for(arma::uword i = 0; i < iter; i++){
      b = 1.0;
      if(i % 100 == 0){
        Rcpp::Rcout << i << std::endl;
      }
      for(arma::uword j = 0; j < thin; j++){
        //updateProjBeta(Lambda, Beta, Eta, Delta, Prec, X, Y, B, Theta, b);
        //updateEtaP(Lambda, Beta, Eta, Delta, Theta, X);
        //updateTau(Beta, Lambda, Tau);
        //updateThetaLambdaP(Lambda, Beta, Eta, Delta, Theta, Tau, X);
        //Prec = updatePrecP(Theta, Y, B);
        //updateDelta(Theta, Beta, Lambda, Eta, Delta, X);
        
        //updateProjBeta(Lambda, Beta, Eta, Delta, Prec, X, Y, B, Theta, b);
        BtY = Ypred * B;
        updateProjBeta2(Ypred, X, Z, B);
        fit = Theta * B.t();
        
        completeY(Ypred, missing_sub, missing_time);
        
        updateAlphaBeta(Ypred);
        updatePhiBeta(Ypred);
        Prec = Phi * alpha;
        updateNu();
        Tausq = updateTausq();
        updateEtaPBeta(X, Z);
        updateTauBeta();
        updateBeta(X, Z);
        updateLambda(X, Z);
        //updateThetaLambdaBeta(X, Z);
        updateDeltaBeta(X, Z);
        if((loglik == 1) || (loglik == 2)){
          log_lik.tube(i, u) = update_log_lik(Ypred, loglik);
        }
      }
      
      LambdaF(u, i) = Lambda;
      ThetaF(u).slice(i) = Theta;
      EtaF(u).slice(i) = Eta;
      PrecF(u).col(i) = Prec;
      TauF(u).col(i) = Tau;
      DeltaF(u).col(i) = Delta;
      BetaF(u).slice(i) = Beta;
      PhiF(u).col(i) = Phi;
      TausqF(u)(i) = Tausq;
      NuF(u)(i) = Nu;
      AlphaF(u)(i) = alpha;
    }
  }
  Rcpp::Rcout << "All done!";
  List mod = List::create(Named("Lambda", LambdaF), Named("Beta", BetaF),
                          Named("Eta", EtaF), Named("Theta", ThetaF),
                          Named("Delta", DeltaF),
                          Named("Prec", PrecF), Named("Tau", TauF),
                          Named("Tausq", TausqF), Named("Y", Y),
                          Named("B", B), Named("Phi", PhiF),
                          Named("Nu", NuF), Named("Alpha", AlphaF),
                          Named("log_lik", log_lik),
                          Named("Time", Time));
  return(mod);
}
