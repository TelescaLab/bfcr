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
arma::mat Tau;
double Tausq;
double alpha;
double b;
arma::mat P;

void updateProjBeta2(arma::mat Y, arma::mat X, arma::mat B){
  arma::uword p = B.n_cols;
  for(arma::uword i = 0; i < Y.n_rows; i++){
    arma::vec g = Prec(i) * BtY.row(i).t() + arma::diagmat(Delta) * Beta * X.row(i).t();
    for(arma::uword k = 0; k < Lambda.n_slices; k++){
      g = g + arma::diagmat(Delta) * Lambda.slice(k) * X.row(i).t() * Eta(i, k);
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
  double nu = 1;
  for(arma::uword i = 0; i < Phi.n_elem; i++){
    Phi(i) = R::rgamma(b * (nu + Y.n_cols) / 2 + 1 - b,
        1.0 / (b * (nu * Tausq + arma::dot(Y.row(i) -
          fit.row(i), Y.row(i) - fit.row(i)) * alpha) / 2));
  }
}
void updateEtaPBeta(arma::mat& X){
  arma::uword n = Theta.n_rows;
  arma::mat Xtilde(Theta.n_cols, Lambda.n_slices);
  for(arma::uword i = 0; i < n; i++){
    for(arma::uword k = 0; k < Lambda.n_slices; k++){
      Xtilde.col(k) = Lambda.slice(k) * X.row(i).t();
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
  arma::uword R = Lambda.n_slices;
  arma::uword D = Lambda.n_cols;
  for(arma::uword d = 0; d < D; d++){
    Tau(0, d) = R::rgamma(b * (t_alpha + (Lambda.n_rows - 1.0) / 2.0) + 1 - b, 1.0 / (b * (t_beta + 1.0 / 2.0 * arma::as_scalar(Beta.col(d).t() * P * Beta.col(d)))));
  }
  /*
  for(arma::uword r = 0; r < R; r++){
  for(arma::uword d = 0; d < D; d++){
  Tau(r + 1, d) = R::rgamma(b * (t_alpha + (Lambda.n_rows - 1.0) / 2.0) + 1 - b, 1.0 / (b * (t_beta + 1.0 / 2.0 * arma::as_scalar(Lambda.slice(r).col(d).t() * P * Lambda.slice(r).col(d)))));
  }
  }
  */
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
  }
}
void updateThetaLambdaBeta(arma::mat& X){
  arma::uword d = X.n_cols;
  arma::uword p = Beta.n_rows;
  arma::uword K = Lambda.n_slices;
  arma::uword n = Theta.n_rows;
  arma::mat X_eta(n, (1 + K) * d);
  arma::mat C, mychol;
  arma::vec g, w, mu, z, v, result;
  X_eta.cols(0, d - 1) = X;
  for(arma::uword k = 0; k < K; k++){
    X_eta.cols(d*(k+1), d*(k+2)-1) = arma::diagmat(Eta.col(k)) * X;
  }
  g = arma::vectorise(arma::diagmat(Delta) * Theta.t() * X_eta);
  C = (arma::kron( X_eta.t() * X_eta, arma::diagmat(Delta)) + arma::kron(arma::diagmat(arma::vectorise(Tau.t())), P));
  mychol = arma::chol(C, "lower");
  w = solve(arma::trimatl(mychol), g);
  mu = solve(arma::trimatu(mychol.t()), w);
  z = 1 / sqrt(b) * arma::randn<arma::vec>(p*d*(K+1));
  v = arma::solve(arma::trimatu(mychol.t()), z);
  result = mu + v;
  
  
  Beta = arma::reshape(result.rows(0, p*d-1), p, d);
  for(arma::uword k = 0; k < K; k++){
    Lambda.slice(k) = arma::reshape(result.rows(p*d*(k+1), p*d*(k+2) - 1), p, d);
  }
  
  //return(C);
}
void updateDeltaBeta(arma::mat& X){
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
      ptm = ptm + Lambda.slice(k) * X.row(i).t() * Eta(i, k); 
    }
    my_sum = my_sum + arma::square(Theta.row(i).t() - ptm);
  }
  for(arma::uword lp = 0; lp < p; lp++){
    Delta(lp) = R::rgamma(b * (a_delta + n/2) + 1 - b, 1.0 / (b * (b_delta + 1.0/2.0 * my_sum(lp))));
  }
}
double updateTausq(){
  double nu = 1;
  return(R::rgamma(b * Phi.n_elem * nu / 2 + 1 - b, 1.0 / (b * (nu / 2.0 * arma::sum(Phi)))));
}
// [[Rcpp::export]]
List MCMC(arma::mat Y, arma::mat X, arma::mat B, int K, arma::uword iter, arma::uword burnin, arma::uword nchains, arma::uword thin){//, arma::mat Beta_init, arma::cube Lambda_init, arma::mat Theta_init, arma::mat Eta_init, double Prec_init){
  int p = B.n_cols;
  int N = Y.n_rows;
  int D = X.n_cols;
  P = getPenalty2(B.n_cols, 2);
  BtB = B.t() * B;
  BtY = Y * B;
  arma::field<arma::cube> LambdaF(nchains, iter);
  arma::field<arma::cube> EtaF(nchains);
  arma::field<arma::cube> TauF(nchains);
  arma::field<arma::cube> cF(nchains);
  arma::field<arma::cube> BetaF(nchains);
  arma::field<arma::mat> SigmaF(nchains);
  arma::field<arma::mat> PrecF(nchains);
  arma::field<arma::cube> ThetaF(nchains);
  arma::field<arma::mat> DeltaF(nchains);
  arma::field<arma::mat> PhiF(nchains);
  arma::field<arma::vec> TausqF(nchains);
  for(arma::uword u = 0; u < nchains; u++){
    for(arma::uword i = 0; i < iter; i++){
      LambdaF(u, i) = arma::cube(p, D, K);
    }
    EtaF(u) = arma::cube(N, K, iter);
    TauF(u) = arma::cube(K + 1, D, iter);
    cF(u) = arma::cube(K, D, iter);
    BetaF(u) = arma::cube(p, D, iter);
    SigmaF(u) = arma::mat(K, iter);
    PrecF(u) = arma::mat(N, iter);
    ThetaF(u) = arma::cube(N, p, iter);
    DeltaF(u) = arma::mat(p, iter);
    PhiF(u) = arma::mat(N, iter);
    TausqF(u) = arma::vec(iter);
  }
  arma::vec Sigma(K);
  //omp_set_num_threads(12);
  //#pragma omp parallel for shared(LambdaF, ThetaF, EtaF, PrecF, TauF) schedule(auto)
  for(arma::uword u = 0; u < nchains; u++){
    // Set initial values
    
    
    //arma::cube Lambda(p, D, K);
    Lambda = arma::cube(p, D, K);
    Eta = arma::mat(N, K);
    Tau = arma::mat(K + 1, D);
    Beta = arma::mat(p, D);
    Delta = arma::vec(p);
    Theta = arma::mat(N, p);
    Prec = arma::vec(N);
    Phi = arma::vec(N);
    alpha = 1;
    Tausq = 1;
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
    b = 0.2;
    for(arma::uword i = 0; i < burnin; i++){
      if(double(i) > double(burnin) / 3.0){
        b = 0.5;
      }
      if(double(i) > double(burnin) * 2.0 / 3.0){
        b = 1.0;
      }
      updateProjBeta2(Y, X, B);
      fit = Theta * B.t();
      
      updateAlphaBeta(Y);
      updatePhiBeta(Y);
      Prec = Phi * alpha;
      
      Tausq = updateTausq();
      updateEtaPBeta(X);
      updateTauBeta();
      updateThetaLambdaBeta(X);
      updateDeltaBeta(X);
      
    }
    Rcpp::Rcout << "Starting MCMC..." << std::endl;
    for(arma::uword i = 0; i < iter; i++){
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
        updateProjBeta2(Y, X, B);
        fit = Theta * B.t();
        
        updateAlphaBeta(Y);
        updatePhiBeta(Y);
        Prec = Phi * alpha;
        
        Tausq = updateTausq();
        updateEtaPBeta(X);
        updateTauBeta();
        updateThetaLambdaBeta(X);
        updateDeltaBeta(X);
      }
      
      LambdaF(u, i) = Lambda;
      ThetaF(u).slice(i) = Theta;
      EtaF(u).slice(i) = Eta;
      PrecF(u).col(i) = Prec;
      TauF(u).slice(i) = Tau;
      DeltaF(u).col(i) = Delta;
      BetaF(u).slice(i) = Beta;
      PhiF(u).col(i) = Phi;
      TausqF(u)(i) = Tausq;
    }
  }

  Rcpp::Rcout << "All done!";
  List mod = List::create(Named("Lambda", LambdaF), Named("Beta", BetaF),
                          Named("Eta", EtaF), Named("Theta", ThetaF),
                          Named("Delta", DeltaF),
                          Named("Prec", PrecF), Named("Tau", TauF),
                          Named("Tausq", TausqF), Named("Y", Y),
                          Named("B", B));
  return(mod);
}
