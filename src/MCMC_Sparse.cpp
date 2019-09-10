#include <RcppArmadillo.h>
#include "updateParam.h"
#include "updateParamSparse.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List MCMC_Sparse(arma::field<arma::vec> Y, arma::mat X, arma::field<arma::mat> B, int K, int iter, int nchains, int thin){
  int p = B(0).n_cols;
  int N = X.n_rows;
  int D = X.n_cols;
  arma::field<arma::cube> LambdaF(nchains, iter);
  arma::field<arma::cube> EtaF(nchains);
  arma::field<arma::cube> TauF(nchains);
  arma::field<arma::cube> cF(nchains);
  arma::field<arma::cube> ThetaF(nchains);
  arma::field<arma::mat> SigmaF(nchains);
  arma::field<arma::vec> PrecF(nchains);
  for(arma::uword u = 0; u < nchains; u++){
    for(arma::uword i = 0; i < iter; i++){
      LambdaF(u, i) = arma::cube(p, D, K);
    }
    EtaF(u) = arma::cube(N, K, iter);
    TauF(u) = arma::cube(K + 1, D, iter);
    cF(u) = arma::cube(K, D, iter);
    ThetaF(u) = arma::cube(p, D, iter);
    SigmaF(u) = arma::mat(K, iter);
    PrecF(u) = arma::vec(iter);
  }
  arma::cube Lambda(p, D, K);
  arma::mat Eta(N, K);
  arma::mat Tau(K + 1, D);
  arma::mat c(K, D);
  arma::mat Theta(p, D);
  arma::vec Sigma(K);
  

  Rcpp::Rcout << "Starting MCMC..." << std::endl;
  for(int u = 0; u < nchains; u++){
    // Set initial values
    Lambda.randn();
    Theta.randn();
    Theta.zeros();
    Sigma.ones();
    Tau.ones();
    c.ones();
    Eta.randn();
    double Prec = 1;
    for(int i = 0; i < iter; i++){
      for(int j = 0; j < thin; j++){
        updateLambdaS(Y, Lambda, Tau, c, Eta, X, B, Prec, Theta);
        //updateLambda2(Ym, L, Tau, c, Eta, X, B(0), Prec, Theta);
        updateThetaS(Y, Lambda, Tau, Eta, X, B, Prec, Theta);
        updateEtaS(Y, Lambda, Sigma, Eta, X, B, Prec, Theta);
        //updateEta(Ym, Lambda, Sigma, Eta, X, B(0), Prec, Theta);
        //updateSigma(Eta, Sigma);
        Prec = updatePrecS(Y, Lambda, Eta, X, B, Theta);
        updateTau(Theta, Lambda, Tau);
        updatec(Lambda, c);
      }
      LambdaF(u, i) = Lambda;
      ThetaF(u).slice(i) = Theta;
      EtaF(u).slice(i) = Eta;
      SigmaF(u).col(i) = Sigma;
      PrecF(u)(i) = Prec;
      TauF(u).slice(i) = Tau;
      cF(u).slice(i) = c;
    }
  }
  
  Rcpp::Rcout << "All done!";
  List mod = List::create(Named("Lambda", LambdaF), Named("Theta", ThetaF),
                          Named("Eta", EtaF), Named("Sigma", SigmaF),
                          Named("Prec", PrecF), Named("Tau", TauF),
                          Named("c", cF));
  return(mod);
}