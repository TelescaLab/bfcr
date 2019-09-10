#include <RcppArmadillo.h>
#include "updateParam.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List MCMC_Impute(arma::field<arma::vec> y, arma::field<arma::vec> observedTimes, arma::vec fullTimes, arma::mat X, arma::mat B, int K, int iter, int nchains, int thin){
  int p = B.n_cols;
  int N = X.n_rows;
  int D = X.n_cols;
  arma::mat ImputedY(X.n_rows, fullTimes.n_elem);
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
  
  arma::field<arma::uvec> observedOrder(N);
  // Initialize complete data
  
  for(arma::uword i = 0; i < X.n_rows; i++){
    ImputedY.row(i) = initializeY(y(i), observedTimes(i), fullTimes);
    observedOrder(i) = getObservedOrder(observedTimes(i), fullTimes);
  }
  
  Rcpp::Rcout << "Starting MCMC..." << std::endl;
  for(int u = 0; u < nchains; u++){
    // Set initial values
    Lambda.randn();
    Theta.randn();
    Sigma.ones();
    Tau.ones();
    c.ones();
    Eta.randn();      
    double Prec = 1;
    for(int i = 0; i < iter; i++){
      for(int j = 0; j < thin; j++){
        updateLambda2(ImputedY, Lambda, Tau, Eta, X, B, Prec, Theta);
        updateTheta(ImputedY, Lambda, Tau, Eta, X, B, Prec, Theta);
        updateEta(ImputedY, Lambda, Sigma, Eta, X, B, Prec, Theta);
        //updateSigma(Eta, Sigma);
        Prec = updatePrec(ImputedY, Lambda, Eta, X, B, Theta);
        updateTau(Theta, Lambda, Tau);
        updatec(Lambda, c);
        //PredictY(ImputedY, X, B, Theta, Eta, Lambda, Prec);
        PredictY2(ImputedY, observedOrder, X, B, Theta, Lambda, Eta, Prec);
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
                          Named("c", cF), Named("ImputedY", ImputedY));
  return(mod);
}