#include <RcppArmadillo.h>
#ifdef _OPENMP
 #include <omp.h>
#endif
#include "updateParam.h"
#include "Utility.h"
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List MCMC_Wrapper(arma::mat Y, arma::mat X, arma::mat B, arma::uword K, arma::uword iter, arma::uword nchains, arma::uword thin, arma::mat Theta_init, arma::cube Lambda_init, arma::mat Eta_init){
  int p = B.n_cols;
  int N = Y.n_rows;
  int D = X.n_cols;
  arma::field<arma::cube> LambdaF(iter);
  arma::cube EtaF(N, K, iter);
  arma::cube TauF(K+1, D, iter);
  arma::cube ThetaF(p, D, iter);
  arma::vec PrecF(iter);
  arma::vec Sigma = arma::ones<arma::vec>(K);
  arma::cube Lambda_temp = Lambda_init;
  arma:: mat Theta_temp = Theta_init;
  arma::mat Eta_temp = Eta_init;
  double Prec_temp = 1;
  arma::mat Tau_temp = arma::ones<arma::mat>(K+1, D);
  for(arma::uword i = 0; i < iter; i++){
    Prec_temp = updatePrec(Y, Lambda_temp, Eta_temp, X, B, Theta_temp);
    updateEta(Y, Lambda_temp, Sigma, Eta_temp, X, B, Prec_temp, Theta_temp);
    updateTau(Theta_temp, Lambda_temp, Tau_temp);
    updateLambda2(Y, Lambda_temp, Tau_temp, Eta_temp, X, B, Prec_temp, Theta_temp);
    updateTheta(Y, Lambda_temp, Tau_temp, Eta_temp, X, B, Prec_temp, Theta_temp);
    LambdaF(i) = Lambda_temp;
    EtaF.slice(i) = Eta_temp;
    ThetaF.slice(i) = Theta_temp;
    TauF.slice(i) = Tau_temp;
    PrecF(i) = Prec_temp;
  }
  List mod = List::create(Named("Prec", PrecF), Named("Lambda", LambdaF), Named("Theta", ThetaF), Named("Tau", TauF), Named("Eta", EtaF));
  return(mod);
}