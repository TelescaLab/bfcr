#include <RcppArmadillo.h>
#include <omp.h>
#include "updateParam.h"
#include "Utility.h"
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List MCMC(arma::mat Y, arma::mat X, arma::mat B, int K, arma::uword iter, arma::uword nchains, arma::uword thin, arma::mat Theta_init, arma::cube Lambda_init, arma::mat Eta_init, double Prec_init){
  int p = B.n_cols;
  int N = Y.n_rows;
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
  //arma::cube Lambda(p, D, K);
  //arma::mat Eta(N, K);
  //arma::mat Tau(K + 1, D);
//  arma::mat c(K, D);
  //arma::mat Theta(p, D);
  arma::vec Sigma(K);
  //double Prec = 1;
  //Theta = Theta_init;
  //Lambda = Lambda_init;
  //Eta = Eta_init;
  //Y_noise = Y + .05*arma::randn<arma::mat>(N, Y.n_cols);
  Rcpp::Rcout << "Starting MCMC..." << std::endl;
  //omp_set_num_threads(12);
  //#pragma omp parallel for shared(LambdaF, ThetaF, EtaF, PrecF, TauF) schedule(auto)
  for(arma::uword u = 0; u < nchains; u++){
    // Set initial values
    
    
    arma::cube Lambda(p, D, K);
    arma::mat Eta(N, K);
    arma::mat Tau(K + 1, D);
    arma::mat Theta(p, D);
    double Prec;
    Tau.ones();
    /*
    Theta.randn();
    Lambda.randn();
    Eta.randn();
    Prec = 1;
    */
    Lambda = Lambda_init;
    Eta = Eta_init;
    Theta = Theta_init;
    Prec = Prec_init;
    
    
    
    
    //double Prec = PrecEM;
    for(arma::uword i = 0; i < iter; i++){
      if(i % 100 == 0){
        Rcpp::Rcout << i << std::endl;
      }
      //Rcpp::Rcout << i << std::endl;
      for(arma::uword j = 0; j < thin; j++){
        
        //updateEta3(Y, Lambda, Eta, X, B, Prec, Theta);
        Prec = updatePrec(Y, Lambda, Eta, X, B, Theta);
        updateThetaLambda(Y, Lambda, Eta, Tau, X, B, Prec, Theta);
        //updateLambda2(Y, Lambda, Tau, Eta, X, B, Prec, Theta);
        //updateTheta(Y, Lambda, Tau, Eta, X, B, Prec, Theta);
        //updateTheta2(Y, Lambda, Tau, X, B, Prec, Theta);
        updateTau(Theta, Lambda, Tau);
      //Tau.zeros();
        //if(i % 10 == 0){
          updateEta(Y, Lambda, Sigma, Eta, X, B, Prec, Theta);
        //}
      }

      LambdaF(u, i) = Lambda;
      ThetaF(u).slice(i) = Theta;
      EtaF(u).slice(i) = Eta;
      PrecF(u)(i) = Prec;
      TauF(u).slice(i) = Tau;
    }
  }

  Rcpp::Rcout << "All done!";
  List mod = List::create(Named("Lambda", LambdaF), Named("Theta", ThetaF),
                          Named("Eta", EtaF), Named("Sigma", SigmaF),
                          Named("Prec", PrecF), Named("Tau", TauF),
                          Named("c", cF));
  return(mod);
}