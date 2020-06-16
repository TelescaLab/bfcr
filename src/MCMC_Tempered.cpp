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
List MCMC_Tempered(arma::mat Y, arma::mat X, arma::mat B, int K, arma::uword iter, arma::uword nchains, arma::uword thin, double beta, arma::mat Theta_init, arma::cube Lambda_init, arma::mat Eta_init, double Prec_init){
  arma::uword coord = 1;
  double p1 = 0, p2 = 0;
  arma::vec p_vec = arma::zeros<arma::vec>(6);
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
  arma::field<arma::cube> ProjF(nchains);
  arma::field<arma::mat> DeltaF(nchains);
  arma::field<arma::vec> ZF(nchains);
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
    ProjF(u) = arma::cube(N, p, iter);
    DeltaF(u) = arma::mat(p, iter);
    ZF(u) = arma::vec(iter);
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
#ifdef _OPENMP
  omp_set_num_threads(12);
  #pragma omp parallel for shared(LambdaF, ThetaF, EtaF, PrecF, TauF) schedule(auto)
#endif
  for(arma::uword u = 0; u < nchains; u++){
    // Set initial values
    
    
    arma::cube Lambda(p, D, K);
    arma::mat Eta(N, K);
    arma::mat Tau(K + 1, D);
    arma::mat Theta(p, D);
    arma::vec Delta(p);
    arma::mat Proj(N, p);
    double Prec;
    Tau.ones();
    Delta.fill(1);
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
        coord = choose_coordinate(p_vec);
        if(coord == 1){
          p1 = updateThetaLambdaPT(Lambda, Theta, Eta,
                              Delta, Proj, Tau, X,
                              beta);
        }
        if(coord == 2){
          p2 = updateEtaPT(Lambda, Theta, Eta,
                             Delta, Proj, X, beta);
        }
        if(coord == 3){
          updateProj(Lambda, Theta, Eta, Delta, Prec, X, Y, B, Proj);
        }
        if(coord == 4){
          Prec = updatePrecP(Proj, Y, B);
        }
        if(coord == 5){
          updateDelta(Proj, Theta, Lambda, Eta, Delta, X);
        }
        if(coord == 6){
          updateTau(Theta, Lambda, Tau);
        }
        p_vec(0) = p1;
        p_vec(1) = p2;

      }
      
      LambdaF(u, i) = Lambda;
      ThetaF(u).slice(i) = Theta;
      EtaF(u).slice(i) = Eta;
      PrecF(u)(i) = Prec;
      TauF(u).slice(i) = Tau;
      DeltaF(u).col(i) = Delta;
      ProjF(u).slice(i) = Proj;
      ZF(u)(i) = arma::sum(arma::exp(p_vec));
    }
  }
  
  Rcpp::Rcout << "All done!";
  List mod = List::create(Named("Lambda", LambdaF), Named("Theta", ThetaF),
                          Named("Eta", EtaF), Named("Proj", ProjF),
                          Named("Delta", DeltaF),
                          Named("Prec", PrecF), Named("Tau", TauF), Named("ZF", ZF));
  return(mod);
}