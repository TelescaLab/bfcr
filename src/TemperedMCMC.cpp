#include <RcppArmadillo.h>
#include <omp.h>
#include "updateParam.h"
#include "Utility.h"
// [[Rcpp::plugins(openmp)]]

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List TemperedMCMC(arma::mat Y, arma::mat X, arma::mat B, int K, arma::uword iter,
                  arma::uword thin, arma::mat Theta_init,
                  arma::cube Lambda_init, arma::mat Eta_init, double Prec_init,
                  arma::vec beta){
  int p = B.n_cols;
  int N = Y.n_rows;
  int D = X.n_cols;
  int u1;
  int u2;
  double lik1, lik2;
  arma::uword nchains = beta.n_elem;
  arma::field<arma::cube> LambdaF(nchains, iter+1);
  arma::field<arma::cube> EtaF(nchains);
  arma::field<arma::cube> TauF(nchains);
  arma::field<arma::cube> cF(nchains);
  arma::field<arma::cube> ThetaF(nchains);
  arma::field<arma::mat> SigmaF(nchains);
  arma::field<arma::vec> PrecF(nchains);
  arma::field<arma::cube> ProjF(nchains);
  arma::field<arma::mat> DeltaF(nchains);
  for(arma::uword u = 0; u < nchains; u++){
    for(arma::uword i = 0; i < iter; i++){
      LambdaF(u, i) = arma::cube(p, D, K);
    }
    EtaF(u) = arma::cube(N, K, iter+1);
    TauF(u) = arma::cube(K + 1, D, iter+1);
    cF(u) = arma::cube(K, D, iter+1);
    ThetaF(u) = arma::cube(p, D, iter+1);
    SigmaF(u) = arma::mat(K, iter+1);
    PrecF(u) = arma::vec(iter+1);
    ProjF(u) = arma::cube(N, p, iter+1);
    DeltaF(u) = arma::mat(p, iter+1);
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
  arma::vec DeltaTemp1, DeltaTemp2;
  arma::mat ThetaTemp1, ThetaTemp2, EtaTemp1, EtaTemp2, TauTemp1, TauTemp2, ProjTemp1, ProjTemp2;
  arma::cube LambdaTemp1, LambdaTemp2;
  double PrecTemp1, PrecTemp2;
  Rcpp::NumericVector switches (2);
  arma::cube Lambda(p, D, K);
  arma::mat Eta(N, K);
  arma::mat Tau(K + 1, D);
  arma::mat Theta(p, D);
  arma::vec Delta(p);
  arma::mat Proj(N, p);
  double Prec;
  /*
  Theta.randn();
  Lambda.randn();
  Eta.randn();
  Prec = 1;
  */
  Rcpp::NumericVector probs (nchains, 1);
  for(arma::uword u = 0; u < nchains; u++){
    LambdaF(u, 0) = Lambda_init;
    EtaF(u).slice(0) = Eta_init;
    ThetaF(u).slice(0) = Theta_init;
    PrecF(u)(0) = Prec_init;
    TauF(u).slice(0).ones();
    DeltaF(u).col(0).ones();
    ProjF(u).slice(0).randn();
  }
  
  for(arma::uword i = 0; i < iter; i++){
    for(arma::uword u = 0; u < nchains; u++){
      Lambda = LambdaF(u, i);
      Theta = ThetaF(u).slice(i);
      Eta = EtaF(u).slice(i);
      Prec = PrecF(u)(i);
      Tau = TauF(u).slice(i);
      Delta = DeltaF(u).col(i);
      Proj = ProjF(u).slice(i);
      updateProjBeta(Lambda, Theta, Eta, Delta, Prec, X, Y, B, Proj, beta(u));
      Prec = updatePrecPBeta(Proj, Y, B, beta(u));
      //updateDeltaBeta(Proj, Theta, Lambda, Eta, Delta, X, beta(u));
      //updateEtaPBeta(Lambda, Theta, Eta, Delta, Proj, X, beta(u));
      //updateThetaLambdaBeta(Lambda, Theta, Eta, Delta, Proj, Tau, X, beta(u));
      updateTauBeta(Theta, Lambda, Tau, beta(u));
      LambdaF(u, i+1) = Lambda;
      ThetaF(u).slice(i+1) = Theta;
      EtaF(u).slice(i+1) = Eta;
      PrecF(u)(i+1) = Prec;
      TauF(u).slice(i+1) = Tau;
      DeltaF(u).col(i+1) = Delta;
      ProjF(u).slice(i+1) = Proj;
    }
    if(i % 10 == 0){
      switches = Rcpp::sample(int(nchains), 1, false);
      u1 = switches(0) - 1;
      if(u1 == 0){
        u2 = 1;
      } else if(u1 == int(nchains) - 1){
        u2 = nchains - 2;
      } else{
        u2 = Rcpp::sample(Rcpp::NumericVector::create(u1-1, u1+1), 1, false)(0);
      }
      DeltaTemp1 = DeltaF(u1).col(i+1);
      DeltaTemp2 = DeltaF(u2).col(i+1);
      ThetaTemp1 = ThetaF(u1).slice(i+1);
      ThetaTemp2 = ThetaF(u2).slice(i+1);
      LambdaTemp1 = LambdaF(u1, i+1);
      LambdaTemp2 = LambdaF(u2, i+1);
      PrecTemp1 = PrecF(u1)(i+1);
      PrecTemp2 = PrecF(u2)(i+1);
      ProjTemp1 = ProjF(u1).slice(i+1);
      ProjTemp2 = ProjF(u2).slice(i+1);
      TauTemp1 = TauF(u1).slice(i+1);
      TauTemp2 = TauF(u2).slice(i+1);
      EtaTemp1 = EtaF(u1).slice(i+1);
      EtaTemp2 = EtaF(u2).slice(i+1);
      lik1 = cpploglik_bayes(ThetaTemp1, LambdaTemp1, PrecTemp1, DeltaTemp1,
                             X, B, Y, 12);
      lik2 = cpploglik_bayes(ThetaTemp2, LambdaTemp2, PrecTemp2, DeltaTemp2,
                             X, B, Y, 12);
      double A = exp((lik1 - lik2) * (beta(u1) - beta(u2)));
      if(R::runif(0.0, 1.0) < A){
        ThetaF(u1).slice(i+1) = ThetaTemp2;
        ThetaF(u2).slice(i+1) = ThetaTemp1;
        LambdaF(u1, i+1) = LambdaTemp2;
        LambdaF(u2, i+1) = LambdaTemp1;
        PrecF(u1)(i+1) = PrecTemp2;
        PrecF(u2)(i+1) = PrecTemp1;
        DeltaF(u1).col(i+1) = DeltaTemp2;
        DeltaF(u2).col(i+1) = DeltaTemp1;
        ProjF(u1).slice(i+1) = ProjTemp2;
        ProjF(u2).slice(i+1) = ProjTemp1;
        TauF(u1).slice(i+1) = TauTemp2;
        TauF(u2).slice(i+1) = TauTemp1;
        EtaF(u1).slice(i+1) = EtaTemp2;
        EtaF(u2).slice(i+1) = EtaTemp1;
      }
    
    
    }
  }

  List mod = List::create(Named("Lambda", LambdaF), Named("Theta", ThetaF),
                          Named("Eta", EtaF), Named("Proj", ProjF),
                          Named("Delta", DeltaF),
                          Named("Prec", PrecF), Named("Tau", TauF));
  return(mod);
}