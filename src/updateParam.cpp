#include <RcppArmadillo.h>
#include "Utility.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
void updateLambda(arma::mat& Y, arma::cube& Lambda, arma::vec& r, arma::mat& Gamma, arma::mat& X, arma::mat& B, double prec){
  arma::mat G;
  arma::mat condMean = arma::zeros<arma::mat>(Y.n_rows, Y.n_cols);
  arma::uword Q = Lambda.n_slices;
  arma::uword counter = 0;
  arma::mat mu;
  arma::mat Z;
  arma::mat D;
  arma::mat Pb;
  
  //arma::mat Dchol;
  //arma::mat Pbchol;
  for(arma::uword outer = 0; outer < Q; outer++){
    for(arma::uword q = 0; q < Q; q++){
      if(q != counter){
        condMean = condMean + arma::diagmat(Gamma.col(q)) * X * arma::trans(Lambda.slice(q)) * B.t();
      }
    }
    G = arma::diagmat(Gamma.col(outer));
    mu = arma::trans(prec * arma::solve(X.t() * arma::square(arma::diagmat(G)) * X, X.t()) * 
      arma::diagmat(G) * (Y - condMean) * B * arma::inv(prec * B.t() * B + r(outer) * getPenalty2(B.n_cols, 2)));
    Z = arma::randn<arma::mat>(B.n_cols, X.n_cols);
    D = (arma::chol(X.t() * arma::square(arma::diagmat(G)) * X));
    Pb = (arma::chol(prec * B.t() * B + r(outer) * getPenalty2(B.n_cols, 2), "lower"));
    //Dchol = arma::chol(D);
    //Pbchol = arma::chol(Pb, "lower");
    Lambda.slice(outer) = mu + arma::solve(arma::trimatl(Pb), arma::solve(arma::trimatl(D.t()), Z.t()).t());
    condMean.zeros();
    counter++;
  }
}

// [[Rcpp::export]]
void updateLambda2(arma::mat& Y, arma::cube& Lambda, arma::mat& Tau, arma::mat& Gamma, arma::mat& X, arma::mat& B, double prec, arma::mat& Theta){
  //arma::mat P = getPenalty(Lambda.n_rows);
  arma::mat P = getPenalty2(Lambda.n_rows, 2);
  arma::mat PriorPrecision;
  arma::mat Precision;
  arma::mat mychol;
  arma::vec mymean;
  arma::vec Ltemp;
  arma::uword Q = Lambda.n_slices;
  arma::uword counter = 0;
  arma::mat condMean = arma::zeros<arma::mat>(Y.n_rows, Y.n_cols);
  for(arma::uword i = 0; i < Q; i++){
    for(arma::uword q = 0; q < Q; q++){
      if(q != counter){
        condMean = condMean + arma::diagmat(Gamma.col(q)) * X * arma::trans(Lambda.slice(q)) * B.t();
      }
    }
    //PriorPrecision = arma::kron(P, arma::diagmat(Tau.row(i+1))) + arma::kron(arma::eye<arma::mat>(Lambda.n_rows, Lambda.n_rows), arma::diagmat(c.row(i)));
    PriorPrecision = arma::kron(P, arma::diagmat(Tau.row(i+1)));
    Precision = prec * arma::kron(B.t() * B, X.t() * arma::diagmat(arma::square(Gamma.col(i))) * X) + PriorPrecision;
    mychol = arma::chol(Precision, "lower");
    mymean = arma::solve(arma::trimatu(mychol.t()), arma::solve(arma::trimatl(mychol), prec * arma::vectorise(X.t() * arma::diagmat(Gamma.col(i)) * (Y - X * Theta.t() * B.t() - condMean) * B)));
    
    Ltemp = mymean + arma::solve(arma::trimatl(mychol), arma::randn<arma::vec>(Lambda.n_rows * Lambda.n_cols));
    Lambda.slice(i) = arma::trans(arma::reshape(Ltemp,Lambda.n_cols,Lambda.n_rows));
    //Lambda.slice(i) = arma::orth(Lambda.slice(i));
    //Ltemp = mymean;
    //Lambda.slice(i) = arma::trans(arma::reshape(Ltemp,Lambda.n_cols,Lambda.n_rows));
    counter++;
    condMean.zeros();
   }

}

// [[Rcpp::export]]
void updateTheta(arma::mat& Y, arma::cube& Lambda, arma::mat& Tau, arma::mat& Gamma, arma::mat& X, arma::mat& B, double prec, arma::mat& Theta){
  arma::mat P = getPenalty2(Theta.n_rows, 2);
  arma::mat PriorPrecision = arma::kron(P, arma::diagmat(Tau.row(0)));
  arma::mat Precision = prec * arma::kron(B.t() * B, X.t() * X) + PriorPrecision;
  arma::mat mychol = arma::chol(Precision, "lower");
  arma::mat partialmean = X.t() * Y * B;
  for(arma::uword i = 0; i < Lambda.n_slices; i++){
    partialmean = partialmean - X.t() * arma::diagmat(Gamma.col(i)) * X * Lambda.slice(i).t() * B.t() * B;
  }
  arma::vec mymean = arma::solve(arma::trimatu(mychol.t()), arma::solve(arma::trimatl(mychol), prec * arma::vectorise(partialmean)));
  Theta = arma::trans(arma::reshape(mymean + arma::solve(arma::trimatl(mychol), arma::randn<arma::vec>(Theta.n_rows * Theta.n_cols)), Theta.n_cols, Theta.n_rows));
}

// [[Rcpp::export]]
void updateEta(arma::mat& Y, arma::cube& Lambda, arma::vec& Sigma, arma::mat& Eta, arma::mat& X, arma::mat& B, double prec, arma::mat& Theta){
  arma::uword nsubj = Y.n_rows;
  arma::mat Xtilde(Y.n_cols, Lambda.n_slices);
  arma::mat PriorPrecision = arma::eye<arma::mat>(Lambda.n_slices, Lambda.n_slices);
  arma::mat Precision(Lambda.n_slices, Lambda.n_slices);
  arma::mat mychol;
  for(arma::uword i = 0; i < nsubj; i++){
    for(arma::uword k = 0; k < Lambda.n_slices; k++){
      Xtilde.col(k) = B * Lambda.slice(k) * X.row(i).t();
    }
    Precision = prec * Xtilde.t() * Xtilde + PriorPrecision;
    Eta.row(i) = arma::mvnrnd(prec * arma::inv_sympd(Precision) * Xtilde.t() * (Y.row(i).t() - B * Theta * X.row(i).t()), arma::inv_sympd(Precision)).t();
    //mychol = arma::chol(Precision, "lower");
    //Eta.row(i) = arma::trans(arma::solve(arma::trimatu(mychol.t()), arma::solve(arma::trimatl(mychol), prec * Xtilde.t() * (Y.row(i).t() - B * Theta * X.row(i).t()))) +
     //+ arma::solve(mychol, arma::randn<arma::vec>(Lambda.n_slices)));
  }
  //for(arma::uword i = 0; i < Eta.n_cols; i++){
  //  Eta.col(i) = Eta.col(i) / arma::stddev(Eta.col(i));
  //}
}

// [[Rcpp::export]]
void updateEta2(arma::mat& Y, arma::cube& Lambda, arma::vec& Sigma, arma::mat& Eta, arma::mat& X, arma::mat& B, double prec, arma::mat& Theta){
  arma::uword n = Y.n_rows;
  arma::uword K = Lambda.n_slices;
  arma::mat Ltilde(Lambda.n_rows, K);
  arma::vec EtaM(K);
  arma::mat EtaV(K, K);
  for(arma::uword i = 0; i < n; i++){
    Ltilde.zeros();
    for(arma::uword k = 0; k < K; k++){
      Ltilde.col(k) = Lambda.slice(k) * arma::trans(X.row(i));
    }
    EtaV = arma::inv_sympd(arma::eye<arma::mat>(K, K) + prec * arma::trans(Ltilde) * arma::trans(B) * B * Ltilde);
    EtaM = prec * EtaV * arma::trans(Ltilde) * arma::trans(B) * (arma::trans(Y.row(i)) - B * Theta * arma::trans(X.row(i)));
    Eta.row(i) = arma::mvnrnd(EtaM, EtaV).t();
  }
}

// [[Rcpp::export]]
void updateEta3(arma::mat& Y, arma::cube& Lambda, arma::mat& Eta, arma::mat& X, arma::mat& B, double prec, arma::mat& Theta){
  arma::uword n = Y.n_rows;
  arma::uword K = Lambda.n_slices;
  arma::vec fullmean;
  double fullprec;
  arma::vec condmean;
  double condprec;
  for(arma::uword i = 0; i < n; i++){
    fullmean = B * Theta * X.row(i).t();
    fullprec = 1;
    for(arma::uword kt = 0; kt < K; kt++){
      fullmean = fullmean + B * Lambda.slice(kt) * X.row(i).t() * Eta(i, kt);
      //fullprec = fullprec + prec * arma::as_scalar(X.row(i) * Lambda.slice(kt).t() * B.t() * B * Lambda.slice(kt) * X.row(i).t());
    }
    for(arma::uword k = 0; k < K; k++){
      condmean = fullmean - B * Lambda.slice(k) * X.row(i).t() * Eta(i, k);
      //condprec = fullprec - prec * arma::as_scalar(X.row(i) * Lambda.slice(k).t() * B.t() * B * Lambda.slice(k) * X.row(i).t());
      condprec = 1 + prec * arma::as_scalar(X.row(i) * Lambda.slice(k).t() * B.t() * B * Lambda.slice(k) * X.row(i).t());
      Eta(i, k) = arma::as_scalar(prec/condprec * X.row(i) * Lambda.slice(k).t() * B.t() * (Y.row(i).t()-condmean)) + 1/std::sqrt(double(condprec)) * R::rnorm(0, 1);
      fullmean = condmean + B * Lambda.slice(k) * X.row(i).t() * Eta(i, k);
    }
  }
}
// [[Rcpp::export]]
double updatePrec(arma::mat& Y, arma::cube& Lambda, arma::mat Gamma, arma::mat& X, arma::mat& B, arma::mat& Theta){
  arma::uword nsubj = Y.n_rows;
  arma::uword k = Lambda.n_slices;
  arma::vec mean;
  //double a = .0001;
  //double b = .0001;
  double a = 0;
  double b = 0;
  double mycumsum = 0;
  double totalN = Y.n_elem;
  for(arma::uword i = 0; i < nsubj; i++){
    mean = B * Theta * X.row(i).t();
    for(arma::uword j = 0; j < k; j++){
      mean = mean + Gamma(i, j) * B * Lambda.slice(j) * X.row(i).t();
    }
    mycumsum = mycumsum + arma::sum((Y.row(i).t() - mean) % (Y.row(i).t() - mean));
  }
  return(R::rgamma(a + totalN/2, 1.0 / (b + 1.0 / 2.0 * mycumsum)));
}



// [[Rcpp::export]]
void updateTau(arma::mat& Theta, arma::cube& Lambda, arma::mat& Tau){
  double t_alpha = .01;
  double t_beta = .01;
  //double t_beta = 1;
  arma::mat P = getPenalty2(Lambda.n_rows, 2);
  arma::uword R = Lambda.n_slices;
  arma::uword D = Lambda.n_cols;
  for(arma::uword d = 0; d < D; d++){
    Tau(0, d) = R::rgamma(t_alpha + (Lambda.n_rows - 1.0) / 2.0, 1.0 / (t_beta + 1.0 / 2.0 * arma::as_scalar(Theta.col(d).t() * P * Theta.col(d))));
  }
  for(arma::uword r = 0; r < R; r++){
    for(arma::uword d = 0; d < D; d++){
      Tau(r + 1, d) = R::rgamma(t_alpha + (Lambda.n_rows - 1.0) / 2.0, 1.0 / (t_beta + 1.0 / 2.0 * arma::as_scalar(Lambda.slice(r).col(d).t() * P * Lambda.slice(r).col(d))));
    }
  }
}

// [[Rcpp::export]]
void updateSigma(arma::mat& Eta, arma::vec& Sigma){
  arma::uword K = Sigma.n_elem;
  double n = Eta.n_rows;
  double a = 1;
  double b = 1;
  for(arma::uword k = 0; k < K; k++){
    //Sigma(k) = gam_trunc_left(a + n / 2.0, 1.0 / (b + 1.0 / 2.0 * arma::sum(arma::square(Eta.col(k)))), Sigma(k - 1));
    Sigma(k) = R::rgamma(a + n / 2.0, 1.0 / (b + 1.0 / 2.0 * arma::sum(arma::square(Eta.col(k)))));
  }
}

void updatec(arma::cube& Lambda, arma::mat& c){
  double a = .1;
  double b = .1;
  arma::uword K = Lambda.n_slices;
  arma::uword D = Lambda.n_cols;
  for(arma::uword k = 0; k < K; k++){
    for(arma::uword d = 0; d < D; d++){
      c(k, d) = R::rgamma(a + Lambda.n_rows / 2.0, 1.0 / (b + 1.0 / 2.0 * arma::sum(arma::square(Lambda.slice(k).col(d)))));
      /*
      if(k == 0){
        c(k, d) = R::rgamma(a + Lambda.n_rows / 2.0, 1.0 / (b + 1.0 / 2.0 * arma::sum(arma::square(Lambda.slice(k).col(d)))));
      } else {
      c(k, d) = gam_trunc_left(a + Lambda.n_rows / 2.0, 1.0 / (b + 1.0 / 2.0 * arma::sum(arma::square(Lambda.slice(k).col(d)))), c(k - 1, d));
      }
      */
    }
  }
}

double gam_trunc_left(double a, double b,  double cut){ 
  double u, pcut, y; 
  pcut = R::pgamma(cut,a,b,1,0);
  if(pcut>0.99){
    return(cut);
  } 
  u = R::runif(0.0,1.0); 
  u = pcut + (1-pcut)*u; 
  y = R::qgamma(u, a, b, 1, 0);
  return y; 
} 