#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat DiffOp(arma::uword n){
  arma::mat D = arma::eye(n, n);
  for(arma::uword i = 1; i < n; i++){
    D(i, i - 1) = -1;
  }
  return(D);
}

// [[Rcpp::export]]
arma::mat getPenalty2(arma::uword n, arma::uword D){
  arma::mat DiffMat = DiffOp(n);
  for(arma::uword d = 1; d < D; d++){
    DiffMat = DiffMat * DiffMat;
  }
  return(arma::trans(DiffMat) * DiffMat);
}

// [[Rcpp::export]]
arma::mat getPenalty(arma::uword n){
  arma::mat Penalty = arma::zeros<arma::mat>(n, n);
  for(arma::uword i = 0; i < n - 1; i++){
    Penalty(i + 1, i) = -1.0;
  }
  for(arma::uword i = 0; i < n - 1; i++){
    Penalty(i, i + 1) = -1.0;
  }
  for(arma::uword i = 1; i < n - 1; i++){
    Penalty(i, i) = 2.0;
  }
  Penalty(0, 0) = 1.0;
  Penalty(n - 1, n - 1) = 1.0;
  
  return(Penalty);
}

// [[Rcpp::export]]
arma::rowvec initializeY(arma::vec y, arma::vec observedTimes, arma::vec fullTimes){
  arma::uvec result(observedTimes.n_elem);
  arma::uword J = 0;
  for(arma::uword i = 0; i < observedTimes.n_elem; i++){
    for(arma::uword j = J; j < fullTimes.n_elem; j++){
      if(observedTimes(i) == fullTimes(j)){
        result(i) = j;
        J = j;
      }
    }
  }
  arma::vec yfull = arma::zeros<arma::vec>(fullTimes.n_elem);
  yfull(result) = y;
  return(yfull.t());
}

// [[Rcpp::export]]
arma::uvec getObservedOrder(arma::vec observedTimes, arma::vec fullTimes){
  arma::uvec result(observedTimes.n_elem);
  arma::uword J = 0;
  for(arma::uword i = 0; i < observedTimes.n_elem; i++){
    for(arma::uword j = J; j < fullTimes.n_elem; j++){
      if(observedTimes(i) == fullTimes(j)){
        result(i) = j;
        J = j;
      }
    }
  }
  return(result);
}

// [[Rcpp::export]]
void PredictY(arma::mat& ImputedY, arma::mat X, arma::mat B, arma::mat Theta, arma::mat Eta, arma::cube Lambda, double Prec){
  arma::mat mymean = X * Theta.t() * B.t();
  Rcpp::Rcout << arma::size(X) << std::endl;
  
  for(arma::uword k = 0; k < Lambda.n_slices; k++){
    mymean = mymean + arma::diagmat(Eta.col(k)) * X * Lambda.slice(k).t() * B.t();
  }
  mymean = mymean + (1.0 / sqrt(Prec)) * arma::randn(mymean.n_rows, mymean.n_rows);
  ImputedY = mymean;
}

// [[Rcpp::export]]
void PredictY2(arma::mat& ImputedY, arma::field<arma::uvec> observedOrder, arma::mat X, arma::mat B, arma::mat Theta, arma::cube Lambda, arma::mat Eta, double Prec){
  
  arma::mat temp = X * Theta.t() * B.t();
  for(arma::uword k = 0; k < Lambda.n_slices; k++){
    temp = temp + arma::diagmat(Eta.col(k)) * X * Lambda.slice(k).t() * B.t();
  }
  temp = temp + (1.0 / sqrt(Prec)) * arma::randn(X.n_rows, B.n_rows);
  for(arma::uword i = 0; i < X.n_rows; i++){
    for(arma::uword j = 0; j < observedOrder(i).n_elem; j++){
      temp(i, observedOrder(i)(j)) = ImputedY(i, observedOrder(i)(j));
    }
  }
  ImputedY = temp;
}
arma::uvec test(arma::field<arma::vec> observedTimes, arma::vec fullTimes){
  return(getObservedOrder(observedTimes(0), fullTimes));
}
