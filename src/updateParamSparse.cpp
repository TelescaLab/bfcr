#include <RcppArmadillo.h>
#include "Utils.h"
using namespace Rcpp;

// [[Rcpp::export]]
void updateLambdaS(arma::field<arma::vec>& Y, arma::cube& Lambda, arma::mat& Tau, arma::mat& c, arma::mat& Gamma, arma::mat& X, arma::field<arma::mat>& B, double prec, arma::mat& Theta){
  arma::mat P = getPenalty(Lambda.n_rows);
  arma::mat PriorPrecision;
  arma::mat Precision;
  arma::mat mychol;
  arma::vec mymean;
  arma::vec Ltemp;
  arma::uword K = Lambda.n_slices;
  arma::uword counter = 0;
  arma::mat condMean = arma::zeros<arma::mat>(Theta.n_rows, Theta.n_cols);
  for(arma::uword k = 0; k < K; k++){
    Precision = arma::kron(P, arma::diagmat(Tau.row(k+1))) + arma::kron(arma::eye<arma::mat>(Lambda.n_rows, Lambda.n_rows),
                           arma::diagmat(c.row(k)));
    for(arma::uword i = 0; i < X.n_rows; i++){
      condMean = condMean + Gamma(i, k) * B(i).t() * (Y(i) - B(i) * Theta * X.row(i).t()) * X.row(i);
      Precision = Precision + prec * Gamma(i, k) * Gamma(i, k) * arma::kron(B(i).t() * B(i), X.row(i).t() * X.row(i));
      //Precision = Precision + prec * Gamma(i, k) * Gamma(i, k) * arma::kron(X.row(i).t() * X.row(i), B(i).t() * B(i));
      for(arma::uword q = 0; q < K; q++){
        if(q != k){
          condMean = condMean - Gamma(i, k) * Gamma(i, q) * B(i).t() * B(i) * Lambda.slice(q) * X.row(i).t() * X.row(i);
        }
      }
    }
    condMean = prec * condMean;
    mychol = arma::chol(Precision, "lower");
    mymean = arma::solve(arma::trimatu(mychol.t()), arma::solve(arma::trimatl(mychol), arma::vectorise(condMean.t())));
    /*
    Lambda.slice(k) = arma::reshape(mymean + arma::solve(arma::trimatl(mychol),
                                    arma::randn<arma::vec>(Lambda.n_rows * Lambda.n_cols)),
    Lambda.n_rows, Lambda.n_cols);
    */
    //Lambda.slice(k) = arma::reshape(mymean,
    //                              Lambda.n_rows, Lambda.n_cols);
    Lambda.slice(k) = arma::trans(arma::reshape(mymean + arma::solve(arma::trimatl(mychol), 
                                                arma::randn<arma::vec>(Lambda.n_rows *
                                                  Lambda.n_cols)),Lambda.n_cols,Lambda.n_rows));    
    condMean.zeros();
    counter++;
  }
  //Rcpp::Rcout << Lambda << std::endl;
}

// [[Rcpp::export]]
void updateThetaS(arma::field<arma::vec>& Y, arma::cube& Lambda, arma::mat& Tau, arma::mat& Gamma,
                  arma::mat& X, arma::field<arma::mat> B, double prec, arma::mat& Theta){
  arma::mat P = getPenalty(Theta.n_rows);
  arma::mat PriorPrecision = arma::kron(P, arma::diagmat(Tau.row(0)));
  arma::mat Precision = PriorPrecision;
  for(arma::uword i = 0; i < X.n_rows; i++){
    Precision = Precision + prec * (arma::kron(X.row(i).t() * X.row(i), B(i).t() * B(i)));
  }
  arma::mat mychol = arma::chol(Precision, "lower");
  arma::mat partialmean = arma::zeros<arma::mat>(Theta.n_rows, Theta.n_cols);
  for(arma::uword i = 0; i < X.n_rows; i++){
    partialmean = partialmean + B(i).t() * Y(i) * X.row(i);
    for(arma::uword k = 0; k < Lambda.n_slices; k++){
      partialmean = partialmean - Gamma(i, k) * B(i).t() * B(i) * Lambda.slice(k) * X.row(i).t() * X.row(i);
    }
  }
  arma::vec mymean = arma::solve(arma::trimatu(mychol.t()), arma::solve(trimatl(mychol), prec * arma::vectorise(partialmean)));
  Theta = arma::reshape(mymean + arma::solve(arma::trimatl(mychol), arma::randn<arma::vec>(Theta.n_rows * Theta.n_cols)), Theta.n_rows, Theta.n_cols);
}

// [[Rcpp::export]]
void updateEtaS(arma::field<arma::vec>& Y, arma::cube& Lambda, arma::vec& Sigma, arma::mat& Eta, arma::mat& X, arma::field<arma::mat>& B, double prec, arma::mat& Theta){
  arma::uword nsubj = X.n_rows;
  arma::field<arma::mat> Xtilde(X.n_rows);
  arma::mat PriorPrecision = arma::diagmat(Sigma);
  arma::mat Precision(Lambda.n_slices, Lambda.n_slices);
  arma::mat mychol;
  for(arma::uword i = 0; i < nsubj; i++){
    Xtilde(i) = arma::mat(B(i).n_rows, Lambda.n_slices);
    for(arma::uword k = 0; k < Lambda.n_slices; k++){
      Xtilde(i).col(k) = B(i) * Lambda.slice(k) * X.row(i).t();
    }
    Precision = prec * Xtilde(i).t() * Xtilde(i) + PriorPrecision;
    mychol = arma::chol(Precision, "lower");
    Eta.row(i) = arma::trans(arma::solve(arma::trimatu(mychol.t()), arma::solve(arma::trimatl(mychol), prec * Xtilde(i).t() * (Y(i) - B(i) * Theta * X.row(i).t()))) +
      + arma::solve(arma::trimatl(mychol), arma::randn<arma::vec>(Lambda.n_slices)));
  }
}

// [[Rcpp::export]]
double updatePrecS(arma::field<arma::vec>& Y, arma::cube& Lambda, arma::mat Gamma, arma::mat& X, arma::field<arma::mat>& B, arma::mat& Theta){
  arma::uword nsubj = X.n_rows;
  arma::uword K = Lambda.n_slices;
  arma::field<arma::vec> mean(nsubj);
  double a = .001;
  double b = .001;
  double mycumsum = 0;
  double totalN = 0;
  for(arma::uword i = 0; i < nsubj; i++){
    mean(i) = B(i) * Theta * X.row(i).t();
    for(arma::uword k = 0; k < K; k++){
      mean(i) = mean(i) + Gamma(i, k) * B(i) * Lambda.slice(k) * X.row(i).t();
    }
    mycumsum = mycumsum + arma::sum((Y(i) - mean(i)) % (Y(i) - mean(i)));
    totalN = totalN + Y(i).n_elem;
  }
  return(R::rgamma(a + totalN / 2.0, 1.0 / (b + 1.0 / 2.0 * mycumsum)));
}