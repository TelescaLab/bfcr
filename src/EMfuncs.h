#ifndef MODSTRING_H2
#define MODSTRING_H2

#include <RcppArmadillo.h>
using namespace Rcpp;

Rcpp::List cpp_EM2(arma::mat X, arma::mat B, arma::mat Y, arma::uword K,
                   double tol, arma::uword max_iter);
void completeY2Means(arma::mat& Y, arma::uvec missing_sub,
                     arma::uvec missing_time);
void cppgetX(arma::mat &EtaM, arma::cube &EtaV, arma::mat &X, arma::mat &newX,
             int cores = 1);
void cppupdateall(arma::mat &Theta, arma::mat &Lambda, arma::vec &precision,
                  arma::mat &newX, arma::mat &B, arma::mat &newY,
                  arma::uword K);
void cppupdateeta(arma::mat &Theta, arma::mat &Lambda, arma::vec &precision,
                  arma::mat &EtaM, arma::cube &EtaV, arma::mat &X,
                  arma::mat &B, arma::mat &Y, double K);

#endif