#ifndef MODSTRING_H1
#define MODSTRING_H1

#include <RcppArmadillo.h>
using namespace Rcpp;

double get_proposal(double old);
arma::uvec armadillo_modulus(arma::uvec indicies, arma::uword n);
List extract_eigenfn(arma::cube& Lambda,
                     arma::mat& Psi, arma::mat& Psi_sqrt,
                     arma::mat& Psi_sqrt_inv, arma::mat& B,
                     arma::uword eigenvals, arma::vec z,
                     arma::vec time);

#endif