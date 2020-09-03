#ifndef GAUSSTRANSFORMATIONS_H
#define GAUSSTRANSFORMATIONS_H

#include <RcppArmadillo.h>

class GaussTransformations {
public:
  arma::mat chol;
  arma::vec w, mu, z, v, result;
  GaussTransformations(arma::uword dim);
  GaussTransformations() {};
  arma::vec bayes_reg(arma::mat& precision, arma::vec& mu);
};

#endif