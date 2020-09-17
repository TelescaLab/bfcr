#include "GaussTransformations.h"

GaussTransformations::GaussTransformations(arma::uword dim) {
  arma::mat chol = arma::mat(dim, dim);
  arma::vec w = arma::vec(dim), mu = arma::vec(dim), v = arma::vec(dim),
    z = arma::vec(dim), result = arma::vec(dim);
}

arma::vec GaussTransformations::bayes_reg(arma::mat& precision, arma::vec& g) {
  chol = arma::chol(precision, "lower");
  w = solve(arma::trimatl(chol), g, arma::solve_opts::fast);
  mu = solve(arma::trimatu(chol.t()), w, arma::solve_opts::fast);
  z = arma::randn<arma::vec>(g.n_elem);
  v = arma::solve(arma::trimatu(chol.t()), z, arma::solve_opts::fast);
  result = mu + v;
  return(result);
}