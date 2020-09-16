#ifndef DATA_H
#define DATA_H
#include <RcppArmadillo.h>
#include "Utils.h"

class Data {
public:
  arma::uword d1, d2, basis_dim, kdim, n;
  arma::uvec missing, missing_sub, missing_time;
  arma::mat response, design_mean;
  arma::mat design_var, basis;
  arma::field<arma::mat> penalties_mean;
  arma::field<arma::mat> penalties_var;
  arma::uvec indices_mean, indices_var;
  arma::uword n_smooths_mean, n_smooths_var;
  arma::uword iter, thin, burnin;
  arma::vec time;
  arma::uvec seq_along_start_delta, seq_along_end_delta;
  arma::ivec rank_mean, rank_var;
  arma::uvec seq_along_elongated;
  Data(arma::mat&, arma::mat&, arma::mat&,
       arma::mat&, arma::vec&,
       arma::field<arma::mat>&,arma::field<arma::mat>&,
       arma::uvec&, arma::uvec&, arma::uword, arma::uword, 
       arma::uword, arma::uword);
  Rcpp::List write_data();
  arma::umat get_seq_along_tau(arma::field<arma::mat>&, arma::uvec&);
  arma::umat seq_along_tau1;
  arma::umat seq_along_tau2;
};

#endif