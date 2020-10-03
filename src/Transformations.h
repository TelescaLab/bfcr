#ifndef TRANSFORMATIONS_H
#define TRANSFORMATIONS_H

#include <RcppArmadillo.h>
#include "Parameters.h"
#include "Data.h"
#include "GaussTransformations.h"

class Parameters;

class Transformations {
public:
  GaussTransformations beta_gauss, lambda_gauss, eta_gauss;
  arma::mat beta_precision;
  arma::vec beta_g, beta_result;
  arma::mat lambda_precision;
  arma::vec lambda_g, lambda_result;
  arma::mat eta_precision;
  arma::vec eta_g, eta_result;
  arma::mat btb, bty, fit, lambda_tilde;
  arma::mat blk_diag_mean_penalties;
  arma::cube blk_diag_var_penalties;
  arma::mat psi_mat;
  arma::mat lin_constr;
  void complete_response(Data&, Parameters&);
  void build_blk_diag_mean(Data&, Parameters&);
  void build_blk_diag_var(Data&, Parameters&);
  void build_blk_diag_phi_delta(Data&, Parameters&);
  Transformations(Data&, Parameters&);
  Transformations() {};
  arma::mat fit_beta_removed;
  arma::mat fit_lambda_removed;
  arma::mat eta_transf;
  arma::mat fit_lambda, fit_beta, lambda_old;
  arma::mat phi_lambda_sum;
  arma::mat delta_cumprod;
  arma::cube blk_diag_phi_delta;
  arma::mat blk_diag_delta_cumprod;
  arma::vec squared_diff;
};

#endif
