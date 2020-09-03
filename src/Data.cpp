#include "Data.h"
Rcpp::List Data::write_data() {
  return(Rcpp::List::create(Rcpp::Named("response", response),
                            Rcpp::Named("basis", basis),
                            Rcpp::Named("latent_dimension", kdim),
                            Rcpp::Named("missing_subjects", missing_sub),
                            Rcpp::Named("missing_time", missing_time),
                            Rcpp::Named("design_mean", design_mean),
                            Rcpp::Named("design_var", design_var)));
}

Data::Data(arma::mat& response, arma::mat& design_mean,
           arma::mat& design_var, arma::mat& basis,
           arma::field<arma::mat>& penalties_mean, 
           arma::field<arma::mat>& penalties_var,
           arma::uvec& indices_mean, arma::uvec& indices_var,
           arma::uword kdim, arma::uword iter,
           arma::uword thin) {
  d1 = design_mean.n_cols;
  d2 = design_var.n_cols;
  basis_dim = basis.n_cols;
  n = response.n_rows;
  this->response = response;
  this->design_mean = design_mean;
  this->design_var = design_var;
  this->basis = basis;
  this->penalties_mean = penalties_mean;
  this->penalties_var = penalties_var;
  this->indices_mean = indices_mean;
  this->indices_var = indices_var;
  this->kdim = kdim;
  this->iter = iter;
  this->thin = thin;
  missing = arma::find_nonfinite(response);
  missing_sub = armadillo_modulus3(missing, response.n_rows);
  missing_time = arma::floor(missing / response.n_rows);
  this->response.elem(missing).fill(0);
}