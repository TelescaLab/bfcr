#include "Data.h"
Rcpp::List Data::write_data() {
  return(Rcpp::List::create(Rcpp::Named("response", response),
                            Rcpp::Named("basis", basis),
                            Rcpp::Named("time", time),
                            Rcpp::Named("latent_dimension", kdim),
                            Rcpp::Named("missing_subjects", missing_sub),
                            Rcpp::Named("missing_time", missing_time),
                            Rcpp::Named("design_mean", design_mean),
                            Rcpp::Named("design_var", design_var)));
}

Data::Data(arma::mat& response, arma::mat& design_mean,
           arma::mat& design_var, arma::mat& basis,
           arma::vec& time,
           arma::field<arma::mat>& penalties_mean, 
           arma::field<arma::mat>& penalties_var,
           arma::uvec& indices_mean, arma::uvec& indices_var,
           arma::uword kdim, arma::uword iter,
           arma::uword burnin, arma::uword thin) {
  d1 = design_mean.n_cols;
  d2 = design_var.n_cols;
  basis_dim = basis.n_cols;
  n = response.n_rows;
  n_smooths_mean = arma::size(arma::find_unique(indices_mean))(0);
  n_smooths_var = arma::size(arma::find_unique(indices_var))(0);
  Rcpp::Rcout << "n_smooths_var: " << n_smooths_var << std::endl;
  this->response = response;
  this->design_mean = design_mean;
  this->design_var = design_var;
  this->basis = basis;
  this->time = time;
  this->penalties_mean = penalties_mean;
  this->penalties_var = penalties_var;
  this->indices_mean = indices_mean;
  this->indices_var = indices_var;
  this->kdim = kdim;
  this->iter = iter;
  this->burnin = burnin;
  this->thin = thin;
  missing = arma::find_nonfinite(response);
  missing_sub = armadillo_modulus3(missing, response.n_rows);
  missing_time = arma::floor(missing / response.n_rows);
  this->response.elem(missing).fill(0);
  seq_along_start = arma::uvec(n_smooths_var);
  seq_along_end = arma::uvec(n_smooths_var);
  int old_index = -1;
  int start = 0;
  int end = -1;
  arma::uword counter = 0;
  seq_along_start(0) = start;
  seq_along_end(0) = end;
  for (arma::uword i = 0; i < penalties_var.n_elem; i++) {
    if (indices_var(i) != old_index) {
      end = end + penalties_var(i).n_rows / basis_dim;
      seq_along_start(counter) = start;
      seq_along_end(counter) = end;
      start = end + 1;
      counter++;
      old_index = indices_var(i);
    }
  }
  Rcpp::Rcout << "in data construction\n";
  Rcpp::Rcout << arma::size(seq_along_start) << "\n";
  Rcpp::Rcout << arma::size(seq_along_end) << "\n";
}