#include <RcppArmadillo.h>
#include <progress.hpp>

static double const log2pi = std::log(2.0 * M_PI);

/* C++ version of the dtrmv BLAS function */
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

// [[Rcpp::export]]
double dmvnrm_arma_fast(arma::mat const &x,  
                           arma::rowvec const &mean,  
                           arma::mat const &sigma, 
                           bool const logd = false) { 
  using arma::uword;
  uword const n = x.n_rows, 
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())), 
    constants = -(double)xdim/2.0 * log2pi, 
    other_terms = rootisum + constants;
  
  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);     
  }  
  
  if (logd)
    return arma::sum(out);
  return exp(arma::sum(out));
}

// [[Rcpp::export]]
void calculate_waic_cpp(Rcpp::List& mcmc_output, Rcpp::List& observed_time,
                        arma::mat& waic_return) {
  Rcpp::List samples = mcmc_output["samples"];
  Rcpp::List control = mcmc_output["control"];
  Rcpp::List data = mcmc_output["data"];
  arma::mat response = data["response"];
  arma::mat design_mean = data["design_mean"];
  arma::mat design_var = data["design_var"];
  arma::mat basis = data["basis"];
  arma::uword kdim = data["latent_dimension"];
  arma::uword iterations = control["iterations"];
  arma::uword burnin = control["burnin"];
  arma::cube beta = samples["beta"];
  arma::field<arma::cube> lambda = samples["lambda"];
  arma::cube eta = samples["eta"];
  arma::mat varphi = samples["varphi"];
  arma::uword num_subjects = response.n_rows;
  arma::field<arma::mat> basis_sub(num_subjects);
  arma::field<arma::vec> mu(num_subjects);
  arma::field<arma::mat> sigma(num_subjects);
  arma::field<arma::uvec> indices(num_subjects);
  arma::field<arma::vec> response_vec(num_subjects);
  arma::uvec num_time(num_subjects);
  Progress progress_bar(num_subjects, true);
  arma::uword counter = 0;
  for (arma::uword s = 0; s < num_subjects; s++) {
    arma::uvec indices_curr = observed_time[s];
    indices(s) = indices_curr;
    basis_sub(s) = basis.rows(indices_curr);
    mu(s) = arma::vec(indices_curr.n_elem);
    sigma(s) = arma::mat(indices_curr.n_elem, indices_curr.n_elem);
    arma::rowvec response_curr = response.row(s);
    response_vec(s) = response_curr.elem(indices_curr);
    num_time(s) = indices_curr.n_elem;
  }
  
  for (arma::uword s = 0; s < num_subjects; s++) {
    for (arma::uword iter = burnin; iter < iterations; iter++) {
      
      if (Progress::check_abort()){
        Rcpp::Rcout << "Aborted\n";
        goto stop;
      }
      mu(s) = basis_sub(s) * beta.slice(iter) * design_mean.row(s).t();
      sigma(s).zeros();
      for (arma::uword k = 0; k < kdim; k++) {
        sigma(s) = sigma(s) + basis_sub(s) * lambda(iter).slice(k) * 
          design_var.row(s).t() * design_var.row(s) * 
          lambda(iter).slice(k).t() * basis_sub(s).t();
      }
      sigma(s) = sigma(s) + 
        1 / varphi(s, iter) * arma::eye<arma::mat>(num_time(s), num_time(s));
      waic_return(counter, s) = dmvnrm_arma_fast(response_vec(s).t(),
                  mu(s).t(), sigma(s), true);
      counter++;
    }
    progress_bar.increment();
    counter = 0;
  }
  
  stop:
    NULL;
}