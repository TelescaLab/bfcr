#include "Utils.h"
#ifdef _OPENMP
 #include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
using namespace Rcpp;

// [[Rcpp::export]]
double get_proposal(double old) {
  double proposal;
  while(true) {
    proposal = old + R::rnorm(0, 1);
    if (proposal > 0) break;
  }
  return proposal;
}

// [[Rcpp::export]]
arma::uvec armadillo_modulus(arma::uvec indicies, arma::uword n){
  return(indicies - n * arma::floor(indicies / n));
}

List extract_eigenfn(arma::cube& Lambda,
                     arma::mat& Psi, arma::mat& Psi_sqrt,
                     arma::mat& Psi_sqrt_inv, arma::mat& B,
                     arma::uword eigenvals, arma::vec z,
                     arma::vec time){
  arma::uword dim_latent = Lambda.n_rows;
  arma::uword dim_spline = B.n_rows;
  arma::mat eigenfn_latent(dim_latent, dim_latent);
  arma::mat eigenfn_spline(dim_spline, dim_latent);
  arma::mat eigenfn_spline_ordered(dim_spline, eigenvals);
  arma::vec eigenval_latent(dim_latent);
  arma::mat eigenfn_latent_ordered(dim_latent, eigenvals);
  arma::vec eigenval_spline(eigenvals);
  arma::vec eigenval_pve(eigenvals);
  arma::mat cov_latent = arma::zeros<arma::mat>(dim_latent, dim_latent);
  arma::mat cov_spline;
  for(arma::uword k = 0; k < Lambda.n_slices; k++){
    cov_latent = cov_latent + Lambda.slice(k) * z * z.t() * Lambda.slice(k).t();
  }
  arma::mat cov_latent_transformed = Psi_sqrt * cov_latent * Psi_sqrt;
  arma::eig_sym(eigenval_latent, eigenfn_latent, cov_latent_transformed);
  eigenfn_spline =  B * Psi_sqrt_inv * eigenfn_latent;
  cov_spline = B * cov_latent * B.t();
  double eval_sum = arma::sum(eigenval_latent);
  for(arma::uword v = 0; v < eigenvals; v++){
    eigenfn_spline_ordered.col(v) = eigenfn_spline.col(dim_latent - v - 1);
    eigenval_spline(v) = eigenval_latent(dim_latent - v - 1);
    eigenval_pve(v) = eigenval_spline(v) / eval_sum;
  }
  arma::vec eval_cumsum = arma::cumsum(eigenval_pve);
  eigenval_pve.elem(arma::find(eval_cumsum > 1)).zeros();
  double magnitude = arma::sum(eigenval_latent);
  return(List::create(Named("eigenfn_spline", eigenfn_spline_ordered),
                      Named("eigenval", eigenval_spline),
                      Named("eigenval_pve", eigenval_pve),
                      Named("cov_spline", cov_spline),
                      Named("magnitude", magnitude)));
}