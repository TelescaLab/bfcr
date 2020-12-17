#include <RcppArmadillo.h>
#include "Utils.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List get_posterior_subject_bands_cpp(List mcmc_output,
                                           double alpha,
                                           std::string mode){
  Rcpp::List samples = mcmc_output["samples"];
  Rcpp::List data = mcmc_output["data"];
  Rcpp::List control = mcmc_output["control"];
  arma::cube beta = samples["beta"];
  arma::mat varphi = samples["varphi"];
  arma::field<arma::cube> lambda = samples["lambda"];
  arma::cube eta = samples["eta"];
  arma::uword kdim = data["latent_dimension"];
  arma::mat basis = data["basis"];
  arma::mat response = data["response"];
  arma::mat design_mean = data["design_mean"];
  arma::mat design_var = data["design_var"];
  arma::uword iterations = control["iterations"];
  arma::uword burnin = control["burnin"];
  arma::uword num_subjects = response.n_rows;
  arma::running_stat_vec<arma::vec> stats;
  arma::mat m_alpha(num_subjects, iterations - burnin);
  arma::vec m_star(num_subjects);
  arma::mat current_lower_dim_fit;
  arma::mat current_fit;

  arma::vec response_vectorized = arma::vectorise(arma::trans(response));
  for (arma::uword iter = burnin; iter < iterations; iter++) {

    current_lower_dim_fit = beta.slice(iter) * design_mean.t();
    for (arma::uword k = 0; k < kdim; k++) {
      current_lower_dim_fit = current_lower_dim_fit +
        lambda(iter).slice(k) * design_var.t() *
        arma::diagmat(eta.slice(iter).col(k));
    }
    current_fit = basis * current_lower_dim_fit;
    stats(arma::vectorise(current_fit));
  }

  arma::mat subj_mean = arma::reshape(stats.mean(), basis.n_rows, num_subjects);
  arma::mat subj_sd = arma::reshape(stats.stddev(), basis.n_rows, num_subjects);
  m_star.fill(R::qnorm5(1 - alpha / 2, 0, 1, 1, 0));
  if (mode == "simultaneous") {
    arma::uword counter = 0;
    for (arma::uword iter = burnin; iter < iterations; iter++) {
      current_lower_dim_fit = beta.slice(iter) * design_mean.t();
      for (arma::uword k = 0; k < kdim; k++) {
        current_lower_dim_fit = current_lower_dim_fit +
          lambda(iter).slice(k) * design_var.t() *
          arma::diagmat(eta.slice(iter).col(k));
      }
      current_fit = basis * current_lower_dim_fit;
      m_alpha.col(counter) =
        arma::max(arma::abs(current_fit - subj_mean) /
          subj_sd, 0).t();
      counter++;
    }
    arma::vec alpha_vec = {1-alpha};
    m_star = arma::quantile(m_alpha, alpha_vec, 1);
  }

  arma::mat mean_overall = subj_mean;
  arma::mat lower = mean_overall - subj_sd * arma::diagmat(m_star);
  arma::mat upper = mean_overall + subj_sd * arma::diagmat(m_star);
  Rcpp::List posterior_bands =
    Rcpp::List::create(Rcpp::Named("subject_means", mean_overall),
                       Rcpp::Named("subject_lower", lower),
                       Rcpp::Named("subject_upper", upper),
                       Rcpp::Named("current_lower_dim_fit",
                                   current_lower_dim_fit));
  return(posterior_bands);
}

// [[Rcpp::export]]
arma::mat get_posterior_means_cpp(List mcmc_results, arma::vec xi,
                                          double alpha,
                                          std::string mode){
  Rcpp::List data = mcmc_results["data"];
  Rcpp::List samples = mcmc_results["samples"];
  Rcpp::List control = mcmc_results["control"];
  arma::cube beta = samples["beta"];
  arma::mat basis = data["basis"];
  arma::uword iterations = control["iterations"];
  arma::uword burnin = control["burnin"];
  arma::uword basis_dim = basis.n_cols;
  arma::vec m_alpha(iterations - burnin);
  arma::running_stat_vec<arma::vec> stats;
  arma::vec m_beta(basis_dim);
  arma::uword num_time_pts = basis.n_rows;
  arma::mat quants(num_time_pts, 3);
  arma::vec current_mean;
  for(arma::uword i = burnin; i < iterations; i++){
    current_mean = basis * beta.slice(i) * xi;
    stats(current_mean);
  }
  double m_star = R::qnorm5(1 - alpha / 2, 0, 1, 1, 0);
  if (mode == "simultaneous") {
    arma::uword counter = 0;
    for(arma::uword i = burnin; i < iterations; i++){
      current_mean = basis * beta.slice(i) * xi;
      m_alpha(counter) =
        arma::max(arma::abs(current_mean - stats.mean()) / stats.stddev());
      counter++;
    }
    arma::vec alpha_vec = {1 - alpha};
    m_star = arma::as_scalar(arma::quantile(m_alpha, alpha_vec));

  }

  quants.col(0) = stats.mean() - m_star * stats.stddev();
  quants.col(1) = stats.mean();
  quants.col(2) = stats.mean() + m_star * stats.stddev();
  arma::vec v = stats.stddev();
  return(quants);
}

// [[Rcpp::export]]
List get_posterior_eigen_cpp(Rcpp::List mcmc_results,
                                     arma::uword eigenvals,
                                     arma::vec zi, double alpha,
                                     std::string mode){
  Rcpp::List control = mcmc_results["control"];
  Rcpp::List data = mcmc_results["data"];
  Rcpp::List samples = mcmc_results["samples"];
  arma::mat basis = data["basis"];
  arma::field<arma::cube> lambda = samples["lambda"];
  arma::cube beta = samples["beta"];
  arma::vec time = data["time"];
  arma::uword iterations = control["iterations"];
  arma::uword burnin = control["burnin"];
  arma::uword num_post_iter = iterations - burnin;
  arma::uword basis_dim = basis.n_cols;
  arma::mat m_alpha(num_post_iter, eigenvals);
  arma::running_stat_vec<arma::vec> stats_vec;
  arma::running_stat_vec<arma::vec> stats_val;
  arma::running_stat_vec<arma::vec> stats_cov;
  arma::mat psi(basis_dim, basis_dim);
  arma::mat psi_sqrt;
  arma::mat psi_sqrt_inv;
  arma::vec q_alpha(eigenvals);
  for(arma::uword j = 0; j < basis_dim; j++){
    for(arma::uword i = 0; i < basis_dim; i++){
      psi(i, j) =
        arma::as_scalar(arma::trapz(time, basis.col(i) % basis.col(j)));
    }
  }
  psi_sqrt = arma::sqrtmat_sympd(psi);
  psi_sqrt_inv = arma::inv_sympd(psi_sqrt);

  arma::mat temp_evec;
  arma::mat eval_mat(num_post_iter, eigenvals);
  arma::mat eval_pve_mat(num_post_iter, eigenvals);
  List eigen_list;
  arma::vec magnitude(num_post_iter);
  arma::mat tempcov;
  arma::uword idx1, idx2;
  arma::uword counter = 0;
  for(arma::uword i = burnin; i < iterations; i++){
    eigen_list = extract_eigenfn(lambda(i), psi, psi_sqrt, psi_sqrt_inv,
                                 basis, eigenvals, zi, time);

    temp_evec = Rcpp::as<arma::mat>(eigen_list["eigenfn_spline"]);
    eval_mat.row(counter) = Rcpp::as<arma::rowvec>(eigen_list["eigenval"]);
    eval_pve_mat.row(counter) = Rcpp::as<arma::rowvec>(eigen_list["eigenval_pve"]);
    arma::rowvec cumsum_vec = arma::cumsum(eval_pve_mat.row(counter));
    stats_cov(arma::vectorise(as<arma::vec>(eigen_list["cov_spline"])));
    magnitude(counter) = eigen_list["magnitude"];
    if (counter == 0) stats_vec(arma::vectorise(temp_evec));
    else {
      // align eigenvectors

      for (arma::uword k = 0; k < eigenvals; k++) {
        idx1 = k * basis.n_rows;
        idx2 = (k + 1) * basis.n_rows - 1;
        if (arma::sum(arma::square(temp_evec.col(k) + stats_vec.mean().subvec(idx1, idx2))) <
          arma::sum(arma::square(temp_evec.col(k) - stats_vec.mean().subvec(idx1, idx2)))) {
          temp_evec.col(k) = -temp_evec.col(k);
        }
      }

      stats_vec(arma::vectorise(temp_evec));
    }
    counter++;
  }
  for (arma::uword k = 0; k < eigenvals; k++) {
    q_alpha(k) = R::qnorm5(1 - alpha / 2, 0, 1, 1, 0);
  }
  if (mode == "simultaneous") {
    counter = 0;
    for (arma::uword i = burnin; i < iterations; i++) {
      eigen_list = extract_eigenfn(lambda(i),
                                   psi,
                                   psi_sqrt,
                                   psi_sqrt_inv,
                                   basis,
                                   eigenvals,
                                   zi,
                                   time);
      temp_evec = Rcpp::as<arma::mat>(eigen_list["eigenfn_spline"]);

      for (arma::uword k = 0; k < eigenvals; k++) {
        idx1 = k * basis.n_rows;
        idx2 = (k + 1) * basis.n_rows - 1;
        if (arma::sum(arma::square(temp_evec.col(k) + stats_vec.mean().subvec(idx1, idx2))) <
          arma::sum(arma::square(temp_evec.col(k) - stats_vec.mean().subvec(idx1, idx2)))) {
          temp_evec.col(k) = -temp_evec.col(k);
        }
        m_alpha(counter, k) = arma::max((temp_evec.col(k) -
          stats_vec.mean().subvec(idx1, idx2)) / stats_vec.stddev().subvec(idx1, idx2));
      }
      counter++;
    }
    arma::vec alpha_vec = {1 - alpha};
    for (arma::uword k = 0; k < eigenvals; k++) {
      q_alpha(k) = arma::as_scalar(arma::quantile(m_alpha.col(k), alpha_vec));
    }
  }


  arma::vec alpha_vec_eval = {alpha / 2.0, 0.5, 1 - alpha / 2.0};
  arma::mat lower(time.n_elem, eigenvals);
  arma::mat mean(time.n_elem, eigenvals);
  arma::mat upper(time.n_elem, eigenvals);
  arma::mat eigenbands(time.n_elem, eigenvals * 3);
  arma::mat eigenval_intervals(3, eigenvals);
  arma::mat eigenval_pve_intervals(3, eigenvals);
  arma::vec magnitude_interval(3);
  for (arma::uword k = 0; k < eigenvals; k++) {
    idx1 = k * basis.n_rows;
    idx2 = (k + 1) * basis.n_rows - 1;
    lower.col(k) = (stats_vec.mean().subvec(idx1, idx2) -
      q_alpha(k) * stats_vec.stddev().subvec(idx1, idx2));
    mean.col(k) = (stats_vec.mean().subvec(idx1, idx2));
    upper.col(k) = (stats_vec.mean().subvec(idx1, idx2) +
      q_alpha(k) * stats_vec.stddev().subvec(idx1, idx2));
    eigenval_intervals.col(k) = arma::quantile(eval_mat.col(k), alpha_vec_eval);
    eigenval_pve_intervals.col(k) = arma::quantile(eval_pve_mat.col(k), alpha_vec_eval);
  }
  magnitude_interval = arma::quantile(magnitude, alpha_vec_eval);
  arma::mat mean_cov =
    arma::reshape(stats_cov.mean(), basis.n_rows, basis.n_rows);
  return(List::create(Named("lower_eigen", lower),
                      Named("mean_eigen", mean),
                      Named("upper_eigen", upper),
                      Named("eigenval_intervals", eigenval_intervals),
                      Named("eigenval_pve_intervals", eigenval_pve_intervals),
                      Named("surface", mean_cov),
                      Named("magnitude", magnitude_interval),
                      Named("raw_magnitude", magnitude),
                      Named("time", time),
                      Named("pve_mat", eval_pve_mat)));
}

// [[Rcpp::export]]
Rcpp::List get_posterior_covariance_cpp(Rcpp::List mcmc_results, arma::vec zi) {
  Rcpp::List control = mcmc_results["control"];
  Rcpp::List data = mcmc_results["data"];
  Rcpp::List samples = mcmc_results["samples"];
  arma::mat basis = data["basis"];
  arma::field<arma::cube> lambda = samples["lambda"];
  arma::cube beta = samples["beta"];
  arma::vec time = data["time"];
  arma::uword iterations = control["iterations"];
  arma::uword burnin = control["burnin"];
  arma::uword kdim = data["latent_dimension"];
  arma::mat est_cov = arma::mat(basis.n_cols, basis.n_cols);
  arma::mat est_cov_full;
  arma::running_stat_vec<arma::vec> stats;
  for (arma::uword i = burnin; i < iterations; i++) {
    est_cov.zeros();
    for (arma::uword k = 0; k < kdim; k++) {
      est_cov = est_cov + lambda(i).slice(k) * zi * 
        zi.t() * lambda(i).slice(k).t();
    }
    est_cov_full = basis * est_cov * basis.t();
    stats(arma::vectorise(est_cov_full));
  }
  return(Rcpp::List::create(Rcpp::Named("mean", stats.mean()),
                            Rcpp::Named("sd", stats.stddev())));
}
