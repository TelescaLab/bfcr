#include "Transformations.h"

Transformations::Transformations(Data& dat, Parameters& pars) {
  beta_gauss = GaussTransformations(dat.basis_dim * dat.d1);
  lambda_gauss = GaussTransformations(dat.basis_dim * dat.d2);
  eta_gauss = GaussTransformations(dat.kdim);
  fit_beta_removed = arma::mat(dat.basis_dim, dat.n, arma::fill::zeros);
  fit_lambda_removed = arma::mat(dat.basis_dim, dat.n, arma::fill::zeros);
  btb = dat.basis.t() * dat.basis;
  bty = dat.basis.t() * dat.response.t();
  fit = dat.design_mean * pars.beta.t() * dat.basis.t();
  for (arma::uword k = 0; k < dat.kdim; k++) {
    fit = fit + arma::diagmat(pars.eta.col(k)) * 
      dat.design_var * pars.lambda.slice(k).t() * dat.basis.t();
  }
  lambda_tilde = arma::mat(dat.basis_dim, dat.kdim);
  blk_diag_mean_penalties = arma::mat(
    dat.basis_dim * dat.d1,
    dat.basis_dim * dat.d1,
    arma::fill::zeros);
  build_blk_diag_mean(dat, pars);
  blk_diag_var_penalties = arma::cube(
    dat.basis_dim * dat.d2,
    dat.basis_dim * dat.d2,
    dat.kdim,
    arma::fill::zeros);
  eta_transf = arma::mat(dat.basis_dim, dat.kdim);
  fit_lambda = arma::mat(dat.basis_dim, dat.n, arma::fill::zeros);
  for (arma::uword k = 0; k < dat.kdim; k++){
    fit_lambda = fit_lambda + 
      pars.lambda.slice(k) * 
      dat.design_var.t() * 
      arma::diagmat(pars.eta.col(k));
  }
  fit_beta = pars.beta * dat.design_mean.t();
  lambda_old = arma::mat(dat.basis_dim, dat.d2);
  phi_lambda_sum = arma::mat(dat.n_smooths_var,
                             dat.kdim, arma::fill::zeros);
  // delta_cumprod = arma::mat(dat.n_smooths_var, dat.kdim, arma::fill::ones);
  delta_cumprod = arma::mat(dat.n_smooths_var, dat.kdim, arma::fill::zeros);
  for (arma::uword i = 0; i < dat.n_smooths_var; i++) {
    delta_cumprod.row(i) = arma::cumprod(pars.delta.row(i));
  }
  blk_diag_phi_delta = arma::cube(dat.basis_dim * dat.d2,
                                  dat.basis_dim * dat.d2,
                                  dat.kdim, arma::fill::zeros);
  psi_mat = arma::mat(dat.basis_dim, dat.basis_dim);
  for(arma::uword j = 0; j < dat.basis_dim; j++){
    for(arma::uword i = 0; i < dat.basis_dim; i++){
      psi_mat(i, j) = 
        arma::as_scalar(
          arma::trapz(dat.time, dat.basis.col(i) % dat.basis.col(j)));
    }
  }
  lin_constr = arma::mat(dat.basis_dim * dat.design_var.n_cols, dat.kdim - 1);
}

void Transformations::build_blk_diag_mean(Data& dat, Parameters& pars) {
  blk_diag_mean_penalties.zeros();
  for (arma::uword i = 0; i < dat.indices_mean.n_elem; i++) {
      blk_diag_mean_penalties.submat(
        dat.basis_dim * dat.seq_along_tau1(i, 0),
        dat.basis_dim * dat.seq_along_tau1(i, 0),
        dat.basis_dim * (dat.seq_along_tau1(i, 1) + 1) - 1,
        dat.basis_dim * (dat.seq_along_tau1(i, 1) + 1) - 1) =
          pars.tau1(i) * dat.penalties_mean(i);
  }
}

void Transformations::build_blk_diag_var(Data& dat,
                                         Parameters& pars) {
  for (arma::uword k = 0; k < dat.kdim; k++) {
    blk_diag_var_penalties.slice(k).zeros();
    for (arma::uword i = 0; i < dat.indices_var.n_elem; i++) {
      blk_diag_var_penalties.slice(k).submat(
          dat.basis_dim * dat.seq_along_tau2(i, 0),
          dat.basis_dim * dat.seq_along_tau2(i, 0),
          dat.basis_dim * (dat.seq_along_tau2(i, 1) + 1) - 1,
          dat.basis_dim * (dat.seq_along_tau2(i, 1) + 1) - 1) =
            pars.tau2(i, k) * dat.penalties_var(i);
    }
  }
}

void Transformations::build_blk_diag_phi_delta(Data& dat, Parameters& pars) {
  for (arma::uword i = 0; i < dat.n_smooths_var; i++) {
    // Rcpp::Rcout << "START: " << dat.seq_along_start_delta(i) << "\n";
    // Rcpp::Rcout << "END: " << dat.seq_along_end_delta(i) << "\n";
    
    for (arma::uword k = 0; k < dat.kdim; k++) {
      blk_diag_phi_delta.slice(k).submat(
          dat.basis_dim * dat.seq_along_start_delta(i), 
          dat.basis_dim * dat.seq_along_start_delta(i),
          dat.basis_dim * (dat.seq_along_end_delta(i) + 1) - 1,
          dat.basis_dim * (dat.seq_along_end_delta(i) + 1) - 1) =
            arma::diagmat(arma::vectorise(
                pars.phi.slice(k).cols(dat.seq_along_start_delta(i),
                               dat.seq_along_end_delta(i)) * delta_cumprod(i, k)));
    }
  }
  /*
  int old_index = 0;
  int start = 0;
  int end = -1;
  int start_phi_delta = 0;
  int end_phi_delta = -1;
  for (arma::uword d = 0; d < dat.indices_var.n_elem; d++) {
    if (dat.indices_var(d) != old_index) {
      start = end + 1;
      end = end_phi_delta + dat.penalties_var(d).n_rows / dat.basis_dim;
      start_phi_delta = end_phi_delta + 1;
      end_phi_delta = end_phi_delta + dat.penalties_var(d).n_rows;
      blk_diag_delta_cumprod.row(d) = arma::cumprod(pars.delta.row(d));
      for (arma::uword k = 0; k < dat.kdim; k++) {
        Rcpp::Rcout << "start: " << start << " end: " << end << std::endl;
        blk_diag_phi_delta.slice(k).submat(
            start_phi_delta, start_phi_delta,
            end_phi_delta, end_phi_delta) =
              blk_diag_delta_cumprod(d, k) * 
              arma::diagmat(arma::vectorise(pars.phi.slice(k).cols(start, end)));
      } 
      old_index = dat.indices_mean(d);
    }
    
  }*/
}

void Transformations::complete_response(Data& dat, Parameters& pars) {
  for(arma::uword i = 0; i < dat.missing_sub.n_elem; i++){
    dat.response(dat.missing_sub(i), dat.missing_time(i)) = 
      arma::as_scalar(dat.basis.row(dat.missing_time(i)) * (fit_beta.col(i) +
      fit_lambda.col(i))) + 
      R::rnorm(0, std::pow(pars.varphi(dat.missing_sub(i)), -.5));
    bty.col(dat.missing_sub(i)) = dat.basis.t() * 
      dat.response.row(dat.missing_sub(i)).t();
  }
}
