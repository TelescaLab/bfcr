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
  //delta_cumprod = arma::vec(dat.indices_var.n_elem);
  delta_cumprod = arma::mat(dat.n_smooths_var, dat.kdim);
  blk_diag_phi_delta = arma::cube(dat.basis_dim * dat.d2,
                                  dat.basis_dim * dat.d2,
                                  dat.kdim, arma::fill::zeros);
  //blk_diag_delta_cumprod = arma::mat(dat.indices_var.n_elem,
                //                     dat.kdim, arma::fill::ones);*
}

void Transformations::build_blk_diag_mean(Data& dat, Parameters& pars) {
  arma::uword old_index = 0;
  arma::uword start;
  arma::uword end = -1;
  for (arma::uword d = 0; d < dat.indices_mean.n_elem; d++) {
    if (dat.indices_mean(d) == old_index) {
      blk_diag_mean_penalties.submat(
        start, start, end, end) =
          pars.tau1(d) * dat.penalties_mean(d) + 
          blk_diag_mean_penalties.submat(start, start, end, end);
    } else {
      start = end + 1;
      end = start + dat.penalties_mean(d).n_rows - 1;
      
      blk_diag_mean_penalties.submat(
        start, start, end, end) =
          pars.tau1(d) * dat.penalties_mean(d);
      
    }
    
    old_index = dat.indices_mean(d);
  }
}

void Transformations::build_blk_diag_var(Data& dat,
                                         Parameters& pars) {
  for (arma::uword k = 0; k < dat.kdim; k++) {
    arma::uword old_index = 0;
    arma::uword start;
    arma::uword end = -1;
    for (arma::uword d = 0; d < dat.indices_var.n_elem; d++) {
      if (dat.indices_var(d) == old_index) {
        blk_diag_var_penalties.slice(k).submat(
            start, start, end, end) =
              pars.tau2(d) * dat.penalties_var(d) + 
              blk_diag_var_penalties.slice(k).submat(
                  start, start, end, end);
      } else {
        start = end + 1;
        end = start + dat.penalties_var(d).n_rows - 1;
        blk_diag_var_penalties.slice(k).submat(
            start, start, end, end) =
              pars.tau2(d, k) * dat.penalties_var(d);
        
      }
      
      old_index = dat.indices_var(d);
    }
  }
}
void Transformations::build_blk_diag_phi_delta(Data& dat, Parameters& pars) {
  for (arma::uword i = 0; i < dat.n_smooths_var; i++) {
    delta_cumprod.row(i) = arma::cumprod(pars.delta.row(i));
    for (arma::uword k = 0; k < dat.kdim; k++) {
      blk_diag_phi_delta.slice(k).submat(
          dat.basis_dim * dat.seq_along_start(i), 
          dat.basis_dim * dat.seq_along_start(i),
          dat.basis_dim * (dat.seq_along_end(i) + 1) - 1,
          dat.basis_dim * (dat.seq_along_end(i) + 1) - 1) =
            arma::diagmat(arma::vectorise(
                pars.phi.slice(k).cols(dat.seq_along_start(i),
                               dat.seq_along_end(i)) * delta_cumprod(i, k)));
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
