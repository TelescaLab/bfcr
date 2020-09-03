#include <RcppArmadillo.h>
#include "Utility.h"
#include <progress.hpp>
#include <iostream>
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]
static double const log2pi = std::log(2.0 * M_PI);

class Parameters;
class Data;
class Transformations;
class GaussTransformations;
class McmcObject;


class Parameters {
  public:
    arma::cube lambda;
    arma::mat beta;
    arma::mat eta;
    arma::vec varphi;
    arma::cube phi;
    double alpha;
    double tausq;
    arma::vec psi;
    arma::vec tau1;
    arma::mat tau2;
    arma::mat delta;
    arma::uword iteration = 0;
    Parameters(Data&);
    Parameters() {};
    void update_beta(Data&, Transformations&);
    void update_lambda(Data&, Transformations&);
    void update_eta(Data&, Transformations&);
    void update_varphi(Data&, Transformations&);
    void update_tau1(Data&, Transformations&);
    void update_tau2(Data&, Transformations&);
    void update_phi(Data&, Transformations&);
    void update_delta(Data&, Transformations&);
    void update_psi(Data&, Transformations&);
    void update_tausq(Data&, Transformations&);
    void update_alpha(Data&, Transformations&);
    void write_parameters();
    double varphi_a = .001;
    double varphi_b = .001;
    double tau_a = 1;
    double tau_b = 0.005;
    double phi_a = 2;
    double phi_b = 2;
    double delta_a1 = 2;
    double delta_a2 = 2;
    double delta_b = 1;
    double nu = 5;
    arma::cube beta_container;
    arma::field<arma::cube> lambda_container;
    arma::cube eta_container;
    arma::mat varphi_container;
    arma::mat tau1_container;
    arma::cube tau2_container;
    arma::field<arma::cube> phi_container;
    arma::cube delta_container;
    arma::vec tausq_container;
};

class GaussTransformations {
  public:
    arma::mat chol;
    arma::vec w, mu, z, v, result;
    GaussTransformations(arma::uword dim);
    GaussTransformations() {};
    arma::vec bayes_reg(arma::mat& precision, arma::vec& mu);
};
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
    arma::vec delta_cumprod;
    arma::cube blk_diag_phi_delta;
    arma::mat blk_diag_delta_cumprod;
    arma::vec squared_diff;
};

class Data {
  public:
    arma::uword d1, d2, basis_dim, kdim, n;
    arma::uvec missing, missing_sub, missing_time;
    arma::mat response, design_mean;
    arma::mat design_var, basis;
    arma::field<arma::mat> penalties_mean;
    arma::field<arma::mat> penalties_var;
    arma::uvec indices_mean, indices_var;
    arma::uword iter, thin;
    Data(arma::mat&, arma::mat&, arma::mat&, arma::mat&, 
         arma::field<arma::mat>&,arma::field<arma::mat>&,
         arma::uvec&, arma::uvec&, arma::uword, arma::uword, 
         arma::uword);
    Data() {};
};

Parameters::Parameters(Data& dat) {
  lambda = arma::cube(dat.basis_dim, dat.d2, dat.kdim, arma::fill::randn);
  beta = arma::mat(dat.basis_dim, dat.d1, arma::fill::randn);
  eta = arma::mat(dat.n, dat.kdim, arma::fill::randn);
  varphi = 100 * arma::vec(dat.n, arma::fill::ones);
  psi = arma::vec(dat.n, arma::fill::ones);
  alpha = 1;
  tausq = 1;
  tau1 = arma::vec(dat.penalties_mean.n_elem, arma::fill::zeros);
  tau2 = arma::mat(dat.penalties_var.n_elem, dat.kdim, arma::fill::zeros);
  phi = arma::cube(dat.basis_dim, dat.d2, dat.kdim, arma::fill::ones);
  delta = arma::mat(dat.penalties_var.n_elem, dat.kdim, arma::fill::ones);
  beta_container = arma::cube(dat.basis_dim, dat.d1, dat.iter);
  lambda_container = arma::field<arma::cube>(dat.iter);
  eta_container = arma::cube(dat.n, dat.kdim, dat.iter);
  varphi_container = arma::mat(dat.n, dat.iter);
  tau1_container = arma::mat(dat.indices_mean.n_elem, dat.iter);
  tau2_container = arma::cube(dat.indices_var.n_elem, dat.kdim, dat.iter);
  phi_container = arma::field<arma::cube>(dat.iter);
  delta_container = arma::cube(dat.indices_var.n_elem, dat.kdim, dat.iter);
  tausq_container = arma::vec(dat.iter);
}

void Parameters::update_beta(Data& dat, Transformations& transf) {
  transf.beta_g = 
    arma::vectorise((transf.bty - transf.btb * transf.fit_lambda) *
    arma::diagmat(varphi) * dat.design_mean);
  transf.beta_precision = arma::kron(dat.design_mean.t() *
    arma::diagmat(varphi) * dat.design_mean, transf.btb) +
    transf.blk_diag_mean_penalties;
  if(!transf.beta_precision.is_symmetric()) {
    transf.beta_precision = arma::symmatu(transf.beta_precision);
  }
  transf.beta_result = 
    transf.beta_gauss.bayes_reg(transf.beta_precision, transf.beta_g);
  beta = arma::reshape(transf.beta_result, dat.basis_dim, dat.d1);
  transf.fit_beta = beta * dat.design_mean.t();
}

void Parameters::update_lambda(Data& dat, Transformations& transf) {
  for (arma::uword k = 0; k < dat.kdim; k++) {
    transf.fit_lambda_removed = transf.fit_beta + 
      transf.fit_lambda - (lambda.slice(k) *
      dat.design_var.t() * arma::diagmat(eta.col(k)));
    
    
    transf.lambda_precision = 
      arma::kron(dat.design_var.t() * arma::diagmat(eta.col(k)) * 
      arma::diagmat(varphi) * 
      (dat.design_var.t() * arma::diagmat(eta.col(k))).t(),
      transf.btb) + transf.blk_diag_var_penalties.slice(k) +
        arma::diagmat(transf.blk_diag_phi_delta.slice(k));
    transf.lambda_g = arma::vectorise((transf.bty - 
      transf.btb * transf.fit_lambda_removed) * arma::diagmat(varphi) *
      arma::diagmat(eta.col(k)) * dat.design_var);
    if(!transf.lambda_precision.is_symmetric()){
      transf.lambda_precision = arma::symmatu(transf.lambda_precision);
    }
    transf.lambda_result = 
      transf.lambda_gauss.bayes_reg(transf.lambda_precision, transf.lambda_g);
    /*
    transf.lambda_chol = arma::chol(transf.lambda_precision, "lower");
    transf.lambda_w = solve(arma::trimatl(transf.lambda_chol), transf.lambda_g);
    transf.lambda_mu = solve(arma::trimatu(transf.lambda_chol.t()),
                             transf.lambda_w);
    transf.lambda_z = arma::randn<arma::vec>(dat.basis_dim * dat.d2);
    transf.lambda_v = arma::solve(arma::trimatu(transf.lambda_chol.t()),
                                  transf.lambda_z);
    transf.lambda_result = transf.lambda_mu + transf.lambda_v;
     */
    transf.lambda_old = lambda.slice(k);
    lambda.slice(k) = arma::reshape(transf.lambda_result, dat.basis_dim, dat.d2);
    transf.fit_lambda = transf.fit_lambda + 
      (lambda.slice(k) - transf.lambda_old) * dat.design_var.t() *
      arma::diagmat(eta.col(k));
    //transf.update_fit(dat, *this);
    
    
  }
}

void Parameters::update_eta(Data& dat, Transformations& transf) {

  for (arma::uword subject = 0; subject < dat.n; subject++) {
    for (arma::uword k = 0; k < dat.kdim; k++) {
      transf.eta_transf.col(k) = lambda.slice(k) * dat.design_var.row(subject).t();
    }
    transf.eta_g = varphi(subject) * transf.eta_transf.t() * 
      (transf.bty.col(subject) - transf.btb * transf.fit_beta.col(subject));

    transf.eta_precision = varphi(subject) * transf.eta_transf.t() * 
      transf.btb * transf.eta_transf + 
      arma::eye<arma::mat>(dat.kdim, dat.kdim);
    
    if(!transf.eta_precision.is_symmetric()){
      transf.eta_precision = arma::symmatu(transf.eta_precision);
    }
    transf.eta_result = 
      transf.eta_gauss.bayes_reg(transf.eta_precision, transf.eta_g);
    eta.row(subject) = (transf.eta_result).t();

  }
  transf.fit_lambda.zeros();
  for (arma::uword k = 0; k < dat.kdim; k++) {
    transf.fit_lambda = transf.fit_lambda +
      lambda.slice(k) * 
      dat.design_var.t() * 
      arma::diagmat(eta.col(k));
  }
  
}

void Parameters::update_tau1(Data& dat, Transformations& transf) {
  double update_a = 0, update_b = 0;
  arma::uword start = 0;
  arma::uword end = static_cast<double>(dat.penalties_mean(0).n_rows) /
    static_cast<double>(dat.basis_dim) - 1;
  arma::uword num_field_elements = dat.penalties_mean.n_elem;
  arma::uword old_index = 1;
  
  for(arma::uword i = 0; i < num_field_elements; i++){
    
    if(dat.indices_mean(i) != old_index){
      start = end + 1;
      end = end + dat.penalties_mean(i).n_rows / dat.basis_dim;
    }
    
    update_a = dat.penalties_mean(i).n_rows;
    update_b = 1.0 / 2.0 *
      arma::as_scalar(arma::vectorise(beta.cols(start, end)).t() *
      dat.penalties_mean(i) *
      arma::vectorise(beta.cols(start, end)));
    tau1(i) = R::rgamma(tau_a + update_a / 2.0, 1.0 / (tau_b + update_b));
  }
  transf.build_blk_diag_mean(dat, *this);
}

void Parameters::update_tau2(Data& dat, Transformations& transf) {
  double update_a = 0, update_b = 0;
  arma::uword start = 0;
  arma::uword end = static_cast<double>(dat.penalties_var(0).n_rows) /
    static_cast<double>(dat.basis_dim) - 1;
  arma::uword num_field_elements = dat.penalties_var.n_elem;
  arma::uword old_index = 1;
  
  for(arma::uword i = 0; i < num_field_elements; i++){
    
    if (dat.indices_var(i) != old_index) {
      start = end + 1;
      end = end + dat.penalties_var(i).n_rows / dat.basis_dim;
    }
    
    for (arma::uword k = 0; k < dat.kdim; k++) {
      update_a = dat.penalties_var(i).n_rows;
      
      update_b = 
        arma::as_scalar(arma::vectorise(lambda.slice(k).cols(start, end)).t() *
        dat.penalties_var(i) * arma::vectorise(lambda.slice(k).cols(start, end)));
      tau2(i, k) = R::rgamma(tau_a + update_a, 1.0 /
        tau_b + update_b);
    }
    update_a = 0;
    update_b = 0;
    old_index = dat.indices_var(i);
    
  }
  transf.build_blk_diag_var(dat, *this);
}
void Parameters::update_varphi(Data& dat, Transformations& transf) {
  transf.fit = transf.fit_beta + transf.fit_lambda;
  double my_sum = arma::accu(arma::square(dat.response - 
                             arma::trans(dat.basis * (transf.fit))));
  double p = R::rgamma(varphi_a + dat.response.n_elem/2, 1.0 / 
                       (varphi_b + 1.0/2.0 * my_sum));
  varphi.fill(p);
}


void Parameters::update_phi(Data& dat, Transformations& transf) {
  for (arma::uword i = 0; i < dat.basis_dim; i++) {
    for (arma::uword j = 0; j < dat.d2; j++) {
      for (arma::uword k = 0; k < dat.kdim; k++) {
        phi(i, j, k) = 
          R::rgamma(phi_a + .5, 1.0 / (phi_b + ::pow(lambda(i, j, k), 2)));
      }
    }
  }
}

void Parameters::update_delta(Data& dat, Transformations& transf) {
  double update_a = 0, update_b = 0;
  arma::uword start = 0;
  arma::uword end = static_cast<double>(dat.penalties_var(0).n_rows) /
    static_cast<double>(dat.basis_dim) - 1;
  arma::uword num_field_elements = dat.penalties_var.n_elem;
  arma::uword old_index = 1;
  for(arma::uword i = 0; i < num_field_elements; i++){
    if (dat.indices_var(i) != old_index) {
      start = end + 1;
      end = end + dat.penalties_var(i).n_rows / dat.basis_dim;
    }
    for (arma::uword k = 0; k < dat.kdim; k++) {
      transf.phi_lambda_sum(i, k) = arma::as_scalar(
        arma::accu(arma::square(lambda.slice(k).cols(start, end)) %
          phi.slice(k).cols(start, end)));
    }
    old_index = dat.indices_var(i);
  }
  start = 0;
  end = static_cast<double>(dat.penalties_var(0).n_rows) /
    static_cast<double>(dat.basis_dim) - 1;
  old_index = 1;
  for(arma::uword i = 0; i < num_field_elements; i++){
    if (dat.indices_var(i) != old_index) {
      start = end + 1;
      end = end + dat.penalties_var(i).n_rows / dat.basis_dim;
    }
    
    for (arma::uword k = 0; k < dat.kdim; k++) {
      
      transf.delta_cumprod = delta.row(i).t();
      transf.delta_cumprod(k) = 1;
      transf.delta_cumprod = arma::cumprod(transf.delta_cumprod);
      for (arma::uword kp = 0; kp < k; kp++) {
        transf.delta_cumprod(kp) = 0;
      }

      update_b = delta_b + .5 * arma::as_scalar(arma::accu(
        transf.delta_cumprod % transf.phi_lambda_sum.row(i).t()));
      if (k == 0) {
        update_a = delta_a1 + dat.basis_dim * (end - start + 1) * (dat.kdim - k) / 2;
        
        delta(i, k) = R::rgamma(update_a, 1.0 / update_b);
      } else {
        update_a = delta_a2 + dat.basis_dim * (end - start + 1) * (dat.kdim - k) / 2;
        
        delta(i, k) = R::rgamma(update_a, 1.0 / update_b);
        
      }
    }
    old_index = dat.indices_var(i);
  }
  transf.build_blk_diag_phi_delta(dat, *this);
  
}
void Parameters::update_psi(Data& dat, Transformations& transf){
  transf.fit = transf.fit_beta + transf.fit_lambda;
  transf.squared_diff = arma::sum(arma::square(dat.response - 
    arma::trans(dat.basis * (transf.fit))), 1);
  double update_a = 0;
  double update_b = 0;
  double psi_a = nu / 2;
  double psi_b = nu * tausq / 2;
  update_a = dat.basis.n_rows / 2;
  for (arma::uword i = 0; i < dat.n; i++) {
    update_b = alpha * transf.squared_diff(i) / 2;
    psi(i) = R::rgamma(update_a + psi_a, 1 / update_b + psi_b);
  }
}

void Parameters::update_tausq(Data& dat, Transformations& transf) {
  tausq = R::rgamma(dat.n * nu / 2, 1.0 / (nu / 2 * arma::accu(psi)));
}

void Parameters::update_alpha(Data& dat, Transformations& transf) {
  double update_a = dat.response.n_elem / 2;
  double update_b = arma::accu(transf.squared_diff % psi) / 2;
  //Rcpp::Rcout << "MEAN: " << update_a / update_b << std::endl << 
    //"UPDATE_B: " << update_b << std::endl << "SQUARED DIFF: " << transf.squared_diff(0) << std::endl;
  alpha = R::rgamma(update_a, 1 / update_b);
}

void Parameters::write_parameters() {
  beta_container.slice(iteration) = beta;
  lambda_container(iteration) = lambda;
  eta_container.slice(iteration) = eta;
  varphi_container.col(iteration) = varphi;
  tau1_container.col(iteration) = tau1;
  tau2_container.slice(iteration) = tau2;
  phi_container(iteration) = phi;
  delta_container.slice(iteration) = delta;
  iteration = iteration + 1;
}

GaussTransformations::GaussTransformations(arma::uword dim) {
  arma::mat chol = arma::mat(dim, dim);
  arma::vec w = arma::vec(dim), mu = arma::vec(dim), v = arma::vec(dim),
    z = arma::vec(dim), result = arma::vec(dim);
}

arma::vec GaussTransformations::bayes_reg(arma::mat& precision, arma::vec& g) {
  chol = arma::chol(precision, "lower");
  w = solve(arma::trimatl(chol), g);
  mu = solve(arma::trimatu(chol.t()), w);
  z = arma::randn<arma::vec>(g.n_elem);
  v = arma::solve(arma::trimatu(chol.t()), z);
  result = mu + v;
  return(result);
}

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
  phi_lambda_sum = arma::mat(dat.indices_var.n_elem,
                             dat.kdim, arma::fill::zeros);
  delta_cumprod = arma::vec(dat.indices_var.n_elem);
  blk_diag_phi_delta = arma::cube(dat.basis_dim * dat.d2,
                                  dat.basis_dim * dat.d2,
                                  dat.kdim, arma::fill::zeros);
  blk_diag_delta_cumprod = arma::mat(dat.indices_var.n_elem,
                                     dat.kdim, arma::fill::ones);
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
        blk_diag_phi_delta.slice(k).submat(
            start_phi_delta, start_phi_delta,
            end_phi_delta, end_phi_delta) =
              blk_diag_delta_cumprod(d, k) * 
              arma::diagmat(arma::vectorise(pars.phi.slice(k).cols(start, end)));
      } 
      old_index = dat.indices_mean(d);
    }
    
  }
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





class McmcObject {
public:
  arma::uword iter, thin;
  Data dat;
  Transformations transf;
  Parameters pars;
  std::string var;
  McmcObject(Data& dat, Transformations& transf, Parameters& pars,
             std::string var) :
    dat(dat), transf(transf), pars(pars), var(var){}
  void sample_parameters();
  Rcpp::List get_samples();
};

void McmcObject::sample_parameters() {
  Progress progress_bar(dat.iter, true);
  for (arma::uword i = 0; i < dat.iter; i++) {
    for (arma::uword j = 0; j < dat.thin; j++) {
      if (Progress::check_abort() ){
        Rcpp::Rcout << "MCMC aborted" << std::endl;
        
        goto stop;
      }
      pars.update_beta(dat, transf);
      pars.update_lambda(dat, transf); 
      pars.update_eta(dat, transf);
      pars.update_tau1(dat, transf);
      pars.update_tau2(dat, transf);
      pars.update_phi(dat, transf);
      pars.update_delta(dat, transf);
      if (var == "pooled") pars.update_varphi(dat, transf); // uses updated beta and lambda
      if (var == "unequal") {
        pars.update_psi(dat, transf);
        pars.update_tausq(dat, transf);
        pars.varphi = pars.psi;
      }
      transf.complete_response(dat, pars);
    }
    progress_bar.increment();
    pars.write_parameters();
  }
  stop:
    NULL;
}

Rcpp::List McmcObject::get_samples() {
  return(Rcpp::List::create(Rcpp::Named("lambda", pars.lambda_container),
                            Rcpp::Named("beta", pars.beta_container),
                            Rcpp::Named("eta", pars.eta_container),
                            Rcpp::Named("varphi", pars.varphi_container),
                            Rcpp::Named("tau1", pars.tau1_container),
                            Rcpp::Named("tau2", pars.tau2_container),
                            Rcpp::Named("phi", pars.phi_container),
                            Rcpp::Named("delta", pars.delta_container),
                            Rcpp::Named("phi_delta", transf.blk_diag_phi_delta),
                            Rcpp::Named("var_penalty", transf.blk_diag_var_penalties)));
}


 
// [[Rcpp::export]]
Rcpp::List mymain(arma::mat response, arma::mat design_mean,
                  arma::mat design_var, arma::mat basis,
                  arma::field<arma::mat> penalties_mean,
                  arma::field<arma::mat> penalties_var,
                  arma::uvec indices_mean, arma::uvec indices_var,
                  arma::uword kdim, arma::uword iter,
                  arma::uword thin=1, std::string var="unequal"){
  Data dat(response, design_mean,
           design_var, basis,
           penalties_mean, penalties_var,
           indices_mean, indices_var, kdim,
           iter, thin);
  
  Parameters pars(dat);
  
  Transformations transf(dat, pars);
  McmcObject mcmc(dat, transf, pars, var);
  mcmc.sample_parameters();
  Rcpp::List output = mcmc.get_samples();
  return(output);
}


class Sampler
{ 
public:
  Data dat;
  Parameters pars;
  Transformations transf;
  // parameterized constructor
  Sampler(Data& dat, Parameters& pars, Transformations& transf) :
    dat(dat), pars(pars), transf(transf)
  { 
    Rcpp::Rcout << "Base Parameterized Constructor\n";
  }
  virtual void print() = 0;
  virtual void sample_parameters() = 0;
  virtual void write_parameters() = 0;
  virtual Rcpp::List get_samples() = 0;
  ~Sampler(){}
};

class SamplerPooled : public Sampler
{ 
public:
  // parameterized constructor
  SamplerPooled(Data& dat, Parameters& pars, Transformations& transf) :
  Sampler(dat, pars, transf)
  { 
    Rcpp::Rcout << "Derived pooled\n";
  }
  void print() {Rcpp::Rcout << "printing pooled\n";}
  void sample_parameters();
  void write_parameters();
  Rcpp::List get_samples();
  ~SamplerPooled(){};
};

class SamplerUnequal : public Sampler{
  public:
    SamplerUnequal(Data& dat, Parameters& pars, Transformations& transf) :
    Sampler(dat, pars, transf) {
      Rcpp::Rcout << "Derived unequal\n";
    }
    void print() {Rcpp::Rcout << "printing unequal\n";}
    void sample_parameters();
    void write_parameters();
    Rcpp::List get_samples();
    ~SamplerUnequal(){}
};



class SamplerFactory {
public:
  static Sampler *new_mcmc(std::string type, Data& dat, Parameters& pars, Transformations& transf) {
    if(type == "pooled") return new SamplerPooled(dat, pars, transf);
    if(type == "unequal") return new SamplerUnequal(dat, pars, transf);
    return nullptr;
  }
};

void SamplerUnequal::sample_parameters() {
  Progress progress_bar(dat.iter, true);
  for (arma::uword i = 0; i < dat.iter; i++) {
    for (arma::uword j = 0; j < dat.thin; j++) {
      if (Progress::check_abort() ){
        Rcpp::Rcout << "MCMC aborted" << std::endl;
        
        goto stop;
      }
      pars.update_beta(dat, transf);
      pars.update_lambda(dat, transf); 
      pars.update_eta(dat, transf);
      pars.update_tau1(dat, transf);
      pars.update_tau2(dat, transf);
      pars.update_phi(dat, transf);
      pars.update_delta(dat, transf);
      pars.update_psi(dat, transf);
      pars.update_tausq(dat, transf);
      pars.varphi = pars.psi;
      transf.complete_response(dat, pars);
    }
    progress_bar.increment();
    write_parameters();
  }
  stop:
    NULL;
}

void SamplerPooled::sample_parameters() {
  Progress progress_bar(dat.iter, true);
  for (arma::uword i = 0; i < dat.iter; i++) {
    for (arma::uword j = 0; j < dat.thin; j++) {
      if (Progress::check_abort() ){
        Rcpp::Rcout << "MCMC aborted" << std::endl;
        
        goto stop;
      }
      pars.update_beta(dat, transf);
      pars.update_lambda(dat, transf); 
      pars.update_eta(dat, transf);
      pars.update_tau1(dat, transf);
      pars.update_tau2(dat, transf);
      pars.update_phi(dat, transf);
      pars.update_delta(dat, transf);
      pars.update_varphi(dat, transf); 
      transf.complete_response(dat, pars);
    }
    progress_bar.increment();
    pars.write_parameters();
  }
  stop:
    NULL;
}

void SamplerUnequal::write_parameters() {
  pars.beta_container.slice(pars.iteration) = pars.beta;
  pars.lambda_container(pars.iteration) = pars.lambda;
  pars.eta_container.slice(pars.iteration) = pars.eta;
  pars.varphi_container.col(pars.iteration) = pars.varphi;
  pars.tau1_container.col(pars.iteration) = pars.tau1;
  pars.tau2_container.slice(pars.iteration) = pars.tau2;
  pars.delta_container.slice(pars.iteration) = pars.delta;
  pars.tausq_container(pars.iteration) = pars.tausq;
  pars.iteration++;
}

void SamplerPooled::write_parameters() {
  pars.beta_container.slice(pars.iteration) = pars.beta;
  pars.lambda_container(pars.iteration) = pars.lambda;
  pars.eta_container.slice(pars.iteration) = pars.eta;
  pars.varphi_container.col(pars.iteration) = pars.varphi;
  pars.tau1_container.col(pars.iteration) = pars.tau1;
  pars.tau2_container.slice(pars.iteration) = pars.tau2;
  pars.delta_container.slice(pars.iteration) = pars.delta;
  pars.iteration++;
}

Rcpp::List SamplerUnequal::get_samples() {
  return(Rcpp::List::create(Rcpp::Named("lambda", pars.lambda_container),
                            Rcpp::Named("beta", pars.beta_container),
                            Rcpp::Named("eta", pars.eta_container),
                            Rcpp::Named("varphi", pars.varphi_container),
                            Rcpp::Named("tau1", pars.tau1_container),
                            Rcpp::Named("tau2", pars.tau2_container),
                            Rcpp::Named("phi", pars.phi_container),
                            Rcpp::Named("delta", pars.delta_container),
                            Rcpp::Named("tausq", pars.tausq_container)));
}

Rcpp::List SamplerPooled::get_samples() {
  return(Rcpp::List::create(Rcpp::Named("lambda", pars.lambda_container),
                            Rcpp::Named("beta", pars.beta_container),
                            Rcpp::Named("eta", pars.eta_container),
                            Rcpp::Named("varphi", pars.varphi_container),
                            Rcpp::Named("tau1", pars.tau1_container),
                            Rcpp::Named("tau2", pars.tau2_container),
                            Rcpp::Named("phi", pars.phi_container),
                            Rcpp::Named("delta", pars.delta_container)));
}


// [[Rcpp::export]]
Rcpp::List run_mcmc(arma::mat response, arma::mat design_mean,
               arma::mat design_var, arma::mat basis,
               arma::field<arma::mat> penalties_mean,
               arma::field<arma::mat> penalties_var,
               arma::uvec indices_mean, arma::uvec indices_var,
               arma::uword kdim, arma::uword iter,
               arma::uword thin=1, std::string var="unequal") {
  
  Data dat(response, design_mean,
           design_var, basis,
           penalties_mean, penalties_var,
           indices_mean, indices_var, kdim,
           iter, thin);
  Parameters pars(dat);
  Transformations transf(dat, pars);
  Sampler* mysampler = SamplerFactory::new_mcmc(var, dat, pars, transf);
  mysampler->sample_parameters();
  return(mysampler->get_samples());
}