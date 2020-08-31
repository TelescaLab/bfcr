#include <RcppArmadillo.h>
#include "Utility.h"

static double const log2pi = std::log(2.0 * M_PI);

class Parameters;
class Data;
class Transformations;

class Parameters {
  public:
    arma::cube lambda;
    arma::mat beta;
    arma::mat eta;
    arma::vec varphi;
    arma::vec phi;
    double alpha;
    double nu;
    arma::vec tau1;
    arma::vec tau2;
    Parameters(Data&);
    void update_beta(Data&, Transformations&);
    void update_lambda(Data&, Transformations&);
    void update_eta(Data&, Transformations&);
    void update_prec(Data&, Transformations&);
    void update_tau1(Data&, Transformations&);
    void update_tau2(Data&, Transformations&);
    
};


class Transformations {
  public:
    arma::mat beta_chol, beta_precision;
    arma::vec beta_g, beta_w, beta_mu, beta_z, beta_v, beta_result;
    arma::mat btb, bty, fit, lambda_tilde;
    arma::mat blk_diag_mean_penalties;
    arma::cube blk_diag_var_penalties;
    void build_blk_diag_mean(Data&, Parameters&);
    void build_blk_diag_var(Data&, Parameters&);
    Transformations(Data&, Parameters&);
    arma::mat fit_beta_removed;
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
    Data(arma::mat&, arma::mat&, arma::mat&, arma::mat&, 
         arma::field<arma::mat>&,arma::field<arma::mat>&,
         arma::uvec&, arma::uvec&, arma::uword);
};

Parameters::Parameters(Data& dat) {
  lambda = arma::cube(dat.basis_dim, dat.d2, dat.kdim, arma::fill::randn);
  beta = arma::mat(dat.basis_dim, dat.d1, arma::fill::randn);
  eta = arma::mat(dat.n, dat.kdim, arma::fill::randn);
  varphi = arma::vec(dat.n);
  phi = arma::vec(dat.n);
  double alpha = 1;
  nu = 1;
  
}

void Parameters::update_beta(Data& dat, Transformations& transf) {
  transf.fit_beta_removed = transf.fit.t() - 
    dat.basis * beta * dat.design_mean.t();
  transf.beta_g = 
    arma::vectorise((transf.bty - dat.basis.t() * transf.fit_beta_removed) *
    arma::diagmat(varphi) * dat.design_mean);
                      
  transf.fit_beta_removed.zeros();
}
Transformations::Transformations(Data& dat, Parameters& pars) {
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
  blk_diag_var_penalties = arma::cube(
    dat.basis_dim * dat.d2,
    dat.basis_dim * dat.d2,
    dat.kdim,
    arma::fill::zeros);
  

}
Data::Data(arma::mat& response, arma::mat& design_mean,
           arma::mat& design_var, arma::mat& basis,
           arma::field<arma::mat>& penalties_mean, 
           arma::field<arma::mat>& penalties_var,
           arma::uvec& indices_mean, arma::uvec& indices_var,
           arma::uword kdim) {
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
  missing = arma::find_nonfinite(response);
  missing_sub = armadillo_modulus3(missing, response.n_rows);
  missing_time = arma::floor(missing / response.n_rows);
  response.elem(missing).fill(0);
}
// [[Rcpp::export]]
void mymain(arma::mat response, arma::mat design_mean,
            arma::mat design_var, arma::mat basis,
            arma::field<arma::mat> penalties_mean,
            arma::field<arma::mat> penalties_var,
            arma::uvec indices_mean, arma::uvec indices_var,
            arma::uword kdim, arma::uword iter,
            arma::uword thin=1){
  Data dat(response, design_mean,
              design_var, basis,
              penalties_mean, penalties_var,
              indices_mean, indices_var, kdim);
  Parameters pars(dat);
  Transformations transf(dat, pars);
  pars.update_beta(dat, transf);
}