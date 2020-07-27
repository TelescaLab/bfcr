#include <RcppArmadillo.h>
#include "Utility.h"

class Parameters;
class Data;
class Transformations;
class Parameters
{

  public:
    arma::cube exp_lambda;
    arma::cube var_lambda;
    arma::mat exp_beta;
    arma::mat var_beta;
    arma::vec exp_prec;
    arma::vec var_prec;
    arma::vec exp_tau1;
    arma::vec var_tau1;
    arma::vec exp_tau2;
    arma::vec var_tau2;
    arma::mat exp_eta;
    arma::cube var_eta;
    arma::cube prec_eta;
    Parameters(Data&);
    void update_lambda(Data&, Transformations&);
    void update_eta(Data&, Transformations&);
    void update_beta(Data&, Transformations&);
    void update_prec(Data&, Transformations&);
};
class Data
{
  public:
    arma::uword D1, D2, p, K, n, n_unique_mean, n_unique_var;
    arma::mat Y, X, Z, B;
    arma::uvec missing, missing_sub, missing_time;
    arma::field<arma::mat> MeanPenalties;
    arma::field<arma::mat> VarPenalties;
    arma::uvec MeanIndices;
    arma::uvec VarIndices;
    Data(arma::mat&, arma::mat&, arma::mat&, arma::mat&, arma::uword,
         arma::field<arma::mat>&,
         arma::field<arma::mat>&,
         arma::uvec&,
         arma::uvec&);
    void impute_response(Transformations&);
};

class Transformations
{
  public:
    arma::cube C_lambda, mychol_lambda, partialZ_eta, ZetaZ;
    arma::mat g_lambda, w_lambda, mu_lambda, z_lambda, v_lambda, result_lambda;
    arma::mat C_beta, mychol_beta;
    arma::vec g_beta, w_beta, mu_beta, z_beta, v_beta, result_beta;
    arma::mat BtY, BtB, Y, fit, mean_only, zzt, lambda_tilde;
    arma::uvec seq_mean, seq_var;
    Transformations(Data&);
    void initialize_transformations(Data&, Parameters&);
};
Parameters::Parameters(Data& dat) {
  exp_lambda = arma::cube(dat.p, dat.D2, dat.K, arma::fill::randn);
  var_lambda = arma::cube(dat.p * dat.D2,
                          dat.p * dat.D2, dat.K, arma::fill::ones);
  exp_beta = arma::mat(dat.p, dat.D1, arma::fill::randn);
  var_beta = arma::mat(dat.p * dat.D1, dat.p * dat.D1, arma::fill::randn);
  exp_prec = arma::vec(dat.n, arma::fill::ones);
  var_prec = arma::vec(dat.n, arma::fill::ones);
  exp_tau1 = arma::vec(dat.MeanPenalties.n_elem);
  var_tau1 = arma::vec(dat.MeanPenalties.n_elem);
  exp_tau2 = arma::vec(dat.VarPenalties.n_elem);
  var_tau2 = arma::vec(dat.VarPenalties.n_elem);
  exp_eta = arma::mat(dat.n, dat.K, arma::fill::randn);
  var_eta = arma::cube(dat.K, dat.K, dat.n, arma::fill::ones);
  prec_eta = arma::cube(dat.K, dat.K, dat.n, arma::fill::ones);
}
Data::Data(arma::mat& y, arma::mat& x, arma::mat& z, arma::mat& b,
           arma::uword dim, arma::field<arma::mat>& mp,
           arma::field<arma::mat>& vp,
           arma::uvec& mi,
           arma::uvec& vi) {
  D1 = x.n_cols;
  D2 = z.n_cols;
  p = b.n_cols;
  K = dim;
  n = y.n_rows;
  Y = y;
  X = x;
  Z = z;
  B = b;
  Y = y;
  missing = arma::find_nonfinite(Y);
  missing_sub = armadillo_modulus3(missing, Y.n_rows);
  missing_time = arma::floor(missing / Y.n_rows);
  Y.elem(missing).fill(0);
  MeanPenalties = mp;
  VarPenalties = vp;
  MeanIndices = mi;
  VarIndices = vi;
}

Transformations::Transformations(Data& dat) {
  // initialize Y
  arma::mat fit(dat.n, dat.B.n_cols, arma::fill::zeros);
  arma::uword n_unique_mean = arma::size(arma::unique(dat.MeanIndices))(0);
  arma::uword n_unique_var = arma::size(arma::unique(dat.VarIndices))(0);
  C_lambda = arma::cube(dat.p * dat.D2, dat.p * dat.D2, dat.K);
  mychol_lambda = arma::cube(dat.p * dat.D2, dat.p * dat.D2, dat.K);
  g_lambda = arma::mat(dat.p * dat.D2, dat.K);
  w_lambda = arma::mat(dat.p * dat.D2, dat.K);
  mu_lambda = arma::mat(dat.p * dat.D2, dat.K);
  partialZ_eta = arma::cube(dat.p, dat.n, dat.K);
  ZetaZ = arma::cube(dat.D2, dat.D2, dat.K);
  C_beta = arma::mat(dat.p * dat.D1, dat.p * dat.D1);
  mychol_beta = arma::mat(dat.p * dat.D1, dat.p * dat.D1);
  g_beta = arma::vec(dat.D1 * dat.p);
  w_beta = arma::vec(dat.D1 * dat.p);
  mu_beta = arma::vec(dat.D1 * dat.p);
  arma::uvec seq_mean(n_unique_mean + 1);
  arma::uvec seq_var(n_unique_var + 1);
  lambda_tilde = arma::mat(dat.p, dat.K);
  arma::uword Meanold_index = 0;
  seq_mean(0) = 0;
  for (arma::uword u = 0; u < dat.MeanPenalties.n_elem; u++) {
    if (dat.MeanIndices(u) != Meanold_index) {
      seq_mean(dat.MeanIndices(u)) = dat.MeanPenalties(u).n_rows;
    }
    Meanold_index = dat.MeanIndices(u);
  }
  
  arma::uword Varold_index = 0;
  seq_var(0) = 0;
  for (arma::uword u = 0; u < dat.VarPenalties.n_elem; u++) {
    if (dat.VarIndices(u) != Varold_index) {
      seq_var(dat.VarIndices(u)) = dat.VarPenalties(u).n_rows;
    }
    Varold_index = dat.VarIndices(u);
  }
  BtB = dat.B.t() * dat.B;
  zzt = dat.Z * dat.Z.t();
}
void Transformations::initialize_transformations(Data& dat, Parameters& pars) {

  BtY = dat.B.t() * dat.Y.t();
  pars.exp_prec.fill(10000);
  // initialize fit
  fit = dat.X * pars.exp_beta.t() * dat.B.t();
  for (arma::uword k = 0; k < dat.K; k++) {
    fit = fit + arma::diagmat(pars.exp_eta.col(k)) * 
      dat.Z * pars.exp_lambda.slice(k).t() * dat.B.t();
  }
  
  // initialize ZetaZ
  arma::vec eta_tube;
  for (arma::uword k = 0; k < dat.K; k++) {
    eta_tube = pars.var_eta.tube(k, k);
    ZetaZ.slice(k) = dat.Z.t() * arma::diagmat(pars.exp_prec) *
      arma::diagmat(arma::square(pars.exp_eta.col(k)) +
      eta_tube) * dat.Z;
  }

}

void Data::impute_response(Transformations& transf) {
  for (arma::uword i = 0; i < missing_sub.n_elem; i++) {
    Y(missing_sub(i), missing_time(i)) = 
      transf.fit(missing_sub(i), missing_time(i));
  }
}

void Parameters::update_lambda(Data& dat, Transformations& transf) {
  
  
  for (arma::uword k = 0; k < dat.K; k++) {
    transf.C_lambda.slice(k) = arma::kron(transf.ZetaZ.slice(k), transf.BtB);

    transf.g_lambda.col(k) = arma::vectorise(dat.B.t() * (dat.Y.t() -
      transf.fit.t() + dat.B * exp_lambda.slice(k) * 
      dat.Z.t() * arma::diagmat(exp_eta.col(k))) * arma::diagmat(exp_prec) * 
      arma::diagmat(exp_eta.col(k)) * dat.Z);

    if (!transf.C_lambda.slice(k).is_symmetric()) {
      transf.C_lambda.slice(k) = arma::symmatu(transf.C_lambda.slice(k));
    }
    var_lambda.slice(k) = arma::inv_sympd(transf.C_lambda.slice(k));
    transf.mu_lambda = var_lambda.slice(k) * transf.g_lambda.col(k);
    
    // updated fit
    transf.fit = transf.fit - arma::diagmat(exp_eta.col(k)) * dat.Z * 
      exp_lambda.slice(k).t() * dat.B.t();
    exp_lambda.slice(k) =
      arma::reshape(transf.mu_lambda, dat.p, dat.D2);
    transf.fit = transf.fit + arma::diagmat(exp_eta.col(k)) * dat.Z *
      exp_lambda.slice(k).t() * dat.B.t();
  }
}

void Parameters::update_eta(Data& dat, Transformations& transf) {
  for (arma::uword i = 0; i < dat.n; i++) {
    for (arma::uword k = 0; k < dat.K; k++) {
      transf.lambda_tilde.col(k) = exp_lambda.slice(k) * dat.Z.row(i).t();
    }
    prec_eta.slice(i) = exp_prec(i) * transf.lambda_tilde.t() * transf.BtB
      * transf.lambda_tilde + arma::eye<arma::mat>(dat.K, dat.K);
    
    for (arma::uword k = 0; k < dat.K; k++) {
      prec_eta(k, k, i) = prec_eta(k, k, i) + 
        arma::trace(arma::kron(dat.Z.row(i).t() *
         dat.Z.row(i), transf.BtB) * var_lambda.slice(k));
    }
     
    prec_eta.slice(i) = arma::symmatu(prec_eta.slice(i));
    
    var_eta.slice(i) = arma::inv_sympd(prec_eta.slice(i));
    for (arma::uword k = 0; k < dat.K; k++) {
      transf.lambda_tilde.col(k) = exp_lambda.slice(k) * dat.Z.row(i).t();
    }
    exp_eta.row(i) = (var_eta.slice(i) *
      transf.lambda_tilde.t() * dat.B.t() *  
      exp_prec(i) * (dat.Y.row(i).t() - dat.B * exp_beta * dat.X.row(i).t())).t();
  }
  
  
  // updated fit
  transf.fit = dat.X * exp_beta.t() * dat.B.t();
  for (arma::uword k = 0; k < dat.K; k++) {
    transf.fit = transf.fit + arma::diagmat(exp_eta.col(k)) * 
      dat.Z * exp_lambda.slice(k).t() * dat.B.t();
  }
  
  // updated ZetaZ
  arma::vec eta_tube;
  for (arma::uword k = 0; k < dat.K; k++) {
    eta_tube = var_eta.tube(k, k);
    transf.ZetaZ.slice(k) = dat.Z.t() * arma::diagmat(exp_prec) *
      arma::diagmat(arma::square(exp_eta.col(k)) +
      eta_tube) * dat.Z;
  }
}

void Parameters::update_beta(Data& dat, Transformations& transf) {

  
  transf.C_beta = arma::kron(dat.X.t() * arma::diagmat(exp_prec) * dat.X,
                             transf.BtB);
  var_beta = arma::inv_sympd(transf.C_beta);
  transf.mu_beta = var_beta * arma::vectorise(dat.B.t() * 
    (dat.Y.t() - transf.fit.t() + dat.B * exp_beta * dat.X.t()) * arma::diagmat(exp_prec) * dat.X);
  
  // updated fit
  transf.fit = transf.fit - dat.X * exp_beta.t() * dat.B.t();
  exp_beta = arma::reshape(transf.mu_beta, dat.p, dat.D1);
  transf.fit = transf.fit + dat.X * exp_beta.t() * dat.B.t();

}

void Parameters::update_prec(Data& dat, Transformations& transf) {
  
}
// [[Rcpp::export]]
Rcpp::List my_main(arma::mat& Y, arma::mat& X, arma::mat& Z,
             arma::mat& B, arma::uword K,
             arma::field<arma::mat>& MeanPenalties,
             arma::field<arma::mat>& VarPenalties,
             arma::uvec& MeanIndices,
             arma::uvec& VarIndices) {
  Data dat(Y, X, Z, B, K, MeanPenalties, VarPenalties, MeanIndices, VarIndices);
  Parameters pars(dat);
  Transformations transf(dat);
  transf.initialize_transformations(dat, pars);
  for (arma::uword i = 0; i < 100; i++) {
    Rcpp::Rcout << i << std::endl;
    pars.update_beta(dat, transf);
    Rcpp::Rcout << "update_beta passed" << std::endl;
    pars.update_lambda(dat, transf);
    Rcpp::Rcout << "update_lambda passed" << std::endl;
    pars.update_eta(dat, transf);
    Rcpp::Rcout << "update_eta passed" << std::endl;
  }
  /*
  transf.update_transformations(dat, pars);
  pars.update_beta(dat, transf);
  pars.update_lambda(dat, transf);
  pars.update_eta(dat, transf);
   */
  return(Rcpp::List::create(Rcpp::Named("Beta", pars.exp_beta),
                     Rcpp::Named("Lambda", pars.exp_lambda),
                     Rcpp::Named("Eta", pars.exp_eta)));
}
