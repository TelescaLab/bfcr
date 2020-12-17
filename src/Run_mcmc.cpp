#include <RcppArmadillo.h>
#include <iostream>
#include "Utils.h"
#include "Data.h"
#include "Parameters.h"
#include "Transformations.h"
#include "GaussTransformations.h"
#include "Sampler.h"
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]


//' Run Markov-Chain Monte-Carlo
//'
//' Generate samples from the posterior distribution. This
//' algorithm pre-dominantly uses Gibbs sampling
//' @param response N x T response matrix, where N is number of subjects and T is number
//' of time points. Values can be NA if there's missing data
//' @param design_mean N x d_{1} design matrix for mean structure
//' @param design_var N x d_{2} design matrix for the covariance structure
//' @param basis User generated basis matrix
//' @param time vector of time points
//' @param penalties_mean List of smoothing penalties for mean structure
//' @param penalties_var List of smoothin penalties for covariance structure
//' @param indices_mean Maps penalties in the mean structure to beta coefficients
//' @param indices_var Maps penalties in the covariance structure to lambda
//' coefficients
//' @param kdim Dimension of latent subspace
//' @param iter Number of iterations to run
//' @param burnin Number of iterations to use as burn-in. This is only relevant when
//' passing the returned object into post-processing functions
//' @param thin Thinning defaulting to 1
//' @param var Can be set to "unequal" to estimate subject-specific measurement
//' errors or "pooled" to estimate a pooled measurement error variance
//' @export run_mcmc
//' @return A List containing 3 lists including data, control, and samples.
// [[Rcpp::export]]
Rcpp::List run_mcmc(arma::mat response, arma::mat design_mean,
                    arma::mat design_var, arma::mat basis, arma::vec time,
                    arma::field<arma::mat> penalties_mean,
                    arma::field<arma::mat> penalties_var,
                    arma::uvec indices_mean, arma::uvec indices_var,
                    arma::uword kdim, arma::uword iter, arma::uword burnin,
                    arma::uword thin=1, std::string var="unequal") {
  Data dat(response, design_mean,
           design_var, basis, time,
           penalties_mean, penalties_var,
           indices_mean, indices_var, kdim,
           iter, burnin, thin);
  Parameters pars(dat);
  Transformations transf(dat, pars);
  Sampler* mysampler = SamplerFactory::new_mcmc(var, dat, pars, transf);
  mysampler->sample_parameters();
  Rcpp::List return_me;
  return_me["data"] = mysampler->write_data();
  return_me["samples"] = mysampler->get_samples();
  return_me["control"] = mysampler->write_control();
  return(return_me);
}
