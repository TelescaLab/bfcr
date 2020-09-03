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



// [[Rcpp::export]]
Rcpp::List run_mcmc(arma::mat response, arma::mat design_mean,
               arma::mat design_var, arma::mat basis,
               arma::field<arma::mat> penalties_mean,
               arma::field<arma::mat> penalties_var,
               arma::uvec indices_mean, arma::uvec indices_var,
               arma::uword kdim, arma::uword iter, arma::uword burnin,
               arma::uword thin=1, std::string var="unequal") {
  
  Data dat(response, design_mean,
           design_var, basis,
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
